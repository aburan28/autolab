#!/usr/bin/env sage -python
"""N611S: materialize the cubic H018 descended principal polarization map."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


P = 103
ROOT = Path("notes")
N609I = ROOT / "poincare_global_class_n609i.json"
N611Q = ROOT / "poincare_cubic_product_quotient_n611q.json"
N611Q_VERIFY = ROOT / "poincare_cubic_product_quotient_n611q_independent_verifier.json"
N611R = ROOT / "poincare_cubic_quotient_polarization_n611r.json"
N611R_VERIFY = ROOT / "poincare_cubic_quotient_polarization_n611r_independent_verifier.json"


@dataclass(frozen=True)
class CM:
    """Element a+b*pi in Z[pi], pi^2+5*pi+103=0."""

    a: int = 0
    b: int = 0

    def __add__(self, other):
        return CM(self.a + other.a, self.b + other.b)

    def __neg__(self):
        return CM(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        return CM(self.a * other.a - 103 * self.b * other.b,
                  self.a * other.b + self.b * other.a - 5 * self.b * other.b)

    def conj(self):
        return CM(self.a - 5 * self.b, -self.b)

    def norm(self):
        value = self * self.conj()
        if value.b:
            raise ValueError("CM norm did not descend to Z")
        return value.a

    def text(self):
        if not self.b:
            return str(self.a)
        if not self.a:
            return "pi" if self.b == 1 else ("-pi" if self.b == -1 else f"{self.b}*pi")
        sign = "+" if self.b > 0 else "-"
        coeff = "" if abs(self.b) == 1 else f"{abs(self.b)}*"
        return f"{self.a}{sign}{coeff}pi"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text(encoding="ascii"))


def receipt_matches(receipt, primary_path):
    if receipt.get("verified") is not True:
        return False
    if receipt.get("primary_sha256") == sha256(primary_path):
        return True
    primary = receipt.get("primary")
    return isinstance(primary, dict) and primary.get("path") == str(primary_path)


def multiply(left, right):
    return tuple(tuple(sum((left[i][k] * right[k][j] for k in range(2)), CM()) for j in range(2)) for i in range(2))


def dagger(matrix):
    return tuple(tuple(matrix[j][i].conj() for j in range(2)) for i in range(2))


def equal(left, right):
    return all(left[i][j] == right[i][j] for i in range(2) for j in range(2))


def matrix_text(matrix):
    return [[entry.text() for entry in row] for row in matrix]


def frobenius(point):
    return point if point.is_zero() else point.curve()(point[0] ** P, point[1] ** P)


def point_text(point):
    return "O" if point.is_zero() else f"({point[0]}, {point[1]})"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    n609i, n611q, n611q_verify = load(N609I), load(N611Q), load(N611Q_VERIFY)
    n611r, n611r_verify = load(N611R), load(N611R_VERIFY)

    field = GF(P**3, name="a")
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    ring = PolynomialRing(field, "T")
    roots = sorted(root for root, _ in (ring.gen()**3 + ring.gen() + 24).roots())
    r = curve(roots[0], 0)
    r1, r2 = frobenius(r), frobenius(frobenius(r))
    phi = curve.isogeny(r)
    phi_dual = phi.dual()

    pi = CM(0, 1)
    h = ((CM(9), CM(2, 1)), (CM(-3, -1), CM(11)))
    t_inverse = ((CM(1), -pi), (CM(), CM(1)))
    k = multiply(multiply(dagger(t_inverse), h), t_inverse)
    expected_k = ((CM(9), CM(2, -8)), (CM(42, 8), CM(742)))
    u = CM(21, 4)
    u_dual = u.conj()
    q = CM(2, -8)
    beta_degree = 2 * u.norm()
    m_determinant = 9 * 371 - beta_degree

    # beta=phi o u and beta^dagger phi=2*u^dagger are tested on a fixed
    # mix of base-field and cubic two-torsion points.
    base_curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    base_point = next(point for point in base_curve.points() if not point.is_zero())
    base_point = curve(base_point)
    samples = {f"base_multiple_{n}": n * base_point for n in range(1, 19)}
    samples.update({"r": r, "r1": r1, "r2": r2, "origin": curve(0)})

    def u_map(point):
        return 21 * point + 4 * frobenius(point)

    def u_dual_map(point):
        return point - 4 * frobenius(point)

    point_checks = {}
    for name, point in samples.items():
        beta_dual_phi = u_dual_map(phi_dual(phi(point)))
        required = 2 * point - 8 * frobenius(point)
        point_checks[name] = {
            "input": point_text(point),
            "passes": beta_dual_phi == required,
        }
    identity_control = all(u_dual_map(phi_dual(phi(point))) == 2 * point - 8 * frobenius(point) for point in samples.values())
    negative_control = any(phi_dual(phi(point)) != 2 * point - 8 * frobenius(point) for point in samples.values())

    gates = {
        "n611q_explicit_g_quotient_receipt_binds": receipt_matches(n611q_verify, N611Q) and n611q["preflight_pass"] is True and next(row for row in n611q["records"]["rows"] if row["line"] == "G_g")["factor_degree"] == 2,
        "n611r_numerical_principal_receipt_binds": receipt_matches(n611r_verify, N611R) and n611r["preflight_pass"] is True and n611r["records"]["descended_self_intersection"] == 2,
        "h018_source_form_binds": n609i["records"]["determinant"] == 2 and n609i["records"]["self_intersection"] == 4,
        "cm_transport_is_exact": equal(k, expected_k),
        "u_dual_and_offdiagonal_match": u_dual == CM(1, -4) and CM(2) * u_dual == q,
        "explicit_dual_isogeny_map_identity_replays": identity_control,
        "identity_negative_control_rejects_wrong_offdiagonal": negative_control,
        "explicit_target_form_is_principal_positive": beta_degree == 3338 and m_determinant == 1 and 9 > 0,
    }
    output = {
        "schema": "ecdlp.h018.cubic-explicit-principal-polarization.n611s.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__)),
        "claim_status": "RESTRICTED THEOREM / EXPLICIT_CUBIC_PRINCIPAL_POLARIZATION_HOMOMORPHISM / STANDARD_DUAL_ISOGENY_IDENTITY_BOUND / MODEL-BOUND / TOY-EVIDENCE / THETA_DIVISOR_OPEN / NO_ECDLP_CLAIM",
        "assumptions": [
            "The standard dual-isogeny identity phi_g^dagger phi_g=[2] holds for the separable degree-two Velu quotient.",
            "The standard Hermitian-matrix description of polarizations on this CM product applies to the displayed maps.",
            "This receipt constructs a polarization homomorphism; it does not choose its associated line bundle or theta divisor.",
        ],
        "bound_inputs": {str(path): sha256(path) for path in (N609I, N611Q, N611Q_VERIFY, N611R, N611R_VERIFY)},
        "records": {
            "curve": "y^2=x^3+x+24 over F_(103^3)",
            "quotient_factor": "phi_g:E -> E_g with kernel <R>",
            "source_form_H": matrix_text(h),
            "transport_inverse_T": matrix_text(t_inverse),
            "transported_form_K": matrix_text(k),
            "target_polarization_M": [["9", "beta^dagger"], ["beta", "371"]],
            "u": u.text(),
            "u_dual": u_dual.text(),
            "beta": "phi_g o ([21]+[4]Frob)",
            "beta_dual_phi": q.text(),
            "u_norm": u.norm(),
            "beta_degree": beta_degree,
            "target_determinant": m_determinant,
            "point_checks": point_checks,
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "In the cubic F_(103^3) fixture, the N611Q quotient factor phi_g and beta=phi_g o ([21]+[4]Frob) explicitly realize the descended Hermitian principal polarization M=[[9,beta^dagger],[beta,371]]: its pullback transport is H018 and its determinant is one, conditional on the stated standard dual-isogeny and CM-polarization identities.",
        "next_requirement": "Construct an associated ample line bundle or theta divisor for M, then pull it back through Psi_g and compute the Frobenius-conjugate divisor data before any scalar H018 section or cryptanalytic relation claim.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611S explicit principal-polarization preflight failed")
    print(json.dumps({"output": str(args.out), "checks": sum(gates.values())}, sort_keys=True))


if __name__ == "__main__":
    main()
