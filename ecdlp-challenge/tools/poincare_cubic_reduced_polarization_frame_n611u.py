#!/usr/bin/env sage -python
"""N611U: reduce the explicit cubic H018 principal-polarization frame."""
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
N611P = ROOT / "poincare_cubic_theta_kernel_n611p.json"
N611P_VERIFY = ROOT / "poincare_cubic_theta_kernel_n611p_independent_verifier.json"
N611Q = ROOT / "poincare_cubic_product_quotient_n611q.json"
N611Q_VERIFY = ROOT / "poincare_cubic_product_quotient_n611q_independent_verifier.json"
N611S = ROOT / "poincare_cubic_explicit_principal_polarization_n611s.json"
N611S_VERIFY = ROOT / "poincare_cubic_explicit_principal_polarization_n611s_independent_verifier.json"


@dataclass(frozen=True)
class CM:
    a: int = 0
    b: int = 0

    def __add__(self, other): return CM(self.a + other.a, self.b + other.b)
    def __neg__(self): return CM(-self.a, -self.b)
    def __sub__(self, other): return self + (-other)
    def __mul__(self, other): return CM(self.a * other.a - 103 * self.b * other.b, self.a * other.b + self.b * other.a - 5 * self.b * other.b)
    def conj(self): return CM(self.a - 5 * self.b, -self.b)
    def norm(self): return (self * self.conj()).a
    def text(self):
        if self.b == 0: return str(self.a)
        if self.a == 0: return "pi" if self.b == 1 else ("-pi" if self.b == -1 else f"{self.b}*pi")
        sign = "+" if self.b > 0 else "-"
        coefficient = "" if abs(self.b) == 1 else f"{abs(self.b)}*"
        return f"{self.a}{sign}{coefficient}pi"


def sha256(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path): return json.loads(path.read_text(encoding="ascii"))
def receipt_matches(receipt, path): return receipt.get("verified") is True and receipt.get("primary_sha256") == sha256(path)
def mul(left, right): return tuple(tuple(sum((left[i][k] * right[k][j] for k in range(2)), CM()) for j in range(2)) for i in range(2))
def dagger(matrix): return tuple(tuple(matrix[j][i].conj() for j in range(2)) for i in range(2))
def equal(left, right): return all(left[i][j] == right[i][j] for i in range(2) for j in range(2))
def text_matrix(matrix): return [[entry.text() for entry in row] for row in matrix]
def frobenius(point): return point if point.is_zero() else point.curve()(point[0] ** P, point[1] ** P)


def lower_coefficient(a, b):
    return 9 * a * a - 45 * a * b + 927 * b * b + a - 196 * b + 11


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    n611p, n611p_v = load(N611P), load(N611P_VERIFY)
    n611q, n611q_v = load(N611Q), load(N611Q_VERIFY)
    n611s, n611s_v = load(N611S), load(N611S_VERIFY)

    pi, d = CM(0, 1), CM(2, 1)
    h = ((CM(9), CM(2, 1)), (CM(-3, -1), CM(11)))
    t_inv = ((CM(1), -d), (CM(), CM(1)))
    k = mul(mul(dagger(t_inv), h), t_inv)
    expected_k = ((CM(9), CM(-16, -8)), (CM(24, 8), CM(690)))
    u, u_dual, q = CM(12, 4), CM(-8, -4), CM(-16, -8)

    field = GF(P**3, name="a")
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    ring = PolynomialRing(field, "T")
    roots = sorted(root for root, _ in (ring.gen()**3 + ring.gen() + 24).roots())
    r = curve(roots[0], 0)
    phi, phi_dual = curve.isogeny(r), curve.isogeny(r).dual()
    generator = (frobenius(r), r)
    reduced_image = (generator[0] + 2 * generator[1] + frobenius(generator[1]), phi(generator[1]))
    wrong_image = (generator[0], phi(generator[1]))
    base_curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    base = curve(next(point for point in base_curve.points() if not point.is_zero()))
    samples = [n * base for n in range(1, 19)] + [r, frobenius(r), frobenius(frobenius(r)), curve(0)]
    point_identity = all((point - 4 * frobenius(point)) * 2 == 2 * point - 8 * frobenius(point) for point in samples)
    beta_identity = all((u_dual_map := (lambda point: -8 * point - 4 * frobenius(point)))(phi_dual(phi(point))) == -16 * point - 8 * frobenius(point) for point in samples)

    # Complete the square in A.  For odd |B|>=3 the real lower bound already
    # exceeds 690; B=-1 is separately larger, and B=1 is minimized at A=2
    # among even A.
    # With A=2k, the B=1 and B=-1 quadratics have real vertices at
    # k=11/9 and k=-23/18, respectively. Their integer minima are k=1 and
    # k=-1, with values 690 and 1078.
    b_one_exact_minimum = lower_coefficient(2, 1)
    b_minus_one_exact_minimum = lower_coefficient(-2, -1)
    # After minimizing over real A, the remaining lower bound is the convex
    # quadratic (31347*B^2-6966*B+395)/36. Its vertex lies in (0,1), so
    # among odd integers with |B|>=3 the minimum is attained at B=3.
    continuous_bound_b_plus_three = (31347 * 3 * 3 - 6966 * 3 + 395) / 36
    continuous_bound_b_minus_three = (31347 * 3 * 3 + 6966 * 3 + 395) / 36
    continuous_bound_abs_b_ge_3 = min(continuous_bound_b_plus_three, continuous_bound_b_minus_three)
    gates = {
        "n611p_and_n611q_kernel_receipts_bind": n611p["preflight_pass"] and receipt_matches(n611p_v, N611P) and n611q["preflight_pass"] and receipt_matches(n611q_v, N611Q),
        "n611s_baseline_receipt_binds": n611s["preflight_pass"] and receipt_matches(n611s_v, N611S) and n611s["records"]["target_polarization_M"][1][1] == "371",
        "d_congruent_pi_mod_two_preserves_kernel_generator": reduced_image[0].is_zero() and reduced_image[1].is_zero(),
        "noncongruent_zero_control_fails_kernel": not (wrong_image[0].is_zero() and wrong_image[1].is_zero()),
        "exact_reduced_cm_transport": equal(k, expected_k),
        "dual_map_offdiagonal_identity_replays": beta_identity and point_identity and u.conj() == u_dual and CM(2) * u_dual == q,
        "reduced_form_is_principal_positive": 2 * u.norm() == 3104 and 9 * 345 - 3104 == 1 and 9 > 0,
        "triangular_kernel_preserving_minimum_is_690": b_one_exact_minimum == 690 and b_minus_one_exact_minimum > 690 and continuous_bound_abs_b_ge_3 > 690,
    }
    output = {
        "schema": "ecdlp.h018.cubic-reduced-principal-polarization-frame.n611u.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__)),
        "claim_status": "RESTRICTED THEOREM / KERNEL_PRESERVING_TRIANGULAR_CUBIC_POLARIZATION_REDUCTION / MODEL-BOUND / TOY-EVIDENCE / THETA_EQUATIONS_OPEN / NO_ECDLP_CLAIM",
        "assumptions": ["The standard dual-isogeny and CM-polarization identities used in N611S apply.", "The reduction minimum is only over triangular source changes T_d with d=pi mod 2."],
        "bound_inputs": {str(path): sha256(path) for path in (N611P, N611P_VERIFY, N611Q, N611Q_VERIFY, N611S, N611S_VERIFY)},
        "records": {"d": d.text(), "transported_form": text_matrix(k), "u": u.text(), "u_dual": u_dual.text(), "beta": "phi_g o ([12]+[4]Frob)", "beta_degree": 2 * u.norm(), "target_form": [["9", "beta^dagger"], ["beta", "345"]], "target_determinant": 1, "old_second_diagonal": 371, "new_second_diagonal": 345, "old_beta_degree": 3338, "new_beta_degree": 3104, "lower_coefficient_formula": "9*A^2-45*A*B+927*B^2+A-196*B+11", "triangular_minimizer": {"A": 2, "B": 1, "lower_coefficient": 690}, "minimum_proof_controls": {"B_plus_one_exact_minimum": b_one_exact_minimum, "B_minus_one_exact_minimum": b_minus_one_exact_minimum, "continuous_lower_bound_B_plus_three": continuous_bound_b_plus_three, "continuous_lower_bound_B_minus_three": continuous_bound_b_minus_three, "continuous_lower_bound_abs_B_ge_3": continuous_bound_abs_b_ge_3}},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "Within all triangular source changes T_d with d=pi modulo 2, d=2+pi is the exact minimizer of the transported lower coefficient. It realizes the same cubic H018 quotient kernel with principal target form [[9,beta^dagger],[beta,345]] and beta=phi_g o ([12]+[4]Frob), reducing the cross-map degree from 3338 to 3104.",
        "next_requirement": "Use this reduced line-bundle frame to compute an evaluable theta generator, or leave the triangular family for a full quotient-lattice reduction before attempting Riemann-Roch.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N611U reduced-frame preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__": main()
