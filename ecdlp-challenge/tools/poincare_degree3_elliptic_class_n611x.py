#!/usr/bin/env sage -python
"""N611X: materialize a degree-three elliptic divisor class for reduced M."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


P = 103
ROOT = Path("notes")
N611U = ROOT / "poincare_cubic_reduced_polarization_frame_n611u.json"
N611U_VERIFY = ROOT / "poincare_cubic_reduced_polarization_frame_n611u_independent_verifier.json"


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
        sign = "+" if self.b > 0 else "-"
        coeff = "" if abs(self.b) == 1 else f"{abs(self.b)}*"
        return f"{self.a}{sign}{coeff}pi"


def sha256(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path): return json.loads(path.read_text(encoding="ascii"))
def receipt_matches(receipt, path): return receipt.get("verified") is True and receipt.get("primary_sha256") == sha256(path)
def frobenius(point): return point if point.is_zero() else point.curve()(point[0] ** P, point[1] ** P)


def intersection(a, b, u, a_other, b_other, u_other):
    # beta=phi*u, so degree(beta)=2N(u); this is det(M+D)-det(M)-det(D).
    return a * b_other + a_other * b - 2 * ((u + u_other).norm() - u.norm() - u_other.norm())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    n611u, n611u_verify = load(N611U), load(N611U_VERIFY)
    u_m, u_d = CM(12, 4), CM(13, 4)
    a_m, b_m, a_d, b_d = 9, 345, 9, 346
    det_m = a_m * b_m - 2 * u_m.norm()
    det_d = a_d * b_d - 2 * u_d.norm()
    m_dot_d = intersection(a_m, b_m, u_m, a_d, b_d, u_d)
    coordinate_degree = a_m
    primitive_gcd = math.gcd(math.gcd(a_d, b_d), math.gcd(u_d.a, u_d.b))

    field = GF(P**3, name="a")
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    ring = PolynomialRing(field, "T")
    root = sorted(value for value, _ in (ring.gen()**3 + ring.gen() + 24).roots())[0]
    phi, phi_dual = curve.isogeny(curve(root, 0)), curve.isogeny(curve(root, 0)).dual()
    base_curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    base = curve(next(point for point in base_curve.points() if not point.is_zero()))
    samples = [n * base for n in range(1, 19)] + [curve(root, 0), frobenius(curve(root, 0)), curve(0)]
    def beta_map(u, point): return phi(u.a * point + u.b * frobenius(point))
    def beta_dual_map(u, point):
        dual = u.conj()
        return dual.a * phi_dual(point) + dual.b * frobenius(phi_dual(point))
    map_replay = all(
        beta_dual_map(u_d, beta_map(u_d, point)) == 2 * u_d.norm() * point
        for point in samples
    )

    # Bounded diagnostic: all rank-one Hermitian classes with a,b<=400 and
    # beta=phi*u that were enumerated have no positive intersection below 3.
    norm_index = {}
    for coeff_b in range(-30, 31):
        for coeff_a in range(-400, 401):
            value = CM(coeff_a, coeff_b)
            norm = value.norm()
            if 0 <= norm <= 100000:
                norm_index.setdefault(norm, []).append(value)
    bounded_minimum = None
    for left in range(1, 401):
        for right in range(1, 401):
            if left * right % 2:
                continue
            for value in norm_index.get(left * right // 2, []):
                candidate = intersection(a_m, b_m, u_m, left, right, value)
                if candidate > 0 and (bounded_minimum is None or candidate < bounded_minimum):
                    bounded_minimum = candidate

    gates = {
        "n611u_reduced_polarization_receipt_binds": n611u["preflight_pass"] is True and receipt_matches(n611u_verify, N611U),
        "m_is_principal": det_m == 1,
        "candidate_is_primitive_rank_one": det_d == 0 and primitive_gcd == 1 and a_d > 0 and b_d > 0,
        "candidate_has_degree_three_against_m": m_dot_d == 3,
        "coordinate_divisor_positive_control": coordinate_degree == 9,
        "cross_map_pointwise_arithmetic_replays": map_replay,
        "bounded_lower_intersection_diagnostic": bounded_minimum == 3,
    }
    output = {
        "schema": "ecdlp.h018.degree3-elliptic-class.n611x.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__)),
        "claim_status": "OBSERVATION / PRIMITIVE_RANK_ONE_DEGREE_THREE_ELLIPTIC_DIVISOR_CLASS / STANDARD_HERMITIAN_NS_CORRESPONDENCE_BOUND / MODEL-BOUND / TOY-EVIDENCE / EFFECTIVITY_AND_THETA_CURVE_OPEN / NO_ECDLP_CLAIM",
        "assumptions": ["The standard Hermitian Neron-Severi description for E x E_g identifies the displayed primitive rank-one positive semidefinite matrix with an elliptic divisor class.", "A degree-three elliptic cover of a theta curve additionally requires effectivity and irreducibility, neither of which is proved here."],
        "bound_inputs": {str(path): sha256(path) for path in (N611U, N611U_VERIFY)},
        "records": {"M": {"diagonal": [a_m, b_m], "u": u_m.text(), "determinant": det_m}, "D": {"diagonal": [a_d, b_d], "u": u_d.text(), "beta": "phi_g o ([13]+[4]Frob)", "determinant": det_d, "primitive_gcd": primitive_gcd}, "intersection_M_D": m_dot_d, "coordinate_elliptic_degree": coordinate_degree, "bounded_rank_one_search": {"a_and_b_upper_bound": 400, "cm_coefficient_box": {"a": 400, "b": 30}, "minimum_positive_intersection": bounded_minimum}},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The reduced principal form M has an explicit primitive rank-one Hermitian class D of M-degree three. Under the stated standard Neron-Severi correspondence, this is a concrete degree-three elliptic-divisor candidate for the theta surface; effectivity, irreducibility, quotient equation, and genus-two curve construction remain open.",
        "next_requirement": "Construct the elliptic quotient associated with D and prove or disprove effectivity/irreducibility of the unique M theta divisor. A positive result should then materialize a degree-three genus-two cover and its evaluator.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N611X degree-three class preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__": main()
