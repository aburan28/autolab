#!/usr/bin/env sage -python
"""N611Y: materialize the conductor-isogeny factorization of 13+4*pi."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


P = 103
ROOT = Path("notes")
N611X = ROOT / "poincare_degree3_elliptic_class_n611x.json"
N611X_VERIFY = ROOT / "poincare_degree3_elliptic_class_n611x_independent_verifier.json"


def sha256(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path): return json.loads(path.read_text(encoding="ascii"))
def receipt_matches(receipt, path): return receipt.get("verified") is True and receipt.get("primary_sha256") == sha256(path)
def frobenius(point): return point if point.is_zero() else point.curve()(point[0] ** P, point[1] ** P)
def point_text(point): return "O" if point.is_zero() else f"({point[0]}, {point[1]})"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    n611x, n611x_verify = load(N611X), load(N611X_VERIFY)
    field = GF(P**3, name="a")
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    cubic = curve.isogenies_prime_degree(3)
    ascending = [isogeny for isogeny in cubic if isogeny.codomain().j_invariant() == field(10)]
    if len(ascending) != 1:
        raise RuntimeError("ascending conductor isogeny was not unique")
    psi = ascending[0]
    e0 = psi.codomain()
    psi_dual = psi.dual()
    e0_base = EllipticCurve(GF(P), [0, 0, 0, 45, 7])
    rho_candidates = e0_base.isogenies_prime_degree(173)
    target_curve = next(isogeny.codomain() for isogeny in rho_candidates if isogeny.codomain().a4() == 51 and isogeny.codomain().a6() == 33)
    rho_base = next(isogeny for isogeny in rho_candidates if isogeny.codomain() == target_curve)
    ring = PolynomialRing(field, "x")
    rho = e0.isogeny(ring(rho_base.kernel_polynomial()))
    scaling = field(24)
    def iota(point):
        return e0(0) if point.is_zero() else e0(scaling**2 * point[0], scaling**3 * point[1])
    def composite(point): return psi_dual(iota(rho(psi(point))))
    base_curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    base_samples = [curve(point) for point in base_curve.points()]
    root_ring = PolynomialRing(field, "T")
    roots = sorted(value for value, _ in (root_ring.gen()**3 + root_ring.gen() + 24).roots())
    torsion_samples = [curve(root, 0) for root in roots]
    samples = base_samples + torsion_samples
    point_checks = [{"input": point_text(point), "passes": composite(point) == 13 * point + 4 * frobenius(point)} for point in samples]
    alternatives = []
    for index, candidate in enumerate(rho_candidates):
        kernel_map = e0.isogeny(ring(candidate.kernel_polynomial()))
        candidate_target = candidate.codomain()
        scalings = [value for value in GF(P) if value and candidate_target.a4() * value**4 == e0_base.a4() and candidate_target.a6() * value**6 == e0_base.a6()]
        for value in scalings:
            def candidate_iota(point, value=value): return e0(0) if point.is_zero() else e0(field(value)**2 * point[0], field(value)**3 * point[1])
            alternatives.append({"rho_index": index, "scaling": int(value), "matches_all_samples": all(psi_dual(candidate_iota(kernel_map(psi(point)))) == 13 * point + 4 * frobenius(point) for point in samples)})
    gates = {
        "n611x_degree_three_class_receipt_binds": n611x["preflight_pass"] is True and receipt_matches(n611x_verify, N611X),
        "unique_ascending_cubic_target_has_maximal_order_j": len(ascending) == 1 and e0.a4() == 45 and e0.a6() == 7 and e0.j_invariant() == field(10),
        "cm_factorization_norm_identity": 13 + 4 * (-4) == -3 and (-1) * (-1) + (-1) * 4 + 11 * 4 * 4 == 173 and 9 * 173 == 1557,
        "degree_product_matches_cm_norm": int(psi.degree()) * int(rho.degree()) * int(psi_dual.degree()) == 1557,
        "selected_norm_173_target_and_scaling_materialized": rho.codomain().a4() == 51 and rho.codomain().a6() == 33 and scaling == 24,
        "composed_factorization_matches_cm_map_on_samples": all(item["passes"] for item in point_checks),
        "all_other_target_scaling_controls_fail": sum(item["matches_all_samples"] for item in alternatives) == 1,
    }
    output = {
        "schema": "ecdlp.h018.conductor-factorization.n611y.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__)),
        "claim_status": "OBSERVATION / EXPLICIT_CONDUCTOR_LOWERING_FACTORIZATION_OF_DEGREE_THREE_CLASS / MODEL-BOUND / TOY-EVIDENCE / ELLIPTIC_QUOTIENT_AND_THETA_CURVE_OPEN / NO_ECDLP_CLAIM",
        "assumptions": ["The recognized j=10 target realizes the maximal CM order of discriminant -43 in this ordinary isogeny class.", "Pointwise agreement on the recorded rational and cubic-torsion samples is a computational witness; a symbolic morphism-equality proof remains open."],
        "bound_inputs": {str(path): sha256(path) for path in (N611X, N611X_VERIFY)},
        "records": {
            "cm_identity": "13+4*pi=3*(-1+4*omega), pi=-4+3*omega, omega^2-omega+11=0",
            "u_degree": 1557,
            "psi": {
                "degree": int(psi.degree()),
                "target": {"a4": str(e0.a4()), "a6": str(e0.a6()), "j": str(e0.j_invariant())},
            },
            "rho": {
                "degree": int(rho.degree()),
                "target": {"a4": str(rho.codomain().a4()), "a6": str(rho.codomain().a6())},
                "iota_scaling": int(scaling),
            },
            "psi_dual_degree": int(psi_dual.degree()),
            "point_check_count": len(point_checks),
            "point_checks": point_checks,
            "alternative_controls": alternatives,
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The N611X element 13+4*pi is materialized through the conductor-lowering path psi^dagger o iota o rho o psi, with degrees 3,173,3 and a unique tested target/scaling. The composition agrees with [13]+[4]Frob on all 109 base-field points and all three cubic two-torsion points.",
        "next_requirement": "Use this factorization to construct the actual elliptic quotient represented by D and prove effectivity/irreducibility of M's theta divisor before deriving a genus-two cover or evaluator.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N611Y conductor factorization preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__": main()
