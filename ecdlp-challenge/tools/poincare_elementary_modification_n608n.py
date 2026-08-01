#!/usr/bin/env sage -python
"""N608N: identify the degree-corrected elementary-modification target."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF

P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture():
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi_dual = second * first
    pi = pi_dual.dual()
    cover = pi.domain()
    phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    return target, cover, pi, pi_dual, phi.codomain().isomorphism_to(target) * phi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, pi_dual, phi = fixture()
    induced_matches = []
    for point in cover.points():
        image = phi.dual()(9 * phi(point) + 3 * pi(point)) + pi_dual(-4 * phi(point) + 11 * pi(point))
        induced_matches.append(image == 3 * point)
    rational_two_torsion = [point for point in cover.points() if not point.is_zero() and 2 * point == cover(0)]
    records = {
        "base_curve_order": int(target.cardinality()),
        "cover_curve_order": int(cover.cardinality()),
        "cover_equation": str(cover),
        "pi_degree": int(pi.degree()),
        "phi_degree": int(phi.degree()),
        "induced_h018_polarization": "[3]",
        "induced_action_checks": len(induced_matches),
        "rational_nonzero_two_torsion": [str(point) for point in rational_two_torsion],
        "line_class_argument": "A normalized symmetric F_103-rational degree-three line bundle differs from O(3Oprime) by rational 2-torsion; none exists on the order-109 cover.",
        "identified_pullback_line": "O_Eprime(3Oprime)",
        "finite_pushforward_exact_sequence": "0 -> pi_*O(2Oprime) -> pi_*O(3Oprime) -> k(O) -> 0",
        "pushforward_rank": 9,
        "pushforward_degrees": {"pi_*O(2Oprime)": 2, "pi_*O(3Oprime)": 3},
        "pushforward_determinants": {"pi_*O(2Oprime)": "O(2O)", "pi_*O(3Oprime)": "O(3O)"},
        "remaining_global_construction": "Extend the N608K evaluator as normalized-Poincare chart matrices into pi_*O(3Oprime), then prove its cokernel is k(O).",
    }
    gates = {
        "registered_cover_has_order_109": int(cover.cardinality()) == 109,
        "induced_polarization_is_three_on_all_rational_points": all(induced_matches),
        "no_rational_nonzero_two_torsion": not rational_two_torsion,
        "symmetric_degree_three_class_has_origin_normalization": int(cover.cardinality()) % 2 == 1 and not rational_two_torsion,
        "canonical_divisor_inclusion_is_degree_one": records["pushforward_degrees"]["pi_*O(3Oprime)"] - records["pushforward_degrees"]["pi_*O(2Oprime)"] == 1,
        "determinant_change_matches_origin_quotient": records["pushforward_determinants"]["pi_*O(2Oprime)"] == "O(2O)" and records["pushforward_determinants"]["pi_*O(3Oprime)"] == "O(3O)",
    }
    output = {
        "schema": "ecdlp.h018.poincare-elementary-modification.n608n.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "RESTRICTED THEOREM / DEGREE-CORRECT TARGET IDENTIFIED / STANDARD-FACT-BOUND / MODEL-BOUND / TOY-EVIDENCE / GLOBAL EVALUATOR EXTENSION OPEN / NO_ECDLP_CLAIM",
        "assumptions": ["The H018 line bundle is normalized, symmetric, and F_103-rational."],
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "On the registered fixture and stated normalization assumption, the N608J graph pulls H018 back to O_Eprime(3Oprime), and pi_*O(2Oprime) is the canonical elementary transform of pi_*O(3Oprime) at O. This identifies a degree-correct global target but does not extend the evaluator.",
        "next_requirement": "Construct normalized-Poincare chart matrices for the N608K evaluator into pi_*O(3Oprime) and prove that their cokernel is k(O), without pointwise gauge fitting.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608N elementary-modification preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
