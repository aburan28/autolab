#!/usr/bin/env sage -python
"""N608M: exact degree screen for the N608J integral correspondence lattice."""
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


def degree(a, b):
    return 99 * a * a + 27 * a * b + 81 * b * b - 195 * a - 9 * b + 99


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi_dual = second * first
    pi = pi_dual.dual()
    cover = pi.domain()
    phi_11 = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    phi_11 = phi_11.codomain().isomorphism_to(target) * phi_11
    points = list(cover.points())
    combination_controls = {17: (1, -1), 23: (-1, -1), 41: (1, -2)}
    control_results = {}
    for prime_degree, (a, b) in combination_controls.items():
        candidate = next(item for item in cover.isogenies_prime_degree(prime_degree) if item.codomain().is_isomorphic(target))
        candidate = candidate.codomain().isomorphism_to(target) * candidate
        direct = [candidate(point) for point in points]
        combined = [a * phi_11(point) + b * pi(point) for point in points]
        same = direct == combined or [-value for value in direct] == combined
        control_results[str(prime_degree)] = {
            "combination": [a, b],
            "h018_pullback_degree": degree(a, b),
            "pointwise_equal_up_to_codomain_sign": same,
        }
    search = [(degree(a, b), a, b) for a in range(-48, 49) for b in range(-48, 49)]
    positive_values = sorted(value for value in search if value[0] > 0)
    records = {
        "cover_order": int(cover.cardinality()),
        "pi_degree": int(pi.degree()),
        "phi_11_degree": int(phi_11.degree()),
        "degree_polynomial": "99*a^2 + 27*a*b + 81*b^2 - 195*a - 9*b + 99",
        "modulo_nine": "degree(a,b) == 3*a (mod 9)",
        "degree_two_congruence_obstruction": "3*a == 2 (mod 9) has no integer solution",
        "positive_control": {"a": 1, "b": 0, "degree": degree(1, 0)},
        "minimum_positive_search_window": {"degree": positive_values[0][0], "a": positive_values[0][1], "b": positive_values[0][2]},
        "degree_two_solutions_search_window": [[a, b] for value, a, b in search if value == 2],
        "known_map_combination_controls": control_results,
    }
    gates = {
        "n608l_positive_control_degree_three": degree(1, 0) == 3,
        "degree_two_rejected_modulo_nine": all((3 * a - 2) % 9 != 0 for a in range(9)),
        "no_degree_two_in_wide_exact_window": not records["degree_two_solutions_search_window"],
        "known_maps_replay_as_lattice_combinations": all(item["pointwise_equal_up_to_codomain_sign"] for item in control_results.values()),
        "positive_degree_minimum_is_three_in_window": records["minimum_positive_search_window"]["degree"] == 3,
    }
    output = {
        "schema": "ecdlp.h018.poincare-correspondence-lattice.n608m.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "NEGATIVE RESULT / DECLARED INTEGRAL CORRESPONDENCE LATTICE ONLY / MODEL-BOUND / TOY-EVIDENCE / OTHER CORRESPONDENCES OR MODIFICATIONS OPEN / NO_ECDLP_CLAIM",
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The declared integral lattice u_(a,b)=((a phi_11+b pi),pi) has no H018 pullback of degree two: its degree is congruent to 3a modulo 9. Its least positive sampled value is the exact degree-three N608J correspondence.",
        "next_requirement": "Search a correspondence outside this declared lattice or build a genuine elementary modification with normalized-Poincare charts; do not cure the mismatch by pointwise weights.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608M lattice degree audit failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
