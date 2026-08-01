#!/usr/bin/env sage -python
"""Independent N608N verifier for the degree-corrected target receipt."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    target = EllipticCurve(GF(103), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi_dual = second * first
    pi = pi_dual.dual()
    cover = pi.domain()
    phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    phi = phi.codomain().isomorphism_to(target) * phi
    induced_matches = [
        phi.dual()(9 * point_phi + 3 * point_pi) + pi_dual(-4 * point_phi + 11 * point_pi) == 3 * point
        for point in cover.points()
        for point_phi, point_pi in [(phi(point), pi(point))]
    ]
    torsion = [point for point in cover.points() if not point.is_zero() and 2 * point == cover(0)]
    records = primary["records"]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "cover_order_recomputes": int(cover.cardinality()) == 109 == records["cover_curve_order"],
        "no_rational_two_torsion_recomputes": not torsion and records["rational_nonzero_two_torsion"] == [],
        "induced_action_recomputes": all(induced_matches) and len(induced_matches) == int(cover.cardinality()) == records["induced_action_checks"],
        "line_target_is_degree_three": records["identified_pullback_line"] == "O_Eprime(3Oprime)",
        "elementary_modification_sequence_is_declared": records["finite_pushforward_exact_sequence"] == "0 -> pi_*O(2Oprime) -> pi_*O(3Oprime) -> k(O) -> 0",
        "rank_degree_delta_matches": records["pushforward_rank"] == 9 and records["pushforward_degrees"]["pi_*O(3Oprime)"] - records["pushforward_degrees"]["pi_*O(2Oprime)"] == 1,
        "determinant_delta_matches": records["pushforward_determinants"] == {"pi_*O(2Oprime)": "O(2O)", "pi_*O(3Oprime)": "O(3O)"},
    }
    output = {
        "schema": "ecdlp.h018.poincare-elementary-modification.n608n.independent-verifier.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N608N independent verification failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
