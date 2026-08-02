#!/usr/bin/env sage -python
"""N610D: measure coefficient-wise exceptional-slope degree growth for H018."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import time
from pathlib import Path

from sage.all import Matrix, vector

from poincare_determinantal_section_n609j import moving_basis
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def lead(series):
    return series[series.valuation()]


def replay_degree(field, samples, degree):
    matrix = Matrix(field, [[slope**power for power in range(degree + 1)] for slope, _ in samples[: degree + 1]])
    coefficients = matrix.solve_right(vector(field, [value for _, value in samples[: degree + 1]]))
    return all(sum(coefficients[power] * slope**power for power in range(degree + 1)) == value for slope, value in samples[degree + 1:])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--max-order", type=int, default=24)
    parser.add_argument("--max-degree", type=int, default=26)
    args = parser.parse_args()
    started = time.perf_counter()

    target, cover, pi, phi = fixture()
    field, ring, u, support, matrix, xvalues = origin_data(target, cover, pi, phi)
    determinant = matrix.det()
    cauchy_slope = lead((-support[0] / support[1]) / u)
    formal = target.formal_group()
    slopes = [field(value) for value in range(2, 32)]
    level_profiles = []

    for level in (0, 1, 17):
        moving = matrix.solve_right(xvalues - ring(level) * vector(ring, [1] * 9))
        samples = []
        for slope in slopes:
            parameter = ring(slope) * u
            point = (ring(formal.x(120)(parameter)), ring(formal.y(120)(parameter)))
            raw = sum(moving[index] * moving_basis(point, support)[index] for index in range(9))
            normal = -determinant * parameter**8 * raw
            samples.append((slope, (slope - cauchy_slope) * normal))
        profile = []
        for order in range(args.max_order + 1):
            coefficient_samples = [(slope, value[order]) for slope, value in samples]
            minimal_degree = next(
                (degree for degree in range(args.max_degree + 1) if replay_degree(field, coefficient_samples, degree)),
                None,
            )
            profile.append({"u_order": order, "minimal_replayed_slope_degree": minimal_degree})
        level_profiles.append({"level": level, "profile": profile})

    first_profile = level_profiles[0]["profile"]
    constant_degree = first_profile[0]["minimal_replayed_slope_degree"]
    higher_degree_nine_rejection = all(
        profile["profile"][10]["minimal_replayed_slope_degree"] is not None
        and profile["profile"][10]["minimal_replayed_slope_degree"] > 9
        for profile in level_profiles
    )
    odd_coefficients_vanish = all(
        row["minimal_replayed_slope_degree"] == 0
        for profile in level_profiles
        for row in profile["profile"]
        if row["u_order"] % 2 == 1
    )
    gates = {
        "constant_term_replays_at_degree_one": constant_degree == 1,
        "every_tested_coefficient_has_a_replayed_degree": all(
            row["minimal_replayed_slope_degree"] is not None for row in first_profile
        ),
        "higher_coefficients_reject_degree_nine": higher_degree_nine_rejection,
        "all_tested_odd_coefficients_vanish": odd_coefficients_vanish,
    }
    output = {
        "schema": "ecdlp.h018.poincare-slope-degree-profile.n610d.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / H018 COEFFICIENT-WISE EXCEPTIONAL-SLOPE DEGREE PROFILE / MODEL-BOUND / TOY-EVIDENCE / EXACT_TWO_VARIABLE_CHART_OPEN / NO_ECDLP_CLAIM",
        "parameters": {"max_u_order": args.max_order, "max_slope_degree": args.max_degree, "sample_slopes": [int(slope) for slope in slopes]},
        "records": {"cauchy_slope": str(cauchy_slope), "profiles": level_profiles, "elapsed_seconds": time.perf_counter() - started},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The replayed degree profile isolates increasing slope complexity at the even formal coefficients while the tested odd coefficients vanish. It distinguishes the degree-one exceptional restriction from the higher-order behavior that invalidated N610C's fixed degree-nine model; level-dependent cancellations may lower individual degrees.",
        "next_requirement": "Use the observed coefficient-degree growth to derive an exact rational two-variable expression or a rigorously bounded truncation. This profile alone gives no branch, base-scheme, Cartier, normalization, Jacobian, relation, rank, descent, cost, or ECDLP conclusion.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N610D slope-degree profile failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
