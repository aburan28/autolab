#!/usr/bin/env sage -python
"""N609G: screen rational removable and origin boundary points for singularities."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import GF

from poincare_global_sections_n608r import fixture, translate_by_fixed
from poincare_linear_projection_n608v import P, coefficient_table
from poincare_local_singularity_probe_n609f import (
    PRECISION,
    coefficient_at,
    p_value_and_derivative,
    q_direction_coefficients,
    q_value_and_derivative,
)
from poincare_origin_boundary_n609c import pole_degree, regular_sections
from poincare_removable_graph_n609e import regularized_value


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def l9_terms(point):
    x, y = point[0], point[1]
    return [1, x, y, x * x, x * y, x**3, x * x * y, x**4, x**3 * y]


def l9_series_terms(point):
    x, y = point
    return [1, x, y, x * x, x * y, x**3, x * x * y, x**4, x**3 * y]


def boundary_p_derivative(target, field, ring, coefficients, point):
    point_k = target.base_extend(field)(field(point[0]), field(point[1]))
    formal = target.formal_group()
    x0, y0 = ring(formal.x(PRECISION)), ring(formal.y(PRECISION))
    moving_point = translate_by_fixed(point_k, x0, y0, ring)
    value = sum(coefficients[index] * l9_series_terms(moving_point)[index] for index in range(9))
    return coefficient_at(value, 0), coefficient_at(value, 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, extension_field, sections = regular_sections()
    graph_target, cover, pi, phi = fixture()
    coefficients = coefficient_table(target)

    graph_rows = []
    for level in (0, 1, 17):
        sources = [(4 * q0, q0) for q0, coefficient in coefficients.items() if regularized_value(target, coefficient, q0) == GF(P)(level)]
        p_zero = q_zero = both_zero = mismatch = 0
        for point, q0 in sources:
            field, ring, support_series, moving = q_direction_coefficients(graph_target, cover, pi, phi, q0)
            q_value, q_derivative = q_value_and_derivative(moving, point, support_series)
            p_value, p_derivative = p_value_and_derivative(target, field, ring, coefficients[q0], point, -4 * q0)
            mismatch += int(GF(P)(q_value) != GF(P)(level) or GF(P)(p_value) != GF(P)(level))
            p_zero += int(p_derivative == 0)
            q_zero += int(q_derivative == 0)
            both_zero += int(p_derivative == 0 and q_derivative == 0)
        graph_rows.append({
            "level": level,
            "restored_graph_source_count": len(sources),
            "p_direction_zero_derivative_count": p_zero,
            "q_direction_zero_derivative_count": q_zero,
            "both_direction_zero_candidate_count": both_zero,
            "value_reconstruction_mismatch_count": mismatch,
        })

    boundary_rows = []
    field = GF(P)
    ring = sections["x"][0].parent()
    for level in (0, 1, 17):
        section_coefficients = [sections["x"][index] - field(level) * sections["1"][index] for index in range(9)]
        finite_sources = [
            point
            for point in target.points()
            if not point.is_zero() and coefficient_at(sum(section_coefficients[index] * ring(field(l9_terms(point)[index])) for index in range(9)), 0) == 0
        ]
        p_zero = q_zero = both_zero = mismatch = 0
        for point in finite_sources:
            q_value = sum(section_coefficients[index] * ring(field(l9_terms(point)[index])) for index in range(9))
            p_value, p_derivative = boundary_p_derivative(target, extension_field, ring, section_coefficients, point)
            mismatch += int(coefficient_at(q_value, 0) != 0 or p_value != 0)
            p_zero += int(p_derivative == 0)
            q_zero += int(coefficient_at(q_value, 1) == 0)
            both_zero += int(p_derivative == 0 and coefficient_at(q_value, 1) == 0)
        constant_coefficients = [field(coefficient_at(value, 0)) for value in section_coefficients]
        boundary_rows.append({
            "level": level,
            "finite_boundary_source_count": len(finite_sources),
            "p_direction_zero_derivative_count": p_zero,
            "q_direction_zero_derivative_count": q_zero,
            "both_direction_zero_candidate_count": both_zero,
            "value_reconstruction_mismatch_count": mismatch,
            "origin_zero_multiplicity_as_l9_section": 9 - pole_degree(constant_coefficients),
        })

    gates = {
        "graph_counts_and_value_reconstruction_exact": [
            (row["level"], row["restored_graph_source_count"], row["value_reconstruction_mismatch_count"])
            for row in graph_rows
        ] == [(0, 8, 0), (1, 0, 0), (17, 2, 0)],
        "no_restored_graph_singularity_candidate": all(row["both_direction_zero_candidate_count"] == 0 for row in graph_rows),
        "finite_q_origin_boundary_counts_and_value_reconstruction_exact": [
            (row["level"], row["finite_boundary_source_count"], row["value_reconstruction_mismatch_count"])
            for row in boundary_rows
        ] == [(0, 2, 0), (1, 0, 0), (17, 0, 0)],
        "no_finite_q_origin_boundary_singularity_candidate": all(row["both_direction_zero_candidate_count"] == 0 for row in boundary_rows),
        "all_q_origin_members_have_simple_p_origin_zero": all(row["origin_zero_multiplicity_as_l9_section"] == 1 for row in boundary_rows),
    }
    output = {
        "schema": "ecdlp.h018.poincare-boundary-singularity-screen.n609g.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / NO RATIONAL SINGULARITY CANDIDATE ON SCREENED H018 BOUNDARY LOCI / MODEL-BOUND / TOY-EVIDENCE / GLOBAL SMOOTHNESS OPEN / NO_ECDLP_CLAIM",
        "records": {
            "restored_graph_rows": graph_rows,
            "q_origin_boundary_rows": boundary_rows,
            "excluded_loci": ["P=-4Q true Cauchy pole", "non-rational geometric points", "global gluing and normalization"],
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The restored rational P=4Q points and finite rational Q=O boundary zeros have no simultaneous formal derivative zero at levels 0,1,17. Every Q=O pencil member has a simple P=O boundary zero in the L(9O) section model.",
        "next_requirement": "Construct a global compact residual section and analyze the true pole and non-rational loci before promoting these local screens to a smooth curve, genus, Jacobian, relation, or ECDLP claim.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N609G boundary-singularity preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
