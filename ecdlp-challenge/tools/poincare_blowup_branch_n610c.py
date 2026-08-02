#!/usr/bin/env sage -python
"""N610C: evaluate H018 pencil differences on the formal level-zero branch."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import Matrix, PolynomialRing, vector

from poincare_determinantal_section_n609j import moving_basis
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def lead(series):
    return series[series.valuation()]


def interpolate_series(field, ring, samples):
    vandermonde = Matrix(
        field, [[field(slope) ** power for power in range(10)] for slope, _ in samples[:10]]
    )
    inverse = vandermonde.inverse()
    coefficients = [
        sum(ring(inverse[power, row]) * samples[row][1] for row in range(10))
        for power in range(10)
    ]
    polynomial_ring = PolynomialRing(ring, "s")
    slope = polynomial_ring.gen()
    return sum(coefficients[power] * slope**power for power in range(10))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    target, cover, pi, phi = fixture()
    field, ring, u, support, matrix, xvalues = origin_data(target, cover, pi, phi)
    determinant = matrix.det()
    cauchy_slope = lead((-support[0] / support[1]) / u)
    formal = target.formal_group()
    slope_values = [field(value) for value in range(2, 32)]
    forms = {}
    reconstruction_rows = []

    for level in (0, 1, 17):
        moving = matrix.solve_right(xvalues - ring(level) * vector(ring, [1] * 9))
        samples = []
        for slope in slope_values:
            parameter = ring(slope) * u
            point = (ring(formal.x(120)(parameter)), ring(formal.y(120)(parameter)))
            raw = sum(moving[index] * moving_basis(point, support)[index] for index in range(9))
            normal = -determinant * parameter**8 * raw
            samples.append((slope, (slope - cauchy_slope) * normal))
        form = interpolate_series(field, ring, samples)
        heldout_matches = all(form(slope) == value for slope, value in samples[10:])
        forms[level] = form
        reconstruction_rows.append(
            {
                "level": level,
                "sample_count": len(samples),
                "train_count": 10,
                "heldout_count": len(samples) - 10,
                "heldout_series_reconstruction_match": heldout_matches,
            }
        )

    boundary_ring = PolynomialRing(field, "s0")
    boundary_slope = boundary_ring.gen()
    boundary = sum(field(forms[0][power][0]) * boundary_slope**power for power in range(10))
    roots = boundary.roots(field)
    gates = {
        "all_three_bivariate_reconstructions_reject_the_degree_nine_model": all(
            not row["heldout_series_reconstruction_match"] for row in reconstruction_rows
        ),
        "boundary_constant_term_still_has_the_expected_simple_slope": len(roots) == 1
        and str(roots[0][0]) == "83"
        and roots[0][1] == 1,
    }
    output = {
        "schema": "ecdlp.h018.poincare-blowup-branch.n610c.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "NEGATIVE RESULT / H018 DEGREE-NINE BIVARIATE SLOPE RECONSTRUCTION REJECTED / MODEL-BOUND / TOY-EVIDENCE / HIGHER-SERIES-STRUCTURE_OPEN / NO_ECDLP_CLAIM",
        "records": {
            "cauchy_slope": str(cauchy_slope),
            "boundary_roots": [{"root": str(root), "multiplicity": int(multiplicity)} for root, multiplicity in roots],
            "reconstruction_rows": reconstruction_rows,
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The degree-nine slope interpolation that succeeds for the exceptional constant term fails to replay the full formal series at every tested pencil level. Therefore it cannot be used to construct a formal branch or to infer a local base-scheme length.",
        "next_requirement": "Characterize the slope-degree growth coefficient by coefficient, or derive an exact rational two-variable chart. Keep the constant-term slope result separate from the rejected full-series reconstruction; global Cartier, normalization, Jacobian, relation law, rank, descent, cost, and ECDLP claims remain open.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N610C formal-branch probe failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
