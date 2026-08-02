#!/usr/bin/env sage -python
"""N610E: solve a bounded H018 branch from coefficient-wise replayed slope forms."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import time
from pathlib import Path

from sage.all import Matrix, PolynomialRing, vector

from poincare_determinantal_section_n609j import moving_basis
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def lead(series):
    return series[series.valuation()]


def interpolate(field, samples, degree):
    matrix = Matrix(field, [[slope**power for power in range(degree + 1)] for slope, _ in samples[: degree + 1]])
    coefficients = matrix.solve_right(vector(field, [value for _, value in samples[: degree + 1]]))
    polynomial_ring = PolynomialRing(field, "s")
    slope = polynomial_ring.gen()
    return sum(coefficients[power] * slope**power for power in range(degree + 1))


def replayed_polynomial(field, samples, max_degree):
    for degree in range(max_degree + 1):
        polynomial = interpolate(field, samples, degree)
        if all(polynomial(slope) == value for slope, value in samples[degree + 1:]):
            return polynomial, degree
    return None, None


def evaluate_polynomial(ring, polynomial, series):
    return sum(ring(polynomial[power]) * series**power for power in range(polynomial.degree() + 1))


def evaluate_form(ring, u, coefficients, series):
    return sum(evaluate_polynomial(ring, polynomial, series) * u**order for order, polynomial in enumerate(coefficients))


def differentiate_coefficients(coefficients):
    return [polynomial.derivative() for polynomial in coefficients]


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
    forms = {}
    reconstruction_rows = []

    for level in (0, 1, 17):
        moving = matrix.solve_right(xvalues - ring(level) * vector(ring, [1] * 9))
        samples = []
        for slope in slopes:
            parameter = ring(slope) * u
            point = (ring(formal.x(120)(parameter)), ring(formal.y(120)(parameter)))
            raw = sum(moving[index] * moving_basis(point, support)[index] for index in range(9))
            normal = -determinant * parameter**8 * raw
            samples.append((slope, (slope - cauchy_slope) * normal))
        coefficients = []
        degrees = []
        for order in range(args.max_order + 1):
            polynomial, degree = replayed_polynomial(
                field, [(slope, value[order]) for slope, value in samples], args.max_degree
            )
            if polynomial is None:
                raise RuntimeError(f"N610E could not replay coefficient u^{order}")
            coefficients.append(polynomial)
            degrees.append(degree)
        forms[level] = coefficients
        reconstruction_rows.append({"level": level, "minimal_replayed_degrees": degrees})

    roots = forms[0][0].roots(field)
    if len(roots) != 1:
        raise RuntimeError("N610E expected one bounded-jet boundary root")
    root, multiplicity = roots[0]
    branch = ring(root)
    derivative_coefficients = differentiate_coefficients(forms[0])
    iterations = []
    for iteration in range(8):
        residual = evaluate_form(ring, u, forms[0], branch)
        slope_derivative = evaluate_form(ring, u, derivative_coefficients, branch)
        if slope_derivative.valuation() != 0:
            raise RuntimeError("N610E bounded slope derivative is not a unit")
        iterations.append({"iteration": iteration, "residual_valuation": None if residual == 0 else int(residual.valuation())})
        if residual == 0 or residual.valuation() >= args.max_order:
            break
        branch -= residual / slope_derivative

    residual = evaluate_form(ring, u, forms[0], branch)
    differences = []
    for level in (1, 17):
        difference = evaluate_form(ring, u, forms[level], branch) - residual
        valuation = difference.valuation()
        differences.append({"level": level, "valuation": int(valuation), "leading_coefficient": str(difference[valuation]), "order_eight_unit": valuation == 8 and difference[8] != 0})

    correction = branch - ring(root)
    gates = {
        "constant_term_has_expected_simple_root": str(root) == "83" and multiplicity == 1,
        "level_zero_slope_derivative_is_a_unit": evaluate_form(ring, u, derivative_coefficients, ring(root)).valuation() == 0,
        "bounded_branch_residual_reaches_order_twenty_four": residual == 0 or residual.valuation() >= args.max_order,
        "both_differences_remain_order_eight_units_on_bounded_branch": all(row["order_eight_unit"] for row in differences),
    }
    output = {
        "schema": "ecdlp.h018.poincare-coefficientwise-branch.n610e.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / BOUNDED H018 COEFFICIENT-WISE FORMAL BRANCH / MODEL-BOUND / TOY-EVIDENCE / EXACT_TWO_VARIABLE_CHART_AND_GLOBAL_RESOLUTION_OPEN / NO_ECDLP_CLAIM",
        "parameters": {"max_u_order": args.max_order, "max_slope_degree": args.max_degree, "sample_slopes": [int(slope) for slope in slopes]},
        "records": {
            "cauchy_slope": str(cauchy_slope), "boundary_root": str(root), "boundary_root_multiplicity": int(multiplicity),
            "branch_correction_valuation": None if correction == 0 else int(correction.valuation()),
            "branch_correction_leading_coefficient": None if correction == 0 else str(correction[correction.valuation()]),
            "residual_valuation": None if residual == 0 else int(residual.valuation()),
            "reconstruction_rows": reconstruction_rows, "newton_iterations": iterations, "differences": differences,
            "elapsed_seconds": time.perf_counter() - started,
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "Each coefficient through the declared jet bound replays independently on held-out slopes. In that bounded reconstructed model, the level-zero branch through slope 83 has a unit slope derivative and the two independent differences remain order-eight units after the branch correction.",
        "next_requirement": "Derive an exact rational two-variable chart or extend the coefficient-wise validation with an explicit degree law. This bounded jet does not itself prove local colength, a global Cartier divisor, normalization, Jacobian, relations, rank, descent, cost, or an ECDLP advantage.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N610E coefficient-wise branch probe failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
