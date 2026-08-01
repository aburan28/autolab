#!/usr/bin/env sage -python
"""Exact ternary elementary-coset admission gate for REP-AUXLINE-017."""
from __future__ import annotations

import argparse
import json
import random
import sys
import time
from math import factorial
from pathlib import Path
from typing import Any

from sage.all import GF, PolynomialRing

sys.path.insert(0, str(Path(__file__).resolve().parent))
import auxiliary_isogeny_rational_chain_solver_probe as n477
import compositional_coset_factor_base_probe as n470


S4_CASES = ((101, 1, 32), (107, 1, 4), (127, 1, 7))
FIBERS = (("0x12345678", 20, 190, 1), ("0x1234567a", 22, 486, 0))


def semaev_s3(x1: Any, x2: Any, x3: Any, a_curve: Any, b_curve: Any) -> Any:
    return (
        (x1 - x2) ** 2 * x3**2
        - 2 * ((x1 + x2) * (x1 * x2 + a_curve) + 2 * b_curve) * x3
        + (x1 * x2 - a_curve) ** 2
        - 4 * b_curve * (x1 + x2)
    )


def support_triples(degree: int) -> list[tuple[int, int, int]]:
    return [
        (count_e1, count_e2, count_e3)
        for count_e3 in range(degree // 3 + 1)
        for count_e2 in range((degree - 3 * count_e3) // 2 + 1)
        for count_e1 in [degree - 2 * count_e2 - 3 * count_e3]
    ]


def girard_coefficient(degree: int, counts: tuple[int, int, int], field: Any) -> Any:
    count_e1, count_e2, count_e3 = counts
    length = count_e1 + count_e2 + count_e3
    numerator = degree * factorial(length - 1)
    denominator = factorial(count_e1) * factorial(count_e2) * factorial(count_e3)
    sign = -1 if (degree - length) % 2 else 1
    return field(sign * numerator) / field(denominator)


def power_sum_value(e1: Any, e2: Any, e3: Any, degree: int) -> Any:
    if degree == 0:
        return e1.parent()(3)
    if degree == 1:
        return e1
    previous3 = e1.parent()(3)
    previous2 = e1
    previous1 = e1 * e1 - 2 * e2
    if degree == 2:
        return previous1
    current = e1 * previous1 - e2 * previous2 + 3 * e3
    if degree == 3:
        return current
    for _ in range(4, degree + 1):
        previous3, previous2, previous1 = previous2, previous1, current
        current = e1 * previous1 - e2 * previous2 + e3 * previous3
    return current


def random_membership_check(field: Any, degree: int, constant: Any, samples: int, seed: int) -> dict[str, Any]:
    rng = random.Random(seed)
    passed = 0
    for _ in range(samples):
        roots = [field(rng.randrange(int(field.order()))) for _ in range(3)]
        e1 = sum(roots, field(0))
        e2 = roots[0] * roots[1] + roots[0] * roots[2] + roots[1] * roots[2]
        e3 = roots[0] * roots[1] * roots[2]
        equations = (
            power_sum_value(e1, e2, e3, degree) == 3 * constant
            and power_sum_value(e1, e2, e3, 2 * degree) == 3 * constant**2
            and e3**degree == constant**3
        )
        roots_in_coset = all(root**degree == constant for root in roots)
        if equations != roots_in_coset:
            raise AssertionError("ternary membership biconditional failed")
        passed += 1
    return {"samples": passed, "biconditional_passed": True}


def symmetric_reduce(coefficient: Any, ring3: Any) -> dict[tuple[int, int, int], Any]:
    x1, x2, x3 = ring3.gens()
    elementary1 = x1 + x2 + x3
    elementary2 = x1 * x2 + x1 * x3 + x2 * x3
    elementary3 = x1 * x2 * x3
    remaining = ring3(coefficient)
    result: dict[tuple[int, int, int], Any] = {}
    while remaining:
        exponents = tuple(int(value) for value in remaining.lm().degrees())
        if not exponents[0] >= exponents[1] >= exponents[2]:
            raise AssertionError(f"non-symmetric leading exponent {exponents}")
        counts = (exponents[0] - exponents[1], exponents[1] - exponents[2], exponents[2])
        leading = remaining.leading_coefficient()
        result[counts] = result.get(counts, ring3.base_ring()(0)) + leading
        remaining -= leading * elementary1**counts[0] * elementary2**counts[1] * elementary3**counts[2]
    return {key: value for key, value in result.items() if value}


def s4_elementary_degree(p: int, a_curve: int, b_curve: int) -> dict[str, Any]:
    field = GF(p)
    ring = PolynomialRing(field, names=("x1", "x2", "x3", "z", "u"), order="lex")
    x1, x2, x3, z, u = ring.gens()
    s4 = semaev_s3(x1, x2, u, field(a_curve), field(b_curve)).resultant(
        semaev_s3(u, x3, z, field(a_curve), field(b_curve)), u
    )
    ring3 = PolynomialRing(field, names=("x1", "x2", "x3"), order="lex")
    by_z: dict[int, Any] = {}
    for exponents, coefficient in s4.dict().items():
        x1_degree, x2_degree, x3_degree, z_degree, u_degree = (int(value) for value in exponents)
        if u_degree:
            raise AssertionError("resultant retained bridge variable")
        by_z[z_degree] = by_z.get(z_degree, ring3.zero()) + ring3(coefficient) * ring3.gen(0)**x1_degree * ring3.gen(1)**x2_degree * ring3.gen(2)**x3_degree
    elementary_by_z = {z_degree: symmetric_reduce(coefficient, ring3) for z_degree, coefficient in by_z.items()}
    reconstruction = ring.zero()
    elementary1 = x1 + x2 + x3
    elementary2 = x1 * x2 + x1 * x3 + x2 * x3
    elementary3 = x1 * x2 * x3
    for z_degree, terms in elementary_by_z.items():
        for counts, coefficient in terms.items():
            reconstruction += coefficient * elementary1**counts[0] * elementary2**counts[1] * elementary3**counts[2] * z**z_degree
    if reconstruction != s4:
        raise AssertionError("symmetric S4 reconstruction failed")
    e3_degree = max(counts[2] for terms in elementary_by_z.values() for counts in terms)
    return {
        "p": p,
        "a": a_curve,
        "b": b_curve,
        "s4_monomial_count": len(s4.dict()),
        "s4_z_degree": int(s4.degree(z)),
        "s4_e3_degree": e3_degree,
        "s4_symmetric_reconstruction_passed": True,
        "elementary_monomial_count": sum(len(terms) for terms in elementary_by_z.values()),
    }


def fiber_row(seed: str, bits: int, degree: int, coset_index: int, samples: int) -> dict[str, Any]:
    instance = n477.n473.load_instance(Path("target/release/gen_instance"), seed, bits)
    p = int(instance["p"])
    factors = n470.factor_integer(p - 1)
    factor_base, metadata = n470.coset_factor_base(instance, degree, factors, coset_index)
    field = GF(p)
    constant = field(int(metadata["coset_constant"]))
    support_d = support_triples(degree)
    support_2d = support_triples(2 * degree)
    if any(girard_coefficient(2 * degree, counts, field) == 0 for counts in support_2d):
        raise AssertionError("unexpected zero Girard coefficient")
    membership = random_membership_check(field, degree, constant, samples, int(instance["seed"]) ^ degree)
    factor_size = len([point for point in factor_base if point is not None])
    pair_state = factor_size * (factor_size + 1) // 2
    return {
        "seed": int(instance["seed"]),
        "bits": bits,
        "p": p,
        "degree": degree,
        "coset_index": coset_index,
        "coset_constant": int(constant),
        "factor_base_size": factor_size,
        "power_d_support": len(support_d),
        "power_2d_support": len(support_2d),
        "pair_state_dimension": pair_state,
        "power_2d_to_pair_state_ratio": len(support_2d) / pair_state,
        "power_2d_support_at_least_pair_state": len(support_2d) >= pair_state,
        "membership": membership,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=128)
    parser.add_argument("--out", type=Path, default=Path("notes/compositional_coset_ternary_preflight_rep_auxline017.json"))
    args = parser.parse_args()
    started = time.perf_counter()
    s4_rows = [s4_elementary_degree(*case) for case in S4_CASES]
    fibers = [fiber_row(*fiber, args.samples) for fiber in FIBERS]
    report = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "claim_status": "HYPOTHESIS / MODEL-BOUND / TOY-EVIDENCE",
        "contract": "notes/compositional_coset_ternary_preflight_contract.md",
        "exact_membership": [
            "P_d(e1,e2,e3)=3c",
            "P_2d(e1,e2,e3)=3c^2",
            "e3^d=c^3",
        ],
        "s4_elementary_rows": s4_rows,
        "fiber_rows": fibers,
        "elapsed_seconds": time.perf_counter() - started,
    }
    report["summary"] = {
        "all_membership_biconditionals_pass": all(row["membership"]["biconditional_passed"] for row in fibers),
        "all_s4_e3_degrees_are_four": all(row["s4_e3_degree"] == 4 for row in s4_rows),
        "all_power_2d_supports_at_least_pair_state": all(row["power_2d_support_at_least_pair_state"] for row in fibers),
        "advance_to_ternary_source_solver": False,
        "wall_clock_not_rho_comparable": True,
    }
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
