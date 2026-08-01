#!/usr/bin/env sage -python
"""Independent replay for the REP-AUXLINE-017 ternary admission gate."""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path
from typing import Any

from sage.all import GF, PolynomialRing

sys.path.insert(0, str(Path(__file__).resolve().parent))
import auxiliary_isogeny_rational_chain_solver_probe as n477
import compositional_coset_factor_base_probe as n470


S4_CASES = ((101, 1, 32), (107, 1, 4), (127, 1, 7))


def s3(x: Any, y: Any, z: Any, a_curve: Any, b_curve: Any) -> Any:
    return (x - y) ** 2 * z**2 - 2 * ((x + y) * (x * y + a_curve) + 2 * b_curve) * z + (x * y - a_curve) ** 2 - 4 * b_curve * (x + y)


def power_sum(e1: Any, e2: Any, e3: Any, degree: int) -> Any:
    field = e1.parent()
    sequence = [field(3), e1, e1**2 - 2 * e2, e1**3 - 3 * e1 * e2 + 3 * e3]
    for index in range(4, degree + 1):
        sequence.append(e1 * sequence[index - 1] - e2 * sequence[index - 2] + e3 * sequence[index - 3])
    return sequence[degree]


def support_count(degree: int) -> int:
    return sum((degree - 3 * count_e3) // 2 + 1 for count_e3 in range(degree // 3 + 1))


def e3_degree_by_leading_reduction(s4: Any, field: Any) -> tuple[int, int]:
    ring = s4.parent()
    x1, x2, x3, z, _u = ring.gens()
    ring3 = PolynomialRing(field, names=("x1", "x2", "x3"), order="lex")
    y1, y2, y3 = ring3.gens()
    e1 = y1 + y2 + y3
    e2 = y1 * y2 + y1 * y3 + y2 * y3
    e3 = y1 * y2 * y3
    max_e3 = 0
    total_terms = 0
    for z_degree in range(int(s4.degree(z)) + 1):
        component = ring3.zero()
        for exponents, coefficient in s4.dict().items():
            if int(exponents[3]) == z_degree:
                component += ring3(coefficient) * y1**int(exponents[0]) * y2**int(exponents[1]) * y3**int(exponents[2])
        while component:
            alpha, beta, gamma = (int(value) for value in component.lm().degrees())
            if alpha < beta or beta < gamma:
                raise AssertionError("S4 z-coefficient is not symmetric")
            counts = alpha - beta, beta - gamma, gamma
            coefficient = component.lc()
            component -= coefficient * e1**counts[0] * e2**counts[1] * e3**counts[2]
            max_e3 = max(max_e3, counts[2])
            total_terms += 1
    return max_e3, total_terms


def verify(report_path: Path, samples: int) -> dict[str, Any]:
    report = json.loads(report_path.read_text())
    s4_checks = []
    for expected, case in zip(report["s4_elementary_rows"], S4_CASES, strict=True):
        p, a_curve, b_curve = case
        field = GF(p)
        ring = PolynomialRing(field, names=("x1", "x2", "x3", "z", "u"), order="lex")
        x1, x2, x3, z, u = ring.gens()
        s4 = s3(x1, x2, u, field(a_curve), field(b_curve)).resultant(s3(u, x3, z, field(a_curve), field(b_curve)), u)
        e3_degree, elementary_terms = e3_degree_by_leading_reduction(s4, field)
        valid = (
            int(s4.degree(z)) == expected["s4_z_degree"]
            and len(s4.dict()) == expected["s4_monomial_count"]
            and e3_degree == expected["s4_e3_degree"]
            and elementary_terms == expected["elementary_monomial_count"]
        )
        s4_checks.append({"p": p, "valid": valid, "e3_degree": e3_degree})
    fiber_checks = []
    for expected in report["fiber_rows"]:
        instance = n477.n473.load_instance(Path("target/release/gen_instance"), hex(int(expected["seed"])), int(expected["bits"]))
        p = int(instance["p"])
        degree = int(expected["degree"])
        factors = n470.factor_integer(p - 1)
        factor_base, metadata = n470.coset_factor_base(instance, degree, factors, int(expected["coset_index"]))
        field = GF(p)
        constant = field(int(metadata["coset_constant"]))
        rng = random.Random(int(expected["seed"]) ^ 0xA017)
        membership_valid = True
        for _ in range(samples):
            roots = [field(rng.randrange(p)) for _ in range(3)]
            e1 = sum(roots, field(0))
            e2 = roots[0] * roots[1] + roots[0] * roots[2] + roots[1] * roots[2]
            e3 = roots[0] * roots[1] * roots[2]
            equations = power_sum(e1, e2, e3, degree) == 3 * constant and power_sum(e1, e2, e3, 2 * degree) == 3 * constant**2 and e3**degree == constant**3
            if equations != all(root**degree == constant for root in roots):
                membership_valid = False
                break
        factor_size = len([point for point in factor_base if point is not None])
        valid = (
            membership_valid
            and support_count(degree) == expected["power_d_support"]
            and support_count(2 * degree) == expected["power_2d_support"]
            and factor_size == expected["factor_base_size"]
            and support_count(2 * degree) >= factor_size * (factor_size + 1) // 2
        )
        fiber_checks.append({"degree": degree, "samples": samples, "valid": valid})
    return {
        "report": str(report_path),
        "s4_checks": s4_checks,
        "fiber_checks": fiber_checks,
        "overall_passed": all(row["valid"] for row in s4_checks + fiber_checks),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("report", type=Path)
    parser.add_argument("--samples", type=int, default=257)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    output = verify(args.report, args.samples)
    if args.out:
        args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"overall_passed": output["overall_passed"]}, sort_keys=True))
    if not output["overall_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
