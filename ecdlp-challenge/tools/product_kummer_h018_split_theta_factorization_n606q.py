#!/usr/bin/env python3
"""Exact N606Q gate for a split type-(1,2) product-theta factorization."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import time
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class CM:
    a: int
    b: int = 0

    def __add__(self, other: "CM") -> "CM":
        return CM(self.a + other.a, self.b + other.b)

    def __mul__(self, other: "CM") -> "CM":
        return CM(
            self.a * other.a - 103 * self.b * other.b,
            self.a * other.b + self.b * other.a - 5 * self.b * other.b,
        )

    def conj(self) -> "CM":
        return CM(self.a - 5 * self.b, -self.b)

    def norm(self) -> int:
        value = self * self.conj()
        if value.b != 0:
            raise AssertionError("CM norm did not descend to Z")
        return value.a

    def text(self) -> str:
        if self.b == 0:
            return str(self.a)
        if self.a == 0:
            return "pi" if self.b == 1 else f"{self.b}*pi"
        sign = "+" if self.b > 0 else "-"
        coefficient = "" if abs(self.b) == 1 else str(abs(self.b)) + "*"
        return f"{self.a}{sign}{coefficient}pi"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return "unknown"


def gate(name: str, passed: bool, evidence: Any) -> dict[str, Any]:
    return {"name": name, "passed": bool(passed), "evidence": evidence}


def norm_lower_bound_for_nonintegral(pi_coefficient: int) -> int:
    """min_a N(a+b*pi) is ceil(387*b^2/4), exact for integral a."""
    if pi_coefficient == 0:
        return 0
    numerator = 387 * pi_coefficient * pi_coefficient
    return (numerator + 3) // 4


def elements_of_norm_at_most(bound: int) -> list[CM]:
    max_b = 0
    while norm_lower_bound_for_nonintegral(max_b + 1) <= bound:
        max_b += 1
    elements: list[CM] = []
    for b in range(-max_b, max_b + 1):
        # The displayed norm is positive definite, so this interval is exhaustive.
        for a in range(-bound - 6 * abs(b), bound + 6 * abs(b) + 1):
            element = CM(a, b)
            if element.norm() <= bound:
                elements.append(element)
    return sorted(set(elements), key=lambda item: (item.norm(), item.a, item.b))


def matrix_from_columns(first: tuple[CM, CM], second: tuple[CM, CM]) -> list[list[CM]]:
    a, c = first
    b, d = second
    return [
        [CM(a.norm() + 2 * c.norm()), a.conj() * b + CM(2) * (c.conj() * d)],
        [b.conj() * a + CM(2) * (d.conj() * c), CM(b.norm() + 2 * d.norm())],
    ]


def matrix_text(matrix: list[list[CM]]) -> list[list[str]]:
    return [[entry.text() for entry in row] for row in matrix]


def factorization_candidates() -> list[dict[str, Any]]:
    small_nine = elements_of_norm_at_most(9)
    small_eleven = elements_of_norm_at_most(11)
    first_columns = [
        (a, c)
        for a in small_nine
        for c in small_nine
        if a.norm() + 2 * c.norm() == 9
    ]
    second_columns = [
        (b, d)
        for b in small_eleven
        for d in small_eleven
        if b.norm() + 2 * d.norm() == 11
    ]
    target = CM(2, 1)
    matches: list[dict[str, Any]] = []
    off_diagonals: set[str] = set()
    for first in first_columns:
        for second in second_columns:
            matrix = matrix_from_columns(first, second)
            off_diagonals.add(matrix[0][1].text())
            if matrix[0][1] == target:
                matches.append(
                    {
                        "first_column": [entry.text() for entry in first],
                        "second_column": [entry.text() for entry in second],
                        "matrix": matrix_text(matrix),
                    }
                )
    return [
        {
            "first_column_count": len(first_columns),
            "second_column_count": len(second_columns),
            "candidate_pair_count": len(first_columns) * len(second_columns),
            "attainable_off_diagonals": sorted(off_diagonals),
            "matches": matches,
        }
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    elements_nine = elements_of_norm_at_most(9)
    elements_eleven = elements_of_norm_at_most(11)
    search = factorization_candidates()[0]
    split_control = matrix_from_columns((CM(1), CM(0)), (CM(0), CM(1)))
    nonintegral_lower_bounds = {
        str(value): norm_lower_bound_for_nonintegral(value) for value in (-2, -1, 1, 2)
    }
    diagonal_forces_integrality = all(item.b == 0 for item in elements_nine + elements_eleven)
    off_diagonal_integral = all("pi" not in value for value in search["attainable_off_diagonals"])
    contradiction = diagonal_forces_integrality and off_diagonal_integral and not search["matches"]
    gates = [
        gate(
            "nonintegral_cm_elements_have_norm_at_least_97",
            min(nonintegral_lower_bounds.values()) == 97,
            nonintegral_lower_bounds,
        ),
        gate(
            "h018_diagonal_budgets_force_integral_matrix_entries",
            diagonal_forces_integrality,
            {
                "norm_at_most_9": [item.text() for item in elements_nine],
                "norm_at_most_11": [item.text() for item in elements_eleven],
            },
        ),
        gate(
            "split_type_12_positive_control_reconstructs_diagonal",
            matrix_text(split_control) == [["1", "0"], ["0", "2"]],
            matrix_text(split_control),
        ),
        gate(
            "all_diagonal_budget_factorizations_have_integral_off_diagonal",
            off_diagonal_integral,
            search,
        ),
        gate(
            "no_split_type_12_factorization_matches_h018_off_diagonal_2_plus_pi",
            not search["matches"],
            search,
        ),
    ]
    passed = all(item["passed"] for item in gates)
    result = {
        "schema": "ecdlp.product-kummer.h018.split-theta-factorization.n606q.v1",
        "candidate": "N606Q split type-(1,2) product-theta factorization of H018",
        "claim_status": (
            "RESTRICTED THEOREM / NEGATIVE SPLIT_THETA_FACTORIZATION_GATE / INTRINSIC_NONSPLIT_THETA_DATA_REQUIRED / MODEL-BOUND / TOY-EVIDENCE / NOVELTY-UNVERIFIED / NO_ECDLP_CLAIM"
            if passed else "REVIEW_REQUIRED / N606Q SPLIT THETA FACTORIZATION GATE FAILED / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM"
        ),
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "command": ["python3", str(Path(__file__).resolve()), "--output", str(args.output)],
        "git_commit": git_commit(),
        "tool_sha256": sha256(Path(__file__).resolve()),
        "cm_order": "Z[pi], pi^2+5*pi+103=0",
        "target": {
            "h018": [["9", "2+pi"], ["-3-pi", "11"]],
            "split_polarization": [["1", "0"], ["0", "2"]],
            "factorization": "H018=A^*diag(1,2)A",
        },
        "norm_bound": {
            "formula": "N(a+b*pi)=(a-5b/2)^2+387*b^2/4",
            "nonintegral_minimum": 97,
            "small_elements_norm_at_most_9": [item.text() for item in elements_nine],
            "small_elements_norm_at_most_11": [item.text() for item in elements_eleven],
        },
        "bounded_factorization_search": search,
        "gates": gates,
        "split_theta_factorization_gate_pass": passed,
        "split_product_theta_evaluator_available": False,
        "still_false_gates": {
            "intrinsic_nonsplit_theta_sections_constructed": False,
            "rigidified_scalar_sections_s0_s1_constructed": False,
            "section_evaluator_on_refined_cover": False,
            "base_locus_computed": False,
            "fiber_degree_audit": False,
            "source_recovery_audit": False,
            "relation_rank": False,
            "target_descent": False,
            "subrho_cost": False,
            "algorithmic_success": False,
        },
        "strongest_valid_statement": (
            "Any split type-(1,2) factorization H018=A^*diag(1,2)A would force the four CM entries of A to have norm at most 11. Every such entry is integral because every nonintegral CM element has norm at least 97; hence the off-diagonal is integral, contradicting H018[0,1]=2+pi. A split product-theta coordinate change cannot supply the N606Q evaluator."
            if passed else "The N606Q split-factorization gate did not verify every exact norm and matrix condition."
        ),
        "next_requirement": (
            "Construct intrinsically non-split theta/Poincare sections rather than a product-theta pullback, with named scalar sections, rigidification anchors, and refined-cover evaluation."
        ),
        "elapsed_seconds": time.monotonic() - started,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not passed:
        failed = ", ".join(item["name"] for item in gates if not item["passed"])
        raise RuntimeError("N606Q split-theta factorization gate failed: " + failed)
    print(json.dumps({"status": result["claim_status"], "checks": f"{sum(item['passed'] for item in gates)}/{len(gates)}", "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
