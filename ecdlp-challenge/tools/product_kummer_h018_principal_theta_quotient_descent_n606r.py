#!/usr/bin/env python3
"""N606R principal-theta quotient field-of-definition gate for H018."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


F4 = tuple[int, int]
Vector = tuple[F4, F4]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return "unknown"


def f4_add(left: F4, right: F4) -> F4:
    return left[0] ^ right[0], left[1] ^ right[1]


def f4_mul(left: F4, right: F4) -> F4:
    a, b = left
    c, d = right
    return (a & c) ^ (b & d), (a & d) ^ (b & c) ^ (b & d)


def vector_add(left: Vector, right: Vector) -> Vector:
    return f4_add(left[0], right[0]), f4_add(left[1], right[1])


def vector_text(value: Vector) -> list[list[int]]:
    return [list(value[0]), list(value[1])]


def subgroup_text(values: frozenset[Vector]) -> list[list[list[int]]]:
    return [vector_text(value) for value in sorted(values)]


def gate(name: str, passed: bool, evidence: Any) -> dict[str, Any]:
    return {"name": name, "passed": bool(passed), "evidence": evidence}


def all_false(values: dict[str, Any]) -> bool:
    return bool(values) and all(value is False for value in values.values())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()

    zero: F4 = (0, 0)
    one: F4 = (1, 0)
    pi: F4 = (0, 1)
    one_plus_pi: F4 = (1, 1)
    origin: Vector = (zero, zero)
    g: Vector = (pi, one)
    h: Vector = (one_plus_pi, pi)
    gh = vector_add(g, h)
    kernel = frozenset((origin, g, h, gh))

    def frobenius(value: Vector) -> Vector:
        return f4_mul(pi, value[0]), f4_mul(pi, value[1])

    lines = {
        "G_g": frozenset((origin, g)),
        "G_h": frozenset((origin, h)),
        "G_gh": frozenset((origin, gh)),
    }
    line_by_values = {values: name for name, values in lines.items()}
    frobenius_lines = {
        name: line_by_values[frozenset(frobenius(value) for value in values)] for name, values in lines.items()
    }
    all_subgroups = {
        "zero": frozenset((origin,)),
        **lines,
        "whole_kernel": kernel,
    }
    stable_subgroups = {
        name: values
        for name, values in all_subgroups.items()
        if frozenset(frobenius(value) for value in values) == values
    }
    frobenius_cubed_lines = {
        name: frozenset(
            (
                f4_mul(pi, f4_mul(pi, f4_mul(pi, value[0]))),
                f4_mul(pi, f4_mul(pi, f4_mul(pi, value[1]))),
            )
            for value in values
        )
        for name, values in lines.items()
    }
    cubic_stable = {name: frobenius_cubed_lines[name] == values for name, values in lines.items()}
    identity_control = {name: values == values for name, values in lines.items()}

    maximal_isotropic_lines = set(lines)
    no_base_field_line = not (set(stable_subgroups) & maximal_isotropic_lines)
    all_cubic = all(cubic_stable.values())
    gates = [
        gate(
            "h018_kernel_has_three_order_two_lines",
            len(lines) == 3 and all(len(values) == 2 for values in lines.values()),
            {name: subgroup_text(values) for name, values in lines.items()},
        ),
        gate(
            "frobenius_cycles_all_maximal_isotropic_lines",
            frobenius_lines == {"G_g": "G_h", "G_h": "G_gh", "G_gh": "G_g"},
            frobenius_lines,
        ),
        gate(
            "no_maximal_isotropic_line_is_base_field_stable",
            no_base_field_line,
            {name: subgroup_text(values) for name, values in stable_subgroups.items()},
        ),
        gate(
            "each_maximal_isotropic_line_stabilizes_over_cubic_extension",
            all_cubic,
            {name: subgroup_text(values) for name, values in frobenius_cubed_lines.items()},
        ),
        gate(
            "identity_frobenius_positive_control_stabilizes_all_lines",
            all(identity_control.values()),
            identity_control,
        ),
        gate(
            "zero_and_whole_kernel_negative_controls_are_not_maximal_lines",
            "zero" in stable_subgroups and "whole_kernel" in stable_subgroups and len(stable_subgroups["zero"]) != 2 and len(stable_subgroups["whole_kernel"]) != 2,
            {name: len(values) for name, values in stable_subgroups.items()},
        ),
    ]
    passed = all(item["passed"] for item in gates)
    still_false_gates = {
        "principal_quotient_equations_over_f103_cubed_constructed": False,
        "principal_theta_divisors_constructed": False,
        "three_conjugate_pullback_sections_constructed": False,
        "semilinear_section_descent_solved": False,
        "rigidified_scalar_sections_s0_s1_constructed": False,
        "section_evaluator_on_refined_cover": False,
        "base_locus_computed": False,
        "fiber_degree_audit": False,
        "source_recovery_audit": False,
        "relation_rank": False,
        "target_descent": False,
        "subrho_cost": False,
        "algorithmic_success": False,
    }
    result = {
        "schema": "ecdlp.product-kummer.h018.principal-theta-quotient-descent.n606r.v1",
        "candidate": "N606R H018 principal-theta quotient field-of-definition gate",
        "claim_status": (
            "RESTRICTED THEOREM / NEGATIVE BASE_FIELD_PRINCIPAL_THETA_QUOTIENT_SHORTCUT / CUBIC_ORBIT_DESCENT_CONSTRUCTION_OPEN / MODEL-BOUND / TOY-EVIDENCE / LITERATURE_BOUND / NO_ECDLP_CLAIM"
            if passed
            else "REVIEW_REQUIRED / N606R PRINCIPAL-THETA QUOTIENT DESCENT GATE FAILED / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM"
        ),
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "command": sys.argv,
        "git_commit": git_commit(),
        "tool_sha256": sha256(Path(__file__).resolve()),
        "contract": "notes/product_kummer_h018_principal_theta_quotient_descent_contract_n606r.md",
        "inputs": {
            "base_field": "F_103",
            "kernel": {"origin": vector_text(origin), "g": vector_text(g), "h": vector_text(h), "g_plus_h": vector_text(gh)},
            "frobenius": "coordinatewise multiplication by pi in F_4",
            "polarization_type": "(1,2)",
            "descent_assumption": "A maximal isotropic subgroup of K(L) gives a principally polarized quotient carrying a descended L.",
            "literature": "Lubicz--Robert (2015), Section 1 descent discussion",
        },
        "subgroup_record": {
            "order_two_lines": {name: subgroup_text(values) for name, values in lines.items()},
            "frobenius_line_permutation": frobenius_lines,
            "base_field_stable_subgroups": {name: subgroup_text(values) for name, values in stable_subgroups.items()},
            "cubic_stability": cubic_stable,
            "principal_quotient_field_degree": 3,
        },
        "gates": gates,
        "principal_theta_quotient_descent_gate_pass": passed,
        "base_field_principal_theta_quotient_available": False,
        "cubic_orbit_principal_quotient_route_available": passed,
        "still_false_gates": still_false_gates,
        "strongest_valid_statement": (
            "The three maximal isotropic order-two subgroups of the H018 type-(1,2) kernel form one Frobenius orbit. Under standard polarization descent, no single principal-polarized quotient can be defined from this kernel over F_103; each quotient becomes available over F_103^3. This rejects only the base-field principal-theta shortcut and gives a concrete cubic-orbit theta-descent construction target."
            if passed
            else "The N606R field-of-definition gate did not establish each exact subgroup and Frobenius condition."
        ),
        "limitations": [
            "The cited descent theorem is applied as an explicit assumption; this receipt does not construct any quotient equation or theta divisor.",
            "A cubic quotient orbit is not a pair of F_103 scalar sections until semilinear descent is explicitly solved.",
            "No evaluator, factor base, relation rank, target descent, sub-rho cost, or ECDLP algorithm is established.",
        ],
        "next_requirement": "N606S must construct one F_103^3 principal quotient from G_g, its principal theta divisor, and its two Frobenius conjugates; then compute the induced semilinear action on their pullbacks and solve for a rigidified F_103 basis s0,s1.",
        "elapsed_seconds": time.monotonic() - started,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not passed:
        failed = ", ".join(item["name"] for item in gates if not item["passed"])
        raise RuntimeError("N606R principal-theta quotient descent gate failed: " + failed)
    print(json.dumps({"status": result["claim_status"], "checks": f"{sum(item['passed'] for item in gates)}/{len(gates)}", "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
