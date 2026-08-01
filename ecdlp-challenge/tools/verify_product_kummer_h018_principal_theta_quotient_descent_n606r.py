#!/usr/bin/env python3
"""Independent finite-kernel replay for the N606R quotient descent gate."""

from __future__ import annotations

import argparse
import json
import time
from datetime import UTC, datetime
from pathlib import Path


def add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return left[0] ^ right[0], left[1] ^ right[1]


def multiply_pi(value: tuple[int, int]) -> tuple[int, int]:
    a, b = value
    return b, a ^ b


def frobenius(value: tuple[tuple[int, int], tuple[int, int]]) -> tuple[tuple[int, int], tuple[int, int]]:
    return multiply_pi(value[0]), multiply_pi(value[1])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    origin = ((0, 0), (0, 0))
    g = ((0, 1), (1, 0))
    h = ((1, 1), (0, 1))
    gh = (add(g[0], h[0]), add(g[1], h[1]))
    lines = {"G_g": frozenset((origin, g)), "G_h": frozenset((origin, h)), "G_gh": frozenset((origin, gh))}
    by_line = {values: name for name, values in lines.items()}
    cycle = {name: by_line[frozenset(frobenius(value) for value in values)] for name, values in lines.items()}
    stable = [name for name, values in lines.items() if frozenset(frobenius(value) for value in values) == values]
    cubic = {
        name: all(frobenius(frobenius(frobenius(value))) == value for value in values)
        for name, values in lines.items()
    }
    result = {
        "schema": "ecdlp.product-kummer.h018.principal-theta-quotient-descent.n606r.verify.v1",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "primary": {"path": str(args.primary), "claim_status": primary["claim_status"]},
        "recomputed": {
            "kernel_vectors": {"g": [list(g[0]), list(g[1])], "h": [list(h[0]), list(h[1])], "g_plus_h": [list(gh[0]), list(gh[1])]},
            "frobenius_line_cycle": cycle,
            "base_field_stable_maximal_lines": stable,
            "cubic_stability": cubic,
        },
    }
    result["verified"] = (
        primary["principal_theta_quotient_descent_gate_pass"] is True
        and cycle == {"G_g": "G_h", "G_h": "G_gh", "G_gh": "G_g"}
        and stable == []
        and all(cubic.values())
        and primary["base_field_principal_theta_quotient_available"] is False
        and primary["cubic_orbit_principal_quotient_route_available"] is True
    )
    result["strongest_valid_statement"] = (
        "Independent finite-field replay confirms that the three H018 maximal isotropic lines form a single Frobenius three-cycle and only stabilize after cubic extension."
        if result["verified"]
        else "The N606R finite-kernel replay did not verify every required quotient-descent condition."
    )
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not result["verified"]:
        raise RuntimeError("N606R independent principal-theta quotient descent replay failed")
    print(json.dumps({"verified": True, "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
