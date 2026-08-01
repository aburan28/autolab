#!/usr/bin/env sage -python
"""Structural replay for the N609K rational true-pole extension receipt."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    records = primary["records"]
    rows = records["rows"]
    formal = records["formal_replay"]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "complete_rational_grid_retained": records["row_count"] == 324
        and len(rows) == 324
        and len({row["q"] for row in rows}) == 108
        and {row["level"] for row in rows} == {0, 1, 17},
        "local_nonvanishing_data_retained": all(
            row["leading_determinant_nonzero"]
            and row["cauchy_coefficient_nonzero"]
            and row["regularized_true_pole_value_nonzero"]
            for row in rows
        ),
        "laurent_replay_retained": formal["q"] == "(1 : 51 : 1)"
        and formal["level"] == 0
        and formal["formula_matches_laurent_replay"]
        and formal["laurent_value_nonzero"],
        "local_trivialization_retained": records["local_parameter"] == "u=x(P)-x(-4Q)"
        and records["regularized_limit_formula"] == "-det(V_Q)*cauchy_coefficient*2*y(-4Q)",
        "global_extension_boundary_retained": "global cartier extension" in primary["claim_status"].lower()
        and "transition functions" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-true-pole-extension.n609k.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N609K structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
