#!/usr/bin/env sage -python
"""Structural replay for the N609L rational P-origin extension receipt."""
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
    replay = records["local_parameter_replay"]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "complete_rational_grid_retained": records["row_count"] == 324
        and len(rows) == 324
        and len({row["q"] for row in rows}) == 108
        and {row["level"] for row in rows} == {0, 1, 17},
        "local_nonvanishing_data_retained": all(
            row["leading_determinant_nonzero"]
            and row["x4_coefficient_nonzero"]
            and row["regularized_p_origin_value_nonzero"]
            for row in rows
        ),
        "local_parameter_replay_retained": replay["q"] == "(1 : 51 : 1)"
        and replay["level"] == 0
        and replay["formula_matches_local_parameter_replay"]
        and replay["local_parameter_value_nonzero"]
        and replay["only_x4_has_nonzero_normalized_constant"],
        "trivialization_data_retained": records["local_parameter"] == "v=-x(P)/y(P)"
        and records["regularized_limit_formula"] == "-det(V_Q)*x4_coefficient"
        and records["zeroed_x4_coefficient_mutation_value"] == 0,
        "global_extension_boundary_retained": "global cartier extension" in primary["claim_status"].lower()
        and "cech transition functions" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-p-origin-extension.n609l.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N609L structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
