#!/usr/bin/env python3
"""N608H: reproducible fresh-seed cost panel for the local generic solver."""
from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import subprocess
from pathlib import Path

DEFAULT_ORDER = 616_882_790_773
SEEDS = (0x1234_5679, 0x1234_567A, 0x1234_567B, 0x1234_567C)
TOKEN_SEEDS = (0xA608_0001, 0xA608_0002, 0xA608_0003, 0xA608_0004)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--results", type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    solver = repo / "src/solver/mod.rs"
    results_path = args.results or repo / "results.tsv"
    target_notes = [f"N608H clean panel row {index}" for index in range(1, 5)]
    found = {}
    for line in reversed(results_path.read_text(encoding="ascii").splitlines()):
        fields = line.split("\t", 7)
        if len(fields) != 8 or fields[7] not in target_notes or fields[7] in found:
            continue
        found[fields[7]] = fields
    if set(found) != set(target_notes):
        raise RuntimeError("N608H clean-panel oracle rows are incomplete")
    rows = []
    for index, (seed, token_seed) in enumerate(zip(SEEDS, TOKEN_SEEDS), start=1):
        fields = found[f"N608H clean panel row {index}"]
        public = json.loads(
            subprocess.check_output(
                ["./target/release/gen_instance", hex(seed), "40"],
                cwd=repo,
                text=True,
            )
        )
        rows.append(
            {
                "instance_seed": hex(seed),
                "token_seed": hex(token_seed),
                "group_ops": int(fields[2]),
                "rho_reference": int(fields[4]),
                "ratio_to_rho": float(fields[5]),
                "order": public["n"],
                "oracle_status": fields[6],
                "result_commit": fields[1],
            }
        )
    ratios = [row["ratio_to_rho"] for row in rows]
    total_ops = sum(row["group_ops"] for row in rows)
    total_rho = sum(row["rho_reference"] for row in rows)
    records = {
        "bits": 40,
        "row_count": len(rows),
        "rows": rows,
        "solver_sha256": sha256_file(solver),
        "aggregate_group_ops": total_ops,
        "aggregate_rho_reference": total_rho,
        "aggregate_ratio_to_rho": total_ops / total_rho,
        "ratio_min": min(ratios),
        "ratio_median": statistics.median(ratios),
        "ratio_max": max(ratios),
    }
    gates = {
        "all_rows_correct": all(row["oracle_status"] == "OK" for row in rows),
        "all_rows_are_nondefault_order": all(row["order"] != DEFAULT_ORDER for row in rows),
        "all_rows_use_distinct_seed_pairs": len({(row["instance_seed"], row["token_seed"]) for row in rows}) == len(rows),
        "aggregate_beats_rho_reference": records["aggregate_ratio_to_rho"] < 1.0,
    }
    output = {
        "schema": "ecdlp.generic-radding.fresh-seed-benchmark.n608h.v1",
        "claim_status": "OBSERVATION / FRESH-SEED GENERIC COST MEASUREMENT / MODEL-BOUND / TOY-EVIDENCE / NO SUB-SQUARE-ROOT CLAIM / NO ECDLP BREAK CLAIM",
        "records": records,
        "gates": gates,
        "preflight_pass": all(value for name, value in gates.items() if name != "aggregate_beats_rho_reference"),
        "hypothesis_supported": gates["aggregate_beats_rho_reference"],
        "strongest_valid_statement": "This fixed fresh-seed panel measures only the current generic solver's finite cost distribution. Passing would not evade the generic-group square-root scale; failure preserves the current r-adding construction as a variance-sensitive baseline.",
        "next_requirement": "Only modify the solver after a candidate is compared on this fixed panel and a larger independent holdout panel with the same oracle accounting.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    print(json.dumps({"rows": len(rows), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
