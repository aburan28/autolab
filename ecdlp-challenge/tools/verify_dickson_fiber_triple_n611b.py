#!/usr/bin/env sage -python
"""Independent structural replay for the N611B Dickson-fiber screen."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import random
import sys
from pathlib import Path

from sage.all import EllipticCurve, GF

sys.path.insert(0, str(Path(__file__).resolve().parent))
import xonly_pair_sum_membership_probe as xprobe


P = 101


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def seed_for(*parts):
    return int.from_bytes(hashlib.sha256("|".join(map(str, parts)).encode("ascii")).digest()[:8], "big")


def evaluate(x_value, degree, parameter, field):
    values = [field(2), field(x_value)]
    while len(values) <= degree:
        values.append(field(x_value) * values[-1] - field(parameter) * values[-2])
    return values[degree]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    field = GF(P)
    curve = EllipticCurve(field, [0, 0, 0, 1, 32])
    by_x = {}
    for point in curve.points():
        if point.is_zero():
            continue
        by_x.setdefault(int(point[0]), []).append((int(point[0]), int(point[1])))
    x_pairs = {x_value: tuple(sorted(points)) for x_value, points in by_x.items() if len(points) == 2}
    all_x = sorted(x_pairs)

    rows = []
    for row in primary["rows"]:
        if not row["admitted_fiber"]:
            continue
        degree, parameter, constant = (row[key] for key in ("degree", "parameter", "constant"))
        selected_x = [x_value for x_value in all_x if evaluate(x_value, degree, parameter, field) == field(constant)]
        points = [point for x_value in selected_x for point in x_pairs[x_value]]
        triple_count = sum(
            xprobe.add(xprobe.add(first, second, P, 1), third, P, 1) is None
            for first, second, third in itertools.combinations(points, 3)
        )
        controls = []
        for control in range(primary["parameters"]["controls"]):
            rng = random.Random(seed_for("n611b", degree, parameter, constant, control))
            sampled_x = rng.sample(all_x, len(selected_x))
            sampled_points = [point for x_value in sampled_x for point in x_pairs[x_value]]
            controls.append(sum(
                xprobe.add(xprobe.add(first, second, P, 1), third, P, 1) is None
                for first, second, third in itertools.combinations(sampled_points, 3)
            ))
        rows.append({
            "degree": degree,
            "parameter": parameter,
            "constant": constant,
            "factor_base_size_matches": len(points) == row["factor_base_size"],
            "membership_replays": all(evaluate(point[0], degree, parameter, field) == field(constant) for point in points),
            "triple_count_matches": triple_count == row["zero_sum_distinct_triples"],
            "control_counts_match": controls == row["random_control_counts"],
            "density_does_not_exceed_control_maximum": triple_count <= max(controls),
            "support_exceeds_pair_state": row["generic_remainder_profile"]["total_monomials"] >= row["unordered_x_pair_state"],
        })
    checks = {
        "primary_schema": primary["schema"] == "ecdlp.dickson-fiber-triple.n611b.v1",
        "primary_script_hash_matches": primary["script_sha256"] == sha256_file(Path(__file__).with_name("dickson_fiber_triple_probe_n611b.py")),
        "all_admitted_rows_replay": bool(rows) and all(all(value for key, value in row.items() if key not in {"degree", "parameter", "constant"}) for row in rows),
        "source_solver_remains_rejected": primary["promotion"]["admit_source_solver"] is False,
    }
    output = {
        "schema": "ecdlp.dickson-fiber-triple.n611b.independent-verifier.v1",
        "primary_sha256": sha256_file(args.primary),
        "verifier_sha256": sha256_file(Path(__file__)),
        "rows": rows,
        "checks": checks,
        "verification_pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verification_pass"]:
        raise RuntimeError("N611B independent verifier failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
