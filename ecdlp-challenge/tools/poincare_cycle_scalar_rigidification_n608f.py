#!/usr/bin/env sage -python
"""N608F: rule out an evaluation-independent scalar cycle strictification."""
import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

from sage.all import EllipticCurve, GF, gcd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import product_kummer_h018_poincare_scalar_cocycle_n606z as transport


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--sample-count", type=int, default=12)
    args = parser.parse_args()
    base = EllipticCurve(GF(103), [0, 0, 0, 1, 24])
    field = GF(103 ** 2, name="a")
    curve = base.base_extend(field)
    base_generator = next(point for point in base.points() if not point.is_zero())
    generator = curve(field(base_generator[0]), field(base_generator[1]))
    orbit = [index * (84 * generator) for index in range(109)]
    products = []
    for point in curve.points():
        if point.is_zero():
            continue
        try:
            factors = [transport.transport(orbit[index], generator, point + (108 - index) * generator) for index in range(109)]
        except ZeroDivisionError:
            continue
        if not all(factors):
            continue
        products.append((point, field.prod(factors)))
        if len(products) == args.sample_count:
            break
    unit_order = field.order() - 1
    values = [value for _point, value in products]
    records = {
        "field_degree": 2,
        "cycle_order": 109,
        "field_unit_order": int(unit_order),
        "common_regular_sample_count": len(products),
        "distinct_full_cycle_products": len(set(values)),
        "full_cycle_products": [str(value) for value in values],
        "sample_points": [str(point) for point, _value in products],
    }
    gates = {
        "enough_common_regular_samples": len(products) == args.sample_count,
        "all_cycle_products_nonzero": all(values),
        "unique_109th_root_per_sample": gcd(109, unit_order) == 1,
        "evaluation_independent_scalar_rigidification_rejected": len(set(values)) > 1,
    }
    output = {
        "schema": "ecdlp.h018.poincare-cycle-scalar-rigidification.n608f.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "NEGATIVE RESULT / EVALUATION-INDEPENDENT-SCALAR-RIGIDIFICATION-REJECTED / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The N607A one-point scalar strictification cannot be promoted to a common-domain scalar rigidification: the 109-edge product varies across sampled regular extension-field evaluations. A point-dependent or matrix-valued geometric bundle map remains untested.",
        "next_requirement": "Construct a point-dependent normalized Poincare frame or an explicit matrix-valued intertwiner to the N606U pushforward bundle; do not reuse a single scalar cycle gauge.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608F scalar rigidification preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
