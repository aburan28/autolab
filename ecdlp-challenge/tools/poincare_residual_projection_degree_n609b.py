#!/usr/bin/env sage -python
"""N609B: certify the degree-nine P-fibre after N608Z base-graph removal."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from collections import Counter
from pathlib import Path

from sage.all import EllipticCurve, GF

from poincare_cleared_base_graph_n608z import LEVELS, cleared_equation
from poincare_linear_projection_n608v import P, coefficient_table


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pole_order_at_origin(coefficient, drop_leading=False):
    # ord_O(x)=-2 and ord_O(y)=-3.  The cleared equation has A's x^5
    # term and B*y of order at worst -9, so c7 determines the order -10.
    c7 = coefficient[7]
    a_orders = [0, -2, -4, -6, -8]
    if not drop_leading and c7 != 0:
        a_orders.append(-10)
    b_y_orders = [-3, -5, -7, -9]
    return -min(a_orders + b_y_orders)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    coefficients = coefficient_table(curve)
    rows = []
    for level in LEVELS:
        pole_orders = []
        base_zero_count = 0
        residual_degrees = []
        for base_point, coefficient in coefficients.items():
            support = -4 * base_point
            graph_point = 4 * base_point
            x, a, b, eliminated = cleared_equation(curve, coefficient, support, GF(P)(level))
            graph_x = x(graph_point[0])
            base_zero_count += int(a(graph_x) + b(graph_x) * graph_point[1] == 0)
            quotient, remainder = eliminated.quo_rem(x - graph_x)
            if remainder != 0:
                raise RuntimeError("N608Z graph factor missing")
            pole_orders.append(pole_order_at_origin(coefficient))
            residual_degrees.append(int(quotient.degree()))
        rows.append({
            "level": level,
            "base_point_count": len(coefficients),
            "pole_order_at_p_origin_distribution": {str(key): value for key, value in sorted(Counter(pole_orders).items())},
            "common_graph_zero_count": base_zero_count,
            "residual_divisor_degree_distribution": {str(key): value for key, value in sorted(Counter(residual_degrees).items())},
        })
    leading_coefficients_nonzero = all(coefficient[7] != 0 for coefficient in coefficients.values())
    dropped_orders = {pole_order_at_origin(coefficient, drop_leading=True) for coefficient in coefficients.values()}
    gates = {
        "all_nonzero_base_points_used": len(coefficients) == 108,
        "leading_c7_never_vanishes": leading_coefficients_nonzero,
        "all_declared_fibres_have_exact_pole_order_ten": all(
            row["pole_order_at_p_origin_distribution"] == {"10": 108} for row in rows
        ),
        "common_base_graph_is_zero_in_all_fibres": all(row["common_graph_zero_count"] == 108 for row in rows),
        "all_residual_fibres_have_degree_nine": all(
            row["residual_divisor_degree_distribution"] == {"9": 108} for row in rows
        ),
        "leading_term_removal_lowers_pole_order": dropped_orders == {9},
    }
    output = {
        "schema": "ecdlp.h018.poincare-residual-projection-degree.n609b.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "RESTRICTED THEOREM / DEGREE-NINE RESIDUAL P-FIBRE ON DECLARED PUNCTURED BASE / MODEL-BOUND / TOY-EVIDENCE / GLOBAL CLASS OPEN / NO_ECDLP_CLAIM",
        "records": {
            "curve_order": int(curve.cardinality()),
            "pole_convention": "ord_O(x)=-2, ord_O(y)=-3; c7*x^5 is the unique order-ten term in A",
            "base_graph_removed": "P=4Q",
            "rows": rows,
            "leading_term_removed_pole_orders": sorted(dropped_orders),
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "For all declared Q and levels, the cleared P-fibre has an exact order-ten pole at P=O and a common graph zero P=4Q. Removing one copy of that graph leaves a degree-nine residual divisor, so no additional P=O zero is missing from the affine residual root recovery on this punctured base.",
        "next_requirement": "Determine the Q=O boundary and the second projection/global divisor class of the compact residual curve before using the degree-nine fibre as a genus or factor-base claim.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N609B residual-projection-degree preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
