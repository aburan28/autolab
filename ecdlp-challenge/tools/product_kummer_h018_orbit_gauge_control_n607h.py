#!/usr/bin/env sage -python
"""N607H: gauge-control the raw N607E finite-orbit transition sequence."""
import argparse
import json
import sys
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, identity_matrix

sys.path.insert(0, str(Path(__file__).resolve().parent))
import product_kummer_h018_poincare_scalar_cocycle_n606z as scalar


def basis(point, support):
    if support.is_zero():
        x, y = point[0], point[1]
        return [1, x, y, x*x, x*y, x**3, x*x*y, x**4, x**3*y]
    x, y = point[0], point[1]
    xr, yr = support[0], support[1]
    return [1, x, y, x*x, x*y, x**3, x*x*y, x**4, (y + yr) / (x - xr)]


def edge_matrix(curve, extension, q, step):
    q_next = q + step
    old_support, new_support = -4*q, -4*q_next
    points = [point for point in curve.points() if not point.is_zero()]
    old_rows, new_rows = [], []
    for point in points:
        try:
            if (not old_support.is_zero() and point in (old_support, -old_support)) or (not new_support.is_zero() and point in (new_support, -new_support)):
                continue
            factor = scalar.transport(q, step, point)
            if factor:
                old_rows.append([factor * value for value in basis(point + step, old_support)])
                new_rows.append(basis(point, new_support))
            if len(old_rows) == 9:
                matrix_new, matrix_old = Matrix(extension, new_rows), Matrix(extension, old_rows)
                if matrix_new.rank() == 9:
                    return matrix_new.solve_right(matrix_old)
        except ZeroDivisionError:
            continue
    raise RuntimeError("no invertible N607E interpolation edge")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    base = GF(103)
    curve0 = EllipticCurve(base, [0, 0, 0, 1, 24])
    extension = GF(103**2, name="a")
    curve = curve0.base_extend(extension)
    generator0 = next(point for point in curve0.points() if not point.is_zero())
    generator = curve(extension(generator0[0]), extension(generator0[1]))
    step = 84 * generator
    qs = [index * step for index in range(109)]
    edges = [edge_matrix(curve, extension, qs[index], step) for index in range(109)]
    gauges = [identity_matrix(extension, 9)]
    for edge in edges[:-1]:
        gauges.append(edge * gauges[-1])
    normalized_prefix = [gauges[index + 1].inverse() * edges[index] * gauges[index] for index in range(108)]
    monodromy = edges[-1] * gauges[-1]
    records = {
        "orbit_length": len(qs),
        "edge_count": len(edges),
        "all_edges_invertible": all(edge.rank() == 9 for edge in edges),
        "prefix_gauge_identity_count": sum(item == identity_matrix(extension, 9) for item in normalized_prefix),
        "closing_monodromy_rank": int(monodromy.rank()),
        "closing_monodromy_is_scalar": all(monodromy[i, j] == (monodromy[0, 0] if i == j else 0) for i in range(9) for j in range(9)),
    }
    gates = {
        "all_edges_invertible": records["all_edges_invertible"],
        "all_nonclosing_edges_gauge_to_identity": records["prefix_gauge_identity_count"] == 108,
        "closing_matrix_remains_full_rank": records["closing_monodromy_rank"] == 9,
    }
    output = {
        "schema": "ecdlp.product-kummer.h018.orbit-gauge-control.n607h.v1",
        "claim_status": "RESTRICTED THEOREM / FINITE_ORBIT_MONODROMY_NOT_BUNDLE_IDENTIFICATION / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The 108 nonclosing N607E edge matrices can be gauged to identity exactly. The remaining closing matrix is a holonomy of the chosen edge transport, not an explicit normalized-Poincare or pushforward-bundle identification.",
        "next_requirement": "Construct a geometric Poincare rigidification or an explicit bundle morphism on a genuine cover; finite-orbit matrix spectra alone are not an admissible substitute.",
    }
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N607H gauge control failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
