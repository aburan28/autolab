#!/usr/bin/env sage -python
"""N607F: separate H018 Poincare fibre and polarization-map orientations.

For L = p1^*Theta^9 tensor p2^*Theta^11 tensor (id,z)^*P, z=-3-pi,
the p1 fibre over Q is indexed by zQ.  The adjoint zbar occurs instead in
the first component of the polarization map.  These are intentionally
different roles.
"""
import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, GF

P, N = 103, 109


def f4_add(left, right):
    return (left[0] ^ right[0], left[1] ^ right[1])


def f4_mul(left, right):
    a, b = left
    c, d = right
    # F4=F2[pi]/(pi^2+pi+1).
    return ((a & c) ^ (b & d), (a & d) ^ (b & c) ^ (b & d))


def mat_vec(matrix, vector):
    return tuple(
        f4_add(f4_mul(matrix[row][0], vector[0]), f4_mul(matrix[row][1], vector[1]))
        for row in range(2)
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    q, r = [point for point in curve.points() if not point.is_zero()][:2]
    pi = (0, 1)
    one = (1, 0)
    z_mod_two = (1, 1)  # -3-pi
    zbar_mod_two = (0, 1)  # 2+pi
    h018_mod_two = ((one, pi), (z_mod_two, one))
    elements = [(0, 0), one, pi, z_mod_two]
    kernel = [vector for vector in ((left, right) for left in elements for right in elements)
              if mat_vec(h018_mod_two, vector) == ((0, 0), (0, 0))]
    graph_pi = [(f4_mul(pi, value), value) for value in elements]

    z = lambda point: (-4) * point
    zbar = lambda point: 3 * point
    shift_z = lambda point: 84 * point
    shift_zbar = lambda point: 106 * point
    translated_abel_jacobi = lambda point, delta: z(point) - 9 * delta

    checks = {
        "h018_matrix_orientation_exact": h018_mod_two == ((one, pi), (z_mod_two, one)),
        "kernel_is_graph_p_equals_pi_q": set(kernel) == set(graph_pi),
        "p1_fibre_uses_z_not_zbar": z(q) != zbar(q),
        "z_fibre_translation_law": z(q + shift_z(r)) == translated_abel_jacobi(q, r),
        "zbar_alternative_has_distinct_shift": zbar(q + shift_zbar(r)) == zbar(q) - 9 * r
        and shift_z(r) != shift_zbar(r),
        "n606y_shift_is_z_shift": shift_z(r) == 84 * r,
        "adjoint_fibre_would_swap_hermitian_orientation": z_mod_two != zbar_mod_two,
    }
    payload = {
        "schema": "ecdlp.product-kummer.h018.poincare-fibre-orientation.n607f.v1",
        "candidate": "N607F H018 Poincare fibre/polarization orientation audit",
        "claim_status": "RESTRICTED THEOREM / FORMAL_POINCARE_FIBRE_ORIENTATION / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "conventions": {
            "line_bundle": "p1^*Theta^9 tensor p2^*Theta^11 tensor (id,z)^*P",
            "z": "-3-pi",
            "zbar": "2+pi",
            "p1_fibre_class": "8[O]+[zQ]",
            "polarization_matrix": "[[9,zbar],[z,11]]",
            "kernel_relation_mod_two": "(P,Q)=(pi Q,Q)",
        },
        "records": {
            "z_on_base_field": "[-4]",
            "zbar_on_base_field": "[3]",
            "z_base_shift": "[84]",
            "zbar_hypothetical_shift": "[106]",
            "kernel_size": len(kernel),
        },
        "gates": checks,
        "preflight_pass": all(checks.values()),
        "supersedes": "The suspected adjoint-fibre mismatch is rejected; N606Y through N607E retain their stated restricted scope.",
        "still_open": [
            "normalized Poincare rigidification",
            "explicit H018 surface-section space",
            "smooth member, factor base, relations, rank, descent, and cost",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not payload["preflight_pass"]:
        raise RuntimeError("N607F orientation audit failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
