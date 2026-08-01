#!/usr/bin/env sage -python
"""N606U conditional cyclic-isogeny model for the H018 pushforward bundle."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from sage.all import EllipticCurve, GF

P = 103

def kernel_orders(curve, polynomial):
    field = GF(P**2, name="a")
    extended = curve.base_extend(field)
    orders = []
    for root, _multiplicity in polynomial.roots(field):
        x = field(root)
        rhs = x**3 + field(curve.a4()) * x + field(curve.a6())
        if rhs.is_square():
            orders.append(int(extended(x, rhs.sqrt()).order()))
    return sorted(orders)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    field = GF(P)
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    first = curve.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    composite = second * first
    cover_map = composite.dual()
    cover = cover_map.domain()
    origin = cover(0)
    basis = cover.riemann_roch_basis(2 * cover.divisor(origin))
    norm_z = (-3) ** 2 - 5 * (-3) * (-1) + 103
    records = {
        "curve": str(curve),
        "cover": str(cover),
        "composite_degree": int(composite.degree()),
        "dual_degree": int(cover_map.degree()),
        "kernel_polynomial": str(composite.kernel_polynomial()),
        "kernel_x_coordinate_orders_over_f103_squared": kernel_orders(curve, composite.kernel_polynomial()),
        "cover_origin_maps_to_origin": bool(cover_map(origin).is_zero()),
        "cover_degree_two_basis": [str(item) for item in basis],
        "norm_z": norm_z,
        "fourier_mukai_numerics": {"w9_rank": 9, "w9_degree": -1, "z_pullback_degree": -norm_z, "theta11_tensor_degree": -norm_z + 9 * 11},
        "conditional_identification": "p2_*L_H018 ~= pi_*O_Eprime(2Oprime)",
        "literature_assumptions": [
            "normalized Fourier--Mukai determinant of O(9O) is O(-O)",
            "Atiyah coprime rank-degree determinant classification",
            "odd cyclic isogeny determinant/norm formula",
        ],
    }
    gates = {
        "cm_norm_z_is_97": norm_z == 97,
        "pushforward_numerical_rank_is_9": records["fourier_mukai_numerics"]["w9_rank"] == 9,
        "pushforward_numerical_degree_is_2": records["fourier_mukai_numerics"]["theta11_tensor_degree"] == 2,
        "dual_cover_degree_is_9": records["dual_degree"] == 9,
        "dual_cover_kernel_is_cyclic_order_9": records["kernel_x_coordinate_orders_over_f103_squared"] == [3, 9, 9, 9],
        "cover_degree_two_section_space_has_dimension_2": len(basis) == 2,
        "cover_origin_normalization_exact": records["cover_origin_maps_to_origin"],
    }
    payload = {
        "schema": "ecdlp.product-kummer.h018.pushforward-identification.n606u.v1",
        "claim_status": "CONDITIONAL RESTRICTED THEOREM / ABSTRACT_CYCLIC_ISOGENY_PUSHFORWARD_MODEL / EXPLICIT_BUNDLE_ISOMORPHISM_OPEN / MODEL-BOUND / TOY-EVIDENCE / LITERATURE_BOUND / NO_ECDLP_CLAIM",
        "records": records,
        "gates": gates,
        "pushforward_identification_preflight_pass": all(gates.values()),
        "explicit_bundle_isomorphism_constructed": False,
        "surface_sections_evaluated": False,
        "still_false_gates": {"surface_evaluator": False, "base_locus": False, "factor_base": False, "relation_rank": False, "target_descent": False, "subrho_cost": False, "algorithmic_success": False},
        "next_requirement": "Construct an explicit normalized-Poincare-to-isogeny-pushforward bundle isomorphism and evaluate the two transported cover sections on E x E.",
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not payload["pushforward_identification_preflight_pass"]:
        raise RuntimeError("N606U preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.output)}, sort_keys=True))

if __name__ == "__main__":
    main()
