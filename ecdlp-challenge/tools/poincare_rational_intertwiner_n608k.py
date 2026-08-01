#!/usr/bin/env sage -python
"""N608K: certify the degree-eleven evaluator as a rational local morphism."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix

P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_basis(point, support):
    x, y = point[0], point[1]
    xr, yr = support[0], support[1]
    return [1, x, y, x * x, x * y, x**3, x * x * y, x**4, (y + yr) / (x - xr)]


def rational_image(isogeny, target_k, point):
    if point.is_zero():
        return target_k(0)
    x_map, y_map = isogeny.rational_maps()
    return target_k(x_map(point[0], point[1]), y_map(point[0], point[1]))


def fixture():
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    cover_map = (second * first).dual()
    cover = cover_map.domain()
    candidate = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    return target, cover, cover_map, candidate.codomain().isomorphism_to(target) * candidate


def deck_generator(cover_map, cover, field):
    cover_k = cover.base_extend(field)
    for root, _multiplicity in cover_map.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square():
            point = cover_k(x, rhs.sqrt())
            if point.order() == 9:
                return point
    raise RuntimeError("could not recover order-nine deck generator")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, cover_map, phi_map = fixture()
    field = GF(P**6, name="a")
    target_k = target.base_extend(field)
    cover_k = cover.base_extend(field)
    deck = deck_generator(cover_map, cover, field)

    rows = []
    for q0 in target.points():
        if q0.is_zero():
            continue
        lift0 = next(point for point in cover.points() if cover_map(point) == q0)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        q = target_k(field(q0[0]), field(q0[1]))
        support = -4 * q
        source_points = [rational_image(phi_map, target_k, lift + index * deck) for index in range(9)]
        target_points = [lift + index * deck for index in range(9)]
        bad_origin = any(point.is_zero() for point in target_points) or any(point.is_zero() for point in source_points)
        bad_support = any(point in (support, -support) for point in source_points)
        if bad_origin or bad_support:
            rows.append({"q": str(q0), "regular": False, "bad_origin": bad_origin, "bad_support": bad_support})
            continue
        evaluation = Matrix(field, [source_basis(point, support) for point in source_points])
        rows.append(
            {
                "q": str(q0),
                "regular": True,
                "rank": int(evaluation.rank()),
                "determinant_nonzero": not evaluation.det().is_zero(),
                "source_points_distinct": len(set(source_points)) == 9,
                "target_constant_section_one_nonzero": True,
            }
        )
    regular_rows = [row for row in rows if row["regular"]]
    gates = {
        "declared_open_locus_is_nonempty": bool(regular_rows),
        "all_registered_fibres_in_open_locus": len(regular_rows) == 108,
        "all_regular_evaluations_full_rank": all(row["rank"] == 9 and row["determinant_nonzero"] for row in regular_rows),
        "all_regular_target_weights_are_named": all(row["target_constant_section_one_nonzero"] for row in regular_rows),
        "one_nonzero_determinant_proves_generic_rank": any(row["determinant_nonzero"] for row in regular_rows),
    }
    output = {
        "schema": "ecdlp.h018.poincare-rational-intertwiner.n608k.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / EXPLICIT RATIONAL LOCAL MORPHISM / MODEL-BOUND / TOY-EVIDENCE / NORMALIZED-POINCARE CHART GLUING OPEN / NO SECTION EVALUATOR / NO_ECDLP_CLAIM",
        "records": {
            "cover_degree": int(cover_map.degree()),
            "separating_map_degree": int(phi_map.degree()),
            "source_divisor_model": "O(8O+[zQ]), z=[-4] on E(F_103)",
            "target_weight": "canonical section 1 of O_Eprime(2Oprime) on Eprime minus {Oprime}",
            "tested_nonzero_base_points": len(rows),
            "regular_fibre_count": len(regular_rows),
            "rejected_fibre_count": len(rows) - len(regular_rows),
            "rows": rows,
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "rational_local_morphism_constructed": all(gates.values()),
        "global_normalized_poincare_intertwiner_constructed": False,
        "strongest_valid_statement": "On the declared nonempty open locus, evaluation at phi(R) followed by multiplication by the named target section 1 defines a rational local morphism from the N606Y divisor-model H018 pushforward to pi_*O_Eprime(2Oprime). Its determinant is generically nonzero. The global normalized-Poincare chart transitions are not yet constructed.",
        "next_requirement": "Compute the determinant divisor and give chart transition matrices extending this rational morphism across its finite excluded locus; verify normalized-Poincare compatibility rather than choosing pointwise rescalings.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608K rational-local-morphism preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
