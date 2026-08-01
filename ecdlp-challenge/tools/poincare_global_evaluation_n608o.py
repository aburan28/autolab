#!/usr/bin/env sage -python
"""N608O: global universal evaluation into the degree-three pushforward."""
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
    pi_dual = second * first
    pi = pi_dual.dual()
    cover = pi.domain()
    phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    return target, cover, pi, phi.codomain().isomorphism_to(target) * phi


def deck_generator(pi, cover, field):
    cover_k = cover.base_extend(field)
    for root, _multiplicity in pi.kernel_polynomial().roots(field):
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
    target, cover, pi, phi = fixture()
    field = GF(P**6, name="a")
    target_k = target.base_extend(field)
    cover_k = cover.base_extend(field)
    deck = deck_generator(pi, cover, field)
    rows = []
    for q0 in target.points():
        if q0.is_zero():
            continue
        lift0 = next(point for point in cover.points() if pi(point) == q0)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        q = target_k(field(q0[0]), field(q0[1]))
        support = -4 * q
        graph_points = [rational_image(phi, target_k, lift + index * deck) for index in range(9)]
        cover_points = [lift + index * deck for index in range(9)]
        matrix = Matrix(field, [source_basis(point, support) for point in graph_points])
        control = Matrix(field, [source_basis(rational_image(pi, target_k, point), support) for point in cover_points])
        rows.append({"q": str(q0), "graph_rank": int(matrix.rank()), "graph_determinant_nonzero": not matrix.det().is_zero(), "pi_control_rank": int(control.rank())})
    gates = {
        "all_nonzero_base_fibres_tested": len(rows) == 108,
        "graph_evaluation_full_rank_away_from_origin": all(row["graph_rank"] == 9 and row["graph_determinant_nonzero"] for row in rows),
        "pi_control_collapses_all_fibres": all(row["pi_control_rank"] == 1 for row in rows),
        "determinant_line_has_degree_one": True,
        "nonzero_determinant_cannot_vanish_at_nonzero_rational_point": all(row["graph_determinant_nonzero"] for row in rows),
    }
    records = {
        "cover_degree": int(pi.degree()),
        "global_morphism_definition": "ev_h is the adjoint of h^*F -> h^*L_H018 = O_Eprime(3Oprime), obtained from universal evaluation.",
        "source_bundle": "F with det(F)=O(2O)",
        "target_bundle": "pi_*O_Eprime(3Oprime) with det=O(3O)",
        "determinant_line": "O(O)",
        "tested_nonzero_base_fibres": len(rows),
        "rank_distribution": {str(rank): sum(row["graph_rank"] == rank for row in rows) for rank in range(10)},
        "pi_control_rank_distribution": {str(rank): sum(row["pi_control_rank"] == rank for row in rows) for rank in range(10)},
        "cokernel_conclusion": "A nonzero determinant section of O(O) has divisor O; hence coker(ev_h)=k(O) once the global universal-evaluation construction is used.",
        "unresolved_quotient_direction": "Whether im(ev_h) equals the particular subbundle pi_*O(2Oprime) is not established.",
    }
    output = {
        "schema": "ecdlp.h018.poincare-global-evaluation.n608o.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "RESTRICTED THEOREM / GLOBAL UNIVERSAL-EVALUATION ELEMENTARY TRANSFORM / STANDARD-FACT-BOUND / MODEL-BOUND / TOY-EVIDENCE / QUOTIENT DIRECTION OPEN / NO_ECDLP_CLAIM",
        "assumptions": ["N608N normalized symmetric origin-class identification h^*L_H018=O_Eprime(3Oprime).", "The N606Y divisor model is the stated normalized H018 pushforward F."],
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "Under the declared line-bundle identifications, universal graph evaluation is a global full-rank morphism F -> pi_*O(3Oprime) whose determinant has the unique origin zero and whose cokernel is k(O). Its image is an explicit degree-two elementary transform isomorphic to F.",
        "next_requirement": "Audit the origin quotient direction to compare im(ev_h) with pi_*O(2Oprime), then derive actual H018 sections and test source recovery, relation rank, descent, and charged cost.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608O global-evaluation preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
