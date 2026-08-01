#!/usr/bin/env sage -python
"""N608I: discharge the N606U stability and determinant prerequisites."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import sys
from pathlib import Path

from sage.all import EllipticCurve, GF, gcd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import product_kummer_h018_graph3pi_pole_budget_n606d as graph3pi

P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cover_fixture():
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    cover_map = (second * first).dual()
    return target, cover_map.domain(), cover_map


def kernel_generator(cover_map, cover):
    field = GF(P**6, name="a")
    cover_k = cover.base_extend(field)
    for root, _multiplicity in cover_map.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square():
            point = cover_k(x, rhs.sqrt())
            if point.order() == 9:
                return field, point
    raise RuntimeError("could not recover an order-nine cover kernel generator")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    base = GF(P)
    alpha_x_num, alpha_x_den, _alpha_y_num, _alpha_y_den = graph3pi.alpha_maps(base)
    alpha_degree = int(alpha_x_num.degree())
    target, cover, cover_map = cover_fixture()
    _field, deck = kernel_generator(cover_map, cover)

    source_rank = 9
    w9_degree = -1
    z_degree = alpha_degree
    source_after_pullback_degree = w9_degree * z_degree
    theta11_degree_addend = source_rank * 11
    source_degree = source_after_pullback_degree + theta11_degree_addend
    source_determinant_origin_coefficient = source_after_pullback_degree + theta11_degree_addend
    cover_rank = int(cover_map.degree())
    cover_degree = 2
    kernel_indices = list(range(9))
    translate_stabilizer_indices = [
        index for index in kernel_indices if (2 * index) * deck == deck.curve()(0)
    ]

    records = {
        "curve": str(target),
        "cover": str(cover),
        "characteristic": P,
        "alpha_name": "3+pi",
        "z_name": "-3-pi=-alpha",
        "alpha_x_degree": alpha_degree,
        "z_is_separable": alpha_degree % P != 0,
        "fourier_mukai_w9": {"rank": source_rank, "degree": w9_degree, "determinant": "O(-O)"},
        "h018_pushforward": {
            "rank": source_rank,
            "degree_after_z_pullback": source_after_pullback_degree,
            "theta11_degree_addend": theta11_degree_addend,
            "degree": source_degree,
            "determinant": f"O({source_determinant_origin_coefficient}O)",
        },
        "cyclic_cover": {
            "degree": int(cover_map.degree()),
            "kernel_generator_order": int(deck.order()),
            "kernel_translation_stabilizer_indices_for_O2": translate_stabilizer_indices,
            "nonzero_kernel_translate_count": 8,
            "all_nonzero_translates_of_O2_distinct": translate_stabilizer_indices == [0],
            "determinant_pi_pushforward_O": "O",
            "norm_O2": "O(2O)",
            "pushforward_rank": cover_rank,
            "pushforward_degree": cover_degree,
            "pushforward_determinant": "O(2O)",
        },
        "coprime_rank_degree": {
            "source_gcd": int(gcd(source_rank, source_degree)),
            "cover_gcd": int(gcd(cover_rank, cover_degree)),
        },
        "standard_theorem_dependencies": [
            "Fourier--Mukai transform of O(9O) is stable of rank 9, degree -1, determinant O(-O).",
            "Separable pullback preserves semistability on a smooth projective curve.",
            "A semistable elliptic bundle with coprime rank and degree is stable.",
            "For a finite etale Galois cover, pullback of pi_*L is the direct sum of kernel translates of L.",
            "Stable elliptic bundles with coprime rank and degree are classified by determinant.",
        ],
    }
    gates = {
        "alpha_degree_97": alpha_degree == 97,
        "z_separable": records["z_is_separable"],
        "source_rank_degree_is_9_2": source_rank == 9 and source_degree == 2,
        "source_determinant_is_O2": source_determinant_origin_coefficient == 2,
        "cover_is_cyclic_degree_9": int(cover_map.degree()) == 9 and int(deck.order()) == 9,
        "cover_odd_kernel_translate_control": translate_stabilizer_indices == [0],
        "cover_rank_degree_is_9_2": cover_rank == 9 and cover_degree == 2,
        "cover_determinant_is_O2": records["cyclic_cover"]["pushforward_determinant"] == "O(2O)",
        "both_rank_degree_pairs_coprime": records["coprime_rank_degree"]["source_gcd"] == 1 and records["coprime_rank_degree"]["cover_gcd"] == 1,
    }
    output = {
        "schema": "ecdlp.h018.poincare-pushforward-stability.n608i.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "RESTRICTED THEOREM / ABSTRACT BUNDLE ISOMORPHISM / STANDARD-THEOREM-BOUND / MODEL-BOUND / TOY-EVIDENCE / NO EXPLICIT INTERTWINER / NO_ECDLP_CLAIM",
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "abstract_bundle_isomorphism_discharged": all(gates.values()),
        "explicit_h018_to_cover_intertwiner_constructed": False,
        "strongest_valid_statement": "Using the named standard elliptic-bundle theorems, the N606U H018 pushforward and pi_*O_Eprime(2Oprime) are stable rank-nine degree-two bundles with determinant O(2O), hence are isomorphic on this F_103 fixture. This is an existence theorem only: it supplies no matrix, section evaluator, source inverse, relation, descent, or ECDLP speedup.",
        "next_requirement": "Construct the isomorphism as an explicit local matrix with declared charts and verify it against the normalized Poincare fibre model; do not substitute a finite-orbit gauge.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608I stability/determinant audit failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
