#!/usr/bin/env sage -python
"""N609I: bind the compact P fibres to the conditional global H018 class."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path


P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cm_mul(left, right):
    a, b = left
    c, d = right
    return (a * c - P * b * d, a * d + b * c - 5 * b * d)


def cm_conjugate(value):
    a, b = value
    return (a - 5 * b, -b)


def cm_norm(value):
    product = cm_mul(value, cm_conjugate(value))
    if product[1] != 0:
        raise RuntimeError("CM norm unexpectedly non-rational")
    return product[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    n608r = json.loads((root / "notes/poincare_global_sections_n608r.json").read_text(encoding="ascii"))
    n609h = json.loads((root / "notes/poincare_residual_pole_compactification_n609h.json").read_text(encoding="ascii"))
    n609h_replay = json.loads((root / "notes/poincare_residual_pole_compactification_n609h_structural_replay.json").read_text(encoding="ascii"))
    z = (-3, -1)
    zbar = cm_conjugate(z)
    determinant = 9 * 11 - cm_norm(z)
    h018 = [[9, zbar], [z, 11]]
    swapped = [[9, z], [zbar, 11]]
    p_fibre = "O(9O) tensor O([-4Q]-O) = O(8O+[-4Q])"
    records = {
        "cm_relation": "F^2+5F+103=0",
        "z": list(z),
        "zbar": list(zbar),
        "h018_matrix": h018,
        "determinant": determinant,
        "self_intersection": 2 * determinant,
        "conditional_arithmetic_genus": 1 + determinant,
        "base_field_z_action": "[-4]",
        "p_fibre_restriction": p_fibre,
        "p_fibre_degree": 9,
        "q_fibre_degree": 11,
        "n608r_selected_cover_subspace": n608r["records"]["selected_cover_subspace"],
        "n608r_preflight_pass": n608r["preflight_pass"],
        "n609h_preflight_pass": n609h["preflight_pass"],
        "n609h_replay_pass": n609h_replay["pass"],
        "n609h_nonzero_q_fibre_count": n609h["records"]["nonzero_q_fibre_count"],
        "n609h_pole_divisor": n609h["records"]["nonzero_q_pole_divisor"],
    }
    gates = {
        "h018_orientation_exact": h018 == [[9, (2, 1)], [(-3, -1), 11]],
        "swapped_orientation_rejected": swapped != h018,
        "determinant_and_intersection_exact": determinant == 2 and records["self_intersection"] == 4,
        "base_field_poincare_action_matches_compact_p_fibres": records["base_field_z_action"] == "[-4]"
        and records["p_fibre_restriction"] == "O(9O) tensor O([-4Q]-O) = O(8O+[-4Q])"
        and records["p_fibre_degree"] == 9,
        "q_projection_degree_obligation_exact": records["q_fibre_degree"] == 11,
        "bound_section_and_compactification_receipts_pass": records["n608r_preflight_pass"]
        and records["n608r_selected_cover_subspace"] == "span{1,x}"
        and records["n609h_preflight_pass"]
        and records["n609h_replay_pass"]
        and records["n609h_nonzero_q_fibre_count"] == 108
        and records["n609h_pole_divisor"] == "8[O]+[-4Q]",
    }
    output = {
        "schema": "ecdlp.h018.poincare-global-class.n609i.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "RESTRICTED THEOREM / CONDITIONAL GLOBAL H018 CLASS OF COMPACT RESIDUAL / STANDARD-FACT-BOUND / MODEL-BOUND / TOY-EVIDENCE / EXPLICIT GLOBAL EQUATION OPEN / NO_ECDLP_CLAIM",
        "assumptions": [
            "N608O's normalized H018 global pushforward model.",
            "N608R's selected global cover section model.",
            "The local N609H evaluator is the restriction of that selected section.",
        ],
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "Under the declared normalized H018 global-section model, the N609H compact P fibres are the restrictions of an H018 divisor class. The class has P-fibre degree 9, required Q-fibre degree 11, self-intersection 4, and conditional arithmetic genus 3 for an effective Cartier representative.",
        "next_requirement": "Construct an explicit global section or product-surface equation realizing this class, then prove Cartier/smoothness and analyze normalization before using the degree/genus target for a Jacobian, factor-base, relation, rank, descent, cost, or ECDLP claim.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N609I global-class preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
