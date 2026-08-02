#!/usr/bin/env sage -python
"""N611W: determine the 2-adic invariant factors of the reduced cross map."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Integers, Matrix, PolynomialRing, identity_matrix


P = 103
ROOT = Path("notes")
N611V = ROOT / "poincare_reduced_level2_gluing_n611v.json"
N611V_VERIFY = ROOT / "poincare_reduced_level2_gluing_n611v_independent_verifier.json"


def sha256(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path): return json.loads(path.read_text(encoding="ascii"))
def receipt_matches(receipt, path): return receipt.get("verified") is True and receipt.get("primary_sha256") == sha256(path)


def frobenius_order(modulus):
    ring = Integers(modulus)
    frobenius = Matrix(ring, [[0, -103], [1, -5]])
    power, identity = identity_matrix(ring, 2), identity_matrix(ring, 2)
    for degree in range(1, 100):
        power *= frobenius
        if power == identity:
            return degree
    raise RuntimeError("Frobenius order not found")


def torsion_image_size(level, extension_degree):
    field = GF(P ** extension_degree, name="a")
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    ring = PolynomialRing(field, "T")
    root = sorted(value for value, _ in (ring.gen() ** 3 + ring.gen() + 24).roots())[0]
    phi = curve.isogeny(curve(root, 0))
    def frobenius(point): return point if point.is_zero() else point.curve()(point[0] ** P, point[1] ** P)
    def beta(point): return phi(12 * point + 4 * frobenius(point))
    torsion = curve(0).division_points(level)
    return len(torsion), len({str(beta(point)) for point in torsion})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    n611v, n611v_verify = load(N611V), load(N611V_VERIFY)
    rows = []
    expected = {4: 1, 8: 2, 16: 8}
    extension_degrees = {4: 6, 8: 12, 16: 24}
    for level in (4, 8, 16):
        torsion_size, image_size = torsion_image_size(level, extension_degrees[level])
        rows.append({"level": level, "extension_degree": extension_degrees[level], "frobenius_order": frobenius_order(level), "torsion_size": torsion_size, "image_size": image_size, "expected_image_size": expected[level]})
    gates = {
        "n611v_level_two_receipt_binds": n611v["preflight_pass"] is True and receipt_matches(n611v_verify, N611V),
        "unit_factor_has_odd_norm": 3 * 3 - 5 * 3 + 103 == 97,
        "frobenius_extension_degrees_are_exact": [row["frobenius_order"] for row in rows] == [6, 12, 24],
        "full_torsion_sizes_replay": [row["torsion_size"] for row in rows] == [16, 64, 256],
        "smith_factor_image_sizes_replay": all(row["image_size"] == row["expected_image_size"] for row in rows),
        "cross_map_is_not_direct_graph_isomorphism": all(row["image_size"] < row["torsion_size"] for row in rows),
    }
    output = {
        "schema": "ecdlp.h018.reduced-2adic-crossmap.n611w.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__)),
        "claim_status": "RESTRICTED THEOREM / NEGATIVE_DIRECT_POWER_OF_TWO_FREY_KANI_GRAPH_BY_REDUCED_CROSS_MAP / MODEL-BOUND / TOY-EVIDENCE / QUOTIENT_LATTICE_AND_HIGHER_THETA_OPEN / NO_ECDLP_CLAIM",
        "assumptions": ["The degree-two separable isogeny phi_g has two-adic Smith factors (1,2).", "The odd-norm CM map 3+pi is a two-adic unit."],
        "bound_inputs": {str(path): sha256(path) for path in (N611V, N611V_VERIFY)},
        "records": {"cross_map_factorization": "beta=phi_g o [4] o ([3]+Frob)", "unit_norm": 97, "predicted_smith_factors": [4, 8], "rows": rows},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The reduced cross map has the predicted 2-adic factors (4,8): exact full-torsion image sizes are 1, 2, and 8 at levels 4, 8, and 16. It supplies nonzero rank-one residue at level 8 and higher residue at level 16, but is never a direct power-of-two torsion isomorphism and cannot itself define a Frey-Kani graph.",
        "next_requirement": "Use the nontrivial 2-adic residue only through an explicit quotient-lattice or theta-group construction; direct cross-map graph gluing is closed before relation work.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N611W 2-adic cross-map preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__": main()
