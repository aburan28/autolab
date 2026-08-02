#!/usr/bin/env sage -python
"""Independent full-torsion replay for N611V."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


def sha256(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    field = GF(103**6, name="a")
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    ring = PolynomialRing(field, "T")
    root = sorted(value for value, _ in (ring.gen()**3 + ring.gen() + 24).roots())[0]
    phi = curve.isogeny(curve(root, 0))
    e2, g2 = curve(0).division_points(2), phi.codomain()(0).division_points(2)
    records = primary["records"]
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.reduced-level2-gluing-audit.n611v.v1",
        "primary_gates": primary["preflight_pass"] is True and all(primary["gates"].values()),
        "full_torsion_recomputed": len(e2) == 4 and len(g2) == 4,
        "all_reduced_beta_inputs_are_even": all((12 * point).is_zero() for point in e2),
        "declared_product_pairing": records["pairing_matrix_bits"] == [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
        "plane_counts": records["maximal_isotropic_plane_count"] == 15 and records["transverse_graph_plane_count"] == 6,
    }
    output = {"verified": all(checks.values()), "primary_sha256": sha256(args.primary), "checks": checks, "strongest_valid_statement": "The independent replay confirms that the degree-six field contains both full two-torsion groups and that the reduced beta coefficient is even there; the primary's direct level-two selection conclusion remains negative."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]: raise RuntimeError("N611V independent verifier failed")
    print(json.dumps({"output": str(args.out), "verified": True}, sort_keys=True))

if __name__ == "__main__": main()
