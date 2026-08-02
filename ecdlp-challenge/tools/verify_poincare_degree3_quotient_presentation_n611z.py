#!/usr/bin/env sage -python
"""Independent replay for N611Z's cubic quotient presentation."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


P = 103


def sha256(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def frobenius(point): return point if point.is_zero() else point.curve()(point[0] ** P, point[1] ** P)
def curve_record(curve): return {"a4": str(curve.a4()), "a6": str(curve.a6()), "j": str(curve.j_invariant())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))

    field = GF(P**6, name="a")
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    torsion_nine = curve(0).division_points(9)
    ring = PolynomialRing(field, "T")
    root = sorted(value for value, _ in (ring.gen()**3 + ring.gen() + 24).roots())[0]
    phi = curve.isogeny(curve(root, 0))
    def beta(point): return phi(13 * point + 4 * frobenius(point))
    def beta_m(point): return phi(12 * point + 4 * frobenius(point))

    paths = []
    for first_index, first in enumerate(curve.isogenies_prime_degree(3)):
        for second_index, second in enumerate(first.codomain().isogenies_prime_degree(3)):
            kernel = [point for point in torsion_nine if second(first(point)).is_zero()]
            paths.append({
                "first_index": first_index,
                "second_index": second_index,
                "first_target": curve_record(first.codomain()),
                "target": curve_record(second.codomain()),
                "first_kernel_polynomial": str(first.kernel_polynomial()),
                "second_kernel_polynomial": str(second.kernel_polynomial()),
                "composite_degree": int(first.degree()) * int(second.degree()),
                "kernel_size_in_recorded_nine_torsion": len(kernel),
                "beta_annihilates_kernel": all(beta(point).is_zero() for point in kernel),
                "beta_m_annihilates_kernel": all(beta_m(point).is_zero() for point in kernel),
            })
    selected = [path for path in paths if path["composite_degree"] == 9 and path["kernel_size_in_recorded_nine_torsion"] == 9 and path["beta_annihilates_kernel"]]
    records = primary["records"]
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.degree3-quotient-presentation.n611z.v1",
        "primary_gates": primary["preflight_pass"] is True and all(primary["gates"].values()),
        "replayed_path_table": paths == records["all_degree_nine_paths"],
        "unique_selected_path": len(selected) == 1 and selected[0] == records["selected_path"],
        "selected_explicit_models": selected[0]["first_target"] == {"a4": "45", "a6": "7", "j": "10"} and selected[0]["target"] == {"a4": "42", "a6": "28", "j": "6"},
        "kernel_and_control": selected[0]["kernel_size_in_recorded_nine_torsion"] == 9 and selected[0]["beta_annihilates_kernel"] and not selected[0]["beta_m_annihilates_kernel"],
        "degree_presentation": records["beta_degree"] == 3114 and records["gamma_degree"] == 346 and records["beta_degree"] == 9 * records["gamma_degree"],
    }
    output = {
        "verified": all(checks.values()),
        "primary_sha256": sha256(args.primary),
        "checks": checks,
        "strongest_valid_statement": "The independent replay reconstructs all seven cubic paths, the unique kernel inclusion for beta_D, the beta_M negative control, and the degree-9-by-346 quotient presentation. Effectivity, theta irreducibility, a genus-two equation, and all ECDLP algorithmic gates remain open.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]: raise RuntimeError("N611Z verifier failed")
    print(json.dumps({"output": str(args.out), "verified": True}, sort_keys=True))


if __name__ == "__main__": main()
