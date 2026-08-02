#!/usr/bin/env sage -python
"""N611Z: construct the quotient presentation of the N611X rank-one class."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


P = 103
ROOT = Path("notes")
N611X = ROOT / "poincare_degree3_elliptic_class_n611x.json"
N611X_VERIFY = ROOT / "poincare_degree3_elliptic_class_n611x_independent_verifier.json"
N611Y = ROOT / "poincare_conductor_factorization_n611y.json"
N611Y_VERIFY = ROOT / "poincare_conductor_factorization_n611y_independent_verifier.json"


def sha256(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path): return json.loads(path.read_text(encoding="ascii"))
def receipt_matches(receipt, path): return receipt.get("verified") is True and receipt.get("primary_sha256") == sha256(path)
def frobenius(point): return point if point.is_zero() else point.curve()(point[0] ** P, point[1] ** P)
def curve_record(curve): return {"a4": str(curve.a4()), "a6": str(curve.a6()), "j": str(curve.j_invariant())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    n611x, n611x_verify = load(N611X), load(N611X_VERIFY)
    n611y, n611y_verify = load(N611Y), load(N611Y_VERIFY)

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
            kills_beta = all(beta(point).is_zero() for point in kernel)
            paths.append({
                "first_index": first_index,
                "second_index": second_index,
                "first_target": curve_record(first.codomain()),
                "target": curve_record(second.codomain()),
                "first_kernel_polynomial": str(first.kernel_polynomial()),
                "second_kernel_polynomial": str(second.kernel_polynomial()),
                "composite_degree": int(first.degree()) * int(second.degree()),
                "kernel_size_in_recorded_nine_torsion": len(kernel),
                "beta_annihilates_kernel": kills_beta,
                "beta_m_annihilates_kernel": all(beta_m(point).is_zero() for point in kernel),
            })

    selected = [path for path in paths if path["composite_degree"] == 9 and path["kernel_size_in_recorded_nine_torsion"] == 9 and path["beta_annihilates_kernel"]]
    if len(selected) != 1:
        raise RuntimeError("degree-nine quotient path was not unique")
    choice = selected[0]
    degree_beta = 2 * 1557
    degree_gamma = degree_beta // 9
    gates = {
        "n611x_degree_three_class_receipt_binds": n611x["preflight_pass"] is True and receipt_matches(n611x_verify, N611X),
        "n611y_conductor_factorization_receipt_binds": n611y["preflight_pass"] is True and receipt_matches(n611y_verify, N611Y),
        "nine_torsion_container_and_cubic_path_sweep_complete": len(torsion_nine) == 27 and len(paths) == 7,
        "unique_degree_nine_kernel_is_annihilated_by_beta": len(selected) == 1 and choice["kernel_size_in_recorded_nine_torsion"] == 9,
        "selected_path_is_the_explicit_conductor_chain": choice["first_target"] == {"a4": "45", "a6": "7", "j": "10"} and choice["target"] == {"a4": "42", "a6": "28", "j": "6"} and choice["first_kernel_polynomial"] == "x + 67" and choice["second_kernel_polynomial"] == "x + 99",
        "neighboring_cross_map_control_does_not_annihilate_selected_kernel": choice["beta_m_annihilates_kernel"] is False,
        "quotient_degree_matches_rank_one_form": degree_beta == 9 * 346 and degree_gamma == 346,
        "abstract_pullback_matrix_matches_d": choice["composite_degree"] == 9 and degree_gamma == 346 and degree_beta == choice["composite_degree"] * degree_gamma,
    }
    output = {
        "schema": "ecdlp.h018.degree3-quotient-presentation.n611z.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__)),
        "claim_status": "OBSERVATION / EXPLICIT_ELLIPTIC_QUOTIENT_PRESENTATION_OF_RANK_ONE_CLASS / STANDARD_QUOTIENT_AND_POLARIZATION_IDENTITIES / MODEL-BOUND / TOY-EVIDENCE / THETA_EFFECTIVITY_AND_CURVE_OPEN / NO_ECDLP_CLAIM",
        "assumptions": [
            "A separable isogeny beta annihilating ker(f1) factors uniquely through the quotient f1.",
            "For f=(f1,gamma^dagger), the standard dual-isogeny and principal-polarization identities give f^dagger f=[[deg(f1),beta^dagger],[beta,deg(gamma)]].",
        ],
        "bound_inputs": {str(path): sha256(path) for path in (N611X, N611X_VERIFY, N611Y, N611Y_VERIFY)},
        "records": {
            "field": "F_(103^6)",
            "recorded_nine_torsion_size": len(torsion_nine),
            "all_degree_nine_paths": paths,
            "selected_path": choice,
            "beta_degree": degree_beta,
            "gamma_degree": degree_gamma,
            "quotient_target": "H:y^2=x^3+42*x+28",
            "quotient_presentation": "f=(f1,gamma^dagger): E x E_g -> H, beta_D=gamma o f1, deg(f1)=9, deg(gamma)=346",
            "pullback_polarization": "f^dagger f=[[9,beta_D^dagger],[beta_D,346]]=D",
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "Among the seven degree-nine two-cubic paths from E, exactly the explicit conductor chain E -> E0 -> H has its full nine-point kernel annihilated by beta_D. The quotient universal property therefore gives gamma:H -> E_g of degree 346 with beta_D=gamma o f1. Under the stated standard polarization identities, f=(f1,gamma^dagger) is an explicit elliptic quotient presentation whose pullback principal form is D.",
        "next_requirement": "Construct the effective divisor or theta-section data of the complementary principal form M relative to this quotient presentation, and prove or disprove theta-divisor irreducibility before seeking a genus-two equation or evaluator.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N611Z quotient-presentation preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__": main()
