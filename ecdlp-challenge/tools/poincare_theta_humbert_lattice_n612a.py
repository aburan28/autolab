#!/usr/bin/env sage -python
"""N612A: certify the full Hom lattice and Humbert irreducibility condition."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QuadraticForm, ZZ, gcd, matrix


P = 103
ROOT = Path("notes")
N611X = ROOT / "poincare_degree3_elliptic_class_n611x.json"
N611X_VERIFY = ROOT / "poincare_degree3_elliptic_class_n611x_independent_verifier.json"
N611Z = ROOT / "poincare_degree3_quotient_presentation_n611z.json"
N611Z_VERIFY = ROOT / "poincare_degree3_quotient_presentation_n611z_independent_verifier.json"


def sha256(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path): return json.loads(path.read_text(encoding="ascii"))
def receipt_matches(receipt, path): return receipt.get("verified") is True and receipt.get("primary_sha256") == sha256(path)
def frobenius(point): return point if point.is_zero() else point.curve()(point[0] ** P, point[1] ** P)


def norm(r, s): return r * r - 5 * r * s + 103 * s * s
def q_ns(vector):
    a, b, r, s = vector.list() if hasattr(vector, "list") else vector
    return a * b - 2 * norm(r, s)
def intersection_with_theta(vector):
    a, b, r, s = vector.list() if hasattr(vector, "list") else vector
    return 345 * a + 9 * b - 8 * r - 1528 * s
def humbert(vector): return intersection_with_theta(vector) ** 2 - 4 * q_ns(vector)


def ternary_coefficients(form):
    basis = [tuple(int(index == coordinate) for index in range(3)) for coordinate in range(3)]
    coefficients = []
    for left in range(3):
        coefficients.append(form(basis[left]))
        for right in range(left + 1, 3):
            summed = tuple(basis[left][index] + basis[right][index] for index in range(3))
            coefficients.append(form(summed) - form(basis[left]) - form(basis[right]))
    return coefficients


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    n611x, n611x_verify = load(N611X), load(N611X_VERIFY)
    n611z, n611z_verify = load(N611Z), load(N611Z_VERIFY)

    field = GF(P**6, name="a")
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    ring = PolynomialRing(field, "T")
    root = sorted(value for value, _ in (ring.gen()**3 + ring.gen() + 24).roots())[0]
    phi = curve.isogeny(curve(root, 0))
    two_torsion = curve(0).division_points(2)
    multipliers = []
    for coeff_a, coeff_b, label in ((0, 0, "0"), (1, 0, "1"), (0, 1, "Frob"), (1, 1, "1+Frob")):
        annihilates = all(phi(coeff_a * point + coeff_b * frobenius(point)).is_zero() for point in two_torsion)
        multipliers.append({"label": label, "coefficients_mod_2": [coeff_a, coeff_b], "annihilates_all_e_two": annihilates})

    theta = matrix(ZZ, 4, 1, [9, 345, 12, 4])
    divisor = matrix(ZZ, 4, 1, [9, 346, 13, 4])
    _smith_diagonal, left, _right = theta.smith_form()
    transport = matrix(ZZ, left.inverse())
    quotient_form = lambda vector: humbert(transport * matrix(ZZ, 4, 1, [0, vector[0], vector[1], vector[2]]))
    coefficients = [int(value) for value in ternary_coefficients(quotient_form)]
    ternary = QuadraticForm(ZZ, 3, coefficients)
    short_to_eight = ternary.short_vector_list_up_to_length(8)
    witness_coordinates = matrix(ZZ, 3, 1, [4, 3, -3])
    witness_vector = transport * matrix(ZZ, 4, 1, [0, *witness_coordinates.list()])
    min_vectors = ternary.short_vector_list_up_to_length(10)[9]

    gates = {
        "n611x_rank_one_class_receipt_binds": n611x["preflight_pass"] is True and receipt_matches(n611x_verify, N611X),
        "n611z_quotient_presentation_receipt_binds": n611z["preflight_pass"] is True and receipt_matches(n611z_verify, N611Z),
        "full_two_torsion_and_phi_kernel_materialized": len(two_torsion) == 4 and sum(phi(point).is_zero() for point in two_torsion) == 2,
        "only_even_cm_multiplier_annihilates_all_e_two": multipliers[0]["annihilates_all_e_two"] is True and all(not item["annihilates_all_e_two"] for item in multipliers[1:]),
        "theta_is_primitive_principal_ns_class": q_ns(theta) == 1 and gcd(theta.list()) == 1,
        "unimodular_humbert_quotient_basis": left * theta == matrix(ZZ, 4, 1, [1, 0, 0, 0]) and matrix(ZZ, 4, 1, transport.column(0)) == theta and abs(transport.det()) == 1,
        "quotient_humbert_form_is_positive_integral": coefficients == [259512, 8136, 699824, 72, 10992, 471817] and ternary.is_positive_definite() and ternary.Gram_matrix().det() == 6192,
        "exact_humbert_minimum_is_nine": all(not vectors for vectors in short_to_eight[1:]) and ternary(witness_coordinates.list()) == 9 and len(min_vectors) == 2,
        "n611x_divisor_is_value_nine_witness": q_ns(divisor) == 0 and intersection_with_theta(divisor) == 3 and humbert(divisor) == 9 and witness_vector == -divisor,
    }
    output = {
        "schema": "ecdlp.h018.theta-humbert-lattice.n612a.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__)),
        "claim_status": "RESTRICTED THEOREM / FULL_CM_HOM_LATTICE_AND_GEOMETRIC_THETA_IRREDUCIBILITY / STANDARD_CM_ENDOMORPHISM_AND_HUMBERT_CRITERION_BOUND / MODEL-BOUND / TOY-EVIDENCE / CURVE_EQUATION_AND_EVALUATOR_OPEN / NO_ECDLP_CLAIM",
        "assumptions": [
            "End(E)=Z[pi] with pi^2+5*pi+103=0, and phi_g^dagger phi_g=[2].",
            "The standard Hermitian Neron-Severi identification and Kani/Weil refined Humbert irreducibility criterion apply over the algebraic closure.",
        ],
        "bound_inputs": {str(path): sha256(path) for path in (N611X, N611X_VERIFY, N611Z, N611Z_VERIFY)},
        "records": {
            "hom_lattice_argument": "For h:E->E_g, put v=phi_g^dagger h. Then phi_g v=[2]h annihilates E[2]. The recorded mod-two action forces v=2u, so h=phi_g u.",
            "two_torsion_size": len(two_torsion),
            "mod_two_cm_multiplier_controls": multipliers,
            "ns_quadratic_form": "q_A(a,b,r,s)=a*b-2*(r^2-5*r*s+103*s^2)",
            "theta_vector": [int(value) for value in theta.list()],
            "humbert_form": "q_M(D)=(D.M)^2-4*q_A(D)",
            "smith_left_matrix": [[int(value) for value in row] for row in left.rows()],
            "quotient_basis_matrix": [[int(value) for value in row] for row in transport.rows()],
            "ternary_humbert_coefficients": coefficients,
            "ternary_humbert_determinant": int(ternary.Gram_matrix().det()),
            "minimum_value": 9,
            "minimum_witness_quotient_coordinates": [int(value) for value in witness_coordinates.list()],
            "minimum_witness_ns_vector": [int(value) for value in witness_vector.list()],
            "n611x_divisor_vector": [int(value) for value in divisor.list()],
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "Under the stated CM endomorphism, Hermitian Neron-Severi, and refined Humbert assumptions, the full Hom lattice is phi_g Z[pi] and the reduced Humbert form of M has exact minimum nine. It therefore never represents one, so the associated principal theta class is geometrically irreducible. This establishes existence of a smooth genus-two theta curve over an algebraic closure, but not an equation, descent datum, evaluator, or ECDLP algorithm.",
        "next_requirement": "Construct an explicit equation and degree-three maps for the irreducible genus-two theta curve, then build and charge its evaluator before any relation-collection or rho-comparison claim.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N612A Humbert-lattice preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__": main()
