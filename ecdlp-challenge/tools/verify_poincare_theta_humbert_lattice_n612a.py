#!/usr/bin/env sage -python
"""Independent replay of N612A's Hom-lattice and Humbert minimum."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QuadraticForm, ZZ, matrix


P = 103


def sha256(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def frobenius(point): return point if point.is_zero() else point.curve()(point[0] ** P, point[1] ** P)
def norm(r, s): return r * r - 5 * r * s + 103 * s * s
def q_ns(vector):
    a, b, r, s = vector.list() if hasattr(vector, "list") else vector
    return a * b - 2 * norm(r, s)
def humbert(vector):
    a, b, r, s = vector.list() if hasattr(vector, "list") else vector
    intersection = 345 * a + 9 * b - 8 * r - 1528 * s
    return intersection * intersection - 4 * q_ns(vector)
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
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))

    field = GF(P**6, name="a")
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    ring = PolynomialRing(field, "T")
    root = sorted(value for value, _ in (ring.gen()**3 + ring.gen() + 24).roots())[0]
    phi = curve.isogeny(curve(root, 0))
    two_torsion = curve(0).division_points(2)
    multipliers = []
    for coeff_a, coeff_b, label in ((0, 0, "0"), (1, 0, "1"), (0, 1, "Frob"), (1, 1, "1+Frob")):
        multipliers.append({"label": label, "coefficients_mod_2": [coeff_a, coeff_b], "annihilates_all_e_two": all(phi(coeff_a * point + coeff_b * frobenius(point)).is_zero() for point in two_torsion)})

    theta = matrix(ZZ, 4, 1, [9, 345, 12, 4])
    divisor = matrix(ZZ, 4, 1, [9, 346, 13, 4])
    _smith_diagonal, left, _right = theta.smith_form()
    transport = matrix(ZZ, left.inverse())
    quotient_form = lambda vector: humbert(transport * matrix(ZZ, 4, 1, [0, vector[0], vector[1], vector[2]]))
    coefficients = [int(value) for value in ternary_coefficients(quotient_form)]
    ternary = QuadraticForm(ZZ, 3, coefficients)
    short_to_eight = ternary.short_vector_list_up_to_length(8)
    witness = matrix(ZZ, 3, 1, [4, 3, -3])
    records = primary["records"]
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.theta-humbert-lattice.n612a.v1",
        "primary_gates": primary["preflight_pass"] is True and all(primary["gates"].values()),
        "two_torsion_replay": len(two_torsion) == 4 and sum(phi(point).is_zero() for point in two_torsion) == 2 and multipliers == records["mod_two_cm_multiplier_controls"],
        "unimodular_lattice_replay": left * theta == matrix(ZZ, 4, 1, [1, 0, 0, 0]) and matrix(ZZ, 4, 1, transport.column(0)) == theta and abs(transport.det()) == 1 and records["quotient_basis_matrix"] == [[int(value) for value in row] for row in transport.rows()],
        "humbert_form_replay": coefficients == [259512, 8136, 699824, 72, 10992, 471817] and records["ternary_humbert_coefficients"] == coefficients and ternary.is_positive_definite() and ternary.Gram_matrix().det() == 6192,
        "exact_minimum_replay": all(not vectors for vectors in short_to_eight[1:]) and ternary(witness.list()) == 9,
        "n611x_witness_replay": q_ns(divisor) == 0 and humbert(divisor) == 9 and transport * matrix(ZZ, 4, 1, [0, 4, 3, -3]) == -divisor,
    }
    output = {
        "verified": all(checks.values()),
        "primary_sha256": sha256(args.primary),
        "checks": checks,
        "strongest_valid_statement": "The independent replay reconstructs the full-Hom mod-two gate and the integral ternary Humbert form with no value through eight and a value-nine N611X witness. The resulting irreducibility conclusion remains conditional on the stated CM, Neron-Severi, and Humbert criteria and does not supply an explicit curve equation or ECDLP algorithm.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]: raise RuntimeError("N612A verifier failed")
    print(json.dumps({"output": str(args.out), "verified": True}, sort_keys=True))


if __name__ == "__main__": main()
