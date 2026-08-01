#!/usr/bin/env sage -python
"""N608L: exact polarization degree of the N608J degree-eleven correspondence."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, gcd

P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add_points(left, right, rhs):
    x_left, r_left = left
    x_right, r_right = right
    slope = (r_right - r_left) / (x_right - x_left)
    x_sum = slope**2 * rhs - x_left - x_right
    r_sum = -r_left + slope * (x_left - x_sum)
    return x_sum, r_sum


def double_point(point, curve_a, rhs):
    x_value, r_value = point
    slope = (3 * x_value**2 + curve_a) / (2 * r_value * rhs)
    x_double = slope**2 * rhs - 2 * x_value
    r_double = -r_value + slope * (x_value - x_double)
    return x_double, r_double


def m_plus_frobenius_x_map(curve, multiplier):
    base = curve.base_field()
    ring = PolynomialRing(base, "X")
    X = ring.gen()
    field = ring.fraction_field()
    x = field(X)
    rhs = x**3 + field(curve.a4()) * x + field(curve.a6())
    point = (x, field.one())
    multiple = point
    bits = bin(multiplier)[3:]
    for bit in bits:
        multiple = double_point(multiple, field(curve.a4()), rhs)
        if bit == "1":
            multiple = add_points(multiple, point, rhs)
    frobenius = (x**P, rhs ** ((P - 1) // 2))
    image = add_points(multiple, frobenius, rhs)[0]
    numerator = ring(image.numerator())
    denominator = ring(image.denominator())
    common = gcd(numerator, denominator)
    return numerator // common, denominator // common


def multiply_cm(left, right):
    """Multiply a+bF using F^2+5F+103=0."""
    a, b = left
    c, d = right
    return (a * c - 103 * b * d, a * d + b * c - 5 * b * d)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi_dual = second * first
    pi = pi_dual.dual()
    cover = pi.domain()
    phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    phi = phi.codomain().isomorphism_to(target) * phi
    alpha = phi.dual() * pi
    candidate_num, candidate_den = m_plus_frobenius_x_map(cover, 4)
    alpha_x, _alpha_y = alpha.rational_maps()
    ring = candidate_num.parent()
    alpha_num = ring({exponent[0]: coefficient for exponent, coefficient in alpha_x.numerator().dict().items()})
    alpha_den = ring({exponent[0]: coefficient for exponent, coefficient in alpha_x.denominator().dict().items()})
    x_maps_equal = alpha_num * candidate_den == candidate_num * alpha_den

    generator = next(point for point in cover.points() if not point.is_zero())
    frobenius_generator = cover(generator[0] ** P, generator[1] ** P)
    sign_check = alpha(generator) == 4 * generator + frobenius_generator
    trace = int(cover.trace_of_frobenius())
    a = (4, 1)
    b = (-1, -1)  # Rosati dual of 4+F, since Fbar=-5-F.
    cross_one = multiply_cm((2, 1), a)
    cross_two = multiply_cm((-3, -1), b)
    induced = (198 + cross_one[0] + cross_two[0], cross_one[1] + cross_two[1])
    direct_scalar = next(index for index in range(int(cover.cardinality())) if index * generator == (phi.dual()(9 * phi(generator) + 3 * pi(generator)) + pi_dual(-4 * phi(generator) + 11 * pi(generator))))

    records = {
        "cover": str(cover),
        "trace_of_frobenius": trace,
        "frobenius_relation": "F^2+5F+103=0",
        "pi_degree": int(pi.degree()),
        "phi_degree": int(phi.degree()),
        "alpha_phi_dual_pi_degree": int(alpha.degree()),
        "alpha_candidate": "[4]+F",
        "alpha_x_maps_equal_exactly": x_maps_equal,
        "alpha_sign_check_on_generator": sign_check,
        "alpha_cm_coordinates": list(a),
        "dual_alpha_cm_coordinates": list(b),
        "h018_cross_terms": {"phi_dual_zbar_pi": list(cross_one), "pi_dual_z_phi": list(cross_two)},
        "h018_induced_polarization_cm_coordinates": list(induced),
        "direct_base_group_scalar_control": direct_scalar,
        "pullback_line_degree": induced[0] if induced[1] == 0 else None,
        "target_line_degree": 2,
    }
    gates = {
        "expected_frobenius_trace": trace == -5,
        "alpha_rational_map_exact": x_maps_equal and sign_check,
        "cross_terms_cancel_frobenius": cross_one[1] + cross_two[1] == 0,
        "induced_polarization_is_three": induced == (3, 0) and direct_scalar == 3,
        "degree_two_direct_extension_rejected": records["pullback_line_degree"] == 3 != records["target_line_degree"],
    }
    output = {
        "schema": "ecdlp.h018.poincare-evaluation-degree.n608l.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "NEGATIVE RESULT / DIRECT GLOBAL EXTENSION OF N608K REJECTED / MODEL-BOUND / TOY-EVIDENCE / OTHER CORRESPONDENCES OR MODIFICATIONS OPEN / NO_ECDLP_CLAIM",
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The N608J correspondence pulls the H018 line bundle back with induced polarization [3], hence degree three, while the N606U target line has degree two. The N608K rational local evaluator cannot extend directly to a global morphism into O_Eprime(2Oprime).",
        "next_requirement": "Search for a different correspondence whose H018 pullback has degree two, or construct and charge a genuine elementary modification with explicit normalized-Poincare charts. Do not repair the degree mismatch by pointwise weights.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608L degree audit failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
