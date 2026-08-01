#!/usr/bin/env sage -python
"""N608D: materialize and replay a degree-three genus-two elliptic component."""
import argparse
import datetime
import hashlib
import json
import subprocess
from pathlib import Path

from sage.all import EllipticCurve, GF, HyperellipticCurve, PolynomialRing, gcd


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_head():
    result = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True)
    return result.stdout.strip() if result.returncode == 0 else None


def transform_at_root(poly, root, degree, ring):
    variable = ring.gen()
    return sum(
        ring.base_ring()(coefficient) * variable ** (degree - index) * (root * variable + 1) ** index
        for index, coefficient in enumerate(poly.list())
    )


def monic(poly):
    return poly / poly.leading_coefficient()


def make_fixture():
    field = GF(101)
    x_ring = PolynomialRing(field, "X")
    x_var = x_ring.gen()
    t_ring = PolynomialRing(field, "t")
    t_var = t_ring.gen()
    a_value = field(0)
    b_value = field(4)
    factor = x_var ** 3 + a_value * x_var ** 2 + b_value * x_var + 1
    companion = 4 * x_var ** 3 + b_value ** 2 * x_var ** 2 + 2 * b_value * x_var + 1
    hyperelliptic = factor * companion
    numerator = x_var ** 3 - b_value * x_var - 2
    root = factor.roots()[0][0]
    discriminant_factor = 4 * a_value ** 3 + 27 - 18 * a_value * b_value - a_value ** 2 * b_value ** 2 + 4 * b_value ** 3
    leading = -discriminant_factor
    scale = leading.sqrt()
    factor_t = transform_at_root(factor, root, 3, t_ring)
    numerator_t = transform_at_root(numerator, root, 3, t_ring)
    hyperelliptic_t = transform_at_root(hyperelliptic, root, 6, t_ring)
    curve = HyperellipticCurve(hyperelliptic_t)
    quotient = EllipticCurve(
        field,
        [
            0,
            (12 * a_value ** 2 - 2 * a_value * b_value ** 2 - 18 * b_value) / leading,
            0,
            (-12 * a_value + b_value ** 2) / leading,
            field(4) / leading,
        ],
    )
    return {
        "field": field,
        "t_ring": t_ring,
        "t": t_var,
        "a": a_value,
        "b": b_value,
        "root": root,
        "factor": factor,
        "companion": companion,
        "hyperelliptic": hyperelliptic,
        "factor_t": factor_t,
        "numerator_t": numerator_t,
        "hyperelliptic_t": hyperelliptic_t,
        "curve": curve,
        "quotient": quotient,
        "scale": scale,
        "leading": leading,
    }


def regular_image(point, fixture):
    if not point[2]:
        return None
    field = fixture["field"]
    t_value = field(point[0])
    y_value = field(point[1])
    factor_t = fixture["factor_t"]
    if factor_t(t_value) == 0:
        return None
    root = fixture["root"]
    u_value = t_value * (root * t_value + 1) ** 2 / factor_t(t_value)
    w_value = y_value * fixture["numerator_t"](t_value) / (factor_t(t_value) ** 2 * fixture["scale"])
    return fixture["quotient"](u_value, w_value)


def reduce_mumford(u_poly, v_poly, fixture):
    curve_poly = fixture["hyperelliptic_t"]
    u_poly = monic(u_poly)
    v_poly = v_poly % u_poly
    reductions = 0
    while u_poly.degree() > 2:
        quotient, remainder = (curve_poly - v_poly ** 2).quo_rem(u_poly)
        if remainder:
            raise RuntimeError("Mumford reduction received a non-divisor pair")
        u_poly = monic(quotient)
        v_poly = (-v_poly) % u_poly
        reductions += 1
    return u_poly, v_poly, reductions


def raw_fiber_mumford(quotient_point, fixture):
    field = fixture["field"]
    t_var = fixture["t"]
    factor_t = fixture["factor_t"]
    numerator_t = fixture["numerator_t"]
    if quotient_point.is_zero():
        u_poly = factor_t
        v_poly = fixture["t_ring"].zero()
    else:
        u_value = field(quotient_point[0])
        w_value = field(quotient_point[1])
        u_poly = t_var * (fixture["root"] * t_var + 1) ** 2 - u_value * factor_t
        inverse = numerator_t.inverse_mod(u_poly)
        v_poly = (w_value * fixture["scale"] * factor_t ** 2 * inverse) % u_poly
    congruence_exact = ((v_poly ** 2 - fixture["hyperelliptic_t"]) % u_poly).is_zero()
    reduced_u, reduced_v, reductions = reduce_mumford(u_poly, v_poly, fixture)
    return {
        "raw_u": u_poly,
        "raw_v": v_poly,
        "congruence_exact": congruence_exact,
        "reduced_u": reduced_u,
        "reduced_v": reduced_v,
        "reductions": reductions,
    }


def run():
    fixture = make_fixture()
    field = fixture["field"]
    curve = fixture["curve"]
    quotient = fixture["quotient"]
    jacobian_points = curve.jacobian()(field)
    curve_points = list(curve.points())
    regular_images = []
    exceptional_points = []
    for point in curve_points:
        try:
            image = regular_image(point, fixture)
        except (TypeError, ValueError):
            image = None
        if image is None:
            exceptional_points.append(str(point))
        else:
            regular_images.append((point, image))

    raw_by_quotient = {}
    image_by_quotient = {}
    for point in quotient.points():
        raw = raw_fiber_mumford(point, fixture)
        raw_by_quotient[str(point)] = raw
        image_by_quotient[str(point)] = jacobian_points([raw["reduced_u"], raw["reduced_v"]])
    origin_key = str(quotient(0))
    origin_image = image_by_quotient[origin_key]
    pullback_by_quotient = {
        key: image - origin_image for key, image in image_by_quotient.items()
    }

    homomorphism_pairs = 0
    homomorphism_failures = 0
    for left in quotient.points():
        for right in quotient.points():
            homomorphism_pairs += 1
            if pullback_by_quotient[str(left + right)] != pullback_by_quotient[str(left)] + pullback_by_quotient[str(right)]:
                homomorphism_failures += 1
    image_frequencies = {}
    for _point, image in regular_images:
        image_frequencies[str(image)] = image_frequencies.get(str(image), 0) + 1
    pullback_examples = []
    for point in quotient.points()[:4]:
        raw = raw_by_quotient[str(point)]
        pullback_examples.append(
            {
                "quotient_point": str(point),
                "raw_u": str(raw["raw_u"]),
                "raw_v": str(raw["raw_v"]),
                "reduced_u": str(raw["reduced_u"]),
                "reduced_v": str(raw["reduced_v"]),
                "reductions": raw["reductions"],
            }
        )
    checks = {
        "smooth_genus_two_model": gcd(fixture["hyperelliptic_t"], fixture["hyperelliptic_t"].derivative()).degree() == 0,
        "prime_order_quotient": int(quotient.cardinality()) == 103,
        "all_raw_fibres_satisfy_curve_congruence": all(raw["congruence_exact"] for raw in raw_by_quotient.values()),
        "all_raw_fibres_reduce_to_mumford_degree_at_most_two": all(raw["reduced_u"].degree() <= 2 for raw in raw_by_quotient.values()),
        "pullback_homomorphism_exact": homomorphism_failures == 0,
        "pullback_image_has_full_quotient_size": len(set(pullback_by_quotient.values())) == int(quotient.cardinality()),
    }
    return {
        "schema": "ecdlp.genus2.degree3-elliptic-component.n608d.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "git_head": git_head(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "HYPOTHESIS / EXPLICIT-COVER-MATERIALIZATION / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "literature": {
            "source": "Tony Shaska, Genus 2 fields with degree 3 elliptic subfields, arXiv:math/0109155",
            "normal_form_equations": ["(2)", "(5)", "(7)"],
        },
        "fixture": {
            "field_order": int(field.order()),
            "normal_form_a": int(fixture["a"]),
            "normal_form_b": int(fixture["b"]),
            "weierstrass_root": int(fixture["root"]),
            "original_curve": f"Y^2 = ({fixture['factor']}) * ({fixture['companion']})",
            "odd_degree_curve": f"y^2 = {fixture['hyperelliptic_t']}",
            "quotient_curve": str(quotient),
            "quotient_order": int(quotient.cardinality()),
            "jacobian_order": int(jacobian_points.cardinality()),
        },
        "map_audit": {
            "curve_rational_point_count": len(curve_points),
            "regular_chart_rational_source_count": len(regular_images),
            "exceptional_chart_source_count": len(exceptional_points),
            "regular_chart_image_fibre_sizes": sorted(set(image_frequencies.values())),
            "all_fibres_have_degree_at_most_three": all(raw["raw_u"].degree() <= 3 for raw in raw_by_quotient.values()),
            "all_nonzero_regular_fibres_have_raw_degree_three": all(
                raw["raw_u"].degree() == 3 for key, raw in raw_by_quotient.items() if key != origin_key
            ),
            "pullback_homomorphism_pairs": homomorphism_pairs,
            "pullback_homomorphism_failures": homomorphism_failures,
            "pullback_image_size": len(set(pullback_by_quotient.values())),
            "examples": pullback_examples,
        },
        "checks": checks,
        "preflight_pass": all(checks.values()),
        "strongest_valid_statement": (
            "N608D materializes a degree-three elliptic quotient and exact Jacobian pullback on a public toy curve. "
            "The pullback is a verified homomorphism on the 103-point quotient, not a relation-density, source-inverse, "
            "target-descent, or ECDLP speedup result."
        ),
        "next_requirement": (
            "Define a target-independent reduced-divisor theta slice, compare triple-sum occupancy against matched random "
            "Jacobian subsets, and charge source recovery and same-form target descent before any promotion."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    output = run()
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608D materialization preflight failed")
    print(json.dumps({"checks": len(output["checks"]), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
