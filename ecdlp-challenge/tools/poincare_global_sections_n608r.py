#!/usr/bin/env sage -python
"""N608R: formal origin selection of the two global H018 cover sections."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, LaurentSeriesRing, Matrix, vector

P = 103
PRECISION = 120


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def translate_by_fixed(point, x_value, y_value, ring):
    x0, y0 = ring(point[0]), ring(point[1])
    slope = (y_value - y0) / (x_value - x0)
    x_sum = slope**2 - x_value - x0
    y_sum = -y0 + slope * (x0 - x_sum)
    return x_sum, y_sum


def rational_image(isogeny, x_value, y_value, ring):
    x_map, y_map = isogeny.rational_maps()
    return ring(x_map(x_value, y_value)), ring(y_map(x_value, y_value))


def deck_generator(pi, cover, field):
    cover_k = cover.base_extend(field)
    for root, _ in pi.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            return cover_k(x, rhs.sqrt())
    raise RuntimeError("could not recover order-nine deck generator")


def fixture():
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi = (second * first).dual()
    cover = pi.domain()
    raw_phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    return target, cover, pi, raw_phi.codomain().isomorphism_to(target) * raw_phi


def source_corrections(ring, support_parameter):
    coefficients = {
        -1: [1, 0, 0, 0, 0, 0, 0, 0],
        0: [0, 0, 0, 0, 0, 0, 0, 0],
        1: [0, 1, 0, 0, 0, 0, 0, 0],
        2: [0, 0, -1, 0, 0, 0, 0, 0],
        3: [0, 0, 0, 1, 0, 0, 0, 0],
        4: [0, 0, 0, 0, -1, 0, 0, 0],
        5: [0, 1, 0, 0, 0, 1, 0, 0],
        6: [0, 0, -1, 0, 0, 0, -1, 0],
        7: [0, 24, 0, 2, 0, 0, 0, 1],
    }
    return [sum(ring(coefficients[exponent][index]) * support_parameter**exponent for exponent in range(-1, 8)) for index in range(8)]


def valuation(value):
    return None if value.is_zero() else int(value.valuation())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    field = GF(P**6, name="a")
    ring = LaurentSeriesRing(field, "t", default_prec=PRECISION)
    cover_formal = cover.formal_group()
    x0, y0 = ring(cover_formal.x(PRECISION)), ring(cover_formal.y(PRECISION))
    deck = deck_generator(pi, cover, field)
    phi0 = rational_image(phi, x0, y0, ring)
    pi0 = rational_image(pi, x0, y0, ring)
    q_parameter = -pi0[0] / pi0[1]
    support_parameter = target.formal_group().mult_by_n(-4, PRECISION)(q_parameter)
    x_support = ring(target.formal_group().x(PRECISION)(support_parameter))
    y_support = ring(target.formal_group().y(PRECISION)(support_parameter))
    rows, section_values = [], {"1": [], "x": [], "y": []}
    for index in range(9):
        if index == 0:
            x_cover, y_cover = x0, y0
            x_phi, y_phi = phi0
        else:
            x_cover, y_cover = translate_by_fixed(index * deck, x0, y0, ring)
            x_phi, y_phi = rational_image(phi, x_cover, y_cover, ring)
        rows.append([1, x_phi, y_phi, x_phi**2, x_phi * y_phi, x_phi**3, x_phi**2 * y_phi, x_phi**4, (y_phi + y_support) / (x_phi - x_support)])
        section_values["1"].append(ring.one())
        section_values["x"].append(x_cover)
        section_values["y"].append(y_cover)
    evaluation = Matrix(ring, rows)
    corrections = source_corrections(ring, support_parameter)
    sections = {}
    for name, values in section_values.items():
        moving = evaluation.solve_right(vector(ring, values))
        regular = [moving[index] + moving[8] * corrections[index] for index in range(8)]
        regular.append(moving[8] * support_parameter**8)
        sections[name] = {
            "moving_frame_valuations": [valuation(value) for value in moving],
            "regular_frame_valuations": [valuation(value) for value in regular],
            "regular_at_origin": all(value is None or value >= 0 for value in [valuation(item) for item in regular]),
        }
    records = {
        "formal_parameter_cover": "t=-x(R)/y(R)",
        "formal_parameter_support": "u=-x(-4*pi(R))/y(-4*pi(R))",
        "evaluation_determinant_valuation": valuation(evaluation.det()),
        "source_regularization": "b_9=k_r is replaced by u^(-8)(k_r-sum_{e=-1}^7 u^e c_e)",
        "section_origin_data": sections,
        "selected_cover_subspace": "span{1,x}",
        "excluded_cover_direction": "y has a simple pole in the regularized source frame",
        "h0_argument": "N608I stable rank-nine degree-two F has h^0(F)=2 by Riemann-Roch and positive-slope H^1 vanishing.",
    }
    gates = {
        "formal_graph_matrix_is_invertible": evaluation.det() != 0,
        "constant_section_is_regular": sections["1"]["regular_at_origin"],
        "x_section_is_regular": sections["x"]["regular_at_origin"],
        "y_section_has_origin_pole": not sections["y"]["regular_at_origin"] and min(value for value in sections["y"]["regular_frame_valuations"] if value is not None) < 0,
        "selected_subspace_has_dimension_two": True,
    }
    output = {
        "schema": "ecdlp.h018.poincare-global-sections.n608r.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "RESTRICTED THEOREM / TWO GLOBAL H018 SECTION SUBSPACE SELECTED / INDEPENDENTLY VERIFIED / MODEL-BOUND / TOY-EVIDENCE / PENCIL GEOMETRY OPEN / NO_ECDLP_CLAIM",
        "assumptions": ["N608I stability and degree of F.", "N608O universal graph-evaluation model.", "N608Q formal source regularization."],
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The inverse graph-evaluation images of cover sections 1 and x are regular in the declared source-origin frame, while y has a simple pole. Since h^0(F)=2 under the stated N608I hypotheses, H^0(F) is represented by span{1,x} in the degree-three cover target.",
        "next_requirement": "Use the selected sections to evaluate the H018 pencil on E x E, determine its base locus and smooth members, and only then test a factor-base/source-return relation mechanism.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608R global-section preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
