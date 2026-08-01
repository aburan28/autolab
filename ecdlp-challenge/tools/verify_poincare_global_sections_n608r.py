#!/usr/bin/env sage -python
"""Independent N608R verifier for formal global-section selection."""
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


def translate(point, x_value, y_value, ring):
    x0, y0 = ring(point[0]), ring(point[1])
    slope = (y_value - y0) / (x_value - x0)
    x_sum = slope**2 - x_value - x0
    return x_sum, -y0 + slope * (x0 - x_sum)


def image(isogeny, x_value, y_value, ring):
    x_map, y_map = isogeny.rational_maps()
    return ring(x_map(x_value, y_value)), ring(y_map(x_value, y_value))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi = (second * first).dual()
    cover = pi.domain()
    raw_phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    phi = raw_phi.codomain().isomorphism_to(target) * raw_phi
    field = GF(P**6, name="b")
    ring = LaurentSeriesRing(field, "t", default_prec=PRECISION)
    x0, y0 = ring(cover.formal_group().x(PRECISION)), ring(cover.formal_group().y(PRECISION))
    deck = None
    cover_k = cover.base_extend(field)
    for root, _ in pi.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            deck = cover_k(x, rhs.sqrt())
            break
    if deck is None:
        raise RuntimeError("could not recover deck generator")
    phi0, pi0 = image(phi, x0, y0, ring), image(pi, x0, y0, ring)
    u = target.formal_group().mult_by_n(-4, PRECISION)(-pi0[0] / pi0[1])
    xu, yu = ring(target.formal_group().x(PRECISION)(u)), ring(target.formal_group().y(PRECISION)(u))
    correction_data = {
        -1: [1, 0, 0, 0, 0, 0, 0, 0], 0: [0] * 8, 1: [0, 1, 0, 0, 0, 0, 0, 0],
        2: [0, 0, -1, 0, 0, 0, 0, 0], 3: [0, 0, 0, 1, 0, 0, 0, 0],
        4: [0, 0, 0, 0, -1, 0, 0, 0], 5: [0, 1, 0, 0, 0, 1, 0, 0],
        6: [0, 0, -1, 0, 0, 0, -1, 0], 7: [0, 24, 0, 2, 0, 0, 0, 1],
    }
    corrections = [sum(ring(correction_data[e][i]) * u**e for e in range(-1, 8)) for i in range(8)]
    rows, values = [], {"1": [], "x": [], "y": []}
    for index in range(9):
        if index == 0:
            xc, yc = x0, y0
            xp, yp = phi0
        else:
            xc, yc = translate(index * deck, x0, y0, ring)
            xp, yp = image(phi, xc, yc, ring)
        rows.append([1, xp, yp, xp**2, xp * yp, xp**3, xp**2 * yp, xp**4, (yp + yu) / (xp - xu)])
        values["1"].append(ring.one())
        values["x"].append(xc)
        values["y"].append(yc)
    matrix = Matrix(ring, rows)
    observed = {}
    for name, raw in values.items():
        moving = matrix.solve_right(vector(ring, raw))
        regular = [moving[i] + moving[8] * corrections[i] for i in range(8)] + [moving[8] * u**8]
        observed[name] = [None if item.is_zero() else int(item.valuation()) for item in regular]
    primary_data = primary["records"]["section_origin_data"]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "regular_valuations_match": all(observed[name] == primary_data[name]["regular_frame_valuations"] for name in ("1", "x", "y")),
        "one_and_x_regular": all(value is None or value >= 0 for name in ("1", "x") for value in observed[name]),
        "y_has_simple_pole": min(value for value in observed["y"] if value is not None) == -1,
        "selected_subspace_retained": primary["records"]["selected_cover_subspace"] == "span{1,x}",
        "determinant_is_nonzero": matrix.det() != 0,
    }
    output = {
        "schema": "ecdlp.h018.poincare-global-sections.n608r.independent-verifier.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N608R independent verification failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
