#!/usr/bin/env sage -python
"""N611P: materialize cubic H018 principal-theta quotient kernels on E[2]^2."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF


P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def point_text(point):
    return "O" if point.is_zero() else "({}, {})".format(point[0], point[1])


def frobenius(point):
    return point.curve()(point[0] ** P, point[1] ** P) if not point.is_zero() else point


def pi_action(point):
    return frobenius(point)


def one_plus_pi_action(point):
    return point + pi_action(point)


def h018_mod_two(vector):
    left, right = vector
    return left + pi_action(right), one_plus_pi_action(left) + right


def line(generator):
    origin = generator[0].curve()(0)
    return (origin, origin), generator


def line_text(values):
    return [[point_text(left), point_text(right)] for left, right in values]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    field = GF(P**3, name="a")
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    roots = sorted(root for root, _multiplicity in (curve.base_field()['X'].gen()**3 + curve.base_field()['X'].gen() + 24).roots(field))
    if len(roots) != 3:
        raise RuntimeError("two-torsion cubic did not split into three roots")
    r = curve(roots[0], 0)
    r1, r2 = frobenius(r), frobenius(frobenius(r))
    g = (r1, r)
    h = (r2, r1)
    gh = (r1 + r2, r + r1)
    lines = {"G_g": line(g), "G_h": line(h), "G_gh": line(gh)}
    frobenius_lines = {
        name: next(target for target, values in lines.items() if line((frobenius(generator[0]), frobenius(generator[1]))) == values)
        for name, (_origin, generator) in lines.items()
    }
    kernel = ((curve(0), curve(0)), g, h, gh)
    identity_graph = (r, r)
    f3_graph = (r2, r)
    gates = {
        "two_torsion_cubic_splits_only_over_cubic": len(roots) == 3 and all(root ** P != root for root in roots),
        "frobenius_orbit_is_exact": r1 != r and r2 != r and frobenius(r2) == r,
        "three_generators_are_distinct_nonzero_order_two": len({g, h, gh}) == 3 and all(2 * component == curve(0) for vector in (g, h, gh) for component in vector) and all(vector != (curve(0), curve(0)) for vector in (g, h, gh)),
        "h018_mod_two_annihilates_full_kernel": all(h018_mod_two(vector) == (curve(0), curve(0)) for vector in kernel),
        "frobenius_cycles_quotient_lines": frobenius_lines == {"G_g": "G_h", "G_h": "G_gh", "G_gh": "G_g"},
        "frobenius_cubed_fixes_each_line": all(line((frobenius(frobenius(frobenius(generator[0]))), frobenius(frobenius(frobenius(generator[1]))))) == values for values in lines.values() for _origin, generator in [values]),
        "identity_and_f3_graph_controls_differ_from_h018_generator": g != identity_graph and g != f3_graph,
    }
    output = {
        "schema": "ecdlp.h018.cubic-principal-theta-kernel.n611p.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / CUBIC_H018_PRINCIPAL_THETA_KERNEL_COORDINATES / MODEL-BOUND / TOY-EVIDENCE / QUOTIENT_EQUATIONS_OPEN / NO_ECDLP_CLAIM",
        "assumptions": [
            "N606R's polarization-descent route is used only as a construction hypothesis.",
            "Frobenius on E[2] realizes multiplication by pi in the N606R F4 model.",
        ],
        "records": {
            "curve": "y^2=x^3+x+24 over F_(103^3)",
            "two_torsion_roots": [str(root) for root in roots],
            "frobenius_orbit": [point_text(point) for point in (r, r1, r2)],
            "generators": {"g": [point_text(value) for value in g], "h": [point_text(value) for value in h], "g_plus_h": [point_text(value) for value in gh]},
            "lines": {name: line_text(values) for name, values in lines.items()},
            "frobenius_line_permutation": frobenius_lines,
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The three N606R maximal isotropic candidate lines are now exact cyclic order-two subgroups of E[2]^2 over F_(103^3), with the declared Frobenius orbit and H018 mod-two kernel equation verified directly on coordinates.",
        "next_requirement": "Construct the quotient A/G_g over F_(103^3), its principal polarization and theta divisor, then compute the Frobenius-semilinear action on all three pullbacks before claiming scalar H018 sections.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611P cubic theta kernel preflight failed")
    print(json.dumps({"output": str(args.out), "checks": sum(gates.values())}, sort_keys=True))


if __name__ == "__main__":
    main()
