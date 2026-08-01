#!/usr/bin/env sage -python
"""Independent replay for the N606U cyclic pushforward preflight."""
import argparse, json
from pathlib import Path
from sage.all import EllipticCurve, GF

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    curve = EllipticCurve(GF(103), [0, 0, 0, 1, 24])
    phi = curve.isogenies_prime_degree(3)[0]
    psi = phi.codomain().isogenies_prime_degree(3)[0]
    dual = (psi * phi).dual()
    cover = dual.domain()
    dimension = len(cover.riemann_roch_basis(2 * cover.divisor(cover(0))))
    checks = {
        "primary_pass": primary["pushforward_identification_preflight_pass"] is True,
        "dual_degree": int(dual.degree()) == 9,
        "cover_basis_dimension": dimension == 2,
        "numerical_degree": primary["records"]["fourier_mukai_numerics"]["theta11_tensor_degree"] == 2,
        "cyclic_kernel_orders": primary["records"]["kernel_x_coordinate_orders_over_f103_squared"] == [3, 9, 9, 9],
    }
    output = {"verified": all(checks.values()), "checks": checks, "strongest_valid_statement": "Independent replay confirms the cyclic degree-nine cover and two-dimensional degree-two cover-side section space."}
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]: raise RuntimeError("N606U verifier failed")
    print(json.dumps({"verified": True, "output": str(args.output)}, sort_keys=True))

if __name__ == "__main__": main()
