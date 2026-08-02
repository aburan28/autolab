#!/usr/bin/env sage -python
"""N611V: audit whether the reduced cubic cross map selects a 2-gluing."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from itertools import combinations
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


P = 103
ROOT = Path("notes")
N611U = ROOT / "poincare_cubic_reduced_polarization_frame_n611u.json"
N611U_VERIFY = ROOT / "poincare_cubic_reduced_polarization_frame_n611u_independent_verifier.json"


def sha256(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path): return json.loads(path.read_text(encoding="ascii"))
def receipt_matches(receipt, path): return receipt.get("verified") is True and receipt.get("primary_sha256") == sha256(path)
def frobenius(point): return point if point.is_zero() else point.curve()(point[0] ** P, point[1] ** P)


def bit_vector(value): return tuple((value >> shift) & 1 for shift in range(4))
def add_bits(left, right): return tuple(a ^ b for a, b in zip(left, right))
def bits_to_int(value): return sum(bit << shift for shift, bit in enumerate(value))
def span(left, right):
    zero = (0, 0, 0, 0)
    return frozenset(bits_to_int(value) for value in (zero, left, right, add_bits(left, right)))


def choose_symplectic_basis(points):
    nonzero = [point for point in points if not point.is_zero()]
    for left, right in combinations(nonzero, 2):
        if left.weil_pairing(right, 2) != 1:
            return left, right
    raise RuntimeError("no symplectic torsion basis found")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    n611u, n611u_verify = load(N611U), load(N611U_VERIFY)
    field = GF(P**6, name="a")
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    ring = PolynomialRing(field, "T")
    roots = sorted(root for root, _ in (ring.gen()**3 + ring.gen() + 24).roots())
    r = curve(roots[0], 0)
    phi, phi_dual = curve.isogeny(r), curve.isogeny(r).dual()
    target = phi.codomain()
    e2, g2 = curve(0).division_points(2), target(0).division_points(2)
    e_first, e_second = choose_symplectic_basis(e2)
    g_first, g_second = choose_symplectic_basis(g2)
    e_lookup = {(0, 0): curve(0), (1, 0): e_first, (0, 1): e_second, (1, 1): e_first + e_second}
    g_lookup = {(0, 0): target(0), (1, 0): g_first, (0, 1): g_second, (1, 1): g_first + g_second}

    def beta(point): return phi(12 * point + 4 * frobenius(point))
    def beta_dual(point): return -8 * phi_dual(point) - 4 * frobenius(phi_dual(point))
    def polarization(point):
        left, right = point
        return 9 * left + beta_dual(right), beta(left) + 345 * right
    def pairing(left, right):
        image_left, image_right = polarization(right)
        return left[0].weil_pairing(image_left, 2) * left[1].weil_pairing(image_right, 2)

    basis_vectors = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
    basis_points = [(e_lookup[value[:2]], g_lookup[value[2:]]) for value in basis_vectors]
    pairing_matrix = [[int(pairing(left, right) != 1) for right in basis_points] for left in basis_points]
    planes = set()
    for left, right in combinations([bit_vector(value) for value in range(1, 16)], 2):
        if left != right:
            planes.add(span(left, right))
    def dot(left, right):
        return sum(left[row] * pairing_matrix[row][column] * right[column] for row in range(4) for column in range(4)) % 2
    isotropic = [plane for plane in planes if all(dot(bit_vector(left), bit_vector(right)) == 0 for left in plane for right in plane)]
    coordinate_e = frozenset((0, 1, 2, 3))
    coordinate_g = frozenset((0, 4, 8, 12))
    graph_planes = [plane for plane in isotropic if plane.intersection(coordinate_e) == {0} and plane.intersection(coordinate_g) == {0}]
    beta_nonzero_control = any(not phi(point).is_zero() for point in e2)
    gates = {
        "n611u_reduced_frame_receipt_binds": n611u["preflight_pass"] is True and receipt_matches(n611u_verify, N611U),
        "degree_six_contains_full_two_torsion": len(e2) == 4 and len(g2) == 4,
        "reduced_cross_map_vanishes_on_e_two_torsion": all(beta(point).is_zero() for point in e2),
        "reduced_dual_cross_map_vanishes_on_eg_two_torsion": all(beta_dual(point).is_zero() for point in g2),
        "polarization_pairing_is_block_product": pairing_matrix == [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
        "all_fifteen_maximal_isotropic_planes_replay": len(isotropic) == 15,
        "six_transverse_graph_planes_replay": len(graph_planes) == 6,
        "unsuppressed_phi_positive_control_is_nonzero": beta_nonzero_control,
    }
    output = {
        "schema": "ecdlp.h018.reduced-level2-gluing-audit.n611v.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__)),
        "claim_status": "NEGATIVE RESULT / REDUCED_CROSS_MAP_DOES_NOT_SELECT_LEVEL2_FREY_KANI_GRAPH / MODEL-BOUND / TOY-EVIDENCE / HIGHER_LEVEL_OR_QUOTIENT_LATTICE_OPEN / NO_ECDLP_CLAIM",
        "assumptions": ["The Weil pairing computed from the reduced principal polarization represents the level-two commutator form.", "The result concerns only direct selection by N611U's cross map on full scalar two-torsion."],
        "bound_inputs": {str(path): sha256(path) for path in (N611U, N611U_VERIFY)},
        "records": {"field": "F_(103^6)", "e_two_torsion_size": len(e2), "eg_two_torsion_size": len(g2), "pairing_matrix_bits": pairing_matrix, "maximal_isotropic_plane_count": len(isotropic), "transverse_graph_plane_count": len(graph_planes), "cross_map": "phi_g o ([12]+[4]Frob)", "cross_map_control": "phi_g"},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "Over F_(103^6), the reduced cross map beta and beta^dagger both vanish on full two-torsion, and M induces the block-product symplectic pairing. All fifteen maximal isotropic planes and six transverse graph planes exist, but beta selects none of them. A direct level-two Frey-Kani theta construction is therefore not supplied by the N611U cross term.",
        "next_requirement": "Search a full quotient-lattice transformation or higher-level theta structure that selects a product presentation, then construct and evaluate the resulting genus-two theta curve before relation work.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N611V level-two gluing audit failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__": main()
