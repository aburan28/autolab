# Result: N611W Reduced 2-Adic Cross-Map Structure

## Status

`RESTRICTED THEOREM / NEGATIVE_DIRECT_POWER_OF_TWO_FREY_KANI_GRAPH_BY_REDUCED_CROSS_MAP /
MODEL-BOUND / TOY-EVIDENCE / QUOTIENT_LATTICE_AND_HIGHER_THETA_OPEN /
NO_ECDLP_CLAIM`

## Exact Result

The N611U cross map factors as

```text
beta = phi_g o [4] o ([3]+Frob).
```

The CM norm `Norm(3+pi)=97` is odd, so the final factor is a two-adic unit.
The degree-two separable isogeny has two-adic factors `(1,2)`. Therefore
`beta` has factors `(4,8)`.

The full torsion computations agree exactly:

| Level | Field degree | `|E[2^m]|` | `|beta(E[2^m])|` |
|---|---:|---:|---:|
| 4 | 6 | 16 | 1 |
| 8 | 12 | 64 | 2 |
| 16 | 24 | 256 | 8 |

Thus the first nonzero cross-map residue occurs at level eight, but it has
image order two. At level sixteen the residue grows to order eight, still far
below the full torsion group. Since the determinant is never a two-adic unit,
`beta` cannot be an isomorphism on any `E[2^m]` and cannot itself define a
direct Frey-Kani graph.

## Evidence

- N611V binding and the degree-six scalar-torsion prerequisite replay.
- Frobenius orders modulo `4`, `8`, and `16` are exactly `6`, `12`, and `24`.
- Each extension contains the full expected torsion group.
- Producer passes all six arithmetic and image-size gates.
- The independent verifier separately rechecks factorization, torsion image
  sizes, Frobenius orders, and nonisomorphism.

## Interpretation

The higher-level cross map contains real, nonzero torsion residue, but it is
not a direct gluing isomorphism. This closes only the attempt to obtain the
theta curve by taking the graph of `beta` at a power-of-two level. It does not
rule out an auxiliary quotient lattice, a theta-group construction using the
residue, an explicit genus-two curve, a theta evaluator, semilinear descent,
relations, rank, target descent, a cost advantage, or an ECDLP improvement.

## Next Concrete Action

Construct an explicit quotient-lattice or theta-group model that absorbs the
`(4,8)` elementary divisors before asking for a product presentation. It must
then produce a named genus-two/theta object with a source evaluator before
relation collection is meaningful.

## Reproduction

```bash
sage -python tools/poincare_reduced_2adic_crossmap_n611w.py \
  --out notes/poincare_reduced_2adic_crossmap_n611w.json
sage -python tools/verify_poincare_reduced_2adic_crossmap_n611w.py \
  --primary notes/poincare_reduced_2adic_crossmap_n611w.json \
  --out notes/poincare_reduced_2adic_crossmap_n611w_independent_verifier.json
```
