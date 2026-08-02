# Result: N611J Slope-Aligned Graph Lattice

## Status

`NEGATIVE COMPACT_INTEGRAL_SLOPE_ALIGNED_GRAPH_CORRECTION / MODEL-BOUND /
TOY-EVIDENCE / MATRIX_VALUED_OR_OUTSIDE_LATTICE_OPEN / NO_ECDLP_CLAIM`

## Exact Result

For the N608M integral correspondence family

```text
h_(a,b) = ([a]phi + [b]pi, pi),
```

N611I supplies `pi=[50]phi` on the rational cover orbit. Requiring the
N606Y transport law gives the exact congruence

```text
84a + 58b = 50 mod 109.
```

The degree-nine deck kernel is separated precisely when `a` is a unit modulo
`9`. Under those two constraints, the exact N608M pullback-degree form

```text
99a^2 + 27ab + 81b^2 - 195a - 9b + 99
```

has unique minimum

```text
(a,b)=(-2,0),   degree=885.
```

The map `chi=[-2]phi` matches the Poincare transport on all 108 nonzero
rational cover translations and has rank nine on all 108 nonzero rational
deck fibres. In contrast, the compact N608O graph `(a,b)=(1,0)` has degree
three and rank nine but matches transport on `0/108` translations.

The transport-aligned factoring control `(a,b)=(0,61)` passes `108/108`
translation checks but has rank one on every fibre and degree `300951`.

## Proof Bound

The positive quadratic part of the degree form is at least
`72(a^2+b^2)`, while the linear part is bounded below by
`-196 sqrt(a^2+b^2)`. Any solution of degree at most `885` must therefore
have `sqrt(a^2+b^2)<6`, reducing the minimum proof to the exact finite search
over `[-5,5]^2`. Its only aligned deck-separating candidate is `(-2,0)`.

## Evidence

- Producer: `5/5` gates pass.
- Independent replay: `6/6` checks pass.
- The rank panels use the full 108 nonzero base fibres and `F_(103^6)` deck
  points, not a sampled interpolation subset.

## Interpretation

Direct integral correction can reconcile translation slope and deck separation,
but not while retaining the compact degree-three graph class that underpins the
N608O elementary modification and two-section construction. This is a scoped
negative result for the direct integral graph lattice, not an impossibility
claim about non-scalar Poincare transitions or external correspondences.

## Next Concrete Action

Search for a genuinely matrix-valued correction between the N607E Cauchy frame
and N608G deck frame, or for a correspondence outside the integral
`([a]phi+[b]pi,pi)` lattice. Any successor must state a low-degree global
class, literal overlap data, full rank, source return, relation rank, target
descent, and charged cost before an ECDLP claim.

## Reproduction

```bash
sage -python tools/poincare_slope_aligned_graph_lattice_n611j.py \
  --out notes/poincare_slope_aligned_graph_lattice_n611j.json
sage -python tools/verify_poincare_slope_aligned_graph_lattice_n611j.py \
  --primary notes/poincare_slope_aligned_graph_lattice_n611j.json \
  --out notes/poincare_slope_aligned_graph_lattice_n611j_independent_verifier.json
```
