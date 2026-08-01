# Experiment Contract: N609J deck-invariant determinantal section

## Hypothesis

The bordered determinant built from the nine deck-evaluation rows gives an
explicit, deck-invariant equation for the selected H018 pencil on the regular
chart. Its `Q=O` formal limit is a fixed nonzero scalar times the N609C
regularized boundary section.

## Construction

For a lift `R` of `Q`, let `V_Q` have rows given by the moving Cauchy basis at
`phi(R+iT)`, and let `v_t` have entries `x(R+iT)-t`. For a test point P define

```text
D_t(P,Q) = det([[V_Q, v_t], [b(P,Q), 0]]).
```

On the regular chart, `D_t=-det(V_Q)*(lambda(P,Q)-t)`. Deck translation cyclically
permutes nine rows, an even permutation, so the determinant descends from the
chosen lift.

## Parameters

- H018 fixture over `F_103`.
- Deck degree nine and selected levels `t=0,1,17`.
- A rank-nine finite rational P interpolation set plus six held-out finite
  rational P values for each formal Q-origin boundary check.

## Metrics

- Schur-complement identity on a nonzero rational Q fibre;
- deck-invariance under one deck shift;
- determinant valuation and regularity at Q origin;
- equality of boundary constants with one scalar times N609C's section.

## Controls

- The leading 9x9 determinant must be nonzero at the selected nonzero Q.
- The determinant must change if the selected x-values are mutated.
- Boundary equality must hold at zeros and nonzeros, not only after division.

## Success criterion

The determinant is deck invariant, satisfies the exact open-chart equation,
and has a regular Q-origin limit whose degree-nine model is proportional to
the regularized boundary section on the interpolation and held-out points.

## Falsification criterion

Any failure leaves the determinant as a lift-dependent or unnormalized chart
expression and blocks its use as an explicit global-section candidate.

## Boundary

This constructs an equation on the declared regular chart and a formal origin
limit. It does not yet prove a globally regular Cartier section across the true
pole or non-rational loci, nor smoothness, normalization, a Jacobian, relations,
rank, descent, cost, or ECDLP improvement.

## Reproduction

```bash
sage -python tools/poincare_determinantal_section_n609j.py --out notes/poincare_determinantal_section_n609j.json
```
