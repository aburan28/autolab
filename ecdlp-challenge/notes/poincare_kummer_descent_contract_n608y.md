# Experiment Contract: N608Y diagonal product-Kummer descent admission

## Hypothesis

The complete `t=0` H018 pencil level may be invariant under simultaneous
negation `(P,Q) -> (-P,-Q)`, and its quotient image in coordinates

```text
(u,x,w) = (x(P), x(Q), y(P)y(Q))
```

may admit a lower-than-random finite interpolation equation. Such an equation
would be an initial, not sufficient, ingredient of a product-Kummer model for
the normalized level curve.

## Null Hypothesis

Diagonal sign symmetry may hold, but the first finite interpolation relation is
forced merely by having 42 quotient samples and is indistinguishable from
matched random 42-point subsets of the same Kummer locus.

## Parameters

- Registered N608U `t=0` complete source set.
- Diagonal product-Kummer coordinates on nonidentity pairs.
- Monomials `u^i x^j` and `w u^i x^j`, scanned by increasing `i+j`.
- Sixteen deterministic random 42-point controls from the same rational
  product-Kummer locus, seed `0x608`.

## Metrics

- closure under diagonal, left, and right sign involutions;
- quotient source count;
- first monomial box with a nonzero interpolation kernel;
- rank and kernel dimension of the first box;
- matched-control first-box distribution.

## Controls

- Positive: every diagonal sign image remains in the complete source set.
- Negative: left and right sign images must not be silently accepted as the
  same symmetry.
- Matched random controls use the same quotient surface, point count, and
  monomial scan.

## Success Criterion

Admission to a finite-model construction requires diagonal closure and a first
interpolation box strictly smaller than every matched random control, or a
non-generic rank defect within the same first box.

## Boundary

A rational-point interpolation relation is not a global defining equation. No
normalization, factor base, relation law, rank, descent, cost, or ECDLP claim
follows unless the admission criterion passes.

## Reproduction

```bash
sage -python tools/poincare_kummer_descent_n608y.py --out notes/poincare_kummer_descent_n608y.json
sage -python tools/verify_poincare_kummer_descent_n608y.py --primary notes/poincare_kummer_descent_n608y.json --out notes/poincare_kummer_descent_n608y_structural_replay.json
```
