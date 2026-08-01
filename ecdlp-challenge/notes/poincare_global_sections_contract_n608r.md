# Experiment Contract: N608R select the global H018 two-section subspace

## Hypothesis

In the N608Q regular source-origin frame, precisely two of the three degree
three cover sections `1`, `x`, and `y` extend through the origin under inverse
graph evaluation. Those two sections are the global H018 section space because
the N608I stable rank-nine degree-two bundle has `h^0=2`.

## Exact Test

Along the formal cover branch through the origin, build the graph evaluation
matrix on the nine deck sheets. Invert the target value vectors for `1`, `x`,
and `y`, then change the resulting coefficients to the N608Q regularized
source frame. Record every Laurent valuation.

## Controls

- `1` is the constant target-section control.
- `x` is admitted only if all regular-frame coefficients have nonnegative
  valuation.
- `y` is the negative control and must retain a negative valuation.

## Success Criterion

Exactly `1` and `x` are regular. The result selects `H^0(F)=span{1,x}` in the
degree-three cover representation, but does not yet identify the full
nine-dimensional origin quotient direction.

## Boundary

This is a two-section construction on the registered toy fixture. It does not
yet establish a smooth pencil member, factor base, source recovery, relation
rank, target descent, cost advantage, or ECDLP speedup.

## Reproduction

```bash
sage -python tools/poincare_global_sections_n608r.py --out notes/poincare_global_sections_n608r.json
sage -python tools/verify_poincare_global_sections_n608r.py --primary notes/poincare_global_sections_n608r.json --out notes/poincare_global_sections_n608r_independent_verifier.json
```
