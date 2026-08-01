# Experiment Contract: N609A degree-nine residual H018 source inverse

## Hypothesis

After N608Z removes the common cleared base graph `P=4Q`, the residual
degree-nine polynomial can replace the N608U degree-ten eliminant while
recovering exactly the same declared open-level sources.

## Parameters

- Registered `F_103` H018 fixture.
- Levels `t=0,1,17` and all 108 nonzero base points.
- N608Z quotient `F_Q,t(x)/(x-x(4Q))`.
- N608U exhaustive open-surface source enumeration.

## Metrics

- residual polynomial degree;
- root-factor calls;
- recovered and exhaustive open-source counts;
- complete set equality;
- residual intersections with the removed base graph.

## Controls

- Positive: each residual inverse set equals exhaustive open enumeration.
- Negative: reversing the Cauchy numerator sign eliminates the required
  `x-x(4Q)` divisibility, so the residual construction must reject it.

## Boundary

This is a degree reduction in one toy source-inversion subroutine. It does not
provide a factor base, relation law, rank, target descent, fully charged cost
comparison, or ECDLP improvement.

## Reproduction

```bash
sage -python tools/poincare_residual_inverse_n609a.py --out notes/poincare_residual_inverse_n609a.json
sage -python tools/verify_poincare_residual_inverse_n609a.py --primary notes/poincare_residual_inverse_n609a.json --out notes/poincare_residual_inverse_n609a_structural_replay.json
```
