# Experiment Contract: N610O H018 true-pole moving crossing

## Hypothesis

In the exact chord chart `P=-4Q+R`, multiplying the determinant by the formal
increment parameter removes the true-pole singularity and exposes its Q-origin
weight.

## Success criterion

All selected levels have nonnegative increment valuation after the declared
factor, with finite inspected Q-origin valuations.

## Reproduction command

```bash
sage -python tools/poincare_true_pole_crossing_n610o.py \
  --out notes/poincare_true_pole_crossing_n610o.json
```
