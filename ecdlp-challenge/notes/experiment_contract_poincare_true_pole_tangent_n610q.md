# Experiment Contract: N610Q H018 true-pole tangent transition

## Hypothesis

The elliptic-curve doubling branch supplies the same weighted transition as
N610P at the omitted tangent direction `r=1`.

## Success criterion

For each selected level, the weighted true-pole representative is regular, the
transition factor is a unit, and the exact transition identity holds.

## Reproduction command

```bash
sage -python tools/poincare_true_pole_tangent_transition_n610q.py \
  --out notes/poincare_true_pole_tangent_transition_n610q.json
```
