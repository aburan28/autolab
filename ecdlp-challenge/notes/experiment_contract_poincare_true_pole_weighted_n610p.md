# Experiment Contract: N610P weighted H018 true-pole transition

## Hypothesis

The true-pole chart has weighted local representative `u^7 zD`, and on
`z=r*u` it transitions exactly to the N610J crossing representative through a
formal-group unit.

## Parameters

- levels: `0`, `1`, and `17`
- directions: `r=2..21`; the tangent direction `r=1` belongs to a separate
  addition chart
- true-pole increment: `z=r*u`

## Success criterion

All directions have a regular weighted true-pole representative, unit
transition factor, and exact transition identity.

## Reproduction command

```bash
sage -python tools/poincare_true_pole_weighted_transition_n610p.py \
  --out notes/poincare_true_pole_weighted_transition_n610p.json
```
