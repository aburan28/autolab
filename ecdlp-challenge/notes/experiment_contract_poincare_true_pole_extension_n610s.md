# Experiment Contract: N610S extension-field H018 true-pole directions

## Hypothesis

The N610P weighted transition is not restricted to base-field direction
parameters.

## Parameters

- coefficient field: `F_(103^6)`
- directions: first sixteen nonexceptional powers of the field generator
- levels: `0`, `1`, and `17`

## Success criterion

Every extension-field row has a regular weighted representative, unit
transition factor, and exact transition identity.

## Reproduction command

```bash
sage -python tools/poincare_true_pole_extension_direction_n610s.py \
  --out notes/poincare_true_pole_extension_direction_n610s.json
```
