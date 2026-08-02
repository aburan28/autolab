# Experiment Contract: N610T generic-direction H018 true-pole transition

## Hypothesis

Over the rational-function direction field `F_(103^6)(r)`, the weighted
true-pole representative and transition unit are regular in the generic
formal direction chart.

## Parameters

- coefficient field: `F_(103^6)(r)`
- formal precision: `32`
- levels: `0`, `1`, and `17`
- exceptional direction charts retained separately: `r=0`, `r=1`, `r=-1`

## Metrics

- weighted true-pole valuation
- transition-unit valuation
- crossing valuation
- truncated-formal exact transition identity

## Positive control

N610P/N610Q/N610R/N610S give regular sampled nonexceptional directions.

## Negative control

Specializing the generic chord formula at `r=1` must be rejected because the
denominator becomes the tangent branch, which N610Q handles by doubling.

## Success criterion

For all selected levels, the generic weighted representative is regular, the
generic transition factor has valuation zero, and the identity holds in the
declared truncated formal ring.

## Falsification criterion

Any negative valuation or failed identity narrows the direction transition to
the already sampled directions and requires a different chart.

## Reproduction command

```bash
sage -python tools/poincare_true_pole_generic_direction_n610t.py \
  --out notes/poincare_true_pole_generic_direction_n610t.json
```
