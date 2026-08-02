# Experiment Contract: N610G H018 formal slope-degree bound

## Hypothesis

For the Cauchy-cleared H018 normal form, every tested odd `u` coefficient
vanishes and the coefficient of `u^r` has slope degree at most `r+1`.

## Null hypothesis

An odd coefficient is nonzero, or an even coefficient through `u^36` fails
held-out replay at degree `r+1`.

## Parameters

- field/curve family: fixed H018 fixture on `E/F_103: y^2=x^3+x+24`
- levels: `0`, `1`, and `17`
- orders: `u^0` through `u^36`
- slopes: deterministic values `2` through `49`
- factor base: not applicable; formal normalization diagnostic
- baseline: N610D degree profile through `u^24`

## Metrics

- formal group evaluations: 144 total level/slope samples
- interpolation degree: `r+1` at each even coefficient
- held-out slope replay: at least ten points at the maximum order
- wall-clock: producer elapsed time

## Positive control

The constant coefficient replays at degree one and has root `83`.

## Negative control

Any violation of parity or the declared degree bound is a direct rejection.

## Success criterion

All levels pass the parity and degree-bound replay through `u^36`.

## Falsification criterion

Any failed held-out replay or nonzero odd coefficient.

## Reproduction command

```bash
sage -python tools/poincare_slope_degree_bound_n610g.py \
  --out notes/poincare_slope_degree_bound_n610g.json
```
