# Experiment Contract: N610D coefficient-wise H018 slope degree

## Hypothesis

After clearing the Cauchy unit, each fixed formal `u` coefficient of the H018
normal form admits a low-degree polynomial in the exceptional slope, but the
required degree grows with the coefficient index.

## Null hypothesis

The degree-nine model replays every tested coefficient, or no bounded
coefficient-wise polynomial model replays held-out slope samples.

## Parameters

- field/curve family: the fixed H018 fixture on `E/F_103: y^2=x^3+x+24`
- sizes: coefficients `u^0` through `u^24`; candidate degrees `0` through `26`
- seeds: deterministic slopes `2` through `31`
- factor base: not applicable; local formal-chart diagnostic
- relation shape: Cauchy-cleared normal-form coefficient versus slope
- baseline: degree-nine whole-series interpolation rejected by N610C

## Metrics

- group operations: formal group expansions per level and slope
- field operations: Vandermonde interpolation/replay per coefficient and degree
- memory: in-process exact Laurent-series samples
- relation probability/rank/solver degree: not applicable
- wall-clock: producer elapsed time

## Positive control

The constant coefficient must replay at degree one, matching N609W/N610B.

## Negative control

At least one higher coefficient must reject degree nine, matching N610C's
whole-series rejection.

## Success criterion

Every tested coefficient has a replayed degree at most 26, with a recorded
minimal degree profile, and both controls pass.  Exact degree agreement between
pencil levels is not required: level-dependent cancellations are admissible.

## Falsification criterion

Any coefficient has no replayed degree at most 26, the constant-term control
fails, or all tested coefficients replay at degree nine.

## Reproduction command

```bash
sage -python tools/poincare_slope_degree_profile_n610d.py \
  --out notes/poincare_slope_degree_profile_n610d.json
```
