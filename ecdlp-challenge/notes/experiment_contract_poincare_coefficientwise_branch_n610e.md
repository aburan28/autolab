# Experiment Contract: N610E coefficient-wise H018 formal branch

## Hypothesis

The coefficient-wise replayed slope models through `u^24` determine a valid
bounded formal branch of the level-zero H018 member through slope `83`, and
the two independent pencil differences remain order-eight units on that branch.

## Null hypothesis

The individually replayed coefficient models do not yield a simple bounded
branch, or a pencil difference changes order after the branch correction.

## Parameters

- field/curve family: fixed H018 fixture on `E/F_103: y^2=x^3+x+24`
- sizes: `u`-jet through order `24`; per-coefficient degree cap `26`
- seeds: deterministic slopes `2` through `31`
- factor base: not applicable; local formal-chart diagnostic
- relation shape: Cauchy-cleared normal-form pencil differences
- baseline: N610A uncorrected order-eight contact and N610D replay profile

## Metrics

- group operations: formal group expansions per level and slope
- field operations: per-coefficient interpolation and Newton updates
- memory: exact series samples and coefficient polynomials
- relation probability/rank/solver degree: not applicable
- wall-clock: producer elapsed time

## Positive control

The `u^0` coefficient has simple root `83` and a nonzero slope derivative.

## Negative control

The producer rejects any coefficient that cannot replay at degree at most `26`.

## Success criterion

All coefficients through `u^24` replay on held-out slopes, the level-zero
residual vanishes through order `24`, and both corrected differences have
order-eight unit leading terms.

## Falsification criterion

Any replay failure, nonsimple root, residual below order `24`, or a corrected
difference whose first nonzero term is not an order-eight unit.

## Reproduction command

```bash
sage -python tools/poincare_coefficientwise_branch_n610e.py \
  --out notes/poincare_coefficientwise_branch_n610e.json
```
