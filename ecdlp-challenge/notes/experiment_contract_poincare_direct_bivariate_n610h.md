# Experiment Contract: N610H direct ordered bivariate H018 numerator

## Hypothesis

Keeping the Q-dependent determinant in its univariate formal ring and clearing
the leading Cauchy factor `v-u` leaves a nonnegative ordered Laurent expansion
in the independent P parameter `v` and Q parameter `u`.

## Null hypothesis

The direct numerator has a negative v valuation or an inspected v coefficient
has a negative u valuation.

## Parameters

- field/curve family: fixed H018 fixture on `E/F_103: y^2=x^3+x+24`
- levels: `0`, `1`, and `17`
- Q ring: Laurent series in `u`, precision 120
- P ring: Laurent series in `v` over the Q ring, precision 48
- relation shape: local determinantal numerator only

## Metrics

- v valuation of each direct numerator
- u valuation of v coefficients zero through fifteen
- direct univariate matrix determinant/solve and bivariate basis evaluation

## Success criterion

All three levels have nonnegative v valuation and nonnegative inspected u
valuations.

## Falsification criterion

Any negative valuation in the declared ordered chart.

## Reproduction command

```bash
sage -python tools/poincare_direct_bivariate_numerator_n610h.py \
  --out notes/poincare_direct_bivariate_numerator_n610h.json
```
