# Experiment Contract: N610J exact H018 Cauchy transition

## Hypothesis

Replacing the tangent factor `v-u` with the literal pole factor
`v-w(u)`, where `w(u)=-x(-4Q)/y(-4Q)`, regularizes the local determinant
numerator in both Laurent orders.

## Parameters

- field/curve family: fixed H018 fixture on `E/F_103: y^2=x^3+x+24`
- levels: `0`, `1`, and `17`
- inspected coefficients: zero through fifteen in each Laurent order

## Success criterion

The support parameter has valuation one with nonzero linear term, and both
orders have nonnegative numerator and inspected coefficient valuations.

## Reproduction command

```bash
sage -python tools/poincare_exact_cauchy_transition_n610j.py \
  --out notes/poincare_exact_cauchy_transition_n610j.json
```
