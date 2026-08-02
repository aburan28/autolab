# Experiment Contract: N610I reverse H018 bivariate numerator

## Hypothesis

The `v-u`-cleared N610H numerator is regular in the reverse iterated Laurent
order, or else exhibits a reproducible reverse pole profile identifying the
missing transition correction.

## Parameters

- field/curve family: fixed H018 fixture on `E/F_103: y^2=x^3+x+24`
- levels: `0`, `1`, and `17`
- outer ring: Laurent series in `u` over Laurent series in `v`
- inspected terms: `u^0` through `u^15`

## Success criterion

All levels either have nonnegative reverse valuations or share one explicit
nonregular valuation profile.

## Reproduction command

```bash
sage -python tools/poincare_reverse_bivariate_numerator_n610i.py \
  --out notes/poincare_reverse_bivariate_numerator_n610i.json
```
