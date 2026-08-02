# Experiment Contract: N610K common H018 bivariate coefficient grid

## Hypothesis

The two exact-Cauchy Laurent expansions represent the same local numerator on
a common finite coefficient grid.

## Parameters

- field/curve family: fixed H018 fixture on `E/F_103: y^2=x^3+x+24`
- levels: `0`, `1`, and `17`
- grid: `v^i u^j` for `0 <= i,j < 12`
- inputs: exact Cauchy parameter from N610J

## Success criterion

All 144 coefficient comparisons agree at every level, with nonzero entries
present on every grid.

## Reproduction command

```bash
sage -python tools/poincare_common_bivariate_grid_n610k.py \
  --out notes/poincare_common_bivariate_grid_n610k.json
```
