# Experiment Contract: N610M exact H018 origin smoothness screen

## Hypothesis

The exact-Cauchy local numerator of each selected H018 member vanishes at the
simultaneous origin with a nonzero linear gradient.

## Parameters

- field/curve family: fixed H018 fixture on `E/F_103: y^2=x^3+x+24`
- levels: `0`, `1`, and `17`
- local chart: N610J exact Cauchy factor `v-w(u)`

## Success criterion

Every member has zero constant coefficient and at least one nonzero linear
coefficient in `v,u`.

## Reproduction command

```bash
sage -python tools/poincare_exact_origin_smoothness_n610m.py \
  --out notes/poincare_exact_origin_smoothness_n610m.json
```
