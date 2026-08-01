# Experiment Contract: N609O coefficient-ratio H018 relation law

## Hypothesis

An explicit nonlinear scalar `a(Q)` formed by a moving-frame coefficient ratio
makes `P+a(Q)Q` determined by a certified complete H018 pencil level.

## Parameters

- H018 `E/F_103`, order 109.
- Certified levels `0,1,17`.
- All ratios `c_i/c_7` and `c_i/c_8` with nonzero denominators.

## Success criterion

One coefficient-ratio map has a singleton image on a certified level.

## Falsification criterion

All admitted maps retain more than one value on every certified level.

## Boundary

This tests a compact, explicit coefficient-ratio family only. It does not
rule out other nonlinear corrections, normalization/Jacobian relations, or an
end-to-end ECDLP improvement.

## Reproduction

```bash
sage -python tools/poincare_coefficient_ratio_relation_preflight_n609o.py --out notes/poincare_coefficient_ratio_relation_n609o.json
```
