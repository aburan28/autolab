# Experiment Contract: N609N additive H018 relation preflight

## Hypothesis

For some `c`, a certified complete H018 level determines `P+cQ`, giving a
direct target-bearing relation `P+cQ=H(t)`.

## Null hypothesis

Every additive map retains multiple values on every certified level.

## Parameters

- H018 `E/F_103`, order 109.
- Certified complete levels `0,1,17`.
- All 109 multipliers `c`.

## Success criterion

One multiplier has a singleton image on a certified level.

## Falsification criterion

No multiplier has a singleton level image.

## Boundary

This rejects direct additive level laws only. Nonlinear relations, compact
correction states, normalization/Jacobian mechanisms, and all end-to-end
ECDLP requirements remain open.

## Reproduction

```bash
sage -python tools/poincare_additive_relation_preflight_n609n.py --out notes/poincare_additive_relation_n609n.json
```
