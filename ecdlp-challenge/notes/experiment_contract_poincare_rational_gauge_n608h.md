# Experiment Contract: N608H low-pole Poincare gauge screen

## Hypothesis

The unique finite-field 109th roots of the N608F full-cycle transport products
are evaluations of one low-pole rational gauge `N/D` on the common regular
domain, rather than unrelated pointwise normalizations.

## Parameters

- `E/F_103`, common regular points in `F_(103^2)`
- 24 deterministic regular sample points
- rational numerator and denominator spaces `L(dO)`, `d=2,4,6`

## Success criterion

One degree has a nonzero homogeneous fit with nonzero numerator and
denominator blocks and with denominator nonzero on every sample.

## Falsification criterion

All bounded-pole fits have zero nullity or only denominator-invalid solutions.

## Reproduction command

```bash
sage -python tools/poincare_rational_gauge_probe_n608h.py \
  --out notes/poincare_rational_gauge_probe_n608h.json
```
