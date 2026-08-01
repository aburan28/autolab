# Experiment Contract: N608V linear-projection factor-base screen

## Hypothesis

For at least one declared complete H018 pencil level, a projective linear map

```text
ell_(a,b)(P,Q) = a P + b Q
```

from its source pairs to the elliptic-curve group has image support materially
smaller than the base-field scale.  Such a map would be a necessary first
ingredient for a linear-projection factor base.

## Null Hypothesis

Every projective direction has image support of order `#E(F_103)`, so linear
projection alone does not produce a sub-base-field factor base.

## Parameters

- Registered `F_103` H018 fixture, with curve order `109`.
- Complete N608U source sets at levels `t=0,1,17`.
- All projective directions: `(1,b)` for `b` in `F_103`, plus `(0,1)`.
- Baseline scale: `#E(F_103)=109`.

## Metrics

- complete source-set size;
- distinct projected images for every direction;
- minimum image size and its direction;
- minimum image fraction of the base curve order.

## Positive Control

The source generator must reproduce N608U's declared complete source counts
before projections are measured.

## Negative Control

An image support at most one quarter of the base curve order for every
declared level would reject the null and retain this linear route for a
relation-law experiment.  Failure of that gate is a scoped negative for these
linear maps only.

## Boundary

This screen does not construct relations or compute factor-base logarithms.
It excludes neither nonlinear correspondences nor a Jacobian, Prym, cover, or
other non-product-linear representation.

## Reproduction

```bash
sage -python tools/poincare_linear_projection_n608v.py --out notes/poincare_linear_projection_n608v.json
sage -python tools/verify_poincare_linear_projection_n608v.py --primary notes/poincare_linear_projection_n608v.json --out notes/poincare_linear_projection_n608v_structural_replay.json
```
