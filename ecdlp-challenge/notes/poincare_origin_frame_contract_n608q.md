# Experiment Contract: N608Q formal origin frame for the moving Cauchy basis

## Hypothesis

As the N606Y support point `r` approaches the origin with formal parameter
`u=-x(r)/y(r)`, the moving kernel

```text
k_r(P) = (y(P)+y(r))/(x(P)-x(r))
```

has a regularized order-eight limit that supplies the missing `x^3 y`
direction of `H^0(E,O(9O))`. This gives a declared source-side origin frame
without a fitted finite-orbit gauge.

## Exact Test

Expand `k_r(P)` through `u^8` using the formal group of
`E/F_103: y^2=x^3+x+24`. Express each coefficient in the standard basis

```text
1, x, y, x^2, xy, x^3, x^2y, x^4, x^3y
```

of `L(9O)`, with held-out point controls.

## Success Criterion

All coefficients through `u^7` lie in the first eight-dimensional subspace,
while the `u^8` coefficient has nonzero `x^3y` component. The resulting
regularized kernel limits to a new independent direction.

## Boundary

This constructs only the source-side formal origin frame. The target-side
evaluation normalization and origin quotient direction remain open.

## Reproduction

```bash
sage -python tools/poincare_origin_frame_n608q.py --out notes/poincare_origin_frame_n608q.json
sage -python tools/verify_poincare_origin_frame_n608q.py --primary notes/poincare_origin_frame_n608q.json --out notes/poincare_origin_frame_n608q_independent_verifier.json
```
