# Experiment Contract: N608X raw Cauchy-coefficient interpolation screen

## Hypothesis

The nine base-dependent coefficients of the N608S raw moving Cauchy frame may
be functions of low pole order on the base elliptic curve. If so, substituting
them into the denominator-cleared level equation could give a low-bidegree
finite model for the complete pencil curve `C_0`.

## Null Hypothesis

The raw-frame coefficients require pole order comparable to the full number of
rational base points. In that case, interpolation of this frame cannot be used
as a low-complexity global curve equation.

## Parameters

- Registered N608S/N608U `F_103` H018 fixture.
- All 108 nonzero rational base points.
- Raw Cauchy coefficient vector of length nine.
- Riemann--Roch spaces `L(d O)` for `1 <= d <= 109`.
- Low-pole admission threshold `d <= 24`.

## Metrics

- the minimum `d` for which each coefficient is the restriction of an element
  of `L(d O)`;
- matrix rank and interpolation residual at that degree;
- low-pole admission for every coefficient.

## Controls

- Positive: the known coordinate functions `1`, `x`, and `y` must first fit in
  `L(O)`, `L(2O)`, and `L(3O)` respectively.
- Negative: each raw coefficient must fail interpolation in every declared
  low-pole space through `L(24O)` before a high-degree fit is accepted.

## Boundary

This tests only the unnormalized raw moving Cauchy basis on rational base
points. It does not rule out a normalized Poincare frame, a non-scalar vector
bundle description, or a different finite model of the pencil curve.

## Reproduction

```bash
sage -python tools/poincare_raw_coefficient_interpolation_n608x.py --out notes/poincare_raw_coefficient_interpolation_n608x.json
sage -python tools/verify_poincare_raw_coefficient_interpolation_n608x.py --primary notes/poincare_raw_coefficient_interpolation_n608x.json --out notes/poincare_raw_coefficient_interpolation_n608x_structural_replay.json
```
