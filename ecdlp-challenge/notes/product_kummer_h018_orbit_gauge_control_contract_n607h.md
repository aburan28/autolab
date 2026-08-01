# Experiment Contract: N607H Finite-Orbit Gauge Control

## Hypothesis

The N607E finite-orbit transition matrices provide an intrinsic monodromy that
can be compared directly to the N606U pushforward bundle.

## Null hypothesis

Without a geometric rigidification, the nonclosing edge matrices are gauge
choices; a finite-orbit matrix product is not itself a bundle identification.

## Parameters

- H018 toy curve `E/F_103: y^2=x^3+x+24`.
- The 109-point N606Y `[84]` base orbit.
- N607E Cauchy basis and explicit Miller transport.

## Metrics

- Invertibility of all 109 reconstructed edge matrices.
- Number of nonclosing edges gauged exactly to identity.
- Rank and scalar status of the closing product.

## Positive control

Recursive gauges `G_(i+1)=C_i G_i` must trivialize each of the first 108
edge matrices exactly.

## Negative control

The closing product must not be promoted to a normalized-Poincare transition,
even when it is invertible or nonscalar.

## Success criterion

The result records the raw finite-orbit data as a control only and requires a
geometric Poincare rigidification or genuine-cover bundle morphism next.

## Reproduction

```bash
sage -python tools/product_kummer_h018_orbit_gauge_control_n607h.py --output notes/product_kummer_h018_orbit_gauge_control_n607h.json
python3 tools/verify_product_kummer_h018_orbit_gauge_control_n607h.py --primary notes/product_kummer_h018_orbit_gauge_control_n607h.json --output notes/product_kummer_h018_orbit_gauge_control_n607h_independent_verifier.json
```
