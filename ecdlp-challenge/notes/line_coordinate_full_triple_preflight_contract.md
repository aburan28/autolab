# Experiment Contract: REP-AUXLINE-018 Full Line-Coordinate Triple Preflight

## Hypothesis

For the direct line-coordinate factor family

```text
F = { P in E(F_p) : (y(P) - m*x(P) - r)^d = c },
```

the condition that a variable affine line has all three intersection points in
`F` can be expressed by a bounded-support quotient-remainder system.  If that
system has nonrandom source density and sub-pair-scale static state, it could
be a non-Cartesian source of zero-sum factor relations.

## Null hypothesis

The full-triple condition is either random-density geometry or its three
remainder coefficients have pair-scale-or-worse support.  In that case it does
not admit a source solver, even when individual toy lines give true triples.

## Parameters

- prime-order control curve: `E/F_101: y^2=x^3+x+32`, order `101`;
- base coordinate: `ell(P)=y-7*x-11`;
- split outer-coset degrees: `5,10,20,25,50`;
- one deterministic outer coset per degree;
- 64 equal-size random direct-point controls per degree.

## Metrics

- exact equality between quotient-remainder source lines and all direct-factor
  triples summing to zero;
- direct-factor triple count and matched-random distribution;
- total degree and monomial count of the three generic remainder coefficients;
- quotient support relative to unordered oriented pair state;
- exhaustive affine-line scan count relative to the rho operation scale.

## Positive control

Every source line must reconstruct three distinct direct-member points whose
sum is the identity, and every such direct-member triple must appear once as a
source line.

## Negative control

For each degree, compare the exact direct-factor triple count with 64 equal
size uniformly sampled point sets on the same prime-order curve.

## Success criterion

Advance only if a registered row has all of the following:

- an exact source/geometry biconditional;
- at least a fourfold triple-density lift over the matched random mean;
- total generic quotient-remainder support strictly below the unordered
  oriented-pair state; and
- a stated source inversion and target-descent plan that does not charge a
  full `p^2` affine-line scan.

## Falsification criterion

If every registered degree has pair-scale-or-worse remainder support, or no
row has a fourfold random-density lift, preserve the scoped negative and do
not implement a Groebner/resultant collector.

## Reproduction command

```bash
sage -python tools/line_coordinate_full_triple_preflight.py \
  --out notes/line_coordinate_full_triple_preflight_rep_auxline018.json
sage -python tools/verify_line_coordinate_full_triple_preflight.py \
  --primary notes/line_coordinate_full_triple_preflight_rep_auxline018.json \
  --out notes/line_coordinate_full_triple_preflight_independent_verifier_rep_auxline018.json
```
