# Experiment Contract: N607D Moving Cauchy-Kernel Section Ansatz

## Hypothesis

For `R=zQ`, the rational function

```text
k_R(P) = (y(P)+y(R)) / (x(P)-x(R))
```

has only the expected poles at `P=R` and `P=O`.  Together with the eight
standard functions of pole order at most eight at `O`, it gives an explicit
fibre basis of `H^0(E,O(8O+[R]))`.  This may make the N606U/N607A Poincare
frame evaluable enough to search for two globally compatible sections.

## Null hypothesis

The kernel has an additional pole, degenerates on the required parameter
domain, or the resulting moving basis cannot satisfy the Poincare transition
law.  Then this ansatz cannot construct the H018 section space.

## Parameters

- `E/F_103: y^2=x^3+x+24`.
- `z=-3-pi`, acting as `[-4]` on `E(F_103)`.
- Fibre divisor `8[O]+[zQ]`.
- Base functions `1,x,y,x^2,xy,x^3,x^2y,x^4` plus `k_(zQ)`.

## Metrics

- fibre-basis rank at multiple rational `Q` values;
- pole and exceptional-point behavior of `k_(zQ)`;
- compatibility with N606Y transport at common regular points;
- coefficient-function degrees needed for global two-section candidates.

## Positive control

For a nonzero `R`, `k_R` must have a simple pole at `R`, be regular at `-R`,
and be independent from `L(8O)`.

## Negative control

The origin parameter and a sign-flipped numerator must not be accepted as the
same regular kernel without separate limiting data.

## Boundary

Passing only constructs a fibre basis. It does not produce globally compatible
H018 sections, a smooth pencil member, a factor base, relations, rank, target
descent, sub-rho cost, or an ECDLP improvement.

## Reproduction

```bash
sage -python tools/product_kummer_h018_cauchy_kernel_section_n607d.py \
  --output notes/product_kummer_h018_cauchy_kernel_section_n607d.json
```
