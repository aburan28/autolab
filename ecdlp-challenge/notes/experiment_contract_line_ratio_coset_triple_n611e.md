# Experiment Contract: N611E Line-Ratio Coset Triple Preflight

## Hypothesis

A factor base defined by a nontrivial ratio of affine line functions

```text
R(P) = L_1(P) / L_2(P),   R(P)^d = c
```

may have a three-point zero-sum membership system whose denominator-cleared
remainder is smaller than the direct pair state and whose relation density is
nonrandom. This differs from the already-rejected single-line and Dickson
families.

## Null hypothesis

For the registered line-ratio maps, the generic all-three-points condition
has pair-scale-or-larger support or the resulting zero-sum triples occur at
matched-random density. In that case no source solver, rank computation, or
descent is admitted.

## Parameters

- Curve: `E/F_101: y^2=x^3+x+32`, prime group order `101`.
- Maps, fixed before execution:
  - `y / (y-x)`;
  - `(y-1) / (y+x+1)`;
  - `(y-2x-1) / (y+x-2)`;
  - `(y-x-1) / (y+2x+1)`.
- Exponents: `d in {5,10,20,25}`.
- Cosets: `c in {1,g^d}`, where `g` is the least primitive generator of
  `F_101^*`; both values are in the image of the `d`th-power map.
- Matched controls: 32 deterministic equal-size samples of defined affine
  nonidentity points per admitted row.

## Exact source relation

For a source line `y=s*x+t`, let

```text
C(x) = x^3-s^2*x^2+(1-2*s*t)*x+(32-t^2)
N(x) = L_1(x,s*x+t)
D(x) = L_2(x,s*x+t).
```

The line supplies a complete affine three-point source exactly when

```text
N(x)^d - c*D(x)^d = 0 mod C(x).
```

The producer must replay this condition for every enumerated distinct
factor-base zero-sum triple and verify that its line is nonvertical.

## Metrics

- factor-base size and denominator exclusions;
- exact triple/source replay;
- generic remainder monomial count and total degree;
- remainder support divided by the unordered pair state;
- distinct zero-sum triples and matched-random density lift.

## Positive control

All direct factor-base triples summing to the identity must reproduce the
quotient-remainder condition with their public affine line witness.

## Negative controls

- A factor-base point with `L_2(P)=0` is rejected, not silently assigned a
  ratio value.
- Every control has the same point count and draws only from points where the
  candidate ratio is defined.
- A candidate must exceed every control and have at least a fourfold density
  lift; a mean-only fluctuation does not promote it.
- A passing preflight still establishes no source solver, factor-log rank,
  target descent, or ECDLP speedup.

## Admission criterion

Advance one fixed row only if its factor base is between 6 and 48 points, all
triple witnesses replay, generic remainder support is below its unordered pair
state, and its triple density is at least four times the control mean and
strictly above the control maximum.

## Reproduction command

```bash
sage -python tools/line_ratio_coset_triple_probe_n611e.py \
  --out notes/line_ratio_coset_triple_probe_n611e.json
```

## Claim boundary

`HYPOTHESIS / LINE_RATIO_COSET_TRIPLE_PREFLIGHT / MODEL-BOUND /
TOY-EVIDENCE / NO_ECDLP_CLAIM`.
