# Experiment Contract: N608D degree-three genus-two elliptic component

## Hypothesis

An explicit degree-three elliptic subcover of a genus-two curve can provide a
non-product Jacobian representation worth testing for ECDLP relation density.

## Null hypothesis

The subcover only yields an explicit isogeny-factor embedding. Until a
target-independent divisor factor base, a nonrandom relation density, a source
inverse, rank, and same-form target descent are demonstrated, it is not an
ECDLP improvement.

## Literature basis

The fixture uses the nondegenerate degree-three elliptic-subcover normal form
in Shaska, *Genus 2 fields with degree 3 elliptic subfields*, arXiv:math/0109155,
equations (2), (5), and (7). That work describes a degree-three elliptic
subfield; this experiment does not treat the cited construction as an attack.

## Parameters

- base field: `F_101`;
- normal-form parameters: `a=0`, `b=4`;
- original curve: `Y^2=(X^3+4X+1)(4X^3+16X^2+8X+1)`;
- rational Weierstrass root: `X=86`, used to move to an odd-degree model;
- elliptic quotient: the resulting `E/F_101` has prime order `103`;
- map chart: `U=X^2/(X^3+4X+1)` and the normalized version of the published
  `V` coordinate.

## Metrics

- smoothness of the genus-two model and quotient;
- exact map replay for every rational source point in the regular chart;
- full degree-three inverse polynomial and exceptional-fibre accounting;
- reduction of every quotient fibre to a Mumford pair;
- exact homomorphism replay for every pair in `E(F_101)^2`;
- cardinality of the induced Jacobian pullback image.

## Positive control

The pullback of every quotient sum must equal the sum of the two pullbacks.
The `103` quotient points must produce `103` distinct Jacobian image points.

## Negative controls

- Treating an exceptional map-chart point as a regular finite source is an
  error.
- A raw degree-three inverse polynomial without reduction and exact Jacobian
  replay is not accepted as a target lift.
- This map must not be called a relation packet, factor base, or ECDLP speedup.

## Success criterion

This contract passes only as a materialization receipt when all map, divisor,
and homomorphism checks pass. Admission to an occupancy experiment further
requires an explicit theta/divisor slice, matched random subsets, source
recovery, relation rank, and target descent.

## Reproduction command

```bash
sage -python tools/genus2_degree3_elliptic_component_n608d.py \
  --out notes/genus2_degree3_elliptic_component_n608d.json
```

## Claim boundary

`HYPOTHESIS / EXPLICIT-COVER-MATERIALIZATION / MODEL-BOUND / TOY-EVIDENCE /
NO_ECDLP_CLAIM`.

## Results

The fixture is smooth and has a prime-order `103` elliptic quotient inside a
`11536=103*112` point Jacobian. Every quotient fibre satisfies the raw Mumford
congruence and reduces to degree at most two. The normalized pullback has 103
distinct images and passes all `103^2=10609` homomorphism pairs exactly.

## Interpretation

This furnishes a reproducible non-product coordinate model with an explicit
elliptic component, but the component is still an explicit toy isogeny factor.
N608E is required to test a concrete factor-base geometry rather than infer a
cryptanalytic benefit from the existence of the cover.
