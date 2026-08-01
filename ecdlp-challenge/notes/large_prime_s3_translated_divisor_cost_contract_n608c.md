# Experiment Contract: N608C translated-divisor backend cost audit

## Hypothesis

For the N608B signed-pair-plus-secondary relation shape, a translated-secondary
divisor and GCD may be a usable relation-membership backend, rather than merely
an exact certificate.

## Null hypothesis

Once construction of the translated divisor, polynomial GCD, root recovery,
source verification, and both backends' precomputation are charged, the
translated-divisor backend is not faster than the existing per-secondary S3
quadratic scan on the same toy instances.

## Parameters

- field/curve family: public generated prime-field short-Weierstrass challenge
  curve at 20 bits, seed `0x12345678`.
- sizes: primary factor-base sizes `B in {8, 16}` and secondary bucket sizes
  `S in {16, 32}`.
- secondary buckets: `low_tail` and `random_x`.
- query suite: four constructed positive targets and 24 deterministic random
  scalar targets per configuration, replayed 16 times for timing stability.
- factor base: low-x points, signed in the primary pair and secondary term.
- relation shape: `T = A + L`, where `A` is a signed primary pair and `L` is a
  signed secondary point.
- baseline: `large_prime_s3_probe.partials_for_target`, which evaluates the
  S3 quadratic for every unsigned secondary point and verifies candidate pair
  hits.

## Metrics

- primary pair-table build time and signed-pair entry count;
- primary divisor construction time and degree;
- per-query translated point additions, divisor construction time, GCD/root /
  source-recovery time, candidate root count, and verified sources;
- baseline S3 scan time, quadratic root tests, square-root candidates, and
  verified partials;
- separately reported local wall-clock ratio. This is an implementation metric,
  not a field-operation equivalence claim.

## Positive control

Each constructed target is formed directly from a sampled signed primary pair
and signed secondary point. Both methods must return at least one verified
relation source for the raw target.

## Negative control

Deterministic random scalar targets are processed by both methods. Any divisor
source is checked by an explicit elliptic-curve addition; no root or x-match is
accepted by itself.

## Success criterion

The backend may be called a local membership signal only if all positive
controls verify, all recovered sources verify, and its fully measured query
wall-clock is no greater than the direct S3 scan on every configuration.
This does not establish a full index-calculus improvement: collision rank,
relation probability at relevant scales, sparse linear algebra, and target
descent remain separate obligations.

## Falsification criterion

Any incorrect source, missed constructed source, or fully charged local
slowdown retains the result as an exact-certificate or negative-backend result,
not an ECDLP speedup.

## Reproduction command

```bash
sage -python tools/large_prime_s3_translated_divisor_cost_n608c.py \
  --gen-instance target/release/gen_instance \
  --out notes/large_prime_s3_translated_divisor_cost_n608c.json
sage -python tools/verify_large_prime_s3_translated_divisor_cost_n608c.py \
  --result notes/large_prime_s3_translated_divisor_cost_n608c.json
```

## Claim boundary

`HYPOTHESIS / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM`. The tested object
is only the large-prime S3 relation-membership backend in the specified raw
x-coordinate representation.

## Results

All eight configurations passed constructed-source recovery and independent
source verification for both backends. Each configuration replayed the same
28-target suite 16 times. The fully charged divisor-to-direct local wall-clock
ratios ranged from `1.1196` to `1.2682`; every row was slower than the direct
S3 scan. The independent verifier also confirms the expected signed-pair
degree, secondary-work count, query shape, and random-hit agreement.

## Interpretation

This is a negative result for the translated-divisor backend in the tested raw
x-coordinate model. The exact GCD certificate remains reusable for correctness
or another representation, but it supplies no local speed signal here. A full
ECDLP assessment would still require relation probability, collision rank,
linear algebra, target descent, and a common field-operation cost model.
