# Experiment Contract: N608B Translated Secondary-Divisor GCD

## Hypothesis

For a large-prime S3 relation `W=A+L`, where `A` is a signed primary pair,
the gcd of the primary-pair x-divisor and the target-translated secondary
x-divisor can return an exact source witness with lower effective cost than
scanning all secondary points.

## Null hypothesis

The gcd is an exact algebraic restatement of the same secondary scan: it needs
a degree-`Theta(S)` translated divisor per target and its decoded roots provide
no relation-density or source-inverse improvement.

## Parameters

- Generated 20-bit prime-order curve, seed `0x12345678`.
- Low-x primary bases `B=8,16`.
- Low-tail and random-x secondary buckets `S=16,32`.
- Signed pair and secondary representations from the existing S3 harness.

## Metrics

- Pair-divisor degree and translated-secondary-divisor degree.
- Exact gcd/source-witness agreement on constructed positives.
- Exact agreement between gcd roots and direct secondary scan on random targets.
- Number of translated roots per target and required source postings.

## Positive control

A target constructed as a signed pair plus signed secondary point must yield a
nonconstant gcd and an exact signed source replay.

## Negative control

Random-target gcd roots must agree exactly with direct scan witnesses; any
additional root is a failure.

## Success criterion

Promotion requires a demonstrated query or inverse cost below the direct scan,
not merely an exact gcd identity.

## Reproduction

```bash
sage -python tools/large_prime_s3_translated_divisor_n608b.py --out notes/large_prime_s3_translated_divisor_n608b.json
```

## Results

All eight configurations (`B=8,16`, `S=16,32`, and low-tail/random-x
secondary buckets) passed the constructed-source and random-source agreement
checks. The primary divisor degree was exactly `B^2` after the x-sign quotient,
and the target-translated secondary divisor degree was exactly `2S`.

## Interpretation

The GCD is an exact source certificate for this raw representation, but it does
not remove the secondary linear degree or the primary quadratic degree. N608C
therefore charges it as a backend before any relation-collection claim.
