# Result: N611M Level-Curve Canonical Elliptic Returns

## Status

`NEGATIVE CANONICAL_LEVEL_CURVE_ELLIPTIC_RETURN_COMPRESSION /
MODEL-BOUND / TOY-EVIDENCE / NONCANONICAL_NORMALIZATION_PACKET_OPEN /
NO_ECDLP_CLAIM`

## Exact Result

For each N608U complete open level `C_t^open`, N611M exhausts all 104
projective fixed returns

```text
(P,Q) -> [a]P + [b]Q,  (a:b) in P^1(F_103).
```

The direct canonical projection supports `(p1,p2)` and best span support are:

| Level | Sources | `p1` | `p2` | Best direction | Best support |
|---:|---:|---:|---:|---|---:|
| 0 | 84 | 62 | 60 | `(1:13)` | 49 |
| 1 | 108 | 74 | 68 | `(1:7)` | 56 |
| 17 | 144 | 76 | 76 | `(1:91)` | 68 |

The declared quarter-order threshold is `floor(109/4)=27`.  No direction at
any level reaches it.  For each level, the degree-ten inverse has no
exceptional branch, and direct `p1,p2` supports exactly equal the `(1:0)` and
`(0:1)` panel entries.

## Evidence

- Producer: four structural gates pass; the quarter-order return hypothesis
  is false.
- Independent verifier: all six checks pass after separately rebuilding the
  degree-nine cover interpolation and every source level from N608U primitives.
- The panel is exhaustive for the stated fixed canonical elliptic span, not a
  sample of directions.

## Interpretation

This does not reject the normalization/Jacobian route.  It rejects a narrow
way for that route to gain compression: a target recovery that is a fixed
linear combination of its two visible elliptic divisor projections.  A viable
normalization packet must instead construct an explicit noncanonical return
and then demonstrate source recovery, a target-bearing relation, rank,
individual-log descent, and fully charged sub-rho cost.

## Next Concrete Action

Derive a canonical projective or normalization model for `C_0`, including its
singular locus and explicit divisor classes.  Before relation collection,
either materialize a noncanonical target-return map or prove that the available
divisor classes factor through the now-screened canonical span.

## Reproduction

```bash
sage -python tools/poincare_level_curve_canonical_returns_n611m.py \
  --out notes/poincare_level_curve_canonical_returns_n611m.json
sage -python tools/verify_poincare_level_curve_canonical_returns_n611m.py \
  --primary notes/poincare_level_curve_canonical_returns_n611m.json \
  --out notes/poincare_level_curve_canonical_returns_n611m_independent_verifier.json
```
