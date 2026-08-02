# Result: N611O Cauchy-Coordinate Labels

## Status

`NEGATIVE BOUNDED_CAUCHY_COORDINATE_FACTOR_LABEL_SCREEN / MODEL-BOUND /
TOY-EVIDENCE / GLOBAL_NORMALIZATION_FUNCTIONS_OR_OTHER_RATIONAL_MAPS_OPEN /
NO_ECDLP_CLAIM`

## Exact Result

N611O evaluates all complete N608U sources using the four fixed labels

```text
d_x,  m,  m^2,  c_8(Q)m,
```

where `m=(y(P)+y(-4Q))/(x(P)-x(-4Q))`.  Every label is defined on every
declared open source and varies within multiple Q-fibres.  No label meets the
joint support threshold `<=27` and `<= 1/2` of the deterministic uniform
control at all three levels.

| Label | Level 0 | Level 1 | Level 17 |
|---|---:|---:|---:|
| `d_x` support | 30 | 44 | 56 |
| `m` support | 55 | 72 | 80 |
| `m^2` support | 28 | 36 | 40 |
| `c_8(Q)m` support | 31 | 46 | 53 |
| Uniform control support | 54 | 71 | 76 |

The natural sign-quotient candidate `m^2` is nearest, with support ratios
`0.5185..., 0.5070..., 0.5263...`, but it is above both the `27` support gate
and the `0.5` ratio gate at every level.  The product-linear Kummer control
`x(P-4Q)` has supports `29,37,44`.

## Evidence

- N608U source counts `84,108,144` and degree-ten, no-exceptional inverse
  conditions replay in the producer.
- The independent verifier separately recomputes all complete sources, all
  Cauchy labels, and the deterministic uniform controls; all five checks pass.
- No source is silently discarded for a Cauchy denominator: the declared open
  condition makes each denominator nonzero.

## Interpretation

This rejects only the four stated elementary Cauchy-coordinate labels.  The
partial compression of `m^2` is consistent with a basic sign quotient, but it
does not yield the required factor-base scale, and there is no target-bearing
relation, rank, individual-log descent, charged cost comparison, or ECDLP
speedup.  Changing exponents or thresholds after this panel would be selection
on the data and is not promoted here.

## Next Concrete Action

Return to the global normalization route.  Any further rational-map candidate
must be derived from an explicit global divisor or normalization function and
must enter with a predeclared target-return and relation contract.

## Reproduction

```bash
sage -python tools/poincare_cauchy_coordinate_labels_n611o.py \
  --out notes/poincare_cauchy_coordinate_labels_n611o.json
sage -python tools/verify_poincare_cauchy_coordinate_labels_n611o.py \
  --primary notes/poincare_cauchy_coordinate_labels_n611o.json \
  --out notes/poincare_cauchy_coordinate_labels_n611o_independent_verifier.json
```
