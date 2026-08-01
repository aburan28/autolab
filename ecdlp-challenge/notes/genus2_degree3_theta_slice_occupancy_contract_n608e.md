# Experiment Contract: N608E degree-three component theta-slice occupancy

## Hypothesis

The N608D degree-three elliptic component may meet triple sums of a public
genus-two theta slice more often than matched random Jacobian subsets, while
retaining exact source labels and the same-form target decomposition.

## Null hypothesis

A theta slice formed from rational curve-point divisor classes has ordinary
occupancy in the embedded 103-point elliptic subgroup. Any apparent relation
count is explained by its finite labelled size or by the explicit isogeny
component, not a new decomposition structure.

## Parameters

- N608D public fixture over `F_101` with `#E(F_101)=103` and
  `#J_C(F_101)=11536=103*112`.
- factor bases: three deterministic theta slices of `B=4,6,8` rational curve
  points, each expanded by its hyperelliptic inverse.
- relation shape: a labelled triple sum in `J_C(F_101)` equals a point in the
  induced 103-point elliptic image.
- controls: 24 deterministic negation-stable pseudo-random Jacobian subsets
  of equal cardinality per slice, generated as random combinations of public
  rational-divisor classes without enumerating the full Jacobian.
- targets: all nonzero points of the induced elliptic image, identified only
  through the explicitly constructed toy pullback map.

## Metrics

- exact labelled triple-source count and target coverage;
- subgroup-hit count among all triple sums;
- coefficient-row rank modulo 103;
- matched-control mean/max occupancy and coverage;
- explicit source replay for every accepted candidate;
- target descent coverage using the same triple-sum enumerator.

## Positive control

The induced elliptic subgroup itself must replay exactly via N608D. A factor
base enlarged to include three known summands of an embedded target must yield
at least one explicitly replayed target decomposition.

## Negative control

Random negation-stable Jacobian subsets of the same labelled size use the same
enumerator, target set, and rank calculation. A toy scalar lookup may label
the known 103-point subgroup only; it must not be treated as an ECDLP oracle.

## Success criterion

A local representation signal requires all source replays, full coefficient
rank for at least one nontrivial slice, target coverage above every matched
control, and a multi-scale lift that survives the exact label cost. It still
does not establish a full ECDLP improvement without charged relation
collection, sparse linear algebra, and blind target descent.

## Falsification criterion

If the slice has no nonrandom occupancy/coverage lift, lacks rank, or requires
enumerating the labelled triple surface for descent, retain the narrow negative
result and move to a different factor-base geometry.

## Reproduction command

```bash
sage -python tools/genus2_degree3_theta_slice_occupancy_n608e.py \
  --out notes/genus2_degree3_theta_slice_occupancy_n608e.json
```

## Claim boundary

`HYPOTHESIS / INTERNAL-JACOBIAN-OCCUPANCY / MODEL-BOUND / TOY-EVIDENCE /
NO_ECDLP_CLAIM`.

## Results

All source replays passed across `B=4,6,8` and 24 matched controls each. The
theta-slice hit ratios against control means were `0.00`, `1.85`, and `2.11`,
but target coverages `0/4/10` never exceeded control maxima `4/6/12` and the
coefficient ranks were only `0/4`, `2/6`, and `4/8`.

## Interpretation

The apparent mean-hit lift at the larger slices is insufficient: it does not
beat the control tail and it does not produce full factor-log rank. This rejects
the tested rational curve-point theta-slice geometry. It does not rule out a
different divisor slice, a non-enumerative packet inverse, or another
internal-Jacobian construction.
