# Experiment Contract: N611A non-split Mumford-coordinate slice

## Hypothesis

A target-independent factor base defined directly by a low-degree reduced
Mumford coordinate condition can give signed triple relations into an embedded
elliptic target subgroup more often, and with better row rank, than matched
random genus-two divisor bases.

## Null hypothesis

The coordinate slice has ordinary finite-group occupancy.  Any hit count is
explained by its factor-base size, and it does not provide a usable
non-split ECDLP relation mechanism.

## Parameters

- N608D's smooth genus-two curve over `F_101`, with Jacobian order `11536`;
- explicit degree-three elliptic quotient of prime order `103` and its exact
  Jacobian pullback embedding;
- factor bases: the canonical `D/-D` representatives with degree-two
  Mumford `u(t)=t^2+a*t+b` and `a` in `{0}` or `{0,1}`;
- relations: signed triples whose Jacobian sum equals an embedded nonzero
  target-subgroup point;
- controls: eight deterministic uniform samples of canonical degree-two
  Mumford representatives of the same size for each coordinate slice.

## Metrics

- factor-base size and intersection with the embedded target subgroup;
- exact source replay for every accepted relation;
- target-subgroup hit count and coverage;
- coefficient-row rank modulo `103`;
- matched-control mean and maximum hit/coverage/rank;
- all charged triple candidates and wall-clock time.

## Positive control

Three declared embedded subgroup points, together with their inverses, must
produce a replayed triple relation.  This checks the enumerator, not the
candidate factor base.

## Negative controls

- Random degree-two Mumford bases use the same signed-triple enumerator and
  embedded-target membership test.
- Scalar target labels are used only to count distinct targets after equality
  testing; they must not select factors or accepted rows.
- A hit or a rank observation is not a return map, factor-log system, blind
  descent, cost comparison, or ECDLP speedup.

## Success criterion

Promotion requires a coordinate slice to exceed every matched control in
target coverage and row rank, with complete source replay and a specified
same-form path back to the original ECDLP target.  It remains only a
preflight until relation collection, sparse linear algebra, blind descent,
and charged rho comparison succeed.

## Falsification criterion

If no slice exceeds the matched-control maxima, or its rank is deficient,
preserve the coordinate-slice negative result and search a different
Mumford-factor geometry or non-enumerative packet inverse.

## Reproduction command

```bash
sage -python tools/mumford_coordinate_slice_probe_n611a.py \
  --out notes/mumford_coordinate_slice_probe_n611a.json
```

## Claim boundary

`HYPOTHESIS / NONSPLIT_MUMFORD_FACTOR_BASE_PREFLIGHT / MODEL-BOUND /
TOY-EVIDENCE / NO_ECDLP_CLAIM`.
