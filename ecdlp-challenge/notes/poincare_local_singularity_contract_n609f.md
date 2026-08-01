# Experiment Contract: N609F open-chart local singularity screen

## Hypothesis

The declared H018 residual level curves have no rational singularity candidate
on their regular open chart: at every rational source, at least one exact
formal derivative in the P or Q direction is nonzero.

## Null hypothesis

Some declared open source has both local derivatives zero, producing a genuine
rational singularity candidate that must be resolved before any curve or
Jacobian claim.

## Parameters

- H018 fixture over `F_103` and all 108 nonzero base points.
- Levels `t=0,1,17`.
- P-direction formal translation on the target curve.
- Q-direction formal translation through a chosen cover lift, its deck orbit,
  and the exact graph-evaluation matrix.

## Metrics

- formal Q evaluator reconstruction on all regular rational `(P,Q)` pairs;
- P and Q derivative-zero counts at all open residual sources;
- simultaneous-zero count at each declared level.

## Controls

- The Q-direction formal evaluator must reproduce the raw Cauchy evaluator on
  all `108*106` regular rational pairs.
- Both local expansions must reproduce the recorded level at every screened
  source.
- The source counts remain N609A's open counts `84,108,144`.

## Success criterion

All reconstruction controls pass and no screened source has both derivatives
zero.

## Falsification criterion

Any simultaneous-zero source is a concrete singularity candidate. It blocks
promotion to a smooth compact curve until its local ring is analyzed.

## Boundary

This excludes only rational candidates on the specified open chart. It does
not cover the true Cauchy pole, the separately regularized `P=4Q` graph,
the `Q=O` boundary, non-rational geometric points, or a global compact model.
It cannot establish smoothness, genus, a Jacobian, factor bases, relations,
rank, descent, cost, or an ECDLP improvement.

## Reproduction

```bash
sage -python tools/poincare_local_singularity_probe_n609f.py --out notes/poincare_local_singularity_n609f.json
```
