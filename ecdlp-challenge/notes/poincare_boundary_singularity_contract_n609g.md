# Experiment Contract: N609G rational boundary singularity screen

## Hypothesis

The locally extended `P=4Q` graph and the rational finite `Q=O` boundary
zeros contain no simultaneous formal-derivative zero for the declared H018
pencil levels.

## Null hypothesis

One of these rational boundary points has both local derivatives zero, giving
a concrete singularity candidate that prevents a smooth-curve interpretation.

## Parameters

- H018 fixture over `F_103`.
- Levels `t=0,1,17`.
- Restored graph sources from N609E.
- Q-origin sections and finite zero sets from N609C/N609D.

## Metrics

- P/Q local derivative-zero counts at restored graph sources;
- P/Q local derivative-zero counts at finite rational Q-origin zeros;
- L(9O) P-origin zero multiplicity for every Q-origin pencil member.

## Controls

- Graph source counts and regularized values must remain `8,0,2`.
- Finite Q-origin source counts must remain `2,0,0`.
- Every local expansion must reproduce its declared level or boundary zero.

## Success criterion

No screened boundary point has both derivatives zero, and every Q-origin
pencil member has the prior simple P-origin zero.

## Falsification criterion

Any simultaneous-zero point is a singularity candidate requiring a local-ring
analysis before a global smoothness claim.

## Boundary

This does not cover the true Cauchy pole `P=-4Q`, non-rational geometric
points, or global gluing. It cannot establish a compact curve, smoothness,
genus, a Jacobian, relations, rank, descent, cost, or an ECDLP improvement.

## Reproduction

```bash
sage -python tools/poincare_boundary_singularity_probe_n609g.py --out notes/poincare_boundary_singularity_n609g.json
```
