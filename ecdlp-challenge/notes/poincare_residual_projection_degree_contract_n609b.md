# Experiment Contract: N609B residual fibre-degree certificate

## Hypothesis

For every declared nonzero base point `Q`, the N608Z cleared equation has an
exact order-ten pole at `P=O`. Its common graph zero `P=4Q` therefore leaves a
complete residual divisor of degree nine in the `P` fibre, with no omitted
point at infinity in the affine root recovery.

## Parameters

- Registered `F_103` H018 fixture.
- All 108 nonzero rational base points.
- Levels `t=0,1,17`.
- Cleared equation `A_Q,t(x)+B_Q(x)y=0`.

## Metrics

- nonvanishing of the `c_7 x^5` leading term in `A_Q,t`;
- leading-pole separation from `B_Q y`;
- exact pole order at `P=O`;
- common graph zero and residual divisor degree.

## Controls

- Positive: every recorded raw coefficient `c_7` is nonzero.
- Negative: removing that leading term lowers the computed pole order, showing
  that the degree certificate is not inferred from the eliminant alone.

## Boundary

The result certifies the degree of the residual `P`-fibre on the declared
punctured base. It does not compactify the `Q=O` fibre, prove a global divisor
class, identify the second projection degree, or establish an ECDLP mechanism.

## Reproduction

```bash
sage -python tools/poincare_residual_projection_degree_n609b.py --out notes/poincare_residual_projection_degree_n609b.json
sage -python tools/verify_poincare_residual_projection_degree_n609b.py --primary notes/poincare_residual_projection_degree_n609b.json --out notes/poincare_residual_projection_degree_n609b_structural_replay.json
```
