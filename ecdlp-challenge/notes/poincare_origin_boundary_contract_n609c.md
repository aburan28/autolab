# Experiment Contract: N609C Q=O boundary of the regularized H018 x section

## Hypothesis

The N608R regular origin limit of the selected global cover section `x` gives
the `Q=O` boundary fibre of the `t=0` pencil member. If its limit lies in
`L(8O)` but is viewed as a section of `L(9O)`, the missing pole order produces
one boundary zero at `P=O`, completing the degree-nine residual P-fibre.

## Parameters

- Registered `F_103` H018 fixture.
- N608R formal cover parameter and N608Q source regularization.
- Regularized cover sections `1`, `x`, and `y`.
- Standard `L(9O)` basis
  `1,x,y,x^2,xy,x^3,x^2y,x^4,x^3y`.

## Metrics

- origin regularity and base-field descent of the x-section limit;
- its `L(9O)` coefficient vector;
- finite pole degree and finite zero divisor degree;
- rational finite zero count;
- inferred section zero at `P=O` in `L(9O)`.

## Controls

- Positive: section `x` is regular, descends to `F_103`, and has the recorded
  nonzero `x^4` coefficient.
- Negative: cover section `y` retains a source-origin pole and is not accepted
  as a boundary section.
- Degree control: deleting the `x^4` term lowers the finite pole degree.

## Boundary

This is a formal `Q=O` boundary computation under the N608R global-section
model. It does not yet prove that the compact residual divisor is globally
identified with a particular H018 member, nor does it establish smoothness,
normalization, factor bases, relations, rank, descent, cost, or ECDLP speed.

## Reproduction

```bash
sage -python tools/poincare_origin_boundary_n609c.py --out notes/poincare_origin_boundary_n609c.json
sage -python tools/verify_poincare_origin_boundary_n609c.py --primary notes/poincare_origin_boundary_n609c.json --out notes/poincare_origin_boundary_n609c_structural_replay.json
```
