# Experiment Contract: N609E removable P=4Q graph extension

## Hypothesis

The N608Z common factor `x-x(4Q)` is a removable Cauchy-chart zero at
`P=4Q=-(-4Q)`, not a genuine component of the level equation.  The degree-nine
residual retains that graph point exactly when the regularized Cauchy limit has
level `t`.

## Null hypothesis

The regularized value disagrees with the residual quotient, or adding graph
points produces source pairs absent from exhaustive regularized evaluation.

## Parameters

- Registered H018 curve over `F_103`.
- All 108 nonzero base points `Q`.
- Declared levels `t=0,1,17`.
- Support `R=-4Q`, chart-removable point `P=-R=4Q`.

## Regularization

For

```text
k_R(P) = (y(P)+y(R))/(x(P)-x(R)),
```

the value at `P=-R` is the tangent limit

```text
k_R(-R) = -(3*x(R)^2+a)/(2*y(R)).
```

The denominator is nonzero for the declared nonzero prime-order points.

## Metrics

- raw-chart failure count at `P=4Q`;
- agreement of the regularized value with the degree-nine residual root;
- extended source recovery versus exhaustive regularized enumeration;
- source counts before and after restoring the removable points.

## Controls

- Negative: direct raw Cauchy evaluation must fail at every graph point.
- Algebraic: retained residual roots must reconstruct `P=4Q` through the
  unsquared cleared equation.
- Completeness: residual factoring plus the regularized graph values must
  equal exhaustive evaluation on the extended chart.

## Success criterion

Every retained graph intersection is and only is a regularized level source,
and the extended inverse matches exhaustive regularized enumeration at all
declared levels.

## Falsification criterion

Any disagreement leaves the graph excluded and narrows the removable-chart
interpretation.

## Boundary

This restores a removable point of the declared Cauchy chart only. It does not
identify a global compact residual divisor, establish its second projection or
normalization, or provide a factor base, relation law, rank, target descent,
cost advantage, or ECDLP improvement.

## Reproduction

```bash
sage -python tools/poincare_removable_graph_n609e.py --out notes/poincare_removable_graph_n609e.json
sage -python tools/verify_poincare_removable_graph_n609e.py --primary notes/poincare_removable_graph_n609e.json --out notes/poincare_removable_graph_n609e_structural_replay.json
```
