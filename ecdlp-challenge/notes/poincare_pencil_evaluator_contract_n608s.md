# Experiment Contract: N608S evaluate the selected H018 pencil on the open surface

## Hypothesis

The N608R sections represented by cover functions `1` and `x` define a
computable pencil coordinate on the open surface:

```text
lambda(P,Q) = f_x,Q(P),
```

where `f_x,Q` is recovered by inverting the N608J graph-evaluation matrix.
The coordinate should satisfy the graph identity

```text
lambda(phi(R), pi(R)) = x(R).
```

## Domain

- `Q` is nonzero in `E(F_103)`.
- `P` avoids `O`, `[-4]Q`, and `-[ -4Q ]`, the declared affine Cauchy chart
  exclusions.

## Metrics

- coefficient descent to `F_103`;
- graph identity on all nonzero rational cover points;
- open-surface coverage and value-fibre distribution;
- explicit exclusion count.

## Controls

- Positive: graph identity for every nonzero rational `R`.
- Negative: origin/support points are excluded rather than assigned by an
  arbitrary affine value.

## Boundary

This is an evaluator for the selected pencil on one toy open surface. It does
not prove a smooth genus-three member, factor-base source inverse, relation
rank, descent, charged cost, or ECDLP improvement.

## Reproduction

```bash
sage -python tools/poincare_pencil_evaluator_n608s.py --out notes/poincare_pencil_evaluator_n608s.json
sage -python tools/verify_poincare_pencil_evaluator_n608s.py --primary notes/poincare_pencil_evaluator_n608s.json --out notes/poincare_pencil_evaluator_n608s_independent_verifier.json
```
