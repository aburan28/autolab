# Experiment Contract: N608Z base graph in denominator-cleared H018 levels

## Hypothesis

For a fixed base point `Q`, clearing the moving Cauchy denominator in
`lambda(P,Q)=t` introduces a point common to every level. If that point is the
graph `P=4Q`, the degree-ten N608U eliminant has a removable common degree-one
factor and the residual fibre has degree nine.

## Parameters

- Registered `F_103` H018 fixture.
- All 108 nonzero base points.
- Levels `t=0,1,17`.
- N608U denominator-cleared equation `A_Q,t(x)+B_Q(x)y=0` and eliminant
  `F_Q,t(x)=A_Q,t(x)^2-B_Q(x)^2(x^3+x+24)`.

## Metrics

- exact equation evaluation at `P=4Q` and at `P=-4Q`;
- divisibility of `F_Q,t` by `x-x(4Q)`;
- degree after quotienting the common factor;
- multiplicity distribution of the common x-root;
- a sign mutation of the Cauchy numerator contribution.

## Controls

- Positive: `P=4Q` satisfies the unsquared cleared equation for every row.
- Negative: `P=-4Q` does not satisfy that equation, even though it has the
  same x-coordinate; this distinguishes the actual graph from an artefact of
  squaring the eliminant.
- Mutation: reverse the Cauchy numerator sign and require the `P=4Q` graph
  test to fail.

## Boundary

The quotient degree is a fibrewise statement on the declared open chart. It
does not prove a global divisor decomposition, smoothness, normalization,
factor-base relation, rank, descent, cost, or ECDLP improvement.

## Reproduction

```bash
sage -python tools/poincare_cleared_base_graph_n608z.py --out notes/poincare_cleared_base_graph_n608z.json
sage -python tools/verify_poincare_cleared_base_graph_n608z.py --primary notes/poincare_cleared_base_graph_n608z.json --out notes/poincare_cleared_base_graph_n608z_structural_replay.json
```
