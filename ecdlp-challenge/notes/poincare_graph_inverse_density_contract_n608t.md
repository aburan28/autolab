# Experiment Contract: N608T graph-restricted pencil inverse density

## Hypothesis

The exact graph inverse from N608S,

```text
x(R)=t -> (phi(R), pi(R)),
```

occupies a large enough fraction of each pencil level `lambda=t` to seed a
surface relation collector without a separate complete-level inverse.

## Null Hypothesis

The graph is only a one-dimensional subset of the two-dimensional pencil
surface, so its exact inverse is too sparse for direct relation collection.

## Metrics

- graph sources per reachable pencil value;
- full open-surface fibre counts from N608S;
- graph fraction per reachable level and total open density.

## Success Criterion

The graph source fraction must remain bounded below by a substantial constant
across its reachable values. Fractions vanishing at toy scale reject this
graph-only source-return mechanism.

## Boundary

This tests only the graph-restricted inverse. A complete projective level
curve, another correspondence, source recovery, rank, descent, and charged
cost remain separate questions.

## Reproduction

```bash
sage -python tools/poincare_graph_inverse_density_n608t.py --out notes/poincare_graph_inverse_density_n608t.json
sage -python tools/verify_poincare_graph_inverse_density_n608t.py --primary notes/poincare_graph_inverse_density_n608t.json --out notes/poincare_graph_inverse_density_n608t_independent_verifier.json
```
