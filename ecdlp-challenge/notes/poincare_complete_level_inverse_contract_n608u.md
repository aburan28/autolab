# Experiment Contract: N608U complete open-level inverse for the H018 pencil

## Hypothesis

For fixed nonzero `Q` and pencil value `t`, the N608S equation
`lambda(P,Q)=t` admits a source generator that does not scan all `P`:

```text
A_Q,t(x) + B_Q(x)y = 0,
A_Q,t(x)^2 - B_Q(x)^2 (x^3+x+24) = 0.
```

The eliminated polynomial has degree at most ten. Factoring it for every `Q`
should recover the complete declared open level, rather than only the sparse
N608T graph subset.

## Parameters

- Registered `F_103` H018 fixture.
- Deterministic level values `t=0,1,17`.
- All 108 nonzero base `Q` values.

## Controls

- Positive: compare every inverse output with exhaustive open-surface
  evaluation.
- Negative: flip the Cauchy numerator contribution in `B_Q` and require a
  set mismatch.

## Metrics

- elimination degree;
- root candidates and recovered points;
- complete-level coverage;
- exceptional `A=B=0` branches;
- charged root-factor calls versus output size.

## Boundary

Passing gives a complete toy source generator for a pencil level. It does not
provide a relation law, factor-base logarithms, matrix rank, descent, or a
sub-rho ECDLP algorithm.

## Reproduction

```bash
sage -python tools/poincare_complete_level_inverse_n608u.py --out notes/poincare_complete_level_inverse_n608u.json
sage -python tools/verify_poincare_complete_level_inverse_n608u.py --primary notes/poincare_complete_level_inverse_n608u.json --out notes/poincare_complete_level_inverse_n608u_independent_verifier.json
```
