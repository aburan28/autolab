# Experiment Contract: REP-AUXLINE-017 Ternary Coset Symmetric Preflight

## Hypothesis

For a coset factor constraint `x^d=c`, a triple can be represented by its
elementary coordinates `(e1,e2,e3)` with exact sparse membership equations.
If the target S4 polynomial has bounded degree in `e3`, a quotient evaluator
may make ternary relations dense enough to improve on REP-AUXLINE-016.

## Null hypothesis

Exact ternary membership necessarily includes a power sum whose elementary
support is quadratic in `d`, recreating pair-scale state before source return.
The bounded S4 degree then does not establish a sub-rho collector.

## Parameters

- symbolic S4 checks on the three prime-order toy curves used by RT1476;
- registered coset fibers: 20-bit `d=190,c=375819` and 22-bit `d=486,c=1`;
- exact membership: `P_d=3c`, `P_2d=3c^2`, `e3^d=c^3`;
- comparison: REP-AUXLINE-016 pair-state dimension and rho lower bound.

## Metrics

- biconditional membership tests on random field triples;
- support of `P_d` and `P_2d` using the exact Girard formula;
- S4 degree in `e3` after exact symmetric reduction;
- support-to-pair-state ratio on registered fibers;
- a pre-implementation relation probability and per-query state gate.

## Success criterion

Advance only if the complete exact membership representation has subquadratic
pre-target support and an explicit collection inequality can plausibly beat
rho after source, rank, descent, and memory costs.

## Falsification criterion

If `P_2d` has quadratic support comparable to or exceeding the pair state,
record the narrow structural negative for this direct ternary elementary-coset
representation and do not implement an S4 source solver.

## Reproduction command

```bash
sage -python tools/compositional_coset_ternary_preflight.py \
  --out notes/compositional_coset_ternary_preflight_rep_auxline017.json
```
