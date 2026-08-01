# N609H Result: exact P-fibre compactification of the H018 residual pencil

## Handoff: compact residual P-fibres

### Claim or task

Identify the natural compact P-fibre line bundle of the residual evaluator and
check that its Q-origin limit gives the matching degree-nine boundary fibre.

### Status

RESTRICTED THEOREM / P-FIBRE COMPACTIFICATION OF H018 RESIDUAL PENCIL /
STRUCTURAL REPLAY / MODEL-BOUND / TOY-EVIDENCE / GLOBAL SURFACE CLASS OPEN /
NO_ECDLP_CLAIM

### Assumptions

- The N608S moving Cauchy evaluator is the declared local representation of
  the selected H018 pencil.
- The N608R/N609C regularized origin frame supplies the Q-origin limit.

### Evidence so far

The primary computation passes four gates and the structural replay passes
four checks. For every one of the 108 nonzero base points, both the `x^4`
coefficient and the Cauchy coefficient are nonzero. Since `x^4` has exact
pole order eight at `P=O` and the Cauchy term has a simple pole at `P=-4Q`,
each member `lambda(P,Q)-t` has exact pole divisor

```text
8[O] + [-4Q].
```

Its zero divisor therefore has degree nine in the P-fibre. At `Q=O`, all 103
selected boundary members have finite pole degree eight, nonzero `x^4`
coefficient `45`, zero `x^3y` coefficient, and one extra `P=O` zero when read
as sections of `L(9O)`. Thus the boundary fibre also has degree nine.

### Failure modes

- The calculation compactifies P-fibres only; it does not construct one global
  section on the product surface or establish the Q-projection degree.
- Global divisor class, gluing, normalization, non-rational geometry,
  smoothness, genus, and Jacobian structure remain open.
- No factor base, relation law, rank, target descent, cost advantage, or ECDLP
  improvement follows yet.

### Next concrete action

Glue the P-fibre line bundles into a global surface section, determine its
Q-projection and H018 class, then analyze normalization before any Jacobian or
relation mechanism is considered.

### Artifact paths

- `notes/poincare_residual_pole_compactification_contract_n609h.md`
- `tools/poincare_residual_pole_compactification_n609h.py`
- `notes/poincare_residual_pole_compactification_n609h.json`
- `tools/verify_poincare_residual_pole_compactification_n609h.py`
- `notes/poincare_residual_pole_compactification_n609h_structural_replay.json`
