# N609L Result: rational P-origin local extension of the H018 determinant

## Handoff: local P-origin representative for the determinant

### Claim or task

Test whether the N609J determinant has a regular nonvanishing representative
at `P=O` on every rational nonzero Q fibre after its eighth-order origin
normalization.

### Status

OBSERVATION / RATIONAL P-ORIGIN LOCAL REPRESENTATIVE FOR N609J DETERMINANT /
STRUCTURAL REPLAY / MODEL-BOUND / TOY-EVIDENCE / GLOBAL CARTIER EXTENSION OPEN /
NO_ECDLP_CLAIM

### Assumptions

- N609J's moving Cauchy determinant is the declared regular-chart model.
- `v=-x(P)/y(P)` is the standard P-origin parameter.
- The moving basis has pole weights at most eight, with only `x^4` at weight
  eight; the Cauchy term has weight one at P-origin.

### Evidence so far

Across all 108 rational nonzero Q values and levels `0,1,17`, both the leading
determinant and the moving `x^4` coefficient are nonzero. Hence each of the
324 normalized values

```text
lim_{P -> O} v^8 D_t(P,Q) = -det(V_Q)c_7
```

is nonzero. The representative local-parameter replay uses
`v^2x -> 1`, lower pole weight for the other polynomial basis elements, and
`v^8(y+y_s)/(x-x_s) -> 0`; it leaves `x^4` as the sole constant contribution
and matches the formula. Zeroing `c_7` forces the expected zero limit.

The producer and structural replay each pass six checks.

### Failure modes

- The replay is algebraic local-parameter accounting, not an independent
  direct Sage Laurent expansion of the Cauchy quotient at infinity.
- The result covers rational nonzero-Q fibres only and does not cover the
  Q-origin intersection or non-rational geometric points.
- There are still no Cech transition functions, global Cartier proof,
  smoothness/normalization/Jacobian result, factor base, relation law, rank,
  descent, cost comparison, or ECDLP improvement.

### Next concrete action

Derive a common transition convention for the rational P-origin and true-pole
frames and compare both to the N608R Q-origin frame; use that to prove or
refute the global Cartier extension before compact-curve analysis.

### Artifact paths

- `notes/poincare_p_origin_extension_contract_n609l.md`
- `tools/poincare_p_origin_extension_n609l.py`
- `notes/poincare_p_origin_extension_n609l.json`
- `tools/verify_poincare_p_origin_extension_n609l.py`
- `notes/poincare_p_origin_extension_n609l_structural_replay.json`
