# N609C Result: regularized Q=O boundary of the selected H018 x section

## Handoff: Q=O boundary fibre candidate

### Claim or task

Compute the origin limit of the N608R-selected cover section `x` and determine
the degree of its boundary divisor in the declared `L(9O)` model.

### Status

OBSERVATION / REGULARIZED Q=O BOUNDARY OF SELECTED H018 X SECTION /
STRUCTURAL REPLAY / MODEL-BOUND / TOY-EVIDENCE / COMPACT RESIDUAL
IDENTIFICATION OPEN / NO_ECDLP_CLAIM

### Assumptions

- The N608R regularized global-section model identifies the selected space with
  the cover span of `1,x` for the registered H018 fixture over `F_103`.
- The formal N608Q source-origin frame is the applicable boundary frame.
- The standard basis of `L(9O)` is
  `1,x,y,x^2,xy,x^3,x^2y,x^4,x^3y`.

### Evidence so far

The primary computation passes nine gates and the structural replay passes six
checks.  The regularized `Q=O` limit of the selected `x` section is

```text
27 + 42x + 94x^2 + 56x^3 + 45x^4.
```

Its coefficient vector in the declared basis is
`[27,42,0,94,0,56,0,45,0]`, so all coefficients descend to `F_103`, the
`x^3y` coefficient is zero, and the `x^4` coefficient is nonzero.  The finite
pole degree is eight; removing the `x^4` term lowers that degree.  Consequently
this element of `L(8O)` has one additional zero at `P=O` when regarded as an
`L(9O)` section.  Its finite zero divisor has degree eight and its rational
finite zeros are `(64 : 32 : 1)` and `(64 : 71 : 1)`.

The cover section `y` remains nonregular at the source origin, providing the
negative control.

### Failure modes

- This is a formal boundary calculation, not an identification of the N609A
  degree-nine residual divisor with a compact global H018 divisor.
- It does not establish the second projection, global divisor class,
  normalization, smoothness, a factor base, relation law, rank, target descent,
  charged cost, or an ECDLP improvement.
- The replay checks preserved the recorded structure but did not independently
  derive the N608R model.

### Next concrete action

Compare the formal `Q=O` degree-nine fibre, including its origin zero, against
the compactification of the N609A residual divisor and determine whether their
P-fibre boundary data agree in a global H018 divisor class.

### Artifact paths

- `notes/poincare_origin_boundary_contract_n609c.md`
- `tools/poincare_origin_boundary_n609c.py`
- `notes/poincare_origin_boundary_n609c.json`
- `tools/verify_poincare_origin_boundary_n609c.py`
- `notes/poincare_origin_boundary_n609c_structural_replay.json`
