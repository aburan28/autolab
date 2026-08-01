# N609D Result: uniform formal Q=O boundary pencil

## Handoff: full selected boundary pencil

### Claim or task

Determine whether the N608R-selected regular pencil `x-t*1` has a stable
formal `Q=O` divisor degree as `t` ranges over `F_103`.

### Status

OBSERVATION / UNIFORM REGULARIZED Q=O BOUNDARY PENCIL / STRUCTURAL REPLAY /
MODEL-BOUND / TOY-EVIDENCE / COMPACT RESIDUAL IDENTIFICATION OPEN /
NO_ECDLP_CLAIM

### Assumptions

- N608R selects the regular cover-section span `1,x` in the stated H018
  model, and the N608Q frame is valid at the source origin.
- Boundary coefficients are read in the standard `L(9O)` basis.

### Evidence so far

The primary computation passes seven gates and the structural replay passes
six checks.  The regular unit boundary is exactly
`[1,0,0,0,0,0,0,0,0]`; the x boundary is
`[27,42,0,94,0,56,0,45,0]`.  Therefore every member `x-t*1`, for all 103
values of `t` in `F_103`, retains its nonzero `x^4` coefficient `45`, finite
pole degree eight, and one zero at `P=O` when regarded as a section of
`L(9O)`.

Rational finite zero support is not uniform: its count histogram across the
pencil is `0:59`, `2:35`, `4:8`, `6:1`.  At the declared open levels the
counts are `2,0,0` for `t=0,1,17`, respectively.  The section `y` remains
nonregular and cannot furnish a third boundary direction.

### Failure modes

- The rational finite-zero counts concern only the formal boundary and do not
  match, model, or explain N609A open source-pair counts.
- Stability of the formal degree does not identify a compact global residual
  divisor, its second projection, class, normalization, or a relation law.
- Nothing here supplies a factor base, rank, target descent, charged cost
  comparison, or ECDLP advantage.

### Next concrete action

Construct or identify a compact closure of the N609A residual divisor and
compare its `Q=O` fibre multiplicities with this full formal boundary pencil.

### Artifact paths

- `notes/poincare_origin_pencil_contract_n609d.md`
- `tools/poincare_origin_pencil_n609d.py`
- `notes/poincare_origin_pencil_n609d.json`
- `tools/verify_poincare_origin_pencil_n609d.py`
- `notes/poincare_origin_pencil_n609d_structural_replay.json`
