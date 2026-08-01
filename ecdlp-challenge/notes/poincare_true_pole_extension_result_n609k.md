# N609K Result: rational true-pole local extension of the H018 determinant

## Handoff: local Cartier representative at the true-pole graph

### Claim or task

Test whether N609J's determinant has a regular nonvanishing local
representative along the rational nonzero graph `P=-4Q` after the natural
Cauchy-coordinate normalization.

### Status

OBSERVATION / RATIONAL TRUE-POLE LOCAL REPRESENTATIVE FOR N609J DETERMINANT /
STRUCTURAL REPLAY / MODEL-BOUND / TOY-EVIDENCE / GLOBAL CARTIER EXTENSION OPEN /
NO_ECDLP_CLAIM

### Assumptions

- N609J's moving Cauchy determinant represents the selected H018 regular
  chart.
- The natural local parameter at the true-pole graph is
  `u=x(P)-x(-4Q)` for rational nonzero `Q`.
- The H018 fixture has prime order 109, so the rational support has nonzero
  y-coordinate away from the origin.

### Evidence so far

For all 108 rational nonzero `Q` values and levels `0,1,17`, the N609J
leading determinant and moving Cauchy coefficient are nonzero. The local
normalization has exact candidate value

```text
lim_{P -> -4Q} u D_t(P,Q) = -det(V_Q) c_8 2y(-4Q).
```

All 324 values are nonzero. A direct Laurent-series expansion at the
representative `(Q,t)=((1:51:1),0)` agrees with the formula. Setting the
Cauchy coefficient to zero forces the expected zero local limit.

The producer and structural replay each pass six checks.

### Failure modes

- This is a rational nonzero true-pole graph calculation. It does not verify
  the chart intersection at `Q=O` or extension at non-rational geometric
  points.
- It identifies a local representative, but it has not derived Cech
  transition functions proving compatibility with the origin chart or a
  global line bundle.
- It establishes no global Cartier divisor, smoothness, normalization,
  Jacobian, factor base, relation law, rank, descent, costs, or ECDLP result.

### Next concrete action

Derive the overlap ratio between this `uD_t` true-pole frame and the N608R
origin regular frame, then verify its Cech cocycle and the `Q=O` meeting before
testing global Cartier extension or compact-curve geometry.

### Artifact paths

- `notes/poincare_true_pole_extension_contract_n609k.md`
- `tools/poincare_true_pole_extension_n609k.py`
- `notes/poincare_true_pole_extension_n609k.json`
- `tools/verify_poincare_true_pole_extension_n609k.py`
- `notes/poincare_true_pole_extension_n609k_structural_replay.json`
