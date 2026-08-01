# N609I Result: conditional global H018 class of the compact residual

## Handoff: global-class target for the residual pencil

### Claim or task

Bind the N609H compact P-fibres to the normalized H018 global-section model
and derive the corresponding Q degree and arithmetic-genus target.

### Status

RESTRICTED THEOREM / CONDITIONAL GLOBAL H018 CLASS OF COMPACT RESIDUAL /
STANDARD-FACT-BOUND / STRUCTURAL REPLAY / MODEL-BOUND / TOY-EVIDENCE /
EXPLICIT GLOBAL EQUATION OPEN / NO_ECDLP_CLAIM

### Assumptions

- N608O's normalized H018 global pushforward model is the intended model.
- N608R's selected cover span `span{1,x}` represents its global sections.
- The N609H local evaluator is the restriction of that selected section.

### Evidence so far

The primary computation passes six gates and the structural replay passes five.
With `z=-3-Frob`, `zbar=2+Frob`, and
`Frob^2+5Frob+103=0`, the normalized line bundle has matrix

```text
[[9, 2+Frob], [-3-Frob, 11]],
```

with determinant two and self-intersection four. On `F_103` rational points,
`z` acts as `[-4]`; hence its P-fibre restriction is

```text
O(9O) tensor O([-4Q]-O) = O(8O+[-4Q]),
```

which is exactly N609H's pole compactification. Under the stated model, the
global residual class therefore has P/Q fibre degrees `9/11`. If an effective
Cartier representative is constructed, its arithmetic genus target is
`1 + 4/2 = 3`.

### Failure modes

- This is conditional on the declared global H018 model; it is not an explicit
  product-surface equation or an independent construction of a Cartier divisor.
- Smoothness, normalization, non-rational geometry, Jacobian structure,
  factor bases, relations, rank, descent, cost, and ECDLP relevance remain open.

### Next concrete action

Construct an explicit global section or equation of this H018 class, prove it
is Cartier, and analyze its true-pole and non-rational loci before using the
genus target for a Jacobian or relation mechanism.

### Artifact paths

- `notes/poincare_global_class_contract_n609i.md`
- `tools/poincare_global_class_n609i.py`
- `notes/poincare_global_class_n609i.json`
- `tools/verify_poincare_global_class_n609i.py`
- `notes/poincare_global_class_n609i_structural_replay.json`
