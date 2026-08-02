# N611E Result: Line-Ratio Coset Triple Preflight

## Result

`NEGATIVE LINE-RATIO COSET TRIPLE PREFLIGHT / MODEL-BOUND / TOY-EVIDENCE /
NONLINEAR_RATIONAL_MAP_FAMILY_OPEN / NO_ECDLP_CLAIM`

N611E tested the fixed line-ratio factor families

```text
R(P) = L_1(P) / L_2(P),     R(P)^d = c,
```

on `E/F_101: y^2=x^3+x+32`, using four distinct numerator/denominator line
pairs, `d in {5,10,20,25}`, and the two valid `d`th-power image cosets
`c in {1,g^d}`. This corrects the initial invalid `c=g` choice before any
result was recorded: for these divisors of 100, `g` itself is not a
`d`th power.

For every one of the 27 admitted factor bases, every enumerated distinct
zero-sum triple has a nonvertical public line witness and exactly replays

```text
N(x)^d - c*D(x)^d = 0 mod C(x).
```

The independent verifier separately reconstructs the curve points, factor
bases, generic remainders, all 32 matched controls per row, and every source
witness.

## Admission Outcome

| Gate | Outcome |
|---|---|
| Prime-order curve and defined-point controls | pass |
| Exact three-point source replay | pass for all 27 rows |
| Generic support below unordered pair state | fail for all rows |
| Fourfold lift strictly above control maximum | fail for all rows |
| Source-solver admission | fail |

The best generic support ratio is `90 / 28 = 3.214...`, so even the least
dense symbolic source representation already exceeds its unordered pair state.
The largest observed density lift is `6.4x`, but it occurs on an
eight-point base with two triples, exactly tying the 32-control maximum rather
than exceeding it. These small-base fluctuations are therefore not promoted.

## Interpretation

This is an exact source-geometry result, not a source-solver or ECDLP result.
The tested line ratios make full zero-sum triples algebraically visible, but
their denominator-cleared generic membership remains pair-scale or worse
before inversion. The finite family is rejected under the stated controls.

## Limits

The result does not rule out a degree-two rational map, a conic or higher
linear-system factor base, a nontrivial extension quotient, a normalized
Poincare/Jacobian packet, a different source inverse, factor-log rank,
target descent, or a non-generic ECDLP improvement.

## Next Concrete Action

Preflight a factor family whose map is not a ratio of affine line functions.
Its contract must first exhibit a generic three-point remainder support below
the direct pair state, then an exact public source inverse and a nonrandom
relation-density signal before a solver is implemented.

## Handoff: Line-Ratio Coset Closeout

### Claim or task

Determine whether fixed line-ratio cosets provide bounded-support,
nonrandom-density triple relations.

### Status

NEGATIVE RESULT

### Assumptions

- The fixed map and coset family represents the intended rational-function
  class.
- Matched point-sampling controls are a useful finite baseline, not an
  asymptotic theorem.

### Evidence so far

- All exact factor and source identities replay.
- All 27 admitted rows fail the support gate.
- No row exceeds the control maximum at the preregistered fourfold threshold.

### Failure modes

- A higher-degree map may have cancellations absent from affine line ratios.
- A useful source inverse may use a different geometry than an affine source
  line.

### Next concrete action

Construct a degree-two/conic-map preflight with symbolic support, exact
witnesses, matched controls, and a target-return plan.

### Artifact paths

- `notes/experiment_contract_line_ratio_coset_triple_n611e.md`
- `tools/line_ratio_coset_triple_probe_n611e.py`
- `notes/line_ratio_coset_triple_probe_n611e.json`
- `tools/verify_line_ratio_coset_triple_n611e.py`
- `notes/line_ratio_coset_triple_probe_n611e_independent_verifier.json`

