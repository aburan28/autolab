# N611F Result: H018 Residual Ramification Audit

## Result

`OBSERVATION / EXACT OPEN-FIBER RAMIFICATION MATCH / MODEL-BOUND /
TOY-EVIDENCE / GLOBAL_BRANCH_DIVISOR_OPEN / NO_ECDLP_CLAIM`

N611F audits all 324 degree-nine residual fibers from the selected H018
pencil. For each nonzero base point `Q) and each level `0,1,17), it removes
the certified `P=4Q) graph factor, checks the residual discriminant and
`gcd(f,f')), factors the residual, and reconstructs every repeated
base-field root through the declared `y=-a(x)/b(x)) formula.

| Level | Discriminant-zero fibers | GCD degree | Repeated open roots | N609F P-critical sources |
|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 0 |
| 1 | 0 | 0 | 0 | 0 |
| 17 | 2 | 2 | 2 | 2 |

At level 17 the two exceptional fibers each contain one repeated linear
factor of multiplicity two. Both roots reconstruct valid curve points at the
declared level; neither is the removed graph point. Factor multiplicity excess
equals the derivative-gcd degree in every row. The independent verifier
recomputes the full table and its exact agreement with N609F's formal
P-direction critical-source counts.

## Interpretation

The exact residual multiplicity and independent formal-derivative calculations
are now two views of the same two rational open critical sources. This is a
useful geometric consistency signal for the compact-normalization track.

It is not a global ramification theorem. A compact H018 member would still
need a global Cartier construction, boundary and non-rational branch analysis,
and normalization before its genus or Jacobian can be asserted. No factor
base, relation, rank, descent, cost advantage, or ECDLP speedup follows.

## Next Concrete Action

Construct the global compact member or derive its all-geometric branch divisor.
Only then apply Riemann-Hurwitz and test whether the normalization has an
explicit map, Jacobian decomposition, or packet relation that survives a full
target-to-cost audit.

## Artifacts

- Contract: `notes/experiment_contract_poincare_residual_ramification_n611f.md`
- Producer: `tools/poincare_residual_ramification_probe_n611f.py`
- Receipt: `notes/poincare_residual_ramification_probe_n611f.json`
- Independent verifier:
  `tools/verify_poincare_residual_ramification_n611f.py`
- Verifier receipt:
  `notes/poincare_residual_ramification_probe_n611f_independent_verifier.json`

