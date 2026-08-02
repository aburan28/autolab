# N611G Result: H018 Quadratic-Subfield Residual Audit

## Result

`OBSERVATION / EXTENSION-FIELD RESIDUAL RAMIFICATION SIGNAL / MODEL-BOUND /
TOY-EVIDENCE / ALL-GEOMETRIC_BRANCH_DIVISOR_OPEN / NO_ECDLP_CLAIM`

N611G constructs 16 deterministic non-base points in the quadratic subfield of
`F_(103^6)` by first constructing cover points and then applying the declared
degree-nine cover map. This avoids assuming that an arbitrary quadratic target
point has a lift in the chosen ambient field.

For every retained sample, all nine deck translates map back to Q, the
nine-by-nine H018 evaluation matrix has rank nine, the coefficient vector is
quadratic-subfield fixed, and each selected residual has degree nine after
removing the `P=4Q` graph factor. Four candidate cover samples that hit poles
of the auxiliary `phi` map were explicitly skipped before the retained
sample was filled.

| Level | Sum of derivative-gcd degrees across 16 samples |
|---:|---:|
| 0 | 1 |
| 1 | 0 |
| 17 | 2 |

Every observed extension-field multiplicity is a single repeated linear factor
in the ambient field. The independent verifier reconstructs the sample
construction and matches the sample count, pole exclusions, and all three
multiplicity totals.

## Interpretation

The rational N611F critical-point signal is not confined to base-field fibers:
the selected H018 family has repeated residual roots at non-base quadratic Q
values as well. This strengthens the case for deriving a real branch divisor.

It remains a finite sample. It does not enumerate all quadratic fibers, prove
a global Cartier correspondence, determine total ramification, establish a
normalization genus or Jacobian, or supply factor-base relations, rank,
descent, cost improvement, or an ECDLP algorithm.

## Next Concrete Action

Construct a symbolic/global discriminant or otherwise certify the
all-geometric branch divisor of a compact H018 member. Only then use
Riemann-Hurwitz and inspect the normalization/Jacobian for an explicit,
target-returning packet mechanism.

## Artifacts

- Contract:
  `notes/experiment_contract_poincare_quadratic_fiber_ramification_n611g.md`
- Producer: `tools/poincare_quadratic_fiber_ramification_n611g.py`
- Receipt: `notes/poincare_quadratic_fiber_ramification_n611g.json`
- Independent verifier:
  `tools/verify_poincare_quadratic_fiber_ramification_n611g.py`
- Verifier receipt:
  `notes/poincare_quadratic_fiber_ramification_n611g_independent_verifier.json`

