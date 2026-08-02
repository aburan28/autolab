# N611H Result: Determinant-Weighted H018 Coefficient Interpolation

## Result

`NEGATIVE DETERMINANT-SCALAR NORMALIZATION RESULT / MODEL-BOUND /
TOY-EVIDENCE / VECTOR_VALUED_TRANSITION_OPEN / NO_ECDLP_CLAIM`

For every nonzero rational Q, N611H reconstructs the N609J nine-by-nine
moving-basis matrix V_Q, checks its determinant is nonzero, base-field fixed,
and invariant under a deck-cycle row shift, then multiplies each N608V raw
coefficient by det(V_Q).

The resulting weighted coefficients have exactly the same first interpolation
orders as the raw frame:

| Coefficient indices | First pole order |
|---|---:|
| 0, 1, 3, 5, 7 | 106 |
| 2, 4, 6, 8 | 109 |

None fits through L(24O), and no coefficient improves on the raw 106-pole
barrier. The independent verifier reconstructs all determinants and the full
weighted interpolation table.

## Interpretation

The N609J determinant is a correct deck-invariant scalar on the rational
chart, but scalar multiplication alone does not normalize the moving Cauchy
coefficient vector into a low-pole global law. This rejects this specific
scalar-gauge route, not the existing non-scalar Poincare transition framework.

## Next Concrete Action

Use the certified vector-valued transition matrices directly: construct a
finite normalized module/equation whose entries transform by those matrices,
rather than attempting another scalar rescaling. It must expose an explicit
global section, branch data, source inverse, target return, rank, descent, and
charged cost before any ECDLP claim.

## Artifacts

- Contract:
  `notes/experiment_contract_poincare_determinant_weighted_interpolation_n611h.md`
- Producer: `tools/poincare_determinant_weighted_interpolation_n611h.py`
- Receipt: `notes/poincare_determinant_weighted_interpolation_n611h.json`
- Independent verifier:
  `tools/verify_poincare_determinant_weighted_interpolation_n611h.py`
- Verifier receipt:
  `notes/poincare_determinant_weighted_interpolation_n611h_independent_verifier.json`

