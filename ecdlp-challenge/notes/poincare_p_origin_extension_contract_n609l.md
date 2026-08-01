# Experiment Contract: N609L P-origin local extension of the determinant

## Hypothesis

For rational nonzero `Q`, multiplying the N609J determinant by
`v^8`, where `v=-x(P)/y(P)` is the formal parameter at `P=O`, gives a regular
nonzero local value. The value is `-det(V_Q)c_7`, with `c_7` the `x^4`
coefficient in the moving Cauchy frame.

## Null hypothesis

The eighth-order normalization retains a pole, has zero limit, or disagrees
with a direct Laurent-series computation.

## Parameters

- H018 fixture over `F_103`.
- All 108 rational nonzero Q fibres.
- Selected levels `t=0,1,17`.
- One formal P-origin Laurent replay.

## Metrics

- leading determinant and `x^4` coefficient;
- `v^8 D_t` local limit;
- algebraic local-parameter basis replay;
- limit after a zeroed `x^4` coefficient mutation.

## Positive control

N609H reports nonzero `x^4` coefficient on each rational nonzero fibre.

## Negative control

Zeroing the `x^4` coefficient must force the normalized local limit to zero.

## Success criterion

All 324 `(Q,t)` rows have nonzero origin-normalized value, and the algebraic
local-parameter replay matches the exact formula at the representative fibre.

## Falsification criterion

Any residual pole, zero, or formula mismatch blocks a rational P-origin
local representative for the determinant.

## Boundary

This is a rational nonzero-Q local calculation. It does not prove Cech
compatibility with the Q-origin frame, a global Cartier section, non-rational
extension, smoothness, normalization, a Jacobian, relation mechanism, rank,
descent, costs, or ECDLP improvement.

## Reproduction

```bash
sage -python tools/poincare_p_origin_extension_n609l.py --out notes/poincare_p_origin_extension_n609l.json
```
