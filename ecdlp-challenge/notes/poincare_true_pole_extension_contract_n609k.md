# Experiment Contract: N609K true-pole local extension of the determinant

## Hypothesis

On each rational nonzero `Q` fibre, multiplying the N609J determinant by the
local parameter `u=x(P)-x(-4Q)` removes its sole Cauchy pole and gives a
nonzero local value at `P=-4Q`. Thus the determinant has a regular nonvanishing
local representative across the rational true-pole graph under the declared
Poincare line-bundle trivialization.

## Null hypothesis

The regularized determinant vanishes, retains a pole, or its explicit local
formula disagrees with a Laurent-series replay.

## Parameters

- H018 fixture over `F_103`.
- All 108 nonzero rational base points `Q`.
- Selected levels `t=0,1,17`.
- One formal P-direction Laurent replay at a representative nonzero Q.

## Metrics

- leading nine-by-nine determinant;
- moving Cauchy coefficient;
- regularized local value `-det(V_Q)*c_8*2*y(-4Q)`;
- zero or pole count after multiplication by `u`;
- direct Laurent-series agreement.

## Positive control

The Cauchy coefficient is nonzero in N609H, so the regularized local value
should be nonzero on the declared prime-order rational fibres.

## Negative control

Mutating the moving Cauchy coefficient to zero must force the local limit to
zero.

## Success criterion

All 324 `(Q,t)` rows have nonzero regularized local value, and the formal
series at the representative fibre agrees with the exact formula.

## Falsification criterion

Any zero, residual pole, or formula mismatch blocks the claim that N609J's
determinant provides even a rational-graph local Cartier representative.

## Boundary

This checks a local representative only on the rational nonzero true-pole
graph. It does not prove global Cech compatibility, extension at the meeting
with `Q=O`, non-rational geometric extension, a global Cartier divisor,
smoothness, normalization, a Jacobian, relations, rank, descent, cost, or
ECDLP improvement.

## Reproduction

```bash
sage -python tools/poincare_true_pole_extension_n609k.py --out notes/poincare_true_pole_extension_n609k.json
```
