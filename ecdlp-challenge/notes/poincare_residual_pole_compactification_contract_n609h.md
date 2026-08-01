# Experiment Contract: N609H P-fibre compactification of the residual pencil

## Hypothesis

For every nonzero base point `Q`, the evaluated pencil member `lambda(P,Q)-t`
has exact pole divisor

```text
8[O] + [-4Q]
```

in the P variable. Its degree-nine zero divisor is therefore a compact
P-fibre, and its regularized `Q=O` limit is a section of `L(9O)` with one
additional origin zero.

## Null hypothesis

Either leading coefficient vanishes on some base fibre, lowering a pole order,
or the Q-origin limit has incompatible degree/multiplicity.

## Parameters

- Registered H018 fixture over `F_103`.
- All 108 nonzero base points.
- All 103 members of the selected `x-t*1` Q-origin pencil.
- Moving Cauchy basis with support `R=-4Q`.

## Metrics

- nonvanishing of the `x^4` and Cauchy coefficients;
- P-pole orders at `O` and `R=-4Q`;
- Q-origin finite pole degree and P-origin zero multiplicity;
- resulting zero-divisor degree in every P-fibre.

## Controls

- The only order-eight P-origin term is the nonzero `c_7 x^4` term.
- The Cauchy kernel has a simple support pole with nonzero coefficient `c_8`.
- Every Q-origin member retains finite pole degree eight and an `L(9O)` origin
  zero of multiplicity one.

## Success criterion

All nonzero Q fibres have exact pole divisor `8[O]+[-4Q]`; all origin-boundary
members have the compatible degree-nine `L(9O)` zero divisor.

## Falsification criterion

Any coefficient vanishing or boundary mismatch rules out this P-fibre
compactification statement for the declared pencil.

## Boundary

This compactifies individual P-fibres only. It does not construct a global
surface section, determine the Q projection or global H018 divisor class,
analyze normalization/non-rational points, or supply a factor base, relations,
rank, descent, cost advantage, or ECDLP improvement.

## Reproduction

```bash
sage -python tools/poincare_residual_pole_compactification_n609h.py --out notes/poincare_residual_pole_compactification_n609h.json
sage -python tools/verify_poincare_residual_pole_compactification_n609h.py --primary notes/poincare_residual_pole_compactification_n609h.json --out notes/poincare_residual_pole_compactification_n609h_structural_replay.json
```
