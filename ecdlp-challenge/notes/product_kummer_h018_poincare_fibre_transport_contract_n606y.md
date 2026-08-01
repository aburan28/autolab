# Experiment Contract: N606Y Discrete Poincare Fibre Transport

## Hypothesis

On the H018 toy fixture, the Poincare fibre class at `Q` is
`D_Q=8[O]+[zQ]`, with `z=-3-pi`.  For rational `R`, translation on the first
factor satisfies `t_R^*D_Q ~ D_(Q+84R)`, because Frobenius is identity on
`E(F_103)` and `z=[-4]` there.

## Tests

- `h0(D_Q)=9` for sampled rational fibres.
- Exact degree-zero divisor-class equality for the translation law.
- Two-step composition of the induced base shift.

## Boundary

This checks the line-bundle classes needed for a Poincare frame. It does not
construct rational trivializing functions, their scalar cocycle, an isomorphism
to N606U's pushforward bundle, surface sections, factor-base relations, or an
ECDLP improvement.
