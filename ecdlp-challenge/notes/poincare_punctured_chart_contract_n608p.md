# Experiment Contract: N608P target-induced chart on the punctured H018 base

## Hypothesis

The N608O global evaluation morphism gives an explicit matrix-valued frame for
the H018 divisor-model bundle away from the origin when the target
`pi_*O(3Oprime)` is trivialized by its deck-labelled cover fibres. The cover
sections `1`, `x`, and `y` should invert to base-field coefficient vectors in
every nonzero H018 fibre.

## Construction

For a base point `Q`, order its cover fibre as `R+jK` and form the N608J graph
evaluation matrix `V_Q`. For a target cover section `s`, define

```text
c_s(Q) = V_Q^(-1) (s(R), s(R+K), ..., s(R+8K))^T.
```

For adjacent nonzero points in the N606Y `[84]` orbit, use a fixed rational
lift of the base step and set `T_Q = V_(Q+step)^(-1) V_Q`.

## Controls

- Each `V_Q` has rank nine.
- `T_Q` is base-field valued and fixes the constant section coefficient.
- The three target sections have rank three on every tested cover fibre.

## Success Criterion

Every `c_s(Q)` and every transition is Frobenius fixed, with no gauge fitted
point-by-point. This produces an explicit punctured-base chart only.

## Boundary

The origin quotient direction and extension across the origin remain open, so
these are not yet global H018 sections. No source return, relation, rank,
descent, cost, or ECDLP speedup claim follows.

## Reproduction

```bash
sage -python tools/poincare_punctured_chart_n608p.py --out notes/poincare_punctured_chart_n608p.json
sage -python tools/verify_poincare_punctured_chart_n608p.py --primary notes/poincare_punctured_chart_n608p.json --out notes/poincare_punctured_chart_n608p_independent_verifier.json
```
