# Experiment Contract: N606W Cyclic-Cover Fibre Evaluator

## Hypothesis

The N606U cyclic cover `pi:E' -> E` has an explicitly evaluable nine-point
fibre over a nonzero `F_103` point.  The cover-side basis `1,x` of
`H^0(E',O(2O'))` evaluates with rank two on that fibre and carries an exact
Frobenius permutation over `F_(103^6)`.

## Null hypothesis

The degree-nine kernel does not generate the selected fibre, rational-map
evaluation disagrees with the isogeny, the fibre is not Frobenius-stable, or
the two cover sections become dependent on the fibre.

## Parameters and controls

- N606U fixture and deterministic dual degree-nine isogeny.
- One nonzero rational target point and its rational cover preimage.
- A degree-nine kernel generator over `F_(103^6)`; the quadratic field has
  only the order-three kernel point and is a field-cost control.
- Cover section basis `1,x`.
- Negative control: the fibre over the origin contains the pole of `x`, so it
  is recorded as unsuitable for this affine evaluator.

## Metrics

- map degree and kernel-generator order;
- nine distinct finite fibre points;
- agreement of explicit rational maps with the isogeny;
- evaluation rank;
- Frobenius multiplier and fibre permutation;
- Frobenius equivariance of the two base-field sections.

## Promotion boundary

This builds only the cover-side evaluation vector.  It does not construct an
isomorphism from `p2_*L_H018`, so it does not evaluate either H018 surface
section or establish a base locus, factor base, relation rank, target descent,
sub-rho cost, or ECDLP improvement.

## Reproduction

```bash
sage -python tools/product_kummer_h018_cover_fiber_evaluator_n606w.py \
  --output notes/product_kummer_h018_cover_fiber_evaluator_n606w.json
sage -python tools/verify_product_kummer_h018_cover_fiber_evaluator_n606w.py \
  --primary notes/product_kummer_h018_cover_fiber_evaluator_n606w.json \
  --output notes/product_kummer_h018_cover_fiber_evaluator_n606w_independent_verifier.json
```
