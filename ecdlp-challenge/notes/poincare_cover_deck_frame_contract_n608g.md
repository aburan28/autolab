# Experiment Contract: N608G H018 cover deck-frame audit

## Hypothesis

The N606U cyclic degree-nine cover has a concrete deck-labelled local frame on
the nonzero H018 base orbit.  That frame supplies a geometric pushforward-side
reference for a future explicit intertwiner; it must not be replaced by a
finite-orbit interpolation gauge.

## Null hypothesis

The selected degree-nine isogeny does not admit the claimed cyclic lift,
kernel labelling, or translation/Frobenius compatibility on the registered
finite-field orbit.  In that case it cannot support the proposed local
pushforward frame.

## Parameters

- H018 target curve `E/F_103: y^2=x^3+x+24`;
- deterministic cyclic degree-nine dual cover from N606U;
- base shift `[84]G` on the order-109 target group;
- a degree-nine kernel generator over `F_(103^6)`;
- the 108 nonzero target orbit points, excluding the origin fibre where the
  constant `O(2O')` trivialization is not used.

## Metrics

- order of the lifted base shift and the deck generator;
- exact cover-map image of all `108 * 9` labelled points;
- injectivity of the labelled frame;
- deck and lifted-base translation permutation laws;
- Frobenius permutation law on each labelled fibre;
- absence of the cover origin on the regular base orbit.

## Positive control

The degree-nine kernel generator cycles each fibre, while the lifted
order-109 base shift advances the fibre index and the cover rational map sends
every labelled point to its intended target point.

## Negative control

The origin target fibre is excluded: it contains the cover origin and therefore
cannot use the asserted constant local trivialization of `O(2O')`.

## Success criterion

All labelled regular-fibre checks pass exactly.  This constructs only a
pushforward-side local frame, not an H018 section frame or a bundle
isomorphism.

## Falsification criterion

Any failed cover-map, order, injectivity, translation, Frobenius, or regularity
check rejects this frame construction at the stated finite-field scope.

## Reproduction

```bash
sage -python tools/poincare_cover_deck_frame_n608g.py \
  --out notes/poincare_cover_deck_frame_n608g.json
sage -python tools/verify_poincare_cover_deck_frame_n608g.py \
  --primary notes/poincare_cover_deck_frame_n608g.json \
  --out notes/poincare_cover_deck_frame_n608g_independent_verifier.json
```

## Claim boundary

`OBSERVATION / COVER_SIDE_LOCAL_FRAME / MODEL-BOUND / TOY-EVIDENCE /
NO_EXPLICIT_H018_INTERTWINER / NO_ECDLP_CLAIM`.
