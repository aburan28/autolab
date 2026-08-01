# Experiment Contract: N607G Theta11 Scalar Transport Gate

## Hypothesis

The missing `p2^*Theta^11` factor in N607E could normalize its raw
nine-dimensional Cauchy monodromy into the transition expected by N606U.

## Null hypothesis

Along a base orbit, `p2^*Theta^11` supplies only scalar fibre transport and
therefore cannot correct a nonscalar matrix monodromy.

## Parameters

- `E/F_103: y^2=x^3+x+24`, rational group order `109`.
- N607E's `[84]` base orbit and recorded nonscalar raw monodromy.
- Local `Theta^11` frame `y*x^4+c*x^5`, selected to be nonzero on the orbit.

## Metrics

- Nonvanishing on all 109 rational orbit points.
- Product of the 109 scalar transition ratios.
- N607E scalar-monodromy status.

## Positive control

The selected degree-11 frame must be nonzero on the full orbit and its scalar
ratios must telescope to one.

## Negative control

The existing N607E monodromy must remain nonscalar after this scalar-only
normalization test.

## Success criterion

The test proves only that Theta11 alone is insufficient; it redirects the
construction to a normalized Poincare/pushforward transition.

## Reproduction

```bash
sage -python tools/product_kummer_h018_theta11_scalar_gate_n607g.py --output notes/product_kummer_h018_theta11_scalar_gate_n607g.json
sage -python tools/verify_product_kummer_h018_theta11_scalar_gate_n607g.py --primary notes/product_kummer_h018_theta11_scalar_gate_n607g.json --output notes/product_kummer_h018_theta11_scalar_gate_n607g_independent_verifier.json
```
