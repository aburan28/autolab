# Experiment Contract: N606U H018 Pushforward Identification

## Hypothesis

For the normalized H018 Poincare bundle `L`, the rank-nine degree-two bundle
`F=p2_*L` is isomorphic to `pi_*O_E'(2O')` for an explicit cyclic degree-nine
isogeny `pi:E' -> E`.  Its cover-side Riemann--Roch basis would then give an
abstract, named basis of `H^0(E x E,L)`.

## Literature-Derived Assumptions

- Normalized Fourier--Mukai on an elliptic curve sends `O(9O)` to a stable
  bundle of rank `9`, degree `-1`, and normalized determinant `O(-O)`.
- Pullback by the separable CM endomorphism `z=-3-pi` and tensoring by
  `O(11O)` preserve the stated rank/degree and normalized determinant.
- A stable elliptic bundle with coprime rank and degree is determined by rank,
  degree, and determinant (Atiyah, 1957).
- For an odd cyclic isogeny, `det(pi_*O)=O`; the norm formula gives
  `det(pi_*O(2O'))=O(2O)`.

## Minimal Test

Construct a rational cyclic degree-nine isogeny on the H018 fixture, verify
its kernel has order-nine points, compute `H^0(E',O(2O'))`, and compare all
numerical and determinant invariants under the assumptions above.

## Promotion Boundary

Passing identifies only the isomorphism class conditionally on the cited
bundle facts.  It does not construct the isomorphism `F -> pi_*O(2O')`,
evaluate a surface section, define a factor base, or establish an ECDLP
speedup.

## Reproduction

```bash
sage -python tools/product_kummer_h018_pushforward_identification_n606u.py \
  --output notes/product_kummer_h018_pushforward_identification_n606u.json
sage -python tools/verify_product_kummer_h018_pushforward_identification_n606u.py \
  --primary notes/product_kummer_h018_pushforward_identification_n606u.json \
  --output notes/product_kummer_h018_pushforward_identification_n606u_independent_verifier.json
```
