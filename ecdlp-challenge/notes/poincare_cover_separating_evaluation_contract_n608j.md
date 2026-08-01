# Experiment Contract: N608J cover-separating Poincare evaluation

## Hypothesis

An explicit degree-eleven isogeny `phi:Eprime -> E` that is independent of the
degree-nine cover `pi` separates every degree-nine deck fibre.  Evaluating the
named N606Y basis of `H^0(E,O(8O+[zQ]))` at `phi(pi^-1(Q))` may therefore give
a full-rank, source-noninterpolated fibre map and a concrete starting point for
the missing H018-to-cover intertwiner.

## Null hypothesis

The degree-eleven map fails to separate the deck kernel, or the nine evaluation
functionals are generically singular/undefined.  Then this isogeny-evaluation
ansatz cannot supply the required fibre map.

## Parameters

- H018 target `E/F_103: y^2=x^3+x+24`;
- N606U cover `pi:Eprime -> E` of degree `9`;
- deterministic prime-degree `11` map from `Eprime` to a codomain isomorphic
  to `E`;
- order-nine cover-kernel generator over `F_(103^6)`;
- all nonzero `Q` in `E(F_103)` for which the named N606Y basis is regular.

## Metrics

- degree of `phi` and exact codomain identification;
- order of `phi(K)` for deck generator `K`;
- count and location of regular versus rejected fibres;
- rank and determinant nonvanishing of the nine-by-nine evaluation matrix;
- rank-one control using `pi`, which collapses every deck fibre.

## Positive control

Because `gcd(11,9)=1`, a degree-eleven map must be injective on the cyclic
order-nine deck group.  The `pi` control instead sends the whole deck fibre to
one point and must have evaluation rank one.

## Falsification criterion

Failure of deck separation or a rank below nine on any regular fibre rejects
this fibre-evaluation ansatz.  A full-rank result alone does not establish a
bundle morphism: regularity on charts, target `O(2Oprime)` weights, and
normalized-Poincare compatibility remain separate gates.

## Reproduction

```bash
sage -python tools/poincare_cover_separating_evaluation_n608j.py \
  --out notes/poincare_cover_separating_evaluation_n608j.json
sage -python tools/verify_poincare_cover_separating_evaluation_n608j.py \
  --primary notes/poincare_cover_separating_evaluation_n608j.json \
  --out notes/poincare_cover_separating_evaluation_n608j_independent_verifier.json
```

## Claim boundary

`HYPOTHESIS / COVER-SEPARATING FIBRE EVALUATION / MODEL-BOUND /
TOY-EVIDENCE / NO BUNDLE MORPHISM / NO_ECDLP_CLAIM`.
