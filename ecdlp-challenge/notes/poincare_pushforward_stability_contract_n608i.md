# Experiment Contract: N608I H018 pushforward stability and determinant audit

## Hypothesis

The N606U rank-nine degree-two H018 pushforward and the explicit degree-nine
cover pushforward have the same determinant and are stable.  Standard elliptic
bundle classification would then make their abstract isomorphism an actual
restricted theorem, rather than a numerical identification conditional on
unchecked stability.

## Null hypothesis

One of the required numerical, separability, determinant, or stability inputs
fails on the registered `F_103` fixture.  Then the N606U target cannot be
upgraded through coprime rank-degree classification.

## Parameters

- H018 curve `E/F_103: y^2=x^3+x+24`;
- explicit CM map `alpha=3+pi`, with `deg(alpha)=97` and `z=-alpha`;
- normalized Fourier--Mukai input `W_9` of rank `9`, degree `-1`, and
  determinant `O(-O)`;
- explicit cyclic degree-nine cover `pi:E' -> E` from N606U;
- cover line bundle `O_E'(2O')`.

## Metrics

- exact degrees, separability, and coprime rank-degree pairs;
- determinant arithmetic for both bundles;
- cyclic-kernel oddness and the translation-stabilizer check for `O(2O')`;
- explicit distinction between an abstract isomorphism and a matrix/evaluator.

## Positive control

The known cover has a cyclic order-nine kernel, and its nine translate line
bundles are distinct because no nonzero kernel point is two-torsion.

## Negative control

Replacing the odd degree-nine kernel by an even-order model would invalidate
the translation-stabilizer argument; that argument must not be silently reused
outside the odd cyclic cover.

## Success criterion

All fixture arithmetic gates hold and the proof dependencies reduce to named
standard theorems: Fourier--Mukai numerics/stability, separable pullback
semistability, and coprime rank-degree classification on an elliptic curve.

## Falsification criterion

Any mismatch in rank, degree, determinant, separability, kernel parity, or
translation stabilizer rejects the theorem upgrade at this fixture.

## Reproduction

```bash
sage -python tools/poincare_pushforward_stability_n608i.py \
  --out notes/poincare_pushforward_stability_n608i.json
sage -python tools/verify_poincare_pushforward_stability_n608i.py \
  --primary notes/poincare_pushforward_stability_n608i.json \
  --out notes/poincare_pushforward_stability_n608i_independent_verifier.json
```

## Claim boundary

`RESTRICTED THEOREM / ABSTRACT BUNDLE ISOMORPHISM / STANDARD-THEOREM-BOUND /
MODEL-BOUND / TOY-EVIDENCE / NO EXPLICIT INTERTWINER / NO_ECDLP_CLAIM`.
