# Experiment Contract: N608L degree audit for the separating evaluator

## Hypothesis

The N608J degree-eleven correspondence `u=(phi,pi):Eprime -> E x E` pulls the
H018 line bundle back to a degree-two line bundle, so the N608K rational map
can extend directly to the N606U target `O_Eprime(2Oprime)`.

## Null hypothesis

The pullback polarization has degree other than two.  Then the N608K evaluator
may be a valid rational local functional but cannot be its desired global
line-bundle morphism without changing the correspondence or adding a genuine
geometric modification.

## Exact Calculation

Use the H018 polarization matrix

```text
H = [[9, zbar], [z, 11]],  z=-3-Frob, zbar=2+Frob,
```

and calculate the induced endomorphism

```text
phi_dual (9 phi + zbar pi) + pi_dual (z phi + 11 pi).
```

The map `A=phi_dual*pi` is compared exactly, as rational maps, with
`[4]+Frob`; duality then gives `pi_dual*phi=[-1]-Frob`.  The Frobenius relation
`Frob^2+5Frob+103=0` determines the induced endomorphism.

## Positive and Negative Controls

- Positive: exact rational-map equality for `A=[4]+Frob`.
- Negative: degree-two target expectation, which must fail if the induced
  polarization is `[3]`.

## Success and Falsification

The direct-global-extension hypothesis passes only when the induced scalar is
`[2]`.  Any other scalar is a scoped rejection of this correspondence as a
direct global extension to `O(2Oprime)`.

## Reproduction

```bash
sage -python tools/poincare_evaluation_degree_n608l.py \
  --out notes/poincare_evaluation_degree_n608l.json
sage -python tools/verify_poincare_evaluation_degree_n608l.py \
  --primary notes/poincare_evaluation_degree_n608l.json \
  --out notes/poincare_evaluation_degree_n608l_independent_verifier.json
```

## Claim Boundary

`NEGATIVE RESULT / DIRECT GLOBAL EXTENSION OF N608K REJECTED / MODEL-BOUND /
TOY-EVIDENCE / OTHER CORRESPONDENCES OR MODIFICATIONS OPEN / NO_ECDLP_CLAIM`.
