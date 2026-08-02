# Experiment Contract: N611Z Degree-Three Quotient Presentation

## Hypothesis

The N611X rank-one class

```text
D = [[9,beta_D^dagger],[beta_D,346]]
```

is the pullback of the principal polarization on an explicit elliptic target
along a quotient map from `E x E_g`. Concretely, a degree-nine isogeny
`f1:E -> H` should have kernel contained in `ker(beta_D)`.

## Parameters

- field: `F_(103^6)`, containing the cubic two-torsion and the tested
  degree-nine kernel;
- source curve: `E:y^2=x^3+x+24`;
- cross map: `beta_D=phi_g o ([13]+[4]Frob)`;
- degree-nine candidates: every two-step degree-three isogeny path from `E`;
- baseline/control: the neighboring `beta_M=phi_g o ([12]+[4]Frob)` and all
  nonselected degree-nine paths.

## Success Criterion

Exactly one tested degree-nine path has a nine-point kernel annihilated by
`beta_D`. The quotient universal property then gives a unique
`gamma:H -> E_g` with `beta_D=gamma o f1`, of degree `346`. With
`f2=gamma^dagger`, the map

```text
f=(f1,f2): E x E_g -> H
```

has pullback matrix `D` under the stated standard dual-isogeny and principal
polarization identities.

## Falsification Criterion

No selected path, multiple incompatible selected paths, an incomplete
degree-nine kernel, or control-map annihilation invalidates this presentation.
Even a passing presentation does not prove theta-divisor effectivity or
irreducibility, give a genus-two equation/evaluator, or establish relation,
rank, descent, cost, or ECDLP improvement.

## Reproduction

```bash
sage -python tools/poincare_degree3_quotient_presentation_n611z.py \
  --out notes/poincare_degree3_quotient_presentation_n611z.json
sage -python tools/verify_poincare_degree3_quotient_presentation_n611z.py \
  --primary notes/poincare_degree3_quotient_presentation_n611z.json \
  --out notes/poincare_degree3_quotient_presentation_n611z_independent_verifier.json
```
