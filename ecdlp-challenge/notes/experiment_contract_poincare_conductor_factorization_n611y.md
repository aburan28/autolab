# Experiment Contract: N611Y Conductor Factorization of the Degree-Three Class

## Hypothesis

The N611X CM element `u_D=13+4*pi` factors through the conductor-lowering
cubic isogeny from the order `Z[pi]` to the maximal order of discriminant
`-43`. This should materialize the quotient-lattice structure behind the
degree-three elliptic class.

## Calculation

Write `pi=-4+3*omega` with `omega^2-omega+11=0`. Then

```text
u_D = 13+4*pi = 3*(-1+4*omega),
Norm(-1+4*omega)=173.
```

Let `psi:E -> E0` be the ascending degree-three isogeny to

```text
E0: y^2=x^3+45*x+7, j(E0)=10.
```

Let `rho:E0 -> C` be the degree-173 isogeny to
`C:y^2=x^3+51*x+33`, and let
`iota:C -> E0` be `(x,y)->(24^2 x,24^3 y)`. Put
`alpha=iota o rho`. The target identity is

```text
psi^dagger o alpha o psi = [13]+[4]Frob.
```

## Success Criterion

The displayed targets, degrees, scaling isomorphism, CM norm identity, and
pointwise composition identity all replay. This only materializes a lattice
factorization; it does not construct the rank-one class quotient or theta
curve.

## Controls

- the remaining cubic isogeny targets must not have `j=10`;
- the other degree-173 target/scalings must fail the target composition;
- map degrees must multiply to `1557`.

## Falsification Criterion

If no composed map matches the CM endomorphism, retain N611X only as a rank-one
class observation. Passing this does not prove effectivity, irreducibility,
the elliptic quotient associated with `D`, a genus-two equation, evaluation,
relations, rank, target descent, cost advantage, or ECDLP speedup.

## Reproduction

```bash
sage -python tools/poincare_conductor_factorization_n611y.py \
  --out notes/poincare_conductor_factorization_n611y.json
sage -python tools/verify_poincare_conductor_factorization_n611y.py \
  --primary notes/poincare_conductor_factorization_n611y.json \
  --out notes/poincare_conductor_factorization_n611y_independent_verifier.json
```
