# Experiment Contract: N611R Cubic Quotient Polarization Type

## Hypothesis

Under the N606R isotropic-kernel polarization-descent assumption, the explicit
degree-two quotient `Psi_g` carries the H018 line bundle to a principal
polarization on `E x E_g` over `F_(103^3)`.

## Inputs

- N606R: `G_g` is a maximal isotropic order-two H018 kernel line under the
  stated standard descent theorem;
- N611P: `G_g` is an actual subgroup of `E[2]^2(F_(103^3))` in the H018
  mod-two kernel;
- N611Q: `Psi_g` is an explicit separable quotient of degree two with kernel
  exactly `G_g`;
- N609I/N606N: `H018^2=4`, `det(H018)=2`.

## Calculation

For a finite degree-two quotient `Psi` and a descended line bundle `M` with
`Psi^* M = L_H018`, intersection pullback gives

```text
L_H018^2 = deg(Psi) M^2.
```

Thus `M^2=4/2=2`; on an abelian surface this is the numerical signature of a
principal polarization (`chi(M)=M^2/2=1`).

## Success Criterion

All source receipts bind, the degree/intersection calculation is integral, and
the conclusion is recorded as a `RESTRICTED THEOREM` under the stated descent
assumption.  This does not construct `M`, its theta divisor, or scalar H018
sections.

## Falsification Criterion

If the explicit quotient degree or kernel disagrees with the descent input, or
the intersection calculation is nonintegral, this principal-type conclusion
is unavailable.  No theta, relation, rank, descent, cost, or ECDLP claim may
then proceed from this quotient.

## Reproduction

```bash
python3 tools/poincare_cubic_quotient_polarization_n611r.py \
  --out notes/poincare_cubic_quotient_polarization_n611r.json
python3 tools/verify_poincare_cubic_quotient_polarization_n611r.py \
  --primary notes/poincare_cubic_quotient_polarization_n611r.json \
  --out notes/poincare_cubic_quotient_polarization_n611r_independent_verifier.json
```
