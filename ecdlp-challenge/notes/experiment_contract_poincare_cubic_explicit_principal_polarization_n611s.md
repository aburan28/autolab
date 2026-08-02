# Experiment Contract: N611S Cubic Explicit Principal Polarization

## Hypothesis

The N611Q cubic quotient `Psi_g:E^2 -> E x E_g` carries the H018 Hermitian
form to an explicit principal polarization homomorphism on `E x E_g`, not
merely to a numerical principal type.  In the cubic toy fixture, its
off-diagonal map is

```text
beta = phi_g o ([21] + [4]Frob): E -> E_g.
```

## Inputs

- N609I: the H018 form has determinant two and self-intersection four;
- N611Q: `phi_g:E -> E_g` is the explicit degree-two factor of `Psi_g`;
- N611R: the quotient has the stated numerical principal-polarization type.

## Model and Calculation

Write `pi` for cubic Frobenius, with
`pi^2 + 5*pi + 103 = 0`.  The source form is

```text
H = [[9, 2+pi], [-3-pi, 11]].
```

For `T=[[1,pi],[0,1]]`, direct Hermitian transport gives

```text
K = (T^-1)^dagger H T^-1
  = [[9, 2-8*pi], [42+8*pi, 742]].
```

Since `phi_g^dagger phi_g=[2]`, take
`u=[21]+[4]pi`, so `u^dagger=[1]-[4]pi`, and define
`beta=phi_g o u`.  Then

```text
beta^dagger phi_g = [2]u^dagger = [2]-[8]pi.
```

Thus the explicit target polarization matrix is

```text
M = [[9, beta^dagger], [beta, 371]],
```

and `D^dagger M D=K` for `D=diag(1,phi_g)`.

## Metrics

- exact CM matrix identities;
- exact degree and norm identities;
- deterministic pointwise replay of `beta^dagger phi_g=[2]-[8]Frob`;
- Hermitian determinant and positivity.

## Positive Control

The dual-isogeny identity `phi_g^dagger phi_g=[2]` must replay on base-field
and cubic two-torsion points.

## Negative Control

Replacing `u` with the identity must fail the required off-diagonal identity
`beta^dagger phi_g=[2]-[8]Frob`.

## Success Criterion

All predecessor receipts bind, the explicit map identity and CM transport
replay, and `M` has positive leading coefficient and determinant one.

## Falsification Criterion

Any failed dual-map, transport, degree, or determinant check leaves only the
N611R numerical statement.  Passing this construction does not construct an
associated line bundle, theta divisor, quotient theta equations, Frobenius
descent datum, scalar sections, relation system, rank argument, target
descent, cost advantage, or ECDLP speedup.

## Reproduction

```bash
sage -python tools/poincare_cubic_explicit_principal_polarization_n611s.py \
  --out notes/poincare_cubic_explicit_principal_polarization_n611s.json
sage -python tools/verify_poincare_cubic_explicit_principal_polarization_n611s.py \
  --primary notes/poincare_cubic_explicit_principal_polarization_n611s.json \
  --out notes/poincare_cubic_explicit_principal_polarization_n611s_independent_verifier.json
```
