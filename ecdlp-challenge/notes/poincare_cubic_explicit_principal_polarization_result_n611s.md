# Result: N611S Cubic Explicit Principal Polarization

## Status

`RESTRICTED THEOREM / EXPLICIT_CUBIC_PRINCIPAL_POLARIZATION_HOMOMORPHISM /
STANDARD_DUAL_ISOGENY_IDENTITY_BOUND / MODEL-BOUND / TOY-EVIDENCE /
THETA_DIVISOR_OPEN / NO_ECDLP_CLAIM`

## Exact Result

Let `pi` be cubic Frobenius on the N611Q source curve, so
`pi^2+5*pi+103=0`.  The H018 form is

```text
H = [[9, 2+pi], [-3-pi, 11]].
```

For the N611Q quotient factor `phi_g:E -> E_g`, define

```text
u    = [21] + [4]Frob,
beta = phi_g o u.
```

The dual identity `phi_g^dagger phi_g=[2]` and
`u^dagger=[1]-[4]Frob` give

```text
beta^dagger phi_g = [2]-[8]Frob.
```

With `T=[[1,pi],[0,1]]`, exact Hermitian transport gives

```text
(T^-1)^dagger H T^-1
  = [[9, 2-8*pi], [42+8*pi, 742]].
```

Therefore the target map

```text
M = [[9, beta^dagger], [beta, 371]]
```

has the required pullback transport through
`Psi_g=diag(1,phi_g) o T`.  Since

```text
deg(beta) = 2*Norm(21+4*pi) = 3338,
det(M) = 9*371 - 3338 = 1,
```

and the leading coefficient is positive, `M` is an explicit principal
polarization homomorphism in the stated CM-product model.

## Evidence

- The producer binds N609I, N611Q, and N611R hashes and passes eight gates.
- It replays `beta^dagger phi_g=[2]-[8]Frob` on 18 deterministic nonzero
  base-field multiples, all three cubic two-torsion points, and the origin.
- The independent Sage verifier recomputes those 22 point identities, rejects
  the identity-map off-diagonal control, and separately checks the CM witness,
  transported form, degree, and determinant.

## Interpretation

This removes a specific ambiguity in N611R: the descended principal object is
now represented by concrete isogeny and Frobenius maps, rather than only by
its intersection number.  It does **not** supply an associated line bundle,
an effective theta divisor, quotient theta equations, the two conjugate
polarizations with compatible divisor data, semilinear descent, scalar H018
sections, or any factor-base, relation, rank, target-descent, cost, or ECDLP
speedup claim.

## Next Concrete Action

Construct an explicit ample line bundle or theta divisor realizing `M` on
`E x E_g`; then pull it back through `Psi_g`, build the two Frobenius
conjugates, and test whether their descent datum yields named scalar H018
sections.

## Reproduction

```bash
sage -python tools/poincare_cubic_explicit_principal_polarization_n611s.py \
  --out notes/poincare_cubic_explicit_principal_polarization_n611s.json
sage -python tools/verify_poincare_cubic_explicit_principal_polarization_n611s.py \
  --primary notes/poincare_cubic_explicit_principal_polarization_n611s.json \
  --out notes/poincare_cubic_explicit_principal_polarization_n611s_independent_verifier.json
```
