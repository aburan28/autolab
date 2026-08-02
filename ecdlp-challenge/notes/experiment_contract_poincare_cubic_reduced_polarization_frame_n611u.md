# Experiment Contract: N611U Reduced Cubic Principal-Polarization Frame

## Hypothesis

The N611S cubic principal polarization has a smaller explicit representative
after a kernel-preserving triangular source change.  This can reduce the
degree of the Poincare cross map before attempting theta evaluation.

## Construction

For `d=2+pi`, set

```text
T_d = [[1,d],[0,1]],
Psi_d(P,Q) = (P+dQ, phi_g(Q)).
```

Because `d=pi (mod 2)`, `Psi_d` kills exactly the same N611P line
`G_g=< (pi R,R) >` as N611Q.  Transporting

```text
H018 = [[9,2+pi],[-3-pi,11]]
```

gives

```text
K_d = (T_d^-1)^dagger H018 T_d^-1
    = [[9,-16-8*pi],[24+8*pi,690]].
```

Let `u=[12]+[4]Frob` and `beta=phi_g o u`.  Then
`beta^dagger phi_g=[-16]-[8]Frob`, so the target form is

```text
M_d = [[9,beta^dagger],[beta,345]].
```

## Optimality Scope

This is an exact reduction only within triangular changes
`T_d` satisfying `d=pi (mod 2)`.  Writing `d=A+B*pi`, that condition means
`A` is even and `B` is odd; the lower transported coefficient is

```text
b_d = 9A^2 - 45AB + 927B^2 + A - 196B + 11.
```

The integer minimum in this family is `690`, attained at `(A,B)=(2,1)`.

## Controls

- N611Q/N611S receipts must bind;
- the old `d=pi` frame must replay lower coefficient `742`;
- a non-kernel-preserving `d=0` control must fail the H018 kernel condition;
- the explicit dual-isogeny identity must replay on deterministic points.

## Success Criterion

The new quotient has the same exact kernel, the map identity and matrix
transport pass, determinant is one, and the complete triangular-family lower
coefficient proof certifies the claimed minimum.

## Falsification Criterion

If any map, kernel, arithmetic, or minimization gate fails, retain N611S
without a reduced-frame claim.  Passing this coordinate reduction does not
construct a theta equation, sections, relations, rank, target descent, a cost
advantage, or an ECDLP speedup.

## Reproduction

```bash
sage -python tools/poincare_cubic_reduced_polarization_frame_n611u.py \
  --out notes/poincare_cubic_reduced_polarization_frame_n611u.json
sage -python tools/verify_poincare_cubic_reduced_polarization_frame_n611u.py \
  --primary notes/poincare_cubic_reduced_polarization_frame_n611u.json \
  --out notes/poincare_cubic_reduced_polarization_frame_n611u_independent_verifier.json
```
