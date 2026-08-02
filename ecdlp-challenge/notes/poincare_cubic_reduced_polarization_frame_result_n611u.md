# Result: N611U Reduced Cubic Principal-Polarization Frame

## Status

`RESTRICTED THEOREM / KERNEL_PRESERVING_TRIANGULAR_CUBIC_POLARIZATION_REDUCTION /
MODEL-BOUND / TOY-EVIDENCE / THETA_EQUATIONS_OPEN / NO_ECDLP_CLAIM`

## Exact Result

N611S used the cubic quotient frame with `d=pi`. N611U instead uses

```text
d = 2+pi,
Psi_d(P,Q) = (P+dQ, phi_g(Q)).
```

Since `d=pi (mod 2)`, this map kills the same H018 kernel generator
`(pi R,R)` as N611Q. Direct transport gives

```text
(T_d^-1)^dagger H018 T_d^-1
  = [[9, -16-8*pi], [24+8*pi, 690]].
```

Taking

```text
u = [12]+[4]Frob,
beta = phi_g o u
```

gives `beta^dagger phi_g=[-16]-[8]Frob`. Therefore the explicit target
principal form becomes

```text
M_d = [[9,beta^dagger],[beta,345]],
deg(beta)=3104,
det(M_d)=9*345-3104=1.
```

This improves the N611S values `371` and `3338` to `345` and `3104`.

## Restricted Optimality Proof

For every triangular source change with `d=A+B*pi=pi (mod 2)`, `A` is even,
`B` is odd, and the transported lower coefficient is

```text
9A^2 - 45AB + 927B^2 + A - 196B + 11.
```

At `B=1`, the even-`A` minimum is `690` at `A=2`; at `B=-1`, the minimum is
`1078`. After real minimization in `A`, every odd `|B|>=3` has lower bound
greater than `7267`. Thus `d=2+pi` is the exact minimizer in this full
kernel-preserving triangular family.

## Evidence

- N611P/N611Q/N611S receipts bind by independently verified hashes.
- The new quotient kills the designated generator; the `d=0` negative control
  does not.
- The dual-isogeny cross-map identity replays on 22 deterministic points.
- Producer passes eight gates; the independent verifier separately checks the
  reduced matrix, map data, strict improvement, and complete restricted
  minimum argument.

## Interpretation

This is a coordinate reduction, not a new ECDLP algorithm. It makes the next
line-bundle/theta calculation slightly smaller, but supplies no theta equation,
section evaluator, semilinear descent, factor base, relation rank, target
descent, charged cost advantage, or ECDLP speedup.

## Next Concrete Action

Replace the N611T `G_g` bundle frame with the reduced N611U Poincare formula
and compute an evaluable Riemann-Roch generator. If that remains infeasible,
search the full quotient-lattice reduction rather than broadening this
triangular optimum into a claim about all principal-polarization frames.

## Reproduction

```bash
sage -python tools/poincare_cubic_reduced_polarization_frame_n611u.py \
  --out notes/poincare_cubic_reduced_polarization_frame_n611u.json
sage -python tools/verify_poincare_cubic_reduced_polarization_frame_n611u.py \
  --primary notes/poincare_cubic_reduced_polarization_frame_n611u.json \
  --out notes/poincare_cubic_reduced_polarization_frame_n611u_independent_verifier.json
```
