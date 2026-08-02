# Result: N611P Cubic Principal-Theta Kernel Coordinates

## Status

`OBSERVATION / CUBIC_H018_PRINCIPAL_THETA_KERNEL_COORDINATES / MODEL-BOUND /
TOY-EVIDENCE / QUOTIENT_EQUATIONS_OPEN / NO_ECDLP_CLAIM`

## Exact Result

Over `K=F_(103^3)`, the three nonzero two-torsion points of

```text
E: y^2=x^3+x+24
```

form a Frobenius orbit `R, R^103, R^(103^2)`.  The three candidate principal
quotient generators are materialized as

```text
g     = (R^103,     R),
h     = (R^(103^2), R^103),
g + h = (R,         R^(103^2)).
```

Each generator has exact order two.  Direct evaluation of the mod-two H018
operator

```text
(P,Q) -> (P + Frob(Q), P + Frob(P) + Q)
```

annihilates all four elements of its kernel.  Frobenius cycles the three
order-two lines

```text
G_g -> G_h -> G_(g+h) -> G_g,
```

and Frobenius cubed fixes every line.  The identity and F3 graph generators
are distinct from `g`.

## Evidence

- Primary construction: all seven coordinate, order, kernel, and Frobenius
  gates pass.
- Independent replay: all five checks pass without importing the producer.
- The serialized coordinates give a reproducible `F_(103^3)` kernel input for
  a quotient implementation.

## Interpretation

This closes the abstract-to-geometric gap for the cubic quotient kernel.  It
does not construct `E^2/G_g`, its principal polarization, any theta divisor,
the three conjugate pullbacks, semilinear descent, scalar H018 sections, a
factor base, a relation matrix, target descent, charged cost, or an ECDLP
speedup.

## Next Concrete Action

Construct the separable quotient of `E^2` by `G_g` over `F_(103^3)` and
materialize a principal polarization or theta divisor on that quotient.  The
result must then be repeated for its two Frobenius conjugates before any
base-field section claim.

## Reproduction

```bash
sage -python tools/poincare_cubic_theta_kernel_n611p.py \
  --out notes/poincare_cubic_theta_kernel_n611p.json
sage -python tools/verify_poincare_cubic_theta_kernel_n611p.py \
  --primary notes/poincare_cubic_theta_kernel_n611p.json \
  --out notes/poincare_cubic_theta_kernel_n611p_independent_verifier.json
```
