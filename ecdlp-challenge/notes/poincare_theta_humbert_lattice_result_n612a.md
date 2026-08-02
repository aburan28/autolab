# Result: N612A Full Hom Lattice and Humbert Irreducibility Gate

## Status

`RESTRICTED THEOREM / FULL_CM_HOM_LATTICE_AND_GEOMETRIC_THETA_IRREDUCIBILITY /
STANDARD_CM_ENDOMORPHISM_AND_HUMBERT_CRITERION_BOUND / MODEL-BOUND /
TOY-EVIDENCE / CURVE_EQUATION_AND_EVALUATOR_OPEN / NO_ECDLP_CLAIM`

## Exact Result

Write `pi=Frob` with `pi^2+5*pi+103=0`, and assume
`End(E)=Z[pi]` and `phi_g^dagger phi_g=[2]`. For any
`h:E -> E_g`, set `v=phi_g^dagger h`. Then

```text
phi_g o v = [2] o h
```

annihilates all of `E[2]`. Over `F_(103^6)`, the producer materializes all
four points of `E[2]` and checks the mod-two CM multipliers `0`, `1`, `Frob`,
and `1+Frob`. Only zero is killed by `phi_g` on all of `E[2]`. Hence `v` is
even and `h=phi_g u`; under the stated endomorphism assumption this gives

```text
Hom(E,E_g) = phi_g Z[pi].
```

Thus the full integral Neron-Severi lattice is represented by `(a,b,r,s)` with

```text
q_A(a,b,r,s) = a*b - 2*(r^2-5*r*s+103*s^2).
```

For the N611U principal theta class `M=(9,345,12,4)`, an explicit unimodular
reduction of the refined Humbert form gives the positive ternary form

```text
259512*x^2 + 8136*x*y + 699824*x*z
+ 72*y^2 + 10992*y*z + 471817*z^2.
```

Its determinant is `6192`. Exact short-vector enumeration finds no nonzero
value through eight; its minimum is nine, attained at `(4,3,-3)`. Transported
back to the Neron-Severi lattice, this vector is `-D`, where N611X's rank-one
degree-three class has `M.D=3` and refined Humbert value nine.

By the refined Humbert irreducibility criterion, the form does not represent
one. Therefore the principal theta class of `M` is geometrically irreducible
under the stated standard assumptions. This uses Proposition 7 and criterion
(4) in Kani, *Jacobians isomorphic to a product of two elliptic curves and
ternary quadratic forms* (2007). The producer passes 9/9 gates and the
independent Sage replay reconstructs the full-Hom control, quotient lattice,
and minimum.

## Interpretation

The quotient geometry has crossed a meaningful threshold: it now has an
irreducible genus-two theta class rather than merely a rank-one degree-three
elliptic class and quotient target. This establishes no explicit curve model.
The curve's field of definition, equation, degree-three maps, evaluator,
factor base, relations, rank, descent, charged cost, and any ECDLP advantage
remain open.

## Next Concrete Action

Construct an explicit genus-two equation with Jacobian polarized by `M` and
materialize its degree-three elliptic maps to `E` and `E_g`. Then test an
evaluator on exact source/target controls before relation collection.

## Reproduction

```bash
sage -python tools/poincare_theta_humbert_lattice_n612a.py \
  --out notes/poincare_theta_humbert_lattice_n612a.json
sage -python tools/verify_poincare_theta_humbert_lattice_n612a.py \
  --primary notes/poincare_theta_humbert_lattice_n612a.json \
  --out notes/poincare_theta_humbert_lattice_n612a_independent_verifier.json
```
