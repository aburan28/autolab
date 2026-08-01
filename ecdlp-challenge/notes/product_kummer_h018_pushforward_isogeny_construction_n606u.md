# N606U: H018 Pushforward Isogeny Construction Preflight

## Status

HYPOTHESIS / POSITIVE CONSTRUCTION INFRASTRUCTURE / BUNDLE_IDENTIFICATION_OPEN /
MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM

## Candidate

Let `L_H018` be the normalized-Poincare line bundle and set

```text
F = p2_* L_H018.
```

The Fourier--Mukai numerical transform of `O_E(9O)` has rank `9`, degree
`-1`.  Pulling it back by `z=-3-pi`, whose norm is `97`, and tensoring by
`O_E(11O)` gives the exact numerical target

```text
rank(F)=9, degree(F)=-97 + 9*11 = 2.
```

The fixture admits rational cyclic degree-nine isogenies.  For the first
two-step composite, its dual isogeny is

```text
pi: E' -> E,
E': y^2 = x^3 + 42x + 28,
deg(pi)=9.
```

On `E'`, `H^0(E',O(2O'))` has Sage-computed basis

```text
1,
(11*y^2 - 50*x*z + z^2)/x^2.
```

Thus a precise construction hypothesis is

```text
F ~= pi_* O_E'(2O').
```

If it holds, the displayed basis gives two named sections of `F`, hence two
candidate H018 sections after an explicit bundle isomorphism.

## Evidence

- `z=-3-pi` has norm `97` in `Z[pi]`.
- Sage found a rational degree-three isogeny from `E`, four further degree-
  three isogenies, and three cyclic degree-nine composites with order-nine
  kernel points; the selected composite's dual has the displayed codomain and
  degree.
- The exact pushforward numerical invariants are consistent with the
  coprime-rank-degree stable-bundle classification on an elliptic curve.

## Open Proof Obligations

1. Prove semihomogeneity/stability and determine `det(F)` with the normalized
   Poincare convention.
2. Compute the determinant of `pi_*O(2O')` and prove it agrees with `F`.
3. Construct an explicit isomorphism, not merely an isomorphism-class
   argument, and evaluate the transported sections on `E x E`.
4. Audit base locus, fibre degrees, source recovery, relation rank, target
   descent, and fully charged cost against rho.

No section evaluator, factor base, relation, descent, or ECDLP improvement is
claimed by this preflight.
