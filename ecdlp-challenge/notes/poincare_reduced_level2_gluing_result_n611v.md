# Result: N611V Reduced Level-Two Gluing Audit

## Status

`NEGATIVE RESULT / REDUCED_CROSS_MAP_DOES_NOT_SELECT_LEVEL2_FREY_KANI_GRAPH /
MODEL-BOUND / TOY-EVIDENCE / HIGHER_LEVEL_OR_QUOTIENT_LATTICE_OPEN /
NO_ECDLP_CLAIM`

## Exact Result

The first field containing both the N611Q cubic quotient and full scalar
two-torsion is `F_(103^6)`. There, both `E[2]` and `E_g[2]` have four points.

For the reduced N611U map

```text
beta = phi_g o ([12]+[4]Frob),
beta^dagger = ([-8]-[4]Frob) o phi_g^dagger,
```

all coefficients acting before the degree-two maps are even. Thus `beta`
vanishes on every point of `E[2]`, and `beta^dagger` vanishes on every point
of `E_g[2]`.

Consequently the `M=[[9,beta^dagger],[beta,345]]` Weil pairing matrix on
symplectic bases is the direct block product

```text
[[0,1,0,0],
 [1,0,0,0],
 [0,0,0,1],
 [0,0,1,0]].
```

The finite symplectic audit finds all 15 maximal isotropic planes and six
planes transverse to both elliptic factors, exactly as for an undecorated
product. The cross map does not distinguish any one of the six transverse
graph planes. The unsuppressed `phi_g` control acts nontrivially, so the zero
action is specific to the reduced even coefficients rather than an artifact
of the quotient map.

## Interpretation

The N611U principal polarization does not directly select a level-two
Frey-Kani gluing datum. Choosing one of the six graph planes would add an
unjustified external choice, not materialize the theta divisor of the current
polarization. This closes only the direct level-two selection shortcut.

Higher-level theta structures, a full quotient-lattice transformation, an
explicit product presentation, a genuine genus-two curve, a theta evaluator,
scalar descent, relation collection, rank, target descent, sub-rho cost, and
an ECDLP speedup all remain open.

## Next Concrete Action

Search the full quotient-lattice transformation for a product presentation or
a higher-level theta datum that selects one. Only then construct and evaluate
the associated genus-two theta curve before testing any relation mechanism.

## Reproduction

```bash
sage -python tools/poincare_reduced_level2_gluing_n611v.py \
  --out notes/poincare_reduced_level2_gluing_n611v.json
sage -python tools/verify_poincare_reduced_level2_gluing_n611v.py \
  --primary notes/poincare_reduced_level2_gluing_n611v.json \
  --out notes/poincare_reduced_level2_gluing_n611v_independent_verifier.json
```
