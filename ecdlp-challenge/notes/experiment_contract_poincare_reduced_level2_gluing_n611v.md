# Experiment Contract: N611V Reduced Level-Two Gluing Audit

## Hypothesis

The reduced N611U Poincare cross map might induce a non-product coupling on
the full `2`-torsion, selecting a Frey-Kani-style graph plane from the reduced
principal polarization.

## Null Hypothesis

The cross map is zero on both `2`-torsion factors. The resulting level-two
pairing is the direct product symplectic form, so its graph planes are only
abstract choices and do not supply a polarization-selected theta curve.

## Parameters

- field: `F_(103^6)`, the first field containing both the cubic quotient and
  full scalar `2`-torsion;
- curves: N611Q's `E` and `E_g`;
- reduced map: `beta=phi_g o ([12]+[4]Frob)`;
- reduced polarization: `M=[[9,beta^dagger],[beta,345]]`.

## Metrics

- cardinalities of `E[2]` and `E_g[2]`;
- action of `beta` and `beta^dagger` on every torsion point;
- full four-by-four Weil-pairing matrix induced by `M`;
- number of maximal isotropic planes and planes transverse to both factors.

## Positive Control

The unsuppressed map `phi_g` must act nontrivially on `E[2]`.

## Success Criterion

The test distinguishes a genuine nonproduct level-two coupling from the
product-pairing null, records every plane count, and keeps a Frey-Kani theta
construction open only if a selected graph plane exists.

## Falsification Criterion

If `beta` and `beta^dagger` vanish on torsion and the pairing is block
product, this closes only direct level-two selection by the N611U cross map.
It does not rule out higher-level theta data, another quotient lattice,
explicit genus-two construction, relations, rank, descent, cost advantage, or
an ECDLP speedup.

## Reproduction

```bash
sage -python tools/poincare_reduced_level2_gluing_n611v.py \
  --out notes/poincare_reduced_level2_gluing_n611v.json
sage -python tools/verify_poincare_reduced_level2_gluing_n611v.py \
  --primary notes/poincare_reduced_level2_gluing_n611v.json \
  --out notes/poincare_reduced_level2_gluing_n611v_independent_verifier.json
```
