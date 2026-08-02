# Experiment Contract: N611W Reduced 2-Adic Cross-Map Structure

## Hypothesis

The N611U cross map may become an isomorphism at a higher two-power level and
therefore directly select a Frey-Kani graph gluing.

## Null Hypothesis

The map has nonunit two-adic determinant. Its Smith factors remain `(4,8)`,
so it can carry higher-level residue information but can never itself be an
isomorphism `E[2^m] -> E_g[2^m]`.

## Calculation

```text
beta = phi_g o ([12]+[4]Frob)
     = phi_g o [4] o ([3]+Frob).
```

`Norm(3+pi)=97` is odd, hence `[3]+Frob` is a two-adic unit. The separable
degree-two factor `phi_g` has two-adic elementary divisors `(1,2)`, giving
`beta` Smith factors `(4,8)`.

The corresponding expected image sizes are `1`, `2`, and `8` on full
`4`-, `8`-, and `16`-torsion.

## Controls

- N611V's full level-two receipt must bind;
- the Frobenius order must be `6,12,24` modulo `4,8,16`;
- all listed extension fields must contain full torsion;
- exact image enumeration must match the Smith-factor prediction.

## Success Criterion

The arithmetic identifies whether `beta` becomes a direct torsion
isomorphism. A nonunit result preserves higher-level residue work but closes
only direct graph selection by this map.

## Falsification Criterion

Any image size larger than the Smith prediction or a unit determinant changes
the gluing branch and requires a new construction. Passing the null does not
rule out quotient-lattice transformations, non-direct higher-level theta
data, theta curves, relations, rank, descent, cost advantage, or ECDLP gains.

## Reproduction

```bash
sage -python tools/poincare_reduced_2adic_crossmap_n611w.py \
  --out notes/poincare_reduced_2adic_crossmap_n611w.json
sage -python tools/verify_poincare_reduced_2adic_crossmap_n611w.py \
  --primary notes/poincare_reduced_2adic_crossmap_n611w.json \
  --out notes/poincare_reduced_2adic_crossmap_n611w_independent_verifier.json
```
