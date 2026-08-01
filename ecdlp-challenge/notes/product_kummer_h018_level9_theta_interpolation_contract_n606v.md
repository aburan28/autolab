# Experiment Contract: N606V Level-Nine Theta Interpolation Preflight

## Hypothesis

For the verified cyclic degree-nine isogeny `c:E -> E'` in N606U, evaluation
of the standard basis of `H^0(E,O(9O))` on a finite coset of `ker(c)` gives an
invertible level-nine interpolation frame.  The cyclic ordering of that coset
is diagonalized by a ninth-root discrete Fourier matrix over a field containing
both the kernel and ninth roots of unity.

## Null hypothesis

The evaluation frame is singular, the selected kernel generator is not of
order nine, or the cyclic shift is not diagonalized by the claimed Fourier
matrix.  Any of these failures blocks the proposed explicit theta-function
route before a Poincare cocycle is considered.

## Parameters

- Fixture: `E/F_103: y^2=x^3+x+24`, trace `-5`, CM polynomial
  `pi^2+5*pi+103=0`.
- Isogeny: the first two degree-three Sage isogenies, composed as in N606U.
- Kernel: a point of exact order nine of `ker(c)` over `F_(103^6)`.
- Section basis: `1,x,y,x^2,xy,x^3,x^2y,x^4,x^3y` for `L(9O)`.
- Coset anchor: the first finite base-field point, lifted to `F_(103^6)`.
- Baseline/control: nine successive points from a non-kernel base-field
  generator, which can also give full rank but does not define the isogeny
  kernel frame.

## Metrics

- exact kernel-point order;
- finite-coset check;
- interpolation-matrix rank and determinant nonzero check;
- cyclic row-shift identity;
- Fourier diagonalization identity;
- control-frame rank.

## Positive control

The kernel coset has nine distinct finite points and its evaluation matrix has
rank nine.

## Negative controls

- The kernel coset through `O` is rejected because the chosen affine basis has
  poles there.
- The non-kernel progression may have rank nine, demonstrating that rank by
  itself does not identify the required theta/Poincare cocycle.

## Success criterion

All exact identities pass over `F_(103^6)`, with the control recorded
separately from the kernel construction.

## Falsification criterion

Failure of the kernel order, coset finiteness, rank, shift, or Fourier identity
falsifies this particular explicit-frame route on the registered toy fixture.

## Scope boundary

This is only a level-nine source-curve interpolation preflight.  It does not
construct the normalized-Poincare bundle isomorphism in N606U, transport either
cover-side section, evaluate a section on `E x E`, identify the H018 base locus,
or establish a factor base, relation rank, target descent, sub-rho cost, or an
ECDLP improvement.

## Reproduction

```bash
sage -python tools/product_kummer_h018_level9_theta_interpolation_n606v.py \
  --output notes/product_kummer_h018_level9_theta_interpolation_n606v.json
sage -python tools/verify_product_kummer_h018_level9_theta_interpolation_n606v.py \
  --primary notes/product_kummer_h018_level9_theta_interpolation_n606v.json \
  --output notes/product_kummer_h018_level9_theta_interpolation_n606v_independent_verifier.json
```

## Handoff: Explicit Theta Frame

### Claim or task
Establish the finite level-nine interpolation and Fourier frame needed before
checking its Poincare translation cocycle.

### Status
HYPOTHESIS

### Assumptions
- Standard Riemann--Roch basis for `L(9O)` on this Weierstrass model.
- N606U's deterministic degree-nine isogeny selection.

### Evidence so far
- N606U verifies a cyclic dual degree-nine cover and two cover-side sections.

### Failure modes
- An interpolation frame can be full rank without carrying the normalized
  Poincare cocycle.
- The field extension may make any future descent cost prohibitive.

### Next concrete action
Run the primary and independent N606V computations, then compare the resulting
translation cocycle with the normalized Poincare rigidification.

### Artifact paths
- `tools/product_kummer_h018_level9_theta_interpolation_n606v.py`
- `tools/verify_product_kummer_h018_level9_theta_interpolation_n606v.py`
