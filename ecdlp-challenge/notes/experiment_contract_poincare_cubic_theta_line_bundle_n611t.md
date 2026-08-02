# Experiment Contract: N611T Cubic Theta Line-Bundle Orbit

## Hypothesis

The explicit N611S polarization maps determine three concrete normalized-
Poincare line-bundle formulas on the cubic quotient surfaces, one for each
Frobenius-conjugate H018 kernel line.  Each formula has a one-dimensional
theta space, although this experiment does not compute a theta equation.

## Construction

For the N611Q factor `phi_i:E -> E_i`, put

```text
u_i = [21] + [4]Frob,
beta_i = phi_i o u_i,
L_i = p1^*Theta_E^9 tensor p2^*Theta_Ei^371
      tensor (id_E, beta_i^dagger)^*P_E.
```

Here `P_E` is the normalized Poincare bundle on `E x E`.  The standard
Poincare matrix formula assigns

```text
[[9, beta_i^dagger], [beta_i, 371]]
```

to `L_i`.  N611S proves this matrix has determinant one for `i=g`; the same
dual-isogeny calculation applies to the Frobenius conjugates.

## Controls and Metrics

- N611S and N611Q independently replayed receipts must bind;
- all three factors must have degree two and form the stated Frobenius orbit;
- the exact Poincare formula must recover the N611S matrix entries;
- the determinant and conditional Riemann-Roch dimension must equal one;
- a product-bundle control without the Poincare term must have zero
  off-diagonal, and therefore must not match the polarization.

## Success Criterion

All three bundle formulas, their Frobenius cycle, and their determinant-one
theta-space budget replay.  The result is a line-bundle-level construction,
not an explicit section or divisor equation.

## Falsification Criterion

If the Poincare cross term, determinant, or Frobenius cycle fails, the theta
target remains only N611S's polarization-homomorphism data.  Passing this
does not compute an effective divisor, a section evaluator, scalar descent,
relations, rank, target descent, cost advantage, or ECDLP speedup.

## Reproduction

```bash
python3 tools/poincare_cubic_theta_line_bundle_n611t.py \
  --out notes/poincare_cubic_theta_line_bundle_n611t.json
python3 tools/verify_poincare_cubic_theta_line_bundle_n611t.py \
  --primary notes/poincare_cubic_theta_line_bundle_n611t.json \
  --out notes/poincare_cubic_theta_line_bundle_n611t_independent_verifier.json
```
