# Result: N611T Cubic Theta Line-Bundle Orbit

## Status

`RESTRICTED THEOREM / EXPLICIT_CUBIC_NORMALIZED_POINCARE_LINE_BUNDLE_ORBIT /
RIEMANN_ROCH_DIMENSION_BOUND / MODEL-BOUND / TOY-EVIDENCE /
THETA_EQUATIONS_OPEN / NO_ECDLP_CLAIM`

## Exact Result

For each N611Q cubic quotient factor `phi_i:E -> E_i`, N611T names

```text
u_i    = [21] + [4]Frob,
beta_i = phi_i o u_i,

L_i = p1^*Theta_E^9 tensor p2^*Theta_Ei^371
      tensor (id_E,beta_i^dagger)^*P_E,
```

where `P_E` is the normalized Poincare bundle on `E x E`.  The normalized
Poincare matrix formula gives the polarization matrix

```text
[[9, beta_i^dagger], [beta_i, 371]].
```

For `i=g`, this is the explicit determinant-one principal polarization from
N611S.  The `G_g -> G_h -> G_gh -> G_g` quotient cycle and the corresponding
target curves replay, so the three named line-bundle formulas form the cubic
Frobenius orbit.  Under the stated Poincare and Riemann-Roch assumptions,
each ample principal bundle has `h0(L_i)=1`.

## Evidence

- N611Q and N611S primary/verifier hash bindings pass.
- All three quotient factors have degree two and their recorded target curves
  satisfy the cubic Frobenius cycle.
- The `G_g` Poincare matrix matches N611S's explicit polarization map.
- The three determinant-one, one-dimensional theta budgets replay.
- The product-bundle control `diag(9,371)` has zero cross term and therefore
  cannot represent the N611S polarization.
- The independent verifier separately checks the formula records, cycle,
  determinant budget, and control.

## Interpretation

The cubic principal-polarization target is now an explicit line-bundle orbit,
not just a Hermitian matrix.  This still does not exhibit a nonzero generator
of `H0(L_g)`, an effective theta divisor equation, an evaluator, pullback
functions on `E^2`, semilinear scalar descent, relations, rank, target
descent, a cost model, or any ECDLP speedup.

## Next Concrete Action

Implement Riemann-Roch on `E x E_g` for `L_g` to compute and evaluate its
unique theta generator.  Then transport the two conjugate generators and
solve the resulting semilinear descent equations before attempting scalar
H018 sections.

## Reproduction

```bash
python3 tools/poincare_cubic_theta_line_bundle_n611t.py \
  --out notes/poincare_cubic_theta_line_bundle_n611t.json
python3 tools/verify_poincare_cubic_theta_line_bundle_n611t.py \
  --primary notes/poincare_cubic_theta_line_bundle_n611t.json \
  --out notes/poincare_cubic_theta_line_bundle_n611t_independent_verifier.json
```
