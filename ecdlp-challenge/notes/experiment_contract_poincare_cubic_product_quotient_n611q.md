# Experiment Contract: N611Q Cubic Product Quotient Models

## Hypothesis

Each cubic H018 quotient line can be realized as the kernel of an explicit
separable degree-two map from `E^2` to a product `E x E_i` over `F_(103^3)`.

## Construction

Let `R,R^103,R^(103^2)` be the nonzero two-torsion orbit.  The unipotent
automorphism

```text
T(P,Q) = (P + Frob(Q), Q)
```

sends the three H018 line generators to `(O,R_i)`.  For the Velu isogeny
`phi_i:E -> E_i` with kernel `<R_i>`, define

```text
Psi_i(P,Q) = (P + Frob(Q), phi_i(Q)).
```

The declared expected kernels are `G_g,G_h,G_(g+h)` in Frobenius order.

## Metrics

- degree and exact factor-kernel points of each `phi_i`;
- exact images of all three candidate H018 generators under each `Psi_i`;
- target-curve coefficient Frobenius orbit;
- separability and non-base-field definition of the three quotients.

## Positive Controls

- `phi_i` has degree two and kernel exactly `{O,R_i}`;
- `Psi_i` kills only its designated H018 line among the three candidates;
- the three target curves are Frobenius conjugate.

## Falsification Criterion

Any mismatched kernel, degree, or Frobenius coefficient destroys this explicit
product quotient construction.  A passing construction does not establish a
principal polarization, theta divisor, semilinear descent, scalar H018
sections, or an ECDLP algorithm.

## Reproduction

```bash
sage -python tools/poincare_cubic_product_quotient_n611q.py \
  --out notes/poincare_cubic_product_quotient_n611q.json
sage -python tools/verify_poincare_cubic_product_quotient_n611q.py \
  --primary notes/poincare_cubic_product_quotient_n611q.json \
  --out notes/poincare_cubic_product_quotient_n611q_independent_verifier.json
```
