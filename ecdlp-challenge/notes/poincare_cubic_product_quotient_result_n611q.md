# Result: N611Q Cubic Product Quotient Models

## Status

`OBSERVATION / EXPLICIT_CUBIC_H018_PRODUCT_QUOTIENT_MODELS / MODEL-BOUND /
TOY-EVIDENCE / PRINCIPAL_POLARIZATION_AND_THETA_OPEN / NO_ECDLP_CLAIM`

## Exact Result

For the cubic two-torsion orbit `R,R^103,R^(103^2)`, N611Q constructs the
three degree-two Velu maps `phi_i:E -> E_i` with kernels `<R_i>`.  Combined
with

```text
T(P,Q) = (P + Frob(Q), Q),
Psi_i(P,Q) = (P + Frob(Q), phi_i(Q)),
```

they give explicit product quotients

```text
Psi_i : E^2 -> E x E_i.
```

`Psi_g`, `Psi_h`, and `Psi_(g+h)` kill exactly `G_g`, `G_h`, and `G_(g+h)`,
respectively, among the three materialized H018 order-two lines.  Each factor
isogeny has exact kernel `{O,R_i}` and degree two.

The first quotient elliptic factor is

```text
E_g: y^2 = x^3 + (42*a^2+16*a+62)x + (30*a^2+47*a+53)
```

over `F_(103^3)`.  The other two target curves are its Frobenius conjugates;
none is base-field defined.

## Evidence

- Producer: five of five degree, kernel, unique-line, conjugacy, and
  non-base-field gates pass.
- Independent replay: separately reconstructs all three Velu maps, product
  kernels, and target coefficient cycle; all four checks pass.

## Interpretation

This is an actual quotient construction, not merely an abstract subgroup
record.  It realizes the abelian-variety quotient required by the cubic theta
route as a product model.  It does **not** yet show that the H018 polarization
descends to a principal polarization on that product, identify its theta
divisor, construct three conjugate pullbacks, solve semilinear descent, or
produce scalar H018 sections.  It supplies no factor base, relation rank,
target descent, cost advantage, or ECDLP speedup.

## Next Concrete Action

Transport the H018 Hermitian/polarization data through `Psi_g` and determine
the descended polarization on `E x E_g`.  The next experiment must prove its
type is principal before choosing a theta divisor or discussing section
descent.

## Reproduction

```bash
sage -python tools/poincare_cubic_product_quotient_n611q.py \
  --out notes/poincare_cubic_product_quotient_n611q.json
sage -python tools/verify_poincare_cubic_product_quotient_n611q.py \
  --primary notes/poincare_cubic_product_quotient_n611q.json \
  --out notes/poincare_cubic_product_quotient_n611q_independent_verifier.json
```
