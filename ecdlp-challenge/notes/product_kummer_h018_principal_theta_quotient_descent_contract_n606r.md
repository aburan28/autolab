# Experiment Contract: N606R Principal-Theta Quotient Descent Gate

## Hypothesis

The intrinsic type-(1,2) H018 polarization may yield computable theta sections
by descending it through a maximal isotropic order-two subgroup of its kernel
to a principally polarized quotient.  If a subgroup is defined over `F_103`,
the quotient supplies a base-field principal theta divisor; otherwise the
three conjugate quotients define a concrete cubic-orbit descent construction.

## Null Hypothesis

No maximal isotropic order-two subgroup of the H018 kernel is Frobenius stable.
Then a single principal-polarized quotient cannot provide an `F_103` theta
section directly.  The intrinsic route remains open only through an explicitly
constructed cubic-orbit descent, not through a base-field quotient shortcut.

## Assumptions

- The N602A H018 mod-two kernel and its nondegenerate alternating theta
  pairing are correct.
- Standard polarization descent applies in characteristic `103 != 2`: an
  ample line bundle descends through an isotropic kernel subgroup, and the
  quotient polarization is principal for a maximal isotropic subgroup.
- Literature source: Damien Lubicz and Damien Robert, *Computing separable
  isogenies in quasi-optimal time*, LMS J. Comput. Math. 18 (2015),
  [Section 1 / descent discussion](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/CB4581162DE6B7CDC53DC6BFC26F42A4/S146115701400045Xa.pdf/computing_separable_isogenies_in_quasioptimal_time.pdf).

## Parameters

- base field: `F_103`;
- H018 kernel: `K(H018) ~= (Z/2)^2`;
- kernel basis from N602A:
  `g=(pi,1)`, `h=(1+pi,pi)` over `F_4`;
- Frobenius acts by multiplication by `pi` on each coordinate;
- candidate quotient subgroups: the three nonzero order-two lines.

## Metrics

- exact Frobenius orbit of the three order-two lines;
- enumeration of all Frobenius-stable subgroups of `K(H018)`;
- field degree required to stabilize each maximal isotropic line;
- principal-quotient eligibility under the stated descent assumption.

## Positive Control

Replace the order-three Frobenius action by the identity action on the same
abstract kernel.  Every order-two line must then be base-field stable.

## Negative Controls

- The zero subgroup and the whole kernel are Frobenius stable but are not
  maximal isotropic lines.
- The actual order-three Frobenius action must leave no nonzero line stable.

## Success Criterion

If a base-field-stable line exists, construct its quotient and pull back its
principal theta divisor as a named H018 section.  Otherwise preserve the
cubic-orbit quotient data and require explicit `F_103^3` quotient equations,
three conjugate theta divisors, and a semilinear descent calculation before
claiming scalar sections.

## Falsification Criterion

If Frobenius cyclically permutes all three nonzero lines, reject only the
single base-field principal-quotient shortcut.  Do not treat the class-level
`h0=2` budget, a finite theta module, or a quotient over an extension as an
already-evaluated `F_103` section basis.

## Reproduction

```bash
python3 tools/product_kummer_h018_principal_theta_quotient_descent_n606r.py \
  --output notes/product_kummer_h018_principal_theta_quotient_descent_n606r.json
python3 tools/verify_product_kummer_h018_principal_theta_quotient_descent_n606r.py \
  --primary notes/product_kummer_h018_principal_theta_quotient_descent_n606r.json \
  --output notes/product_kummer_h018_principal_theta_quotient_descent_n606r_independent_verifier.json
```
