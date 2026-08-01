# Experiment Contract: N606S Theta-Eigenline Field Gate

## Hypothesis

The three `F_103^3` principal-quotient candidates from N606R may already have
scalar theta eigendirections over the cubic field.  If so, cubic quotient
equations could suffice for named theta sections.

## Null Hypothesis

Each quotient-line lift has irreducible characteristic polynomial over
`F_103`.  Then its scalar theta eigendirections require the quadratic phase
field, and actual quotient-plus-theta data must be handled over `F_103^6`.

## Parameters

- quaternionic H018 lifts `I,J,K` and Frobenius descent matrix `D` from N602B;
- `I^2=J^2=K^2=-1`, `D^3=1`, and `D` cycles `I,J,K` by conjugation;
- base field `F_103`, cubic kernel field `F_103^3`.

## Metrics

- characteristic polynomials and base-field roots of `I,J,K`;
- explicit `+i` and `-i` eigenlines over `F_103(i)`;
- exact transport of the `+i` eigenspace under `D`;
- field-degree compositum requirement.

## Controls

- `D^3=1` and all three conjugacies must replay exactly.
- A base-field root of `T^2+1` would falsify the claimed quadratic phase need.

## Reproduction

```bash
sage -python tools/product_kummer_h018_theta_eigenline_field_n606s.py \
  --output notes/product_kummer_h018_theta_eigenline_field_n606s.json
sage -python tools/verify_product_kummer_h018_theta_eigenline_field_n606s.py \
  --primary notes/product_kummer_h018_theta_eigenline_field_n606s.json \
  --output notes/product_kummer_h018_theta_eigenline_field_n606s_independent_verifier.json
```
