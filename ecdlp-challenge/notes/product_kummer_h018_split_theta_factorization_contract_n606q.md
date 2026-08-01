# Experiment Contract: N606Q Split Theta-Factorization Gate

## Hypothesis

The positive Hermitian H018 polarization might secretly be the pullback of the
split type-(1,2) product polarization through a CM-coordinate change. Such a
factorization would give a concrete product-theta starting point for the two
sections required by N606P:

```text
H018 = A^* diag(1,2) A,  A in M_2(Z[pi]).
```

## Null Hypothesis

No such factorization exists. Then a product-theta evaluator cannot be
obtained by a CM coordinate change, and any N606Q section construction must
use intrinsically non-split Poincare/theta data.

## Parameters

- CM order: `Z[pi]`, with `pi^2+5*pi+103=0`;
- norm: `N(a+b*pi)=a^2-5ab+103b^2`;
- target matrix:

```text
H018 = [[9, 2+pi], [-3-pi, 11]].
```

## Metrics

- exact lower bound for the norm of every element with nonzero `pi`
  coefficient;
- all possible diagonal factorizations of `9` and `11` through
  `diag(1,2)`;
- exact off-diagonal comparison;
- a split positive control `diag(1,2)=I^*diag(1,2)I`.

## Success Criterion

If a factorization is found, materialize its product-theta sections and begin
the N606P base-locus, fiber-degree, and source-recovery audits.

## Falsification Criterion

If the diagonal norm budgets force all entries of `A` into `Z`, then the
off-diagonal of `A^*diag(1,2)A` is integral and cannot equal `2+pi`. Record
the restricted obstruction and do not treat a product-theta basis change as
an H018 evaluator.

## Reproduction

```bash
python3 tools/product_kummer_h018_split_theta_factorization_n606q.py \
  --output notes/product_kummer_h018_split_theta_factorization_n606q.json
python3 tools/verify_product_kummer_h018_split_theta_factorization_n606q.py \
  --primary notes/product_kummer_h018_split_theta_factorization_n606q.json \
  --output notes/product_kummer_h018_split_theta_factorization_n606q_independent_verifier.json
```
