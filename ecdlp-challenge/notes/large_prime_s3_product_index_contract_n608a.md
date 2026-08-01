# Experiment Contract: N608A Large-Prime S3 Product Index

## Hypothesis

The pair-plus-large-prime S3 image has an exact univariate product-tree or
resultant membership polynomial whose degree or coefficient support is much
smaller than its explicit image set.

## Null hypothesis

For this raw x-image representation, the square-free membership polynomial has
one root per distinct image x-value and therefore retains essentially the full
`B^2 S` image scale.

## Parameters

- Generated public prime-order curves at 20 and 22 bits, seed `0x12345678`.
- Low-x primary bases of sizes `8,16`.
- Independent low-tail and random-x secondary buckets of sizes `16,32`.
- Signed primary pairs and signed secondary additions, as in the prior S3
  large-prime collector.

## Metrics

- Raw image-entry count and distinct x-image cardinality.
- Degree and nonzero-coefficient count of the square-free product polynomial.
- Exact membership on constructed positive x-values and sampled nonimage x-values.
- Image cardinality relative to the raw entry count and to `p`.

## Positive control

Every constructed signed pair-plus-large image must be a zero of the product
polynomial.

## Negative control

Sampled base-field x-values outside the enumerated image must be nonzeros.

## Success criterion

Promotion requires a repeatable substantial degree or support compression while
preserving exact membership; otherwise this rejects only the obvious univariate
product-index representation.

## Reproduction

```bash
sage -python tools/large_prime_s3_product_index_n608a.py --out notes/large_prime_s3_product_index_n608a.json
```
