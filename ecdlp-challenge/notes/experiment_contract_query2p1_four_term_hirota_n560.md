# Experiment Contract: N560 four-term quadratic Fay/Hirota preflight

## Hypothesis

The Query2P1 raw product admits a centered four-term bilinear identity with
coefficients in the affine monomial basis of `L(12O)`, beyond N550's
three-term, pole-at-most-six family.

## Null hypothesis

Even after post-fitting this larger coefficient space, no identity survives
independent source decks on the same curve.

## Parameters

- prime-order curves over `F_251`, `F_257`, `F_269`
- deck widths `2`, `3`
- centered terms at shifts `4`, `3`, `2`, `0`
- coefficient space `L(12O)`
- affine-template, coordinate-hash, and three shuffled decks

## Success criterion

An all-term null vector must survive the concatenated independent-deck matrix,
while direct source-root and constant-product controls pass and shuffled data
rejects it.

## Reproduction command

```bash
python3 -B tools/query2p1_four_term_hirota_probe_n560.py \
  --out notes/query2p1_four_term_hirota_probe_n560.json
```
