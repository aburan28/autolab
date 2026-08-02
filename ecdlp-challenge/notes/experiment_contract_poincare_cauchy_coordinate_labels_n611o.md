# Experiment Contract: N611O Cauchy-Coordinate Labels

## Hypothesis

On complete H018 pencil levels, a rational Cauchy coordinate associated with
the moving support `S=-4Q` has a correction alphabet substantially smaller
than the base group, unlike the product-linear returns screened by N608V and
N611M.

## Candidate Labels

For a complete open source `(P,Q)`, let `S=-4Q` and

```text
d_x = x(P)-x(S),        d_y = y(P)+y(S),
m = d_y/d_x,
w = c_8(Q) m,
```

where `c_8(Q)` is the Cauchy-basis coefficient in the N608U evaluator.  The
screen tests `d_x`, `m`, `m^2`, and `w`.  Open-source exclusions ensure
`d_x != 0`; a zero `w` is a valid label.

## Parameters

- registered H018 `E/F_103` fixture;
- N608U complete source levels `0,1,17`;
- all complete sources and no source thinning;
- linear/Kummer control `x(P-4Q)` when finite;
- deterministic uniform `F_103` source-label control.

## Metrics

- defined-label count and support;
- Q-fibres with more than one label;
- support ratio to the matched deterministic uniform-label control;
- best candidate support against `floor(#E/4)=27`.

## Positive Controls

- N608U source counts remain `84,108,144`;
- every Cauchy denominator is nonzero on the declared open sources;
- the linear/Kummer control is reported separately, not used as a candidate;
- deterministic uniform labels use every retained source once.

## Success Criterion

At every level, one fixed candidate has complete defined labels, is not
Q-only, has support at most `27`, and has at most half the uniform-control
support.  Passing is only a factor-label signal.  Target relation, rank,
descent, and charged cost remain mandatory.

## Falsification Criterion

If no candidate meets the joint threshold, reject this bounded Cauchy
coordinate family on the registered fixture.  Other rational maps, global
normalization functions, covers, and relation laws remain open.

## Reproduction

```bash
sage -python tools/poincare_cauchy_coordinate_labels_n611o.py \
  --out notes/poincare_cauchy_coordinate_labels_n611o.json
sage -python tools/verify_poincare_cauchy_coordinate_labels_n611o.py \
  --primary notes/poincare_cauchy_coordinate_labels_n611o.json \
  --out notes/poincare_cauchy_coordinate_labels_n611o_independent_verifier.json
```
