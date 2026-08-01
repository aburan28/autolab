# Experiment Contract: N609D formal Q=O boundary pencil

## Hypothesis

Within the N608R regularized model, every pencil member `x-t*1` has a uniform
degree-nine `Q=O` boundary divisor: finite pole degree eight plus one zero at
`P=O` as a section of `L(9O)`.

## Null hypothesis

Changing `t` cancels the leading `x^4` term, changes the finite pole degree,
or makes the selected regularized pencil members fail base-field descent.

## Parameters

- Registered H018 fixture over `F_103`.
- N608R selected regular cover sections `1` and `x`.
- All 103 values of `t` in `F_103`, with declared open-level controls
  `t=0,1,17`.
- `L(9O)` basis `1,x,y,x^2,xy,x^3,x^2y,x^4,x^3y`.

## Metrics

- boundary coefficient vectors and their base-field descent;
- finite pole degree and inferred origin-zero multiplicity;
- rational finite zero count by level;
- histogram of rational finite zero counts across the full pencil.

## Controls

- Positive: the regularized `1` section has boundary vector
  `[1,0,0,0,0,0,0,0,0]`.
- Uniform-leading-term: the `x^4` coefficient is nonzero and independent of
  `t`.
- Negative: the `y` section remains nonregular at the same source origin and
  cannot define an additional boundary pencil direction.

## Success criterion

Every `t` produces a base-field boundary coefficient vector with finite pole
degree eight and exactly one origin zero as an `L(9O)` section.

## Falsification criterion

Any `t` cancels the leading term, lowers the finite pole degree, or fails
descent. Such a result would narrow the claimed uniform boundary pencil.

## Boundary

This measures only formal `Q=O` fibres under N608R. It neither identifies the
N609A residual compactification nor establishes a global divisor class,
second projection, normalization, factor base, relations, rank, target
descent, cost advantage, or ECDLP speedup.

## Reproduction

```bash
sage -python tools/poincare_origin_pencil_n609d.py --out notes/poincare_origin_pencil_n609d.json
sage -python tools/verify_poincare_origin_pencil_n609d.py --primary notes/poincare_origin_pencil_n609d.json --out notes/poincare_origin_pencil_n609d_structural_replay.json
```
