# Experiment Contract: N611M Level-Curve Canonical Elliptic Returns

## Hypothesis

The complete H018 level curves can support a normalization/Jacobian packet whose
eventual target return does not collapse through either canonical elliptic
projection or their fixed linear span.

## Model

For each declared pencil level `t`, use N608U's complete open correspondence

```text
C_t^open = { (P,Q) : lambda(P,Q)=t }.
```

The two canonical maps are `p1(P,Q)=P` and `p2(P,Q)=Q`.  A divisor packet that
returns through their fixed elliptic span has finite-field label
`[a]P+[b]Q`; projective pairs `(a:b)` enumerate the distinct nonzero spans.
This is only a necessary screen for canonical returns, not a theorem about all
maps from a future normalization or Jacobian.

## Parameters

- registered H018 `E/F_103` fixture;
- complete N608U source levels `0, 1, 17`;
- all 104 projective pairs `(a:b)` over `F_103`;
- exact N608U degree-ten inverse for every nonzero rational `Q`.

## Metrics

- complete source count and elimination degree;
- support of `p1` and `p2` returns;
- minimum support in the canonical return span;
- number of span directions below the quarter-order threshold;
- whether the two projections appear in the projective span panel.

## Positive Controls

- N608U source counts remain `84,108,144`;
- every inverse polynomial has degree ten and no exceptional branch;
- the `(1:0)` and `(0:1)` panel entries equal direct `p1` and `p2` supports.

## Success Criterion

At least one canonical return direction has support at most `#E/4`.  Such a
signal is still not an ECDLP result: it would require a normalized curve,
explicit divisor classes, source return, target descent, rank, and charged
cost.

## Falsification Criterion

If every canonical return direction has support above `#E/4`, the
normalization/Jacobian route may not use a target recovery that factors through
this fixed projection span as its compression mechanism.  Other Jacobian
factors, non-fixed descent, and other curve families remain open.

## Reproduction

```bash
sage -python tools/poincare_level_curve_canonical_returns_n611m.py \
  --out notes/poincare_level_curve_canonical_returns_n611m.json
sage -python tools/verify_poincare_level_curve_canonical_returns_n611m.py \
  --primary notes/poincare_level_curve_canonical_returns_n611m.json \
  --out notes/poincare_level_curve_canonical_returns_n611m_independent_verifier.json
```
