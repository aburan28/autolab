# Experiment Contract: N610U generic P-origin H018 transition

## Hypothesis

The P-origin frame `p^8 D` remains regular on the generic simultaneous
P/Q-origin direction and joins the exact crossing frame by its literal factor
`p-u`.

## Parameters

- coefficient field: `F_(103^6)(c)`
- P parameter: `p=c*u`
- Q/support parameter: `u`
- formal precision: `32`
- levels: `0`, `1`, and `17`

## Metrics

- P-origin valuation of `p^8 D`
- crossing valuation of `(p-u)p^8 D`
- literal transition identity

## Positive control

N609L gives the rational-Q P-origin local representative, and N610J/N610N
give the exact crossing factor.

## Negative control

The unnormalized evaluator must retain its P-origin pole; a nonnegative
valuation after omitting `p^8` would signal that the test is measuring the
wrong frame.

## Success criterion

All selected levels have regular `p^8 D`, crossing weight one, and an exact
identity in the declared truncated formal ring.

## Falsification criterion

A pole in `p^8 D` or failed identity rejects this generic overlap convention
and leaves the inverse-support transition unresolved.

## Reproduction command

```bash
sage -python tools/poincare_p_origin_generic_transition_n610u.py \
  --out notes/poincare_p_origin_generic_transition_n610u.json
```
