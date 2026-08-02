# Experiment Contract: N610V exceptional P-origin specialization

## Hypothesis

After cancelling the P-origin pole with `p^8`, the generic simultaneous-frame
coefficients extend to the exceptional direction `c=0` through the declared
formal order, and the crossing factor specializes to `-u`.

## Parameters

- coefficient field before specialization: `F_(103^6)(c)`
- P parameter before specialization: `p=c*u`
- exceptional specialization: `c=0`
- inspected coefficient orders: `0..11`
- formal precision: `32`
- levels: `0`, `1`, and `17`

## Success criterion

At each level, every inspected coefficient of `p^8 D` has denominator
nonzero at `c=0`, its specialized leading coefficient is nonzero, and the
specialized crossing frame equals `-u` times the specialized P-origin frame.

## Falsification criterion

A denominator vanishing at `c=0`, a zero leading term, or a failed specialized
transition rejects this exceptional-direction extension at the inspected order.

## Reproduction command

```bash
sage -python tools/poincare_p_origin_c_zero_specialization_n610v.py \
  --out notes/poincare_p_origin_c_zero_specialization_n610v.json
```
