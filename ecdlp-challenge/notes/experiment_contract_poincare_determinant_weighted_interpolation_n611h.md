# Experiment Contract: N611H Determinant-Weighted H018 Coefficient Interpolation

## Hypothesis

The deck-invariant determinant factor det(V_Q) from N609J cancels the
high-pole gauge of the raw moving Cauchy coefficients. The weighted
coefficients det(V_Q)*c_i(Q) may descend to low-pole functions of Q, giving
a finite normalized equation candidate.

## Null hypothesis

The determinant-weighted coefficients fail base-field descent or retain
near-full pole order. Then this scalar determinant normalization does not
supply a global coefficient law.

## Parameters

- All 108 nonzero Q values on the H018 F_103 fixture;
- N608V coefficient vectors and N609J moving-basis determinant;
- Riemann--Roch interpolation spaces L(dO) through d=109;
- low-pole admission threshold d<=24.

## Controls

- Every det(V_Q) must be nonzero and base-field fixed.
- A single deck shift preserves the determinant.
- The raw N608X first-fit orders remain the negative baseline.
- A low-pole fit is only a normalized-chart observation, not a Cartier,
  branch, normalization, or ECDLP result.

## Success criterion

All weighted coefficients descend and at least one nonconstant coefficient
fits below the raw 106 pole barrier; promote to global-equation work only if
the full selected coefficient vector fits through L(24O).

## Reproduction

```bash
sage -python tools/poincare_determinant_weighted_interpolation_n611h.py \
  --out notes/poincare_determinant_weighted_interpolation_n611h.json
```

