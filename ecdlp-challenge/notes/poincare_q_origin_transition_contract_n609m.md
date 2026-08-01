# Experiment Contract: N609M Q-origin frame transition for the determinant

## Hypothesis

N608R's source correction is a literal Cech frame transition for N609J's
determinant: the raw moving Cauchy evaluation and the corrected Q-origin
evaluation agree exactly, and the corrected form is regular at `Q=O` for every
finite rational P at levels `0,1,17`.

## Null hypothesis

The two evaluations differ, the corrected expression retains a Q-origin pole,
or its boundary constant does not match the N609C degree-nine pencil.

## Parameters

- H018 fixture over `F_103`.
- All 108 finite rational P values.
- Levels `t=0,1,17`.
- Formal Q-origin series precision 120.

## Metrics

- raw and corrected evaluator equality;
- valuation of the corrected expression;
- boundary-constant match with the N609C section.

## Positive control

The raw-to-corrected change of basis is an algebraic identity before taking a
boundary constant.

## Negative control

Perturbing one correction term must change at least one finite-P formal
expression.

## Success criterion

All 324 finite `(P,t)` rows agree exactly, are Q-origin regular, and match the
N609C boundary pencil.

## Falsification criterion

Any mismatch or remaining pole blocks treating the correction as a usable
Q-origin frame transition for the determinant.

## Boundary

This is a finite-P Q-origin transition test. It does not cover the simultaneous
P-origin/Q-origin point, non-rational geometry, global cocycle compatibility,
a global Cartier divisor, smoothness, normalization, Jacobians, relations,
rank, descent, costs, or ECDLP improvement.

## Reproduction

```bash
sage -python tools/poincare_q_origin_transition_n609m.py --out notes/poincare_q_origin_transition_n609m.json
```
