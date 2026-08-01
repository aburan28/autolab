# Experiment Contract: N609I conditional global H018 class of the residual

## Hypothesis

Under the declared N608O/N608R normalized-H018 model, the N609H compact
P-fibres are restrictions of the selected global section of

```text
L_H = p1^*O(9O) tensor p2^*O(11O) tensor (id,z)^*P,
z = -3-Frob.
```

Consequently the residual divisor has P-fibre degree 9, Q-fibre degree 11,
self-intersection 4, and arithmetic genus 3 if it is an effective Cartier
divisor on the product surface.

## Null hypothesis

The H018 matrix or its P-fibre restriction fails to match N609H's exact pole
divisor, so the conditional global-class identification is invalid.

## Parameters

- `E/F_103: y^2=x^3+x+24`.
- CM relation `F^2+5F+103=0`.
- H018 form `[[9,2+F],[-3-F,11]]`.
- N608R selected section space and N609H compact P-fibre receipt.

## Metrics

- exact H018 Hermitian matrix and determinant;
- base-field restriction `z=[-4]`;
- P-fibre class `8O+[-4Q]` agreement;
- Q-fibre degree, self-intersection, and conditional arithmetic genus.

## Controls

- Swapped Poincare orientation must not equal H018.
- The N609H primary and replay receipts must pass and retain all 108 P-fibres.
- N608R must retain selected cover space `span{1,x}`.

## Success criterion

All formal-class and bound-receipt checks pass, yielding the stated
model-conditional degrees and arithmetic genus.

## Falsification criterion

Any orientation, restriction, or receipt mismatch leaves the global class
unidentified and blocks Q-projection/genus use.

## Boundary

This identifies a class only under the stated normalized global model. It does
not give an explicit global equation, prove the compact divisor is Cartier or
smooth, analyze normalization, construct a Jacobian, define relations, or
provide rank, descent, cost, or an ECDLP improvement.

## Reproduction

```bash
sage -python tools/poincare_global_class_n609i.py --out notes/poincare_global_class_n609i.json
sage -python tools/verify_poincare_global_class_n609i.py --primary notes/poincare_global_class_n609i.json --out notes/poincare_global_class_n609i_structural_replay.json
```
