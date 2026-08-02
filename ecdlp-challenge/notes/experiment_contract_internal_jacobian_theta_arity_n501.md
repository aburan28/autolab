# Experiment Contract: N501 higher-arity internal theta occupancy

## Hypothesis

Four or five signed theta factors can overcome the split-cover Prym occupancy
loss observed for three factors, yielding rank-bearing relations and blind
target descents above matched random odd controls.

## Parameters

- covers: existing `F_59`, `F_83`, and `F_103` rows with `c=1,2,3`
- signed theta base: eight target-independent split-fiber points
- arities: `3`, `4`, and `5`
- relation and descent targets: four each
- matched random odd controls: four per row

## Metrics

- signed witness count and materialized factor/sign assignments
- factor-log and augmented rank
- blind descent targets with witnesses
- matched-control relation and descent lifts
- wrong-target rejection and true-log consistency

## Positive control

The attack-invalid zero-Prym selector must yield relation witnesses while
retaining the same E-component and signed-factor representation.

## Success criterion

Promotion requires an at-least-eightfold relation and descent lift across the
field sweep, increasing rank, correct target rows, and a charged collector
below rho and BSGS. This preflight does not claim the final cost gate.

## Falsification criterion

Control-like relation/descent supply, fixed low rank, or tuple growth that
dominates the density gain rejects this split higher-arity realization.

## Reproduction command

```bash
HOME=/private/tmp sage -python tools/internal_jacobian_theta_arity_probe_n501.py \
  --out notes/internal_jacobian_theta_arity_probe_n501.json
```
