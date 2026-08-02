# Experiment Contract: N610F H018 Q-only target relation obstruction

## Hypothesis

At least one certified complete H018 level has no relation of the form
`P + G(Q) = H(t)`, even when `G` is an arbitrary group-valued function of `Q`.

## Null hypothesis

Every fixed-`Q` fibre at each certified level has at most one `P`, leaving the
Q-only correction family unblocked by fibre multiplicity.

## Parameters

- field/curve family: fixed H018 fixture on `E/F_103: y^2=x^3+x+24`
- levels: complete certified levels `0`, `1`, and `17`
- factor base: not applicable; target-relation preflight
- relation shape: `P + G(Q) = H(t)` with unrestricted Q-only `G`
- baseline: N609N additive and N609O coefficient-ratio corrections

## Metrics

- source count per level
- Q-fibre count and maximum P-fibre size
- number of Q fibres containing distinct P values
- explicit witness fibre

## Positive control

The complete-source generator must reproduce the known source counts.

## Negative control

A witness fibre must contain two unequal P values at the same Q.

## Success criterion

Every certified level has a repeated-Q, distinct-P witness.  Such a witness
proves the declared Q-only relation family cannot be single-valued at that
level.

## Falsification criterion

Any certified level has no repeated-Q distinct-P fibre.

## Reproduction command

```bash
sage -python tools/poincare_q_only_relation_preflight_n610f.py \
  --out notes/poincare_q_only_relation_n610f.json
```
