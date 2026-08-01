# Experiment Contract: N608K rational H018-to-cover intertwiner

## Hypothesis

The N608J degree-eleven evaluation matrices arise from a regular morphism on a
nonempty declared open locus, not merely unrelated fibre maps.  On the cover
open set where `phi(R)` avoids the source divisor support and `R` avoids the
cover origin, the canonical section `1` of `O_Eprime(2Oprime)` supplies the
target weight.

## Null hypothesis

The evaluation construction has no regular open locus, or its determinant is
identically zero.  Then N608J's finite full ranks do not even define a rational
intertwiner in the stated divisor model.

## Parameters

- N608I abstract target `F ~= pi_*O_Eprime(2Oprime)`;
- N608J degree-nine `pi`, degree-eleven `phi`, and `z=-3-pi`;
- source divisor model `O(8O+[zQ])` from N606Y;
- cover points `R` with `Q=pi(R)` and evaluation point `phi(R)`.

## Declared Open Locus

Remove from `Eprime` the finite conditions

- `R=Oprime`, where the target section `1` does not trivialize `O(2Oprime)`;
- `phi(R)=O`, where the affine source evaluator has a pole;
- `phi(R)=z*pi(R)` or `phi(R)=-z*pi(R)`, where the named source basis has a
  possible support pole.

The test verifies that this locus contains every registered rational fibre and
that one nonzero determinant proves generic full rank.

## Positive and Negative Controls

- Positive: degree-eleven `phi` separates the cyclic degree-nine deck group.
- Negative: `pi` collapses the same deck group and has rank one.

## Success Criterion

The source and target weights are algebraically named on the open locus, all
registered rational fibres lie in it, and the determinant is nonzero at least
once (here, on every registered rational fibre).

## Falsification Criterion

An empty regular locus, rank failure, or failure of the rank-one collapse
control rejects the rational-intertwiner claim.

## Reproduction

```bash
sage -python tools/poincare_rational_intertwiner_n608k.py \
  --out notes/poincare_rational_intertwiner_n608k.json
sage -python tools/verify_poincare_rational_intertwiner_n608k.py \
  --primary notes/poincare_rational_intertwiner_n608k.json \
  --out notes/poincare_rational_intertwiner_n608k_independent_verifier.json
```

## Claim Boundary

`OBSERVATION / EXPLICIT RATIONAL LOCAL MORPHISM / MODEL-BOUND / TOY-EVIDENCE /
NORMALIZED-POINCARE CHART GLUING OPEN / NO SECTION EVALUATOR / NO_ECDLP_CLAIM`.
