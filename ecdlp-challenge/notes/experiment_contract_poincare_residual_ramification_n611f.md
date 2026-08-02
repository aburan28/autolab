# Experiment Contract: N611F H018 Residual Ramification Audit

## Hypothesis

The selected H018 degree-nine residual correspondence has a sparse, exact
open-fiber ramification profile compatible with the conditional arithmetic
genus-three class from N609I. In particular, repeated residual roots should
agree with the formal P-direction critical sources from N609F rather than
appear at base-field scale.

## Null hypothesis

Residual multiplicities are absent where the formal critical calculation
predicts them, occur at base-field scale, or fail to reconstruct the declared
open source geometry. In any case, no normalization or Jacobian construction
is promoted without a genuine global Cartier curve.

## Parameters

- `E/F_103: y^2=x^3+x+24`, group order `109`;
- H018 residuals from all 108 nonzero `Q) values;
- levels `0,1,17);
- residual polynomial: the degree-nine quotient after removing the certified
  `P=4Q) graph factor.

## Metrics

- residual degree, discriminant-zero count, and derivative-gcd degree;
- irreducible repeated-factor degrees and multiplicities;
- recovery of every repeated base-field root as a valid curve point through
  the declared residual `y=-a(x)/b(x)) formula;
- comparison of repeated-root locations with N609F's formal P-direction
  critical-source counts.

## Positive control

For every residual, the factorization multiplicity excess must equal the
degree of `gcd(f,f')). Every repeated base-field root must satisfy the curve
equation and the selected level equation.

## Negative controls

- The removed graph factor is divided exactly once.
- A result is not interpreted as global ramification unless it includes the
  true pole, `Q=O), non-rational fibers, and a Cartier extension.
- The N609F derivative count is a cross-check, not evidence of a global
  genus theorem.

## Admission criterion

Advance only if exact residual controls pass and the repeated-root profile
matches the formal P-direction profile at every selected level. This remains
a geometry observation, not a factor base, relation, rank, descent, or ECDLP
result.

## Reproduction command

```bash
sage -python tools/poincare_residual_ramification_probe_n611f.py \
  --out notes/poincare_residual_ramification_probe_n611f.json
```

## Claim boundary

`HYPOTHESIS / H018_RESIDUAL_RAMIFICATION_AUDIT / MODEL-BOUND /
TOY-EVIDENCE / NO_ECDLP_CLAIM`.

