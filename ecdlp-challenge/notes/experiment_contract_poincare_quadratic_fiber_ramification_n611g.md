# Experiment Contract: N611G H018 Quadratic-Subfield Residual Audit

## Hypothesis

The H018 residual multiplicity signal from N611F is not an artifact of
restricting Q to `E(F_103)). Reconstructing the selected coefficient vector
from actual degree-nine cover lifts of deterministic points in
`E(F_(103^2)) \\ E(F_103)` may expose additional geometric branch candidates
or confirm square-free behavior on a genuinely non-rational sample.

## Null hypothesis

The cover lift, coefficient reconstruction, or residual identity fails outside
the base field; or the finite extension sample has no interpretable relation
to global branch geometry. Neither outcome supports a normalization claim.

## Parameters

- H018 target curve `E/F_103: y^2=x^3+x+24`;
- its declared degree-nine cover, deck generator, and cover-to-target map;
- 16 deterministic non-base points selected from the quadratic subfield of
  `F_(103^6)`;
- levels `0,1,17`;
- residual after exact removal of the `P=4Q` graph factor.

## Metrics

- quadratic-subfield and non-base membership of Q;
- cover-lift and nine-point deck-orbit checks;
- exact coefficient interpolation against cover-side x values;
- residual degree, discriminant-zero status, derivative-gcd degree, and
  repeated-factor profile.

## Positive controls

- `pi(R)=Q`, all nine `R+iT` map to Q, and their image-evaluation matrix
  has rank nine;
- each reconstructed coefficient vector lies in the quadratic subfield;
- the residual has degree nine after graph removal.

## Negative controls

- A base-field Q is rejected from the extension sample.
- Mutating one deck x value breaks the interpolation control.
- This finite sample is not an all-geometric branch divisor, genus, or
  ECDLP result.

## Admission criterion

Advance only if all lift/interpolation/residual controls pass. Any observed
multiplicity is a candidate for later branch-divisor construction; no
normalization, Jacobian, relation, rank, descent, or cost claim is allowed.

## Reproduction command

```bash
sage -python tools/poincare_quadratic_fiber_ramification_n611g.py \
  --out notes/poincare_quadratic_fiber_ramification_n611g.json
```

## Claim boundary

`HYPOTHESIS / H018_QUADRATIC_SUBFIELD_RESIDUAL_AUDIT / MODEL-BOUND /
TOY-EVIDENCE / NO_ECDLP_CLAIM`.

