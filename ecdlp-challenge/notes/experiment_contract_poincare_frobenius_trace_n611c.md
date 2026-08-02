# Experiment Contract: N611C H018 Frobenius-Trace Factor Relations

## Hypothesis

Factoring the degree-nine H018 residual fibers into irreducible x-polynomials
produces Frobenius-orbit divisor packets whose trace points form a compressed
factor base on `E(F_103)`. The exact divisor trace can then yield a complete
relation, rank, and individual-descent route unlike pointwise `P+cQ` laws.

## Null hypothesis

The trace packets occupy a base-field-scale factor base or require at least
linear-in-group-order fiber work, so exact relation/rank/descent does not
beat Pollard rho even on the registered toy curve.

## Parameters

- `E/F_103: y^2=x^3+x+24`, prime group order `109`;
- H018 residual fibers from N609A at levels `0,1,17` and all 108 nonzero Q;
- residual x-polynomial degree nine after removing the declared extraneous
  graph factor;
- one factor variable for each nonzero base-field trace point of an
  irreducible residual factor;
- held-out target: `37G`, excluded from collection rows.

## Exact relation

For a residual factor `f` of multiplicity `m`, lift one full root orbit to
`F_(103^deg(f))`, form its point trace `T_f`, and use

```text
sum_f m*T_f = -4Q.
```

This is the group-law image of the degree-nine P-zero divisor with pole
divisor `8[O]+[-4Q]` certified by N609H.

## Metrics

- factor degrees, lifted root/curve checks, and Frobenius descent of traces;
- exact fiber trace replay;
- trace-factor-base support and relation-matrix rank;
- held-out factor-log solve and fresh target descent;
- fiber-factor calls, the rank-imposed minimum independent relation rows,
  extension-root evaluations, matrix state, and rho group-addition estimate.

## Positive control

Every full residual fiber must reproduce `-4Q`; after excluding `37G` from
collection, the solved trace logs must recover its scalar through a fresh
level fiber.

## Negative controls

- The extraneous graph factor is removed exactly once before factorization.
- Every trace must be Frobenius fixed and lie on the base curve.
- Rows involving the held-out target are excluded from the factor-log system.
- A successful toy solve is not an ECDLP improvement unless the complete
  charged source work and factor-base support beat rho under a common model.

## Admission criterion

Promote only if exact controls pass, trace support is below one quarter of the
base group, and even the rank-imposed minimum number of independent relation
fibers is below the toy rho estimate before charging extension arithmetic or
linear algebra. A result failing either cost or compression is a scoped
negative result.

## Reproduction command

```bash
sage -python tools/poincare_frobenius_trace_probe_n611c.py \
  --out notes/poincare_frobenius_trace_probe_n611c.json
```

## Claim boundary

`HYPOTHESIS / H018_FROBENIUS_TRACE_FACTOR_RELATIONS / MODEL-BOUND /
TOY-EVIDENCE / NO_ECDLP_CLAIM`.
