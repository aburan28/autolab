# Experiment Contract: N611D H018 Oriented Kummer Trace Packets

## Hypothesis

The exact N611C Frobenius-trace factors can be compressed through the Kummer
quotient `E/{+/-1}`: use one label for each nonzero `x(T)` and attach the
orientation sign determined by a canonical representative of `{T,-T}`. This
could lower the factor-base rank while preserving an exact target relation and
held-out descent.

## Null hypothesis

The quotient only yields the involution's constant-factor compression, so its
independent label count and rank lower bound remain above rho even before
extension arithmetic, factorization, or linear algebra are charged.

## Parameters

- `E/F_103: y^2=x^3+x+24`, prime group order `109`;
- H018 degree-nine residual fibers at levels `0,1,17`, for all nonzero `Q`;
- every residual factor is lifted and traced exactly as in N611C;
- Kummer label: `x(T)` for nonzero trace point `T`;
- canonical orientation: the member of `{T,-T}` with smaller generator scalar;
- held-out target: `37G`, excluded from collection rows.

## Exact relation

For each traced factor point `T`, write `T = epsilon(T) C[x(T)]`, where
`C[x(T)]` is the canonical orientation and `epsilon(T)` is `+1` or `-1`.
Then each residual fiber must satisfy

```text
sum_f multiplicity(f) * epsilon(T_f) * C[x(T_f)] = -4Q.
```

The unknown for a Kummer label is the generator scalar of its canonical point;
the signs are public, reconstructed from the exact lifted trace, and belong in
the matrix coefficients rather than as extra unknowns.

## Metrics

- exact root lifting, trace descent, orbit pairing, and orientation replay;
- Kummer-label support, compression ratio, relation rank, and matrix state;
- held-out factor-log solve and fresh target descent;
- independent-row lower bound compared with rho.

## Positive control

All 324 fibers must replay exactly. The Kummer relation matrix, trained after
excluding every `Q=37G` row, must solve the held-out level-zero fiber.

## Negative controls

- Every nonzero label must have exactly the two-point orbit `{T,-T}`.
- The chosen sign must reconstruct the full base-field trace point.
- The held-out target is absent from collection rows.
- A quotient result is not an ECDLP improvement unless support and its
  rank-imposed minimum collection both beat rho after all costs are charged.

## Admission criterion

Promote only if all exact controls pass, support is at most one quarter of the
group order, the minimum independent rows is below the rho estimate, and the
compression is materially greater than the involution-only factor two. Any
failure is a scoped result about this direct oriented Kummer trace route.

## Reproduction command

```bash
sage -python tools/poincare_oriented_kummer_trace_probe_n611d.py \
  --out notes/poincare_oriented_kummer_trace_probe_n611d.json
```

## Claim boundary

`HYPOTHESIS / H018_ORIENTED_KUMMER_TRACE_PACKETS / MODEL-BOUND /
TOY-EVIDENCE / NO_ECDLP_CLAIM`.

