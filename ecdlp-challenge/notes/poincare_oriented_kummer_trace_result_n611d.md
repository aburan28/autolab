# N611D Result: H018 Oriented Kummer Trace Packets

## Result

`NEGATIVE END-TO-END ORIENTED-KUMMER TRACE-PACKET RESULT / MODEL-BOUND /
TOY-EVIDENCE / NONTRACE_PACKET_ALGEBRA_OPEN / NO_ECDLP_CLAIM`

N611D converts every N611C base-field trace factor `T` to an oriented Kummer
label. The label is `x(T)`, and its signed canonical representative is the
member of `{T,-T}` with smaller scalar against the declared generator. The
sign is reconstructed from the exact lifted trace, so it appears as a public
matrix coefficient rather than another unknown.

Every one of the 324 degree-nine H018 residual fibers satisfies the exact
signed relation

```text
sum_f multiplicity(f) * epsilon(T_f) * C[x(T_f)] = -4Q.
```

The 321 collection rows exclude all `Q=37G` rows. Their 54-column matrix has
rank 54, and its factor-log solution recovers the held-out scalar `37` from a
fresh level-zero fiber. The independent verifier recomputes all roots, traces,
inverse-orbit pairings, rank, and descent without importing the producer.

## Charged Boundary

| Metric | Direct trace | Oriented Kummer |
|---|---:|---:|
| Factor labels | 108 | 54 |
| Compression | 1.0x | 2.0x |
| Independent-row lower bound | 108 | 54 |
| Relation rank | 108 | 54 |
| Matrix state for 321 rows | 34,668 | 17,334 |
| Extension-root evaluations | 2,914 | 2,914 |
| Rho group-addition estimate | 13.08 | 13.08 |

The quotient exactly removes the unavoidable inverse-pair duplication, but it
does not change the collection source or the extension arithmetic. Its 54
independent factor logs still require at least 54 independent relation rows,
over four times the registered rho estimate before the remaining costs.

## Interpretation

This is a successful exact Kummer representation test and a negative ECDLP
cost result. It establishes an involution-only constant-factor compression for
this direct trace factor base. It does not establish a generic or
cryptographic-scale speedup, and it does not rule out a packet construction
whose variables are not individually base-group trace logs.

## Red-Team Checks

- Every nonzero traced point occurs with its inverse, and each `x` label has
  exactly that two-point orbit.
- The sign reconstructs the original trace point before any linear solve.
- Held-out target rows are excluded from the collection matrix.
- The admission gates reject promotion: support is not below one quarter of the
  group, rank is not below rho, and compression is not greater than two.

## Three Follow-On Theories

1. Conservative: test whether a nontrivial public endomorphism orbit can give
   a larger signed packet quotient on an explicitly special curve, while
   separately charging the curve restriction and rejecting any conclusion for
   generic prime-field ECDLP.
2. Representation change: construct a compact normalized-Poincare/Jacobian
   packet whose relation variables remain geometric divisors rather than their
   individual traces, with an explicit target-return and source-inverse map.
3. High risk: search for a non-graph H018 correspondence whose projection
   fibers yield a bounded relation kernel and a nonlinear packet-log system
   that avoids one independent original-group log per orbit.

## Next Concrete Action

Specify the target-return and source-inverse maps for a non-trace divisor
packet before collecting another relation matrix. The first admission test must
show super-involution compression and a rank lower bound below rho on a
matched toy family.

## Artifacts

- Contract: `notes/experiment_contract_poincare_oriented_kummer_trace_n611d.md`
- Producer: `tools/poincare_oriented_kummer_trace_probe_n611d.py`
- Receipt: `notes/poincare_oriented_kummer_trace_probe_n611d.json`
- Independent verifier: `tools/verify_poincare_oriented_kummer_trace_n611d.py`
- Verifier receipt:
  `notes/poincare_oriented_kummer_trace_probe_n611d_independent_verifier.json`

## Reproduction

```bash
sage -python tools/poincare_oriented_kummer_trace_probe_n611d.py \
  --out notes/poincare_oriented_kummer_trace_probe_n611d.json
sage -python tools/verify_poincare_oriented_kummer_trace_n611d.py \
  --primary notes/poincare_oriented_kummer_trace_probe_n611d.json \
  --out notes/poincare_oriented_kummer_trace_probe_n611d_independent_verifier.json
```

