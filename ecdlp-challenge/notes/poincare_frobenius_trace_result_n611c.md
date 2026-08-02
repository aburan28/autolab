# N611C Result: H018 Frobenius-Trace Factor Relations

## Result

`NEGATIVE END-TO-END TRACE-FACTOR RELATION RESULT / MODEL-BOUND /
TOY-EVIDENCE / NO_ECDLP_CLAIM`

The degree-nine H018 residual divisor supplies an exact relation mechanism.
For every irreducible residual x-factor, lifting its root orbit to the
corresponding extension field and summing its curve points gives a
Frobenius-fixed trace point in `E(F_103)`. Across all 324 fibers at levels
`0,1,17`, every relation replays exactly:

```text
sum_f multiplicity(f) * Trace(f) = -4Q.
```

The relation matrix has full rank `108`. Excluding every row for the held-out
target `37G`, the solved trace logs recover `37` through a fresh level-zero
fiber. The independent verifier recomputes all extension roots, trace sums,
factor-degree counts, rank, and held-out descent.

## Charged Boundary

| Metric | Value |
|---|---:|
| H018 residual fibers replayed | 324 |
| Collection fibers after holding out `37G` | 321 |
| Extension-root evaluations | 2914 |
| Trace-factor-base size | 108 of 109 group points |
| Relation rank / minimum independent rows | 108 / 108 |
| Matrix state | 34,668 field elements |
| Rho group-addition estimate | 13.08 |

Even an optimal collector needs at least 108 independent relation fibers to
solve 108 trace-factor logs, before factorization or extension arithmetic.
That rank lower bound is already over eight times the registered rho estimate,
and the trace-factor base occupies essentially the whole base group.

## Interpretation

This is a complete toy relation/rank/descent receipt, but it is not an ECDLP
speedup. Passing to Frobenius traces collapses each geometric packet back into
an almost full-size base-field factor base. The result rejects this direct
trace-factor route, not a normalization/Jacobian representation that retains
additional non-base-field structure with a proven target-return map.

## Artifacts

- Contract: `notes/experiment_contract_poincare_frobenius_trace_n611c.md`
- Probe: `tools/poincare_frobenius_trace_probe_n611c.py`
- Receipt: `notes/poincare_frobenius_trace_probe_n611c.json`
- Independent verifier: `tools/verify_poincare_frobenius_trace_n611c.py`
- Verifier receipt:
  `notes/poincare_frobenius_trace_probe_n611c_independent_verifier.json`

## Next Concrete Action

Search for a compact normalization/Jacobian factor whose labels do not trace
directly to an almost full `E(F_p)` factor base, and require an explicit
target return, source inverse, rank, descent, and charged rho comparison.
