# N611B Result: Dickson-Fiber Triple Screen

## Result

`NEGATIVE RESULT / DICKSON-FIBER TRIPLE MEMBERSHIP PREFLIGHT / MODEL-BOUND /
TOY-EVIDENCE / NO_ECDLP_CLAIM`

The declared Dickson family does not pass the source-solver admission gate.
For each factor triple with elementary polynomial
`U(T)=T^3-e1*T^2+e2*T-e3`, membership in a fiber is exactly the three
coefficient equations of

```text
D_d(T,a)-c == 0 mod U(T).
```

Only two of the twelve pre-registered fibers had 6 to 48 affine factor
points.  Their generic membership remainder is already larger than the
unordered x-pair state, so it does not supply the bounded-support algebra
required to improve on a pair representation.

| `(d,a,c)` | Factor points | Remainder monomials | x-pair state | Support ratio | Triple count | Control mean / max |
|---|---:|---:|---:|---:|---:|---:|
| `(10,2,0)` | 10 | 66 | 15 | 4.40 | 4 | 0.875 / 4 |
| `(25,1,0)` | 20 | 807 | 55 | 14.67 | 8 | 10.125 / 18 |

The first row has a 4.57x lift over the *mean*, but it merely ties the
strongest one of 16 matched controls and fails the static-support gate. The
second row is below the control mean. All membership and zero-sum triple
replays pass; the independent verifier reproduces both rows and all controls.

## Interpretation

This is a negative result for the finite `D_d(x,a)=c` family at the stated
degrees, parameters, constants, and toy curve. It does not rule out a
rational-map factor family with a different pre-symmetrization algebra, a
Frobenius construction with an explicit nonalgebraic correspondence, or the
separate product-Kummer line-bundle route. No source inverse, factor-log rank,
target descent, charged rho comparison, or ECDLP speedup is established.

## Next Candidate Set

1. Search maps with a bounded-support *fiber ideal* before symmetrization,
   rather than another recurrence-derived univariate polynomial.
2. Permit a Frobenius-conjugate construction only with an explicit
   nonalgebraic target-label path that clears the existing trace-zero audit.
3. Continue the product-Kummer normalized-Poincare section construction, then
   test a complete factor-base/relation/rank/descent path.

## Artifacts

- Contract: `notes/experiment_contract_dickson_fiber_triple_n611b.md`
- Probe: `tools/dickson_fiber_triple_probe_n611b.py`
- Receipt: `notes/dickson_fiber_triple_probe_n611b.json`
- Independent verifier: `tools/verify_dickson_fiber_triple_n611b.py`
- Verifier receipt: `notes/dickson_fiber_triple_probe_n611b_independent_verifier.json`
