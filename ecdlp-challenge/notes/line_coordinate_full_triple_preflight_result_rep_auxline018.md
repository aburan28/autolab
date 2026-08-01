# Result: REP-AUXLINE-018 Full Line-Coordinate Triple Preflight

## Status

NEGATIVE FULL-TRIPLE QUOTIENT-REMAINDER RESULT / MODEL-BOUND / TOY-EVIDENCE /
NO_ECDLP_CLAIM

## Candidate

For the direct oriented line-coordinate family

```text
F = { P : (y(P)-m*x(P)-r)^d=c },
```

let `C_{s,t}(x)` be the intersection cubic of `E` and `y=s*x+t`, and let
`D_{s,t}(x)=((s-m)*x+(t-r))^d-c`. The exact full-triple condition is
`D mod C = 0`: every intersection point belongs to the direct factor set, so
the three points sum to zero.

## Evidence

On `E/F_101: y^2=x^3+x+32`, using `m=7`, `r=11`, and every split degree
`5,10,20,25,50`, the quotient-remainder sources agree exactly with every
distinct direct-factor triple that sums to zero. The independent replay
reconstructed every source count, every random control, and every static
support count.

| d | oriented base | direct points | exact triples | random lift | remainder terms | pair state | terms / pair state |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 3 | 3 | 0 | n/a | 104 | 6 | 17.33 |
| 10 | 7 | 8 | 1 | 1.56 | 446 | 28 | 15.93 |
| 20 | 15 | 15 | 4 | 0.84 | 1,799 | 120 | 14.99 |
| 25 | 23 | 27 | 26 | 0.87 | 2,804 | 276 | 10.16 |
| 50 | 38 | 51 | 211 | 1.03 | 10,733 | 741 | 14.48 |

The exhaustive `p^2=10,201` line scan is about `809.88` times the toy rho
group-addition scale, and is not a comparable solver cost. More importantly,
each degree has a generic three-coefficient remainder whose static support
already exceeds the oriented pair state by more than a factor of ten. No row
has the preregistered fourfold density lift over 64 equal-size random point
controls.

## Interpretation

The line quotient provides a correct ternary relation identity, not a compact
relation source. It makes real zero-sum triples visible, but on this direct
factor family they are random-density at the tested scale and their defining
system is pair-scale-or-worse before source inversion. A Groebner/resultant
collector, rank computation, and target descent were therefore not
implemented.

## Limits

This is a negative result for the unmodified direct line-coordinate coset
family and the all-three-points affine-line quotient. It does not rule out a
new signed or multi-line factor representation, a different rational map with
bounded denominator-cleared support, a Frobenius/extension quotient,
product-Kummer geometry, another source inversion mechanism, or ECDLP
generally.

## Artifacts

- Contract: `notes/line_coordinate_full_triple_preflight_contract.md`
  (`122b3770784ad53125c821c786db72f7e869ee13719900def181de8f0d7f72a5`)
- Probe: `tools/line_coordinate_full_triple_preflight.py`
  (`51e33631b83ced72c4dd52a2bee5ccd2e8288c85c071ae860d5aacea57df9dfb`)
- Independent verifier: `tools/verify_line_coordinate_full_triple_preflight.py`
  (`c107261132a4d73c94b781785ba1562106e4e8b73642ef97f470bc11f16f96d2`)
- Primary artifact:
  `notes/line_coordinate_full_triple_preflight_rep_auxline018.json`
  (`f4ba1e0565d02479ffd026cb8979091d9c0e25b10f8288df5f37d12047b0dd42`)
- Independent replay:
  `notes/line_coordinate_full_triple_preflight_independent_verifier_rep_auxline018.json`
  (`6557804a7897729b88e48c00919297a99c4d9d666d1f8777ef61254eee3f84c7`)

## Next Concrete Action

Preflight a rational-map family only when its triple membership has bounded
support before denominator clearing, a public factor witness, and a non-scan
source inversion with a credible rank and blind-descent path.
