# Result: REP-AUXLINE-017 Ternary Coset Symmetric Preflight

## Status

NEGATIVE DIRECT-TERNARY-COSSET REPRESENTATION RESULT / MODEL-BOUND /
TOY-EVIDENCE / NO_ECDLP_CLAIM

## Exact membership identity

For three roots `x1,x2,x3` of `x^d=c`, let `(e1,e2,e3)` be their elementary
symmetric coordinates and let `P_k=x1^k+x2^k+x3^k`. The following are
equivalent in odd characteristic when `d` is below the field characteristic:

```text
x1^d=x2^d=x3^d=c

P_d(e1,e2,e3)=3c
P_2d(e1,e2,e3)=3c^2
e3^d=c^3.
```

The reverse implication follows because the three values `ui=xi^d` have
elementary symmetric coordinates `(3c,3c^2,c^3)`, hence characteristic
polynomial `(T-c)^3`.

## Preflight result

The S4 target polynomial is quartic in the third elementary coordinate `e3`:
three independently reconstructed prime-order toy curves each have 439 raw
S4 monomials, 118 elementary-coordinate monomials, and exact `deg_e3(S4)=4`.
That bounded target quotient does not help the direct ternary coset route,
because `P_2d` has quadratic support.

| fiber | d | factor base B | `P_d` terms | `P_2d` terms | pair state | ratio |
|---|---:|---:|---:|---:|---:|---:|
| 20-bit, coset 1 | 190 | 88 | 3,104 | 12,224 | 3,916 | 3.12 |
| 22-bit, coset 0 | 486 | 218 | 19,927 | 79,219 | 23,871 | 3.32 |

For `P_k`, every weighted partition `a+2b+3c=k` gives a nonzero Girard
coefficient in the registered fields. Therefore the exact `P_2d` support is
Theta(`d^2`) before target specialization, already larger than the complete
unordered pair-coordinate state in both fibers. This fails the stated
admission gate; no S4 ternary source solver was promoted.

## Evidence

- Contract: `notes/compositional_coset_ternary_preflight_contract.md`
  (`31e6cd5a9e4d4d2d237950c9fe412b32f93db7a7613fd4b1f15c786e6b8f355b`)
- Probe: `tools/compositional_coset_ternary_preflight.py`
  (`82e4365fed4e7ea7977a2ad1457fabad9bc726156abab367d2dcb64e18408d51`)
- Independent verifier: `tools/verify_compositional_coset_ternary_preflight.py`
  (`7f280c711b559988bc78b8ceb4f236f38fea971b4b6cf35c55f0057ef5c20553`)
- Primary report and replay:
  `ead19c31f8e110bc9a43dc3967672ec7f33fbd89c0d89e71f43f40ce5e85c18e`,
  `5de79fdd87b4e7a802313092a2aed49b7984ec76c5f74cb18368fd3e70874933`

## Limits

This closes only the direct elementary-symmetric ternary representation of a
multiplicative-coset factor base. It does not rule out a factor family with a
different low-support triple identity, an extension-field/Frobenius quotient,
an endpoint-producing S4 inverse, product-Kummer geometry, relation-rank
amortization, target descent, or ECDLP generally.

## Next concrete action

The next candidate must avoid the `P_2d` dense invariant itself. Test factor
families whose triple membership is specified by a bounded-support algebraic
map before symmetrization, beginning with a public rational-map search that
requires an exact factor witness, nonrandom triple density, and a complete
source-to-rank-to-descent cost model.
