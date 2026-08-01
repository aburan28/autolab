# Result: N608Y finds diagonal symmetry but rejects finite-interpolation admission

## Status

`NEGATIVE RESULT / DIAGONAL KUMMER DESCENT WITHOUT FINITE-INTERPOLATION
ADVANTAGE / STRUCTURAL REPLAY / MODEL-BOUND / TOY-EVIDENCE / NORMALIZED GLOBAL
MODEL OPEN / NO_ECDLP_CLAIM`

## Result

The complete N608U `t=0` source set has 84 points. It is exactly closed under
simultaneous negation:

```text
(P,Q) -> (-P,-Q).
```

Neither one-sided sign involution preserves any source. The quotient therefore
has 42 distinct labels in diagonal product-Kummer coordinates
`(x(P),x(Q),y(P)y(Q))`.

This symmetry alone does not identify the level curve. Scanning monomials
`u^i x^j` and `w u^i x^j` finds the first interpolation kernel in the
`(u,x)=(3,5)` box: 48 monomials, evaluation rank 42, and kernel dimension 6.
All sixteen deterministic random 42-point subsets of the same 5,832-point
rational product-Kummer locus produce exactly that same first box and rank.

## Interpretation

Diagonal Kummer descent is a genuine structural property of this toy level,
but a rational-point interpolant in the first available kernel is no evidence
of a global equation or a special low-complexity curve. The observed kernel is
fully explained by the number of samples.

The next viable path is to derive a normalized global product-Kummer or
Poincare model from explicit transition data, then independently verify the
resulting divisor, source return, and fixed-return boundary. This receipt does
not address factor bases, relations, rank, descent, cost, or ECDLP speed.
