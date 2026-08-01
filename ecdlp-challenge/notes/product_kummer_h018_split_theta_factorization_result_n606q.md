# Result: N606Q Split Theta-Factorization Gate

## Status

RESTRICTED THEOREM / NEGATIVE SPLIT_THETA_FACTORIZATION_GATE /
INTRINSIC_NONSPLIT_THETA_DATA_REQUIRED / MODEL-BOUND / TOY-EVIDENCE /
NO_ECDLP_CLAIM

## Claim

Test whether the H018 polarization is a CM-coordinate pullback of the split
type-(1,2) product polarization:

```text
H018 = A^* diag(1,2) A,
H018 = [[9, 2+pi], [-3-pi, 11]].
```

Such a factorization would have provided an explicit product-theta starting
point for the two scalar sections required by N606P.

## Exact Obstruction

For `u=a+b*pi` in the fixture CM order,

```text
N(u) = a^2 - 5ab + 103b^2
     = (a - 5b/2)^2 + 387b^2/4.
```

Every element with nonzero `pi` coefficient therefore has norm at least `97`.
If `A=[[a,b],[c,d]]` factored H018 as above, its diagonal entries would give

```text
N(a)+2N(c)=9,
N(b)+2N(d)=11.
```

All four entries would then be integral. Its off-diagonal
`conj(a)*b + 2*conj(c)*d` would be integral too, contradicting `2+pi`.

The explicit bounded search agrees: the six possible first columns and four
possible second columns yield 24 candidates, whose off-diagonals are only
`-9,-7,-1,1,7,9`; none is `2+pi`. The split positive control
`diag(1,2)=I^*diag(1,2)I` passes.

## Interpretation

This rules out only a CM-coordinate change from a split product-theta basis.
It does not rule out the positive H018 class, its two class-level sections, an
intrinsically non-split theta/Poincare construction, a refined-cover evaluator,
or any later factor-base, rank, descent, or ECDLP result.

## Next Concrete Action

The next N606Q phase must construct intrinsically non-split theta/Poincare
sections with named scalar sections `s0,s1`, explicit rigidification anchors,
and an evaluator on the refined cover. It cannot substitute a product-theta
pullback.
