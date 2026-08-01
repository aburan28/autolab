# Result: N608Z extracts a common base graph from cleared H018 levels

## Status

`OBSERVATION / COMMON BASE GRAPH EXTRACTED FROM CLEARED H018 LEVELS /
STRUCTURAL REPLAY / MODEL-BOUND / TOY-EVIDENCE / GLOBAL DIVISOR MODEL OPEN /
NO_ECDLP_CLAIM`

## Result

For every nonzero `Q` and each declared level `t=0,1,17`, the N608U
denominator-cleared equation has the common solution

```text
P = 4Q.
```

It is an actual solution of the unsquared equation `A_Q,t(x)+B_Q(x)y=0`, not
merely an x-coordinate artefact. The point `P=-4Q`, which has the same x
coordinate, never satisfies the unsquared equation. Reversing the Cauchy
numerator sign rejects the graph in all 108 `t=0` rows.

The degree-ten eliminant has an exact `x-x(4Q)` factor in every row. Dividing
it out leaves a degree-nine quotient for every `Q` and every declared level.
The x-root multiplicity is usually one; it is two in 8 of 108 rows at `t=0`
and 2 of 108 rows at `t=17`, recording residual intersection with the base
graph rather than erasing it from the global geometry.

## Interpretation

N608U's degree-ten inverse includes a common graph created by clearing the
moving Cauchy denominator. A residual degree-nine inverse can now be stated
and checked on the open source set. This is a local source-generation
simplification only. It neither constructs the compact residual divisor nor
supplies relations, factor logs, rank, descent, a charged cost comparison, or
an ECDLP speedup.

## Next Action

Implement the degree-nine residual inverse, compare it against exhaustive
open-source enumeration, and then derive the compact residual divisor and its
global class rather than extrapolating from the affine factorization.
