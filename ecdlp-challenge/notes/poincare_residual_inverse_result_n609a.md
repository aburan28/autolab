# Result: N609A upgrades the complete open-level inverse to degree nine

## Status

`OBSERVATION / DEGREE-NINE COMPLETE OPEN-LEVEL SOURCE INVERSE /
STRUCTURAL REPLAY / MODEL-BOUND / TOY-EVIDENCE / RELATION MECHANISM OPEN /
NO_ECDLP_CLAIM`

## Result

N608Z identifies the common cleared base graph `P=4Q`. N609A divides its
factor `x-x(4Q)` from every N608U eliminant and factors the degree-nine
residual instead. For each of 108 nonzero base points, the residual inverse
matches exhaustive open-surface enumeration exactly:

| Level | Residual inverse sources | Exhaustive open sources | Degree | Factor calls |
|---:|---:|---:|---:|---:|
| 0 | 84 | 84 | 9 | 108 |
| 1 | 108 | 108 | 9 | 108 |
| 17 | 144 | 144 | 9 | 108 |

The residual still contains an intersection with the removed graph in 8 rows
at level `0` and 2 rows at level `17`; those graph points are explicitly
excluded from the declared open source set. Reversing the Cauchy numerator
sign makes the unsquared graph identity fail in all `t=0` rows.

## Interpretation

This is a measured reduction in the algebraic source-return subroutine, from
one degree-ten to one degree-nine factorization per nonzero `Q`. It does not
yet change the base-point scan, define a factor base, create relations, or
provide an ECDLP complexity improvement.

## Next Action

Compactify the residual degree-nine divisor after base-graph removal, determine
its global class and boundary, and test whether it has a target-bearing
relation mechanism that survives rank, descent, and fully charged cost checks.
