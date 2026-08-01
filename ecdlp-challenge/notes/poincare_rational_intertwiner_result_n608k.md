# Result: N608K Rational Local H018-to-Cover Morphism

## Named Map

For a cover point `R` set `Q=pi(R)`.  On the open locus where

```text
R != Oprime,
phi(R) != O,
phi(R) != zQ,
phi(R) != -zQ,
```

the N606Y divisor model identifies the H018 fibre with

```text
F_Q = H0(E, O(8O+[zQ])).
```

The section `1` trivializes the source line away from its divisor and also
trivializes `O_Eprime(2Oprime)` away from `Oprime`.  Therefore the formula

```text
f |-> f(phi(R)) * 1|R
```

is a named regular functional from `pi^*F` to `O_Eprime(2Oprime)` on that
open locus.  Adjunction gives a rational local morphism

```text
F --> pi_* O_Eprime(2Oprime).
```

No interpolation coefficients or per-point rescalings are used.

## Exact Test

On all 108 nonzero `F_103` base points, every one of the nine cover positions
lies in the declared open locus.  The corresponding nine-by-nine matrix has
nonzero determinant at every fibre, so the rational local morphism has generic
rank nine.  The independently recomputed result agrees.

The negative control is essential: replacing `phi` by `pi` maps every deck
fibre to one point and yields rank one, so the full rank depends on the
degree-eleven cover-separating map rather than on arbitrary matrix fitting.

## What Is Still Open

This is not yet a global normalized-Poincare intertwiner.  The construction
uses the `1` trivialization only away from the finite source and target divisors.
The next proof must compute the determinant divisor of this morphism and give
regular transition matrices on complementary charts.  In particular, it must
not repair zeros or poles by fitting independent weights at the 108 rational
fibres.

## Boundary

`OBSERVATION / EXPLICIT RATIONAL LOCAL MORPHISM / MODEL-BOUND / TOY-EVIDENCE /
GLOBAL NORMALIZED-POINCARE GLUING OPEN / NO GLOBAL SECTION EVALUATOR /
NO FACTOR BASE / NO RELATION RANK / NO DESCENT / NO_ECDLP_CLAIM`.
