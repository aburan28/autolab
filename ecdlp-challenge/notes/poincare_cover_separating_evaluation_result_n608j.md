# Result: N608J Cover-Separating H018 Fibre Evaluation

## Construction

Let `pi:Eprime -> E` be the N606U cyclic degree-nine cover and let `K` generate
its deck kernel.  The cover curve has a deterministic degree-eleven isogeny

```text
phi:Eprime -> E.
```

For each nonzero `Q` in `E(F_103)`, choose a lift `Qprime` in `pi^-1(Q)` and
evaluate the explicit N606Y fibre basis of

```text
H0(E, O(8O + [zQ]))
```

at the nine points

```text
phi(Qprime + jK),  j=0,...,8.
```

Changing the lift only permutes the labelled deck coordinates.  Since
`deg(phi)=11` is coprime to `|K|=9`, `phi(K)` has order nine and these nine
evaluation points are distinct.

## Result

Every one of the 108 nonzero rational base fibres is regular for this evaluator
and has a nonzero determinant:

| Map used for nine evaluation points | Rank distribution |
|---|---:|
| `phi` of degree 11 | `rank 9: 108` |
| `pi` of degree 9, negative control | `rank 1: 108` |

The independent Sage replay recomputes the same full panel, including deck
separation, matrix ranks, and the rank-one collapse control.

## Why This Matters

This is an explicit fibrewise isomorphism from the H018 fibre to nine scalar
evaluation coordinates.  It is not a finite interpolation gauge: both the
source basis and the nine evaluation points are named algebraically by `z`,
`pi`, `K`, and `phi`.

## Missing Step

The values currently live in scalar coordinates, whereas the N608I target is
`pi_*O_Eprime(2Oprime)`.  To obtain an actual bundle map, the next construction
must provide nonzero local weights

```text
w_j(Q) in O_Eprime(2Oprime)|_(Qprime+jK)
```

on declared charts and prove that the weighted matrices are regular and have
the normalized-Poincare transition law.  A determinant nonzero at all rational
base points does not prove this: its divisor can still have nonrational support,
and arbitrary pointwise rescaling would only recreate the rejected finite-gauge
escape.

## Boundary

`OBSERVATION / COVER-SEPARATING FIBRE EVALUATION / MODEL-BOUND / TOY-EVIDENCE /
NO BUNDLE MORPHISM / NO SECTION EVALUATOR / NO FACTOR BASE / NO_ECDLP_CLAIM`.
