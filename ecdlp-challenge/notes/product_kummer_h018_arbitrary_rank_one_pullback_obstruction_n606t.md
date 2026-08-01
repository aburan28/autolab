# N606T: Arbitrary Rank-One Pullback Obstruction

## Status

RESTRICTED THEOREM / NEGATIVE FINITE_RANK_ONE_PULLBACK_DECOMPOSITION /
INTRINSIC_POINCARE_COUPLING_REQUIRED / MODEL-BOUND / TOY-EVIDENCE /
NO_ECDLP_CLAIM

## Claim

No finite product of theta sections pulled back along homomorphisms
`E x E -> E` can have H018 as its polarization class.

## Model

Let `R=Z[pi]`, with `pi^2+5*pi+103=0`.  A finite product of such pullbacks
has Hermitian class

```text
H = sum_j v_j^* v_j,    v_j=(a_j,b_j) in R^2.
```

Each summand is positive rank one.  H018 has diagonal entries `9,11` and
off-diagonal `2+pi`.

## Proof

For `u=a+b*pi`,

```text
N(u) = (a-5b/2)^2 + 387b^2/4.
```

Every nonintegral `u` therefore has `N(u) >= 97`.  In a positive sum equal to
H018, the diagonal equations are

```text
sum_j N(a_j)=9,    sum_j N(b_j)=11.
```

Thus every `a_j` and `b_j` has norm at most `11`, hence lies in `Z`.  Every
summand's off-diagonal `conj(a_j)b_j` is then integral, and so is their sum.
It cannot equal `2+pi`.

This argument is independent of the number of factors.  It extends N601A and
N606Q from two-coordinate pullbacks to every finite positive rank-one
pullback decomposition.

## Boundary

The theorem excludes only sections assembled from finite products of
homomorphism-pullback theta factors.  It does not exclude a genuinely
non-factorizable Poincare/theta section, a non-split principally polarized
quotient over `F_103^6`, a refined-cover evaluator, relation rank, target
descent, or an ECDLP improvement.

## Next Action

Any positive H018 construction must provide an explicit non-factorizable
Poincare coupling or non-split genus-two quotient model, then test its section
evaluator against base-locus, relation-rank, descent, and fully charged rho
baselines.
