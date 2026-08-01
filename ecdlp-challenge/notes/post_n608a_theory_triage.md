# Post-N608A Candidate Triage

`OPEN / HYPOTHESIS / NOVELTY-UNVERIFIED / MODEL-BOUND.` This note does not
claim an ECDLP improvement. It selects distinct next tests after N608A showed
that an exact univariate product index preserves the raw pair-plus-large image
degree.

## Baseline and admission rule

The relevant baseline is Pollard rho on a prime-order base-field subgroup.
Every candidate below must provide all of: target-independent factor-base
construction, exact source replay, a relation matrix with rank, the same-form
target descent, and fully charged time/memory below rho. A toy section space,
membership filter, or verifier receipt alone is not an attack.

## 1. Conservative extension: bilinear large-prime quotient certificate

### Core mechanism

Replace N608A's univariate x-image index by a two-coordinate certificate
`(u,v)` derived from the pair sum and large prime, where the relation
`P_i+P_j+L=W` becomes a bidegree-bounded intersection rather than a product
polynomial over all images. The certificate must return the signed pair and
large-prime witness without a `B^2 S` posting table.

### Why it is distinct

N608A only rejects a square-free univariate product over target x-values. A
bounded bivariate ideal with a source inverse would be a different operation.
The existing low-x trace/norm work is nearby but does not supply this
pair-plus-large witness or large-prime collision descent.

### First falsifying test

On a 20-bit generated curve, define `u=x(P_i+P_j)` and
`v=x(P_i+P_j+L)`, eliminate only the two pair coordinates, and measure:

- bidegrees and monomial count of the exact ideal;
- number of source witnesses per accepted target;
- storage needed to recover the signed pair and `L`;
- random-control occupancy and collision-row rank.

Reject this candidate if eliminating the pair variables has degree or support
`Omega(B^2 S)`, or if the inverse stores one posting per image witness.

### Likely failure

The second coordinate makes source recovery explicit but merely converts the
same product image into a bivariate table.

## 2. Representation change: non-split genus-two elliptic component

### Core mechanism

Use a genus-two Jacobian with an embedded elliptic component that is not the
split `E x E'` coordinate model screened by REP-N500. Define a theta-slice
factor base in reduced Mumford coordinates and test whether sums meeting the
elliptic component have a non-product-group occupancy or a low-degree Cantor
decomposition certificate.

### Why it is distinct

REP-N500 rejects the explicit split Prym cancellation equation. It does not
settle a non-split Jacobian whose theta divisor, embedding, and descent do not
separate into one base coordinate and one independent Prym coordinate.

### First falsifying test

Materialize one public genus-two toy curve with an explicit elliptic isogeny
factor and enumerate a small reduced-divisor theta slice. Compare its exact
triple intersection with the embedded elliptic target coset to matched random
Jacobians/subsets of equal size. Record Cantor operation counts, Mumford
degrees, fibre multiplicities, and a target-shaped inverse.

Reject it if the target intersection matches the product occupancy model, or
if target descent requires enumerating a Prym coordinate or a quadratic-size
divisor table.

### Likely failure

The non-split model becomes split after choosing a computational basis, and
the missing Prym state reappears in source recovery.

## 3. High-risk speculative idea: type-(1,2) Poincare pencil packet

### Core mechanism

Complete the H018 type-(1,2) Poincare construction to a genuine two-section
pencil, then ask whether a smooth member gives an explicit degree-three packet
correspondence whose packet sums recover base-group relations without a
product-state inverse.

### Why it is distinct

This is a surface-level theta construction, not an x-image filter, a finite
orbit frame, or a split-Jacobian Prym cancellation. Current H018 work has only
formal class and local-frame evidence; it has no section, curve, packet map,
factor base, rank, or descent.

### First falsifying test

Construct either an explicit normalized-Poincare/pushforward bundle morphism
or prove that the selected local reflexive modules cannot realize the H018
class on a genuine cover. If a pencil is materialized, immediately require a
smooth-member check, packet quotient maps, a public factor base, and a
source-to-target cost inequality before collection.

### Likely failure

The surface construction exists but its packet inverse is a disguised
quadratic table, or its smooth members do not carry an addition-compatible
factorization.

## Decision

Run candidate 1 first because it has the smallest exact falsification test and
directly targets the only unclosed large-prime escape hatch. Candidate 2 is
the most promising representation change but needs an explicit non-split toy
model before an occupancy experiment. Candidate 3 remains the highest-upside
geometric route, with its current Poincare rigidification prerequisite still
active.

## Handoff: post-N608A theory slate

### Claim or task

Choose mathematically distinct successor tests after the raw large-prime
product-index representation fails to compress.

### Status

OPEN

### Assumptions

- All work remains on generated or toy public curves.
- No candidate is novel or an ECDLP improvement without the full route.

### Evidence so far

- N608A rejects the exact univariate product index.
- REP-N500 rejects the stated split-Jacobian theta slice.
- H018 lacks a realized Poincare bundle morphism and any relation gate.

### Failure modes

- A certificate or theta construction can hide a full source table.
- An internal Jacobian can reintroduce a product-state descent coordinate.

### Next concrete action

Write the N608B contract for the bilinear large-prime quotient certificate.

### Artifact paths

- `notes/post_n608a_theory_triage.md`
