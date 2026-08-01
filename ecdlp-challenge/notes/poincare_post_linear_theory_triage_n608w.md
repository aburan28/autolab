# Post-N608V theory triage: complete source levels after linear compression fails

## Status

`OPEN THEORY SLATE / NOVELTY-UNVERIFIED / MODEL-BOUND / NO_ECDLP_CLAIM`

N608U gives complete source sets for declared H018 pencil levels. N608V shows
that every product-linear map `E x E -> E` tested on those sources retains a
base-field-scale image. The next candidates must add structure rather than
repeat a linear projection.

## Theory 1: Kummer quotient control

### Core mechanism

Apply the `R ~ -R` quotient after a product-linear map and use the x-coordinate
as the factor label. This is the conservative extension of N608V: it tests
whether sign identification gives more than its unavoidable factor-of-two
compression.

### Why it might evade the current barrier

The product image might be highly negation-paired on a special pencil level,
so a Kummer factor base could be substantially smaller than the group image.

### Minimal test

For every N608V direction and each declared level, compare the number of
distinct x-labels with the exact sign-orbit lower bound for that direction's
group image. A viable result needs collisions materially beyond sign pairing
and a signed relation-return rule.

### Preliminary exploration

For the N608V best directions, the Kummer supports are `25`, `28`, and `34`
for levels `0`, `1`, and `17`, exactly the ceilings of half the corresponding
group supports `49`, `56`, and `68`. This is not yet a contract result; it
motivates a full-direction control rather than a factor-base claim.

### Likely failure and reusable lesson

The quotient may only identify signs, leaving support of order `q`. If so,
the useful lesson is that a scalar Kummer label cannot provide asymptotic
compression of the complete-source levels.

## Theory 2: Normalize a complete pencil level and use its Jacobian internally

### Core mechanism

Construct the normalization `C_t` of a declared complete level
`lambda(P,Q)=t` and make its projection/divisor structure explicit. Search for
a factor system and relations inside `Jac(C_t)`, rather than returning source
pairs through a map `E x E -> E`.

### Why it might evade the current barrier

N608V constrains product-linear maps only. A non-product curve can have
different divisor geometry and may expose a relation representation not visible
in either elliptic factor separately.

### Minimal test

Derive a reproducible plane or correspondence model for `C_0`, calculate its
projection degrees and singularities, and identify the induced maps on divisor
classes. Reject the route immediately if its target recovery factors through a
fixed algebraic map from an abelian variety back to `E`: the existing
Jacobian-return theorem makes that component a public scalar or zero on the
prime subgroup.

### Likely failure and reusable lesson

The source curve may merely lift the original target through a fixed return or
require a factor base of order `q`. Either outcome specifies which extra
divisor geometry a viable level curve would need.

## Theory 3: deck-character packet with Frobenius descent

### Core mechanism

Lift a complete level source through the degree-nine cover, form character
packets of the deck orbit over a field containing ninth roots of unity, and
seek a Frobenius-stable packet whose addition defect has a smaller correction
alphabet than the raw source pairs.

### Why it might evade the current barrier

The H018 cover supplies a genuine cyclic deck action and nontrivial Frobenius
multiplier. Character packets retain data discarded by product-linear maps and
could, in principle, organize source lifts into lower-dimensional correction
classes.

### Minimal test

Specify an exact lift rule for every N608U complete source, enumerate all deck
packets on the registered fixture, and compare additive-defect support and
collision probability against matched randomized deck labels. Promotion needs
an exact target-bearing relation rule and a correction alphabet smaller than
the proposed factor base.

### Likely failure and reusable lesson

The packet may be only a relabeling of the cover state, with defects comparable
to randomized labels. That would narrow deck character routes while leaving a
different cover or a non-character correspondence open.

## Selected Next Action

Start Theory 2 by deriving a canonical finite model of `C_0` from the N608U
Cauchy-cleared equations. Its contract must state projection degrees, an exact
source-return map, a fixed-return check, and the required relation/rank/descent
measurements before any algorithmic claim.
