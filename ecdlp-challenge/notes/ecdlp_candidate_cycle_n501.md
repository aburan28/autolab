# Candidate Cycle N501: post-selector ECDLP relation mechanisms

## Candidate A: higher-arity internal theta relations

### Status

HYPOTHESIS / MODEL-BOUND / TOY-EVIDENCE / NOVELTY-UNVERIFIED

### Claim or task

Use four or more signed theta factors in the genus-two cover so that internal
Prym cancellation has enough combinatorial supply to support factor logs and
blind target descent.

### Assumptions

- The theta factors admit exact Cantor composition and separate signed labels.
- A meet-in-the-middle collector can be fully charged, including product-state
  occupancy, memory, relation rank, and descent.

### Minimal falsification

On the nine existing split covers, compare arities three, four, and five with
matched random odd controls. Promote only if the relation and descent lift is
at least eightfold, factor-log rank grows with the signed base, and the
charged meet-in-the-middle cost is below rho.

### Likely failure mode

The Prym-zero condition remains an independent product-group constraint, so
higher arity buys only materialized tuple cost rather than a useful density
gain.

## Candidate B: non-split Mumford-coordinate factor base

### Status

CONJECTURE / MODEL-BOUND / NOVELTY-UNVERIFIED

### Claim or task

Define a factor base directly on reduced Mumford divisors of a non-split
genus-two Jacobian, aiming for a low-degree membership surface that avoids a
fixed selected sheet and retains an explicit target-label path.

### Assumptions

- The target embedding is specified before relation collection.
- The map back to the original ECDLP is neither a hidden fixed algebraic return
  nor an uncharged discrete-log oracle.

### Minimal falsification

Construct the smallest non-split toy cover and audit: target embedding,
composition, membership degree, relation rank, blind descent, and the full
cost of the target-label map. Reject immediately if the target path factors
through an invertible isogeny or loses its scalar label.

### Likely failure mode

The label path collapses to the algebraic-return obstruction or the factor-base
membership surface has product-group occupancy comparable to random sampling.

## Candidate C: compact joint selector correction

### Status

HYPOTHESIS / MODEL-BOUND / NOVELTY-UNVERIFIED

### Claim or task

Choose sheet signs only at relation time using a bounded public correction
function of the complete tuple, rather than a pointwise selector or the
zero-Prym balance rejected by REP-N499.

### Assumptions

- The correction has a fixed, public description independent of target and
  relation row identity.
- It is evaluable before offset opening and does not encode a posting list.

### Minimal falsification

Enumerate all bounded corrections from a fixed low-complexity basis on the
existing covers, validate train/holdout and wrong-target behavior, then charge
candidate generation, rank, and descent against matched rho.

### Likely failure mode

Every compact correction is either control-like or becomes a row-specific
payload when it gains true-row recall.

## Next action

Implement Candidate A first: it has the existing exact cover machinery and a
clear product-occupancy null hypothesis. Candidates B and C remain independent
proof/disproof tracks.
