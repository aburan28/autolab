# Result: N611L Deck-Character Packet

## Status

`NEGATIVE FIXED_DECK_CHARACTER_PACKET / MODEL-BOUND / TOY-EVIDENCE /
NORMALIZATION_OR_JACOBIAN_PACKET_OPEN / NO_ECDLP_CLAIM`

## Exact Result

The declared packet

```text
L(P,Q) = Norm_F(103^6)/F(103)(sum zeta^(-j) f(phi(R+jK)-P))^9
```

with fixed `f(S)=1/(x(S)-1)` is invariant under all nine changes of the deck
lift and is not a function of `Q` alone.  It does not meet its contract.

At complete source levels `0, 1, 17`, respectively, its label supports are
`22, 18, 21`, versus sourcewise randomized-order controls `30, 33, 34`.
The candidate/control support ratios are therefore
`0.73333..., 0.54545..., 0.61764...`, all above the required `0.5` gate.
The level-one panel also has two explicit packet poles, so only `106/108`
complete source pairs have valid labels there.  Levels zero and seventeen are
valid on `84/84` and `144/144` source pairs.  The packet has more than one
label on `14, 26, 40` Q-fibres, respectively.

## Evidence

- The producer retains the receipt and then fails preflight because the fixed
  probe has two poles; that failure is required by the contract rather than a
  dropped observation.
- The independent replay passes all seven checks: schema, complete-source
  counts, packet counts, validity, lift invariance, support statistics, and
  Q-fibre dependence.
- The character mechanics remain exact: zero deck-lift invariance failures
  occur on every valid packet and the absolute norm is Frobenius fixed.

## Interpretation

This rejects only the declared character packet and fixed rational probe.  Its
partial alphabet compression does not beat a matched random ordering, and its
two poles prevent complete-source use.  Tuning a different probe after this
panel would be selection on the observed data, so this result leaves that
route closed as stated.  It provides no target-bearing relation, relation
rank, individual-log descent, charged cost comparison, or ECDLP speedup.

## Next Concrete Action

Construct a canonical normalization or Jacobian packet whose labels arise
from a globally specified divisor or reflexive module, then require complete
source return, a target relation, rank, descent, and charged sub-rho cost
before treating any representation signal as ECDLP progress.

## Reproduction

```bash
sage -python tools/poincare_deck_character_packet_n611l.py \
  --out notes/poincare_deck_character_packet_n611l.json || true
sage -python tools/verify_poincare_deck_character_packet_n611l.py \
  --primary notes/poincare_deck_character_packet_n611l.json \
  --out notes/poincare_deck_character_packet_n611l_independent_verifier.json
```
