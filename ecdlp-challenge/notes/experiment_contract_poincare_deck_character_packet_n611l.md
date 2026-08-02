# Experiment Contract: N611L Deck-Character Packet

## Hypothesis

For a complete H018 source pair `(P,Q)`, a Frobenius-descended nontrivial
character packet of the nine `pi`-deck lifts of `Q`, translated by `-P`, has a
smaller correction alphabet than a sourcewise randomized deck ordering while
still depending on more than `Q` alone.

## Packet

For any lift `R` of `Q`, an order-nine deck generator `K`, a primitive ninth
root `zeta`, and `f(S)=1/(x(S)-1)` with `f(O)=0`, define

```text
A(P,Q;R) = sum_{j=0}^8 zeta^(-j) f(phi(R+jK)-P),
L(P,Q) = Norm_F103^6/F103(A(P,Q;R)^9).
```

The ninth power removes a deck-origin shift; the absolute norm is the declared
Frobenius descent. A packet with a pole is retained as an explicit failure,
not silently dropped.

## Parameters

- N608U complete sources at levels `0,1,17` on the registered H018 fixture;
- all nonzero rational `Q` and all nine deck lifts over `F_(103^6)`;
- fixed rational probe `f=1/(x-1)`; no probe selection after observing data;
- a deterministic sourcewise Fisher-Yates permutation of the same nine values
  as the matched control.

## Metrics

- complete-packet validity and deck-lift invariance;
- base-field descent of `L`;
- label support and maximum collision multiplicity;
- number of Q fibres with more than one packet label;
- candidate/control support ratio.

## Positive controls

- each candidate label is invariant under all nine changes of lift;
- the norm is Frobenius fixed;
- sourcewise randomized orderings use exactly the same packet values.

## Success criterion

All complete sources have valid lift-invariant labels, at least one Q fibre has
more than one label, and the candidate support is at most half the matched
random-ordering support at every level. Passing is only a representation
signal; a target-bearing relation, rank, and descent remain required.

## Falsification criterion

Invalid packets, Q-only collapse, or failure to beat the matched random
ordering rejects this fixed character packet. It does not rule out other
characters, probes, covers, or a Jacobian-normalization route.

## Reproduction

```bash
sage -python tools/poincare_deck_character_packet_n611l.py \
  --out notes/poincare_deck_character_packet_n611l.json
sage -python tools/verify_poincare_deck_character_packet_n611l.py \
  --primary notes/poincare_deck_character_packet_n611l.json \
  --out notes/poincare_deck_character_packet_n611l_independent_verifier.json
```
