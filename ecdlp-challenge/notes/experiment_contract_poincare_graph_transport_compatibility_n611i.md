# Experiment Contract: N611I Graph/Transport Compatibility

## Hypothesis

The degree-eleven graph evaluator used in N608O can intertwine the explicit
N606Y first-factor Poincare transport with the N608G degree-nine deck frame.
For a cover translation `L`, that requires

```text
pi(L) = [84] phi(L),
```

because the Poincare transport shifts the base fibre by `[84]R` when its
first-factor argument shifts by `R`.

## Null hypothesis

The graph map has a different translation slope. Then it may separate deck
sheets fibrewise but cannot serve as this normalized-Poincare/deck
intertwiner without a further non-scalar geometric correction.

## Parameters

- target: `E/F_103: y^2 = x^3 + x + 24`, of order `109`;
- cover: the cyclic degree-nine cover `pi:Eprime -> E` from N606U;
- graph map: the degree-eleven `phi:Eprime -> E` from N608O;
- Poincare base-shift multiplier: `84` from N606Y/N607F;
- all 108 nonzero `Eprime(F_103)` translations;
- all 108 nonzero rational base fibres, with sheet computations over
  `F_(103^6)`.

## Metrics

- scalar slope `s` defined by `pi = [s] phi` on `Eprime(F_103)`;
- number of nonzero translations satisfying `pi(L)=[84]phi(L)`;
- graph evaluation rank across all regular deck fibres;
- rank of the transport-aligned control `psi=[84]^-1 pi` across the same
  fibres.

## Positive control

`psi=[84]^-1 pi` satisfies `pi(L)=[84]psi(L)` identically. It factors through
`pi`, so its nine deck-sheet values must coincide and its moving-basis
evaluation matrix must have rank one.

## Negative control

The cover projection `pi` itself is deck-invariant and cannot provide a
rank-nine sheet evaluator.

## Success criterion

The graph map has slope `84` and passes the translation identity on every
nonzero rational cover translation while retaining rank nine on every regular
deck fibre. Only then may it advance to a literal vector-valued Poincare/deck
transition construction.

## Falsification criterion

Any non-`84` slope or failed translation identity rejects this graph map as a
direct intertwiner, even if its fibre matrices are full rank.

## Reproduction

```bash
sage -python tools/poincare_graph_transport_compatibility_n611i.py \
  --out notes/poincare_graph_transport_compatibility_n611i.json
sage -python tools/verify_poincare_graph_transport_compatibility_n611i.py \
  --primary notes/poincare_graph_transport_compatibility_n611i.json \
  --out notes/poincare_graph_transport_compatibility_n611i_independent_verifier.json
```

## Claim boundary

`HYPOTHESIS / GRAPH_TO_POINCARE_DECK_INTERTWINER / MODEL-BOUND /
TOY-EVIDENCE / NO_ECDLP_CLAIM`.
