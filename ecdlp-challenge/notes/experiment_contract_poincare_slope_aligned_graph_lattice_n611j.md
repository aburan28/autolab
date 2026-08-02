# Experiment Contract: N611J Slope-Aligned Graph Lattice

## Hypothesis

An integral graph map

```text
chi_(a,b) = [a]phi + [b]pi
```

can match the N606Y Poincare translation slope `84`, retain rank-nine deck
separation, and remain in the compact degree-three graph regime needed by the
N608O elementary modification.

## Parameters

- `E/F_103` and `Eprime` as in N608M;
- `h_(a,b)=(chi_(a,b),pi)`;
- exact N608M pullback-degree form
  `99a^2 + 27ab + 81b^2 - 195a - 9b + 99`;
- N611I rational translation slope `pi=[50]phi`;
- Poincare slope `84` and all 108 nonzero rational cover translations;
- all 108 nonzero base fibres over `F_(103^6)`.

## Exact constraint

Transport compatibility requires

```text
84a + 58b = 50 mod 109.
```

Deck separation requires `a` to be a unit modulo `9`, since `pi` vanishes on
the deck kernel and `phi` is injective there.

## Controls

- `(a,b)=(1,0)` is the N608O degree-three rank-nine graph control and must
  fail the slope-84 transport identity.
- `(a,b)=(0,61)` is transport-aligned but factors through `pi`, so it must
  have rank one on every deck fibre.
- The candidate minimizer must match the transport identity and retain rank
  nine on every regular fibre.

## Success criterion

At least one slope-aligned deck-separating pair has pullback degree at most
three. Only then can it replace the compact N608O graph directly.

## Falsification criterion

If every aligned deck-separating integral graph has degree greater than three,
the compact direct-integral correction route is rejected. A higher-degree map
is recorded as an algebraic observation, not promoted to a section or ECDLP
mechanism.

## Reproduction

```bash
sage -python tools/poincare_slope_aligned_graph_lattice_n611j.py \
  --out notes/poincare_slope_aligned_graph_lattice_n611j.json
sage -python tools/verify_poincare_slope_aligned_graph_lattice_n611j.py \
  --primary notes/poincare_slope_aligned_graph_lattice_n611j.json \
  --out notes/poincare_slope_aligned_graph_lattice_n611j_independent_verifier.json
```
