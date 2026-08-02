# Experiment Contract: N611K Low-Degree Graph Map Screen

## Hypothesis

A low-degree prime isogeny `psi:Eprime -> E` outside the previously replayed
`11,17,23,41` controls aligns the N606Y Poincare transport multiplier `84`
while retaining all-sheet deck separation.

## Parameters

- target `E/F_103: y^2=x^3+x+24` and N606U cover `pi:Eprime -> E`;
- every prime degree below `100` for which Sage returns an isogeny from
  `Eprime` to a curve isomorphic to `E`;
- N608M basis maps `phi` (degree 11) and `pi` (degree 9);
- rational-orbit norm form `11a^2+3ab+9b^2` used to identify the observed
  low-degree maps as `[a]phi+[b]pi` on `Eprime(F_103)`;
- all 108 nonzero rational cover translations and all 108 nonzero base fibres;
- deck computations over `F_(103^6)`.

## Metrics

- discovered prime-degree maps and their integral-combination identifiers;
- transport matches `pi(L)=[84]psi(L)` on the full rational cover orbit;
- rank of the nine-sheet moving-basis evaluator on every nonzero base fibre;
- H018 pullback degree of the corresponding graph pair.

## Positive controls

- degree 11 `phi` must replay as `(a,b)=(1,0)` and retain rank nine;
- degrees 17, 23, and 41 must replay as the N608M combinations;
- every map with `gcd(a,9)=1` must act injectively on the order-nine deck
  kernel.

## Success criterion

At least one discovered low-degree prime isogeny matches all 108 transport
translations and keeps rank nine on every regular fibre. Only then may it
advance as a compact graph candidate.

## Falsification criterion

If all discovered maps are transport-misaligned, this bounds only the
prime-degree-below-100 graph-map family. It does not rule out higher-degree,
non-isogeny, vector-valued, or non-graph corrections.

## Reproduction

```bash
sage -python tools/poincare_low_degree_graph_map_screen_n611k.py \
  --out notes/poincare_low_degree_graph_map_screen_n611k.json
sage -python tools/verify_poincare_low_degree_graph_map_screen_n611k.py \
  --primary notes/poincare_low_degree_graph_map_screen_n611k.json \
  --out notes/poincare_low_degree_graph_map_screen_n611k_independent_verifier.json
```
