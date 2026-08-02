# Result: N611I Graph/Transport Compatibility

## Status

`NEGATIVE DIRECT_GRAPH_POINCARE_DECK_INTERTWINER / MODEL-BOUND /
TOY-EVIDENCE / VECTOR_VALUED_CORRECTION_OPEN / NO_ECDLP_CLAIM`

## Exact Result

The N608O degree-eleven graph map `phi:Eprime -> E` and the N606U degree-nine
cover `pi:Eprime -> E` obey the rational-point scalar relation

```text
pi = [50] phi on Eprime(F_103).
```

The N606Y/N607F first-factor Poincare transport instead requires a graph map
with translation slope `84`, because a first-factor shift by `R` sends the
base fibre from `Q` to `Q + [84]R`. Consequently, none of the 108 nonzero
rational cover translations satisfies

```text
pi(L) = [84] phi(L).
```

This is not a sheet-separation failure. On all 108 nonzero rational base
fibres, the moving-basis graph evaluation at the nine deck sheets has exact
rank nine.

The transport-aligned control

```text
psi = [84]^-1 pi = [61] pi
```

satisfies the Poincare translation identity on all 108 cover translations,
but it factors through `pi`. Its values are constant on every deck fibre, and
the same moving-basis evaluation has rank one on all 108 nonzero fibres.

## Evidence

- Producer controls: `4/5` pass; the sole failed gate is the intended direct
  graph/Poincare transport identity.
- Independent replay: `7/7` checks pass.
- The producer exhausts the 108 nonzero `Eprime(F_103)` translations and the
  108 nonzero rational deck fibres; fibre matrices use `F_(103^6)` so the
  order-nine deck generator is present.

## Interpretation

The known graph evaluator has exactly the desired fibrewise separation but the
wrong translation slope for the known Poincare transport. The obvious
transport-aligned replacement destroys that separation. Therefore no scalar
relabeling of this graph can supply the missing normalized-Poincare/deck
intertwiner.

This is a restricted negative result for the named graph map and transport
law, not an obstruction to a genuinely vector-valued transition, a different
correspondence, or a non-scalar correction that changes the slope while
preserving rank nine.

## Next Concrete Action

Construct a matrix-valued correction between the N607E Cauchy transport and
the N608G deck frame. It must change the graph slope from `50` to `84`, remain
evaluation-independent on a declared open cover, preserve rank nine on all
regular fibres, and then be tested as literal Cech data.

## Reproduction

```bash
sage -python tools/poincare_graph_transport_compatibility_n611i.py \
  --out notes/poincare_graph_transport_compatibility_n611i.json
sage -python tools/verify_poincare_graph_transport_compatibility_n611i.py \
  --primary notes/poincare_graph_transport_compatibility_n611i.json \
  --out notes/poincare_graph_transport_compatibility_n611i_independent_verifier.json
```
