# Experiment Contract: N606X Actual Challenge-Curve Special-Structure Audit

## Hypothesis

The published 40-bit benchmark curve may expose a rational small-prime-degree
self-isogeny or a same-`j` base-field near miss not visible in the prior
20--22 bit generated controls.  A self-map would be only a candidate
GLV/GLS-style representation route, requiring a separate subgroup action and
fully charged scalar-decomposition comparison.

## Null hypothesis

No rational self-isogeny or same-`j` base-field near miss occurs through
prime degree 43 on the published 40-bit fixture.  This is a bounded negative
result, not a claim about higher-degree, non-rational, extension-field, or
other non-generic algorithms.

## Parameters

- Public fixture: `p=616883774851`, `a=24569641080`, `b=109798974509`,
  prime subgroup order `n=616882790773`, seed `0x12345678`.
- Isogeny degrees: primes 2 through 43 used by REP-N479.
- Embedding-degree guard: exponents 1 through 200.

## Controls

- The registered 20--22 bit REP-N479 rows are the negative family controls.
- Cardinality, primality, ordinary status, non-anomalous order, and
  `j not in {0,1728}` must be recomputed rather than inherited from JSON.

## Metrics

- Exact cardinality and matching published order.
- Embedding-degree guard, `j`, supersingularity, and endomorphism-order
  discriminant.
- Rational self-isogenies and same-`j` nonisomorphic codomains at each degree.
- Runtime charged separately as structural reconnaissance, not solver work.

## Success criterion

Find a reproducible rational self-isogeny and then open a new contract for its
action on the published prime-order subgroup and its cost against Pollard rho.

## Falsification criterion

No map in the registered degree range rejects only this low-degree rational
self-map family for the actual fixture.

## Reproduction

```bash
HOME=/private/tmp/sage-home sage -python tools/challenge_curve_endomorphism_audit_40bit_n606x.py \
  --out notes/challenge_curve_endomorphism_audit_40bit_n606x.json
```
