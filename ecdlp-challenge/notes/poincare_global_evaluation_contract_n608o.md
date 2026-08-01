# Experiment Contract: N608O global universal-evaluation elementary transform

## Hypothesis

Let `F` be the normalized H018 divisor-model pushforward with fibre
`H^0(E,O(8O+[-4]Q))`, and let

```text
h:Eprime -> E x E,  R |-> (phi_11(R), pi(R)).
```

Because N608N identifies `h^*L_H018` with `O(3Oprime)`, universal evaluation
defines a global morphism

```text
ev_h: F -> pi_*O(3Oprime).
```

Its determinant should vanish only at the base origin, giving a length-one
cokernel and an explicit degree-two elementary transform isomorphic to `F`.

## Null Hypothesis

The graph evaluation is not generically rank nine, or its determinant data are
incompatible with a single origin zero. Either result rejects this global
elementary-transform route on the registered fixture.

## Exact Inputs

- N608I: `det(F)=O(2O)`.
- N608N: `h^*L_H018=O(3Oprime)` under the stated normalization and
  `det(pi_*O(3Oprime))=O(3O)`.
- N608J: all 108 nonzero rational base fibres have rank-nine graph evaluation.

## Controls

- Positive: recompute all 108 rank-nine graph evaluations over `F_103^6`.
- Negative: the `pi` control has rank one on every corresponding deck fibre.

## Success Criterion

The universal-evaluation map is globally defined by pullback and evaluation,
is generically full rank, and its determinant is a nonzero section of `O(O)`.
The resulting cokernel is `k(O)`.

## Boundary

This does not identify the resulting elementary transform with the specific
canonical subbundle `pi_*O(2Oprime)`; that requires a quotient-direction audit
at the origin. It gives no surface sections, source inverse, relation matrix,
descent, cost comparison, or ECDLP speedup.

## Reproduction

```bash
sage -python tools/poincare_global_evaluation_n608o.py --out notes/poincare_global_evaluation_n608o.json
sage -python tools/verify_poincare_global_evaluation_n608o.py --primary notes/poincare_global_evaluation_n608o.json --out notes/poincare_global_evaluation_n608o_independent_verifier.json
```
