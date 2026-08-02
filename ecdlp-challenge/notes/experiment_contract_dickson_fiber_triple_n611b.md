# Experiment Contract: N611B Dickson-Fiber Triple Screen

## Hypothesis

A Dickson-polynomial factor fiber `D_d(x,a)=c` can retain the sparse
recurrence of a monomial coset while giving a lower-support exact
three-factor membership condition and nonrandom zero-sum triple density.

## Null hypothesis

Reducing `D_d(T,a)-c` modulo the generic cubic with roots `x1,x2,x3` has
pair-scale-or-worse support, or its curve-point triple density is matched by
negation-stable random bases of the same size.

## Parameters

- toy curve: `E/F_101: y^2=x^3+x+32`, of prime order `101`;
- maps: Dickson recurrence `D_0=2`, `D_1=x`,
  `D_k=xD_(k-1)-aD_(k-2)`;
- degrees: `10,20,25`; parameters `a=1,2`; fiber constants `c=0,1`;
- factor base: all affine curve points satisfying `D_d(x(P),a)=c`;
- controls: 16 deterministic negation-stable random x-pair bases of equal
  cardinality.

## Exact source condition

For `U(T)=T^3-e1*T^2+e2*T-e3`, a triple of distinct affine x-coordinates is
in the declared fiber exactly when

```text
D_d(T,a)-c == 0 mod U(T).
```

The three remainder coefficients are the pre-solver membership equations.

## Metrics

- factor-base size and exact membership replay;
- generic remainder monomial support and degree;
- complete distinct zero-sum triple count;
- matched-control mean and maximum triple count;
- support-to-pair-state ratio and candidate enumeration count.

## Positive control

The Dickson recurrence must agree with direct finite-field recurrence values
for every point in the toy curve scan, and every accepted triple must sum to
the identity.

## Negative controls

- Matched random bases preserve x-pair and negation structure.
- Incomplete or non-paired fibers are rejected rather than silently padded.
- A support or density observation is not a source inverse, factor-log rank,
  target descent, or ECDLP speedup.

## Admission criterion

A follow-on source solver requires one declared fiber with exact controls,
generic remainder support strictly below the unordered x-pair state, and at
least a fourfold density lift above every matched control. It must then pass a
separate rank, individual-descent, and fully charged rho comparison.

## Reproduction command

```bash
sage -python tools/dickson_fiber_triple_probe_n611b.py \
  --out notes/dickson_fiber_triple_probe_n611b.json
```

## Claim boundary

`HYPOTHESIS / DICKSON_FIBER_TRIPLE_PREFLIGHT / MODEL-BOUND / TOY-EVIDENCE /
NO_ECDLP_CLAIM`.
