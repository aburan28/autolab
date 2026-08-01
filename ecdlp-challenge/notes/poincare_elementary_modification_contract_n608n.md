# Experiment Contract: N608N degree-corrected H018 evaluator target

## Hypothesis

The N608J correspondence pulls the normalized symmetric H018 line bundle back
to `O_Eprime(3Oprime)`. Its degree-three pushforward therefore has the
canonical elementary modification

```text
0 -> pi_* O(2Oprime) -> pi_* O(3Oprime) -> k(O) -> 0.
```

This supplies a degree-correct target for extending the N608K evaluator; a
global extension with a length-one cokernel at `O` would recover the desired
degree-two bundle without ad hoc fibre weights.

## Null Hypothesis

The pullback has a nontrivial rational two-torsion translate, or the divisor
inclusion does not push forward to the stated elementary modification. Either
failure rejects this target construction on the registered fixture.

## Assumptions

- The H018 line bundle is the normalized symmetric, `F_103`-rational line
  bundle defined by its recorded Hermitian polarization.
- `pi:Eprime -> E` and `phi_11:Eprime -> E` are the N608J maps.

## Controls

- Positive: the N608L induced polarization action is `[3]` on every one of
  the 109 rational cover points.
- Negative: any rational nonzero two-torsion point would leave the line-class
  normalization ambiguous.

## Success Criterion

The degree-three symmetric pullback has origin class `O(3Oprime)` and the
finite pushforward of `O(2Oprime) -> O(3Oprime)` is the stated exact sequence.

## Next Construction

Construct normalized-Poincare chart matrices for the evaluator into
`pi_*O(3Oprime)`, then prove its cokernel is exactly the origin skyscraper.
No section, relation, rank, descent, or ECDLP speedup claim is permitted here.

## Reproduction

```bash
sage -python tools/poincare_elementary_modification_n608n.py --out notes/poincare_elementary_modification_n608n.json
sage -python tools/verify_poincare_elementary_modification_n608n.py --primary notes/poincare_elementary_modification_n608n.json --out notes/poincare_elementary_modification_n608n_independent_verifier.json
```

## Claim Boundary

`RESTRICTED THEOREM / DEGREE-CORRECT TARGET IDENTIFIED / STANDARD-FACT-BOUND /
MODEL-BOUND / TOY-EVIDENCE / GLOBAL EVALUATOR EXTENSION OPEN / NO_ECDLP_CLAIM`.
