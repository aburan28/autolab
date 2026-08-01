# Experiment Contract: N608M integral correspondence-lattice degree screen

## Hypothesis

Within the explicit integral correspondence lattice `u_(a,b) = ((a phi_11 + b pi), pi)`, some pair of integers `(a,b)` pulls the H018 polarization back to degree two.

## Null Hypothesis

The degree-two target is not represented in this lattice. This rejects the direct-global-extension route for this declared family only; another correspondence lattice or a geometric modification remains a separate task.

## Exact Model

Use the N608J maps `pi` of degree nine and `phi_11` of degree eleven, with `phi_11_dual pi = 4 + Frob`, `Frob^2 + 5 Frob + 103 = 0`, `z = -3-Frob`, and `zbar = 2+Frob`. The tool derives the pullback degree polynomial, proves the degree-two nonrepresentation modulo nine, and checks the named degree `17`, `23`, and `41` maps against lattice combinations on every rational cover point.

## Controls

- Positive: `(a,b)=(1,0)` reproduces the N608L degree-three correspondence.
- Negative: the degree-two target must fail by the exact congruence proof.

## Success and Falsification

The hypothesis passes only if the polynomial equals two for integral `(a,b)`. A modular obstruction is a scoped negative result for this lattice.

## Reproduction

```bash
sage -python tools/poincare_correspondence_lattice_n608m.py --out notes/poincare_correspondence_lattice_n608m.json
sage -python tools/verify_poincare_correspondence_lattice_n608m.py --primary notes/poincare_correspondence_lattice_n608m.json --out notes/poincare_correspondence_lattice_n608m_independent_verifier.json
```

## Claim Boundary

`NEGATIVE RESULT / DECLARED INTEGRAL CORRESPONDENCE LATTICE ONLY / MODEL-BOUND / TOY-EVIDENCE / OTHER CORRESPONDENCES OR MODIFICATIONS OPEN / NO_ECDLP_CLAIM`.
