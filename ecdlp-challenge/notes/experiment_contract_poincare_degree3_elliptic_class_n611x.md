# Experiment Contract: N611X Degree-Three Elliptic Class in the Reduced PPAV

## Hypothesis

The N611U principal polarization `M` contains a primitive rank-one divisor
class of `M`-degree three. If its unique theta divisor is irreducible and the
standard Hermitian Neron-Severi correspondence is effective in this fixture,
this class yields a degree-three elliptic component of a genus-two theta
curve.

## Candidate

With

```text
u_M = 12+4*pi,
u_D = 13+4*pi,
beta_M = phi_g o u_M,
beta_D = phi_g o u_D,
```

define

```text
M = [[9,beta_M^dagger],[beta_M,345]],
D = [[9,beta_D^dagger],[beta_D,346]].
```

The claimed identities are

```text
det(D)=9*346-2*Norm(13+4*pi)=0,
M.D=3.
```

## Metrics

- exact CM norms, determinant, and intersection;
- primitiveness of the displayed Hermitian class;
- deterministic pointwise replay of the two cross maps;
- bounded search for lower positive intersections as a diagnostic only.

## Positive Control

The coordinate elliptic divisor `diag(0,1)` has `M`-degree nine, confirming
the intersection convention.

## Negative Control

The zero class is not accepted as an elliptic divisor, and diagonal coordinate
classes do not give degree three.

## Success Criterion

The exact displayed class is primitive, rank one, has positive `M`-degree
three, and all map arithmetic replays. This is a class-level positive signal,
not a constructed curve or quotient.

## Falsification Criterion

If any identity fails, do not carry a degree-three cover claim forward.
Passing this does not prove effectivity, irreducibility of the theta divisor,
the quotient elliptic curve, a genus-two equation, a section evaluator,
relations, rank, target descent, cost advantage, or ECDLP speedup.

## Reproduction

```bash
sage -python tools/poincare_degree3_elliptic_class_n611x.py \
  --out notes/poincare_degree3_elliptic_class_n611x.json
sage -python tools/verify_poincare_degree3_elliptic_class_n611x.py \
  --primary notes/poincare_degree3_elliptic_class_n611x.json \
  --out notes/poincare_degree3_elliptic_class_n611x_independent_verifier.json
```
