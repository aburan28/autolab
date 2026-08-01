# Result: N608L Reject Direct Global Extension of the Degree-11 Evaluator

## Exact Endomorphism Calculation

For the N608J maps `phi:Eprime -> E` and `pi:Eprime -> E`, set

```text
A = phi_dual * pi.
```

The producer compares the x-coordinate rational maps exactly and uses a
generator sign check to establish

```text
A = [4] + Frob,
A_dual = [-1] - Frob,
Frob^2 + 5 Frob + 103 = 0.
```

With `z=-3-Frob` and `zbar=2+Frob`, the H018 induced polarization is

```text
phi_dual (9 phi + zbar pi) + pi_dual (z phi + 11 pi)
  = [198] + (2+Frob)(4+Frob) + (-3-Frob)(-1-Frob)
  = [3].
```

The independent verifier also evaluates the induced map on every point of the
109-point cover group and obtains `[3]` throughout.

## Consequence

The pullback of the H018 line bundle along

```text
R |-> (phi(R), pi(R))
```

has degree three.  N606U's target is `O_Eprime(2Oprime)`, of degree two.
Therefore the N608K rational local evaluator cannot extend directly to the
target as a global line-bundle morphism.  This is a line-bundle degree mismatch,
not a failed finite-field rank test.

## Boundary

This rejects only this direct degree-eleven correspondence as a global
extension route. It does not rule out another correspondence with pullback
degree two, a rigorously constructed elementary modification, a different
H018 frame, source recovery, relation rank, target descent, or an ECDLP
improvement.
