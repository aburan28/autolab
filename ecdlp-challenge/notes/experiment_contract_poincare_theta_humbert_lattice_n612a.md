# Experiment Contract: N612A Full Hom Lattice and Humbert Irreducibility Gate

## Hypothesis

The N611U principal class `M` is geometrically irreducible as a theta class on
`E x E_g`. The test must use the full Neron-Severi lattice, not only the
subfamily of classes already written as `phi_g o u`.

## Formal Model

Assume the ordinary CM endomorphism identity

```text
End(E) = Z[pi],  pi^2+5*pi+103=0,
```

and the standard dual-isogeny identity `phi_g^dagger phi_g=[2]`. For any
`h:E -> E_g`, set `v=phi_g^dagger h`. If `2h=phi_g v` annihilates `E[2]`,
the explicit `E[2]` action decides whether `v` is even. This tests
`Hom(E,E_g)=phi_g Z[pi]` in the stated CM model.

On the resulting integral lattice write a class as `(a,b,r,s)` with
`u=r+s*pi` and

```text
q_A(a,b,r,s) = a*b - 2*Norm(u),
M = (9,345,12,4),
q_M(D) = (D.M)^2 - 4*q_A(D).
```

Kani's refined Humbert criterion says that a principal theta class is
irreducible exactly when `q_M` does not represent one.

## Success Criterion

The full Hom lattice equality is supported by the `E[2]` divisibility gate,
and an integral unimodular reduction of `q_M` has no nonzero vector of value
at most eight while retaining an exact value-nine witness. This proves the
criterion's `q_M(D) != 1` condition under the stated standard assumptions.

## Controls

- nonzero `1`, `Frob`, and `1+Frob` mod-two CM multipliers must not be
  annihilated by `phi_g` on all of `E[2]`;
- the N611X rank-one class must remain a value-nine positive control;
- the quotient reduction must have a primitive null vector exactly equal to
  the theta class.

## Falsification Criterion

A missing full-Hom reduction, any nonzero multiplier annihilating all
`E[2]`, a nonunimodular quotient basis, or a Humbert vector of value one
invalidates the irreducibility gate. Passing does not give an equation,
section evaluator, factor base, relation, rank, descent, cost advantage, or
ECDLP break.

## Reproduction

```bash
sage -python tools/poincare_theta_humbert_lattice_n612a.py \
  --out notes/poincare_theta_humbert_lattice_n612a.json
sage -python tools/verify_poincare_theta_humbert_lattice_n612a.py \
  --primary notes/poincare_theta_humbert_lattice_n612a.json \
  --out notes/poincare_theta_humbert_lattice_n612a_independent_verifier.json
```
