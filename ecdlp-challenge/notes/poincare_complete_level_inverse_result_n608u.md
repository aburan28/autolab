# Result: N608U constructs a complete open-level source generator

## Status

`OBSERVATION / COMPLETE OPEN-LEVEL SOURCE GENERATOR / STRUCTURAL REPLAY /
MODEL-BOUND / TOY-EVIDENCE / RELATION MECHANISM OPEN / NO_ECDLP_CLAIM`

## Result

For the registered `F_103` H018 fixture, fix a nonzero base point `Q` and a
pencil value `t`.  Clearing the moving Cauchy denominator in the N608S
equation `lambda(P,Q)=t` gives an equation

```text
A_Q,t(x) + B_Q(x)y = 0.
```

Eliminating `y` with the curve equation produces

```text
A_Q,t(x)^2 - B_Q(x)^2 (x^3+x+24) = 0.
```

The eliminated polynomial has degree exactly ten in every declared trial.
Factoring one such polynomial for each of the 108 nonzero base points, then
recovering `y`, reproduces every enumerated point in the declared open level.

| Level | Algebraic inverse sources | Exhaustive open sources | Root-factor calls |
|---:|---:|---:|---:|
| 0 | 84 | 84 | 108 |
| 1 | 108 | 108 | 108 |
| 17 | 144 | 144 | 108 |

No exceptional `A=B=0` branch occurred.  A mutation of the `x` coefficient in
the Cauchy-cleared equation changes the recovered set, so the stated
construction is not accepted merely because the test compares a formula with
itself.

## Interpretation

N608U removes the N608T graph-density obstruction for this source-generation
subproblem: it recovers complete open pencil levels rather than the sparse
graph slice.  It does not yet supply a relation law, a factor base, factor-base
logarithms, matrix rank, target descent, or an operation-counted comparison to
Pollard rho.  It is therefore not evidence of an ECDLP speedup.

## Next Action

Test whether complete level sources admit a projection with a factor base
strictly smaller than the base-field scale, while charging source generation
and preserving any target-bearing relation needed by the projected attack.
