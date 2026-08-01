# Result: N608V rejects product-linear factor-base compression

## Status

`NEGATIVE RESULT / LINEAR-PROJECTION FACTOR BASE REMAINS BASE-FIELD SCALE /
STRUCTURAL REPLAY / MODEL-BOUND / TOY-EVIDENCE / NONLINEAR REPRESENTATIONS
OPEN / NO_ECDLP_CLAIM`

## Result

N608U supplies complete open source sets for the H018 pencil levels `0`, `1`,
and `17`.  N608V evaluated every projective product-linear map

```text
ell_(a,b)(P,Q) = a P + b Q
```

with directions `(1,b)` for all `b` in `F_103` and `(0,1)`.  The minimum
distinct image supports were:

| Level | Complete sources | Best direction | Minimum image | Fraction of `#E=109` |
|---:|---:|---|---:|---:|
| 0 | 84 | `(1,13)` | 49 | `49/109` |
| 1 | 108 | `(1,7)` | 56 | `56/109` |
| 17 | 144 | `(1,91)` | 68 | `68/109` |

None satisfies the declared one-quarter-order compression gate (`<=27`).

## Interpretation

The complete source-return mechanism does not become a sub-base-field factor
base merely by applying a linear homomorphism from `E x E` to `E`.  This is a
scoped negative for the declared product-linear maps over this toy fixture. It
does not rule out a nonlinear target-bearing correspondence, a map to a
Jacobian or Prym, a cover-side representation, or a different curve family.

## Next Positive Question

Can a nonlinear correspondence preserve enough of the complete level geometry
to make its target image genuinely sparse, while still providing a relation
law, rank, individual-log descent, and fully charged comparison to Pollard
rho?
