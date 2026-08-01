# Experiment Contract: N607F H018 Poincare Fibre Orientation Audit

## Hypothesis

For `L = p1^*Theta^9 tensor p2^*Theta^11 tensor (id,z)^*P` with
`z=-3-pi`, the `p1` fibre over `Q` is represented by `8[O]+[zQ]`.  The
adjoint `zbar=2+pi` belongs in the first component of the polarization map,
not in that fibre index.

## Null Hypothesis

The mod-two kernel graph forces the fibre to use `zbar Q`, invalidating the
N606Y--N607E transport convention.

## Parameters

- Toy curve: `E/F_103: y^2=x^3+x+24`, order `109`.
- CM relation: `pi^2+5pi+103=0`.
- Hermitian form: `[[9,zbar],[z,11]]`.
- Controls: actual fibre action `z=[-4]`; hypothetical adjoint action
  `zbar=[3]` on `E(F_103)`.

## Metrics

- Exact mod-two kernel of `H018`.
- `z` and `zbar` base-field translation shifts.
- Equality of the N606Y `z`-fibre translation law.

## Positive control

The kernel must be the graph `(P,Q)=(pi Q,Q)` and the `z` fibre must use
shift `[84]`.

## Negative control

The adjoint-fibre alternative must require its distinct shift `[106]`, so it
cannot silently equal the N606Y transport model.

## Success criterion

All orientation, kernel, and translation checks pass; the result preserves the
limited scope of N606Y--N607E without asserting a global theta construction.

## Falsification criterion

Any failure requires reclassifying the affected N606Y--N607E receipts and
rebuilding the moving frame under the corrected action.

## Reproduction commands

```bash
sage -python tools/product_kummer_h018_poincare_fibre_orientation_n607f.py --output notes/product_kummer_h018_poincare_fibre_orientation_n607f.json
sage -python tools/verify_product_kummer_h018_poincare_fibre_orientation_n607f.py --primary notes/product_kummer_h018_poincare_fibre_orientation_n607f.json --output notes/product_kummer_h018_poincare_fibre_orientation_n607f_independent_verifier.json
```

## Handoff: N607F fibre/polarization orientation

### Claim or task

Resolve whether the mod-two kernel requires replacing the N606Y fibre action.

### Status

HYPOTHESIS

### Assumptions

- The normalized-Poincare convention already audited by A1 is used.

### Evidence so far

- A1 fixes `H018=[[9,zbar],[z,11]]` for `(id,z)^*P`.
- N606Y uses `8[O]+[zQ]` and shift `[84]`.

### Failure modes

- Confusing the Poincare fibre index with a polarization-map component.

### Next concrete action

Run the primary audit and its independent replay.

### Artifact paths

- `tools/product_kummer_h018_poincare_fibre_orientation_n607f.py`
- `tools/verify_product_kummer_h018_poincare_fibre_orientation_n607f.py`
