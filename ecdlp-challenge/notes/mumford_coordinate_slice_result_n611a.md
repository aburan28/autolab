# N611A Result: Non-Split Mumford-Coordinate Factor-Base Screen

## Result

`NEGATIVE RESULT / NONSPLIT MUMFORD-COORDINATE FACTOR-BASE OCCUPANCY /
MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM`

Two target-independent factor bases on the N608D genus-two fixture were
defined directly by the linear coefficient of the reduced degree-two Mumford
polynomial `u(t)=t^2+a*t+b`: `a=0` and `a in {0,1}`.  Signed triple sums were
tested for exact equality to the public embedded 103-point elliptic subgroup,
with eight matched random degree-two Mumford bases for each size.

The initial run has full target coverage and row rank for both coordinate
slices, but so does every matched-control maximum.  Its smaller slice contains
four direct embedded-target factors, compared with a control maximum of two,
so the observed hit counts are not interpretable as a coordinate advantage.

The target-free companion removes every embedded-target divisor from candidate
and control bases.  It preserves full coverage and rank on both sides while
the coordinate-slice hit counts fall below the matched-control means:

| `a` values | Factor base | Coordinate hits | Control mean | Control maximum |
|---|---:|---:|---:|---:|
| `{0}` | 52 | 1510 | 1619.0 | 1712 |
| `{0,1}` | 102 | 12288 | 12387.5 | 12696 |

All accepted sources replay exactly.  The largest target-free slice has
`1,435,820` logical signed triples, a `20,910`-entry pair index, and `20,808`
target-minus-third queries; its measured toy collection time is about 3.14
seconds.  Those are collector metrics, not charged cryptanalytic complexity.

## Interpretation

This rejects these two low-degree Mumford `u`-coefficient factor-base
geometries as a nonrandom source of triple-relation coverage or rank in this
fixture.  The result does not rule out a different Mumford surface, a
non-enumerative packet inverse, or another non-split internal-Jacobian route.
It supplies no factor-log projection back to `E`, blind target descent,
relation collection analysis, sparse linear algebra, or sub-rho ECDLP
algorithm.

## Artifacts

- Contract: `notes/experiment_contract_mumford_coordinate_slice_n611a.md`
- Probe: `tools/mumford_coordinate_slice_probe_n611a.py`
- Primary receipt: `notes/mumford_coordinate_slice_probe_n611a.json`
- Target-free red-team receipt:
  `notes/mumford_coordinate_slice_probe_n611a_target_free.json`

## Next Concrete Action

Define a Mumford factor surface whose membership yields a non-enumerative
packet inverse and specify an explicit same-form factor-log projection back to
the original elliptic target before repeating occupancy collection.
