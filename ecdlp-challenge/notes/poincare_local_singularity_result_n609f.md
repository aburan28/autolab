# N609F Result: no rational singularity candidate on the open residual chart

## Handoff: formal two-direction local screen

### Claim or task

Screen every rational source of the declared open H018 residual levels for a
simultaneous zero of the exact P- and Q-direction local derivatives.

### Status

OBSERVATION / NO RATIONAL SINGULARITY CANDIDATE ON DECLARED OPEN H018 CHART /
STRUCTURAL REPLAY / MODEL-BOUND / TOY-EVIDENCE / GLOBAL SMOOTHNESS OPEN /
NO_ECDLP_CLAIM

### Assumptions

- The N608S raw evaluator and the N608R graph/cover formal model describe the
  same regular open chart.
- The cover formal parameter maps separably to the Q direction at the declared
  nonzero base points.
- The true pole and boundary loci remain outside this screen.

### Evidence so far

The primary computation passes five gates and the structural replay passes six.
The exact Q-direction Laurent construction reproduces the raw evaluator at all
11,448 regular rational `(P,Q)` pairs.  At every open source, both local
expansions reproduce its level value.

| Level | Open sources | P derivative zero | Q derivative zero | Both zero |
|---:|---:|---:|---:|---:|
| 0 | 84 | 0 | 0 | 0 |
| 1 | 108 | 0 | 0 | 0 |
| 17 | 144 | 2 | 4 | 0 |

Thus the observed separate projection ramification at level `17` does not
produce a rational singularity candidate on this open chart.

### Failure modes

- This does not cover `P=O`, the true Cauchy pole `P=-4Q`, the separately
  regularized `P=4Q` locus, the `Q=O` boundary, or non-rational geometric
  points.
- No global compact equation, smoothness theorem, genus, normalization,
  Jacobian, factor base, relation law, rank, descent, cost, or ECDLP advantage
  follows from this screen.

### Next concrete action

Extend the exact derivative screen across the removable graph and the finite
`Q=O` boundary, then construct a global compact residual section before any
smoothness or Jacobian claim.

### Artifact paths

- `notes/poincare_local_singularity_contract_n609f.md`
- `tools/poincare_local_singularity_probe_n609f.py`
- `notes/poincare_local_singularity_n609f.json`
- `tools/verify_poincare_local_singularity_n609f.py`
- `notes/poincare_local_singularity_n609f_structural_replay.json`
