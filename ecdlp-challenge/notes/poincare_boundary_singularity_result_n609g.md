# N609G Result: no rational singularity candidate on screened boundary loci

## Handoff: removable graph and Q-origin boundary screen

### Claim or task

Extend N609F's exact local derivative screen to the rational removable graph
and finite Q-origin boundary points of the declared H018 pencil.

### Status

OBSERVATION / NO RATIONAL SINGULARITY CANDIDATE ON SCREENED H018 BOUNDARY LOCI
/ STRUCTURAL REPLAY / MODEL-BOUND / TOY-EVIDENCE / GLOBAL SMOOTHNESS OPEN /
NO_ECDLP_CLAIM

### Assumptions

- N609E's Cauchy tangent regularization is the correct local extension at
  `P=4Q`.
- N608R/N609C's regularized `Q=O` section model is the applicable boundary
  frame.
- The true Cauchy pole and non-rational points are outside this rational screen.

### Evidence so far

The primary computation passes five gates and the structural replay passes
five checks.  All restored graph sources reconstruct their level under both
local expansions and have nonzero derivatives in both directions:

| Level | Restored graph sources | Both derivatives zero |
|---:|---:|---:|
| 0 | 8 | 0 |
| 1 | 0 | 0 |
| 17 | 2 | 0 |

At `Q=O`, the finite rational boundary source counts are `2,0,0` at levels
`0,1,17`; every tested derivative is nonzero.  Independently, every member of
the selected boundary pencil retains one simple `P=O` zero as an `L(9O)`
section.

### Failure modes

- The true pole `P=-4Q`, non-rational geometric points, global gluing, and
  normalization remain untested.
- This local evidence does not establish smoothness, genus, a Jacobian,
  factor bases, relations, rank, target descent, charged cost, or ECDLP speed.

### Next concrete action

Construct a global compact residual section and analyze the true pole and
non-rational loci before promoting any local screens to a smooth curve or
Jacobian claim.

### Artifact paths

- `notes/poincare_boundary_singularity_contract_n609g.md`
- `tools/poincare_boundary_singularity_probe_n609g.py`
- `notes/poincare_boundary_singularity_n609g.json`
- `tools/verify_poincare_boundary_singularity_n609g.py`
- `notes/poincare_boundary_singularity_n609g_structural_replay.json`
