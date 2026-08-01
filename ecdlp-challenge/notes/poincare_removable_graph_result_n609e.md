# N609E Result: removable P=4Q extension of the residual source inverse

## Handoff: restore the Cauchy chart's removable graph point

### Claim or task

Determine whether the common cleared graph `P=4Q` can be restored by the
regular Cauchy limit, and whether doing so makes the degree-nine residual
inverse complete on the extended chart.

### Status

OBSERVATION / REMOVABLE P=4Q CHART EXTENSION OF DEGREE-NINE RESIDUAL INVERSE /
STRUCTURAL REPLAY / MODEL-BOUND / TOY-EVIDENCE / GLOBAL COMPACTIFICATION OPEN /
NO_ECDLP_CLAIM

### Assumptions

- The N608S/N608U Cauchy-coordinate evaluator and N608Z clearing convention
  are the declared local chart.
- `Q` is nonzero on the prime-order `F_103` fixture, so `y(-4Q)` is nonzero.
- The tangent limit of the Cauchy coordinate is the applicable regularization
  at `P=4Q`.

### Evidence so far

The primary computation passes seven gates and the structural replay passes
six.  Direct raw Cauchy evaluation fails at `P=4Q` for all 108 nonzero base
points, as expected from its `0/0` representation.  Its tangent limit

```text
-(3*x(-4Q)^2 + a) / (2*y(-4Q))
```

agrees exactly with the degree-nine residual's retained graph intersections.
The unsquared cleared equation reconstructs each retained graph point.

| Level | Recovered open chart | Restored graph points | Extended sources |
|---:|---:|---:|---:|
| 0 | 84 | 8 | 92 |
| 1 | 108 | 0 | 108 |
| 17 | 144 | 2 | 146 |

At every level, degree-nine factoring plus the regularized graph values equals
an exhaustive evaluation of the extended chart.

### Failure modes

- This is a local chart extension. It does not yet prove a compact global
  residual divisor or determine its second projection, divisor class, or
  normalization.
- The source-count increase is not a factor-base reduction, relation law,
  rank result, target descent, charged cost advantage, or ECDLP improvement.
- The true pole point `P=-4Q` remains excluded; only its inverse `P=4Q` was
  removable.

### Next concrete action

Use the extended chart as local evidence while constructing a global residual
divisor, then compare its Q-boundary and second projection with the formal
H018 degree requirements before testing smoothness or any relation mechanism.

### Artifact paths

- `notes/poincare_removable_graph_contract_n609e.md`
- `tools/poincare_removable_graph_n609e.py`
- `notes/poincare_removable_graph_n609e.json`
- `tools/verify_poincare_removable_graph_n609e.py`
- `notes/poincare_removable_graph_n609e_structural_replay.json`
