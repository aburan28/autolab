# REP-AUXLINE-017 Follow-On Theories

## Theory 1: Bounded-support rational-map triples

### Core mechanism

Search public rational maps whose target-curve fibers have a triple membership
ideal with bounded support before elementary symmetrization, avoiding the
coset power-sum `P_2d` expansion.

### Minimal test

Enumerate low-degree maps on generated prime-field curves, certify the exact
fiber, derive the symmetric triple ideal, and reject any candidate whose static
support is Omega(`B^2`) or whose triple density matches a random base.

### Likely failure

Clearing denominators produces a dense degree-`B` fiber polynomial, as in the
existing auxiliary-isogeny map audit.

### Learning if it fails

It identifies which map invariants force dense triple membership and narrows
the needed structure to a genuine bounded-support fiber algebra.

## Theory 2: Frobenius-conjugate triple quotient

### Core mechanism

Move to an explicitly bounded extension representation where a factor triple
is one Frobenius orbit. Its trace, second trace, and norm are base-field data,
potentially replacing the three independent coset power sums.

### Minimal test

Construct a public extension-field orbit family with an exact base-field
target witness; count relation density, base-field descent cost, and all
extension arithmetic against rho.

### Likely failure

The factor points either do not descend to the target prime-order group or
the Weil-restriction/Semaev system restores the same quadratic state.

### Learning if it fails

It separates a true Frobenius-invariant loophole from merely changing
coordinates around a generic relation search.

### Triage

This is not a fresh implementation branch in the current ledger. REP-N489
already tested fixed trace-zero-section norm fibers through full rank and
descent; REP-N490 found their static S3 eliminant pair-scale-or-worse; and
REP-N493 gives a restricted obstruction to nonconstant base-field algebraic
trace-zero sections on the registered ordinary curves. A successor must name
an explicit non-algebraic or non-invertible correspondence and satisfy those
existing source, rank, and blind-descent gates before another experiment.

## Theory 3: Product-Kummer pencil factor base

### Core mechanism

Complete an actual product-Kummer line-bundle/pencil construction whose fibers
carry correlated factor labels and a target-local incidence relation not
equivalent to ordinary pair or triple enumeration.

### Minimal test

Materialize the pending global gluing/section gate, then certify a public base
scheme, factor map, relation source, rank, and independent descent on toy
prime-field instances.

### Likely failure

The section geometry reduces to reducible divisors or an ordinary Semaev
relation, as the current reducible-divisor audit suggests.

### Learning if it fails

It gives a precise geometric obstruction to using the product-Kummer surface
as a non-generic relation source.
