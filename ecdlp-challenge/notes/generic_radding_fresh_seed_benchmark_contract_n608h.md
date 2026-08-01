# Experiment Contract: N608H fresh-seed generic r-adding benchmark

## Hypothesis

The current opaque-token r-adding walk has an aggregate 40-bit cost at or
below the Pollard-rho reference when evaluated on independently derived,
nondefault instances and independently fixed token encodings.

## Null hypothesis

Across the registered fresh seeds, the aggregate counted group operations
exceed the aggregate rho reference. Individual low scores are then treated as
collision variance, not as a reusable generic improvement.

## Parameters

- benchmark: `ECDLP_BITS=40`;
- instance seeds: `0x12345679`, `0x1234567a`, `0x1234567b`, `0x1234567c`;
- token seeds: `0xa6080001`, `0xa6080002`, `0xa6080003`, `0xa6080004`;
- solver: the current local `src/solver/mod.rs`, recorded by SHA-256 and left
  unmodified by this experiment;
- baseline: trusted oracle `rho_reference` for each generated instance.

## Metrics

- oracle correctness;
- charged group operations;
- per-instance rho reference and ratio;
- aggregate operations divided by aggregate rho reference;
- min, median, and max per-instance ratios.

## Positive control

Each row must return a correct scalar under a nondefault instance seed and a
fixed independent token encoding.

## Negative control

The default development fingerprint is excluded: each registered instance must
have an order different from the locally hardcoded default order.

## Success criterion

All rows are correct and aggregate cost is below the aggregate rho reference.
This would still be a finite generic constant-factor observation, not a
sub-square-root result.

## Falsification criterion

Any incorrect row, default-order match, or aggregate ratio above one rejects
the claimed expected-cost improvement for this seed panel.

## Reproduction

```bash
ECDLP_BITS=40 ECDLP_SEED=0x12345679 ECDLP_TOKEN_SEED=0xa6080001 ./benchmark.sh --note "N608H clean panel row 1"
ECDLP_BITS=40 ECDLP_SEED=0x1234567a ECDLP_TOKEN_SEED=0xa6080002 ./benchmark.sh --note "N608H clean panel row 2"
ECDLP_BITS=40 ECDLP_SEED=0x1234567b ECDLP_TOKEN_SEED=0xa6080003 ./benchmark.sh --note "N608H clean panel row 3"
ECDLP_BITS=40 ECDLP_SEED=0x1234567c ECDLP_TOKEN_SEED=0xa6080004 ./benchmark.sh --note "N608H clean panel row 4"
python3 tools/generic_radding_fresh_seed_benchmark_n608h.py \
  --out notes/generic_radding_fresh_seed_benchmark_n608h.json
python3 tools/verify_generic_radding_fresh_seed_benchmark_n608h.py \
  --primary notes/generic_radding_fresh_seed_benchmark_n608h.json \
  --out notes/generic_radding_fresh_seed_benchmark_n608h_independent_verifier.json
```

## Claim boundary

`OBSERVATION / FRESH-SEED GENERIC COST MEASUREMENT / MODEL-BOUND /
TOY-EVIDENCE / NO SUB-SQUARE-ROOT CLAIM / NO ECDLP BREAK CLAIM`.
