# Decision Records: fs-recompute

## DR-RECOMPUTE-001: Semantic Computation Key, Tolerance Role, and Determinism Class Partitioning

- **Status**: Ratified (2026-08-25)
- **Bead**: `frankensim-mkfvu.1`
- **Authors**: FrankenSim Core Team

### Context

FrankenSim's incremental-recompute store (`fs-recompute`) operates as a content-addressed Merkle DAG with first-class slack certificates and a self-policing determinism trip-wire (`StoreError::DeterminismViolation`).

In the legacy v1/v2 schema, every node record stored a single 7-field tuple:
`(op_id, input_hashes, params, code_version_hash, rng_seed, achieved_error, required_tolerance)`.

However, across different scientific and engineering operation families:
1. Some algorithms (e.g. adaptive meshing, Krylov solvers with residual stop criteria, interval Newton roots) legitimately branch and produce different output artifacts when `required_tolerance` changes. Treating `required_tolerance` as merely a query threshold in these algorithms caused either false-positive determinism trip-wire violations or incorrect memo-hits across different discretization fidelities.
2. Other algorithms (e.g. discrete topology evaluations, exact reductions, closed-form formulas) compute a fixed artifact regardless of downstream query thresholds.
3. Stochastic or exploratory algorithms may run in nondeterministic fast modes, which must be explicitly tracked without laundering unverified results into certified authority.

### Decision

1. **`ComputationKey` vs `OutputObservation` Partitioning**:
   - `ComputationKey` captures everything that determines the execution of the kernel: `op_id`, `input_hashes`, `params`, `rng_seed`, `code_version_hash`, `max_iterations`, and `effective_tolerance_bits`.
   - `OutputObservation` captures post-execution metrics: `artifact_hash`, `achieved_error`, `wall_time_s`, `peak_memory_bytes`.

2. **Explicit `ToleranceRole`**:
   - `ToleranceRole::InputParameter`: Tolerance is a direct geometric or grid parameter (e.g. mesh edge length). Affects computation; bound into `ComputationKey`.
   - `ToleranceRole::StoppingCriterion`: Tolerance is an iterative residual threshold. Affects computation; bound into `ComputationKey`.
   - `ToleranceRole::QueryThreshold`: Tolerance is only a post-hoc query threshold. Does NOT affect computation; `effective_tolerance_bits` is 0.
   - `ToleranceRole::None`: Operation does not use tolerances.

3. **Explicit `DeterminismClass`**:
   - `ExactDeterministic`: Bit-identical outputs required across all runs, worker thread counts, and completion permutations on the same ISA.
   - `ToleranceDependentDeterministic`: Bit-identical output guaranteed for fixed `(inputs, params, policy, seed, tolerance)`.
   - `Nondeterministic`: Relaxed / fast mode; cannot mint verified or certified evidence.

4. **Legacy 7-Field Migration**:
   - Legacy records map to `ToleranceRole::StoppingCriterion` by default when `required_tolerance > 0`, preserving conservative semantics without silent reinterpretation.
