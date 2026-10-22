# Restarting feasibility-first differential evolution

`ConstrainedDifferentialEvolutionState::checkpoint_bytes` and
`from_checkpoint_bytes` preserve a drained search across application restarts.
They retain the original objective and maximum scaled constraint violation,
not a penalty objective. Candidate generation and selection remain in the same
serial/parallel DE implementation.

## Retained work

A checkpoint contains the finite search box, initial point, configuration,
spent callback attempts, original lifetime allowance, complete parent
population, partial replacement population, pending failed/unlaunched candidate,
queued successful siblings and earliest best-ranked observation. The candidate
coordinates and scores are stored as exact binary64 words. Restoration performs
no objective calls. A spent allowance stays spent; extending it is explicit.

A successful callback cached immediately before cancellation is restored as
completed work. A failed member keeps its original candidate coordinates and
must pay for another attempt; successful sibling results are not recomputed.
A restored search may continue through either the serial or parallel API.

## Binding the model and selection policy

The caller supplies a 32-byte immutable model identity when saving and the
same expected identity when restoring. This identity must describe the objective,
constraint definitions and row order/scaling, physical inputs, and numerical
implementation. The codec cannot determine whether two arbitrary callbacks
implement the same physical model. Run the same build and numerical mode when
bitwise continuation is required.

Restore also takes the **expected constraint count and exact feasibility
threshold**, separately from the serialized data. Changing either refuses.
Each finite score stores objective and nonnegative maximum scaled violation;
its feasibility flag is recomputed using the admitted threshold. The source
format cannot independently assert that an infeasible score is feasible.

The new format has magic `FSDECF01`, version 1, and checksum domain
`fs-dfo/constrained-differential-evolution-checkpoint-v1`. It includes the
original `fs-rand` semantics version and DE kernel identity. Scalar DE's
`FSDECP01` version-1 format and public methods are unchanged. The two decoders
refuse each other's snapshots.

Byte, dimension and population-shape limits come from the restoring caller,
not from the checkpoint. The common decoder reuses original DE configuration
admission and checks retained counters, ranked best consistency and cached
candidate coordinates. Best consistency uses the original feasibility-first
ordering: the best feasible score need not have the smallest scalar objective,
and equal-objective feasible replacements need not retain its raw violation.

## API

```rust,ignore
let bytes = search.checkpoint_bytes(cx, model_identity, 1 << 20)?;
let limits = DifferentialEvolutionRestoreLimits {
    max_bytes: 1 << 20,
    max_dimension: 16,
    max_population_entries: 4096,
};
let mut restored = ConstrainedDifferentialEvolutionState::from_checkpoint_bytes(
    cx, &bytes, model_identity, expected_rows, expected_tolerance, limits,
)?;
restored.extend_evaluation_budget(larger_lifetime_allowance)?;
let report = restored.advance(cx, next_callback_allowance, &mut evaluate)?;
```

The types are available from
`fs_dfo::differential_evolution::constrained` and its `snapshot` module.
After cancellation, use a fresh noncancelled execution context to capture the
already-drained search. The returned bytes still require application-owned
crash-safe storage. Loading an older checkpoint cannot account for work done
after it was saved, and the checksum is not authentication against a deliberate
rewrite.

This is search-engine restart support. It does not serialize a thermal factory,
mesh or material cards, reconstruct a `PopulationScanStudy`, or restore that
owner's separately maintained physical-evaluation counter. Such owners must
bind and restore their own model/counters before claiming whole-study recovery.
No constrained-stationarity or global-optimality claim follows from a restored
best-so-far candidate.

## Focused checks

Seven Rust regressions cover early-prefix continuation, cancellation retention,
parallel-failure recovery, malformed inputs and policy mismatch, feasible ties,
infeasible target stopping, and the original scalar wire format. The scalar
payload fixture is independently encoded from the pre-refactor field grammar.
Run the new and existing snapshot tests with:

```bash
cargo test -p fs-dfo --lib differential_evolution::constrained::snapshot
cargo test -p fs-dfo --lib differential_evolution::snapshot
```

These commands have not been executed in the delivery environment: Rust and the
repository's DSR/RCH runners are unavailable. Patch applicability and fixture
encoding checks do not establish compilation or runtime correctness.
