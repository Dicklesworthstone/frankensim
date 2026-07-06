# SAFETY: <crate>/<module path>

> One file per registered unsafe capsule (patch Rev P). Every section is
> mandatory; "N/A" requires a one-line justification. The capsule registry
> (unsafe-capsules.json) records this file's existence; xtask check-unsafe
> fails CI if the capsule drifts (unregistered unsafe, >300 lines, missing
> SAFETY.md).

## Invariants
What must ALWAYS hold inside this capsule, stated precisely.

## Aliasing assumptions
Which pointers/references may alias; how exclusivity is guaranteed.

## Alignment assumptions
Required alignments and who establishes them (e.g., 128-byte arena policy).

## Lifetime assumptions
What outlives what, and which lifetimes are erased/reconstructed.

## Panic behavior
What happens if code inside the capsule panics; unwind safety across the
facade (project policy: panics are captured at tile scopes and converted to
structured diagnostics — docs/CONVENTIONS.md).

## Cancellation behavior
Interaction with asupersync scopes: what state exists at cancellation points;
why cancellation cannot tear capsule state.

## Concurrency behavior
Thread-safety claims (Send/Sync justifications), memory-ordering choices.

## Miri coverage
Which tests run under Miri; known Miri limitations for this capsule
(e.g., inline asm/intrinsics) and the compensating checks.

## Model-checking coverage
For concurrency capsules: which interleavings are enumerated, by what harness
(G4). "N/A (single-threaded)" is acceptable with justification.

## Fuzz/property coverage
Fuzz targets on the safe facade; property tests; corpus location.

## Proof obligations discharged by callers
What the SAFE facade's callers must uphold (should be: nothing — if callers
carry obligations, the facade is not safe; document why it still qualifies).
