# SAFETY: fs-substrate/src/morton/bmi2/mod.rs

> BMI2 PDEP/PEXT Morton interleave capsule. Registered in
> unsafe-capsules.json; enforced by `cargo run -p xtask -- check-unsafe`
> (registration, <300 lines, this file).

## Invariants
All `unsafe` is confined to two `#[target_feature(enable = "bmi2")]`
functions whose bodies call the register-only intrinsics `_pdep_u64` /
`_pext_u64` on integer arguments. No memory is read or written; no pointer
exists in this capsule.

## Aliasing assumptions
N/A — no references or pointers; pure integer-register computation.

## Alignment assumptions
N/A — no memory access.

## Lifetime assumptions
N/A — all values are `Copy` integers.

## Panic behavior
No panics are possible (no asserts, no arithmetic that can overflow-check,
no allocation). Nothing unwinds across the capsule.

## Cancellation behavior
No poll points: each call is a handful of register instructions. Tile-level
cancellation lives in callers (fs-exec discipline).

## Concurrency behavior
Stateless pure functions; trivially `Send`/`Sync`-compatible. Nothing to
order, nothing to race.

## Miri coverage
Miri cannot interpret BMI2 intrinsics. The dispatch layer selects this
capsule only via `is_x86_feature_detected!`, which reports false under Miri,
so Miri runs route to the magic-bits twin. Compensating checks: the
`dispatched_backend_matches_the_magic_reference` battery (100k seeded cases)
and the in-capsule known-answer sweep run on every native x86-64 test run
(the project's Threadripper runner has BMI2).

## Model-checking coverage
N/A (no concurrency).

## Fuzz/property coverage
`morton` module G0 batteries: exhaustive 16-cube bijection, 100k-case seeded
LCG bijection and backend-equivalence sweeps covering coordinate extremes
(0 and 2^21 - 1); the conformance suite repeats the equivalence law.

## Proof obligations discharged by callers
The dispatch table (`morton::fns`) must select this capsule only after
`is_x86_feature_detected!("bmi2")` returns true — the standard runtime-
detection contract, identical to fs-simd's x86 capsule, satisfied at the
single construction site of the `OnceLock` table. Facade callers of
`morton3_encode`/`morton3_decode` carry no obligations.
