# CONTRACT: frankensim-apple

## Purpose and layer

Small C ABI adapter between the native Apple application and `fs-wasm`. It is
an L6 interface surface and deliberately lives in an isolated workspace.

## Public types and semantics

`frankensim_apple_run` executes one bounded catalog entry and replaces the
calling thread's result packet. `frankensim_apple_result_len` and
`frankensim_apple_result_value` copy that packet without transferring Rust
ownership. `frankensim_apple_last_error` reports the calling thread's status.

## Invariants

1. All computations call existing `fs-wasm` functions.
2. No allocator-owned pointer crosses the ABI.
3. Every successful result has the documented six-value packet header.
4. Panics are contained and never unwind through Swift.
5. Work sizes are bounded before invoking a kernel.

## Error model

Error `0` is success, `1` is an unknown experiment, `2` is a caught panic, and
`3` is an invalid result. Failed runs publish an empty packet.

## Determinism class

The adapter adds no nondeterminism. Individual kernel contracts remain
authoritative.

## Cancellation behavior

Runs are currently bounded but not externally cancellable. The Swift actor
prevents overlapping calls and can discard a result after task cancellation.

## Unsafe boundary

There are no unsafe blocks or raw-memory operations. Rust 2024 nevertheless
requires the `unsafe(no_mangle)` attribute for stable C symbol names; that
linkage declaration is the only reason the crate permits the `unsafe_code`
lint.

## Feature flags

None.

## Conformance tests

Unit tests pin structured packet shapes, unknown-id refusal, thread-local
replacement, and a bounded real-kernel smoke run for every public catalog ID.

## No-claim boundaries

This adapter does not upgrade scientific evidence, authenticate results, expose
the full FrankenSim API, or make a general end-user simulator claim.
