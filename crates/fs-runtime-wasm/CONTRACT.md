# fs-runtime-wasm Contract

## Purpose and layer
`fs-runtime-wasm` (layer `L6`, HELM; standalone nested workspace) exposes `fs_runtime::AssemblyRuntime` across `wasm-bindgen` for browser, AudioWorklet, and Node hosts while keeping `wasm-bindgen` out of the native workspace dependency graph.

## Public types and semantics
- `WasmRuntime::load(fsim_bytes, pack_bytes, pack_lengths, seed_override)`: decodes canonical `.fsim` source and length-framed `FSMCDPK\0`/`FSINTPK\0` card packs through `fs-fsim` and lowers through `fs-runtime`.
- `inputs()`, `set_input()`, `set_input_by_name()`, `probes()`, `snapshot_len()`, `read_probes()`, `step()`, `render_f32()`, `ledger()`, `take_refusal()`, `is_halted()`, `time_s()`, `steps()`.

## Invariants
Pure Rust (`#![forbid(unsafe_code)]`); delegates all physics, lowering, energy accounting, and admission to `fs-runtime`, `fs-fsim`, and `fs-system`. `slice_packs` validates pack counts (`<= MAX_PACKS = 128`), byte ranges, and exact consumption before forwarding to `AssemblyRuntime::load`.

## Error model
Refusals cross JS as canonical `RuntimeRefusal::to_json()` strings (`Err(String)` on fallible calls, `take_refusal()` once after `render_f32` halts with exact silence).

## Determinism class
Deterministic (`ExecMode::Deterministic`); `tests/runtime.test.mjs` verifies 0-ULP IEEE-754 parity against native `fs-runtime` on `forced_oscillator.fsim` (50 stepped frames, mid-run `set_input_by_name`, and `render_f32`).

## Cancellation behavior
Delegates to `AssemblyRuntime`'s private `CancelGate::new_clock_free()`; no `WasmRuntime` API requests cancellation, so every step ends in completion or a JSON refusal.

## Unsafe boundary
`#![forbid(unsafe_code)]`. Zero `unsafe` blocks.

## Feature flags
None. `wasm-bindgen` is confined to `cfg(target_arch = "wasm32")`.

## Conformance tests
Native unit test in `src/lib.rs` and `node --test tests/runtime.test.mjs` verifying dynamic descriptor lookup, bit-for-bit parity on `forced_oscillator.fsim`, `render_f32`, ledger closure, and halt/refusal once-delivery.

## No-claim boundaries
Runtime transport only; scientific validity belongs to the admitted assembly and cards.
