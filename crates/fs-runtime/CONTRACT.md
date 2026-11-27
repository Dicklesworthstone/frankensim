# fs-runtime Contract

## Purpose and layer
`fs-runtime` (layer `L6`, HELM) is the shared runtime crate over `fs-system` for native Rust and browser WebAssembly front ends (`frankensim-rc-1109-root-3gc12.11.9`). It contains no physics or integrator of its own; every step executes `fs_system::System::step` directly.

## Public types and semantics
- `AssemblyRuntime`: owns a sealed kernel `System`, pre-allocated `ArenaPool`, lock-free `AtomicU64` interactive input cells (initialized from each `InputDescriptor::default_value` in `from_system`), pre-resolved probe readers, generated `SnapshotDescriptor`, and once-delivered refusal state.
- `InputBinding`, `ProbeBinding`, `CustomProbeFn`: typed input and probe bindings for `AssemblyRuntime::from_system`.
- `LiveSource`, `LiveInputSpec`, `RUNTIME_INPUT_FAMILY`: lock-free interactive power source (`Component + PortTangent + StageResidual`), declarative live input channel, and power-port lowering hook.
- `SnapshotDescriptor` (`SNAPSHOT_DESCRIPTOR_SCHEMA = "frankensim.runtime.snapshot_descriptor.v1"`), `ProbeDescriptor`, `InputDescriptor`, `RuntimeRefusal`, `inputs_to_json`, `MAX_PACKS`, `UNIT_DIMS`, `dims_to_unit` (re-exported from `fs-fsim`): canonical descriptors, refusals, and JSON encoders.

## Invariants
- On the success path of `set_input`, `step`, `read_probes`, and `render_f32`, the runtime wrapper adds no heap allocation (refusal paths format diagnostics, and any allocation inside a component belongs to the sealed `System`); `tests/no_alloc.rs` measures 0 heap allocations for a 1 kg `Mass` + `ControlledSource` system after one warm-up step.
- Every interactive power input enters the kernel through `LiveSource` (`Component + PortTangent + StageResidual`), preserving `BlockSolveMode::MonolithicNewtonKrylov` parity with static lowering and closing first-law energy conservation in `System::ledger()`.
- On refusal or non-finite/overflowing probe output during `render_f32`, the runtime halts, writes exact `0.0` silence for remaining frames, and returns the refusal once from `take_refusal()`.
- Canonical JSON encoders escape strings deterministically and emit `null` for non-finite `default_value`.

## Error model
Every invalid binding, step refusal, or buffer mismatch returns `RuntimeRefusal { code, stage, message }` with canonical JSON encoding; nothing is silently clamped or ignored.

## Determinism class
`ExecMode::Deterministic` with explicit caller `seed` and fixed `StreamKey`; identical inputs and step sizes produce bit-identical native probe outputs.

## Cancellation behavior
Each runtime owns a private `CancelGate::new_clock_free()` passed into `Cx` on every `System::step`; no R1c API requests cancellation, so every step ends in completion or a refusal.

## Unsafe boundary
`#![deny(unsafe_code)]` across the crate. Zero `unsafe` blocks in library code (only the test-only `GlobalAlloc` counter in `tests/no_alloc.rs`).

## Feature flags
None. All dependencies are pure Franken workspace crates; `wasm32-unknown-unknown` compilation is checked.

## Conformance tests
- `descriptor::tests::canonical_json_writers_pin_exact_bytes`: byte-pin JSON checks.
- `live_source::tests`: power-port role/dimension guard, `LiveSource` vs static `Force` bit-for-bit & `BlockSolveMode::MonolithicNewtonKrylov` parity (including no-source graph mode check), and first-law work booking.
- `tests/no_alloc.rs`: atomic `GlobalAlloc` proof of 0 heap allocations on a 1 kg `Mass` + `ControlledSource` system after warm-up, `from_system` default-cell initialization check, and `render_f32` `f32`-overflow silence halt.

## No-claim boundaries
`fs-runtime` adds no physical model, integrator, material provenance, or audio post-processing. Scientific claims remain owned by the admitted assembly's cards and `fs-system` components.

## Planned (R1d-R3)
- Cut R1d: `AssemblyRuntime::from_draft`, `AssemblyRuntime::load`, observer probe resolution, and `tests/runtime_core.rs`.
- Cut R3: `fs-runtime-wasm` (`wasm-bindgen` wrapper + cross-target native-vs-`wasm32` parity tests).
