# fs-runtime Contract

## Purpose and layer
`fs-runtime` (layer `L6`, HELM) is the shared runtime crate over `fs-system` for native Rust and browser WebAssembly front ends (`frankensim-rc-1109-root-3gc12.11.9`). In Cut R1a, it ships the canonical snapshot/input/refusal descriptor types, dimension-to-unit helper, and zero-dependency canonical JSON encoders.

## Public types and semantics
- `SnapshotDescriptor` (`SNAPSHOT_DESCRIPTOR_SCHEMA = "frankensim.runtime.snapshot_descriptor.v1"`), `ProbeDescriptor`, and `InputDescriptor`: assembly-derived descriptors carrying probe offsets, lengths, shapes, SI dimensions, and defaults.
- `RuntimeRefusal`: typed refusal `{ code, stage, message }` with canonical JSON encoding (`RuntimeRefusal::to_json`).
- `inputs_to_json`: canonical JSON array encoder for `&[InputDescriptor]`.
- `dims_to_unit`: canonical SI unit label for common `fs_qty::Dims` (unified with `fs-cli`'s `UNIT_DIMS` in `fs-fsim` in Cut R2b).

## Invariants
- `SnapshotDescriptor::to_json`, `inputs_to_json`, and `RuntimeRefusal::to_json` emit deterministic canonical JSON with escaped string literals and emit `null` for non-finite `default_value` entries.
- `SnapshotDescriptor::probe(name)` resolves descriptors by declared name without hard-coded indices.

## Error model
Refusals are represented by `RuntimeRefusal { code, stage, message }` (implementing `Display` and `std::error::Error`) and serialize to canonical JSON without panicking or silently clamping values.

## Determinism class
Pure deterministic data formatting and lookup: identical descriptor values always produce identical UTF-8 JSON bytes natively and in `wasm32-unknown-unknown`.

## Cancellation behavior
Synchronous, bounded O(n) formatting over caller-provided slices; no unbounded loops or background work.

## Unsafe boundary
`#![deny(unsafe_code)]` across the crate. Zero `unsafe` blocks.

## Feature flags
None. All dependencies are pure Franken workspace crates compatible with `wasm32-unknown-unknown`.

## Conformance tests
- `descriptor::tests::canonical_json_writers_pin_exact_bytes`: unit test pinning exact canonical JSON bytes for `SnapshotDescriptor`, `inputs_to_json` (including non-finite `NaN` -> `null`), `RuntimeRefusal`, string escaping, and `dims_to_unit`.

## No-claim boundaries
`fs-runtime` adds no physical model, integrator, material provenance, or audio post-processing. Scientific claims remain owned by the admitted assembly's cards and `fs-system` components.

## Planned (R1b-R3)
- Cut R1b: `AssemblyRuntime::from_system`, `InputBinding`, `ProbeBinding`, zero-allocation `set_input`/`step`/`read_probes`/`render_f32`, and `tests/no_alloc.rs`.
- Cut R2a/R2b: `fs-fsim` split and unified `.fsim` codec + `UNIT_DIMS`/`dims_to_unit`.
- Cut R1c/R1d: `LiveSource` (`StageResidual`), `LiveInputSpec`, `from_draft`/`load`, and `tests/runtime_core.rs`.
- Cut R3: `fs-runtime-wasm` (`wasm-bindgen` wrapper + Node/WASM parity tests).
