# fs-tesla-wasm — CONTRACT

Layer: L6 (boundary crate). Standalone `[workspace]` following the
`fs-wasm` / `fs-flyer-wasm` pattern: browser builds never depend on the
state of unrelated in-progress native crates.

## Purpose and layer

The single boundary between the Tesla Coil museum page (browser JS) and
the simulation plane (`fs-flux` lumped LC step), compiled to
`wasm32-unknown-unknown`. It exposes one step function over the
`fs_flux::lc` kernel and renders its result as one flat JSON record.
Native builds compile the same pure function for tests; `wasm-bindgen` is
confined to wasm32 via a target-specific dependency so native builds stay
dependency-clean.

## Public types and semantics

| Entry | Kind | Contract |
|---|---|---|
| `tesla_coil_step(resonant_freq_khz, input_kv, spark_gap_mm, q_factor) -> String` | wasm + native | Executes one deterministic Tesla-coil LC step through `fs_flux::lc::step_tesla_coil` and returns `{"ok":{...}}` with the result fields rendered as decimal numbers. The JSON shape is part of the page contract: field names are stable, values are plain JSON numbers. |

## Invariants

1. **Shape stability.** The emitted record keeps exactly the fields of
   `TeslaCoilResult` under their struct names inside an `"ok"` object; no
   field is dropped, renamed, or reordered semantically by the renderer.
2. **Dependency purity.** Runtime dependencies are the workspace crate
   `fs-flux` only (Decalogue P1); serialization is explicit formatting,
   not a serde edge.
3. **Deterministic content.** Identical inputs produce byte-identical
   output on the same platform; no wall clock, randomness, host identity,
   or allocation-order dependence reaches the string.

## Error model

v0 performs no admission: non-finite or out-of-domain inputs propagate
into the arithmetic and surface as `null` in JSON (JSON has no NaN). Typed
refusal envelopes (`{"refusal": ...}` per the fs-flyer-wasm pattern) are
the v1 upgrade path and are NOT yet implemented; callers must not rely on
refusal semantics from this crate.

## Determinism class

Bit-identical replay on the same ISA for identical inputs. Cross-ISA
identity is not claimed (floating-point rendering follows the platform's
shortest-roundtrip display).

## Cancellation behavior

Synchronous single-step entry; there is nothing to cancel. A JS caller
abandons the promise-like return by simply not using it.

## Unsafe boundary

None: `#![forbid(unsafe_code)]` is declared in Cargo.toml lints and the
crate contains no unsafe block.

## Feature flags

None. The wasm32 binding is selected by target architecture
(`cfg_attr(target_arch = "wasm32", wasm_bindgen)`), not by a feature.

## Conformance tests

Native unit coverage lives with the owning lane's battery plan; today the
crate is covered indirectly by the fs-flux kernel tests plus the
dependency-policy and contract gates (`cargo run -p xtask -- check-deps`,
`check-contracts`). A dedicated round-trip golden for the rendered string
is the next planned slice.

## No-claim boundaries

This crate is a presentation seam. It makes no claim that the embedded
LC-step model is validated, that its lumped constants match any measured
coil, or that the rendered numbers carry uncertainty bounds. Physics
claims belong to `fs-flux`'s own contract and evidence artifacts — never
to this renderer.
