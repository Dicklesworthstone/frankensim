# fs-goddard-wasm — CONTRACT

Layer: L6 (boundary crate). Standalone `[workspace]` following the
`fs-wasm` / `fs-flyer-wasm` pattern: browser builds never depend on the
state of unrelated in-progress native crates.

## Purpose and layer

The single boundary between the Goddard Rocket museum page (browser JS)
and the simulation plane (`fs-mbd` rigid-body/rocket physics), compiled to
`wasm32-unknown-unknown`. It exposes one step function over the
`fs_mbd::goddard` kernel and renders its result as one flat JSON record.
Native builds compile the same pure function for tests; `wasm-bindgen` is
confined to wasm32 via a target-specific dependency so native builds stay
dependency-clean.

## Public types and semantics

| Entry | Kind | Contract |
|---|---|---|
| `goddard_rocket_step(chamber_pressure_psi, fuel_flow_kg_per_sec, throat_area_cm2, expansion_ratio) -> String` | wasm + native | Executes one deterministic Goddard rocket step through `fs_mbd::goddard::step_goddard_rocket` and returns `{"ok":{...}}` with the result fields rendered as decimal numbers. The JSON shape is part of the page contract: field names are stable, values are plain JSON numbers. |

## Invariants

1. **Shape stability.** The emitted record keeps exactly the fields of
   `GoddardResult` under their struct names inside an `"ok"` object; no
   field is dropped, renamed, or reordered semantically by the renderer.
2. **Dependency purity.** Runtime dependencies are the workspace crate
   `fs-mbd` only (Decalogue P1); serialization is explicit formatting, not
   a serde edge.
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
crate is covered indirectly by the fs-mbd kernel tests plus the
dependency-policy and contract gates (`cargo run -p xtask -- check-deps`,
`check-contracts`). A dedicated round-trip golden for the rendered string
is the next planned slice.

## No-claim boundaries

This crate is a presentation seam. It makes no claim that the embedded
Goddard model is validated, that its constants are experimentally
identified, or that the rendered numbers carry uncertainty bounds. The
physics claims, if any ever, belong to `fs-mbd`'s own contract and
evidence artifacts — never to this renderer.
