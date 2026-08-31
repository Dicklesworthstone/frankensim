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
| `tesla_coil_step(resonant_freq_khz, input_kv, spark_gap_mm, q_factor) -> String` | wasm + native | Admits one bounded positive finite input tuple, executes one deterministic Tesla-coil LC step through `fs_flux::lc::step_tesla_coil`, and returns either `{"ok":{...}}` or a typed `{"refusal":{...}}` envelope. |

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
4. **Fail-closed admission.** Non-finite, non-positive, or above-cap inputs
   refuse before the owning kernel runs. A non-finite or non-positive result
   also refuses rather than being rendered as JSON `null` or a misleading
   zero. Each boundary is exercised by native tests.

## Error model

The boundary returns stable refusal codes with a human-readable message and
ranked repairs. `non-finite-input` covers NaN and infinity;
`input-outside-domain` covers non-positive values and the documented safety
caps; `non-finite-output` catches a broken downstream result;
`output-outside-domain` catches zero or negative output. The admitted
maximums are 10,000 kHz, 1,000 kV, 1,000 mm spark gap, and Q of 1,000,000.

## Determinism class

Bit-identical replay on the same ISA for identical inputs. Cross-ISA
identity is not claimed (floating-point rendering follows the platform's
shortest-roundtrip display).

## Cancellation behavior

Synchronous single-step entry; there is nothing to cancel. A JS caller can
discard the returned string if the result is no longer needed.

## Unsafe boundary

None: `#![forbid(unsafe_code)]` is declared in Cargo.toml lints and the
crate contains no unsafe block.

## Feature flags

None. The wasm32 binding is selected by target architecture
(`cfg_attr(target_arch = "wasm32", wasm_bindgen)`), not by a feature.

## Conformance tests

Native unit tests exercise a valid result, non-finite and non-positive input
refusals, lower and upper admission boundaries, underflowed output refusal,
and ranked repairs. The fs-flux kernel tests and dependency-policy/contract
gates remain complementary evidence.

## No-claim boundaries

This crate is a presentation seam. It makes no claim that the embedded
LC-step model is validated, that its lumped constants match any measured
coil, or that the rendered numbers carry uncertainty bounds. Physics
claims belong to `fs-flux`'s own contract and evidence artifacts — never
to this renderer.
