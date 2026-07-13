# CONTRACT: fs-wasm

> Status: PARTIAL. This crate is a standalone nested workspace that exposes
> browser/WASM demos over existing FrankenSim kernels; it is not part of the
> native Cargo workspace build graph.

## Purpose and layer

Browser surface over FrankenSim numerical leaves and end-to-end campaign
crates. Layer: **L6 HELM / interface surface**. The crate compiles as an
`rlib` for native smoke checks and as a `cdylib` for
`wasm32-unknown-unknown`.

## Public types and semantics

- Root functions in `src/lib.rs` expose deterministic numerical demos:
  sparse Poisson/heat kernels, Chebyshev/Orr-Sommerfeld probes, QMC,
  interval/Taylor arithmetic, forward-mode AD, FFT, eigensolves, and
  randomized NLA summaries.
- `deep`, `geom`, `pde`, `dynamics`, `certified`, and `campaigns`
  modules re-export broader browser demos over upper-stack crates.
- The `#[wasm_bindgen]` JavaScript boundary is compiled only for
  `wasm32`; native builds exercise the same pure Rust functions.

## Invariants

1. Browser demos call the real crate kernels, not mocks or rewritten JS
   numerics.
2. Public demo inputs are clamped or bounded before allocating or
   iterating so browser calls cannot request unbounded work.
3. Fallible demo paths return `NaN`, empty vectors, or bounded fallback
   values rather than trapping across the WASM boundary. The vessel CVaR
   surface maps canonical `fs-robust` validation errors to `NaN` here instead
   of reintroducing a panic-only risk implementation.
4. The nested workspace isolates browser-only dependencies from the
   native workspace dependency policy.

## Error model

No structured error type in v0. Browser-facing functions use bounded
fallback outputs for invalid or failed kernel calls. Native helper
functions may still use ordinary Rust assertions from their upstream
crates when called outside the clamped public surface.

## Determinism class

Deterministic for fixed inputs on one target/ISA, subject to the
determinism contracts of the underlying crates. Cross-browser and
cross-ISA bit identity is not claimed for floating-point visual demos.

## Cancellation behavior

None in v0. Work is bounded by input clamps and fixed iteration caps.

## Unsafe boundary

None in this crate. `unsafe_code = "forbid"` is set locally.

## Feature flags

None. WASM-only dependencies are target-gated under `cfg(target_arch =
"wasm32")`.

## Conformance tests

Native unit tests in the nested workspace exercise root demos, campaign
defaults, geometry/PDE/deep modules, and flagship headline/determinism cases.
Current verification is native cargo test/clippy of the nested workspace plus
any wasm32 build lane provided by DSR or site automation. The wasm32 browser
surface itself remains a build/smoke lane rather than a browser-E2E test suite.

## No-claim boundaries

- Not a packaged public simulator API.
- Not a certification surface; campaign functions surface summaries and
  visualizable traces from lower crates.
- No browser performance claim without wasm32 benchmark artifacts.
- No guarantee that every native crate feature is available in WASM.
