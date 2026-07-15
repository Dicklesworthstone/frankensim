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
- `geom::marching_cubes` retains its historical Three.js triangle/normal wire
  layout but samples the selected analytic demo field into `fs-viz::Grid3` and
  delegates all polygonization, indexing, budget, and winding semantics to the
  shared native marching-tetrahedra implementation.
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
5. TrussPath optimality uses the same `fs-truss-e2e` promotion gate as the
   native campaign. The transcribed sparse arrays and PDHG iterates are bound
   into the private `fs-truss` receipt; only a matching outward certificate
   serializes rank `2`, a verified flag, and finite optimum endpoints. A hard
   error or numerical unavailability serializes Estimated/no-bound fields.
   TrussPath load-path promotion likewise calls the exact shared
   `fs-truss-e2e::certify_load_path` implementation: its six wire fields carry
   rank/flag/outward path endpoints and two exact u32 words for the 64-bit
   replay golden. The golden is a drift sentinel only; promotion authority
   remains the private exact receipt and lower-layer BLAKE3 identities.
6. FRAME wire version 2 applies that same gate to its normalized layout LP,
   then outward-divides verified endpoints by the physical yield stress. The
   four claim fields are appended to the layout block immediately before the
   existing sizing offset, so all earlier layout fields retain their positions.
7. The browser isosurface is a serialization adapter, not a second polygonizer:
   each indexed `fs-viz` triangle is expanded to three positions plus three
   analytic gradient normals in the documented 18-value block.

## Error model

No structured error type in v0. Browser-facing functions use bounded
fallback outputs for invalid or failed kernel calls. Native helper
functions may still use ordinary Rust assertions from their upstream
crates when called outside the clamped public surface.
`marching_cubes` maps non-finite input, sampling/allocation refusal, or a
surface exceeding 60,000 triangles to the one-value `[0]` sentinel; it never
serializes a silently truncated/open prefix.

## Determinism class

Deterministic for fixed inputs on one target/ISA, subject to the
determinism contracts of the underlying crates. Cross-browser and
cross-ISA bit identity is not claimed for floating-point visual demos.

## Cancellation behavior

Most browser work is bounded by input clamps and fixed iteration caps.
TrussPath optimum and load-path certificate construction additionally run under
a deterministic `fs-exec::Cx` and poll cancellation through their cold proof
stages. That scoped context uses `CancelGate::new_clock_free`, so constructing
it never reads the unsupported `wasm32-unknown-unknown` platform time source;
its private sentinel request marker is omitted from timestamp accessors and
latency reports. The browser surface does not yet expose an external
cancellation handle.

## Unsafe boundary

None in this crate. `unsafe_code = "forbid"` is set locally.

## Feature flags

None. WASM-only dependencies are target-gated under `cfg(target_arch =
"wasm32")`.

## Conformance tests

Native unit tests in the nested workspace exercise root demos, campaign
defaults, geometry/PDE/deep modules, flagship headline/determinism cases, and
the exact clock-free TrussPath certificate context. The native-host TrussPath
transcription test compares both serialized claim ranks/flags/outward endpoints
and reconstructs the load-path replay golden from its two wire words for exact
comparison against the native campaign; it is not browser execution or
cross-target bit-identity evidence.
The geometry module additionally checks the shared isosurface's triangle-count
wire length, finite unit normals, exact replay, and non-finite sentinel.
Current verification is native cargo test/clippy of the nested workspace plus
any wasm32 build lane provided by DSR or site automation. The wasm32 browser
surface itself remains a build/smoke lane rather than a browser-E2E test suite.

## No-claim boundaries

- Not a packaged public simulator API.
- Not a general certification API; campaign functions surface summaries and
  visualizable traces from lower crates. TrussPath's serialized optimum and
  material-volume path intervals are narrow exceptions and carry only their
  lower-layer receipts' declared graph/LP claims.
- The shared promotion gate gives native and browser code the same claim-strength
  rules, but cross-target endpoint bit identity remains unclaimed until a retained
  browser runner or WASM golden exists.
- No browser performance claim without wasm32 benchmark artifacts.
- No guarantee that every native crate feature is available in WASM.
