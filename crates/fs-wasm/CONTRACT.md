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
   of reintroducing a panic-only risk implementation. A fallback must not be a
   value the payload's own documented meaning would read as a successful
   result: `dynamics::ga_motor_orbit` therefore folds a refused PGA sandwich
   product to `(NaN, NaN, NaN)` rather than republishing the untransformed seed
   point, which would be bit-indistinguishable from a genuine identity image.
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
8. GrammarForge uses the exact shared `fs-grammar-e2e` simplification assessment
   and summary, not a browser-local interpretation of ShapeProg certificates.
   Its soundness bit additionally requires one assessment per serialized elite;
   a sound strict subset cannot promote the browser headline.
   Its 32-value header preserves the historical fields at `[0..=20]` and adds:
   local radius threshold `[21]`, maximum admitted outward finite-sample check `[22]`,
   typed aggregate status code `[23]`, observed assessment count `[24]`, and
   refusal/non-finite/finite-negative/evidence-refusal/structural-empty/
   conservative-check-exceedance/threshold-mismatch counts `[25..=31]`. The
   `0.02` offset admitted at local threshold `0.03` serializes the sound global
   certificate `0.04`; status remains independently explicit.
9. NeuroShape preserves its 24-value header and field offsets while separating
   enclosed-component existence from exact component counting. Its closed
   interval boundary frame can certify only a component-count lower bound of
   one. Header slot `[17]` is therefore always the JSON-safe `-1` exact-count
   unknown sentinel in this tranche; `[20]` is `0` for unknown evidence or `1`
   for a certified enclosed-component lower bound, and `[21]` carries that
   lower bound. Reserved slot `[22]` now carries
   `NEUROSHAPE_COMPONENT_EVIDENCE_SCHEMA_VERSION = 1`; version-aware consumers
   must refuse unsupported values before interpreting slots `[16]`, `[17]`,
   `[20]`, or `[21]`. Slot `[23]` remains reserved zero. Legacy consumers that
   ignored both reserved slots are not thereby made version-aware and must be
   migrated explicitly. The positive-definite finite-difference Hessian at the
   origin is curvature corroboration only; without a zero-gradient certificate
   it is not a critical-point or minimum theorem. It and the sampled contour
   crossings cannot promote these fields to an exact count.
10. FlowCert is the `fs-flowcert-e2e` campaign, not a browser-local
    transcription of it. `campaigns::flowcert` calls `run_campaign`, so every
    published headline — the MAP-Elites atlas statistics, `all_accurate`, and
    the `map_color_rank` in slot `[7]` — is the report's own, and each point's
    `accurate` bit is the campaign's `converged && profile_error <= tol`. A
    point whose chunked steady-state march exhausts the step cap is published
    as unresolved (`converged = 0`, `accurate = 0`), never as an accurate or
    cleanly-inaccurate point, and it drops the map's color out of `Verified`.
    The payload is `FLOWCERT_SCHEMA_VERSION = 2` in header slot `[8]`; version
    1 was the fixed-step-budget layout with 8-value header and point blocks and
    no way to express an unresolved point. Version 2 keeps every version-1
    field at its old position within the header and within each point block,
    appends `converged`/`steps_run` to each block, and adds `all_converged` at
    `[9]`. Consumers must refuse an unrecognized `[8]` before reading any other
    slot. Only the two spotlight velocity profiles are recomputed in this
    crate, because the native `OperatingPoint` does not carry a profile; that
    re-march is the same chunked loop and mints no claim.

## Error model

No structured error type in v0. Browser-facing functions use bounded
fallback outputs for invalid or failed kernel calls. Native helper
functions may still use ordinary Rust assertions from their upstream
crates when called outside the clamped public surface.
`marching_cubes` maps non-finite input, sampling/allocation refusal, or a
surface exceeding 60,000 triangles to the one-value `[0]` sentinel; it never
serializes a silently truncated/open prefix.
GrammarForge clamps infinite local radius thresholds to the documented public
range. A NaN threshold remains a `NaN` numeric sentinel while separately
serializing `SimplifierRefused` and its count. Non-finite and finite-negative
nominal certificates have distinct typed counts. A transactional core rollback
therefore cannot be laundered as a zero-error successful simplification.

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
GrammarForge tests compare every native/WASM simplification-summary value and
status bit, assert the exact `0.02/0.03 -> 0.04` envelope, and verify that a NaN
threshold serializes typed refusals and cannot promote the headline.
NeuroShape tests assert that the default closed-frame certificate serializes a
lower bound of one but an unknown exact count, while an unenclosed case retains
the same wire shape and claims neither an exact zero nor a positive lower
bound. Both cases pin schema version 1 in slot `[22]` and retain zero in the
remaining reserved slot `[23]`. A decoder-shaped conformance fixture accepts
only the exact version-1 bit pattern and refuses legacy zero, future version 2,
fractional, non-finite, and truncated headers before reading topology evidence.
FlowCert tests pin the whole payload field-for-field against
`fs_flowcert_e2e::run_campaign` (including that slot `[7]` is the rank of the
report's own `credibility_color` and slot `[8]` the schema version), assert that
no point claims `accurate` while `converged` is 0, name the default-budget
point whose profile error is inside tolerance but whose march never reaches
steady state and require it to serialize as unresolved, check that the minimum
step budget publishes nine unresolved points and a non-`Verified` map, and pin
the spotlight re-march's `converged`/`steps_run` to the campaign's.
A dynamics test builds a degenerate (non-versor) motor whose sandwich product is
refused and requires the `(NaN, NaN, NaN)` fold rather than the seed point.
Current verification is native cargo test/clippy of the nested workspace plus
any wasm32 build lane provided by DSR or site automation. The wasm32 browser
surface itself remains a build/smoke lane rather than a browser-E2E test suite.

## No-claim boundaries

- Not a packaged public simulator API.
- Not a general certification API; campaign functions surface summaries and
  visualizable traces from lower crates. TrussPath's serialized optimum and
  material-volume path intervals are narrow exceptions and carry only their
  lower-layer receipts' declared graph/LP claims.
- GrammarForge's radius threshold controls local rewrite admission only. It is
  not a global error budget, and no browser consumer may infer
  `max_certified_error <= radius_threshold`. The sampled discrepancy field is a
  conservative outward finite-grid admission check, not a continuum proof or a
  downward-rounded lower bound; the compositional ShapeProg certificate carries
  the declared global algebraic authority.
- NeuroShape's certified-inside witness plus closed positive boundary frame
  proves that at least one enclosed component exists. It does not prove there
  are no additional components inside or outside the frame. The finite-
  difference origin-Hessian check and sampled contour crossings are
  corroborating evidence only. No zero-gradient certificate is present, so the
  Hessian check does not establish a critical point or minimum, and the browser
  surface makes no exact component-count or full homeomorphism claim.
- The shared promotion gate gives native and browser code the same claim-strength
  rules, but cross-target endpoint bit identity remains unclaimed until a retained
  browser runner or WASM golden exists.
- No browser performance claim without wasm32 benchmark artifacts.
- No guarantee that every native crate feature is available in WASM.
