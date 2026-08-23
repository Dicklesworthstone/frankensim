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
9. NeuroShape preserves the field offsets `[0..=21]` of its historical 24-value
   header while separating enclosed-component existence from exact component
   counting. Its closed interval boundary frame can certify only a
   component-count lower bound of one. Header slot `[17]` is therefore always
   the JSON-safe `-1` exact-count unknown sentinel in this tranche; `[20]` is
   `0` for unknown evidence or `1` for a certified enclosed-component lower
   bound, and `[21]` carries that lower bound. The positive-definite
   finite-difference Hessian at the origin is curvature corroboration only;
   without a zero-gradient certificate it is not a critical-point or minimum
   theorem. It and the sampled contour crossings cannot promote these fields to
   an exact count.
10. NeuroShape wire version `NEUROSHAPE_SCHEMA_VERSION = 3` (header slot `[22]`,
   header length 30) publishes a no-tunnel step whose authority is an INTERVAL
   sign margin. Slot `[5]` is `fs_rep_neural::derive_safe_step`'s
   downward-rounded `magnitude_lower_bound / L`, where the margin at `[23]` is
   the inward endpoint of the degenerate IBP enclosure at the origin — a
   certified lower bound on `|f(0)|` — and `[24]` is the derivation's typed
   status (`1` sign-separated, `0` no finite sign margin, `2` malformed
   enclosure, `3` invalid Lipschitz bound). `[5]` is `0` whenever `[24] != 1`.
   The nominal forward pass stays at `[4]` as display data only: version `1`
   published `|f(0)_nominal|/L` at `[5]` as a proven step, which overstated the
   true `|f(0)|/L` because the forward pass's own evaluation error was
   unaccounted for. Version `1` also carried
   `NEUROSHAPE_COMPONENT_EVIDENCE_SCHEMA_VERSION` in `[22]`; that value moved to
   `[25]`. Version `3` keeps every version-2 slot at its position and extends
   the header from 27 to 30 so the report's typed zero-set localization outcome
   (`fs_neuroshape_e2e::SurfaceLocalization`, bead frankensim-o33vo) crosses
   the ABI: `[26]` is its stable status code (`1` localized, `2` valid-empty,
   `3` invalid input, `4` unrepresentable, `5` resource refused,
   `6` cancelled, `7` allocation refused, `8` internal fault; `0` stays the
   reserved no-claim code), `[27]` is the producing stage (`0` success/no
   single stage, `1` grid construction, `2` isocontour extraction). Legacy
   slots `[6]`, `[7]`, and `[15]` are derived, non-authoritative views of the
   typed outcome: every non-finite legacy slot is folded to `NaN` by the
   JSON-safe `fon` encoder, so a valid-empty payload carries `[6] = NaN`
   (the folded `+inf`) and `[7] = +0`, while the native report fields keep
   `+inf`; the typed `[26]` status, never a NaN sentinel, distinguishes
   valid-empty from every refusal. `[28]` carries
   `NEUROSHAPE_LOCALIZATION_SCHEMA_VERSION = 1`, which freezes those code
   tables and gates interpretation of `[26]`/`[27]`; `[29]` remains a
   reserved zero. A consumer
   that gated on `[22] == 2` refuses this payload rather than silently
   ignoring the typed localization record; version-aware consumers must refuse
   an unrecognized `[22]` before interpreting any slot, an unrecognized `[25]`
   before interpreting `[16]`, `[17]`, `[20]`, or `[21]`, and an unrecognized
   `[28]`, `[26]`, or `[27]` code before trusting any crossing view. The export
   also runs the campaign through `try_run_campaign`, so an inadmissible net or
   geometry serializes an empty vector instead of trapping or publishing a
   partial header.
11. FlutterCert slot `[8]` is the flag of `witness_decay_rate_color`, which names
   exactly one quantity: the LARGEST eigenvalue real part of `A(witness_mu)`.
   Its endpoints are `fs_flutter_e2e::spectral_abscissa_interval`'s
   outward-rounded ones, not a round-to-nearest `−1 + √(μ−1)`, and they are not
   an enclosure of the operator's spectrum — for `μ > 1` the second eigenvalue's
   real part lies strictly below them. The endpoints are not serialized in this
   layout, so the wire format is unchanged.
12. FlowCert is the `fs-flowcert-e2e` campaign, not a browser-local
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
bound. Both cases pin payload schema version 2 in slot `[22]`, the
component-evidence version in `[25]`, and zero in the remaining reserved slot
`[26]`. A decoder-shaped conformance fixture accepts only the exact version-2
bit pattern and refuses the legacy version-1 value, zero, future version 3,
fractional, non-finite, and truncated headers before reading any evidence, and
refuses an unrecognized component-evidence version in `[25]`. A safe-step test
pins slots `[5]`/`[23]`/`[24]` bit-for-bit against the native campaign's
`SafeStepDerivation` and requires the published step to be STRICTLY below the
nominal `|origin_value|/L` that version 1 published, while still
under-estimating the sampled nearest-surface radius; a companion case pins the
`NoFiniteSignMargin` status of a zero-straddling field to the zero radius and
the `0` wire code. A FlutterCert test drives the reaching input
`fluttercert(1.2, 1.9, 8)` (whose witness has `μ > 1`, unlike the graceful
default sweep) and pins slot `[8]` to a `witness_decay_rate_color` whose
endpoints are `spectral_abscissa_interval`'s, strictly outward of the
round-to-nearest abscissa, with the operator's second eigenvalue outside the
claim.
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
- NeuroShape's safe step is a no-tunnel bound at the ORIGIN derived from that
  point's enclosure and the global Lipschitz upper bound. It is not a Euclidean
  distance to the zero set, not a bound at any other point, and slot `[6]`'s
  sampled `nearest_surface_radius` is localization evidence, never the
  certificate. The derivation is arithmetic conditional on its inputs: it
  carries no field identity, issuer, or portable receipt on the wire.
- The shared promotion gate gives native and browser code the same claim-strength
  rules, but cross-target endpoint bit identity remains unclaimed until a retained
  browser runner or WASM golden exists.
- No browser performance claim without wasm32 benchmark artifacts.
- No guarantee that every native crate feature is available in WASM.
