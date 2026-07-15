# fs-scenario — CONTRACT

The boundary-condition and load-case ALGEBRA (plan patch Rev D): a
`Scenario` is a typed value answering "what is being done to the domain?"
— with dimensional analysis on every value, provenance (seed + canonical
IR), and admission-time validity checks that catch the class of mistakes
no solver can fix.

Ambition tags: typed BCs/frames/signals/combos [S]; seeded ensembles
(Dryden, Kanai–Tajimi, Carreau bands) [S]; canonical IR [S].

## Purpose and layer

Layer **L3** (FLUX support). Runtime deps: `std`, fs-blake3, fs-qty,
fs-rand, fs-cheb, fs-exec, fs-ga, fs-math. The Design Ledger stores scenarios
as canonical IR artifacts; that integration lives ABOVE L3 (fs-ledger is L6) and is
exercised here through a dev-dependency in conformance tests only.
Consumers: fs-solid, fs-flux, fs-lbm, fs-uq, fs-regime, the milestone
flagships.

## Public types and semantics

- `signal::TimeSignal` — `Constant`, `Ramp` (finite strict interval, clamped;
  stable convex interpolation; the vessel tilt
  `(ramp 0deg 65deg 3s)`), `Table` (strictly increasing times + declared
  `Interp` contract: Linear/Hold, clamped ends), `Chebfun` (fs-cheb
  function object). Every signal knows its `Dims`.
- `frame::FrameTree` — frames with `Fixed`, `Rotating`, and `Tilt`
  motions; poses are fs-ga MOTORS composed down the parent chain
  (`world_pose`), cycle/dangling-parent checked. Rotation about an
  off-origin axis is `T(c)·R·T(−c)`.
- `bc` — `BoundaryCondition { region, physics, kind, value, compatibility,
  frame }`; `expectation(physics, kind)` is the dimensional contract
  table (velocity for flow Dirichlet, kg/s for mass-flow inlets, Pa for
  pressure/traction, K / W/m² / W/(m²K) for thermal, m for elastic
  Dirichlet; no-value kinds; everything else structurally Unsupported).
  Flux-carrying inlets MUST declare `Compat::Incompressible`. A declared
  total mass flow is uniform or time-varying kg/s; a spatial profile is refused
  until a geometry-bound layer can retain and certify its surface integral.
- `ensemble::StochasticEnsemble` — seeded generators: Dryden gust PSD,
  Kanai–Tajimi ground-acceleration PSD (spectral representation with
  Gaussian coefficients — a genuine Gaussian process), Carreau parameter
  bands. `realize(member)` is a deterministic function of the complete
  canonical ensemble spec (including duration and dt), member, and versioned
  stream/synthesis semantics; its Philox `StreamKey` is replayable bitwise.
  The typed receipt that exposes those algorithm-version fields is tracked by
  `frankensim-sj31i.39`. `realize` uses the
  documented `DEFAULT_REALIZATION_BUDGET`; `realize_with_budget` admits work
  under explicit `max_samples` and deterministic `max_work` limits before any
  output allocation. Both entry points independently validate public structs.
  Duration/dt are non-sampling placeholders for band models, but remain finite
  positive time quantities because they are retained in canonical IR.
  `SpectrumModel::try_psd` likewise refuses invalid frequency, dimensions,
  parameter domains, and non-finite derived values instead of leaking raw
  NaN/∞ through a direct query. Its one-sided spectral domain is ω ≥ 0;
  non-spectral Carreau bands are a typed refusal rather than a fictitious zero
  PSD.
- `scenario::Scenario` — root value: frames, base BCs, `LoadCase`s,
  factored `Combination`s (`1.2D + 1.6L`), ensembles, `ContactLaw`s
  (Frictionless/Coulomb/Tied), explicit `Environment` (gravity, ambient
  temperature/pressure — REQUIRED constructor argument, never silently
  defaulted). `validate()` returns structured `Violation { code, what,
  fix }` values under `DEFAULT_VALIDATION_BUDGET`; a resource refusal becomes a
  non-green structured violation. `validation_plan` exposes the checked
  collection/signal/checkpoint/identity/work shape for later ledger binding,
  while `validate_with_budget` accepts an explicit `ValidationBudget` and
  fs-exec `Cx`, returning typed `ValidationError` refusals. Scenario, frame,
  case, combination, ensemble, and region
  identities are exact UTF-8 strings: non-ASCII is admitted, empty identities
  and duplicates within each identity role are refused, and no normalization
  is implicit. Repeated combination terms are refused rather than silently
  summed. Contact pairs are unordered; a repeated equal model is a duplicate
  and a repeated different model is a conflict, both with declaration-row
  provenance.
- `ir::write_ir`/`ir::parse_ir` — canonical byte-stable, explicitly versioned
  s-expression encoding. v2 writes six-base `[m, kg, s, K, A, mol]` vectors;
  explicit v1 and historical unversioned five-vector forms decode by appending
  `mol = 0` and return a `DecodedScenario` carrying the immutable
  `DimensionCrosswalkReceipt`. The receipt BLAKE3-binds the exact supplied v1
  bytes to the exact canonical v2 re-emission and exposes a verifier; explicit
  and implicit v1 spellings therefore retain distinct source identities.
  Floats use shortest-round-trip form and u64 seeds remain exact integers.
  Physically irrelevant signed zero is canonicalized to `0`, matching the
  scenario types' semantic `PartialEq` rather than creating two content
  identities for equal values.
  Strings use exactly the writer's quote and backslash escapes; every other
  backslash sequence is rejected so distinct authority bytes cannot alias one
  decoded identity through an undocumented escape rule.
  `parse_ir` applies `DEFAULT_IR_PARSE_BUDGET`; `parse_ir_with_budget` exposes
  explicit byte, recursive-depth, total-node, atom/string-byte, and per-list
  child limits. Every limit is checked before recursive descent or syntactic
  tree growth; caller-selected depth may tighten but cannot exceed the
  recursive implementation's hard safety ceiling. Non-finite wire numbers and
  invalid Chebyshev constructor inputs are structured parse refusals, never
  panics. The byte/node limits also bound decoder collection sizes, but are not
  a separately metered exact heap-byte budget.

## Invariants

1. **Round-trip losslessness**: `parse_ir(write_ir(s)).scenario() == &s` for every
   representable scenario; `write_ir` is byte-stable canonical v2. Historical
   v1 bytes are accepted without mutation and their migration context is not
   discarded; every legacy receipt records and verifies exact `old_hash →
   new_hash` evidence.
2. **Dimensional soundness**: `validate()` rejects any BC/frame/ensemble/
   environment value whose SI exponents disagree with the contract table.
3. **Net-flux compatibility**: if any condition declares
   `incompressible`, either declared mass flows balance to 1e−9 relative
   or a pressure outlet exists — otherwise `flux-imbalance` with the
   imbalance quantified in the message and an actionable fix. Evaluation
   failure and non-finite aggregation are explicit violations; neither is
   silently reinterpreted as zero flow, and a pressure outlet cannot mask a
   malformed inlet declaration.
4. **Frame chains terminate**: cycles and dangling parents are violations;
   `world_pose` refuses cyclic chains at runtime too (hop budget).
5. **Bitwise ensemble replay (G5)**: identical complete canonical ensemble
   specifications, member identities, and implementation stream/synthesis
   semantics produce identical realization bits. Distinct member/seed streams
   are domain-separated; retained non-degenerate fixtures verify that they
   produce distinct draws. The current API does not yet expose the versioned
   recipe/seed-tree receipt tracked by `frankensim-sj31i.39`. Spectral grids
   with fewer than two samples are refused instead of publishing the
   seed-independent zero trace produced by an empty harmonic basis.
6. **Statistical spectrum match**: the ensemble-averaged periodogram of
   Kanai–Tajimi realizations converges to the target PSD (conformance
   holds band-mean relative error < 15% at 48 members with fixed seed).
7. **Nothing defaulted silently**: `Scenario::new` requires an
   `Environment`; `Environment::earth_lab()` exists but must be named at
   the call site.
8. **Finite authority values**: validation rejects non-finite signals, frame
   components, environment values, BC values, contact coefficients, spectral
   parameters, and Carreau bounds. Absolute ambient temperature and pressure
   are nonnegative; frame axes/quaternions are finite and unit length; Carreau
   viscosity/time bounds are positive, its shear-thinning index lies in
   `(0, 1]`, and every admitted independent draw satisfies
   `eta_zero >= eta_inf`.
9. **Evaluation APIs fail closed**: `TimeSignal::eval`,
   `BoundaryCondition::mass_flow_at`, `FrameTree::{local_pose,world_pose}`, and
   `StochasticEnsemble::{realize,realize_with_budget}`
   revalidate the public data they consume. Missing, wrong-physics, and
   geometry-unintegrated mass-flow values are typed refusals. Calling these
   APIs directly on a malformed public struct cannot bypass scenario-level validation. IR parsing
   independently enforces wire/resource/numeric-constructor safety, but is not
   semantic scenario admission; callers still invoke `Scenario::validate`.
10. **Identity integrity**: every scenario/frame/case/combination/ensemble/
    region identity is nonempty; names are unique within their role;
    combination terms do not repeat a case; and an unordered contact pair has
    exactly one declared model. Diagnostics retain first and repeated rows.
11. **Indexed structural validation**: frame identity indexes are built once;
    the tri-color parent traversal visits each storage row at most once, and
    case/frame/combination/ensemble/contact reference checks use deterministic
    ordered indexes rather than repeated prefix or whole-collection scans.
12. **Semantic preflight before validation**: top-level collection caps precede
    nested traversal; checked plans account for aggregate case BCs, combination
    terms, dynamic signal scalars, raw flux checkpoints, exact identity bytes,
    ordered-index comparisons, checkpoint sorting, and flux evaluation work.
    Exact requested limits admit and one-unit-short limits refuse for every
    budget field.

## Error model

`ScenarioError`: `Dimensions { context, expected, got }`, `Frame`,
`Evaluate`, `Parse { at, what }`. Parse/evaluation errors include deterministic
budget refusals and allocation-refusal context. `ValidationError` distinguishes
named limit refusal, work-plan overflow, total-work refusal, scratch-allocation
refusal, and cancellation with phase/completed/planned work. Admitted validation produces
`Vec<Violation>` (code + what + fix) rather than failing fast — agents get the
whole repair list at once.

## Determinism class

**D0**: signal evaluation, frame poses (fs-ga + fs_math::det trig), and
ensemble realizations (Philox + det trig, fixed draw/summation order) are
bit-identical across runs and platforms. IR text is canonical.

## Cancellation behavior

All parsing is bounded by explicit byte/depth/node/atom/list limits. Spectral
realization is O(samples × harmonics), admitted against explicit sample/work
budgets before allocation; collection reservations are fallible. Semantic
validation preflights all public collection families, dynamic signal payload,
identity bytes, raw flux-checkpoint allocation shape, and deterministic work
before executing. Net-flux validation streams base and optional-case slices
without materializing a vector for every effective set; its exact raw
checkpoint capacity is fallibly reserved before append/sort. Its
identity/reference phase uses deterministic O(N log N) indexes and its
frame-cycle traversal is linear after indexing. The explicit `Cx` lane polls
before preflight, at every top-level and nested record visited while constructing
the semantic plan, after planning, after fixed phases, at every frame-index row,
frame-cycle traversal/finalization step, and frame validation row, at
BC/case/combination term/ensemble/contact boundaries, before and after each
net-flux provider evaluation, at every tabulated signal scalar and Chebyshev
coefficient, and after private validation before publication. Tabulated signal
validity and ordering are accumulated in one pass while retaining diagnostic
order. A request observed at any checkpoint publishes no partial findings.
Fallible index/output reservation remains active work under
`frankensim-sj31i.24`. No loop is admitted from an unchecked float-to-size
conversion.

## Unsafe boundary

Zero `unsafe`.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs` (JSON verdicts, suite `fs-scenario/conformance`):

- **sc-001** rich representative vessel-pour fixture round-trips
  memory ↔ versioned v2 IR ↔ fs-ledger artifact losslessly; canonical text
  byte-stable. A retained legacy fixture verifies implicit/explicit v1
  five-vector decoding, `mol = 0`, pinned distinct old hashes, a shared pinned
  canonical hash, and receipt verification/tamper refusal.
- **sc-002** seeded violations caught with structured fixes: flux
  imbalance (repaired by adding an outlet), cyclic + dangling frames,
  wrong-dimension BC, undeclared inlet compatibility, unknown combo case,
  kinetic > static friction.
- **sc-003** KT/Dryden/Carreau members bitwise-identical across repeated
  realization; members differ from each other; seed matters; Carreau
  draws stay inside declared bands.
- **sc-004** 48-member KT ensemble periodogram vs target PSD: per-bin
  and band-mean tolerances hold (metrics logged as JSON).
- **sc-005** G0 frame laws: chain composition equals the manual motor
  product; the tilt ramp matches a directly-built motor at five times,
  clamps past its end; points on a rotation axis stay fixed.
- **sc-006** G3 unit coherence: deg/rad and mm/m spellings (via
  fs-qty parsing) converge to the same SI values and the same canonical
  IR; validation is spelling-blind.
- **sc-007** adversarial parser budgets hit byte/depth/node/atom/list limits at
  deterministic boundaries; malformed Chebyshev domains/coefficients and
  non-finite wire numbers return `Parse` without unwinding.
- **sc-008** direct signal/frame/ensemble result APIs and whole-scenario
  validation refuse non-finite/domain-invalid public structs. Realization
  sample/work budgets are checked at exact boundary and boundary+1, and valid
  Carreau draws remain finite, in-band, and physically ordered. A one-sample
  spectral grid is rejected before it can masquerade as a stochastic trace.
- **sc-009** exact non-ASCII identities remain admissible, while empty and
  duplicate role identities, repeated combination terms, self-contact, and
  duplicate/conflicting unordered contact pairs fail closed with declaration
  provenance.
- **sc-010** the retained semantic plan replays exactly; frames, base/case BCs,
  cases, combinations/terms, ensembles, contacts, signal scalars, raw flux
  checkpoints, identity bytes, and total work each admit at the exact boundary
  and refuse one unit short; pre-requested cancellation publishes no findings.
- A focused `scenario` unit regression forces checkpoint-capacity overflow and
  proves a typed `AllocationRefused` with the scratch vector left empty.
- A deterministic injected-checkpoint regression cancels at the environment
  phase and proves the private finding buffer is not returned; the public `Cx`
  path uses the same checkpoint route.
- A deterministic injected preflight regression cancels while counting nested
  case boundary conditions and proves no semantic plan is published.
- A focused frame regression injects cancellation inside the tri-color cycle
  walk and proves frame findings remain private at that boundary.
- A focused table-signal regression proves the checkpointed scalar traversal
  is one pass, preserves public diagnostic order, and observes injected
  cancellation before findings escape.

## No-claim boundaries

- **Physics vocabulary is v0**: IncompressibleFlow / Thermal /
  Elasticity kinds only. New physics extend `expectation` — adding a
  (physics, kind) pair is a table row plus tests, not a redesign.
- **Region names are strings here**: binding to fs-geom `Region` objects
  (existence, patch-measure integration for velocity-inlet flux) happens
  in the consuming solver layer; net-flux checking covers DECLARED
  mass-flow values, not velocity-profile surface integrals.
- **Recorded ground-motion suites (PEER-class) are not bundled**: the
  `Table` signal is the container for imported records; curation of suites is
  data, not code, and lives with fs-uq.
- **No load-combination EVALUATION**: combinations are typed references
  with factors; assembling factored response quantities is solver-side.
- **The ledger `scenarios` integration is a thin artifact row** (canonical
  IR + seed); a dedicated relational table is deferred to the ledger's
  next schema migration if queries demand it.
- **IR resource budgets are syntactic, not semantic or byte-exact heap
  admission**: parsing may intentionally return a finite but dimensionally
  invalid `Scenario` so migration/diagnostic tooling can inspect it; call
  `Scenario::validate` before solver admission. Input/node/list limits bound
  decoder cardinality. Semantic validation has its own explicit work plan and
  must not reuse syntax limits as an admission receipt. Exact decoded-heap
  accounting plus fallible semantic-index/finding reservation remain active
  work under `frankensim-sj31i.24`.
