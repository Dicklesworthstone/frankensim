# fs-scenario — CONTRACT

The boundary-condition and load-case ALGEBRA (plan patch Rev D): a
`Scenario` is a typed value answering "what is being done to the domain?"
— with dimensional analysis on every value, provenance (seed + canonical
IR), and admission-time validity checks that catch the class of mistakes
no solver can fix.

Ambition tags: typed BCs/frames/signals/combos [S]; seeded ensembles
(Dryden, Kanai–Tajimi, Carreau bands) [S]; canonical IR [S]; persistent entity
identity with rebinding receipts [S].

## Purpose and layer

Layer **L3** (FLUX support). Runtime deps: `std`, fs-blake3, fs-qty,
fs-rand, fs-cheb, fs-exec, fs-ga, fs-ivl, fs-motion, fs-math. The Design Ledger stores scenarios
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
  off-origin axis is `T(c)·R·T(−c)`. `rotating_motor_path(frame)` is the
  one-way L3-to-L2 admission boundary: it lowers a rotating target with only
  fixed ancestors into an owned `FrameTreeMotorPath` implementing
  `fs_motion::LowerToMotorTube`. It refuses a non-rotating target or any
  dynamic ancestor instead of sampling a general composed path into a
  constant-screw claim. `fs-motion` does not import scenario types.
- `payload` — a closed, versioned algebra for scalar, fixed-width vector,
  rectangular tensor, paired complex-phasor, canonical species-bundle,
  directed characteristic-state, field-trace-reference, and component-port
  payloads. Every payload declares a six-base dimension or sealed `fs-qty`
  semantic kind, canonical basis id, frame, orientation parity, and either
  continuous or named-reset reference semantics. Numeric carriers admit fixed
  values, strictly increasing time tables with explicit interpolation and
  outside-domain policy, or explicitly tagged distribution parameters/support.
  `Continuous` means the reference origin is not reset; it is not a claim that
  a `StepLeft` sample path is mathematically continuous. Complex phasors always
  carry one phasor-capable semantic kind and one Peak/RMS convention across all
  samples. Reference variants are typed indirections and therefore do not
  duplicate the referenced value's sample source.
  Constructors validate finiteness, dimensions, semantic scalar domains,
  stable shapes/axes, distribution order, and aggregate item bounds. The
  canonical payload V1 decoder enforces byte, aggregate-item (including the
  retained species-axis cache), and identifier limits before allocation and
  reconstructs only through those constructors. Its conservative default wire
  ceiling is derived from the closed item/id ceilings and remains large enough
  for every admitted V1 payload. Canonical scenario v2 IR embeds those exact
  bytes as lowercase hex under an independent payload-version tag. Scenario v1
  refuses the typed form, preserving its immutable five-to-six-base dimension
  crosswalk instead of guessing how a typed payload should migrate. The
  default parser's per-atom ceiling equals its already-bounded 16 MiB total
  input ceiling, so a canonical typed atom cannot trip an unrelated historical
  1 MiB limit while the complete artifact remains admitted. The payload
  encoder exposes an allocation-free exact byte counter and a fallible
  exact-reservation materializer; the historical infallible convenience
  delegates to that single codec grammar.
- `bc` — `BoundaryCondition { region, physics, kind, value, compatibility,
  frame }`; `expectation(physics, kind)` is the closed dimensional and carrier
  contract table. Existing rows retain velocity for flow Dirichlet, kg/s for
  mass-flow inlets, Pa for pressure/traction, K / W/m² / W/(m²K) for thermal,
  m for elastic Dirichlet, and no-value wall kinds. Typed-only rows add magnetic
  vector potential (vector Wb/m), normal magnetic flux density (scalar T),
  electric potential (scalar V), normal current density (scalar A/m²), species
  amount/mass flux bundles, and heterogeneous incoming/outgoing gas
  characteristic states. Every other pair is structurally Unsupported; generic
  Dirichlet/Neumann aliases do not smuggle these meanings. Payload and boundary
  frame ids must agree. Flux-carrying total-flow inlets MUST declare
  `Compat::Incompressible`; typed payloads never stand in for that legacy
  uniform/time-signal kg/s declaration, and all new typed rows forbid the
  incompressible compatibility tag. A spatial total-flow profile is refused
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
- `entity` — persistent ENTITY IDENTITY: `Assembly -> (Assembly | Part)`,
  `Part -> (Region | Surface)`, and `Interface` under the assembly or part that
  contains both of its sides. `EntityId { kind, digest }` is derived by BLAKE3
  over a length-prefixed, domain-separated preimage of the kind tag, the parent
  identity (which transitively encodes the whole parent path), the declared
  name, the optional `GeometryFingerprint`, and the optional `InterfacePair`.
  `EntityDeclaration::identity` is pure — no catalog is consulted — so a caller
  can compute an identity before declaring it. DISPLAY NAMES ARE NOT IN THE
  PREIMAGE: `EntityCatalog::rename` moves the display name and appends a
  `Rename` receipt while every binding stays on the same identity.
  `InterfacePairing::Ordered` keeps declaration order (`applied_side` is the
  side a one-sided treatment such as a TIM is applied to) and makes `(a, b)`
  and `(b, a)` distinct identities; `Unordered` canonicalizes the pair into one
  identity and `applied_side` returns `None` rather than defaulting to a side.
  `EntityCatalog` stores entities in declaration order with a sorted identity
  index, refuses a repeated declaration as `DuplicateDeclaration` and a
  different declaration landing on a present digest as `IdentityCollision`,
  and appends every identity event to a hash-chained `IdentityReceipt` log
  (`receipt_root`, `verify_receipts`, per-receipt `verifies`).
  `ImportRevision { label, event, scope, entities }` applies an import,
  re-mesh, or revision migration ATOMICALLY: `Correspondence::Auto` matches on
  geometry fingerprint first and declared name plus corresponding parent path
  second, refuses a tie as `AmbiguousImportMatch` instead of taking a first
  match, and `ImportScope::Complete { root }` retires active strict descendants
  of `root` the revision omits (the root itself is never auto-retired) while
  `ImportScope::Partial` leaves them untouched. `resolve` follows supersession
  links under an explicit hop budget and reports the WEAKEST `EvidenceTier`
  along the chain. `EntityRef`/`ReferenceSite`/`BindingTable` carry scenario
  references BY IDENTITY, and `validate_bindings` resolves all of them into the
  same `Violation { code, what, fix }` shape. `migrate_legacy_scenario` turns
  every distinct region string into a declared-name `Surface` carrying the
  legacy marker under one synthetic part; the marker is metadata, not identity,
  so re-running the migration re-derives the same identities and appends no
  receipts. `Datum`/`Tolerance`/`Placement` are typed declarations: datum
  hierarchies are acyclic by construction (a datum's identity depends on the
  identities it references), tolerance magnitudes must be finite positive
  lengths, form controls forbid a datum frame while orientation/location
  controls require one of at most `MAX_DATUM_FRAME_LEN` distinct datums, and a
  placement names a scenario `FrameId` rather than introducing a parallel frame
  system (its existence and uniqueness in the scenario's `FrameTree` are
  checked by `validate_bindings`, not at declaration time, because the catalog
  does not hold the scenario).
- `ir::write_ir`/`ir::parse_ir` — canonical byte-stable, explicitly versioned
  s-expression encoding. v2 writes six-base `[m, kg, s, K, A, mol]` vectors;
  explicit v1 and historical unversioned five-vector forms decode by appending
  `mol = 0` and return a `DecodedScenario` carrying the immutable
  `DimensionCrosswalkReceipt`. The receipt BLAKE3-binds the exact supplied v1
  bytes to the exact canonical v2 re-emission and exposes a verifier; explicit
  and implicit v1 spellings therefore retain distinct source identities.
  Accepted current-version bytes that differ from the canonical writer output
  return a separate `SourceCanonicalizationReceipt`. It BLAKE3-binds the exact
  supplied v2 bytes to the exact canonical v2 re-emission, names the versioned
  re-emission rule, and refuses source or target tampering. Exact canonical v2
  writer output carries no redundant receipt. Thus alternate whitespace,
  integer/version spellings, finite float spellings, and signed zero may be
  inspected without laundering distinct authority bytes into one silent
  identity. Unsupported string escapes remain refusals rather than aliases.
  Floats use shortest-round-trip form and u64 seeds remain exact integers.
  Physically irrelevant signed zero is canonicalized to `0`, matching the
  scenario types' semantic `PartialEq` rather than creating two content
  identities for equal values.
  Strings use exactly the writer's quote and backslash escapes; every other
  backslash sequence is rejected so distinct authority bytes cannot alias one
  decoded identity through an undocumented escape rule.
  `parse_ir` applies `DEFAULT_IR_PARSE_BUDGET` and
  `DEFAULT_IR_DECODE_BUDGET`; `parse_ir_with_budget` exposes
  explicit byte, recursive-depth, total-node, atom/string-byte, and per-list
  child limits while retaining the default end-to-end decode authority.
  `plan_ir_decode` derives a checked conservative syntax/semantic/output/work
  plan directly from source cardinality, and
  `parse_ir_with_resource_budget` admits that plan before syntax-tree
  allocation. Every syntax limit is checked before the corresponding recursive or syntactic
  tree growth; caller-selected depth may tighten but cannot exceed the
  recursive implementation's hard safety ceiling. Non-finite wire numbers and
  invalid Chebyshev constructor inputs are structured parse refusals, never
  panics. Syntax-tree growth, atom buffers, every decoded collection vector,
  and every decoded string-field copy reserve fallibly before population. The
  decoded plan charges retained syntax, semantic slots/text, typed-payload
  expansion/scratch, canonical receipt output, and deterministic work under
  separate caps. Resource refusals report operation, phase, requested units,
  cap, and completed/planned work; no partial scenario or receipt is returned.
  `write_ir_plan` computes exact canonical text, payload, peak-logical-heap,
  and work cardinalities without materializing bytes;
  `write_ir_with_budget` reserves the exact output before emission and encodes
  payload scratch sequentially under that admitted peak. Every syntax node retains its exact half-open
  source span. Semantic decoding adds a deterministic `$`-rooted structural
  path only while unwinding a refusal, so nested field/index diagnostics do not
  allocate on the green path. Parse and reserved-machine-role refusals expose
  both the path and span; typed-payload hex nibble offsets are translated back
  to absolute scenario-source bytes rather than reported relative to the
  embedded string.

## Invariants

1. **Round-trip losslessness and source identity**:
   `parse_ir(write_ir(s)).scenario() == &s` for every representable scenario;
   `write_ir` is byte-stable canonical v2. Source evidence is mutually
   exclusive: exact canonical v2 has no receipt, accepted noncanonical v2 has
   one source-canonicalization receipt and no migration receipt, and legacy v1
   has one dimension-crosswalk receipt and no source-canonicalization receipt.
   `DecodedScenario` equality includes this provenance; callers that mean only
   semantic equality compare the decoded `Scenario` values.
   Historical v1 bytes are accepted without mutation and their migration
   context is not discarded; every legacy receipt records and verifies exact
   `old_hash → new_hash` evidence.
2. **Dimensional soundness**: `validate()` rejects any BC/frame/ensemble/
   environment value whose SI exponents disagree with the contract table.
3. **Net-flux compatibility**: if any condition declares `incompressible`,
   either declared mass flows balance to 1e−9 relative at every instant or a
   pressure outlet exists. This layer certifies the all-time branch only for
   uniform and `TimeSignal::Constant` total flows. Ramp, Linear/Hold Table, and
   Chebfun checkpoints remain deterministic falsifiers: one sampled point may
   produce a concrete `flux-imbalance`, but a green finite screen refuses as
   `flux-certification-unavailable` instead of impersonating an all-time
   floating-point proof. Evaluation failure and non-finite aggregation are
   explicit violations; neither is silently reinterpreted as zero flow, and a
   pressure outlet cannot mask a malformed inlet declaration.
4. **Frame chains terminate**: cycles and dangling parents are violations;
   a parent reference to a duplicated frame id is explicitly ambiguous and is
   never resolved through an arbitrary storage row;
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
    exactly one declared model. A pair group with multiple models classifies
    every repeated row as a conflict regardless of declaration permutation; an
    all-equal group classifies repeats as duplicates. Diagnostics retain first
    and repeated rows. Combination-term diagnostics copy at most 128 source
    bytes from each combination or case identity, append the exact original
    byte length when truncated, and retain combination-row/term coordinates so
    equal UTF-8-safe previews do not erase exact declaration provenance.
    Net-flux set labels allocate no eager case-name copy; findings render the
    same bounded identity preview and retain the exact case row.
11. **Indexed structural validation**: frame identity indexes are built once;
    the tri-color parent traversal visits each storage row at most once, and
    follows only uniquely resolved parent ids so duplicate declarations cannot
    make cycle findings depend on declaration order. The frame structural
    checker is crate-private and runs only after whole-scenario budget/work
    preflight; there is no public unadmitted `FrameTree::check` escape hatch;
    case/frame/combination/ensemble/contact reference checks use deterministic
    ordered indexes rather than repeated prefix or whole-collection scans.
12. **Semantic preflight before validation**: top-level collection caps precede
    nested traversal; checked plans account for aggregate case BCs, combination
    terms, dynamic signal and typed-payload scalar slots, typed-payload identity
    bytes, raw flux checkpoints, all other exact identity bytes, worst-case
    finding slots, ordered-index comparisons, checkpoint sorting, and flux
    evaluation work. String comparison work is charged per identity role using
    its maximum key width: twice the checkpointed heap-sort envelope plus both
    operands of every subsequent ordered lookup. Combination
    references additionally charge the sum of the case/reference maxima for
    their cross-lookup into the case index; contact keys charge both
    canonicalization passes, self-pair comparison, and grouped adjacency.
    Empty identities retain a one-unit comparison width. Each identity or
    reference component has a hard 4,096-byte cap that callers may tighten but
    cannot raise. Numeric frame/parent and BC-frame index lookups are charged
    separately. Exact requested limits admit and one-unit-short limits refuse
    for every budget field.
    Valid boundary-condition, frame, and ensemble validation carry diagnostic
    context as borrowed formatters, so the green path does not allocate or copy
    identity names merely in case a finding is needed.

13. **Identity is a pure function of the declaration preimage**: equal
    `(kind, parent, declared name, fingerprint, interface pair)` tuples derive
    equal `EntityId`s and nothing else does within the checked catalog. Every
    variable-length component is length-prefixed, so no two distinct
    declarations share a preimage by concatenation. Display names and the
    legacy marker are excluded, which is exactly why a rename cannot orphan a
    binding and why mechanical legacy migration is idempotent.
14. **One successor per identity**: a superseded identity records exactly one
    successor (`AmbiguousSupersession` otherwise), a superseded or retired
    identity can never be redeclared or revived
    (`InactiveIdentityRedeclared`), and resolution therefore walks a chain, not
    a graph. Resolution reports the WEAKEST `EvidenceTier` on that chain: a
    content-matched hop after an asserted hop still resolves as `Asserted`.
15. **Revisions are atomic**: `apply_import` either applies every step with its
    receipts or leaves the catalog exactly as it was; a refused revision
    mutates nothing, including the receipt log.
16. **The receipt log is an append-only hash chain**: receipt `i` has sequence
    `i` and binds the previous root, so an interior edit, a reorder, or a
    dropped receipt fails `verify_receipts`.
17. **Referential integrity**: `validate_bindings` resolves EVERY binding and
    enumerates every base-BC, case-BC, and contact-side site. Unresolved,
    retired, over-deep, wrong-kind, duplicated, orphaned, unbound, and
    string-ambiguous references are structured findings with fix hints; there
    is no first-match fallback anywhere on the path from a reference to an
    entity.

## Error model

`ScenarioError`: `Dimensions { context, expected, got }`, `Frame`,
`Evaluate`, `Parse { span, path, what }`, and
`ReservedBoundaryRole { role, span, path }` for an attempt to encode a
Machine-IR joint, terminal, controller, or reset as a BC kind. `span` is an
exact half-open byte range in the supplied source and `path` is a deterministic
structural address such as `$.cases[0].bcs[1].kind`. Parse/evaluation errors
include deterministic budget refusals and allocation-refusal context.
`ValidationError` distinguishes
named limit refusal, work-plan overflow, total-work refusal, scratch-allocation
refusal, and cancellation with phase/completed/planned work. Admitted validation produces
`Vec<Violation>` (code + what + fix) rather than failing fast — agents get the
whole repair list at once.

`EntityError` is the entity layer's typed refusal set (name shape, unknown or
inactive parent, containment, interface pair shape/side/scope, duplicate
declaration, digest collision, inactive-identity redeclaration, unknown entity,
unknown or duplicate datum, datum-owner kind, placement-occurrence kind,
ambiguous import match, ambiguous supersession, inactive predecessor,
non-revision event, ambiguous or unknown declared name, capacity, allocation,
hop budget). Every variant exposes a stable `code()`, a ranked `fix()`, and
`into_violation()`, so an entity refusal enters the same diagnostic shape as
the rest of the crate rather than a second one. `ResolutionFault`
(`Dangling`, `Retired`, `DepthExceeded`, `KindMismatch`) is the resolution-time
subset and carries the same `code`/`fix` pair.

## Determinism class

**D0** for the entity layer: identity derivation, datum identity, and receipt
digests are pure BLAKE3 over byte-exact preimages with no floating point, no
hashing-order dependence (entities, receipts, datums, tolerances, placements,
and bindings are all declaration-ordered `Vec`s with sorted lookup indexes),
and no wall-clock or address-derived input. Two independent builds of the same
model produce identical identities and identical receipt roots.

**D0**: signal evaluation, frame poses (fs-ga + fs_math::det trig), and
ensemble realizations (Philox + det trig, fixed draw/summation order) are
bit-identical across runs and platforms. Reconstructing the admitted
constant-screw lowering with the same frame declaration and tube parameters is
bit-identical. IR text is canonical. Payload V1 bytes are canonical for every
admitted value: signed zero is normalized on write and rejected when supplied
noncanonically on the wire; all other finite IEEE-754 bits are retained.

## Cancellation behavior

All parsing is bounded by explicit byte/depth/node/atom/list limits. Spectral
realization is O(samples × harmonics), admitted against explicit sample/work
budgets before allocation; collection reservations are fallible. Semantic
validation preflights all public collection families, dynamic signal and typed
payload scalar slots, typed-payload and scenario identity bytes, raw flux-
checkpoint allocation shape, and deterministic work before executing. Net-flux
validation streams base and optional-case slices
without materializing a vector for every effective set; its exact raw
checkpoint capacity is fallibly reserved before append/sort. Its
checkpoint list sorts in place without hidden scratch allocation, using the
same deterministic checkpointed heap-sort skeleton as semantic indexes and
`f64::total_cmp` for a total ordering. The preflight charges its conservative
comparison/swap envelope before validation starts. A checkpointed in-place
deduplication pass canonicalizes both signed-zero encodings to `+0.0`. Its
identity/reference phase uses deterministic O(N log N) indexes. Frame ID/name,
frame-membership, case, combination, term, ensemble, and
unordered-contact indexes plus linear cycle traversal scratch are exactly and
fallibly reserved before population. Contact conflict flags are also fallibly
reserved and populated by one grouped pass over the pair index. Every index
sorts in place with row index as the total-order tiebreaker using a deterministic
checkpointed heap sort. The preflight work total counts each index population
item and a conservative heap comparison/swap envelope, including
per-combination term entries and every net-flux checkpoint set. The explicit
per-component cap bounds each opaque string comparison to 4 KiB (8 KiB across
both components of an unordered contact key). One heap-sift checkpoint covers
at most two such comparisons; an opaque ordered lookup covers at most
`usize::BITS + 1` comparisons including final equality. The widest opaque
comparison envelope between surrounding record checkpoints is one two-part
contact lookup plus two direct component comparisons, conservatively 528 KiB
on a 64-bit target; bounded diagnostic rendering is additional work outside
that comparison envelope. Two single-component combination-term lookups fit
below the same comparison envelope. The
explicit `Cx` lane polls
before preflight, at every top-level and nested record visited while constructing
the semantic plan, at every typed-payload source sample and retained identity
component traversed during that plan, after planning, after fixed phases, at every frame-index row,
every frame-cycle scratch-initialization, traversal, and finalization step, and
frame validation row, at
BC/case/combination term/ensemble/contact boundaries, at every unordered
contact-conflict scratch-initialization step and classification group (including
all-unique inputs), before and
after each net-flux provider evaluation, at every tabulated signal scalar and
Chebyshev coefficient during structural scanning, before every nonconstant
Chebyshev recurrence step during net-flux evaluation, throughout index
population/sort steps, and after private
validation before publication. Net-flux checkpoint sorting polls throughout
the in-place heap sort rather than only after an opaque library sort; set
classification, checkpoint counting/materialization, and deduplication likewise
poll at every provider or checkpoint. Tabulated signal validity and ordering are
accumulated in one pass while retaining diagnostic order. A request observed at
any checkpoint publishes no partial findings. Preflight proves a conservative
finding bound of 12 fixed slots plus 13/frame, 8/BC, 3/case, 2/combination,
4/term, 16/ensemble, and 4/contact; `max_findings` admits that heap authority and
the private result vector fallibly reserves it before the first finding. No loop
is admitted from an unchecked float-to-size conversion.

After that single checkpointed structural scan, net-flux compatibility uses a
crate-private prevalidated signal evaluator. Public `TimeSignal::eval` and
`BoundaryCondition::mass_flow_at` remain independently fail-closed, while the
whole-scenario path retains O(log N) table lookup without rescanning all N
samples at every one of N validation checkpoints. The prevalidated evaluator
still checks finite time, table shape, lookup bounds, and finite results. The
work plan charges each table provider by its binary-search height, each
Chebyshev provider by coefficient count, and each materialization/deduplication
and set-scan pass explicitly.

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
- **sc-001f** accepted current-version whitespace, version, integer, finite
  float, and signed-zero aliases re-emit to one canonical v2 artifact while
  retaining distinct source hashes in versioned tamper-checked receipts;
  exact canonical writer output remains receipt-free.
- **sc-001g** nested ensemble scalars, frame-motion arity, cross-field profile
  domains, reserved/unknown BC roles, and embedded typed-payload hex mutations
  pin exact structural paths and absolute half-open source spans.
- **sc-002** seeded violations caught with structured fixes: flux
  imbalance (repaired by adding an outlet), cyclic + dangling frames,
  wrong-dimension BC, undeclared inlet compatibility, unknown combo case,
  kinetic > static friction.
- **sc-002a, sc-002e–i** deterministic ramp/table/Chebfun screens catch concrete
  time-varying flux counterexamples, while adversarial green Hold,
  relative-tolerance, and smooth-grid fixtures refuse without an all-time
  certificate. A pressure outlet is the explicit alternative; uniform and
  constant totals retain the exactly checkable path.
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
- **sc-007a** exact writer and conservative decoder heap/output/work plans
  admit the retained collection-rich scenario with UTF-8, quote/backslash
  escapes, and a typed payload; every one-short cap returns a typed preflight
  refusal before publication. The same exact/one-short boundary covers legacy
  v1 migration and its receipt-producing canonical re-emission.
- **sc-008** direct signal/frame/ensemble result APIs and whole-scenario
  validation refuse non-finite/domain-invalid public structs. Realization
  sample/work budgets are checked at exact boundary and boundary+1, and valid
  Carreau draws remain finite, in-band, and physically ordered. A one-sample
  spectral grid is rejected before it can masquerade as a stochastic trace.
- **sc-009** exact non-ASCII identities remain admissible, while empty and
  duplicate role identities, repeated combination terms, self-contact, and
  duplicate/conflicting unordered contact pairs fail closed with declaration
  provenance; mixed-model pair classification is invariant under model/order
  permutations.
- **sc-010** the retained semantic plan replays exactly; frames, base/case BCs,
  cases, combinations/terms, ensembles, contacts, signal scalars, raw flux
  checkpoints, aggregate and per-component identity bytes, finding slots, and
  total work each admit at the exact boundary and refuse one unit short;
  pre-requested cancellation publishes no findings. A focused common-prefix
  regression admits exactly 4,096 bytes and refuses both a 4,097-byte component
  and any caller attempt to widen the hard comparison tile.
- **sc-011** all 1,296 four-frame parent graphs match an independent exhaustive
  chain oracle, including world roots, dangling parents, self-cycles, and
  multi-node cycles.
- **sc-012** 100k-frame deep, wide, and cyclic adversaries complete under the
  admitted `Cx`; deep/wide graphs remain green, every cyclic row is diagnosed,
  wall times are emitted, and deterministic planned work grows by less than 3x
  for each 25k -> 50k -> 100k doubling.
- A focused `scenario` unit regression forces checkpoint-capacity overflow and
  proves a typed `AllocationRefused` with the scratch vector left empty.
- A deterministic injected-checkpoint regression cancels at the environment
  phase and proves the private finding buffer is not returned; the public `Cx`
  path uses the same checkpoint route.
- A deterministic injected preflight regression cancels while counting nested
  case boundary conditions and proves no semantic plan is published.
- A representative phase-matrix regression discovers every reached frame,
  signal, BC/case/combination/ensemble/contact/index/net-flux checkpoint and
  injects cancellation at each distinct phase, proving private findings never
  escape.
- A focused frame regression injects cancellation inside the tri-color cycle
  walk and proves frame findings remain private at that boundary.
- A focused frame permutation regression proves duplicated parent ids produce
  the same explicit ambiguity finding without declaration-order-dependent cycle
  findings.
- A focused frame scratch regression forces a capacity overflow and proves a
  typed allocation refusal without partial scratch state.
- A focused table-signal regression proves the checkpointed scalar traversal
  is one pass, preserves public diagnostic order, and observes injected
  cancellation before findings escape.
- A focused evaluator regression proves the prevalidated constant, ramp, table,
  and Chebyshev paths match public evaluation for valid signals while retaining
  malformed-table and non-finite-time guards.
- A focused Chebyshev flux regression injects cancellation inside the Clenshaw
  recurrence and proves the private finding buffer is not published.
- A focused checkpoint-deduplication regression proves per-element polling and
  canonical `+0.0` retention when both signed-zero encodings occur.
`tests/entity_identity.rs` (JSON verdicts, suite `fs-scenario/entity`):

- **ec-001** the representative cooling stack (assembly, three parts, two
  regions, seven surfaces, three interfaces, three datums, three tolerances,
  two placements) is expressed by identity and resolves every one of its eight
  reference sites; ordered interfaces name an applied side and the unordered
  one refuses to.
- **ec-002** six renames move display names only: every binding resolves to the
  same identity with the same hop count and evidence tier, and the declared
  names do not follow the display names.
- **ec-003** a full CAD-rename re-import supersedes 16 identities — 13 bodies on
  content-matched receipts and 3 interfaces on declared-path receipts — and
  every binding still resolves, one hop away, onto an entity with the same
  geometry fingerprint and kind but a different declared name.
- **ec-004** re-applying the identical revision mints no identity, appends 16
  `Unchanged` receipts, and leaves the rendered binding table byte-equal.
- **ec-005** dangling, retired (via a `Complete` revision), caller-kind,
  site-kind, orphan-site, duplicate-binding, unbound-site, and
  ambiguous-declared-name references are all reported with nonempty fix hints.
- **ec-006** a string-only scenario migrates mechanically: seven distinct
  region strings become seven legacy-marked declared-name surfaces with no
  invented fingerprints, re-running appends zero receipts, and a drifted region
  string is caught as `entity-legacy-name-drift`.
- **ec-007** datum hierarchy, tolerance dimensional typing, form/location datum
  rules, datum-frame arity and repetition, empty tolerance sources, placement
  round-trip through `FrameTree::world_pose`, and an undeclared placement frame.
- **ec-008** the lifecycle: author from strings, import-bind onto real
  geometry (path-matched), rename, re-mesh re-import, then validate — all eight
  bindings still resolve, the re-meshed face through two hops, and a reserved
  `Sensor` site binds through the same shape.
- **ec-009** identity derivation, receipt roots, and site enumeration are
  identical across two independent builds of the same model.
- In-crate `entity::identity_tests` cover preimage injectivity and
  length-prefixing, kind/parent/fingerprint sensitivity, interface ordering,
  containment rules, duplicate-versus-collision refusals (the collision branch
  is exercised by forcing a colliding row on both the declare and the import
  path), datum-owner and placement-occurrence kind checks, rename receipts, chain
  verification against field tampering/reordering/truncation, weakest-link
  resolution, hop budgets, ambiguous automatic matching, revision atomicity,
  and the declared-name ambiguity that motivates the whole module.
- `tests/typed_bc.rs` locks all eight typed expectation rows, carrier/kind/
  six-base-dimension/frame refusals, heterogeneous characteristic admission,
  typed-total-flow refusal, exact compatibility semantics, typed payload
  scalar/identity/work budget charging in both base and case paths, checkpointed
  preflight cancellation, and canonical v2 payload-IR round trips (including a
  typed atom above 1 MiB) with legacy/version/case/trailing-byte refusals.

## No-claim boundaries

- **Payload sources are declarations, not evaluators or stochastic proofs**:
  this crate validates table structure and distribution parameter/support
  shape, units, finiteness, and canonical ordering. It does not claim a table
  interpolator, sampler, correlation model, probability law certificate, or
  field/port existence proof. Reference identifiers are canonical syntax; the
  owning ledger/component layer must resolve and authorize them. Complex
  phasor distributions are empirical-only in V1 because independent
  Normal/Uniform parameters would not honestly encode paired covariance.

- **General frame-path lowering is not claimed**: the current one-way adapter
  covers one rotating target below fixed ancestors. Scheduled tilts, multiple
  rotating ancestors, and other composed dynamic chains require a certified
  general-path constructor or certified tube composition in `fs-motion`; they
  are structured refusals here. The current `LowerToMotorTube` constructor has
  no `Cx`, so cancellation-correct tube construction is also not claimed by
  this adapter.

- **Typed multiphysics rows are declarations, not solver claims**: Magnetics,
  Electrics, and GasExchange establish exact structural payload kinds,
  six-base dimensions (or a heterogeneous characteristic contract), frames,
  and wire semantics only. `fs-qty` does not yet expose sealed electromagnetic
  semantic kinds, so these rows make no stronger semantic-kind claim. Field
  equations, orientation transforms, conservation checks, and solver support
  remain in their owning FLUX crates. New pairs still require an explicit
  `expectation` table row plus tests; there is no open-ended fallback.
- **A content-derived identity proves byte equality of its preimage, and
  nothing about the physical world**: equal `EntityId`s prove the two
  declarations supplied the same kind, parent identity, declared name,
  fingerprint, and pairing. They are NOT a claim that two revisions describe
  the same physical part, that the entity exists in any geometry kernel, or
  that anyone measured anything. A `GeometryFingerprint` is BLAKE3 over exactly
  the bytes the importer chose to hash: this crate never parses, canonicalizes,
  meshes, or inspects geometry. Equal fingerprints therefore prove equal
  supplied bytes; UNEQUAL fingerprints prove nothing at all, because a
  re-export of an unchanged part routinely produces different bytes. That
  asymmetry is why `MatchBasis::proves_geometry_bytes_matched` is true only for
  `GeometryFingerprint` and why a fingerprint mismatch falls back to the weaker
  `DeclaredPath` tier rather than to "different part".
- **"Collision-checked" is a fail-closed check, not a collision claim**: on an
  identity hit the catalog compares the full stored declaration preimage and
  refuses with `IdentityCollision` when it differs, so two distinct
  declarations can never alias onto one entity. This is not a claim that
  BLAKE3-256 collisions are impossible, and the crate does not search for them.
- **A receipt records a claimed correspondence, not a witnessed operation**:
  the chain proves the retained log is internally consistent and detects any
  interior edit, reorder, or dropped entry, because the root is stored beside
  it. Truncating the log AND the retained root together is only detectable
  against a root pinned outside the catalog (the Design Ledger's job).
  `MatchBasis::Asserted` records that a caller asserted a correspondence; this
  crate verified nothing about it, and `EvidenceTier::Asserted` is the value
  every downstream reader gets for such a chain — including chains whose other
  hops are content-matched.
- **Entity bindings are not in canonical IR**: `Scenario` still carries its
  string `region`/`region_a`/`region_b` fields, and `write_ir`/`parse_ir`
  remain byte-identical to their pre-entity behaviour. The `EntityCatalog` and
  `BindingTable` are sibling values that a scenario is validated AGAINST; they
  do not round-trip through scenario IR yet, so a scenario reconstructed from
  canonical bytes alone has no bindings and `validate_bindings` will report
  every site as `entity-unbound-site`. Persisting the catalog and bindings is
  the versioned `.fsim` project schema's job (`frankensim-extreal-program-f85xj.6.1`).
- **`Scenario::validate` does not run entity validation**: `validate_bindings`
  is a separate explicit entry point. It does not participate in
  `ValidationPlan`/`ValidationBudget` work accounting, polls no fs-exec `Cx`,
  and is therefore not cancellation-correct. Entity admission has its own
  explicit `EntityBudget` caps (entities, receipts, name bytes, supersession
  hops, hierarchy depth, datums, tolerances, placements, bindings) and reserves
  fallibly before every catalog mutation, but diagnostic vectors are not
  byte-metered — the same gap the finding-capacity boundary below describes.
- **Reserved reference sites are not existence-checked**: `ReferenceSite::Load`,
  `MaterialBinding`, `Sensor`, and `Requirement` exist so the later scenario
  objects reuse ONE diagnostic shape. This crate cannot enumerate those
  collections, so their rows are not checked for existence and they are never
  reported as unbound; only their entity references are resolved.
- **Datums, tolerances, and placements are declarations, not evaluations**:
  this crate checks structure, dimensions, datum-frame arity, and source
  presence. It does not construct a datum reference frame geometrically,
  perform tolerance stack-up, decide whether a part conforms, or evaluate a
  drawing. `PlacementBasis::Nominal` names the frame that is DECLARED to carry
  an occurrence's placement; the crate does not verify that the frame's
  transform is where the part actually is, and there is deliberately no
  as-built variant until the as-built layer lands.
- **Legacy migration infers no structure**: a bare string carries no part
  hierarchy and no geometry, so migration puts every migrated surface under one
  synthetic legacy part with no fingerprint. It is a mechanical renaming of the
  reference layer, not a reconstruction of the assembly.
- **Region names are strings here**: binding to fs-geom `Region` objects
  (existence, patch-measure integration for velocity-inlet flux) happens
  in the consuming solver layer; net-flux checking covers DECLARED
  mass-flow values, not velocity-profile surface integrals. The `entity` module
  gives those references a stable IDENTITY; it does not resolve them to
  geometry, and an `EntityId` that resolves in the catalog says nothing about
  whether a chart, mesh, or patch exists for it.
- **Recorded ground-motion suites (PEER-class) are not bundled**: the
  `Table` signal is the container for imported records; curation of suites is
  data, not code, and lives with fs-uq.
- **No load-combination EVALUATION**: combinations are typed references
  with factors; assembling factored response quantities is solver-side.
- **The ledger `scenarios` integration is a thin artifact row** (canonical
  IR + seed); a dedicated relational table is deferred to the ledger's
  next schema migration if queries demand it.
- **IR heap accounting is portable logical admission, not allocator
  introspection**: parsing may intentionally return a finite but dimensionally
  invalid `Scenario` so migration/diagnostic tooling can inspect it; call
  `Scenario::validate` before solver admission. Decode admission conservatively
  counts requested vector/string capacities, retained typed values, canonical
  output, and sequential payload scratch; it does not claim allocator metadata,
  page commitment, fragmentation, or implementation-defined capacity rounding.
  The convenience `parse_ir` retains a 16 MiB total-input ceiling and the
  conservative default logical plan; an already-materialized larger canonical
  artifact needs explicit syntax and resource authorities, while the L6 Machine
  projection separately refuses scenario artifacts above its own 16 MiB
  transport bound. The
  semantic `check_round_trip` helper derives only byte/atom authority from the
  exact writer output and retains default structural ceilings. Semantic
  validation has its own explicit work plan and must not reuse syntax limits as
  a solver-admission receipt. `write_ir` preserves its historical infallible
  convenience signature by deriving exact caps internally; callers that must
  handle allocator refusal use `write_ir_with_budget`. Source and migration
  receipts remain identity evidence, not claims about physical resident memory.
- **Finding capacity is not exact diagnostic-heap admission**: bounded identity
  previews prevent one long combination name from being copied in full into
  every term finding, but other diagnostic fields and each final `String`
  allocation are not yet byte-metered. A global diagnostic-byte budget remains
  active work under `frankensim-sj31i.24`.
