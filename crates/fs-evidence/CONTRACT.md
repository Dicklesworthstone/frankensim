# CONTRACT: fs-evidence

## Purpose and layer
`Evidence<T>` / `Certified<T>` (patch Rev B): numerical, statistical, and
MODEL-FORM certificates carried inside values, with a conservative
composition algebra, model cards + the registration lint, two-fidelity
discrepancy models with out-of-distribution refusal, model bracketing, and
decision-aware escalation. The reason this crate exists: without model
evidence the system can produce beautifully certified WRONG answers (mesh
error 0.7%, closure discrepancy 10%) — being able to SAY that is the
product. Layer: UTIL (usable by every layer; the bead label said L6, but
the bead scope explicitly demands a low-layer home — this is it). Depends
on fs-blake3 for typed canonical identities and fs-obs for deterministic
telemetry/legacy correlation.

## Public types and semantics
- `Evidence<T> { value, qoi, numerical, statistical, model, sensitivity,
  provenance, adjoint_ref }` — the traveling noun. `Certified<T>` is the
  opaque newtype whose constructor discipline is `Evidence::certified()`: rigorous
  numerics (Exact/Enclosure), valid statistical parameters, and in-domain
  non-negative model discrepancy; pure math
  certifies with `ModelEvidence::none()` (the explicit "no model involved"
  statement); refusals are structured `CertifyError`s.
- `NumericalCertificate { kind: Exact|Enclosure|Estimate|NoClaim, lo, hi }`
  — the plan Appendix B `Certified` fields (value + interval bound +
  provenance + adjoint hook) kept intact as the numerical slice. Severity
  is ordered; float composition never claims Exact; NoClaim absorbs.
- `StatisticalCertificate { None | EValue{e, alpha} | HalfWidth{...} }` —
  finite non-negative e-values and widths with levels/confidences strictly in
  `(0,1)`; statistical v1 composition is conservative-weakest (see no-claims).
- `ModelEvidence { cards, assumptions, validity, discrepancy_rel,
  in_domain }`; `ValidityDomain` — named-parameter boxes with
  intersection/containment; `SensitivitySummary` — d(qoi)/d(param)
  headlines, merged by magnitude.
- `Evidence::combine(op, a, b, value)` — Add/Sub/Mul/Min/Max on the QoI
  with certificates composed conservatively and provenance chained
  (`ProvenanceHash::chain`, order-sensitive, FNV until the ledger hash).
- `Evidence::assess(threshold_rel) -> DecisionStatus` and
  `UncertaintyBreakdown` — per-source relative bands, first-order-sum
  total, dominant source with the declaration-order tie law (ModelForm
  first: ties escalate the model, the band cheap refinement cannot fix).
  `escalation_advice` maps dominance to RefineNumerics /
  GatherMoreSamples / EscalateModelFidelity (the HELM governor hook).
- `ModelCard` (name, version, ambition tag, assumptions, validity, known
  failures, calibration provenance, discrepancy band) and `ModelRegistry`
  — `register_solver` REFUSES without a registered card (the lint).
- `DiscrepancyModel::fit(&[FidelityPair])` — observed parameter box +
  mean/max relative discrepancy; `query`/`evidence_at` refuse
  out-of-distribution points and any query whose key set differs from the exact
  training schema, with the first missing, unexpected, non-finite, or
  out-of-range parameter named (`OutOfDomain`).
- `ModelBracket` — N plausible models; evidence = midrange value, an
  enclosure spanning every member, spread as the model band, and a
  bracket-spread sensitivity entry (the vessel flagship's contact-line
  mitigation).
- `to_ledger_row_json` on evidence and cards — the `evidence` /
  `model_cards` table rows (canonical order, no clocks, no addresses).
- `identity` module (sj31i.52.2 tranche 1) —
  `ColorEvidenceSourceIdV1` and `ColorEvidenceNodeIdV1` are nominal
  `SourceId` / `EvidenceNodeId` roles over static v1 framing schemas. Direct use
  of those low-level aliases with `CanonicalEncoder` proves only schema-shaped
  framing; it is not semantic admission. The opaque `ColorEvidenceSourceV1`
  helper result additionally enforces an explicit source-schema domain,
  nonzero source-schema version, and exact retained source bytes. The opaque
  `ColorEvidenceNodeV1` keeps a color attached to its receipt and is the type
  that carries the helper invariants: its source uses a complete child
  role/schema descriptor, while recursive parent rows bind the evidence-node
  role, complete node `SchemaId`, and full typed root. Add/Mul/Hull composition
  sorts parent IDs, preserves duplicates, recomputes the output with the current
  `IntervalOp` algebra, and then binds kind, operation, parent law,
  `COLOR_ALGEBRA_VERSION`, source/parents, and exact canonical output. Color
  encoding is hard-capped at 1 MiB, uses a fallible exact buffer reservation,
  and polls cancellation at the caller's byte stride while copying. Binary
  parent rows remain in fixed-size array storage rather than a fallible heap
  collection. Construction returns only unanchored receipts: it neither
  authenticates origin nor changes `Color`, `Certified<T>`, or `AdmittedColor`
  trust state.
- `identity` module (sj31i.52.2 tranche 2) —
  `ValidityDomainIdV1` is a low-level schema-shaped `SemanticId`; only the
  opaque `IdentifiedValidityDomainV1` proves helper validation of a normalized
  declared parameter box. The helper consumes and retains the exact domain
  while binding sorted, length-framed axis UTF-8 bytes without normalization
  plus both finite, ordered IEEE-754 endpoint bit patterns. Unconstrained is the
  canonical empty row sequence. Axis count and the exact ordered-field payload
  are admitted before allocation-free row streaming; the payload is hard-capped
  at 1 MiB and the complete frame separately obeys caller `CanonicalLimits`.
  The current scatter/gather producer preflights exactly four non-semantic
  stream chunks per axis against the caller's shared collection/chunk budget.
  Raw frame aliases remain low-level framing, not a semantic-admission bypass.
- `identity` module (sj31i.52.2 tranche 3) —
  `CertifiedF64EvidenceIdV1` is a strong `SemanticId` for the helper-defined
  semantic projection of one opaque `Certified<f64>`. The helper consumes and
  retains that certified record and binds its carried scalar and QoI,
  numerical kind and finite bounds, statistical variant payload, canonical
  model-card-name and assumption sets, the typed validity child, exact raw
  discrepancy bits, in-domain state, ordered sensitivity-name/derivative-bit
  rows, and whether an adjoint correlation is present. Card and assumption
  vectors must already be strictly byte-sorted and duplicate-free; malformed
  public vectors refuse rather than being silently reordered. The raw legacy
  FNV provenance and adjoint-correlation values are deliberately excluded so a
  weak collision cannot be laundered into apparent strong equality. Adjoint
  `None` versus `Some` remains semantic claim state, while the original token
  remains inspectable only through the attached certified record. Low-level
  identity/receipt aliases remain schema-shaped framing and do not prove the
  opaque helper relationship.

- `color` module (bead qmao.1): the THREE-COLOR epistemic schema —
  `Color::{Verified{lo,hi}, Validated{regime: ValidityDomain, dataset},
  Estimated{estimator, dispersion}}` with the `ColorRank` lattice
  (verified > validated > estimated), the TOTAL conservative pairwise
  `compose` (the result never outranks the weaker operand; verified intervals
  combine per `IntervalOp` with outward-rounded arithmetic; validated regimes
  INTERSECT, and a disjoint intersection demotes further with
  infinite/no-claim dispersion; estimated absorbs everything with additive
  dispersion),
  `check_regime` (validated is a REGIONAL property: exiting, failing to
  report a regime axis, supplying a non-finite state, or declaring an
  empty/non-finite/inverted regime AUTO-DEMOTES to estimated with a
  `Demotion` flag), `regime_demotion` (the borrowed form used for bounded
  multi-parent admission preflight), `verified_from`
  (the only door to a verified color — non-enclosure certificates
  refuse with the laundering teaching error), and `color_of` (the
  honest bridge from existing Evidence receipts: model-free enclosures may be
  `Verified`, while plain `ModelEvidence` is always `Estimated`). Model cards,
  simulated fidelity pairs, discrepancy bands, and declared validity boxes do
  not authenticate experimental membership and cannot mint `Validated`;
  promotion awaits a typed, independently checkable anchored-source receipt.
  `Color::payload_json`
  escapes caller-controlled strings and represents non-finite floats as
  tagged JSON strings, never invalid bare numeric tokens. The distinct
  `Color::canonical_bytes` identity encoding is versioned (v2), structurally
  length-prefixed, deterministically ordered, and preserves every IEEE-754 bit;
  display rounding therefore never aliases color identity or authorization.
  Color algebra v2 also reserves every `derived:` leaf identity and emits
  bounded composition identities below `derived:v2:` so a caller cannot re-root
  computed evidence as an independently anchored source. Readable derived
  identities are emitted only when every component already satisfies the shared
  grammar; invalid/sentinel components use the domain-separated compact form,
  so even empty-regime demotion returns a structurally valid Estimated payload.
  `color_of` recognizes model absence only for the exact
  `ModelEvidence::none()` shape; an empty card list cannot hide discrepancy,
  assumptions, validity restrictions, or an out-of-domain model behind a
  `Verified` numerical interval.
  Write-time enforcement lives HELM-side in fs-ledger over these types.

- `falsify` module (bead qmao.4): bounded declaration catalog and diagnostic
  telemetry — `FalsifierRegistry` refuses a class with zero declarations and
  `standard()` provides seven intended checker families, separating retained
  sampled-interface replay from the much stronger continuum-watertightness
  proposal (certified oriented intersections, winding/degree, and
  coverage-complete subdivision). `catalog_gate` reports missing declarations
  only; it is explicitly not executable or release authority.
  `FalsifierHistory` ingests bounded, source-referencing `FalsifierAttempt`
  values carrying an idempotency ID, class/regime/falsifier, caller-asserted
  claim-revision and retained-artifact references, seed, positive compute
  charge, and typed outcome. It
  rejects undeclared class/falsifier pairs, treats byte-identical retries as
  no-op replays, and rejects conflicting attempt-ID reuse. `doubt` is the
  empirical discrepancy rate widened by a per-row time-uniform union-bound
  Hoeffding heuristic, with cold-start maximum and a never-zero floor.
  Discrepancy attempts emit opaque escaped, fixed-order
  `fs-evidence/falsifier-candidate` schema-version-1
  tombstone/estimator-bug *candidate* projections, correlated by attempt and
  exact seed/compute bit strings; ingestion neither authenticates the caller's
  references nor adjudicates or persists the candidates.
  `allocate_budget` validates
  and max-rescales consequence × doubt × class-review-share weights before
  normalization. `rent_review` is a preliminary class-level diagnostic over
  fixed ordered `RENT_VOLUME`-attempt windows whose zero-discrepancy windows
  decay toward a nonzero floor. Window closure is independent of review-call
  cadence; subthreshold reviews do not erase observations.

## Invariants
1. Conservativeness (G0, evd-001): composed enclosures contain every
   propagation of operand-enclosed true values (300k seeded samples);
   composed validity is exactly the per-parameter intersection;
   assumptions union sorted; discrepancy bands add; `in_domain` is a
   conjunction. Indeterminate IEEE endpoint arithmetic widens to the whole
   real line instead of discarding NaN corners.
2. Severity monotonicity: composition kind = max operand severity, floored
   at Enclosure for float ops; NoClaim absorbs to infinite bounds; an
   Estimate anywhere poisons `certified()` downstream (evd-006).
3. No card, no solver: `ModelRegistry::register_solver` refuses unknown
   cards with teaching text (evd-002).
4. Out-of-distribution discrepancy queries refuse with the violated
   parameter named — never silent extrapolation. Query keys must equal the
   training schema exactly, so an untrained physical dimension cannot be
   supplied and silently ignored; non-finite training or query coordinates are
   unusable, not a way to synthesize or enter a trained box (evd-004).
5. Dominance ties break in declaration order (ModelForm, Statistical,
   Numerical) — deterministic verdicts.
6. Ledger rows and provenance chains are deterministic (repeat-identical).
   Evidence and model-card rows stay valid JSON under hostile metadata and
   tagged non-finite no-claim values; evidence rows retain model assumptions
   and sensitivity headlines rather than dropping those semantic slices.
   Public set-like card, assumption, and known-failure vectors are sorted and
   deduplicated again at the durable rendering boundary, so caller mutation or
   insertion order cannot change row identity.
   Provenance chaining is order-sensitive.
7. Color regime checks fail closed: no empty, inverted, or non-finite regime
   and no non-finite current state can retain `Validated`; disjoint validated
   composition and regime exit both carry infinite/no-claim dispersion.
8. Certified trust boundary (gp3.2.1, evd-012): `Certified<T>` is an
   OPAQUE newtype over `Evidence<T>` — private inner, no `DerefMut`, no
   field access for writing. The ONLY constructor is
   `Evidence::certified()`, which validates the ACTUAL numbers, not the
  constructor that claimed them: scalar evidence requires bit-identical carried
  value and QoI; Exact requires a finite QoI with
   bit-identical bounds; Enclosure requires finite ordered bounds that
   CONTAIN the QoI; statistical e-values, levels, widths, and confidence
   parameters must satisfy their finite domains; model discrepancy must be
   non-negative (positive infinity is the explicit unbounded claim); empty,
   inverted, or non-finite model-validity domains refuse even when a public
   literal asserts `in_domain: true`; Estimate/NoClaim and out-of-domain models
   refuse. Decision breakdown applies the same validity check and assigns an
   infinite model band to an impossible domain.
   Reads flow through `Deref<Target = Evidence<T>>`;
   `Certified::into_evidence()` is the explicit downgrade — the mark is
   lost and any reconstruction must re-enter `certified()`
   (re-validated round trip). Escape hatches are plain `Evidence<T>` or
   `NumericalCertificate::no_claim()`, never a `Certified<T>`. Certification
   requires an owned/`'static` payload so the boundary can detect scalar `f64`
   values and bind them bit-exactly to their QoI; borrowed payloads remain
   plain evidence until promoted to owned values.
   MIGRATION: `Certified<T>` was a type alias for `Evidence<T>`;
   callers that mutated or moved fields of a certified value now call
   `.into_evidence()` first (one workspace site: fs-geom conformance),
   and callers that only read fields are unchanged via Deref. This
   crate has no serializer; persisted evidence re-enters through
   `certified()` on ingest by construction.
9. Decision assessment fails closed (evd-013): malformed or negative
   uncertainty becomes an infinite band; infinite totals and malformed
   thresholds cannot become `DecisionGrade`.
10. Color provenance identities are bounded, byte-length-framed, and
    domain-separated (evd-014). Every generated identity is in the reserved
    `derived:v2:` namespace, while all `derived:` identities are forbidden at
    source admission so computed evidence cannot be re-rooted as a leaf.
11. The shared public color validator and `color_of` bridge fail closed on
    NaN/inverted intervals, malformed identities and regimes, negative/NaN
    dispersion, and malformed evidence/model inputs (evd-015). Plain model
    evidence never self-promotes to `Validated`, including in-domain
    two-fidelity discrepancy evidence (evd-004). Ordered infinite Verified
    endpoints remain sound but vacuous enclosures.
12. `demotion_estimator_identity` is total over arbitrary strings and always
    emits a bounded identity accepted by `color_identity_reason`; invalid
    readable inputs are hash-compacted rather than interpolated. Malformed
    model-card diagnostics bind the complete sorted/deduplicated card set and
    each entry's validation reason in a v2 domain-separated streaming hash;
    distinct residual cards cannot alias merely because the first invalid card
    is the same.
13. Opaque helper-built color-evidence nodes (sj31i.52.2, G0/G3/G4) are
    same-schema replay-stable and publish no partial root on malformed color,
    resource, or cancellation refusal. Source nodes require exactly one
    schema-bound typed source and zero parents. Composition nodes carry no
    source and exactly two descriptor-bound parent rows. Add/Mul/Hull normalize
    parent order before recomputing both the color and identity, without
    deduplicating multiplicity. Helper output agrees with an independent
    canonical construction for every color variant. Exact buffer-reservation
    fault injection and cancellation after an absorbed byte prefix both refuse
    without publishing a root. Raw frame-ID/receipt aliases do not assert these
    helper-level semantic invariants. Legacy `ProvenanceHash` has no conversion
    or rehash bridge into either typed role.
14. Opaque helper-built validity domains (sj31i.52.2, G0/G3/G4) normalize
    insertion order through `BTreeMap`, preserve arbitrary exact UTF-8 axis
    bytes and signed-zero endpoint bits, and change identity for every
    normalized axis or endpoint-bit mutation. Non-finite/inverted bounds,
    row-count/field/chunk/frame resource overflow, and entry or mid-stream
    cancellation refuse without publishing an identity. Helper output agrees
    with independently framed canonical rows.
15. Opaque helper-built certified-f64 semantics (sj31i.52.2, G0/G3/G4) retain
    the exact `Certified<f64>` while agreeing with an independently framed
    strong projection. Every governed scalar, certificate, model, validity,
    sensitivity, or adjoint-presence mutation moves the root; signed-zero and
    sensitivity NaN payload bits remain distinct where their wire types permit
    them. Raw legacy provenance changes and raw `Some` adjoint-token changes do
    not move the root, while `None` versus `Some` does. Non-canonical string
    sets, validity refusal, exact field/row/chunk/frame resource overflow, and
    entry or late cancellation refuse without publishing an identity. The raw
    semantic-ID/receipt aliases do not assert helper-level consistency with an
    attached certified record.

## Error model
Structured teaching errors throughout: `CertifyError`, `RegistryError`,
`OutOfDomain`, `FitError`, `FalsifyError` — all `core::error::Error` with actionable
Display text. Constructors are total (enclosure bounds normalize by
swapping); no panics cross the boundary.

## Determinism class
Deterministic: pure values and mutable diagnostic state machines produce the
same results for the same ordered call sequence; renderings use sorted
(`BTreeMap`) order; there are no clocks, addresses, or hidden randomness.
Bit-stable across runs and platforms up to fs-math-class scalar-arithmetic
divergence.

## Cancellation behavior
Core certificate/color algebra is bounded small synchronous work. Typed color,
validity-domain, and certified-f64 identity helpers accept an explicit
cancellation probe. Color payload copies poll at the configured byte stride;
validity and sensitivity rows poll at stream boundaries; certified-f64 set and
row preflights also poll while traversing caller data. Every refusal consumes
the in-flight encoder and publishes no partial identity. Falsifier allocation
and history review iterate caller data without a `Cx`; allocation length and
distinct history rows are defensively capped, but these diagnostic APIs are not
P7 hot-kernel or cancellation-authoritative paths. Callers must not place large
reviews inside latency-bounded tile loops.

## Unsafe boundary
None. `unsafe_code` denied workspace-wide.

## Feature flags
None. The mechanisms are `[S]`-grade bookkeeping; the models they DESCRIBE
carry their own ambition tags in their cards.

## Conformance tests
tests/conformance.rs, cases evd-001..evd-015 (JSON-line verdicts; seeded
cases carry seeds): the G0 conservativeness battery, the registration
lint, the worked model-discrepancy-dominates example (10% closure vs 0.7%
mesh at a 5% threshold → NotDecisionGrade{ModelForm} + escalation advice,
flipping to DecisionGrade with a 2% calibrated closure), the
out-of-distribution refusal on a synthetic two-fidelity corpus, bracketing
spread reporting with deterministic schema-valid rows, and certification
poisoning. The color cases cover disjoint-regime demotion with infinite
dispersion, outward-rounded verified arithmetic, non-finite state/regime
demotion, and deterministic escaping/tagging of hostile JSON payloads.
In-module suites cover the certificate algebra, validity laws, tie-breaking,
provenance chaining, and card rendering.
`evd-013` exercises the public `Evidence` layer against indeterminate interval
arithmetic and malformed statistical/model uncertainty.
`evd-014` locks bounded collision-resistant derived identities and source
namespace separation; `evd-015` locks the shared color structural invariant
across every public variant and bridge input, including the malformed empty-
regime demotion sentinel and complete-set malformed-card identity stability.
`evd-016` uses `fs-propcheck` seed `0xE71D_4A48_0001` for 512 shrinkable color
composition recipes: all 27 operand-kind/operation combinations, both operand
orders, overlapping and disjoint validated regimes, and structural validation
of every input and output. It locks the no-laundering rank law while preserving
the fixed cases above as regression pins.

- Falsifier registration is total: malformed, duplicate, and empty declaration
  sets refuse; the catalog lint names each distinct undeclared class once in
  canonical sorted order without making a release claim.
- Fixed and adversarial falsifier cases check that budget allocation is
  monotone in consequence AND doubt, with honest
  boundaries (cold-start max, asymptotic doubt floor, dependent-free floor,
  empty-job zero), finite max-rescaled extreme weights, and validation even
  when the requested total budget is zero.
- Every accepted discrepancy attempt produces BOTH correlated pending candidate
  payloads; neither candidate is an active tombstone or an adjudicated
  estimator bug. Exact retries are idempotent, conflicting ID reuse refuses,
  malformed/undeclared attempts leave telemetry unchanged, and two clean rent
  windows decay identically under batched versus incremental review calls.

## Operational V&V artifact schemas (bead frankensim-ext-vv-artifact-schemas-x68z)

The `vv` module is the machine-readable boundary for operational verification,
validation, calibration, and prediction assessment. Its seven top-level schemas
are `ContextOfUse`, `ValidationPlan`, `ExperimentArtifact`, `CalibrationSplit`,
`SolutionVerificationReceipt`, `PredictionAssessment`, and
`AssumptionsLedger`.

- Every individual artifact carries a stable artifact identity and the Five Explicits:
  units, a fixed seed or an explicit not-applicable reason, accuracy/time/memory
  budgets, component/schema versions, and required capabilities. References
  bind both the named artifact and its content identity; a display label is
  never lineage authority. A finite zero-valued header accuracy limit is stored
  as canonical mathematical zero, so IEEE `-0.0` cannot create a second header
  or artifact identity. The complete `VvCase` has its own separately
  versioned identity domain (`org.frankensim.fs-evidence.vv-case.v3`), distinct
  from `…vv-artifact.v3`; this exact case identity is what schema-admission
  receipts bind. The receipt digest itself is governed independently under
  `…vv-schema-admission-receipt.v2`; its v2 preimage includes the admitted
  wire-schema/ruleset versions rather than conflating those versions with the
  receipt-identity era. Each artifact row encodes the artifact family's stable
  one-byte wire tag as well as its slug, id, and content hash; the tag is not
  merely an implementation detail used for sorting.
- `ContextOfUse` binds the decision, named QoIs, acceptance bands, and the
  applicability domain. `ValidationPlan` binds each QoI to its experiment,
  solution-verification, diagnostic, and assessment dependencies. Dependency
  closure is checked per QoI: a missing or mismatched edge refuses, and evidence
  for one QoI cannot silently support another QoI.
- `CalibrationSplit` requires non-empty, pairwise-disjoint calibration,
  validation, and blind-holdout dataset sets. Blind data is inadmissible before
  a non-zero preregistration identity and a release record bind the prior
  commitment and released dataset identity;
  calibration data can never be relabeled as validation or blind evidence.
  The sealed commitment has its own registered v2 identity
  (`…vv-blind-holdout.v2`) over the preregistration hash and canonically ordered,
  length-framed `(observation id, immutable source-locator hash)` rows. Split
  headers and non-blind partitions remain outside this narrow commitment but
  inside the complete split artifact identity.
  Selection minting recomputes the split's canonical content hash and requires
  the exact kind/id/hash reference; a same-id forged reference cannot mint a
  validation or released-blind capability even if later whole-case admission
  is never invoked. Canonical decoding enforces the same invariant locally:
  every validation or blind selection must carry a non-zero split content
  identity, and the split reference enclosing a `BlindHoldout` selection must
  exactly equal
  the release receipt's split reference in kind, id, and content hash before a
  typed selection is returned. A release constructor and decoder also refuse
  all-zero sentinels for the exact split content hash, blind commitment, or
  authority receipt; mutually consistent zero-valued references cannot become
  typed authority. The blind commitment and non-zero release-authority receipt
  remain independently identity-bearing; whole-case admission additionally
  checks the commitment against the admitted split.
- `ExperimentArtifact` records physical dataset lineage, instrument-calibration
  receipts, clock-synchronization evidence and tolerance, repeatability and
  covariance, and an authenticity reference. Covariance rows and columns bind
  an explicit, unique QoI order. Admission first interprets the supplied matrix
  in that declared order, then permutes the tensor into sorted QoI order for
  transport and identity. Simultaneous axis/matrix permutations that pass the
  deterministic covariance screen therefore canonicalize identically, while
  relabeling axes without moving the corresponding matrix entries changes
  tensor semantics. The current finite-precision, unpivoted LDL^T screen runs
  both before and after canonical axis permutation so no admitted artifact can
  fail its own canonical decode. It is deliberately fail-closed, but it is not
  an exact or interval PSD certificate: near-singular tensors can be refused in
  one declared ordering even when another ordering reaches the same canonical
  tensor. Scientific PSD certification and fully presentation-invariant
  admission remain an explicit stronger follow-up, not an implicit claim of
  this structural schema.
  IEEE signed zero is normalized before covariance validation and encoding, so
  `-0.0` cannot mint an identity distinct from mathematical zero.
  Constructor input order is never discarded into a set and silently reused as
  covariance meaning. Every observation row binds an
  `ObservationSourceRef`: the exact dataset source-bytes hash, bounded locator
  domain, positive locator-contract version, non-zero immutable locator hash,
  and non-zero extraction-receipt hash. It also binds the row's QoI,
  instrument, acquisition channel, and clock. Admission requires every row's
  receipt-independent `(dataset, locator domain/version, locator)` identity to
  be unique: changing only extraction evidence cannot manufacture a second raw
  observation or move one locator across calibration/validation partitions.
  The receipt remains part of the complete typed source and artifact identity.
  Admission also requires every row's
  dataset hash to equal the experiment's exact source-bytes hash; both that
  source hash and the custody-receipt hash must be non-zero. It further requires
  exact manifest/declaration QoI equality, exactly one current non-zero
  calibration receipt for each referenced instrument, and membership of each
  referenced clock in the declared synchronization topology. The schema checks
  presence, dimensions, finite values, typed reference closure, and these
  internal bindings; it does not itself authenticate a laboratory, instrument,
  clock, signer, extraction receipt, or dataset.
- Raw-lineage carriers and their top-level case/receipt wrappers use deliberately
  redacted `Debug` implementations. Formatting observation sources/manifests,
  experiments, authenticity records, calibration splits, observation
  selections, blind releases, `VvCase`, `SchemaAdmissionReceipt`, or
  `AdmittedVvCase`—and formatting the `VvArtifact` sum wrapper—may expose
  bounded counts, stage/variant tags, contract/schema versions, and
  authentication/binding-presence booleans, but never raw row
  membership, dataset/locator/extraction/custody hashes, metrology bindings,
  preregistration material, artifact maps, or release-authority material. This
  is defense against accidental log and panic disclosure, not an authorization
  boundary: explicit typed getters still expose the exact values to code that
  already holds the evidence object.
- A validation plan must name observability, identifiability, confounding, and
  inverse-crime diagnostics. A synthetic oracle, high-fidelity code, reduced
  model, or second implementation is code/solution verification or discrepancy
  evidence. It cannot mint physical validation or a `Validated` color without
  an admitted physical experiment.
- `SolutionVerificationReceipt` reports the QoI-specific mesh, time,
  nonlinear-solve, and iterative-solve uncertainty components. Its reported
  numerical envelope may not understate their conservative composition, and a
  code-comparison receipt is kept distinct from the numerical enclosure.
- The uncertainty waterfall has exactly the separately named model-form,
  parameter, numerical, data, aleatory, and epistemic sources. It is either a
  conservative bound composition or an explicitly probabilistic composition
  with confidence/dependence semantics and supporting independence evidence;
  neither mode may omit, duplicate, or silently average sources.
- `PredictionAssessment` binds validation metrics to both experimental and
  numerical uncertainty, posterior-predictive checks, the QoI waterfall, and
  an applicability decision. Its seven evidence axes are categorical: code
  verification, solution verification, numerical uncertainty, parameter/data
  uncertainty, model-form validation, prediction-domain relevance, and
  comparison to experiment. Colors are never converted to percentages or
  averaged into a numeric confidence score.
- Applicability is evaluated against the declared domain. Domain exit, a
  missing axis, or a non-finite state follows the declared `Demote` or `Refuse`
  policy and is retained in the assessment; silent extrapolation is forbidden.
  Process-standard conformance receipts are recorded separately and cannot
  satisfy model-validation or experiment-comparison axes.
- `AssumptionsLedger` uses `AssumptionId | predicate | scope | evidence |
  runtime monitor | violation effect | owner | expiry/review gate` and seeds the
  plan's runtime assumptions: A-001 reduced/rigid-body adequacy; A-002 MQS
  regime; A-003 mixed-cylinder/section-averaged-duct adequacy; A-004 smooth
  contact and continuum-scale adequacy; A-005 symmetry preservation; A-006
  material/process/query validity; A-007 closure/correlation/turbulence/
  lubrication applicability; and A-008 probability/dependence/population
  representativeness. Missing monitors or violated predicates demote, escalate,
  or refuse according to the row; they never disappear from the lineage.

Preregistered validation metrics are ENFORCED, not merely listed (bead
gt1k3): every `QoiValidationPlan` spec is evaluated against the typed
artifact numbers — interval agreement derives from
`|observed - predicted|` vs the combined uncertainty, normalized
discrepancy from the preregistered maximum (zero combined uncertainty
admits only exact agreement), and every posterior-predictive check must
carry EXACTLY the preregistered minimum tail probability and pass it.
Outcomes are derived, never caller-asserted; failed outcomes stay
visible as violations and additionally refuse any POSITIVE
model-form-validation or comparison-to-experiment evidence axis.

Observation identity and interpretation are authority-bound (beads xl3yi and
i94v.3.3.1, schema v3): every experiment carries a canonical
`ObservationManifest` mapping each observation id to a typed
`ObservationSourceRef` plus the exact QoI, instrument, acquisition channel, and
clock identities used to interpret it. Each source reference binds exact
dataset bytes, locator domain and positive contract version, non-zero locator,
and non-zero extraction receipt; experiment admission proves the row's dataset
hash equals its `DataAuthenticity` source-bytes hash. Receipt-independent raw
locator identities are INJECTIVE, so two ids can never alias one immutable
source row by relabelling only its extraction receipt, while genuinely distinct
locators with equal values remain distinct. The complete source, including its
receipt, remains identity-bearing.
Every one of those bindings enters the derived aggregate observations hash
under domain `…vv-observation-manifest.v3`; none is transported as a
caller-supplied aggregate. Blind-holdout partitions retain their hash-only
`(id, source)` bindings and sealed v2 commitment (domain
`…vv-blind-holdout.v2`), while case closure cross-checks each source against the
richer experiment manifest in addition to the exact id-coverage rule. The wire
schema and artifact content-identity domain deliberately advance together to
v3 (`…vv-artifact.v3`); v1/v2 artifacts do not decode under v3 and cannot share
content identities with current artifacts.

Two no-claim boundaries remain explicit in schema v3. First,
`ObservationLocatorIdentity` is receipt-independent but is still scoped by its
caller-declared locator domain and contract version. It does **not** normalize
equivalent raw records expressed under different locator domains or versions,
so it cannot by itself prevent cross-contract relabelling. Promotion work must
add a stable schema-owned `RawRecordIdentity`, recomputable from an exact
dataset-scoped byte-span or canonical Merkle-leaf proof, while retaining the
locator contract as separate interpretive authority. A supplied leaf hash is
not authentication; typed verification must check the byte-span/inclusion
proof against the admitted dataset bytes.

Second, a `VvCase` does not yet prove disjoint raw-locator use across two
different `ExperimentArtifact` wrappers and does not yet encode explicit
cross-experiment sharing or joint-likelihood groups. A blanket overlap ban
would incorrectly reject legitimate shared-data analyses. Successor work must
make the policy explicit, require exact group coverage and likelihood
authority, and test both undeclared reuse refusal and declared reuse admission.

V&V artifacts have a versioned canonical bounded binary encoding. Fields are
fixed-order and length-framed; maps/sets use canonical ordering; floating-point
values preserve their exact IEEE-754 bits except at schema-declared
mathematical-zero seams (currently header accuracy and covariance entries),
where `-0.0` is normalized to `+0.0`. Decoding caps total bytes, counts,
nesting, and string sizes, and refuses unknown tags, invalid UTF-8, duplicate or
out-of-order keys, non-canonical encodings, and trailing bytes. An accepted
decode must re-encode byte-identically. The 4 MiB transport ceiling is exact:
an N-byte encoder/decoder input is accepted by the size gate and N+1 is refused
with the canonical rule and byte offset. Caller-retained applicability decisions
must also match recomputed numeric violations bit-for-bit, and signed-zero
metric thresholds share one semantic ordering key so aliases cannot evade
duplicate detection. The admission receipt binds the schema version, dedicated
complete-case content identity, complete individual-artifact identity map,
QoI/context identity, and validation-rule version under the registered
`…vv-schema-admission-receipt.v2` identity; it proves that those
structural rules ran, not that referenced
experimental evidence is authentic or scientifically sufficient.
Receipt artifact-map rows are canonically sorted by the artifact family's
explicit stable wire tag and then artifact id. They never inherit ordering from
Rust enum declaration position or caller/map insertion order; the explicit tag
mapping and artifact-kind variants are governed inputs to the receipt identity
schema.

Stable rule-coded errors identify the failed invariant and artifact path. G0
tests cover all seven schema round trips, byte/identity determinism, split
disjointness, numerical/waterfall arithmetic, rule-code stability, and
applicability demotion/refusal. G3 cases mutate dataset roles, dependency QoIs,
comparison source kinds, evidence-axis order, domain coordinates, and binary
framing to demonstrate fail-closed behavior. G5 repeats canonical identity
construction across insertion orders. These tests prove schema mechanics only:
they do not prove experimental authenticity, independence, model adequacy,
physical validation, process-standard conformance, or decision fitness.

## Admitted scientific color (bead 6pf9, stage S1)

- `Color` is a DECLARATION: publicly constructible, structurally validated,
  never authority. `AdmittedColor` is the opaque positive-evidence handle:
  private fields, single constructor `AdmittedColor::from_receipt(color,
  receipt, verifier)`.
- Local gates fire before the injected capability: malformed payloads
  (`validate_color_payload`), non-positive ranks (Estimated), and
  stale-algebra receipts (`color_algebra_version !=
  COLOR_ALGEBRA_VERSION`) refuse even under an accept-everything verifier.
- An accepting `AdmissionDecision` must name exactly the policy fingerprint
  committed by the receipt. A mismatch returns structured `PolicyMismatch`
  rather than silently detaching the admitted value from its lineage; retained
  tests cover both a matching authority and an accepting policy substitution.
- `AdmissionReceipt` is plain data (node provenance hash, row schema
  version, algebra version, policy fingerprint). Authority lives in the
  `AdmissionVerifier` capability; the default `NoAdmissionVerifier` is
  deny-all, so at this layer NOTHING admits. The authenticating verifier is
  HELM-side (`fs_ledger::LedgerColorAdmissionVerifier`), keeping this crate
  free of upward dependencies.
- No-claim: `AdmittedColor` is capability-gated, not cryptographically
  unforgeable. A lying verifier at the composition root can admit anything —
  the same trust model as `WaiverVerifier` and `SourceOriginVerifier`, and
  visible at the same audit surface. Stage S1 is additive: positive-evidence
  consumer APIs migrate to require `AdmittedColor` in later stages of bead
  6pf9.

## No-claim boundaries (colors)

- Verified-interval composition here covers outward-rounded Add/Mul and exact
  endpoint Hull; the full ledger operation algebra composes through fs-ivl
  when wired.
- Estimated dispersion combines additively (conservative); calibrated
  dispersion algebra joins the color-probes bead.
- `ModelEvidence` carries model-form diagnostics, not experimental authority.
  Even an in-domain discrepancy model trained on paired simulations remains
  `Estimated`; this crate has no authenticated anchor type that can admit a new
  `Validated` leaf.
- Estimator identity chains are human-readable up to
  `MAX_COLOR_IDENTITY_BYTES` and then collapse to a domain-separated BLAKE3
  composition identity. Human-readable pairs are byte-length-framed, so an
  identity containing `+` or `@` cannot alias a different composition tree.
  Every binary composition retains both framed operands, including repeated
  identities, and mixed `Verified`/`Validated`/`Estimated` operand classes use
  distinct domain labels. A source literally named `verified` therefore cannot
  impersonate a Verified operand. Pass-through operands receive a bounded
  `derived:v2:` identity rather than retaining a source-leaf identity. Single
  model-card identities remain source identities until an operation composes
  them.
  This bounds provenance payload growth without changing the Estimated rank or
  its accumulated dispersion. Compact identities are deterministic provenance
  labels, not signatures or source authority.

## No-claim boundaries
- Statistical composition is CONSERVATIVE-WEAKEST v1 (half-widths add,
  confidences min, mixed kinds keep the width-bearing certificate);
  proper e-value arithmetic (products under independence, e-BH) is
  fs-eproc's contract and will replace this composition rule behind the
  same API.
- Discrepancy models are honest BOOKKEEPING (observed box + mean/max
  band), not learning; trained/learned discrepancy models (FrankenTorch)
  arrive with fs-surrogate and will implement the same query/refusal
  surface.
- The adjoint hook is carried, never composed here — composed tapes are
  fs-ad's contract.
- First-order band addition is conservative for small relative bands; it
  is NOT a rigorous product-form bound for large ones — fs-ivl composition
  should be used for the numerical slice when bands are large (the
  numerical slice already does interval arithmetic; the TOTAL across
  sources is the first-order sum).
- Ledger persistence is row RENDERING only; the `evidence` /
  `model_cards` tables land with fs-ledger (rows are shaped for that
  migration).
- `ProvenanceHash` is FNV-1a until the BLAKE3-class ledger hash supersedes
  it (same upgrade path as fs-obs). It remains legacy deterministic
  correlation only and cannot construct `ColorEvidenceSourceIdV1` or
  `ColorEvidenceNodeIdV1`.
- `IdentifiedValidityDomainV1` proves helper validation and exact normalized
  declared bounds. A raw `ValidityDomainIdV1`/receipt proves only schema-shaped
  framing. Neither proves that observations occupy the box, that a model is
  valid there, or that any source/calibration authority admitted the
  declaration.
- `IdentifiedCertifiedF64EvidenceV1` proves only the helper-defined strong
  semantic projection of an already-local `Certified<f64>` and keeps that
  record attached. It does not add units, a quantity kind, source or model-card
  content identity, model-card version/calibration/known-failure binding,
  gradient authority, seeds, budgets, versions, capabilities, or external
  trust. The statistical fields bind only the numeric local variant payload;
  `Evidence` currently carries no null/hypothesis/estimand, method, sample or
  dataset, or dependence-context identity, so none is implied or bound here.
  Model-card names are identifiers only. The projection is not a scientific
  certificate, an origin signature, a ledger admission, or a commutativity
  claim. Raw legacy FNV provenance and adjoint values remain inspectable
  correlation metadata but are intentionally outside the strong root; only
  adjoint presence is bound. A raw
  `CertifiedF64EvidenceIdV1`/receipt is merely schema-shaped framing and can
  encode a frame that disagrees with any purported attached record.
- A successful opaque helper build proves canonical source framing or exact
  replay of the named Add/Mul/Hull color operation. It proves
  nothing about source origin, experimental membership, model correctness, or
  admission. External trust must travel via an independent `AuthorityRef`;
  scientific rank remains in `Color`, numeric structural consistency remains in
  `Certified<T>`, and color admission remains in `AdmittedColor`. These axes do
  not promote one another.
- G5 promotion is intentionally pending a retained known-answer vector replayed
  across the reference ISAs and thread counts. Same-process recomputation and
  independent G3 construction are not reported as that stronger proof.

## No-claim boundaries (falsifiers)

- The registry stores falsifier IDENTITIES and stated methods; executing
  a falsifier (running rays, FD probes, full solves) is each consumer
  kernel's code. A public method string cannot prove executable binding or
  checker independence, and `catalog_gate` cannot authorize release.
- Catalog class/spec identifiers are bounded ASCII slugs; human prose belongs
  in the separately bounded method/detail fields. Class count, declarations per
  class, per-class and total declaration bytes, lint request/output size,
  allocation length, history rows/key bytes, attempt count, and retained
  attempt bytes all have explicit ceilings. These are defensive synchronous
  caps, not a substitute for a budgeted `Cx` API.
- Exact-instance authority requires a typed retained receipt bound to the claim,
  specification, policy, implementation/TCB, seed, stopping rule, budget,
  outcome, and artifacts, admitted by fs-package/fs-checker. That successor
  architecture is not implemented by this telemetry module.
- `consequence` is supplied by the caller (ledger-DAG dependent-weight
  traversal is HELM-side); the allocator's contract is what it does with
  the number, floors included.
- Candidate payloads are opaque escaped
  `fs-evidence/falsifier-candidate` schema-version-1 fixed-order JSON
  projections, with the full-width seed and compute charge represented as exact
  hexadecimal bit strings; they are not authenticated content addresses.
  History mutation is not transactional with an external
  sink; downstream adjudication and atomic persistence are mandatory before a
  candidate can invalidate a claim or count as an estimator defect.
- Stable attempt IDs prevent identical retry/replay from double-counting and
  conflicting reuse fails closed. They remain caller assertions rather than
  authenticated ledger identities, so history still cannot be trusted against
  fabricated or selectively omitted attempts. `doubt` assumes
  an admitted stationary Bernoulli process and controls time multiplicity only
  within one row, not the family of class/regime rows; here it is planning
  telemetry only. Candidate discrepancy reports conservatively raise this
  telemetry but cannot activate authority.
- Rent decay is per-class rather than per-falsifier. Fixed-volume windows are
  invariant to review-call cadence, but are not authenticated calendar/policy
  windows and have no restoration policy. Per-falsifier
  receipts, independence threat graphs, cadence enforcement, and restoration
  policy remain successor work (xpck.6).

## Conservation-defect accounting (bead frankensim-leapfrog-2026-program-i94v.1.3.1, I04.1)

`balance` defines the typed accounting language for balance defects:
`BalanceDefectReceipt` pins quantity kind (with per-kind admissible
accounting rule and canonical unit), extensive-versus-rate semantics,
sign/orientation, chain/cochain/region spatial support on one named complex,
instant/window time support on one named logical clock, producer identity,
an exact `ChartRef` (id, version, content digest), owned per-account terms
with `ColorRank` evidence, an expected-closure interval, a five-state
`DefectState` (exact zero / bounded / unknown / unowned remainder /
inapplicable — distinct, never collapsed), split numerical/material/model
uncertainty, and content-addressed lineage. Receipts compose over provably
disjoint spatial partitions and adjacent extensive windows only; every
incompatibility (chart, clock, sign, scale, quantity, semantics, support,
window adjacency, duplicate ownership, accounting-rule violation) refuses
with a stable named `BalanceRule` slug. Canonical transport is bounded,
version-gated, and canonical (decode re-encodes bit-for-bit or refuses);
identity is the domain-separated BLAKE3 hash of the canonical bytes.

### No-claim boundaries (conservation-defect accounting)

- The module accounts for defect values it is HANDED; it does not detect,
  localize, or attribute physical events (those are I04.2+ scopes).
- No unit-scale conversion: `scale_pow10` is identity metadata and scale
  mismatches refuse; cross-scale rebasing is a consumer decision.
- Split uncertainties add linearly (perfect-correlation conservative bound);
  no independence, distribution, or coverage-probability claim.
- Interval sums widen one ulp outward per composition; enclosures are
  conservative, not tight.
- Mass accounting is non-relativistic; species ledgers are open (reactions
  may produce/consume); entropy/exergy sign restrictions apply only to
  chart accounts with the matching declared role.
- `ColorRank` is carried and composed weakest-wins, but no payload-level
  color algebra (regimes, datasets, estimators) is evaluated here.

## Evidence-action vocabulary (bead frankensim-leapfrog-2026-program-i94v.1.4.1, I08.1)

`action` defines the decision-facing vocabulary for buying uncertainty
reduction: `ActionProposal` pins a closed ten-kind taxonomy (solver
tolerance, mesh/time refinement, representation escalation, UQ samples,
material/coupon test, sensor campaign, falsification, standards obligation,
explicit refusal — coupon tests and sensor campaigns are the physical
kinds), one targeted claim slice (claim id × six-component uncertainty
decomposition), a planning-only expected-response factor in (0, 1] or
explicit Unknown, per-axis cost states (money/compute/memory/lead-time,
each a `NonNegRange` envelope in one unit or `Unknown` with a NAMED
authority gap — never a silent zero), capability demands, dependencies on
content ids, exclusivity and correlation group identities, optional expiry
on a named clock, an evidence-color CEILING, and the proposer.
`Portfolio::admit` refuses duplicate content identity (idempotent
proposals), missing/self dependencies and cycles (content addressing makes
true cycles inexpressible), exclusivity violations, cross-clock expiry, and
expired proposals; per-axis totals are `Known` only when every contribution
shares one unit, else `Incomparable` listing every gap and unit — cross-
currency sums refuse because price normalization is a policy input. The
planning/execution split is structural: proposals cannot carry outcome
colors; only `ExecutionReceipt::admit` (bound to the exact proposal id,
outcome never above the ceiling) can, so planned physical tests cannot
raise evidence until executed. Canonical transport is bounded,
version-gated, canonical (re-encode bit-for-bit or refuse), with
domain-separated BLAKE3 identity.

### No-claim boundaries (evidence-action vocabulary)

- Cost and response "distributions" are interval envelopes; no parametric,
  moment, or coverage claim.
- Correlation is a group identity only; no joint-distribution claim.
- Portfolio admission proves feasibility structure, not optimality; no
  decision, ranking, or value-of-information computation lives here.
- No currency, time, or unit conversion — ever; mismatches refuse and
  normalization decisions belong to policy layers.
- An `ExecutionReceipt` binds outcome color to a proposal identity but does
  not itself verify the execution happened; execution attestation is
  ledger/checker scope.
