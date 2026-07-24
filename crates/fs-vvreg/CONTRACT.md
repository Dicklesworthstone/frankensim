# CONTRACT: fs-vvreg

The Gauntlet G1/G2 benchmark & V&V registry (bead
frankensim-ext-benchmark-vv-registry-f1gq, epic E0c): the single place
where a benchmark family name becomes an executable claim target — or is
refused until it is one.

## Purpose and layer

Layer UTIL (versioned registry data + fail-closed citation, corpus-query, and
standards-source gates). Depends on `fs-blake3` (domain-separated content
identity), `fs-evidence` (`Evidence<T>` plus claim-color caps), and `fs-qty`
(runtime dimensions for loaded measurement and context data). A family name (TEAM, NAFEMS,
CFR, IFToMM, ECN) is NOT an executable benchmark: every G1/G2 entry needs
exact version/edition, source, license, input-deck identity, oracle binding,
QoIs, and acceptance envelopes before any solver claims against it.

The `standards` module is the metadata-only source boundary for
standards-derived rules. It represents an exact standard, part, edition,
ordered amendment or corrigendum chain, jurisdiction/profile, lifecycle state,
external locator, source hash, license/access policy, supersession edge, and
explicit reference state without accepting or storing a protected standards
body.

## Public types and semantics

- `RegistryEntry` — one row: id, tier (`RegistryTier::G1Analytic` /
  `G2Benchmark`), family, title, `Edition`, source, `LicenseState`,
  `DeckPin`, `OracleBinding`, QoIs (`Qoi` + `AcceptanceEnvelope`), notes.
  Incomplete rows may EXIST (recording a known target is honest) but
  refuse citation.
- `Registry::check_acceptance_envelope` — arithmetic-only executable QoI gate
  on the seeded registry. Callers name an entry and QoI; the gate uniquely
  selects the stored unit/envelope, refuses caller-built registries, and binds
  the entry + registry digests. A tolerance observation supplies independent
  reference + computed values and uses
  `abs(computed-reference) <= atol + rtol*abs(reference)`; an interval
  observation supplies the computed value and uses inclusive `[lo, hi]`.
  Passing calls return a sealed `EnvelopeVerdict`; violations retain the same
  full verdict, while pre-verdict arithmetic refusals retain a sealed,
  replay-complete `EnvelopeAttempt`. All diagnostics have canonical exact-bit
  JSON. Modes cannot be mixed, non-finite inputs and derived overflow fail
  closed, and zero signed margin passes exactly on the boundary.
- `validate_entry` — the citation gates as a validation-only probe: it
  returns `Ok(())` or a typed `CitationRefusal` naming the first failing
  gate, in documented order: id shape/size, QoI-count cap, blank text
  fields (family, title, source, notes, QoI names/units), edition, license,
  deck, oracle binding, QoI presence, duplicate QoI names, per-QoI
  envelope pin/validity. It can never mint a receipt.
- `Registry::cite` — the ONLY receipt-minting path. It refuses
  caller-built registries (`UnauthoritativeRegistry`), refuses ambiguous
  (duplicated) ids rather than picking one of the conflicting rows, runs
  the full gate chain, and binds the resulting `CitationReceipt` to the
  seeded registry's content digest.
- `OracleBinding` — `Unpinned` targets have no pinned oracle identity or
  comparison procedure; `SelfContained` decks carry their complete closed
  form/procedure; `DerivationRequired` decks deliberately delegate a
  load-bearing derivation and stay NON-CITABLE (typed `UnboundOracle`
  refusal) until a derivation receipt mechanism binds the obligation.
- `CitationReceipt` — SEALED: private fields, no public constructor;
  holding one proves admission ran. Accessors expose entry id, tier,
  exact edition, deck digest, entry digest, registry digest, and registry
  version. Color rule lives here: `numerical_claim_cap()` is at most
  `Verified` for the exact edition and scope; `physical_claim_cap()` is
  unconditionally `Estimated` in this slice (the `Validated` upgrade
  requires a typed held-out-evidence binding, not a caller-asserted flag).
  No color is inherited from a publisher's name.
- `ConsumptionStatus` / `ConsumptionRecord` — Appendix-D discipline:
  consuming beads record unread/read/derived/reproduced/
  independently_falsified and pin the exact artifact version (the entry
  digest).
- `PrimaryReference` — the 30-entry seed of definition/provenance anchors.
  Deliberately mints no color and exposes no authority path.
- `Registry` — sorted rows + references, `lint()` (citability partition +
  seed-integrity findings: duplicate ids/keys/indices — including same key
  at different indices — and blank reference fields; duplicated ids are
  never citable), `canonical_rows()` (deterministic serialization; floats
  as IEEE-754 bit tokens), `digest()` (domain-separated, length-framed
  BLAKE3 identity), `entry_digest`.
- `registry()` — the seeded workspace registry: 12 citable G1 analytic
  entries (authored canonical specs), 6 derivation-required G1 targets
  (Geneva, Atkinson, Bennett mobility, isentropic nozzle, Sod, Lax —
  registered, non-citable until their delegated oracle/deck content is
  pinned), and 15 deliberately unpinned G2 targets.
- `standards::StandardEditionKey` — exact `(standard_id, part, edition)`
  identity. `part=None` means the standard is not partitioned; it never means
  "choose any part." Lookups never inherit support from another edition.
- `standards::StandardSourceDraft` / `StandardManifest` — untrusted metadata
  input and its validated, sorted, sealed form. Admission rejects exact-key
  collisions, duplicate typed changes, missing/self/cyclic supersession
  targets, zero source hashes, invalid text, and explicit resource-cap
  violations. The canonical `FSMF` schema is versioned, length-framed, and
  rejects semantic-but-noncanonical encodings on decode.
- `standards::ProtectedTextReference` — publisher/repository catalog plus an
  external locator only. There is deliberately no protected-text/body field;
  authorized bytes remain out of band and are presented only as an observed
  hash at the binding gate.
- `standards::RuleBindingRequest` / `RuleProvenance` — untrusted request and
  sealed exact provenance for a derived rule. Binding requires an exact
  admitted edition, a nonzero pinned source hash matching the observed bytes,
  currently available access, a nonzero derived-rule hash, and an explicit
  read/derived/reproduced state. Historical editions additionally require
  `RuleUsePolicy::HistoricalReplay`.
- `corpus::DatasetDraft` / `admit_dataset` / `CorpusDataset` — the version-3
  validation-dataset boundary. Fifteen top-level fields are mandatory and
  missing fields return `CorpusError::MissingField { field: DatasetField }`:
  id, title, immutable earliest-retained payload, sensor roster, geometry, acquisition
  environment, partition, preprocessing lineage, final artifact,
  context-of-use ranges, license/redistribution terms, acquisition provenance,
  retention policy, acceptance-envelope records, and a legacy A-E evidence
  tag interpreted as non-ranked portfolio coordinates.
  Sensors bind raw channel and measurement dimensions. Instrument identity,
  calibration, placement, geometry, acquisition environment/date, license, and
  complete preprocessing lineage each distinguish `Available` from a reason-bearing
  unavailable/unreplayable state. Those states are mandatory to record but
  cap use at `Estimated`; they are never replaced with invented metadata.
  Measurement uncertainty is bounded, covariance-diagonal, or explicit
  `Unstated`; covariance carries squared `fs_qty::Dims`.
- `partition::PartitionLedger` / `CorpusRegistry::query` — the
  seeded-registry-only, purpose-typed evidence boundary. The request names
  calibration, validation, or blind-evaluation purpose independently of the
  governed partition and supplies every context coordinate with matching
  dimensions and an in-range finite value. Training and calibration partitions
  both admit only calibration use. Validation admits only validation use, and
  blind holdout additionally requires a nonzero preregistration plus blind-
  manifest release bound to the current generation. Missing, duplicate,
  unknown, mismatched, and out-of-range coordinates remain distinct typed
  refusals. Success returns
  `fs_evidence::Evidence<&CorpusDataset>` with a numerical `NoClaim`, unknown
  statistical width, unbounded model-form discrepancy, and the declared
  context box copied into model validity, wrapped with a sealed access receipt
  binding dataset content, generation, partition, purpose, exact canonical
  query context, preceding repartition event, and blind release; the
  dataset's separately exposed physical claim cap is `Validated` only for C/D
  data carrying the controlled-experimental-validation coordinate, with
  original raw retention, stated uncertainty, complete instrument,
  calibration, placement, geometry, environment, acquisition-date, license,
  replayable lineage, and pinned-envelope authority. Every explicit gap demotes
  the cap to `Estimated`. Level E field monitoring alone is always
  `Estimated`.
- `partition::{RepartitionReceipt,BlindReleaseReceipt}` — canonical versioned
  state-transition records. Repartition increments a per-dataset generation,
  chains the previous event identity, clears any blind release, and explicitly
  records whether an evaluation-to-calibration move stales earlier validation
  claims. A blind release is idempotent only for identical inputs; conflicting
  replacements refuse.
- `partition::PartitionLedger::register_model` / `ModelTaint` — the only public
  model-taint constructor. It requires fresh calibration-purpose access
  receipts and computes a bounded, order-independent transitive closure over
  direct datasets and parent model artifacts. Each source retains a
  deterministic artifact path to the direct consumer.
- `partition::PartitionLedger::validate_model` / `TaintValidationReceipt` —
  requires fresh validation or released-blind accesses and refuses any exact
  dataset-content intersection with the model's full calibration closure. A
  success receipt binds the exact model-taint identity and every distinct
  access identity; it is a disjointness record, not scientific authority.
- `partition::PartitionReceiptRecord::{encode,decode}` — one closed `FSVVREC`
  version-1 wire envelope for dataset-access, repartition, blind-release,
  model-taint, and validation receipts. Decode bounds every string, source
  count, access count, model path, and total record; validates closed tags and
  ordering; re-derives the existing semantic identity from exact fields; and
  refuses future versions, truncation, extension, malformed UTF-8, semantic
  inconsistency, identity substitution, and noncanonical bytes. The envelope
  is stored under generic fs-ledger artifact kind
  `vv-partition-receipt-v1`; wire transport does not mint ledger authority.
- `corpus::CorpusDataset::{encode,decode,digest}` — bounded canonical `FSVVCRP`
  version-3 binary round trip and domain-separated content identity. Version 3
  keeps the v2 `BlindHoldout` partition tag but rotates dataset and registry
  identity domains because A-E are coordinate labels rather than ranks and
  field evidence no longer inherits a validation cap. Decode
  validates, canonicalizes, and byte-compares; noncanonical, future-version,
  truncated, trailing, invalid-tag, and invalid-UTF-8 inputs refuse.
- `corpus::CorpusRegistry::audit` / binary `corpus-audit` — deterministic
  per-dataset mandatory/optional completeness table followed by the fixed-scope
  seven-axis cooling-QoI coverage map. Every QoI retains a count for numerical
  verification, cross-code agreement, controlled experiments, blind
  prediction, field monitoring, transferability, and independent reproduction;
  no maximum or average is computed. Zero controlled-experiment counts emit
  structured `level=WARN qoi_gap=... evidence_axis=...` diagnostics for
  downstream planning. Optional as-built geometry and calibration-valid-through gaps,
  plus mandatory reason-bearing no-claim states, produce structured
  `level=WARN` rows; validation defects produce `level=ERROR` and a nonzero CLI
  exit.
- `portfolio::{EvidenceAxis,EvidencePortfolio,PortfolioClaimClass,
  PortfolioAdmission}` — the schema-v1 claim-scoped portfolio algebra. The
  seven axes are independent categorical coordinates, not an epistemic order.
  Legacy A-E tags map exactly to coordinates: A numerical verification; B
  cross-code agreement; C controlled experiment; D controlled experiment plus
  blind prediction; E field monitoring. Claim classes name non-substitutable
  required axes and admission matches exact QoI and regime. Exact duplicate
  observations are idempotent. Independent reproduction additionally requires
  a source and group both distinct from the controlled experiment. Portfolio
  and admission identities use versioned, domain-separated canonical
  encodings; the opaque admission proves only that this structural rule ran.
- `adversarial::{AdversarialRegistry,AdversarialCase,AdversarialOutcome,
  AdversarialAssessment}` — the schema-v1 honesty registry for thermal
  challenge regimes. Eight cases attack attached-flow, known-contact,
  convection-dominance, stable-fan, known-leakage, fixed-property,
  lumped-temperature, and forced-convection assumptions. Each names one
  dominant uncertainty and distinguishes an exact retained corpus binding
  from planned evidence with a reason and tracking Bead. An accepted
  prediction passes only inside its inclusive challenge envelope; an
  out-of-envelope accepted prediction is counted as a false acceptance.
  Refusal or demotion passes only when it attributes the registered dominant
  uncertainty. `render_regime_limitations` emits every case in canonical order
  for the future public scorecard and renders missing execution or unretained
  planned evidence as `NO-DATA`, never zero.
- `thermal_level_a::thermal_level_a_cases` — 19 frozen Level-A thermal
  definitions backed by one retained TSV manifest: 12 analytic values across
  planar/2-D/axisymmetric/spherical conduction, fin efficiency, a
  lumped-capacitance limit, fully developed duct Nusselt limits, a parallel-
  plate view factor, and series contact resistance; plus seven two-sided G1
  targets spanning P1/P2 primal and adjoint order, mixed Neumann/Robin
  boundaries, and anisotropic nonlinear conductivity. Each case declares
  dimensions, formula semantics, context, and an acceptance rule.
- `thermal_level_b::thermal_level_b_cases` — 4 typed cross-code case
  definitions (convective block with a bilinear source, declared rotated-frame
  tensor, element-mean `k(T)` slab, 3-D film fin) bound to one retained
  external-reference TSV manifest frozen by the pinned `tools/vvref`
  scikit-fem deck runner. `parse_thermal_level_b_manifest` is fail-closed;
  `verify_spec_echo` requires every load-bearing echoed spec number to be
  BIT-IDENTICAL to the catalog constants; `verify_probe_grid` recomputes
  probe positions with the fixture arithmetic `index * (extent / count)` and
  requires bit-identity, witnessing cross-language Kuhn-mesh parity without
  running the external stack; `thermal_level_b_deck_bytes` exposes committed
  deck bytes whose BLAKE3 must match the manifest's recorded deck hash. The
  external run's own known-answer self-check must have recorded `pass` or
  parsing refuses. Declared per-case probe-temperature envelopes are
  investigation triggers, never silently widened.
- `scorecard::build_scorecard` — the public V&V scorecard as a pure,
  fail-closed projection of (corpus registry, adversarial registry,
  ledgered `ScorecardRunRecord`s, executed `AdversarialAssessment`s). Cells
  are keyed by (QoI metric, canonical regime label derived from the
  acceptance record's `ContextRange`s); each cell reports reference counts
  per portfolio axis, external-reference coverage (`EXTERNAL_AXES` =
  cross-code, controlled-experimental, blind-predictive, field), sorted
  per-run signed-error rows with the declared `ReferenceUncertainty`
  rendered next to the model error, envelope arithmetic
  (`EnvelopeCheck::{Pass, Fail, Unpinned}`), interval coverage (always
  `NO-DATA` until e07 machinery exists), and a `FalseAcceptanceCell` that
  is `NoData` whenever zero bound challenges executed. Run records naming
  unknown datasets or undeclared metrics refuse the whole build; foreign
  or duplicate assessments refuse through the adversarial renderer's own
  validation. `render_markdown`/`render_json` are pure functions of the
  typed content; the scorecard identity binds typed fields (never rendered
  bytes) plus sorted assessment identities. The honesty header carries
  corpus schema/digest, seeded-vs-caller-built authority, per-level and
  per-axis dataset counts, and a sorted KNOWN GAPS list naming every
  (QoI, regime) cell — plus every scoped cooling QoI with no dataset at
  all — that has zero external coverage.

## Invariants

- FAIL-CLOSED CITATION: an entry missing any load-bearing field (edition,
  license, deck hash, oracle, QoI, envelope) cannot be cited; the refusal
  is typed and names the field. An unpinned family name never acts as an oracle.
  Ambiguous ids and duplicate QoI names refuse; a deck that delegates its
  oracle refuses.
- SEALED RECEIPTS AND AUTHORITY: `CitationReceipt` cannot be constructed
  outside `Registry::cite`, and `cite` refuses every registry except the
  seeded one behind `registry()` (a private authority marker that public
  constructors cannot set). Synthetic rows and caller-built registries
  can be validated and linted but can NEVER reach the receipt/color-cap
  API. There is no caller-asserted evidence flag that upgrades a cap.
- MUTATION-SENSITIVE IDENTITY: every semantic field of an entry moves
  `entry_digest` (length framing, variant tags, exact IEEE-754 float
  bits). Registry input order cannot move `Registry::digest()`: rows sort
  by (id, content identity) and references by their full field tuple, so
  even conflicting duplicate-key rows land in one canonical arrangement.
- ROW/IDENTITY AGREEMENT: `canonical_row` preserves the deck variant and
  state (authored / external / malformed-external / unpinned) and uses one
  canonical hex spelling; a valid external digest is normalized to its raw
  32 bytes in the identity, so hex case cannot fork either surface. Oracle
  state (unpinned / self-contained / derivation-required) is likewise
  distinct in both the row and identity.
- NO AUTHORITY-BY-CITATION: `PrimaryReference` has no color/authority API;
  receipts cap colors, they never mint them; composition cannot upgrade
  the physical-prediction cap without independent held-out evidence.
- DERIVATION-REQUIRED DECKS: the Bennett mobility, Geneva, Sod, Lax,
  isentropic-nozzle, and Atkinson G1 targets pin their parameterization and
  mandatory limit checks while delegating a load-bearing derivation; they
  are not mnemonic-formula oracles.
- G2 seeds stay uncitable until edition/license/deck/oracle/QoIs/envelopes
  are pinned; pinning them is downstream work, not a tolerance relaxation.
- EXECUTABLE ENVELOPES ARE FAILURE GATES, NOT AUTHORITY: every finite
  completed comparison retains the seeded entry/QoI/unit, exact entry and
  registry identities, reference or interval, computed value, derived
  tolerance/deviation, signed margin, and pass state. The attempt/verdict fields
  are private, so callers cannot forge a passing registry record. Negative
  margins return `EnvelopeGateError::Violation`; malformed or unpinned stored
  definitions, mismatched modes, NaN/infinity, and arithmetic overflow never
  produce a passing verdict. Pre-verdict arithmetic failures retain the exact
  definition and observation needed for replay. The gate cannot mint
  `CitationReceipt`, `ColorRank`, or trust.
- EXACT STANDARDS SOURCE BINDING: standard family, part, edition, and ordered
  amendment/corrigendum chain are semantic identity. Unknown editions never
  fall back; current-use rules cannot bind withdrawn or superseded rows; and a
  supersession edge must resolve to an exact in-manifest key without cycles.
- SOURCE AUTHENTICATION AND ACCESS ARE FAIL-CLOSED: unpinned, zero-hash,
  revoked, mismatched, or unread sources cannot mint `RuleProvenance`.
  Historical rows remain inspectable and may bind only under an explicit replay
  policy; they never regain current support through a successor.
- PROTECTED TEXT STAYS OUT OF THE MODEL: source bytes are supplied out of band.
  The manifest retains bibliographic metadata, an external locator, a complete
  source hash, policy identifiers, and derived-rule identity. Redacted rows
  additionally omit title, license terms, revocation reason, and no-claim prose
  while preserving exact keys, locator, source/row hashes, status, and coarse
  access state.
- CANONICAL STANDARDS IDENTITY: source rows sort by exact key before encoding;
  every semantic row field, ordered change, status edge, and policy variant is
  length-framed into the row and manifest identities. Rule provenance binds the
  schema, manifest, source row, authenticated source hash, exact clause,
  derived-rule id/hash, reference state, use policy, and historical bit.
- EARLIEST-RETAINED CORPUS BINDING: an admitted corpus row has a nonzero retained
  payload digest and positive byte length. `OriginalRaw` is required for a
  physical cap above `Estimated`; a `DerivedOnly` row must carry the reason the
  original is unavailable. Complete preprocessing ordinals are contiguous from
  zero, the first input equals the retained digest, every subsequent input
  equals the preceding output, and `final_artifact` equals the last output. An
  unreplayable historical lineage instead binds its retained input/output
  hashes and states why replay is impossible. Retention policy must preserve
  payloads and calibration evidence together.
- DIMENSIONED UNCERTAINTY AND CONTEXT: sensor placement is finite length data,
  bounded uncertainty matches the measured quantity dimensions, covariance
  diagonal entries use the squared dimensions, acquisition uncertainties are
  nonnegative and dimension-compatible, and context/envelope regimes are
  finite, ordered, and nested. `Unstated` is explicit and admitted only as a
  color demotion; absence of the uncertainty field is still a missing-field
  refusal.
- PORTFOLIO AXES DO NOT LAUNDER: `EvidenceLevel` deliberately has no ordering.
  Field monitoring, cross-code agreement, and numerical verification cannot
  substitute for controlled experimental validation. Repeated observations on
  one axis are idempotent, different QoIs or regimes never satisfy one another,
  and a claimed independent reproduction must use a different source and
  declared independence group. Corpus audit rows aggregate each axis
  separately.
- ADVERSARIAL HONESTY IS NOT ACCURACY AUTHORITY: case declarations are
  immutable, bounded, identity-sensitive challenge specifications. Retained
  dataset ids must resolve in the authoritative corpus and their A-E
  coordinate must match the declared evidence basis; planned cross-code or rig
  evidence remains visibly unretained. A prediction outside its envelope
  cannot pass, while refusals and demotions cannot pass by naming an unrelated
  uncertainty. Duplicate/foreign assessment rows refuse scorecard assembly,
  and absent rows remain `NO-DATA`.
- CORPUS AUTHORITY SEPARATION: caller-built registries can validate, serialize,
  hash, and audit rows but cannot return evidence-bearing query results. The
  workspace seed behind `corpus()` is the only query authority. It includes 19
  Level-A thermal reference/target definitions, an explicitly synthetic
  Level-B migration of `fs-benchmark` CHT query `cht-q3`, four external
  cross-code Level-B thermal references (seeded only after their retained
  manifest passes fail-closed verification; a corrupted manifest panics the
  seed rather than shrinking the registry), the retained
  Martin-Moyce 1952 digitized square-column curve, and three published
  electronics-cooling Level-C records: Pires-Fonseca flat/strip fins, Nunes
  HFE-7100 micro-pin fins, and Markal-Kul fin distributions. The new records
  bind retained publisher source bytes, licenses, nominal geometry,
  measurement-uncertainty statements, validation partitions, context ranges,
  QoI declarations, and complete retained-source-to-final lineage. Pires-
  Fonseca and Nunes retain explicit plot-digitization half-widths; Markal-Kul
  retains the publisher supplementary archive. The acquisition log also pins
  five rejected candidates and their byte-access, metrology, geometry, or
  payload failures. Original instrument histories, calibration certificates,
  acquisition dates, as-built metrology, and governed acceptance envelopes are
  not invented: the affected rows remain physically `Estimated`, and every
  evidence query remains numerically `NoClaim`.
- PURPOSE BEFORE DATA: the raw declared-partition query is crate-private. Public
  callers must cross `PartitionLedger` with a semantic purpose; repeating the
  stored partition name cannot bypass purpose checks. Blind data cannot be
  queried through ordinary validation, and every repartition stales all older
  access receipts by generation.
- CALIBRATION TAINT IS TRANSITIVE AND CONSERVATIVE: model registration accepts
  only fresh calibration-purpose accesses. Training data enter the same taint
  class as calibration data. Validation compares dataset content identities,
  so two contexts from the same dataset are treated as overlapping rather than
  laundering reuse through a coordinate change. Parent paths are complete up
  to the explicit depth/item caps; cycles and empty lineage refuse.
- RECEIPT IDENTITY IS REPLAY-COMPLETE FOR THIS LAYER: access identity binds the
  exact query context after sorting axes by name and encoding SI values as
  IEEE-754 bits plus signed dimension exponents. Validation retains every
  distinct access identity, including multiple coordinates from one held-out
  dataset. Caller ordering cannot change access, taint, or validation identity.

## Error model

Total functions; no panics in library paths. Citation failures are
`CitationRefusal` values (recoverable, typed, actionable); executable-envelope
failures are `EnvelopeGateError` values retaining either the complete sealed
failing verdict or the exact registry-bound definition, observation, and
non-finite-input/arithmetic-overflow refusal; consumption binding failures are
`ConsumptionRefusal`; seed-authoring defects surface as `IntegrityFinding`s
from `lint()`, not as crashes.

Standards-manifest admission and decoding return typed `ManifestError` values
for validation, graph, cap, allocation, framing, UTF-8, tag, version, and
canonicality refusals. Rule admission returns `RuleBindingError`, retaining the
exact edition and expected/observed hashes where relevant. No partially
admitted manifest or provenance object is published on error.

Corpus admission and decoding return `CorpusError`; query failures return
`CorpusQueryRefusal` wrapped by `PartitionRefusal`. Purpose, freshness,
repartition, blind-release, lineage, taint-cycle, resource, and taint-
intersection failures are direct `PartitionRefusal` variants and retain the
dataset, generation, attempted purpose, or model path needed to diagnose the
boundary. Both preserve the exact failing field or query boundary.
No partially admitted dataset or successful evidence wrapper is returned on a
refusal.

Portfolio construction and claim admission return `PortfolioRefusal`. Missing
axes name the exact claim class, axis, QoI, and regime; malformed identifiers,
zero content identities, resource overflow, and shared-group reproduction each
remain distinct. No partial `EvidencePortfolio` or `PortfolioAdmission` is
published.

## Determinism class

Fully deterministic: seed metadata and tracked fixture bytes are fixed inputs;
serialization renders
floats as bit tokens (never locale/formatting dependent); digests are
domain-separated BLAKE3 over length-framed canonical bytes. Bitwise
reproducible across runs, thread counts, and ISAs. Envelope-verdict JSON uses a
fixed field order and exact IEEE-754 bit tokens for every floating-point value.
The standards manifest is likewise fully deterministic: input order is erased
by exact-key sorting; its integer/framing encoding is fixed little-endian; row,
manifest, and rule identities are domain-separated; and canonical decode
re-encodes and byte-compares before admission.

Corpus datasets sort sensor, environment, context, acceptance, and regime rows
by stable keys; preprocessing sorts by ordinal and is then chain-validated.
The versioned binary wire length-frames all strings and collections, stores
floating-point values as exact IEEE-754 bits, and is byte-stable across runs,
thread counts, and ISAs. Corpus registries sort by dataset id before hashing and
audit rendering.

Partition state uses sorted dataset maps, monotonically increasing integer
generations, fixed little-endian framing, domain-separated identities, and
canonical sets/maps for query contexts, taint sources, and validation accesses.
Equivalent caller order therefore cannot move a receipt identity. Exact float
bits are identity-bearing; numerically equal but bit-distinct coordinates such
as signed zero remain distinct provenance.

Portfolio observations sort by axis, exact QoI/regime bytes, source, and
independence group before versioned domain-separated hashing. Exact duplicates
are idempotent; input order cannot move portfolio or admission identity.

The scorecard is byte-identical for equal inputs: cells live in a sorted map
keyed by (metric, canonical regime label), run rows sort by (dataset id, run
identity, exact value bits), gap rows sort lexicographically, and both renders
are pure. Markdown/JSON float display uses Rust's shortest-round-trip
formatting (value-injective) and the JSON additionally carries exact IEEE-754
bit tokens; the scorecard identity is a domain-separated BLAKE3 over
length-framed typed fields, not over rendered bytes.

## Cancellation behavior

None; operations are synchronous with no cancellation points. Honest
cost model for caller-supplied data: `Registry::build` sorts
(`O(n log n)` comparisons over rows and references, with one content
digest per row for the canonical tie-break); `lint`/`canonical_rows`/
`digest` are linear in rows plus per-entry gate cost; the
executable-envelope gate performs a linear entry/QoI lookup and one registry
digest before constant-cost scalar arithmetic; the
`MAX_QOIS_PER_ENTRY` cap is checked before any QoI traversal, so QoI-count
work and the quadratic name-comparison count are capped at 64 and 64².
String byte lengths outside registry ids remain uncapped.
Enforced input caps: `MAX_QOIS_PER_ENTRY` on gate checks,
`MAX_LOOKUP_ID_LEN` plus lowercase-ASCII-slug validation on
entry validation, `Registry::cite`, and `Registry::check_acceptance_envelope`,
and `MAX_BEAD_ID_LEN` on
`ConsumptionRecord::bind` (validated before any trim or copy).
Row/reference COUNTS are uncapped — see no-claim boundaries.

Standards-manifest work is synchronous and non-preemptible, but explicitly
bounded before traversal/allocation by `ManifestLimits`: record count,
changes-per-record, bytes per string, and total canonical bytes. Default hard
caps are 4,096 rows, 64 changes per row, 4,096 bytes per string, and 16 MiB per
manifest. Construction sorts rows and validates the functional supersession
graph in `O(n log n)`; exact lookup is `O(log n)`.

Corpus admission/query/audit are synchronous and non-preemptible. Admission is
bounded before allocation/traversal by 4,096 datasets, 4,096 sensors per
dataset, 4,096 entries in every other collection, 4,096 UTF-8 bytes per string,
and 16 MiB canonical bytes per dataset. Registry construction is
`O(n log n)`; dataset admission sorts bounded collections; query is one binary
dataset lookup plus linear bounded context matching; audit is linear in
datasets plus admission cost.

Partition access, repartition, release, registration, and validation are also
synchronous and non-preemptible. Each bounded operation admits at most 4,096
direct items; taint explanations admit depth at most 256; justifications admit
at most 4,096 UTF-8 bytes. Transitive closure and validation use ordered maps or
sets and are `O(n log n)` within those caps.

Portfolio construction is synchronous and non-preemptible, capped at 4,096
observations and 256 bytes per QoI/regime identifier. Canonicalization is
`O(n log n)`. Ordinary claim admission is one bounded scan per required axis;
independent-reproduction admission additionally performs a deterministic
worst-case `O(n^2)` pair search within the same hard cap.

## Unsafe boundary

None. Workspace `deny(unsafe_code)` applies.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs`: G1 seeds split into 12 citable receipts + 6
derivation-blocked refusals; the synthetic-forge regression (a fully
pinned synthetic row validates but no caller-built registry can mint a
receipt, and seeded receipts bind the registry digest); the bead's named
fixtures — registry lint
partition and the unpinned family-name citation refusal (TEAM 10) — plus
unknown-id refusal, duplicate-id fail-closed citation and lint exclusion,
bounded/malformed lookup-id refusal before input copies,
duplicate reference keys at non-adjacent indices, gate ordering probes
(including oracle-before-QoI and dedup-before-envelope), invalid-envelope
reasons, executable tolerance and interval gates with inclusive endpoints,
disclosed seeded boundary-plus-ULP corruption, exact-bit JSON verdicts, and
fail-closed mode-mismatch/non-finite/overflow/unpinned probes, seeded
entry/QoI/digest binding, unknown-entry/QoI refusal, and caller-built-registry
forgery refusal, the canonical-row golden for the unpinned TEAM 10 row,
sorted/deterministic serialization with bit-token floats and variant-
tagged 64-hex deck digests, deck row/identity agreement (case
normalization, malformed-vs-unpinned distinction, authored-vs-external
separation), order-canonicalized registry digest, per-field mutation
sensitivity of `entry_digest` (oracle included), Appendix-D
consumption-record round trips, color-cap laws against the `ColorRank`
lattice, the complete 1..=30 primary-reference seed pinned by ordered
(key, locator) identity, and per-field reference digest mutation locks
(index, key, citation, locator, anchors, boundary each move the registry
digest).

`tests/standards_manifest.rs`: G0 exact-edition, source-pin, hash-match,
access-revocation, explicit-reference-state, zero-hash, collision, change-chain,
resource-cap, and closed-acyclic-supersession refusals; explicit historical
replay admission; G3 input-order invariance, ordered-change sensitivity, and a
protected-text leakage/redacted-JSON fixture; G5 every-field source identity
mutation, rule-provenance mutation, canonical v1 round trip, frozen `FSMF`
header, future-schema refusal, truncation/trailing-byte refusal, and
semantic-but-unsorted wire rejection.

`tests/corpus.rs`: G0 typed refusal for every mandatory top-level field;
calibration, placement, bounded-uncertainty, covariance-dimension, retention,
lineage, canonical-codec, tamper, and caller-authority probes; G3 input-order
invariance for registry identity and explicit authority-gap demotion; exact
partition plus missing/unknown/duplicate/dimension/out-of-range context
refusals; retained synthetic, Martin-Moyce, and three electronics-thermal
source/final artifact hash bindings; independent Pires-Fonseca and Nunes
re-digitized subsamples checked against declared half-widths; acquisition-log
admission/rejection counts; and deterministic audit rendering with warn-level
optional/claim-authority gaps plus the complete seven-axis QoI coverage map.

`tests/portfolio.rs`: G0 exact A-E-to-axis mappings, typed missing-axis
refusals for every claim class, field-only validation refusal, exact QoI/regime
matching, duplicate idempotence, and distinct-group independent reproduction;
G3 order-invariant portfolio/admission identities and source/regime mutation
sensitivity.

`tests/partition.rs`: G0 purpose/partition matrix, direct and transitive data-
laundering refusals with exact model paths, stale access refusal during model
registration and validation, blind release gating, and wrong-purpose model
input refusal; G3 order-independent repartition/taint identities, query-context
identity sensitivity, retention of multiple held-out accesses from one dataset,
and disjoint held-out success constrained to a non-certifying taint-check
receipt; canonical round trips for all five receipt variants, future/truncated/
extended/tampered wire refusal, and exact-byte fs-ledger persistence across
reopen with dedupe plus corruption detection.

`tests/thermal_level_a.rs`: G0 manifest/catalog identity and family coverage;
independent recomputation of all 12 closed-form scalar values; exact retained-
artifact binding; seeded query/no-claim checks over every context axis; and G1
target coverage for two element degrees, Neumann/Robin boundaries, nonlinear
anisotropy, and primal/adjoint order. These tests verify reference definitions
and targets, not thermal-kernel convergence.

`tests/thermal_level_b.rs`: G0 catalog well-formedness (unique ids, positive
extents/counts, symmetric-PD tensor via leading minors, increasing positive
`k(T)` knots, picard-iff-nonlinear, in-grid probes); fail-closed verification
of the committed manifest; BLAKE3 deck-byte binding; typed parser refusals for
unknown kinds, wrong column counts, non-finite/unparsable numbers, duplicate
keys, non-dense probe indices, missing mandatory metadata, and a recorded
external self-check other than `pass`; ULP-scale spec-echo and probe-position
tampers refused bit-exactly; nonlinear-case iteration-count sanity; and corpus
registration with `CrossCode` level, `Estimated` physical cap, retained-locator
identity, and pinned per-case envelopes. The executing cross-code comparison
itself lives in `fs-conduction/tests/level_b_crosscode.rs`.

`tests/scorecard.rs`: G0 known-error cell arithmetic on a synthetic corpus
(exact regime label, axis/external counts, pass/fail/unpinned split, exact
signed errors); loud NO-DATA for empty cells and zero executed challenges
(never rendered as zero); reference-uncertainty rendering in both renders;
G5 byte-identity of Markdown, JSON, and identity across independent builds;
the real-corpus e2e lane logging per-cell inputs and the gap list while
pinning the seeded honesty header; a seeded false acceptance surfacing in
exactly the cells of its retained challenge dataset and nowhere else; and
fail-closed refusals for non-finite values, negative half-widths, blank ids,
unknown datasets/metrics, the run-record resource limit, and foreign or
duplicate assessments.

## No-claim boundaries

- Admission proves the ENTRY is fully pinned; it does not prove any solver
  result, and it does not verify that an external deck's bytes exist or
  match their registered digest — artifact retrieval/verification is the
  consuming lane's job.
- The scorecard is a diagnostic projection, never authority: supplying run
  records or assessments cannot upgrade any corpus claim cap, an in-envelope
  run row is envelope arithmetic (not an evidence colour), and cell
  false-acceptance attribution follows the retained challenge dataset only —
  planned-evidence cases stay in the adversarial section as `NO-DATA`. The
  scorecard does not verify that a `ScorecardRunRecord`'s ledger identity
  exists in any ledger or that its values were actually produced by the
  named run; binding run provenance is the persisting lane's job (e05
  scope). No retained ledgered run-result store exists yet, so the
  committed repository artifact renders every prediction-error, coverage,
  and false-acceptance cell as `NO-DATA`; the embedded adversarial
  section's `false_acceptance_count` line counts supplied assessments and
  is subordinate to the header's `NO-DATA` marker when nothing executed.
- An executable envelope proves only that a caller-supplied scalar satisfies
  the exact arithmetic rule stored on the named seeded-registry QoI. It does not
  bind a reference to an oracle, bind a computed value to a deck/run, establish
  that the whole entry is citable, or grant evidence authority; the consuming
  lane must bind those identities and provenance.
- The 15 G2 seeds are targets, not benchmarks: no claim may cite them
  until their decks are pinned (exact edition, license, deck hash, oracle,
  QoIs, acceptance data).
- G1 acceptance envelopes bound agreement with the authored analytic
  oracle under its stated assumptions; they say nothing about physical
  validity (the physical cap stays `Estimated` without held-out evidence).
- The whole-registry BLAKE3 digest golden constant is deliberately NOT
  frozen in this crate yet: per `docs/GOLDEN_POLICY.md` a golden pin
  requires committed-tree, two-mode (and where claimed, two-ISA)
  reproduction, which is scheduled with the batch-verify lane. The
  serialization goldens pin exact row strings instead.
- `ConsumptionRecord` records discipline; it does not enforce that the
  consuming bead actually read/derived/reproduced anything.
- The `Validated` physical-cap upgrade is deliberately absent: it requires
  a typed, non-forgeable binding of independent held-out evidence (future
  work tracked on the bead), not a boolean argument.
- Derivation-required entries (Bennett mobility, Geneva, Sod, Lax,
  isentropic nozzle, and Atkinson) stay non-citable until a derivation
  receipt mechanism exists; their registration is a target declaration,
  not an oracle.
- No registry-size caps: rows are compiled-in seed data here, and a
  caller-built `Registry` is the caller's resource decision —
  `build`/`canonical_rows`/`digest` do not police hostile row counts.
- A blank authored spec has no well-formed deck digest (`DeckPin::digest`
  is `None`) and renders as `{"authored":null}` — a distinct visible
  state; the entry identity still covers the raw state bytes.
- A standards manifest proves only that metadata is structurally valid and
  content-addressed. It does not prove that a publisher locator resolves, that
  a hash is publisher-authoritative, that access is legally authorized, that a
  declared edition is actually current in a jurisdiction, or that a derived
  rule is semantically correct.
- `ProtectedTextReference` excludes a body from the API and redacted rows omit
  sensitive metadata, but arbitrary public metadata strings cannot be
  cryptographically classified as copyright-safe. Callers remain responsible
  for supplying only identifiers/locators/policy codes and for keeping licensed
  source bytes out of manifest fields.
- `ReferenceState::{Read, Derived, Reproduced}` is explicit provenance state,
  not independent proof that the human or agent performed the declared work.
  It affects identity and makes unread use impossible; downstream evidence must
  substantiate stronger claims.
- Historical replay is provenance preservation, not standards support: a
  withdrawn or superseded edition remains historical even when its exact bytes
  are still pinned and accessible.
- Corpus admission proves schema completeness, internal hash references,
  dimension compatibility, lineage closure, and policy declarations. It does
  not retrieve external artifacts, authenticate a lab, verify calibration
  certificate signatures, establish population representativeness, prove a
  solver result, or establish that an acceptance envelope is scientifically
  appropriate.
- `EvidenceLevel` is a compact legacy source tag, not a rank.
  `EvidenceLevel::PublishedExperiment` contributes the controlled-experiment
  coordinate but does not mean original instrument data survived. Paper figures and publisher
  tables reduced from unretained sensor histories are `DerivedOnly` even when
  the article calls them data. Digitization half-widths bound coordinate
  extraction only and do not replace the paper's experimental uncertainty or
  missing calibration authority. Level E contributes field monitoring only and
  cannot mint a `Validated` cap without separate controlled-experiment
  evidence. The per-axis QoI map counts declarations, not validated solver
  comparisons, population coverage, or acceptance authority.
- A corpus query proves only that the seeded row was revalidated and that the
  caller requested an admitted purpose inside its declared context under one
  captured in-memory partition generation. The
  returned `Evidence<&CorpusDataset>` deliberately has numerical `NoClaim` and
  cannot be certified. Its statistical and model-form components are likewise
  unknown/unbounded rather than silently absent, while model validity retains
  the dataset context box for conservative downstream intersection. The
  dataset-level color is a maximum support cap, not a minted `Color` or a
  validated prediction.
- Partition/repartition/blind-release/access/model-taint/validation receipts
  are exact canonical records. Their v1 wire records are self-consistent, and
  the integration battery proves that exact bytes survive fs-ledger's generic
  content-addressed artifact write, reopen, dedupe, bounded read, and
  corruption scan. The UTIL runtime still retains partition state and ordered
  events only in memory. A generic artifact row is not dedicated receipt
  membership, a query index, an atomically coupled audit event, a refusal log,
  HELM orchestration, authorization, secrecy of blind bytes, crash-fault proof,
  or partition-state restoration.
  The nonzero preregistration and manifest hashes are caller-supplied bindings,
  not proof that either artifact exists, predates model freeze, was independently
  witnessed, or was honestly concealed. Full blind-protocol orchestration and
  adversarial end-to-end audit logging remain downstream work.
- A clean taint intersection proves only that the exact registered dataset
  identities are disjoint. It cannot detect unregistered copies, transformed
  leaks, common upstream source material, human prior exposure, fabricated
  corpus metadata, or semantic dependence. It never upgrades numerical,
  statistical, physical, or model evidence color.
- `EvidencePortfolio` validates bounded structural coordinates and required-axis
  presence only. Source and independence-group hashes are caller-supplied
  identities, not signatures or authenticated independence. A
  `PortfolioAdmission` does not mint `Color`, `AdmittedColor`, a corpus access
  receipt, or runtime authority. Durable authentication and exact-instance
  policy remain `fs-package`/`fs-checker`/`fs-ledger` scope.
- `AdversarialAssessment` is a sealed diagnostic honesty receipt, not a solver
  result, evidence color, validation claim, or ledger authority. Registration
  does not execute a challenge. The retained analytic cases are reference
  definitions, the retained fin case has the corpus row's existing metrology
  gaps, and planned fan, vent, material-lot, and cross-code cases have no
  acquired result. The deterministic regime-limitation Markdown is an input
  surface for the downstream scorecard; it is not the DSR artifact lane or the
  public scorecard itself.
- The seed CHT row is synthetic Level B and therefore physically `Estimated`.
  Its CSV is an authored tabulation of a hard-coded query, not raw sensor data;
  the stored `1.0 K` value is an acceptance tolerance, not measurement
  uncertainty. Instrument, calibration, placement, geometry, environment,
  acquisition-date, and replayable-export authority are explicitly unavailable
  and emitted as WARN gaps. It is a worked schema migration, not external
  validation.
- The Martin-Moyce row binds the exact retained digitized JSONL bytes and the
  live consumer path, not the lost experiment or digitization workflow. The
  original raw records, instruments, calibration, placement, acquisition
  environment/date, digitizer lineage, redistribution authority, and a
  defensible scalar acceptance envelope are unavailable or unresolved. The
  row is Level C by source family but remains physically `Estimated`, carries a
  numerical `NoClaim`, and cannot establish L4 experimental validation.
- The 19 Level-A thermal rows are reference definitions and theoretical order
  targets. Their tests independently recompute closed-form scalars and enforce
  the target matrix, but no row binds a thermal solver output, mesh/refinement
  ladder, adjoint run, or machine fingerprint. P2 primal/adjoint rows and all
  unimplemented domain families are intentionally target-only. Registration is
  therefore not a G1 pass, not solution verification, and not evidence that a
  thermal kernel exists; a consuming crate must retain its own comparison or
  `fs-mms::OrderVerdict` receipt before making that claim.
- The four cross-code Level-B rows are SAME-DISCRETIZATION parity references:
  by construction the external code assembles the identical discrete system
  (bit-identical Kuhn mesh, P1, element-mean `k(T)`, consistent source/Robin
  mass), so within-envelope agreement checks independent assembly/boundary/
  solver implementations and says NOTHING about discretization error, mesh
  convergence, or physical validity. Two codes agreeing is not truth: every
  row stays `Estimated` with a numerical no-claim, and the external code
  (scikit-fem/scipy) is a development-only oracle that never enters the
  runtime dependency graph. This corpus does not yet contain any
  CalculiX/Elmer-class independent-discretization reference; a case solved on
  a DIFFERENT mesh or element family would need new case ids, new envelopes,
  and its own corpus version bump.
