# CONTRACT: fs-vvreg

The Gauntlet G1/G2 benchmark & V&V registry (bead
frankensim-ext-benchmark-vv-registry-f1gq, epic E0c): the single place
where a benchmark family name becomes an executable claim target — or is
refused until it is one.

## Purpose and layer

Layer UTIL (versioned registry data + fail-closed citation and standards-source
gates). Depends
only on `fs-blake3` (domain-separated content identity) and `fs-evidence`
(`ColorRank` for the citation color caps). A family name (TEAM, NAFEMS,
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

## Determinism class

Fully deterministic: all seed data is `const`; serialization renders
floats as bit tokens (never locale/formatting dependent); digests are
domain-separated BLAKE3 over length-framed canonical bytes. Bitwise
reproducible across runs, thread counts, and ISAs. Envelope-verdict JSON uses a
fixed field order and exact IEEE-754 bit tokens for every floating-point value.
The standards manifest is likewise fully deterministic: input order is erased
by exact-key sorting; its integer/framing encoding is fixed little-endian; row,
manifest, and rule identities are domain-separated; and canonical decode
re-encodes and byte-compares before admission.

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

## No-claim boundaries

- Admission proves the ENTRY is fully pinned; it does not prove any solver
  result, and it does not verify that an external deck's bytes exist or
  match their registered digest — artifact retrieval/verification is the
  consuming lane's job.
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
