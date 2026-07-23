# CONTRACT: fs-project

The versioned `.fsim` project schema (bead f85xj.6.1): the user-facing
contract for the ratified thermal-design-assurance vertical
(`frankensim-vertical-ratification-v1`). One semantic model, two spellings,
canonical bytes, typed violations, receipted defaults, and explicit migration
receipts. The `bind` module (bead f85xj.6.4) adds validate-time resolution of
material/interface bindings against matdb cards: envelope-vs-domain refusal,
surfaced uncertainty, retained usage receipts, and pin-only resolution of
conflicting claims. The `assignment` module (bead f85xj.6.3) compiles
persisted mesh-index-free selectors to persistent `EntityId` tokens, resolves
them against caller-supplied promoted meshes through `fs-io`, and retains the
exact lower-layer reports.

## Purpose and layer

Layer L6 (HELM). Depends on `fs-ir` (AST + both concrete syntaxes),
`fs-scenario` (entity identities, `Violation`), `fs-matdb` (cards, receipted
queries), `fs-io` (deterministic mesh assignments), `fs-rep-mesh` (promoted
finite meshes), `fs-exec` (cancellation context), `fs-qty`, and `fs-blake3`.
This crate persists project intent and resolves its card and geometry
bindings; it runs no solves and admits no scenarios itself.

## Public types and semantics

- `ProjectSpec` carries every section of the cooling vertical's contract:
  `Metadata` (name, created, context of use, intended decision), the Five
  Explicits (`Versions`, `Seeds`, `Budgets`, capabilities list,
  `UnitsDoctrine`), `GeometryArtifact` references (imported quarantine-receipt
  ids — geometry never lives inline), `EntityDecl`
  assembly/part/region/interface declarations with persistent identities,
  mandatory `GeometryAssignment` rows (one declared artifact, persistent
  region/interface target, explicit coordinate unit, one `MeshSelector`, and
  overlap policy),
  `MaterialBinding` matdb card refs (state, admitted temperature range,
  source, optional `claim` pin), `InterfaceCardBinding` TIM/contact card refs
  (optional `claim` pin), `PowerDissipation`
  rows, `Cooling` (fans/vents/leakage; declared-empty lists are facts,
  omission is a violation), `Envelope`, `ThermalLimit` requirements with
  margins, `SolverSettings`, and `OutputRequest`s. Sections are `Option`s so
  recognition stays lenient and validation can name every omission at once.
- `ProjectSpec::validate() -> Vec<Violation>` (the fs-scenario
  code/what/fix triple): empty output is the definition of admissible. Every
  mandatory section omission has a stable `project-*-missing` code; empty
  load-bearing collections, dimension mismatches (checked against the six-base
  SI vectors in `spec::dims`), inverted ranges, malformed card hashes, duty
  range, unknown entity references, duplicate names, and out-of-order parents
  all carry named codes with fix hints.
- Entity identities: `resolve_entities` recomputes `fs_scenario::EntityId`s
  from the declarations (parents before children; interfaces reference two
  regions under an assembly). An optional `:id` pin is verified against the
  recomputation and drift is `project-entity-id-mismatch` — a pin proves
  byte-equal derivation inputs, never physical sameness. Display names are
  outside identity, exactly as in fs-scenario.
- Wire: `lower`/`recognize` map `ProjectSpec` to and from the `fs_ir::Node`
  envelope `(fsim-project :version 1 ...)`. `print_sexpr`/`parse_sexpr` and
  `print_json`/`parse_json` are the two spellings; `parse_sexpr_lenient`
  accepts noncanonical bytes and omitted defaultable fields, issuing a
  `CanonicalizationReceipt` (both hashes, `verifies()`) and `DefaultReceipt`s.
  The strict parsers refuse noncanonical input (`fsim-non-canonical`, first
  difference position) and applied defaults (`fsim-default-in-strict-mode`).
  The canonical `assignments` section spells named-group, half-space, box,
  cylinder, nearest-datum, and explicit-face-set selectors without depending
  on imported face ordinals except where the explicitly fragile
  `explicit-face-set` variant is intentionally chosen and acknowledged.
- Canonical bytes are the checked s-expression render; `canonical_hash`
  hashes them under `org.frankensim.fs-project.canonical.v1`. The JSON
  spelling parses to the same AST, so both spellings reach one hash.
- `FSIM_VERSION = 1`. Readers refuse other versions
  (`fsim-unsupported-version`); `migrate_envelope` is the only path from an
  older envelope, applies a registered `MigrationRule`, and returns a
  `ProjectMigrationReceipt` (old/new hashes + rule, `verifies()`). The only
  registered rule is the synthetic version-0 proof rule, named as such.
- The ONE defaultable field is power-row `duty` (`1.0`, continuous
  dissipation is the conservative assumption); everything else is mandatory,
  and the default is always receipted.
- `bind::resolve_bindings(&ProjectSpec, &CardLibrary, &BindingRequirements)
  -> MaterialResolution` (bead f85xj.6.4): every declared region must bind
  exactly one `fs_matdb::MaterialCard` (manufactured state spelled exactly as
  the card's `MaterialStateId`) and every declared interface exactly one
  `InterfaceSystemCard`. Per required property, the resolver queries matdb at
  BOTH endpoints of the binding's admitted range (validity boxes are per-axis
  intervals, so endpoint containment implies range containment) and requires
  the SAME claim selected at both. Every value arrives as matdb
  `Evidence` + `PropertyUsageReceipt`; receipts are retained
  (`RetainedReceipt`: receipt, canonical bytes, content hash, context) for
  ledger retention and later `verify_receipt` replay. `render_table()` is the
  deterministic region -> card -> property -> uncertainty -> receipts log.
  `BindingRequirements::thermal_steady_v1()` pins the cooling vertical's
  property set (`thermal-conductivity`,
  `area-specific-thermal-contact-resistance`, axis `T`), drift-tested against
  `fs-conduction`'s constants.
- Conflicting-claims doctrine: when a card carries coexisting conflicting
  claims, resolution refuses (`project-binding-claims-conflict`, listing every
  candidate hash) unless the project file records an explicit `:claim <hex>`
  pin, which resolves through matdb's receipted pinned-claim query. There is
  no auto-pick path, and `PreferObservationBacked` is deliberately never used
  by the resolver. A pin never bypasses validity.
- The per-region required range is the envelope's ambient range with its
  ceiling raised to every declared thermal limit on that region (data must be
  valid up to the limit being judged); an interface must be covered wherever
  either side region can operate. The binding's admitted range must cover the
  required range (`project-material-envelope-uncovered`), and the card's
  claims must cover the admitted range (`project-binding-domain-uncovered` —
  the extrapolation refusal, surfaced at validate time).
- Uncertainty surfacing: `Unstated` uncertainty is an `Advisory`
  (`binding-uncertainty-unstated`) up front, never a refusal and never
  laundered — it caps downstream evidence at Estimated and the table says so.
- `assignment::resolve_geometry_assignments(&ProjectSpec,
  &ImportedMeshLibrary, AssignmentLimits, &Cx) -> GeometryResolution`
  resolves declarations first, compiles project-local targets to the actual
  `EntityId::token()` strings, then delegates every artifact's selector plan
  to `fs_io::resolve_mesh_assignments`. `ImportedMeshLibrary::insert`
  computes its own domain-separated key from the exact `GeometryArtifact`
  row; callers cannot supply a mismatched map key. Every successful
  `ResolvedGeometryArtifact` retains fs-io's report unchanged, its exact
  canonical JSON bytes, a domain-separated report hash, and the declared-name
  to `EntityId` correspondence. `render_table()` surfaces entity, artifact,
  source identity, unit, selector fingerprint, selected face count, area,
  optional enclosed volume, bounds, and report hash.

## Invariants

- Parse-render-parse is byte-idempotent for canonical documents in both
  spellings, and `print_sexpr` is deterministic.
- The strict parsers accept exactly the canonical bytes; the lenient parser
  never accepts silently (receipt or violation for every relaxation).
- Unknown sections and unknown keywords are refused by name
  (`project-unknown-field`), so typos cannot silently drop intent.
- Recognition never panics on malformed clause shapes; it degrades to named
  violations with placeholder values that validation then rejects.
- `Violation.fix` is non-empty for every emitted code (tested via the broken
  corpus).
- Every declared region/interface has exactly one geometry assignment; every
  declared geometry artifact has at least one assignment; roles and targets
  are unique where required. Geometry reports publish atomically only after
  every artifact succeeds, so refusal or cancellation exposes no partial
  assignment result.
- The adapter preserves fs-io assignment order and rejects any returned
  subject/order mismatch before retention. Project-wide request and selected
  face counts are checked against the explicit `AssignmentLimits`.

## Error model

Two tiers, exactly as in fs-scenario/fs-ir: wire-layer refusals are
`ProjectError { code, detail, hint }` (syntax, non-project shape, unsupported
version, non-canonical bytes, strict-mode defaults, quantity spelling,
migration refusals); semantic findings are `Violation` triples collected by
`validate()`/`DecodedProject::findings()`. Binding resolution reports the
same triples (`project-binding-*` / `project-material-*` /
`project-interface-*` codes, each with an actionable fix; matdb's own
refusal text travels in `what`) plus non-refusing `Advisory` rows; it is
total over garbage input, naming preconditions instead of panicking.
Internal `expect` is limited to infallible writes to `String`.
Geometry-assignment preflight and adapter failures use the same
`Violation` triple (`project-assignment-*`); fs-io refusal codes, `what`, and
`fix` propagate without changing their lower-layer meaning, with the geometry
role added as context. Refusal and cancellation leave `artifacts` empty.

## Determinism class

Fully deterministic: pure functions of the input bytes/spec and explicitly
supplied card/mesh libraries; canonical rendering has one spelling; hashing
is domain-separated BLAKE3 over exact bytes. Mesh resolution inherits
fs-io's deterministic face-order and selector-fingerprint contract. No
clocks, no RNG, no environment reads.

## Cancellation behavior

Static document recognition and validation do not take a cancellation
context; their descriptions are bounded by the caller's document admission
layer. Geometry assignment resolution takes an explicit `Cx`, polls before
entry, at artifact boundaries, and inside fs-io's bounded face/selector tile
loops. Cancellation is a typed `mesh-assignment-cancelled` violation and
publishes no partial report. A validation budget/plan triad (fs-scenario
style) remains deferred until real projects show the structural loader can be
large.

## Unsafe boundary

None. Workspace `deny(unsafe_code)` lint applies.

## Feature flags

None.

## Conformance tests

`tests/bind.rs` (f85xj.6.4): the reference project binds fully — three
resolved bindings, six retained receipts, each replayed via
`verify_receipt` and byte-decoded via `from_bytes_verified`, with the full
chain logged in the table; envelope-outside-card-domain refuses at validate;
an admitted range narrower than the envelope refuses with both ranges named;
Unstated uncertainty warns without refusing; conflicting claims refuse
without a pin (no resolved row exists — auto-pick impossible) and resolve to
exactly the pinned claim with the `pinned-claim` receipt policy; pin
refusals (unknown, out-of-domain, malformed) are typed; a range stitched
from two claims refuses; card-unknown, state-mismatch, non-region target,
duplicate, unbound-region/interface, wrong-dims, and missing-property each
carry their named code with a fix; missing sections name the precondition;
and the property constants drift-test against `fs-conduction`.

`tests/project.rs`: the reference cooling project renders, parses strictly,
is admissible, and hash-stable; canonical bytes are identical across both
spellings; parse-render-parse idempotence across variants (empty cooling
lists, identity pins); all seventeen mandatory-section omissions surface their
named violations; unknown section and keyword refusal; the receipted duty
default with canonicalization receipt (and strict-mode refusal of the same
bytes); noncanonical whitespace refusal/receipt; the synthetic v0 migration
round trip with verifying receipt plus refusal of no-op and unknown-version
migrations; entity identity pins accepting the recomputed token and refusing
a stale one; and the broken-project corpus (14 rows) logging every violation
with its fix as the error-message quality bar.

`tests/assignment.rs` (f85xj.6.3): a promoted cube resolves named groups to
the exact persistent region/interface identities; retained fs-io JSON and its
domain hash replay exactly; re-tessellation changes source/report material
while preserving the project entity and physical area; dangling and
wrong-kind targets, unit mismatch, empty selection, and unacknowledged
explicit-face fragility refuse without partial publication; pre-cancelled
resolution is atomic. `tests/project.rs` additionally round-trips all six
selector variants in both canonical spellings, including the maximum `u32`
explicit face index.

## No-claim boundaries

- An admissible `.fsim` document proves the project is well-formed and
  internally consistent. It does NOT prove the referenced geometry artifacts,
  matdb cards, or capabilities exist, are compatible, or are validated —
  geometry existence and selector checks require
  `assignment::resolve_geometry_assignments`, not structural schema
  validation. That resolver checks only the caller-supplied promoted
  `ImportedMeshLibrary`: its computed keys bind an entry to one exact project
  geometry row, but WHO supplied the mesh is the caller's trust channel.
  Card existence and coverage ARE checked by `bind::resolve_bindings`, but
  only against the
  caller-supplied `CardLibrary`: the library keys cards by their own content
  hashes (a key cannot lie about its card), while WHO supplied the collection
  is the caller's trust channel — resolution proves binding coverage, never
  card authenticity or scientific truth of the underlying claims.
- Binding resolution is temperature-axis-only in v1: cards whose claims
  depend on further validity axes refuse (`project-binding-axis`) rather
  than being partially resolved. The required ceiling uses DECLARED thermal
  limits, not solved temperatures — a solve exceeding its limit fails its
  requirement; the resolver only guarantees the material data covers the
  range over which that judgment is made. The declared `source` channel on a
  binding is recorded into the table, not authenticated against claim
  provenance.
- A pinned-claim receipt proves the pin resolves and evaluates identically
  on replay; WHO pinned it is proven by the project file's canonical bytes
  (the pin is part of the hashed document), so a verifier must cross-check
  `receipt.selected` against the project's pin, exactly as matdb's CONTRACT
  states.
- The optional `claim` pin was added to version 1 before the schema freeze
  (.16.5): version 1 is still pre-freeze (see below), absence of the field
  is the canonical spelling of "no pin", and pre-pin documents are
  byte-identical under the extended grammar.
- The mandatory `assignments` section was also added to version 1 before the
  schema freeze. Unlike the optional claim pin, this intentionally makes old
  pre-freeze documents incomplete: they must declare one selector per
  region/interface before becoming admissible again. No compatibility shim
  invents geometric intent.
- Entity identity pins prove byte-equal derivation inputs at declaration
  time; unequal pins prove nothing about physical sameness (re-exports change
  bytes), which is why drift is a violation to adjudicate, not an automatic
  rebind.
- The schema does not persist fs-scenario's full `EntityCatalog`/
  `BindingTable` receipts (rebind chains, correspondences, evidence tiers);
  the fs-scenario CONTRACT's delegation of catalog+binding persistence to
  this crate is met for declarations and will extend to receipt chains when
  the import path (.6.3) produces them.
- `source_hash` on geometry artifacts reproduces `fs_io::ImportReceipt`'s
  current 64-bit FNV field: deterministic correlation, not
  collision-resistant authentication; the HELM-side BLAKE3-class upgrade is
  tracked upstream. The adapter's BLAKE3 source identity binds that field and
  the complete geometry row exactly, but cannot strengthen the importer's
  underlying FNV collision guarantee or authenticate the supplied mesh.
- The schema freeze (e16, bead .16.5) has NOT happened: version 1 is the
  first implemented version, not yet a frozen public promise. The migration
  machinery exists so the freeze can be honest when it lands.
- Validation is structural and dimensional; it makes no physics claim (a
  well-formed project can still describe an unsolvable or nonsensical study).
