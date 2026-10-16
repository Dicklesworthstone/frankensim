# fs-matdb-store — Contract

## Purpose and layer

Layer L6 (HELM). FrankenSQLite-backed queryable store over compiled
fs-matdb material packs (beads frankensim-oecdy and
frankensim-material-reality-ui2if.2.4): SQL answers the
discovery questions (which materials carry a property, which
properties a material carries, what is valid at an ambient point);
every actual evaluation goes through the canonical hash-verified pack
bytes and the in-memory fs-matdb evaluator. The pack corpus remains
the only source of truth; the store is a derived, regenerable index
plus a canonical-bytes vault.

## Public types and semantics

- `MaterialStore::open(path)` — file or `":memory:"`; DDL v2 applied
  transactionally (`STORE_SCHEMA_VERSION` in `PRAGMA user_version`).
  V1 stores gain the explicit `properties` family on existing pack rows;
  canonical bytes and existing claim/validity rows are preserved. The v2
  corpus digest binds family, so a migrated v1 seal requires deliberate
  resealing. Newer/negative schema versions refuse before any DDL.
- `CatalogPack` / `PackKind` — the five existing canonical families:
  property claims, named material cards, ordered interfaces, constitutive
  models, and species associations. `from_bytes_verified` dispatches to
  the exact L1 decoder, including its hash, version, and canonicality checks.
  Material-card v2 carries an explicitly associated, verified model pack;
  ordinary ingest/reopen/content lookup preserves its exact member cards and
  normalization records without changing the SQL schema. Its parameters remain
  models, not indexed scalar properties. A same-state v1 card gains no implicit
  membership from a separately stored model pack.
- `ingest_bundle(&[CatalogPack])` — store related artifacts in one transaction.
  A failure in any member rolls back every earlier member and its index rows;
  commit failures also attempt rollback and report rollback failure explicitly.
  Empty bundles are no-ops. Pack ids remain globally unique across families.
  Bundle bytes are prepared in caller order before the transaction; callers
  bound bundle size to their memory budget. This synchronous API is not a
  streaming ingest or a new concurrent-writer contract.
- `ingest_pack(&NormalizedPack)` — canonical bytes + derived index
  rows (claims, validity axes) in the ClaimSet's canonical order;
  refuses duplicates and empty license/redistribution (the license
  gate survives the store; the license check is a PRE-PASS and the row
  writes are one `BEGIN IMMEDIATE` transaction, so a mid-ingest
  failure leaves no residue and a retry never hits `DuplicatePack` on
  a half-ingested id).
- `seal_corpus` / `require_sealed` — a domain-separated BLAKE3 digest
  folded over pack content hashes in pack-id order; EVERY discovery
  and evaluation surface recomputes and compares, refusing
  (`FS-MATDB-STORE-CORPUS-CHANGED`) on drift and
  (`FS-MATDB-STORE-NOT-SEALED`) before the first seal.
- Discovery: `packs(optional_kind)` returns family/name/content identities;
  `properties_of(pack_id)`, `materials_with(property,
  value_range)`, `valid_at(property, axis, value)` (missing axis =
  unconstrained, matching `ValidityDomain` semantics). Every `PropertyRow`
  includes `pack_kind`, so an interface claim is not mislabeled as a bulk
  property. Equal-property rows use claim hash as a deterministic tie-break.
  Value-range discovery inspects canonical claims and uses the pinned evaluator
  at admitted scalar/sample points and clipped linear-segment endpoints. It
  includes curves with overlapping supported values, without filling exact-only
  sample gaps or extending beyond validity/knot support. Range endpoints must
  be finite and ordered (`InvalidValueRange`). This asks whether some supported
  state matches; it does not admit a complete state envelope or select a claim.
- `evaluate(pack_id, property, &QueryPoint, policy)` — decodes the
  stored bytes via the family's hash-verified decoder and delegates
  to `ClaimSet::query` + `verify_receipt`: the SAME evaluator,
  receipts, and refusal set as direct pack use, passed through
  unchanged (`StoreError::MatDb`).
  Named material and ordered-interface packs expose their original nested
  claims; model parameters and species metadata never become synthetic scalar
  claims (`NoPropertyClaims` on property evaluation).
- `evaluate_typed(pack_id, &PropertyKey, &QueryPoint, policy)` delegates to
  the same canonical decoder and `ClaimSet::query_typed`, then verifies the
  receipt. It retains complete quantity, hardness/tensor context and typed
  axis requirements; no property name or equal dimensions can erase them.
- `discover(&DiscoveryRequest, &LawRegistry)` evaluates a complete typed
  property/model bundle
  separately on each canonical candidate. Targets distinguish named materials,
  explicitly unbound property packs, and interfaces filtered by ordered A/B
  material-state identities. An interface hit still requires exact texture,
  medium, environment and history matching when bound to a simulation.
  Results include every requested property and its evaluated evidence or named
  gap; all globally unknown property names are retained together. Known names
  with no matching target remain distinguishable from unknown names.
  `Complete`/`Partial`/`Unavailable` describe requested property support and
  model admission, without promoting evidence strength or solver qualification.
  The report retains the request: `LocalState` success is conditional on that
  one state, while `Envelope` requires finite ordered corners with identical
  axes and quantity descriptors. The selected claim must remain the same over
  the box, with continuous support. Exact-only curve knots never fill an
  interval; linear curves cannot extend past their knot span. Every competing
  same-key claim's box intersection is queried under the original selection
  policy, detecting conflicts or source changes confined to the interior.
  Queries and receipt replay use the canonical evaluator. No property is borrowed
  from a different material condition to make a candidate complete. Unrelated
  weak evidence does not alter the requested property's status or receipts.
  `ModelRequirement` names an exact law/version and optional model-card pin.
  Only models explicitly associated with that candidate may satisfy it; a
  separately stored model makes the identity known, never associated by name.
  `unknown_models` retains all globally absent law/version requests, distinct
  from candidate-local `Missing` gaps listing actually associated versions.
  Model-only requests are permitted; empty bundles, duplicate requirements,
  blank law names and zero implementation versions refuse before discovery.
  Each returned `DiscoveredModel` retains either the complete original card
  or a `ModelDiscoveryGap`: absent/mismatched pin, missing version, unsupported
  domain, competing calibration, original registry refusal, or narrower built-
  node support. Both card and implementation boxes must contain the local point
  or both envelope corners with exact axis semantics. Overlapping model sources,
  including interior-only overlaps, require an explicit pin; pins never grant
  missing support. Disjoint irrelevant models cannot create ambiguity. No
  parameters are fused or models stitched together across an unsupported range.
  Successful selection calls the supplied registry's real factory through
  `LawRegistry::instantiate`, preserving card/law/state-schema/admission checks.
  Merely having metadata in the store or an implementation registered is
  insufficient. Returned cards preserve `RequiresDeclaredState` when applicable;
  discovery neither invents initial state nor connects or executes a solver graph.
- `verify_index(pack_id)` — cross-checks every derived row against the
  decoded pack, claims table AND validity table (claim hash, axis,
  bitwise bounds); `FS-MATDB-STORE-INDEX-MISMATCH` names the first
  disagreement.
- `load_pack(pack_id)` — hash-verified decode of the stored canonical
  property-pack bytes; other families return `WrongPackKind`.
- `load_catalog_pack(pack_id)` / `load_by_hash(kind, hash)` — typed canonical
  loads for all families, behind `require_sealed`. Content lookup pins the
  whole artifact, not its nested card/claim hash; it never chooses by a similar
  name. The decoded pack id must equal its stored key. Ordered surfaces,
  material-state identities, law versions, state conventions, sources,
  normalization, and species associations remain inside the verified bytes.
- `canonical_dump` — fixed-order render of the derived tables for the
  bitwise-rebuild proof, including each pack's family and content identity.
- `StoreError` — stable `FS-MATDB-STORE-*` codes.

## Invariants

1. Evaluation parity BY CONSTRUCTION, and asserted: the store's answer
   is bitwise the in-memory answer, receipt content-hash included.
2. Index tampering cannot poison an answer: evaluation never reads the
   index (proven by the tamper test — a doctored `scalar_value` row is
   caught by `verify_index` while `evaluate` still returns the exact
   pack value); tampered canonical bytes die at the content hash.
3. Staleness fails closed: after any post-seal ingest, every surface
   refuses until an explicit re-seal.
4. Rebuild determinism (G5): two stores built from the same packs have
   identical corpus digests and identical canonical dumps.
5. Absent data is a named refusal — never a fabricated row or a
   silent empty set: a property name the corpus has never seen is
   `FS-MATDB-STORE-UNKNOWN-PROPERTY` from `materials_with`/`valid_at`
   (an empty result on a KNOWN property remains a legitimate empty
   set), and fs-matdb's own `UnknownProperty`/extrapolation refusals
   pass through evaluation unchanged, per the population strategy in
   `docs/MATERIAL_PROPERTY_TAXONOMY.md`.
6. Ingest is atomic: an induced mid-transaction failure rolls back every
   pack, claim, and validity row in the bundle together. Tests exercise
   both a missing validity table and an existing-id conflict after all
   five families have been written, then successfully retry.

## Error model

Typed `StoreError`; fs-matdb refusals pass through unchanged; SQL
driver failures carry their context. No silent degradation anywhere.

## Determinism class

Deterministic ingest order (ClaimSet canonical iteration), fixed DDL,
canonical digest fold; bitwise-identical rebuild asserted.
Model discovery additionally depends on the caller's exact registry and its
factory determinism. Reports do not serialize or authenticate that registry;
consumers must repeat admission when binding the selected content-pinned card.

## Cancellation behavior

Synchronous short-running statements via the fsqlite sync API; bulk
ingest is caller-chunkable per pack or atomic bundle. No `Cx` integration (workspace
`frankensim-ccmn` effort).
Model discovery uses synchronous law construction/admission under the existing
law-node contract, not a time-stepping workload or a new parallel executor.

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None (fsqlite's `async-api` feature is a dependency detail).

## Conformance tests

`tests/store.rs`: discovery surfaces (per-material, per-property with
range, validity-window); evaluation parity + refusal passthrough;
staleness fail-closed + deliberate re-seal; index tamper detection
(claims AND validity rows) + poison-proof evaluation +
corrupted-bytes refusal; bitwise rebuild; duplicate refusal; atomic
rollback of a failed ingest (plus the executed fact that an
unlicensed claim is UNREPRESENTABLE — `ClaimSet::insert_claim`
refuses it upstream, making the store's admission pre-pass
defense-in-depth); unknown-property named refusal. Five-family file-backed
round trips include a named body and its ordered interface, exact hash lookup,
and direct/store evaluation receipt parity. Additional G0/G4/G5 tests cover
wrong family/hash/wire version, atomic bundle rollback with a preserved seal,
ingest-order-independent rebuilds, and v1 migration with deliberate resealing.
Synthetic fixtures prove storage semantics, not physical dataset accuracy.
G0 compound-discovery tests cover all missing heating requirements in one
report, room-temperature versus high-temperature support, local versus envelope
meaning, sparse curve holes, clipped knot spans, interior-only conflicts, typed
axis/quantity refusals, ordered counter-material filtering, request validation,
and stale corpus refusal. G3 adds unrelated weak evidence and checks unchanged
requested-property support and receipts; existing observation precedence is
also exercised over a competing claim's interior intersection.
G0 model-discovery cases use real canonical stores and the public executable
registry with an explicitly synthetic spring adapter. They check mixed/model-
only bundles, all unknown and candidate-local gaps, exact member provenance,
required-initial-state retention, consumer force evaluation after explicit
initialization, unknown implementations, factory parameter refusal, exact
versions/pins, interior ambiguity, irrelevant disjoint models, typed-axis
refusal and narrower implementation domains. Existing property-only checks use
an explicitly empty registry and retain their assertions.

## No-claim boundaries

- The store never compiles TSV sources: pack compilation stays in the
  fail-closed `xtask matdb-pack` path, and ingest accepts only
  already-admitted values of the five existing canonical pack types. The corpus-wide e2e
  (compile all `data/matdb/seed-v1` packs and ingest) lives with the
  xtask tests where the compiler binary exists — recorded follow-up.
- Existing family wire versions are preserved. Material/interface v1 packs
  do not embed model cards, and the store does not infer model/species
  associations from similar names. Material-card v2 uses the L1-owned explicit
  association described above. Model discovery uses only actual associated
  cards and the caller's executable registry; interface v1 still has no model
  association payload. Exact whole-artifact identities resolve through `load_by_hash`.
- Value-range and compound discovery cover the existing scalar and
  one-dimensional curve payloads. Compound envelopes prove declared data
  support and stable source selection, not a physical trajectory, numeric
  admissibility at every future solver state, or physical validity. Model results
  prove only the stated registry admission and declared card/node box support;
  factory correctness, graph wiring, initialization and future execution remain
  separate obligations. Scenario/CLI material selection and actual lead-heating
  missing-input acceptance remain MR13 work.
  Exact state queries are still required when physics consumes a candidate;
  evolving-state domain exits and rollback belong to the runtime coupling.
- The seal is an integrity mechanism, not authentication: it detects
  drift and tampering against the sealed identity, but a hostile party
  who can rewrite BOTH packs and seal defeats it — authenticity needs
  the fs-package/fs-checker trust channel, out of scope here.
- Concurrent writers are out of contract (single-writer usage; the
  underlying FrankenSQLite locking applies but is not part of this
  crate's claims).
