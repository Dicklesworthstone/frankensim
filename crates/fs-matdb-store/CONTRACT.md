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
  scalar_range)`, `valid_at(property, axis, value)` (missing axis =
  unconstrained, matching `ValidityDomain` semantics). Every `PropertyRow`
  includes `pack_kind`, so an interface claim is not mislabeled as a bulk
  property. Equal-property rows use claim hash as a deterministic tie-break.
- `evaluate(pack_id, property, &QueryPoint, policy)` — decodes the
  stored bytes via the family's hash-verified decoder and delegates
  to `ClaimSet::query` + `verify_receipt`: the SAME evaluator,
  receipts, and refusal set as direct pack use, passed through
  unchanged (`StoreError::MatDb`).
  Named material and ordered-interface packs expose their original nested
  claims; model parameters and species metadata never become synthetic scalar
  claims (`NoPropertyClaims` on property evaluation).
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

## Cancellation behavior

Synchronous short-running statements via the fsqlite sync API; bulk
ingest is caller-chunkable per pack or atomic bundle. No `Cx` integration (workspace
`frankensim-ccmn` effort).

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

## No-claim boundaries

- The store never compiles TSV sources: pack compilation stays in the
  fail-closed `xtask matdb-pack` path, and ingest accepts only
  already-admitted values of the five existing canonical pack types. The corpus-wide e2e
  (compile all `data/matdb/seed-v1` packs and ingest) lives with the
  xtask tests where the compiler binary exists — recorded follow-up.
- Existing family wire versions are preserved. Material/interface v1 packs
  do not embed model cards, and the store does not infer model/species
  associations from similar names. Compound discovery and new cross-pack
  binding formats remain separate work; callers can already resolve exact
  whole-artifact identities with `load_by_hash`.
- Discovery indexes scalar claims' values; curve claims are listed
  (kind = "curve") but range-filtered discovery over curve knots is a
  follow-up.
- The seal is an integrity mechanism, not authentication: it detects
  drift and tampering against the sealed identity, but a hostile party
  who can rewrite BOTH packs and seal defeats it — authenticity needs
  the fs-package/fs-checker trust channel, out of scope here.
- Concurrent writers are out of contract (single-writer usage; the
  underlying FrankenSQLite locking applies but is not part of this
  crate's claims).
