# fs-matdb-store — Contract

## Purpose and layer

Layer L6 (HELM). FrankenSQLite-backed queryable store over compiled
fs-matdb material packs (bead frankensim-oecdy): SQL answers the
discovery questions (which materials carry a property, which
properties a material carries, what is valid at an ambient point);
every actual evaluation goes through the canonical hash-verified pack
bytes and the in-memory fs-matdb evaluator. The pack corpus remains
the only source of truth; the store is a derived, regenerable index
plus a canonical-bytes vault.

## Public types and semantics

- `MaterialStore::open(path)` — file or `":memory:"`; DDL v1 applied
  idempotently (`STORE_SCHEMA_VERSION` in `PRAGMA user_version`).
- `ingest_pack(&NormalizedPack)` — canonical bytes + derived index
  rows (claims, validity axes) in the ClaimSet's canonical order;
  refuses duplicates and empty license/redistribution (the license
  gate survives the store).
- `seal_corpus` / `require_sealed` — a domain-separated BLAKE3 digest
  folded over pack content hashes in pack-id order; EVERY discovery
  and evaluation surface recomputes and compares, refusing
  (`FS-MATDB-STORE-CORPUS-CHANGED`) on drift and
  (`FS-MATDB-STORE-NOT-SEALED`) before the first seal.
- Discovery: `properties_of(pack_id)`, `materials_with(property,
  scalar_range)`, `valid_at(property, axis, value)` (missing axis =
  unconstrained, matching `ValidityDomain` semantics).
- `evaluate(pack_id, property, &QueryPoint, policy)` — decodes the
  stored bytes via `NormalizedPack::from_bytes_verified` and delegates
  to `ClaimSet::query` + `verify_receipt`: the SAME evaluator,
  receipts, and refusal set as direct pack use, passed through
  unchanged (`StoreError::MatDb`).
- `verify_index(pack_id)` — cross-checks every derived row against the
  decoded pack (`FS-MATDB-STORE-INDEX-MISMATCH` names the first
  disagreement).
- `canonical_dump` — fixed-order render of the derived tables for the
  bitwise-rebuild proof.
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
5. Absent data is a named refusal (`UnknownProperty`,
   `NoClaimInDomain` extrapolation) — never a fabricated row, per the
   population strategy in `docs/MATERIAL_PROPERTY_TAXONOMY.md`.

## Error model

Typed `StoreError`; fs-matdb refusals pass through unchanged; SQL
driver failures carry their context. No silent degradation anywhere.

## Determinism class

Deterministic ingest order (ClaimSet canonical iteration), fixed DDL,
canonical digest fold; bitwise-identical rebuild asserted.

## Cancellation behavior

Synchronous short-running statements via the fsqlite sync API; bulk
ingest is caller-chunkable per pack. No `Cx` integration (workspace
`frankensim-ccmn` effort).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None (fsqlite's `async-api` feature is a dependency detail).

## Conformance tests

`tests/store.rs`: discovery surfaces (per-material, per-property with
range, validity-window); evaluation parity + refusal passthrough;
staleness fail-closed + deliberate re-seal; index tamper detection +
poison-proof evaluation + corrupted-bytes refusal; bitwise rebuild;
duplicate refusal.

## No-claim boundaries

- The store never compiles TSV sources: pack compilation stays in the
  fail-closed `xtask matdb-pack` path, and ingest accepts only
  already-verified `NormalizedPack` values. The corpus-wide e2e
  (compile all `data/matdb/seed-v1` packs and ingest) lives with the
  xtask tests where the compiler binary exists — recorded follow-up.
- Material-claims packs (`NormalizedPack`) only in v1: interface,
  model, species, and material-card pack kinds are recorded follow-ups
  (the same vault pattern applies; the dispatcher enum in
  `xtask/src/matdb_pack.rs` is the template).
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
