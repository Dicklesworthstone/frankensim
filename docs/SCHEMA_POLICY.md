# Schema Policy — what FrankenSim promises not to break

Bead: `frankensim-extreal-program-f85xj.16.5`. Registry: [`schema-policy.json`](../schema-policy.json).
Gate: `cargo run -p xtask -- check-schemas` (also runs inside `check-all`).

A product's trust compounds through stable interfaces. A research workspace's
velocity comes from breaking them. FrankenSim resolves that tension the only way
that survives contact with a large tree: by naming a **small** set of serialized
schemas that are public promises, and declaring that everything else is
**internal and breakable**.

The list is the promise. There is no implicit public schema. If a serialized
format is not in the frozen set below, no compatibility is offered for it, no
matter how stable it looks or how many consumers have started reading it.

## The frozen set

Eight schemas are promised. Each has a policy record in `schema-policy.json`
carrying its version location, compatibility promise, migration obligation, and
deprecation horizon.

| Id | Owner | Version constant | Current |
|----|-------|------------------|---------|
| `project.fsim` | `fs-project` | `FSIM_VERSION` | 1 |
| `package.format` | `fs-package` | `FORMAT_VERSION` | 9 |
| `checker.protocol` | `fs-checker` | `CHECKER_PROTOCOL_VERSION` | 7 |
| `scenario.ir` | `fs-scenario` | `SCENARIO_IR_VERSION` | 2 |
| `matdb.pack` | `fs-matdb` | `MATDB_PACK_SCHEMA_VERSION` | 1 |
| `ledger.schema` | `fs-ledger` | `SCHEMA_VERSION` | 20 |
| `sbom.source_manifest` | `xtask` | `SCHEMA` | `frankensim-source-manifest-v1` |
| `maturity.registry` | `xtask` | `REGISTRY_SCHEMA` | `frankensim-capability-maturity-v1` |

The registry is authoritative; this table is a reading aid. `check-schemas`
reads each constant out of its declared source file and refuses if the recorded
value has drifted, so the promise cannot describe a format that no longer
exists.

Two of these carry a **lockstep set**: constants in other crates that must move
in the same commit. The evidence-package format binds
`CHECKER_SUPPORTED_PACKAGE_FORMAT`, `CHECKER_DECISION_IDENTITY_VERSION`, and the
crosswalk's `SUPPORTED_PACKAGE_FORMAT`; the FSMATPK family binds the species,
interface, and model pack versions. Lockstep values are verified against source
in the same way.

## Migration obligations

Every frozen schema declares exactly one obligation, and `check-schemas`
enforces the matching evidence.

- **`auto-migration-receipt`** — an older document loads and migrates, emitting
  a receipt that names each applied rule. Requires at least one named
  cross-version migration test, verified to exist. Held by `project.fsim`,
  `scenario.ir`, and `ledger.schema`.
- **`refuse-unmigratable`** — an older document is refused by name with an
  honest error rather than reinterpreted approximately. Requires a named
  refusal/advertisement test. Held by `package.format` and `checker.protocol`.
  This is a real obligation, not an absence of one: silently accepting a stale
  package would be worse than refusing it.
- **`no-predecessor`** — the schema has never been bumped, so there is nothing
  to migrate from. Held by `matdb.pack`, `sbom.source_manifest`, and
  `maturity.registry`.

`no-predecessor` is the load-bearing case. It may cite no evidence, and the
check additionally refuses any `no-predecessor` record whose current version has
moved past its first. A first bump therefore **cannot** land quietly: it fails
`check-all` until the record is rewritten to declare `auto-migration-receipt` or
`refuse-unmigratable` and to cite a real test.

### The honest boundary

Three frozen schemas have no cross-version migration test because they have no
predecessor, and two have a refusal test rather than a migration path. This
document does not claim that all eight schemas have exercised migrations. It
claims that each declares how its predecessor is handled, that the declaration
matches a test that exists, and that the declaration cannot silently expire.

## Compatibility promise

Each record states what a minor and a major bump may change. The general shape:

- A **minor** bump may add optional fields that an older reader can ignore. It
  may not change an existing value, verdict, or canonical hash for unchanged
  input. Some records have no minor channel at all — the evidence-package format
  and the ledger schema each carry a single integer that a consumer binds
  exactly — and say so.
- A **major** bump may remove or retype fields, change canonical bytes, or
  change identity construction. It requires a migration path, a lockstep sweep,
  and a new **deprecation horizon**.

The **deprecation horizon** is per-record and explicit: how long a superseded
version stays supported. For `project.fsim` and `scenario.ir` a predecessor
stays accepted for one major version. For `package.format` and
`checker.protocol` the horizon is deliberately none — a superseded version is
refused immediately, because a mis-read evidence package is a false certificate.
For `ledger.schema` every rung from v1 forward stays migratable, and dropping a
rung requires a decision record.

## Accretion control

The failure mode this gate exists to stop is not a deliberate break. It is
**accretion**: a new serialized format appears in a product-boundary crate, some
consumer starts depending on it, and it becomes accidentally public without
anyone deciding that it should be.

`check-schemas` scans every public integer version constant declared in the
crates named by `accretion_scope` — currently `fs-checker`, `fs-crosswalk`,
`fs-ledger`, `fs-matdb`, `fs-package`, `fs-project`, and `fs-scenario` — and
requires each to be classified exactly once, either in `frozen` (a promise) or
in `not_promised` (internal and breakable, with a stated reason). A new constant
that is neither fails the gate.

`not_promised` reason classes:

| Class | Meaning |
|-------|---------|
| `identity-domain` | a content-hash domain separator, not an interchange format |
| `internal-row` | a database row shape carried by the ledger schema ladder |
| `internal-receipt` | a receipt returned through the API, not exchanged as a file |
| `internal-wire` | an encoding internal to its crate |
| `derived-lockstep` | travels inside a frozen schema and is promised through it |

The check also refuses **stale rows**: a `not_promised` entry naming a constant
that no longer exists is removed, so the list cannot look complete while
covering nothing. Entries are sorted by `(location, constant)` so an addition
reads as a one-line diff.

### What accretion control does not cover

The scan is scoped to the product boundary, by design. The other workspace
crates are research surfaces whose serialized shapes are internal and breakable
under the blanket exclusion above; they are not scanned, and their absence from
the registry is not an oversight. Private constants, non-integer constants, and
computed or aliased constants are also outside the scan — a version that is not
a public integer literal is not a format this policy can bind, and
`check-schemas` refuses rather than guesses when a declared constant fails to
resolve.

## Changing a frozen schema

1. Change the schema and bump its version constant.
2. Update its record in `schema-policy.json` in the **same commit**: the new
   version, any lockstep values, and the migration obligation with its evidence.
3. Add or extend the migration test the obligation names.
4. If the obligation was `no-predecessor`, rewrite it — the gate will refuse it
   otherwise.
5. Run `cargo run -p xtask -- check-schemas`, then `check-all`.

## Adding a new serialized format

If it is internal, add a `not_promised` row with a reason class. That is the
common case and it is cheap.

If it is a promise, add a `frozen` record with a full policy record and expect
to keep it. The set is meant to stay small: every entry is a commitment that
constrains future work, and the value of the frozen set comes from its size.
