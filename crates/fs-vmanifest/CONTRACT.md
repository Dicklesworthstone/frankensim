# CONTRACT: fs-vmanifest

The typed VerificationManifest schema for leapfrog G1 claim/evidence
freezes, plus the frozen I01 (multi-field equation/compiler) instance
(bead frankensim-leapfrog-2026-program-i94v.1.1.8.1). Sibling G1-freeze
beads (I02…I15, CP/EM/RL/PD) reuse this schema and add their own instance
modules.

## Purpose and layer

Layer UTIL (schema + seed data; depends only on `fs-blake3`). A G1 freeze
preregisters claims, fixtures, obligations, and waivers BEFORE any
implementation result is inspected, so tolerances, claim wording,
capability scope, and failure policy cannot drift toward favorable
outcomes. Preregistration is not proof; nothing here mints evidence
colors or promotion authority.

## Public types and semantics

- `ClaimSpec` — one preregistered claim: ambition lattice element
  (`Ambition::Solid/Frontier/Moonshot`), polarity (affirmative or
  refutation/falsifier lane), statement, explicit hypotheses, QoI + unit,
  `ToleranceSemantics` acceptance arithmetic, targeted `GauntletTier`,
  `OracleRoute` (identity, declared independence, declared TCB overlap),
  activation/kill criteria, fallback, and the Unknown/no-claim boundary.
- `FixturePin` — a corpus element: authored generator spec (digest
  computed from the exact bytes) or external artifact by 64-hex digest,
  each in a `Partition` (development or held-out; the split is frozen).
- `ObligationRow` — one execution leaf's complete mapping: covered claim
  ids, unit-case classes, G0 generators/laws/shrinkers/seeds, deck ids,
  G3 relations, G4 cancellation/fault schedule, G5 determinism matrix,
  `scripts/e2e/leapfrog/*.sh` entry point, smoke/core/max tier, DSR lane,
  fs-obs event kinds, and the exact replay command.
- `Waiver` — a NAMED skip: subject, narrow reason, owner, retirement
  predicate, expiry/review point, and explicit promotion effect. There
  are no unnamed skips: a claim must be covered by an obligation or a
  waiver, and an obligation deck must resolve to a pin or a waiver.
- `ManifestDraft::freeze` — fail-closed gates in documented order:
  collection caps (BEFORE any deep scan), version, top-level blanks,
  per-component blanks/list caps, duplicate ids, oracle independence
  (production-oracle reuse is a refusal, not a style issue), tolerance
  validity, fixture well-formedness, orphan claim references, orphan
  decks, uncovered claims. Refusals are typed `FreezeRefusal` values.
- `FrozenManifest` — SEALED: no public constructor, no mutating API;
  producers are `freeze` and `amend` only, so holding one proves the
  gates ran on exactly this content. `digest()` is the canonical
  identity the freeze is bound to.
- `FrozenManifest::amend` — the only change path: the successor must
  carry `version + 1` and pass every gate; the `AmendmentRecord` names
  exactly the invalidated descendants (claims and obligation leaves whose
  content identity changed or vanished).
- `claim_digest` / `fixture_digest` / `obligation_digest` /
  `waiver_digest` — per-component canonical identities (length-framed,
  variant-tagged, exact IEEE-754 float bits, domain-separated BLAKE3).
- `i01_draft()` — the I01 (multi-field equation/compiler) instance: 9
  claims (5 baseline [S], 2 [F], 2 [M], including a refutation-polarity
  falsifier lane for the completeness moonshot), 6 fixture pins (2
  held-out), 6 obligation rows covering every claim, 1 waiver for the
  not-yet-licensed external benchmark deck slot. Exposed as a draft so
  consumers freeze it themselves — no panic path hides in a static
  initializer.
- `i02_draft()` — the I02 (machine causalization and structural-index
  compiler) instance on the same schema: 9 claims (5 [S]: incidence,
  deterministic matching/DM/SCC/BLT, scoped index reduction, consistent
  initialization, block plans + repair witnesses; 2 [F]: causal-witness
  minimality/presentation-invariance, hybrid-mode structural
  completeness; 2 [M]: certified hidden-constraint discovery and a
  refutation-polarity globally-optimal-tearing falsifier lane), 7
  fixture pins (2 held-out), 6 obligation rows, 1 waiver.
- `i04_draft()` — the I04 (conservation-defect microscope) instance:
  8 claims (4/2/2 with a counterfactual-completeness falsifier), 8
  fixture pins (2 held-out), 5 obligation rows, 1 waiver.
- `i08_draft()` — the I08 (evidence-budget co-design planner) instance:
  8 claims (5/2/1 with a robust multi-horizon optimality falsifier), 9
  fixture pins (2 held-out), 6 obligation rows, 1 waiver.
- `i15_draft()` — the I15 (executable standards compiler) instance:
  9 claims (5/2/2 with a transitive-impact completeness falsifier and a
  human-authority preservation moonshot), 8 fixture pins (2 held-out), 6
  obligation rows, 1 waiver. LICENSING BOUNDARY: all I15 fixtures are
  synthetic standard-shaped packs; no licensed standard text is embedded;
  real editions enter only through the waived external slot pinned via
  fs-vvreg.

## Invariants

- FAIL-CLOSED FREEZE: a manifest missing any load-bearing field, with
  duplicate ids, a non-independent oracle, an invalid tolerance, a
  malformed fixture, an orphan reference, or an uncovered claim cannot
  freeze; the refusal names the gate.
- SEALED AUTHORITY: `FrozenManifest` is immutable by construction;
  "alter flags after freeze" has no code path. Change is amendment,
  amendment is a new version, and the record names the invalidated
  descendants — nothing else is invalidated.
- CANONICAL IDENTITY: components sort into one total order with content
  tie-breaks; assembly/input order can never move `digest()`; valid
  external hex normalizes to raw bytes so case cannot fork identity;
  every semantic field of every component is mutation-sensitive.
- LATTICE SEPARATION: solid/frontier/moonshot claims are distinct
  elements; a weaker receipt closes its own element and is never
  relabeled as the stronger theorem (the I01 maximal lanes activate only
  after the baseline closes, per their activation fields).
- BOUNDED GATES: collection and per-row list caps are checked before any
  content scan (`MAX_CLAIMS`/`MAX_FIXTURES`/`MAX_OBLIGATIONS`/
  `MAX_WAIVERS`/`MAX_ROW_ITEMS`).

## Error model

Total functions; no panics in library paths. Freeze failures are typed
`FreezeRefusal` values; amendment failures are `AmendmentRefusal` values
(wrong version, or the successor's own refusal). Seed defects surface as
freeze refusals in the conformance battery, not as crashes.

## Determinism class

Fully deterministic: seed data is `const`/static text; digests are
domain-separated BLAKE3 over length-framed canonical bytes with exact
float bit patterns. Byte-stable across runs and thread counts on the
same ISA (the G5 test); cross-ISA stability of the digest is expected
but not yet claimed — see no-claim boundaries.

## Cancellation behavior

None; freezing is synchronous and pure. Cost is `O(n log n)` sorting
plus one content digest per component and `O(rows × list-items)`
reference checks, all bounded by the caps above. Interruption/resume of
manifest ASSEMBLY is covered by the staged-assembly equivalence test:
identity depends only on frozen content, never on assembly history.

## Unsafe boundary

None. Workspace `deny(unsafe_code)` applies.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs`: the I01 draft freezes with the declared lattice
split (5/2/2), held-out partitions, full coverage, and a refutation-
polarity falsifier lane; the freeze-gate battery in documented order
(caps before deep scans, version, blanks, duplicates, oracle reuse,
tolerance cases, malformed fixtures, orphan refs/decks, uncovered claims,
waiver-covered claims accepted); G5 input-order invariance of the
digest; G4 staged-assembly equivalence; the G3 mutation battery
(loosened bands, weakened hypotheses, swapped held-out fixtures move the
identity; production-oracle reuse fails closed; post-freeze alteration
has no code path); per-component mutation sensitivity including external
hex case normalization; amendment semantics (successor version enforced,
defective successors refused, invalidated descendants named exactly,
dropped obligations named); the I02 draft freezes with its own 5/2/2
lattice, 2 held-out partitions, full coverage, refutation-polarity
tearing falsifier, an identity distinct from I01's, and I02 input-order
invariance.

## No-claim boundaries

- A frozen manifest asserts NOTHING about implementation correctness:
  preregistration is not proof, and no evidence color, receipt, or
  promotion authority is minted here.
- The named `scripts/e2e/leapfrog/*.sh` entry points and DSR lanes are
  preregistered locations; their existence and behavior are verified at
  campaign time, not freeze time.
- The manifest digest golden constant is deliberately NOT frozen in this
  crate yet: per `docs/GOLDEN_POLICY.md` a golden pin requires
  committed-tree, two-mode reproduction, scheduled with the batch-verify
  lane. Cross-ISA digest stability is likewise expected but unproven
  until the two-host campaign runs.
- Waivers record discipline; they do not verify that owners discharge
  them by expiry — that policing belongs to the governance lane.
- Amendment invalidation names descendants by claim/obligation identity;
  it does not revoke any downstream artifact itself.
- The I01 held-out fixture SPECS are frozen here, but held-out seed
  DISCIPLINE (that development never peeks) is enforced by the campaign
  runner, not by this crate.
