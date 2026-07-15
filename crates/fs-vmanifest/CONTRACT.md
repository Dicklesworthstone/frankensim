# CONTRACT: fs-vmanifest

The typed VerificationManifest schema for leapfrog G1 claim/evidence freezes,
with authored draft constructors for I01, I02, I03, I04, I05, I08, I12, and
I15.
Sibling G1-freeze beads reuse this schema and add their own instance modules.

## Purpose and layer

Layer UTIL (schema + seed data; depends only on `fs-blake3`). A G1 freeze
preregisters claims, fixtures, obligations, and waivers BEFORE any
implementation result is inspected, so tolerances, claim wording,
capability scope, and failure policy cannot drift toward favorable
outcomes. Preregistration is not proof; nothing here mints evidence
colors or promotion authority.

## Public types and semantics

Canonical manifest/obligation encoding is schema v2. V2 introduces unordered,
duplicate-free obligation-set identity and therefore uses new manifest and
obligation hash domains; unchanged claim, fixture, waiver, and authored-spec
component encodings retain their v1 subdomains.

- `ClaimSpec` — one preregistered claim: ambition lattice element
  (`Ambition::Solid/Frontier/Moonshot`), polarity (affirmative or
  refutation/falsifier lane), statement, explicit hypotheses, QoI + unit,
  `ToleranceSemantics` acceptance arithmetic, targeted `GauntletTier`,
  `OracleRoute` (identity, declared independence, declared TCB overlap),
  activation/kill criteria, fallback, and the Unknown/no-claim boundary.
- `ManifestDraft` — mutable assembly authority. Its numeric `version` is
  the sole machine-interpreted manifest-instance revision.
  `FiveExplicits::versions` is opaque, identity-bearing provenance for
  schemas, toolchains, dependencies, and data contracts; freeze hashes but
  deliberately does not semantically parse that prose. Authored instance
  constructors have an additional conformance lint against mirroring the
  numeric revision through a small, explicitly non-exhaustive set of known
  legacy semicolon-field spellings. A future structured-pin schema may remove
  this remaining human-authoring ambiguity without inventing a prose parser.
- `FixturePin` — a corpus element: authored generator spec (digest
  computed from the exact bytes) or external artifact by 64-hex digest,
  each in a `Partition` (development or held-out; the split is frozen).
- `ObligationRow` — one execution leaf's complete mapping: covered claim
  ids, unit-case classes, G0 generators/laws/shrinkers/seeds, deck ids,
  G3 relations, G4 cancellation/fault schedule, G5 determinism matrix,
  `scripts/e2e/leapfrog/*.sh` entry point, smoke/core/max tier, DSR lane,
  fs-obs event kinds, and the exact replay command.
  `claims_covered`, `unit_cases`, `decks`, `g3_relations`, and `obs_events`
  are nonempty duplicate-free sets for identity purposes, so
  presentation-only reordering cannot revoke an otherwise identical leaf.
  The five `canonical_*` accessors expose the prospective lexical view on a
  draft. They sort but do not repair duplicates in an unfrozen raw row,
  because duplicate input is a freeze refusal.
- `FrozenObligationRow` — the immutable accepted projection returned by a
  frozen manifest. Its five owned set fields are lexically sorted during
  freeze after validation; every field is private and exposed only through
  read-only getters, so external code cannot forge or mutate the accepted
  projection. `digest()` returns the exact authored-row component identity.
  Canonical execution/serialization is therefore enforced by the frozen API
  rather than left as a caller convention. This projection changes no
  schema-v2 digest bytes.
- `Waiver` — a NAMED skip: subject, narrow reason, owner, retirement
  predicate, expiry/review point, and explicit promotion effect. There
  are no unnamed skips: a claim must be covered by an obligation or a
  waiver, and an obligation deck must resolve to a pin or a waiver.
- `ManifestDraft::freeze` — fail-closed gates in documented order: all
  collection/per-row-list and cumulative UTF-8 text-byte caps (BEFORE any
  semantic scan), version, top-level blanks, required nonempty claim
  authority, per-component blanks, duplicate ids, oracle independence
  (production-oracle reuse is a refusal, not a style issue), tolerance
  validity, fixture well-formedness, orphan claim references, orphan decks,
  uncovered claims, and unused/misspelled waiver subjects. Each category is a
  whole-manifest pass, so a later-row failure in an earlier category cannot be
  preempted by an earlier-row failure in a later category. Claim
  ids and obligation-leaf ids also share one evidence-id namespace, so
  string-valued amendment invalidation can never collapse two authorities;
  an untyped waiver subject must resolve in exactly one of the claim,
  referenced-deck, or leaf namespaces. Refusals are typed
  `FreezeRefusal` values.
- `FrozenManifest` — SEALED: no public constructor, no mutating API;
  producers are `freeze` and `amend` only, so holding one proves the
  gates ran on exactly this content. `digest()` is the canonical
  identity the freeze is bound to, and `obligations()` exposes only canonical
  `FrozenObligationRow` projections.
- `FrozenManifest::amend` — the only change path: the successor must
  preserve the initiative, carry representable `version + 1`, and pass
  every gate. A predecessor claim id may not become an immediately succeeding
  obligation-leaf id, or vice versa; that adjacent-version authority-kind
  alias is a typed
  `EvidenceKindChanged` refusal. `AmendmentRecord` names exactly the predecessor
  claim and obligation-leaf authority that must be re-earned or explicitly
  rebound: direct claim edits
  reach their producer leaves; fixture/deck, obligation, and waiver
  changes follow reverse dependencies; title, Five Explicits, or
  amendment-rule changes invalidate all predecessor evidence. Filling a
  formerly waived deck slot also propagates. A leaf receipt binds both its
  execution payload and claim mapping: a mapping-only edit invalidates the
  predecessor leaf authority even though its byte-identical execution payload
  may be reused through authenticated amendment lineage. An unrelated sibling
  claim sharing that payload is not invalidated by a claim content edit or
  mapping-only removal/rename, while changing execution semantics or a deck
  invalidates every claim fed by that evidence.
- `claim_digest` / `fixture_digest` / `obligation_digest` /
  `waiver_digest` — per-component canonical identities (length-framed,
  variant-tagged, exact IEEE-754 float bits, domain-separated BLAKE3).
  `ToleranceSemantics` equality follows the same float-bit rule (including
  distinct valid `-0.0`/`+0.0` bounds); valid external digest hex is compared
  by decoded bytes, so `FixtureSource` and `FixturePin` equality normalize
  hex case exactly as their digest does.
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
- `i03_draft()` — the I03 electrostatic/EQS gate: 16 claims split 8
  baseline [S], 2 [F], and 6 [M]; 22 detailed fixture/theorem-card pins
  with 7 stage-local heldouts; 8 obligation rows split 4 Core/4 Max; and
  1 licensed-industrial-pack waiver. Exact cochain algebra is separate
  from numerical convergence. The baseline also freezes the closed-
  system versus grounded capacitance distinction, locally conservative
  current, class-specific dielectric passivity, quantitative EQS-to-
  Maxwell escalation, total-current/work ownership, and constrained
  held-variable force formulas. For a fixed port with dielectric-outward
  normal and current positive into the field, total current is primary and
  satisfies `I=-integral(J_f+partial_t D).n=dot(Q_free)-dot(Q_transfer)`;
  it equals `dot(Q_free)` only for a blocking carrier boundary. Exact
  reciprocal capacitance structure is
  attached to self-adjoint assembly/terminal-adjointness and oriented incidence
  receipts before interval PSD checks. EQS monitors carry explicit norms,
  positive gross-input denominators, and normalized QoI bounds. Fixed-charge
  force uses an affine charge space and gauge quotient, while radial coenergy
  is allowed only on a certified in-domain segment. Independent maximal lanes
  cover governed candidate-first, sealed-simultaneous raw-block commit/reveal
  IID-lot discharge calibration with exact dyadic-atom sampling and
  compatible-set width normalization; support-typed space-charge/aging;
  global-homeomorphism versus a.e.-injectivity electrostriction routes crossed
  with reduced-Schur versus mixed inf-sup stability; and machine-bound
  stationary condensation or complete filtered pronilpotent cyclic-
  L-infinity/BV force naturality, certified signed refinement defects,
  regularized topology-event jumps, and a refutation-polarity counterexample
  search. The latter targets a cardinality-proved exhaustive microgrammar,
  symbolic full-domain theorems, and non-exhaustive adversarial supergrammar,
  with canonicalization of the full decorated object. Manifest version 1
  deliberately grants neither theorem nor exhaustive-search authority from
  prose: pre-candidate successors must freeze complete machine proposition and
  grammar ASTs, definitions, translations, predicates, and rank/unrank proofs.
  Bare chain or cohomology equivalence is explicitly not claimed to preserve
  equilibria or force.
- `i04_draft()` — the I04 (conservation-defect microscope) instance:
  8 claims (4/2/2 with a counterfactual-completeness falsifier), 8
  fixture pins (2 held-out), 5 obligation rows, 1 waiver.
- `i05_draft()` — the I05 HIL/WCET/fixed-point deployment-twin compiler
  gate: 12 claims split 7 baseline `[S]`, 2 `[F]`, and 3 `[M]`; 15 fixture
  and theorem/exhaustiveness-card pins with 6 stage-local held-outs; 7
  execution leaves; and 1 production-target characterization waiver. The
  baseline freezes bounded `DeploymentTwinIR` admission, exact-rational
  fixed-point range/error containment, reproducible safe `no_std` images,
  non-confusable static/compositional/measured/unavailable timing evidence,
  interval-valued multi-clock HIL causality, typed fault/safe-state behavior,
  and directed bounded-horizon target refinement. Frontier lanes attempt
  complete interference ownership and compositional deadline/age bounds while
  reconciling authenticated measurements without relabeling samples as WCET.
  Maximal lanes separately target machine-checked floating-to-binary trace
  refinement and exact worst-case cycles for a complete finite target model,
  while a refutation lane attacks false certificates. Version 1 prose cards
  mint neither theorem nor exhaustive-search authority: successors must first
  freeze executable proposition/definition/target-semantics or finite-state
  grammar/transition/cost/coverage artifacts and independent checkers.
  Generated controller code is safe Rust; unavoidable boot/vector/HAL/
  interrupt target support is a separate content-addressed, capability-
  allowlisted, audited runtime capsule and remains an explicit TCB rather than
  being hidden inside code generation.
- `i08_draft()` — the I08 (evidence-budget co-design planner) instance:
  8 claims (5/2/1 with a robust multi-horizon optimality falsifier), 9
  fixture pins (2 held-out), 6 obligation rows, 1 waiver.
- `i12_draft()` — the I12 hybrid mode-automaton compiler instance: 10
  claims (6/2/2, including the Solid grazing false-certificate
  refutation tripwire), 7 fixture pins (3 held-out), 6 obligation rows,
  and 1 waiver.
- `i15_draft()` — the I15 (executable standards compiler) instance:
  9 claims (5/2/2 with a transitive-impact completeness falsifier and a
  human-authority preservation moonshot), 8 fixture pins (2 held-out), 6
  obligation rows, 1 waiver. LICENSING BOUNDARY: all I15 fixtures are
  synthetic standard-shaped packs; no licensed standard text is embedded;
  real editions enter only through the waived external slot pinned via
  fs-vvreg.

## Invariants

- FAIL-CLOSED FREEZE: a manifest missing any load-bearing field or required
  list item, carrying no claims, with
  duplicate ids, a non-independent oracle, an invalid tolerance, a
  malformed fixture, a duplicate claim mapping, an orphan reference or
  waiver, an uncovered claim, or a claim/leaf evidence-id collision
  cannot freeze; the refusal names the gate.
- SEALED AUTHORITY: `FrozenManifest` is immutable by construction;
  "alter flags after freeze" has no code path. Change is amendment,
  amendment is a same-initiative successor version, exhaustion is a
  typed refusal rather than integer wrap, and the record names the exact
  reverse-dependency invalidation set — nothing else is invalidated.
- VERSION AUTHORITY: the numeric manifest revision advances exactly once
  per amendment and is the only revision field interpreted by the schema.
  Authored constructors are checked only against the documented,
  non-exhaustive known legacy semicolon-field spellings; arbitrary public
  drafts receive no semantic natural-language parsing claim.
  A version-only successor can carry an empty invalidation set; an
  identity-bearing campaign-authority title or campaign-policy change cannot
  masquerade as version-only.
- CANONICAL IDENTITY: components sort into one total order with content
  tie-breaks; assembly/input order and presentation order of declared set
  fields can never move `digest()`; valid external hex normalizes to raw bytes
  so case cannot fork identity; exact floating-point bits govern both digest
  and component equality; `FrozenManifest` equality is digest equality; every
  semantic field of every component is mutation-sensitive. Frozen manifests
  expose accepted obligation rows only through the canonical owned-set
  projection, never raw authored order.
- LATTICE SEPARATION: solid/frontier/moonshot claims are distinct
  elements; a weaker receipt closes its own element and is never
  relabeled as the stronger theorem (the I01 maximal lanes activate only
  after the baseline closes, per their activation fields).
- BOUNDED GATES: collection, per-row list, and cumulative UTF-8 text-byte caps
  are checked before any semantic content scan (`MAX_CLAIMS`/`MAX_FIXTURES`/
  `MAX_OBLIGATIONS`/`MAX_WAIVERS`/`MAX_ROW_ITEMS`/
  `MAX_MANIFEST_TEXT_BYTES`). The byte cap covers text payload only, not an
  overclaimed bound on allocator or canonical-framing overhead.

## Error model

Total functions; no panics in library paths. Freeze failures are typed
`FreezeRefusal` values; amendment failures are `AmendmentRefusal` values
(initiative change, exhausted/wrong version, adjacent-version evidence-kind
alias, or the successor's own refusal). Generic structural seed defects are
freeze refusals; semantic defects specific to an initiative, such as an
incorrect Philox alias/range, are focused conformance or campaign-integrity
failures, never silently normalized and never crash a library path.

## Determinism class

Fully deterministic: seed data and governed-sampling protocol text are
`const`/static; digests are
domain-separated BLAKE3 over length-framed canonical bytes with exact
float bit patterns. Byte-stable across runs and thread counts on the
same ISA (the G5 test); cross-ISA stability of the digest is expected
but not yet claimed — see no-claim boundaries.

## Cancellation behavior

None; freezing is synchronous and pure. Cost is `O(n log n)` sorting plus at
most two content-digest evaluations per claim/fixture/waiver and three per
obligation (manifest grouping, content-tie ordering, then sealed-projection
digest caching), with bounded cross-reference scans that are
worst-case quadratic in capped component/list counts. Amendment propagation
uses a deterministic set of predecessor authority ids, so fan-out cannot
amplify duplicate owned strings before finalization; its reverse-dependency
scans remain within the same caps. Chunked in-memory manifest assembly
identity is covered by clone-boundary equivalence tests: identity depends only
on frozen content, never on append history. This is a G4 precursor, not
durable checkpoint encoding, process restart, corruption recovery, or runtime
request-drain-finalize proof.

## Unsafe boundary

None. Workspace `deny(unsafe_code)` applies.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs`: the I01 draft freezes with the declared lattice
split (5/2/2), held-out partitions, full coverage, and a refutation-
polarity falsifier lane; the freeze-gate battery in documented order
(all four collection caps, all six per-row/list caps, and cumulative text cap
before semantic scans, version, blanks, required
nonempty claims, duplicates, category-global oracle and reference precedence,
tolerance cases, malformed fixtures, orphan refs/decks, uncovered claims,
waiver-covered claims accepted); G5 input-order invariance of the
digest; chunked in-memory assembly identity as a G4 precursor; the G3 mutation battery
(loosened bands, weakened hypotheses, swapped held-out fixtures move the
identity; production-oracle reuse fails closed; post-freeze alteration
has no code path); compile-time-maintained, every-field component mutation
sensitivity for claims/oracles, fixtures, obligations, waivers, and every
tolerance payload/tag, including exact-bit signed-zero distinction; direct
external-hex case-normalized equality; canonical draft accessors plus sealed
canonical frozen obligation projection/digest; amendment semantics (initiative/version/
overflow/evidence-kind gates, defective successors, claim/fixture/filled-waiver-slot/
waiver discharge/obligation/waiver/title/global-policy propagation,
version-only and new-authority preservation, and mapping-set order
invariance); the I02 draft freezes with its own 5/2/2 lattice, 2 held-out
partitions, full coverage,
refutation-polarity tearing falsifier, an identity distinct from I01's,
and I02 input-order invariance; and sibling instance smoke checks cover
I04, I08, I12, and I15.

`tests/i05.rs`: exact I05 7/2/3 lattice, single maximal refutation lane,
15-fixture/6-held-out corpus and once-only seven-leaf claim map; all nine unit
case classes, campaign-policy participation, structured lifecycle events,
entry/replay/DSR bindings, sealed canonical obligation projections, and
stage-specific tiers; fixed-point signedness/scale/rounding/overflow semantics;
the four non-confusable timing-evidence kinds; whole-interval clock alignment;
fault-fidelity/safe-state and regulatory no-claim boundaries; directed target
refinement; theorem/proposition and finite-target/exhaustiveness successor
ratchets; G3 assumption/oracle/tolerance/holdout/partition/policy mutations;
G4 chunked assembly equivalence; G5 input-order invariance; and targeted versus
global amendment invalidation.

`tests/i03.rs`: exact I03 8/2/6 claim lattice and once-only 4-Core/
4-Max leaf mapping; exact nine unit-case classes; content-bound campaign
policy, FailureBundle, and independent-adjudication obligations; exact
22-fixture/7-heldout stage-local partitions with exact deterministic fixture
aliases, inclusive development/Core/Max Philox ranges, collision-free derived
keys, and one public governed statistical-holdout protocol whose realized
nonce/lots/features/predictions/labels remain non-public until adjudication and
whose exact IID authority is conditional on candidate-first semantics, an
authenticated peer-withheld 1024-block join-DAG transcript, an audited honest
uniform component, exact dyadic sampling, candidate-input isolation, and
order/external-id invariance; separately pinned theorem-axiom and formal-
projection no-authority gates; target-only M0 grammar formalization gate;
orthogonal execution/predicate/claim/completeness/integrity/support/observable/
promotion axes, including `DeclaredDivergent` only as an observable disposition;
exact per-leaf deck/event/entry/replay/DSR maps, and one-consumer holdouts;
explicit common cancellation/checkpoint lifecycle events, two declared ISA
families with bitwise comparison confined to identical ISA fingerprints, and
I03-specific chunked clone-boundary assembly identity; per-leaf G4 fault,
no-partial-publication, drain/resume, and durable-retention authority; per-leaf
G5 exact-output authority; G5 order stability; G3 hypothesis/oracle/tolerance/
holdout/policy mutations; and targeted versus global amendment invalidation.

## No-claim boundaries

- A frozen manifest asserts NOTHING about implementation correctness:
  preregistration is not proof, and no evidence color, receipt, or
  promotion authority is minted here.
- The named `scripts/e2e/leapfrog/*.sh` entry points and DSR lanes are
  preregistered locations; their existence and behavior are verified at
  campaign time, not freeze time.
- Authored fixture text is the immutable generator/theorem contract, not
  proof that an executable generator or dataset already exists. Campaign
  receipts must additionally bind implementation/toolchain identity and
  the exact generated artifact bytes before evidence can promote.
- Named independent BEM, Maxwell, interval, and Lean oracle routes are
  preregistered identities and separation contracts. Freeze does not prove
  that those implementations exist, are correct, or satisfy the claimed
  independence; campaign receipts and certifier-adversary lanes must do so.
- The manifest digest golden constant and fixed stream-key/Philox known-answer
  vectors are deliberately NOT frozen in this crate yet: per
  `docs/GOLDEN_POLICY.md` those pins require committed-tree, two-mode
  reproduction scheduled with the batch-verify lane. Current tests lock seed
  syntax, ranges, derivation domain/endianness prose, and collision structure,
  not a known-answer value. Cross-ISA digest stability is likewise expected
  but unproven until the two-host campaign runs.
- Waivers record discipline; they do not verify that owners discharge
  them by expiry — that policing belongs to the governance lane.
- Amendment invalidation names predecessor claim/obligation authority;
  the downstream ledger/governance layer must actually revoke or refuse
  stale receipts carrying those identities. Unchanged component evidence
  may be rebound across a version-only or targeted successor only through
  authenticated amendment lineage and byte-identical component digests;
  whole-manifest digest inequality alone neither revokes everything nor
  authorizes blind reuse.
- An invalidated obligation leaf is a mapping-bound authority whose receipt
  must be reissued or explicitly rebound. `AmendmentRecord` is exact at
  claim/leaf authority granularity; it is not by itself a numerical-kernel
  rerun schedule, and authenticated lineage may reuse byte-identical execution
  payloads where the surviving claim authority permits it.
- Every authored held-out fixture SPEC is frozen here. Deterministic holdouts
  also freeze public replay ranges and therefore carry no statistical
  untouched/IID authority. The discharge max holdout deliberately has no
  public seed. Candidate/model/toolchain and every campaign semantic freeze
  first; only then do three named custodians enter a peer-withheld,
  all-three-fixed commitment phase over exact lot-major 1024x256-bit vectors.
  Strict RFC8032-authenticated transcripts, an audited at-least-one-honest IID
  uniform vector independent of the already-fixed candidate and adversarial
  masks, exact coordinate-wise XOR, 257-bit half-open cumulative intervals for
  the pinned dyadic atom law, a later 256-bit order-challenge beacon that never
  enters sampling, and complete
  candidate-input isolation are all authority-bearing. Exactness does not rely
  on computational hiding or a short-seed PRG. Its IID receipt, external-id
  relabeling/order invariance, access control, one-shot adjudication, and later
  raw-draw reveal must be established by campaign integrity receipts. This
  crate preregisters that protocol but cannot enforce the external governance.
- I03 synthetic discharge holdouts can support calibration evidence only;
  physical/industrial authority remains blocked until the waived,
  independently governed experimental pack is admitted. I03's exact
  topology/force theorems apply only under the complete frozen
  variational or regularized-event premises; generic cohomology
  preservation, bare chain equivalence, refinement, or topology change
  receives no free invariance claim.
- I03 theorem target cards in manifest version 1 are not machine propositions.
  They mint no theorem color until a pre-proof successor freezes canonical
  proposition AST and definition bytes/digests, a total runtime-premise map,
  and a deterministic AST-to-Lean translation with structural round-trip
  checks. The exact axiom policy admits only `propext`, `Quot.sound`, and
  `Classical.choice`; its digest and the complete transitive axiom closure are
  receipt-bound, while `sorryAx`, custom postulates, and unsafe/native-oracle
  shortcuts fail integrity.
- I03's version-1 `M0` prose similarly freezes the ambitious
  16x16x16x16/N=65536 target but has no exhaustive authority. A pre-search
  successor must freeze the full decorated-record grammar, canonical
  encodings/domains, validity/stratum/tag/parameter semantics, explicitly
  encoded event primitives, total enumeration and exclusion order,
  rank/unrank/sharding algorithms, source digests, independent decoder and
  bijection proofs, cost preflight, and Merkle completeness root.
- I05 timing evidence never promotes a `MeasuredSampleMaximum` into a WCET
  upper bound. Static and compositional bounds carry authority only for their
  exact binary, target, clock/power state, task context, microarchitectural and
  interference assumptions. HIL clock containment, injected safe-state
  behavior, and bounded trace agreement are likewise not regulatory safety,
  physical-plant validation, unbounded stability, or authority on another
  target/rig profile.
- I05 version-1 theorem and finite-target cards are ambitious targets, not
  machine propositions or executable exhaustive grammars. A pre-proof
  successor must freeze canonical proposition/definition/decoded-target ASTs,
  total runtime premises, translation and axiom closure. A pre-search successor
  must freeze complete instruction/microstate/input grammar, initial/valid/
  transition/cost semantics, unsupported-state closure, enumeration or sound
  symbolic coverage, quotient obligations, independent decoding, preflight and
  completeness root. Even a model-exact maximum becomes silicon WCET only
  after independent model-fidelity qualification on the exact hardware.
- The adjacent-version `EvidenceKindChanged` guard prevents an immediate
  claim/leaf kind swap. `FrozenManifest` carries no lineage-wide tombstone set;
  a ledger spanning nonadjacent versions must key authority by typed kind plus
  manifest/version digest (or add an authenticated tombstone schema) rather
  than treating a raw string id as globally unique forever.
- The clone-boundary assembly tests do not encode or decode durable
  checkpoints, restart a process, detect corruption, or exercise runtime
  cancellation. Those remain executable G4 campaign obligations.
