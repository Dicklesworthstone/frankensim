# fs-flywheel-e2e — CONTRACT

FLYWHEEL CLOSES (bead lmp4.18): the whole-loop e2e harness testing the
addendum's central claim — that speculation (9), incremental recompute
(2), the sheaf-adjudicated merge (10), and tombstones (E) COMPOUND.

## Purpose and layer
Layer L6 (integration harness over the flywheel crates). [F]-gated per
the Ambition-Tag rule.

## Public types and semantics
- `LoopConfig`: per-proposal toggles + the G4 cancel point; `baseline()`
  and `composed()` are the measurement anchors.
- `run_loop(config, iterations, seed) -> LoopReport`: the modeled
  two-agent design-iteration workload on the corpus's CHT-wedge edit
  trace — the exact replay config/requested iterations/seed, total modeled
  cost, completed design iterations (including terminal tombstone/dead
  outcomes), per-stage event trace, accept rate, skips, merge verdicts,
  tombstone blocks, the end-to-end headline node, its retained
  `ColorGraph`, and `trace_hash()` for G5.
- `speedups(iterations, seed)`: isolated per-proposal speedups + the
  composed loop (baseline cost ratios); a zero-work comparison is the neutral,
  finite ratio `1.0`.

## Invariants
- Colors: the modeled baseline is a Verified source minted from a typed
  enclosure origin and admitted by an exact-request fixture capability.
  Every accepted speculation retains distinct agent-specific proposer and
  modeled holdout-dataset Estimated source leaves in the SAME append-only
  `ColorGraph`; those heterogeneous leaves first form an accepted-evidence
  derived node, and the headline advances only from that node. Derived
  identities occupy the reserved `derived:v2:` namespace and cannot be
  re-rooted as sources. Downstream upgrades are refused by the graph write
  gate (laundering-across-the-loop). `LoopReport::headline()` resolves only
  `ColorNode::scientific_color()`: a waived or unresolved node fails closed
  instead of exposing its unverified declaration as the scientific headline.
  `LoopReport::admitted_headline()` (bead 6pf9) goes further: the retained
  lineage itself mints an `fs_evidence::AdmissionReceipt` for the headline
  node and re-verifies it through `LedgerColorAdmissionVerifier`, so a
  positive replayable headline converts to an opaque
  `fs_evidence::AdmittedColor` carrying the exact node lineage, while
  estimated, waived, and replay-divergent headlines refuse with the exact
  gate (`HeadlineAdmissionRefusal`). Positive campaign reporting should
  consume the admitted form; the raw `headline()` color remains a declared
  candidate.
- The compounding measurement follows the review-round-3 protocol:
  isolated AND composed over 5 seeded replays, composed >
  max(isolated) by the stated 1.15x margin, coefficient of variation
  reported and bounded.
- Member-crate refusals (merge conflicts, tombstone blocks, skip
  misses) are OUTCOMES the report records, not errors.

## Error model
The harness has no recoverable error channel for its admitted fixture-scale
campaigns; a fixed corpus/source invariant or member-crate refusal contract
being violated is a bug (panic), not a runtime outcome. Arbitrarily large
public iteration counts are not resource-admitted and can exhaust memory;
this feature-gated measurement harness makes no totality claim under resource
exhaustion.

## Determinism class
Deterministic in `seed`: domain-separated counter-addressed draws are keyed by
logical iteration/agent/operation identity, so toggling one proposal cannot
shift the candidate, speculation, or merge draws owned by another stage. Stage
order and accumulation order are fixed — reports, retained color graphs, and
trace hashes replay bit-for-bit (G5, tested). `trace_hash()` is a versioned,
domain-separated, length-prefixed commitment to every public semantic
report field, every event, every canonical graph row, and every field of
every color node (including origins, admission-policy fingerprints,
demotions, waivers, transitive waiver dependencies, and node hashes).
The replay input config, requested iteration count, and logical seed are
part of that commitment even when cancellation or a zero-work request
would otherwise produce identical outputs.
Member crates are deterministic per their contracts.

## Cancellation behavior
`cancel_after_stages` models the G4 storm at every stage boundary; the
loop unwinds BEFORE any state mutation for the interrupted stage, so a
partial event trace and retained evidence graph are replayable clean
prefixes (tested at five cancel points). Work completed by one or both branches
before the cancellation boundary remains charged to `total_cost`: parallel
merge runs settle the maximum completed branch cost, while serialized runs
settle their sum.

## Unsafe boundary
No `unsafe` anywhere in this crate.

## Feature flags
`flywheel-e2e` ([F], default OFF): gates the entire harness until its
Gauntlet tier + kill metric are green.

## Conformance tests
tests/e2e.rs — fw-001 compounding (margin + CV over 5 replays), fw-002
retained-lineage laundering/re-rooting refusals, fw-003 G5 whole-report
and graph determinism plus field-by-field hash sensitivity, fw-004
G4 cancellation storm, fw-005 telemetry completeness, fw-006 toggle-invariant
logical workload draws, and fw-007 cancellation-prefix cost settlement;
fw-008 defines finite neutral speedups for a zero-work request;
tests/dbg.rs —
config-sweep smoke; tests/phase1_gate.rs — THE PHASE-1 MILESTONE GATE
(xpck.3): skip-yield dashboard live, accept-rate telemetry stratified
by proposer × regime, the merge candidate-remainder diagnostic (<25% in its
synthetic corpus and seeded live fixtures; explicitly not the full unresolved-
merge kill criterion, which also needs retained escalation, refusal, and
type-conflict counts), and
Proposal 9's checkpoint FIXTURE (accept rate > 30% AND median
warm-start savings ≥ 1.5× at the calibrated realistic tolerance, with
a hostile control proving the measurement can fail) — machinery
CONFORMANCE on a deliberately favorable linear θ-family, never
activation evidence (bead sj31i.18): the six-month ACTIVATION claim is
adjudicated exclusively by `activation::adjudicate` against a
preregistered `SamplingFrame` (identity-retained strata, minimum
counts, seed, thresholds) over independent holdout rows with exact
denominators — development rows, duplicate identities, zero warm
denominators, and warm-start-free (unmeasured-savings) strata are typed
refusals, and each gate is an anytime-valid betting e-process at the
preregistered δ, per stratum, so optional stopping and favorable
pooling cannot manufacture a pass (`tests/activation_stats.rs`);
tests/cal.rs — the tolerance-calibration probe behind the gate's
0.05/0.02 choices. `tests/acceptance.rs` exercises a stronger receipt path for
ac_003 through ac_005: the harness retains the verifier's original serialized
receipt, resolves its presented artifact root, independently reruns exact
verifier admission for the bound problem, candidate, tolerance, query/QoI,
units, flux reconstruction, hypotheses, work disposition, and verifier family,
then requires a separate detached-format producer-process attestation fixture
before exposing the lower-owned claim or constructing the solver-free
source-certificate verifier. The attestation subject canonically binds its
purpose, exact receipt root (the complete ordered receipt-sequence root in
ac_004), current-process executable schema/byte count/raw-content hash,
producer crate/version/features, source and dependency cones, workspace
manifest and lock, and toolchain. The fixture sidecar first crosses concrete
exact verifier and policy capabilities; the
domain-owned `PromotionTrustRoot` then re-adjudicates exact pinned verifier and
key-policy byte observations and alone can mint the opaque `PromotionWitness`
required by the claim-bearing type. Missing or truncated sidecars, executable
drift, foreign receipt scopes, verifier/policy substitution, fixed/fake hashes,
producer relabeling, changed interval endpoints, and cross-problem, cross-QoI,
or cross-tolerance replay fail closed. ac_004 binds the executable root,
attestation audit/anchor/policy roots, retained receipt roots, and source/build
inputs into its version-2 whole-path replay commitment.
`tests/phase3_gate.rs` retains the horizon activation
ledger: Proposal A's single typed coverage battery is recomputed and stored as
`Estimated` with exact value bits. Its domain-separated identity binds the
schema and algorithm versions, model, truth dimension, exact range bits, RB
dimensions, concept flag, and exact parameter/tolerance bits. The numeric
kill-floor is observed, but the certification/activation trigger remains
unmet; no disconnected fixture certificate promotes it. The package-root
fixture signature authenticates the five-claim holding-pen record, not their
scientific authority. No wall-time claim is attached to this feature-gated
battery; retained lineage work and shared-host load are part of the observed
test cost.

## No-claim boundaries
- COSTS ARE MODELED UNITS from the corpus's op counts: the loop
  MECHANICS are measured; wall-clock physics compounding lands when
  the wedge's real solvers (CHT vertical) replace the cost model.
- The baseline source capability authenticates one exact, compiled-in
  modeled enclosure request. It proves that the loop exercises typed source
  admission, retained origin/policy provenance, and replay; it is not a
  cryptographic signature, an external trust root, or retained physical
  experiment evidence. Production claims require a real injected origin
  verifier and retained certificate artifact resolution.
- The Proposal-8 query stage is a soft edge (Phase 2 per the polish
  note); the headline-color tail models its admission check.
- The two-agent concurrency model is synchronous round-based; a live
  multi-process swarm trial is the xpck.3 milestone's territory.
- The G4 storm is a deterministic stage-boundary cancel model, not the
  base plan's thread-storm harness (fs-exec owns that).
- ac_003 through ac_005 hash the exact on-disk `current_exe()` bytes through a
  fixed 64-KiB streaming buffer under a 1-GiB cap before producer work and
  recapture them before promotion. Path spelling and read chunking are not
  identity inputs. This proves stable executable-content binding for Cargo's
  `acceptance` test process; it does not prove that those bytes equal the
  already mapped memory image, that a separately deployed `fs-verify` binary
  ran, or that the operating system prevented replacement between reads.
- The detached producer-attestation sidecar, verifier, and key policy used by
  the acceptance battery are deterministic, publicly derivable fixtures. The
  opaque `PromotionWitness` makes the trust-root/type boundary and transplant
  refusals executable, but the fixture is not a cryptographic signature,
  transparency-log inclusion, revocation check, hardware/process attestation,
  or vendor-independent trust root. Production promotion still requires an
  externally produced retained sidecar and independently provisioned verifier
  and key-policy observations while preserving this fail-closed subject and
  promotion-root protocol.
- Other acceptance/phase-gate fixture certificate, signature, and waiver
  verifiers still accept exact in-memory declarations or deterministic,
  publicly derivable test strings. They prove typed capability plumbing,
  policy/request binding, waiver taint, color separation, and root-tamper
  refusal, but do not provide cryptographic authentication, retained physical
  experiment evidence, or vendor-independent third-party review.
