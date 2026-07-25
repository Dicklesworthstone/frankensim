# FrankenSim Constellation Governance

Policy version: 1
Owning Bead: `frankensim-extreal-program-f85xj.13.6`
Policy owner: FrankenSim release owner (`jemanuel` in the owning Bead)
Effective date: 2026-07-23
Next scheduled review: 2026-10-23

## Scope and authority

This policy governs the seven repositories pinned by `constellation.lock`:
asupersync, FrankenNetworkX, FrankenNumpy, FrankenPandas, FrankenScipy,
FrankenSQLite, and FrankenTorch. It governs maintainership assumptions,
incident handling, support, and source retention. It does not prove that any
sibling is correct.

The generated [trust-cone assessment](CONSTELLATION_TRUST_CONE.md) is the
authority for measured consumers, API references, risk classes, and verification
gaps. [Claim integrity](CLAIM_INTEGRITY.md) is the authority for deciding when a
sibling defect invalidates a FrankenSim claim. [CI gates](CI_GATES.md) describe
the checks that exist today. In a conflict, executable evidence and an explicit
no-claim boundary take precedence over optimistic policy prose.

`constellation.lock` proves selected content identity. It does not establish
maintainer availability, independent review, vulnerability response, long-term
hosting, or recoverability if a remote disappears.

## Maintainer reality and decision rights

The operational planning assumption is a bus factor of one for each sibling
unless a repository-local, current maintainer registry proves otherwise. No
machine-readable cross-repository maintainer registry exists today, so this
document does not invent additional reviewers or a guaranteed escalation
contact.

The FrankenSim release owner is accountable for:

- deciding whether a sibling revision is admissible to a candidate pin set;
- recording the evidence and no-claim boundaries for that decision;
- stopping a release when required evidence is absent;
- coordinating an upstream report or candidate repair; and
- ensuring that a pin change, its compatibility evidence, and affected
  semantic goldens move together.

The release owner may delegate review, but delegation does not transfer
accountability. A sibling maintainer may propose a fix; that proposal is not
FrankenSim admission evidence by itself. For correctness-critical asupersync and
FrankenSQLite changes, independent review is preferred and its absence remains
an explicit residual risk.

FrankenPandas is `pinned-unused`. It receives lock and availability governance,
but no FrankenSim runtime correctness claim. Its first consumer requires a new
admission decision, measured usage, compatibility coverage, and an updated
trust-cone assessment before use can be described as supported.

## Review cadence

The release owner reviews this policy and the trust-cone assessment:

1. at least once every three calendar months;
2. before each planned release train;
3. after every P0 or P1 sibling incident;
4. when a sibling changes ownership, archival status, license, security posture,
   or distribution location; and
5. before a pinned-unused sibling gains its first consumer.

A review may conclude that no pin should move. Missed review dates do not make a
stale assessment current; they create a governance finding that must be
recorded. This cadence is a review obligation, not an uptime or response-time
service-level agreement.

## Incident classification

Classify the reachable FrankenSim consequence, not the upstream label:

| Class | FrankenSim consequence | Claim-integrity mapping | Default action |
| --- | --- | --- | --- |
| C0 | credible data corruption, security compromise, false certificate, or default-path result that may be wrong | P0 when a public/default claim can be false | stop affected releases and claims; preserve evidence |
| C1 | correctness-critical or availability-critical default path refuses, hangs, leaks, or loses bounded cancellation without evidence of a wrong accepted result | P1 unless the refusal concealed a false public claim | contain the path; open an emergency candidate |
| C2 | optional, feature-gated, or narrow interop regression with an honest no-claim boundary | P1 or P2 according to reachability | disable or retain the boundary; schedule a tested bump |
| C3 | planned-only or pinned-unused surface | P2 documentation/governance finding | do not promote the planned surface |

Uncertainty moves classification upward until evidence narrows it. A false
certificate or false durability statement is more serious than an explicit
refusal. An upstream correctness bug is not automatically a FrankenSim
claim-integrity incident: the defect must be reachable from a cited FrankenSim
surface or must have contaminated retained evidence.

## Incident-response protocol

Every suspected C0-C2 sibling incident follows this order:

1. **Report and preserve.** Record the reporter, observed version and Git head,
   exact FrankenSim root identity, affected commands/artifacts, and the smallest
   retained reproducer. Preserve the original database, logs, receipts, and
   hashes. Do not repair the only copy.
2. **Triage reachability.** Map the defect to measured consumers in the
   trust-cone assessment. Determine whether it affects runtime, test-oracle, or
   pinned-unused surfaces and whether any public claim or release evidence cited
   the affected result.
3. **Contain.** Stop affected release or evidence lanes. Prefer refusal,
   feature disablement, pin retention, evidence revocation, and explicit
   no-claim language. Never weaken a trigger, checksum, certificate, or
   cancellation contract merely to regain green status.
4. **Coordinate upstream.** Open a sibling issue or private security report with
   the minimal reproducer and impact statement. Public records must be redacted
   when disclosure would expose a live exploit or sensitive user data.
5. **Prepare a candidate pin.** Keep `main` as the only branch. A candidate is a
   proposed `constellation.lock` state plus the upstream change and release
   note; it is not admitted merely because focused upstream tests pass.
6. **Run the emergency compatibility train.** Apply the E13.4 compatibility
   matrix to the candidate, including the affected load-bearing surface,
   same-snapshot DSR evidence, and all golden/identity consequences. If E13.4 is
   not live, the release remains stopped or explicitly degraded; manual probes
   cannot be promoted into a complete-train claim.
7. **Land or reject atomically.** Land the pin, compatibility evidence,
   changelog/severity note, and required golden movements in one reviewed
   commit, or retain the old pin and rejection record. Emergency status does not
   waive tests.
8. **Publish and archive.** Once E13.3 exists, produce and verify the vendored
   source bundle and attach the incident disposition. Until then, the ordinary
   DSR artifact and sibling-layout bootstrap are not called a self-contained
   archival bundle.
9. **Correct the record.** Tombstone or supersede contaminated claims,
   certificates, ledgers, or releases without erasing the original evidence.
   Record who was notified and which later artifact is authoritative.
10. **Review.** Retain a post-incident review or tabletop transcript, update
    regression coverage, and revise risk/review priorities if the event exposed
    a broader boundary.

There is no guaranteed response or fix time. The control is fail-closed release
authority: absent evidence blocks or narrows a claim.

## Release trains and emergency updates

Routine pin movement is event-driven, not automatic dependency churn. E13.4
owns the executable compatibility suite and full release-train protocol. Until
that Bead lands, no document may claim a current pin bump passed the complete
cross-repository train.

### The compatibility registry and the bump gate

`fs_govern::compatibility` registers, per sibling, the load-bearing claims
FrankenSim depends on and the FrankenSim tests that exercise them, and
adjudicates a proposed pin change. `cargo run -p xtask -- compatibility-report`
names which siblings are off-pin and prints the exact selector each one
requires. It is a report rather than a `check-all` gate, matching the existing
decision to keep `check-constellation` out of `check-all`: a pin drift is a
release-train event, not a source-tree defect.

The registry deliberately records **uncovered** surfaces as well as covered
ones, each with the reason it carries no test. An uncovered surface is a
standing gap that the report shows every run; it is never a blank.

### The bump protocol

A pin change is adjudicated by `evaluate_bump`, which is fail-closed. A bump is
admitted only when all of the following hold, and refusal reports every failing
condition at once rather than one per attempt:

1. Something actually moved. An empty train cannot be recorded as a successful
   one.
2. Every moved sibling has a registered compatibility surface. A sibling nobody
   registered cannot be adjudicated at all.
3. Every moved sibling's surface carries tests. A surface with no coverage
   cannot supply evidence, so a bump that moves it is refused rather than waved
   through.
4. Every `P1`/`P2` surface reports a green execution, whether or not that
   sibling moved, because one sibling's move can break another's surface.
5. Each required result is EXECUTED. `NotRun` is never a pass, and an executed
   result reporting zero tests is refused as well: a selector that matched
   nothing is not evidence.
6. The golden implication is declared, and this is now ENFORCED rather than
   advisory. A semantic golden is identified as `crate:surface`, so it is
   coupled to a sibling exactly when its owning crate consumes that sibling at
   runtime (`coupled_golden_surfaces`). Declaring "no golden surface" while
   coupled goldens exist is a refusal that names them. Ownership is matched on
   the whole segment before the colon, so `fs-ad` does not capture
   `fs-adjoint`.

Emergency justification is a closed classification, so an out-of-train bump
must name one of: reachable security defect, credible corruption, false
scientific or certificate result, cancellation/durability contract violation,
or a critical sibling becoming unavailable. Convenience, new upstream features,
and version freshness are not merely discouraged — they are unrepresentable in
the type, so they cannot be argued into an emergency.

### Current state (2026-07-24): one rehearsed bump, REJECTED

The protocol has been exercised once, end to end, against the live pins. It was
not a synthetic drill — the constellation had genuinely drifted, and the
transcript below is the real outcome. It is pinned as an executable test
(`rehearsed_bump_2026_07_24_is_refused`) so the verdict cannot quietly change.

Stage 1, candidate. All SEVEN siblings were off their recorded pins, including
a `franken_numpy` MINOR move (`0.1.0` to `0.2.0`). The previously understood
single-sibling drift was an artifact of the aggregate lock hash, which reports
that *a* repository moved without naming it.

Stage 2, suite. Two P1 surfaces were executed against the live pins:

- `asupersync` at `0.3.9@054cff23` — **GREEN, 25/25** (`fs-exec` conformance 14,
  `constellation_smoke` 1, `lease_battery` 10, 0 failed). FrankenSim's bounded
  cancellation, drain, budget-propagation and latency-lane contracts hold on
  the new commit.
- `frankensqlite` at `31fc4a3b` — **DID NOT BUILD**. `fsqlite-btree` is mid
  async-pager migration: seven errors where synchronous callers `?` a future.
  The surface therefore recorded `NotRun`. A dependency that fails to compile
  is an absence of evidence, never an absence of a problem.

Stage 3, adjudication. **REFUSED**, for two independent reasons: the
`frankensqlite` P1 durability surface never executed, and `franken_numpy` moved
across a minor version into a surface with NO compatibility coverage, so no
evidence exists that could admit it.

Stage 4, disposition. The bump is rejected and the recorded lock stands as the
rollback reference. No pin was moved and no train is claimed to have passed.

This is the protocol working as intended: a green surface does not launder the
refusal for the others, and the outcome of a real train may legitimately be
"reject".

Adjudicating a PROPOSED pin set rather than whatever happens to be checked out
uses `compatibility-report --candidate <lock>`. The candidate is parsed by the
same canonical reader as the tracked lock, so a malformed or hash-inconsistent
proposal is refused rather than half-read.

Known limitation: there is no tool that MINTS a candidate lock.
`lock-constellation` writes the lock from the live checkouts and refuses when
any tree is dirty, so today a candidate is produced either by checking the
proposed pins out first or by hand-editing the lock and recomputing its
FNV-1a-64 `lock_hash`. Until a minting command exists, "candidate" and "live"
coincide more often than the protocol intends.

Note also that because every `P1` surface must be green on every train, no bump
can be admitted at all while FrankenSQLite does not compile — including a
proposal that touches only asupersync, for which green evidence already exists.
That is deliberate: a tree whose durability surface cannot be verified is not a
tree in which any pin change should be accepted.

Known limitation: `fs-govern`'s dev-dependencies include `fs-ledger`, so
`cargo test -p fs-govern` cannot run during a FrankenSQLite outage even though
the library itself stays sibling-free and `xtask` keeps building. The adjudicator
is available when it matters; its own test target is not.

## Archival, escrow, and retention

E13.3 owns the vendored, deterministic, content-addressed source bundle. Once
that mechanism is live, every published release must retain:

- FrankenSim and all seven exact sibling source trees;
- `constellation.lock`, the
  [unified structural source manifest](../frankensim-source-manifest.json), and
  toolchain identity;
- bundle and per-tree content roots plus independent verification instructions;
- the relevant DSR logs, before/after snapshots, compatibility verdicts, and
  claim corrections; and
- enough metadata to build without sibling directories or network access.

Published source bundles and incident correction records are retained
indefinitely. At least two independently administered storage locations should
hold each published bundle, and a restore/verify drill should run during the
quarterly review. Loss of one copy is an incident; a hash without retrievable
bytes is not escrow.

Current no-claim boundary: E13.3 is open. Git remotes, local sibling checkouts,
`constellation.lock`, bootstrap provenance, and the configured DSR release
artifact are useful inputs, but they are not an independently escrowed,
self-contained source bundle. FrankenSim therefore does not yet claim
reproducibility forever after sibling-remotes disappear.

## Support horizon

Only the newest admitted release train is presumed eligible for ordinary fixes.
Older pins and bundles are immutable historical evidence; they are not silently
patched. A backport requires an explicit new supported release and the same
admission evidence as any other pin change.

When E13.3 is live, a retained bundle is intended to preserve source
reproducibility, not ongoing platform support, security maintenance, hosted
service availability, or performance on future hardware. Without a verified
bundle, support is limited to the presently retrievable pinned sources and the
honest no-claim boundaries above.

## Retained exercise and policy evidence

The initial retained exercise is
[the July 2026 synthetic FrankenSQLite corruption tabletop](CONSTELLATION_INCIDENT_TABLETOP_2026-07.md).
It deliberately stops at missing E13.4 and E13.3 controls rather than
pretending that a focused repair, emergency train, or archival bundle exists.

`cargo run -p xtask -- check-constellation-assessment` checks this policy's
required sections, its cross-references, the tabletop's required stages, and
the generated trust-cone artifacts. That check proves document presence and
consistency only; it does not prove maintainers are available or that an
incident was operationally resolved.
