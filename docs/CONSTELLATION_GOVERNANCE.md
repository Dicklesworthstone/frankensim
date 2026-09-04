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

Six of the seven siblings carry boundary coverage. Several sit behind
NON-DEFAULT features (`fnp-interop`, `fnx-interop`, `torch-bridge`) and several
are unit tests inside `src/` rather than integration tests, so each registered
test records its kind and its required features — a selector that omits the
feature compiles the boundary out and reports a vacuous pass. Only
`frankenpandas` is uncovered, correctly: it is pinned-unused and FrankenSim
makes no claim that depends on it.

An unmoved sibling is required on every train only when it is high priority AND
sits in the runtime graph. That rule exists because one sibling's move can
break another's surface, which can only propagate through runtime consumption;
a dev-only oracle such as `frankenscipy` is outside it, and is required when it
moves itself.

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

### Current state: four trains run; all seven pins are current

The protocol has now been exercised end to end against live pins three times, and the fourth train settled the sibling the third had to hold. Neither
was a synthetic drill, and both transcripts are pinned as executable tests so
their verdicts cannot quietly change.

**Train 1 (2026-07-24) — refused, evidence missing.** All seven siblings were
off-pin. `asupersync 0.3.9@054cff23` ran GREEN 25/25, but `frankensqlite` could
not build: `fsqlite-btree` failed under the `async-api` feature that
`fs-ledger`'s dev-dependency enables, so the P1 durability surface recorded
`NotRun` and the bump was refused. A dependency that fails to compile is an
absence of evidence, never an absence of a problem.

**Train 2 (2026-07-25) — every surface GREEN, refused on the golden
obligation.** Once the sibling's async-pager migration landed, every registered
surface executed and passed against the live pins:

| sibling | measured |
|---|---|
| asupersync | 25/25 |
| franken_networkx | 8/8 |
| franken_numpy | 2/2 |
| frankenscipy | 2/2 |
| frankensqlite | 5/5 |
| frankentorch | 4/4 |
| frankenpandas | no claim, no runtime consumer — no evidence required |

So on suite evidence the whole seven-sibling drift is compatible. The bump is
nonetheless REFUSED, on the one obligation nobody had discharged: **24 semantic
golden surfaces are owned by crates that consume the movers at runtime**
(`fs-exec:*`, `fs-ledger:*`, `fs-plan:*`, `fs-vskeleton:*`), so the attempt may
not declare `NoGoldenSurface`. Someone must state what happened to those
goldens before the pins move under them. That is the rand_nla mis-pin lesson
made executable rather than remembered.

Two boundaries this exposed and fixed:

- A surface selector must run EXACTLY the registered claim set (`-- --exact`).
  Running whole targets meant an unregistered test could fail a surface for a
  reason the registry never claimed — during train 2 a throughput floor
  (`ledger_008_event_throughput_ledgered`, a perf smoke test) failed on a
  contended shared host while all five registered frankensqlite tests passed.
  Machine contention must not be indistinguishable from sibling
  incompatibility.
- A moved sibling that declares no claim AND has no runtime consumer requires
  no evidence, because a pin move cannot break something nothing depends on.
  Both conditions are required, so "no claims" cannot become a way to dodge
  evidence for something that is actually consumed.

**Train 3 (2026-09-02) — six of seven advanced; frankentorch REFUSED on
absent evidence; golden implication declared.** Owner decision rc-o3 (bead
frankensim-rc-root-q61wp.23) chose to advance the lock to the seven
fast-forwarded sibling heads; bead frankensim-rc-root-q61wp.39 executed it.
The verification host (yto, x86-64, release) had been bootstrapped at the old
pins with the drift gate reporting all seven `on-pin`; its sibling clones were
moved forward to the exact candidate heads (forward-only; the Mac's shared
checkouts were never touched) and the registered selectors ran exactly as
`compatibility-report` prints them:

| sibling | head | measured |
|---|---|---|
| asupersync | 03a0a298 (0.4.9 → 0.4.10) | 6/6 |
| franken_networkx | ab730335 | 4/4 (after the selector fix; 0/0 before) |
| franken_numpy | 9b6b5828 | 2/2 (after the selector fix; 0/0 before) |
| frankenscipy | a75ad6ed | 3/3 |
| frankensqlite | d5c68ea3 (fsqlite 0.3.8 → 0.3.15) | 5/5 |
| frankentorch | HELD at 9627f39c | candidate 74df606b does not compile under nightly-2026-07-06 (kernel-cpu `Mask::select`, bead frankensim-r55qa): 0 tests executed, refused |
| frankenpandas | 38a4b26f | no claim, no runtime consumer — no evidence required |

The golden obligation that refused train 2 is discharged the way train 2's
transcript said it must be: the 24 coupled semantic golden surfaces (the same
registry-derived list) are declared `Unaffected`, backed by the owner-crate
identity batteries executing green at the new heads (fs-exec, fs-ledger,
fs-plan, fs-vskeleton). The sibling-review drills re-ran at the new heads
(`sibling_review_cancellation` 7/7, `sibling_review_durability` 5/5),
and the Journey A spine crates were run there as the bead's step 4: fs-conduction
13/13 across its first three targets (analytic 8/8 once the new Level-A
`thermal-a-heatsink-fin-array-ntu` row got its declared gap entry in the
binding matrix), while fs-project's `fansystem` target and fs-cli did not
execute because the train's own debug-profile builds filled the host's disk
(ENOSPC, linker bus error) — an environment refusal recorded as NotRun, to be
re-run after the host is reclaimed, not a sibling incompatibility. The transcript is pinned as
`train_2026_09_02_exact_heads_are_admitted_after_golden_disposition` in
`crates/fs-govern/tests/compatibility.rs`.

One audit fact this train closes: the committed `Cargo.lock` had already
recorded the drifted sibling versions (asupersync 0.4.10, fsqlite 0.3.15,
frankentorch-kernel-cpu 0.1.1) against pins that named older heads, so `HEAD`
was building against a trusted computing base the lock did not declare. After
this train the declared TCB and the built TCB agree on six of the seven
heads, with frankentorch the one declared exception. The
four fs-ledger GC tests carried in `suite-known-red.json` under bead
frankensim-fsqlite-fk-cascade-ordering-f2jag pass at the new frankensqlite
head and leave that registry in the same change.

Two boundaries this train exposed and fixed:

- **Unit-test selectors matched nothing.** The registry rendered
  `--exact <bare name>` for `#[cfg(test)]` tests inside `src/`, but libtest's
  exact match is against the full path `<module>::tests::<name>`, so
  franken_networkx and franken_numpy each reported `0 passed; 47 filtered
  out` — precisely the vacuous pass the renderer exists to prevent, and it
  had been printing that selector since train 2 without anyone executing it
  exactly. The renderer now emits the qualified path for unit surfaces
  (integration surfaces keep bare names, which ARE their libtest path), with
  a regression test that refuses the bare form.
- **A mover that cannot compile is refused, not skipped.** frankentorch's
  candidate head builds on the sibling's own newer nightly but not on the
  workspace's pinned one. The train moved the other six and left
  frankentorch at its recorded pin; the Mac's frankentorch checkout is
  therefore the one remaining off-pin sibling, visible as `stale-lock` in
  the drift gate until bead frankensim-r55qa resolves it.

**Train 4 (2026-09-04) — the held sibling advanced, and what holding cost.**
frankentorch moved from `9627f39c` to `74df606b`; the other six heads are
unchanged and the lock hash went `4a7ffa243c0cf006` → `0730b0a5295eeab7`.
The advance was forced by a cost train 3 did not price. Cargo records a PATH
dependency's version as found on disk, so holding the pin left the
workstation checkout (kernel-cpu `0.1.1`) and the Linux verification clone
(`0.1.0`) disagreeing, and one shared `Cargo.lock` cannot satisfy `--locked`
on both machines at once. Whichever machine ran cargo last rewrote the line
and broke the other; it flipped four times in two days. Retreating the
workstation was refused by this document's own guardrail — the pin was 619
commits behind and pins never move backwards — so both checkouts were
aligned forward and the pin followed them.

The evidence discipline here is worth stating because the surface test still
cannot run under the pinned nightly. What was verified is what the advance
actually touches, not a substitute claim: nothing in frankensim compiles
this sibling (it enters only as an OPTIONAL dependency of `fs-ad`, off in
the default feature set), `cargo metadata --locked` resolves at the new
head, and on the verification host `check-constellation` reports OK with
all sibling trees clean while `check-constellation-drift` reports all seven
on-pin. A pin advance whose sibling is never built is a bookkeeping change,
and it is recorded as one rather than dressed as a compatibility result.

One property of the drift gate this exposed, sound but easy to misread: its
`no-data` verdict is only reachable when HEAD already differs from the pin
(`pin_relation` returns `OnPin` on equality before any ancestry test), and
`no-data` is not charged as a violation. So `check-constellation-drift`
alone reporting `policy OK` does not establish that every sibling is
on-pin — a shallow or grafted clone reports `no-data` for a genuinely
off-pin sibling. Equality is `check-constellation`'s job, and it does fire:
measured exit 1 against exactly this case before the re-lock.

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
