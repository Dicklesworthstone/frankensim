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

An out-of-train candidate is justified by a reachable security defect, credible
corruption, false scientific/certificate result, cancellation or durability
contract violation, or a critical sibling becoming unavailable. Convenience,
new upstream features, or version freshness are not emergency criteria.

The previous verified lock is the rollback reference. Rollback is a new,
recorded `main` commit; it is never an unrecorded filesystem or history rewrite.
Rollback is insufficient when already-published artifacts are contaminated:
those artifacts also require correction or tombstoning.

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
