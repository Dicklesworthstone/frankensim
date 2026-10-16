# Constellation Incident Tabletop: Synthetic FrankenSQLite Corruption

Exercise ID: `FS-CONSTELLATION-TT-2026-07`
Owning Bead: `frankensim-extreal-program-f85xj.13.6`
Exercise date: 2026-07-23
Scenario status: synthetic; no actual corruption was observed
Readiness result: failed closed at missing E13.4 and E13.3 controls

## Purpose and non-evidence boundary

This retained tabletop walks a hypothetical FrankenSQLite corruption report
through report -> triage -> emergency train -> bundle. It tests the decision
path in [Constellation Governance](CONSTELLATION_GOVERNANCE.md); it is not an
incident report and is not evidence that FrankenSQLite corrupted a FrankenSim
database.

The scenario borrows the shape of a real integration hazard—storage semantics
can disagree with fs-ledger assumptions—but injects a fictional recovery defect:
after an acknowledged commit and process restart, one artifact page contains
older bytes while its lineage row retains the newer content hash. No such event
was observed in this exercise.

## Participants and authorities

| Tabletop role | Exercised authority |
| --- | --- |
| Reporter | retained database copy, command transcript, artifact and lineage hashes |
| Incident lead | FrankenSim release owner; classification, containment, and release stop |
| Storage investigator | minimal reproducer and sibling-impact analysis |
| Evidence reviewer | [Claim Integrity](CLAIM_INTEGRITY.md) classification and correction scope |
| Release operator | candidate pin, E13.4 compatibility train, DSR evidence, E13.3 bundle |

The current governance is single-maintainer in practice, so one person may hold
several roles. The transcript records the missing independent review rather than
inventing separation of duties.

## Inject

A release-candidate fs-ledger run reports success, is restarted, and then fails
an integrity re-hash:

- the ledger row names artifact hash `H_new`;
- the reopened blob hashes to `H_old`;
- the retained pre-restart receipt says the transaction committed;
- the database and logs are still available; and
- a release note already describes the candidate as durability-verified, but no
  public release has shipped.

The facilitator additionally states that the candidate uses the currently
pinned FrankenSQLite revision. Whether the fault is in FrankenSQLite,
fs-ledger, storage media, or the test harness is initially unknown.

## Transcript

### T+00 — report and preserve

Decision: stop the release-candidate lane and preserve the database, WAL or
journal files, logs, receipts, root HEAD, `constellation.lock`, sibling heads,
toolchain identity, and exact re-hash command. Work continues on copies. No
vacuum, recovery rewrite, migration, or regenerated golden may touch the
original.

Evidence produced in a real incident: a redacted Bead, immutable artifact
inventory, hashes of every preserved file, and a minimal access list.

### T+15 — triage reachability

The trust-cone assessment identifies FrankenSQLite as an active,
correctness-critical dependency beneath fs-ledger durability and lineage
claims. Because accepted bytes disagree with the acknowledged content hash,
this is classified C0. The candidate durability sentence would become a P0
claim-integrity defect if published without correction.

The team does not yet attribute root cause. “FrankenSQLite corruption” remains
the scenario hypothesis, not a finding. If independent replay instead showed an
fs-ledger hashing defect, ownership would move but the release stop would
remain.

### T+30 — containment and claim correction

Decision: keep the candidate unpublished, mark its durability verdict invalid,
and disable reuse of its receipts. Preserve the prior verified pin. Do not
weaken fs-ledger integrity checks or relabel the mismatch as an availability
failure.

If the candidate had shipped, the same decision would require a public
correction that names affected versions and artifacts, supersedes rather than
deletes old evidence, and explains whether user databases require inspection.
No confidentiality-sensitive database bytes belong in a public issue.

### T+45 — upstream coordination and candidate repair

The storage investigator reduces the mismatch to a deterministic
commit-restart-rehash reproducer, tests it against both the pinned revision and
a proposed upstream repair, and sends the minimal case to the sibling
maintainer. A focused upstream pass is necessary but not sufficient.

Decision: describe the proposed revision as a candidate pin only. Record
upstream test command, revision, machine/ISA, and before/after heads. Do not
change `constellation.lock` yet.

### T+75 — emergency compatibility train

Expected control: run E13.4 against the candidate lock, including
FrankenSQLite transaction, reopen, migration, WAL/checkpoint, crash/durability,
large-blob, fs-ledger time-travel, and integrity-guard surfaces. Bind the run to
one stable root/constellation snapshot and retain exact DSR logs. Check identity
and golden consequences before an atomic pin/evidence/changelog commit.

Observed tabletop readiness: E13.4 is open, so this complete train cannot be
claimed today. Focused probes cannot substitute for it. The correct exercised
decision is to keep the release stopped or explicitly remove the affected
durability claim; emergency status does not authorize an untested pin bump.

### T+105 — release and archival bundle

Expected control after a green train: land the candidate on `main`, build the
release, then use E13.3 to assemble and independently verify the deterministic
vendored source bundle containing FrankenSim, all seven sibling trees,
`constellation.lock`, SBOM/source manifest, toolchain identity, compatibility
evidence, incident disposition, and correction metadata.

Observed tabletop readiness: E13.3 is open. The sibling-layout bootstrap,
ordinary Git remotes, and configured DSR artifact are not that bundle. The
exercise therefore cannot demonstrate offline restoration or independent
escrow, and records a failed-readiness result rather than a fictional bundle
hash.

### T+120 — closeout decision

The incident would remain open until root cause, the candidate disposition,
affected-claim correction, regression coverage, and archival state are
recorded. A rejected candidate is a valid outcome. Closure would cite the
specific compatibility and bundle receipts; this tabletop has none and cannot
serve as those receipts.

## Observations

What the policy handled:

- reachability-based severity instead of accepting an upstream label;
- evidence preservation before repair;
- release and claim containment without deleting history;
- separation of upstream focused proof from FrankenSim admission; and
- explicit refusal when required controls are absent.

Readiness gaps retained by the exercise:

1. E13.4 must make candidate-pin compatibility a runnable, refusing gate.
2. E13.3 must make the source bundle buildable, verifiable, reproducible, and
   restorable without sibling directories.
3. Independent review and private security intake are not guaranteed.
4. No live archive inventory currently proves two retrievable bundle copies.

## Follow-up and review trigger

This transcript must be revisited when either E13.4 or E13.3 lands. Replace the
corresponding hypothetical stage with a retained rehearsal receipt, but preserve
this original failed-readiness record. A real incident requires its own record
and must not overwrite this exercise.
