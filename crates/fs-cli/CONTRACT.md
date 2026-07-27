# CONTRACT: fs-cli

The stable command-line membrane for the Cooling 0.1 product workflow (bead
`frankensim-extreal-program-f85xj.6.2`). The binary is named `frankensim`;
the package is `fs-cli` so the workspace retains its flat `fs-*` crate
convention.

## Purpose and layer

Layer L6 (HELM). `fs-cli` turns command-line arguments and project bytes into
deterministic result records and structured diagnostics. It owns presentation
and exit semantics, not project-schema, solver, report, or package authority.
Those remain with `fs-project`, `fs-session` and the cooling pipeline,
`fs-report`, and `fs-package` respectively.

It also owns the product-level geometry-import orchestration checkpoint for
bead `frankensim-extreal-program-f85xj.6.3`: exact caller-supplied raw bytes
flow through `fs-io` quarantine/promotion, `fs-project` persistent assignment,
and one atomic `fs-ledger` operation. Lower layers remain free of filesystem
and L6 policy.

As of bead `frankensim-extreal-program-f85xj.6.5` (slice 1) it additionally
owns the deterministic solve-orchestration driver: content-derived run
identity, session budgets derived from the project, staged ledgered
execution, cancellation at stage boundaries, durable stage checkpoints
through the fs-exec legacy v1 snapshot envelope, and honest
budget-exceeded/stage-gap terminals.

## Public types and semantics

The v0 grammar is intentionally small:

```text
frankensim [--json] validate <project.fsim|project.json>
frankensim [--json] import <project> <source> <ledger.db> --unit <unit> --max-hole-edges <n>
frankensim [--json] import <project> <source> <ledger.db> --unit <unit> --step-root <id> --target-h <spacing>
frankensim [--json] solve <project.fsim|project.json> <ledger.db>
frankensim [--json] solve --resume <run-id> <ledger.db>
frankensim [--json] report <run-id>
frankensim [--json] package <run-id>
```

`--json` may appear once at any position. Unknown flags, duplicate/missing
operands, mixed mesh/STEP policies, non-integer repair/root values, and
non-finite or non-positive STEP spacing are refused. Project inputs are capped
at 16 MiB before parsing.
`.fsim` selects the canonical s-expression spelling and `.json` the canonical
JSON spelling; unknown extensions are refused rather than guessed.

`validate` invokes the strict `fs-project` reader and all of its recognition
and semantic checks. A successful result reports the canonical project hash,
schema version, zero findings, and the exact authority class
`structural-project-admission`.

`report` and `package` are present in the parser but currently return the
stable `cli-stage-unavailable` refusal naming the producer Bead that must land
before the verb can execute:

- report: `frankensim-extreal-program-f85xj.6.9`;
- package: `frankensim-extreal-program-f85xj.6.10`.

This is a deliberate fail-closed integration seam. Reusing the photovoltaic
skeleton or emitting placeholder artifacts would turn a CLI-shaped mock into
a product claim.

### Geometry import

The library surface exposes `RawGeometryLibrary` and
`import_project_geometry`. A raw-source library binds one caller path/label,
byte payload, length unit, repair or STEP root/sampling policy, and optional
named face groups to the strong identity of an exact canonical project geometry
row. Physical labels are provenance only: `.fsim` continues to identify
imported receipt/content rows, never machine-local paths. STL, OBJ, and PLY use
the quarantine/promote route. Strict triangular faceted STEP uses the
caller-selected positive `FACETED_BREP` root, validated unit ID, and finite
positive target spacing, then assigns on the exact repaired soup returned by
the lower topology/SDF handoff. Its retained wrapper preserves the separate
native decoder and tessellation-import receipts without merging their
authority.

The `import` command executes this path for the product reference shape: exactly
one canonical project geometry row, one raw source file, and one SQLite ledger
destination. Mesh callers must explicitly select a hole-repair cap; STEP
callers must explicitly select root and spacing; both supply the source
coordinate unit. Multi-source and named-group adapter construction remain on
the library surface until a bounded source-manifest grammar is ratified.
Project validation and bounded source-file reading finish before the ledger is
opened; format admission and import refusals are then retained in that ledger.
The command derives deterministic execution seed from the project and caps raw
source bytes by both the import default and declared project memory budget.

On success the library atomically retains:

- exact hostile source bytes as input lineage;
- the exact `fs-io` promoted receipt in both an artifact and the Imports
  extension table;
- a deterministic lossless PLY spelling of the promoted finite mesh;
- each `fs-project` assignment report;
- a complete orchestration receipt and terminal successful operation.

Parse, source-hash, parser-version, promotion, and post-promotion assignment
refusals retain all evidence available at the refusal stage and finish one
terminal error operation. Project-admission, resource-envelope, and
pre-cancellation failures occur before ledger side effects.

### Solve orchestration (slice 1)

The library surface exposes `run_solve`, `resume_solve`, `SolveRunId`,
`SolveStage`, `SolveDriverState`, `SolveOutcome`, and `SolveRefusal`. The
`solve` verb wraps them with bounded project I/O and a real monotonic clock.

Run identity is content-derived before any side effect:
`hash_domain("org.frankensim.fs-cli.solve-run.v1", project canonical hash ||
constellation || workspace || root seed || driver version)`. Budgets travel
inside the project hash, so raising a budget starts a new run whose completed
artifacts still deduplicate by content. Driver semantics version 2 binds the
strict lineage/evidence rules in this contract; version-1 run ids do not
collide with or silently resume under version 2. Every solve operation carries the
32-byte run identity as its ledger `session` value; the run's own operations
are its index — resume and downstream consumers locate evidence through
`visible_op_ids` + session filtering + bounded edge reads, with no ledger
schema extension.

The session capability token is derived from project budgets: wall seconds
from `:solve-time`, memory bytes from `:memory-bytes`, cores pinned to 1
(the slice-1 driver is single-threaded), core-seconds = wall × cores, and one
verb per stage. The driver measures its own wall time (the governor meters
what it is told) and charges the `fs-session` governor after every stage.
Warnings fire at declared consumption fractions (0.5 and 0.9 of the wall
grant); any `Throttled`/`Paused` enforcement stops the run after the current
stage with the honest `budget-exceeded` terminal (exit 6) carrying the
completed stages and their durable artifacts. The partial is never presented
as complete.

Stages execute in pinned order: `import-verify`, `assign`,
`material-resolve`, `flow-network`, `conduction`, `qoi`. Each completed stage
is one atomic ledger operation: frozen stage IR, the Five Explicits, the
stage receipt artifact, lineage links, and the sealed driver state. The first
stage additionally retains the exact canonical project source as input
lineage. In slice 1 `import-verify` streams every retained import artifact
through the ledger's content-hash verifier (row presence is not authority),
checks the frozen import byte/count envelope against retained raw, canonical
PLY, and canonical assignment-report bytes, and `assign` binds declared
targets to that verified evidence; `material-resolve`
(frankensim-hp7tb), `flow-network` (frankensim-frn2i), `conduction`
(frankensim-s93ej), and `qoi` (frankensim-s2l9v) refuse with
`cli-solve-stage-gap` (exit 5) naming their producer bead, and the refusal is
itself retained as a terminal error operation.

Stage checkpoints use the fs-exec **legacy v1** snapshot envelope
deliberately: the v2 envelope's expectation token is in-process-only and its
post-restart authorized resume path is explicitly unfinished (fs-exec
CONTRACT), while the v1 `LegacySnapshotExpectationV1` is reconstructible from
durable fields alone. The sealed driver state (run, project hash, consumption
totals, completed stages) is retained as a `solve-stage-state` artifact; its
ledger content hash and `open_expected` admission prove only bounded codec
integrity. Public `SolveDriverState`/`CompletedStage` construction and a valid
legacy envelope grant no resume authority.

Before opening a governor or publishing progress, resume independently
re-attests every candidate prefix that can tie or extend the longest verified
prefix. After bounded checkpoint decode and shape validation, a candidate
strictly shorter than an already verified prefix is skipped because it cannot
affect longest-prefix selection. Resume discovery requires the exact supported
solve-stage schema, field order, run, driver version, stage/ordinal pairing,
diagnostic state, and logical clocks before reading that operation's bounded
edge set; unrelated successful same-session operations therefore cannot
obscure a valid checkpoint by carrying a wide output set. A stage-shaped
candidate must then have positive, distinct, increasing operation ids; finite
nonnegative one-core totals (`core_s == wall_s`); exact canonical stage rows,
Five Explicits, deterministic main-branch execution context, and typed
operation-content identity sidecars; complete nontruncated edge sets with the
writer's exact edge-count seal; exact canonical stage receipts; one checkpoint
per stage; and a direct predecessor-checkpoint chain. Stage zero must consume
the exact strictly decoded canonical project and one exact import summary. Its
import operation must carry the exact project, Five Explicits, canonical
versioned IR field order and complete limits/policy shapes (including
named-group arrays), deterministic main-branch context, typed content sidecar,
exact typed edges, and strict whole-input-consumed receipt grammar. Unrelated
successful import operations are rejected by the cheap
row/IR/Five-Explicits prefilter before any bounded edge scan. Import evidence
is capped at 255 sources so the complete worst-case `4*N+1` edge set fits the
ledger's 1024-edge scan. Raw artifacts are streamed under their frozen
per-source/aggregate caps; opaque promotion receipts are streamed under a
4 MiB per-artifact cap. PLY and assignment-report evidence that solve must
parse has a 64 MiB per-artifact cap. Before the general PLY parser runs, solve
requires the exact nine-line ASCII header emitted by `fs_io::ply::write_ply`,
canonical unsigned vertex/face counts within the frozen import limits, and a
matching number of newline-terminated body records. Before reading any source
payload, solve preflights the checked total of the summary and every raw,
promotion, PLY, and report artifact against the smaller of the project memory
budget value and a 512 MiB hard admitted-input/work envelope. This is a
deterministic per-candidate resource-admission bound, not a claim that peak
allocations fit within that value: canonical PLY checking can simultaneously
retain the input bytes, decoded soup, and canonical re-encoding. Independently
valid equal-length or longer parallel histories can repeat this bounded
attestation within one resume invocation; no invocation-wide cumulative work
cap or memoization claim is made. A genuine candidate outside these solve
envelopes produces `cli-solve-import-envelope`. Competing independently valid
longest checkpoints are ambiguous and refuse rather than winning by row order.
Only after this complete attestation does resume re-charge recorded
consumption, so the budget continues instead of resetting.

## Output and exit contract

- stdout carries final result records only;
- stderr carries diagnostics and, for solve in JSON mode, one
  `frankensim.cli.solve-progress.v1` JSON-line per completed stage and per
  budget warning;
- JSON mode emits one complete object per line in deterministic field order;
- text mode emits stable `key=value` result rows and `ERROR`/`FIX` diagnostic
  pairs;
- exit `0` is success, `2` usage, `3` input I/O/encoding/size, `4` project or
  run refusal, `5` unavailable product stage or solve stage gap, and `6`
  budget-exceeded honest partial.

Diagnostic codes and fix text are machine-facing compatibility surface.
Human prose may improve without changing a code or exit class.

## Error model

Command-line failures are `Diagnostic` records with stable codes on stderr
plus a result record naming the status. Import failures are
`GeometryImportRefusal` values; solve failures are `SolveRefusal` values
(code, optional stage, what, fix, optional owning-bead dependency, optional
run identity, optional recorded ledger operation). Both refusal families are
retained in the ledger when a ledger is open and the refusal happened inside
a stage; recording failure degrades loudly into the refusal text rather than
silently. Budget exhaustion is not an error: it is a successful
`SolveOutcome` with the `budget-exceeded` status, because the work performed
is real and durable.

Solve refusal codes: `cli-solve-project-invalid`, `cli-solve-budget`,
`cli-solve-session`, `cli-solve-ledger`, `cli-solve-ledger-transaction`,
`cli-solve-import-evidence`, `cli-solve-import-envelope`,
`cli-solve-assignment`, `cli-solve-capability`,
`cli-solve-stage-gap`, `cli-solve-cancelled`, `cli-solve-run-id`,
`cli-solve-unknown-run`, `cli-solve-resume-identity`,
`cli-solve-resume-complete`, `cli-solve-resume-budget`,
`cli-solve-ledger-path`, `cli-solve-ledger-open`, `cli-solve-budget-exceeded`.

## Invariants

- Argument order never changes semantic output except for the documented
  position-independent `--json` flag.
- A successful validation has exactly zero `DecodedProject::findings()` and
  no lenient default or canonicalization receipt, because the CLI uses strict
  readers.
- Every refusal has a non-empty code, message, and suggested fix.
- User-controlled strings are escaped before JSON emission; every JSON record
  is one line.
- No unavailable stage (`report`, `package`) writes a run, report, package,
  checkpoint, or ledger artifact. A solve run writes only through its staged
  operations, and a stage the driver cannot execute refuses with retained
  evidence instead of substituting a skeleton stage.
- Solve run identity is derived before any side effect, printed in every
  result record, and bound as the `session` of every operation the run
  writes; the same project in the same workspace always derives the same run.
- Each completed solve stage commits atomically: stage receipt, sealed driver
  state, exact sealed lineage links, and terminal outcome land together or not at all
  (fs-ledger's crash battery is the durability authority).
- The driver state sealed after stage N lists exactly stages 0..=N with their
  positive, distinct, increasing operation ids and receipt hashes. Resume uses
  the unique longest fully re-attested retained prefix; a malformed candidate
  cannot lend authority to its decoded fields, and equally long valid
  competitors refuse as ambiguous.
- The slice-1 single-core driver records finite nonnegative consumption with
  `consumed_core_s == consumed_wall_s` in every checkpoint. A non-finite or
  backward fresh clock refuses before a checkpoint or lineage seal is
  published.
- A pre-cancelled solve publishes nothing; cancellation between stages leaves
  the completed prefix durable and resumable with bit-identical stage
  evidence relative to an uninterrupted run.
- Budget continuity survives resume: recorded consumption is re-charged
  before any new stage runs, and a resume that would already be exhausted
  refuses rather than resetting the meter.
- Import accepts exactly one source for every exact project geometry row and
  no extras; declaration order, not insertion order, determines retained rows.
- Raw bytes must reproduce the project row's FNV source hook and exact
  `fs-io` parser version before promotion.
- Import refuses a caller-owned ledger transaction so its artifacts,
  extension rows, lineage, and terminal outcome commit or roll back together;
  solve refuses a caller-owned transaction for the same reason.
- Every import and solve operation freezes project-derived units, seed,
  budgets, versions, and capabilities in the ledger Five Explicits. Frozen
  import IR also binds every import/assignment resource limit and, in project
  declaration order, exact source-row identity, source unit, repair cap or
  STEP root and target-spacing bits, and ordered named-group mappings. Caller
  path labels do not enter semantic identity. Frozen solve stage IR binds the
  stage name, ordinal, run identity, project hash, and driver version.
- Solve admission checks frozen raw-source byte caps against streamed actual
  lengths (including checked aggregate addition), assignment request and
  predicate-work caps against the exact project and promoted PLY face counts,
  mesh vertex/face and named-group face-range caps against canonical retained
  PLY, and selected-face totals against canonical retained assignment-report
  rows.
- A STEP success retains both lower receipts and writes/assigns the exact
  repaired soup whose counts and fingerprint appear in the import receipt.
- Caller-supplied named-group face ordinals are never laundered across
  face-removing repair. Duplicate/degenerate removal with non-empty groups
  refuses until an adapter supplies an explicit remap or callers use geometric
  selectors. Orientation-only repair, vertex compaction, and appended hole
  faces preserve existing face ordinals.
- Successful import never truncates assignment results: resolver/report count
  must equal the prepared geometry count.

## Determinism class

Argument parsing, validation formatting, unavailable-stage refusals, and
geometry import identities are pure functions of arguments and input bytes
except for the explicit file/ledger boundaries. They read no clock, RNG,
network, or machine state.

Solve stage receipts, driver-state payloads, run identities, and operation IR
contain no wall-clock values; operations use logical times, so deterministic
stages reproduce identical artifact content hashes across independent runs
and ledgers (conformance-tested). Wall time enters only budget accounting
(the meter and the `wall_s` reporting fields) and is explicitly not identity.
Whether a run ends `completed` or `budget-exceeded` is wall-clock dependent
by design — a budget is a wall phenomenon — but every artifact the run
retains on either path is deterministic.

## Cancellation behavior

Validation is bounded by the 16 MiB CLI input cap but has no asynchronous
cancellation surface.

Geometry import has explicit source-count, per-source-byte, aggregate-byte,
and assignment-work caps. It polls the supplied `fs-exec::Cx` before source
work, per source, before and after promotion, and before ledger publication.
A pre-cancelled attempt publishes nothing. Once the atomic SQLite transaction
begins, the bounded ledger calls finish or roll back; cancellation does not
leave a partial successful operation.

The solve driver owns a caller-supplied `fs-exec::CancelGate` and polls it at
every stage boundary before starting the next stage. Cancellation returns a
`cli-solve-cancelled` refusal naming the resume command; the completed prefix
is already durable, so nothing extra is written. Cancellation inside a stage
body is not yet claimed: slice-1 stage bodies are short bounded verifications,
and the physics stages that need intra-stage polling arrive with their
producer beads. The binary does not yet install an OS signal handler; gate
wiring exists at the library seam and is conformance-tested there.

## Unsafe boundary

No unsafe code.

## Feature flags

No feature flags. Runtime dependencies remain Franken-only.

## Conformance tests

`tests/cli.rs` covers the grammar and all v0 verbs, stable exit classes,
strict validation success, structural findings with fixes, noncanonical input
refusal, JSON escaping/line discipline, import-policy conflict/numeric
refusals, routing of both admitted import policy shapes into bounded project
I/O, the exact producer-Bead refusal for report/package, and the solve
grammar's ledger-operand requirement.

`tests/import.rs` supplies a closed reference tetrahedron and covers G0 retained
lineage, repair of deterministic duplicate/degenerate STL facets, strict
faceted-STEP decoding through topology/SDF handoff, separate nested receipts,
and exact repaired-mesh retention. G3 covers changed source identity, open-mesh
promotion refusal, mis-scaled unit refusal, dangling assignment refusal, and a
clean/dirty re-tessellation pair with identical selector statistics. It also
drills the fail-closed named-group behavior when repair removes faces. G4
covers pre-cancellation with zero publication. G5 covers content-identity
equivalence across independent ledgers and proves that changing exact STEP
sampling bits moves the frozen operation IR and retained summary. Every
recorded case runs the ledger linter.

`tests/solve.rs` covers G0 run-identity determinism and input sensitivity,
the pinned stage order and gap owners, prefix execution with a ledgered
`cli-solve-stage-gap` refusal, and the honest budget-exceeded partial with
warning fractions and refused re-resume. G3 covers missing/mismatched import
evidence, an exact-schema import-summary decoy, a caller-sealed but
non-authoritative state, a canonical higher checkpoint with one substituted
predecessor, `{}` replacements for otherwise exact import limits/policy
objects, isolated per-source and aggregate byte-cap forgeries, a genuine import
outside the project-memory solve envelope, a tiny PLY header with an unknown
zero-property element and a huge declared count, unrelated successful import
and same-run operations wider than the complete edge scan, competing valid
longest checkpoints, unknown run ids, and a corrupted retained project pin
(refused by the ledger's own read-integrity gate before the driver's identity
check — both fail closed). Identity refusals occur before progress or ledger
publication. G4 covers pre-cancellation with zero
publication and between-stage cancellation whose resumed evidence is
bit-identical to an uninterrupted run. G5 covers identical stage-receipt
identities across independent fresh ledgers.

## No-claim boundaries

- `validate` proves only canonical structural and dimensional admissibility.
  It does not prove referenced artifacts or material cards exist, a requested
  capability is installed, the project is solvable, or any physical model is
  valid.
- Geometry import binds exact raw bytes, lower-layer receipts, one promoted
  finite tessellation, assignment reports, and their lineage. The legacy FNV
  hook and caller path/label do not authenticate custody, physical/CAD
  sameness, continuum coverage, units, or topology beyond the retained
  lower-layer claims. The summary's `assignment_table` is strictly encoded and
  bounded as part of the retained bytes but is diagnostic rendering, not
  resume authority; typed assignment-report artifacts carry the retained
  assignment evidence.
- Faceted STEP support is limited to fs-io's pinned triangular root-reachable
  resource subset and estimated SDF handoff. It is not full EXPRESS/AP
  interpretation, representation/unit-context discovery, NURBS/surface
  tessellation, component nesting, self-intersection certification, or
  physical/CAD sameness. Named face groups are caller-supplied labels on the
  promoted soup, not independently certified CAD product-structure identity.
- The presence of report/package in help and parsing is not an implementation
  claim. Until their named authorities land, execution fails before side
  effects.
- The solve verb executes only its slice-1 prefix. A run cannot currently
  complete: `material-resolve`, `flow-network`, `conduction`, and `qoi` are
  typed gaps owned by their named beads, and no end-to-end solve, physics
  answer, QoI, or product-determinism claim is made. `import-verify` proves
  canonical internal row/receipt/lineage structure, retained content hashes,
  and the documented resource-count consistency only. PLY/report parsing
  checks exact writer grammar, project source/unit/subject ordering, and
  bounded counts; it does not recompute the source FNV/parser admission,
  cross-check or replay promotion/STEP repair-root-spacing policy, authenticate
  lower-layer producers, verify mesh/named-group/request/assignment
  fingerprints, re-run named-group or selector classification, or recompute
  reported areas, volumes, and bounds. It therefore does not prove the
  geometry is watertight, meshable, or physically meaningful. `assign` binds
  declared targets to verified import evidence; it does not re-run selector
  resolution.
- The session opened by the driver uses caller-declared capability evidence
  (`open_session_declared`); no external issuer, grant authentication, or
  revocation policy is claimed. Memory consumption is charged as zero — the
  governor does not sample RSS — so governor/runtime-memory enforcement is
  currently inert while its core-seconds and wall-seconds enforcement is real.
  Separately, the solve import byte-total preflight genuinely refuses evidence
  above the declared memory value used as an admitted-input/work cap; it
  neither measures nor proves peak memory use.
- The governor is in-memory per process; exactly-once submission and durable
  session terminals (`Governor::new_durable`) are not yet wired. Budget
  continuity across resume relies on the recorded consumption in the sealed
  driver state, not on a durable governor.
- Stage-level checkpoints are driver state only. Intra-stage solver
  checkpoints (fs-exec v2 envelopes, fs-ledger solver-checkpoint receipts,
  pause acknowledgement) arrive with the physics stage producers. Migration
  of the driver state to the v2 envelope waits on fs-exec's post-restart
  authorized-resume seam, which its CONTRACT names as unfinished.
- Resume attestation proves consistency with this driver's retained canonical
  row/receipt/lineage format. Operation-content sidecars and artifact-edge
  seals detect divergence and post-seal mutation; they are not signatures,
  external producer authentication, or authorization against an actor that can
  manufacture an entirely canonical ledger history through the public ledger
  API. That stronger threat model requires a separately authenticated producer
  capability or signed receipt scheme.
- Deterministic replay is claimed for retained stage artifacts, not for
  enforcement outcomes; kill -9 crash recovery of the ledger is fs-ledger's
  proven claim, and the driver's own kill -9 drill (resume to identical
  content root after a mid-transaction kill) is owed with the first
  physics-bearing stage.
