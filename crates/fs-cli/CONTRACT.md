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
execution, cancellation at stage boundaries and fixed solve-owned evidence
checkpoints, durable stage checkpoints
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
frankensim [--json] compare <left-run> <right-run> [<ledger.db>]
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

`report` and `package` execute against a run whose seven stages are all
sealed (an incomplete run refuses `cli-report-run-incomplete` and names the
next stage). Before either verb reads a byte of the retained report or
package, it proves the run by **sealed evidence**: the driver-state chain is
re-attested exactly as resume does (row identity, driver version, exact
lineage edges, intact edge seals, recovered card packs reproducing the run
identity), every retained receipt is re-hashed against the hash the sealed
driver state recorded, and the conduction operation's retained outputs are
taken from the sealed operation itself with their kinds checked. Exports never
replay physics; only `solve --resume` re-executes retained stages, because it
continues computing from the verified state. The report result discloses this
as `"verification":"sealed-evidence"` next to
`"authority":"projection-of-retained-receipts"`. MEASURED 2026-09-02: with the
replay, `report`/`package` on the heatsink example cost twice the solve
(3.4 s of a 4.2 s validation was re-solving conduction).

`compare <left-run> <right-run> <ledger.db>` diffs two completed runs of the
same project by their retained receipts and, like the other export verbs,
never replays physics. Both runs go through the resume loader and the report
receipt binding above; the rows are read from the retained QoI receipt (values,
colours, requirement outcomes, the eight budget terms), the material-resolve
receipt (pack-set root and per-kind card identities), and the seven stage
receipts. A stage row is `unchanged (same receipt)` when the hashes agree,
`unchanged (same inputs; differs only by binding keys)` when the two receipts
differ in no top-level key but `run`, `project_hash`, or `import_op` (which
run, project, and import op the receipt belongs to, not what the stage
computed), and `changed` otherwise; the differing keys are listed so the
reader can see why. A design change (another card pack, limit, fidelity, or
geometry import) is a different canonical project, so the two project hashes
are reported side by side with `same_project` and a changed hash is one of the
rows. The result carries `"authority":"projection-of-retained-receipts"` and a
no-claim: a changed value is a change between two Estimated candidates, not a
verified effect. Refusals: no ledger operand (`cli-export-ledger-required`), a
missing ledger (`cli-export-ledger-missing`), an unknown or incomplete run (the
solve loader's own codes), two runs with no common requirement to diff
against because their QoI names or requirement identities differ
(`cli-compare-project-mismatch`, naming both project hashes), and a receipt
whose shape cannot be read (`cli-compare-receipt-shape`). Comparing a run with
itself reports `changed:false` and no row marked changed. MEASURED 2026-09-02:
the hostile foam twin of the reference project (rebound to a 0.04 W/mK card
under the same chemistry key) shows exactly one material card change, a hotter
maximum with a smaller nominal margin, both verdicts `indeterminate`, and five
of seven stages changed (`import-verify` and `assign` differ only by binding
keys); a merely renamed project shows the hash change and a bit-identical QoI.

### Euler cinematic static admission

The `cinematic` command is the stable user membrane for bead
`frankensim-h7xu5.8.1`:

```text
frankensim [--json] cinematic inspect <config.fscine> (--trajectory <artifact>|--run-reduced) [--dry-run]
frankensim [--json] cinematic storyboard <config.fscine> (--trajectory <artifact>|--run-reduced) [--dry-run] <host-resources>
frankensim [--json] cinematic daily <config.fscine> (--trajectory <artifact>|--run-reduced) [--dry-run] <host-resources>
frankensim [--json] cinematic representative-4k-frame <config.fscine> (--trajectory <artifact>|--run-reduced) [--dry-run] <host-resources>
frankensim [--json] cinematic final <config.fscine> (--trajectory <artifact>|--run-reduced) [--dry-run] <host-resources>
frankensim [--json] cinematic resume <config.fscine> (--trajectory <artifact>|--run-reduced) [--dry-run] <host-resources>
frankensim [--json] cinematic verify <config.fscine> --trajectory <artifact> [--dry-run]
frankensim [--json] cinematic mux <config.fscine> --trajectory <artifact> [--dry-run]
```

`<host-resources>` is the complete tuple `--memory-bytes`,
`--free-storage-bytes`, `--wall-time-s`, `--workers`, and
`--paths-per-second`; a partial tuple refuses. Launch modes require it because
their named profile estimate is admitted against caller-declared host
availability. `inspect`, `verify`, and `mux` may omit it and then report
`resource_admission=not-requested`; that is not a host-capability claim.

The command decodes the complete 1 MiB-capped `.fscine` document, binds both
budget references to the exact named profile identity, resolves regular-file
assets relative to the canonical config directory, rejects path escape and
stale bytes, and hashes admitted assets under per-file and aggregate byte
ceilings with cancellation checkpoints. `--trajectory` additionally decodes
the canonical bounded Euler trajectory artifact, checks the configured
trajectory schema version and exact artifact identity, and leaves the artifact
reader to validate its distinct wire-codec version. `--run-reduced` records an
unverified request in the write-free plan; it does not schedule, persist, or
pretend that the campaign has run.

Every successful inspect/plan admission emits deterministic
composition/partition identities, the exact profile and estimate (if
admitted), trajectory verification facts, planned stages, authority, and
no-claim text without physical locator paths. Static admission never creates
the configured artifact root. Until their producer Beads land, non-dry
launch/resume, whole-bundle verify, and mux calls fail closed with the owning
dependency; storyboard/daily execution belongs to `frankensim-h7xu5.8.3`,
representative/final/resume to `.8.2`, verification to `.8.4`, and quarantined
mux to `.8.5`.

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

### Solve orchestration

The library surface exposes `run_solve`, `resume_solve`, `SolveRunId`,
`SolveStage`, `SolveDriverState`, `SolveOutcome`, and `SolveRefusal`. The
`solve` verb wraps them with bounded project I/O and a real monotonic clock.

Run identity is content-derived before any side effect:
`hash_domain("org.frankensim.fs-cli.solve-run.v1", project canonical hash ||
constellation || workspace || root seed || driver version)`. Budgets travel
inside the project hash, so raising a budget starts a new run whose completed
artifacts still deduplicate by content. Driver semantics version 7 binds the
current executable prefix: version 3 added retained card packs and material
resolution, version 4 added flow-network execution, version 5 conditionally
added conduction when the project carries its explicit conduction setup, and
version 6 added exact matching-P1 finite-contact lowering. Version 7
length-frames every variable-width interface-resolution identity field so
distinct declared names cannot collide by shifting delimiters. Older run ids do
not collide with or silently resume under version 7. Every solve
operation carries the 32-byte run identity as its ledger `session` value; the
run's own operations are its index. Resume and import discovery walk globally descending
`visible_op_ids_page_controlled` pages under a first-page high-water mark,
then reconstruct each bounded candidate through
`read_op_fields_controlled`. Retained field hashes are compared with the
fixed operation-content sidecar only after cheap row/IR/Five-Explicits
prefilters; a malformed sidecar on an unrelated row therefore cannot block
discovery. No ledger schema extension is required.

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
stage additionally retains the exact canonical project source and every
admitted card pack as input lineage. `import-verify` streams every retained
import artifact through the ledger's content-hash verifier (row presence is
not authority), checks the frozen import byte/count envelope against retained
raw, canonical PLY, and canonical assignment-report bytes; `assign` binds
declared targets to that verified evidence; and `material-resolve` builds the
complete card library from the admitted packs, resolves every declared
material and interface binding through `fs_project::resolve_bindings` under
`BindingRequirements::thermal_steady_v1()`, and retains each selected claim's
replayable usage receipt as a `solve-material-usage-receipt` artifact.
`flow-network` lowers the declared fan system and enclosure network to its
interval-certified operating point. `conduction` executes only when
`cooling.conduction` explicitly supplies one interior seed per region and a
partition of physical boundary laws. It replays the retained promoted meshes
and assignments, constructs and audits a labeled tetrahedral volume, binds
region labels to material-card conductivity, distributes declared regional
power, resolves declared matching-P1 contact faces against their interface
cards, solves the heterogeneous steady problem, and retains the canonical
stage receipt, a `solve-conduction-solution` temperature-field artifact, and
interface-resolution evidence when contact is present.
Projects without that optional declaration refuse with
`cli-solve-conduction-undeclared` (exit 4), and a project with no
`temperature-max` requirement refuses at `qoi` with `cli-solve-qoi-undeclared`
(exit 4): project defects named at the stage, not stage gaps. `cli-solve-stage-gap`
(exit 5) is reserved for a stage with no producer at all (none today); such
refusals name their producer bead and are retained as
terminal error operations.
The conduction receipt (schema v6) discloses the PUBLISHED mesh's quality
census (tets, vertices, min dihedral, max radius-edge, sliver and flat-tet
counts), the base mesh's flat-tet repair (found / repaired / unrepaired /
rounds — fs-mesh removes zero-volume tets by edge removal before the audit) and
the recovery budgets and statistics that produced it (memory budget, Steiner
cap and depth, segment/facet recovery rows). Radius-edge is disclosed, not
enforced: constrained refinement is bridge plan B2c. The stage REFUSES
(`cli-solve-conduction-mesh-quality`, exit 4) when a zero-volume tet survives
the repair or the smallest dihedral of any solved rung is below 1° (measured:
the reference heatsink sits at 7.3° after repair, 6.4° on its second uniform
rung); slivers between 1° and 5° solve and are disclosed. The Steiner cap
derives from the declared memory budget with the fixture default as a floor,
so identical inputs at the fixture budget mesh identically.

**The uniform h-ladder** (receipt block `ladder`, driver version 12). When the
project declares `solver.fidelity = "ladder"`, the stage solves the audited base
and then up to two further rungs, each one uniform 1→8 refinement of the
labeled complex (fs-mesh CONTRACT item 16: walls split in place with their
parent facet, labels replicated, volume preserved), taken while the next rung
still fits the declared memory budget (`memory_bytes / 256` tets) and the
project declares no interface pairs (they bind to base faces: `stop` says so).
Every other fidelity solves the base once and the block says
`"stop":"fidelity-single-rung"`. The block carries one row per rung (tets,
vertices, `h_m = (volume/tets)^(1/3)`, min dihedral, the QoI stage's functional
`t_max_k` — the nodal maximum over the ThermalLimit region — nonlinear and
Krylov iterations, final residual) and a `richardson` estimate over the last
three rungs by the Eça–Hoekstra procedure: `observed-order` (monotone
differences, observed order in [0.5, 2]: GCI half-width `1.25|f3−f2|/(2^p−1)`
and the extrapolated value), `data-range` (anything else three rungs can show:
`3 × (max − min)` of the QoI over the rungs, a bound on the observed variation
rather than an asymptotic estimate), `converged-exactly`, or `single-rung` /
`two-rungs` (no estimate). The published field, temperature range, energy
balance and quality census are the FINEST rung's. Cost is disclosed, not
hidden: three rungs are about seventy base solves (MEASURED 2026-09-02, debug
build, heatsink example: 685 → 5,480 → 43,840 tets in 0.5 s, 7.8 s and 130 s;
T_max 301.99578 → 301.99615 → 301.99610 K, i.e. `data-range`, half-width
1.1 mK), which is why the ladder is a declared study and not the default. The
linear solve's true-residual gate is two decades below the declared
`tolerance-rel` (floor 1e-13); the crate default of 1e-12 refused the second
rung at 1.37e-12 after 1,256 Krylov iterations.

The QoI receipt (schema v2) carries the ladder's estimate as the
**Discretization term** — the first measured term of the eight — through
fs-airflow's discretization-receipt seam: `"state":"interval"` with the
half-width in kelvin, the term's provenance role
`thermal-qoi-discretization-receipt`, and a `derivation` naming the method
(`richardson-gci`, `eca-hoekstra-data-range` or `bitwise-agreement`), the
ladder status, order, rung count, refinement ratio, safety factor and the
conduction receipt it was read from. The stage still refuses
(`cli-solve-qoi-budget-authority`) if any OTHER term, or this term without a
ladder estimate, carries measured authority; the seven remaining terms are
explicit NO-DATA and the `no_claim` text counts them. The verdict stays
Estimated / indeterminate: one measured term does not close a budget. The
report projects the interval term's magnitude into the uncertainty table and
the QoI's `discretization_error`, and renders the ladder rows and estimate as
its convergence section (fs-ladder's vocabulary via fs-report); a single-rung
run shows no study.
A run id whose ops exist in the ledger but none of which is admissible under the
running driver (another driver version, a stage-receipt schema this driver does
not read, a non-deterministic or off-branch op) refuses with
`cli-solve-run-incompatible` (exit 4), listing how many ops were seen and why
each class was rejected; `cli-solve-unknown-run` is reserved for a run id with
no ops at all.

Card packs reach the run through the repeatable `--materials <pack>` and
`--interfaces <pack>` flags. Admission decodes each pack through its own
canonical `fs-matdb` envelope, orders the set by pack root, collapses
byte-identical duplicates, and refuses when two distinct packs reconstruct
one card — the card library is keyed by card identity, so there is no
last-one-wins path. Caller order and caller paths are not semantic: they
never reach the run identity or any receipt. The canonical set root **is**
bound into the run identity, so a different pack set is a different run
rather than a second answer for one run. Because the packs are retained
against the run's first operation, resume recovers them from the ledger
whichever stage was interrupted, and re-deriving the run identity from the
recovered set is what makes that recovery an attestation rather than a
restatement.

Stage checkpoints use the fs-exec **legacy v1** snapshot envelope
deliberately: the v2 envelope's expectation token is in-process-only and its
post-restart authorized resume path is explicitly unfinished (fs-exec
CONTRACT), while the v1 `LegacySnapshotExpectationV1` is reconstructible from
durable fields alone. The sealed driver state (run, project hash, consumption
totals, completed stages) is retained as a `solve-stage-state` artifact; its
ledger content hash and `open_expected` admission prove only bounded codec
integrity. Public `SolveDriverState`/`CompletedStage` construction and a valid
legacy envelope grant no resume authority.

Before opening a governor or publishing progress, resume collects every
retained driver state of the run (bounded checkpoint decode plus shape
validation), then independently re-attests candidates longest prefix first:
the first fully verified prefix wins, an equally long valid competitor refuses
as ambiguous, and strictly shorter candidates are re-attested only when every
longer one failed, because they cannot otherwise affect longest-prefix
selection. `solve --resume` re-attests by replay (each retained stage is
re-executed and its receipt must be bitwise the retained one); the export
verbs re-attest the same chain by sealed evidence (re-hashed receipts, sealed
outputs) and say so. Resume discovery requires the exact supported
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
parse has a 64 MiB per-artifact cap. Before payload verification, solve
requires the exact nine-line ASCII header emitted by `fs_io::ply::write_ply`
within a 256-byte header window, canonical unsigned vertex/face counts within
the frozen import limits, and a matching number of newline-terminated body
records. Solve then makes two bounded passes over that exact writer grammar:
the first decodes finite vertex coordinates and in-range `u32` triangle
indices; the second requires every numeric token to equal Rust's canonical
writer spelling with exact single spaces and newlines. No generic
`fs_io::ply::read_ply` or full-buffer `write_ply` call remains in solve
verification. Before reading any source payload, solve preflights the checked
total of the summary and every raw,
promotion, PLY, and report artifact against the smaller of the project memory
budget value and a 512 MiB hard admitted-input/work envelope. This is a
deterministic per-candidate resource-admission bound, not a claim that peak
allocations fit within that value. One invocation examines at most 8,192
visible operation IDs in at most 128 64-row pages and charges a single 1 GiB
cumulative work ledger. This is a deterministic conservative byte-equivalent
admission/work measure, not exact CPU accounting: controlled
operation/artifact input bytes are charged once, accepted UTF-8 copies and
direct byte-comparison passes add their processed bytes, derived render/receipt
outputs add their exact output bytes, and fixed page rows, sidecars, edge
items, descriptors, and project/entity/assignment items add documented proxy
charges. If captured history continues beyond either cap, or cumulative work
would exceed 1 GiB, solve returns
`cli-solve-work-envelope`; it does not report a false not-found result or
publish a refusal row that would worsen the same history bound. A genuine
candidate outside the separate per-candidate solve envelope produces
`cli-solve-import-envelope`. Competing independently valid longest checkpoints
are ambiguous and refuse rather than winning by row order. Only after this
complete attestation does resume re-charge recorded consumption, so the budget
continues instead of resetting.

## Output and exit contract

- stdout carries final result records only;
- stderr carries diagnostics and, for solve in JSON mode, one
  `frankensim.cli.solve-progress.v1` JSON-line per completed stage and per
  budget warning;
- JSON mode emits one complete object per line in deterministic field order;
- text mode emits stable `key=value` result rows and `ERROR`/`FIX` diagnostic
  pairs;
- exit `0` is success, `2` usage, `3` input I/O/encoding/size, `4` project or
  run refusal, `5` unavailable product stage or solve stage gap, `6`
  budget-exceeded honest partial, and `130` explicit cancellation observed
  before publishing a cinematic result. The older solve driver retains its
  existing refusal-class cancellation surface; this bead does not silently
  redefine that separate compatibility boundary.

Diagnostic codes and fix text are machine-facing compatibility surface.
Human prose may improve without changing a code or exit class.

Each machine code maps to exactly ONE exit class regardless of which layer
catches it. In particular the card-pack resource caps
(`cli-solve-card-pack-count`, `cli-solve-card-pack-size`) always exit `3`
(input/resource cap): such invocations match the documented grammar — the
card-pack flags are unbounded repetition — so `2` (usage) would misreport
correct syntax, and the cap is enforced at parse, read, and admission layers
without changing the class (bead p63op).

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
`cli-solve-work-envelope`, `cli-solve-assignment`, `cli-solve-capability`,
`cli-solve-stage-gap`, `cli-solve-conduction-undeclared`, `cli-solve-qoi-undeclared`,
`cli-solve-run-incompatible`, `cli-solve-conduction-mesh-quality`,
`cli-solve-cancelled`, `cli-solve-run-id`,
`cli-solve-unknown-run`, `cli-solve-resume-identity`,
`cli-solve-resume-complete`, `cli-solve-resume-budget`,
`cli-solve-ledger-path`, `cli-solve-ledger-open`, `cli-solve-budget-exceeded`,
`cli-solve-usage`, `cli-solve-card-pack-read`, `cli-solve-card-pack-size`,
`cli-solve-card-pack-count`, `cli-solve-card-pack-source`,
`cli-solve-card-pack-decode`, `cli-solve-card-pack-identity`,
`cli-solve-card-pack-conflict`, `cli-solve-material-nonfinite`,
`cli-solve-material-receipt-envelope`.

The `material-resolve` stage additionally propagates `fs_project::bind`'s own
violation codes verbatim (`project-material-card-unknown`,
`project-material-state-mismatch`, `project-material-unbound-region`,
`project-binding-domain-uncovered`, and the rest of that vocabulary). This is
deliberate: flattening every binding cause into one CLI code would destroy
the diagnostic identity the binding layer already defines.

## Invariants

- Argument order never changes semantic output except for the documented
  position-independent `--json` flag.
- A successful validation has exactly zero `DecodedProject::findings()` and
  no lenient default or canonicalization receipt, because the CLI uses strict
  readers.
- Every refusal has a non-empty code, message, and suggested fix.
- User-controlled strings are escaped before JSON emission; every JSON record
  is one line.
- No unavailable `report`, `package`, or `compare` invocation writes a run, report, package,
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
  evidence relative to an uninterrupted run. Cancellation during fresh or
  resume evidence verification returns `cli-solve-cancelled` directly with no
  refusal operation, stage receipt, checkpoint, seal, progress line, or
  terminal-success publication. A zero-stage prefix recommends a fresh retry;
  `--resume` guidance is emitted only after at least one durable stage exists.
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

The solve driver owns a caller-supplied `fs-exec::CancelGate`. It polls at
every stage boundary and throughout solve-owned evidence materialization,
stream verification, incremental UTF-8 validation/copy, JSON cursor
advancement, PLY body scanning, canonical payload decoding, writer-token
verification, named-group face-reference and duplicate validation, retained
project comparison, import-IR canonical comparison, and reconstructed
stage-receipt comparison. Solve-owned artifact readers poll at phase entry
before ledger dispatch, before and after every shared controlled tile, and
again after the controlled reader returns before classifying a missing,
stopped, corrupt, or complete result. Byte scans and comparisons use an
input-byte checkpoint cadence no wider than 64 KiB and repeat the current
checkpoint after inspected work before returning a mismatch or parse error.
It also polls after every successful stage body immediately before clock,
charge, or ledger publication. Cancellation returns `cli-solve-cancelled`
directly; it never passes through the retained-refusal writer, so the already
durable prefix is unchanged.

The 64 KiB value is an input-byte checkpoint cadence, not a wall-clock latency
claim. The controlled fs-ledger reader hashes and delivers one shared
at-most-64-KiB tile before invoking solve; a callback break prevents later
tile hashing/delivery. Its SQL engine may already have materialized the
containing storage row of up to 4 MiB. JSON cursors poll at phase entry, at
each crossed fixed boundary, and at the final current position before their
wrapper classifies success or syntax failure. UTF-8 conversion and output
reservation results, byte-comparison mismatches, named-face range failures,
and PLY line/header/payload failures are likewise classified only after the
post-work checkpoint. Import named-group duplicate checks preserve the
writer's face and group order: they reserve each complete duplicate set once,
scan face references in tiles representing at most 64 KiB of `u32` values,
and bound every solve evidence label at 4 KiB. The enclosing operation IR
inherits fs-ledger's guarded `MAX_OP_IR_BYTES` 1 MiB ceiling.

Allocation and the named duplicate-set insertion itself remain bounded opaque
intervals. Retained project parsing through `fs_project::parse_sexpr` is one
opaque call over at most `MAX_PROJECT_BYTES` (16 MiB), bracketed by phase-entry
and completion checkpoints; solve makes no intra-call latency claim for that
lower-layer parser. Canonical import-project rendering and canonical stage
receipt construction are also captured as bounded opaque calls bracketed by
entry/completion checkpoints, followed by tiled byte comparison. State-envelope
decoding is bounded by 4 MiB and bracketed before its result is classified.
The PLY preflight, payload walk, and writer-spelling wrappers each repeat a
bounded post-result checkpoint, including on short malformed input. Canonical
JSON integer and finite-float token widths are capped before standard-library
numeric conversion; the finite writer bound is 327 bytes, including the
negative smallest subnormal. The binary does not yet install an OS signal
handler; gate wiring exists at the library seam and is conformance-tested
there.

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
publication, every fresh solve-owned evidence phase, a storage-row-spanning
raw artifact, an exact-4-MiB inline promotion receipt stopped after its first
shared 64-KiB controlled tile, extreme finite PLY writer spellings, incremental
UTF-8 and JSON-parser phases, named-group face-range checks, chunked retained
project comparison, resume checkpoint/project/receipt re-attestation, and the
final pre-publication boundary with exact zero publication beyond any already
durable prefix. It also proves a cancelled resume retry reproduces the same stage-receipt
identities as an uninterrupted run. G5 covers identical stage-receipt
identities as an uninterrupted run. G0 additionally lowers the committed
two-region imported fixture through its declared finite-resistance contact,
checks nonzero interface heat transfer and energy closure, and replays the same
receipts on resume. G5 covers identical stage-receipt identities across
independent fresh ledgers. The typed G4 matrices enumerate
61 fresh and 79 resume checkpoints, including entry/intermediate/completion
points around visible-ID pages, controlled candidate fields and text
conversion, operation sidecars, artifact descriptors, edge seals/pages/set
comparisons, the invocation-cached canonical project JSON and Five Explicits,
project validation/identity, entity resolution, source-indexed and stage
assignment derivation, canonical receipt construction, retained-project
parsing, state decoding, both available reconstructed stage receipts, and
order-preserving duplicate work.
Focused unit regressions cover immediate JSON phase entry, final-current
syntax priority, empty/missing evidence dispatch entry, mismatch and invalid
UTF-8 result priority, malformed PLY completion polling, named-face range
priority, exact-4-KiB/cap-plus-one labels, and the exact-1-GiB cumulative-work
boundary whose refused plus-one charge does not advance the meter. A focused
accounting regression proves accepted UTF-8 copies and direct byte comparisons
charge each processed tile once while invalid UTF-8 is not charged as accepted
copy work. Integration regressions cover multi-page unrelated history, the
8,192/8,193 visible-ID boundary, and cancellation/retry through exact-1-MiB
operation fields. A G3
import-IR regression keeps
an intentionally unsorted named-face list unchanged while refusing repeated
faces, repeated group names, and an over-ceiling label with zero stage
publication.

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
- `report`, `package`, and `compare` execute against completed runs and prove
  them by sealed evidence (see the verb section); none of them replays
  physics, and `compare` writes nothing.
- All seven stages execute; a project without an explicit conduction setup or
  without a `temperature-max` requirement refuses at the stage that needs it
  (`cli-solve-conduction-undeclared`, `cli-solve-qoi-undeclared`). The
  conduction path proves only the declared finite-mesh, static
  Dirichlet/Neumann/Robin solve, exact matching-P1 finite contact, and the
  reported algebraic residual and energy balance. It does not authenticate
  source geometry or material claims, establish mesh convergence beyond what a
  declared uniform h-ladder measured on the QoI's nodal maximum (no
  energy-norm or DWR bound), or lower nonmatching or temperature-dependent
  contact. Unsupported interface models refuse rather
  than silently assuming perfect contact. When the project declares a
  schema-v4 `(airflow-convection ...)` law, the stage performs ONE branch's
  conjugate airflow exchange (driver v10, receipt v3 `conjugate` object):
  the coefficient is derived from the named `fs-convection` card at the
  flow-network branch Reynolds number with dry-air properties frozen at
  300 K, the reference temperature is the exponential-law marched air
  temperature over the target's retained exterior area, and the solid/air
  fixed point must pass a kelvin convergence test, a scale-free watt
  balance gate, and an `fs-conduction` decomposition cross-check or the
  stage refuses by name (`cli-solve-conduction-airflow-*`). It does not
  close more than one branch, re-evaluate properties or the coefficient at
  film temperature, model recirculation, buoyancy, redistribution, or
  momentum feedback, propagate the flow bracket through the card, or claim
  experimental validation or a maturity level for the exchange.
- `material-resolve` proves that every declared region and interface resolves
  to an admitted card whose selected claim answers the required property at
  both endpoints of the declared temperature range, and it retains that
  claim's replayable usage receipt. It does **not** authenticate the pack
  producer, validate a claim against any external corpus, narrow or replace a
  claim's stated uncertainty, or turn an `Unstated` uncertainty into a bound.
  A pack that decodes is admissible input, not validated data. Perfect
  contact remains refused by the binding layer, so declaring it still fails
  closed. `import-verify` proves
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
- The solve evidence cadence covers controlled visible-ID page progress,
  controlled operation/artifact tiles, guarded text conversion, prehashed
  sidecar checks, bounded descriptors and edge seals/pages, UTF-8 copy, JSON
  cursors, canonical PLY decoding and writer-token comparison, Five-Explicits
  and canonical byte comparison, project/source identity and validation,
  entity resolution, assignment filtering/derivation, named-face ranges,
  duplicate work, and exact edge-set comparison. SQL may materialize one
  guarded operation field of at most 1 MiB before its at-most-64-KiB callback
  slices, one bounded descriptor, or the bounded rows selected for a 64-ID
  page/1,024-edge read; solve claims cancellation and accounting at the API
  entry/callback/result boundaries, not inside the storage engine.
  JSON/PLY cursor advancement is checkpointed, but its already-admitted input
  bytes are not charged again per parser instruction. Standard-library
  allocation, hash-table insertion, guarded text conversion, `parse_sexpr`
  over at most 16 MiB, `findings`, canonical project/Five-Explicits rendering,
  `resolve_entities`, and each bounded receipt-row formatting call remain
  opaque between their pre/post checkpoints; derived render/receipt bytes and
  accepted UTF-8/direct-comparison passes receive the conservative charges
  described above, but no intra-call wall-time latency or peak-allocation
  claim is made. Canonical project JSON and Five Explicits render once per
  fresh invocation; resume caches them by fully attested retained-project
  source hash and fully attests any distinct source before using its own
  result. Lower-layer work outside these calls is not owned by this driver.
- The order-preserving duplicate sets use a fixed standard-library hasher.
  Their membership result is deterministic and their iteration order never
  enters retained output; the implementation does not request a randomized
  set seed.
- Resume attestation proves consistency with this driver's retained canonical
  row/receipt/lineage format. Operation-content sidecars and artifact-edge
  seals detect divergence and post-seal mutation; they are not signatures,
  external producer authentication, or authorization against an actor that can
  manufacture an entirely canonical ledger history through the public ledger
  API. That stronger authority model requires a separately authenticated
  producer capability or signed receipt scheme.
- Deterministic replay is claimed for retained stage artifacts, not for
  enforcement outcomes; kill -9 crash recovery of the ledger is fs-ledger's
  proven claim, and the driver's own kill -9 drill (resume to identical
  content root after a mid-transaction kill) is owed with the first
  physics-bearing stage.
