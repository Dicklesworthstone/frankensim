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

## Public surface

The v0 grammar is intentionally small:

```text
frankensim [--json] validate <project.fsim|project.json>
frankensim [--json] solve <project.fsim|project.json>
frankensim [--json] solve --resume <run-id>
frankensim [--json] report <run-id>
frankensim [--json] package <run-id>
```

`--json` may appear once at any position. Unknown flags, missing operands, and
extra operands are refused. Project inputs are capped at 16 MiB before parsing.
`.fsim` selects the canonical s-expression spelling and `.json` the canonical
JSON spelling; unknown extensions are refused rather than guessed.

`validate` invokes the strict `fs-project` reader and all of its recognition
and semantic checks. A successful result reports the canonical project hash,
schema version, zero findings, and the exact authority class
`structural-project-admission`.

The remaining verbs are present in the parser but currently return the stable
`cli-stage-unavailable` refusal naming the producer Bead that must land before
the verb can execute:

- solve/resume: `frankensim-extreal-program-f85xj.6.5`;
- report: `frankensim-extreal-program-f85xj.6.9`;
- package: `frankensim-extreal-program-f85xj.6.10`.

This is a deliberate fail-closed integration seam. Reusing the photovoltaic
skeleton or emitting placeholder artifacts would turn a CLI-shaped mock into
a product claim.

The library surface exposes `RawGeometryLibrary` and
`import_project_geometry`. A raw-source library binds one caller path/label,
byte payload, length unit, repair cap, and optional named face groups to the
strong identity of an exact canonical project geometry row. Physical labels
are provenance only: `.fsim` continues to identify imported receipt/content
rows, never machine-local paths. The current checkpoint admits STL, OBJ, and
PLY. Faceted STEP fails with `cli-import-format-unavailable` until the product
surface carries an explicit root entity, length unit, and target sampling
spacing. Import is not yet a v0 command-line verb.

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

## Output and exit contract

- stdout carries final result records only;
- stderr carries diagnostics (and will carry solve progress JSON-lines once
  solve orchestration exists);
- JSON mode emits one complete object per line in deterministic field order;
- text mode emits stable `key=value` result rows and `ERROR`/`FIX` diagnostic
  pairs;
- exit `0` is success, `2` usage, `3` input I/O/encoding/size, `4` project
  refusal, and `5` unavailable product stage.

Diagnostic codes and fix text are machine-facing compatibility surface.
Human prose may improve without changing a code or exit class.

## Invariants

- Argument order never changes semantic output except for the documented
  position-independent `--json` flag.
- A successful validation has exactly zero `DecodedProject::findings()` and
  no lenient default or canonicalization receipt, because the CLI uses strict
  readers.
- Every refusal has a non-empty code, message, and suggested fix.
- User-controlled strings are escaped before JSON emission; every JSON record
  is one line.
- No unavailable stage writes a run, report, package, checkpoint, or ledger
  artifact.
- Import accepts exactly one source for every exact project geometry row and
  no extras; declaration order, not insertion order, determines retained rows.
- Raw bytes must reproduce the project row's FNV source hook and exact
  `fs-io` parser version before promotion.
- Import refuses a caller-owned ledger transaction so its artifacts,
  extension rows, lineage, and terminal outcome commit or roll back together.
- Every import operation freezes project-derived units, seed, budgets,
  versions, and capabilities in the ledger Five Explicits.
- Successful import never truncates assignment results: resolver/report count
  must equal the prepared geometry count.

## Determinism and cancellation

Argument parsing, validation formatting, unavailable-stage refusals, and
geometry import identities are pure functions of arguments and input bytes
except for the explicit file/ledger boundaries. They read no clock, RNG,
network, or machine state. Validation is bounded by the 16 MiB CLI input cap
but has no asynchronous cancellation surface.

Geometry import has explicit source-count, per-source-byte, aggregate-byte,
and assignment-work caps. It polls the supplied `fs-exec::Cx` before source
work, per source, before and after promotion, and before ledger publication.
A pre-cancelled attempt publishes nothing. Once the atomic SQLite transaction
begins, the bounded ledger calls finish or roll back; cancellation does not
leave a partial successful operation.

Solve cancellation is not implemented by this checkpoint. It must use the
`fs-session` request -> drain -> finalize protocol, checkpoint on cancellation,
and prove resume equivalence before the solve verb stops returning
`cli-stage-unavailable`.

## Unsafe boundary and features

No unsafe code. No feature flags. Runtime dependencies remain Franken-only.

## Conformance tests

`tests/cli.rs` covers the grammar and all four v0 verbs, stable exit classes,
strict validation success, structural findings with fixes, noncanonical input
refusal, JSON escaping/line discipline, and the exact producer-Bead refusal
for each not-yet-integrated stage.

`tests/import.rs` supplies a closed reference tetrahedron and covers G0 retained
lineage, repair of deterministic duplicate/degenerate STL facets, G3 changed
source identity, open-mesh promotion refusal, mis-scaled unit refusal, dangling
assignment refusal, G4 pre-cancellation with zero publication, and G5
content-identity equivalence across independent ledgers. Every recorded case
runs the ledger linter.

## No-claim boundaries

- `validate` proves only canonical structural and dimensional admissibility.
  It does not prove referenced artifacts or material cards exist, a requested
  capability is installed, the project is solvable, or any physical model is
  valid.
- Geometry import binds exact raw bytes, lower-layer receipts, one promoted
  finite tessellation, assignment reports, and their lineage. The legacy FNV
  hook and caller path/label do not authenticate custody, physical/CAD
  sameness, continuum coverage, units, or topology beyond the retained
  lower-layer claims.
- Importing STL/OBJ/PLY is not a faceted-STEP support claim. Named face groups
  are caller-supplied labels on the promoted soup, not independently certified
  CAD product-structure identity.
- The presence of solve/report/package in help and parsing is not an
  implementation claim. Until their named authorities land, execution fails
  before side effects.
- No solve cancellation, checkpoint/restart, run identity, report rendering,
  evidence packaging, continuum enclosure, or end-to-end product determinism
  claim is made by this checkpoint.
