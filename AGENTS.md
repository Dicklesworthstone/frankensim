# AGENTS.md - FrankenSim

> Guidelines for AI coding agents working in this Rust simulation, geometry,
> optimization, and rendering codebase.

---

## RULE 0 - THE FUNDAMENTAL OVERRIDE PREROGATIVE

If I tell you to do something, even if it goes against what follows below, YOU
MUST LISTEN TO ME. I AM IN CHARGE, NOT YOU.

---

## RULE 0.1 - VALUE DELIVERY OVER PROCESS: NO PROCESS PORN

FrankenSim is missing large amounts of real product functionality. Agent time,
user time, tokens, review attention, and repository complexity are scarce. The
default job is to implement the requested capability or fix the concrete defect,
not to elaborate the machinery surrounding the work.

**Process is never the product unless the user explicitly asks for process
work.** Beads, plans, Agent Mail, audits, manifests, provenance, logging, CI,
test harnesses, dashboards, status reports, and agent coordination exist only to
support delivery. They must not become a self-perpetuating substitute for
delivery.

### The value test

Before doing non-product work, answer all three questions:

1. What concrete user-visible capability, correctness defect, or immediate
   implementation blocker does this work address?
2. Is this the smallest direct action that addresses it?
3. Will its likely value exceed its implementation, maintenance, review, and
   delay cost right now?

If any answer is unclear, **do not do the work**. Return to implementing the
requested functionality. Speculative future usefulness, elegance, completeness,
or “more confidence” is not enough.

### Hard scope limits

- A review request authorizes reading the relevant code, reproducing concrete
  defects, making the smallest sound fixes, and adding focused regressions. It
  does **not** authorize redesigning adjacent infrastructure, inventing a new
  analyzer, exhaustively hardening hypothetical cases, or chasing unrelated
  pre-existing failures.
- Never turn a small review or repair into hundreds or thousands of lines of CI,
  harness, analyzer, schema, logging, provenance, or planning code without the
  user's explicit approval for that expansion.
- If incidental support work is becoming comparable to or larger than the
  requested implementation, stop before expanding it. Report the concrete
  blocker and ask whether the user wants the support work. Do not rationalize
  continuing because time has already been spent.
- Do not build a validator for a validator, a harness for a harness, or an
  analyzer whose principal purpose is proving internal process artifacts unless
  that exact system is the requested deliverable.
- Do not chase a failure already present on `HEAD` unless it blocks the requested
  deliverable and the user authorizes broadening the scope. Record it briefly and
  continue or stop at the real boundary.
- Do not repeatedly re-audit, re-plan, re-hash, re-seal, poll, or add reviewers
  after the requested behavior has focused evidence. One coherent implementation
  plus proportionate verification is better than layers of ceremonial assurance.
- Do not spawn an agent swarm for a narrow task. Delegate only concrete,
  independent implementation or bounded verification that materially shortens
  the path to the requested result. Never create recursive review loops or wait
  repeatedly for agents with no deliverable in hand.

### Implementation must dominate

Unless the user explicitly requests planning, governance, CI, or tooling:

- Spend the dominant share of effort on working product code and direct tests of
  that product code.
- Prefer an existing test seam over creating a new framework. Tests should prove
  changed behavior, important boundaries, and named claims; exhaustive testing
  of internal ceremony is negative value.
- Keep logging actionable and proportional. “Great logging” means enough context
  to diagnose a real failure, not recording every intermediate state or building
  a second product around evidence collection.
- Treat Beads as concise execution bookkeeping. Once a Bead is sufficiently
  clear to implement safely, implement it. Do not spend hours optimizing issue
  prose, graph metrics, dependencies, or acceptance wording while its actual
  functionality remains absent.
- Run the narrowest relevant checks while iterating. Run broad DSR/repository
  gates only when there is a coherent implementation ready for that proof or the
  user explicitly asks for them.
- When choosing between missing P0 functionality and optional meta-infrastructure,
  implement the P0 functionality. Tooling wins only when it is a demonstrated
  blocker to that implementation.

### Mandatory checkpoint and stop rule

At the first sign of scope expansion, pause and state in plain language:

- the user-facing outcome being delivered;
- the files and approximate size of the proposed expansion;
- why the expansion is strictly necessary; and
- the smaller alternatives considered.

If that explanation cannot establish immediate net value in a few sentences,
**do not proceed**. Ask the user before expanding. A “fresh-eyes” instruction is
not permission for an unbounded expedition.

When work has drifted into process porn, stop immediately. Do not add more tests,
proof layers, or cleanup to justify the sunk cost. Freeze the tree, disclose the
exact state candidly, and wait for direction.

### Concrete anti-pattern

Adding roughly 1,500 lines to an internal Beads/CI harness during a fresh-eyes
review, then spending hours extending static analysis and evidence machinery
while user-facing FrankenSim functionality remains missing, is a canonical
failure—not thoroughness. Never repeat it.

---

## RULE NUMBER 1: NO FILE DELETION

**YOU ARE NEVER ALLOWED TO DELETE A FILE WITHOUT EXPRESS PERMISSION.** Even a
new file that you yourself created, such as a test file. You must always ask and
receive clear, written permission before deleting a file or folder of any kind.

---

## Irreversible Git & Filesystem Actions - DO NOT EVER BREAK GLASS

1. **Absolutely forbidden commands:** `git reset --hard`, `git clean -fd`,
   `rm -rf`, or any command that can delete or overwrite code/data must never be
   run unless the user explicitly provides the exact command and states, in the
   same message, that they understand and want the irreversible consequences.
2. **No guessing:** If there is any uncertainty about what a command might
   delete or overwrite, stop immediately and ask the user for specific approval.
3. **Safer alternatives first:** When cleanup or rollbacks are needed, use
   non-destructive inspection first: `git status`, `git diff`, backups, or
   explicit hand-written patches.
4. **Mandatory explicit plan:** Even after explicit user authorization, restate
   the command verbatim, list exactly what will be affected, and wait for
   confirmation that your understanding is correct.
5. **Document the confirmation:** When running any approved destructive command,
   record the user text that authorized it, the command actually run, and the
   execution time in your final response.

---

## Git Branch: ONLY Use `main`, NEVER `master`

When this directory is a git repository, the default branch is `main`.

- All work happens on `main`.
- Never create, switch to, or push feature branches unless the user explicitly
  overrides this file.
- Never reference `master` in code or docs. If you see it, treat it as a bug.
- If the remote also needs a legacy `master` ref, synchronize it from `main`
  only when the user or project automation asks for that exact operation.

---

## RULE 2: NO GIT BRANCHES. NO GIT WORKTREES. EVER.

`main` is the one and only branch. There is no "temporary" branch, no per-agent
branch, no per-task branch, and no scratch worktree.

### FORBIDDEN

- `git branch <anything-other-than-main>`
- `git checkout -b <foo>` or `git switch -c <foo>`
- `git worktree add ...`
- Pushing non-main refs to `origin`
- Creating pull requests or draft PRs from feature branches
- Working in scratch clones at paths like `/tmp/frankensim-*`,
  `/data/projects/frankensim-*`, or `~/projects/frankensim-*` to isolate work
- Using any tool or harness that creates branches or worktrees as a side effect

### WHAT YOU DO INSTEAD

- Commit directly to `main` when the user asks for commits and the work is ready.
- Keep unfinished work in the working tree.
- Coordinate through Agent Mail reservations when multiple agents are active.
- Use Beads issue IDs and file reservations as the isolation mechanism, not git
  branches.
- If another agent changed files, do not revert or stash their work. Work with
  the current tree.

---

## Project Truth Sources

The code and its tests say what exists. Documents say what was intended. Where
they disagree, the code is right about the present and the document is either a
goal or stale; say which when you report.

- `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md` is the technical constitution: the
  Decalogue, architecture, crate atlas, Gauntlet verification program and
  performance targets. Read the relevant sections before broad design work.
- The other root `COMPREHENSIVE_*.md` files are plans for one program each
  (bridge plan, new domains, instruments, building blocks, Wright Flyer).
- `README.md` and each crate's `CONTRACT.md` state what is claimed and what is
  refused. Read the crate's contract before a nontrivial edit to that crate.
- `.beads/` is the work tracker. The current convergence program is the epic
  `frankensim-rc-1109-root-3gc12` (label `reality-check-2026-11-09`). Its
  operating rules are the next section; pick work from it with
  `bv --robot-plan --label reality-check-2026-11-09`.

This `AGENTS.md` is the operating contract for agents.

---

## Operating Rules From The 2026-11-09 Review

The owner agreed these on 2026-11-09 after a whole-project review. They bind
every agent that writes to `main`.

### 1. Compose, do not enumerate

Coupled-physics capability lands in one of three forms only:

- a **Component**: one physics owner's state, ports and step behind the shared
  contract in `fs-ports`, with the conformance tests;
- a **connection or interface law**: a junction, a transformer, a mapped field
  transfer, or a dissipative law built from a material or interface card;
- a **description plus an oracle test**: an assembly file that the generic path
  runs, with an analytic, manufactured or measured reference.

Do not add a new pairing module, and do not add a new `*_with_*` entry point
(for example `advance_liquid_regions_with_moving_gas_and_heat`). Each one is a
cell of a physics x integrator x controller matrix that has no end, and every
cell must later be migrated or removed. Between 2026-10-26 and 2026-11-09
`crates/fs-couple/src` gained 822,563 lines and lost 12,841 that way.

Until the Component contract is frozen at v1, the accepted work is the menu on
bead `frankensim-rc-1109-root-3gc12.8.3`: compile repairs, Component
implementations for existing physics owners, adapters, the kernel slices, card
data, and oracle tests for existing leaf laws.

### 2. Behaviour is computed, never put in

An application, demo or instrument gets its behaviour from conservation laws
and interface laws acting on a described assembly of parts and cards.

- No rate, amplitude, threshold, event time or loss chosen from the expected
  result. A value that was fitted to the answer is named as fitted, with its
  source, in the card that holds it.
- An app describes the device, sets inputs, steps the engine and reads probes.
  Do not add new integration code to an app crate; the app engines that
  integrate their own state today move to this form under the `.13` beads of
  the epic.
- A demo that integrates nothing says so on its page ("kinematic
  illustration"); it does not claim to use a physics crate it does not call.
- Each demo has at least one test of a quantity that was not an input.

### 3. Source that was never compiled does not land on `main`

- Before pushing a change to Rust source, compile every crate the change
  touches, on the feature sets it touches, and run the tests you added or
  changed. Put the exact commands and their results in the commit message.
- If you cannot compile (no Cargo, missing sibling checkouts, no disk, no build
  host), do not push. Keep the change in the working tree, or send it to the
  coordinator as a patch (`git format-patch` output in a mail thread) marked
  UNCOMPILED. The landing lane applies it on the current tip, compiles and
  tests it, and then pushes it or returns it with the compiler output.
- "Compile was unavailable" in a commit message on `main` describes a defect.
  It is not an accepted disclosure.
- Documentation-only changes are not covered. A data file that a test reads is
  covered by that test.

### 4. Large files: no length limit, and no write-back of a partial read

There is no size limit on any file (owner ruling, 2026-11-09). Do not split a
file for its size and do not add size checks.

Measured on 2026-11-09: 27 times, a commit replaced a `CONTRACT.md` larger than
1 MiB with only that commit's new section (`fs-couple` 18, `fs-material` 8,
`fs-lbm` 1); megabytes of contract text were absent for days. The cause is not
identified. So:

- If a tool shows you a tracked file as empty or cut short, or refuses to read
  it, do not write the file back. Add your lines with a patch.
- Before every commit, read `git diff --stat` for the files you touched. A
  deletion you did not intend means stop and restore the file from `HEAD`.
- The `fs-couple`, `fs-material` and `fs-lbm` contracts are an index
  (`CONTRACT.md`) plus `contracts/part-NN.md`. New sections go in the index,
  where its layout note says.

### 5. Deletion during migration

RULE NUMBER 1 stands without exception. A migration bead that replaces a
pairing by a composition posts its deletion list (paths and line counts) with
the equivalence evidence. The files are removed only after the owner approves
that list in writing.

### 6. Byte controls for instrument renders

Compare renders only between runs of the same binary. A change that alters
arithmetic order, such as a kernel migration, needs a recorded A/B on the same
inputs before and after, and regenerated controls committed with both commit
hashes named.

---

## FrankenSim - This Project

FrankenSim is intended to be a single memory-safe Rust continuum for:

- computational geometry
- certified representation conversion
- physics simulation
- adjoint and derivative-free optimization
- uncertainty quantification
- rendering and scientific visualization
- replayable design ledgers and agent-native orchestration

The mission is not to wrap legacy CAD/FEM/CFD/optimization tools. The mission is
to build one typed algebra where geometry, fields, operators, derivatives, error
bounds, budgets, provenance, and cancellation travel together through every
layer.

### Load-Bearing Principles

- **Pure Rust first:** runtime dependencies are limited to `std`, asupersync,
  FrankenSQLite, FrankenNumpy, FrankenTorch, FrankenScipy, FrankenPandas, and
  FrankenNetworkx.
- **Memory safety:** default code must be safe Rust. `unsafe` is allowed only in
  narrow audited leaf modules where there is no safe alternative.
- **Determinism:** deterministic mode must be bit-stable across runs and thread
  counts on the same ISA wherever the plan claims it.
- **Cancellation-correct compute:** every unit of compute runs under explicit
  asupersync scopes and checks cancellation at bounded tile boundaries.
- **Budgets first:** accuracy, time, memory, and capability budgets are explicit
  values, not comments.
- **Certificates over vibes:** error bounds, adjoints, convergence orders,
  watertightness, topology, and stochastic claims need executable evidence.
- **One data model:** complexes and cochains are the shared substrate for
  geometry and physics.
- **Provenance-complete:** artifacts are content-addressed; operations are
  replayable through the Design Ledger.
- **Agent-first:** units, seeds, budgets, versions, and capabilities are always
  explicit.

---

## Toolchain: Rust & Cargo

Use Cargo for Rust builds. Do not introduce another package manager unless the
user explicitly asks for a non-Rust subproject.

Baseline:

- Rust 2024 edition.
- Cargo workspace with flat `fs-*` crates.
- `#![deny(unsafe_code)]` at crate or module level wherever practical.
- Explicit feature flags for frontier/moonshot capabilities.
- Release profile changes must justify the performance, determinism, and
  cancellation tradeoffs.

### Dependency Policy

Production runtime dependencies must stay Franken-only:

- `std`
- `asupersync`
- `FrankenSQLite`
- `FrankenNumpy`
- `FrankenTorch`
- `FrankenScipy`
- `FrankenPandas`
- `FrankenNetworkx`

No BLAS, LAPACK, C, C++, Fortran, OpenCASCADE, gmsh, FEniCS, MFEM, OpenFOAM,
SU2, ParaView SDKs, Python runtime requirements, or FFI-backed numerical
shortcuts in production paths.

Development-only references and comparison oracles are acceptable only when they
are isolated, documented, and not in the production dependency graph.

### External-Adapter Ruling (ADPT-2026-07)

Decision record: `fs-govern::adapter_policy` (bead f85xj.11.1, ratified under
explicit owner delegation on 2026-07-23). The tension between the Franken-only
production graph and real-world CAD/CAE adoption is resolved as OFFICIAL
QUARANTINED ADAPTERS:

- FrankenSim MAY ship optional official adapter binaries wrapping foreign
  tools (OpenCASCADE/gmsh-class kernels) as separate out-of-process
  executables with separate distribution.
- Adapters are never build-time or link-time dependencies of any workspace
  crate. The production runtime dependency graph stays Franken-only and
  `check-deps` enforcement is unchanged.
- Adapter output enters the workspace only through the fs-io quarantine
  boundary: evidence-bearing receipts, Estimate-only authority, and the
  adapter identity and version recorded in the receipt. An adapter can never
  mint authority by itself.
- In-process foreign kernels (FFI plugins) remain forbidden in production
  paths. Re-opening that option requires a full backend trust-contract
  design through a new decision record, not ad-hoc FFI.

### Unsafe Code

Unsafe is a last resort and must be treated as an auditable boundary:

- Prefer safe Rust, const generics, ownership, and explicit layouts first.
- Keep unsafe leaves small, local, and behind safe facades.
- Candidate unsafe zones: SIMD microkernels, arena allocation internals, memory
  mapping, architecture dispatch, exact low-level layout handling.
- Each unsafe exception needs a documented invariant, tests, and preferably a
  ledger/contract entry once the repo has the relevant artifacts.
- Never use unsafe to paper over lifetime, aliasing, cancellation, or
  synchronization design problems.

---

## Architecture

FrankenSim is organized as layers. Respect the dependency direction: lower
layers know nothing about higher layers.

| Layer | Name | Responsibility |
| --- | --- | --- |
| L0 | SUBSTRATE | hardware topology, arenas, SIMD dispatch, two-lane execution, determinism |
| L1 | BEDROCK | dense/sparse/FFT math, certified arithmetic, AD, RNG/QMC |
| L2 | MORPH | Regions, charts, representation routing, meshing, geometry certificates |
| L3 | FLUX | FEEC/DEC physics, CutFEM, LBM, BEM/FMM, solvers, adjoints |
| L4 | ASCENT | shape/topology/global/Bayesian/multi-objective optimization |
| L5 | LUMEN | path tracing, direct chart rendering, scientific visualization |
| L6 | HELM | FrankenScript IR, sessions, capabilities, ledgers, planner, reports |

### Dependency Direction

- `fs-substrate` and `fs-exec` may depend on asupersync and low-level utility
  crates.
- BEDROCK crates can depend on SUBSTRATE but not on MORPH, FLUX, ASCENT, LUMEN,
  or HELM.
- MORPH can depend on BEDROCK and SUBSTRATE.
- FLUX can depend on MORPH, BEDROCK, and SUBSTRATE.
- ASCENT can depend on FLUX, MORPH, BEDROCK, and SUBSTRATE.
- LUMEN can depend on MORPH, FLUX field abstractions, BEDROCK, and SUBSTRATE.
- HELM can orchestrate all layers, but lower layers must not call HELM.

If you need to violate this, stop and redesign the API.

---

## Workspace Map

As of 2026-11-09:

- `crates/` holds 193 flat `fs-*` crates. Each declares its layer in
  `[package.metadata.frankensim] layer` in its `Cargo.toml` (L0 to L6, plus
  `UTIL` and `TOOL`); `cargo run -p xtask -- check-layers` enforces the
  dependency direction. Appendix A of the plan lists the crates the plan
  foresaw; the manifests are the list of what exists.
- `xtask/` holds the repository checks (`check-all`, `check-layers`,
  `check-deps`, `check-contracts`, `check-docs` and the rest).
- `apps/` holds the browser applications: `instrument-lab`,
  `instrument-explorer`, `observatory`, `pocket-watch`, `wright-flyer`.
- `data/` holds material database packs (`data/matdb`), instrument
  measurements (`data/instruments`), and reference and V&V corpora.
- `docs/` holds design and decision records.
- Sibling repositories (the constellation, for example `../asupersync`) are
  path dependencies pinned in `constellation.lock`.

The repository contains three programs that share leaf crates but do not yet
share a composition layer:

| Program | Entry point | Reaches |
| --- | --- | --- |
| Evidence and thermal design | `frankensim` (`crates/fs-cli`) | projects, studies, cooling and electrothermal commands, reports, the design ledger |
| Coupled multiphysics | library API of `fs-couple`, `fs-material` and the physics owners (`fs-mbd`, `fs-thermal`, `fs-gas`, `fs-flux`, `fs-solid`, `fs-em`, `fs-tribo`) | tests and examples only: no `frankensim` command and no assembly file runs a coupled pairing |
| Instruments and demos | `music_render` (`crates/fs-couple/src/bin/music_render.rs`) and `apps/` | offline instrument renders graded against recordings; browser demos |

The convergence program joins them: one composition kernel (Component
contract in `fs-ports`, coupled-step driver in `fs-system`), one assembly
description that lowers to it (`fs-scenario`, `.fsim`), and one run path
through the ledger (`frankensim run <assembly>`).

Every crate must eventually ship:

- `CONTRACT.md` with invariants, error models, determinism class, and no-claim
  boundaries.
- Unit tests for local behavior.
- Conformance tests for cross-crate semantics.
- Golden or ledger-backed evidence for claims that need replay.

---

## Core Invariants

### Execution and Cancellation

- Every hot kernel is a tile program under an explicit `Cx`.
- Tiles are the unit of scheduling, cancellation, determinism, and NUMA
  placement.
- Kernels must poll cancellation at tile boundaries.
- Cancellation is request -> drain -> finalize, not silent drop.
- Speculative races must cancel and fully drain losing branches.
- Resumable solvers must serialize enough state to pause, migrate, resume, and
  fork deterministically.

### Determinism

- Deterministic reductions use fixed-shape trees keyed by logical tile identity.
- RNG streams are counter-based and keyed by logical identity, not worker thread.
- Tie-breaking must be deterministic.
- Fast mode may relax determinism only when the mode is part of the ledgered
  provenance.

### Geometry

- A `Region` is abstract; concrete representations are charts.
- No chart type is privileged globally.
- Conversions carry cost, error bounds, certificate availability, and
  provenance.
- Rep Router decisions must be reproducible and explainable.
- Watertightness, manifoldness, self-intersection freedom, topology, and
  chart-consistency claims require certificates or explicit no-claim language.

### Physics

- Fields are cochains on complexes where possible.
- Preserve exact discrete sequences: grad/curl/div identities must hold by
  construction for FEEC paths.
- Prefer matrix-free operators and roofline-honest kernels.
- Solvers must expose error estimates, residuals, cancellation points, and
  adjoint hooks when the layer claims gradient support.
- Goal-oriented error estimates should target the actual objective, not generic
  field accuracy.

### Optimization

- Optimization problems are data: typed objective/constraint graphs with
  manifold variables, budgets, seeds, and provenance.
- Gradients must be verified against independent checks before they are trusted.
- Surrogates operate under certify-or-escalate policy.
- Stochastic decisions use anytime-valid confidence/e-process machinery when
  optional stopping is possible.
- Returned optima need KKT/residual/certificate context or honest no-claim
  boundaries.

### Ledger and Agent Interface

- The Five Explicits are mandatory: units, seeds, budgets, versions, and
  capabilities.
- Mutating calls need idempotency keys.
- Failed admission should return structured diagnostics and ranked fixes.
- Every artifact should be content-addressed and explainable through lineage.
- Generated catalogs must come from code/contracts and must not drift by hand.

---

## Ambition Tags

The plan uses three risk tags. Preserve them in docs, feature gates, and status
reports:

- `[S] Solid` - established mathematics and engineering.
- `[F] Frontier` - research-backed, high-upside engineering risk.
- `[M] Moonshot` - novel synthesis; must stay behind feature flags until proven.

Do not promote `[F]` or `[M]` functionality into the default path without
Gauntlet evidence and explicit no-claim boundaries.

---

## The Gauntlet - Correctness Program

The Gauntlet is the definition of "done" for technical claims.

| Tier | Meaning |
| --- | --- |
| G0 | property tests and algebraic laws |
| G1 | manufactured solutions and convergence-order verification |
| G2 | canonical benchmarks |
| G3 | metamorphic tests |
| G4 | chaos, cancellation storms, leak/deadlock checks |
| G5 | determinism audits |

When implementing a feature, name the relevant Gauntlet tier in the tests or
docs. If a feature cannot yet satisfy the intended tier, document exactly what
is proven and what is not.

### Certifying the Certifiers

Certificate machinery is security-critical for science:

- Interval/Taylor enclosures need high-precision or independent spot checks.
- e-process validity needs null simulations with adversarial stopping.
- Conformal prediction needs coverage tests under drift scenarios.
- Sheaf/watertightness certificates need adversarial seam tests.
- A false certificate is worse than an ordinary wrong answer.

---

## Performance Program

Performance claims must be roofline-aware and measurable. Do not write "fast"
unless there is a benchmark, target, machine fingerprint, and acceptance band.

Reference target families from the plan:

- Apple Silicon, especially M-series unified-memory machines.
- Many-core x86, especially high-core-count Threadripper/EPYC class systems.

Required habits:

- Track cache-line differences: 128 bytes on Apple aarch64, 64 bytes on x86-64.
- Design for bandwidth-rich and bandwidth-starved schedules.
- Keep CCD/L3 locality in mind on high-core-count x86.
- Resolve SIMD dispatch once, not in hot loops.
- Store performance evidence in ledger/artifacts when that infrastructure
  exists.
- Treat performance regressions as test failures once baselines exist.

---

## Code Editing Discipline

### No Script-Based Code Changes

Do not run broad regex or script-based code rewrites over source files. Make
code changes manually with focused patches. Use structured tools such as
`ast-grep` only when the pattern is genuinely syntactic and the diff can be
reviewed.

### No File Proliferation

Revise existing files in place unless a new file represents genuinely new
functionality or a required contract/test artifact.

Forbidden naming patterns:

- `main_v2.rs`
- `improved.rs`
- `new_version.rs`
- `final_final.rs`
- duplicate experimental copies of existing modules

### Backwards Compatibility

This project is early-stage. Prefer the correct design over compatibility shims.
Do not preserve bad APIs through wrappers unless the user explicitly asks for a
migration layer.

### Comments and Docs

- Document invariants, error models, determinism class, and no-claim boundaries.
- Avoid comments that merely restate code.
- In math-heavy code, include enough references or derivation notes for the next
  agent to verify signs, units, and assumptions.

---

## Output Style

Core library code should not print casually to stdout/stderr.

- Use structured tracing or ledger events for observability.
- CLI output, when added, must be deterministic and documented.
- Errors intended for agents should be structured and actionable.
- Diagnostics should include units, budgets, capability context, and suggested
  fixes when possible.

---

## Compiler Checks

After substantive Rust code changes, verify that the relevant checks pass.
DSR is the first choice for repo-level gates and release builds. Prefer RCH
only for narrow ad hoc Cargo probes or when DSR itself is unavailable.

Typical lanes:

```bash
dsr quality --tool frankensim
dsr build frankensim --target darwin/arm64
```

For the crates a change touches, the narrow check is
`cargo check -p <crate> --all-targets` and `cargo test -p <crate>` on the
feature sets the change touches. Operating rule 3 requires it before a push.

---

## DSR - Required CI and Release Runner

GitHub Actions is not the CI source of truth for this repository. The account is
throttled/cut off, so agents must always use DSR in preference to GitHub
Actions for repo-level verification, release builds, and fallback release work.

- Use `dsr` if it is on `PATH`; otherwise use
  `/Users/jemanuel/projects/doodlestein_self_releaser/dsr`.
- Run `dsr repos info frankensim` when checking the registry wiring.
- Run `dsr quality --tool frankensim` for the configured quality gate.
- Run `dsr quality --tool frankensim --dry-run` when you need to inspect the
  gate without executing the Cargo workload.
- Run `dsr build frankensim --target darwin/arm64` for the configured native
  release artifact lane. Use `--allow-dirty` only when the user explicitly
  wants a build from the current dirty tree.
- Run `dsr fallback frankensim --version <version>` only for intentional
  release fallback work.
- Use `dsr doctor` and `dsr health all` for DSR/host diagnostics.

Do not wait on, poll, or cite GitHub Actions as required proof unless the user
explicitly asks for that. The workflow files are retained as manual specs and
historical gate documentation, not as automatic merge/release criteria.

When reporting verification, include the exact DSR command, pass/fail status,
and any run log or artifact path DSR prints. If DSR is unavailable, report the
exact blocker and then use RCH or local Cargo only as a clearly labeled fallback.

### Cargo-Independent Constellation Preflight

Before attributing a workspace-wide Cargo failure to FrankenSim source, run:

```bash
scripts/ci/checkout_constellation.sh --verify-only
```

This validates the canonical `constellation.lock` and compares every sibling's
live Git head and cleanliness without invoking Cargo or building `xtask`. A
nonzero result is an authoritative drift/dirty-checkout refusal; report the
expected and actual heads it prints. Never move a shared sibling checkout back
to its pin to make this preflight green. Pin movement must instead follow the
constellation governance and compatibility-train process.

The preflight tells you a sibling is off-pin. It does **not** tell you whether
that is a fast-forward or a divergence, and only the second is an incident.
The standing gate `cargo run -p xtask -- check-constellation-drift` (a
check-all member since bead frankensim-es6pt) now makes that classification
continuously: per sibling it reports `on-pin`, `stale-lock` (fast-forward;
visible but not red), `drift` (dirty, retreated, or diverged; red, naming the
sibling with expected and actual heads), or NO-DATA. Before proposing any
remedy, classify: if `git merge-base --is-ancestor <pin>
HEAD` succeeds in the sibling, the checkout is strictly newer and the *lock* is
stale, so the safe direction is forward. Measured 2026-07-29: all seven
siblings were clean strict fast-forwards, 25 to 611 commits ahead. Treating
that as drift and "restoring" the pins would have discarded every one of those
commits.

### Regenerating Derived Artifacts On A Shared Tree

`cargo run -p xtask -- generate-source-manifest` with **no arguments** is safe
while other agents have uncommitted work. It is a pure function of `HEAD`:
every tracked path that is dirty and not named on the command line is read
from `HEAD` via `git show`, so another lane's edits cannot leak into an
artifact that claims to describe a commit
(`xtask/src/source_manifest.rs::capture_source_files_scoped`). Pass paths as
arguments only for the files you are about to commit yourself.

There is therefore **no "quiet tree" precondition**. Do not decline to
regenerate on the belief that it will absorb someone else's WIP; that belief
appeared in several commit messages, was copied forward as precedent, and is
false. If a derived artifact is stale, it is stale at `HEAD` and that is always
a real violation.

---

## Testing Policy

Tests scale with risk.

### Unit Tests

Every module should include focused tests for:

- happy paths
- empty/boundary/max cases
- error conditions
- unit/dimension correctness where relevant
- deterministic tie-breaking

### Property and Metamorphic Tests

Use property tests for algebraic laws and invariants:

- chart conversion round trips within certified bounds
- adjoint identity checks
- exact-sequence identities
- rigid transform invariance
- unit-rescaling invariance
- refinement monotonicity
- interval containment under equivalent rewrites

### Concurrency and Cancellation Tests

Concurrency-sensitive code needs deterministic lab-runtime or model-checking
coverage:

- no task leaks
- no arena leaks
- loser branches drained
- cancellation latency bounded at tile boundaries
- pause/resume deterministic equivalence
- panic/fault containment propagates structured errors

### Golden Evidence

Golden artifacts should be deterministic, reproducible, and tied to contracts.
Do not regenerate golden files casually. If a golden changes, explain the
semantic reason and run the relevant verifier.

---

## Documentation and Contracts

Each crate should have a `CONTRACT.md` before it becomes a dependency target for
other crates. A contract should state:

- purpose and layer
- public types and semantics
- invariants
- error model
- determinism class
- cancellation behavior
- unsafe boundary, if any
- feature flags
- conformance tests
- no-claim boundaries

The `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md` should stay high-level. Do not turn
it into a dumping ground for implementation notes that belong in crate docs,
contracts, or issue records.

---

## Agent Workflow

### Start of Work

1. Read this file.
2. Read `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md` sections relevant to the task.
3. Read the crate `CONTRACT.md` and the relevant Beads issue.
4. Inspect the tree before editing.
5. Reserve files through Agent Mail if multiple agents are active and the tools
   are available.

### While Working

- Keep changes tightly scoped.
- Do not disturb unrelated files.
- Do not revert changes you did not make.
- Keep technical claims attached to tests, contracts, or no-claim language.
- Prefer small, reviewable patches.

### End of Work

Before ending a session:

1. Run applicable format/check/test lanes.
2. Report exactly what passed and what did not run.
3. If Beads is in use, update issue status and run `br sync --flush-only`.
4. Release Agent Mail file reservations if you made any.
5. Leave clear handoff notes for any blocker.

---

## MCP Agent Mail - Multi-Agent Coordination

When Agent Mail tools are available, use them for multi-agent coordination and
file reservations.

Typical same-repository flow:

1. Register identity for this project path.
2. Reserve files before editing.
3. Use the Beads issue ID or task name as the thread ID.
4. Send start/progress/completion messages for shared work.
5. Release reservations when finished.

Reservations are advisory, but in this project they are the isolation mechanism.
They replace branches and worktrees.

---

## Beads (`br`) - Issue Tracking

If `.beads/` exists, use Beads for task state and dependency tracking.

Useful commands:

```bash
br ready
br list --status=open
br show <id>
br update <id> --status=in_progress
br close <id> --reason "Completed"
br sync --flush-only
```

Conventions:

- Use Beads IDs in Agent Mail thread IDs and commit messages.
- Do not run bare interactive tools in automated sessions if robot/non-TUI
  modes exist.
- `br sync --flush-only` does not commit; stage and commit intentionally only
  when the user asks for commits or the workflow requires it.

---

## `bv` - Graph-Aware Triage

If `bv` is available and `.beads/` exists, use robot modes only. Bare `bv`
launches an interactive TUI and can block the session.

```bash
bv --robot-triage
bv --robot-next
bv --robot-plan
bv --robot-insights
```

Use `bv` for work selection and dependency insight. Use Agent Mail for
coordination and file reservations.

---

## `ubs` - Bug Scanner

Before committing code, run `ubs` on changed files when available:

```bash
ubs <changed-files>
ubs $(git diff --name-only --cached)
```

Fix true positives at the root cause and rerun on the affected files.

---

## RCH - Remote Compilation Helper

Use DSR first for repo-level gates. Use RCH for CPU-heavy Cargo probes when DSR
does not cover the needed check or when you are intentionally running a narrow
diagnostic:

```bash
rch exec -- env CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_check" cargo check --all-targets
rch exec -- env CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" cargo test --all-targets
```

**ALWAYS use `${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}` as the base — never `${TMPDIR:-/tmp}`
or a bare `/tmp` — when you invent a new per-task target dir name.** These isolated
target dirs are 5-15G *each* and are never cleaned up by the task that made them.
On the Mac, `$TMPDIR` and `/tmp` both live on the internal Data volume: 600G+ of
`rch_target_frankensim_*` accumulated there and wedged the machine at 100% full
(2026-07-11). Since 2026-10-11 `RCH_TARGET_BASE` is unset on the Mac as well
(the external NVMe is retired for builds), so the base falls back to `$TMPDIR`
on the internal disk; on Linux workers `$TMPDIR`=/data/tmp.

Prefer reusing one of the two names above. If you truly need a distinct target dir
for an isolated task, still base it on `$RCH_TARGET_BASE` and delete it when done.

**Disk preflight (machine law, bead frankensim-z4igc):** before any heavy
build/test lane on this Mac, run `sbh check --need 20G` (adjust `--need` for
the lane's expected footprint). A nonzero exit means the Data volume is under
pressure: STOP and reclaim first. A full disk silently corrupts shared SQLite
writers (`br`, cass, agent mail) mid-commit. The `sbh` daemon already scans
and auto-cleans `/tmp`; this preflight is the agent-side gate for everything
else.

Do not rely on local fallback for heavy builds in a shared-agent environment
unless the user explicitly authorizes it.

Quick diagnostics:

```bash
rch doctor
rch status
rch queue
```

---

## Search Tools

- Use `rg` for fast targeted text search.
- Use `rg --files` to inspect file sets.
- Use `ast-grep` when syntax matters.
- Use project-specific AI search tools only for exploratory architecture
  questions, not for exact symbol searches.

Do not use broad scripted rewrites where a hand patch is safer.

---

## Current Repository State

As of 2026-11-09 (replace the dated facts when they change):

- The workspace is described in "Workspace Map" above.
- The workspace has not compiled as a whole with default features since
  2026-10-23 (bead `frankensim-g7u99`). Single crates do compile; check the
  ones you touch.
- `suite-receipt.json` is from 2026-08-01 and is not green. No retained
  receipt describes the current tree. Do not cite it as current evidence.
- `README.md` was last revised on 2026-10-16 and is behind the code.

---

## What To Build First

The plan's P0 crates (substrate, execution, dense and sparse algebra, FFT,
intervals, random numbers, ledger) exist. The work that matters most now is
convergence, in this order:

1. Keep `main` compiling: the compile census and landing rule (beads `.1.x` of
   the epic).
2. The composition kernel: the Component contract, connections, the system
   state store and the coupled-step driver (beads `.5.x`), proved on five
   slices of existing behaviour (beads `.6.x`).
3. The description path: one assembly model, its lowering to the kernel, and
   `frankensim run <assembly>` (beads `.11.x`).
4. Migration of the existing pairings, instruments and demos onto the kernel
   (beads `.9.x`, `.12.x`, `.13.x`), with oracle tests (beads `.14.x`).

Do not start with moonshot features. Keep `[M]` work behind flags and contracts.

<!-- bv-agent-instructions-v2 -->

---

## Beads Workflow Integration

This project uses [beads_rust](https://github.com/Dicklesworthstone/beads_rust) (`br`) for issue tracking and [beads_viewer](https://github.com/Dicklesworthstone/beads_viewer) (`bv`) for graph-aware triage. Issues are stored in `.beads/` and tracked in git.

### Using bv as an AI sidecar

bv is a graph-aware triage engine for Beads projects (.beads/beads.jsonl). Instead of parsing JSONL or hallucinating graph traversal, use robot flags for deterministic, dependency-aware outputs with precomputed metrics (PageRank, betweenness, critical path, cycles, HITS, eigenvector, k-core).

**Scope boundary:** bv handles *what to work on* (triage, priority, planning). `br` handles creating, modifying, and closing beads.

**CRITICAL: Use ONLY --robot-* flags. Bare bv launches an interactive TUI that blocks your session.**

#### The Workflow: Start With Triage

**`bv --robot-triage` is your single entry point.** It returns everything you need in one call:
- `quick_ref`: at-a-glance counts + top 3 picks
- `recommendations`: ranked actionable items with scores, reasons, unblock info
- `quick_wins`: low-effort high-impact items
- `blockers_to_clear`: items that unblock the most downstream work
- `project_health`: status/type/priority distributions, graph metrics
- `commands`: copy-paste shell commands for next steps

```bash
bv --robot-triage        # THE MEGA-COMMAND: start here
bv --robot-next          # Minimal: just the single top pick + claim command

# Token-optimized output (TOON) for lower LLM context usage:
bv --robot-triage --format toon
```

#### Other bv Commands

| Command | Returns |
|---------|---------|
| `--robot-plan` | Parallel execution tracks with unblocks lists |
| `--robot-priority` | Priority misalignment detection with confidence |
| `--robot-insights` | Full metrics: PageRank, betweenness, HITS, eigenvector, critical path, cycles, k-core |
| `--robot-alerts` | Stale issues, blocking cascades, priority mismatches |
| `--robot-suggest` | Hygiene: duplicates, missing deps, label suggestions, cycle breaks |
| `--robot-diff --diff-since <ref>` | Changes since ref: new/closed/modified issues |
| `--robot-graph [--graph-format=json\|dot\|mermaid]` | Dependency graph export |

#### Scoping & Filtering

```bash
bv --robot-plan --label backend              # Scope to label's subgraph
bv --robot-insights --as-of HEAD~30          # Historical point-in-time
bv --recipe actionable --robot-plan          # Pre-filter: ready to work (no blockers)
bv --recipe high-impact --robot-triage       # Pre-filter: top PageRank scores
```

### br Commands for Issue Management

```bash
br ready              # Show issues ready to work (no blockers)
br list --status=open # All open issues
br show <id>          # Full issue details with dependencies
br create --title="..." --type=task --priority=2
br update <id> --status=in_progress
br close <id> --reason="Completed"
br close <id1> <id2>  # Close multiple issues at once
br sync --flush-only  # Export DB to JSONL
```

### Workflow Pattern

1. **Triage**: Run `bv --robot-triage` to find the highest-impact actionable work
2. **Claim**: Use `br update <id> --status=in_progress`
3. **Work**: Implement the task
4. **Complete**: Use `br close <id>`
5. **Sync**: Always run `br sync --flush-only` at session end

### Key Concepts

- **Dependencies**: Issues can block other issues. `br ready` shows only unblocked work.
- **Priority**: P0=critical, P1=high, P2=medium, P3=low, P4=backlog (use numbers 0-4, not words)
- **Types**: task, bug, feature, epic, chore, docs, question
- **Blocking**: `br dep add <issue> <depends-on>` to add dependencies

### Session Protocol

```bash
git status              # Check what changed
git add <files>         # Stage code changes
br sync --flush-only    # Export beads changes to JSONL
git commit -m "..."     # Commit everything
git push                # Push to remote
```

<!-- end-bv-agent-instructions -->

For any web requests you must make with curl or otherwise, identify yourself honestly: use a User-Agent that names this project, e.g. "frankensim-agent (+https://github.com/Dicklesworthstone/frankensim)". Never borrow another product's or company's crawler or client User-Agent, and never spoof a client to get past a block.
