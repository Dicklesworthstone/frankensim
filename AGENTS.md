# AGENTS.md - FrankenSim

> Guidelines for AI coding agents working in this Rust simulation, geometry,
> optimization, and rendering codebase.

---

## RULE 0 - THE FUNDAMENTAL OVERRIDE PREROGATIVE

If I tell you to do something, even if it goes against what follows below, YOU
MUST LISTEN TO ME. I AM IN CHARGE, NOT YOU.

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

This repository is currently plan-first. The authoritative design document is:

- `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md`

Read it before broad design work. It defines the Decalogue, architecture,
roadmap, crate atlas, Gauntlet verification program, performance targets, and
flagship pipelines. This `AGENTS.md` is the operating contract for agents; the
plan is the technical constitution.

If a `README.md`, crate `CONTRACT.md`, `TESTING.md`, or `.beads/` directory is
added later, read the relevant files before making nontrivial edits.

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

Expected baseline when the workspace exists:

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

## Planned Workspace Crates

Appendix A of the plan is the authoritative crate atlas. The rough intended
shape is:

### L0 / L1

- `fs-substrate` - architecture detection, topology fingerprints, dispatch
- `fs-simd` - portable SIMD, NEON, AVX-512, exploratory SME2 tiers
- `fs-alloc` - scope arenas, huge pages, pools
- `fs-exec` - asupersync-backed two-lane executor and deterministic reductions
- `fs-la` - GEMM, small dense batches, TSQR, LOBPCG/Lanczos, mixed precision
- `fs-sparse` - BSR, SELL-C-sigma, SpMV/SpMM, AMG
- `fs-fft` - Stockham FFTs, real transforms, 3D pencils
- `fs-ivl` - interval, affine, Taylor models, exact predicates
- `fs-cheb` - Chebyshev/Fourier function objects and spectral eigenproblems
- `fs-ad` - duals, FrankenTorch tape bridge, adjoint support
- `fs-rand` - Philox, Sobol/Owen, distributions
- `fs-ga` - projective and conformal geometric algebra

### L2 / L3

- `fs-geom` - Regions, charts, Rep Router, sheaf certificates
- `fs-rep-sdf`, `fs-rep-frep`, `fs-rep-nurbs`, `fs-rep-mesh`,
  `fs-rep-voxel`, `fs-rep-neural`
- `fs-xform` - FFD, RBF, level-set velocity, manifold harmonics, SIMP fields
- `fs-mesh` - dual contouring, Delaunay/Ruppert, anisotropic remeshing
- `fs-feec` - FEEC families, exact sequences, DWR
- `fs-cutfem` - ghost penalty and cut quadrature
- `fs-iga` - spline elements and shells
- `fs-solid`, `fs-lbm`, `fs-bem`, `fs-fmm`, `fs-vpm`
- `fs-time`, `fs-couple`, `fs-adjoint`, `fs-uq`

### L4 / L5 / L6

- `fs-eproc` - e-processes, e-BH, conformal e-prediction
- `fs-opt` - L-BFGS, TR-Newton-Krylov, CMA-ES, BO, NSGA, PDHG
- `fs-topo` - SIMP, level-set, ground structure, homogenization, persistence
- `fs-sos` - moment/SOS, Burer-Monteiro SDP, Lyapunov certificates
- `fs-surrogate` - neural operators, POD-DEIM, Koopman, certify-or-escalate
- `fs-render`, `fs-viz`, `fs-img`
- `fs-ir` - FrankenScript, catalog, admission
- `fs-session` - capabilities, governor, idempotency
- `fs-ledger` - schema, content addressing, time travel, explain
- `fs-plan` - cost/error models, budget allocator, tropical analytics
- `fs-report` - notebooks and reports via FrankenPandas

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
Prefer RCH for compute-heavy Cargo work in shared-agent environments.

Typical lanes once the workspace exists:

```bash
rch exec -- env CARGO_TARGET_DIR="${TMPDIR:-/tmp}/rch_target_frankensim_check" cargo check --all-targets
rch exec -- env CARGO_TARGET_DIR="${TMPDIR:-/tmp}/rch_target_frankensim_clippy" cargo clippy --all-targets -- -D warnings
rch exec -- env CARGO_TARGET_DIR="${TMPDIR:-/tmp}/rch_target_frankensim_fmt" cargo fmt --check
```

If the workspace is still plan-only, do not invent checks. State that no Cargo
workspace exists yet and validate markdown or file presence only.

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
3. If present, read `README.md`, crate `CONTRACT.md`, and the relevant Beads
   issue.
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

Use RCH for CPU-heavy Cargo work when available:

```bash
rch exec -- env CARGO_TARGET_DIR="${TMPDIR:-/tmp}/rch_target_frankensim_check" cargo check --all-targets
rch exec -- env CARGO_TARGET_DIR="${TMPDIR:-/tmp}/rch_target_frankensim_test" cargo test --all-targets
```

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

At the time this file was written, `/Users/jemanuel/projects/frankensim`
contained only:

- `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md`
- `AGENTS.md`

There was no git repository and no Cargo workspace yet. When those are added,
update this section or replace it with the actual workspace map.

---

## Building P0 First

The coda of the plan says: build P0. Agents should interpret that literally.

Initial implementation should prioritize:

- `fs-substrate`
- `fs-exec`
- `fs-alloc`
- `fs-la`
- `fs-sparse`
- `fs-fft`
- `fs-ivl`
- `fs-rand`
- ledger v0

P0 exit criteria from the plan:

- G0 and G4 green.
- GEMM, SpMV, and FFT within the stated target bands on both reference ISA
  families.
- deterministic mode bit-stable.

Do not start with moonshot features. Keep `[M]` work behind flags and contracts
until the solid spine exists.

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
