<div align="center">
  <img src="frankensim_illustration.webp" alt="FrankenSim illustration" width="720">
</div>

# FrankenSim

<div align="center">

[![Status](https://img.shields.io/badge/status-active%20Rust%20workspace-2ea44f)](#implemented-workspace)
[![Rust](https://img.shields.io/badge/rust-nightly%202024-b7410e)](rust-toolchain.toml)
[![Crates](https://img.shields.io/badge/workspace-29%20library%20crates-0969da)](#implemented-workspace)
[![Contracts](https://img.shields.io/badge/contracts-28%20of%2029%20crates-8250df)](#contracts-and-verification)
[![Tests](https://img.shields.io/badge/tests-42%20crate%20test%20files-1f883d)](#contracts-and-verification)
[![License](https://img.shields.io/badge/license-MIT%20%2B%20AI%20rider-yellow)](LICENSE)

</div>

FrankenSim is a working Rust workspace for deterministic geometry, certified numerics, meshing, execution, evidence, and design-ledger infrastructure for simulation and design optimization.

The tree contains a real Cargo workspace with 29 library crates, repository policy tooling, conformance contracts, integration tests, and working implementations across substrate/runtime, numerical kernels, geometry representations, meshing, transformation, imaging, evidence, and ledger layers.

There is not yet a packaged end-user simulation application or crates.io release. Today, FrankenSim is usable as a source workspace and library substrate.

## TL;DR

**The problem:** simulation systems often split physical units, numerical error, runtime behavior, geometry validity, evidence, and reproducibility across separate tools. That makes it too easy for an optimization run to produce an answer without a durable explanation of which assumptions, approximations, kernels, and machine conditions made the answer valid.

**The solution:** FrankenSim builds those concerns into the workspace architecture. Units are represented explicitly, kernels have deterministic contracts, geometry conversions carry evidence, runtime behavior is structured around cancellable work contexts, and the ledger records artifacts, operations, events, roofline measurements, and time-travelable design state.

### What Exists Now

| Area | Current implementation |
|------|------------------------|
| Workspace | Rust 2024 nightly Cargo workspace with 29 `fs-*` library crates plus `xtask` |
| Contracts | 28 of 29 library crates have `CONTRACT.md`; `fs-img` is the current contract gap |
| Runtime substrate | Capability probing, SIMD facades, aligned arenas, two-lane execution, cancellation contexts, tile pools, tuner and race scaffolding |
| Numerics | Deterministic elementary math, dense/sparse linear algebra, FFT/DCT, interval/affine/Taylor arithmetic, Chebyshev collocation, random/QMC streams, AD/adjoint infrastructure, e-process inference |
| Geometry | Region/chart abstraction, SDF, mesh and F-rep charts, representation conversion hooks, transformations, tet meshing, remeshing, quality audits |
| Evidence and ledger | Composable `Evidence<T>`/`Certified<T>`, model cards, bracketing, FrankenSQLite-backed design ledger, artifact hashes, event streams, tune cache, roofline recording |
| Policy tooling | `xtask` checks for layer direction, Franken-only runtime dependencies, contracts, unsafe capsules, and constellation lock verification |
| Tests | 42 crate-level conformance and integration test files exercising the implemented contracts |

### What You Can Use Today

| Task | Implemented path |
|------|------------------|
| Build source crates and run conformance tests | `cargo test --workspace` |
| Enforce repo architecture rules | `cargo run -p xtask -- check-all` |
| Assemble deterministic sparse matrices and multiply them | `fs_sparse::{Coo, Csr, Bsr, Sell}` |
| Run Chebyshev collocation and Orr-Sommerfeld stability probes | `fs_cheb::Cheb1` and `fs_cheb::orr_sommerfeld` |
| Attach numerical/model/statistical evidence to values | `fs_evidence::{Evidence, Certified}` |
| Use physical quantities with dimensional checks | `fs_qty` compile-time quantities and runtime `QtyAny` |
| Work with SDF, mesh, and F-rep geometry charts | `fs_geom`, `fs_rep_sdf`, `fs_rep_mesh`, `fs_rep_frep` |
| Apply design parameterizations and detect foldover | `fs_xform` FFD, RBF, velocity band, density, and composition types |
| Encode image artifacts and apply deterministic film/denoising transforms | `fs_img` PNG/OpenEXR subset plumbing, film transforms, and bias-labeled denoising |
| Record artifacts, operations, metrics, events, and tune rows | `fs_ledger` on FrankenSQLite |
| Probe machine axes and record roofline measurements | `fs_roofline` plus `fs_substrate` |

## Quick Start

FrankenSim currently builds from source. The workspace expects sibling Franken projects for path dependencies, especially `~/projects/asupersync` and `~/projects/frankensqlite`.

```bash
git clone https://github.com/Dicklesworthstone/frankensim.git
cd frankensim

# Use the pinned nightly toolchain from rust-toolchain.toml.
rustup toolchain install nightly
rustup component add rustfmt clippy --toolchain nightly

# Build and test the workspace.
cargo test --workspace
```

Run the repository policy checks:

```bash
cargo run -p xtask -- check-all
```

Run the project validation lane through DSR, which is the preferred verification path for this repo:

```bash
DSR_BIN="$(command -v dsr || printf '%s' /Users/jemanuel/projects/doodlestein_self_releaser/dsr)"
"$DSR_BIN" repos info frankensim
"$DSR_BIN" quality --tool frankensim
```

GitHub Actions is not the source of truth for this project. Use DSR for automation and verification.

## Command Reference

FrankenSim does not yet ship an end-user simulator CLI. The commands below are the implemented repository and validation entrypoints.

| Command | Purpose |
|---------|---------|
| `cargo test --workspace` | Build the workspace and run crate-level unit, integration, and conformance tests |
| `cargo fmt --check` | Check formatting under the pinned nightly toolchain |
| `cargo clippy --workspace --all-targets -- -D warnings` | Run the strict local lint lane |
| `cargo run -p xtask -- check-layers` | Enforce the layer dependency direction |
| `cargo run -p xtask -- check-deps` | Enforce the Franken-only runtime dependency policy |
| `cargo run -p xtask -- check-contracts` | Verify `fs-*` crate contracts; currently expected to flag `fs-img` until its contract lands |
| `cargo run -p xtask -- check-unsafe` | Check unsafe code against registered capsules |
| `cargo run -p xtask -- check-constellation` | Verify sibling Franken repository pins |
| `cargo run -p xtask -- check-all` | Run the implemented `xtask` policy checks together |
| `dsr quality --tool frankensim` | Run the configured repo-level DSR quality gate when DSR is available |
| `dsr build frankensim --target darwin/arm64` | Run the configured native DSR build lane |

## Implemented Workspace

FrankenSim is organized as layered Rust crates. The names below are crates present under `crates/`.

### Units, Observability, and Evidence

| Crate | What is implemented |
|-------|---------------------|
| `fs-qty` | Compile-time dimensional quantities, runtime `QtyAny`, SI parsing, dynamic dimension checks, and serialization support |
| `fs-obs` | Structured event schema, event kinds, emitters, content hashing, and shared observability types for kernels and ledgers |
| `fs-evidence` | `Evidence<T>` and `Certified<T>` wrappers, numerical/statistical/model evidence, model cards, discrepancy records, bracketing, composition, and conservative assessment |

### Substrate and Execution

| Crate | What is implemented |
|-------|---------------------|
| `fs-substrate` | Machine capability probes, fingerprints, topology maps, dispatch tiers, bandwidth measurement, Morton/tile identifiers, tiled field helpers, and CCD topology primitives |
| `fs-simd` | Scalar baseline operations plus registered NEON/x86 unsafe capsules behind safe facades and alignment contracts |
| `fs-alloc` | 128-byte aligned allocation helpers, cache padding, scoped arenas, arena pools, site stats, hugepage policy/outcome types, and sharded pools |
| `fs-exec` | `Cx` execution context, stream keys, cancellation gates, tile kernels, tile pools, reductions, latency lane hooks, work distribution, race/tuner scaffolding, solver state, and kill registry |

### Numerical Kernels

| Crate | What is implemented |
|-------|---------------------|
| `fs-math` | Deterministic `exp`, `expm1`, `ln`, `sin`, `cos`, `tanh`, `sqrt`, complex numbers, error-free transforms, double-double helpers, and conformance boundaries |
| `fs-la` | Dense GEMM, batched small dense operations, factorizations, mixed-precision helpers, real and complex eigensolver scaffolding, randomized NLA range/trace estimators, and conformance tests |
| `fs-sparse` | COO, CSR, BSR, and SELL formats, deterministic assembly, transpose/symmetrize, SpMV/SpMM, Chebyshev, ILU0, PCG, and smoothed-aggregation AMG pieces |
| `fs-fft` | Complex FFT, real FFT, DCT-II/DCT-III, Stockham-style transform structure, and transform conformance tests |
| `fs-ivl` | Outward-rounded intervals, affine arithmetic, first-order Taylor models, Newton/Krawczyk root helpers, expansions, and exact predicates (`orient`, `incircle`, `insphere`) |
| `fs-cheb` | Adaptive 1D Chebyshev functions, Lobatto points, differentiation matrices, Dirichlet Laplacian eigen checks, and Orr-Sommerfeld growth-rate utilities |
| `fs-rand` | Philox counter-based streams keyed by logical identity, distributions, Sobol sequences, Owen/QMC-style components, and lattice helpers |
| `fs-ad` | `Real` trait, forward-mode dual numbers, gradient checks, implicit-function adjoint hooks, and checkpointed/full adjoint scaffolding |
| `fs-eproc` | Betting e-processes, pairwise races, Gaussian mixture confidence sequences, and e-Benjamini-Hochberg support |

### Geometry, Representations, and Meshing

| Crate | What is implemented |
|-------|---------------------|
| `fs-geom` | Points, vectors, AABBs, Betti-bound summaries, `Chart` and `Region` traits, conversion records, sampled SDF fixtures, and representation-router scaffolding |
| `fs-rep-sdf` | Dense tiled SDF grids, FrankenVDB-style sparse tree structure, adaptive octree SDF, narrow-band helpers, and SDF chart implementations |
| `fs-rep-mesh` | Half-edge surfaces, oriented tet complexes, soup repair, generalized winding-number support, point-triangle distance, mesh charts, shapes, dual contouring, mesh-to-SDF, and quality assessment |
| `fs-rep-frep` | CSG/F-rep builders, differentiable Boolean operations, chart implementation, finite-difference parameter derivatives, and interval/Lipschitz/gradient-style evaluation hooks |
| `fs-xform` | Free-form deformation lattices, RBF morphs, level-set velocity bands, SIMP-style density fields, composed parameterizations, Jacobian actions, and foldover detection |
| `fs-mesh` | Incremental Delaunay/tetrahedralization scaffolding, exact audits, quality refinement, ghost/refinement helpers, metric fields, remeshing, and execution-aware meshing hooks |

### Ledger, Planning, Roofline, and Vertical Skeleton

| Crate | What is implemented |
|-------|---------------------|
| `fs-ledger` | Design Ledger schema v2 on FrankenSQLite, content-addressed artifacts, operations, event streams, lineage, metrics, tune cache, extension rows, integrity/lint checks, and time-travel queries |
| `fs-ir` | FrankenScript typed AST with spans, isomorphic s-expression and JSON syntaxes, shape comparison, study recognition, lowering, and structured IR errors |
| `fs-plan` | Cost model types, error and time ledgers, plan-cost oracle, and cost-model construction from tune records |
| `fs-roofline` | Machine-axis probing, kernel specs, roofline kernel registry, measurement harness, section-target constants, ledger recording, and staleness checks |
| `fs-img` | LUMEN image pipeline crate with deterministic PNG/OpenEXR subset plumbing, film/display transforms, bias-labeled denoising, and conformance tests |
| `fs-vskeleton` | Photovoltaic vertical skeleton tying SDF/PDE/objective/adjoint/optimization/ledger concepts into a narrow demonstrator path |
| `fs-opt` | ASCENT optimization problem IR crate with typed objective/constraint graphs, dimensional validation, differentiability-class routing, canonical serialization, manifold metadata, and a conformance contract |

## Examples

These examples are library-level examples. FrankenSim does not yet expose a stable end-user CLI.

### Deterministic Sparse Assembly

```rust
use fs_sparse::{Coo, Csr};

let mut coo = Coo::new(3, 3);
coo.push(0, 0, 2.0);
coo.push(0, 1, -1.0);
coo.push(1, 0, -1.0);
coo.push(1, 1, 2.0);
coo.push(1, 2, -1.0);
coo.push(2, 1, -1.0);
coo.push(2, 2, 2.0);

let csr: Csr = coo.assemble();
let mut y = vec![0.0; 3];
csr.spmv(&[1.0, 2.0, 3.0], &mut y);
assert_eq!(y, vec![0.0, 0.0, 4.0]);
```

### Chebyshev Collocation and Stability Probes

```rust
use fs_cheb::orr_sommerfeld::{critical_reynolds, max_growth};

let growth = max_growth(5772.0, 1.02, 48).expect("eigen solve");
let critical = critical_reynolds(1.02, 48, 4000.0, 8000.0).expect("bracketed solve");

assert!(growth.is_finite());
assert!(critical > 5000.0);
```

### Evidence-Carrying Values

```rust
use fs_evidence::{Evidence, ProvenanceHash};

let provenance = ProvenanceHash::of_bytes(b"example kernel output");
let value = Evidence::exact(42.0, provenance)
    .certified()
    .expect("exact pure-math evidence is certifiable");

assert_eq!(value.value, 42.0);
```

### Repository Policy Checks

`xtask` is part of the workspace and emits JSON-lines verdicts so agents and DSR can parse results without scraping prose.

```bash
cargo run -p xtask -- check-layers
cargo run -p xtask -- check-deps
cargo run -p xtask -- check-contracts
cargo run -p xtask -- check-unsafe
cargo run -p xtask -- check-all
cargo run -p xtask -- check-constellation
```

## Architecture

The implemented crates follow a layered architecture. Lower layers do not depend on higher layers; `xtask check-layers` enforces that direction.

```text
                         fs-vskeleton
                              |
                +-------------+-------------+
                |                           |
             fs-plan                    fs-roofline
                |                           |
                +-------------+-------------+
                              |
                          fs-ledger
                              |
        +---------------------+---------------------+
        |                                           |
      fs-ir                                      fs-mesh
                                                   |
       +----------------------+--------------------+----------------------+
       |                      |                    |                      |
    fs-geom              fs-rep-sdf           fs-rep-mesh           fs-rep-frep
       |                      |                    |                      |
       +----------------------+---------+----------+----------------------+
                                        |
                                   fs-xform
                                        |
       +----------------------+---------+----------+----------------------+
       |                      |                    |                      |
    fs-math                fs-la              fs-sparse              fs-fft
       |                      |                    |                      |
    fs-ivl                fs-cheb              fs-rand                fs-ad
       |                                           |
       +----------------------+--------------------+
                              |
           +------------------+------------------+
           |                  |                  |
        fs-exec            fs-alloc           fs-simd
           |                  |                  |
           +------------------+------------------+
                              |
                         fs-substrate
                              |
                 fs-qty   fs-obs   fs-evidence
```

The diagram is simplified: some crates are siblings rather than strict parents, and `fs-qty`, `fs-obs`, and `fs-evidence` provide cross-cutting utility contracts. The authoritative dependency rule is the layer metadata in each crate manifest plus `xtask check-layers`.

## Design Principles

### Determinism Is a Contract

FrankenSim treats determinism as something each crate must state. Crate contracts identify determinism class, cancellation behavior, unsafe boundaries, and no-claim boundaries. The goal is not to pretend every floating-point path is identical on every machine; it is to say exactly what is stable, what is measured, and what evidence travels with a result.

### Evidence Travels With Results

The workspace uses explicit evidence wrappers instead of leaving proof obligations in comments. `fs-evidence` provides numerical, statistical, and model-form evidence; geometry conversion crates and ledger code can carry these records forward instead of discarding them at API boundaries.

### Geometry Is Representation-Aware

`fs-geom` defines chart and region contracts, while SDF, mesh, and F-rep crates implement specific representation behavior. This lets conversion, meshing, and transformation code reason about representation validity instead of treating geometry as an opaque blob.

### Runtime Behavior Is Inspectable

Execution contexts, cancellation gates, stream keys, tile pools, event schemas, and roofline measurement are regular workspace concepts. Long-running kernels are expected to be cancellable, observable, and attributable to ledgered artifacts and machine fingerprints.

### Repository Policy Is Code

The repo does not rely on prose alone for architecture rules. `xtask` enforces layer direction, dependency policy, contract presence, unsafe-capsule registration, and constellation lock consistency.

## Implemented Algorithms and Mechanics

The important point is not just that the crate names exist. The workspace already contains concrete algorithms and data structures with tests and contracts around them.

| Area | Implemented mechanics |
|------|-----------------------|
| Sparse assembly | `fs_sparse::Coo` stages triplets, canonicalizes them with deterministic duplicate accumulation, and emits CSR. CSR, BSR, and SELL kernels accumulate rows in a fixed column order so cross-format SpMV behavior is testable. |
| Randomized NLA | `fs-la::rand_nla` provides seeded range finding, randomized SVD, Nyström PSD approximation, sketch-and-precondition least squares, Hutchinson trace estimation, and Hutch++ trace estimation. Its golden-hash sentinel is still proof-pending. |
| Certified arithmetic | `fs-ivl` provides outward-rounded intervals, affine arithmetic, Taylor models, Newton/Krawczyk helpers, expansions, and exact predicates for orientation and in-sphere/incircle decisions. |
| Spectral numerics | `fs-cheb` implements adaptive 1D Chebyshev expansions, Lobatto grids, differentiation matrices, Dirichlet Laplacian checks, and an Orr-Sommerfeld stability battery. |
| Randomness and sampling | `fs-rand` uses counter-based Philox streams keyed by logical identity, with Sobol/QMC and lattice helpers so replay does not depend on thread arrival order. |
| AD and adjoints | `fs-ad` includes forward-mode dual numbers, a `Real` scalar contract, gradient checking, implicit-function hooks, and checkpointed/full-adjoint scaffolding. |
| Anytime inference | `fs-eproc` implements betting e-processes, pairwise races, Gaussian mixture confidence sequences, and e-BH support for conservative stopping decisions. |
| Geometry representations | `fs-geom` defines chart/region contracts; SDF, mesh, and F-rep crates implement concrete chart families, conversion records, quality checks, and representation-specific query behavior. |
| Transformations | `fs-xform` implements FFD lattices, RBF morphs, velocity bands, density fields, composed parameterizations, Jacobian actions, and foldover detection. |
| Meshing | `fs-mesh` contains tetrahedralization/refinement scaffolding, exact audits, ghost/refinement helpers, metric fields, remeshing, and execution-aware hooks. |
| Imaging | `fs-img` contains deterministic image artifact plumbing for the LUMEN layer: PNG/OpenEXR subset handling, film/display transforms, and explicitly bias-labeled denoising. |
| Execution | `fs-exec` provides cancellation contexts, stream keys, tile kernels, tile pools, deterministic reductions, racing/tuning scaffolding, solver state, and kill registries over `asupersync`. |
| Ledger | `fs-ledger` records content-addressed artifacts, operations, events, metrics, tune rows, extension rows, integrity checks, lineage, and time-travel queries on FrankenSQLite. |
| Roofline | `fs-roofline` probes machine axes, registers measurable kernels, records runs to the ledger, and checks staleness against machine fingerprints. |

## Contracts and Verification

The workspace currently has 28 `CONTRACT.md` files for 29 `fs-*` crates. The known exception is `fs-img`, which is present in the workspace but still needs a crate contract.

Existing contracts use these required sections:

- purpose and layer
- public types and semantics
- invariants
- error model
- determinism class
- cancellation behavior
- unsafe boundary
- feature flags
- conformance tests
- no-claim boundaries

Use these commands for local checks:

```bash
cargo fmt --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo run -p xtask -- check-all
```

Use DSR for the project validation path:

```bash
DSR_BIN="$(command -v dsr || printf '%s' /Users/jemanuel/projects/doodlestein_self_releaser/dsr)"
"$DSR_BIN" repos info frankensim
"$DSR_BIN" quality --tool frankensim
```

The current DSR setup is preferred over GitHub Actions for this repository. If a future GitHub workflow appears, do not treat it as authoritative unless `AGENTS.md` and DSR policy are changed first.

## Repository Layout

```text
.
|-- Cargo.toml                         # Workspace manifest
|-- Cargo.lock                         # Committed lockfile
|-- rust-toolchain.toml                # Nightly toolchain and components
|-- crates/
|   |-- fs-qty/                        # Dimensional quantities
|   |-- fs-obs/                        # Structured observability
|   |-- fs-evidence/                   # Evidence and certification wrappers
|   |-- fs-substrate/                  # Hardware/substrate probes
|   |-- fs-simd/                       # SIMD facades and capsules
|   |-- fs-alloc/                      # Aligned arenas and pools
|   |-- fs-exec/                       # Execution context, cancellation, tile work
|   |-- fs-math/                       # Deterministic scalar math
|   |-- fs-la/                         # Dense linear algebra
|   |-- fs-sparse/                     # Sparse formats and solvers
|   |-- fs-fft/                        # FFT and DCT transforms
|   |-- fs-ivl/                        # Certified arithmetic and predicates
|   |-- fs-cheb/                       # Chebyshev functions and collocation
|   |-- fs-rand/                       # Counter-based RNG and QMC
|   |-- fs-ad/                         # Automatic differentiation
|   |-- fs-eproc/                      # E-process inference
|   |-- fs-geom/                       # Chart/region abstractions
|   |-- fs-rep-sdf/                    # Signed-distance-field charts
|   |-- fs-rep-mesh/                   # Mesh charts and repair
|   |-- fs-rep-frep/                   # F-rep/CSG charts
|   |-- fs-xform/                      # Design parameterizations
|   |-- fs-img/                        # LUMEN image artifact plumbing
|   |-- fs-mesh/                       # Tetrahedralization and remeshing
|   |-- fs-ledger/                     # Design ledger on FrankenSQLite
|   |-- fs-ir/                         # FrankenScript IR
|   |-- fs-plan/                       # Cost and error planning
|   |-- fs-opt/                        # Optimization problem IR scaffold
|   |-- fs-roofline/                   # Roofline measurement harness
|   `-- fs-vskeleton/                  # Narrow PV vertical skeleton
|-- xtask/                             # Repository policy checks
|-- docs/                              # Conventions, CI/DSR gates, templates
|-- COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md
|                                       # Long-form architecture reference
|-- AGENTS.md                          # Agent operating rules
|-- unsafe-capsules.json               # Registered unsafe capsules
`-- constellation.lock                 # Franken sibling repository pin file
```

## How FrankenSim Compares

| Capability | FrankenSim today | Mature solver stacks | General scientific Python |
|------------|------------------|----------------------|---------------------------|
| End-user simulation app | Not yet | Yes, depending on stack | Usually assembled by user |
| Library substrate | Yes, source workspace | Yes | Yes |
| Determinism and cancellation contracts | Explicit crate contracts | Varies by project | Varies by package |
| Evidence-carrying values | Implemented core wrappers | Usually external/ad hoc | Usually external/ad hoc |
| Geometry representation contracts | Implemented SDF/mesh/F-rep chart layer | Often specialized | Usually library-specific |
| Design ledger | Implemented FrankenSQLite-backed crate | Often separate tooling | Usually separate notebooks/files |
| Roofline measurement | Implemented harness crate | Varies | Usually external tooling |
| Packaged distribution | Not yet | Often yes | Yes |

Use FrankenSim today when you want to work inside a deterministic, evidence-oriented Rust substrate for simulation infrastructure. Use a mature solver stack when you need a ready-made production physics solver or GUI today.

## Troubleshooting

### `failed to load source for dependency asupersync`

FrankenSim uses local path dependencies for the Franken constellation. Keep sibling repositories next to this checkout:

```text
~/projects/frankensim
~/projects/asupersync
~/projects/frankensqlite
```

Then rerun:

```bash
cargo test --workspace
```

### `the option Z is only accepted on the nightly compiler`

Use the pinned nightly toolchain:

```bash
rustup toolchain install nightly
rustup component add rustfmt clippy --toolchain nightly
cargo +nightly test --workspace
```

### `xtask check-deps` rejects a dependency

Runtime dependencies are intentionally restricted to the Franken constellation and workspace crates. Dev-only oracle dependencies may be allowed in test scope, but runtime additions need an architecture decision and matching contract update.

```bash
cargo run -p xtask -- check-deps
```

### `xtask check-contracts` fails

Every completed `fs-*` crate must ship a `CONTRACT.md` with the required sections. The current tree has one known gap, `fs-img`, because that crate is still an in-progress scaffold. Use an existing crate contract as the local template and keep the no-claim boundaries honest.

```bash
cargo run -p xtask -- check-contracts
```

### DSR and local commands disagree

Treat DSR as the project automation source of truth. Local commands are useful for fast iteration; DSR is the lane to cite when reporting repository health.

```bash
DSR_BIN="$(command -v dsr || printf '%s' /Users/jemanuel/projects/doodlestein_self_releaser/dsr)"
"$DSR_BIN" quality --tool frankensim
```

## Limitations

FrankenSim has substantial working code, but it is still early infrastructure.

| Capability | Current state |
|------------|---------------|
| Stable public API | Not promised yet; contracts exist, but APIs may still change |
| End-user CLI/application | Not implemented as a packaged simulation app |
| crates.io distribution | Not published |
| GitHub Actions | Not authoritative for this repo; use DSR |
| Full multiphysics solver suite | Not complete in the current workspace |
| NURBS, voxel, and neural representations | Not implemented as first-class crates in the current workspace |
| Randomized NLA golden sentinel | `fs-la::rand_nla` is implemented and mostly covered, but `rand_nla_golden_hash` still carries a placeholder hash until the observed value is deliberately recorded and cross-ISA proof is completed |
| Production validation corpus | In progress through contracts, tests, ledger records, and roofline harnesses |
| Performance claims | Must be backed by `fs-roofline`/ledger evidence; do not infer claims from architecture text alone |

The long-form architecture reference remains in `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md`, but this README describes the code that is already present in the repository.

## FAQ

### Is FrankenSim usable today?

Yes, as a Rust source workspace and simulation-infrastructure library substrate. It is not yet a polished end-user simulator.

### What should I run first?

Start with the workspace tests and policy checks:

```bash
cargo test --workspace
cargo run -p xtask -- check-all
```

For the project validation lane, use:

```bash
DSR_BIN="$(command -v dsr || printf '%s' /Users/jemanuel/projects/doodlestein_self_releaser/dsr)"
"$DSR_BIN" quality --tool frankensim
```

### Why nightly Rust?

The repository pins nightly in `rust-toolchain.toml` for Rust 2024 plus narrow nightly features documented in the repo, including const-generic dimension arithmetic and optional portable-SIMD experiments. Most code is written so future de-nightlying should be localized rather than a rewrite.

### Why are local Franken dependencies required?

The workspace intentionally depends on other Franken projects by path for runtime pieces such as `asupersync` and `frankensqlite`. `xtask check-constellation` and `constellation.lock` make that relationship explicit.

### Are the crate contracts normative?

Yes. Existing `CONTRACT.md` files are part of how the workspace communicates implemented behavior. If code and contract disagree, that is a bug to resolve, not a documentation detail to ignore. The current exception is `fs-img`, which still needs its first contract.

### Should I use GitHub Actions?

No. Use DSR for this repository unless the repository policy changes. GitHub Actions may be absent, stale, or non-authoritative.

## About Contributions

Please don't take this the wrong way, but I do not accept outside contributions for any of my projects. I simply don't have the mental bandwidth to review anything, and it's my name on the thing, so I'm responsible for any problems it causes; thus, the risk-reward is highly asymmetric from my perspective. I'd also have to worry about other "stakeholders," which seems unwise for tools I mostly make for myself for free. Feel free to submit issues, and even PRs if you want to illustrate a proposed fix, but know I won't merge them directly. Instead, I'll have Codex or Codex review submissions via `gh` and independently decide whether and how to address them. Bug reports in particular are welcome. Sorry if this offends, but I want to avoid wasted time and hurt feelings. I understand this isn't in sync with the prevailing open-source ethos that seeks community contributions, but it's the only way I can move at this velocity and keep my sanity.

## License

FrankenSim is licensed under the MIT License with the repository's Anthropic/OpenAI rider. See [LICENSE](LICENSE) for the exact terms.
