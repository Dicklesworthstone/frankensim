<div align="center">
  <img src="frankensim_illustration.webp" alt="FrankenSim illustration" width="720">
</div>

# FrankenSim

<div align="center">

[![Status](https://img.shields.io/badge/status-active%20Rust%20workspace-2ea44f)](#implemented-workspace)
[![Rust](https://img.shields.io/badge/rust-nightly%202024-b7410e)](rust-toolchain.toml)
[![Crates](https://img.shields.io/badge/workspace-126%20fs--%2A%20crates-0969da)](#implemented-workspace)
[![Contracts](https://img.shields.io/badge/contracts-127%20of%20127%20crates-8250df)](#contracts-and-verification)
[![Tests](https://img.shields.io/badge/tests-263%20crate%20test%20files-1f883d)](#contracts-and-verification)
[![License](https://img.shields.io/badge/license-MIT%20%2B%20AI%20rider-yellow)](LICENSE)

</div>

FrankenSim is a working Rust workspace for deterministic geometry, certified numerics, meshing, execution, evidence, and design-ledger infrastructure for simulation and design optimization.

The tree contains 127 `fs-*` crate directories: 126 in the native Cargo workspace plus the standalone nested `fs-wasm` workspace. They include repository policy tooling, conformance contracts, integration tests, and working implementations across substrate/runtime, numerical kernels, geometry representations, meshing, physics, solvers, adjoints, optimization, imaging, evidence, packaging, and ledger layers.

There is not yet a packaged end-user simulation application or crates.io release. Today, FrankenSim is usable as a source workspace and library substrate.

## TL;DR

**The problem:** simulation systems often split physical units, numerical error, runtime behavior, geometry validity, evidence, and reproducibility across separate tools. That makes it too easy for an optimization run to produce an answer without a durable explanation of which assumptions, approximations, kernels, and machine conditions made the answer valid.

**The solution:** FrankenSim builds those concerns into the workspace architecture. Units are represented explicitly, kernels have deterministic contracts, geometry conversions carry evidence, runtime behavior is structured around cancellable work contexts, and the ledger records artifacts, operations, events, roofline measurements, and time-travelable design state.

### What Exists Now

| Area | Current implementation |
|------|------------------------|
| Workspace | Rust 2024 nightly Cargo workspace with 126 native `fs-*` workspace crates plus `xtask`; `fs-wasm` is a standalone nested workspace |
| Contracts | 127 of 127 `fs-*` crate directories have `CONTRACT.md` files |
| Runtime substrate | Capability probing, SIMD facades, aligned arenas, two-lane execution, cancellation contexts, tile pools, tuner and race scaffolding |
| Numerics | Deterministic elementary math, dense/sparse linear algebra, FFT/DCT, interval/affine/Taylor arithmetic, Chebyshev collocation, random/QMC streams, AD/adjoint infrastructure, e-process inference |
| Geometry | Region/chart abstraction, SDF, mesh and F-rep charts, representation conversion hooks, transformations, tet meshing, remeshing, quality audits |
| Evidence and ledger | Composable `Evidence<T>`/`Certified<T>`, model cards, bracketing, FrankenSQLite-backed design ledger, artifact hashes, event streams, tune cache, roofline recording |
| Policy tooling | `xtask` checks for layer direction, Franken-only runtime dependencies, contracts, unsafe capsules, and constellation lock verification |
| Tests | 263 crate-level conformance and integration test files in the intended snapshot, exercising the implemented contracts |

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
| Work with NURBS, voxel, point-cloud, and lattice representations | `fs_rep_nurbs` and `fs_rep_voxel` |
| Assemble exterior-calculus operators and cohomology probes | `fs_feec` cochains, incidence, Hodge stars, Betti checks, and Hodge decomposition |
| Run matrix-free/resumable linear solves | `fs_solver` CG, MINRES, GMRES, transposed solves, and p-multigrid support |
| Build CutFEM and elasticity fixtures | `fs_cutfem`, `fs_solid`, `fs_material`, and `fs_scenario` |
| Verify gradients and build adjoint pipelines | `fs_ad`, `fs_adjoint`, and the gradient verification gate |
| Track spectral health signals | `fs_spectral` sheaf-Laplacian gap monitoring and low-confidence propagation |
| Harvest cyclic symmetry in solves | `fs_symmetry` cyclic residuals, DFT block-diagonalized circulant solves, and perturbation bounds |
| Run gradient-based, derivative-free, or density-based optimizers | `fs_ascent`, `fs_dfo`, `fs_constraint`, `fs_topopt`, and `fs_opt` |
| Keep objective uncertainty honest | `fs_robust` CVaR, weakest-input colors, robust optima, and colored fragility curves |
| Apply design parameterizations and detect foldover | `fs_xform` FFD, RBF, velocity band, density, and composition types |
| Encode image artifacts and apply deterministic film/denoising transforms | `fs_img` PNG/OpenEXR subset plumbing, film transforms, and bias-labeled denoising |
| Record artifacts, operations, metrics, events, and tune rows | `fs_ledger` on FrankenSQLite |
| Package solver evidence for independent review | `fs_package`, `fs_checker`, and `fs_crosswalk` |
| Check plugin conformance and restriction-map behavior | `fs_conform` restriction-map compatibility checks and diagnostics |
| Keep the go-to-market wedge explicit | `fs_wedge` vertical ranking, proposal mapping, and cycle-time kill-criterion data |
| Probe machine axes and record roofline measurements | `fs_roofline` plus `fs_substrate` |

## Quick Start

FrankenSim currently builds from source. The workspace expects exact pinned
sibling Franken projects for path dependencies, especially
`~/projects/asupersync` and `~/projects/frankensqlite`. On a fresh checkout,
materialize those siblings first with the standalone, zero-dependency bootstrap
at `tools/bootstrap`; it is deliberately outside the root workspace so Cargo
does not have to resolve the missing path dependencies before building it. The
current FrankenNumpy pin is case-collision-safe on default macOS filesystems.
See [`docs/BOOTSTRAP.md`](docs/BOOTSTRAP.md) for offline and mirror modes, the
fail-closed verification contract, and the remaining blank-machine evidence
boundary.

```bash
git clone https://github.com/Dicklesworthstone/frankensim.git
cd frankensim

# Use the pinned nightly toolchain from rust-toolchain.toml.
rustup toolchain install nightly
rustup component add rustfmt clippy --toolchain nightly

# Materialize every sibling at the exact constellation.lock pin.
cargo run --manifest-path tools/bootstrap/Cargo.toml

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
| `cargo run --manifest-path tools/bootstrap/Cargo.toml` | From a fresh checkout, materialize and verify the pinned sibling constellation before the root workspace resolves |
| `cargo test --workspace` | Build the workspace and run crate-level unit, integration, and conformance tests |
| `cargo fmt --check` | Check formatting under the pinned nightly toolchain |
| `cargo clippy --workspace --all-targets -- -D warnings` | Run the strict local lint lane |
| `cargo run -p xtask -- check-layers` | Enforce the layer dependency direction |
| `cargo run -p xtask -- check-deps` | Enforce the Franken-only runtime dependency policy |
| `cargo run -p xtask -- check-contracts` | Verify `fs-*` crate contracts and required sections |
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
| `fs-ledger` | Versioned Design Ledger schema (currently v5) on FrankenSQLite, immutable persisted ledger identity, content-addressed artifacts, operations, event streams, lineage, metrics, tune cache, extension rows, integrity/lint checks, and time-travel queries |
| `fs-ir` | FrankenScript typed AST with spans, isomorphic s-expression and JSON syntaxes, shape comparison, study recognition, lowering, and structured IR errors |
| `fs-plan` | Cost model types, error and time ledgers, plan-cost oracle, and cost-model construction from tune records |
| `fs-roofline` | Machine-axis probing, protected historical baselines, kernel registry, measurement receipts, baseline-bound ledger recording, and staleness checks |
| `fs-img` | LUMEN image pipeline crate with deterministic PNG/OpenEXR subset plumbing, film/display transforms, bias-labeled denoising, and conformance tests |
| `fs-vskeleton` | Photovoltaic vertical skeleton tying SDF/PDE/objective/adjoint/optimization/ledger concepts into a narrow demonstrator path |
| `fs-opt` | ASCENT optimization problem IR crate with typed objective/constraint graphs, dimensional validation, differentiability-class routing, canonical serialization, manifold metadata, and a conformance contract |

### Expanded Crate Map

The workspace has grown beyond the first substrate and geometry layer. These crates are also implemented today and fill in the physics, optimization, audit, and orchestration layers.

| Crate family | Implemented role |
|--------------|------------------|
| `fs-query`, `fs-topo`, `fs-geocon` | Geometry query and certificate layer: closest-point, raycast, offset, clearance, thickness, curvature, topology validity, and geometric constraints |
| `fs-rep-nurbs`, `fs-rep-voxel` | Additional MORPH representations: rational B-spline charts with exact spline algebra and certified closest-point brackets; voxel, point-cloud, and lattice/strut representations |
| `fs-feec` | Exterior-calculus core: cochains, Whitney forms, exact incidence operators, Hodge stars, Betti checks, high-order tensor/simplex spaces, and cohomology/Hodge-decomposition utilities |
| `fs-opdsl`, `fs-tilelang`, `fs-tilelang-macros`, `fs-soa`, `fs-soa-derive` | Operator and layout tooling: DSL scaffolding for operators, deterministic tile-language lowering, and structured data layouts for stable kernel memory behavior |
| `fs-solver` | Resumable Krylov and multigrid stack: CG, MINRES, GMRES, matrix-free operators, transposed solves, mixed-precision refinement, and p-multigrid over the FEEC hierarchy |
| `fs-cutfem`, `fs-solid`, `fs-material`, `fs-scenario`, `fs-time` | FLUX physics layer: CutFEM on certified SDF cuts, elasticity/hyperelasticity, constitutive laws, typed boundary/load scenarios, and time-support primitives |
| `fs-adjoint` | Gradient truth layer: IFT adjoints, density/SIMP pullbacks, Sobolev smoothing, Hadamard shape gradients, revolve-style time adjoints, gradient verification, and feature-gated ledger transposition |
| `fs-spectral`, `fs-symmetry` | Spectral and symmetry health tools: sheaf-Laplacian gap tracking, hysteresis, confidence propagation, cyclic residuals, DFT block-diagonalized circulant solves, and perturbation bounds |
| `fs-constraint`, `fs-ascent`, `fs-dfo`, `fs-dimine`, `fs-regime`, `fs-topopt`, `fs-robust` | ASCENT layer: typed constraints with unsat cores and repairs, L-BFGS/TR-Newton/augmented-Lagrangian/Riemannian optimizers, derivative-free CMA-ES/BIPOP/Nelder-Mead, dimensional law mining, regime/validity-domain machinery, density-based topology optimization, and objective epistemics |
| `fs-contract`, `fs-verify`, `fs-bisect`, `fs-conform` | Verification layer: assume-guarantee component contracts, certified speculation verification, bisect-style proof/debug infrastructure, and restriction-map plugin conformance |
| `fs-package`, `fs-checker`, `fs-crosswalk` | Evidence distribution layer: content-addressed evidence packages, a solver-free standalone checker, and a standards-language crosswalk for ASME V&V and FAA/EASA certification-by-analysis concepts |
| `fs-recompute`, `fs-session`, `fs-govern`, `fs-probe`, `fs-iface`, `fs-ladder`, `fs-spececo`, `fs-wedge`, `fs-io` | HELM and utility layer: incremental recomputation, sessions/resource governance, machine-readable proposal governance, probes, interfaces, degradation ladders, specialization ecology, wedge-selection data, and I/O boundaries |

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

### Evidence Packages and Standalone Checking

```rust
use fs_checker::check_against_root;
use fs_package::{Claim, EvidencePackage, Provenance};

let claim = Claim::estimated(
    "lift-preview",
    "surrogate lift estimate before wind-tunnel admission",
    "pod-surrogate-v3",
    0.04,
);

let pkg = EvidencePackage::new(Provenance::new("commit-abc", "lock-def"))
    .with_claim(claim);

let root = pkg.try_merkle_root().expect("bounded example package");
let report = check_against_root(&pkg, root);

assert!(report.passed());
assert_eq!(report.merkle_root(), root);
assert!(report.render_pie().contains("estimated"));
```

This deny-all example checks an honestly Estimated claim and the expected
content root. Verified certificates, anchoring datasets, waivers, signatures,
and falsifier artifacts never authorize themselves: use the corresponding
typed verifier capabilities and retain the returned verification receipt.
Scientific release admission additionally requires a purpose-bound approval
signature over the exact checker protocol, expected root, scientific policy
fingerprints, waiver clock, and per-claim admission context, at least one
admitted scientific Verified or Validated claim, and authenticated falsifier
evidence for every certificate-class claim. It does not infer signer identity,
role, or authorship from detached signature bytes.

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
| Randomized NLA | `fs-la::rand_nla` provides seeded range finding, randomized SVD, Nyström PSD approximation, sketch-and-precondition least squares, Hutchinson trace estimation, and Hutch++ trace estimation. Its golden-hash sentinel is recorded and verified identical across both reference ISAs and both build modes. |
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
| Roofline | `fs-roofline` probes machine axes, checks pre/post probes against a protected historical baseline, records citable runs under the machine-plus-baseline identity, and reports either kind of drift. |

## Why This Is Useful

FrankenSim is not trying to be another wrapper around an external solver. The useful part is that the same workspace owns the path from representation to evidence. A result can carry dimensional meaning, numerical method, solver status, derivative provenance, validity regime, machine fingerprint, and package-level audit data without crossing an untyped boundary.

| Common failure mode | FrankenSim's implemented answer |
|---------------------|----------------------------------|
| A quantity loses units once it enters a kernel | `fs-qty`, `fs-scenario`, and `fs-opt` keep dimensional information explicit at API boundaries |
| A geometry conversion silently changes validity | `fs-geom`, representation crates, `fs-query`, and `fs-topo` expose chart contracts, conversion records, and validity/topology checks |
| A solve returns a number but not a reason to trust it | `fs-solver`, `fs-adjoint`, `fs-verify`, and `fs-evidence` return residuals, certificates, brackets, validation domains, and evidence colors |
| An optimizer stalls and nobody knows why | `fs-ascent`, `fs-dfo`, and `fs-constraint` report stop reasons, KKT residuals, unsat cores, repair suggestions, and restart schedules |
| A benchmark only works on one machine | `fs-substrate`, `fs-roofline`, `fs-ledger`, and DSR make machine axes, kernels, measurements, and validation lanes explicit |
| An auditor wants the evidence, not the solver | `fs-package`, `fs-checker`, and `fs-crosswalk` package claims, re-check them without the solver stack, and describe them in existing certification vocabulary |

The design target is a workflow where "what did we compute?" and "why should anyone believe it?" are not separate artifacts.

## End-to-End Flow

The crates are deliberately small, but the intended flow is continuous. A design study can move from typed intent to geometry, physics, derivatives, optimization, and audit without dropping provenance.

```text
typed intent / study IR
        |
        v
fs-ir, fs-opt, fs-scenario
        |
        v
geometry and representation layer
fs-geom -> fs-rep-sdf / fs-rep-mesh / fs-rep-nurbs / fs-rep-voxel / fs-rep-frep
        |
        v
validity, queries, topology, and meshing
fs-query -> fs-topo -> fs-mesh -> fs-feec / fs-cutfem / fs-solid
        |
        v
operators and solves
fs-opdsl / fs-tilelang -> fs-sparse / fs-la -> fs-solver
        |
        v
derivatives and optimization
fs-ad / fs-adjoint -> fs-ascent / fs-dfo / fs-constraint
        |
        v
evidence, replay, and audit
fs-evidence -> fs-ledger -> fs-package -> fs-checker -> fs-crosswalk
```

That flow is still library-level rather than an end-user app, but the important seams already exist in code: geometry knows its chart, physics has typed scenarios and solver reports, adjoints are checked against finite differences, optimizers return certificates, and evidence packages can be verified without running the solver.

## Algorithmic Spine

### Geometry and Topology

FrankenSim uses representation-specific charts instead of flattening all geometry into one mesh too early. SDF grids, triangle/tet meshes, rational NURBS patches, voxel fields, point clouds, lattices, and F-reps each keep their own query and validity contracts. Conversion is treated as an operation with a record, not as a silent cast.

Important mechanics already implemented include exact rational spline refinement in `fs-rep-nurbs`, sparse VDB-style fields and exact Euclidean distance transforms in `fs-rep-voxel`, mesh repair and winding support in `fs-rep-mesh`, CSG/F-rep evaluation in `fs-rep-frep`, and topology/validity checks in `fs-topo`.

### Discretization and Physics

The FLUX layer is built around discrete exterior calculus and explicit weak forms. `fs-feec` owns cochains, incidence, Hodge stars, high-order FEEC spaces, Betti checks, and Hodge decomposition. `fs-cutfem` works directly on certified SDF cuts with ghost penalty, Nitsche embedded boundary conditions, and cut quadrature. `fs-solid` provides small-strain and finite-strain elasticity, B-bar locking relief, CutFEM/body-fitted frontends, and Newton globalization against constitutive laws from `fs-material`.

The principle is that numerical stabilization should either follow from the complex/operator structure or be recorded as a named method with tests. Spurious modes, checkerboarding, cut-cell conditioning, locking, and buckling-adjacent tangents are treated as design constraints rather than after-the-fact patches.

### Solvers, Adjoints, and Optimization

The solver stack is resumable and transposition-aware. `fs-solver` implements CG, MINRES, GMRES, matrix-free operators, mixed-precision refinement, and p-multigrid with deterministic reductions. `fs-adjoint` builds gradients by transposed solves and IFT machinery rather than differentiating through Krylov iterations. The feature-gated ledger-transposition module records a DAG and pulls cotangents backward through registered VJPs, refusing missing or declared non-differentiable seams loudly instead of returning a silent zero.

On top of that, `fs-ascent` provides L-BFGS with strong-Wolfe search, trust-region Newton-Krylov with Steihaug negative-curvature handling, augmented-Lagrangian constraints with KKT certificates, Riemannian L-BFGS over manifold metadata, and a stop-rule algebra that attributes why a run ended. `fs-dfo` adds deterministic CMA-ES, BIPOP restarts, and Nelder-Mead for cases where gradients are unavailable or untrusted.

### Evidence, Packages, and Standards

`fs-evidence` gives values a claim type. `fs-package` groups color-typed claims and provenance into a content-addressed Merkle bundle. `fs-checker` re-verifies package completeness, content address, signature status, and budget breakdown without depending on solvers, geometry, or license gates. Its distinct release gate additionally requires an injected signature verifier, falsifier pairing, and matching anchors; ordinary integrity checking does not authenticate a signature merely because text is present. `fs-crosswalk` maps the same package concepts onto ASME V&V 10/20/40 and FAA/EASA certification-by-analysis vocabulary.

That division matters: the expensive solver can produce evidence once, while a reviewer can later check the package structure and provenance without re-running the solver stack.

## Design Rules Encoded In The Code

The README-level principles are backed by concrete implementation patterns:

| Rule | Where it shows up |
|------|-------------------|
| Make invalid states loud | Typed `Error`/refusal enums in geometry, packages, constraints, solvers, and scenarios; missing VJPs and non-differentiable seams block gradients |
| Keep replay independent of arrival order | Philox stream keys, deterministic reductions, sorted assembly, fixed traversal order, and content-addressed artifacts |
| Put budgets in data | Stop rules, solver iteration reports, cut quadrature depth, roofline runs, package color breakdowns, and session resource estimates |
| Separate production code from oracle code | Franken-only runtime dependency policy plus isolated test/oracle paths checked by `xtask` |
| Prefer certificates to comments | `CONTRACT.md` files, `fs-evidence`, gradient verification gates, package checks, topology checks, and ledger receipts |
| Treat performance as a measured claim | `fs-roofline`, substrate machine fingerprints, tune rows, and DSR validation lanes |
| Avoid hidden global state | Explicit contexts, stream keys, scenarios, manifests, packages, and ledger records |

## Why The Pieces Live Together

Most simulation systems split their real work across several unrelated layers: a
CAD kernel owns geometry, a mesher owns discretization, a solver owns residuals,
an optimizer owns design changes, a notebook owns provenance, and a separate CI
system tries to prove the run still works. FrankenSim keeps those concerns in one
typed Rust workspace because the hard bugs usually happen at the boundaries.

| Boundary that usually leaks | Why it matters | FrankenSim's approach |
|-----------------------------|----------------|-----------------------|
| Geometry to mesh | A conversion can change topology, watertightness, or feature size before physics sees it | Chart and representation crates expose conversion records, validity checks, topology summaries, and no-claim boundaries |
| Mesh/operator to solver | A residual can be small for the wrong norm, wrong operator, or wrong boundary condition | Operator crates, solver reports, and contracts name the algebraic object being solved and the residual semantics |
| Solver to adjoint | A gradient can silently differentiate the wrong fixed point or ignore a transpose seam | `fs-adjoint` uses IFT/transposed-solve machinery and feature-gated ledger transposition that refuses missing VJPs |
| Optimizer to evidence | An optimizer can return a design without explaining feasibility, stopping reason, or uncertainty | `fs-ascent`, `fs-dfo`, `fs-constraint`, and `fs-robust` report KKT residuals, stop reasons, CVaR/DRO quantities, unsat cores, and evidence colors |
| Experiment to artifact | A report can preserve plots but lose code version, machine, budgets, or assumptions | `fs-ledger`, `fs-report`, `fs-package`, and `fs-checker` keep content hashes, provenance, budget pies, package roots, and standalone checks |

The result is not just "many crates." The result is a set of small interfaces
that make the expensive parts of simulation inspectable: what representation was
used, what numerical contract was in force, what was measured, what was
estimated, and which claim boundary still needs proof.

## What "Certified" Means Here

FrankenSim uses certification language carefully. A value is not "certified"
because a test passed somewhere; it is certified only when the code path returns
an explicit claim with the right evidence color, provenance, and budget context.

| Term | Meaning in this repo | Typical source |
|------|----------------------|----------------|
| `Verified` | A property was checked by a deterministic proof, exact computation, residual certificate, invariant, or bounded numerical argument inside the modeled assumptions | Exact predicates, solver residual certificates, KKT residuals, package integrity checks, impulse-conservation bands |
| `Validated` | A model or result is supported against an external or empirical reference over a stated validity domain | Model cards, benchmark fixtures, calibration data, crosswalk records |
| `Estimated` | The claim is useful but relies on statistical, model-form, discretization, or approximation evidence that is not upgraded to a hard bound | DWR estimates, stochastic confidence sequences, surrogate bands, model-form evidence |

The important design rule is **no laundering**: composing a `Verified` value with
an `Estimated` value cannot produce a stronger result. That shows up in the
evidence wrappers, the robust/fragility tooling, and the end-to-end smoke
campaigns. For example, `fs-thrust-e2e` can place verified bands around
full-fidelity point-vortex runs that conserve impulse, while surrogate-served
elites remain estimated. The campaign-level rank is the weakest elite color, not
the most flattering one.

## Representative Vertical Slices

Several crates are intentionally small capstones. They are not full products,
but they prove that lower-level pieces compose into end-to-end workflows.

| Slice | What it demonstrates | Why it matters |
|-------|----------------------|----------------|
| `fs-marquee` | A feature-gated raw-SDF/CutFEM/DWR smoke runner for a topology-style layout study | Geometry, embedded-boundary physics, gradient checks, and certificate fields can live on the same path without meshing as the first operation |
| `fs-frame` | A seismic-minimal frame smoke tier: truss layout, sizing, fiber hysteresis, e-stopped fragility, and CVaR sizing | Structural optimization can carry layout diagnostics, nonlinear response evidence, anytime-valid probability intervals, and robust tail metrics together |
| `fs-thrust-e2e` | CertQD-Thrust: point-vortex thruster simulation, certify-or-escalate surrogate use, MAP-Elites illumination, evidence colors, and a content-addressed notebook | Quality-diversity search can produce an atlas with per-elite trust labels instead of a single unqualified best design |
| `fs-vskeleton` | A photovoltaic vertical skeleton connecting SDF/PDE/objective/adjoint/optimization/ledger ideas | A narrow domain slice is useful for testing orchestration seams before every physics model is mature |

These slices are deliberately honest about scope. They are smoke tiers and
integration fixtures, not claims that every flagship has reached its full
production resolution. Their value is that they exercise the architecture under
real composition pressure: dependencies point downward, evidence has to survive
multiple layers, and tests can catch mismatches between the story and the code.

## Algorithm Families In Practice

FrankenSim combines conservative numerical machinery with explicit
engineering-grade fallbacks. The algorithms are implemented where they give the
workspace a contract it can enforce, not just because they are impressive names.

| Family | Implemented examples | Design principle |
|--------|----------------------|------------------|
| Strict elementary math | `fs-math` deterministic `exp`, `ln`, `sin`, `cos`, `tan`, `atan2`, `erf`, `pow`, double-double helpers, and Payne-Hanek large-argument trig reduction | Low-level math must not delegate determinism to platform libm when higher layers depend on bit-stable replay |
| Exact and bounded predicates | `fs-ivl` intervals, affine/Taylor arithmetic, expansions, `orient`, `incircle`, and `insphere` | Geometry and topology decisions should be made with explicit uncertainty instead of hoping floating-point signs are stable |
| Spectral numerics | `fs-cheb` adaptive Chebyshev functions, collocation, differentiation matrices, and Orr-Sommerfeld probes | Smooth 1D and stability problems get compact, high-order representations with conformance batteries |
| Sparse and structured linear algebra | CSR/BSR/SELL sparse formats, deterministic assembly, randomized NLA, eigensolver scaffolding, Krylov solvers, and p-multigrid | Large solves need reproducible assembly order, inspectable residuals, and operator structure that adjoints can reuse |
| Optimization | L-BFGS, trust-region Newton-Krylov, augmented Lagrangian, Riemannian L-BFGS, CMA-ES, BIPOP, Nelder-Mead, NSGA-II, Pareto tracing, CVaR, and discrete-support Wasserstein DRO | Different design regimes need different optimizers, but all of them should report why they stopped and what assumptions their answer uses |
| Anytime and robust inference | Betting e-processes, pairwise races, e-BH, Gaussian mixture confidence sequences, conformal bands, CVaR, and DRO inner suprema | A study should stop when the decision is good enough, not when a preselected sample count happens to be exhausted |
| Evidence packaging | Claim colors, validity domains, Merkle roots, package checks, crosswalk mappings, and budget pies | Solver output should be reviewable as a standalone artifact after the expensive run is over |

Two recent examples show the pattern. `fs-math::payne` extends trig reduction
beyond the old Cody-Waite domain with self-verifying `2/pi` constants generated
from an all-integer Machin bignum in tests. `fs-dfo::dro` implements the exact
dual for a discrete-support Wasserstein ambiguity set and tests it against
closed-form endpoints, kink cases, and tiny LP oracles. In both cases, the code
does not merely add an algorithm; it adds the tests and contract language that
make the algorithm usable elsewhere in the workspace.

## Core Abstractions That Repeat Across Crates

FrankenSim is easier to understand if you look for the same few shapes
reappearing in different layers. The crate names are numerous, but the design is
not arbitrary: values, fields, operators, budgets, and evidence all have
standard carriers.

| Abstraction | Where it appears | What it buys |
|-------------|------------------|--------------|
| `Qty` and dimensional records | `fs-qty`, `fs-scenario`, `fs-opt`, material and constraint crates | A boundary condition, objective, or material parameter can reject the wrong unit before it contaminates a solve |
| `Cx`, cancellation gates, and stream keys | `fs-exec`, solver drivers, stochastic crates, e2e campaigns | Work can be interrupted, replayed, and attributed without depending on thread arrival order |
| Chart and region traits | `fs-geom`, `fs-rep-*`, `fs-query`, `fs-topo` | Geometry keeps representation-specific meaning instead of becoming an untyped point cloud too early |
| Linear operators and transpose actions | `fs-solver`, `fs-adjoint`, `fs-bem`, `fs-opdsl` | Primal solves, adjoint solves, and gradient checks can share the same operator contract |
| Evidence colors | `fs-evidence`, `fs-package`, campaign crates, robust/validation crates | A result says whether it is verified, validated, or estimated, and composition cannot upgrade the weakest ingredient |
| Contract files | Every `fs-*` crate | Each crate states public semantics, invariants, determinism, feature flags, tests, and no-claim boundaries in one place |
| Content roots and package checks | `fs-ledger`, `fs-report`, `fs-package`, `fs-checker` | A later reader can inspect what was claimed without reconstructing the whole run from a notebook |

This repetition is intentional. A new physics module, optimizer, or campaign
does not need to invent its own story for cancellation, uncertainty, provenance,
or validation. It plugs into the same carriers and either satisfies their
contracts or refuses loudly.

## Evidence-Bearing Campaign Catalog

The e2e crates are small, but they are the best way to see the workspace acting
as one system. Each campaign composes several lower-level crates and forces a
real answer to carry evidence instead of just a number.

| Campaign crate | Composition | Decision or artifact produced |
|----------------|-------------|-------------------------------|
| `fs-thrust-e2e` | Vortex particle dynamics, conformal surrogate gating, MAP-Elites, evidence colors, report generation | A quality-diversity atlas of self-propelling vortex-thruster designs with per-elite trust labels |
| `fs-robustopt-e2e` | SOS proof machinery, robust objective epistemics, evidence wrappers | A design choice that is both globally certified in the nominal model and robust under perturbation |
| `fs-schedule-e2e` | Tropical critical paths, value-of-information, evidence colors | A campaign schedule with named bottlenecks and an act/stop recommendation |
| `fs-flutter-e2e` | Coupling, spectral health, SOS Lyapunov checks, evidence | A flutter boundary and stability certificate for a partitioned added-mass model |
| `fs-neuroshape-e2e` | Neural implicit charts, Lipschitz bounds, interval topology checks, visualization evidence | A certified neural-SDF shape with no-tunnel tracing and bounded-component evidence |
| `fs-grammar-e2e` | Shape-program synthesis, archive illumination, fabrication checks, evidence | A diverse family of simplified CSG programs that match a target while retaining certificates |
| `fs-oed-e2e` | Kalman assimilation, value-of-information, tolerance allocation, evidence | A sensor placement plan that stops when the design decision is robust |
| `fs-metamat-e2e` | Lattice homogenization, SOS PSD checks, evidence | A stiffness-density frontier whose points are stable and Voigt-admissible |
| `fs-truss-e2e` | Ground-structure truss optimization, endpoint-checked tropical load paths, evidence | A deterministic approximate truss iterate plus an advisory connected load path and conditionally unique bottleneck member |
| `fs-adaptbo-e2e` | Bayesian optimization, e-process stopping, confidence sequences, evidence | An anytime BO run that can stop early under optional-stopping-valid evidence |
| `fs-flowcert-e2e` | LBM, analytic Poiseuille checks, MAP-Elites, evidence | A CFD credibility map over Reynolds number and resolution |
| `fs-vessel` | Chebyshev stability objective, free-surface LBM, e-racing, robust CVaR, volume rendering | A laminar-pour vessel study with stability, mass-ledger, robust-family, and render evidence |
| `fs-ornith` | BEM/VPM screening, e-racing, LBM refinement, SOS stability certificates, conformal surrogate bands, NSGA-II atlas construction | A smoke-tier ornithoid aircraft study that carries lineage, certificates, model-form honesty, and replay evidence through a full design campaign |
| `fs-topopt` `cutfem-marquee` | Density fields as CutFEM SDFs, DWR-guided octree refinement, zero-remesh logs, topology/thickness witnesses | A feature-gated topology-optimization lane that keeps the background grid stable while recording exactly where the proof work still lives |
| `fs-flagship-e2e` | Smoke/mid/full stage wiring, metric-only golden hashes, cross-flagship audits, failure drills, forensic JSON rows, notebook replay | A replay suite with frozen smoke-stage goldens, shared-core drift checks, and ignored mid/full fidelity lanes for later perf cadence |

The common pattern is: generate candidates, evaluate them with a physics or
numerical kernel, attach evidence to the result, and stop only when the decision
or artifact is good enough under the declared budget. This is why the campaigns
are useful even when they are smoke-tier: they stress the boundaries between
representations, solvers, optimizers, reports, and evidence.

## Recent Working Surfaces Worth Reading

The fastest way to understand the codebase is to read a few recently landed
surfaces that each compress a larger design rule into a small amount of code.
They show how FrankenSim tries to turn "a solver produced a number" into "a
claim survived a particular chain of assumptions."

| Surface | What was built | Why it is useful |
|---------|----------------|------------------|
| `fs-bo::sparse` | Inducing-point sparse Gaussian processes with deterministic farthest-point inducing locations, DTC/SoR prediction, Titsias ELBO accounting, and sparse/exact conformance tests | Bayesian optimization can use large-data approximations without hiding what the approximation discarded |
| `fs-uq::adaptive_mlmc` | Adaptive MLMC admission checks, rate recovery, CVaR and anytime probability tests | Uncertainty estimates reject invalid budgets and tolerances before they produce misleading evidence |
| `fs-flutter-e2e` | A partitioned added-mass flutter campaign whose Lyapunov certificate is checked against the actual eigenvalue boundary and a separate numerical-abscissa implementation | Stability certificates are not accepted merely because two equivalent code paths agree; the README now distinguishes a true independent boundary check from an implementation cross-check |
| `fs-neuroshape-e2e` | A neural-SDF campaign with a closed interval boundary frame, no-tunnel Lipschitz tracing, Morse evidence, and visualization checks | Neural geometry gets a bounded-region proof instead of relying on sampled ring points that could leave escape gaps |
| `fs-ornith` | A smoke-tier ornithoid aircraft campaign: parameter decoding, BEM/VPM screening, e-raced candidate elimination, LBM refinement, SOS trim certificates, conformal L/D bands, and Pareto atlas rows | A flagship-style study can carry certificates, surrogate uncertainty, model-form honesty, replay, and graceful degradation through one pipeline |
| `fs-topopt::marquee` | A feature-gated CutFEM-octree topology lane with density-as-SDF geometry, DWR-guided refinement, zero-remesh iteration logs, and thickness/topology witnesses | Topology optimization can be exercised without making remeshing the first operation, while the feature gate keeps the proof state explicit |

The details matter. The sparse GP surface does not just say "approximate GP";
it records the ELBO slack that measures how much covariance mass was thrown
away. The flutter surface does not just say "stable"; it separates the
sufficient Lyapunov condition, the actual eigenvalue criterion, and the
implementation cross-check. The ornith surface does not claim high-fidelity
aerodynamics; it labels panel-vs-LBM agreement as model-form evidence and keeps
surrogate fallbacks inside conformal bands.

## Current Snapshot: What Is Newly Concrete

The current `main` snapshot is no longer just a crate atlas with isolated
building blocks. Several formerly aspirational paths now have concrete, tested
surfaces that show how the layers are meant to compose.

| Surface | What exists now | Why it matters |
|---------|-----------------|----------------|
| `fs-wasm::flagships` | Tier V browser-facing `run_ornithoid`, `run_vessel`, and `run_frame` pipelines | The browser surface now exercises real reduced flagship workflows instead of only leaf numerical kernels |
| `fs-mesh` v2/v3 | Conforming PLC recovery, non-convex facet triangulation, exact audits, measured 10-million-point perf lane, and a measured boundary-layer decision | Meshing claims now include recovery quality, determinism, and honest continuation notes for residual near-coplanar slivers |
| `fs-sparse::CsrCompact` | Compact u32 column-index CSR, sharded deterministic SpMV, tiled deterministic parallel assembly, and NUMA first-touch helpers | Sparse performance work is framed as a bandwidth and page-placement problem without sacrificing bitwise equality |
| `fs-rand` fast paths | Ziggurat normal fast path, bulk Philox fills, and dev-only stream statistics | Random performance can improve while the strict deterministic Box-Muller path remains the certification default |
| `fs-fft` N-D transforms | c2r inverse plus 2D/3D separable pencil transforms | Downstream physics and visualization crates get a broader deterministic transform surface before SIMD/roofline optimization lands |
| IO/risk/probe hardening | PLY face-list integer validation, non-finite CVaR/objective/probe rejection, and package provenance checks | Failure modes that used to look like ordinary inputs now fail closed before they contaminate evidence |

That pattern is important: a new feature is not considered mature simply
because the algorithm appears in code. It becomes useful when the contract names
its determinism class, tests pin representative behavior, the no-claim boundary
is explicit, and any performance claim is tied to a measured lane.

## Latest Implementation Deep Dives

Several recent implementation slices are useful to read because they show the
project's design rules becoming ordinary code rather than project vocabulary.

| Slice | What landed | Why it matters |
|-------|-------------|----------------|
| Deterministic integer powers | `fs-math` now owns pinned integer-power semantics and `xtask` checks for drifting `.powi` usage across dependent crates | Golden values no longer depend on libm or build-mode accidents for a common scalar operation |
| Declared run identity | `fs-exec` and `fs-rand` bind random streams to a declared run identity instead of pool history | Stochastic replay follows the logical study, not whichever worker happened to execute first |
| Caller-owned cancellation gates | `fs-race` and session pressure handling now require gates supplied by the caller/session owner | A race, pause, or memory-pressure response cannot manufacture private cancellation state that the owner cannot observe |
| Versioned solver snapshots | `fs-exec` solver state is wrapped in a versioned, self-authenticating envelope | Pause/resume/fork support gets a concrete artifact boundary instead of a raw struct dump |
| GEMM performance evidence | `fs-la` added packed f32 and mixed-precision paths, transposed/strided op-form GEMM, batched perf lanes, and roofline regression hooks | Dense-kernel speedups are tied to workload shape, denominator, and regression checks instead of prose claims |
| Risk and certificate hardening | CVaR uses fractional boundary weighting, adjoint certificates fail closed on vacuous evidence, and explain/DWR regressions stay executable | Evidence objects stop accepting plausible-but-empty proofs or biased tail estimates |
| Import-order correctness | `fs-io` accepts legal face-before-vertex PLY element order while still validating face payloads | The importer distinguishes format legality from malformed data, which is the right failure boundary for quarantined IO |

Together, these slices describe the current engineering center of gravity:
determinism is being pushed down into reusable primitives, cancellation is
treated as an ownership contract, evidence is allowed to weaken or refuse, and
performance work is expected to leave behind a measurable lane.

## Algorithms And Design Patterns Now Visible

The implemented workspace has enough surface area that recurring algorithms and
design patterns are visible across crates.

| Pattern | Where it appears | Design rule |
|---------|------------------|-------------|
| Fixed logical identity | Philox stream keys, tile IDs, solver snapshot envelopes, package roots | The identity of work must come from the problem and run record, not scheduler timing |
| Deterministic accumulation | COO/CSR assembly, sharded SpMV, GEMM tests, FFT goldens, evidence composition | Parallelism may change when work runs, but not the semantic order of the claim being checked |
| Fail-closed certification | Evidence colors, adjoint/DWR certificates, CVaR admission, PLY quarantine, package checking | Missing, non-finite, vacuous, or malformed evidence should refuse or downgrade, not pass as a stronger claim |
| Measurement-backed performance | Roofline axes, tune rows, perf lanes, CUSUM/regression gates, machine fingerprints | A performance statement needs a denominator, workload, hardware context, and acceptance band |
| Representation-preserving conversion | Region/chart records, SDF/mesh/NURBS/voxel/F-rep contracts, sheaf and topology witnesses | Geometry conversion is not just data translation; it must carry validity, topology, error, and no-claim context |
| Agent-readable governance | Beads, `CONTRACT.md`, DSR policy, `xtask` gates, changelog research notes | Automated work needs durable intent, exact acceptance criteria, and replayable evidence rather than conversational memory |

This is the useful part of FrankenSim as a source workspace today. The crates
are not merely a collection of numerical routines; they are converging on a
shared protocol for scientific software: state the claim, encode the assumptions,
run under explicit budgets, produce evidence with a color, and retain enough
lineage for another process to check the result later.

## Internal Data Model

FrankenSim repeatedly moves the same kinds of records across crate boundaries.
The implementation is easier to understand if you track those records rather
than memorizing every crate name.

| Record | Carried by | Used for |
|--------|------------|----------|
| Dimensions and units | `fs-qty`, `fs-scenario`, material cards, objectives, constraints | Rejecting unit drift before geometry or physics consumes an invalid scalar |
| Geometry representation | `fs-geom`, `fs-rep-*`, `fs-query`, `fs-topo`, `fs-mesh` | Preserving chart identity, conversion cost, validity, topology, and no-claim boundaries |
| Execution context | `fs-exec::Cx`, stream keys, cancellation gates, tile identities | Making cancellation, deterministic reduction, RNG replay, and scheduling explicit |
| Operator semantics | `fs-opdsl`, `fs-feec`, `fs-cutfem`, `fs-sparse`, `fs-solver` | Naming the algebraic object, residual norm, transpose action, and matrix-free/materialized boundary |
| Sensitivity path | `fs-ad`, `fs-adjoint`, VJP registries, gradient checks | Preventing a design update from using a gradient whose transpose or differentiability claim is missing |
| Evidence color | `fs-evidence`, `fs-package`, campaign crates | Distinguishing verified, validated, and estimated claims without laundering weaker evidence into stronger language |
| Provenance root | `fs-ledger`, `fs-package`, `fs-checker`, `fs-report` | Letting a later reader inspect the commit, lock state, artifacts, package root, and budget pie |

In practice, a run should not produce just `f64`. It should produce a value plus
the surrounding records that explain its units, representation, operator,
budget, solver status, evidence color, and provenance. The project is useful
because those records are normal Rust types instead of comments in a notebook.

## Execution Model

The runtime model is deliberately plain: work is split into tiles, tiles run
under an explicit context, and every stochastic or parallel decision is keyed by
logical identity rather than by worker timing.

```text
study / scenario / seed
        |
        v
Cx + stream keys + budgets
        |
        v
tile programs
        |
        +--> deterministic assembly / operator application
        +--> cancellable solver or simulator step
        +--> stochastic draw keyed by logical identity
        |
        v
reports, evidence colors, tune rows, ledger artifacts
```

This is why the code spends so much effort on things that look mundane:
canonical sort order, row-range sharding, fixed reduction trees, Philox stream
keys, idempotency keys, and content hashes. They make it possible to run a
parallel or stochastic computation and still answer: which logical tile did
this work, which seed did it use, what budget did it consume, and can the same
claim be replayed?

## Validation Strategy

FrankenSim uses the Gauntlet tiers from the project plan as a shared vocabulary
for proof strength. Not every crate is at every tier, but the README and
contracts should make the intended tier visible.

| Tier | What it proves | Examples in the workspace |
|------|----------------|---------------------------|
| G0 | Algebraic laws and local invariants | FFT round trips, sparse format equality, package root checks, quantity/unit laws |
| G1 | Manufactured solutions or convergence order | PDE and operator fixtures where closed-form or manufactured references exist |
| G2 | Canonical benchmarks | Chebyshev/Orr-Sommerfeld probes, Poiseuille-style flow checks, topology/optimization benchmark envelopes |
| G3 | Metamorphic behavior | Translation, scaling, relabeling, refinement, cross-format, and representation round-trip checks |
| G4 | Chaos and cancellation | Cancellation storms, budget exhaustion drills, ledger crash recovery, losing-branch drain tests |
| G5 | Determinism audits | Bitwise replay, fixed golden hashes, cross-thread sparse assembly, Philox stream identity, flagship smoke-stage hashes |

The strongest tests are usually paired. A deterministic golden catches drift,
but a golden alone can freeze a bug. A closed-form oracle catches meaning, but
an oracle alone may not exercise parallel replay. The codebase therefore tends
to pair fixtures: oracle plus determinism, dense reference plus sparse path,
feature gate plus no-claim boundary, or smoke campaign plus failure drill.

## Performance Model

Performance work in FrankenSim is framed as an evidence problem, not as a style
preference. A kernel is not "fast" because it uses the right algorithmic buzzword;
it is fast only when the repo records the machine, denominator, workload, and
acceptance band.

| Concern | Implementation direction |
|---------|--------------------------|
| Memory bandwidth | `fs-roofline` and `fs-substrate` measure machine axes and STREAM-like baselines before interpreting kernel throughput |
| Sparse SpMV | `CsrCompact` reduces index traffic, `spmv_sharded` balances contiguous row ranges by nnz, and NUMA first-touch places pages with the shard that will read them |
| Assembly | `Coo::assemble_parallel` tiles by row range and preserves duplicate accumulation order so thread count cannot change the numerical result |
| SIMD | Unsafe or architecture-specific code is kept behind registered capsules and safe facades; exploratory SME2 stays gated rather than becoming a default path |
| FFTs | The current N-D implementation is correctness-first and separable; higher-radix, SIMD, cache-blocked transposes, and executor-tiled pencils remain explicit follow-up work |
| Random generation | Bulk Philox and ziggurat normals are performance paths, while the strict path stays cross-ISA deterministic until the faster path earns the same proof |

The recurring tradeoff is intentional: first make the semantics deterministic
and testable, then optimize the hot path while proving it is still the same
computation. When the optimized path is not yet equally proven, it stays behind
a feature flag, fast-mode API, or no-claim boundary.

## How To Read A Crate Contract

Every `fs-*` crate has a `CONTRACT.md`. These files are not decorative. They are
the quickest way to tell whether a crate's README-level description is backed by
code and what the code refuses to claim.

| Section | Question to ask |
|---------|-----------------|
| Purpose and layer | Which architectural layer owns this behavior, and which layers may depend on it? |
| Public types and semantics | What values does this crate expose, and what does each value mean? |
| Invariants | What should remain true across refactors, thread counts, seeds, or representations? |
| Error model | Which failures are recoverable results, which are programmer errors, and which are explicit refusals? |
| Determinism class | Is replay bitwise, same-ISA, cross-ISA, statistical, or intentionally fast-mode-only? |
| Cancellation behavior | Where can long-running work stop, drain, and report structured cancellation? |
| Unsafe boundary | Is there unsafe code, and if so what invariant makes it valid? |
| Feature flags | Which parts are default, frontier, moonshot, or intentionally off the critical path? |
| Conformance tests | Which tests define the reimplementation contract? |
| No-claim boundaries | What would be dishonest to infer from the current implementation? |

If you only have time to read one file inside a crate, read the contract first.
Then read the tests named by the contract. The implementation usually makes much
more sense after those two files explain the claim surface.

## Concrete First Dives

These are good first code paths for understanding how the architecture works in
practice.

| If you want to understand... | Read these files first | What to look for |
|------------------------------|------------------------|------------------|
| Evidence packages | `crates/fs-package/src/lib.rs`, `crates/fs-checker/src/lib.rs`, their tests | How claims, provenance, Merkle roots, signatures, and budget pies become independently checkable |
| Deterministic sparse work | `crates/fs-sparse/src/lib.rs`, `crates/fs-sparse/src/perf.rs`, `crates/fs-sparse/tests/conformance.rs` | Sorted assembly, duplicate accumulation order, compact CSR, sharded SpMV, and cross-thread equality |
| Mesh honesty | `crates/fs-mesh/CONTRACT.md`, `crates/fs-mesh/src/recovery.rs`, `crates/fs-mesh/tests/conformance.rs` | Exact audits, recovery counters, non-convex facet handling, and measured boundary-layer decisions |
| Browser flagships | `crates/fs-wasm/src/flagships.rs`, `crates/fs-wasm/src/lib.rs` | How reduced aircraft/vessel/frame campaigns compose existing crates into trap-free browser exports |
| Deterministic transforms | `crates/fs-fft/src/lib.rs`, `crates/fs-fft/CONTRACT.md` | Fixed twiddle generation, oracle tests, c2r inverse, N-D pencil decomposition, and no performance overclaiming |
| Random replay | `crates/fs-rand/src/lib.rs`, `crates/fs-rand/CONTRACT.md`, `crates/fs-rand/tests/ziggurat.rs` | Logical stream identity, random access, bulk fills, fast-mode ziggurat, and strict-mode boundaries |

These paths are intentionally small enough to read in one sitting. Together
they show the project style: data structures first, contracts second, tests that
define the claim, and implementation details that avoid hidden global state.

## How The Algorithms Compose

FrankenSim is useful because its algorithms are not isolated demos. A typical
workflow crosses several crates, and each crossing has a contract.

```text
geometry or design variables
        |
        v
representation-specific chart / SDF / mesh / field
        |
        v
operator assembly or matrix-free operator
        |
        v
solver, verifier, or simulator
        |
        v
adjoint / sensitivity / uncertainty wrapper
        |
        v
optimizer, race, or value-of-information planner
        |
        v
evidence color + ledger/package/report artifact
```

| Crossing | Implementation shape | Failure it prevents |
|----------|----------------------|---------------------|
| Design to geometry | `fs-xform`, representation crates, interval enclosures, topology checks | A parameter update silently leaves the geometry validity domain |
| Geometry to physics | CutFEM, FEEC, BEM, LBM, material/scenario records | A solver receives an untyped shape with unknown boundary meaning |
| Physics to solver | Linear-operator traits, Krylov reports, p-multigrid, block preconditioners, deterministic sparse formats | A residual or matrix format looks fine locally but cannot be replayed or checked elsewhere |
| Solver to sensitivity | AD, IFT adjoints, VJP registries, transposed solves, gradient verification | An optimizer follows a gradient that differentiated the wrong computation |
| Sensitivity to decision | L-BFGS/TR-Newton, CMA-ES, NSGA-II, races, MLMC, CVaR, DRO | A campaign optimizes a noisy proxy without recording stop reasons or uncertainty |
| Decision to artifact | Evidence colors, content roots, package checks, crosswalk records, reports | A reviewer sees only a final plot and cannot inspect the chain of assumptions |

This is why the workspace has many small crates instead of one monolith. The
crate boundary is where the system can ask a useful question: did this layer
preserve dimensions, topology, determinism, residual meaning, evidence color,
or provenance?

## Design Principles For Agentic Simulation Work

FrankenSim is written with automated agents as a real workload. Agents are good
at running many narrow checks, but they are also prone to losing context,
over-trusting a local success, or optimizing for the easiest measurable target.
The architecture tries to make the correct behavior the path of least
resistance.

| Principle | What it means in practice |
|-----------|---------------------------|
| Make claims local | A crate should state exactly what it can prove, what it estimates, and what it refuses to claim in `CONTRACT.md` |
| Prefer typed refusal to heroic guessing | Invalid tolerances, missing VJPs, unsupported representation conversions, and incomplete evidence should fail loudly or weaken the evidence color |
| Keep replay independent of scheduling | Stream keys, deterministic reductions, sorted assembly, and cancellation contracts make parallel or stochastic work auditable |
| Separate model-form honesty from numerical accuracy | A low-Re LBM channel, inviscid panel method, surrogate, and SOS proof can coexist only when each says what kind of evidence it provides |
| Use smoke tiers as integration pressure | Small e2e campaigns are valuable because they force geometry, solvers, uncertainty, optimization, and reports to share one evidence language |
| Treat performance as an evidence object | Roofline numbers, tune rows, machine fingerprints, and DSR logs are the claim; architectural intent is not enough |

The result is a codebase where an agent can make progress without inventing a
new project philosophy for every task. The same rules keep appearing: check the
contract, run the conformance test, preserve the proof state, and do not upgrade
a claim just because one path passed.

## How To Evaluate A FrankenSim Result

When a result looks interesting, inspect it in this order:

1. Find the crate contract and read the no-claim boundaries.
2. Identify the evidence color or certificate type attached to the output.
3. Check whether the computation is deterministic, stochastic with a valid
   stopping rule, or model-form estimated.
4. Look for an independent oracle: dense reference, closed-form endpoint,
   symmetry law, LP enumeration, interval enclosure, conservation law, or
   cross-ISA golden.
5. Check the stop reason, residual norm, budget, and validity domain.
6. Confirm that any generated report or package contains enough provenance to
   rerun or independently check the claim.
7. If performance is part of the claim, require a roofline/tune/DSR artifact
   rather than accepting timing prose.

This checklist is intentionally mundane. The project is trying to make advanced
simulation work feel less like a notebook archaeology exercise and more like a
chain of inspectable, typed claims.

## Failure Modes The Architecture Is Built To Catch

Many simulation bugs do not look like panics. They look like plausible numbers
that lost their assumptions. FrankenSim tries to turn those cases into typed
refusals, weaker evidence colors, or failing conformance tests.

| Failure mode | How it can happen | Mechanism in FrankenSim |
|--------------|-------------------|-------------------------|
| Unit drift | A pressure, length, time step, or material constant enters as an unlabelled `f64` | `fs-qty`, scenario contracts, and objective/constraint records keep dimensions explicit |
| Topology drift | A representation conversion closes a hole, flips orientation, or hides a thin feature | Chart contracts, topology checks, conversion records, and no-claim rows make representation changes auditable |
| False solver confidence | A residual is small in the wrong norm or under the wrong boundary rows | Solver reports, operator contracts, conformance batteries, and explicit stall diagnoses name what was solved |
| Silent gradient loss | A non-differentiable operation returns zero cotangent or a missing transpose is ignored | VJP registries, IFT adjoints, gradient certificates, and merge gates refuse missing derivative paths |
| Stochastic overconfidence | A campaign keeps sampling after seeing lucky early results and treats the final mean as fixed-sample evidence | E-processes, confidence sequences, e-racing, and value-of-information logic keep stopping rules visible |
| Platform-dependent replay | Libm, thread order, or unordered assembly changes a golden value | Deterministic math, Philox stream keys, sorted assembly, fixed reductions, and content roots pin the replay path |
| Evidence laundering | A validated model, estimated surrogate, and verified invariant are collapsed into one flattering status | Evidence colors compose conservatively; the weakest required ingredient limits the claim |
| Performance folklore | A kernel is called fast because the algorithm sounds fast | Roofline records, machine fingerprints, tune rows, DSR lanes, and no-claim boundaries separate measured facts from ambitions |

This is the practical design principle behind the project: failures should be
local, named, and inspectable. A user should be able to tell whether a result is
blocked by geometry validity, algebraic convergence, model-form uncertainty,
budget exhaustion, missing derivative support, or simply an unmeasured
performance claim.

## How A Claim Becomes Auditable

The audit path is designed to be mechanical:

1. A crate declares what it can and cannot claim in `CONTRACT.md`.
2. A computation returns values with residuals, stop reasons, evidence colors,
   validity domains, or budget records instead of a naked scalar.
3. The ledger records artifacts, events, metrics, tune rows, machine fingerprints,
   and lineage.
4. A report or package bundles the result with its provenance and budget pie.
5. A checker can later validate package structure, content roots, signature
   status, and claim completeness without needing to rerun the full solver;
   authenticated release admission requires a caller-supplied verifier.

That path is useful for research code, but it is also useful for day-to-day
engineering. If a golden hash changes, if a solver starts returning a weaker
certificate, if a new representation loses topology, or if a stochastic estimate
stops being decision-grade, the failure has a place to surface. The project tries
to make those failures boring and local instead of late and forensic.

## Evidence Package Lifecycle

The evidence-package crates are small, but they define an important workflow:

1. A solver, verifier, optimizer, or study runner creates claims with evidence colors: verified, validated, or estimated.
2. `fs-package` bundles those claims with provenance and a deterministic content address.
3. `fs-checker` re-verifies the bundle's structural completeness, expected root, signature status, and budget pie without pulling in the solver stack. Its explicit release API additionally requires a caller-supplied authenticator, certificate-class falsifiers, and matching validation anchors.
4. `fs-crosswalk` gives the package a standards-facing vocabulary so the same artifact can be discussed as validation metrics, domain of validation, provenance, certified bounds, and content integrity.
5. `fs-ledger` can retain artifacts, operations, events, tune rows, roofline measurements, and lineage so a later report can connect package claims back to the actual run history.

This is why the package/checker/crosswalk layer is useful even before FrankenSim has a polished application: it turns internal numerical claims into an inspectable artifact boundary.

## Reading Paths

Different readers should start in different places.

| Goal | Start here |
|------|------------|
| Understand the architecture | `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md`, then this README's architecture and crate-map sections |
| Check repository policy | `AGENTS.md`, `docs/CI_GATES.md`, `xtask`, and `constellation.lock` |
| Study numerical kernels | `crates/fs-math`, `crates/fs-la`, `crates/fs-sparse`, `crates/fs-fft`, `crates/fs-ivl`, `crates/fs-cheb` |
| Study geometry | `crates/fs-geom`, `crates/fs-rep-*`, `crates/fs-query`, `crates/fs-topo`, `crates/fs-mesh` |
| Study physics and solvers | `crates/fs-feec`, `crates/fs-cutfem`, `crates/fs-solid`, `crates/fs-material`, `crates/fs-solver` |
| Study gradients and optimization | `crates/fs-ad`, `crates/fs-adjoint`, `crates/fs-ascent`, `crates/fs-dfo`, `crates/fs-constraint` |
| Study audit and packaging | `crates/fs-evidence`, `crates/fs-ledger`, `crates/fs-package`, `crates/fs-checker`, `crates/fs-crosswalk` |

## Contracts and Verification

The tree currently has 127 `CONTRACT.md` files for 127 `fs-*` crate directories.
The contract count is meant to be checkable, not aspirational.

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
|-- crates/                            # 127 fs-* crates; selected entries shown below
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

The complete crate list is the `members` array in [Cargo.toml](Cargo.toml). Use `find crates -mindepth 2 -maxdepth 2 -name Cargo.toml` for a filesystem inventory.

## How FrankenSim Compares

| Capability | FrankenSim today | Mature solver stacks | General scientific Python |
|------------|------------------|----------------------|---------------------------|
| End-user simulation app | Not yet | Yes, depending on stack | Usually assembled by user |
| Library substrate | Yes, source workspace | Yes | Yes |
| Determinism and cancellation contracts | Explicit crate contracts | Varies by project | Varies by package |
| Evidence-carrying values | Implemented core wrappers | Usually external/ad hoc | Usually external/ad hoc |
| Geometry representation contracts | Implemented SDF/mesh/NURBS/voxel/F-rep chart layer | Often specialized | Usually library-specific |
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

From a fresh checkout, the supported way to create and verify that layout is:

```bash
cargo run --manifest-path tools/bootstrap/Cargo.toml
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

Every `fs-*` crate must ship a `CONTRACT.md` with the required sections. If this check fails, either a new crate was added without a contract or an existing contract lost one of the required sections. Use an existing crate contract as the local template and keep the no-claim boundaries honest.

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
| Neural representations | Not implemented as a first-class representation crate in the current workspace |
| Randomized NLA golden sentinel | Resolved: `rand_nla_golden_hash` is deliberately recorded (`0xeef1_0550_7daf_c0d5`) and verified identical on arm64 and x86-64 in both debug and release, after fixing a build-mode-dependent `powi` fixture; the workspace-wide `powi` sweep is tracked in bead `frankensim-powi-build-mode-determinism-4xnt` |
| Ascent golden sentinel | `fs-ascent` implements the optimizer stack, but its trajectory golden hash lane still needs the recorded hash to be frozen in the test source |
| Long-running stability fixtures | Some structural stability and snap-through tests are active proof lanes and may need targeted runtime/threshold work rather than being treated as cheap smoke tests |
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

Yes. Existing `CONTRACT.md` files are part of how the workspace communicates implemented behavior. If code and contract disagree, that is a bug to resolve, not a documentation detail to ignore.

### Should I use GitHub Actions?

No. Use DSR for this repository unless the repository policy changes. GitHub Actions may be absent, stale, or non-authoritative.

## About Contributions

Please don't take this the wrong way, but I do not accept outside contributions for any of my projects. I simply don't have the mental bandwidth to review anything, and it's my name on the thing, so I'm responsible for any problems it causes; thus, the risk-reward is highly asymmetric from my perspective. I'd also have to worry about other "stakeholders," which seems unwise for tools I mostly make for myself for free. Feel free to submit issues, and even PRs if you want to illustrate a proposed fix, but know I won't merge them directly. Instead, I'll have Codex or Codex review submissions via `gh` and independently decide whether and how to address them. Bug reports in particular are welcome. Sorry if this offends, but I want to avoid wasted time and hurt feelings. I understand this isn't in sync with the prevailing open-source ethos that seeks community contributions, but it's the only way I can move at this velocity and keep my sanity.

## License

FrankenSim is licensed under the MIT License with the repository's Anthropic/OpenAI rider. See [LICENSE](LICENSE) for the exact terms.
