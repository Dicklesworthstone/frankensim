<div align="center">
  <img src="frankensim_illustration.webp" alt="FrankenSim - Rust continuum for geometry, physics, optimization, and rendering" width="900">
</div>

<h1 align="center">FrankenSim</h1>

<p align="center">
  <strong>A plan-first Rust continuum for certified geometry, physics, optimization, and rendering.</strong>
</p>

<div align="center">

![Status](https://img.shields.io/badge/status-design%20plan-blue)
![Language](https://img.shields.io/badge/language-Rust-orange)
![Runtime deps](https://img.shields.io/badge/runtime%20deps-Franken--only-lightgrey)

</div>

FrankenSim is designed as one continuous pipeline from geometry to physics to
optimization to rendering, with derivatives, error bounds, budgets, provenance,
and cancellation kept inside the values that move through the system.

```bash
git clone https://github.com/Dicklesworthstone/frankensim.git
cd frankensim
sed -n '1,180p' COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md
```

This repository currently contains the public design plan and agent operating
rules. It is not yet an installable crate or executable.

## TL;DR

**The problem:** Shape optimization against real physics is usually split across
CAD kernels, meshers, FEM/CFD solvers, optimizers, plotting tools, and notebooks.
Derivatives disappear at those boundaries. Error budgets are not composed.
Provenance is fragile. Cancellation is usually a process kill.

**The solution:** FrankenSim is designed as one typed Rust continuum where a
geometry, field, mesh, solver result, or Pareto point can carry its derivative
hooks, certified error bounds, budget state, provenance hash, and cancellation
scope through the whole pipeline.

### Why FrankenSim?

| Capability | What it is meant to provide |
| --- | --- |
| Certified representation routing | Move between SDFs, F-reps, NURBS, meshes, voxels, point clouds, and neural implicits with explicit cost and error receipts. |
| Structure-preserving physics | FEEC, CutFEM-on-SDF, variational integrators, port-Hamiltonian coupling, and matrix-free solvers. |
| Optimization with evidence | Adjoint-native gradients, CMA-ES and BO paths, topology optimization, e-process racing, and certificate-aware surrogate use. |
| Deterministic cancellation | asupersync scopes, tile-level cancellation checks, fixed reduction trees, logical RNG streams, and replayable solver states. |
| Design ledger | FrankenSQLite-backed artifacts, ops, metrics, tuning data, lineage, time travel, and `explain()` queries. |

## Quick Example

The current repository is a plan, so this is the intended FrankenScript shape,
not runnable code yet:

```lisp
(study "frame-seismic-cvar-v9"
  (seed 0xF00D0002) (versions (constellation :lock "2026-07"))
  (capability :cores 96 :mem 384GiB :wall 36h :ops (flux.* ascent.* uq.*))
  (budget (qoi "P(drift>2e-2)" :rel-error 0.15 :confidence 0.95))
  (let ground (topo.ground-structure (grid 8 x 5 x 24m)
                                      :knn 14
                                      :rules "AISC-cat.json"))
  (let layout (ascent.solve-lp (min (member-volume ground)) :method pdhg))
  (let frame  (topo.size layout :method tr-newton-krylov))
  (let resp   (flux.fiber-frame frame site :integrator variational))
  (let frag   (uq.probability (exceeds (peak-drift resp) 2e-2)
                :stop (e-process :alpha 0.05)))
  (ascent.optimize (min (mass frame)) :over (sections frame)
    :subject-to ((cvar frag :beta 0.9 :le 0.02))
    :method augmented-lagrangian
    :emit (frame frag report ledger)))
```

The syntax is provisional. The invariant is that seeds, versions, capabilities,
budgets, physical models, optimization state, and output artifacts are all
first-class data.

## Design Principles

| Principle | Meaning |
| --- | --- |
| Pure memory-safe Rust | Runtime dependencies are limited to `std` plus the Franken constellation: asupersync, FrankenSQLite, FrankenNumpy, FrankenTorch, FrankenScipy, FrankenPandas, and FrankenNetworkx. |
| Determinism by contract | Deterministic mode uses logical RNG streams, fixed reduction trees, and deterministic tie-breaking so studies can be replayed. |
| Differentiable or certifiable | Operators expose adjoints, interval/Taylor enclosures, convergence evidence, or explicit no-claim boundaries. |
| Budgets first | Accuracy, wall time, memory, and capability grants are values in the IR and ledgers. |
| One data model | Complexes and cochains connect geometry and physics instead of leaving every solver to invent its own representation. |
| Provenance complete | Artifacts are content-addressed and every operation lands in the Design Ledger. |

## Architecture

```text
+-----------------------------------------------------------------------+
| L6 HELM      FrankenScript IR, sessions, budgets, ledger, planner     |
+-----------------------------------------------------------------------+
| L5 LUMEN     spectral rendering, direct chart tracing, visualization  |
+-----------------------------------------------------------------------+
| L4 ASCENT    adjoints, topology optimization, CMA-ES, BO, SOS, UQ     |
+-----------------------------------------------------------------------+
| L3 FLUX      FEEC, CutFEM, LBM, BEM/FMM, solvers, coupling, adjoints   |
+-----------------------------------------------------------------------+
| L2 MORPH     Regions, charts, Rep Router, meshing, validity proofs    |
+-----------------------------------------------------------------------+
| L1 BEDROCK   dense/sparse/FFT math, intervals, AD, RNG, GA            |
+-----------------------------------------------------------------------+
| L0 SUBSTRATE asupersync execution, arenas, SIMD, NUMA, determinism     |
+-----------------------------------------------------------------------+
| Franken constellation: asupersync, SQLite, Numpy, Torch, Scipy,        |
| Pandas, Networkx                                                       |
+-----------------------------------------------------------------------+
```

The workspace is intended to be a flat set of `fs-*` crates with a strict
acyclic dependency order. Each crate should have a `CONTRACT.md`, executable
conformance tests, clear invariants, an error model, and no-claim boundaries.

## Planned Crates

| Layer | Crates |
| --- | --- |
| L0/L1 | `fs-substrate`, `fs-simd`, `fs-alloc`, `fs-exec`, `fs-la`, `fs-sparse`, `fs-fft`, `fs-ivl`, `fs-cheb`, `fs-ad`, `fs-rand`, `fs-ga` |
| L2 | `fs-geom`, `fs-rep-sdf`, `fs-rep-frep`, `fs-rep-nurbs`, `fs-rep-mesh`, `fs-rep-voxel`, `fs-rep-neural`, `fs-xform`, `fs-mesh` |
| L3 | `fs-feec`, `fs-cutfem`, `fs-iga`, `fs-solid`, `fs-lbm`, `fs-bem`, `fs-fmm`, `fs-vpm`, `fs-time`, `fs-couple`, `fs-adjoint`, `fs-uq` |
| L4 | `fs-eproc`, `fs-opt`, `fs-topo`, `fs-sos`, `fs-surrogate` |
| L5/L6 | `fs-render`, `fs-viz`, `fs-img`, `fs-ir`, `fs-session`, `fs-ledger`, `fs-plan`, `fs-report` |

## Comparison

| Question | FrankenSim target | Conventional CAD + FEM/CFD stack | Optimizer around black-box solver |
| --- | --- | --- | --- |
| Do derivatives cross geometry, mesh, and solver boundaries? | Yes, by design | Usually no | Usually no |
| Are error sources composed? | Yes, through Error Ledger models | Rarely | Rarely |
| Can bad candidates be cancelled mid-solve? | Yes, through asupersync scopes and tile checkpoints | Usually no | Often only by killing a process |
| Is provenance queryable? | Yes, through content-addressed artifacts and ops | Often manual | Often manual |
| Does it depend on native BLAS/Fortran/C++ kernels? | No production dependency | Common | Common |
| Can it render the same chart used by the solver? | Yes, through LUMEN and direct chart backends | Usually requires export | Usually out of scope |

## Installation

There is no installable release yet.

### Read the Plan

```bash
git clone https://github.com/Dicklesworthstone/frankensim.git
cd frankensim
less COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md
```

### Watch for the First Workspace

When P0 exists, the expected source workflow will be:

```bash
git clone https://github.com/Dicklesworthstone/frankensim.git
cd frankensim
cargo test --workspace
```

### Package Managers

No Homebrew, crates.io, or binary packages exist yet.

## Current Repository Contents

```text
frankensim/
|-- AGENTS.md
|-- COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md
|-- README.md
|-- frankensim_illustration.webp
`-- .gitignore
```

## Roadmap

| Phase | Scope | Exit criterion |
| --- | --- | --- |
| P0 Bedrock | `fs-substrate`, `fs-exec`, `fs-alloc`, `fs-la`, `fs-sparse`, `fs-fft`, `fs-ivl`, `fs-rand`, ledger v0 | G0 and G4 green; GEMM, SpMV, and FFT within target bands; deterministic mode bit-stable |
| P1 Geometry + eyes | Regions, SDF/F-rep/mesh charts, Rep Router v1, meshing, preview tracer | Certified chart round trips and watertightness checks |
| P2 Elasticity + first optimization | FEEC elasticity, CutFEM-on-SDF, matrix-free multigrid, adjoints, SIMP | Topology optimization on a raw SDF with composed error certificate |
| P3 Fluids I | Sparse/free-surface LBM, scaling assistant, thermal and non-Newtonian paths | Cavity, Taylor-Green, and cylinder benchmarks green |
| P4 Structures at scale | IGA shells, fiber beams, ground-structure PDHG, MLMC, e-stop | Seismic frame flagship with anytime-valid fragility |
| P5 Aero stack | BEM/FMM, vortex particles, coupling, SE(3), Koopman surrogates | Ornithoid Pareto run with e-raced generations |
| P6 Certificates and planning | SOS/Lasserre, sheaf certificates, conformal e-prediction, planner, diff-rendering | Moonshot features pass certifier tests or stay flagged off |

## The Gauntlet

FrankenSim treats technical claims as obligations. The planned verification
program has six tiers:

| Tier | Purpose |
| --- | --- |
| G0 | Property tests and algebraic laws |
| G1 | Manufactured solutions and convergence-order checks |
| G2 | Canonical benchmarks |
| G3 | Metamorphic tests |
| G4 | Chaos, cancellation storms, leak checks, deadlock checks |
| G5 | Determinism audits |

Features marked `[F]` or `[M]` in the plan must stay out of default paths until
the relevant Gauntlet evidence exists.

## Performance Targets

The plan is written around Apple Silicon and many-core x86. Targets are meant to
be measured and failed, not treated as marketing copy.

| Kernel family | Example target |
| --- | --- |
| LBM sparse D3Q19 | 1.0 GLUP/s class on Apple M-series; 0.6 GLUP/s class on 96-core Threadripper |
| GEMM f64 | 75% of measured peak for the selected SIMD tier |
| SpMV / SELL-C-sigma | 85% of STREAM-class bandwidth |
| Matrix-free FEEC apply | 30% of peak FLOPs for p=4 sum-factorized paths |
| Sphere-traced SDF rays | 80 to 120 Mray/s class, depending on machine |

Every real target should eventually live beside a machine fingerprint,
benchmark command, acceptance band, and ledger record.

## Command Reference

There is no FrankenSim CLI yet. Current useful commands are repository commands:

```bash
# Read the plan
sed -n '1,180p' COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md

# Inspect the agent contract
sed -n '1,220p' AGENTS.md

# Once the Rust workspace exists
cargo test --workspace
```

Planned command surfaces will likely come from `fs-ir`, `fs-session`,
`fs-ledger`, and `fs-report`.

## Configuration

No runtime configuration file exists yet. The intended system configuration is
explicit in the FrankenScript IR:

```lisp
(study "example"
  (seed 0x5EED0001)
  (versions (constellation :lock "2026-07"))
  (capability :cores 16 :mem 64GiB :wall 2h)
  (budget (qoi-rel-error 2e-2))
  ...)
```

The Five Explicits are required for every real study:

- units
- seeds
- budgets
- versions
- capabilities

## Limitations

| Area | Current state |
| --- | --- |
| Code | No Rust workspace has been created yet. |
| Installation | No release, package, or installer exists. |
| Claims | Performance, correctness, and certificate claims are design targets until implemented and tested. |
| Dependencies | The production dependency policy is intentionally narrow and will make some implementation work harder. |
| Moonshot items | Sheaf certificates, e-raced optimization, and self-optimizing planners must stay behind feature flags until validated. |

## Troubleshooting

### `cargo test` says there is no `Cargo.toml`

That is expected right now. The repository is still in the design-plan stage.

### The README mentions packages that do not exist

The `fs-*` crate list is the intended workspace map from the plan. It is not a
published Cargo workspace yet.

### There is no installer

Correct. Clone the repo and read the plan for now.

### A claimed feature looks too ambitious

Check the tag in the plan. `[S]` is solid engineering, `[F]` is frontier work,
and `[M]` is moonshot work that must be feature-gated until the Gauntlet validates
it.

### Can I rely on benchmark numbers now?

No. The numbers in the plan are targets. Real claims require benchmark artifacts,
machine fingerprints, and acceptance bands.

## FAQ

### Is this usable today?

No. Today it is a public plan and coordination repo.

### Why build this instead of binding to existing CAD and solver libraries?

The central design goal is to keep derivatives, error bounds, budgets,
provenance, and cancellation together across every layer. Wrapping a pile of
separate tools would preserve the boundaries that the project is trying to
remove.

### Why Rust?

Rust gives the project ownership, lifetimes, const generics, zero-cost
abstractions, fearless concurrency, and a practical path to high-performance
safe code with narrow audited unsafe leaves.

### Why avoid BLAS, LAPACK, C, and C++ in production paths?

FrankenSim needs kernels shaped around its own layouts, tile scheduler,
determinism model, and cancellation protocol. External native kernels can still
serve as development or conformance references when they are isolated.

### What should be built first?

P0: `fs-substrate`, `fs-exec`, `fs-alloc`, `fs-la`, `fs-sparse`, `fs-fft`,
`fs-ivl`, `fs-rand`, and ledger v0.

### Is this open to outside contributions?

Bug reports are welcome. The contribution policy below is intentionally strict.

## About Contributions

*About Contributions:* Please don't take this the wrong way, but I do not accept outside contributions for any of my projects. I simply don't have the mental bandwidth to review anything, and it's my name on the thing, so I'm responsible for any problems it causes; thus, the risk-reward is highly asymmetric from my perspective. I'd also have to worry about other "stakeholders," which seems unwise for tools I mostly make for myself for free. Feel free to submit issues, and even PRs if you want to illustrate a proposed fix, but know I won't merge them directly. Instead, I'll have Codex or Codex review submissions via `gh` and independently decide whether and how to address them. Bug reports in particular are welcome. Sorry if this offends, but I want to avoid wasted time and hurt feelings. I understand this isn't in sync with the prevailing open-source ethos that seeks community contributions, but it's the only way I can move at this velocity and keep my sanity.

## License

No license file is present yet. Treat this repository as public planning
material until a license is added.
