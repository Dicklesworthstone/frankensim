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
![License](https://img.shields.io/badge/license-MIT%20%2B%20OpenAI%2FAnthropic%20Rider-blue)

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

## What Exists Now

This repo currently publishes the design, operating rules, and public project
shell. That matters because FrankenSim is large enough that the contracts have
to come before the crates.

| Asset | Purpose |
| --- | --- |
| `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md` | Full technical plan: principles, architecture, crate atlas, Gauntlet, roadmap, risks, and flagship studies. |
| `AGENTS.md` | Rules for coding agents: no destructive actions, no branches or worktrees, P0-first scope, Rust dependency policy, contracts, and verification expectations. |
| `README.md` | Public project overview and onboarding path. |
| `frankensim_illustration.webp` | Compressed overview image for the public README. |
| `.gitignore` | Rust-oriented ignore rules, with the large source PNG excluded from git. |

The next meaningful milestone is P0: creating the initial Rust workspace and
landing the substrate, execution, allocator, numerical, interval, random, and
ledger foundations.

## Why This Is Useful

FrankenSim is aimed at design loops where a shape must be optimized against
physics rather than merely simulated once. Those loops are expensive today
because each tool boundary loses information.

| User question | Conventional workflow | FrankenSim target |
| --- | --- | --- |
| "If I move this control point, how does drag change?" | CAD changes, mesher changes, solver runs, optimizer estimates with finite differences. | Geometry parameters expose Jacobian actions, solvers expose adjoints, and the optimizer gets a gradient with provenance. |
| "How much should I trust this result?" | Tolerances are scattered across CAD, mesh, solver, and scripts. | The Error Ledger attributes uncertainty to geometry, discretization, algebraic residual, surrogate error, and statistical noise. |
| "Can we stop evaluating this bad candidate?" | The job usually runs to completion or gets killed from outside. | Candidate evaluations live in cancellation scopes, and statistical tests can stop losing branches mid-solve. |
| "Where did this Pareto point come from?" | Reconstructing the answer means reading old scripts and filenames. | `explain(artifact)` walks the content-addressed lineage graph. |
| "Can a mesh-free design loop still get FEM-grade evidence?" | Often no; meshing becomes the bottleneck. | CutFEM-on-SDF is a first-class physics path, with body-fitted meshes reserved for final verification or export. |

## How the Pipeline Works

A study is intended to move through the system as typed data rather than as a
pile of side effects:

```text
FrankenScript study
    |
    v
Static admission
    units, budgets, versions, capabilities, chart support
    |
    v
HELM task DAG
    budgeted nodes, idempotency keys, ledger scope
    |
    v
SUBSTRATE execution
    asupersync scopes, tile kernels, arenas, deterministic reductions
    |
    v
MORPH geometry
    Region charts, representation routing, validity certificates
    |
    v
FLUX physics
    complexes, cochains, solvers, residuals, adjoints, error estimates
    |
    v
ASCENT optimization
    gradients, CMA-ES, BO, topology updates, e-process stopping
    |
    v
LUMEN and reports
    renders, scientific visualization, Pareto tables, lab notebooks
    |
    v
Design Ledger
    artifacts, ops, edges, metrics, tune rows, events
```

The central rule is simple: values carry the information needed by downstream
layers. A field should know its units, budget slice, provenance, derivative
support, error bound, and cancellation context. A chart conversion should return
a receipt, not merely a mesh. A solver result should carry enough state to
resume, audit, and differentiate.

## Core Data Model

The plan organizes geometry and physics around a small set of recurring nouns.

| Concept | Role |
| --- | --- |
| `Region` | Abstract shape or domain. It can be presented by multiple concrete charts. |
| `Chart` | Concrete representation of a `Region`: SDF, F-rep, NURBS, mesh, voxel, point cloud, lattice, or neural implicit. |
| `Certified<T>` | Value plus error bound, provenance hash, and the hooks needed by downstream layers. |
| `Complex` | Cell complex used as the common substrate for geometry and physics. |
| `Cochain` | Field value living on cells, edges, faces, or volumes. |
| `Budget` | Accuracy, time, memory, and capability allocation for an operation. |
| `Cx` | Execution context carrying cancellation, arena, RNG identity, budget slice, and ledger handle. |
| `LedgerHandle` | Append path for artifacts, ops, metrics, tune rows, and events. |

Two sketches from the plan capture the intended style:

```rust
pub trait Chart: Send + Sync {
    type Param;

    fn eval(&self, x: Point3, cx: &Cx) -> ChartSample;
    fn support(&self) -> Aabb;
    fn topology_hint(&self) -> BettiBounds;
}
```

```rust
pub trait TileKernel: Sync {
    type Out: Reduce;

    fn tiles(&self) -> TilePlan;
    fn run(&self, tile: TileId, cx: &Cx) -> ControlFlow<Cancelled, Self::Out>;
}
```

These are not final APIs. They describe the pressure the final APIs must bear:
parallel execution, cancellation, determinism, error accounting, and provenance.

## Algorithmic Spine

FrankenSim is intentionally not a wrapper over one solver family. It is a set of
compatible algorithmic choices that share data layouts, cancellation semantics,
and ledgered evidence.

| Area | Algorithms and design choices |
| --- | --- |
| Representation routing | Directed chart graph, Pareto shortest paths, composed error bounds, cost models, certificate-preferred conversions. |
| Geometry robustness | Exact predicates, generalized winding numbers, interval broad/narrow phases, dual contouring, NURBS distance by Bezier bounds plus interval Newton where feasible. |
| Certified arithmetic | Interval arithmetic, affine arithmetic, Taylor models, adaptive floating-point expansions, double-double and quad-double escalation paths. |
| Execution | Two-lane asupersync runtime, work-stealing tile pool, scope arenas, deterministic reductions, logical RNG streams, checkpointable solver states. |
| Dense math | BLIS-style GEMM, batched small dense kernels, TSQR, LOBPCG, Lanczos, mixed precision with refinement. |
| Sparse math | CSR, BSR, SELL-C-sigma, matrix-free operators, Chebyshev smoothers, smoothed aggregation AMG, p-multigrid. |
| Geometry to physics | FEEC on cell complexes, CutFEM directly on SDFs, IGA on spline charts, DWR adaptivity driven by the objective. |
| Fluids | Sparse tiled LBM, free-surface LBM, BEM plus FMM, vortex particles, pressure-robust incompressible formulations. |
| Structures | Hyperelasticity, rods, shells, fiber beam-columns, buckling eigenproblems, arclength continuation, IPC-style contact. |
| Optimization | L-BFGS, trust-region Newton-Krylov, augmented Lagrangian, CMA-ES, differential evolution, Bayesian optimization, PDHG layout optimization. |
| Stochastic decisions | QMC, MLMC, e-process confidence sequences, e-BH, conformal e-prediction, certify-or-escalate surrogate policy. |
| Rendering | Spectral path tracing, direct SDF/F-rep sphere tracing, NURBS tracing, volume rendering, line integral convolution, differentiable rendering. |

## Representation Router

The Representation Router is the planned mechanism that decides how to satisfy a
request such as "give FLUX a chart that can evaluate signed distance, curvature,
and boundary traces within this error budget."

```text
SDF grid ---- dual contouring ----> mesh
   |                                ^
   |                                |
   v                                |
F-rep CSG ---- sampling ----------> sparse voxel chart
   |
   v
NURBS refit ---- certificate -----> spline chart
```

Each edge in the graph has:

- a cost model
- an error model
- a certificate status
- a provenance record
- a deterministic replay path

The Router solves a multi-objective path problem over cost and composed error.
If no path satisfies the caller's budget, admission fails early with a structured
diagnosis instead of running a doomed simulation.

## Physics Model

The physics layer is built around complexes and cochains because that makes
structure preservation a default constraint rather than an afterthought.

| Physics concern | Planned treatment |
| --- | --- |
| Conservation laws | Discrete exterior derivative is combinatorial, so identities such as curl-grad and div-curl are exact on supported complexes. |
| Embedded geometry | CutFEM runs on SDF and sparse grid charts, avoiding body-fitted remeshing inside topology loops. |
| Thin structures | IGA and Kirchhoff-Love shells preserve spline geometry where that matters. |
| Fluids | LBM handles many-core tiled free-surface and non-Newtonian flows; BEM/FMM and vortex methods cover screening and exterior aerodynamics. |
| Coupling | Port-Hamiltonian interfaces exchange power-conjugate effort/flow pairs to avoid energy-creating coupling artifacts. |
| Adaptivity | Dual-weighted residual estimators spend resolution where the objective is sensitive. |

The system should be able to say not only "the answer is 12.4" but "the answer
is 12.4 with this geometric tolerance, this discretization contribution, this
solver residual, this surrogate band, and this statistical half-width."

## Optimization Model

ASCENT treats optimization problems as data:

```text
variables:
  shape parameters, density fields, section choices, manifold coordinates

constraints:
  physics equations, volume, stress, stability, manufacturability, topology

objectives:
  drag, compliance, fragility, mass, robustness, render-derived losses

evidence:
  gradients, KKT residuals, confidence sequences, certificates, lineage
```

Gradient-based optimization uses discrete adjoints and Sobolev/Riesz smoothing.
Derivative-free optimization uses CMA-ES, differential evolution, DIRECT, and
trust-region interpolation models where gradients are unavailable or dishonest.
Bayesian and multi-fidelity paths decide whether to trust a surrogate, escalate
to a higher-fidelity solve, or stop early with an anytime-valid decision.

## What Counts as a Valid Result

FrankenSim should make result quality visible. A successful run should be able to
answer these questions without manual archaeology:

| Question | Expected evidence |
| --- | --- |
| What exact input produced this artifact? | FrankenScript IR, seed, versions, capability grant, and ledger op row. |
| Which representation path was used? | Router path and conversion receipts. |
| What did the error budget buy? | Error Ledger attribution by geometry, discretization, solver, surrogate, and statistics. |
| Was the run deterministic? | Execution mode, reduction shape, RNG stream identity, and replay metadata. |
| Were candidates cancelled safely? | Scope tree events, drain/finalize records, arena accounting. |
| Why should a gradient be trusted? | Adjoint identity tests, finite-difference or dual checks, and residual tolerances. |
| Why should a stochastic decision be trusted? | e-process or confidence-sequence record with stopping rule. |

## Flagship Studies

The plan uses three end-to-end studies to keep the architecture honest.

| Study | What it exercises |
| --- | --- |
| Ornithoid multi-inlet aircraft | F-rep and spline geometry, BEM/FMM screening, vortex wakes, LBM refinement, Koopman surrogates, SOS stability certificates, multi-objective Pareto search. |
| Seismic-minimal building frame | FrankenNetworkx ground structures, PDHG layout optimization, nonlinear sizing, fiber-section dynamics, MLMC ground motion ensembles, anytime-valid fragility estimates. |
| Laminar-pour vessel | F-rep vessel geometry, level-set lip optimization, Chebyshev stability models, free-surface LBM validation, viscosity-band robustness, LUMEN render output. |

These are forcing functions, not demos. A feature that cannot help one of these
studies needs a strong reason to enter the roadmap.

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

## Development Discipline

FrankenSim is meant to be built by multiple coding agents without letting the
repository turn into a pile of incompatible experiments. The operating discipline
is part of the design:

- one branch: `main`
- one crate, one contract
- no production FFI shortcuts
- no unledgered performance claims
- no unbounded `unsafe`
- no solver or conversion claim without a conformance path
- no `[M]` feature in the default path without certifier evidence

Every crate should eventually answer four questions in its `CONTRACT.md`:

1. What invariants does this crate own?
2. What errors can it bound, estimate, or refuse to claim?
3. What determinism and cancellation behavior does it guarantee?
4. Which Gauntlet tiers prove those claims?

## Near-Term Build Order

The first useful code should make later claims cheap to prove.

| Step | Deliverable | Why it comes early |
| --- | --- | --- |
| 1 | Workspace skeleton and contracts | Agents need stable crate boundaries before parallel implementation. |
| 2 | `fs-substrate` | Hardware fingerprints, cache geometry, SIMD tier selection, and topology map feed every hot path. |
| 3 | `fs-exec` and `fs-alloc` | Tile execution, cancellation scopes, deterministic reductions, and scoped arenas are cross-cutting. |
| 4 | `fs-ivl` and `fs-rand` | Certified arithmetic and deterministic random streams underpin tests, geometry, UQ, and rendering. |
| 5 | `fs-la`, `fs-sparse`, `fs-fft` | These kernels are used by almost every higher layer and give early roofline evidence. |
| 6 | ledger v0 | Provenance, metrics, and replay evidence should exist before complex studies begin. |

The planned P0 exit is deliberately concrete: G0 and G4 green, core numerical
kernels inside target bands, and deterministic mode bit-stable.

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

FrankenSim is licensed under the MIT License with the OpenAI/Anthropic rider.
See [LICENSE](LICENSE) for the exact terms.
