# COMPREHENSIVE PLAN FOR FRANKENSIM

*A single, memory-safe Rust continuum for computational geometry, physics, optimization, and rendering — designed from a blank slate for Apple Silicon and many-core x86, built on the Franken constellation (asupersync, FrankenSQLite, FrankenNumpy, FrankenTorch, FrankenScipy, FrankenPandas, FrankenNetworkx), with zero other runtime dependencies.*

---

## 0. How to read this document

This is a design plan, not a survey. Every mechanism described here is chosen because it is load-bearing for the mission: **given a physics-based objective and constraints, synthesize the geometry that optimizes it — faster, more correctly, and more verifiably than any existing system, on commodity many-core CPUs, in pure safe Rust.**

Ambition is calibrated with three tags used throughout:

- **[S] Solid** — established mathematics and engineering; the work is implementation excellence, not research risk.
- **[F] Frontier** — published research from roughly the last decade that no mainstream system has productized; real engineering risk, enormous payoff.
- **[M] Moonshot** — novel synthesis proposed here for the first time (to my knowledge); prototyped behind feature flags, promoted only after the Gauntlet (§13) validates it.

The mix is deliberate: the spine of FrankenSim is [S], the leapfrog features are [F], and a handful of [M] bets are what make the system unlike anything else. Nothing tagged [M] sits on the critical path of the roadmap.

---

## 1. Thesis: why a blank slate wins

Every existing pipeline for "optimize a shape against physics" is an **archipelago**: OpenCASCADE or a B-rep kernel for geometry, gmsh or a proprietary mesher for discretization, an FEM/CFD code (FEniCS, MFEM, deal.II, OpenFOAM, SU2, COMSOL, Abaqus) for physics, SciPy/NLopt/Dakota for optimization, ParaView for looking at the wreckage. Each island is excellent. The water between them is where everything drowns:

1. **Derivatives don't cross boundaries.** The CAD kernel doesn't know the mesh's sensitivity to a control point; the solver's adjoint dies at the mesher; the optimizer sees a noisy black box and falls back to finite differences or pure evolution.
2. **Error bounds don't cross boundaries.** Geometry tolerance, meshing error, discretization error, solver tolerance, statistical noise — nobody composes them. You get a number with no pedigree.
3. **Provenance doesn't exist.** Six weeks into a design study, nobody can reconstruct which mesh, which solver settings, and which random seed produced the Pareto point the client liked.
4. **Cancellation doesn't exist.** An optimizer that discovers a candidate is hopeless after 10% of its simulation cannot claw back the other 90% of the compute — the process model is "run to completion or `kill -9`."
5. **The hardware is wasted.** These codebases predate 96-core CCD-partitioned CPUs and 546 GB/s unified-memory laptops; they are MPI-shaped or single-threaded-with-OpenMP-sprinkles, allergic to work stealing, and dependent on BLAS binaries tuned for a different decade.

FrankenSim's founding move is to make **one typed algebra** in which geometry, fields, operators, derivatives, error bounds, budgets, provenance, and cancellation are all first-class values that travel *together* through every layer. When the value that flows through the system is `Certified<Field>` — a field plus its derivative hooks plus its interval-verified error bound plus its provenance hash plus its cancellation scope — the archipelago collapses into a continuum. That is not achievable by wrapping existing libraries; it requires the blank slate, and Rust's type system is the only mainstream language substrate strong enough (affine types for resource scopes, const generics for units and dimensions, traits for the representation algebra, fearless concurrency for the executor) to hold it together without a garbage collector or a C ABI in the hot path.

---

## 2. The Decalogue: ten non-negotiable principles

**P1 — Pure memory-safe Rust, Franken-only dependencies.** The runtime dependency set is exactly: `std`, asupersync, FrankenSQLite, FrankenNumpy, FrankenTorch, FrankenScipy, FrankenPandas, FrankenNetworkx. No BLAS, no LAPACK, no C, no C++, no FFI. `unsafe` is permitted only inside a small set of audited leaf modules (SIMD microkernels, arena allocators, memory mapping), each under 300 lines, each with a safe façade, each exhaustively tested including under an in-house model-checking harness that enumerates interleavings of the concurrency primitives.

**P2 — Determinism is a feature, not an accident.** A `deterministic` execution mode guarantees bit-identical results across runs, across thread counts, and (best effort, documented exceptions) across ISAs — achieved via fixed-shape reduction trees, counter-based RNG streams keyed by logical work identity rather than thread identity, and compensated summation. Optimization studies are replayable from a seed and a ledger entry, always.

**P3 — Everything is differentiable or certifiable, ideally both.** Every operator either exposes an adjoint (for gradients) or an interval/Taylor-model enclosure (for certificates), and the flagship paths expose both. "We can't differentiate that" is treated as a bug with a mitigation plan, never a permanent state.

**P4 — Budgets first.** Every operation accepts an explicit budget — accuracy, time, memory — and the system maintains an **Error Ledger** and a **Time Ledger** that compose budgets end-to-end across the pipeline (§11.4). "How accurate is this drag number and where did the error come from" is a query, not a research project.

**P5 — Structure preservation over brute force.** Discretizations preserve the exact sequences of exterior calculus (div∘curl = 0 holds *exactly* on the discrete level); time integrators are symplectic/variational where the physics is Hamiltonian; multiphysics coupling goes through power-conserving ports. Structure preservation is the cheapest accuracy money can't buy: it eliminates whole classes of spurious modes, drift, and instability before a single flop is tuned.

**P6 — Matrix-free and roofline-honest.** We never assemble what we can apply, never store what we can recompute cheaper than a cache miss, and every kernel ships with its arithmetic-intensity analysis and a target expressed as a percentage of the roofline for the machine it's running on (§14).

**P7 — Cancellation-correct compute.** Every unit of work runs inside an asupersync scope; every kernel polls cancellation at tile boundaries with a bounded latency-to-cancel; every iterative solver is a serializable state machine that can pause, migrate, and resume. Cancellation is a *numerical primitive*: the optimizer's sequential tests (§9.6) kill dominated candidates mid-solve and the executor reclaims their cores within a millisecond. This is where asupersync's cancel-correctness stops being plumbing and becomes algorithmic leverage.

**P8 — One data model: complexes and cochains.** The shared substrate for geometry and physics is the cell complex; fields are cochains on it; representations (SDF, NURBS, mesh, voxel) are *charts* of an abstract region; conversions are functors with certified error. There is exactly one canonical in-memory layout for each of these, shared by all kernels.

**P9 — Provenance-complete.** Every artifact — a geometry, a field, a mesh, a Pareto front — is content-addressed, and every operation that produced it is an event in the FrankenSQLite ledger. `explain(artifact)` reconstructs the full causal tree, and any study can be replayed, forked, or audited from the ledger alone.

**P10 — Agent-first ergonomics, human-compatible.** The canonical interface is a typed intermediate representation ("FrankenScript IR") with a machine-readable catalog, structured errors that carry suggested fixes, dimensional analysis enforced at compile time, and the **Five Explicits**: units, seeds, budgets, versions, and capabilities are never implicit, ever (§11).

---

## 3. The Twelve Big Bets

These are the theses that make FrankenSim a leapfrog rather than a re-implementation. Each is one paragraph here and elaborated in its home section.

**Bet 1 — Geometry as a category of charts with a certified-error routing plane. [F/M]** There is no privileged representation. An abstract `Region` is presented through charts (SDF grid, F-rep CSG tree, NURBS patchwork, half-edge mesh, sparse voxel tree, point cloud, lattice graph), conversions between charts are functors carrying *(cost, error-bound)* labels, and a **Rep Router** solves a Pareto shortest-path problem over the chart graph to pick the cheapest conversion chain that keeps composed error inside the caller's budget. Error bounds are interval-verified, not estimated, wherever the chart pair permits. This turns "which representation should I use" — the perennial unforced error of every CAD/simulation handoff — into a solved planning problem. (§7.3)

**Bet 2 — Geometric algebra as the coordinate-free core. [S/F]** Rigid motions, screws, joints, and incidence live in projective geometric algebra Cl(3,0,1); sphere/tangency-rich constructions live in conformal geometric algebra Cl(4,1). Points, lines, planes, circles, and spheres are all blades; intersections are meets; motions are versor sandwiches; twists and wrenches are bivectors. One algebra replaces the usual zoo of quaternion/matrix/Plücker special cases, kills gimbal-class bugs by construction, and makes kinematic constraints for mechanisms and flapping-wing linkages compositional. (§7.7)

**Bet 3 — Structure-preserving physics: FEEC + CutFEM + variational integrators + port-Hamiltonian coupling. [F]** The physics kernel is built on Finite Element Exterior Calculus from day one — fields are differential forms, discrete de Rham complexes are exact, Hodge theory is native — with **cut-cell embedded discretizations (CutFEM with ghost-penalty stabilization) directly on signed distance fields**, so a level-set geometry can be simulated *without ever generating a body-fitted mesh*. Time integration is variational/symplectic where applicable; multiphysics subsystems compose through power-conjugate ports with Dirac-structure interconnection, so partitioned coupling is passive (provably non-energy-creating) by construction. This trio is the single biggest accuracy-and-robustness leapfrog over legacy FEM/CFD. (§8)

**Bet 4 — Adjoint-native, goal-driven everything. [S/F]** Discrete adjoints are designed in from the first line of every solver (differentiate *through the solution*, via the implicit function theorem, never through solver iterations), and mesh/grid adaptivity is driven by dual-weighted-residual estimates of the *design objective* — the simulator spends resolution only where it changes the answer the optimizer cares about. (§8.7, §8.6)

**Bet 5 — The anytime-valid stochastic layer: e-processes and conformal e-prediction. [F/M]** All stochastic estimation (Monte Carlo aero coefficients, seismic response statistics, noisy objective comparisons) runs under test-martingale/e-process confidence sequences, giving *anytime-valid* inference: the optimizer may peek continuously and stop the moment a comparison is decided, with rigorous error control under optional stopping. Surrogate models are wrapped in conformal e-prediction to produce distribution-free, anytime-valid error bands, and the optimizer follows a **certify-or-escalate** policy: trust the surrogate inside its certified band, escalate to higher fidelity outside it. This is the honest, load-bearing realization of "conformal e-martingales": game-theoretic statistics fused with conformal prediction, deployed as the safety layer of the whole optimization stack. (§9.6, §9.7)

**Bet 6 — Optimization as geometry. [S/F]** The optimizer suite is unified under information geometry and Riemannian geometry: CMA-ES and natural evolution strategies are implemented as what they are — natural-gradient flows on the Gaussian family; shape gradients are preconditioned by Sobolev metrics (steepest descent in the right inner product, which is the difference between mesh noise and clean convergence); design spaces that are manifolds (rotations, subspaces, volume-constrained level sets) get genuinely Riemannian optimizers; and the space of shapes itself carries elastic metrics so the system can compute *geodesics between designs* for interpolation, morphing, and trust regions that respect shape geometry. Optimal transport supplies Wasserstein metrics and barycenters for comparing and blending designs and for distributionally robust objectives. (§9)

**Bet 7 — Certificates everywhere. [F/M]** Interval and Taylor-model arithmetic certify geometric predicates, watertightness, Lipschitz bounds for sphere tracing, and enclosure of critical quantities. Sum-of-squares/Lasserre relaxations (via an in-house first-order SDP solver) produce *proofs* of global optimality or infeasibility for low-dimensional polynomial subproblems — including SOS Lyapunov certificates for the stability regions of flight dynamics. Persistent homology of density fields turns topology (number of tunnels, enclosed voids) into constraints the optimizer can actually enforce. FrankenSim doesn't just return answers; where it matters, it returns proofs. (§9.8, §7.8)

**Bet 8 — Cancellation as a numerical primitive; determinism as a contract. [S/M]** See P2/P7. The novel part [M] is the integration: sequential statistical tests (Bet 5) *drive* structured cancellation (asupersync) inside evolutionary and Bayesian optimizers — "e-raced CMA-ES" kills the losing half of a generation the instant the e-process crosses threshold, reallocating cores to survivors, with the whole tournament bit-reproducible from its seed. No system I'm aware of closes this loop. (§5.2, §9.6)

**Bet 9 — Learned accelerators with guarantees. [F]** Neural operators (Fourier neural operators, DeepONets, via FrankenTorch) and classical reduced-order models (POD-Galerkin with DEIM, balanced truncation, Koopman/DMD for unsteady flow) act as 100–10,000× cheap physics inside the optimization loop — but only ever inside their conformal-e-certified validity bands (Bet 5). Machine learning proposes; certified numerics disposes. (§9.7)

**Bet 10 — The Design Ledger. [S]** FrankenSQLite is the system of record: content-addressed artifact store, event-sourced operation log, metric time series, autotuning cache, and lineage graph — WAL-crash-safe, snapshot-isolated for concurrent readers during sweeps, and the substrate for time travel, world forking, and `explain()`. The Ledger is what makes a six-month design campaign a database you can query instead of a directory you fear. (§11.2)

**Bet 11 — Sheaf-theoretic consistency. [F/M]** Cellular sheaves over the patch/tile adjacency complex give the mathematically correct language for "locally represented, globally consistent": stalk = local function space or local representation, restriction map = conversion/trace operator, H⁰ = globally consistent assemblies, and **H¹ = the obstruction**. Concretely: watertightness of a multi-representation model is certified by verifying the interface mismatch cocycle vanishes to tolerance (with interval bounds); domain-decomposition coarse spaces (BDDC-style) are computed as approximate harmonic sections of the interface sheaf; distributed constraint propagation across an assembly is sheaf diffusion. Sheaf cohomology here is not decoration — it is the datatype of gluing, which is precisely the operation multi-representation geometry and domain-decomposed physics perform constantly. (§7.3, §8.9)

**Bet 12 — The pipeline optimizes itself. [F/M]** Cost models learned from the autotuning cache predict runtime and error for every operator on the current machine; the HELM planner co-optimizes *accuracy per second* across the whole pipeline (allocating the error budget between meshing, discretization, solver tolerance, surrogate error, and Monte Carlo noise); and the task DAG's timing is analyzed in the max-plus (tropical) semiring, where the critical path is a tropical eigenvector — so the scheduler knows analytically which kernel to accelerate next. Meta-optimization of the toolchain, running continuously, on your own hardware. (§11.4, §14.3)

---

## 4. Architecture at a glance

```
┌─────────────────────────────────────────────────────────────────────────┐
│  L6  HELM        orchestration · FrankenScript IR · sessions · ledger   │
│                  agent API · budgets · planner · capability governor    │
├─────────────────────────────────────────────────────────────────────────┤
│  L5  LUMEN       spectral path tracing · sphere tracing · sci-viz       │
│                  differentiable rendering                               │
├─────────────────────────────────────────────────────────────────────────┤
│  L4  ASCENT      adjoint/shape/topology optimization · CMA-ES/NES ·     │
│                  Bayesian & multi-fidelity · e-process racing ·         │
│                  SOS certificates · multi-objective                     │
├─────────────────────────────────────────────────────────────────────────┤
│  L3  FLUX        FEEC/DEC · CutFEM-on-SDF · IGA · LBM · BEM+FMM ·       │
│                  vortex methods · contact · variational integrators ·   │
│                  port-Hamiltonian coupling · adjoints · UQ              │
├─────────────────────────────────────────────────────────────────────────┤
│  L2  MORPH       Region/chart algebra · SDF/F-rep/NURBS/mesh/voxel ·    │
│                  Rep Router · meshing · parameterizations · GA ·        │
│                  validity certificates                                  │
├─────────────────────────────────────────────────────────────────────────┤
│  L1  BEDROCK     fs-la GEMM/sparse/eigen · FFT · intervals/Taylor ·     │
│                  autodiff · RNG/QMC · extended precision                │
├─────────────────────────────────────────────────────────────────────────┤
│  L0  SUBSTRATE   asupersync two-lane executor · work stealing · arenas  │
│                  SIMD dispatch · NUMA/topology · determinism · autotune │
├─────────────────────────────────────────────────────────────────────────┤
│  FRANKEN CONSTELLATION   asupersync · FrankenSQLite · FrankenNumpy ·    │
│  FrankenTorch · FrankenScipy · FrankenPandas · FrankenNetworkx          │
└─────────────────────────────────────────────────────────────────────────┘
```

Dataflow in one sentence: HELM lowers an agent's FrankenScript program to a task DAG; the DAG's nodes are budgeted, cancellable kernel invocations scheduled by SUBSTRATE; MORPH supplies geometry as charts of Regions; FLUX turns charts into complexes, cochains, and solved fields with adjoints and error estimates; ASCENT consumes objectives, gradients, and confidence sequences to propose new designs; LUMEN renders anything on demand; and every artifact and event lands in the Ledger.

The Rust workspace is a flat set of `fs-*` crates with a strict acyclic dependency order (full inventory in Appendix A). The discipline that keeps a system like this from rotting: **each crate ships a CONTRACT.md stating its invariants, a conformance test suite that any reimplementation must pass, and zero access to any layer above it.**

### 4.1 The Franken constellation: division of labor

The rule is simple: **FrankenSim implements its own HPC-critical inner loops; the constellation provides substrate, orchestration, ML, graphs, tables, and storage.**

- **asupersync** is the execution substrate — structured concurrency, cancellation tokens, scopes — extended by `fs-exec` with a throughput-oriented work-stealing lane (§5.2). Cancel-correctness is load-bearing for Bets 5 and 8.
- **FrankenSQLite** is the Design Ledger (§11.2): artifacts, events, metrics, lineage, autotune cache. WAL mode, snapshot-isolated readers, single file per project.
- **FrankenNumpy** provides the ndarray interchange type: every FrankenSim buffer exposes zero-copy views as FrankenNumpy arrays, which is how agents, notebooks, and the other Franken libraries touch the data without copies.
- **FrankenTorch** provides the reverse-mode tape used for surrogate training and for differentiating learned components, plus the NN machinery for neural operators, DeepSDF-style implicit shapes, and GP hyperparameter fitting. (Solver adjoints are hand-derived in FLUX, not taped — §8.7 explains why.)
- **FrankenScipy** provides special functions, quadrature baselines, ODE reference integrators, and reference optimizers used as cross-checks in the Gauntlet; ASCENT's production optimizers are native to FrankenSim but validate against FrankenScipy on shared problems.
- **FrankenPandas** is the results plane: metric tables, Pareto front tables, study reports materialize as dataframes for analysis and for auto-generated lab-notebook documents.
- **FrankenNetworkx** powers everything graph-shaped: ground-structure generation for truss layout (§9.5), lattice/TPMS infill graphs, mesh adjacency and partitioning experiments, constraint-propagation graphs, and max-plus critical-path analysis of the task DAG (§14.3).

---

## 5. L0 — SUBSTRATE: hardware, execution, memory, determinism

### 5.1 The two target microarchitectures, honestly

FrankenSim treats "Apple Silicon" and "high-core-count x86" not as compile targets but as *design constraints* with numbers attached. Illustrative reference machines (numbers are order-of-magnitude design inputs, re-measured by the autotuner on first run, never trusted from a spec sheet):

| Property | Apple M4 Max (16c: 12P+4E) | TR PRO 7995WX / 9995WX (96c) | EPYC 9755-class (128c) |
|---|---|---|---|
| Vector ISA | NEON 128-bit ×4 pipes; **SME/SME2** matrix unit | AVX-512 (Zen4 double-pumped / Zen5 full 512b) | AVX-512 (Zen5c) |
| Cache line | **128 B** | 64 B | 64 B |
| Core topology | Unified, P/E asymmetric | 12 CCDs × 8 cores, per-CCD 32 MB L3 islands | 16 CCDs |
| Memory | Unified, ~546 GB/s | 8-ch DDR5, ~330 GB/s | 12-ch DDR5, ~575 GB/s |
| NUMA | UMA | NPS1–NPS4 modes, CCD-local L3 | ditto, more so |
| FP64 peak (order) | ~0.5 TF (NEON) + SME | ~6–10 TF | ~8–12 TF |

Consequences baked into the substrate:

1. **Cache-line-size is a `const`**, resolved per-target (`128` on aarch64-apple, `64` on x86-64), used for padding, false-sharing avoidance, and tile geometry. Getting this wrong on M-series silently halves effective bandwidth on contended structures.
2. **Bandwidth-per-core inverts between the machines.** M4 Max has ~34 GB/s per core; a 96-core TR has ~3.4 GB/s per core. The same algorithm must therefore ship with a bandwidth-rich schedule (fewer, fatter tiles; streaming-friendly) and a bandwidth-starved schedule (aggressive cache blocking, fusion, recomputation-over-reload). The autotuner selects; kernels are written to be schedule-polymorphic.
3. **CCD L3 islands are soft NUMA.** The work-stealing scheduler steals within a CCD first, then within the socket; memory is first-touch initialized by the worker that will own the tile; long-lived field storage is sharded per-CCD.
4. **P/E asymmetry on Apple** is handled by per-core-class throughput measurement: the autotuner benchmarks each kernel per core class and sets stealing weights so E-cores take proportionally smaller tile quanta rather than being ignored or stalling the join.
5. **SIMD dispatch tiers**: a portable-SIMD baseline (compiles everywhere), a NEON path, an AVX-512 path, and an exploratory SME2/streaming-mode matrix path for GEMM-shaped work on M4 [F — intrinsics maturity risk, see §17]. Dispatch is resolved once at startup into function tables; no per-call branching.
6. **Huge pages** (2 MB via `madvise`/THP on Linux; Apple's 16 KB base pages noted) for field arenas; explicit software prefetch in stencil kernels with autotuned prefetch distance.

### 5.2 Execution: the two-lane executor on asupersync

`fs-exec` provides two lanes sharing one topology map:

- **Latency lane** — asupersync's native async scheduling for orchestration: ledger I/O, IR interpretation, progress streaming, watchdogs, inter-study coordination. Milliseconds matter, throughput doesn't.
- **Throughput lane** — a work-stealing fork-join pool (per-core Chase–Lev-style deques, CCD-aware steal order) whose *units of work are tiles* and whose *lifetimes are asupersync scopes*. Every compute task is spawned inside a cancellation scope; the tile protocol requires kernels to poll `ctx.checkpoint()` at tile boundaries, giving a **bounded latency-to-cancel target of ≤ 200 µs** on reference hardware.

The contract that makes P7 real:

```rust
/// Every hot kernel is a TileKernel. Tiles are the unit of scheduling,
/// cancellation, determinism, and NUMA placement.
pub trait TileKernel: Sync {
    type Out: Reduce;                       // fixed-shape reduction tree
    fn tiles(&self) -> TilePlan;            // geometry decided by autotuner
    fn run(&self, t: TileId, ctx: &Cx) -> ControlFlow<Cancelled, Self::Out>;
}
```

`Cx` carries: the cancellation token, the tile-scoped bump arena, the counter-based RNG stream keyed by `(seed, kernel_id, tile_id, iteration)`, the budget accountant, and the ledger handle. Because RNG streams are keyed by *logical* identity, results are independent of which worker ran which tile — the foundation of P2.

Three executor behaviors that most HPC runtimes lack and FrankenSim treats as core:

1. **Speculative races with loser-cancellation.** "Try Delaunay refinement and octree CutFEM discretization concurrently; first to meet the error budget wins; the loser's scope is cancelled and its arenas reclaimed in O(1)."
2. **Resumable solvers.** Krylov, multigrid, and time-stepping loops are explicit state machines (`SolverState: Serialize`); a paused solve can be checkpointed to the Ledger, migrated, resumed, or *forked* (two optimizer branches continuing one warm-started solve).
3. **Statistical preemption.** The e-process layer (§9.6) holds kill-handles to candidate evaluations; crossing an elimination threshold cancels the candidate's entire scope tree mid-flight. Cancel-correctness (asupersync) guarantees no torn state, and arena scoping guarantees no leaks. This is Bet 8's machinery.

Failure containment: a panicking tile poisons its scope, siblings are cancelled, the structured error propagates with full tile provenance, and the study continues or halts per policy — never a process abort mid-campaign.

### 5.3 Memory discipline

- **Scope arenas**: bump allocators tied 1:1 to asupersync scopes; cancellation or completion reclaims in O(1); cross-scope escapes are a compile error (lifetimes do the enforcement). Long-lived artifacts graduate to the content-addressed store.
- **Structure-of-arrays everywhere** hot, via an in-house `#[derive(Soa)]` proc-macro (ours, hence P1-compliant). Arrays are 128-byte aligned unconditionally (superset of both targets' needs).
- **Morton/tile-major field layout**: fields live in Z-order 8³ (or autotuned 4³/16³) tiles — the layout is simultaneously stencil-friendly, sparse-friendly (§7.2's FrankenVDB), cache-oblivious-ish, and gives spatial tasks natural tile identities for the executor.
- **First-touch + shard-by-CCD** for large fields; allocation-site tracking feeds the Ledger so memory regressions are diffable between runs.

### 5.4 The determinism contract (P2)

`ExecMode::Deterministic` guarantees: fixed-shape pairwise reduction trees keyed by tile index (never by arrival order); Neumaier compensated accumulation for global reductions; counter-based RNG (Philox-class, implemented in-house) keyed by logical identity; deterministic tie-breaking in all selection/partition steps; and a cross-ISA report in CI documenting exactly which operations may differ across NEON/AVX-512 (FMA contraction differences) with switches to force strict orderings where bit-equality across ISAs is required. `ExecMode::Fast` relaxes reduction shape for ~5–15% throughput when reproducibility across core counts isn't needed; the mode is part of every ledger event, so a result always knows how it was made.

### 5.5 The autotuner

A first-run (and periodically refreshed) calibration pass measures, per machine: sustained bandwidth per core class, GEMM microkernel throughput across tile shapes, stencil throughput vs. tile size and prefetch distance, steal latencies, and reduction costs. Results persist in the Ledger's `tune` table keyed by a hardware topology hash. Every `TilePlan` and every cost model in the HELM planner (Bet 12) reads from this table. Nothing performance-critical is hardcoded to either microarchitecture.

---

## 6. L1 — BEDROCK: the numerical foundations

Everything here is written in-house, in safe Rust with audited SIMD leaves, because the entire point is kernels shaped for *our* data layouts and *our* two machines rather than for a 1979 Fortran calling convention.

### 6.1 Dense linear algebra (`fs-la`) [S]

- **GEMM in the BLIS style**: an arch-specific register microkernel (e.g., 2×… — shapes chosen by the autotuner, typically 8–14 rows × 2–4 vector columns) wrapped in autotuned KC/MC/NC cache blocking with packing into arena buffers. Targets in §14. f64, f32, and f32-storage/f64-accumulate mixed mode.
- **Batched small dense** is a first-class citizen, because FEM element matrices (6×6 … 48×48) dominate assembly-adjacent work: the layout interleaves *across* the batch so SIMD lanes run across elements, not within one tiny matrix — routinely 5–15× over loop-over-LAPACK-calls patterns.
- Blocked Cholesky/LU/QR; **TSQR** (tree QR) for tall-skinny blocks (Krylov basis orthogonalization, POD snapshots); one-sided Jacobi SVD for accuracy-critical small problems; randomized range finders + randomized SVD/Nyström for big low-rank work [F]; **LOBPCG** and Lanczos (with selective reorthogonalization) as the eigen-solvers — chosen precisely because they are matrix-free-native (modal analysis, buckling, stability eigenproblems all reduce to operator applications).
- **Mixed precision with iterative refinement** as a default policy: factor/precondition in f32, refine to f64, escalate to double-double (in-house error-free-transformation arithmetic) only where the Error Ledger demands it.
- Randomized NLA utilities: sketch-and-precondition least squares, Hutchinson/Hutch++ stochastic trace estimation (log-det, sensitivity aggregation in topology optimization) [F].

### 6.2 Sparse (`fs-sparse`) [S]

Formats: CSR for generality, **BSR** (block CSR) matched to FEM block structure, and **SELL-C-σ** for SIMD-regular SpMV on both ISAs. Assembly path is COO → parallel sort → dedupe → target format, all tiled. Preconditioner components: Chebyshev polynomial smoothers (the many-core-friendly choice — no sequential dependencies), IC(0)/ILU(0) with level scheduling where serial chains are tolerable, and **smoothed-aggregation AMG** built in-house — primarily employed as the coarse solver under p-multigrid (§8.9). Direct sparse factorization is deliberately minimized (P6): where unavoidable (small coarse grids, 1D-ish problems), a supernodal Cholesky exists but is not a load-bearing wall.

### 6.3 FFT (`fs-fft`) [S]

Iterative Stockham autosort, mixed radix 4/8, real-input split paths, cache-blocked transposes, 3D via pencil decomposition across the tile executor. Consumers: spectral Poisson solves on boxes, Fourier neural operators, optimal-transport convolutions, rendering image ops, and Chebyshev transforms (via DCT) for the function layer below.

### 6.4 Certified arithmetic (`fs-ivl`) [S/F]

Interval arithmetic with rigorous outward rounding via `next_up`/`next_down` directed nudging (no global rounding-mode fiddling — Rust-safe and thread-safe); affine arithmetic for correlation-aware tightness; **Taylor models** (multivariate polynomial + interval remainder) for high-order enclosures; interval Newton and the Krawczyk operator for certified root isolation. Shewchuk-style adaptive-precision floating-point expansions provide **exact geometric predicates** (`orient3d`, `insphere`) — the difference between a mesher that works and a mesher that works until Tuesday. Double-double/quad-double types round out the escalation ladder. This crate is what the word "Certified" means everywhere else in this document.

### 6.5 Function objects (`fs-cheb`) — compute with functions [F]

A chebfun-inspired layer: 1D/2D/3D functions represented by adaptively truncated Chebyshev/Fourier expansions with automatic degree selection to near machine precision; algebra, calculus, rootfinding on *functions* as values. Used for: airfoil sections and spline profiles, boundary condition specification, 1D stability eigenproblems (Orr–Sommerfeld in §15.3), and anywhere an agent should say "the inlet velocity profile is *this function*" instead of shipping a lookup table.

### 6.6 Automatic differentiation (`fs-ad`) [S]

Three regimes, chosen per situation and composable:

1. **Forward-mode duals** `Dual<const N: usize>` with SIMD-lane packing — for few-parameter sensitivities and for stress-testing hand adjoints.
2. **Reverse-mode tape** — FrankenTorch's tape for learned components and surrogate training.
3. **Adjoint-by-construction** — the workhorse: solvers differentiate through *the solution*, not the iterations, via the implicit function theorem (solve the transposed system with the same preconditioner infrastructure); time-dependent adjoints use binomial (revolve-style) checkpointing schedules integrated with the Ledger so checkpoints spill to storage gracefully. §8.7 owns the details.

### 6.7 Randomness and quasi-randomness (`fs-rand`) [S]

Counter-based Philox-class generator (in-house) for reproducible parallel streams; Sobol' sequences with Joe–Kuo direction numbers and Owen scrambling; rank-1 lattice rules constructed component-by-component — because UQ workloads (§8.8) get 1–2 orders of magnitude from QMC over MC at these dimensionalities, and low-discrepancy points also drive sampling in Bayesian optimization and light transport.

---

## 7. L2 — MORPH: the geometry kernel

MORPH's job: represent, construct, interrogate, convert, mesh, and parameterize geometry — with certified validity, native derivatives with respect to design parameters, and no privileged representation.

### 7.1 The core abstraction: Regions and charts [F/M]

The canonical object is an abstract `Region` — semantically, a measurable subset of ℝ³ with (piecewise) smooth boundary — that is never stored directly. It is *presented* through **charts**, each a concrete representation implementing the `Chart` trait: query interfaces (inside/outside, distance, closest point, ray intersection, curvature, integrals over the region/boundary) plus a declared error model relative to the abstract region. A `Region` holds one or more charts plus the provenance of how each was obtained; agreement between charts is a checkable proposition, not an assumption.

This framing does three things no conventional kernel does: it makes "the same shape held three ways" a normal, coherent state instead of a synchronization bug factory; it makes every conversion's error explicit and composable (feeding the Error Ledger, P4); and it lets the Rep Router (§7.3) choose representations per *operation* rather than per project.

### 7.2 The chart zoo

- **SDF charts** [S]: dense tiled grids (Morton 8³ tiles, f32 storage/f64 evaluation), **FrankenVDB** — an in-house sparse hierarchical tile tree (hash root → 32³ internal → 8³ bitmasked leaves) for huge sparse domains — plus adaptively-sampled SDFs (octree with per-cell polynomial fits, error-controlled) and hash-grid variants. Narrow-band mode for level-set evolution.
- **F-rep charts** [S]: CSG DAGs over implicit primitives and blends (R-functions for differentiable Booleans — critical: hard min/max Booleans have derivative discontinuities that poison shape optimization; R-function blends are the smooth, principled alternative), with interval and Lipschitz-bound evaluators derived automatically from the DAG (feeding certified sphere tracing, §10.2, and certified inside/outside tests).
- **Neural implicit charts** [F]: small MLPs (FrankenTorch) as shapes, DeepSDF-style, with Lipschitz-constrained layers so certified bounds remain available; used as ultra-compact, inherently differentiable design parameterizations.
- **NURBS/B-rep charts** [S/F]: full rational B-spline surfaces/curves with exact knot insertion and degree elevation; trimmed-patch B-reps with trim curves held in exact rational form plus interval boxes. Watertight trimmed-NURBS Booleans are the graveyard of CAD kernels; our position (honest): Boolean ops route through F-rep/SDF charts by default and re-fit splines afterward when a spline chart is required [F], while direct B-rep Booleans remain a constrained, certificate-gated capability rather than a promise. IGA (§8.1) is the payoff for keeping first-class splines: geometry basis = analysis basis.
- **Mesh charts** [S]: half-edge for surfaces, oriented tet/hex complexes for volumes; a repair suite (polygon-soup healing via **generalized winding numbers** for robust inside/outside on broken input — the single most effective modern trick in mesh robustness, accelerated by an octree dipole approximation to run fast on millions of triangles).
- **Voxel/point/lattice charts** [S]: occupancy and multi-material voxel fields on FrankenVDB; point clouds with normals (fitting targets); explicit lattice/strut graphs via FrankenNetworkx (ground structures, TPMS approximations, infill).

### 7.3 The Rep Router and the sheaf of representations [F/M]

Conversions form a directed multigraph: nodes are chart types, edges are converters annotated with *(cost model, error model, certificate availability)*. Requests arrive as "give me a chart supporting operation X with error ≤ ε under budget B"; the Router solves a **Pareto shortest-path** problem (cost × composed error, with certified edges preferred) over this graph — FrankenNetworkx territory — and executes the winning chain, recording actual cost/error back into the tune tables so the models improve. Representative certified edges:

- mesh → SDF: exact point-triangle distance + winding-number sign, per-tile error = 0 up to fp rounding [S];
- SDF → mesh: **dual contouring** with QEF vertex placement (sharp-feature capable) and a manifold variant, with interval verification that the extracted surface brackets the zero set within tolerance [S/F];
- NURBS → SDF: certified distance via convex-hull/Bézier bounding + interval Newton polishing [F];
- SDF → NURBS: MLS/spline re-fitting with thin-plate smoothing and a *reported, measured* Hausdorff bound (estimated, promoted to certified where the interval machinery can close it) [F].

**Gluing as sheaf cohomology (Bet 11).** A model whose patches live in different charts is a cellular sheaf over the patch-adjacency complex: stalks are local chart spaces, restriction maps are trace/conversion operators to shared interfaces. Global consistency = existence of a section (H⁰); the **watertightness certificate** is an interval-verified bound ‖δs‖ ≤ tol on the interface mismatch cocycle — if the coboundary can't be driven below tolerance, the H¹ obstruction is reported *with the offending interface cells attached*, which is exactly the diagnostic an agent needs to fix a leaky model. The same machinery later powers domain-decomposition coarse spaces (§8.9). This is the most "galaxy-brained" component that is nonetheless straightforwardly implementable: cellular sheaves are finite linear algebra over an adjacency complex.

### 7.4 Queries and constructions [S]

Closest point, signed distance, ray casting, clearance fields, offsets (exact-ish via distance-field composition; spline offsets via re-fit), Minkowski sums via SDF convolution/max-plus composition, medial axis via filtered Voronoi (λ- and angle-criteria) — the medial axis doubles as the *thickness oracle* for manufacturability constraints in ASCENT. Curvature (principal, mean, Gaussian) available on every chart type, by analytic differentiation where the chart permits and by certified stencils otherwise.

### 7.5 Meshing — when we must [S/F]

CutFEM-on-SDF (§8.1) exists precisely to make meshing optional inside optimization loops. When body-fitted meshes are wanted (final verification, shells, export):

- **Tet pipeline** [S]: BRIO-ordered incremental Delaunay on exact predicates → constrained boundary recovery → Ruppert/Shewchuk-style quality refinement → sliver exudation via weighted-Delaunay perturbation. Fully tiled and parallel via domain coloring.
- **Surface remeshing** [S]: isotropic (incremental edge ops toward target length) and **anisotropic under a metric tensor field** — the metric is *supplied by FLUX's error estimator*, closing the goal-oriented adaptivity loop (§8.6).
- **Hex-dominant** [F]: octahedral frame-field computation (spherical-harmonic SH9 parameterization with MBO-style smoothing) → parameterization → hex extraction, with polycube fallback; honestly flagged as the hardest meshing deliverable in the plan (§17), mitigated by the fact that IGA and CutFEM cover most hex use-cases at higher accuracy anyway.

### 7.6 Design parameterizations — the levers ASCENT pulls [S/F]

A parameterization is a differentiable map θ → Region, shipped with its Jacobian action (`dRegion/dθ` as a boundary velocity field or density perturbation): free-form deformation lattices; RBF morphs with compact support; direct spline control points; level-set velocity fields (the topology-optimization workhorse); density fields (SIMP); **manifold-harmonic bases** — eigenfunctions of the Laplace–Beltrami operator on the current shape, giving an automatically smooth, hierarchical, low-dimensional "shape spectrum" ideal for CMA-ES-scale global search [F]; neural implicits (§7.2); and constructive/procedural programs (F-rep DAGs with exposed parameters). Constraint primitives are first-class and differentiable: minimum thickness (via medial oracle), draft angle, symmetry groups (enforced by construction — quotient the parameterization by the group), bounding envelopes, volume/mass.

### 7.7 The geometric-algebra layer (`fs-ga`) [S/F]

Cl(3,0,1) (projective GA) is the kinematics substrate: motors for rigid motions, join/meet for incidence, screws/twists/wrenches as bivectors, clean interpolation (motor logarithms) with no gimbal pathologies; Cl(4,1) (conformal GA) handles sphere/circle/tangency-rich construction and blending problems where it is genuinely superior. Implementation: types generated by an in-house const-evaluated code generator over the algebra's multiplication table — fully monomorphized, SIMD-friendly, no runtime blade bookkeeping. GA is deliberately an *internal excellence* layer with conventional Vec3/quaternion façades at the API boundary, so agents never pay a formalism tax to move a part.

### 7.8 Validity and topology certificates [F]

Beyond watertightness (§7.3): manifoldness checks; self-intersection freedom via interval-arithmetic broad/narrow phases (a *certificate*, not a sampling heuristic); and **persistent homology** of density/occupancy fields (in-house cubical complex reduction, chunked and parallel) reporting Betti numbers with persistence — consumed by ASCENT as differentiable-ish topology constraints ("exactly one tunnel; no enclosed voids" — castability, drainability, cable routing) via persistence-diagram penalties [F/M].

---

## 8. L3 — FLUX: the physics kernel

FLUX's job: turn charts into discrete complexes and cochains, apply and invert physical operators at roofline speed, integrate in time with preserved structure, couple physics passively, and hand ASCENT exact discrete gradients and honest error bars.

### 8.1 The unifying formalism: exterior calculus on complexes [F]

Fields are cochains (discrete differential forms) on cell complexes; operators are built from the exterior derivative d (exact, purely combinatorial — an incidence matrix) and discrete Hodge stars (where all metric information and all approximation lives). Element families follow **Finite Element Exterior Calculus**: Whitney forms at lowest order, hierarchical high-order P_rΛ^k / P_r⁻Λ^k families on simplices, tensor-product families on cubes/octree cells. Why this is the spine and not a garnish:

- The discrete de Rham sequence is **exact**: grad→curl→div identities hold to machine precision, killing spurious pressure modes, spurious EM modes, and checkerboarding *by construction* rather than by stabilization folklore.
- **Cohomology is computable and used**: harmonic cochains give the correct treatment of multiply-connected domains (circulation around a wing → lift; pressure fields in domains with handles), where naive FEM silently produces the wrong answer or an underdetermined system.
- One assembly/apply machinery serves elasticity, flow, heat, and (later) electromagnetics — the "generic, abstract, maximally reusable mechanics" requirement is satisfied by choosing the formalism in which those physics *are already the same code*.

Three discretization frontends share this machinery:

1. **Body-fitted simplicial FEEC** [S] on MORPH's tet meshes.
2. **CutFEM on SDFs** [F] — octree background grids (FrankenVDB-aligned) with cut cells at the zero level set, ghost-penalty stabilization for small-cut robustness, and aggregated-element fallback. **This is the marquee bridge**: level-set/SDF geometry is simulated with FEM-grade accuracy and *zero meshing*, which is what makes topology optimization loops (§9.5) run at interactive cadence. Non-uniform resolution comes free from the octree; the error estimator (§8.6) drives refinement.
3. **Isogeometric analysis (IGA)** [S/F] — Galerkin directly on NURBS/B-spline spaces, multi-patch coupling via mortar methods; Kirchhoff–Love shells on spline surfaces (the right tool for thin aero skins and vessel walls). Geometry basis = analysis basis: the CAD→mesh information massacre simply never happens.

### 8.2 Solid mechanics [S/F]

Linear elasticity through finite-strain hyperelasticity (Neo-Hookean, Mooney–Rivlin, Ogden), with locking-free mixed formulations (TDNNS-style / weakly imposed symmetry) where near-incompressibility or bending-domination demands it. Reduced models: geometrically exact **Cosserat rods** (Lie-group state on SE(3)) for members and struts; Kirchhoff–Love IGA shells. Inelasticity: J2 plasticity with return mapping; fiber-section beam elements with Mander concrete + Menegotto–Pinto steel laws for the reinforced-concrete frame use-case. Stability: geometric-stiffness buckling eigenproblems (LOBPCG), **pseudo-arclength continuation** through limit points, and Koiter-style post-buckling asymptotics [F] — because the steel-minimizing frame that ignores buckling is a paperclip. Contact: SDF-native barrier formulation in the Incremental Potential Contact family [F] — smooth, robust, intersection-free by construction, and unusually easy for us because penetration depth queries against SDF charts are MORPH one-liners.

### 8.3 Fluids [S/F]

Multiple solvers, deliberately, because the flagship problems span regimes:

- **Incompressible Navier–Stokes, FEEC-native** [F]: H(div)-conforming velocity spaces with pressure-robust discretization (velocity error independent of pressure — a modern correctness property most production codes lack), upwinded DG convection, projection and coupled variants.
- **Lattice Boltzmann (`fs-lbm`)** [S] — the many-core queen: D3Q19/D3Q27 on sparse FrankenVDB tiles, **cumulant/central-moment collision** for high-Re stability, interpolated (Bouzidi-type) curved boundaries sampled *directly from SDF charts*, octree grid refinement with rescaling, thermal double-population, non-Newtonian via local relaxation (power-law/Carreau — "a liquid of a given viscosity" is a parameter, including shear-dependent ones), and **free-surface LBM** (mass-tracking VOF-style) for pouring and sloshing. LBM is embarrassingly tile-parallel, bandwidth-bound with a precisely computable roofline (§14), and pairs perfectly with the substrate's Morton tiles. A *lattice-scaling assistant* automates the dx/dt/τ/Mach bookkeeping under stability constraints — a chronic source of human LBM error that agents should never be exposed to raw.
- **Boundary elements + fast multipole (`fs-bem`, `fs-fmm`)** [F]: kernel-independent black-box FMM (Chebyshev interpolation based) driving Laplace BEM — potential-flow panel methods with Kutta condition and free wakes for O(N) exterior aerodynamics screening; elastostatic BEM later.
- **Vortex particle methods** [F]: Lagrangian vorticity with FMM Biot–Savart, particle-strength-exchange viscosity, hybridized with the BEM surface solution — the natural tool for unsteady, wake-dominated flapping flight (§15.1).
- Turbulence honesty: LES (Smagorinsky/WALE) on LBM and NS for resolved-eddy regimes; RANS-style closures only as clearly-labeled low-fidelity screens; the multi-fidelity optimizer (§9.7) is the system's answer to "LES everywhere is unaffordable," not wishful meshing.

### 8.4 Multiphysics composition: port-Hamiltonian Dirac structures [F]

Subsystems (structure, fluid, thermal, later EM) expose **power-conjugate ports** (effort/flow pairs: force/velocity, pressure/flux, temperature/entropy-flow); composition is interconnection through a Dirac structure, which conserves power *exactly* at the discrete level. Consequences: partitioned coupling schemes are **passive by construction** — the added-mass instabilities and energy-drift pathologies of ad-hoc FSI staggering are structurally excluded — and multirate co-simulation (fast fluid, slow structure) inherits stability guarantees. Interface acceleration via Aitken and interface quasi-Newton (IQN-ILS) when strong coupling is needed. Compositionally, this is the "wiring diagram of open systems" picture made executable: FrankenSim assemblies are literally morphism composition in the category of open port-Hamiltonian systems — the category theory is not decorative, it is the coupling API's type discipline.

### 8.5 Time integration [S/F]

Variational/symplectic integrators for conservative dynamics (energy behavior over 10⁶ steps that dissipative-by-accident schemes cannot match); Lie-group integrators on SE(3) for rods and rigid bodies (no quaternion drift, no renormalization hacks); Newmark/generalized-α for structural dynamics with controllable numerical dissipation; IMEX splittings for stiff-nonstiff mixtures; exponential integrators where the stiff part is linear; embedded-pair adaptivity with PI step-size controllers. Every integrator is a resumable state machine (P7) and every one ships its adjoint (§8.7).

### 8.6 Adaptivity: spend resolution only where the objective feels it [F]

Non-uniform resolution appears in four coordinated forms — octree h-refinement (CutFEM/LBM), anisotropic metric-driven remeshing (body-fitted), p-enrichment (high-order FEEC), and wavelet-style coefficient thresholding on tiles — all driven by one signal: **dual-weighted-residual (goal-oriented) error estimates**, where the "goal" is ASCENT's actual objective functional. The adjoint solution weights local residuals by their influence on the quantity being optimized; refinement follows the product. The philosophical point: in a design-optimization system, "accurate simulation" is not the goal — *accurate objective and gradient* is — and DWR is the estimator that knows the difference.

### 8.7 Adjoints and shape derivatives — gradient truth [S/F]

- **Discrete adjoints everywhere**: hand-derived transposed operators sharing the primal's matrix-free infrastructure and preconditioners; differentiation through solvers via the implicit function theorem (never through Krylov iterations — that path is both wasteful and wrong); time-dependent adjoints under binomial checkpointing with Ledger-backed spill.
- **Shape gradients** via Hadamard boundary formulas evaluated on FEEC traces, and via level-set/density chain rules for CutFEM/SIMP paths; **Sobolev (H¹) gradient smoothing** applied as the Riesz representation step — the single most important practical trick in shape optimization, turning raw boundary gradients (mesh-noise amplifiers) into smooth descent directions, and slotting neatly into Bet 6's "right inner product" doctrine.
- Gradient verification is a Gauntlet gate: every adjoint is checked against forward-mode duals and complex-step-quality finite differences on randomized problems in CI; a solver without a passing gradient check cannot merge.
- Honesty note: free-surface LBM adjoints are a research swamp; the plan's position is continuous-adjoint NS for gradient-based fluid shape steps, LBM under gradient-free/surrogate-gradient optimization, revisited later [M].

### 8.8 Uncertainty quantification [S/F]

Random inputs as Karhunen–Loève expansions of random fields (material scatter, geometric tolerance, gust spectra); propagation via polynomial chaos for smooth low-dimensional cases, QMC (§6.7) as the general workhorse, and **multilevel Monte Carlo** coupling coarse/fine discretizations for order-of-magnitude variance-per-cost wins [F]; seismic specifics: stochastic ground motion (Kanai–Tajimi-class spectra), response-spectrum CQC fast path, nonlinear time-history slow path, fragility via incremental dynamic analysis. Every stochastic estimate runs under the e-process layer (§9.6), so UQ inherits anytime-valid stopping: sample until the confidence sequence is tight enough for the decision at hand, then stop — automatically, validly.

### 8.9 Solvers and preconditioning [S/F]

The default stack is **matrix-free p-multigrid** (high-order FEEC operators applied via sum-factorization) with smoothed-aggregation AMG on the assembled lowest-order coarse problem, Chebyshev smoothing, wrapped in CG/GMRES/MINRES with mixed-precision inner work + f64 refinement. Saddle-point systems (Stokes, contact, incompressibility) get block preconditioners with Schur-complement approximations. Domain decomposition (BDDC/FETI-DP family) provides the extreme-core-count path; its coarse space is computed as approximate **harmonic sections of the interface sheaf** — the same cellular-sheaf machinery as §7.3, now earning its keep twice (Bet 11). All solvers: resumable, cancellable, deterministic-mode capable, adjoint-equipped.

---

## 9. L4 — ASCENT: the optimization kernel

ASCENT's job: given `minimize J(u(θ), θ) subject to physics(u, θ) = 0, constraints(θ) ≤ 0, θ ∈ M`, find the best designs — with gradients when they exist, guarantees when they're possible, statistical validity when there's noise, and maximal exploitation of 100+ cores when there's a population.

### 9.1 The problem IR [S]

Optimization problems are data: typed objective/constraint graphs over manifold-valued variables, with PDE constraints referencing FLUX studies, stochastic operators (expectations, CVaR, quantiles) referencing UQ configurations, and budgets attached per P4. Multi-objective, chance-constrained, and bilevel structures are representable. Because problems are data, they are storable, diffable, replayable, and — crucially for agents — *constructible incrementally with validation at each step*.

### 9.2 The gradient stack [S]

FLUX adjoint → shape/density gradient → **Sobolev/Riesz smoothing** (§8.7) → optional Riemannian projection (tangent space of the constraint manifold) → optimizer. First-order methods: L-BFGS (the default), Adam-family for stochastic gradients; second-order: matrix-free trust-region Newton–Krylov using Hessian-vector products assembled from second-order adjoints; constraints via augmented Lagrangian (robust default) and an interior-point option; SQP for tightly-constrained small-dimension polish. Manifold variants (Stiefel, sphere, SO(3), fixed-volume level sets) implement retraction + vector transport so "optimize an orientation" never becomes "optimize 9 numbers and renormalize when it explodes."

### 9.3 Derivative-free and global search [S/F]

**CMA-ES** as the flagship (full covariance, rank-µ update, BIPOP restarts; separable and low-rank variants for dim > ~200), implemented explicitly as **natural-gradient flow on the Gaussian family** per the information-geometric-optimization view — which is not aesthetic name-dropping: it dictates the correct step-size/covariance couplings and unifies CMA-ES with natural evolution strategies in one codebase (Bet 6). Alongside: differential evolution, DIRECT for Lipschitz-ish global structure, Nelder–Mead for cheap polish, and trust-region DFO with quadratic interpolation models for expensive smooth-ish objectives. Population methods are where the 96–128 core machines shine: a CMA-ES generation is an embarrassingly parallel wave of FLUX evaluations, scheduled as sibling scopes, raced under §9.6.

### 9.4 Bayesian and multi-fidelity optimization [F]

In-house Gaussian processes (Matérn kernels, exact to ~10⁴ points, inducing-point sparse beyond), batched q-EI/q-NEI acquisition via MC with reparameterized gradients (FrankenTorch), trust-region BO for the 30–300-dim regime, and **multi-fidelity acquisition** where fidelity is a controllable knob (LBM resolution, mesh level, panel vs. LES) with cost-aware selection — the mathematically principled version of "screen cheap, verify expensive." GP hyperparameters fit by marginal likelihood with QMC-seeded multistart.

### 9.5 Topology and layout optimization — generative design with adult supervision [S/F]

- **Density (SIMP)** [S]: Helmholtz PDE filtering (reusing FLUX's Poisson machinery), Heaviside projection with continuation, robust erode/dilate formulations for guaranteed length scale, stress constraints via adaptive p-norm aggregation, compliance/volume/eigenfrequency objectives. Runs on CutFEM octrees → no remeshing per iteration.
- **Level-set** [S/F]: shape-gradient velocity advection (WENO on octree narrow bands), fast-iterative-method redistancing (the parallel-friendly cousin of fast marching), **topological derivatives** to nucleate holes — genuine topology changes with mathematical justification rather than heuristic hole-punching.
- **Truss/frame layout** [S/F]: ground structures generated by FrankenNetworkx (k-NN candidate members under fabrication rules); member-force **layout optimization as a linear/second-order-cone program**, solved by an in-house first-order primal-dual (PDHG-class) solver that scales to millions of candidate members; Michell-continuum limits as analytic sanity oracles; then a nonlinear sizing/geometry pass with Euler and code-based buckling constraints, joint-cost regularization, and constructability rules (member catalogs, angle limits) as graph constraints. This is the steel-and-concrete flagship's engine (§15.2).
- **Lattice/infill** [F]: numerical homogenization of unit cells (TPMS, strut families) → graded-property macro-optimization → de-homogenization back to printable/buildable micro-geometry, with FrankenNetworkx carrying the conforming lattice graphs.
- **Topology constraints via persistence** [F/M]: Betti-number targets enforced through persistence-diagram penalties (§7.8) — the only principled way to tell an optimizer "no enclosed voids" without hand-crafted heuristics.

### 9.6 The anytime-valid racing layer — "conformal e-martingales," cashed out [F/M]

The statistical core of Bet 5. Every noisy comparison ("is design A's CVaR-drag lower than B's?") is monitored by an **e-process / test-martingale confidence sequence** (betting-style construction in the Waudby-Smith–Ramdas lineage): valid at every sample size simultaneously, immune to optional stopping, composable across comparisons via e-value arithmetic, with family-wise control over thousands of candidates via **e-Benjamini–Hochberg**. Operationally:

- **e-raced generations** [M]: within a CMA-ES/NES generation or a successive-halving bracket, per-candidate Monte Carlo evaluation streams feed the e-processes; the moment a candidate's elimination e-value crosses threshold, ASCENT cancels its entire evaluation scope (asupersync, ≤ 200 µs to reclaim) and reassigns cores to survivors. Expected effect on noisy objectives: 2–5× more generations per core-hour at *identical* statistical guarantees, and bit-reproducible tournaments in deterministic mode.
- **Anytime-valid UQ stopping**: MLMC/QMC estimators stop exactly when the confidence sequence is tight enough for the pending decision — no conservative pre-registered sample sizes, no invalid peeking.

### 9.7 Surrogates under certificates: certify-or-escalate [F/M]

Neural operators (FNO on box-resampled fields, DeepONet for parametric families) and classical ROMs (POD-Galerkin+DEIM, balanced truncation, Koopman/DMD for unsteady flows) are trained continuously on the Ledger's accumulating (design → field/objective) corpus. Each surrogate is wrapped in **conformal e-prediction**: distribution-free prediction bands with anytime validity under the e-value formulation, recalibrated online as new ground-truth solves arrive. The optimizer's policy is mechanical: *inside* the certified band and band-width below decision-relevance → use the surrogate (10²–10⁴× cheaper); *outside* or band too wide → escalate fidelity, and the new solve both decides the point and tightens the band. ML proposes, certified numerics disposes; the system gets faster every week it runs without ever silently trusting a hallucinated flow field.

### 9.8 Proof-carrying optimization [F/M]

For low-dimensional polynomial subproblems (≤ ~8–10 variables after parameterization: linkage synthesis, section-shape polish, safety-margin verification), **Lasserre/SOS moment relaxations** solved by an in-house first-order SDP solver (Burer–Monteiro low-rank with dual certificate extraction) return *global* optimality or infeasibility certificates. The show-stealer application: **SOS Lyapunov functions certifying regions of attraction for trim states of the aircraft's flight dynamics** — "stability" in the objective becomes a proven set, not an eigenvalue vibe (§15.1). Complementary certificate machinery: interval branch-and-bound for small global problems, and KKT-residual certificates attached to every returned local optimum so downstream consumers can distinguish "converged" from "stalled."

### 9.9 Multi-objective and robust formulations [S/F]

NSGA-II/III and MOEA/D for evolutionary Pareto exploration (population = core count, naturally); hypervolume (WFG algorithm) and knee-point detection; scalarization sweeps for gradient-based Pareto tracing; interactive steering — an agent can re-weight objectives mid-campaign and the study forks (P9) rather than mutates. Robustness: CVaR via the Rockafellar–Uryasev reformulation (turning tail-risk minimization into a smooth problem FLUX gradients can serve); **distributionally robust optimization over Wasserstein balls** [F] — optimal transport (Sinkhorn-scaled, FFT-accelerated where grids allow) supplying both the ambiguity sets here and the shape-difference metrics/barycenters used for design-space diversity maintenance and semantic diffs (§11.5).

---

## 10. L5 — LUMEN: the rendering kernel

LUMEN's job: physically-based still/animated imagery of designs and physics, direct scientific visualization of every FLUX field type, and differentiable rendering as an optimization sensor — all CPU-native and tile-parallel.

### 10.1 Physically-based light transport [S/F]

Unbiased spectral path tracing in the Maxwell-render tradition: hero-wavelength spectral sampling with MIS, next-event estimation, light-BVH sampling for many emitters, layered measured-spectrum materials, Beer–Lambert and heterogeneous (Woodcock-tracked) volumetrics — the latter pointed directly at LBM density/temperature fields, so "render the pour" is the same machinery as "render the product shot." Wide (8-way) BVH with SIMD traversal, watertight ray-triangle tests, ray streams sorted for coherence; progressive tile streaming to HELM so agents watch convergence live. Optional polarization (Mueller calculus) behind a flag for optics-adjacent studies [F].

### 10.2 Geometry backends without conversion [S]

Sphere tracing renders SDF/F-rep charts *directly* using certified Lipschitz bounds from `fs-ivl` — no meshing for visualization, and step sizes that provably never tunnel. NURBS patches trace via Bézier-clipping-seeded Newton iteration; meshes trace natively. The Rep Router means "render whatever chart exists" is the default, not a pipeline stage.

### 10.3 Scientific visualization [S]

Dual-contoured isosurfaces with sharp features; direct volume rendering with preintegrated transfer functions; line-integral convolution and streamline/streakline transport for flow; tensor glyphs (stress ellipsoids) and principal-direction ribbons for load paths; Morse–Smale/Reeb-graph structural overlays (vortex skeletons, load-path skeletons) [F] — the topological summaries that make a 10⁸-cell field legible to an agent in one image. All viz shares LUMEN's tiling and streams progressively.

### 10.4 Differentiable rendering [F]

Reparameterized/edge-aware gradient estimators for silhouette-and-shading derivatives, enabling inverse problems — match a target photograph, a sketch, or a legacy artifact's appearance as an objective term — and closing the loop where aesthetic constraints join physical ones in ASCENT.

### 10.5 Image plumbing [S]

In-house PNG and OpenEXR writers, à-trous denoiser (optional, clearly labeled as biased), film response curves, and deterministic-mode rendering (per-tile Philox streams) so an image is as replayable as a solve.

---

## 11. L6 — HELM: orchestration, sessions, and the agent interface

Everything below this layer is *callable*. HELM makes it **programmable, replayable, and governable** — for AI agents first, humans equally. HELM is where the Decalogue's P4 (budgets), P9 (provenance), and P10 (agent-first) become running code.

### 11.1 FrankenScript IR — the canonical interface [S/F]

The system's one true interface is a **typed, versioned intermediate representation**, not a pile of function signatures. Two isomorphic concrete syntaxes — canonical s-expressions and a lossless JSON mapping — so agents emit whichever their tooling prefers; both parse to the same typed AST. Values in the IR are the system's real nouns: dimensioned quantities, Regions and charts, fields, operators, budgets, studies, capability tokens.

Before anything executes, the IR passes **static admission**: dimensional analysis (a pressure cannot be added to a stress *in the type system*, via const-generic `Qty` types — Appendix B), chart-compatibility checking through the Rep Router, budget-feasibility screening against learned cost models, and capability sufficiency. An ill-typed study is rejected in milliseconds with a structured diagnosis, not discovered at hour six.

A flavor of the surface (the spout study, abbreviated; fuller example in Appendix C):

```lisp
(study "spout-laminar-v3"
  (seed 0x5EED0001) (versions (constellation :lock "2026-07"))
  (budget (wall 2h) (mem 96GiB) (qoi-rel-error 2e-2))
  (let vessel (frep (revolve (cheb-profile "body.chb")) (fillet :edge lip :r 3mm)))
  (let lever  (xform.level-set-velocity vessel :band 12mm :dof 4096))
  (let pour   (flux.free-surface-lbm vessel
                (fluid :model (carreau :mu0 0.12Pa*s :n 0.8) :sigma 0.061N/m)
                (schedule :rate 0.5L/s :tilt (ramp 0deg 65deg 3s))))
  (let J (min (perturbation-growth pour :at lip :modes (1 .. 8))))
  (ascent.optimize J :over lever :method (lbfgs :m 17)
    :until (any (grad-norm 1e-5) (e-value 20) (budget-exhausted))
    :emit (pareto ledger report)))
```

Every operator lives in a **machine-readable catalog**: full type signature, units, cost model, error model, determinism class, capability requirements, worked examples, and semantic version. `catalog(query)` is how an agent discovers the system; it is generated from the code, so it cannot drift. High-level verbs (`optimize_shape`, `simulate_pour`) exist for terseness but **lower to explicit IR**, and the lowering is inspectable — progressive disclosure with nothing hidden.

### 11.2 The Design Ledger [S]

FrankenSQLite is the system of record (Bet 10). One file per project; WAL mode; snapshot-isolated readers so dashboards and agents query freely while a 96-core sweep writes. Six core tables (DDL in Appendix D):

- **artifacts** — content-addressed (in-house BLAKE3-class tree hash) typed blobs: charts, fields, meshes, fronts, images. Identical results dedupe to one row, across forks, forever.
- **ops** — the event-sourced log: every executed IR node with its full inputs (seeds, versions, budgets, capabilities — the Five Explicits, frozen), timing, and outcome.
- **edges** — the lineage DAG connecting ops to the artifacts they consumed and produced.
- **metrics** — time series: residuals, e-values, objective traces, ledger balances, roofline attainment.
- **tune** — the autotuner's cache (§5.5), keyed by kernel × shape-class × machine fingerprint.
- **events** — the fine-grained session stream (tile completions, cancellations, preemptions) that powers live dashboards and post-hoc flame graphs.

Because sessions are event-sourced over content-addressed state, three expensive-sounding features are nearly free: **time travel** (`at(t)` replays the log to any instant), **forkable worlds** (a fork is a new branch of the op log sharing every artifact by hash — an agent speculatively explores three formulations for the cost of the deltas), and **`explain(artifact)`** (walk `edges` backward and render the full causal tree: which chart, which mesh, which seed, which solver tolerances, which library versions produced this number). A six-month campaign is a database you query, not a directory you fear.

### 11.3 Sessions, capabilities, and the resource governor [S/F]

Every IR program executes inside a **session** holding a **capability token**: an explicit grant of operators, core-seconds, resident memory, wall time, and ledger scope. Budgets are *enforced*, not advisory — the resource governor does admission control onto the two-lane executor, meters per-study core quotas, and on memory pressure degrades gracefully in a declared order: spill coldest arenas, coarsen adaptively, or pause–serialize–resume (P7 makes every solver checkpointable, so "pause the LES, run the urgent trim study, resume" is routine, not heroic).

Agent-proofing details that matter in practice: **idempotency keys** on every mutating call, so a retried request cannot double-spend a budget or duplicate a study; **dry runs** — `estimate(program)` returns predicted wall time, memory, error, and energy from the cost models *without executing*, so agents plan before they spend; and **errors as guidance** — every failure is a structured value carrying a machine-readable diagnosis plus candidate fixes ranked by the cost model:

```json
{"error":"BudgetInfeasible","stage":"flux.lbm","need":{"wall":"5.1h"},"have":{"wall":"2h"},
 "fixes":[{"action":"relax qoi-rel-error to 4e-2","est_wall":"1.7h","est_qoi_impact":"+1.8e-2"},
          {"action":"surrogate screen (§9.7), certify top-4 only","est_wall":"1.9h"},
          {"action":"coarsen lattice to dx=1.4mm outside lip band","est_wall":"1.8h"}]}
```

A refusal that teaches is worth ten silent successes.

### 11.4 The Error Ledger, the Time Ledger, and the self-optimizing pipeline (Bet 12) [F/M]

Every operator publishes two small models: an **error model** (a-priori rates plus a-posteriori estimators — DWR from §8, interval widths from `fs-ivl`, conformal bands from §9.7, MC half-widths) and a **cost model** (quantile-regression fits over the `tune` table for *this* machine). Composition of operators composes the models: the **Error Ledger** attributes the end-to-end uncertainty of a quantity of interest to its sources — geometry tolerance, discretization, algebraic residual, surrogate band, statistical noise — and the **Time Ledger** does the same for seconds.

That turns budget allocation into an optimization problem HELM solves *about itself*: given "drag to 2% in 2 hours," the planner allocates error and time across mesh size, polynomial order, solver tolerance, surrogate usage, and sample counts to maximize expected accuracy-per-second — convex where the models are (rate-based bounds usually are), CMA-ES on itself where they aren't, re-planned online as a-posteriori estimates replace a-priori guesses. Version 1 ships greedy-plus-lookup [F]; the full co-optimizer is an [M] bet promoted only after the Gauntlet shows its plans beat hand allocation on the flagships. Scheduling analytics for the same planner live in §14.3.

### 11.5 Agent ergonomics, the last mile [S/F]

The **Five Explicits** (units, seeds, budgets, versions, capabilities) are mandatory fields of every op — an agent literally cannot express an ambiguous request. Compile-time dimensional analysis via `Qty<const M:i8, const KG:i8, const S:i8, …>` (Appendix B) extends into the IR checker so unit errors die at admission. **Semantic design diffs**: `diff(a, b)` between two shapes returns not a file diff but a geometric one — varifold/OT distance, transport plan visualization, per-region attribution ("Pareto-31 moved lip curvature −18%, thinned wall 0.4mm") using the §9.9 machinery. **Automatic lab notebook**: every study emits, via FrankenPandas over the ledger, a human-readable report — objective traces, ledger attributions, LUMEN renders, convergence tables, the exact IR to reproduce — because reproducibility should be a side effect, not a virtue. And everything conversational (catalog, estimates, ledger queries, `explain`) is served from the latency lane in ≤100ms even while a campaign saturates the throughput lane — the agent never waits behind its own physics.

---

## 12. The constellation, deepened: contracts with the Franken libraries [S]

§4.1 assigned roles; this section states the *contracts* — what FrankenSim demands of each library, and what flows back upstream, since FrankenSim will be the constellation's most demanding consumer.

- **asupersync** — Demands: structured-concurrency scopes with sub-millisecond latency-to-cancel at tile granularity, cancellation trees deep enough for speculative races over solver hierarchies, and stable semantics under the G4 chaos harness. Returns upstream: the executor's model-checking harness, cancellation-latency benchmarks, and the two-lane scheduler as a reusable pattern.
- **FrankenSQLite** — Demands: WAL crash safety under append-heavy event logging with concurrent snapshot readers, multi-GiB blob handling for content-addressed artifacts, and predictable checkpoint behavior mid-sweep. Returns: a punishing OLTP+blob stress suite and the ledger schema as a reference workload.
- **FrankenNumpy** — Contract: zero-copy interop at the membrane. Cochain and lattice storage expose FrankenNumpy-compatible strided views so agents drop to raw arrays without copies; FrankenNumpy ufuncs run on the throughput lane.
- **FrankenTorch** — Contract: its tape is the reverse-mode backend for surrogate training (§9.7) and neural-implicit charts (§7.2). FrankenSim contributes custom differentiable ops upstream — "adjoint linear solve," "implicit-function-theorem node," "conformal calibration layer" — turning solver adjoints into first-class autograd citizens.
- **FrankenScipy** — Contract: the conformance oracle. Where semantics overlap (special functions, reference optimizers, sparse baselines), FrankenSim's kernels must match FrankenScipy in the Gauntlet and beat it on the roofline; disagreements become bug reports for one side or the other.
- **FrankenPandas** — Contract: first-class frames over every ledger table; powers §11.5 notebooks and §14.4 performance analytics.
- **FrankenNetworkx** — Contract: graph algorithms at scale for ground structures (§9.5), mesh/partitioning experiments, constraint propagation, and the max-plus critical-path analysis (§14.3); FrankenSim contributes the tropical-semiring iteration kernels.

One **constellation lockfile** pins exact versions per project; the lock hash is part of every op's Five Explicits, so `explain()` can always answer "which FrankenTorch produced this gradient."

---

## 13. The Gauntlet: the correctness program

A system that sells *certified* error had better be paranoid about its own. The Gauntlet is a permanent, tiered verification program; nothing tagged [F] or [M] ships to the default path without passing every tier below it, on both target ISAs, nightly.

### 13.1 The six tiers

- **G0 — Properties and laws.** In-house property-based testing (proptest-class shrinking) of algebraic laws: chart-conversion functor laws (round-trips within certified bounds), adjoint consistency ⟨Av, w⟩ = ⟨v, Aᵀw⟩ to roundoff for every operator with an adjoint, GA identities, exact-sequence identities (δδ = 0, div∘curl = 0 *bitwise* on discrete cochains), conservation checks (mass in free-surface LBM, energy drift bounds for symplectic integrators), interval-arithmetic containment under random rewrites.
- **G1 — Manufactured solutions and order verification.** Every PDE frontend × element family × BC type carries an MMS battery; CI fits the observed convergence order and **fails the build if the slope deviates from theory by more than 0.2**. Adjoint order is verified the same way. This is the tier that catches the quiet death of accuracy.
- **G2 — Canonical benchmarks.** Cook's membrane; lid-driven cavity at Re 1000 against Ghia et al.; Taylor–Green vortex decay; circular-cylinder Cd/Strouhal at Re 100; NAFEMS shell and locking tests; plate buckling loads; NACA 0012 polars against experiment (with LES honesty labels); an El Centro-class frame time-history against published fiber-model results. Each stores its acceptance envelope in the repo, results in the ledger.
- **G3 — Metamorphic tests.** Frame invariance (rigidly move every input; outputs must be invariant to tolerance), unit-rescaling invariance through runtime paths (compile-time `Qty` can't cover data files), estimator monotonicity under refinement, adjoint-vs-finite-difference agreement along random directions at randomized points, Rep-Router path independence (any admissible conversion path must land within the composed certificate).
- **G4 — Chaos and cancellation storms.** Kill random scopes at random times during real solves; assert zero leaks via arena accounting, no deadlocks, ledger consistency, and that pause–serialize–resume reproduces the uninterrupted trajectory in deterministic mode. Fault-inject tiles, poison allocators, torture the two-lane executor with adversarial preemption; run asupersync's model checker over the executor's core protocols. Cancellation-correctness (P7) is a *claim*; G4 is the proof engine.
- **G5 — Determinism audits.** Bit-identical results across runs and thread counts on the same ISA in deterministic mode; a cross-ISA divergence report classifying every difference (FMA contraction, libm ULP, reduction shape) with documented envelopes. Artifact hashes in the ledger make silent drift impossible.

### 13.2 Certifying the certifiers [F/M]

The exotic machinery gets its own trials: interval/Taylor-model enclosures spot-checked against a slow in-house high-precision (double-double/quad-double) oracle; **e-process validity verified empirically** — under simulated nulls, the anytime-valid type-I error at *adversarially chosen* stopping times must stay ≤ α; conformal e-prediction coverage tested under distribution shift scenarios; the sheaf watertightness certificate (§7.3) cross-examined by brute-force ray-parity queries on adversarial multi-chart seams. If a certificate can be fooled, that is a Sev-0 bug — worse than a wrong answer, it is a wrong answer wearing a badge.

### 13.3 Conformance as the inter-agent contract [S]

Every crate ships `CONTRACT.md` (semantics, invariants, error models, determinism class) plus an executable conformance suite; "done" is defined as *conformance green on both ISAs*, nothing else. Because the suites speak FrankenScript IR, they double as the integration language between the AI agents building the system (§16.2): an agent implementing `fs-fft` never negotiates with the agent implementing `fs-lbm` — both negotiate with the contract.

---

## 14. The performance program

Performance claims in design documents are usually vibes. Ours are rooflines with names attached, re-baselined continuously by the autotuner on the actual machines (§5.5); the numbers below are **illustrative targets** for the two reference machines, stated so they can be *failed*.

### 14.1 Roofline targets [S]

| Kernel | Arithmetic intensity | M4 Max target | TR PRO 96c target |
|---|---|---|---|
| LBM D3Q19 stream+collide (sparse tiles) | ~152 B/LUP, bandwidth-bound | ≥ 1.0 GLUP/s | ≥ 0.6 GLUP/s (≥ 1.0 on 12-CCD EPYC-class) |
| GEMM f64 (fs-la, BLIS-style) | compute-bound | ≥ 75% of ~0.5 TF NEON peak (SME2 upside exploratory) | ≥ 75% of ~6–7 TF all-core AVX-512 peak |
| SpMV / SELL-C-σ | bandwidth-bound | ≥ 85% of STREAM | ≥ 85% of STREAM (per-CCD placement) |
| Matrix-free FEEC apply, p = 4 (sum-factorized) | ~mixed | ≥ 30% of peak FLOPs | ≥ 30% of peak FLOPs |
| Batched 8×8–32×32 dense (SIMD-across-elements) | compute-bound | ≥ 60% peak | ≥ 60% peak |
| Stockham FFT 3D, in-cache pencils | mixed | ≥ 40% of memory-bound limit | ≥ 40% |
| Sphere-traced SDF primary rays | latency-bound | ≥ 80 Mray/s | ≥ 120 Mray/s |

Two honest caveats: the M4's 546 GB/s unified memory makes it a *bandwidth monster per core* (LBM and SpMV love it), while the Threadripper's strength is aggregate compute and cache — so kernel-level winners flip by workload, and the Rep Router's cost models are machine-specific on purpose. And every target is stored in the ledger with the machine fingerprint; a "target" that was never re-measured is a lie waiting to happen.

### 14.2 Machine-adaptive execution [S/F]

The `tune` table keys on kernel × shape-class × machine fingerprint (ISA, core topology, cache geometry, measured bandwidths). On Threadripper: CCD-aware tile ownership so working sets stay inside their 32 MiB L3 island, first-touch page placement, and reduction trees shaped to the CCD hierarchy. On M4: asymmetric scheduling that keeps latency-lane work on P-cores, streams throughput tiles across P+E with E-core-sized tile variants, and (exploratory, [F]) SME2 streaming-mode GEMM. Hugepage-backed arenas everywhere; software prefetch distances autotuned per machine; SIMD tier chosen once at startup.

### 14.3 Tropical critical-path analytics (Bet 12) [F/M]

Model the task DAG's timing in the **max-plus (tropical) semiring**: composition of stage delays is ⊕ = max, ⊗ = +; steady-state throughput of the pipelined DAG is the tropical eigenvalue of its timing matrix, and the critical circuit is its eigenvector's support. FrankenNetworkx runs the tropical power iteration over the live event stream, and the output is *actionable*: the makespan's sensitivity to each kernel's latency tells the autotuner **which kernel to tune next**, and tells the §11.4 planner where a relaxed tolerance actually buys wall-clock (slack elsewhere buys nothing). Scheduling theory as an online instrument, not a whiteboard aphorism.

### 14.4 Performance regressions are test failures [S]

Nightly roofline runs on both reference machines write to `metrics`; any kernel dropping below its tolerance band fails CI exactly like a wrong answer. The event stream reconstructs flame-graph-equivalents post hoc, so a regression arrives with its own diagnosis. FrankenPandas dashboards over the ledger make "what got slower this month, and why" a one-liner.

---

## 15. Three flagship pipelines, end to end

The flagships are not demos; they are the forcing functions that keep every layer honest. Each is written as the IR-level story an agent would actually run, with an explicit "what breaks first" clause.

### 15.1 Ornithoid multi-inlet aircraft: L/D × stability × maneuverability

1. **Parameterize.** Baseline bird-like body as an F-rep (fuselage blend + wing surfaces as IGA shells); inlets as smooth R-function Booleans whose positions/lip shapes are design variables; global form deformed by manifold harmonics (~200 coefficients) + local FFD near inlets. All levers expose Jacobian actions (§7.6).
2. **Screen wide.** BEM panel method with Kutta condition + kernel-independent FMM for induced velocities; vortex-particle wake shed from trailing edges for *flapping* gaits (§8.3); hundreds of candidates per generation raced under **e-raced generations** (§9.6) — dominated wings die at 10–20% of their solve and their cores are reassigned within a millisecond (P7).
3. **Refine survivors.** Cumulant LBM/LES on sparse VDB lattices with SDF-sampled curved boundaries; forces via momentum exchange; LES honesty labels attached (§8.3); DWR adaptivity keyed to the L/D functional (§8.6).
4. **Trim and certify stability.** Fit a reduced flight-dynamics model (Koopman/DMD surrogate with conformal e-bands, §9.7) around each trim state; then **SOS Lyapunov synthesis** (§9.8) produces a certified region of attraction — "stable" becomes a proven basin volume, not an eigenvalue vibe. Maneuverability scored as certified reachable-set proxies under bounded control effort.
5. **Pareto.** NSGA-III over (L/D, ROA volume, maneuver metric, inlet mass-flow satisfaction) with gradient polish (adjoint LBM being deferred [M], polish uses BEM-level adjoints + surrogate gradients, honestly labeled); deliverable is a **certified Pareto atlas**: every point carries its Error Ledger, its stability certificate, and `explain()` lineage.

*What breaks first:* LES cost per candidate — mitigated by the multi-fidelity stack (panel → surrogate → LBM) and e-racing; if SOS scaling disappoints, ROA falls back to sampled Lyapunov with conformal bounds [F].

### 15.2 Seismic-minimal building frame: minimum material, certified fragility

1. **Layout.** FrankenNetworkx ground structure (millions of candidate members under fabrication rules); member-force layout LP/SOCP via the in-house PDHG solver (§9.5); Michell-continuum limits as analytic sanity oracles on subproblems.
2. **Size nonlinearly.** Trust-region Newton–Krylov sizing with Euler + code-based buckling constraints, joint-cost regularization; arclength continuation spot-checks global buckling modes (§8.2).
3. **Nonlinear time history.** Fiber-section beam-columns (Mander concrete, Menegotto–Pinto steel) under suites of recorded motions + **Kanai–Tajimi** stochastic ground-motion ensembles (§8.5); variational integrators keep long-duration response drift-free.
4. **Quantify risk anytime-validly.** MLMC over motion ensembles with **e-process stopping** on P(peak drift > limit): the study stops itself the moment the fragility estimate is decision-grade at the requested confidence — valid at that data-dependent stopping time by construction (Bet 5).
5. **Optimize under tail risk.** CVaR-constrained mass minimization (Rockafellar–Uryasev smooth form, §9.9); constructability pass snaps members to catalog sections and enforces angle/joint rules as graph constraints.

Deliverable: frame + fragility curves + an **anytime-valid certificate** that the collapse-proxy probability meets target. *What breaks first:* nonlinear time-history volume — mitigated by MLMC level design, surrogate screening of motion suites, and the fiber-beam kernels' batched-small-dense hot loop (§6.1) which is exactly what the substrate is tuned for.

### 15.3 Laminar-pour vessel: the spout that never dribbles

1. **Parameterize.** F-rep vessel of revolution + free-form lip; spout region driven by level-set velocity parameterization (§7.6) — topology of the lip channel may change and that is fine.
2. **Stability objective, cheaply.** Quasi-steady thin-film flow along the lip streamlines; **Orr–Sommerfeld eigenproblems solved by `fs-cheb`** spectral collocation (§6.5) give modal growth rates σ₁…σ₈ along the pour path — the objective is min–max certified modal growth, a differentiable proxy for "laminar."
3. **Validate with physics.** Free-surface cumulant LBM pours across the tilt schedule and viscosity band (Carreau parameters spanning the target fluids); a Plateau–Rayleigh breakup model scores the free jet; spill/dribble indicator from contact-line behavior on the lip SDF.
4. **Robustify.** CVaR over the viscosity/rate band so the spout is laminar for *the family* of liquids, not the nominal one; e-raced candidates under the LBM validator.
5. **Deliver.** Optimized vessel + LUMEN spectral render of the pour (Woodcock-tracked LBM density field, §10.1) — the marketing shot and the physics are the same bytes.

*What breaks first:* contact-line modeling fidelity (a genuinely open problem) — handled by bracketing contact-angle models and reporting the objective's sensitivity band in the Error Ledger rather than pretending certainty.

---

## 16. Roadmap: phases, exit criteria, and how it's actually built

### 16.1 Phases [S]

| Phase | Weeks | Scope | Exit criterion (Gauntlet-enforced) |
|---|---|---|---|
| P0 — Bedrock | 0–6 | fs-substrate, fs-exec (two-lane), fs-alloc, fs-la, fs-sparse, fs-fft, fs-ivl, fs-rand, ledger v0 | G0+G4 green; GEMM/SpMV/FFT within 80% of §14.1 targets on both ISAs; deterministic mode bit-stable |
| P1 — Geometry + eyes | 6–14 | fs-geom, SDF/F-rep/mesh charts, Rep Router v1, dual contouring, Delaunay tet, LUMEN preview tracer | Chart round-trips certified; watertightness certificate vs ray-parity oracle; sphere-traced turntable of an F-rep at target ray rates |
| P2 — Elasticity + first optimization | 14–24 | FEEC elasticity, CutFEM-on-SDF, matrix-free p-MG + AMG, fs-ad adjoints, SIMP + Helmholtz filter | **Marquee demo: topology optimization on a raw SDF — no mesh in the loop — with a composed error certificate**; G1 orders verified; adjoint FD-checks green |
| P3 — Fluids I | 24–34 | fs-lbm (cumulant, sparse, free-surface), lattice-scaling assistant, thermal, non-Newtonian | Cavity/TGV/cylinder G2 green; ≥ targets on GLUP/s; first spout pours end-to-end |
| P4 — Structures at scale | 34–44 | IGA + Kirchhoff–Love shells, fiber beams, ground-structure PDHG, Kanai–Tajimi + MLMC + e-stop | Frame flagship v1: fragility with anytime-valid stopping; NAFEMS shell suite green |
| P5 — Aero stack | 44–56 | BEM+FMM+Kutta, vortex particles, Dirac coupling, SE(3) integrators, Koopman surrogates | Ornithopter flagship v1: screened→refined Pareto with e-raced generations live |
| P6 — Certificates & self-optimization | 56–68 | SOS/Lasserre SDP, sheaf certificates promoted, conformal e-prediction hardened, §11.4 planner, diff-rendering | [M] features pass §13.2 or ship flagged-off; planner beats hand-allocated budgets on all three flagships |

Phases overlap deliberately (geometry hardens while physics starts); each phase gate is a Gauntlet state, not a date. Nothing [M] gates anything [S].

### 16.2 Team-of-agents methodology [F]

FrankenSim is sized for a swarm of AI coding agents with human architectural review. The load-bearing practices: **one crate = one contract** (`CONTRACT.md` + executable conformance suite, §13.3); **IR as the integration language** — agents integrate against frozen IR semantics, never against each other's internals; **golden ledgers** — every merged feature lands with a replayable ledger of its acceptance run; **the Decalogue as tie-breaker** — disputes between agents resolve by principle number, not seniority. The repository is organized so the maximum context an agent needs is one crate + its contracts + the IR spec — deliberately smaller than a frontier context window.

---

## 17. Risk register

| Risk | Sev | Mitigation |
|---|---|---|
| Pure-Rust kernels miss BLAS-class performance | High | BLIS-style microkernel discipline + autotuner + §14.4 CI; targets stated to be *failable*; batched-small-dense (the actual hot loop of FEM) favors our design over generic BLAS anyway |
| SME2/streaming intrinsics immaturity on M4 | Med | NEON path is the committed baseline; SME2 is flagged [F] exploratory upside only |
| Hex meshing is famously hard | Med | IGA + CutFEM cover most hex use-cases better; frame-field hexing stays [F], never on critical path |
| Trimmed-NURBS watertightness | Med | Booleans route through SDF + certified re-fit (§7.2); sheaf certificate catches seam failures |
| LES cost explodes optimization budgets | High | Multi-fidelity ladder (panel → surrogate → LBM), e-racing, DWR adaptivity, honesty labels prevent silent under-resolution |
| First-order SDP (SOS) scaling/conditioning | Med | Scope to ≤ ~10 design vars; Burer–Monteiro + dual certificates; fall back to sampled Lyapunov with conformal bounds |
| No LBM adjoint (yet) | Med | Honest [M] deferral; gradients via BEM-level adjoints + certified surrogates; black-box lane (CMA-ES) always available |
| [M] math (sheaf/e-race/planner) fails validation | Med | All [M] behind flags off the critical path; §13.2 decides; the [S] spine ships regardless |
| Single-machine bias blocks future scale-out | Low | Serializable solver states + content-addressed artifacts + event-sourced sessions are exactly the primitives distribution needs later |
| Scope creep (this document is the evidence) | High | Phase gates + Decalogue + ambition tags; anything new must displace something or wait |

---

## 18. Appendices

### Appendix A — The crate atlas

`fs-substrate` (arch detect, fingerprints, dispatch tables) · `fs-simd` (portable/NEON/AVX-512/SME2 tiers) · `fs-alloc` (scope arenas, hugepages, pools) · `fs-exec` (two-lane executor, speculative races, deterministic reductions, on asupersync) · `fs-la` (GEMM, batched small dense, TSQR, LOBPCG/Lanczos, mixed precision, Hutch++) · `fs-sparse` (BSR/SELL-C-σ, SpMV/SpMM, AMG) · `fs-fft` (Stockham, real transforms, pencils) · `fs-ivl` (interval/affine/Taylor models, exact predicates, double-double oracle) · `fs-cheb` (function objects, spectral 1D–3D, eigenproblems) · `fs-ad` (duals, tape bridge to FrankenTorch, adjoint-by-construction, revolve checkpointing) · `fs-rand` (Philox, Sobol/Owen, distributions) · `fs-ga` (PGA/CGA, const-evaluated ops) · `fs-geom` (Regions, charts, Rep Router, sheaf certificates) · `fs-rep-sdf` / `fs-rep-frep` / `fs-rep-nurbs` / `fs-rep-mesh` / `fs-rep-voxel` / `fs-rep-neural` · `fs-xform` (FFD, RBF, level-set velocity, manifold harmonics, SIMP fields) · `fs-mesh` (dual contouring, Delaunay/Ruppert, anisotropic remesh, frame-field hex [F]) · `fs-feec` (Whitney/Nédélec/RT families, sum-factorized apply, DWR) · `fs-cutfem` (ghost penalty, cut quadrature) · `fs-iga` (splines-as-elements, shells) · `fs-solid` (hyperelastic, TDNNS, rods, fiber beams, contact/IPC) · `fs-lbm` (cumulant, sparse VDB, free surface, thermal, scaling assistant) · `fs-bem` / `fs-fmm` (panels, Kutta, bbFMM) · `fs-vpm` (vortex particles) · `fs-time` (variational/symplectic/SE(3), arclength) · `fs-couple` (Dirac ports, IQN-ILS) · `fs-adjoint` (discrete adjoints, shape gradients, Sobolev smoothing) · `fs-uq` (KL, PCE, QMC, MLMC, Kanai–Tajimi, CVaR) · `fs-eproc` (e-processes, e-BH, conformal e-prediction) · `fs-opt` (L-BFGS, TR-Newton–Krylov, AL, CMA-ES/BIPOP, DIRECT, GP/TuRBO, NSGA, PDHG LP/SOCP) · `fs-topo` (SIMP, level-set, ground structure, homogenization, persistence constraints) · `fs-sos` (moment/SOS, Burer–Monteiro SDP, Lyapunov) · `fs-surrogate` (FNO/DeepONet/POD-DEIM/Koopman, certify-or-escalate) · `fs-render` (spectral PT, sphere/NURBS tracing, diff-render) · `fs-viz` (isosurface, DVR, LIC, glyphs, Morse–Smale) · `fs-img` (PNG/EXR) · `fs-ir` (FrankenScript, catalog, admission) · `fs-session` (capabilities, governor, idempotency) · `fs-ledger` (schema, content addressing, time travel, explain) · `fs-plan` (cost/error models, budget allocator, tropical analytics) · `fs-report` (lab notebooks via FrankenPandas).

### Appendix B — Load-bearing trait sketches

```rust
/// A chart presents a Region through one representation. (§7.1)
pub trait Chart: Send + Sync {
    type Param;                      // design-lever handle
    fn eval(&self, x: Point3, cx: &Cx) -> ChartSample;      // value + gradient + certified Lipschitz data
    fn support(&self) -> Aabb;
    fn topology_hint(&self) -> BettiBounds;
}

/// Certified conversion between representations: a functor with a receipt. (§7.3)
pub trait Convert<Dst: Chart>: Chart + Sized {
    fn convert(&self, budget: ErrBudget, cx: &Cx) -> Result<Certified<Dst>, Diag>;
    // Certified<T> = value + interval error bound + provenance hash + adjoint hook
}

/// Compile-time dimensional analysis: the Five Explicits, unit edition. (§11.5)
#[derive(Clone, Copy)]
pub struct Qty<const M: i8, const KG: i8, const S: i8, const K: i8, const A: i8>(pub f64);
pub type Length   = Qty<1, 0, 0, 0, 0>;
pub type Stress   = Qty<-1, 1, -2, 0, 0>;   // Pa = kg·m⁻¹·s⁻²
impl<const M: i8, const KG: i8, const S: i8, const K: i8, const A: i8>
    core::ops::Add for Qty<M, KG, S, K, A> { /* same dims or it does not compile */ }

/// Every kernel is a tile program under a context. (§5.2; TileKernel shown there.)
pub struct Cx<'scope> {
    pub cancel: CancelToken<'scope>,   // asupersync scope; poll at tile edges (P7)
    pub arena:  &'scope Arena,
    pub rng:    PhiloxKey,             // keyed by logical identity, not thread (P2)
    pub budget: BudgetSlice,
    pub ledger: LedgerHandle,
}
```

### Appendix C — A fuller FrankenScript study (seismic frame)

```lisp
(study "frame-seismic-cvar-v9"
  (seed 0xF00D0002) (versions (constellation :lock "2026-07"))
  (capability :cores 96 :mem 384GiB :wall 36h :ops (flux.* ascent.* uq.*))
  (budget (qoi "P(drift>2e-2)" :rel-error 0.15 :confidence 0.95))
  (let site   (uq.ground-motion (kanai-tajimi :S0 0.03m2/s3 :wg 15rad/s :zg 0.6)
                                (records "PEER-set-A") (mlmc :levels 4)))
  (let ground (topo.ground-structure (grid 8 x 5 x 24m) :knn 14 :rules "AISC-cat.json"))
  (let layout (ascent.solve-lp (min (member-volume ground)) :method pdhg
                               :oracle (michell :tol 0.08)))
  (let frame  (topo.size layout :method tr-newton-krylov
                :constraints ((buckling :code "AISC-E3") (drift-elastic 5e-3))))
  (let resp   (flux.fiber-frame frame site :integrator variational :dt-adapt true))
  (let frag   (uq.probability (exceeds (peak-drift resp) 2e-2)
                :stop (e-process :alpha 0.05)))          ; anytime-valid (Bet 5)
  (ascent.optimize (min (mass frame)) :over (sections frame)
    :subject-to ((cvar frag :beta 0.9 :le 0.02) (constructable :catalog "AISC"))
    :method augmented-lagrangian
    :emit (frame frag report ledger)))
```

### Appendix D — Design Ledger schema (v0)

```sql
CREATE TABLE artifacts(hash BLOB PRIMARY KEY, kind TEXT NOT NULL, bytes BLOB,
                       meta JSON, created_at INTEGER) STRICT;
CREATE TABLE ops(id INTEGER PRIMARY KEY, session BLOB, ir JSON NOT NULL,
                 seed BLOB, versions JSON, budget JSON, capability JSON,
                 t_start INTEGER, t_end INTEGER, outcome TEXT, diag JSON) STRICT;
CREATE TABLE edges(op INTEGER REFERENCES ops(id), artifact BLOB REFERENCES artifacts(hash),
                   role TEXT CHECK(role IN ('in','out')), PRIMARY KEY(op, artifact, role)) STRICT;
CREATE TABLE metrics(op INTEGER, t INTEGER, name TEXT, value REAL,
                     PRIMARY KEY(op, t, name)) STRICT;
CREATE TABLE tune(kernel TEXT, shape_class TEXT, machine BLOB, params JSON,
                  measured JSON, PRIMARY KEY(kernel, shape_class, machine)) STRICT;
CREATE TABLE events(id INTEGER PRIMARY KEY, session BLOB, t INTEGER,
                    kind TEXT, payload JSON) STRICT;
```

### Appendix E — Reading list, keyed to the bets

Arnold–Falk–Winther, *Finite Element Exterior Calculus* (Bet 2) · Hansen–Ghrist, cellular sheaves and sheaf Laplacians (Bet 3) · Ramdas et al., game-theoretic statistics and e-processes (Bet 5) · Vovk et al., conformal e-prediction (Bet 6) · Waudby-Smith–Ramdas, betting confidence sequences (Bet 5) · Ollivier et al., information-geometric optimization (Bet 7: CMA-ES as natural gradient) · Hughes et al., isogeometric analysis (Bet 2) · Burman et al., CutFEM/ghost penalty (Bet 1) · Barill et al., generalized winding numbers (§7.2) · Shewchuk, robust adaptive predicates (§6.4) · Allaire–Jouve–Toader, level-set shape/topology optimization (§9.5) · He–Gilbert, rationalized layout optimization (§9.5) · Marsden–West, variational integrators (Bet 2) · van der Schaft–Jeltsema, port-Hamiltonian systems (Bet 4) · Lasserre, moment-SOS hierarchies (Bet 8) · Peyré–Cuturi, computational optimal transport (§9.9) · Michor–Mumford, Riemannian shape spaces (§7.6) · Li et al., Incremental Potential Contact (§8.2) · Kiendl et al., Kirchhoff–Love IGA shells (§8.2) · Geier et al., cumulant LBM (§8.3) · Körner et al., free-surface LBM (§8.3) · Ying–Biros–Zorin, kernel-independent FMM (§8.3) · Eriksson et al., TuRBO trust-region BO (§9.4) · Meyer–Musco–Musco–Woodruff, Hutch++ (§6.1) · Halko–Martinsson–Tropp, randomized numerical linear algebra (§6.1) · Griewank–Walther, *Evaluating Derivatives* + revolve (§6.6) · Dorier et al./Trefethen, chebfun-style computing with functions (§6.5) · Ju–Losasso–Schaefer–Warren, dual contouring (§7.5) · Nieser–Reitebuch–Polthier and Ray et al., frame fields for hexing (§7.5) · Baydin et al., AD survey (§6.6).

---

## Coda

The archipelago dies the day derivatives, error bounds, budgets, provenance, and cancellation ride *inside the values* instead of living in six incompatible tools' heads. Everything above is in service of that one move: a single typed continuum from "describe a region" to "here is the certified Pareto atlas," saturating whatever cores it is given, honest about every approximation it makes, and legible — natively — to the agents who will do most of the designing. Build P0. The rest is compounding.
