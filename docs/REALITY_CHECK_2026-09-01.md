# FrankenSim Reality Check — 2026-09-01

HEAD `b46da547`. Evidence gathered by twelve read-only investigation lanes plus one executed compile/test run on the yto verification host. Code is the ground truth; the five plan documents (880 KB) and the README are the measuring stick.

---

## 0. The one-paragraph answer

FrankenSim is a very large, mostly honest library substrate: 1.33 M lines of Rust source and 0.57 M lines of tests across 171 `fs-*` crate directories, 5,433 commits in 58 days, 171 of 171 contracts present, layer direction enforced. Its plan documents promise a continuum — typed intent → geometry → certified physics → adjoints → optimization → render → ledger — plus a product wedge (Cooling 0.1), a browser Wright Flyer, musical-instrument physics, an Euler-disc flagship, five machine-domain flagships, and a fifteen-initiative "leapfrog" program. What exists end to end today: **one** production-grade 3-D physics solver (steady conduction) reachable through a CLI that stops before any user-visible answer; **one** shipped product (the Wright Flyer, live, real wasm, reduced physics, self-referential validation); a library of instrument physics nobody has listened to; a real spectral path tracer that has produced a 32×18-pixel film; solid L0/L1 numerics that miss the plan's own performance bars on every measured kernel; and a bead graph of 3,289 items that is 83 % blocked, has 2 ready items, and was seeded with ~2,200 beads in two July weeks that target physics that does not exist yet. The central thesis — a design study running from typed intent through certified physics to an optimized geometry — **has never executed once on a production PDE.** Libraries compile at HEAD; `cargo test --workspace` cannot build (four test targets committed red between Aug 23 and Aug 27).

---

## 1. Headline numbers

| Dimension | Value | Note |
|---|---|---|
| Vision documents | 5 docs, 880 KB, ~170 extracted testable goals | main plan + addendum + new-domains + Wright Flyer + music |
| Source | 1,329,114 src lines; 568,622 test lines; 171 crate dirs | 3 of the 9 "standalone fs-* workspaces" contain zero FrankenSim code |
| Commits | 5,433 since 2026-07-05; 2,216 in the last 30 days (995 touch source) | all authored as the owner; swarm output |
| Beads | 3,289 total · 1,699 open · 204 in_progress · 1,300 closed · **2 ready** · 303 open epics · 368 open P0 | 152 of 204 in_progress untouched > 14 days (103 > 30 days) |
| Bead creation | W28 452 · **W29 1,522** · W30 223 · **W31 693** · W32–W36 399 | closure runs ~150/week; the two July waves alone exceed everything closed to date |
| Build at HEAD | libs/bins ✓ (162/163 members checked on yto) · **workspace test build ✗** (4 red test targets) | fs-scenario, fs-feec, fs-session, fs-evidence-runner |
| xtask gates at HEAD | 7 OK · 5 violating | check-deps 2, check-unsafe 1, check-suite-receipt 4, check-program-metrics 2, check-source-manifest 2 |
| Suite receipt | 2026-08-01, HEAD 88b18883: 7,621 pass / 95 unexpected red / 38 ignored | 31 days and ~2,200 commits stale |
| Capability registry | L1 = 3 · L2 = 11 · **L3 = 1** · L4 = 0 · L5 = 0 | the L3 entry is not backed by its cited lane (see §4.1) |
| Product CLI | `solve` executes 4 of 6 stages (5 with a declared conduction section); QoI · report · package · run · compare all refuse | the QoI bead is closed while the code still refuses |
| Constellation | all 7 siblings off-pin on the Mac (fast-forwards, 25–288 commits ahead); committed Cargo.lock records drifted versions | Mac builds are not against the declared constellation |

---

## 2. Vision checklist (core plan + addendum), with status

Status vocabulary: WORKING · PARTIAL · STUB · UNPROVEN · NOT_STARTED · REGRESSED · WRONG_APPROACH. Bead-coverage column counts open/in-progress beads whose title or description mentions the item (grep over 1,902 live beads).

### 2.1 Architecture and substrate (L0/L1)

| # | Goal (plan source) | Status | Evidence | Beads |
|---|---|---|---|---|
| 1 | Flat acyclic `fs-*` workspace, CONTRACT per crate, layer direction enforced (MP §4, §13.3) | WORKING | check-layers 0, check-contracts 0, 171/171 | n/a |
| 2 | Franken-only runtime deps; unsafe confined to audited capsules (P1) | WORKING, 3 violations | check-deps: fs-g1-train ships serde/serde_json at runtime; check-unsafe: `fs-simd/src/wasm/mod.rs` unregistered | 0 |
| 3 | Deterministic mode bit-identical across runs/threads on one ISA (P2, G5) | PARTIAL | G5 cross-ISA report clean on 6 artifacts; 27 goldens in golden-couplings.json; music/WF/audio chains one-host only | many |
| 4 | Two-lane executor, tile cancellation, ≤ 200 µs latency-to-cancel (§5.2) | PARTIAL | executor, Cx, gates, races exist; fs-exec CONTRACT: "NO 200 µs cancel-latency CLAIM yet" | 45 |
| 5 | G0 + G4 green (P0 exit) | PARTIAL | G0 harness closed; G4 chaos harness bead 6nb.5 in_progress since 07-17 | 4 |
| 6 | GEMM ≥ 75 % peak, SpMV ≥ 85 % STREAM, FFT ≥ 40 %, batched dense ≥ 60 % (§14.1; P0 bar = 80 % of these) | **NOT MET** | all-core GEMM 0.39 (M4) / 0.32 (x86) vs bar 0.60; 3-D FFT 0.22–0.25 / 0.17–0.31 vs 0.32; batched dense 10–29 % vs 60 %; SpMV x86 row invalidated by a STREAM-denominator defect | 136 |
| 7 | LBM ≥ 1.0 / 0.6 GLUP/s (§14.1) | UNPROVEN | fs-lbm CONTRACT carries no GLUP/s claim | 7 |
| 8 | Nightly roofline to `metrics`; kernel below band fails CI (§14.4) | PARTIAL | perf-baselines dated 2026-07-11, 90-day policy → expire 2026-10-09; gates opt-in via env var | — |
| 9 | Autotuner persists to `tune` keyed by fingerprint (§5.5) | PARTIAL | GEMM tune rows in fs-session; nothing reads them into TilePlans | 10 |
| 10 | fs-ivl intervals/affine/Taylor/Krawczyk/exact predicates (§6.4) | WORKING | L2 registered | — |
| 11 | fs-la/fs-sparse/fs-fft/fs-rand/fs-cheb/fs-ad/fs-eproc as specified (§6) | WORKING | L2 for sparse assembly; rand_nla goldens; Orr–Sommerfeld; e-processes | — |

### 2.2 Geometry (L2 MORPH)

| # | Goal | Status | Evidence | Beads |
|---|---|---|---|---|
| 12 | Region/Chart/Convert traits + Rep Router as Pareto shortest path (§7.1, §7.3) | PARTIAL | traits, conversion records, sheaf seam (sampled-agreement only) exist; router is scaffolding | 12 |
| 13 | Certified round trips; watertightness faced by oriented-intersection/winding/interval oracles (P1 exit) | PARTIAL | fs-topo certificates L2; continuum successors "explicit planned work" | many (sj31i.43/.44) |
| 14 | Dual contouring with interval bracket verification (§7.3) | UNPROVEN | dual contouring exists; interval bracket proof not found | — |
| 15 | BRIO Delaunay + Ruppert + sliver exudation (§7.5) | PARTIAL | fs-mesh v2/v3 conforming PLC recovery, exact audits, perf ladder; residual near-coplanar slivers noted | — |
| 16 | FrankenVDB sparse SDF; F-rep R-functions with interval/Lipschitz; NURBS exact i128 refinement (§7.2) | WORKING | fs-rep-sdf, fs-rep-frep, fs-rep-nurbs | — |
| 17 | Persistent homology of density fields for ASCENT constraints (§7.8) | NOT_STARTED | no cubical persistence found | 0 |

### 2.3 Physics (L3 FLUX) — what you can actually simulate

| # | Goal | Status | Evidence | Beads |
|---|---|---|---|---|
| 18 | 3-D steady conduction, nonlinear k(T), D/N/R, contact, radiation (EXTREAL E05) | **WORKING** | fs-conduction: MMS ladders, analytic rows, 220 tests green on yto at HEAD; the only production-grade 3-D solver | — |
| 19 | FEEC exact sequences bitwise; matrix-free apply p = 4 ≥ 30 % peak (§8.1) | PARTIAL | derham battery green; vector families are interpolation ladders; no 3-D curl-curl/Darcy solve; **78 % of fs-feec's 19 K lines is EM winding-topology schema for leapfrog I13** | — |
| 20 | CutFEM ghost penalty + cut quadrature on octree/VDB (§8.1, P2) | PARTIAL | **2-D quadtree Q1 only**; 3-D octree NOT_STARTED per CONTRACT | 3 |
| 21 | IGA Kirchhoff–Love shells, NAFEMS suite (§8.1, P4) | **STUB** | fs-iga = 287 lines, 1-D clamped B-spline Poisson | 9 (all EM/leapfrog mentions, no shell bead) |
| 22 | fs-solid: hyperelastic 3-D, TDNNS, rods, fiber beams, buckling (§8.2) | PARTIAL | 2-D nonlinear; 3-D linear K/M only; Cook's membrane, snap-through green | — |
| 23 | LBM D3Q19/D3Q27 cumulant on sparse VDB; cavity Re 1000 vs Ghia, TGV, cylinder Re 100 (§8.3, G2, P3) | PARTIAL | BGK D2Q9/D3Q19 + thermal; cylinder/large-duct G2 lanes are `#[ignore]`; no cumulant, no D3Q27, no free-surface 3-D | 7 |
| 24 | 3-D incompressible Navier–Stokes with scalable saddle preconditioning | NOT_STARTED | fs-flux is 2-D BDM1-P0 on triangles with dense LU | 1 (battery only) |
| 25 | Turbulence (LES/RANS), compressible flow | NOT_STARTED | fs-scenario "RANS solver" is a 282-line 1-D channel behind a feature flag | 0 implementation beads |
| 26 | BEM/FMM bbFMM + Kutta, free wakes; VPM with PSE and remeshing (§8.3, P5) | PARTIAL / STUB | 2-D panel + wake; 3-D non-lifting FMM panels on spheres; VPM 2-D O(N²) direct | 22 FMM mentions (Maxwell-BEM only), 3 VPM |
| 27 | fs-time variational/symplectic/SE(3), all resumable + adjoint (§8.5) | WORKING (library) | frankenscipy ODE oracle casebook | — |
| 28 | fs-couple Dirac ports + IQN-ILS, field-level FSI (§8.4) | **WRONG_SCOPE** | 90 % of the 26 K-line "coupling" crate is musical-instrument synthesis; FSI is a 2-DOF toy | — |
| 29 | Transient multiphysics end to end (conjugate, FSI, thermo-mechanical) | NOT_STARTED | only steady conduction ↔ 1-D air-path fixed point exists | 1 (Euler-specific) |
| 30 | fs-adjoint IFT/Hadamard/Sobolev; FD gate blocks merges (§8.7) | WORKING (library) | 73 % of the crate is DWR-acceptance/explain machinery; ~1.9 K lines of adjoint numerics | — |
| 31 | fs-uq KL/PCE/MLMC/Kanai–Tajimi/CVaR (§8.8) | PARTIAL (toy) | 1.8 K lines; dense KL, fixed-ladder MLMC | — |

### 2.4 Optimization, rendering, orchestration (L4/L5/L6)

| # | Goal | Status | Evidence | Beads |
|---|---|---|---|---|
| 32 | **P2 marquee: topology optimization on a raw SDF, no mesh, composed error certificate** (§16.1) | **STUB** | fs-marquee `study::run_study` = 369 lines behind an off-by-default feature nobody builds; fs-topopt's real fixture is a 3×3×3 Kuhn cube, 12 OC iterations | 20 |
| 33 | Optimizer suite: L-BFGS, TR-NK, AL, CMA-ES BIPOP, NSGA, PDHG (§9.2–9.3) | WORKING (library) | ZDT/DTLZ/WFG, CUTEst-scale, oracle casebooks | — |
| 34 | Any geometry → physics → adjoint → optimizer → new geometry loop on a production PDE | **NOT_STARTED** | only two crates depend on both fs-adjoint and an optimizer; both are toys | 0 |
| 35 | fs-topo SIMP/level-set/homogenization/ground-structure PDHG (§9.5) | PARTIAL | SIMP + Helmholtz + Heaviside on fixed tets; no level-set, no homogenization; PDHG only in fs-truss-e2e | 117 (mostly EM topopt) |
| 36 | fs-sos Lasserre / Burer–Monteiro SOS Lyapunov (§9.8) | STUB | 391 lines; verifies supplied certificates | 17 |
| 37 | fs-surrogate certify-or-escalate over FNO/DeepONet/POD/Koopman (§9.7) | PARTIAL | POD + conformal band + 1-D reduced-basis ladder; no neural operator anywhere | 4 |
| 38 | fs-render spectral PT, sphere tracing ≥ 80/120 Mray/s, differentiable (§10) | WORKING (correctness-first) | real BDPT/MIS/hero-wavelength tracer, Cornell golden 24×24×8 spp; no ray-rate claim; "4K attainment explicitly disclaimed" | — |
| 39 | fs-viz isosurfaces, DVR, LIC, Morse–Smale (§10.3) | PARTIAL | streamlines, marching tets, VTU; DVR/LIC/Morse–Smale "staged" | 23 |
| 40 | FrankenScript: parse + admit + **execute** studies (P10, §11.1) | **STUB as a runtime** | parser/AST/admission/catalog WORKING; zero `execute`/interpreter; 75 % of fs-ir's 53 K lines is leapfrog Machine-IR | 46 |
| 41 | Design Ledger STRICT tables, explain, at(t), forks (§11.2) | WORKING | v20; 4 known-red GC tests on an upstream FrankenSQLite cascade-order bug | — |
| 42 | Sessions: capability tokens, governor, idempotency, `estimate()` (§11.3) | WORKING | consumed by fs-cli `solve` | — |
| 43 | Error Ledger + Time Ledger + planner beating hand budgets (§11.4, P6) | UNPROVEN | fs-plan drives no real work | — |
| 44 | Three-color evidence, laundering refusal, budget pie (AD P3) | WORKING | L2 | — |
| 45 | Falsifier pairing, deny-by-default release gate (AD P6) | WORKING (0A) / PARTIAL (0B) | fs-checker Phase-0A; 0B descriptive only | — |
| 46 | Three flagships: ornithoid, seismic frame, laminar-pour vessel (§15) | PARTIAL (smoke tiers) | honest smoke-tier crates; none at plan resolution | 21 |

### 2.5 Programs added after the core plan

| Program | Plan promise | Status | Evidence | Beads live / closed |
|---|---|---|---|---|
| **Cooling 0.1 (EXTREAL E06)** | one-command run + validate/import/solve/report/package at L3 | PARTIAL | 4 (5) of 6 stages execute; QoI never wired although `fs_airflow::qoi::extract_thermal_qois` exists; report/package/run/compare refuse; the only on-disk report/package implementations are **orphaned files with hard-coded 342.15 K "Verified" output** (wired 08-25 14:05, unwired the same day) | 57 / 150 |
| **Wright Flyer** | real-time browser sim, wake vortices visible, four-flight historical scoring, V-14 perf contract | WORKING product, PARTIAL fidelity | live on Vercel, real 209 KB wasm at 120 Hz, real `fs-wing`/`fs-airscrew`/`fs-atmo`; **fs-vpm not in the browser cone**, physical wake and prop slipstream "declared unsupported"; all 6 vvreg rows self-referential; V-14 NO-DATA; live bundle older than `dist/` | 5 / 164 |
| **Music building blocks** | "hear physics": menus of instrument images, claims registry, listening receipts | WORKING library, UNPROVEN product | 66/72 beads closed; 16 of 45 claim rows green, 29 ungated, 45/45 `live_default=no`; **7 listening receipts, 0 adjudicated**; 1 budget row; no instrument reaches wasm/Apple/CLI; the doctrine's chord-editor consumer does not exist in-repo | 6 / 66 |
| **Euler-disc flagship (t6314 + h7xu5 + b8bxd + jmh21)** | blinded physical prediction; 4K simulation-driven film | PARTIAL / sink | largest crate (107 K lines): 63 % film+audio, 22 % physics, 15 % protocol; the campaign uses the crude coefficient laws, not the 10 K lines of finite-patch/tribology/gas-film adapters; **no retained campaign result anywhere**; film proven at **32×18 px, 2 frames**; 88 open t6314 beads, ~half require a physical laboratory | 145 / 60 |
| **New domains (Geneva, gear, motor, Wankel, genset…)** | five machine flagships in `fs-flagship-e2e` | NOT_STARTED (flagships) | E0 pieces partially landed (six-base dims, fs-matdb, machine IR, weighted operators); **0 of 66 flagship beads closed, 0 of 79 theorem-lane beads closed** | 66+79 / 0 |
| **Leapfrog-2026 (i94v)** | fifteen platform leaps + five theorem ratchets, twenty verification journeys | PLAN | **1,009 live / 38 closed**; 124 open P0s; targets EM, semantic manufacturing, HIL twins, Floquet — on physics that is 2-D or absent | 1,009 / 38 |
| **Apple app** | native SwiftUI studio, 43 kernels, builds for 3 targets, visually inspected | WORKING build, UNPROVEN runtime | 43/43 route to fs-wasm; sim + Catalyst Debug builds succeeded 08-30; **never launched or tested per DerivedData logs**; commit 297ebcf2 "preserve App Store metadata" deleted it; **0 beads** | 0 / 0 |
| **Website frankensim.org** | 30 lab kernels + 10 campaigns + 3 flagships | WORKING, 8 weeks stale | deployed wasm sha == the 2026-07-09 build | 0 |

---

## 3. What specifically IS working right now

1. **The substrate.** Deterministic scalar math, dense/sparse/FFT kernels, intervals and exact predicates, Philox streams, e-processes, the two-lane executor with cancellation, the ledger on FrankenSQLite, evidence colours, packages and the standalone checker. All L2-registered items cite independent oracles and passed on yto at HEAD.
2. **Steady conduction.** A real 3-D P1 FEM with nonlinear conductivity, contact, surface and enclosure radiation, transient (linear) and adjoint, backed by MMS ladders, analytic rows, and 19 Level-A corpus rows. 220 tests green at HEAD.
3. **The product prefix.** `frankensim validate | import | solve` run for real: quarantined STL/STEP import into a ledger, material-card resolution with usage receipts, an interval-bracketed fan/loss-network operating point, and (with a declared section) the conduction solve with min/max temperature, residual and energy closure in the receipt. The governor enforces budgets on this path. The Python client over JSON lines is real.
4. **The Wright Flyer.** A live, replayable, deterministic (one golden across native and wasm lanes) browser sim built on workspace crates, with historical presets, human and pilot control modes, results card, field service, sound, and a 53-file node test suite. The most finished thing in the repository.
5. **Instrument physics.** Viscothermal bores (Ernoult-validated), scattering-port delay lines, port-Hamiltonian ladders, modal plates (Leissa-validated), reeds, brass lip loops, felt hammers, bowed strings, glottis/tract, pickups, jet lab (Brown 1937), selector and PCM render path — all with bake-off receipts.
6. **Rendering.** A deterministic unbiased spectral path tracer with BDPT, AOVs, checkpoint/resume, bit-exact tiling, EXR output, and a cinematic mux pipeline that has run end to end at smoke scale.
7. **Governance that works as built.** check-layers, check-contracts, check-docs, check-maturity, check-schemas, claim-integrity — 7 gates green at HEAD; contracts' no-claim sections are accurate and candid throughout. **The contracts are honest; the README's implied breadth and the registry's one L3 are the problem.**

## 4. What is NOT working or not implemented

### 4.1 Truth defects (claims the repo makes that the code contradicts)

| Defect | Evidence | Severity |
|---|---|---|
| `thermal.conduction-solve` registered **L3** via `scripts/e2e/cooling_01.sh` | the lane drives heatsink-fan, which has no conduction section, so it never runs conduction; it asserts report/package/run exit 5; promotion bead 4gh8t's clause 1 ("runs through to package with ZERO stage gaps") is false; bead 6.11 was reopened for this the same day; the entry stands | Sev-0 (false certificate) |
| QoI bead **s2l9v closed** with all four children | fs-cli still `unreachable!`s at the QoI stage naming s2l9v as the gap; no open bead owns the wiring | Sev-0 (false close on the critical path) |
| `src/report.rs` / `src/package.rs` in fs-cli | tracked, not declared in `lib.rs`, contain hard-coded "Verified 342.15 K" and six literal Estimated claims with `_ledger` ignored; produced the retained `target/cooling-01/report.json` and the repo-root `.fspkg` files | Sev-0 (fabricated evidence on disk) |
| `g0_all_product_workflow_stages_are_integrated` | passes on a USAGE refusal because it calls report/package with no operands | vacuous test |
| Known-red owner bead **f2jag CLOSED** while its three tests are still red | check-suite-receipt 4 violations at HEAD | gate drift |
| README: "no capability at L3" (Current Snapshot, FAQ) vs maturity table "L3 = 1"; ".fsim schema v2" vs `FSIM_VERSION = 3`; Quickstart step 4 says cooling-reference refuses at conduction (it declares conduction and refuses at QoI); steps 5–6 run report/package/run (all exit 5) | README + docs/QUICKSTART.md | doc truth |
| program-metrics.{md,json}, source manifest, SPDX all stale at HEAD | xtask check-program-metrics 2, check-source-manifest 2 | gate drift |
| `frankensim` CLI test `g0_conduction_stage_executes_declared_card_backed_contact` **fails at HEAD** | refuses with `project-conduction-interface-undeclared` (TEMP-DIAG); solve.rs is modified-uncommitted by another lane | the README's "declared conduction executes" is unproven at HEAD |
| Four test targets committed red (08-23 → 08-27), never announced | fs-scenario `rans_card_gates` (feature-gated import), fs-feec `differential_characters` (missing imports), fs-session `snapshot_freeze_gate` (API moved), fs-evidence-runner `value.rs:4158` (undefined `domain`) | workspace test build broken ≥ 9 days |
| Constellation | Mac siblings 25–288 commits ahead of `constellation.lock`; committed Cargo.lock records asupersync 0.4.10 / fsqlite 0.3.13 vs pins 0.4.9 / 0.3.8 | declared TCB ≠ built TCB |
| Apple app metadata | commit "preserve App Store metadata" is 21 deletions; no App Store metadata remains | shipped-surface regression |
| Tracked debug junk | fs-cmaes-viz-wasm: `campaign_err.txt`, `dbg.txt`, `col_test.txt`, `test_out.txt` (~215 KB) | hygiene |

### 4.2 Structural findings

- **Describing outnumbers doing ~3:1 in HELM.** fs-evidence-runner (106.6 K lines, **zero dependents**, grew in 8 commits), fs-vmanifest (43.5 K, authored "G1 drafts" for fifteen initiatives), fs-govern (24.9 K descriptive), Machine-IR assurance/codecs (~40 K) versus ~70–80 K lines that admit or execute user work.
- **The biggest "physics" crates are plumbing for other programs.** fs-feec 78 %, fs-material 72 %, fs-adjoint 73 %, fs-scenario ~60 % non-physics; fs-couple 90 % music. Line counts overstate physics depth by roughly 3×.
- **Effort inversion in the largest crate.** fs-euler-disc-e2e: 68 K lines of film/audio + 16 K of protocol versus 24 K of physics, in a project whose mission is certified simulation.
- **The original plan's own milestones are open.** epic-foundations-huq: 163 live / 11 closed, **156 open P0s**; P0 (bedrock), P1 (geometry + eyes), P2 (marquee), P3 (fluids) milestones all open. The plan's coda — "Build P0. The rest is compounding." — was not followed; the compounding started anyway.

## 5. What is blocking

1. **No single owner-visible journey has a gate.** Nothing forces "a user gets an answer" before scope expands. The one L3 was promoted on a lane that never ran the capability; nothing caught it.
2. **Vision fan-out outran delivery.** Two beadification waves (W29: 1,522; W31: 693) created more beads than have ever been closed, most for leapfrog/new-domain/theorem work sitting on 2-D or absent physics. 464 live beads are tagged moonshot; 283 are theorem beads; 66 flagship + 79 theorem-lane beads have zero closures.
3. **The graph cannot dispense work.** 83 % blocked, 2 ready, 303 open epics, 368 open P0. Priority has no signal. 152 of 204 in_progress beads are dead sessions.
4. **Truth decay.** False closes and an unbacked L3 mean "closed" and "L3" are no longer trustworthy inputs to planning.
5. **Build hygiene under swarm load.** Tests committed red for nine days; the Mac constellation drifted; rch root-workspace offload is refused (hz3 alias) and clients die with exit 144; test suites take 26 minutes for one crate on the quiet box.
6. **Physical reality.** Half the Euler flagship and all of L4 validation need a laboratory, instruments, and licensed corpora that software agents cannot produce.
7. **Performance bars.** Every measured P0 kernel misses; baselines expire 2026-10-09; the fleet has no citable x86 GEMM row.

## 6. If we implemented every open and in-progress bead, would the gap close?

**No.** Four reasons, each independently sufficient:

1. **The open beads are mostly about the wrong layer.** 1,009 leapfrog + 66 machine-flagship + 79 theorem-lane + 99 Euler-lab beads (≈ 1,250 of 1,900 live) build schema, compilers, assurance overlays and theorem cards on top of physics that is 2-D, toy-scale, or absent. Completing them yields more descriptions of studies that still cannot run.
2. **The core plan's physics gaps have almost no implementation beads.** 3-D CutFEM octree (3 mentions, no implementation bead), 3-D incompressible NS (1 battery bead, no solver bead), IGA shells (0), FMM-accelerated VPM (0), turbulence (0), compressible (0), field-level FSI (0 generic), 3-D lifting BEM (0). The P2 marquee has 20 mentions and no bead that un-gates and ships it.
3. **The product critical path has a false close.** "All open beads" excludes s2l9v (QoI) because it is closed; the product would still stop at the QoI gap.
4. **~90 beads cannot be closed by software.** Specimens, vacuum chambers, launch rigs, DAQ, custody chains, blind unlocks.

## 7. Vision goals with NO bead coverage (or none that implements them)

- IGA Kirchhoff–Love shells / NAFEMS shell suite (P4 exit criterion)
- 3-D CutFEM on octree/VDB (P2 marquee prerequisite at 3-D)
- 3-D incompressible Navier–Stokes solver + scalable saddle preconditioner
- FMM-accelerated VPM with PSE viscous diffusion and remeshing (P5)
- Turbulence modelling on a mesh (LES/RANS); compressible flow
- Field-level FSI / conjugate CFD (any transient multiphysics)
- 3-D lifting-surface BEM with wake sheets
- Persistent-homology density constraints
- FrankenScript executor (46 mentions, none builds an interpreter)
- Direct volume rendering, LIC, Morse–Smale (fs-viz staged items)
- Neural-operator / DEIM / Koopman surrogates
- Sphere-traced ray-rate target; LBM GLUP/s target; ≤ 200 µs cancel-latency claim
- QoI stage wiring in fs-cli (the bead exists but is closed)
- Apple app (any bead); website redeploy; WF live-bundle refresh; WF resolved wake in the browser cone
- Music: owner listening adjudication; any product exposure of an instrument
- Cooling: a physical example that declares conduction; an L3 lane that actually runs conduction

---

## 8. Bridge plan

### Reading of the gap

The distance between vision and code is not a list of missing features; it is a **missing spine**. The workspace has more than enough substrate to run one complete design study, and it has never done so. Every hour spent on layer-6 descriptions of studies that cannot execute, or on flagships whose physics is not built, widens the gap. The bridge therefore has five moves, in strict order: make the repo tell the truth; ship one real journey to a legitimate L3; ship the plan's own marquee (the continuum thesis, demonstrated once); park everything that is not on those two paths behind an explicit, owner-signed decision; then reopen the frontier with kill criteria. Beads for this plan are few and mostly product code. Where a process item appears it names the exact defect it prevents.

### Move 0 — Truth (days, not weeks)

- T1. Fix the four red test targets; make `cargo test --workspace --no-run` part of the DSR quality gate so a red test target can never again sit unannounced for nine days.
- T2. Demote `thermal.conduction-solve` to L2 in `capability-maturity.json` **or** make `cooling_01.sh` run a conduction-declaring project through every stage. Until it does, L3 = 0. Add one rule to check-maturity: an L3 entry must cite a retained receipt whose executed-stage list includes the capability (the receipt schema already exists in spine-e2e-summary).
- T3. Reopen s2l9v (QoI) with the code condition as its DONE-WHEN; reopen f2jag or file the successor that owns the three still-red ledger tests.
- T4. Replace the orphaned `report.rs`/`package.rs` bodies with ledger-traced renderers (Move 1 does the work); until then they must be declared as gaps, not left as fabricated evidence. Remove the vacuous `g0_all_product_workflow_stages_are_integrated` assertion or make it call the verbs with operands.
- T5. README/Quickstart truth pass: L3 count, `.fsim` v3, step 4/5/6 expectations, "9 standalone workspaces" (three have no FrankenSim code), Apple "43 kernels" is fine; add an explicit "what stops where" line.
- T6. Regenerate program-metrics and the source manifest at HEAD; restore Apple Info.plist metadata; untrack the four debug files in fs-cmaes-viz-wasm; register or remove the fs-simd wasm capsule; move fs-g1-train's serde behind a dev/feature boundary.
- T7. Constellation decision: advance `constellation.lock` to the fast-forwarded sibling heads under the compatibility-train process (all seven are strict fast-forwards) so the declared TCB equals the built TCB.

### Move 1 — One real journey: Cooling 0.1 to a legitimate L3 (weeks 1–4)

Target: `frankensim run examples/heatsink-fan/heatsink-fan.fsim ledger.db --materials aa6061.fsmcdpk` completes six stages and the user reads a QoI verdict with an evidence colour, a report rendered from the ledger, and a package the standalone checker accepts.

- J1. Give heatsink-fan (and cooling-enclosure) a real `cooling.conduction` section: solid regions, contact traces, convection boundaries bound to the flow-network operating point.
- J2. Conjugate coupling (s93ej.3): Robin coefficients from `fs_airflow` branch velocities → `fs_convection` cards → conduction boundary rows; partitioned fixed point with Anderson/IQN-ILS relaxation (fs-couple already has IQN-ILS), energy balance gate.
- J3. QoI stage: wire `fs_airflow::qoi::extract_thermal_qois` into `solve.rs`; delete the `unreachable!`; requirement verdict with composed colour; **first Verified-colour QoI bound** via the existing heat adjoint + a DWR-style residual estimate on the P1 solve (fs-adjoint heat adjoint + fs-dwr pattern lifted from 2-D CutFEM to P1 tets).
- J4. Report: deterministic HTML + JSON twin rendered from ledger artifacts and receipts (no literals); package: format-9 evidence package from retained receipts with the source-manifest identity; `run` completes; `compare` diffs two runs' QoIs.
- J5. Make `cooling_01.sh` the honest lane: run J1's project through all six stages, retain the receipt with stage list and digests, fail on any gap; regenerate spine-ratchet/spine-e2e-summary; then promote to L3 through check-maturity with the receipt rule from T2.
- J6. Fix the failing `g0_conduction_stage_executes_declared_card_backed_contact` (coincident boundary-slot orientation) as part of J1.

### Move 2 — The marquee: the continuum demonstrated once (weeks 3–10)

Target: the plan's P2 exit criterion, un-gated and shipped: topology optimization on a raw SDF with CutFEM, no mesh in the loop, composed error certificate, as a tracked example/CLI verb with a retained receipt.

- M1. Un-gate fs-marquee: default feature, built and tested by the workspace suite; `frankensim study marquee.fsim` (or a `--study` verb) runs SDF → 2-D CutFEM elasticity → self-adjoint gradient → SIMP/level-set update → DWR certificate → ledger, with a retained golden and a receipt. This is the first geometry → physics → adjoint → optimizer → geometry loop in the repo.
- M2. 3-D cut quadrature on tets/hexes from SDF (Saye-style high-order quadrature on implicitly defined regions, or moment fitting) in fs-cutfem; 3-D ghost penalty; assembly through fs-sparse with p-MG/AMG from fs-solver.
- M3. 3-D marquee on an fs-rep-sdf VDB grid: cantilever/bracket compliance minimization with volume constraint, DWR-driven octree refinement, certificate composed from cut-quadrature error + solver residual + adjoint FD gate. Registry entry `optimization.marquee-topopt` at L2 then L3.

### Move 3 — Portfolio decision (owner-signed, week 0)

- D1. Park (status `deferred`, labels kept) until Moves 1–2 close: leapfrog-2026 (1,009 beads), ext-theorem lanes (79), ext machine flagships (66), Euler t6314 physical-lab beads (~50), b8bxd/jmh21 enablers not consumed by WF/marquee. Parking is reversible and preserves intent; it removes 1,200+ items from the swarm's field of view.
- D2. Finish-out lanes with explicit done lines: Wright Flyer (refresh live bundle; measure V-14 on one qualified device; either put fs-vpm's wake in the browser cone or delete the "watch the vortex sheets" promise from the plan); Music (owner adjudicates the 7 listening receipts; expose one instrument in fs-wasm and the Apple catalog; retain the 3-D jet receipts in-tree); Euler cinematic (one 1080p daily frame or park); Apple (restore metadata, run XCTests, launch on 3 form factors with screenshots); website (redeploy current fs-wasm).
- D3. P0 performance decision: measure SpMV x86 with the corrected STREAM denominator and record the cancel-latency claim; then either fund GEMM/FFT/batched kernel work to the 80 % bars or re-baseline the plan's targets to measured values with a dated receipt. Either way, the P0 milestone bead closes with a true statement.
- D4. Bead-graph hygiene: audit the 152 stale in_progress beads (finished-but-unclosed vs abandoned), release the abandoned ones; re-tier P0 to the Moves 0–2 critical path only (tropical critical path over the graph already exists — use it); collapse the 303 open epics to the ones with a live child.

### Move 4 — Reopen the frontier with kill criteria (after Move 2)

- F1. One epic per core-plan physics gap, each stating the plan section, the current state ("not started"), the gate it waits on, and a retire criterion: 3-D incompressible NS + saddle preconditioning; IGA shells or formal retirement of the promise; FMM-VPM or retirement; turbulence retirement from v1; field-level FSI via fs-couple ports (and move the instrument code out of fs-couple into its own crate so the coupling crate is a coupling crate).
- F2. FrankenScript executor v0: lower an admitted study to the fs-cli stage pipeline (Cooling and marquee) so P10 stops being parse-only; or amend the README to "IR and admission only".
- F3. Un-park leapfrog/new-domain work initiative by initiative, only where the underlying solver exists at 3-D fidelity with a G1 ladder.
