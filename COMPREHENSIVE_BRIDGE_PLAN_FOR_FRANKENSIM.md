# Comprehensive Bridge Plan for FrankenSim

**Status:** round 3, 2026-09-02 (early). Revised in place; rounds are logged in §18.
**Inputs:** `docs/REALITY_CHECK_2026-09-01.md` (Phase 1: where the code really is), the steering epic `frankensim-rc-root-q61wp` (.1–.43, label `reality-check-2026-09`), the five owner decisions of 2026-09-01, and the landings of 2026-09-01/02: the report stage and `run`; the finned heatsink solving seven stages; the fs-mesh facet-recovery fix; the A8 gap-table truth fix.
**Purpose:** the Phase 2 document of the reality check — a plan that closes **every** gap between `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md` (plus its addendum) and the code, granular enough that each item becomes a self-contained bead and no reader needs this document afterwards. The plan is the measuring stick; the code is the ground truth; every item names the code it starts from.

---

## 0. How to read this plan

- **Item IDs** are stable (`A1`, `B4`, …). A bead that implements an item names it. Existing steering beads are cited as `q61wp.N`.
- Every item carries: **goal** (plan section served), **state** (what the code does today, with the file or test that proves it), **work** (what changes), **seams** (what already exists and is reused rather than rebuilt — verified names, not guesses), **falsifier** (the test that fails if the item were faked), **receipt** (what is retained), **logging** (what a failing run must print so it is diagnosable without a debugger), **kill** (when to stop and say so), **depends**, **size** (S = a session, M = days, L = weeks).
- **No claim advances on prose.** An L3, a README fact, a critical-path bead close each resolve to a retained, digest-pinned receipt from an executed lane (§2). That rule is the plan's spine; everything else is work.
- The plan **does not proliferate process**: no new dashboards, analyzers or harnesses beyond the receipt rule and the gates that already exist (`xtask check-*`, the two shell lanes, the G0/G1 CLI battery). Where a lane already exists it is repaired, not duplicated.

## 1. The gap, in one page (state at 2026-09-02 early)

**What the vision promises.** A certified simulation product: a user declares a design study and the product either answers with an evidence colour and an error budget or refuses by name; the continuum geometry → physics → adjoint → optimizer → geometry is demonstrated once, for real; the fluid frontier (NS, LBM, VPM, compressible, turbulence) and the structural frontier (IGA shells) exist at 3-D with G1 ladders; performance meets the §14.1 bars; every claim traces to a receipt.

**What the code delivers, measured.**

| Journey / promise | 2026-09-01 morning | 2026-09-02 early |
|---|---|---|
| `frankensim run` on a real body | stopped at the QoI gap; report/package were fabricated literals | **all seven stages execute on the finned heatsink** (`examples/heatsink-fan`, 108-facet single shell) and on the cooling reference; report HTML + JSON twin + format-9 package sealed in the ledger; `report`/`package` export the retained bytes; `scripts/ci/solve_stage_producers_e2e.sh --profile full` 65/65, `scripts/ci/examples_freshness_e2e.sh` 22/22 |
| Conduction on real geometry | proven only on 4-facet tetrahedra; the heatsink STL was five glued boxes; every comb refused at facet recovery at any budget | fs-mesh facet recovery accepts any coplanar tiling and iterates to a fixed point (`crates/fs-mesh/tests/comb_prism.rs`: one/two/four fins, exact volume, default budget); conduction mesh on the heatsink = 246 vertices in 4.6 s |
| Conjugate solid/air coupling | declared coefficients only | derived Robin rows from the flow-network operating point through `fs-convection` cards with domain gating (Hausen in domain at Re ≈ 1.5 × 10³; energy imbalance 9 × 10⁻¹¹ W; 12 fixed-point iterations) |
| QoI with the 8-term budget | unreachable | executes; all eight terms are honest NO-DATA → Estimated / indeterminate |
| Gap table truth | conduction and QoI were "typed gaps" even when they executed | `gap_dependency` is `None` for every stage; undeclared inputs refuse by name (`cli-solve-conduction-undeclared`, `cli-solve-qoi-undeclared`, exit 4) |
| L3 capabilities | one, falsely certified | **one, legitimately**: `thermal.conduction-solve` at L3 on the retained seven-stage lane receipt (65/65 at 503205c8; commit 974b1cc6), stated in the registry and README as "executes end to end with retained receipts, not converged or validated"; five gates at 0 violations on yto |
| Marquee (P2) | gated, library did not compile | compiles; `mq_004` shows the gradient sign is wrong (q61wp.16 owner) |
| Frontier (NS, 3-D LBM, VPM, compressible, turbulence, IGA) | absent or 2-D | unchanged; IGA / turbulence / compressible / FMM-VPM **retired from v1** by owner decision (q61wp.40–.43); NS gated behind the LBM wedge (q61wp.35) |
| Graph | 2 ready, 83 % blocked, 204 in_progress | 41 ready, 164 in_progress, 1,173 beads deferred to 2026-12-01 |

**Why the gap persisted.** Claims were self-certified prose, and nothing gated scope on a delivered journey. Both mechanisms are being removed: R1 below, and the ordering of this plan (journeys before frontier).

**What remains, by weight.** (i) Journey A *runs* but is not *verified*: the conduction mesh is the recovered surface with no refinement, the budget terms are NO-DATA, there is no independent consumer, no physical anchor, no promotion receipt. (ii) Journey B has a compiling library and a falsified gradient. (iii) Truth repairs are mostly landed and must be closed with evidence, and two new ones were found tonight (A9 fabricated attribution, A12 quality slivers). (iv) Shipped surfaces need done-lines. (v) Performance and constellation decisions need receipts. (vi) The frontier waits for the LBM wedge.

## 2. The one rule and the standing doctrines

**R1 — Claims read receipts.** A capability level, a README fact, or a critical-path close resolves to a retained receipt from an executed lane that names the stages it ran, with digests, on a stated commit. Parts that exist: `spine-e2e-summary.json` (schema `frankensim-spine-e2e-receipt-v1`), `spine-ratchet.json`, `xtask check-maturity` (`receipt_stage_run` requires a top-level `stages` array of `{capability, stage, status:"executed"}` rows), `check-docs` + `doc-facts-inventory.json`, `check-spine-metrics`, `check-suite-receipt`. Delta: the lane emits the `stages` array (`--retain-receipt PATH`, q61wp.13), the DSR quality gate includes a no-run build of every workspace test target (q61wp.1).

**R2 — Falsifier before feature.** Each item names the test that fails if the item is faked: a hostile twin, a permutation injection, a domain-boundary probe, a cross-ISA replay. A feature bead without a falsifier does not close. Tonight's precedent: the tiling criterion looked complete by area and by row count, and the independent winding audit caught the leaked flood anyway — the audit *is* the falsifier.

**R3 — Measured tolerances carry their measurement.** A gate set above a measured floor states the measurement and keeps ~5× headroom (the comb test pins ≤ 600 facet Steiner points against a measured 139); when a fixture is re-dimensioned the measurement is redone in the same commit.

**R4 — Refuse by name, never approximate.** Outside a card's domain, without a declared input, on an open shell: refuse with a code, a message and a fix. Tonight's examples: Gnielinski refusing L/Dh = 5 < 10; the volumetricizer refusing the five-box soup; `report` refusing an unknown run without writing; the undeclared-conduction refusal that replaced a false "waiting for a producer" message.

**R5 — Delete nothing; rewrite in place; park reversibly.** `br defer` with a written un-park criterion is the only portfolio instrument.

**R6 — No process porn.** No new harness, dashboard or registry unless it replaces a Sev-0 mechanism observed in the audit. Existing lanes are repaired.

**R7 — Permutation and order are inputs.** Determinism is bit-for-bit for identical bytes. A permutation of the input (region order, vertex order) is a *different* input and may yield a different mesh; any statement meant to be permutation-invariant is made on geometry or identity, never on indices (fs-mesh CONTRACT §14, tonight's g3 test rework).

### 2.1 Mathematics that buys alpha on the two journeys (each entry names its item and its falsifier)

The substrate already carries exact predicates, intervals, e-processes and deterministic streams; the items below use mathematics from the last sixty years to turn "runs" into "verified" with fewer tets, fewer runs and stronger claims than the textbook path.

1. **Goal-oriented (dual-weighted-residual) mesh adaptation for the QoI — B2.** Becker–Rannacher (2001): refine where the adjoint of the junction-temperature functional weights the residual, not uniformly. fs-conduction has the heat adjoint; fs-dwr exists. Alpha: a converged junction temperature at a fraction of the tets a uniform ladder needs, and an error *indicator* per element that the report can show. Falsifier: on the heatsink, DWR-adapted meshes reach the uniform ladder's fine-rung QoI within the stated tolerance with ≤ ¼ of its tets; the effectivity index is reported and stays in [0.5, 2].
2. **Equilibrated-flux, constant-free error bound — B6.** Prager–Synge hypercircle through Braess–Schöberl (2008) and Ern–Vohralík (2015): reconstruct an H(div) flux on the RT0 patches (fs-feec `whitney`), and the energy-norm error is bounded above *without an unknown constant*; the QoI bound follows by the adjoint. This is what lets the discretization term carry the Verified colour. Falsifier: the bound brackets the true error on every rung of the B2 ladder; a mesh coarsened on purpose widens it.
3. **Conforming-Delaunay termination theory for recovery — B4.** Murphy–Mount–Gable (2001) and Cohen-Steiner–de Verdière–Yvinec (2002): midpoint refinement of segments and facets terminates when the local feature size is bounded below; tonight's tiling acceptance and fixed-point passes are the practical form. Alpha: a stated, checkable reason recovery ends, and a Steiner-count bound proportional to facet area over the square of the local feature size that becomes the budget formula in B3. Falsifier: the corpus bodies' Steiner counts stay under the formula's bound.
4. **Simulation of simplicity for order-canonical meshes — B12.** Edelsbrunner–Mücke (1990): resolve degeneracies by a symbolic perturbation keyed on a canonical *geometric* order (lexicographic coordinates) instead of insertion index, so identical point sets in any order tie-break identically. Alpha: R7 collapses to "identical geometry, identical mesh"; the g3 slot-swap invariant returns. Kill: if the exact audit or the kernel's cavity logic cannot take the keyed perturbation in a week, keep R7.
5. **Sliver exudation and radius-edge refinement — B2.** Cheng–Dey–Edelsbrunner–Facello–Teng (2000) and Shewchuk's Delaunay refinement with protected constraints: `fs_mesh::exude` and `refine` already implement the mechanics; the plan makes them mandatory after carve with a dihedral floor. Falsifier: a forced 1° sliver refuses.
6. **Anytime-valid (e-process) acceptance for replay and falsifier lanes — B7.** Ramdas–Grünwald–Vovk–Shafer (2020–23): the determinism and hostile-twin lanes accumulate evidence across runs and hosts with an e-value, so "the twin flips the verdict" and "two ISAs agree" are sequential tests that never need a fixed sample size and never inflate error under optional stopping. The substrate already ships e-processes. Falsifier: a deliberately injected 1-in-20 nondeterministic run drives the e-value below the threshold within the stated number of replays.
7. **Certified cut quadrature by interval root isolation — D2.** Saye (2015) dimension reduction with roots enclosed by `fs_ivl::newton_roots_bounded` (Krawczyk / Hansen–Sengupta interval Newton): the geometric quadrature error is an enclosure, not an estimate, so the cut-cell term of the marquee certificate is rigorous. Falsifier: sphere-volume enclosure tightens with level; a wrong root box is caught.
8. **Ghost-penalty stabilization — D2/D3.** Burman (2010), Burman–Hansbo (2012): condition-number control on cut cells independent of how the interface cuts them; without it the 3-D marquee's solver counts explode. Falsifier: the condition number stays bounded across a sweep of interface positions on one cell.
9. **Topological derivative with Sobolev-smoothed descent — D1.** Sokołowski–Żochowski (1999), Allaire–Jouve–Toader (2004) level-set method; `fs_adjoint::sobolev_smooth` exists. Falsifier: the adjoint-FD gate per iteration; monotone compliance decrease on the cantilever.
10. **Tropical (max-plus) critical path over the bead graph — E2.** `xtask generate-tropical-path` exists: the P0 re-tier uses zero-slack paths to Journeys A/B rather than centrality. Falsifier: after re-tiering, every P0 lies on a zero-slack path or carries a written reason.
11. **Quasi-Newton coupling acceleration if the fixed point stalls — B5.** Degroote–Bathe–Vierendeels (2009) IQN-ILS; today Aitken (`fs_couple::iterate_aitken`) suffices at 12 iterations. Kill: implement only if a corpus body exceeds the iteration cap.

## 3. Definition of done for v1 (the acceptance test of the whole plan)

An outsider with the repository, the constellation at pins, and one Mac or Linux box can, in order:

1. Build the workspace with zero red test targets (`check-suite-receipt` 0 violations) — A1.
2. Run `frankensim run examples/heatsink-fan/heatsink-fan.fsim ledger.db --materials aa6061.fsmcdpk` and read a report whose junction-temperature QoI carries a colour, an 8-term budget with the discretization and coupling terms **measured**, a mesh-convergence statement over three rungs, and a Verified-colour discretization bound — B2, B5, B6.
3. Re-derive that verdict from the package with the solver-free checker (`fs_checker::check_against_root`) and from the Python client (`python/frankensim/client.py`), and get the same colour — B7.
4. Swap the aluminium card for the foam card and watch the verdict flip; omit `--materials` and watch material-resolve refuse; strip `cooling.conduction` and watch conduction refuse by name — B7, A8.
5. Find `thermal.conduction-solve` at L3 in the maturity registry citing that lane's receipt (B10) and nothing else at L3 without a receipt (A2).
6. Run `frankensim study` on the tracked marquee project and get a topology-optimized design with a composed certificate and a Goodhart rung-climb delta — D1–D3.
7. Read a README whose every quantitative statement is generated from or checked against `doc-facts-inventory.json`, whose quickstart runs the heatsink end to end, and whose capability table matches the registry — A5.
8. Open `br ready` and find ≥ 10 items, every one on a journey path; find no in_progress bead older than 7 days without a live claim — E2.
9. Find §14.1 either met or re-baselined with dated receipts on the same page — F1, F2.
10. Find every shipped surface (Wright Flyer, instruments, Apple app, website, Euler cinematic) with a done-line and a receipt, or a dated park — H1–H5.

### 3.1 Coverage matrix: every uncovered vision goal from the reality check (§7) → an item or a decision

| Vision goal with no implementing bead (2026-09-01) | Disposition in this plan |
|---|---|
| IGA Kirchhoff–Love shells / NAFEMS suite | retired from v1 (q61wp.40); I5 |
| 3-D CutFEM on octree/VDB | D2, D3 |
| 3-D incompressible NS + saddle preconditioner | gated epic (q61wp.35); I2 after I1 |
| FMM-accelerated VPM | retired (q61wp.43); I5 |
| Turbulence on a mesh; compressible flow | retired (q61wp.41, .42); I5 |
| Field-level FSI / conjugate CFD | I1 (LBM conjugate channel) then I2; B5 is the reduced form that ships now |
| 3-D lifting-surface BEM with wake sheets | H1 decision (Wright Flyer wake promise) |
| Persistent-homology density constraints | D3 stretch (§7 D4), not v1 |
| FrankenScript executor | I4 |
| Direct volume rendering, LIC, Morse–Smale | not v1; README states "staged" (A5) |
| Neural-operator / DEIM / Koopman surrogates | not v1; parked with un-park criterion (E1) |
| Sphere-traced ray-rate, LBM GLUP/s, ≤ 200 µs cancel latency | F1 measurements; F2 restatement |
| QoI stage wiring | landed (A3 closes s2l9v on the receipt) |
| Apple app, website redeploy, WF bundle, WF wake | H3, H4, H1 |
| Music: listening adjudication; product exposure | H2 |
| Cooling: physical example with conduction; an L3 lane that runs conduction | landed (B1); B10 in flight |

## 4. Workstream A — Truth and gates

**A1. Four committed-red test targets + DSR no-run build.** *Goal:* R1. *State:* q61wp.1 in progress; the four targets named in the audit were red 08-23…08-27 and several have since been fixed by other lanes — re-verify at HEAD with `cargo test --workspace --no-run` on the rch lane (root-workspace builds are refused locally; see G2). *Work:* every target builds; the DSR quality gate gains the no-run build; `docs/CI_GATES.md` records the gate line and its expected wall-clock. *Falsifier:* revert one fix locally; the gate refuses. *Receipt:* suite receipt at HEAD, 0 violations. *Logging:* the failing target name and the first compiler error line. *Size:* S–M.

**A2. Receipt rule in check-maturity.** *Goal:* R1. *State:* `xtask/src/maturity.rs::receipt_stage_run` exists and demands a `stages` array; the lane does not yet emit it; L3 count is 0. *Work:* (with q61wp.13) the lane emits `stages` rows from the per-stage status lines it already checks; check-maturity verifies receipt existence, the stage row `executed`, and receipt digest newer than the lane script; refusal text names the missing field. *Falsifier:* an L3 entry citing the old cooling_01 receipt (never ran conduction) fails; an entry citing tonight's receipt without `stages` fails. *Size:* S. *Owner:* NobleLion (in flight).

**A3. Reopened beads closed on code conditions.** *State:* s2l9v reopened; QoI executes; f2jag ready with a live test. *Work:* close s2l9v citing the 65-check lane (stage `qoi` executed) and the driver-version note; leave f2jag to its owner until its three FrankenSQLite tests are green. *Size:* S.

**A4. Fabricated renderers replaced — close with evidence.** *State:* `report.rs`/`package.rs` are ledger-traced exports; the workflow test calls the verbs with operands; `run` chains seven stages; g1 CLI tests pin it. *Work:* close q61wp.3 and q61wp.12 citing the lane and the tests. The stray `cooling_demo_run.fspkg`, `profile_run.fspkg`, `feedface….fspkg` at the repo root are droppings of the old stub verbs: **owner permission required to delete** (R5). *Size:* S.

**A5. README / QUICKSTART / plan §16 truth pass.** *State:* q61wp.4 open. Tonight fixed the command-reference paragraph and three example READMEs; the quickstart steps, the capability table narrative, the schema version (v4), the workspace inventory line and plan §16 remain. *Work:* rewrite "what happens where" around the seven stages and the two named refusals; quickstart runs the heatsink end to end and shows the export paths and the unknown-run refusal; `.fsim` v4; L3 count from the registry; "9 standalone workspaces" → the measured split; plan §16 status column (P0–P3 open, dated); regenerate `doc-facts-inventory.json` and pass `check-docs`. *Falsifier:* `check-docs` 0; a hand run of the quickstart on a clean clone by someone who did not write it. *Depends:* B10 for the L3 line (write "as of <commit>"). *Size:* M.

**A6. Derived registries at HEAD.** *State:* program-metrics and the source manifest regenerate cleanly with the native xtask (the tracked `target/debug/xtask` is a Linux ELF; G4); spine files are NobleLion's for q61wp.13. *Work:* regenerate source manifest + SPDX + program-metrics as the last step of each commit that adds tracked files (they are pure functions of HEAD; stage first because the generators read the index); fix the Info.plist metadata; untrack the fs-cmaes-viz-wasm debug files (owner permission); gate fs-g1-train's serde behind a feature. *Falsifier:* `check-source-manifest`, `check-program-metrics`, `check-spine-metrics` all 0 at HEAD. *Size:* S–M.

**A7. Contact test at HEAD.** *State:* q61wp.6 ready; the g3 undeclared-interface pair is green again (message format fixed tonight); `g0_conduction_stage_executes_declared_card_backed_contact` status must be re-read from tonight's full suite run. *Work:* if red, fix the coincident boundary-slot orientation or the fixture's interface declaration; the test stays the only executable proof of declared contact. *Size:* S–M.

**A8. Gap-table truth — landed 2026-09-02.** `gap_dependency` is `None` for all seven stages; `declaration_gap` refuses undeclared conduction (`cli-solve-conduction-undeclared`) and undeclared QoI (`cli-solve-qoi-undeclared`) at exit 4; fs-cli CONTRACT and the spine ratchet read 7 executing stages. *Close* with the suites and the ratchet regeneration.

**A9. Fabricated attribution in the vertical-profile lane.** *State:* `scripts/ci/e2e_extreal_vertical_profile.sh` synthesizes a "PipelineAttributionReceipt" whose top-three kernels are literal fractions of the solve time (`t_solve * 0.45`, `0.35`, …) with invented arithmetic intensities and roofline regimes. That is a fabricated receipt of exactly the kind the audit found in `report.rs`. *Work:* keep the per-stage wall-clock rows (they can be read from the retained stage receipts' `wall_s`); delete the kernel rows or produce them from a real profiler run with the tool named; the receipt schema gains `measured: true|false` per row and any promotion-context consumer refuses `false`. *Falsifier:* a run with the solve stage stubbed to zero time must not produce kernel rows. *Size:* S.

**A10. Freshness lanes have owners.** *State:* `examples_freshness_e2e.sh` (22 checks) and the vertical-profile script are run by no xtask gate; the G0/G1 CLI battery is the executable proof. *Work:* register the freshness lane in `docs/CI_GATES.md` as a native-binary lane with its expected check count; state that the CLI battery is the gate. No new harness. *Size:* S.

**A11. Marquee gradient sign.** *State:* `fs-marquee` `mq_004` falsifies the gradient sign. *Work:* fix against a finite-difference gate per iteration (D1 prerequisite); the test stays. *Owner:* q61wp.16. *Size:* S.

**A12. Mesh-quality truth in the conduction receipt.** *State:* the conduction receipt records vertices and iterations, not quality; recovery can leave slivers; tonight's ulp-twin adoption removed one class of them. *Work:* the receipt records tets, min dihedral, max radius-edge, recovery statistics; the README's "conduction executes" line states the mesh is unrefined until B2 lands. *Falsifier:* a receipt without the quality fields fails `check-schemas` for the conduction receipt schema. *Size:* S. *Depends:* B3 (same receipt edit).

## 5. Workstream B — Journey A: Cooling 0.1 from "runs" to a legitimate L3

**Target.** `frankensim run examples/heatsink-fan/heatsink-fan.fsim` yields a junction-temperature QoI with a colour, a measured 8-term budget, a report and a checker-accepted package — and an outsider can falsify it in the four ways of §3.

**B1. Heatsink declares conduction — landed.** *State:* v4 project; single-shell 108-facet STL from `examples/heatsink-fan/generate_heatsink_stl.py`; corrected vent arithmetic (0.00072 m²); Hausen developing-flow card in domain at 70 % fan speed; seven stages; freshness lane 22/22. *Remaining:* the cooling-enclosure example gets a real body (its `plate.stl` is a tetrahedron) and the derived law; heated-plate stays the schema tour. *Falsifier:* the freshness lane. *Size:* S–M.

**B2. Mesh size control and quality after recovery (load-bearing).** *Goal:* the QoI is a converged discretization, not a coarse estimate; the refinement is goal-oriented (§2.1 item 1) with a uniform h-ladder as the control. *State:* the conduction mesh is the recovered PLC itself (246 vertices for the heatsink); `fs_mesh::refine` (`RefineOptions { max_radius_edge, max_steiner, split_hull_facets, min_edge_factor }`, worst-first circumcenter insertion, hull-escaping offenders skipped and counted) and `fs_mesh::exude` (sliver perturbation) exist and are not called by the conduction stage. *Work:* after carve-and-label: constrained Delaunay refinement to a target edge length derived from the declared accuracy budget and the geometry's smallest declared feature (default feature/4), with recovered walls protected by the encroachment rule the contract already sketches, then exudation; then a DWR loop (solve → adjoint of the QoI → element indicators → local refinement of the top fraction, walls protected) until the indicator sum is below the accuracy budget or the memory budget stops it; a hard quality floor (min dihedral ≥ 5°, max radius-edge ≤ 2) below which the stage refuses rather than solves; mesh statistics in the receipt (A12). *Seams:* `refine`, `exude`, `VolumetricPolicy`, the audited complex, `fs_conduction::ConductionMesh::new_region_owned`. *Falsifier:* a three-rung h-ladder on the heatsink whose QoI converges monotonically at the P1 rate against a fine reference; a fixture with a forced 1° sliver refuses. *Receipt:* per-rung vertices/tets/min-dihedral/QoI in the lane's JSONL. *Logging:* rung, target h, achieved h, Steiner count, worst dihedral, QoI. *Kill:* if refinement cannot reach the quality floor on the heatsink within the memory budget, the stage refuses and the plan says so — no silent coarse solve. *Depends:* B1. *Size:* M–L.

**B3. Recovery budget from the declared memory budget; receipt fields — landed 2026-09-02 (q61wp.45).** *State:* `recovery_budget(memory_bytes)` derives `max_steiner` with the fixture default as the floor (64 MiB → 32 768); the conduction receipt is schema v4 with `mesh.quality` (tets, vertices, min dihedral, max radius-edge, slivers below 5°, flat tets — `fs_mesh::QualityCensus`, det-routed) and a `recovery` object (memory budget, depth, Steiner cap, segment and facet statistics — `fs_mesh::RecoveryEvidence`); disclosure, not enforcement. **Measured on the four-fin comb: min dihedral 0.000°, 9 flat tets of 722, 301 tets above radius-edge 2** — the input to B2. *Previously:* `solve.rs` derived `max_tets` from `budgets.memory_bytes` and hard-coded `RecoveryOptions::default()`. *Work:* derive `max_steiner` from the same budget with the default as a **floor** (never below today's behaviour); write `max_steiner`, `max_depth`, the memory budget, `FacetRecoveryStats` and segment stats into the conduction receipt so the mesh row is reproducible from the receipt alone. *Falsifier:* a project with a tiny memory budget still recovers the heatsink (floor); the receipt names the budget. *Size:* S.

**B4. Facet recovery for the rest of real geometry.** *State:* triangular facets accept any coplanar tiling; non-triangular loops keep the old path; near-duplicate vertices adopt within the chord tolerance; each pass recomputes the face set (O(tets) per facet-round — fine at 10² facets, not at 10⁴). *Work:* (a) polygonal facet loops with holes through the same tiling test and an exact point-in-polygon on the drop-axis projection; (b) incremental face bookkeeping (maintain the face set across insertions from the kernel's cavity delta); (c) corpus: comb prism (done), plate with a through-hole, enclosure box with vents, a rotated (non-axis-aligned) comb that pins the tolerance path, a 10⁴-facet tessellated cylinder for the complexity gate; (d) recovery statistics in the receipt (B3). *Falsifier:* every corpus body volumetricizes with analytic volume under the default budget; the 10⁴-facet body finishes under a stated wall-clock on the quiet host. *Kill:* a body that needs Steiner points beyond the memory-derived budget refuses with the count. *Size:* M.

**B5. Conjugate coupling falsifiers (extends s93ej.3).** *State:* derived Robin rows, fixed point (12 iterations on the heatsink), energy-balance gate, card-domain refusal with a named negative case. *Work:* the coupling residual becomes its own budget term; a metamorphic pair (fan speed × 1.3 → Re × 1.3 → refusal at the Hausen ceiling or the coefficient change the card's exponent predicts); a permutation injection with two channels fed by two branches (swapping branch names must move per-region heat rates — sum tests are blind to this); cross-ISA replay of the receipt on yto once its disk is freed (G3). If the fixed point ever needs acceleration, `fs_couple::AitkenRelaxation` / `iterate_aitken` exist; IQN-ILS does not and is not promised. *Size:* M.

**B6. Guaranteed discretization bound (q61wp.11).** *State:* open; fs-feec's `whitney` module carries the discrete spaces; fs-conduction has the P1 solve and adjoint; `fs-dwr` exists (API to verify at execution). *Work:* equilibrated-flux a posteriori bound (Braess–Schöberl / Ern–Vohralík) with an RT0 flux reconstruction on the B2 mesh; the QoI's discretization term becomes a Verified-colour bound under the modelled assumptions. *Falsifier:* the bound brackets the fine-reference error on the B2 ladder at every rung; an under-resolved mesh gives a wide bound, never a narrow lie; the bound's efficiency index is reported. *Depends:* B2. *Size:* L.

**B7. Falsifier battery (q61wp.14).** *Work:* (1) independent consumer — `python/frankensim/client.py` runs `run` and reads the verdict; `fs_checker::check_against_root` re-derives the colour from the package without the solver; (2) hostile twins — foam card flips the verdict; missing `--materials` refuses at material-resolve; a duty of 2.0 refuses at validate; stripped conduction refuses at conduction by name; (3) physical anchor — junction temperature vs the Level-A hand calculation (fin efficiency × Hausen h; B9) within the stated budget; (4) determinism — same inputs → same run identity and byte-identical report twin across two ledgers and two ISAs. *Receipt:* one JSONL per falsifier retained by the lane. *Depends:* B2 (anchor tolerance), B5. *Size:* M.

**B8. `compare` verb and recompute skipping.** *State:* `run`, `report`, `package` real; `compare` absent; `fs-recompute` exists. *Work:* `frankensim compare <run-a> <run-b> <ledger>` diffs QoIs, budgets and stage receipts by identity; fs-recompute skips unchanged stages. *Falsifier:* comparing a run with itself is empty; comparing the foam twin lists the flipped verdict and the changed material identity and nothing else. *Size:* M.

**B9. Level-A hand-calc row for the heatsink — landed 2026-09-02 (q61wp.48).** *State:* `thermal-a-heatsink-fin-array-ntu` (fs-vvreg Level-A, family fin): the lumped NTU chain from the retained operating point — C_air = ρ V A c_p, ε = 1 − exp(−h A_w / C_air), T_s = T_in + Q / (C_air ε) — gives **301.982 K**; the solver's retained maximum is 301.996 K (min 301.966 K); the straight-fin tip deficit (η = 0.9947 at k = 167 W/(m K)) is 0.047 K; acceptance atol 0.1 K, rtol 0. This is a same-inputs consistency anchor (same card, same h); the card's own uncertainty (developing-flow correlations ±20 % in h → ±1.6 K) is the budget term B5 measures, not this row's envelope. Cited by B7(3). README counts 20 Level-A rows.

**A13. A run whose ops exist but are inadmissible was reported as "unknown run" — landed 2026-09-02 (q61wp.63).** *State:* the loader now names every discovery predicate and refuses with `cli-solve-run-incompatible` ("8 retained operation(s) … none is admissible under driver version 11: 7 × driver version unsupported; 1 × unknown stage") when ops exist; `cli-solve-unknown-run` only when no op carries the id. The falsifier showed the reproduction ledger had been retained under driver version 10; a same-driver additive receipt bump (B3's v4) does not make old runs inadmissible, so no compatibility rule beyond "a driver bump retires old runs, and the message says so" is needed today. *Open:* a checked-in old-driver fixture for the test suite (the scratch ledger is the only one).

**B8a. Export verbs take 60–130 s on a 7-stage ledger (new, q61wp.62).** *State:* measured in both shell lanes (65 s / 133 s at moderate load; 107 s / 123 s at load 860). *Work:* instrument `load_completed_run`, make export O(run) not O(ledger) without skipping validation, re-measure on a quiet host, pin with headroom. *Size:* M.

**B10. The honest lane and the L3 promotion (q61wp.13, in flight).** *Work:* the solve-stage lane retains the receipt with the `stages` array at a stated HEAD, built on yto; `capability-maturity.json` moves `thermal.conduction-solve` to L3 citing that lane + receipt; spine ratchet/metrics and program-metrics regenerate; program-metrics rows (L3 count, error-budget completeness, decision-turnaround) flip from NO-DATA to measured. *Depends:* A2, A8. *Note:* this L3 asserts "the pipeline executes conduction end to end with retained receipts", not "the answer is converged" — B2/B6/B7 are what turn the colour. *Size:* S.

**B11. Transient and multi-region follow-through (after L3).** *State:* fs-conduction has linear transient and adjoint; the CLI has no transient stage; multi-region contact is declared and tested (g3 pair). *Work:* staged transient behind an explicit project flag (extreal 5.11–5.13); the enclosure example with declared contact. *Gate:* only after B10. *Size:* L.

**B12. Order-canonical meshes (research item, not on the promotion path).** *State:* R7 — tilings follow the kernel's index-ordered symbolic perturbation, so region/vertex order changes the mesh. *Work (optional):* a geometric tie-break (Edelsbrunner–Mücke simulation of simplicity keyed on a canonical coordinate order rather than insertion index) so identical geometry in any declaration order yields the same complex. *Falsifier:* the g3 reverse-order test regains the exact slot swap. *Kill:* if the kernel change costs more than a week or breaks the exact audit, keep R7 and the contract sentence. *Size:* M–L.

## 6. Workstream C — Geometry import and the physical corpus

**C1. STEP path parity with STL.** *State:* STEP import exists behind the AP242 minimum slice; the heatsink goes through STL. *Work:* one STEP body (the heatsink as a B-rep) through import → assignment → volumetricization with the same receipt; refuse non-manifold or unit mismatch by name. *Falsifier:* the STEP heatsink's QoI matches the STL heatsink's within the discretization term. *Size:* M.

**C2. Assignment selectors on real bodies.** *State:* half-space and named-group selectors; a chip footprint cannot yet be a region. *Work:* face-set selectors from STEP named faces; overlapping assignments keep refusing. *Size:* M.

**C3. Corpus of tracked bodies.** heatsink (done), enclosure with vents, plate with hole, PCB slab (anisotropic card from f85xj.5.6), rotated comb; each with a generator script, a fingerprint pinned in its `.fsim`, and a row in the freshness lane. *Size:* M (spread).

## 7. Workstream D — Journey B: the P2 marquee, once, for real

**D1. Un-gate and make `run_study` a loop (q61wp.16).** *State:* `fs_marquee::study::run_study` exists; lib compiles; `mq_004` falsifies the gradient sign (A11). *Work:* fix the sign against a finite-difference gate; N iterations of compliance decrease under a volume constraint on the 2-D quadtree CutFEM; level-set update with topological derivative and `fs_adjoint::sobolev_smooth`; DWR certificate per iteration; ledger + golden. Goodhart guard: re-solve the final design one rung finer and on a body-fitted mesh; the certificate carries both deltas. *Falsifier:* the adjoint-FD gate; a run with the sign flipped fails to descend and the test notices. *Kill:* if compliance does not decrease monotonically on the cantilever after the sign fix, stop and report the mechanism. *Size:* M–L.

**D2. Certified 3-D cut quadrature (q61wp.18).** *Work:* Saye-style dimension reduction with roots isolated by `fs_ivl::newton_roots_bounded` (certified `RootBox`es); ghost penalty; assembly via fs-sparse; p-MG/AMG from fs-solver. *Falsifier:* the enclosure of a sphere's volume tightens with the level; a deliberately wrong root box is caught by the interval check. *Size:* L.

**D3. 3-D marquee on `fs_rep_sdf::VdbGrid` (q61wp.19), `frankensim study` verb (q61wp.20), marquee falsifiers (q61wp.17).** Sequence D1 → D2 → D3. **D4 (stretch, not v1):** persistent-homology density constraint on the level set.

## 8. Workstream E — Portfolio and graph

**E1. Parking executed (q61wp.26, done).** Un-park criterion on each root; review date 2026-12-01.
**E2. Graph hygiene (q61wp.27).** 164 in_progress: finished-but-unclosed vs abandoned (close with evidence or release); P0 re-tier to zero-slack beads on the tropical critical path to Journeys A/B (`xtask generate-tropical-path` / `check-tropical-path`); close epics with no live child. *Done when:* `br ready` ≥ 10, every item on a journey path; no in_progress bead > 7 days without a live claim (Agent Mail or a note).
**E3. Bead template hygiene.** Every implementation bead carries goal, falsifier, seams, receipt, logging; beads without a falsifier return to plan space.

## 9. Workstream F — Performance (O2 executed as re-baseline)

**F1. Measurements (q61wp.28).** SpMV with the corrected STREAM denominator on x86; cancel latency measured by the executor's own timer; a citable x86 GEMM row; all on a quiet host with fingerprints.
**F2. §14.1 page.** Bars restated as measured values with dated receipts; original bars kept as stretch with rationale.
**F3. Kernel work (after Journeys A/B are L3).** Batched dense → SIMD-across-elements on fs-soa layouts; all-core GEMM → per-core L2 blocking and two-level packing; 3-D FFT → cache-blocked pencil transposes; each with a roofline receipt and a kill criterion (no more than one week without a measured gain).

## 10. Workstream G — Constellation, build and fleet hygiene

**G1. Constellation advance (q61wp.39).** Seven siblings are strict fast-forwards; advance through the compatibility train, relock, rerun the trust-cone gates; record the per-sibling delta.
**G2. Shared build hygiene.** The global `CARGO_TARGET_DIR` lock serializes every agent; document per-agent private targets; the rch lane stays the default for root-workspace builds; document the exit-144 client kills and server-side survival.
**G3. yto disk.** The verification host is at 100 %; the owner decides deletions (a 66 GB disposable target of mine is listed); until then no yto verification.
**G4. Native tooling.** The tracked `target/debug/xtask` is a Linux ELF on the Mac; document the native build (`cargo lbuild -p xtask` in a private target, 8 s incremental) and the macOS first-exec wedge diagnostic.

## 11. Workstream H — Shipped surfaces: done-lines, not drift

**H1. Wright Flyer (q61wp.29).** Redeploy current dist; measure V-14 on one qualified device; decide the wake promise.
**H2. Music (q61wp.24, .30).** Owner adjudicates the seven listening receipts; one instrument exposed in fs-wasm and the Apple catalog; jet receipts in-tree.
**H3. Apple (q61wp.32).** Metadata restored, XCTests run, three form factors with screenshots.
**H4. Website (q61wp.33).** Redeploy with current fs-wasm.
**H5. Euler cinematic (q61wp.31).** One 1080p daily frame with receipt, or park.

## 12. Workstream I — The frontier, wedge first, with kill criteria

**I1. LBM for the wedge (q61wp.34).** D3Q19 + D3Q7 thermal: promote the ignored release lanes (`cylinder_re100`, `d3q19_battery` bead 84hv, `d3q19_boundaries` bead 40p2) to run in the release profile; Boussinesq conjugate channel as the Level-B reference for the correlation cards Journey A uses (the Hausen Nusselt at Re ≈ 1.5 × 10³ checked against the LBM channel). *Kill:* if the D3Q7 thermal port cannot hold the analytic Poiseuille–Graetz Nusselt within 5 % at a stated resolution, the card stays Level-A-only and the plan says so.
**I2. 3-D incompressible NS (q61wp.35, gated).** Opens only after I1 produces a converged conjugate channel; retire criterion on the epic.
**I3. fs-instrument split (q61wp.37).** fs-couple becomes a coupling crate.
**I4. FrankenScript executor v0 (q61wp.38).** Lower an admitted study to the stage pipeline; otherwise the README says "IR and admission only".
**I5. Retired epics (q61wp.40–.43).** Dated notes in plan and README; un-retire only with a G1 ladder and a consuming journey.

## 13. Sequencing and critical path

```
week 0   A2* A3 A4 A6 A9 A10 A12     (* NobleLion, in flight)   B10*
week 1   B2 B3 B4(a,c)   A5 A7 A11        D1               E2   G1 G4
week 2   B5 B7 B9        C1               D1 loop+golden   F1
week 3   B6  B8          C2               D2               H1 H3 H4
week 4+  B11 B12 C3      D3  I1           F2 F3 (gated)    H2 H5 I3 I4
```
Critical path to the first legitimate L3 (execution claim): A2 → A8 → B10 — days. Critical path to a *verified* answer (the §3 step 2): B2 → B6 → B7 → A5. Journey B: A11 → D1 → D2 → D3. E/F/G/H are parallel and mechanical once decided.

## 14. Risk register

| Risk | Signal | Mitigation |
|---|---|---|
| Sweeper commits half-verified work (it committed tonight's slice mid-verification) | commits under the owner's name minutes after an edit | run suites before editing shared files; note "verified by" in bead closes, not in commit messages |
| Two agents edit the same seam (solve.rs, spine files) | Agent Mail/SendMessage claims | claim before editing; pathspec commits; the receipt-head rule forces sequencing |
| Refinement (B2) blows the memory budget on real bodies | Steiner count vs `max_tets` | refuse with the count; the budget is the project's declaration |
| Verified-colour bound (B6) is not tight enough to be useful | efficiency index ≫ 10 | report the index; the colour is still Verified, the width is honest |
| Marquee never descends after the sign fix | compliance non-monotone | D1 kill criterion; report the mechanism |
| LBM thermal port misses the Graetz check | > 5 % at stated resolution | I1 kill criterion; cards stay Level-A |
| Docs drift again after A5 | `check-docs` red | doc-facts inventory is the only source of numbers |

## 15. Falsification and logging standard (applies to every bead)

- **Falsifier named** in the bead body: the exact command and the expected refusal or delta.
- **Receipt retained**: path, schema, digest; the receipt names the commit and host fingerprint.
- **Logging**: stage, input digests, verdict, budget terms, wall-clock — the lane's JSONL alone must locate a failure.
- **Tolerances measured**: value, source run, headroom factor.
- **No `#[ignore]` on a true defect**; a red test gets an owner bead with the failure text.
- **Kill criterion** for every research-grade item (B12, D1, D2, F3, I1).

## 16. Bead mapping

Existing steering beads carry A1–A7 (q61wp.1–.6), B1/B5/B6/B7/B10 (q61wp.8 closed, .9, .11, .14, .13 closed at 974b1cc6), D1–D3 (q61wp.16–.20), E1–E2 (q61wp.26 closed, .27), F1 (q61wp.28), G1 (q61wp.39), H (q61wp.29–.33), I (q61wp.34–.43; .36 closed). **Beads created from this plan (Phase 3a, 2026-09-02, label `bridge-plan`):** B2 = q61wp.44, B3 (+A12) = .45, B4 = .46, B8 = .47, B9 = .48, B11 = .49, B12 = .50, C1 = .51, C2 = .52, C3 = .53, D4 = .54, E3 = .55, F2 = .56, G2 = .57, G3 = .58 (owner decision), G4 = .59, A9 = .60, A10 = .61; A11 is a note on .16. Dependencies as in §13 (C2→C1, C3→B4, F2→.28, .14→B9 and B2, .11→B2, B11→.13 and B2, D4→.19); no cycles. Refinement rounds R1–R4 live as dated notes on the beads.

## 17. What this plan deliberately does not do

- It does not re-beadify the 880 KB of plans; parked work stays parked until a journey consumes its solver.
- It adds no harness beyond R1 and the gates that exist.
- It deletes nothing; fabricated files are rewritten in place; deletions wait for the owner.
- It does not promise physics that has no solver: retired items stay retired until a G1 ladder exists.

## 18. Round log

- **Round 1 (2026-09-01, late):** first tracked version; superseded §8 of `docs/REALITY_CHECK_2026-09-01.md`; added workstream C and the B2/B3/B4 mesh items; added A8/A9.
- **Landings after round 3 (2026-09-02, 00:00–02:30):** A8 (gap-table truth), A9 (vertical-profile fabrication removed), A10 (lane registration), B3 (+A12: receipt v4 with mesh quality and recovery evidence; measured 0.000° dihedrals and 9 flat tets on the comb), B9 (heatsink NTU anchor 301.982 K vs solver 301.996 K), G2/G4 (build-hygiene section), A5 partial (README command reference, three example READMEs, the quickstart rewritten around the seven stages, the run-id capture, the exports and the three named refusals); q61wp.13 (first receipt-backed L3, NobleLion); closes on evidence: .6, .8, .10, .12 (blocked by 6.9/6.10 owners), .26, .36, .45, .48, .57, .59, .60, .61, s2l9v. New findings filed: A13 (q61wp.63, schema-drift runs reported as unknown), B8a (q61wp.62, 60–130 s exports). Not committed by me and not mine: a locally resolved `Cargo.lock` recording frankensqlite 0.3.15 against the 0.3.14 pin (G1 territory).
- **Round 3 (2026-09-02, early):** added §2.1 — eleven mathematical instruments, each bound to an item and a falsifier (DWR goal adaptation, equilibrated flux, conforming-Delaunay termination, simulation of simplicity, sliver exudation, e-process acceptance, interval-isolated Saye quadrature, ghost penalty, topological derivative, tropical critical path, IQN-ILS as a gated fallback); B2 upgraded to goal-oriented adaptation with the uniform ladder as control; B2's falsifier gains the ≤ ¼-tets and effectivity-index checks.
- **Round 2 (2026-09-02, early):** A8 landed and recorded; verified every cited seam by name (Aitken not IQN-ILS; `newton_roots_bounded`; `VdbGrid`; the Python client path; the LBM release lanes and their beads); added R7 and B12 (order dependence), A11, A12, D4, the coverage matrix (§3.1), kill criteria, the risk register (§14), the receipt-head sequencing with q61wp.13.
