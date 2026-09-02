# Comprehensive Bridge Plan for FrankenSim

**Status:** round 1, 2026-09-01 (late). Revised in place; rounds are logged in §12.
**Inputs:** `docs/REALITY_CHECK_2026-09-01.md` (Phase 1: where the code really is), the steering epic `frankensim-rc-root-q61wp` (.1–.43, label `reality-check-2026-09`), the five owner decisions taken on 2026-09-01, and the landings of the same evening (report stage, `run`, the finned heatsink solving seven stages, the fs-mesh facet-recovery fix).
**Purpose:** the Phase 2 document of the reality check — a plan that closes **every** gap between `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md` (plus its addendum) and the code, granular enough that each item becomes a self-contained bead and no reader needs this document afterwards. It is the measuring stick; the code is the ground truth.

---

## 0. How to read this plan

- **Item IDs** are stable (`A1`, `B4`, …). A bead that implements an item names it. Existing steering beads are cited as `q61wp.N`.
- Every item carries: **goal served** (plan section), **state** (what the code does today, with evidence), **work** (what changes), **seams** (what already exists and is reused rather than rebuilt), **falsifier** (the test that would fail if the item were faked), **receipt** (what is retained), **logging** (what a failing run must print so it is diagnosable without a debugger), **depends**, **size** (S = a session, M = days, L = weeks).
- **No claim advances on prose.** An L3, a README fact, a critical-path bead close each resolve to a retained, digest-pinned receipt from an executed lane (§2). That rule is the plan's spine; everything else is work.
- The plan **does not proliferate process**: no new dashboards, analyzers or harnesses beyond the receipt rule and the gates that already exist. Where a lane already exists it is repaired, not duplicated.

## 1. The gap, in one page (state at the end of 2026-09-01)

**What the vision promises.** A certified simulation product: a user declares a design study, the product either answers with an evidence colour and an error budget or refuses by name; the continuum geometry → physics → adjoint → optimizer → geometry is demonstrated; the fluid frontier (NS, LBM, VPM, compressible, turbulence) and structural frontier (IGA shells) exist at 3-D with G1 ladders; performance meets §14.1 bars; every claim traces to a receipt.

**What the code delivers today (measured, not described).**

| Journey / promise | 2026-09-01 morning | 2026-09-01 night |
|---|---|---|
| `frankensim run` on a real body | stopped at the QoI gap; report/package were fabricated literals | **all seven stages execute on the finned heatsink** (`examples/heatsink-fan`); report HTML + JSON twin + format-9 package sealed in the ledger; `report`/`package` export the retained bytes; 65-check lane green, 22-check examples lane green |
| Conduction on real geometry | proven only on 4-facet tetrahedra; the heatsink STL was five glued boxes; any comb refused at facet recovery | fs-mesh facet recovery accepts any coplanar tiling and iterates to a fixed point; one/two/four-fin combs volumetricize with exact volume under the default budget |
| Conjugate solid/air coupling | declared coefficients only | derived Robin rows from the flow-network operating point through `fs-convection` cards with domain gating (Hausen in domain at Re ≈ 1.5 × 10³ on the heatsink; balance 1e-10 W) |
| QoI with 8-term budget | unreachable | executes; all eight terms are honest NO-DATA → Estimated / indeterminate |
| L3 capabilities | one, falsely certified | zero (demoted); promotion now waits on the receipt rule (§2) and the falsifier battery (§4) |
| Marquee (P2) | gated, lib did not compile | lib compiles; one test shows the gradient sign is wrong (q61wp.16 owner) |
| Frontier (NS, LBM 3-D, VPM, compressible, turbulence, IGA) | absent or 2-D | unchanged; IGA/turbulence/compressible/FMM-VPM **retired from v1** by owner decision; NS gated behind the LBM wedge |
| Graph | 2 ready, 83 % blocked, 204 in_progress | 41 ready, 164 in_progress, 1,173 beads deferred to 2026-12-01 |

**Why the gap persisted** (both mechanisms are still partly live and the plan removes them): claims were self-certified prose, and nothing gated scope on a delivered journey. The remedy is not more beads; it is two journeys taken to L3 under a receipt rule, then the frontier reopened wedge-first with kill criteria.

**What remains after tonight, by weight.** (i) Journey A is *running*, not *verified*: mesh quality and convergence, falsifiers, an independent consumer, a physical anchor, the honest lane and the L3 promotion are open. (ii) Journey B has a compiling library and a falsified gradient. (iii) Truth repairs T1–T7 are partly landed and must be closed with evidence. (iv) The shipped surfaces need done-lines. (v) Performance and constellation decisions need their receipts. (vi) The frontier needs the LBM wedge before anything else.

## 2. The one rule and the standing doctrines

**R1 — Claims read receipts.** A capability level, a README fact, or a critical-path close resolves to a retained receipt from an executed lane that names the stages it ran, with digests, on a stated commit. Parts that exist: `spine-e2e-summary.json`, `spine-ratchet.json`, `check-maturity`, `check-docs` + `doc-facts-inventory.json`, `check-spine-metrics`, `check-suite-receipt`. Delta: check-maturity requires an L3 entry to cite a receipt whose stage list contains the capability's stage and whose digest is newer than the lane script (q61wp.2); the DSR quality gate includes a no-run build of every workspace test target (q61wp.1).

**R2 — Falsifier before feature.** Each item names the test that fails if the item is faked (a hostile twin, a permutation injection, a domain-boundary probe, a cross-ISA replay). A feature bead without a falsifier does not close.

**R3 — Measured tolerances carry their measurement.** A gate set above a measured floor states the measurement and keeps ~5× headroom; when a fixture is re-dimensioned the measurement is redone in the same commit.

**R4 — Refuse by name, never approximate.** Outside a card's domain, without a declared input, on an open shell: refuse with a code, a message and a fix. Tonight's examples: Gnielinski refusing L/Dh < 10, the volumetricizer refusing the five-box soup, `report` refusing an unknown run without writing.

**R5 — Delete nothing; rewrite in place; park reversibly.** `br defer` with a written un-park criterion is the only portfolio instrument.

**R6 — No process porn.** No new harness, dashboard or registry unless it replaces a Sev-0 mechanism observed in the audit. Existing lanes are repaired.

## 3. Definition of done for v1 (the acceptance test of the whole plan)

An outsider with the repository, the constellation at pins, and one Mac or Linux box can, in order:

1. Build the workspace with zero red test targets (`check-suite-receipt` 0 violations).
2. Run `frankensim run examples/heatsink-fan/heatsink-fan.fsim ledger.db --materials aa6061.fsmcdpk` and read a report whose junction-temperature QoI carries a colour, an 8-term budget with at least the discretization and coupling terms **measured** (not NO-DATA), a mesh-convergence statement, and a Verified-colour discretization bound (B6).
3. Re-derive that verdict from the package with the solver-free checker, and from the Python client, and get the same colour (B7).
4. Swap the aluminium card for the foam card and watch the verdict flip; omit `--materials` and watch material-resolve refuse (B7).
5. Find `thermal.conduction-solve` at L3 in the maturity registry citing that lane's receipt (B10) and nothing else at L3 without a receipt.
6. Run `frankensim study` on the tracked marquee project and get a topology-optimized design with a composed certificate and a Goodhart rung-climb delta (D-workstream).
7. Read a README whose every quantitative statement is generated from or checked against `doc-facts-inventory.json`, whose quickstart runs the heatsink end to end, and whose capability table matches the registry.
8. Open `br ready` and find ≥ 10 items, every one on a journey path; find no in_progress bead older than 7 days without a live claim.
9. Find §14.1 either met or re-baselined with dated receipts on the same page.
10. Find every shipped surface (Wright Flyer, instruments, Apple app, website, Euler cinematic) with a done-line and a receipt, or a dated park.

## 4. Workstream A — Truth and gates (close T1–T7 with evidence; add what tonight exposed)

**A1. Four committed-red test targets + DSR no-run build.** *Goal:* R1. *State:* q61wp.1 in progress; fs-scenario `rans_card_gates`, fs-feec `differential_characters`, fs-session `snapshot_freeze_gate`, fs-evidence-runner `value.rs` were red on 08-23…08-27; several are since fixed by other lanes (verify at HEAD). *Work:* make each target build; add `cargo test --workspace --no-run` (via the rch lane) to the DSR quality gate; record the gate line in `docs/CI_GATES.md`. *Falsifier:* revert one fix locally and watch the gate refuse. *Receipt:* suite receipt at HEAD with 0 violations. *Size:* S–M.

**A2. Receipt rule in check-maturity.** *Goal:* R1. *State:* L3 demoted (q61wp.2); the rule itself may be partially landed — verify the checker reads a receipt path per L3 entry. *Work:* L3 entries carry `receipt`, `stage`, `receipt_digest`; the checker verifies the receipt exists, lists the stage as executed, and its digest is newer than the lane script; refusal text names the missing field. *Falsifier:* an L3 entry citing the old cooling_01 receipt (which never ran conduction) must fail. *Depends:* none. *Size:* S.

**A3. Reopened beads with code conditions.** *State:* s2l9v reopened, QoI now executes (f47411c7) — close s2l9v on the lane receipt listing `qoi` executed (tonight's 65-check lane does). f2jag (FrankenSQLite cascade) is ready with a live test; keep open until its three tests are green. *Work:* close s2l9v with the receipt path in the close reason; leave f2jag to its owner. *Size:* S.

**A4. Fabricated renderers replaced — done tonight; close with evidence.** *State:* `report.rs`/`package.rs` are ledger-traced exports; the workflow test calls the verbs with operands; `run` chains seven stages. *Work:* close q61wp.3 and q61wp.12 citing `scripts/ci/solve_stage_producers_e2e.sh` (65 checks, cross-ledger byte-identical twins) and the two `g1_*` CLI tests. Remove the stray `cooling_demo_run.fspkg`, `profile_run.fspkg`, `feedface…fspkg` at the repo root **only with the owner's written permission** (they are droppings of the old stub verbs). *Size:* S.

**A5. README / QUICKSTART / plan §16 truth pass.** *State:* q61wp.4 open; the README still says the pipeline stops at conduction for the heatsink; `.fsim` is now v4; L3 count 0; the quickstart's steps 4–6 must now expect completion. *Work:* rewrite "what stops where", the capability table, the schema version, the quickstart steps (heatsink end to end, report/package export paths, unknown-run refusal), the "9 standalone workspaces" line; add a status column to plan §16; regenerate `doc-facts-inventory.json`. *Falsifier:* `check-docs` at HEAD; a hand run of the quickstart on a clean clone. *Depends:* A4, B-items landing dates (write "as of <commit>"). *Size:* M.

**A6. Derived registries at HEAD.** *State:* program-metrics, source manifest and SPDX were stale; tonight adds tracked files (`crates/fs-mesh/tests/comb_prism.rs`, `examples/heatsink-fan/generate_heatsink_stl.py`, `docs/REALITY_CHECK_2026-09-01.md`, this plan). *Work:* stage first, regenerate both registries (`xtask regenerate-source-manifest`, program-metrics), regenerate `spine-e2e-summary.json`/`spine-ratchet.json` from tonight's lane (7/7, first gap none), fix the Info.plist metadata, untrack the fs-cmaes-viz-wasm debug files (owner permission), gate fs-g1-train's serde. *Falsifier:* `check-source-manifest`, `check-program-metrics`, `check-spine-metrics` all 0. *Size:* S–M.

**A7. Contact test at HEAD.** *State:* q61wp.6 ready; `g0_conduction_stage_executes_declared_card_backed_contact` refused with `project-conduction-interface-undeclared`; tonight's suite run will say whether it still fails. *Work:* fix the coincident boundary-slot orientation or the fixture's interface declaration; keep the test as the only executable proof of declared contact. *Size:* S–M.

**A8. Gap message honesty (new).** *State:* a project with no `conduction` section refuses with `cli-solve-stage-gap … cannot execute until frankensim-s93ej supplies its authoritative producer` (exit 5). That message is now false: conduction executes when declared. *Work:* split the code: undeclared conduction → `cli-solve-conduction-undeclared` (exit 4, fix: "declare cooling.conduction with a seeded region and a boundary"); keep `stage-gap` for stages that truly have no producer. Update the CLI test `g0_run_stops_at_the_conduction_gap_when_the_project_declares_no_conduction` to the new code. *Falsifier:* the stripped-heatsink project. *Owner:* the solve.rs seam (NobleLion's lane). *Size:* S.

**A9. Fabricated attribution in the vertical-profile lane (new).** *State:* `scripts/ci/e2e_extreal_vertical_profile.sh` synthesizes a "PipelineAttributionReceipt" whose top-three kernels are literal fractions of the solve time (`t_solve * 0.45`, `0.35`…) with invented arithmetic intensities. That is a fabricated receipt. *Work:* either measure (per-stage wall from the retained stage receipts' `wall_s`, which exist; kernel rows only from a real profiler run) or delete the kernel rows and label the receipt "stage wall-clock only". *Falsifier:* the receipt schema gains `measured: true|false` per row and the checker refuses `false` rows in a promotion context. *Size:* S.

**A10. Freshness lanes have owners.** *State:* `examples_freshness_e2e.sh` (22 checks) and the vertical-profile script are not run by any xtask gate; the G0/G1 CLI battery is the executable proof. *Work:* register the freshness lane in `docs/CI_GATES.md` as a native-binary lane with its expected check count and make `check-suite-receipt` aware of it, or state plainly that the CLI battery is the gate. No new harness. *Size:* S.

## 5. Workstream B — Journey A: Cooling 0.1 from "runs" to a legitimate L3

**Target (unchanged):** `frankensim run examples/heatsink-fan/heatsink-fan.fsim` yields a junction-temperature QoI with a colour, a measured 8-term budget, a report and a checker-accepted package — and an outsider can falsify it.

**B1. Heatsink example declares conduction — landed tonight (q61wp.8).** *State:* v4 project, single-shell 108-facet STL from a tracked generator, corrected vent arithmetic (0.00072 m²), Hausen developing-flow card in domain at 70 % fan speed, seven stages, 22-check lane. *Remaining work:* cooling-enclosure example gets the same treatment (its `plate.stl` is a tetrahedron; give it a real enclosure body and a declared conduction section); heated-plate stays the schema tour. *Falsifier:* the freshness lane. *Size:* S–M.

**B2. Mesh quality and size control after recovery (new, load-bearing).** *State:* the conduction mesh is the recovered PLC itself: 246 vertices for the heatsink. That is a valid volume, not a converged discretization; the QoI is a coarse estimate. `fs_mesh::refine` exists (radius-edge, circumcenter insertion) but skips hull-escaping offenders and is not called by the conduction stage; recovery produces some slivers. *Work:* after carve-and-label, run constrained refinement to a target edge length derived from the declared accuracy budget and the geometry's feature size (the seed's region declares it; default = smallest declared feature / 4), with the recovered walls protected (Ruppert-style encroachment rule already sketched in the contract); retain mesh statistics (vertices, tets, min dihedral, max radius-edge) in the conduction receipt. *Seams:* `fs_mesh::refine`, `exude`, `RecoveryOptions`, the audited complex. *Falsifier:* a three-rung h-ladder on the heatsink whose QoI converges monotonically with the expected P1 rate against a fine reference; a quality gate that refuses a mesh with a dihedral below a stated floor. *Logging:* per-rung vertices/tets/min-dihedral/QoI. *Depends:* B1. *Size:* M–L.

**B3. Recovery budget from the declared memory budget; receipt fields (peer asks).** *State:* `solve.rs` derives `max_tets` from `budgets.memory_bytes` but hard-codes `RecoveryOptions::default()`. *Work:* derive `max_steiner` from the same budget with the default as a floor (never below today's behaviour); write `max_steiner`, `max_depth`, the memory budget and the recovery statistics (`FacetRecoveryStats`, segment stats) into the conduction receipt so the mesh row is reproducible from the receipt alone. *Falsifier:* a project with a tiny memory budget still recovers the heatsink (floor); the receipt names the budget. *Size:* S.

**B4. Facet recovery for the rest of real geometry (new).** *State:* tonight's fix covers triangular facets with any coplanar tiling; non-triangular loops keep the old path; near-duplicate vertices are adopted within the chord tolerance; each pass recomputes the face set (O(tets) per facet-round). *Work:* (a) polygonal facet loops (holes, notches) through the same tiling test with an exact point-in-polygon; (b) incremental face bookkeeping so a 10⁴-facet STL is tractable; (c) a corpus: comb prism (done), plate with a through-hole, enclosure box with vents, a non-axis-aligned (rotated) comb; (d) retain recovery statistics in the receipt (B3). *Falsifier:* each corpus body volumetricizes with analytic volume under the default budget; the rotated comb pins the tolerance path. *Size:* M.

**B5. Conjugate coupling falsifiers (extends s93ej.3, landed by NobleLion).** *State:* derived Robin rows, fixed point, energy-balance gate, card-domain refusal (negative case landed tonight). *Work:* coupling residual as its own budget term; a metamorphic pair (double the fan speed → Re doubles → refusal or coefficient change per the card's exponent); permutation injection on branch order (two channels fed by two branches: swapping branches must change per-region heat rates); cross-ISA replay of the receipt on yto. *Seams:* `fs_airflow::conjugate`, `fs_couple` (IQN-ILS, if the fixed point needs acceleration), `fs-convection` cards. *Size:* M.

**B6. Guaranteed discretization bound (q61wp.11).** *State:* open; fs-feec has H(div)/RT0, fs-conduction has the P1 solve and adjoint. *Work:* equilibrated-flux a posteriori bound (Braess–Schöberl / Ern–Vohralík) on the conduction mesh; the QoI's discretization term becomes a Verified-colour bound under the modelled assumptions. *Falsifier:* the bound brackets the fine-reference error on the B2 ladder for every rung; a deliberately under-resolved mesh produces a wide bound, never a narrow lie. *Depends:* B2. *Size:* L.

**B7. Falsifier battery (q61wp.14).** *Work:* (1) independent consumer: the Python client runs `run` and reads the verdict; `fs-checker` re-derives the colour from the package without the solver; (2) hostile twins: foam card flips the verdict; missing `--materials` refuses at material-resolve; a duty of 2.0 refuses at validate; (3) physical anchor: junction temperature vs a Level-A hand calculation (fin efficiency × Hausen h) within the stated budget; (4) determinism: same inputs → same run identity and byte-identical report twin across two ledgers and two ISAs. *Receipt:* one JSONL per falsifier retained by the lane. *Depends:* B2 (anchor tolerance), B5. *Size:* M.

**B8. `compare` verb and recompute skipping.** *State:* `run`, `report`, `package` real; `compare` absent. *Work:* `frankensim compare <run-a> <run-b> <ledger>` diffs QoIs, budgets and stage receipts; fs-recompute skips unchanged stages by identity. *Falsifier:* comparing a run with itself is empty; comparing the foam twin lists the flipped verdict and the changed material identity. *Size:* M.

**B9. Level-A hand-calc corpus row for the heatsink.** *Work:* fin-efficiency + developing-flow Hausen hand calculation as a `fs-vvreg` Level-A row with stated assumptions and tolerance; cited by B7(3). *Size:* S.

**B10. The honest lane and the L3 promotion (q61wp.13).** *Work:* `scripts/e2e/cooling_01.sh` (or the solve-stage lane, whichever the registry cites) runs the heatsink through seven stages and retains the receipt; regenerate `spine-e2e-summary.json` / `spine-ratchet.json`; promote `thermal.conduction-solve` to L3 through A2's rule; flip program-metrics rows (L3 count, error-budget completeness, decision-turnaround) from NO-DATA to measured. *Depends:* A2, B2, B7. *Size:* S.

**B11. Transient and multi-region follow-through (after L3).** *State:* fs-conduction has linear transient and adjoint; the CLI has no transient stage. *Work:* a staged transient (5.11–5.13 in the extreal program) behind an explicit project flag; multi-region enclosure with declared contact (A7). *Size:* L. *Gate:* only after B10.

## 6. Workstream C — Geometry import and the physical corpus (the wedge needs bodies)

**C1. STEP path parity with STL.** *State:* STEP import exists behind the AP242 minimum slice; the heatsink goes through STL. *Work:* one STEP body (the heatsink as a B-rep) through import → assignment → volumetricization with the same receipt; refuse on non-manifold or unit mismatch by name. *Size:* M.

**C2. Assignment selectors on real bodies.** *State:* half-space selector covers whole bodies; named groups exist. *Work:* selector by face-set from the STEP (named faces) so a chip footprint can be a region; refusal for overlapping assignments stays. *Size:* M.

**C3. Corpus of tracked bodies.** heatsink (done), enclosure with vents, plate with hole, PCB slab (anisotropic card from f85xj.5.6), rotated comb; each with a generator script, a fingerprint pinned in its `.fsim`, and a row in the freshness lane. *Size:* M (spread).

## 7. Workstream D — Journey B: the P2 marquee, once, for real

**D1. Un-gate and make `run_study` a loop (q61wp.16, in progress).** *State:* lib compiles (fixed tonight); mq_004 shows the gradient sign is falsified — the loop cannot descend until the adjoint sign is right. *Work:* fix the sign against a finite-difference gate per iteration; N iterations of compliance decrease under a volume constraint on the 2-D quadtree CutFEM; level-set update with topological derivative and Sobolev-smoothed gradient; DWR certificate per iteration; ledger + golden. Goodhart guard: re-solve the final design one rung finer and on a body-fitted mesh; certificate carries both deltas. *Falsifier:* adjoint-FD gate; a run with the sign flipped must fail to descend and the test must notice. *Size:* M–L.

**D2. Certified 3-D cut quadrature (q61wp.18, in progress).** Saye-style dimension reduction with interval-isolated roots (fs-ivl); ghost penalty; assembly through fs-sparse; p-MG/AMG from fs-solver. *Falsifier:* enclosure of a sphere's volume tightens with the level; a wrong root isolation is caught by the interval check. *Size:* L.

**D3. 3-D marquee on VDB (q61wp.19); `frankensim study` verb (q61wp.20); marquee falsifiers (q61wp.17).** As already beadified; sequencing D1 → D2 → D3.

## 8. Workstream E — Portfolio and graph (mechanical after the decisions)

**E1. Parking executed (q61wp.26, done).** Un-park criterion on each root; review date 2026-12-01.
**E2. Graph hygiene (q61wp.27).** 164 in_progress: audit finished-but-unclosed vs abandoned (close with evidence or release); P0 re-tier to zero-slack beads on the tropical critical path to Journeys A/B (`xtask tropical-path`); close epics with no live child. *Done when:* `br ready` ≥ 10, every item on a journey path; no in_progress bead > 7 days without a live claim.
**E3. Bead template hygiene.** Every implementation bead carries the goal served, falsifier, seams, receipt, logging; beads without a falsifier are returned to plan space.

## 9. Workstream F — Performance (O2 executed as re-baseline)

**F1. Measurements (q61wp.28).** SpMV with the corrected STREAM denominator on x86; cancel latency (fs-exec) with the executor's own timer; a citable x86 GEMM row from the fleet; all on a quiet host (ts1/ts2), retained with host fingerprints.
**F2. §14.1 page.** Bars restated as measured values with dated receipts; the original bars kept as stretch with a rationale.
**F3. Kernel work (after Journeys A/B are L3).** Batched dense → SIMD-across-elements on fs-soa layouts; all-core GEMM → per-core L2 blocking and two-level packing; 3-D FFT → cache-blocked pencil transposes; each with a roofline receipt.

## 10. Workstream G — Constellation, build and fleet hygiene

**G1. Constellation advance (q61wp.39).** Seven siblings are strict fast-forwards; advance through the compatibility train, relock, rerun the trust-cone gates; record the per-sibling delta (the aggregate lock hash hides which drifted).
**G2. Shared build hygiene.** The global `CARGO_TARGET_DIR` lock serializes every agent; document per-agent private targets; keep the rch lane the default for root-workspace builds; fix the exit-144 client kills or document the server-side survival.
**G3. yto disk.** The verification host is at 100 %; the owner decides what to delete (my 66 GB `rc-check-target-20260901` is disposable). Until then, no yto verification.
**G4. macOS first-exec wedge** documented with the sample/ps diagnostic; prebuilt-xtask workaround retained.

## 11. Workstream H — Shipped surfaces: done-lines, not drift

**H1. Wright Flyer (q61wp.29).** Redeploy current dist; measure V-14 on one qualified device; decide the wake promise (fs-vpm in the browser cone or delete the sentence).
**H2. Music (q61wp.24, .30).** Owner adjudicates the seven listening receipts; one instrument exposed in fs-wasm and the Apple catalog; jet receipts in-tree.
**H3. Apple (q61wp.32).** Metadata restored, XCTests run, three form factors with screenshots.
**H4. Website (q61wp.33).** Redeploy with current fs-wasm.
**H5. Euler cinematic (q61wp.31).** One 1080p daily frame with receipt, or park.

## 12. Workstream I — The frontier, wedge first, with kill criteria

**I1. LBM for the wedge (q61wp.34).** D3Q19 + D3Q7 thermal: promote the ignored G2 lanes; Boussinesq conjugate channel as the Level-B reference for the correlation cards Journey A uses (the Hausen card's Nusselt at Re 1.5 × 10³ can be checked against the LBM channel).
**I2. 3-D incompressible NS (q61wp.35, gated).** Opens only after I1 produces a converged conjugate channel; retire criterion stated on the epic.
**I3. fs-instrument split (q61wp.37).** fs-couple becomes a coupling crate.
**I4. FrankenScript executor v0 (q61wp.38).** Lower an admitted study to the stage pipeline; otherwise the README says "IR and admission only".
**I5. Retired epics (q61wp.40–.43).** Dated notes in plan and README; un-retire only with a G1 ladder and a consuming journey.

## 13. Sequencing and critical path

```
week 0   A2 A3 A4 A6 A8 A9   (truth closes; registries at HEAD)
week 1   B2 B3 B4(a,c)  A5 A7          D1 (sign fix)        E2  G1
week 2   B5 B7 B9        C1            D1 loop + golden      F1
week 3   B6 (bound)  B8  B10 (L3)      D2                    H1 H3 H4
week 4+  B11 C2 C3       D3 I1         F2 F3 (gated)         H2 H5 I3 I4
```
Critical path to the first legitimate L3: A2 → B2 → B7 → B10 (B6 upgrades the colour but is not on the promotion path). Journey B's critical path: D1 → D2 → D3. Everything in E/F/G/H is parallel and mechanical once decided.

## 14. Falsification and logging standard (applies to every bead)

- **Falsifier named** in the bead body: the exact command and the expected refusal or delta.
- **Receipt retained**: path, schema, digest; the receipt names the commit and host fingerprint.
- **Logging**: stage, input digests, verdict, budget terms, wall-clock — the lane's JSONL alone must locate a failure.
- **Tolerances measured**: value, source run, headroom factor.
- **No `#[ignore]` on a true defect**; a red test gets an owner bead with the failure text.

## 15. Bead mapping

Existing steering beads carry A1–A7 (q61wp.1–.6), B1/B5/B6/B7/B10 (q61wp.8, .9, .11, .14, .13), B8 (q61wp.12 remainder), D1–D3 (q61wp.16–.20), E (q61wp.26–.27), F (q61wp.28), G (q61wp.39), H (q61wp.29–.33), I (q61wp.34–.43). **New beads from this round:** A8, A9, A10, B2, B3, B4, B9, B11, C1–C3, E3, F2, G2–G4, plus the fs-mesh follow-through inside B4. Phase 3a creates them self-contained (§14 template) under the steering epic with dependencies as in §13.

## 16. What this plan deliberately does not do

- It does not re-beadify the 880 KB of plans; parked work stays parked until a journey consumes its solver.
- It adds no harness beyond R1 and the gates that exist.
- It deletes nothing; fabricated files are rewritten in place.
- It does not promise physics that has no solver: retired items stay retired until a G1 ladder exists.

## 17. Round log

- **Round 1 (2026-09-01, late):** first tracked version. Supersedes §8 of `docs/REALITY_CHECK_2026-09-01.md` (which stays as the Phase 1 record). Adds workstreams C (bodies) and the B2/B3/B4 mesh items from tonight's meshing work; adds A8/A9 truth defects found tonight.
