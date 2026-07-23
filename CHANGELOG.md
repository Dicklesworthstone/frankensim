# Changelog

This is a synthesized, agent-facing changelog for FrankenSim.

Scope window: historical reconstruction from project inception on 2026-07-05 through
[`main@291c7db`](https://github.com/Dicklesworthstone/frankensim/commit/291c7dbfdc6fd366f6ee55a6dbc39f137e05afd6)
on 2026-07-10.

The unreleased documentation entry below is intentionally narrower. It records
the verified workspace and capability-maturity facts at
[`main@55379ba`](https://github.com/Dicklesworthstone/frankensim/commit/55379baf41cf1d3584c87a81ae1dde20508a4c8f)
on 2026-07-22; it is not a claim that all 1,676 intervening commits have already
been reconstructed into capability waves here.

This document was rebuilt from git history, tag/release metadata, the checked-in
Beads tracker, and the committed README/contract surface at each endpoint. It
is organized by landed capabilities rather than raw diff order.

## Unreleased — Documentation Truth Refresh (2026-07-22)

This small-scope update reconciles the public README with the current tree and
the new machine-readable maturity registry. Its evidence sources are the root
workspace manifest, the live `crates/` inventory, crate test and contract paths,
[`capability-maturity.json`](capability-maturity.json),
[`docs/MATURITY_LEVELS.md`](docs/MATURITY_LEVELS.md), and the checked-in Beads
workstream.

### Delivered documentation correction

- Corrected the workspace inventory from 126 native `fs-*` members plus
  standalone `fs-wasm` to 135 native members plus `fs-wasm`: 136 `fs-*` crate
  directories in total.
- Corrected contract coverage from 127/127 to 136/136 and the tracked crate
  integration-test inventory from 276 to 509. The README now says explicitly that a
  file count is not a green-suite claim.
- Added the registry-backed capability boundary: 14 registered capabilities,
  with 10 at L2 and 4 at L1; none is registered at L3, L4, or L5.
- Documented the registry's current lowest-tier caveat: an L1 entry can still
  carry an operative note that the capability is unbuilt or temporarily
  build-broken, so the level label alone is not a green-build assertion.
- Replaced stale prose saying neural representations were absent and external
  validation was merely "in progress." The README now records that
  `fs-rep-neural` exists but has no registered maturity claim, while an external
  production validation corpus does not yet exist.
- Qualified design-ledger and browser-flagship prose with their recorded L1
  status and known failing/build-broken boundaries instead of presenting crate
  presence as integrated capability proof.
- Added generated README repository facts and a generated capability matrix,
  both enforced by the named `xtask check-docs` gate and composed into
  `xtask check-all`.
- Added [`doc-facts-inventory.json`](doc-facts-inventory.json), a portable
  tracked-input registry that is checked exactly against `git ls-files` in a
  worktree and remains usable in RCH/archive snapshots where `.git` is
  intentionally absent. This keeps untracked working-tree files out of public
  counts without making remote verification depend on transferred Git internals.

### Active workstream

- [`frankensim-extreal-program-f85xj.16.6`](https://github.com/Dicklesworthstone/frankensim/blob/55379baf41cf1d3584c87a81ae1dde20508a4c8f/.beads/issues.jsonl)
  tracks the wider documentation-truth program. The implementation now covers
  the visible narrative, generated facts, the registry-backed capability
  matrix, a seeded stale-count regression, and `xtask check-all` composition;
  its final proof and tracker closeout are recorded on the Bead rather than
  overstated here as a repository-wide green-suite claim.

### Representative commit

- [`55379ba`](https://github.com/Dicklesworthstone/frankensim/commit/55379baf41cf1d3584c87a81ae1dde20508a4c8f)
  — add the capability-maturity registry and governed promotion events that
  supply the README's current evidence boundary.

### Post-checkpoint validation-corpus authority slice

- Added a versioned, canonical `fs-vvreg` corpus schema covering retained
  payload identity, sensors and uncertainty, geometry, environment, dataset
  partition, preprocessing lineage, context of use, license, provenance,
  retention, acceptance regimes, and A-E evidence levels.
- Added exact-partition and dimensioned context query gates. Successful seeded
  queries deliberately return numerical `NoClaim`; caller-built registries
  cannot acquire query authority.
- Registered two worked rows: an explicitly synthetic Level-B CHT fixture and
  the existing Martin-Moyce 1952 digitized square-column curve as a Level-C
  retained artifact.
- Preserved the Martin-Moyce row's evidentiary gaps instead of fabricating
  metadata: original raw records, instrument/calibration/placement authority,
  acquisition environment/date, digitization replay, redistribution authority,
  and a defensible scalar envelope are recorded as unavailable or unresolved.
  Those gaps force the physical support cap to `Estimated`.
- Added deterministic audit output that reports optional omissions and
  reason-bearing claim-authority gaps as `WARN`, while structural defects remain
  `ERROR` failures.
- This registry slice does not establish a production validation corpus or an
  L4 capability. It records what bytes and metadata exist, and makes the
  no-claim boundary executable.
- Added a retained Level-A thermal manifest with 19 versioned rows: 12
  independently recomputed analytic references across conduction, fin,
  transient, convection-limit, radiation, and contact families, plus seven G1
  primal/adjoint order targets covering P1/P2, Neumann/Robin, and nonlinear
  anisotropic cases.
- Every Level-A row is explicitly reference-only or target-only. Registry
  queries remain numerical `NoClaim`; no solver output, mesh/refinement ladder,
  adjoint run, or machine fingerprint is bound yet, so registration is not
  reported as a G1 pass.

### Post-checkpoint Level-A execution crosswalk

- Bound six canonical Level-A analytic rows directly into the executing
  `fs-conduction` fixtures: two slab fluxes, uniform-source center rise,
  rectangular affine temperature, cylindrical-shell conductance, and `mL=1`
  adiabatic-tip fin efficiency.
- Bound three canonical P1 MMS rows directly into their executing L2 ladders:
  isotropic Dirichlet, mixed Neumann, and Robin. The tests now take the
  theoretical slope and order-gate tolerance from `fs-vvreg` instead of
  repeating those target values locally.
- Added complete analytic and MMS crosswalks. The ten unexecuted rows retain
  specific gaps: missing spherical/transient/radiation/contact physics,
  convection ownership, P2 elements, a combined anisotropic nonlinear case,
  and dual-order ladders.
- Kept authority boundaries explicit. These are test-time execution bindings,
  not retained corpus receipts; solver/model envelopes remain separate from
  the catalog's formula-reproduction tolerances, no ladder or machine
  fingerprint is persisted, and `fs-vvreg` queries remain numerical
  `NoClaim` with physical support capped at `Estimated`.
- Proved the bound battery remotely at stable HEAD `bff2be3f` with
  `RCH_REQUIRE_REMOTE=1 rch exec --no-self-healing -- env
  CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test"
  cargo test -p fs-conduction --all-targets`: 39/39 tests passed (2 adjoint,
  7 analytic, 22 conformance, 8 MMS) on `vmi1227854`, RCH job
  `29943194691043348`. Cargo execution took 193.346 seconds after the
  dependency transfer; the run produced no `fs-conduction` warning.

### Post-checkpoint Level-A convection execution crosswalk

- Bound the two `fs-vvreg` `ConvectionLimit` rows directly into the executing
  `fs-convection` circular-duct limiting-case test. The 3.66 constant-wall-
  temperature and 4.36 constant-heat-flux results now use the registry values,
  context, metric, and tolerance instead of duplicate local reference literals.
- Added a complete two-row family partition check and JSON comparison verdicts
  carrying the computed value, reference, absolute error, envelope, and
  `executed-formula-limit-not-registry-receipt` authority label.
- Raised aggregate Level-A execution coverage from 9/19 to 11/19. Eight gaps
  remain: spherical conduction, lumped transient conduction, radiation,
  contact resistance, P2 primal, combined anisotropic-nonlinear MMS, and the
  P1/P2 adjoint-order ladders. No retained comparison receipt or machine
  fingerprint is added, so registry authority remains `NoClaim`/`Estimated`.
- Remote proof passed all 11 `fs-convection` targets tests with
  `RCH_REQUIRE_REMOTE=1 rch exec --no-self-healing -- env
  CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test"
  cargo test -p fs-convection --all-targets` on worker `vmi1149989`, job
  `j-29943194691043354`. The shared HEAD advanced from `d8c2c422` to
  `877790ab` during the lane without touching these owned files, so this is a
  passing dirty-worktree proof, not a same-HEAD receipt.

### Post-checkpoint Level-A spherical-shell execution binding

- Added a deterministic pole-free spherical-shell patch to `fs-conduction`.
  The radial faces carry catalog-derived Dirichlet data; the conical and
  azimuthal faces are adiabatic for the radial `1/r` solution, avoiding pole
  and seam degeneracy while still exercising the 3-D conduction assembly.
- Bound `thermal-a-sphere-shell` directly into the analytic solver battery.
  The test checks second-order L2 refinement, a separate nodal envelope, and
  heat-rate convergence after normalizing the patch by its analytic solid
  angle. The focused remote probe observed a `3.962` L2 refinement ratio,
  `0.005565 K` fine L2 error, `0.5942 K` fine nodal error, and `0.373%`
  full-shell conductance error.
- Raised aggregate Level-A execution coverage from 11/19 to 12/19. Seven gaps
  remain: lumped transient conduction, radiation, contact resistance, P2
  primal, combined anisotropic-nonlinear MMS, and the P1/P2 adjoint-order
  ladders. This adds no retained comparison receipt or machine fingerprint;
  registry authority remains `NoClaim`/`Estimated`.
- Remote proof passed with
  `RCH_REQUIRE_REMOTE=1 rch exec --no-self-healing -- env CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" cargo test -p fs-conduction --all-targets`:
  40/40 tests (`2` adjoint, `8` analytic, `22` conformance, `8` MMS) on worker
  `vmi1227854`, job `j-29943194691043361` (614.2 seconds total; 356.1 sync,
  189.3 remote command, 2.1 artifact retrieval). Shared `main` advanced from
  `22f2ab44` to `6a63d601` during the lane without touching these reserved
  files, so this is passing dirty-worktree evidence, not a same-HEAD receipt.

### Post-checkpoint Level-A lumped-transient execution binding

- Bound `thermal-a-lumped-transient` directly into `fs-time`'s existing dense
  and operator-backed first-order generalized-alpha convergence test. The
  executed model is the normalized scalar decay `u' + u = 0` at `t/tau = 1`;
  the test checks both second-order paths against the canonical `exp(-1)`
  reference and emits structured verdicts with separate formula and numerical
  envelopes. At `h = 0.05`, the dense and operator paths both delivered
  approximately `1.998043` observed order and `7.649e-5` absolute endpoint
  error against a `5e-4` solver envelope.
- Remote regression proof passed with
  `RCH_REQUIRE_REMOTE=1 rch exec --no-self-healing -- env CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" cargo test --locked -p fs-time --all-targets`
  on worker `vmi1149989`: 44/44 tests passed (2 FrankenScipy oracle, 2
  generalized-alpha probe, 11 hybrid, 5 operator-integrator, 8 SE(3), and 16
  time-battery tests; the library target contained 0 tests). RCH reported
  615.2 seconds total, including 342.169 seconds of sync, 203.191 seconds for
  the remote command, and 2.819 seconds for artifact retrieval. Shared `HEAD`
  moved from `ac4a1d252c4afb30a53003357c73c5eaec20cb13` to
  `7458f7dea9d65e01a3f17188145c855006301afe` during the lane, so this is passing
  dirty-worktree evidence, not a same-HEAD receipt.
- Raised aggregate Level-A execution coverage from 12/19 to 13/19. Six gaps
  remain: radiation, contact resistance, P2 primal, combined
  anisotropic-nonlinear MMS, and the P1/P2 adjoint-order ladders. This is not a
  spatial transient-conduction model and does not prove that an arbitrary body
  satisfies the catalog's `Bi <= 0.1` lumped-regime condition. It adds no
  retained comparison receipt or machine fingerprint; registry authority
  remains `NoClaim`/`Estimated`.

### Post-checkpoint Level-A contact-resistance execution binding

- Added an internal Robin-jump operator for exact duplicated matching-P1
  traces. Every coincident boundary-face pair must be named and bound exactly
  once; an undeclared pair refuses instead of becoming implicit perfect
  contact or an unnoticed adiabatic gap.
- Resolved positive area-specific resistance `R''` in canonical m² K/W from an
  ordered `fs-matdb::InterfaceSystemCard`. The card identity and exact
  property-use receipt travel into the integrated per-interface heat-rate
  report and replay against the source claim set.
- Added deterministic series-resistance budgets. Stable term ordering makes
  input permutations bit-identical, stated symmetric half-widths combine by a
  conservative sum with a confidence floor, finite-sum overflow refuses, and
  any unstated term keeps the complete uncertainty band unavailable.
- Bound `thermal-a-contact-series` to a G1 two-slab solve. The focused remote
  fixture reproduced the `0.3 K/W` Level-A network as
  `0.30000000000000004 K/W`, with `1.393e-11 K` maximum nodal error and
  approximately `100 W` A-to-B interface heat rate. Missing cards, missing
  bindings, changed-mesh reuse, and aggregate uncertainty overflow have typed
  refusal tests.
- Remote regression proof passed with
  `RCH_REQUIRE_REMOTE=1 rch exec --no-self-healing -- env
  CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test"
  cargo test --locked -p fs-conduction --all-targets` on worker `vmi1149989`:
  46/46 tests passed (2 adjoint, 8 analytic, 22 conformance, 6 contact, and 8
  MMS; the library target contained 0 tests). RCH reported 581.0 seconds total,
  including 346.395 seconds of sync, 163.212 seconds for the remote command,
  and 1.836 seconds for artifact retrieval. Shared `HEAD` moved from
  `90137712` to `34c57c5b` during the lane without touching the reserved
  contact files, so this is passing dirty-worktree evidence rather than a
  same-HEAD retained receipt.
- Raised aggregate Level-A execution coverage from 13/19 to 14/19. Five gaps
  remain: radiation, P2 primal, combined anisotropic-nonlinear MMS, and the
  P1/P2 adjoint-order ladders. This rung has no nonmatching/mortar projection,
  pressure/gap/temperature-dependent contact law, changing contact area,
  scenario-interface auto-lowering, contact adjoint, dominance-attribution
  corpus case, retained registry comparison receipt, or L4 authority.

### Post-checkpoint Level-A surface-radiation execution binding

- Added card-backed hemispherical-total emissivity with an explicit material
  state, immutable card identity, exact property-use receipt, validity-domain
  query temperature, and source uncertainty statement. Missing, ambiguous,
  dimensionally wrong, non-physical, or out-of-domain properties refuse.
- Added a validity-bounded surface-to-ambient linearization
  `h_rad = 4 ε σ T_mean³`. Each evaluation returns the exact Robin row consumed
  by conduction, the full `ε σ (T_s⁴ - T_a⁴)` comparison, their pointwise
  discrepancy, and first-order emissivity-width propagation. Evaluation beyond
  the declared temperature departure refuses instead of stretching the model.
- Added admitted closed-enclosure view-factor matrices with finite/range and row
  closure checks, area-weighted reciprocity, stable content identity, and
  explicit analytic or external-QMC evidence. External-QMC evidence requires a
  seed, positive sample count, and generator/version label; this crate does not
  generate or validate the rays themselves.
- Added deterministic opaque gray-diffuse radiosity and a budgeted,
  under-relaxed conduction coupling. Surface bindings retain an exact mesh
  identity, flux overlays may replace only explicit adiabatic-remainder faces,
  and the reports expose radiosity residuals, enclosure heat closure, fixed-point
  updates, and convergence thresholds. The two-surface analytic fixture
  produced `517.7298382565217 W/m²` with zero measured enclosure-energy
  residual in the focused remote run.
- Bound `thermal-a-parallel-plate-view-factor` to the executing radiation
  battery and raised aggregate Level-A execution coverage from 14/19 to 15/19.
  Four gaps remain: P2 primal, combined anisotropic-nonlinear MMS, and the P1/P2
  adjoint-order ladders. This is a test-time binding, not a retained comparison
  receipt or L4 authority.
- Remote regression proof passed at stable `HEAD`
  `b3b5f2c1c809eec06cde1e40cbc916d6995469b5` with
  `RCH_REQUIRE_REMOTE=1 rch exec --no-self-healing -- env
  CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test"
  cargo test --locked -p fs-conduction --all-targets` on worker `vmi1149989`:
  53/53 tests passed (1 library, 2 adjoint, 8 analytic, 22 conformance, 6
  contact, 8 MMS, and 6 radiation). RCH reported 632.1 seconds total, including
  356.637 seconds of sync, 199.358 seconds for the remote command, and 2.292
  seconds for artifact retrieval; the start and end SHAs were identical.
- No general-geometry view-factor/QMC generator, participating-medium or
  spectral/specular/semitransparent model, spatially varying emissivity/trace
  flux, additive radiation on an occupied boundary row, monolithic nonlinear
  operator, nonlinear uncertainty certificate, radiation adjoint,
  radiation-dominated retained corpus case, or registered maturity promotion is
  claimed by this rung.

### Post-checkpoint Level-A P2 primal MMS binding

- Bound `thermal-a-mms-p2-dirichlet` to a four-grid L2 ladder that executes
  the existing `fs_feec::highorder::simplex::SimplexSpace` order-2
  tetrahedral kernel. The fixture solves
  `-Delta T = 3*pi^2*sin(pi*x)*sin(pi*y)*sin(pi*z)` with homogeneous
  Dirichlet temperature-excess data and evaluates the discrete field through
  the kernel's high-order L2 path.
- Focused remote proof passed 1/1 on `hz1` (job
  `j-29943194691043463`) with
  `RCH_REQUIRE_REMOTE=1 rch exec --no-self-healing -- env
  CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test"
  cargo test --locked -p fs-conduction --test mms
  mms_p2_isotropic_dirichlet_order -- --nocapture`. The four L2 errors were
  `0.00566462`, `0.00167318`, `0.000704082`, and `0.000360002`, giving
  observed order `3.007794` against theoretical order `3.0`; the worst linear
  residual was `9.87e-13`. The lane started and ended at
  `7a2a63ffc1535e62452550624cb8bc719b29e836`, so this is stable-HEAD
  dirty-worktree execution evidence.
- The first complete MMS-test run on `vmi1227854` (job
  `j-29943194691043465`) passed all nine numerical and matrix checks that it
  reached, then exposed the crosswalk guard's stale expected executing-binding
  count (`5` rather than `6`). After correcting that count, the focused guard
  retry on `hz1` (job `j-29943194691043470`) exposed a second stale assumption
  that every executing target lived in `tests/mms.rs`; the already-bound P1
  adjoint correctly lives in `tests/adjoint.rs`. The guard now validates a
  canonical crate-relative integration-test function path, and its next retry
  passed 1/1 on `hz1` (job `j-29943194691043475`). These were registry-test
  bookkeeping failures, not numerical ladder failures, and remain recorded
  here rather than being folded into the passing result.
- Full-crate remote proof then passed all 56 tests on `vmi1227854` (job
  `j-29943194691043477`) with `RCH_REQUIRE_REMOTE=1 rch exec
  --no-self-healing -- env
  CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test"
  cargo test --locked -p fs-conduction --all-targets`: 1 library, 3 adjoint,
  8 analytic, 22 conformance, 6 contact, 10 MMS, and 6 radiation tests. The
  lane started and ended at `c942e7100989e3fa6bff641317350d0b9d9c74c1`;
  RCH reported the artifact retrieval in its terminal output and a successful
  remote verdict.
- Aggregate Level-A execution coverage rises from 17/19 to 18/19. Only the P2
  heat-adjoint order target remains without an executing kernel. This rung
  does not turn the P1 `fs-conduction` material/boundary frontend into a
  general P2 thermal frontend, persist a corpus receipt or ladder, or add a
  machine fingerprint or L4 authority.

### Post-checkpoint Level-A P1 heat-adjoint MMS binding

- Bound `thermal-a-mms-p1-adjoint` to a four-grid P1 dual L2 ladder. The
  fixture manufactures `z(s) = s - s^4/4` with `-Delta z = 3s^2`, homogeneous
  Dirichlet data at `s=0`, and a natural boundary elsewhere, then solves the
  transposed conduction operator through `fs_adjoint::ift_gradient_matfree`.
- Every rung also checks the discrete defining identity
  `(dJ/du)^T u_h = lambda_h^T b_h`, so an operator-transpose,
  boundary-elimination, or objective-weight mismatch cannot hide behind the
  observed-order fit. This is dual-consistency evidence for the declared
  fixture, not a DWR estimate or goal-oriented error bound.
- Focused remote proof on `hz1` (job `j-29943194691043456`) passed 1/1 with
  `RCH_REQUIRE_REMOTE=1 rch exec --no-self-healing -- env
  CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test"
  cargo test --locked -p fs-conduction --test adjoint mms_p1_adjoint_order --
  --nocapture`. The L2 errors were `0.00717563`, `0.00322824`, `0.00182664`,
  and `0.00117296`, giving observed order `1.976427` against theoretical order
  `2.0`; the worst primal/dual identity error was `5.97e-13`, and the worst
  adjoint residual was `9.88e-13`. The initial probe refused before forming a
  ladder because a primal residual of `1.725e-14` narrowly missed a `1e-14`
  solver threshold; the dedicated ladder threshold is now `1e-12`, still 100x
  tighter than the asserted `1e-10` evidence bound. The shared-tree HEAD moved
  from `e58ce6ca` to `d50e1004` during proof, so this is dirty-worktree
  execution evidence rather than a same-commit certificate.
- The broader `fs-conduction` regression command,
  `RCH_REQUIRE_REMOTE=1 rch exec --no-self-healing -- env
  CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test"
  cargo test --locked -p fs-conduction --all-targets`, did not reach the
  conduction tests on `vmi1149989` (job `j-29943194691043460`). Compilation
  stopped in a concurrent dirty `fs-matdb` edit because `query.rs:1103-1104`
  referenced `PropertyClaim` without importing it. The shared-tree HEAD moved
  from `d50e1004` to `3bee6d8b` during the lane. This is an external
  shared-snapshot blocker, not crate-wide regression evidence for this rung.
- Aggregate Level-A execution coverage rises from 16/19 to 17/19. Only the P2
  primal and P2 adjoint targets remain without executing kernels; no retained
  corpus receipt, full ladder artifact, or machine fingerprint is added here.

### Post-checkpoint Level-A combined anisotropic-nonlinear MMS binding

- Bound `thermal-a-mms-p1-anisotropic-nonlinear` to an executing P1 L2
  refinement ladder rather than treating the existing constant-anisotropic and
  isotropic-`k(T)` ladders as a substitute for the combined target. The new
  fixture uses a rotated orthotropic frame, three distinct linear principal
  conductivity laws, and the full manufactured chain rule
  `-div(K(T) grad(T)) = -(K(T):H(T) + grad(T)^T K'(T) grad(T))`.
- The fixture explicitly checks that the tensor has an off-diagonal component,
  that multiple tensor components move with temperature, and that the nonlinear
  driver takes at least two iterations. Its order gate is resolved from the
  canonical `fs-vvreg` row at runtime; the separate isotropic-`k(T)` regression
  remains in the battery.
- Focused remote proof on `vmi1149989` passed 1/1 test with
  `RCH_REQUIRE_REMOTE=1 rch exec --no-self-healing -- env
  CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test"
  cargo test --locked -p fs-conduction --test mms
  mms_anisotropic_temperature_dependent_order -- --nocapture`. The four-grid
  ladder reported L2 errors `0.0685325`, `0.0304509`, `0.0171268`, and
  `0.0109606`, for observed order `2.000463` against theoretical order `2.0`.
  The shared-tree HEAD advanced from `8550df13` to `4795745c` during the lane,
  so this is dirty-worktree execution evidence rather than same-commit proof.
- The crate-wide remote follow-up passed 54/54 tests across seven test binaries
  with `RCH_REQUIRE_REMOTE=1 rch exec --no-self-healing -- env
  CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test"
  cargo test --locked -p fs-conduction --all-targets` on `vmi1149989`; this
  includes all nine MMS tests. HEAD again moved underneath the shared-tree lane,
  from `4795745c` to `f467b119`, so the result remains explicitly
  dirty-worktree evidence rather than a same-commit certificate.
- Raised aggregate Level-A execution coverage from 15/19 to 16/19. Three gaps
  remain: P2 primal and the P1/P2 adjoint-order ladders. This observed-order
  ladder is test-time execution evidence only: it persists no registry receipt,
  machine fingerprint, or universal convergence certificate.

### Post-checkpoint convection-correlation rung

- Added `fs-convection`, an L3 library with 11 Nusselt relations spanning
  analytic duct limits, developing and turbulent ducts, rectangular ducts,
  flat plates, isolated-cylinder crossflow, and vertical-plate natural
  convection.
- Each relation carries a source, assumptions, known failures, a shared
  machine-readable validity domain, and an explicit discrepancy basis.
  Missing, non-finite, non-positive, or out-of-domain inputs refuse instead of
  extrapolating silently.
- The public dimensional conversion returns typed W/(m² K) with attached model
  evidence. The paired Robin-boundary adapter keeps that evidence beside the
  exact boundary row consumed by `fs-conduction`.
- Added Level-A limit and frozen-formula spot checks, domain-edge refusals, a G3
  coherent-unit-rescaling check, evidence-retention checks, and a
  three-flow-rate conduction seam test. The focused remote Cargo lane passed 10
  tests.
- The empirical discrepancy values are declared engineering allowances, not
  source-published confidence intervals. Source-table corpora, automatic model
  selection or blending, interrupted-fin and array-interference models,
  conjugate CFD, experimental validation, and capability-maturity registration
  remain outside this slice.

### Post-checkpoint fan/network operating-point rung

- Added `fs-airflow`, an L3 library for typed, monotone fan curves and
  series/parallel quadratic-loss networks. Admission makes pressure tolerance,
  tolerance authority, speed-scaling bounds, and the non-admissible stall
  region explicit; malformed or out-of-domain inputs refuse structurally.
- Added identical-fan series/parallel laws, exact quadratic resistance
  composition, and an enclosure model that requires leakage to be represented
  as a distinct branch rather than buried in a fitted system coefficient.
- Added a bounded interval-Newton solve for the unique nominal intersection of
  the retained fan and system curves. Only that declared mathematical root is
  `Certified`; fan-data tolerance, loss uncertainty, model discrepancy,
  pressure, branch flow, and leakage conclusions remain `Estimate` evidence.
- Added deterministic terminal-branch allocation and a typed velocity/Reynolds
  handoff into `fs-convection`, retaining the weakest evidence authority and
  refusing missing branches or invalid handoff quantities. Operating-point
  provenance binds the full curve, source and tolerance authority, fan-bank
  configuration, recursive loss-network semantics, and leakage identity.
- Added nine conformance tests covering monotone-curve refusal, fan and loss
  composition identities, a unique sign-changing nominal bracket, stall
  refusal, three speed points, leakage sensitivity and flow balance, the
  convection handoff, and provenance separation for equal-nominal models with
  different uncertainty authority. The focused remote Cargo lane passed all
  nine tests on RCH worker `ovh-a` (job `j-29943194691043331`).
- The retained fan fixture and uncertainty allowances are synthetic. No
  manufacturer curve/tolerance corpus, unequal parallel-fan model,
  installation-effects model, CFD comparison, experimental validation, or
  capability-maturity registration is claimed by this slice.

## Version Timeline

There are no git tags and no GitHub Releases as of
[`main@291c7db`](https://github.com/Dicklesworthstone/frankensim/commit/291c7dbfdc6fd366f6ee55a6dbc39f137e05afd6).

| Version | Kind | Date | Summary |
|---------|------|------|---------|
| [`main@55379ba`](https://github.com/Dicklesworthstone/frankensim/commit/55379baf41cf1d3584c87a81ae1dde20508a4c8f) | Current factual checkpoint; not a full history reconstruction | 2026-07-22 | 2,675 commits; 136 `fs-*` crate directories, 136 contracts, 509 tracked crate integration-test files, and the first machine-readable 14-capability maturity registry. |
| [`main@291c7db`](https://github.com/Dicklesworthstone/frankensim/commit/291c7dbfdc6fd366f6ee55a6dbc39f137e05afd6) | Public mainline snapshot | 2026-07-10 | 999 commits, refining the FEEC register-accumulator contractions into const-bound index loops while preserving their golden and pending both-ISA gate. |
| [`main@2568c82`](https://github.com/Dicklesworthstone/frankensim/commit/2568c8262bf50789d20be528265b5c96d575fb1b) | Public mainline snapshot | 2026-07-10 | 998 commits, adding a measured NEON complex-transpose capsule, SLP-friendly FEEC contraction accumulators, and a standalone zero-dependency clean-machine constellation bootstrap. |
| [`main@fb34ef6`](https://github.com/Dicklesworthstone/frankensim/commit/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8) | Public mainline snapshot | 2026-07-10 | 992 commits, adding governed historical roofline baselines, documenting why global x86 FMA contraction violates G5 determinism, and re-pinning FrankenNumpy past an unrepresentable case-colliding corpus revision. |
| [`main@f93c9de`](https://github.com/Dicklesworthstone/frankensim/commit/f93c9de48367ec43a1ac4ab3df7e24a1050882f0) | Public mainline snapshot | 2026-07-10 | 988 commits, adding evidence-package schema v4 and shared BLAKE3 roots, failure-compounding and repository claim/closure gates, machine-adaptive execution, fail-closed numerical repairs, and the first production GEMM autotune/roofline integration. |
| [`main@e993e76`](https://github.com/Dicklesworthstone/frankensim/commit/e993e7640a547ca9b11ded6d580a3ce6846a4c82) | Public mainline snapshot | 2026-07-10 | 857 commits, adding deterministic `.powi` policy, run-identity RNG binding, caller-owned race/session cancellation gates, solver snapshot envelopes, GEMM perf evidence lanes, legal face-first PLY import, standard empirical CVaR weighting, and adjoint certificate fail-closed regressions. |
| [`main@6725739`](https://github.com/Dicklesworthstone/frankensim/commit/6725739f42b878b310e7e7e318fb8c980cef71f3) | Public mainline snapshot | 2026-07-09 | 738 commits, adding live browser flagship pipelines, mesh v2/v3 closure, sparse roofline and NUMA lanes, rand/FFT perf work, and fail-closed IO/risk/probe hardening. |
| [`main@e08e302`](https://github.com/Dicklesworthstone/frankensim/commit/e08e30280bcd7af05ae55e990b129d6f75192ead) | Public mainline snapshot | 2026-07-09 | 696 commits, adding the gated SME2 exploratory capsule, SME2 battery, mesh hull-regression guard, and topopt proof-hygiene cleanups. |
| [`main@438128d`](https://github.com/Dicklesworthstone/frankensim/commit/438128d988082f2183f406d75883b766fd6b7324) | Public mainline snapshot | 2026-07-09 | 688 commits, adding hull-facet encroachment protection to mesh refinement and tightening the topology-optimization marquee evidence fixtures. |
| [`main@1fe4ef5`](https://github.com/Dicklesworthstone/frankensim/commit/1fe4ef5c21737a92aa58dd02f901052de5b0a88b) | Public mainline snapshot | 2026-07-09 | 686 commits, adding the flagship replay-suite scaffold, topopt marquee evidence hardening, exact certificate payload refinements, and tracker state for the fs-mesh v2 follow-up. |
| [`main@9dc7417`](https://github.com/Dicklesworthstone/frankensim/commit/9dc7417090e62e7ab2756719034e2ab0d3b32875) | Mainline checkpoint | 2026-07-09 | 684 commits, registering the staged flagship replay suite and its contract. |
| [`main@d5873bf`](https://github.com/Dicklesworthstone/frankensim/commit/d5873bfd82a3c2dbe359c11aef4947a5def8cdba) | Public mainline snapshot | 2026-07-09 | 681 commits, adding the ornithoid flagship contract, CutFEM-octree topology lane, sparse-GP and MLMC decision loops, browser campaign tiers, and proof-hardening fixes. |
| [`main@319cb64`](https://github.com/Dicklesworthstone/frankensim/commit/319cb64f052d76e15882ee53ace41092881c7fa8) | Public mainline snapshot | 2026-07-08 | 633 commits, extending the campaign suite, inverse-trig/AD surface, and p-MG smoother evidence. |
| [`main@fb08842`](https://github.com/Dicklesworthstone/frankensim/commit/fb088428ae810ee3ddd893712588c1e64f6cc0c4) | Prior changelog checkpoint | 2026-07-08 | 622 commits, including the proof-robust, schedule, and flutter end-to-end campaign capstones. |
| [`origin/main@941a67e`](https://github.com/Dicklesworthstone/frankensim/commit/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f) | Earlier public checkpoint | 2026-07-08 | 621 commits, covering the large 2026-07-08 implementation wave through value-of-information query planning. |
| [`main@43d52f2`](https://github.com/Dicklesworthstone/frankensim/commit/43d52f2ed960ede2067eec353359063c86d9dbb3) | Prior changelog baseline | 2026-07-07 | Working Rust workspace with 333 commits, no formal release tag, and the Phase 0 spine closed. |
| [`8e4c0a5`](https://github.com/Dicklesworthstone/frankensim/commit/8e4c0a5c4f18aa7d0bd0add47e407b835c7a3b86) | Inception commit | 2026-07-05 | Initial FrankenSim plan. |

## 1. Project Foundation And Policy Spine

FrankenSim started as a plan-first repository and quickly became a real Rust
workspace with project-specific agent policy, Beads issue tracking, licensing,
README/illustration assets, DSR-first validation policy, and workspace policy
checks.

### Delivered capability

- Initial architecture plan and README explaining the deterministic,
  evidence-oriented simulation substrate.
- Beads tracker and `bv` workflow documentation for project workstreams.
- MIT license with the Anthropic/OpenAI rider.
- Rust workspace scaffold with flat `fs-*` crates.
- `xtask` policy gates for layer direction, dependency policy, contracts,
  unsafe capsules, and constellation locks.
- DSR documented as the preferred validation lane over GitHub Actions.

### Closed workstreams

- [`frankensim-epic-foundations-huq.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - workspace scaffold, layer enforcement, lint/format rails.
- [`frankensim-epic-foundations-huq.2`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - pinned Franken constellation integration.
- [`frankensim-epic-foundations-huq.4`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - CI/Gauntlet gate policy.

### Representative commits

- [`8e4c0a5`](https://github.com/Dicklesworthstone/frankensim/commit/8e4c0a5c4f18aa7d0bd0add47e407b835c7a3b86) - initial FrankenSim plan.
- [`48e45c7`](https://github.com/Dicklesworthstone/frankensim/commit/48e45c7e2965db98c29e85b7bf06cf2f021eb22e) - initialize Beads project tracker.
- [`7554599`](https://github.com/Dicklesworthstone/frankensim/commit/755459986a056368685250c932407da019af7a8b) - add the MIT license with AI rider.
- [`a7e4d54`](https://github.com/Dicklesworthstone/frankensim/commit/a7e4d54d904d131850ef238f4b8e32b366c60736) - add the initial Rust workspace scaffold.
- [`4b31ce8`](https://github.com/Dicklesworthstone/frankensim/commit/4b31ce842848f40c9548a21a6513bbc050f38aa6) - document DSR as the primary CI runner.

## 2. Substrate, Execution, Determinism, And Performance

The first implementation wave established the lower layers that later crates
depend on: hardware capability discovery, explicit execution contexts,
deterministic reductions, aligned allocation, tile identity, SIMD capsules,
autotuning, and roofline measurement.

### Delivered capability

- `fs-substrate` capability probes, machine fingerprints, tiled fields, Morton
  identities, bandwidth probes, and CCD-aware topology helpers.
- `fs-simd` scalar and capsule-backed SIMD facades.
- `fs-alloc` scoped arenas, aligned allocation, hugepage policy, arena pools,
  and site tracking.
- `fs-exec` execution contexts, cancellation gates, tile pools, deterministic
  reductions, speculative races, resumable solver state, and autotuner hooks.
- `fs-roofline` machine-axis probing, kernel registry, ledgered measurements,
  and staleness checks.

### Closed workstreams

- [`frankensim-epic-substrate-wf9.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - substrate capability probes and machine fingerprints.
- [`frankensim-epic-substrate-wf9.2`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - SIMD facades and capsules.
- [`frankensim-epic-substrate-wf9.4`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - scoped arenas and allocation pools.
- [`frankensim-epic-substrate-wf9.7`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - two-lane executor and `Cx` contract.
- [`frankensim-epic-perf-fz2.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - roofline harness.

### Representative commits

- [`47c7719`](https://github.com/Dicklesworthstone/frankensim/commit/47c771994b3294d1dda1949640e4c46fc920203a) - implement `fs-substrate` capability probes.
- [`741979d`](https://github.com/Dicklesworthstone/frankensim/commit/741979d7cc8f11b9ca1e0f556d28fc7095476ca0) - implement `fs-simd` scalar and capsule tiers.
- [`59b85cb`](https://github.com/Dicklesworthstone/frankensim/commit/59b85cb47f69f2fffbf2258ff119288de05c58a8) - implement `fs-alloc` scope arenas and pools.
- [`39bd2f8`](https://github.com/Dicklesworthstone/frankensim/commit/39bd2f8e686cb97623d322fd0bd8ab96e2e8db80) - add `fs-exec` tile context and pool modules.
- [`05e3c3d`](https://github.com/Dicklesworthstone/frankensim/commit/05e3c3d6554ca094f7a8d526834f0142107ff4fb) - add deterministic reduction support.

## 3. Bedrock Numerics And Certified Arithmetic

The BEDROCK layer filled in deterministic math, random streams, dense and sparse
linear algebra, FFTs, interval/affine/Taylor models, Chebyshev tools,
automatic differentiation, e-process inference, and mixed-precision ladders.

### Delivered capability

- `fs-math` deterministic elementary functions, strict complex arithmetic, EFT,
  and double-double helpers.
- `fs-rand` logical Philox streams, deterministic distributions, QMC/Sobol
  components, and replayable sampling.
- `fs-la` GEMM/factorization kernels, eigensolver scaffolding, mixed-precision
  refinement ladders, randomized NLA, and batched small-dense kernels.
- `fs-sparse` deterministic COO/CSR/BSR/SELL assembly and sparse operations.
- `fs-fft` FFT/DCT transform kernels.
- `fs-ivl` intervals, affine arithmetic, exact predicates, Taylor models, and
  Newton/Krawczyk certification helpers.
- `fs-cheb` Chebyshev function objects and Orr-Sommerfeld spectral probes.
- `fs-ad` forward duals and adjoint infrastructure.
- `fs-eproc` anytime-valid e-process machinery.

### Closed workstreams

- [`frankensim-epic-bedrock-6ys.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - logical Philox streams.
- [`frankensim-epic-bedrock-6ys.3`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - GEMM and dense kernels.
- [`frankensim-epic-bedrock-6ys.9`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - sparse formats and deterministic assembly.
- [`frankensim-epic-bedrock-6ys.12`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - interval and affine arithmetic.
- [`frankensim-epic-bedrock-6ys.13`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - Taylor models and certified roots.
- [`frankensim-epic-ascent-7tv.6`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - e-process inference core.

### Representative commits

- [`089cc72`](https://github.com/Dicklesworthstone/frankensim/commit/089cc72cd6d00db69e272fcfae13af42b1629f19) - implement logical Philox streams.
- [`7456302`](https://github.com/Dicklesworthstone/frankensim/commit/7456302d9ae1b6aaa96be2f7f04a8aeb8a12209f) - add FFT correctness-first transform core.
- [`44a883b`](https://github.com/Dicklesworthstone/frankensim/commit/44a883b17167dc0f848b13710208f9a606c35d67) - add deterministic CSR core.
- [`e7dc872`](https://github.com/Dicklesworthstone/frankensim/commit/e7dc872f211fa66bc3d41969ebe21c4c886dc97b) - add interval arithmetic core.
- [`a433365`](https://github.com/Dicklesworthstone/frankensim/commit/a433365ab46ce18726c20cb80ff25f527507e432) - add GEMM and factorization core.
- [`a8f88f8`](https://github.com/Dicklesworthstone/frankensim/commit/a8f88f82c3aaebc0864926e9b6125d7d27f1c14c) - add evidence and model-form certification core.

## 4. Geometry, Representation, Meshing, And Shape Control

The MORPH layer became a working representation stack rather than a plan:
regions/charts, SDFs, meshes, F-reps, NURBS, voxel fields, exact topology,
representation routing, geometric constraints, IO quarantine, and
parameterization levers all landed as crate-level surfaces with contracts.

### Delivered capability

- `fs-geom` chart/region contracts, conversion records, representation router,
  semantic diff, sheaf watertightness certificates, and sheaf repair.
- `fs-rep-sdf`, `fs-rep-mesh`, `fs-rep-frep`, `fs-rep-nurbs`, and
  `fs-rep-voxel` representation crates.
- `fs-mesh` Delaunay/refinement/remeshing scaffolding.
- `fs-xform` FFD, RBF, level-set, density, and composition parameterizations
  with Jacobian actions and foldover checks.
- `fs-topo` validity and topology certificates.
- `fs-query` geometry interrogation layer and `fs-geocon` geometric constraint
  primitives.
- `fs-io` import/export with quarantine receipts.

### Closed workstreams

- [`frankensim-epic-morph-wqd.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - region/chart core.
- [`frankensim-epic-morph-wqd.2`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - SDF representation crate.
- [`frankensim-epic-morph-wqd.4`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - mesh representation crate.
- [`frankensim-epic-morph-wqd.5`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - NURBS spline algebra and trims.
- [`frankensim-epic-morph-wqd.19`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - design parameterizations.
- [`frankensim-epic-morph-wqd.23`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - topology certificates.
- [`frankensim-epic-morph-wqd.25`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - import/export quarantine.

### Representative commits

- [`397e325`](https://github.com/Dicklesworthstone/frankensim/commit/397e3255e28599bbd1ebc99dee9fd12923b26837) - add `fs-geom` region and chart crate.
- [`ebf97fb`](https://github.com/Dicklesworthstone/frankensim/commit/ebf97fb63dc4b3c921b6b08e2d7db9f8ad0c1834) - add SDF representation crate.
- [`124ed81`](https://github.com/Dicklesworthstone/frankensim/commit/124ed81ceb7a8b89283d504177b9c5739e997208) - add mesh representation crate.
- [`a295b6b`](https://github.com/Dicklesworthstone/frankensim/commit/a295b6b3f6a380e79dd4cf4ec4f9fb6e4a47b4df) - add F-rep CSG DAGs and evaluators.
- [`4a17114`](https://github.com/Dicklesworthstone/frankensim/commit/4a17114c15002faf5d400706ef8971187da26614) - add design parameterizations.
- [`3948a82`](https://github.com/Dicklesworthstone/frankensim/commit/3948a82fb0793064c3428052335c59a46d70f8fe) - add validity and topology certificates.
- [`f8a3cf7`](https://github.com/Dicklesworthstone/frankensim/commit/f8a3cf7d60c56afdc121f53f85a9cbbb9c10d9eb) - implement NURBS algebra, trims, and honest Boolean policy.
- [`3432e6e`](https://github.com/Dicklesworthstone/frankensim/commit/3432e6e6d5f6a847cab52d00eaa47afff8eba2b0) - implement quarantined import/export.

## 5. Physics, Operators, Solvers, And Time Integration

The FLUX layer gained typed scenario/regime/material foundations, exterior
calculus infrastructure, operator generation, time integration, solver
state/reporting, p-multigrid, and mixed-precision Krylov refinement.

### Delivered capability

- `fs-scenario` typed boundary/load-case algebra.
- `fs-regime` Buckingham-Pi scaling, nondimensionalization, and admission
  predicates.
- `fs-material` constitutive kernels with consistent tangents.
- `fs-feec` exterior calculus complexes and high-order tensor/simplex families.
- `fs-opdsl` typed operator IR that generates primal/JVP/VJP/adjoint/DWR hooks.
- `fs-time` structure-preserving and resumable integrator batteries.
- `fs-solver` Krylov solvers, p-multigrid preconditioning, structured stall
  diagnosis, and mixed-precision refinement.
- `fs-iface` coupling-graph static checker.

### Closed workstreams

- [`frankensim-epic-flux-tfz.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - scenario algebra.
- [`frankensim-epic-flux-tfz.2`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - material kernels.
- [`frankensim-epic-flux-tfz.3`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - regime/admission gates.
- [`frankensim-epic-flux-tfz.4`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - operator DSL/compiler.
- [`frankensim-epic-flux-tfz.5`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - FEEC complexes.
- [`frankensim-epic-flux-tfz.6`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - high-order FEEC.
- [`frankensim-epic-flux-tfz.12`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - time integrators.

### Representative commits

- [`88c1b46`](https://github.com/Dicklesworthstone/frankensim/commit/88c1b46c85f0895f25cc84f8f85a60e1e28c95b5) - implement regime/admission gates.
- [`052d38b`](https://github.com/Dicklesworthstone/frankensim/commit/052d38bc5524e47f8ec931fcb74f9b852bd937e8) - implement material kernels.
- [`a968527`](https://github.com/Dicklesworthstone/frankensim/commit/a96852737b948c61233e6f8c6541233d3c9d10e2) - complete FEEC exterior calculus contract.
- [`89c1f82`](https://github.com/Dicklesworthstone/frankensim/commit/89c1f827e72a6e22a2a19b3d9a9b3737c5bf224c) - add typed operator IR compiler.
- [`de11ad1`](https://github.com/Dicklesworthstone/frankensim/commit/de11ad1caf20be601d165916685903ba45d0dfe8) - lock time-integrator evidence.
- [`7aec6a5`](https://github.com/Dicklesworthstone/frankensim/commit/7aec6a5b1d603d8fe96040dabdb752fb0baac976) - add Krylov and p-multigrid battery.
- [`43d52f2`](https://github.com/Dicklesworthstone/frankensim/commit/43d52f2ed960ede2067eec353359063c86d9dbb3) - add mixed-precision Krylov refinement.

## 6. ASCENT Optimization And Constraint Machinery

The optimization layer moved from problem sketches to typed optimization data,
constraints, geometric constraints, derivative-free engines, Goodhart checks,
and problem serialization.

### Delivered capability

- `fs-opt` typed objective/constraint graphs, manifold variables, PDE and
  stochastic nodes, differentiability routing, canonical serialization, and
  the Goodhart guard.
- `fs-constraint` typed hard/soft/chance/robust/cert/fabrication/code
  constraints, interval certificates, unsat cores, and ranked repairs.
- `fs-dfo` CMA-ES/IGO core with deterministic restart behavior and DFO
  no-claim boundaries.
- `fs-geocon` geometric constraint primitives for min-thickness, draft,
  symmetry, envelopes, and volume/mass.

### Closed workstreams

- [`frankensim-epic-ascent-7tv.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - optimization problem IR.
- [`frankensim-epic-ascent-7tv.2`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - constraint calculus.
- [`frankensim-epic-ascent-7tv.4`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - CMA-ES/DFO engines.
- [`frankensim-epic-epistype-qmao.5`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - Goodhart guard.
- [`frankensim-epic-morph-wqd.21`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - geometric constraint primitives.

### Representative commits

- [`5b86049`](https://github.com/Dicklesworthstone/frankensim/commit/5b860492632414045326fa7a260701eae515cc16) - scaffold `fs-opt`.
- [`64dbee1`](https://github.com/Dicklesworthstone/frankensim/commit/64dbee1cf45d28669e36e9eb658e08dab5827f78) - add the typed constraint calculus.
- [`f536e98`](https://github.com/Dicklesworthstone/frankensim/commit/f536e98f8254a443988830b2365e6a29d865dd3a) - add the Goodhart guard.
- [`b5f1210`](https://github.com/Dicklesworthstone/frankensim/commit/b5f12106da8c28f0a97e12ba2df968790537c8ac) - harden DFO/CMA behavior.
- [`19a4193`](https://github.com/Dicklesworthstone/frankensim/commit/19a4193b5356f8c910f3b3166ae52cdc8a068584) - add geometric constraint primitives.

## 7. HELM, Ledger, IR, Sessions, And Evidence

The orchestration layer gained evidence-carrying values, ledger schema,
time-travel queries, FrankenScript IR, admission checks, sessions/governor
support, planning models, image artifacts, and a narrow PV vertical skeleton.

### Delivered capability

- `fs-evidence` `Evidence<T>`/`Certified<T>` wrappers, model cards,
  discrepancy/bracketing records, falsifier pairing, and evidence gates.
- `fs-ledger` content-addressed artifacts, operations, event streams, tune rows,
  lineage, time travel, physics-VCS primitives, tombstone memory, and telemetry.
- `fs-ir` typed FrankenScript AST, JSON/s-expression syntax, admission checks,
  declarative query IR, and structured errors.
- `fs-plan` cost/error ledgers and planning models.
- `fs-session` capability tokens, governor checks, idempotency, and dry-run
  estimates.
- `fs-vskeleton` a photovoltaic vertical skeleton tying geometry, PDE,
  objective, adjoint, optimization, and ledger concepts together.
- `fs-img` deterministic PNG/EXR-style artifact plumbing and bias-labeled
  denoising.

### Closed workstreams

- [`frankensim-epic-helm-gp3.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - ledger core.
- [`frankensim-epic-helm-gp3.2`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - evidence wrappers.
- [`frankensim-epic-helm-gp3.3`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - ledger time travel.
- [`frankensim-epic-helm-gp3.4`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - IR core.
- [`frankensim-epic-helm-gp3.5`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - IR admission.
- [`frankensim-epic-helm-gp3.7`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - session/governor.
- [`frankensim-epic-foundations-huq.8`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - PV vertical skeleton milestone.

### Representative commits

- [`1e20bda`](https://github.com/Dicklesworthstone/frankensim/commit/1e20bdac6b18431dae1507c874721709cefc1fee) - add Design Ledger v0.
- [`a8f88f8`](https://github.com/Dicklesworthstone/frankensim/commit/a8f88f82c3aaebc0864926e9b6125d7d27f1c14c) - add evidence and certification core.
- [`33c1ed4`](https://github.com/Dicklesworthstone/frankensim/commit/33c1ed4c6ba49af09740e4b3c1cc9ec57449a30d) - complete ledger time travel.
- [`21369ba`](https://github.com/Dicklesworthstone/frankensim/commit/21369ba96892ba136bf8630a221005107f698be9) - add FrankenScript typed AST.
- [`df6e8dc`](https://github.com/Dicklesworthstone/frankensim/commit/df6e8dc93478dcf151de77b773ef1b3fd3429763) - implement IR static admission.
- [`c8e5ca9`](https://github.com/Dicklesworthstone/frankensim/commit/c8e5ca9523f7586d219799cbefbaf2d4e1b9dbfe) - implement session governor.
- [`8f71102`](https://github.com/Dicklesworthstone/frankensim/commit/8f71102ca9f82c68d69848dc5d9d493e164c435d) - complete `fs-img`.

## 8. Addendum Flywheel, Certified Speculation, And Governance

The addendum wave added a machine-readable governance layer and a speculative
execution/evidence flywheel: three-color epistemic schema, falsifier gates,
incremental recompute, physics VCS, semantic diff, proposer zoo, speculation
economics, risk registers, contract registries, and an executable Phase 0 spine
gate.

### Delivered capability

- Three-color ledger schema for verified/validated/estimated evidence.
- Falsifier pairing, consequence-by-doubt budgets, and no-falsifier-no-ship
  checks.
- Slack-bearing recompute store and tolerance-aware invalidation.
- Physics-VCS commit/checkout primitives and semantic diff.
- Certified speculation verifier, proposer zoo, and accept/reject economics.
- `fs-spececo`, `fs-verify`, `fs-probe`, `fs-ladder`, `fs-recompute`,
  `fs-govern`, and `fs-contract` addendum crates.
- Phase 0 spine gate: laundering refusal, falsifier pairing, tombstone protocol,
  interface-type checks, Goodhart guard, and budget-pie reporting.

### Closed workstreams

- [`frankensim-epic-epistype-qmao.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - three-color schema.
- [`frankensim-epic-epistype-qmao.4`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - falsifier pairing.
- [`frankensim-epic-flywheel-lmp4.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - certified speculation verifier.
- [`frankensim-epic-flywheel-lmp4.2`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - proposer zoo.
- [`frankensim-epic-flywheel-lmp4.3`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - accept/reject economics.
- [`frankensim-epic-flywheel-lmp4.6`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - recompute store.
- [`frankensim-epic-flywheel-lmp4.10`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - semantic diff.
- [`frankensim-epic-addendum-xpck.2`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - Phase 0 spine complete.
- [`frankensim-epic-flywheel-lmp4.14`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - assume-guarantee contracts.

### Representative commits

- [`e43e3b1`](https://github.com/Dicklesworthstone/frankensim/commit/e43e3b14c8db88bbc5263a98cfc9cfde19fcc2e3) - add the three-color epistemic schema and ledger gate.
- [`39fd1a5`](https://github.com/Dicklesworthstone/frankensim/commit/39fd1a5144ae211073fbf8100774f1b4366304b3) - add falsifier pairing.
- [`ea102b5`](https://github.com/Dicklesworthstone/frankensim/commit/ea102b596b725fe77808ee6a6c16473a8d2e5f11) - add the slack-bearing recompute store.
- [`772d975`](https://github.com/Dicklesworthstone/frankensim/commit/772d97508f1e19bf5d551c1147bd83b65f4e470f) - add physics-VCS base verbs.
- [`2f2fe56`](https://github.com/Dicklesworthstone/frankensim/commit/2f2fe56339aa811649d1ebd0702d0f1cf43e09b3) - add certified-speculation verifier v0.
- [`6c6f633`](https://github.com/Dicklesworthstone/frankensim/commit/6c6f6338e04fd412ba29c630f4d4bfd86a56d1d3) - add proposer zoo.
- [`92636d0`](https://github.com/Dicklesworthstone/frankensim/commit/92636d0bb5071a503ed872bb8337cd90b2676350) - land integrated speculation economics and ledger v3 telemetry.
- [`544eaee`](https://github.com/Dicklesworthstone/frankensim/commit/544eaee3c0fd27ab672c9f6cf36fbcbd7a99204a) - close the Phase 0 spine milestone.
- [`b33496e`](https://github.com/Dicklesworthstone/frankensim/commit/b33496e291b928f5113bd65c651b2d5b3b97b044) - add assume-guarantee component contracts.

## 9. Certified Geometry Conversion, CutFEM, Solids, And Topology Optimization

The post-baseline wave turned geometry and structural mechanics into a more
direct design loop: exact and certified representation conversion, CutFEM on
SDFs, elasticity and contact on cut geometry, structural stability, and both
density and level-set topology optimization landed as working crate surfaces.

### Delivered capability

- `fs-rep-nurbs` gained a certified NURBS-to-SDF converter and an SDF-to-NURBS
  refit path, giving the workspace a round-trip bridge between spline CAD,
  fields, F-reps, and sheaf watertightness checks.
- `fs-cutfem` added certified cut classification, depth-controlled cut
  quadrature, ghost penalties, Nitsche embedded boundary handling, quadtree
  adaptivity, and zero-meshing FEM on SDFs.
- `fs-solid` landed linear elasticity, finite-strain hyperelasticity through
  material cards, B-bar near-incompressibility mitigation, CutFEM frontends,
  stability continuation, buckling pencils, Koiter scaffolding, Cosserat rods,
  fiber sections, force-based beams, and SDF-native contact.
- `fs-topopt` added density topology optimization with robust three-field
  filtering, stress aggregation, eigenfrequency objectives, and golden
  batteries.
- `fs-topols` added level-set topology optimization over the CutFEM stack with
  WENO advection, fast-marching redistancing, velocity extension, topological
  derivatives, and deterministic snapshots.

### Closed workstreams

- [`frankensim-epic-morph-wqd.11`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - certified NURBS-to-SDF conversion.
- [`frankensim-epic-morph-wqd.12`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - SDF-to-NURBS refit.
- [`frankensim-epic-flux-tfz.13`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - `fs-solid` elasticity core.
- [`frankensim-epic-flux-tfz.14`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - advanced structural element kernels.
- [`frankensim-epic-flux-tfz.15`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - stability, continuation, and Koiter path.
- [`frankensim-epic-flux-tfz.16`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - SDF-native IPC contact.

### Representative commits

- [`279d611`](https://github.com/Dicklesworthstone/frankensim/commit/279d611817f854930a8dc20cd7605702a0f14c2f) - add the certified NURBS-to-SDF converter.
- [`c9e1a4b`](https://github.com/Dicklesworthstone/frankensim/commit/c9e1a4b984773dfed3fcd953fd251394567867ec) - close the certified SDF-to-NURBS round trip.
- [`f781085`](https://github.com/Dicklesworthstone/frankensim/commit/f781085a538063c1d6711777030efa3ca5817d19) - ship CutFEM on SDFs.
- [`15eb757`](https://github.com/Dicklesworthstone/frankensim/commit/15eb75700d69ca73d6fec1cc72b18f3f65666a01) - add the `fs-solid` elasticity core.
- [`e396160`](https://github.com/Dicklesworthstone/frankensim/commit/e3961608169ab16ea19beb69b368857ba5636221) - add stability, continuation, and Koiter scaffolding.
- [`bcc71d4`](https://github.com/Dicklesworthstone/frankensim/commit/bcc71d484b001d7188b920308971b7800c4cf543) - add structural rods, fiber sections, and force-based beams.
- [`4b4c24e`](https://github.com/Dicklesworthstone/frankensim/commit/4b4c24ef48d93be96cb95c4b9b362bb414d6c8d7) - add the `fs-topols` level-set topology optimization crate.
- [`5015e57`](https://github.com/Dicklesworthstone/frankensim/commit/5015e5795a805e37fda6f58b3f46e736c0428557) - add stress aggregate objective plumbing.
- [`0c66ce7`](https://github.com/Dicklesworthstone/frankensim/commit/0c66ce7741e3e48f9289f4a8e975d2ad12e74084) - add eigenfrequency objective support.
- [`7616efe`](https://github.com/Dicklesworthstone/frankensim/commit/7616efedd5e76fb5394a284dd7dab4bc867462b1) - close the SDF-native contact lane.

## 10. Differentiation, Adaptivity, Planning, And Self-Knowledge

The second addendum wave made derivatives, goal-oriented accuracy, and economic
planning first-class rather than incidental: gradients now carry evidence,
planners refuse when budgets are inadequate, and value-of-information work
turns ignorance into a priced queue.

### Delivered capability

- `fs-adjoint` gained a ledger-DAG VJP registry, fail-loud missing-gradient
  behavior, differentiability mitigations, interval residual gradient
  certificates, explanation objects, DWR accept gates, and exact Hessian-vector
  products for density misfit.
- `fs-dwr` added goal-oriented error estimation, anisotropic metric synthesis,
  marking, tile-level adaptivity scaffolding, and contract-level accuracy
  boundaries.
- `fs-ir` added a greedy fidelity-ladder planner, anytime query semantics,
  colored intervals, priced tightening hints, and explicit refusal reports.
- `fs-plan` added value-of-information query planning that ranks surrogate rung
  climbs and physical validation probes by decision impact per dollar.
- `fs-flywheel-e2e` added whole-loop, phase-1, and phase-2 gates that measure
  the compound effect of speculation, recompute, sheaf merges, tombstones,
  adjoint gradients, and planner ladders.

### Closed workstreams

- [`frankensim-epic-coupling-bk0o.1`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - ledger-DAG transposition.
- [`frankensim-epic-coupling-bk0o.3`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - gradient certificates and merge gates.
- [`frankensim-epic-flywheel-lmp4.16`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - fidelity-ladder planner.
- [`frankensim-epic-flywheel-lmp4.17`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - anytime and refusal semantics.
- [`frankensim-epic-flywheel-lmp4.18`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - whole-loop flywheel proof harness.
- [`frankensim-epic-addendum-xpck.4`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - Phase 2 leverage gate.

### Representative commits

- [`41cadcc`](https://github.com/Dicklesworthstone/frankensim/commit/41cadccb11349a1ab5b48374b5e5f875906cdfd7) - add ledger-DAG transposition.
- [`3fab970`](https://github.com/Dicklesworthstone/frankensim/commit/3fab970ad002fb16223020cd8dfde362c117e004) - add gradient certificates.
- [`e907bc6`](https://github.com/Dicklesworthstone/frankensim/commit/e907bc6796b5e8b840b52bc618eae1657a969bf2) - add the DWR goal-oriented accept test.
- [`9ee7227`](https://github.com/Dicklesworthstone/frankensim/commit/9ee7227858c7f09c86bd694b600a1bb263a013b1) - add the greedy fidelity-ladder planner.
- [`1ade44a`](https://github.com/Dicklesworthstone/frankensim/commit/1ade44a99cf9a67a6be8f6fa10b914209f9622c6) - add anytime and refusal semantics.
- [`9ab427e`](https://github.com/Dicklesworthstone/frankensim/commit/9ab427e3ca95b3e0e0556547c84234ce4ae37321) - land the whole-loop flywheel harness.
- [`3f9a714`](https://github.com/Dicklesworthstone/frankensim/commit/3f9a714868c2de4a8ab62359f903a7990c272b4c) - add the Phase 2 leverage gate.
- [`941a67e`](https://github.com/Dicklesworthstone/frankensim/commit/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f) - add value-of-information queries to `fs-plan`.

## 11. Extended Physics, Numerics, And Solver Kernels

The FLUX and BEDROCK layers broadened from core abstractions into usable
simulation kernels: harmonic FEEC analysis, domain decomposition, pressure-
robust flow, FMM/BEM, LBM, vortex particles, time slabs, and large-domain
deterministic trigonometry.

### Delivered capability

- `fs-feec` added cohomology, harmonic cochains, Hodge splitting, circulation
  extraction, and high-order vector-family spaces.
- `fs-dd` added BDDC domain decomposition with sheaf-framed edge coarse spaces
  and measured conditioning tables.
- `fs-flux` added a BDM1-P0 pressure-robust Navier-Stokes scaffold with
  boundary moments, conformance batteries, Picard stepping, and discrete
  adjoint hooks.
- `fs-fmm` and `fs-bem` added black-box Chebyshev FMM, Laplace panel methods,
  Hess-Smith airfoil screening, wake roll-up, and transpose/FMM validation.
- `fs-lbm` added D2Q9 core dynamics, rheology, thermal coupling, free-surface
  VOF, contact-line bracketing, and refinement.
- `fs-vpm`, `fs-iga`, and `fs-couple` added vortex-particle, isogeometric, and
  passive port-Hamiltonian coupling surfaces.
- `fs-time` gained time slabs as ledger cells, while `fs-math` extended
  deterministic trig with Payne-Hanek reduction across all finite `f64`.

### Closed workstreams

- [`frankensim-epic-flux-tfz.7`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - FEEC cohomology and harmonic decomposition.
- [`frankensim-epic-flux-tfz.11`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - domain decomposition.
- [`frankensim-epic-flux-tfz.17`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - pressure-robust Navier-Stokes.
- [`frankensim-epic-flux-tfz.19`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - LBM extensions.
- [`frankensim-r6r5`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - Payne-Hanek trig reduction.

### Representative commits

- [`7d95552`](https://github.com/Dicklesworthstone/frankensim/commit/7d95552cf6bbb1055296991fb43f119b27916635) - add FEEC cohomology.
- [`57d775e`](https://github.com/Dicklesworthstone/frankensim/commit/57d775e47f3d38f94b7c1d017c3c74888e551ec0) - add BDDC domain decomposition.
- [`8f031ad`](https://github.com/Dicklesworthstone/frankensim/commit/8f031ad48a56a3ac5e91023cd99782138918bcd4) - complete pressure-robust Navier-Stokes.
- [`8c5b0e5`](https://github.com/Dicklesworthstone/frankensim/commit/8c5b0e555e9d3d0f1a7840e941fed43890519008) - add `fs-fmm` and `fs-bem`.
- [`eeb919a`](https://github.com/Dicklesworthstone/frankensim/commit/eeb919a98ddb7dc8acb8aff87fb6fa85138eeb3b) - implement the FMM transpose panel operator.
- [`714e9e6`](https://github.com/Dicklesworthstone/frankensim/commit/714e9e6ef39d590ddfd4fa42254fec89a136cc4b) - add the D2Q9 LBM core.
- [`ded5b78`](https://github.com/Dicklesworthstone/frankensim/commit/ded5b78909a40b81b493289bea2cd6742aa99eb1) - complete LBM rheology, thermal, free-surface, and refinement extensions.
- [`43ea2dd`](https://github.com/Dicklesworthstone/frankensim/commit/43ea2dd57b29c33e221250cec35b5d8db814a6d5) - add the vortex-particle method core.
- [`8a90167`](https://github.com/Dicklesworthstone/frankensim/commit/8a9016744727cdeb0b9788df6999e86e9bdc38ab) - add time slabs as cells.
- [`7049ca3`](https://github.com/Dicklesworthstone/frankensim/commit/7049ca38995f62664d0715580ae7b8fe987f22ae) - ship Payne-Hanek trig reduction for all finite `f64`.

## 12. Optimization, UQ, Robust Design, And Archives

The ASCENT layer became much broader than the initial DFO core. It now includes
gradient optimization, Bayesian optimization, uncertainty quantification,
surrogates, proof-carrying optimization, value-of-information, fabrication
constraints, quality-diversity archives, robust design, and multiple explicit
design-space representations.

### Delivered capability

- `fs-ascent` added L-BFGS, trust-region Newton, Riemannian routing, augmented
  Lagrangians, Wolfe search, stopping policy, and gradient Pareto tracing.
- `fs-bo` added Gaussian-process Bayesian optimization, TuRBO trust regions,
  and multi-fidelity BO.
- `fs-dfo` expanded with log-domain Sinkhorn optimal transport, deterministic
  multi-objective utilities, and a finite-support Wasserstein DRO oracle.
- `fs-uq`, `fs-surrogate`, `fs-voi`, `fs-race`, and `fs-sos` added UQ,
  certify-or-escalate ROMs, active validation, anytime-valid e-racing, and SOS
  certificates.
- `fs-archive` added MAP-Elites/CVT quality-diversity archives, while
  `fs-lattice`, `fs-truss`, and `fs-fab` added homogenized lattice design,
  ground-structure sizing, and manufacturing/code-compliance constraints.
- `fs-shapeprog`, `fs-rep-neural`, `fs-toleralloc`, `fs-assimilate`, and
  `fs-robust` added program synthesis, Lipschitz-certified neural implicits,
  tolerance allocation, as-built/data assimilation, and objective epistemics.

### Closed workstreams

- [`frankensim-epic-ascent-7tv.18`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - value-of-information and active validation.
- [`frankensim-epic-ascent-7tv.20`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - fabrication and code-compliance constraints.
- [`frankensim-epic-ascent-7tv.7`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - e-racing with anytime-valid cancellation.
- [`frankensim-epic-ascent-7tv.10`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - proof-carrying SOS optimization.
- [`frankensim-epic-ascent-7tv.14`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - lattice homogenization and graded optimization.
- [`frankensim-vcia`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - Sinkhorn OT and DRO evidence lane.

### Representative commits

- [`f1d2028`](https://github.com/Dicklesworthstone/frankensim/commit/f1d202841d6cb4c29fb60a1fdbfd8624df1c8c69) - add the gradient optimizer stack.
- [`b93bd21`](https://github.com/Dicklesworthstone/frankensim/commit/b93bd21afce726eb3fddbd3b8c20db3d428a796d) - add gradient Pareto tracing.
- [`a76cbea`](https://github.com/Dicklesworthstone/frankensim/commit/a76cbea571cb0cb9b324a4734a200e40a63a8ca2) - add the value-of-information crate.
- [`dcfc2ac`](https://github.com/Dicklesworthstone/frankensim/commit/dcfc2ac3503004c5cf227e05282d47fdbb659ac2) - add TuRBO trust-region BO.
- [`2d3a3c4`](https://github.com/Dicklesworthstone/frankensim/commit/2d3a3c41de507156d9a22eac5385f649dc217f2c) - add multi-fidelity BO.
- [`961e2d6`](https://github.com/Dicklesworthstone/frankensim/commit/961e2d60cad948a0b48ac52256dbc89ee154bf8c) - add log-domain Sinkhorn optimal transport.
- [`382e616`](https://github.com/Dicklesworthstone/frankensim/commit/382e6169cfcfb951006ba53dfe9cd6b175454ee1) - add deterministic multi-objective utilities.
- [`d947001`](https://github.com/Dicklesworthstone/frankensim/commit/d9470015316f0e587cc797ff53b4a0dbf359d7fc) - add the discrete Wasserstein DRO oracle.
- [`c0f1a5c`](https://github.com/Dicklesworthstone/frankensim/commit/c0f1a5c2ef63582e173efc87be255e3af1242e6c) - add SOS proof-carrying optimization.
- [`3c56418`](https://github.com/Dicklesworthstone/frankensim/commit/3c56418c549420e1f62886af9cd1c31a78c38a38) - add MAP-Elites/CVT quality-diversity archives.
- [`4c3ec84`](https://github.com/Dicklesworthstone/frankensim/commit/4c3ec84cb5cdf5dcbf53e301accb3fc9e908be70) - ship lattice homogenization and graded optimization.
- [`f9fad53`](https://github.com/Dicklesworthstone/frankensim/commit/f9fad5309f23a18de2d9356cecb01d21938e8aa7) - add certified ground-structure truss sizing.

## 13. Rendering, Reporting, Browser Surface, And End-To-End Campaigns

The project now has several capstone paths where many crates compose into a
visible workflow: scientific rendering, browser demos, lab notebooks, flagship
optimization studies, and cross-layer end-to-end certification harnesses.

### Delivered capability

- `fs-render` added unbiased spectral path tracing, certified chart tracing,
  mixed-scene backends, and Woodcock volume rendering over live LBM fields.
- `fs-img`, `fs-viz`, and `fs-report` hardened PNG validation, added analytic
  scientific visualization primitives, and generated automatic lab notebooks
  with semantic diffs.
- `fs-wasm` added browser-oriented numerical kernel surfaces and tiered demo
  modules.
- `fs-marquee` added the raw-SDF, zero-meshing marquee study runner.
- `fs-frame` added the seismic-minimal frame flagship with truss layout,
  fiber-hinge time history, anytime-valid fragility, and CVaR sizing.
- `fs-thrust-e2e` added a certified quality-diversity vortex-thruster campaign
  composing VPM, evidence, surrogates, archives, and reports.
- [`fb08842`](https://github.com/Dicklesworthstone/frankensim/commit/fb088428ae810ee3ddd893712588c1e64f6cc0c4)
  added three additional certified end-to-end campaigns:
  proof-robust optimization, campaign scheduling, and flutter certification.

### Closed workstreams

- [`frankensim-epic-lumen-qfx.1`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - spectral path tracing.
- [`frankensim-epic-lumen-qfx.2`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - chart render backends.
- [`frankensim-epic-lumen-qfx.3`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - volume rendering.
- [`frankensim-epic-flagships-mye.1`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - marquee raw-SDF optimization study.
- [`frankensim-epic-flagships-mye.3`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - seismic-minimal frame flagship.

### Representative commits

- [`b089421`](https://github.com/Dicklesworthstone/frankensim/commit/b089421c5494e231d3623973251e91a30169efb0) - add the spectral path-tracing core.
- [`f8a927e`](https://github.com/Dicklesworthstone/frankensim/commit/f8a927e601372affc84484e041f8f99d513fe88a) - add certified chart render backends.
- [`0ccee78`](https://github.com/Dicklesworthstone/frankensim/commit/0ccee78a87ec4e81ca38554142e7506aea356bf5) - ship Woodcock volume rendering over live LBM fields.
- [`8aaecbc`](https://github.com/Dicklesworthstone/frankensim/commit/8aaecbcd1881349484bdabd3f71a5f38d250ee28) - add analytically verified scientific visualization primitives.
- [`6e5ae66`](https://github.com/Dicklesworthstone/frankensim/commit/6e5ae66563ec1267711c591ef23bbf47a4568596) - add automatic lab notebooks and semantic diffs.
- [`68768ba`](https://github.com/Dicklesworthstone/frankensim/commit/68768ba496da9ab6698041429f3b9c23a5d291f6) - add the standalone browser numerical-kernel surface.
- [`6df2c03`](https://github.com/Dicklesworthstone/frankensim/commit/6df2c03eb3a2ea5c8c731faad7cabda19475aabb) - ship the raw-SDF marquee study runner.
- [`8d27622`](https://github.com/Dicklesworthstone/frankensim/commit/8d27622657d0f4ace6bdc40991f295e9b958e34c) - ship the seismic-minimal frame flagship.
- [`dc1bf7f`](https://github.com/Dicklesworthstone/frankensim/commit/dc1bf7faf7ed6feafcda0f7c3c3eb78fe6d1d0e0) - add the certified quality-diversity vortex-thruster campaign.
- [`fb08842`](https://github.com/Dicklesworthstone/frankensim/commit/fb088428ae810ee3ddd893712588c1e64f6cc0c4) - adds proof-robust, schedule, and flutter certified end-to-end campaigns.

## 14. Latest Campaigns, AD Bridges, And p-Multigrid Smoothing

The latest 2026-07-08 wave adds more proof-carrying end-to-end examples and
tightens two core numerical surfaces: inverse-trig AD support and p-independent
multigrid smoothing.

### Delivered capability

- `fs-neuroshape-e2e` and `fs-grammar-e2e` add certified neural-shape topology
  and geometric-grammar campaigns, composing neural implicit charts,
  visualization/evidence, shape programs, archives, and fabrication checks.
- `fs-oed-e2e` adds SensorForge, a value-of-information sensor-placement
  campaign that stops when the posterior design choice is robust.
- `fs-vessel` adds the laminar-pour vessel flagship, combining Chebyshev
  stability objectives, free-surface LBM, CVaR robustification, e-racing, and
  spectral volume rendering.
- `fs-metamat-e2e`, `fs-truss-e2e`, `fs-adaptbo-e2e`, and `fs-flowcert-e2e`
  add certified metamaterial frontiers, optimal truss load paths, anytime BO,
  and CFD credibility maps.
- `fs-math` adds deterministic `asin`/`acos`, and `fs-ad` extends the `Real`
  surface with inverse-trig operations and chain-rule tests.
- `fs-ad` also gains the FrankenTorch tape bridge, binomial Revolve schedule,
  snapshot-store seam, and matrix-free IFT tangent route.
- `fs-solver` replaces p-MG Jacobi scaling with PU-symmetrized
  vertex-centered additive Schwarz smoothing and an exact r=1 Pavarino coarse
  term.

### Closed workstreams

- [`frankensim-epic-flagships-mye.4`](https://github.com/Dicklesworthstone/frankensim/blob/319cb64f052d76e15882ee53ace41092881c7fa8/.beads/issues.jsonl) - laminar-pour vessel flagship.
- [`frankensim-t88x`](https://github.com/Dicklesworthstone/frankensim/blob/319cb64f052d76e15882ee53ace41092881c7fa8/.beads/issues.jsonl) - inverse-trig `Real` operations for AD consumers.
- [`frankensim-o3ui`](https://github.com/Dicklesworthstone/frankensim/blob/319cb64f052d76e15882ee53ace41092881c7fa8/.beads/issues.jsonl) - `fs-ad` bridge, Revolve, spill, and IFT integrations.
- [`frankensim-x08j`](https://github.com/Dicklesworthstone/frankensim/blob/319cb64f052d76e15882ee53ace41092881c7fa8/.beads/issues.jsonl) - p-independent p-MG smoothing.

### Representative commits

- [`5cbdd90`](https://github.com/Dicklesworthstone/frankensim/commit/5cbdd903e8caf389e98eafe5ebc2436ac7a0093c) - add neural-shape topology and geometric grammar e2e campaigns.
- [`94404c4`](https://github.com/Dicklesworthstone/frankensim/commit/94404c4d8ef5f4fb7a33d9e46fddb32f51cb0aaa) - add SensorForge OED e2e.
- [`b95e00f`](https://github.com/Dicklesworthstone/frankensim/commit/b95e00ff7ee8140d0c1427824a23953af1353963) - ship the laminar-pour vessel flagship.
- [`8028bbc`](https://github.com/Dicklesworthstone/frankensim/commit/8028bbc29de12aff9fcb46849e688d0cd28bbab5) - add the metamaterial stiffness-density frontier campaign.
- [`4c573d6`](https://github.com/Dicklesworthstone/frankensim/commit/4c573d648bafd6d779def16e4ca7303fdb99f870) - add the optimal truss critical-load-path campaign.
- [`76e9c89`](https://github.com/Dicklesworthstone/frankensim/commit/76e9c89f62e91612b8886115c06ba8260c76cc94) - add AnytimeBO.
- [`922c835`](https://github.com/Dicklesworthstone/frankensim/commit/922c83597d2a79c1a862d9b8e37f5240b8ec9c8e) - add deterministic inverse trig and AD `Real` operations.
- [`3eb480c`](https://github.com/Dicklesworthstone/frankensim/commit/3eb480c8f4c1faa0d2b92a3b31e0d9db46256f0e) - add the FlowCert CFD credibility map.
- [`7575cdd`](https://github.com/Dicklesworthstone/frankensim/commit/7575cdd6bb0728794ef0fce4eae3512f9667b93d) - ship `fs-ad` adjoint integrations.
- [`319cb64`](https://github.com/Dicklesworthstone/frankensim/commit/319cb64f052d76e15882ee53ace41092881c7fa8) - add vertex-patch Schwarz smoothing to p-MG.

## 15. 2026-07-09 Mainline Expansion: Flagships, Decision Loops, And Proof Hygiene

The window from
[`319cb64`](https://github.com/Dicklesworthstone/frankensim/commit/319cb64f052d76e15882ee53ace41092881c7fa8)
through
[`e08e302`](https://github.com/Dicklesworthstone/frankensim/commit/e08e30280bcd7af05ae55e990b129d6f75192ead)
is another implementation-heavy mainline slice, not a formal release. It adds
63 commits across 116 files, with 16,163 insertions and 369 deletions. The work
broadens flagship coverage, strengthens decision-making loops for optimization
and uncertainty, and fixes campaign checks that had been too correlated or too
underspecified to support the claims around them.

### Delivered capability

- `fs-ornith` lands as a smoke-tier ornithoid aircraft flagship with parameter
  decoding, BEM/VPM screening, e-raced candidate elimination, LBM refinement,
  SOS-style trim certificates, conformal lift/drag bands, Pareto-atlas rows,
  replay evidence, and a checked-in contract.
- `fs-topopt` gains a feature-gated `cutfem-marquee` lane for topology
  optimization over density-as-SDF fields, DWR-guided octree refinement,
  zero-remesh iteration logs, and topology/thickness witnesses.
- `fs-solver` adds Stokes-class saddle-point machinery: Schur approximations,
  block preconditioners, and PMINRES coverage over the FEEC mixed-space surface.
- `fs-dfo` expands from many-objective hypervolume accounting into NSGA-III,
  MOEA/D, contribution archiving, and world-fork steering. `fs-uq` adds chance
  constraints, seismic uncertainty helpers, anytime stopping, and adaptive MLMC
  admission checks.
- `fs-bo` adds inducing-point sparse Gaussian-process support and then hardens
  the sparse/exact comparison so approximate covariance loss is visible instead
  of smoothed over.
- `fs-ascent` adds interior-point and SQP constrained polish engines, a
  FrankenScipy oracle lane, Rosenbrock basin checks, and a Problem-IR study
  runner.
- `fs-render` adds edge-aware differentiable rendering and an inverse-rendering
  objective surface, while `fs-wasm` exposes deeper upper-stack kernels and a
  certified campaign browser tier.
- `fs-gen` lands as a proposal-only generation crate with compiler-enforced
  epistemic boundaries, making generated candidates explicit about what they
  are not allowed to certify.
- `fs-flagship-e2e` lands as an L6 replay suite for the flagship family:
  smoke/mid/full stage wiring, frozen smoke-stage content hashes,
  cross-flagship LBM and e-race audits, structured failure drills, forensic
  JSON rows, and ignored mid/full fidelity lanes for later perf cadence.
- `fs-mesh` hardens full-Ruppert refinement by checking candidate Steiner
  points against hull-facet diametral balls before insertion, while `fs-topopt`
  tightens the marquee tests around faster deterministic topology, benchmark,
  and thickness witnesses.
- `fs-simd` adds a `frontier-sme2` exploratory capsule for a runtime-gated
  streaming-mode 16x16 f32 GEMM tile, registers the unsafe boundary, documents
  the SME2 safety invariants, and adds a hardware-optional battery that keeps
  unsupported machines inert.
- End-to-end proof hygiene improves across `fs-flutter-e2e`,
  `fs-neuroshape-e2e`, `fs-flowcert-e2e`, `fs-truss-e2e`, and
  `fs-metamat-e2e`, including independent flutter stability checks, closed
  neural-shape boundary frames, steady-state credibility maps, corrected
  duality-gap enclosure direction, and PSD certificate symmetrization.
- HELM and flywheel surfaces gain a budget allocator, moonshot allocation
  tests, self-knowledge e2e battery, epistemic-engine acceptance gate, and
  Phase 3 horizon gate.

### Closed and active workstreams

- [`frankensim-epic-flagships-mye.2`](https://github.com/Dicklesworthstone/frankensim/blob/d5873bfd82a3c2dbe359c11aef4947a5def8cdba/.beads/issues.jsonl) - ornithoid multi-inlet aircraft flagship.
- [`frankensim-b7d0`](https://github.com/Dicklesworthstone/frankensim/blob/d5873bfd82a3c2dbe359c11aef4947a5def8cdba/.beads/issues.jsonl) - active CutFEM-octree topology-optimization marquee lane.
- [`frankensim-avuw`](https://github.com/Dicklesworthstone/frankensim/blob/d5873bfd82a3c2dbe359c11aef4947a5def8cdba/.beads/issues.jsonl) - Stokes-class saddle-point block preconditioners.
- [`frankensim-vcia`](https://github.com/Dicklesworthstone/frankensim/blob/d5873bfd82a3c2dbe359c11aef4947a5def8cdba/.beads/issues.jsonl) - many-objective, steering, and Wasserstein/DRO decision surfaces.
- [`frankensim-ijil`](https://github.com/Dicklesworthstone/frankensim/blob/d5873bfd82a3c2dbe359c11aef4947a5def8cdba/.beads/issues.jsonl) - constrained polish engines, FrankenScipy oracle, and Problem-IR runner.
- [`frankensim-epic-helm-gp3.9`](https://github.com/Dicklesworthstone/frankensim/blob/d5873bfd82a3c2dbe359c11aef4947a5def8cdba/.beads/issues.jsonl) - `fs-plan` allocation and replanning machinery.
- [`frankensim-epic-lumen-qfx.5`](https://github.com/Dicklesworthstone/frankensim/blob/d5873bfd82a3c2dbe359c11aef4947a5def8cdba/.beads/issues.jsonl) - differentiable rendering objective lane.
- [`frankensim-epic-flagships-mye.5`](https://github.com/Dicklesworthstone/frankensim/blob/20e9825018aeb30b5b2b12c122c7c95c84b3a646/.beads/issues.jsonl) - flagship e2e suite with frozen smoke-stage goldens.
- [`frankensim-uee3`](https://github.com/Dicklesworthstone/frankensim/blob/1fe4ef5c21737a92aa58dd02f901052de5b0a88b/.beads/issues.jsonl) - active fs-mesh v2 follow-up.

### Representative commits

- [`7624964`](https://github.com/Dicklesworthstone/frankensim/commit/7624964781185c5d82879ba0e63aba9da304021d) - add the self-knowledge e2e battery.
- [`3dcbfaa`](https://github.com/Dicklesworthstone/frankensim/commit/3dcbfaaf3d3e90496719bb6a879500ca825e61b7) - add Stokes saddle-point block preconditioners and PMINRES.
- [`8c20bab`](https://github.com/Dicklesworthstone/frankensim/commit/8c20bab6da5ef694d4243ded896ac540d87502c8) - add the NSGA-III reference-direction lane.
- [`61e9348`](https://github.com/Dicklesworthstone/frankensim/commit/61e93489985b6cdbcea768535e6a12d12885f471) - add edge-aware differentiable rendering.
- [`fe8a2c2`](https://github.com/Dicklesworthstone/frankensim/commit/fe8a2c2061377e47293a3c638eb152337a97adcd) - add the `fs-plan` budget allocator.
- [`6d074b8`](https://github.com/Dicklesworthstone/frankensim/commit/6d074b8c92c396ac320815074640c63410f6df95) - add seismic uncertainty, anytime stopping, and adaptive MLMC.
- [`e2b3dff`](https://github.com/Dicklesworthstone/frankensim/commit/e2b3dffb1b4b58b7378e337e3134f4e4d95fb7ac) - add the sparse Gaussian-process lane in `fs-bo`.
- [`73c1c6b`](https://github.com/Dicklesworthstone/frankensim/commit/73c1c6b46f0f36673cc2bbda70164c0b11eee23c) - add `fs-gen` with proposal-only epistemics.
- [`51adc98`](https://github.com/Dicklesworthstone/frankensim/commit/51adc987414491b75ad4543fae7925e83804471a) - add steering and chance-constraint integration.
- [`2f29772`](https://github.com/Dicklesworthstone/frankensim/commit/2f29772b8536171929f2b4bae4ece8aa849caa3c) - add interior-point, SQP, oracle, and IR runner machinery.
- [`2199248`](https://github.com/Dicklesworthstone/frankensim/commit/219924842e0e0ebc71ec975b7d17bddf193db403) - expose the certified campaign tier in `fs-wasm`.
- [`ffcb4a9`](https://github.com/Dicklesworthstone/frankensim/commit/ffcb4a91e047b5f660d28083194da6daeb47b190) - add the smoke-tier `fs-ornith` crate.
- [`8c06c93`](https://github.com/Dicklesworthstone/frankensim/commit/8c06c93e7404b7bbd528eebc14d4af43914a8906) - add the feature-gated CutFEM-octree topology lane.
- [`d5873bf`](https://github.com/Dicklesworthstone/frankensim/commit/d5873bfd82a3c2dbe359c11aef4947a5def8cdba) - add the `fs-ornith` contract and close the ornithoid flagship bead.
- [`fbba704`](https://github.com/Dicklesworthstone/frankensim/commit/fbba7049d3f92714185f26ed1473503b78903281) - refine flutter and neuroshape certificate payloads.
- [`1092e94`](https://github.com/Dicklesworthstone/frankensim/commit/1092e946be32b211bf20f4748be2f44683feaeaf) - harden CutFEM-octree marquee boundary evidence.
- [`9dc7417`](https://github.com/Dicklesworthstone/frankensim/commit/9dc7417090e62e7ab2756719034e2ab0d3b32875) - add the staged flagship e2e replay-suite scaffold.
- [`1337058`](https://github.com/Dicklesworthstone/frankensim/commit/1337058ac6bc4b8f1dd9b7cdeb85db2fd0713bfd) - tighten the flagship replay-suite battery lint posture.
- [`1fe4ef5`](https://github.com/Dicklesworthstone/frankensim/commit/1fe4ef5c21737a92aa58dd02f901052de5b0a88b) - mark the fs-mesh v2 follow-up in progress.
- [`438128d`](https://github.com/Dicklesworthstone/frankensim/commit/438128d988082f2183f406d75883b766fd6b7324) - protect hull refinement and tighten marquee witnesses.
- [`70e1d24`](https://github.com/Dicklesworthstone/frankensim/commit/70e1d2437087b8b3dd04237cd56a8daefedbac2b) - refine boundary splitting and grid traversal.
- [`a8e98b2`](https://github.com/Dicklesworthstone/frankensim/commit/a8e98b21b8eb793069be9153d366618af172b801) - ledger the hull encroachment regression bound.
- [`d93ca59`](https://github.com/Dicklesworthstone/frankensim/commit/d93ca592bb6f408095eb1d72ecb25e2a09dd7674) - add the gated SME2 exploratory capsule.
- [`e08e302`](https://github.com/Dicklesworthstone/frankensim/commit/e08e30280bcd7af05ae55e990b129d6f75192ead) - add the SME2 exploratory battery.

## 16. Post-SME2 Follow-Up: Browser Flagships, Mesh Decisions, And Perf Lanes

The window from
[`e08e302`](https://github.com/Dicklesworthstone/frankensim/commit/e08e30280bcd7af05ae55e990b129d6f75192ead)
through
[`6725739`](https://github.com/Dicklesworthstone/frankensim/commit/6725739f42b878b310e7e7e318fb8c980cef71f3)
adds another 42 committed changes across 77 files, with 6,195 insertions and
219 deletions. This is still a moving `main` snapshot rather than a package
release, but it closes several proof lanes that were active in the prior
changelog endpoint and opens clearer continuation surfaces for performance
work.

### Delivered capability

- `fs-wasm` now exposes live flagship browser pipelines rather than only
  lower-level browser kernels. The new `flagships` module exports reduced but
  real `run_ornithoid`, `run_vessel`, and `run_frame` paths that compose the
  corresponding aircraft, laminar-pour vessel, and seismic-frame stacks into
  browser-facing headline reports.
- `fs-mesh` closes the v2 constrained-recovery lane and the v3 interior-facet
  recovery decision. The stack now includes conforming PLC segment recovery,
  exact-predicate ear clipping for simple non-convex planar facets, a
  measured 10-million-point performance lane, and a recorded boundary-layer
  decision: hull-edge protection was measured counterproductive, while
  weighted exudation remains the honest continuation for residual near-coplanar
  slivers.
- `fs-mesh` also adds a feature-gated hex-dominant meshing surface with SH9
  frame fields, MBO smoothing, and a polycube fallback, keeping the frontier
  capability behind `frontier-hexmesh` while preserving the default test lane.
- `fs-sparse` and `fs-substrate` gain a compact CSR performance path, sharded
  deterministic SpMV, tiled deterministic parallel assembly, real shared STREAM
  bandwidth slices, and NUMA first-touch helpers. The work records all-core
  roofline progress while keeping the broader `fs-sparse-perf` lane open.
- `fs-rand` adds fast-path deterministic normal generation through a ziggurat
  sampler, bitwise-equivalent bulk Philox fills, and a dev-only statistical
  stream battery. The strict distribution contract remains explicit about which
  pieces are performance paths and which pieces are certification gates.
- `fs-fft` adds complex-to-real inverse support plus 2D/3D N-D pencil
  transforms, expanding the transform surface that downstream physics and
  visualization crates can use without leaving the deterministic math layer.
- IO, risk, and packaging code paths were made more fail-closed: `fs-io`
  rejects fractional and negative PLY face-list values, conformance tests stop
  using direct panic arms, CVaR/objective/probe inputs reject non-finite values,
  `fs-robust` removes unwrap/tail-slice panic surfaces, and `fs-package`
  includes provenance in package roots while removing unreachable string-write
  panics.
- `fs-topopt` records the CutFEM-octree marquee as a zero-remeshing topology
  evolution contract and adds invalid-density-shape rejection, while the SME2
  prototype remains documented as an exploratory, runtime-gated capsule rather
  than a portable default lane.

### Closed and active workstreams

- [`frankensim-uee3`](https://github.com/Dicklesworthstone/frankensim/blob/6725739f42b878b310e7e7e318fb8c980cef71f3/.beads/issues.jsonl) - closed `fs-mesh` v2: constrained recovery, protected Ruppert refinement, sliver exudation, deterministic parallel coloring, and the measured 10-million-point perf lane.
- [`frankensim-iw3l`](https://github.com/Dicklesworthstone/frankensim/blob/6725739f42b878b310e7e7e318fb8c980cef71f3/.beads/issues.jsonl) - closed `fs-mesh` v3: non-convex interior facet recovery plus a measurement-backed boundary-layer quality decision.
- [`frankensim-b7d0`](https://github.com/Dicklesworthstone/frankensim/blob/6725739f42b878b310e7e7e318fb8c980cef71f3/.beads/issues.jsonl) - closed `fs-topopt` CutFEM-octree marquee lane, including the zero-remesh contract and thickness/boundary evidence.
- [`frankensim-wsbf`](https://github.com/Dicklesworthstone/frankensim/blob/6725739f42b878b310e7e7e318fb8c980cef71f3/.beads/issues.jsonl) - in progress `fs-sparse-perf` lane for STREAM-relative SpMV, parallel assembly, SIMD kernels, CCD sharding, and Franken interop.
- [`frankensim-1za9`](https://github.com/Dicklesworthstone/frankensim/blob/6725739f42b878b310e7e7e318fb8c980cef71f3/.beads/issues.jsonl) - open `fs-rand` performance lane for ziggurat normals, bulk generation, and statistical batteries.
- [`frankensim-27d3`](https://github.com/Dicklesworthstone/frankensim/blob/6725739f42b878b310e7e7e318fb8c980cef71f3/.beads/issues.jsonl) - open `fs-fft` performance lane for radix kernels, SIMD lanes, N-D pencils, and roofline gates.

### Representative commits

- [`e8c496e`](https://github.com/Dicklesworthstone/frankensim/commit/e8c496e4ddf2486a6bdcf3cba9c605adc8e4e5ef) - record the CutFEM-octree topology marquee as zero-remeshing topology evolution.
- [`29fb1ea`](https://github.com/Dicklesworthstone/frankensim/commit/29fb1ea4f4fbe84d9181e804edead61f8925e4dc) - add the SME2 streaming-mode GEMM prototype contract surface.
- [`f759916`](https://github.com/Dicklesworthstone/frankensim/commit/f7599163ddc5a7309e3095a651c7b18454541e75) - add deterministic ziggurat normal generation for the fast-mode rand path.
- [`cb6dbde`](https://github.com/Dicklesworthstone/frankensim/commit/cb6dbdec9a8892a1c09a57ad939bf184801aa0a4) - add the dev-only Philox statistical stream battery.
- [`8254164`](https://github.com/Dicklesworthstone/frankensim/commit/8254164a7dea53dc02343de01c9543a23d99dcf5) - add bitwise bulk Philox fills.
- [`648a891`](https://github.com/Dicklesworthstone/frankensim/commit/648a891c0f7554121bc281ff291ca3045f42a9f1) - add complex-to-real inverse FFT and 2D/3D pencil transforms.
- [`9a6d303`](https://github.com/Dicklesworthstone/frankensim/commit/9a6d303960d638e1f91e5b243bb29f679cc6eb53) - add conforming interior facet recovery to `fs-mesh`.
- [`6510b02`](https://github.com/Dicklesworthstone/frankensim/commit/6510b0266933ae410926b2d572fc7b32efa47fb0) - triangulate simple non-convex recovery facets.
- [`fa3d93d`](https://github.com/Dicklesworthstone/frankensim/commit/fa3d93d2e75acd6aa91216802da489bda3faf50c) - record the 10-million-point `fs-mesh` perf lane closeout.
- [`151a643`](https://github.com/Dicklesworthstone/frankensim/commit/151a6431f338b8e045f1ea7dd17b6f2664bf766e) - add the hex-dominant meshing module.
- [`0cd9de2`](https://github.com/Dicklesworthstone/frankensim/commit/0cd9de26e383f8a804e32287994d8f4b22fbcd97) - add compact/sharded SpMV, deterministic parallel assembly, and the sparse roofline lane.
- [`dccc078`](https://github.com/Dicklesworthstone/frankensim/commit/dccc0784e9250557328671c584e17df8b5d7e6a1) - add NUMA first-touch to the STREAM triad baseline and compact CSR path.
- [`b4f849f`](https://github.com/Dicklesworthstone/frankensim/commit/b4f849fdfb9ef12373a4720b84274eaf7a46ac75) - reject invalid density lattice shapes in `fs-topopt`.
- [`220d28c`](https://github.com/Dicklesworthstone/frankensim/commit/220d28c0d64f61a7607ad17d751a46a85039c0a9) - reject fractional PLY face-list indices.
- [`3842369`](https://github.com/Dicklesworthstone/frankensim/commit/3842369190ee421d6d707037f72258c6d0e9513f) - lock non-convex recovery and boundary-layer gates.
- [`4b77149`](https://github.com/Dicklesworthstone/frankensim/commit/4b771496e1b68bedddb6d5ecd292aca34f332cbf) - resolve the `fs-mesh` v3 boundary-layer decision by measurement.
- [`ff9ef32`](https://github.com/Dicklesworthstone/frankensim/commit/ff9ef327a049e5acc22b404dc7653d370dcd77a5) - reject invalid CVaR and objective samples across risk consumers.
- [`38f871f`](https://github.com/Dicklesworthstone/frankensim/commit/38f871ff568f9e84700c7c08606466a4d911e8e3) - fail closed on non-finite probe budget inputs.
- [`fdffa73`](https://github.com/Dicklesworthstone/frankensim/commit/fdffa73531ab7b00e72e1e5535ca4a8a5bdf194b) - include provenance in evidence package roots.
- [`91478e2`](https://github.com/Dicklesworthstone/frankensim/commit/91478e28d420746aca0f835aff109096daa60fd5) - remove unreachable string-write panics from `fs-package`.
- [`462f0f6`](https://github.com/Dicklesworthstone/frankensim/commit/462f0f638dbb348b23e9edb39b87a4a60eba7787) - remove panic arms from `fs-io` conformance checks.
- [`6725739`](https://github.com/Dicklesworthstone/frankensim/commit/6725739f42b878b310e7e7e318fb8c980cef71f3) - expose live flagship browser pipelines in `fs-wasm`.

## 17. Determinism, Cancellation Ownership, GEMM Evidence, And Fail-Closed Fixes

The window from
[`6725739`](https://github.com/Dicklesworthstone/frankensim/commit/6725739f42b878b310e7e7e318fb8c980cef71f3)
through
[`e993e76`](https://github.com/Dicklesworthstone/frankensim/commit/e993e7640a547ca9b11ded6d580a3ce6846a4c82)
adds 119 committed changes across 169 files, with 11,967 insertions and 1,317
deletions. This update is still a moving `main` snapshot, not a tagged release,
but it materially hardens the lower-level semantics that later simulation
claims rely on.

### Delivered capability

- `fs-math` now pins deterministic integer-power behavior and extends `xtask`
  policy checks so dependent crates do not drift back to platform-sensitive
  `.powi` calls without an explicit review point.
- `fs-exec` and `fs-rand` bind stream replay to declared run identity rather
  than execution-pool history, preserving stochastic replay when scheduling or
  worker counts change.
- `fs-race` and `fs-session` now require caller-owned cancellation gates for
  race losers and memory-pressure pauses. Cancellation and pause/resume
  authority is visible to the owner instead of hidden inside convenience
  registrations.
- `fs-exec` wraps solver snapshots in a versioned, self-authenticating envelope,
  giving resumable solver state an artifact boundary suitable for pause,
  migration, fork, and replay checks.
- `fs-la` extends the GEMM performance program with packed f32 and
  mixed-precision BLIS-style paths, transposed/strided op-form GEMM, measured
  blocking defaults, batched small-dense perf tests, and roofline regression
  evidence lanes.
- `fs-uq` now computes empirical CVaR with fractional boundary weighting, which
  avoids under-reporting risk at quantile boundaries and locks the behavior with
  regression tests.
- `fs-adjoint` fails closed on vacuous certificate evidence and retains
  executable regressions for explain and DWR acceptance paths.
- `fs-io` accepts legal face-before-vertex PLY element order while keeping
  payload validation in place, separating valid file ordering from malformed
  face data.
- `fs-render`, `fs-geom`, and small lint/test cleanups remove dependency,
  Lipschitz-bound, and boundary underflow traps that could otherwise surface as
  confusing downstream evidence failures.

### Closed and active workstreams

- [`frankensim-xlvx`](https://github.com/Dicklesworthstone/frankensim/blob/e993e7640a547ca9b11ded6d580a3ce6846a4c82/.beads/issues.jsonl) - closed GEMM perf slice covering AVX/NEON capsules, autotuned blocking, f32/mixed packing, strided forms, and parallel tiling.
- [`frankensim-zsvk`](https://github.com/Dicklesworthstone/frankensim/blob/e993e7640a547ca9b11ded6d580a3ce6846a4c82/.beads/issues.jsonl) - closed CVaR risk under-reporting bug by switching to standard empirical CVaR with fractional boundary weight.
- [`frankensim-epic-morph-wqd.25.1`](https://github.com/Dicklesworthstone/frankensim/blob/e993e7640a547ca9b11ded6d580a3ce6846a4c82/.beads/issues.jsonl) - closed legal PLY element-order import bug.
- [`frankensim-9sf6`](https://github.com/Dicklesworthstone/frankensim/blob/e993e7640a547ca9b11ded6d580a3ce6846a4c82/.beads/issues.jsonl) - closed adjoint certification soundness work so verified color honors finite-difference falsifiers.
- [`frankensim-epic-helm-gp3.13`](https://github.com/Dicklesworthstone/frankensim/blob/e993e7640a547ca9b11ded6d580a3ce6846a4c82/.beads/issues.jsonl) - in progress session memory-pressure gate ownership.
- [`frankensim-1hhx`](https://github.com/Dicklesworthstone/frankensim/blob/e993e7640a547ca9b11ded6d580a3ce6846a4c82/.beads/issues.jsonl) - in progress race cancellation-gate ownership.

### Representative commits

- [`8c6db27`](https://github.com/Dicklesworthstone/frankensim/commit/8c6db27329488ff61ae3be9a33c1ba3591ea5e95) - pin deterministic integer powers and extend the policy surface that catches drift.
- [`921c486`](https://github.com/Dicklesworthstone/frankensim/commit/921c486e400ad918b9853e2ffd5a9a6f681d4658) - bind RNG streams to declared run identity rather than pool history.
- [`d31573b`](https://github.com/Dicklesworthstone/frankensim/commit/d31573bd7ddcb82a1cac0906d471af303ca8c38f) - make races panic-total and require cancellation wiring.
- [`98ea5db`](https://github.com/Dicklesworthstone/frankensim/commit/98ea5dbc2fc0f783b4a1753e1d074f5eb4dc0098) - require caller-owned cancellation gates in `fs-race`.
- [`361bb36`](https://github.com/Dicklesworthstone/frankensim/commit/361bb360e49b82f59b43b220bd6170cd9fb04d0d) - bind session pressure pauses to owned gates.
- [`b2cb2c2`](https://github.com/Dicklesworthstone/frankensim/commit/b2cb2c24796f1f43be972c1b59b456fee6f5cd3d) - wrap solver snapshots in a versioned self-authenticating envelope.
- [`af0339e`](https://github.com/Dicklesworthstone/frankensim/commit/af0339e2f0f6356faf4a901f788579351f2402aa) - add packed f32 and mixed-precision GEMM paths.
- [`dbbffa8`](https://github.com/Dicklesworthstone/frankensim/commit/dbbffa83901a6554478e3a850824c9ca231eed16) - add transposed and strided op-form GEMM.
- [`5b8aeb7`](https://github.com/Dicklesworthstone/frankensim/commit/5b8aeb71737a758a0eedb22d9e9bf1702f295acc) - wire GEMM and batched perf evidence lanes.
- [`7f6420f`](https://github.com/Dicklesworthstone/frankensim/commit/7f6420f585d4a2baf8a5cca4f9112ca5d1c0ca7c) - weight CVaR boundary samples fractionally.
- [`4fbdefc`](https://github.com/Dicklesworthstone/frankensim/commit/4fbdefc124c4b7e98617c467be24757fd6c51391) - accept legal face-before-vertex PLY element order.
- [`e993e76`](https://github.com/Dicklesworthstone/frankensim/commit/e993e7640a547ca9b11ded6d580a3ce6846a4c82) - cover adjoint certificate fail-closed regressions.

## 18. Golden Governance, Certificate Soundness, And Failure Compounding

The 26-commit window from
[`e993e76`](https://github.com/Dicklesworthstone/frankensim/commit/e993e7640a547ca9b11ded6d580a3ce6846a4c82)
through
[`81c11df`](https://github.com/Dicklesworthstone/frankensim/commit/81c11df36fb421d0e7bbfede4c094c96602dce6f)
turns several determinism and certification incidents into governed repository
surfaces. It spans 87 files, with 4,913 insertions and 660 deletions. The
important change is not just more tests: semantic dependencies, replay
identities, package payloads, and minimized failures became inspectable data.

### Delivered capability

- `golden-couplings.json`, semantic-surface version constants, and `xtask
  check-goldens` now identify which downstream sentinels must be deliberately
  re-frozen when a bit-producing surface changes. `docs/GOLDEN_POLICY.md`
  requires committed-tree reproduction, both build modes, a plausible
  bit-moving cause, and a same-commit registry update.
- `xtask check-claims` ties README-cited hashes, `fs-*` crate names, and named
  sentinel functions to the tracked tree instead of trusting manually
  maintained capability prose.
- Deterministic-power enforcement expanded across workspace surfaces, while
  dimension exponent arithmetic now refuses overflow instead of narrowing,
  wrapping, or saturating silently.
- `fs-evidence::Certified<T>` became an opaque validated trust boundary, and
  `fs-sos` separated algebraic decomposition checks from rigorous global or
  radius-scoped value bounds. Under-calibrated conformal bands now return an
  honest unbounded interval rather than finite under-coverage.
- Evidence-package schema v2 began carrying complete color payloads,
  provenance, magnitude budgets, and roots through a strict parser. The
  solver-free checker re-derives budgets and roots, rejects hostile or
  incomplete transport, and treats signature validity as an injected
  capability rather than inferring it from signature presence.
- `fs-obs::ident` introduced typed, versioned, length-prefixed replay identity
  with float-bit encoding, explicit exclusions, dependency roots, and
  collision-boundary tests. Flagship metric identities migrated onto that
  surface and their expected golden changes were re-frozen deliberately.
- Topology fixes aligned closed-cube component counting with 26-connectivity
  and made persistence penalties count components alive at the requested
  level, removing phantom tunnels and historical-bar over-counting.
- `fs_bisect::compound` now turns a counterexample into a deterministic
  regression family: capture, greedy shrinking, bounded neighborhood probes,
  content-addressed manifests, tracking references, replay, and stale-family
  detection. Its first fixture reduces the real `.powi` incident to the exact
  `k = 7` divergence boundary and reproduces the manifest on both reference
  ISAs in debug and release.

Package checking in this slice proves structural and content-root integrity; it
does not rerun solvers or establish scientific truth. The schema-v2 in-tree
default also makes no cryptographic signature-validity claim. Failure
compounding returns a minimum under the caller's shrink order, not a global
minimum, and it records recommended admission rules without enacting them.

### Closed and active workstreams

- [`frankensim-golden-coupling-discipline-y4pt`](https://github.com/Dicklesworthstone/frankensim/blob/81c11df36fb421d0e7bbfede4c094c96602dce6f/.beads/issues.jsonl) - closed golden dependency registry and justified-bump protocol.
- [`frankensim-claim-state-lint-06yc`](https://github.com/Dicklesworthstone/frankensim/blob/81c11df36fb421d0e7bbfede4c094c96602dce6f/.beads/issues.jsonl) - closed README/source claim-state coupling.
- [`frankensim-powi-build-mode-determinism-4xnt`](https://github.com/Dicklesworthstone/frankensim/blob/81c11df36fb421d0e7bbfede4c094c96602dce6f/.beads/issues.jsonl) - closed deterministic integer-power repair and four-quadrant sentinel.
- [`frankensim-epic-epistype-qmao.6.1`](https://github.com/Dicklesworthstone/frankensim/blob/81c11df36fb421d0e7bbfede4c094c96602dce6f/.beads/issues.jsonl) - closed schema-v2 package round-trip and solver-free reverification.
- [`frankensim-epic-gauntlet-6nb.9`](https://github.com/Dicklesworthstone/frankensim/blob/81c11df36fb421d0e7bbfede4c094c96602dce6f/.beads/issues.jsonl) - closed failure-compounding workflow with cross-ISA manifest evidence.
- [`frankensim-epic-helm-gp3.14`](https://github.com/Dicklesworthstone/frankensim/blob/81c11df36fb421d0e7bbfede4c094c96602dce6f/.beads/issues.jsonl) - remained in progress because canonical replay identity had landed as a shared primitive and selected adoption, not a completed workspace-wide migration.

### Representative commits

- [`10bdac0`](https://github.com/Dicklesworthstone/frankensim/commit/10bdac0a123543d8c14512a5244470447af02cb0) - add the golden-coupling registry and justified-bump protocol.
- [`2ca177e`](https://github.com/Dicklesworthstone/frankensim/commit/2ca177ee8bd6157048e6d3cd3a06cfb858cc65cf) - add README claim-state drift checks.
- [`3712602`](https://github.com/Dicklesworthstone/frankensim/commit/3712602e520a8de14276773a4dc2954bd5127f4d) - refuse overflowing dimensional exponent arithmetic.
- [`d0c1e09`](https://github.com/Dicklesworthstone/frankensim/commit/d0c1e092c6d0b48ff95adf34dd2f2e2b2beeca07) - make certified evidence opaque and validated.
- [`049538f`](https://github.com/Dicklesworthstone/frankensim/commit/049538f05bc718738824ee55e68a5965542a1bfb) - separate SOS algebra verification from certified value bounds.
- [`914ff55`](https://github.com/Dicklesworthstone/frankensim/commit/914ff55fe51452c50c04b78fdd89f05139843826) - introduce canonical typed replay identities.
- [`f882987`](https://github.com/Dicklesworthstone/frankensim/commit/f882987b5912a3513af2aa5e036c9e803c3f26e2) - complete schema-v2 package/checker transport and robust-SOS integration.
- [`81c11df`](https://github.com/Dicklesworthstone/frankensim/commit/81c11df36fb421d0e7bbfede4c094c96602dce6f) - close failure compounding with four-quadrant evidence.

## 19. Reproducible Bootstrap, Canonical Replay, And Honest Performance

The 43 commits from
[`81c11df`](https://github.com/Dicklesworthstone/frankensim/commit/81c11df36fb421d0e7bbfede4c094c96602dce6f)
through
[`cdb62ee`](https://github.com/Dicklesworthstone/frankensim/commit/cdb62eea326e523fed26048d4d53c717bf43d455)
span 111 files, 6,443 insertions, and 836 deletions. This wave made
bootstrap, replay identity, package semantics, hardware placement, and
performance admission more reproducible while retaining measured misses as
explicit evidence.

### Delivered capability

- `xtask bootstrap-constellation` now materializes locked dependencies into
  detached, transform-free checkouts; verifies existing heads and cleanliness;
  supports offline caches and mirror transport; distinguishes provable
  case-collision artifacts from unrelated dirt; and writes bootstrap
  provenance. A layer inversion found during this work moved package coverage
  out of `fs-crosswalk` and into the legal L6 dependency direction.
- DSR quality coverage gained Cargo-metadata-derived feature lanes, visible
  `required-features` skips, standalone `fs-wasm` native/build lanes, complete
  per-lane logs with tested-tree identity, and derived repository inventories.
  The inventory is also checked against README crate, contract, and test counts.
- Canonical typed replay identity spread into capability fingerprints,
  reports, flagship metric streams, and steered-study fingerprints. Vessel,
  LBM, and BEM paths replaced platform `sin`, `cos`, `powf`, `atan2`, `ln`,
  and related calls with deterministic math, restoring the poured-mass and
  ornithoid ROA metrics to their full cross-mode hash streams.
- Evidence-package schema v3 added backward-only composition receipts,
  falsifier records, and content-hashed dataset anchors. The checker recomposes
  parent evidence, refuses color laundering and refuted claims, and keeps those
  semantics bound into the package root without importing solver code.
- Package admission was tightened across in-memory verification, JSON parsing,
  and standalone checker diagnostics: complete provenance, unique claim IDs,
  finite/ordered payloads, exact integer handling, magnitude-overflow refusal,
  and bounded bytes, depth, nodes, containers, tokens, and strings.
- `fs-sparse` gained optional FrankenNumpy dense-array conversion interop with
  explicit copy/densification semantics, shape and overflow refusal,
  non-finite rejection, and documented treatment of stored and signed zeros.
- Linux affinity syscalls were isolated in a registered capsule; measured L3
  groups, first-touch audits, advisory TilePool CCD pinning, hugepage/THP
  provenance, and A/B harnesses made machine placement measurable. Per-kernel
  SIMD tier audits found and repaired a slower NEON scale path and an x86
  fused-multiply-add path that had fallen back to per-element libm calls.
- Autotuner persistence began refusing foreign probes, invalid pins,
  noncanonical records, and refresh-counter exhaustion. Per-worker throughput
  rows also demonstrated that bandwidth-rich and bandwidth-starved machines
  can select different schedule families without changing result bits.
- The roofline gate gained an `EnvironmentInvalid` result after a saturated
  host collapsed both the measured machine roof and the kernel by roughly the
  same factor, producing a vacuous passing ratio. Absolute plausibility floors,
  over-roof detection, and run-wide poisoning now prevent that false
  performance certificate.
- FFT moved to mixed radix-8/4/2 Stockham stages and improved raw M4 throughput
  by 11–18%, but the corrected traffic model left attainment around 16–25%; the
  40% gate remained open. Batched f32 and mixed f32-storage/f64-accumulation
  GEMM also landed without claiming the still-unmet 60% small-dense target.

Affinity, hugepages, prefetch distance, and schedule selection are
machine-specific measurement surfaces, not portable speed guarantees. The
bootstrap fetch path was exercised through mirror transport rather than a new
clean remote host. Package verification remains solver-free structural and
semantic checking, and at this point package roots were still the in-house FNV
content address rather than the later BLAKE3 schema.

### Closed and active workstreams

- [`frankensim-epic-foundations-huq.17`](https://github.com/Dicklesworthstone/frankensim/blob/cdb62eea326e523fed26048d4d53c717bf43d455/.beads/issues.jsonl) - closed lock-driven constellation bootstrap.
- [`frankensim-epic-foundations-huq.18`](https://github.com/Dicklesworthstone/frankensim/blob/cdb62eea326e523fed26048d4d53c717bf43d455/.beads/issues.jsonl) - closed derived feature-matrix, standalone-workspace, and inventory coverage.
- [`frankensim-epic-helm-gp3.14`](https://github.com/Dicklesworthstone/frankensim/blob/cdb62eea326e523fed26048d4d53c717bf43d455/.beads/issues.jsonl) - closed canonical replay identity adoption.
- [`frankensim-xfxq`](https://github.com/Dicklesworthstone/frankensim/blob/cdb62eea326e523fed26048d4d53c717bf43d455/.beads/issues.jsonl) - closed evidence-package schema v3 receipts, falsifiers, and anchors.
- [`frankensim-gtql`](https://github.com/Dicklesworthstone/frankensim/blob/cdb62eea326e523fed26048d4d53c717bf43d455/.beads/issues.jsonl) - closed sparse FrankenNumpy and graph-constellation interop with explicit conversion no-claims.
- [`frankensim-epic-perf-fz2.2`](https://github.com/Dicklesworthstone/frankensim/blob/cdb62eea326e523fed26048d4d53c717bf43d455/.beads/issues.jsonl) - remained active at this endpoint while the machine-adaptive scoreboard was reconciled; it closes later in the next window.
- [`frankensim-27d3`](https://github.com/Dicklesworthstone/frankensim/blob/cdb62eea326e523fed26048d4d53c717bf43d455/.beads/issues.jsonl) - remained active because the FFT roofline target was not met.
- [`frankensim-9ekv`](https://github.com/Dicklesworthstone/frankensim/blob/cdb62eea326e523fed26048d4d53c717bf43d455/.beads/issues.jsonl) - remained open because batched small-dense GEMM had not met its 60% target.

### Representative commits

- [`6fb8513`](https://github.com/Dicklesworthstone/frankensim/commit/6fb85130c893b2e9de652adec124f6fd38d48083) - add lock-driven constellation bootstrap and repair the exposed layer inversion.
- [`f96156d`](https://github.com/Dicklesworthstone/frankensim/commit/f96156dd91662856971bb2494cca7b1f190a1bc3) - derive DSR feature lanes, visible skips, standalone WASM coverage, and inventories.
- [`3a263c2`](https://github.com/Dicklesworthstone/frankensim/commit/3a263c2fb14d94eeda511d1a6e6a4c86ddffbc21) - canonicalize replay identities and deterministic flagship trigonometry.
- [`b24b0dc`](https://github.com/Dicklesworthstone/frankensim/commit/b24b0dc914c74b39f95a8ef2c20693a99dd87edc) - finish deterministic BEM routing and four-quadrant flagship goldens.
- [`e00858c`](https://github.com/Dicklesworthstone/frankensim/commit/e00858c2b97451baef410681e422199cc0963cb4) - add schema-v3 composition receipts, falsifiers, and anchors.
- [`cd185a6`](https://github.com/Dicklesworthstone/frankensim/commit/cd185a6e5808394edafc9884f525c2a999fdbbff) - harden standalone evidence-verification boundaries.
- [`dc99094`](https://github.com/Dicklesworthstone/frankensim/commit/dc990942c6aca864644df5b5860061482e964ab4) - add explicit FrankenNumpy dense-array conversion for sparse matrices.
- [`fb3bfc8`](https://github.com/Dicklesworthstone/frankensim/commit/fb3bfc83320ae26806d9f8e007c38618d9f50057) - add OS affinity, measured L3 topology, and CCD/first-touch A/B evidence.
- [`ceccb00`](https://github.com/Dicklesworthstone/frankensim/commit/ceccb003ba0914b488474a2dfaeb57cd3eedb080) - add advisory TilePool CCD pinning with bit-equivalence gates.
- [`535aa83`](https://github.com/Dicklesworthstone/frankensim/commit/535aa8393aa1484cdef4df4edd36dfe21bf951ac) - improve FFT throughput while recording the unmet roofline target.
- [`cdb62ee`](https://github.com/Dicklesworthstone/frankensim/commit/cdb62eea326e523fed26048d4d53c717bf43d455) - refuse self-normalizing performance gates on saturated machines.

## 20. Evidence-Bound Autotuning, BLAKE3 Packages, And Closure Enforcement

The 66-commit window from
[`cdb62ee`](https://github.com/Dicklesworthstone/frankensim/commit/cdb62eea326e523fed26048d4d53c717bf43d455)
through
[`fb34ef6`](https://github.com/Dicklesworthstone/frankensim/commit/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8)
touch 119 files, with 13,684 insertions and 1,863 deletions. They establish the
first production-facing GEMM measure-cache-model-dispatch path, move package
integrity onto domain-separated BLAKE3 roots, and extend the repository's
failure-compounding doctrine to issue closure itself. Several commits retain
explicit `wip` boundaries; this section preserves those distinctions.

### Delivered capability

- `fs-roofline` moved its citation boundary from plausible-looking attainment
  ratios to exact timed receipts bound to pre-run machine axes, kernel identity,
  work count, raw sample-duration bits, derived statistics, and verdict.
  Admission re-derives every row, corroborates post-run axes, and makes analytic,
  tampered, duplicate, invalid-axis, or drifted measurements non-citable.
- A governed historical-baseline API models the remaining sustained-contention
  hole: trusted axes require at least three mutually agreeing, floor-plausible
  runs plus operator, justification, fingerprint, topology,
  OS/architecture/firmware identity, and age policy. First-run measurements are
  candidates rather than self-authorizing evidence; the API marks degraded,
  suspicious, stale, or identity-drifted environments ineligible. Binding that
  admission result into every citable GEMM receipt remains an active blocker.
- `fs-exec` tuning evidence became versioned, unit-safe, canonical, bounded,
  and explicitly non-statistical. Exact wall-time samples determine a
  deterministic ranked argmin; summaries and candidate separation are
  re-derived rather than trusted as independent claims.
- A bounded GEMM block-plan lattice and `fs-session` orchestration now implement
  cold measurement, validated cache adoption, pinned replay, ledger
  read/write-through, and dispatch. Tune rows bind bit semantics, shape,
  requested threads and thread budget, exact probe dimensions, ISA tier,
  placement, and implementation identity through one canonical key.
- Tuning publication was split into prepare, durable persistence, in-memory row
  commit, transactional dispatch, and decision commit. Cancelled or failed
  dispatches do not publish a success decision, and cache imports must match the
  requested execution identity and the ranked-evidence winner.
- `fs-la` added cancellation-aware parallel f64 GEMM with private output
  staging. Cancellation stops new acquisition, drains scoped workers, and
  preserves caller-visible `C` bit-for-bit; only a successful operation crosses
  the documented non-cancellable final-copy boundary.
- Core and batched GEMM paths added checked extent, stride, packing, allocation,
  and tuning-quantum admission before mutation or unsafe SIMD entry. Batched
  f64, f32, and mixed paths also implement `alpha == 0` without reading
  poisoned inputs. The final endpoint still carries open SIMD-pointer and
  batched-semantics workstreams, so this is landed hardening rather than a
  claim that every safe-facade audit is closed.
- The shipped roofline registry now invokes the session-backed GEMM tuner and
  derives §14.1 implementation coverage from registered kernels. A `landed`
  target row means the registry contains the production path; it does not mean
  the declared 75% performance target was met.
- `fs-substrate` added architecture-specific software-prefetch hints and a
  distance sweep that found different optima on the M4 Pro and Threadripper,
  while proving the hint cannot affect result bits. A repository Cargo note
  also records why globally enabling x86 FMA is rejected: LLVM contraction
  changed a Chebyshev golden, so the libm-FMA trap must be repaired in narrow
  pure-`mul_add` target-feature capsules instead of changing workspace-wide
  arithmetic.
- A dependency-free `fs-blake3` crate became the shared safe-Rust content-hash
  owner. Evidence-package schema v4 moved from 64-bit FNV roots to 32-byte
  BLAKE3 roots with distinct derive-key domains for headers, claims, internal
  nodes, and raw artifacts while keeping the checker solver-free.
- A clean-bootstrap audit narrowed the earlier `xtask` claim: Cargo cannot run
  that command before unresolved path dependencies have been materialized.
  Bootstrap now refuses unverifiable substitutions and the docs distinguish
  prebuilt-cache verification from true clean-clone setup. The FrankenNumpy pin
  was advanced past its case-colliding fuzz-corpus filenames, restoring clean
  checkouts on default macOS, but a strict pre-Cargo bootstrap remains open.
- Evidence waivers now bind to exact canonical color payloads rather than only
  a color rank. Package coverage depends on valid authenticated records, and a
  refused package exposes no positive evidence summary that could be mistaken
  for a partial pass.
- Interval Newton gained bounded work and completeness receipts; high-order IGA
  gained degree-aware Gauss-Legendre quadrature; constrained optimization began
  validating complete KKT state; chart agreement became
  `Agreed`/`Disagreed`/`Unknown`; and voxel queries refuse mixed frames,
  unrepresentable coordinates, and unbudgeted dense transforms.
- The cache-blocked six-step FFT passed its correctness and golden batteries but
  measured roughly twice as slow as the existing stage walk, so it remains
  default-off behind `frontier-sixstep`. FEEC's high-order x86 path removed the
  per-element libm `fma` hole, but scalar fused code generation still left the
  30% cross-ISA vectorization target open.
- `xtask check-closures` became the eighth `check-all` lane. A closed bug Bead
  must name retained regression evidence or an explicit disposition, connecting
  tracker closure to the previously landed failure-compounding workflow.

This snapshot does not close every production claim. The exhaustive GEMM
performance oracle was still manual and exercised the legacy parallel path; the
later hardening commits carried no current x86 production-oracle result; GEMM
staging was not memory-budgeted; and worker scopes still used scoped standard
threads rather than a full asupersync `Cx`. The roofline wrapper could also
persist a reusable session tune row during warmup before post-run roofline
admission, so rejected axes withheld citable roofline output without yet
guaranteeing that no tuning state survived. Content roots prove transit
integrity, not scientific truth, authorship, or claim origin. The closure gate
checks that a close reason names regression evidence or a disposition; it does
not resolve the cited artifact and independently prove that artifact passes.

### Closed and active workstreams

- [`frankensim-yqug`](https://github.com/Dicklesworthstone/frankensim/blob/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8/.beads/issues.jsonl) - closed at [`f93c9de`](https://github.com/Dicklesworthstone/frankensim/commit/f93c9de48367ec43a1ac4ab3df7e24a1050882f0) after the first real GEMM autotune loop, then reopened by the final endpoint for receipt/decision binding, real `TilePool`/`Cx` traversal, exact build/effective-tier identity, and no-write-before-admission blockers.
- [`frankensim-epic-perf-fz2.2`](https://github.com/Dicklesworthstone/frankensim/blob/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8/.beads/issues.jsonl) - closed the machine-adaptive hardening scoreboard with measured per-ISA outcomes and remaining gap beads.
- [`frankensim-hx4p`](https://github.com/Dicklesworthstone/frankensim/blob/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8/.beads/issues.jsonl) - closed the self-expanding Gauntlet slice by reconciling failure compounding with universal closure-evidence enforcement.
- [`frankensim-7uq9`](https://github.com/Dicklesworthstone/frankensim/blob/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8/.beads/issues.jsonl) and [`frankensim-t7x3`](https://github.com/Dicklesworthstone/frankensim/blob/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8/.beads/issues.jsonl) - remained in progress in the tracker despite the shared BLAKE3 owner and schema-v4 root migration landing materially.
- [`frankensim-krym`](https://github.com/Dicklesworthstone/frankensim/blob/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8/.beads/issues.jsonl) - remained open for schema-v5 non-forgeable claim origins and authenticated waiver transport.
- [`frankensim-27d3`](https://github.com/Dicklesworthstone/frankensim/blob/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8/.beads/issues.jsonl) - remained open after the six-step FFT was honestly rejected for default promotion.
- [`frankensim-a55x`](https://github.com/Dicklesworthstone/frankensim/blob/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8/.beads/issues.jsonl) - closed the narrow FEEC libm-call defect, while [`frankensim-cwjn`](https://github.com/Dicklesworthstone/frankensim/blob/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8/.beads/issues.jsonl) remained open for packed vectorization and the 30% target.
- [`frankensim-9ekv`](https://github.com/Dicklesworthstone/frankensim/blob/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8/.beads/issues.jsonl) - remained open for the batched small-dense 60% target.
- [`frankensim-rpgc`](https://github.com/Dicklesworthstone/frankensim/blob/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8/.beads/issues.jsonl) - closed core GEMM extent-overflow admission, while [`frankensim-zevq`](https://github.com/Dicklesworthstone/frankensim/blob/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8/.beads/issues.jsonl) and [`frankensim-9ry9`](https://github.com/Dicklesworthstone/frankensim/blob/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8/.beads/issues.jsonl) remained in progress for SIMD pointer-safety and batched no-read/overflow closure.
- [`frankensim-dfh3`](https://github.com/Dicklesworthstone/frankensim/blob/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8/.beads/issues.jsonl) - closed the governed historical machine-axis baseline layer; receipt binding remains under the reopened GEMM blockers.
- [`frankensim-1t8i`](https://github.com/Dicklesworthstone/frankensim/blob/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8/.beads/issues.jsonl) - remained open for a strict pre-Cargo clean-clone bootstrap; [`frankensim-7n2n`](https://github.com/Dicklesworthstone/frankensim/blob/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8/.beads/issues.jsonl) also remained open despite the case-collision re-pin landing in [`fb34ef6`](https://github.com/Dicklesworthstone/frankensim/commit/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8).

### Representative commits

- [`2ecb0cd`](https://github.com/Dicklesworthstone/frankensim/commit/2ecb0cd5396a5c632c64441f6527844a94cc526e) - bind citable roofline runs to re-derivable timed receipts.
- [`dc859e2`](https://github.com/Dicklesworthstone/frankensim/commit/dc859e2504aa0249937d74c8a10f0f41993f0702) - introduce bounded GEMM plans and the session autotune scaffold.
- [`666967d`](https://github.com/Dicklesworthstone/frankensim/commit/666967def6b9827ff1862b2458b66a72afb70632) - make tuner evidence unit-safe, canonical, bounded, and non-statistical.
- [`f91b775`](https://github.com/Dicklesworthstone/frankensim/commit/f91b7753f9c665f468ee6edd571f66f3996311a7) - add transactional cancellation-aware parallel GEMM.
- [`3096ee5`](https://github.com/Dicklesworthstone/frankensim/commit/3096ee5bc110504ed5b04b56352d2eabffa3ead8) - bind tune rows and decisions to completed execution identity.
- [`da8f438`](https://github.com/Dicklesworthstone/frankensim/commit/da8f438fb7c442500bf62e94bb14380686dbfddf) - wire the session-backed tuner into the shipped roofline registry.
- [`f93c9de`](https://github.com/Dicklesworthstone/frankensim/commit/f93c9de48367ec43a1ac4ab3df7e24a1050882f0) - derive target coverage from the production registry while separating coverage from performance proof.
- [`3a4257a`](https://github.com/Dicklesworthstone/frankensim/commit/3a4257a76c031eaa93f4763ed9e39c6de061d2f7) - add bit-neutral software prefetch and per-machine distance sweeps.
- [`2d6d3f9`](https://github.com/Dicklesworthstone/frankensim/commit/2d6d3f99bc103684d33c183589d6d7f19a6f3828) - add governed historical roofline baselines and sustained-contention admission semantics.
- [`005f0e2`](https://github.com/Dicklesworthstone/frankensim/commit/005f0e2b4570f200006ba4d1805fca61d86b7789) - document why global x86 FMA contraction violates G5 bit determinism.
- [`9711705`](https://github.com/Dicklesworthstone/frankensim/commit/971170581e555ec44ede17f6bec3612aa023ed78) - extract the shared dependency-free BLAKE3 owner.
- [`10d1e2c`](https://github.com/Dicklesworthstone/frankensim/commit/10d1e2cfb80dc1bed6df237cc1695364d912d61f) - migrate package/checker roots to domain-separated 32-byte BLAKE3.
- [`d4b9c04`](https://github.com/Dicklesworthstone/frankensim/commit/d4b9c04c74dc173bb6f47d735f9f999a9483ec0f) - retain the correct-but-slower six-step FFT as gated frontier evidence.
- [`acedb1b`](https://github.com/Dicklesworthstone/frankensim/commit/acedb1b87e686013173081683e8de79068bf1e1f) - remove the FEEC x86 libm-FMA call path while preserving an explicit performance no-claim.
- [`504d2a8`](https://github.com/Dicklesworthstone/frankensim/commit/504d2a89003913e04efad13a6e9d0c464f3dd04e) - require closed bugs to cite regression evidence or an explicit disposition.
- [`fb34ef6`](https://github.com/Dicklesworthstone/frankensim/commit/fb34ef6b6fb2b57f6c55e2ec8c4cc33653a958f8) - re-pin FrankenNumpy after removing case-colliding corpus filenames upstream.

## 21. SIMD Transpose, FEEC Accumulation, And Pre-Cargo Bootstrap

The seven-commit follow-up from
[`f7e9787`](https://github.com/Dicklesworthstone/frankensim/commit/f7e9787f7dcf412b35642b7fffb7e943d1c64045)
through
[`291c7db`](https://github.com/Dicklesworthstone/frankensim/commit/291c7dbfdc6fd366f6ee55a6dbc39f137e05afd6)
touches 19 files, with 1,043 insertions and 104 deletions. It advances two
measured numerical kernels without erasing their remaining no-claim boundaries,
and it moves clean-machine constellation materialization out of the Cargo
workspace so it can run before sibling path dependencies exist.

### Delivered capability

- `fs-simd` gained an 8-by-8-tiled NEON transpose for interleaved complex f64
  matrices. One complex value occupies one 128-bit vector register, so the
  capsule performs exact loads and stores without floating-point arithmetic or
  shuffle-induced bit changes. The scalar twin and tier battery cover irregular
  and tile-aligned sizes, while `fs-fft` routes its six-step transpose through
  the shared operation table.
- The capsule improved the large-transform six-step/stage-walk ratio from the
  earlier approximately 0.44 to measured values of 0.731 and 0.679 on the
  macOS-aarch64 lane. That halves the deficit but does not produce a win, so
  `frontier-sixstep` remains default-off and the existing
  `0x79aa_108f_a517_012f` six-step golden remains unchanged. The tracker released
  `frankensim-27d3` for its next measured lever rather than declaring the
  performance target complete.
- A post-landing audit found that the safe scalar and NEON transpose facades
  compute `2 * n1 * n1` with unchecked `usize` arithmetic before validating
  slice lengths. The endpoint therefore records a measured optimization, not a
  completed safe-facade proof; overflow-safe geometry admission and a retained
  regression are still required before treating the capsule as fully closed.
- `fs-feec` rewrote all three sum-factorized hexahedral contractions to keep
  each output row in a local `[f64; P]` accumulator and store once. The
  per-output operation chain remains `+0.0` followed by ascending-`l`
  `mul_add`, preserving the `0xaaf1_076a_196c_6902` golden while removing repeated
  destination loads/stores and making the const-`P` loops friendlier to SLP
  vectorization.
- The FEEC change measured 0.56-0.70 attainment across `r = 1..5` on the M4
  lane, up from the prior 0.44 row. The commit deliberately leaves
  `frankensim-cwjn` in progress: the x86 gate verdict and the declared both-ISA
  30% target were not complete at this endpoint.
- A follow-up replaced iterator/zip chains inside those register accumulators
  with explicit `0..P` index loops. This exposes constant loop bounds more
  directly to x86 SLP while retaining the same ascending-`l` fused chain; the
  release battery remained 10/10 with the frozen golden unchanged. M4
  remeasurement and the quiet x86 verdict were still pending.
- A standalone `tools/bootstrap` package now lives outside the root workspace,
  has its own empty workspace and lockfile, and uses no dependencies. It can
  parse `constellation.lock`, clone missing siblings at detached pinned commits,
  verify existing sibling HEAD and cleanliness, honor a network-free offline
  mode, use a mirror base, and write deterministic schema-v1 provenance without
  first asking Cargo to resolve the very path dependencies it is meant to
  materialize.
- Hermetic replay drills cover first clone, idempotent offline verification,
  dirty and drift refusal, and offline-missing refusal. The committed close
  evidence also records a real seven-sibling offline verification including the
  case-collision-safe FrankenNumpy pin. `frankensim-1t8i` and
  `frankensim-7n2n` are closed at this endpoint.

The bootstrap closure is narrower than a claim that every hostile invocation
is already fail-closed. The endpoint README and parts of the bootstrap/CI docs
still retain earlier language saying clean-clone bootstrap is unavailable or
open. Code review also found that the newly cloned path re-checks the detached
HEAD but does not apply the existing-sibling clean-tree check after checkout,
and a bare `--root` or `--from` silently falls back instead of refusing a
missing value. These are follow-up documentation and admission defects, not
reasons to erase the landed standalone bootstrap and replay evidence.

### Closed and active workstreams

- [`frankensim-1t8i`](https://github.com/Dicklesworthstone/frankensim/blob/291c7dbfdc6fd366f6ee55a6dbc39f137e05afd6/.beads/issues.jsonl) - closed with the standalone zero-dependency bootstrap, hermetic replay drills, and real-lock offline verification; a literal blank-machine fetch from the public remotes remained outside the recorded session.
- [`frankensim-7n2n`](https://github.com/Dicklesworthstone/frankensim/blob/291c7dbfdc6fd366f6ee55a6dbc39f137e05afd6/.beads/issues.jsonl) - closed after the FrankenNumpy corpus rename and pinned clean checkouts on case-insensitive macOS and case-sensitive Linux.
- [`frankensim-27d3`](https://github.com/Dicklesworthstone/frankensim/blob/291c7dbfdc6fd366f6ee55a6dbc39f137e05afd6/.beads/issues.jsonl) - open after the NEON transpose improvement because the six-step FFT still loses to the stage walk and remains gated off.
- [`frankensim-cwjn`](https://github.com/Dicklesworthstone/frankensim/blob/291c7dbfdc6fd366f6ee55a6dbc39f137e05afd6/.beads/issues.jsonl) - in progress after the register-accumulator and const-bound-loop rewrites; the committed M4 row improved, while the x86 verdict remained pending.
- [`frankensim-dfh3`](https://github.com/Dicklesworthstone/frankensim/blob/291c7dbfdc6fd366f6ee55a6dbc39f137e05afd6/.beads/issues.jsonl) - open again with a coordination note that production historical-baseline wiring was in flight on another lane.

### Representative commits

- [`f7e9787`](https://github.com/Dicklesworthstone/frankensim/commit/f7e9787f7dcf412b35642b7fffb7e943d1c64045) - add the NEON interleaved-complex transpose capsule and record its measured, still-insufficient six-step improvement.
- [`ea03f01`](https://github.com/Dicklesworthstone/frankensim/commit/ea03f0125e3100e14dd93449ab0109adf9d4906c) - record the 27d3 capsule slice and release the still-open workstream.
- [`2168810`](https://github.com/Dicklesworthstone/frankensim/commit/2168810f3544652cb6e53dc46c8c488930a0b741) - make FEEC sum-factorized contractions register-accumulator friendly without moving their golden.
- [`a3edf0e`](https://github.com/Dicklesworthstone/frankensim/commit/a3edf0ec681d679865d5e5796f0802e1a4c3f48b) - add the standalone zero-dependency constellation bootstrap and its replay harness.
- [`2568c82`](https://github.com/Dicklesworthstone/frankensim/commit/2568c8262bf50789d20be528265b5c96d575fb1b) - close the clean-machine bootstrap tracker slice with explicit evidence and a residual blank-machine no-claim.
- [`291c7db`](https://github.com/Dicklesworthstone/frankensim/commit/291c7dbfdc6fd366f6ee55a6dbc39f137e05afd6) - expose const-bound FEEC accumulation loops to x86 SLP while preserving the golden and open cross-ISA verdict.

## Current Non-Release Status

- No package or crate release has been tagged yet.
- The canonical project state is `main`, not a versioned artifact; the latest
  factual checkpoint covered here is
  [`55379ba`](https://github.com/Dicklesworthstone/frankensim/commit/55379baf41cf1d3584c87a81ae1dde20508a4c8f),
  while the complete capability-wave reconstruction currently ends at
  [`291c7db`](https://github.com/Dicklesworthstone/frankensim/commit/291c7dbfdc6fd366f6ee55a6dbc39f137e05afd6).
- The repository is actively changing; use crate `CONTRACT.md` files and Beads
  close reasons for detailed no-claim boundaries.
- The README describes the implemented workspace; the long-form plan remains
  `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md`.
- Uncommitted work in the local working tree is not covered by this changelog.

## Notes For Agents

- Start with the version timeline for chronology.
- Use thematic sections for architectural orientation.
- Use commit links for implementation evidence and Beads links for intent and
  acceptance criteria.
- Do not cite GitHub Releases for this repo until a real release exists.
- Treat DSR as the project validation source of truth unless `AGENTS.md`
  changes.
