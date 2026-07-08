# fs-topopt CONTRACT

## Purpose and layer

Layer: **L4 ASCENT** (deps: fs-adjoint/fs-solver/fs-feec L3,
fs-ascent L4, fs-material L3, fs-la/fs-sparse L1, fs-math L0,
fs-rep-mesh L2). Density-based topology optimization (plan §9.5 [S]):
SIMP with the modern hygiene stack — Helmholtz PDE filtering,
Heaviside projection with continuation, exact chain-rule
sensitivities, and the classical optimality-criteria driver.

NAMING: the plan's atlas used "fs-topo" for this stack; that crate
name carries the L2 topology-CERTIFICATE machinery (persistence,
cubical homology), so the optimization stack lives here.

## Public types and semantics

- `DensityElasticity` — matrix-free K(ρ̄) = Σ_c E_c·K_c with per-cell
  UNIT-modulus 12×12 blocks (fs-material `IsotropicElastic` tangent ×
  fs-feec barycentric-gradient B-matrices) kept separate for the
  exact chain rule; Dirichlet dofs handled by identity-on-fixed
  masking (SPD on the full vector space); `cell_energies` =
  uᵀK_c u per cell (the compliance sensitivity kernel).
- `DensityFilter` — Helmholtz filter: volume-weighted cell→vertex
  scatter, (M + r²K)⁻¹M solve on the FULL vertex space (natural BCs —
  the correct filter behavior: no boundary droop), vertex→cell
  gather. Linear; `apply_transpose` is the exact chain-rule pullback
  (adjointness ⟨Fx, w⟩ = ⟨x, Fᵀw⟩ G0-tested; constants preserved to
  solver tolerance). One assembled operator, built once.
- `heaviside`/`heaviside_derivative` — tanh projection with β/η,
  exact endpoints, monotone, closed-form slope; tanh through the
  strict exp kernel (no platform libm in the pipeline).
- `DesignPipeline` — ρ → filter → projection → SIMP
  (E_min + (1−E_min)·ρ̄^p) with `pullback` reversing the chain
  exactly, and `compliance_and_gradient` exploiting self-adjointness
  (λ = u ⇒ dc/dE_c = −u_cᵀK_cu_c: ZERO extra solves — stated and
  FD-verified).
- `optimality_criteria` — the classical OC driver (documented choice
  for compliance/volume; fs-ascent's augmented Lagrangian is the
  general path): multiplicative update with move limits, volume
  multiplier by fixed 80-step bisection on the PROJECTED volume —
  fully deterministic, whole runs replay bitwise.
- `robust::{RobustPipeline, robust_optimality_criteria}` — the
  erode/dilate three-field formulation: one filter, three projections
  (η ± δ), POINTWISE-ORDERED realizations (tested); minimize ERODED
  compliance s.t. volume on the DILATED field whose target is adapted
  (DAMPED, 0.3 blend) so the NOMINAL design carries the budget — the
  undamped adaptation measured a period-2 limit cycle on cold starts
  (kept as a regression probe); reports carry the erosion-retention
  ratio vol(eroded)/vol(nominal), the measured minimum-length-scale
  signal.

## Invariants

- Every stage of the density chain has an exact derivative; the
  composed sensitivity is FD-verified at MULTIPLE continuation
  stages (p = 1 → 3, β = 1 → 8) per the acceptance.
- The filter preserves constants and is symmetric in the
  volume-weighted pairing (mesh-independent length scale r).
- OC keeps designs in [1e−3, 1] with move limits; the volume
  constraint tracks the projected design.

## Error model

Structured panics on solver failures and invalid materials
(modeling errors). Optimization outcomes are reported traces
(compliance, volume, final change), never silent.

## Determinism class

Bit-deterministic: fixed bisection schedules, deterministic solves
throughout; a WHOLE topology-optimization run replays bitwise
(G5-tested). Golden FNV-64 over pipeline stages, compliance
gradient, and a short OC run: `0x772a_2f8c_a720_dd64`; robust
three-field golden `0x519a_41e3_466e_4b7d`. Recorded on
Apple M4 Pro, verified on Threadripper (x86_64).

## Cancellation behavior

Iteration-granular through the resumable fs-solver states; OC
iterations are bounded and the driver can stop between them. Cx
wiring is driver scope.

## Unsafe boundary

None. `unsafe_code = "deny"`.

## Feature flags

None.

## Conformance tests

`tests/topopt_battery.rs` (8 cases): filter G0 laws (linearity ≤
1e−9, transpose adjointness ≤ 1e−9, constants preserved); projection
G0 (exact endpoints, monotone on a 100-point sweep, slope vs FD ≤
1e−8); FULL-CHAIN sensitivity vs FD at three continuation stages
(rel ≤ 2e−4 through solve+SIMP+projection+filter); OC cantilever
(kuhn(3), fixed face + edge load): compliance reduced ≥ 20%, volume
within 0.03 of the 0.4 target, design range > 0.5 (not gray), and
the ENTIRE run replaying bitwise; three-field pointwise ordering on
random designs + eroded-compliance sensitivity FD gate (rel ≤ 2e−4);
robust OC vs the non-robust baseline AUDITED WITH THE SAME
three-field probe — eroded compliance descends, volumes ordered,
nominal volume on budget, and erosion retention at least matching
the baseline (the min-length-scale claim, measured); two cross-ISA
golden hashes (slice-1 pipeline + robust). Plus
`tests/probe_robust.rs`: the limit-cycle regression.

## No-claim boundaries

- Slices 1–2 scope: compliance + volume (+ robust three-field) on
  FIXED kuhn meshes. Eigenfrequency objectives (clustered-eigenvalue
  handling), stress p-norm aggregation (singularity-trap treatment),
  the medial-axis thickness oracle (the geometry-layer audit beyond
  the erosion-retention signal), and the CutFEM-octree marquee (zero
  remeshing + DWR adaptivity + literature-benchmark envelopes on
  MBB/L-bracket class fixtures) are the bead's later slices /
  recorded splits.
- OC is the compliance/volume driver; MMA is not implemented
  (fs-ascent AL is the general constrained path — documented
  choice).
- No multi-load/worst-case formulations, no continuation SCHEDULER
  (drivers own β/p ramps; the primitives take fixed parameters).
