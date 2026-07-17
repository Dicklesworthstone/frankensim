# fs-topols CONTRACT

## Purpose and layer

Layer: L4 (ASCENT). Level-set topology optimization (plan §9.5 [S/F],
bead 7tv.12): shape-gradient velocity advection with topological
derivatives for hole nucleation — genuine, mathematically justified
topology changes. The level set IS the geometry: physics is evaluated
by fs-cutfem's canonical elasticity operator directly on the evolving discrete SDF
with zero meshing anywhere in the loop.

## Public types and semantics

- `GridSdf`: nodal φ on a uniform `[0,1]²` lattice, bilinear between
  nodes, implementing fs-cutfem's `CutSdf`. The enclosure is EXACT up
  to an outward roundoff pad — a bilinear attains its extrema at the
  corners of any axis-aligned rectangle, so a box's range is the hull
  of per-cell clipped-corner evaluations (certified classification
  holds on the moving geometry). Horizontal, vertical, and point boxes
  lying exactly on lattice lines still visit an adjacent cell and return a
  finite enclosure; they do not degrade to `Interval::WHOLE`.
- `weno::advect` + `Velocity::{Linear, Normal}`: WENO5 (Jiang–Shu)
  spatial stencils with TVD-RK3 time stepping under a CFL number;
  linear advection for the order battery, Godunov/Rouy–Tourin normal
  flow `φ_t + v_n|∇φ| = 0` for the optimizer; band-masked updates.
- `fim::redistance` → `RedistanceAudit`: fast-iterative-method eikonal
  relaxation (Godunov upwind updates, alternating-order Gauss–Seidel
  to an order-independent fixed point) with the interface FROZEN at
  values reconstructed from zero crossings; audits carry the sampled
  zero-set Hausdorff drift (units of h) and |∇φ|−1 statistics — the
  redistancing-frequency policy's inputs. `zero_crossings`/`hausdorff`
  are the audit primitives.
- `veloext::extend_velocity`: interface-normal extension by
  ascending-|φ| upwind sweeps (deterministic ordering).
- `topder::topological_derivative`: the compliance sensitivity to
  infinitesimal traction-free hole insertion (plane-strain
  Amstutz-form coefficients), NUMERICALLY GATED — the battery punches
  a real hole and checks the measured compliance change against
  `DT·πρ²` within a documented first-order band, so wrong signs or
  scales cannot ship. `nucleate`: greedy best-gain hole punching with
  spacing, box-edge margins (clamps/loads excluded), and per-event
  ledger rows (`NucleationEvent`).
- `optimize::optimize_compliance` + `OptimizeSettings`/`Cantilever` →
  `OptimizeReport`: the descent loop — CutFEM solve (traction-free Γ,
  strong box-edge clamps, a checked right-box-edge traction band), nodal
  strain-energy densities sampled inside material, normal extension,
  fs-adjoint Sobolev H¹ smoothing, `v_n = w − ℓ`, one interface move
  per iteration, redistance + audit, augmented-Lagrangian volume
  multiplier NORMALIZED BY THE MEAN ENERGY SCALE (an O(1) multiplier
  against O(J) energies shrinks the structure to nothing at full
  speed — measured failure mode), scheduled nucleation; ledger rows
  with compliance, volume, ℓ, drift, and FNV snapshot hashes. The load is
  definitionally zero outside the checked `EdgeBand`; unrelated SDF cuts on
  the same edge are skipped, while a caller-supplied cut through supported load
  refuses before mutation. After every optimizer-owned advection,
  redistancing, or nucleation update, a deterministic two-cell-deep
  non-design material pad is retained around the loaded band with one cell of
  endpoint clearance. Per-iteration ledger rows and `OptimizeReport` retain
  the number of nodal assignments made by that policy.
  fs-topols constructs `IsotropicElastic` and calls fs-cutfem directly; the
  historical Nitsche and ghost coefficients are multiplied by
  `(lambda + 2*mu)/mu` so removing the fs-solid compatibility facade does not
  change the discrete operator. The solver retains its `1e-12` relative gate
  and 60,000-iteration cap.
  Reported compliance is canonical assembled-load `b^T u`, not the former
  node-mask/trapezoid proxy. With zero body force and embedded displacement
  data here, it is the exact discrete external work.
- `optimize::material_volume`: certified cut-quadrature area of
  `{φ < 0}`.

## Invariants

1. Enclosure containment holds on the DISCRETE moving geometry
   (bilinear corner-hull law, dense-sample battery).
2. WENO order on smooth advection: measured order > 1.6 at the finest
   pair, absolute error < 5e-4 after half a revolution (tls-001;
   design order 3 with RK3 dominates asymptotically).
3. Redistancing moves the sampled zero set < 0.2h, restores
   |∇φ| = 1 to mean deviation < 0.05, and is idempotent (< 0.05h on
   repeat) — tls-002.
4. Extension is constant along interface normals (radial fixture
   deviation < 0.15 across a 0.13-wide annulus) — tls-003.
5. Volume is conserved under rigid rotation through five
   advect+redistance cycles within 1%, with per-cycle drift audits
   < 0.25h (tls-004 — the drift POLICY evidence).
6. Sensitivities are numerically gated: topological-derivative
   prediction vs a real punched hole and shape-velocity prediction vs
   a real uniform boundary motion both land in [0.25, 4]× with the
   right SIGN (tls-005) — the bead's adjoint-vs-FD gate.
7. The cantilever descent converges to the volume target within 0.05,
   stabilizes (tail variation < 20%), fires nucleation with positive
   predicted gain producing a GENUINE interior hole (flood-fill
   verified), keeps every redistancing drift < 0.5h, and BEATS the
   trivial uniform-band design at equal volume (tls-006). The actual final
   optimizer field also crosses the L4→L2 honesty boundary: cell-center
   sampling into fs-topo's public `VoxelField` must yield exact cubical
   `b0 = 1`, `b1 >= 1`, and `b2 = 0`; the H0 persistence filtration's
   face-connected alive-component count is retained separately for
   localization, and an independent repeat reproduces the complete audit
   identity bitwise.
8. Determinism (P2): two descent runs produce bitwise-identical FNV
   snapshot sequences.
9. Typed right-edge load support (G0/G3): aligned or non-aligned checked bands
   wholly inside material solve deterministically; an SDF crossing through
   the supported band refuses. Non-finite/non-positive load magnitudes, band
   half-widths outside `[0, 0.5]`, and material settings outside the canonical
   finite coercive plane-strain regime refuse before the level set mutates.
10. LOAD SUPPORT IS NON-DESIGN: once the initial canonical solve admits the
    caller geometry, every optimizer-owned geometry update restores a
    two-cell-deep right-edge material pad extending one cell past both support
    endpoints. The complete supported segment remains strictly material on the
    next solve. The policy scans only the last three lattice columns, records
    every changed-node count, and replays exactly; it never weakens canonical
    cut-edge refusal or fills unrelated columns.

## Error model

fs-cutfem's `CutFemError` propagates unchanged from physics solves.
`InvalidElasticityInput` names
invalid cantilever load/band data, invalid optimizer material data, or the
canonical typed-support refusal when the initial SDF cuts a loaded segment.
Optimizer-owned erosion of an already-admitted load pad is repaired by the
explicit bounded non-design policy and recorded, not passed to the solver or
silently accepted. The
inherited certified plane-strain bound is `(lambda + 2*mu)/mu <= 4`,
equivalent for the isotropic card to `nu <= 1/3`; larger values refuse rather
than entering the unsupported near-incompressible regime. Structured asserts
(panics) guard programmer contracts: lattice/grid alignment, nodal-array
lengths, non-uniform grids where v1 requires uniform. Audits never silently
degrade — drift and gradient deviations are returned, not clamped.

## Determinism class

Bit-deterministic across runs on a fixed platform: BTree/index-ordered
sweeps, deterministic solvers underneath, snapshot hashes asserted
bitwise in the battery. Cross-ISA goldens not yet recorded.

## Cancellation behavior

Bounded synchronous loops (fixed iteration counts, FIM sweep caps,
CFL-bounded step counts). Load-pad retention adds at most three lattice-column
visits per geometry publication. Chunked Cx polling belongs to the fs-exec
driver (L4 consumes fs-cutfem's canonical L3 kernel directly).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None. §9.5's [S] core (evolution, redistancing, coupling) and the [F]
nucleation operator ship together as a standalone L4 crate per the
crate-granular gating rule.

## Conformance tests

`tests/battery.rs`: tls-001 WENO order; tls-002 FIM audits +
idempotency; tls-003 normal extension; tls-004 volume conservation +
drift policy; tls-005 numerical sensitivity gates (DT vs punched
hole, shape velocity vs FD); tls-006 the cantilever descent (volume,
stabilization, nucleation with interior-hole flood-fill proof,
determinism, recorded/replayed load-pad intervention, final typed-support
validity, beats-trivial-at-equal-volume). tls-006 additionally sends that
same final field—not a synthetic replacement—through fs-topo's production H0
persistence and exact cubical Betti oracle after a one-cell extrusion. A
versioned fs-obs receipt binds both crate versions, the optimizer's final
snapshot, every sampled level-set bit, grid metadata, threshold, Betti tuple,
and persistent material-component count. A repeat must reproduce the receipt
bitwise; the receipt explicitly records whether face-connected persistence
localization agrees with the closed-cube Betti component authority. Disclosed
filled-void and two-island mutations must respectively erase the tunnel and
produce exactly two components before the passing `ConformanceCase` can be
emitted. Unit tests: enclosure
containment including degenerate lattice-aligned boxes, planar extension,
typed aligned-band success/replay, supported-cut refusal, deterministic strict
load-pad retention without unrelated-column fill, and invalid fixture/material
no-mutation refusals.

## No-claim boundaries

- SIMP cross-validation on a SHARED fixture: fs-topopt's SIMP is
  3D-tet-based while this crate is 2D plane strain — a
  dimension-matched pairing (2D density pipeline or 3D level set)
  closes it; recorded follow-up. The shipped quality gate is
  beats-trivial-at-equal-volume plus the sensitivity gates.
- Octree narrow bands (uniform-lattice bands ship; fs-rep-sdf tile
  storage wiring is the consumer path per the plan's shared-substrate
  design).
- 3D, multiple load cases, stress constraints, compliant mechanisms.
- Certified clipping of an SDF-cut loaded edge segment is not claimed; typed
  support refuses that case. The optimizer's non-design pad prevents its own
  evolution from creating that unsupported case; it is not a cut-edge
  quadrature implementation. The two-point rule is exact for the shipped
  constant traction times Q1 edge shapes, but no quadrature-error claim is made
  for arbitrary traction callbacks.
- fs-adjoint Hadamard boundary-form velocities (the volumetric
  energy-density form ships; the Hadamard trace form composes when
  the FEEC trace machinery lands).
- Perimeter/curvature regularization beyond Sobolev smoothing;
  velocity CFL coupling to BAND MAINTENANCE beyond the shipped
  one-move-per-iteration policy.
- The vessel flagship's lip-channel fixture (§15.3) — the cantilever
  classic ships; the flagship composes downstream.
- The output audit is exact for the sampled cubical solid only. Cell-center
  sampling can miss sub-cell features of the continuous bilinear `GridSdf`, so
  it is a resolution-scoped cross-layer audit rather than a continuous-topology
  theorem. fs-topo's shipped persistence pairing is H0; the nucleated planar
  hole is checked by exact cubical `b1`, not laundered into an H1 persistence
  claim.
