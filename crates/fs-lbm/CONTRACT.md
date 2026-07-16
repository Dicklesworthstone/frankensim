# CONTRACT: fs-lbm

Lattice Boltzmann core with a D2Q9 BGK path, a tile-major D3Q19 BGK/Guo path,
link-wise D3Q19 wall/open-boundary topology, the lattice-scaling assistant,
and dense-grid extension scaffolding for vector forcing, local rheology,
thermal double-population fixtures, and a free-surface mass-ledger prototype.

## Purpose and layer

Layer L3 (FLUX). Depends on `fs-evidence` (the `Color` for the Evidence-typed
scaling plan) and deterministic `fs-math` primitives. Pure, deterministic
(fixed tile/cell/link order).

## Public types and semantics

- D2Q9 constants: `Q` (9), `CS2` (`1/3`), `MACH_LIMIT` (`0.3`).
- `equilibrium(rho, ux, uy) -> [f64; 9]` — the D2Q9 equilibrium distribution
  (recovers `Σf = ρ`, `Σeᵢfᵢ = ρu`).
- `Lbm::channel(nx, ny, tau, gx)` — a body-force-driven channel (periodic x,
  bounce-back y-walls). `step`/`run` (collide + Guo forcing + stream +
  bounce-back); `density`, `velocity` (Guo-corrected), `total_mass`,
  `viscosity` (`ν = (τ−½)/3`), `x_velocity_profile`.
- `core2::Grid` / re-exported `Grid` — a general dense D2Q9 grid with cell
  flags, vector gravity, per-cell relaxation time, per-cell external force,
  periodicity flags, deterministic collide/stream steps, and gas/wall
  boundary bounce handling for the plain non-free-surface step.
  `core2::MomentumExchange2`, `stream_from_with_wall_momentum`, and
  `step_with_wall_momentum` add an opt-in raw lattice-impulse receipt for an
  explicitly selected subset of stationary wall cells without changing the
  legacy step path. `core2::MovingWallMomentumExchange2` plus the paired
  moving-wall stream/step methods apply per-wall-cell halfway-bounce velocity
  corrections and return boundary-relative impulse, torque, work, resolved
  population momentum, and moving-mass transfer terms for an explicitly
  selected wall subset. `core2::VelocityPressureX2` plus the paired Grid step
  methods provide a low-Mach regularized velocity inlet at x-min and density
  outlet at x-max with periodic y closure; density/velocity and non-equilibrium
  stress are extrapolated from the respective first-interior columns.
- `rheology::Rheology`, `rheology::update_tau`, and
  `rheology::channel_flow` — local apparent-viscosity laws and explicit
  τ updates with floor/cap counts for cells outside the representable
  relaxation window.
- `thermal::ThermalLbm` and `thermal::gbeta_for_rayleigh` — D2Q9 flow plus
  D2Q5 temperature populations for Rayleigh-Bénard-style onset fixtures with
  fixed-temperature wall rows.
- `freesurface::FreeSurface`, `ContactModel`, `dam_break`, and `surge_front`
  — dense-grid VOF-style mass tracking with interface/gas/fluid conversion
  bookkeeping, conservative carry redistribution, contact-model bracketing,
  and qualitative dam-break / jet-fragment fixtures.
- `d3q19::{Duct, equilibrium3, duct_analytic}` — the frozen D3Q19 BGK + Guo
  body-force core on aligned 4x4x4 SoA tiles, with stationary halfway
  bounce-back x/y walls and periodic z for the rectangular-duct fixture. Its
  axial-z BGK collision dispatches once to a bit-neutral NEON/AVX2 tile capsule
  or the scalar twin. Its separate pull-stream dispatcher moves four-cell x
  rows through a bit-neutral NEON/AVX2 capsule or the scalar stencil; SIMD
  lanes are independent cells, and streaming performs no floating arithmetic.
- `d3q19_bgk_simd_tier` / `D3q19BgkSimdTier` — truthful operation-specific
  receipt for the Duct collision kernel. AVX-512 hosts report AVX2 because this
  capsule deliberately uses four-lane AVX2; x86 without AVX2+FMA reports scalar.
- `d3q19_stream_simd_tier` / `D3q19StreamSimdTier` — separate receipt for the
  Duct pull-stream kernel. AArch64 reports NEON; x86 requires only AVX2 for this
  pure-move capsule and otherwise reports scalar.
- `d3q19::{CollisionModel3, collide_cell3, CollisionError3}` — the checked
  per-cell collision authority used by `BoundaryGrid3` and public cell calls.
  `Bgk { tau }` retains the frozen general-force projection; the Duct tile
  scalar/SIMD twin is separately bit-locked to its test-only axial oracle. The unforced
  `CentralMoment` reference rung relaxes a full-rank centered D3Q19 monomial
  basis with independent second/higher-order rates. `ReducedCumulant` applies
  the nonlinear `C220`/`C202`/`C022` corrections to the three independent
  fourth-order rows represented by that basis. All paths reject invalid
  parameters or state before publication.
- `D3Q19_MOMENT_COLLISION_SEMANTICS_VERSION` — independent version for the
  optional moment-space basis, rate grouping, `mul_add` chains, cumulant
  projection, and deterministic back-transform. It does not refreeze the
  established BGK grid golden.
- `Face3`, `FaceBoundary3`, and `BoundarySpec3` — six-face axis-aligned
  boundary declarations: paired periodic faces, tangential moving/stationary
  halfway walls, regularized velocity faces, and isothermal pressure/density
  faces. Open faces currently require zero body force.
- `LinkMaskTile3` and `BoundaryLink3` — aligned per-tile D3Q19 effective-wall
  masks and canonical `(tile, lane, direction)` enumeration. Separate open
  masks distinguish populations reconstructed after streaming.
- `BoundaryGrid3` — tile-major D3Q19 grid whose compiled masks drive the
  runtime pull stencil. The legacy `new(..., tau, ...)` constructor remains
  BGK, while `with_collision_model` selects an explicit checked model for the
  whole grid. `voxelize_sdf` samples a scalar field at cell centers, commits
  solid occupancy/masks atomically, and then locks topology. Solid and planar
  walls use link-wise halfway bounce-back. Moving walls add the standard
  `2 w_i rho (c_i dot u_wall) / c_s^2` incoming-link correction.
  Velocity/pressure face-interior cells use second-order Hermite regularized
  non-equilibrium-stress reconstruction from the first interior cell.
- `plan_scaling(reynolds, char_length_lu, u_lattice) -> ScalingPlan { tau,
  viscosity, u_lattice, mach, tau_margin, stable }` — the lattice-scaling
  assistant. `ScalingPlan::color()` (verified when comfortably stable, else
  estimated). Panics on non-positive Reynolds / length.
- `poiseuille_analytic(gx, viscosity, ny, y)` — the analytic reference profile.

## Invariants

- The equilibrium recovers its density + momentum moments exactly.
- Unforced collision preserves density and all three momentum components to
  roundoff. The boundary grid delegates to the public cell authority; the Duct
  tile scalar/SIMD twin is differential-tested against its frozen axial oracle.
- Equal-rate central-moment collision reduces to BGK within deterministic
  transform/solve roundoff. Split higher-order relaxation changes
  nonequilibrium modes without relaxing degree-zero/one invariants.
- Reduced-cumulant collision preserves degree-zero/one invariants, fixes the
  discrete equilibrium to solve roundoff, and is covariant under Cartesian
  axis permutation when rates are equal within each represented order. Its
  nonlinear fourth-order relaxation does not reduce to BGK at equal rates.
- MASS is conserved by a closed-domain step (collision, forcing, streaming,
  and bounce-back all conserve mass). Prescribed velocity/pressure faces are
  open-system flux boundaries and do not claim global mass conservation.
- Duct SIMD streaming is a pure population-bit permutation. Its x/y wall
  crossings read the opposite population at the destination cell, including
  corners and simultaneous z wraps; every non-wall z source is periodic.
- Steady Poiseuille channel flow matches the analytic parabola to a few percent
  (halfway bounce-back resolves the quadratic profile).
- `plan_scaling` derives `τ = 3ν + ½`, flags `stable` iff `τ > ½` AND
  `Mach < MACH_LIMIT`.
- General dense-grid constructors reject zero dimensions and nonphysical
  relaxation times before arithmetic can produce NaNs.
- Gas cells do not act as fluid population sources in the plain dense-grid
  stream step; absent gas-side populations bounce at the fluid boundary until
  explicit free-surface bookkeeping lands.
- D2Q9 wall momentum exchange sums `-2 f_post c_q` in fixed row-major/link
  order for selected stationary halfway-bounce-back links. Gas and exterior
  bounces are excluded; an isolated wall in equilibrium fluid has zero net
  impulse, and selecting one obstacle cannot silently include another.
- D2Q9 moving-wall streaming adds
  `2 w_q rho_post (c_q dot u_wall) / c_s^2` to each incoming wall link. Its
  selected-wall impulse uses the boundary-relative Wen-style exchange
  `(c_out - u_wall) f_out - (c_in - u_wall) f_in`. The receipt independently
  retains the resolved fluid-population impulse and `u_wall (f_in - f_out)`;
  their linkwise balance closes in fixed row-major/direction order. Torque uses
  destination-local halfway-link midpoints about the caller's finite origin.
  The force convention is from Wen et al., *Galilean Invariant Fluid-Solid
  Interfacial Dynamics in Lattice Boltzmann Simulations* (2014,
  <https://arxiv.org/abs/1303.0625>).
- D2Q9 regularized x faces impose the declared inlet velocity and outlet
  density to roundoff while copying the complementary moment and independently
  measured non-equilibrium stress from the first interior column. The measured
  and unmeasured open-step paths are bit-identical.
- Rheology laws reject non-finite or non-positive physical parameters, and
  every update reports floor/cap counts when viscosity leaves the representable
  τ window.
- Thermal wall populations encode the declared wall temperatures, so the
  public `temperature` query is consistent on wall and fluid rows.
- Free-surface steps conserve the tracked ledger mass (fluid `Σf` plus
  interface mass plus carry) to the test tolerance, and gas/interface/fluid
  conversions are counted rather than hidden.
- D3Q19 wall and open masks are disjoint and compiled in tile/lane/direction
  order. Where an open face meets a wall, directions crossing both faces are
  stationary wall-owned while pure open-face directions are reconstructed;
  where moving exterior walls disagree at an edge/corner, the whole cell is
  stationary so diagonal corrections cannot create an unpaired mass source.
- Voxelized topology is initialization-only and two-phase: failed/non-finite,
  all-solid, or open-neighbor-invalid proposals leave occupancy and masks
  unchanged. A committed topology cannot be silently changed after
  perturbation or stepping.
- Stationary planar and voxel bounce-back conserve mass to roundoff. Tangential
  moving-wall corrections sum to zero over the represented link set.
- On face-interior cells, regularized velocity/pressure reconstruction
  preserves the declared zeroth and first moments and copies the independently
  measured non-equilibrium second moment from the first interior cell. Mixed
  wall/open rim cells preserve their per-link wall ownership instead of
  claiming exact target moments.

## Error model

Most operations are total over physically admissible inputs. Constructors and
parameter helpers panic on nonsensical requests: zero dimensions, non-finite
forces/relaxation times, non-positive viscosities/rheology indices, non-positive
Rayleigh height, or non-positive Reynolds/length in the scaling assistant.
D2Q9 momentum measurement additionally requires a full-grid boolean mask whose
selected entries are all `Cell::Wall`; this is checked before a measured step
can mutate populations.
D2Q9 moving-wall calls additionally require one finite velocity per grid cell,
zero velocity on every non-wall cell, speed squared below 0.03, and a finite
moment origin. Every moving wall adjacent to fluid requires positive finite
post-collision density, and every outgoing wall-link population must be finite.
Request fields and post-collision state are admitted before streaming mutates
the grid.
D2Q9 regularized x flow requires at least three columns, non-periodic x,
periodic y, fluid face/first-interior columns, zero gravity/external forcing,
a positive finite outlet density, and a finite inlet speed squared below 0.03.
Every condition is checked before collision.
D3Q19 boundary construction additionally rejects non-4-multiple dimensions,
tile-count overflow, invalid collision parameters, moment-space collision with
nonzero force, unpaired periodic faces, non-tangential/non-finite or
outside-low-Mach wall velocities, non-finite or outside-low-Mach inlet
velocities, non-positive pressure density, open faces on more than one axis,
and body force combined with the current regularized open-face closure.
Voxelization rejects non-finite samples, an all-solid domain, an obstructed
first-interior open-face layer, and topology mutation after initialization.
Boundary-grid perturbation rejects non-finite amplitudes or magnitudes at least
one before changing populations or locking topology.
`collide_cell3` returns typed errors rather than publishing non-finite or
non-positive cell states. Both moment-space rungs additionally refuse rates
outside `(0, 2)`, nonzero body forcing, or a numerically rank-deficient moment
transform.

## Determinism class

Fully deterministic: fixed cell iteration order for D2Q9 and fixed logical
tile/cell/direction identities plus face/z/y/x reconstruction order for D3Q19;
no RNG. Duct collision retains the frozen direction-ascending reductions and
scalar operation tree in fixed lane blocks. Duct pull-streaming uses a pinned
direction/z/y/x-tile row schedule, but only moves bits into slot-exclusive
destinations, so schedule changes cannot create arithmetic reassociation.

## Cancellation behavior

None here (a step is synchronous); polling at tile boundaries under `Cx` is the
production kernel's concern.

## Unsafe boundary

Four registered leaf capsules opt out of the workspace default at module scope:
the collision pair under `src/d3q19/simd/{neon,x86}/mod.rs` and the pull-stream
pair under `src/d3q19/simd/stream/{neon,x86}/mod.rs`. Each stays below the
capsule line cap with an adjacent `SAFETY.md`. Safe façades validate tile shapes,
Rust borrows prove aliasing and lifetimes, and x86 one-shot selectors admit
private thunks only after their operation-specific feature checks. Miri selects
safe scalar twins.

## Feature flags

None.

## Conformance tests

`tests/lbm.rs` covers the v0 core: equilibrium moments; mass conservation;
Poiseuille flow matches the analytic parabola (symmetric, centered); the
scaling assistant derives τ + flags stability + colors the plan; it rejects a
high-Mach plan and nonsense inputs; determinism; exact one-link wall-impulse
sign/magnitude, obstacle selection, equilibrium cancellation, replay
determinism, and pre-step mask refusal for D2Q9 momentum exchange; an
independently enumerable moving-wall link that pins bounce correction,
boundary-relative force, torque, work, and moving-mass balance; exact
zero-velocity compatibility with the stationary API; moving-step bit replay;
and fail-closed moving-field admission. It also covers regularized x-face
moment/stress reconstruction, measured-path bit equivalence, and pre-step
topology/forcing refusal.

`tests/extensions.rs` covers the current extension scaffolding: power-law and
Newtonian-limit channel profiles, Carreau plateaus, Rayleigh-Bénard onset
bracketing and Nusselt heat transport, gas-neighbor streaming behavior, thermal
wall-temperature queries, invalid-parameter rejection, free-surface mass-ledger
conservation, qualitative dam-break front advance, rotation equivariance,
contact-model bracketing, and qualitative jet fragmentation.

`tests/extensions.rs` (tfz.19): lbm-101 power-law/Carreau profiles vs
analytic (0.12% with τ floor+cap ledger); lbm-102 Rayleigh–Bénard onset
bracket (decay Ra=1200 / growth Ra=2500, Nu>1); lbm-104 STRICT free-
surface mass ledger (5e-14 over 600 dam-break steps, conversions
counted); lbm-105 dam-break front envelope (coarse honesty band);
lbm-106 G3 rotation equivariance (3e-14); lbm-107 contact-model
bracket band + Plateau–Rayleigh jet fragmentation with strict ledger;
lbm-108 level-jump refinement (Poiseuille through the interface +
shear-wave decay-rate transmission, first-order labels).

`tests/d3q19_battery.rs` covers the frozen D3Q19 core: exact integer lattice
moments/opposites, equilibrium moments, mass conservation, analytic
rectangular-duct flow, replay determinism, the registered core golden, shared
kernel bit-equivalence, BGK/central-moment equal-rate equivalence, split-rate
conservation, reduced-cumulant equilibrium/nonlinearity/axis-covariance laws,
and fail-closed inputs.

The D3Q19 SIMD module adds G0 seeded differential batteries for both operations.
Collision compares every lane to its scalar tile twin and directly locks that
twin to the frozen axial cell authority. Pull-streaming compares every
population bit on a single-tile fixture and an asymmetric 2x3x4-tile fixture;
independent route anchors cover inter-tile x/y/z moves, periodic z, corner
bounce, and wall precedence over simultaneous z wrap. Both batteries retain
first-divergence JSONL receipts, hostile dispatch-selector cases, and one-shot
table identity. The existing Duct golden remains the end-to-end bit-surface
gate.

`tests/d3q19_boundaries.rs` (bead 40p2) covers all six hand-enumerated planar
link masks, aligned deterministic mask ordering, the exact 18 links around one
voxel, atomic immutable SDF topology and rejection paths, axis-generic
regularized face density/velocity plus independent inlet/outlet stress
reconstruction, grid-level collision selection and equal-rate BGK/central
agreement, split-rate cumulant conservation, open/wall and moving/open rim
ownership, stationary planar and voxel mass conservation, exact one-step
moving-lid momentum oracles, qualitative primary cavity circulation,
pressure-driven Poiseuille shape, and a boundary replay-hash candidate. Ignored
release fixtures carry the full 10,000-step leak and 32x32 full-rim 3%
pressure-Poiseuille gates.

`tests/cylinder_re100.rs` (frankensim-wghy) keeps the inexpensive evidence
machinery in the default suite: a checked linear-detrended, Hann-window
`fs-fft` lift-frequency estimator with peak-prominence, exact-bin,
affine-trend, replay, flat/non-finite/ambiguous/Nyquist, and
resolution-refusal oracles; and an exactly projected stair-step cylinder mask
pinned by occupied cells, centroid, reflection symmetry, fluid-wall link count,
and FNV-1a fingerprint. The ignored release-scale `lbm-109` fixture runs Re=100
at diameter 10, nominal inlet speed 0.1, 32D streamwise extent, 12D and 16D
periodic lateral spans, 8,192 warm-up steps, and 32,768 retained force samples.
It normalizes drag by measured mean inlet density, requires raw-domain and
split-window Cd guards before admitting a bounded empirical two-width
zero-blockage Cd intercept, and requires the full-window and both half-window
intercepts to stay in `[1.25, 1.45]` with at most 2% split drift. It does
not extrapolate spectral bins: both widths and both half-windows must put St in
`[0.155, 0.175]`, split windows must agree within `0.0062`, and the reported 16D
St must agree with the 12D sensitivity run within `0.01`. The deterministic
symmetry-breaking seed is disclosed as offset `(10, 6)` at transverse speed
`1e-4`. Roshko NACA TR-1191, Posdziech-Grundmann 2007, Behr et al. 1995, and
Maskell ARC R&M 3400 are cited at the executable gate; no publisher artifact is
redistributed.

## No-claim boundaries

- D3Q19 grids remain BGK + Guo on a dense set of aligned SoA tiles. The
  selectable central-moment and `ReducedCumulant` cell operators are
  deterministic `O(Q^3)` correctness references: they are unforced and have no
  performance or high-Re stability claim. D3Q27, sparse active-tile
  storage/sweeps, a production cumulant collision, normalized drag/lift
  histories, and bandwidth roofline / fs-tilelang kernels remain staged. The
  raw D2Q9 stationary-wall momentum receipt is not a normalized aerodynamic
  coefficient. Geier
  et al.'s primary derivation (doi:10.1016/j.camwa.2015.05.001) explicitly
  restricts itself to D3Q27 after identifying non-refining D3Q19 anisotropy.
  `ReducedCumulant` therefore implements only the general cumulant definitions
  for the independent D3Q19 `C220`/`C202`/`C022` projection; it is not a
  "Geier D3Q19" operator and cannot borrow the paper's high-Re evidence.
- The current SIMD slice accelerates the all-fluid Duct axial-BGK collision and
  its pure-move pull stream. `BoundaryGrid3` mixed masks/general forcing and the
  central-moment/cumulant paths remain scalar; collision and streaming are not
  fused. Therefore this slice alone makes no stream+collide GLUP/s claim. Both
  native ISA batteries, the unchanged Duct golden, and the four-quadrant
  performance lane remain required before promotion evidence is green.
- SDF voxelization is midpoint classification followed by stair-step halfway
  bounce-back. It is second-order for the represented flat, lattice-aligned
  halfway wall, not a second-order certificate against the original continuous
  curved SDF. Interpolated Bouzidi-type curved boundaries remain staged.
- D2Q9 `MomentumExchange2` is a raw stationary-wall lattice impulse over the
  caller's exact cell mask. `MovingWallMomentumExchange2` adds a static-topology
  halfway-wall velocity correction and boundary-relative force evaluation, but
  it does not move cell topology or initialize newly uncovered fluid. Neither
  receipt applies physical-unit conversion, reference-area normalization,
  curved-boundary interpolation, blockage correction, averaging, or
  shedding-frequency estimation; therefore this is not yet the Re=100 cylinder
  Cd/St validation. Boundary-relative exchange improves the force receipt's
  frame behavior; it is not by itself a proof that the full bounce-back solver
  is Galilean invariant.
- The separate `lbm-109` release fixture encodes the intended normalization,
  warm-up, detrended FFT, raw/split-window guards, primary-source envelopes,
  and empirical two-width Cd sensitivity treatment. Maskell's closed-tunnel
  law is not claimed to transfer to a periodic-y fixture; the linear-beta Cd
  intercept merely eliminates one assumed leading coefficient from two
  disclosed widths, and its magnitude is bounded. Strouhal uses the 16D result
  plus 12D sensitivity rather than a frequency-bin extrapolation. Until that
  explicitly ignored release lane is executed green on the combined batch
  snapshot, the crate still makes no completed Re=100 cylinder validation
  claim.
- `VelocityPressureX2` is a low-Mach regularized fixture boundary with periodic
  lateral closure. It is not a characteristic or non-reflecting far-field
  condition, accepts no body force, and makes no unbounded-domain or blockage
  claim. The cylinder battery must disclose its domain, resolution, warm-up,
  averaging window, and lateral-domain sensitivity separately.
- The cylinder fixture's diameter-10 stair step and `tau = 0.53` have no grid-
  convergence certificate, curved-wall accuracy claim, or seed-independent
  shedding-phase claim. Split-window and two-width agreement detect selected
  transients and domain sensitivity; they do not replace a resolution family.
- The face-generic regularized closure was selected instead of six
  face-specialized Zou-He tables because its Hermite stress projection has an
  independent second-moment oracle and preserves arbitrary tangential target
  components under one implementation. This is not evidence that a correct
  Zou-He/Hecht-Harting implementation is unstable. The current closure is a
  low-Mach fixture, not a high-Re characteristic/non-reflecting boundary; it
  accepts only constant per-face targets and refuses simultaneous Guo force.
- Open/wall intersections use an explicit per-link mixed policy: wall-crossing
  links bounce stationary, pure open-normal links are regularized, and links
  crossing both faces are wall-owned. This is not a geometric corner
  reconstruction, and exact target moments are claimed only on face-interior
  cells. Moving-wall cavity evidence is qualitative until the separate G2
  benchmark bead lands.
- Solid topology is static after initialization. Moving/deforming SDF
  boundaries require a future explicit population initialization, mass delta,
  and topology-transition receipt rather than reusing `voxelize_sdf`.
- Interface and gas flags exist so the plain core can share the future data
  model, but free-surface mass/VOF bookkeeping is not implemented in
  `Grid::step`; gas-side pulls currently bounce rather than reconstructing
  missing free-surface populations.
- Thermal and rheology fixtures are dense-grid correctness scaffolding, not
  validated LES, cumulant, sparse-tile, or production multiphase solvers.
- The free-surface implementation is a dense prototype with ledger and
  metamorphic gates. It does not yet claim quantitative dam-break agreement,
  contact-angle calibration, production wetting physics, or validated
  surface-tension breakup rates.
- The scaling assistant covers the `τ`/`ν`/`Mach` core; consuming fs-regime's
  dimensionless groups and emitting a full `dx`/`dt` unit conversion with
  Evidence provenance is the fuller deliverable.
- Grid refinement is the TWO-LEVEL 1:2 channel coupling with Dupuis–Chopard
  non-equilibrium rescaling and a FIRST-ORDER interface handoff: measured
  2.5% steady Poiseuille deviation at the level jump and 5.6% extra shear-
  wave decay rate — honesty-labeled in lbm-108; the post-collision
  (Filippova–Hänel-style) second-order transfer, general octree topologies,
  and dwr-adaptivity-driven refinement signals are recorded successors.
- Contact-line physics is MODEL-BRACKETED (neutral vs wetting fill ghosts,
  lbm-107 reports the sensitivity band), never pretended-certain — the
  §15.3 caveat is a design decision here, not an omission.
- ADJOINT HONESTY (plan §8.7 [M]): free-surface LBM adjoints are NOT
  promised — cell-conversion events make the map non-differentiable;
  gradients for free-surface objectives go through surrogate or
  gradient-free lanes. The model card is this paragraph.
- Pouring scenarios: tilt schedules enter as the rotating gravity vector
  (lbm-106 pins 90-degree equivariance at 3e-14); full fs-scenario moving-
  frame integration and Plateau–Rayleigh breakup SCORING (beyond the
  qualitative fragment gate) are staged with the vessel flagship.
