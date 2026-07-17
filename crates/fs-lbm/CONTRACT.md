# CONTRACT: fs-lbm

Lattice Boltzmann core with a D2Q9 BGK path, a tile-major D3Q19 BGK/Guo path,
link-wise D3Q19 wall/open-boundary topology, the lattice-scaling assistant,
and dense-grid extension scaffolding for vector forcing, local rheology,
thermal double-population fixtures, and a free-surface mass-ledger prototype.

## Purpose and layer

Layer L3 (FLUX). Depends on `fs-evidence` (the `Color` for the Evidence-typed
scaling plan), deterministic `fs-math` primitives, `fs-matdb` for immutable
interface-system queries and usage receipts, and `fs-qty` for runtime
dimension checks. Pure, deterministic (fixed tile/cell/link order).

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
  selected wall subset. `core2::CurvedMovingWallLink2` plus the paired curved
  moving-wall stream/step methods consume an exhaustive per-fluid-cell,
  per-direction wall-link map and apply the linear Bouzidi-Firdaouss-Lallemand
  off-lattice interpolation. Link-local velocity and the declared wall
  fraction determine both the incoming population and the exact intersection
  used for torque. `core2::WallTopologyTransition2` and
  `Grid::transition_wall_topology` atomically replace a fluid/wall cell mask,
  initialize newly uncovered cells from unique surviving one-ring donors, and
  receipt the exact active-population mass/momentum removal and insertion.
  `core2::PartialSaturationCell2`, `core2::PartialSaturationExchange2`, and
  the paired Grid collision/step methods add a force-free dense-D2Q9
  Noble-Torczynski partially saturated cell rung. Caller-supplied solid
  fractions and velocities blend fluid and non-equilibrium-bounce-back solid
  collisions while the receipt retains aggregate equal-and-opposite
  fluid/solid impulse, angular impulse, solid work, and coupling-weight sum.
  `core2::UnderResolvedCylinderPair2` and
  `core2::UnderResolvedLubricationExchange2` add a separate pair-local,
  per-unit-depth normal-resistance correction for smooth circular cylinders
  below a caller-declared cutoff. The immutable descriptor owns checked gap,
  reduced radius, viscosity, common normal, line-of-action point, and surface
  velocities; its receipt returns exact equal-and-opposite solid impulse and
  angular impulse, dissipative pair work, the effective gap, and the applied
  resistance increment without mutating the lattice or a rigid-body state.
  `core2::VelocityPressureX2` plus the paired Grid step methods provide a
  low-Mach regularized velocity inlet at x-min and density outlet at x-max with
  periodic y closure; density/velocity and non-equilibrium stress are
  extrapolated from the respective first-interior columns.
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
- `freesurface::query_dynamic_wetting2`, `DynamicWettingAnswer2`,
  `ContactLineRegime2`, and `DynamicWettingError2` — a read-only adapter from
  an ordered, history-bearing `InterfaceSystemCard` to a receipt-bearing
  receding/advancing contact-angle bracket. Both speed inputs are `QtyAny`
  values carrying coherent-SI velocity dimensions. It queries both endpoints
  at the absolute contact-line speed, selects by signed speed outside a
  symmetric pinning deadband, and returns no selected angle while pinned. The
  answer retains both `MaterialAnswer` receipts and the interface-system
  content identity.
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

### Sparse active-tile infrastructure (`d3q19::sparse`, bead sjro)

- `morton3(tx, ty, tz)` / `demorton3(key)` — 3×21-bit Morton (Z-order) tile
  keys; `morton3` panics past the 21-bit range (a compile-shape bound, not a
  runtime input).
- `SparseGrid3::new(nx, ny, nz, tau, force)` — an EMPTY sparse D3Q19 grid
  over a domain of positive 4-multiples; state exists only for active tiles.
  `activate_tiles(&[(tx,ty,tz)])` is idempotent, atomic (refuses whole-batch
  on any out-of-domain tile), preserves existing tile state bit-exactly, and
  re-sorts the active list into ascending Morton order — iteration, sweep,
  and reduction order depend only on the active SET, never on activation
  history.
- `step_serial()` / `step_pooled(&runner, &gate)` — one sweep step as two
  passes (collide from the shared per-cell `collide_cell3` kernel, then
  pull-stream) over three retained buffers. The pull-stream pass writes a
  private transactional destination; only a completely drained successful
  pass swaps that destination into the published state. A link whose source
  lies in an inactive tile or outside the domain bounces back (inactive space
  is solid wall in this increment). `step_pooled` runs both passes through any
  `fs_exec::KernelRunner` and is bitwise-identical to `step_serial` for
  every worker count: kernel tiles own disjoint Morton-contiguous slot
  groups sized by a worker-count-independent constant, and all inputs are
  read-only within a pass.
- `perturb(seed, amplitude)` — splitmix64 hash-seeded density perturbation
  (the hash IS the seed schedule); `state_bits()` — canonical-order exact
  bits for golden hashing; `total_mass()`, `cell_macros(slot, lane)`,
  `active_tiles()`, `allocated_state_bytes()` (exactly three population
  buffers per active tile), `state_bytes_per_tile()`.
- `SparseError3` — typed refusals: `Dims`, `TileOutOfDomain`, `Collision`
  (Morton key + lane + the underlying `CollisionError3`), `Cancelled`
  (pre-step state intact; re-issue is deterministic), `Pool`.

### 3-D free surface on sparse tiles (`d3q19::freesurface3`, bead sxnm)

- `FreeSurface3::new(nx, ny, nz, tau, gravity, sigma, contact, fluid)` /
  `with_solid_tiles(..., solid_tiles, fluid)` — the 2-D Körner VOF
  machinery ported onto the sparse grid: per-lane `Cell3`
  (fluid/interface/gas) and tracked mass keyed by MORTON KEY
  (activation-reshuffle-proof), interface inserted between fluid and gas,
  closure asserted. Solid tiles are permanent walls that advance-driven
  activation refuses to open (fluid declared inside one panics).
- `step()` — refresh smoothed fill, collide wet lanes via the shared
  `collide_cell3`, pairwise-antisymmetric mass exchange, pull-stream with
  gas reconstruction at the `σκ/cs²` reference density, conversion
  cascade with the conservation carry. ACTIVATION INVARIANT: every newly
  created interface cell immediately activates every inactive in-domain
  tile in its 18-neighborhood (activate-on-fluid-arrival deadlocks: an
  interface backed by the activation frontier pins as fluid-against-wall).
- Probes: `fill`, `ledger_mass` (fluid Σf + interface m + carry),
  `region_mass(pred)` + `carry()` (ledger partitioning), `cell`,
  `wet_extent`, `surge_front_x(z_below)`, `conversions()`,
  `worst_step_violation()` (per-step ledger violation, logged not
  assumed), `retire_gas_tiles()` (ledger-neutral memory maintenance).
- `FreeSurfaceError3` typed refusals: `Grid`, `Collision` (Morton key +
  lane + source), `NoFluid`.
- Goldens: `fs-lbm:d3q19-pour` frozen four-quadrant (see
  golden-couplings.json). Batteries: G0 ledger algebra, G2 Martin-Moyce
  envelope (reference data in `data/reference/martin-moyce-1952.jsonl`,
  report-only curve comparison), G5 bitwise replay.
- No-claims: conservation and determinism only — no throughput claim
  (712t), no turbulence or contact-angle physics beyond the 2-D
  neutral/wetting ghost, quantitative Martin-Moyce central band is
  fine-lattice scope, solid geometry is tile-granular.

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
- D2Q9 curved moving-wall streaming uses the linear BFL rule with fraction
  `r = distance(fluid center, wall intersection) / link length`. Writing
  `delta = 2 w_q rho_post (c_q dot u_wall) / c_s^2`, the short-link branch is
  `f_in = 2r f_out + (1-2r) f_far + delta` for `r < 1/2`; the long-link branch
  is `f_in = (f_out + delta)/(2r) + (2r-1) f_opposite/(2r)` for `r > 1/2`.
  At `r = 1/2` the established moving-halfway operation tree is used exactly.
  The complete geometry map is deterministic and exhaustive: every actual
  fluid-wall pull link has one entry and every other slot is empty. Short links
  require the next lattice node away from the wall to be fluid or interface.
  Receipt torque uses the declared intersection `x_fluid - r c_q`, while its
  impulse/mass/work balance retains the existing boundary-relative convention.
  The interpolation rule follows Bouzidi, Firdaouss, and Lallemand (2001,
  <https://doi.org/10.1063/1.1399290>).
- D2Q9 partial saturation uses the Noble-Torczynski non-equilibrium-bounce-back
  variant with lattice time step one. For solid fraction `epsilon_s`,
  `B = epsilon_s (tau - 1/2) / ((1 - epsilon_s) + (tau - 1/2))`, and
  `f_post = f + (1-B) Omega_f + B Omega_s`, where
  `Omega_s[q] = f[opp(q)] - f[q] + feq[q](rho,u_s) - feq[opp(q)](rho,u)`.
  The solid-coupling correction conserves population mass to roundoff; its
  fluid momentum and angular-momentum changes are equal and opposite to the
  solid receipt. Cells and directions reduce in fixed row-major order. The
  construction follows Noble and Torczynski (1998,
  <https://doi.org/10.1142/S0129183198001084>).
- D2Q9 under-resolved cylinder lubrication is an **Estimated** local subgrid
  closure. With unit normal `n` directed from solid 1 to solid 2, signed gap
  rate `s = (u_2 - u_1) dot n`, reduced radius
  `R* = (1/R_1 + 1/R_2)^-1`, supplied gap `h`, positive floor `h_min`, and
  cutoff `h_c`, it uses `h_eff = max(h, h_min)` and
  `lambda = 3 pi sqrt(2) mu (R*)^(3/2)` times
  `(h_eff^(-3/2) - h_c^(-3/2))` when `h < h_c`; otherwise `lambda = 0`.
  The one-step correction is `J_1 = lambda s n`, `J_2 = -J_1`, with pair
  work `-lambda s^2 <= 0`. Kromkamp et al. derive the leading equal-cylinder
  force per unit length (2005, <https://doi.org/10.1017/S0022112005003551>);
  the unequal-radius expression is the local parabolic-gap Reynolds inference,
  not a separately validated unequal-cylinder result. Subtracting the
  resistance evaluated at `h_c` follows Nguyen and Ladd's general cutoff rule
  (2002, equations 34-35,
  <https://doi.org/10.1103/PhysRevE.66.046708>); their three-dimensional sphere
  exponent and calibrated cutoff values are not transferred. Powers of
  three-halves use `x * sqrt(x)` and pair evaluation has no reduction-order
  dependence.
- D2Q9 wall-topology replacement is a two-phase publication. Fresh cells use
  an equal-weight average of unique, surviving pre-transition one-ring fluid
  populations, relaxation times, and external forces; newly covered
  populations are cleared. The receipt's fresh-minus-removed mass and momentum
  equal the active-grid change to roundoff. Donor discovery and all receipt
  reductions use fixed row-major/direction order.
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
- Dynamic wetting admission requires both endpoint claims to answer under the
  same explicit `SelectionPolicy`. Each endpoint must be dimensionless,
  finite, and strictly inside `(0, pi)` radians, with receding no greater than
  advancing. Both caller inputs and any selected curve whose abscissa is
  `contact_line_speed` must carry SI velocity dimensions. The cloned query
  point binds that axis to the absolute coherent-SI speed; the sign only
  chooses the endpoint. Exact replay reproduces both usage receipts, and the
  retained system hash moves when ordered surfaces, history, claims, or models
  move.
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
D2Q9 curved moving-wall calls require a full-grid link map that is exactly
populated on fluid/interface-to-wall pulls, fractions in `(0, 1]`, finite
link-intersection velocities below the same low-Mach envelope, and a finite
moment origin. Short-link interpolation requires an in-domain fluid/interface
far donor with a finite outgoing population; long-link interpolation requires
the finite local opposite population. Moving links require positive finite
local post-collision density. Geometry is admitted before collision, and all
required post-collision values are admitted before grid publication.
D2Q9 partial-saturation descriptors require a finite solid fraction in
`[0,1]`, exactly zero velocity at zero fraction, and finite solid velocity
inside the same low-Mach envelope. The full-grid field may assign positive
fraction only to `Cell::Fluid`; all populations must be finite, every active
cell needs positive finite density, finite momentum, and `tau > 0.5`, forcing
must be exactly zero, and the moment origin must be finite. Collision is built
in a private proposal and the caller's output/scratch is published only after
every proposed population and the aggregate receipt are finite.
D2Q9 under-resolved cylinder lubrication descriptors require a finite
nonnegative gap; positive finite cutoff, floor, reduced radius, and dynamic
viscosity; `h_min < h_c`; a finite unit normal within `128 * f64::EPSILON`;
finite surface velocities individually inside the squared-speed `< 0.03`
low-Mach envelope; and a finite line-of-action point. Evaluation requires a
finite moment origin and refuses non-finite derived resistance, gap rate,
impulse, work, offset, or angular impulse. Negative gaps belong to contact and
are refused; exactly zero gap is admitted only through the explicit positive
floor. The descriptor and evaluation mutate no state, so refusal is atomic.
D2Q9 wall-topology replacement currently requires a fluid/wall-only domain, a
full-grid target mask that leaves at least one fluid cell, positive finite
covered/donor population mass, finite donor populations/momentum/external
force, and donor relaxation times above 0.5. Every fresh cell needs at least one
unique surviving one-ring fluid donor. The complete proposal and receipt are
validated before flags, populations, relaxation times, or forces are
published.
D2Q9 regularized x flow requires at least three columns, non-periodic x,
periodic y, fluid face/first-interior columns, zero gravity/external forcing,
a positive finite outlet density, and a finite inlet speed squared below 0.03.
Every condition is checked before collision.
Dynamic wetting queries return typed errors when either `QtyAny` input is not a
velocity; for a non-finite signed speed or negative/non-finite pinning speed;
when a selected `contact_line_speed` curve declares non-velocity abscissa
dimensions; for either material-query refusal, non-dimensionless or
out-of-range angles, or a reversed hysteresis bracket. They query both
endpoints before selecting one, so an incomplete card refuses even when the
current direction would consume the endpoint that is present. The adapter
mutates neither the card nor free-surface state.
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

Dense-grid steps are synchronous (no polling). `SparseGrid3::step_pooled`
polls `Cx::checkpoint` once per kernel tile: a tripped gate drains the pass,
surfaces `SparseError3::Cancelled`, and leaves the pre-step state intact. The
collide pass writes its retained scratch buffer and the stream pass writes a
third retained transactional destination; neither becomes published state
until both passes complete. A re-issued step is bitwise-identical to an
uncancelled one.

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
and fail-closed moving-field admission. A periodic-x, two-moving-wall Couette
deck checks the manufactured linear relative-velocity profile, exact replay,
separate top/bottom impulse and receipt balance, and a measured common-boost
profile/force residual inside one low-Mach lattice envelope. A separately
enumerated partial-saturation cell pins every population of the published
blend, mass and equal-and-opposite impulse/torque/work receipts, origin
translation, exact zero-fraction compatibility, the full-fraction stationary
fixed point, replay, descriptor admission, and atomic refusal. Linear
curved-wall fixtures independently pin both BFL branches, off-lattice torque
arms, full receipt balance, exact halfway compatibility, replay, checked link
construction, exhaustive geometry, far-donor requirements, finite-population
admission, and atomic refusal. A separate topology-transition fixture
pins unique-donor fresh-cell initialization, exact covered/fresh counts,
mass/momentum delta closure, replay, idempotence, and atomic no-donor/mixed-domain
refusal. The file also covers regularized x-face moment/stress reconstruction,
measured-path bit equivalence, and pre-step topology/forcing refusal.

`tests/extensions.rs` covers the current extension scaffolding: power-law and
Newtonian-limit channel profiles, Carreau plateaus, Rayleigh-Bénard onset
bracketing and Nusselt heat transport, gas-neighbor streaming behavior, thermal
wall-temperature queries, invalid-parameter rejection, free-surface mass-ledger
conservation, qualitative dam-break front advance, rotation equivariance,
contact-model bracketing, qualitative jet fragmentation, and dynamic wetting
selection, receipt replay, interface-system identity binding, and typed
refusal paths.

`tests/extensions.rs` (tfz.19): lbm-101 power-law/Carreau profiles vs
analytic (0.12% with τ floor+cap ledger); lbm-102 Rayleigh–Bénard onset
bracket (decay Ra=1200 / growth Ra=2500, Nu>1); lbm-104 STRICT free-
surface mass ledger (5e-14 over 600 dam-break steps, conversions
counted); lbm-105 dam-break front envelope (coarse honesty band);
lbm-106 G3 rotation equivariance (3e-14); lbm-107 contact-model
bracket band + Plateau–Rayleigh jet fragmentation with strict ledger;
lbm-108 level-jump refinement (Poiseuille through the interface +
shear-wave decay-rate transmission, first-order labels); lbm-110 through
lbm-112 dynamic wetting interpolation/pinning, G5 receipt replay and ordered
system identity, and typed incomplete/ambiguous/out-of-domain/unit/bracket
refusals.

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

`tests/d3q19_sparse.rs` drives both entry and mid-stream cancellation. Its
deterministic one-worker runner lets collision finish, requests cancellation
after the first stream kernel tile, and pins unchanged published bytes, zero
completed steps, drained structured refusal, and bitwise serial-equivalent
fresh-gate reissue.

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
  performance or high-Re stability claim. D3Q27, a production cumulant
  collision, normalized drag/lift histories, and bandwidth roofline /
  fs-tilelang kernels remain staged. The sparse active-tile layer
  (`d3q19::sparse`) claims determinism and active-set memory proportionality
  only — NOT throughput (the GLUP/s perf lane is bead 712t): its collision is
  the scalar per-cell kernel, occupancy is whole-tile (per-cell masks and the
  free-surface activation/deactivation rules are WS1-E's increment), tile
  size is pinned at 4³ (the 4³→8³ autotune sweep is follow-on work through
  the fs-exec tune path). The frozen serial sparse golden covers only its
  pinned L-shaped active set and step count; a separate bitwise test covers
  the worker-count sweep. The retained transactional stream destination raises
  population storage from two to three buffers per active tile; this is the
  explicit memory cost of failure-atomic mid-stream cancellation, not a
  throughput optimization. The raw D2Q9 stationary-wall momentum receipt is
  not a normalized aerodynamic
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
  curved SDF. D3Q19 interpolation remains staged; the separate D2Q9 linear BFL
  path consumes caller-supplied link fractions and does not derive intersections
  from an SDF or certify geometric/convergence order.
- D2Q9 `MomentumExchange2` is a raw stationary-wall lattice impulse over the
  caller's exact cell mask. `MovingWallMomentumExchange2` adds a static-topology
  boundary-relative force evaluation for both halfway and caller-described
  linear-BFL links, but it does not move cell topology or initialize newly
  uncovered fluid. Neither
  receipt applies physical-unit conversion, reference-area normalization,
  blockage correction, averaging, or shedding-frequency estimation; the BFL
  path additionally has no quadratic/multireflection option or measured curved-
  geometry convergence deck. Therefore this is not yet the Re=100 cylinder
  Cd/St validation. Boundary-relative exchange improves the force receipt's
  frame behavior. The moving-Couette conformance deck now measures relative-
  profile and shear-impulse residuals for one flat, periodic-x, steady, low-Mach
  common-boost pair; it does not prove Galilean invariance for curved walls,
  topology changes, transients, finite-Re bodies, or the full coupled solver.
- `PartialSaturationCell2` consumes a caller-supplied cell volume fraction and
  one aggregate solid velocity; it does not rasterize geometry, resolve
  overlapping bodies, or attribute the aggregate receipt among bodies. This
  first rung is force-free dense D2Q9 and selects the published
  non-equilibrium-bounce-back variant. It has no stresslet-accuracy,
  convergence-order, free-surface coexistence, rigid-body integration,
  high-Re stability, or FSI claim.
- `UnderResolvedCylinderPair2` is a leading normal squeeze-film correction for
  smooth convex circular or osculating-cylinder geometry, Newtonian fully
  wetted creeping-gap flow, and raw lattice quantities per unit out-of-plane
  depth. Its caller supplies geometry, the gap and cutoff calibration, and
  enumerates each unordered pair exactly once. It neither estimates the
  already-resolved lattice resistance nor mutates fluid populations, so using
  an uncalibrated cutoff or applying the pair twice will double count force.
  The positive gap floor is a disclosed numerical regularization, not a
  roughness, contact, or continuum-validity model. There is no tangential or
  rotational lubrication, many-body hydrodynamics, cavitation, asperity/contact,
  elastic-hydrodynamic, thermal/non-Newtonian, gap-generation, rigid-body
  integration, cutoff-calibration, convergence, or validated-bearing claim.
  Full Reynolds/Elrod-Adams/EHL films remain owned by `fs-tribo`.
- D2Q9 `transition_wall_topology` consumes an already-discretized next wall
  mask; it does not integrate geometry motion, infer covered cells, or couple a
  rigid-body state. Equal-weight one-ring population averaging is a
  deterministic first fresh-cell rung, not a conservative remap or an
  accuracy-order certificate. The receipt deliberately exposes any net active
  mass/momentum change; callers must ledger or correct it rather than assuming
  conservation. Gas/interface transitions and newly covered-body impulse/work
  coupling remain staged.
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
- The dynamic wetting adapter proves system- and endpoint-provenance-bound
  hysteresis admission and directional selection only. Its caller-supplied
  pinning deadband is recorded but has no independent calibration evidence or
  MatDB usage receipt. The adapter does not impose a returned angle in
  `phi_at`, evolve roughness or contamination, implement Cox-Voinov,
  Hoffman-Tanner, or another contact-line law, or claim capillary-rise or
  contact-line validation. The existing binary neutral/wetting solver path is
  behaviorally unchanged.
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

## WS1 3-D feature preservation: thermal + non-Newtonian (bead eg62)

`d3q19::coupled` ports the two VALIDATED 2-D capabilities to the 3-D
path without forking constitutive code. `PlatesGrid3` is a plain-storage
D3Q19 plate channel (rigid z-plates via halfway bounce-back, periodic
x/y) with per-cell relaxation time and per-cell external force, colliding
through the SHARED `collide_cell3` BGK + Guo kernel; the tiled SoA
`Duct`/`BoundaryGrid3` fast paths and their frozen bit-semantics goldens
are untouched. `ThermalLbm3` adds a D3Q7 advection-diffusion temperature
population (c_s² = ¼) with anti-bounce-back fixed-temperature plates and
Boussinesq coupling into the Guo term. `update_tau3` adapts τ from the
local D3Q19 non-equilibrium second moment through the SAME
`rheology::Rheology` law, `TAU_FLOOR`/`TAU_CAP`, and `RheologyStats`
ledger as 2-D; the plate-channel profile is certified against the SHARED
`powerlaw_poiseuille_analytic` (the geometry is 1-D in the wall-normal
axis). Battery (tests/d3q19_coupled.rs, verdict-JSON): Rayleigh–Bénard
onset bracket (decay at Ra = 1200, growth at Ra = 2500, Nu > 1), n ∈
{0.5, 1.0, 1.5} profiles to the 2-D 3% bar, 1e-11 mass ledger, exact
plate temperatures, linear conduction, ledger honesty, and a
replay-stable seeded determinism CANDIDATE.

### No-claim boundaries (3-D coupled ports)

- Single-threaded correctness-first reference rung: no SIMD/tile/parallel
  claim, no performance claim against the tiled fast path.
- No turbulence model and no claim beyond the fixture-scale onset band
  and profile bars actually gated in the battery.
- The seeded determinism hash is a replay-checked CANDIDATE, not a frozen
  golden: freezing awaits the GOLDEN_POLICY four-quadrant ceremony
  (committed tree, M4 + x86, debug + release), per the 40p2 precedent.
- Wall layers hold equilibrium state and are excluded from physics
  claims; only the fluid interior is certified.

## D3Q19 sparse-sweep performance evidence model (bead 712t)

`perf` is the normative semantic substrate for the ignored GLUP/s runner. It
freezes these workload rows before any timing result is inspected:

- one 128³ sparse-grid domain with every 4³ tile active (`dense-active`);
- one 128³ sparse-grid domain whose whole-tile active set rounds upward to at
  least ten percent (`sparse-ten-percent`: 3,277 tiles / 209,728 cells);
- a serial row with exactly one worker and an all-core row retaining the exact
  normalized worker/placement identity.

The version-2 logical traffic model is source-counted, not a compulsory-miss
fantasy. One active BGK update performs two population reads and two population
writes across collide and pull-stream: `19 * 8 * 4 = 608 B/cell`. Splitting
same-tile pulls from cross-tile/domain halo pulls removes sparse-map lookup from
784 of the 1,216 links in one 4³ tile. The remaining 432 halo links are exactly
6.75 lookups/cell, adding `6.75 * 16 = 108 B/cell` of logical source key/slot
traffic and `16 / 64 = 0.25 B/cell` of amortized active-tile identity, for
exactly `716.25 B/cell`. BTree node/cache-line traffic can only increase real
traffic; it is not invented as a portable constant. The generic three-axis
BGK+Guo path counts 142 FLOPs for density/momentum/forced velocity, 271 for
equilibrium, and 611 for force projection plus relaxation, totaling 1,024
FLOPs/update and approximately 1.4293 FLOP/B. Any kernel/model change must bump
`D3Q19_PERF_MODEL_VERSION` and restate this arithmetic.

The plan's 1.0 GLUP/s Apple-M-class and 0.6 GLUP/s Threadripper-class values are
retained as informational `plan_target_glups` fields. They are never silently
used as initial gates. A positive/negative gate requires `EvidenceClass::Citable`,
an exact admission-receipt hash, at least three timed repetitions, stable
critical-path signatures, and a separately authorized positive `floor_glups`.
Before the first honest quiet-window calibration, rows are explicitly
`report_only` with a named refusal and a null floor. Contaminated pre/post axes
are `environment_invalid`, neither floor pass nor floor miss.

Per-repetition tile/halo timing forms a bounded task DAG.
`attribute_critical_path` computes the longest path in the max-plus semiring,
independent of task and predecessor enumeration. Exact path ties choose the
smaller stable task id; exact dominant-class ties choose the stable class order.
Receipts retain the source-to-sink path, makespan, per-class critical wall,
dominant class from `activation|collide|halo|stream`, and its parts-per-million
share. Citable rows require the same path and dominant class across
deterministic repetitions; varying wall values remain telemetry, while
topology/class drift fails closed.

`SparseGrid3::step_pooled_observed` is the opt-in measurement seam for that
driver. It retains the exact executor reports and whole-pass wall envelopes for
the existing collide and combined stream/halo pool invocations. The stream
kernel uses the same worker-count-independent groups as the ordinary path and
partitions its writes into two disjoint loops: same-tile pulls first, then
cross-tile/domain halo lookup and pulls. The ordinary path runs those same
loops without clocks or timing buffers; the observed path adds one clock pair
per completed group and returns canonical `(group, first_tile_slot, tiles,
local_stream_wall_ns, halo_wall_ns)` rows. Group walls overlap across workers
and therefore MUST enter max-plus task-DAG composition; summing them as serial
time is invalid. A missing group observation refuses before the private
destination is published, and cancellation retains the pre-step state exactly.

### No-claim boundaries (D3Q19 performance model)

- `perf` performs no timing, machine probe, authority validation, baseline
  promotion, durable-ledger write, SIMD dispatch, scheduler change, or tuning
  publication. It validates the record the ignored runner must emit.
- Logical byte and FLOP accounting defines a portable roofline denominator; it
  does not claim measured DRAM traffic, cache behavior, instruction count, or
  target attainment on any host.
- No floor is authorized by source code alone. Both reference-ISA rows,
  pre/post admitted axes, quiet-host dispersion, retained external receipts,
  and the sealed tune-row path remain required measurement evidence.
- The observed sweep does not yet time active-set construction, authenticate
  machine axes, turn overlapping group walls into the final critical circuit,
  emit a GLUP/s row, or write the ledger. It exposes truthful primitive timing
  inputs for the ignored driver; it is not itself citable performance evidence.
