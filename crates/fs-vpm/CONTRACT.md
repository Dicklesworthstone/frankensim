# CONTRACT: fs-vpm

Vortex particle methods: 2-D advection and viscous diffusion, plus finite-core
connected 3-D wakes.

## Purpose and layer

Layer L3 (FLUX). Safe Rust with `fs-exec` for the production `Cx` and admitted
budget accountant. Tests use `fs-alloc` to construct scoped contexts.

## Public types and semantics

- `VortexParticle { pos, circulation }` — a point vortex carrying circulation
  `Γ`.
- `VpmBudget` — explicit caps for particle count, RK4 steps, exact attempted
  pair evaluations, and logical peak live payload bytes.
- `simulate_with_cx(cx, particles, dt, steps, core, budget)` — the checked,
  cancellation-aware production transaction. It returns a complete `VpmRun`
  with exact work and ambient-budget consumption, or a typed `VpmError` with no
  partial particle state.
- `induced_velocity(&particles, point, core)` — the desingularized 2-D
  Biot–Savart velocity `Σ (Γⱼ/2π)·perp(r)/(|r|²+core²)`; `core = 0` skips
  coincident particles.
- `advect(&particles, dt, core)` / `simulate(&particles, dt, steps, core)` —
  original unchecked compatibility helpers for analytic fixtures and existing
  callers. They have no explicit work/memory/cancellation authority and are not
  the production entry point.
- `total_circulation(&particles)` — `Σ Γᵢ` (conserved).
- `vorticity_centroid(&particles)` — `Σ Γᵢxᵢ / Σ Γᵢ` (invariant), `None` when
  the total circulation is zero.

## Invariants

- A single vortex `Γ` induces a purely TANGENTIAL field of magnitude `Γ/(2πr)`
  (`u·r = 0`) and does not advect itself.
- A counter-rotating PAIR separated by `d` self-propels in a straight line at
  exactly `Γ/(2πd)`, preserving its separation (the 2-D analog of a vortex
  ring's translation).
- Total circulation is conserved to roundoff; a co-rotating pair conserves its
  vorticity centroid.
- The desingularized kernel keeps the velocity finite at/near a particle.
- Deterministic (fixed particle order, no RNG).
- The checked path attempts exactly `4 S N²` source-target contributions for
  `N` particles and `S` steps, including coincident/self visits that the kernel
  skips numerically.
- Its exact logical-work charge is `N + S(4N² + 4N + 1)`: one combined input
  validation/copy visit per particle, four all-pairs sweeps, three stage
  displacements plus one final combine, and one private step commit.
- For `S > 0`, the reusable-buffer implementation admits a logical peak of
  `N(5P + 4V)` bytes, where `P = size_of::<VortexParticle>()` and
  `V = size_of::<[f64; 2]>()`; for `S = 0`, it admits `NP` bytes. This is
  operation-owned vector payload, not allocator-realized resident memory.
- Input order, target order, source order, and RK4 arithmetic association match
  the compatibility path. Valid admitted fixtures therefore have bitwise-equal
  results on the same build and ISA.


The 2-D kernel has no absolute distance cutoff. Only exact coincidence gives
zero self velocity; a finite nonzero separation keeps its signed tangential
field at any admitted scale. Checked and compatibility paths use the same
kernel and retain the original arithmetic on ordinary finite intermediates.
The exceptional-range path normalizes the squared distance with the original
core radius, then scales each velocity component by binary exponents. This
retains finite fields when a squared length, scalar vortex strength, or
normalized transverse separation alone would overflow or underflow. It adds no
buffers, pair visits, or logarithms. The checked owner still refuses coordinate
differences, actual velocity sums, or RK stage positions that do not fit f64;
it does not turn an unrepresentable physical field into zero.

The velocity and checked-dipole regressions in `tests/vpm.rs` cover signed
scaling from `1e-180` through `1e300`, core-radius retention, subnormal source
geometry, exceptional component ranges, the original ordinary arithmetic, and
transactional refusal. These checks concern the represented algebraic vortex
model, not the continuum validity of fluid mechanics at arbitrary lengths.

## Error model

The checked path validates finite particle coordinates/circulations and finite
`dt`; signed finite `dt` is allowed. `core` must be finite and nonnegative;
its square need not be representable. Count/work/byte arithmetic is checked
before allocation.
Local cap, ambient budget, allocation, and non-finite-intermediate failures are
typed. All owned vectors use fallible reservation, and an error publishes no
particle vector.

The compatibility helpers retain their original infallible signatures and
ordinary `Vec` allocation behavior. `vorticity_centroid` returns `None` on zero
total circulation.

## Determinism class

Same-build, same-ISA deterministic: fixed input/target/source order, fixed RK4
association, and no RNG or worker-count dependence. No cross-ISA bit-stability
claim is made.

## Cancellation behavior

`simulate_with_cx` admits the ambient cost quota against its complete checked
work plan, then checkpoints before allocation, before each chunk of at most
`VPM_WORK_CHECKPOINT_STRIDE = 256` logical units, at RK-stage and step
boundaries, and immediately before publication. Completed units are charged
before the next checkpoint. Cancellation/deadline/poll refusal latches through
`AdmittedBudget`; synchronous private buffers drain by stack unwinding, and no
partial state escapes.

The cancellation bound is in logical work units, not microseconds. The
compatibility `advect` and `simulate` helpers do not accept a `Cx` and make no
cancellation claim.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/vpm.rs` contains the seven original analytic/determinism fixtures plus G0
checked-path evidence for bitwise legacy equivalence, exact work/byte receipts,
zero-step identity, cap and arithmetic-overflow refusal, invalid inputs, and
non-finite intermediates; G4 evidence covers pre-cancellation and ambient
cost/poll refusal without partial publication.

## No-claim boundaries

- The original 2-D advection core uses DIRECT `O(N²)` Biot–Savart. The FMM /
  treecode acceleration for large particle counts (the whole point of the
  method at scale) is the fuller deliverable, staged.
- Three-dimensional particle-strength-exchange viscosity and dynamically stretched 3-D
  filaments/rings remain deferred. The connected 3-D wake's prescribed
  geometry and reduced finite-core model are specified below.
- The HYBRID BEM+VPM airfoil (pitching/plunging unsteady lift vs Theodorsen —
  the flapping credential) composes this core with fs-bem downstream.
- Spatial adaptivity (particle insertion/merging, remeshing) is a follow-on.
- Logical live bytes exclude the borrowed input, allocator metadata, capacity
  rounding, allocator fragmentation, and process RSS. Admission cannot promise
  allocator success; reservation failure remains typed.
- No integration-error, convergence-order, stability/CFL, conservation
  certificate, wall-clock deadline-latency, resumable-partial-state,
  cross-ISA-bit-stability, or performance claim is made by the checked receipt.
- G4 cancellation evidence is synchronous and deterministic. It does not yet
  establish cancellation under multi-threaded signal races or a wall-time
  latency bound.
- Retired from v1 (owner decision 2026-09-01, bead frankensim-rc-root-q61wp.25):
  FMM-accelerated Biot–Savart, 3-D particle-strength-exchange viscosity, remeshing,
  and the BEM+VPM hybrid are NOT scheduled; the promise is preserved in deferred
  epic frankensim-rc-root-q61wp.43 and un-parks only behind a consumer journey.

## Momentum-preserving 2-D inviscid advection

`pairwise::advance_pair([a, b], dt, core)` solves an isolated pair of the existing
algebraic-core Biot–Savart model without allocation. Both vortices use the
same declared core. Their separation is invariant, giving constant angular
rate `(Gamma_a + Gamma_b)/(2 pi (d² + core²))`. Equal and opposite strengths
produce exact dipole translation. A sinc/cosc displacement form remains
continuous as the total circulation approaches zero; no circulation floor or
division by the cloud's net circulation is used. Strict `fs_math::det` trig,
normalized lengths and an exponent-based exceptional-range fallback avoid
platform-libm variation and spurious overflow of squared lengths or summed
like-signed strengths. Invalid or nonrepresentable pairs refuse.

`conservative::simulate_with_cx` composes a lexicographic forward half sweep
with its exact reverse half sweep. The result is an explicit second-order,
time-reversible alternative to the unchanged fourth-order RK4 owner. It
preserves each circulation, the first moments `sum Gamma*x` and
`sum Gamma*y`, and the quadratic moment `sum Gamma*|x|²`, to floating-point
roundoff. This remains valid for signed clouds with zero total circulation.
These vortex impulse invariants do not assert the momentum of fluid in an
arbitrary finite control volume or a body/wall reaction force. Many-vortex
energy and finite-step relabeling invariance are not exact; trajectory and
energy accuracy still need time refinement. Zero-core coincident vortices
retain the existing zero-induction convention, not a resolved collision model.

Each nonzero step costs `N(N-1)` pair visits plus one step unit; input admission
costs `N`. Zero time skips pair visits. Only one particle vector is owned,
with no per-step allocations. The same `VpmBudget`/`Cx` caps, cancellation
checkpoints and final publication checkpoint apply; failures return no partial
state. A lower pair count is not a measured wall-clock speed guarantee.

The exact pair-Hamiltonian construction and symmetric composition follow
J. Frank, *Symplectic and energy conserving discretizations of Hamiltonian
PDEs*, pp. 11–12,
<https://webspace.science.uu.nl/~zegel101/NUMpde2025/Lecture8Extra2.pdf>.
The radial pair flow here is derived for this crate's shared regularized core.
`tests/pairwise.rs` checks analytic rotations/dipoles, signed impulse closure,
reversibility, dimensional/frame covariance, independent RK4 references and
second-order three-vortex refinement (G0/G1/G3).
`tests/conservative.rs` exercises the checked transaction, one-buffer/work
receipts, identity/refusal cases and long-time signed-cloud impulse.

## Conservative 2-D viscous diffusion

The diffusion-focused implementation adds `diffusion::PseDiffusion` and two
checked entries without changing the existing inviscid model:

- `diffuse_with_cx(cx, particles, areas_m2, dt, steps, law, budget)` evolves
  circulation on fixed quadrature sites.
- `simulate_viscous_with_cx(cx, particles, areas_m2, dt, steps, core, law, budget)`
  composes two diffusion half steps around the original checked RK4 advection.

`PseDiffusion::new(nu, epsilon)` declares a constant nonnegative kinematic
viscosity in m²/s and a positive Gaussian kernel width in metres. Every original
particle has an independently supplied positive quadrature area in m²;
`omega_i = Gamma_i/A_i`. The area stays with that particle under the existing
incompressible advection. Areas, the diffusion width and the algebraic velocity
core are distinct inputs. No material property or quadrature support is inferred
from a particle's circulation.

The spatial law is

```text
dGamma_i/dt = sum_j k_ij (A_i Gamma_j - A_j Gamma_i)
k_ij = 4 nu exp(-|x_i-x_j|²/epsilon²) / (pi epsilon^4).
```

The Gaussian's second moment is 2 along each axis, giving the Laplacian in the
resolved quadrature/kernel limit. This uses the positive integral-operator PSE
construction of Degond and Mas-Gallic (1989), DOI
`10.1090/S0025-5718-1989-0983559-9`. Each unordered pair is solved analytically:
its vorticity difference decays as `exp(-k_ij (A_i+A_j) dt)`. A forward half
sweep followed by the fully reversed half sweep gives second-order splitting
at fixed sites. The moving entry uses `D(dt/2) A(dt) D(dt/2)` with the same
circulation inventory, checked advection buffers, work meter and final acceptance.

Each pair transfers equal and opposite circulation and contracts the
vorticity range and area-weighted enstrophy. These properties hold for signed
vorticity and for stiff diffusion without an explicit diffusion CFL or clipping,
up to floating-point roundoff. `expm1` retains small exchanges; scaled area
weights and a logarithmic fallback avoid spurious overflow/underflow in stiff
or rescaled inputs. This is an integrator property, not permission to skip time
refinement. Zero viscosity recovers the original advection bits and work counts.
The Gaussian, exact relaxation and logarithmic fallback use the strict
`fs_math::det::{exp, expm1, ln}` functions, avoiding platform-libm variation.
The finite-step splitting error depends on particle order and converges at
second order; relabeling particles is not an exact discrete symmetry.

Successful pure-diffusion runs use `S N(N-1)` unordered pair visits; the moving
entry adds `2 S N(N-1)` visits to the existing `4 S N²` advection visits. Identity
diffusion skips these visits. Each visit is one logical work unit under the
existing `VpmBudget` and ambient `Cx` budget. Pure diffusion owns one particle
vector; coupled advection reuses the original scratch without another buffer.
There is no dense pair matrix, pair cache, tail cutoff or per-step allocation.
Cancellation is checked within the existing 256-work-unit bound, and failed
input/domain/budget/finite-output checks publish no state.

`tests/diffusion.rs` checks unequal-area analytical pair relaxation, signed
maximum principle and weighted enstrophy, Gaussian/Lamb–Oseen spatial
refinement, independent three-site temporal modes, dimensional rescaling, and
zero/positive-viscosity coupled trajectories. These are G0/G1/G3 implementation
checks, not measured-fluid validation or a machine-performance claim.
`tests/viscosity2d.rs` additionally exercises the complete moving-cloud
advection/diffusion transaction against Gaussian spreading, unequal-area
analytical relaxation, translation and dimensional scaling, permutation-error
time refinement, subnormal zero-time identity, exact costs and atomic refusals.
`cargo run -p fs-vpm --example viscous_vortex` runs an illustrative vortex with
declared SI kinematic viscosity and reports circulation, mean-square radius,
the Lamb–Oseen continuum width, enstrophy, and admitted work/memory consumption.
Callers starting from dynamic viscosity must divide by their actual fluid
density; this example does not infer a material specification or remesh areas.

Scope: constant-viscosity 2-D bulk diffusion on an explicitly supplied finite
cloud. Continuum accuracy requires overlapping, sufficiently regular quadrature
and domain padding; advection can degrade that regularity. Diffusion cannot
populate absent sites. Remeshing, physical walls, particle insertion,
variable-viscosity stress, 3-D stretching and viscous boundary shedding remain
outside this entry. Pair conservation does not give exact linear/angular impulse
on arbitrary finite clouds. The velocity core remains algebraic, not a Gaussian
heat kernel. This bounded 2-D addition does not close the deferred 3-D/FMM/BEM
epic or change the existing inviscid thrust-screening consumer.

## Viscous evolution with conservative advection

`conservative::simulate_viscous_with_cx` composes the original PSE diffusion
half steps around the original exact-pair conservative advection sweep:
D(dt/2), A(dt), D(dt/2). Quadrature areas remain explicit and distinct from
the Biot–Savart core. Time is nonnegative; every area is validated even for
zero motion or zero viscosity. Within this domain, zero viscosity recovers
the existing conservative trajectory, work counts and live-memory bound.

Positive time and viscosity require `3*S*N*(N-1)` pair visits and
`N + S + 3*S*N*(N-1)` logical work units. One particle vector is owned, with
no RK4 buffers or per-step allocation. Original local/ambient budgets,
bounded checkpoints and atomic publication apply.

Advection preserves its vortex-impulse moments; PSE preserves circulation
and contracts its vorticity range and weighted enstrophy. PSE does not
exactly preserve linear/angular impulse on an arbitrary finite cloud, and
the composition does not hide that quadrature error with a projection.
Temporal and cloud/kernel refinement remain distinct requirements.
No complete fluid/body momentum balance, viscous heat, wall condition,
evolving-core model or measured performance claim is added.
`tests/conservative_viscous.rs` checks zero-viscosity equivalence, analytical
dipole decay and travel, refinement, signed-cloud behavior and refusals.

## Three-dimensional impulse and conservative coarsening

`FilamentWake::hydrodynamic_impulse_per_density` reports the complete closed
wake's geometric hydrodynamic impulse divided by uniform fluid density:
`sum Gamma * oriented_vector_area`, in m^4/s. Multiplication by density gives
kg m/s. Quad orientation matches the published physical-edge convention,
positive along the newer row's left-to-right span. Local coordinate differences
and compensated accumulation reduce translation and cancellation error.
The checked accessor rejects invalid geometry, unrepresentable results, and
an attached near/far interface: the retained near sources alone are open.
This is neither angular momentum nor momentum inside an arbitrary finite
fluid control volume.

`coarsen_oldest` now uses the actual three-dimensional quad areas at every
station to solve `Gamma_new A_new = Gamma_a A_a + Gamma_b A_b`. It checks all
three components against a vector-wide floating-point roundoff scale.
Area and circulation are normalized separately before combining their products.
A representable merged circulation is retained even when the dimensional
`Gamma * area` intermediate would overflow, underflow, or lose accuracy in the
subnormal range. The coarse representation still refuses incompatible vectors.
Rotated or transverse wakes and small length units no longer depend on an
x-spacing assumption or a fixed dimensional spacing cutoff. A single merged
circulation cannot reproduce a vector outside its new quad's area direction;
such pairs refuse with `coarsen-impulse-unrepresentable` and retain the
resolved wake. Invalid counts, geometry, ages, cores and degenerate candidate
areas also refuse before any mutation. Attached near/far boundaries retain
their existing coarsening refusal.

Merged node positions use absolute-circulation weights. The row-wide core
smoothing scale uses the middle node's separation in all three coordinates.
That geometric smoothing model does not establish conservation of every
viscous second moment on warped rows. The legacy `station_impulse` remains
streamwise-integrated circulation; only suitable parallel span lines make it
proportional to a physical impulse component. No angular-impulse or energy
conservation claim follows from this coarsener.

The focused momentum/coarsening tests use independent polygon edge integrals,
an analytic rectangular ring, proper rotations and translations, transverse
and rescaled wakes, near/far source accounting, and atomic refusal of warped
or malformed inputs. Existing topology, viscosity and field tests retain
their own source-model checks.
Additional physical-edge integrals isolate each span on sheared, nonuniform-width
wakes through repeated coarsening and length rescaling. A late incompatible pair
checks atomic refusal after an earlier pair is prepared. Analytical rectangular
pairs check independently scaled area/circulation, including signed cancellation.

## Self-induced attached near-wake motion

`FilamentWake::shed_self_induced` convects the existing connected near wake
using its own finite-core Biot–Savart velocity plus a uniform imposed flow,
then emits the supplied circulation at the fixed lifting line. Each shedding
interval is `1 / shed_hz` seconds with an explicit positive substep count.
Every substep evaluates the old source sheet, moves a private midpoint sheet,
rebuilds its physical edges, and evaluates that midpoint before updating nodes.
The captured old circulations remain constant during the interval; the supplied
new circulation enters at its endpoint. The midpoint scheme is second order
for the smooth held-circulation problem, without an order claim across arbitrary
discontinuous changes of boundary circulation.

`SelfInducedShed` supplies the emitted positive core radius and the existing
`WakeCoreEvolutionMode`. Existing base cores remain intact; initial and
midpoint evaluations use their elapsed age plus the stage fraction. On
publication existing row ages advance once and the emitted row has age zero.
Moving shared material vertices transports the filament geometry, including
line stretching, while preserving every retained quad circulation and nodal
Kelvin balance. The fixed line can exchange hydrodynamic impulse with the
wake; full vector impulse conservation is not imposed on this supported system.

Admission includes a conservative direct-work bound covering both RK stages,
source assembly and node updates even for zero circulation. Every active source
needs a finite positive core. Both stages limit each full-substep displacement
against the declared fraction of active core radii and resolved edge lengths.
Edges exactly zero-length at interval start are newborn edges: the core guard
applies immediately, and their length joins the edge guard in later intervals.
This permits finite-core release from the line. Excessive motion refuses and
can be resolved by increasing the fixed substep count. The guard is a resolution
condition, not a certified integration-error estimate. Cancellation is checked
within 256 source-target contributions and before publishing the private
candidate; any refusal preserves the original state.

This operation is a direct all-pairs near-wake solver. An attached far-field
interface refuses rather than moving only one end of shared filaments. Images,
nonuniform external fields, remeshing and resolved viscous tube-strain dynamics
remain separate. Its focused tests use the independently integrated velocity
of a rectangular ring: with one span edge fixed, the opposite edge rotates at
an analytic angular rate. They check second-order temporal refinement,
circulation reversal, support-relative vector impulse, stage-dependent core
spreading, three-dimensional frame covariance, uniform transport, repeated
positive-circulation emission, exact age accounting and transactional refusals.

## Finite-core connected 3-D wake

`FilamentWake` stores captured **quad strengths**, positive along the current
lifting line's left-to-right span. Velocity evaluation first assembles unique
oriented physical edges: a row's spanwise circulation is the older quad strength
minus the newer one; a trailing edge carries the right-minus-left spanwise
difference; the current line closes the sheet with the newest strength. This is
the boundary of an oriented vortex-ring sheet, not a separate full-strength
spanwise filament for every time sample. In particular a steady wake with a
common core evaluates to one outer ring, independently of its row subdivision.
The current line's closing segments are included; downstream consumers must
not add those same bound segments again.

`visit_segments_with_core_growth` exposes that same allocation-free source
traversal to exact evaluation and far aggregation. Shared spans cancel **before**
kernel evaluation, so distinct adjacent core widths/ages cannot create fake
transverse vorticity for steady circulation. A shared span uses the newer row's
core and age, and trailing edges use their emitting row. This is a declared
reduced edge-core model, not a resolved viscous vortex sheet.

`induced_kinematics_with_core_growth` sums the analytic finite-segment velocity
Jacobians through this same traversal. `WakeKinematics` carries physical velocity
in m/s, the spatial Jacobian in 1/s, and exposes its curl (vorticity) and trace
(divergence). There is no probe-spacing parameter or finite-difference stencil.
Invalid probes/growth, singular zero-core encounters, or unrepresentable outputs
return `None`. A positive-core axis has finite nonzero vorticity even when its
segment velocity vanishes. G2 tests compare derivatives to the independently
differenced whole field and check zero off-source curl for a closed inviscid
sheet. This exact near-source API does not differentiate multipole far cells.

`node_circulation_residual` measures incoming minus outgoing circulation on the
actual assembled edge strengths, including starting, tip, and closing-line
nodes. The older `cell_net_circulation` entry point now reports this incident-edge
balance at the named cell's newer right-hand vertex. Small roundoff residuals
can occur; it no longer reports a tautological subtraction unrelated to the
velocity sources. `tests/topology3d.rs` independently assembles endpoint
incidence and compares velocities to separately integrated closed rings, with
nonuniform time-varying strengths and the hydrodynamic impulse of the actual
filaments as sign/closure checks (G0/G2/G3).

After aggregation, `downstream_gamma` retains the last transferred quad's
strengths at the shared near/far interface. The far field owns its trailing
edges and the near wake owns the single shared span. Repeated conversion
preserves the exact combined field, including mixed cores and image reflection.
Coarsening an already attached oldest row refuses with
`coarsen-interface-attached`: moving it alone would disconnect the retained far
segments. Coarsen before aggregation. As before, conversion is a snapshot;
continued prescribed convection uses `FarField::shed_prescribed` to keep the
near/far interface connected. It translates existing near and far sources by
the same displacement, advances old cores over `shed_dt_s`, and adds a zero-age
row at the current lifting line. Geometry, age, and diffusion are preflighted;
a refusal leaves both input states unchanged without cloning the wake.

`segment_kinematics_with_core` returns the finite-segment velocity and its
exact spatial Jacobian, with no finite-difference step or allocation. The
exterior derivative preserves the rationalized endpoint expression.
Positive-core axial and endpoint probes retain a finite, generally nonzero
curl even when their velocity is zero. Singular on-segment zero-core probes
and unrepresentable results return `None`; regular extended-axis probes
remain available. Scaling lengths by s gives velocity /s and Jacobian /s²
at fixed circulation. The six `tests/segment_kinematics.rs` G0/G1 checks use
independent differentiated line quadrature, the long-filament core profile,
Stokes circulation, closed-loop off-filament curl and geometry scaling.
The kernel convention follows Lindsay and Krasny (JCP 172, 2001), equations
(1)-(2), DOI 10.1006/jcph.2001.6862. Regularization gives algebraic vorticity
tails; a visualization core mask is not a claim that curl vanishes outside it.

The connected filament, aggregation, and image lanes use the Rosenhead--Moore
kernel `dl cross r / [4 pi (|r|^2 + a^2)^(3/2)]`. Straight-segment integration
is analytic. Its unit-tangent form avoids comparing incompatible length powers,
and its exterior endpoint difference is rationalized to retain the small field
near the segment's extended axis. There is no absolute length cutoff. Exactly
axial or degenerate segments return zero, including the singular zero-core
convention. Geometry and core radius scaled by s give velocity scaled by 1/s
at fixed circulation.

`WakeCoreEvolutionMode::CoreSpreading` admits finite effective kinematic
viscosity in `[0, MAX_NU_EFF]`, in m^2/s. For a fixed declared viscosity,
`a^2 = a0^2 + 4 nu t`; zero viscosity recovers Frozen velocities. This is a
reduced algebraic core-spreading model, not exact Gaussian Lamb--Oseen
diffusion, particle-strength exchange, or resolved Navier--Stokes diffusion.
The core width is a regularization scale; it is not a finite physical second
moment of the algebraic profile. Coarsening also contributes a geometric
smoothing scale without changing that interpretation.

Each shed row retains elapsed `age_ticks`, starting at zero and increasing
only on a successful shed. Coarsening preserves circulation-weighted age
alongside position and core width. Relative weights are normalized without a
dimensional circulation floor, so changing circulation units cannot alter
geometry, core width, or age. Aggregation uses retained age and its declared
`shed_dt_s`; removing rows never resets the ages of those left near. Public
hybrid evaluation and image snapshots apply the same viscosity and shedding
interval to the remaining near wake and the aggregated sources.

Every far cell retains each original source segment and its own squared core.
Exact evaluation, near-field fallback, and prune audits use those individual
cores. The far expansion uses a circulation-length-weighted common core,
exact line monopole/dipole/quadrupole moments (including the intrinsic
line quadrupole), and the leading core-variation moment. Spatial/core mixed
terms remain in the truncation. Probes inside the geometric bounding sphere,
or with maximum squared-core deviation greater than
`(|r|^2 + common_core^2)/16`, use exact retained sources. This is an expansion
domain guard, not a certified total-velocity error bound. Separated mixed-core
cells retain constant work per cell rather than always scanning their sources.
The existing singular contribution bound remains conservative for regularized
sources; pruning still independently audits the exact contribution.

`FarField::advance_core(dt)` adds continued diffusion to already aggregated
cells in O(number of cells), without scanning their segments or rebuilding
moments. A common offset and cached maximum base core provide transactional
finite checks; zero time, Frozen, and zero viscosity preserve state exactly.
The method does not convect geometry. Callers keep near shedding and far
diffusion on the same elapsed-time schedule, and advance existing far cells
before appending newly aggregated cells to avoid aging those new cells twice.
`FarField::translate(displacement_m)` moves cell centers and retained exact
endpoints transactionally, preserving relative multipole moments and source
cores. Its cost is linear in retained sources, with no reconstruction or
allocation. The coupled `shed_prescribed` operation performs this translation
and diffusion on the same interval as near shedding; it does not model
free-wake self-induction. Images are evaluation snapshots and must be rebuilt after geometry, diffusion,
or mode changes. Mirroring preserves source cores, diffusion offsets and the
axial-vector core moment.

The filament and far-field digest domains are v3 because circulation assembly
and edge orientation now describe the closed sheet rather than the former
nonclosing row sources. Filament identity includes the closing-line geometry
and retained interface circulation as well as row geometry, strengths, cores,
and ages. The updated replay goldens cover this intentional physical change.
`tests/viscosity3d.rs` checks independent line quadrature, the analytic long
filament, length scaling, thin exterior fields, exact/near/far/image agreement,
mixed-core convergence and fallback, chronology through conversion, continued
diffusion and transactional refusals. These checks do not establish free-wake
stretching, a viscous particle solver, cross-ISA bit stability, or a wall-clock
performance guarantee.

