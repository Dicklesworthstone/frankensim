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
- Automatic spatial adaptivity and 3-D remeshing remain follow-ons; explicit
  moment-preserving 2-D remeshing is specified below.
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
  FMM-accelerated Biot–Savart, 3-D particle-strength-exchange viscosity, general
  adaptive remeshing,
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

## Moment-preserving 2-D remeshing

`remesh::remesh` redistributes signed particle circulation onto an explicitly
authored `CartesianGrid2` with the tensor-product M4' interpolation kernel.
The stencil reproduces `sum Gamma`, `sum Gamma*x`, `sum Gamma*y` and all three
quadratic position moments to floating-point roundoff, including zero-net
circulation clouds. No net-circulation division or global correction is used.
The unchanged shared velocity core stays separate from grid spacing and the
returned uniform quadrature area `dx*dy`. These are point-vortex moments;
they do not assert momentum in an arbitrary finite fluid control volume.

A nonzero stencil lobe outside the authored grid refuses the complete operation.
No tail clipping, circulation threshold, inferred wall or periodic wrap is
applied. Exact existing grid coordinates reproduce their original node even
with decimal spacing and at the boundary. All grid nodes, including zero
circulation, remain in the result to supply diffusion quadrature. Signed outer
lobes are required for second-moment reproduction: this remap does not promise
positive weights, a maximum principle or enstrophy contraction.

There are exactly sixteen stencil visits per input and
`17*N_input + 2*N_grid` logical work units. `RemeshBudget` admits input count,
grid count, work and logical live bytes before fallible allocation. Only a
particle buffer and an f64 buffer are owned; compensated deposition uses the
latter as scratch before replacing it with uniform areas. Input is borrowed
and immutable. Cancellation polls at entry, within 256 logical visits and
before final publication; a supplied callback can poll the caller's execution
scope. Invalid geometry, inadequate support, arithmetic/resource refusal and
cancellation return no partial cloud.

The interpolation follows Siemaszko et al., *International Journal for
Numerical Methods in Fluids* (2025), equations 34–35,
<https://doi.org/10.1002/fld.70002>. `tests/remesh2.rs` checks all quadratic
monomials, signed cancellation, repeated-grid identity, translation and length
rescaling, induced-velocity refinement, full-support/resource/cancellation
refusal, and an actual remesh–PSE–remesh evolution with explicit new padding.

The returned particles and areas directly feed either checked viscous owner.
Remeshing restores regular quadrature at an explicit time, not an automatic
adaptation schedule. Finite-cloud PSE retains its own spatial, temporal and
domain-padding errors; its circulation conservation does not become exact
impulse conservation after remeshing. The whole-plane diffusion equation has
`d(integral |x|²*omega dxdy)/dt = 4*nu*Gamma_total` for decaying finite moments,
rather than an invariant second moment at nonzero circulation. This bounded 2-D
operation does not implement 3-D remeshing or close the deferred FMM epic.

## Scheduled remeshed viscous evolution

`conservative::remeshed::simulate_with_cx` evolves the complete 2-D trajectory
with the original conservative viscous step and M4' remap owners. The borrowed
`RemeshSchedule` supplies a positive cadence `K` and exactly
`1 + ceil(S/K)` target grids: initial, each complete interval, and the final
partial interval. Zero steps performs the initial remap only; zero time still
executes the authored grid changes. Each resulting grid supplies its own
positive `dx*dy` quadrature before the next PSE step. Viscosity and the shared
velocity core remain independent physical inputs.

Targets may widen, move or change resolution. PSE generally fills old boundary
nodes, so further remapping needs explicit padding for their nonzero stencils.
The driver never clips those tails or infers a wall, strength threshold, moving
window or adaptive-grid law. Decimal-grid changes recover a representable tiny
offset relative to the nearest physical node if the origin-relative quotient
rounds to an integer; offsets that actually underflow still refuse. This is
local arithmetic, with no coordinate snapping or moment projection.

One ambient account admits and charges the whole trajectory. If `N` is the
borrowed input count and `G_k` the successive grid sizes, remap work is
`R = 17*N + 2*G_0 + sum_{k>0}(17*G_(k-1) + 2*G_k)`. Each positive-time
evolution interval adds its original `S_k*G_k*(G_k-1)` pair visits, multiplied
by three for positive viscosity. Total work is `R + pair_visits + S`.
`RemeshedVpmBudget` bounds input/grid counts, steps, total pairs/work and peak
reserved payload before any allocation. At most two particle-plus-f64 banks
are reserved once to the largest declared grid; a single initial-only remap
needs one bank. All later remaps reuse them, with no intermediate trajectory
copies or further allocations. Borrowed inputs/grids and allocator overhead
are excluded from the stated payload bound. All stages poll the same account
within 256 logical visits and publish only the final regular cloud.

The composite inherits signed M4' interpolation and finite-cloud PSE errors.
It does not make diffusion impulse-exact or positive, or claim second-order
convergence for arbitrary changing remesh schedules. Cadence, spatial/kernel
resolution, padding and physical timestep need separate refinement. The G0/G3/G4
checks in `tests/remeshed_viscous.rs` compare explicit original-owner composition,
zero-net inviscid impulse, changing quadrature, whole-run work and combined-bank
admission, and late support/cancellation refusal with deterministic retry.

## Three-dimensional impulse and conservative coarsening

`FilamentWake::hydrodynamic_impulse_per_density` reports the complete closed
wake's geometric hydrodynamic impulse divided by uniform fluid density:
`sum Gamma * oriented_vector_area`, in m^4/s. Multiplication by density gives
kg m/s. Quad orientation matches the physical-edge convention, positive along
the newer row's left-to-right span. Local coordinate differences and compensated
accumulation reduce translation and cancellation error.

`hydrodynamic_angular_impulse_per_density(origin)` reports the corresponding
geometric angular impulse about an explicit origin, in m^5/s (kg m^2/s after
multiplication by density). For localized vorticity this is
`J = -1/2 integral |x-origin|^2 omega dV`; each closed filament contributes
`-Gamma/2 integral |x-origin|^2 dx`. See Carpenter and Guha, *Physical Review
Fluids* **5**, 064702 (2020), [doi:10.1103/PhysRevFluids.5.064702](https://doi.org/10.1103/PhysRevFluids.5.064702).
The exact straight-edge contribution is
`-Gamma/6 * (b-a) * (|a-origin|^2 + (a-origin).(b-origin) + |b-origin|^2)`.
The implementation evaluates an equivalent symmetric triangle decomposition
near each ring, then transports its moment using the existing vector-area
owner. A small ring far from world zero therefore avoids cancellation between
squared world coordinates. Moving the origin by `a` obeys
`J(origin+a) = J(origin) - a cross I`.

Both checked accessors reject invalid geometry, unrepresentable results, and
an attached near/far interface: the retained near sources alone are open.
The angular query also rejects a nonfinite origin. These are localized-vorticity
geometric moments, not momentum or angular momentum inside an arbitrary finite
fluid control volume. Their definitions do not depend on the reduced core
smoothing scale.

`coarsen_oldest` retains all three components of BOTH linear and angular
impulse separately at every span station. It first tries retaining each pair's
oldest node line, with the actual successor boundary. Every original pair and
station must then be representable by the actual returned scalar circulation
under both moment gates. This fast construction telescopes a constant-strength
straight filament exactly, including two rows merged into one; an already
conservative common circulation is retained bit for bit. Equal strengths alone
do not admit warped geometry.

If that bounded pairwise construction fails, node positions use the original
absolute-circulation weights. Each station starts with its pairs' closest
linear-impulse scalar representations, then applies the native rectangular SVD
to find the minimum Euclidean circulation change satisfying its six moments.
The station's aggregate moments are preserved; circulation may redistribute
among retained rows, but impulse is not exchanged with other stations.
`CoarsenMetric::used_stationwise_correction` identifies this correction.
The matrix has six fixed columns, padded when fewer rows are retained; its
factorization does not form a Gram matrix or square the condition number.

Coordinates are differenced about a common station-local origin and normalized
by a local length before area or cubic products. Linear and angular moment
blocks and circulation are then scaled separately. Representable circulation
can therefore survive even when a dimensional `Gamma * area` or angular-moment
product would overflow or underflow. Translation, transverse convection, and
length/circulation units do not introduce a world-axis or dimensional cutoff.

The actual dimensional strengths must pass two independent residual gates:
each moment's vector defect is bounded by `512 * f64::EPSILON` times the sum
of that station's original unsigned panel moments of the same kind. Angular
moments use the common local origin. Neither gate has a dimensional floor, and
a vanishing signed target does not erase its unsigned budget. Singular modes
below `64 * EPSILON` times the largest singular value are omitted only if both
original-target checks still pass.

Failure of these bounded oldest-pairwise/weighted-stationwise constructions
returns `coarsen-impulse-unrepresentable` and retains the original wake. This
is not a search over arbitrary replacement geometries: even a planar unequal-
strength pair can require more retained rows on the tested candidates. A pair
with zero net linear impulse may still carry nonzero angular impulse and cannot
be discarded. Zero-circulation degenerate panels carry neither moment.
Invalid counts, geometry, ages, cores, or nonfinite candidate metrics likewise
refuse without publishing a changed wake. Untouched node/circulation buffers
and the lifting line are reused; metric refusal restores original buffers.
Attached near/far boundaries retain their existing refusal.

The row-wide core smoothing scale includes all three components of middle-node
separation, and age retains the original absolute-circulation weighting.
This smoothing model does not establish conservation of every viscous second
moment or energy. The legacy `station_impulse` remains streamwise-integrated
circulation; only suitable parallel span lines make it proportional to a
physical impulse component. Fixed world-space probe metrics remain sampled
velocity differences, not rotation-invariant error norms.

`tests/coarsen_angular.rs` independently integrates physical edges with exact
quadratic Gauss quadrature. An analytic rectangle detects the old 25% angular
impulse defect while checking useful one-row coarsening. Twisted six-direction
geometry exercises the full station solve, rotations and translations; a ring
at a `2^40` world offset checks local-origin accuracy. Refusal tests include a
pure-angular source with zero linear impulse. Straight refinement and combined
near/far physical edges retain both moments. Existing coarsening, topology,
viscosity, image and field tests retain their own Kelvin, symmetry, unequal-
strength age/core weighting, unit-scaling, deterministic and atomicity oracles,
with enough retained rows to represent the additional angular constraints.

## Connected three-dimensional streamwise refinement

`FilamentWake::refine_streamwise(max_edge_length_m, max_new_rows, cancelled)`
inserts complete node lines until every streamwise edge meets the positive
length target within 64 machine epsilons. Each original row-to-successor
interval uses uniform subdivisions selected by its longest edge, including
zero-strength edges. This retains a conforming connected lattice without
T-junctions. `StreamwiseRefinement` reports row counts and longest edges before
and after the operation. Caller and certificate row limits must admit the
entire change; insufficient capacity refuses instead of partly refining.

Each child retains the parent quad circulation, base squared core and elapsed
age exactly. Inserted spanwise edges have zero circulation after cancellation;
original trailing edges are subdivided collinearly with the same circulation
and effective core. Consequently the represented finite-core velocity,
velocity gradient/vorticity, and geometric linear and angular filament impulse
remain unchanged to roundoff, including for warped quads. Neither dividing
circulation among the children nor interpolating their ages would preserve
that representation. Existing row buffers, the lifting line, shed count and
attached far-field seam remain intact; no physical time elapses.

Work and storage scale linearly with original plus inserted lattice size.
Only inserted row buffers are constructed privately, with cancellation checks
per row and before publication. Failed allocation, budget, geometry, age/core,
coordinate-resolution or cancellation checks retain the original state.
This operation resolves straight represented filaments for subsequent motion;
it does not recover unresolved curvature, coarsen the wake, change the core
model or establish a viscous tube-strain or energy-conservation claim.
`tests/remesh3d.rs` checks independent physical-edge linear and angular
integrals per station, finite-core velocity/Jacobian/curl at multiple core
growth rates, retained material age and buffers, attached near/far continuity,
row limits, idempotence and cancellation/late-refusal atomicity.

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
source assembly, node updates even for zero circulation, and final source
admission: `(2 * substeps * (targets + 1) + 1) * (source_bound + 1)`. Every
active source needs a finite positive core. Both stages limit each full-substep displacement
against the declared fraction of active core radii and resolved edge lengths.
Edges exactly zero-length at interval start are newborn edges: the core guard
applies immediately, and their length joins the edge guard in later intervals.
This permits finite-core release from the line. Excessive motion refuses and
can be resolved by increasing the fixed substep count. The guard is a resolution
condition, not a certified integration-error estimate. The accepted endpoint
sources are assembled and admitted once more, without another velocity sweep:
a resolved edge that collapses only in the final RK update, or an invalid final
evolved core, refuses before publication. Cancellation is checked within 256
source-target contributions and before publishing the private candidate; any
refusal preserves the original state.

The default near-only operation requires a complete near wake; an attached
far-field interface uses the coupled operation below. Its original focused
tests use the independently integrated velocity of a rectangular ring: with
one span edge fixed, the opposite edge rotates at an analytic angular rate.
They check second-order temporal refinement, circulation reversal,
support-relative vector impulse, stage-dependent core spreading,
three-dimensional frame covariance, uniform transport, repeated
positive-circulation emission, exact age accounting and transactional refusals.

### Moving line, incident field and image boundary

`FilamentWake::shed_self_induced_with_environment` adds a `WakeBoundary` and
an incident-velocity callback to the same integrator. Optional endpoint line
nodes prescribe a linear per-node trajectory over the interval. Both RK stages
assemble sources at the corresponding line position, and the new row is
emitted at the endpoint. Free nodes initially coincident with the line remain
free. This represents translation or prescribed linear deformation; it does
not claim exact finite rigid-body rotation. Prescribed line motion obeys the
same per-substep core/edge displacement limit as the material vertices.

The callback receives the actual stage position [m] and elapsed time [s]
within the interval, and returns an additive velocity [m/s] beyond the uniform
ambient value. It must exclude this wake's own induction and images. Spatially
and temporally varying fields therefore participate in both midpoint stages.
Callback errors and non-finite outputs refuse without publishing physical
state; callback side effects themselves are caller-owned.

An optional stationary `ImagePlane` mirrors every assembled physical source
and reverses its circulation, retaining its stage core. Images are rebuilt in
the source scratch each stage, not materialized in physical rows, far cells,
impulse ledgers or digests. The induced normal velocity cancels at the plane.
The caller supplies a boundary-compatible incident field; wake images do not
cancel imposed normal flow. Physical nodes and line trajectories must stay
strictly above the plane at initial, midpoint and final positions. The work
budget counts the additional image sources. These checks do not certify
clearance continuously between discrete stages.

`tests/advection_environment.rs` compares the default environment with the
existing entry point, checks independent space/time-dependent transport and
moving-frame covariance, verifies image dynamics against direct finite-core
line integration, and exercises boundary/callback/work-limit rollback.

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

