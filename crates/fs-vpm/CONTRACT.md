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

### Optional fourth-order conservative time composition

`conservative::simulate_with_order_cx` and
`conservative::ambient::simulate_with_order_cx` accept
`CompositionOrder::{Second, Fourth}`. Existing entry points retain Second,
including their arithmetic, resource counts and checkpoint sequence. Fourth
composes the original reversible second-order step `S` as
`S(w1*dt) S(w0*dt) S(w1*dt)`, with `w1=1/(2-cbrt(2))` and `w0=1-2*w1`.
The conditions `2*w1+w0=1` and `2*w1^3+w0^3=0` cancel the leading splitting
error. The construction follows Yoshida, *Physics Letters A* 150 (1990),
equations (4.1), (4.5)-(4.6),
<https://doi.org/10.1016/0375-9601(90)90092-3>.

All stages use the original equal-and-opposite pair flow; no state projection,
circulation correction or division by net circulation is introduced. Thus
linear and angular point-vortex impulse remain conserved to roundoff in the
unforced lane, including signed zero-net clouds and passive tracers. The held
affine lane composes the complete source/pair step at fourth order. Adjacent
exact maps of the same source combine into four passes: outer weights
`w1/2` and inner weights `(w1+w0)/2`. The source moment receipt includes every
signed source update in its original double-double accumulator. A prescribed
source can still change impulse; no body reaction or closed-system momentum
claim is inferred from that receipt.

The fixed positive and negative weights are kept separate from `dt`. Pair
updates use the original scaled time-fraction seam with the sign in `dt`;
affine generators apply each weight to the original double-double half-product.
This avoids erasing resolvable motion when `dt` is subnormal but circulation or
source velocity is large. Both distinct signed affine maps are prepared before
evolution. Unrepresentable internal maps or updates refuse atomically, including
a nonzero source coefficient that underflows to zero. Negative internal stages
are admitted only by the reversible inviscid/source path. Viscous/PSE evolution
and remapping retain their original nonnegative-time methods; a fourth-order
request cannot carry a diffusion law. Scheduled inviscid evolution can select
this composition between its original positive-time remap boundaries, as
specified below.

For `N` particles and `S` nonzero-time physical steps, Fourth admits exactly
`3*S*N*(N-1)` pair visits, plus `4*N*S` source visits when the held affine field
is active. Input validation/copy costs `N`; physical step boundaries cost `S`.
One original particle vector is retained, with no stage or per-step allocation.
These are deterministic operation/payload counts, not a wall-clock speed claim.
Zero-time, zero-step, empty-cloud and zero-source paths preserve their original
inactive accounting. Fourth-order accuracy requires smooth finite-core or
noncolliding evolution and a held source; it is not an adaptive error bound,
exact many-vortex energy conservation or an order claim across changing sources,
viscosity or remapping.

`tests/conservative_order.rs` (G0/G1/G3/G4) checks independent full-ODE fourth-order
refinement with and without a noncommuting affine field, improved trajectory
accuracy at equal pair work, 3,000-step signed impulse balance, reversal,
analytic first-moment evolution and source accounting, frame/unit covariance,
subnormal-time transport, exact resource admission and cancellation inside
signed pair and coalesced source stages. Original affine/viscous/remeshed
regressions remain applicable to the unchanged Second path.

### Adaptive positive-time inviscid windows

`conservative::adaptive::simulate_with_cx` and
`simulate_with_affine_flow_cx` complete an explicit finite nonnegative duration
with the original selected Second/Fourth pair maps and, optionally, one held
affine source. `AdaptiveOptions` declares initial/minimum/maximum physical
timesteps, a positive absolute position tolerance in metres and a nonnegative
relative tolerance. Every trial compares one selected-order full step with two
half steps. For order `p`, its local estimate is

`max_i ||x_fine_i-x_coarse_i|| / ((2^p-1)*(atol+rtol*max(R_old,R_fine)))`.

`R` is the maximum Euclidean radius about the unweighted particle centroid.
The centroid uses double-double accumulation; the norm uses deterministic
`fs_math::det::hypot`. The scale is objective under translations and rotations,
transforms with physical length units, admits signed zero-net-circulation
clouds and does not relax when a uniform source moves the whole cloud far from
the origin. The estimate includes passive tracers. A ratio at most one accepts
the fine state unchanged: there is no Richardson state extrapolation, impulse
projection or circulation correction. Thus accepted unforced steps retain
the original circulation and linear/angular impulse invariants to roundoff.
Rejected trials halve the physical step; accepted ratios at most `1/8` or
`1/32` permit doubling at Second or Fourth order. Limits bound all positive
trials, except that the final remainder may be shorter than the minimum.
Failure to meet tolerance at that floor refuses the whole window.

Original composition weights remain separate from the physical timestep, so
the two half steps retain resolvable subnormal-time motion. Endpoint time uses
a double-double accumulator and a downward-rounded final remainder when
needed; a remaining low part is advanced explicitly rather than silently
discarded. Original signed Fourth stages are internal inviscid/source maps.
No negative step reaches diffusion or remapping. Original numerical refusals
propagate; the controller does not hide unrepresentable source maps or geometry
by returning a partial trajectory.

The affine path prepares both full/fine source-map sets once per trial and
combines adjacent exact source maps at the fine midpoint. Only accepted fine
source increments enter the original double-double moment account. Rejected
and coarse source visits still count as computational work. A trial uses
`3*q*N*(N-1)` pair visits with `q=1` or `3`, plus `5*N` or `11*N` source visits
when active. Initial validation/centered geometry costs `2*N`; each trial adds
`3*N+4` copy/geometry/boundary units beyond pair and source work. Three particle
banks are fallibly reserved once and reused; no trial allocates. Empty,
zero-duration and unforced singleton windows are identities with one bank and
`N` validation units. Zero source takes the exact unforced path.

The owner reuses the original adaptive budget, result and error vocabulary,
including a zero `fraction_evaluations` receipt because this lane evaluates no
PSE kernels. It shares the dyadic step controller with the fixed-site and
moving-viscous lanes; their original Second-order threshold is unchanged.
Whole-window particle, attempt, pair, source-visit, work and live-payload caps
are checked before each candidate's work; one admitted `Cx` account charges
all accepted and rejected candidates. Original bounded work checkpoints and
the final publication checkpoint make cancellation atomic. Deterministic
operation counts and the three-bank payload are stated; no wall-clock speed
claim is inferred. The local estimate is neither a certified/global trajectory
bound nor a spatial/core/continuum-error bound. Variable accepted steps do not
establish time reversibility or symplecticity. Held prescribed sources retain
their original external-impulse interpretation; no body reaction, moving-source
adaptation, adaptive remapping or broad 3-D VPM completion is claimed.
Viscous trajectories retain the separate original positive diffusion owner
and its independent position/vorticity tolerances described below.

`tests/adaptive_conservative.rs` (G0/G1/G3/G4) compares both orders with an
independent refined full-ODE reference for close vortices and noncommuting
strain, checks actual error improvement under tighter local tolerances and
lower admitted work for Fourth on that fixture, signed zero-net impulse and
analytic source first moments, equivalence to original fine maps, frame/unit
covariance, inactive/subnormal/short-endpoint motion, exact cap admission,
deterministic retry, tolerance/numerical refusals and cancellation throughout
pair/source trials and final publication. Existing adaptive fixed-site PSE,
moving-viscous and conservative/affine/viscous/remeshed tests remain applicable.

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

`PseDiffusion::compact_wendland(nu, radius)` selects a positive C2 compact
kernel in the same diffusion owner. Its pair coefficient is
`1008*nu/(5*pi*radius^4) * (1-q)^4*(1+4*q)` for `q=r/radius < 1`, exactly
zero at and beyond the radius. The radial moment integral is `5/504`, so the
normalized Cartesian second moment is 2 as for the Gaussian. This uses the
Wendland C2 shape (1995, DOI `10.1007/BF02123482`) with its diffusion
normalization derived here. `kernel()` and `support_radius_m()` expose that
physical choice. The original constructor continues to select Gaussian.

The compact law uses the same conservative analytical pair map, bounded
execution and Strang compositions, including the existing conservative,
ambient-flow and remeshed evolution entries. Kernel width denotes support
radius for this law; the same numeric width in another kernel family does
not imply the same regularization. Continuum accuracy still requires support
overlap, resolved quadrature and radius/cloud refinement. Near the support
edge, `(1-r²)/(1+r)` preserves the small polynomial factor; absent interactions
are exact identity, including signed-zero circulations.

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
`tests/compact_diffusion.rs` recovers the kernel's second moment from actual
pair relaxation and independent polar quadrature, checks Gaussian heat
refinement, signed/disconnected invariants, analytical viscous dipole travel,
physical unit scaling and resolvable near-support transfers.
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

## Prepared compact diffusion on fixed sites

`diffusion::sparse::PreparedPse::prepare_with_cx` binds the original particle
positions, positive quadrature areas and compact law into an immutable neighbor
graph. Radius-sized cells are sorted deterministically; bounded searches over a
two-cell halo identify candidates, and the original dense support metric admits
edges. Accepted unordered edges are sorted by original particle indices, so
forward/reverse evolution retains exactly the dense compact pair composition.
Gaussian laws refuse; no Gaussian tails are truncated.

`PreparationBudget` independently limits input particles, edges, actual
candidate visits, preparation work and live payload. Coordinates divided by
support radius must be finite and within +/-2^50; this keeps the rounding-safe
cell halo valid. Two passes count then reserve the exact edge inventory, and
cancellable in-place heap sorting requires no auxiliary sorting allocation.
Preparation uses O(N log N + C + E log E) work, where C counts candidates and E
is the actual edge count. Dense clusters can still require quadratic C and E;
the caps refuse that case rather than promising linear work for every cloud.
Peak logical preparation payload on 64-bit targets is 48N+16E bytes, with
24N+16E retained. Allocation metadata and capacity rounding remain excluded.

`PreparedPse::diffuse_with_cx` accepts arbitrary finite circulation on the
same bitwise-pinned geometry. It computes each half-step relaxation fraction
once and reuses those values for all steps. The existing equal-and-opposite
transfer, exceptional-range arithmetic and execution accountant remain the
owners. Active runs use 2ES pair visits, N+E+S(2E+1) logical work units and
N particle values plus E scalar fractions. There is no per-step allocation or
repeated pair-kernel evaluation. Identity runs skip fraction/pair work and
charge N+S. The immutable prepared graph is borrowed during a run and has its
own preparation receipt; moved/remeshed clouds require new preparation.

The remaining compact regressions compare independent exhaustive edge
membership, bitwise dense/sparse evolution, reused geometry after refusals,
cancellation, and 64/256-particle neighbor-work scaling.
`cargo run -p fs-vpm --example sparse_diffusion` prepares one Gaussian-vortex
cloud, reuses its graph for two successive runs, compares the heat field with
the analytical continuum solution and reports actual sparse/dense pair counts.
These counts establish algorithmic work, not measured machine speed. The
finite-cloud, support-resolution and splitting-error limits still apply.

### Adaptive time accuracy on the prepared graph

`PreparedPse::advance_adaptive_with_cx` completes an explicit duration using
one full symmetric pair step versus two half steps. It accepts the fine state
when every particle's local estimate
`|Gamma_fine-Gamma_coarse| / [3(A*atol + rtol*max(|Gamma_old|,|Gamma_fine|))]`
is at most one. The positive absolute tolerance is in vorticity [1/s]; the
relative tolerance is dimensionless. The division by three is the second-order
step-doubling estimate. No extrapolation, clipping, negative diffusion time or
circulation correction is applied, so accepted steps keep the original pair
maps' conservation and vorticity-range/enstrophy contraction.

Rejected trials halve the time step. Accepted ratios at most 1/8 allow
doubling, within explicit minimum/maximum steps; the last endpoint remainder
may be shorter than the minimum. A failed tolerance floor, stalled represented
time, aggregate work/pair/attempt cap, cancellation or original physical
refusal returns no partial window. Geometry remains pinned to the original
graph. Exact zero-time, zero-viscosity and edgeless identities still validate
the input and preserve its bits.

Three particle buffers and two fractions per supported edge are reserved once.
An active trial has six E pair visits and 2 N + 6 E + 1 logical work units;
initial validation costs N. Both fractions are cached while the exact trial
step is unchanged. Each changed step adds E logical units and two E actual
kernel evaluations, retained in `fraction_evaluations`. Rejected trials are
included in these counts. No trial allocates. The ambient accountant admits the explicit work ceiling
and reports actual consumption; inactive runs admit only N units. The final
result includes exact requested elapsed time, accepted/rejected counts and the
largest accepted local estimate.

This is automatic control of temporal splitting on a fixed spatial graph.
It is not a global time-error certificate or an estimate of finite-cloud,
quadrature, kernel, wall or continuum-model error. `tests/adaptive_diffusion.rs`
uses independent three-site diffusion eigenmodes, an unequal-area analytic
pair, tolerance refinement, conservation, identities and transactional refusal.

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

### Adaptive coupled motion and diffusion

`conservative::adaptive::simulate_viscous_with_cx` completes an explicit
duration while controlling both sampled vorticity and particle trajectories.
The original `D(h/2) A(h) D(h/2)` step is compared with two half steps;
only the fine positive-time state is accepted. Position error uses the
Euclidean norm divided by three times an explicit absolute metre tolerance.
Circulation uses the existing area-weighted absolute-vorticity/relative
estimate. A small circulation error cannot hide inaccurate particle motion.

For positive viscosity, each trial also compares its fine result with an
independent `A(h/2) D(h) A(h/2)` composition. The witness reuses the coarse
buffer after the ordinary comparison and calls the original positive pair
owners. Position and circulation differences use the same physical tolerances
without the step-doubling division by three. This detects initial transport
that rapid diffusion can extinguish in both ordinary candidates. All estimates
must pass, and only the original two-half-step state is accepted. These remain
local estimators, not a certificate for arbitrary stiff transients or global
trajectory error.

The existing viscous owner now accepts an internal positive time fraction,
so the full/half trials share the original pair maps, ordering, validation
and checkpoints. Keeping the fraction separate from the physical time avoids
erasing a representable subnormal interval before multiplication by a large
rate. The old fixed-step and scheduled callers still use exactly fraction one.
The same dyadic controller, robust scaled-error arithmetic, public budgets
and complete-window result are shared with fixed-site adaptive diffusion.

Each diffusion substage evaluates its actual moved geometry. Both Gaussian
and compact laws are supported; no pinned neighbor graph is used after motion.
With P=N(N-1), each positive-viscosity trial charges 12P pair visits, 7P
diffusion-kernel evaluations and 4N+12P+5 logical work units. Zero viscosity
skips the witness and retains 3P visits and 2N+3P+4 work units. Initial
validation costs N. Three particle buffers are reserved once; the witness
requires no additional storage or per-trial allocation.
Zero duration and fewer than two particles preserve input bits using one
particle buffer, while still validating all areas and physical controls.

Rejected trials count toward aggregate attempts/work/pair budgets. A tolerance
floor, invalid physical stage, cancellation or exhausted budget publishes no
partial trajectory. The accepted fine maps preserve the original circulation
and vorticity-range/enstrophy properties. Advection's impulse conservation
does not repair PSE's finite-cloud impulse error.

This is temporal error estimation for the original direct 2-D particle model,
not a global error certificate, spatial-remeshing controller, wall treatment,
variable-viscosity stress law or three-dimensional viscous solver.
`tests/adaptive_viscous.rs` checks the independent decaying-dipole
translation/strength solution including rapid initial viscous decay,
geometry-changing clouds against the original
refined owner, identities, resource counts and transactional refusal.

### Compact neighbors for adaptive moving diffusion

`conservative::adaptive::simulate_viscous_compact_with_cx` uses the same
adaptive full/half-step and reverse-splitting witness with a bounded compact
neighbor workspace. The original cell coordinates, two-cell halo, exact
compact-support predicate, and canonical pair order are shared with the
fixed-site graph owner. Every advection or candidate reset invalidates the
graph; only adjacent diffusion stages at unchanged positions reuse it.
Each supported edge evaluates its original fraction once per sweep and
reuses it for both forward/reverse conservative exchanges.

`CompactNeighborBudget` bounds live edges at every geometry and candidate
distance tests over the entire window. Three original particle buffers plus
N cell keys and at most min(max_edges,N(N-1)/2) edges/fractions are reserved
once. Graph probes, searches, sorts, candidates, fractions, exchanges and
discarded trials consume the original aggregate work account. Edge, memory,
pair, search and work limits are checked before exceeding their allowances;
failure returns no partial trajectory. Zero-duration, small-inventory and
zero-viscosity paths retain the original result and need no neighbor graph.

The dense Gaussian/compact path and existing fixed-site preparation remain
available. No Gaussian tail truncation, new diffusion coefficient, or change
to the accepted adaptive composition is introduced. Diffusion visits depend
on supported edges; direct advection still requires 5N(N-1) visits per active
viscous trial, and graph construction has its own charged cost. Subquadratic
whole-solver scaling and wall-clock speedups are not claimed.

Four public-API regression groups cover bitwise dense-path parity on separated
clouds, particles entering support during motion and rejected trials, exact
aggregate caps and late refusal/replay, and identity/zero-viscosity behavior.
These groups pass the bounded native suite together with the dense and
source-coupled adaptive owners.

### Adaptive viscous motion in a held affine source

`conservative::adaptive::simulate_viscous_with_affine_flow_cx` joins the held
incompressible affine owner to the same adaptive viscous clock and positive
PSE/advection stages. A full trial is
`F(h/2) D(h/2) A(h) D(h/2) F(h/2)`; its accepted candidate is two original
half steps. The adjacent fine quarter-source maps coalesce into one half map.
Both circulation/vorticity and Euclidean position estimates must pass. The
positive-viscosity reverse-order witness remains
`F(h/2) A(h/2) D(h) A(h/2) F(h/2)`, preserving the independent check for
initial induced motion erased by rapid diffusion. No negative diffusion time,
extrapolation, moment projection or frozen post-motion neighbor graph enters.

`AmbientViscousAdaptiveRun` separates accepted affine moment increments from
accepted PSE increments. The former reuse `AmbientMomentChange`. The latter
report actual changes in `sum Gamma*x`, `sum Gamma*y` and `sum Gamma*|x|^2`,
plus circulation/first/quadratic representation defects relative to each
original signed pair transfer. A finite irregular PSE cloud generally changes
vortex impulse even in exact arithmetic; that change is not mislabeled as an
external source or floating-point loss. The pair observer derives commanded
increments from relative-position products and defects from local before/after
circulations, using compensated accumulation. It does not subtract large
whole-cloud moments. Lost sub-ULP transfers therefore appear as representation
defects without inventing an actual moment change. Original pair advection
retains its separate roundoff-level impulse defect.

Only accepted fine paths contribute to the physical accounts. Rejected,
coarse and stiffness-witness work contributes to execution receipts and caps.
With P=N(N-1), an active field and positive viscosity with N>=2 use 12P pair
visits, 7P diffusion-kernel evaluations, 7N source visits and 12P+11N+5 total
work units per trial. Zero viscosity or a singleton uses 3P pair visits,
no diffusion evaluations, 5N source visits and 3P+7N+4 work units. Initial
validation costs N. Zero fields add no source visits and retain the native
viscous trajectory, adaptive decisions and work counts exactly.

The existing `AmbientAdaptiveBudget` caps source visits over the entire
window. All attempts/pairs/source/work limits are checked before each complete
trial. The same three particle banks are reserved once; the witness reuses
the coarse bank. Generic observer calls occur inside the original metered PSE
pair unit with constant stack storage, and the unobserved path retains the
original pair arithmetic and order. There are no per-trial allocations or
extra moment-scanning cloud passes. Source maps are prepared before candidate
copy work. Active-source windows share the original inviscid compensated time
helper, preventing accumulated-time roundoff from inventing a tiny final map;
the zero-source viscous clock remains unchanged. Empty/zero-duration paths
remain identities, while a held source transports a singleton. All borrowed
inputs remain unchanged on invalid data, numerical/moment representation,
tolerance, execution-budget or cancellation refusal; no partial cloud escapes.

These are temporal error estimators for the stated direct finite-cloud model,
not a global error certificate or a guarantee for every stiff source transient.
Gaussian/compact spatial quadrature still needs overlap, domain padding and
refinement. Fixed incompressible particle areas, held source parameters and
constant viscosity remain assumptions. No finite carrier recoil, solved wall
reaction, viscous heating or finite-volume mechanical-work account is claimed.
`tests/adaptive_affine_viscous.rs` independently integrates the simultaneous
nonlinear position/strength ODE and six source/PSE moment derivatives for both
kernels. It also checks a stiff rotating decaying dipole, exact native
zero-source parity, original fine-map composition, frame/unit covariance,
singleton/subnormal motion, sub-ULP transfer loss, whole-window caps and
cancellation inside sources, diffusion and final publication (G0/G1/G3/G4).

`simulate_viscous_compact_with_affine_flow_cx` selects that same moving compact
backend in the same adaptive source owner. Every active affine stage also
invalidates cached geometry; source-driven entry into support cannot be missed
between the adjacent fine diffusion stages. The original per-edge fractions,
canonical ordering and signed-transfer observer remain shared. The dense and
compact physical states and accepted source/PSE accounts match bit for bit.
The report additionally exposes all graph rebuilds and candidate distance
tests. In an active positive-viscosity source trial, seven geometry rebuilds
replace the unforced compact path's six; direct advection still costs 5P visits.
If K is the sum of supported edges over its seven diffusion stages and G the
charged graph work, the trial costs 5P+2K pair visits, K kernel evaluations,
7N source visits and 5P+3K+G+11N+5 work units. These costs include rejected
trials and witnesses. The complete graph payload is reserved before evolution;
data-dependent work is admitted before each bounded operation. A zero source
retains the original compact driver's complete physical/execution result.
The focused source/compact tests check strain-induced support entry, exact
dense state/account parity, rejected-trial counts, graph/source cap exhaustion
and deterministic reuse after a late edge or work refusal.


## Prescribed affine flow on conservative 2-D vortices

`conservative::ambient::{simulate_with_cx, simulate_viscous_with_cx}` adds one
held `AffineFlow2`, `u(x)=U+A(x-o)`, with `A=[[a,b],[c,-a]]`. Its three authored
gradient entries make incompressibility structural. Translation, rigid
rotation, simple or rotated shear, and planar extension use an exact affine
half-map around the original native step. The viscous composition is
`F(dt/2) D(dt/2) A(dt) D(dt/2) F(dt/2)`; here the central `A(dt)` denotes
original conservative pair advection. All particle areas remain explicit and
unchanged under this area-preserving source. Negative time is supported only
by the inviscid entry.

The half-map is prepared once, using Cayley–Hamilton, the original deterministic
elementary functions and double-double arithmetic. Entire-series evaluation
handles the zero-eigenvalue/shear limit without division by a small eigenvalue;
spectral projectors preserve expanding and contracting strain components.
Unrepresentable coefficients, endpoints or moment accounts refuse atomically.
An exact affine stagnation point is retained; a nonzero endpoint displacement
that disappears completely at the stored position refuses rather than silently
turning the source off. No step subdivision, position snapping, core deformation
or arbitrary spatially varying source is inferred.

`AmbientVpmRun` retains the original `VpmRun` plus actual source increments in
`sum Gamma*x`, `sum Gamma*y` and `sum Gamma*|x|²`, measured during the source
substeps with compensated endpoint differences. No division by total
circulation is used, so signed zero-net clouds remain valid. For inviscid motion
the first moment obeys `M'=A(M-Gamma_total*o)+Gamma_total*U`; the source map
integrates that balance exactly up to roundoff, independently of pair splitting.
The quadratic moment can change under strain. The source accounts for those
changes while pair advection retains its original impulse invariants. PSE still
has its separate finite-cloud impulse error. These vortex moments do not claim
finite-domain fluid/body momentum, reaction force, pressure work or a closed
energy budget for a background extending to infinity.

An active source adds exactly `2*N*S` logical visits to the original plan and
no heap payload beyond its one particle vector. Map preparation and receipts
use fixed stack storage. All phases share one original ambient account and
bounded work checkpoints. Zero source, zero time, zero steps and empty-cloud
paths preserve the original run and resource counts; within nonnegative time,
zero viscosity also recovers the new inviscid result exactly. Fixed-source
second-order time refinement is checked against an independent full ODE,
alongside analytic affine motions, signed moment balance, reversal, world-frame
and dimensional covariance, resource admission and atomic refusals in
`tests/ambient_flow.rs` (G0/G1/G3/G4). This extension does not close the deferred
3-D/FMM/BEM program or provide automatic remeshing under strain.

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

## Scheduled remeshed evolution

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

`conservative::remeshed::simulate_with_affine_flows_cx` composes the same
trajectory with exactly `ceil(S/K)` authored `AffineFlow2` sources. Source `k`
is held after remapping to grid `k` until remapping to grid `k+1`; the last
partial interval has its own source. Each physical step uses the original
`F(dt/2) D(dt/2) A(dt) D(dt/2) F(dt/2)` owners. Flow changes occur only at
these authored interval boundaries, without inferred time interpolation.
Zero steps requires no sources. Zero source or zero time preserves the original
scheduled result and resource receipt exactly, including all authored remaps.

Each active interval adds `2*S_k*G_k` source visits to the complete work plan,
including zero-strength quadrature nodes. The unchanged two-bank payload bound
still applies. One fixed-size affine map is prepared on entering each interval;
there are no per-step cloud copies, map allocations or further heap payload.
Source, PSE, pair advection and remapping share one execution account, including
source checkpoints and complete work admission before any allocation.
`AmbientRemeshedVpmRun` returns the original regular-cloud receipt plus actual
first/quadratic source moment increments accumulated in double-double arithmetic
across all intervals. PSE and remap defects remain separate from that source
account. This does not imply a fluid/body reaction force or close the momentum
balance of an infinite prescribed background. Arbitrary changing source/remap
schedules do not acquire a global second-order convergence guarantee.

G0/G3/G4 checks in `tests/remeshed_affine_flow.rs` compare the complete run with
explicit original-owner composition for positive and zero viscosity, moving
windows and changing anisotropic resolution. Independent piecewise affine moment
laws cover signed zero-net circulation, while whole-source work, exact combined
payload, inactive-source parity and late source/support/cancellation refusals
exercise the shared transaction and deterministic retry.

### Selected-order inviscid trajectories

`conservative::remeshed::simulate_with_order_cx` and
`simulate_with_affine_flows_order_cx` accept `CompositionOrder` in place of a
diffusion law. They use the same scheduled-remesh transaction, original two
particle/area banks, exact pair owner, signed composition coefficients and
affine source maps. `Second` preserves the original zero-viscosity particles,
moment changes, work/memory receipts and checkpoint sequence exactly. Existing
viscous entry points continue to call the unchanged D/2-A-D/2 owner.

`Fourth` completes all three signed pair sweeps within each positive physical
step. For an active held field, adjacent source maps combine into the same
four passes as the standalone affine owner; both distinct weighted maps are
prepared once on entering the interval. There are no added cloud buffers,
per-step maps, allocations or dependencies. The original M4' remap runs once
after each complete authored interval, including the final partial interval.
Neither a negative substep nor a partial composition can trigger remapping,
change the held source or invoke irreversible PSE diffusion.

For interval `k`, Fourth admits `3*S_k*G_k*(G_k-1)` pair visits and, for an
active field, `4*S_k*G_k` source visits. Total work remains the exact sum of
those visits, the original remap work `R` and `S` physical step units. The
largest target still determines the original two-bank payload cap. Inactive
fields or zero time add no source visits; zero time still performs all authored
remaps and zero steps only the initial remap. Signed source moment increments
accumulate in the original double-double account over all intervals; remap
interpolation error remains distinct. Whole-run budget admission, bounded
cancellation and late failure with no partial publication cover every signed
stage and endpoint remap.

Fourth-order temporal convergence applies to the reversible evolution between
fixed remap endpoints with each field held. Changing spatial resolution or the
number/timing of remaps changes the discretization separately; an arbitrary
changing schedule has no automatic fourth-order global claim. This extension
does not infer adaptive timesteps or meshes, a boundary reaction, 3-D diffusion
or a closed-system momentum account for prescribed ambient flow.

`tests/remeshed_order.rs` (G0/G1/G3/G4) compares changing-grid trajectories with
explicit existing-owner composition, checks signed zero-net impulse and an
independent piecewise affine first-moment law, and demonstrates fourth-order
time refinement with fixed initial/final remaps against an independently
refined full-ODE RK4 reference. Equal-pair-work comparisons check the practical
accuracy benefit. Tests also cover exact original/inactive receipts, changing
quadrature, complete signed-stage work and two-bank admission, late source-map
and support refusal, and cancellation inside a coalesced source bridge, the
negative-time pair stage and the final remap, followed by deterministic retry.

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

## Complete physical wake impulse and interval exchange

`FarField::hydrodynamic_impulse_per_density(near, origin_m, max_sources,
cancelled)` returns a `WakeImpulse` containing all three linear and angular
geometric impulse components of the supplied closed physical filament set.
It uses the near sheet's physical shared edges, including its closing lifting
line, and every retained exact far segment. Core growth, multipole caches and
image devices do not enter these moments. Signed endpoint circulation must
close within roundoff, using the same deterministic checker as hybrid
advection. An open pruned or mismatched interface refuses; removal of an entire
closed component requires separate caller accounting because closure alone
cannot detect that missing history.

A nearby deterministic anchor and compensated sums retain both moments before
origin transport. `about_origin` re-expresses angular impulse without revisiting
sources, preserving `J(O+a)=J(O)-a cross I`. `exchange_since` first differences
the local compensated moments and then transports the difference, so a small
change need not disappear between two large rounded world-origin readouts.
Both endpoints use the later snapshot's origin as one fixed spatial reference;
changing only a query origin generates no exchange.

`WakeImpulseExchange` reports positive wake linear/angular impulse changes in
N s and N m s for one declared constant uniform density, together with their
negative mean force/moment equivalents over a positive elapsed interval. These
are contributions from the represented wake, not a complete aircraft load
model or arbitrary control-volume momentum. No load is automatically applied.
Omitted body/bound vorticity, potential/added-mass impulse, imposed flow,
supports, outer boundaries and numerical conversion/pruning changes retain
their separate physical owners.

The source cap precedes the near-geometry scan and bounds conservative
near-edge capacity plus retained far records; far-cell count is separately
capped, including empty cells. Closure uses O(N) scratch and O(N log N)
indexing; geometric accumulation is linear. Source sweeps checkpoint every
256 visits and before returning, and typed numerical/cap/cancellation refusals
return no partial diagnostic. `tests/hybrid_impulse3d.rs` covers independent
closed-edge integrals, partition equality, translation/rotation and origin
transport, pure angular impulse, a small exchange hidden by a distant origin,
supported-ring dynamics, density/time/sign conventions, image/core/cache
exclusion, source bounds and refusal/cancellation.

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

### Error-controlled tree dynamics

Both near and combined near/far systems expose `shed_self_induced_with_tree`.
The caller supplies `WakeTreeSpec`, with an absolute stage-velocity error
budget in m/s and an exact-leaf capacity. Zero tolerance evaluates the original
source order exactly. The default existing entry points retain direct
summation. Freeze the tree settings per run and include `WakeTreeSpec::digest`
in its numerical identity; these settings are not a runtime autotuner.

Each velocity stage builds a deterministic balanced source tree from the same
prepared continuous finite segments, including optional image segments. A total
midpoint-coordinate/source-index ordering fixes subdivision. Nodes retain the
shared Cartesian monopole, dipole, continuous-line quadrupole and leading
core-variation moments. They do not clone physical segment lists. Scratch
storage is reused across stages. Free near/far vertices still evolve once,
source circulation and age/core semantics remain unchanged, and final source
admission precedes atomic publication.

An internal node is accepted only when its full expansion error bound fits its
share of the absolute budget, allocated by unsigned circulation-length measure.
A half-budget margin accommodates accumulation/partition roundoff. Near or
uncertifiable nodes open to smaller cells; leaves use the original exact
finite-segment kernel. A cheap squared-distance rejection gate inverts a lower
bound on the required Taylor error bound, so it can save bound evaluations but
can never accept an uncertified approximation. Images use the same budget;
their approximate normal-velocity residual is bounded by the local error budget,
while physical stage/endpoint nodes still obey strict plane clearance checks.

`WakeAdvectionReport` reports exact pair evaluations, accepted/visited cells,
peak node count, metered work and the largest summed local velocity error
bound. `SelfInducedShed::max_work` meters this mode as follows: source assembly
uses the conservative source bound plus one per stage, construction charges
bounded source passes and sorting, each target charges one, each tree visit
eight, each accepted expansion sixteen, and each exact pair one. Construction
and source/target work are included even when all circulations vanish. Exhaustion
or cancellation rolls back the entire interval; accelerated admission does not
require reserving the worst-case direct all-pairs workload.

The report controls local source-evaluation error, with the analytic and
floating-point qualifications of the expansion bound below. It does not
certify accumulated trajectory, impulse or energy error, impose mutual
pairwise approximate forces, or certify a real-time aircraft schedule. The
whole-step battery checks exact zero-tolerance parity, accuracy refinement,
nonzero accepted-cell work, image/core evolution, bit-identical physical
partitioning, deterministic measured work caps and transactional cancellation.

### Coupled material motion across near/far representations

`FarField::shed_self_induced` and its environment form advance retained exact
far sources together with the attached near sheet. Near and far must share
the same shedding interval and core mode. Signed circulation must balance at
every combined physical endpoint; incomplete/pruned open vortex lines refuse.
A deterministic coordinate index binds shared free endpoints once per interval
(signed zeros denote the same point). Prescribed line nodes have separate
identities even when a newborn free row occupies the same coordinates.

The common stage integrator updates each unique free vertex once. Exact far
sources precede near sources in the direct sweep, preserving the physical
ordering across aggregation. Both representations see the same current and
midpoint geometry, incident flow, line trajectory, images and elapsed diffusion.
Far source circulations and stored base cores remain unchanged. After all
stages succeed, each cell's geometry and multipole moments are rebuilt from
its moved exact sources once; the accumulated diffusion offset remains separate
from those stored cores. Near and far states publish together only after the
last cancellation check. Separately stored image snapshots must be rebuilt
after successful physical motion.

This is a bounded direct reference path, with cost proportional to unique
nodes × retained sources × substeps. It avoids repeatedly convecting shared
endpoints and preparing segment geometry per target, but does not accelerate
the dynamics using far multipoles. Source assembly is a linear pass between
substep checkpoints; direct target sweeps check cancellation within every 256
source contributions. Circulation transport is conservative; exact impulse or
energy conservation is not asserted for a discrete wake interacting with its
prescribed lifting line or image boundary. Core spreading remains the declared
algebraic model rather than a resolved viscous vortex tube. Streamwise
refinement is an explicit separate operation.

`tests/hybrid_advection3d.rs` checks the independent translating-square-ring
solution and time-dependent diffusion integral, second-order convergence,
bit-identical frozen-core partitioned/unpartitioned motion over multiple
intervals, rebuilt cell geometry and core moments, shared seam continuity,
environment equivalence, and atomic topology/time/core/work/cancellation
refusals. The tests do not infer an accurate far expansion inside its source
cell from the exact direct evolution used here.

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

`FarCell::exact_kinematics`, `FarField::exact_kinematics` and
`hybrid_kinematics` extend the same analytic velocity/Jacobian through retained
far sources. Every segment keeps its individual source core plus continued
diffusion offset; the near part uses the far owner's admitted core mode,
shedding interval and retained row ages. The complete hybrid sheet includes
its closing bound line, so a consumer must not add that line again.
`ImageSystem::exact_kinematics` and `grounded_hybrid_kinematics` evaluate the
same physical sources and mirrored boundary devices. Velocity is a polar
vector, vorticity is axial, and the reflected Jacobian is `R grad(u) R`.
Images never enter physical circulation or impulse ledgers. Image construction
checks both near and far clearance and finite reflected coordinates; a snapshot
must be rebuilt after source geometry, age/core or mode changes.

These queries allocate no sampling buffers and cost O(retained sources) per
point. They use exact segment derivatives throughout, not derivatives of
velocity-only multipole expansions or their switching conditions. There is no
derivative-error claim for approximated or pruned sources. Invalid public source
shapes, non-finite data/results, incompatible shedding intervals and singular
active zero-core encounters refuse; a zero-strength segment has no singular
physical contribution. `FarField::validate_kinematics_sources` exposes the
successful allocation-free source preflight for downstream batch consumers.
It does not establish nodal closure or completeness after pruning, and exact
sampling deliberately ignores cached expansion moments. G0/G1/G3 regressions
in `tests/hybrid_kinematics.rs` check mixed-core partition invariance, continued
diffusion, independent velocity finite differences, axial reflection, plane
derivatives, finite-core axis curl, malformed source refusal and far clearance.

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

## Absolute three-dimensional expansion error

`FarCell::expansion_error_bound(p)` bounds the Euclidean velocity error of the
stored Cartesian Taylor polynomial, in m/s. It returns infinity inside the
existing geometry/core fallback domains, for invalid caches or probes, and
when the polynomial's arithmetic cannot be represented safely. The existing
`eval` and source-moment construction retain their arithmetic and fallback
behavior. Both now share a fixed-size expansion cache with no owned source
vectors; its source-iterator constructor includes the exact continuous-line
quadrupole `dl_l * dl_m / 12`.

Write `R=|p-centroid|`, `a=radius`, `d=R-a>0`, `c=mean core²`,
`Delta=max|source core²-c|`, and `W=sum |Gamma|*|dl|`. For
`F(x,c)=x/(|x|²+c)^(3/2)`, the retained polynomial is
`F - DxF*delta + DxxF[delta,delta]/2 + DcF*eta`. Its full real-arithmetic
remainder is bounded by

```text
W/(4*pi) * [(4*a³ + 6*a*Delta)/(d²+c)^(5/2)
           + (15/8)*Delta²/(d²+max(0,c-Delta))³].
```

The terms cover spatial order three, omitted space/core coupling, and core
order two, respectively. The implementation tightens the first two terms by
also integrating their derivative bounds along the Taylor path and taking the
smaller independently valid bound. With `rho=sqrt(R²+c)` and `theta=a/rho`,
the integrated spatial factor is
`theta³*(4-3*theta)/(rho²*(1-theta)²)` and the mixed factor is
`(3/2)*Delta/rho⁴*((1-theta)^(-4)-1)`. Cancellation cannot reduce `W`.

A separate conservative floating-point allowance includes source-count
accumulation and world-coordinate conditioning. This is not a directed-rounding
enclosure, a temporal integration-error certificate, or an audit of arbitrary
external edits to public cached moments. Continued common core spreading is
included through the evolved mean core and unchanged core deviations.
`tests/expansion_bound3d.rs` checks heterogeneous cores, continued diffusion,
signed cancellation, an independent analytic long-line oracle, scale and
far-distance behavior, and mandatory exact-fallback domains.
