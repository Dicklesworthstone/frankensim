<!-- fs-lbm contract, part 1 of 3; a higher number is newer. Index: ../CONTRACT.md. -->
## Rigid wall momentum on the original thermal VOF trajectory

The immutable thermal rigid-wall candidate methods use the existing link-resolved rigid reflection with current cellwise material relaxation and the original Guo force frame. That same native tick records the actual mass routes consumed by extensive-enthalpy transport. Final source/rheology, wetting, wall mass, population/Mach, material speed, SI conversion and cancellation checks precede return.

A private evolved-state tick separates preparation from publication without cloning prepared source or electrical operators. Each thermal rigid-wall trial clones the original native fluid once and retains the same material/caloric bindings and physical clock. The existing electrical/conduction ownership gate remains operative. Optional molecular heating reuses the original held initial-state stress-power source; actuator work and numerical loss are not additional caloric inputs. Fixed flooded-plane and single-phase limits remain. Focused regression definitions exercise isothermal parity, heat and buoyancy, actual rigid replay, original molecular source reuse, and refusal/retry. Formatting, parsing and scoped whitespace checks passed; native compilation and execution remain outstanding.


### Finite guided plane sliders on the physical D2Q9 disk

The native dry-disk owner now accepts FiniteLubricationSlider2 through
DryDisk2::step_with_lubrication_sliders. Each source supplies a positive
physical mass per out-of-plane depth, a distinct center of mass, a signed
tangential velocity and the original stationary plane/film geometry.
The guide fixes normal position and orientation. Infinite-plane geometry
therefore stays exactly unchanged under travel; no finite face extent,
rotary inertia, resolved-fluid wall or free normal recoil is inferred.

Every original three-unknown dry-fluid candidate condenses the finite
slider midpoint with k_eff = k/(1+k/(2m)) and
u_free = u_before + J_authored/(2m). The stiff branch avoids overflow of
k/(2m). The disk receives the film load inside its existing native solve;
the slider receives the opposite tangent impulse and actual kinetic work.
Normal traction and the complete distributed-traction couple go to a
separate ideal guide port, whose work is zero. The original cylinder spin
resistance is retained; translating recoil does not silently remove it.

Slider centers undergo both free half-drifts around the same impulse.
Receipts separately retain authored work, true endpoint momentum/kinetic
change, guide impulse/couple and their measured balance. The combined
exterior-fluid/dry-disk/slider momentum account includes guide reactions
as external inputs. Every slider endpoint, work and momentum check is
prepared before the original sole body/fluid publication, followed only
by infallible assignments of the already-admitted slider states.

Focused native-source checks cover simultaneous oblique sliders with
Guo acceleration and authored loads, independent fixed-origin momentum
and kinetic-work reconstruction, stiff recoil, original inactive/empty
bit parity, late refusal and cancellation with exact retry. The bounded
stiff limit is also checked without forming the overflowing raw ratio.
Implicit midpoint remains passive but is not L-stable: unresolved stiff
films may reverse relative velocity during a tick, so physical relaxation
still needs time resolution. Existing model cutoffs, positive-gap
admission, 64-plane limit, no contact and no film-heat destination remain.

## Native curved D3Q19 wall geometry and momentum

`CurvedWall3` binds canonically ordered original wall links to their actual
fluid-side fractions and caches the opposite fluid donor for fractions below
one half. `from_level_set` checks all native occupancy signs and brackets every
voxel-wall intersection with explicit root tolerance/iteration limits. The
callback's continuity and single crossing on each link remain caller geometry
obligations; no SDF/topology certificate is inferred. An unresolved short-link
gap, stale donor, mismatched mask or exhausted root budget refuses.

`BoundaryGrid3::step_with_curved_walls` uses the original four collision laws,
Guo source frame and pull table, then applies linear Bouzidi-Firdaouss-Lallemand
interpolation only at selected links. Surface velocity is sampled once per
link at its actual intersection. Exactly-halfway geometry retains the original
population and wall-exchange arithmetic. Sparse geometry/cursors and reusable
post/next fields avoid a dense distance/velocity field or fresh per-step
population allocation after warmup.

The SAME incoming population is used for both streaming and the shared wall
momentum convention. Actual population impulse, boundary-relative impulse and
velocity-times-mass correction remain separate; torque uses the off-lattice
intersection, and work uses that same wall velocity. Unselected walls, held
fluid forcing and native periodic angular transport are separately counted.
Interpolation mass drift is measured, never renormalized into a false swept
volume. Every candidate population and load passes finite-state/cancellation
checks before one final swap. Late refusal publishes neither fluid nor loads.

This is fixed-occupancy curved geometry with prescribed surface kinematics,
including tangential rotation of invariant spheres/cylinders. It neither
advances a solid pose nor remaps newly covered/uncovered fluid. Existing
moving-body owners remain authoritative for evolving immersed geometry.
Regularized reservoir faces have their own measured exchange below. There is no general
GCL, full-solver Galilean-invariance, calibrated body-drag, complete fluid-energy
or global second-order convergence claim.

Six native regressions pass, covering original halfway bit compatibility and
global momentum accounts for BGK/TRT/central-moment/reduced-cumulant collision,
both interpolation branches, actual generated sphere torque/work, periodic
transport, malformed/missing donors and late cancellation. Fractional planar
Couette recovers its physical-wall linear profile within 2e-10. Resolved
cylindrical Couette relative velocity errors are 0.01472427871037479 on 24x24x4
and 0.007795637138649035 on 48x48x4, under respective 2%/1% budgets; this measured
pair does not establish asymptotic order. A 2.3-cell inner radius is outside
that accuracy envelope (about 12.4% error). The ten existing native finite
wall-rotor regressions also pass after sharing the incoming-population momentum
helper. These are direct original numerical-owner checks with only unused
runtime/pool and unrelated module wrappers omitted; full workspace Cargo,
DSR, cross-ISA and performance gates were not executed.

### Local material viscosity on the original curved-wall step

`BoundaryGrid3::step_with_curved_walls_and_relaxation` accepts a step-local
x-fastest shear-relaxation field. It uses the SAME native collision, held Guo
force, interpolated pull populations and wall/periodic momentum accounts as
the uniform curved-wall entry. BGK/TRT/central-moment/reduced-cumulant families
retain their original non-shear controls, including TRT's magic product and
D3Q19's existing bulk/shear tie. Unchanged rates retain original bits; a later
ordinary step still uses its constructor model. Solid padding is ignored.
Invalid field dimensions/rates or native endpoint/load/cancellation refusal
leave accepted populations and topology unchanged.

`material::MaterialScaling::step_curved_walls_with_material` samples
each original `FluidTransport` response once, retaining that concrete value
and its exact applied tau. Borrowed source values permit deliberately frozen
fields without repeating property interpolation. The existing momentum
scaling owner converts nu=mu/rho0 using the reference dx/dt. Every local source
must retain the same physical reference density; mass, inertia and fluid
inventories are not remapped. Cp/k remain source-admitted but impose no
thermal-lattice restriction on this momentum-only operation. The explicit
cell budget precedes sampling/allocation. All source and conversion checks
precede native stepping, with no fallible publication after native acceptance.

Temperature/composition evaluation, source epochs and any thermal evolution
remain with the caller. The method does not reconstruct source pressure,
solve implicit shear rheology, deposit viscous heat or substitute actuator
work for heat. Geometry, wall velocity, moment origin, held force and returned
momentum/work accounts retain the native lattice-unit convention. The caller
must use the pinned density/dx/dt mapping and maintain its scenario speed
bound. Curved-wall interpolation mass drift and fixed-occupancy restrictions
are unchanged.

Focused regressions exercise constant-source population/report identity for
all four collision families, independent local collision/pull reconstruction,
two-temperature piecewise Couette fluidity and actual wall work, original
source retention and single evaluation, and budget/source/density/rate/late
wall/cancellation rollback. These new regressions await native execution;
no new measured convergence or performance claim is made here.

### Source shear-dependent viscosity at curved moving walls

`MaterialScaling::step_curved_walls_with_material_rheology` resolves each
incoming fluid cell through the original Guo-corrected D3Q19 stress estimator
and shared safeguarded `RheologySolve`. The held physical force supplies
`A = gamma*tau`; the same source query closes
`tau = 1/2 + 3 (mu(gamma)/rho0) dt/dx^2` with `gamma = A/(tau*dt)`.
The constructor rate is only an initial iterate. All six symmetric stress
components participate, and a held acceleration is not inserted as strain.

The complete stress field is prepared once before source work. Each source
response is read once per bounded constitutive evaluation; its accepted
concrete value, physical query, iteration count and tau defect remain in
`MaterialBoundaryRheology`. Cached transport then passes through the SAME
material boundary sampler and native curved step. No second collision,
streaming, wall law or publication path is introduced. BGK/TRT/central-moment/
reduced-cumulant non-shear controls and open-reservoir accounts remain intact.

Every trial keeps the original reference density and momentum scaling.
Temperature, composition and claim selection stay with the original source
producer. The callback must not advance external state. Source-domain refusal,
exhausted root/cell budgets, invalid native stress or late cancellation preserve
accepted flow; the original sole population swap remains final. No arbitrary
viscosity floor/cap or unsupported source extrapolation is supplied.

This is low-Mach incoming-time generalized-Newtonian constitutive splitting,
not implicit time integration or a thermal, yield-history or viscoelastic
model. Next-time material admission belongs to the next call. Physical-rate
and SI-source queries coexist with the original lattice-unit wall/momentum
accounts; fixed-occupancy geometry, interpolation mass drift and the native
bulk/shear tie are unchanged.

Four focused regressions cover constant-source identity for every collision
family, independent full Guo-frame stress and a two-axis material-law quadratic,
Carreau Couette velocity/wall work, and source/density/root/cancellation refusal
with exact retry. Product-source parsing passed before the executor outage;
the reconstructed final payload and new tests remain uncompiled/unexecuted,
and the tests have not passed a parser check. No measured nonlinear-flow
accuracy, performance or arbitrary constitutive convergence claim is made.

## Constitutive shear and volumetric strain

Molecular shear-viscosity queries use `gamma=sqrt(2 dev(D):dev(D))`, with
`D=sym(grad u)`. Three-dimensional gradient, population-stress, sparse,
immersed, free-surface and tetrahedral-transfer paths use the same invariant
as `fs-material::DeformationRate`. Isotropic dilation and rigid rotation have
zero shear rate; extension retains its deviatoric contribution. D2Q9 retains
its existing dimensional strain and bulk conventions, described below.
Independent bulk collision and explicitly sourced bulk heat retain the trace.
BGK/TRT liquid heating retains its tied `zeta=2mu/3` contribution. Source
viscosity observes deviatoric shear while `2mu D:D` includes volume-change work;
the thermal feedback cache retains both strain components in one stress pass.
Native MAC rate observation removes spherical normal strain before squaring;
its existing full-stress momentum/power operator retains its original response.
Analytic regressions cover dilation, extension, shear, rotation, force-frame
changes, source viscosity, and the separate nonzero heat of pure dilation.

## Source rheology on moving immersed bodies

`MaterialScaling::step_moving_body_with_material_rheology` and
`step_moving_bodies_with_material_rheology` bind original
`FluidTransport` sources to the native moving partial-saturation owner.
The existing physical half-force velocity and wall/open/periodic gradient
use the shared material deformation-rate kernel to supply
`sqrt(2 dev(D):dev(D)) / dt`. Rigid rotation and isotropic dilation contribute
no material shear. Each non-solid cell resolves one original source and
reads its transport once, independently of body count. A mixed solid/fluid
cell is sampled through this velocity gradient; its blended nonequilibrium
stress is not treated as a pure-fluid constitutive estimator.

The source's dynamic viscosity and retained reference density determine the
same local relaxation used by fluid collision, the Noble--Torczynski weight,
and the implicit six-velocity rigid-body response. Multiple bodies sharing
disjoint subvolumes of a cell use one total-occupancy denominator. BGK, TRT,
central-moment and reduced-cumulant non-shear controls remain native. An
unchanged relaxation preserves the authored collision model exactly. A
one-body assembly delegates to the same single-body coupling.

The shared source field is prepared once before native rigid half-drift,
geometry sampling, simultaneous exchange and final half-drift. Existing
pair, wall/halo, endpoint-force and cancellation checks precede the original
sole population publication; no second fluid step or sequential body update
is introduced. Source, allocation, density, geometry or cancellation refusal
preserves accepted populations, rigid states and force history. Scratch
workspaces may change. Explicit cell and kinematics byte budgets bound the
preparation; no fallible source work follows native publication.

Sources retain the original rho0/dx/dt conversion and exact reference density.
Temperature, composition and material history remain with the caller, whose
callback must not publish tentative history. Cp and conductivity impose no
thermal lattice constraint on this momentum-only path. Coefficients are
frozen at the incoming time level and resampled on the next call. This is
not an implicit endpoint rheology solve, viscous-heat integrator or subgrid
lubrication law. Native lattice units, low-Mach restrictions and additional
body inertia over retained interior fluid remain unchanged.

Focused regressions cover all-family constant-source/native identity,
independent Carreau collision/body response, all deviatoric strain components,
rigid rotation and dilation invariance of source viscosity and coupling weights,
source and geometry refusal, late cancellation with exact retry, shared-cell
simultaneous momentum exchange and one-body-array identity. These sources
passed Rust parsing and scoped whitespace checks; native compilation and
execution remain pending. No measured accuracy or performance result is
claimed for this extension.

## Source-backed gas convection on moving-liquid interfaces

`MaterialFreeSurfaceEnthalpy3::step_with_surface_convection` and
`step_with_finite_gas_surface_convection` evaluate a call-local
`FreeSurfaceGasConvection` inside the original zero-storage surface solve.
Each trial receives its actual interface temperature, current original gas
bulk state, and the same selected gas source claims at the mean interface/bulk
temperature and current pressure/composition. A positive finite coefficient
replaces only `Gg = kg A/dg` with `Gg = h A`; it does not add another gas film or
remove the liquid half-cell resistance or gray radiation. The exact coefficient
evidence and original film source remain in `FreeSurfaceInterfaceSample`.
`None` is accepted only at exact interface/gas thermal equilibrium. A failed
correlation, unsupported source, invalid coefficient or unrepresentable
conductance refuses the complete transaction.

The same borrowed provider remains attached to every cloned thermal observer,
nonlinear trial and endpoint evaluation. Existing positive pairwise surface
elimination supplies the signed liquid/gas/radiant heat and the explicit
capacity-based interval bound. The finite atmosphere stores exactly the
opposite gas heat on each original thermal substep; its new temperature, EOS
pressure and source properties feed subsequent trials. Optional pressure loading
uses this same parcel through the existing native pressure owner. Accepted
liquid, wall and gas states still publish together. The opt-in implicit residual
calls the provider at each trial endpoint for prescribed or single finite-gas
surfaces. Finite solid cells and contacts participate through their original
thermal owner; electrical and separate liquid environmental paths retain their
implicit exclusions.

This entry admits the original Cartesian, height-graph and local PLIC surface
geometry. Its provider receives the actual unit outward normal, physical area
and current patch owner. A higher owner must still admit its physical surface
model, characteristic length and held relative gas flow; a grid patch is not a
new boundary-layer entrance. The existing plate adapter is
`fs-airflow::graph::thermal::medium::convection::PlanarFreeSurfaceConvection`.
It admits a planar surface in world coordinates, including reconstructed planar
geometry, but does not supply an empirical law for a curved interface.
The held-geometry thermal/native-flow split, isothermal-patch and held-property
assumptions still require time, mesh and model assessment. No gas momentum,
evaporation, arbitrary curved-surface law or total mechanical-energy closure is
introduced. Focused regressions in
`tests/material_free_surface_enthalpy/conduction/surface/convection.rs` exercise
the actual surface root and liquid resistance, finite-gas reciprocal energy and
step admission, exact equilibrium, late refusal with identical accepted retry,
and the original implicit endpoint residual.

## Curved D3Q19 walls with native open reservoirs

`step_with_curved_walls` also admits the original velocity/pressure open faces.
It streams the actual curved populations first, then invokes the SAME cached,
regularized open-face owner used by immersed-body coupling. The first-interior
stress and prescribed velocity/density use their native Guo force frame,
including spatially varying held forces and existing inlet profiles. At mixed
wall/open nodes, the native mask retains every curved or halfway wall-owned
population; only the open directions reconstruct.

`CurvedWallStep3::open` reports reservoir mass, linear impulse and angular
impulse TO the fluid. The shared compensated accumulator counts the temporary
open-link reflection and every actual reconstruction delta, including the
replacement of tangential populations on face interiors. These are direct
operator contributions, not a residual inferred from a whole-domain balance.
They remain distinct from curved-body reaction, other walls, held volume-force
input and periodic coordinate transport. BFL population mass drift remains
measured without renormalization. Reservoir exchange is not attributed to the
body, and is not labeled pressure work or a complete fluid-energy balance.

The candidate reservoir populations and accounts finish before the existing
single publication swap. A late reconstruction, arithmetic or cancellation
refusal preserves accepted populations and geometry. The open-face cache,
post/next population buffers and compensated boundary owner are reused; no
second collision law, reservoir law or whole-volume balance scan is added.

Four focused G0/G1/G3 native regressions pass: exact halfway recovery for both
velocity/pressure and pressure/pressure ends on all three axes and all four
collision families; rotating curved-body mass/linear/angular balance with
separate reservoir, periodic, stationary-wall, uniform and spatially varying
force contributions; wall-corner population preservation and steady
pressure-driven flow between fractional walls; and late reservoir-account
overflow without publication. The resolved 96x16x4 duct with a 14.6-cell gap
has relative axial-profile error 0.009887809049680325 under a 3% budget. A
6.6-cell gap gave 0.046973215227596746 and is outside that envelope; no global
order claim follows. The six original curved-boundary regressions also pass.
These are direct original numerical-owner checks with unused runtime/pool,
body-motion and unrelated API wrappers omitted. Full Cargo, DSR, cross-ISA
and performance gates were not run. The fixed-occupancy and prescribed-wall-
kinematics limits of the curved owner remain unchanged.


### Finite translating plane carriers on the physical D2Q9 disk

FiniteLubricationPlane2 and
DryDisk2::step_with_translating_lubrication_planes advance finite normal
AND tangential plane momentum within the original dry-disk/fluid transaction.
The supplied mass is per unit out-of-plane depth; center and surface point
are distinct. Orientation remains fixed. The guide receives only a pure
reaction couple, with zero work and no supported linear impulse.

The original UnderResolvedCylinderPair2 owner supplies the normal
resistance, cutoff and gap floor at the actual first-half-drift geometry.
The shared finite-mass elimination condenses normal and tangent resistance
independently with k/(1+k/(2m)), using one carrier mass. Both enter every
existing three-unknown dry-fluid candidate, with no later impulse kick or
extra fluid collision. Original cylinder spin resistance is unchanged.

Each accepted carrier moves its center AND plane point through both free
half-drifts. The second uses its solved endpoint velocity. Admission checks
the actual carrier/disk endpoint gap, never a prescribed-wall extrapolation.
With fixed normal, the relative gap is affine within each drift; positive
start/midpoint/endpoint gaps exclude crossing along these two trajectories.
Complete vector low-Mach admission applies before, at midpoint and after.

Actual carrier momentum and kinetic increments, authored loads/work and
the orientation-guide couple are independently checked before native
publication. The aggregate exterior-fluid/dry-disk/carrier momentum account
includes guide couples as external sources. Refusal/cancellation preserves
all accepted states. Only infallible assignments publish the already
admitted finite endpoints.

Focused native tests reconstruct three successive two-plane steps from
actual changing geometry, coefficient laws, momentum and work. They check
recoil that rescues a falsely predicted held-wall crossing, true first/last
drift crossings, complete vector speed refusal, exact retry and inactive
parity. This supplies missing-film physics, not resolved moving walls,
finite face extents, contact or a thermal destination. The inherited
implicit midpoint method is not L-stable; physical stiff-film relaxation
still needs time resolution.

## Coupled implicit internal energy of one finite atmosphere

Selecting `FreeSurfaceInterfaceSolve::implicit` on the existing finite-gas
surface or surface-convection entry couples one gas energy coordinate to the
original liquid/solid backward-Euler residual and bounded Newton-Krylov driver.
Every trial reconstructs a private parcel through its original isochoric caloric
inversion. Its actual temperature, EOS pressure, composition and selected source
claims feed the same interface root as the trial liquid endpoint. No physical
history advances in a residual or directional-difference evaluation.

The added gas row is original-law `m_g Delta u_g + dt sum(Q_from_gas) = 0`.
It uses actual represented storage, with formation references canceled by the
original temperature-anchored caloric law. Liquid and solid rows retain the
original `BackwardEulerStep` heat-rate and local/total energy gates. The gas row
includes both solved liquid surfaces and every original direct solid/gas film.
A separate union gate compares actual liquid-plus-solid-plus-gas storage with
independently accumulated authored/molecular and liquid/solid radiant heat. Gas exchange is internal to that union.
The accepted observer validates the already prepared gas endpoint and accounts
for its surface heat once; it never applies that heat to the parcel a second
time. Native flow, source, pressure and final endpoint failures still discard
all candidates before publication.

`maximum_unknowns` includes the gas coordinate. Existing residual, Newton,
Krylov, interval, source-inversion and interface-root budgets remain operative.
`cell_energy_tolerance_j` bounds each liquid row, each solid row and the gas row
per accepted interval. `GasEnergySolve` retains its separate original caloric-inversion
controls. `FreeSurfaceImplicitReport::maximum_gas_residual_j` exposes the actual
accepted gas defect; `FreeSurfaceFiniteGasExchange::implicit_energy_tolerance_j`
identifies the accumulated declared gas-row allowance within the reported
combined gas budget. Conditioning scales do not replace any raw joule gate.

This removes the explicit gas-capacity CFL restriction; it does not establish
timestep accuracy or unconditional nonlinear convergence. The original
first-order thermal/native-flow split and fixed gas mass, composition and volume
remain. Current and separate liquid environment/radiation paths remain excluded
from this implicit rung. Focused public-API regressions compare the
stiff coupled endpoint with an independent caloric/series-resistance equation,
check large formation-reference invariance and evolving gas/film source states,
close the radiating liquid-plus-gas union, and verify late-domain, cancellation
and work-budget refusal with an identical accepted retry.


### Shared finite translating film carriers in D2Q9

TranslatingFilmCarrier2, CarrierLubricationPatch2 and
DryDisk2::step_with_lubrication_carriers attach several original cylinder/
plane films to one positive finite mass. Each stable carrier index owns
one center, velocity, authored WORLD impulse and kinetic store, independent
of its patch count. Patches retain immutable COM-relative points and fixed
world normals; the orientation guide supplies only a zero-work couple.

For a carrier, R sums every attached normal/tangent resistance tensor.
Its midpoint solves (2m I+R) u_mid = 2m u_before + J_authored + R u_disk,mid.
The shared response retains off-diagonal coupling between oblique faces;
independent mass elimination per patch would produce different forces.
Original native normal resistance supplies its cutoff/floor. Cylinder
tangent/spin coefficients and viscosity semantics remain unchanged.

One scaled native 2x2 Cholesky response is prepared per carrier. The
existing maximum_condition_number also admits this response. Every
retained-fluid Newton candidate evaluates all carrier midpoints once in
bounded stack storage; no carrier vector or factorization is allocated
there. Each patch uses the common motion for its force, couple and work.
Aggregate carrier momentum/kinetic storage is counted exactly once.

Midpoint and endpoint patches are reconstructed from the actual common
carrier COM. Both relative-gap half-drifts are admitted. COM/local
cancellation that prevents a patch from following its shared translation
at the local lattice geometry scale refuses before publication. All
carrier trajectories, aggregate guide torques and independent momentum/
work checks finish before the original sole native disk/fluid commit,
followed only by infallible carrier assignments.

Focused tests independently assemble/invert the unscaled coupled matrix
for unequal oblique faces across three moving steps, demonstrate failure
of independent per-face masses, and reconstruct actual system momentum,
guide torque and unique kinetic work. Single-face equivalence, inactive/
empty bit parity, one authored impulse per carrier, conditioning/index/
unused-carrier/speed/gap/representability refusals and exact retry pass.
At most 64 carriers and 64 patches are admitted, with every carrier used.
These remain local infinite missing-film patches: no resolved moving wall,
finite extent, contact, GCL, rotation or thermal destination is inferred.
The inherited midpoint method is not L-stable; stiff relaxation still
requires adequate time resolution.

## Sparse liquid curl and full-gradient WALE/Vreman

`FreeSurface3::kinematics` prepares an explicitly byte-budgeted
`FreeSurfaceKinematicsWorkspace3` over active tiles, using actual retained
Guo forces and native population density. Its borrowed view returns the shared
`FlowKinematics3` (full gradient, curl, divergence, enstrophy, Q and stretching)
at occupied liquid centers. Queries allocate nothing. The cache is rebuilt
after each new preparation, including sparse slot reshuffles and piston motion;
the view borrows the owner so geometry cannot change underneath a query.

Centered and contiguous one-sided three-point differences are second order.
When only one inward liquid neighbor exists, diagnostics expose their first
order explicitly. No adjacent liquid on either side of an axis is a local
`MissingDerivative` refusal. Gas, zero-fill interfaces, solid/inactive tiles,
out-of-domain space and the current piston slab never supply fake stationary
fluid. These are reconstructions from actual interior liquid velocities, not
imposed wall-law gradients: transient prescribed wall velocities are not
retained by this owner and are not guessed or extrapolated from wall nodes.
The same rule avoids assuming an unmodeled free-surface tangential stress.

`prepare_gradient_sgs` and `step_with_gradient_sgs` separately opt into the
shared WALE/Vreman formulas using this full, force-centered gradient. Every
occupied liquid cell needs second-order support; unresolved thin sheets refuse
the whole step. Zero-fill interface ghosts retain the original molecular
relaxation: they carry populations for gas reconstruction and mass exchange,
but no occupied liquid supports an inferred gradient or modeled dissipation.
Diagnostics still return None there. Reports count these native wet collision
cells while assigning their zero liquid volume zero molecular and SGS power.
The original stress-based Smagorinsky arithmetic is unchanged. BGK/TRT retain
their molecular family and authored TRT product, with effective shear
`tau = tau_molecular + 3*nu_t`. Eddy/molecular limits refuse rather than clamp.
Reports weight deviatoric shear dissipation by native fill and density; modeled
viscosity sees the full gradient, including antisymmetric rotation. No reported
dissipation is added as heat. The original VOF owner performs collision, gas
reconstruction, mass transport and conversion on one candidate before commit.
Zero coefficient uses the original native step bit-for-bit on admitted support.

Preparation and sparse scans checkpoint cancellation; rejected queries and
steps never mutate accepted flow. The kinematics and native SGS tests cover analytic
affine curl/stretching with unequal density and forcing, quadratic interface
gradients in all orientations, solid/piston masks, under-resolved refusal,
budget/cancel retry, full-gradient viscosity and fill-weighted power, and native
BGK/TRT evolution with zero-coefficient controls, including a finite liquid
block and its empty interface shell. These are exploration closures
and resolved diagnostics, not turbulent-flow validation, entrainment physics,
wall modeling, baroclinic generation or a complete vorticity evolution equation.

## Finite solids in the implicit moving-liquid surface solve

The opt-in surface residual retains each `FreeSurfaceThermalBody` as the same
original `DepositionSimulation` caloric cell used by explicit conduction.
Liquid/solid contacts and solid/solid links enter its original pair operator;
there is no duplicate wall store or second heat integrator. Each row uses its
own original caloric slope for numerical conditioning. Source conductivity,
ordered interface resistance and both contact states are re-evaluated at the
actual trial endpoint, without freezing their physical coefficients.

A finite atmosphere's one energy equation includes the opposite of every
liquid-surface and direct solid-film transfer. Original solid radiation remains
external to the liquid/solid/gas union. Raw per-cell, gas-row, aggregate thermal
and independent union joule gates remain in force. Only the accepted endpoint
advances the original observer, and the native flow acceptance still precedes
publication of liquid, all bodies and gas together.

This extension retains stationary solid geometry, the declared scalar bulk and
contact/film models, and each body's uniform-cell interpretation. A source or
fully-solid phase-domain exit refuses the complete transaction; it does not
invent a melting or mechanical handoff. Focused regressions check independent
five-store backward-Euler equations with distinct material capacities, evolving
contact/film sources, radiation and one shared gas; a separately solved
prescribed-reservoir contact; and late source/cancellation rollback followed by
an identical accepted retry.

## Convection on original reconstructed liquid interfaces

Thermal substeps and implicit trials cache normals once from the original held
height-graph or PLIC reconstruction, after existing complete-owner and area
admission. Native endpoint queries use the newly reconstructed owners, areas
and normals. The accepted `FreeSurfaceInterfaceSample` retains the exact
`convection_outward_normal` used at that epoch. Reusing a prepared thermal
sample against changed native orientation refuses, alongside existing stale
owner and area checks. Native endpoint evaluation updates surface-temperature
and capillary samples; it supplies no second heat debit or gas advance.

Local `UniformReservoir` evolution additionally requires the provider's
explicit `admits_uniform_boundary_evolution` capability. The default refuses.
The capability declares spatially uniform authored flow, dimensions and source
boundary data; a coefficient may still depend on actual temperatures and local
geometry. Newly exposed owners use original boundary index zero and their real
face. Existing uniform-source/context admission remains authoritative, and
fixed-owner evolution retains its original owner-change refusal.

Manufactured positive normal-dependent coefficients test actual curved
height-graph and PLIC paths without asserting a flat-plate model for curvature.
The cases independently check the liquid-half-cell/surface-convection root,
mean-film source queries, held and live geometry, finite-gas reciprocal energy,
changing PLIC ownership and a late constitutive refusal followed by accepted
retry. Their coefficient evidence carries an estimate and no calibrated model.

## Moving-wall WALE/Vreman with material SI receipts

`planar_wall_velocity_field_gradient_sgs_candidate` and its publishing entry
prepare the full-gradient closure from accepted force-centered liquid velocities
and the SAME prescribed halfway wall-face velocity used by native reflection.
Quadratic wall/center/interior reconstruction supplies second-order support at
the wall-adjacent cell and both centers of a two-layer liquid film. Gas, empty
interfaces, inactive solids and piston slabs remain absent velocity samples;
a single occupied wall/gas center still refuses. The public liquid-only
kinematics observer keeps its existing stencil semantics.

The resulting WALE/Vreman viscosity directly supplies the existing local
collision relaxation map. That same collision produces reflected populations
and actual wall impulse, torque, residual couple and actuator work. Empty
interface carriers retain molecular relaxation, and the viscosity-ratio gates
remain active. Zero coefficient passes no local tau map, preserving original
wall trajectories and receipts on admitted support. The native collision family,
independent molecular bulk rate and source-rheology integration are retained.

The corresponding `MaterialFreeSurface3` SI entries retain the original
material, force frame, capillarity, pressure datum and source/scenario limits.
Filter width uses `Delta_lu = Delta_m/dx`; viscosity uses `dx²/dt`, volume
`dx³`, and integrated deviatoric dissipation `rho0*dx⁵/dt³`. All native checks
and SI conversions finish before publication of the same candidate. Positive
quantities lost to conversion underflow and nonfinite receipts refuse.
One native clone is used; no thermal route history is allocated.

The scope remains a fixed, fully flooded tangential wall, with gradients
reconstructed from liquid and the prescribed halfway face. No empirical wall
function, heat deposition or turbulence validation is inferred. Existing
material regression groups cover disabled-coefficient BGK/TRT parity, enabled
physical evolution and reciprocal wrench/work with independent SI conversions,
stencil/budget/cancellation refusal, and late SI-conversion rollback.
`sgs::gradient::wall_tests` adds all six wall orientations, quadratic shear and
rigid rotation, force/Galilean invariance, two-layer films, and actual BGK/TRT
and independent-bulk collision populations and wall impulses compared with an
independently prescribed analytic-gradient relaxation map (G0/G3).

## Source-resolved Newtonian films

`film::material::MaterialFilmLiquid` borrows the original
`ResolvedLiquidSurfaceState` used by material VOF, retaining its exact query,
viscosity basis, receipts and capillarity decision. It supplies native
Nusselt/Benney coefficients without source re-selection, cloning or caloric
prerequisites. Evolving shear/pressure rheology refuses; only explicitly
neglected capillarity maps to zero. Repeated film queries allocate nothing;
temperature/composition changes require a newly admitted source state.
`tests/material_film.rs` covers both viscosity bases, pinned temperature and
receipt preservation, dimensional growth, explicit capillarity omission and
rheology refusal using synthetic source inputs, without a calibration claim.

## TRT on coupled laminar plate flow and heat

`PlatesGrid3::{set_trt_magic,set_bgk_collision,collision_model_at}` selects
the original D3Q19 BGK/TRT kernel while retaining each cell's current
material/rheology tau and held Guo force. `try_step` admits original storage
and completes every collision in reusable scratch before streaming; a late
cell refusal leaves accepted populations, forces, tau and selection unchanged.
No new collision law, per-step allocation or extra collision pass is introduced.
BGK remains the default.

Original D2Q9/D2Q5 and D3Q19/D3Q7 source transactions admit TRT, preserve its
product through snapshots, and validate every implied odd relaxation before
accepting a local viscosity field. Heat/force units, thermal balance, source
history and rollback stay in the common thermal/material owners. Other
D2Q9 collision operators retain their previous refusal.

Regressions in `d3q19_coupled`, `thermal_sources` and `local_thermal_transport`
cover planar flow with a BGK slip control, constitutive retuning, refusal/retry,
native/source thermal history, physical source conversion and spatial relaxation.
The 3/16 wall-location claim is restricted to uniformly forced, fully developed,
constant-viscosity planar flow. It does not certify variable-viscosity wall
accuracy, interface flow, nonlinear convection or global stability.

## Continuous long-wave film bands

`NusseltFilm::long_wave_band` maximizes the native Benney growth over an
explicit closed nonnegative wavenumber interval. The quartic dispersion needs
only its endpoints and, when present, its interior stationary point: at most
three native mode evaluations, without allocation or an arbitrary sampling
resolution. The entire interval must pass both long-wave budgets. The result
retains the maximizing wavenumber and actual band-wide kh/Re*kh maxima.
Including zero retains the neutral mass-conservation mode; zero capillarity,
point intervals and boundary maxima retain their physical meanings.

This is an Estimated maximum of the declared local model, not a floating-point
enclosure, arbitrary-wavelength stability proof or full free-surface eigensolve.
`tests/film_band.rs` checks independent dimensional extrema, an unstable
interior missed by endpoint samples, monotone/neutral limits and tail refusals.

## Evolving Mann inflow across native characteristic intervals

`NsfEvolvingMannInflow::preview_boundary_interval` evolves the existing inlet
at every original fixed-dt NSF step boundary. The profile at time n is held
through native step n to n+1; a final endpoint refresh aligns both accepted
clocks without adding a fluid transfer. Point samples, face/cell averages,
mean/frame/gas and planar filtering reuse the existing sampler.
`NsfEvolvingMannIntervalBudget` bounds total cell steps, OU transitions
(including initial catch-up and endpoint), and N+1 tile-mode profile evaluations.

The existing `NsfCharacteristicTrial` owns the single fluid candidate and
unchanged chronological mass, energy, pressure-work and covariance aggregation.
The inlet reuses its source trial and profile buffers. `NsfEvolvingMannTrial::fluid()`
exposes those interval receipts before acceptance; the endpoint rebind retains
the original control-update semantics for the candidate's last-step queries.
Acceptance publishes both owners. Drop, refusal or cancellation retains both
accepted states for exact retry. Standalone `apply` semantics remain unchanged.

This is native stepwise prescribed synthetic inflow, without interpolation
within predictor passes or a higher-order, nonreflecting or eddy-turnover
qualification. Regressions in `mann_inflow/evolving/interval_tests.rs` compare
all six point/face/cell and homogeneous/plane observations with explicit
apply/step/endpoint sequences, demonstrate changed resolved velocity and energy
exchange relative to held inflow, and cover budgets, catch-up, cancelled/dropped
retry and native fluid refusal.

## Live source rheology under sparse SGS

`d3q19::sparse::sgs::material::rheology::SparseRheologySgsWorkspace3` and
`MaterialScaling::{sparse_rheology_sgs_workspace,
step_sparse_with_material_rheology_and_sgs,
step_sparse_pooled_with_material_rheology_and_sgs}` resample original
`FluidTransport` responses at incoming resolved shear on the native sparse
trajectory. Existing source points/samples, robust `sqrt(2 dev(D):dev(D))/dt`,
momentum scaling and the force-centered wall-aware sparse gradient are shared.
One reconstructed gradient supplies both rheology and SGS. Each response is
read once per active cell; Morton coordinates decode once per tile.

The local molecular response supplies collision viscosity, molecular
dissipation and the SGS ratio. Zero SGS retains local rheology. Topology,
dx/dt/reference density and the native collision law bind the workspace.
Two reserved source buffers introduce no adapter step allocation; storage
owned by returned material responses and native pooled scheduling remain
with their original owners. Accepted populations, step count, source epoch
and report survive source, ratio, native-step or cancellation refusal.
Callback side effects remain caller-owned and cannot be rolled back.

This is explicit incoming-time resolved-shear generalized-Newtonian LES.
Evaluating nonlinear molecular stress at filtered strain does not close its
unresolved correlations. This adds no constitutive memory, implicit coupling,
thermal/composition evolution or turbulence-validation authority. Four groups
in `material/rheology/tests.rs` cover frozen-source and zero-SGS equivalence,
physical shear and live viscosity feedback, source/ratio refusals, serial/pooled
parity, late cancellation and preserved accepted source epochs.

## Prescribed moving surfaces on shared D3Q19/D3Q27 transport

`Grid3::{step_serial_with_moving_wall_momentum, step_pooled_with_moving_wall_momentum}` consume a borrowed wall-velocity field in `WallMomentumProbe3::wall_links()` order. The original pull stencil applies the standard low-Mach moving halfway-bounce correction using post-collision density before admitting the final streamed state. Optional material/SGS collision fields retain their existing owner and step-local meaning. The fixed topology and sparse link selection are reused; no additional population buffer or per-step wall-field allocation is introduced. An all-zero velocity field selects the original stationary stream path and preserves stationary population and wall-load arithmetic.

Exterior surface velocity must be tangential to every physical face crossed by a selected link. Each fluid cell's selected moving-wall corrections must sum to zero mass source to roundoff, including interior solid links. This tests the weighted sum of link-direction dot products, not each dot product separately, so tangential shear remains possible. Incomplete stencils and normal/swept motion refuse rather than leak mass or project the prescribed velocity. Geometry is fixed: no fresh-cell initialization, remap or moving-body geometry is implied.

The receipt reuses `WallMomentumExchange3` and reads actual outgoing and corrected incoming populations. Boundary-relative impulse follows Wen et al., https://arxiv.org/abs/1303.0625. Separate fluid-population, wall-velocity mass, and angular terms retain `wall + fluid_population = wall_velocity_mass`; measured population mass change is also reported. Torque acts at the halfway intersection. Wall work is the sum of actual link impulse dotted with prescribed wall velocity; its negative is actuator work supplied to the fluid. This is boundary mechanical work, not a thermal balance or a claim that fluid kinetic-energy change equals input work. Canonical compensated reductions are independent of worker count. Collision, motion, final-state, receipt and cancellation checks all finish before the original population/step-count publication.

`tests/lattice3_moving_walls.rs` covers analytical Couette velocity, viscous drag and wall power, initial half-space momentum, zero-slip entrainment, changing drives with per-cell viscosity, native momentum/mass and relative angular identities, zero-wall bit parity, worker counts 1/2/7, refusal and post-stream cancellation/retry.

## Native CFD intervals driven by OU timestep means

`NsfEvolvingMannInflow::preview_time_averaged_boundary_interval` is an opt-in
consumer of `fs_atmo::mann_field::fixed_step::MannStepAverage`. Each original
native CFD step receives the exact linear velocity mean over its corresponding
advected OU interval, including the inlet's existing point/face/cell observation,
proper coordinate rotation and optional harmonic slip-plane correction. The
original source gas, actual ghost locations and reservoir profile rebind remain
authoritative. No second flow solve or extra convection displacement is added.

An external prepared mean workspace is reused with the inlet's temporal trial and
profile buffers and the original single characteristic fluid candidate. Whole-
interval bounds cover native cell steps, normal draws (four per mode for prior
catch-up, eight per mode per averaged step) and all `steps + 1` spatial profiles.
Catch-up runs only after native interval eligibility and cell work are admitted.
The final instantaneous profile is restored at the common fluid/source endpoint;
ordinary same-time apply is idempotent and the endpoint checkpoint path is
unchanged. Trial acceptance publishes both fluid and atmosphere; refusal,
cancellation or dropping the trial publishes neither. The external mean can
retain a rejected candidate observation and is explicitly scratch, not a third
accepted owner.

The original interval retains its chronological mass, energy, pressure-work and
exchange-variation accounts. These are the native fluxes produced by a velocity
held at its exact linear step mean, not exact averages of nonlinear momentum or
energy fluxes. Spatial/advective Nyquist admissions remain active; OU noise is not
band-limited. No added SGS kinetic energy, acoustic source or higher-order CFD
accuracy is asserted.

Three regression groups cover independently constructed ghost controls through
all six wall/spatial policies and native transfer accounts, bound-clock catch-up,
whole-interval budget admission before OU work, cancellation/drop/retry, final
profile idempotence and changed native flow with unchanged source endpoints.
The full native fs-lbm library and unit tests with `turbulent-inflow` and the
boundary-SGS integration targets typecheck; these new CFD tests were not executed
in the resource-limited session. The underlying atmosphere's 33 Mann-field tests
were executed successfully.

## Finite-inertia planar drag on shared D3Q19/D3Q27 transport

`lattice3::finite_wall::{FinitePlanarWall3,PlanarWallState3}` binds a fixed
planar coupling window to two tangential material velocities and spin about
its positive Cartesian normal axis. Positive mass and axial inertia are
physical inputs. Complete face-cell stencils, a center on the wall plane,
nonperiodic normal and unchanged voxel topology are required. Corners crossing
another physical surface refuse. The window and inertia are frozen: a moving
finite footprint, changing geometry, pose integration and conservative
moving-volume remap remain separate work.

`Grid3::{step_serial_with_finite_wall,step_pooled_with_finite_wall}` accepts
external tangential impulses and an axial couple, permitting driven motion
and coasting under actual fluid traction. One original collision pass supplies
`b` and the positive-semidefinite resistance `K` in `J(v)=b-Kv`. The existing
prepared LU workspace solves `(M+K/2) delta_v=b-Kv_before+J_external`, after
normalizing by the square root of the physical inertia. Streaming uses the
solved midpoint field through the original moving-wall owner. Complete planar
stencils make its moving-mass terms cancel from generalized translation and
axial spin; actual Wen-relative impulse, torque and wall work are nevertheless
measured and checked before acceptance. No drag coefficient, force
extrapolation, repeated collision or additional population buffer is used.

The receipt separates finite-wall traction, other stationary-wall loads,
zero-work normal/transverse support reactions, external drive and held Guo
fluid forcing. It reports the polarized wall kinetic-energy increment and
external work. Declared relative impulse/work budgets gate the actual measured
exchange. Cancellation or any admission failure publishes neither fluid nor
mechanical state; accepted updates publish together. Sparse selections,
velocity storage and the three-by-three factor workspace are prepared once.
Optional per-cell collision fields retain their original material/SGS meaning.
Isothermal fluid kinetic-energy loss is not wall heat, and this local receipt
does not supply missing periodic global angular transport.

`tests/lattice3_finite_wall.rs` covers independent half-space resistance and
midpoint inertia on all six faces for both lattices, BGK/TRT variants,
off-center translation-spin response, driven/coasting evolution with varying
viscosity and prescribed-motion population parity, linear momentum and work,
worker parity, and post-stream cancellation/refusal with exact retry.

## Dense boundary material rheology with native SGS and spatial wall loads

`material::boundary::sgs` adds a prepared `BoundaryRheologySgsWorkspace3` and
serial/pooled `MaterialScaling` steps on the original `BoundaryGrid3`. One incoming
force-centered velocity gradient supplies both the physical resolved material
shear `sqrt(2 dev(D):dev(D))/dt` and the selected existing SGS closure. Isotropic
dilation does not change the molecular shear query; the original SGS gradient
and eddy-viscosity law retain their own constitutive convention. Each original material
response is queried once and its transport state read once. The original scaling
converts that local molecular viscosity to collision relaxation; its reference
density must remain unchanged. Eddy/molecular ratio admission and resolved
viscous dissipation use the local molecular value, and zero SGS still applies
heterogeneous material transport. Collision family, TRT magic product and
non-shear controls retain their original owners.

The workspace retains two preallocated source-response epochs, molecular tau and
the last accepted native SGS/wall-load report. Sampling, unit-binding, source,
closure, cancellation and native-flow failures preserve accepted samples and
populations. Source callbacks must not publish external history. Material
response heap storage remains caller-owned. Native pooled chunk allocations and
the existing initial stream-buffer allocation retain their own limits.

`WallPatchProbe3` partitions a unique canonical wall-link selection into authored,
nonempty spatial patches with finite moment origins and an explicit payload byte
budget. It retains two exchange buffers. Ordinary and SGS patch steps measure
the actual collision that streams that step, through the same link traversal and
reflection arithmetic as aggregate wall loads. Every patch receives its force,
moment, mass exchange and moving-wall work; fluid publication precedes the
infallible exchange-buffer swap. Native topology and periodic moment-coordinate
qualifications remain unchanged. Measurement needs no additional collision,
full-domain force scan or per-step patch allocation. The shared serial stream
helper preserves the existing transactional publication path.

This material law is evaluated on resolved incoming shear; it does not close the
nonlinear correlation of unresolved strain with molecular stress. No thermal
integration, constitutive-memory evolution, immersed-fluid filtering, new wall
model, stochastic SGS force or turbulent calibration is supplied. Retained
fictitious interior fluid makes immersed-body LES a separate unresolved physical
question; this capability is restricted to the original dense wall/open domain.

Material integration regressions cover all four collision families and
wall/open boundaries, independent heterogeneous collision/transport with zero
and nonzero SGS, dilation-invariant molecular shear with active SGS,
analytic forced Couette shear and local dissipation, retained
source epochs on refusal, and serial/pooled replay. Three spatial-load groups
cover all-axis momentum and aggregate agreement, zero-SGS parity, memory and
partition refusal, and late cancellation. Full native fs-lbm production/unit
sources with `turbulent-inflow`, plus `material_boundary_sgs` and `boundary_sgs`
integration targets, typecheck after merging concurrent boundary changes.
These runtime tests were not executed in this session.

## Prescribed thermal plates with source-owned viscous self-heating

`PlatesGrid3::set_plate_motion` retains uniform tangential velocities on both
fixed halfway z-plates and an explicit torque origin. The original BGK/TRT
collision and pull stream use `f_in = f_out + 2 w rho (e.u_wall)/cs²` with
post-collision density. Complete tangential planes introduce no modeled mass
source. A boundary-only preflight uses that same reflection, admits finite
positive-density/Mach-bounded endpoint populations, and accumulates compensated
per-plane momentum, angular impulse and work before publication. It adds no
population buffer or wall map. Zero drives preserve stationary population bits;
invalid drives and rejected native steps retain accepted state and receipts.

The original thermal/source candidate snapshot carries this drive and its last
accepted receipt. `MaterialSlab3::set_plate_motion_si` uses the retained dx/dt;
`last_plate_exchange_si` converts the actual native exchange with the retained
rho0, cell volume and time unit. Negative wall work is energy supplied by the
prescribed actuator. The existing `step_with_material_rheology_and_viscous_heating`
resolves the same actual-temperature, force-corrected shear material for momentum
and `mu*gamma_dot²` heat. Wall work is not added again as thermal power. Both
source epochs, thermal wall exchange, populations and mechanical plate receipt
publish together through the existing checked material transaction, including
late native failure, endpoint-source refusal and cancellation. Low-level thermal
steps retain their original nontransactional status; native failures now return
through `ThermalForcingError::Flow3` to the checked owner instead of panicking.

Focused regression definitions cover the analytic Couette velocity and viscous
heating temperature profile under refinement, independent physical wall work and
heat conversion, retained thermal Carreau-Yasuda feedback, stationary bit parity,
and native/late-cancellation refusal followed by exact retry. They use declared
synthetic material calibrations, not experimental validation. This fixed-geometry,
low-Mach, generalized-Newtonian path adds no swept-cell transport, finite actuator
reservoir, elastic memory, numerical-loss heating or transient total-energy proof.

## Independent bulk viscosity in the dense D2Q9 central-moment solver

`CollisionModel2::CentralMomentBulk { bulk_tau }` separates the second-moment
trace from the deviatoric modes. The shear modes retain each cell's `Grid::tau`;
the trace uses the finite material relaxation time `bulk_tau > 1/2`. Third and
fourth central moments keep their original equilibrium treatment. The original
`CentralMoment` option remains exactly the `bulk_tau = 1` arithmetic, including
forcing and streaming. The bulk reciprocal is prepared once per collision sweep.

In dimension D, `nu = cs² (tau_s - 1/2)` and
`zeta/rho = (2/D) cs² (tau_b - 1/2)` in lattice units. Therefore D2Q9 uses
`zeta/rho = cs² (bulk_tau - 1/2)`, not the D3 factor `2/3`. For planar flow in a
three-dimensional material, the corresponding two-dimensional coefficient is
`zeta_2/rho = zeta_3/rho + nu/3`. Scale the appropriate dimensional kinematic
coefficient by `dt/dx²` before inverting the D2 relaxation law. Linear isothermal
longitudinal waves decay at `Gamma = (nu + zeta_2/rho) k²/2`; transverse waves
retain `nu k²` attenuation. The planar 3-D conversion therefore recovers
`Gamma = (4nu/3 + zeta_3/rho) k²/2` without confusing dimensional conventions.

The same physical-population Guo transform applies each force moment at its own
relaxation rate. Uniform and spatial rheology remove isotropic stress from the
shear-tau-dependent deviatoric part, reconstruct dilation using the configured
bulk tau, and retain the existing full strain invariant `sqrt(2 D:D)` for the
constitutive solve. Changing bulk tau must not be misinterpreted as material
shear thinning. The checked endpoint-force step rejects invalid bulk relaxation
before numerical proposal work; the low-level collide API retains its assert
contract for invalid relaxation parameters.

Focused regressions cover independently manufactured trace/shear relaxation
under Guo forcing, compatibility at unity, two-wavelength continuum sound decay,
shear-wave and force-driven channel viscosity, source ramps/reversals, and
uniform/spatial mixed-strain rheology. This is a weakly compressible isothermal
bulk-viscosity capability; it supplies no thermal acoustic attenuation or
frequency-dependent molecular relaxation. D3Q19 and D3Q27 retain their existing
collision contracts.

## Independent bulk viscosity in the D3Q19 central-moment solver

`CollisionModel3::CentralMomentBulk { shear_rate, bulk_rate, higher_order_rate }`
relaxes the five deviatoric second moments independently of their trace. All
three rates must be finite and strictly between zero and two. In lattice units,
`nu = (1/shear_rate - 1/2)/3` and
`zeta/rho = (2/9) (1/bulk_rate - 1/2)`; the latter is also exposed by
`kinematic_bulk_viscosity`. The existing models retain their original bulk/shear
relation. Equal bulk and shear rates preserve the original central-moment
population arithmetic, including the selected higher-order rate and Guo force
transform. The extra relaxation corrects only the isotropic second moment.

Native material and lattice SGS shear rebinding preserve the independently
configured bulk and higher-order rates. Curved-wall source rheology extracts
the deviatoric stress numerator from the physical Guo stress for every collision
family, then queries the original constitutive source at `sqrt(2 D_dev:D_dev)`.
The bounded local solve retains source-domain, iteration and cancellation
refusals. Isotropic stress changes neither this shear query nor its apparent
viscosity, while independently relaxed dilation remains available for physical
bulk heating. Native boundary and immersed-solid relaxation
continue to use the shear rate for their viscous response. The separate
free-surface stress SGS path also admits independent central moment while
preserving the original bulk and higher-order rates.

This is a weakly compressible, isothermal collision capability. Longitudinal
continuum waves decay at `Gamma = (4nu/3 + zeta/rho) k²/2`, while transverse waves
decay at `nu k²`. The collision operator alone adds no thermal equation,
molecular relaxation or temperature-dependent bulk-viscosity sampling. The
separate thermal-plate coupling below explicitly supplies physical bulk heat.

`MaterialForceGrid3::try_new_with_bulk_viscosity` binds the existing
`BulkViscousMaterial` snapshots to transient 3-D material forcing. Supported
independent-trace collision families refresh shear and bulk from one accepted
material epoch while preserving their other relaxation rates. Physical bulk
uses the pinned Boussinesq reference density and original dx/dt mapping; local
density departures remain buoyancy inputs. An accepted endpoint supplies the
next collision's viscosities, preserving the existing explicit lag and atomic
source/population/force publication. Families without an independent trace mode
refuse this constructor. This momentum owner adds no thermal equation or heat
reservoir; the separate thermal owners retain those responsibilities.

## Physical shear and bulk heating on 3-D thermal plates

`PlatesGrid3::set_central_moment_bulk_collision` retains independent bulk and
higher-order rates while local material/rheology fields change only shear.
Checked thermal snapshots preserve that selector, the held Guo force frame,
plate state and directional buoyancy. Selecting numerical damping alone does
not create a heat source.

`MaterialSlab3::step_with_material_rheology_and_bulk_viscous_heating` takes one
explicit physical `zeta` in Pa s, spatially uniform and held for that step. The
same coefficient sets `tau_bulk=1/2+(9/2)(zeta/rho_ref)dt/dx²` and volumetric
heat `zeta(div u)²`. The retained source callback supplies local
`mu(T,sqrt(2 D_dev:D_dev))` and conductivity. Its bounded constitutive solve separates
the Guo-corrected deviatoric numerator from independently relaxed dilation.
Cached components at the actual accepted query give
`q=2mu D_dev:D_dev+zeta(div u)²`; isotropic shear power is not added a second time.

The existing source half-kicks, material endpoint checks and independently
measured wall/sensible-heat balance publish the collision selection and thermal
state together. Cancellation or any source/conversion/temperature refusal keeps
the accepted slab. Higher-order relaxation and numerical kinetic-energy
residuals supply no heat. The earlier heated API refuses independent bulk unless
this explicit physical binding is provided. The ordinary BGK/TRT heated path
uses the same deviatoric material query while retaining its collision-implied
`zeta=2mu/3` and the complete `2mu D:D` work. Pure dilation therefore supplies
bulk heat even when the material remains at its zero-shear viscosity.

This remains explicit low-Mach, fixed-density/Cp sensible-heat transport. It
does not assert compressible total-energy conservation or an implicit heating
solve. This entry point holds bulk viscosity uniform for a step; the source-local
extension below supplies temperature/strain-dependent bulk viscosity. Positive bulk
viscosity must remain representable with `0<bulk_rate<2`; zero and underresolved
bulk coefficients refuse rather than acquiring an artificial positive floor.
The native `material_bulk_viscous_heating` tests compare forced pure-shear,
pure-dilation and mixed-strain tensor power against actual thermal transport,
exercise successive temperature/material queries and changes in physical zeta,
and verify endpoint/cancellation refusal without publishing state or selectors.

## Source-local shear and bulk viscosity on 3-D thermal plates

`BulkViscousMaterial<F>` retains one original `FluidTransport` response and an
explicit positive physical `zeta` in Pa s from the same material query. It reads
the original transport owner once; no bulk/shear ratio or numerical damping law
is inferred. `MaterialSlab3::step_with_local_bulk_rheology_and_viscous_heating`
queries this owner at each cell's actual temperature and deviatoric invariant
`gamma=sqrt(2 D_dev:D_dev)`. Shear viscosity, conductivity and bulk viscosity may all
vary by cell and through successive accepted temperature/strain queries.

For Guo-corrected stress numerators `Adev` and `Abulk`, the local root solves
`gamma=Adev/(tau_s dt)`. The reported relaxation defect belongs to this physical
shear-tau root. The same accepted source returns the independent bulk viscosity
and hence `tau_b`; `Abulk/(tau_b dt)` reconstructs the separate dilation for heat.
Collision retains both physical relaxation times. Pure compression queries
zero shear while preserving its bulk relaxation and positive volumetric heat.
The reconstructed deviatoric strain agrees with the queried invariant to the
admitted root tolerance, and heat remains
`q=2mu D_dev:D_dev+zeta(div u)^2` without counting isotropic shear twice.

`PlatesGrid3::set_central_moment_bulk_field` validates and retains a complete
bulk-rate field with a uniform higher-order rate. Checked thermal snapshots
and local shear retuning preserve this field. The heated material step publishes
its applied bulk rates, populations and sensible heat atomically; endpoint source
responses are retained for diagnostics and admission, not applied a second time.
Source-domain, convergence, conversion, heat-balance and cancellation refusals
leave the accepted slab and collision field unchanged. The legacy heated API
refuses both uniform and spatial independent-bulk selections without a physical
binding. A constant field preserves the uniform selector's population arithmetic.

This uses the existing fixed reference density/Cp and explicit thermal update.
Zero or unrepresentable bulk viscosity still refuses; no floor, molecular
relaxation memory or compressible total-energy claim is introduced. Native tests
cover forced pure shear, pure dilation and mixed strain against independent
tensor power and the original thermal solver, spatial/temperature feedback,
constant-field bit parity, and source/cancellation rollback.

## Source-local bulk viscosity in compressible NSF

`NsfGrid3::with_bulk_viscosity` retains a cited `NsfBulkViscositySource` on the
original gas law before the first step. The source sees each query's actual
absolute temperature, thermodynamic pressure, density and molecular shear
viscosity in SI units. It returns positive physical bulk viscosity `zeta` in
Pa s. A cloned characteristic trial shares the same original source. Sources
must be deterministic functions of the supplied state; predictor queries can
repeat and may run in parallel. This is an instantaneous molecular law, with
no inferred eddy bulk viscosity or frequency-dependent molecular memory.

The native D3Q27 trace relation is
`zeta_lu=(2/3)rho theta0(1/omega_b-1/2)`. The material snapshot retains the
source response and its converted trace rate. Collision and the existing
conservative total-energy stress work use that same rate. Shear SGS modifies
only the shear transport: it is not fed to the molecular source and does not
augment its returned bulk viscosity. No second viscous heat source is added.
The dispersion reference uses the source's physical reference-state bulk
viscosity while retaining the original gas, shear and thermal transport.

Accepted endpoint diagnostics retain the final corrected temperature and
pressure, including the last caloric energy corrector. Source, range or rate
refusal leaves the accepted state and material epoch unchanged. The consuming
cold builder rebuilds only the material, stress and conservative flux history
required by its source; it preserves the original populations and total
energy. Without this opt-in source, the existing bulk/shear-ratio route retains
its original arithmetic and allocation behavior. Source values are never
floored or clamped, and the existing gas, Mach and temperature envelopes remain
in force. This extension does not widen the solver's qualified stability range.

Focused native tests exercise local temperature/pressure feedback, accepted
material epochs, molecular-ratio recovery, nonzero shear SGS, source refusal,
and characteristic trial accept/drop/cancel. A two-rung acoustic regression
compares the incremental attenuation caused by added bulk viscosity against
the exact linearized NSF modes, with relative error below `2e-4` on the finer
rung and decreasing under refinement. A transverse shear wave checks that
changing only bulk viscosity does not introduce shear loss. These focused
checks supplement, rather than replace, the broader NSF qualification matrix.

## Molecular viscous heat with thermal feedback on the original free surface

`MaterialFreeSurfaceEnthalpy3::step_with_conduction_and_viscous_feedback`
refreshes molecular stress power at every source evaluation on the original
thermal timeline. It captures the accepted Guo-corrected nonequilibrium stress
numerator once, alongside the existing source-owned mechanical pressure when
pressure-dependent rheology is selected. The original `RheologyModel::sample`
then closes `gamma*tau(T,gamma)=A/dt` at each actual caloric state. Newtonian
sources use the same numerator with their newly resolved molecular tau. Thus
`mu(T,gamma)*gamma^2` uses one consistent material/strain response; it never
combines an old physical strain rate with a new viscosity. The source coordinate
is `sqrt(2 dev(D):dev(D))`. The orthogonal volumetric contribution is retained
separately, with the same local tau for tied collision or its independent bulk
rate when explicitly selected.

Within one private thermal candidate, an identical enthalpy bit pattern can
reuse its complete resolved sample for the same original source, held stress,
pressure and physical scaling. Cancellation still polls before each cache hit.
Only two bounded enthalpy-key vectors accompany the existing diagnostic
samples; the cache never survives a different flow interval or source binding.

Every accepted thermal substep integrates its applied power once over the
located liquid volume and actual duration. Zero-mass interface ghosts and
unlocated conversion carry receive no inferred strain or heat. The report
retains the last applied and thermal-endpoint source samples, including the
physical rate, molecular tau and power. Endpoint samples provide admission and
diagnostics without a second deposit. Authored heat, internal conduction and
molecular heat retain distinct accounting, and the original thermal audit must
close their sum. The original VOF collision, mass routes and conservative
`mass*h` transport then advance once using the accepted thermal endpoint.
Source, thermal-work-budget, native-flow or cancellation refusal preserves the
accepted liquid state, material/rheology observations and physical clock.

This opt-in is explicit thermal subcycling with the original first-order
thermal/flow split. It holds the kinetic stress numerator, not physical strain
or physical Cauchy stress, during thermal evolution. Nonlinear thermal time
refinement is still required; no unconditional stability, simultaneous momentum
feedback, kinetic-plus-thermal conservation, phase heat,
pressure work or numerical-loss heating is asserted. Existing frozen-source
viscous-heating APIs retain their original arithmetic and meaning.

Focused regression sources compare thermal/rheological feedback with an
independent dimensional RK4 heat ODE and closed-form constitutive root, check
substep refinement and original source identity, exercise Newtonian conduction and
constant-viscosity native-population parity, and verify complete rollback and
exact retry after source and late thermal-budget refusals. These are synthetic
constitutive coupling checks, not experimental validation of a liquid.

### Independent physical bulk viscosity on the sparse free surface

`FreeSurface3::set_central_moment_bulk` selects the existing D3Q19 central-moment
operator, retaining the molecular shear rate and authored independent trace and
higher rates. Prepared analytic rheology replaces only shear relaxation.
Material plans select `with_bulk_viscosity(zeta_pa_s, higher_rate)` for a
uniform coefficient;
the shared dense/free-surface conversion is
`tau_bulk = 1/2 + (9/2) zeta dt/(rho0 dx^2)`. This is an explicitly authored
uniform physical scalar, without a coefficient inferred from numerical damping.
The separate source selection below retains a material law.
Zero/unrepresentable viscosity or
rates outside `(0,2)` refuse without a floor. D2Q9 and sparse SGS consumers
explicitly refuse this selection.

Native/P1 initialization uses `tau_s dev(D) + tau_bulk tr(D) I/3` for the retained
nonequilibrium second moment. Frozen and thermally refreshed viscous heating use
`mu gamma^2 + zeta (div u)^2`, with separate diagnostic powers and absolute
divergence; neither carry nor zero-mass ghosts receive heat. Thermal feedback
holds the independent trace rate for a uniform authored bulk coefficient.
The default still retains the original `zeta=2mu/3` arithmetic.
Runtime coefficient setters preserve populations and therefore change the
inferred trace strain; selecting the plan before initialization preserves an
authored initial gradient. This supplies no reversible pressure work, surface
energy, compressible EOS, physical acoustic speed or discrete kinetic-energy
thermalization claim. Focused regressions cover forced moments, source shear
independence, actual collision, transfer heat and thermal bulk/shear separation.

`with_bulk_viscosity_source(property, higher_rate)` instead binds one already
selected positive Pa s property in the original material bundle. The original
card, complete scalar requirement and selected claim remain pinned. The same
temperature, deviatoric-shear and supported pressure query supplies both mu and
zeta; all retained properties pass the existing source-axis admission. Frozen
consumers permit only source support appropriate to their held coordinates.
Missing sources, unsupported dynamics and invalid relaxation refuse.

The sparse native owner applies a complete local bulk-rate field alongside the
local shear rates. Every wet collision cell, including empty interface
carriers, requires its own admitted rate. New wet cells first expose raw stress
numerators to the material query, and the complete source field is installed
before publication. Isothermal, moving-wall and thermal endpoints retain this
same rule. Failed or cancelled coefficient replacement leaves accepted state
unchanged.

P1 transfer validates original vertices and pinned source support lines, then
initializes shear and trace moments with coefficients from each receiver's same
material sample. Thermal feedback retains raw trace stress and recomputes both
zeta and tau_bulk at every fresh material evaluation; it never reuses a stale
trace rate with a new coefficient. Reported bulk viscosity, bulk relaxation and
separate shear/bulk powers identify the applied response. Focused regressions
cover spatially varying source coefficients, the first transferred heat step,
nonlinear thermal feedback, strict local fields and refusal without mutation.

## Implicit conduction inside phase-resistance heat/volume coupling

With `implicit-enthalpy-coupling`,
`FlowField::advance_isotropic_implicit_conductive_phase_resistance_thermal_volume_momentum_with_baths`
solves source-owned Fourier/bath heat within every existing bounded physical
drag/phase/pressure trial. Each trial starts from the same original M, h and P.
Authored net heat plus that trial's physical Darcy-heat guess drives one native
backward-Euler caloric interval; its independently reconstructed endpoint face
heat then drives the original chart, pressure, volume and material transport.
The returned phase resolves molecular viscosity and permeability through the
original material source. Physical-versus-used drag heat and endpoint
constitutive impulse/heat must pass the original coupling gates together.

Only the accepted trial's conductive increments, bath account, nonlinear report,
phase sources, applied resistance and pressure readback are returned. Its native
pressure/drag impulse is scattered once and the accepted mass routes carry
momentum once. Numerical mechanical damping and pressure/source work retain
their separate accounts. Implicit endpoint/assembly defects, native caloric
transport and drag/constitutive heat residuals share both the original enthalpy
and mechanical-energy defect budgets. Initial, incoming, accepted pre-transport
and final transported states retain original source/chart and mobile-phase
admission; fixed-solid voxels own no fluid inventory.

The existing implicit residual-evaluation and total Krylov-iteration ceilings
cover all outer trials, including rejected guesses. Per-direction Krylov cycles
and per-trial Newton/backtrack bounds retain their native meanings, while the
original resistance controls bound the number of trials. The result reports
aggregate thermal work separately from the accepted physical interval. There
is no hidden thermal subcycling or accumulated physical time across trials.
Source, budget, native-solver or cancellation refusal publishes no candidate.
Legacy resistance entries use an identity caloric transform and retain their
existing arithmetic and physical accounts.

This remains the original first-order caloric/flow split: k(T,f) is implicit at
the pre-advection heat endpoint, and phase resistance uses the transported
endpoint. It does not assert simultaneous advection-diffusion, finite bath
storage, exact-solid locking, or a complete fluid-plus-bath entropy audit.
The focused regression source manufactures a mobile latent cell with equal
phase densities, k(f), a nonzero throughflow, physical drag and a one-sided bath.
Independent Fourier, Darcy, material-energy and momentum formulas distinguish
heat inside the implicit endpoint from a late heat deposit; the fixture also
checks whole-call work exhaustion, late source/cancellation refusal and retry.

## Finite solid storage in implicit phase-volume conduction

With `implicit-enthalpy-coupling`, the conjugate thermal-volume and native
momentum entries accept explicit `ThermalVolumeSolidCell` inventories alongside
the original fluid inventory. Every selected fixed-solid voxel retains its
actual held mass, original enthalpy/phase state, authored heat and borrowed
`PreparedPhaseConductivity`. Its own caloric and k(T,f) laws, including analytic
enthalpy tangents and latent plateaus, participate in the SAME sparse implicit
graph as the fluid. Repeated source identities borrowing the same original
caloric object share a native law binding; separately allocated caloric curves
retain distinct bindings even when their immutable contents are identical.

Every face between participating fluid/solid or solid/solid cells uses the
original area and both actual half distances. Optional `ThermalVolumeContact`
selects a fixed `Resistance` or an original constitutive `Source`, adding one
nonnegative area-specific resistance in series; unlisted active interfaces
are perfect contacts. Duplicate/reversed, nonadjacent,
inactive and fluid-only contact assignments refuse. A prescribed bath cannot
also own any face between active thermal cells. Unselected fixed solids and
unselected exposed faces remain insulated; authored solid heat and prescribed
exposed-face baths provide explicit external input. Source contact resistance
and both enthalpy derivatives enter the existing native Newton system.

One backward-Euler interval evaluates both temperatures and source conductivities
at the coupled heat endpoint. Independently reconstructed face heat has opposite
signs in the two receiving inventories. Local storage/physical-heat mismatches,
addition roundoff and whole-graph balance share the original enthalpy-defect
budget with subsequent fluid transport. `ImplicitConjugateThermalVolumeReport`
retains the native whole-graph account, authored solid heat, actual solid energy
change, prescribed solid-boundary heat and signed solid-to-fluid transfer.
Solid endpoint states and source
identities are returned separately, in increasing voxel order. They never
occupy the fluid mass, phase, heat-source or material-advection arrays.

The native momentum variant applies only its original accepted pressure impulse
and fluid mass routes. Both returned fluid and solid states are atomic on
source, thermal, flow, budget or cancellation refusal. Empty solid/contact
inputs retain the existing insulated-solid path. Workspace counts all fluid
cells plus explicitly selected finite solids. Solid geometry, mass and occupancy
remain fixed even if their caloric law crosses a latent interval; this is not
wall melting motion, solid advection or thermal expansion. The original
first-order heat/flow split and its time-refinement requirement remain.

Focused regression sources cover an independently manufactured graded fluid/PCM
wall with unequal half distances, finite contact resistance and k(f), a heated
solid strip conducting through both solid/solid and fluid/solid faces, total
fluid-plus-solid energy, unchanged solid mass, native momentum composition,
empty-wall parity, exclusive face ownership and source/cancellation rollback.

### Finite-solid conduction inside molecular and phase-resistance coupling

The implicit conjugate viscosity entry prepares the original molecular kick
once, then places its physical dissipation and authored fluid/solid heat in the
same finite-wall graph before native volume and momentum transport. Each wall
retains its separate caloric inventory and source conductivity. The accepted
wall endpoint therefore responds to molecular heat during that heat interval.
Molecular backward-Euler/scatter losses and mechanical pressure work retain
their existing separate accounts; they do not become additional heat. The
combined caloric, heat-addition and transport defects use the original shared
enthalpy budget. Source viscosity remains frozen over the original split.

The implicit conjugate phase-resistance entry uses that same fluid/wall graph
inside each native drag/phase/pressure trial. Both inventories start from their
original mass and enthalpy on every trial. Authored wall heat is applied for
one physical interval each time; rejected trials never advance wall history.
The accepted fluid phase, constitutive drag, physical heat and pressure must
meet the existing outer gates before either fluid or wall endpoint is returned.
The accepted wall result travels with that trial's heat payload, while the
whole-call residual/Krylov ceilings account for all attempted thermal solves.

Both entries reuse the same fixed geometry, fixed or source-owned contact
resistance, source bounds and finite-wall energy account as the base solver.
Only fluid heat enters native advection; no solid mass or momentum is routed.
Existing fluid-only methods use empty wall/contact inputs through the same
owners. These paths retain the original first-order heat/flow split and do not
add evolving wall geometry, wall mechanical recoil or viscosity iteration.

Focused regression sources independently close a one-fluid/one-wall molecular
case with Q_mu=8/25 J, Q_wall=46/225 J and wall h=652/225 J/kg, and a latent
fluid/finite-wall phase-drag case with Q_wall=1 J and Q_drag=1/27 J. Nonzero
fluid throughflow leaves the manufactured caloric endpoint unchanged, while
solid transport stays absent. The latter requires multiple coupling trials
and checks late cancellation followed by exact retry from original wall state.

### Prescribed heating and cooling of finite solid walls

`ThermalVolumeBathFace` can select an exposed face of an explicitly bound
`ThermalVolumeSolidCell` in every implicit conjugate entry: thermal-volume,
native momentum, molecular viscosity and phase resistance. The shared admission
owner uses the actual participating thermal cells; an unbound solid receiver,
duplicate face or bath facing any other active cell refuses. Original fluid-only
entries retain their fluid-only receiver rule and empty-wall numerical path.

Each bath uses the RECEIVING cell's original k(T,f), caloric law and analytic
enthalpy tangent. Its conductance is `A/(d/k + R)` with actual face area and
receiving half width, including on graded walls. Distinct hot and cold faces
retain their individual temperatures and series resistances. Their powers and
`d[G(h)(Tb-T(h))]/dh` enter the same sparse endpoint residual/Jacobian as bulk
conduction and finite storage. There is no additional timestep or thermal solve.

`solid_boundary_heat_j` is signed inward heat on finite solids;
`absolute_solid_boundary_heat_j` is the non-cancelling sum of its bath-channel magnitudes.
Both are subaccounts of the existing total boundary report, never extra deposits.
Solid storage closes authored solid heat plus solid boundary heat minus
`solid_to_fluid_heat_j`, within the shared energy-defect budget. Wall endpoint
conductive increments already include their bath heat. Each coupled resistance
trial resets original wall history and retains only its accepted boundary account.
Numerical damping, pressure work, mass advection and momentum ownership are unchanged.

The focused regression source independently manufactures a graded latent wall
with k(f), finite contact resistance, simultaneous +12/-3 J wall baths and a
separate +1 J fluid bath. Its fluid/wall storage changes are 7/8 J with 5 J
authored wall heat. Native momentum composition, original source bounds,
internal-face/duplicate/unbound refusal and exact retry supplement this balance.
An equal-content/distinct-caloric-object case also closes an independent
two-cell balance while preserving the native bath-law identity requirement.

### Original interface sources in finite-wall conduction

All four implicit conjugate entries take `ThermalVolumeContact` selections.
`Resistance(FaceContactResistance)` retains the fixed-resistance path; an
existing fixed selection can be converted with `.into()`. `Source { cells, law }`
borrows an original `ImplicitContactLaw`, including `PreparedMaterialContact`
from the existing material/interface owner. Fluid-only entry signatures remain
unchanged. Every selected face must touch a participating finite solid.

The source retains its original caloric objects, selected response, validity
domain, material claim and fixed external conditions. Its two enthalpy
derivatives include the original temperature/phase chain rule. Native
`with_contacts` attaches it to the SAME implicit conduction graph; geometric
contact resistance is zero on that link so the source contributes exactly once.
There is no second thermal solve or additional contact iteration. If voxel
numbering reverses source A/B order, caloric identities, enthalpy domains,
trial values and both derivative components are reordered together. The
original interface context and source receipt remain unchanged. External
pressure remains the prepared source's pressure; fluid projection pressure
does not silently replace it.

Accepted heat is independently reconstructed with the same endpoint contact
law and actual half-cell geometry, then deposited with opposite signs in the
two inventories. Source bounds and state admission are checked even when a
bulk side has zero conductivity. Rejected trials publish no wall history.
The existing shared energy-defect and whole-call work budgets remain in force.
Fixed-only and empty contact paths skip source attachment and preserve their
existing endpoint resistance arithmetic.

The focused regression source uses an original two-coordinate interface
surface `R = 0.5 + f_wall - p/400000`, with wall A at voxel 1 and fluid B at
voxel 0. At fixed pressure 100000 Pa, unequal masses, half widths and
conductivities transfer 20 J into a latent wall: fluid/wall specific enthalpies
change from 51/25 to 41/30 J/kg while the wall stays at 320 K. Freezing the
initial resistance instead transfers 31/1.5 J. A second fixed-pressure root,
both derivative components, original object order, conservation, native
momentum forwarding, mixed contact duplicates, source-domain refusal and
cancellation/retry supplement this physical distinction. These are synthetic
constitutive regression sources, not experimental validation. Geometry,
interfacial storage and external mechanical conditions do not evolve here;
the original first-order heat/flow split remains authoritative.

### Resolved exposed-surface conduction in phase-volume flow

All four implicit conjugate entries accept `ThermalVolumeBoundaries`, borrowing
fixed `baths` and original per-face `sources`. Existing fixed-bath slices and
arrays convert with `.into()` without allocation. Fluid-only public signatures
remain unchanged. Each `ThermalVolumeBoundaryFace::Law` borrows a native
`ImplicitBoundaryLaw`; `Exchange` borrows an `ImplicitBoundaryExchangeLaw` for
two distinct baths sharing one physical surface. The receiving cell retains
its exact original caloric object, actual graded face area and half-width.
Mixed fixed/source duplicates, occupied faces and unbound receivers refuse.

This admits the existing `ImplicitLayerBoundary` source-integrated pure-phase
k(T) half-cell/film balance, and native `GrayAmbientBoundary` or
`GrayTwoBathBoundary` surface solves. Their physical assumptions remain distinct:
the integrated layer enforces its material-temperature support through the
half-cell; gray owners use cellwise k(h), its derivative and optionally original
enthalpy-dependent emissivity/HTC while solving the actual exterior temperature.
Film resistance belongs to the original boundary law. For gray owners the
descriptor's resistance is a coating/contact before the exterior surface;
the convective HTC owns its own film and is never counted again there.

The same per-cell source aggregation supplies every selected face's complete
inward heat and analytic enthalpy derivative to the existing sparse implicit
graph. Original source callbacks retain their bounded work and cancellation
checks. The bridge adds no new timestepper, face storage or thermal iteration.
Accepted heat is independently reconstructed by evaluating the same face laws
at the native endpoint. Each phase-resistance trial starts from the original
wall/fluid inventories, and only the accepted face exchanges are returned.

`ImplicitConductiveThermalVolumeReport::bath_exchanges` retains cell, face,
original channel order, actual reservoir temperature and signed heat into the
body. Directional and solid boundary totals remain net physical heat;
absolute totals sum individual channel magnitudes. Opposing convection and
radiation on one surface therefore remain visible even at zero net wall heat.
These are subaccounts of the existing energy balance, never extra deposits.
The native local-source entropy account remains explicitly unclassified; no
complete body-plus-bath entropy claim is made by this adapter.

Independent regression sources manufacture: a graded nonlinear material wall
whose integrated X-face removes 87.5 W while a separate Y-face supplies 20 W;
a 400 K bulk cooling through a 350 K physical surface with 200 W total gray
exchange; and one shared surface exchanging +10/-10 J with distinct 400/300 K
baths at zero net storage change. The mixed-face case also checks original
caloric identity, unsupported material-face temperature, exclusive ownership
and cancellation/retry. Held geometry, fixed external conditions, original
source domains and the first-order heat/flow split remain authoritative.

## Fitted dynamic-SGS coefficients

`PeriodicDynamicSgs` and `ChannelDynamicSgs` expose the admitted nonnegative
`C_s²` independently of cell-centered viscosity, through an indexed accessor
and an allocation-free iterator. Outputs are consumable only after successful
preparation; an out-of-grid index returns `None`. The completed stress-basis
scratch stores these coefficients, so numerical storage budgets are unchanged.
Spatial and Lagrangian fits retain their original viscosity arithmetic; a
uniform-flow preparation clears the coefficient cache. Regression sources
compare the coefficients against the existing independent full-tensor Germano
oracle and check that scalar preparation preserves them without extra storage.
The coefficient remains an exploratory model estimate, not DNS validation.

`PeriodicDynamicScalar`, `ChannelDynamicScalar` and their common-component
wrappers likewise expose fitted `C_phi` through `coefficient` and
`coefficients`. The completed scalar-basis scratch retains the coefficient
without allocating another field or changing the reported cell diffusivity.
This lets conservative face consumers apply their own resolved strain even
when the adjacent centered cell strains vanish. Uniform-flow or uniform-scalar
preparation clears the coefficient output, including when reusing a workspace.
The existing successful-preparation, bounds and exploratory-use rules apply.

`sgs3::smagorinsky_viscosity` accepts `C_s` and the physical `sqrt(2 S:S)`
invariant, allowing conservative faces to reuse one strain evaluation across
independent fitted coefficients. A fitted `C_s²` or `C_phi` is supplied as its
square root. The existing `(C_s Delta) * ((C_s Delta) * |S|)` arithmetic and
input validation are unchanged, including validation at exactly zero transport.

## Resolved inviscid vorticity production

`FlowKinematics3::inviscid_vorticity_budget` retains the separate stretching
`G omega`, compressible dilatation `-omega div(u)`, and baroclinic
`grad(rho) cross grad(p) / rho^2` contributions. Pressure is thermodynamic
pressure in the same coherent units as density and velocity, never p/rho.
`VorticityBudget3` exposes their sum and the corresponding local inviscid
enstrophy rate. Physical conversion uses the actual timestep: curl scales
with 1/dt, the three rates with 1/dt^2, and enstrophy production with 1/dt^3.
No reference density or guessed cell spacing enters that conversion.

`NsfGrid3::vorticity_budget` samples the original accepted populations,
current half-force velocity and accepted p=rho*theta. Seven fixed-size
population reductions supply the centered fluid stencil without allocation,
source requery or state mutation. Periodic neighbors wrap through the native
grid. A solid or nonperiodic face interrupting the scalar stencil refuses;
the ordinary boundary-aware kinematics field remains available there. This
query supplies no invented pressure/density wall condition.

These are resolved inviscid material-rate terms, not an additional momentum
force or a complete vorticity transport discretization. Viscous stress curl,
nonconservative applied-force curl, explicit source injection, and subgrid
terms are not included. Local bulk enstrophy is not a global boundary flux
balance or an uncertainty certificate.

The module tests cover expansion/rotation, barotropic cancellation,
pressure-density misalignment, mass-unit invariance, SI conversion, invalid
data, periodic agreement with native force-centered curl, and incomplete
thermodynamic stencils (G0/G1/G3). Independent analytic checks were performed;
native Rust execution was unavailable when this extension was introduced.

## Native lattice channel trajectories with turbulence memory

`Grid3::channel_lagrangian_sgs_workspace` binds the unique fixed halfway-wall
normal and reuses the original byte-admitted local collision fields and two
history banks. The existing serial/pooled Lagrangian step methods now accept
stationary all-fluid channels as well as periodic domains. Explicit
`step_serial_with_channel_lagrangian_sgs_and_moving_walls` and pooled
counterparts bind the complete native reflected-link probe and one uniform
tangential velocity per wall to the same accepted step.

The current wall velocities enter the shared Favre contractions, native
halfway-wall strain diagnostics and original moving bounce-back. Only
tangential filtering/backtracking wraps periodically; normal histories never
connect opposite walls. Fits and pathline departures use the incoming
Guo-centered physical velocity. Original molecular collision fields and
D3Q27 independent trace/higher rates remain authoritative.

The native owner must publish populations, lattice clock and any wall
impulse/work receipt before the infallible history commit. Preparation,
collision, streaming, wall-drive, Mach or cancellation refusals retain the
accepted history and native state. A workspace cannot silently change wall
normal, become periodic or admit interior solids. No additional field
arrays or per-preparation allocations are introduced.

Six native regression sources cover both stencils, all wall normals,
independent shared-history/native-transport references, molecular fields,
serial/pooled identity, moving-wall impulse/work, Couette/comotion,
storage/topology refusals and late-step rollback. Changed sources passed
pinned Rust formatting/parsing and bounded production review; these native
consumer tests were not compiled or executed here. The shared core's prior
74-test run is not a native-consumer test claim. Authored memory times,
first-order trajectories and exploratory-model limits remain.

## Lagrangian SGS with spatial and changing body forces

`ForcedGrid3` composes the original Lagrangian SGS workspace with spatial
and time-dependent force fields in periodic domains and stationary channels.
`channel_lagrangian_sgs_workspace` binds the physical wall normal, and the
serial/pooled `step_with_channel_lagrangian_sgs_and_moving_walls` and
`advance_with_channel_lagrangian_sgs_and_moving_walls` families also compose
moving walls with held or next force fields. Fits, pathline departures and
pre-step dissipation use the incoming per-cell Guo half-force velocity and
the prescribed physical wall velocities.

The original forced-wall transaction admits and publishes any next force
together with populations, clock and wall impulse/work. Turbulence history
commits only after that native step succeeds. Refused drive, wall data,
departure or late cancellation leaves the accepted history retryable.
Molecular collision fields and D3Q27 nonshear rates remain under their
original owners. The adapters add no preparation arrays or history banks
and retain the original workspace byte gate.

Explicit memory time, departure limits, resolved-wall treatment and
exploratory-use restrictions remain. This adds no open-boundary history,
scalar-flux or turbulence-validation claim. `tests/lagrangian_forcing.rs`
retains the periodic regression sources; its `channel.rs` child adds four
regression groups covering independent pathline/native owners, held and
changing force epochs, moving-wall receipts, worker identity, endpoint
refusal and post-stream cancellation/retry. The new channel sources were
parsed and reviewed, but have not been compiled or executed.

## Compressible gas gravity and mechanical work

`NsfGrid3::with_gravity_si` and `set_gravity_si` bind uniform signed SI gravity
to the original compressible momentum and total-energy equations. Every force
epoch uses the evolving lattice density in `rho*g`, including characteristic
and thermal correctors; no Boussinesq reference-density subtraction or
caller-authored expansion coefficient is inserted. The exact endpoint
trapezoid of `rho*g.u` supplies total energy and the separate
`EnergyLedger::gravity_work_j` account. This work is not a heat source.
Characteristic trials retain the binding and aggregate its work before their
ordinary energy and cancellation checks.

Changing gravity rebases the complete Guo half-source at fixed physical
velocity and stress, refreshes the original conservative flux history and
publishes only after admission. It advances no clock and supplies no impulse
or heat. Zero gravity on an unforced grid is an exact no-op. Conversion is
performed once at binding; force and energy evaluation add constant work per
cell without another field allocation or material query.

Seven regression definitions cover uniform acceleration and no artificial
heating, density scaling, evolving-density gravity with mass supply, separate
heat/work accounting, evolved stress-preserving rebinding, exact zero parity,
and refused/cancelled characteristic publication. Native compilation and Rust
execution remain outstanding; these definitions are not an executed
qualification. Existing Mach, source/material, boundary and timestep limits
remain. The gravity binding alone claims no hydrostatic well-balancing,
exact discrete fluid-plus-potential energy conservation, self-gravity or
measured performance.

## Curved-wall local SGS on native D3Q19

`BoundaryGrid3::prepare_curved_sgs` and
`step_with_curved_walls_and_sgs` admit explicit Smagorinsky, WALE or Vreman on
frozen BFL geometry. Actual axial intersection distances and frozen per-link
wall velocities supply nonuniform quadratic velocity gradients. The same
velocities drive the native curved population and load transaction. Physical
Guo-centered fluid velocity, original molecular collision, non-shear collision
controls, wall mass/force/torque/work and external transport retain their
existing owners. Sparse ordered traversal reuses the existing SGS workspace.

Closure and load refusal precede population publication; zero coefficient
retains native curved-step population and load bits. This supplies resolved
gradients at the admitted fixed geometry, without a wall law, a nonlocal
curved-domain dynamic filter or evolving geometry. The focused
`tests/curved_sgs.rs` cases cover displaced linear/quadratic walls, rotation,
force centering, five collision families, load accounting and refusal/retry.

## Native staggered-edge SGS stress

The transient finite-volume SGS owner evaluates reciprocal shear-edge
viscosity from the actual two staggered velocity differences at that edge.
Other full-gradient entries and the authored scalar or directional filter
width fields are interpolated using physical graded-grid distances. A single
cached kinematic viscosity supplies both reciprocal stresses; material jumps
retain the original sorted harmonic transmission of each adjacent TOTAL
viscosity. Smagorinsky, WALE and Vreman keep their shared constitutive kernels.

Prepared full-gradient and edge arrays are reused on every nonlinear sweep.
Normal stress and scalar heat retain their cell-interior closure values.
Four-fluid-cell interior edges use the new stencil; walls and incomplete
stencils retain the original boundary treatment. Each adjacent material
independently constrains the edge eddy/molecular ratio during both refresh
and endpoint admission. Zero coefficient keeps the original disabled path.
The SGS unit regressions cover a centered-strain null with nonzero native
momentum stress, full affine tensors on graded grids, reciprocal coefficients,
material endpoint refusal and cache reuse. These are discretization checks,
without a developed turbulent-spectrum or wall-resolution claim.

## Nonlinear periodic capillary jets

`jet::nonlinear::PeriodicSlenderJet` evolves a uniform periodic grid of
cell-centered squared radius and staggered axial momentum using the
Newtonian Eggers-Dupont slender equations with full axisymmetric curvature.
Conservative mass and capillary/viscous/momentum flux differences preserve
periodic liquid volume and axial momentum to floating-point roundoff.
Physical cross-sectional area and volume flux retain their factors of pi.

Minmod-limited MUSCL momentum reconstruction improves smooth translating
disturbances while retaining mass-consistent dual-cell fluxes and nonnegative
numerical kinetic loss. The moving-wave G1 test checks 64-to-128 refinement
against independent pressureless transport; limiting remains active at extrema.

The owner uses three-stage SSP Runge-Kutta, adaptive advection/viscous/capillary
step limits and reusable O(N) buffers. Each caller supplies a maximum timestep,
CFL factor, attempted-substep budget, neck radius, maximum slope and minimum
radial resolution. A model-domain stop publishes only admitted substeps.
Cancellation, invalid state or exhausted work budget leaves the entire
pre-call physical state and time unchanged. Polling is bounded by 256 cells.

The result exposes the evolved radii/velocities, volume, momentum, surface and
kinetic energy, integrated physical viscous loss, separate upwind numerical
loss and a signed energy-balance remainder. The full spatial discretization
does not claim exact energy conservation or a certified error bound.
Focused G1/G3 tests check the small-amplitude Rayleigh limit and refinement,
nonlinear conservation, Newtonian extensional damping and SI similarity;
G0/G4 cover model stops, refusal, rollback and deterministic retry.

`jet::material::nonlinear::MaterialPeriodicSlenderJet` retains the exact
resolved liquid source beside the evolving native state. Density, kinematic
viscosity and tension are read from that same point with its original receipts.
Both original viscosity bases are supported. Omitted active surface data and
sources requiring evolving non-Newtonian rheology refuse construction.
The shared `SlenderJetEvolution` trait lets downstream studies retain this
source-bound owner without duplicating the PDE or exposing coefficient mutation.

This is a constant-property incompressible Newtonian periodic cylinder with
a passive exterior. It does not supply nozzle/end conditions, gravity-induced
stretching, variable surface tension, non-Newtonian stress, active gas dynamics
or post-pinch topology. Reaching an authored neck/resolution limit is not a
resolved breakup or a measured droplet-size prediction.

## Isothermal hydrostatic initialization of compressible gas

`NsfGrid3::new_isothermal_hydrostatic` creates a resting column with signed
SI gravity along one declared grid axis. The reference pressure belongs to
the declared reference cell's center. At the uniform reference temperature,
adjacent densities obey
`rho[j+1]/rho[j] = (1+b)/(1-b)`, `b = g_lattice/(2 theta0)`.
Equivalently the pressure jump equals gravity times the trapezoidal density
times the cell width. Combined with the existing full Guo half-source chart,
this is a fixed point of the original D3Q27 collision, streaming and stationary
halfway bounce-back operators up to floating-point roundoff. It approaches
the ideal-gas barometric exponential with second-order spatial error.

The gravity axis must be nonperiodic unless gravity is zero. Transverse
periodicity, stationary domain walls and embedded solids are supported;
thermal walls must be adiabatic or at the reference temperature. The ordinary
grid admission and original gas/material pressure and temperature domains
remain authoritative. Nonfinite gravity, nonpositive or overflowing profiles,
and `|b| >= 1` refuse. No clipping or altered equation of state is introduced.
The temporary layer profile adds no persistent field or step-time overhead.
Zero gravity delegates to the original uniform constructor with identical
state bits. Subsequent changes in gravity or boundary forcing evolve normally.

Six regression definitions cover every axis and gravity sign, pressure datum
and force balance, resting walls/solids and worker identity, second-order
barometric convergence, release after gravity removal, source-gas pressure
support, zero parity and invalid columns. Rust syntax/formatting checks and
independent discrete algebra checks were run; full crate execution remains
outstanding. This initialization is not a general well-balanced perturbation
scheme for arbitrary gravity directions or nonuniform temperature, and it
adds no fluid-plus-gravitational-potential energy conservation claim.

`NsfGrid3::new_isothermal_hydrostatic_profile` extends that resting state to
one signed SI acceleration per layer, including spatial variation and force
reversal. Adjacent layers obey
`rho[j+1]/rho[j] = (1+b[j])/(1-b[j+1])`, so the pressure jump balances the
trapezoidal **force density** `rho*a`, with each endpoint's own acceleration.
The existing body-force owner supplies the physical Guo half-source and
reports its work in `body_force_work_j`; uniform gravity stays disabled.
No independent pressure evolution, forcing operator or profile sampling is
added to the time step. The original cell cap bounds the retained body field;
the all-zero profile delegates to the original uniform constructor. A constant
profile reproduces the uniform-gravity pressure and evolution, with the work
assigned to the corresponding body-force account.

A periodic force axis is supported when its discrete pressure increments
close around the entire loop. Compensated summation of
`log1p(b[j])-log1p(-b[j])` checks that compatibility at floating-point tolerance;
zero arithmetic mean acceleration alone does not establish it. An incompatible
loop, nonfinite/unrepresentable layer force, `|b[j]| >= 1`, an invalid profile
shape or an unsupported material state refuses. Six additional physical
regression definitions cover varying/reversing forces along every axis with
stationary walls and embedded solids, periodic seam balance, second-order
convergence to an independently integrated varying-force potential, pressure
release after removing support, constant/zero parity and invalid profiles.
The same isothermal, axis-aligned scope and potential-energy limits apply.

### Held spatial body acceleration and vorticity production

`NsfGrid3::with_body_acceleration_si` and `set_body_acceleration_si`
hold one prescribed SI acceleration at each native cell center. The original
compressible density supplies `rho*a` at every thermal/characteristic/SGS
force epoch; the original energy corrector adds the same endpoint trapezoid
of `rho*a.u` reported by `EnergyLedger::body_force_work_j`. Gravity and
manufactured sources remain separate supplies. Mechanical work cancels its
kinetic contribution before the internal-energy remainder; it is not added
as heat or subtracted from internal energy a second time. Characteristic
trials retain the field and its work through copy, reduction, refusal and
acceptance.

The held field is x-fastest, bounded by the original grid's cell cap, and
requires 48 numerical bytes per cell when enabled. Solid entries are ignored.
An empty slice or all-zero fluid field disables it; no field storage or
additional floating-point force/energy arithmetic is used on the disabled
route. Per-step evaluation performs no field allocation or callback. The
setter shares the gravity owner's full Guo half-source rebinding at fixed
accepted velocity, density and physical stress, then refreshes the original
energy-flux history. SI conversion, field-length, finite-state and cancellation
refusals precede publication. The setter polls every 64 input/rebinding cells,
around native history reconstruction and immediately before publication; that
existing reconstruction is not internally interruptible. The accepted SI
field participates in the state hash.

`body_acceleration_curl_si` returns the held field's resolved `curl(a)`
in s^-2 using the native fluid-neighbor topology and physical spacing.
It does not use `curl(rho*a)/rho`, which would invent a density-gradient
source for uniform acceleration. Missing fluid neighbors and two-cell
derivative axes refuse; a one-cell periodic axis denotes a homogeneous
extrusion. The existing inviscid budget stays separate. Manufactured
mass/momentum source curls, viscous/SGS stress curls and a full discrete
vorticity balance are outside this point query.

This is an externally prescribed acceleration held during each native step,
not a velocity-dependent Coriolis/drag integrator, a potential-energy
conservation method or a new turbulence closure. The body-force regressions
exercise uniform equivalence to native gravity, separate source-work
composition, evolving density under mass injection, actual vortical onset
and spatial refinement, stress-preserving rebinding, cancellation/invalid
inputs, stencil refusals, and characteristic trial/drop/accept behavior.

### Joint Lagrangian momentum and scalar memories

sgs3::dynamic::lagrangian::scalar::LagrangianDynamicScalar owns one momentum
trajectory and an independent scalar Germano history. Both fits reuse the same
current Favre fields, but transport separate signed numerator/denominator pairs
with independently authored memory times. Scalar physical contraction scaling
uses rho^2 u^2 phi^2, distinct from momentum rho^2 u^4. One trial commits both
histories and their common sample time, or drops both; uniform current sources
still advect and decay accepted history. The original scalar contraction kernel
and a shared pathline/relaxation/normalization implementation supply the math.
Fitted coefficients remain available independently of centered strain, reusing
the consumed scalar vector-basis storage. Prepared numerical payload is 264
bytes/cell, with checked required bytes and actual capacity admission.

Six new focused tests and eight existing momentum-history tests compiled and
passed in a bounded native source probe (14 passed, 0 failed, 0 ignored). That
probe used the actual dynamic/filter/history modules and required pinned fs-math
deterministic functions; it is not a Cargo crate build or DSR gate. Rustfmt and
bounded source review also completed. The original periodic entry remains
first-order with explicit authored timescales, nonnegative emitted transport
and signed history. The resolved-channel extension below adds physical-wall
histories; adaptive-timescale calibration, inflow memory, scalar variance,
reaction-rate, SGS energy and DNS/turbulence-accuracy qualification remain absent.

### Joint Lagrangian stress, heat and component histories

`sgs3::dynamic::lagrangian::scalar::components::LagrangianDynamicComponents`
adds a common same-unit component history to the original paired momentum/heat
owner. Every component, including closure and absent species, contributes to
one summed Germano objective with a common variation scale. Separate physical
`rho^2 u^2 phi_common^2` history scaling prevents trace components or changing
field amplitudes from acquiring spurious weight. All three fits reuse one
current Favre preparation and the original scalar and pathline kernels.

`LagrangianComponentsConfig` supplies an explicit positive memory time in the
sample-time units, a per-axis departure limit in cells, a coefficient bound
and exploratory-use policy. Each fit retains its own memory time and signed
numerator/denominator history. Zero current sources still advect and decay
accepted memory; emitted coefficients remain nonnegative and independent of
centered strain. The native face owner supplies its actual current strain and
physical diffusion admission.

One exclusive trial exposes all three coefficient sets and reports. Its
infallible commit publishes every history, logarithmic scale and the common
sample time together. Refusal, cancellation and dropped trials preserve all
accepted banks. Numerical payload is 376 bytes/cell, independent of component
count; filtering work is linear in cells times components, with bounded
polling and no transposed component field or per-preparation field allocation.

Six added regressions cover scalar equivalence, component order/units/trace
weighting, independently reconstructed physical-history relaxation, joint
rollback/retry, storage/admission and centered-null/uniform-source memory.
All 68 real-source SGS tests passed with pinned nightly-2026-07-06 and the
actual compiled fs-math dependency. This bounded probe is not a full crate
build or DSR gate. First-order pathline integration and explicit authored
timescales remain. The resolved-channel extension below adds physical-wall
domains; inflow memory, adaptive calibration, reaction-rate and SGS-energy
closure are not supplied.

### Resolved-channel Lagrangian turbulence histories

The momentum, paired scalar and common-component owners now provide
`new_channel`, `prepare_with_walls` and `normal_axis`. They reuse the
original channel Favre filter, Cartesian spacings and physical wall
derivatives. Current no-slip velocities and scalar Value/ZeroGradient/Interior
conditions enter normalization and both physical and filtered gradients.
Common species retain interior-only normal reconstruction; impermeability
does not invent a mass-fraction wall gradient.

Startup averages signed contractions in a local tangential 3-by-3 stencil.
It does not reuse instantaneous whole-plane averaging: each history retains
its local trajectory. Common-component contractions remain raw until the
chosen instantaneous or Lagrangian averaging owner consumes them.

History interpolation remains periodic only in the two tangential axes.
Normal backtracking must remain inside the actual half-cell walls. Positive
trilinear weights interpolate fluid-center histories; within either outer
half cell, constant extension uses the nearest fluid history. A departure
through a physical wall refuses the trial instead of wrapping or reflecting
history across the channel. This is an explicit first-order boundary
approximation, not a calibrated wall model.

All three histories share the original exclusive publication transaction,
with independent signed numerator/denominator scales and memory times.
Cancellation, wall crossing or coefficient refusal preserves every accepted
bank and sample time. Periodic entry points and their arithmetic remain.
Numerical payload remains 168/264/376 bytes per cell for momentum, paired
scalar and common-component owners; channel preparation adds no field
allocation or species-count-dependent storage.

All 74 real-source SGS tests passed using pinned nightly-2026-07-06 and
the actual fs-math dependency, including existing periodic and instantaneous
channel cases. Six added groups cover local tangential initialization,
all normal orientations, wall values and physical scaling, no opposite-wall
history, signed decay/translation, storage admission and joint rollback.
This is a bounded native source probe, not a full crate/DSR gate, turbulent
channel DNS qualification or adaptive memory calibration.

## Dense source-resolved shear rheology

`MaterialFreeSurfaceRheology2` binds the original `IsothermalRheologyLaw`
and bounded constitutive root to an existing dense BGK `FreeSurface`.
Incoming Guo-corrected stress closes physical shear against the selected
source viscosity. Temperature, reference density and uniform sourced
capillarity remain fixed. Source selections are pinned by the same resolver
used by sparse isothermal material flow.

Two reusable source buffers retain the most recent accepted collision input,
before streaming and interface conversion. Their samples are not relabeled
as final-field constitutive evaluations. The admitted tau field enters the
original conservative momentum candidate through one preparation hook,
without a second lattice clone or an alternative collision/VOF integrator.

Source, root-convergence, declared-speed, native-flow and cancellation
refusals publish neither the supplied liquid nor the accepted receipt epoch.
Pressure-evolving sources and non-BGK collision models refuse. There is no
viscosity clipping, analytic material fit, temperature evolution or extra
endpoint source sweep. The original law retains sourced surface tension;
its internal shear-only reference omits that field, so dense scaling uses
the same explicit capillary override as the sparse owner.

Focused regressions cover constant-source native parity for both contact
models, an independently derived quadratic rate/viscosity root that changes
the native collision, and late source-domain refusal preserving flow and
prior receipts. These are synthetic constitutive checks, not measurements
or a whole-scheme nonlinear stability certificate.

### Time-dependent spatial body acceleration on native NSF steps

`NsfGrid3::step_with_body_acceleration_si` accepts the prescribed SI field at
the NEXT native time level. Collision, OLD velocity and incoming population
history retain the accepted force. Every NEW thermal/characteristic/SGS
corrector instead uses the endpoint field multiplied by its actual NEW
density. Acceptance retains that field for ordinary subsequent steps. This
is distinct from `set_body_acceleration_si`, which rebinds the force chart
at fixed physical state without advancing time or applying an impulse.

For physical momentum `m=rho*u`, the original total-energy corrector and
`body_force_work_j` both use the second-order centered work
`(a_old+a_new).(m_old+m_new)/4` in lattice units. This continuously reduces
to the held quadrature. For a spatially uniform drive in uniform constant-density flow,
the original Guo update supplies `Delta(m)=rho*(a_old+a_new)/2`, so the work equals the kinetic
change even for unequal reversal, enabling or disabling. A trapezoid of
endpoint powers would instead add `Delta(a).Delta(m)/4` as artificial heat.
Variable-density transport and manufactured mass/momentum sources retain the
original unresolved mechanical/internal partition; no exact local kinetic
partition is claimed for those cases. Gravity and source energy remain
separate, included once with the same accepted endpoint velocities.

Input conversion/admission occurs once. The existing force and energy owners
receive borrowed OLD/NEW slices; no solver copy, acceleration callback or
extra manufactured-source sample is introduced. One admitted endpoint vector
is staged until publication; the old field is then replaced. Identical
endpoints use the original held arithmetic, and two disabled fields retain
the original zero-field path. Cancellation polls input chunks, native step
stages and corrector cells, with a final poll after all native refusals and
before any state/control publication. External callback side effects are
not rollback state. The caller must resolve temporal variation on the native
step; this is not an implicit velocity-dependent or rotating-frame solver.

G0/G1/G3 regressions cover per-step kinetic work/no heat, force reversal and
enable/disable, unchanged held bits, quadratic-time second-order refinement,
actual spatial vortical onset and retained vorticity after disabling,
evolving density and source sampling, enthalpy/SGS closure and worker parity,
and endpoint refusal/cancellation followed by deterministic retry.


## Molecular-kinetic planar contact-line response

The existing source-bound 2-D, sparse 3-D and thermal planar wetting owners
can select `MolecularKineticLaw` as an alternative to the Cox closure.
It inverts the full Blake-Haynes relation, with explicit calibrated jump
length/frequency and explicit speed, temperature and angular domains.
The selected static advancing/receding endpoint remains a zero-speed source
query. The computed dynamic angle is a model prediction, not a new measurement.

The named absolute-temperature coordinate and sourced surface tension belong
to the original liquid. Thermal preparation and endpoint sampling use the
actual local liquid temperature and the same sigma as native capillary
pressure; genuinely dry numerical ghosts retain the reference context.
Jump kinetics stay constant only over their declared temperature interval.
No Arrhenius response, viscosity-based calibration or gas-wedge matching is
inferred. Measured rate curves, an existing Cox law and duplicate molecular
laws refuse composition. The observed speed must remain in the admitted
domain even inside the numerical pinning deadband.

The shared history owner supplies the signed geometric contact speed and
installs the predicted angle in the original wall ghosts, curvature and
pressure reconstruction. Static-source and native-step rollback remain
unchanged. Local rate-dependent excess force relative to the chosen endpoint
and its nonnegative force-times-speed loss are reported per contact-line
length; neither is an evolved heat source, a total pinning-loss account, or
a closed surface/flow energy balance. A pinned branch consumes model speed
zero while retaining the measured geometric speed separately.

The allocation-free scalar tests check the independent forward sinh law,
linear-friction limit, advancing/receding reflection, exact stationary
endpoints, temperature/tension/calibration changes and domain refusal.
Integration regressions exercise actual original-source native wetting and
thermal context. These are model/numerical checks, not experimental calibration
or a resolved microscopic contact-line claim.


## Optional implicit-viscosity capillary-jet evolution

The original `PeriodicSlenderJet` and `MaterialPeriodicSlenderJet` can select
`SlenderJetIntegrationMethod::ImplicitViscosity`. The explicit SSPRK3 policy
remains the default. Selection preserves physical state, time, density, tension,
viscosity and source identity; cancellation refuses without changing policy.

The four-stage third-order IMEX-SSP3(4,3,3) method advances the SAME conservative
area/momentum operator. Newtonian extensional viscosity is solved implicitly
in face velocity with an O(N) cyclic tridiagonal LDL/rank-one solve. Twelve
additional N-element work arrays are allocated only when first selecting this
policy and are reused across substeps. The solve retains the uniform-velocity
null mode, restores momentum to roundoff and checks the componentwise residual.

Removing the explicit viscous restriction does not remove capillary, advection,
drainage, physical-domain or authored maximum-step limits. The implicit tableau
is L-stable but not stiffly accurate; unresolved stiff transients do not acquire
uniform accuracy. The optional temporal report gives the maximum accepted-step
radius/velocity difference from an embedded second-order estimate in SI units.
It is not a global bound or an automatic tolerance controller.

Positive final stage weights integrate physical and numerical loss rates
separately. The signed energy residual retains temporal/spatial truncation;
implicit damping is not reported as molecular heat or an unconditional
nonlinear energy theorem. Whole-call cancellation, numerical errors and work
exhaustion roll back the accepted state; physical-domain stops keep admitted
substeps. Existing sourced jets and generic resumable vessel flight use this
policy without a duplicate equation or detached material trajectory.

Native regressions check a manufactured nonuniform periodic solve, resolved
third-order viscous decay, high-viscosity agreement/work against an explicit
reference, the inviscid limit, conservation, policy/call rollback and domain
stops. Material regressions keep the original dynamic/kinematic viscosity
sources and compare the same native trajectory. Model restrictions remain
axisymmetric, periodic, slender, Newtonian and pre-pinch-off.

## Native free-surface loads on submerged solid objects

`FreeSurface3::step_with_submerged_solid_loads` advances the original VOF
owner once and returns per-object pressure and viscous reaction from its
actual outgoing/reflected populations. `submerged_solid_load_candidate`
retains the same result unpublished for coupled acceptance. Each object is
a selected group of original permanently solid tiles; several objects share
one fluid step. Native `StationaryWallMomentum3` carries the impulse on each
object and its moment at actual halfway links about the declared origin.
The liquid receives the opposite impulse and angular impulse. Stationary
boundaries do zero mechanical work; no caloric deposition is inferred.

Every exposed D3Q19 link must remain surrounded by full liquid before and
after the step. Gas, interface cells, another touching object, domain-wall
contact, duplicate tiles and missing solid geometry refuse. Complete closed
surfaces cancel a uniform pressure datum, so resultant force and torque need
no invented absolute pressure. Canonical compensated sums retain small
buoyancy and drag contributions against the opposing background pressures.
Explicit tile/link and relative mass-defect bounds, physical low-Mach checks,
finite loads and cancellation complete before candidate publication.

The material entries retain the original source, dx/dt scaling and scenario
speed bound, and finish SI force/torque and mass admission before committing
the same native candidate. Auxiliary reaction storage is proportional to
selected tiles and objects; no second collision or population buffer is added
by the measurement. The ordinary unobserved stream keeps its arithmetic.

Five native and two material regression definitions cover Archimedes force,
center-of-buoyancy torque, uniform-pressure cancellation, analytical initial
drag, continued free-surface trajectory parity, multiple objects and joined
tile assemblies, shifted origins, SI scaling/source identity, and refused
or cancelled publication. All seven changed Rust files passed rustfmt parsing;
the four new source/test files were formatted. Full crate compilation and
execution remain outstanding because DSR/RCH are unavailable and the shared
filesystem is full. These are reaction ports for fixed, fully flooded voxel
geometry, not complete MR69 moving-boundary coupling, partial-wetting geometry,
body inertia, finite-speed piston work or a fluid kinetic/thermal energy certificate.

## Isothermal initialization of material thermal flow

Both material thermal-slab constructors also admit equal positive wall
temperatures. They initialize every thermal population and both walls at
`theta = 0`, retaining `T = T_initial + theta * 1 K` and
`g_beta = g beta (1 K) dt^2/dx`. The one-kelvin coordinate scale supplies only
unit conversion; it imposes no wall contrast or initial thermal load. Nonzero
initial contrasts retain their original signed normalization. Source-owned
constructors normalize the actual material reference temperature in this same
coordinate.

The retained original wall coordinates travel with all material snapshots and
remain the admission authority for source, wall, local-transport, tensor and
electrothermal steps. Conservative material refresh checks physical cell and
wall temperatures without treating the coordinate scale as a boundary
difference. Existing heat-source, wall-flux and film conversions therefore
remain valid from an isothermal start. Nusselt remains singular at zero wall
contrast; Rayleigh initially vanishes for equal prescribed temperatures. The
`material_buoyancy` regressions cover isothermal 2-D/3-D startup, adiabatic source
energy and physical temperature rise, and tight-budget conservative material
refresh.

### Source-resolved free-surface rheology with bound SGS

`MaterialFreeSurfaceRheology3::with_sgs_si` binds one immutable exploratory
closure in `MaterialFreeSurfaceRheologySgs3`. The original material source and
physical units remain authoritative. Smagorinsky trial-rate viscosity or
WALE/Vreman full native velocity-gradient viscosity joins the molecular
relaxation inside the existing bounded constitutive root. The source query is
still the canonical Guo-corrected deviatoric stress invariant. Active gradient
closures retain the original occupied-liquid second-order stencil and workspace
admission. Zero coefficient delegates the original molecular trajectory and
needs no gradient workspace.

The same accepted material response supplies independent molecular bulk
viscosity. Source-local bulk is installed before collision and rebuilt after
native interface conversion, including every newly wet cell. Original
collision, mass transport, pressure reconstruction and conversion remain the
only evolution owners; ordinary and prescribed-gas-pressure steps share this
path. Molecular shear power, molecular bulk power and modeled deviatoric power
are reported separately in SI without thermal deposition. One native candidate
publishes only after endpoint source, speed, viscosity-ratio and cancellation
gates. The successful scratch cache exposes endpoint total tau; source samples
retain molecular tau. This explicit wrapper prevents a later unbound step from
consuming a source query fitted under modeled relaxation. Integration regression
definitions in `tests/material_free_surface_rheology/sgs.rs` cover nonlinear
all-closure forcing, exact disabled-source trajectory and newly wet bulk
coverage. SGS remains exploratory, without turbulence qualification or a
subgrid free-interface model.

### Moving walls with the bound source-rheology and SGS owner

`MaterialFreeSurfaceRheology3::with_sgs_and_planar_wall_velocity_field_si`
binds the original joint molecular/eddy constitutive root with the actual SI
wall field. WALE/Vreman use the prescribed halfway-wall gradient, including
supported two-layer liquid films; Smagorinsky retains its kinetic-stress root.
Binding does not advance populations, time or mass, and does not implicitly
retain wall actuation for a future step.

`MaterialFreeSurfaceRheologySgs3::step_with_planar_wall_velocity_field_si`
recloses the input source/SGS response under the current wall field before
collision. The same fresh source response supplies independent bulk viscosity
inside the original single native wall candidate. This prevents a changed
wall gradient from leaving the previously fitted molecular or bulk response
in the applied collision. The accepted endpoint reuses its actual pressure,
Guo-corrected stress and the held wall field for source admission. Source
samples retain molecular tau while scratch retains endpoint total tau.

Native wall reflection, gas reconstruction, interface conversion, mass
transport and SI wall impulse/work remain authoritative. The input state's
molecular shear, independent bulk and modeled powers remain separate and do
not deposit heat. Source or convergence failure, invalid wall/SI data and
cancellation preserve accepted flow and source samples. Zero SGS coefficient
uses the original molecular epoch and endpoint root, retaining population,
source and wall-receipt parity. Existing fully flooded, tangent-only fixed-wall
limits and exploratory isothermal LES scope remain.

Four integration regression groups in
`tests/material_free_surface_rheology/sgs/motion.rs` cover analytic rotating
quadratic two-layer binding, changed prior wall/source-bulk epochs, independent
nonlinear endpoint-root and wall-work checks, constant-source native wall-SGS
parity, all-three-closure disabled parity, and late source refusal/retry.
All changed Rust files passed parser/source and whitespace review; these
integration tests were not compiled or executed in the constrained environment.

## Original liquid enthalpy volume demand at internal gas pressure faces

`FlowField::project_density_volume_source_with_internal_pressure_faces`
extends the original physical-density pressure graph with explicitly authored
fluid-to-inactive full-face pressures. Each opening uses only the occupied
liquid half-cell mass; the inactive side supplies no density, velocity unknown
or transported material. Unlisted inactive faces remain walls. Component
compatibility, gauges, Krylov solve and independent volume/work accounting
remain with the existing owner. Sealed components must balance independently.
The companion velocity entry supplies an unprojected physical predictor.
Relative source and opening pressure work retain their common component datum.

`MaterialFreeSurface3::project_enthalpy_volume_demand` reads the original
located native inventory as fixed physical M, advances its original caloric
chart by the supplied NET extensive enthalpy increment, and computes
`S = M (1/rho_new - 1/rho_old) / dt`. Current M/voxel-volume supplies projection
inertia; numerical lattice pressure-density excursions never become a guessed
thermodynamic expansion law. Initial chart volume, represented enthalpy and
source/chart density agreement have separate explicit physical accuracy
budgets and measured defects. Original material/claim identities and external
query conditions remain required. Identical temperatures share source snapshots.

This entry admits full occupied Cartesian voxels and empty ghosts. Partial or
overfilled interfaces, conversion carry, piston-remapped geometry and nonliquid
states refuse. Gas correction pressures and signed physical fluxes are exposed
with the original pressure-work report. The result proposes compatible fluxes;
it does not advance native populations, VOF geometry, mass, enthalpy or time.
A conservative moving-interface owner must still consume these fluxes with
material/energy transport. Mechanical correction work is not silently heat.
Ordinary conjugate heat transport refuses nonzero internal solid-face flux,
including a field accepted by this explicitly opened pressure graph.

Nine regression definitions cover graded unequal-density inertia on every
axis, signed openings and large pressure offsets, component compatibility,
original-source heating/cooling and interior partial increments, enthalpy datum
invariance, explicit source-density roundoff budgets, incompatible geometry,
ordinary thermal refusal and cancellation. Changed Rust sources passed parser
and whitespace review; these tests were not compiled or executed in the
constrained environment.

### Localized linear and angular impulses on dense native lattices

`Grid3::apply_wrench_impulse` and `ForcedGrid3::apply_wrench_impulse`
deposit one authored total `WrenchImpulse3` onto an increasing, unique stencil
of actual fluid cells. `origin` fixes the angular-impulse reference in lattice
lengths; native positions are cell centers. The stencil's original population
masses determine its center and inertia. The rotation solves
`I * delta_omega = angular - (center-origin) cross linear`, using the original
rank-aware angular correction. Together with uniform translation, this is the
minimum mass-weighted squared velocity change that supplies the wrench.

Collinear stencils carry transverse angular impulse but cannot invent axial
spin. A single cell carries its point-force moment but no independent couple.
Unresolved angular impulse outside the explicit budget is refused. The API
requires nonperiodic topology; it does not infer an unwrapped moment chart.

The existing Guo-frame impulse retains physical nonequilibrium stress and
density, supplies signed kinetic work, and preserves the original uniform or
cellwise force epoch. The dependent receipt observes original and staged
physical momentum at each selected position. Deposition error and measured
local/global linear and angular defects must pass their explicit bounds before
the original population swap. Rank failure, excess speed, a late torque defect,
or cancellation preserves all accepted fields and completed-step counts.

There is one compact stencil vector plus the native bounded 3-by-3 factorization
workspace; the original next-population storage is reused without a dense
impulse field. No time, collision, pressure projection, heat or receiving-body
reaction is added. The coupled body owner retains its reciprocal dynamics and
work account. Tests cover unequal-mass couples on D3Q19 and D3Q27, off-center
loads and origin covariance, minimum-change deposition, rank deficiency,
spatial-force torque, retained stress, Galilean work, and late refusal/cancellation.

## Native planar liquid expansion from original thermal density

`MaterialFreeSurfaceEnthalpy3::advance_uniform_slab_volume` advances uniformly
heated or cooled planar liquid into actual native full/partial cells. The
liquid occupies the complete transverse section from either anchored domain
wall, with one optional common partial terminal layer; the whole slab may be
thinner than one cell. An actual unoccupied native gas-side layer must remain
at both endpoints. A partial last domain cell beside the opposite wall is
not silently assigned a different pressure boundary.

Initial states require canonical density-one rest populations after the
same density-preserving normalization as native ingress. Moving continuation
requires the exact retained per-layer native population/mass profile produced
by this bridge, unchanged uniform enthalpy and the same anchored face.
Storage is proportional to slab depth, not the entire three-dimensional grid.
Accepted ordinary native ticks invalidate the binding; exact readback also
refuses unmatched changes through other paths. Arbitrary incoming motion or
nonequilibrium, holes, disconnected geometry, internal solids, piston history
and unsupported coupled force or surface owners remain outside this rung.

The existing pressure/source owner receives an explicitly graded metric:
terminal wet width is `fill*dx`, with its missing width added to the adjacent
gas-side cell. Original indices, total extent and other occupied cells remain.
The private metric admission independently checks native
`mass/population_density*dx^3` and original `M/rho(h)` against actual wet
volume. The public voxel proposal retains its original full-cell admission.
The retained affine predictor reaches the actual gas face, including the
occupied half-cell inertia; terminal-center velocity is not copied there.
Pressure impulse is measured relative to actual incoming native momentum.

Exact affine overlaps carry physical M, uniform h and momentum through the
moved interface. Contraction obtains no mass or enthalpy from gas. Endpoint
source density consistently sets physical mass/pressure units; population
density remains the numerical pressure coordinate, not an imposed thermal
EOS ratio. The existing strain initializer supplies the new known affine
molecular response before native mass/momentum readback. The terminal native
velocity stores the liquid-volume mean; its centroid offset from the native
Cartesian sample point is explicitly reported and bounded in metres.

The original caloric chart, datum, source claims, external conditions and
clock remain authoritative. A physical source/chart/reference density
agreement allowance is retained and reported. Ordinary constructors retain
their exact admission; subsequent constant-reference stepping refuses larger
unmodeled density changes. Input is NET caloric increment. Required pressure
mechanical work is separate from actual native kinetic storage and heat.
Incoming native-to-face reconstruction and both outgoing kinetic mapping
changes consume one noncancelling absolute accuracy budget. The pressure
owner's compensated kinetic increment closes the mechanical account. No
projection, reconstruction or remap loss becomes heat. Extensive enthalpy
uses retained Dd mass sums and h1-h0, including the explicit mass-defect datum.

Repeated calls compose this thermal source split. They do not perform native
viscous collision, wall relaxation or mass-route advection and therefore do
not replace a general time integrator. Ordinary native continuation remains
available; general arbitrary-moving thermal VOF and phase change still need
their own physical evolution and conservative transport owners.

Numerical-leaf regressions cover exact affine mass/momentum/kinetic integrals,
noninteger initial depth, subcell slabs, repeated density maps, tiny changes,
extreme finite products, invalid inputs and cancellation. Native integration
definitions cover all oriented axes, partial gas-face geometry, heat/cool/heat
sequences, actual incoming pressure impulse and work, zero-source projection
without new expansion/heat, centroid and kinetic refusals, missing gas-side
capacity, ordinary-tick invalidation, original datum and atomic cancellation.


## Localized wrench on retained free-surface liquid

`FreeSurface3::apply_wrench_impulse` deposits a total linear/angular impulse
on a strictly increasing, unique x-fastest stencil of original located wet
cells. It reuses the native rank-aware `WrenchPlan3` at actual cell centers.
Fluid inertia is population density; interface inertia is its retained liquid
mass. Empty cell volume and unlocated conversion carry supply no inertia or
lever arm. Gas, inactive/solid cells, piston-occupied cells and inadmissible
interface inventories refuse before publication. Collinear stencils retain
their represented transverse couples and refuse unresolved axial spin.

The minimum mass-weighted squared velocity increment passes through the same
Guo-frame population impulse as the original density-field entry. Canonical
density, physical nonequilibrium stress, fill, tracked mass, classification,
force history, geometry epoch and time remain fixed. Independent angular
readback integrates the actual before/after physical velocities with the SAME
retained liquid masses. Local/global torque, linear momentum, signed kinetic
work and low-Mach limits are admitted together before any population publishes.
A dependent receipt and final cancellation poll precede publication as well.

Storage and work follow the selected stencil: two compact vectors, native
sparse tile lookups and the bounded 3-by-3 rank solve. There is no full-domain
scan, dense impulse array, liquid-owner clone, collision, pressure projection,
heat assignment or implicit receiving-body reaction. The coupling owner must
supply reciprocal body dynamics and its own work account.

`MaterialFreeSurface3`, `MaterialFreeSurfaceRheology3` and
`MaterialFreeSurfaceEnthalpy3` expose the same entry with
`MaterialFreeSurfaceWrenchImpulse3` in kg m/s and kg m²/s. Original cell size,
density and velocity units convert the total wrench; no extra timestep is
applied. Linear/angular SI defects are re-evaluated against the authored
physical request inside the native transaction, including origin-conversion
transport and representability. Rheology and thermal histories remain attached
to the original source state; signed mechanical work changes resolved kinetic
energy and is not deposited again as enthalpy.

Native regressions check unequal partial-fill inertia, actual torque/work,
retained Guo stress, untouched nonselected cells, origin/boost covariance,
rank refusal, invalid inventory, speed limits and dependent-receipt/cancellation
atomicity. SI receipt regressions exercise distinct mass/momentum/angular/work
units, authored-origin transport and unrepresentable readback. The material
integration regression compares the same physical wrench at two timesteps.

## Original-reference density through deformed tetrahedral liquid ingress

`free_surface3_from_reference_tetrahedra` and
`free_surface3_from_reference_tetrahedra_with_enthalpy` retain original
reference positions and per-tetrahedron P1 corner density on the current
snapshot's original connectivity. Positive affine deformation transports
density as rho_reference/J. Original mass remains the reference integral;
first and second mass moments use its affine push-forward. Current
tetrahedron/cell intersections use the existing exact density-moment owner.
Original P1 velocity and enthalpy retain their cubic kinetic and quadratic
thermal mass integrals.

These entries use the existing component momentum correction, native
initializer, thermal attachment and projection-error accounts. The receiver
retains its original uniform physical density and supplied geometric fill.
Explicit kg/m3 density and m3 occupied-volume bounds admit the measured
reference/current compatibility. Original overlap mass partition and actual
per-cell native mass defects share the existing absolute mass bound.
Numerical pressure cannot replace physical density, and native fill is not
rescaled to conceal an incompatible source inventory.

The report retains reference and deformed density moments, maximum measured
current-density discrepancy, mass partition/readback defects, and actual
native occupied volume from fill_fraction*dx^3. Thermal ingress repeats this
audit after final native stress/constitutive attachment. Signed kinetic and
gravitational projection changes remain numerical accounts and never become
heat. The input must already be liquid with an independently admitted
disposition of stored solid stress/history; this does not supply melting,
general variable-density VOF or a phase-density transition.

The production kinetic regressions cover the independent graded-density
cubic form, exact uniform-density compatibility, and retained variance beside
large common transport, permutations and mass scaling. Native integration
regressions specify deformed reference inventory, moment/enthalpy
conservation, refusal and continuing flow behavior. MR28 remains open for
its broader phase and history obligations.

## Supplied heat drives conservative liquid slab expansion

`MaterialFreeSurfaceEnthalpy3::prepare_uniform_slab_heat` closes an authored
heat interval against original caloric storage and the existing pressure
owner's mechanical correction work: Q = M*delta_h + W. Fixed-pressure enthalpy
already contains ambient p*dV; correction pressure does not change the
thermodynamic source context. Projection loss and all native/face/affine
mapping changes retain their separate numerical accounts and never become heat.

On the admitted original slab metric, affine normal motion gives
W = B*u*(u-u0). B is measured once by the original pressure owner's energy
audit on an already-compatible unit-speed field. Scalar trials use the SAME
original chart and analytic specific-volume tangent, with a caller-selected
finite sign bracket, evaluation limit and joule residual. No repeated 3-D
pressure solve or material-property refresh occurs inside the scalar solve.
The owned original-source caloric endpoint keeps the accepted density and
expansion law available for the continuing native liquid.

`FreeSurfaceSlabHeatPlan::advance` delegates to the existing conservative
native slab step and independently gates its measured pressure work and
represented enthalpy against supplied heat. The plan borrows the original
trajectory; failed preparation or advancement changes no accepted state.
All six orientations, a partial terminal layer and subcell slabs retain the
existing mass, momentum, geometry, source and cancellation gates. Repeated
calls remain thermal source splitting; arbitrary flowing geometries, viscous
integration, capillarity and melting require their existing physical owners.


### Conservative dyadic resolution transfer of sparse liquid momentum

`FreeSurface3::refine_dyadic` and frozen
`MaterialFreeSurface3::refine_dyadic` return a new original owner on the same
bounded physical domain. Every spatial axis is refined by two and every active
tile supplies eight child tiles, with dx and dt both halved. This P0 transfer
gives each symmetric child the parent's actual intrinsic density and retained
liquid fill; physical child mass is one eighth. Explicit budgets admit actual
stored mass, first mass moment, linear/angular momentum, kinetic storage,
speed, six-mode physical viscous stress, total domain cells and active tiles.
Native receipts use original coarse units; the material receipt uses kg, m, s,
J and Pa and independently checks rounded SI totals before returning the candidate.

Uniform BGK viscosity is retained through tau_f=2*tau_c-1/2. Physical
Guo-frame nonequilibrium is scaled by tau_f/(2*tau_c), preserving physical
viscous stress with its (tau-1/2)/tau factor. Native velocity/pressure units stay
fixed, held force density halves and native surface tension doubles, subject
to the original fine capillary admission. Uniform/affine acceleration remains
the original force owner and is evaluated at actual child centers with a
zero-time force reframe. Fixed solid tiles, uniform planar contact angle and
integer piston location retain their physical geometry. Geometry epoch
increments because cell identities change; source state stays immutable and
no elapsed time, work, heat or conversion event is authored.

This owner refuses unlocated carry, non-BGK collision, local surface-tension/
bulk-relaxation fields, spatial or dynamic contact-angle history and
unrepresentable scales/populations/balances. The SI route retains the original
frozen source and pressure datum; evolving thermal/rheology and other histories
need their own transfer owners. It reconstructs no hidden subcell free-surface
shape or higher velocity modes, performs no local AMR/coarsening, and does not
promise identical later trajectories on the finer discretization.

Eight focused checks passed in this session: five tests compile the complete
new native refinement module unchanged with actual sparse/native storage and
original constructor, read/write, sample, force and closure source spans,
compensated sums, cancellation, scalar BGK collision and Guo frame/stress
routines. Three tests execute the unchanged SI conversion/receipt helpers with
a minimal error facade. The focused harness does not include full streaming,
conversion timesteps or material-owner construction. Two original material
integration regressions in `tests/material_free_surface3/refinement.rs` specify
coarse evolution, independent physical inventories, source/pressure/unit
retention, fine evolution and atomic refusal. Those complete integration tests
were added but not executed under the current full-Cargo resource limit.
