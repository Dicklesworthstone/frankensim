# CONTRACT: fs-couple

## Adaptive time refinement of sourced directional friction and heat

`friction_heat::directional::trajectory::march::advance_source_directional_heat`
advances a held-geometry span using the original coupled phase. Every comparison
starts one full step and two half steps from the same accepted checkpoint. It
compares original linear/angular momenta, each finite cell's energy, unsigned
path, net slip and physical friction work in explicit physical error scales.
Relative scales use increments from that checkpoint, preserving independence
from background momentum and thermal reference energy. Only the actual two-half
result is accepted; no extrapolated momentum, heat or material state is created.
Time refinement updates the original source coefficients and thermal endpoints
as bodies slow down, stop and continue exchanging heat.

An explicit mean-speed domain is intersected with the supplied phase-distance
domain using each represented duration. Smaller timesteps cannot accidentally
turn a fixed trial distance into an out-of-domain source speed, and no bound
based on initial Euclidean speed is assumed. Invalid/nonfinite domains refuse.
Source/path/thermal/geometry refusals propagate rather than being reclassified
as successful refinement. Native force, quadrature, source and Newton limits
remain active; separate comparison and accepted-phase budgets bound total work
and retained accepted payloads. Cell comparisons poll cancellation in bounded
chunks. Source callbacks evaluate private candidates and must not mutate
external physical state.

The complete requested span is atomic, including a failure after private
successful steps. Only the chosen two bodies, original source history and
original thermal timeline publish at the end. Compensated physical-transfer
totals and direct origin-to-end kinetic/cell-energy changes check whole-span
mechanical, thermal and combined closure against an absolute energy tolerance.
Every owned accumulated-path coordinate is independently compared with total
curved travel, limiting accumulated source-root and rounding defects. Accepted
phase reports retain their exact start/end times and source transport payloads.

This is local step-doubling error estimation for a fixed-pose, fixed-load,
endpoint thermal/source discretization. It does not certify global error,
source-root uniqueness, spatial thermal convergence or discontinuous-law
convergence; no external drive, normal dynamics or thermal event scheduler is
added. G1/G3/G4 tests include three time rungs against an independent continuous
temperature-dependent friction ODE, adaptive arrest followed by conduction,
source-speed domain admission and whole-span rollback/retry after late source,
cancellation and work-budget failures.

## Sourced directional arrest coupled to finite heat

`friction_heat::directional::trajectory::apply_source_directional_trajectory_heat`
joins the original source-owned continuous directional relaxation, full-tensor
rigid-body impulse owner and existing `FrictionHeatStep` residual. The native
trajectory provides integrated impulse, curved net displacement, unsigned path,
finite-arrest time and physical traction-slip work as separate quantities.
Average force times net displacement is not substituted for curved work.

Each thermal trial starts from the same accepted source and two body states.
Original pinned principal claims are queried at the solved finite surface
temperatures, optional common interface temperature, load/actual-area pressure,
phase-average unsigned speed and each accepted source path plus actual travel.
The same original caloric objects and located contact face are required. The
common zero-capacity interface derives signed heat into the existing finite
cells from their geometry/conductivity; no heat fraction or parallel store is
created. Physical friction enters once, while hot-to-cold exchange continues
after sliding has stopped.

Exactly the selected two bodies, source history and original thermal timeline
publish together after native impulse/momentum/work checks, finite-cell/interface
residuals, independent Delta(K+H)-Q_external admission and cancellation. Source,
domain, numerical and budget refusals leave all accepted owners unchanged.
Transport payloads remain attached to the accepted candidate.

This is a fixed-geometry, held-normal-load phase with constant sampled coefficients
inside its native curve. Endpoint thermal storage and phase-average coefficient
sampling require time refinement, including the timing of heat released before
arrest. Source regularity and physical branch selection remain necessary;
nonlinear convergence does not establish unique source roots or spatial thermal
accuracy. G1/G3/G4 tests cover independent finite-cell balances through arrest and
rest, temperature-dependent stopping, full-tensor/frame-covariant curved travel,
original query receipts and atomic refusal/retry.

## Spatial shock loads on original structural pressure surfaces

The opt-in `experimental-shock-loading` feature supplies `shock_loading`.
Frozen and species-bearing one-dimensional Euler adapters consume the accepted
SSP RK wall momentum exchange, retaining its actual accepted duration. They
require reflecting, stationary, adiabatic ends with zero mass/energy transport.
An explicit global axis fixes orientation; held exterior absolute pressure is
subtracted once from the structural load, never from the gas flux. The caller
supplies the most recent report; a mismatch with current extensive state refuses.
The optional `experimental-compressible` feature adds the original 2-D and 3-D
Cartesian Euler owners' latest per-face wall impulses, with orthonormal frame
embedding. `WallPressureImpulse::from_euler3d` retains all six boundary sides
and both tangential face coordinates. A nonuniform volumetric pressure footprint
therefore retains both transverse moments instead of becoming an averaged load.

Cold surface preparation reuses `fs-plate::ShellPressurePort` and retains sparse
original-node P1 pressure rows. Each explicitly selected planar patch belongs
to one fluid face. Duplicate facets, area/normal mismatch, differing accepted
intervals, work-budget failure and cancellation refuse without a partial load.
Unloaded and supported nodes keep their original identities. Returned impulse,
interval-mean force, resultant, torque and spatial virtual work preserve this
discrete pressure map. Per-face Cartesian loads retain nonuniform shock moments.

This is a stationary-wall load transfer and a one-way prescribed structural
forcing path. It does not move the gas boundary, assert joint fluid/solid energy
closure, resolve within-step load variation, model fracture or validate material
shock response. G1 source-Hugoniot reflection and independent fluid/support
momentum, G0/G3 affine surface work/moment and native modal impulse-consumer
tests accompany explicit boundary/geometry/refusal and cancellation checks.
The 3-D tests independently reconstruct all-six-wall fluid/support momentum,
rotated per-face resultants and both moments, and affine P1 virtual work with
exterior pressure. They also check retained loads after cancelled source steps
and refusal of periodic/open faces, bad embeddings and invalid face indices.

`WallTractionImpulse::from_molecular_euler2d` and `from_molecular_euler3d`
consume the source-owned molecular Euler operators' complete accepted loads,
including both molecular halves and the middle advective step. They never
substitute the borrowed Euler owner's partial report. The gas-side normal
traction is signed `p-tau_nn`; tensile viscous loads remain negative rather
than being clipped or represented as absolute pressure. `WallPressureImpulse`
retains its positive Euler-pressure meaning. Both types share one private
normal-load account and the same original P1 projection kernel, exposed as
`project` for pressure and `project_tractions` for signed traction. No temporary
pressure-vector allocation or extra force distribution is introduced.

The molecular adapters still require stationary adiabatic slip walls with zero
mass, energy and tangential momentum transfer. Explicit exterior pressure is
subtracted over the same accepted interval. Actual transport-law tests exercise
a tensile 2-D wall load through nodal force and affine virtual work, and a full
3-D six-wall momentum balance with retained loads after late cancellation.
Both cases differ from their Euler-only middle-stage impulses. They transfer
normal traction only; moving walls, tangential no-slip stress and reciprocal
fluid/solid energy closure remain separate capabilities.

## Constitutive feedback through public material equations

`constitutive_graph::MaterialGraphResidual` closes caller-authored, dimensioned
balances of original law inputs and outputs on the existing Newton--Krylov
solver. Scaled unknowns and residuals keep physical units explicit. Jacobian
actions use the material graph's exact chain rule at the same committed
history as each residual; graph-internal edges remain acyclic. A finite work
budget, latched original source errors, cancellation and independent final
equation tolerances precede returning a complete candidate. No physical
history mutates during a trial or failed solve. The native registered
hyperelastic material-card consumer and a closed two-variable electrothermal
reference exercise this API. Convergence of arbitrary caller balances does
not confer conservation, passivity, calibration or a time-discretization law.

## Source-integrated conductivity in composed heat transfer

`PreparedHeatNetwork::with_integrated_conduction` routes the existing material
owner's integrated k(T) layer fluxes through the same circuit, thermal and
observer callbacks. Each selected link retains its original conductivity
property/claim, material card, caloric object, fixed conditions and distances.
The native conductance map and reported link secants both contain the applied
flux. Internal heat still transfers once between the original cells; no extra
source or time integrator is introduced. Immutable correspondence is admitted
once, while the original layer owner checks actual phase, support and geometry
on every trial. A separate thermal-contact owner on the same link refuses,
and topology changes require explicit rebinding. The admitted layer model is
quasi-steady, homogeneous, isotropic and pure-phase; it adds no within-layer
storage, tensor transport or moving phase boundary. Focused synthetic tests
compare the nonlinear flux with an independent linear-k integral and verify
that data-only conductivity changes reach finite heat and electrical history.

M4.1-W1 extends `unilateral_contact` with explicit physical span measures [m],
per-length stiffness, static/kinetic coefficients and a declared friction ramp
speed (mandatory for nonzero kinetic friction). No 1/n normalization or hidden
1e-3 m/s scale remains in the adapter; constitutive refusals retain owner/field
in a typed DContactError. Loss-port Results propagate. `binding::bound_point_obstacle`
lowers SourcePointContact verbatim and is also used by MaterialContact.
`BoundNormalStorage` composes the existing point or finite-gap potential with
the pHS Storage contract, preserving original U/F/dF and source errors in its
fallible API; the infallible Storage view signals invalid trials as nonfinite
values for the numerical owner's refusal. Whole finite-gap footprints already
integrate their areas and are not multiplied by another area weight. Sources
must have an unloaded zero-approach origin. Point Hunt--Crossley loss is a
separate port evaluated at actual discrete effort and explicit dt; its envelope
is checked before projection. No time integration, cancellation framework,
material data, thermal evolution or interpolant authority is invented.

Multiphysics composition through a versioned port-thermodynamic vocabulary: a
lossless interface relation, explicit storage/dissipation/source primitives,
a complete stream-flux bundle, and evidence-bound closed-window physical
accounting without false passivity.

## Reacting atmosphere on the native moving-liquid timeline

`GasReaction::liquid_candidate` in `liquid_gas` composes the original closed gas
reaction owner with `fs-lbm`'s finite-atmosphere/VOF enthalpy owner. The native
liquid timestep supplies one interval to chemistry, thermal exchange and flow.
The reaction changes composition at conserved total internal energy, retaining
species formation references. Its accepted original gas sources then supply
film conductivity, caloric inversion and pressure during heat exchange. The
same liquid source laws feed native viscosity, buoyancy, conductivity and any
attached current/heat response. Subsequent calls use the actual cooled gas and
remaining species inventory.

The candidate clones native populations, carry, source/contact histories and
finite walls; no caller-owned state changes on source, thermal, reaction,
combined-energy or cancellation refusal. An independent final account measures
gas internal energy, post-advection liquid enthalpy and finite wall enthalpy
against authored/electrical/physical viscous heat. Interfacial heat cancels once;
formation energy is not credited again and numerical kinetic losses are not
thermal sources. Cancellation is checked at existing source and solver stages
and immediately before candidate publication.

This is a first-order reaction/thermal/flow split, with well-mixed fixed-volume
gas, the original constant-density liquid approximation and explicitly exposed
film geometry. It adds no gas traction, flame or shock propagation, evaporation,
phase conversion, moving wall or global mechanical/entropy theorem. Time and
space refinement remain necessary. `tests/reactive_liquid_gas.rs` supplies G1
independent isomer/finite-capacity balances, G3 data-only changes that reach
native VOF populations, continuation and G4 rollback checks; synthetic fixtures
do not constitute empirical combustion or liquid validation.

`GasReaction::liquid_surface_candidate` shares that chemical staging and
whole-interval audit with the native solved liquid/gas/radiation interface.
Product composition, NASA-9 energy and transport feed its surface-temperature
solve; the original liquid `sigma(T_surface)` then drives native capillary
pressure and Marangoni traction while bulk transport retains its enthalpy
temperature. Radiation enters the external account once, and gas/liquid
transfer cancels within the combined finite inventory. The complete-surface
geometry and source restrictions, fixed gas mass/volume and first-order split
remain binding; this adds no surface-energy or mechanical-energy theorem.
`tests/reactive_liquid_gas/surface.rs` supplies independent zero-storage
radiative and finite-capacity references, data-only kinetics and surface-tension
changes reaching native populations, continuation and refusal rollback.
These are manufactured model checks, not experimental validation.

## Assembled shell adhesive dissipation

`shear_shell::assembly::MaterialShellAssembly::reduce` projects the original
mortar glue damping through the same assembled mass-normalized modes as
stiffness, component viscosity, bridge forces and pressure surfaces. It sums
the complete modal tensor before the existing positive-semidefinite damping
owner realizes passive channels, preserving cross-mode terms and undamped
subspaces. With fixed stiffness and geometry, changing glue viscosity changes
neither the elastic basis nor the radiation map. Zero supplied viscosity adds
no channel; no loss is inferred from a material name.

`has_viscous_glue` reports supplied positive edge or face viscosity. Projection
and connection budgets include this work. Only rigid mortar ties eliminate
slave coordinates; supports on compliant slaves remain active. Original tie
owners still validate coefficients and geometry. Source adapters own the
validity of a constant-coefficient glue law across the retained state and
frequency range.

## Purpose and layer

Layer L3 (multiphysics coupling). Domain assembly depends on physical
solvers. The neutral vocabulary is owned by `fs-ports` and reexported here
without conversion or duplicate types. That lower boundary depends only on
`fs-qty` and `fs-iface`; FEEC imports it without depending on domain assembly.

## Public types and semantics

### Shared cavity and distributed downstream duct

`SharedCavity::with_downstream_duct` binds the unchanged compact cavity/neck
to an additional source-bound `DgPortLoad`. Admission requires an
`ExteriorMouth` neck with no direct exterior load, equal inlet gas and area
(roundoff-only area allowance), an ordinary non-driven inlet, cold accepted
clocks and identical step duration. The original neck retains its bore, inner
end inertia and selected viscous law. DG owns additional distributed length,
wall memory, branches and outlet terminals.

A scalar displacement mobility is prepared once. Each step solves mean inlet
pressure against current complete DG history; the original modal and DG
candidates publish together only after finite pressure, swept-volume,
equal/opposite port work and joint budget checks. The original energy ceiling
covers both owners' represented storage; the pressure ceiling also bounds mean
DG inlet pressure. Rational terminal work remains signed boundary exchange.
Dense-DG loss inferred from its storage deficit is explicitly reported as such,
not as an independently integrated heat channel.

The held-pressure/held-mean-flow connection approximates the continuous
interconnection; it is not the joint continuous exponential. Coupling-clock,
spatial and constitutive-band convergence remain separate. Tests replay both
original propagators with accepted interface quantities, retain actual viscous
and thermal memories, refine the coupling clock at fixed geometry, and refuse
late candidates without publishing either owner. The returned cavity accessor's
harmonic response describes the cavity subsystem and excludes the attached DG
load. A broadband 5 cm/order-2 fixture had nonmonotone 24-to-48 kHz
endpoint error. Its retained DG spectrum reached 47,649 rad/s; the unchanged
error gates passed at 96/192/384 kHz against a 1.536 MHz fixed-spatial reference.
This is fixture convergence evidence, not a universal timestep certificate or
a continuum validation.

### Continuous loaded modal and cavity response

`CoupledModalSystem::harmonic_response` solves the original continuous linear
network in peak phasors with `exp(-i omega t)`. Original component damping,
complete signed connection columns, Maxwell material relaxation and attached
multiport radiation all enter the same complex system. Maxwell histories are
eliminated exactly; acoustic oscillators remain explicit, preserving a finite
joint antiresonance at a bare lossless acoustic pole. The existing cancellable
complex LU supplies the solve, with an original-equation backward-residual
gate. Returned original modal and memory coordinates include cycle-mean storage,
source power, each physical loss and their power residual. These are the AC
particular solution about equilibrium, independent of accepted transients and
DC rest lengths; analysis never advances the clock or changes a ledger.
Prepared constant-loss connections reuse the original `g(K)b` damping columns,
including their full cross-mode reaction, while springs keep their original
columns. Returned separate spring/damping efforts and the summed generalized
modal reaction retain those distinct directions. The legacy scalar effort sum
is a channel diagnostic when the columns differ. This is the continuous response
of the admitted prepared image; no new loss law is evaluated at the drive
frequency. G1 tests compare resonant amplitudes with the anchored one-mode law,
and a three-clock ladder converges to the same response through the original
prepared time owner with a nonproportional two-mode connection.

`CoupledModalSystem::harmonic_responses` shares one dynamic-stiffness assembly
and factorization across 1..64 peak-force right-hand sides at a common forcing
frequency. It reuses the original cancellable complex LU solves, original-
equation residual gate, observer and physical accounting for every column.
For d mechanical-plus-air coordinates, G constant-loss links, C connections,
M Maxwell branches and r inputs, the admitted work is d^2*[d+G+r*(C+M+3)];
r=1 retains the scalar allowance. Retained response arrays and struct headers
are bounded at 8 MiB with checked arithmetic before allocation. No coordinate
is truncated and no accepted state changes on success or refusal. Both friction
owners expose the same released batch path. Four production batch tests pass
natively for scalar parity, coherent cross terms, reciprocity, phase cancellation,
work/storage bounds and unchanged state. The complete harmonic module reports
13/14 passing; its unchanged two-body test fails during initial time stepping
before any harmonic call. The same heat-resistance setup failure was reproduced
with pre-batch source c72b13d1881ace6718aefc3c4937ef2a5358b25f.

`ModalPreparation::HarmonicOnly(PreparedModalBudget)` constructs these same
continuous equations without a coupled time propagator, quadratic heat-integral
matrices or held-step friction compliance. `Ordinary` and `Prepared` preserve
the existing playback choices; conversion from `Option<PreparedModalBudget>`
preserves existing callers. The shared preparation owner retains every original
component, connection, material history and radiation pole, and prepares the
same constant-loss `g(K)b` directions using the original stiffness/eigen owners.
Coefficient, material, radiation and internal-component rebuilds preserve the
chosen purpose. Time staging, prediction, compliance and joint-DG time
preparation refuse for harmonic-only images before accepted state changes.
Explicit mode, byte and work allowances still apply; the construction allowance
includes conservative harmonic workspace and coefficient-preparation estimates.
Original physical component initialization and shell reduction still run.

Connection and material rank forms assemble on their exact nonzero support,
using one reused ordered scratch vector. Only exact zero coefficients are
skipped; there is no numerical sparsity threshold or physical truncation. All
retained coordinates, the complex LU, original-equation residual and conservative
work bound are unchanged. Two native support checks compare independent dense
outer products; two native construction checks establish bit-for-bit response
and constant-loss parity with prepared playback, including coefficient updates,
and state-preserving refusal/cancellation. Full-bank and friction construction
regressions remain authored but unrun.

`SharedCavity::harmonic_response` maps original surface forces and incident
neck pressure into that same solve and returns physical internal pressure,
surface volume rate, outward neck flow, mutual exterior reaction and gas
compression availability. Internal and radiated pressures remain reactions,
not extra external drives. Thermal pressure uses the admitted compression law;
viscous and thermal queries refuse outside their sampled admission bands.
Radiation validity remains with its producer, including the declared band of
geometry-derived fits.

A caller work allowance bounds `d^3+d^2*(connections+constant_loss_links+material_arms+3)` before
dense allocation, with `d=mechanical_modes+acoustic_poles <= 1024`. No mode or
pole is truncated. Nonpositive frequency, singular or unrepresentable response,
and missing prepared constant-loss columns explicitly refuse.
This continuous response is not the sampled held-force transfer. Independent
physical pressure/flow condensation, pressure-diffusion availability, separate
loss channels, signed basis invariance and unchanged accepted histories are
tested alongside analytic loaded resonances.

### Continuous finite side-branch resonance

`CoupledModalSystem::harmonic_response_with_load` and the shared-cavity wrapper
retain the main DG bore, every original finite `DgBranchLoad`, and the complete
original MIMO terminal in one continuous frequency solve. Each junction adds an
explicit flow unknown with equal original pressure traces and opposite inflows.
No isolated duct impedance is inverted, including at a bare branch pole.
Terminal channel IDs retain their authored order; a branch channel observes its
actual distal outlet. Original wall memories, material bands, mutual terminal
entries and strictly proper terminal states are preserved. The complete state
and junction count consumes the existing dense work allowance.

Returned branch responses carry original fields, storage, wall losses and both
junction powers. Whole-network accounting sums every owner once. Independent
original pressure-trace checks accompany the complete-equation residual; queries
never advance accepted histories. Analytic cylinder-network peak and linewidth
refinement, a bare branch pole, mutual/strictly proper terminals, branch and
channel permutations, viscothermal bands and shared-cavity accounting are covered
by the native focused tests. These are zero-volume point junctions joining finite
plane-wave ducts, not a resolved three-dimensional entrance-flow model. Finite
gap channels and driven exterior mouths remain explicit refusals.

### Continuous finite side-branch evolution

The explicit exact downstream backend also retains every original branch state
on its common clock. In original DG energy coordinates, junction columns are
`H=[-B_bore,+B_branch]`; the complete normalized Gram solve constructs
`P=I-H(H^T H)^-1 H^T`. The generator `P L P` and projected external columns
preserve the original pressure constraints from the admitted cold state. This
eliminates only zero-storage junction flows; it does not reduce the acoustic or
wall-state basis. The same cached exponential and quadratic-integral owner
evolves the coupled modal, cavity, duct and original rational terminal states.
Terminals admit conditioned invertible D or exactly zero D. For zero D,
original states obey `x'=A x+F U+b u` and `H x=0`; the conditioned solve
`(-H F) U=H A x+H b u` removes only zero-storage flow multipliers. Full
projected inlet coupling, original pressure/work observations and all DG/rational
coordinates are retained. Initial and endpoint constraints and their
differentiated identities must pass; no artificial resistance or state reset
is applied.

Each original duct's unprojected loss form is integrated separately. Junction
flow includes its acoustic-state, modal-velocity and complete terminal-feedback
terms. Both original pressure impulses, both actual integrals of `p(t) Q(t)`,
and swept volume are independently observed. Original per-duct energy closure,
endpoint and integrated pressure continuity, junction work, terminal work and
whole-system budgets must pass before the existing `DgPortTrial` publishes any
history. Terminal work remains signed external exchange because a generic
rational realization supplies no physical storage certificate. Tests compare
24/48/96 kHz constant-drive trajectories, reversed branch ownership, original
wall memory, an independent physical mutual-terminal heat ladder, and late
refusal/cancellation with identical retry. Preparation retains the existing
256-dimensional complete-state allowance; finite gaps, driven mouths,
mixed-rank singular D and unresolved higher-index pressure constraints remain
outside this backend.

### Shared cavity exterior loading

`cavity_modal::shared::SharedCavity::with_exterior_radiation` connects the
existing passive multiport model to the original surface/cavity/neck exact
propagator. A selected neck column is in physical pressure/volume-flow units;
its power-preserving pullback is `Q = S/sqrt(M) q_dot`. Other selected columns
retain the original surface modal velocities, including signed mutual coupling.
The neck must declare `ExteriorMouth`, so exterior mass and resistance are not
counted twice. Additional prescribed exterior pressure remains an incident
source. `neck_exterior_pressure_pa` observes the accepted internal reaction;
feeding that pressure back as an incident source would double-count it.
Analytic coupled-pole motion, exact time partition, energy balance and rejected
step preservation are tested in `cavity_modal::shared::tests`. Geometry/band
validity of the supplied impedance remains with its producer; this composition
does not claim a new BEM fit, distributed neck field or experimental validation.

`SharedCavity::new_with_viscous_neck` replaces the frozen internal neck
resistance with the existing positive `ZkDiffusive` circular-bore law:
Poiseuille drag plus passive `R s/(s+p)` arms, represented by free auxiliary
masses in the same coupled exponential. Only actual bore length contributes
viscous drag and storage; inviscid bore and inner-end inertia remain counted
once, and `ExteriorMouth` radiation can be attached separately. Auxiliary
states consume the original mode, connection and work budgets, retain their
kinetic storage, receive no external port or radiation column, and participate
in accepted-state rollback. `ViscousNeckOptions` declares the sampled Bessel
accuracy band and error budget. This does not claim thermal wall compliance,
distributed neck propagation, nonlinear flow, measured fit accuracy or realtime
throughput. Tests compare independent duct Bessel impedance, a physical-flow
RK4 reference, energy balance and exact timestep partition.

`SharedCavity::with_baffled_exterior` derives the complete exterior reaction
from supplied original signed normal-velocity traces and the actual circular
neck opening in one infinite rigid baffle. The generated annular-sector area
quadrature retains the original mouth radius and area; body traces must retain
every original surface coordinate and reproduce its signed cavity volume
projection. The existing power-consistent Rayleigh owner supplies all mutual
surface/mouth loads, and the existing positive-residue fitter retains its
unchanged training/held-out complex and resistance error gates. Fitting uses a
common mechanical mass-normalized basis before the neck is converted to the
physical pressure/volume-flow convention of the existing attaching method.
Explicit cell, frequency, sweep and complete-pole budgets apply without field
truncation. Preparation is cold; ordinary stepping keeps the existing loaded
propagator and accepted-state accounting. A three-rung mouth quadrature checks
the independent flanged inertance and compact radiation-power limits; signed
mutual loads, energy and exact time partition are also tested. This is a
sampled, spatially discretized infinite-baffle reference, not a curved/unbaffled
body solution, quadrature-error certificate or out-of-band guarantee.

### Cylindrical cavity thermal compression

`SharedCavity::new_with_cylindrical_thermal_wall` uses an explicit physical
cavity radius and the original volume, with an isothermal cylindrical sidewall,
insulated endcaps and the existing uniform-pressure assumption. Positive ZK
thermal compliance is inverted through the existing Jacobi owner into
isothermal equilibrium compression plus positive Maxwell arms, retaining the
adiabatic instantaneous limit. Thermal history, acoustic availability and
diffusion dissipation remain in the original constitutive-memory/exact account.
Pressure and compression-energy observations include accepted thermal branch
reactions once. Relaxed versus instantaneous initial compression is explicit;
viscous-neck and exterior radiation loading compose without changing physical
port count. Three-rung independent Bessel and physical pressure/diffusion
time-domain tests cover this compact reference. No arbitrary-cavity geometry,
finite wall conductivity, distributed acoustics or full wall-temperature
evolution is inferred.

`SharedCavity::new_with_cylindrical_molecular_wall` and
`CylindricalThermalCompression::new_with_molecular_relaxation` bind the
original gas reservoirs to those same radial heat coordinates before lowering
compression. The complete material gas snapshot must match the cavity source.
A positive capacity/conductance pencil yields positive Maxwell arms in one
joint symmetric eigensolve; the existing coupled propagator, harmonic solve,
pressure observation and material-memory work account consume them unchanged.
There is no second heat deposit. Setup explicitly admits the joint coordinate
count (thermal poles times one plus nonzero molecular reservoirs; ceiling 128).
The frozen molecular instantaneous stiffness and equilibrium-gas adiabatic
reference have distinct accessors; the fully relaxed limit stays isothermal.
A 33-frequency check bounds sampled complex and dissipative realization error
to 1e-8 within the original wall band, separately from equilibrium Bessel
admission. Joint spatial error remains a thermal-pole refinement question.
Independent tests cover a quadratic constitutive oracle, three-rung complex
Bessel convergence, reservoir splitting, physical temperature ODE storage/heat,
the original loaded two-coordinate solve, and accepted-time work closure.
This preserves local Debye reservoirs and the original constant Fourier
conductivity; separate molecular-energy diffusion, finite wall conductivity,
nonuniform pressure, flow and changing background state are outside this lane.

### Molecular cavity compression

`SharedCavity::new_with_molecular_relaxation` maps the original material
`RelaxationState::bulk_relaxation` spectrum through the actual cavity volume
and neck area into the existing reciprocal modal solve. Internal pressure,
compression storage, molecular loss, harmonic response and accepted histories
use the same source gas. Original surfaces, neck viscous memory, downstream DG
and exterior radiation retain their owners and work accounts. The supplied
material law must match the complete gas snapshot; missing/extra history,
setup-budget failures and cancellation refuse publication. Independent
internal-temperature dynamics and thermodynamic availability check the
time-domain pressure/energy; direct caloric response checks harmonic pressure
and dissipation in `cavity_modal::shared::molecular::tests`.

The medium is frozen, the compact pressure spatially uniform, and the cavity
wall adiabatic. Molecular loss is acoustic availability dissipation, without a
second gas-heat deposit. For a conducting cylindrical wall, use
`SharedCavity::new_with_cylindrical_molecular_wall` above: its joint caloric law
replaces the adiabatic specialization. The two compression owners cannot be
added as independent shunts. No measured instrument fidelity or evolving gas
chemistry is inferred.

### Accelerating material-fiber supports

`MultiBowPerformance::with_fiber_support_program` supplies declared continuous
piecewise-linear axial velocity to the existing `with_fiber_support_motion`
owner. That owner projects `Gamma_n = integral(rho A phi_n dx)`, supplies
`-Gamma_n a` to the original modal/contact/nonlinear solve, uses step-average
support speed, and migrates the crossing along the original material hair.
There is one inertia law and one accepted support state. The program rejects
competing tangential gesture speeds or another support-motion control. Velocity,
acceleration and complete translation extrema are bounded before playback;
the original finite span independently bounds the moving crossing. Original
materials, pressure, state and sample clock survive; callback partition cannot
change the motion.

The native mechanics stores fiber-relative kinetic energy. `FiberSupportFrame`
adds `V Gamma.qdot + M V^2/2` to recover laboratory mechanical energy. Physical
support work is native bow drive work plus inertial-coordinate work plus the
change of that kinetic correction. Both direct and programmed support motion
use this account. Component-state and work bounds admit representability before
any physical step; publishing an accepted frame is infallible. Inertial work is
never heat. Prescribed normal load and dry friction remain this moving-hair
rung's assumptions; normal/rosin/thermal history transport and stick flexure
are not supplied. Focused regressions compare free acceleration with the
independent sine-wave solution, work closure, width refinement, material
substitution, callback partitions and the original moving-contact owner.

### Coupled extension and normal deflection of material fibers

`TensionedFiberRibbon::with_axial_transverse_strain` uses the original
`fs_nlmodal::longitudinal::LongitudinalString` storage, with each quadrature
strip's original mass, EA and tensile preload. Normal deflection changes
`u_x + y_x²/2`, generating axial waves and changing local tension and contact
reaction. All original fiber coordinates and linear material losses remain in
the shared exact-linear/contact/discrete-gradient owner. Its existing axial
quadratic storage is subtracted once from the nonlinear correction. Disjoint
cold material potentials compose under one identical solve configuration;
overlapping coordinates or incompatible budgets refuse. Existing string, body,
cavity and pressure states are retained. Continuous-span slope/strain bounds
and positive local tension gate every trial. This fixed-end moderate-slope
rung retains selected axial inertia; omitted axial harmonics equilibrate
statically. Axial/normal modes and time require refinement. It adds no stick
camber, moving contact, nonlinear viscous cross terms or real-instrument claim.
Focused G1/G3 checks compare independent continuum strain energy and gradient,
quadrature-width invariance, normal-to-axial energy transfer, and rollback.

### Material gas loading in the file-driven plate renderer

The optional `radiation-load none|baffled-piston` record after `observer` in a
plate performance selects the existing compact piston reaction. The one
resolved ambient gas and output clock supply both mechanical loading and the
retarded pressure observer. Physical face area, modal mass and signed volume
and force projections remain those of the original chart reduction. The input
hash and CLI sidecar retain the load selection. Omission preserves the original
unloaded image; no material loss or radiation model is inferred.

`CompactPlateVoice::from_radiators_with_piston_loads` exposes the same path to
code consumers. A retained filter owns its original gas snapshot and sample
clock: mismatched consumers or step durations refuse before filter mutation;
reattachment of the identical load preserves its history. Gas snapshots pass
the material owner's coherence admission. This remains the existing sampled,
per-mode compact piston approximation, without mutual radiation, atmospheric
absorption, a whole-fluid energy claim or continuum/experimental validation.
G3 regressions distinguish a pressure-induced change in mechanical response
from density-only output scaling, check the actual retained body's kinetic
energy, and compare loaded file/WAV renders across callback partitions. The original piston
decay/frequency oracle retains its tolerance and uses a coherent low-pressure
gas for its vanishing-load limit.

### Material-driven finite-body sliding

`tribothermal::dynamics::SlidingPair` adds force-driven tangent inertia to the
original material sliding contact and finite enthalpy network. Each complete
uniform body's original thermal-cell mass also supplies its inertial mass.
Both velocities and applied forces have two components in the original fixed
tangent frame, ordered by body A then body B. Centre-of-mass momentum and two
Newton--Krylov relative-velocity residuals determine the complete planar motion;
every trial restarts the original contact/thermal candidate,
including its actual endpoint temperature, pressure, speed and slip-history
source queries. Separate impulse and whole-system energy gates precede joint
publication of velocities, original slip/work history and finite temperatures.
Original friction work is an internal transfer; external force work is measured
as the dot product with each body's endpoint displacement. Physical friction alone is partitioned
as heat, while backward-Euler kinetic loss remains an explicit numerical account.
The fixed patch/normal-load reduction passes released motion and inertial
compliance `dt²(1/mA + 1/mB)` to the original tangential return map. That owner
selects exact sticking or vector sliding; the material source queries the
Euclidean speed, and accumulated distance follows actual slip arc length.
The source motion converges to that same two-component native slip. A stuck
branch queries exact zero speed
and leaves slip/work history unchanged. Bodies can accelerate together, stop
and restart in the opposite direction, under the native endpoint static-first
convention. It does not locate the physical stopping instant within a coarse
interval: a backward-Euler stop may have numerical kinetic loss and zero endpoint
friction heat. Refine time to resolve stopping heat; numerical loss is never
thermalized. Zero-speed material support is required. No normal, rotational or
contact-geometry evolution is added. Focused synthetic regressions cover static
capacity, signed source feedback, oblique motion, rotation of the tangent frame,
stopping/restart, continuation and rollback;
native execution remains required for qualification. The internal same-trial
friction callback allows other energy owners to share this original solve;
compliant consumers must gate queried motion against the native solved slip.

The `SlidingPair::from_material_bodies` constructor also binds the original
`fs_mbd_material::ResolvedBody` objects to their ordered interface materials
and finite caloric cells. It requires uniform same-card bodies with matching
mass and total reference volume. Motion and kinetic storage consume the
retained canonical body inertias; accepted candidates share those immutable
sources. A density change at fixed geometry therefore changes both inertial
response and the matching finite heat inventory through data resolution.
Mismatched identities, masses or geometry refuse before motion advances.

`SlidingPair::step_planar_candidate` advances turning and stopping motion in
that same fixed tangent frame through the original material/contact/enthalpy
candidates. Each segment uses midpoint displacement, with a scalar radial
speed solve for the current isotropic friction law. Physical traction work
heats the original finite bodies; no backward-Euler kinetic defect is booked
as heat. Static reaction is considered only at an accepted zero relative
velocity. A stopping-time root must also satisfy the full vector impulse
gate; unresolved turning is subdivided instead of accepting a mirror-velocity
midpoint root. Both bodies' accumulated impulse defects and the whole-call
energy account must pass before any segment is published.

Constant collinear Coulomb motion is exact at arbitrary tangent orientation,
including arrest and restart. General curved paths require refinement of the
declared maximum segment duration, independently of root accuracy. An oblique
limiting stop is approximate within its vector impulse/time budgets. G1/G3
tests compare curved-stop integrals and a separate continuous RK4 turning
reference on three time rungs, and exercise rotated frames, distinct static
and kinetic capacity, source temperature feedback, and transactional refusal.
Temperature-dependent coefficients still use endpoint thermal states and
segment mean speed; thermal transport remains first order. No normal motion,
anisotropic material law, within-stick thermal event localization, or real
material calibration is inferred.

### Directional material friction with finite heat and inertia

`tribothermal::directional` retains the original four directional coefficient
claims, ordered material/phase identities and finite caloric cells. Each thermal
trial refreshes declared surface/mean temperatures, nominal pressure, sliding
speed and unsigned path; undeclared source conditions stay fixed. The native
associated ellipse supplies actual traction/slip work, partitioned explicitly
between the two bodies and any named outside destination.
`DirectionalPair::from_material_bodies` also retains the original source mass,
density and geometry receipts.

Force-driven planar candidates use midpoint compliance `h²/(2 m_reduced)` and
the native scalar ellipse projection. Only speed/path-dependent sources need
an additional bounded scalar consistency root. Accepted source path never
exceeds actual accumulated motion. Static capacity applies only from accepted
rest; moving segments use kinetic coefficients. Full-vector arrest, turning
subdivision and rank-one partial arrest retain motion along a frictionless
material axis, including widely separated component scales. Independent
per-body impulse and whole-system work/heat gates precede joint publication;
cancellation or exhausted budgets leave all original owners unchanged.

Fixed normal load, nominal area and tangent/material frame are required.
General turning and temperature/speed/path dependence require timestep
refinement; endpoint thermal transport is first order, and arbitrary weakening
is not claimed to have a unique solution. No normal, rotational or
evolving-patch dynamics or experimental calibration is inferred. Nine focused
synthetic test groups cover source ownership, finite inertia/heat, independent
continuous turning refinement, static/kinetic separation, temperature and motion
feedback, partial arrest, monotone continuation and atomic refusal.

### Material contact loads on a field actuator

`tribothermal::field_motion` closes the existing dry-contact active-set/enthalpy
owner with `MaterialFieldMotion`, including optional sourced axial support.
Two displacement residuals identify actual field gap change, original solved
contact slip, and the motion used to query interface material data. The native
traction supplies the opposing field force. Field/circuit/body/support history
and contact/thermal history publish together; inertia is not integrated twice.

The same thermal trial resolves bulk electrical transport, pinned electrical
contacts, and pinned mechanical friction. Original conductor/contact losses and
the original friction heat partition enter explicit receivers exactly once.
External thermal supply is independently evaluated before either internal loss.
Field, kinetic, support elastic, interface and thermal storage plus exported
friction work must close original electrical, mechanical and thermal supplies.
A single descriptor budget includes all discarded outer and thermal trials.

Surface A must match the single-region field body's source card, mass and
reference volume; surface B is a stationary finite thermal guide. Normal force,
area and tangent frame are fixed. Static reaction at rest uses an exact zero
speed/history source query; uninterrupted sliding in either direction is
supported. An unresolved stop/reversal refuses;
there is no velocity clamp, impact impulse, moving thermal mesh, or new normal
contact law. Focused `tribothermal::field_motion` regressions cover independent
constant-charge trajectories, both static reaction directions, body/source
inventory matching, authored/callback external heat, caloric/resistance/friction
feedback including electrical contacts, cancellation and global-budget refusal.
Syntax parsing and independent numeric checks passed; native execution is pending.

### Typed gesture consumer bindings

`gesture_path::bind_paths` binds each description-validated parameter path to
an explicit existing pressure voice or actual material-string bow member.
Every source path requires a destination; incompatible/unconsumed paths return
`GestureError::NoConsumer { path }`. `BoundPaths` retains the original resolved
v2 schedule, alias-migration identity, description and exact destination map for
render provenance. V2 schedules remain typed throughout lowering and call the
original pressure/friction sampling loops without changing time rounding, ramps, station
checks, integration, cancellation, force transfer or stored physical state.
Numeric u8 alias decoding is isolated in fs-scenario's compatibility module;
typed destination members use usize without a legacy u8 conversion.
The existing single-string, material-bank and simultaneous-bow entry points have
path-bound counterparts returning that provenance alongside their renderer.
Legacy file entry points remain v1 compatibility adapters; a typed schedule
requires description admission and explicit path bindings. No new consumer or
silent fallback is invented for topology events, frets, valves, pitch, prestress
or terminations. G0 `gesture_schedule::unconsumed_targets_reach_a_typed_consumer_or_refusal`
checks every alias and bitwise pressure-control equivalence; original render
equivalence tests remain unchanged. Addressing alone establishes no physical
validation or conservative topology-change claim.

### Declared per-note control resolution

`description::per_note` resolves `MinimumThresholdPressure` using a model-only
monotone onset predicate and explicit SI domain, pressure tolerance and maximum
evaluations. Its API takes no audio/reference input. The non-oscillating lower
and oscillating upper bracket, every probe, model identity, objective, target,
note, domain and budget enter a deterministic content-addressed receipt.
`record_resolution` binds the receipt to the corresponding Performance control;
`fs_scenario::performance::lower` retains it in render provenance. Other declared
objectives remain representable but refuse this pressure-specific resolver.
Model identities must bind initial conditions and the predeclared observable;
evaluation cancellation is the model's responsibility. The search is bounded
and synchronous and publishes only a completed result. Monotonicity is an
explicit model obligation; a sampled bracket is not a proof of global
monotonicity or asymptotic bifurcation onset. G1 `per_note_rule::pn_001` uses
the existing reed fixture's 0.25 s sustained-pressure observable with 18.75 Pa
bisection width; G0 `pn_002`–`pn_004` cover API shape, provenance and refusals.
The authored four-note `music_lane.sh --fixture midi-reed` exercises MIDI,
lowering, model-only resolution, schedule rendering, replay and named refusal.
This fixture is plumbing evidence and makes no reed realism claim.

### Port-schema adoption contract

`PortedModule` exposes stable local schemas, coordinate charts, internal N-way
junctions and the last accepted physical-time transfers for CoupledModalSystem,
SharedCavity, GasPiston, DoubleActingGasPiston, DynamicAperture, ModalAperture and
the ReedBoreVoice characteristic inlet. All schemas pass PortSchema::try_new.
Graph embedding must namespace local identities. The chart timestamp names the
clock origin; each audit binds its own exact accepted interval. Held effort work
uses the original conjugate displacement increment, never logical ticks as time.

Adapters observe the existing solvers; they do not feed back into their arithmetic
or acceptance. Typed transfers and window evidence are materialized lazily on
request. Exact-network connection work is minus the spring-energy change and
the independently integrated damper loss along the same full trajectory, never
endpoint reaction times displacement. Multiple damped columns require bounded
per-column integration scratch; this is not an allocation-free observer.
Rejected trials preserve the preceding account. An independent observation
refusal is retained as `ModuleObservationFailed`, and observational overflow as
a typed transfer error; neither publishes an invented zero record or rewinds
an already accepted step.

SharedCavity retains its mass-normalized modal chart and signed volume projection.
The canonical schema is the explicitly transformed reference-coordinate pair:
F = f_q sqrt(M_ref), v = q_dot/sqrt(M_ref), with M_ref=1 kg, not a specimen modal
mass or a point-force measurement. Raw owner units remain N/sqrt(kg) and
m sqrt(kg)/s. Its compression junction uses the original signed volume pullbacks
to pressure/volume-flow coordinates, including neck and compression storage.
Each CoupledModalSystem column similarly names its B_i^T q_i work charts and
one connection-storage endpoint. DoubleActingGasPiston's two gas faces share a
single common-flow wall junction; internal face work is not an exterior supply.

Available energy windows use the owner's actual storage change, prescribed work,
loss and existing joule tolerance. No balance is reconstructed as a leftover.
Piston heat ports use the source-backed discrete gas-midpoint temperature and
Q/T_midpoint entropy quadrature; gas entropy change remains independently
measured. This energy projection claims neither gas species/element balances,
exact isentropy nor external-reservoir thermodynamic validation. Negative
numerical entropy generation remains a visible WindowAuditReport violation.

MagentaBadger's ruling in Agent Mail 2636 explicitly leaves two audits unavailable:
DynamicAperture has no owner-declared energy tolerance, and ReedBoreVoice /
DelayedFilter has no physical line storage/loss owner or joule tolerance. Their
real transfers remain available, but ModuleWindowAudit::Unavailable names the
missing inputs. Those gaps keep IRP .162 open; unavailable is never a passing
window. No stepping or music-lane golden is intentionally changed.

G0/G3 tests in tests/port_schema_adoption.rs cover admitted schemas for all seven
owners, available/unavailable accounts, sign-flip refusal, two-port parity,
three-component/shared-cavity N-way descriptors, the common-wall descriptor,
and accepted-record preservation on refusal. Each module emits a native-run
JSONL row with residual/tolerance or explicit null/unavailable fields. Native
verification and reed/string baseline hashes remain central-run requirements.

### Integrator doctrine and current audit boundary (IRP M5.2-W1)

Render-path linear parts must use exact held-input propagation or explicitly
frequency-matched energy-consistent propagation. Nonlinear storage retains its
discrete-gradient owner; generalized-alpha is not admitted on render paths.
Held effort exchanges use the same displacement work at both ports, rate
changes require antialias filtering, and phase/contact/solver budgets must be
admitted explicitly. A midpoint path that remains is admitted only inside the
fs-phs one-cent band (`admit_midpoint_band`, pole-rate bound times dt <= 0.04);
the table below says which paths are exact and which are banded.

`tests/integrator_warping_audit.rs` contains the current path inventory and G1
one-step-map measurements at 1/5/10 kHz and 48/96/192 kHz. The common fs-phs
linear owner is measured for its callers; that is not a complete nonlinear or
connected render. Public ImpactSystem and Dynamic/Paired/ModalAperture seams
are exercised directly, including the uniform scalar condensation. The exact
modal owner is the control. Midpoint phase is compared to
`2*atan(omega*dt/2)/dt` within 1e-6 cents. The closure test
`every_render_path_is_listed_with_integrator_and_measured_warping` accepts
a row only if it is exact, admits no more than the one-cent band, is the
first-order network row with its measured order ladder, or is a PianoAudio
factor whose companion measurement stays within 1 cent. A path outside its
band refuses with a typed `TimeResolution` error; nothing is retuned.

| Current path | Linear owner | Current band / audit boundary |
| --- | --- | --- |
| KC string (`realize_coupled_kc`) and nonlinear plate (`thin_plate::VkBody`) through fs-phs `step_render` | Exact ZOH when the storage is exactly quadratic, otherwise Gonzalez midpoint | Midpoint admits only the one-cent band on the split's pole-rate bound and refuses with `TimeResolution` |
| Dirac joins (`realize_dirac_join` / `step_descriptor`) | Gonzalez midpoint | One-cent band on the descriptor pole-rate bound; typed `TimeResolution` refusal |
| Glottis fold island, coupled ODE chains, blown and reed ODEs, brass lip held-force island, plate-cavity transformer | Exact ZOH (fs-phs `linear_exact`) | No warping; nonlinear forces are held over the step |
| PlateBank complete body/material/cavity/exterior state | Prepared full exact ZOH for quadratic storage, otherwise Gonzalez midpoint | Reciprocal loads share one accepted solve; the original nonlinear pole-rate band and complete-state capacity apply |
| TsDriver speaker coil/cone/panel | Full coupled exact ZOH | Held voltage, independently integrated current/work/loss; dense reference step, not an allocation-free callback |
| ImpactSystem | Gonzalez midpoint | One-cent band on the complete retained block (modal stiffness, volume-spring dyads and resistance), typed `ImpactError::TimeResolution`; the per-mode 0.9*pi guard remains a separate validity check |
| ModalAperture, uniform condensation, and Dynamic/PairedAperture with `BladePropagator::Midpoint` | Midpoint response | One-cent band (pole-rate bound times dt <= 0.04, DEF-04) refuses at construction and at every step; the material-modal adapter additionally caps its caller's bound at 1 |
| PairedAperture with `BladePropagator::ExactZoh` | Exact ZOH per elastic blade (scaled exponential; Van Loan viscous loss) | No warping at any omega*dt (audit rows ~1e-11 cents up to omega*dt 1.31). Contact/Bernoulli forces are held over the step and solved at the midpoint opening; blade material memory refuses |
| DynamicAperture with `BladePropagator::ExactZoh`, authored or plate-bound (closure, Maxwell arms), and every plate-valve performance file | Exact ZOH of the opening plus its Maxwell arms (scaled exponential; Van Loan viscous and arm loss), through `reed_bore::contact::step_with_opening` | No warping (audit rows <= 3e-10 cents up to omega*dt 1.31); arms checked against an independent RK4 integration. Contact/Bernoulli forces held over the step and solved at the midpoint opening |
| CoupledModalSystem | Exact ZOH of the full retained bilateral pencil (component modes, connection springs and dashpots), split into connected K/C islands (cd7f65bc7); nonlinear contacts keep held discrete-gradient reactions around the same affine response | Coupled algebraic Nyquist bound; measured against the exact two-body eigenproblem; the two-body sticking traction matches an independent exact-flow integration within 1e-10 N |
| ModalAcousticTimeModel | Exact ZOH | Caller-declared Nyquist guard; no frequency-warping claim for its observer transfer |
| Source-backed GasPiston | Midpoint gas/momentum | One-cent band on the gas-spring rate sqrt(rho c^2 A^2/(m V)) at the start, midpoint and end states; typed `PistonError::TimeResolution` |
| Circuit/thermal/network piston | First order in general | Observed order 1 +- 0.1 on the four-rung companion ladder (`FS_IRP_NETWORK_AUDIT_LOG`); no exact-linear claim |
| Live PianoAudio at 1/2/4x | Exact local blocks plus finite-step coupled solves | Companion-measured coupled-board frequency within 1 cent, window discrepancy included, at each factor (`FS_IRP_PIANO_AUDIT_LOG`) |

The coupled-mode test uses the exact symmetric stiffness eigenproblem and a
three-rate phase-error ladder. The coupled pencil steps exactly, so the scalar
midpoint atan oracle does not apply to it. Piston tests reuse the actual
source-backed thermochemical and thermal/electrical fixtures. The speaker
migration has a held-voltage time-partition regression and a public-state pole
audit. The audit and its companions pass natively (M5.2-W1 closed on ts1 b46 at
a4b6c4c22); consumer fixtures that request clocks outside a midpoint band
refuse with TimeResolution and must move to a resolved clock or an exact block.

- `PORT_SCHEMA_VERSION = 2`; `StableId` admits a canonical transport-safe ID.
- `PortValueShape` distinguishes scalar, non-empty vector/tensor, and field
  values; field shapes carry separate neutral `fs-iface::SpaceType` roles for
  effort and flow. `PowerPairing` distinguishes scalar product, Euclidean dot
  product, and field duality with explicit integration-measure dimensions plus
  `FieldMeasureSide` (effort-density or flow-density).
- `PortKind` covers mechanical translation, rotational torque/angular velocity,
  fluid pressure/volume flow, thermal temperature/entropy flow, electrical
  voltage/current, magnetic mmf/flux rate, and chemical electrochemical
  potential/amount flow. Each kind owns canonical generalized effort/flow
  dimensions; semantic kind identity remains distinct even when dimensions
  coincide (for example torque and energy).
- `CoordinateBinding` makes basis, frame, and positive orientation explicit;
  `PortTimestamp` carries a named logical clock and tick;
  `ConservationRole` is canonicalized into a sorted, duplicate-free list.
- `PortSchema::try_new` binds stable ID, legacy/current `PortKind`, six-base
  effort/flow dimensions, value/field shape, coordinates, power pairing,
  timestamp, and conservation roles. V2 admission checks that the contraction
  matches the shape, effort × flow has watt dimensions without exponent
  overflow, the measure-adjusted generalized dimensions match the declared
  `PortKind`, and PR-2 schema-only kinds name their non-energy conserved flow
  where one exists. The original three seed schemas retain their PR-1
  Energy-only role vectors exactly.
- `PortKind::scalar_seed_schema` constructs one canonical scalar coordinate of
  any kind without inventing identity, frame/basis, or clock data; retained
  goldens use it to migrate the existing mechanical, fluid, and thermal seeds.
- `STREAM_PORT_VERSION = 1`; `StreamPort` is not a scalar `PortKind` or a fifth
  thermodynamic relation. It bundles signed mass (`kg/s`), canonically ordered
  species/element amount (`mol/s`), three momentum-rate (`N`), energy (`W`),
  and entropy-rate (`W/K`) values under one `StreamChartBinding`. Its fixed
  roles are Energy, Mass, Amount, LinearMomentum, and Entropy.
- `StreamChartBinding` separates the constituent-basis artifact and its
  explicit, canonically ordered species/element axis from the spatial
  basis/frame/orientation. It also binds the state schema, chemical reference
  state, logical clock/tick, gravity datum, and the closed
  `StreamStressWorkConvention`. Stream admission is owner-outward only until a
  public pullback exists.
- `StreamEnergyChart` structurally selects exactly one of: the canonical
  moving-stream enthalpy chart
  `mdot * (h + |u|^2/2 + g*z) + W_deviatoric`; internal energy plus pressure and
  deviatoric Cauchy work; or one coordinate from an exact mixture
  Euler/Legendre family. Internal-energy, enthalpy, Helmholtz, and Gibbs
  selections all reconstruct canonical enthalpy from the retained conjugate
  terms before transported energy is formed. The caller-declared stream energy
  rate must equal the selected chart bit-for-bit.
- `PressureWorkCrosswalk` recomputes `h = e + p/rho`,
  `volume_flow = mdot/rho`, and `mdot*(h-e) = p*volume_flow` exactly.
  `EulerLegendreCrosswalk` recomputes the mixture Euler identity plus enthalpy,
  Helmholtz, and Gibbs transforms. Both also require an
  `ExactIdentityProofRef` bound to the complete stream context.
- `ChemicalEnergyAccounting` admits either energy embedded in the selected
  state potential or one separately proved species-potential power term. It
  refuses dual ownership, a foreign chemical reference, and a foreign or
  wrong-kind partition proof. Explicit species-potential and Euler/Legendre
  modes require a species-only constituent axis.
- `ConservativeJunction::iterate_added_mass_fixed` /
  `iterate_added_mass_aitken` are schema-bound migration bridges for the
  retained mechanical scalar fixture; they must remain bitwise equal to the
  legacy results and are not general FSI operators.
- `ConservativeJunction` admits two conjugate schemas and evaluates the finite
  scalar seed as `SchemaInterconnection`/`SchemaPort` (shared effort, opposite
  flow). It localizes metadata mismatches and refuses non-scalar/non-finite
  evaluation.
- Four closed primitive variants are distinct in `PortPrimitive`:
  `ConservativeJunction`; `StorageElement` with Hamiltonian/free-energy state
  schema and constitutive-gradient operator; `DissipativeRelation` with a
  constitutive family plus mandatory monotonicity/nonnegative-production
  evidence reference; and `SourceOrReservoir` with an explicit signed
  `AccountingBoundary` carrying the same basis/frame/orientation as its port
  plus included-source/external-reservoir treatment.
- `ACCOUNTING_WINDOW_VERSION = 1`; `AccountingWindowInterval` names an exact
  nonempty interval in one logical clock without treating tick differences as
  seconds. `WindowEvidenceRef` binds a receipt/verifier/statement digest to
  that interval and to one exact `WindowEvidenceRole`: manifest closure,
  manifest row and local port, initial/final inventory, integrated row and
  local port, boundary temperature, element projection, or charge projection.
- `WindowElementSchema` either declares a nonchemical no-claim or retains a
  nonempty canonical element axis plus interval-bound projection evidence.
  Stream and chemical-power manifests require the audited form.
  `WindowChargeSchema` requires direct-coulomb, species-to-coulomb, or
  proven-neutral evidence; neutrality admits only exact zero charge values.
- `WindowManifestEntry` retains contribution, exchange, local port,
  counterparty, topology role, port kind, logical timestamp, optional external
  `AccountingBoundary`, and row-binding evidence. Stream entries additionally
  retain the complete `StreamChartBinding`; every stream in a window must use
  the same state schema, constituent basis/axis, chemical reference, spatial
  coordinates, logical clock, gravity datum, and stress convention. Internal
  exchanges are admitted only as two reciprocal rows of one kind and timestamp.
  Manifest-closure evidence carries the complete canonical row subjects,
  including topology, boundary, and retained stream context, so it cannot be
  rebound to a different or truncated same-interval manifest.
- `WindowBalance` and `WindowInventorySnapshot` carry finite typed energy,
  mass, canonical elemental amounts, charge, and entropy at an exact endpoint.
  `IntegratedWindowTransfer` uses one convention for every signed value:
  positive into the audited system. `BoundaryEntropyBreakdown` makes advected,
  diffusive/chemical, heat, and radiation terms mandatory and retains a
  contribution-bound constant/profile temperature reference whenever heat or
  radiation is nonzero.
- `AccountingWindowSpec` owns the exact manifest, a common stored/stream energy
  reference contract, element/charge/entropy conventions, typed nonnegative
  tolerances, optional reference environment, and closure evidence. An empty
  manifest is legal only as an explicitly evidenced isolated window. Filled
  rows must cover the declared manifest exactly; ordinary power rows may carry
  only balance axes compatible with their `PortKind`, while stream rows carry
  the full multi-balance vector.
- `WindowAuditReport::audit` returns `Err` only for malformed, incomplete,
  foreign, nonfinite, or semantically incompatible input. A physically failed
  audit returns `Ok(report)` with retained `WindowAuditViolation`s. It computes
  `R = (final - initial) - sum(transfer_into)` for energy, mass, each element,
  and charge; computes `S_gen = delta(S) - sum(S_boundary,into)`; applies
  absolute typed tolerances to equalities and a one-sided lower bound to
  entropy generation; and retains canonically ordered per-port and per-exchange
  ledgers. Every internal exchange is independently gated for zero net energy,
  mass, each element, charge, and entropy transfer, so defects cannot cancel
  between pairs.
- `ExergyEnvironment` and `ExergyLedger` optionally report the narrow
  Gouy-Stodola diagnostic `T0 * S_gen` for a finite positive explicit `T0`.
- The raw scalar `Port` remains a backwards-compatible, non-admitting numeric
  container. The legacy `conjugate_to`/`interconnect` migration oracle composes
  only the original mechanical, fluid, and thermal seed kinds and refuses all
  schema-only kinds; raw construction or arithmetic alone carries no
  dimensional, coordinate, clock, identity, or conservation certificate.
- `interconnect(kind_a, kind_b, effort, flow) -> Result<Interconnection,
  CoupleError>` — a Dirac structure (shared effort, opposite flow) whose
  `interface_power` is `0` exactly (power-conserving by construction); refuses
  incompatible ports. `interface_power(&[Port])` = `Σ effort·flow`.
- `EnergyAudit` — `record`, `max_generation`, `is_passive(tol)`: the legacy
  `is_passive` name checks only caller-supplied scalar interface imbalance at
  each recorded exchange. A nonzero balance is a bug alarm, not a proof of
  whole-system or closed-window passivity.
- `AitkenRelaxation::new(omega_init, omega_max)` + `next_omega(residual)` — the
  scalar Δ² dynamic relaxation factor, magnitude-capped.
- `iterate_fixed_relaxation` / `iterate_aitken` — the added-mass interface
  fixed point under fixed vs Aitken relaxation → `FsiResult { converged, steps,
  solution, final_residual }`.

### `vibroacoustic` (bead frankensim-fsim-vibroacoustic-wgkq7)

The crate's first numerical engine: modal structure x rigid-wall-cavity
x exterior-radiation coupling, frequency domain under `e^{-i omega t}`
(matching `fs_bem::helmholtz`), hysteretic stiffness `k (1 - i eta)`.
Sign convention: positive structural deflection points AWAY from the
cavity (composes with the BEM's outward panel velocities).

- `StructuralModes` (mass-normalized, `fs_modal` convention) /
  `CavityModes` / `AcousticMedium` — the data carriers; any cavity
  basis plugs in, `rectangular_cavity_modes` (analytic, with a
  sufficiency-PROVEN enumeration box — an isotropic heuristic silently
  dropped axial modes of elongated cavities, executed regression) and
  `helmholtz_resonator_mode` (lumped lowest mode, `2 (8/3 pi) a` end
  correction; the psi = 1, Lambda = V presentation is EXACT, proven by
  the 3-state reduction) ship.
- `assemble_coupling` — area-integral `C_rq = INT phi_r psi_q dA`,
  pinned against closed-form sin x cos overlap integrals.
- `project_radiation_impedance` — panel-space `Z` (`p = Z v`) reduced
  to modal `Zm` with area weighting on the force side only.
- `VibroacousticModel::with_structural_viscosity` — attaches the existing
  `ModalPlateViscosity` PSD owner in the same mass-normalized modal basis.
  The full signed matrix enters the coupled structural dynamic stiffness as
  `-i omega C`; passive-channel squares supply `1/2 Re(v^H C v)` to the
  structural power account. Cross-mode damping and undamped nullspaces are
  retained, including each principal block used by the truncation solve.
  Scalar structural hysteretic loss must be zero when this matrix owner is
  attached; cavity loss and exterior radiation retain separate accounts.
  No basis reprojection or frequency-dependent relaxation law is inferred.
  `src/vibroacoustic/viscosity_tests.rs` checks an independent analytic
  two-mode/cavity response, integrated cycle work with radiation, a relative
  dashpot's common-motion nullspace, signed-basis invariance, and truncation.
- `VibroacousticModel::with_cavity_relaxation` — binds the original gas's
  molecular reservoirs to the cavity through `fs-material::relaxation`'s
  complex bulk modulus `K(omega)`. Density and equilibrium sound speed must
  match the cavity basis exactly; scalar cavity loss must be zero. The pressure
  mass term becomes `omega^2 rho0 c0^2/K(omega)`, with the original reciprocal
  inertial boundary forcing retained. Molecular heat is independently evaluated
  as `omega Lambda Im(1/K) |p|^2 / 2`, including uniform sealed-cavity compression.
  This is a homogeneous, quiescent molecular constitutive law. Shear viscosity,
  thermal diffusion, wall boundary layers and mean flow are not inferred from a
  propagation attenuation coefficient. Per-frequency caloric work is bounded by
  the original reservoir count and performed once outside the modal loops.
  `src/vibroacoustic/gas_tests.rs` checks a sealed gas spring, independent Schur
  elimination, integrated pressure/strain cycle work, reciprocity, modal
  rescaling/truncation, explicit humid-air species and stale-gas/double-loss
  refusal. The undamped eigenfrequency diagnostic remains the equilibrium-gas
  reference; loaded dispersive resonances require the FRF.

- `VibroacousticModel::frf` / `frf_truncated` — dense complex coupled
  solve; every response carries the complete power breakdown with BOTH
  audit residuals (`input = structural + interface + radiated`,
  `interface = cavity`) — algebraically exact identities, measured at
  1e-13.
- `frf_with_convergence` — truncation error OBSERVED (full vs halved
  bases) and refused (`FS-COUPLE-VIBRO-TRUNCATION-NOT-CONVERGED`)
  above the caller tolerance; near-resonance rows are genuinely
  truncation-sensitive (measured 0.28 on the casebook fixture) and the
  casebook exposes rather than hides that.
- `undamped_natural_frequencies` — exact `x = omega^2` linearization
  (unit-triangular `B`), pinned by the two-oscillator closed form at
  ~1e-16 across a five-decade coupling sweep, a first-order added-mass
  falsifier (quartering 4.00), and a multi-mode independent-determinant
  oracle.
- `VibroError` — stable `FS-COUPLE-VIBRO-*` refusals.

### `modal_acoustic_time`

A re-export of `fs_modal::time` (IRP M1.1-W1a, owner decision OD10): the exact zero-order-hold modal realization, its types, caps and
semantics are specified in `crates/fs-modal/CONTRACT.md` (`### time`), and every former path `fs_couple::modal_acoustic_time::*`
still resolves. fs-couple keeps `modal_acoustic_time::broadband` (the opt-in broadband receiver observation; it needs `fs-vfit`,
`fs-phs` and `broadband_radiation`, so it cannot live in L1). Consumer facts that stay here: `AcousticRealizeError::Modal` wraps
`ModalAcousticTimeError`; the bowed runtime's set-valued contact uses the free-motion plus unit-force endpoint response, cached per
bow station (see the Bow bullet under the string paths); the plate and bowed-string paths call the single-mode steps
`advance_exact_zoh` and `ModalAcousticTimeModel::advance_exact_zoh_step`, which are public in fs-modal for that reason.

### `render`

Block render API (music bead `frankensim-music-v8-root-3ez8g.2.1`): the
callback-shaped hosting layer for performance images. `RenderContext`
owns admitted voices (`ReedBoreVoice` — the `realize_reed_bore` physics
restructured for blockwise stepping, with the one-shot realizer now a
thin wrapper over it so the two paths share one loop body;
`ModalStringVoice` — the exact-ZOH modal runtime hosted verbatim) and
advances them one block at a time into a pre-sized scratch buffer.
INVARIANTS: a block boundary performs no arithmetic, so block size is a
pure loop-partition choice — the one-shot path and every block partition
are bitwise-identical (tested at 64/571/full); typed `ControlDelta`s
apply only BETWEEN blocks, transactionally.
`SetModalForce` changes one mass-normalized held force [N/sqrt(kg)] in an
admitted modal voice, preserving displacement, velocity and its physical clock.
Zero force releases the input without silencing or resetting the resonator.
All voice/mode indices and force values in a batch are checked before any
input or control log changes. This exposes generalized excitation, not a
spatial pluck/contact law or an automatic score scheduler. Controls are logged as
`ControlRecord`s carrying their D17 lift description (empty = pure input
move); cancellation (`render_under_gate` over `fs_exec::CancelGate`) is
polled only at block boundaries with drain semantics — an in-flight
block completes, the rendered prefix is a whole number of blocks, and a
resumed context continues bitwise-identically (tested). ALLOCATION:
the massless-reed voice with an empty plate bank is allocation-free per
block (counting-allocator gate); the massive-reed lay path
(`dissipative_modal_forces` returns a `Vec` per sample) and the modal
voice (`step` builds a per-sample energy frame) are DISCLOSED allocating
voices — fusion candidates (bead 3ez8g.15), not silently admitted.
Determinism: one-host bitwise, inherited from the hosted kernels.
No-claims: no device/OS audio host (Franken-only policy; live output
would be a quarantined adapter decision); no image hopping yet (D17
lifts beyond input parameters land with the articulation/track beads);
refusal mid-block permanently poisons the context (mid-sample state is not
rewound). The failed block's output is partial and must be discarded.
Subsequent block or control calls return `RenderError::Poisoned` before
touching output, voices, or control logs; reconstruction is required.
Request sizing and control-admission refusals do not poison a healthy context.
`ScheduledRenderer` admits immutable `ScheduledControl` schedules.
`ScheduledRenderer::from_pressure_gestures` connects complete pressure tracks
directly to hosted voices, requiring their actual sample periods to match the
declared audio rate. It validates destination bindings even for an empty horizon,
requires context sample zero, and preserves supplied physical state. Compilation
and all admission finish before controls or simulation state advance.
The renderer
splits caller blocks exactly before their integer sample indices. Equal-index
events retain input order; end-boundary events wait for the next nonempty block.
All controls are validated and log capacity reserved before stepping. The
constructor refuses invalid schedules or budgets; no resampling or
score-to-excitation compilation is claimed. Its cancellation gate polls before each
caller block and consumes no events on cancellation. G3 tests compare scheduled
force changes, release and re-excitation against direct modal stepping across
multiple block partitions, with cancellation/resume at callback boundaries.
`pressure_gesture_controls` binds a selected blowing-pressure gesture track to
that scheduler. It samples existing gesture ramps at the declared control rate,
holds each value, and places ticks at the first audio sample at or after their
time using integer clock conversion. Admission bounds evaluated ticks and storage;
the caller supplies the destination voice's audio rate and a half-open horizon.
G3 checks a nonintegral clock ratio and pressure ramp against per-sample control
application across callback partitions. Other gesture targets and continuous-time
interpolation between control ticks are not implemented by this binding.
Boundary facts a consumer must
know, both discovered by this bead's battery: the characteristic-line
realization refuses `UnflangedOpen`/`FlangedOpen` terminations whose
Nyquist `ka` exceeds 1 (at 48 kHz that is every bore wider than
~2.27 mm radius — the measured-load `MouthLoad::Tabulated` lane
is the lift), and the Zwikker–Kosten Bessel path's former
`r_v ≈ 160` continued-fraction ceiling was an iteration-budget bug fixed
in `fs-phs` (budget now scales with `|z|`).

The scalar characteristic-line constructors, `realize_line`, `impedance_line`,
regional/material-wall lines and all duct/field reed-voice constructors accept
the shared `MouthLoad` selection. Existing `Termination` arguments retain the
analytic path. A tabulated load requires an explicit `TabulatedMouth`: the
mouth plane/radius and exterior gas must match, and its real DC impedance is
supplied separately rather than extrapolated. Every positive-frequency bin
comes from the table, with no analytic splice or out-of-band fallback. The
complete declared FFT/ladder band is admitted before realization. Raw terminal
data alone cannot supply the bore's static behavior: the gas/wall propagation
must also admit the existing 1e-3 rad/s quasi-static probe. A material wall
whose source band excludes that probe explicitly refuses, even when its
positive FFT band is covered; no static wall response is extrapolated. Raw terminal
passivity is checked before propagation, so bore loss cannot conceal an active
mouth; characteristic FIRs also retain the 4x-oversampled passivity refusal.
The existing Bessel propagation, FIR/ladder and reed stepper remain the numerical
owners. Scalar lines refuse `Matched`, which is a multimodal continuation.
The full-band regressions compare characteristic reflection and mouth-pressure
transfer with frequency-domain TMM at 48/96 kHz, plus uniform-regional identity,
changed-load sensitivity, physical-reference/band refusals and active-load
mutations. They use an explicit resistive control, not a validated exterior
radiation bake; table convergence, complete-instrument realism and multimodal
runtime adoption are separate work.

### Aperture-junction solver modes (bead frankensim-2s4i5)

IRP M5.5 adds `render::ladder::RealizedResolution` and
`ReedBoreVoice::{new_with_fft_length, realized_resolution,
set_junction_step_tolerance, junction_residual_max, initial_ir_tail_energy}`.
The constructor shares the existing `with_line` admission/initialization; it
only passes the explicit request to `characteristic_line_dense`. Actual FIR
history length exposes the 8192 clamp. Strict mode has **no** pressure-step
tolerance: it still terminates on an exact zero or adjacent binary64 endpoints.
The massless FastNewton seam admits an explicit relative tolerance in (0,1],
keeps 1e-9 as the existing default, and preserves every corner/fallback rule.
Massive reeds refuse this tolerance seam rather than reporting a tolerance the
coupled structural solver did not use. The observed volume-flow defect is
recorded in m3/s, and the initialized periodic FIR's last-quarter tap norm is a
dimensionless diagnostic. Neither bounds pitch, onset, or partial-level error.

`tests/convergence_ladder.rs` binds the existing M-CL0 fixture to fs-ladder,
with at least three measured sample-rate, FFT, and Newton-tolerance rungs and
the joint group at fixed physical IR duration. The legacy observer's backward
jet-force difference supplies the declared smooth-signal first-order hypothesis;
this is not a global-order assertion for the oscillator. The fixture measures
raw pitch/partials and a **finite-horizon onset** (0.25s, 1e-3 Pa RMS, nondecaying
last-quarter RMS, 8 bisections in the declared 0..4800 Pa bracket). An unbracketed
onset refuses. This is not an asymptotic bifurcation/extinction measurement.
G1/G3/G5 report construction is exercised, but Q2 accuracy authority remains
unavailable where no conditioning map supplies a bound in the acoustic QoI's
units. Spreads and signed finest-reference joint interactions remain data.
`scripts/e2e/music_ladder.sh` runs a prebuilt native test executable and retains
JSONL, a raw-scaled PCM render/provenance sidecar and its fs-blake3-bound report
in a fresh directory. It performs no build and never overwrites evidence.

Both modes share the closed-slit characteristic constraint
`(1-r0) p_plus - p_minus_history = Zc Ubody`. The reflection denominator
is not clamped, and closure is evaluated at the body-loaded bore pressure;
body motion can reopen the slit. Singular or nonfinite candidates do not
qualify for this shortcut. G1 junction tests cover signed body flow,
reflected incoming waves, and reopening.

The reed aperture junction has two solver modes on
`ReedBoreVoice::set_solver_mode`. `Strict` (default) is the
deterministic path: a local sign-changing bracket, up to 32 symmetric
bracket expansions if needed, then up to 2100 bisections (covering the
binary64 exponent range, including subnormals). The former 96-step bound
could refuse near-zero roots solely because their pressure scale was far
below the initial bracket width. It stops at an
exact residual zero or adjacent floating-point pressure endpoints,
returning the endpoint with smaller flow residual. Nonfinite residuals,
failure to bracket, or exhausted refinement return errors; a grid guess
is never returned as a solution. This removes the former dimensionally
invalid pressure-scaled flow tolerance. Root selection is deterministic,
but multiple quasistatic roots and the Bernoulli dead zone remain model
limitations; this is not a uniqueness or exact-arithmetic certificate.
G1 tests independently reconstruct volume-flow balance across body flow,
reflection, and large incoming pressure, and test refusal without a root.
`FastNewton` is a DECLARED FAST MODE (the
fs-rand ziggurat precedent): DAMPED island Newton on the analytic
Jacobian of the shared Bernoulli-plus-wave-junction residual — each
step is halved up to four times until the residual magnitude strictly
decreases, and a sample that never improves, hits a sqrt kink or
non-finite state, exhausts its iteration budget, or shows a vanishing
slope is handed untouched to the strict path and counted in
`FastSolveStats`. Convergence requires an UNDOUBLED full step inside
`NEWTON_STEP_TOL`, so a still-descending line search never reads as
a converged root. It is NOT bitwise-equal to the strict path and must
never become the default until it earns the same proof. Receipts (bead
2s4i5): the junction-level battery lives inline (`fast_mode_tests`) —
a 64-point open/interior/closing stimulus sweep checks both modes against
a volume-flow scale `w H sqrt(2 Pc/rho)` and checks pressure agreement
within `1e-8 * (1 Pa + |p_strict|)`. The integration lane
(`tests/reed_newton_fusion.rs`) owns the fallback-hit-rate receipt on
a full nominal render (bounded structurally at 0.5; measured 12-22%
warm-seeded) and before/after budget rows stamped with build profile —
debug rows are diagnostics forever, per the budget-lane doctrine.
Recorded finding: end-to-end rendered waveforms are NOT comparable
between modes beyond onset — the self-oscillating loop amplifies
per-sample root differences into phase-divergent (O(kPa)) waveforms,
so waveform equality is never a claim either mode makes.

### Material-driven viscous slit flow

`bernoulli_aperture::viscous::ViscousSlit` owns the shared reduced channel law
`dp = 12 mu L U/(w h^3) + rho U|U|/(2 w^2 h^2)`. `from_gas` retains one
coherent `GasState`, derives viscosity from it, and refuses a different density
at evaluation. `PairedAperture` and `ModalAperture` can select that law from
their retained material gas plus the physical channel length. `DynamicAperture`
accepts the same source-bound law; a plate-valve performance's optional
`slit_length` record uses its actual regional inlet gas, not the exterior gas.

All three mechanical owners use this flow inside their same-sample pressure,
contact and swept-volume solve and account for the channel loss as `dp * U`.
Independent blades and signed modes retain their own conjugate pressure areas;
the shared slit width is counted once. Spatial single/symmetric/independent
plate closures integrate the same law analytically along each positive P1 edge
segment, with an analytic pressure tangent. Closed portions contribute zero.
An area-equivalent height is never substituted into the cubic viscous term.
Spatial pressure condensation uses safeguarded Newton with a passive bracket,
96 evaluations and a 32-epsilon relative pressure-balance target; failure
refuses the trial before any mechanical, material-memory or wave state advances.
Law attachment is allowed only before the first accepted step.

G1/G3 tests in `viscous_aperture`, `modal_aperture::viscosity`, the paired spatial
module and the constitutive owner cover flow/work closure, material-gas swaps,
the cubic-gap laminar limit, independent point quadrature, partial closure,
edge subdivision, paired/modal equivalence and full modal damping. This remains
a fixed-state, parallel-channel, quasi-steady model with one common streamwise
length. It does not resolve transverse/entrance flow, channel inertance,
pressure recovery, gas evolution, compressible flow or measured instrument
accuracy; it is not the plan's resolved-flow acceptance rung.

### `music_render` (binary)

`ReedPerformance` accepts optional `slit_length LENGTH_M` immediately after its
reed record. `ReedBoreVoice::set_viscous_slit_length` binds the shared
`bernoulli_aperture::viscous::ViscousSlit` to the exact source inlet `GasState`,
including regional mixtures and axial fields. Both quasistatic and massive
reed/contact branches solve this constitutive law in the original same-sample
characteristic junction; no output filter, alternate contact owner or new
integrator is introduced. Omission retains the inviscid arithmetic. Selection
after rendering refuses. The massless viscous branch uses strict refinement;
fast mode counts this fallback and exposes no Newton stopping-tolerance claim.
G1/G3/G5 regressions in `render::scattering`, `reed_bore::fast_mode_tests` and
`render::schedule::reed::viscous_tests` cover the signed pressure/flow equation,
Poiseuille/closure limits, ten beating/reopening cycles with a pressure-work
account, inlet-material binding, and callback invariance. Quasisteady planar
flow does not establish unsteady channel, turbulence, or instrument realism.

The music lane's render CLI (bead ib15w): `music_render <reed|string>
<out.wav> [--seconds S] [--block N] [--full-scale-pa P]`. Pinned fixture
compositions render through the block API and encode through the ONE
pascals→PCM owner, `pcm_wav::encode_pcm16_wav` — the recorded seam
decision (beads ib15w + h7xu5.7.8): the cinematic stereo/receipt-hashed
encoder stays cinematic, no third RIFF writer exists in the music lane.
Reed runs accept `--temperature-k K` and `--ambient-pressure-pa PA`; these
feed the shared `GasState` into the existing duct/reed calculation.
Optional `--relative-humidity RH` supplies a dimensionless fraction in `[0, 1]`
to the shared `GasState::try_new_moist_air` mixture owner. Default zero follows
the exact dry-air path. Nonzero humidity changes the gas properties consumed
by the duct and reed calculation; it is not an output pitch/gain adjustment.
For nonzero humidity the owner refuses temperatures outside its saturation-fit
window and excessive vapor fractions. This does not model condensation or
moisture uptake in solids. The ambient values, relative humidity and gas
specification are recorded in the sidecar. Other fixtures refuse these
overrides. This is a stationary gas-state change, not
solid thermal expansion, heat transfer or phase change. Reed runs also accept
`--duct-length-m M` and `--duct-radius-m M` for the cylindrical duct, feeding
the declared geometry to the existing acoustic solver and recording it in the
sidecar. These dimensions must be finite and positive; other fixtures refuse
them. Optional `--duct-outlet-radius-m M` selects a linear cone when different
from `--duct-radius-m` (the inlet); omission or equal radii retains the cylinder.
Both cone endpoint radii are recorded for replay. Other fixtures refuse
them. Existing acoustic admission also applies: a positive radius can still
exceed the selected radiation model's `ka` limit and must refuse.
Geometry remains fixed during the performance. No automatic wall-material
response is implied. Gas-model validity
limits do not certify the fixed solid specimen at extreme temperatures.

The file-driven path also accepts `frankensim-flue-performance-v1`:
`music_render wind crates/fs-couple/examples/flue.performance OUT.wav --block 512`.
The fixed-order SI records explicitly supply mechanical/output clocks, the output
sample window and pascal scale/ceiling, bore and exterior gas, cylinder dimensions
and prescribed terminal reflection, mouth dimensions, jet card/gain/onset seed, an
optional opt-in jet response in place of the gain chain (`jet_response forced-bickley
SPACING WINDOW` or `jet_transfer nolle-price-2015 CHANNEL CHAMFER_DEPTH_M`, exclusive),
jet limits, wave memory, fractional phase band/budget, optional viscothermal cells
and band, receiver position/resolution/band, and control compilation limits. The
suffix is the existing canonical, typed blowing-pressure gesture schedule.
`jet_card tonal-interim` explicitly selects the retained stage-I source card;
`jet_card json {CANONICAL_JSON}` supplies a complete one-line card. No filename or
instrument name selects physical dimensions, material properties or calibration.
`gas` and `exterior` accept the same legacy, named-mixture and source-backed NASA-9
records as other wind sources; `exterior same` shares the supplied bore gas. A
separate exterior must have the same static pressure. Gas changes enter jet speed,
Reynolds number, propagation, selected viscothermal loss and physical radiation.
The optional `windway none`, `windway slit LENGTH_M` or
`windway unsteady LENGTH_M TRANSVERSE_CELLS CONVECTION_CELLS` record after
`flue` selects the shared ideal, quasi-steady or unsteady source-gas windway
before performance setup. Actual channel length, model, viscosity and resolved
cell counts remain in the output sidecar. Omitted/`none` retains the ideal jet.
A finite channel is explicit physical geometry; no lip-channel length is
invented for a transverse flute.

`windway slit-developing LENGTH_M` (opt-in, pexfp.135.4) holds the jet by
`bernoulli_aperture::developing`: a straight parallel-plate slit fed through a
uniform inlet, `P = rho U^2/2 [1 + 4 f_app L/D_h]` with Shah's (1978) apparent
friction factor for developing laminar flow, as stated openly by McCreery et al.
(NURETH-15 2013, eq. 2). No constant is tuned; at large `x+` it is the
fully developed slit law plus `0.674 rho U^2/2`. The held note's pressure (the
gesture's highest blowing pressure) must lie inside the correlation's validity
or the input is refused: `x+ >= 4e-4`; Chen's low-Re inlet term
`0.64 + 38/Re` may move the held drop by at most the stated 2.4 % (Chen's term
and Shah and London's agreement statement behind the `x+` floor are secondary,
search snippet only; primary not read); an
undeveloped exit keeps `Re_delta* <= 520` (core speed and displacement thickness;
the Blasius layer's linear instability, Manneville 2016); and a developed exit
keeps `Re_Dh <= 1760` (`Re_cl = 660` of the laminar centreline speed and half gap,
above which plane Poiseuille flow sustains turbulence, Kashyap et al. 2020). The
sidecar records `x+`, `Re_Dh`, `Re_L`, the core
speed, displacement thickness and these measures at the held pressure. The
attack and release pass through the same quasi-steady law and are not covered.
The jet transfers refuse it (their exit profiles assume Poiseuille or uniform
flow; pexfp.135.5 owns the core speed). A converging windway declared as its
exit-height slit most likely overstates the loss.

`jet_speed card-centreline` (opt-in, pexfp.135.5) after the `jet` record feeds the
jet card the exit's centre speed. By default the card reads the windway's bulk
speed, but its Strouhal number and band are defined on the centre (peak) speed of
a near-top-hat jet (`MeanJetProfile::u_centerline`). The rule is:

| windway | card speed |
| --- | --- |
| none | the bulk speed `sqrt(2P/rho)`; nothing changes |
| `slit-developing` | the inviscid core `sqrt(2P/rho)` while the layers have not met, then the Poiseuille centre `1.5 U` |
| `slit` | the Poiseuille centre `1.5 U`, refused where it exceeds `sqrt(2P/rho)` at the held pressure |

The card speed sets the convection `c_p = (c_p/U) U`, the growth band with its
band-lag tap, the onset seed and the card Reynolds number. The bulk speed alone
carries the flux, so the held flow is unchanged. The unsteady channel and the jet
transfers, which own their exit speed, refuse it. The sidecar records the law and
the held bulk and card speeds. Omission keeps every default result bit-identical.
No card selects it: under it the island's loop-phase window moves every flue card
off its pitch (pexfp.135.5 A/B), and that window is the open defect.

`bernoulli_aperture::slit_channel` resolves fixed, fully developed planar
Newtonian momentum from the retained `GasState` density and viscosity. Its
2..=4096 cell-centred finite volumes share a no-slip wall diffusion operator;
backward Euler evolves the transverse profile from rest under signed pressure.
An optional uniform kinetic-head exit is solved with the mean flow in the same
implicit update. Kinetic storage, pressure work, physical viscous heat, exit
work and numerical dissipation remain separate. Standalone refusal leaves the
accepted channel state intact. No fitted relaxation poles or alternate gas
formulas enter this owner.

`FlueVoice::with_unsteady_windway` connects that accepted flow to jet speed,
source volume, Reynolds coverage and the original card's growth law. Pressure
controls accelerate the fluid and release retains its momentum. A separate
2..=4096 cell semi-Lagrangian convection field transports the existing feedback
across the actual mouth cut-up at the evolving speed, including stop and
startup, without a minimum-speed history clamp. The labium feedback is read
at the frozen path's own delay (pexfp.320, 1718d79fb): the transit less the
band filters' low-frequency lag `2/omega_c + 1/omega_h = 16 theta / U`, with
its one-sample floor and sub-sample rule, as the distance `c_p tau` along the
field, so `16 theta c_p / U` upstream of the cut-up, independent of `U`. The
read interpolates linearly between the bracketing nodes (between the current
mouth velocity and the first node below it), the transport step's own
interpolant. A sub-sample transit and a jet at rest read at the labium. Before
this the transported loop was longer by `16 theta / U`, which is 88 to 114
degrees of jet phase at f0 on recorder a1/f1 and whistle D5. The voice's
original failure poisoning applies to the combined sample. Channel and acoustic mouth work
ledgers retain their distinct ports; numerical damping is not material heat.
Characteristics entering during a sample interpolate the accepted old boundary
velocity and the current solved mouth velocity at their physical entry time.
The affine current coefficient participates in the existing joint jet/mouth
equation; interior characteristics use the preceding spatial field. Zero speed
preserves interior parcels and retains the accepted boundary trace for restart.
Channel momentum, convection and bulk speed remain private candidates until
the mouth and wave step accepts. A rejected trial preserves their preceding
accepted state under the existing poisoned-voice policy. Convection uses a
piecewise-constant accepted bulk speed during each source sample. Focused
characteristic, source-history and gas-to-PCM regressions are authored; syntax
parsing and an independent characteristic/refinement oracle passed. Native
execution remains required.
Single-tube and regional-gas branched file sources use this same builder and
the same physical receiver. Independent steady Poiseuille and transient
sine-series regressions exercise space/time refinement, signed release,
reversal, gas interventions and callback partitioning; native execution is
required before treating them as passing evidence.

These are explicit numerical axes for a fixed, incompressible, fully developed
planar channel and the original estimated jet card. Side-wall/entrance/axial
profile development, compressibility, thermal transport, moving walls and
turbulence remain outside the channel. The exit uses uniform bulk kinetic
head despite the resolved internal profile. This does not establish a resolved
jet, moving-reed channel or experimentally validated instrument claim.

This loader constructs the existing `FlueVoice::from_fractional_network` and
`FlueOutletPerformance::with_multirate_radiation`. Every mechanical outlet-flow
sample reaches the shared differentiating Rayleigh receiver before the sole PCM
owner. The file declares 48 kHz CLI output and an integer mechanical/output ratio;
`--decimate` refuses a second observer. The sidecar retains both clocks, original
input/card identities, source gas coefficients and the actual receiver delay.
This is an **EstimateOnly** composition of a tonal jet and stationary rigid
cylindrical pipe with a prescribed reflection and one-way infinite-baffle outlet
radiation. The selected card keeps its narrower laboratory scope. It supplies no
resolved turbulence, matched radiation feedback, measured flute geometry, wall
material response or calibrated SPL claim; `--acceptance` refuses this source.

### `bakeoff`

Bake-off receipt schema and protocol for the music-program claims registry
(bead `frankensim-music-v8-root-3ez8g.1.4`). `BakeoffReceipt` records two
contender images run on one shared-card fixture against a caller-supplied
reference QoI map, with logical budgets (states, steps, solver iterations —
never wall-clock, so bytes are host-stable), observed failure modes, and a
reviewed outcome. `BakeoffOutcome` has exactly three variants —
keep-both, keep-for-subset, refuse-newcomer — and structurally cannot
express deleting an admitted image (doctrine D21: menus, not winners).
Receipts serialize to a canonical line-oriented byte encoding
(`frankensim-bakeoff-receipt-v1`, strict round-trip decoder) and are
content-addressed via `fs-blake3` under
`org.frankensim.fs-couple.bakeoff-receipt.v1`. The harness MEASURES;
the outcome is a reviewed judgment passed in by the caller — residual
tables cannot adjudicate claim scopes alone. The executed
`tests/bakeoff_string.rs` fixture (modal-ZOH vs the Gonzalez-stepped
pHS modal bank, analytic damped-oscillator reference) golden-pins its
receipt at `tests/receipts/string-modal-zoh-vs-phs-bank.bakeoff`; the
pHS contender's phase error is gated against the implicit-midpoint
dispersion model `(2/dt)·atan(omega·dt/2)` rather than a blanket state
band, because pointwise residuals near zero-crossings are unboundedly
phase-sensitive. No-claim: a receipt proves the comparison ran and what
it measured, never that a claim scope is correct — scope promotion is
the per-track gates beads' review job through the instrument-claims
registry.

### `broadband_radiation`

`broadband_radiation::modal` and `modal_acoustic_time::broadband` add an opt-in
one-way pressure receiver over complete mass-normalized modal velocity fields.
Computed training/independent held-out samples use exp(-i omega t), with explicit
producer identity, SI units, pure delay and producer discretization estimate.
Continuous vector fits are conjugated for the fs-vfit owner; proper stable
state-space filters join each original mode in an augmented held-force ZOH
exponential, reusing fs-phs's checked deterministic exponential. There is no
Tustin warp or velocity difference. The original modal owner still advances
mechanical states and work/loss records, preserving those bits exactly. Its
legacy pressure coefficients are zeroed only in the observer's private clone,
so the admitted receiver owns the pressure ceiling. A failed step publishes no
filter/delay/mechanical state. Block render failures poison their stream wrapper.

The pure delay uses exact integer taps or an explicit centered Blackman/sinc
fractional kernel; all taps must be causal within the physical delay. History
starts at zero. This approximates a bandlimited sampled pressure history, not
continuous retiming of all held-force spectral images. Held-out error includes
the actual kernel and retains EstimateOnly authority; it is not a continuum
certificate. Producer error remains separate. Constant samples with zero delay
use the exact legacy equilibrium-relative quadrature operation order.

Caps: all 4096 modal inputs are retained in <=256-input cold batches, each grid
has <=512 frequencies, <=4194304 transfer cells total, per-mode order<=128,
iterations<=128, delay-kernel half-length<=4096, delay history<=1048576 samples,
and total augmented matrices<=4194304 entries. No mode/pole truncation or implicit
budget enlargement is allowed. Host cancellation is at bounded pressure-block
boundaries; no hard-real-time, continuum passivity or load-model claim follows.

The v1 modal-file optional `broadband` block reuses its original parser and
ProjectedForceProgram; v2-v5 coupled/nonlinear trajectories refuse this
independent-mode receiver. CLI sidecars retain every producer/fit/error/delay
field, and the legacy path remains the default. The native G0/G1 battery is
`tests/broadband_modal_transfer.rs`: independent analytic held-input sphere
endpoint response at five off-resonance drives (1e-6 amplitude/radian bands),
constant compatibility, bitwise mechanical energy/work/loss, all-mode cap,
held-out refusals, integer delay, and three-rung delay/pole/grid ladders. The
pole/grid axes use the independent nonrational Bessel/Struve piston owner:
the delay-stripped sphere is already first-order rational and cannot establish
pole-count liveness. Numerical/native qualification is pending. The audible
replay/refusal sub-lane is `scripts/e2e/music_modal_broadband.sh`, requiring the
native qualification log and a supplied binary; it makes no instrument-realism
or general source-mortar/loading claim. M-BAR-Q2/M-GT-Q2/M-VNP-Q2 must supply their
own computed mode-basis samples explicitly before adopting this observer.

Generic offline bridge from solver-neutral sampled scalar-input radiation transfers to
real-tesseral SH filters. Radiation solvers supply complex-SH training rows, direct
held-out fields, source semantics, and diagnostics; `fs-couple` does not depend on them.
It converts `exp(-i omega t)` to `fs-vfit`'s `s=+i omega`, fits prewarped abscissae, and
applies unprewarped Tustin. Direct disjoint held-out fields gate artifacts; stable
proper fits may retain caller-admitted constant feedthrough but always forbid `s*e`.
The runtime owns visible per-filter state, superposes in fixed input order, and refuses
whole samples transactionally. EstimateOnly covers the fixed reference before `1/r`;
moving-source, FW-H, near-field, feedback, room, head, and passivity remain unclaimed.

### Layer-specific material memory in the reciprocal string bank

`with_material_layered_string_bending_relaxation` selects an independent
equilibrium-modulus role and ordered Maxwell-arm property pairs for every
original core/coating layer. Layers may have different names and arm counts;
an explicitly empty spectrum represents an elastic layer without invented loss
data. The existing uniform-selector API lowers through this same path.
`fs-material::section::relaxation` now owns source admission and section
integration instead of duplicating them in the acoustic coupling layer.

The bank retains its existing mass-normalized curvature projection, per-arm
history, reciprocal body/bridge loading and exact propagation. Explicit
Kelvin-Voigt viscosity remains a separate parallel dashpot; supplied Maxwell
increments must describe additional relaxing stiffness. Source frequency bands,
memory-state budgets, cold binding and sample-clock admission still apply.
Source images expose normalized `layered_selections()` plus the original
sections, layer spectra and channels. Focused G3 regressions cover bitwise
uniform/explicit layer equality and a two-arm coating/elastic-core substitution
through actual pressure and sympathetic body response, retaining the original
energy-balance assertions. These are numerical coupling checks, not physical
instrument qualification. Helical winding, axial/torsional memory, evolving
temperature and finite-strain material history remain outside this binding.

### Material string bank rotational bridges

`bowed_string::runtime::reciprocal::bank` admits explicit small-angle
`StringBodyRotationConnection` laws between a sourced torsional string mode
and a body rotation. The same polar inertia and GJ/shear-viscosity sources that
drive bow surface motion supply the string angle map. A body rotation can come
from two retained bridge footprints in the original material plate basis and
their physical separation; this is a secant rotation, exact for rigid rocking,
with a spatial approximation for a deforming plate. It is not an independently
tuned body mode. Raw angular modal maps retain their producer's authority.

The existing coupled-modal owner receives the exact work-conjugate conversion
`x = r theta`, `k_x = k_theta/r²`, `c_x = c_theta/r²`, using the sourced outer
string radius. Reciprocal torque, rotational spring storage and dashpot loss
therefore enter the same solve as bending, all other strings and the shared
body/cavity. Its force limit bounds `|torque|/r`. Body acceleration supplies the
existing pressure observer after bow release; twist receives no direct audio
gain. The new binding must precede playback and string motion. Original source,
body, gas, pressure and interface bindings survive multiple-bow scheduling.
G2 checks an independent two-oscillator angular-spring solution; G3 exercises
torque-only release, shear and body-density substitutions, original spatial
footprints and joint bowing; G0 covers rotational preload energy and refusals.
This does not establish nonlinear bridge contact, moving endpoints, local
plate-slope convergence, large-angle torsion, exterior radiation loading, or
measured violin/cello/bass fidelity.

### Material string bank independently directed translational bridges

`bowed_string::runtime::reciprocal::bank::polarization::StringBodyTranslationConnection`
and `ReciprocalStringBank::with_translational_body_connections` add independent
centreline translation laws to the existing material string bank. Each law
retains an original string index, interior station, signed unit
`TransverseDirection`, original-body `ModalAttachment`, nonnegative stiffness
and viscous resistance, and signed rest extension. The body map must measure
displacement along the SAME physical signed direction as the string map.
Body indices address the supplied bodies before strings are prefixed.

The relative coordinate is `direction dot y_string - x_body - rest`.
The original material source supplies the string projection with its actual
gas-loaded transverse mass. Its transpose supplies reciprocal force in the
same `CoupledModalSystem`; spring storage and dashpot loss enter the existing
work ledger. Original scalar bridges remain active, so distinct directions
and body maps can load both transverse planes without duplicating the
material string or prescribing modal detuning. Centreline translation has
zero direct torsional projection. Neither a complete rigid bridge nor its
contact law is inferred.

The retained source exposes `translational_body_connections()`. Cold gas,
cavity and rotational-bridge rebuilds retain those declarations and lower
their string maps against the current sourced realization. Original body
maps, material memory, gas state, exterior radiation, pressure observation
and dry bow interfaces stay with their existing owners. Binding accepts
one complete nonempty set before motion or playback and before nonlinear
storage; an empty set is a no-op. A one-plane source refuses any nonzero
second direction component. Scalar omitted-mode residual bridge flexibility
refuses this composition. Invalid stations, body indices, nonfinite or zero
body maps, nonpassive laws and owner budget violations refuse without
publishing a replacement system.

The `material_strings` frontend accepts an optional
`translational-bridges COUNT` block after thermal-support and rotational-bridge
records, before cylinder gas loading and `cavity`:

`translational-bridge STRING STATION SX SY SHELL_COMPONENT NODE GX GY GZ K C REST`

It requires an original material shell or assembled-shell body. The adapter
dots a global unit direction with the selected physical node's XYZ
translation in every retained `physical_modes` vector, preserving the original
mass normalization and eigenvector signs. Assembly components share ONE
assembled modal basis; their physical component index selects geometry,
while the modal attachment uses body zero. Recompiling after a body-material
substitution regenerates this projection from the new material solve.

String-local and body-global directions declare an orientation mapping.
The input has no string-to-world frame, so consistency between those
directions, including across connections, remains the geometry author's
responsibility. Unit-vector admission does not verify that mapping. The
frontend admits at most 16 additional rows and requires a retained second
plane whenever `SY` is nonzero.

Core regressions cover transverse-basis covariance, reciprocal second-plane
pressure, material substitution, zero direct torque, energy closure, cold
admission and gas/rotational attachment order. Frontend regressions project
the actual shell node, retain gas/cavity bindings, compare callback sizes
and change the body material. `directional-bridge.material-strings` is a
synthetic integration fixture, not measured instrument or bridge
qualification. Full three-dimensional unilateral bridge contact, moving
geometry and measured violin/cello/bass fidelity remain separate requirements.

### `string_specimen`

`acoustic_realize::LinearMaterialStringRuntime` incrementally advances an explicitly
linearized pinned string from a `ResolvedStringSpecimen`, using the existing
exact-ZOH modes, loss admission and compact strip acceleration/jerk observer.
Material/operator state, modal displacement/velocity, observer history and sample
position publish at one epoch after all admission checks and a final `Cx` poll.
`step` returns a `MaterialStringFrame` pairing that epoch and the acoustic frame
with `StringDissipation`: material Kelvin-Voigt heat, oscillating-cylinder air
loss, and authored damping without an identified thermal destination. Coefficients
are retained when the modes are constructed and rebuilt with each rebind.
For constant additive viscous coefficients during a sample, each mechanism's
energy is its coefficient fraction times the mode-local viscous loss: all share
the same integral of squared velocity. This follows the viscous work law in
[MIT 8.01, section 23.5](https://ocw.mit.edu/courses/8-01sc-classical-mechanics-fall-2016/mit8_01scs22_chapter23.pdf).
The existing Rayleigh option still replaces the entire damping model; its loss
is authored, with no silently added air or material contribution. Tiny negative
mode-local differences and undamped drift remain in a signed roundoff residual,
not in heat. That residual must fit the combined modal/aggregate/reduction
roundoff allowance before the runtime commits. Diagnostic microphone pressure
adds no second radiation-energy debit. The returned partition does not itself
debit reservoirs or establish a closed thermal feedback/entropy model.
Refused steps and rebindings preserve the complete accepted runtime. Cloning an
accepted runtime retains the acceleration history needed for identical continuation.

`project_point_forces` maps prescribed transverse loads in newtons at material
stations to the current mass-normalized sine basis. `point_motion` uses the same
participation to return accepted displacement [m] and velocity [m/s]. The transpose
pair preserves power and held-load work; endpoint participation is exactly zero
for pinned supports. The resulting forces feed the existing mechanical and thermal
steps, including after refinement and material updates. Projections must be used
before the next epoch change. This does not supply contact laws, opposite-body
reactions, moving-load work or compliant-support dynamics. G1 checks an independent
physical oscillator and work; G3 checks multipoint power/work through thermal
updates; G4 checks endpoints, invalid inputs, overflow, cancellation and retry.

`initialize_point_load_equilibrium` connects the same physical force projection
to the modal owner's static solve `q_k = f_k / omega_k²`. It initializes only an
unexcited runtime before its first sample, establishes zero prior acceleration,
and advances the epoch without changing material state or the sample clock.
Holding those forces retains equilibrium; releasing them excites the same live
mechanics and pressure observer. Stored initial energy assumes a settled loading
history before the simulation window; it is not a simulated actuator-work or
contact receipt. Existing motion cannot be overwritten. G1 checks static compliance
and release response; G4 checks repeat, invalid-force and cancellation refusals.

`rebind_ambient` resolves a new prescribed gas through the existing gas owner
and rebuilds air drag and compact observer weights atomically. Material receipts,
modal coordinates, previous acceleration and accepted sample count are retained;
the epoch advances with zero parameter work because this model has no gas-added
mass or gas-dependent elastic storage. The thermal string exposes the same update
with expected-epoch checking and unchanged material enthalpy. Boundary heat must
still enter through thermal transport. This is a piecewise prescribed ambient
change, not an atmosphere transient, fluid-loading work or propagation-delay solve.
G1/G3 checks use cylinder drag and compact pressure formulas including retained
jerk history; G4 checks invalid gas, stale epochs and cancellation preserve retry.

Rebinding requires the exact same immutable card, length, linear density and
sine-mode count. This makes the mass-normalized basis identical, with no sorted-mode
matching, displacement reset or sound-level rescaling. Tension, EI, damping and
observer weights rebuild together; the source state remains inspectable by other
consumers. Signed parameter work is the change in modal elastic energy at fixed
displacement/momentum, checked against an explicit caller budget. It is returned
for the thermal/external parameter owner's opposite transfer, not booked as damping.
The axial prestress reservoir, heat transport and closed thermomechanical balance
are not supplied here. Coefficients are externally scheduled and piecewise constant.
The acoustic observer retains the preceding accepted acceleration across updates;
abrupt coefficient changes can produce a transient, without a crossfade.

`substitute_material` additionally accepts a different source card and line mass
at exactly the same length, circular radius and pinned sine-mode count. It uses
the same coefficient/loss/observer rebuild, with an explicit canonical transfer:
`q_new = sqrt(mu_new / mu_old) q_old` and
`p_new = sqrt(mu_old / mu_new) p_old`. Physical displacement and transverse
momentum remain continuous; velocity generally changes. The reported external
parameter work includes both kinetic and elastic energy changes and must fit the
caller's budget. Previous acceleration is rescaled into the new mass normalization
so the pressure observer retains its physical history. Physical forces must be
projected again after the epoch changes. The compact observer includes endpoint
acceleration changes but does not reconstruct distributional radiation from the
instantaneous velocity jump; that needs a time-resolved transition.
G1 checks an independent forced oscillator,
canonical momentum and pressure including the even-mode jerk term; G3 checks a
substitution round trip and identical-source continuation; G4 checks complete
rollback after geometry/basis, work-budget and cancellation refusals. This is a
prescribed material intervention, not chemical conversion, melting, mass transport
or heat deposited in either material. The thermal wrapper retains its stricter
same-card binding and does not expose this substitution.

Explicit `refine_modes` increases the retained sine basis at fixed geometry without
advancing time. Existing coordinates and observer acceleration history are retained;
added modes start at rest with zero previous acceleration. Shared admission enforces
mode, Nyquist and material-frequency limits before publication. The thermal wrapper
checks the expected epoch and preserves enthalpy. G1/G3 checks continuation energy,
pressure and independent forced response; G4 checks atomic refusals and cancellation.
This does not reconstruct omitted history or estimate truncation error, and does not
support coarsening or transfer between different geometries or eigenbases.

The rung omits transverse-motion-induced axial stretching, moving supports, extra
polarizations, Prony memory, basis/geometry transfer, phase change and radiation
backreaction. It uses prescribed piecewise ambient gas and fixed listener distance, with no atmospheric
path filter or common propagation delay. G1 uses a synthetic thermal-expansion card
and an independent damped-oscillator/pressure calculation. G3 cycles a same-card
density/temperature update at prescribed mass through changed cross-section, bending
stiffness and pressure without changing the retained sine basis. G4 exercises material-work
and late observer refusals, cancellation and exact suffix replay. This is focused
implementation evidence, not measured string fidelity or full MR05/MR34 completion.
The loss-partition G1 test independently integrates squared analytic modal
velocity and feeds only material heat into the existing equilibrium enthalpy
owner. Its synthetic constant-coefficient solid absorber checks the combined
mechanical/thermal/air balance and rejects duplicate heat, air credited to the
solid, and a missing mechanical term. It stays wholly on the solid branch;
thermoelastic feedback and phase dynamics are not tested. G3 changes loss laws
without resetting vibration and checks that the returned destinations change.
Existing G4 refusal/retry tests compare the complete frame, including loss fields.

`acoustic_realize::ThermalMaterialStringRuntime` owns that incremental string,
its immutable material card and a uniform equilibrium enthalpy state. It closes
the dissipative temperature/viscosity feedback. Its ordinary constructor admits
the prescribed-radius or prescribed-mass, fixed-tension, constant-density,
constant-Young's-modulus solid rung. Those
restrictions keep elastic and axial prestress storage independent of temperature;
there is no omitted stiffness-change or expansion work. Density must agree with
the solid branch of the thermal chart. The initial absolute temperature and card
must agree across both owners. Original query-axis schemas and selected claims
are retained; subsequent queries pin those claims, with no automatic substitution.
Constant properties and linearly interpolated temperature-only laws are admitted;
exact-only temperature tables do not authorize an evolving path between samples.
Viscosity must be nonnegative throughout its supplied source curve, including
interior knots whose values cannot be inferred from positive endpoints.
For `FixedMass`, the prescribed mass and `m/L` remain exact through thermal
updates; constant density and loaded length keep the resolved radius and basis
unchanged. Thermal area and volume use `ResolvedStringSpecimen::radius_m()` and
its cross-section, rather than reconstructing mass through a rounded radius.
This ordinary constructor does not admit changing density, expansion or geometric
state transfer. The source-bound thermoelastic constructor below also keeps mass,
section and loaded length fixed while evolving its thermal eigenstrain.

`try_new_thermoelastic` admits the original sourced `FixedThermalExtension`, a
scalar expansion coefficient and scalar or globally affine temperature-dependent
Young's modulus. Every original modulus knot must belong to that same affine law;
derivative jumps require a different caloric model. The original material and
expansion owners resolve the complete temperature path, retaining source claims,
reference temperature and reference stress-free length. At each accepted sample,
thermal expansion changes the actual tensile force and modal operator without
resetting displacement, velocity, observer history or the accepted sample count.
Slack, changing density/geometry, melting and unsupported source paths refuse.

The axial Helmholtz storage is `Psi_a = A L_ref E(T) e(T)^2 / 2`, with
`e = (L-L_ref)/L_ref - integral(alpha dT)`. The retained native transverse
storage is `Psi_v = sum(omega_n(T)^2 q_n^2 / 2)`. Mechanical entropy is
`s = -partial_T(Psi_a + Psi_v)` at fixed supports and modal coordinates, and
mechanical internal energy is `K + Psi_a + Psi_v + T*s`. The public
`mechanical_internal_energy_j` reports this same storage for coupled consumers.
Positive total fixed-coordinate heat capacity is checked on every original
caloric segment, including the amplitude-dependent modal term; the exact
quadratic extrema of `T*Psi_TT` are checked for the admitted affine law.

Mechanics advances once at the accepted temperature. At its resulting modal
coordinates the existing bounded scalar feedback solver closes reversible heat
`Q_rev = -delta(Psi_a) - modal_parameter_work - delta(T*s)` through the original
caloric/transport owner. The heat includes entropy change due to modal motion,
not just the endpoint temperature change. Thermal callbacks are repeatable
`FnMut` proposals from the same accepted state; only the final report/reservoir
candidate survives. The ordinary constant-stiffness path still calls transport
once. The final frame exposes axial, modal-parameter and entropy-energy changes
separately, updates the actual modal energies and pressure, and checks the full
first law. Material dissipation remains distinct from signed reversible heat;
neither numerical loss nor diagnostic acoustic pressure is deposited as heat.
The electrical and finite-gas paths share this owner, including reciprocal
conduction at the actual solved solid and gas endpoint temperatures.

This remains a first-order partitioned small-strain model with a fixed section,
uniform temperature and diagnostic compact radiation. Energy closure does not
prove exact entropy conservation, timestep convergence or instrument fidelity.
`tests/string_specimen/thermoelastic.rs` defines synthetic live-pitch, independent
internal-energy, electrical-source and refusal/retry checks. Native compilation
and execution were unavailable in the implementation environment; syntax and
energy-algebra checks are not reported as native test passes.

Each sample advances mechanics at its initial viscosity, credits material loss
and explicit external heat to total enthalpy, then resolves the new material
temperature and rebuilds damping. Final pressure and retained observer history
use the updated acceleration. This is a first-order partitioned coupling; exact
mechanical substeps do not make the coupled solution exact. A positive caller
temperature-increment budget constrains every accepted step, and temporal
refinement must assess coupling error. Air loss and mechanical roundoff remain
outside material heat. The returned first-law residual includes work, both storage
changes, outgoing losses and numerical loss, with a scale-aware arithmetic
allowance including the enthalpy reference. It is not a physical-error certificate.
All candidates publish at one epoch after the final cancellation poll. Reusing
an accepted input epoch refuses before consuming another heat input or sample;
thermal/material refusal after the mechanical substep leaves both owners intact.
Cloning preserves exact continuation, including thermal state and source pins.
The direct `step` takes external heat in joules. `step_with_current` instead
resolves the existing `fs-material::conductor` from the specimen's pinned
electrical-resistivity claim, current temperature and actual length/area. An
ideal source deposits `I² R(T_initial) dt` once in the external heat account;
material vibration loss remains separate. The returned `ElectrothermalStringFrame`
retains initial/final sourced resistance and the same accepted thermal/mechanical
epoch. Final electrical resolution precedes publication; any refusal preserves
the accepted state and retry position. This is first-order partitioned uniform
Ohmic heating with insulated boundaries, not circuit dynamics, skin/contact
resistance, a charge-field solve, Lorentz forces or motional emf. Synthetic
G1/G3/G4 tests check analytical temperature refinement, current-sign symmetry,
resistance/damping/pressure response and atomic retry, not experimental fidelity.
`step_with_current_and_thermal_transport` composes the same current source
with the existing transport callback. The thermal owner receives vibration
plus Joule heat as internal deposition and returns actual boundary heat.
`ElectrothermalStringFrame<R>` retains Joule heat, boundary heat and the
concrete transport report separately; its coupled external-energy account is
their sum, so electrical work is neither lost nor counted twice. Both paths
share current admission, and final electrical resolution still precedes
publication. Zero current reduces to the ordinary thermal transport trajectory.
The independently prescribed fluid/enclosure reservoirs retain the existing
uniform-body, Biot and infinite-reservoir restrictions.
`step_with_voltage` and `step_with_voltage_and_thermal_transport` share these
same transactions under `OhmicDrive::Voltage`. Each sample freezes the current
`V/R(T_initial)` and deposits its Joule energy; the next sample uses the updated
resistance. The frame retains the authored drive, initial sample current and
current demanded at the final accepted resistance. Final electrical evaluation
must succeed before publication. Positive temperature coefficients therefore
reduce fixed-voltage heating while increasing fixed-current heating. The
partition remains first order; no within-sample circuit solve is claimed.
`step_with_thermal_transport`
accepts a thermal proposal through the shared `fs-material::phase` input/result
carriers. Passing `|input| environment.advance(cx, input)` with
`fs-conduction::lumped::LumpedThermalEnvironment` uses the existing
`solve_lumped_enthalpy` for one
backward-Euler step over the same audio sample. It derives the exposed cylinder
area `2 pi r (L+r)`, volume and Biot length `V/A` from the resolved specimen,
with no independently authored thermal mass or geometry. Convection and surface
radiation act on this same complete boundary; supports are thermally insulating.
The environment is a prescribed infinite, uniform reservoir and enclosing
radiation field with independent fluid and enclosure temperatures. Opposing
convection and radiation fluxes contribute their signed heat once. It neither changes nor debits the separately
declared acoustic gas. Transport may be declared or card-backed; the conduction
owner enforces its existing card identity and temperature-coverage rules.
The callback must not mutate or debit an external source before acceptance.
`step_with_thermal_state` also accepts caller-owned state by exclusive borrow.
Its callback reads that accepted state and returns an owned candidate alongside
the normal thermal proposal. Only after all coupled checks succeed does the
method replace the caller's state. A finite support therefore advances with
the string without a separate caller commit. G3/G4 finite-support integration
checks heat debit, temperature/damping/pressure feedback, late foreign-chart
refusal and exact retry of both owners. The transport owner validates its
opaque state; callbacks must not mutate shared interior state or external
resources. This is synchronous publication through exclusive borrows, not a
cross-thread or durable transaction.
`step_with_drive_and_thermal_state` provides the same publication boundary
for either `OhmicDrive::Current` or `OhmicDrive::Voltage`, reusing the existing
electrothermal path. External thermal state remains a candidate until final
resistance/current resolution also succeeds. The callback receives Joule plus
material heat; its boundary transfer excludes both. Finite-support tests check
this separation, current/voltage callback parity, cancellation and retry.
The coupling owner verifies that the proposal uses its original phase chart
and finite heat/residual budgets, then checks the coupled energy balance. This
keeps fs-conduction out of fs-couple's production dependencies: FEEC's optional
topology feature already depends on fs-couple. The cross-crate test has a
development-only conduction dependency.

Material loss enters once as mean internal power for that sample. Boundary heat
uses the actual endpoint convection/radiation powers, not a subtraction from
the computed enthalpy change that could hide thermal solve error. The returned
transport trajectory retains those powers, its Biot value and energy residual.
The caller's thermal energy-residual budget is checked separately and appears
as `thermal_solve_tolerance_j`, zero on the direct-heat path. Coupled first-law
acceptance uses that budget plus the reported arithmetic allowance; neither
certifies time-discretization error. A late thermal, phase or material refusal
preserves mechanics, temperature, source receipts, epoch and observer history.

The fixed corpus Biot ceiling is applied, including the conservative radiative
coefficient over the chart. This is a necessary uniform-body check, not proof
of initially uniform temperature or homogenized internal heating. Rapid local
sources and nonuniform support transfer require spatial transport regardless
of step size. The analytical convection reference and lumped assumptions are
described in [MIT's transient heat-transfer notes](https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node129.html).
The thermal chart retains its existing source authority. Finite-support contact
uses the conduction adapter's fixed resistance and uniform-body limits; spatial
conduction, entropy production certification, thermal
expansion, changing stiffness/density, and solid-liquid transfer are not claimed.
G1 compares the coupled trajectory with an independently integrated continuous
ODE and checks temporal convergence. Synthetic heat capacity and viscosity make
feedback observable in a bounded test; they are not measured material data.
G1 additionally checks analytical heating/cooling and first-order refinement
for two radii, plus separate and combined convection/radiation inputs changing
the material damping and diagnostic pressure. G4 tests Biot, thermal residual,
source-domain, phase and cancellation refusals, duplicate epochs, retry and
cloned continuation. These exercise the existing solver directly, not a mock.

`with_uniform_circular_material_state` lowers one homogeneous isotropic circular
specimen to the existing `PrestressedString` realizer. It consumes the shared
`ResolvedMaterialStatePoint` API, requiring only positive SI density and Young's
modulus. Area, area second moment, linear density, mass, EA, EI, and diameter
follow from the same radius and length. All source receipts and the complete
material state remain available in `ResolvedStringSpecimen`; its specimen hash
binds geometry and material resolution, excluding mechanical loading and solver
options. It is not the complete scenario identity.

The original binding holds the template's geometry and tension fixed.
`with_uniform_circular_material_and_prestress` additionally takes an explicit
`StringPrestress`: fixed tension, fixed extension from a stress-free length, or
target primary undamped fundamental frequency. Fixed extension computes `T=EA*strain`
and enforces the caller's declared linear strain limit; that assumption is not a
measured yield boundary. Target frequency inverts the existing pinned-pinned
Euler–Bernoulli dispersion relation, including EI. Slack/compressive solutions and
extension/target prescriptions on the moving-end basis refuse. The prescription
remains separate from the material receipts, available through `prestress()`.

`with_uniform_circular_material_and_constraints` also accepts a
`StringGeometryConstraint`: fixed radius or fixed total mass at the template's
current loaded length. Fixed mass derives `r=sqrt(m/(rho L pi))`, retains the
prescribed mass, and sets `mu=m/L` directly. It never reconstructs either mass
quantity through the rounded radius. Thus changing density at fixed mass/length
keeps exactly the same numerical modal basis while changing EA, EI and observer
width. The radius-derived area and `rho A` agree with these mass constraints to
floating-point roundoff. Geometry and prestress constraints are independent and both remain
available on the result, separate from material authority. The specimen hash
binds resolved geometry, mass and linear density, not the authored comparison
policy. The v2 identity distinguishes prescribed masses even when rounding maps
them to the same radius.

`with_uniform_circular_thermal_extension` consumes the existing
`IntegratedIsotropicThermalExpansion` and current elastic resolution from the
same exact card, query coordinates and axis conventions. It computes the
small-strain eigenstrain relation `T = EA * ((L-L_ref)/L_ref - integral(alpha dT))`,
using one reference configuration, as in
[MIT 16.20, section 3.7](https://web.mit.edu/16.20/homepage/3_Constitutive/Constitutive_files/module_3_no_solutions.pdf).
The current geometry/mass constraint remains explicit; both elastic and expansion
receipts survive in the result. `FixedThermalExtension` records the mechanical
prescription, while `thermal_expansion()` distinguishes the sourced path from
raw authored eigenstrain. A subsequent ordinary binding requires the thermal
path to be supplied again; it cannot silently retain old thermal evidence.
Total, thermal and positive elastic strains must stay within the authored limit.
Slack, compression and domain exit refuse; no buckling solver is implied.
Fixed tension remains a separate force-controlled prescription. G1/G3 check a
piecewise-linear expansion curve across temperature knots, heating/cooling force
changes, mass and receipt retention, and pitch from the real pressure solver.
G0 checks current-state/source mismatch, invalid reference length, unsupported
moving supports and strain-domain refusals. This initializes independent uniform
specimens; it does not evolve temperature, exchange heat, transport strain or
solve lateral contraction during playback.

All policies preserve the selected loss law, supports, and modal budget. These compare
specimens at resolved states; thermal evolution, transverse contraction,
yielding, anisotropy, or melting are not solved. Target pitch describes an undamped
linear eigenfrequency, not a damped or finite-amplitude pressure frequency. Diameter
supplies the existing compact strip-radiator width, without a new claim of accurate
circular-wire radiation. The density and modulus requirements use exact
dimension-only quantity schemas; equal-dimension semantic kinds or waveform
forms cannot be erased when lowering to the solver. Admission refuses missing
or quantity-schema-mismatched properties,
invalid mechanical inputs, and nonfinite/zero derived products before publication.
The finite scalar calculation is deterministic; no inner-loop database access is
added. `tests/string_specimen.rs` checks G0 admission, G1 mass,
EA/EI and beam modes, and G3 geometric/material scaling plus pressure-measured pitch
through the actual acoustic realizer. G1/G3 cases also cover extension-derived
tension and material changes under explicit extension/target-frequency constraints,
plus fixed-mass diameter/stiffness changes and their bending-sensitive pressure pitch.
Synthetic data check numerical behavior, not physical validation of a measured material.

`ResolvedStringSpecimen::with_kelvin_voigt_bending_loss` selects a constant
Kelvin–Voigt bending law using the retained `kelvin_voigt_bending_viscosity`
property [Pa s]. Density, Young's modulus and viscosity must be dimension-only
validity-wide scalar constants or temperature-only curves at fixed positive
absolute T. Frequency samples are not promoted to a broadband law. The
viscosity must declare a frequency domain, intersected with any
frequency restrictions on density and modulus. Other state coordinates stay at
their admitted query point. The descriptor carries `eta I`, the intersected band,
and the material bundle identity. Rebinding an already viscous specimen recomputes
all three from the new geometry and material. No source uncertainty is improved.

All three material loss binders (Kelvin–Voigt string/plate and Prony string)
use `ResolvedScalarProperty::angular_frequency_band_rad_s`. One static typed
frequency coordinate can have any name and be expressed in Hz or rad/s;
legacy untyped `omega` remains rad/s. The normalized band intersects every
participating coefficient's applicability before the existing modal checks.
Multiple frequency coordinates refuse because they may be independent physical
variables. No source receipts or interpolation decisions are rewritten. G3
tests compare actual pressure bit-for-bit across the two unit conventions and
legacy inputs; narrowing only density support rejects an out-of-band string.

The shared linear/Kirchhoff–Carrier/moving-end modal path projects
`mu y_tt - T y_xx + EI y_xxxx + eta I y_xxxxt = 0`, giving
`zeta_k = eta I k^4/(2 mu omega_k)`. Thus tension contributes stiffness but no
viscous loss. This is the constant-coefficient Kelvin–Voigt term in
[Sakthivel et al. (2023)](https://arxiv.org/abs/2301.07931), with tension added
and projected onto the existing sine/cosine bases. The explicit law replaces
the authored constant modal ratio; simultaneous Rayleigh
or nonzero authored damping refuses. Air resistance is still added separately.
Every retained reference frequency, including the secondary polarization, must
lie in the descriptor band. G1 checks `eta I` and rebinding; G3 checks actual
pressure decay against an independent diameter-form expression in all three
solver paths, including zero-viscosity air-only decay. G0 covers unavailable
data, wrong quantities, sampled frequency curves, duplicate losses and out-of-band modes.

This Kelvin–Voigt choice models linear bending strain-rate loss only. It does not supply axial
viscosity, material memory, amplitude-dependent loss, time-evolving material
state, thermoelastic/joint/radiation dissipation, or measured material validation.
The modal applicability check concerns reference frequencies; nonlinear
frequency shifts and generated harmonics are not an in-band certificate. The
existing moving-end cosine approximation and compact observer remain unchanged.

`with_prony_bending_loss` binds a causal generalized Maxwell bending spectrum
from explicit modulus [Pa] / relaxation-time [s] property pairs and
`equilibrium_young_modulus` [Pa]. `with_standard_linear_solid_bending_loss`
selects the single `relaxing_bending_modulus` / `bending_relaxation_time` pair
through that same binding. The explicit equilibrium modulus must equal the
Young's modulus used to derive EA/EI; a sampled dynamic modulus cannot silently
be treated as an equilibrium value. All selected coefficients, including density,
must be dimension-only validity-wide constants or temperature-only curves at
fixed positive absolute T. Every relaxation-time claim must
declare a frequency band, intersected with all other coefficient bands. An empty
spectrum is explicitly elastic and requires an equilibrium-modulus band. Geometry
or material rebinding resolves every retained source pair again, recomputing
each delta EI and replacing the source identity and time constants. Missing or
mismatched selectors refuse rebinding; duplicate source pairs refuse binding.
These are resolved inputs, not measured-material validation.

Each mass-normalized modal coordinate receives generic `fs-phs` relaxation arms
with `a_j = delta EI_j k^4/mu`: `q'' + c q' + omega0^2 q + sum_j a_j(q-v_j)=force`,
`v_j'=(q-v_j)/tau_j`. This is the internal-strain generalized Maxwell formulation
described in [Bleyer's viscoelasticity tutorial](https://bleyerj.github.io/comet-fenicsx/tours/nonlinear_problems/linear_viscoelasticity_jax/linear_viscoelasticity_jax.html),
projected onto the existing beam modes. Nonlinear axial storage, contact ports,
body joins and pressure observation use the same augmented pHS state. A held
pluck initializes every `v_j=q`, representing full relaxation before release.
Positive-strength branches add memory states in branch-major, mode-major order;
zero-strength branches retain source selection but add no state or active pole.
The owned spectrum is borrowed throughout the audio loop. These
material arms produce storage dispersion and dissipation together; no PCM decay
envelope or substitute modal damping ratio implements it.

Admission checks both relaxed and instantaneous reference frequencies against
the material band and Nyquist guard, using the sum of all branch stiffnesses
and including fixed-end detuned polarizations.
Moving-end secondary polarization refuses on this path. Kelvin–Voigt, Rayleigh
and nonzero authored internal loss cannot be combined with the Prony descriptor;
air drag remains separate at the relaxed reference frequency. G0/G1 tests cover
sample-clock admission (`dt/tau_j <= 0.1` for every nonzero arm and
`dt*omega_instant <= 0.1`, otherwise increase sample rate),
source/quantity/band refusals, geometry rebinding and actual pressure pitch/decay
against independent characteristic-polynomial roots in fixed linear, nonlinear
and moving-end realizations. Multi-arm G1 pressure is compared with independently
integrated hereditary-strain dynamics and checked for second-order refinement;
G3 checks equivalent branch splitting/permutation, zero-arm omission, detuned
polarizations, and body/contact composition in fixed and moving-end paths.
This remains linear bending memory at a fixed state. No measured spectrum fit,
axial creep, evolving temperature, phase change, arbitrary pre-release loading,
or in-band certificate for nonlinear harmonics is supplied.

### `driving_point::regional::field`

`RegionalDuct::with_exterior_gas` now binds side-hole radiation as well as the
terminal load and observer. Every regional run retains its own chimney gas
and receives the same explicit exterior through `fs_duct::exterior`; pressure
and volume flow join the two media without a second end correction. Open and
vented holes may therefore contain a different mixture from the atmosphere.
This also reaches the continuous-field response/reflection, baked terminal
loads and scheduled reed renderer. A closed main outlet does not disable
radiation through open side holes. An omitted exterior retains the original
local-gas convention and identical media retain their original bits. Static
pressure mismatch and unavailable varying-medium mutual-hole laws still
refuse. Focused G1/G3 tests cover physical side-hole transfers, upstream gas
regions, continuous fields and changes to both reed motion and microphone
pressure. This is frozen pressure-matched acoustics, not a mixing/mean-flow
model or experimental validation of a complete instrument.

`RegionalDuct::from_gas_field` discretizes the shared `fs_material::gas_field::AxialGasField`
at explicitly refinable midpoint cells (at most 64 total), retaining the
existing cylinder/cone TMM and pressure/volume-flow interconnection. It adds no
propagation or loss law. The exact x=0 state feeds source density and the inlet
characteristic impedance, and x=L feeds the existing terminal radiation load;
an independent exterior can still be supplied explicitly. Uniform fields keep
the original complete geometry and arithmetic. Nonuniform tone holes and smooth
radius profiles refuse pending the corresponding continuous-field duct support.
`ReedBoreVoice::new_with_gas_field` feeds this same regional line into the existing
reed runtime and compact-jet observer. No instrument constants, tuning law,
evolving breath or calibrated exterior observation are added.

The focused field tests compare uniform responses bitwise, exercise cone
radiation and endpoint states, and refine a supplied 1 K axial temperature
gradient over 4/8/16/32/64 cells. Fundamental resonance convergence is second
order; the third-to-first register interval is checked against independent
energy-weighted perturbation. A rendered-pressure test checks uniform legacy
identity, exact source density and the gradient's nonzero audio effect. These
are numerical tests of prescribed fields, not measured instrument qualification;
native results are recorded separately.

`RegionalDuct::from_continuous_gas_field` uses the existing fourth-order
Magnus operator for simultaneous smooth geometry and primitive gas fields.
Its optional `WallPin` now reaches the same local wall-admittance owner as
uniform cylinder/cone ducts at every gas/geometry quadrature point. Responses,
characteristic-line feedback and baffled outlet pressure retain that common
wall binding. The wall's explicitly supplied per-area mass/stiffness/resistance
remain constant; actual radius and slope set the surface coupling. This is a
locally reacting wall, not an axial shell, evolving wall temperature or a
nonlocal wall-radiation model. The continuous binding now admits tone holes
at their exact local gas state through the original homogeneous-chimney and
junction laws. Closed holes and a single open/vented hole may lie in a
gradient. Multiple open/vented holes retain every mutual correction and
require a constant primitive field throughout their connecting span;
unsupported varying inter-hole media still refuse. This also reaches the
existing scheduled reed renderer without a separate propagation path.
Open/vented chimney flows now cross this continuous-field binding into the
same bell-and-hole microphone realization. Axial delays use the original
gas-field quadrature and each chimney's local sound speed. Feedback, chimney
mouth loads and exterior microphones use one homogeneous exterior: the explicit
binding, or the actual outlet gas when omitted. The same choice applies to
regional interfaces, arbitrary terminal loads and characteristic-line shortcuts. Selecting a continuous field cannot
silently turn an open hole into a non-radiating source. Closed holes remain
loads only. This retains the existing sampled-transfer and infinite-baffle
approximations, rather than adding a new radiation or transport law.
A focused impulse regression covers the combined wall/field changes to
both incoming feedback and causal external microphone pressure.

### `acoustic_realize` / `pcm_wav`

`realize_assembly_with_gas_and_absorption` accepts a caller-resolved gas and
an explicit `GasAbsorptionModel`. Its observer always uses the selected law;
invalid ISO-air conditions refuse without a classical fallback. `RealizeOptions`
carries the same resolved gas and observer law alongside the bore discretization.
Ordinary assemblies, regional plate charts, `CompiledMaterialAssembly` and
`CompiledValveMaterialAssembly` all forward those options to one dispatcher.
Material compilation therefore preserves the caller's explicit absorption and
reed-bore resolution, including their admission failures. The gas is validated
once before stepping and supplies excitation, duct coefficients, body/cavity
loading and observed pressure; solid card temperature and prestress stay
independent. Missing absorption selection on a supplied gas chooses classical
Stokes-Kirchhoff by input semantics, never by coefficient equality with a preset
or an invalid humidity sentinel. The scenario-air entry point preserves its
legacy absorption arithmetic. `RayleighMedium::from`
copies density/sound-speed bits without reevaluating an atmosphere. G3 tests
cover explicit ISO/classical selection through source-bound uniform/regional
solid bodies and material valves, pressure/composition changes reaching actual
pressure histories, unchanged solid identities, retained bore refusals, invalid
ISO humidity, exact medium views and legacy history equality. Ideal-gas pressure
changes density/loading without an invented sound-speed or pitch retuning.
GasState itself carries no source receipts; arbitrary-mixture molecular
relaxation, changing composition during a render and experimental validation
remain outside these frozen-state paths.

File-driven bowed, reed and plate-valve performances expose
`from_bytes_with_gas_sources` and retain their authored `GasRecord`s. Ambient,
regional and exterior declarations use the same L3 reader. The CLI loads
bounded source artifacts itself: repeated `--gas-component SPECIES NASA9
TRANSPORT-CLAIMS`, explicit `--gas-water none|ID` and `--gas-phase
frozen|liquid-water`. A primitive reed render additionally supplies
`--gas-record FILE`; wind/ensemble gas records live in their performance files.
No thermochemical route falls back to a preset. Source-aware sidecars retain
composition, temperature, pressure and all three pack identities per species.
Omitting the new options preserves legacy numerical and sidecar paths.

The native `percussion` example uses this same bounded source loader for
drum, cymbal-profile, explicit shell-mesh and paired hi-hat `ambient` rows,
including both shell files referenced by a pair and `--gas-record` overrides.
The selected `Ambient` owns the original `GasRecord`; its exact admitted
density and sound speed supply exterior radiation and distributed cavities,
and `rho*c^2` supplies the existing drum volume spring. Composition, phase
admission and species/NASA-9/transport identities remain with the run and its
ambient receipt. The same registry resolves a hi-hat squeeze film's separate
`gas` row; source viscosity and p/T/R enter the original compressible film
owner, which retains the complete record and its existing isothermal storage.
No new pressure force, damping term, time loop or solid retuning is added.
Missing, ambiguous and out-of-domain sources refuse, while legacy declarations
preserve their numerical state. Focused native tests cover actual source bytes
and nested file routes, NASA-9-dependent drum pressure feedback with exact
replay, and source species mass in compressed real-skin gas pockets. This
extends M0.1-W1 input-to-consumer coverage; it does not qualify a full-band
percussion instrument or add evolving gas temperature/composition.

The same native percussion route accepts `head_material` in a drum specification:
side, thickness, installed tension, explicit loss token and the shared
`MaterialRecord` (`PACK STATE [AXIS=VALUE ...]`). The existing native 8 MiB card
loader is shared with plate inputs; each named pack is loaded once per drum.
`material_source::plate_section` resolves isotropic E, nu and density at the
same physical query before the original `TensionedDisk` assembly. The full pack,
source record and resolved state remain with the executing drum. Numerical
coefficient overwrites and missing, ambiguous, foreign-state or out-of-domain
properties refuse. Only static properties at fixed external conditions are
admitted; a frequency/strain/rate sample cannot masquerade as an elastic law.
Geometry and installed tension remain explicit, with no inferred thermal
expansion or retuning. Existing head relaxation may compose only through its
unchanged equilibrium-modulus and loss-ownership checks. Focused regressions
cover exact numeric-source parity, independent mass/bending scaling, actual
contact/nonlinear-head/cavity response and measured-pack file loading. This
extends M0.5-W1 without claiming measured drum-film damping or instrument
qualification.

`realize_assembly` turns an `fs_scenario::AcousticAssembly` into observer
pressure by composing nameless primitives: a prestressed Euler–Bernoulli
waveguide (`fs-nlmodal::prestressed_beam_omega` / Kirchhoff–Carrier), a
`bernoulli_aperture`, a `stribeck_friction` contact, a
`traveling_wave_line`, `ModalAcousticTime`, and compact modal radiators.
There is no instrument crate and no instrument algebra. A guitar or
clarinet is one filling of those objects.

- Ambient state comes from `GasState::try_new_moist_air(T, p, RH)` once,
  before realization. Excitation, duct impedance/resonance, radiation loads,
  and observer absorption share that state. RH=0 uses the original dry-air
  constructor bit for bit; RH>0 requires 253.15–323.15 K and water mole
  fraction at most 0.15, with model-domain failures returned as `Ambient`
  errors. Invalid RH refuses admission. This is an estimated ideal mixture:
  viscosity and conductivity retain the dry-air fits, the saturation
  enhancement factor is omitted, and no exact humid-medium prediction or
  experimental validation is claimed. G0/G3 assembly tests check the dry
  boundary, refusals, independent mixture density/speed/impedance, deterministic
  pressure, and waveform-derived duct pitch. The fixed eight-cell Cauer duct
  retains its first-order spatial boundary error; the pitch check uses its
  analytic LC frequency scale and does not establish 5% continuum accuracy.
  ISO absorption remains the sole
  bulk absorption inside its window; no extra classical term is added.
- String path, `EA = 0`: triangular-pluck (or bow) modal ICs marched by
  `ModalAcousticTimeModel`. The shared fixed-end string observer reads actual
  acceleration: odd modes give `p = rho A_odd q_ddot / (4 pi r)` in physical
  displacement coordinates. Even modes have zero monopole area and signed
  moment `M_even = -w L²/(n pi)`, giving `rho M_even q_dddot/(4 pi r c)`.
  This follows the spatial expansion of distributed retarded monopoles
  ([Wood, eq. 6](https://euphonics.org/11-7-1-compact-sound-sources-monopoles-dipoles-quadrupoles/)).
  Linear acceleration includes held forcing and viscous loss; nonlinear
  acceleration uses the admitted pHS storage and `(J-R) grad H + G u`, so
  stretching, conservative contact, and applied forces reach emitted pressure.
  Jerk is a backward difference of accepted accelerations: first-order at the
  reported endpoint, second-order at the interval midpoint for smooth motion.
  The initial acceleration is evaluated from the initial state and first held
  force; no zero-history startup spike or ideal release impulse is synthesized.
  The observer retains the linear path's 1e7 Pa pressure ceiling.
  This is an uncoupled compact strip approximation with the even-mode observer
  on the positive string axis. Common propagation delay, angular listener
  geometry, cylindrical-wire radiation, resolved impact impulses, and radiation
  back-reaction are not supplied by this observer. G1 tests check signed SI
  pressure against a damped oscillator, even-mode jerk and refinement, nonlinear
  stretching-force pressure, and silence under a held static load. Fletcher
  inharmonicity `ω_n = n ω_1 √(1+B n²)` when `EI > 0`. Stokes air
  drag from `GasState` plus the authored constant modal ratio. The reusable
  `air_path::oscillating_cylinder_air_resistance_per_length` computes the
  Desvages (2018) section 3.2.2 cylinder resistance approximation
  `R=2 pi mu (1+r sqrt(2 rho omega/mu))` [N s/m²]; both string paths use
  dimensionless `zeta=R/(2 mu_linear omega)`. Width is interpreted as circular
  diameter. The approximation retains a constant low-frequency continuation;
  added mass, end effects, confinement, rarefaction, turbulence and nonlinear
  fluid drag remain unclaimed. Explicit Rayleigh damping still replaces the
  complete loss sum when no material bending law is selected. The explicit
  Kelvin–Voigt path above instead combines only bending viscosity and air drag.
  Without a selected constitutive law, `damping_ratio` is a finite nonnegative
  constant per-mode reduced-model input. Zero explicitly assumes no internal
  loss; elastic EI does not create dissipation. No hidden `2e-7*omega` addend
  or relaxation spectrum inferred from a single ratio remains. The authored
  ratio is not a material measurement; physical loss requires admitted data
  and a selected Kelvin–Voigt or causal Prony law. G1 checks the third mode's
  actual SI pressure against the independently damped oscillator, and G3 checks
  equivalence to explicitly zero bending viscosity across linear, KC and Dirac
  paths. Invalid authored ratios refuse instead of being clamped.
  G1/G3 tests check the resistance formula and
  pressure ring-down under changed gas pressure. A second
  polarization at `1+detune` is a second member on the same clock
  and shares every plate (it is not an independent unused body).
- String path, `EA > 0`: `fs-nlmodal::kirchhoff_carrier_string`
  (fixed-fixed) or `kirchhoff_carrier_moving_end` (Dirac join)
  + Gonzalez
  pHS step, with the same stiff-string frequencies overwriting the
  linear ω_n. Bridge force `T_eff y'(0)` drives every listed body mode
  (Carcagno-band Sitka pair is a constructor, not a named guitar).
- Bow: a set-valued Stribeck contact (`stribeck_friction::StribeckLoadLine`).
  Each substep solves the exact-ZOH endpoint load line `v = v_free - m F`
  against the unregularized kinetic law of `fs-tribo::FrictionLaw::Stribeck`
  and its static cone `mu_s N`. The residual has at most one local maximum and
  one local minimum, so its stable roots are bracketed exactly. Where a stick
  reaction and sliding roots coexist, the contact keeps its branch
  (McIntyre–Schumacher–Woodhouse hysteresis) and changes only at a fold, a
  slip reversal or a lift. No capture window or ramp scale enters the bow
  path. G1 checks root selection against a dense-scan oracle, the hysteresis
  loop against an independent fold scan, and the Coulomb closed form. The
  bowed gates measure a 4..64 substep ladder at the playable point (playable
  on every rung, steady pitch within 12 cents). The former capture window
  (`max(stiction, 0.02)` and a 0.15 m/s flip basin) is gone; at 3.9 N, its old
  gate point, the regime was not converged in substeps under either contact.
  Helmholtz motion is possible with enough modes; it is not guaranteed and
  not a measured rosin curve. For explicit ODE embeddings elsewhere,
  `StribeckFriction` delegates traction to
  `fs-tribo::FrictionLaw::regularized_traction_1d`, reversing the port's
  driver-minus-body velocity to the owner's body-minus-driver convention. The
  shared owner supplies deterministic exponential evaluation and refuses
  malformed coefficients, negative/non-finite normal loads, and
  unrepresentable force. G1/G3 checks compare the force with an independent
  Stribeck formula and the work of equal-and-opposite contact forces. This
  wiring admits no sourced interface data or authored modal losses.
  `BowedStringCard::from_material_specimen` lowers a resolved circular specimen
  using the shared pinned-beam dispersion, Kelvin–Voigt loss and cylinder-air
  drag owners. It requires a matching loss/material identity and frequency
  band, rejects memory, moving supports and polarization splitting, and applies
  runtime bandwidth admission. The caller retains the source receipts; this
  fixed-state, small-amplitude card grants no new material authority.
  Elastic and viscous support shear are evaluated at the bridge, `x = 0`:
  `T y'(0) - EI y'''(0) - eta I y'''_t(0)`. Explicit cards must include the
  viscous bending contribution in their modal damping ratios. Bow velocity
  and bridge-force histories observe the same completed sample endpoint.
  G1 tests compare independent held-load motion, energy and support shear;
  G3 checks shared pluck/bow losses and sensitivity to material viscosity and
  ambient pressure. This does not evolve material state, route damping heat,
  model nonlinear stretch, or establish bow-contact passivity; substep
  convergence is measured at the playable point only.
  `FrictionIsland::InterfaceStribeck` consumes an ordered dry-interface card
  through `fs-tribo`, retaining its identity, history and caller authority.
  The complete sliding law, including velocity strengthening, drives force;
  other friction families refuse. Fixed temperature, load/nominal-area pressure
  and actual slip must remain inside the card domain at admission and both
  endpoints of each substep, including sticking. These sampled checks do not
  bound interior extrema. G1/G3 compare signed, load-dependent motion to an
  independent oscillator and check runtime domain departure. This path supplies
  neither measured-data admission nor heating, wear or evolving contact area.
  Bowed plate outputs report signed volume velocity in m³/s, not point velocity.
  Plate observers admit positive finite listener distance before stepping and
  refuse non-finite volume flow or pressure. G0/G3 checks cover malformed
  distances and inverse-distance pressure scaling with identical motion.
  This compact observer adds no propagation delay or general near-field claim.
  The material string bank can retain both transverse displacement planes
  through `MaterialStringSource::with_transverse_polarizations(TransverseDirections::try_new(bow, bridge)?)`.
  Both planes use the same original circular section, installed tension, density,
  bending stiffness, material loss and gas drag. Their physical directions project
  displacement and force with the same work-conjugate map; shared body compliance
  determines polarization mixing without an authored detuning. Coordinate order is
  first transverse plane, second plane, then optional material torsion. Surface
  bow/pluck ports include the original radius-times-angle contribution; centreline
  loads and unilateral stops carry no invented radial torque. Stops follow the
  source's declared centreline/bow axis. Existing explicit torque/rotation bridge
  connections compose with both planes and retain their own source-derived angular
  maps and passive laws. Material body/cavity pressure and scheduling remain in the
  original joint solver.

  `frankensim-material-strings-v3` extends the v2 section/support grammar with one
  required `polarizations single` or `polarizations two BOW_X BOW_Y BRIDGE_X BRIDGE_Y`
  record after each string section. Both directions must be finite unit vectors.
  Existing v1/v2 files and optional `rotational-bridges` records retain their
  semantics. `polarized-pluck.material-strings` exercises the data-to-PCM route.
  This circular, frozen-state linear model does not supply arbitrary
  six-degree-of-freedom bridge motion, independent stop-normal geometry, shared
  nonlinear two-plane tension, helical winding/slip, or a measured-instrument timbre
  guarantee. Nonlinear material source conversion uses the existing shared
  two-plane Kirchhoff–Carrier stretching law.

  Seven focused physical tests cover material-mode identity, rotational
  covariance, directional sound changes, force work conjugacy, constructor
  ordering/refusals, combined rotational-bridge response and two-plane unilateral
  contact. Three example tests cover source input, material/angle swaps,
  PCM/callback identity and admission. Native verification is recorded separately.
  `ReciprocalStringBank::with_cylinder_gas_loading` selects a causal gas
  operator for the same circular material strings. The shared
  `air_path::cylinder` owner derives the unsteady incompressible Stokes impedance
  `Z(s)=pi*rho*a^2*s+4*pi*mu*z*K1(z)/K0(z)`, `z=a*sqrt(rho*s/mu)`, from
  actual outer radius and the retained complete gas. The continuum follows
  [NASA/CR-2007-215050, equation (7)](https://ntrs.nasa.gov/api/citations/20080006726/downloads/20080006726.pdf).
  Its positive spectrum per log relaxation rate,
  `8*mu/[pi*(J0^2+Y0^2)]`, supplies physical free inertias joined by passive
  dashpots in the original exact reciprocal system. Both transverse planes
  carry instantaneous added gas mass and viscous memory; original solid
  density, bending, installed tension, material viscosity and torsional inertia
  remain distinct. The active transverse card replaces sampled gas resistance
  with intrinsic material loss plus the explicitly collapsed slow tail, so
  fluid damping is not counted twice. The same intrinsic-loss owner supplies
  Kelvin–Voigt bending and the construction's independent tensile loss factor
  in the loaded mass basis; selecting causal fluid loading must not erase
  either. Prony bending keeps its internal strains alongside that tensile
  loss, without inventing a second bending-viscosity law.

  All material-string performance versions accept the optional
  `cylinder-gas-loading X_MIN X_MAX INTERVALS` row before an optional
  `prepared-linear` row and then `cavity`, after any section, polarization and
  rotational/translational-bridge records. The same
  original geometry/gas supplies the operator; no extra material coefficients
  are authored. This binding precedes additional bows, stop profiles, static
  preload and the physical pressure/PCM schedule. File construction counts all
  mechanical/fluid coordinates and string/body/torque/cavity/tensor-loss ports
  against the explicit mode, connection and complete setup-work allowances.
  The source fixtures and regressions compose material/mixture substitutions,
  regional body loss, both transverse planes, torsion, cavity, stops and release.

  These fluid components append after physical bodies. Full tensor body loss,
  cavity connections, pressure observation, dry interface and torque maps
  survive binding, either cold torque/fluid order, and multiple-bow scheduling.
  Pressure and force maps follow the changed transverse mass normalization;
  fluid coordinates receive no fictitious direct acoustic pressure weights.
  Added kinetic storage, reciprocal reactions and viscous heat enter the same
  complete work balance. Static initialization admits exactly zero-stiffness
  positions only when their applied static force is exactly zero. Their
  already-zero positions remain gauges while their velocities remain physical
  free inertias; no artificial spring or general singular preload solve enters.

  `CylinderFlowResolution` declares independent slow/fast dimensionless
  wavenumber cutoffs and log-grid intervals: `0<X_MIN<=0.5`,
  `X_MIN<X_MAX<=3000`, and `1..=4096` intervals. The exact integrated slow
  tail is collapsed to resistance with dynamic error bounded by
  `R_low*lambda_min/omega` at positive frequency; the fast tail is truncated.
  Refine all three axes and the actual coupled response. Every branch adds one
  coordinate and one connection per transverse mode. The shared mechanical
  connection ceiling is 4096, with caller mode/connection/setup budgets still
  enforced before allocation; that port capacity alone does not raise the
  ordinary exact connected-island limit. An explicit prepared linear budget
  admits the complete larger system as described below. Repeated/live bindings and loaded elastic
  modes outside their original material bands refuse.

  Operator regressions compare complex force with an independent modified
  Bessel reference and exercise every cutoff/grid axis. Material-bank tests
  exercise pressure-dependent inertia, viscosity-only changes in actual body
  pressure, both transverse planes, unchanged torsion, fluid kinetic energy,
  torque binding order and preload/release. Source review is not native test
  execution. This local long-cylinder, small-motion Newtonian incompressible
  model supplies no finite-amplitude convection/wakes, axial/end/confinement
  or neighbouring-string interaction, rotary drag, rarefaction, compressible
  radiation, evolving atmosphere or measured violin/cello/bass fidelity claim.

- Reed path: quasistatic or massive Bernoulli valve. Reeds with explicit
  positive stiffness derive their effective pressure area from the
  declared closing relation `A = k H / Pc`. Characteristic-line motion and
  ODE pressure force/swept-face flow share that area; a fixed 25 mm face
  length no longer overrides explicit stiffness and closing pressure.
  This is an effective modal area, not measured geometry. The zero-stiffness
  legacy reduction still derives stiffness using a 25 mm face length, and
  damping is now an explicit finite nonnegative ratio, consumed as
  `c = 2 zeta sqrt(k m)`. Reduced valve cards pass `zeta = eta/2`, preserving
  their existing `c = eta sqrt(k m)` model. This removes the runtime's fixed
  0.35 substitution; authored fixture ratios remain explicit assumptions,
  not measurements. G1 covers actual velocity decay; G4 rejects invalid ratios.
  G1 tests exercise the force balance and actual moving-reed acceleration.
  The massive characteristic-line junction now uses the current opening and
  swept-face flow `-A v_mid`, plus body volume flow, in the pressure solve. Positive
  velocity widens the opening: the ODE and characteristic paths share this
  sign, conjugate to `F = -A (Pmouth - Pbore)`. It corrects the ODE's former
  positive swept-flow sign; see Silva et al., arXiv:0810.2870, Eq. 27.
  Its fixed-opening residual is monotone, bracketed between zero jet flow and
  zero pressure drop, and uses the same signed Bernoulli primitive as the ODE.
  The massless solver and its Newton Jacobian are not used for massive reeds;
  fast-mode massive samples count as strict fallbacks. G1 tests cover closed
  apertures, reverse flow, and instantaneous pressure/face work cancellation;
  G3 covers mass-dependent bore pressure at equal static compliance.
  Characteristic-line linear mechanics now uses midpoint position and velocity,
  with that same average velocity in swept flow. Eliminating the linear reed
  response produces a positive effective load for the existing aperture solve.
  G1 checks the discrete spring/kinetic energy change against pressure work,
  Bernoulli loss and viscous loss without contact, plus energy preservation and
  second-order phase convergence of the isolated oscillator.
  Aperture geometry and elastic lay force are held at the start of the step.
  Single-slit Hunt–Crossley loss uses the same midpoint velocity as mechanics:
  the active reaction is `elastic - elastic*chi*v_mid`; a second zero-reaction
  solve handles the nonadhesive unloading clamp. The coefficients come from
  the admitted fs-dcontact obstacle, restricted to one unit opening coordinate.
  G1 exercises loading, unloading and within-step reversal; a lagged-velocity
  negative control injects energy where the accepted loss removes it.
  This proves the sign of loss work against the held elastic reaction, not
  that the latter matches the change in nonlinear contact potential.
  The isolated ODE reed retains its own pHS stepping scheme. This is not a
  nonlinear contact or whole-bore energy certificate: contact work, coupled
  geometry/time refinement and onset validation remain MR68 obligations.
  Isolated
  cylindrical bores use the `acoustic_chain` ODE with a
  `ViscothermalPin` (massive
  reed = `mass_spring_damper` + jet + face flow). Blow/reed on a
  moving-end string×plate×duct is a leftover `join_port` inlet.
  Moving-end string×duct without a plate is
  `transformer(waveguide, chain, A_inlet)`. Fixed-fixed
  string×duct without a plate shares the ODE clock without
  a force port (`φ(0)=0`).
  Reed lay
  is `fs-dcontact` (not a private Hunt–Crossley).
- `ThinPlate` modes come from `fs-plate` + `fs-modal` (certified
  eigenpairs), not caller Hertz. Mouth pressure × area drives the
  plate; plate volume velocity returns to the waveguide
  (structure–bore).
- String obstacles are `fs-dcontact` power-law potentials on sine
  collocation. On the Kirchhoff–Carrier path the conservative
  potential is `ContactStorage` wrapping the KC Hamiltonian, so
  Gonzalez sees the contact energy. Optional `mu_kinetic` adds
  `fs-tribo` Coulomb traction as a modal port force (not a gradient
  of `H`). The linear modal path still applies both as a force.
  Hunt–Crossley `χ` is a scenario field on the obstacle and a
  dissipative port force (`modal_hunt_crossley_forces`), never a
  term in `H`.
- Reed and structure–bore time port: TMM `R(ω)` sampled on the DFT
  grid and inverse-transformed to an FIR. Isolated linear blow (no
  body) is the same `DelayedFilter` object filled with `IFFT[Z_in]`
  (or the mouth transfer `p/u` when the end is open). A vented
  reflectance FIR does not ring a measurable period; the impedance
  FIR does, so tone-hole shortening is TMM-emergent on one port
  type. There is no one-pole `TravelingWaveLine` fallback on those
  paths. Unflanged-open Nyquist `ka > 1` is refused; a
  flanged mouth uses the Rayleigh piston and is not.
- Card-driven acoustic compilation: `material_assembly::compile_material_assembly`
  resolves exact `MaterialCard`/`QueryPoint` inputs with explicit claim selection,
  then calls the existing circular-string, uniform-plate and regional-chart
  specimen binders. Only consumer-needed scalar properties are resolved: rho/E
  for elastic strings, the selected in-plane plate constants, and additional
  viscosity or equilibrium/modulus/time fields for explicitly selected sourced
  string/plate bending loss and explicitly selected isotropic plate thermal loss.
  Inactive property families are not queried. Loss selectors
  are retained; old numerical loss coefficients and citations are rebuilt.
  `StringPrestressBinding::ThermalExtension` selects the existing alpha(T)
  integrator and sourced thermal string binder. Its reference point and alpha
  claim selection are explicit; the immutable card and current point come from
  the same string source as density, elasticity and selected loss properties.
  Mechanical/loss pins and expansion pins cover their separate requirements.
  Piecewise-linear expansion curves are integrated across all crossed knots;
  current alpha times temperature change is not substituted for that integral.
  Fixed supports use the reference stress-free length to derive elastic strain
  and tension. Fixed mass or current radius remains an independent constraint.
  The specimen retains both mechanical and integrated-expansion receipts, and
  the ordinary acoustic realizer uses the resulting tension and rebuilt loss.
  `Prescribed` keeps explicit force/extension/tuning constraints and raw authored
  eigenstrain; it does not query or claim sourced expansion. Invalid paths,
  missing/ambiguous data and slack/compressive states refuse through the existing
  owners. No thermal transport, finite deformation or melting is added.
  `CompiledMaterialAssembly` keeps the numeric assembly and source-bound specimens
  immutable, and its `realize` method calls the existing ordinary/chart dispatcher.
  Unbound components remain authored. Material conditions remain separate from
  ambient-gas conditions; this compiler does not infer thermal equilibrium or
  upgrade the authority of the whole assembly. Geometry, constraints, forcing,
  supports, damping choices and solver controls retain their existing meanings.
  Resolution/binding errors preserve component context and the owning diagnostic;
  failure leaves the input and prior compiled specimens intact. A thermal plate
  selects the existing six-property thermoelastic resolver; optional bending
  viscosity is requested in that same atomic bundle, including exact claim pins.
  The resulting state rebuilds the fixed-mass/thickness specimen and both loss
  mechanisms. Stale numerical coefficients and citations are replaced together.
  Regional bindings explicitly select `RegionalPlateLoss::Authored` or
  `KelvinVoigt`. The latter resolves every region's loss with its elastic state:
  isotropic regions require uniaxial bending viscosity, while orthotropic regions
  require all six reduced plane-stress viscosity coefficients, including signed
  or explicitly zero couplings. The original tensor owner handles rotation and
  admissibility. `realize_assembly_with_material_plate_chart_and_options` keeps
  the resolved chart through the common acoustic dispatcher, so its existing
  modal curvature-energy loss projection and frequency-domain refusals cannot
  disappear at a bare-chart conversion. Authored damping must be zero alongside
  sourced loss. This retains the existing diagonal weak-loss image, not full
  nonproportional damping. G1/G3 compare uniform eta/E against its analytic modal
  decay, substitute one region's loss source, and observe changed coupled
  pressure and free decay without changing elastic frequencies or the gas state.
  An orthotropic material selection with this isotropic thermal law refuses.
  Specimen temperature must be an explicit positive absolute T coordinate; the
  surrounding gas does not supply it implicitly. Compilation is
  synchronous with no cancellation guarantee; ordinary runtime admission runs at
  realization. This is a Rust description API, without a new wire format, global
  scenario identity, optical binding or conservative material evolution.
  G0/G1/G3 exercise explicit-binding pressure parity, consumer-minimal exact-source
  retention, simultaneous string/plate card substitution at fixed mass, regional
  chart substitution through moving-end string/cavity pressure, active Kelvin–Voigt
  and Prony coefficient rebuilding, and missing-member/property/domain/law refusals.
  Thermal compiler tests vary specimen temperature through source curves and
  check changed pressure in linear, sine von Karman and sampled-FE plate paths,
  with and without separately selected viscosity. A fixed-mass card swap checks
  rebuilt thickness and both damping addends against independent Zener and
  Kelvin–Voigt arithmetic; conflicting heat-capacity claims require the selected
  pin to survive compilation. These fixtures validate implementation, not an
  experiment or a transient heating simulation.
  String thermal compiler tests independently calculate integrated strain,
  fixed-mass area/inertia, tension and stiff-string frequency, then measure
  emitted-pressure pitch at two temperatures with and without sourced viscosity.
  The viscosity fixtures supply scalar coefficients at explicitly supported
  temperature states, preserving the existing frequency-wide constant-law
  requirement; sampled coefficients are not promoted to broadband loss laws.
  They also check exact alpha selection independent of mechanical pins, force
  control without alpha resolution, and unchanged inputs after path/slack refusals.
- Uniform plate material binding: `with_uniform_plate_material_state` consumes
  the same `ResolvedMaterialStatePoint` carrier as the circular-string and
  disc adapters. Isotropic elasticity requires rho/E/nu; principal-axis
  orthotropic plane stress requires rho/E1/E2/nu12/G12. It retains every original
  property receipt and checks exact scalar quantity schemas. Card axes map to
  section material axes, with reciprocal nu21 on exchange. The explicit
  `ThinPlate::material_angle_rad` rotates that uniform section counterclockwise
  from rectangle x, including D16/D26 bending/twisting coupling. Geometry,
  supports and forcing remain fixed. The angle participates in the v2 specimen
  identity. No out-of-plane constants are invented in either
  linear or FE-sampled acoustic plate construction. A fixed-thickness or
  fixed-mass prescription derives thickness and mass coherently; geometry and
  material/model identity are retained separately from supports, pretension,
  authored damping and solver controls. These are independent frozen specimens,
  not a conservative thermal evolution. The purely elastic binder refuses stale
  thermoelastic addends. `with_uniform_isotropic_thermoelastic_material_state`
  instead rebuilds the mechanics and loss inputs from the existing resolved
  thermal bundle, retaining its complete receipts and optional viscosity data.
  The source-bound orthotropic adapter also refuses geometric
  nonlinearity: the existing sampled nonlinear plate's isotropic membrane
  approximation is not a resolved orthotropic membrane law. Anisotropic thermal
  loss and phase changes remain unavailable here. The raw
  nonlinear path checks E, nu and G before selecting isotropic sine modes.
  Its rotated FE bending still uses the documented isotropic membrane
  approximation. The modal search uses a sine trial Rayleigh quotient;
  it is not an exact off-axis frequency or a completeness certificate.
  G0/G1/G3 cover receipts, reciprocity, fixed-mass geometry,
  semantic aliases, model boundaries and real pressure changes from material
  and axis swaps, with aligned eigenfrequencies checked against Navier theory.
  An arbitrary angle changes actual assembly pressure; a quarter-turn recovers
  the pressure and eigenfrequency of reciprocal axis exchange.
- Uniform material-bound plate viscosity: `ResolvedPlateSpecimen::with_kelvin_voigt_bending_loss`
  selects proportional isotropic Kelvin–Voigt bending. The dimension-only
  `kelvin_voigt_bending_viscosity` [Pa s] claim must be nonnegative and carry an
  explicit frequency band. Rho, E, nu and eta must be validity-wide scalar
  constants or temperature-only curves at fixed positive absolute T; their
  frequency bands intersect. Temperature, moisture and other
  coordinates remain at the admitted state point, with original uncertainty
  and observation status retained. The assumed viscous plane-stress tensor has
  the same Poisson ratio as the elastic tensor: the moment law is
  `M = D (curvature + (eta/E) curvature_rate)`, giving `C_b = (eta/E) K_b`.
  This is a proportional constitutive specialization of
  [Kelvin–Voigt strain-rate damping](https://doc.comsol.com/6.3/doc/com.comsol.help.sme/sme_ug_theory.06.029.html#3441080),
  not a measured general bulk/shear viscosity tensor.
  The ordinary linear, sine von Karman and sampled-FE plate paths consume it.
  Each retained reference-mode frequency must lie in the admitted band, even
  for a zero-strength viscosity claim. Only bending energy receives this loss;
  geometric prestress does not. For FE modes the diagonal coefficient is
  `(eta/E) phi^T K_b phi / (phi^T M phi)`. Off-diagonal damping coupling is
  omitted, a modal approximation when prestress makes damping nonproportional.
  The nonlinear path retains bending-rate damping but supplies no nonlinear
  membrane viscosity or guarantee that finite-amplitude frequency shifts and
  generated harmonics remain in band. Authored scalar damping alongside this
  law refuses. Explicit thermal loss remains separately additive; its inputs
  must exclude this viscosity's contribution to avoid double counting.
  `compile_material_assembly` resolves the selected viscosity with the plate's
  elasticity and density. A material rebind recomputes its coefficients and
  source identity rather than carrying the old material's viscosity forward;
  missing source data refuses. Thermal rebinding of an already viscous raw
  plate refuses because it would otherwise retain stale viscosity inputs.
  G1 checks the actual sine-mode damping operator and ringdown against the
  analytic Kelvin–Voigt plate equation with timestep refinement. Independent
  geometric-energy integration checks prestressed DKT and sampled-FE damping,
  including additive thermal loss. G3 changes only viscosity through the card
  compiler, retaining elastic eigenfrequencies while changing actual assembly
  pressure in all three plate paths. Fixed-temperature curve tests retain the
  source interpolation receipts and compare plate frequencies/damping, every
  Prony string branch and actual pressure against independently evaluated
  coefficients. The thermal string compiler also exercises interior T values
  with simultaneous integrated expansion and bending viscosity. These are
  frozen-state calculations, not constitutive evolution during heating.
  These synthetic implementation fixtures
  are not experimental material/ringdown validation.
- Regional linear plates: `ResolvedPlateSpecimen::section` exposes the rotated
  numerical section for assignment to `PlateChart` elements. The caller retains
  the resolved specimens for source receipts; a numerical chart does not mint
  material authority or a regional provenance identity. The full-rectangle
  specimen mass is not reused as regional mass: fs-plate integrates local rho*h
  and rotary inertia from each actual section.
  `compile_plate_material_chart` additionally binds disjoint named regions to
  exact resolved material states, constitutive/axis choices and fixed-thickness
  or fixed-region-mass prescriptions. The immutable `ResolvedPlateChart` retains
  these inputs alongside its numerical chart. It integrates actual triangle
  areas to derive h = m/(rho A_region); `with_region_material` rebuilds sections,
  mass and identity using the retained constraints. Failure leaves the original
  binding unchanged. Empty/duplicate names, overlapping or missing assignments
  and out-of-range triangles return `PlateAssignment`; geometry and section
  refusals preserve the owning `PlateError`. No fallback material is inserted.
  The specimen identity includes mesh geometry/connectivity, boundary-node set
  and per-element source/model/thickness/angle in mesh order. Names, assignment
  order and repartitioning with identical resolved per-element inputs do not
  alter it. It is not a complete acoustic-scenario or authored-constraint identity.
  Compilation is synchronous, with no cancellation guarantee or material-state
  evolution; cloned numerical charts do not automatically retain source authority.
  `certified_chart_radiators` consumes that chart with `PlateChartRadiation`
  (explicit eigenvalue search interval, mode budget, support/prestress, authored
  damping and nodal unit-force footprint). All force weights must be finite
  and sum to one; supported forces react into the support without redistribution.
  `realize_assembly_with_plate_chart` feeds the resulting modes into the same
  assembly dispatcher, body ports, piston loads, string/cavity/duct coupling
  and pressure observer as uniform plates. A simultaneous `assembly.plate`
  refuses. Extra authored lumped bodies still compose. The modal interval
  need not include the fundamental; certification applies only in that window.
  It does not certify stability or exclude buckling modes outside that window.
  `ResolvedPlateChart::with_kelvin_voigt_bending_loss` selects regional scalar
  loss (`RegionalBendingViscosity`): isotropic regions use the proportional law
  above; orthotropic regions use their six explicit `plane_stress_viscosity_*`
  coefficients, resolved and rotated by the same owner as the full tensor image
  (`plate_viscosity::resolve_region_viscosity`; no coefficient is inferred). The
  one diagonal evaluator, `plate_viscosity::DiagonalPlateViscosity`, gives
  `zeta_i = sum_e w_r q_e^T K_e q_e / phi^T K phi`: isotropic regions use the
  elastic DKT matrix with `w = tau omega / 2` (the legacy arithmetic, moved bit
  for bit); orthotropic regions use the DKT matrix of the rotated viscous tensor
  with `w = omega / 2`, equal to `C_ii / (2 omega_i)` of the full projection for
  an exact eigenpair. This image omits off-diagonal modal damping;
  `PlateViscosity` projects the full tensor, and a chart with scalar loss
  selected cannot also select it. Tests:
  `isotropic_regional_loss_routed_through_plate_viscosity_is_bitwise_legacy`
  (zeta and radiator bits against a verbatim pre-move copy) and
  `orthotropic_region_now_admitted_and_matches_full_matrix_diagonal` (1e-13 plus
  twice the measured eigenpair residual).
  No spatial thermal-loss law, nonlinear membrane, remeshing, phase evolution
  or full-wave radiation is claimed. Interfaces are perfectly bonded, with
  thickness centered on one flat mid-surface as specified by fs-plate.
  G3 checks uniform/element partition pressure equivalence, legacy uniform
  path parity including moving-end/cavity ports, and pressure/mode changes
  from local density, stiffness, thickness and grain on an unstructured chart.
  G0/G1/G3 additionally check assignment refusals, exact source retention,
  integrated regional mass, region-order/partition pressure equivalence and
  actual pressure changes after named material replacement under both mass
  and thickness constraints. These use synthetic materials, not experimentally
  validated material data or measured sound.
  G0 checks force footprints, modal controls and conflicting descriptions.
- Linear plate projection: the eigenvector, its modal mass `phi^T M phi`,
  signed radiating area `integral(phi dA)`, and signed bridge-force projection
  remain in one basis. Surface integrals use triangle-area P1 nodal quadrature;
  the declared drive footprint is uniform traction over the generated
  rectangle's first cell strip, normalized to one unit of total force.
  There is no mass inferred from drive amplitude, absolute-value radiation,
  or artificial area floor for cancelling modes. Generalized modal force is
  `b_drive F + A_mode p_face`. The pHS ports divide both coefficients by
  `sqrt(modal_mass)`, and its observers convert acceleration back by the same
  factor. Thus eigenvector sign/scale cannot change physical pressure or work.
  This follows modal virtual work and the compact limit of the
  [Rayleigh integral](https://euphonics.org/4-3-2-the-rayleigh-integral-and-the-baffled-piston/).
  The optional diagonal piston self-load uses the physical panel area to set
  its equivalent circular radius; input velocity is `A_mode q_dot/A_panel`
  and generalized reaction is `-A_mode p_face`. Modal normalization cannot
  change that geometry. The approximation omits mutual radiation impedances,
  higher multipoles and higher-order DKT trace reconstruction; eigenvalue
  certification does not certify radiation accuracy. G1/G3 check manufactured
  quadrature/mass, zero-monopole cancellation, mode sign/scale invariance with
  and without the actual fitted self-load, conjugate pHS ports, and density
  scaling through the actual assembly pressure output.
- Plate damping: a single authored viscous ratio is the two-point
  Rayleigh fit (`fs-material::visco`) through that ratio at ω₀ and
  4ω₀ so higher certified modes sit on the stiffness limb. Both linear and
  von Karman banks add Zener thermoelastic loss only from explicit specimen
  inputs. `with_isotropic_thermoelastic_state` binds a resolved six-property
  bundle to the plate, moving elasticity, density, thermal properties, and
  temperature together and retaining its state identity. There is no
  density-threshold material guess or fixed-temperature substitute. Missing
  thermal inputs leave that loss unavailable/unclaimed; requesting the scalar
  law on an orthotropic plate refuses (including a mismatched shear modulus
  when E1=E2). The law estimates the first through-thickness thermal relaxation
  mode of a homogeneous isotropic thin beam/plate; it does not solve a thermal
  field or establish an anisotropic loss law. Callers must exclude this addend
  when authored damping or resolved thermal transport already includes it.
  With prestress, the thermal addend is weighted by the mode's bending/total
  strain-energy ratio: `phi^T K_b phi / phi^T (K_b + K_G) phi` for DKT modes,
  and `D k^2 / (D k^2 + T)` for the isotropic sine modes. Authored loss is
  unchanged. The same projection reaches the sampled-FE von Karman bank.
  This is weak-loss, diagonal modal bending damping about the prestressed
  equilibrium; it does not model finite-amplitude membrane thermal loss or
  off-diagonal thermal coupling. The energy separation follows
  [Fedorov et al., Eq. (1) and Sec. II](https://arxiv.org/abs/1807.07086).
  G1 tests independently integrate geometric energy and compare actual
  unforced sine-mode ringdown to the analytic damped oscillator under refinement.
  G0/G1 tests inspect both assembled damping operators against independent
  Zener arithmetic, vary equal-density specimens and temperature, and cross
  the former density threshold without switching laws. Radiation
  reaction on linear compact radiators uses the sole Bessel–Struve piston
  owner in `fs-phs`, converting acoustic p/U to face p/v by the physical
  area, then fitting and passivity-repairing with `fs-vfit`. Evaluation and
  fitting refusals propagate; there is no silent low-ka fallback or unloaded
  assembly success. `fs-bem` is not a production dependency of this crate.
  G0 `thin_plate_refuses_instead_of_low_ka_fallback` exercises this refusal.
  The existing compact one-branch load keeps ka<=0.5 and dt*fmax<=0.1;
  its analytic leading resistance and inertance coefficients are tested at
  1e-14, with a 512-frequency grid targeting <2.2% complex and <4.2%
  resistance error versus the closed form (native measurement pending).
  These are sampled discrepancies, not continuum enclosures. Linear plates
  accept in-plane pretension and clamped edges. Von Karman is the
  isotropic SS sine path when `e1 = e2` and edges are free to
  rotate; clamped or orthotropic bending uses DKT-sampled
  displacement with the same sine Airy membrane channel
  (`von_karman_sampled_plate`). The compact far-field observer is
  the baffled on-axis piston `p = ρ A ÿ / (2 π r)`, the same
  half-space as the self-load. A full BEM observer is the
  solver-neutral `DirectionalFarFieldTable` /
  `far_field_observer_pressure` (`p = F e^{ikr}/r`); `fs-bem`
  produces `FarFieldTable` as a **dev-dependency** only (the
  production cycle `couple → bem → solver → feec → couple` is
  refused).
  `CompactBody::drive` reuses the shared exact-ZOH modal transition under
  held force, converting through mass-normalized coordinates. Acceleration,
  volume velocity and pressure describe the completed endpoint. G1/G3 tests
  compare analytical damped/undamped motion and held-force subdivision,
  including a step beyond the former Euler stability limit. Radiation-filter
  reaction remains sampled and held, not an exact continuous coupled solution.
  If motion refuses after the filter advances, the caller must discard the run;
  the filter memory is not rolled back.
  The finished pressure history is passed through ISO 9613-1
  absorption (`air_path`) with the assembly's explicit humidity;
  Stokes–Kirchhoff is only the fallback outside the ISO window.
  ISO is applied only with an `fs_material::iso9613::TerrestrialAir`
  admission bound to the gas state. The scenario's legacy moist air makes its
  own. An explicit `GasAbsorptionModel::Iso9613Air` carries one made from the
  gas's composition, so heliox, oxygen-poor or CO2-rich gas refuses with
  `AcousticRealizeError::NonTerrestrialAir`, naming the species.
  The receipt's `relaxation_omitted` is
  `OmissionBound::StandardAirConvention`: standard air's worst ISO absorption.
  It is a convention and not a bound for a gas that relaxes more strongly than
  air (CO2-rich, for example).
  `air_path::omitted_observer_absorption(record, range, band_top)` is the
  figure a render that applies no law (or Stokes-Kirchhoff only) leaves out,
  read from the listener gas's record, never inferred from its state: ISO
  9613-1 for a legacy record or terrestrial-air composition; Stokes-Kirchhoff,
  complete, for monatomic species only; the relaxation law with ISO's rates
  held (declared) for N2 and O2 off ISO's fractions; otherwise `Unbounded`,
  naming the species (heliox's He, CO2 above ISO's ceiling, SF6). The band top
  is `observer_band_top_hz`: the declared band, Nyquist and 20 kHz
  (`OBSERVER_BAND_TOP_CAP_HZ`), whichever is lowest.
  `ObserverOmission::admit` refuses with
  `AcousticRealizeError::ObserverOmission` when the omitted figure exceeds
  `OBSERVER_OMISSION_TOLERANCE_DB` (1.0 dB) or is unbounded, unless the input
  declares a law (`observer_law none` or `stokes-kirchhoff`), which accepts the
  omission with its receipt. In air at 20 °C, 50 % RH and 1 atm, at 20 kHz, a
  1 m listener omits 0.524 dB and passes; a 2 m listener omits 1.048 dB and
  refuses. A gas without a record naming its composition (a resolved field
  cell) is `Unrecorded` and is treated like `Unbounded`.
  `GasAbsorptionModel::MolecularRelaxation` applies the gas's own relaxation
  law (`fs_material::relaxation`). It is evaluated from the gas's composition
  and its species' data, and supplied as `RealizeOptions::observer_relaxation`
  at exactly the realized gas state. It is applied as the causal transfer
  `air_path::propagate_pressure_history_relaxing`: absorption and dispersion
  come from one wavenumber, Stokes-Kirchhoff is counted once, and the bulk
  delay is not applied. The receipt names it, with nothing omitted.
  These refuse rather than being ignored:
  - the law without its state;
  - a state of another gas;
  - a state supplied with any other law.
  `xtask/tests/observer_relaxation_amount.rs` gates the amount on a rendered
  plucked string at 100 m. The humidity change is gated against ISO 9613-1,
  and the oxygen change against the law with declared rates.
  The shared `pcm_wav::baffled::BaffledPressure::new_multirate_gas_differentiating`
  composes these same `ObserverFilter` laws on every original surface-to-receiver
  path before signed modal summation. Only bit-equal computed distances share a
  prepared filter; there is no common-distance filter on the summed pressure.
  Each path retains its `rho*dA/(2*pi*r)` gain, source filter receipt and residual
  delay `r/c - filter.low_frequency_delay_s`. A negative residual receives the
  smallest common whole-sample causal padding across ALL paths; that added
  numerical latency is reported separately and preserves relative path timing.
  The existing 65-tap sinc contributes its original 32-sample interpolation
  latency. The source response's own delay is counted once, including dispersion
  and its numerical padding; FIR support is not physical travel time.
  Explicit distinct-path, retained-response and convolution/projection budgets
  bound cold preparation, with cancellation between bounded designs and surface
  tiles. The actual pressure loop and clone retain the original receiver owner.
  The reference cold budgets are 16,384 distinct paths, 4,000,000 source-filter
  taps and 500,000,000 convolution/projection terms. These are admission limits,
  not a measured real-time guarantee. The receiver exposes its actual multiply
  count per sample; the live-piano clock reports that count and the physical and
  numerical delays separately. Larger responses, including some reduced-pressure
  or higher-rate configurations, can exceed the declared limits and refuse;
  the caller may select larger explicit cold budgets in the native API.
  The law is explicit: classical transport does not infer molecular relaxation,
  ISO requires its state-bound terrestrial-air admission, and a supplied
  molecular state must belong to the same gas. Prepared live-piano atmospheres
  now use this path with the generic classical law; native configurations can
  supply another admitted law or retain the lossless comparison receiver.
  This remains fixed-gas, infinite-baffle, one-way observation, with original
  source-law, filter-tail, temporal and spatial approximation limits.

  Streaming renderers apply the same laws through `air_path::ObserverFilter`,
  a causal FIR built once per law, range and clock. The reed performance
  (`render::schedule::reed`) applies it to both observers, the legacy
  `listener_m` proxy and the `baffled_outlet` receiver. It chooses as above:
  ISO 9613-1 for a legacy (T, p, RH) medium, Stokes-Kirchhoff (with the
  omission receipted) for a supplied mixture or thermochemical medium, or the
  `observer_law none|iso9613|stokes-kirchhoff|relaxation
  iso9613-terrestrial-air` record. `none` (`GasAbsorptionModel::Lossless`) is
  a receipted lossless choice for comparisons with lossless analytic oracles,
  never a default. `ReedPerformance::observer` returns the receipt.
  How each law is made causal:
  - Molecular relaxation: its own causal transfer.
  - ISO 9613-1: each O2 and N2 term by its exact complex Debye
    (Kramers-Kronig) pair, so the magnitude is ISO's exactly and the dispersion
    is its own. The classical term is the in-band minimum phase.
  - Stokes-Kirchhoff: the in-band minimum phase. The receipt names the band
    where `omega tau_v << 1` fails.
  A minimum phase built from the in-band magnitude is not used where
  dispersion exists, because it cannot see absorption above Nyquist. For air
  at RH 50 % and 30 m, against the direct causal relaxation transfer:

  | Grid Nyquist | Worst in-band complex error | LF delay (physical 34.5 us) |
  |---|---|---|
  | 24 kHz | 2.0 | 20.2 us |
  | 48 kHz | 1.07 | 27.4 us |
  | 192 kHz | 0.11 | 33.8 us |

  A dispersive transfer sampled only to Nyquist rings before `t = 0`. Its
  excess phase at Nyquist is nonzero, so the periodic DFT spectrum jumps there,
  and no zero-latency causal FIR drops the precursor below a fixed fraction.
  That is a property of band-limited dispersive transfers, not of one law. Its
  FIR is therefore delayed by a sub-sample alignment `d` (making the Nyquist
  bin real) plus the fewest whole samples `D` for the precursor to meet the
  certificate. Each receipt's certificate (`budget_energy`) is
  `min(OBSERVER_TRUNCATION_ENERGY = 1e-9, OBSERVER_UNCERTAINTY_MARGIN = 1e-3 x
  alpha_uncertainty_energy)` of the response energy. `alpha_uncertainty_energy`
  is what ISO 9613-1's stated +-10 % on alpha is in the same measure. The
  dropped energy therefore sits three orders below the law's own uncertainty
  at every range, and short ranges pay more padding. The receipt prints `d`,
  `D`, the budget, the latency and the realized dropped energy. There is no
  fixed range floor: a baffled receiver whose own delay cannot absorb the
  filter's refuses (`ObserverLatency`, naming `observer_law none`), wherever the
  physics puts that for its law, gas, clock and range. A 5 cm receiver at
  48 kHz in air refuses.
  Precursor energy left before `-D`, ISO companion, air, 48 kHz (prototype):

  | Range | Alignment | D = 16 | D = 64 | D = 256 | D = 1024 |
  |---|---|---|---|---|---|
  | 0.5 m | none | 3.1e-5 | 8.0e-6 | 1.9e-6 | 3.7e-7 |
  | 0.5 m | d | 1.7e-9 | 3.3e-11 | 5.5e-13 | 9.4e-15 |
  | 30 m | none | 2.7e-4 | 7.1e-5 | 1.7e-5 | 3.4e-6 |
  | 30 m | d | 3.0e-7 | 6.3e-9 | 1.1e-10 | 1.8e-12 |
  The filter is referenced to its front. The baffled receiver's propagation
  delays use the receipt's `delay_speed_m_s`
  (`RegionalDuct::baffled_outlet_delayed_at`), so receiver delay plus the
  filter's low-frequency delay (dispersion and latency) is `r / c_0` exactly
  once. A receiver too close to absorb it refuses with
  `AcousticRealizeError::ObserverLatency` (range, available and needed
  samples). The law is applied at the receiver's distance from the outlet
  centre; `range_spread_m` bounds the difference over the Rayleigh elements.
  The proxy has no propagation delay, so its filter delay is uncompensated.
  `ReedObserverReceipt::consumer_delay_s` is the delay a consumer must apply
  to place the proxy's radiated part as the baffled path is placed: both
  total `r / c_0`. The proxy's sample is bore pressure plus radiation (plates
  and jet at `listener_m`). The law acts on the radiation only
  (`ReedBoreVoice::last_observer_parts`); the bore pressure travels no
  distance and passes unfiltered and undelayed
  (`internal_part_unpropagated`). Under `observer_law none` the proxy is the
  voice's own sum, bit for bit.
  `ReedPerformance::observer_admission` reports the listener gas's figure
  (`air_path::omitted_observer_absorption`, from its record; a field cell is
  unrecorded) at the observer's range and band top (the receiver's declared
  band, or Nyquist for the proxy, capped at 20 kHz), and what the applied law
  leaves out: nothing under ISO or the relaxation law, the relaxation part
  under Stokes-Kirchhoff, everything under `none`. It refuses an undeclared
  Stokes-Kirchhoff law whose omission is unbounded, unrecorded or above
  1.0 dB. The library builds regardless; `music_render wind` calls it before
  any output exists and writes the receipt as `reed_input.observer`. The five
  reed CO2 fixtures declare `observer_law stokes-kirchhoff`, the law they
  already had by default, so their sound is unchanged.
  Authored string `ζ` stays a constant modal ratio; a relaxation spectrum
  is present only when the caller explicitly selects `PronyBending`.
  An optional `HelmholtzCavity` faces the plate monopoles as a
  flow-driven pHS whose damper is compact-mouth `Re Z_rad(ω₀)`.
- Bow roughness: an optional `ContactTexture` is a declared
  self-affine height spectrum. `fs-tribo::surface_excitation` samples
  it; the height perturbs the normal load that Stribeck sees. No
  measured rosin is invented.
- Three-way clock: when `EA > 0` the string member is Kirchhoff–
  Carrier with contact inside `H`. Polarizations share the plate and
  the duct FIR. The certified three-pHS Dirac interconnection is
  `common_flow_dirac(moving_end_waveguide, transformer(plate, load))`
  in `fs-phs` (free-end attachment port, area transformer, Kirchhoff
  star). The load is a Helmholtz cavity or an `acoustic_chain`
  duct — the same object. `PrestressedString.moving_end` selects
  that clock. Bow, obstacle stations, and blow/reed are leftover
  ports of `join_port` (Stribeck and Hunt–Crossley stay port
  forces). Von Karman is the same plate pHS with quartic storage.
  `EA > 0` is `kirchhoff_carrier_moving_end` on the same
  free-fixed port. Fixed-fixed sines keep the one-way bridge
  force; with a plate the duct is still the ODE chain
  (bridge force on the leftover plate port), not the FIR line.
  Without a plate the same cylindrical chain is still the ODE
  clock. A cylindrical bore (uniform or stepped) is the
  `acoustic_chain` LC ladder with an all-regime
  `ViscothermalPin` from `GasState` (`μ`, `γ`, `Pr`;
  series `R` and thermal `G` are three-term Foster networks
  collocated to Bessel Zwikker–Kosten `F(r_v)`, all shear
  numbers). A single widening taper is `spherical_cone`
  (`ψ = x p` plus the one inlet Euler near-field shunt).
  Any other taper, alone or in a multi-section chain, is
  frustum LC cells (exact frustum inertance and volume),
  frankensim-d1j6s. A `Closed` termination is the chain's
  rigid outlet (`ChainOutlet::Rigid`, `U = 0`), matching
  the characteristic line's `Termination::Closed`; before
  frankensim-gs17m it was the pressure-release end, an
  octave away.
  Open tone holes are `AcousticTap` side branches on
  the ODE path: a compact neck is one inductor; a long
  neck (`kℓ > 0.2` or `ℓ > 4b`) is the same 2-cell LC
  line plus flanged mouth the TMM chimney already is.
  The same `ViscothermalPin` sits on a lumped neck or
  on each chimney cell of a line.
  A scenario `LocallyReactingWall` is the ODE `WallPin`
  (per-cell shunt LC, `A_w = 2π a dx` with slant on a
  taper) and the TMM `Y' = 2π a slant / Z'` addend on the
  FIR `Z(ω)` / `R(ω)` path, including tone-hole chimneys.
  Same numbers both clocks. The one-pole `TravelingWaveLine`
  is not that path.
  A closed pad is the TMM cavity compliance on the station
  cell; `foster_branches > 0` adds the bore's thermal
  Foster on that remaining `C`. A scenario `ToneHole.open_fraction`
  σ is the same admittance mix
  `Y = σ Y_open + (1−σ) Y_closed` the TMM `HoleState::Vent`
  already is.
  A quasistatic reed is the Bernoulli port on that
  inlet; a massive reed is `mass_spring_damper` plus the same
  jet, face-velocity flow, and a Hunt–Crossley lay as a
  dissipative port force (not a term in H). With a linear plate
  the mouth is a transformer. An unflanged ODE mouth is the
  same compact `(R, X)` as the TMM: `Re Z_rad` plus mass
  `Δℓ = 0.6133 a` on the last flux. A flanged mouth is
  `0.8216 a` (Rayleigh piston above `ka = 1` on the TMM
  path; the unflanged Nyquist ceiling does not apply).
  `foster_branches > 0` replaces the pin `Re Z_rad` with
  a Foster match of that load across a band (unflanged
  samples stay under `ka = 1`; flanged may use the
  piston). Tabulated `R(ω)` remains the FIR path.
  An open chimney with
  `foster_branches > 0` carries the bore's Foster series
  on a lumped neck, or series plus thermal Foster on
  each cell of a long-neck line (extra states). Frequency-by-frequency
  Bessel TMM (`LossModel::Bessel`, spherical-wave cones
  with local-radius lossy substations) is the FIR path. 3-D jet broadband remains a
  no-claim. Plate ×
  Helmholtz in `cavity_phs` is the same transformer, not a
  staggered pair of steps.
  Its modal damping inputs are complete: it never selects aluminum/steel from
  density or takes the solid's temperature from the gas. `PlateCavitySpec::bind_plate`
  uses the ordinary `certified_radiators` owner to rebuild frequency, material
  damping, force participation and monopole area together. Force and flow ports
  are normalized by the same square root of modal mass. Source-bound callers
  pass their `ResolvedPlateSpecimen::plate()` and retain that specimen's receipts;
  this numeric adapter does not mint material authority. Thermal and viscous loss,
  including prestress dilution, are included once by the plate owner. Nonlinear
  plates refuse at this linear adapter. Failed admission preserves previous modes.
  G1/G3 tests carry compiled source curves through this pHS pressure path at two
  independent solid temperatures on both sides of the former density threshold,
  compare damping with Zener plus Kelvin–Voigt arithmetic, and exercise repeat
  determinism and failed rebinding. They are synthetic implementation checks;
  this adds neither thermal evolution nor exterior plate-radiation back-reaction.
- Von Karman geometric nonlinearity is the isotropic simply-supported
  analytic primitive when those hypotheses hold, and the FE-sampled
  Airy construction otherwise. The membrane channel remains the
  sine Airy basis (in-plane movable).
- Reed observer adds the compact far-field dipole of the slit force
  `F = Δp · w · y`. That is not 3-D jet broadband and not the
  fs-aeroac 2-D Curle Hankel kernel. 3-D jet broadband remains a
  no-claim. There is no instrument crate.
- `encode_pcm16_wav` maps pascals through a declared full-scale. It never
  peak-normalizes.

## Invariants

- The Dirac interconnection conserves interface power EXACTLY (to roundoff) —
  the G0 law; incompatible ports are refused.
- Vibroacoustic power identities are exact by derivation and asserted
  at solver roundoff; the input-balance identity is TAUTOLOGICAL for
  any consistently assembled system, so the one-sided interface-normal
  mutation is caught by the `interface == cavity` cross-row residual
  (measured alarm 2.0 vs 1e-16). Repeat solves are bitwise identical.
- Every admitted v2 schema has a shape-compatible power contraction, a checked
  watt-dimensional effort/flow product, and energy. Rotational, electrical,
  and chemical PR-2 kinds additionally require angular momentum, electric
  charge, and amount respectively; legacy seed role vectors remain unchanged.
  Field duality records separate effort/flow function-space roles,
  includes its integration-measure dimensions in the checked product, and
  explicitly assigns the measure to the density side before kind-specific
  dimension comparison.
- Stable relation IDs cannot alias their owned port/boundary IDs. A
  conservative pair requires distinct port IDs, matching physical metadata and
  clock/timestamp, plus outward-from-owner conventions on both sides. PR-1
  refuses common-frame orientations until a public pullback can prove them.
- Storage, dissipation, and source ownership remain separate typed claims. A
  conservative junction cannot stand in for any of them; a dissipative
  relation cannot exist without an evidence reference; a source cannot exist
  without a named accounting boundary whose coordinate/sign binding exactly
  matches the port.
- Every admitted stream has a non-empty duplicate-free constituent axis, one
  shared outward sign convention, finite fixed-dimension rates, fixed bundled
  conservation roles, and exactly one energy chart. Alternate charts cannot
  enter without exact scalar identities and context-bound durable evidence;
  the actual rate axis must equal the proof-bound axis, chemical power is owned
  exactly once, and pressure/deviatoric work uses the single normative
  Cauchy-tension-positive integrated outward-power convention.
- Every accounting report is tied to one exact logical interval, subject-bound
  evidence roles, exact initial/final endpoints, one canonical expected-row
  set, one element axis, one charge policy, and one entropy convention.
  Caller order cannot change arithmetic: rows are sorted by stable contribution
  ID and reduced with a fixed-order compensated sum.
- Global first-law, mass, element, and charge closure use two-sided typed
  tolerances. Entropy generation is never absolute-valued or clamped: only
  `S_gen >= -entropy_tolerance` passes. Individual boundary contributions may
  be negative. Internal pair closure is a separate equality gate on every
  audited axis, and violations retain the exact pair, ports, and rows rather
  than guessing a single cause from a global residual.
- Inventory snapshots and filled contributions retain the exact energy
  reference, element schema, charge policy/basis, and entropy convention under
  which they were admitted; audit refuses same-interval context rebinding.
- A constant-temperature thermal power row must carry only heat entropy and
  satisfy `energy = T * integrated_heat_entropy` bit-for-bit. Radiation on a
  thermal power row requires profile-bound external evidence. Profile-bound
  thermal and stream heat/radiation integrals remain external-evidence claims.
- The energy audit reports an interface-balance failure exactly when some
  caller-supplied exchange has absolute imbalance above `tol` or is non-finite.
- On the added-mass fixture (`μ ≥ 1`): naive staggering (`ω = 1`) diverges while
  Aitken-relaxed coupling converges to `x* = c/(1+μ)`; Aitken never takes more
  steps than a stable fixed under-relaxation.

## Error model

DEF-05 modal apertures now admit relative tolerances strictly below 1/128.
Work and flow closure use the actual declared relative ceiling, with no 128x
post-solve allowance. `ModalApertureFrame` records the work and flow
`SolveCertificate`s (joules and m^3/s), including the observed residual, scale,
requested/effective budgets, gamma-n arithmetic floor and the actual scalar
residual-evaluation count. These explicit budgets never auto-raise below their
floor; a refusal retains `fs_phs::PhsError` through
`AcousticRealizeError::AlgebraicSolve`. No state or sample clock is published.
Work's floor uses `gamma_(4*modes+12)` times the absolute energy/work operands;
flow uses `gamma_(2*modes+8)` times its absolute flow operands. These are
conditional arithmetic diagnostics, not bounds on material, flow-model or
time-discretization error. G0/G4 `modal::refusal_tests` tests the 1/128 boundary
and a 1% unaccounted external-force work defect through the production gate.

Structured `CoupleError`; no panics. In addition to legacy incompatible kinds,
errors cover invalid IDs/empty shapes, pointwise or field-measure dimension
overflow, non-power products, kind-specific or shape/pairing mismatch, missing
kind-required conservation roles, schema-only kinds sent through the raw
interconnection oracle, localized schema conjugacy mismatch, identity aliasing,
accounting-boundary coordinate mismatch, scalar/non-scalar misuse, wrong
added-mass port kind, and non-finite schema-bound values. Stream errors localize
non-finite fields/components, empty or duplicate constituent axes, non-outward
orientation, proof-bound/actual axis mismatch, species-potential accounting on
an element axis, chart/proof binding or identity mismatches, non-positive
density, chemical double counting, and declared/chart energy-rate disagreement.
Window errors additionally localize clock/interval/endpoint mismatches,
evidence-role/subject mismatches, malformed or one-sided internal exchanges,
duplicate rows/ports/elements, missing or unexpected manifest rows, foreign
stream context, missing stream/chemical element schema, unproved charge policy,
negative inventory/tolerance, nonpositive temperatures, incompatible
power-port transfer axes, absent temperature evidence, nonfinite integrated
values, and element/charge axis violations. Physical conservation or
second-law failures are report violations, not construction errors.

## Determinism class

Fully deterministic: schema admission canonicalizes roles by enum order,
stream admission canonicalizes constituents by typed ID, and window admission
canonicalizes manifest, per-exchange, port, and element ledgers by stable typed
identity. Exact crosswalks use a pinned operation order and bit equality;
window reductions use a fixed-order compensated sum. Interconnection, relation
construction, audit, and iterations are pure functions of their inputs.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

- `mm_line` module (bead 3ez8g.4.2): the multimodal characteristic-line
  runtime for brass. `MmLineBank` realizes one exact-FIR plane-mode
  reflectance line PER VALVE COMBINATION from the fs-duct multimodal
  authority. Its length comes from the shared `LineLadder`
  convergence ladder (DEF-02, IRP .82), not from a fixed round-trip
  multiple; an executed lesson is that a four-trip tail left the
  lowest peak 92 cents flat from wrap-around aliasing. The ladder
  starts at four round trips and doubles until every combo's
  fourth-quarter tail energy is below `tail_tol` and its first three
  peaks move less than `peak_tol_cents`. It refuses with
  `LadderUnconverged` at the ceiling. The line is terminated into the
  selected mouth load. `MmLineBank::new` and `MatrixFirLine::{new,with_ladder}`
  accept the shared `fs_duct::MouthLoad`. Its `Tabulated` arm requires a
  matching mouth plane/radius/gas and complete positive-frequency support
  for the requested ladder, including its ceiling's lowest bin and Nyquist.
  One bank uses one common table, so every combo must match its declared
  inlet-relative axial position; differing path-length combos refuse.
  Per-combo tables remain outside this interface.
  Every retained positive bin comes from the table; the bank receipt reports
  `table_bins == n_fft/2` and zero fallback bins. DC uses the table's declared
  finite static impedance through the same multimodal bore at the existing
  1e-3 rad/s quasi-static probe, with no table extrapolation, default -1, or
  clipping. Strict-load errors propagate, and plane reflectance is checked
  before FIR realization with the existing 1e-12 rounding allowance.
  `MouthLoad::Matched` refuses explicitly on these staircase time-domain lines.
  The existing FIR Nyquist rolloff and precursor folding remain in effect;
  full-band source coverage is not a new high-band realization-error claim.
  `MmLoad::TabulatedWithFallback` remains an explicit diagnostic policy.
  A spliced analytic fallback outside the
  table's support is a spectral jump, so the ladder refuses it by
  name. Passivity (`|R| <= 1 + 1e-12` on the 4x-oversampled grid) is
  checked on the stored taps without clipping and refuses otherwise.
  A valve switch rebuilds from those taps. Valve
  switches at block boundaries use the carry-outgoing-history lift
  (in-flight waves persist; the new taps govern future reflections)
  plus an optional linear CROSSFADE from the old line — an
  instantaneous reflectance swap is genuinely discontinuous physics,
  and the crossfade is the no-click lane; every switch is logged with
  its lift name and carried sample count. `MatrixFirLine` is the
  full N^2-FIR bake-off arm with per-mode incoming-energy logging.
  It is sized by the same ladder:
  - DC is the matrix's own low-frequency limit;
  - the convergence tail is the Frobenius sum over all N^2 entries;
  - peak shift and passivity are checked on the plane entry;
  - its `MatrixLineReceipt` records each entry's tail.
  `cup_junction` + `CupState` is the D18 mouthpiece-cup lumped
  compliance, TRAPEZOIDAL so the runtime matches the analytic shunt
  `Z/(1 + i omega C Z)` to `O((omega dt)^2)`. Realization doctrine
  (the recorded bake-off): at a brass input plane the higher modes
  are deeply evanescent, so the MM-derived `R_00` — which already
  contains every interior mode-conversion round trip — is what a
  plane-driven source sees; the receipt (ml-002) measured the matrix
  arm within 1 cent of the dominant line with a plane-mode energy
  fraction above 0.999 at 9x the cost.

- `brass_loop` module (bead 3ez8g.4.3): the composed brass playing
  loop — no `Trumpet` type. `BrassVoice` = a 1-DOF
  `fs_phs::mass_spring_damper` lip island (outward-striking; parameters
  at reduce-lab card scale with provenance labels on every input) x the
  per-valve multimodal characteristic lines x the D18 cup junction.
  PITCH IS NEVER ASSIGNED: the control surface is blowing pressure,
  lip tension, valve combination — structurally enforced by the
  battery's source grep. THE NUMERICS LESSON (executed): an explicit
  Bernoulli flow computed from the previous sample's bore pressure
  carries a one-sample delay whose phase flips the aperture's
  small-signal conductance ANTI-dissipative above fs/4 — a crook combo
  locked onto a 7.3 kHz parasitic scream; the fix is the closed-form
  IMPLICIT flow/junction solve (the linear cup relation substituted
  into the Bernoulli law gives a quadratic in the flow), leaving no
  delay in the flow path while the lip gap stays explicit (lip mass
  filters it). Passivity on the line taps is additionally enforced on
  a 4x-oversampled DTFT grid (inter-bin Gibbs overshoot makes a
  grid-only |R| <= 1 claim insufficient). Per-block diagnostics carry
  the lock estimate, gap statistics, blow/bore work, and the lip
  island's supply defect.

- `flue_loop` module (bead 3ez8g.10.3): the composed flue playing loop
  — no `Flute` type. `FlueVoice` = the jet-drive island x the wind
  line. NO CARD -> REFUSE: the island is parameterized only from a jet
  card minted by fs-aeroac's lab (`JetIsland::from_card`), and only a
  card whose claim class can drive a tone is admitted (the tonal
  interim card; the broadband-refusal-boundary card refuses with
  `CardClassCannotDrive`). Card-backed numbers: convection ratio
  `c_p/U = St (h/delta) / (stage + 1/4)` from the stage-I lock, the
  profile ratio `theta/b`, the Reynolds validity band, the claim class;
  authored and labeled: the spatial deflection gain, the onset seed,
  the legacy mouth radiation resistance (unless physical reaction is selected).
  PITCH IS NEVER ASSIGNED: blowing
  pressure and geometry gestures (cut-up, labium offset, jet angle)
  are the controls. Model: Bernoulli jet speed; receptivity = the
  leaky INTEGRAL of the delayed mouth velocity (its 1/omega cancels the
  drive derivative's omega, so the pipe's resonance alone selects the
  mode); the sinuous mode's amplification BAND from the card's theta
  (two low-pass sections at Michalke's tanh neutral `omega_c = U/(4
  theta)` and a high-pass at `omega_c/2` because growth vanishes at
  low Strouhal); labium cut `Q_j = Q_max (1 + tanh((eta - y0)/b))/2`;
  series drive source `rho (Delta_m/S_m) dQ_j/dt` in an implicit mouth
  loop with the mouth inertance, radiation resistance, and the line;
  an EXACT per-block energy ledger (source = stored + numerical +
  radiation + pipe; closure ~1e-15). THREE EXECUTED TRAPS on the way
  to a note: (1) a frequency-flat deflection with the omega of the
  drive derivative locks the highest stage the delay allows (9.4 kHz);
  (2) a continuous onset seed at 1e-3 U, integrated, flaps the jet by
  its own half-width and the labium slope averages to nothing (no
  lock; the seed-independence falsifier in the fixture is the guard:
  a limit cycle's RMS does not move when the seed drops tenfold);
  (3) an amplification band with no low edge lets the integral gain
  win at 6 Hz where the open pipe has no impedance (a relaxation cycle,
  not a note). Historical observations on the flute-class fixture (6 mm bore, 12 kHz
  line for the radiation law's `ka < 1`, 10 x 1 mm flue, 4 mm cut-up,
  tonal interim card): 10-50 Pa lock the fundamental (0.90-1.01 f1),
  110-240 Pa the octave (1.85-2.02 f1), 360 Pa the twelfth (2.84 f1);
  a 5 mm cut-up locks one mode below a 3 mm one at equal pressure; four
  rank voices are bitwise their solo selves in any stepping order;
  replay is bitwise and keyed by Philox voice identity. Those register-pressure
  observations predate the retarded-time correction below and do not qualify
  its corrected feedback phase; they require a fresh sweep. The card's
  Reynolds band (144-264) does not cover playing pressures (Re ~1500),
  so `claim_check` refuses every claim from this fixture: that boundary
  is inherited by the flue-gates bead, not papered over. This chain is
  the DEFAULT, and its gain is AUTHORED: the cards' gain 90 is about
  7-18x the forced jet's response at each note's speaking mode (mail
  7491). It stands until the forced jet passes fg_003, fg_004, fj_001 and
  gas_direction_amount (ruling 7518, mail 7859, bead pexfp.135.1).

  OPT-IN FORCED JET (`FlueVoice::with_forced_jet`, card row `jet_response
  forced-bickley SPACING WINDOW`; rulings 7376/7518): the forced Bickley
  jet's complete causal response replaces the band chain and the authored
  gain. The labium displacement is `eta_n = sum_k g_k v_(n-k)` of the exit
  transverse velocity (mouth velocity plus onset seed). `g_k` are the hat
  weights of one step response integrated by `fs_aeroac::forced_jet` at
  the station `cut-up / b` and the Reynolds number of the held pressure
  (the gesture's highest blowing pressure, including a declared source's
  pressure path). The exit is the Poiseuille slit (`b = 0.4 h`, `U0 =
  1.5 U`) behind a slit windway attached first; otherwise it is a uniform
  lip jet (`b = 0.75 h`, `U0 = U`). At variable speed, accepted velocity
  intervals retain their own advances of `tau = integral U0/b dt`; a new
  speed never retroactively retimes the history. The explicit fixed-Re
  extension is `dz/dt = (U0(t)/b) A z + B v`, with held interval speed and
  linear physical-time velocity. Local derivative-knot integration avoids
  cancellation of accumulated response integrals. The current velocity's
  affine coefficient still enters the same-sample mouth solve; the response
  is held after its window. Constant positive speed from rest retains the
  old cached hat-weight arithmetic. Zero speed disables input/output and
  freezes this approximate clock, with no claim of molecular decay at rest.
  Trial preparation never commits the clock or velocity; accepted publication
  alone appends an interval. The declared minimum positive speed bounds
  the same number of retained samples (three scalars per interval).
  G0/G1/G3 source tests in `flue_loop::forced_jet::tests` cover original
  uniform arithmetic at three clocks, independent response-integral
  correspondence after speed changes/ring wrap, the instantaneous-retiming
  falsifier, and transactional rest/refusal behavior.
  Register physics is graded on the
  production row 0.1 b and 90 b/U0 (`PRODUCTION_SPEC`). Numerical fixtures
  declare `COARSE_SPEC` 0.25/40, which renders the whistle D5 card within
  0.3 cents and 0.3 dB of the production row. A process memo retains 32
  integrated responses keyed by their exact inputs; the integration is
  deterministic, so reuse is bit-identical.
  REFUSALS: current-sample feedback beyond the solve's uniqueness bound
  asks for refinement (a jet that crosses the cut-up and grows within one
  sample; at 48 kHz over flue 0.1-1 mm, cut-up 0.06-0.8 mm and
  100-4000 Pa, only the 5.3 b and 10.7 b stations at 4 kPa reach it); a
  cut-up shorter than one grid cell of the response refuses at attach; a
  retained history at the declared minimum speed beyond the declared budget
  refuses; so do a non-positive reference pressure or one whose jet is below
  that minimum, the unsteady windway
  (its jet starts from rest, where this frozen-time response needs
  unbounded history), a windway attached after it, and the Nolle/Price
  transfer beside it.
  NOT CLAIMED: the base jet's own spreading over the cut-up, the labium's
  hydrodynamic back-reaction, the Reynolds number of attack and release
  (frozen at the held pressure), and memory beyond the window.
  Of the four gates of the switch, its fg_004 twin
  (`fg_004_forced_jet_*`) passes unchanged, but it LANDS HIGH on the three
  register gates. Their ignored forced-jet twins
  (`fj_001_forced_jet_*`, `fg_003_forced_jet_*`,
  `gas_direction_amount_forced_jet.rs`) carry the measured registers. The
  default switches only when all four pass with their assertions unchanged.

  JET SPEED EDGE (ruling 8132, bead pexfp.135.2): the declared minimum jet
  speed U_min (`jet_limits`, or 1 m/s on the legacy TMM image) sizes the
  retained history, so every lookup at or above it is retained. A slower
  positive jet is admitted and its delay is never clamped. In the onset,
  before the jet first reaches U_min, a lookup beyond the history that
  precedes the voice's first motion reads the rest state. That is exactly
  zero at any delay: unbounded history at zero memory, bit for bit. The
  forced jet's response reads input older than its first accepted interval
  the same way. Its intervals keep their own advances, so a release usually
  still reaches the window; once a jet below U_min fills the ring with
  intervals that do not, the forced clock freezes as at rest and the voice
  takes the declared state. Any other lookup beyond the history below
  U_min (a release, a dip mid-note, an onset slower than its own transit)
  takes the declared state: no deflection, so the windway's steady flow
  splits at zero deflection, `Q = q_max/2 (1 + tanh(-y0/b))`, and its ramp
  source stays.
  `FlueBlockDiag::below_speed_edge_samples` counts those samples. U_min is
  a declared validity edge, not a stability limit: below it the jet carries
  at most U_min/U of the held flow, and the parallel Bickley jet stays
  unstable down to Re ~ 4 (Tatsumi & Kakutani 1958). Zero blowing keeps its
  own path; the Nolle/Price transfer and the unsteady windway keep their
  own histories.

  Quasi-steady convection reads mouth velocity at the actual retarded time
  `t_n - tau`: the newest committed history entry is `t_(n-1)`, not `t_n`.
  For transit shorter than one sample, the current-velocity interpolation
  term is solved simultaneously with the jet cut-flow and acoustic mouth,
  including any selected contraction loss. Trial evaluations preserve the
  filter memories and consume no extra onset-noise draws. A sufficient global
  uniqueness bound is required before committing the mouth wave; unresolved
  feedback refuses rather than acquiring an artificial extra sample of delay.
  The interpolation regression establishes second-order retarded interpolation,
  not second-order accuracy of the complete jet integrator or measured flute
  fidelity. Gas and windway material changes cross the one-sample boundary
  through the same physical transit calculation.

  Optional `with_viscous_windway(length_m)` binds the same material-gas
  `ViscousSlit` used by reeds to the actual flue exit width/height. Its solved
  bulk speed reaches volume injection, Reynolds number, convective delay,
  amplification and the physical outlet receiver. Schedule admission and
  accepted pressure controls check that actual speed against history/Mach
  bounds; cold attachment also rechecks a previously held pressure. Omission
  preserves the ideal jet. `flue_windway` tests cover viscosity-only material
  swaps at unchanged density/sound speed, actual microphone pressure, control
  refusals and callback partition invariance with simultaneous mouth loss.
  The supplied rectangular channel is a quasi-steady parallel-plate reduction,
  retaining the card profile; exit kinetic head is not acoustic mouth loss.
  No entrance/profile development, windway inertia, turbulent jet, transverse
  flute lip geometry or measured instrument fidelity follows from this option.

  `flue_loop::network::performance::file` also admits explicit branched tube
  graphs, section-local material gas records, addressed viscothermal sections,
  passive endpoint/interior R-L-C loads, and multiple coplanar radiating
  openings. The jet and windway bind the actual inlet section's complete gas,
  independent of source order. The same physical graph, fractional/loss owner
  and multirate receiver retain the original node addresses through loss
  subdivision. CLI provenance describes every supplied node, section gas and
  selected opening. The single-tube grammar retains its original construction.
  `music_render_wind::flue::network` supplies G0/G3 source/direct-construction,
  gas/branch/receiver sensitivity, callback/PCM and refusal regressions. This
  adds data access to existing physics; it does not add finger-pad mechanics,
  junction end corrections, mutual radiation feedback, gas evolution or
  measured-instrument authority. See `examples/FLUE_NETWORK_PERFORMANCES.md`.

  `FlueOutletPerformance::with_mouth_reaction` selects compact self reaction of
  the actual rectangular mouth in the same exterior gas and area basis as its
  pressure receiver. Optional `mouth_aperture rectangle WIDTH_M LENGTH_M`
  before `mouth_radiation` supplies acoustic dimensions independently of the
  jet's width and cut-up; omission retains the original rectangular window.
  Network observation and compact reaction share this declared aperture.
  Independent aperture dimensions are not yet supported by the DG route.
  This choice does not alter the aerodynamic jet geometry,
  source coefficient or independently selected separation law. The exact
  rectangular static Rayleigh integral and compact
  monopole resistance determine one positive-real `R*s/(s+p)` branch in the
  existing `fs-vfit` midpoint owner. Current-flow pressure and retained history
  enter the original same-sample jet/mouth/duct solve. Original windway, jet
  feedback, gas sources and separately declared contraction loss remain. The
  branch replaces the legacy end correction and authored linear resistance;
  near-field storage and escaped radiation energy have separate ledger terms.
  The file selection is `mouth_reaction compact-baffled` immediately after
  `mouth_radiation`, with zero `mouth_loss_rel`.

  Construction bounds cell work, requires matched static pressure and a compact
  band, checks static quadrature against the exact rectangle, and checks both
  analogue and realized complex impedance and resistance within 5% at 32
  frequencies. These sampled checks are not a continuum or broadband
  certificate. Mouth-to-terminal mutual reaction, free-field body scattering,
  resolved turbulence and measured flute qualification remain outside this
  selection. Physical storage/work and source-to-PCM regressions are authored;
  syntax checks pass, while native compilation and execution remain required.

- `wind_line` module (bead 3ez8g.6.3): the wind articulation runtime —
  `WindLineBank` holds per-fingering exact-FIR characteristic lines
  over a typed `FingeringTable` with the carry-history + crossfade
  switch (the MM bank's D17 lift, plane-wave edition), and the
  char -> VFIT hop: a settled note may leave the exact-FIR phrase
  image for the cheaper rational hold image via the REPLAY-PRIME lift
  (the incoming IIR line is fed the recent outgoing history, outputs
  discarded, before the crossfade), hopping back on gesture
  resumption. Hop-readiness is GATED: entering an image whose
  registry gate is not green refuses structurally
  (`ImageNotGated`). "When to hop" is MEASURED DATA
  (`data/claims/wind-hop-policy.tsv`: settle-detector parameters +
  click-vs-hop-timing table, minted by the fixture, gated on every
  run), not doctrine. Fundamental estimation in the fixtures is
  AUTOCORRELATION — three estimator traps were executed en route:
  a global FFT peak returns the dominant harmonic; even-factor HPS is
  structurally blind on odd-harmonic (closed-pipe) spectra (the
  |X(2f)| factor lands on the absent even harmonic); odd-factor
  products cannot disambiguate f0 from 3f0 (the odd series is
  self-similar under 3x). The period domain has neither ambiguity.

- `piano_vertical` module (bead 3ez8g.5.1): the piano vertical
  composition — still no `Piano` type. Three unison strings (exact-ZOH
  modal images whose mode series carries the stiff-string
  inharmonicity law `f_n = n f0 sqrt(1 + B n^2)`; tension is state) +
  a duplex segment + a small modal board coupled through the
  POWER-CONJUGATE bridge (`F = sum c_k v_k` in, `g_k = -c_k v_b`
  back) + the 87zbd felt hammer island + pedals as coupling states
  (dampers = viscous drag through the strings' own force ports; una
  corda = struck-string count). The Weinreich aftersound EMERGES:
  in-phase unison motion pumps the lossy board (fast early decay);
  the surviving quasi-antisymmetric configuration barely couples
  (slow aftersound). MEASUREMENT LAWS learned executing the gates
  (each an executed misread first): the board-velocity spectrum
  carries the BOARD's own lines — the string law must be read from
  the string state; the board envelope is structurally BLIND to the
  antisymmetric survivors — the two-stage fit runs on the string
  ensemble energy; strong detune SUPPRESSES energy exchange (real
  physics) — the beat gate reads the AUDIO envelope, whose chorus
  modulation persists; and envelope smoothers must be scaled between
  the coupling slosh (~33 Hz, measured) and the detune beat (2-4 Hz).

- vibration gates (bead 3ez8g.7.1): the string/plate/nonlinear
  registry review. The TRUNCATION GATE makes modal honesty executable:
  a row's retained series must have converged (top-retained-mode energy
  share below the authored gate), and the deliberately under-truncated
  falsifier (N=3 at 220 Hz, whose top mode sits mid-band where the felt
  hammer still delivers ~1e-1 of the energy) FAILS the same gate the
  disclosed N=12 fixture passes at 6e-20. The linear<->nonlinear
  SELECTOR THRESHOLD is committed data with fixture provenance
  (`data/claims/vibration-selector-thresholds.tsv`), re-derived every
  run from the KC glide law the fs-nlmodal battery pins at 1e-12 — a
  hand-edited threshold fails. THE GONZALEZ CAVEAT travels with the
  KC/von-Karman promotion: conservation under discrete-gradient
  stepping is structurally blind to gradient errors, so nonlinear
  evidence cites FD-gradient/trajectory oracles, never energy
  conservation.

- Dispersive-waveguide bake-off (bead 3ez8g.7.2): the [F] string
  alternative measured against modal ZOH on ONE shared music-wire card
  with the analytic stiff law `f_n = n f0 sqrt(1 + B n^2)`,
  `B = pi^2 E I/(T L^2)` derived from the card, as the oracle. The
  executed receipt decided KeepBoth: modal keeps every certified claim
  (D21 asserted against the live registry in the gate test); the
  waveguide earns the budget-constrained hero-string scope
  (B_hat/B = 1.014, worst partial 9.8 cents with the linear
  fractional-interpolation phase as the named residual owner, ~10x
  fewer flops per sample, state cost flat in partial count — while
  modal holds FEWER states on this card, recorded honestly). The
  DISPERSIONLESS control on the same loop misses the law by 41 cents:
  the allpass stage, not the loop, carries the stiffness.

- Xylophone bar (bead 3ez8g.12.1): the cheapest FULL filling — a pure
  composition of live parts (analytic free-free bar modes + exact-ZOH
  runtime + a Hertz mallet island; no new physics, only binding and
  gates). THE BOUNDARY CONDITION IS LOAD-BEARING: the free-free
  `beta L` values are the `cosh cos = 1` roots computed in-fixture by
  deterministic Newton (self-verified at machine zero, never
  transcription-trust), giving the non-harmonic ladder
  f2/f1 = 2.7565, f3/f1 = 5.4039; the executed falsifier proves the
  gate discriminates (the PINNED n^2 family measures f2/f1 = 4.000 and
  fails). Material: the Indian-rosewood matdb pack (E_L = 13.53 GPa;
  density = SG*1000 as an authored Estimate, disclosed). The named
  successor is the ARCHED bar (ratio tuning toward 3-4, a
  varying-section solver) — v1 claims the uniform prismatic bar
  exactly and only that.

- Cane reed card (bead 3ez8g.6.2, completing T-Wind): "shaped like a
  reed, made of cane" — a tapered vamp mesh + the CC-BY cane packs
  through the fs-solid reduce lab (f1 = 2709 Hz, the measured
  clarinet-reed class) driving the EXISTING realize_reed_bore loop;
  the authored-parameter path remains the Estimate debug image. The
  LAY is a measured-geometry chart (tip opening + facing length ->
  circular-arc radius) consumed as the rest gap, rolling-stiffness
  hardening the named successor; contact stays the fs-dcontact
  Obstacle pattern. Provenance chain in the card: mesh digest,
  material resolution identity (heel E_L dry; the blank-zone density
  mismatch and the 10x E_T anisotropy DISCLOSED; eta = delta/pi from
  the damping pack), reduction options, retained-compliance
  disclosure. Measured sensitivities: more-open lay +46% threshold;
  harder cane brighter at matched breath (+1.9%) with the ordering
  inverting near the soft threshold (both regimes recorded); the
  stiffness threshold shift is below the ladder rung (logged, never
  overclaimed).

- `glottis` module (bead 3ez8g.8.2): the glottis as the SAME object
  as a lip reed and relief valve — two islands on the menu (D21): the
  1-DOF image with the TITZE SURFACE-WAVE mechanism (the flow sees the
  gap one mucosal propagation time late; a bare quasi-steady one-mass
  valve measured NO phonation to 6 kPa — the delay IS the mechanism;
  wave speed authored Estimate) and the two-mass image (4-state
  coupled pHS with the coupling INSIDE the Hamiltonian, corrected
  Steinecke-Herzel driving, asymmetric seed escaping the unstable
  symmetric rest, dcontact collision). FIVE energy/mechanism defects
  were EXECUTED en route and their causes are recorded in-code where
  each was fixed: the short low-loss tract's FIR wraps non-passively
  (4x-oversampled enforcement, the MM-bank lesson); the junction
  Picard iteration does not contract when Zc U rivals the drive
  (closed-form quadratic solve, the brass lesson); a stale one-sample
  line push adds junction phase slip; the S-H recovery factor's SIGN
  is load-bearing (inverted, the cycle locks antiphase onto the tract
  resonance); and the symmetric rest state is a driving-pressure
  deadlock (S-H's own displacement seed). Fold parameters:
  masses/stiffnesses/length RECEIPTED from the licensed two-mass pack;
  damping, rest gap, collision, mucosal speed authored Estimate — the
  honesty ladder's permitted rung, card-driven either way.

- `tract` module (bead 3ez8g.8.1): the tract as A(x,t) with LICENSED
  charts (Assaneo 2011 PLOS CC-BY Table 2 — the authors' own 10-tube
  /a/ and /u/; Spanish, acoustically-inverted, formant-valid register
  disclosed; Story 1996 is unlicensable and the Dresden CC0 meshes are
  the named MRI-grade upgrade vein). D22 as a WHOLE-CURVE gate: char
  impulse-spectrum envelope vs the TMM |Z_in| envelope of the SAME
  walled chart. STRUCTURAL FINDS the lane earned: characteristic_line
  consumers must pass VOLUME-normalized zc (rho c/S — an area-free
  rho*c flattens R to +1 and silently erases the tract shape; the
  glottis module carried this defect and its gates were re-based under
  the corrected load); a rigid-wall short tract realizes ACTIVELY (its
  ringing outlasts the FIR window and time-aliasing makes |R| > 1 —
  the tissue walls are required for the REALIZATION, not just the
  physics); a uniform passivity clip that fires flattens the shape
  (backstop only); dense realization via characteristic_line_dense
  (n_fft override). The morph carry needs p-minus CONTINUITY or every
  rebuild micro-clicks.

- Vowel gates (bead 3ez8g.8.3, completing T-Voice): the first voice
  claims — the vowel from the licensed area function, the source from
  the self-oscillating folds, their INTERACTION through the load at
  the glottis. The corpus loop closes (rows
  acoustic-assaneo-2011-tract-{a,u}): the charts' TMM F1 lands in the
  published Spanish vowel CLASS bands (827 / 340 Hz; a class check,
  never a per-subject match) and the length-shifted falsifier fails
  the same gate at 1104 Hz. THE FLAGSHIP EMERGENT GATE: sweeping the
  fold frequency toward F1(/u/), the phonated f0 rides ~35% above the
  free-fold frequency (load stiffening) and REGISTER-BREAKS at the F1
  ceiling (deviation -147 Hz at the top of the sweep) — the physics a
  source-filter vocoder structurally cannot have, measured. The
  1-DOF-vs-two-mass bake-off receipt is EXECUTED (KeepBoth:
  cheap-ensemble vs solo-realism; slope -6.6 vs -11.5 dB/oct on the
  identical fixture); the ooh/aah ladder is minted for the owner's
  ear.

- `pickup` module (bead 3ez8g.9.2): the Faraday pickup as a
  PURE-READER PORT — D20 as structure (pose + per-mode window gains
  only; `emf()` reads the string owner's modal velocities; no step, no
  state, no re-discretization anywhere). The raised-cosine aperture
  window is unit-normalized so the point-pickup limit is EXACT; the
  window is what makes bridge-vs-neck voicing physical. v1 provenance:
  linearized flux with an authored dipole-like height profile
  (Estimate, disclosed); the MAGNETOSTATIC BAKE (offline B-map lab,
  the lips/jet authority pattern) is the recorded follow-up that
  upgrades provenance without changing the runtime shape.

- `speaker` module (bead 3ez8g.9.4): the Thiele-Small driver + sealed
  cabinet as ONE port-Hamiltonian — the Bl gyrator is an antisymmetric
  J pair (power-exact by construction), the sealed box is
  `Sd^2 rho c^2 / V` extra cone stiffness (exactly why the classic
  resonance-shift algebra holds), and one fs-plate-class panel mode
  couples through the shared cavity term (the cabinet sound as
  physics; X-Consist with ONE retained mode, disclosed). Excursion
  past the card's x_max REFUSES — the suspension nonlinearity is not
  modeled and silent linear forte is the lie the refusal prevents.
  Card = authored datasheet-class Estimate; low-ka piston load
  disclosed; cabinet IR packs are FORBIDDEN and absent.
  The opt-in `MagneticTsDriver` retains that original mechanical assembly and
  admits caller-supplied quadratic Bl(x)/L(x) profiles over an explicit closed
  displacement interval. L must remain positive throughout, including interior
  extrema; the interval must fit inside the suspension's original travel limit.
  Total linkage lambda is stored with H_m=(lambda-psi_pm(x))^2/(2L(x)), where
  psi_pm is the analytic integral of Bl. The same gradient supplies current,
  Bl(x)i + L'(x)i^2/2 cone force, and reciprocal voltage including L'(x)i*x_dot.
  No additional linear gyrator duplicates the interaction. The existing
  descriptor discrete-gradient solver owns the common step and work/heat ledger;
  rejected trials do not change accepted state. This nonlinear path is second
  order and needs timestep refinement; the linear TsDriver keeps ExactZoh.
  Tests in `speaker/magnetic/tests.rs` cover constitutive energy derivatives,
  independent voltage/power identities, finite-domain admission, constant-model
  convergence, nonlinear cabinet/panel evolution, and source/heat balance.
  These bounded profiles do not add data authority, current-dependent core
  saturation, hysteresis, eddy losses, coil heating or nonlinear suspension.

- Electric chain (bead 3ez8g.9.5): `tests/electric_chain.rs`
  composes string -> pickup -> tone RC -> triode island -> TS driver
  with every stage a named, receipted object. The composition laws:
  ONE string owner (a bare string stepped alongside stays BITWISE
  identical — the chain reads, never owns); voicing survives the
  chain (per-harmonic bridge/neck output ratios = pickup window-gain
  ratios on string-carried lines); the harmonic ladder is
  DEVICE-made (single-mode string, so 2f0/3f0 at the speaker is
  constitutive curvature); the box composes transparently
  (chain-measured excursion-per-volt = driver-alone). The
  split-vs-full-DAE bake-off receipt
  (`tests/receipts/electric-split-vs-full-dae.bakeoff`, KeepBoth)
  pins the split image to the analytic gain and to the monolithic
  1e-14 authority. MEASUREMENT LAW recorded: harmonic line
  projections need INTEGER-CYCLE windows — at a non-integer cycle
  count the fundamental leaks ~1.7% into every harmonic line,
  drive-independent, and read as a 27% violation of a ratio law
  that holds.

- Composed guitar body (bead 3ez8g.7.4): `tests/guitar_body_gate.rs`
  gates the braced-top + ladder-braced-back + Helmholtz model's
  BRIDGE-ADMITTANCE triad against the Carcagno 2018 BR-guitar corpus
  rows (F1 +5.9%, F2 +9.9% in 15% envelopes; F3 -32.6% in an
  authored 40% envelope with the mechanism DISCLOSED — the flat
  rectangle surrogate has no X-bracing/arching, so its third
  breathing family sits low: the named weakest claim). Executed
  falsifiers: the midline-brace order flip (frequency-order breathing
  selection picks a 6e-15-volume mode; select by max |INT phi dA|)
  and the MEASURED refutation of nodal-line folklore (an ECCENTRIC
  brace on the nodal line shifts the mode 11.6% via torsion +
  eccentric-membrane coupling). The wolf note is a pHS pair: tune
  EFFECTIVE frequencies (k = m w^2 - k_c; the coupling spring detunes
  a naive pair by 10% and no hybridization happens), then the string
  on the breather beats at DETRENDED depth 0.046 (the estimator
  lesson: divide out decay before measuring modulation) vs 0.990
  clean off-resonance. Truncation sensitivity is PER-FREQUENCY
  (spread 23.6x) and its direction is fixture-specific: here
  between-peaks suffers more than on-peak because low-index resonant
  modes survive basis halving.

## Conformance tests

`src/broadband_radiation.rs`: G0 complex/real-SH round trip/reconstruction against
`fs-bem`; admission, neutral replay, properness/direct-feedthrough retention, and
held-out mutation; causal transactional runtime superposition and bitwise replay.

`src/vibroacoustic.rs` unit tests: closed-form overlaps, two-oscillator
split + dropped-coupling mutation, added-mass first-order falsifier,
exact energy balance + one-sided normal-flip mutation, reciprocity,
truncation refusal, cavity/resonator pins (incl. the elongated-cavity
regression), multi-mode pencil vs independent determinant, refusals,
determinism. `tests/vibroacoustic_casebook.rs`: the fs-plate ->
fs-couple -> fs-bem composition — box-with-flexible-top vs independent
perturbation with the sealed-cavity stiffening direction, coupled FRFs
with observed per-frequency truncation deltas, and the radiation solve
whose per-period energies audit GREEN through a real
`WindowAuditReport` (with the `entropy_on_nonthermal_power_port`
refusal retained as doctrine: dissipation leaves through THERMAL
ports); plus the guitar T(1,1) comparison against Carcagno et al. 2018
with the exact two-DOF product invariant and the material-vs-measured
damping inequality.

`tests/couple.rs` (38 cases): v2 scalar-seed migration goldens for all three
legacy kinds and bitwise whole-result migration of the added-mass fixture;
rotational/electrical/magnetic/chemical watt-dimensional goldens, required-role
admission, raw-oracle refusal, and semantic kind mismatch refusal; schema
fail-closed metadata; localized junction mismatch; non-scalar refusal by the
scalar evaluator; field-duality measure dimensions, density-side assignment,
measure-application overflow, and distinct effort/flow spaces; all four
primitive descriptors and identity-alias refusals;
complete moving-enthalpy stream admission and canonical constituent ordering;
one proved explicit-chemical contribution;
bit-exact enthalpy/internal-energy pressure-work equivalence; one-ULP and
foreign-context crosswalk/stress-evidence refusal; wrong-identity and
non-positive-density refusal; all four exact Euler/Legendre coordinates
reconstructing the same canonical enthalpy and chemical double-count refusal;
proof-bound axis mismatch and species-mode/element-axis refusal;
empty/duplicate/non-finite/wrong-orientation stream refusals and
declared-energy mismatch;
green first-law/mass/element/charge/entropy closure with explicit signed
advected/diffusive/heat/radiation rows and `T0*S_gen`; retained red first-law
and negative-entropy reports; exact missing-row and foreign-evidence-role
refusal; exact manifest/element/charge evidence-subject refusal; same-interval
energy/accounting-context rebinding refusal; stream/chemical element-schema
and typed power-axis refusal, including cancelling forbidden entropy slots,
inexact constant-temperature energy, and constant-temperature radiation;
boundary-temperature and nonfinite fail-closed admission; all six
permutations of a cancellation-sensitive compensated sum; an explicitly closed
isolated entropy window; one-sided internal-pair refusal; and a hidden internal
entropy source localized to the exact two rows, ports, and exchange while an
external row masks its pair-local energy defect from the global first law;
legacy power-conjugate ports; exact interface power conservation and
incompatible-port refusal; energy-audit imbalance and non-finite alarms; the
Aitken Δ² factor; naive staggering diverges where Aitken stays stable; Aitken
accelerates over stable fixed relaxation; light added mass converges naively;
determinism.

`src/mm_line.rs` `mm_line_tests`, cases ml-001..ml-006 — held-fingering
line fidelity vs the MM authority (0.39 cents worst over three valve
combinations; the explicit scattering LOOP additionally carries its
one-sample junction skew, measured 38 cents at a ~127-sample period and
gated at a disclosed 45), the realization bake-off receipt, the valve
no-click gate (carry+crossfade within 3x the control's slew; the
hard-cut falsifier is at least 2x worse), per-line passivity on the
stored taps, the cup shunt vs its analytic transform, and
causality/splice/refusal arms (the causality gate is the executable
witness of the e^{-i omega t} conjugation discipline).

`src/brass_loop.rs` `brass_gate_tests` + `brass_observer_tests`, cases
bl-001..bl-008 — the EMERGENT gates: lock 238 cents above the
cup-loaded column peak (outward-striking band) and 281 cents from the
lip alone at 10.4 kPa; tension walks the lock 208.9 -> 445.1 Hz across
distinct column peaks; valves drop the whole series (-125/-200 cents,
monotone, no pitch anywhere in the inputs); slotting (+78 cents of
continuous lipping inside a slot, 1310-cent jumps between); the cold
horn plays -50.4 cents flat vs the -48.6 sound-speed prediction with
temperature entering ONLY through the gas card; the no-click valve
change (envelope held through the switch; the fresh-voice hard-cut
falsifier collapses to 40% RMS); structural honesty + refusal arms;
and the radiated-observer demonstration from the SAME bake's SH
directivity table at the lock frequency (compact-mouth near-omni is
the correct physics at ka 0.09).

`src/wind_line.rs` `wind_articulation_tests` + `wind_hop_tests`,
cases wa-001..wa-004 + the committed-policy gate — the sigma slur
(+265 cents with the envelope held through the change), the register
vent (MEASURED, not the textbook sketch: with an INERTIVE 1.8 mm /
6 mm-chimney vent the twelfth-class flip happens at 0.70L, x3.36,
18 cents from the TMM authority's own peak, while the 0.30L vent
leaves the lock in the low regime — position matters AND the lock
follows the impedance peaks), the structural gate refusal, and the
hop policy (settled hop click 0.504 vs the early-hop falsifier 0.864;
the note survives on the hold image; the committed policy artifact
re-gated every run with the measured nuance that the settle detector
is conservatively late — the right side to err).

`src/piano_vertical.rs` `piano_vertical_tests`, cases pv-001..pv-005 —
THE TILT CONTRAST on the composed vertical (felt trend 159.7 vs the
matched linear spring's 1.000: the hysteresis is audible, not
asserted); the inharmonicity law measured in the audio (worst 4.23
cents over partials {1..5,7} with the B=0 falsifier missing by 15.04;
partial 6 disclosed as the duplex collision line, partial 8 as
weakly-excited at the strike point); the two-stage decay (string
ensemble 3.28/s early vs 1.15/s late) with audio-envelope beats
scaling with detune (0.370 s -> 0.239 s at doubled detune); pedal
topology (dampers collapse the note; una corda energy ratio 0.67 with
the third string still fed through the bridge); and the dissipation
ledger + bitwise replay (worst per-window growth 0.0).

`tests/vibration_gates.rs`, cases vg-001..vg-003 — the truncation gate
with its executed under-truncation falsifier (6.2e-20 vs 1.0e-1 across
a 2e-3 gate), the selector-threshold re-derivation pin (relative 1e-8
against the committed artifact), and the gate-summary enumeration
(every vibration-filling registry row must appear in
`data/claims/vibration-gate-summary.tsv` with a status MATCHING the
registry; the plate 1-port's honestly-ungated row is asserted by
name).

`tests/bakeoff_dispersion.rs` — the committed-receipt gate (bytes
decode, hash, KeepBoth outcome, D21 modal-row-untouched against the
live registry, the newcomer's own row present) and the
passivity/control gate (loop bounded and decaying under g < 1 with the
unit-magnitude stage; constructor refusals; the dispersionless
control's partial 8 sits within 20% of harmonic while the law demands
+18 cents).

`tests/xylophone_bar.rs`, cases xb-001..xb-004 — root self-verification
(residual at machine zero, classical table at 1e-8, asymptote check) +
named refusals; the rendered ratio gate within 1% with the executed
pinned-family falsifier; strike-position physics emergent (mode-1 node
found in-fixture at 0.2242 = the cord-mount point; striking it kills
the fundamental ~8 orders, the center strike kills partial 2, and the
non-noded partials keep speaking); and the listening chain gate
(WAV -> sidecar -> receipt, Unadjudicated).

`tests/cane_reed_card.rs`, cases cr-001..cr-003 — the provenance-chain
card mint (identity reproducibility, retained-compliance floor, the
reed-class f1 band, the lay's facing-scale radius); the sensitivity
sweeps with a SUSTAINED-oscillation speak detector (tail at 2% of
blowing pressure and non-decaying — the absolute-threshold detector
fired on ringing at every rung, twice, before this); and refusals by
name (a card without a modulus cannot resolve; lay-chart nonsense;
an empty orifice face).

`src/glottis.rs` `glottis_tests`, cases gl-001..gl-005 — both islands
self-oscillate against the live tract with onset curves in the
physiological range (1-DOF 188-276 Pa, two-mass 404-865 Pa, rising
with stiffness, none on the ladder floor — the detector demands
DETRENDED sustained output after a settling ramp read as phonation on
the first run); the collision/damping ledger (undriven kick decays,
per-step dissipation non-negative); the two-mass vertical phase
(upper lags lower 0.306 rad at the fundamental) with faster-closing
flow; the bake-off FIXTURE QoIs on the identical tract + card (f0
152/173 Hz, OQ 0.78/0.82, slope -6.6/-10.0 dB/oct — the receipt is
vowel-gates' to mint); refusals by name.

`src/tract.rs` `tract_tests`, cases vt-001..vt-004 — the whole-curve
D22 gate (envelope correlation 0.87//a/ and 0.90//u/ with dominant
resonances within 50 Hz; identity /a/ >> /u/); the yielding-wall F1
raise (340 -> 401 Hz, nearest-peak comparison past the wall's own
~100 Hz mode); the u->a morph (100 Hz control rate, carried, clicks
18x less than the five-cut waypoint control on the REFLECTED wave,
F1 track 353 -> 806 Hz logged); licensing-first refusals (unlicensed
data refuses by name). The glottis module's gates were RE-MEASURED
under the corrected load in the same change (two-mass onsets
404-715 Pa physiological and monotone; the 1-DOF's ~5.8 kPa
non-monotone onsets recorded as the fixed-delay mechanism's real
behavior; flow skew logged, not asserted).

`tests/vowel_gates.rs`, cases vw-001..vw-004 — the formant class gate
with its executed length-shift falsifier (a uniform AREA scale would
not move formants; the falsifier must shift LENGTH); the sung
interaction sweep (ten stiffness points crossing F1 with per-point
f_free/f0/deviation logging); the committed bake-off receipt gate
(KeepBoth scopes + both rows alive in the registry); and the
vowel-ladder listening chain (Unadjudicated).

`src/pickup.rs` `pickup_tests`, cases pk-001..pk-005 — station
voicing with every per-mode EMF line matching an INDEPENDENT
Simpson-integrated analytic gain within 5% (bridge centroid 452 Hz vs
neck 328 Hz); two mixed pickups combing per the gain-sum prediction
with an executed notch; the exact point-pickup limit (5e-4) plus the
pure-reader law (identical repeated reads; EMF recomputable from the
public gains and the same state slice); refusals by name; and the EMF
driving the .9.1 circuit DAE end-to-end (RC-load rolloff matching the
first-order analytic response at 0.25x and 4x cutoff with supply
defects at 1e-21).

`src/speaker.rs` `speaker_tests`, cases ts-001..ts-005 — the sealed
box shifts resonance per fc = fs sqrt(1+Vas/Vb) (measured 64.7 vs
65.1 Hz analytic); the electrical impedance matches the analytic TS
motional curve to the printed digit at four frequencies; the panel
colors the response near its resonance and NOT far away (executed
contrast against the rigid box); a 60 V resonance drive trips the
excursion refusal while a gentle drive never does; refusals by name
plus the supply-defect ledger at 3e-19.

`tests/electric_chain.rs` (7 + 2 mints): ec-001 compose + logging +
the one-owner witness; ec-002 voicing ratios within 6% (gain-floored
to string-carried lines); ec-003 H2 0.07% -> 0.86% across grid drive
0.03 -> 0.40 V RMS with H2 >> H3; ec-004 composed-vs-direct box
transfer within 0.04% at 65/260 Hz; ec-006 committed bake-off
receipt re-asserted (split gain 59.950 = authority 59.950, analytic
59.923; pluck L2 4.6e-12); ec-008 listening digest chain
(electric-clean-vs-driven, Unadjudicated); ec-009 the electric
gate-summary artifact enumerates every registry row.

`tests/guitar_body_gate.rs` (5 + 1 mint): gt-001 triad + product
invariant (defect 0.0) + damping honesty (material-only Q 39/54/78
above measured 34/18/36 as it must be); gt-002 the executed
midline-brace negative control; gt-003 wolf note; gt-004 truncation
spread; gt-006 guitar-strum listening digest chain (Unadjudicated).

## No-claim boundaries

- TOMBSTONE-3ez8g.13.5: `traveling_wave_line` is dead code with zero
  in-tree consumers, superseded by the exact-FIR characteristic line;
  every public item is deprecated + doc(hidden) and the module header
  denies new consumers. No music image may cite it (standing law).
  Deletion stays available on explicit written owner permission.


- The electric chain claims NOTHING about the string (one owner; the
  vibration epic's rows carry it) and nothing about perceived tone
  until the listening receipt is adjudicated. The preamp input gain
  is a measured-calibration volume knob (disclosed, not physical);
  budget rows wait on the budget lane covering this chain (D25);
  live-default stays no.

- `acoustic_realize` is a composed description→waveform, not a named
  instrument product and not a 3D jet. Bow is regularized friction, not
  measured rosin or Helmholtz-guaranteed. A linear `ThinPlate` is DKT +
  certified modes + compact monopoles. `geometric_nonlinearity` selects
  the von Karman modal pHS (SS sine or FE-sampled DKT `w`), not a
  full shell or an in-process BEM body.
  String+duct+plate share one sample clock and, when `EA > 0`, a
  Gonzalez string whose contact lives in `H`. The three-pHS Dirac
  join lives in `fs-phs`; this realize path still uses the FIR
  duct and one-way bridge force. Isolated blow is an impedance
  FIR, not a reflectance FIR. 3D broadband jet noise remains a
  no-claim. WAV encoding is a physical-scale dump, not mastering.
- `vibroacoustic` is frequency-domain only: time-domain realization is
  the vector-fitting bead's scope; non-rectangular cavity bases beyond
  the lumped Helmholtz mode arrive through the same `CavityModes`
  carrier (numeric Laplacian producers are not shipped here); no
  modal-density/SEA regime; the casebook's box fixture is validated by
  independent perturbation and exact power audits; the guitar fixture
  compares the coupled T(1,1) pair against Carcagno et al. 2018
  (CC-BY, Table I) through a DECLARED shape-surrogate top (deviations
  +10-17% recorded, product invariant exact) — measured mobility
  MAGNITUDES exist only as figures in the license-compatible
  literature and are not compared; `condition`-free dense solves inherit `fs_la::eigen_complex`
  boundaries.

- The FSI fixture is the classic LINEARIZED added-mass interface map
  (`H(x) = −μx + c`) — enough to reproduce the instability and its fix; a full
  nonlinear FSI solve over real fluid/structure subsystems is the consumer.
- `AitkenRelaxation` is the scalar Δ² relaxer; the vector INTERFACE
  QUASI-NEWTON (IQN-ILS) accelerator and MULTIRATE co-simulation are staged.
- PR-2 admits the rotational/electrical/magnetic/chemical vocabulary and
  kind-specific generalized dimensions. It does not by itself prove causality,
  DAE index, source closure, a constitutive law, or a physical port adapter.
- `ChemicalPotentialAmountFlow` describes species electrochemical potential
  and species amount flow. Reaction affinity/extent rate requires an explicit
  stoichiometric crosswalk and is not represented by that kind.
- The scalar evaluator does not execute vector/tensor/field or general
  multi-port Stokes–Dirac operators. `StorageElement` and
  `DissipativeRelation` carry durable public operator/evidence references; this
  crate does not execute or validate the referenced domain law in PR-1.
- Field-duality admission checks dimension arithmetic including the declared
  measure; it does not prove quadrature, trace pullback, orientation, or the
  numerical duality operator.
- `AlongFrame`/`AgainstFrame` can describe schemas and audit boundaries, but
  PR-1's conservative junction refuses them rather than assuming an unproved
  effort/flow pullback. A later neutral transfer API must make that transform
  explicit before admission.
- PR-3 admits `StreamPort`, one selected energy-accounting chart, exact numeric
  crosswalk identities, and context-bound durable proof references. It does
  not execute the referenced verifier, equation of state, stress operator,
  constituent-map artifact, or chemical partition proof.
- PR-3 freezes a stationary or co-moving accounting surface. It does not model
  boundary velocity, Reynolds-transport terms for a moving/deforming control
  boundary, or cross-frame kinetic-energy transforms.
- PR-4 audits direct mass, exact caller-projected element amounts, direct or
  externally projected coulombs, first-law energy, boundary entropy, internal
  pair closure, and optional `T0*S_gen`. It does not audit momentum, species
  production/reaction progress, mass-versus-molar consistency, or derive
  element/charge values from stoichiometry, valence, molar mass, or Faraday's
  constant.
- `fs-couple` checks evidence interval, semantic role, and subject identity but
  does not execute the referenced verifier or discover omitted sources. A green
  report means the supplied, closure-evidenced scalar algebra closes; it does
  not independently prove quadrature, physical-time integration, source
  discovery, EOS validity, species-to-element/charge projection, or the truth
  of the retained evidence statement.
- Logical clock ticks are ordering coordinates, not seconds. Every window row
  must arrive pre-integrated in physical time with an integrated-transfer
  receipt; endpoint sampling is not promoted into an integral.
- Constant-temperature thermal power rows check `Q = T*DeltaS_heat` exactly
  and therefore refuse radiation; radiative thermal rows require the profile
  evidence lane. Profile temperature, stream heat/radiation, and
  nonequilibrium entropy decompositions are retained external-verifier claims;
  this crate does not execute their quadrature or infer an average temperature.
- PR-4 currently admits a conservative homogeneous stream-context window:
  external stream rows must share state schema, constituent basis/axis,
  chemical reference, spatial coordinate binding, clock, gravity datum, and
  stress convention. Heterogeneous inlet/outlet contexts require an explicit
  upstream projection/crosswalk into one common audit context; this PR does not
  define or execute that crosswalk.
- Included-source entropy means entropy transferred into the audited system.
  Internally generated entropy must remain in `S_gen`; supplying it as an input
  would subtract it and is outside the admitted convention.
- The exergy ledger is only Gouy-Stodola destruction relative to the named
  environment, not a complete open-system availability-flow balance.
- Neither `WindowAuditReport::is_green`, `StreamPort`, nor an
  `AccountingBoundary` is a subsystem, discretization, co-simulation, or
  whole-machine passivity certificate.
- The two `fs-iface::SpaceType` entries record neutral effort/flow field roles;
  they do not supply a mortar/Nitsche/harmonic transfer operator, certify that
  the declared pair is dual, or prove inf-sup compatibility for a particular
  domain adapter.
- The energy audit's balances are supplied by the caller each exchange; wiring
  them onto the ledger is the coupling driver's integration.
- Dirac interface losslessness does not establish passivity of component
  storage/dissipation/source laws, spatial or temporal discretization,
  interface transfer, nonlinear iteration, multirate windows, or the coupled
  system. Those obligations require a signed, closed-window energy audit.
- The multimodal line inherits the fs-duct modal image's disclosed
  boundaries (matched-mouth higher modes, no tone holes, no mean flow)
  and adds its own: the explicit scattering junction carries a
  one-sample loop skew shared with the reed voice (disclosed, measured,
  not a line error); a hard (fade = 0) valve swap is genuinely
  discontinuous and is the caller's choice; the crossfade lift is an
  engineering lane, not a claim about real valve port geometry; and
  the realization cost is O(n_fft) MACs per sample per line (the
  matrix arm is N^2 of that) with no real-time claim until a budget
  row exists.
- The brass loop is Estimate end to end (authored card-scale lip; the
  claim decision is the brass-gates bead). Its own v1 boundaries: the
  lip pressure-collection face area is authored ABOVE the
  orifice-strip area (mouth-side lip face; disclosed in provenance);
  no WallPin lane through the multimodal image yet; the radiated
  observer is demonstrated per-direction at the lock frequency, not as
  a broadband convolution (the broadband lane is
  `broadband_radiation`); valve changes are discrete combos with the
  crossfade lift — a trombone's continuous slide needs the wind
  epic's fractional-delay bead, a dependency of claim, not of code.
- The wind articulation runtime inherits the char image's plane-wave
  scope; the hold-image realization here is a fixture-grade 12-pole
  fit (the certified hold lane is the clarinet casebook's); the
  register-vent finding is geometry-specific (an inertive chimney
  vent inverts the classic pressure-release node rule — recorded as a
  measured hypothesis, not a general law).
- The piano vertical's board is an authored spruce-scale modal set (the
  full fs-plate orthotropic chart with rib stiffeners is the
  piano-gates upgrade lane); bridge geometry is authored v1 (the
  ingest lane upgrades later); the explicit board coupling is
  semi-implicit at audio rate (the ledger gate bounds it); cabinet,
  room, and JCA linings stay out by the bead's own boundary. Claim
  promotion happens in piano-gates (.5.3), not here.
- The vibration gate review promotes claims on EXISTING evidence (the
  .1.4 bake-off receipt, the fs-plate Olson-Hazell battery, the
  fs-nlmodal analytic pins); it mints no new physics. The plate
  vfit-driving-point row stays ungated with the missing item named (no
  executed plate-1-port fit artifact); bowed and vibroacoustic rows
  carry evidence but await their own mechanism/corpus gates (.7.3-.7.5,
  .7.4).
- The glottal islands claim PHONATION MECHANISMS, not vowels: no
  formant filters, no vocoder, no glottal-waveform playback (the
  waveform emerges from the valve + tract loop or not at all). The
  uniform 17.5 cm tract is a fixture, not an articulation model
  (.8.1's lane); the bake-off receipt between the two islands is
  minted in vowel-gates (.8.3), not here.
- The tract module claims chart realization and morph mechanics, not
  vowels-with-a-glottis (.8.3 composes them); no fricative jet, no
  nasal branch (named refusals); mouth radiation is IdealOpen v1 (the
  compact-mouth Z bake is the recorded refinement — the radiation-fit
  ka ceiling forbids UnflangedOpen at full audio band).

### IRP radiated observation admission (M2.6-W0)

`pcm_wav::observation::ObservationPolicy` distinguishes `Acceptance` from
`Diagnostic`. `ReedPerformance::from_bytes_with_policy` refuses `listener_m`
under Acceptance; its physical outlet still requires the unchanged flanged
feedback load and uses one exterior gas for feedback and reception. The raw
voice, positional reed CLI fixture and legacy baseline remain diagnostics.
`music_render wind --acceptance` selects the flanged-development path explicitly.

`TwoSidedPressure::with_observation_policy` and
`AperturePerformance::new_with_observation_policy` require the actual installed
reaction binding under Acceptance. They check the load's identity, radius,
clock, band and exterior density/sound speed (relative tolerance 1e-12).
`from_baffled_environment` uses that checked path. Observation adds no second
reaction load, changes no bore gas and never changes a termination to admit it.
Existing one-way builders retain Diagnostic scope.

`TwoSidedAperture` binds the original `BoundaryApertures` on either physical
duct, including material-modal valves. Every driven node retains its complete
local gas, represented chimney and material/contact/Maxwell owner on the same
valve clock and horizon. Both graph/key previews and the whole-system work
gate complete before either graph, key, valve or supply state publishes; paired
driven acceptance retains MIMO, fractional, reactive and supply histories.
Storage/loss includes each graph and key once, with key-reservoir work and
mechanical actuation separate from series-pressure/reed work and signed MIMO
terminal exchange. `ScheduledTwoSidedPressure` reuses `BoundaryProgram` and
`BoundarySchedule` per side at mechanical sample addresses, including audio
substeps; both sides and source-force programs share the auxiliary control/work
budgets. A keyed side without an explicit program refuses before playback.
G0/G3 regressions exercise material substitution, complete work, atomic refusal
and cancellation, source/receiver parity and material-modal retention; native
execution remains required. The approximation remains a material cantilever
moving slit, not a resolved circular pad, key linkage, moving-opening mutual
exterior, or full M2.1-W4b qualification.

Native material-modal `boundary_groups` and mechanical-sample
`boundary_control` records bind these same keys on either or both sides.
The post-blade declarations use the existing shared `material_boundary`
compiler and the same bounded immutable cache as the reed plates. Explicit
held closing-force and exterior-reservoir assignments use the original
boundary scheduler at mechanical substeps. Existing typed force gestures
remain available; combining typed and explicit sample ownership refuses.
Keys derive their complete local gas, impedance, represented chimney, clock
and horizon from the realized source. Compact, imported-MIMO and
computed-baffled radiation prepare through their original owners before key
attachment and final receiver admission. Missing driven-node bindings,
duplicate physical identities and incomplete programs refuse. The native CLI
retains source records and controls; `paired-keyed-sampled.performance`
demonstrates both sides. Focused native/direct/PCM and material-substitution
regressions are included; their execution remains required. The moving
cantilever-slit approximation does not add resolved circular pads, linkages
or radiation from evolving key openings.

`TwoSidedAperture::with_mimo_environment` attaches the same source-bound
`MimoEnvironment` and original `WaveguideNetwork` terminal host used by the
flue to single, paired and material-modal valves. Each side admits one ordered
group before stepping; valve materials, full modal viscosity, regional gases,
fractional travel, distributed losses, material walls and supply states retain
their original owners. The terminal host returns mutually coupled pressures
inside the valve feedback loop. `TwoSidedFrame::coupled_terminal_work_j` sums
both sides' signed exterior work exactly once in the complete work admission;
it is neither dissipation nor an inferred rational stored-energy coordinate.

The two-sided microphone checks installed groups even under Diagnostic. It
preserves each side's port order while remapping global network 0/1 to the
source environment's local 0, then requires the exact source gas, actual outlet
geometry and mechanics clock. Omitted/reordered ports and mismatches refuse.
Output reconstruction/rate may change within the shared receiver's admission.
The supplied models are proper, undelayed, reciprocal and bound to the existing
coplanar infinite-baffle observer. Separate side groups do not represent
cross-side mutual radiation. The sampled passivity report remains GridOnly;
this path does not promote observation to Acceptance or close M2.5-W3c's
full-band bake, clock-refinement and closed-loop qualification requirements.

Native flute and material-modal files now select original radiation artifacts
through one bounded `multiport::input` decoder and their existing performers.
The artifact independently encodes its exterior, ordered aperture basis,
physical band and full pole/residue/direct matrices. Exact imported-byte hashes
remain distinct from declared geometry/material and realization identities;
import retains those declarations without attesting an offline bake. A
thermochemical exterior additionally binds the original canonical species,
NASA-9 and transport pack hashes and all twelve resolved gas-state fields by
bits. Replacing a registry pack under the same species id, changing a snapshot
or relabeling coefficients with the current caller gas cannot silently pass.
The existing runtime owns stability, reciprocity, properness, full source-clock
band, passivity, resource and geometry admission. Selected groups replace only
their observed memoryless endpoints. Unselected outlets/sides retain authored
loads, material/gas histories and losses; side groups have no cross-side
exterior coupling. Both source callbacks and CLI artifact reads are bounded,
legacy selectors are preserved, and failures precede performer/output
publication. This adds no fitting engine or instrument qualification. Grammar,
source scope and synthetic integration fixtures are in
`examples/MIMO_RADIATION_ARTIFACTS.md`; native execution of their new regression
sources remains required.

`multiport::bake::bake_baffled_mimo` connects the original `fs-bem` circular
piston/Rayleigh array to the shared `fs-vfit` proper MIMO fit and existing
terminal host. Its source is the actual ordered geometry and coherent exterior
gas, with force/velocity converted by both aperture areas and conjugated once
to the causal phasor convention. Native `computed-baffled` selectors on flute
and either material-modal side rebuild that source after data-only geometry or
atmosphere substitutions, without an imported matrix or a parallel runtime.
The full mechanical Nyquist band, bounded source sweeps, three explicit orders,
separate working/runtime memory, training and withheld matrix errors, source
causality assumptions, continuous passivity and interval-mean ZOH passivity are
admitted before attachment. All three unchanged fits and refusals are retained;
the smallest passing order is selected, without repair, clipping, damping or
gain changes. The gas record retains original thermochemical receipts.
Cancellation is polled between samples, causality and orders; the existing
individual quadrature/fitting/eigen kernels have no mid-call poll. This remains
a finite-band proper approximation of stationary coplanar uniform apertures in
an infinite rigid baffle and homogeneous inviscid gas. Sampled errors and
declared-tail causality consistency do not certify unsampled frequencies or
closed-loop accuracy; sampled passivity is GridOnly, not stored energy. No
delay extraction, measured-instrument claim or full M2.5-W3 qualification is
added. G0/G3 source/unit/gas and native integration regressions are included;
native execution remains required. An independent development-only SciPy
oracle found admissible order-four two-port fits at 8 and 48 kHz, but does not
substitute for the unavailable Rust checks.

G1/G2 analytic tests in `baffled_pressure_qualification.rs` check on-axis piston
pressure at z/a=2,10,100 and signed 2J1(x)/x directivity at ka=0.5,2,5 within
1e-3, with radial/angular refinement. The on-axis fixture uses a=25 mm to
respect the existing z>=50 mm domain; the far-field fixture uses a=7 mm.
The band parameter is admission metadata: its ladder checks invariant kernels
and live spatial/temporal refusals, not waveform convergence as the band changes.
G0 checks the finite-difference 0.1-rate and differentiating 0.45-rate caps at
cap and cap+1 Hz. G3/G5 adapter tests compare the same accepted terminal flow
against the shared owner and repeat pressure bitwise across callback sizes.
These are authored infinite-baffle fixtures, not real-instrument fidelity or a
full-band reaction match. Native qualification remains required; rustfmt alone
does not establish the numerical claims.

## Frequency-domain branched chart and network port (IRP M2.1-W8)

`tract::network::TractBranch` and `TractChart::network_with_branches` turn a
chart and explicit branch geometry/source identities into `fs_duct::network`
data. The resulting ports are glottis, oral and the supplied branch labels.
All propagation, junction equations and graph elimination stay in fs-duct;
there is no nasal-specific solver. The former no-nasal-branch documentation is
lifted for this frequency-domain API. TractVoice's time-domain series image is
unchanged; dynamic branch-line realization remains M2.2-W3c.

`driving_point::network::NetworkDrivingPoint` is the explicit FD port seam for
player-airway/bore consumers. Construction requires the graph, active label,
every other port's load, and a caller-supplied scalar `FluidPressureFlux`
PortSchema, with owner-outward orientation. Gas, omega and loss remain explicit
at each evaluation. `harmonic_port` performs the actual complex transfer
`p=Z U_in`, `U_out=-U_in` with that schema. It does not claim a time-window
energy audit, glottal dynamics, or a completed M3.1-W7 player-airway realization.

`NetworkDrivingPoint::realize_line` realizes that driving point as a
characteristic line. It uses the same ladder, Fourier convention and
4x-oversampled passivity gate as `realize_line`, and samples
`R = (Z - zc)/(Z + zc)` from the network solve. The nominal round trip is twice
the longest shortest axial path from the active port, computed with the
single-duct arithmetic. A network solve may refuse the realizer's
1e-3 rad/s DC probe as ill-conditioned. Below the largest rung's first bin,
that probe is answered from the lowest admitted decade. A fallback of
`R = -1` would ignore a lossy bore's Poiseuille resistance (`R(0+) = -0.945`
on a 2 mm bore) and leave a tail that only halves per rung. `new_coupled`
declares ports that one
`fs_duct::network::CoupledTermination` terminates together, for example a
multi-aperture radiation matrix. The caller supplies that matrix at every
sampled frequency, including the near-DC point, and a source is required
exactly when coupled ports are declared. `ReedBoreVoice::new_with_network` runs
the existing reed and compact-jet observer on that line. The active port's edge
radius sets the source impedance. This is a frozen-geometry FIR line. It is not
M2.2-W3c's DG network, a moving-geometry line, or a radiated far field from the
openings.

Tests in `driving_point::ladder_tests`:
- `network_line_matches_the_duct_line_for_one_edge`: same rung, taps within
  1e-9 of the duct line.
- `network_line_coupled_constant_load_is_the_siso_line`: bitwise equal, plus
  both declaration refusals.
- `network_line_keeps_the_branched_bore_resonances`: an open mid-bore chimney
  and a radiating bell. Line and network `|Z|` maxima agree within 0.05 cent.

`render::inlet_geometry_tests::reed_voice_on_a_one_edge_network_matches_the_duct_voice`
renders within 1e-6 of the duct voice.

### Closed-loop radiation audit (IRP M2.5-W3b G3)

`tests/closed_loop_radiation.rs::closed_loop_audit_reed_bore_radiation` closes
reed, bore and multiport radiation in one loop. Its fixture is declared
synthetic: a 7.5 mm bore with an open 4 mm chimney, whose bell and chimney
mouths end flush in one infinite baffle 45 mm apart. There the coplanar
baffled-piston matrix (fs-bem) is the exact `Z_rad` at every frequency.

The in-loop load is a passive, delay-extracted fs-vfit MIMO model built in
four steps:
- the exact bake passes the dispersion-relation causality test with its gap
  delays removed;
- it is fitted from 200 Hz to 8 kHz;
- grid passivity is repaired from 1 Hz, so it holds at every line bin;
- the repair's change is gated at 1e-4 over the fit band.

The loop runs at 16 kHz. Each check compares the model loop with the exact one:
- input-impedance peaks: within 1 cent and 0.5 dB;
- loaded modes and decay rates of the realized line below 0.5 Nyquist: within
  1 cent and 2 %;
- the reed's oscillation threshold: within 1 %;
- the playing frequency: within 1 cent;
- the harmonic spectrum up to 8 kHz: within 0.5 dB rms.

Measured (yto g3n):
- order-6 fit error 1.07e-7, repair change 1.3e-8;
- 15 peaks within 3e-6 cent and 6e-7 dB;
- 7 modes within 0.0025 cent and 7e-5 decay;
- threshold 2242.7 Pa identical; pitch 298.5 Hz within 1e-8 cent; 11
  harmonics within 3e-6 dB rms.

The exact loop is measured once, and every candidate is audited against it.
A planted negative must fail: independent mouths, the exact self terms with
the mutual terms removed. It is passive and exact on the diagonal, so only the
coupling separates it from the reference. It fails six of the seven gates:
- peaks shift 2.59 cent and 2.30 dB;
- modes shift 1.99 cent, with 23 % decay;
- the pitch shifts 1.24 cent;
- the harmonics differ by 1.42 dB rms.

The threshold alone does not see it (0.17 %). A NaN departure fails its gate.

No-claim: this is a well-fitted model on one fixture, with one planted
negative. Relaxed vector fitting collapsed at every order when the band
started at 10 Hz.

`tract::network::tests::network_branched_tract_nasal` uses Assaneo's licensed
/a/ chart plus a DECLARED SYNTHETIC 0.13 m x 4 mm-radius branch after section 2.
The mouths are explicitly ideal pressure release; this is a pole-zero diagnostic,
not a measured nose or an acceptance voice render. G3 checks the lower first
formant and added oral-transfer antiformant, with 8/4/2 Hz frequency grids
(last-pair feature positions within 4 Hz); all sampled N-port matrices must be
reciprocal and have PSD real Hermitian part to 1e-12 normalized tolerance.
Qualitative side-cavity context: Dang, Honda & Suzuki, JASA 96,2088 (1994),
<https://pubmed.ncbi.nlm.nih.gov/7963023/>. No human-specimen frequency oracle is
inferred from that citation. Tests print formants/antiformants and provenance.
Additional cross-crate tests exercise the explicit typed transfer, missing-load
and wrong-schema refusals, and unbranched graph/chain equivalence (1e-10).
All new numerical assertions remain pending central native verification.

## Assembly model-domain admission (IRP M6.1-W1b)

`description::validate(&AssemblyDescription, &ClaimDomain)` is the mandatory
admission seam for description-render consumers (M5.4-W2). Its private
`Admission` token borrows the exact immutable description and claim, exposes the
unchanged requested band, and records each selected registry row, module,
configuration, evaluated axis, refinement route and resolved material-state
identity. It does not render or narrow a band. Existing raw component render
APIs are unchanged; the description renderer must accept Admission directly.

`ClaimDomain` declares gas, dynamics, an existing `selector::RegistrySnapshot`
and explicit `ModuleDomain` selections by EntityId. This reuses
`fs_evidence::ValidityDomain`; every required or extra constrained axis must be
provided. No missing-axis-as-unconstrained rule applies here. A refused or absent
registry row cannot admit a model. Ungated rows remain visibly ungated in the
limitations; domain admission is not a gate promotion or Q3/Q4 acceptance.

Primitive parametric geometry vocabulary (all dimensioned SI; a
`test-fixture:` prefix marks authored fixtures without changing semantics):

- `cylinder.v1`: radius, length, evanescent_relative_error.
- `tone-hole.v1`: hole_radius, bore_radius, chimney_height, edge_radius.
- `circular-opening.v1`: radius.
- `lumped-reed.v1`: rest_opening, width, closing_pressure, attack, mass,
  stiffness, damping_ratio.
- `circular-string.v1`: radius.

Construction and every named ConfigurationState are checked. Assignments must
name actual parameters with matching dimensions. Unsupported imported geometry,
other part/model combinations, incomplete coverage, duplicate part/module
selections and mismatched material card/query bindings refuse. Required valve
and string elasticity/density resolve through MaterialSource and the existing
state-point owner; explicit Unknown refuses before card resolution. No property
is defaulted. Geometry-source citations and authored bounds are declarations,
not a proof of the bound or of real-specimen equivalence.

Domain policies explicitly supply frequency_hz for every module, ka for compact
radiation, edge_displacement_ratio for compact holes,
length_over_shock_distance and evanescent_relative_error for ducts,
mouth_pressure_pa for reeds, bow_force_n/pluck_force_n for strings. These
policies intersect owner limits. The §1.2 rounded audit cutoffs are
1.841*c/(2*pi*a) and 3.832*c/(2*pi*a), with no symmetric-only escape. Segment
geometry, Bessel/wide-tube scope, compact radiation, reed primitives and gas
coherence reuse their existing owner admission functions. A strict hole/bore
boundary is recorded with upper_inclusive=false.

Edge displacement uses the declared lowest excited frequency and conservative
pressure envelope. Inverse shock distance uses beta=(gamma+1)/2; linear string
drag requires KC=2*pi*displacement/diameter <= 1. A zero pressure gives zero
inverse shock distance without an infinity sentinel. The numerical policies for
“much less than” and evanescent error must come from the caller's declared
evidence; this validator does not invent them. These algebraic evaluations are
deterministic for identical inputs and perform no stochastic work or hot kernel.

G0/G1 tests `description_validate::dv_001` through `dv_010` cover analytic
cutoff agreement at 1e-12 relative, exact/next-value boundaries, ka=1 and declared
ka=0.5 policies, strict hole radius, unknown materials, slit refusal, mandatory
axes, retained limitations, dynamics-driven edge/shock/KC and configuration
coverage. Native execution is pending central verification. No audio golden or
wire schema changes; no specimen-fidelity, resolved anatomy/flow or model-form
equivalence claim follows from this token.

## Material-derived circuit storage

`MaterialCircuit::new_with_storage` binds every R/C/L branch to its material
owner, lowering dielectric/core coefficients through `fs-material::field_storage`
into the existing Kirchhoff descriptor. Original branch ordering, conserved
charge/flux linkage, source receipts and per-resistor Joule heat are retained.
Changing dielectric or core data changes the actual electrical evolution.

`rebind_storage_candidate` proposes a zero-time material/geometry change at the
exact same charge and flux linkage. It returns each branch's field-energy change
and the total signed parameter work. The coupled physical owner must supply or
receive that work when accepting the candidate; it is not counted as resistor
heat or source-port work. Incompatible ideal-source constraints, incomplete or
duplicate bindings, unrepresentable energy and cancellation refuse atomically.
C/L remain fixed within an ordinary circuit step. This linear lossless update
does not supply continuous thermal/mechanical feedback or switching impulses.

## Material conductivity at solved field and current

`MaterialCircuit::conductivity_candidate` and its thermally refreshed-seed
variant bind explicitly typed electric-field and current-density coordinates to
the original bulk conductor cards. Each bounded constitutive iteration calls
the existing Kirchhoff descriptor from the SAME accepted charge/flux history,
obtains actual bulk branch voltage/current, and refreshes the complete original
material bundle through `ResolvedMaterialStatePoint::at_query_point`. Geometry,
source basis, selected claims, requirements and other physical coordinates
remain explicit. Signed and magnitude axes are separate authored choices; a
source law that requires an electrical operating axis cannot silently freeze it.

Admission requires resistance consistency and an independent constitutive
power gate. The original native circuit step remains the sole owner of charge,
flux, source work and deposited Joule heat. The material evaluation at the
solved operating point and its resistance/voltage/power residuals are retained
separately; iteration discrepancy never becomes heat. Every reported iteration
is one native circuit solve. Unsupported source states, changed source/geometry
bindings, exhausted work and cancellation refuse without advancing history.
This is isotropic, uniform, quasi-static bulk conduction with fixed C/L and
geometry. It does not add switching, hysteresis, tensor/AC transport or contact
laws; contact-bound circuits require their separate explicit contact owner.
Five regression sources cover constant-law parity, analytic nonlinear RC
discharge, current-density/polarity response, complete temperature refresh and
atomic refusal. Native execution remains pending.

## Joint nonlinear conduction and finite enthalpy

`circuit_enthalpy::conductivity::step_candidate` places the same field/current
material iteration inside the existing finite-enthalpy residual. Each thermal
trial first refreshes the original source bundle at that proposed endpoint;
the original circuit then resolves the electrical operating point from the
same accepted charge/flux. Only its native per-branch Joule account reaches the
existing receiving cells, conduction and authored boundary/source owners.
There is no additional thermal integrator or duplicate heat calculation.

The shared driver counts every inner native circuit solve, including discarded
thermal stencils and fresh final admission. Source errors, operating-domain
failures, nonconvergence, exhausted work and cancellation publish neither
history. Authored conduction controls are admitted before the remaining outer
budget narrows them. The report retains both actual circuit work and material
evaluation at the joint temperature/field/current state. Thermal-only circuit
members preserve their original query schema without dummy field coordinates.

The existing first-order endpoint-thermal/midpoint-electrical discretization,
fixed geometry/C/L and lumped bulk heating limitations still apply. Scalar
response surfaces may depend simultaneously on temperature and an explicitly
owned field or current coordinate, using the canonical material evaluator and
only its source-declared interpolation cells. A two-resistor finite-capacitor fixture independently
eliminates the coupled field and temperature laws and two-cell enthalpy
equations; regressions cover data-only substitutions, energy balance, replay
and atomic refusal. Equation checks do not substitute for native test execution.

## Prepared material networks for nonlinear conductivity

`MaterialNetwork::prepare_conductivity` binds explicit electrical property
bases and typed field/current coordinates to the existing cell-resistor data.
`PreparedMaterialNetwork::resolve_conductivity_seeds` initializes the original
circuit from complete thermal/electrical material bundles.
`PreparedMaterialNetwork::conductivity_candidate` then advances the original
joint conduction/enthalpy owner without a caller-written material callback.
Cards, geometry, destinations, caloric instances and selected source claims
come from the same retained network; accepted continuation preserves those
pins while re-evaluating the actual thermal and electrical state.

The explicit operating axes apply to each resistor's complete source bundle.
After the original electrical owner resolves its actual E=V/L and J=I/A, the
thermal network reuses that same evaluated bundle to lower its half-cell series
conductances. Those coefficients enter the existing backward-Euler residual
before acceptance, including every nonlinear stencil and final admission.
Reported thermal samples and links describe the same coefficients used by the
residual. `step_candidate_with_transport` exposes this post-electrical update
through the original `StepCoefficients` type; the ordinary generic entry point
retains the identity update for independent thermal coefficients.

Ordinary thermal-only resolution stays strict. Passive thermal cells have no
electrical field owner and cannot admit k(E/J) by inventing coordinates. Missing,
duplicate or incorrectly typed electrical axis plans refuse. Native Joule heat,
external sources, work limits and atomic publication retain their existing
owners. Focused regressions cover k(E) and k(J) substitutions, independently
eliminated coupled equations over two accepted intervals, exact source-sample
equality, series conductance, native Joule accounting and solved-field source
refusal. These synthetic numerical cases do not establish measured material
accuracy or extend the fixed-geometry, scalar conduction model. A two-coordinate
regression additionally changes both resistivity and thermal conductivity with
the same actual temperature and field, retaining exact thermal/electrical
source-sample equality through two accepted intervals.

The same `PropertyValue::axes()` metadata admission applies to the existing
tribothermal, electrical-contact and plastic-contact material owners. Every
source coordinate must still have an explicitly owned, correctly dimensioned
physical state. A temperature/sliding-speed friction surface therefore feeds
the original coupled momentum/finite-enthalpy residual; it does not create a
second friction law or convert numerical loss into heat. The force-driven
surface regression independently checks the resulting velocity, heating and
momentum under data-only speed-slope changes.

### Nonlinear conduction through actual field motion and finite heat

`material_circuit::motion::conductivity` reuses the original conductivity
iteration with the original reciprocal field/body descriptor. Every trial starts
from the same charge/flux, gap, momentum and support history. Explicit source E/J
axes are refreshed from that trajectory's actual bulk voltage/current; native
Joule work stays authoritative and constitutive residuals never become heat.
Resistance and power gates both precede joint acceptance.

`circuit_enthalpy::motion::conductivity` deposits those native losses through the
existing finite-enthalpy driver. Its global work cap includes all inner native
descriptor solves, discarded thermal stencils and final admission.
`PreparedMaterialNetwork::field_motion_conductivity_candidate` uses existing
cards, cells, source pins, geometry and heat destinations directly. The same
solved E/J bundle also refreshes thermal conductivity and original link transport;
temperature, operating point, charge/flux and actual motion close together.

The original uniform-gap, constant field coefficient, guided inertia and optional
conservative-support scope remains. Conductor geometry is stationary; no moving
conductor transport, Hall law, hysteresis, field entropy or series-contact
response is added. The coupled thermal time level remains generally first order.
Focused synthetic regressions cover native constant-law parity, E/J-source
substitution, both electrical polarities, magnetic linkage, refinement against
an independent RC/Maxwell-force ODE, closed-form external work/heating, complete
prepared data, continuation, budget refusal and cancellation. Syntax parsing and
independent equation checks passed; native Rust execution remains unverified.

## Reciprocal material electric and magnetic gap motion

`material_circuit::motion::MaterialFieldMotion` composes the existing source-bound
RLC circuit, uniform dielectric or magnetic gap and regional rigid-body inertia.
The moving element's original quadratic term is replaced exactly once by
q²d/(2 epsilon A) or lambda²d/(2 mu N²A), with lambda the winding flux linkage.
All other C/L storage and original J/R/G topology are retained. The same descriptor
solve advances charge, flux linkage, gap and momentum. Terminal effort and
attractive force come from the same energy; body mass comes only from `ResolvedBody`.
Changing supported gap-medium or moving-body material data changes the coupled
motion and electrical history without changing a physics law.

`new_with_support` adds an optional layered axial return member to those SAME
gap/momentum coordinates. `motion::support::AxialSupport` delegates force,
elastic storage and per-layer strain admission to the original
`ConcentricThermoelasticSection`; reference length and gap map displacement to
axial strain. Canonical section density supplies consistent member inertia
sum(rho A) L/3 for the declared linear displacement shape, added to body mass.
Changing support layers changes return force, preload, motion and electrical
back-reaction. Elastic storage participates in the native descriptor, ordinary
enthalpy driver and thermoelectric-motion balance; it is never heat. Modulus and
density require validity-wide constant claims. Temperatures/eigenstrains remain
frozen; support damping, heating, live attachment, waves and buckling are outside
this conservative reduction. `material_field_support` covers both field types
against independent midpoint solutions, source substitutions through actual
charge/Joule/finite-enthalpy evolution, and atomic source/strain/cancellation
refusal. These regressions are added but native execution remains pending.

Original electrical ports keep their order and an external force port appends.
Bulk conductor and electrical contact losses reuse the existing native-effort
accounting and remain separately sourced heat destinations. Explicit refreshed
contacts must retain their IDs/branch assignments. Their losses are counted once.
The combined field/kinetic/elastic storage, original source work, force times actual
displacement, and heat must close before any joint candidate is returned.

The endpoint retains exact proposed charge/flux and rechecks the original circuit
constraints at the new gap. Incompatible ideal clamps, unsupported source data,
open gap bounds, unrepresentable values, balance failures and cancellation refuse
without partial publication or charge projection. Native domain-aware finite
differences shrink unsupported probes rather than expanding material validity.
A validity-wide constant field coefficient, fixed area and rigid translational
guide are required. Electric geometry assumes ideal electrodes; magnetic geometry
is one uniform gap with an ideal, negligible-reluctance return path. The resolved
permeability belongs to the gap medium, not an effective core material. This does
not model a length-changing solenoid, pole saturation/hysteresis, fringing, leakage,
impact, breakdown, dielectric loss or deformation. Thermal coefficient curves
require field entropy before reversible thermal coupling and are refused here.

Focused G1/G3 regressions cover data-only density/permittivity/permeability swaps,
both signs and zero of charge/flux linkage, actual winding count, driven RLC and
coil-contact loss, a 100 µm electric gap at 100 V, a 1 mm magnetic gap and atomic
constraint/domain refusal. Fixed-charge/fixed-linkage motion is checked against
the independently integrated constant Maxwell-pressure force. These are numerical
model checks, not experimental actuator or ferromagnetic-material validation.
Native execution is pending; source review and independent numerical oracles
are not a Rust test pass.

### Complete tetrahedral material valve resonance

`ModalAperture::from_material_solid` consumes the complete source-owned
`fs_solid::reduce::visco::MaterialModalValve` through the original full-C
`ExactZoh` pressure/contact/swept-volume solve. It returns the ordinary
`ModalAperture`, so existing two-sided and duct couplings use the same runtime.
No selected mode is removed, including modes orthogonal to mean slit opening.

The explicit coordinate map is `q = sqrt(total_mass_kg) x`, where q is the
source mass-normalized coordinate and x is metre-valued. Mass, stiffness,
complete viscosity, pressure, optional nodal-force patterns and nodal
observations all use this same congruence. No mode divides by its opening.
The source, full shapes, approximation report and gas remain attached; authored
diagonal damping is zero and the source viscosity is counted once.

`material_modal_valve` regressions cover pressure/force virtual work, complete
matrix scaling, source retention, material-driven excitation of a dark mode,
free-decay energy closure and time-step subdivision, and incomplete-state and
cancellation refusal. The reduction remains a frozen material image with one
effective mean slit; this does not add spatial solid contact, evolving material
memory or broadband accuracy beyond its original admission.

### Hereditary tetrahedral material aperture resonance

`ModalAperture::from_hereditary_material_solid` consumes the complete
equilibrium-mode `HereditaryModalValve` in the same existing aperture,
pressure/contact transaction and exact structural propagator. Original
tensor-Prony times and coupled factor rows become actual `ModalRelaxation`
states; the complete Kelvin matrix is present in that SAME exact block.
No source mode is divided by mean opening, and no frequency-frozen C or
independent damping ratio replaces hereditary loss. Pressure, optional nodal
forces, swept volume and `solid_nodal_motion` use the source's one
q=sqrt(total_mass)*x coordinate map.

Explicit relaxed/unrelaxed initialization selects the original Maxwell
viscous deformation, including its initial recoverable energy. The existing
exact owner independently integrates Kelvin and Maxwell heat along the
actual trajectory and returns their restoring reactions to the same momentum
and work equations. Failed trials publish neither motion nor memory.
The source remains available through `hereditary_solid()`; scalar plate
spectra are empty because the original tensor laws remain on that source.

All channels must fit the existing exact-state capacity. Every active source
time requires dt/tau<=2, and the complete instantaneous frequency must satisfy
the caller's angular-step ceiling in (0,pi). Original material strain limits
are checked at initial and candidate endpoints. These checks do not establish
within-sample strain bounds, arbitrary forcing-band validity, nonlinear
contact accuracy, modal/spatial convergence or thermally evolving material.
Heat reports total constitutive dissipation, without invented per-element deposition.

Regressions cover relaxed/unrelaxed storage, full-C plus tensor-memory
momentum and independent power integrals, original material stress/history,
dark-mode excitation, free decay, held-force sample subdivision, work-conjugate
nodal observation and transactional clock/strain/cancellation refusals.

### Semidefinite modal reed viscosity

`ModalAperture::with_viscosity` admits passive relative-motion losses with a
coupled undamped subspace, including a dashpot between two modal coordinates.
The existing shared pivoted PSD factor acts on mass-normalized C; its channels
must also reconstruct every physical entry within `1e-10` of the corresponding
diagonal geometric scale. The original C remains in the simultaneous force and
response equations. Midpoint loss uses nonnegative channel squares for singular
laws; exact propagation retains its integrated full-C loss. No damping floor
fills null directions. Unresolved reconstructions and indefinite laws refuse.
Previously accepted SPD and exact-zero laws retain their original arithmetic;
factorization stays outside the stepping loop. Authored raw matrices confer no
material identity. Regressions cover a lossless common mode, analytic relative
decay, simultaneous pressure/flow/work, signed coordinate scaling and refusals.

### Material modal reed relaxation

`ModalAperture::with_material_relaxation` binds generalized-Maxwell spectra from
the original regional plate materials to every retained reed mode. Equal
normalized spectra within each blade project exactly onto its original modal
basis, including oriented, tapered or heterogeneous equilibrium plates. Each
mode carries separate internal strains; the existing `ExactBlade` exponential
and independent loss quadratures supply storage, dissipation and restoring
force to the same pressure/contact/acoustic transaction. Rejected trials leave
material histories unchanged. Source substitution reaches returning ducts and
exterior PCM through the existing wind file and CLI, selected with
`loss proportional-maxwell relaxed|unrelaxed COUNT MODULUS TIME ...`.

The owner requires exact propagation, zero prestress, no separate material
viscosity, source coverage of equilibrium and instantaneous frequencies, the
original angular-step cap and `dt/tau <= 2`. Nonproportional regional spectra
require a full-matrix memory projection and refuse here. Source data remain
frozen during playback. Manufactured material regressions check stored energy,
independent force work, momentum, positive loss, trial rollback, data-only
substitution and PCM callback identity; they do not establish measured cane
behavior or complete clarinet/oboe fidelity.

`ModalAperture::with_regional_material_relaxation` removes the blade-wide
proportional-spectrum restriction while retaining each region's original
plane-stress constitutive approximation. Every region supplies its own modulus
ratios and times; its full rotated bending form is projected using the original
DKT displacement and rotation eigenvectors, retained by the same eigensolve
that supplied pressure, gap and force ports. Positive factors retain all
cross-mode reactions and explicit internal strains in the existing exact modal
owner. Equilibrium and geometric-prestress stiffness remain unrelaxed. Loss and
mean material force are independently integrated, and rejected pressure/contact
trials publish no memory. `loss regional-maxwell relaxed|unrelaxed COUNT
MODULUS TIME ...` reaches the existing wind file, CLI, metadata and PCM.

The complete exact loss block must satisfy `2*(3*modes+channels) <= 256` before
allocation. Original projection budgets, instantaneous frequency/source bands,
angular caps and `dt/tau <= 2` apply. Exact-zero strain rows carry no memory.
Localized regions can also have nonzero modal rows with a null strain subspace:
these use the same pivoted positive-semidefinite factor as localized interface
viscosity. Both reed and body memory retain the original cross terms within the
existing `1e-10` dimensionless reconstruction allowance, without added stiffness
or artificial damping in null directions. Successful positive-definite factors
retain their original arithmetic and channel order. Indefinite or unresolved
forms still refuse; every resolved channel counts against the original budget.
The old proportional selector retains its original diagonal propagation.
This is regional proportional-to-local-tensor bending memory, not arbitrary
three-dimensional Prony plane-stress condensation, nonlinear wet-cane physics
or measured instrument qualification. Physical-coordinate propagation/power
oracles and material-to-microphone/PCM regressions are authored; native source
compilation and execution remain required.

### Material two-plane nonlinear string composition

`MaterialStringSource::into_nonlinear` admits both transverse planes to one existing fs-nlmodal Kirchhoff–Carrier stress channel, with one source-derived EA, elongation, live tension and cross-plane stretching work. Original EI, prestress, mass, material/gas loss, torsion and directional centreline/surface/torque maps retain their ordering. `NonlinearMaterialString::prepare` retains that source while composing explicitly supplied receiver/contact/felt/volume participants in the original `PreparedImpactSystem`; `MaterialStringRun` implements the existing `ImpactSource` for the existing causal pressure renderer. The continuous vector-slope bound is invariant to transverse basis rotation. This is a frozen, moderate-slope reference model. Full rod/longitudinal dynamics, transfer from a nonlinear bow-friction bank, acoustic reaction, full instrument validation and exact-linear/SAV bandwidth work remain open for this impact route.

`SupportedMaterialStrings` also retains both source-selected transverse planes
on one physically directed translating support. Its exact kinetic transformation
counts each string's mass once; both planes share material EA and the affine
fixed-nut/moving-bridge strain. Bow and bridge directions may differ. Directed
centreline/contact rows, original torsion and each plane's stationary-gas loss
remain in the same mechanical solve. Independent spatial integrals check
kinetic energy, strain, drag and force maps; rotations and material/gas swaps
check the source binding. This is one bridge translation, not a six-DOF bridge,
resolved winding slip, dynamic longitudinal inertia or measured instrument
qualification.

### Source-bound axial loss at moving string terminations

`SupportedMaterialStrings::in_stationary_gas_with_longitudinal_waves_and_loss`
uses the same `LongitudinalMaterialLoss` and per-layer
`ConcentricAxialViscosity` as the reciprocal fixed-end bank. The original
supported dynamics retain both the relative axial resistance and the affine
endpoint-strain resistance. The existing exact kinetic pullback supplies their
support cross terms; no separate damping step or transverse/gas viscosity is
substituted. Equal endpoint translation has no affine axial strain loss.

Every original layer receipt remains available, including explicit zero loss.
The source band must cover DC through the complete coupled small-signal clock
bound. The linear uniaxial strain-rate law excludes finite-slope viscous terms,
heat evolution and axial fluid drag. Solid axial mass stays distinct from the
gas-loaded transverse mass. Existing states, support reactions and stationary
gas channels are retained in the original physical owner.

Focused G0/G3 regression sources cover the independent spatial strain-rate
power, both moving endpoints, rigid translation, explicit-zero parity, layered
material substitution reaching support/neighbor motion, and source refusal.
Rust syntax checks and an independent spatial quadrature check passed in the
implementation environment; native Rust execution remains required.

### Spatial closure on retained multimode material reeds

`ModalAperture::with_plate_closure` uses the original plate displacement
columns at every physical contact quadrature point and slit edge. Higher modes
with zero mean opening still change partial closure and contact. The same
triangle-area weights and transpose modal maps supply storage, force and work.
The existing contact owner's stable power-law secants are cached per point and
sample. Positive P1 gap integrals retain partially open edges; the shared
viscous-slit owner supplies local flow when selected.

The full modal Newton island uses the existing midpoint/full-viscosity or exact
material-memory response. Trials do not advance accepted structural or material
history. One or two blades share a matching planar mesh; their original
pressure and swept-volume projections remain authoritative. Moving pressure
footprints, nonmatching contact meshes and resolved moving-wall flow are absent.
Initial states still pass the original constructor's effective-slit admission
before the spatial closure is attached.

The optional `closure spatial` material-modal performance record reaches the
ordinary pressure renderer, music CLI and streamed PCM writer. Focused G0/G3
regression sources cover partial opening, zero-mean contact, local work, full
viscosity, regional material memory, source substitution, callback partitions,
CLI PCM parity and rejected-step rollback. Parse-only Rust and whitespace checks
passed; native execution and measured-instrument qualification remain open.

### Constitutive memory in the shared coupled modal owner

`CoupledModalSystem::with_material_relaxation` appends the existing fs-phs
Maxwell internal-strain coordinates to the original coupled exact islands.
Displacement projections connect those islands; mechanical q/v ordering and
all force ports remain unchanged. Each positive arm stores
`(sqrt(k) b.q - z)^2/2` and dissipates that mismatch squared divided by its
relaxation time. No artificial inertia or independent material time stepper is
introduced. Explicit initial history or fully relaxed initialization is required.

The owner independently integrates modal, connection and material loss on one
trajectory. Material restoring loads enter endpoint acceleration and pressure.
Scalar/multiple contact, scalar/multiple friction and nonlinear storage use the
complete account. Candidate q/v/z publish together only after acceptance;
rejected pressure, nonlinear or cancellation trials preserve accepted history.
Cold internal-component extension zero-extends original projections and keeps
memory; static equilibrium publishes the corresponding relaxed memory.

The original dense augmented exponential ceiling, instantaneous stiffness
Nyquist bound, total energy and explicit setup-work budget remain binding.
The existing prepared backend also advances material histories, in either
binding order. Its input-augmented dimension is `d = 3*N + M` for N mechanical
coordinates and M histories; the original Van Loan ceiling `2*d <= 6144`
and the explicit mode, byte and multiply-work limits remain binding. The same
scaled exponential and integral Gramians retain separate mechanical loss,
connection loss, material heat and external work. Instantaneous stiffness
balances the calculation while the mechanical baseline remains equilibrium
stiffness. Coefficient, constant-loss connection and radiation rebuilds retain
the histories. Empty memory retains the original dimensions and arithmetic.

Four focused native tests of this original prepared owner check complete
state/work/heat against an independent full fs-phs realization, whole-step
versus two-half-step agreement, 125 histories beyond the small backend ceiling,
and release of a stiff isolated history without invented mechanical work.
Additional contact/air rollback regression sources are present; those wider
integration tests were not executed in that focused native check. These tests
do not establish general continuum/modal or whole-instrument accuracy.

### Regional body-material memory through the original string performer

`ReciprocalStringBank::with_material_plate_relaxation` binds selected excess
modulus/time property pairs on every original plate-chart material region. The
chart constructor retains its complete DKT displacement/rotation eigenvectors
from the same solve as the bridge and pressure maps. `PlateBendingRelaxation`
projects each region's actual rotated equilibrium bending form into that basis;
all cross-mode terms and independently sourced regional histories survive. The
shared positive factor is also used by regional material reeds. A region's
excess modulus scales its existing tensor; prestress remains conservative. This
is fixed-state regional proportional bending memory, not a general tensor-Prony
plane-stress reduction or an assembled/carved-shell memory model.

The existing exact coupled owner advances each real internal strain alongside
strings, body, cavity, contact and friction. Restoring loads change both emitted
pressure and reciprocal bridge reactions on neighboring strings. Cold cylinder
gas and torsion rebuilds preserve the original body basis, source image and
accepted histories; internal axial-component extension preserves and extends
the same projections. Scalar/tensor body viscosity and authored body decay are
mutually exclusive with this selection. Source support covers DC through the
complete instantaneous coupled frequency bound, with positive-arm `dt/tau <= 2`.
Original exact-state, retained-mode and projection-work limits remain binding;
singular nonzero strain forms refuse without rank truncation.

The ordinary `material_strings --body-profile` v1/v2 plate route accepts
`loss regional-maxwell relaxed|unrelaxed MAX_STATES MAX_PROJECTION_TERMS COUNT
MODULUS TIME ...`. It resolves only the chosen material properties and attaches
the law after selected gas/torsion/axial construction, before the existing
multibow performer and PCM writer. An omitted record keeps the prior tensor
selection. Source bytes and identities remain inspectable. The mechanical-only
prepared backend refuses this selection before material-pack I/O, and augmented
setup work is reserved in the original budget. `regional-relaxing.body` is an
explicitly synthetic runnable example, not a measured instrument calibration.

Focused G0/G3 regression sources check regional data-only substitutions reaching
actual pressure, PCM and sympathetic-string motion; unchanged equilibrium
modes; gas/axial attachment-order equality; torsion/cavity composition; explicit
zero parity; callback partitions; source/budget refusals; and cancellation
rollback. Rust parsing and scoped whitespace checks passed. Native compilation
and execution remain required; no measured violin/cello/bass fidelity is claimed.

### Original string-material spectra in the shared instrument owner

An original homogeneous `ResolvedStringSpecimen` selected with
`with_prony_bending_loss`, or a bonded `MaterialStringSource::ConcentricProny`,
retains its exclusive constitutive law when lowered into `ReciprocalStringBank`. A private elastic-backbone lowering path
checks the original material identity, property selection and loss exclusivity;
public memoryless card and impact consumers still refuse the same law. Each
positive source arm contributes one internal strain per retained transverse
plane and sine mode, with stiffness `DeltaEI*(n*pi/L)^4/mu_loaded`. The actual
gas-loaded mass supplies `mu_loaded`; tension, axial stiffness and torsional
stiffness do not acquire an invented relaxation law.

Automatic source histories, regional body histories and the existing explicit
per-layer string selection share the same exact mechanical/material owner, in
that order. Explicit selection on another string preserves all prior history;
selecting a second bending-memory mechanism on an already exclusive Prony
string refuses, even when its selected strengths are zero. Cold gas changes
reproject automatic source histories with the new mass while preserving their
accepted energy-normalized strain. Explicit layered selections retain their
existing requirement to bind after gas loading. Torque reconstruction and
axial extension preserve all history groups. Every attachment validates source
support through the complete instantaneous coupled frequency bound.

The ordinary performer accepts `bending-loss prony COUNT MODULUS TIME ...`
after each string's versioned section/polarization rows. It requires original
equilibrium Young modulus and source excess-modulus/time pairs, without
requiring unselected bending viscosity. Each bonded layer resolves those roles on its own original material state.
The common `ConcentricPronySection` supplies every `DeltaE * I_layer` and
independent relaxation time, with no inferred bending viscosity or averaged
time constant. Its geometry and source receipts also feed torsion, both
transverse planes, support laws, cylinder gas loading and longitudinal dynamics.
The section's 128-arm bound and the complete owner's 256-history/setup ceilings
remain enforced. The existing `bending-memory` block remains the additional
per-layer path for homogeneous or bonded sections with selected Kelvin viscosity.
Automatic, explicit and regional body history allowances share one checked
setup budget before material I/O. Static preload initializes relaxed history;
material reactions reach the common body, neighboring strings and PCM writer.
`relaxing-strings.material-strings` and its source pack are synthetic examples.
Focused regression sources cover pressure/PCM substitution, coexistence,
gas/torsion/axial composition, zero parity and duplicate/band/budget refusals.
Syntax and whitespace checks passed; native compilation and execution remain
required. Frozen small-curvature bending does not establish axial/torsional
creep, interlayer slip, measured wound-string behavior or instrument fidelity.

`SupportedMaterialStrings::in_stationary_gas_with_longitudinal_waves` extends
that same support owner with material-solid axial inertia and selected
longitudinal sine modes. Actual 3D span endpoints, a declared transverse frame
and support direction determine every axial and transverse projection. The
exact triangular kinetic transform retains each directional affine-lift mass
once. `fs_nlmodal::LongitudinalString` supplies positive strain storage and its
derivatives, including actual moving-end extension; transposed forces drive
the same support and neighboring spans. Original transverse material/gas loss,
torsion and source states remain. The reference is already strung and
equilibrated, so initial static tension reactions are not reapplied. Axial
modes are conservative; axial viscosity and gas drag are not inferred.

This extension still has one translating support, without a full 6-DOF bridge,
moving torsional endpoints, winding slip or an exact linear propagation
upgrade. Existing `ImpactSystem` owns admission, stepping and work balance.
Independent spatial mass/strain/gas quadrature, gradient/Hessian,
support/neighbor work and material/gas/torsion regressions are authored. Syntax
and whitespace checks pass; native compilation and execution remain required.
No measured instrument realism is claimed.

`ReciprocalStringBank::with_longitudinal_waves` selects the existing
`fs_nlmodal::LongitudinalString` energy on the playable material bank. Original
homogeneous/bonded material supplies EA and SOLID axial mass; retained
transverse coordinates keep their own gas-loaded normalization. Explicit axial
sine modes append after every original string/body/cavity/torsion/fluid state.
Mapped modal storage couples them to both transverse planes in the same exact
linear/friction/DG transaction. The full positive strain energy replaces each
axial oscillator's linear stiffness once; the reported nonlinear correction can
therefore be signed. Contiguous original bindings keep borrowed coordinate
slices. Body states, connections, complete tensor losses, interface source and
pressure weights are preserved rather than rebuilding an independent voice.

`longitudinal-waves SLOPE STRAIN ITERATIONS MAX_FORCE ABS_TOL REL_TOL COUNT INDEX...`
selects this path in the ordinary material-string file through scheduling and
PCM. It replaces instantaneous geometric stretching. Requested indices, full
mode budgets, material elasticity/density bands, clock, rotation-invariant slope
bounds, strain, positive local tension and nonlinear closure remain explicit
admission conditions. Axial modes are conservatively elastic; original bending,
torsion, bridge, body and fluid losses remain active. Omitted axial channels are
statically equilibrated. No axial viscosity is inferred, and no independent
longitudinal end-force port, moving endpoint, helical slip or temperature evolution
is added. The shared nonlinear force preload below also equilibrates retained
axial coordinates. Focused regressions cover exact axial linear
limits, material frequency scaling, distinct fluid/solid mass, sympathetic motion,
total work, source-band refusal and callback-identical PCM. They do not qualify
measured violin/cello/bass fidelity or close the full string-unification program.

`FrictionModalSystem::initialize_released_equilibrium` and its simultaneous-bow
wrapper admit conservative nonlinear storage before sample zero. The existing
bilateral static response is factored once and preconditions a bounded energy
descent with safeguarded secant step lengths. Original storage energies and
exact static gradients include mapped replacement stiffness exactly once.
Every generalized force must close its existing absolute/relative tolerance;
the existing nonlinear iteration cap includes rejected backtracking trials.
No trial publishes state. Final component, connection-force and total-energy
checks stage the same original coordinates and relaxed Maxwell histories, then
publish them together without advancing a clock. Initialization storage is not
pre-window actuator work. `NonlinearPreloadDiagnostics` records the accepted
residual and storage; no stability or uniqueness certificate is implied.
The native material-string `forces preload` path therefore supports geometric
stretching, retained longitudinal waves and compatible conservative bow fibers
through their original coupled body/cavity/fluid network. Compliant normal
obstacles join the same stationary energy and exact gradient through their
original `ContactStorage`, including reciprocal receiver forces, original
accepted attachment maps and normal offsets. Dynamic contact loss is absent at
equilibrium. Each actual candidate normal force, penetration and material
kinematic envelope passes before publication; no source law is averaged or
replaced. Contact storage contributes to the complete preload energy while its
first dynamic report remains absent until the first accepted physical sample.
Retained shear films, obstacle residual compliance and nonlinear-only support
of otherwise free coordinates still refuse this static path; their dynamic
owners remain unchanged. Focused
regressions cover independent cubic and mapped-strain equilibria, strong
nonlinear loading, held/released motion, relaxed material histories, refusal
atomicity, source substitution and callback-identical native pressure/PCM.
The joint contact regressions independently solve a quartic-plus-normal force
equation and a two-body Hertz/strain equilibrium, preserve relaxed Maxwell
history, exercise held silence and physical release, and check a source-bound
receiver substitution in the native material performer. The manufactured
`nonlinear-stopped-preload.material-strings` uses the existing stop-profile
grammar. Source/format checks and independent numerical oracles do not replace
native Cargo execution or measured instrument qualification.

`ReciprocalStringBank::with_geometric_stretching` adds those same source-owned
quartic channels to the existing `FrictionModalSystem`, preserving every original
string/body/bridge/cavity coordinate, optional torsion, complete tensor damping,
and causal cylinder-fluid mass/memory. The actual loaded string mass normalizes
the geometric channel; material EA is neither fitted nor counted as linear
stiffness again. Each nonlinear trial reuses the original exact complete linear
response and original scalar or simultaneous friction branch/root solve, with
`fs_phs::discrete_gradient` supplying held constitutive reactions. Both reaction residual and total
storage/work balance must close before the original transaction publishes.
Rejected trials preserve accepted components, fluid state, clock and diagnostics.
The acceleration pressure path evaluates the endpoint geometric force. Existing
single- and multiple-bow scheduling retain this state, including the exact
one-port forwarding route. Several prescribed bow contacts share every DG trial
and one final publication; each contact retains its own traction, actual slip,
work and loss. The joint energy tolerance uses the sum of absolute bow supplies,
so opposing surface work cannot hide scale. Joint reports retain nonlinear
storage and reaction convergence. The material performance file selects it with
the optional
`geometric-stretching SLOPE ITERATIONS MAX_FORCE ABS_TOL REL_TOL` record, so
material and gas changes propagate through the bow, shared body and PCM writer.
Independent compliant stops and contact-derived bow loads use that same original
joint normal/tangential solve inside each DG trial. Initial and final contact
potential, original normal loss and signed support work `-R*delta_s` enter the
combined energy balance and energy budget; contact storage remains separately
reported from linear-plus-geometric network storage. The energy tolerance uses
the sum of absolute individual support supplies. Contact-loaded tangential laws
use the jointly solved current reaction and release on separation. Prescribed
support motion translates along fixed normals; optional source-derived station
motion below admits independent contact relocation. Rejected trials preserve accepted
contact reports and support offsets as well as the complete mechanical state.
Cold attachment is explicit; subsequent gas/torque rebuilds refuse rather than
dropping storage. Static force loading dispatches to the joint conservative
preload above when geometric storage is present; the linear-only contact
equilibrium owner cannot silently omit it.
This bounded iterated DG composition is not SAV, a full rod model, a unique
nonlinear-branch guarantee or measured instrument/trajectory qualification.

Independent compliant stops admit prescribed longitudinal migration along the
original full supported spans, composable with normal support translation.
Source-derived endpoint maps retain material/gas-loaded mass, declared
transverse directions and the complete string/body/cavity/torsion/fluid state.
The original joint normal/friction solver applies held reactions through mean
endpoint maps; unilateral potential and damping use the complete endpoint
closure change. Attachment work `R*delta_N^T*mean_q` and normal translation
work `-R*delta_s` enter the same energy acceptance with geometric stretching.
Acceleration pressure uses endpoint maps and endpoint geometric force.
Accepted stations, offsets, mechanics, reports and clock publish atomically;
constant station paths retain fixed-position arithmetic. Receiver footprints
and normals remain fixed. A stop supplying a fixed bow's normal load cannot
migrate independently; its bow-position path owns the joint normal/tangential
endpoint described below. Nonpositive
endpoint/mean-map compliance, trajectory limits and existing material, force,
penetration, energy and pressure limits refuse. The native v1 stop profile
carries bounded station trajectories on the original performance clock, with
a pressing/migrating/releasing finger example. G0/G1/G3 regressions cover the
independent work identity, endpoint pressure, refusal/retry/hold, complete
source-state retention and file/direct PCM identity; native execution remains
required. This is not a solved finger, evolving contact patch, rigid stopping
constraint or measured instrument qualification.

Prescribed bow-position migration uses the original `MultiBowStringBank` and
`MultiBowPerformance` transaction. The source string's surface map is evaluated
at each step's midpoint station, including its actual material/gas-loaded mass,
declared transverse directions and retained surface torque. This held map is a
power-conjugate velocity/force transformer: slip is `u*dt - B^T*delta_q`, with no
invented `Bdot*q` sliding term, geometric storage or frequency assignment. Every
traction shares the original simultaneous friction, normal-contact and optional
nonlinear solve. The original modes, material/fluid histories, body/cavity
reaction and pressure observer remain live. Independent stop station/normal
translation paths retain their own attachment and actuator work in that SAME
accepted sample. Geometry, physical state, reports and source endpoint cursors
publish together; a refused trial permits an unchanged replay. A later held
sample settles the tangent map onto the accepted endpoint, rather than leaving
the preceding interval's midpoint in place.

Finite ribbons translate all original quadrature cells together at their source
width; the complete footprint must stay inside the original supported string
span. Explicit speed, per-sample distance and complete response-rebuild budgets
bound admission. Receiver footprints and contact material stations along any
resolved fibers stay fixed, preserving their free length and mechanical state.
A contact-loaded normal follows its bow's string station through the original
endpoint displacement-map product rule, including attachment work
`R*delta_N^T*mean_q`. Its held reaction supplies that same port's tangential
friction and pressure checks. The material bank checks exact physical
co-location at each endpoint and retains normal direction, clearance, receiver
and constitutive law. Independent stop programs retain original placeholders
for these mapped contacts; the bow program owns their actual paths. Binding
normal contacts precedes migration admission, and native loading preserves
that order. A failed joint update restores both accepted maps and histories.
Retained shear/thermal layers still refuse migration until their finite patch
history transport is implemented; the generic joint entry point also refuses
moving their normal alone even if the tangent map happens to be unchanged.
Native material-string files expose bounded station knots before the original
gesture schedule, and use the same renderer and callback-independent clock.
Regression coverage includes independent midpoint sine maps and work/slip,
fixed-position parity, material/fiber/gas state retention, combined independent
stop work, refusal/replay, finite-width geometry and temporal refinement.
Joint-patch regressions additionally cover independent normal actuator work,
material/interface substitution, resolved material fibers, explicit stop-path
composition and native callback/PCM identity. This
bounded prescribed longitudinal migration does not resolve motion along bow
hair, changing supports, contact area, transported film, or measured instrument
accuracy; spatial and temporal refinement remain separate qualifications.

`ReedPerformance::from_bytes_with_gas_field` accepts the same source-retaining `AxialGasField` for authored geometry and gestures. Its exact inlet state drives reed density and impedance; the resolved midpoint cells supply one existing regional TMM for both feedback and terminal pressure. The exact outlet field state supplies the default exterior, and an explicit exterior record replaces that single binding for both flanged loading and baffled observation. The field remains frozen during rendering. Common pressure, full length, geometry/refinement and conflicting per-segment source declarations refuse through the existing owners. Standalone and ensemble sidecars can inspect the original primitive stations and source-pack identities.

The flue performance file accepts optional `walls COUNT FMIN_HZ FMAX_HZ` plus addressed `wall NODE isotropic PACK T_K OMEGA_RAD_S THICKNESS_M stress-free|strain-fixed elastic|kelvin-voigt` records after the original duct geometry, followed by `radiation prescribed|compact-baffled`. The CLI loads caller-relative NormalizedMaterialCardPack artifacts (8 MiB each, 64 MiB aggregate), while the library resolves original cards at independent solid T/omega states with SingleClaimOnly property selection and retained source/usage receipts. Shared artifact paths load once; geometry, branches, regional gases, existing losses and outlet addresses remain with their original owners. Each material wall occupies the existing degree-two node's physical dual strip; unselected and endpoint half-strips stay rigid, and actual area/coverage is reported. Compact radiation binds every declared observed outlet through FlueVoice::from_baffled_environment, then validates each installed binding against the lowered graph and the same exterior/receiver declaration. Material substitutions reach wall storage, jet-returning bore pressure and the common physical microphone; exterior substitutions affect both reaction and propagation when selected. This remains the existing frozen small-motion isotropic local-ring image plus independent low-ka baffled loads: no full shell coupling/exterior wall radiation, full-band unbaffled or mutual reaction, resolved turbulence or measured Boehm claim. The composition remains EstimateOnly and these additions do not close full M0.5 or M2.5-W3c/W4b qualification beads.

### Shared material moving keys on retained wind networks

`bernoulli_aperture::network::boundaries` attaches the original dynamic
aperture/contact/material owner to every explicit driven chimney endpoint.
The actual local gas, impedance, clock and integer or fractional transit must
match; closure retains propagating waves and material memory. Exterior reservoir
pressure and opening-positive mechanical force supply separately reported work.
All keys and the original graph, including any coupled radiation operator,
preview before any accepted physical state changes. This is the existing moving
slit approximation, with prescribed exterior pressure; it does not resolve a
circular pad or directly observe moving-slit radiation.

`ApertureNetwork` and `PairedApertureNetwork` bind those same keys before the
first sample and advance them through `step_with_boundaries`. Single or
independent double reeds retain their original nonlinear solvers. The nested
`ApertureNetworkFrame::boundaries` reports key supplies separately, while the
complete balance includes primary, key and duct storage and each supply exactly
once. Missing key drives refuse. Every participant shares a sample horizon;
extension validates all keys before changing any budget. The paired owner keeps
its original strict whole-system energy gate, including key residuals, and its
gated keyed step checks cancellation before publication. These APIs establish
shared physical coupling, not measured clarinet or oboe qualification.

`AperturePerformance` and `PairedAperturePerformance` accept the same
`performance::boundaries::BoundaryProgram`: complete ordered initial key drives
plus changes addressed by mechanical sample and original node. The existing
pressure renderer and baffled receiver consume the original coupled step;
opening-positive key forces and exterior pressure are data inputs. A keyed
source without its program refuses rendering before changing output. Cold
admission bounds pressure, source-force and key assignments together regardless
of attachment order. Preallocated control candidates publish only after the
physical step and receiver succeed, so callback size cannot shift key timing.
Key actuator work remains separate from primary reed work. This sample-addressed
library API adds no second interpolation owner or instrument-specific PCM path.

The flue file path resolves cantilever closure mass, stiffness and pressure
area through the existing material chart and DKT reduction, with elastic or
sourced Maxwell memory, contact and optional local-gas slit viscosity. Material
walls and closures share the bounded source cache. Typed pressure, jet-angle
and `ForceN` closing controls use the shared gesture sampler and reach the
existing pressure/PCM renderer. Geometry and material source bands still bound
the admitted approximation; substituting a material does not widen them.

### Fixed supply-passage momentum in material reed instruments

The shared `bernoulli_aperture::slit_channel` exposes its retained no-slip
profile as the exact next-step discrete affine boundary `dp = Z*Q - history`
for a fixed `StaticPressure` passage. Its source gas, dimensions and clock
supply the impedance; history is fluid momentum, not a second authored pressure
source. `TwoSidedAperture::with_supply_channel` adds this impedance and history
to the existing single, independent-pair or modal reed load. The same accepted
through-flow drives the channel and both acoustic ducts; physical blade-face
pressure subtracts the actual channel pressure drop. All participant trials
and the complete flow/energy account pass before channel, reed and paired
waveguide states are published.

The supply profile retains initial kinetic energy, release and signed reverse
flow. `TwoSidedFrame::supply` separately exposes storage, pressure work, wall
heat and backward-Euler numerical dissipation; both losses enter the numerical
account once. The material-wall and shared baffled-environment constructors
compose the same supply with original regional wall sources, Maxwell histories,
fractional/viscothermal propagation, radiation reaction and the physical
microphone. Existing schedules consume the resulting source directly. The
`two_sided_reed` example selects explicit passage geometry and mesh with
`--supply-slit LENGTH_M WIDTH_M HEIGHT_M TRANSVERSE_CELLS`.

The passage is additional fixed, planar, fully developed incompressible
Newtonian geometry, distinct from an existing duct segment or the moving reed
clearance. Closing-wall remapping, squeeze work, compressibility, axial flow
and evolving atmosphere remain outside it. A kinetic-head terminal refuses
affine condensation. Source/clock/load mismatches and complete-step failures
leave accepted state unchanged. Added regression sources cover pressure/flow
condensation, material single/paired reeds, initial fluid momentum, viscosity
interventions, physical microphone output and transactional refusals; native
execution and measured clarinet/oboe validation remain separate requirements.

`SteadyAxialGas::solve` in fs-material can feed its source-retaining `AxialGasField` directly to the existing regional duct, reed-bore and scheduled ReedPerformance field constructors. Wall cooling and explicit gas injections therefore change the same EOS/transport medium used for bore propagation, reed inlet loading and physical termination. The shared thermal owner supplies stationary enthalpy/species balances; the emitted field is frozen for acoustics. Its grid, caloric, conservation and low-Mach limits remain explicit and do not imply resolved transient breath evolution or measured clarinet/oboe equivalence.

### Explicit prepared propagation for complete modal networks

`PreparedModalBudget` selects cold preparation of the complete, fixed-period
linear mechanical operator. `CoupledModalSystem::new_prepared` and
`new_with_columns_prepared` admit initially larger networks; consuming
`with_prepared_linear_propagation` on the coupled, friction or reciprocal
string-bank owner preserves its original state and sources. The budget has
positive `max_modes`, `max_dense_bytes` and `max_multiply_terms` fields, with
a shared ceiling of 512 mechanical modes. Checked conservative additional
dense-storage and cold multiply-work counts precede prepared allocations.
Original setup-work, physical, material-band and transaction checks remain.

The prepared state retains every displacement, velocity and held force.
Original signed spring/damper columns, rest loads and material damping form
the same mass-normalized stiffness and damping matrices. A scaled full-state
exponential gives the held-force transition; separate component and individual
connection loss Gramians retain physical resistance ownership. Short-interval
Van Loan integration and Gramian doubling avoid a stiff inverse exponential
over the full sample period. Independent transition refinement, symmetry and
the complete quadratic energy identity gate cold acceptance; original sample
energy, per-port work and force/pressure ceilings gate publication. Preparation
does not truncate modes, replace the nonlinear friction solve or change clocks.

The bank retains the selected budget across cold torque and gas reconstruction,
and its complete owner continues through additional bows, stops, preload and
pressure/PCM rendering. All material-string performance versions accept
`prepared-linear MAX_MODES MAX_DENSE_BYTES MAX_MULTIPLY_TERMS` after the optional
cylinder record and before `cavity`. The file retains its 128-coordinate cap
and selects preparation before attaching torque/fluid. Its initial bank must
still fit the ordinary constructor; use the direct prepared constructors for
an initially larger generic network. Omission retains the ordinary path and
its existing 42-mode connected-island limit.

The supplied two-plane layered-string fixture contains 54 total mechanical
coordinates, including torsion, body, cavity and 48 free gas-memory inertias.
Its explicit 64 MiB / 200-billion-multiply-term budget authorizes costly cold
preparation, not a runtime or real-time guarantee. Regression sources cover
ordinary/prepared equivalence, noncommuting matrices and rest loads, individual
loss/work integrals, timestep partition, an analytic 43-mode connected system,
stiff free inertia, complete resource refusal and actual prepared material-string
preload/pressure. Native execution, numerical convergence of the chosen physical
basis and measured violin/cello/bass fidelity remain separate requirements.

`ReedPerformance` accepts an optional authored `gas_field` block after segments, selecting the existing continuous Magnus or midpoint field owner with explicit resolution, station positions, temperatures, canonical source IDs and mole fractions. Numeric resolution remains midpoint shorthand. Only caller-loaded species/NASA-9/transport packs supply the material model; declared components may be exactly zero at the inlet or throughout the field. A positive dyadic composition is used solely to admit the complete source catalogue and is never evaluated physically. Original knots and pack identities remain retained, with active knot receipts on the existing gas-record/provenance surface. The same field supplies acoustic propagation, exact inlet reed coupling and default exterior loading/pressure observation. Existing exact-length, common-pressure, source-domain, phase, geometry and aggregate-work refusals remain, and explicit segment or injected-field coexistence is refused. Wind and ensemble sidecars distinguish continuous steps per smooth span from midpoint cells and authored segment samples. Regression coverage includes source/order/zero/seed invariance through actual midpoint and continuous pressure, callback parity and the unchanged absent-block path. This remains a prescribed frozen field through the existing compact reed/duct/radiation model, not a transport/reaction solution or measured instrument qualification. The checked-in N2/CO2 example includes the existing source-pack construction and rendering commands.

### Material fluid films and finite heat

`render::plate::impact::squeeze::material::MaterialSqueezeFilm` resolves the same `FluidTransport` snapshot used by flow and heat, retains the concrete source and its receipts, and compiles viscosity into the original `fs-tribo::ResistiveFilm`. Its pressure and analytic tangent remain inside the existing joint mechanical equation. `HeatedImpactSystem::new_with_films` routes selected contact and incompressible-film work to the same finite caloric stores. Film heat is evaluated by the original pressure owner at the actual discrete-gradient velocity and midpoint geometry, then deposited once before mechanical acceptance. Existing thermal-network subcycling, conduction, external sources and rollback apply unchanged. Network reports distinguish `deposited_heat_j` from true combined-system `external_heat_j`.

The film is frozen, incompressible, Newtonian, no-slip and quasistatic. Caller-declared thin-gap, pressure and applicable inertia/cavitation/rarefaction restrictions remain; no continuous viscosity-temperature feedback, added fluid inertia, phase-changing geometry or implicit wall/fluid heat partition is inferred. Compressible-film free-energy loss and unrelated mechanical/acoustic losses are not relabelled finite heat. The source binding remains inspectable while mechanics owns its immutable compiled pressure operator. Regression sources cover data-only liquid viscosity/caloric swaps, analytic pressure/work, unchanged mechanical stepping, finite temperature and network conduction, and atomic refusals. Native execution was unavailable in this session; the independent numerical reduction is not a Rust test pass.

### Material relaxation through finite heat and transport

`HeatedImpactSystem::new_with_losses` also binds explicit `RelaxationHeatRoute`
destinations for the original Maxwell arms, including arms projected from shell
bending spectra. Each selected arm deposits `dt * e_z * (e_z / tau)` from the
successful `StepWorkspace` discrete effort. That effort retains the complete
nonlinear Gonzalez correction and is unavailable after a failed numerical
attempt. Endpoint power, total energy deficits, and stored branch energy are
not substituted for the arm's irreversible work.

Branches are routed once in canonical order into existing `CaloricReservoir`
stores; unselected losses remain outside that thermal boundary. The same route
feeds uniform transport and the original subcycled enthalpy network. Material
history, motion, heat, and accepted time publish together, with no second
mechanical loss debit. Caloric-domain or solver refusal retains all accepted
states. Constitutive coefficients stay frozen over this coupling: it does not
infer temperature-dependent modulus, relaxation-time feedback, spatial heat
allocation, or measured material validity. `material_relaxation_heat` exercises
an independent nonlinear polynomial-work oracle, selective routing, retained
history, physical conduction, common clocks, and transactional refusal.

### Source-driven thermoelectric transport

`circuit_enthalpy::thermoelectric::step_candidate` composes the existing prepared
material network, Kirchhoff descriptor and enthalpy residual. Each absolute-Seebeck
solid leg binds an original voltage-source branch in exclusive oriented series
with its material resistor; bulk distances and areas must match its two
source-derived conductive links. One pinned scalar or piecewise-linear S(T)
claim covers both junctions and the actual Joule-heated bulk. Source integrals
on a-to-bulk and bulk-to-b thermal edges supply their mean coefficients. The
native midpoint current drives heat `[-I*S_left*T_a,
I*(S_left-S_right)*T_bulk, I*S_right*T_b]` into those same three caloric stores.
The existing resistor owner deposits Joule heat once in the bulk. Each thermal
edge closes electrical/thermal work and reversible entropy transfer; a constant
S retains zero bulk transfer. For variable S, the bulk receives the interior-node
Thomson allocation and the terminal edge averages contain the remaining spatial
allocation. These discrete nodal transfers require spatial refinement; they are
not exact continuum pointwise Peltier fluxes or a whole-leg Thomson source.
Shared junctions sum the signed material contributions. Internal Seebeck work
and ALL reversible heat are removed from external-supply accounting before
the total electrical/thermal balance gate. All source, geometry, numerical and
cancellation refusals preserve both histories. The zero-field, fixed-geometry,
single-phase leg still excludes tensor/AC and electrical-contact extensions;
other thermal cells retain their original phase and boundary behavior.
Endpoint coefficient coupling and lumped Joule deposition require time/spatial
refinement. Focused G1/G3 cases cover integrated affine and piecewise S(T),
current/gradient reversal, signed cancellation, actual bulk enthalpy feedback
and out-of-support atomic refusal. These numerical fixtures are not measured
material validation. Physical convention reference:
https://doc.comsol.com/6.4/doc/com.comsol.help.heat/heat_ug_theory.07.087.html .

### Nozzle flow through finite electrical and thermal networks

`cavity_phs::piston::exchange::network` stages one original pressure-driven
species/enthalpy transfer, then invokes the existing circuit/enthalpy gas-wall
network solve from that same transferred composition on every nonlinear trial.
`nozzle_material_network_candidate` retains prepared electrical/conductive source
claims and the actual gas midpoint; the material-bound piston entry point also
retains its original regional body and inertia. One candidate returns thermal,
electrical, piston and finite-reservoir endpoints together. Full-endpoint mass
and energy checks enclose both gases, electrical storage, cell enthalpy and
kinetic energy, with external work/heat separated from internal nozzle enthalpy,
Joule loss and reciprocal wall heat. A shut path preserves the closed-network
state and port history. Nonzero flow retains the explicit missing-stream audit
boundary of the chamber's two-port observation; the outer paired balance is
authoritative for this finite system.

This is a first-order transfer-then-network split at one physical duration,
with a uniform fixed-volume reservoir, original ideal-gas/species sources,
fixed thermal geometry and rigid piston. Existing donor, pressure-crossing,
source-domain, clock, event, work and cancellation checks remain. It does not
supply combustion, jet thrust, moving-mesh conduction or a simultaneous nozzle
residual. Regression sources cover material/mixture interventions, both flow
directions, one withdrawal during repeated trials, full storage balance,
closed-path parity and complete rollback; native execution remains unverified.

### Molecular diffusion through the complete material assembly

`MaterialExchangeState::diffusion_step_candidate` and `advance_with_diffusion`
connect the original source-backed Maxwell–Stefan network to the same reacting
gas, finite thermal cells, electrical circuit and moving pistons. Every face
reads the original accepted composition/T/p, resolves explicit binary-pair
diffusivity and transfers species plus their partial enthalpies once. All
outgoing inventories are admitted jointly before incoming species can replenish
them. Each chamber retains its existing caloric solver. Changed composition
therefore enters actual chemistry, wall heat, charge/flux and momentum feedback.

The whole-interval mass/energy balance starts before diffusion. A gross inlet or
outlet makes the chamber's closed two-port audit unavailable even if net mass
is zero; the enclosing account retains actual storage and external work/heat.
All owners publish atomically after the original source, entropy, withdrawal,
nonlinear and conservation checks. Chemical formation energy and carried
enthalpy are never added again as heat.

This is a first-order diffusion/reaction/coupled-heat-motion split over one
interval, retaining the original ideal-gas, molar-average, approximately
isothermal/isobaric face domain and explicit temperature/pressure limits.
It supplies no thermal/pressure diffusion, membrane forces, resolved bulk flow
or simultaneous nozzle/diffusion solve. Arrivals may react now but cannot cross
a second diffusion face until the next interval. Focused regressions cover
the binary limit, diffusivity-only physical response, post-mixing reaction,
zero-net-mass stream accounting, closed faces and transactional refusal.
Native execution and experimental qualification remain separate requirements.

### Combined species transport with selected material heat laws

`MaterialExchangeState::{transport_step_with_gas_law_candidate,
advance_with_transport_and_gas_law}` lets the original combined nozzle and
Maxwell–Stefan transport use the same typed gas-heat callback as the closed and
sequential-flow assemblies. Transport samples one accepted snapshot and checks
one gross species budget, then chemistry runs once and the existing joint
electrical/enthalpy/heat/momentum solve receives those actual material gases.
The heat law sees each trial gas midpoint and exact cell material endpoint;
its source payload and refusal type survive unchanged. Empty path lists select
either transport mechanism alone, including a closed assembly.

All transport variants and heat laws now share the original whole-interval
mass/energy audit and changed-inventory port accounting. There is one accepted
clock and one publication; carried partial enthalpy or reaction formation
energy is not deposited again as heat. Failed transport, heat-law resolution,
nonlinear solve or cancellation leaves every owner intact. Existing
first-order splitting, distinct physical transport paths, molar-average
diffusion domains and material/geometry restrictions remain; this is not a
resolved advection-diffusion or boundary-layer solve. The existing material
exchange regressions exercise conductive-law parity, actual post-transport
chemistry once, typed late heat errors, atomic publication and cancellation.

### Curved material shell bodies in the reciprocal string network

`shear_shell::compile_material_shear_shell` binds original
`ShellRegionMaterial` selections to one admitted 3D midsurface. Each selection
retains its material state, ordered plane-stress model, thickness constraint,
global grain axis and shear correction. The grain axis is projected onto each
physical facet before the additional in-plane angle is applied. Orthotropic
transverse shear uses the source's `shear_modulus_13` and
`shear_modulus_23`; isotropic shear comes from the same source E/nu.
Missing transverse elasticity or any of the six in-plane and three transverse
viscosity entries refuses. Signed cross terms remain explicit and validated.

The material image uses the existing MITC3+ stiffness, support transforms,
mass choice, eigenvalue owner and physical mode expansion. Viscous membrane,
bending and shear energies use the same elastic bubble condensation; numerical
drilling stabilization receives no invented material loss. Full modal damping
is factored by the existing passive dashpot owner and retained in the common
mechanical solve. Original six-DOF force/moment rows and actual facet-normal
volume integrals are projected on those same mass-normalized modes. Independent
assembly, projection, modal and complete-network budgets remain explicit.

A source that also states `constant_loss_reference_omega` [rad/s] declares that
its viscosities are its loss factors at that angular frequency (viscosity =
loss factor x modulus / omega_ref). Its elements are projected on their own and
read at every retained mode's own frequency: the projected matrix C becomes
S C S with s_a = sqrt(omega_ref G(omega_a) / omega_a), where G is a fixed
five-branch Debye (generalized Maxwell) spectrum fitted to unit loss over
50 Hz-10 kHz (`shear_shell::constant_loss`; flat to 0.990-1.013 over
100 Hz-2 kHz). Each mode's damping then carries the source's loss factors at its
own frequency instead of growing with it. The anchor's claim domain and that
band clip the region's frequency support. Sources without the anchor keep the
single unchanged Kelvin-Voigt projection, bit for bit. This is a per-mode loss
law: the stiffness dispersion that a causal constant-loss solid implies is not
applied.

`ResolvedShearShell::with_coatings` adds ordered, source-bound films on either
side of each substrate region without moving its midsurface or attachment
skins. The original laminate/MITC owners retain full A/B/D, I0/I1/I2 and each
film's offset viscous coupling; the regional shear correction remains explicit.
By default, storage is frozen at an authored reference, with every retained
mode checked against both source support and an authored relative
storage-drift bound.
Each film is projected separately and reads its own E''(omega)/omega; positive
diagonal congruence gives exact modal-diagonal loss for that basis and a passive
geometric-mean cross-mode approximation. Wood's anchored law never rescales
film loss. Complete projection budgets count every material-law group.
Replacing or removing the stack and rebinding wood preserve the original
substrate mass/thickness policy, interfaces and pressure ownership. Films add
mass. Source, geometry, budget and cancellation refusals are transactional.
Perfect bonding and linear small strains are assumptions; this frozen option
supplies no hereditary memory. Neither coating option supplies penetration,
cure stress, slip, an evolving local strain-amplitude gate or a displaced
acoustic boundary.

An explicit `ResolvedShearShell::with_coating_relaxation` instead rebuilds each
film at its source equilibrium coefficient E0, retaining the same signed
offsets, A/B/D, I0/I1/I2, bridge maps and pressure geometry. The original MITC
viscosity projection yields each film's complete unit-modulus PSD form H in
the retained equilibrium basis, including membrane, bending, shear and their
cross terms. Positive scalar Prony arms multiply H; equal relaxation times
are summed across all regions/components before the shared PSD factorization.
This preserves the transfer while bounding history count by retained rank
per unique time, rather than by number of coated faces. This aggregation
covers relaxed or uniformly unrelaxed initialization; it does not identify
independently authored histories for every source film.

`ShellCoatingRelaxation` retains projected films, original source identities,
the intersected support band, separate sampled storage/loss fit certificates,
and complete coupling rows/time constants. Wood/glue and interface viscosity
remain in their original owners; no frozen film damping is also added.
Projection work, sampled approximation tolerances, history count and
cancellation are explicit transactional gates. Coated parts in an assembly
must select one consistent memory policy. The finite-band
`prony-continuation` selection explicitly accepts the globally passive fitted
continuation for a transient; it does not widen the source's measured domain.
The fixed modal/bubble transformation is a Galerkin approximation, not
frequency-dependent local condensation or a convergence certificate.

`ReciprocalStringBank::with_material_shell_coating_relaxation` explicitly
binds this retained image to the original coupled memory owner. Library callers
must bind before playback; shell constructors alone return the equilibrium
image. The material_strings frontend performs this binding automatically when
`coating-memory prony-continuation TERMS STORAGE_TOL LOSS_TOL MAX_HISTORY`
is selected in a V3/V4 body profile or assembly options. The same body rows
enter the shared string/body/cavity state and survive source rebuilds and
additional string-memory selection. Prepared rendering uses the existing
resource declaration. The Cannone O1 recipe's 16 candidates retain 13 positive
times; 30 body modes need at most 390 histories. Its conservative admission
bound is 3,954,407,664 bytes, within the existing 4 GiB declaration. The
independent 512-history cap does not guarantee that a 512-history case fits
that byte limit.

Both shell image types expose `coating_harmonic_response`. For bare images
or equilibrium/hereditary coating images it solves the unloaded structural
response with the original complex LU. With exp(-i*omega*t), the dynamic
matrix is K0 - omega^2*I - i*omega*C_background plus the sum of
[(E'(omega)-E0)-i*E''(omega)]*H; interface viscosity is added independently.
Source support is checked at the common drive frequency. Frozen-film images
refuse this exact-source analysis. Returned bridge and signed exterior/cavity
volume velocities, supplied power, separate wood/glue, interface and per-film
heat, and original-equation solve/power residuals permit physical comparison
without adding a sound-pressure gain. The frontend
`--body-response NEW.tsv MIN_HZ MAX_HZ POINTS BRIDGE_INDEX` exports 2..1024
log-spaced frequencies within a complete work limit and 8 MiB report limit,
validating all samples before a create-new write. Unit-force normalization
does not certify finite-amplitude strain. The structural FRF has no string
or air-pressure feedback; the ordinary coupled WAV route supplies those
separate dynamics. Focused shell, harmonic and frontend native checks do not
constitute a completed coated-violin render or mesh/mode/time convergence.

The additional `--body-response-only` selection omits the WAV positional
argument and projects the same body through the original shell/assembly
resolver without constructing strings, cavity dynamics, the bow performer or
the prepared time propagator. The performance still supplies body dimensions,
material state, mode count and input limits; only body/interface packs are
loaded. Playback-only options refuse. This remains the unloaded exact-source
structural transfer and opens no WAV output.

`ReciprocalStringBank::material_shell_bridge_harmonic_response` and the same
method on `MultiBowStringBank` instead force one original physical shell bridge
inside the complete connected performer. They reuse the original coupled
harmonic owner with every retained string, cavity, bilateral connection,
Maxwell history and selected exterior radiation pole. Contacts are explicitly
released; this is not a linearization of nonlinear bowed playing. The pressure
readout is exactly the original time-owner observation: acceleration rows use
-omega^2*q and real modal transfer rows use Re(T)*v + Im(T)*omega_mode*q.
No extra receiver transfer, normalization or pressure gain is introduced into
this native pressure. An explicitly selected exterior receiver is separate.

`material_shell_bridge_harmonic_matrix` on both banks returns one unit-peak-
force column per physical bridge and all observed bridge velocities, sharing
the original harmonic factor. `material_shell_bridge_forces_harmonic_response`
combines one coherent complex force per bridge before solving and accounting;
pressure superposes, while power, branch heat and stored energy retain cross
terms and cannot be obtained by adding isolated column losses.

The harmonic owner retains branch ordering in per-Maxwell mean storage and
heat arrays, respectively |c^T*q-z|^2/4 and |c^T*q-z|^2/(2*tau). The bank
attributes coating heat and storage through its original bound history range;
the complement is other material heat. Equal-time factors can mix films, so
this loaded query does not invent per-film allocation. It evaluates the actual
positive Prony playback realization, whereas `coating_harmonic_response`
evaluates the original film law. Source/fitted-operator frequency supports
are intersected across the retained strings, shell, films, interfaces and air
fit, with bound fluid impedance admission at the drive frequency. Bare images
are supported. Frozen films, unbound projected coating memory, unsupported
nonlinear/adhesive/contact/local-strain attachments and residual flexibility
without complete direct bridge-force columns refuse. Queries change no accepted
state, history or sample clock. Source admission is not an amplitude or
mesh/modal convergence certificate.

`ShellHarmonicReceiver` borrows the actual multi-bow bank and owns the original
receiver geometry and source traces. Its prepared queries perform the bank's
harmonic solve internally; foreign same-sized response arrays cannot be bound.
A selected `shell_observer` uses the original reciprocal P1/P0 observation owner
with the same complete closed-shell load, gas, basis and surface/mouth map.
A selected `one_way_observer` uses the profile's original ruled Burton-Miller
skin producer, every retained body mode and the actual aperture-flow shares.
Playback and harmonic observation share source-coordinate conversion, including
spatial-mouth inverse normalization and split neck-flow shares. The producer is
evaluated at the forcing frequency with raw complex pressure and gradient;
no FIR latency, taper, normalization, second mechanical reaction or extra loss
debit enters. One BEM frequency sample serves all bridge matrix columns. The
report records receiver/gas/formulation/operator/loading/skin identities,
source counts, pass band, observation work and loaded versus one-way body-mode
scope. Omitted feedback is not certified negligible. Boundary, quadrature and
projection work are reserved before mechanics; skin work counts every actual
256-source batch. Three focused production-helper native tests pass for phase,
source mapping, invalid inputs and combined/batched work; they do not execute
a BEM solve. Full bank/microphone integration tests remain authored and unrun.

The frontend `--loaded-body-response` exports bridge mobilities, signed shell
volume velocities, native and selected receiver pressure in Pa/N, supplied and
component/connection/coating/other-material/radiation power, stored energy and
solve/power residuals before the ordinary WAV render. The v2 TSV identifies
`drive_bridge`; a final bridge argument of `all` exports the complete matrix in
frequency-then-bridge order. A selected exterior receiver adds three complex
pressure gradients; without one, receiver pressure equals native pressure.
The existing prepared multiply-work declaration bounds the entire sweep,
including all columns and receiver work, divided across frequencies; small
nonprepared runs retain their allowance. All samples pass before a create-new
write, with the same 2..1024 point and 8 MiB limits. Receiver sample allowances
and pass bands remain binding alongside all source supports. This selection
still builds the complete performer/prepared propagator. Plate microphones
remain unsupported by this shell-specific report. Focused native containers
do not validate complete bank construction. Normal-Cargo bank and projection
integration regressions remain unexecuted.

`--loaded-body-response-only` takes the same report arguments without a WAV
positional path. The frontend shares physical bank construction with playback,
selects `HarmonicOnly` from the required explicit `prepared-linear` resource
record, then writes the loaded response before gesture scheduling. It creates
no coupled time-transition or heat-integral matrices and performs no held-step
bow-compliance solve. The multiply-work allowance bounds the complete sweep.
The parser retains computed BEM/Rayleigh source declarations; actual prepared
playback still refuses those backends before material I/O, while harmonic-only
construction can retain admitted shell BEM or skin loading. The library's
Rayleigh binding likewise distinguishes preparation purpose; this shell-only
CLI does not admit a planar body. Source identity, geometry, gas, fit-band and
resource gates remain unchanged. Playback-only options and plate microphone
selections refuse. Existing shell/skin microphone declarations instead select
their original continuous frequency producers without building playback FIRs.
Three production option-parser tests pass natively, including `all` selection;
three exact-production v2 writer tests retain original coupled numerics behind
small bank/receiver transport facades and verify native/selected pressure,
gradient units, matrix row order and create-new/refusal behavior. These do not
validate a BEM solve or complete CLI build. Closed-shell BEM and skin bank/receiver
regressions remain authored and unrun.

`SharedMaterialShellBody` and
`ReciprocalStringBank::from_material_shear_shell` connect that image to every
original material string, its own physical bridge map, shared gas and optional
body cavity. Shell source and image accessors retain geometry, region packs,
physical modes, force maps, volume shapes, eigenvalue intervals and damping.
Existing torsion joints, two-plane strings, gas-memory inertias, additional bows,
stops, preload, scheduled performance and pressure continue through the same
reciprocal owner.

The model is an unloaded, frozen, linear midsurface with homogeneous supports
and explicit numerical drilling stabilization. Its compact baffled observer
does not imply full exterior radiation feedback or an assembled, prestressed
violin/cello/bass corpus. Mesh/mode/stabilization/time refinement and measured
instrument qualification remain separate from source-level composition tests.

`shear_shell::pressure::ShellPressureSurfaces` independently selects exterior
and outward-from-cavity original facets. Both rows use the shared
`fs-plate::shell::pressure` owner and the same material eigenbasis. Exterior
acceleration pressure uses only the exterior row; the existing single
`HelmholtzModalCavity` compression spring uses only the cavity row. Internal
structural members retain their complete mechanics without being counted as
fluid boundaries. Explicit empty rows mean zero participation; omitted source
selection preserves the historical whole-midsurface path. Material rebinding
retains authored fluid roles and geometry. The uniform cavity, midsurface
pressure and compact observation remain model-form rungs: this does not imply
closed finite-thickness fluid skins, resolved interior acoustics, or exterior
radiation back-loading. Source tests exercise independent physical pressure and
cavity rows, identical mechanical histories under observer-only changes, one
collective air spring, material rebinding and zero exterior participation.

### Prepared material conduction on the shared thermal timeline

`PreparedMaterialNetwork::observe` composes the original pinned transport
resolver with another observer on the existing enthalpy substeps. Networks
with resistor branches refuse; circuit storage and Joule heating retain their
circuit owner. The shared observer pair owns heat-source addition,
competing-link refusal, timestep composition and atomic history acceptance.
Material snapshots are refreshed at actual thermal endpoints and events,
while authored link coefficients remain unchanged. Conduction remains the
declared uniform-cell, fixed-geometry reduction; temperature-dependent source
coefficients do not by themselves establish time or spatial convergence.

### Material interface compliance in the shared string shell

`shear_shell::connections` retains complete ordered interface states and claim
identities, exact surface material roles, property names, actual contact area
and weighted footprints with explicit global offsets. Source stiffness per area
[N/m³] and viscous resistance per area [N s/m³] require finite nonnegative
fixed-state laws with intersecting frequency applicability. Every retained mode
must lie in that band.

Regional membership, positive normalized weights, physical moment arms and
coincident application points for two moving sides are checked. Interface
stiffness enters the monolithic shell pencil before the shared eigensolve.
Interface viscosity remains its original passive connection channel, supporting
semidefinite interface-only loss alongside the separately owned full bulk
material viscosity. Elastic stiffness is never duplicated as a runtime spring.

The reciprocal string bank budgets these additional channels and retains their
sources through existing extensions. Bulk material rebinds recheck the ordered
interface surface roles. Fixed counterparts have zero motion and no modeled
mass. Regressions cover changed material support stiffness and viscosity reaching
actual microphone pressure, exact projected moment arms, source rebinding,
interface-only loss and source/geometry/frequency/resource refusals.

This finite bilateral model does not establish continuum mortar, glue coverage,
opening contact, complete assembled instrument anatomy, thermal/history
evolution, a prestressed equilibrium or measured timbre.

### Flue jet-direction controls and gas-bound mouth separation

`FlueControl::SetJetAngle` changes only the prescribed direction in the
existing fixed-window tilted-cut image. Its effective offset is the original
`labium_offset + cut_up * tan(angle)`; the resulting injected-flow change
enters the same source-pressure derivative and mouth-work balance. Finite
direction and transverse-coordinate admission precede held-input mutation.
Acoustic, convection, windway, material-wall and radiation histories persist.

`FlueOutletPerformance` accepts one blowing-pressure track and an optional
typed jet-angle track. The original pressure compiler and `GestureSchedule`
sampler share exact ceiling-rounded mechanical timestamps; both controls apply
before that sample. Whole-phrase admission bounds every track/tick assignment
and source-sampling work before stepping. The existing multirate observer
retains partial frames across controls and callbacks.

The file's optional `mouth_loss contraction CC MACH SOURCE` row binds the
existing `ContractionSpec` to the actual mouth window and actual inlet gas.
Its separation loss enters the same implicit acoustic-mouth solve, remains
distinct from radiation/viscous loss, and checks the solved contracted Mach.
Omission and `mouth_loss none` select the original linear mouth.

Focused G0/G3 regressions compare static geometry against angle controls,
manually sampled pressure/direction against the actual shared microphone,
retained moving-gas histories, callback partition invariance, data-only gas
changes and file/CLI construction. These tests do not establish measured
instrument fidelity. The direction remains an instantaneous geometric
approximation, not resolved transient jet turning, moving lips, coverage,
time-varying mouth area/cut-up, or a material-state evolution law. The
separation closure does not resolve vortex memory or turbulence.

### Retained material inertia in contact and electrical feedback

`solid_enthalpy::plastic::dynamics::ThermoplasticDynamics` now composes with the
original pressure-dependent thermal supports, finite circuits and electrical
support contacts. Its `advance_with_thermal_contacts` and series-conduction
variant use the same contact Newton residual as the quasistatic solid. Its
`advance_electrical_candidate`, support variant and pressure-contact variant use
the same original circuit/enthalpy residual. Existing quasistatic entry points
retain their signatures and omit inertia.

Every thermal trial calls the original native small-strain mechanical owner
with the same accepted displacement, velocity, J2 history and source-bound
reference mass. Endpoint stress supplies support compression; that pressure and
the actual two-sided temperatures resolve the original pinned interface laws.
Source resistivity determines current, native current supplies Joule heat, and
plastic plus reversible material heat reaches the original caloric cells once.
Bulk conductivity and contact resistance retain their original series owner.

The complete coupled balance includes kinetic change and separately reported
backward-Euler inertial loss, alongside original constitutive numerical loss.
Neither numerical loss becomes physical heat. Dynamic candidates retain velocity
and its native kinetic account; their `into_parts` returns thermal and dynamic
owners together, plus the original circuit where present. Continuation must use
all returned histories. Refusal, cancellation and source exit publish none.

Reference geometry, mass, electrical topology and transport areas remain fixed.
Support contacts retain their original complete, prescribed, compressive face
model; this adds no opening, impact, friction, evolving contact area, motional EMF,
Lorentz force or finite-strain mechanics. Loads are endpoint values and time
accuracy is first order. Focused regressions add independent consistent-mass
displacement/contact-pressure and two-cell heat oracles, source-only resistivity
changes, circuit and motion continuation, physical versus numerical loss,
pressure-dependent current and atomic source refusal. Native execution remains
pending in the resource-constrained build environment.

### Deforming bulk conductor paths on the native solid endpoint

`solid_enthalpy::plastic::electrical::geometry` binds explicit
`DeformingConductorPath` branch directions through
`ThermoplasticSolid::with_deforming_conductors`. The prepared heat network
supplies the sole branch-to-cell mapping and original whole-cell length/area;
the adapter accepts no competing element, material or Joule-heat route. Each
endpoint reuses the P1 displacement gradients already computed and admitted by
the native mechanical owner.

For reference unit direction n, the reduction forms F = I + grad(u), J = det(F)
and lambda = |F n|, then uses L = lambda L0 and A = J A0 / lambda. Positive
volume, the explicit gradient bound and the native material's small-strain and
rotation domain must all hold. The relative residual of
F^T F n = lambda^2 n must meet the declared orthogonality tolerance, so the
deformed path stays normal to its deformed end planes. A nonprincipal sheared
path needs a spatial conduction solve and refuses here. Directions are
normalized once at binding; branches must be unique and counts bounded.

The original `ResolvedConductorState` adopts the complete endpoint source
bundle and lowers the changed dimensions through `into_uniform_conductor`.
This occurs inside the existing electrical/enthalpy callback after the same
native mechanical proposal, including its inertial variant. Reference
dimensions always come from the prepared network, never from a previously
stretched resistor. Material source pins, charge/flux, plastic history and
enthalpy remain with their original owners. Original circuit Joule heat is
deposited once; a changing purely resistive coefficient adds no separate field
storage, parameter work or heat source.

This is a bounded affine geometry correction to the existing fully solid,
isotropic DC, small-strain uniform-path reduction. Reference inventory and
caloric geometry, other thermal/contact transport areas, topology and C/L stay
fixed. It does not infer piezoresistivity, motional EMF, Lorentz forces, melting
or finite-strain constitutive mechanics. Endpoint material/mechanical evaluation
with midpoint current retains the existing generally first-order coupled
accuracy. Refusal publishes no partial electrical, thermal or mechanical state.

Four focused synthetic regressions cover independently evaluated affine
resistance/current/heat with unchanged inventory, native consistent-mass axial
motion over two steps without accumulating geometric stretch, thermal-expansion
data feeding shape and current at one endpoint, and nonprincipal or
out-of-domain atomic refusal.

### Sourced thermal expansion in the native thermoplastic solid

`ThermoplasticSolid` accepts optional source-integrated isotropic expansion
through `SourceQuery`. Each native element receives its actual temperature and
`integral(alpha(T) dT)` stress-free strain. The original P1/J2 solve owns
compatible displacement, support reactions and plastic history. The material
owner supplies physical yield heat and signed reversible exchange to the same
enthalpy graph, including the existing coupled contact and electrical paths.
Reports retain actual elastic/hardening free energy, `T s_mechanical` internal
energy, plastic heat, reversible exchange and numerical loss separately.
Globally affine sourced E(T)/H(T) use the material owner's entropy
`s=alpha(T) tr(sigma)-(E'/E) psi_elastic-H' p²/2`. The separately exposed
`coefficient_storage_change_j` corrects native endpoint-modulus element accounts
and is already included in `stored_mechanical_change_j`; it is never counted
twice or turned into plastic heat. Combined energy and actual nodal-work gates
must pass before acceptance, including the existing support-contact path.

Expansion claims remain pinned over the fully solid caloric span. Curves and
fixed-context surfaces may vary with temperature, including negative/zero
coefficients; the shared material owner admits a conservative positive
effective-capacity bound for the complete strain and temperature domain.
Reference mass, geometry and Poisson ratio remain fixed. E(T) and native H(T)
may be globally affine over the supported solid span; source-cell derivative
jumps, finite deformation and melting remain outside this model. Focused
native-mesh tests exercise data-only expansion and modulus swaps, free and
constrained heating, two successive loaded plastic increments, independent
shape/heat/free-energy oracles and unchanged inventory on capacity refusal.
Native execution of the new modulus regressions is pending in the configured
build environment.

### Shared structural and conductive material assignments

`solid_enthalpy::plastic::ThermoplasticSolid::thermal_material_cells` lowers
the existing per-tetrahedron bindings into the original material transport
network. It retains cell/material identity, the original caloric object and
fixed pressure/composition context, selects conductivity at accepted enthalpy,
and pins the admitted claim. The existing prepared observer evaluates that
claim at actual heat states. Source support loss refuses the complete coupled
interval; mechanically or frequency-bounded conductivity requires an evolving
response. Additional thermal cells need their own complete network bindings.
The original source solid remains the immutable borrowed context anchor while
accepted candidate histories evolve.

G1/G3 regressions compare two-tet mixed-shear equilibrium and conservative heat
exchange under data-only conductivity/Cp changes. They exercise a free
mechanical displacement, local plastic history through unloading, source
identity and context continuity, physical versus numerical heat, and G4 atomic
source-domain refusal. These synthetic numerical fixtures do not extend the
existing small-strain, fixed-density, fixed-geometry, single-solid model.

### Solved support pressure, electrical current and heating

`solid_enthalpy::plastic::electrical::ElectricalSupportContact` adds a pinned
electrical contact claim and explicit two-body heat split to an existing
whole-face thermal/mechanical support. The accepted-stress source initializer
and joint `advance_electrical_with_pressure_contacts_candidate` use the same
ordered interface card, reference face area, actual compressive traction and
both live material temperatures. The shared static response admission retains
declared environmental data; source-domain failure refuses the complete step.

The original circuit reduces R'' by area, solves the current and supplies
contact Joule heat to the original caloric destinations once. The original
thermal contact, native plastic/reversible heat and closed energy checks stay
in the same endpoint residual. Support bindings may not duplicate a prepared
interface correction or change the circuit topology during a step. Explicit
initial circuit rebinding preserves the native conserved charge/flux state.

The model remains a fixed-footprint, prescribed, compressive support with
zero interfacial energy storage. It does not infer opening, arcing, evolving
area, spreading temperature or electromagnetic force. Focused synthetic G1/G3
regressions compare independent RC/current, two-cell heat and scalar plasticity
equations, data-only interface swaps, source-pressure refusal and cancellation;
native execution is required before claiming those regressions pass.

## Assembled material shell bodies in the reciprocal string bank

`shear_shell::assembly::MaterialShellAssembly` joins independently sourced
component midsurfaces through the original shell mortar owner. Complete
component stiffness and mass are tied before one shared eigensolve; no part
is independently truncated or fitted to a desired body spectrum. Each mode
is recovered through the original component support/director mapping. Physical
bridge forces, full tensor viscosity, source-bound interface damping and the
independently selected exterior/cavity pressure rows use this same basis.

`ReciprocalStringBank::from_material_shell_assembly` retains the original
parts, material records, ties, physical ports and projected image in the
existing bank. All strings share one body and optional Helmholtz cavity.
Changing one component's material recomputes its actual mass, stiffness, loss
and thickness-dependent tie geometry while keeping the other sources intact.
Internal parts must explicitly select empty fluid surfaces when appropriate.

Current junction admission requires global support frames and refuses
prescribed slave coordinates that the original mortar elimination cannot
compose. Rigid ties do not provide compliant glue, assembled prestress or
nonlinear contact. Pressure remains the inherited compact observer and single
Helmholtz closure; a closed skin or resolved interior/exterior air is not
inferred. Focused regression sources cover common modes and work, independent
fluid roles, full summed viscous cross terms, material substitution and original
performer retention. Their native numerical execution remains required.

### Source data for assembled string bodies

The existing `material_strings --body-profile` route accepts v4 component
assemblies through the shared v3 component reader and material cache. It
retains explicit directors, regional packs, support frames, material interfaces
and independently selected fluid surfaces per part, then supplies physical
edge/face ties and per-string bridge footprints to the common assembled owner.
One whole-body window and projection allowance apply after component assembly.
Empty fluid roles are explicit, individually free parts may be tied, and at
least one region retains the performance's original body source. Source-bound
bow interfaces, torsion, physical rocking levers, preloads and the original
receiver remain in the existing performer. The bounded grammar, synthetic
two-component fixture and unqualified native regression sources are documented
in `examples/material_strings/BODY_PROFILES.md`.

### Shared material keys in native reed performance inputs

The material-modal valve file consumes the same material-boundary grammar and
DKT/contact/Maxwell/slit builder as flue files. Explicit upstream/downstream
driven chimney nodes are completely matched to source-card keys using the
actual lowered regional gas, impedance, clock and lifetime. Blades and keys
share one immutable path cache and the original artifact byte allowances.

`control_source` binds all-v2 pressure and nonnegative closing-force tracks
through the original gesture sampler. The legacy pressure suffix retains
explicit initial key drives. The canonical two-sided observer and boundary
programs apply assignments at every mechanical substep; cursors publish after
physical acceptance and receiver ingestion. Both key sides, source forces and
pressure share stored-assignment limits, with combined auxiliary compilation
work. Releasing a force or changing callback partition keeps the original
material, contact and acoustic histories. Imported/computed MIMO and compact
outlet reactions remain installed in the original owners.

The CLI retains each key's side, physical node, original material artifact,
typed consumer and local gas. `examples/material_modal_valve/paired-keyed.performance`
demonstrates the shared path. Focused regressions cover source-to-pressure/PCM
changes, mechanical substep controls, one/two-blade sources and the three
runtime source families. Native execution remains required. Keys remain
one-mode material cantilever slits with prescribed exterior reservoirs; this
does not add circular pads, linkages, direct key radiation or instrument
qualification.


### Shared material-card relaxation for flute and reed duct walls

Both wind-performance readers use one bounded wall-loss selector.
`maxwell COUNT MAX_DT_OVER_TAU MAX_ANGULAR_STEP MODULUS_PROPERTY TIME_PROPERTY ...`
selects 1..=8 original modulus/time claim pairs from each wall's normalized
material card. Density, equilibrium Young modulus and Poisson ratio come from
that same resolved state. Every coefficient must have its physical units,
fixed-temperature constancy and frequency coverage; duplicate identical
branches and unresolved source times or instantaneous wall frequencies refuse.

`WallRelaxation::from_state_point` lowers those claims through the existing
generalized-Maxwell and passive pressure/volume-flow owners. The original
two-sided reed and flute network equations carry wall memory, returning
pressure and radiation without an extra wall integrator or Kelvin-Voigt loss.
Distinct branches may share a time claim; each selected property is resolved
once and all source branches remain visible.

Regression sources change only relaxation-time material data while holding
geometry, gases, equilibrium stiffness and performance bytes fixed, and check
wall-memory energy, source feedback, actual microphone pressure, callback
partition equality and the flute CLI's PCM output. These are synthetic
coupling checks. The inherited fixed-state, isotropic thin-cylinder membrane
limits remain; this does not add orthotropic bore shells, thermal evolution,
exterior wall radiation or measured instrument qualification.


### Reciprocal tangential receivers for material bow contacts

`MultiFrictionModalSystem::with_surface_attachments` binds an existing
receiving component through its physical displacement row. Relative slip and
equal/opposite traction use the same signed column in the original coupled
solver. The complete contact compliance matrix is rebuilt, including reactions
between contacts sharing a receiver; the single-contact path uses that same
map. Initial motion, material histories, source interfaces, pressure observation
and existing budgets remain with their original owners.

`MultiBowStringBank::with_bow_body_surfaces` accepts original body indices and
retains the source declarations through scheduled playback. Body maps must come
from the same geometry/material basis as the receiver dynamics. Bind before
playback and normal obstacles; a contact cannot be rebound. A resolved normal
and tangential counterpart must identify the same receiving body. Prescribed
gesture speed remains additional actuator motion and contributes external work;
zero speed makes the tangential interface passive.

Focused regression sources cover an analytic two-mass momentum/work balance,
shared-receiver cross-reactions, material-density changes reaching traction and
actual body pressure, and direct/scheduled callback equality. This is an
explicit mechanical connection, not a supplied bow-hair construction, rosin
thermal model or measured violin/cello/bass qualification.


### Independent material meshes on a common double-reed interface

Opposing material-modal blades may retain different node counts, triangles and
slit subdivisions when they cover the same simple planar footprint and complete
physical slit in the same coordinate frame. The first blade remains the
declared clearance, contact and slit integration surface. The second blade's
P1 displacement trace uses the shared `fs_plate::loading::PlateTraceStencil`;
the identical transpose returns contact forces into its original modal basis.
Original pressure faces, swept-volume rows, full viscosity and material
histories keep their existing owners. Identical meshes preserve the original
columns exactly.

Complete oriented boundary and slit coverage is checked before interpolation.
Missing coverage, holes, disconnected sheets, duplicate or degenerate cells,
nonmanifold topology and crossing boundaries refuse. Source budgets also bound
trace construction; there is no nearest-node snap, extrapolation or inferred
frame alignment.

Regression sources cover rigid/affine traces, virtual-work and dissipation
balance, independent mesh refinement against analytic contact energy in both
blade orders, and independently refined material blades through the scheduled
gas-loaded microphone. The P1 trace and chosen quadrature remain finite
approximations: a coarse interface is not invariant to swapping blade order,
and source/interface refinement must be assessed separately. Curved or
noncoincident blades, moving contact footprints and measured reed fidelity
remain outside this connection.

### Material-surface feedback in hereditary vibration

The existing `thermovisco_dynamics::Oscillator` consumes source-backed Maxwell
laws built from complete fixed-context temperature lines, including sampled
temperature/pressure surfaces. Changing the pressure query on an identical
material card changes the actual hereditary motion and physical dashpot heat;
the same heat reaches the original lumped enthalpy-feedback owner. Focused
regressions check the thermal endpoint and independent combined energy balance.
Source continuity, role pins and shared material identity remain enforced.
Non-temperature coordinates stay fixed external conditions; this adds no
evolving pressure/composition work, phase transition, variable reference mass
or alternate caloric model. Numerical dissipation is not converted to heat.

### Environmental material surfaces in thermoplastic conduction

The existing native tetrahedral plasticity solve and prepared material heat
network consume the same fixed-context source surfaces. `ThermoplasticSource`
retains every crossed temperature-line breakpoint, original role pin, material
identity and typed external condition. At fixed pressure, Poisson ratio and
reference density remain constant while native E/H may be globally affine in
temperature with explicit coefficient work and mechanical entropy. Yield and
heat capacity may vary over their supported source paths; conductivity resolves
at the actual thermal endpoint.

The spatial regression changes only the pressure query on identical material
cards and checks the resulting free-node displacement, plastic history,
temperatures and conduction against independent mechanical and thermal
reductions. Physical plastic heat, stored elastic/hardening energy, boundary
heat and numerical loss remain separately balanced. Proposals retain atomic
mechanical/thermal acceptance. This remains the existing small-strain,
single-solid-phase, fixed-reference-mass and fixed-geometry model; it does not
supply evolving-pressure work or a new constitutive or heat-transfer solver.

### Source-bound torsional material memory

Original string sources may bind independent per-layer shear spectra through
`with_torsional_relaxation` after sourced torsion. The shared material section
owner supplies each `delta_G * J` and relaxation time; the existing exact owner
projects `delta_GJ * (n*pi/L)^2 / sum(rho*J)` into the original twist coordinates.
These histories follow original bending histories and precede body and explicit
bending additions. They survive gas, torque and axial rebuilds and act through
the existing surface-friction and rotational-bridge maps, with no direct twist
acoustic gain. Existing shear Kelvin viscosity remains a separate parallel loss.

Complete instantaneous frequency support, positive-arm `dt/tau <= 2`, original
memory/setup ceilings and cancellation gates apply. Unsupported scalar/impact
conversion refuses instead of dropping history. This is frozen small-angle
torsion of bonded circular sections; it does not add axial Prony loss, winding
slip, heating or measured-alloy accuracy. Focused source and runtime regressions
cover layer subdivision, independent complex shear response, body/pressure
coupling, history preservation and explicit zero-arm reduction.

### Shared material walls on regional duct paths

`RegionalDuct::with_material_wall` binds one `LocalMaterialWall` to input
impedance, stable reflection, characteristic feedback and physical baffled-outlet
pressure. The local owner retains the resolved material identity, thickness,
axial restraint and full elastic, Kelvin–Voigt or Maxwell law. The existing graph
material-wall reduction uses this same owner. A source-resolved relaxation law
must match the exact material snapshot; matching equilibrium modulus alone
cannot admit a foreign source. Explicitly authored spectra retain their declared
unbound semantics.

Every actual cone, profile and continuous-gas quadrature radius and admitted
chimney is reduced through the original thin isotropic membrane map. Maxwell
frequency response calls the existing `GeneralizedMaxwell::modulus`, retaining
positive-arm geometry admission and the graph RC image. Every sampled frequency
must lie inside material support; unsupported data never becomes a rigid wall
or a constant-loss substitute. Existing constant-pin wrappers retain their
original arithmetic and compatibility path.

This remains a frozen-state, local thin-membrane image with thickness/radius
at most 0.05. It adds no axial shell coupling, wall thermal evolution, pressure
prestress or exterior wall radiation. FIR sampling and history admission remain
with the existing characteristic/outlet realizers; graph midpoint-state
timestep limits are not reinterpreted as frequency-domain physics. Synthetic
source-card regressions compare graph/continuum impedance, the independent
elastic-cylinder limit, joint gas/geometry/material refinement, source identity,
and material-only relaxation changes in feedback and microphone impulse pressure.
They do not establish measured-material or experimental-wall validation.

## Live strings in a finite source-backed gas

`LinearMaterialStringRuntime::try_new_with_gas` and both runtime `rebind_gas`
methods admit the shared complete `GasState` through its material-owner
coherence check. Source-backed gas composition, pressure and transport affect
actual modal aerodynamic drag and the compact pressure observer. Rebinding
retains material receipts, modal displacement/momentum, accepted acceleration
history and sample position; it advances one epoch and transfers no inferred
solid heat or elastic work. Existing `AmbientGas` entry points use the same path.

`acoustic_realize::finite_gas` composes the existing thermal string and
`fs_material::thermochemical_gas::energy::IsochoricGasParcel`. One sample evolves
mechanics at the initial gas/material laws, deposits material viscous loss into
the solid and aerodynamic loss into the gas exactly once, and reconstructs both
states with their existing source owners. Final damping and observer
coefficients use the accepted endpoint gas, with the previous accepted
acceleration retained for the observer difference. There is one accepted sample
and epoch. Caller-prescribed solid and gas heats remain separate from the two
internal dissipation transfers.

The returned combined residual measures actual modal, solid-enthalpy and
gas-internal-energy changes, applied mechanical work, both external heats,
authored loss and the explicit damping-partition roundoff term. It does not
define loss as an energy deficit. Absolute original gas/solid solve allowances
and arithmetic roundoff bound admission. Gas source-domain, timestep-temperature,
pressure-observer, stale-epoch and cancellation refusals publish neither owner.
The initial runtime medium must equal the parcel's complete snapshot; gas mass,
volume, composition and source receipts remain with the immutable source-backed
parcel model.

This is first-order partitioned feedback with a well-mixed fixed-volume local
gas inventory. Its volume does not create acoustic chamber boundaries.
The original no-exchange step omits convective/conductive solid–gas transfer.
Gas motion/pressure work, radiation loading and propagation of a changing
atmosphere remain outside both finite-gas paths. The
compact pressure observer is diagnostic, not an additional heat debit.
Temperature-increment limits require temporal-refinement assessment and do not
certify acoustic or thermodynamic accuracy. The same thermal-string admission
applies: either constant stiffness/fixed tension or sourced fixed-support
thermoelasticity, with constant density, solid phase and a fixed modal basis.
The combined account includes axial storage and mechanical entropy for the latter.

G1/G3/G4 regression definitions in `tests/string_specimen/finite_gas.rs` check
source-dependent drag without a vibration reset, measured combined energy,
independent endpoint pressure with the new sound-speed coefficient,
frozen-initial-medium staging, and atomic refusal. Synthetic source data test
numerical coupling, not experimental gas/material fidelity.

### Reacting finite gas on a continuing material string

`acoustic_realize::finite_gas::reactive` feeds the native string sample's
aerodynamic dissipation and prescribed gas heat into the existing
`GasReaction::heated_implicit_candidate`. Every declared reaction rate sees
the same source-backed endpoint temperature and composition. The chemical
interval is exactly one native audio sample; there is no second physical step,
reaction-only subcycle or inferred mechanism. The original implicit mechanism,
caloric, species, consumption, entropy and work limits remain authoritative.

Formation energies stay inside the gas internal-energy datum and are not
added as an external heat release. The original finite-gas publisher measures
actual mechanical, solid and gas storage, rebuilds drag and compact pressure
from the final mixture and preserves the preceding accepted acceleration.
Gas, thermal string, modal history and sample position publish at one epoch.
Source, nonlinear, temperature-limit, energy, observer and cancellation failures
leave both accepted owners unchanged and permit an identical retry.

The insulated reacting entry remains a first-order frozen-initial-mechanics split
around simultaneous heated chemistry in a uniform fixed-volume inventory. It does not resolve flame
or chamber propagation, pressure work, reciprocal conductive exchange or
acoustic radiation loading. The microphone remains a diagnostic observer.
G1/G3/G4 definitions in `tests/string_specimen/finite_gas/reactive.rs` check
independent isomer extent/caloric/pressure balances, formation-energy accounting,
continued mixture-dependent drag and sound, repeated species evolution and
atomic chemical/late-publication refusal. Synthetic sources establish numerical
integration behavior, not experimental chemical or acoustic fidelity.

The `reactive::exchange` entry points also close reciprocal thermal exchange
with the reacting gas. They use the SAME bounded exchange/thermal-proposal
owner as ordinary finite-gas conduction, replacing only the endpoint gas
proposal with the original heated implicit reaction owner. Each proposed Q
restarts from accepted species, supplies +Q to gas and -Q to solid, and evaluates
the interface using the final reacting temperature/composition. Chemical heat
is never a second source. The concrete conduction path shares solid source,
condition, geometry, conductivity-selection and pointwise Biot admission.
The original nested reaction/exchange/thermoelastic work ceilings remain explicit;
there is no additional time advance per trial. Original chemical and thermal
reports publish with final drag, pressure and the one accepted physical state.
Focused `reactive/exchange.rs` regressions independently solve both endpoint
equations, check composition-dependent conductance, insulated trajectory parity,
late publication refusal, identical retry and cancellation after gas proposal.

The `reactive::flow` entries compose those same heated chemical endpoints with
conservative nozzle exchange on the same audio interval. The ordinary flow
owner transports species and upstream enthalpy once, measures both finite gas
storage changes, and publishes the string, local gas and reservoir together.
The reaction/thermal report describes the pre-flow endpoint; final drag and
compact pressure use the transported mixture. The original reaction catalog
admits final local species, caloric/transport sources and explicit water identity
before publication, including newly incoming inert species. Unknown inflow
cannot strand the next chemical sample in an unsupported accepted state.
Formation energy and transported enthalpy are never additional heat credits.
This retains first-order heat/reaction-then-flow splitting and diagnostic sound;
it does not resolve a flame, jet momentum or reciprocal acoustic loading.

### Reciprocal thermal exchange with finite surrounding gas

The original `step_with_isochoric_gas` remains the insulated solid/gas API.
`step_with_isochoric_gas_exchange` adds one signed transfer Q from solid to
gas, using the existing bounded scalar feedback owner and concrete
`ConductanceSample` source evaluations. Each proposal advances the original
solid enthalpy chart by material loss + prescribed solid heat - Q, and the
original isochoric gas owner by aerodynamic loss + prescribed gas heat + Q.
The interface residual is Q - dt G(T_s,end - T_g,end). Prescribed heating is
included before choosing transfer direction; constant heat capacity and
averaged temperatures are not substituted for either source law.

The source callback receives actual endpoint states and must be repeatable,
proposal-only, and justified by its own geometry and validity. The concrete
`step_with_isochoric_gas_conduction` composes the existing
`MaterialConductionLayers` evaluator: G = A/(d_g/k_g + d_s/k_s), using the
accepted trial gas conductivity and the same immutable solid card/conditions
at its trial temperature. The declared area cannot exceed the actual string
surface. Every trial applies the existing 0.1 lumped Biot ceiling with the
declared transfer area, Bi = G(V/A)/(A k_s). Zero area or explicitly zero
conductivity means insulation; missing or dynamically supported source data
still refuses. Effective distances, initially uniform temperature and local
spreading remain explicit caller/model limitations. This is a pointwise Biot
gate, not a validity-wide spatial-gradient certificate.

The bracket contains zero and all evaluated proposals, including the
source-only probe, must stay within both original charts. Bounded iterations
admit actual source and thermal storage residuals. A final transfer that
crosses endpoint temperatures refuses. There is no uniqueness claim for an
arbitrary nonmonotone conductance, source extrapolation or hidden substepping.

Final material/medium rebinding, pressure history, one sample/epoch, combined
measured energy and atomic paired publication use the original finite-gas
path. Internal Q cancels from the combined account; the solid subaccount's
external boundary heat is prescribed heat - Q, while the frame separately
exposes prescribed solid heat, signed Q, concrete interface source and both
exchange residuals. Gas/solid source, geometry, Biot, iteration, energy,
temperature-increment, epoch, observer or cancellation refusals publish
neither state. Mechanics still freezes the initial laws for the sample:
this first-order split requires temporal refinement and does not add a
spatial gas layer, acoustic chamber, pressure work or radiation loading.

G1/G3/G4 definitions in `tests/string_specimen/finite_gas/exchange.rs` compare
the independent two-capacity backward-Euler solution in both directions,
including external heating that reverses the initial ordering; source/geometry
dependent conduction and damping; insulation parity with the original
vibration/observer path; independent reciprocal energy; and late refusal and
cancellation rollback. Synthetic data establish coupling tests, not measured
material fidelity.

### Shared environmental response of electrical and thermal contacts

`ElectricalContactNetwork` and the existing finite-cell thermal-contact adapter
use the same `ThermalInterfacePoint` response-axis admission. Original interface
curves and sampled surfaces may depend on declared evolving endpoint thermal
coordinates and explicitly typed fixed external pressure/composition conditions.
Electrical and thermal properties retain their own claims, dimensional meanings
and source support. Their common ordered interface does not equate the two
body temperatures or either body's full material query point.

The existing prepared electrical owner pins sources; original circuit and
enthalpy owners solve current, bulk/contact Joule heat and thermal exchange.
Pressure-dependent `R(T_A,p)` and `R''(T_B,p)` therefore affect the same trial,
with each Joule loss deposited once. Geometry, interface storage assumptions
and explicit heat partitions are unchanged. Focused regressions in
`tests/electrical_contacts/conditions.rs` compare pressure and composition
substitutions with analytic midpoint discharge and both body heat budgets,
exercise both live source surfaces through finite history, and refuse mechanical
or untyped response axes and surface-support holes without partial publication.
These synthetic source tests do not establish measured-material validation.

`PlateValvePerformance` also admits explicit `dg_tube` source data through the original exact DG exponential and its gas-derived Bessel memory, coupled to the retained ExactZoh plate/contact/Maxwell aperture by the exact step-mean port law. The same ambient gas drives the column, aperture density and slit viscosity; the existing wind/ensemble renderer observes internal inlet pressure only. Full-state DOF, finite-band Bessel and exponential admissions remain separate from spatial/nonlinear convergence, and inferred duct loss is not an independent heat certificate. See `examples/PLATE_VALVE_PERFORMANCES.md` for the bounded cylinder scope and source grammar.

### Fixed-volume chambers in the shared reacting gas assembly

`GasPiston::with_fixed_support` and `MaterialGasPiston::with_fixed_support`
attach an ideal stationary support to an existing wall at rest. The same owner
retains physical wall mass, area, material-backed gas and accepted volume;
fixed support never approximates a constraint with a large inertial mass.
The existing gas midpoint heat/volume solve runs at exactly unchanged volume.
`PistonStep` records the balancing signed support impulse, zero support work,
and the independently closed momentum and energy residuals. Attaching a support
to a moving wall refuses, since its stopping impulse and energy need an actual
accounted mechanism. Double-acting construction requires a free primary wall.

The supported owner is retained by existing reaction, nozzle/species, finite
heat, radiation/conduction and circuit compositions. Formation energy and
transport enthalpy still pass through their original owners exactly once.
Fixed support eliminates the wall-compressibility oscillator and its clock
restriction, while all gas, reaction, transport and thermal controls remain.
A fixed chamber accepts an exact `[V,V]` volume bracket. This provides
stationary reactors and reservoirs through the same assembly as moving walls.

Focused synthetic regressions in `piston/support/tests.rs`, reactive tests and
`tests/material_exchange/flow.rs` cover changing temperature/pressure/viscosity/
sound at fixed volume, independence from physical wall mass, formation-inclusive
reaction energy, a complete flow/reaction/thermal/circuit history, support work,
cancellation and atomic refusal. Source parsing was checked; native execution
remains a separate verification requirement.

### Spatial enthalpy with directional source-backed solids

`solid_enthalpy::TetCellMaterial` selects `ThermalStiffnessModel` and
`ThermalExpansionModel` independently. These are the same types and
implementations used by `fs-mbd-material`, owned by
`fs-solid::linear3::source`. Complete isotropic E/nu or nine orthotropic
engineering constants lower to the existing solid tensor. Scalar or three-axis
source expansion integrates from the original stress-free state, then rotates
once into the declared nonzero reference mesh frame. The constructor now
requires that frame; an equilibrium solve cannot relabel it.

Reference density and initial stiffness resolve atomically before source pins
are retained. During heating, only current elastic constants and expansion
are refreshed. The exact delivered mass divided by reference volume supplies
inertia and body loads; a density claim at the heated temperature is neither
required nor substituted. `TetMaterialSample` retains the current elastic
record, original reference record, both model choices, and the exact integrated
expansion/transform receipt. The identical accepted material and strain arrays
feed equilibrium, stress recovery and the existing prestress/vibration pencil.

The original substep/event observer admits every directional source path.
Heating outside one selected coefficient's support refuses the complete
candidate even if later cooling would return to a supported endpoint.
`tests/solid_enthalpy/directional.rs` covers rotated free motion, restrained
stress, data-only modulus/expansion changes in the native prestress tangent,
independent elastic and expansion axes, fixed reference mass with density
supported only at the original temperature, frame mismatch and heat/cool
rollback. These are synthetic analytic fixtures. The existing sequential
small-strain, fully-solid, one-way quasistatic approximation remains explicit;
caloric and elastic energy are reported separately, with no reciprocal
thermoelastic heat, changing transport geometry or experimental claim.

`SolidEnthalpyState::advance_transient_candidate` also feeds those same local
material/thermal arrays into the existing sparse 3-D Newmark trajectory. The
previous state must match the run's exact material and thermal fields and
clock. Both endpoints retain the same source-pinned reference records, cell
inventory, axes, frames and strain bounds. The next problem is checked against
the resolved endpoint before its source identities are borrowed by the native
candidate. Changed source stiffness rebuilds the native operator at fixed
reference mass; temperature-dependent softening and stiffening continue the
same displacement, velocity and acceleration history. Density changes and
changing stiffness with stiffness-proportional damping refuse. Constant mass
damping remains admissible.

Expansion work and constitutive-stiffness work are independently integrated
from the native stress/strain endpoints, including fully clamped storage, and
kept separate from caloric heat. The latter is the signed tensor term
`sum V/2 * elastic_strain_old : (C_next-C_old) : elastic_strain_next`.
This is sequential one-way thermal-to-inertial loading; no reciprocal heat
correction or whole caloric/mechanical energy closure is claimed. Tests in
`tests/solid_enthalpy.rs` change material expansion, modulus and heat capacity
as source data, compare actual heated inertial endpoints with a one-DOF
reference, continue the same motion through softening and stiffening, and
refuse missing source/clock bindings or unsupported damping. A nonzero initial
elastic strain makes omitted stiffness work exceed the energy allowance.
The endpoint adapter does not observe intervening thermal events: callers
must admit that path, resolve heating events and check temporal refinement.

`SolidEnthalpyState::rebind_transient_candidate` connects explicit source-card
replacement to the same retained 3-D motion. Both snapshots must originate from
their material/caloric owners and match the accepted clock, reference mesh/frame,
zero supports and persistent cell correspondence. The old snapshot must exactly
match the run; the new problem must exactly match the replacement snapshot. The
native owner transfers displacement and reduced canonical momentum, rebuilds
mass/stiffness/thermal forces, re-equilibrates acceleration and reports bounded
signed mechanical work. Subsequent ordinary heating retains the new sources and
delivered mass.

`SolidMaterialRebind` separately retains each cell's old/new mass, exact
phase/card/chart identity and extensive `m*h` inventory. Different caloric charts
may use unrelated energy zeros, so these values are not subtracted into heat or a
purported closed total-energy balance. The caller accepts the new thermal owner,
source binding, resolved state and run together and owns the external matter/energy
replacement. This is neither continuous mass flow nor a melting closure. Synthetic
G1/G4 tests in `tests/solid_enthalpy/rebind.rs` exercise changed density, modulus,
expansion and heat capacity through nonzero-time motion, reverse replacement,
subsequent heating, and exact binding/work/cancellation rollback.

### Finite source gas on uniform thermal feedback

`finite_gas_thermal::FiniteGasThermalEnvironment` adapts the original
`LumpedThermalEnvironment` to an `IsochoricGasParcel` through the existing
`UniformEnthalpyStep` callback. Each proposal retains the accepted body, gas and
hereditary history. The original gas energy inversion receives prescribed gas
heat + Q; the original lumped body solves against that trial gas endpoint. The
shared bounded feedback owner enforces Q + dt*convection_into_body_w = 0. No
heat-capacity, EOS, Biot or enthalpy solver is replaced.

Accepted body boundary heat is separately integrated enclosure radiation - Q.
The same Q enters the gas, while the original body's endpoint powers and the
source mismatch remain explicit in the report. Actual body and combined body/gas
storage residuals must each meet the declared absolute allowance. Radiation is
an external enclosure transfer, not a second gas heat credit. Discarded
material/thermal trials never debit the finite reservoir; callers publish the
returned gas only after their entire coupled candidate accepts.

The original source, whole-curve Biot, chart and cancellation refusals remain.
Film coefficient applicability and a well-mixed, fixed-volume gas inventory are
explicit assumptions. No flow, pressure work, reaction, spatial gas gradient or
acoustic loading is implied. Source data and numerical energy closure are not
experimental material validation.

Focused G1/G3/G4 regression definitions in
`tests/thermovisco_dynamics/finite_gas.rs` cover independent two-capacity exchange,
finite gas inventory changing actual hereditary vibration, exact insulation
parity, independent radiation ownership, measured combined energy and
failed/cancelled proposals preserving both accepted owners.

`advance_with_convection` uses the same bounded residual and enthalpy owner
while resolving a new whole-surface film coefficient from the complete gas
endpoint at each trial. `advance_with_convection_feedback` additionally supplies
the signed proposed heat to select a correlation's heating/cooling branch;
the accepted temperature ordering and energy gates still apply. The exact
source sample accompanies the final gas and body
reports. The declared-coefficient `advance` delegates to this path with a
constant callback. This is a bulk-gas coefficient law, not a wall-temperature
film iteration or resolved flow. The original Biot admission remains active;
callers must not dilute a local coefficient by averaging partial coverage.

### Shared DG causal terminal boundary

The source-bound DgPortLoad now also admits a separately prepared proper
rational outlet through the common fs-vfit MIMO admission. The inlet source
sees the Schur complement of the original two-port DG mean-pressure matrix
and the terminal's original mean-output map. This preserves the signed
cross-port terms and both accepted histories without a delayed feedback
sample or another bore/filter solver. Existing ideal-open/closed/matched
loads retain their original propagation path.

Both aperture and flue coupling validate inlet work independently of the
duct's net all-port work. DgTerminalFrame reports actual outward held flow,
mean pressure, signed DG/rational work and their compatibility residual.
Terminal exchange is not labeled as heat or stored energy. The existing
dense DG inferred-loss limitation and the sampled rational GridOnly report
remain explicit. Each owner is exact for held inputs; clock refinement of
the coupled mean solve remains necessary and is not a joint continuous
exponential certificate. Failed or abandoned trials publish neither owner.

Focused tests compare original DG and rational trajectories/work over
retained history in both propagator kernels, exercise energy returned by a
reactive terminal, preserve every history on refusal, and account for
terminal exchange in the moving-aperture balance.

### Distributed material-shell cavity in the reciprocal string bank

ReciprocalStringBank::with_material_shell_cavity attaches the original
tetrahedral Helmholtz pressure basis to every retained displacement mode of
a material shell or assembled shell body. Every declared cavity facet must
exactly match one boundary triangle with outward-from-gas ownership. The
original FE quadrature integrates the complete linear structural trace and
pressure shape. No hole filling or fitted body resonances are introduced.

The existing CavityCoupling::modal_network supplies all reciprocal storage
and connections in the original bank solve. The uniform mode is the exact
normalized constant and contributes compression stiffness without a
fictitious inertial coordinate. Every standing-wave mode in the declared
window must fit the explicit pressure, physical-state and work budgets.
The shared pressure ceiling is now 256, separate from the total mechanical
budget and unchanged compact-neck allowance; excess modes refuse.

The source retains geometry, gas, body basis identity, numerical declarations,
pressure shapes, measured modal report and overlap rows. Material histories,
body observation, bridge and torque links, unsteady string-air components
and original exterior body loading survive the cold attachment and existing
rebuild paths. Distributed air participates in material-interface frequency
admission. Interior pressure is exposed separately from microphone pressure.

This path covers one connected closed conforming midsurface cavity with
frozen coherent gas, linear structural kinematics and lossless interior
acoustics. It refuses open/unmatched geometry, stale gas and duplicated
lumped air compliance. It does not supply f-hole/exterior exchange or
thermoviscous cavity wall loss, and does not establish measured instrument
accuracy. Mesh, pressure-order, modal-window and time convergence remain
necessary. Focused source tests cover body/string reaction, exact uniform
compression, material and atmosphere changes, physical refusals and both
orders of cavity/unsteady-string-air attachment.

## Distributed winding material, circuit and enthalpy coupling

`circuit_enthalpy::winding::WindingThermalNetwork` connects the existing
`fs-em::winding::ohmic::OhmicWindingCircuit` to the same finite enthalpy
Newton/backward-Euler owner used by material resistor and field-motion circuits.
It reuses `MaterialCell`, `ConductionLink` and `PreparedMaterialNetwork`;
there is no second material evaluator, circuit integrator or thermal solver.

Every original magnetic triangle explicitly names an existing thermal cell or
a `MagneticRegion::Fixed` prepared isothermal state. A current-carrying triangle
must name a real receiving cell. Mapped triangle volumes cover each whole cell,
including the original extrusion thickness; the source density and cell mass
must agree with that fixed volume within the declared absolute allowance.
Mapped magnetic cells retain the actual borrowed caloric law and its material
card, and admit fully solid, constant-density charts. Thermal cells outside the
magnetic geometry retain the existing latent storage, source and boundary
behavior. This is an authored geometric correspondence, not an inferred remap.

Conductor resistivity or conductivity, caloric temperature and source thermal
conductivity come from the same immutable cards and actual proposed endpoint.
Electrical claim identities persist with the accepted circuit; thermal
transport selections stay pinned through each interval's nonlinear trials.
Canonical complete state bundles are refreshed through their existing
`at_query_point` owner. All nonthermal coordinates remain fixed. Thermally
mapped nonconducting cores also refresh their original magnetic sources.
Changed permeability or a temperature-dependent magnetic source whose
temperature moves refuses, even when endpoint values coincide. Explicit
Fixed regions retain their prepared full sources. All magnetic source
snapshots remain available in the accepted report.

Each nonlinear trial starts from identical accepted flux/charge/time. The
native descriptor midpoint supplies the mutually coupled winding currents and
all spatial and graph-resistor losses. Every original numeric graph resistor
requires its own explicit heat destination; its authored resistance remains
unchanged. The private thermal-driver routing boundary deposits both kinds of
loss exactly once into the existing local residual. Source-derived link
conduction and authored boundary heat enter that same backward-Euler solve.
Cell residuals, complete deposited Joule heat, native electrical work and
combined thermal/capacitor/magnetic energy undergo fresh final admission.
Budget, source, topology or cancellation failure publishes neither endpoint.

This is generally first-order coupled accuracy: endpoint material coefficients
with electrical midpoint currents and thermal backward Euler. Source and
topology events require interval cuts; target events remain pending. Geometry,
inductance and externally fixed magnetic regions remain fixed. No moving
winding, thermal expansion, melting conductor, magnetic parameter work,
hysteresis, skin/proximity loss or unrepresented winding-end resistance is
inferred. Regressions cover actual data-only conductor swaps, nonlinear
temperature feedback and continuation, mutual current evolution, source
conduction, all graph/coil heat, caloric/geometry identity, complete passive-core
mapping and fixed-state/source refusals, native budgets/events and late
cancellation. Independent numerical energy checks supplement source review;
native build/regression execution remains a separate requirement.

## Electrical heating of a live string in finite surrounding gas

`acoustic_realize::finite_gas::electrical` connects the existing source-bound
Ohmic drive to both finite-gas string steps. The original conductor resolves
initial resistance, current and Joule work; that work enters prescribed solid
heat exactly once. The existing owners advance vibration, material loss,
aerodynamic drag, solid enthalpy, gas energy/EOS, optional reciprocal layer
conduction and the compact pressure observer. Endpoint resistance and current
must also resolve before either accepted owner publishes. Late source refusal,
stale epoch or cancellation preserves both states and acoustic history.

`FiniteGasElectrothermalStringFrame` retains both source-resolved conductors,
initial and endpoint currents, Joule work and separately prescribed solid heat.
Its nested combined energy account already includes Joule work; callers must
not add it again. Material viscous loss heats the solid, aerodynamic drag heats
the finite gas, and signed interfacial heat cancels between their inventories.
Numerical loss and diagnostic radiation are not converted to extra heat.

The existing first-order frozen-coefficient split, fixed geometry/density,
fully solid caloric state, source validity, layer geometry and Biot restrictions
remain in force. The sourced thermoelastic constructor also propagates heating
into tension/stiffness and its reversible work through this same path. This
supplies no circuit dynamics, contact/skin effects, Lorentz force, gas pressure
work or radiation loading. Gas volume specifies a thermal inventory, not an
acoustic chamber geometry.

The three focused `string_specimen::finite_gas::electrical` regressions exercise
independent current/voltage work and signed-drive parity, source resistance and
damping feedback into actual pressure, cumulative mechanical/solid/gas storage,
and late refusal with deterministic retry and cancellation. Their manufactured
coefficients establish integration behavior, not measured material validation.

### Endpoint electrical work on the same material/gas trajectory

`finite_gas::electrical::feedback` adds bounded backward-Euler source feedback
to the original conduction transaction. Each private trial starts from the
same accepted string/gas history, deposits work J through the shared electrical
proposal owner, and resolves the original conductor at the thermal endpoint.
The existing safeguarded scalar solver closes J - dt I_end^2 R_end = 0 with an
explicit nonnegative work bracket, residual tolerance and iteration cap. Held
voltage uses the original source's I_end = V/R_end. Inner thermal/thermoelastic
work limits and material domains remain authoritative at every trial.

`FiniteGasElectricalFeedbackFrame` distinguishes deposited heat from actual
endpoint boundary work and retains their signed residual. The outer energy
audit uses actual source work and includes the declared residual allowance;
the nested original thermal account uses deposited J. No state or source is
renormalized, and neither quantity is an extra heat credit beside the other.
The nozzle variant transports species/enthalpy once AFTER convergence, over
the same sample interval. Its common two-volume balance includes the admitted
source-work correction before all three owners publish atomically. A failed
bracket, work cap, source domain, final flow/observer check or cancellation
preserves the accepted trajectory and input inventory for retry.

The initial-law APIs remain available with their original source quadrature.
Endpoint feedback changes electrical/thermal time discretization, not the
frozen-initial mechanics or first-order heat-then-flow split. It supplies no
finite circuit storage, induction, uniqueness or global-stability theorem.
Focused regressions in `tests/string_specimen/finite_gas/electrical/feedback.rs`
cover independent affine-resistance current/voltage roots, two-capacity thermal
exchange, single-transfer composition, zero-flow parity and atomic refusal.

## Finite enthalpy heating through a pressure-loaded solid valve

`solid_enthalpy::valve` binds each accepted P1 material/caloric state to the
existing thermal pressure-valve reducer. `advance_pressure_valve_to` advances
the original source-observed thermal timeline, refreshes pinned elastic and
expansion claims, solves compatible equilibrium once, then admits displaced
aperture geometry, complete tangent stability and pressure/modal compliance.
It returns the thermal state, caloric audit, resolved source state and valve as
one candidate. Reference mass, mesh frame and the original dead-pressure
footprint remain fixed. Source, phase, geometry, reduction or cancellation
refusal leaves the accepted input intact.

`examples/heated_valve_sound` accepts normalized material packs or explicitly
manufactured coefficients and carries finite heat through ordinary material
acoustic assembly to physical-pressure CSV and stereo WAV. The caller declares
a fully solid temperature interval; the example does not infer a melting law
from heat capacity. Its cold and heated outputs are frozen-state experiments.
They do not claim live modal-state transfer, reciprocal thermoelastic heat,
moving load geometry, measured damping or instrument validation. The thermal
audit remains caloric rather than a claim of reciprocal thermomechanical work.

Three core regressions compare an independent six-tetrahedron tangent and
pressure compliance, source/aperture refusal and cancellation. Two native
example regressions check finite-energy temperature, fixed reference mass,
data-only heat-capacity changes reaching actual microphone pressure,
deterministic replay and refusal without a sound comparison.

### Native radiative heat on the same compatible solid state

`solid_enthalpy::spatial::SpatialSolidEnthalpy::advance_with_radiation` takes
the native borrowed `EnthalpyRadiation` selection. Ambient reservoirs and
closed reflected exchange, using either mean patches or pointwise P1 emission,
therefore change the original nodal enthalpy, pinned material response,
integrated thermal strain and compatible native stress/deformation together.
The result retains `EnthalpyHeatStep`, including complete physical heat audits
and separate convection/radiation accounts. The original `advance` result and
its arithmetic remain unchanged; the candidate's defaulted thermal-result type
also carries the complete radiative variant without discarding its evidence.

All corners must remain fully solid and within source support before centroid
mechanical quadrature. A failed thermal, source, strain or mechanical admission
returns no candidate and leaves old enthalpy untouched. Fixed reference mass,
geometry and one-way thermal loading remain the model; this adds no reciprocal
elastic heat or acoustic radiation reaction. Focused G1/G3 tests compare both
ambient discretizations to a scalar Stefan-Boltzmann backward-Euler balance,
both closed-enclosure discretizations to two finite-body exchange equations,
then verify integrated expansion and clamped stress from those same states.
They also check default equivalence and late geometric/source refusal.

### Native nodal heat on the same continuing modal trajectory

`solid_enthalpy::spatial::motion::SpatialHeatedValveMotion` joins the native
heterogeneous P1 enthalpy owner to the same `HeatedValveReduction`,
`ThermalValveModalMotion` and `fs-modal` physical state transfer used by the
finite-cell path. The shared reduction helper preserves the original support,
initial-stress, aperture, pressure-compliance and modal admissions. Local
enthalpies remain on the original heat mesh, and every vertex is admitted
before the existing enthalpy-centroid mechanical quadrature. Reference mass
does not change with the queried hot density.

Each interval advances held-pressure motion, native backward-Euler heat,
source-dependent equilibrium/tangent, shape-based displacement/momentum
transfer and the new basis's pressure observer. Heat and motion publish
together only after all admissions and a final cancellation poll. Late source,
reduction, observation, transfer or budget refusal retains the accepted heat
field, vibration and sample clock. No independent conductivity, thermal store,
modal integrator or mode-number state matching is introduced.

`step_with_radiation` admits the same explicit native ambient/enclosure and
mean/pointwise heat owners as the quasistatic spatial path. It shares the
original continuing-motion transaction and returns the complete selected heat
report. Radiative heating refreshes source stiffness, integrated strain,
equilibrium and tangent modes before shape-based state transfer and the new
pressure observation. Both accepted heat and the modal sample clock remain
unchanged after a late observation or transfer-budget refusal. The focused
motion radiation regression independently integrates P1 fourth-power emission,
checks native mass/storage, follows the analytic changing-stiffness oscillator
and velocity-derived pressure, varies only emissivity, and checks exact retry.

This is a first-order one-way thermal composition on fixed reference geometry,
supports and pressure footprint. Signed parameter work and discarded modal
state are reported separately; neither becomes heat. Spatial/time/modal
convergence, reciprocal thermoelastic heat, moving geometry, melting and
acoustic radiation back-reaction remain separate obligations. Focused manufactured
tests in `spatial_solid_enthalpy::motion` compare an independent affine-tet
oscillator, unchanged-temperature continuation, heat-capacity sensitivity,
and late refusal with exact retry. These are implementation checks, not
measured material or instrument validation.

## Native continuously heated string and finite-gas journey

The `finite_gas_string` binary and example share one implementation. A bounded
scenario declares geometry, force/release, observation, current or voltage,
prescribed heat, thermal inventory and solve controls. A normalized solid claim
pack supplies constant density/modulus, temperature-dependent Kelvin-Voigt loss,
heat capacity, conductivity and optional electrical resistivity. The shared
native gas loader retains explicit species, NASA-9 and transport associations;
the finite parcel requires a declared frozen-gas policy.

The same modal state and acceleration history continue through all accepted
samples. Material loss heats the solid, drag heats the gas, and source-bound
layer conduction exchanges equal and opposite heat. Resistance, damping, gas
EOS and the compact pressure observer follow their accepted material states.
Every CSV row measures the cumulative mechanical/solid/gas storage change from
actual endpoints and checks it against source work and explicitly accounted
loss. Metadata retains card/curve/source identities and initial mechanical
energy. The whole physical trajectory is admitted before output; WAV creation
requires a new destination and uses the declared pressure scale without peak
normalization. Clipping is reported. I/O errors are not transactional rollback.

The single-phase interval is a caller declaration, not discovered melting
physics. All original fixed-geometry/tension, constant-E/density, first-order,
well-mixed gas and diagnostic-radiation limitations remain. The supplied solid
fixture is manufactured. Native regressions cover replay, independent storage
integrals, data-only material/gas/volume changes, strict input and source or
cancellation refusal; they do not establish measured instrument accuracy.

## Thermal material field storage

`circuit_enthalpy::field_storage::step_candidate` closes sourced C(T) and L(T), native Kirchhoff charge/flux, resistor heating and the original finite enthalpy/Newton residual on one candidate interval. Every field branch maps to one whole material cell with the same immutable card, source context, caloric law and geometry. The prepared material network keeps resistor/conduction source choices pinned. Field coefficients are refreshed at each proposed temperature through `fs-material::field_storage::thermal`; the original circuit owns rebinding and midpoint evolution. No replacement circuit or heat integrator is introduced.

For each field, the source law supplies psi=z²/(2K), s=-partial_T psi and U=psi+Ts. The frozen-initial-coordinate coefficient change supplies explicit parameter work to the native electrical step. Its negative, together with minus delta(Ts), becomes signed reversible transfer into the same material's caloric store; Joule dissipation is routed separately, exactly once. The final independent enclosing balance uses actual U at both endpoints plus caloric storage, against only authored electrical work and external heat. Original native reports retain their isothermal free-energy convention and expose the separate physical account explicitly.

The bounded reduction requires fixed geometry, complete storage mapping, one field branch per thermal material cell, fully solid constant-density caloric charts for mapped fields, positive affine coefficients and positive total caloric plus field heat capacity across each declared temperature/coordinate envelope. Both endpoint charges/fluxes must satisfy that envelope. Other thermal cells keep their original phase and boundary behavior. Source constraints, material-domain errors, unresolved thermal events, failed residuals, budget exhaustion and cancellation return no live-state mutation. Ideal-source constraint incompatibilities during frozen-coordinate rebinding refuse; no impulsive source exchange is inferred. Hysteresis, dielectric relaxation, moving geometry, field-dependent constitutive response, spatial field reassembly and phase transitions of mapped fields are separate owners. This first-order endpoint-temperature split does not claim unconditional passivity or temporal accuracy.

Focused authored RC/RL fixtures cover both thermal coefficient signs, exact constant-coefficient recovery, independent total-U scalar roots, repeated accepted history, Joule versus reversible versus external heat, positive closed-system entropy, source/volume mismatch, capacity envelope, budget and late cancellation. These are numerical regressions, not measured material validation. Native Cargo/DSR execution was unavailable for this increment.

### Implicit conservative homogeneous reaction extent

`cavity_phs::piston::reactive::implicit::ImplicitReactionSolve` selects
backward Euler for one reversible source-bound reaction, including inerts and
initially absent products. It solves `xi - dt V r_end = 0` in moles with the
same bounded scalar feedback algorithm used by thermal exchange. Each trial
starts from the same accepted gas mass, volume and TOTAL internal energy,
reconstructs the original NASA-9 mixture, and re-evaluates canonical partial
activities, the actual temperature-dependent material rate and thermodynamic
reverse closure. No formation-energy heat is added. Source/model identities
and explicit/implicit paths share the original inventory reconstruction and
final species, mass, energy and entropy admission.

The caller's signed molar bracket is intersected with current donor allowances;
accepted species are never clipped. Search endpoints retain strictly positive
donors. Zero-rate states accept before querying the nonzero edge. Every queried
state must remain inside the declared caloric/transport/kinetic source domain.
Rejected bracket samples need not pass the accepted-trajectory entropy gate;
the converged physical endpoint must. Reports retain endpoint rates, actual
nonlinear residual, admitted extent bracket and bounded trial count. Thermal
feedback's joule API and numerical operation order remain unchanged; the shared
internal scalar core attaches units only in its physical wrappers.

Existing parcel, piston, flow/diffusion and finite electrical/thermal consumers
select this method through the same `GasReaction` value. Existing dyadic
subcycling also refines implicit trajectories and may retry donor, entropy or
nonlinear-root failures, while source, arithmetic and cancellation errors
remain hard refusals. The root cap is at most 128 interior trials plus two
endpoints per attempted substep; every gas inversion retains its separate cap.
An aggregate subcycling report does not mislabel its summed extent as one root:
individual endpoint equations remain in the actual accepted substeps.

This is a first-order homogeneous reaction-then-motion/heat split, not a
multi-reaction implicit mechanism, chemical-equilibrium acoustic derivative or
validated ignition/engine model. Multi-channel implicit construction refuses.
Temporal refinement and nonlinear-source branch admission remain necessary.
Focused regressions use synthetic balanced isomers with an independent analytic
backward-Euler limit and exponential relaxation, actual temperature-dependent
source kinetics, formation-energy pressure/sound response, material piston
motion, absent reagents, source/work limits and cancellation rollback. Native
workspace execution and experimental mechanism qualification remain separate.

## Finite moisture transport into native plate mechanics

The material-reality MR46 slice in `moisture` keeps local water inventory
separate from ambient relative humidity. `MoistureSlab` is a homogeneous,
isothermal, one-dimensional fixed-reference specimen with kg water/kg dry solid
as its sole evolving field. Its source query explicitly distinguishes absolute
temperature and dry-basis moisture. Reference density at the declared reference
moisture gives conserved dry skeleton mass; later density queries never create
or destroy that inventory. Every interval measures the water received at each
end and checks total water change against those signed transfers.

The effective diffusivity is required material data; no wood/material name
selects a number. Its original claim is pinned and re-evaluated locally.
Backward Euler uses harmonic interior-face diffusivities evaluated at the
accepted start and half-cell resistance at prescribed surfaces. This lagged
nonlinear coefficient is a first-order approximation, with explicit interval,
content-increment, true-residual, conservation and work limits. The tridiagonal
SPD operator uses the original `fs-solver::checked_cg` and original
`fs-sparse::Precond` interface, not a second linear solver. Its identity
term and nonnegative graph conductances establish SPD. Source, numerical
maximum-principle, increment and water balance admissions precede publication.
No water is clipped or silently renormalized. Sealed ends carry zero transfer;
prescribed dry-basis values are infinite reservoirs. The explicit FPL wood
sorption boundary calls the existing coefficient-card owner at matching
temperature and actual RH; it changes the surface condition, never every
interior cell at once.

`moisture::plate::MoisturePlate` maps that same inventory to the existing
rectangular DKT mesh: diffusion is along in-plane x, while width and thickness
are well mixed and their exterior faces are sealed in the transport reduction.
Each x strip re-resolves its original isotropic or orthotropic elastic claims at
its actual local moisture. The existing plane-stress/axis lowering and native
DKT assembler rebuild stiffness and translational/rotary inertia. Reference
dry density times (1+M) supplies wet mass on this fixed geometry. The assembled
mass is independently compared with conserved dry mass plus transported water.
Ordinary modal/force/pressure consumers receive the native pencil; no modal
frequency multiplier or instrument preset is introduced. Coupled advance
returns transport and mechanical endpoints together only after both admit.

This plate rung does not model hygroscopic swelling, humidity prestress, sorption
heat, hysteresis, free liquid, stress-assisted
diffusion, changing porosity/geometry, dynamic modal-state remapping, or
moisture-dependent damping. It makes no complete thermodynamic energy or
chemical-entropy claim. Its no-swelling/fixed-geometry approximation must be
admitted by the caller; it is not the full MR46 completion. Drying/wetting
response is memory through the diffusion field, not a sorption scanning law.
No physical wood diffusivity or validation is inferred from authored fixtures.

Focused tests cover a closed-form one-cell update, sealed water conservation,
Fourier-series slab refinement, squared-length/time scaling, data-only
diffusivity changes, local nonlinear source refresh, the shipped FPL sorption
pack, native modal stiffness/inertia scaling, spatial anisotropic material
consumption, complete mechanical-source refusal, budgets, quantity-basis
refusal and cancellation. Native Cargo execution is a separate gate; parsing
and independent numerical oracles do not establish executed native test results.

### Source-backed finite surface moisture exchange

`moisture::surface::DryBasisSurfaceTransfer` admits the canonical
`dry_basis_surface_emission_coefficient` source in m/s. Its convention is
`j_water,in = rho_dry * h_s * (M_equilibrium - M_surface)`, with both moisture
contents in kg water/kg dry solid. A coefficient normalized to vapor
concentration or vapor pressure is not interchangeable with this source.
The reference dry density and inventory remain owned by the original slab.
This surface-emission boundary follows the explicit surface-resistance form
in [USDA FPL's diffusion treatment](https://www.fpl.fs.usda.gov/documnts/pdf1993/simps93a.pdf);
the implementation imports no material coefficients from that paper.

`DryBasisSorption::with_surface_transfer` attaches the pinned coefficient to
the existing isotherm. Both sources accept only the same fixed absolute
temperature and actual relative humidity; additional dynamic surface axes
refuse. `MoistureFace::SourceSorption` exposes this law to prescribed ambient
RH histories. Existing finite-atmosphere, native plate and native solid
entry points consume the same attached law without a separate solver.

The boundary coefficient combines the half-cell diffusion resistance
`dx/(2*D)` and surface resistance `1/h_s` in series. The native backward-Euler
operator remains SPD with nonnegative exchange coefficients. A zero `h_s`
has the sealed diffusion limit; omitting the surface law preserves the
original instantaneous-equilibrium arithmetic. Surface content is recovered
from that same flux and admitted against the original transport source.
Reports distinguish actual surface content from reservoir equilibrium and
retain sorption and transfer source evaluations. No surface water is stored,
clipped or counted twice.

For a finite atmosphere, every conserved-water residual trial resolves both
the equilibrium and transfer coefficient at its endpoint RH. Thus a changing
humidity-dependent surface law changes uptake, vapor depletion, native wet
inertia, stiffness and compatible swelling through the existing owners.
This adds no spatial vapor diffusion, airflow solve, coating storage,
surface-moisture-dependent kinetics, hysteresis or sorption heat. A source's
calibration and the original isothermal dilute bound-water limits remain
required. Broad MR46 stays open.

Focused authored regressions cover independent one-cell uptake/desorption
and surface-flux formulas, sealed and instantaneous limits, a closed-form
finite-vapor balance with humidity-dependent surface transfer, continuation
and deterministic replay, actual native plate response, and late source,
mechanical and cancellation refusals. Native test execution and physical
validation remain distinct from these authored verification cases.

### Finite atmospheric water inventory and native material response

`moisture::atmosphere::FiniteMoistureAtmosphere` replaces the infinite surface
reservoir with the actual water-vapor inventory of the original
`ClosedGasParcel`. The gas is well mixed at prescribed fixed temperature and
volume. At least one noncondensing carrier and an explicitly identified water
species are required; the original Buck liquid-water saturation law retains
its stated temperature domain. The solid keeps its original isothermal
dry-basis finite-volume field and conserved dry skeleton.

`DryBasisSorption` pins an equilibrium surface-content property with semantic
kg-water/kg-dry-solid output. Its two query coordinates are typed absolute
temperature and relative humidity; duplicate aliases, wrong moisture bases,
additional conditions and unsupported source states refuse. Its equilibrium
content and optional finite surface resistance enter the existing
backward-Euler Fick boundary, not the entire interior field. Both ends may
share the same atmosphere or be sealed.

A bounded scalar solve chooses endpoint RH so endpoint vapor plus endpoint
bound water equals their initial total. Every trial uses the original slab
diffusion/source solver with its existing accepted-start diffusivity lag and
work limits. The accepted vapor mass and unchanged carrier masses determine
the new mole composition and ideal-gas pressure at the same T and V. The
original thermochemical gas model then re-resolves caloric, transport and
acoustic properties from the same retained source catalog. Exact vapor
depletion removes the inactive species without a trace floor while preserving
its source catalog for later desorption. An already balanced initial RH,
including sealed diffusion, retains the complete gas state unchanged.

Publication checks actual reconstructed total water, total gas-plus-solid
mass and each noncondensing carrier independently. Source, native diffusion,
bracket, iteration, conservation or cancellation refusal changes neither
accepted owner. Reports retain both sorption source evaluations, actual gas
and solid water changes, signed surface transfers and the gas internal-energy
change. No new mixture-property evaluator, diffusion solver or species store
is introduced.

`MoisturePlate::advance_finite_atmosphere` rebuilds native DKT stiffness and
wet inertia from that same moisture endpoint. The corresponding
`MoistureSolid::advance_finite_atmosphere` resolves the original compatible
swelling, stress and native prestressed tangent. Each wrapper publishes gas,
moisture and mechanical endpoints together only after the original mechanical
source and assembly admissions succeed. A changed source isotherm,
diffusivity or finite atmospheric volume therefore changes the coupled native
response without a material-name switch or a frequency multiplier.

This is isothermal water-mass coupling. A prescribed thermostat maintains T;
bound-water caloric storage and sorption heat are absent, so reported gas
`delta(m*u)` is not the thermostat's total heat or a closed total-energy
account. Gas volume stays fixed and its changed pressure is not automatically
applied as a mechanical surface traction. Condensation, hysteresis, gas
momentum, geometry-dependent reservoir volume, pressure-dependent sorption
and thermal moisture coupling remain outside this slice. It advances MR46
(`frankensim-material-reality-ui2if.7.4`) without closing the broader bead.

Focused synthetic regressions use independent one-cell adsorption/desorption
and continuation formulas, the native ideal-gas pressure and NASA-9 storage
limits, conserved carrier and water inventories, actual gas transport/sound
changes, data-only native plate response, compatible solid deformation and
prestress, exact sealed/equilibrium retention, and late source or mechanical
refusal/cancellation. Rust parsing and these authored tests do not establish
executed native results or empirical sorption accuracy.

### Simultaneous implicit conservative reaction mechanisms

`GasReaction::with_implicit_mechanism(implicit::coupled::ImplicitMechanismSolve)`
solves every canonical reaction extent together with the existing
`fs-solver::NewtonKrylovState`. The residual is
`F_r = xi_r - dt V rate_r(n0 + S xi, U0, V)`. A trial reconstructs one original
`ClosedGasParcel` at conserved mass, volume and total internal energy, then
evaluates every original material rate and detailed-balance reverse closure at
that same actual temperature, pressure and composition. NASA-9 formation energy
is included once by the existing gas chart. No reaction heat source or new
caloric/rate engine is introduced. Existing piston, flow/diffusion and finite
thermal/electrical consumers select the method through their original
`parcel_candidate` calls.

Each channel has an explicit signed extent range keyed by its original reaction
ID. Newton uses dimensionless extents/residuals normalized by the channel's
largest absolute range endpoint. Physical acceptance separately checks every
unscaled molar residual, source support, stoichiometric species balance, mass,
total-U balance and entropy. A converged native solver report alone is
insufficient. Canonical accumulation uses the explicit mechanism's compensated
sum owner. The numerical Jacobian admits feasible centered/one-sided columns
inside the stoichiometric positivity domain and caller ranges; it supplies the
existing native Krylov/preconditioner interfaces and has no analytic tangent or
adjoint authority. Negative-inventory/range trials request native globalization;
missing source data, arithmetic failure, exhausted work or cancellation latch a
hard refusal. Accepted species are never clipped.

Initially absent intermediates can be produced and consumed in the same
implicit interval. A conservative interior *initial guess* is built only on the
reaction support reachable from the original species. Inaccessible reaction
channels remain at zero extent, with their original rates and residuals still
checked at acceptance. The seed is never published or counted as a physical
burn. For this simultaneous method, the existing consumed-fraction control
bounds gross withdrawals divided by original stock plus simultaneous channel
production. Explicit and scalar integration retain their original-stock
denominator. Source domains and finite representability can refuse an interior
seed or derivative stencil; there is no guaranteed solve at every valid state.

The original 32-channel bound remains. Outer attempts are capped at 128, total
residual evaluations (including stencils/rejected trials) at 65536, and native
restart/cycle counts at 32 each. Every gas inversion retains its own work bound.
Cancellation is polled through physical work and native solver callbacks.
Accepted gas and enclosing motion/heat history remain unchanged on any refusal.
The report retains every channel's endpoint residual and the native nonlinear
report. Existing dyadic refinement replays complete coupled trajectories; an
aggregate report retains solved equations in its actual substeps instead of
mislabeling summed extents as one nonlinear root.

This is a first-order homogeneous backward-Euler reaction stage with the existing
reaction-then-motion/heat split, not a validated combustion mechanism, flame or
shock solver, equilibrium acoustic derivative, or nonlinear uniqueness/stability
certificate. The scalar `with_implicit` API still accepts one reaction; it and
the coupled method cannot be bound together. Direct manufactured regressions
cover a stiff A/B/C chain against an independent matrix inverse, absent
intermediate production/consumption, source-driven temperature/pressure/sound,
inactive species, source/budget/cancellation refusals, native piston energy,
replay and actual subcycling. Native Cargo/DSR execution was unavailable in this
environment; source parsing and independent numerical references do not replace
those gates.


## Transported moisture to compatible solid deformation and prestress

`moisture::solid::MoistureSolid` maps each original dry-basis finite-volume cell
to a reference rectangular-prism strip with six native P1 tetrahedra. Tetrahedra
never straddle moisture cells; the well-mixed transverse faces remain sealed.
The source-bound solid independently selects isotropic or oriented orthotropic
elasticity and fixed principal hygroscopic directions. Every accepted moisture
cell re-evaluates the pinned elastic bundle and integrates its three hygroscopic
strain coefficients from the original explicit source length/moisture datum.

The existing solid owner solves compatible quasistatic displacement, recovers
actual stress/support reactions/elastic energy, admits the deformed geometry,
and assembles the original signed reference-configuration initial-stress pencil.
Dry skeleton and transported water own inertia. Reference element density is
wet cell mass/reference strip volume, with current density obtained from wet
element mass/admitted deformed volume; no fresh density curve creates matter.
Native mass and transported inventory must agree under an explicit tolerance.
Fully prescribed bodies retain stress and reactions and have no modal pencil.

`advance` returns proposed transport, geometry and all mechanical projections
together after source/domain, equilibrium, strain, geometry, mass and cancellation
gates. Failure retains the original accepted slab and does not debit reservoir
water. The source-driven outputs expose the native modal/stability pencil;
negative eigenvalues remain instability, never converted to audible frequencies.
Diffusion stays on its original reference chart. This is slow isothermal open-
reservoir hygroscopic mechanics: chemical reservoir work and sorption heat are
not closed by the mechanical energy report. Stress-assisted diffusion, hysteresis,
finite deformation, a deformed-configuration tangent and dynamic mode-state
transfer remain separate models. Analytic fixtures verify transport-to-free
swelling, constrained stress/energy, rotated directions, orthotropic data swaps,
spatial variation, native wet inertia, changed modal pencil and atomic refusal.
All fixture coefficients are authored verification data, not wood validation.

## Spatial magnetic thermal feedback

`circuit_enthalpy::eddy::thermal` extends the same `MagneticThermalNetwork` with `step_eddy_thermal_candidate`. Every thermally mapped triangle, conducting or insulating, has a complete source-pinned `ThermalMagneticMaterial` law and explicit |B| envelope. Original material/caloric context, immutable card and selected claims, whole-volume receiving cells, density, geometry, electrical regions and accepted clocks stay bound. The summed field heat-capacity lower bounds plus each receiving cell caloric capacity must remain positive over a common supported temperature interval. Both accepted endpoint fields must respect their envelopes. Fixed external regions require temperature-independent permeability; changing field entropy cannot disappear into an unmodeled reservoir.

Every original enthalpy/Newton trial starts from the same accepted A/H, resolves local sigma(T), mu(T) and thermal conduction, rebuilds the original planar operator, and calls the original endpoint-effort induction solve. Its native per-triangle sigma E² integral is deposited once through the same heat router as the fixed-mu path. Separately, Q_rev=−P−delta(T S_field) enters that triangle receiving cell. Native magnetic parameter work P and signed reversible heat remain distinct from Ohmic heat. The returned report independently reconstructs delta U_field+delta U_caloric−impressed work−external heat, with U_field=psi+T S; reversible heat is internal, not external supply. Both retained histories publish together, and cancellation/source/budget/nonlinear failures publish neither.

Focused synthetic tests cover both permeability-temperature slope signs, the zero-slope existing limit, independent scalar elimination/root references, repeated accepted-history evolution and replay, spatially unequal permeability/conductivity with local temperatures, summed capacity bounds, missing/duplicate/foreign sources, field bounds, and interrupted/budget-exhausted proposals. This is a fixed-geometry linear reversible affine-permeability rung with first-order thermal endpoint coupling, not saturation, hysteresis, phase evolution, a time-error certificate, unconditional entropy/passivity, or measured magnetic-material qualification.

## Spatial solid cells exchanging heat with one finite source gas

`solid_enthalpy::gas::step_with_isochoric_gas` closes the existing thermal
backward-Euler cell balances and the existing fixed-volume NASA-9 gas energy
balance in one bounded Newton--Krylov solve. Every exposed patch resolves the
original gas and solid conductivity at the common trial endpoint. Its signed
heat is removed from its own cell and added to the single gas inventory;
opposite patch transfers remain distinct even when their net gas heat is zero.
Existing internal links and authored external heat remain on the original
thermal owner. No mean solid temperature or substitute heat capacity is used.

The source-bound structural field is resolved at every trial and at admission.
Conductivity selection is pinned at the original stress-free reference, and
the accepted start must still support it. Local cell, gas and combined caloric
energy allowances are independent of the nonlinear stopping rule. Contact
IDs, aggregate exposed area, cell/patch/trial budgets, source domains and
cancellation are admitted before any endpoint can be returned. The immutable
candidate owns thermal, gas, material-field and transfer reports together;
callers perform any required equilibrium/valve admission before publishing it.

This slice assumes a well-mixed gas at fixed volume and composition, fixed
thermal geometry and authored internal conductances. It is sequential
quasistatic thermoelastic loading: reciprocal elastic heat, pressure work,
gas flow, moving surfaces and phase change are not supplied by this coupling.
Numerical Jacobian probes must remain inside the declared material domains.
Independent three-reservoir regressions cover local heat, source conductivity
and gas-volume changes reaching displacement, zero-net spatial exchange,
source continuity, cancellation and transactional refusal.

### Computed DG radiation and source-bound exterior pressure

DgPortLoad::new_with_computed_baffled derives the single circular outlet from
original DG area and an explicit shared GasRecord. ComputedBaffledInput now
exposes its existing aperture/gas bake independently of waveguide topology;
its original BEM samples, rational fit, full mechanical Nyquist admission and
source/held-out/continuous/sampled reports remain unchanged. The retained
ComputedBaffledRecord owns source geometry, gas receipts and numerical results,
not another copy of the accepted physical terminal history.

The shared dg_tube/dg_bore reader admits OUTLET baffled followed, after all
geometry and gas-station rows, by dg_radiation computed-baffled with the original
fourteen controls and dg_exterior with a complete gas record. Static pressure
must match the quiescent bore; mixture and temperature may differ. Cold source
preparation rebuilds the physical load once. The same source serves moving
apertures and original jet-driven flues through the shared DG endpoint solve.

FlueDgPerformance::with_baffled_outlet and AperturePerformance's existing
TubeBaffled selector consume each accepted outward DG terminal flow once in
the original same-clock differentiating BaffledPressure owner. Their opening
and exterior come from that computed source; observing adds no second load.
The default original inlet observation remains available with the same reaction
active. Geometry, source/wall frequency admissions, receiver limits and clock
remain checked; a conflicting explicit exterior refuses. A failed pressure
callback poisons its finite wrapper rather than claiming complete acceptance.

Both authored formats reach the existing pressure renderer and explicit CLI
decimation. Metadata retains the computed source even for inlet observation.
An active ensemble microphone uses its actual bound exterior gas, checked
against the numerical receiver medium and existing common-listener policy.
Focused tests cover gas/geometry source records, original thermochemical
receipts, accepted-flow observer equality, physical-state preservation,
callback partitioning, receiver refusals and pressure-ceiling poisoning.

This is one uniform circular opening in an infinite rigid baffle, without
mouth mutual loading, unflanged bell diffraction, body scattering or room
physics. The original flue jet remains an exploration source. Wall-law,
spatial, rational and nonlinear-time convergence stay separate; sampled
passivity remains GridOnly. No measured flute, clarinet or oboe qualification
is asserted. Source tests and illustrative fit declarations require native
execution before numerical or acoustic accuracy claims.

## Native finite-gas heating of a spatial material valve

The `heated_valve_sound` example accepts an explicit finite gas volume, layer
distance and step count alongside the existing species/NASA9/transport packs.
The original 48-cell solid receives gas heat on actual exterior triangles,
with an insulated root and conservative shared-face internal links. The native
slice pins the original conductivity source and resolves it at each cell's
actual nonlinear endpoint within the declared fully solid 280–400 K interval.
Internal half-cell links and gas-facing layers consume the same material
snapshot; constant claims, temperature curves and supported fixed-context
response surfaces use the original material-network owner. Its centroid-normal
two-point network remains a coarse transport model, not a continuum-consistent
tetrahedral heat discretization.

The optional manufactured `--demo-conductivity K_AT_280 K_AT_300 K_AT_400`
source permits conductivity comparisons with identical initial coefficients.
A regression follows the changed temperature dependence through spatial heat,
quasistatic deformation, reduced valve frequency and actual pressure, with
combined storage closure and deterministic replay checks.

Each accepted temperature field uses the same source stiffness and expansion
owner as the final equilibrium, prestress, aperture and valve reduction.
Initial and final gas snapshots feed the existing common-medium acoustic
realizer. Pressure CSV and stereo WAV preserve physical pressure scaling; an
additional thermal CSV reports spatial temperature extrema, gas T/p and
independently measured solid/gas storage and cumulative heat. Both complete
frozen-state acoustic experiments precede output creation. The gas volume is
thermal storage, not a modeled acoustic enclosure; pressure work, live modal
state transfer and reciprocal thermoelastic heat remain outside this slice.

## Material valves on the original DG bore

`bernoulli_aperture::dg::LoadedApertureDg` binds the existing
`LoadedAperture` single, independent-paired or material-modal owner to the same
DG mean-port relation used by the scalar valve and flue source. It requires cold
ExactZoh mechanics, matching physical clock, inlet density and mean impedance;
retained paired/modal gas states must match the complete inlet material state.
No blade, mode, contact coordinate or constitutive memory is collapsed to the
diagnostic opening. Both wrappers share one acoustic/work acceptance path.
Additional source-mesh pressure work appears in the frame and is already
credited once in mechanical work. Rational terminal exchange remains signed
port work, not an invented heat or storage certificate. Failed trials leave
mechanical, wall and terminal histories unchanged.

`dg::loaded::performance::LoadedDgPerformance` adds a finite canonical
blowing-pressure phrase and the existing source-mesh force/pressure compilers.
All auxiliary controls share one remaining storage/work allowance. Each accepted
mechanical sample advances its source once and publishes its controls; a failed
callback poisons further output while retaining the accepted prefix. Inlet
pressure and an explicitly selected computed baffled outlet are distinct
observations. `DgPortLoad::baffled_receiver` is the shared scalar/flue/loaded
admission owner for original outlet geometry, exterior gas, physical clock and
radiation/wall bands. It feeds the existing retarded pressure observer without
rebaking or adding a reaction.

Regressions exercise original scalar arithmetic, independent material blades,
all retained modal/Maxwell coordinates, prescribed traction work, rollback and
scheduled playback. These are coupling invariants, not measured
clarinet/oboe fidelity or convergence certificates. Native test execution and
instrument-specific time/space/material validation remain separate requirements.

## Authored material-modal DG sources

`frankensim-material-modal-dg-valve-v1` makes the loaded DG composition
reachable through the existing material-card/mesh readers and `music_render
wind`. One or two original regional plates retain all requested signed modes,
selected tensor viscosity or Maxwell histories, wetted geometry, spatial
pressure, contact closure and source-mesh force controls. The material source
binds the actual DG inlet gas and mean impedance. Shared DG input owns local
gas fields, wall memory and any computed exterior reaction; material inputs
do not replace those owners with density-only tuning.

The source declares equal mechanical/output clocks and ExactZoh; downstream
PCM decimation is explicit. A named inlet trace remains internal pressure.
Baffled output requires the original computed terminal and uses its exterior
gas, outlet geometry and accepted outward flow. The supplied paired-DG
fixture and regression sources are manufactured integration data, not an
instrument qualification. Clock/structure, missing material claims, absent
pressure ports and unbound exterior observations refuse without substitution.

## Authored distributed interior for material strings

An optional `interior-fe-v1` trailer on v3/v4 body profiles supplies the
original tetrahedral volume, pressure FE order, complete frequency window and
explicit numerical allowances to the existing Helmholtz and
`with_material_shell_cavity` owners. Gas properties come only from the
performance's coherent material gas or its explicit `cavity_gas` source. The declared volume must conform to the
original shell's complete cavity surface; missing walls and open apertures
refuse rather than acquire invented closing faces. A separate lumped cavity
cannot be supplied simultaneously.

Preparation installs this reciprocal air load before subsequent cylinder-fluid,
torsional, nonlinear, material-memory, radiation and bow/performance attachment.
Whole-network admission reserves every possible pressure inertia and compression
connection before material reads. The exact uniform pressure mode contributes
compression but no fictitious air inertia. Body radiation keeps its original
structural basis while the same accepted solve carries strings, shell and air.

The supplied closed-shell-cavity source contains its exact shared surface and
P2 tetrahedral geometry. Its compact body observation is an integration
diagnostic; the fixture does not claim open violin/cello/bass geometry,
aperture radiation, interior viscothermal loss or measured-instrument fidelity.
Regression source covers physical motion/PCM propagation, pressure/composition
and solid-material substitutions, retained fluid augmentation, zero-inertia
uniform compression and strict resource/geometry admission. Fixture topology
was checked independently; native Rust execution remains required.

`with_material_shell_cavity_in_gas` binds an independent frozen interior to the
same FE, compression, compact-neck and optional thermal-diffusion owners. The
bank's exterior gas remains authoritative for string drag, exterior loading and
pressure observation. The interior FE density/sound speed must exactly match its
declared interior GasState. Both regions must share the static reference
pressure to roundoff; mean pressure preload, heat exchange and species diffusion
are not supplied. Native `cavity_gas` follows `gas` before `body`, retains the
original GasRecord and refuses when no actual cavity consumes it. Omission
preserves the original equal-medium trajectory. Distinct lumped media require
an exterior-owned mouth and its native computed skin load; one-medium compact
outer-end models refuse substitution instead of silently using interior density
for exterior inertia. No new gas law or time integrator is introduced.

## Material-wind ensembles and complete numerical observation latency

The existing ensemble `--valve` input selects plate-valve, material-modal and
material-modal-DG schemas without replacing their physical renderers. Material
modes, constitutive histories, controls, reaction loads and microphone histories
stay source-owned. Pressure-producer clocks govern outer conversion; actual
mechanical clocks/counts are recorded separately. Listener compatibility uses
the microphone's bound complete exterior gas. Internal DG inlet traces are
explicitly identified and do not infer a listener gas.

`PressureEnsemble::new_with_source_delays` accepts the numerical receiver
latency already present in each producer, in final-output samples, and aligns
that plus the existing outer decimator delay using the original pure-delay
owner. Every total must be finite, nonnegative, integral and at most 4096
samples. Fractional totals refuse. Physical propagation remains within each
receiver; no source state, source event or sound-travel time is resampled to
perform alignment. The original `new` constructor retains zero-source-delay
behavior. Focused regression source compares independently filtered modal
trajectories across callback partitions and covers exact zero-delay parity and
invalid latency admission; native execution is still required.

## Material-driven sliding heat allocation (MR20)

`tribothermal::SlidingHeatPartition` makes the heat model explicit:
`Prescribed` retains supplied fractions; `SemiInfiniteEffusivity` uses each
surface's original material and accepted caloric state. Conductivity is queried
at each body's own endpoint temperature, preserving its external coordinates
and pinning the initially selected claim. Density and differential heat
capacity are read from the same original enthalpy chart used by the thermal
solver. There is no second Cp law, replacement material, constant heat split
after a swap, or separate heat integrator.

The selected split participates in every nonlinear endpoint trial, including
the finite-body dynamics and field-driven motion paths. Friction work is still
partitioned and deposited exactly once by the existing thermal owner. The
report retains the conductivity bundles, original phase states, differential
Cp, density, diffusivity, heat fractions and current contact age. Failure to
resolve either material or to fit the model envelope returns no new thermal,
slip, mechanical or work state.

This is an explicit local half-space estimate for incremental friction heat:
homogeneous isotropic media, equal exposure times, negligible thermal
resistance at the generation plane. It does not predict flash temperatures or
moving-source Peclet effects. Background heat exchange remains the separately
declared thermal network. A caller-specified diffusion-depth envelope compares
`2 sqrt(alpha_max age)` with the distance to competing boundaries (including
footprint edges); age and the maximum sampled diffusivity survive accepted
steps and cloned restarts. The envelope is an approximation guard, not a
finite-domain error certificate or proof of property behavior between sampled
states. Contact resistance, third-body layers, melting and inadequate spatial
resolution require other models.

Focused G1/G3 regressions in `tests/tribothermal/partition.rs` exercise changes
to conductivity and Cp through actual finite heat stores, independently solve
a temperature-dependent conductivity feedback case, and verify cumulative
exposure refusal, ordered-source identity and unchanged accepted state.

## Source-dependent internal transport in the spatial solid/finite-gas solve

`solid_enthalpy::gas::step_with_isochoric_gas_transport` composes the original
prepared material network with the same backward-Euler solid/gas residual as
`step_with_isochoric_gas`. A resistor-free network pins each original cell's
conductivity source once, refreshes it at the actual endpoint enthalpy/phase,
and lowers existing two-layer links through their original conductance owner.
No dummy circuit, alternate thermal state or additional heat integrator is
introduced. Networks with resistor work refuse this thermal-only entry point.

The structural and transport bindings must share their original material,
borrowed caloric chart and complete physical condition point. Gas patches
consume the network's existing material snapshot through the original layer
adapter, including the exact conductivity key and reference source pin. The
accepted report retains those cell/link source snapshots alongside local,
gas and combined energy accounts. Failed source refresh, nonlinear admission,
geometry or cancellation leaves both inventories intact. The original entry
point retains authored internal coefficients; geometry and all original
quasistatic, single-phase and fixed-volume limitations remain unchanged.

Three nonlinear-transport regressions exercise different material cards,
independent three-reservoir equations and compatible nonuniform displacement,
subsequent source refresh, conductivity changes reaching gas/motion, and
wrong-context/chart or failed-trial refusal.


## DG flue mouth observation and compact reaction

`FlueDgPerformance::with_mouth_radiation` observes the source's actual flue
width and cut-up through the existing rectangular Rayleigh owner. Placement
and surface resolution are explicit; the already computed terminal supplies
the complete exterior gas and common microphone. Each accepted negative DG
inlet flow enters that receiver once, then its physical pressure is added to
the terminal pressure. Each aperture retains its own propagation and numerical
filter support. Geometry, equal static pressure, quadrature, source band and
cold attachment remain original-owner admission checks.

Optional `with_mouth_reaction` attaches the existing `RectangularMouthLoad`
through `FlueVoice::attach_mouth_radiation` before playback, replacing the legacy
mouth load and feeding its retained passive state into the current jet/DG
solve. The original source equations, controls, pressure ceiling and callback
failure behavior remain in effect. The file grammar adds optional
`mouth_radiation X Y WIDTH_CELLS CUT_UP_CELLS` and
`mouth_reaction none|compact-baffled` after `observation baffled ...`; a physical
reaction requires zero legacy `mouth_loss_rel`.

This is a disjoint rectangular mouth and circular terminal in an infinite
rigid baffle, with independent self reactions. Mutual mouth/terminal loading,
free-field body scattering and a resolved compressible jet are not provided.
Focused tests compare the pressure against the original Rayleigh owners fed
with accepted signed flows, source-state change under mouth reaction, block
partition equality and physical admission. Native execution and measured flute
qualification remain pending.


## Complete DG terminal groups

`DgPortLoad::new_with_terminal_ports` binds every non-inlet port exactly once
in explicit MIMO channel order. The complete mean relation is
`(Zdd + L) u = Hd - Hl + Zd0 Qin`, with outward flow into the original
`PreparedMimoLoad`. Its full cross entries and signed exchanges are retained.
The one-channel branch preserves the original divide; larger constant systems
use the existing normalized prepared LU with infinity-norm condition <=1e10.
No terminal pressure, work or accepted memory is published until every original
DG/load candidate and the nonlinear source accept. An abandoned trial advances
nothing. The coupling is exact for each owner's held-flow image, not the joint
continuous exponential; coupling-clock refinement remains separate.

`terminal_ports` and `last_terminal_frames` retain declared matrix order.
`last_terminal_frame` remains only the actual axial-outlet convenience, including
None for side-only groups. Aperture and flue frames carry every terminal;
aperture balance credits their signed work once, with absolute per-port scales
for numerical cancellation checks. Tests cover independent 2x2 elimination,
raw original DG/MimoZoh state and work, both kernels, wall memory, rollback,
side-only ownership and bitwise legacy scalar arithmetic. These regressions
are authored; native execution remains pending.


## Authored DG openings with shared exterior reaction and observation

`DgPortLoad::new_with_computed_baffled_ports` binds explicit ordered port/circle
pairs to one original computed BEM/VF source and one complete terminal group.
Every non-inlet port is required. Axial geometry retains the original outlet
area; side geometry is independent and uses the actual face's gas for common
static-pressure admission. All mutual source entries, local bore/wall physics,
full mechanical Nyquist fit and original passivity/error refusals are retained.
No duplicate propagation or load owner is introduced. A microphone consumes the
same original geometry, exterior gas and accepted ordered terminal-flow vector.
Scalar, paired/material-modal and DG flue performances now pass every accepted
flow exactly once, with buffers allocated at construction.

The shared `dg_tube`/`dg_bore` input admits optional `dg_side_openings COUNT`
before `dg_radiation`, followed by `dg_side_opening PORT_ID DOWNSTREAM_ELEMENT
X_M Y_M RADIUS_M`. IDs exceed 1; faces and IDs are unique. The computed source's
existing eight-opening cap permits seven side openings plus the axial outlet.
The original element face is explicit: no station snapping occurs. These are
zero-volume flush plane-wave junctions in one infinite rigid baffle, without
finite chimneys, moving pads or free-field instrument-body scattering. Flue
mouth reaction remains its separate compact self load; terminal mutual terms
do not imply mouth-to-terminal mutual reaction.

`examples/flue-dg-side-openings.performance` exercises the complete authored
selection. Focused source tests retain original face/material/area identities,
reverse matrix order independently of DG order, and compare the common
microphone with original individual receivers driven by actual accepted flows.
The fit, native tests and example audio remain unexecuted in this environment;
no instrument-accuracy or measured-audio qualification is claimed.


## Whole-interval reactive material heat and piston refinement

`cavity_phs::piston::reactive::thermal::subcycling` composes the existing conservative reaction/thermal/motion candidate on real subintervals. `GasPiston::reactive_thermal_subcycled_candidate` and `reactive_thermal_layers_subcycled_candidate`, with the corresponding `MaterialGasPiston` forwarding methods, feed each complete accepted gas, finite-body enthalpy and wall-motion endpoint into the next chemical stage. Wall conductivity, mixture transport, phase/enthalpy and material rate sources therefore affect chemical progress within the requested macro interval. The original rigid material binding remains the sole piston inertia authority.

`ReactiveThermalSubcycling` compares complete N and 2N trajectories from immutable initial owners, using independent absolute scales plus a relative tolerance for each species amount, each canonical reaction-channel extent, gas temperature/pressure/volume, wall velocity and body specific enthalpy. The finer actual trajectory is accepted without extrapolation. Dyadic meshes are bounded at 256 substeps and all attempted coupled steps at 511; every nested chemistry, heat and motion owner retains its original separate work limits. Donor/entropy/iteration and resolved-clock refusals may request refinement; missing material support, invalid data and cancellation remain hard failures. The callback is repeatable and cannot debit live external state. This is controlled refinement of the existing first-order split, not a monolithic implicit chemical/thermal closure or detailed combustion/engine validation.

Integrated gas and body source heats are apportioned at constant power with cumulative fractions and a final remainder; exterior pressure and applied force remain constant. The returned `UniformEnthalpyStep` covers the complete interval, including the sum of boundary heat and the original source-bearing reports of every accepted substep. Independent macro checks enforce the original species, isolated chemical energy, finite-body energy, total energy and momentum tolerances. Momentum uses actual endpoint MΔv against separately accumulated midpoint gas-pressure, exterior-pressure, applied-force and support impulses. Chemical formation energy is never credited again as heat. One macro port clock/account replaces tentative substep accounts and retains actual accepted work, gas endpoint entropy, summed heat and the full heat entropy integral Σ(Q/T); heat reversals do not create a fictitious single boundary temperature.

Focused regressions in `cavity_phs::piston::reactive::tests::thermal` cover material-only conductivity swaps changing intrainterval reaction progress, piston compression feeding source-backed rates, independent reversible-isomer convergence under tighter tolerances, whole-interval source/heat/impulse accounting, and bounded-work/late-source/late-cancellation rollback. Synthetic source cards establish numerical and integration behavior; they are not empirical combustion calibration.


## Explicit compact mouths on the material shell cavity

`ShellCavityOpening` binds a connected planar disk of original interior FE
boundary faces, internal neck length and nonnegative internal resistance.
`with_material_shell_cavity_openings` integrates its actual area and every
retained pressure-shape average through `HelmholtzFe::boundary_load`, then
uses the original `CavityCoupling` neck inertia and collective volume springs.
Each face has exactly one solid or declared fluid owner. Existing solid caps,
undeclared holes, duplicate/disconnected patches and moving retained rims
refuse. The maximum patch extent and internal length must remain within the
original compact acoustic chart across the declared band. No exterior end mass
is included in internal neck data. The sealed constructor keeps its old path.

A joint exterior source map contains all structural modes followed by actual
neck coordinates, excluding interior standing-wave inertias. Existing shell
BEM receives solid panels and explicit fluid panels with normal trace W/area;
its complete passive Model acts through the original exact mechanical network.
The existing finite-distance observer instead labels physical aperture volume
flow Q=W*v_neck and uses the inverse-scaled 1/area trace. Body-mouth and
mouth-mouth terms, material memory, cold rebuild mappings, signed work and
accepted physical state remain with their original owners. New readouts expose
accepted neck pressure, flow and the exterior reaction without another solve.

The unloaded constructor explicitly means zero exterior pressure. The atomic
`with_material_shell_cavity_and_bem` and the existing later BEM attachment both
require the unchanged exterior admissions. In particular, dense collocation
does not guarantee the modal fitter's 1e-12 reciprocity bound for independent
localized body/mouth columns. Such data refuse; the matrix is never repaired or
its tolerance relaxed. Source-derived general open-body radiation is therefore
not qualified by this integration. Native regressions and the manufactured
joint-BEM example remain pending and may expose this existing numerical limit.
The direct FE/neck/passive-network tests cover original normalization, pressure
reaction and work, atmosphere substitution, geometry refusal and retained source
maps. This is a fixed compact-mouth model, not distributed f-hole flow, interior
thermoviscous loss, moving aperture rims or measured violin/cello/bass accuracy.


The material-string v3/v4 body profile exposes these same cavity mouths after
its `interior-fe-v1` tetrahedra with `air-openings COUNT`, then
`air-opening INTERNAL_LENGTH_M INTERNAL_RESISTANCE_PA_S_M3 FACE_COUNT` and
that many `air-opening-face A B C` records. Original FE indices and geometry
are retained. Each opening reserves one fluid mode and one possible internal
damper connection before material-pack reads. Opening-bearing profiles require
BEM and a covering band at preflight. The cavity binds at the original early
assembly point; mandatory joint BEM binds later, preserving existing nonlinear
attachments. No incomplete performance is returned if exterior preparation
refuses. Sealed source grammar is unchanged.

The open-shell-cavity profile/performance and focused actual-flow microphone
regression are unexecuted integration inputs. Their 80-panel source may refuse
at the unchanged dense-BEM reciprocity gate described above. No native test,
render, material calibration or instrument-realism parity is claimed.

## Thermally live layered supports in reciprocal field motion

`circuit_enthalpy::motion::support` connects the original
`MaterialFieldMotion` actuator to the shared
`ConcentricThermoelasticLaw`. `AxialSupport::from_thermoelastic_law`
retains the original source-bound section, reference gap/length and consistent
inertia. `MaterialNetwork::field_motion_support_candidate` binds each layer to
one whole finite caloric cell with the same card, state context, reference
volume and mass. Independent layer temperatures update expansion, stiffness,
restoring force and continuing charge/flux, gap and momentum.

Every trial starts from the same accepted histories. At unchanged initial gap,
the member is rebound to its proposed layer temperatures and the original
field-motion descriptor advances once. Each layer receives
`Q_rev = -delta(Psi_at_initial_gap) - delta(T*s)`; its source-backed Joule heat
is deposited independently, once. The native mechanical potential remains
Helmholtz energy. A second physical account closes the change in caloric,
support `Psi+T*s`, field and kinetic energy against electrical/mechanical work
and authored heat. Reversible transfer is removed from the reported external
heat. The existing thermal network owns conduction, reservoirs, events and
bounded nonlinear work; no new field, motion or thermal solver is introduced.

Scalar reference density/expansion and globally affine modulus use the shared
material admission. A conservative bound requires positive combined caloric
and thermoelastic heat capacity over the caloric interval and declared elastic
strain envelope. Original cards, selected claims, thermal references, geometry
and inertia survive rebinding. Ordinary frozen-support stepping refuses a
thermally live support. Enclosed electrical storage coefficients must remain
source-constant because this lane does not supply their field entropy.
Unsupported sources, inventory mismatches, cancellation and failed balances
publish no partial physical history.

This is the existing small-strain axial reduction with fixed reference geometry
and electrical paths. It does not add radial work, interlayer slip, changing
mass, strain-dependent resistance, distributed Lorentz forces, phase change,
magnetic hysteresis or heating of the moving-gap medium. The endpoint thermal
coupling is generally first order; nonlinear tolerances are not time-error
bounds. Bulk electrical heat routes are supported here; source-bound electrical
contact heating remains on the existing separate contact-aware motion path.

Focused tests in `tests/circuit_enthalpy_motion/support.rs` cover unequal layer
heating, a source-only expansion swap affecting actual gap/velocity, independently
reconstructed internal energy, original inertia/flux, positive-capacity refusal,
cancellation and exact retry. The changed Rust files passed syntax parsing;
native compilation and these tests were not executed in the connector-only
session because DSR, RCH, Cargo and rustc were unavailable.

## Continuing layered-string heat, electrical drive and sound

`acoustic_realize::concentric::ThermalConcentricStringRuntime` connects the
shared source-bound thermoelastic section to continuing native string modes,
gas resistance, pressure history and finite thermal cells. Each perfectly
bonded annulus retains its own material, temperature, expansion, modulus,
caloric law and Kelvin–Voigt bending viscosity. Reference cards, complete
property requirements, source selections and original annular inventories
survive every step; there is no equivalent homogeneous material card.

The existing exact modal propagator advances once at the accepted temperature.
The existing backward-Euler network and Newton–Krylov owner then solve layer
temperatures at frozen new modal displacement and momentum. The original
source section rebuilds tension, bending stiffness, damping and pressure
coefficients without resetting physical motion or observer acceleration history.
Each layer receives its actual share of bending dissipation, proportional to
its own `eta_i*I_i`, and reversible heat
`-delta(Psi_axial_i)-delta(Psi_modal_i at fixed new q)-delta(T_i*s_i)`.
The entropy term includes the change in physical modal displacement, not just
temperature. Independent combined-energy acceptance includes full mechanical
internal energy, caloric storage, supplied force work and the native air-loss
account. Local and global heat residuals must also pass before publication.

`step_with_ohmic_heating` reuses the original parallel-conductor model on the
SAME already-resolved trial section. Current or voltage determines each layer's
current and Joule heat from its actual resistivity/conductivity. Total current
is shared by conductance; it is not applied independently to every layer.
The payload separates electrical supply and per-layer Joule energy from other
external heat. Joule heat, physical bending heat and reversible exchange each
enter the original thermal cells once. Independent sources and thermal links
continue through the existing callback and network; no second electrical,
mechanical or thermal integrator is introduced.

This is a fixed-section, fixed-density, small-strain bonded-member model with
scalar expansion, globally affine modulus, solid caloric charts and sourced
Kelvin–Voigt bending. A per-layer capacity check includes axial AND modal
thermoelastic curvature over every original caloric segment. Slack, source
exits, excessive strain/temperature increments, changed inventory, cancellation
and solver/balance failures publish neither thermal nor acoustic history.
The thermal/mechanical split is first order and needs timestep/modal refinement.
The thermal-only and Ohmic methods above retain their supplied fixed gas state;
the finite-gas transaction below evolves it explicitly. Microphone pressure is
the existing diagnostic radiation approximation. This does not thermally evolve a complete
bow/string/body network, resolve helical winding slip, radial gradients inside
an annulus, skin effects, changing geometry or phase transitions.

Four focused tests cover unequal temperatures with continuing motion, independent
internal-energy accounting and layer-loss ownership; source-only expansion and
resistivity changes affecting tension and pressure; implicit `R(T)`/current/
Joule energy; and refusal with exact retry. The changed Rust files passed syntax
parsing and the thermal identities were checked by independent numerical
reductions. Native compilation/tests and audio rendering/listening were not
performed in this connector-only session. No measured instrument-fidelity claim
is made by these synthetic constitutive tests.


## Finite DG side branches

`DgPortLoad::new_with_terminal_branches` joins original two-port `DgDuct` branches to distinct existing bore side ports. Each branch retains its actual geometry, material gas, viscothermal memory, original propagator and state. Junction mean pressure is common and held volume flows are opposite. The constant junction matrix is eliminated through the existing prepared LU; the resulting complete external mean-port response enters the unchanged terminal MIMO coupling with every cross term retained.

External IDs remain the authored bore-side IDs, but terminal pressure, work and microphone flow belong to the actual branch outlets. `DgBranchFrame` separately exposes junction flow, both independently integrated pressure/work observations and the original branch step. Compressible branch storage permits different junction and outlet flows; neither is substituted for the other.

Bore, branch and rational histories publish together only after every original candidate and signed junction/terminal work check accepts. Main `state()` retains its original meaning; `branches()` exposes the independent states. `stored_energy_j`, cumulative loss and coupled-source ledgers sum every original DG owner once. Rational terminal exchange is reported separately without fabricating a storage function.

Cold admission requires complete terminal ownership, two-port branches without nested side ports, common equilibrium junction pressure, a common clock, resolved native wall bands and a combined acoustic-plus-memory DOF allowance before exponentials. Different coherent gas compositions at the junction are permitted. An empty branch list preserves existing scalar and multiport arithmetic. This is a finite one-dimensional duct attached at a zero-volume point junction: it does not add three-dimensional entrance/jet flow, moving pads, or a continuous-time joint exponential. Coupling-clock, spatial and wall-fit refinement remain separate qualifications.


## Finite material-gas chimneys in authored DG instruments

The shared DG source grammar now allows `dg_side_chimney PORT_ID
DOWNSTREAM_ELEMENT X_M Y_M RADIUS_M LENGTH_M ELEMENTS ORDER` among the
existing `dg_side_openings` rows, followed by `dg_chimney_gas` and one
complete shared material gas record. The original uniform-cylinder DG and
Bessel wall owners retain each chimney's physical propagation, storage,
viscothermal memory and gas dependence. Parent wall poles, band and tolerance
apply to each branch; the common state allowance covers the main bore and all
chimneys before any propagator or radiation fit is constructed. Full
thermochemical receipts remain in the existing source record collection.

`new_with_computed_baffled_branches` binds the actual branch outlet area and
gas to the original complete radiation matrix. The source retains every
terminal in its authored order, with accepted OUTER flows feeding the existing
microphone once. An internal junction is not a radiating aperture. Bore,
chimney and exterior compositions may differ at common static pressure; a
mean pressure jump refuses. The microphone band is limited by every retained
wall owner. Empty branch lists preserve the existing flush-opening route.
The same source input serves scalar and material-modal apertures and DG flues.

The new flue and paired material-modal chimney performance inputs demonstrate
data-only selection. Focused regressions check material receipts, original
wall-operator identity, complete memory limits, actual outlet area, equilibrium
admission and the original microphone driven by accepted branch outlet flows.
Native Rust tests and audio rendering remain unexecuted here. This is a
fixed plane-wave junction/finite cylindrical branch reduction, without moving
pads, resolved three-dimensional junction corrections, mean flow, or empirical
instrument qualification.


## One physical Rayleigh source for rectangular mouths and circular terminals

The computed baffled source can now prepare a complete mixed opening image
through `ComputedBaffledRecord::with_driven_mouth`. Channel zero is the actual
rectangular mouth; following channels retain the original circular-terminal
order. Shared receiver surface quadratures provide normal velocity per SI
volume flow to the existing power-consistent Rayleigh owner. Full signed
mutual entries pass unchanged through the original causal fitting ladder,
withheld errors, continuous passivity and interval-mean sampled passivity.
The bound full gas record, mechanical clock, complete band, opening geometry
and source quadrature accompany the actual replacement fit and identities.

Original rectangular and circular static-inertia references independently
refuse insufficient source quadrature without replacing computed entries.
The rectangular reference formula is now shared from the BEM owner.
Complete source-sweep work and coexisting source/fit memory are admitted
before surface allocation. Circular-only inputs retain their original path.
Focused tests cover physical flow normalization, compact mutual phase,
signed interference, gas changes and source refusal; native execution remains
pending. This is a declared infinite-baffle surface reduction, without resolved
jet flow, free-field instrument scattering or empirical acoustic qualification.


## Shared reciprocal material-shell exterior and observer

Material-string `with_material_shell_bem` and
`prepare_material_shell_bem` now select the existing symmetric P1/P0 reciprocal
boundary owner through the shared complete modal-loading adapter. Existing
`ShellBemOptions` fields are unchanged. The corresponding
`*_with_angular_refinement` methods expose the original 0..=2 angular refinement;
the existing methods select level zero. Every original signed structural mode
and explicitly admitted compact mouth volume-flow column enters the same
work-conjugate matrix. The strict existing reciprocity/passivity assessment and
full causal fit precede attachment to the original complete mechanical model.
No computed collocation matrix is repaired or silently substituted.

`ShellBemReport::reciprocal` records actual original unknown/cache/work and
backward-residual/flux diagnostics. `maximum_sweep_terms` admits the complete
original per-frequency source work times the declared grid, including singular
moment preparation, angular Gram assembly, dense factorization and all source
columns. The retained legacy condition-number lower-bound field is the
universal lower bound one; it is not a condition probe or enclosure.

The material-shell observer reconstructs the same original boundary, gas,
complete source order and load angular refinement, then observes its actual
P1/P0 solution at the declared receiver. It consumes each accepted body
generalized velocity and compact mouth's actual SI volume flow exactly once on
the existing mechanical clock. The default receiver uses the original adaptive
near-field quadrature defaults. `with_material_shell_observer_with_quadrature`
exposes that owner's independent tolerance, depth and per-frequency kernel
allowance. The full boundary sweep allowance and the receiver kernel allowance
are separate units; neither replaces the unchanged FIR tail, wrap and DTFT
admission or replay-work cap.

The authored material-string format exposes optional
`shell_angular_refinement LEVEL` immediately after `body-radiation bem`, and
`shell_quadrature REL_TOL MAX_DEPTH MAX_KERNEL_EVALUATIONS` immediately after
`shell_observer`. Missing records retain level zero and the original receiver
defaults; malformed or incompatible records refuse before material-pack I/O.
The manufactured 80-panel/42-vertex, two-source, 65-frequency level-zero example
requires 10,535,395,000 conservative source-work units and declares 12 billion.
Its observer separately declares at most 80 billion boundary-work units under
the unchanged 100-billion hard ceiling, plus the original two-million
adaptive-kernel allowance per frequency. These declarations change numerical
work capacity, not physical materials or fit/FIR tolerances.

Original geometry, resonance-exclusion, angular-resolution, compact-mouth,
gas-coherence, complete-source, cold-history and exact-state gates remain
authoritative. This joins the previously admitted compact uniform opening
velocity to its exterior reaction and receiver; it does not resolve an f-hole
flow field or interior viscothermal losses. Native compilation and the new
integration regressions were unavailable in the implementation environment;
the source-defined gates remain required before qualification of these
manufactured examples.


### Joint DG mouth and terminal exterior

The explicit DG flue mutual-mouth route replaces independent mouth and terminal loads with one original rational MIMO owner. Its rectangular mouth is channel zero with outward flow exactly minus accepted inlet flow; all circular terminal channels retain their physical order, including finite branch outlets. The same-step source effort is inlet pressure minus exterior mouth pressure. The existing Schur solve retains every mutual term, and bore, branch and original rational candidates publish together only after their independent pressure/work checks.

Source preparation reuses the original receiver surface quadratures and power-consistent Rayleigh matrix in one infinite rigid baffle, with the complete retained exterior GasRecord. It keeps the common full-mechanical-band fit, withheld errors, source causality assumptions and continuous/sampled passivity admissions. Sampled admission remains GridOnly; finite surface/frequency checks are not continuum or rational-storage certificates. Signed joint mouth work is distinct from heat. This does not resolve the jet or free-field scattering by an instrument body. Cold source reuse preserves the joint load and does not restore legacy mouth inertia.

The shared DG flute file grammar selects this path with
`mouth_reaction mutual-baffled` after its explicit `mouth_radiation`
rectangle. Existing `none` and `compact-baffled` selections remain distinct.
The same joint load composes with the original finite DG chimney states;
`flue-dg-mutual-chimney.performance` supplies the combined data-only selection.
The receiver retains the source's actual surface quadratures on cold reuse.
Generic aperture adapters refuse this pressure-jump load at construction,
because their existing source law requires absolute bore inlet pressure.

Source syntax parsing and an independent 90-case numerical elimination/work
oracle passed during integration; neither executes the Rust implementation.
Native tests and example audio remain pending. The declared geometry, source
fit and material validity gates remain mandatory.

## Finite source gas on continuing bonded layered strings

`acoustic_realize::concentric::finite_gas` couples the existing layered-string
thermal step to the original `IsochoricGasParcel`. All annuli keep their own
material cards, caloric cells, temperatures and source-derived electrical
paths. The gas retains its original species, thermochemistry, transport
sources, mass and volume. Its initial complete snapshot must agree with the
string's accepted ambient.

The existing bounded scalar feedback owner solves signed heat exchange with
the outermost annulus around the existing backward-Euler layer network. Each
proposal begins from the same accepted histories. The conductance callback
receives the candidate gas and exposed material states, so source-driven
conduction or convection can share this transaction. The original parallel
Ohmic owner optionally supplies current/voltage heating at each trial
temperature. Solid bending loss and per-layer reversible thermoelastic heat
remain in their original cells; the actual native aerodynamic loss enters
the gas once. Interface heat is debited from the exposed solid and credited
to the gas.

Acceptance independently measures caloric solid storage, mechanical internal
energy (kinetic plus Helmholtz plus `T*s`) and gas internal energy. The
external account includes force work, electrical supply and other declared
solid/gas heat; signed interface heat and drag cancel internally. Damping
roundoff stays an explicit numerical account. Endpoint gas thermochemistry,
EOS and transport rebuild drag and pressure observation while preserving
modal displacement, momentum and observer acceleration history. Final
publication advances the string, thermal network and gas together;
source-domain exits, failed balances, work limits and cancellation leave all
three accepted owners unchanged.

The heat-only path keeps species fixed. The optional
`ThermalConcentricStringRuntime::step_with_isochoric_gas_nozzle` composition
then exchanges species and upstream enthalpy with a second finite gas parcel
using the original source-resolved nozzle. Both directions preserve species
inventories and total extensive gas storage; transported flow work and jet
energy are not deposited again as heat. Final drag and pressure use the
conserved mixture on the continuing layered motion, with one sample and the
previous accepted observer acceleration. Heat-stage gas/source receipts remain
explicitly separate from the subsequent flow endpoint. All four owners publish
together.

These are well-mixed fixed-volume gas inventories. Heat followed by flow is
an explicit first-order partition and needs physical timestep refinement.
No chamber acoustics, resolved jet force, pressure-volume work, reactions,
evolving geometry or complete bowed-instrument body/radiation feedback is
supplied. The inherited fixed-section, solid-caloric, small-strain and source-law
limits still apply. The compact pressure observer is diagnostic; it is not a
second heat source.

The source-bound conduction adapter also checks the actual exposed outer area,
original conductivity claim and conditions, and the native lumped Biot limit.
A pre-window point-load initializer uses the original modal equilibrium owner;
existing motion or pressure history cannot be overwritten by initialization.

Three focused regression definitions independently reconstruct mixture heat
capacity, layer Ohmic current and complete stored energy; compare composition/
pressure changes at identical geometry; and check late refusal, cancellation
and exact retry. Production and regression source passed Rust syntax parsing.
Native compilation, numerical execution and audio rendering/listening were not
available in this session; these tests do not establish instrument fidelity.


## Original finite gap channels on complete DG exterior groups

`DgPortLoad::with_terminal_channels` binds explicit fixed planar slit channels
or rigid closures to actual external ports, including finite chimney outlets.
Each slit is the existing `SlitChannel` transverse no-slip momentum owner with
the original local `GasState` and common clock. Density, viscosity, inertial
history and resolved wall layers enter the simultaneous mean-pressure solve;
geometry or material substitution changes the physical source feedback.
The admitted source remains affine. All exterior mutual entries and rational
histories remain intact. A closed channel has exactly zero held flow and an
unconstrained pressure reaction, including when every terminal is closed;
there is no infinite resistance, deleted exterior channel or reset history.

Each original gap owner advances privately at the shared held flow. Its actual
FV flow, constitutive pressure, storage, viscous heat and backward-Euler
numerical dissipation participate in the unchanged relative work gates before
bore, branch, gap and rational states publish together. Channel kinetic
storage and both explicitly distinguished losses join the aggregate source
ledger once. `DgTerminalFrame::channel_work_j` is signed exchange into that
owner; exterior work subtracts it from the independently integrated DG work.
Negative channel work is permitted when stored fluid energy returns.
`last_terminal_channel_frames` retains the complete original channel step;
numerical dissipation is never relabeled as physical heat. Existing actual
outer flow drives the exterior observer once. Cold joint-mouth attachment
retains and recomputes the same gap condensation.

Admission checks fixed thin geometry, unique declared terminal ownership,
local duct area, the complete explicitly budgeted transverse field count
(hard limit 4096), and a supplied resolved-cell Mach ceiling in `(0,0.3]`.
Runtime Mach refusal or any late work/flow failure leaves all accepted
histories and observations intact. Empty channel declarations retain the
original direct-terminal path. This is a fixed, laminar, incompressible planar
port reduction with original backward-Euler channel stepping; it does not
resolve axial entrance/exit jets, lateral walls, annular pad radiation,
turbulence, moving-pad swept volume/contact or full instrument qualification.
The exterior surface geometry is unchanged and remains independently explicit.

Focused regression definitions compare one/two finite gaps, both original DG
kernels, a finite material-gas chimney and the full exterior history against
independently assembled uncondensed pressure equations. They also cover exact
closure, all-closed source response, signed kinetic-energy return, material
and gap substitutions, joint-mouth construction order, aggregate source work,
late refusal and exact retry. Rust syntax parsing passed; native compilation,
numerical test execution and acoustic qualification remain pending.


## Continuing finite-channel momentum through source gas changes

`SlitChannel`, `CircularChannel` and their shared `FixedChannel` owner accept
`rebind_gas_controlled` at fixed geometry and sample period. The original
finite-volume momentum owner reassembles density-dependent cell masses and
source-viscosity face conductances, preserving each cell's axial momentum.
Consequently a uniform density change scales velocity by `rho_old/rho_new`
and kinetic storage by the same ratio. A viscosity-only change alters the
subsequent diffusion while preserving the instantaneous field and storage.

`ChannelGasUpdate` reports original and updated kinetic energy and their
signed difference as parameter work. Work limits, representability and
cancellation are checked before publishing the gas or any physical history.
This zero-duration intervention consumes no flow sample and creates no
pressure work, wall heat or numerical dissipation. `last_step()` remains the
historical accepted exchange; `stored_energy_j()` reads the actual current
kinetic storage, including an accepted gas intervention.

Planar and annular cells use the same transfer and accounting implementation.
The operation does not model background mass/heat transport, moving walls,
compressible channel storage or turbulence. The external gas owner remains
responsible for its background balance and the reported kinetic work.
Focused native regressions use independent cell momentum and kinetic-energy
calculations, changed-viscosity continuation and transactional refusal.

## Original cavity gas diffusion in the reciprocal string bank

`with_material_shell_cavity_thermal` projects scalar heat diffusion on the same
original tetrahedral Helmholtz FE into the existing complete mechanical and
material-memory propagator. `HelmholtzFe::isothermal_diffusion_modes` restricts
the original consistent K/M pencil to homogeneous Dirichlet traces on selected
solid walls; actual compact mouths remain natural zero-diffusive-flux faces.
The original sparse inertia/slicing owner retains the entire declared positive
eigenwindow. A source ceiling refuses, never truncates. Its bounded eigensolver
call retains the original cancellation limitations; assembly/projection polls.

The coherent source gas supplies kappa/(rho*cv), cv=cp/gamma, adiabatic modulus
rho*c² and isothermal modulus rho*c²/gamma. Complete L2 pressure/temperature
cross projections form reciprocal compression-memory rows. Positive equilibrium
compression is factored by the existing Cholesky owner, and one first-order
Maxwell state per active heat eigenmode enters the existing exact network.
All modal cross terms and actual neck volume coordinates remain coupled. The
instantaneous stiffness equals the original adiabatic cavity. Missing thermal
projection remains explicitly adiabatic, with a declared infinity-norm bound;
no modal-Q fit, eigenvalue clipping or duplicated compliance is introduced.

Existing source/body/layered history stays a prefix, followed by gas history.
New gas states start unrelaxed so initial mechanical compression is preserved.
Cold gas/torque/axial rebuilds retain thermal projections and states. Accepted
single/multiple-bow pressure and neck observations include that exact memory.
Coordinate-only physical pressure observations refuse on a memory-bound cavity;
the separately named adiabatic reference is not a physical pressure readout.
The original work/energy owner counts thermal availability and dissipative loss
exactly once, together with material memories and exterior acoustic loading.

This is linear gas conduction to fixed-temperature walls, not complete
viscothermal boundary-layer physics: no no-slip shear, wall effusivity/finite
thermal inventory, aperture enthalpy advection, molecular relaxation or measured
instrument agreement is claimed. The gas/solid state remains frozen. The FE,
heat-mode window, compression-projector remainder and pressure basis require
convergence; a coarse FE does not resolve audio-frequency thermal skin depths.
Focused tests cover analytic star-tet heat rates/projections, coherent gas and
length substitutions, reciprocal pressure/force/storage/loss, strict budgets,
original bank history and reconstruction. Native Rust execution remains
required; syntax and independent matrix identities are not runtime proof.

## Refined sealed multi-chamber material assemblies

`cavity_phs::piston::thermal::boundary::ports::network::material::state::subcycling` extends the existing co-owned `MaterialExchangeState` with `step_subcycled_candidate`/`advance_subcycled` and their source-bearing `with_gas_law` counterparts. Each real subinterval advances the original circuit charge/flux, every finite cell enthalpy, declared chamber reactions and independent piston walls together. Subsequent chemical stages therefore see actual material wall heat and compression within the same requested interval. Existing thermal contact resistance, electrical contact losses, radiation and caller-selected gas convection/conduction continue through their original source owners and nonlinear equations.

`MaterialExchangeSubcycling` reuses `ReactiveThermalSubcycling`'s gas/species/enthalpy scales, scalar comparison and 256-substep/511-evaluation bounds, adding independently dimensioned charge and flux tolerances. It compares complete N and 2N trajectories and accepts the finer actual trajectory without extrapolating physical state. Each species and reaction channel is compared separately; every finite cell uses its own enthalpy coordinate, including latent plateaus. Node potentials and branch currents remain the original circuit's resolved algebraic efforts, not new differential stores. Nested chemistry/network/thermal/momentum solvers retain their original explicit work limits. All accepted source reports are retained within the bounded fine trajectory; discarded trajectories publish nothing.

Electrical port inputs, exterior pressures and applied forces are held at their caller values. Each chamber's `heat_j` is integrated external heat, apportioned at constant power over actual representable subinterval boundaries, with a final remainder. Thermal `source_w` and existing boundary/source integrals are not divided again. Numerical cell-energy guesses are scaled to the actual subinterval without changing a physical source. Chemical formation energy, Joule heat, wall exchange, contact conduction and closed-enclosure radiation are never credited twice.

The macro audit independently reconstructs all finite-cell enthalpy change, circuit-owner initial/final Hamiltonians, each sealed gas internal-energy change, piston kinetic changes and accepted source/mechanical work. Original electrical, thermal, circuit/thermal and complete-network energy tolerances remain enforced. The sum of accepted maximum cell residuals conservatively bounds each fixed cell's accumulated residual at the original cell tolerance. Original chemical energy/entropy/species and independently reconstructed pressure/applied/support impulse bounds also remain in force. Each gas receives one macro port tick with its actual endpoint entropy and integrated heat entropy Σ(Q/T), including sign-reversing heat histories. Final publication moves the entire assembly only after all checks and cancellation gates succeed.

This is bounded temporal refinement of the existing first-order sealed-chamber split, not a monolithic implicit chemistry solve or a reacting compressible-flow model. Cell masses, thermal/circuit topology and chamber identities remain fixed; the `state::flow` transport stages are not composed by this API. Source-domain refusal, failed numerical work and cancellation preserve the complete original assembly and its port/history state.

Focused regressions in `tests/material_exchange/subcycling.rs` exercise a data-only emissivity change reaching reaction progress within the interval, real circuit/contact/radiation source histories, independent total physical energy and impulse balances, integrated heat versus held electrical-current/force semantics, nonzero-clock resume, charge-controlled refinement, and late source/work-limit/cancellation rollback. Fixtures are synthetic and establish integration behavior rather than empirical combustion calibration.

## Consistent-inertia solids with original tensor material memory

`tensor_material::TensorMaterialBody::from_mass_matrix` admits the complete
symmetric positive-definite inertia in the caller's displacement coordinates.
Momentum is `M v`, velocity is `M^-1 p`, and the kinetic Hamiltonian retains
the full quadratic block. The existing `fs-la` Cholesky owner supplies the
inverse; a dimension-scaled floating-point identity-residual check rejects
unrepresentable inertia, including dependent bases that leave a tiny positive
pivot. It never regularizes mass or removes coordinates. Existing diagonal
inertia and plane-stress paths keep their own original arithmetic.

`tensor_material::solid::TensorSolidBody` connects
`fs-solid::linear3::visco::ProjectedViscoSolid` to this same full-memory
`TensorMaterialBody` and the existing `fs-phs` exact-ZOH propagator. The
source tetrahedral basis owns consistent mass, strain maps, supports and
force/velocity transformations together. Every material point retains its
equilibrium stiffness, Kelvin matrix and all Maxwell branches. Initial
constitutive prehistory remains explicit. Held nodal forces are projected by
the same basis transpose; the reported mean nodal velocity expands the same
work-conjugate port output used for integrated supplied work.

`TensorMaterialHistory::advance_with_material_heat` independently integrates
each original point's Kelvin and Maxwell resistance channels along the same
complete coupled trajectory using `fs-phs::integrated_resistance_loss`.
`TensorSolidHistory::advance` aggregates those channels to their original
elements. No work-minus-storage residual is assigned as material heat.
Both point and element totals must agree with the original total dissipated
energy within the caller's explicit positive absolute attribution tolerance.
Failed propagation, heat attribution, nodal expansion or cancellation preserves
the entire original history, clock and accumulated heat.

Existing material-point/state-coordinate and dense-propagation capacities remain
in force; oversized requests refuse without dropping modes or memories.
Geometry, basis and material temperature/composition are fixed. The returned
element heat is available for a thermal owner but this API does not itself
advance caloric state. Separate point heat for plane-stress histories is
explicitly refused; their existing total-loss dynamics remain available.
There is no claim of nonlinear mechanics, acoustic radiation, spatial/modal
convergence or empirical material calibration. Five focused G0/G3 regressions
in `tests/tensor_solid.rs` exercise off-diagonal inertia and analytic force work,
Tet4/Tet10 shear/mass, Maxwell prehistory, regional loss attribution, constrained
nonorthogonal work maps, dependent-basis refusal and transactional failure.

## Retained inertial motion during native finite-gas heating

`heated_valve_sound --inertial-history MAX_STEP_S MASS_DAMPING_PER_S` retains
each accepted finite-gas/source state and advances one original sparse Newmark
trajectory from stress-free rest on those exact times. Nodal displacement,
velocity and acceleration persist while source stiffness and expansion refresh
at fixed reference mass, mesh and supports. `PREFIX.motion.csv` records actual
nodal motion and global kinetic/shifted elastic energy, signed external,
expansion and stiffness work, damping and cumulative mechanical residual.

The explicit step ceiling is an admission limit, not a temporal-error
certificate. This is prescribed one-way thermal loading: parameter work and
damping are not silently exchanged with caloric storage. Pressure CSV/WAV retain
the existing frozen-state acoustic experiments. All mechanical admission and
sound calculations precede output creation. Focused native tests check the
retained Newmark kinematic identities, shared heat/motion clock, actual inertia,
source-stiffness work, energy balance, replay and refused/cancelled inputs.

`LocalMaterialWall::bind_dg` carries the same source-resolved thin-cylinder
mass, elasticity, Kelvin-Voigt viscosity or complete generalized-Maxwell
relaxation into the native DG duct time operator. `MaterialWallDg` retains the
original material receipts, geometry
reduction and declared frequency applicability beside that operator. Material
data changes affect pressure, reciprocal wall motion and stored energy through
the existing gas/port/time owners. Every positive Maxwell source arm retains
its radial reaction pressure, elastic storage and nonnegative loss in that
same operator. The exact exponential advances all memories together, so the
network midpoint clock restrictions do not apply to this route. The source
band remains an applicability statement, not a
filter or broadband certificate; uniform gas/cylinder, small displacement and
frozen source state are the bounded model. The source-data regression changes
only Young modulus and compares actual time-domain pressure and wall motion
with the independent sealed-tube oscillator while checking total energy.
A further source-state regression changes only a relaxation-time claim and
checks actual DG pressure, retained branch energy and material loss; the
unused Kelvin property is not counted again. Native execution is pending.

## Nonlinear induction coupled to finite material enthalpy

`MagneticThermalNetwork::step_nonlinear_eddy_candidate` advances the original
`PlanarNonlinearEddy` field and `DepositionSimulation` on one accepted time
interval. Each thermal Newton trial starts from the same accepted magnetic
potential, resolves the same material card's conductivity at the trial endpoint
temperature, and solves the source-defined reversible H(B) induction problem.
The original conductor-region projection supplies exact triangle-integrated
sigma*E^2 losses to the existing thermal cells; the shared material network
supplies conduction, authored sources, boundaries and finite caloric storage.
The coupled candidate publishes both histories together. Cancellation, missing
routes, source mismatch or nonlinear failure leaves both inputs unchanged.

Only physical Joule heat enters material enthalpy. The field step's nonnegative
backward-Euler magnetic convexity loss remains separately visible as
`report.electrical.numerical_dissipation_j`. The shared driver audit checks
magnetic storage change + thermal storage change + numerical dissipation
against impressed-current work and other external heat. Numerical loss is not
an empirical core-loss model and is never deposited as heating.

Source identity, property basis, requirements, nonthermal coordinates, fixed
geometry, original triangle volumes and complete conductor/thermal routing are
re-admitted at every trial. The complete H(B) source path must remain unchanged;
a temperature-dependent H(B,T) source refuses thermal feedback even if selected
endpoint values happen to agree. This path does not supply magnetic entropy,
thermal parameter work, hysteresis, remanence, moving conductors or terminal
circuit feedback. Caloric charts are fixed and fully solid in this bounded
coupling; unsupported phase changes refuse through the existing material owner.

`tests/nonlinear_eddy_enthalpy.rs` specifies an independent scalar-bisection
oracle for source-only changes to saturation and temperature-dependent
resistivity, continuation from the accepted field, physical heat versus
numerical loss, complete energy balance, and transactional failure. Full native
integration execution is required in a build environment that can compile the
pinned constellation; the production-leaf checks alone do not establish it.

## Reacting finite gas and spatial material-solid heating

`SolidEnthalpy::step_with_reacting_isochoric_gas` advances the original declared
`GasReaction::parcel_candidate` at fixed gas mass, volume and formation-inclusive
total internal energy, then delegates the existing joint gas/spatial-cell
backward-Euler owner over the same physical interval. Reaction products and their
material transport enter the thermal solve. The optional original
`PreparedMaterialNetwork` keeps internal conductivity and gas-facing material
layers at the actual endpoint temperatures. The returned native solid state
feeds the existing deformation, stress, prestress and transient-solid owners.

Original-to-final gas energy plus cell enthalpy minus authored external heat
must close before either stage publishes. Chemical formation energy is already
in gas storage and is never added as a separate heat source. The gas temperature
increment bound covers both the chemical intermediate and final state. Failed
chemistry, thermal iteration, source admission, events or cancellation preserve
the original gas and thermal histories.

This is a first-order chemistry-to-thermal operator split for a well-mixed,
fixed-mass, fixed-volume gas. The isolated chemical stage retains its own entropy
admission; this is not a whole-system entropy-production claim. Flame propagation,
gas traction, moving boundaries, melt mechanics and reciprocal thermoelastic
heat remain outside this path's scope.

`tests/solid_reactive_gas.rs` specifies independent spatial heat and native 3-D
solid displacement/stress references under data-only kinetic, caloric and product
transport changes, cooled-species continuation with material k(T), and late
thermal failure/replay/cancellation. Source parsing and review do not substitute
for executing these full integration regressions.

## Circuit-driven nonlinear magnetic fields and finite material heat

`WindingThermalNetwork::step_nonlinear_candidate` closes the original native
circuit, nonlinear magnetic field and finite enthalpy histories over one
accepted interval. Source-defined H(B) determines the actual field and winding
flux linkage; the circuit voltage determines current. The same material cards
supply winding and core resistivity at each proposed thermal endpoint.
Every thermal trial restarts both accepted electrical histories and uses the
existing material-network conduction and caloric residual.

Linear and nonlinear winding consumers share admission of original magnetic
triangle volumes, explicit heat destinations and native graph-resistor routes.
Core eddy regions and stranded winding regions must be disjoint. Their physical
Joule losses and every original graph resistor's loss are deposited exactly
once. Private winding-port exchanges cancel inside the combined field/circuit
owner: only external source work enters the global thermal energy account.
Backward-Euler field and explicit end-inductance numerical dissipation is
reported separately and never becomes material heat.

The complete H(B) law, fixed geometry, source identity and original caloric
instances remain pinned. Thermally mapped magnetic cells require wholly solid,
constant-density caloric charts. These entry points retain fixed magnetic
storage; the thermal-magnetic entry below supplies its entropy/work owner.
Existing thermal
events, conduction, source domains, work budgets and late cancellation remain
authoritative. All three accepted histories are unchanged on refusal.
This is a first-order, reversible isotropic magnetoquasistatic coupling; it
does not add hysteresis, moving conductors or thermal expansion.

G1/G3/G4 regressions in `nonlinear_winding_enthalpy.rs` compare two accepted
steps against an independent scalar voltage/current/enthalpy reduction,
change only source saturation data, verify endpoint resistivity and every
physical heat contribution, and distinguish external work from private port
work and numerical loss. They also exercise incomplete routes, trial budgets,
late cancellation, clock mismatch and deterministic retry. The integration
tests are added and parsed; full native execution is not yet claimed.

### Ordinary resistor material feedback in nonlinear field circuits

`WindingThermalNetwork::step_nonlinear_material_candidate` extends that SAME
residual to source-bound ordinary graph resistors. Complete `CellResistor`
bindings and explicit conductivity/resistivity bases use the original
`PreparedMaterialNetwork`; its card, caloric instance, phase, geometry and
source-selection checks remain authoritative. Accepted native resistance and
any retained source/geometry snapshot must match the actual starting thermal
state. Each trial refreshes the original selected sources at its proposed
endpoint, then invokes the native resistance-only circuit reassembly before
solving the existing nonlinear circuit/field equations. The accepted circuit
and report retain the resulting source bundles across continuation steps.

Every whole-cell resistor binding must match its existing graph heat route and
occupy a separate body from the cells covered by magnetic triangles. This
avoids representing the same physical conductor volume with both bulk graph
loss and winding/core loss. Existing material-conduction links couple those
separate thermal bodies. Native graph-local Joule heat is deposited once through
the original routing; source resolution adds no second heat input. Resistance
changes store no parameter energy, and field backward-Euler loss remains
numerical dissipation. The legacy numeric-resistor entry refuses omission of
already attached material bindings, rather than freezing them while heating.

Focused manufactured regressions discharge a finite capacitor into the original
nonlinear field and two distinct thermal bodies. An independent scalar current
solve and exact two-body thermal equations check temperature-dependent bulk
resistance, heat transfer, charge drain, field energy and physical/numerical
loss separation. Changing only the resistor's source curve changes both actual
current and body temperatures. Source-domain, missing-binding, work-budget and
late-cancellation refusals preserve all accepted histories. These are discrete
coupling checks, not measured resistor or magnetic-device validation.

### Voltage-driven magnetic temperature and finite caloric feedback

`WindingThermalNetwork::step_nonlinear_thermal_candidate` closes sourced mu(T),
winding/core resistivity, optional material graph resistors and finite caloric
storage on that same native field/circuit/enthalpy residual. It retains actual
nodal field, terminal current, linkage, remaining circuit charge/flux and clock;
every trial starts from the accepted histories. The native field owner rebuilds
the original mesh at the proposed material state and retains its new operator
for continuation. Complete source cards, claim selections, fixed coordinates,
material volumes, thermal receivers and original topology remain bound.

The shared `magnetic::caloric::PreparedMagneticCaloric` implementation now serves
both this winding consumer and the original induction thermal consumer. It uses
the existing `ThermalMagneticMaterial` for field Helmholtz energy psi, entropy S
and internal energy U=psi+T*S. Each native cell's signed parameter work P supplies
reversible heat Q_rev=-P-delta(T*S) to its original caloric receiver. Physical
Joule loss is routed once through existing winding/core/graph destinations;
backward-Euler numerical loss is never deposited as heat.

The ordinary report retains native free-energy/parameter-work closure. The
additional `caloric_balance()` independently checks field internal energy,
native circuit storage, declared winding end storage and caloric energy against
external electrical work, external heat and numerical loss. Both source samples
and the accepted physical balance survive `into_parts` in the returned report.
No constitutive work is silently assigned to an external thermal reservoir.

This path requires linear isotropic permeability with an admitted positive
affine mu(T) over each declared interval, fully solid fixed-density magnetic
volumes and positive combined caloric/field heat capacity throughout the declared
temperature and flux-density envelope. Fixed external magnetic regions require
temperature-independent permeability. H(B,T), hysteresis, remanence evolution,
moving geometry and thermal expansion require their own thermodynamic owners.
Source, capacity, field-envelope, solver and cancellation refusal publish none
of the proposed histories. This is first-order time integration, not an
unconditional passivity or experimental material-validation claim.

`nonlinear_winding_enthalpy::magnetic_thermal` adds independent voltage/current
and physical-energy references for data-only permeability changes, nonzero end
inductance, continuing intervals, constant-law equivalence, replay and atomic
refusal. Parser/source checks and independent reference-equation evaluations
are reported separately from native execution of these regressions.

## Source gas fields and source walls in the native DG operator

The authored material-modal DG valve and DG flue loaders accept the same
optional `dg_wall` record after the complete `dg_tube`/`dg_bore` block:
`dg_wall isotropic PACK T OMEGA H AXIAL FMIN FMAX LOSS`. `AXIAL` is
`stress-free` or `strain-fixed`; `LOSS` uses the existing elastic,
Kelvin-Voigt or source-property Maxwell grammar. The declared material band
must cover the propagation band. `dg_wall rigid` and omission retain the
original operator. Missing source data, unsupported geometry, thin-wall
violations and complete-state budget overflow refuse before rendering.

A complete physical source partition uses `dg_wall regions N`, followed by
`dg_wall_region X_BEGIN X_END isotropic PACK T OMEGA H AXIAL FMIN FMAX LOSS`
for each region in inlet-to-outlet order. Coordinates are in metres from the
original bore inlet. Regions must be contiguous, cover the full bore and
coincide with actual DG element faces (only floating-point roundoff is
tolerated). Refinement repeats the same source assignments; cutting a material
interface refuses instead of averaging the two solids. The common applicability
band is the intersection of the source bands. Up to 256 regions share an 8 MiB
unique-artifact budget and the original complete-state limit. Repeated paths
reuse their immutable artifacts, including paths shared with reed blades.
The same CLI and pressure/PCM renderer retain every region's geometry, source
identity and local material memory; no second propagation law is introduced.

Both loaders bind the original `LocalMaterialWall` or its
`RegionalMaterialWallDg` composition before preparing the same DG exponential. Density, stiffness, viscosity and every positive source
Maxwell arm affect the actual pressure feedback, radial wall motion and
retained storage. Prescribed gas fields and independent exterior gas remain
their original owners. The reed loader shares a single normalized artifact
cache between its blades and wall; CLI metadata retains the solid source and
model scope. This remains a frozen, isotropic thin cylindrical main-bore wall;
side chimneys stay rigid, and axial shell waves and exterior wall radiation
are absent. Source-swap, gas-field, callback, original CLI/PCM and refusal
regressions are added; full native verification remains required.

`LocalMaterialWall::bind_dg` composes its original elastic, Kelvin or
source-Maxwell material with a prescribed nonuniform `AxialGasField` on
identical circular segments. The DG owner uses the already assembled pointwise
gas pressure masses and independent exact geometric wall masses. Gas
composition affects the actual pressure/flow operator and reciprocal wall
motion; the source solid still owns inertia, equilibrium stiffness and every
positive relaxation arm. Original gas/source receipts, pressure ports, stored
energy, dissipation and state-size limits are retained.

`bind_dg_controlled` supports cancellation before candidate publication and
reports the existing typed cancellation error. Original source/domain and
thin-wall refusals remain authoritative. This is cold assembly of a frozen gas
field and frozen wall, with the same uniform-cylinder and material-band limits;
it does not transfer an already vibrating material or solve gas
advection/thermal evolution.

Source-pack regressions change only downstream composition while keeping inlet
gas, geometry and wall source fixed, then check propagated pressure, wall motion,
Maxwell storage/loss and total passive energy. A separate source-solid
substitution uses the same gas field. Missing transport, out-of-domain
temperature and cancellation refuse. These full integration tests are added and
parsed but not natively executed. Native exact-production matrix tests cover
independent heterogeneous sealed-pocket dynamics, physical energy/force for
noncommuting mass roots, and refusal/cancellation; they do not certify the
complete material-query/DG exponential path.

## Finite spatial enthalpy drives continuing reduced solid motion

`solid_enthalpy::valve::motion::HeatedValveMotion` binds one finite thermal
inventory, pinned solid source/reduction and existing `ThermalValveModalMotion`.
Each interval advances the native held-pressure motion at the accepted source,
advances the original spatial thermal owner over the same physical duration,
rebuilds compatible thermal equilibrium and tangent modes, and uses the existing
mass-inner-product transfer for the same displacement and conjugate momentum.
The changed equilibrium contributes its actual displacement offset. Mode order,
sign or count never substitutes for physical state identity. The new observer is
formed on the actual endpoint basis and evaluated after transfer.

All candidates remain private until source admission at every thermal substep
and event, equilibrium/reduction, native modal limits, projection/work limits,
endpoint observation and cancellation succeed. Heat, motion, basis, observation
and sample count publish together. The absolute thermal clock is the stored
initial origin plus the native modal elapsed clock; incompatible intervals
refuse. The original caloric heat audit, pressure-step work/loss and signed
incremental modal parameter work/discarded-state accounts remain separate.
Projection loss and coefficient-change energy are never physical heating.

This is a first-order prescribed-thermal split on fixed reference mass, mesh,
supports and pressure footprint. It does not close reciprocal thermoelastic
heat, dissipative heat return, bore pressure feedback or radiation reaction;
the supplied narrow-band pressure transfer is an observation. G1/G3/G4 tests in
`tests/solid_enthalpy/valve/motion.rs` compare an independent heated one-DOF
oscillator, a heat-capacity-only material swap in actual motion, and late
observer/transfer/cancellation refusals followed by deterministic retry. Full
native runtime execution remains required for numerical qualification.

## Finite source gas drives an elastic pressure boundary

`cavity_phs::modal::GasValveMotion` joins the original source-bearing
`ClosedGasParcel` to `ThermalValveModalMotion`. The source mesh, consistent
mass, prestressed tangent modes and pressure footprint remain authoritative.
The initial physical gas volume anchors an affine volume chart in the actual
retained modal coordinates; volume is not reconstructed from a separate
effective mass or piston area.

One native sample solves finite gas heat and pressure work together with the
original exact held-pressure modal transition. The existing Newton-Krylov
solver sees a dimensionless volume residual and an analytic Jacobian built
from the native held-force map and gas midpoint pressure derivative. Every
trial starts from the same accepted histories. The original mixture's caloric,
EOS and transport laws are resolved at the gas stages; no nominal air constants
or fixed gamma are substituted.

Publication requires independent physical-volume, pressure-transfer and
enclosing energy bounds. Exterior pressure work uses the actual wall sweep.
The enclosing balance is gas internal-energy change plus modal-energy change
plus native modal damping loss plus exterior work minus supplied gas heat.
Structural loss is not credited as gas heat. Source refusal, failed nonlinear
iteration, work limits, cancellation and underresolved gas-loaded stiffness
leave both participants and the sample clock unchanged. The original modal
Nyquist guard and existing midpoint resolved-band limit apply to a conservative
bound that includes the gas rank-one stiffness, at the old, midpoint and final
gas states.

The base `step` uses a closed, spatially uniform, frozen-composition gas and a
fixed structural material/geometry epoch. Its pressure/volume feedback requires
time refinement; exact structural ZOH does not make the full nonlinear
interconnection exact. There is no new wall heat transfer, reciprocal
thermoelastic heating, reactive/nozzle flow, finite-deformation geometry or
automatic exterior radiation-transfer refresh. Changing the structural source
requires its native state transfer and explicit geometry/parameter work before
binding a new volume chart.

## Material-driven longitudinal piezoelectric sensor/actuator

`piezoelectric` composes the original `fs-material::piezoelectric` source law
with a rod's reference area/length and the shared `fs-phs` stepping and port
owners. The fixed-root axial shape `u(z)=x*z/L` gives consistent mass
`rho*A*L/3`. Combined energy is `A*L*w(x/L,Q/A)+p^2/(2m)`; its derivatives
provide force, voltage and velocity. Geometry and material data therefore
change inertia, open/short stiffness, capacitance, direct sensing and converse
actuation together, without independent sensor/actuator gains.

`PiezoelectricBody::step_candidate` advances extension, momentum and free
charge in one native descriptor solve. A held current includes open-circuit
zero current. A held ideal voltage uses an algebraic source-current coordinate
and includes a short circuit at zero voltage. Initial charge must already be
compatible with a voltage clamp: native zero-time admission may not project
any conserved coordinate. Switching impulses need a separate circuit model.
Accepted endpoints satisfy the declared small-signal strain/field envelope,
explicit voltage tolerance and independently reconstructed mechanical/electric
work balance. Cancellation, source/domain refusal or failed work balance
publishes no partial state. Cancellation is polled around the bounded 3/4-state
native solve; no interruption inside its factorization is promised.

`PiezoelectricRod::port_hamiltonian` exposes the same energy as an ordinary
two-port component: inputs [force,current], outputs [velocity,voltage], state
[extension,momentum,charge]. Existing native interconnections can compose it
with physical electrical/mechanical storage. G1/G3/G4 synthetic regressions in
`tests/piezoelectric.rs` exercise solved data-only material substitutions,
direct/converse responses, distinct open/short mechanical trajectories, actual
R-L circuit exchange, work conservation, and atomic refusal/cancellation.

This single-element longitudinal model is fixed-state and reversible. It adds
no empirical coefficients, resolved axial waves, transverse mechanics,
fringing, material hysteresis, heating, or measured material validation.
Amplitude limits are caller-admitted model restrictions, not sourced strength
or breakdown measurements. It is a bounded MR59 implementation, not closure
of that bead's full tensor/electrostatic/structural qualification.

## Regional material walls in native duct propagation

`bernoulli_aperture::wall::material::RegionalMaterialWallDg` retains one original
`LocalMaterialWall` per DG element, including its material state, thickness,
restraint and explicit loss law. It compiles independent surface mass, stiffness
and Kelvin--Voigt resistance into `fs-duct`'s existing reciprocal gas/wall
operator and intersects the original source bands. Different wall materials
therefore change both pressure propagation and wall motion in the same native
time solve, including spatial gas composition and pre-existing gas memories.

Exactly uniform assignments retain the original full Maxwell-spectrum route.
Heterogeneous assignments carry every original positive Maxwell arm on its own
element's nodes, alongside the original local mass and equilibrium stiffness.
Elastic, Kelvin--Voigt and hereditary regions may coexist; no common spectrum
replaces their sources. Gas memory, ports and wall histories advance in one
native operator, and physical pressure observations are zero outside each arm's
support. Material boundaries align with element boundaries. The existing thin cylindrical
local-membrane model supplies no axial shell waves or stress transfer across
material joints. Cold attachment and the live elastic/Kelvin--Voigt transfer
described below are distinct operations. Focused G1/G3 tests reconstruct local force/energy relations and
change only one wall assignment through the actual pressure/motion evolution.

## Spatial Ohmic current and finite material-temperature feedback

`circuit_enthalpy::conduction::ConductionThermalNetwork::step_candidate`
connects `fs-em::conduction::PlanarConduction` to the existing material-bound
enthalpy/Newton owner. Each electrical triangle has one explicit heat receiver;
mapped triangles must cover that receiver's entire physical volume. The same
card, caloric law instance, source selections and fixed external coordinates
remain bound. All electrically mapped caloric charts have constant density
compatible with the retained mass and geometry. Confined non-solid charts use
the explicit phase-coordinate admission below. Nonlinear caloric storage and
temperature/phase-dependent resistivity/conductivity remain supported.

Every thermal trial re-queries each original electrical source at its proposed
temperature and declared liquid mass fraction, rebuilds the native FEEC
conductivity operator and solves the held
electrode voltages. Actual spatial current redistribution and volume-integrated
`E.J` determine local heat; no uniform-resistor reduction supplies it. Original
thermal links resolve their conductivity from the same endpoint material state.
The shared backward-Euler residual accounts for these local losses, thermal
conduction and authored boundaries together. Independently recovered terminal
work must close against the deposited heat and combined stored-energy change.

The candidate returns both accepted thermal history and the rebuilt endpoint
electrical problem, together with the exact field observations and original
linear/nonlinear reports. Initial electrical/thermal state mismatch, unsupported
sources, incomplete correspondence, unaccepted balance, budgets and cancellation
return no candidate. All nonlinear trials leave original inputs unchanged.
`maximum_circuit_trials` bounds field evaluations on this shared driver;
`ConductionControls` governs the field solve, and the unused circuit-DAE controls
do not govern it. Events obey the original enthalpy interval rules.

This is fixed-geometry quasistatic DC with admitted scalar/tensor coefficients,
held electrode voltages and
first-order thermal time integration. It does not add induction, dielectric
storage, thermoelectric terms, contact resistance, density-changing phase or moving-mesh
work. G1/G3/G4 synthetic tests in `tests/conduction_enthalpy.rs` exercise layered
current/heat, an independent thermal reference, nonlinear source-only feedback,
continuation, material substitution, gauge/current reversal and atomic refusal.
They verify the declared numerical model, not empirical material qualification.


### Spatial Ohmic current-drive feedback (MR43)

`circuit_enthalpy::conduction::ConductionThermalNetwork::step_driven_candidate`
accepts a held `ConductionDrive::Voltage` or a held total two-terminal
`ConductionDrive::Current`. The existing `step_candidate` voltage API delegates
to the same complete thermal/material residual. Current-drive electrodes are
disjoint complete exterior traces; the negative electrode is the potential
reference. Each trial re-extracts the native spatial resistance from the
actual source-pinned conductivity field, solves its electrode voltage, and
checks the recovered total current independently. This drives local E.J heating,
finite caloric storage and material refresh without a separate heat/current
integrator or a guessed uniform wire. Both current directions and zero current
retain source, geometry and trace admission. Electrical storage, induction and
source compliance limits remain outside this ideal DC boundary model.

G1/G3/G4 regressions in `tests/conduction_enthalpy.rs` cover the independent
current-heated slab equation, data-only resistivity changes, continued operator
updates, reversal/zero-current limits, matching voltage-driven BE endpoints,
invalid electrodes/current and late transactional cancellation. Source review
and parsing are not substitutes for executing these coupled tests.

## Continuing vibration in an exchanging material atmosphere (MR34)

`ThermalMaterialStringRuntime::step_with_isochoric_gas_nozzle` advances one
existing mechanical/solid/gas thermal sample, then uses the native source-backed
pressure-driven nozzle to exchange species and upstream enthalpy with a second
finite gas volume. The final conserved mixture supplies both aerodynamic drag
and compact pressure observation. Modal displacement/momentum and the previous
accepted observer acceleration survive; one epoch and one audio sample advance.

The three owners publish atomically after an independently measured combined
mechanical, solid and two-volume internal-energy balance. Transported enthalpy
is internal to the pair; flow work and thermalized jet energy are not added as
separate heating. Source conflicts, excessive withdrawal, pressure reversal,
final temperature bounds and cancellation leave every accepted state intact.
The G1/G3/G4 regressions in `tests/string_specimen/finite_gas/flow.rs` cover native
transfer endpoints, independent energy/pressure expressions, data-only species
substitution, retained motion/history, reverse flow and transactional failure.

`step_with_isochoric_gas_nozzle_conduction` supplies the original source-bound
reciprocal thermal proposal to that SAME flow/accounting/publication owner.
`FiniteGasStringFlowFrame<R>::thermal` retains the original typed thermal report,
including its conductance/source receipt and pre-transfer gas heat endpoint;
its contained string observer is refreshed to the final transported mixture.
The insulated route uses `R = FiniteGasStringFrame`. Internal conduction heat
and transported enthalpy both cancel from the enclosing energy balance.

`step_with_drive_and_isochoric_gas_nozzle_conduction` also uses the original
conductor's current and Joule work. Its final resistance/current are admitted
after the complete thermal/flow proposal; a late electrical refusal preserves
the reservoir as well as the string and local gas. All three old owners remain
alive until their replacements have been published, with no fallible work
between them. The electrical and flow adapters share proposal owners instead
of duplicating conductor, thermal, nozzle or caloric laws. Focused regressions
in `tests/string_specimen/finite_gas/flow/conduction.rs` cover powered conduction
with bidirectional transport, complete storage accounting, zero-area parity,
retained sound history and late three-owner rollback/retry.

This first-order heat-then-flow, well-mixed fixed-volume split requires timestep
refinement; it does not solve simultaneous transport and conduction. It adds no
resolved jet or wind force, chamber radiation loading or propagation through a
varying medium. The compact pressure observer remains diagnostic; no experimental
material qualification is implied.

### Distributed current injection on the common electrothermal drive

`ConductionDrive::CurrentDensity` adds native prescribed-current-density edges
and voltage return traces to the same `step_driven_candidate` residual used by
voltage and total-electrode-current drives. Each material-temperature trial
rebuilds the source-bound spatial operator and solves its original mixed
boundary problem. Source work includes the solved voltage at current traces;
the actual per-cell E.J heats the existing finite enthalpy cells exactly once.
The accepted report retains the native nodal applied currents and total boundary
currents beside the material operator and thermal history. Current edges need
not be equipotential; their prescribed density differs physically from the
existing equipotential two-terminal total-current drive. The same solid,
constant-density, fixed-geometry, first-order DC limits and atomic failure rules
apply. Focused tests add independent heterogeneous current-heating equations,
temperature feedback, continuation, material substitution and late rollback.

## Source-backed modified Arrhenius gas kinetics

`cavity_phs::piston::reactive::material_rate::MaterialRateScale::modified_arrhenius` binds explicit material-card prefactor, temperature exponent, signed molar activation energy and temperature normalization to the existing reactive gas owner. All coefficients resolve at the actual gas T/p and retain original usage receipts. The activity convention uses `A (T/T0)^b exp(-Ea/(R T))` in mol/(m3 s). The alternative normalized-concentration convention additionally resolves `c0` in mol/m3 and converts to the existing pressure-activity scale using `(p0/(R T c0))^sum(alpha_forward)`; reverse rates still come from the original NASA-9 equilibrium closure. Units and normalizations must be authored explicitly; concentration-based literature prefactors cannot be copied without conversion.

The new binding pins coefficient claim identities before evolution. `SingleClaimOnly` requires globally unique claims for every required coefficient and becomes an explicit pin set; ambiguous/observation-preferred selection requires caller-supplied pins. Source validity crossing cannot silently choose another claim, including within implicit endpoint iterations. The legacy directly resolved scalar-rate path is unchanged. Existing explicit, scalar/coupled implicit, subcycled, piston, species-flow and electrical/thermal consumers use the same channel evaluator. Formation energy is never credited again as heat; source/range/cancellation refusal publishes no proposed physical state. Positive rates that underflow/overflow refuse rather than become zero or infinity.

Focused tests in `cavity_phs::piston::reactive::tests::arrhenius` cover independent dimensional rate scaling, first-/second-order concentration conversion and reverse equilibrium, normalization invariance, actual material-piston motion and energy, implicit endpoint temperature feedback, and bad-unit/range/source/cancellation rollback including disjoint source domains. Coefficients are synthetic. The Arrhenius coefficient binding alone does not supply collider/falloff or surface chemistry; the optional material collider extension below supplies bounded third-body, Lindemann and Troe laws. Neither slice supplies a validated fuel mechanism, flame/shock model, or measured engine/ignition prediction. Program anchor: `frankensim-ext-gas-chemistry-ladder-paqh`, remains open.

### Spatial conductor, transient circuit storage, and finite heat

`circuit_enthalpy::conduction::circuit` connects a source-resolved planar
conductor to an original resistor slot in the bounded native Kirchhoff graph.
The complete positive/negative electrodes follow the branch orientation. The
native spatial operator supplies the resistance; the original circuit descriptor
advances capacitor charge and inductor flux. At each proposed thermal endpoint,
the existing material resolver refreshes conductivity from the same retained
sources, the descriptor restarts from the accepted electrical storage, and its
actual midpoint voltage drives the native spatial E/J field.

Every triangle retains its full-volume thermal receiver. Every other graph
resistor must name exactly one active receiving cell. The existing enthalpy
candidate owns the common nonlinear residual and publication. Field terminal
work is internal to this boundary; only the native circuit ports supply external
electrical work. Independent checks compare field/circuit current, local and
graph Joule heat, native dissipation, electrical storage, and the combined
thermal/electrical balance. Missing, duplicate, foreign, or double-counted heat
routes refuse. Initial source snapshots and electrical/thermal clocks must agree.
Cancellation or any late refusal leaves all accepted fields, charges, fluxes,
source selections, and thermal stores unchanged.

`ConductionCircuit::new_with_material_resistors` also binds complete material
conductors to every nonspatial resistor slot. Their geometry and source basis
come from the original `ResolvedConductor`; placeholder graph resistances do not
enter the calculation. `step_material_circuit_candidate` uses the existing
`PreparedMaterialNetwork` to re-evaluate these same sources at every proposed
thermal endpoint, alongside the spatial conductivity and heat-transport links.
Conductivity and resistivity remain distinct source bases. The accepted branch
reports retain complete conductor snapshots and their original Joule-power law.
Each material resistor occupies its whole caloric cell; it cannot also represent
the spatial body's already covered volume. Its initial geometry, complete source
bundle and state must agree with the accepted electrical model. The ideal-resistor
step refuses material-bound circuits rather than freezing their coefficients.

This rung represents one heterogeneous two-electrode body in a bounded graph.
It holds geometry, solid density, circuit topology, C/L values, explicitly ideal
graph resistances, and source-port inputs fixed over the interval. Source-bound
graph resistances evolve through the same material/thermal residual. Electrical midpoint
stepping with endpoint material coefficients and backward-Euler heat is generally
first-order as a coupled method; energy admission is not a discretization-error
bound. It does not add conductor induction/displacement current, moving geometry,
phase transitions, contact resistance, or thermoelectric transport.

The added `tests/conduction_enthalpy/circuit.rs` regressions use independent RC
and RL decay formulae, layered local heating, a nonlinear scalar temperature/
resistance balance, an external source with another heat-producing resistor,
source/charge continuation, and refusal/retry checks (G1/G3/G4). Its material-load
regressions additionally compare closed-form current-driven heating for both
source bases, an independent RC/thermal balance after a data-only material swap,
continuing charge/temperature feedback, whole-volume overlap and source-geometry
refusals, and late cancellation without partial publication. These tests were
authored and source-reviewed; Rust parsing/format checks passed. They were not
executed in the implementation environment because the workspace disk was full
and DSR/RCH were unavailable. No runtime conformance or empirical calibration
claim follows from their presence.



## Material collider density and pressure-dependent falloff

`MaterialRateScale::with_colliders` binds a declared third-body, Lindemann or Troe law
to an existing reaction channel. Every species in the original source catalog,
including inert gas, has an explicit reaction-card efficiency property; multiple
species may intentionally share one property. Claim identities for the base
kinetic law and all added coefficients are pinned before evolution. Each
chemical stage resolves the complete coefficient set at the actual parcel
condition point and forms `q = sum(epsilon_i n_i)/(V c0)` from the live gas
inventory. Neither a substitute pressure nor a second species store is accepted.

Third-body kinetics apply `q` to both directional rates. Lindemann uses the base
activity scale as `k_inf`, the explicitly normalized material property
`reaction_low_pressure_rate_scale` as `k0`, and the common factor
`Pr/(1+Pr)`, with `Pr = k0 q/k_inf`. Both scales use mol/(m3 s) because the
pressure activities and `q` are dimensionless; imported concentration-based
low-pressure constants require the explicit activity and `c0` conversion.
The existing `fs-thermochem` evaluator retains real species activities, reverse
equilibrium and the common multiplier. Exactly zero effective colliders give
zero rates; a positive rate lost to numerical range refuses coupled publication.

The shared channel evaluator carries this law into existing explicit,
implicit, mechanism, subcycled, piston, flow and electrical/thermal reaction
paths. Formation energy, mass, element conservation and accepted-state entropy
remain owned by the original conservative gas update. No extra combustion heat
is injected. Focused synthetic regressions cover inert-pressure changes at fixed
reactant partial pressure, independent low/high-pressure limits and reverse
balance, absent colliders, data-only efficiency changes through conservative
piston motion and acoustic state, implicit endpoint temperature/composition,
source-domain refusal and cancellation.

This is bounded ideal-gas third-body/Lindemann/Troe chemistry. It does not provide
SRI broadening, chemically activated reactions, empirical fuel mechanisms,
flames, shock/ignition prediction or validated engine fidelity. Formula
conventions follow the [Cantera reaction-rate reference](https://www.cantera.org/stable/reference/kinetics/reaction-rates.html)
and [falloff reference](https://www.cantera.org/stable/reference/kinetics/rate-constants.html).
Program anchor `frankensim-ext-gas-chemistry-ladder-paqh` remains open.

### Source-backed Troe broadening on the same reaction channel

`ColliderRegime::Troe3` and `Troe4` extend the existing collider owner with the
common factor `Pr/(1+Pr) * F(T, Pr)`. Both pin and resolve the dimensionless
`reaction_troe_alpha` and kelvin-valued `reaction_troe_temperature_3` and
`reaction_troe_temperature_1`. `Troe4` additionally requires
`reaction_troe_temperature_2`; `Troe3` deliberately omits that last term.
There are no hidden coefficient defaults or species-name special cases.

The center is `(1-alpha) exp(-T/T3) + alpha exp(-T/T1)`, with `exp(-T2/T)`
added only for `Troe4`. This bounded model requires `alpha` in `[0,1]`, positive
characteristic temperatures and `0 < Fcent <= 1` at the actual gas temperature.
Source coefficients retain their existing units, validity domains and pinned
claim identities during implicit iterations. Invalid selected coefficients
still refuse when effective colliders are exactly absent.

The implementation evaluates the center by log-sum-exp and combines the
Lindemann and Troe factors in logarithmic form. For the canonical coordinates
`x = log10(Pr) - 0.4 - 0.67 log10(Fcent)` and
`y = 0.75 - 1.27 log10(Fcent) - 0.14 x`, the broadening exponent uses
`y^2/(x^2+y^2)` with scaled ratios. This admits the continuous `y=0` value
without dividing by zero or overflowing either square. Positive final factors
lost to numerical range refuse; physical collider absence remains exact zero.
The specified finite Troe algebra, including its pressure tails, is preserved
without replacing the broadening factor by an assumed exact limiting value.

Focused synthetic regressions compare independent Troe equations across
pressure regimes, exercise source-only changes through implicit piston motion
and conservative energy, and cover coefficient/source/cancellation refusal.
Supplemental shared-law regressions cover concentration normalization, exact
absence before extreme normalization, native range statuses and precedence of
`ZeroActivity`. These tests are authored and source-reviewed; native Cargo
execution remains unverified in the resource-exhausted implementation workspace.
The original channel evaluator supplies all reaction consumers, thermodynamic
reverse closure and conservation; this extension creates no second rate or
combustion-heat owner.

## Directional electrical transport on the shared enthalpy timeline

`ConductionThermalNetwork` and its jointly accepted candidate now use the same
generic native field owner for scalar, oriented orthotropic and mixed regions.
Each nonlinear thermal trial refreshes every original principal electrical
claim at its actual temperature, retains the original basis and material/mesh
frame, solves current, and deposits the original per-cell `E.J` exactly once.
Voltage, equipotential total-current and distributed current-density drives
all retain their existing boundary and source-work accounting. Continuation
publishes the accepted thermal history and endpoint tensor operator together.
Source-domain, solve-budget and cancellation failures retain the old pair.
The fixed-density, fixed-geometry DC limits remain; an explicitly owned liquid
mass fraction also admits the confined phase-changing path below. Principal
axes do not rotate with temperature. Independent tests derive the
nonlinear heating equation for a rotated resistivity law, continue its accepted
state, verify transverse current, and exercise source and late-cancellation
rollback. Coefficients are synthetic, not a material-accuracy claim. Program
anchor: MR43 (`frankensim-material-reality-ui2if.7.1`), remains open.

Oriented planar conductor materials also pass through this native circuit
connection. Each thermal trial refreshes both principal source properties at
its endpoint temperature, retains the fixed material frame, and uses the same
field resistance to stamp the circuit branch. The original circuit midpoint
supplies the spatial voltage; actual per-triangle Joule losses enter the common
enthalpy residual. Scalar defaults, topology, fixed C/L, source work, storage
and heat routing remain shared. A rotated anisotropic RC fixture checks three
continuing discharge/heating steps against independently solved nonlinear
midpoint circuit and finite-capacity heat equations, including original pins,
transverse coordinate components and combined energy closure.

### Confined phase-changing spatial conductors and circuits

The same spatial/circuit/enthalpy owner admits solid, mushy and fully liquid
caloric states when every chart knot's density matches the represented fixed
volume within its existing volume tolerance. Non-solid charts require an
explicit liquid MASS fraction coordinate. The electrical material retains that
declared coordinate, and the thermal binding must drive the exact same axis
from its actual enthalpy state on every nonlinear trial and accepted endpoint.
This remains active on an isothermal latent plateau, where temperature alone
cannot describe changing transport.

Both scalar and oriented principal electrical laws use the existing source
resolver, original claims, electrical basis and fixed tensor frame. The supplied
material law determines phase-dependent conductivity/resistivity; no effective
mixture rule, phase-interface impedance or extra latent-heat source is invented.
Current redistribution, terminal resistance, capacitor/inductor evolution and
local Joule heat remain in their original native owners and one energy balance.
Moving geometry, shrinkage, fluid motion and density-changing phase charts
remain outside this stationary confined model.

Three synthetic regressions in `conduction_phase_enthalpy` independently check
current/voltage-driven latent evolution, continuation into liquid sensible heat,
capacitor storage and phase-dependent discharge, and source/axis/density/late
cancellation refusal. Focused parsing/formatting and source review passed;
native Cargo execution remains unverified. MR43 remains open for broader EM
and moving-material coupling.

## Source-backed updates of an already vibrating duct wall

Uniform and regional material-wall owners now expose `rebind_controlled` on
the continuing native DG state. New source-resolved density, stiffness and
viscosity, including supported temperature changes, retain radial displacement,
surface momentum, gas waves and every viscothermal memory coordinate. The native
transfer reports actual changed wall storage as signed parameter work once;
the external thermal/material owner owes that work and cannot also call it heat.
Source receipts, state, reciprocal operator and a newly certified propagator
publish together at the same physical clock after applicability-band, work,
state-identity and cancellation checks. No vibration restart or amplitude
normalization occurs. Failed source or propagator preparation retains the old
complete state. Changing gas/mesh, thermal expansion and mass transfer remain
outside this operation; cold regional hereditary construction is unchanged.
Focused tests continue an independent
sealed-wall oscillator, preserve regional diffusion history, change supported
source temperature, and exercise source/work/history/cancellation refusals.
Program anchors: MR08/MR34 (`frankensim-material-reality-ui2if.1.8`, `.5.4`),
remain open for broader changing-geometry and reduced-history behavior.

The same live operation now retains sourced Maxwell internal variables. Every
arm keeps its elastic extension `z_j/K_j`, hence its viscous offset relative to
the retained wall displacement. Its modulus/time property pair and exact
selected claim IDs establish identity, with the same immutable material card,
base claims, wall thickness, axial restraint and positive-arm support. Reordered
declarations and uniform/regional storage layouts are handled by complete
per-element correspondence. Unlabelled authored spectra, changed sources,
activated/extinguished arms and incomplete maps refuse atomically.

Changed Maxwell spring storage joins the existing parameter-work ledger once;
a relaxation-time-only update has zero instantaneous work and changes future
loss in the native DG propagator. Five focused regressions cover an independent
analytic gas/wall/Maxwell trajectory, physical storage quadrature, source-driven
temperature updates and rollback with nonzero history. Parsing/formatting and
bounded source review passed; native test execution remains unverified. This
does not solve heat, infer branch lineage across different materials, or turn
instantaneous parameter work into dissipated heat.

## Material-driven biased magnetostrictive coil/core

`magnetostrictive` composes the original source-bound incremental magnetic
material with the existing `fs-phs` descriptor and port interconnections. One
fixed-root rod uses axial displacement `u(z)=x*z/L`, reference volume `A*L`,
and consistent mass `rho*A*L/3`. Its ideal N-turn coil uniformly links the core:
`lambda=N*A*Delta B`, `Delta i=L*Delta H/N`. The native state is
`[x, momentum, Delta lambda]`. Differentiating the same integrated material
energy supplies axial force and current; no independently fitted sensing or
actuating gains, material-name dispatch, or additional time integrator exists.

The complete material, typed nonzero DC bias and reference geometry remain
fixed. The reference shape is already equilibrated at that bias. Native ports
are force/velocity and incremental voltage/current, and the reported account is
reversible perturbation energy plus kinetic energy. It excludes the DC-bias
source's preparation/absolute energy. This incremental passivity/accounting
boundary must not be relabeled a closed energy account for an entire biased
device. A coupled caller may connect the native component to electrical and
mechanical loads while retaining that same perturbation convention.

The coil reduction gives `L_c=N^2*A*mu_axial/L`, `theta=N*A*(d/s_H)/L`,
`lambda=L_c*i+theta*x`, and `F=k_H*x-theta*i`. A held current gives
`k_H=A/(L*s_H)`; a held flux gives `k_B=A*mu_T/[L*(s_H*mu_T-d^2)]`.
Zero incremental voltage conserves flux and is not the zero-current boundary.
Ideal current clamps require compatible accepted flux; no switching impulse
or conserved-state projection is introduced. Candidate source, amplitude,
native solve, work, or late cancellation refusal leaves the body unchanged.

The caller's strain/H increment envelope is checked through the material
storage and at acceptance; its H limit must be below the absolute bias.
The ideal fixed magnetic path equals the rod length. Distributed waves,
off-axis deformation, magnetic fringing/demagnetizing fields, Maxwell forces
from changing magnetic geometry, coil resistance, eddy/hysteresis losses,
thermal evolution and acoustic radiation are not supplied by this component.
It computes axial vibration and sensing; it does not synthesize a hum signal
or claim microphone pressure or empirical transducer accuracy.

Synthetic `tests/magnetostrictive.rs` regressions independently predict the
fixed-current/fixed-flux midpoint oscillators, source-sign/density/turn-count
effects, current- and voltage-driven motion/work, and both coupled normal
modes of a native coil/capacitor interconnection. They also check incompatible
clamps, amplitude/domain refusal and late-cancellation retry. This advances
MR60 (`frankensim-material-reality-ui2if.7.8`); the broader thermal/history,
spatial field and acoustic-observation requirements remain open.

## Time-dependent plastic flow on the common solid/enthalpy owner

`ThermoplasticSolid::with_linear_overstress` selects the canonical positive
viscosity from each element's original mechanical/caloric material before
loading. Selection is in existing connectivity order; `None` retains the
rate-independent law. No second mesh, history, constitutive solver, caloric
inventory or integration timeline is introduced.

The existing shared mechanical candidate passes its actual physical dt to the
native J2 stress, tangent, history and work account. This same endpoint is used
by quasistatic loading, consistent inertia, solved-pressure thermal contacts
and electrical/contact heating, including their inertial wrappers. Each thermal
or Newton trial restarts from the accepted displacement and plastic history;
thermal substeps do not advance the constitutive clock. Persistent positive
overstress is admitted on continuation. A held load now causes finite-time
creep rather than an immediate rate-independent return.

Yield plus viscous dissipation enters the original finite caloric cells exactly
once. Viscosity, yield and the other selected coefficients are refreshed from
the same complete source at actual endpoint temperature, so material data and
the generated heat alter subsequent deformation. Hardening storage, reversible
heat, constitutive numerical loss, kinetic storage and inertial numerical loss
keep their existing distinct accounts. The common source/heat/work gates still
publish all candidates together or leave accepted state and clock untouched.

Synthetic `tests/solid_plastic_enthalpy/overstress.rs` independently solves the
uniaxial Perzyna and finite-capacity heat equations over successive actual
loaded-tetrahedron steps in tension and compression. It checks temperature-
dependent viscosity substitution, material/time scaling, thermal-substep
independence, shared inertial response and losses, missing/unsupported sources,
duration-free refusal and failed-proposal retry. This is the existing fixed-
reference, single-solid, small-strain/rotation model with linear overstress;
it does not claim finite-strain creep, phase change, damage, general rate
exponents or measured material accuracy. MR22
(`frankensim-material-reality-ui2if.3.4`) remains open for broader coverage.

## Live gas propagation through source-bound compliant ducts

`MaterialWallDg::rebind_gas_field_controlled` and its regional counterpart
compose the original native uniform gas intervention with the retained
`AxialGasField`. Actual compression and gas momentum survive; pressure,
volume flow, native propagation, matched boundaries, pressure ports and
reciprocal wall loading use the same new medium. Material data determines the
change, without a pitch multiplier or replacement oscillator.

Uniform/regional wall motion and local Maxwell histories retain their
physical values. The original signed parameter-work ledger records the
actual acoustic storage change once. The adapter prepares the new source
field, state, operator and same-timestep/tolerance native propagator before
publishing any of them. Subsequent existing wall-material updates remain
bound to the refreshed gas.

The operation is a prescribed temporal material intervention at fixed
geometry. The enclosing gas/thermal owner must separately account for
background heat, species transfer and the stated parameter work. It adds
neither mean flow nor static wall pressure traction. Nonuniform gas and
native viscothermal gas closures refuse until their state transfer is
represented. Source erasure, foreign owners, excessive work or late
cancellation preserve the prior accepted trajectory and propagator.

Five native/integration regressions cover independent wave continuation,
regional Maxwell motion, port and cold-assembly equivalence, actual
source-data gas changes followed by wall rebinding, and atomic refusal.
This advances MR08/MR33/MR34 while their wider spatial/history requirements
remain open. Native execution and measured acoustic validation remain
separate gates.

## Live gas changes through the shared acoustic port

`DgPortLoad::rebind_gas_field_controlled` extends the original material-gas
transfer to a continuing bore with finite branches and causal exterior loads.
`DgGasRefresh` selects the existing uniform lossless, uniform viscothermal or
spatial viscothermal transfer. Source-bound uniform/regional walls and their
finite gas evolution use the same preparation path. No replacement storage
law or memory fit is introduced.

Compression, momentum, wall motion and admitted diffusion/Maxwell histories
survive. Pressure, volume flow and the exact mean inlet impedance follow the
new gas. The original native storage difference enters parameter work once;
the enclosing background owner owes its heat/species balance and that signed
work. A refresh consumes no sample and preserves accepted output frames.

Independently sourced branches and exterior loads retain their gas and causal
history. Their junction and terminal eliminations are rebuilt against the new
bore. Common static-pressure and computed-opening admissions remain active.
Attached planar and circular channels refresh their actual local gas through
the original native momentum transfer. Each resolved cell retains momentum,
new source viscosity drives its next implicit response, and the retained Mach
ceiling applies to the remapped field. Independent branch channels retain
their original gas. Rigid closures update their source without fluid storage.
All source, state, propagator and algebraic coupling changes publish together
after work, source and cancellation checks.

`DgPortGasUpdate` retains the original main-bore report and each identified
channel's native kinetic transfer, plus aggregate represented storage and
signed parameter work. The work allowance bounds the sum of absolute native
interventions, including every channel; opposite transfers cannot cancel an
excess. Each channel's separate cumulative parameter-work account enters the
aggregate ledger once and is subtracted from its storage/work residual.
Pressure work, wall heat, numerical dissipation and rational exterior exchange
remain distinct. Historical channel frames survive a zero-duration refresh;
subsequent steps start from current kinetic storage.

This prescribed fixed-geometry intervention advances MR33/MR67 and the MR05
atomic-state requirement. It does not solve mean flow, heat/species transport,
static pressure traction or an exterior gas change. Independent-source
branch/terminal histories and late refusal are covered by authored focused
regressions; native execution remains a separate verification gate.

`ApertureDg::rebind_gas_field_controlled` publishes that continuing load with
the original scalar, independent-pair or material-modal source. The source's
density and exact mean impedance follow the same inlet gas; source-bound slit
viscosity follows its complete gas state. Explicit numerical viscosity remains
explicit. Structural propagators, positions, velocities, signed pressure
projections, contact and solid relaxation histories are retained. The
quasisteady aperture carries no gas inertia to reset or charge as storage.

For a valve between two continuing ducts,
`rebind_series_gas_fields_controlled` takes explicit upstream/downstream fields
on their own geometries and prepares both before publishing either. The
complete inlet gas must agree, and the valve receives the sum of the new exact
mean impedances. The work allowance bounds all absolute interventions in both
ducts and their finite channels; opposite signed transfers cannot hide an
exceeded budget. Independent pressure,
viscous/Bernoulli flow, continuum storage, two-blade oscillator and exact-retry
regressions cover source changes after nonzero mechanical/acoustic histories.
This does not add a density-changing mass-flow junction or a new gas-energy law.

Focused regressions also exercise a moving aperture, both finite channel
geometries and retained causal terminal history in one continuing source-gas
transaction, with independent cell momentum/kinetic work and atomic retry.
The exact pinned compiler checked the complete production graph and the
`channel_gas` and `dg_gas_rebind` integration targets successfully. Runtime
test execution remains unverified: dependency code generation was interrupted
under sustained shared memory pressure before a test binary ran. This is not
complete constellation qualification or measured acoustic validation.

## Finite thermochemical feedback through the complete acoustic load

`DgPortLoad::step_with_finite_gas_controlled` connects the original finite
thermochemical inventory to the complete continuing bore, finite terminal
passages, branch coupling and causal exterior load. One original held-flow
sample is followed by the existing scalar work solve. Each trial starts from
the same original gas inventory and the same marched acoustic histories;
the native gas transfer supplies the actual bore AND channel parameter work.
The existing heat convenience uses the original caloric/transport source
owner, with no replacement gas equation or integration scheme.

The inventory's physical volume must equal the uniform circular main bore
plus all of its main-gas slit/circular terminal passages. Independently
sourced branches and their passages remain separate inventories. Their
original equilibrium-pressure requirements, computed exterior source gates,
channel Mach limits and sum-of-absolute-work allowance remain active.

`FiniteGasDgPortStep` retains the original background report, gas energy change,
heat, separately transported enthalpy, signed work and actual native transfer.
The whole endpoint balance includes represented bore/branch/channel storage,
original loss and net acoustic port work. Rational terminal and driven-mouth
exchange remain signed boundary work; no rational storage is invented.
Historical interval frames and impulses survive the material refresh. Returned
endpoint pressures and the full current storage/parameter-work ledger use the
accepted gas. All source, work, volume and cancellation checks precede joint
gas/load publication.

This is a bounded first-order frozen-medium partition at fixed geometry.
Acoustic volume flow does not imply thermodynamic mass transport, and original
wall/viscous/numerical loss is not silently deposited into the background.
The callback owns any actual heat/species process and returns uncommitted
companions for publication only after success. Focused synthetic regressions
cover independent caloric/compressive-energy references, native nozzle/species
transfer with per-cell momentum and kinetic work in both channel geometries,
and refusal/cancellation/exact retry. Native execution and measured physical
validation remain separate gates. This advances MR05/MR06/MR33/MR67.

## Finite thermochemical feedback through the moving aperture

`ApertureDg::step_with_finite_gas_controlled` composes the same complete-load
transaction with the original scalar aperture, or `LoadedAperture` single,
paired and material-modal mechanics. The source and acoustic load advance
once in the initial medium. The finite inventory then pays the actual bore
and channel parameter work, and the original source consumes the accepted
density, viscosity and exact load impedance without losing mechanical or
material memory. Independent blade forces and modal pressure work remain in
their original owners. The heat convenience shares the existing caloric owner.

`FiniteGasApertureDgStep` reports one accepted interval. Its `gas.acoustic`
and `coupled.duct` are the same acoustic account and must not be summed.
Held-force/flow integrals retain the medium in which they were advanced;
endpoint pressure and storage use the final gas. `DgApertureFrame` exposes
the signed `parameter_work_j` separately, and its balance subtracts that
intervention. The combined gas/mechanical/acoustic residual includes actual
drives, original losses, signed exterior exchange, heat and transported
enthalpy exactly once. All fallible checks and the final cancellation poll
precede joint publication of gas, mechanics and the complete load.

An attached upstream tract is refused before any trial: two independent
finite gas inventories require a joint solve. The first-order fixed-volume
partition does not infer mean-pressure prestress, thermodynamic mass flow
from acoustic flow, or heating from mechanical/acoustic dissipation.
Focused synthetic regressions cover scalar caloric/wave-work feedback,
independent paired-blade motion, exact retry after source refusal or final
cancellation, and rejection of an unowned second inventory. Native execution
and measured physical validation remain separate gates.

## One finite thermal body on the moving two-gas valve

`DoubleActingGasPiston::thermal_valve_candidate` composes the metering valve
with the existing two-face finite-enthalpy thermal owner. Native compressible
species/enthalpy transport is sampled once at the original opening. Every
bounded heat/momentum trial then starts from that same post-transfer pair and
the same original body; no trial resamples mass flow or spends received gas.

The original thermal solver resolves both conductances at the actual two gas
midpoints and one body enthalpy midpoint. It retains the source evaluations,
separate heat brackets and aggregate joint-trial budget. Sensible or latent body
storage uses the existing material phase chart. The accepted face heats enter
the original valve finalization, including its post-motion pressure-order gate
and independent native caloric, kinetic and force-work measurements.

`DoubleThermalValveCandidate` returns both gases, shared wall and the one body
endpoint together. Its report contains the native mechanical account once,
inside the valve report, plus the original thermal/interface observations.
The enclosing energy residual adds the body's independently measured enthalpy
balance to the reference-stable valve balance: both face heats and transported
enthalpy cancel internally, while all prescribed heat and force work stay
explicit. Source failure, nested-work exhaustion or final cancellation returns
no partially advanced owner. Callers publish the complete candidate without
applying its accepted heat a second time.

The metering geometry remains frozen over the first-order transport stage.
Finite-body negligible-pV-work and fixed-geometry assumptions remain those of
the original enthalpy owner; no wall melting-induced shape, jet recoil or
heterogeneous acoustic junction is inferred.

Three focused regression groups cover independent caloric/species/body energy
with sensible and latent storage, exact reductions to the original insulated
valve and closed-valve thermal owners, and late source/budget/cancellation
refusals with exact retry and continuation. They are authored but unexecuted in
the unavailable pinned native build environment; syntax/format checks do not
establish numerical or experimental validation.

## Moving metering valves between independent finite gases

`cavity_phs::piston::double_acting::valve` connects the original two finite
gas inventories through a rectangular metering window on their existing shared
wall coordinate. `RectangularMeteringValve` derives edge displacement from
the primary volume and its original face area; the exposed area is the exact
window overlap. It creates no second position state and does not clamp motion.

`DoubleActingGasPiston::valve_candidate` and `step_with_valve` select the
donor from actual initial pressures, call the original compressible nozzle and
species/total-enthalpy transfer once, then call the original two-gas midpoint
pressure-work and single-wall momentum solve. Different gas compositions,
source caloric laws, temperatures and densities retain their original owners.
Native throat/choking, source-domain, gross donor-fraction and transfer
conservation gates remain active. Final pressure ordering is checked again
after the actual heat and motion.

The endpoint audit independently re-evaluates post-transfer caloric increments,
actual wall kinetic change and force work. The original transfer's conserved
species datums cancel before accumulation, preserving small thermal signals.
Transferred enthalpy is internal to the pair; jet kinetic energy and pressure
flow work are not deposited again. A late source, geometry, pressure, balance
or cancellation refusal leaves both gases and the shared wall unchanged.

This remains a first-order frozen-opening, transfer-then-motion partition.
It supplies no jet recoil, resolved valve-fluid inertia, metering-edge friction
or heterogeneous acoustic junction law. Actual mechanical/heat port observations
are retained; for positive mass transfer, the full typed stream/mixing-entropy
audit is explicitly unavailable because the native pair report does not supply
that entropy account. Closed flow preserves the original native port accounts.

Four focused regression groups cover independent species/caloric balances,
both donor directions, geometry-driven opening/closure, analytical constant-
caloric choking, changed shared-wall response, and late refusal/cancellation
with exact retry. These tests are authored but not natively executed here;
formatting/parser checks are separate from numerical or physical validation.

## Finite gas feedback retains spatial and molecular histories

`FiniteGasDgLoss::SpatialViscothermal` and
`FiniteGasDgLoss::UniformMolecularDiffusion` connect the existing native memory
transfers to finite-gas stepping. One conversion selects the same native owner
for source-wall loads, complete port loads, transport networks, and one- or
two-tract moving apertures. Spatial quadrature and radial sections remain bound
to their original geometry. Molecular transfers retain source laws, species
fractions and reservoir support while preserving physical entropy perturbations,
compression and momentum; their actual storage work is paid by the finite gas.

The background remains uniform and the original geometry-volume admission
remains active. Molecular composition changes refuse while the original
acoustic species histories are retained. Same-composition inventory transport
and heating use the existing caloric and species/enthalpy owners. These choices
do not create a memory realization or discard an existing one.

After native source, memory and work admission, `prepare_dg_gas` retains the
certified propagator when the complete native operator identity is unchanged.
That identity includes physical mass transforms, geometry and typed boundaries.
Changed operators rebuild at the original clock and tolerance.

Focused public regressions cover spatial and molecular finite heating against
the original native transfer, unchanged-medium continuation, same-composition
transport, source/control refusals and final cancellation with exact retry.
They are authored but not executed in this environment; newer-toolchain
formatting/parser checks do not establish native compilation or numerical proof.

## One finite gas inventory through the two-sided moving aperture

`LoadedApertureDg::step_with_shared_finite_gas_controlled` advances the original
two-sided moving source and both complete acoustic tracts against one explicitly
shared, well-mixed `ClosedGasParcel`. Its volume must equal both main bores plus
their main-gas finite terminal passages, and both original complete inlet states
must match that parcel. Independently sourced side branches retain their own
inventories. The heat convenience method uses the existing caloric gas owner.

The original mechanics and both acoustic tracts advance once per accepted sample.
A shared scalar feedback solve then pays the sum of both native gas-parameter
work accounts with one debit from the original parcel. Each trial starts from
the same original gas and the same marched acoustic histories; it does not
repeat the mechanical step. The same native admission, background energy,
memory transfer and endpoint-report helpers serve the existing single-tract
finite-gas path and this two-tract path.

Both refreshed loads retain their bore and finite-channel momentum, terminal
pole histories and interval force/flow/work observations. One absolute-work
allowance covers both sides. Endpoint storage, gas density/viscosity and the
source's exact summed load impedance follow the accepted shared gas.
`FiniteGasSeriesApertureDgStep` reports the ordered upstream/downstream native
refreshes, one background result, one gas-work debit, and the enclosing
gas/mechanics/acoustics/loss/exterior-work residual.

Admission, source evolution, native transfers, whole-system energy closure and
final cancellation all precede publication. A refusal leaves the gas, original
moving blades and both complete tracts retryable. This remains a fixed-volume,
first-order partition with one common gas background; it does not infer mean
mass transport from acoustic volume flow or admit independently evolving inlet
densities through the existing equal-gas valve.

Three public regression groups in `tests/dg_gas_rebind/finite_series.rs` cover
an independent caloric/two-bore-work reference and complete energy account,
actual summed-volume/source admission, and late refusal/cancellation with exact
retry. These tests are authored but have not been executed in this environment:
the pinned native compiler/build route and required sibling Git checkouts are
unavailable. Formatting/parser and whitespace checks are not native test results.

## Live regional gas fields through branches and moving apertures

`DgPortLoad::rebind_regional_gas_fields_controlled` changes the optional main
field and explicitly identified finite branches in one transaction. Branch
inputs use their original external port IDs; unknown and duplicate IDs refuse
before native work. Reports follow original attachment order. The main-only
entry point delegates to this same owner, and unselected native propagators
are retained instead of rebuilt.

Each selected duct uses the existing uniform lossless, uniform diffusion or
spatial diffusion transfer, preserving its original compression, momentum,
wall and admitted memory coordinates. Junction and computed-exterior pressure
gates inspect the complete proposed endpoint. Coordinated regional pressure
changes can therefore be admitted without publishing an inconsistent
intermediate junction. Finite slit/circular terminal passages follow the
accepted gas at their actual branch or main outlet and retain cell momentum,
history and original Mach limits. Branch and terminal elimination are rebuilt
once from the accepted regional sources; exterior pole states do not reset.

`DgRegionalGasUpdate` identifies the optional main transfer, selected branch
transfers and every native terminal-channel transfer. One allowance bounds
their sum of absolute parameter work; the signed sum equals the represented
storage change within native arithmetic. Subsequent aggregate step ledgers
include each branch's cumulative parameter work once, separately from pressure
work, dissipation and rational exterior exchange. Source, state, propagator
and coupling changes publish only after all checks and final cancellation.

Scalar, independent-pair and material-modal `ApertureDg` sources consume this
same regional transaction. A branch-only change updates the source's exact
mean impedance while retaining inlet gas, viscosity and mechanical history.
The series regional entry point prepares both tracts, checks their common
complete inlet gas and summed impedance, and applies one absolute-work
allowance across both. An explicitly unselected side retains its state and
reports `None`. Existing main-only APIs share these candidate helpers.

This is a prescribed, fixed-geometry regional material intervention. The
background owner owes its thermal/species account and signed work; it does not
add a density-changing valve mass-flow law, static pressure force, mean
transport or exterior-medium update. Six focused synthetic regression groups
cover independent branch/channel momentum and energy, native diffusion memory,
coordinated pressure, scalar and paired source response, series coupling and
transactional refusal/exact retry. Native compilation, execution and measured
physical validation remain separate gates. This advances MR05/MR06/MR33/MR67.

## Finite transport networks drive complete acoustic loads

`dg_port::finite_gas::transport` connects the original multi-parcel material
transport transaction to `BTreeMap<u64, DgPortLoad>`. The original wall-DG and
complete-load entry points share one private, statically dispatched transaction
kernel. No second transport, caloric, work-feedback or conservation algorithm
is introduced. Existing wall-node structs and their entry points remain usable.

Nozzle, diffusion, porous and membrane paths consume one unchanged original
species/enthalpy snapshot. Each complete load then advances once from its
original medium and pays its actual bore/channel parameter work from the
transported inventory. Finite reservoirs share the same source map and caloric
owner. Final driving-force gates run after every endpoint's heat and acoustic
work. Compensated sums independently close total mass and enclosing gas plus
bore/branch/channel energy; native losses and net acoustic port work retain
their original signs. Exterior exchange is already in net port work and is
not charged again.

Complete source, actual bore-plus-passage volume, controls, clock and the one
original inlet flow are admitted before transport. Explicit candidates expose
all proposed parcels, complete loads and original reports. The advance APIs
publish both maps only after all endpoint, conservation and cancellation checks.
Failed attempts retain every accepted source, channel momentum and causal pole.

The gas volumes remain separate fixed-volume, well-mixed owners connected by
native material transport. A restriction does not thereby gain a resolved
acoustic scattering model, jet noise or static pressure force. Independent
branch gases remain their original sources, and native composition/pressure,
transport-domain and first-order partition limits remain active. Focused
synthetic regressions cover actual transport into complete warm loads and
transactional refusal/retry. Native execution and measured validation remain
separate gates. This advances MR05/MR06/MR33/MR67.

## Continuing spatial dielectric charge through thermal phase change

`ConductionThermalNetwork::step_dielectric_candidate` connects the native
`PlanarLeakyDielectric` history to the original cell enthalpy residual. Every
thermal trial resolves the same complete material source bundle at actual
temperature and liquid mass fraction, refreshes scalar or oriented conductivity,
and advances from the same accepted electric charge and potentials. Local
midpoint E.J power is routed once through the existing original-triangle to
receiving-cell correspondence. Electric storage change remains a separate
physical ledger term; it is not deposited as additional heat.

The jointly accepted candidate owns the continued native field, finite thermal
state and original combined energy report. Source work includes displacement
current during electrode-voltage ramps. Thermal conduction, prescribed heat,
Joule deposition and dielectric storage close in the same acceptance path.
Confined equal-density phase change follows the existing enthalpy admission;
changed geometry or phase density is refused. Permittivity must be explicitly
constant over its source validity domain and unchanged along the trajectory,
so variable dielectric parameter work is not silently discarded.

Thermal rates and conductivity use the endpoint; overall coupled temporal
accuracy is first order in general even though the held-coefficient electrical
interval uses midpoint. The original material identities, frames, source
receipts, clocks and trace topology must agree. Cancellation, missing sources,
changed storage or late thermal refusal leave all accepted owners unchanged.
Dispersion, induction, breakdown, charge-carrier chemistry and moving interfaces
remain outside this coupling.

Four regression cases in `tests/conduction_phase_enthalpy/dielectric.rs`
cover independent epsilon/sigma relaxation and local heat, heterogeneous
interface charge and driven electrode work, latent melting with phase-dependent
conductivity, and source/storage/cancellation rollback. Scalar and rotated
tensor conductivity both use the native path. These cases are authored and
source/parser checked, but were not executed in this environment.

## Joint endpoint chemistry, finite gas heat and spatial solid enthalpy

`GasReaction::heated_implicit_candidate` solves the original implicit mechanism
at the actual heated or cooled endpoint `U_end = U_start + Q`, with unchanged
closed-parcel mass and volume. Every existing channel rate, material source,
activity, collider and Troe evaluation therefore sees that same endpoint
temperature, pressure and composition. The original stoichiometry, consumption,
source-support, caloric inversion and work limits remain authoritative. Chemical
formation energies remain inside gas storage; no additional reaction-heat source
is introduced. Zero supplied heat preserves the existing adiabatic path.

`SolidEnthalpy::step_with_coupled_reacting_isochoric_gas` uses that capability
inside the original spatial cell/gas thermal residual. The unknowns remain cell
energy increments and total gas energy increment. For each outer thermal trial,
the original inner mechanism solves its extents at that proposed gas energy;
the resulting gas state determines the original film transport and thermal
residual. Wall cooling consequently changes the reaction extent within the
same accepted interval. This is nonlinear elimination of the existing chemical
equations, with no intermediate adiabatic-reaction-then-cooling state.

The common gas endpoint callback carries the identical accepted parcel, native
heat receipt and chemical report through the original final thermal acceptance.
Ordinary nonreacting calls use the same path with their existing heat-only
callback. Finite wall heat cancels between participants once, authored external
heat remains supply, nonlinear spatial k(T) uses the original prepared material
network, and actual solid endpoints re-admit original expansion/elasticity
sources before mechanical recovery. A single gas inversion policy must be
shared by chemistry and heat. No candidate changes accepted gas, thermal or
mechanical state on numerical, source, budget or late cancellation refusal.

The heat-aware reaction report distinguishes signed endpoint heat-entropy
transport `Q/T_end` from `delta S - Q/T_end`. Cooling may lower gas entropy;
the latter residual is admitted nonnegative within the original tolerance.
It includes backward-Euler numerical entropy and is neither isolated chemical
entropy production nor the entropy of an unspecified external reservoir.

This path requires `with_implicit_mechanism`, including for one channel;
chemistry-only subcycling is refused because the complete coupled interval must
be refined together. Inner and outer residuals have their original independent
caps. Nested work can be expensive, and no uniqueness or unconditional
convergence claim is made. Temporal accuracy remains first order. The gas is
homogeneous and fixed-volume, solid thermal geometry is fixed, and mechanical
recovery retains the original sequential thermoelastic interpretation. Flame
fronts, shocks, pressure traction, reciprocal thermoelastic heat and moving
interfaces remain outside this entry point.

Four regression cases in `tests/solid_reactive_gas/coupled.rs` cover an
independently derived heated/cooled isomer endpoint and zero-heat parity, an
independent joint reaction/wall-heat reference with native displacement and
stress, continued species/source history with nonlinear spatial conduction,
and source/control/budget/cancellation rollback with deterministic retry.
Changing only wall-conductivity data changes same-interval reaction extent.
The cases and runtime were parser/source reviewed, including independent
numerical checks of the analytic references, but native tests were not executed
because the pinned sibling dependencies and configured remote runners are
unavailable in this environment.

## Source-defined reactions drive the elastic gas boundary

`GasValveMotion::reactive_step_candidate` and `step_reactive` connect the original
`GasReaction` to the existing mesh/material-derived modal pressure boundary.
One native modal sample supplies the chemical duration. Source-backed explicit,
implicit, multichannel and subcycled kinetics remain owned by `GasReaction`.

The accepted owner is cloned with its original reference volume, reference modal
coordinates, equilibrium, material identity, modal state and observation transfers.
A fixed-volume chemical proposal conserves total gas internal energy, including
formation references. Every subsequent gas/motion trial starts from that same
reacted gas and accepted mechanical history; chemical substeps never advance the
structural clock. Native gas-loaded time-resolution gates are retained. A complete
candidate publishes once, or both histories remain unchanged.

The split is first-order chemistry then joint pressure/work/motion. Compression
and external heat affect the next chemical sample. The enclosing account starts
before chemistry: gas internal-energy change plus modal storage change, physical
modal loss and exterior pressure work equals supplied external heat within the
declared tolerance. Formation energy and modal loss are not added as gas heating.
Original chemical and post-reaction mechanical reports remain separate.

G1/G3/G4 regressions in `tests/solid_enthalpy/valve/gas_motion/reactive.rs` use the
existing six-tet cube's independent one-DOF mass/stiffness integrals and a
constant-cv energy-midpoint quadratic with the changed formation offset. They
cover source-rate/modulus substitution, forward/reverse excitation, inert
equality with the native step, continued anchored volume, physical damping
balance and rejected chemical/late-mechanical proposals with deterministic retry.
Coefficients are synthetic. Native regression execution is still required; only
source/API review and independent arithmetic were available in this session.

This remains a spatially uniform closed gas and fixed structural
material/geometry epoch. It adds no flame/shock propagation, chemical acoustic
relaxation, wall thermal transfer, nozzle flow, moving geometry/basis or automatic
exterior-radiation refresh. The wider gas-chemistry ladder remains open.


### Source-bound finite-strain dynamics (MR58 bounded increment)

`finite_solid::FiniteSolidMaterial::resolve` binds an explicit model index and reference
density from the same immutable `MaterialCard` at the exact
`ResolvedHyperelastic` context, retaining the original model, typed source
query, material identity and density receipt. The existing fixed-external
source admission rejects unsupported stress, strain, rate and history axes;
a numerically matching query alone cannot freeze them. Held temperature,
pressure and composition describe reference preparation/environment
conditions. This path does not refresh them during a trajectory.

`FiniteSolidBody` uses `hyper3`'s total-Lagrangian P1 internal force, exact
consistent tangent and strain energy, `linear3`'s consistent reference-mass
kernel, `fs-time`'s bounded `second_order_acceleration` initialization and its
`OperatorGeneralizedAlpha` at `rho_inf=1`. Uniform and per-element source fields
use the same operators. Reference mass is immutable during finite motion;
material density is never multiplied by deformed volume. Every free DOF
requires positive finite mass, including refusal of unused free nodes.
Prescribed displacements are stationary and attained initially.

Candidate intervals use held referential nodal/body loads and return actual
kinetic and strain storage, force-displacement source work, midpoint support
reactions and the signed energy integration defect. Nonlinear midpoint is not
exactly energy conserving: that signed defect must pass the requested
absolute gate and is never relabeled as heat. The gate is not a trajectory
accuracy certificate. No accepted q/v/a/history changes until the caller
publishes the returned candidate; domain, solver, workspace, cancellation and
late energy refusals leave the original body intact. Immutable matrices are
shared between candidates; the accepted-step budget bounds retained history.
Native geometry and solver errors retain their original owners; material
binding preserves the original source refusal text.

`FiniteSolidBody::new_with_pressure` additionally binds an immutable
`hyper3::pressure::ClosedTetPressureSurface` to exactly the original coordinates
and tetrahedra. Its entire exterior must be closed and consistently oriented.
Positive signed pressure acts inward on the current triangle areas. The existing
initial-acceleration solve includes this load, and each
`step_candidate_with_pressure` explicitly supplies the uniform pressure held
over that interval. The original dead-load step refuses a pressure-bound body;
a body initialized without pressure cannot acquire a binding after its initial
acceleration was prepared. Zero pressure skips the follower arithmetic and
recovers the original native state and replay exactly.

The adapter subtracts the native pressure force and analytic force Jacobian
from the material residual and tangent. Reference mass and the native midpoint
integrator remain unchanged. Pressure work is the endpoint expression
`pressure_work_j = -p (V_after - V_before)` and is included in `source_work_j`.
It is not replaced by midpoint force dotted with displacement: P1 volume is
cubic in position, so their difference belongs to the signed integration
defect. Stationary support reactions include the same current-area pressure.
Surface identity, finite pressure, current element orientation, source-domain,
solver, energy and cancellation refusals leave the accepted state unchanged.

The focused `finite_solid` integration tests specify an independent scalar
finite-stretch oscillator and its physical work/storage, density and
stiffness substitutions, rigid rotation with translating acceleration,
heterogeneous reference masses, source-axis refusal, inversion/stretch
limits, massless free-node refusal and late energy/cancellation rollback.
`tests/finite_solid/pressure.rs` adds a three-axis breathing fixture whose
independent scalar midpoint equation includes the analytic volumetric material
force and current pressure area. It checks initialized acceleration, changing
held pressures, actual motion, exact volume work, support inertia/reactions,
fixed mass, zero-pressure parity, deterministic replay and atomic refusal.
`examples/pressurized_material_solid.rs` exposes the same bounded motion using an
explicit normalized material card/model/temperature, or a labelled synthetic
`--demo`, and reports the computed volume, motion, storage and pressure work.
These are numerical fixtures, not experimental validation or measured material
calibration. The pressure production sources passed native parsing; full native
compilation and test execution remain pending because the shared executor
killed the real `asupersync` prerequisite even during a metadata-only check.

MR58 remains open. This increment is compressible, isothermal,
rate-independent hyperelasticity on conforming P1 tetrahedra. Evolving thermal
reference state, finite-strain inelastic evolution, moving supports,
contact, selected open pressure patches, nonuniform pressure fields, a coupled
fluid equation of state and a mixed/stabilized nearly-incompressible formulation
remain outside this path. No continuous-time inversion or
self-intersection guarantee, mesh-convergence certificate, or new adjoint
claim is made.


### Material piezoelectric circuit and finite load heating (MR59 increment)

`material_circuit::piezoelectric::PiezoelectricCircuit` places the original
longitudinal piezoelectric rod inside the native Kirchhoff descriptor. The
selected capacitor branch defines the electrode order. Its clamped capacitance
comes from the rod material; its charge-only energy stamp is removed before
inserting the original charge/strain energy and consistent inertia. Other
capacitors, inductors, source ports, transformers and resistor bindings retain
their native operators. Initial and endpoint algebraic admission use the true
coupled voltage and never project charge, flux, extension or momentum to repair
an incompatible ideal clamp.

One candidate solves rod motion and circuit current together. Existing
material-circuit accounting evaluates resistor heat at the actual algebraic
midpoint. Electrical port work and held force times extension change close
against combined storage and heat. Piezoelectric field, strain and kinetic
storage are counted exactly once.

`circuit_enthalpy::piezoelectric::step_candidate` uses the existing thermal-driver
boundary and original backward-Euler/Newton enthalpy solve. Load resistance and
conduction are resolved at the proposed endpoint temperatures. Each trial
restarts from the accepted electrical/mechanical prefix; each resistor has one
finite thermal receiver. Any source, domain, budget, work or cancellation
refusal returns no joint endpoint.

Five authored numerical fixtures cover independently eliminated loaded-rod
equations, source density/coupling/resistance and polarity substitutions,
additional C/L storage and Kirchhoff exchange, open/short limits, finite-load
thermal feedback, and late refusal/retry. Independent arithmetic is not native
Rust validation: parsing, compilation and tests remain unrun because the
executor is unavailable. These fixtures are not measured material calibration.

MR59 remains open. The existing uniform, small-signal longitudinal rod and
bounded dense circuit domains apply. Ceramic temperature, thermal strain,
dielectric loss, depoling, nonlinear/anisotropic 3-D fields and electrical
contacts are not evolved. The rod context and reference mass stay fixed;
explicit resistor/thermal-network states evolve. Coupled thermal accuracy is
generally first order; work closure is not a time-accuracy certificate.


### Joint implicit chemistry and native elastic gas-boundary motion

`cavity_phs::modal::reactive::coupled` adds
`GasValveMotion::{coupled_reactive_step_candidate, step_coupled_reactive}`.
It closes the original implicit mechanism, material gas energy inversion and
mesh-derived elastic wall response at one accepted endpoint. Each pressure
trial starts from the original accepted gas and modal states. The native
held-pressure response supplies the anchored final volume and actual wall
sweep; the gas receives the matching debit `W = p_trial * (V1 - V0)`.
The original species equations use that same `V1`, composition and
`U1 = U0 + Q - W`. Formation energy remains inside the original caloric
storage and is never supplied again as reaction heat.

The original Newton--Krylov owner solves the endpoint EOS pressure residual.
The adapter has explicit pressure/volume intervals, a bounded pressure
difference stencil, outer probe count and native solver budget. Every inner
mechanism retains its source, species, consumed-fraction, caloric and work
limits. Native gas-loaded clock admission still applies; it is a
frozen-composition restriction, not an error estimate for the whole reacting
system. Final admission independently checks pressure, anchored volume,
reciprocal native input work, EOS/work mismatch and the complete first law:
gas and modal storage changes plus native structural loss and exterior
pressure work equal external heat. Structural damping remains an outgoing
loss; it is not silently deposited in the gas.

`piston::reactive::boundary` exposes the same closed-parcel implicit mechanism
with explicit final volume, heat and work. It reuses the original composition
reconstruction and material energy inversion, including the zero-exchange
path. Heat entropy transport is `Q/T_end`; mechanical work is excluded from
that transport. Arbitrary enclosing mechanical trial probes defer only the
entropy inequality, which the native modal owner must enforce after pressure
and work close and before publication. The public standalone boundary
candidate always enforces it. All other source and conservation gates apply
to every trial. Refusal or cancellation leaves both original owners and their
histories available for deterministic retry.

The focused child `tests/solid_enthalpy/valve/gas_motion/reactive/coupled.rs`
specifies independent equal-mass isomer and one-DOF six-tet-cube references,
source-rate and stiffness substitutions, explicit volume/heat/work and native
zero-exchange parity, continued anchored motion/species/energy, and
budget/domain/caloric/cancellation rollback. The scalar oracle eliminates
the species and exact held-force wall response algebraically; it does not
reuse the runtime nonlinear solver. Coefficients are synthetic.
An earlier runtime draft passed the pinned Rust syntax parser before the
executor became unavailable. The final sources and new tests received
source/API review and independent reference arithmetic only. Native
compilation and numerical execution remain pending.

This path is first order: chemistry and gas pressure use backward Euler,
while the structural held-force response retains its native ZOH integration.
The original reaction-then-motion split remains available. The joint path
requires one original implicit mechanism and a shared gas inversion policy;
chemistry-only subcycling is refused so refinement covers the complete
coupled sample. It does not add flame/shock propagation, thermal wall
exchange, evolving structural material/basis, moving mesh reconstruction,
or composition-dependent radiation transfers. No nonlinear uniqueness or
unconditional convergence claim is made. The broader material-reality and
gas-chemistry beads remain open.


### Source-declared finite thermal reference (MR58/MR21 increment)

`finite_solid::thermal::MeanThermalReference` uses the existing mean engineering
expansion resolver with explicit reference and endpoint temperature roles.
The coefficient must carry that source convention; instantaneous small-strain
coefficients are not silently reinterpreted. Both contexts retain all remaining
coordinates and quantity descriptors. The selected law's actual source axes
and validity axes must be supported fixed external conditions, including when
a curve declares a dependence absent from its applicability bounds.

`FiniteSolidMaterial::resolve_mean_thermal` binds the thermal reference,
reference density and current-temperature constitutive model to one immutable
material card. With source length ratio `s=L(T)/L(T_ref)`, `Fe=F/s` and the
response per original reference volume is `W0=s^3 W(Fe)`, `P0=s^2 Pe`,
`A0=s Ae`. The original material checks the elastic gradient; `hyper3` also
checks the actual total Jacobian. Density and mass stay at the original
reference temperature. Uniform and per-element references use the same native
force, energy, tangent and dynamic operators.

Five new fixtures cover independent energy/stress/tangent formulas, free
expansion under rigid rotation, source-driven finite motion, heterogeneous
held temperatures, original mass, and source/domain/Jacobian refusals.
Source/API review and independent arithmetic were performed; native Rust
parsing, compilation and execution remain unavailable. These fixtures do not
supply experimental calibration.

This is thermal preparation followed by mechanics at held temperature per
element. It does not evolve a heat/entropy equation or account for work during
a temperature change. Isotropic reference preparation is supported;
anisotropic finite thermal mapping, coupled evolving thermomechanics,
inelasticity and nearly-incompressible stabilization remain open. MR58 and
MR21 remain open.


### Directional dielectric storage through the original phase/heat trajectory

The spatial conduction/enthalpy coupling now accepts the explicit dielectric
laws retained by `PlanarLeakyDielectric`. Its existing complete-bundle
refresh evaluates actual temperature/liquid fraction and conductivity while
preserving both principal permittivity claims, their independent material
frame and the original geometry frame. The original local midpoint `E.J`
is deposited once into the original enthalpy owner; field storage remains
`E.epsilon.E/2`. No extra electric state, material epoch or heat credit is
introduced.

The focused anisotropic dielectric child tests specify two continuing
charge/phase/heat steps against independent scalar energy/charge elimination,
scalar and tensor conductivity, independent dielectric-axis substitution,
source identity and phase-coordinate retention, wrong-frame/missing-source/
varying-permittivity refusal and late cancellation. The synthetic fixture
is not experimental calibration. Source/API review and independent scalar
reference arithmetic completed; reconstructed Rust tests remain unexecuted
and unparsed in the unavailable executor.

Permittivity and its frame remain source-declared constant through these
trajectories. Variable-permittivity work and reversible thermal entropy,
dielectric dispersion and structural motion remain outside this increment.


### Source-retaining live viscothermal atmosphere updates

`MaterialWallDg` and `RegionalMaterialWallDg` expose
`rebind_viscothermal_gas_field_controlled`. The existing prepare/publish
boundary now supports the native uniform circular diffusion transfer with an
explicit new-gas admission budget. It retains the complete source catalog,
query and composition, every wall material receipt and local Maxwell
history, and the original physical clock and propagator tolerance. Updated
gas field, native state and newly certified propagator publish together.
The prior lossless entry keeps its original admission contract.

The source integration regression changes temperature, pressure and
composition with nonzero acoustic and wall histories. It specifies preserved
physical memory and receipts, changed continuing pressure trajectories,
complete parameter-work accounting, and source/band/foreign-operator/late
cancellation rollback. Its narrow-duct wall fixture explicitly satisfies
the original thin-wall limit. These are synthetic numerical references;
the reconstructed final test has not been parsed or executed.

This path is a prescribed uniform background material intervention. Actual
wave and diffusion storage changes are external parameter work, not an
unmodeled heat credit. It adds no spatial gas mixing, mass/species injection,
mean-pressure traction or thermochemical reservoir solve. Geometry, radius,
native diffusion sections and thermal-arm support stay fixed.


### Same-endpoint reacting gas and material-bound inertial wall

`GasPiston::coupled_reactive_candidate` closes the existing implicit mechanism
at the actual moving gas endpoint: each pressure trial supplies the original
wall kinematics and `W = p (V1 - V0)`, then the original reaction owner solves
`xi_r = dt V1 rate_r(T1, p1, n1)` with `U1 = U0 + Q - W`. The shared bounded
scalar solver closes applied pressure against that same endpoint EOS. Native
wall kinematics and impulse accounting are shared with the nonreacting piston.
Fixed supports retain their original volume and balance pressure with a
zero-work support impulse. `MaterialGasPiston` keeps the original resolved
rigid-body inertia and source receipts through the same candidate.

Gas composition, temperature, transport, acoustics, volume, momentum and port
history publish together. Complete formation energy is already in gas U;
no reaction-enthalpy credit or numerical residual is added as heat. Independent
EOS-pressure, pressure-work, impulse, mechanical/total energy and signed
`Delta S - Q/T1` admission precede return. The latter includes the temporal
method's entropy and is not an isolated chemical-production claim. Every
pressure and chemical probe has explicit work/source-domain bounds.

Three focused synthetic tests specify independent endpoint caloric/EOS/rate
relations, fixed support, zero-rate nonzero-heat evolution, two continuing
exothermic compression intervals, source/bracket refusals and final-checkpoint
cancellation. Scoped Rust parsing/formatting, whitespace and source/API review
passed; native compilation and test execution remain unavailable because the
pinned constellation siblings and remote runners are absent and disk is full.
Independent reduced-equation arithmetic is not Rust execution or experiment.

This is homogeneous closed-gas kinetics with endpoint pressure/chemistry and
midpoint wall kinematics, generally first order. It requires the original
coupled implicit mechanism and shared gas inversion policy; chemistry-only
subcycling refuses. It adds no flame, shock, valve flow, shared crank, wall
material heating or guaranteed nonlinear uniqueness. The gas-control-volume
and gas-chemistry ladder beads remain open.


### Simultaneous reactive gas, finite material walls and signed heat contacts

`GasPiston::coupled_reactive_thermal_candidate` and its material-layer entry
close reciprocal wall heat against the same endpoint chemistry/pressure/work
solve. Each proposed Q starts from the original gas inventory and body
enthalpy: gas receives `Q_external + Q`, body receives `Q_body - Q`. The
actual endpoint gas composition/T/p and body enthalpy/phase determine the
original conductance law, and the shared scalar residual closes
`Q = dt G (T_body,end - T_gas,end)`. The material-layer entry reuses the
same-card enthalpy/conductivity condition owner and the reacted gas's original
mixture transport. The exposed patch cannot exceed the body's surface area.
`MaterialGasPiston` retains the original source-derived translating inertia.

`coupled_reactive_thermal_ports_candidate` uses that same heat-root owner for
multiple contacts. It retains individual signed heats and concrete source
records in stable contact-ID order, including opposing nonzero transfers at
zero net gas heat. Every contact sees one actual reacted gas endpoint.
Finite counterpart callers must debit those returned transfers exactly once
and publish their histories with the gas. The existing split assembly entry
points keep their declared split semantics; these are explicit coupled entries.

Independent gas/kinetic/body energy, applied/evaluated boundary heat, source
identity and nested pressure/chemical work counts accompany joint candidates.
Chemical formation and body latent energy remain in their original storage;
neither numerical residual nor a duplicate reaction/phase heat source is added.
The same endpoint kinetic policy, source envelopes, whole-interval clock and
cancellation boundary apply throughout all discarded and accepted trials.

Three new focused tests specify two continuing source-conductivity substitutions
with independent caloric/endpoint-rate/transport/acoustic equations, actual
source retention, cancellation and foreign-card refusal, and opposing signed
thermal contacts. Scoped Rust parsing/formatting, whitespace and source/API
review passed. Native compilation and execution remain blocked by the missing
pinned siblings/runners and full filesystem; independent scalar reference
calculations do not replace those tests or experimental validation.

The uniform finite body retains fixed geometry and negligible pV-work scope;
its phase chart does not move a melt boundary or alter piston inertia. Heat
and kinetics use endpoint evaluation with midpoint mechanical kinematics,
generally first order. Arbitrary spatial flow/flames, shared crank mechanics,
material-network wiring of this coupled policy and a complete wall-entropy
model remain outside this increment. Broad gas/thermal beads remain open.


### Regional atmosphere updates on continuing material walls

Both `MaterialWallDg` and `RegionalMaterialWallDg` expose
`rebind_spatial_viscothermal_gas_field_controlled`. The shared preparation path
uses the original spatial gas-memory transfer, retains all wall material
assignments and physical histories, rebuilds the original propagator at its
existing timestep and tolerance, then publishes field, state and propagator
together. Complete gas catalogs, actual local composition/temperature/pressure
and new port coefficients survive; the inlet does not stand in for the field.
The native weak-coordinate/arm conservation and external parameter-work scope
are specified in the fs-duct contract. No background thermal/species transport
or mean-pressure wall-prestress update is inferred.

A focused source regression definition changes the outlet atmosphere while
holding the inlet fixed, continues native pressure and regional wall motion,
checks the same clock and retained source catalogs, and exercises both wall
source layouts plus late cancellation and stale-propagator refusal. Syntax
parsing and source/API review passed; native Rust execution was unavailable in
this environment, so this is not a passing native or experimental result.
Focused native command:
`cargo test -p fs-couple --test material_duct_wall spatial::gas_regional_live`.


### Joint finite material heat, reaction and native elastic motion

`GasValveMotion::coupled_reactive_thermal_candidate` closes the original
implicit reaction/pressure/native-modal endpoint with a finite enthalpy body.
The existing bounded scalar feedback owner finds one signed transfer Q into
the gas. Every trial restarts ALL owners from their accepted histories:
`U_gas,1=U_gas,0+Q_gas+Q-W` and
`H_body,1=H_body,0+Q_body-Q`. Rates, gas pressure/work, source transport and
interface conductance use those same proposed endpoints. No reaction-first
heat lag, separate piston, constant-Cp surrogate or additional formation-heat
credit is introduced.

The existing source-bound series-conduction adapter also has a direct entry
point, `coupled_reactive_thermal_layers_candidate`. Body conductivity and
caloric state must retain the same material card and actual query point.
The original reacted gas supplies its own conductivity. The returned candidate
contains the original motion/gas owner, `UniformEnthalpyStep` and full account;
they publish together. Cancellation or source, bracket, pressure/chemistry,
clock, work or caloric refusal leaves every accepted history unchanged.

Independent measured balances enclose gas storage, body enthalpy, modal
storage, original structural damping and exterior pressure work. Q cancels;
damping leaves this enclosure and is not deposited again. Body entropy uses
the original chart's exact piecewise `integral(dh/T)`, including latent heat.
The native gas entropy gate remains active; final body, interface
`Q(1/T_gas,1-1/T_body,1)`, and enclosing endpoint entropy inequalities are
checked separately. These discrete inequalities include backward-Euler
integration entropy and are not a continuous entropy-production certificate.

This is a first-order endpoint reaction/conduction and native held-pressure
structural method, with explicit nested work/residual limits. The thermal body
uses a justified uniform-temperature, fixed-geometry, negligible-pV-work
enthalpy model and is a separate participant from the source-fixed elastic wall.
It does not silently update that wall's basis, expansion, strength or phase.
No flame, shock, phase-density motion, mechanical/fluid handoff, unconditional
convergence or continuum combustion accuracy is claimed.

Four synthetic native regression definitions cover independent scalar
elimination of heat, chemistry and six-tet elastic motion under conductivity
and caloric DATA changes; continuing latent-body/native-mode trajectories and
whole-system accounting; late source/clock/budget failures and deterministic
retry; and cancellation after actual chemical/mechanical/source evaluation.
Rust syntax parsing, bounded independent source/API review and reference
arithmetic passed. Native compilation/tests were unavailable (no Rust/Cargo,
DSR or RCH); no native or experimental success is asserted.
Focused command:
`cargo test -p fs-couple --test solid_enthalpy valve::gas_motion::reactive::coupled::thermal`.


### Spatial electrocaloric feedback through the original charge/enthalpy solve

`ConductionThermalNetwork::step_dielectric_thermal_candidate` resolves endpoint
conductivity and affine scalar/orthotropic epsilon(T) inside the existing
enthalpy residual and retained `PlanarLeakyDielectric`. It reuses the original
Joule heat router. Each original triangle additionally deposits
`Qrev=-P-delta(TS)`, where P is its native coefficient work and S is the
source derivative at fixed D. Reversible heat and Ohmic heat remain separate;
no numerical residual is reclassified as physical heating.

The final enclosing account is
`delta_U_field+delta_H_caloric-W_electrode-Q_external`.
Original complete source identities, claims, nonthermal axes, fixed material
frames, whole-volume thermal receivers and accepted clocks remain pinned.
Capacity bounds are summed over ALL triangles assigned to each caloric cell;
caloric plus field capacity must stay positive throughout the common source
interval and explicit |D| envelope. Actual accepted and proposed displacement
magnitudes both satisfy that envelope. The full field/thermal candidate
publishes together after source, charge, energy, capacity and cancellation gates.

This increment admits fixed geometry and axes, fully solid constant-density
mapped caloric charts and positive affine nondispersive linear epsilon.
The separate constant-epsilon trajectory retains its existing phase support.
Nonlinear polarization, phase-dependent field entropy, dispersion,
electrostriction, moving interfaces, unconditional entropy/passivity and
temporal-accuracy claims remain outside this model.

Three focused regression definitions cover independent oriented fixed-D
thermodynamic derivatives; both slope signs and the zero-slope limit over
continued spatial charge/heat and driven-electrode trajectories; and the
weak-leakage conserved-charge limit plus source, capacity, phase, budget,
late-cancellation and retry refusals. Independent scalar elimination checks
covered nine slope/orientation/drive combinations over two continued steps.
The weak-leakage reference endpoint is 300.50484323440674 K, versus
300.50484323440355 K in the exactly conserved-charge limit, with approximately
7.95e-12 J of physical leakage heat. These are independent reference
calculations, not executed native outputs. Final changed files passed Rust
syntax parsing and bounded independent source/accounting review; native
compilation/tests were unavailable (no Rust/Cargo, DSR or RCH).
Focused command:
`cargo test -p fs-couple --test conduction_phase_enthalpy dielectric::thermal`.


### Same-endpoint reactive chambers on the shared material exchange network

The shared circuit/enthalpy network now accepts native midpoint and coupled
reactive boundary drivers through one internal boundary interface. Electrical
route admission, receiving-cell heat aggregation, circuit advancement, the
outer thermal residual and final all-owner energy check remain single owners.
Existing public midpoint candidate/error names specialize the same generic
carriers and retain their original behavior and source records.

`network::reactive::step_with_contacts_candidate` binds optional original
reaction mechanisms to stable chamber IDs. Listed chambers close reaction
progress, pressure work, wall momentum and all thermal contacts against the
same proposed endpoint; unlisted chambers retain their native midpoint law.
Every trial restarts the original accepted gas inventory. The result retains
an explicit native-or-reactive boundary and its actual detailed receipt;
endpoint chemistry is never represented as a fabricated midpoint step.

`PreparedMaterialExchangeNetwork::coupled_reactive_step_candidate` and its
custom gas-law counterpart connect that policy to the original material bulk,
electrical contacts, thermal contacts, finite enthalpy cells and closed
radiation enclosure. Conductivity/resistivity, source/phase snapshots, heat
routes, surface-area admission and radiation closure use their existing owners.
Each signed gas contact is debited at its physical receiving cell, even when
other contacts cancel its net gas heat. Formation and latent energy remain in
original storage; reaction energy and electrical heat are not credited twice.
No accepted circuit, thermal, gas or moving-wall history changes on refusal.

Six focused regression definitions cover mixed reacting/ordinary chambers,
independent per-cell and all-owner energy equations, endpoint kinetic rates,
two continuing material-data substitutions, exact source retention, custom
heat laws, no-reaction equivalence, mapping/area refusals and late cancellation.
Scoped Rust syntax/formatting, whitespace and bounded source/API review passed.
Native compilation and test execution remain unavailable because pinned
constellation siblings and remote runners are absent and the filesystem is
full. These definitions and source checks are not executed native validation.
Focused command: `cargo test -p fs-couple --test material_exchange_network coupled_`.

This is a first-order homogeneous closed-gas model with independent translating
walls and fixed circuit topology/C/L and thermal geometry. The original coupled
implicit mechanism and explicit source/work limits remain required. Whole-
assembly interval refinement, rather than chemistry-only substeps, governs
coupled time resolution. It adds no spatial flame, valve mass flow, shared
crank or unconditional convergence guarantee. The broad physics ladder beads
remain open; this increment closes the previously separate material-network
connection for the coupled chemistry/heat/pressure policy.


### Continuing simultaneous reactions in the persistent material assembly

`MaterialExchangeState::coupled_reactive_candidate` and
`advance_coupled_reactive`, with corresponding custom gas-law entries, move
all endpoints from the prepared coupled network into the original live state.
Gas species and complete source catalogs, volume/momentum, port histories,
finite-cell enthalpy/event cursor and electrical charge/flux advance together.
The same clock/drive admission is shared with existing staged and flow routes.
There is no additional reaction stage, thermal deposition or constitutive law.

`CoupledMaterialExchangeAdvance` retains the original circuit/thermal report,
material/contact samples and each actual `CoupledChamberStep` variant. Chemical
progress remains inside the original reacting boundary account rather than
being copied into a second state or an incompatible midpoint report. Empty
reaction bindings preserve ordinary chamber integration. A final cancellation
checkpoint precedes the single live-state replacement; all prior work remains
a candidate, including source validation and collection of every chamber.

Two focused regression definitions cover two continued intervals, direct vs
preview vs live-state equality, exact source/report and port-history retention,
dismantle/resume, ordinary-path reduction, outer-work refusal and late source
cancellation followed by deterministic retry in a fresh execution scope.
Scoped Rust syntax/formatting, whitespace and bounded source/API review passed.
Native compilation and tests remain unavailable for the environment reasons
recorded above. Focused command:
`cargo test -p fs-couple --test material_exchange_network state_tests::coupled`.
The fixed-topology, independent-wall, closed homogeneous-gas and first-order
limits of the prepared coupled solver are unchanged.


### Finite plumbing and the native elastic-cavity state (MR05/MR06/MR07)

`GasValveMotion::{network_candidate,reactive_thermal_network_candidate,reactive_thermal_network_layers_candidate}`
compose the original shared-snapshot nozzle/species/enthalpy owner with the
original mesh/material-derived elastic motion. `ModalGasNetwork` assigns one
chamber identity and any bounded finite-reservoir plumbing; every inlet and
outlet samples the same accepted inventories. All transported enthalpy enters
storage once, including formation reference and flow work. A single endpoint
inversion policy is required for the chamber's flow and mechanical stages.

The reactive thermal path transfers once, then restarts every original coupled
heat/chemistry/pressure/motion probe from that same transferred chamber. Heat,
reaction rates, EOS and structural pressure work use the original simultaneous
endpoint solver. Thermal source selection remains the original complete gas/body
law. Anchored volume, modal basis/equilibrium, structural history, observers and
accepted clock are retained. The candidate returns the chamber, every finite
reservoir, optional finite body and all original stage reports together.

A fresh extensive balance encloses all final gas inventories, body enthalpy,
modal storage, original modal loss, actual exterior-pressure work and external
heat. Stream enthalpy cancels internally; chemical formation energy is never
added as heat. Source, aggregate donor withdrawal, pressure crossing, local
nonlinear, final mass/energy and cancellation failures publish no endpoint.

This is a first-order transport-then-coupled-motion composition over one native
sample, not an implicit flow network. The thermal body's geometry is fixed and
separate from the source-fixed elastic wall. No jet momentum, flame/shock,
thermal basis update, phase handoff, exact continuum accuracy or whole-network
entropy theorem is claimed. The existing thermal entropy report covers only
the post-flow evolution window. Refine the complete sample at flow events.

Six focused G1/G3/G4 regression definitions compare intake/exhaust against
independent ideal-nozzle, mixture caloric, isomer-rate, finite heat and six-tet
modal equations; they cover material/composition changes, reporting-ID
invariance, closed-flow equality, anchored continuation, late refusal and
cancellation/retry. Independent scalar references and Rust syntax parsing were
run in the editing session. Native Rust compilation/tests were not run because
Rust/Cargo and DSR/RCH were unavailable.


### Finite thermochemical background and continuing material-wall acoustics

`MaterialWallDg` and `RegionalMaterialWallDg` expose
`step_with_finite_gas_controlled` and the closed-inventory
`step_with_finite_gas_heat_controlled` convenience. The original native
propagator advances one sample in the initial medium. A repeatable original-owner
background proposal receives that original finite parcel, the same dt and trial
work out of gas; it returns a complete endpoint with separately identified heat
and transported enthalpy plus original reports/companion candidates. Existing
scalar feedback closes gas work against the actual native wave/diffusion
parameter-work increment. The heat convenience reuses the original source gas's
caloric/EOS/model-and-work update; no new thermal or flow solver is introduced.

Measured extensive gas storage, native storage, original port work and
unpartitioned loss independently close before the source wrapper, gas, state and
same-clock propagator publish together. Complete initial source-state equality,
actual circular volume, source/phase policy, fixed geometry, temperature/work/
error budgets and cancellation are admitted. Final nodal and native port
observations and ledger use the accepted gas, while held-port integrals retain
their actual substep. Kelvin/Maxwell wall histories and supported uniform
diffusion histories survive, and parameter work appears once across later ledgers.

Native dissipation remains exported energy: there is no inferred gas/wall
thermal partition. Scope is uniform well-mixed, fixed-volume, first-order
partitioned acoustics; no mean-flow momentum, pressure prestress, new material
kinetics, wall heating or spatial background update. External callbacks must
return immutable companion proposals and publish those only after acceptance.

Three integration regression definitions cover independent variable-Cp
heating/storage, original nozzle species transfer with continuing wall/diffusion
memory, and late transactional refusals. All eight changed Rust files passed
syntax parsing. Independent manufactured-profile arithmetic was evaluated;
it is not a native acoustic trajectory execution. Native Rust compilation/tests
were not run because Rust/Cargo and DSR/RCH were unavailable.


### Ceramic material temperature in the original piezoelectric/enthalpy solve

`circuit_enthalpy::piezoelectric::thermal` advances the existing
`PiezoelectricCircuit` and original finite enthalpy/Newton timeline atomically.
Each nonlinear trial refreshes the original rod at trial endpoint T, retaining
accepted charge, flux, extension, momentum, geometry and source identity, then
solves the original Kirchhoff/piezo descriptor. Native rebind exposes
P=Psi(T1,x0)-Psi(T0,x0); Q_rev=-P-[T1*S1-T0*S0] enters the same whole-rod
caloric cell separately from native resistor Joule heat. The final account
closes ceramic U, rod kinetic energy, other circuit storage and caloric enthalpy
against electrical work, force-displacement work and external heat.

The coupled owner requires a fully solid constant-density chart covering the
whole source interval, matching entire rod volume/mass, and positive
caloric-plus-field capacity over explicit strain/D bounds. Original small-signal
and ideal-source compatibility gates remain authoritative. Ordinary fixed-ceramic
stepping has P=0 and retains the existing native trajectory. Unsupported
source/density/phase, capacity/coordinate bounds, native failure, budget or
cancellation return no joint endpoint. This endpoint-T coefficient split is
generally first order; energy closure is not an unconditional entropy/stability
or time-error certificate. Thermal strain, depoling, pyroelectric offsets and
irreversible ceramic losses require their own admitted laws.

Four focused regression definitions cover independent full-block potentials and
derivatives, signed poling reversal/zero coupling, continued material-load motion
and physical energy, an isolated fixed-charge caloric limit, fixed-source native
equality, and late source/budget/cancellation refusal with deterministic retry.
All nine changed Rust files passed syntax parsing; independent scalar equations
were evaluated. Native Rust compilation/tests and rustfmt were not run because
the required Rust/Cargo and DSR/RCH environment was unavailable.


### Shared finite gas transport and continuing acoustic chambers

`advance_finite_gas_dg_transport` composes the original combined nozzle and
Maxwell-Stefan transport with multiple existing source-backed uniform or
regional-wall DG owners and finite reservoirs. All rates, gross species
withdrawals and partial enthalpies use one immutable original gas snapshot over
the common native sample duration. Transport reconstructs each inventory once;
repeated local gas/acoustic work feedback reads that transported proposal
without resampling or withdrawing again. Each duct's complete original source
state and clock must match its native parcel. Unbound reservoirs use the
original caloric heat owner.

`GasTransportStep::admit_endpoint_nozzle_pressures` delegates the original
active-flow ordering gate against the complete final heat/work-adjusted map.
Actual global gas and wave/wall storage, native loss, external heat, port work
and mass close before both maps publish. Transported enthalpy and gas-to-wave
work are internal and never counted again. Per-duct source fields, ports,
diffusion memories and physical Maxwell histories continue.

These are distinct well-mixed fixed-volume chambers: no continuous spatial
pressure-equilibrated field, acoustic scattering through restrictions, jet
momentum/noise, static pressure prestress, wall thermal evolution or guessed
gas/wall loss partition. Transport entropy applies to its original stage;
later heat/work does not inherit an entropy certificate. Original source,
near-isothermal diffusion and gross-withdrawal gates remain; the composed
evolution is first order and requires physical timestep refinement.

Three focused regressions cover independent binary flux and source substitution,
multi-owner continuation, global energy/species conservation and late final
pressure/work/cancellation rollback. All six changed Rust files passed syntax
parsing. Independent diffusion and caloric reference arithmetic was executed.
Native compilation, Cargo tests and rustfmt remain unexecuted in this environment.


### Reciprocal spatial elastic work and finite enthalpy

`solid_enthalpy::thermoelastic::ThermoelasticSolid` retains the original
`SolidEnthalpy` binding, native compatible equilibrium/stress, actual assembled
nodal loads and finite thermal state. One declared source interval per reference
tetrahedron uses the source-backed affine isotropic thermoelastic potential:
E(T), fixed Poisson ratio and constant instantaneous expansion. The same
source/card/eigenstrain and delivered reference mass remain pinned. It admits
a fully solid constant-density interval within the original complete caloric
chart, preserving later latent/liquid segments. Any latent segment intersecting
the admitted interval, or any actual phase/source-domain exit, refuses. The
interval's actual caloric capacity plus the conservative fixed-total-strain
field capacity must stay positive.

Each bounded heat fixed-point trial starts from accepted thermal and observer
histories, runs the original thermal substeps/events with additional reversible
powers, then resolves local sources and solves the original P1 load,
displacement and stress owner. With unchanged homogeneous supports, local work
W_e = V_e (sigma_0 + sigma_1) : delta(epsilon) / 2 sums to actual averaged
endpoint nodal work (f_0 + f_1) dot delta(u) / 2. Q_rev,e = W_e - delta(U_e)
enters that element's original cell exactly once.

The independent free-energy identity separates expansion work
-V (sigma_0 + sigma_1) : delta(epsilon_free) / 2 and stiffness work
V e_0 : delta(C) : e_1 / 2, while T*S remains part of physical U. Numerical
residuals are never converted to heat. Local source/native stress and energy
parity, per-element applied/source heat mismatch, mechanical work, thermal
storage and the combined internal/caloric/external balance all gate publication.
Original observer source, conductance, transfer, interval and event hooks are
forwarded; only the final accepted observer clone is returned.

This is quasistatic continuation with declared trapezoidal endpoint mechanical
work and the original explicit thermal substeps. Changing-load/temperature
coupling is generally first order. No exact continuous-path work, unconditional
entropy/stability, inertia, plasticity, moving supports, follower loads, changed
heat-transfer geometry or phase-flow claim is added. Source, phase, capacity,
native solve, iteration-budget, observer or cancellation refusal leaves accepted
history unchanged.

Three synthetic independent regression definitions exercise free/constrained
finite heating, heterogeneous connected tetrahedra with changing loads and
material-only substitutions, full phase-chart retention and late refusal/retry.
All four changed Rust files passed syntax parsing; independent two-element
reference equations and bounded source/API/physics reviews were completed.
Native compilation, Cargo tests and rustfmt remain unexecuted in this environment.


### Shared unsteady planar and circular fluid passages

`slit_channel::FixedChannel` retains either the original planar slit or an
explicit circular passage. Both use one finite-volume momentum owner with
geometry-supplied cell areas, masses and no-slip shear conductances. Circular
annuli implement radial symmetry at the axis and the half-cell outer-wall
condition. The same pressure/flow condensation, continuing velocity field and
signed work account distinguish kinetic storage, physical viscous loss and
backward-Euler numerical loss. Failed proposals publish no accepted history.
The original public `SlitChannel` constructors and methods remain available;
coupled consumers select `.slit()` or `.circular()` for geometry-specific data.

Material-modal reed files accept `supply circular LENGTH_M RADIUS_M CELLS`.
The shared DG grammar accepts
`dg_terminal_circular PORT_ID RADIUS_M LENGTH_M CELLS MAX_MACH`. Each uses
the actual retained local gas and mechanical clock. Circular and slit terminal
fields share the declared cell budget and original duct/radiation transaction;
microphone pressure follows the accepted physical flow. CLI provenance retains
the actual radius and radial resolution. No instrument-name dispatch, fitted
oscillator or independent material constants are introduced.

This is a fixed, fully developed, incompressible laminar passage reduction.
Additional passage length must not duplicate either acoustic duct. Entrance
development, taper, turbulent transport, compressible passage storage and moving
walls require separate admitted models. Regression definitions cover radial
Poiseuille and temporal Fourier–Bessel refinement, signed work and rollback,
material/atmosphere changes through reed and flute mechanics and pressure,
and the original uncondensed multiport equations with mixed passage geometries.

### Moving material fibers retain solved normal contact

The existing accelerating-fiber owner now carries canonical normal receivers
along the same material hair crossing as tangential contact. It reprojects the
original normal basis at each accepted endpoint and supplies that geometry to
the original simultaneous unilateral-contact/friction solve. Indentation,
release, material-domain checks and traction therefore share one physical step.
String-side station motion, retained body/cavity waves, and the source-bound
normal law remain with their original owners. Noncanonical receiver maps refuse
because they do not identify a physical material station.

`FiberSupportFrame::normal_attachment_work_j` is the fiber-side contribution of
the existing endpoint product rule. It is included in laboratory support work
and already appears in native obstacle work; the laboratory balance subtracts
that share from the native obstacle account exactly once. String-side migration
and independent normal indentation retain their separate accounts. Support
kinematics and accounts publish only after the complete mechanical step accepts.

The material-strings executable composes its existing fiber-contact profile with
`bow_fiber_program`, including continuous reversal and simultaneous bow station
motion. No force envelope replaces solved normal reaction. Normal translation
of hair supports, stick flexure, and transport of rosin/shear/thermal patch
histories remain outside this axial-support model. Regression definitions cover
canonical projection, normal loading/release, laboratory work, refusal/retry,
callback partition identity and source input through the physical PCM path.

### Spatial current with sourced thermal contacts and radiation

`PreparedHeatNetwork::step_conduction_candidate` binds the existing spatial
Ohmic field to the original complete heat assembly. Every implicit trial uses
the same cell enthalpy endpoint to reevaluate material conductivity, pinned
thermal-interface resistance and pinned surface emissivity. The original
spatial field deposits its actual per-triangle E.J integral into the admitted
cells exactly once. The bulk-only and complete-heat APIs share one private
current/enthalpy solver and retain their own complete source samples.

The original geometry/volume, material/caloric identity and phase-density
admission remains mandatory. Electrical graph bindings refuse on this path
because no graph circuit is being advanced. Thermal contact resistance is not
an electrical contact model. Closed-enclosure radiosity remains internal heat
exchange; its signed numerical net is checked against the original total-energy
allowance. Geometry, contact areas, emitting areas and view factors are fixed.
The composition does not imply volume-changing melting or moving interfaces.
Regressions define an independent two-body nonlinear balance, source-only
emissivity/contact substitutions changing current and temperature, continuation,
bulk-path equivalence and atomic domain/cancellation refusal.

### Source-bound regional acoustic networks and their exterior

`NetworkDrivingPoint::response_regional` and `response_coupled_regional` use
`fs_duct::network::NetworkMedia` on the original nodal/chain solve. Each exact edge ID
retains its original gas record; an independently declared exterior supplies
the external loads. No mixture or wave impedance is averaged at a junction.
`realize_line_regional` shares the existing operator-realization ladder and
passivity gate. Its nominal delay uses the longest shortest travel path with
each edge's own L/c, rather than assigning one speed to the complete graph.

`RadiatingNetwork::new_regional` retains those same records. Its full BEM
opening matrix and coherent microphone fields use the same exterior medium,
surface, flow basis and solved opening flows. Changing only an interior gas
therefore changes excitation of the exterior without changing that exterior's
operator. Changing the exterior changes both reaction and observation. Raw BEM
matrices are not clipped, diagonalized or repaired for passivity.

The `network_exterior` executable accepts one `--edge-gas EDGE FILE` per edge
and `--exterior-gas FILE`. It reuses the original gas grammar and native
thermochemical pack loader, including explicit water/phase policy. This is
static, no-mean-flow, rigid-wall acoustics at one common background pressure;
composition and temperature may vary by edge. It does not evolve mixing, heat
transport, reaction or a moving gas interface. The exterior remains homogeneous
with one uniform-flow coordinate per declared opening. Regression definitions
cover uniform equivalence, local travel times, source-only gas substitutions,
full-matrix reaction/pressure and lossless exterior-work closure. Compiler
checks do not establish numerical convergence or instrument realism.

## Prescribed-volume finite gas networks

`prescribed_gas::step_candidate` composes the original same-snapshot finite-nozzle network with original NASA-9 midpoint heat/volume proposals for explicitly driven chambers and finite reservoirs. All drives and endpoint solvers cover the same parcel identities; concurrent inlet/outlet enthalpy is counted once. Gross gas `p_mid * dV` work, exterior-pressure work and net work to the prescribed drive remain separate. Actual local storage balances include transported enthalpy; the complete finite network closes `sum(dU) + sum(W_gross) - sum(Q_external) = 0` before all original owners can be returned atomically. The original graph-size, source, inventory, withdrawal and pressure-crossing limits remain active.

This flow-then-motion split is first order. Prescribed-drive mechanics, timestep/event resolution, combustion and empirical valve applicability remain with callers. G0/G1/G3 tests cover independently computed constant-Cv pressure work, analytic Otto efficiency under refinement, variable-Cp compression against an entropy oracle, simultaneous inlet/outlet flow with motion, species conservation and link-order invariance, held-volume heat, and late-refusal/cancellation rollback.

## Stateful tangential contact on the shared modal network

`friction::tangential::TangentialModalSystem<L>` accepts the existing
`fs_tribo::TangentialLaw` directly. It uses the complete original bilateral
network's released increment and held-force compliance, then applies that law's
work-conjugate effort to the same network. Both body reactions, material memory
and admitted acoustic loads remain in the original evolution. Reversible
interface storage is included in the total-energy ceiling and work audit;
only the law's irreversible loss is exported as friction heat. Constitutive
state, mechanics and observations publish atomically. Actual displacement,
constitutive storage, physical heat and external work are checked together.

The current adapter has one fixed tangential coordinate and a prescribed normal
load. Explicit pressure and temperature inputs pass unchanged to the law; this
adapter alone does not evolve thermal state or establish material authority.
Transverse constitutive response refuses rather than losing an unassembled work
coordinate. G1/G3 regressions compare elastic shear against an independent
two-oscillator load line, close reversal hysteresis work, reproduce the original
Coulomb consumer and preserve loaded history on refusal/cancellation. The
same public consumer runs the analytic Cattaneo--Mindlin law: its accepted
motion reproduces the Hertz tangential backbone and its path-average effort
closes body work, reversible contact storage and microslip heat together.
Loaded reversals retain the original contact memory and reject a changed normal
load without modifying mechanics or history; the law's fixed-load scope remains.
Each accepted sample represents one constitutive displacement segment, so time
refinement must resolve reversals within a sample. A three-rung elastic-limit
test checks second-order convergence to the independent stiffness-shifted
oscillator frequency; energy closure alone is not a time-accuracy claim.

### Two-axis retained shear on the original modal response

`friction::tangential::planar::PlanarTangentialModalSystem<L>` attaches both
physical tangent directions to the original modal bodies. It computes the full
2-by-2 response `C_ij = b_i^T D_network b_j` through the existing exact network,
including unequal principal responses, cross-axis terms, bilateral connections
and retained material/acoustic loads. Both directions on each side retain the
same original component and modal basis. The caller supplies their projections
in an unchanged orthonormal physical tangent frame; no shape is normalized or
replaced with an effective isotropic mass.

The shared `fs-tribo` matrix adapter owns the bounded constitutive solve. Every
trial starts from the accepted law state. The two held efforts then stage one
original network update, and actual driver-minus-body displacement must agree
with both returned constitutive increments. Reversible shear storage, physical
friction heat, original network losses, external work and driver work close
together before the network, history and sample clock publish. Source refusal,
nonconvergence and cancellation leave all accepted state unchanged. Response
columns and force workspace are retained between samples. Mechanical-only
steps take explicit normal load and temperature and export heat. The optional
`planar::heat` composition solves temperature with the original finite thermal
network; normal geometry and tangent-frame evolution remain caller-owned.

`tests/modal_planar_tangential.rs` checks an independent coupled elastic matrix
solution, nonproportional two-axis viscoplastic histories, basis covariance,
anisotropic static/sliding Coulomb behavior, singular response and exact retry
after refusal. This unlocks the original vector layer law in modal motion; it
does not turn the scalar Hertz--Mindlin law into a vector partial-slip solution.

### Simultaneous two-axis retained shear on shared modal bodies

`friction::tangential::multiple::MultiTangentialModalSystem<L>` advances a
bounded set of original `TangentialLaw` contacts in one mechanical transaction.
Each contact has its own law, retained history, physical limits, prescribed
normal load, drive and fixed pair of tangent projections. All response entries
`C_ij = b_i^T D_network b_j` are computed through the original modal response,
including cross-contact and cross-axis terms. Shared bodies, bilateral links,
material memory and admitted acoustic loads retain their existing evolution.
Setup caches the complete matrix and its shared constitutive admission; it
checks contact, matrix-storage and projection-work budgets before allocation.

The shared `fs-tribo` simultaneous solver queries every law from its accepted
history. Its held efforts are summed before one original network update.
Actual driver-minus-body motion is then checked independently for every
contact, including local stored-energy change and dissipated work. Original
network storage/loss, all contact storage/heat, external work and all drive
work must close together. Only then do the original network, all contact
histories, sample clock and observations publish. A late physical-budget
failure or source refusal preserves the complete accepted state and permits
an identical retry. Reversible contact energy is included in the network's
total-storage ceiling; irreversible heat is exported once per contact.

This retains the original time integrator and constitutive displacement-segment
meaning. It is a frozen-geometry, prescribed-normal-load contact interval,
without contact discovery, normal/shear evolution, tangent-frame transport or
a multiple-interface thermal solver. Redundant response rows are admitted
without physical regularization when the constitutive problem resolves.
Nonunique friction roots and iteration limits retain the shared solver's
explicit selection/refusal; no universal convergence or ordering-independence
claim is made. Scalar Mindlin keeps its transverse-motion refusal. Setup is
reused, but candidate constitutive histories and reports may allocate.

`tests/modal_multi_tangential.rs` checks an independent four-axis elastic
solution and all original body recoils, nonproportional plastic histories,
independent tangent rotations and contact permutation on a resolved problem,
simultaneous static/sliding behavior, redundant contact rows, loaded-history
rollback/retry, and a three-rung comparison with the analytic elastic frequency
of a body loaded by two retained contacts.

### Stateful shear memory and finite-body thermal evolution

`friction::tangential::heat` connects `TangentialModalSystem<L>` to the original
`FrictionHeatRoute` and finite `DepositionSimulation` in one accepted sample.
The existing nonlinear thermal owner solves original cell enthalpy increments
and a common interface temperature. Each trial supplies that temperature to the
same immutable `TangentialLaw`, repeats mechanics from accepted history, and
routes only its irreversible loss into the original cells. Geometry and actual
endpoint conductivities determine the heat split; no fraction or thermal store
is introduced. Reversible shear storage remains mechanical and still counts
toward the total-energy ceiling. Friction heat cancels once in the combined
mechanical/caloric work balance.

Prediction and publication are separated in the original tangential consumer.
The final accepted root is explicitly restaged after numerical probes; all
thermal admission, source callbacks, energy checks and cancellation polls finish
before mechanics, constitutive history and caloric state publish together.
Pressure and normal load remain explicit. A prescribed interface temperature
is refused because this operation solves it. Source callbacks may supply
transport and receipts but cannot replace the retained law or checkpoint.

Temperature-dependent dissipation is supported when the law admits the solved
temperature. Frozen-temperature material laws retain their refusal instead of
being reconstructed with unloaded history. Variable elastic coefficients need
explicit free-energy/entropy work; this adapter does not infer that coupling.
The fixed-property native Mindlin law therefore supplies a declared constant-
property heating rung. Moving footprints and thermal geometry evolution remain
with their existing dedicated owners.

Five groups in `tests/modal_tangential/heat.rs` cover native Mindlin heating and
reversal storage, an independently reduced thermal-softening fixed point using
the original Coulomb solver, ordinary hot-to-cold conduction at zero friction,
loaded-history refusal/cancellation/retry, and frozen-temperature law refusal.
The tests use original mechanical and thermal owners; their complete native
Cargo execution remains pending the upstream runtime build resource limit.

## Shared finite-network friction heating

`friction_heat::FrictionHeatRoute` binds an existing perfect-contact link between
adjoining orthogonal, uniform prismatic cells. Its massless common-temperature
interface derives each conductance from the original face area, that body's
half-cell distance and its endpoint material conductivity. The old direct link
is suppressed only in the thermal residual while surface coverage is retained.
Signed body heat includes ordinary hot-to-cold exchange; only the sum equals
physical friction work. Original cell storage, other sources and the thermal
energy account remain authoritative. No prescribed heat fraction or extra heat
store is introduced.

`FrictionModalSystem::step_with_friction_heat` solves the original scalar modal
friction proposal and this thermal interface in the same bounded Newton
transaction. Every temperature trial repeats mechanics from accepted history,
including existing nonlinear storage, and regenerates a temperature-specific
Stribeck load line. Retained dry-interface cards preserve their original law and
applicability; unbound laws retain explicit caller authority. Original physical
friction work enters caloric cells once. Independent cell, interface, thermal
and combined energy gates pass before either owner publishes. Source refusal,
iteration budget, cancellation or mismatched clocks preserve accepted motion,
law, branch history and thermal state.

Focused G1/G3 tests cover an independent two-free-body/finite-capacity solution
with temperature-dependent friction, equal mechanical loss and thermal gain,
the original Stribeck root at the solved temperature, rollback and exact retry.
General thermal tests additionally cover unequal half cells, ordinary
conduction, analytic steady temperature and three-level spatial convergence.
This backward-Euler finite-volume route requires time and spatial refinement.
Moving footprints, finite third-body storage, wear, thermal expansion and
flash-temperature closures remain outside this implementation.

## Physical stopping and reversal with finite thermal feedback

`SlidingPair::step_collinear_candidate` advances relative motion on one declared fixed tangent axis while retaining common transverse motion. It refuses transverse relative velocity/acceleration. Each sliding segment solves the original source/enthalpy feedback at midpoint slip with a bounded safeguarded secant/bisection solve of the reduced impulse. A stop is localized before the original static-capacity test admits sticking or reversal. Actual traction times sliding distance supplies friction heat; numerical kinetic loss is never deposited as heat. All segments, body momenta, original slip/work histories and finite caloric stores publish atomically after independent local and cumulative momentum/energy gates.

Constant held Coulomb coefficients recover analytic stopping/reversal. Endpoint temperature/speed-dependent coefficients still require time refinement; a sign bracket proves neither uniqueness nor a continuous first-event theorem for arbitrary nonmonotone curves. Thermal weakening within an initially stuck segment is not separately localized. G0/G1/G3 tests use source-resolved masses and verify complete conversion of relative kinetic energy at a stop, same-call forced reversal, independent finite-temperature feedback, static/common transverse force work, and late-refusal rollback. The existing general planar endpoint integrator remains available.


## Finite caloric feedback for retained partial-slip contact

`tribothermal::partial_slip` closes the original loaded two-region contact
history against the original finite enthalpy cells. Its immutable, explicitly
declared affine stiffness `k(T_a,T_b)` remains positive throughout its admitted
temperature rectangle. The model identity is an input; it does not establish
measured or source-derived contact stiffness. The original caloric objects and
transport laws retain their source domains.

For both native regions, `psi = sum(k z^2/2)`, each surface entropy contribution
is `s_i = -partial psi/partial T_i`, and contact internal energy is
`U = psi + sum(T_i s_i)`. Each cell receives its declared share of physical
microslip heat plus reversible exchange `-T_i,new * delta s_i`. Parameter work
and backward-Euler numerical loss are never relabeled as heat. Native contact
work, cell residuals, thermal balance, and the combined balance
`delta H + delta U + D_numerical - W_wrench - Q_external = 0` pass independent
gates before the coupled candidate can publish. Retained strains, slip
branches, frame transport, source/caloric object identity, and the physical
clock continue together; late failures and cancellation publish neither owner.

The fixed-coefficient ingress holds normal patch/load and Coulomb coefficients
fixed and uses an explicit constant physical-heat partition. The material
friction ingress below re-evaluates pinned Coulomb sources. Neither claims normal work,
expansion, wear, flash temperature, sourced stiffness admission, Euler-body
integration, total entropy production, uniqueness, or unconditional convergence.
Endpoint coupling is first order and needs independent timestep refinement.
Six focused integration references cover reversible elastic cooling, stationary
loaded heating, analytical Jenkins heat/numerical-loss separation, unequal
surface temperatures with conduction and frame transport, and domain,
rollback, cancellation and deterministic replay behavior.

## Finite third-body friction heating

`friction_heat::third_body` places physical friction loss in an existing finite
enthalpy layer cell. Its original mass and caloric law retain heat; both body
conductances come from actual half-cell geometry, endpoint conductivity and
source contact resistance. Independent layer and whole-network accounts keep
layer storage, signed body conduction, outside heat and friction work distinct.
An insulated layer retains heat. Reversible elastic and coefficient work remain
with their constitutive owners.

`FrictionModalSystem::step_with_third_body_heat` accepts original mechanical
motion, friction history and the finite layer/body thermal network together.
The layer's own caloric law supplies the same-step friction temperature, with
no extra temperature store. Five thermal and two native modal regressions cover
analytical finite-capacity/resistance balances, actual viscoplastic heat,
insulation, cancellation/retry and an 8/16/32 timestep ladder with approximately
first-order error reduction. This fixed-geometry rung uses one uniform,
orthogonal layer cell with two opposite conductive faces; spatial uniformity
and loaded thermoelastic applicability remain explicit model requirements.

## Resolved-conduction flash comparison

`friction_heat::flash::FlashHeatHistory` observes accepted, conduction-derived
signed body heat as bounded piecewise-constant flux pulses. Heating and cooldown
reuse the original `fs_tribo::flash_temperature_candidate` half-space kernel;
`FlashTemperatureInput::for_duration` supplies physical exposure time without
inventing traverse geometry. Stable pulse superposition avoids cancellation of
nearly equal square roots. Original route/material/geometry and initially
uniform body temperatures bind the comparison. Refused receipts, exhausted
history and cancellation leave the observer and thermal owner unchanged.

The oracle is read-only: it does not route production heat, create storage,
prescribe body fractions or replace an accepted temperature. Its declared limit
is stationary, one-dimensional, constant-property half-space conduction.
Three resolved spatial rungs give maximum relative temperature-rise errors
1.5781%, 0.3687% and 0.07344%; variable heating/cooldown and exact legacy-duration
equivalence are also checked. These are numerical model comparisons with
synthetic analytical materials, not experimental qualification or a validity
claim at arbitrary moving-source Peclet number.

### Reacting chambers with prescribed volume

`prescribed_gas::reactive::step_candidate` composes the original shared-snapshot
nozzle transport, `GasReaction::parcel_candidate`, and midpoint heat/volume
owners over one physical interval. Incoming species react in the actual mixed
chamber before prescribed boundary motion. Held reservoirs may also react.
Reaction catalogs retain the original caloric/transport/rate source identities,
exact stoichiometric conservation certificates, reverse equilibrium closure,
consumption and entropy admission, and selected explicit, implicit or chemical
subcycling policy. There is no duplicated chemical inventory, EOS or heat law.

All chamber endpoints publish atomically. Original-to-final local and complete
network balances include total species formation internal energy, transported
upstream enthalpy, external heat, and actual midpoint pressure work. Chemical
formation energy and enthalpy transport receive no second heat credit. Exterior
pressure work and net work to the prescribed drive remain separately reported;
the mechanical consumer must accept the reciprocal drive work. An independent
network molar audit measures each actual species change minus its summed
stoichiometric production across all chambers and every reaction channel.

This flow/chemistry/motion split is first order even if the chemical substep is
implicit or adaptively refined. Whole-step refinement is required; chemistry-only
refinement does not bound splitting error. Products leave in the next flow
interval. No simultaneous chemical/mechanical endpoint, resolved flame, ignition
time accuracy, engine efficiency validation, phase change or measured reaction
mechanism is claimed. The original nonreacting API retains its exact path.

G1/G3 tests in `cavity_phs::piston::reactive::tests::prescribed` independently
reconstruct formation-energy-driven expansion, simultaneous inlet/outlet
enthalpy and live mixed-state kinetic scales, backward-Euler reaction extent,
and the adiabatic reversible-isomer solution under complete-interval refinement.
They also check exact empty-chemistry equivalence, deterministic link reordering,
late source/identity failures, cancellation and replay without partial state.


### Reciprocal inertial thermoelasticity

`solid_enthalpy::thermoelastic::dynamic` launches the original sparse
average-acceleration tetrahedral trajectory from admitted loaded equilibrium
and supplied nodal velocity. Later states evolve actual inertia. It retains the
existing source-pinned affine E(T), fixed Poisson ratio, constant instantaneous
expansion potential, finite caloric charts, stationary supports, fixed reference
mass and positive total-capacity/strain/source admissions.

Each cell receives reversible `Q = W_stress - delta(U_elastic)`, with
`U_elastic = Psi + T*S` and source-derived `S = -dPsi/dT` at fixed total strain.
Native kinetic energy stays separate. The accepted balance is
`delta(K + U_elastic + H_caloric) + D_Rayleigh = W_external + Q_external + E_incoming`;
both native mechanical and thermal gates remain active. Rayleigh loss is
explicitly exported without invented spatial heating. Supplied damping must
exclude mechanisms already resolved by coupled transport.

Every nonlinear trial restarts both owners from accepted history; only the
balanced thermal state, native checkpoint, source state and original observer
are returned. Matrices and step count survive changing endpoint material
lifetimes without a second solve. Cancellation, source/strain exit, iteration
or energy refusal leave both owners unchanged. Accepted nodal motion and native
stress remain exposed.

G1 compares heated inertia and multi-step adiabatic oscillation against
independent consistent-mass, momentum, free-energy and finite-capacity equations.
G3 covers uncoupled/damping limits and native continuation; G4 covers late
refusal, state pairing and replay. Five new dynamic/checkpoint references and
three existing reciprocal quasistatic references pass focused native execution.
Fixtures are authored, not measured-material validation. Endpoint coupling
generally has first-order path error and needs time refinement; no exact entropy
conservation, finite-strain/phase dynamics, moving supports, downstream contact
recoupling or arbitrary initial free-body state is claimed.

## Moving finite-network friction heat and local modal feedback

`friction_heat::moving` binds original orthogonal prismatic strip faces and a
fixed footprint width and length. Linear material/footprint paths are integrated
at every cell-edge crossing, with one solved interface temperature for each
overlapping cell pair. Original cells retain the thermal wake, longitudinal and
depth conduction, and all storage. Ambient loss excludes the mean covered area
without persisting fictitious contacts. A sweep broadphase and one-pass exposed
area gathering avoid a quadratic traversal of the original thermal links.

Full footprint coverage over the entire interval, unchanged original geometry,
material and mass, explicit overlap budgets, and absence of competing cross-body
links are required. Heat flows are derived from local conductivities, actual
half-cell distances and common interface temperatures. There are no prescribed
body heat fractions or extra thermal stores. Endpoint temperature and flux data
use backward Euler; exact geometric overlap does not remove the need for
independent time refinement.

`FrictionModalSystem::step_with_moving_friction_heat` resolves each local
Coulomb law at its own interface temperature with a supplied mean normal force.
Normal loads must sum to the original drive load. Static and kinetic capacities
sum in the existing common-slip mechanical solve; local heat is local traction
times the same actual solved slip. Every nonlinear trial starts from accepted
mechanical/thermal state, and whole-system energy admission precedes one atomic
publication. A retained single-interface source card and non-Coulomb scalar
composition refuse. The geometry path remains an explicit input from the motion
owner; it is not inferred from a provisional friction iterate. Curvature,
changing width, transverse gradients and finite-patch partial slip are outside
this common-slip strip rung.

Three native thermal regression groups check multi-cell crossing conservation,
path reversal, refusal/rollback, thermal wake, stationary-limit recovery and a
Jaeger moving-band comparison at Peclet numbers 0.1, 1 and 10. Across three joint
meshes every local temperature-field error decreases; finest relative field
errors are 0.949%, 1.641% and 1.353%, with mean-temperature errors below 0.26%.
Independent depth and time ladders show observed orders 1.945 and 0.978; axial
refinement also reduces the temperature increment. Two native modal regression
groups check local temperature-dependent friction and heat, original mechanical
energy, and bitwise-identical retry after a staged refusal.


### Source-resolved friction on a loaded caloric contact

`CaloricPartialSlip::with_material_friction` optionally attaches the original
`MaterialSlidingContact` resolver and its pinned friction claims. Admission
requires the same ordered interface/history/provenance, original caloric
objects and cell order, accepted initial coefficients, and declared heat split.
The affine contact free-energy model retains its own independent identity.
Rebinding an attached source is refused.

Each thermal endpoint re-evaluates the original Coulomb sources at both actual
surface temperatures, held normal load divided by the binding's explicit
nominal area, and actual translational speed `norm(creepage * rolling_speed)`.
Those coefficients update the original native capacity while preserving all
retained elastic/kinetic history. Source receipts are returned with the accepted
report. A temperature-driven capacity reduction may relax loaded strain and
produce physical heat even with zero imposed motion; numerical loss remains
separate from both physical and reversible heat.

This material ingress admits translational contact and a fixed prescribed heat
split. It refuses retained twist, relative torsional motion, distance-dependent
source axes and semi-infinite effusivity partitioning; it invents no torsional
speed or preattachment sliding history. Diagnostic relative travel starts at
attachment and advances only with accepted physical motion. Original source
range failures and all other late refusals leave both accepted owners unchanged.

Three additional G1/G3/G4 references independently solve the coupled
mu(T)/finite-capacity two-slider equations, stationary thermally triggered
relaxation, and pinned-source/domain/unsupported-motion behavior. Together with
the six fixed-mode references, all nine focused tests pass. These use declared
synthetic material claims, not measured friction validation. The fixed normal
patch and first-order endpoint/parameter-update limits remain in force.

### Material voice-coil self-heating and power compression

`speaker::thermal::ThermalMagneticTsDriver` binds one whole, wholly solid,
fixed-density coil cell and a sourced cooling network to the original nonlinear
magnetic speaker. The existing `MaterialNetwork` retains pinned electrical,
caloric and thermal-transport sources. Endpoint R(T), discrete-gradient magnetic
and mechanical motion, and native thermal backward Euler are solved jointly
through `cell_driver_candidate`. Every trial starts from the same accepted
magnetic and thermal history, with matching clocks and source snapshots.

`ThermalMagneticTsDriver::from_material` builds a new specimen from the same
whole-cell conductor and caloric source: its initial resistance and coil mass
are resolved together. `SpeakerMechanics` supplies only non-coil moving mass,
suspension, mechanical loss, piston area and travel; the original magnetic
profile supplies Bl and Le, and the native gas load remains separate. A data-only
electrical or density substitution therefore reaches native current, motion,
pressure and heating without a duplicated resistance or total-mass card.
Construction starts at rest at the actual thermal clock; `bind` preserves an
existing driver's state. Accepted coil mass remains fixed and is checked during
continuation. Source-swap, conductivity-equivalence, retained-state, advanced-clock
and changed-inventory regressions exercise this path.

Only the native coil-current I²R loss enters the coil cell. Suspension, compact
radiation and panel dissipation remain separately reported external losses.
Electrical work plus external heat therefore balances magnetic/mechanical
storage, finite thermal storage and the mechanical/acoustic sink without
depositing any loss twice. Source/domain refusal, cancellation, work exhaustion,
excursion failure and balance failure preserve both accepted histories.

Magnetic profiles and mechanical geometry remain fixed. General thermal coupling
is first order and requires timestep refinement. No copper coefficient, heat
capacity or cooling constant is inferred; melting, expansion, skin/eddy losses,
motion-dependent convection and temperature-dependent magnet strength require
their own models. Focused tests exercise self-heating/power compression,
constant-conductivity reduction to the original five-state cabinet, analytic
finite-store cooling and joint refusal behavior.

### Critical and overdamped acoustic radiation memory

The shared `render::plate::impact::radiation::Model` admits all finite
nonnegative pole damping ratios, including critical and overdamped acoustic
memory, with representable stiffness and viscous resistance. The same signed
positive-real impedance and reciprocal mechanical reaction remain in the impact,
coupled-modal, prepared-air and piano owners. Their generic passive/exact
propagators retain both real decay poles; no underdamped phase formula or new
time integrator is introduced. Existing modal-frequency, coupling and work
limits remain; a fast real decay is not treated as a Nyquist frequency.

The piano regressions compare actual prepared free motion and independently
integrated dissipation against critical and two-real-pole solutions, including
a fast decay far above the sample rate. Three mechanical rates recover the
analytic coupled mass/radiator response and dissipated work at second order.
The fixed fit dictionaries and their geometric/frequency qualification are
unchanged; admitting a pole regime alone does not validate a fitted exterior.

## Source-owned directional rigid friction with resolved heat partition

`friction_heat::directional::apply_source_directional_heat` closes original
directional coefficient claims, full-tensor rigid recoil and finite enthalpy
in one thermal endpoint solve. The existing rigid matrix owner supplies actual
slip and physical traction-slip work. Its source callback uses
`advance_motion_consistent` to close every declared speed/path coordinate.
Each trial starts from the same accepted body momenta and source history.

`DirectionalHeatMaterial` retains both original caloric objects and ordered
material/interface identities. Both surface temperatures follow their own
finite cells; an optional distinct Kelvin axis follows the solved common
interface temperature. Nominal pressure is the held load divided by the
actual thermal contact area. Dry contact refuses liquid surfaces.

The existing `FrictionHeatStep` derives signed body transfers from each
material conductivity, actual cell geometry and the common zero-capacity
interface. It replaces the occupied direct link only in the trial, preserving
coverage and preventing duplicate conduction. No heat fraction or additional
thermal store is introduced. Other heat/conductances retain their original
callback receipt; internal friction heat cancels once in the combined
kinetic/contact/caloric account. Only the two selected bodies, original source
history and finite thermal network publish, after independent mechanical,
local thermal, interface and total-energy gates and final cancellation.

The route retains its orthogonal prismatic, perfect-contact discretization.
Thermal cells may be local contact volumes; they need not carry whole-body
mass. Load, poses, arms and frame are held; the midpoint mechanical phase
refuses unresolved slip events. Endpoint thermal/source coupling needs time
refinement, and a nonmonotone motion-query root is not a uniqueness claim.
Five G1/G3/G4 groups cover closed-form recoil/partition, solved-temperature
weakening, full-tensor/frame-covariant motion queries, unequal-temperature
exchange without friction power, and atomic failure/cancellation/retry.

## Native rigid-body partial-slip friction with finite thermal feedback

`friction_heat::rigid::apply_rigid_friction_heat` jointly accepts one supplied
positive-load finite patch, native rigid momenta, retained regional shear/spin
memory, and the original finite enthalpy network. Every Newton trial resolves
the original same-source `PartialSlipLaw` at the common interface temperature
and calls native `apply_rigid_friction` from accepted state. Real masses,
inertias and contact arms determine mobility. Route A/B binds mechanical A/B,
with checks on original ellipse area, shared world face centre and normal.

Only original physical plastic heat enters the common thermal interface.
Elastic storage and both numerical losses remain separate in an independently
reconstructed combined energy balance. Static/kinetic coefficient continuation
retains memory and requires unchanged stiffness, regional decomposition, model
and source. Stiffness changes or any nonzero parameter work refuse because this
consumer has no corresponding physical supply port. Final publication of bodies,
contact history and thermal state follows every gate and cancellation check.

Four native physical tests cover sticking with zero heat, partial slip with only
yielded-region heat, gross 3D shear/spin versus an independent inertia/caloric
reduction, retained temperature feedback, actual momentum/work, staged refusal
and exact retry. The supported phase has one supplied loaded patch, frozen
poses/geometry and a uniform common interface. Normal/impact admission, pose
advance, moving footprints and variable-stiffness thermodynamics remain separate
owners; no instantaneous impact loss is automatically reclassified as heat.


### Permeability paths in continuing finite-gas assemblies

`MaterialExchangeState::transport_paths_step_candidate` and
`transport_paths_step_with_gas_law_candidate` accept the material owner's
`GasTransportPaths`: nozzles, molecular diffusion, porous Darcy leakage and
selective dense-membrane permeation. The original shared species/enthalpy stage
runs once before reaction and the complete electrical/thermal/moving-wall solve.
The transported composition supplies each actual gas heat-law trial. Final
transport driving-force order is checked again after heat, reaction and piston
work. The corresponding `advance_with_transport_paths` methods publish all
physical owners and histories atomically. Existing nozzle/diffusion entry points
delegate with empty additional paths and receive the same final endpoint gate.

`bernoulli_aperture::wall::material::finite_gas_dg_transport_candidate_with_paths`
previews those same paths through the existing uniform/regional source-wall
acoustic owners. Its candidate retains the complete gas and acoustic endpoint
maps plus the shared energy account; `advance_finite_gas_dg_transport_with_paths`
publishes them together. Original species enthalpy, gas-to-wave parameter work
and declared losses remain owned exactly once. All active total/partial-pressure
orders are rechecked on the final heated/acoustic endpoint. Separate finite
ducts do not acquire an invented common pressure field, acoustic continuity,
jet force or loss-to-wall heat partition.

These adapters retain each constitutive source domain and the original explicit
first-order split. A rejected final ordering requires a smaller interval or a
more suitable physical model, not per-edge flux clipping. Added real-fixture
regressions exercise porous transport in the moving assembly and selective
membrane transport through continuing acoustic owners, including immutable
preview, atomic publication and enclosing species/energy accounts.

## Conduction-derived heat for retained thermoelastic partial slip

`tribothermal::partial_slip::resolved` adds an explicit resolved-interface
route to the original `CaloricPartialSlip` owner. Its constructor needs no
prescribed body heat fraction. Existing loaded material-bound contacts can
upgrade their routing without resetting regional strains, sliding branches,
accepted time or pinned friction claims. The original prescribed route remains
available; a resolved route cannot accidentally use that route's step.

The existing finite thermal contact supplies matching orthogonal geometry and
half-cell distances. Same-endpoint original body conductivities determine a
common massless interface temperature and signed heat flows. Only the native
contact's physical microslip loss enters that node. Each body separately retains
its reversible thermoelastic exchange `-T_i Delta s_i` from the original
two-temperature affine free-energy model. Numerical return-map loss stays
outside heat. Reports distinguish conductivity-derived generated-heat
allocation, signed ordinary conduction, reversible exchange and the independent
thermal/interface/whole-system residuals.

Pinned material queries retain their declared body-temperature axes, pressure,
actual translational speed, source domains and loaded memory. The common node
does not fabricate a different temperature axis for an existing claim. Legacy
partition metadata in an upgraded material binding is ignored for heat routing
and replaced by the derived allocation in accepted reports.

Six new native physical tests and nine original retained-contact/material
regressions pass. Independent references cover two-capacity conduction,
reversible entropy/enthalpy exchange, zero-slip hot-to-cold transfer, and
temperature-dependent Coulomb feedback. Source-domain, competing heat-owner,
cancellation and other refused trials leave accepted contact/thermal histories
unchanged. Scope remains a fixed normal patch/load, prescribed relative motion,
matching finite-volume geometry, negligible interface resistance/storage and
first-order endpoint coupling.

### Voltage-driven hysteretic cores with distinct finite heating

`material_circuit::hysteretic::HystereticCircuit` replaces one original
inductor branch by `UniformHystereticCore` through a private native current port.
The existing Newton solve closes linkage change against native terminal voltage
using the exact work-conjugate mean current evaluated from frozen accepted
magnetic memory. One original native DAE assembly serves all scalar magnetic
trials. Ordinary L/C/source branches and material-derived graph resistance retain
their original owners. Private-port magnetic work cancels exactly once; core
pinning heat, sourced graph/copper heat and external work remain separate.
No numerical residual or fitted core-loss map is deposited as physical heat.

The `thermal::ThermalHystereticCircuit` child uses the prepared `MaterialNetwork`
and existing `cell_driver_candidate`, with distinct whole-core and winding
receivers, sourced R(T), fixed core geometry/density and an explicit admitted
core temperature interval for the unchanged magnetic parameters. It checks
source snapshots, volume, clocks, heat deposition and the combined energy
account before returning both candidate histories. Refusal publishes neither.

An increment assumes monotone linkage. Actual magnetic reversals and source
jumps require interval cuts. Start/midpoint/endpoint native voltage probes
refuse detectable turns; they do not certify unresolved interior oscillations.
The linear limit matches native midpoint stepping; a smooth active branch has
second-order time convergence. Pinning/saturation crossings require refinement,
and the general coupled thermal scheme is first order. Current-only constraint
topologies may refuse the admitted private-port native DAE. No vector hysteresis,
eddy/excess loss, frequency-dependent fit or invented material-temperature
coefficient is introduced.

Six original-source physical tests pass: native linear RLC reduction, independent
major-loop pinning area/remanence, analytic active-set RL refinement, reversal/
source/budget/cancellation refusal, distinct finite copper/core heating with R(T)
feedback and global closure, and thermal-domain rollback.


### Reciprocal finite-caloric normal contact

`tribothermal::normal::CaloricNormalContact` advances one prescribed linear indentation segment and the original finite thermal network together. The native sphere or elliptic Hertz/Hunt-Crossley owner supplies force, patch geometry, pressure, free energy and applicability. Normal load and patch evolve with indentation. The segment mean force, rather than endpoint force multiplied by travel, is conjugate to actual normal displacement. Hunt-Crossley work is integrated over that segment and deposited once through the declared heat partition.

The constitutive model supplies the reduced modulus and its two-temperature gradient/Hessian at fixed geometry. Surface entropy is `s_i=-dpsi/dT_i`; contact internal energy is `U=psi+sum(T_i*s_i)`. A symmetric coordinate-path split accounts for nonlinear coefficient changes; reversible exchange is `-C_i-delta(T_i*s_i)`. The original caloric capacities plus the full contact thermal Hessian must admit a positive-definite two-temperature stability matrix. Every trial starts from accepted normal and caloric state. Original caloric objects, cell masses, source domains and native contact applicability remain bound. Independent native-potential, thermal-cell and whole-system energy gates precede atomic publication; failure/cancellation leaves both owners unchanged.

Geometry/curvature and declared damping remain fixed over the segment. This first-order composition prescribes normal motion; it does not implement flight/impact event localization, thermal expansion, phase change, moving curvature, simultaneous tangential/body motion or complete Euler thermal closure. Seven focused native G1/G3 tests cover analytic reciprocal cooling, loaded heating, integrated damping on compression/release, circular degeneration and tiny increments, first-order refinement, capacity refusal and transactional failure/replay. Sources are authored references, not measured calibration.


### Evolving normal contact with directional friction and finite heat

`tribothermal::normal::directional::NormalDirectionalContact` joins the original
Hertz/Hunt--Crossley normal owner to the source-owned associated Coulomb ellipse
inside the same finite-enthalpy solve. Each thermal trial uses the normal
owner's actual endpoint load and circular or elliptic patch area. Their mean
pressure, the two original surface temperatures, declared mean temperature,
prescribed sliding speed and unsigned path refresh the same selected friction
claims. Normal reversible heat/storage and physical frictional dissipation
remain separate.

The original normal callback supplies its already evaluated patch and work;
there is no second normal or thermal solver. Actual traction--slip work is
partitioned explicitly between the two original caloric cells and a named
outside destination. The outer account excludes that internal friction-heat
exchange from external supply and independently balances normal internal
energy, finite heat, normal/tangential mechanical work and exported work.
Normal, friction, path and caloric candidates publish together; source-domain,
work-budget or cancellation failures leave all accepted owners unchanged.

Indentation and actual tangent displacement are prescribed. Tangential load
and pressure use the endpoint normal patch, so changing load/material response
requires time refinement and has first-order temporal accuracy. Curvature and
the tangent/material frame stay fixed. This does not locate contact/separation
events or solve simultaneous normal/tangential inertia, rotation or moving
footprints. Five focused G1/G3 groups exercise the held-patch limit, actual
Hertz area/pressure feedback, thermal/normal reciprocity, continuous-work
refinement, rotated path continuation and transactional refusal.

### Original material sources in reciprocal normal contact

`tribothermal::normal::material::MaterialNormalModulus` binds the original `ThermalHertzPair` to `CaloricNormalContact`. Every thermal iterate re-queries each solid at its own temperature and uses the original nonlinear reduced modulus, gradient and Hessian in native normal response and reversible heat. Both ordered caloric charts must share the exact original material cards, external query context and density claims; the normal owner then retains those exact caloric objects. The reference normal request must match the native reduced modulus and limiting yield strength. No representative bulk temperature, reset of normal storage, coefficient-work heat credit or changed reference mass is introduced.

Three focused native G1/G3 integrations, run alongside the seven original normal cases, verify source-backed loading against an independent rational-compliance caloric root, heating at held indentation with changing contact INTERNAL energy, and original-card/context/density/domain refusal with exact replay. In the authored loading reference, initially 300 K / 360 K bodies reach approximately 288.238949147 K / 360 K; mechanical work is 6.832367020 J, contact internal energy increases 8.008472106 J and physical damping heat is zero. The difference is reversible caloric cooling. These are synthetic implementation references. The original affine per-solid E, fixed nu/yield, fixed geometry, finite source-domain and first-order coupling limits remain in force.


### Inertial normal motion with reciprocal finite heat

`tribothermal::normal::dynamics::InertialNormalPair` determines normal compression from two complete finite-body inertias, constant exterior forces and the original caloric normal-contact owner. A single endpoint indentation unknown closes actual midpoint relative travel. The same native segment-mean force supplies equal/opposite impulses and mechanical work, while the inner normal/thermal solve evolves the Hertz patch, reciprocal source stiffness and Hunt-Crossley heat. Original masses and accepted caloric inventory remain fixed.

Independent momentum, midpoint geometry, native force/work and `delta(K+U_contact+H)-W_external-Q_external` gates precede atomic publication. The discrete-gradient normal work creates no artificial kinetic heat. Five native reference groups verify an independent Hertz midpoint root with unequal masses/drives, conservative compression/rebound and time reversal, integrated damping heat, thermal feedback changing actual motion, and source/budget/cancellation/separation refusal with replay.

This is the fixed-normal, nonnegative-indentation contact branch. An interval crossing separation or original no-tension/source bounds refuses; callers must shorten or locate that event. No free-flight/impact event locator, changing normal direction, tangential coupling, moving thermal footprint or full Euler mechanics is claimed. Complete-body single-cell caloric models need justified uniform-temperature applicability. The inherited thermal composition is first order; zero-slope undamped reversibility does not promote the general thermomechanical scheme to second order. References are authored implementation fixtures, not measured validation.


### Rolling and retained microslip in one finite thermal solve

`tribothermal::partial_slip::rolling::CaloricRollingContact` composes original contour/rolling loss and retained partial slip through one original finite-enthalpy endpoint solve. Original rolling law/state/checkpoint, complete patch/interface identity, independent generalized-work coordinate, explicitly selected actual surface temperature and heat split are retained. Native rolling work is independently reconstructed as `-F_contour*v*dt-M_roll*omega*dt` and deposited once; that same endpoint immediately updates source friction and reciprocal contact stiffness. A combined `delta(H+U_contact)+D_numerical-W_microslip-W_rolling-Q_other` gate precedes atomic publication.

Five focused native reference groups cover all three original loss alternatives, unequal-temperature conduction, loaded elastic thermal response, immediate source mu(T) feedback, separate torsional work, rollback and replay. The normal patch is fixed; rolling kinematics are prescribed, with original constant coefficients and their domains. No mean surface temperature or new empirical law is invented. Endpoint composition is first order; synthetic fixtures do not establish measured calibration or full Euler coupling.

This wrapper uses the prescribed-partition contact entrypoint. The separate conductance-resolved routing API is preserved; it is not automatically substituted for either declared rolling or microslip partition.


### Inertial motion with retained finite-caloric microslip

`tribothermal::partial_slip::dynamics::InertialContactPair` advances two translating bodies in a fixed contact tangent plane. Its two outer unknowns are contact impulses [N s]. Every trial reuses the original loaded partial-slip/finite-enthalpy candidate at actual midpoint relative travel, preserving strains, kinetic branches, pinned friction sources and thermal history. Equal and opposite impulses update both momenta; exterior forces do work on each actual displacement. Midpoint kinetic work introduces no fictitious kinetic heat. Native backward-Euler spring dissipation remains a separate numerical account.

Each moving body occupies one complete finite caloric cell. Optional original `ResolvedBody` bindings retain source density/geometry and must match ordered material, total mass and volume; heating never replaces mass or momentum. Independent momentum, native wrench/travel, kinetic-work and full `delta(K+U_contact+H)+D_numerical-W_external-Q_external` gates precede atomic publication. Late source refusal, work exhaustion and cancellation preserve both histories.

The normal patch/load and tangent frame remain fixed. Rotations, moving thermal footprints and changing normal geometry are outside this rung. Elastic sticking retains finite shear compliance and can have relative motion; event-exact Coulomb arrest is not claimed. Endpoint contact/thermal composition is first order and requires refinement. Five focused analytical reference groups cover planar inertia/force work, source mass and mu(T) feedback, reciprocal stiffness/heat/inertia, loaded elastic reversal with 24/48/96-step refinement, and refusal/replay. Sources are authored references, not measured validation.

This motion wrapper uses the prescribed-partition contact entrypoint. The separate conductance-resolved route remains available through its original API and is not implicitly substituted here.


### Simultaneous retained rigid friction and finite thermal interfaces

`friction_heat::network::FrictionHeatNetworkStep` assembles one backward-Euler residual over the original cell energy increments and one massless common-temperature coordinate per distinct original contact link. It preserves original geometry, material, mass and surface coverage, suppresses each replaced direct series link once, and accumulates signed interface powers into shared cells in canonical link-ID order. No interface capacity or prescribed heat fraction is introduced. Per-interface balance, cell balance, the original thermal audit and total energy have independent admission gates.

`friction_heat::rigid::network::apply_rigid_friction_heat_network` resolves every original contact law at that same candidate thermal endpoint, calls the native rigid multicontact friction owner once, and routes each regional return map's physical microslip heat to its own interface. Shared-body inertia and shared-cell thermal feedback are simultaneous. Bodies, retained regional history and the original thermal simulation publish together after all gates; source refusal, exhaustion and cancellation preserve accepted state. Elastic storage, elastic/mechanical numerical loss and outside heat retain separate accounts. Nonzero coefficient work is refused.

The scope remains frozen poses, supplied positive normal patches, matching finite elliptical contact and orthogonal thermal faces, negligible interface resistance/storage, and uniform temperature within each patch. Temperature continuation may change static/kinetic Coulomb coefficients only; stiffness, region decomposition and source/model identities remain fixed. Normal solving, moving footprints, impact-loss partition and variable-stiffness thermoelastic work remain separate owners.

Four G0/G1/G5 tests in `rigid_friction_heat_network.rs` cover an independently eliminated three-body/two-interface reference with shared middle-body mobility and shared finite heat capacity, retained second-step input permutation, source/budget/cancellation/competing-owner rollback and exact retry, and zero-friction hot-to-cold conduction. The reference interface temperatures are 304.69405993744004 K and 303.18187043010846 K; physical heats are 0.17343547971344292 J and 0.088803662273099 J. Tests admit 2e-7 K and 1e-9 J respectively. Four existing single-interface regressions also pass in the bounded native source-selected lane; this is not full workspace/DSR verification.


### Prepared articulated contacts and same-step finite-cell thermal feedback

`friction_heat::articulated::propose_articulated_friction_heat` borrows the
accepted native articulated systems, original framed constitutive checkpoints
and original enthalpy simulation. It returns one complete candidate containing
every body response, every new framed history and the thermal endpoint; the
surrounding integrator publishes these together. No accepted body or thermal
owner is mutated during preparation, root iterations or final admission.

The original `PreparedArticulatedFramedTangentialContacts` transports histories
once and retains the full native ABA contact response and constitutive factor.
Thermal trials reuse those bodies, histories and immutable laws while solving
each distinct contact temperature. Dense contact matrices stay borrowed until
the final report. Each native irreversible heat enters its explicit original
route in one `FrictionHeatNetworkStep`; shared cells retain one capacity, signed
ordinary conduction remains active, and replaced direct links are suppressed
once. Accepted-state shear energy, native recoil and prescribed-support work
remain in the combined balance. Optional final-root entropy uses the original
caloric integrals without another physical solve.

Each mechanical site must have one distinct heat route. Callback/report order
is canonical ascending original link ID, with an explicit map back to original
mechanical order. Route A belongs to the point, B to its dynamic counterpart or
prescribed support; the caller owns material/body identity. The prepared common
world point must meet the actual shared face centre and the tangent normal must
be parallel to both cell centrelines. Contact temperatures are solved, so
prepopulated drive temperatures refuse. Original source laws own temperature,
pressure and history applicability; source callbacks cannot replace laws.
Frozen elastic laws retain their refusal of unsupported temperature changes;
this adapter supplies no variable-potential thermoelastic heat/entropy port.

Cell/contact/unknown counts, trials and nonlinear iterations are bounded.
Source identity, original per-contact work, native momentum/energy, cell,
interface, thermal, combined-energy and cancellation gates precede return.
`tests/articulated_friction_heat.rs` covers an independently eliminated joint
thermal-softening/rigid-recoil reference on two contacts sharing a body and a
cell, original pinned coefficient receipts, world/frame covariance, zero-slip
conduction/entropy, reversible storage release and frozen-temperature refusal,
and whole-candidate late source/cancellation/budget retry. Its curves are
synthetic fixtures, not measured-material validation. Geometry, normal loading
and material spin remain caller-owned frozen-interval inputs; spatial/time
refinement and source applicability remain necessary.

### Source-backed Hertz recoil and finite-cell thermal feedback

`friction_heat::hertz::propose_rigid_hertz_heat` borrows the original rigid
bodies, their mass properties, a pinned `SourceThermalCoulomb`, and the original
finite enthalpy network. It returns both contacted body candidates, the complete
thermal endpoint, the original Hertz work report, and final coefficient-use
receipts together. Accepted mechanics, source ownership and heat remain
unchanged through preparation, nonlinear trials and final admission.

One `PreparedRigidHertzPair` retains the full point/spin mobility from both
native inertia tensors. Each thermal trial refreshes only the original selected
static/kinetic friction claims at its solved common-face temperature and explicit
caller-supplied pressure, then solves the original joint Hertz force/spin law.
Normal load, circular radius, poses and contact frame stay held. Unsupported
additional evolving query coordinates refuse through the original source owner.
The mechanics and thermal network share the actual representable interval
`(accepted_time_s + dt_s) - accepted_time_s`.

Actual native dissipated work enters one original contact route through
`FrictionHeatNetworkStep`. Endpoint material conductivities and the original
cell-centre distances determine signed interface transport; no heat fraction or
new capacity is supplied. Other links, boundaries and signed outside heat stay
active. The replaced direct contact link is suppressed once. The original
circular face, shared centre and normal must agree with the mechanical patch.
Caller-owned body/material correspondence remains explicit. Per-cell,
interface, thermal, native momentum/work and combined kinetic-plus-caloric
energy gates precede the returned candidate. Optional entropy uses the original
caloric law at the final root.

Five new native integration groups in the existing
`tests/articulated_friction_heat.rs` target cover independent scalar thermal
references for translation and spin, both-body recoil, source receipts,
zero-slip ordinary conduction and entropy, physical geometry/clock/source
admission, and late source/cancellation/budget refusal with exact retry. The
clock case independently checks recoil and heat at a nonzero clock where the
requested and representable durations differ. All five and the six existing
articulated groups pass against the current original-source full-inertia
implementation. These are synthetic analytical fixtures, not measured-material
validation or a full Cargo/DSR workspace gate.

This is the stateless static/gross-sliding/spin Hertz rung; it does not consume
elastic partial-slip history or rolling-deformation work. Normal evolution,
moving footprints and pose integration remain separate owners. Finite cells,
a common patch temperature and backward-Euler heat require spatial/time
refinement; nonlinear convergence alone does not prove root uniqueness.

### Finite-body rolling with original whole-sphere caloric storage

`friction_heat::rolling::propose_rigid_rolling_heat` stages the original native
two-body rolling response, retained rolling history, and finite enthalpy endpoint
as one complete candidate. The original deforming sphere's viscoelastic moment
work enters its own bulk caloric cell exactly once. It is not routed to a
zero-capacity contact interface or divided by an authored heat fraction.
Ordinary network links, boundaries and signed other heat remain active.

The explicit binding matches the native deforming A/B body to an existing
thermal body/material ID and a cell whose mass equals the whole sphere mass.
A second active cell with that thermal body ID refuses, so this one-cell model
cannot duplicate the sphere's caloric inventory. Unborn scheduled cells have no
accepted energy; original event-boundary and activation rules remain in charge.
The caller declares correspondence of the rolling card's material because native
`MassProperties` contains no material ID to compare independently.

The original rolling card has constant coefficients within its authored
temperature/frequency domain. One cached native solve supplies actual work while
every thermal trial admits its actual sphere temperature through that same card.
Final mechanics and history are restaged at the converged endpoint temperature;
exact physical agreement with the cached state, work and impulse is required.
Final source-use receipts are retained. Initial card temperature must equal the
accepted sphere temperature. Both owners share the actual representable clock
interval. No temperature-dependent viscosity or material refresh is invented.

The inner thermal audit includes injected rolling work; a separate outer account
retains only authored/caller outside heat. Actual caloric change minus outside
heat minus native rolling dissipation, and actual contacted-pair kinetic change
plus caloric change minus outside heat, have independent joule gates. Other
bodies are neither copied nor included in pair energy. Accepted states remain
borrowed through source failures, nonlinear budget exhaustion and cancellation.

Seven native groups in `tests/rolling_friction_heat.rs` pass against the current
original-source full-inertia body implementation. They cover independently
derived two-axis recoil and finite-cell conduction, latent enthalpy, signed
other heat and endpoint conductance, zero-speed/viscosity/load conduction,
deforming surface B and remapped body indices with an energetic spectator,
whole-sphere/source/temperature admission including duplicate active inventory,
and loaded-history late failure/cancellation/budget with exact retry. Synthetic
fixtures establish analytical behavior, not measured-material validation or a
full Cargo/DSR workspace gate.

This is uniform whole-sphere bulk temperature with the original held normal
patch and no-slip rolling law. Spatial loss density, flash heating, sliding or
spin loss, thermoelastic storage and temperature-dependent rheology remain
separate physics. Endpoint heat transport is first order and needs time/spatial
refinement independently of nonlinear convergence.

### Simultaneous Hertz surface heat and sphere rolling bulk heat

`friction_heat::hertz_rolling::propose_rigid_hertz_rolling_heat` combines the
original five-coordinate native contact response with one original finite
enthalpy network. It returns both body endpoints, the original rolling history,
the thermal endpoint and final source-use receipts together. All accepted
owners stay borrowed through preparation, nonlinear trials and admission.

One prepared response is reused throughout. The original pinned friction
claims resolve at the solved common interface temperature. The original
constant-coefficient viscoelastic card instead uses the actual sphere bulk
temperature; its before/midpoint/endpoint, source, frequency and normal-patch
limits remain active. Both physical losses change together with native recoil,
including rolling induced by surface traction from an initially zero rate.
The original patch supplies the load and circular area/radius; nominal pressure
remains an explicitly supplied original source coordinate.

Hertz surface loss enters the existing massless interface through
`FrictionHeatNetworkStep`; endpoint conductivities and original distances
determine signed heat into both bodies. Bulk rolling loss enters the explicitly
bound whole-sphere cell once. That cell must occupy the correct A/B route side,
match native mass and declared thermal body/material IDs, and have its centre
at the actual sphere COM. Duplicate active sphere inventory refuses. A thin
surface cell cannot stand in for the whole sphere's heat capacity. The shared
face, normal and circular radius retain the original Hertz geometry checks.
Other caller-owned body/material correspondence remains explicit.

Outside heat is measured before injecting bulk work. The inner interface report
therefore retains an explicitly documented bulk-source contribution in its
external account; the returned combined candidate separately exposes actual
outside heat. Independent gates check
`delta H-Qoutside-Dsurface-Dbulk` and
`delta K_pair+delta H-Qoutside`, alongside the unchanged native work/momentum,
cell, interface and thermal gates. No residual or numerical damping becomes
heat. Mechanics and thermal transport share the actual representable
endpoint-minus-start interval. Final source/temperature/mechanics/history
receipts all come from the final restaged root.

Five new native groups in the existing articulated thermal target use
independently derived planar force/torque equations and an eliminated finite
caloric network to bracket the source-temperature root. They check both losses
and recoil, induced rolling, distinct interface/bulk temperatures and receipts,
signed outside heat, zero-motion/viscosity conduction, original geometry and
whole-sphere ownership, source/clock/effort/reversal admission, and loaded late
source/cancellation/budget with exact retry. All five plus the eleven existing
Hertz/articulated groups pass against the original full-inertia native build.

The native normal/shear decoupling and held-geometry limits remain unchanged.
Uniform sphere temperature and a common interface temperature are finite-cell
approximations; backward-Euler transport requires independent refinement.
Released zero-area patches, spatial loss density, variable viscosity,
thermoelastic storage and entropy classification are not supplied by this
adapter. Analytical synthetic fixtures are not measured validation or a full
Cargo/DSR workspace gate.

## Finite-gas rocket propulsion and translation

`rocket_blowdown::GasRocketState` binds an immutable source-bearing finite
reservoir to the native `fs_mbd::rocket::RocketState`. Propellant mass is derived
from the reservoir, while dry mass includes the tank/container. Callers cannot
author or mutate the two inventories independently. Each step accepts the
original material owner's source-resolved adiabatic blowdown and derives the
mechanical forcing from its accepted extensive quantities:
`mass_flow = expelled_mass / duration` and
`effective_exhaust_velocity = axial_impulse / expelled_mass` along an explicit
finite unit thrust direction in the matching caller-declared inertial frame.
The frame identity is an admission check, not a coordinate transform or proof
that the caller's frame is inertial. Endpoint rate resampling cannot substitute
for the accepted mass and impulse, and nozzle pressure thrust is included once.

The adapter retains the original gas, source/nozzle and dynamics reports. Both
candidate endpoints publish together only after remaining/consumed propellant
mass checks, a directed thrust-impulse check, and a final cancellation poll.
Invalid state or forcing, source-domain failures, unsupported nozzle regimes,
failed balances and cancellation publish neither endpoint. The existing gas
mass tolerance bounds cross-owner discrepancies; the actual residuals remain
visible rather than being repaired or hidden.

The remaining reservoir follows the exact frozen-composition adiabatic
thermodynamic path. The finite-time mass-flow/thrust quadrature is second order
for smooth solutions; native mechanics integrates the resulting step-held
forcing analytically. Refinement tests cover both velocity and position. This
is comoving, spatially mixed, fixed-volume gas storage on a point vehicle, not
an attitude, slosh, combustion, real-gas, wall-heating, rotating-frame or full
world-frame exhaust-energy model. The `rocket_flight` example runs the complete
source -> nozzle -> finite tank -> variable-mass translation journey using the
existing explicitly manufactured gas-source fixture and emits time-resolved
mass, thrust, position, velocity and coupling residuals.

## Signed adhesive normal storage and native approach dynamics

`BoundNormalContact` and `BoundNormalStorage` accept the original conservative
DMT source through `NormalLawSource`, preserve negative geometric approach, and
project signed force through the same work-conjugate mechanical map. Finite-gap
attraction can accelerate bodies before solid contact. Unilateral Hertz and
finite-gap sources keep their unloaded-origin and clipping behavior. DMT
attraction adds no internal damping and is not frictional heat.

The existing native Gonzalez integration consumes a cancellation-stable complete
storage numerator when every constituent supplies one (or the inner owner
explicitly declares an exact quadratic matrix). Missing constituent support
returns `None`; no partial numerator is substituted. Affine-map and supplied
midpoint roundoff corrections are included. Source validity remains active at
every trial, including corner crossings. This is frozen-source elastic
adhesive dynamics; evolving adhesion, radiation, plasticity, and JKR snap
energy require their explicit physical owners.

Focused regressions cover the real pHS reference and prepared steppers,
acceleration from rest at finite gap, contact onset, equal/opposite momentum,
energy conservation, and coordinate translation. Native finite-difference
probes include the current residual scale so internal force at a zero state
remains numerically resolvable; solver tolerances and budgets are unchanged.

## Adhesive junctions on the native modal network

`friction::adhesive::AdhesiveModalBinding` attaches an explicit real area,
repulsive normal load and complete junction shear law to `FrictionModalSystem`.
An optional `SourceAdhesiveTangential` retains the original normal and shear
owners. The existing full finite-body compliance selects static reactions or
compliant sliding; the existing native propagation, state publication, observer
and energy transaction remain in charge. The actual applied force times solved
slip supplies friction loss. Equal/opposite body reactions retain momentum,
and zero repulsion with real contact area does not release the junction.

The `FrictionDrive` normal input must equal the repulsive force held over the
current tangential sample. It is
not the signed net adhesive reaction. Coulomb and affine pressure contributions
use the original held-effort solve. Stribeck and hyperbolic pressure contributions
use the same cached `StribeckLoadLine` continuation as dry friction, with actual
`tau_kinetic A` added to sliding force and independent `tau_static A` in the
static cone. This offset changes the residual without changing its critical
slips; no effective normal load or division by repulsion is needed. A junction
strength drop can create four intersections (the cone plus three sliding roots).
The retained low/high sliding branch persists inside the static cone until its
root ends; a finite-slip fold jumps branches, whereas unloading to zero slip or
reversal allows sticking. Only accepted native steps publish branch history.
Binding adhesion invalidates a former dry load line without erasing body motion.

`step_with_adhesive_contact` and its cancellation-aware variant accept updated
area/load from the original normal owner in the SAME native transaction.
`AdhesiveModalContactUpdate::Source` requires a response-derived
`SourceAdhesiveContactUpdate` with the original normal identity; `Declared`
geometry is accepted only for an originally unsourced binding. Original material
claims, selected strengths, fixed query/temperature and source allocations are
retained. Only the small constitutive leaf and load-line cache are staged.
Changed area invalidates the cached junction offset even at unchanged load.
An attached-to-attached update preserves the occupied low/high branch; an
accepted zero-area detach or subsequent reattachment resets it. A caller must
resolve and submit separation events rather than skip detached intervals.
Actual current geometry enters both scalar solve and accepted-force check.
Late source, constitutive, mechanical, energy or cancellation refusal restores
the previous geometry, load-line cache and branch along with native mechanics.
Normal evolution and work remain external; the held geometry update does not
add a normal energy account, shear peeling, or a temporal accuracy claim.

G1/G3 load-line tests compare all four intersections with an independent physical
residual scan, test fold hysteresis and zero-repulsion sliding, and recover dry
continuation exactly at zero junction strength. `tests/adhesive_modal.rs` exercises
both weakening laws on two real finite bodies, including static-cone overlap,
equal/opposite momentum, actual friction-work balance and failed-trial rollback.
The external normal owner supplies area evolution; mixed-mode peeling and tangential
elastic storage are not inferred. Dry source replacement, simultaneous-port
promotion, and all dry-law thermal refresh entry points refuse to discard or
reinterpret an attached adhesive source. Friction loss remains available for
explicit external heat routing; an evolving thermal adhesive model is not
claimed by this fixed-query binding.

Cold internal-coordinate extension preserves the same adhesive binding alongside
existing source interfaces and nonlinear storage.

## Native mixed-mode cohesive debonding

`contact::cohesive::CohesiveModalSystem` explicitly bonds the current sample-zero
modal configuration at zero local jump. Three complete work columns map opening
and both shear components; their transposes apply the same physical forces to
every retained body coordinate. All nine held-force compliance entries come
from the original `CoupledModalSystem` response, including bilateral springs,
dashpots and admitted material/acoustic memory. No diagonal compliance surrogate,
prescribed reaction history or replacement mechanical propagator is introduced.
The original area and fixed local frame are caller declarations. The optional
`from_source` constructor checks the two ordered receiving material identities
and borrows the complete original `ResolvedBilinearCohesive` and its receipts at
fixed external conditions. Bond creation, changing texture frames, source
refresh and rebonding after rupture are not inferred from geometry.

`CohesiveLoadLine` solves the complete jump using the material owner's analytic
three-component path mean. The native wrapper stages the original mechanics
under those held reciprocal forces, recomputes the actual three-component jump,
and independently reevaluates the timed material path on that motion. Separate
jump and componentwise absolute-plus-relative force gates precede the whole
energy gate. Network and recoverable interface storage changes plus original
network loss, intrinsic fracture and physical surviving-bond viscous loss must
equal authored modal-force work plus applied force dotted with prescribed jump
motion. The original total-energy ceiling remains active. One material
evaluation is reserved for the final actual-state audit within the supplied
evaluation budget. Only after all gates and cancellation checks pass are native
motion, irreversible cohesive history and the complete frame published together.
No mechanical or source owner is cloned per step.

The original common-mode or frozen-onset BK law owns initiation, irreversible
damage, complete failure and undamaged compression on reclosure. Endpoint
traction is reported separately from the mean force that performs body work.
Fracture and viscosity are distinct losses; neither is silently relabeled as
frictional heat. Local softening can admit multiple roots, so bounded Newton
convergence does not certify uniqueness or resolve unreported intrastep path
reversals. Physical timestep convergence and fixed-frame applicability remain
the caller's responsibility. `tests/cohesive_modal.rs` exercises original
six-coordinate finite bodies, three momenta and complete energy, a native
oblique spring with off-diagonal compliance, prescribed-port work, mixed rupture
and reclosure, source identity, and late-refusal replay.

## Material magnetostrictive rod inside a Kirchhoff RLC circuit

`material_circuit::magnetostrictive::MagnetostrictiveCircuit` binds one original
inductor branch to `magnetostrictive::MagnetostrictiveRod`. The selected numeric
inductance is replaced by the source-derived clamped inductance; its uncoupled
flux energy is removed before inserting the original rod's coupled strain/flux
and kinetic energy. Native charges, other fluxes, source indices (including
gaps), resistor stamps and transformer constraints keep their original owner.
Appending the actual extension and momentum coordinates to that descriptor
gives reciprocal actuation and motional voltage from the same material
derivatives. No mechanical time integrator, independent sensing gain or fitted
damping is introduced.

The complete network uses perturbations about the rod's explicitly held
DC-biased, stress-free reference. The original material snapshot, reference
geometry, uniform magnetic path and caller-declared strain/H envelope stay
fixed. Electrical and mechanical source work close against combined incremental
storage and source-resolved `R*(Delta i)^2` resistance dissipation. The reported
branch records are explicitly incremental: this loss is neither the total
`R*(I_bias+Delta i)^2` winding heat nor generally its change from a DC baseline.
Bias preparation/supply and cross terms are outside this account; no thermal
deposition or temperature-dependent magnetic law is inferred.

The original bounded descriptor and discrete-gradient stepper own evolution.
For this fixed quadratic tangent and held inputs the method reduces to
second-order midpoint; source discontinuities require interval cuts and motion
requires time refinement. Candidates preserve both histories on source/domain,
constraint, work or cancellation refusal. The same source energy supplies an
exact quadratic Hessian to the native residual roundoff bound. Small energy
residuals do not certify trajectory accuracy.

Focused synthetic regressions compare an independent coupled RL/mechanical
midpoint solution, source-free mechanical-to-resistive energy transfer, a
continuous damped eigenmode and second-order refinement, zero-coupling reduction
to an interleaved native RLC graph, reversed coil orientation, and true coupled
current clamps with rollback. This extends MR60 with loaded axial motion and
sensing; it does not close its spatial-wave, hysteretic, thermal, acoustic-pressure
or experimental-validation scope.

### Native rigid partial slip with reversible thermal stiffness feedback

`tribothermal::partial_slip::CaloricPartialSlip::step_rigid_resolved_candidate` closes original native rigid shear/spin impulses around the existing conduction-resolved caloric partial-slip owner. Its solve/report/candidate/error types live in `tribothermal::partial_slip::rigid`. Actual material arms, poses, masses and principal inertias determine midpoint slip and spin. Each thermal trial retains the original two-temperature affine contact free energy, regional strain history, pinned friction claims and finite caloric network. Geometry must match the original elliptical patch and shared thermal face. Partition follows original conductivity/geometry; no caller heat split or effective replacement body is introduced.

Independent native point-wrench, linear/angular momentum, kinetic work and `Delta(K+U_contact+H_caloric)+D_numerical-Q_external` checks gate one atomic candidate. Physical microslip heat, reversible per-surface entropy exchange, stored contact energy and numerical loss remain separate. A translational material source solves only its two admitted shear coordinates; genuine spin still reaches and fails the original applicability gate. Fixed poses and positive supplied normal load remain explicit; the designated caloric cells need not represent either complete rigid body. This does not provide normal contact evolution, moving footprints, nonlinear stiffness thermodynamics or measured calibration.

Four native physical groups pass: loaded two-axis shear/spin with variable stiffness against an independent two-temperature entropy solution; gross sliding with original arm and spin moments and heat/storage/loss separation; pinned mu(Ta) with affine stiffness against an independent native/caloric reference (including a deliberately reversed legacy heat-fraction hint that cannot override the derived 3:1 partition); and late source/budget/cancellation/geometry refusal with exact retry. The source fixture also verifies actual-spin refusal and retained pinned claims/history. These are focused source-selected original-owner builds with the pinned nightly, not full workspace Cargo/DSR verification.


## Air-relative rocket loads and explicit gravity

`rocket_blowdown::environment` consumes the existing coherent `GasState`
and native `fs_mbd` quaternion/vector owners. `FlightAirSample` binds the
ambient pressure, thermodynamic/transport snapshot and world wind velocity
to a declared inertial frame. Flow at a body-fixed application point uses
`v_center + R(omega cross arm) - v_wind`, including rotational point velocity.
Mach number, Reynolds number and dynamic pressure are derived from that same
sample and retained immutably. The explicit drag law applies `-q Cd A v/|v|`
and its matching moment about the center of mass. Its relative mechanical
power is nonpositive; this does not allocate dissipation to vehicle heating.
A common Galilean velocity boost leaves the air-relative load unchanged.

The caller supplies the coefficient, area and reference length and remains
responsible for the coefficient's Mach/Reynolds/attitude/configuration domain.
No hidden coefficient, lift law, stability derivative or rarefied-flow model
is supplied. Source authority remains with the original atmospheric/property
producer. `CentralGravity` separately accepts an explicit positive G M,
source center and frame, returning consistent acceleration and specific
potential; it does not infer Earth rotation, J2 or a terrestrial constant.
Frame mismatches, invalid data and unrepresentable loads refuse. Focused
regressions cover co-motion, a common velocity boost, rotational drag work,
inverse-square gravity and the potential gradient. These are instantaneous
loads; trajectory accuracy remains the time integrator's responsibility.


## Finite source gas coupled to native rigid rocket flight

`rocket_blowdown::rigid::GasRigidRocketState` binds one immutable source-bearing
`ClosedGasParcel` to the original `fs_mbd::rocket::rigid` vehicle. All propellant
mass is derived from the actual tank; callers cannot author a second inventory.
A powered interval advances the original gas owner once, derives `mdot` and
effective body exhaust velocity from accepted expelled mass and scalar axial
impulse, and drives the original native dynamics with those quantities.

The adapter checks actual remaining and consumed mass and native bookkeeping
against the gas owner's explicit mass tolerance. Scalar axial mechanical impulse
is checked against nozzle impulse independently of directed world impulse,
whose norm can be smaller as thrust rotates. Pressure thrust contributes once
to total force and its application-arm torque. Original source/nozzle reports,
angular-momentum accounting and actual cross-owner residuals are retained;
cancellation or either owner's refusal publishes neither endpoint.

`coast` explicitly closes the valve. It preserves the exact tank/source state,
uses zero engine mass flow/thrust/exhaust spin, and advances native mechanics
under gravity and other supplied forces/torques. A failed gas/nozzle solve never
implicitly selects coast. Both powered and coast operations retain the native
fixed-COM/shared-principal-axis mass model and separately declared diagonal
rotational exhaust flux; no slosh, general outlet tensor, internal-flow momentum
or complete world exhaust-energy audit is inferred.

Native actual-source verification passed seven rigid-coupling tests: rotating
and off-axis thrust, scalar/world impulse distinction, pressure torque once,
second-order coupled position/velocity/attitude refinement, cross-owner mass,
late cancellation, source/mechanical budgets and exact-gas closed-valve coast.
The existing point-gas coupling and frame-correct environment suite also passed
all nine tests, including wind-relative drag, rotational application-point
power and central-gravity potential. Full workspace validation remains separate.


## Finite winding temperature in the biased magnetostrictive circuit

`material_circuit::magnetostrictive::thermal` composes the original
`MagnetostrictiveCircuit`, source-pinned conductor resolver and
`circuit_enthalpy::cell_driver_candidate`. A fixed three-branch series circuit
contains one stationary winding resistor, the reciprocal rod/coil storage and
a source specified by TOTAL voltage. The bias current is a linearization
reference, not an implicitly maintained current source. Every thermal endpoint
trial restarts the original native descriptor from identical accepted linkage,
extension and momentum, with incremental voltage
`V_total - R(T_endpoint) * I_bias`.

The winding receives the actual nonnegative
`R(T_endpoint) * (I_bias + delta_i_mid)^2 * dt`, including the DC and cross
terms. The total magnetic/mechanical storage CHANGE includes
`I_bias * delta_lambda`; this is not an invented absolute biased energy datum.
Independent work and power gates check total terminal work, reciprocal mechanical
work, stored-energy change and this heat before the thermal owner deposits it
once. Conductivity/resistivity source identity, accepted winding state, clocks,
source envelopes, actual solver budgets and late cancellation remain active.
Only a jointly accepted candidate advances either history.

The winding cell must have its original wholly solid, constant-density caloric
chart and stationary geometry. Other bound cells/reservoirs are passive thermal
recipients. Magnetic/core coefficients and the rod's prepared biased reference
remain at their independently held preparation temperature; no magnetic
temperature law, core heating, hysteresis, eddy loss, thermal strain or expansion
is inferred. Endpoint electrical coefficients plus backward-Euler enthalpy give
first-order coupled trajectories, despite the native circuit's midpoint substep.

The focused `thermal::tests` qualify actual DC heat with zero incremental loss,
independent reciprocal midpoint/current/total-work equations and regeneration, both conductivity
and resistivity thermal feedback at the same accepted endpoint, refinement
against independent continuous coupled equations, simultaneous bias/winding
orientation reversal, and source/budget/cancellation rollback. Synthetic cards
qualify this numerical/physical composition, not empirical device accuracy.


## Atmospheric source-gas rigid flight with explicit propulsion mode

`rocket_blowdown::environment::spherical::SphericalStandardAtmosphere` caches
the original lower USSA1976 profile and an independently declared central
gravity source. It samples geometric altitude at the actual inertial point
and assigns local bulk air velocity `Omega cross (position-center) + wind`.
Radius, center, gravitational parameter, rotation and background wind are
explicit. The profile's effective hydrostatic reference gravity is not silently
used as the vehicle gravity law. No fictitious Coriolis/centrifugal forces are
added to inertial dynamics. Air outside the source's 0..80 km geometric domain
refuses; central gravity remains independently usable outside that air domain.
This spherical rigidly rotating reference is not weather, terrain, oblateness
or a self-consistent rotating hydrostatic atmosphere.

`rocket_blowdown::flight::advance_atmospheric_rocket` composes that environment
with the original source-gas/rigid vehicle. It resolves separate actual body
points for aerodynamic loading and nozzle ambient pressure, native angular
velocity at the aerodynamic point, and gravity/potential at the COM. The caller's
pure/repeatable drag-coefficient law receives actual Mach, Reynolds and body-
relative flow; the same force supplies its matching body torque. No empirical
coefficient, lift derivative or source extrapolation is invented.

One immutable half-step predicts the midpoint environment. The full accepted
interval then restarts from the original gas and vehicle under those midpoint
loads; trial propellant is never subtracted twice. The original gas midpoint
quadrature and bounded native attitude substeps retain their own controls and
reports. Final atmosphere and coefficient-domain admission precede publication.
The candidate retains original accepted gas/mechanical reports, midpoint loads
and independent endpoint diagnostics. Powered and explicit closed-valve modes
are distinct, and source/nozzle failure never silently turns into coasting.

This smooth composition is second order, with physical outer timestep and
native substep controls; it does not locate impact, cutoff, staging or source-
domain events. The caller splits intervals at known mode changes and refines
changing flow regimes. The underlying fixed common COM/principal axes, declared
diagonal exhaust angular flux, ideal frozen nozzle and source-domain limits
remain in force. General internal-flow/transverse outlet linear momentum,
slosh, aerodynamic heating and a full exhaust/world-energy audit are not inferred.

The combined native environment/rigid/flight suite passed all 19 tests, including
rotating-air frame covariance, point-velocity drag/torque, mass and scalar-impulse
closure, source/cancellation atomicity and measured second-order coast
position/velocity with changing drag and gravity. Both native atmospheric
example tests passed. The complete 2 s manufactured example ran at 0.01 s and
0.005 s outer timesteps, with 1 s powered then 1 s closed valve. Coarse-run total
mass residual was -1.084e-18 kg, and refined final altitude differed by about
4.55e-5 m. Those fixture results are not calibrated vehicle performance or a
full-workspace validation claim.

`cargo run -p fs-couple --example atmospheric_rocket --
crates/fs-material/examples/material_gas_spring/mixture-a.tsv 1000 1000000`
exercises the complete path and writes time-labelled CSV diagnostics. Its
README specifies manufactured Cd, spherical gravity/rotation choices, uniform
tank-volume/inertia geometry, prescribed valve schedule, source validity and
model boundaries. The final gas inventory is unchanged throughout closed-valve
coast while native motion, gravity and drag continue.

## Native retained multi-contact caloric feedback

`tribothermal::partial_slip::rigid::network::step_rigid_caloric_network` resolves all original shear/spin impulses, finite-cell energy increments and one massless common temperature per original interface in one Newton system. Contacts may share native bodies and thermal cells. Canonical original link order fixes evaluation and reduction order. Every trial applies all impulses before evaluating any midpoint contact motion; the transport callback sees the whole endpoint once per trial. Generic fs-solver Newton/Krylov and finite-difference actions own the numerical algorithm. No replacement inertia or thermal store is introduced.

The existing `CaloricPartialSlip::mechanical_endpoint` retains each contact's original regional strain, affine two-temperature stiffness/free energy, pinned friction source and applicability, and accepted history. Physical microslip loss enters each original common node; each surface retains its own reversible `-T_i Delta s_i`. Conductances derived from original geometry and endpoint conductivity determine generated-heat partition. Signed interbody conduction, contact internal energy and numerical return-map loss remain separate. Shared finite-cell storage and outside supplies are counted once.

Acceptance checks every original constitutive impulse/couple, each body's point-wrench momentum response, total linear and orbital-plus-spin angular momentum, original cell/interface balances, native kinetic work, and `Delta(K + sum U_contact + H_caloric) + sum D_numerical - Q_external`. A candidate returns all original bodies, contacts and thermal state together. Source, topology, applicability, budget, cancellation or balance refusal leaves accepted inputs unchanged.

Four focused G0/G1/G3/G4/G5 groups in `tests/rigid_caloric_network.rs` cover two contacts on three actual bodies and three finite caloric cells against independent 2-by-2 mobility and tridiagonal entropy/heat references (2e-7 K and 2e-9 N s); simultaneous gross sliding against independent regional plastic loss and shared-cell deposition (1e-9 J); 8/16/32 retained-step continuation converging at first order to the coupled elastic oscillator (relative-velocity errors 0.0292071/0.0157971/0.00821031 m/s); and permutation, late source/cancel/budget refusal and exact retry. Four existing native single-pair caloric groups also pass.

Scope: frozen poses, supplied positive normal elliptical patches, coincident original orthogonal thermal faces, retained two-region law and declared affine temperature-stiffness domain. Nonlinear convergence does not establish uniqueness or temporal accuracy. Endpoint constitutive integration is first order. No evolving normal load, moving footprint, arbitrary geometry, nonlinear stiffness thermodynamics, calibration, numerical-loss thermalization, full Cargo/DSR workspace validation or blanket bead closure is claimed.

## Native finite-thickness shear-layer mechanics and heat

`friction_heat::third_body::rigid` adds an explicit first-order, quasi-static shear-layer director interpolation between the actual separated native faces. With `d = x_B - x_A`, shear rate is the tangential part of `v_B - v_A - (omega_A + omega_B) cross d / 2`. Virtual work derives opposing face forces `+F/-F` and equal face couples `d cross F / 2`. These are declared interpolation terms, not unexplained point-contact corrections. They close original orbital-plus-spin angular momentum and vanish in work under common rigid motion.

Original masses, principal inertias, orientations and material arms supply the coupled two-coordinate mobility, including off-diagonal terms. Geometry/inertial mappings are cached per solve. The retained `ThirdBodyLayerContact` law owns midpoint elastic storage and physical Bingham plastic work. Independent audits check each actual face wrench and `Delta(K + U_layer) + D_plastic`.

`step_heated_rigid_layer_candidate` jointly solves two shear forces and all ORIGINAL finite-cell energy increments, without nested mechanical solves, replacement stores or authored heat fractions. Physical plastic work enters the original layer cell once. Resolved conductivity/resistance determines exchange with both bodies; the real layer mass/caloric law retains heat, and the same endpoint temperature resolves yield and viscosity. Only unchanged elastic potential/domain may continue; changed stiffness refuses without a thermoelastic work port. Actual native faces must match the original thermal face centres, area and thickness. Accepted inputs remain immutable, and one candidate returns bodies, history and thermal network together.

Five focused G0/G1/G3/G4 groups in `tests/rigid_layer_contact.rs` cover independent anisotropic elastic/plastic 2-by-2 references, actual face force and moment, common rigid-motion null response, frame rotation, a 16/32/64 midpoint-elastic refinement ladder, finite thermal feedback, retained continuation and entropy, and late source/cancel/budget/geometry/elastic-potential refusal with exact retry. Elastic displacement errors are 1.02619e-8/2.58197e-9/6.46526e-10 m, approximately fourfold reduction per refinement. The nonlinear thermal reference gives layer temperature rise 0.03176733394 K, force 0.29474358175 N and plastic heat 0.01086889051 J.

Scope: frozen geometry and explicitly massless mechanical, shear-only layer interpolation. Thermal layer mass remains physical; neglect of mechanical inertia is a declared quasi-static approximation. No normal response, bending/twist stiffness, finite-layer mechanical inertia, finite-rotation state transport, wear/opening/advection, continuum thermal convergence, measured calibration, full Cargo/DSR workspace validation or blanket bead closure is claimed.

## Moving material partial-slip history and finite thermal feedback

`friction_heat::moving::partial_slip::MovingPartialSlip` attaches original two-region shear strain and kinetic-branch history to surface A material coordinates. A declared uniform-pressure brush moves over the original planar strips. Surviving parcels retain complete states without interpolation or averaging; entering material is explicitly unloaded. Departing elastic storage is exported as mechanical boundary flux, never silently erased or relabelled as frictional heat.

`step_candidate` performs transport followed by the original regional return maps, a first-order Lie split. Regional extensive force, work and storage scale with physical area; original strains remain intensive. The reference ellipse supplies total area/load and the constitutive law, not an inferred elliptical pressure distribution. Body offsets determine actual longitudinal slip independently of footprint motion and material residence.

One original `MovingFrictionHeatStep` owns the complete finite-cell thermal transaction. Original cell-edge trajectories and each material region's swept overlap weights are prepared once outside Newton iterations. Only physical microslip work enters the common interface nodes; each cell retains its original thermal mass, material law and history. The callback sees the proposed endpoint, original material-region identities, swept cell pairs and local common temperatures. It can resolve regional friction coefficients from actual material temperatures. Source/model identities and the original elastic potential/domain must remain unchanged; unsupported stiffness changes refuse.

The report retains total shear force and its distributed moment about the strip origin, original retained storage, exported elastic energy, physical heat and numerical return-map loss. Before atomic publication, original thermal balances and `Delta(U_contact + H_thermal) + D_numerical + U_out - W_contact - Q_other` must close. Source, geometry, budget, cancellation or energy refusal publishes neither transported history nor thermal state.

Six focused G0/G1/G3/G4/G5 groups in `tests/friction_heat_moving_partial_slip.rs` verify exact preload transport/outflow, swept-cell heat deposition, common translation including distributed moment, independent nonlinear local temperature feedback, 16/32/64 brush-integral refinement, and late refusal with exact retry. A 0.011 J preload retains 0.00825 J and exports 0.00275 J with no fictitious heat. A sliding case deposits 0.00625 J of physical work once. Finest force/storage/heat errors are 0.00761719 N, 5.37618e-5 J and 3.90625e-6 J. Force, storage and history transport are first order. The integrated heat's second-order envelope is limited to the constant-coefficient fixture.

Scope: prescribed planar translation, fixed uniform pressure, fixed reference stiffness, surface A brush memory, longitudinal imposed slip and original finite thermal cells. No torsional/spatial pressure-field law, native body dynamics, evolving normal load, recirculated outflow memory, curved surfaces, large rotations, continuum thermal convergence, calibration, Cargo/DSR workspace validation or blanket bead closure is claimed.

## Finite-temperature three-dimensional magnetic induction

`circuit_enthalpy::spatial_induction::SpatialMagneticThermalNetwork::step_candidate`
composes the original retained Whitney A–phi induction state with the original
finite-enthalpy Newton owner. Every temperature trial starts from the same
accepted field and thermal histories, resolving the same pinned conductivity
and permeability sources at the proposed local temperatures. Immutable geometry,
incidence and Whitney basis allocations are shared across material refreshes.

Complete original tetrahedron-to-thermal correspondence routes actual native
Ohmic heat to each whole-volume solid receiver exactly once. Original
conductivity/resistivity basis and claim identities, constant density, common
clocks, supported temperature/B envelopes and positive combined caloric plus
magnetic heat capacity are required. The shared original magnetic caloric owner
now accepts two- or three-component fields and physical cell volumes; planar
entry points retain their existing behavior.

Signed reversible heat Q_rev=-P-delta(T*S_field) is separate from Joule heat.
The independent physical balance is
delta(U_field)+delta(H_caloric)+W_mechanical-W_impressed-Q_external=0,
using magnetic internal energy rather than duplicating free-energy storage.
Prescribed motion keeps explicit external Lorentz work. The accepted thermal and
field histories are returned together; source, time, volume, work-budget or
cancellation refusal returns neither. Already deposited heat must not be added
again by a caller.

Five new physical tests compare same-step field/temperature solutions with
independent scalar elimination across both electrical source bases and both
permeability-slope signs, six distinct unordered local receivers, independently
derived motion work, and first-order time refinement against continuous
equations. Atomic refusal covers capacity, source/domain, physical volume,
clock/event, actual trial budget and late cancellation. Six existing planar
eddy/enthalpy tests also pass; their thermal fixtures explicitly declare the
absolute-temperature source type required by the original entropy owner.
These are manufactured equation/discretization checks, not measured devices.

The thermal path admits scalar conductivity and induced-linear globally affine
source-backed mu(T), fixed geometry and prescribed velocity. Recoil, saturation,
hysteresis, thermal expansion, moving interfaces and mechanical integration
require compatible additional physical owners; they are not inferred here.

## Nonlinear plate damping and physical port normalization

`VkBody::from_modal_storage_with_damping` binds a source-reduced quartic
storage, the complete mass-normalized viscous matrix, physical force
projections, and signed surface integrals to the original `fs-nlmodal`/`fs-phs`
stepper. The same basis must supply all four inputs. Off-diagonal loss and
coupled lossless directions survive; no extra damping ratio or loss floor is
added. The original storage owner admits frequencies and quartic channels,
and the original pHS owner admits the resistance matrix. This does not infer a
constitutive model, a material source, or radiation back-reaction.

Analytic simply-supported plate force and face ports now use the same
`2/sqrt(rho*h*Lx*Ly)` sine-product normalization as nonlinear storage. The
independent two-dimensional modal mass integral is one, and the conjugate
point-force-to-volume-flow gain scales inversely with surface density.

G1/G3 tests retain the exact common/differential motion of a rank-one damper,
held-force subdivision in the linear limit, signed compact pressure/flow,
and physical-input refusals. The compact pressure observation uses step-mean
acceleration; it is not a retarded or finite-wavelength acoustic field. The
existing nonlinear time-resolution gate and source/reduction limits remain.

### Shared material storage on spatial electrothermal circuits

`ConductionCircuit::new_with_material_storage` applies the original material
capacitor/inductor reductions before ordinary circuit admission. The spatial
conductor, other material resistors, and field storage therefore share the
existing circuit/finite-heat trajectory without independent C/L overrides.
Complete dielectric/core source records and geometry survive continuation.

Storage contexts remain fixed during electrothermal intervals.
`rebind_storage_candidate` changes source/geometry at unchanged charge and
flux linkage and reports required external work from the SAME material energy
law, delta(q²/(2C) + lambda²/(2L)), shared with the lumped material circuit.
The spatial field owner is retained unchanged. Incompatible ideal-source
constraints, incomplete sources and cancellation refuse transactionally.
Parameter work is not Joule heat. Automatic temperature-dependent dielectric
or magnetic storage requires compatible reversible field-entropy ownership.
The `conduction_enthalpy::circuit::storage` regressions cover independent
RC/RL source swaps, actual spatial fields/heat, work reversal and refusal.

### Material frame viscosity through the common modal owner

`FrameVibrationBasis` retains the native `FrameMaterialDamping` and its original
section sources. `coupled_component` projects each native positive resistance
port through the unchanged mass-normalized basis. The existing
`CoupledModalSystem` receives every complete column, so its ordinary and
prepared propagators preserve the full Galerkin damping operator Phi^T C Phi,
including nonproportional cross terms. There is no diagonal damping fit,
second integrator or invented pressure transfer.

`FrameVibrationComponent` keeps source basis, oscillator state and all material
connections together. `into_system` prepares the original coupled owner;
`into_parts` remaps every port for a larger contact/radiation/modal assembly.
The returned basis retains physical force/observation maps and provenance.
Original nonzero modal states survive preparation. Material loss enters
`connection_dissipation_j`; separately authored modal loss enters
`component_dissipation_j`. Neither is implicitly deposited into thermal cells.

The scalar `time_model` and its independent-mode runtime refuse source-viscous
frames instead of dropping their coupling. Complete native-spectrum and
material-band admission remains in force after explicit modal truncation.
Projection and native setup work are bounded. Truncation does not reconstruct
omitted motion/loss. `tests/frame_vibration_material.rs` specifies continuous
physical references, heterogeneous off-diagonal damping, larger-assembly
work/state continuity and rejection behavior. The existing `frame_vibration`
example now accepts manufactured density, stiffness and viscosity data and
uses this complete shared path for its SI response and loss output.

## Directional material damping in existing modal consumers

The canonical `ModalPlateViscosity::factor` retains its original full-rank
Cholesky path and falls back to the existing PSD channel owner for directional
loss and undamped modal subspaces. Plate and both shell consumers therefore
retain the same complete source-projected matrix without another realization.
The zero-loss shortcut requires an exactly zero matrix: nonzero zero-diagonal,
wholly negative and asymmetric zero-diagonal matrices refuse, and negative
diagonal checks scale with the actual tensor. Native regressions verify a
rank-one damper's lossless common motion, complete differential work, retained
full-rank arithmetic and active off-diagonal refusals across coefficient scales.

## Nonlinear plate material memory

`VkBody::with_relaxation_branches` composes the existing Maxwell owner with the
body's original nonlinear storage and full viscous damping matrix. Branches
project only displacement entries in the same mass-normalized basis; one
energy-normalized material history per branch follows the mechanical state.
Complete cross-mode hereditary stiffness and loss are retained, with no
replacement by scalar damping ratios. The equilibrium storage excludes the
added relaxing stiffness. Attachment requires zero original state, no prior
memory, and an explicit memory-state ceiling; shared pHS/render capacity and
full-rate time-resolution checks remain in force. Empty branches are an exact
no-op. Read-only state and total mechanical-plus-memory energy are exposed.

`nonlinear_plate_relaxation` checks the linear exact limit against independent
continuous equations and held-input subdivision, a physically unstrained
lossless combination, and nonlinear convergence over 200/400/800 steps against
independent RK4. It also checks complete energy decay after forcing stops,
compact step-mean pressure/volume-flow consistency, and admission/rollback.
Native tests pass without relaxing the original time-resolution bound. Branch
coefficients and applicability must come from the caller's source reduction;
this API does not automatically derive a hereditary material tensor or add
finite-wavelength acoustic back-reaction.

## Spatial electrothermal field storage and finite frame heat feedback

### One spatial electrical/thermal state

`ConductionThermalNetwork::step_field_circuit_candidate` composes the original
spatial conduction circuit with the same field-storage thermodynamics already
used by the lumped circuit. The original charge/flux coordinates, source
conductivity, surrounding material resistance, capacitance/inductance and finite
caloric cells are evaluated inside one thermal residual. Geometry, caloric
ownership, source claims and state-dependent heat-capacity admission remain
mandatory. The existing fixed-state field API and default report types remain
valid.

The common field law retains free energy `psi=z²/(2K(T))`, entropy
`S=-partial_T psi`, and physical internal energy `U=psi+T*S`. Its reversible
transfer `-delta_psi_parameter-delta(T*S)` is internal heat, distinct from the
native local Joule integrals and external thermal sources. The enclosing
balance uses physical field internal energy plus caloric storage. Source/geometry,
heat-capacity, budget or cancellation refusal leaves both original histories
unchanged. The spatial RC/RL regressions use independent scalar electrical and
caloric equations, including positive, negative and zero coefficient slopes,
thermal resistance feedback, source refusal and continuation.

### Original frame losses heat their own material

`FrameVibrationComponent::into_thermal_system` binds each original member layer
to exactly one whole finite material cell through `FrameLayerCell`. The material
card, source context and caloric object remain original; receiving volume is
`A_i*L` and mass is `rho_i*A_i*L`. Every participating cell must stay in its
fully solid, constant-density caloric domain. Original scalar density and
Young's modulus must explicitly declare constant-within-validity behavior.
This rung retains the reference geometry, mass and elastic basis.

`FrameThermalSystem::step` stages the original ordinary or prepared modal step,
then deposits its independently integrated axial and both curvature-port work
integrals. Axial heat partitions by the original `eta_i*A_i`; bending heat
partitions by `eta_i*I_i`. Neither an energy defect nor an authored heat fraction
defines physical dissipation. The same `PreparedMaterialNetwork` and finite
enthalpy/Newton residual advance source conductivity, conduction, reservoirs and
caloric storage. Authored additional modal losses remain an explicit external
sink. Reports expose material heat, source evaluations, combined energy balance,
actual caloric entropy, conductive production, material `D_i/T_i`, and separate
temporal/algebraic entropy defects.

The complete original section bundles are refreshed with
`ConcentricSectionTemperatureLaw`. Original axial and bending viscosity
properties may vary with temperature. Each interval uses its accepted starting
viscosities and layer heat shares; its resulting finite-cell temperatures select
the NEXT interval's viscosities and shares. This is a first-order
thermal/coefficient split, not a simultaneous within-interval nonlinear
mechanical integration. Refine the shared time step for feedback accuracy.
Every endpoint viscosity must support DC through the actual sampled Nyquist
band; the complete fixed elastic spectrum was admitted when constructing the
frame. Changing an ideal dashpot has no stored energy or parameter-work impulse.
Temperature-dependent elasticity, changing density/geometry, phase change and
thermoelastic work are outside this fixed-reference rung and refuse.

The live section accessor `viscous_sections()` and report's
`interval_viscosity`/`next_viscosity` distinguish current source state from
the fixed `basis()` construction records. Independently altered thermal
temperatures cannot be paired silently with old accepted viscosity.
Source evaluation, independent heat/energy admission, endpoint reaction checks
and operator preparation all finish before either accepted trajectory advances.
A refused step preserves the modal clock, states, port ledger, finite cells,
current coefficients and heat shares, so the same input can be retried.

### Resistance refresh on the complete coupled propagator

`CoupledModalSystem::set_connection_damping` replaces a sorted, unique subset
of the original nonnegative physical resistances. It rebuilds the SAME
ordinary/prepared operator with all original attachments, springs, source port
identities, Maxwell relaxation branches and acoustic radiation bindings. Modal
coordinates, stored material/acoustic histories, sample clock and previous
frame/transfer ledger remain unchanged. Constant-loss anchors retain their
positive-coefficient requirement. Existing setup budgets and instantaneous
reaction bounds apply; failed preparation or cancellation is atomic. Empty or
unchanged updates allocate and rebuild nothing. The next ordinary step retains
its full-interval force and energy checks.

The finite-frame owner prepares its resistance change against the staged
mechanical endpoint inside the native admission callback. The subsequent
publication consists only of moves/assignments. Regressions cover analytical
axial/bending motion and original layer integrals; independent two-cell thermal
conduction; heat-capacity data swaps; temperature-dependent viscosity driving
both subsequent motion and heat partition; ordinary/prepared continuation;
independent Maxwell/acoustic physical ODEs; coefficient/refusal/cancellation
atomicity and deterministic retry. These are native test targets; source review
or formatting checks alone are not evidence of numerical execution.

## Finite cohesive patch load line

`cohesive_contact::patch::CohesivePatchLoadLine` solves all opening and two shear coordinates of up to 64 original bonded quadrature sites together. It consumes the full caller mechanical cross-compliance, heterogeneous original `CohesivePoint` histories and reference areas. The original analytic mean traction supplies physical work. Finite lever arms transmit moments; local failure redistributes force through the retained off-diagonal response. No average damage or independent-site solution replaces the patch.

Prepared scratch and LU storage allocate only at construction. Each Jacobian perturbation evaluates only its local point law, giving linear constitutive-evaluation work per Jacobian while retaining the bounded dense solve. Explicit total point-law evaluation, iteration, backtracking, per-site force/jump and compatibility budgets apply. Refusal/cancellation mutates scratch only; callers must independently audit actual mechanical motion before publishing all histories. Softening continuation does not certify root uniqueness, spatial resolution, finite-rotation objectivity or timestep accuracy. Tests cover heterogeneous full cross-compliance, force redistribution after local failure, exact finite-body force/torque/work balance with progressive edge peeling, the one-site limit, and atomic refusal.

## Native finite cohesive patches

`contact::cohesive::patch::CohesivePatchModalSystem` attaches heterogeneous `CohesiveModalSite` quadrature to the original `CoupledModalSystem`. Each site retains a borrowed material law or the original ordered source, reference area, and three complete modal work columns. Setup computes every cross-site/cross-axis compliance term through the original held-force response under an explicit projection-work budget. Translational, rotational and elastic participation therefore retain their reciprocal forces and moments. Construction explicitly bonds the current sample-zero configuration; it does not infer bonding from proximity.

The prepared patch solve stages original native mechanics with the common held force field. It independently reevaluates every material point on actual body motion and gates each jump and force component, then the complete native-plus-bond energy window. Material calls reserve one final audit per site. Only after all checks and the final cancellation poll do all body states, irreversible histories and per-site observations publish together. Original fracture, surviving-bond viscosity, native losses, external work and prescribed jump work remain separate. Contact scratch allocates at construction; the original propagator retains its existing allocation rules.

The native integration regressions cover quadrature-area partition, finite lever-arm force/torque and edge peeling, and multi-site refusal/replay. The same finite-body fixture passes the independently executable original-source core tests; the full native executable still requires the repository runtime build. Fixed frames, linearized geometry, declared spatial quadrature and straight local jump paths remain explicit; no finite-rotation contact search, remeshing/history transfer, spontaneous rebonding, thermal coefficient refresh or automatic fracture heat partition is inferred.

## Adaptive constant-pressure reacting gas with prescribed heat

`GasReaction::isobaric_adaptive_candidate` advances the original closed gas
against its initial constant pressure. Signed external heat is constant power
over the requested interval. Every accepted substep uses the existing coupled
implicit mechanism, with composition-dependent source enthalpy, temperature
and actual EOS volume solved together. Chemical formation energy stays in the
original total enthalpy and is counted once.

The bounded owner compares complete N and 2N backward-Euler trajectories using
every species amount, every channel extent, temperature and actual volume.
Only the finer physical trajectory can be accepted; there is no extrapolated
composition or concentration floor. The reported discrepancy is an asymptotic
first-order estimate, not a rigorous error guarantee for arbitrary stiffness.

Whole-window gates retain actual mass, pressure and reconstructed species
inventory; compensated accepted source increments; H-H0=Q and U-U0+p Delta V=Q;
and accepted heat entropy sum(Q_k/T_k). Absolute conservation budgets are not
multiplied by substep count. The history retains every accepted original
endpoint report, while discarded paths contribute only to work accounting.
All attempted endpoints and original residual/Jacobian evaluations, including
refused roots, spend the caller's global budgets. Cancellation and refusal
publish no partial state. Original public single-step APIs are unchanged.

Independent tests cover reversible three-species exact kinetics, stiff and
transient trajectories, exothermic heating/adiabatic/cooling histories against
a scalar reference ODE, mole-changing volume-controlled refinement, exact
failed-attempt work accounting, source limits and late cancellation.

## Complete multimodal flanged loads in the existing acoustic lines

`mm_line::MmLoad::FlangedModal(&exterior)` selects the complete physical
flanged-aperture matrix for both `MmLineBank` and `MatrixFirLine`. The same
signed mutual/higher-mode load is used at every positive-frequency bin and
at the existing near-DC probe. Exterior gas and full Nyquist-band coverage
are admitted before the FIR sweep, with common interior/exterior static
pressure. This strict selector does not substitute a scalar termination,
clip reflectance or default the DC response. Existing causal, tail, peak,
passivity and modal/axial refinement gates remain unchanged. The original
scalar selectors and shared `MouthLoad` API retain their existing behavior.

The new `mm_line/flanged_tests.rs` reconstructs physical mouth pressure and
volume flow from the complete reflected modal waves and checks their full
matrix boundary and work, plus exterior-pressure and Nyquist admission.

## Consistent inertia in exact tensor plane-stress dynamics

`TensorMaterialBody::from_plane_stress_mass_matrix` uses the existing validated
full-mass `Inertia` owner. Off-diagonal inertia participates in canonical
momentum `p=Mv`, kinetic storage, force-driven motion, complete Maxwell memory,
and independent work/heat accounting. Complete state capacity is checked before
mass factorization; no inertia coupling or constitutive coordinate is removed.
The original lumped constructor and exact plane-stress condensation remain.

`TensorMembrane::new_with_consistent_mass` assembles the exact P1 triangle mass
`rho*h*A/12 * [[2,1,1],[1,2,1],[1,1,2]]` for each Cartesian component and
restricts both matrix axes by the original support map. Full-mesh row sums
remain physical load weights; they do not replace the constrained kinetic
matrix. Four new native tests cover independent physical-velocity RK4 with
anisotropic memory and full mass, initial momentum/kinetic cross terms,
independent source work and heat, zero normal traction, diagonal compatibility,
invalid inertia, exact velocity quadrature and a three-rung continuum-frequency
ladder. Existing membrane and distributed-loading tests also pass. This remains
a linear-kinematics membrane with the original exact tensor-memory limitations;
it does not add bending, prestress or finite deformation.


## Prepared complete plate vibration on a fixed clock

`VkBody::into_prepared(dt_s)` consumes an accepted body and preserves every
original displacement, momentum and Maxwell-history component. The opt-in
`PreparedVkBody` keeps the fixed positive physical clock and the complete
original force, surface, damping and storage owners. It does not reduce the
state or replace a coupled operator with a diagonal approximation.

Quadratic bodies use the owned `PreparedExactZoh` cache, including independent
integrated heat and source work. Nonlinear bodies reuse the canonical Gonzalez
equation and `StepWorkspace`, with the existing midpoint pole-rate admission
checked for the fixed clock. Both paths reuse scratch during accepted warm
steps. This is an allocation-reuse guarantee for these plate storage owners,
not a hard real-time timing certificate.

`drive_and_radiate[_controlled]` returns a `PreparedVkStep` containing the
original mechanical/material work record, endpoint volume velocity, and compact
baffled-monopole pressure from interval-mean acceleration. Density and listener
distance remain observation parameters. Preparation retains explicitly attached
acoustic loading without adding another radiation reaction or a finite-wavelength
receiver approximation. Invalid inputs, solver or
observation refusal, and cancellation publish no physical-state or
constitutive-history change.

The three prepared-plate regressions compare 150 forced steps against the
original quadratic and nonlinear bodies with full damping and Maxwell memory,
independently check physical force work and the complete energy balance, and
measure zero allocations over 128 accepted warm steps in each path. They also
cover preserving a moving initial body, the nonlinear phase bound, invalid
clock/observer admission, early and late cancellation, and retry after refusal.
These checks exercise the production numerical owners through the focused
native target; they do not assert full workspace qualification.

## Nonlinear electrical material response in the thermal piston

`GasPiston::circuit_conductivity_thermal_candidate` composes the existing
`MaterialCircuit::conductivity_with_materials_candidate` with the existing
powered body/gas/piston feedback. The body enthalpy midpoint updates temperature
and any declared phase coordinate; the native electrical midpoint updates each
explicitly bound electric-field/current-density coordinate. The original full
material bundle, selected source claims, geometry and conductivity/resistivity
basis remain bound throughout both iterations.

Acceptance requires the original electrical constitutive and power residuals,
thermal/source residuals, and complete circuit/body/gas/kinetic energy balance.
Only native circuit Joule heat enters the body. Evaluated material power defects
remain separate residuals; capacitor/inductor depletion is not counted again as
heat. Every trial starts from the same accepted charge/flux and thermal history.
The global work cap counts discarded native electrical iterates as well as the
last iterate. Refusal, nonconvergence, budget exhaustion and cancellation return
no partially accepted physical state.

The existing circuit thermal tests exercise a manufactured two-coordinate
resistivity surface, independent RC midpoint identities, actual body/gas/motion
sensitivity, unbound/out-of-domain operating axes, electrical nonconvergence,
global work exhaustion and final cancellation. This retains the fixed C/L and
geometry, one-whole-resistor, uniform thermal body and prescribed nonthermal
context restrictions. It adds no contact refresh, carrier transport, hysteresis,
spatial temperature field or uniqueness/time-accuracy certificate.


### Frozen-composition caloric limit of reactive owners

When the complete original thermodynamic models and active compositions are
identical before and after a proposed reaction boundary step, the shared
reaction admission uses the original source's specific internal-energy and
entropy differences. The isobaric report likewise uses the source enthalpy
difference. Actual initial/final pressures remain in the entropy integral.
This preserves the material parcel's stable heating/cooling and expansion
limits when a large common formation reference would erase an absolute-total
subtraction. The strict model guard includes all retained species sources,
conventions, molar masses, composition and caloric-continuity policy.

Changed compositions continue to use their original formation-inclusive
energy and entropy accounts. Tests retain identical physics under positive
and negative common reference shifts, replay the original isobaric and
fixed-volume heat owners, resolve free-expansion entropy, and refuse
compression with no supplied heat or work.

### Adaptive reciprocal finite-body heat and chemistry

`GasReaction::isobaric_thermal_body_adaptive_candidate` advances the reacting
gas and original finite solid/liquid body together through the existing
simultaneous chemical/thermal endpoint. The corresponding
`isobaric_thermal_layers_adaptive_candidate` retains original material-card
conductivity selection and source receipts at every actual reacted gas/body
endpoint. Both continuing states are returned together, with contact heat
already applied once to each side with opposite signs.

Temporal comparison extends the shared gas species/extents/temperature/volume
controls with body temperature and body specific enthalpy. The body's relative
enthalpy scale is its change from the initial state, so an arbitrary datum does
not loosen accuracy and latent progress can force refinement at constant T.
Only the finer physical history is retained, with every joint source report
stored once. Prescribed gas/body heats are uniform powers; subinterval heat
brackets retain the declared exchange-power bounds.

Whole-window gates independently check gas Delta H, gas Delta U plus actual
p Delta V, body storage against actual final m Delta h, reciprocal/evaluated
contact heat, and combined first/second-law accounts. Both prescribed heat
entropy contributions use the actual local endpoint temperatures. Source
tolerances divide over accepted substeps; absolute conservation and entropy
budgets never grow with substep count. Failed thermal roots and subsequent
source refusals spend their actual nested chemical evaluations before retry.

Tests compare hot/cold exchange with its exact solution and nonlinear-Cp
exothermic feedback with an independent coupled ODE, resolve sensible-to-latent
entry and body-enthalpy-only refinement, retain original material-layer sources,
and independently count refused heat-root work and final cancellation. This
remains homogeneous ideal gas plus a uniform fixed-geometry solid/liquid body
with negligible body pressure-volume work; water-boiling expansion is not
folded into that body assumption.

## Complete tensor relaxation in material-modal valves

`PlaneStressPoint::relaxation_image` changes coordinates in the existing
complete Maxwell plane-stress energy. For its retained block energy
`H=[q,a]^T[A B;B^T D][q,a]/2`, the shared symmetric eigensolver diagonalizes
`sqrt(R) D sqrt(R)` with `R=diag(1/tau_source)`. The resulting independent
passive channels reproduce the complete energy and original dissipation;
all six history coordinates per source branch survive, including zero or
dependent source factors. The inverse map reconstructs the original strain
and normal-traction constraint. No branchwise plane-stress approximation,
fitted relaxation times or mode truncation is introduced. Orthogonality and
energy reconstruction have explicit numerical gates, with cancellation around
the bounded dense transform.

`ModalAperture::with_tensor_material_relaxation` applies that image at the
original three DKT area points per element, with original region, ordered
material axes, angle, thickness and retained modal coordinates. The equilibrium
section must reproduce the existing chart to `1e-9` relative diagonal scale;
its stiffness is counted once. Exact channels enter the existing mechanical
propagator, pressure/contact solve, swept volume and independent dissipation
audit. The full 3D outer-face strain is checked before accepted state publication.
Initial histories are explicitly relaxed or unrelaxed. Original source bands,
instantaneous-frequency/relaxation clock limits and complete dense-state
capacity remain mandatory; capacity refusal never drops a material history.

Both native material-modal formats share `loss tensor-maxwell` and the generic
`fs-material` source binding. Their data identify the small-strain claim and
each complete Mandel branch, with material-point, memory and projection-work
budgets. Source cards retain all coefficients and evidence; metadata identifies
full tensor memory separately from proportional/regional scalar spectra. The
canonical paired DG example is reachable from Instrument Lab with its original
material pack. Scope remains fixed-state, homogeneous-through-thickness linear
plate bending; nonzero 3D Kelvin normal-rate terms refuse. Measured wet-cane,
nonlinear solid, resolved moving-flow and complete instrument accuracy are not
established. The new complex-Schur, energy, source substitution and rollback
regressions are accompanied by independent formulation checks; native Rust
execution of these additions remains unverified under the shared resource limit.

## Live material elasticity, expansion and finite heat in the original frame

The original coupled mechanical owner now accepts one atomic batch of
`ConnectionCoefficients`: stiffness, rest extension and damping. It rebuilds
its existing ordinary or prepared propagator, complete frequency bound and
endpoint reactions without resetting q/v, Maxwell memory, radiation history,
ports or sample time. At unchanged coordinates it reports each connection's
before/after free energy and parameter work; this work is neither a heat pulse
nor a mechanical impulse. The damping-only API delegates to the same update.

`FrameVibrationComponent::into_thermoelastic_system` binds the complete original
material frame to its same whole-layer caloric cells. The original complete
mass-normalized basis is a fixed coordinate chart, with native free-coordinate
kinetics and three positive physical springs per source layer. Each axial
spring has that layer's integrated thermal rest extension; the two bending
factors preserve the original full stiffness projection. Separate layer springs
retain mismatch energy for unequal temperatures/expansion coefficients.
Off-diagonal K(T) is retained, not replaced by scalar modal frequency retuning.

Each accepted-temperature interval uses the SAME native coupled propagator.
Original spring-port work is checked against the source potential, and original
material viscosity-port losses are partitioned by their actual layer
rigidities. Each layer retains the original thermoelastic Helmholtz law:
S = -dPsi/dT and U = Psi + T*S. The finite thermal residual receives separate
viscous heat and signed reversible exchange
Q_rev = -delta(Psi_at_fixed_endpoint) - delta(T*S).
That expression is independently checked against stress work minus delta(U);
native spring-work roundoff is not deposited as heat. Individual native
coefficient rows must reproduce each layer's complete before/after storage and
parameter work, followed by the combined mechanical/caloric energy account.

Thermal trials query the same pinned sources for E(T), alpha and eta(T). Their
jointly accepted endpoint prepares stiffness/rest/damping together before any
history publishes. Source, thermal, strain, native budget or cancellation
refusals retain motion, coefficients, finite heat and the shared clock for
deterministic retry. Reports retain original caloric entropy and mechanical
source entropy, with thermal and thermoelastic temporal defects distinguished
from physical material-loss entropy.

This first-order freeze-then-refresh rung requires complete source-bound
members, all original free modes, fixed reference mass/geometry, no massless
condensation or frozen preload, no additional authored modal damping, scalar
expansion and scalar or globally affine E(T). All layers remain solid with
positive combined caloric/thermoelastic capacity and source support through the
sampled frequency band. Fibre-strain gates apply to accepted/proposed endpoints;
there is no continuous intersample strain certificate. Radial constraints,
buckling/geometric stiffness, changing mass and melting/topology need their
corresponding physical owners and are not inferred.

Manufactured regressions use independent axial/bending oscillator solutions
and a heterogeneous two-DOF consistent-mass physical reference, with four
unequal-temperature/alpha layers and simultaneous E(T)/eta(T). Both native
backends, actual continued motion, off-diagonal stiffness, layer U/S, caloric
closure, and nonzero-prefix refusal/retry are covered.

## Two-axis modal shear and same-step finite-body thermal feedback

`friction::tangential::planar::heat` connects the original
`PlanarTangentialModalSystem<L>` to its original `FrictionHeatRoute` and finite
`DepositionSimulation`. The existing enthalpy/Newton owner solves original
cell-energy increments and a common interface temperature. At each trial,
that temperature reaches the same retained constitutive law and its accepted
two-axis history. Both force components use the full original 2-by-2 modal
response, including cross-axis compliance; one complete network proposal
provides both actual displacements and the mechanical energy account.

Only native irreversible friction work enters the original caloric cells.
Reversible shear energy in both directions stays mechanical, and original
network damping remains in its own loss account. Geometry and actual endpoint
conductivity determine signed heat flows, including ordinary hot-to-cold
conduction at exactly zero slip. No heat fraction, replacement thermal store or
extra friction deposition is introduced. Friction transfer cancels once in the
combined mechanical/caloric balance.

The final nonlinear root is restaged after finite-difference probes. Original
mechanical compatibility, constitutive work, thermal balances, source admission
and cancellation checks all precede publication of bodies, shear history,
sample clock and caloric state. Source callbacks may return transport receipts
but cannot replace the law or publish trial history. Frozen-temperature elastic
sources retain their refusal; variable elastic free energy requires an explicit
reversible heat/entropy port and is not inferred by this dissipative adapter.

Five focused groups in `tests/modal_planar_tangential/heat.rs` check independent
analytic thermal softening with rotated anisotropic inertia, both force axes
and both bodies' recoil, transverse elastic storage and return work without
heat, conduction at zero slip, loaded-history late refusal/cancellation/retry,
and refusal to rebind a frozen elastic source. The analytic free-body fixture
retains a nonzero off-diagonal response and finite thermal capacities, so an
isotropic mass replacement or prescribed-temperature friction cannot satisfy
its reference. Validation uses the original source/dependency implementations.

Scope: one fixed tangent frame and prescribed normal load, original held-effort
modal mechanics and backward-Euler conduction. Both time integrations require
refinement. No moving footprint, evolving normal geometry, temperature-dependent
elastic entropy closure, measured calibration or full workspace build is claimed.

## Nonlinear spatial conduction and finite material heat share one endpoint

`ConductionThermalNetwork::step_response_candidate` places the original
`PlanarConductionResponse` inside the SAME material network and finite-cell
enthalpy residual. Each proposed thermal endpoint refreshes the original
temperature-dependent conductor bundle, then closes native |E|/|J| through
the original FEEC solve. Actual local E.J power enters the original receiving
cells exactly once. Terminal work, constitutive power defects and outer
thermal residuals remain separate; no iteration defect becomes physical heat.

The original whole-volume routes, caloric instances, material identities,
selected source roles, fixed external conditions and thermal links remain
bound. Continuation retains accepted spatial operating coordinates even when
the original thermal reference query used zero field. Every inner native
field solve consumes the existing shared work budget, including discarded
constitutive trials and fresh final thermal admission. A refused source,
budget, stencil, terminal gate or cancellation publishes neither field/source
nor thermal history.

This scalar, memoryless, fixed-geometry rung supports native voltage and
distributed current-density traces. A nonlinear equipotential total-current
terminal requires its own constraint and refuses here. Thermal conductivity
response/applicability must be independent of E/J. The original caloric chart's
retained Cp, density and melting-transition receipts are inspected through their
EXACT selected claims: E/J-dependent value or applicability refuses before
using an incorrectly frozen heat-storage law. Categorical phase-branch
selectors remain legal, and field-name fallback is scoped to the actual
receiving cell. Explicitly authored h(T,phase) remains the declared independent
caloric law; electrocaloric storage/work is not inferred. The source adapter's
existing joint phase/E/J admission limitation also remains explicit.

Independent two-capacity balances exercise nonlinear voltage and current-density
laws, nonuniform local heat, data-only material swaps, thermal feedback and
continued accepted state. Regressions include nested work limits, late
cancellation, source-domain exits, mixed old/new endpoint refusal, conductivity
applicability, and Cp/rho/latent-heat response/applicability guards with a valid
ordinary melting-chart control.

## Complete Kelvin-Maxwell plane-stress dynamics

`PlaneStressPoint` admits the original Kelvin tensor when its normal block is
positive definite, retaining three massless normal-strain coordinates after all
six coordinates of each Maxwell branch. The normal traction constraint is
(D_{nn}\dot e_n=-H_{e_n}-D_{nu}\dot u). Its exact elimination contributes a
skew momentum/normal-strain coupling, normal resistance (D_{nn}^{-1}), and
in-plane resistance (D_{uu}-D_{un}D_{nn}^{-1}D_{nu}). The latter is assembled
as a congruence through the existing material factor projection, never by
subtracting large loss matrices or fitting modal damping. The existing complete
port-Hamiltonian exact-ZOH owner integrates mechanics, internal history, source
work and heat with either diagonal or consistent physical inertia.

Purely in-plane Kelvin viscosity keeps the existing algebraic normal-strain
condensation and history layout. Zero Kelvin viscosity is unchanged. A nonzero
singular normal-viscosity block explicitly refuses; no normal inertia, penalty,
pseudoinverse, eigenvalue clipping, branch deletion or substituted material law
is introduced. Membrane capacity includes every retained normal strain before
point preparation. A nonzero instantaneous displacement jump with normal
Kelvin viscosity has no finite-work prehistory and refuses; relaxed
initialization and a zero-displacement virgin history are admitted.

Accepted-history traction uses actual velocity and includes the original elastic,
Maxwell and Kelvin stress. The velocity-free point traction API refuses when
normal Kelvin viscosity is active. Strain-domain, cancellation and accepted-state
publication remain owned by the existing material history.
The pure-Maxwell `PlaneStressRelaxation` image explicitly refuses nonzero Kelvin
viscosity; complete Kelvin evolution remains available through the original body.

Three new regressions cover full-mass dynamics against independently integrated
physical traction-rate equations, force work and Kelvin/Maxwell heat; the
analytic purely in-plane damped oscillator; and quadrature subdivision,
initialization/capacity refusals and cancellation. Together with the existing
plane-stress, exact relaxation-image and membrane checks, the focused native
target passes 20 tests.
This establishes these bounded linear-kinematic material equations, not a
general singular DAE, finite-strain shell or full-workspace qualification.

## Complete radiation reaction on nonlinear plate bodies

`thin_plate::VkBody::with_radiation_load` binds the existing shared passive acoustic `Model` to every original mass-normalized modal velocity. Every mode has a model column, including zeros; all signed mutual terms survive. The original point-force and area projections remain independent. The existing acoustic modal bank and explicit physical-port interconnection append the complete acoustic storage and its reciprocal force to the same exact/Gonzalez equation. The existing prepared-body API applies unchanged.

Accepted mechanical, material and other appended host histories stay an identical prefix; only the new acoustic coordinates begin at rest. A moving attachment explicitly starts a new quiescent exterior at that instant, rather than inferring past acoustic memory. Total loss includes acoustic damping once and must not be credited again as material heat. Every physical exterior needs a single ownership declaration. The original complete dense-render cap and caller pole allowance refuse oversized requests without discarding modes or histories. Cancellation is checked around the bounded original setup owners.

The producer remains responsible for geometry, gas, basis and fitted-band validity. A passive authored pole model is not acoustic-fidelity evidence. Existing compact pressure observations remain distinct from a finite-wavelength receiver transfer. Three focused native G1/G3 tests compare complete linear/nonlinear, full-damping, Maxwell-memory and two-pole mutual-radiation trajectories with an independently assembled velocity-space RK4 system, independently check stored energy and force work/heat, measure second-order time convergence over three clocks, and cover moving-history preservation, zero loading, capacities, invalid physics and cancellation. Native tests use original numerical owners and a narrow source-extracted plate/model/cancellation seam; they do not establish full workspace or instrument qualification.

## Native modal mechanics driven by original liquid-bridge work

`render::schedule::force::coupled::contact::capillary::CapillaryModalSystem`
couples the original planar fixed-volume liquid-bridge source to the complete
native mechanical network. The gap column and its transpose map actual
mechanical motion and work-conjugate mean force reciprocally. The original
network's held-force compliance enters the bounded bridge load-line solve;
all native modes, connections, force limits and propagation remain active.
The two receiving substrate identities must match the ordered original
wetting sources, and the liquid source is borrowed throughout.

Every step stages the actual native motion and reevaluates the original
bridge on that motion. Gap and constitutive-force gates precede an independent
combined mechanical/surface energy account. The signed surface potential
retains its dry-substrate reference; prescribed opening work, external work,
native dissipation and topology release are separate. A final cancellation
check precedes joint publication of body, bridge and observation successors.
The evaluation budget reserves this actual-motion constitutive check;
runtime scratch is preallocated. A refused proposal retains all accepted
states and observations. The fixed-direction bridge remains the declared
planar quasistatic model with its own rate/geometry screens; native body
inertia does not add resolved liquid inertia, 3-D end curvature or pinch-off.
No topology release is automatically deposited as heat.

## Objective cohesive transport and finite-rotation work

`cohesive_contact::objective::ObjectiveCohesivePoint` retains the original
cohesive point, bonded area, stress-free material offset and explicit material
frame. Common rigid transport preserves the complete material history. A
physical trial uses the original timed analytic local-path mean traction;
an explicit local-jump route represents zero local deformation under rotation
without clamping a computed strain or altering the constitutive law.

The spatial discrete force pair and balancing frame couple obey the finite
work identity using the arithmetic average of endpoint rotation matrices,
the spatial Cayley increment and midpoint separation. Vector kinematic and
scalar work checks independently refuse an inconsistent proposal. Endpoint
wrenches remain distinct from these discrete work efforts. The mechanics
owner must apply the additional couple to the frame-carrying side and use
compatible rotation increments and moment arms, or audit its actual work;
this adapter does not certify an arbitrary rigid-body integrator. Physical
half-turn/chart-limit steps require subdivision. Pure common-frame transport
does not advance material damage or viscosity.

`objective::source::SourceBoundObjectiveCohesivePoint` retains the original
resolved interface by reference, checks ordered receiving materials and exact
constitutive-law pointer identity, and retains that binding in successors.
Fixed source conditions and caller-supplied attachment/frame geometry remain
explicit. Numerical and source refusals leave accepted history unchanged;
fracture work and viscous loss keep their original separate accounts.

## Reciprocal cavity loading of complete nonlinear plates

`VkBody::with_cavity_load` connects the original `HelmholtzModalCavity`
pressure/flow owner to the complete nonlinear plate, full damping and material
history. The supplied signed cavity-side shapes are independent of both point
drive and exterior observation: outward surface motion gives U_in = -A.v,
while cavity pressure applies +A.p. The shared explicit-column interconnection
therefore closes internal work in the same exact/Gonzalez step.

Every accepted host coordinate remains unchanged. Two appended coordinates retain
gas compression volume and inward acoustic momentum, initialized from explicit
pressure and outward neck flow. Compliance, neck inertia, end corrections and
selected resistance come from the existing source owner once. The returned
`VkCavityLoad` observes accepted ordinary or prepared-body state without
advancing it or applying a second reaction. Gas/neck storage is already part
of the body's total energy; neck loss is already part of total dissipation.

Three focused G1/G3 regressions compare complete linear/nonlinear,
full-damping/Maxwell/cavity trajectories with independent pressure/flow RK4
equations and a three-clock nonlinear convergence ladder; check independent
work/heat, signed surface cancellation and pressure acting on an undriven mode;
and check exact linear subdivision, state preservation, capacity and rollback.

The original uniform small-signal cavity, fixed geometry/gas and selected neck
loss approximation remain explicit. This adds neither cavity standing waves
nor broadband neck/exterior radiation. Native checks use the actual gas and
Helmholtz coefficient owners for authored-resistance fixtures; unused Bessel
loss branches refuse in the narrow native adapter. They do not establish full
workspace, source-material or measured-instrument qualification.

## Resolved fluid pressure controls through the original NSF/DG path

`NsfDg::step_with_boundary` and the typed `NsfDgPerformance` constructors
accept an original NSF incoming characteristic or prescribed pressure offset.
The noncoupled face holds that physical control while the existing scalar
interface solve matches the other face to the live DG pressure/flow response.
No additional flow law, jet source, force integrator or material snapshot is
introduced. The same trial owns fluid populations, boundary work and acoustic
state; failed steps publish neither subsystem nor new boundary controls.

Native `boundary pressure COUNT` records hold Pa relative to the original face
reference pressure. `boundary characteristic COUNT` and legacy `boundary COUNT`
hold m/s. Source-clock events and conservative observation-clock projection
retain the existing timing semantics. `dg_molecular` now reaches the resolved
reader through the shared source parser and material DG owner; it adds no
molecular internal-energy model to the NSF volume. Species/source-band refusals
remain original material refusals.

`resolved-dg-pressure.performance` composes the prescribed pressure, stationary
resolved constriction, DG bore, computed radiation and exterior gas observer.
Focused tests cover physical flow sensitivity at fixed pressure, returning-bore
feedback, direct original-owner replay, typed/legacy characteristic identity,
callback partition, source molecular memory and transactional refusal. This is
an exploratory fixed-composition, local-linear boundary capability; it does
not qualify a measured flute, a finite player airway or nonlinear transparency.

## Simultaneous retained modal contacts and shared finite-cell thermal feedback

`friction::tangential::multiple::heat` connects the original
`MultiTangentialModalSystem<L>` to one `FrictionHeatNetworkStep`. Every physical
interface has its own massless common temperature; all original finite cells
share one enthalpy residual. Each thermal trial solves all retained contacts
against their complete original `2N`-coordinate modal response, then stages one
sum of body efforts. Contact-to-contact compliance, both tangent directions,
original network losses and original constitutive histories remain active.

`ModalTangentialHeatRoute` explicitly binds a mechanical construction index to
an original thermal link and its numerical temperature guess. Bindings can be
enumerated in any order. Duplicate, missing or out-of-range contact indices and
duplicate link ownership refuse before solving. The thermal owner and source
callback use ascending original link IDs; source conductivities must carry
those same identities. Mechanical frames retain construction order, and the
returned index map relates each thermal report to its original contact.
The caller retains authority for the physical modal-body/cell correspondence.

Only native irreversible contact work heats the interface nodes. Geometry and
actual endpoint conductivities determine signed body heat. A shared cell keeps
one original capacity and receives each adjacent interface contribution once;
ordinary interbody conduction continues at zero friction. Reversible shear
storage and original network damping retain their mechanical accounts. The
combined balance cancels friction transfer internally. An optional final-root
audit reuses original exact caloric entropy and reports each interface's
physical friction and conduction generation without adding another solve.

Cell count, interface count, joint unknown count, nonlinear iterations and
total source/mechanical trials have explicit bounds. Finite-difference probes
start from the same accepted bodies and contact histories. After convergence,
all contacts and their summed body force are restaged at the final temperature
vector. Source admission, original mechanical/thermal balances and cancellation
checks precede joint publication. Frozen-temperature laws retain their refusal;
changing elastic free energy requires an explicit reversible heat/entropy port
and cannot be admitted by replacing a law or resetting its retained history.

Five focused groups in `tests/modal_multi_tangential/heat.rs` and seven existing
multi-contact mechanical groups pass through the original native source and
dependency implementations. An independent two-equation thermal-softening
reference eliminates three finite cells and retains shared body compliance;
it verifies distinct interface temperatures, both force axes, all three body
recoils, shared-cell heat once and exact caloric entropy. Further groups verify
mechanical/contact-route permutation, zero-slip tridiagonal conduction,
duplicate/missing/source-identity refusal and loaded vector-history late
refusal/cancellation with exact retry.

Scope: a retained contact set, fixed tangent frames and prescribed normal loads
per sample. Original held-effort mechanics and backward-Euler conduction require
time refinement. Root convergence establishes neither uniqueness nor continuous
event accuracy. Moving footprints, evolving normal geometry, variable-elastic
entropy coupling, measured calibration and full workspace validation are outside
this composition.

## Native three-dimensional material valves on the shared DG renderer

The material-modal DG grammar accepts `solid tet4 SOURCE_TOKEN` in place of
plate loss/blade records. Original material cards provide equilibrium elasticity,
density, complete symmetric Mandel branch coefficients, relaxation times and
strain limits at declared external states. Explicit material axes, regional
element assignments, supports and original exterior pressure faces reach the
existing solid owners; no effective modulus, loss ratio or pressure area is
invented by this reader. The signed reference pressure map is work-conjugate to
swept volume; positive closing pressure acts against it. Nodal force histories
use the same source-port scheduler and original shape projection.

`solid_memory` selects the existing complete hereditary reduction and aperture,
including actual Maxwell history and relaxed/unrelaxed initialization.
`solid_approximation` explicitly selects the existing frozen-frequency K/full-C
image and its retained-frequency dispersion gates. Every selected mode, including
zero-opening modes, survives either route. Both use the same pressure/contact,
DG gas, radiation, exterior observer and source clocks. Metadata distinguishes
the two physical approximations and retains original source pack identities.

The frozen solid aperture now validates original volume strain at cold and
actual candidate endpoints through the source element maps; modal displacement
or mean opening alone cannot admit an invalid material strain. Failure precedes
accepted motion/clock publication. The hereditary owner retains its existing
source-domain check. These are endpoint checks, not temporal strain enclosures.

The native surface-load owner is P1, so Tet10 input refuses rather than silently
using linear pressure faces. Contact remains the existing effective mean slit;
plate spatial traces, follower pressure, arbitrary authored material memory,
large strain and thermal material evolution are not implied. The supplied
`solid-dg.performance` and `solid-hereditary-dg.performance` are manufactured
small-solid integration fixtures, not measured clarinet/oboe specimens.
Focused regressions exercise an independent tetrahedral compliance/work identity,
material-only changes reaching actual DG/microphone state, complete viscosity
and causal memory, partition identity, source/clock/capacity refusals and atomic
strain rejection. Native execution remains a separate verification requirement.

## Source-bound planar layer thermoelasticity

`friction::tangential::planar::material_layer::PlanarMaterialLayerSystem`
continues the original resolved third-body layer in both tangent directions
on the existing modal network and one explicitly bound finite film cell.
The same `ResolvedShearLayer`, source/card/caloric admission and native
`ThirdBodyLayerContact` own the material and elastic/plastic history. The
original full two-axis network response and combined body effort are reused.
The thermal cell must retain its actual mass rho*A*h, reference density,
original caloric object and solid phase. No receiving capacity is invented.

The original source owner admits one affine G(T) over the complete caloric
interval and checks static and sample Nyquist coverage. Its stiffness slope
is cached once; every endpoint refresh pins the same selected claims, geometry
and elastic domain. The initial source is canonicalized to pinned selection
without mistaking a changed policy receipt for changed physical coefficients.
Caller-supplied temperature refuses, and any supplied pressure must equal N/A.
Compression is retained; the native law remains pressure independent.

For psi=k(T)|e|^2/2, mechanical entropy is s=-k'|e|^2/2 and internal energy is
U=psi+T*s. The original isothermal mechanical trial uses endpoint coefficients
and accepted strain. The original cell receives plastic loss minus coefficient
free-energy work minus delta(T*s). This signed reversible exchange is never
classified as frictional heat. The returned report retains plastic heat,
reversible exchange, coefficient work, entropy energy and actual delta(U)
separately. Original other heat and conductances retain their thermal owners.

Each bounded nonlinear or Jacobian trial starts from the same accepted two-axis
history and original finite cells. Final restaging, cell/network acceptance,
a direct endpoint internal-energy balance, the reconstructed mechanical balance
and cancellation checks all precede joint source/law/history/body/thermal
publication. No additional accepted history or thermal reservoir is created.

Five focused groups in `tests/modal_planar_tangential/material_layer.rs`
verify an independent loaded two-axis thermoelastic root under heating and
cooling, nonproportional plastic motion and complete physical energy balance,
constant-modulus recovery of original mechanics, source/geometry/domain refusal,
and loaded late-failure/cancellation with exact retry. All five and the eleven
prior planar mechanics/heat groups pass through the original native source
and material, numerical and thermal dependency implementations.

Scope: a fixed homogeneous film and tangent frame, small elastic shear,
prescribed compressive load and affine temperature-dependent stiffness.
Endpoint thermal feedback is first order and needs temporal/spatial refinement;
root convergence does not establish uniqueness. Expansion, phase change,
moving footprints, film advection, measured calibration and a full workspace
validation gate remain outside this composition.

## Complete source tensor memory in the shared string/body network

`TensorPlateRelaxation` binds the same original `TensorPronyProperties` used by
material valves to every retained DKT body material point. Original equilibrium
elasticity, density, material axes, complete 21-coefficient Mandel branches,
relaxation times and strain limits remain source-owned. The equilibrium section
must reproduce the existing plate stiffness. No additional modal damping or
frozen viscosity may duplicate this material storage and loss.

The original plane-stress reduction retains all six history directions per
branch and area point, including histories with zero mechanical coupling. Its
complete relaxation coordinates enter the existing exact string/body network;
bridge reaction, friction, cavity and radiation use the same actual state.
Uncoupled histories retain their own energy and heat and belong to exactly one
connected mechanical island. Explicit mode, material-point, history, projection,
coupled-spectrum and clock limits refuse unsupported work without truncation.

`ModalPlaneStressPoint` retains original displacement/history ranges through
cold construction, component extensions, static preload and cold restart.
Original complete surface strain and normal traction are reconstructed from the
actual candidate coordinates and histories. Source admission follows the final
contact, friction, DG or thermal solve, immediately before publication; numerical
free-motion predictors and nonlinear trial states are not accepted physical
endpoints. Refusal leaves accepted mechanics, histories and clocks unchanged.
The surface mapping follows the original homogeneous bending quadrature; these
are endpoint and area-quadrature checks, not continuous space/time enclosures.

The native `loss tensor-maxwell` body record names original strain, branch and
time properties and selects relaxed or unrelaxed history. `clamped-x-min`
supports an explicit one-edge cantilever without changing perimeter-supported
profiles. `tensor-body.body` and `tensor-body.material-strings` retain two body
modes, six original area points and 36 material histories within existing exact
owner capacities. They are manufactured integration fixtures, not calibrated
violin, cello or bass bodies. Source thermal evolution, shell thickness modes,
large strain and arbitrary initial point histories are outside this adapter.

Authored regressions cover tensor-only source changes reaching bridge feedback
and PCM, complete normal-history reconstruction, block partition identity,
source/clock/capacity refusal, zero-coupling decay and unique island ownership,
and atomic strain admission after contact, preload and restart. Native Cargo
execution remains a separate verification requirement.

## Native three-dimensional capillary adhesion

`render::schedule::force::coupled::contact::capillary` admits the original
planar and axisymmetric bridge models through one sealed, statically dispatched
native mechanical transaction. Existing planar construction remains available;
`AxisymmetricCapillaryModalSystem::new_axisymmetric` binds the original 3-D
bridge source and its work/accuracy budgets to actual modal bodies.

The owner stores the actual mechanical gap separately from the numerical
equilibrium profile gap. Root residuals therefore do not accumulate as physical
geometric drift. After the original network is staged, the original bridge is
reevaluated at actual displacement. Compatibility, estimated gap/force error,
actual constitutive force, total energy and cancellation must all pass before
body state, bridge state, physical gap and observation publish together.
One constitutive evaluation is reserved for this final actual-motion audit.

Signed surface energy, external and prescribed work, network loss and planar
topology release remain distinct. The axisymmetric model retains its original
quasistatic, fixed-topology applicability and does not manufacture liquid
inertia, rupture or heat. The six focused original-source native tests cover
both models, moving bodies, spring/damper connections, prescribed/external
work, atomic refusals and retry; they are not full workspace validation.

## Native DMT adhesion with two-dimensional junction sliding

`render::schedule::force::coupled::contact::adhesive_dmt::AdhesiveDmtModalSystem`
projects all three original reciprocal contact columns through the actual
mechanical network and couples them to the original sourced DMT and junction
shear laws. Normal and both tangent components share the full compliance;
normal motion changes the actual area that carries static/kinetic shear.

The native owner audits the original constitutive law on actual staged motion
and separately accounts for normal storage, friction loss, external work,
prescribed work and network dissipation. Initial and staged endpoint velocities,
including prescribed support motion, determine occupied kinetic state. Zero
net slip cannot silently restore static strength while the endpoint is moving.
Only after force, work, energy, material and cancellation gates pass do body,
contact, clock and frame advance together. Source mismatch or refused trials
leave accepted history unchanged.

Seven focused original-source native tests cover coupled body motion, momentum,
complete work, off-diagonal compliance, prescribed motion, gap-only attraction,
occupied kinetic continuation, a zero-net-slip reversal, timestep refinement
and atomic retry. This adapter retains the straight-path constitutive model;
it does not resolve intrastep stick/slip events, evolving contact frames,
partial-slip elasticity, JKR peeling, heat destination or bond healing.

## Frictional reclosure of completely failed cohesive interfaces

`cohesive_contact::reclosure` continues the original fully failed material
point with its undamaged compression penalty and separately sourced residual
static/kinetic Coulomb coefficients. Ordered interface identity, original law
and irreversible fracture history remain attached. A wholly open path carries
no friction, and reclosure does not restore cohesion.

All entries of the reciprocal positive-definite 3-by-3 compliance enter the
normal/two-tangent solve. Exact straight-path mean compression supplies both
normal work and residual friction capacity. A conservative sufficient feedback
screen rejects configurations whose scalar sliding monotonicity has not been
established. Fixed-size arrays and bounded continuation avoid runtime allocation;
floating-point residuals do not claim interval or global uniqueness proofs.

The existing `CohesiveModalSystem::enable_frictional_reclosure` performs an
in-place source-checked handoff after complete failure. Original bodies, source,
sample clock and prior observation survive. The native owner reevaluates the
original material at actual staged motion and publishes jointly after force,
work, total-energy, material and cancellation checks. It keeps residual friction
loss separate from fracture dissipation and uses actual endpoint velocity to
retain kinetic occupancy when zero net slip conceals a moving endpoint.

Nine focused original-source tests cover full coupling, manufactured work,
entry/exit, tangent rotation, source/history checks, actual rupture-to-reclosure
handoff, velocity occupancy and atomic retry. Partial-damage friction, rebonding,
evolving sources, evolving frames, heat recipients and exact stopping events
remain outside this owner; a full Cargo/workspace check is separate.


## Live thermoelastic frames in complete native assemblies

`into_thermoelastic_parts` supplies the original free-coordinate frame, material connections and a source/caloric binding for an explicit component and contiguous connection range. The existing standalone constructor delegates to this path. Binding admits a cold, unchanged source frame within the original fully configured `CoupledModalSystem`; signed columns, coefficients, initial state, thermal clock and source temperatures must agree. Native companions, reciprocal acoustic loads, initialized Maxwell memory, propagation and budgets remain in that same owner. Source K/rest/C updates use their original global addresses. Only frame material dashpot loss heats its finite cells; foreign connector loss is recovered from its native work and spring-storage change and checked against the complete loss account. Focused physical-coordinate ODE, acoustic refinement, finite-heat, memory and rollback regressions cover this coupling. This slice supports one live thermoelastic frame, fixed mass/geometry and its complete basis; warm rebinding and changing phase remain unsupported. Prepared acoustic propagation retains its existing second-order averaged interconnection.


## Phase/field source laws under total-current electrothermal control

The original finite-enthalpy spatial response candidate now closes actual E/J and caloric T/liquid mass fraction under voltage, current-density or equipotential total-current drive. Both ordinary and nonlinear held-current paths delegate the same native electrode solver; every thermal trial enforces its authored current inside the source closure. Supplied work uses authored I times solved electrode voltage independently of recovered reaction, while local native E.J remains the deposited heat. All trials share the existing field-solve budget, pinned caloric/source admission, final fresh acceptance and atomic histories. Zero current retains conducting-terminal admission and never invents a potential across disconnected islands. Independent finite-capacity and latent-balance regressions cover both signs, data-only changes, partial/full-melt continuation, zero current and budget/cancellation refusal. Fixed geometry, field-independent caloric storage and bounded memoryless scalar response remain explicit limits.



### Directional Kelvin–Maxwell plane-stress constraints

`PlaneStressPoint` admits Kelvin viscosity acting on any subset of the local
normal tensor-strain directions [zz,yz,zx]. An eliminated direction must have
an exactly zero row in the original Kelvin tensor and its projected matrix;
the retained principal normal block must be positive definite and numerically
resolved. Every nonzero direction is retained, including small positive
diagonal entries that can be resolved. No rank tolerance, pseudoinverse,
penalty stiffness or dropped Maxwell factor row is used.

Undamped normal strains minimize the complete instantaneous equilibrium plus
Maxwell energy, including coupling to retained normal strains. Damped normals
satisfy the original massless traction-rate constraint. Their mobility and
momentum cross-coupling retain symmetric resistance and skew reciprocal work
in the existing fs-phs runtime. Physical inertia, external force ports, all
six coordinates per original Maxwell branch, reconstructed normal tractions,
and one accepted work/heat clock remain intact. Capacity counts the actual
zero-to-three dynamic normal coordinates at point and membrane admission.
Initial histories retain the prior finite-work rules; pure-Maxwell relaxation
images still refuse Kelvin viscosity.

The focused 22-test native target passes, including all six one-/two-direction
cases against an independent full-material DAE reconstruction with LU
constraints and RK4 in physical displacement, velocity, original Maxwell
memory and normal strains. Tests compare consistent inertia, signed forcing,
all normal tractions, supplied work and Kelvin/Maxwell heat. A 1e-24 positive
viscosity stays dynamic, and exact state-cap boundaries are exercised.
Existing plane-stress, membrane, loading and relaxation tests also pass.
This is not full-workspace Cargo/DSR qualification. Singular oblique normal
blocks without exact coordinate-aligned zero rows continue to refuse when
their retained block cannot be resolved; no general rank-revealing DAE solver
is claimed.

### Resolved source with a finite upstream acoustic load

`NsfDg::with_upstream_load` attaches a second original `DgPortLoad` to the
upstream face of the same `NsfZone`. Both duct inlet orientations point out of
the zone; a prescribed series pressure rise satisfies
`p_zone - p_upstream = rise`, with external work `-rise * upstream_flow * dt`.
Neither face flow nor zone pressure is independently prescribed. The original
NSF face maps, downstream observer, material wall and molecular histories remain
their owners. Matching gas, face area, clock, zero reference flow and cold state
are admission conditions; duplicate or independently driven attachments refuse.

The two characteristic solves share one evaluation/cell-step allowance and
reserve a final evaluation to restage both ducts and the original source from
the same accepted state. Pressure residuals, covariance, mass projection and
both original port-work accounts must pass before all three states publish.
Failure or cancellation leaves their histories and clocks unchanged. The
one-load boundary-pressure and boundary-flow paths retain their arithmetic.

The native `upstream_dg`, optional `upstream_molecular`/`upstream_wall`, and
`boundary upstream-pressure` records reach this owner without a surrogate
renderer. Frame and CLI provenance retain the drive kind, original upstream
duct geometry/material identity and separate source work. The 96 kHz Lab
starter is an acoustic feedback experiment: the DG interfaces do not advect
mean mass or enthalpy, and the finite duct is not a qualified anatomical airway.

### Continuous passive piston loading with resolved weak resistance

`thin_plate::BaffledPistonRadiation` uses the original Bessel–Struve circular
piston impedance in an infinite rigid baffle and the evaluated coherent gas.
Its port is pressure/volume flow; the host supplies the same signed volume
shape for flow and reciprocal reaction. The complete shared acoustic modal
bank retains every fitted q/p pair. This is a source/composition adapter, not
another time integrator or a conversion of legacy sampled-filter histories.

The opt-in `radiation::fit::fit_resolved_relative` keeps the existing PSD block
solver, adds frequency-local complex/resistance weights and a logarithmic
enrichment of the first resolved dictionary interval, and balances dictionary
columns without changing the objective or cone. Even frequency samples alone
determine coefficients and convergence. Training and held-out samples must each
satisfy a 2% pointwise complex and 2% Hermitian-part relative limit, in addition
to the original global peak/RMS gates. Legacy fit APIs and their objective are
unchanged. Errors are sampled Frobenius norms, not continuum or directional
multiport power certificates. Zero local references or unresolved dynamic
range refuse instead of being hidden under a global peak normalization.

The piston adapter uses an explicit positive frequency band, an odd 33–257
point quadratic grid, a geometry-resolving dictionary and a complete 1–128 pole
capacity. Outside-band queries, incoherent gas, unresolved sampling, excess
history and cancellation refuse. Density/geometry remain physical inputs; no
coefficient is fitted to sound recordings. Bodies with multiple independent
faces need separate sources or an actual mutual-impedance producer. There is
no claim for unbaffled, nonuniform-face, off-band, or full instrument radiation.

Nine focused native tests pass: four new piston physics groups and five
existing shared-fit regressions. Independent linear/log frequency grids for
0.01/0.05/0.15 m radii over the tested 40 rad/s–0.4*pi*48000 rad/s band show
complex errors below 0.7% and resistance errors below 1.4%, including the weak
low-frequency resistance. Tests compare complete coupled state, work and heat
with independently assembled physical-velocity RK4 equations under both normal
signs and zero coupling, check a loaded small-ka pole/decay limit, and exercise
density scaling, held-out weak-loss corruption, capacities and refusals. These
use actual numerical/gas/piston/fit owners with narrow source extraction; the
large BEM fixture and whole-workspace Cargo/DSR qualification were not run.

## Native moving-body adhesion with dynamic wetting and pinning

`capillary::ContactLineCapillaryModalSystem::new_contact_line` binds the
original dynamic-footprint bridge to the existing native modal transaction.
Actual bodies, springs, dampers, external forces and prescribed support motion
share the original network compliance. Contact-line relaxation changes the
mechanical force, while pinning retains the physical footprint until a sourced
advancing or receding threshold is reached.

The original wetting law is replayed on actual staged gap motion before
publication. Body motion, footprint history, physical gap, clock and frame
advance together only after source, compatibility, force, energy and
cancellation checks. Source geometry/work uncertainty must fit inside the
combined energy tolerance; a small algebraic residual alone is insufficient.
This energy-accuracy admission also retains the original axisymmetric bridge's
reported work uncertainty. Surface storage, external/prescribed work, native
network loss, physical contact-line loss and topology release remain separate.

Five focused original-source tests cover reciprocal moving-body mechanics,
physical pinning, sourced rate dependence, spring/damper and external/prescribed
work, numerical energy uncertainty, and atomic refusal/retry. The six existing
planar/axisymmetric native tests also pass. The new path retains the original
symmetric two-dimensional, connected, contact-line-dominated liquid regime;
it does not resolve bulk liquid motion, pinch-off or thermal deposition.

### Reciprocal PlateBank evolution with complete material and acoustic history

`thin_plate::PlateBank::advance_controlled` composes the original compact-body
mass, signed force and volume projections, optional full nonlinear `VkBody`,
uniform-pressure Helmholtz cavity and declared exterior loads into the existing
`PreparedVkBody` owner. It retains every material and acoustic coordinate, full
modal damping and nonlinear potential. The original exact quadratic or Gonzalez
equation owns the single accepted step; this adapter introduces no integrator.
`PlateBankStep` reports interval-mean compact pressure, endpoint volume flow,
complete storage change, physical input work and total physical dissipation.
The total includes acoustic loss and must not be deposited again as solid heat.

Positive surface area denotes outward motion. Cavity input is `-sum(A_i v_i)`;
its pressure reaction is `+A_i p`. Exterior input is `+sum(A_i v_i)` with the
opposite mechanical pressure reaction. Each exchange uses the same signed map
in both directions, independently of the external point-force projection.
Modes harvested from one plate carry producer-owned face identity and share
ONE uniform-piston history, preserving cross-mode loading. Equal face areas
alone do not identify a shared face; separate faces remain independent loads.
The complete original Vk material/load history is retained without reattaching
or duplicating its existing loads. Automatic exterior attachment applies to the
bank's compact faces; Vk exterior loading remains explicit in its original API.

Continuous piston attachment starts from explicitly quiescent acoustic history
and refuses migration from sampled digital filters. The original fit window
is retained, each active structural frequency must be in that source band, and
the source's pointwise/global sampled fit limits apply unchanged. Gas and face
definitions stay pinned. A finer numerical clock reuses the same continuous
source and complete accepted history; a coarser-than-declared clock refuses.
Repeated same-source attachment cannot reset history. Independent changes to
public body coefficients, topology or accepted mirrors refuse after preparation.
Solver, invalid observation and late cancellation publish no physical state.

Preparation and changed clocks allocate the existing numerical caches; fixed
clocks reuse them. Immutable constitutive sharing preserves stable nonlinear
energy callbacks and does not add wrapper depth on repeated preparation. The
existing dense exponential capacity bounds the COMPLETE joint state (126
coordinates for this single-force adapter); no history or mode is truncated to
fit. A bare bank of independent compact oscillators retains its O(mode-count)
scalar exact path without that dense capacity. Assembly dispatch declares bank
radiation lazily, so paths selecting another physical owner do not prepare an
unused bank model. This does not migrate those other renderers' own loading.

Focused G0/G1/G3/G4 native regressions compare joint motion, Maxwell history,
cavity exchange, external work and heat with independent physical equations;
exercise nonlinear time refinement, signed projections, shared-face loading,
clock subdivision and atomic refusal; and retain the large bare-bank path.
The sharing seam also preserves complete storage and step results across 1,024
cold copies. Full workspace Cargo/DSR qualification remains blocked by the
unavailable constellation siblings. No measured real-time performance, spatial
cavity standing waves, mutual radiation between separate faces, or continuum
frequency-error certificate is claimed by this adapter.

### Molecular material state in characteristic pressure feedback

`driving_point::molecular::MolecularDrivingPoint` connects the original
`fs-duct::molecular::MolecularDuct` to the shared characteristic-line FIR
realizer. The same gas reservoirs and wall diffusion determine storage,
attenuation and dispersion at every sampled frequency; they are not reapplied
as an observer filter. The original frozen molecular sound speed sets the
wavefront-time estimate, while the inlet's equilibrium gas remains the real
pressure/volume-flow reference. The explicitly supplied exterior gas enters
the existing terminal and side-mouth radiation laws. Scattering reflection remains finite at
lossless impedance poles. Source and near-DC band refusals propagate.

`ReedBoreVoice::new_with_molecular_driving_point` feeds that line into the same
reed/contact/pressure stepper and retains its original material and geometry.
Changing only reservoir data therefore changes the internal bore pressure,
before any listener absorption. G1/G3 tests cover the independent analytical
cylinder response, the existing lossless limit and its impedance pole, actual
reed pressure sensitivity, exact block replay and source-domain refusals.
The adapter retains frozen homogeneous linear plane-wave acoustics, explicit
gas shear/thermal resolution, optional source-owned local membrane walls, the
sampled FIR approximation and the existing compact-jet observer. It adds no
nonlinear gas evolution, measured molecular data or experimental instrument claim.

The molecular path also supplies actual flanged outlet flow to the existing
`OutletPressure` observer. Its original chain and radiation load determine
`Qout/p+`; equilibrium and molecular paths share one flow differentiation,
axial-delay factoring and Rayleigh realization. Feedback and receiver retain
the same explicitly bound, pressure-matched exterior. Tests compare molecular
outlet flow to an independent loaded-cylinder solution and exercise reservoir
data changes in physical receiver pressure after the frozen wavefront. The
scheduled adapter primes the observer with the same initial characteristic
wave as the reed. A diagnostic compact proxy requires identical inlet and
exterior gases; the physical baffled receiver can use a distinct exterior.

Finite side chimneys retain the original molecular reservoirs and wall law.
Each open or vented mouth is bound to the flanged compact load in the declared
exterior; the original admittance-mixture owner applies its opening fraction
once. Closed chimneys remain acoustic cavities. A closed terminal with an open
side mouth uses the actual network's near-DC reflection, not a sealed `+1`.
Earlier side reflections stay in the full sampled reflection history.

`baffled_outlet_and_open_holes_with_delay_speeds` obtains terminal and mouth
flows in one original chain solve per frequency and uses the shared outlet
realizer. Mouth flow is not replaced by chimney inlet flow. The retained sweep
has a public 1,048,576-bin budget; it refuses rather than dropping an opening.
Each mouth has its own axial-plus-chimney delay and exterior delay. Native
molecular playback primes every opening from the identical source sample -1,
applies the selected exterior law at each mouth's own receiver range, then sums
pressure. Individual FIR delay admission and omission policy remain binding.
`molecular_open_hole_observers` and native sidecars expose those actual paths.

The observer retains the existing parallel-baffle geometry convention: each
side receiver is translated by `L - x_hole` along x. Compact independent mouth
loads and separate Rayleigh disks are approximations, not an exact radiation
power match. Fractional opening is the original admittance mixture, not a
resolved pad aperture. No nonlinear hole vortex loss or mutual exterior
scattering is added. G1 tests compare an independent finite-line/side-branch
network, and G3 tests exercise source-rate sensitivity, separate wavefronts,
closed-mouth omission and exact block replay.

`ReedPerformance` exposes this same owner through optional `molecular`,
`molecular_resolution` and `molecular_wall` records before `compile_limits`.
The shared material reader evaluates the law against the original ambient
species record. Cone slices, composition blocks, diffusion poles and the
angular-frequency band are explicit data. The same resolved inlet drives the
reed flow, and the same exterior drives terminal loading and observation.
Native sidecars retain the actual reservoirs, frozen sound speed and source
wavefront time. Observer absorption remains a separate declared exterior law.
Homogeneous molecular propagation requires explicit source species. Per-segment
overrides and midpoint field images refuse. Explicitly continuous axial and
steady-generated fields use the field-bound molecular route described below;
the field's own source admits reservoirs even with a legacy ambient record. Existing files without these records keep their prior path.
The runnable `examples/reed-molecular-synthetic.performance` uses authored
synthetic coefficients, not measured N2/O2 material data. Focused file tests
exercise data-only rate changes, physical pressure, exact block replay and
unsupported-source admission. The companion
`examples/reed-molecular-openings-synthetic.performance` supplies fully open and
partly vented chimneys through the same native geometry. Standalone integration
tests in `tests/molecular_reed_performance.rs` cover that physical consumer,
per-mouth exterior paths and observer changes that preserve reed feedback.

### Material-solid surface contact in the shared aperture solve

`ModalAperture::with_solid_surface_closure` maps each retained original solid
mode onto declared exterior Tet4 faces using their signed outward normal.
Original node clearances, area weights and physical boundary-edge lengths
enter the same `ModalPlateClosure`, contact potential, dissipative contact and
positive P1 slit-flow integral as plate valves. Zero-mean modes remain local
gap coordinates. Reported opening and opening velocity use that same spatial
trace. Pressure traction and swept-volume projections remain the original
solid integrals; contact selection does not silently redefine pressure work.

Geometry is prepared once on a cold source and admits parallel fixed-reference
face normals only. Unsupported Tet10 traces, missing or repeated exterior
faces, non-boundary slit edges, incomplete clearances, exhausted projection
budgets and inadmissible penetration refuse. Every trial retains the original
source-strain checks, full viscosity and hereditary material states; failed
contact, source-domain or cancellation checks publish no state. This is
small-displacement surface contact, not resolved moving-channel gas flow.
Native `solid_closure spatial` data and the Lab starter select this exact owner.

## Native event-resolved adhesive stopping and restart

`contact::adhesive_events::AdhesiveEventModalSystem` applies the original
sourced DMT junction events to actual mass-normalized free body coordinates.
A complete fixed tangent column supplies the reciprocal force and exact
effective mass. Both friction impulse and its first time moment propagate
original positions and velocities, preserving common and uncontacted motion.
The original event owner selects actual stops, static support, area-dependent
breakaway and reversals; a native velocity tolerance never selects sticking.

An independent impulse integral supplies friction's actual work on the native
bodies and is checked against the original nonnegative junction loss. The
body subsystem records this signed boundary work. Prescribed normal work is
retained separately in the combined signed surface-energy account. Zero net
slip cannot erase positive friction loss or its effect on body position.

Preallocated candidates use the original modal state budgets and joint
publication methods. Body states, contact history, sample clock, native work
account and observation publish together only after source, compatibility,
work, energy and cancellation checks. Refused trials preserve accepted state
for deterministic retry. Runtime body propagation is linear in the retained
mode count and allocates no new state buffers.

This path admits sample-zero networks of unconnected free coordinates with
zero acoustic transfers and no material memory or radiation. Normal approach
is prescribed linearly within a sample. Elastic networks, normal inertia,
arbitrary vector slip, changing frames and thermal deposition require their
own trajectory/transfer owners and are not inferred by this adapter.

### Tensor-body refinement with exact shared history coordinates

`loss tensor-maxwell-shared` binds the same complete original tensor laws,
source point count, strain limits and relaxed/unrelaxed initialization as
`loss tensor-maxwell`. It selects the material owner's exact equal-time
coordinate map with a separate factor-work budget. Original area and thickness
quadrature, full branch factors and every selected mechanical mode remain.
After static preload, coupled contact and each accepted candidate, the complete
original point histories are reconstructed before source-strain and normal
traction admission. Cancellation or source-domain refusal remains atomic.
The unshared input and complete-history route retain their original behavior.

`tensor-body-shared.body` uses eight triangles, 24 original material points and
144 reconstructible histories under the existing 36-state allowance. The same
canonical material-string performance and source packs select it through the
ordinary body-profile path, including in the Lab. This makes finer source
meshes admissible without raising the exact structural owner capacity or
discarding material history. Manufactured cantilever data, fixed thermal state,
flat-plate geometry and finite modes remain explicit qualification limits.

## Native elastic connections with resolved adhesive stick/slip

`contact::adhesive_elastic_events::AdhesiveElasticEventModalSystem` binds the
original elastic event law to one actual native spring/damper connection and
its complete signed column. Original free body modes, rest extension, stiffness
and Kelvin-Voigt resistance are retained. This exact scalar reduction preserves
common motion while actual spring loading changes the adhesive stopping and
restart history. Normal indentation and the original anchor remain fixed;
external modal forces may change each sample.

The original free and connected event adapters share preallocated native
impulse/time-moment propagation and native account staging. Actual represented
body travel, spring extension and endpoint velocity are checked against the
original event path. Independent integrated total-reaction work closes spring
storage, physical damper loss and junction loss. Original connection ports
retain their actual spring/damper work; friction crosses the body-network
boundary once. Numerical propagation uncertainty must fit inside the combined
energy budget. Full-sample passive force bounds preserve original connection
limits, and rejected steps publish neither body nor contact history.

Six focused native tests check independent connected-oscillator and dashpot
solutions, original ports, stop/hold/restart, center-of-mass motion, numerical
work accuracy and atomic
retry; all six prior free-body event tests also pass. Arbitrary connection
networks, elastic component modes, moving anchors in this native adapter,
normal inertia and vector-slip events remain outside this reduction.

## Native liquid adhesion between unlike wetting surfaces

`capillary::AsymmetricContactLineCapillaryModalSystem::new_asymmetric_contact_line`
uses the same original native mechanical transaction as the symmetric and
axisymmetric bridge adapters. Actual moving bodies and prescribed/external
loads exchange work with both original footprint histories. Ordered source
identities, each surface's viscous/pinning losses and the source's numerical
force/work uncertainty remain in the accepted frame.

Final constitutive replay uses actual mechanical gap and the original accepted
two-footprint state. Body state, both footprints, physical gap, clock and frame
publish jointly after all existing source, force, compatibility, energy and
cancellation gates. Six new native tests cover asymmetric response, independent
pinning, actual source exchange, reciprocal work/momentum and late accuracy
rollback. The five symmetric-wetting and six existing planar/axisymmetric
native tests also pass. The liquid owner remains a connected planar bridge;
topology release and a thermal recipient are not inferred.

## Friction while adhesive cohesion is only partially damaged

`cohesive_contact::frictional_damage::DamageFrictionLoadLine` adds an explicit
damaged-area rigid-Coulomb specialization: surviving bonds retain the original
cohesive/viscous traction, while damaged area carries independently sourced
compression-only friction. The original material path analytically integrates
damage times compression pressure through initiation, reversals and failure.
It is queried through `trial_over_time_with_damage_weighted_compression`; damage
is not frozen at an endpoint or advanced a second time.

This declared mixture adds no elastic friction storage. Original compression,
bond energy, fracture dissipation and physical viscosity keep their existing
ownership; friction adds a separate nonnegative loss. The full reciprocal
normal/two-tangent compliance and original cohesive Newton solve remain active
under one aggregate material-evaluation budget. This is a rigid-friction member
of a damage/friction model family, not the full elastic-slip Alfano-Sacco law.

`CohesiveModalSystem::enable_frictional_damage` retains original bodies, sources,
clock and canonical damage history. The final audit uses actual body motion
and endpoint velocity before joint publication. Partial-damage friction and
the older full-failure reclosure mode are mutually exclusive. Three material
and eight core/native tests cover independent quadrature/work, evolving damage,
intact/open/failed limits, source effects, 5/10/20/40-kHz refinement, kinetic
occupancy and atomic retry. Nine original reclosure tests also pass. Sources
and frames remain fixed; global softening uniqueness and intrastep event
resolution are not asserted by this held-force path.

### Complete PlateBank images for external physical composition

`PlateBank::physical_image` exposes the original complete mechanical, Maxwell,
nonlinear, cavity and continuous-radiation state to another assembly owner.
Cold physical assembly is shared with the direct bank renderer; export does not
first prepare an unused exponential or rebuild a plate from scalar frequencies.
An already moving bank exports every accepted state coordinate unchanged.
The immutable original constitutive callbacks remain shared; the image owns its
explicit state and its caller owns the subsequent accepted clock.

Port 0 retains the original force/drive-velocity projection. Port 1 is pressure
along the signed outward surface map / outward volume flow. These projections
are independent. An interior flow domain and an exterior opposing-pressure
impedance require their respective orientation in the existing transformer.
The new fs-phs external-column seam retains every original J/R/storage entry.
No geometry, material, port map or quadratic law is inferred from probes.

`PlateBankObservation` retains the complete source extent and the original
mass-normalized surface weights. It observes finite complete state prefixes,
including states followed by an appended subsystem, and returns endpoint
volume flow or compact pressure from interval-mean surface acceleration.
It refuses missing histories and invalid observations instead of filling them
with zeros. It advances no state and adds no acoustic loss.

Only an already declared radiation source is realized during export; no exterior
model is silently added to an undeclared moving source. Source gas/clock and
geometry remain checked. Empty banks return no image; a cavity-only island has
no mechanical interface and refuses this export. Cancellation changes no
evolving state. The complete-state capacity and physical model limits of the
direct bank remain, and the receiving assembly must admit its own full state.

Focused G0/G1/G3/G4 native tests cover warmed mixed compact/Vk/full-damping/
Maxwell/cavity history, independent force and volume ports, exact continuation,
independent physical equations for an appended flow load, shared piston memory,
and observation/cancellation refusal. This is a reusable physical composition
API; each renderer must explicitly consume it. Full workspace qualification
remains blocked by missing constellation siblings.

### Complete original plate histories in coupled ODE assemblies

`acoustic_realize::realize_coupled_ode` now consumes the bank's original
`PlateBankImage`. It retains all compact modes alongside the optional Vk body,
the original full damping and nonlinear storage, complete material/cavity/
exterior histories, and accepted initial state. It does not rebuild the plate
from a descriptive constructor or infer its full state extent from mode count.
The ODE duct starts after the complete original prefix. Its positive inlet
flow compresses the fluid: the existing ratio-one transformer receives the
bank's outward-volume port, applies `+A p` on the structure and `-sum(A v)` on
the duct, and leaves the original point-force port external and independent.
The observation uses the same signed normalized areas for all compact and Vk
modes, validates the complete prefix, and introduces no second acoustic debit.
Any cavity or exterior load already owned by the bank remains counted once.

Quadratic bodies prepare the shared `PreparedExactZoh` once at the rendering
clock. Every sample reuses its complete transition and independent heat/work
integrals. Nonlinear storage retains the original Gonzalez equation and
temporal admission; no unit-vector effort probing establishes a linear model.
The string's existing held bridge-force coupling is unchanged. This extension
does not migrate the separate Dirac, blown-duct or reed-specific routes and
does not enlarge the shared dense-state capacity. Candidate-vector allocation
remains, and no measured real-time speed or zero-allocation claim is made.

Three focused G0/G1/G3 native tests cover a warm mixed compact/Vk/Maxwell bank
coupled to an independent physical fluid reference, retention of every warmed
continuous radiation coordinate and its energy, and exact empty-bank duct
reduction. State, work, heat, signed pressure observation, complete prefix and
prepared-clock checks pass using the original numerical and physics owners.

### Rank-free Kelvin-Maxwell plane-stress descriptor

`DescriptorPlaneStressPoint` and `DescriptorPlaneStressBody` expose a complete
alternative for oblique singular Kelvin normal blocks. All three normal tensor
strains remain stored physical coordinates; three algebraic rate efforts per
point enforce the original zero-traction equations. The original full Kelvin
matrix is stamped once as the PSD resistance over physical velocity and normal
strain-rate efforts. It is never inverted, ranked, clipped or regularized.
The original equilibrium tensor, all six spring-memory coordinates per Maxwell
branch, complete SPD mass, signed cross terms and force/velocity ports survive.
This adds neither normal inertia nor a penalty stiffness or another integrator.

The existing `fs_phs::step_descriptor` owns the second-order Gonzalez/midpoint
step. The explicit `history_relaxed_at_rest` initialization solves the original
equilibrium normal block and starts relaxed spring memory and rates at zero.
Moving preloads and warm continuation require `history_from_state` with every
coordinate present and original constitutive normal traction consistent. The
algebraic representatives carry no stored energy and are not an instantaneous
rate interpolant in undamped constraint directions.

State capacity includes displacement, momentum, every Maxwell coordinate,
three normal strains and three algebraic efforts per point. Original strain
domain, finite work/heat, nonnegative computed physical dissipation, normal
traction and independent energy balance are admitted before atomic publication.
Cancellation is polled around the existing bounded descriptor solve and before
commit; interruption inside its Newton solve and allocation-free callbacks are
not claimed. Energy and traction acceptance do not certify frequency accuracy.
Callers must refine time; this material API is not automatically admitted to an
audio render clock and does not use the effort-join-only descriptor pole bound.

Three focused G0/G1/G3/G4 native tests pass. An independent physical DAE/RK4
reference exercises normal viscosity ranks zero through three, including both
oblique singular ranks, a full off-diagonal mass and two nonproportional Maxwell
branches. Motion, external work and heat converge at second order on 20/40/80
steps: finest state error is below 7.5e-10, work error below 9.3e-13 J and heat
error below 3.4e-13 J in the authored fixture. Original normal traction remains
below 2e-10 Pa. Complete restart and changed-clock continuation are reproducible;
late cancellation, inconsistent moving restart and cap22/cap21 refusal preserve
accepted history. These are focused actual-owner native builds, not full
workspace Cargo/DSR qualification or a measured material/real-time certificate.

### Prescribed-flow ducts retain the complete vibrating source

`acoustic_realize::realize_blown_ode` now consumes the same complete
`PlateBankImage` as the coupled string/duct path. Compact, Vk-only and mixed
banks retain their original storage, full damping and material/cavity/exterior
histories, including accepted initial motion. Cylindrical-chain dispatch no
longer excludes a Vk source. A ratio-one transformer joins the bank's signed
outward-volume port to the second duct inlet. The independent structural-force
port remains open with exactly zero applied force, and only the surviving
fluid inlet receives the prescribed flow. Thus a pressure-driven body's
response depends on its original signed surface areas, independently of an
unforced point-drive projection. Observed inlet pressure comes from that same
fluid port, with the complete original plate observer added once.

The quadratic route prepares one shared exact propagator for the full source
and viscothermal duct and reuses it across samples. Nonlinear storage keeps the
shared temporal admission and refuses unresolved clocks; there is no inferred
linearization or new integration owner. Shared dense capacity and existing
per-step candidate allocation remain. Reed/scattering and Dirac routes are
unchanged.

Three focused tests execute the actual complete blown renderer with the actual
coherent gas, transport and Foster duct owners. They check point-drive
independence with nonzero pressure-driven motion, Vk-only and mixed source
participation with signed areas and positive-definite, non-diagonally-dominant
material damping, and audible decay from already accepted Maxwell history
without a new blow. This exercises general material loss beside the duct's
singular passive Foster blocks without removing their physical branches.

## Native mechanics with dynamic axisymmetric liquid footprints

`capillary::AxisymmetricContactLineCapillaryModalSystem::new_axisymmetric_contact_line` binds the original 3-D circular-footprint source to the shared native mechanical transaction. Original bodies, connections, external forces and prescribed gap velocity determine actual gap motion. Source surface energy, circular-line viscous and pinning losses, numerical uncertainty and both ordered material identities remain attached.

The final source replay starts from the old accepted radius and uses actual represented body travel. Mechanics, circular footprint, angle, gap, clock and observation publish together after the existing force/work/source gates. Connected topology contributes zero release; no rupture or bulk liquid dynamics is inferred.

Six native tests cover reciprocal impulse/momentum and complete work, pinned-radius angle change, source-friction sensitivity, original connection/external/prescribed work, cancellation/source/range refusal with deterministic retry, and a late numerical-energy refusal that preserves both bodies and liquid history. All six prior planar and axisymmetric capillary tests pass. The underlying branch remains explicitly symmetric and quasistatic in the liquid.


## Adhesive events through complete linear elastic networks

`render::schedule::force::coupled::contact::adhesive_events::linear_network`
propagates one scalar dry junction through all mass-normalized elastic modes,
modal damping, and arbitrary admitted rank-one spring/damper connections with
their original rest extensions. Sliding uses the held kinetic resistance;
sticking constrains only the contact coordinate while the internal network
continues to move. The changing constraint reaction determines release.

The matrix-free operator uses preallocated workspaces and bounded Taylor cells.
Dense-polynomial event isolation can resolve intervening stops even when the
sample endpoint velocities have the same sign. Positive quadrature accounts
separately for modal damping, each connection loss, and junction friction;
external work, contact impulse, and its time moment are retained independently.
The result carries explicit propagation, projection, and event-location work
uncertainty, which must fit the caller's physical energy budget. This is a
bounded numerical estimate, not a formal certification of all possible roots.

Fresh nonzero velocities are not captured by a sticking tolerance. The caller
owns accepted phase and state; cancellation, unresolved grazing events, or
exhausted work/accuracy budgets cannot publish a partial trajectory. Exactly
stationary constrained states retain the required impulse without traversing
irrelevant high-frequency cells. The law assumes fixed sources, normal
capacity, frame, coefficients, and external loads during a step; nonlinear
material memory, radiation dynamics, vector friction, and simultaneous contact
constraints are outside this owner.


## Sourced adhesive events in the original native network

`render::schedule::force::coupled::contact::adhesive_network_events` binds the
linear-network event owner to the original `CoupledModalSystem`. Admission
checks the ordered material identities, original DMT area and junction shear
strength, complete contact column, and every component/connection budget. It
retains every admitted original mode, modal damping term, spring rest extension,
connection column, and read-only pressure transfer; no reduced surrogate body
or spring replaces the supplied network.

Each candidate replays through the original component state guards and checks
the resulting original storage, reactions, pressure, connection port work, and
loss accounts independently. External mechanical work and junction loss enter
the native boundary account while modal and connection losses remain in their
original channels. Contact impulse and time moment retain momentum and position
effects through stop, hold, and restart. Original energy tolerances include the
event solver's numerical work uncertainty, and conservative interval energy
bounds enforce connection-force and tangent-speed limits between endpoints.

Accepted bodies, phase, frame, connection ports, and clock publish together only
after all checks and the final gate; late refusal and cancellation permit an
identical retry. Normal indentation/capacity is held during a step. Material
memory, acoustic back-reaction, and constant-loss connections with distinct
damping columns require other propagation owners and are refused here.


## Joint finite-patch damage and friction

`cohesive_contact::patch::frictional_damage` solves all three traction components
at all admitted quadrature sites in one load-line residual using the complete
cross-site compliance. The Coulomb disk uses each site's original irreversible
damage and compressive pressure; site areas enter physical force and work once.
The shared local constitutive owner preserves original material identity,
opening/shear history, viscosity, and the analytic damaged-compression integral.
The projection scale conditions the residual and introduces no elastic friction
storage. Rank-deficient compliance is admitted for coincident quadrature or
fewer independent body coordinates; numerical least-squares regularization
does not modify the accepted compliance or constitutive equations.

The native cohesive-patch owner can enable this law on the original network.
It re-audits material/friction work against actual body displacements, preserves
per-site endpoint kinetic occupancy even after zero net slip, and exposes both
site and aggregate friction diagnostics. Original forces, torque, component
budgets, energy tolerances, and source identities remain mandatory. All point
histories, body states, diagnostics, and the clock publish jointly; late energy
refusal, cancellation, and point/iteration/branch budget exhaustion leave the
accepted patch retryable.

This is a bounded local branch solve with fixed frames, source queries and
quadrature during a step, and a straight constitutive path. It supports
progressive peeling, compressed-edge sliding, complete failure and frictional
reclosure without healing history. It does not resolve intrastep friction
events or establish global uniqueness across softening branches. Workspaces
and dense dimensions are explicitly bounded (at most 64 sites); local material
probes are shared across the joint Jacobian rather than replaced by independent
site mechanics.

## Nonlinear spatial current with complete material heat exchange

`PreparedHeatNetwork::step_conduction_response_candidate` closes the original
scalar conductivity/resistivity response in temperature, phase and local E/J
inside the same finite-cell enthalpy residual as sourced thermal contacts,
optional source-integrated k(T) layers and total-emissivity radiosity. Voltage,
distributed current-density and equipotential total-current drives retain their
native spatial solves; every inner field solve consumes the shared trial budget.
Native per-triangle E.J enters the mapped original cells once. The source-generic
candidate retains complete heat samples and publishes the field, material-source
and caloric histories together after the original energy and source gates pass.

Contact and emissivity response AND source applicability must be independent of
local E/J, since these thermal boundaries own no resolved interface electric
field. This admission applies to the new nonlinear consumer; older prescribed
context heat paths retain their existing semantics. Caloric storage and thermal
conductivity retain the original field-independence admission. Closed-enclosure
radiation is internal exchange; its numerical net power is excluded by the
original combined-energy gate rather than treated as external supplied heat.

Fixed geometry and view factors, original whole-volume source routes, independent
caloric laws, bounded memoryless scalar response and first-order endpoint
coupling remain the limits. This adds no moving phase boundary, interface
electrical-field model, unconditional convergence or measured-material claim.
Three focused `conduction_enthalpy::response::heat` regressions cover independent
two-body Fourier/contact/radiation and nonlinear-current balances, current
reversal, data-only electrical/interface/emissivity changes, continued history,
exact no-surface-path parity and budget/cancellation/source-support refusal.


## Source-resolved Hertz islands on shared finite thermal cells — v0

`friction_heat::hertz::multiple::propose_rigid_hertz_island_heat` couples
multiple original Hertz patches to one original native body island and one
finite enthalpy network. Each stable patch identity binds its pinned original
coefficient source, unchanged normal load/radius, one distinct original thermal
route and its own solved common-interface temperature. Explicit pressure and
the T/p-only source-axis refusals remain intact. Conductivity inputs retain
canonical original link identities; a precomputed permutation maps thermal
link order to native patch order without repeated trial searches.

One native island is prepared over the actual representable interval
`(accepted_time + requested_dt) - accepted_time`. Every nonlinear trial
refreshes all original source claims at their respective interface temperatures,
solves all mechanical patches together from unchanged accepted momenta, and
routes their actual dissipation through `FrictionHeatNetworkStep`. That original
network sums powers into shared cells once, suppresses the replaced direct
contact links once, and retains all ordinary outside heat and other conduction.
No fixed heat fraction, duplicate capacity, prescribed face temperature or
numerical-loss heating is introduced. Geometry, circular area, centre and
normal must match the original native contacts.

Only a complete body/source/thermal candidate leaves the owner. Per-patch
mechanical and interface audits remain in force, and both direct
`Delta K + Delta H - Q_outside` and the sum of the original mechanical/network
residuals must pass independently. Final-candidate caloric entropy can be
requested from the same original network. Source failures, exhausted budgets,
late audits and cancellation preserve all accepted owners and exact retry.
Allocation, unknown count, nonlinear iterations and total trials are explicit
bounds; native controls are admitted before trial allocations.

Frozen geometry/normal data, the native response-rank and endpoint-event limits,
original orthogonal finite cells, common massless faces and backward Euler
remain the approximation scope. This does not consume retained elastic
partial-slip memory or rolling loss and does not establish root uniqueness or
experimental validity. Seven focused G1/G3/G4 groups in
`tests/articulated_friction_heat/hertz_multiple.rs` use an independent
three-body/two-patch mechanical response and eliminated finite-cell thermal
root, check source receipts and entropy, canonical permutations and signed
outside heat, stationary two-route conduction, source/geometry/ownership
refusals, loaded late failure/cancellation retry and large-clock synchronization.

### Complete vibrating history at moving-end Dirac attachments

The moving-end renderer consumes the original `PlateBankImage` through
`dirac_body::join_complete_body`. Mixed compact/Vk coordinates, full material
damping, Maxwell memory and attached acoustic histories remain in their
original bank order after the complete string state; any explicit fluid load
follows them, and the existing join multiplier remains last. The body force
port participates in the original common-velocity `join_port`; its independent
signed outward-volume port couples through the existing ratio-one transformer.
The surviving fluid source port is mapped explicitly after the original string
force ports. An unloaded body's pressure port stays open with zero applied
pressure and is not mistaken for a volume-flow source. A reed or volume source
without an independent duct inlet refuses instead of being silently ignored.

Joining accepted histories requires compatible initial attachment velocities.
Admission uses the original efforts and force columns with a relative
roundoff scale from their absolute products. Incompatible moving histories
refuse; no velocity projection, hidden impact, state reset or material-memory
reconstruction occurs. The multiplier has no storage and starts as a zero
representative. Exactly declared quadratic joins prepare the shared
`PreparedCommonFlowZoh` once, retain every physical coordinate, and return the
instantaneous endpoint reaction. Original material heat and source work are
integrated independently by that owner. Its complete-state and dense limits
remain; the renderer reuses candidate-state and mean-output buffers.

Nonlinear storage keeps the original descriptor/Gonzalez update and the
original complete pole-rate bound and midpoint admission, including retained
acoustic poles. Exact held-input evolution does not make state-dependent
exciters monolithic or supply a source-bandwidth/anti-aliasing guarantee.
The existing held-input, contact and observer conventions remain explicit.

The original bank observer reads its complete physical prefix at the declared
string offset, so compact and Vk modes cannot overlap or double count in the
pressure calculation. Focused native tests cover independently written
mechanical/fluid DAE evolution with 20/40/80-step refinement, work/heat and
constraint closure, nonlinear histories, malformed/incompatible joins, and
128 exact steps of accepted warm radiation memory at 48 kHz. Actual renderer
checks cover mixed and Vk-only bodies, signed areas, declared radiation at
48 kHz, render-length prefix identity, an exact quadratic 10 kHz common mode,
and retention of the nonlinear clock refusal. The native target uses the
actual exercised physics; unavailable bow/texture and cavity-source setup
are explicitly outside its tested scope. Full workspace and measured
real-time qualification are not asserted.
A dependent two-body/two-patch fixture also couples a shared native support
to two distinct original finite thermal pads carrying that same body's identity;
its independently derived full mechanical response and finite-cell root check
both body momenta, each physical heat contribution and all temperatures.

## Native circular wetting between unlike surfaces

`contact::capillary::AsymmetricCircularContactLineCapillaryModalSystem` binds
the original axisymmetric bridge with two independently evolving circular
footprints to the shared native mechanical transaction. Both ordered source
identities, surface histories, wetting thresholds and line resistances remain
with the original liquid owner. The actual body gap drives the joint footprint
solve; one line may pin while the other moves. Geometry and line perimeter are
three-dimensional, and each line exposes its own viscous and pinning losses.

Mechanical probes all start from the same accepted liquid state. A reserved
final replay evaluates both footprints on actual represented body motion,
with original component, connection, external-force and prescribed-gap work
checks. Geometric and quadrature uncertainty must fit the native compatibility,
force and energy gates before bodies, liquid history, diagnostics and clock
publish together. Refusal and cancellation leave both original footprint
histories and the native states retryable.

The owner inherits the original connected, quasistatic, acute-angle,
nonpositive-pressure, interior-neck and low-speed applicability limits. It
does not add bulk liquid inertia, evaporation, topology change or a full 3-D
stability claim. G1/G3/G4 native tests cover unlike-source pin/move dynamics,
individual friction sensitivity, surface exchange, reciprocal force and work,
original connections, numerical-accuracy refusal and exact retry.

## Constitutive graphs on the finite thermal timeline

`constitutive_graph::enthalpy::step_candidate` closes public material-law feedback on the existing `DepositionSimulation` and `BackwardEulerStep`. Each endpoint trial binds actual caloric temperatures and the original physical interval duration to declared graph inputs, solves the original dimensioned material equations using their consistent tangents, and deposits explicitly bound watt-valued graph outputs once into the receiving cells. An outer energy-coordinate Newton solve includes the converged material response in its numerical derivative. All trials borrow the same committed graph history; the returned joint candidate carries the thermal endpoint and proposed law history together. Fresh local and total joule gates and an independent external-supply evaluation precede return. Original material-source failures, events, cancellation, and aggregate numerical work limits preserve both committed histories.

`PreparedHeatNetwork::graph_step_candidate` composes this feedback with the existing material heat-network resolution, including source-integrated conductivity attachments, interfaces, and radiation. Its thermal-only admission refuses electrical circuit bindings whose losses would otherwise be omitted or counted twice. Closed-enclosure radiative roundoff is separately reported and excluded from the independent external-supply account. The adapter retains fixed mass, geometry, and topology within the interval, original caloric and latent storage, and first-order backward-Euler rates; it does not certify passivity of arbitrary authored graph laws, supply phase kinetics, or remove the need for independent time-step refinement.

## Continuing ambient-fluid changes on native rod trajectories

`RodMotion::rebind_submerged_fluid_candidate` updates an already bound material-fluid environment at the accepted clock without changing pose, linear momentum, centroid spin, mass, inertia, supports, or step count. The candidate reevaluates the native response at the actual moving pose. `RodFluidRefresh` reports both coherent transport snapshots and the change in hydrostatic potential relative to the original reference geometry; that change is external parameter work, with zero instantaneous impulse and no invented heat. It does not represent the thermodynamic cost of replacing a finite reservoir. Subsequent native kicks include gravity and buoyancy forces and their actual rotating attachment couples in the external linear and angular impulse accounts; conservative potential work is not counted again as applied-wrench work. Failed admission or cancellation returns no changed trajectory. Focused environment tests cover data-driven sinking, neutral buoyancy and floating; analytic piecewise acceleration across a density change; eccentric strands with different mass and displacement centroids; gas-pressure changes; stationary reactions; and rejected or cancelled refreshes.


## Simultaneous adhesive events on complete elastic networks

`contact::adhesive_events::linear_network::multiple` resolves up to eight
independent scalar junctions on the same complete mass-normalized network.
Every original mode, modal damping term, rank-one spring/damper connection and
rest extension remains in the shared matrix-free cell propagator. The scalar
owner uses that same physical propagator and resistance quadrature.

All sticking reactions solve one normalized Gram system; its orthogonal
projection constrains only the active contact coordinates. Sliding uses each
held kinetic capacity. Dense-polynomial roots locate the earliest stop or
release across every junction, including hidden intervening stops. Events
within the declared time resolution are applied together. Zero-velocity
continuation retains the largest feasible sticking set; equally sized
admissible alternatives refuse, rather than depending on contact ordering.
Dependent or unresolved contact columns are refused without regularization.

Per-junction slip, friction loss, impulse and impulse time moment remain
separate from original modal/link damping and external work. Work uncertainty
includes event timing, constraint projection, propagation and normalized-Gram
conditioning; it must fit the original energy gate. These are numerical
estimates, not interval event certificates. Active-set factors are cached,
workspaces are preallocated, and aggregate contact, cell, operator, root and
phase-trial budgets bound computation. Accepted state and phases belong to the
outer physical transaction; refusal and cancellation preserve identical retry.

The lane holds source capacities, local work columns, coefficients and external
loads during each step. It does not supply vector Coulomb cones, redundant
contact force distributions, changing normal geometry, nonlinear material
memory, or acoustic back-reaction. G1/G3/G4 tests exercise analytic nonorthogonal
stops and release, independent dense-ODE work, simultaneous/hidden events,
contact-unit and time-partition covariance, conditioning and ambiguity refusal,
limits, cancellation and deterministic retry.


## Simultaneous sourced adhesive events in the original network

`contact::adhesive_network_events::multiple::MultipleAdhesiveNetworkEventModalSystem`
retains up to eight independently sourced scalar junctions on one complete
original modal network. Each site supplies its own ordered normal/shear source,
held indentation and complete work column. Original DMT real area and junction
strength determine its static and kinetic capacities. Coupled event propagation
selects holding reactions, stops and releases jointly; contacts are not advanced
as separate pairwise mechanical steps. Per-contact phase, accumulated slip,
friction loss, impulse and time moment remain available alongside global work.

Single and multiple owners share one original native-state transaction. It
restores every original component under its budgets, independently evaluates
pressure and storage, retains every spring/damper/rest extension, and computes
each original connection's force and work port. All individual source, slip and
interval-speed limits and the original energy budget, including propagated
numerical work uncertainty, must pass before any accepted body, phase, frame or
clock changes. Late original-body or energy refusal and cancellation preserve
all junction histories for identical retry.

The fixed scalar work columns must be independent under the admitted normalized
Gram pivot. Dependent constraints are refused without regularization or silent
merging. Two scalar contacts are not a vector Coulomb disk. Normal geometry,
source conditions and external force are held during a step; evolving frames,
normal contact, nonlinear memories and radiation feedback remain separate
physical owners. G1/G3/G4 tests include simultaneous three-body stopping and
restart, constrained internal elastic motion, per-site order covariance,
zero-strength recovery of the original complete prepared network, reduction
to the single-junction owner, source/rank admission and joint rollback.


## Finite-rotation native rigid cohesive patches

`cohesive_contact::rigid::RigidCohesiveSystem` evolves two original
`fs_mbd::RigidBodyState` owners and their original full `MassProperties` through
a single finite cohesive patch. Every site retains its source-bound original
point, irreversible history, reference area, body attachment arms and material
frame fixed in body A. Ordered receiving material identities and actual initial
geometry are mandatory; construction does not infer or reset a bond.

A joint twelve-variable wrench solve replays both original native Cayley steps
and all original objective constitutive trials. Its fixed dense dimension does
not grow with quadrature; material work is linear in the site count per residual
evaluation, with preallocated site scratch, prepared LU and explicit evaluation
budgets. The actual midpoint arms transmit patch moments. The original frame
couple acts on body A, and spatial moments use the inverse arithmetic rotation
average, tau_body = R0^T (M - c cross M / 2). These conventions retain reciprocal
force, total angular impulse and the original finite-rotation work identity.

Actual native storage, original cohesive storage, intrinsic fracture loss, bond
viscosity, external force/couple work and gravity are checked independently.
Body states, every point, accepted observations and the clock publish jointly
after force/moment, energy, momentum, source, geometry and final cancellation
gates. Late refusals leave an identical retry possible. Failed bonds retain
fragment motion and damage; the original compression penalty provides
frictionless reclosure within its declared domain.

Source conditions, attachments and topology are fixed. This is a bounded local
branch solve over straight local constitutive paths, not a proof of global
softening-branch uniqueness or intrastep trajectory accuracy. Rotation-chart
and original material/rate limits remain mandatory. Bond creation, healing,
collision discovery and damaged-area friction are outside this owner.
`tests/rigid_cohesive.rs` checks an independent axial oscillator, global frame
covariance, finite-patch peeling with separate fracture/viscous losses, full
native free motion after failure, reclosure, gravity/external work, source and
geometry refusal, work-budget/cancellation refusal, strict energy refusal, and
late body-budget rollback with bit-identical retry.

## Molecular gas refresh on continuing DG ports

`DgGasRefresh::UniformMolecularDiffusion` selects the native uniform-background molecular/thermal entropy-and-compression transfer in the existing `prepare_dg_gas` transaction. It supports the original exact-mass and spatial wall routes, including the explicit spatial quadrature used by `DgSourceInput` gas fields; it never flattens spatial histories. The original molecular source law, species fractions, geometry and reservoir support remain bound. The same option is forwarded by existing one- and two-sided Bernoulli/DG gas-refresh APIs.

The complete port transaction rebuilds the native propagator and algebraic port relations while preserving independently sourced branches, terminal channels and every original exterior MIMO pole state. Existing equilibrium-pressure and total absolute parameter-work checks still apply. Interior material changes do not reidentify the exterior medium. Parameter work is charged once, separately from molecular/thermal dissipation; this operation does not itself solve background heating or species transport. A regression in `dg_port::gas::tests::molecular` exercises nonzero molecular histories with an independently sourced branch/channel and full mutual rational terminal. Its full-crate integration execution remains separate from the focused native DG constitutive/continuation tests.

### Source-material walls in molecular pressure feedback

`MolecularDrivingPoint::new_with_material_wall` retains one
`Arc<LocalMaterialWall>` and passes that same immutable source law into
`MolecularDuct::new_with_compliant_wall`. The material identity, thickness,
axial constraint, supported anisotropy and elastic/viscous/Maxwell response
remain with the existing wall owner. Feedback, terminal flow and all finite
chimney observations therefore consume the same jointly loaded gas operator.
`material_wall()` exposes the original binding for inspection.

The optional gas shear/thermal resolution remains independent of mechanical
wall displacement. The existing characteristic FIR ladder and pressure
reflection, finite-hole loads, opening fractions, wavefront factoring and
Rayleigh outlet observation are unchanged. The positive material surface
mass makes its high-frequency admittance decay as 1/omega, so the leading
wavefront remains the original frozen molecular gas wavefront.

Every actual radius and queried frequency must satisfy the material wall's
domain. The existing near-DC probe also requires material support; missing
data cannot fall back to a rigid wall. The caller still declares the exterior
gas independently. This is a single local membrane material over the admitted
duct and chimneys, without nonlocal shell coupling, exterior wall radiation,
nonlinear gas evolution or a new performance-file grammar.

The focused integration tests independently assemble a Debye caloric law and
a membrane-loaded cylinder, compare pressure reflection and physical outlet
flow, vary gas and solid relaxation independently in the receiver response,
and check common source-band refusals. These tests are retained for the
canonical coupling build; the session's executed native slice verifies the
underlying molecular/compliant-wall operator but does not execute these
coupling/FIR tests.

## Finite source-material cavity walls

`SharedCavity::new_with_cylindrical_material_wall` couples an explicit cylindrical gas region to a finite solid wall using the original `PreparedThermalTransport` source, evaluated at the reference gas temperature. Radial finite-volume capacities and Fourier conductances include gas, wall, and any original local molecular reservoirs in one caloric operator before pressure elimination. Its positive compression realization joins the existing exact mechanical runtime and harmonic response, with one accepted state, pressure, availability, dissipation, and rollback account. Gas and wall cell counts are independent spatial refinement axes; the cold realization is limited to 128 total thermal coordinates and rejects unresolved or unrepresentable spectra. The exterior wall is held at the reference temperature and endcaps are insulated. The scope is small perturbations, uniform acoustic pressure, fixed geometry, and frozen reference coefficients; this path does not infer mean heating, phase evolution, thermal expansion, interface resistance, exterior convection, or spatial accuracy from the declared frequency band.

### Complete vibrating bodies on the reed-driven ODE duct

The reed ODE renderer consumes the original `PlateBankImage`, retaining compact
and Vk modes together, full material damping and every accepted internal state.
The signed surface-volume port supplies reciprocal pressure reaction independently
of the external point-force port, which remains unforced. The complete original
plate observer uses the same state prefix. The existing inlet-pressure plus plate
observation convention remains; this change does not supply a new exterior outlet
microphone model.

For an exactly quadratic body/duct, a quasistatic reed uses the same complete
`PreparedExactZoh` held-volume-flow owner as the massive reed. The original strict
aperture root closes against the SAME interval's affine mean-pressure response,
with no delayed-pressure surrogate, hidden 5 Pa primer, extra initialization
step or render-length-dependent model selection. The original static opening
clamp and flow dead zone remain. Its physical flow is reconstructed from that
constitutive law; positive finite mean body impedance and the actual accepted
body pressure/work must agree before publication. Complete state/output buffers
are reused.

This connection gives exactly `body_work = dt * mean_pressure * held_flow`.
The massless valve's independently evaluated jet loss and mouth/body work obey
`jet_heat + body_work = mouth_work`; zero mechanical storage remains the declared
quasistatic model. The prior held-wave Robin connection included a within-window
pressure variance term absent from the valve's mean-flow work account. The
held-flow connection removes that coupling loss without modifying original
material/acoustic resistance or smoothing the aperture law. It does not certify
the physical adequacy of a massless model near a reed resonance.

For a massive reed and exactly quadratic body, the prepared complete body now
supplies its current-window held-volume-flow response `p_bar = free + B*q`.
`reed_bore::loaded::PreparedReedLoad` closes the original moving aperture, swept
flow and lay reaction against that response through the existing scalar contact
solver. The original linear reed mechanics evolve exactly under the interval's
held pressure/contact force; the aperture uses the midpoint gap. The body's
actual accepted mean pressure and work must agree with the condensed junction
before both candidate states are published. No previous-sample pressure drives
this quadratic path. Nonlinear bodies retain the original held-force Gonzalez
path and temporal admission; they are not inferred quadratic from probes.
Existing dense-state limits remain. Source histories seed a private one-shot
render snapshot.

The loaded reed retains original contact potential, exact integrated viscous
heat, nonadhesive contact loss, and Bernoulli jet loss as separate quantities.
It checks `deltaH_reed + viscous + contact + jet + body_work = mouth_work`.
The exact held-flow body's original storage, material/acoustic loss and port work
close the other side of that same interval. A finite strictly positive mean
body impedance is required; there is no impedance clipping, state projection
or substituted numerical damping. The source scenario's existing effective
area and lay calibration are retained, without implying measured specimen
authority. Exact linear propagation does not imply exact nonlinear continuous
motion or root uniqueness; the nonlinear aperture/contact step needs temporal
refinement.

Preparation owns the existing exact maps and contact storage. Warm reed trials
reuse the existing scalar solve and allocation-free original contact response;
exact blade velocity and viscous-heat contractions no longer allocate temporary
scaled states. A 1,000-candidate probe spanning open and contacting states
observes zero allocations. This is an allocation result, not a real-time timing
qualification.

Focused native tests execute the complete renderer with original gas,
viscothermal Foster, modal material, contact and numerical owners. They check both
reed selections with mixed modes and non-diagonally-dominant damping, independent
drive/surface maps, zero-blow silence and warm-memory response, render-length prefix
identity, same-window aperture pressure/flow against independent physical RK4,
and preservation of the nonlinear clock refusal. These are numerical/structural
checks, not measured-instrument, complete workspace or real-time qualification.

## Simultaneous Hertz surface and material rolling bulk heat islands — v0

`friction_heat::hertz_rolling::multiple::propose_rigid_hertz_rolling_island_heat`
joins one immutable prepared native five-coordinate contact island to the
original finite enthalpy network. Each original surface-friction source is
resolved at its own solved common-interface temperature and explicit pressure.
Each original constant-coefficient rolling material card is admitted at its
deforming sphere's actual bulk temperature. Interface identity must agree
between the original rolling and surface-friction owners.

Every patch routes its surface dissipation once through its distinct original
thermal link. Bulk losses from multiple patches on the same sphere are summed
into one original sphere cell, retaining one heat capacity. Different spheres
retain distinct bulk temperatures, inventories and heat destinations. Admission
requires the native whole-sphere mass, declared body/material identities, the
correct original route side and a caloric centre coincident with the native COM.
Conflicting mappings, a shared inventory for different spheres, duplicate
active sphere cells and duplicate thermal links refuse. Canonical permutations
between original patch IDs, thermal links and unique sphere bodies are
prepared once and retained throughout the nonlinear solve.

The original thermal network retains ordinary conduction and authored outside
heat. Signed outside heat is measured before injecting internal bulk rolling
loss, preserving its physical boundary and small signed contributions. Every
trial starts from unchanged accepted native momenta, original rolling histories
and caloric state over the same actual representable interval. The callback
supplies original endpoint conductivities and other heat only; interface work
and bulk work have distinct original owners and destinations.

Acceptance independently checks actual
`Delta H - Q_outside - sum(D_surface) - sum(D_bulk)`,
`Delta K + Delta H - Q_outside`, and the original native/network residual sum.
Only a complete mechanical/history/source/thermal candidate is returned.
Source failures, late audits, exhausted trials and cancellation leave accepted
owners unchanged and permit exact retry. Original cell/unknown, trial,
nonlinear, mechanical and physical-residual budgets remain explicit.

This retains frozen normal loads and poses, independent nonoverlapping local
sphere-deformation fields, the original rolling-card applicability box, one
uniform caloric sphere cell, and common-face backward-Euler transport.
Temperature admission does not create temperature-dependent viscosity.
It adds no spatial thermal validity, retained partial-slip memory, entropy
classification of injected bulk work, root-uniqueness or experimental claim.

Six focused G1/G3/G4 groups in
`tests/articulated_friction_heat/hertz_rolling_multiple.rs` use independent
coupled torque equations and a two-face thermal root, check one shared sphere
inventory, Hertz-induced rolling from rest, distinct original source
temperatures, signed outside heat and exact patch permutations. They also
check stationary conduction, zero viscosity, invalid inventories/interfaces,
late source failure/cancellation/budget retry, and two distinct deforming
spheres with nontrivial native/patch/link ordering and independent per-cell
heat destinations. All 23 existing groups in the original articulated-friction
heat test target continue to pass, for 29 groups in that target.

### Source gas in spatial thermoviscous interiors

`cavity_modal::interior::thermoviscous_interior` binds the original spatial
Helmholtz pressure solver and its viscous/thermal wall form to one coherent
`GasState`. Density and sound speed must match the FE volume medium bitwise;
viscosity, thermal conductivity, isobaric heat capacity and gamma are copied
from that same material evaluation. No dry-air constants or alternate gas law
are substituted. The original FE owner retains geometry, reduction, frequency
band and boundary-work semantics, including raw wall power up to roundoff.
This is a quiescent homogeneous frequency-domain interior with stationary
isothermal thin boundary layers. It does not add molecular relaxation,
finite-conductivity walls or narrow-gap linearized Navier–Stokes.

Six focused new native cases passed against actual FE/sparse/modal/material
source owners: four boundary-layer cases and two source-gas integration cases.
The complete Cargo workspace was not built in this environment.

### Spatial molecular gas changes in continuing port networks

`DgGasRefresh::SpatialMolecularDiffusion` routes the existing main/selected-branch
transaction through `DgDuct::rebind_spatial_molecular_gas_field`. The original
DG owner preserves local thermal/molecular entropy and weak compression/momentum
on unchanged spatial quadrature, while the complete port transaction retains
branch/rational-terminal histories and charges the aggregate absolute parameter
work. Existing common-pressure junction/exterior checks and all-or-nothing
publication remain in force. The background may change temperature and pressure
but must retain station positions and complete species fractions. This adds no
transport law, interpolation of histories or new exterior gas assumption.

Four new native physical-transfer tests and five existing uniform-transfer
regressions passed using actual DG/material owners. The added complete branch,
channel and mutual-terminal integration test is present but was not executed
because the full coupling dependency graph is unavailable in this checkout.

### Continuous molecular fields in feedback and opening sound

`MolecularDrivingPoint::from_gas_field` and its material-wall counterpart bind
the original primitive `AxialGasField` and molecular law to the same characteristic
pressure feedback, terminal flow and finite-chimney mouth observations. Axial
Magnus resolution, radial diffusion resolution and FIR length remain separate
axes. Inlet loading uses the exact inlet gas, exterior radiation uses the explicit
exterior gas, and terminal/side-mouth causal delays use the core's frozen-slowness
quadrature and each chimney's local gas. Uniform field arithmetic is preserved.
The source-resolved local material wall is shared with the propagation owner;
its positive surface mass preserves the gas-frozen wavefront.

The existing whole-field coverage, source band, common pressure, internal mutual
hole and work-cap refusals remain. The existing reed constructor accepts this
prepared driving point. The field is prescribed and frozen during playback;
this does not solve species transport or add nonlocal structural/exterior loading.
The two core wavefront regressions ran natively. Three added full feedback/FIR
consumer regressions are present but unexecuted because the complete coupling
and material-wall dependency graph is unavailable in this environment.

### Independent massive reed/body coupling checks

The actual full Foster-bore and massive-reed equations are compared with an
independent fine RK4 reference. Relative interval-mean pressure errors at
24/48/96 kHz are 2.215118e-4, 5.541245e-5 and 1.385526e-5, approximately fourfold
reduction per refinement for the smooth constant-pressure fixture. Fourfold
reed mass at unchanged static compliance changes actual bore pressure by 2.14463%
under that test's relative norm. Independent analytical damped-oscillator and
positive Simpson heat references cover contact entry, penetration and release,
with original contact potential, jet heat, mouth/body work and refusal replay.
Existing contact regressions remain active. These bounded numerical checks do
not establish measured-instrument fidelity or arbitrary-clock nonlinear accuracy.

### Molecular gas fields in authored reed performances

`ReedPerformance` accepts the existing `molecular`, `molecular_resolution` and `molecular_wall` records after an explicitly continuous `gas_field` or `gas_steady` block, or with a field injected through `from_bytes_with_field_resolution`. The parser retains the complete `RelaxationAbsorption` declaration and admits it against the actual field inlet. The continuous molecular owner then reevaluates the same law at every original field/geometry quadrature point, including species absent at the inlet but present downstream. The legacy homogeneous path retains its original evaluated-state constructor and arithmetic.

The exact field inlet supplies excitation and feedback; the exact outlet supplies the default baffled exterior, with explicit exterior records retaining precedence. The same molecular operator drives the ordinary reed voice, terminal/chimney volume flows and pressure observation. Authored-segment inlet samples remain read-only metadata; they do not replace continuous quadrature. Field/source/steady-solve receipts and observer-law selection retain their existing meanings. Injected fields retain their complete source owner while `gas_records` continues to mean authored records.

Molecular fields require `GasFieldResolution::Continuous`; midpoint input refuses instead of being silently reinterpreted. Explicit segment gases remain unsupported by this molecular grammar. Existing source domains, composition limits, common-pressure admission, finite work caps, internal mutual-hole refusals, FIR approximation and frozen-field scope remain in force. No new grammar, transport model, source fitting or audio effect is introduced. The new `molecular_input::field_tests` cover actual authored/injected field-to-pressure/block replay, law sensitivity, downstream source/domain refusals and the existing stationary thermal/species generator feeding the molecular performance. These full coupling/render tests require the complete workspace test lane.

All three new full coupling/render regressions are unexecuted in this environment.
Rust formatting/parser checks passed; the constituent molecular field and wavefront
operators have separate native evidence, without claiming a full render run.

### Original work through a quasistatic valve and a warm complete body

The first-window reference integrates the original body J/R/G under constant
physical volume flow independently of the prepared maps. A mixed compact/Vk body
with off-diagonal damping and warm Maxwell history is checked at zero and 100 Pa
mouth pressure by independently integrating its full states, pressure and original
R heat. Both body storage plus heat equals bore work and body storage plus heat
plus jet loss equals mouth work. At zero mouth pressure the original body motion
drives reverse valve flow and pays for both loss channels. Direct helper checks
also cover signed reverse flow and complete static closure. Existing massive
contact/work, refinement, state-prefix and nonlinear-refusal checks still pass.


## Joint two-tangent adhesive friction on complete linear networks

`adhesive_events::linear_network::vector::VectorLinearNetworkFlow` advances
the original mass-normalized modal coordinates and every original passive
spring/dashpot link with simultaneous Euclidean Coulomb disks. Each contact
has two independent tangent columns and one capacity; the monotone rung has
equal static and kinetic strengths. It has genuine set-valued sticking and
direction-changing vector slip, without per-axis force clips or smoothing.

The complete backward-Euler operator is condensed to at most sixteen contact
unknowns. A prepared mechanical factor and mobility are reused at a fixed cell
size; analytic disk projection derivatives drive a bounded semismooth Newton
solve, with a projected-gradient fallback. Workspace allocation occurs at
admission. Unit-column Gram pivots require independent, conditioned contact
coordinates; no dependent reaction choice is assigned physical authority.

Each cell measures original modal/link damping, contact impulse and its first
time moment, per-contact slip and friction work, and original stored energy.
The exact discrete loss `0.5*|delta_v|^2 + 0.5*delta_q^T*K*delta_q` is reported
separately as numerical dissipation, including each connection's spring share.
It is never material damping, friction heat or physical fracture loss. An
independent uncondensed equation check and Coulomb work gap contribute to the
numerical work estimate. Explicit cell, nonlinear, cancellation, physical-state,
energy-accuracy and cumulative numerical-loss gates precede return. Failed calls
retain no accepted history and retry deterministically from caller-owned state.

This is first-order discrete mechanics; cell refinement remains necessary.
Work estimates are not interval or trajectory-error certificates, and link-force
maxima cover discrete endpoints. No exact event time, unequal static/kinetic
strength, evolving normal geometry, moving tangent frame, material memory or
nonlinear network is supplied. Nine G1/G3/G4 tests cover analytic oblique stops,
disk-versus-box behavior, joint momentum/sticking, independent full-network RK4
and first-order convergence of state and numerical loss, basis/contact-order
covariance, internal elastic motion, partitioned impulses and moments, and
strict physical/numerical/cancellation refusal with bit-identical retry.


## Original native mechanics with sourced vector adhesive friction

`adhesive_network_events::vector::VectorAdhesiveNetworkModalSystem` binds every
two-tangent disk to its original ordered DMT/shear sources. Real contact area
and the original junction strength determine capacity. This owner explicitly
requires equal static and kinetic source strengths; it does not flatten a
velocity-weakening law into a monotone one. All native modes, modal damping,
Kelvin-Voigt links, rest extensions, pressure observations and physical limits
remain active. Each contact must supply both complete independent work columns.

The vector kernel and scalar event owners share original native state staging.
Actual component restoration, independent storage/work checks, each original
connection's force/work, source/speed/slip limits and propagated numerical work
uncertainty precede joint publication. The vector owner's backward-Euler loss
has its own `numerical_dissipation_j` field in both native and adhesive frames,
and a separately named numerical-loss entry in the native window audit. Physical
modal damping, link damping, friction, external work and numerical loss remain
distinct; `CoupledModalFrame::total_dissipation_j` contains physical losses only.
Each connection's discrete spring loss also enters its actual work transfer.
No numerical loss is assigned to heat, and no residual is erased as a loss.

Eight G1/G3/G4 tests cover diagonal loading outside the disk, joint sticking
with continuing internal elastic motion, tangent-chart rotation and contact
permutation, individual source causality, first-order refinement toward the
original complete prepared network, separate numerical/work-accuracy limits,
source and rank admission, cancellation and late original-body refusal with
bit-identical retry. The 26 existing scalar native event tests also pass.
This lane retains fixed normal indentation, fixed tangent geometry and source
conditions, with bounded internal timesteps. It does not add rotating frames,
normal impacts, velocity weakening, nonlinear constitutive memory, dependent
reaction selection or exact continuous event times.


## Simultaneous circular wetting bridges on original native mechanics

`capillary::multiple::MultipleCircularWettingModalSystem` binds up to 32 original
asymmetric circular contact-line bridges to one original modal network. Every
bridge retains its two independently sourced surfaces, physical volume,
footprints, pinning/sliding state and loss channels. Complete reciprocal gap
columns condense the full native displacement response, including all
cross-bridge terms, original component modes and bilateral connections.

A bounded simultaneous Newton solve uses the original analytic force tangents,
prepared dense LU and backtracking. Probes always start from accepted footprint
states over the same physical timestep. One aggregate source-evaluation ceiling
includes rejected probes and reserves the final actual-motion replay of every
bridge. All bridge/Newton/candidate buffers are allocated at admission; the
original native propagator retains its own allocation behavior. Inverse-Jacobian
gain screens ill-conditioned local continuation and propagates original force
uncertainty through every cross-bridge compliance. These remain numerical
estimates, not interval, global uniqueness or stability certificates.

One native stage applies the sum of solved forces. Every original bridge is
reevaluated at its actual staged gap, and body motion, all footprint histories,
gap coordinates, diagnostics and sample clock publish together after source,
force, compatibility, complete work/energy and cancellation checks. Original
surface storage, native loss, per-line viscosity/pinning and external/prescribed
work remain distinct. No heat recipient is invented.

Seven focused G1/G3/G4 tests exercise analytical three-body impulse/motion,
coupled two-bridge work, independently pinned footprints, six original elastic
modes with pressure and connection loss, reduction to the existing scalar owner,
bridge permutation, cross-bridge source causality, aggregate work/conditioning
limits, source identity, cancellation and strict/late-native atomic retry. The
18 existing circular and planar/axisymmetric native tests also pass. Each bridge
retains the original connected, acute-angle, nonpositive-pressure, quasistatic
liquid scope. This owner does not model liquid exchange/coalescence between
bridges, changing contact geometry, bridge rupture, liquid inertia, evaporation,
thermal deposition or exact within-sample trajectories.


## Damaged-area friction on original rotating rigid cohesive patches

`cohesive_contact::rigid::RigidCohesiveSystem::enable_damaged_area_friction`
attaches one original `ResidualCoulombSource` to each existing bonded site,
without replacing either native body, changing its full inertia, resetting
damage/history or advancing the physical clock. Original ordered surface cards,
query conditions and material assignments must agree. The explicit constitutive
selection is the shared damaged-area rigid Coulomb mixture: each tangent disk
has capacity `mu * area * integral(damage * compression_pressure)`, evaluated
analytically on the same original material path. Tensile cohesive traction,
surviving viscosity and endpoint damage are not substitutes for that measure.

Two local tangent-force unknowns per site join the original twelve body-wrench
unknowns. Every disk, body pose and material point is solved simultaneously.
Static capture uses the shared law's static coefficient, or kinetic capture
when already sliding; a yielded capture solution resolves the entire patch at
its kinetic coefficients. Numerical cone scaling introduces no elastic spring.
Prepared damped least-squares Newton corrections admit redundant sticking sites
without changing physical compliance or the final force equations. Individual
reactions may be a numerical representative of an indeterminate distribution;
neither their uniqueness nor global softening-branch stability is claimed.

Friction force uses the original arithmetic rotation-average pullback and its
balancing frame couple. Actual native linear/angular momentum and kinetic plus
gravity work therefore use the same finite-rotation geometry as each original
cohesive point. The shared constitutive audit independently checks the actual
local displacement, full two-dimensional Coulomb disk and dissipation; a second
spatial audit checks its objective work. Fracture, surviving bond viscosity and
friction losses remain separate, and friction is not automatically called heat.
The original surface compression energy is counted once. Complete failure
retains physical fragments and residual friction on reclosure, with no healing
or invented bond creation.

Site count bounds dimension at 140. Dense Jacobian, factorization and site
scratch are allocated at admission; no per-step friction buffer allocation is
required. Physical evaluations are linear in site count, while the joint
numerical differentiation and dense solve have their explicit bounded costs.
All probes, branch switches and final material audits share the original work
ceilings and cancellation gate. Bodies, all bond/friction states, observations
and clock publish only after the complete original audit; a refusal preserves
all accepted state and deterministic retry.

Eight focused G1/G3/G4 tests cover static holding between the static and kinetic
thresholds, an analytical two-tangent sliding trajectory, partial-damage area
scaling, independent quadrature through progressive fracture, distinct fracture/
viscous/friction losses, redundant quadrature partition, rotating spread patches
with full inertia, source/cancellation/late-energy rollback, and open/zero-friction
reduction to the original native dynamics. A resolved incoming slip cannot be
silently reflected by a zero-displacement midpoint capture: that case refuses
with a step-subdivision request. Exact stopping times, general intrastep reversals,
contact discovery, changing attachment topology, wear/heat deposition and
arbitrary moving fracture-surface material reselection remain outside this owner.

## Prescribed-input graphs and finite thermoelectric stores

A `MaterialGraphResidual` may have zero algebraic unknowns and zero equations when every graph input is prescribed. The original zero-dimensional Newton owner still performs source admission and the final fresh graph evaluation; both evaluations count against the work allowance, and the proposed graph history advances once from the same committed history. Empty residuals cannot suppress source, missing-input, budget or cancellation failures. The original graph-to-enthalpy adapter supports native Kelvin and Ohmic nodes on the same finite thermal cells: terminal Peltier and bulk Thomson heat are signed reversible transfers, Joule heat remains native nonnegative dissipation, and the original material heat network retains Fourier transport. Current-driven cooling and voltage-constrained generation share the same graph/enthalpy solver; there is no thermoelectric-specific time integrator or invented storage certificate for arbitrary graphs. Focused public integration tests compare both paths with the original circuit thermoelectric owner, an independent constant-property four-cell heat equation, and the native energy/entropy balances under material-only changes and continued thermal history.
