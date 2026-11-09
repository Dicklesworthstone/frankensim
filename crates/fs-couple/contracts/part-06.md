<!-- fs-couple contract, part 6 of 8; a higher number is newer. Index: ../CONTRACT.md. -->
## Complete material and acoustic memory in exponential impact

`ExponentialImpactSystem::with_relaxation_branches` and
`with_radiation_load` cold-compose the existing ImpactSystem material and
passive acoustic owners before the first accepted step. The complete assembled
pHS system is prepared through its native `share_storage` seam, so acceptance
and propagation retain the same original storage callbacks, stable nonlinear
remainder, structural matrices and physical resistance witness.

Each Maxwell arm keeps its original signed displacement projection,
`H=(sqrt(k)*projection.x-z)^2/2`, and `R_zz=1/tau`. Initial relaxed or
unrelaxed history is explicit. Radiation keeps the existing complete signed
mode mapping, passive pole model, original port interconnection and initial
rest state. Material and acoustic memory may be appended in either order;
their native positional metadata still observes the correct physical states.
No modes, off-diagonal terms, poles or histories are dropped, and no source
parameter is fitted or interpreted as another damping ratio.

The entire coupled quadratic block is propagated by the same existing exact
linear owner while nonlinear contact and WoolFelt remain in its implicit
discrete-gradient solve. The original positive loss calculation, conditioning
acceptance, source-domain checks, rollback and accepted clock remain shared.
Every material/acoustic coordinate counts against the unchanged complete
exponential capacity, and its combined instantaneous spectrum must pass the
original sampled-band gate. Original radiation coupling/clock limits still
apply. Empty branches and admitted zero-coupling loads preserve the prepared
image. This is cold composition; live material replacement and its parameter
work are not supplied.

Four focused regression groups specify a manufactured Maxwell eigenmode with
independent positive heat, full cross-mode Maxwell/acoustic response and work
against an independent ODE in both attachment orders, joint nonlinear felt
conditioning with refusal/retry, and complete-state/source/clock/lifecycle
admission. Pinned Rust parsing and formatting passed. Full fs-couple execution
of these four groups remains outstanding; no numerical or replacement-host
substitute is claimed as runtime verification.

Shell/support body potentials, undeclared or nonquadratic strain maps, nonlinear
contact drag and thermal feedback remain outside this entry. A passive
authored acoustic pole model still requires separate geometry validation.
Exact linear propagation does not remove nonlinear contact, source/model,
modal-truncation or temporal-refinement obligations.

### Automatic physical power and cycle-energy adjoints

`PhsPowerObjective` selects mean supplied/resistance power or supplied/resistance
energy per actual cycle. The forced and autonomous sensitivity owners expose
`prepare_power_adjoint` and `prepare_power_adjoint_analytic`: they sample the
original value and complete coefficient (plus autonomous frequency) covector,
then reuse the existing transpose factor. Source-bound
`PreparedPhsForcedPowerAdjoint` and `PreparedPhsAutonomousPowerAdjoint` expose
`gradient_parameter`, accepting the existing complete physical direction and
automatically composing its direct objective and residual chains. These four
objectives require no manually assembled explicit objective partial.

Forced directions retain complete J/R/G/storage/Fourier-input changes and the
prescribed clock. Per-cycle work includes `T*(P_p-Omega_p*P/Omega)` as its direct
partial. Autonomous directions retain the original `G_p/H_p -> y_p -> u_p`
feedback chain once; cycle `q_Omega=-q/Omega` enters the unchanged phase border,
while the explicit parameter partial is `T*P_p` at fixed unknown frequency.
Both grids contract their own direct mean with the cached dual and independently
check direct-partial and total-gradient sampling. Supplied work and resistance
loss retain their independent original laws and signed physical contractions.

Original storage/source/feedback callbacks still evaluate cached physical
states. Parameter queries perform no objective Hessian/feedback-Jacobian
reevaluation, new orbit/factor, or waveform-direction solve. One private source
sample evaluator shares the original direction assembly and scratch between
residual and objective callbacks. Positive objective, parameter and HB scales
condition acceptance only; selected power/energy units apply to the complete
sampled value and coefficient/frequency covector. Refusal and cancellation
publish no partial result and preserve the original prepared source.

The 23 focused physical tests include three new autonomous/source-owner groups
and three expanded forced-direction groups. They cover independent complex
mechanics and separately resolved nonlinear memory systems for all four
objectives; all twelve physical parameter changes of an analytic autonomous
fixture; separately evaluated supply/loss laws; a nonzero period contribution;
conditioning-scale invariance; cached tangent reuse; cancellation/refusal;
signed raw R; and weak net work amid opposing `2^54` terms. Conditioning-scale
variation is not a physical coordinate-unit conversion. These fixed-order
gradients do not establish harmonic-tail convergence, stability, passivity or
a heat destination.

### Complete forced structural/interior resonance parameter jets

`StructuralInterior::prepare` now exposes
`PreparedStructuralInterior::solve_with_frequency_parameter_derivative`
through the existing resonance module. One original coupled factor supplies
all four derivative orders for structural displacement/velocity, pressure
potential/pressure, reciprocal fluid force and displaced flow. Original
structural/acoustic M/K, viscous/gyroscopic C, wall damping and both opposite
interface-direction blocks are supplied through the same physical direction
assembly as the pole owner. Complete complex force/flow source jets retain
frequency, curvature, parameter and mixed phase/amplitude laws.

Pressure remains the velocity of pressure potential: every derivative of
p=-i*omega*psi is retained. Interface recovery includes both B and B_theta,
including the mixed B_theta times frequency-tangent terms. Separate structural
input, acoustic input, signed fluid-to-structure exchange, structural heat and
wall heat expose all five power orders. Symmetric damping and its direction
supply heat; gyroscopic and reciprocal skew work never becomes dissipation.
Borrowed field jets feed the shared coherent receiver, peak/notch and linewidth
sensitivity queries without a new resonance solver.

The adapter reserves 24*n*n+1536*n+2048 additional bytes and work score
2048*n*n+4096*n before physical assembly/recovery. Coordinates and topology
stay fixed; the caller supplies constitutive/geometry directions. This entry
covers the memoryless complete interior pencil; additional exterior loads
refuse unless their complete physical law is supplied through another owner.
No finite parameter interval, continuum convergence or forward error near a
pole is certified. Focused native tests execute the actual recovery and modal
owners against independent complete equations, changing interfaces, DC,
dilation and five-order balances. The original tetrahedral FE regression is
also defined in `resonance_tests.rs`; the bounded native verification does not
execute its full geometry/assembly dependency closure.

### Source-bound shell radiation-load curvature

`PreparedAssemblyMimoHarmonicResponse::solve_with_frequency_curvature`
binds the original admitted continuous MIMO source to the loaded shell's
complete physical curvature path. Source-owned analytic Z_s and Z_ss become
Z_omega=-i*Z_s and Z_omegaomega=-Z_ss on the structural exp(-i*omega*t)
axis. Source phase/curvature, real and conjugate poles and all mutual entries
remain; no impedance is frozen, refitted or interpolated. The original source
band, port ordering, pole history, clock and admission reports are unchanged.

Three source matrices require 48*p*p additional bytes and work score
1024*p*p*(expanded_order+1), reserved before evaluation and deducted from
the native curvature allowance. Physical support/tie/port recovery retains
its own complete budget. The source's independent validity/passivity class
does not gain broadband authority from a local derivative.
Focused native tests combine the actual admitted source and modal owners with
an independently cleared Maxwell memory equation, its asymmetric physical
power curve and maximum at sqrt(2), coherent source-phase invariance,
history/report retention and band/resource/cancellation refusals.

## Original support mechanics behind anisotropic solid coatings

`elastic_reflection::NormalOscillatorBackedElasticCoating` binds original
incident and finite-coating `NormalElasticCharacteristic` owners to the same
three principal `WallPin` laws and physical support rotation. Both original
full tensors retain their native quasi-modes, material frames and strain
domains. The material owner resolves normal-wave propagation, coherent
conversion, interference and interior fields; the adapter supplies the same
original support tensor at the same carrier frequency.

Only actual rear motion feeds the shared principal-coordinate displacement
and kinetic/spring storage evaluator. Incident-face motion is not used for
support storage. Original positive backing work remains the sole dissipation
in these two lossless solids; coating elastic energy is not counted as pin
storage. Exactly zero thickness recovers direct anisotropic supported
reflection. Modal and Cartesian drives reuse prepared maps without allocation,
factorization or constitutive source queries.

This is stationary harmonic reflection with normal wavevectors on homogeneous
planes. Original anisotropic energy can walk laterally; finite-aperture
transport, oblique anisotropic conversion, attenuation and transient layer
history remain outside this owner. Five focused G1/G3/G4 regression definitions
cover exact absence, an independent Cartesian continuum-transfer reference
for two noncommuting tensors and an oriented support, parity with the original
isotropic coating at normal incidence, reactive storage without heat, and
original coating-domain/cancellation admission. Native execution status is
reported separately.

## Initialized source storage in the live nonlinear valve and both acoustic ducts

`TwoSidedAperture::with_stored_pressure_source` binds the original
`StoredPressureSource` to one explicitly authored Driven endpoint at the remote
end of the upstream duct. Its wave reference and physical sample period must
match the actual adjacent section and valve clock. A cold attachment retains
initialized compliance pressure/inertive flow and admits their original
quadratic storage against the existing complete energy bound. It does not
overwrite a passive, fitted or already owned termination.

`step_with_stored_source_under_gate` and the corresponding downstream-boundary/
spatial-pressure entry point solve the source, original nonlinear valve and
both retained networks in one accepted transaction. Every endpoint reads the
old traveling state. The remote source uses the native causal R-L-C reflection;
the valve still uses its original combined characteristic load, Bernoulli and
swept-flow solve. There is no extra source/valve solve, wave-history clone,
replacement FIR termination or implicit zero-delay shortcut across the duct.

The remote ideal pressure supply works against its own outgoing flow, which
can differ from valve flow while the duct stores energy. This work and the
separate local series-source work enter once. Initial source energy is storage,
not supplied work. Original source resistor heat, source storage change,
mechanical/contact/material dissipation, both graph accounts and any fixed
supply-passage account participate in the original whole-system gates. Internal
source-to-duct transfer cancels once. The original paired graph acceptance
remains the last fallible physical operation, followed by infallible publication
of valve/source/optional passage histories. Refusal and cancellation preserve
all accepted clocks, source state, waves, motion and observations.

The upstream graph may contain passive branches and its admitted distributed
loads, but must have exactly this one Driven endpoint and no moving-key owner.
Downstream keys and original spatial pressure ports retain their explicit
entry point. Existing source-free/local-pressure entry points and arithmetic
are unchanged; a bound remote source requires an explicit remote pressure
input. Coefficients and geometry remain fixed. This is a declared linear R-L-C
source, not finite thermodynamic gas inventory, anatomical tissue properties,
moving lung volume or measured player validation.

Six integration regression groups cover initialized pressure traveling to
independently moving nonlinear blades, bore reflections returning to the
reservoir only after the full round trip, independent local/remote source work
with initialized R-L-C storage and resistor heat, the zero-storage
pressure-release limit, late whole-system energy rejection/cancellation with
exact retry, and source/clock/termination/storage admission. All six passed in
a focused native build executing the actual source and transaction methods,
original uniform-slit paired midpoint solver and native fs-vfit wave/load
owners. Its open-contact/scalar-gas fixture explicitly excluded closed contact,
exact/modal/Maxwell blades, moving keys, source geometry/material lowering,
supply passages and coupled radiation; those combinations were not runtime
verified. Seven existing stored-source regressions also passed after sharing
its immutable boundary preview. New files passed pinned formatting; full
workspace Cargo/remote-runner validation remains unavailable.

### Original oscillator-backed coatings on the native thermal timeline

`OscillatorBackedViscoelasticLayers::prepare_heating` binds the original
supported coating and authored peak pressure to native
`SupportedViscoelasticHeating`. Exact original bulk and shear work enters the
existing physical depth cells; accepted thermal intervals and explicit power
envelopes own each joule. The adapter adds no observer, heat account or time
integrator.

A distinct optional backing cell receives the original resistance work once;
otherwise that work is reported as exported. Its caloric mass remains an
explicit thermal property and is not inferred from the pin's reduced inertia.
The wrapper retains the original `WallPin` and exposes `peak_response` for the
actual rear displacement and cycle-mean kinetic and spring storage. Those
stored mechanical energies never enter thermal source credits. The peak
response remains the authored carrier response after an envelope changes or
ends; it does not represent an instantaneous transient state.

Original coating and support coefficients remain frozen inside their declared
temperature bands. The native thermal engine retains conduction, boundary
events, material storage, accepted-step decisions and atomic refusal. The
adapter does not infer temperature-dependent mechanics or oscillator ring-up.

Integration fixtures cover depth-resolved Maxwell-layer heat, a distinct
backing receiver with independent caloric mass, reactive storage without heat,
affine envelopes and dark intervals, split versus whole advances with native
boundary events, original source bindings, temperature support, budgets and
late cancellation without mutation of the input histories.

## Internal finite ducts between modal acoustic compartments

`vibroacoustic::CavityCompartmentBasis` combines 2..=64 original sealed
compartment bases. Every uniform compression mode is retained first, in
compartment order; higher modes follow a deterministic frequency order.
Original pressure normalization and modal volume norms are unchanged.
Interface samples concatenate, with each pressure trace supported only on
its own compartment. A caller-specified `maximum_trace_entries` bounds the
complete zero-padded trace payload before allocation; reservations are
fallible. The controlled constructor polls between compartments and output
modes and can consume the existing Cx cancellation predicate.

`InternalCavityDucts` connects declared pressure faces of different
compartments with finite native cylindrical TwoPorts. Both end flows remain
unknowns in the SAME structural/cavity/exterior solve. Positive left flow
leaves the left compartment; positive right flow enters the right compartment.
The two flows need not agree because a finite duct retains compressibility.
Native scattering constraints introduce no transmission, pressure-node,
duct-admittance or cavity-resonance division. The original compact and
distributed exterior opening owners receive the full simultaneous matrix
stride and retain their complete mutual wall/mouth and mouth/mouth radiation.
A sealed exterior and a structural-only radiation attachment also remain
available.

Each duct face retains the full original local pressure projection, including
its unit-volume-flow uniform-mode normalization. Optional unresolved entrance
inertia is admitted independently at each end only for that compartment's
single uniform pressure mode and reference radius ka <= 0.25. A resolved
entrance uses no additional mass. Native right entrance composition is on the
right side of the original cylinder. Cavity norms/volumes exclude internal
duct volume; the exterior does not own an internal connection's reaction.

The graph reuses the SAME cylinder preparation and face-work observation
helpers as distributed exterior openings. Equilibrium propagation retains
native Bessel viscous/thermal response. Molecular propagation retains the
original MolecularDuct and complete RelaxationState, with the selected native
wall spectrum or explicit adiabatic/slip choice. All compartments and ducts
share the original gas and molecular selection. Both attachment orders refuse
equilibrium substitution or a changed molecular law; exact source-frequency
selection is preserved. Existing common scalar cavity loss and molecular-loss
exclusion remain in force. Original geometry, actual pressure-trace identity,
outward normals and structural-coupling identity remain caller-owned:
matching spectra and norms alone cannot establish that a different mode shape
is the same original basis.

Responses expose each end pressure/flow, native TwoPortPower observations,
signed work leaving each physical compartment, and internal duct absorption.
Independent two-ended pressure work retains transfer and cancellation.
Internal absorption is debited ONCE as `PowerBreakdown::internal_ducts`;
`cavity_connections` is the signed internal transfer used in the cavity
interface audit. The total balance includes structural loss, cavity loss,
external opening-duct absorption, internal duct absorption and complete joint
radiation. No recovered transfer, reactive storage or numerical work defect
is converted into an extra heat source.

Truncated solves must keep every physical compartment's uniform mode.
Convergence keeps that minimum pressure basis and compares both internal end
flow vectors independently of the existing modal and exterior-opening
metrics. This observes the actual graph under reduction instead of replacing
an omitted compartment with an unintended fixed-pressure reservoir.

Scope is frozen, linear, quiescent plane-pressure propagation through straight
mechanically rigid ducts, with original modal cavity and structural bases.
Different compartment gas states, mean flow, nonlinear aperture exchange,
thermal-inventory evolution, compliant duct walls and transverse duct modes
are not inferred. Frequency-domain absorption is not a transient caloric
update or an independently resolved local heat-density field.

Verification: 12 native physical tests pass with the pinned July 2026
compiler. The existing eight opening tests remain green. Four new tests
cover independent finite-duct continuum transfer between sealed compartments
at quarter/half/full waves with and without attenuation; simultaneous internal
ducts with BOTH production compact/distributed exterior dispatchers and
complete mutual work; duct-orientation and signed pressure-basis covariance;
opaque internal transmission; and original combined pressure support/norms,
trace-payload boundary and cancellation. The native source slice executes the
actual basis, numeric operators, production opening dispatcher, original
TwoPort wave/composition/power methods, fs-math and fs-la complex LU. Narrow
error/data-carrier scaffolding does not replace any tested numerical owner.

Three public integration tests specify native Bessel/Rayleigh composition with
both exterior routes, exact molecular cavity/duct binding in both attachment
orders, original complex-compliance transfer, uniform-mode truncation, and
physical/source admission. Those public tests were syntax parsed and source/API
reviewed but were NOT Cargo-executed. Missing sibling workspace sources and
writable dependency/build storage still prevent a complete Cargo/DSR gate;
the focused native result is not a whole-workspace build claim.

## Original spatial Maxwell memory with finite caloric storage

`rigid_damping_heat::relaxation` connects the native simultaneous spatial
Maxwell graph to the existing rigid/thermal transaction. The original rigid
and Kelvin graph, complete branch factors/times, initial elastic histories and
held efforts are explicit inputs. No relaxed spring geometry, equilibrium
stiffness, source calibration or temperature dependence of those histories is
inferred. Empty branch banks retain the existing driven viscous behavior.

Every native moving source integrates original Kelvin and Maxwell heat along
the same trajectory. Its complete candidate memory stays in the pending
mechanical receipt until the matching caloric interval and observer accept.
Only then does the private mechanical acceptance hook move that bank into the
continuing history. Final source/thermal/energy refusal or cancellation still
discards the entire candidate; borrowed original histories never change.

The shared report adds `memory_energy_change_j`, zero for memoryless paths.
The final independent energy account includes the original-to-final Maxwell
storage change, retained before subtracting large absolute energies. Elastic
recovery may accelerate the original bodies while its positive dashpot loss
heats original finite cells. Neither released elastic storage, signed actuator
recovery nor numerical drift is reclassified as heat. Shared receivers,
partial routing, unallocated outflow, nonlinear calorics and actual native
event clocks retain their existing semantics.

The complete body and memory output banks are moved from the native receipt;
initial memory is borrowed until an interval accepts. Branch, source, receiver,
native work and dense-memory limits remain explicit. Original factors/times
and instantaneous viscosity stay fixed through this entry. Native geometric
composition and caloric transport require independent timestep refinement.

Four focused integration groups cover preloaded memory recovery with finite
nonlinear calorics, driven continuation against an independent physical ODE,
empty-memory parity with the existing driven graph, and joint history/observer/
caloric cancellation and refusal. All four passed against the original native
MBD/pHS and finite-caloric/thermal event owners in a bounded native build, as
did all four existing driven-thermal regressions after the shared acceptance
hook changed. Numerical kernels and tested source were unchanged. A complete
workspace build remains outstanding.

## Continuous force-driven tensor solids with finite-cell self-heating

`tensor_material::solid::thermal::tones` connects
`PreparedTensorSolidBody::advance_tones_thermal_candidate` to the existing
continuous-tone mechanical flow and original finite-cell enthalpy timeline.
All supplied nodal tones act simultaneously on the original consistent inertia
and retained Kelvin/Maxwell system. Each mechanical interval advances once,
independently of caloric event/CFL cuts. Positive original element heat enters
its explicit whole-reference-mass/volume receiver exactly once.

The accepted history must borrow the exact prepared body, share the thermal
clock and retain one complete unique caloric route per element. Exact caloric
law identity, reference inventory and explicit fully-solid frozen-coefficient
temperature intervals remain binding at thermal substeps/events. All native
observer source, transport, interval, transition, event and timestep hooks
survive through the shared heat handoff. Mechanical history, thermal state and
observer history publish together; late refusal and cancellation preserve both
accepted endpoints.

The report retains original nodal force-work moments, kinetic/elastic/Maxwell
storage and independently integrated physical heat. Stable original mechanical,
native caloric and actual heat-handoff residuals form the joint account;
recovered storage and numerical defects never become heat. Mechanical
coefficients, supports and geometry remain fixed; the original full/reduced
basis scope and endpoint-only strain/startup-band limitations remain. Mean heat
over a fixed mechanical interval preserves the original explicit thermal
order/CFL; mechanical duration and thermal cuts need physical time refinement.

Three focused tests specify independent original Tet4 inertia/strain/Maxwell
ODE work and heat, startup interference and unforced decay, original multi-cell
routing with transport/events, and source/clock/late refusal/retry/cancellation
isolation. Native execution and shared-helper type checking were unavailable
when this extension was authored; no runtime validation is inferred.

## Original relaxing exterior in the loaded molecular network

`MolecularRadiatingNetwork::prepare_fields_with_exterior_material` uses an
explicit `MolecularExteriorMaterial` with the original borrowed RelaxationState,
declared angular-frequency band, native BEM formulation and native assembly/
response limits. The complete GasState must match the existing exterior record.
A sound-speed/density pair or equal source labels cannot stand in for that law.
The original regional interior gases, material-wall Arc owners, source clock,
generator impedances, probes, geometry and opening reference planes are retained.

The existing `MolecularBemObservation::prepare_response` computes the full
self/mutual opening impedance and finite receiver p/u bank together. Unit-integral
1/area normal-velocity traces retain pressure/volume-flow coordinates. Native
complex K(omega), the decaying Green function and q = i*rho*omega*v remain in
one original solver; k times equilibrium sound speed never replaces omega in
the velocity trace. At most 32 network openings fit one original 256-field BEM
batch. There is one native assembly/factorization for this exterior basis and
one original loaded molecular network factorization for all source coordinates.

The paired interior/exterior result retains the original native response and
an exact bounded clone of the already solved port bank. Coherent drives,
independent source statistics, directional microphones, source control and
recorded Welch bins reuse its matched harmonic observer. No second boundary
solve, additional attenuation filter, material query or opening-power
independence assumption is introduced. The existing field allowance admits
three additional opening-work matrices and their original panel projection;
ComplexOptions and MolecularResponseLimits independently bound native work.

`PreparedRadiatingNetwork::exterior_material` exposes the exact borrowed source,
band and formulation; `molecular_exterior` exposes the same native response.
Its source inputs are actual outward opening flows in boundary-label order,
not generator amplitudes. Existing `exterior_surface_power_w` is signed opening
boundary work already included in the graph external-load account. In a relaxing
exterior it includes molecular absorption and outgoing energy; it is neither
a second dissipation debit nor power on an arbitrarily distant sphere. No
automatic thermal receiver is inferred. Lossless `far_field` stays unavailable,
and the native harmonic owner's radiation-momentum refusal remains intact.

Stationary Welch-band composition additionally requires the exact same original
exterior source object, declared band and boundary formulation. It refuses
ordinary/relaxing mixtures or numerically identical but distinct source objects,
including at silent bins. Domain, storage, work, geometry, loading or cancellation
refusal publishes no partial paired result and leaves the original owners usable.
The default nondispersive exterior path keeps its existing owners and semantics.

This remains a held homogeneous molecular-compression exterior with inviscid
momentum, rigid panels outside the declared openings and original spatial
quadrature. Exterior shear/thermal diffusion, mean flow, distributed wall
radiation, startup and evolving materials are not added. Three focused native
regressions in `driving_point/network/radiation/tests/prepared/molecular/fields/exterior.rs`
cover independent combined-velocity BEM loading/fields and material-wall work,
empty-reservoir and measured-coherence limits, and gas/domain/budget/cancellation
and stationary-source identity. Synthetic source fixtures are numerical tests,
not empirical validation. Whitespace checks reported no errors. Native tests,
Rust compilation and rustfmt parsing were not executed because the original
compiler/toolchain was unavailable; source review is not runtime proof.

## Molecular duct fields and sound under the same full exterior load

`RadiatingNetwork::bind_molecular_with_compliant_walls` binds original molecular
edge laws and exact edge-id `MolecularCompliantWall` Arc owners to the retained
regional network and exterior. Complete material state, source identities,
Maxwell response, thickness/restraint and radius/frequency support remain in
the original wall law, including the existing `LocalMaterialWall` implementation.
Native media pins remain live on unlisted edges; unknown/duplicate ids and
overlapping pin/custom declarations refuse through the same duct owner.

`MolecularRadiatingNetwork::prepare_fields` accepts a
`MolecularRadiatingFieldRequest` with the physical frequency, source clock,
complete scalar/coupled generators, non-opening scalar loads, original probes
and `NetworkFieldBudget`. It computes the original full BEM opening basis
once, then supplies that exact mutual radiation matrix to the native molecular
field preparation and its single loaded network factorization. The existing
exterior observation code consumes an exact clone of that solved native port
bank; it does not independently reconstruct the interior or choose new sources.
`prepare_fields_with_radiation` derives its optional angular observation from
the SAME unit-opening BEM solutions.

The returned `PreparedMolecularRadiatingFields` exposes `interior()` and
`radiating()`. Internal pressure, directed volume flow, coherent moments and
source-owned cylindrical losses share source coordinates and phases with the
actual opening flows, directional microphone fields and far-field response.
Source-amplitude and covariance-factor queries are contractions only. No new
material evaluation, radial eigensolve, duct factorization or BEM solve occurs
during these observations.

The original exterior field allowance additionally admits four
port-by-source complex arrays and one source vector for the retained port-bank
clone before the BEM solve. Native interior probe/transfer/chain caps remain
separately binding, as do exterior geometry, angular work, loaded conditioning
and cancellation. Refusal returns no partial paired preparation and leaves the
original molecular, material and exterior owners usable.

Mechanical compliance changes the loaded opening flows. The exterior remains
rigid outside its original opening patches; distributed wall radiation needs
its own source/trace owner. Gas absorption, mechanical wall work and actual
opening/exterior work stay separate. Only the native viscous channel is
declared direct gas caloric heat; thermal and molecular availability and wall
work do not specify an automatic thermal receiver. Far-field power observes
the existing exterior load work and is never added again as an energy debit.

Three native regressions in
`driving_point/network/radiation/tests/prepared/molecular/fields.rs` use the
original resolved material and Maxwell membrane APIs. They compare independent
complex line loading, native interior gas/wall work and a separately assembled
physical BEM velocity boundary; coherent loss, microphone and angular source
factors; and original material support, retained-bank/field budgets,
cancellation and retry. The source claims are synthetic implementation
fixtures, not empirical material validation. These tests and Rust
parser/formatter checks were not executed because the workspace was
unavailable; the production source received a bounded independent review.

## Finite thermal bodies in moving molecular networks

- `molecular_exchange::network::thermal::{step_candidate, step_candidate_with_convection}` couples finite enthalpy and phase-change bodies to the original moving-gas molecular network. Bodies sharing a gas use the existing simultaneous implicit thermal solve. Independent gas groups start from the same snapshot; the original molecular transaction then samples the complete heated snapshot. Both operators cover the same physical duration, with first-order Lie splitting.
- Each gas has at most one thermal group, and physical body IDs are globally unique. Bodies retain their original material curves, fixed mass and volume. Constant-density curves and consistent mass/density/volume are required so enthalpy change represents fixed-volume stored energy. Heating preserves gas mass, volume, composition, auxiliary conditions and momentum; actual represented momentum and kinetic changes are checked in both the common and laboratory frames.
- Original aperture sources and drift bounds are admitted before heating. Native gross withdrawal, species, momentum, caloric and endpoint support checks remain authoritative. State-dependent convection receives the actual implicit gas state and current velocity. The same repeatable law must also admit the final gas state, post-exchange velocity and accepted wall state. For this final admission only, a copied body descriptor carries its accepted material state; its original curve, geometry, sources and duration remain unchanged. Final receipts apply no extra heat.
- Original source-relative gas caloric integrals and actual body storage close the heat stage. Authored gas/body heat and enclosing radiation enter once; body/gas heat and reciprocal molecular streams cancel internally. Native composition-aware accounts and independently reconstructed thermal corrections close both energy frames. One fixed noncancelling energy allowance and per-component momentum allowance bound the complete proposal. Separate thermal-stage debits remain available for adaptive accumulation.
- Gas, body, heated-snapshot and source receipts publish atomically after the final cancellation checkpoint. Caller work limits and native iteration limits bound evaluation. Supplied heat-transfer laws must qualify the actual contact/flow regime. This fixed-volume, homogeneous-reservoir model does not infer convection coefficients, body expansion/work, pressure or drag forces, container mechanics, or finite-channel flight laws.
- Focused checks cover an independent simultaneous multi-body solution, native heated-network replay, opposing wall transfers, changing composition and large formation datums, common boosts, final source admission, radiation, refusal/cancellation, work limits, and empty or insulated native parity.

## Adaptive moving molecular networks with finite thermal bodies

`molecular_exchange::network::thermal::trajectory::{advance_candidate, advance_candidate_with_convection}` evolves the original moving gas network and finite material bodies over one physical window. Every full and half trial calls the native thermal/network transaction; only accepted half steps supply the next state, chronological history and integrated accounts. The original simultaneous shared-gas thermal solve and simultaneous molecular exchange remain the physical owners, with refined first-order operator splitting.

Body input durations specify the complete window. Authored gas heat and body internal heat represent held powers, apportioned to trials with compensated arithmetic. The final account uses the authored whole-window request and actual accepted native radiation. Source partition error and absolute proration roundoff remain explicit; reciprocal body/gas heat is internal. Native fixed-volume material, source, implicit-solve and final convection admission remain binding.

Local estimates include gas temperature, pressure, every admitted species inventory, velocity and retained momentum increments. Each body independently controls temperature, actual stored-energy increment and liquid fraction, including isothermal latent plateaus. Signed, positive and negative convection/radiation transfers are compared separately per body, so opposing contacts or temporal reversals cannot hide transfer error. Each molecular link retains both gross beams and its authoritative signed net, with vector momentum and source-relative energy evaluated from the exact heated snapshot sampled by each trial.

Existing molecular trajectory helpers own source references, state and stream estimates, accepted molecular totals and original-to-final gas storage/kinetic energy. Opt-in thermal reference registration includes isolated heated gases while preserving the molecular-only path. Gas IDs, body IDs and authored source reductions are canonical. Native thermal energy and momentum debits are added once beside native molecular debits; they are never reconstructed by subtracting rounded combined accounts.

One fixed window allowance bounds actual species/total mass closure, momentum components, common-frame energy and laboratory-frame energy. The energy audit independently combines original-to-final gas and body storage, gas kinetic change, authored heat and accepted radiation. Accumulated absolute defects and source reserves cannot cancel or expand into a per-step allowance. Reports retain native transactions, integrated molecular streams, body heat/storage totals, residuals and consumed budgets.

Caller attempt and accepted-interval limits, native source/body limits, cancellation checkpoints and two retained transactions per interval bound work/history. Only typed molecular gross-withdrawal refusal permits physical step-size retry. Source, phase, convection, inversion and conservation failures remain hard errors; late refusal publishes no gas or body endpoint. Convection callbacks must be repeatable, and their external side effects remain caller-owned.

These are local temporal estimates, not certified global error or material-source uncertainty bounds. The original homogeneous fixed-volume reservoirs, finite stationary bodies, prescribed apertures and qualified heat-transfer laws remain the model. Container mechanics, wall pressure/drag, rotational transport, flight delay and spatial gas resolution remain outside its scope. Six focused regression groups cover independent two-capacity convergence, opposing heat and canonical ordering, latent-plateau refinement, entering species and exact native-history replay, whole-window prescribed heat, fixed budgets, withdrawal retry and atomic callback/cancellation refusal.

## Nonuniform nonlinear contact heating through original finite thermal owners

The planar/axisymmetric nonlinear contact quadrature policy composes with the
existing original interface-temperature/phase/pressure refresh and explicit
contact heat partitions. Every electrical trial closes each source sample on
its refreshed external context; the accepted native integral enters original
thermal receivers once. Constitutive residuals are not heat credits, and a
representative contact seed is not used as a face-average nonlinear law.

The manufactured regression in
`tests/conduction_enthalpy/response/contacts/quadrature.rs` has an affine,
continuum-consistent nonuniform first endpoint with bulk losses 2 W and 1 W
and contact loss 1.5 W. The original interface has nonzero temperature response.
It checks applied/evaluated point sources, independent finite-capacity storage
and source work, the explicit heat partition, a second accepted thermal
interval, work budgets and late cancellation. Original caloric support,
phase/source ownership and transactional publication remain binding. These
checks establish numerical model/accounting targets, not measured-material
authority or a spatial/time/quadrature error certificate.

## Complete thermal assemblies for thermally live actuator supports

`circuit_enthalpy::motion::support::step_heat_assembly_candidate` and
`PreparedHeatNetwork::field_motion_support_candidate` connect the original
thermally live axial support to the complete existing `HeatAssembly`. The plain
material-network and electrical-contact entrypoints retain their existing API
and share the same implementation. Defaulted report/sample type parameters
preserve their original payload while complete providers retain their own
source observations and typed errors.

At every proposed finite-enthalpy endpoint, the original prepared assembly
resolves its bulk links, source-bound thermal and electrical contacts, closed
radiation, integrated/directional transport, and any composed coolant or
prescribed-reservoir sources. Its original material/caloric inventory and
electrical heat destinations admit the support; no independent heat map or
constitutive state is invented. That same endpoint feeds the existing section
Helmholtz/entropy law and native reciprocal field/body descriptor. Per-layer
reversible heat, bulk Joule heat and each interface's A/B heat partition enter
the original finite cells once.

The physical balance encloses caloric, support internal, field and kinetic
energy. Authored source/reservoir exchange remains external. Finite links and
closed radiation remain internal; the enclosure's numerical net power is
explicitly excluded from external heat and remains visible in the enclosing
energy residual. No balancing heat, stored-energy correction or second
mechanical/thermal solver is introduced. Source refusal, circuit-budget
exhaustion and cancellation return no thermal, electrical or mechanical
endpoint; retries start from the same accepted owners.

The original restrictions persist: wholly solid source-compatible layers,
fixed reference inventory/inertia, constant enclosed field coefficients,
affine modulus/declared expansion, fixed electrical geometry and generally
first-order endpoint-temperature coupling. Thermal link geometry, interface
pressure, exposed patch area and radiative view factors remain authored over
the interval. Support motion does not infer remeshing, changing contact,
coolant volume, phase transfer, arcing or nonlinear E/J constitutive response.

The focused `circuit_enthalpy_motion::support` regression definitions exercise
source-dependent thermal contact plus native radiosity, per-cell heat and an
independent total internal-energy account, transport changes reaching native
gap/momentum, and exact reduction to the original contact-support path when
extra heat mechanisms are absent. A complete-provider test composes the native
`fs-convection` coefficient with a prescribed fluid reservoir, checks its
external heat separately from closed radiation, retains the provider's sample
and error types, and injects late source/budget/cancellation refusals. These new
regressions were source-reviewed but not executed in this implementation
session because the execution environment became unavailable; compilation,
runtime acceptance and experimental validation are not claimed here.

## Actual-cell temperature continuation of loaded adhesive angular contacts

The rolling-dissipation and torsion heat owners expose attached, same-time
native temperature continuations. The exact original receiver identities,
material caloric objects, accepted time, source history and caller stamp are
admitted before either bulk temperature is read from its actual finite cell.
The independently supplied interface temperature drives its original adhesion
query. Original ordered cards, source applicability and angular calibration
remain bound through the sealed normal, angular and lifecycle continuations.

The native state and heat/source binding publish together after the final
caller checkpoint and cancellation check. Full spring history and all three
normal-release, angular-release and plastic-loss ledgers survive exactly, as
do their partitions, outside destinations, receiver bindings and clocks.
The next native motion uses the refreshed source and retained pair geometry.
No old loss is deposited again, and this coefficient update changes neither
body momentum nor caloric enthalpy.

Reports retain the original source temperature report, both observed bulk
temperatures, prescribed interface temperature and a named external driver for
signed normal-plus-angular parameter work. A finite nonnegative bound limits
its absolute total. The mechanical JKR potential is not a caloric entropy
model, so storage work is not silently removed from the receivers or declared
irreversible heat. Thermal branch changes, coupled area stability, changing
yield calibration and closed thermoelastic feedback remain outside this
attached-endpoint operation.

Four new focused native/thermal cases cover independent E/G/W expectations,
all three nonzero prior ledgers, retained memory and actual thermal/body state,
subsequent native advancement, stale history, wrong cards/clock, invalid or
exceeded work controls and final-publication cancellation. Pinned formatting,
recursive parsing and independent source review passed. These new cross-crate
tests were not compiled or executed because the session's dependency build
was blocked by exhausted memory/disk; the separate 20-test production angular
leaf harness passed and does not substitute for integrated execution.

### Complete sensitivities and pump-work objectives for modulated Hamiltonians

`ParametricHamiltonian::prepare_sensitivity` binds the original
`H0(x)+x^T Q(theta) x/2`, drive, complete H0 tangent and physical clock to the
shared fs-orbit factor and independent cached AFT trajectories.
`PreparedParametricSensitivity::solve_parameter` accepts simultaneous original
J/R/G/H0, symmetric Fourier Q, input Fourier and clock directions. Its effort
partial is `H0_xp+Q_p*x` at fixed phase; fs-phs owns the physical contractions,
and the existing sI clock column contributes `i*n*Omega_p*X_n` once. Pump/input
direction dimensions, original clock and retained bandwidth are admitted.
Directions may introduce pump harmonics above the base pump bandwidth when
retained by the orbit. Maximum direction matrix/phase/state scratch is reserved
before generic numerical preparation; existing source/callback storage/work
exclusions and hard modulation bounds remain explicit. No finite-difference
or quadratic-storage replacement is introduced.

`power_objective` and `prepare_power_adjoint` expose six independent original
accounts: mean pump/port/resistance power and corresponding energy per actual
cycle. Pump is `x^T Qdot x/2`, with full state gradient `Qdot*x`; port/loss use
`e=grad(H0)+Q*x` and the full Hessian `H0_xx+Q`. Generic two-grid checks cover
value and complete covector in the selected power/energy units, then each grid
combines its own direct objective mean with the cached dual. Automatic complete
directions retain signed pump, port and resistance laws separately. The mean-pump
explicit clock derivative is `(Omega_p/Omega)*P`; for pump energy per cycle it
cancels the explicit period derivative analytically. The solved waveform clock
chain remains. Port/loss work includes its ordinary period derivative.
Parameter queries reuse the original factor/dual without fresh objective
Hessians, orbit solves or waveform-direction solves; original source/parameter
laws still evaluate the cached states.

Three new focused groups and all seven original parametric tests pass natively.
The new evidence compares complete nonlinear-memory and off-diagonal phased-pump
directions with independently rebuilt neighboring HB sources and an independent
physical ODE integrator; evaluates six separate power/work neighbor oracles;
checks the isolated-clock identity and changing-period counterexample; admits a
new retained pump harmonic; and checks original-clock/bandwidth rejection,
bounded scratch refusal, incomplete material actions and cancellation/retry.
These fixed-order derivatives and accounts do not certify infinite harmonic
resolution, stability, passivity or a heat destination.

## Tetrahedral sampled contact response in finite thermal feedback

The same explicit contact quadrature policy is retained by
`SpatialConductionResponse` through the existing thermal source refresh and
frozen-source handoff. The native triangular contact response carries its
integrated loss and original point observations; existing
`step_response_contact_candidate` heat routing consumes that integral once
through the original receiving cells and specified contact partitions. The
representative contact source remains the external-context seed, while
constitutive closure and acceptance use every point's applied/evaluated source.

Bulk/contact source support, caloric phase ownership, independent supplied
electrical work, field-trial budgets and atomic continuation keep their
existing meanings. No face-average nonlinear resistance or additional heat
source is introduced.

A focused three-dimensional thermal regression reuses the original material
and contact fixtures with varying exterior voltage. Two accepted intervals
check each point's original source temperature and independent quadratic
current law, each body's finite-capacity storage and Fourier exchange, contact
heat partitions and total supplied work. Native triangular regressions and
the planar analytic thermal endpoint provide complementary checks. Integration
execution and mesh/time/quadrature error bounds remain separate validation.

## Calibrated spatial relaxation with finite temperature feedback

`rigid_damping_heat::relaxation::temperature` feeds actual finite-cell
temperatures into the original complete Kelvin/Maxwell graph. Its borrowed
`CalibratedRelaxationDamping` binds the original connector forms, original
branch factors/times/order and one existing material `TemperatureShift` per
connector. WLF and Arrhenius coefficients, reference temperatures and support
are supplied explicitly. Reusing shifted outputs as original calibration
refuses when the bound source differs.

The shared material convention is a_T=tau(T)/tau_ref. All original branch times
and the original Kelvin viscosity of each thermorheologically simple connector
receive that same factor. The low-level native held-kinetics bank keeps all
elastic factors and memory coordinates unchanged. A change in rates therefore
creates no stored-energy change, parameter work or independent heat impulse.
Native physical dissipation still integrates Kelvin and Maxwell losses along
the same complete coupled trajectory; original signed actuation remains
separate. No history is reset or reinterpreted.

Each connector names its original caloric CellId independently of heat routes.
Shared cells are sampled once per synchronization, and unchanged connector
temperatures reuse their admitted material coefficients. Same-time births
complete before source sampling. Every actual native thermal interval freezes
kinetics at its initial temperature; original material histories advance only
when their matching mechanical/caloric/observer candidate accepts. After all
final thermal events, the source admits both the calibration domain and all
effective Kelvin forms, branch times and reciprocal rates without preparing a
dense flow. Endpoint refusal discards the complete candidate.

The candidate exposes original accepted bodies/history, finite thermal state,
the shared energy report, actual endpoint temperatures and a_T values. Original
reference laws are retained for continuation. Cached additional storage is
linear in the already admitted connector and branch counts. There is no new
constitutive law, eigensolve, thermodynamic state transfer or source fit.

Temperature feedback is explicit and first order. Its max_dt_s refinement is
independent of mechanical geometry refinement and tight energy balance.
Vertical modulus shifts, thermal expansion, moisture calibration, entropy
models and measured source qualification remain separate physical inputs.

Verification: All five focused temperature-feedback groups passed natively:
independent WLF/Arrhenius isothermal body/history/heat/work response;
16/32/64-step self-heating refinement against an independent coupled physical
ODE with nonlinear finite calorics; exact continuation across original thermal
events; empty-memory parity with the existing calibrated Kelvin feedback;
and final-endpoint calibration/reference/cancellation/observer refusals with
unchanged original inputs and deterministic retry. All four fixed-history
thermal regressions and four driven thermal regressions passed after the shared
source dispatch changed. These runs execute the original TemperatureShift,
MBD/pHS, caloric timeline and observer numerical owners in a bounded build.
No full workspace or complete impact-host execution is claimed.


### Adaptive complete nonlinear heat assemblies

`advance_heat_assembly_adaptive_candidate` connects the original thermal-only
heat transaction to `AdaptiveConductionStepper`. Coarse and two-half-step
proposals use the same complete nonlinear assembly; only accepted fine
intervals contribute physical storage, energy and original caloric entropy.
Native events and explicit continuous-source transition times split the clock.
Cell inventories, caloric objects, topology and contact coverage remain fixed.
One source-evaluation cap includes failed, coarse, rejected and final queries.
Original source/material/energy refusals and cancellation remain typed terminal
failures; native numerical nonconvergence and temporal error can reduce the
interval. A late refusal publishes no partial window. Final source evidence is
refreshed after target-time events. Entropy retains unclassified source terms;
local temporal error control does not establish global or model accuracy.

### Adaptive original continuum momentum and caloric evolution

The opt-in fs-couple laminar-enthalpy feature exposes
laminar_enthalpy::adaptive through LaminarEnthalpy::advance_adaptive and
advance_boussinesq_adaptive. Both compose the original coupled native BDF1
step without replacing its convective momentum, source-rheology, same-endpoint
density-buoyancy, sign-split BDM transport, caloric inventory, or local/global
thermal acceptance owners. The original force, physical boundary and slip
routes remain held for the whole window. Material/exchange callbacks must be
pure and smooth throughout the window and valid on native nonlinear probes.
No thermal event strictly inside the window is consumed: original
BackwardEulerStep admission refuses before a native step. Events exactly at
the target remain pending in the original left-limit state.

Each comparison solves one coarse interval and two original fine intervals
from the same accepted joint history. The second fine interval uses the first
fine thermal state and actual BDM velocity. Only both fine reports and the
final original field/material state survive acceptance; no Richardson field,
interpolated material, or extrapolated heat is substituted. The existing
ConductionTimeOrder::First owner supplies the fine-error correction for the
actual two representable spans and the accepted/rejected step multipliers.
The final suggested_step_s can seed the next window without resetting the
controller; it remains a proposal subject to new admission, not an error bound.
A computed endpoint exactly one representable floating-point value below the
target rounds to that target. Actual spans, including this cap-rounding
effect, still enter all original solves, work and error arithmetic; no broad
epsilon drops an unresolved physical interval.

The joint local-error maximum compares original per-cell stored-energy
changes in joules, original caloric temperature, and physical affine BDM
velocity at every original triangle vertex. Energy differences are formed
before mass multiplication and relative tolerances scale interval changes,
so a formation-enthalpy datum, absolute Kelvin value or uniform velocity
translation cannot inflate an error allowance. The Euclidean affine error
attains its maximum at a triangle vertex; raw edge-moment orientation and
centroid-only sampling do not define the momentum norm. Original energy
coordinates still detect latent storage when temperature remains flat.
Separate coarse-versus-summed-fine discrepancies also bound integrated bulk
strain heat, routed slip heat, signed external stream enthalpy, other external
heat and mechanical buoyancy work. Cancellation between these accounts
cannot make inaccurate transfer reports pass merely because total thermal
storage is unchanged. These are local estimates; finite sampling does not
certify unresolved source oscillations, a global trajectory error, pressure
accuracy, or whole-system passivity.

Global native-step and endpoint-evaluation ceilings include coarse and
rejected comparisons. Every native call receives only the remaining
endpoint allowance, and a terminal failed call reports that upper allowance
beside completed native evaluation counts. Original source, nonlinear,
momentum, thermal, event and cancellation failures are terminal; only a
successful but inaccurate comparison may refine. The original thermal
inventory's substep ceiling and the explicit accepted-report cap both bound
retained fine intervals. The minimum step bounds full comparisons; accepted
halves and a target-clipped final comparison can be shorter. An unresolved
comparison at its allowed floor refuses rather than claiming convergence.

Only accepted fine reports enter compensated physical accounts. Final
storage is independently measured as the sum of original mass times each
cell's final-minus-initial specific stored energy. Its difference from
accepted strain/slip/stream/other heat must satisfy the separate whole-window
joule gate. Buoyancy work remains an independently integrated mechanical
observation and never enters that thermal balance. Coarse/rejected heat,
SIP/upwind/backflow/BDF numerical loss and partial accepted prefixes are never
published as physical heat or replacement histories. All input histories
remain immutable through late failure.

Five native regression groups are authored: reversible two-cell transport
against an independent backward-Euler recurrence and continuous exact
solution under two caloric datums; naturally loaded uniform Boussinesq
acceleration with independent accepted-interval velocity/work quadrature and
24 J original caloric heating; fixed-temperature flow whose time-varying
stream heat and explicit cooling cancel; bounded work/source/late-cancellation
rollback; and original pending-event left limits, decimal endpoint rounding
and accuracy-floor refusal. Bounded source/API review found no concrete blocker. All five changed
coupled Rust files passed the pinned Rust formatter and a subsequent
rustfmt --check with edition 2024 and skip_children=true. Those checks parse
and format the source; they do not type-check or execute it. Native FEM/thermal
compilation and integration tests remain unexecuted in this session. No passing
integration-test claim is made.

## Nonlinear membrane damping in the original exponential impact solve

ExponentialImpactSystem constructors now admit the existing MembranePotential
with its original bounded proportional material preparation.
The same fs-plate condensed film supplies the complete projected linear
tension/bending matrix and original nonlinear stretching energy. Uniform
proportional material spectra automatically retain every original mixed-pair
Maxwell history through ImpactSystem::with_body_relaxation. Existing felt,
elastic contact, reciprocal volumes, spatial damping and passive radiation
continue to use their original owners and the same accepted state transaction.

The quadratic block contains the complete film pencil and every normalized
memory's z²/2. The original nonlinear strain/history coupling remains implicit
in PreparedNonlinearZoh; no lagged strain or separately advanced material state
is substituted. A stable original facet-work callback preserves tiny stretching
work without subtracting large total energies. Original mapped callbacks still
use their bounded temporary derivative vectors; allocation-free complete
impact execution is not claimed.

All film histories count against the unchanged dense capacity before assembly.
Source-bound films retain their complete strain/history tangent envelope and
material band, additionally checked against the original 0.9*pi sampled limit.
Bare films without this original complete envelope refuse; the existing
explicit zero-arm spectrum admits bounded lossless films without dummy states.
Original dt/tau<=2, slope/history/energy, duplicate intrinsic loss, radiation
coupling and transactional cancellation gates remain. Exact linear propagation
does not certify nonlinear transients or permit coarser source data.

Three complete film integration groups are authored: independent
Duffing-Maxwell response/heat refinement, signed mixed-mode acoustic loading
with rollback, and complete-state/source admission. They remain unrun because
the full original impact dependency closure is unavailable in the bounded
native build. The separately exercised original film-work and native mapped
pHS seams do not stand in for that full integration run.

The physical scope remains the original moderate-slope, statically relaxed
in-plane film with a common proportional spectrum. This does not supply
heterogeneous viscoelastic condensation, finite-strain 3-D tissue, wrinkling,
in-plane inertia, exact shell or transverse-only support propagation, or
thermal feedback.


### Adaptive finite-gas windows with nonlinear wall feedback

The original finite-gas thermal window controller can resolve convection at
the actual wall and common gas endpoint in every coarse, fine and retried
interval. It retains source-energy apportionment, local h/u/T/phase error
estimation, accepted fine states and whole-window energy/entropy accounts.
Source-aware callbacks receive the actual interval input. Initial/final
admission hooks can refuse source/domain exits but cannot replace a native
physical candidate; the interval controller remains private. Existing outer
heat-bracket and gas-increment retry rules remain unchanged. Wall, material
and source failures propagate without publishing a partial window. Temporal
error control remains local and does not establish physical-model accuracy.

### Complete material-memory field and physical-power curvature

`PreparedStructuralInteriorRelaxationHarmonic::solve_with_frequency_curvature`
reuses the full original mechanical/pHS factor for every physical field's first
two frequency derivatives. The fixed original M/C/relaxed-K, interface and
Maxwell laws retain every viscous-strain state. Explicit source jets supply both
force and inward-flow first/second derivatives, including phase and amplitude
variation. Material impedance is not frozen and no modal truncation, sampled
loss surrogate or additional factorization is introduced.

The response exposes displacement, velocity, pressure potential, pressure,
interface force/flow and every branch stress curvature. Original kinematic,
momentum and branch equations independently check the second solve, including
twice the first-order state in the frequency action. Complete mean-storage
slopes include the original spatial storage and every material spring.

Each original arm retains its own spring-energy value/slope/curvature and
independent heat, structural work and acoustic work. Positive branch resistance
stays inside the normalized quadratic contraction of the canonical branch
effort; the recovery does not form a potentially overflowing resistance-times-
field product. Structural and instantaneous wall heat use the original
symmetric damping blocks. Gyroscopic and equal/opposite interface exchange
remain outside material heat.

Original structural/acoustic source work, interface exchange, instantaneous
heat, branch heat and both branch-work channels expose physical power value,
slope and curvature. Separate structural, acoustic, relaxation and total
balance methods retain signed independent accounts with no unity floor.
The existing complete preparation bounds cover the added O(n^2+b*n) recovery;
the native curvature owner separately admits its three solves. Cancellation
or original-equation/arithmetic refusal publishes no partial response.

Focused native tests use independently assembled original consistent P1
matrices and the actual complete mechanical/Maxwell owner. They compare an
analytic coupled cubic across frequency, source phase/amplitude jets and delay
invariance, relaxing wall admittance with separate instantaneous losses, and
manufactured mutual-memory/gyroscopic original rows. The numerical recovery
was exercised directly; these tests do not execute the outer FE mesh assembly.

## Sparse full structural/interior acoustic transients

`cavity_modal::interior::structure::sparse::PreparedSparseStructuralInterior`
advances the original heterogeneous pressure-potential/pressure finite-element
field together with every original structural displacement/velocity coordinate.
It accepts the original sparse consistent M/K/C forms and integrates the same
native normal traces as `StructuralInterior`. The conforming tetrahedral entry
shares the existing exact vertex, owning-volume normal, sidedness and actual
free-coordinate lowering; it never constructs a dense coupled model first.

Pressure remains `p = psi_dot`, particle velocity remains `-grad(psi)/rho`,
and B has acoustic rows and structural columns. At midpoint the complete
simultaneous unknown is [structural velocity, pressure]. The native bounded
sparse LDLT owner factors

```text
[ M_s + h*C_s/2 + h^2*K_s/4     -h*B^T/2 ]
[ -h*B/2                       -M_a - h*C_a/2 - h^2*K_a/4 ].
```

Only the acoustic equation's sign is changed to obtain this symmetric factor.
Both original physical equations, kinematics, source work and port work are
checked afterward. Neither a bare cavity stiffness nor a structural stiffness
is inverted. Original lateral elasticity, all off-diagonal inertia, every
pressure coordinate, cavity compression, constant potential and structural
rigid modes survive. The feature supplies an actual sparse path beyond the
existing dense exponential owner's 42-coordinate time limit.

Each candidate retains independent structural and acoustic kinetic/elastic
storage; original generalized-force and projected inward-flow work; intrinsic
structural damping heat; relative-wall admittance heat; and both separately
contracted B/B^T interface works. Signed exchange is positive from fluid to
structure and cancels from the total account. Polarized original forms compute
energy changes without subtracting large stored totals. Structural, acoustic,
interface and complete-system balances all gate publication, alongside the
uncondensed componentwise physical equations. No heat comes from an energy
deficit, and a failed/cancelled step leaves the accepted state and clock intact.
The returned state belongs to its exact prepared spatial/material/time owner.

This lane requires exactly symmetric original forms and memoryless passive
structural damping. Native scaled inertia admits SPD mass. PSD stiffness and
damping admission use a 4096*epsilon diagonal shift in numerical congruence
coordinates ONLY for the admission probe, allowing true null spaces without
altering any retained coefficient. Resolved negative directions and coupled
zero diagonals refuse. Actual storage and heat retain their signed contraction
with evaluated roundoff; no negative value is clipped into a passivity claim.
The existing exact owner remains available for gyroscopic structures and its
explicit hereditary owners retain their original internal histories. A
frequency-local impedance or damping sample does not authorize this transient
law. Molecular relaxation, thermal diffusion, moving material/geometry, ideal
pressure release and exterior radiation are not inferred by this entry.

Preparation admits source/trace counts, all retained sparse entries, every
staging extent, native symbolic/fill/frontal work and numeric factor payload.
Native sparse fill may still be large and refuses its authored limits. Caller
sources, geometry, O(n) candidate/observation vectors and allocator/runtime
overhead are additional to factor bytes. A step reuses the cached factor and
performs sparse solves/contractions without a dense mass inverse, new
factorization, material query or modal reduction. It allocates O(n) candidate
and audit vectors; it is not an allocation-free propagator. Cancellation polls
during original-form work and around the synchronous native solve.

Implicit midpoint is second order and passive for the admitted frozen linear
forms; it is not an exact wave propagator. Time-step phase refinement remains
separate from spatial h/p, material quadrature, source bandwidth and physical
model validity.

Focused native evidence uses the unchanged production sparse dynamics with
the actual fs-sparse CSR/bounded LDLT and fs-math arithmetic. Eight tests pass:
G1 sealed-piston/cavity analytic frequency and second-order refinement;
G0/G1 a 112-coordinate coupled field with both domain balances and separate
positive losses; retained rigid and constant acoustic modes; independent
constant-force and volume-injection work; full off-diagonal blocks compared
with independently solved original nonsymmetric endpoint equations; G3
coordinate-unit congruence; and physical/resource/cancellation/identity/retry
refusals. Three additional public FE/Cx integration tests define refinement
against the existing exact owner, an actual free tetrahedral structure plus
p4 pressure field beyond 42 coordinates, and original-source/symmetry limits.
Those three integration tests are authored and API-reviewed but unexecuted:
the complete Cargo dependency closure remains unavailable in this environment.
The native slice does not establish whole-crate/workspace compilation, a
measured performance ratio, continuum convergence or acoustic-program closure.

## Mixed gas-vapor nuclei through the shared thermal acoustic adapter

With `caloric-cavitation`,
`MaterialConductingVaporBubble::with_noncondensable_gas` attaches one
immutable original gas inventory and its explicit wet-mixture applicability
bounds to the existing native resolved-liquid trajectory. `gas()` exposes
that original inventory. All direct, compact-source and receiver advances
continue to use the same native proposal, liquid-displacement source and
original water acoustic medium. No separate gas/phase update or second
damping model is introduced.

The new `gas_vapor_cavitation` example accepts the existing explicitly
synthetic gas-source TSV and requires wet-mixture conductivity and molecular
diffusivity bounds as command-line inputs. It initializes a gas-supported
10 micrometre nucleus in 353.15 K water and authors a smooth rarefaction and
compression cycle. Output distinguishes water partial pressure, dry partial
pressure and total pressure, retains the exact dry mass, and reports changing
vapor inventory, liquid heat, emitted sound, gas storage/work and independently
reconstructed cumulative thermal/mass/source-energy balances. Physical
transport domains and the selected fast-mixing approximation are printed;
no parameter is silently inferred from unrelated dry-gas transport.

Four new coupled groups are provided for equilibrium support by original gas,
native/adapted mixed-trajectory and causal-source agreement, complete
observer refusal and final cancellation. They reuse original source-bearing
synthetic fixture cards and explicit wet-mixture bounds. The original
pure-vapor tests remain present.

Fresh focused original-source validation passed: material 33/33, native cavitation 105/105, coupled conducting 10/10 (including four new mixed groups and two existing source-parser tests), and coupled vapor 5/5, for 153 freshly executed tests. The native/coupled `cavitation`-only configuration also compiled with the mixed feature disabled. DSR, RCH and the documented DSR fallback are absent; these local one-job Cargo lanes use retained facades for the actual product sources and original numerical dependencies. No full-workspace or release-gate pass is claimed.

The actual `gas_vapor_cavitation` example completed 200 samples with `mixture-a.tsv`, authored wet-mixture conductivity lower bound 0.02 W/(m K), and diffusivity lower bound 2e-5 m2/s. It recorded 100 evaporation and 100 condensation samples; radius 8.48781..14.29779 micrometres; vapor mass 7.52400e-16..3.59147e-15 kg; and exactly unchanged dry mass 2.0853512985751484e-15 kg. Peak sampled receiver pressure at 0.1 m was 1.40411e-4 Pa. Maximum cumulative residuals were 7.25308e-20 J for combined thermal energy, 1.45545e-28 kg for vapor mass, and 7.23157e-24 J for original-source gas energy. These are numerical conservation results within the admitted common-temperature, rapid-mixing model using explicitly synthetic gas and authored mixture transport; they are not empirical water/air calibration or a spatial/time-convergence claim.

## Axisymmetric nonlinear current through complete heat assemblies

`conduction::axisymmetric::response::step_axisymmetric_conduction_response_candidate`
closes original meridional bulk/interface laws, their swept Joule work and any
complete `HeatAssembly` at one finite-enthalpy endpoint. Original thermal
contacts, integrated conduction layers, radiation and finite coolant links can
therefore feed back into the actual field. The returned candidate retains the
assembly's complete sample and typed failure; its default sample preserves the
existing material-network entry points through the same numerical driver.

Voltage, current-density, total-current, regulated-power and grouped
voltage/current/mixed electrodes retain their native nonlinear solve and work
account. Every actual field solve consumes the outer trial budget. Original
bulk and contact operating axes enter thermal boundary admission; graph
electrical bindings cannot disappear. Swept volume routes, source pins, phase
charts and contact heat partitions retain their existing admission. Joule loss
enters its receiving cells once, and reciprocal finite links and enclosure
radiation remain internal exchange.

A focused manufactured regression couples original nonlinear cylinder
transport to a sourced thermal interface, compares two accepted intervals
against independently eliminated body balances, and exercises boundary and
grouped-terminal drive. Source failures and exhausted/cancelled trials retain
both input histories. This composition adds no moving-interface, field-dependent
caloric or reversible electrical-storage law; independent time, field and
contact-quadrature refinement remains required.


## Scenario reservoirs with native spatial force histories

`scenario_open_boundary::lattice3::{step_forced_serial, advance_forced_serial}`
connects existing typed velocity and pressure sampling to `ForcedOpenBoundaryGrid3`.
Held or explicitly supplied endpoint force densities retain the native lattice-unit
and x-fast convention and remain inseparable from their populations. Both adapters
reuse source/frame admission, material scaling, full-vector speed bounds and
pressure-chirp resolution. Selected velocity or pressure profiles refuse even
when the sampled uniform target equals the dormant fallback; unselected profiles
remain active. The same velocity-profile ownership check now applies to generic
velocity-only rebinds and mixed serial/pooled steps.

The native step publishes reservoirs, populations, accepted force epoch, wall
exchange and step count together. Three focused regression groups specify
D3Q19/D3Q27 parity against original explicit native composition, independent
half-force momentum/reservoir closure and SI boundary moments, late refusal/retry,
and profile preservation. Pinned Rust parsing and focused diff checks pass;
Cargo typechecking/runtime tests have not run. No new force model, integrator,
thermal-work or physical-validation claim is introduced.

## Coupled normal, vector shear and twist at finite thermal temperatures

`friction_heat::adhesive::combined_temperature::propose_jkr_combined_temperature_heat`
continues the original attached, full-stick I/II/III JKR source and sealed
native four-coordinate history at actual finite-cell temperatures. The two
bulk temperatures, massless-interface temperature and common relaxed contact
area close in the same bounded thermal endpoint solve. Each trial requeries
the original selected bulk and interface claims at their independent Kelvin
coordinates, then reruns the original combined source and native admission.
Normal, both shear components and twist share one disk; the torsional-only
radius cannot substitute for a shear-loaded contact.

A private statically dispatched backend now serves both this entry and
`propose_jkr_torsional_temperature_heat`. It retains the existing finite
difference/Newton driver, trial budgets, cancellation, variable-area transport,
receiver/coverage checks, caloric registration and final-trial entropy audit.
The original public torsional error enum, variant imports, conversions and
display text remain intact. Neither adapter introduces a second contact model.

The returned candidate owns the successor source, native bodies and history,
thermal state, refreshed route and complete typed reports together. Coordinates,
rates, complete native mobility, pending COM/angular travel, geometry, mechanical
clock and step count survive the temperature phase. The caloric clock advances
to the already staged mechanical endpoint. Original incoming bulk source
temperatures must agree with the accepted receiving cells, and original ordered
body/card/registered caloric identities and the reaction plane must agree.

Normal elastic, vector-shear elastic, torsional elastic and surface-energy
changes remain separate in the original source report. Their signed sum stays
on an explicitly named, bounded external coefficient-work port. It is not
deposited as heat or inferred as entropy: this mechanical potential does not
supply a caloric constitutive law. The driver checks both the complete contact
plus finite-cell balance and absolute independent stage residuals. Conductivity
and other heat supplies remain with the original callback and thermal owners.

Three new combined integration groups specify independently solved simultaneous
thermal/area feedback and all four storage terms, unchanged native history and
subsequent motion, full-pipeline zero-shear reduction to the original torsional
owner, and atomic receiver/domain/work/area/trial/cancellation refusals with
retry. Together with the four new source and three native groups, this extension
adds ten authored integration groups. Pinned formatting, recursive Rust parsing,
scoped whitespace and independent source/API review pass. These integration
tests have not been compiled or executed: the complete dependency build is
blocked by the current workspace's exhausted disk and memory limits. No
whole-crate/workspace or numerical execution pass is claimed.

This entry holds admitted native geometry and attached full-stick topology.
It does not consume pending travel, correct a stale gap, supply moving thermal
faces, infer finite-rotation kinematics, cross a thermal instability or invent
slip/debonding heat. Endpoint area limits do not establish stability along an
unobserved temperature path.

## Exact-vacuum moving molecular receiver startup

`molecular_exchange::filling::step_candidate` joins an original finite moving donor to a receiver represented explicitly as evacuated geometry or an ordinary moving gas. An evacuated receiver supplies no return Maxwellian, temperature, pressure, velocity or seed inventory. One original drifting-Maxwellian sample supplies the same accepted species mass, vector momentum and complete crossing energy to both endpoints. Gross withdrawal uses only the original donor inventory.

The donor retains its original moving species-boundary solve. `MovingGasParcel::fill_evacuated_candidate` retains captured momentum and removes actual bulk kinetic energy through the original material owner's admitted-mass callback before its single caloric inversion. Thermal and common-boost terms are formed directly in the declared aperture frame, with original species formation references removed before arithmetic. No extra receiver-mixing heat is added.

Actual final species/total mass, inertial momentum and source-relative internal-plus-bulk energy are independently checked in aperture and laboratory frames. Native caloric defects, reconstruction corrections and frame conversion consume a fixed noncancelling allowance. Both positive endpoints must satisfy original source and drift support before atomic publication. Closed geometry and active zero transmission retain exact vacuum; filled receivers delegate exactly to the original reciprocal pair.

The result is a frozen-rate, first-order homogeneous mean-drift endpoint. It supplies no collisionless filling front, finite flight delay, container trajectory, remaining wall-pressure force or complete-vessel passivity claim. Original source/species and inversion limits bound work; cancellation and every late source, inventory or balance refusal return no accepted endpoints.

Five coupled regression groups specify independent monatomic/diatomic capture temperature and momentum, species separation with tangential drift and formation/common-boost invariance, exact reciprocal continuation, inactive/blocked vacuum preservation, and source/withdrawal/defect-budget/cancellation refusal. Two material groups specify the admitted-mass energy callback.

## Exponential material-memory loss in finite caloric stores

material::heating::exponential::HeatedExponentialImpactSystem routes selected
original Maxwell-arm losses from ExponentialImpactSystem into the existing
CaloricReservoir charts. The original branch resistance and normalized-memory
coordinate generate a checked positive witness in the complete current state,
including later acoustic coordinates. The native nonlinear held propagator
integrates each source square on its accepted complete trajectory. No midpoint
power, endpoint velocity, fitted damping ratio or energy-deficit heat estimate
is substituted.

Each branch has one explicit destination; multiple branches may heat one finite
specimen. Canonical branch order fixes accumulation and unknown/duplicate
routes refuse. Unselected modal loss, other material arms, radiation and felt
conditioning stay in the mechanical account. Branch stored energy is never
deposited as heat. Source observations are bounded by the original total
dissipation up to declared arithmetic roundoff, with the full mechanical balance
still owned by the existing acceptance gate.

Thermal trials reuse the same compensated cumulative deposit and original
enthalpy datum as the existing prepared-impact wrapper. Finite chart/phase
admission, cumulative inventory and step energy checks happen before mechanical
publication; only infallible caloric assignments follow it. A chart refusal,
Newton failure, source/budget refusal or cancellation therefore advances no
accepted mechanical state, history, sample or reservoir. Attachment and
continuation preserve current history, clock and Newton limits without replaying
prior heat. The wrapper implements the existing ImpactSource audio interface.

Cold preparation admits complete byte/work/dimension allowances for every
retained source integral. Warm host buffers and native integral factors are
reused; constitutive callbacks retain their existing allocation limitations.
This is fixed-coefficient one-way thermalization into finite uniform specimens.
It does not supply spatial conduction or temperature-dependent material
feedback. Nonlinear motion and heat still require temporal refinement; exactness
is limited to the quadratic held-input case.

Three focused original caloric-owner groups cover sensible/latent storage,
compensated tiny deposits at a large enthalpy datum, finite-chart refusal and
retry. Three full-host groups are authored for independent accepted-trajectory
heat, motion parity, continued attachment, canonical routes, rollback and retry.
The full impact-host groups remain unrun in this bounded environment; native
pHS and caloric source checks do not substitute for that integration run.

## Nonlinear electrical response in thermally live actuator supports

`PreparedMaterialNetwork::field_motion_support_conductivity_candidate` connects
the original prepared E/J material-response owner to the same live-support
enthalpy transaction. Initial source selection and circuit construction use
`prepare_conductivity` and `resolve_conductivity_seeds`. The existing
`solve_native` electrical iteration resolves bulk E/J and the actual reciprocal
field/body trajectory at each proposed thermal endpoint. The original prepared
transport owner then refreshes scalar or directional Fourier coefficients from
those same solved-field bundles before the shared heat residual is evaluated.
A thermal link cannot retain its seed-field coefficient after electrical
response changes.

`motion::support::step_conductivity_candidate` provides the same closure for
explicitly resolved original interface sources and A/B heat routes. Each
contact retains its own V/I/J coordinate and source orientation; the containing
branch does not replace the interface voltage. The transport callback receives
the admitted native operating point and may refresh its independently owned
heat paths and source observations. It must not advance history or add native
Joule/reversible heat a second time. Closed-enclosure numerical net heat retains
its explicit internal classification.

All entrypoints share the existing support Helmholtz/entropy law, reversible
heat, native moving descriptor and finite-enthalpy residual. Native
`FieldConductivityStep` reports retain applied/evaluated material and interface
residuals. Every inner descriptor trial consumes the outer circuit budget,
including discarded thermal stencils and final admission; narrowing remaining
work preserves authored fixed/adaptive source-coordinate relaxation. Source
refusal, nonconvergence, budget exhaustion and cancellation publish no thermal,
electrical or support/mechanical history. Default report parameters and the
previous plain-network/complete-heat APIs retain their existing reduction.

The support potential and fixed caloric charts still require electrically
independent consumed properties and applicability. A zero-field-valid modulus
or heat-capacity claim cannot silently remain frozen after E/J changes.
Stationary conductor geometry, fixed reference mass/inertia, constant field
coefficients and the original endpoint-temperature timing remain restrictions.
This does not infer field-dependent mechanical thermodynamics, motional
conductor work, arcs or interface storage. Complete heat providers retain their
own coordinate-admission rules; the new prepared frontend reuses the existing
bulk conductivity owner rather than weakening static heat-assembly admission.

Four focused regressions cover E/T-dependent resistance and thermal links through
the public prepared frontend; original nonlinear contact sources with independent
series laws, per-cell heat and total internal energy; exact accounting of all
inner native trials, late refusal/nonconvergence/budget/cancellation/retry;
reduction to the prior complete heat-support path; and refusal of electrical
applicability hidden in mechanical/caloric sources. Pinned rustfmt parsed the
changed Rust sources and the scoped diff check passed. These regressions were
not executed in this implementation session; compilation, runtime acceptance
and experimental validation are not claimed.

### Finite-amplitude internal nonlinear dissipation

`phs_orbit::dissipative::DissipativeHamiltonian` evaluates the complete original
ODE `xdot=(J-R) grad H-D(x,grad H)+G u(t)` through the existing fs-phs physical
field and nonlinear-resistance owners. Its `NonlinearDissipation` is the same
read-only law used by original time stepping and equilibrium linearization.
All material-memory coordinates remain present, including resisting flows
outside the external input map's range. No artificial physical ports, fitted
damping ratio, resistance secant, or discarded generated harmonic is introduced.

Native scaled seeded HB, adaptive retained/omitted-equation and waveform
refinement, and exact discrete variational shooting consume the same complete
tangent `(J-R) Hxx-D_x-D_e Hxx`. Independent state and effort directions enter
the authored resisting action exactly once. A passive nonlinear law can have
a negative differential slope at a driven bias; that slope is retained.
An incomplete/refusing flow, tangent or original storage action refuses the
result, with no finite-difference replacement of an authored tangent.

The independent instantaneous accounts are original stored energy, external
port power, native linear-resistance power, nonlinear `e dot D`, and
`grad H dot xdot`. Nonlinear loss is never inferred from a residual balance
or reported as external supplied power. The shared pointwise work admission
preserves roundoff-sized signed values; it is not a global passivity proof.
Raw signed R keeps its original semantics. Heat destinations remain with the
physical composition.

Preparation retains six state vectors, one physical-input vector and the
existing storage-tangent owner, with hard limits of 256 states/ports and an
authored action allowance covering every state direction. Adapter evaluation
does not allocate; source callbacks retain their own work/storage behavior.
Native solver budgets and cancellation boundaries are unchanged. Fixed
harmonic order and finite shooting grids supply numerical evidence, not
infinite-tail, primitive-period, stability or constitutive-validity certificates.

Five focused tests compare independent original nonlinear ODE integration on
three grids, generated third/fifth harmonics and explicit memory, separate
power means, complete state/effort chain rules, variational shooting, passive
negative differential slope, weak physical contractions, malformed data,
source refusal and cancellation/reuse.

### Original material-memory parameter response and physical work

`PreparedStructuralInteriorRelaxation::prepare_parameter_direction` binds one
complete original spatial M/C/relaxed-K/interface direction and the signed
stiffness, strain-projection and relaxation-time direction of every Maxwell arm.
The shared spatial assembler inserts both opposite B_theta blocks. The retained
original `LinearMechanicalSystem` supplies its cached inverse mass and canonical
arm map; direction preparation is reused across frequencies. A factor belonging
to another original preparation refuses the direction.

`PreparedStructuralInteriorRelaxationHarmonic::solve_with_frequency_parameter_derivative`
accepts the existing complete force/inward-flow source jet and reuses the SAME
factor for all five native state solves. Physical displacement, velocity,
pressure potential, pressure, interface force/flow and arm stress retain both
parameter and mixed derivatives alongside the complete frequency curvature.
Stress recovery includes the changing sqrt(k) multiplying canonical branch
effort. The original kinematic, momentum and branch equations are independently
checked, including M_theta, reciprocal B_theta, k_theta, b_theta, rate_theta
and every additional -i frequency-action term.

Each branch independently exposes spring-storage parameter/mixed derivatives,
complete heat and signed structural/acoustic work jets. The original decay-rate
direction contributes directly to heat; changing strain maps contribute to
branch work. Structural and instantaneous wall heat include their original
symmetric damping directions, while gyro/interface exchange stays outside heat.
Eight power channels and four normalized balance methods retain all source,
field, coefficient and mixed product rules. A zero-projection arm can develop
a nonzero stress derivative without inventing first-order heat.

The response borrows complete original field jets through
`field_parameter_jet` and scalar arm stress through
`branch_stress_parameter_jet`, for existing coherent receiver, extremum and
bandwidth queries. Coordinates, topology and arm count stay fixed. Missing
material/geometry laws, finite parameter changes and transient parameter-change
work are not inferred.

Existing adapter ceilings admit a complete parameter payload
`8*(64*d^2+40*n^2+16*b*n+160*d+160*n+64*b)` bytes and work score
`2048*d^3+512*n^2+1024*b*n`, with d=2*n+b. Native direction and harmonic
bounds remain separately enforced. Recovery allocates no additional dense
factor and performs no extra state solve; cancellation or source, binding,
capacity, arithmetic or physical-row refusal publishes no partial response.

Focused native tests compare an independently differentiated original
uniform-pressure mixed-material law and coherent source delay, complete
manufactured nonuniform mutual-memory/gyroscopic equations, every independent
storage/work channel and the zero-memory limit. The actual mechanical and
Maxwell owners execute; outer FE mesh assembly is outside this focused runtime
verification. Existing physical field/curvature regressions remain.

## Original longitudinal support damping in exponential impact

ExponentialImpactSystem admits SupportedStrings prepared with the original
longitudinal motion and explicit support geometry for every span. The original
canonical mass transformation, signed moving-end chord, both transverse
planes, all retained axial modes and every mean/cosine Maxwell history survive
in one state. Support reactions, absolute/relative drag, ground damping and
the selected constant axial viscosity retain their original complete
resistance matrix. Configuration-dependent finite-strain axial viscosity and
transverse-only support images refuse this fixed-resistance entry.

The exact quadratic block contains the original support spring, transverse
tension/bending and affine lift, plus the exact zero-constant linear part of
every original equilibrium and relaxing strain. All reciprocal support and
mechanical/history cross terms remain. This is a split of the original
polynomial, not a tangent fitted at a trial state. Cubic/quartic geometric
coupling stays in the existing implicit discrete-gradient equation. Original
storage, force, accepted material history and the work/heat ledger are not
replaced. The shared native quadratic-strain remainder computes equilibrium
work using the original energy coefficient directly, without square-root
renormalization or subtraction of large stored energies.

Cold preparation augments each existing squared support clock bound by
sum_j k_j * SupportedStrainMap::stiffness_bound(maximum_strain_j).
The existing map bounds the complete canonical strain-gradient norm squared
plus maximum_strain times the strain-Hessian norm on the original declared
slope/strain domain. The preexisting small-signal and Maxwell stress/history
bounds remain; no norm is subtracted to claim a sharper estimate. This is a
conservative frequency envelope, never a replacement physical stiffness.
Preparation is idempotent and does not widen any source frequency band,
change the source's exact dt, or relax its dt/tau limit.

Only exponential-prepared support images retain this domain metadata when no
Maxwell arm is selected. No dummy physical history or material law is added.
The original joint material owner checks the complete support/volume/radiation
envelope and original source authority; a zero-state selection cannot discard
those gates. Ordinary midpoint support preparation is unchanged. All actual
body, Kelvin, material and acoustic coordinates count against the existing
dense capacities. Original strain, slope, travel, tension, tensile-source,
history, energy and cancellation checks still own acceptance and rollback.

Four focused source/integration groups are authored beside this split:
signed original Q and stable nonlinear work against an independent two-word
polynomial; storage/resistance/history and envelope parity; complete material
response/heat refinement; and zero-history metadata, capacity, source clock,
volume/radiation and transactional refusal. They passed native parser and
formatter checks but remain unrun because the complete original supported
host closure is unavailable in the bounded build. No replacement physical
host was used. The separate nine passing native pHS groups establish the
shared work/split algebra and original Maxwell behavior, not a full supported
ImpactSystem integration result.

Physical scope remains the original moderate-slope, statically balanced
reference with its explicitly selected axial Fourier channels. Unretained
strain harmonics remain statically condensed by the original model; this is
not a claim of fully resolved axial waves, exact finite geometry or thermal
material evolution. Exact propagation applies to the complete declared linear
block. Nonlinear/contact accuracy still requires temporal and modal refinement.

### Pumped responses over multiple original source periods

`ParametricHamiltonian::over_periods(m)` borrows the unchanged original source
on an integer multiple of its current response period. Nested views multiply
the total authored-period count; `angular_frequency`, `period_multiplier`
and `required_harmonics` expose the response chart. Original Fourier input,
Q/Qdot coefficients, physical frequencies/phases, J/R/G, nonlinear H0, all
memory states and native/authored H0 tangents remain unchanged and borrowed.
Physical-time field, power and tangent evaluation stay identical even for
nonbinary multipliers. Only the response clock and required numerical
bandwidth change. Total period multiplier and mapped source bandwidth are
bounded by 256; invalid clocks, periods or bandwidth refuse before scratch
allocation or source evaluation.

Existing full seeded/adaptive HB, variational shooting, signed sidebands,
sensitivities and objective owners use the longer response clock and complete
response duration. There is no new solver or source-coefficient reconstruction.
`require_subharmonic_separation` requires at least two original source
periods and a complete retained waveform on the response clock. It compares
every proper synchronized repeat through the existing exact-integer Parseval
`repeat_rms` owner and positive per-original-state physical RMS tolerances.
Primitive authored-clock gcd includes nonzero input AND pump harmonics,
admitting pump-only excitation while refusing constant/padded source clocks.
Memory states count; DC does not inflate tolerances and tiny numerical tails
cannot establish a longer period. This is neither automatic branch finding
nor a general primitive-period, existence, stability or infinite-tail proof.

Physical input/modulation parameter directions retain the authored SOURCE
clock and have bandwidth multiplied by total periods.
`ForcedSensitivityDirection.angular_frequency_derivative` denotes RESPONSE
omega_p; source-clock change is m*omega_p, so the explicit pump relative clock
partial is omega_p/omega_response. Per-cycle objectives use the full response
period. Pump explicit period/clock cancellation remains analytic, with its
waveform chain retained.

Three focused tests exercise a biased saturated nonlinear half-frequency
response with relaxation memory against independently refined original ODEs,
full-cycle variational shooting and Liouville's determinant, six independent
physical power/work values, rejection of the ordinary shorter HB branch,
unchanged physical-time laws and authored Hessian under nested 3x2 views,
response-clock neighboring-source gradients, mapped-bandwidth refusal and
cancellation/reuse. Existing parametric and parameter-objective tests remain
applicable to the original one-period chart.


## Adaptive convection of multiple bodies sharing one finite gas

`FiniteGasBodiesEnvironment::advance_adaptive` and its source/wall-feedback
variants advance distinct original enthalpy bodies through one common gas
inventory. Every coarse, fine and retried proposal uses the original
simultaneous native solve. Temporal error checks each body's enthalpy,
temperature and liquid fraction independently and checks the gas's original
fixed-density internal energy and temperature once. Opposing transfers cannot
cancel component error. Only accepted fine intervals contribute physical heat,
elapsed time and retained entropy/source history.

Constant-power total window heats are apportioned through cumulative source
endpoints. Original body IDs, caloric objects, masses, geometry, environments
and gas source remain fixed. The checked wall API exposes actual interval
descriptors and admits initial/final states of every native proposal, including
coarse and rejected ones. Callbacks can refuse unsupported source states but
cannot mint physical candidates. Original wall-law work allowances remain
shared across bodies.

The native-attempt budget includes rejected and refused calls. An independent
body-step cap bounds retained component history and refuses before further
proposals when exhausted. Final component storage, gas storage, source
apportionment, matching and summed absolute local defects pass the common
whole-window energy allowance. Independent endpoint absolute defects also
prevent cancellation between component residuals. No partial physical window
is published on failure or cancellation.

This remains a fixed-volume, fixed-composition, well-mixed gas coupled to
independent native lumped bodies. It introduces no direct body/body contact,
spatial gas field, alternate caloric law, global time-error guarantee or
experimental qualification. Four G1/G3/G4 regression groups are authored.
The changed sources passed the pinned parser and an independent continuous
cooling oracle was checked; native compilation and execution remain pending.


### Nonlinear moving conduction and complete finite heat

Planar conduction::response::motion and tetrahedral
conduction::spatial::response::motion now close the original scalar electrical
sources at T, phase, conductor-frame E' = E + u cross B and actual J inside the
existing finite-enthalpy endpoint. Every thermal proposal refreshes the same
complete pinned bulk and stationary-contact sources; the native response owner
closes their electrical coordinates. Original caloric charts, volume routes,
source applicability, contact heat partitions, native constitutive/contact
refinement and total field-trial budgets remain binding. Rejected stencils do
not advance either history or turn numerical residuals into physical heat.

The accepted moving reports retain native laboratory and conductor-frame
fields, current, Lorentz force, original source-closure diagnostics and separate
work accounts. Only native local conductor-frame Joule loss and integrated
stationary-contact loss are deposited into their original finite receivers,
once. Terminal input and mechanical output remain signed and independently
recovered: net_supplied_work_j is electrical input minus Lorentz output. It
closes Joule heat during motoring, generation and braking, including negative
terminal input or zero terminal supply. Positive mechanical work leaves the
electrical/thermal system. The original thermal driver receives this net work,
not laboratory E.J as an additional heat source.

step_moving_conduction_response_candidate and
step_spatial_moving_conduction_response_candidate accept any complete
HeatAssembly and retain its typed source sample and errors. Thermal contacts,
conduction layers, radiation, reservoirs and finite links evaluate the same
endpoint. Both bulk and nonlinear-contact operating axes participate in
field-independent thermal-source admission. Enclosure radiation and reciprocal
finite links remain internal exchange, and their existing closed-balance gate
is retained. Other graph electrical bindings refuse without an owner for their
work. MaterialNetwork conveniences use the same respective endpoint paths.

Planar and 3D motion accept held boundary voltage/distributed-current traces
and native direct grouped voltage/current/floating/mixed electrodes.
Original electrical-contact adjacent cells must remain stationary, even with
zero motional EMF. Motion elsewhere remains a prescribed fixed-domain snapshot:
occupancy and geometry must stay represented by the original cells. These
entry points supply no mesh advection, induction, magnetic storage, momentum
evolution or sliding-interface work law. Continuation accepts the native source
owner and finite thermal history together; source failure, exhausted budgets,
failed energy gates and cancellation accept neither.

Focused regressions cover nonlinear E'/J laws, native force and per-body energy
across two intervals, source/speed sensitivity, generation and braking,
stationary contact heat/velocity guards, complete radiation/contact provider
samples and typed failures, shared budgets, cancellation and reusable input
owners. Pinned rustfmt parsed/formatted the changed Rust files, and focused
cross-module source/API reviews found no blocker. Extracted independent
manufactured-reference functions were compiled and executed for six planar
complete-assembly endpoints, ten 3D moving-source endpoints and six 3D complete
heat endpoints. Those arithmetic references passed their stated gates; the
worst planar two-body balance defect was below 1.4e-13 J. Full crate compilation
and the newly authored integration regressions have not been executed in this
environment. These checks do not establish mesh/time accuracy, uniqueness,
global convergence or empirical material validation.

### Endpoint mechanical sources in original laminar enthalpy coupling

LaminarEnthalpy::advance_driven and advance_boussinesq_driven resolve an owned
LaminarMechanicalDrive from each actual native caloric candidate and explicit
LaminarDriveQuery. The same immutable snapshot supplies spatial force,
prescribed fluid/wall/reservoir velocity and literal supplied traction to the
original BDF1 momentum owner; its material and thermal callbacks observe that
same source. Original boundary kinds, SI slip coefficients, constitutive
projection, pressure datums, objective strain, quadrature and sign-split
enthalpy flux remain native. No spatial interpolation, second momentum solve,
lagged thermal actuation, or replacement caloric law is introduced. Boundary
kind changes and source discontinuities are caller-authored events; source
callbacks must be pure and valid on nonlinear derivative probes.

The prior held-input and Boussinesq APIs share the extracted original thermal
residual/acceptance core and original constitutive helpers. The accepted
LaminarDrivenExchange retains the exact mechanical source object, query and
thermal payload. Changed routes are checked against the actual endpoint
boundary kinds and original finite receivers. Every source, native solve,
heat, route, work or cancellation refusal leaves both input histories intact.

The optional mechanical receipt independently evaluates additional body-force
work with the original six-point triangle quadrature and accepted BDM basis.
Accepted native reactions supply work conjugate to the actual prescribed
velocity/Navier-slip wall data, supplied open-traction work and separate
dissipative-backflow stress work. These retained powers use the same original
boundary samples and endpoint material; no boundary callback is queried again.
Their sign is positive into the fluid. All-velocity components use the returned
zero pressure gauge, while traction components preserve their native datum.
A velocity condition may carry throughflow and is not inferred to be a wall.
Native retained transport supplies signed outward carried kinetic energy and
separate numerical upwind loss. Extrusion depth and actual duration convert
the SI per-depth observations to joules. Buoyancy keeps its independent
source-density receipt. Physical Navier-slip heat follows its existing explicit
fluid/wall partition exactly once; none of the mechanical observations or
numerical losses is deposited as heat. These are discrete source exchanges,
not pointwise traction estimates or complete mechanical energy/passivity
certificates.

Focused native regressions are authored for an endpoint-temperature actuator
with an independent uniform-acceleration and linear-pressure solution, both
ordinary and original bounded-density Boussinesq coupling, independent signed
body/open-traction/buoyancy work, separate 24 J original caloric heating,
retained original source identity, and source/late-cancellation rollback.
The changed Rust files pass the pinned parser/formatter. Native FEM/thermal
compilation and integration test execution remain unverified in this session.

## Sparse original structural and acoustic Maxwell histories

`cavity_modal::interior::structure::sparse::PreparedSparseStructuralInterior::
new_relaxing` and `from_tet_elastic_relaxing` extend the existing complete
sparse owner with the original `fs_modal::relaxation_modes::
MaxwellDampingBranch` data. Each source retains its full authored
[u, pressure-potential] strain direction, positive stiffness k and positive
physical time tau. The constitutive equation remains

`sigma_dot = k * b^T * [v,p] - sigma/tau`.

The spatial stiffness is RELAXED; adding the branch stiffness to K as well
would double-count it. A structural-only map supplies hereditary material
stress. A pressure-potential-only map supplies the causal acoustic wall
admittance `k*tau/(1+s*tau)*b*b^T`. Original zero maps retain their actual
stress, storage and free relaxation rather than disappearing. Branch order
is source order. A mixed structural/acoustic direction refuses intact
because it does not preserve this real symmetric signed factorization;
the existing general exact and harmonic owners remain its implementation.
Gyroscopic or nonsymmetric spatial forms also retain their existing scope.

The same physical M/K/C/B forms, geometry adapter, source projection and
full spatial coordinates are used. The memoryless requirement on the
original C blocks still applies; the new stress states carry the hereditary
response separately. No material law is inferred from a frequency sample,
damping ratio, scalar loss fit or post-hoc energy deficit.

For half-step t=h/2 and nu=tau/(tau+t), a branch adds the scaled unknown
sigma_mid/sqrt(k*nu). Its diagonal is -1 for structural memory and +1 for
acoustic memory; its sparse column contains the original direction times
t*sqrt(k*nu), with the corresponding physical sign. Thus one native bounded
LDLT factor advances all velocities, pressures and memories simultaneously.
The branch is not expanded into a dense b*b^T update. This is significant for
nonlocal original strain maps: the 112-spatial-coordinate test admits
48-entry and 64-entry directions with a 3000-entry staging ceiling, while
their dense outer products alone would require 6400 entries.

`state_with_stresses` imports every actual initial stress in source order.
The existing `state` and `zero_state` explicitly initialize all admitted
branches at zero stress. Returned states preserve original fields, source
clock and stress histories under one owner identity. State reuse across
another owner, source or timestep refuses. No relaxation law, factorization
or material query is reconstructed during a step.

Storage includes sigma^2/(2k). Each step separately reports every branch's
initial/final storage, polarized storage change, structural or acoustic
work, independently evaluated h*sigma_mid^2/(k*tau) heat and original
constitutive-equation residual. Structural, acoustic, interface, aggregate
memory and total work accounts all gate publication; every individual
branch account also gates. Momentum equations include the original
b*sigma_mid reaction, and audits reconstruct midpoints from accepted
candidate endpoints. No negative heat is clipped, and no physical loss is
manufactured from a balance residual. Cancellation or refusal leaves the
initial fields, histories and clock unchanged.

`SparseStructuralRelaxationLimits` bounds branch count and the full source
coefficient scan before sparse lowering. Retained strain entries and
augmented staging also obey the existing matrix-entry limit; the full
augmented factor obeys the native symbolic/fill/frontal/work and payload
bounds. Dense caller source arrays, allocator metadata and accepted state
vectors are additional. Congruence or nonzero material coupling that cannot
be represented refuses. Cancellation polls through source directions,
augmented assembly, branch RHS/recovery and physical audits, as well as
around the original synchronous native solve.

The method remains passive second-order implicit midpoint on admitted
frozen linear forms. It is not an exact or L-stable relaxation integrator;
short relaxation times need their own temporal refinement. Existing
spatial h/p, source bandwidth and constitutive validity requirements remain.
This entry does not add thermoelasticity, molecular/thermal volume memory,
nonlinear or changing material laws, mixed cross-domain Maxwell channels,
exterior radiation, arbitrary positive-real wall realizations or measured
material validation.

Seven new G1/G3 numerical source tests cover dark original stresses;
independently prescribed continuous fields and work with both solid and
wall memory; full off-diagonal original endpoint equations solved by an
independent nonsymmetric Gaussian elimination; 112-coordinate full fields
with nonlocal maps; original material strain-unit congruence; analytic
sinusoidal wall admittance and retained potential gauge; and source,
budget, cancellation, foreign-owner and exact-retry refusals. Temporal
continuous-field errors decrease 1.58353e-5, 3.96088e-6, 9.90350e-7; the
wall-admittance refinement decreases 1.64801e-4, 4.12133e-5, 1.03041e-5.

The final focused native run passes those seven tests, all eight existing
sparse-field tests and one comparison with the pinned pre-change owner:
100 forced no-memory steps retain identical field, interface and physical
work bits. The build uses full actual fs-sparse and fs-math with the final
production sparse numerical owners. One additional public FE/Cx integration
test compares the original Maxwell trajectories against the existing exact
coupled owner. It is authored, formatted and API-reviewed but unexecuted
because full Cargo dependency closure is unavailable; the three previously
authored FE/Cx tests also remain outside the native numerical slice. This is
not whole-crate/workspace Cargo success or a measured performance ratio.

## Full original-coordinate structural reaction through FMM/Robin radiation

`vibroacoustic::iterative::StructuralRadiation` couples a real frozen original
mechanical M/C/K pencil to a closed exterior surface through the existing
Burton-Miller `IterativeRadiation` owner. Original full consistent inertia,
off-diagonal stiffness, viscous damping and skew gyroscopic velocity coupling
remain in the native `PreparedHarmonicResponse`. One real panel-normal trace
per supplied coordinate retains zero/dark rows as well as visible motion.

Under peak exp(-i omega t), the original equations are
`H q = F - T^T W p`, `H = K - omega^2 M - i omega C`, and
`v_n = T (-i omega q)`. Native structural mobility `Y = -i omega H^-1`
therefore supplies exactly the nonlocal normalized Robin action
`beta p = rho*c*T*Y*T^T*W*p` and free drive `T*Y*F`. Pressure force and velocity
use the same original trace and panel area map. No modal projection,
independent-radiator power addition, diagonal impedance replacement or
resymmetrization of gyroscopic response is performed by this adapter.

Preparation factors the mechanical pencil once and constructs its full
coordinate mobility using structural unit-force solves. It performs no
acoustic unit-coordinate column solves and forms no acoustic N-by-N matrix.
Each forcing reuses one prepared native FMM/Robin operator, obtains one
acoustic solution, then solves the original mechanical pencil against the
actual external force minus recovered pressure traction. Native pressure and
outward velocity remain available as the original `RadiationSolution` for
downstream field observation.

The complete returned work account distinguishes original external input,
structural viscous dissipation, signed structure-to-fluid interface work and
signed radiation work. Independent contractions check the reconstructed
physical boundary velocity, reciprocal interface work and total input =
structural heat + outward sound. A strictly negative native outward-rounded
radiation dot-product enclosure refuses rather than being clipped into
passivity. That enclosure covers floating-point dot evaluation only; it is
not a discretization or FMM error bound, and a pointwise check is not a
global passivity certificate.

Caller controls separately bound coordinate-by-panel trace entries, the full
coordinate mobility and added numerical payload
`16*q*q + 192*q + 160*N + 4096` bytes. Borrowed source traces, allocator
metadata and the original mechanical/FMM owners are excluded from that added
bound. A nonlocal action costs O(q*N + q*q) with O(q) scratch; it allocates
bounded scratch, not a reusable allocation-free workspace. Native panel
quadrature, directional/W5a, near-block and GMRES resource limits remain
unchanged. Cx is polled around synchronous preparation and inside the
existing controlled solve; native CancelGate polling remains active within
FMM stages. Source owners remain immutable after cancellation or refusal.

This route requires an invertible bare mechanical pencil. An exactly
singular lossless bare resonance refuses even if radiation could regularize
the jointly loaded problem; existing loaded dense owners retain that scope.
The adapter does not turn a frequency-dependent material sample into a
frozen M/C/K law, build a tensor-Prony shell material reduction, add incident
acoustic illumination, infer structural traces, or certify structural/modal
truncation, continuum geometry or broadband/causal validity.

Five focused native groups execute the complete new coupling source against
the original mechanical PSD/LU/refinement, FMM, triangle kernels, Robin,
GMRES, cancellation gate and geometry numerics. The independent pulsating
sphere at ka=0.8 refines from 80 to 320 panels with the admitted TriangleNear
rule: velocity error falls from 4.3718% to 1.0934%, and radiation-power error
from 2.7385% to 0.6183%. A separated-body directional full-Triangle fixture
with gyroscopic coupling agrees with dense acoustic columns composed with
the original loaded mechanical owner: relative velocity error 5.104e-10,
pressure error 3.671e-10, and total work defect 1.027e-15. Other groups check
coordinate congruence and structurally coupled dark motion, equality of
nonlocal surface work and original structural dissipation, and singular,
resource, invalid-force, cancellation, deterministic retry and zero-drive
limits. All five groups pass. The native slice replaces only Cx arena/RNG
construction in these new tests and omits unrelated dependency construction
surfaces/modules; full Cargo/workspace compilation remains unexecuted.


### Scheduled mechanical events in coupled laminar evolution

`LaminarEnthalpy::advance_scheduled` and
`advance_boussinesq_scheduled` drive the existing joint BDF1 momentum/enthalpy
step-doubling controller with the original `LaminarMechanicalDrive` snapshot.
Every coarse and fine thermal endpoint, including nonlinear probes, resolves
its own source against the actual candidate thermal inventory. The resulting
native force, velocity/slip/open boundary callbacks, source rheology and heat
exchange all see that same owned snapshot. The accepted report retains the
actual snapshot, query and original thermal payload; equal endpoint times do
not authorize reusing a snapshot from a different candidate.

A caller-authored finite, nonnegative, strictly increasing mechanical event
list clips comparisons at the exact event. The interval closing an event
uses its preceding segment, and the next interval uses the following segment.
Past/future events may remain in the list when resuming from accepted native
histories. Smooth data and fixed boundary kinds within each authored segment
are a source contract; discontinuities and kind changes require explicit
events. One-ULP event/target tails snap only to that exact stop, while all native
steps, work and error corrections use the actual represented durations.
An event-clipped short interval does not permanently collapse the next
controller proposal. Original thermal inventory events remain subject to
their independent native admission: mechanical events do not activate,
consume or reinterpret thermal events.

The shared physical comparison checks body work, generalized prescribed-fluid/
slip-wall work, supplied traction work, dissipative-backflow stress work and
signed carried kinetic energy separately from the existing bulk/slip/stream/
other heat and buoyancy channels. Each channel checks both its signed
coarse/fine defect and the defect between coarse magnitude and the sum of fine
magnitudes, so opposing accounts or reversals cannot hide unresolved transfer
behind cancellation. Only accepted fine reports enter compensated totals.
Numerical upwind loss remains a separate observation and does not enter the
physical work norm. No mechanical observation creates a heat source. Final
thermal storage is independently reconciled against the accepted original
heat and stream exchanges.

All coarse, fine, rejected and terminal native calls share the original work
ceilings. Material/drive/heat refusal, native physical failure and cancellation
are terminal rather than invitations to retry invalid physics with a smaller
step. Failure publishes no accepted prefix and never mutates either input
history. Existing held-drive adaptive APIs use the same owner with an empty
mechanical schedule and preserve their force/boundary semantics.

Focused regressions exercise steady velocity/temperature with opposing
time-dependent body/pressure work that still demands temporal refinement;
exact decimal-time velocity reversals and original two-cell enthalpy
recurrences, including resumed-window parity; a temperature-driven actuator
with original-density buoyancy and independent BDF kinetic-work accounting;
and source refusal, cancellation and global work exhaustion after accepted
internal prefixes. The six changed Rust sources passed the pinned parser,
formatter and `rustfmt --check`. These new native FEM/thermal integration
tests remain unexecuted; parsing and static API/oracle review are not runtime
validation.

## Simultaneous structural and FMM radiation at bare resonances

`vibroacoustic::iterative::joint::PreparedJointStructuralRadiation` retains
every original structural velocity and surface pressure in one native GMRES
system. With peak exp(-i omega t), its equations are
`A p - B(i omega rho T v) = 0` and
`[C + i(K/omega - omega M)] v + T^T W p = F`.
The original native Neumann FMM applies A and B, including its selected
triangle/directional/near corrections. The original all-prescribed mechanical
owner evaluates Hq and viscous work without factoring H; these trial coordinates
are Krylov unknowns, not added physical supports. No bare mechanical inverse,
dense acoustic column matrix, modal replacement or artificial loss is introduced.
Actual radiation can therefore resolve a bare structural resonance. A truly
unresolved forced dark resonance still refuses; convergence is not a general
uniqueness or rank certificate, including for an unforced nullspace.

Original M/C/K admission, full gyroscopic entries, exact outward closed-triangle
admission, reciprocal trace/traction maps and the native near-block pressure
preconditioner remain authoritative. Positive solver-coordinate scaling uses
`d_j^2 = omega M_jj + abs(K_jj)/omega + abs(C_jj) + rho*c*sum_i(area_i*T_ji^2)`
and pressure scale `sqrt(rho*c/area)`. It changes no physical coefficient.
The solve returns the original pressure/velocity field, all structural motion,
signed interface and radiation work, independent viscous heat, a true Euclidean
residual of the scaled joint system, and independently recomputed force/acoustic
block defects and work balances. The force-block defect is in the same D^-1
physical coordinate metric and is explicitly not a coefficientwise backward
error. The acoustic defect is `norm(Ap-Bq)/(norm(Ap)+norm(Bq))`; neither defect
certifies spatial quadrature, FMM error or forward accuracy near a loaded pole.

Added payload admits `16*(N+n)*(restart+20) + 8*restart*(restart+4) + 512*n
+ 160*N + 48*(max_iters+1) + 4096` bytes, separately from the original
mechanical/FMM bounds, borrowed source data and allocator metadata. GMRES counts
actual Arnoldi updates, including partial final cycles. Native Cx/CancelGate
checks and a per-composed-matvec predicate preserve immutable retry after refusal.
The source is still a frozen real M/C/K law at one positive frequency; this
owner does not infer dispersive materials, incident illumination, mechanical
traces, startup transients or a coupled adjoint.

Four focused native groups pass with the actual original mechanical
PSD/LU/refinement/action, full FMM, triangle kernels, geometry, near-block and
GMRES owners. At an exactly singular bare scalar resonance, the independent
pulsating-sphere velocity error falls from 8.8931% to 2.1580% under 80-to-320-panel
refinement; signed work closes within 7.3e-14. A separated-body directional
full-Triangle gyroscopic fixture agrees with dense acoustic columns composed
with the original loaded mechanical owner to 5.104e-10 in velocity and 3.671e-10
in pressure, with 7.91e-15 total-work defect. Other groups cover original-unit
congruence and coupled dark motion, genuine forced dark singularity, invalid
mass/force, open geometry, exact added-memory and iteration limits, cancellation,
bitwise retry and zero forcing. Test execution narrows only executor context
construction and unrelated dependency/module surfaces; no physical numerical
owner is substituted. Full Cargo/workspace compilation remains unexecuted.

## Original molecular sources in full pressure-volume evolution

`cavity_modal::interior::heterogeneous::molecular::MolecularBulkInterior`
connects original `RelaxationState` sources to the complete sparse pressure
space. A complete disjoint partition of original tetrahedra retains every
borrowed source and its exact positive compressibility realization. Original
source density constructs `K_rho`; frozen `beta_inf` constructs `M_inf`, and
the pressure histories supply the relaxed increment once. There is no
homogeneous replacement, equilibrium-mass double count, attenuation fit,
extra Stokes--Kirchhoff term or inferred bulk viscosity.

Each region has frozen temperature/composition and the same stationary
background pressure. The source's coherence and material domain remain
binding. Up to 256 original regions are admitted; the FE preparation supplies
explicit sparse and full-history resource limits. Original spatial receiver
and source projections remain available through the retained heterogeneous
model. Subsequent propagation queries no material and retains all histories,
with separate physical molecular storage/heat and boundary/source work.
This is a linear stationary molecular-pressure model; resolved shear and
Fourier heat conduction remain on their existing thermoviscous owner.

Two focused native integration tests pass using original public molecular
law construction and actual FE assembly. They verify original source identity,
independent frozen-volume mass and variable-density stiffness, exact partition
and pressure compatibility, and full-field relaxation of a closed gas against
the independent single-reservoir caloric solution. Pressure errors at 20, 40
and 80 time steps were `2.420e-4`, `6.048e-5` and `1.512e-5` Pa, with separately
integrated nonnegative heat closing physical energy. The tests use the bounded
native dependency closure documented by the FE bulk owner, not a full
workspace/Cargo or experimental validation claim.

## Two finite rigid walls sharing the original liquid history

`free_surface_wall::planar::pair::FreeSurfacePlanarWallPair` retains one original
source-bound liquid and two finite rigid bodies on distinct opposing flooded
planes. One shared native pilot and one final replay provide both three-mode
wall responses. The original mass-whitened midpoint solve and native rigid
integrator retain each body's two tangent translations, principal normal spin,
raw and physical guide reactions, orbital-plus-spin momentum and mechanical
work. Both bodies and the one liquid clock publish together. Cell work counts
two complete-domain trials total, not two per body; per-wall receipts contain
small mechanical records and no duplicate liquid populations.

Held applied tangent forces and axial couples, physical tolerances and work
bounds remain explicit. Each original wall's final native replay and actual
body endpoint must pass before combined momentum, common-origin angular
momentum and signed-work checks admit the pair. Cross-wall influence travels
through the evolving native liquid over subsequent ticks. Fixed fully flooded
planes, isothermal flow and original source/VOF limits remain; this adds no
normal displacement, tilt, finite edges, cover/vacate or whole-fluid energy
conservation claim. Guide reactions and unrelated bulk forces retain separate
physical ownership.

Validation: four new pair tests and four original isothermal planar tests pass
against the selected original liquid, material, cancellation, rigid-body and
linear-algebra sources. They include an independent six-mode analytic solve on
all axes, actual body linear/spin/orbital momentum and kinetic-work balances,
pressure and ordering covariance, 80-step fluid-mediated momentum transfer and
atomic second-wall refusal/retry. The one-tick oracle compares every population
against one original native replay and checks the coupled owner's physical time
and accepted-step count; SparseGrid3's separate standalone sweep counter is not
a free-surface transport clock. Full Cargo remains blocked in asupersync by the
8 GiB environment limit. Thermal regressions and the full workspace gate were
not executed; their numerical expressions are unchanged by the helper split.


### Self-sustained regimes with original nonlinear internal loss

`phs_orbit::DissipativeAutonomousHamiltonian` composes the original pHS,
physical feedback and shared `NonlinearDissipation` law as
`xdot=(J-R)e-D(x,e)+G*u(x,G^T e)`. All material-memory coordinates remain,
including internal resisting flows outside G's range. The complete tangent
`(J-R)Hxx-D_x-D_e Hxx+G du_total/dx` uses the authored TOTAL feedback
derivative and independent state/effort loss directions exactly once. The
shared compensated physical-field owner performs each final row rounding.
Pointwise passive D may have a negative differential slope, which is retained.
Stored energy, physical feedback supply, signed native R work, nonlinear
`e dot D` and independently contracted energy rate are separate accounts.

Native seeded/scaled autonomous HB solves amplitude and natural frequency
together with a full-waveform phase condition. Native adaptive HB refines
retained/omitted equations and successive waveform/frequency. Variational
shooting solves the period against the original ODE and differentiates the
same RK4 map, retaining every memory and phase multiplier. Constructor limits
are 256 states/ports and all n declared loss tangent actions; scratch is
`6*n+3*m+n*m` scalars plus the existing storage-tangent owner. Warm adapter
source evaluation allocates nothing; callbacks own their resource behavior.
Solver controls bound additional work. Refusal/cancellation publishes no
partial cycle. The source-bound sensitivity and adjoint preparation below
retains the complete nonlinear-loss parameter law and unknown-frequency chain.

Three G1/G3 tests compare the complete tangent and independent work, a
nonlinear self-sustained mechanical/Maxwell cycle against independent RK4
equations on three grids, generated third/fifth harmonics, independently
determined natural frequency, separate mean supply/heat and variational
period refinement on three grids. They also cover incomplete tangents,
active loss, action limits, cancellation/reuse and equilibrium refusal.
These checks do not certify infinite harmonic tails, primitive periods,
global attraction, material validity or hybrid/nonsmooth source laws.


### Nonlinear internal loss on the original pumped source

`ParametricHamiltonian::with_dissipation` consumes a prepared source and
attaches the same `NonlinearDissipation` used by original pHS time stepping
and equilibrium linearization. Its field is `(J-R)e-D(x,e)+Gu`, with complete
pumped effort `e=grad H0+Q(t)x`. Original pointwise work admission evaluates
actual `e dot D`, retaining signed roundoff and negative incremental slopes;
this is neither global nor incremental passivity certification. Every memory
coordinate remains, and D need not lie in the external input map's range.
The complete tangent `(J-R)(H0_xx+Q)-D_x-D_e(H0_xx+Q)` reaches existing HB,
adaptive HB, variational shooting and signed-sideband owners without an
effective damping, reduced state or finite-difference replacement.

`ParametricPower::nonlinear_power` separately reports nonlinear loss; pump,
port and native resistance remain independent accounts. Longer-period views
borrow the same law and unchanged physical-time source. Preparation requires
a complete native/authored storage tangent and n permitted loss actions,
adds 3*n scalar scratch entries and introduces no warm adapter allocations;
callback-internal behavior remains its author's responsibility. No-D
arithmetic and native refusal semantics are preserved. Structured nonlinear
law refusals survive solver callback boundaries. Ordinary parameter queries
hold the nonlinear-loss law fixed while retaining its full effort chain;
explicit loss-direction methods below additionally supply independent D_p.

Four new focused G0/G1 groups and thirteen existing parametric tests pass
natively: complete pumped state/effort chains; a saturated subharmonic branch
with independent original-ODE grids of 256/512/1024 steps per response period;
generated third/fifth harmonics, separate energy accounts, full-cycle
variational/Liouville checks, signed sidebands, refusal and cancellation/reuse.
The branch needs an actual adaptive harmonic ladder up to order 32 at the
stated tolerances. These checks establish bounded numerical behavior, not
infinite tails, stability, constitutive validity or measured-instrument fidelity.


## Combined adhesive failure, finite caloric release and retained thermal history

`friction_heat::adhesive::snap::combined_loading` regenerates ordered proportional
shear/twist segments and explicit separating-fold release through the original
`JkrCombinedLoading`. Reaching the first boundary and detaching are distinct
operations. Only newly released total contact potential enters the caller's
explicit two-surface/outside partition. Signed prescribed work and the four
reversible stores remain mechanical; previously accumulated release imported
with source history is never heated retroactively. A nonzero outside share
requires a named physical destination.

Original ordered material cards, exact registered caloric objects and receiver
cell, segment, body, material, mass and geometry remain pinned. Mechanical
history, revision, cumulative release ledger and both finite cells publish
together after bounded source and caller cancellation checks. A failed second
receiver, late cancellation, stale history or inadmissible separation publishes
neither owner. Cost is linear in the bounded segment count plus original source
solves and one thermal candidate clone; receiver checks are constant-size.

The reports retain each original accepted/clipped segment and close
`Delta U_contact + Delta H_A + Delta H_B + E_outside - W_prescribed`.
Compensated signed and absolute accounts keep constitutive, telescoping,
partition, ledger-accumulation and each receiver's relative caloric defects
separate. Explicit transfer and total-energy budgets apply independently of the
original source and individual-cell tolerances. No subtraction of large absolute
thermal reference energies is used to recover a small heat increment.

`propose_temperature_continuation` adopts the sealed original loading successor
with unchanged receiver bindings and accumulated heat ledger. Held coordinates,
source step count and prior release remain authoritative, including detached
history after all attached equilibria disappear. Independently transported
thermal time may advance but cannot go backwards. The new revision invalidates
the old stamp. This prescribed source update does not deposit its signed
parameter work or equate source temperatures with finite-cell temperatures.
Adopt the source continuation and its explicit external-work account together.

This same-time, frozen-source impulse composition supplies no loading clock,
conducting link, body recoil, automatic fold detachment, rebonding, partial slip,
or implicit thermal feedback. Eight focused integration groups cover reversible
zero-heat loading, independently reconstructed shared-fold release and partition,
complete balance, once-only transfer, receiver and boundary refusals, atomic
caloric/cancellation failure, and temperature continuation of loaded and released
histories. Formatting/parser and whitespace checks passed. These groups were not
compiled or executed: required constellation siblings and DSR/RCH runners were
absent in the working environment.

## Both finite rigid walls on the original thermal liquid timeline

`free_surface_wall::planar::pair::thermal::FreeSurfaceThermalPlanarWallPair`
retains one original thermal liquid and two finite bodies. Existing pair
mechanics share a sealed liquid view and one closure-driven pilot/replay seam;
isothermal aliases retain their original numerical path. One common thermal
pilot and one final replay use source-derived local viscosity/rheology and
optional buoyancy, original native mass-route enthalpy transport and both
original body/guide/spin/orbital/work audits. Body mass and inertia remain the
authored finite values, and the two per-wall mechanical records own no copy
of the thermal liquid populations.

Authored cell heat and optional held initial-state molecular heat belong to
one accepted liquid interval, never one copy per wall or per trial. Only the
final original thermal report and jointly admitted bodies/liquid escape.
Body elapsed time starts at binding, while a pre-advanced liquid retains its
own physical epoch and increments its accepted clock once. Boundary work is
not inferred caloric heat. Original fixed flooded-plane geometry, explicit
thermal source splitting, material and prepared-boundary limits remain.

Four new pair tests and four original single-wall thermal tests pass. The
new tests cover all-axis cold compatibility, binding an already advanced
thermal state, actual summed body linear/spin/orbital momentum and kinetic-work
balances under heated/buoyant flow, original molecular-source equivalence and
late second-body/source refusal with deterministic retry. Initial rest-fluid
wall reflection may produce no intercell mass routes; continuation requires
positive native routes and every tick retains exact original replay and
extensive-enthalpy checks. Four original isothermal pair tests also pass.

Across this extension's focused runs, 31 unique tests pass (eight new and
23 existing). One untouched partially filled liquid fixture refuses its ordinary
stationary mass-transport step before any changed wall entry; that result is
not claimed green. Tests execute original numerical/material/cancellation
sources with optional conduction/electrical feature registrations excluded.
The normal fs-couple Cargo feature graph and full workspace gate remain
unverified because the async runtime dependency exceeds the environment's
8 GiB memory limit. This is not a whole-fluid total-energy certificate.


### Planar moving terminal controls through finite heat

MovingConductionResponseDrive extends the existing planar thermal endpoint to
direct grouped voltage, total-current, floating and mixed electrodes alongside
original traces. Every thermal proposal refreshes the original bulk/contact
sources and solves the actual native moving electrode problem. Ordered terminal
observations stay in the accepted moving source report; no accepted nonlinear
field is rescaled or independently reconstructed.

Terminal supply uses independently authored native electrode work. Signed
Lorentz output is subtracted from that supply, and only local conductor-frame
Joule/contact heat enters the original finite receivers, once. Complete heat
provider samples/errors, source/phase admission, contact partitions, cumulative
native work limits and atomic continuation remain shared. Current controls
retain motional voltage offsets and permit generation; zero total current may
sustain open-circuit voltage without heat. Focused regressions define two-step
current/mixed electrothermal feedback, trace/voltage parity, floating operation
and budget/cancellation refusals. Formatting, targeted API/accounting review
and independent fixture endpoint arithmetic passed; full repository integration
execution remains outstanding. The fixed-domain prescribed-motion limitations
of the existing planar thermal endpoint remain unchanged.


### Original generalized boundary work in coupled laminar evolution

The driven endpoint now obtains supplied-traction, prescribed-velocity/
Navier-slip wall and dissipative-backflow powers directly from the accepted
native reaction owner. The original boundary quadrature, pressure reaction,
endpoint viscosity, prescribed source trace and convective advector remain
authoritative. Boundary callbacks are no longer resampled for reporting.
The material owner supplies W/m with its original reference density; the
coupler multiplies only by declared depth and actual interval duration.
Mean-zero components explicitly retain their returned zero gauge and
traction components retain their imposed datum.

Both previously unknown work channels are present in native driven reports.
They participate independently in the original adaptive signed and magnitude
defects and accepted-only compensated totals. Numerical upwind loss remains
outside the physical work norm. This observation occurs after private native
thermal acceptance and never alters caloric coefficients, deposited energy
or slip-heat routing. A failed power observation still rejects the complete
candidate or adaptive window without publishing a prefix.

A new prescribed-acceleration regression holds coarse/fine endpoint velocity
and temperature equal while boundary source work alone requests refinement.
It independently verifies accepted wall input equals kinetic gain plus the
positive BDF damping defect, with no new thermal heat. The original decimal
velocity-reversal regression now checks its independent 0.078 J prescribed
input. A finite-slip coupling regression checks moving-wall source work,
original bulk and relative-slip heat, explicit 25/75 percent fluid/wall
partition and physical depth/time conversion; wall work is not deposited a
second time. All six changed coupled Rust sources passed the pinned parser,
formatter and rustfmt --check. Native FEM/thermal compilation and these
integration test executions remain unverified; no numerical runtime pass is
inferred from formatting or static review.


### Nonlinear conduction, native rigid momentum and finite heat

conduction::spatial::mechanical::response closes the original scalar E'/J
material laws and all six rigid midpoint velocities through the existing
bounded Newton-Krylov and fallible finite-difference services. Every trial uses
the same initial native mechanical state and supplied external material
context. Final current is freshly solved; original-cell J cross B loads use the
shared P1 force/wrench integration and native free-floating impulse owners.
Independent midpoint, Lorentz/kinetic and authored-source/heat/storage gates
precede acceptance. Native voltage/current-density traces and direct grouped
voltage/current/floating/mixed electrodes retain their actual controls. All
numerical probes, constitutive trials and final admission consume one native
field-work cap; no nonlinear field is replaced by an affine drag response.

The original exact tetrahedral mass geometry, material card and density claim
pins remain bound. The original density bundle is refreshed at the supplied
external T/phase context and must retain its exact inertial values. Density
response or applicability depending on owned electrical coordinates refuses.
Continuation may therefore update electrical source coordinates while retaining
the original inertia without admitting an unowned field-dependent mass law.

conduction::spatial::response::mechanical additionally closes that native
momentum/current candidate inside the original finite-enthalpy endpoint for any
complete HeatAssembly. Every thermal proposal refreshes the original electrical
sources and restarts from the same mechanical and thermal histories. Original
local Joule heat is deposited once, original finite links/contact/radiation
providers remain in the shared heat residual, and all nested electrical work
consumes the existing thermal circuit-trial allowance. The complete source and
heat samples, candidate source owner and native mechanical result are retained.

The final full energy defect is delta_H + delta_K - terminal_work -
external_mechanical_work - other_external_heat. It is formed by adding the
independent native kinetic-minus-external-minus-Lorentz defect to the existing
thermal defect, so Lorentz transfer cancels once and is not deposited as heat.
The final global and closed-enclosure heat gates are reapplied. Original
caloric, constant-density, routing, typed-source-error and atomic acceptance
requirements remain binding. Refusal or cancellation accepts no input owner.

This is a contact-free electrical whole-body force substep at fixed pose,
geometry and inertia with prescribed B. Thermal contacts are still supported
by their existing heat owners. It adds no pose/advection transport, induction,
gyroscopic evolution or sliding electrical-interface law. Independent nonlinear
bar, full rigid and two-body heat/momentum arithmetic checks passed, as did
pinned formatting and bounded exact-production typechecks with dependency
stubs. Six native and two complete-heat integration regressions were authored;
the full dependency-backed repository tests remain unexecuted.

## Scenario reservoir controls with native flow and heat

`scenario_open_boundary::thermal::{step_serial, step_forced_serial}` samples the
existing typed velocity/pressure controls into the original joint momentum/D3Q7
thermal transaction. Uniform force, retained spatial force and authored endpoint
force use their native physical velocity epochs. Scalar relaxation and full tensor
diffusion share this route. Selected spatial profiles still refuse uniform updates;
unselected profiles remain active. No preliminary control rebind is published.

Targets, populations, force buffers, heat/flux histories, finite walls, receipts
and counters publish only after the closing thermal gate. Endpoint force buffers
remain unchanged on refusal. Heat inputs retain the original thermal carrier's
explicit lattice/reference-capacity units; momentum material scaling does not
infer thermal data, mechanical heating, buoyancy or compressible enthalpy.

Focused regressions cover D3Q19/D3Q27 original-owner reduction, scalar/tensor and
held/changing-force paths, independent momentum/heat balances, late heat refusal,
retry and profile ownership. Pinned Rust parsing and scoped diff checks passed;
Cargo typechecking and executable regressions were unavailable in this session.

## Prepared nonlinear contacts with complete live-support heat

`MaterialHeatNetwork::prepare_conductivity` and its prepared
`field_motion_support_conductivity_candidate` join original E/J material and
V/I/J interface response to the same finite-enthalpy/live-support transaction.
`ContactConductionResponse` declares coordinates and source orientation; the
existing electrical network retains interface cards, geometry, source pins and
A/B heat partitions. The prepared seed API supplies original circuit inputs.

Solved-field bundles refresh scalar or directional Fourier conductance before
the original thermal-contact owner applies its series resistance once. The
same-trial closed-enclosure radiation result is retained once and remains
internal energy exchange. Native electrical iterations consume the shared work
budget; source refusal, cancellation and failed admission publish no histories.
The existing static preparation and support entrypoints retain their reduction.

Auxiliary thermal-contact/emissivity claims must be electrically independent in
value and applicability. Integrated-layer and externally supplied additional
link owners explicitly refuse this frontend. Original fixed geometry/inertia,
constant field-storage and support/caloric restrictions remain. Focused tests
cover independent series/Fourier/radiation and per-cell/total energy laws,
constant-response reduction, late source failure, budget/cancellation and retry.
Pinned Rust parsing and scoped diff checks passed; Cargo compilation and these
runtime regressions were not executed in this implementation session.

### Physical heat and wave-flux resonance queries

`fs-couple::phs_power_resonance::PhsPowerQuery` binds an original `HarmonicModel`
or a fully feedback-closed `ScatteringModel` to the shared `fs-modal` signed-power
query owner. It resolves actual extrema, half-power crossings, linewidth and
observed Q, including their local physical-parameter derivatives and bounded
grid surveys. Channels are independently evaluated original R/T dissipation,
aggregate conjugate source work, individual original-port work and selected
incident/outgoing wave flux. Boundary reference loading never becomes internal
material heat. Equilibrium linearizations retain signed incremental-work
semantics. Negative work is not made positive or squared to manufacture peaks.

The source law lends reusable buffers with complete frequency and parameter
jets in original input order. Selected wave slots carry incident output waves;
ordinary slots remain conjugate inputs. Original J/R/Q/G/T and positive-reference
directions enter the native physical equations. References are constant in
frequency; physical dispersion remains in retained states. The first-order and
full wave observations share the same normalized paired-frequency contraction,
so near-stationary slopes reproduce the same original sample.

Every frequency uses one original factorization. Refinement requests two state
solves; only final extrema/edges request three-order curvature or five-order
parameter responses. Native memory, work and residual admission remains active,
and shared query controls bound all callbacks and extrema. The actual P_peak/2
crossings define a positive power peak's linewidth, without a Lorentzian fit.
Q is peak frequency divided by measured width, not a free-pole decay rate.
Sources, models and channels must agree across peak and bandwidth queries.
Grid resolution, local feature isolation and constitutive validity remain
caller obligations. Cancellation, missing/flat crossings, inconsistent samples
and source/physical refusals publish no partial query result.

### Complete material-memory power and receiver resonance queries

`PreparedStructuralInteriorRelaxation` now resolves original signed-power
extrema and coherent field peaks/notches, their independently bracketed
half-power crossings, linewidth/Q, local physical-parameter derivatives and
grid-resolved surveys. All queries reuse the shared `fs-modal` numerical owner;
no reduced pole fit, sampled damping, invented loss floor or finite model
perturbation replaces the original coupled response.

Power selections retain structural/acoustic source work, reciprocal interface
exchange, structural and wall heat, material-memory heat and separately owned
structural/acoustic work into each original Maxwell branch. Total dissipation
adds structural, wall and relaxation heat once. Negative source/interface work
keeps its sign; heat is observed from the original constitutive channels.
The first-frequency response now carries independent first-order power and
per-arm work/heat accounts. The normalized quadratic-work owner excludes gyro
exchange and retains original weak and subnormal coefficients.

Receiver selections cover original displacement, velocity, pressure potential,
physical pressure, fluid force and displaced flow. Complex receiver weights
act as coherent linear rows without conjugation. Pressure retains every
kinematic derivative of `-i*omega*psi`. A pressure-point receiver may supply its
original basis evaluation row. Authored frequency/parameter receiver jets
include changing geometry or observation laws explicitly; no such derivative
is inferred from a label. Simple cancellation notches retain their frequency
and total complex-amplitude derivatives, with absent phase/log-amplitude at
an exact zero.

Source and receiver spectrum traits lend reusable buffers. Frequency requests
supply complete values/slopes; final parameter requests supply all five local
orders. One native factor per frequency serves two state solves during
refinement, three for a final power curvature, or five for final physical
parameter jets. Original material/geometry directions are prepared once and
remain bound to the same immutable preparation. Full response jets are not
computed at every refinement trial. Provider domain/work limits remain with
the authored laws; cancellation is checked around provider calls and native
solves. Existing complete-state and adapter admission remains in force.

Power widths cross `P_peak/2`; receiver widths cross `|H_peak|^2/2`. Both report
actual measured peak frequency divided by width, rather than equating a driven
response Q with pole decay or a material loss factor. Source, receiver, channel
and model must match the supplied peak. Final jets must reproduce the same
first-order sample; inconsistent observations refuse. Strictly local feature
isolation, positive cavity frequencies, grid resolution, fixed coordinates and
topology remain caller obligations. Cancellation, missing/flat crossings,
capacity failure or physical refusal publishes no partial query. Native tests
execute original coupled recovery and query physics; they do not establish
complete outer FE mesh-construction or full-workspace Cargo validation.

## Finite thermal wall motion under prescribed gas loading

The original single-body and opposing-body thermal wall adapters expose
`step_candidate_with_sources`, consuming the lower owner's borrowed
`FreeSurfaceWallSources`. Authored heat, optional original molecular heating
and prescribed cell/link gas pressure share the same accepted liquid history.
The existing heat-only and molecular methods delegate to this path.

The original mechanics are unchanged: one pilot, the existing mass-whitened
midpoint solve for each body, one final native replay, and independent body,
guide, spin/orbital momentum and kinetic-work checks. Both bodies share one
liquid owner. Only the final replay's thermal receipt, molecular samples and
clock advance are retained. The caller must hold pressure values fixed for the
same geometry across both trials; arbitrary interior-mutable callback behavior
is not certified by the function type or affine wall-response check.

Gas pressure acts through the original liquid/gas populations and reaches the
finite walls through subsequent native evolution. It is not an instantaneous
extra body force, and neither pressure work nor wall work is deposited as
caloric heat. The fixed pressure reference, source applicability, low-Mach
bounds, flooded planar geometry, first-order thermal splitting and absence of
a gas storage/pressure-volume energy model remain explicit.

Four new original-source tests pass: reference-pressure identity for both body
adapters with actual gas links; five-tick pressure-driven recoil of both
initially resting bodies with independently reduced momentum and work; single
body molecular-source equivalence; and final-replay pressure/cancellation
refusal with exact retry. All eight existing single/pair thermal-body tests
also pass. Across the related lower pressure and conversion repair checks,
25 unique focused tests pass (10 new, 15 existing), with no failures. This
executes the original numerical kernels and does not claim a full Cargo build
or the optional solid-conduction/electrical enthalpy-coupling graph.


## Event-timed combined adhesion, footprint evolution and separation

`friction_heat::adhesive::snap::combined_loading::timed` consumes the original
sealed affine loading ramp with exactly matching accepted source history and
thermal/heat clocks. The explicit boundary policy either stops at the first
boundary or releases an admitted separating radial fold and follows the
remaining prescription detached. Opening and nonseparating boundaries retain
the stopped prefix and expose untraversed time. Original source steps are
regenerated and compared in full; geometry samples never advance accepted
mechanical history or manufacture release.

The same native enthalpy timeline owns conduction, ambient/reservoir supplies,
births and scheduled events. Events at each endpoint execute before footprint
retargeting or release heat. Original released potential is partitioned once;
its exact owned patch is then removed, restoring actual ambient exposure before
any detached-tail transport. The footprint retains its link identity after
separation. Original ordered receivers, caloric laws, cell geometry/mass and
exact expected patch/absence are checked before native fluxes and after every
individual event and step. Scheduled takeover or revival of that link refuses.
Unrelated thermal links retain their own ownership.

Attached coverage and conductance use explicit bounded left-endpoint footprints.
Source area samples retarget actual native coverage at interval endpoints;
conductance uses the fixed original reference G/A, preserving material/geometry
resistance without cumulative rescaling drift. The complete source radius
envelope is conservatively reserved against both physical cell faces and
competing coverage, including the maximum at a reversal through zero load.
One native contact sweep checks both receiving surfaces with bounded cancellation.
Reversals split at a represented, separately gated zero-load clock. Internal
footprint cuts may round downward to respect their duration ceiling; physical
source/boundary clocks are never moved.

One shared native substep allowance covers every attached interval and detached
tail. Independent footprint-count/duration and original mechanical-stage budgets
bound work. Signed and absolute accounts combine each original mechanical/heat
balance with native transport defects, counting contact heat once and retaining
ordinary thermal supply and incoming material energy separately. Loading history,
caloric state, clock, topology and release ledger publish together after the final
checkpoint. Source, geometry, caloric, topology, budget or cancellation refusal
leaves all accepted owners unchanged.

This is a first-order footprint approximation requiring independent footprint
and thermal timestep refinement, with frozen constitutive properties and
conductance per area. It supplies no exact swept-area integral, automatic
source-temperature feedback, temperature-dependent conductivity update, inertia,
recoil or rebonding. Five focused integration groups cover independently computed
post-event cooling and held/reversing conduction, footprint/clock ownership,
shared budgets, scheduled contact conflicts and atomic retry. Two additional
clock tests cover large origins and unchanged physical event times. Formatting
and whitespace checks passed; these tests were authored but not compiled or
executed because constellation siblings and DSR/RCH runners remain unavailable.



## Implicit laminar enthalpy with external RCR storage

With `laminar-enthalpy`, `LaminarEnthalpy::advance_windkessel` resolves the
original nonlinear BDM momentum/rheology and RCR/Windkessel ports, including
optional series inertance, inside each native backward-Euler caloric residual.
The original material is queried at the actual trial temperature and shear;
its accepted source snapshot remains owned by the returned flow. Original
viscosity, convection, pressure response, final momentum equations and both
unreduced circuit laws retain their native admission and residual gates.
Entirely circuit-connected components require no artificial traction datum.

Every thermal/Newton probe starts from the SAME accepted velocity and compliance
pressures. Previous port Q is integrated from that native velocity, never
supplied independently or advanced by a provisional thermal trial. Compliance
coordinates use the same reservoir references and physical clock; geometry,
port ordering, density and C/L/R remain fixed, and reservoir increments are held
through the interval. BDF1 is used consistently for momentum, external storage
and caloric transport. This entry adds no adaptive circuit window, moving mesh,
variable density or Boussinesq/load combination.

The heat callback borrows the complete native circuit/flow result. Exact
sign-split BDM rates drive the original enthalpy donor selection, with explicit
external inlet caloric states required on every incoming branch, including
reversal. Bulk strain and explicitly routed relative-slip heat are deposited
once. External resistor losses, compliance/inertance storage, reservoir/port
work and numerical BDF defects remain separate; no thermal receiver is inferred
for an external load. Original additional heat links and sources remain usable.

Acceptance moves the original resolved field into the standard coupled
candidate and native pressure/Q/storage/power receipts into its typed exchange
payload without cloning or reconstructing the field. The receipt exposes
accepted compliance coordinates for continuation with that same candidate's
velocity and thermal clock. Native energies/powers remain J/m and W/m; declared
depth and actual interval convert them to joules. Port power covers added
pressure increments relative to original traction, not an inferred full-system
energy or passivity certificate. Source, stream, numerical-budget, cancellation
or final thermal refusal publishes none of the three physical histories.

Focused native regression definitions cover independently eliminated two-port
RCR/inertance momentum and storage equations, signed reversal with two-cell
caloric transport, external energy versus BDF damping, original Arrhenius/
Carreau viscosity with a manufactured finite-slip thermal endpoint, and late
source/stream/budget/cancellation rollback. All four changed Rust files passed
the pinned parser, formatter and `rustfmt --check`. The full native FEM/thermal
dependency build and these integration test executions remain unverified in
the current partial restored environment; no numerical pass is inferred from
formatting or static review.


## Autonomous nonlinear-loss parameter responses and physical work adjoints

`DissipativeAutonomousHamiltonian::prepare_sensitivity` reuses the original
full-state bordered autonomous factor, fixed waveform phase reference, and both
AFT grids. Every storage, feedback, material-memory, and nonlinear-loss tangent
remains in the original physical residual. It returns phase-aligned Fourier
waveform and natural-frequency parameter responses, arbitrary cycle-mean
objectives and adjoints, and separate supplied, native-resistance, and nonlinear
internal-loss power or work over the solved period.

`DissipativeAutonomousParameterDirection` combines the existing physical
J/R/G, storage, and TOTAL feedback partials with an optional original D_p at fixed
state and effort. `None` explicitly holds that constitutive law fixed while
retaining D_e e_p. The shared fs-phs owner keeps direct loss and effort-chain
terms separate until the final compensated physical field contraction.
Nonlinear heat differentiates as e_p dot D + e dot (D_e e_p + D_p); its state
covector includes the full storage Hessian, D_x, and D_e. Cycle work retains the
solved-frequency contribution exactly once through the bordered adjoint.

Preparation reserves 6n + 4m additional query scalars and retains the original
source, callback, factor, and two-grid resource limits. Failed or incomplete
derivative callbacks, cancellation, unresolved finite-grid evidence, singular
cycles, and excess resources refuse explicitly and leave reusable prepared
factors intact. Finite-order and finite-grid admission does not certify an
infinite harmonic tail or orbital stability.

Three new independent-cycle tests cover a live loss direction that changes
frequency and waveform, simultaneous J/R/G/storage/feedback/loss perturbations,
a fixed D law with nonzero D_e e_p, all six physical power/work objectives,
forward-adjoint agreement, callback refusal, cancellation, and scratch budgets.
All 18 relevant autonomous tests pass, including the existing dissipative
forward and ordinary autonomous sensitivity, power, and adjoint regressions.


## Pumped nonlinear-loss parameter responses and physical work adjoints

`ParametricHamiltonian::prepare_sensitivity` accepts the attached original
nonlinear-loss source and factors its complete state/effort tangent.
`solve_parameter_with_dissipation` and `gradient_parameter_with_dissipation`
add an independent D_p callback at fixed state and effort. Existing parameter
methods explicitly hold that constitutive law fixed and still retain D_e e_p.
An explicit D_p without an attached nominal loss source refuses.

Original J/R/G, storage, input, pump, response-period, and response-clock
parameter chains remain in the residual. The shared fs-phs parameter owner keeps
D_e e_p and D_p separate until the final compensated physical contraction.
`MeanNonlinearPower` and `NonlinearEnergyPerCycle` independently measure
e dot D. Their state covectors include D_x, D_e, and the complete H0_xx + Q
storage tangent; direct heat is e_p dot D + e dot (D_e e_p + D_p). Cycle work
retains the response-period derivative, including multi-source-period views.
All eight pump/port/native-resistance/nonlinear-loss mean and cycle objectives
use this same complete physical parameter direction.

Preparation reserves 4n additional query scalars when nonlinear loss is
attached. Nonlinear-work observation uses 2n reusable scratch entries within
that bound. Parameter adjoint queries reuse the original factor and cached
coarse/fine trajectories without additional Hessian evaluations or orbit
solves. Other direct power accounts avoid unnecessary loss-parameter calls.
Original source admission, callback completeness, cancellation, finite-grid
checks, and resource limits remain explicit.

Three new tests compare combined physical/loss/clock waveform and all eight
objective gradients with neighboring original harmonic models and independently
written RK4 equations, cover held-fixed constitutive loss with a live effort
chain on two source periods, and exercise loss-only directions, missing nominal
sources, incomplete callbacks, cancellation, and reuse. All 20 relevant
parametric tests pass. The existing converged order-32 subharmonic branch now
requests 256 sensitivity samples to satisfy the native minimum of 132; physical
tolerances and other limits are unchanged. These finite-order/grid checks do
not certify infinite tails, orbital stability, or nonsmooth constitutive laws.


## Driven friction events in adaptive thermal continuation

The existing Hertz full/two-half thermal comparison checks whether the coarse
and fine driven histories both contain an arrest, including reversals that leave
both endpoints sliding. A disagreement triggers the existing bounded refinement
or refusal. Available inner travel uncertainty also contributes
`|W^-1 a| * travel_error` to the supplied-mechanical-work comparison; it cannot
vanish merely because the coarse and fine nominal work agree.

Source-resolved affine invariant stopping, static hold and reversal retain the
original native event receipts. Cooling continues after frictional heat ceases,
and changed source temperature can later release a held contact through the
original static surface. An analytic event clock is conditional on the held
source parameters of that phase; source evolution still requires thermal-window
refinement. Actual friction heat, external mechanical work and external heat
remain separate, and accepted bodies, mode and thermal state publish together.

Focused physical regressions compare independent signed exponential motion and
traction-power integrals with the original finite-cell thermal roots. They also
cover reversal-history disagreement with matching endpoint modes, original
Stribeck travel uncertainty in supplied work, stopping followed by thermal
release under the same motor, and late source/cancellation rollback. These are
local numerical comparisons and implementation checks, not globally certified
trajectory bounds or experimental validation.



### Nonlinear moving planar conductors on the native circuit and heat timeline

ConductionResponseCircuit::new_with_motion and step_with_motion_candidate use
the original RLC descriptor, authored charge/flux history and original E'/J
material laws. Each frozen conductor enters through its exact resistance and
motional current offset. Shared internal Norton ports preserve non-ground
references and conserved offsets; accepted owners retain original authored
graph topology and ports. Every moving constitutive trial reserves two native
field solves and one descriptor solve, including discarded/refinement trials.
Actual-voltage recovery performs no third sparse solve.

Physical work closes independently as circuit storage change plus graph heat,
conductor-frame field heat and signed Lorentz output equals authored port work.
The descriptor's homogeneous loss and internal offset-port work are algebraic
accounts, never additional material heat. The complete HeatAssembly endpoint
supports fixed C/L and original source R(T)/C(T)/L(T), preserving reversible
field entropy transfer, nonlinear contact sources, local heat receivers and
complete typed provider samples/errors. Its finite-heat supply is authored
port work minus Lorentz output. All thermal/constitutive probes start from one
accepted electrical and thermal history; a refusal publishes none of them.

Focused regressions cover nonlinear RC/RLC motoring, generation into capacitor
storage, zero-motion parity, continuation, complete contact/radiation heat,
source-dependent storage, typed refusal, budgets and cancellation. Pinned
formatting, bounded exact-source type checks and independent midpoint/work/
heat reference arithmetic passed. Full dependency-backed compilation and
integration execution remain outstanding. Geometry, velocity and imposed B
are prescribed; induction and mechanical history evolution are outside this
entrypoint's scope.

## Reciprocal heat in finite-plastic sheet forming

`solid_enthalpy::finite_plastic::sheet::FiniteThermoplasticSheetCoupling` binds
the original P1/Q1 thermoplastic sheet to one active caloric cell per element.
It uses the existing source-admitted `FiniteThermoplasticMaterial`: mechanical
and caloric laws retain the same material card, original caloric object,
reference density and positive total heat-capacity bound. Reference area times
original thickness fixes cell volume and mass. Admission checks every accepted
quadrature temperature against its cell, including after a nonuniform native
sheet path; a diagnostic mean cannot stand in for the retained point state.

The sheet and tetrahedral couplings share one bounded feedback algorithm and
thermal observer adapter. Every proposal starts from the same accepted native
mechanical history and thermal timeline. Original quadrature volumes integrate
plastic and creep dissipation plus signed reversible heat exactly once. The
thermal owner retains its substeps, source transitions, conductance observers,
events and conservative exchanges. Its actual endpoint temperatures then drive
the mechanical endpoint; relaxation operates on caloric enthalpy. An unrelaxed
iteration reuses the preceding actual-temperature mechanical proposal.

Acceptance separately checks temperature agreement, the sum of absolute cell
heat discrepancies, caloric balance and combined mechanical/caloric balance.
Hardening storage, rate potential and signed endpoint work defects remain
separate from physical heat. The returned candidate contains the complete
sheet, thermal simulation, observer history and native equilibrium/energy
reports. Original cell identities, reference geometry, mass, accepted enthalpy
and clocks remain binding; rejected iterations change none of these owners.

This is quasistatic membrane mechanics with element-uniform caloric cells and
independently specified reference-network conductances. Evolving thickness is
mechanical output, not an implicit update to Fourier or convection geometry.
The existing finite-plastic enthalpy regression seam adds native heated sheet
forming and self-heating with independent physical-loss and storage checks.
Native Rust execution of the new coupling remains pending available build
resources.

## Source-bound thermal memory on the original shear shell

MaterialShellImage::with_thermoelastic_memory selects explicit instantaneous
expansion, constant-stress specific heat and normal-conductivity roles from
the SAME held regional material bundles. Positive absolute temperature must
be explicitly typed. The shared coupled slab law derives fixed-strain heat
capacity from the original principal plane-stress Q, expansion and density;
every selected source band is retained. Orthotropic21 exchanges the normal
expansion entries without reversing in-plane engineering shear. Coatings,
asymmetric/effective sections, other substrate memory and frequency-anchored
loss refuse; independent original Kelvin and interface loss is retained once.

ThermoelasticShellImage keeps the original isothermal eigenbasis, physical
mode signs and mass normalization, force footprints, source receipts and both
exterior/cavity volume maps. The original MITC thermal-factor owner projects
all facets once, including the elastic director-bubble map and seven-point
quadrature. Compensated covector rotation precedes compensated physical
strain contraction. The existing power-scaled QR owner compresses complete
rows without an eigenvalue cutoff. No dense thermal nodal matrix or second
mechanical eigensolve is formed.

Under insulated faces and absent lateral conduction, the uniform thickness
temperature mode adds reversible adiabatic membrane stiffness, with no heat.
Odd diffusion modes add positive bending memory, stiffness dispersion and
physical dissipation. Every retained cross-mode term remains. Exactly
identical represented regional spectra share one compression; different
thermal rates retain separate histories. Each regional spectrum retains its
own finite-band response and truncation-error report.

prepare_coupled uses the original complete prepared held-force owner at zero
mechanical/thermal perturbation. Original bulk/interface dashpots, conservative
membrane columns and bending memory share its work, storage, heat and clock.
Interface elastic stiffness remains in the original pencil and is never
reapplied. Membrane reference stiffness is only a coordinate factorization
scale; its reported connection forces are virtual thermal-coordinate forces,
and the original force ceiling applies to that scale. Pressure transfers
remain zero until an explicit acoustic owner is attached. fill_force uses
the unchanged physical port maps without allocation and refuses atomically.

Checked projection, factor storage, memory and original runtime ceilings
remain binding. A single conservative instantaneous-stiffness bound must fit
the intersection of ALL regional and interface applicability bands. The
final augmented exponential extent is checked before avoidable preparation.
Thermal diffusion rates impose no artificial sample-period/relaxation-time
restriction. The local homogeneous plane-stress slab approximation excludes
lateral/boundary heat flow, finite strain, thermal initial fields and evolving
material temperature; static elastic director condensation remains the spatial
approximation. Constitutive bounds do not certify modal/spatial truncation,
projection roundoff, arbitrary forcing or later external couplings.

Source and runtime regressions cover typed held data, signed anisotropic axes,
original curved-shell strain/Gram integration, complete complex stiffness and
separate thermal/Kelvin/interface heat, pulse decay, null expansion, regional
thermal-rate ownership, and transactional resource/band/cancellation refusal.
Independent constitutive and exact-state numerical oracles were checked;
native Rust execution remains unverified in this environment.

## Executable certified damping blocks

The explicit plate-viscosity certificate consumer binds the existing
retained-matrix free-pole decisions to the original modal time owner.
ModalPlateViscosity::certify_for_time preserves the original energy-quotient
singleton owner when available; WindowedPlateViscosity::certify_for_time uses
the existing near entries and far-row bounds. Only complete required Full
principal blocks are materialized. Their original passive factorization
retains every signed coefficient without a new fitting or truncation rule.

CertifiedPlateViscosity::prepare_coupled keeps the original coordinate order,
frequencies, state, clock, pressure transfers and admission budgets.
Singleton physical rates come directly from the certified pole centers.
Full clusters carry their entire intrinsic loss through coupled channels
and receive no parallel scalar damping. Ordinary preparation uses the
existing exact islands/scalar fast path; Prepared uses the existing bounded
dense owner. HarmonicOnly, prior scalar damping and mismatched frequencies
refuse. Frequency equality alone cannot identify a rotated degenerate
eigenspace: the caller must supply the original basis.

The opt-in approximation removes coupling between certified blocks. The
certificate bounds free-decay poles of the retained source operator, not
arbitrary forcing, pressure or heat error, nor poles after adding another
support, contact, radiation or memory coupling. Runtime work and nonnegative
physical heat belong to the realized operator; omitted coupling is never
recorded as dissipation. Compact image storage is linear in mode count plus
the sum of squared Full-block sizes; windowed compilation creates no global
dense damping matrix. Original FE and runtime capacity policies remain.

Four authored regressions cover independent rotated scalar dynamics/heat,
a complete 96-mode dense/windowed source requiring only one coupled pair,
original singleton rates and exact lossless modes, plus binding/resource/
cancellation refusal and retry. Independent spectral and exact-state oracles
were checked; native Rust execution remains unverified in this environment.

## Oscillating axial wall ports on the original channel field

`ChannelContinuousOperator::solve_harmonic_with_walls` solves the original
tridiagonal continuous momentum equations for pressure and independent
axial wall peak phasors under `exp(-i*omega*t)`. It retains every velocity
cell in the planar, circular or annular gas/material channel. The solve is
linear in cell count, uses two bounded work vectors and introduces no
dense matrix, time step, fitted impedance or numerical heat.

The response retains pressure work and both mechanical wall powers
separately. Wall forces act on the fluid; their physical opposite acts on
the mechanical wall. Complex wall power is `F_on_fluid*conj(U)/2`.
Independently accumulated relative-shear heat and kinetic inventory obey
`P_pressure + sum(P_wall) = heat - i*2*omega*mean(E)`.
A source may recover real work. Raw real/reactive closure residuals remain
visible, with roundoff scaled by the actual uncancelled cell and wall work.
`observe_harmonic_with_walls` checks a coupled solver's supplied full field
against those same original equations and accounts.

A declared work bound is checked before allocation, and cancellation is
polled throughout elimination and observation. No accepted time-domain
field or source changes. The circular axis cannot supply a nonzero wall
velocity. Frequency must be positive for peak-phasor cycle means, and the
original nonlinear kinetic-head terminal refuses this linear operator.

G1/G3 `channel_wall_harmonic` checks oscillatory planar Couette refinement,
pressure-only agreement with the original observer, wall-force/work
balance, phase and work recovery, geometry-independent ownership and
invalid/work-budget/cancellation refusal. These are fixed-geometry,
incompressible Newtonian laminar equations; mesh refinement remains the
caller's accuracy control. Normal wall motion, rotation, thermal evolution
and transition prediction are not inferred.

## Material-source resolved laminar channels

`slit_channel::material::MaterialChannel<M: FluidTransport>` admits original
liquid, gas or explicitly raw transport owners into the shared planar,
circular and annular no-slip momentum solver. The source is evaluated once
at construction and retained with its receipts and conditions. Geometry,
cell inertia and face shear come from the same constructors as the gas
owners; no synthetic gas or independently fitted resistance is introduced.
Ordinary stepping does not resample the source or require it to be cloneable.

The original transverse field supplies pressure condensation, harmonic
response, wall reactions, signed mechanical work and physical/numerical
energy accounts. `with_source_controlled` is an immutable, cancellable
proposal at fixed geometry and exactly unchanged density. It preserves
every velocity, kinetic inventory and historical exchange, while the new
viscosity changes future diffusion. Material errors remain typed; unequal
density refuses because this owner cannot supply a changing-mass balance
or density parameter work. Flow previews clone only when the source is
cloneable; a borrowed source can preserve the original source identity.

G1/G3 `channel_material` covers original-gas parity, liquid Poiseuille
behavior, retained source identity and sampling count, non-cloneable
sources, viscosity continuation and source/density/cancellation refusals.
The input describes a spatially uniform Newtonian transport snapshot.
Temperature, species, phase and thermal evolution, non-Newtonian shear
dependence, compressibility and transition prediction remain with their
respective physical owners.

## Resolved axial laminar flow in concentric annuli

`slit_channel::annular::AnnularChannel` retains a coherent source gas and
the original shared finite-volume momentum field between two material walls.
Each cell uses its actual annular area. Face shear uses integrated cylindrical
resistance, `2*pi*mu*L / ln(r_right/r_left)`, including the half-cell wall
intervals. Gap-product areas and compensated logarithms retain narrow-gap
resolution without subtracting nearly equal squared radii. No hydraulic
diameter or fitted annular drag replaces the field.

Both prescribed wall velocities are axial, in inner/outer order. The existing
pressure port, implicit kinetic-head terminal, signed wall work, impulses on
the fluid, relative-shear heat, kinetic storage and numerical dissipation
remain on the shared momentum owner. The factorization and hot-step buffers
are reused. Flow previews and gas changes publish only after acceptance;
gas changes preserve each cell's momentum and expose their signed parameter
work through the original material-transfer receipt.

G1/G3 `annular_channel` covers the independent logarithmic Couette solution,
wall traction and heat, annular Poiseuille refinement, thin-gap behavior,
startup/reversal, pressure condensation and material/refusal behavior.
This is fixed concentric geometry, incompressible fully developed Newtonian
flow. It does not supply cylinder rotation, eccentricity, normal wall motion,
entrance flow, thermal evolution or transition prediction. A zero inner
radius must use the existing circular symmetry owner.

## Native rigid JKR motion and event-timed snap heat

`friction_heat::adhesive::snap::native` composes the original loaded finite-body
JKR phase with the original finite-cell thermal timeline. The heat owner retains
the accepted full source/mass/policy/time/step token, so later native phases
continue their actual event chronology. Changed mobility, stale heat stamps or
inconsistent refreshed marker geometry/rates refuse without resetting history.
Only sealed original snap release enters the declared surface/outside shares;
native conduction and authored thermal events run at their original times.

Actual body endpoints, source history and caloric state publish as one candidate.
The combined balance eliminates internal snap transfer and retains actual body-
wrench work, thermal boundary/source heat and incoming material energy. A separate
sum-of-absolute-defects gate prevents mechanical and thermal errors cancelling.
The native phase keeps its frictionless normal, held-pose/Jacobian assumptions;
the surrounding pose owner consumes travel before refreshing geometry. This
does not add tangential traction or temperature-dependent normal feedback.
The original timed-heat target includes loaded multi-phase attachment/pull-off,
native history mismatches and rollback after late caloric/cancellation refusal.

## Native sparse solid material heat into finite enthalpy cells

`tensor_material::solid::thermal::spatial::step_thermal_candidate` advances the
original fixed-duration sparse Kelvin/Maxwell trajectory and deposits its
independently integrated element losses exactly once into complete reference-
mass/volume caloric receivers. Original inertia, quadrature histories, moving
support and imposed-eigenstrain work survive; work and solve error never become
heat. The existing enthalpy owner retains all source/transport/event hooks and
may subdivide the interval. Its thermal power is the mechanical interval mean;
accuracy requires refining both intervals. Geometry/material coefficients stay
frozen within explicitly admitted fully-solid temperature bounds. Mechanical,
thermal and observer state publish atomically after their separate energy gates.
The source API and four focused integration regressions are syntax/format
checked; full Cargo execution remains blocked by unavailable sibling dependencies.


## Native duct outlet harmonics and coherent aperture radiation

`acoustic_realize::duct_harmonic::observe_duct_harmonics` contracts the original
bell and active tone-hole `AcousticOutletFlow` rows directly against every full
`HbOrbit` coordinate. These internal output rows need not be external power
ports. Loaded-terminal inertia, material memory, native outward signs, sparse
tap IDs, stations and effective open radii remain owned by the original assembly.
DD contraction preserves complex phase; DC flow stays separate and positive
lines use peak `2 conj(c_h)` amplitudes. No inlet-flow substitution is made.

Explicit outward normals and mouth-relative receiver geometry define local
aperture frames. Existing circular quadrature and harmonic Rayleigh transfer
produce each mouth's pressure and velocity; velocity rotates back to one common
Cartesian frame before DD coherent summation. The returned existing
`PreparedHarmonicField<3>` feeds the established harmonic pressure renderer.
State, opening, harmonic, quadrature, storage and cancellation limits apply.

The exterior is explicitly supplied, homogeneous and lossless; each prescribed
uniform circular aperture has its own oriented infinite-baffle observation.
This adds no radiation load, mutual aperture scattering or atmospheric model,
and does not establish candidate-orbit convergence. Five physical regressions
cover the actual loaded bell versus inlet, finite-aperture rendered pressure,
sparse tap identity and phase cancellation, rigid bell/lateral-hole rotation,
closed-mouth silence and atomic admission/cancellation refusals.

## Native continuous tangents in periodic and variational solves

The original pHS orbit adapter discovers `Storage::continuous_tangent()`
after its existing exact-quadratic cached path. Ordinary harmonic balance,
variational shooting and Floquet adapters can therefore differentiate assembled,
boxed nonlinear storage without a separately maintained constitutive callback.
Native SOS, direct sums, shared storage and complete Maxwell memory retain their
original effort law, geometric prestress curvature and memory cross terms.
Explicit caller-provided actions remain available as overrides.

A selected native action must fill every physical derivative row with finite
values. Refusal or incomplete output stops the analytic solve; no hidden
finite-difference columns replace a declared action. Unsupported storage retains
its existing solver support and limitations. Feedback derivatives and harmonic
resolution remain the responsibility of their existing owners. Focused tests
exercise ordinary phased nonlinear HB, original cycle work, variational shooting
against the independent Liouville determinant, and strict native refusal.

## Physical harmonics of original nonlinear port outputs

`phs_orbit::spectrum::observe_port_spectrum` evaluates `G^T grad H(x)` from
every original HB state coordinate, including explicit material memory. Output
order is independent of state order: nonlinear storage can create DC, even
harmonics and higher harmonics absent from the state coefficients. The shared
FFT owner synthesizes the full state and projects the physical output, with
one original-gradient evaluation per sample and O((states + ports) N log N)
transform work. No small-signal linearization, added harmonic envelope or
constitutive replacement is used.

Three or more successively doubled, shifted sample grids preserve the original
phase origin. Two successive comparisons must fit each port's absolute RMS
budget, shared between retained-coefficient change and observed omitted-band
RMS. Explicit source-evaluation, sample and numeric-storage bounds apply;
cancellation or refusal publishes no partial spectrum. DC remains a separate
real physical load. Positive output lines use peak `exp(-i n omega t)` phasors,
`2 conj(c_n)`, compatible with existing mechanical and acoustic consumers.

G1/G3/G4 cases in `tests/port_spectrum.rs` check analytic biased cubic output,
phase and amplitude conversion, weak physical scales, unresolved output order,
DC, budgets and cancellation. An integration case propagates the actual force
harmonics through the existing compact acoustic observer and pressure renderer
against an independent retarded-time force/force-rate oracle. Finite-grid
agreement is not an infinite-tail bound, equation-solve or stability certificate,
continuum validation, reciprocal acoustic load, or measured-instrument claim.

## Molecular absorption on original finite gas storage

`molecular_heating::MolecularGasHeating` binds an original prepared molecular BEM
volume-absorption operator to the original `IsochoricGasParcel` and
`RelaxationAbsorption`. The complete re-evaluated gas/reservoir snapshot must
match, and positive quadrature must represent the parcel's actual volume to
roundoff. The immutable parcel retains its full thermochemical model, source
packs, composition and density; the acoustic snapshot alone does not certify
historical pack identity, so the caller must supply the original caloric source.
The lumped parcel has no independent shape; geometric coverage and cubature
convergence remain the quadrature author's responsibility.

`heat_load` retains original source order, frequency and receiver geometry while
evaluating coherent, independent-source or coherence-factor molecular heat.
It multiplies actual integrated watts by an explicit physical duration exactly
once. No solid heat capacity, extra RMS factor, wall work, outgoing flux or
second attenuation term is injected. Queries reuse the existing bounded BEM
readout without a boundary re-solve or dense source covariance.

The returned source-bound load can apply `heat_gas` through the original
nonlinear NASA-9 internal-energy inversion on a clone. Fixed-volume storage
uses Cv and preserves mass/density; the same source supplies EOS pressure and
transport. The candidate's next molecular state is re-evaluated through the
original relaxation law. Native caloric refusal, endpoint source-domain
failure, a caller-admitted maximum temperature increment or cancellation
publishes no candidate and never changes the accepted gas or acoustic bank.
The native joule balance remains separately available.

Alternatively, the load's original `gas()` and `absorbed_energy_j()` can supply
the existing `FiniteGasThermalEnvironment` gas and prescribed gas heat, so its
native reciprocal body/gas exchange solves the combined balance. That path must
start from the original parcel, not a candidate already heated with the same
load. Harmonic heat is frozen and cycle averaged over the interval; accepting a
candidate requires a new acoustic solve before claiming temperature feedback.
This is uniformly mixed finite gas storage, without spatial gas conduction,
mean flow or transient acoustic onset.

Four G1/G3 integration tests in `tests/molecular_heating.rs` cover independent
positive-volume power, nonlinear Cv/u(T) energy and EOS closure, original source
retention, coherent cancellation/partial coherence/lossless limits, stale
source and volume refusal, caloric and relaxation endpoint domains, and
cancellation after the candidate solve. These synthetic source fixtures are
numerical conservation tests; native execution remains a separate gate.

## Shared affine dry-friction branch continuation

`stribeck_friction::StribeckLoadLine` also admits the original Coulomb and
clipped-linear `VelocityDependent` laws. Static capacity and the one-sided
kinetic intercept remain distinct. An occupied sliding branch continues inside
the larger static cone until its own zero-speed limit or finite-slip fold;
fresh contact and reversal retain the existing static-first initialization.
Exact affine roots consume no bisections and are checked against the unchanged
constitutive law. Strong weakening counts its unstable root without selecting
it. The zero-traction root at the clamp remains present at the fold itself.
A coincident flat sliding interval refuses as nonisolated, rather than claiming
a finite root count. Actual adhesive junction forces remain independent of
repulsive load, including zero repulsion. Arrest is checked on the actual
velocity line, so a rounded force quotient cannot invent a zero-speed slide.
Seven native G0/G1/G3 regression groups
cover these limits, both slope signs, reversal and adhesive offsets. This adds
no state-ageing law, thermal evolution or continuous-time passivity claim.

Single and simultaneous modal friction ports use this continuation for every
active dry law. Each accepted step publishes its own branch and root disclosure;
release retains the original zero-force path. Stateful shear layers and Mindlin
microslip remain with their original owners. Direct and interface-bound bow
declarations now admit Coulomb and clipped-linear laws through the same solver,
without translating them to Stribeck curves. Focused native integration fixtures
exercise kinetic continuation inside the static cone, opposite port directions,
actual mechanical work and the independent exact-oscillator load line.

The simultaneous dry-friction active-set rescue retains zero-load contacts as
zero-traction released rows, including on the condensed normal-obstacle line.
Their relative motion is unconstrained; loaded neighbours retain the original
cross-contact compliance. Enumeration counts each released row only once,
independent of its former sliding direction. Zero kinetic friction permits
zero-traction sliding, and no fixed small-speed cutoff erases a nonzero slip.
G1 regressions compare the coupled released row with analytical elimination and
the actual modal rescue with the converged joint sweep, including its work
balance. Adhesive and stateful shear laws keep their separate owners.

## Source-backed ionic phases during finite moisture exchange

`FiniteMoistureSolution::with_electrolyte_solubility` retains the original
`ElectrolyteSolubilitySolution` while the existing finite solution and native
slab or Cartesian specimens exchange solvent. Every candidate volume updates
the declared source-activity basis, solves pH and acid-family partitioning with
the selected finite salt, and obtains osmotic pressure/water activity from
the actual dissolved particles. Analytical cation/family totals, fixed
spectators, species identities and original source selection remain bound.
Concentration can create a solid; dilution can exhaust its finite stock.

The existing shared scalar water closure and native material/sorption
queries remain the sole moisture integrator. This path serves single and
multiple specimens, their adaptive continuations and existing one-way solid
mechanics consumers. Cancellation checkpoints also enter the nested chemistry
solve; refused candidates leave the original solution and specimens unchanged.
The accepted ionic owner and native fields publish together, with exact
analytical-stock/layout checks and the original free-water conservation check.

Only anhydrous selected salt reactions without a solvent term are admitted.
Hydrate precipitation or hydrolysis with chemical water transfer requires an
additional reaction-water account and refuses at construction. The original
fixed-solvent dilute acid/base approximation remains explicit, including the
small pure-water autoionization particle baseline retained by its actual H+/OH-
inventories; the solvent activity is not exactly pure-water-normalized. This is not a
concentrated-electrolyte activity model or an ion-permeating boundary. No
volatile solute transfer, nucleation kinetics or dissolution heat is inferred.
Focused native tests manufacture both precipitation on uptake and exhaustion
on desorption using the original source-qualified 1:2 salt, then check the
same endpoint pH, component totals, source RH and finite water balance.

## Native converter conduction and finite material temperatures

`circuit_enthalpy::converter::ConverterThermalNetwork` joins the native
three-state buck/boost stage to the existing finite-enthalpy/Newton owner.
Its five explicit receivers retain the original winding, input-high,
input-low, output-high and output-low channel order. Every channel has a
distinct original caloric/material body, including inactive channels; other
finite thermal cells and source-derived conduction links remain in the same
thermal inventory.

`ConverterResistance` separates caller-admitted fixed values, whole-body
uniform bulk resistivity/conductivity paths through `MaterialOhmicBody`, and
explicit static scalar resistance keys in ohms. The latter consume the
original material resolver, claim selection and complete source domain. They
do not reinterpret resistivity, efficiency or switching energy as resistance.
No semiconductor card or device law is fabricated. Unowned electrical-field,
current, gate, rate and frequency conditions refuse even when they appear
only in source applicability. Sources can depend on the actually owned
temperature and explicitly declared liquid mass fraction; no phase mixture
rule is introduced.

`bind` retains the original L/C and exact flux/bus-charge coordinates while
attaching those source laws to actual caloric state. Each `step_candidate`
trial refreshes the same original claims at its proposed thermal endpoint and
restarts the native electrical trajectory from the same accepted coordinates.
`ConverterElectricalDrive` selects held-duty averaging or a complete prescribed
switch sequence. Both feed the same thermal residual and native work ledger.
Each exact electrical interval counts against the shared circuit-trial budget;
a switch cycle must span the thermal interval within f64 duration and absolute
endpoint-clock roundoff. A period unable to advance the original clock refuses;
the authored switching durations are never rescaled to fit the thermal step.

Only the five native integrated conduction heats enter the cells, once, as
interval-average thermal rates. They are checked against native dissipation,
while both signed terminal works independently check native supply. Electrical
storage change uses the original small-transfer-aware native record, never
the difference of separately rounded large endpoint energies. The final shared
gate closes retained electrical storage, finite thermal storage, terminal work
and other original thermal exchange before either history publishes.

L/C and geometry remain fixed; no reversible parameter heat is inferred.
The averaged lane omits unresolved PWM RMS loss; the switched lane resolves
its actual current ripple and each channel's conducting intervals. Neither
adds dead time, recovery, diode clamping, gate drive or switching-transition
heat. A resistance source is a bounded on-state law, not a complete device
model or safe-operating-area claim. Thermal endpoint feedback is first-order
in general and needs refinement independently of the energy gates.

`tests/converter_enthalpy.rs` includes an independent scalar nonlinear energy
oracle for temperature-dependent freewheel decay, continued material/current
history, all-five-channel forward/reverse steady conversion, and actual
switched-ripple deposition with endpoint source feedback. It checks source
units, unambiguous receivers, event-period agreement, work budgets, source
domain and cancellation before joint publication. Its material data are
manufactured numerical fixtures, not measured-device validation.

## Source-bound moving-current feedback into native rigid momentum

`circuit_enthalpy::conduction::spatial::mechanical::PreparedSpatialConductionKick`
closes the original moving-conductor current and six rigid velocity components
at one implicit-midpoint force substep. The exact `TetrahedralMass` mesh and its
retained material-card/state identities must match the complete electrical mesh.
Its original `ResolvedBody::inertia()` enters the native articulated body without
an authored replacement mass, diagonal inertia or fitted motor constant.

Fixed conductivity, imposed B and native linear voltage/current traces make the
Lorentz wrench affine in rigid velocity. Preparation obtains the source wrench
and six homogeneous response columns through seven complete native field solves.
It factors `M - dt*A/2` without deleting couplings or modifying drag. A warm kick
uses one complete electrical solve at the coupled midpoint, the original P1
Lorentz-load/wrench map, and `fs_mbd::apply_free_floating_impulses` for actual
linear/angular momentum and polarized kinetic work. Source and response fields
are released promptly; the prepared cache retains only the small coupling and
the original immutable owners, plus held boundary/B snapshots.

The actual endpoint midpoint must agree with the field velocity. Independent
Lorentz/kinetic and `delta_K + Q_Joule - W_terminal - W_external` gates precede
publication. Only conductor-frame Joule heat belongs in a thermal receiver;
mechanical output changes kinetic storage, and a held external wrench has its
own signed work. Neither Maxwell/induction loads nor the native impulse may be
applied again. Refusal/cancellation changes no accepted state and permits retry.

This is explicitly a fixed-configuration momentum substep. It does not advance
pose, free rigid drift, gyroscopic transport, induction, thermal state or material
occupancy. Changing geometry, sources, B, boundary law, pose or duration requires
new preparation; a changed incoming pose refuses. The caller owns transport and
the physical validity of frozen geometry. Finite electrical contacts and joint
constraints require their own coupling laws and are outside this lane.

`tests/spatial_conduction_mechanical.rs` specifies independent uniform-bar
braking/generation/motoring and external-work oracles, source-density and mesh
variation, actual current-driven voltage feedback, a non-diagonal rotational
response with world-frame momentum transfer, open-circuit screening, second-order
fixed-operator time refinement, and source/budget/cancellation/retry refusals.
These manufactured-law fixtures do not qualify an empirical electric machine.

## Physical harmonic pressure reaches the existing WAV stream

pcm_wav::observation::harmonic::HarmonicPressureRenderer borrows an original
PreparedHarmonicField and its explicit source-clock epoch/output rate. It
implements the existing PressureRenderer contract and uses the sole PCM/WAV
writer, including authored pascal/full-scale conversion and clip reporting.
No gain normalization, propagation delay, attack envelope or source time
integration is introduced. Source and propagation validity remain the original
field caller's responsibility, including the 2D force model's scaling limit.

Construction checks every physical line against output Nyquist and reserves one
bounded callback buffer. Sampling uses the native aeroac/FFT DD phase owner and
absolute integer indices. Successful callbacks allocate nothing and are
bit-identical under block partition changes. A refused callback leaves both
the caller's output and completed sample counter unchanged. Existing streaming
cancellation and renderer/writer clock agreement remain authoritative. Native
regressions compare complete WAV bytes with the original encoder, verify
partition invariance, and exercise rate, size, alias and overflow refusals.



## Coherent winding heat on the native enthalpy timeline

The opt-in `harmonic-winding-heating` feature adds
`MagneticThermalNetwork::step_harmonic_winding_candidate`. The original
`HarmonicWindingExposure` integrates its resolved AC fields on exactly
[thermal.time_s(), target_s]; carrier/beat interference and same-frequency
cancellation survive into the receiving cells. Eddy, driven-solid and stranded
heat use original source/triangle identities and whole-volume routes. Physical
Joule heat is deposited once. Terminal work and incremental magnetic energy
change enter the same local/global native thermal audit.

The existing material network, finite caloric stores, source-based conduction,
boundary handling and backward-Euler/Newton solver remain authoritative. The AC
observation is computed once and borrowed across material/Jacobian trials; no
field re-solve or spectral reintegration occurs there. `circuit_trials` counts
material/thermal evaluations in this adapter. Consecutive intervals retain the
original physical phase clock.

This is a frozen-electromagnetic-law receiver: thermally mapped electrical
sources must be independent of evolving temperature, refreshed coefficients and
selected sources must agree, and fixed-magnetic-law admission remains mandatory.
Thermal transport may depend on actual trial temperatures. Whole solid volumes,
constant geometry/density, original caloric/card identities, source support,
events, work budgets and cancellation retain their refusals. Temperature-dependent
electrical/magnetic feedback, DC operating heat, switch-on history, hysteresis,
moving geometry and finite-amplitude saturation require other owners. Thermal
transport accuracy remains first order and requires time refinement.

Four G1/G3/G4 cases in `tests/harmonic_winding_heat.rs` cover independent two-body
temperatures and conduction, physical source work, consecutive interval
additivity, complete cancellation and transactional source/geometry/event/
budget/interruption refusal. Native workspace execution remains a separate gate.

## Molecular graph propagation reaches loaded exterior sound

`driving_point::network::radiation::RadiatingNetwork::bind_molecular` binds original
molecular reservoirs and optional joint thermal/shear diffusion to every edge of
an existing regional radiation owner. The graph and complete gas records come
from that same immutable owner; each relaxation state must match its original
edge gas. Explicit NetworkMedia records are required even for uniform gas.
Independent exterior gas, common static pressure, block limits and wall
resolution survive, and wall spectra are prepared once.

`MolecularRadiatingNetwork::prepare_sources` and `prepare_sources_with_radiation`
evaluate those laws at the physical angular frequency and reuse the fully loaded
source solve with the complete mutual BEM radiation matrix. The same opening
solutions and solved outward flows supply coherent pressure/velocity, source
statistics, directional microphone cross spectra and optional far-field
powers/moments. Source-bank coupling, signed work, retained-wave pole handling,
cancellation and allocation-free amplitude application keep their original
owners. Frequency sweeps reuse molecular spectra; amplitude/statistics sweeps
require no additional field solve. The result retains `molecular_network()` and
reports `loss_model() == None`; conventional results return `Some(original_loss)`.

The exterior remains homogeneous and nondispersive. This is frozen, stationary,
linear harmonic propagation, without mean flow, exterior molecular absorption,
transient onset or causal time realization. Original plane-wave, mesh/port-basis,
diffusion-band and material limits apply; embedded compact tone holes refuse.
G1/G3/G4 cases in `radiation/tests/prepared/molecular.rs` cover empty-reservoir
parity, analytic dark-mode absorption/work, collision-rate intervention,
independently assembled receiver fields, coherence through directional/angular
observers and original-record/budget/band/cancellation refusals. These numerical
fixtures do not qualify an empirical material or instrument.

## Declared isotropic sources in bonded string input

The layered `material_strings` frontend accepts an optional `elastic-sources`
block with one `elastic-source DENSITY PAIR FIRST SECOND` row per bonded
core/sleeve, in spatial order. The explicit pair is `young-poisson`,
`young-shear` or `bulk-shear`; original property keys, admissibility, material
queries and claim receipts feed `LayerIsotropicParameters` and the shared
concentric-section projection. The resulting mass, EA, EI, rotary inertia and
GJ enter the original coupled performer and pressure output.
`MaterialStringSource::with_isotropic_torsion` uses that same section's derived
or directly selected G without a fabricated shear claim; viscosity remains
independently sourced. Missing, incompatible, repeated-role or incomplete
declarations refuse. Omission retains the original canonical source path.
Separate contact, caloric and memory bindings retain their own source requirements
and validity limits; this input block adds no material law or measured-property
claim.

## Exact-ramp work and heat in deposited layers

`strip_enthalpy::step_candidate` uses the solid owner's native work-conjugate
resultants at prescribed strain/curvature constraints and the caller's explicit
mean force/moment for mean-load controls. Endpoint reactions stay available but
are not substituted for interval-mean work. The existing material histories,
source identities, reversible heat and physical dashpot loss feed the same
thermal cells and nonlinear transaction. Backward-Euler loads retain their
original work and numerical-loss convention.

G1/G3 coupled tests compare exact Maxwell ramp heat against an independent
two-cell heat-network solve and check mean-load work and the combined first law.
Exact-ramp material coefficients are still frozen at the proposed temperature;
spatial heat transfer retains its original backward-Euler discretization.

## Original elastic sources across live material coupling

Modal wet-frame temperature refresh and native rod enthalpy/viscosity preparation
retain each layer's explicit density and elastic source roles (`E/nu`, `E/G` or
`K/G`), original claim selections, geometry and sourced-or-derived shear binding.
Section enrichment and refresh use the original section constructors' retained
roles; derived Young/shear moduli do not become measured claims or canonical-name
substitutes. Rate-only Maxwell updates require constant density and every actual
elastic storage parent, including both K and G when they determine E, and preserve
native history, storage, fluid state and clock at the accepted boundary. Rod
viscosity preparation retains these same mechanical sources while adding only the
selected loss property and depositing its independently accounted heat. Concentric
string caloric admission joins the section's selected density receipt to the
original caloric source, mass, volume and fixed query conditions, including when
the density property has an explicitly selected noncanonical name. Existing
source-support, affine thermoelastic, caloric and energy limits remain binding;
refused source changes publish no coupled state.

## Loaded quadratic Maxwell scenes through native force, contact and drift

`trajectory::quadratic::sustained::mixed_visco_sustained_sphere_trajectory_candidate`
advances original curved P2/P1 Maxwell motion with a complete prepared sphere
scene through native force kick, simultaneous contact and material drift. Every
active solid and declared rigid-pair row uses its original law, witness, solid
trace and full rigid inertia. Restitution targets are computed from accepted
pre-force velocities and fixed across retries, so a resting loaded contact does
not turn each gravity kick into a new restitution bounce.

The original finite-compliance mixed owner determines pressure at the actual
post-contact drift coordinates and advances each retained quadrature history.
Physical storage includes `U(J(p))`; reactions from that pressure solve's fixed
coordinates never become physical support work. Authored stationary supports,
signed nodal/reference-volume loads, wrenches, gravity and actual support torque
retain separate momentum accounts. The original Tet4 and new TET10 sustained
drivers share one accounting implementation through generic path receipts;
existing Tet4 public receipt names remain concrete aliases.

Positive reference and current whole-material maps are required before the first
force phase. Every accepted interval then admits the complete TET10 volume path,
curved boundary paths and declared sphere-pair paths. Each rejected interval
recomputes the entire force/contact/drift sequence from the same accepted state.
No rejected pressure solve, contact trial or material endpoint publishes time,
stress or history. Accepted coordinates and geometry are reused. The absolute
splitting allowance is distributed in proportion to trial duration; signed and
noncancelling momentum/work checks also apply to the complete requested call.

Relaxation, contact-map loss, source/force splitting and rigid drift defects
remain separate. Loaded support contact can contain refinable numerical impact
loss; it is not declared friction heat. Scope is first-order finite-compliance
motion, stationary held velocities and one COM-centred sphere per rigid body.
Exact incompressibility explicitly refuses the unconstrained kinetic subflow.
Geometric admission covers the declared affine numerical path and local material
orientation, not continuum-flow CCD or global solid self-intersection.

`trajectory::quadratic::mixed_visco_sphere_scene_trajectory_candidate` also adds
the same complete scene to the existing native transient-impact driver. Declared
pair encounters localize before entering the simultaneous solve and retain their
own restitution/friction law. Accepted interval keys align with the original
contact impulses. The surface-only entry retains its undeclared-pair refusal;
both entries share the same native motion and physical accounting worker.

Five focused native regressions specify loaded curved two-sphere dwell with
pre-force restitution, complete history replay and timestep-refined support loss;
adaptive/late/cancelled rollback; sliding with original off-diagonal inertia and
actual activation-gap angular residual; unequal-mass pair impact with independent
impulse/loss/trajectory values; and single-sphere parity between transient entries.
Pinned parser/format checks, bounded owner-integration review and independent
mechanics arithmetic passed. Native Cargo compilation and Rust test execution
remain unverified under the workspace resource blocker.

## Complete curved solid and rigid-pair scenes

`inelastic_impact::surface::quadratic::scene` joins original curved TET10
material contacts and explicitly declared rigid sphere pairs through the same
prepared scene owner. It retains every signed native solid trace, actual paired
witness, full rigid inertia and original interface law. Canonical solid IDs
remain even and pair IDs odd. The shared scene query now carries its original
surface receipt type, so quadratic geometry never becomes a fabricated linear
surface. Maxwell and conservative native geometry both use this adapter.

Complete affine scene paths combine the existing curved-boundary query with
swept-box pair pruning and the original directed sphere-path owner. Visits,
candidate storage and time windows share the enclosing scene allowances.
Clearance bounds include culled pairs; earliest unresolved windows retain their
original face/pair key with deterministic ties. Initial actual pair penetration,
missing pair laws and exhausted work refuse. Local volume regularity remains a
separate admission, and the affine path is not a continuum-flow certificate.

Three focused native regressions cover a simultaneous curved-support/pair impact
with independently derived impulses and support torque, native geometry and
history identity, canonical permutation, hidden pair crossings, shared work,
exact touching recession and one-ULP penetration. Pinned parser/format checks,
source review and independent impulse/path arithmetic passed. Native compilation
and execution remain unverified under the current disk/dependency-memory block.

## Native quadratic Maxwell trajectories through curved impacts

`inelastic_impact::trajectory::quadratic::mixed_visco_sphere_trajectory_candidate`
advances the original finite-compliance P2/P1 Maxwell motion and COM-centred
rigid spheres through simultaneous curved-surface impulses, native material
steps and native rigid free flight. The original consistent mass, complete
inertia tensor, pressure gauge, material laws and quadrature histories remain
with their existing owners. Stationary held velocity DOFs retain their initial
prescribed displacement; external nodal/body loads, wrenches and gravity are
explicit and fixed for the call.

The complete reference and initial current TET10 maps must be admitted before
the first impulse. Each trial then requires whole-material volume-path and
curved boundary-path admission, plus independent sphere-pair sweeps sharing
the boundary work allowances. Unresolved volume/boundary paths or ambiguous
endpoint activation shorten the interval and recompute the native solve from
the same accepted post-impact state. Rejected trials contribute no history,
loss or accepted time. Refusal retains the unresolved evidence, and caller
owners remain unchanged unless the complete returned candidate is adopted.

Both Tet4 and TET10 trajectories share one physical accounting implementation.
The existing Tet4 receipt names remain concrete aliases. Actual linear/angular
momentum changes, applied/support impulses, total material storage including
bulk, physical relaxation and contact loss stay separate from backward-Euler
kinetic loss, angular drift and endpoint-work defects. Signed final balances
and accumulated residual magnitudes both obey complete-call limits. Accepted
coordinates and endpoint geometry are reused between intervals.

This driver supports finite compliance, stationary supports and one
COM-centred sphere per rigid body. It explicitly refuses exact-incompressible
backward-Euler flight because the endpoint displacement chord does not prove
the instantaneous volume-tangent velocity needed by constrained contact.
Undeclared sphere-pair encounters, rotating offset centres, unresolved
persistent/grazing contact and exhausted budgets refuse. Geometric admission
covers the declared affine coordinate path and local material regularity,
not the continuum trajectory or global solid self-intersection.

Native regression source covers a curved supported bounce with history replay,
free recoil and continued pressure/history evolution, original signed TET10
body loads with full rigid wrench/gravity, independent numerical-defect gates,
and transactional late/pair/cancellation refusals. Pinned parser/format checks,
native-owner source review and independent exact-rational impulse/loss
calculations passed. Native compilation and test execution remain unverified
under the workspace disk/dependency-memory blocker.

## Whole-material quadratic volume admission

`inelastic_impact::surface::quadratic::volume` binds a prepared volume check to
the original immutable TET10 surface owner and its complete element incidence.
Preparation admits every complete reference map once. Subsequent affine path
queries retain every original node and check every element through the original
`fs-query` quadratic-volume owner, with one cell-work allowance shared across
the mesh and bounded per-element frontier storage. Preparation clones neither
coordinates nor material histories.

A positive result covers closed material domains, including edges and vertices;
positive corner geometry, face samples or quadrature samples are insufficient.
An unresolved result retains the earliest uncertified interval, canonical
element ID, material-coordinate bounds and scaled determinant bounds. Those
determinants use each element's own normalization and are not physical cubic
metres. Cancellation, invalid complete coordinate images, exhausted work and
unresolved reference geometry cannot return admission. Local orientation alone
does not establish global absence of material self-intersection.

Focused native regression source covers a hidden internal fold with unchanged
corners, a reference fold missed by the existing sampled checks, and a
two-element isochoric curved path with shared work limits. Source/format and
independent analytic checks passed. Native compilation and execution remain
unverified under the current workspace resource blocker.

## Material carriers shared by moving film patches

`MaterialFilmCarrierSi` binds each physical registry body to its original
`ResolvedBody` mass/full inertia and accepted SI state. `MaterialFilmPatchSi`
separately describes prescribed reservoirs or body-coordinate patches attached
by explicit carrier index. `ResolvedImmersedBodies::with_shared_surfaces`
retains the existing dry-sphere/fluid history and creates
`ResolvedImmersedSharedSurfaces`. Several patches may reference one carrier;
equal material references in different registry entries remain independent.

One shared native dry-body correction and resistance solve advances spheres,
fluid, registry bodies and prescribed patch positions on the same accepted
clock. Each carrier's complete external force/couple is converted once with
the original shared SI/lattice helper. Its mass, full inertia, free half-drifts,
film impulse and kinetic-storage account are never repeated per patch. All
attached faces use the common accepted pose and solved midpoint rigid velocity.
Source definitions remain immutable and observers read accepted state.

Reports separate per-patch reaction/work from per-carrier aggregate momentum,
kinetic storage, authored work and free-rotation defect. Existing assembly
converters interpret the unchanged native lattice units. The wrapper forwards
the original initial-flow, held-fluid-force, live open-boundary target and
density-preserving viscosity updates; every normal/tangential film reads the same fluid snapshot. Invalid
indices, unused registry entries, bad load counts, cancellation or rejected
native probes preserve the accepted states and clock. There are at most 64
registry carriers and 64 patches, with original native dry/body/film budgets.

This remains a local planar missing-film model. Carrier material supplies
inertia, not certification of patch extent or a resolved no-slip wall. There
is no inferred gravity/buoyancy, carrier sampled fluid interior, contact, GCL
or finite heat destination. Complete carrier external loads are explicit.
Focused material regressions use the original source resolver and native
solver to check multi-face response at the shared midpoint velocity, unique
storage and actual momentum balance, legacy equivalence, coefficient refresh,
source identity and transactional retry.


## Sealed finite gas and full nonlinear membrane equilibrium

`cavity_phs::membrane::GasMembrane` joins the original finite-deformation
`MembranePatch` to an original `ClosedGasParcel` or `ClosedRealGasParcel`.
The caller anchors the parcel's actual volume to a reference shape, supplies
absolute exterior pressure and a positive volume interval, and fixes the
complete membrane perimeter through the original incidence admission.
Orientation points out of the gas. Every trial retains those prescribed
coordinates and obtains its swept volume from the membrane's exact
unit-pressure work.

The original material isentrope supplies pressure, temperature, gas energy
and entropy at conserved mass and frozen composition. Ideal NASA calorics
and homogeneous Peng–Robinson departures retain their original source
domains. Native `GasMechanicalLoad` pulls this pressure through the current
volume gradient and Hessian. The force tangent contains both the geometric
pressure term and `-(rho*c^2/V) grad(V) grad(V)^T`, including nonzero gas
stiffness at zero pressure jump. The original gas owner limits the chart to
1..=256 free Cartesian coordinates; membrane mesh/material admission remains
separate.

The conservative potential is the membrane storage change plus native gas
Delta U, exterior-pressure work and negative dead-load work. The held-pressure
`p_final*Delta V` expression is not added again. The original Newton–Krylov
owner closes the force residual with that same tangent. Repeated Krylov
actions reuse an unchanged physical trial; rejected source/geometry proposals
can backtrack. Cancellation, exhausted budgets and unresolved final residuals
publish no coupled state. Full pressure loads remain on support-reaction rows.

Focused synthetic tests cover independent pyramid equilibrium/work,
variable-Cp entropy, native PR pressure/bulk/work, a six-coordinate geometric
tangent/potential, gas stiffness at zero pressure jump, backtracking and
unchanged reference inventory on refusal. Their presence does not assert
native execution. This is reversible static equilibrium, not a heat/mass-transfer
or finite-time integrator, gas-wave solver, vibration/stability certificate,
or proof against enclosure contact and global surface intersection.

### Small vibrations of the same sealed gas and membrane

`GasMembrane::solve_vibration` first solves the original coupled equilibrium,
then reuses that one accepted source-bearing gas state for every tangent column.
The complete pencil is
`K_membrane - (p-pe) Hess(V) - (dp/dV) grad(V) grad(V)^T`; the original ideal/PR
acoustic bulk modulus supplies `dp/dV`, including nonzero reservoir stiffness
at zero pressure jump. Original membrane reference inertia, material/geometry
admission, full-spectrum positivity, residual gates and nodal projection/
reconstruction remain with `MembraneModes` and `fs-modal`. The dense leaf
retains its 128-free-coordinate cap.

The result exposes the same immutable equilibrium and native modal basis.
Volume perturbation is `B.q`; gas-pressure perturbation is `(dp/dV) B.q`.
This observes restoring pressure already included in the pencil, rather than
adding another gas spring. Independent signed pressure-jump loads use the same
work-conjugate port.

Accepted volume and temperature must be strictly inside the declared geometric
interval and entropy-solve bracket. A continuous homogeneous source branch
around equilibrium remains caller authority; local checks do not certify
unqueried source support or global phase stability. There is no added gas
inertia, wave propagation, damping, thermal exchange or finite-amplitude time
integration. Refusals publish no modal result or source mutation. Focused
synthetic regressions cover analytic zero-jump stiffening with variable Cp,
all six modes of a curved pressure-loaded shape, native PR departures,
reciprocal pressure observations, endpoint/budget/cancellation refusal, and
unstable/nonconservative pencils. Test presence does not claim execution.


### Sealed-gas membrane motion

`cavity_phs::membrane::GasMembraneMotion` advances the full admitted nonlinear
membrane through the original `fs_time::galpha::OperatorGeneralizedAlpha` owner
at `rho_inf = 1`. It retains the original fixed perimeter, reference nodal
masses, membrane force/tangent and source-owned isentropic gas response,
including variable caloric data and the admitted Peng–Robinson branch.
Every actual shape determines its enclosed volume and absolute gas pressure.
Original gas compressibility and follower-pressure geometry enter the same
native tangent; pressure is not added again as external forcing.

Initial, generalized-alpha stage and accepted endpoint shapes must pass
original membrane, volume, temperature and material/source admission. One
original gas trial is cached per coordinate point for force/tangent reuse.
Physical line-search refusals remain recoverable; cancellation and caller
evaluation limits stop the proposal. Mechanics, gas, clock and native history
publish together only after all endpoint and energy checks. A refused proposal
preserves the complete accepted state. Caller controls bound physical
evaluations, retained steps, Krylov restart and the permitted energy defect.

Each accepted interval reports kinetic, membrane and native gas storage
changes, gas work, exact geometric swept volume, exterior-pressure work and
held dead-load work. The measured defect is
`Delta K + Delta E_membrane + Delta U_gas + p_exterior Delta V - W_dead`;
gas work is `-Delta U_gas` and is not debited twice. The caller must select
an acceptable defect and timestep: nonlinear midpoint integration is not
exactly energy conserving, and its defect is not physical heat or damping.
Stored endpoint acceleration is the native algorithmic quantity; endpoint
force residual does not certify static or instantaneous acceleration equilibrium.

This bounded path retains homogeneous sealed gas, stationary supports and the
existing free-coordinate cap. It does not introduce gas inertia/waves, heat or
mass exchange, contact, material damping, prescribed support motion or topology
changes. Five focused native regressions cover small motion with a changed held
load, finite-amplitude timestep/energy refinement, curved full-coordinate
midpoint pressure balance, original variable-caloric/PR endpoint states, and
complete rollback after energy, source, cancellation or work-budget refusal.
Test presence does not claim native execution.


### Heat input during continuing membrane motion

`GasMembraneMotion::into_thermal` moves the complete accepted native mechanical
state, clock, telemetry history, reference masses and current parcels into
`GasMembraneThermalMotion`. It does not restart motion or recompute algorithmic
acceleration. Each thermal interval supplies signed heat to the first parcel
and, for a sealed opposing chamber, separately to that parcel with independent
`GasEnergySolve` controls. Fixed exterior pressure retains the existing
reservoir-work account.

The original ideal/PR `midpoint_heat_volume` owners solve heat and volume work
together. Their midpoint pressures drive the original nonlinear membrane at
the native generalized-alpha midpoint, and their complete endpoints become the
next accepted gas states. The tangent differentiates this discrete law at fixed
interval heat: endpoint coordinates are `q1=2*q_stage-q_old`, so the native
`dp_mid/dV1` contributes a factor-two cross-gradient action using endpoint and
stage volume gradients. Opposing-gas stiffness and follower geometry occur
once. This action is generally nonsymmetric; no isentropic tangent is used
for heated trials.

For curved membrane geometry, midpoint traction work need not equal pressure
times the exact swept volume. The runtime reports and separately gates this
pressure-work mismatch and the enclosing
`Delta K + Delta E_membrane + Delta U1 + Delta U2 + W_reservoir - W_dead - Q1 - Q2`
defect. Both allowances must be finite and nonnegative. Original gas receipts
retain actual heat, work, internal-energy change, entropy change and discrete
pressure derivatives. Neither numerical defect is a heat source or damping
model, and finite midpoint steps do not certify isentropy or nonnegative
entropy production.

Geometry and source admission, bounded native Newton/Krylov work, both energy
gates and final cancellation precede atomic publication of all participants.
Repeated stage/tangent calls reuse their native gas proposal. Heat is
prescribed interval energy; this API does not infer heat-transfer coefficients,
finite heat-source dynamics, gas waves, composition change or phase equilibrium.
Time refinement remains necessary. Focused heated-motion regressions are
authored; this checkpoint does not claim native execution.

### Two finite gas chambers sharing a membrane

`GasMembrane::with_sealed_exterior` replaces the held exterior-pressure reservoir
with a second original ideal or Peng–Robinson parcel. `SealedMembraneExterior`
supplies that parcel's independent positive volume interval and native
entropy-solve controls; existing per-call controls continue to govern the first
parcel. The two source-temperature domains may be disjoint. Both chambers share
the same fixed-perimeter membrane geometry: its native swept volume `s` gives
`V1 = V1_reference + s` and `V2 = V2_reference - s`. Neither parcel exchanges
heat or mass, and both retain their source, composition, mass and reference
entropy.

The shared trial supplies actual `p1 - p2` traction and the tangent
`K_membrane - (p1-p2) Hess(V1) - (dp1/dV1 + dp2/dV2) grad(V1) grad(V1)^T`.
Native gas-load actions provide both compressibilities. The geometric pressure
term is counted once using the pressure difference, including near pressure
balance. Both source states are cached at the same shape and reused by the
existing equilibrium, modal and generalized-alpha owners. Reported gas-trial
counts denote coupled geometry proposals, each reconstructing one or two
parcels with their independently bounded original entropy inversions.

Evaluations expose both native gas states and volume-work reports. Potentials
and motion defects include both internal-energy changes once; separate
opposing-gas energy/work fields distinguish them from fixed-reservoir pressure
work, which is zero for the sealed exterior. Both accepted volumes and their
independent temperature brackets must admit local modal perturbations.
`pressure_perturbation_pa` retains its first-parcel meaning; explicit
exterior-pressure and pressure-jump observations include the equal-and-opposite
volume changes. Their restoring forces already belong to the modal pencil.

All existing geometry/source, time-step, energy, budget and cancellation gates
precede atomic publication of mechanics and both gases. A refusal from either
side preserves the accepted state/history. Existing homogeneous-gas,
stationary-support, coordinate-cap and numerical-energy-defect boundaries
remain unchanged. Four focused native regression cases cover the new physical
behavior; their presence does not assert execution.

## Original spatial gas sources in three-dimensional pressure acoustics

`cavity_modal::interior::heterogeneous::AxialGasInterior` binds the original
`AxialGasField` to the original tetrahedral pressure space through a physical
axis in metres. At each volume quadrature point the material owner reevaluates
temperature and mole fractions together, preserving its pressure, thermochemical
cards and phase policy. The FE owner assembles `1/rho` stiffness and
`1/(rho*c*c)` mass from that one coherent state. It does not reuse the geometry
owner's uniform reference medium or interpolate derived density and sound speed.
Prepared sparse operators serve full harmonic pressure solves and rigid-wall
modes without resampling material data during frequency/drive solves.

The mesh must lie inside the supplied source interval. True cellwise jumps must
conform to tetrahedron faces; unresolved jumps refuse instead of depending on
which side a quadrature point happens to sample. Continuous primitive gradients
retain explicit quadrature and mesh refinement requirements. This binding is a
frozen stationary, isotropic, lossless bulk fluid without mean flow. Local passive
boundary admittance and imposed inward velocity use the FE owner's original
pressure traces and work balance. It adds no heterogeneous thermoviscous closure
or thermal evolution. The `heterogeneous_acoustics` example accepts original
gas source packs, drives a spatial pressure field, and exports pressure and
boundary work; its local outlet impedance is not exact multimode radiation.

`heterogeneous::regions::FluidRegionsInterior` binds original `FluidAcoustics`
providers to explicit conforming tetrahedral regions. Gas, liquid, dense-fluid
and material-card sources retain their own phase/relaxation policies and original
receipts. Each used region is resolved once into a coherent acoustic snapshot;
unused sources are not queried. The original heterogeneous FE owner assembles
both physical matrices, including inverse-density interface coupling, with no
homogeneous reference-fluid substitution. Region maps and source failures remain
explicit. Changing an original source requires preparing a new frozen binding.

Both spatial source bindings expose the original operator's sparse transient
pressure-potential evolution as well as its harmonic and modal consumers.
`fluid_interface_pulse` launches a compact oblique packet across a conforming
original cold/hot gas or IF97-water interface, exports actual nodal pressure and
potential, and reports physical acoustic energy and boundary work at every step.
Its exterior local impedances are exact only for normal plane waves. The example
requires independent spatial and time-step refinement; source thermodynamics,
interfaces and bulk coefficients remain frozen throughout the pulse.

`heterogeneous::mechanical::MaterialMechanicalTime` binds original
`LocalMaterialWall` sources to the same full spatial time operator. Each patch
explicitly supplies its cylindrical axis, radius and relative geometry tolerance;
the actual boundary vertices and outward normals must approximate that
cylinder's lateral surface. Caps, incompatible radii/normals and non-thin source
walls refuse. The source's axial restraint and hoop/axial material directions
remain binding. Facet and time refinement remain caller responsibilities; this
is locally reacting cylindrical membrane motion, not arbitrary plate bending,
lateral shell waves, pressure prestress or evolving geometry/materials.

The existing source reduction supplies surface mass, equilibrium stiffness,
viscosity and every positive original Maxwell arm without fitting. Actual FE
face area enters once through the consistent trace mass. Original source
identities, exact-zero source records and the intersection of applicability
bands remain attached; a declared band is not an input filter. The original
Maxwell relaxation-time and instantaneous-frequency clock limits apply before
preparation. All states, memory storage, dissipation and pressure feedback use
the native implicit owner, including its resource, cancellation and source-state
binding refusals. `material_mechanical_acoustics` exercises source impedance and
reflection, native transient parity, physical memory work and source/geometry
admission with synthetic constitutive records.

## Actual curved TET10 contact through native momentum owners

`inelastic_impact::surface::quadratic::QuadraticSolidContactSurface` binds the
complete original reference coordinates and TET10 incidence. The existing solid
topology owner admits shared midsides; the original closed-tetrahedron owner
supplies exterior incidence. Actual six-node quadratic patches, rather than
corner planes, supply current distance, outward-side bounds and material points.
The prepared query BVH is refitted once per coordinate image and reused across
spheres. Sphere IDs, shared vertex/edge features and contact ordering are canonical.

Distance and witness-localization enclosures must resolve nonpenetration and the
authored activation band before mechanical admission. A potentially active patch
is never silently dropped because its enclosure straddles the threshold. Distinct
unresolved minimizers, backside ambiguity, native point Jacobian failure and
aggregate geometry-work exhaustion refuse. Raw negative gaps remain negative for
an external event locator. Signed native TET10 shapes reconstruct the real solid
point, and the radial rigid material arm retains its actual torque and rounding
error; neither witness is moved to manufacture coincident contact.

The native Maxwell and conservative-solid helpers borrow reference geometry and
accepted displacement without copying constitutive history. Detected material
points feed the existing simultaneous circular-friction contact owner and its
original mass response. Exact-incompressible Maxwell and conservative formulations
reuse the existing pressure-constrained response and complete physical impulse,
angular-momentum, support-work and energy accounts. Successful impulses retain
physical pressure, Maxwell memory, displacement and time through the original
complete-state candidate. No replacement mass, constitutive or contact law is
introduced.

`QuadraticSolidContactSurface::affine_sphere_paths_at_positions` reuses that
prepared curved topology for a complete declared path of every original node
and sphere center. The initial curved query admits exterior containment and
preserves negative-gap refusal. An uncertain nonnegative initial gap can be
resolved only by the subsequent whole-path proof, including exact touching and
receding support cases. Endpoint poses locate sphere centers; interpolating
those centers affinely does not assert an underlying rigid rotational path.
Endpoint Bernstein hull unions supply a conservative swept broad phase. Every
remaining face uses the original directed quadratic sweep, with aggregate
traversal, candidate, time-window and distance-cell limits across the scene.
The earliest unresolved window retains its sphere and face identity; a clear
receipt carries a finite lower clearance bound, including culled pairs.

The input solid must be a consistently oriented, conforming nonoverlapping mesh.
Positive native Jacobians at selected points are local admission checks; they do
not prove global injectivity or a valid material trajectory between states.
Focused native regressions specify a quadratic isochoric shear with unchanged
corners but a displaced contact patch, actual off-center rigid torque, conservative
pressure closure, retained relaxed Maxwell memory and atomic overlap refusal.

## Finite gas species and momentum in an explicit energy frame

`MovingGasParcel::exchange_species_momentum_relative` retains the original
finite gas inventory and momentum owner while applying integrated species
streams and their complete impulse. Energies are specified relative to original
partial enthalpies in one declared inertial frame; impulse uses the retained
momentum frame. A translating energy frame can preserve small thermal changes
beside large common drift when callers construct streams directly in that frame.

The original gas owner supplies its exact admitted final mass to the additional
energy callback. Actual endpoint mass and momentum determine kinetic change via
the changing-mass midpoint identity with compensated arithmetic. That kinetic
change is removed before ONE original caloric inversion. No absolute kinetic
totals are subtracted, and no extra receiver-mixing heat may be added afterward.
Per-component momentum budgets use exchanged impulse, not common drift or an
unrelated large pressure component. Vanished impulse, invalid mass, caloric
refusal and cancellation expose no partial endpoint. The enclosing phase/body
owner must admit gross withdrawals and close reciprocal energy and momentum.
Fixed gas volume, original thermochemistry and original auxiliary axes retain
their existing meaning. Focused tests specify changing mass, common inertial
boosts, native source reconstruction, exact idle state and atomic refusal.

## Stationary native gravity in evolving finite solids

`FiniteSolidBody::new_with_gravity` and `new_with_gravity_and_sources`
bind an explicit stationary native `PointGravityField` and inertial frame to
the original cellwise or P1 reference mass, optionally alongside current-area
pressure. Initialization and ordinary/energy-momentum motion reuse the original
generalized-alpha/Newton operator. Native gravity force and its analytic spatial
derivative are evaluated on the evolving solver stage; the material endpoint
chain factor is not applied a second time to the gravity derivative.

Gravity potential decrease enters source work exactly once. The full nodal
force retained from the accepted native stage feeds stationary support reactions
and mixed solid/rigid free-flight impulse and moment accounts. The binding
survives cloned and instantaneous-impulse candidates, retaining the field and
acceleration history through continued motion. Affine gravity has matching
stage-force and potential work to quadrature roundoff. Nonlinear fields retain
a separate signed quadrature defect and the existing mechanical energy gate;
that numerical defect is never physical heat.

The existing first-order force/contact/drift composition retains the same bound
field. Native gravity at the initial configuration supplies its consistent-mass
kick, external impulse/moment and mean-velocity work. Contact impulses may then
change velocity without changing the field or reference mass. The original
kinetic subflow advances configuration and time, and the native force/mass owner
restores endpoint acceleration with gravity at the new configuration. Original
potential decrease supplies gravity displacement work exactly once. Its mismatch
with kick work has a separate signed account and contributes its absolute value
to the existing split gate; opposing dead-load, pressure and elastic errors cannot
cancel it. Private kick bindings retain this accounting even if public diagnostic
fields are edited. These are first-order subflows requiring timestep refinement,
not an exact forced trajectory or a new contact law.

Field-domain, frame, cancellation and accounting failures publish no candidate.
Callers must not duplicate the bound field in held nodal or body loads.
Moving-source power, field rebinding and self-gravity remain outside this binding.

Focused native regressions cover free fall and global impulse/moment balance,
stiff tidal motion with original graded mass and pressure, impulse and split
continuation, independently computed numerical work defects, noncancelling split
gates, source-domain exits, edited diagnostic receipts and cancellation. Ordinary
midpoint material work defects are checked against independent scalar energy
accounting; the energy/momentum path retains its tighter energy gate.


## Finite gas and rigid wall under molecular shear

`molecular_shear::step_candidate` exchanges tangential momentum between the
existing `MovingGasParcel` and an original finite-mass rigid wall bound to its
enthalpy/phase chart. An explicit Maxwell reflection model supplies the same
species encounter masses and caloric reflected heat through the original
`ThermalAccommodation::evaluate_sealed` owner. There is no volatile placeholder,
mass transfer, replacement gas inventory or independent heat-capacity model.

For total accommodated encounter mass C and reduced mass mu, the original
relative velocity relaxes by exp(-C/mu). The reciprocal impulse and integrated
wall displacement follow that same exact frozen-coefficient solution. Relative
kinetic loss D enters each original caloric store as D/2, in addition to the
reciprocal original reflected heat. This split follows the Maxwell crossing
energy moments; neither kinetic loss nor source heat is counted twice. Actual
native wall impulse, represented gas momentum and each original thermal store
are independently checked using exchange-scale budgets. Large common velocity,
energy datum or normal pressure cannot enlarge the conservation allowance.

The plane is fixed normally and the wall cannot spin. Complete molecular normal
pressure, separate gas/wall confinement impulses and the patch's rotational
guide reaction are explicit and do zero work. The material patch arm is body
fixed and translates with the wall. All crossed caloric segments must retain
constant density, so latent plateaus are supported when they require no missing
expansion/normal-boundary work. Duration, endpoint temperature changes, species
and original chart-knot work are bounded. One original gas heat reconstruction
produces a joint gas/wall candidate; any late refusal or cancellation publishes
neither endpoint.

Temperature-dependent encounters and reflected heat are frozen over each
interval, so trajectory accuracy requires time refinement. The homogeneous
gas owner supplies no spatial transport, gas spin or unrestricted angular
closure. Focused regressions specify unequal-mass relaxation, variable calorics
and latent storage, full support/guide loads, sealed-owner parity, specular
limits, stiff coupling, inertial covariance and atomic source/budget refusal.

## Forced superfluid momentum and original caloric heating

`superfluid::advance_forced_closed_friction_heat` composes the authored component
body forces and full closed-channel momentum account with the same original
thermal timeline used by `advance_closed_friction_heat`. Shared private routing
and thermal-preparation helpers preserve the original caloric laws, sources,
contacts, events, clocks and per-cell heat admission.

External component-force work appears once in the original channel account.
Only the interval's actual mutual-friction loss is deposited in the original
thermal cells; it is an internal transfer in the combined account. A source
force is not replaced by fake heat power, and cumulative heat is not replayed.
Thermal refusal or final cancellation publishes neither source-driven momentum
nor thermal storage. The returned reports retain physical body and wall
impulses, signed source work, friction heat and the unchanged thermal audits.

The fixed-reference wave model and first-order endpoint heat split remain
explicit. Temperature evolution does not silently invent a new helium
equation of state or relinearize the wave coefficients. The example
`cargo run -p fs-couple --features superfluid --example superfluid_piston`
separately demonstrates a sourced 1.8 K helium-II channel driving a native
2 g piston, with independently recomputed combined momentum and energy.
Its geometry, excitation and friction coefficient are authored inputs;
exported heat is accounted for without claiming calibrated apparatus data.

## Reciprocal native superfluid piston

The opt-in `superfluid::piston::advance_balanced_piston` connects the original
two-fluid channel to the original `fs_mbd::RigidBodyIntegrator`. A centered,
axially translating native body drives an impermeable end; the other end is
already stationary. Normal and superfluid endpoint velocities share the
mechanical midpoint velocity. Reference pressure is explicitly balanced by
exterior pressure or preload; the computed loads are pressure perturbations.

The channel's exact two-mode endpoint Schur response solves pressure and piston
velocity reciprocally. Existing tridiagonal factors and cached endpoint inverse
diagonals are reused. The fluid half-cell mass rho*A*dx/2 is tracked separately
from the unchanged native body mass: the actual native body and that half cell
share the adjacent-cell pressure impulse and its work. Native position advances
through the same rigid-body integrator as native momentum. The opposite fixed
support's impulse remains explicit although its work is zero.

The combined wave, native body, half-cell and exported friction-heat balance
must pass authored tolerances. The adapter preserves both owners on native,
physical-domain, transfer or cancellation refusal, including the final
checkpoint. A previously moving opposite end or counterflow at the piston
refuses instead of silently deleting endpoint momentum. Axial motion, no spin,
authored equilibrium position and small displacement are explicit admissions.

This is a fixed-reference linear acoustic/mechanical composition. It does not
supply moving geometry, arbitrary rotations, an unbalanced absolute-pressure
load, or nonlinear thermal feedback. Friction heat remains available on the
original channel and is included once in the combined account; this adapter
does not silently deposit it in an unspecified thermal owner.

## Spatial plate and finite rigid-body normal impacts

`plate_impact::PlateRigidImpact` localizes one frictionless normal event between
an original `fs-plate-transient::PlateDynamics` point and an original
`fs-mbd::RigidBodyState`. The plate point's native work-conjugate transpose,
complete scaled sparse mass and shared bounded mass solver determine its
mobility; the body contributes its actual orientation, off-centre application
arm and full tensor inertia. The original `fs-contact::impact` restitution law
selects a single reciprocal impulse. No mass lumping or modal truncation occurs.

Preparation borrows the original mesh and full reference mass before support
elimination, checks every retained nonzero against the trajectory's reduced
mass, and solves once for point mobility. Retaining constrained rows is required
to recover actual dynamic support impulses: `M_full delta-v - J_full` includes
both constrained-row cross inertia and the applied point load. The represented
resultants are transverse Pz and bending Lx/Ly, including their external support
transfers; stationary supports do no work. A reduced mass alone cannot establish
those complete support reactions. In-plane momentum and drilling rotation
require their corresponding shell/solid owner and are not inferred here.

The event stages both original states, gates represented momentum changes,
actual restitution speed, linear/angular balance including supports, and combined
kinetic change plus the original nonnegative restitution loss. Polarized work
from both bodies is accumulated together in double-double arithmetic before
rounding, retaining loss against large common motion. Relative normal velocity
likewise combines original COM, rotational and plate motion before rounding;
a small approaching rotation cannot disappear into common translation and be
misclassified as a no-op. Native structural impulse
admission independently checks generalized transfer/work and recomputes the
post-event acceleration with the existing finite load, stiffness and damping.
No duration, heat deposition, pose repair or elastic-storage jump is invented.

Only plate velocity/acceleration and rigid momenta publish, after the final
cancellation poll. Physical time, plate displacement, all accepted telemetry and
rigid pose remain unchanged; both native smooth integrators continue from those
states. Telemetry is retained in place rather than copied. Point-location,
full-mass-nonzero and live-workspace caps bound preparation/event work; inactive
or separating contacts return without either mass solve or a state update.

The normal is the fixed +z direction of the original small-deflection plate.
Caller-supplied body/plate points must coincide in xy within an explicit
tolerance; their signed z separation enters the existing gap law unchanged.
This is one already-localized normal event, with stationary supports. Collision
search, moving supports, finite plate-normal rotation, simultaneous contacts,
friction and event adjoints remain outside this interface. Splitting the smooth
trajectory at the event is the caller's responsibility.

`tests/plate_impact.rs` independently compares consistent P1 plate inertia and
off-centre full-tensor rigid response with dense mass/tensor inverses at plastic,
partially elastic and elastic limits. It checks original damped and rigid
continuation, constrained-row support inertia and the fixed-point limit,
compensated physical loss, inactive/separating no-ops, mass/budget refusal and
final-poll atomic cancellation. The bridge advances the flexible-body coupling
gap in `frankensim-b8bxd.9.3`; it does not close its broader program.

## Simultaneous plate impacts against multiple finite rigid bodies

`fs_couple::plate_impact::bodies::PlateRigidImpactNetwork` extends the continuing
spatial plate adapter to several explicitly indexed finite rigid bodies. The
existing isolated and one-body simultaneous public APIs remain available.

Preparation takes the original `PlateDynamics`, mesh, complete reference mass,
`PlateImpactSite` slice, `PlateRigidImpactSetControls`, and
`PlateRigidImpactNetworkBudget`. It reuses the original point stencils, support
admission, consistent-mass solves, and cached `B M^-1 B^T` of
`PlateRigidImpactSet`. The cache is solved once per site. Each event receives a
mutable `PlateImpactBody` slice and one `PlateImpactBodyWitness` per prepared
site, declaring a participant index and a COM-relative material arm in that
body's axes. Equal poses or inertias never merge distinct participants.

The original plate contribution couples all sites, including sites on different
bodies. Rigid translational and full-tensor rotational cross terms are added
only when two sites name the same finite body. The resulting complete matrix
goes to the original `fs-contact` simultaneous normal-law owner. Initially
separating active rows remain in this law. `UniqueImpulses` and `UniqueVelocity`
determinations are returned unchanged; redundant impulses and support entries
remain representative loads, not unique local pressure or heat allocations.

An accepted nonzero event stages one original structural impulse and native
acceleration restart, plus one resultant native momentum impulse per affected
body. Point velocities are remeasured from the actual candidate native states
at every site, including geometrically inactive sites. Unrounded compensated
body resultants remain available through the per-body requested/actual transfer
checks. The shared audit reconstructs full constrained-row inertia reactions,
the represented transverse `Pz`, `Lx`, and `Ly`, each finite body's spin and
orbital increments, and collective kinetic work in one compensated accumulator.
The original native impulse-work bounds and separate physical loss/work bounds
remain authoritative. A failed numerical bound is not reported as dissipation.

Only after all original laws, native updates, endpoint checks, transfer balances,
and the final cancellation poll succeed are the plate velocity/acceleration and
affected rigid states published. Time, displacement, step history, and untouched
participants retain their original bits. No event clones an entire body array.
All-zero events are exact no-ops; cancellation or any early/late refusal leaves
every supplied participant unchanged.

Contact, mobility-product, mass, point-location, and workspace budgets retain
their original owners. `max_bodies` is a separate explicit positive bound;
checked workspace arithmetic reserves the actual storage size of all native
body receipts and compensated resultants before preparation allocates. Event
assembly is quadratic in contact count, with a linear number of original rigid
inverse-inertia actions. There is one plate impulse solve per event and no
per-body plate solve, sequential impulse sweep, or new numerical solver.

Independent tests compare a dense consistent plate mass and per-body full
tensor oracle against mixed same-body and distinct-body rows, all three elastic
regimes, contact/body permutations, translated geometry, and native trajectory
continuation with damping. Analytic shared-plate tests exercise active separating
and inactive rows. Supported and redundant fixtures independently reconstruct
full support reactions. Additional tests falsify sequential isolated impacts and
naive subtraction of large kinetic energies, preserve the strict native work
refusal at unresolved scales, and verify cancellation, bounded participants,
late imbalance refusal, retry equivalence, and one-body API compatibility.

This remains a fixed-configuration frictionless normal event for one flat plate
with its original +z normal and stationary supports. It does not perform collision
discovery, rigid-rigid contact, curved or finite-rotation shell kinematics,
frictional plate contact, moving-support work, or local thermal partition.

## Simultaneous normal impacts between independently meshed spatial plates

`plate_impact::plates::PlatePlateImpact` exchanges reciprocal normal momentum
between two continuing native `PlateDynamics` trajectories. The participants
can have different meshes, consistent masses, stationary support maps, damping,
and concrete load-owner types. `PlateImpactSurface` supplies each original
dynamics, mesh, complete uneliminated mass, and explicit flat-face reference
height. A common `PlateImpactSite` slice names world-x/y locations present on
both meshes; face offsets or material thickness are never inferred.

Positive contact impulse acts along world +z on the positive plate and -z on
the negative plate. At each site the original point observations give
`gap = height_positive + B_positive*q_positive - height_negative -
B_negative*q_negative` and the corresponding relative normal velocity. The
same original stencils transfer the transpose impulses. Both states must have
the same physical time. Their displacements, clocks, step counts and original
telemetry remain unchanged by this frozen-configuration event.

Preparation reuses both original `PlateRigidImpactSet` owners, including their
complete-mass/support admission and cached consistent-mass response solves.
The full contact response is the sum of their original
`B M^-1 B^T` matrices and is assembled once. No plate is lumped, rigidified,
modally truncated, or copied into a surrogate motion owner. The original
`fs_contact::multicontact` law resolves all rows together, including active
initially separating rows, and retains its global work, passivity, and
unique-impulse/unique-velocity determinations. Redundant coordinates remain
representative loads rather than unique local pressure or heat allocations.

A nonzero event accumulates opposite full nodal requests in double-double
arithmetic, then performs one original structural impulse and native
acceleration restart per participant. The shared plate-only audit reconstructs
constrained-row inertia reactions from each full mass, compares actual free-row
transfers with the unrounded requests, and retains plate/support resultants and
polarized kinetic increments unrounded until both participants enter the event
balance. This same extracted owner also serves the existing plate/rigid APIs.
Measured native endpoint speeds, original contact work, support-inclusive Pz,
Lx and Ly, total kinetic work, and the original native numerical work allowance
must all pass. Failed numerical work is not relabeled as physical dissipation.

Only a final successful cancellation poll permits publishing both velocity and
acceleration pairs. A failed load, solve, transfer, geometry, time, work or
balance check leaves both original states intact. All-zero contact solutions
invoke no native impulse/restart and preserve every original state bit.

`PlatePlateImpactBudget` retains both original point/mass/contact/workspace
budgets and adds checked combined live-storage admission for the two owner
ceilings, compensated requests and shared contact arrays. The original
64-contact law bound remains. Cached response preparation is quadratic in
contact count; event transpose/observation work is linear in contact count
times the original state sizes, with one native solve/restart per plate and
no per-contact trajectory clone. Contact observations and transpose accumulation
poll cancellation at bounded array tiles.

Six G0/G1/G3/G4 groups independently compare a coarse triangle and a four-triangle
refinement with exact P1 area integration and dense mass inversion across
inelastic, partially elastic and elastic laws, contact permutations and world
translations. They continue both original damped trajectories; prove active
separating/inactive cross response and falsify sequential pair sweeps; recover
full fixed-slope/free-translation inertia reactions and redundant-load motion;
retain a known tiny physical loss during large common translation where total
energy subtraction fails; refuse unresolvable nodal momentum updates and
unsynchronized times; and verify distinct load types, second-owner load refusal,
native acceleration restart, budget/cancellation boundaries and atomic retry.

This is a frictionless normal event for two flat plates in their shared fixed
world frame. It does not detect collision times, certify finite-thickness
nonintersection, rotate normals, evolve sustained or tangential contact,
advance moving supports, integrate a network of more than two plates, or assign
the admitted mechanical loss to heat or other material destinations.

## Simultaneous plate and finite rigid-body normal impacts

`plate_impact::multiple::PlateRigidImpactSet` resolves an entire localized normal
contact set between the same original spatial plate and one finite rigid body.
Preparation reuses the isolated owner's full/reduced mass admission and solves
the original consistent mass once per point. It caches all `B M^-1 B^T` cross
responses, retaining both independently solved off-diagonal entries for the
native contact owner's symmetry admission. There is no modal truncation, mass
lumping, or dense matrix in the plate's full degrees of freedom. The contact
matrix has at most 64 rows and columns.

Each fixed-pose event adds the shared rigid body's full tensor contribution
`1/m + (r_i cross n)^T I^-1 (r_j cross n)` in its original body frame. The existing
`fs_contact::multicontact::solve_normal_impacts` owns the coupled complementarity
and declared per-site restitution/proximity laws. Geometrically active sites
participate even if initially separating: another contact can make them close.
Inactive sites retain their complete predicted and measured cross response.
Mutually incompatible rebound targets that violate assembled passivity refuse;
neither sequential pair impulses nor a hidden restitution reduction is used.

All point loads accumulate through the original work-conjugate transpose into
one native structural impulse/restart. Their original central moments accumulate
into one native rigid resultant impulse. The same shared admission as the
isolated path checks actual stored transfer, complete constrained-row support
reactions, represented transverse linear/bending-angular momentum, and combined
compensated kinetic work. Every original point's post-event relative velocity
is checked against the coupled law. Both original states publish only after the
last cancellation check, preserving displacement/pose, physical time, accepted
step counts and telemetry. An all-zero contact solution leaves both inputs
unchanged and invokes neither native impulse update.

The original contact determinacy is retained. `UniqueVelocity` identifies a
redundant PSD contact set: its impulse vector and local support entries describe
one checked representative, not unique local pressure, work or heat allocations.
The original range admission is required for total loss. No thermal partition
is inferred from the normal impulses. Numerical passivity residuals remain
separate from nonnegative physical loss.

Point-count, mass-nonzero, mobility-product and live-workspace caps are checked
before preparation. The adapter's workspace bound includes cached rows, native
mass-solve scratch and event arrays; the separately bounded original contact
solver's scratch, borrowed model data and existing trajectory history are
excluded. Repeated events reuse the prepared plate cross mobility and perform
one structural impulse solve/restart irrespective of contact count.

`tests/plate_impact/multiple.rs` compares original plate/tensor response with an
independent dense coupled solve and a falsifying independent-scalar treatment.
It checks contact-order covariance, plastic/partial/elastic impacts, continued
damped plate and rigid motion, initially separating contact activation, complete
support cross inertia, redundant-load determinacy, exact no-ops, budget and
penetration refusals, and final-poll cancellation followed by a valid retry.
The five isolated impact tests retain their original behavior after the shared
preparation/transfer extraction.

The represented model remains the original flat plate's fixed +z normal and
transverse Pz/Lx/Ly dynamics with stationary supports. Collision discovery,
moving sites/supports, finite normal rotation, friction, shell in-plane dynamics,
and event adjoints require their corresponding owners. This advances the actual
flexible-body transfer gap in `frankensim-b8bxd.9.3` without closing that broader
program.

## Rigid impacts against solids with material memory

`inelastic_impact::{plastic_rigid_impact_candidate,visco_rigid_impact_candidate}`
compose one frictionless rigid/material-point impact using the original
`fs-contact` Newton restitution law. The native finite plastic or viscoelastic
owner supplies its current barycentric witness/velocity, supported directional
mobility and work-conjugate impulse distribution. Graded reference density and
the full consistent mass remain active. The rigid owner supplies its current
material-point kinematics, complete body inertia tensor, effective mass and
actual represented linear/angular momentum jump.

The normal points from solid to rigid, with equal and opposite impulses through
the original witnesses. Negative gaps refuse, even inside the activation
tolerance. Off-normal witness separation refuses rather than relocating a
contact point. Both candidates retain configuration, source bindings, material
history and physical time. A failed geometry, solver, physical balance or
cancellation gate leaves both accepted states intact.

The combined report checks actual restitution speed, total linear and angular
momentum about the requested origin, body-local defects, conjugate contact work
and kinetic change plus contact loss minus signed moving-support work. Contact
loss is reported separately and is not deposited as heat. Independent native
regressions in `tests/inelastic_impact.rs` exercise graded edge-point mobility,
full off-diagonal rigid inertia, moving supports, elastic/inelastic restitution,
and retained plastic/Maxwell histories. This is a bounded single-contact event;
collision detection, event-time localization, simultaneous manifolds, friction
and persistent-contact evolution retain their existing owners.

`mixed_visco_rigid_impact_candidate` connects the same impact law and coupled
momentum/work audit to finite-compliance P2/P1 viscoelastic motion. Its
`InelasticRigidContact<Tet10ImpulsePoint>` uses all ten signed quadratic shape
weights for the curved witness, velocity, supported mobility and opposite
impulse. The default contact payload remains the original P1 material point.
Continuous pressure, all Maxwell histories, physical storage and time survive
unchanged. Exact incompressibility refuses this mass-only contact path; the
pressure-aware `multiple::volume` entries below supply its constrained mobility.
`tests/inelastic_impact/mixed.rs` couples a curved, relaxed quadratic solid to
off-diagonal rigid inertia. Independent P2 reference integrals determine its
restitution impulse, moving-support work and kinetic change; the existing audit
checks actual linear/angular momentum and conjugate work at both witnesses.

`inelastic_impact::multiple` extends this to simultaneous circular Coulomb
contact for native plastic, viscoelastic, thermal Tet4 and thermal Hex8 owners.
It shares the existing mixed rigid/solid mechanical core: original paired
witnesses and stable IDs build the complete solid cross-point mobility and
signed shared-rigid response, including the full inertia tensor. One existing
`fs-contact` solve determines the coupled impulses; one summed native update
per body retains their cross-contact work. Material history is borrowed while
building mobility and copied only into the complete candidate. Native source
errors, mass-solve reports and the full `SolidImpulseReport` remain available.

The paired witnesses, caller contact/row/workspace budgets, original contact
determinacy and endpoint momentum/work gates remain binding. Held solid DOFs
retain their accepted velocity. The combined report separates external
`support_work_j` and checks kinetic change minus this work against contact
dissipation. Velocity-only determinacy with moving solid supports returns
`UnresolvedSupportWork` with the original contact receipt before staging an
update; it does not choose an unestablished support-work partition. Stationary
supports and rigid-only redundant manifolds retain their existing admission.
Temperature, density, plastic/Maxwell/quadrature history, storage and clock stay
unchanged. Contact loss remains a mechanical account, without inferred heat.

The existing inelastic regression fixtures now exercise coupled two-point
sticking and sliding, independent graded-mass/full-tensor cross terms, contact
order invariance, moving-support restitution/work, retained histories and
native budget/cancellation refusals. This fixed-pose event does not supply
collision detection, event localization or persistent-contact time integration.

`multiple::mixed_visco_mixed_contact_candidate` and
`multiple::elliptic::mixed_visco_mixed_elliptic_contact_candidate` extend this
same coupled owner to finite-compliance P2/P1 Maxwell solids. The complete
quadratic cross-point matrix supplies circular or directional Coulomb response;
one summed native event preserves physical pressure and all material histories.
The original ellipse owner retains exact free slip on a zero-coefficient axis
and exact circular equivalence for equal principal coefficients. No principal
direction or contact point is solved independently.

The curved quadratic regression adds two simultaneous impacts with independently
integrated P2 mass and full rigid tensors. Its exact normal block predicts both
positive impulses, signed moving-support work and kinetic loss. It also covers
the directional adapter's circular equivalence and a free-slip principal axis.
These Rust regressions remain unexecuted under the existing native dependency/
resource blocker; rational checks and formatting do not substitute for Cargo.

`inelastic_impact::multiple::volume` connects exactly incompressible Maxwell
motion and all three conservative quadratic formulations to this same coupled
contact owner. Circular and directional Coulomb entries retain every signed
P2 cross-point/cross-axis response and the full rigid inertia. Exact mode
assembles and factors the original mass/pressure KKT once per event, then reuses
it for all homogeneous mobility columns and the final summed impulse. Accepted
velocities must already meet the original weak-volume-rate allowance. The
conservative finite-compliance modes retain their original mass-only response.

`VolumeContactControls` carries the native mixed pressure scales, sparse budgets
and balance controls. The returned `MixedImpulseReport` retains physical
external impulse, internal pressure impulse, actual rounded momentum, signed
support work and weak volume rates separately. Pressure resultant, torque and
work defects are added in magnitude to the coupled absolute residuals; they
cannot cancel a body's native endpoint error. Pa s pressure impulses never
replace stored Pa pressure. Maxwell histories remain unchanged, and conservative
event work enters the original compensated cumulative-energy account without
resetting its initial energy. Contact loss is not assigned to material heat.

`tests/inelastic_impact/volume.rs` exercises two actual boundary witnesses against
an independent rational 30-displacement/4-pressure elimination and a two-contact
restitution solve. It checks pressure-induced cross-axis response, paired
momentum/work, nonzero pressure impulses with zero pressure work, both friction
adapters, retained pressure/history and conservative exact/full-material
dispatch. Native Rust execution remains pending the established build-resource
block; these are regression definitions backed by separate algebraic checks.

`solid_enthalpy::finite_plastic::FiniteThermoplasticSolid::mixed_contact_candidate`
binds the thermal Tet4 event to the original source-backed caloric trajectory.
The accepted clock, cell enthalpies, reference geometry/mass and exact caloric
objects must match before the native event is proposed. The returned source
owner retains its cells, source families and clock while adopting the new
mechanical velocity; the borrowed thermal simulation remains unchanged. The
original coupled `advance` then continues from that complete accepted state.
Impact loss stays separately reported, without an invented heat partition.
The existing finite-plastic enthalpy seam checks a hot, plastically deformed
equal-mass collision, exact source/caloric retention, stale and independently
advanced heat refusal, and continuation on the original heat timeline with
the original consistent-mass inertial load.

`inelastic_impact::surface` prepares and queries the original deformed Tet4
boundary for native plastic, viscoelastic and thermoplastic owners. It shares
the existing closed-surface preparation, reference identity check, BVH refit,
orientation admission and canonical face/edge/vertex witnesses. Queries borrow
material history; the same barycentric witnesses feed one simultaneous native
contact candidate with original graded mass, full rigid inertia and support
reactions. Penetration and reference mismatch refuse without repairing geometry.
This discrete adapter makes no swept-path claim and does not approximate Hex8
or curved P2 boundaries with Tet4 faces.

`inelastic_impact::surface::scene` extends the same native Tet4 queries to the
complete prepared sphere scene. Plastic, Maxwell and thermoplastic bodies use
the original sphere-pair broad phase, explicit pair laws and disjoint canonical
surface/pair IDs. A shared scene composition changes only the solid point
payload to its unchanged native tetrahedron and barycentric coordinates. All
rigid/solid and rigid/rigid rows can therefore enter one original native mixed
solve, retaining shared-body cross mobility and reciprocal momentum transfer.
Pair penetration, scene mismatch, reference mismatch and original budgets
refuse before any state changes. The existing relaxed-Maxwell fixture adds an
independent unequal-mass solid/two-sphere chain, separate pair/surface restitution,
contact permutation and refusal checks. This query is instantaneous; its native
consumer owns the subsequent trajectory and any thermal allocation.

`inelastic_impact::directional` samples the original interface sources for all
four native manifold owners before invoking the shared elliptic mechanical
solve. One original history and explicit static/kinetic coefficient role bind
each input contact. Active contacts sample the actual complete incoming tangent
speed, including any principal free-slip direction; inactive rows do not query
sources. Canonical contact IDs retain their original input indices. Local
temperature, pressure and other source conditions come from the caller's typed
callback rather than from an invented impulse duration or pressure. A source
failure, mechanical refusal or cancellation publishes neither mechanical state
nor sampled interface state. Source claims and accepted slip histories survive
unchanged, and physical loss remains separate from any thermal allocation.

The source-caloric `FiniteThermoplasticSolid` also provides
`mixed_elliptic_contact_candidate` and `mixed_source_elliptic_contact_candidate`.
Both admit the original caloric objects, enthalpies and clock before contact
or source evaluation, then retain the complete source-bound mechanical owner.
They share the circular entry's admission/adoption helpers and continue through
the same coupled thermal advance. They do not rebind material sources or heat
the body from contact loss without a declared partition.

The existing native fixtures now exercise current deformed-surface geometry,
independent graded-mass/full-tensor impact response, source speed and explicit
coefficient roles, exact principal free slip, input permutation, late source
failure and retained plastic/Maxwell/thermal histories. The original caloric
seam adds directional parity and source admission ordering. Formatting and
source/API checks pass; native test execution remains unverified while the
required constellation sibling checkouts are unavailable.

`inelastic_impact::trajectory::visco_sphere_trajectory_candidate` continues the
native fixed-temperature Maxwell body and COM-centred rigid spheres through
transient impacts. It composes the original implicit material advance, full
consistent inertia, simultaneous contact solve and rigid endpoint-kick/free-spin
integrator. Stationary held DOFs retain their accepted displacement. Dead solid
loads, rigid wrenches and rigid gravity remain explicit external sources.

The existing interval sweeps enclose the declared affine interpolation of each
trial's native endpoint coordinates and preserve tetrahedron orientation. An
unresolved path shortens the interval and recomputes the original material
advance from the same accepted history. Rejected trials add no time, relaxation
loss or contact impulse to the accepted account. Accepted endpoint geometry and
contacts carry into the next interval without repeating its BVH query. A positive
activation band localizes events without witness relocation; contact first
entering that band at the requested final endpoint resolves on continuation.
Accepted endpoints remain nonpenetrating.

Each report retains exact accepted durations, native material/rigid receipts and
path bounds. Physical Maxwell relaxation and contact loss remain separate from
backward-Euler kinetic loss, signed material endpoint error and rigid drift.
The native actual angular-momentum increment and support/external torque use the
original consistent mass and current coordinates. Backward-Euler angular drift
is reported and bounded independently, never booked as a physical torque. Final
balances and sums of individual residual magnitudes must both pass the caller's
whole-call tolerances, preventing cancellation of large numerical defects.

Every failure preserves all input owners. This bounded isothermal, first-order
trajectory does not certify the unknown continuum path, provide sustained
contact, or assign contact heat. Offset spheres, unsupported sphere-pair
encounters, unresolved grazing/Zeno behavior and exhausted budgets refuse.
G1/G3 cases specify the analytic supported bounce (impulse 3 N s and loss
0.75 J), accepted-duration Maxwell replay, free recoil/deformation, native
continuation and late rollback. These new tests share the execution limitation
stated above; formatting and bounded source review do not establish runtime
convergence or a performance benchmark.

`inelastic_impact::trajectory::sustained::visco_sustained_sphere_trajectory_candidate`
continues loaded resting/sliding contact through the original native Maxwell
body and a complete prepared sphere scene. The native accepted Piola field and
held loads supply a force kick, one simultaneous native contact solve couples
every active solid/sphere and sphere/sphere row, and the original constitutive
endpoint map continues Maxwell history during kinetic drift. Rigid kick/drift
arithmetic and phase gates are shared with the existing elastic trajectory.
Newton targets come from accepted PRE-FORCE normal velocity, so a held load at
rest does not acquire a fresh restitution bounce at every interval.

Each shortened trial recomputes all phases from the same accepted material,
velocity, geometry and time. Rejected trials contribute no accepted relaxation,
contact loss, source transfer or clock advance. Shared interval sweeps admit
the complete declared affine numerical path, including every sphere pair and
tetrahedron orientation. Endpoint queries carry forward without repeating the
BVH work. Exact accepted durations and original canonical contact keys remain
in the report. Final-endpoint activation resolves on the next continuation.

The whole account distinguishes native physical relaxation, unallocated contact
map loss, source displacement work, signed support work, material endpoint error,
solid force/source splitting and rigid force/drift defects. In particular, the
normal loss produced by supporting a resting load in this first-order split is
not automatically physical heat. It accompanies an explicitly bounded force
splitting term and must decrease under timestep refinement. Each trial receives
its duration fraction of the whole-call absolute splitting allowance; native
phase momentum/work gates and accumulated residual-magnitude gates also apply.
Opposing errors cannot cancel their way through admission.

The scope is fixed-temperature native Tet4 material with stationary held DOFs
and COM-centred spheres. This explicit first-order split is not a general
Moreau-Jean existence theorem or an exact-flow collision certificate. Source
friction histories, moving supports and contact heat partition require their
corresponding continuing owners. Focused regressions specify a two-sphere loaded
chain with restitution one and no spurious bounce, independent contact impulses
and support torque, halved numerical support loss under timestep halving, exact
accepted-duration Maxwell replay and late whole-call refusal. Native execution
remains unverified while the required constellation siblings are unavailable.

## Physical pressure from an explicitly selected lumped duct

`acoustic_realize::duct_observer::BaffledDuctObserver` admits the existing
`BoreRealization::LumpedLadder` for an unbranched flanged cylinder/frustum bore,
including a single widening cone, with one coherent gas per authored axial
segment, an independent pressure-matched exterior, and either its original
numerical wall law or source-bound local material walls. Both massive and quasistatic reeds retain their implicit pressure/work
junction. The original `fs-phs` terminal-flow row observes the accepted state;
that flow enters the existing finite-aperture Rayleigh receiver and causal
atmospheric filter. Internal source pressure is not added to physical microphone
pressure.

The explicitly selected lumped path also admits original tone-hole branches
with one uniform interior/exterior gas and rigid walls. Each open or vented
mouth uses its original assembler-selected flow row: the compact neck flux or
the long chimney's actual terminal flux after its compression cells. Native
effective aperture already includes open fraction. Closed pads keep their
acoustic load and emit no receiver flow. Regional tap materials, a distinct
exterior and chimney wall/source partitions retain explicit refusals.

The existing Rayleigh receiver and causal atmosphere filter are prepared
individually for every active mouth. They process its actual flow and physical
range before pressures are summed. `BaffledDuctObservation::openings` retains
each original axial station, effective radius, translated receiver, filter and
delay support. Geometry uses the existing parallel-baffle convention: every
mouth faces +z, with receiver.x shifted by bore_length - opening_station.
No hole azimuth, arbitrary normal or mutual radiation is inferred. Axial and
chimney propagation stay in the native state, without a second propagation
delay, source-pressure proxy or radiation load. The original authored stations
remain independent of the ladder's numerical cell-length correction.

The physical image selects the original mass-conserving frustum owner before
passive-pitch admission. Its prepared chain carries the original `RegionalDuct`
reference, is borrowed by the passive impedance and pitch gate, then moves
unchanged into reed playback with its own outlet-flow row and receiver.
Qualification and rendering use the exact same model, without a second cold
assembly or a spherical-to-frustum substitution after admission. Spatial
screening, measured pitch tolerance, original source-band admission and complete
material state remain binding. Legacy source-pressure comparisons retain their
automatic spherical-cone behavior.

Each region's density, sound speed and transport determine its exact frustum
compression/inertia and local viscothermal memory. Numerical walls keep their
original owner; source walls retain each section's original material and full
local constitutive law. The reed and Bernoulli source use the actual inlet gas.
The exterior determines only terminal radiation loading, Rayleigh pressure and
atmospheric propagation. Passive admission uses the original regional TMM with
the same interface and exterior assignment; its scan follows the original
regional one-way travel time. Frozen regions require coherent gas states and
equal static pressure; mean flow, mixing and thermal equilibration are not
modeled. Identical regional profiles retain the original uniform allocation,
scan arithmetic, physical output and receipts bit for bit.

`MaterialWallLayout::Uniform` and `Segments` bind their original
`LocalMaterialWall` sources to this same spatial chain. The canonical frustum
query supplies authored region, midpoint radius and actual slanted strip area;
axial gas geometry remains separate. The existing material reduction supplies
equilibrium mass/stiffness, Kelvin resistance and every positive Maxwell arm,
with original source errors retained. The existing cell junction and Maxwell
storage carry wall pressure, flow and hereditary work together. The same
material partition enters the original regional TMM before passive admission,
and its source layout remains in the physical observation result.

This is frozen small-displacement local membrane motion. It adds no
static-pressure prestress, axial shell interaction, mean thermal/material
evolution or exterior wall radiation. Every TMM query retains material frequency
support, and cell geometry retains the original thin-wall admission. The
existing exact propagator keeps its frequency gate; graph midpoint-memory
limits are not applied to that different original propagator. Material
applicability is not an antialiasing filter or a broadband excitation certificate.

Nonuniform profiles allocate exactly the requested total with at least two cells
per authored region, refining the greatest local cell travel time first with
deterministic authored-order ties. Their reported spatial screening estimate
uses the largest cell travel time relative to total regional travel time. The
minimum count is obtained by independent monotone integer searches without
building trial chains. This estimate does not certify arbitrary material-interface
eigenmode error: the measured regional TMM pitch gate remains decisive.

The original compact/Foster radiation-load approximation remains. Rayleigh
observation does not certify its continuum loading accuracy. Axial propagation
stays in the spatial ODE, so the receipt reports zero separately factored
axial-delay samples; finite spatial resolution does not claim an exact continuum
wavefront. Receiver geometry, exterior delay, quadrature band and all media remain
explicit. Side openings retain the uniform-gas/rigid-wall scope above; separate
body receivers require their own supported owner and refuse on this path.

Authored G1/G3 regressions cover retarded inverse-range pressure for massive and
quasistatic reeds on cylinders and single cones, absolute-pressure substitution,
warm-bore/cool-room assignment, actual regional cylinder/frustum microphone
output, identical-profile output/receipt bits, unequal-pressure and insufficient
regional-cell refusal, exact cell allocation and minimum screening work,
independent original-TMM composition across a gas interface, the uniform
exterior-reference bits, and a cone tolerance that rejects the physical frustum
even when the legacy psi image would pass. Test presence does not claim native
execution.


Authored source-wall regressions additionally compare a sourced Kelvin cylinder
with the equivalent original numerical-wall pressure and pitch bits; exercise
regional Maxwell wall changes in actual frustum microphone pressure for massive
and quasistatic reeds; retain gas, exterior and source-wall assignments; and
check the original out-of-band source refusal and duplicate-wall rejection.
These are synthetic constitutive checks, not measured-material or instrument
validation. Their presence does not claim native execution.

## Superfluid mutual-friction heat on the original thermal timeline

The opt-in `superfluid::advance_closed_friction_heat` composes the original
closed `fs_duct::superfluid::TwoFluidChannel` with the original
`fs_thermal::enthalpy::DepositionSimulation`. Both clocks must match. Ordered,
unique, active `CellId` destinations route each new interval's dissipated
joules exactly once; unrelated thermal cells retain their original contacts.

The channel stages its wave endpoint and physical friction heat privately.
Native thermal sources, contacts and events advance to that endpoint, then
native `apply_heat_j` applies the new joules using the unchanged borrowed
caloric laws and their material-domain checks. Caller-selected local,
aggregate-transfer and combined-storage tolerances admit the represented
energy increments. The original wave, thermal-march and heat-intake reports
remain separate: the friction deposit is internal to the composed system,
and the thermal march's supplied sources are not counted a second time.
The combined residual uses the owners' original cell-relative residuals,
preserving small heat increments against large thermal energy datums.

Thermal refusal or cancellation, including the final checkpoint after the
thermal candidate is prepared, preserves both accepted histories. An
infallible assignment publishes the thermal candidate only after channel
acceptance; retrying deposits no rejected or cumulative heat.

This is a first-order endpoint heat split with a frozen wave reference.
Finite thermal temperatures evolve through the caller's original caloric
law; they do not relinearize the channel or replace its wave temperature
perturbation. Geometry correspondence and heat destinations are explicitly
authored. Closed ends exclude unhandled mass and entropy imports. No open-flow
thermal transport or interface law is inferred. For a thermal law storing
enthalpy, the account remains represented wave-plus-thermal storage rather
than an unconditional total-internal-energy theorem.

## Impulses on continuing coupled modal systems

`CoupledModalSystem::apply_modal_impulse` applies a simultaneous generalized
momentum jump through the original mass-normalized mechanical coordinates.
The same physical port columns map N s into modal momentum as map N into
modal force. Finite bilateral spring/dashpot forces have no instantaneous
impulse. Mechanical displacement, Maxwell internal strains, proper passive
radiation states, prepared operators, sample clock and preceding interval
frame/port records remain unchanged. Only mechanical velocities publish.

Each represented coordinate increment must satisfy the original relative
solve tolerance against its requested/represented impulse, independently of
background velocity and other coordinates. Stable kinetic work and prescribed
impulse work are accumulated separately; their defect subtracts momentum
before multiplication and uses the existing energy tolerances on the jump.
Original component and whole-system energy, pressure and endpoint connection
force ceilings remain binding. The returned receipt describes external work,
without invented heat, damping or elapsed time. Native staging uses reserved
scratch and publishes all components only after admission and cancellation.

The existing `ImpulseRenderer` now routes bilateral `CoupledModal` voices
through that native transaction and keeps same-time mixed-voice publication
atomic. `from_coupled_modal_forces_and_impulses` shares the original physical
force projector, preload and scheduler. All component impulses at a sample
share one network jump; held forces and boundary-event ownership survive.
Independent voices also refuse a nonzero impulse lost entirely to velocity
rounding and retain small signed work between large cancelling modal works.
The direct native API polls during the jump and just before publication;
gated rendering retains its original whole-callback cancellation boundaries.
Nonlinear contact/friction, rigid constraints, inertance feedthrough and an
empirical strike law are outside this prescribed linear-network reset.
`tests/coupled_impulse.rs` covers analytic two-mass continuation at two sample
rates, ordinary/prepared paths, retained material/acoustic history, reciprocal
momentum and work, refused/cancelled jumps, mixed-voice atomicity and callback
partition invariance.

## Native rigid-body acoustic histories

`pcm_wav::observation::moving::rigid` observes original `fs-mbd` prepared
rigid-point motion at explicitly supplied physical timestamps. Body-fixed source
and ideal microphone attachments retain native position, velocity and full-inertia
acceleration, including tangential and centripetal terms. The mechanical world
frame must be the fluid rest frame. Observation does not advance either owner.

`BodyFrameCompactSourceStrength` carries physical volume injection and force ON
THE FLUID. Force transport uses `F_world = R F_body` and
`Fdot_world = R(Fdot_body + omega_body cross F_body)`; the supplied derivative is
of body-coordinate components. Scalars remain unchanged. Total body wrench,
gravity, support loads and net torque do not infer acoustic source strengths.
Declared source bandwidth must include the exported force's rotational variation.

`RigidAcousticHistory` reserves bounded paired motion/strength arrays once and
commits both only after validation and the final cancellation poll. Append neither
allocates nor advances mechanics; capacity, source, clock and cancellation refusals
leave both arrays unchanged. Borrowed accepted slices feed the original C2
trajectory, moving-microphone and coherent physical-pressure renderer owners.
The caller retains body identity, physical clock, smooth load resolution and
finite-source validity; this does not smooth impacts, model deformation or supply
radiation reaction. Native accepted-rotor and attached-microphone regressions
exercise the complete pressure path and independent physical reference limits.

## Moving microphones on accepted physical trajectories

`pcm_wav::observation::moving::MovingCompactRenderer::with_moving_receiver`
observes physical compact volume-flow/force histories through borrowed
`fs-aeroac::curle3d::moving::trajectory::PreparedMovingReceiverSource` pairs.
Every source must share the same original microphone trajectory owner. Each
pair retains its source path, original fluid snapshot, finite receiver window,
continuous retarded-clearance admission and root-work limits. Source-strength
interpolation still evaluates values and derivatives at the same returned
emission time. Neither mechanical trajectory is advanced by observation.

The existing renderer owns one global receiver clock, deterministic coherent
pressure accumulation and preallocated transactional callback scratch. Microphone
motion changes the retarded source clock and the whole-window Doppler bandwidth
bound; it does not replace source velocity in the pressure law or multiply the
Eulerian pressure by an extra receiver-motion gain. The associated receiver type
keeps fixed microphones and moving trajectory owners distinct. Complete window
and residual-policy admission occurs before exact-endpoint reconciliation.

Pressure remains in Pa through the existing `PressureRenderer`, causal
`DecimatedRenderer` and explicitly scaled PCM writer. Missing source/receiver
history, insufficient bandwidth, cancellation and numerical refusal publish no
partial callback. Hermite interpolation is not an antialias filter; the producer
still declares source bandwidth and geometric amplitude modulation. This is an
ideal scalar microphone in the original homogeneous quiescent lossless medium,
with no microphone-body scattering, wind/self-noise, feedback or experimental
source qualification. Five authored G1/G3/G4 regression groups cover analytic
receiver motion, independently inverted curved pairs, coherent source sums,
partition invariance, finite domains, cancellation and the original PCM path.
Native compilation and execution remain separate verification obligations.

## Fixed laboratory histories in a uniform fluid flow

`MovingCompactRenderer::with_uniform_flow` observes original accepted volume
and force histories through `UniformFlowCompactObserver` paths bound to one
`UniformFlowReceiver`. Every source must retain the same fluid, laboratory
velocity, microphone and limits. The existing Hermite source interpolant,
finite-history and bandwidth admission, coherent pressure sum, resource limits,
global receiver clock and transactional callbacks remain the only render path.

Uniform flow changes the physical pressure coefficients and propagation delay;
fixed laboratory endpoints have observed `d tau/d t = 1`. The fluid-frame
Eulerian `1/kappa` is never used to inflate their required receiver bandwidth.
Exact retained endpoints use cached relative geometry and are accepted only
when the rounded arrival equals the requested timestamp and the original
residual budget holds. Large epochs do not cause the renderer to subtract two
large Galilean positions. Explicit `observer_at` queries do return fluid-frame
coordinates and may refuse an unrepresentable transform.

Pressure remains in Pa for the original observation/decimation/PCM interfaces.
This supplies homogeneous lossless convective propagation, not wind-generated
microphone noise, shear/refraction, scattering, absorption or source feedback.
Four focused regression groups cover coherent cubic histories without Doppler
inflation, common-receiver and history/bandwidth refusals, exact endpoints at
large epochs, and bitwise callback partition/cancellation/retry behavior.

## Stationary compact fields from accepted force and volume histories

`pcm_wav::observation::moving::stationary::StationaryCompactRenderer` reuses
the original finite-history admission, physical receiver clock, source cubic
Hermite coefficients and retarded-time endpoint reconciliation. For fixed
sources and a fixed receiver it observes pressure, all Cartesian particle
velocities, signed instantaneous intensity and the original `FirstOrder`
microphone's Pa-equivalent output. Source pressure and velocity are summed
coherently before intensity or microphone projection; a pressure node can
still contain transverse motion and a directional-microphone signal.

Each source supplies its force impulse at its FIRST retained endpoint. The
observer integrates the exact same force cubics analytically into cached
double-double prefix impulses and evaluates their partial-interval integrals
at emission time. It does not reset history at callbacks, integrate on the
audio clock, infer missing prehistory, advance mechanics or peak-normalize
the signal. Both the monopole `Q,Qdot` terms and the complete compact-force
`impulse,F,Fdot` terms retain their induction and near fields.

The original work, history, bandwidth and payload limits also charge the
added transfer records, impulse endpoints and full-field scratch. Callbacks
allocate nothing and publish atomically; full-field and `PressureRenderer`
callbacks share one finite cursor. The ordinary decimator and explicitly
scaled PCM writer can consume microphone output without changing the physical
source. G1/G3/G4 regression cases cover an independent polynomial primitive,
scalar-pressure parity, explicit initial history, coherent reinforcement and
cancellation, quadratic flux, callback partitioning, PCM output and refusal.

`StationaryCompactRenderer::with_ideal_reflections` extends those accepted
histories to the complete image group of an original `IdealReflectionObserver`:
one to three perpendicular rigid or pressure-release planes. Every image keeps
its own retarded time and reads the same physical cubic and force-impulse
prefix. Cached virtual receivers retain the original geometry. Their complete
particle velocities are mapped back through the transpose image Jacobian and
wall parity before coherent Cartesian field summation. Intensity and ideal
directional microphone output are formed only after all source/image fields
combine; a pressure-release node can retain normal motion and a dipole signal.

All-path finite-history coverage and physical half-spaces are admitted before
publication. Work charges sources times paths times output samples; retained
payload charges every path transfer while source histories and impulse
prefixes remain shared. Existing transactional field/microphone callbacks and
global sample clocks are retained. Regressions in
`tests/moving_pressure_history/stationary.rs` compare rotated mixed corners
against independent source images and analytic polynomial primitives, test
wall boundary conditions and directional pressure nodes, and exercise path
history/work refusals and callback partition/cancellation.

This remains a one-way stationary compact-source model in a homogeneous
lossless linear fluid. It refuses moving sources and inherits each source
field owner's representable range, including individual flux. It supplies no
source qualification, nonlinear propagation, finite-impedance or general
boundary scattering, microphone-body response, radiation reaction or
experimental amplitude validation. Native Cargo execution is a separate
obligation from the authored regressions.

## Full spatial circuits coupled to complete heat assemblies

`conduction::heat::step_spatial_conduction_circuit_candidate` and the
`PreparedHeatNetwork` method of the same name connect the native multi-electrode
circuit descriptor, fixed capacitors/inductors, and original material graph
resistors to a complete `HeatAssembly` at one finite-enthalpy endpoint. The
shared spatial trial driver retains the original material, electrical-contact,
thermal-interface, layer, radiation, bath, and reciprocal finite-link owners.
Every trial starts from the accepted charge and flux state.

Bulk, contact, and graph-resistor Joule work enter their declared destinations
once; whole-cell graph-resistor receivers must be disjoint from spatial field
volumes. The assembly's original `closed_balance` accounts for enclosure
radiation internally and distinguishes external bath work. The returned sample
retains the full assembly state.

Source, geometry, domain, budget, or cancellation failure leaves both accepted
histories unchanged. This route retains the existing fixed-geometry,
field-independent thermal-boundary, confined-phase, and fixed authored
capacitor/inductor restrictions. Graph series contacts and source-dependent
electrical storage remain unsupported and refuse explicitly. The first-order
coupled update requires temporal and spatial refinement for accuracy.

## Planar pressure and viscous force/torque driving native rigid motion

`stokes_motion::planar::PlanarStokesMotion`, under the existing
`experimental-stokes-motion` feature, advances both translations and rotation
through the original material BDF1 Stokes and `fs-mbd` owners. A loaded base
solve and three homogeneous rigid unit motions supply the complete 3-by-3
fluid impedance, including translation/rotation cross terms. A reusable native
`fs-la::LuWorkspace` solves the mass/inertia-scaled coupled equation. A fifth
native flow solve uses the actual represented rigid endpoint rates. No fitted
drag, prescribed body trajectory or assumed rotational added inertia is used.

Depth converts original pressure/viscous force and moment to whole-body N and
N m; pressure datums remain explicit. Moments and rigid boundary velocities
share the reference COM. Native mass properties admit a COM pose only; this
lane additionally requires a world-aligned principal body z axis and planar
momenta. It retains both BDM normal moments, including rotation's nonzero
linear moment. Native quaternion drift supplies the accepted orientation.

Actual represented linear/angular momentum residuals have separate force and
torque tolerances. Spring storage uses represented pose coordinates; physical
dashpot loss and numerical backward-Euler loss remain distinct. Body work and
fluid kinetic change and bulk viscous heating are separate accounts, not a
complete fluid/body energy or passivity claim. Body, fluid history, translation,
angle and time publish
atomically after endpoint gates and cancellation checks. Transport is frozen
once; five-solve source/boundary callbacks must be pure and fixed.

Reference mesh geometry and torque arms remain fixed. Independently authored
small-translation and small-angle bounds are enforced; there is no ALE,
transported volume, convection or finite-motion claim. Every rigid unit mode
must satisfy the native boundary/flux admission with the exterior owners.
Four authored regression groups cover independent pressure-driven translation
and torque with a matching affine inertial body force, two depths, analytic
off-center reciprocal cross impedances, viscosity response, native pose/work,
rotation traces, axis admission and atomic late-cancellation continuation.
Native compilation/test execution remains a separate verification obligation.

## Finite guided translation from the original nonlinear moving-domain fluid

`stokes_motion::translating::GuidedTranslationMotion`, under the existing
`experimental-stokes-motion` feature, couples a finite native body to
`fs-flux::ale::translating::TranslatingAleFlow`. The same reference mesh moves
through world space as x = X(t) + y; there is no small-displacement bound.
Original BDM1/P0 convection, pressure, viscosity, source density and frame
acceleration remain with the nonlinear ALE owner. No fitted drag, prescribed
body trajectory, added-mass coefficient or second fluid model is introduced.

For each trial endpoint speed the original flow solves its own nonlinear
equations from the same accepted relative-velocity history. Its retained
discrete pressure/viscous reaction, including strong-normal constraints and
full weak tangential SIP flux, drives the original native rigid impulse owner.
Depth converts N/m into whole-body N. The actual stored momentum increment
closes m(v1-v0)/dt = F_fluid + F_external - k*q1 - c*v1. Outward advective
momentum is separate and is never added to the impermeable body's wall load.
Pressure datums retain their native component-local meanings. The native COM
pose supplies the moving load origin; transverse guide force and guide couple
are reported explicitly.

The caller supplies an ascending physical speed bracket and a maximum number
of native flow solves, including rejected trials. Safeguarded secants with
periodic bisection require an actual sign bracket; this does not assert that
every nonlinear problem has a unique root. Every trial retains the native
fluid's nonlinear acceptance. Represented body/mesh/drift rates and the maximum
unrelaxed normal-velocity trace correction use the same authored m/s tolerance;
both nonlinear equation residual gates remain active. Body momentum, pose,
fluid history,
world translation and clock publish only after momentum, mechanical work and
final cancellation admission. Unbracketed roots, budgets, native failures and
cancellation preserve both continuing physical owners.

Spring storage, signed external/interface work, physical mechanical dashpot
loss and numerical backward-Euler body/spring loss remain separate. The native
fluid report retains world kinetic energy, momentum and moving-surface flux.
The body work identity is not a complete fluid/body energy or passivity proof;
bulk viscous heat is evaluated from the accepted native strain without
inferring thermal-pressure feedback. The entire outer
mesh boundary also translates: endpoint callbacks supply WORLD conditions on
that moving control surface, not stationary obstacles left behind by the
body. Rotation, deformation, changing topology, variable density and general
mesh remapping remain outside this guided rung.

Four focused regression groups cover signed motion farther than one domain
width with independently integrated fluid/body mass and momentum, physical
forces sampled at the translated endpoint with spring/dashpot work, nonlinear
viscous recoil and continuation, and rejected-root/budget/datum/late-cancellation
atomicity. A native ALE regression additionally compares discard and commit
with the unchanged direct-step path. All four coupling groups and nine
ALE/absolute-correction cases passed in the selected-original-source native
build. Full Cargo/dependency integration remains separate; unrelated source
adapters and runtime wrappers were excluded from this bounded build.

## Implicit pressure and viscous loads driving guided rigid motion

The opt-in `experimental-stokes-motion` feature exposes `GuidedStokesMotion`.
It binds the original material `FluxSystem` to one `fs-mbd::RigidBodyState`
with a translation guide, physical extrusion depth, spring and dashpot. Native
mixed BDF1 Stokes solves compute zero-motion loading and the homogeneous unit
velocity response. Their scalar Schur complement closes the body equation;
a final native solve uses the actual representable body velocity. The response
is solved homogeneously rather than subtracting large absolute-pressure loads.
No fitted added mass, drag law or separate flow solver is introduced. Three
native solves are required per step; trial fields are released after extracting
their response. Material transport is bound once and never requeried by trials.

Selected-edge pressure/viscous forces use the retained discrete reaction
functional, explicit pressure datums, and depth conversion from N/m to N.
Mean-zero components require an authored pressure offset; traction components
retain their datum. External force must explicitly include any pressure on the
body's other side. Ideal-guide transverse force and reference-COM torque are
reported separately. Initial off-guide/rotating states are refused, strong
normal velocity is checked, and the inherited weak SIP tangent treatment remains.

Fluid history, original body momentum, translated pose and time publish together
after force, body-work, finite-arithmetic, displacement-domain and cancellation
gates. Spring storage uses the actual represented displacement increment.
Dashpot loss and backward-Euler mechanical loss are separate; the latter is not
heat. Fluid kinetic change and conjugate interface work are exposed, but the
body-work identity is not a global fluid/body energy or passivity proof. Native
fluid viscous and SIP terms have not acquired an independent whole-system audit.

All five fixed-reference, translating and rigid-ALE body/flow owners expose
`bulk_viscous_dissipation_j`. It uses the native material solution's physical
`2 mu D:D` integrated over cell areas, multiplied by the original extrusion
depth and timestep. Source-backed rheology supplies the accepted endpoint
viscosities; a rigid frame translation or spin does not add strain heat.
This is BDF1 endpoint quadrature, not an exact integral along the unobserved
within-step trajectory. Physical wall-slip loss, mechanical dashpots, SIP,
upwind and backward-Euler numerical losses remain separate. Ordinary motion
steps do not update thermal state; this diagnostic does not establish a coupled
energy theorem.
Nonfinite physical heat refuses before any body or fluid history is committed.

`GuidedStokesMotion::step_with_heating` deposits that same native endpoint
strain heat into original finite thermal cells over the identical interval.
The explicit triangle-to-cell map must match the original x/y triangle prisms
extruded along +z, including centrelines, depth, volume, surface and frozen-density
mass. Native thermal subcycling preserves caloric laws, authored sources,
contacts, reservoirs and scheduled events; generated heat never replaces a
boundary source. Both clocks must agree, and motion plus thermal history publish
only after the native flow/body and thermal gates and final cancellation pass.
This first-order one-way split retains frozen fluid transport and fixed geometry;
it adds no advection, thermoviscosity or global energy/passivity theorem. Caloric
identity remains caller-authored because a transport snapshot cannot supply it.

`GuidedStokesMotion::step_with_slip_heating` additionally deposits the accepted
native Navier-slip dissipation `depth * integral(beta * |(u-wall)_t|^2)` into
the same thermal transaction. Every exterior Navier-slip edge has one explicit
`StokesSlipHeatRoute`, including perfect-slip edges with zero dissipation. The
adjacent mapped fluid triangle receives the authored fluid fraction; an active
wall cell outside the fluid map receives the remainder. Multiple edges may heat
one wall cell. Partitions, wall caloric identity and physical wall correspondence
are caller-authored; this path infers no flash-temperature model or conductive
link. Duplicate/missing routes and other boundary kinds refuse before flow work.
Retained endpoint reactions supply heat without repeating material or boundary
callbacks. Separate deposited bulk and slip energy, and their combined source
integral, are checked against the original physical depth and motion interval.
Authored thermal sources and events remain additive. The original bulk-only
method retains its behavior. Analytic slip-Couette regressions cover fluid/wall
partition endpoints, shared wall receivers, frame and depth changes, zero
perfect-slip loss, source retention, invalid routing and atomic cancellation.

The existing motion regressions now independently check zero strain heat in
an accelerating uniform translation, analytic Couette heating under rigid
spin/Galilean boosts, shear-dependent endpoint viscosity, and the steady
moving-wall work split into bulk heat plus interfacial slip.

`GuidedStokesMotion` remains fixed-reference, small-displacement linear
incompressible Stokes, distinct from the finite-translation owner above. The
caller authors a geometrically appropriate displacement bound, and every
endpoint enforces it. The mesh does not move; no ALE, transported volume,
convection, finite-motion, thermal pressure-work or rheological evolution is
claimed. Fixed source/boundary callbacks must be pure across trials. Polling
between native solves and before publication does not make the synchronous
native linear solve cancellable internally. Four authored G1/G3/G4 regressions
cover exact closed-box added mass and full-box backward-Euler work, open-pressure
driving, transient viscous back-reaction, gauges/depth, guided-state admission,
late cancellation and rejected-step continuation. Native test execution remains
a separate verification obligation.

## Source-backed evolving viscosity across film layers

`LiquidFilmImpactSystem::new_with_layers` accepts explicit source-backed
`(film, channel, layer)` bindings on the original normal-pressure
mechanical/enthalpy clock. Every channel/layer of a selected island is covered
exactly once at its existing finite heat destination; uniform/channel
assignments cannot overlap. Original normalized interfaces, finite inventories
and nonthermal source context remain fixed.

Every endpoint resolves each original liquid law at its actual accepted local
temperature and checks that layer's own relative viscosity change. A small
effect on the hydraulic equivalent cannot hide excessive change in a thin
layer. One group compilation stages the native hydraulic coefficients and
the exact profile work fractions, retaining common Barus separation. Profiles,
sources, mechanical coefficients, reservoirs and any thermal-network state
publish together only after the existing accepted-step gates. The completed
interval uses its old profile's work; endpoint viscosity never reprices heat.

The same implementation serves adiabatic, coupled-transport and network
entries. No second clock, extra heat source, layer advection/mixing or inferred
inventory is introduced. Separate sliding-film routes cannot alias these
normal-pressure bindings. Six focused regressions cover independent old/new
profile work and actual motion, one-layer recovery, network conservation,
source/change-limit/cancellation rollback and mixed source segments.

## Sliding Reynolds films in native mechanical time integration

`ImpactSystem::with_sliding_film` attaches the existing `ShearedFilm`
pressure/entrainment/wall-traction operator to the original joint nonlinear
mechanical equation. Its prescribed wall rows use the same mass-normalized
basis and actual discrete-gradient rates as contacts and normal squeeze films.
The complete resisting force and analytic directional tangent enter both
reference and prepared execution; the fs-phs dissipation receipt counts their
work once. No force lag, prescribed driving waveform or second clock is added.

Pressure-only and sliding islands have separate stable insertion indices and
share the existing four-island work ceiling. Geometry and source profiles
remain immutable during each solve. Endpoint gap admission, compressible-gas
domain exclusion, mapped support-memory recomposition restrictions and atomic
motion/history/cancellation acceptance remain enforced. This binding retains
the native sliding law's prescribed Newtonian profile and zero-Barus scope.

Three focused G1/G3 regressions exercise the independent two-mass Couette
relaxation solution, reference/prepared force and analytic-tangent reduction,
momentum and physical work balance, rejected-step retry and domain/budget
admission. Runtime execution remains a separate verification obligation.

## Layer-resolved incompressible film heat

`render::plate::impact::material::heating::SqueezeFilmLayerHeatRoute`
and `HeatedImpactSystem::new_with_layer_losses` locate an original
passage's pressure-flow work within its retained `GapViscosityProfile`.
Every prescribed viscosity layer explicitly names one existing finite caloric
store. Several layers may share a store. The source layer interfaces,
viscosities and positive integrals of `(z-z0)^2/mu` remain queryable;
no wall heat split, layer inventory, density, mixing or advection is inferred.

The actual channel coefficient must exactly match the retained profile's
hydraulic-equivalent viscosity, both at attachment and before each deposit.
A common native Barus factor preserves these pressure-work fractions.
Whole-film, whole-channel and layered routes cannot overlap or duplicate
physical work. Existing uniform/channel liquid feedback refuses replacement
of a retained layered film. Dedicated source-backed layer bindings update
viscosity and physical work fractions together at accepted endpoints. Couette work is outside this pressure-only
route: combined pressure/shear layer cross terms require their own partition.

The original pressure solve is evaluated once per selected island at the
accepted discrete-gradient velocity and midpoint geometry. Its physical
channel work is counted once in the mechanical loss gate, irrespective of
the number of representational layers. Nonlargest layer energies are summed
with compensation; only the rounding remainder goes to the largest physical
contributor, preserving tiny positive deposits. All deposits remain scratch
until the existing adiabatic, transport or thermal-network owner accepts
the same mechanical step. Any chart, profile or cancellation refusal leaves
motion, stores and clock unchanged. Five focused regressions cover exact
uniform reduction, independent two-layer work, tiny contributions, invalid
coverage, late rollback, network conservation and feedback admission.

## Physical force envelopes in modal performance files

The v1 broadband and v2 coupled continuous-tone file paths accept optional
`envelopes COUNT`, unique `envelope TONE_INDEX POINT_COUNT`, and
`gain SAMPLE GAIN` records after the original physical tone rows. Original
tone-row indices survive per-voice grouping and exact-frequency coalescing.
Each nonempty curve begins at sample zero, increases strictly within the
declared horizon, interpolates finite signed gains linearly, and holds its
last gain. Undeclared curves retain unit gain; zero envelopes keep the old
constant preparation and numerical path. Held DC event assignments remain
independent physical inputs.

The shared parser bounds aggregate knots and reserves payload before allocation.
Per-driver endpoint buffers and, for the independent path, additional scaling
work are charged against the existing source allowances. Envelope-capable
preparation admits complete affine source coordinates and original physical
limits. Carrier coefficients remain unscaled accepted state, including through
silence and re-attack. Both gain endpoints multiply that same interval-start
carrier basis before the exact physical step. No output pressure is faded and
no resonant, material, receiver or delay history is reset. Source phase, held
events and gain cursors publish only after the physical/observation step accepts.

Both existing CLI children expose the number of physical force envelopes and
retain the original fit/delay, decimation and PCM semantics. The complete
examples are `modal_enveloped_tones.performance` and
`coupled_enveloped_tones.performance`. Focused file checks cover independently
evaluated whole ramp segments and physical work, original tone indexing,
silence/release/re-attack phase, bitwise callback partitions, residual physical
ringdown, malformed knots, resource admission and rejected-step cursor
preservation. Native executions of both actual CLI children produce 4800-sample
48 kHz PCM16 WAVs with zero clipping through the unchanged encoding path.

These are prescribed physical force histories. Knot refinement is independent
of the mechanics clock; sharp envelope changes are not certified bandlimited
by carrier admission. This does not infer a nonlinear excitation law or alter
the estimate-only authority of sampled receiver transfers.

## Continuous physical tones in coupled performance files

Continuous coupled force files (`render::schedule::force::file::coupled_tones`):
  a trailing `coupled_tone_limits`/`tones` section on the original v2 mechanical
  schema projects physical cosine/sine newtons through its unchanged signed
  actuator columns. The existing parser, held-force compiler, complete spring
  preload and `CoupledModalTones` owner supply one coupled trajectory, exact
  continuous force work and independent original heat channels. Held events
  remain interval DC; source phase starts at performance time zero and advances
  only with an accepted network sample. Rejected steps retain accepted phase,
  event cursor, state and work/heat. File-wide source payload and all phase arrays
  are reserved before the remaining dense byte budget is admitted; projection,
  clock, guard-frequency and original physical limits remain binding. Without an explicit sampled receiver, the CLI observes real instantaneous
  pressure/velocity maps; imaginary natural-mode quadrature and contact/friction
  schemas refuse. The optional complete-network causal receiver is specified below. This is an authored reduced assembly path,
  not a radiation or instrument-validation claim. `coupled_tones/tests.rs`
  compares actual file callbacks with an independent coupled normal-coordinate
  solution, exact preload, physical work/heat and partition-invariant phase,
  and checks transactional and admission refusals. See
  `examples/COUPLED_CONTINUOUS_FORCES.md` for units, limits and a complete file.

## Coupled causal receivers in performance files

The v2 continuous-tone adapter accepts an optional `coupled_broadband` suffix
followed by one existing `transfer`/frequency/withheld block spanning the complete
voice-major modal basis. The original mechanical parser, signed port compiler,
static preload, events and continuous phases are retained. Supplied sampled
transfers replace legacy observation maps before tone preparation; they do not
change mechanical coefficients or add reaction loads. The shared receiver owns
the network pressure ceiling. Phase and held-event publication wait for the
source/filter/delay transaction to succeed.

Transfer samples use pressure per generalized modal velocity in the original
signed normalization and exp(-i omega t), with the declared propagation delay.
The file exposes a separate bounded receiver preparation allowance and reuses
the existing complete transfer-row parser, fit controls and actual sampled-delay
held-out gate. `music_render modal` writes fit/delay/producer evidence as
EstimateOnly before the usual explicit observation decimation and PCM scale.
The acoustic histories start at zero; mechanical preload supplies no acoustic
prehistory. Fractional delay and held-out samples retain their approximation
scope, and these authored transfers do not establish instrument validity.

`examples/coupled_broadband_tones.performance` supplies two signed first-order
pressure transfers, including negative direct gain, with an explicit three-sample
delay. File-to-callback regressions bind every modal transfer row, preserve
bitwise phase/output across callback partitions, compare every accepted physical
state and work/heat frame with the neutral-map source, and refuse incomplete
bases, bad clocks, held-out misfit and insufficient preparation budgets. See
`examples/COUPLED_CONTINUOUS_FORCES.md` for units and complete record semantics.

## Continuous physical tones in modal performance files

The v1 broadband performance-file path now accepts optional `tone_limits`,
`tones`, and `tone VOICE PORT OMEGA_RAD_S COSINE_N SINE_N` source records before
the receiver block. These prescribe continuous physical newtons through the
original signed work-conjugate port columns; they do not author acoustic output
partials. Existing held assignments remain the DC channel, and static preload
retains its pre-window held-load meaning. Frequencies and phase remain on the
original mechanical clock, referenced to performance time zero across callbacks.
The existing continuous multi-tone mechanical/receiver owner retains all cross
terms in work, viscous heat and pressure. Candidate phase advances only after
an accepted complete observation step. Tone/driver storage and setup/step work
are explicitly bounded across affected voices, with conservative equal shares;
projection is charged against the existing force compiler's work allowance.
The original clock/Nyquist guard and physical model/receiver limits still apply.
Held-only files retain their original stepping path. Complex constant legacy
quadrature refuses continuous forcing. Native source regressions compare with
an independent continuous-force response, expose held-sampling bias, and cover
signed-port cancellation, rejected-step phase/work/heat preservation and invalid
source/clock/resource refusal. The example and record semantics are documented
in `examples/MODAL_CONTINUOUS_FORCES.md`; they make no instrument-realism claim.

## Nonlinear resonance branches through the original pHS

`phs_orbit::continue_forced_frequency_response` connects the complete original
port-Hamiltonian ODE to fs-orbit's existing forced pseudo-arclength continuation.
The continuation parameter and forcing fundamental are the same angular
frequency. Original DC inputs, complex multiport Fourier amplitudes and relative
phases remain fixed; no modal coupling, resistance or material-memory coordinate
is removed. This traverses frequency folds rather than restarting isolated
frequency solves. It preserves native branch order, which can reverse frequency,
and the terminal accepted point may pass the requested end frequency.

Each `PhsFrequencyResponsePoint` carries the unmodified native continuation
point, including original-state Fourier coefficients, residual, AFT sampling
defect and step size, plus an independent `PhsOrbitAudit` at that point's accepted
frequency. The audit uses the original Hamiltonian, port work and full resistance
matrix. Invalid endpoints, work settings, physical port shape, unresolved
authored forcing, source refusal or continuation exhaustion return an error
without a partial branch. Authored bandwidth is admitted before source sampling.

The retained harmonic order is fixed along one branch. AFT refinement checks
retained-row sampling, not missing response harmonics; callers must repeat
retained orders and physical-audit quadratures to measure convergence. Cycle
power balance does not establish waveform completeness or Floquet stability.
The branch can include unstable periodic orbits; no physical jump threshold or
hysteresis width is certified by continuation alone. This synchronous path uses
the native step/Newton/AFT budgets and 32..=65536 audit samples, and has no
cooperative cancellation hook or time integrator.

Focused regressions cover analytic complex multiport response across linear
resonance, DC and higher forcing harmonics, original work, Maxwell memory and
frequency-dependent loss, both folds of a cubic-storage response with generated
higher harmonics, admission before source evaluation, and source refusal.

## Numerical Floquet stability of forced pHS resonance branches

`phs_orbit::analyze_forced_orbit_stability` exposes the complete original ODE to
fs-orbit's existing shooting and monodromy owner. The forcing frequency and phase
must match the supplied HB or continuation point; `PeriodicPortInput::at_frequency`
preserves authored Fourier amplitudes and phases at an accepted branch frequency.
All damping and material-memory coordinates remain. Each native solve fixes the
forcing period, and corrections exceeding caller-authored absolute tolerances in
any original state coordinate return `SeedDisagreement`.

The report retains native shooting states, complete multipliers and residuals on
N, 2N and 4N time grids. Stable/Unstable numerical status requires both adjacent
state/radius comparisons to pass and all radii to lie beyond the explicit
unit-circle guard; otherwise it is Unresolved. No forced multiplier is removed
as an autonomous phase. Dimensions, finite data, harmonic bandwidth and budgets
are admitted before source callbacks. Source and native solver refusals remain
errors, including the native suspected-torus refusal.

This is local numerical evidence for smooth ODEs. Native forward-difference
monodromy and three finite grids do not provide an eigenvalue enclosure,
existence theorem, global stability claim, hybrid-event map or cooperative Cx
cancellation. Regressions cover analytic complex-forced oscillators, all three
Maxwell-memory multipliers, stable and unstable Duffing branches with an
independent Liouville determinant, and unresolved/refusal boundaries.

`analyze_forced_orbit_stability_variational` and its `_analytic` entry use the
same candidate admission, physical seed-agreement checks and three-grid
classifier, with the existing original-state variational shooting solver at
each level. The automatic entry requires complete declared quadratic storage;
the analytic entry takes the complete continuous Hessian action of nonlinear
storage. Both retain material memory, physical source phases, the entire
scaled monodromy, all multipliers, normalized residuals and evaluation counts.
They never silently replace an unavailable or refused tangent with finite
differences. Cancellation is checked before admission, between levels and at
the shared solver's stage/product boundaries, with no partial assessment.

Physical scales and numerical work/payload/map caps apply separately to each
native solve. Exactly three solves run, and previously returned level data are
retained outside the next solve's payload cap. This is a measured temporal
ladder of exact discrete tangents, not an enclosure of continuous Floquet
multipliers, an omitted-harmonic bound or a hybrid-event stability analysis.
Focused regressions compare constant-system multipliers with an independent
RK4 polynomial across excitation scales, retain all Maxwell poles, and check
an exactly forced cubic cycle's parametric instability and Liouville determinant.

## Continuous physical impedance on accepted loaded assemblies

`assembly_harmonic::PreparedAssemblyMimoHarmonicResponse` borrows the original
admitted `PreparedMimoLoad`, enforces its continuous source band, and evaluates
the original proper rational matrix at `s=-i*omega`. It attaches the full self
and mutual load through `AssemblyHarmonicPorts`, retaining the accepted body's
original M/K, glue C, physical force/flow maps and static loaded mean shape.
Reused frequency factors return full component perturbations and external
reaction forces. Independently reconstructed physical input, glue loss and
external load work must balance.

`solve_with_frequency_derivative` consumes the rational owner's analytic
`dZ/ds` with the exact `-i` chain rule at `s=-i*omega`. Every original real and
conjugate pole, mutual entry and physical source slope enters the same-factor
loaded tangent. The source derivative matrix and evaluation work are admitted
before allocation and deducted from the shared frequency controls. Original
physical motion, port reactions and work slopes support actual receiver
resonance refinement. Source band, history and both admission reports retain
their original meaning; frequency slopes add no broadband certificate.

The original source clock, accepted pole history, port ordering and both
passivity report classes remain unchanged. No fitting, sampled-transfer
substitution, acoustic surface, external storage/heat or broadband certificate
is inferred. The focused MIMO test covers the analytic two-port Maxwell source,
Fourier axis, original band, resource/cancellation refusals and unchanged
history/reports. Loaded-shell port regressions separately cover resonance,
nonmatching ties, physical reaction work and mutual-load falsifiers.

## Exact joint assembly and rational-load resonance transients

`assembly_harmonic::transient::PreparedAssemblyMimoTransient` evolves the original
accepted body's complete retained M/K, glue C, physical attachment map and every
original `PreparedMimoLoad` memory coordinate in one continuous linear generator.
Applied incremental component force is held over the original source clock;
attachment flow evolves with the body and load memory. The connection uses no
inverse of the load's direct term, no held-flow split and no discarded pole or
mutual entry. Zero-direct-term Maxwell and complete nonreciprocal MIMO loads are
supported through their existing admitted source owner.

The structural source comes from `AssemblyResponseSource`, shared with the
original standalone response. The existing fs-phs prepared quadratic-flow owner
caches the joint exponential and independent trajectory integrals once. States
remain [q,p,original load coordinates], with p=M*v and the source's original
complex-pair real chart. Initialization preserves load history; checkpoint
restoration publishes complete structural and load histories atomically.

Each accepted step independently reconstructs physical applied work, attachment
flows, endpoint reactions and reaction impulses. Structural energy change, the
integral of original glue dissipation and signed external load work must balance
physical applied work. The load-work integral may be negative when existing
memory releases energy; no load storage or heat certificate is invented. Static
loaded mean shape and all original component/support/tie maps remain unchanged.
Failures of source, physical maps, numerical accounts, budgets or cancellation
publish neither structural state nor load history.

Connection dimensions, allocation envelope, construction/contraction work and
the original exact-owner cache budget are explicit. Source clock, declared band,
port ordering and continuous/sampled passivity report classes remain unchanged.
Exact evolution of the supplied rational equations does not certify broadband
source validity, missing structural modes, nonlinear motion or a changed contact
branch. Focused numerical checks cover analytical loaded decay/frequency through
large time steps, history release, full coupled mass/damping/mutual loading,
independent q/v integration, exact subdivision and original complex-pair charts.

## Original exterior BEM geometry to complete modal mouth loading

`modal_mouth::BemModalMouth` binds a caller-selected circular planar cap on the
original closed outward-oriented exterior triangle mesh to fs-duct's signed
axisymmetric mouth basis. A unit modal volume flow imposes the physical normal
velocity trace psi/S. Native aperture work projection returns the complete modal
pressure/volume-flow impedance directly, including mutual terms and higher modes
whose net volume flow is zero. The adapter never applies the per-port net-flow
conversion of `ApertureImpedance::acoustic` to these modal traces.

Original geometry validation, cap normals/radius, panel uniqueness, area and the
unmodified discrete basis Gram matrix are checked before solving. A caller-owned
Gram-defect limit refuses an under-resolved cap; no trace renormalization or mesh
repair is performed. Finite Gram moments do not certify complete disk coverage
or absence of self-intersection; the selected cap remains the caller's physical
surface. Geometry is not remeshed and no missing cap is generated.

`bake` uses the explicitly selected native collocation or Galerkin aperture route,
shares each frequency factor across modal loads, and preserves the original raw
BEM artifact and diagnostics. Every complete complex matrix must pass fs-duct's
pointwise Hermitian passivity admission, with its original frequency, radius,
exterior gas and mode count. The batch admits trace, matrix and projection work
ceilings and cancellation before publishing all requested loads together.

The harmonic bake is finite-panel, axisymmetric loading of the supplied rigid
exterior geometry; its separate causal realization is described below. Regressions cover
signed modal traces and mesh refinement, independent surface/port work equality,
native Galerkin radiation connected to a loaded duct at two frequencies, mutual
loading liveness, and geometry/resolution/budget/cancellation refusals.

## Continuous multimodal duct with complete rational radiation memory

`mm_line::continuous::rational::RationalModalLine` joins the original lossless
Galerkin duct's full radial-mode state to every original `PreparedMimoLoad`
coordinate. `ModalRationalLoad::from_bem` retains the actual modal-mouth producer's
geometry, signed basis, independent source samples and fit reports. The generic
constructor requires an explicit producer-declared pressure/volume-flow basis.
Mouth radius and coherent static pressure must match the duct attachment.

The sole interval is the original load clock. All incident pressure-wave modes
are held at the inlet, while duct pressure, outward mouth flow and the complete
mutual radiation state evolve jointly throughout the interval. The existing
fs-phs exact exponential/Gramian owner prepares one full finite generator and
independent linear/quadratic observation maps. Warm steps use bounded matrix
contractions and private candidate allocations; no allocation-free claim is made.
Held incident coordinates condition the numerical preparation only and never
enter the public physical history or receive assigned stored energy.

Invertible direct D uses the existing terminal solve for D*U=P_duct-C*z. Exactly
zero D uses the existing differentiated pressure-constraint owner, including
every held input in the original generator. Neither path inserts a resistance,
diagonalizes mutual terms, advances a separately held-flow load or discards pole
history. Mixed-rank D, singular higher-index constraints and inconsistent
accepted zero-D history refuse without projection. Complete original-state,
workspace, work and cancellation admission precedes preparation/publication.

Independent exact integrals report inlet pressure-flow work, incident and
reflected wave energies, duct-mouth work, and original (C*z+D*U)*U load work.
Original duct storage change closes against these exchanges. Endpoint pressure,
constraint and work agreement must pass before any physical state, original load
history or output changes. The mean reflection relation retains every mutual
entry and changes no history. Interval-mean waves remain distinct from integrated
wave energy. Rational work is signed exchange, not inferred heat or storage.

The original load's band, clock and continuous/sampled classes remain unchanged;
component admission does not qualify the whole closed loop. Exact propagation
concerns this admitted finite spatial operator. Spatial/radial convergence, inlet
wave sampling, physical source fit and omitted duct-loss physics remain separate.

## Self-excited resonance through original physical ports

`phs_orbit::AutonomousPortFeedback` is an explicitly authored instantaneous law
producing physical port efforts from the complete original state and conjugate
flows G^T*grad(H). Native autonomous HB and adaptive HB solve the resulting full
ODE without changing its Hamiltonian, native resistance or material-memory
states. Stateful excitation, transport and controller memory must remain in
explicit physical states; deterministic callbacks may not hide evolving history.
The source can supply energy. Component passivity does not establish closed-loop
passivity, and the adapter does not validate the authored excitation physics.

Autonomous stability selects a resolved transverse coordinate at the supplied
HB waveform's phase zero and temporarily permutes it into native shooting's
coordinate-zero section. N/2N/4N native shooting solves for the period, restores
original state order and reports the actual section coordinate and value. Each
shot must retain section orientation and transversality and satisfy explicit
physical-state and period agreement with the seed. An accurate period alone
cannot admit an under-resolved waveform.

Every multiplier remains in the returned spectrum. Exactly one isolated
multiplier near +1 may identify the phase direction for the transverse numerical
classification; multiple neutral directions, including conservative 1:1 families,
remain Unresolved or retain the native solver's refusal. State, period and
transverse-radius agreement across all three grids and the unit-circle guard
control classification. Existing forced-response continuation is unchanged.

The independent cycle audit evaluates original H, physical source work and full
R loss on the supplied HB waveform; it is not a shooting-trajectory quadrature.
These smooth-ODE, finite-difference numerical results do not certify global
stability, a fundamental period, hybrid events or material/source validity. The
bounded native solves have no cooperative Cx cancellation. Focused checks cover
stable/unstable cycles, phase changes, Maxwell storage and heat, a resolved
Van der Pol cycle, conservative degeneracy and transactional source refusal.

## Loaded resonance poles from the original joint transient

`assembly_harmonic::transient::spectrum` analyzes the exact retained generator
already used by `PreparedAssemblyMimoTransient`. The accepted body's complete
M/C/K, physical port coupling and every original continuous MIMO memory state
survive, including mutual/nonreciprocal entries and port-dark states. No second
equation assembly, equivalent-viscous fit, pole pruning, source-history reset
or sample-clock-derived pole estimate is introduced.

Declared time scaling and the existing numerical coordinate scales prepare the
complete nonsymmetric spectrum for fs-spectral. Every eigenstate returns to the
original [q,p,load] chart, with p=M*v and the original pair-memory convention.
Independent residuals check kinematic, momentum and load equations. Explicit
payload/conversion and spectral budgets, complete-basis conditioning admission
and cancellation govern publication.

Poles use x(t)=v*exp(s*t) and omega=i*s under exp(-i*omega*t). Signed decay is
-Re(s), oscillation frequency is |Im(s)|/(2*pi), and optional cycle-decay Q is
|Im(s)|/(-2*Re(s)) for finite decaying oscillation. A pole need not be a peak in
an observed transfer. No damping floor or forced-response linewidth is inferred;
mixed-unit eigenstate normalization does not define stored energy.

All roots remain regardless of source-band membership. The original band and
continuous/sampled classes remain unchanged; a frequency-membership flag does
not certify off-axis source continuation. Results describe the complete retained
finite generator, not omitted modes or continuum completeness. Focused checks
cover analytical Maxwell poles, clock/unit invariance, full mutual/nonreciprocal
eigenstate evolution, dark memory, and an analytical memory-only high-Q case.

## Complete geometric modal radiation memory

`modal_mouth::BemModalMouth::realize` connects every original signed modal mouth
impedance entry to the existing physical-source MIMO fitting owner. The shared
owner supplies the three-order ladder, independent training/held-out/causality
BEM sweeps, dispersion consistency, explicitly reported bounded residue repair,
continuous passivity and actual interval-mean ZOH passivity. The original Pa s/m3
matrices are conjugated once into s=+i*omega; higher-mode net-flow gains remain
diagnostics and never rescaling factors. The raw reciprocity defect is bounded
before reciprocal identification and remains in the reported fit error.

The result retains the original `PreparedMimoLoad`, physical clock, exact BEM
route, circular patch, coherent exterior gas, radial mode count, authored source
and approximation bands, source budgets, all three unchanged physical bakes and
complete fit reports. Every physical query, including the causality grid's lower
probes, must fit the source support. The approximation band must reach mechanics
Nyquist, and the complete three-sweep source budget is checked before sampling.
Unsupported coverage, native refusal or cancellation publishes no realization.

These coincident modal coordinates share one physical mouth reference plane;
no inter-aperture propagation delay is inferred. Geometric diffraction memory
and mutual terms remain in the shared rational state. The sampled report stays
GridOnly. Mesh/basis accuracy, unsampled frequencies and closed-loop resonance
or decay qualification remain distinct from realization error and passivity.
This advances the modal-mouth path without closing the broader multi-aperture
delay and instrument-audit requirements of M2.5-W3b. Focused G1/G0 checks retain
native two-mode BEM samples, an independently loaded duct, higher-mode mutual
pressure and persistent pole memory, plus early band/clock/budget refusals.

## Matched modal radiation loading and receiver fields

`modal_mouth::observation` prepares the original Galerkin exterior geometry once
and reuses its complete signed modal source batch. At a frequency, one boundary
factor supplies both the full `ModalMouthImpedance` and the original P1-pressure/
P0-auxiliary-density solutions. The native adaptive receiver owner integrates
those same densities to physical pressure and particle velocity. Panel pressure
means remain the work projection only; they do not replace the P1 receiver field.

Receiver positions and actual triangle clearances are admitted once. Each
frequency prepares one common adaptive row bank, with separate componentwise
gradient discrepancy checks. Coherent modal-flow cases reuse the resulting
field bank without new source solves or quadrature. Signed higher modes with
zero net volume flow remain present, and original complex phases determine both
interference and the returned mean intensity. The output also reports P=Z*U
and complex mouth work (1/2)*sum(P*conj(U)) at the same frequency and reference
plane. It can feed the native ideal pressure/velocity microphone projection.

Native source geometry, exclusion-band, wavelength, angular-resolution and
quadrature limits remain active. Separate source/receiver payload and contraction
caps cover the complete requested bank. Cancellation or any finite/admission
failure returns no partial response; immutable source/field preparations remain
reusable. Propagated quadrature discrepancies do not cover BEM mesh, solve or
physical-model error. This is an axisymmetric harmonic receiver path, with no
causal observation filter or microphone hardware model implied.

Three G1/G4 regressions compare coherent pressure, particle velocity and mouth
work against a separately solved combined physical source; close a lossless
loaded duct's incident/reflected power against the same radiation work; retain
directional microphone response and higher-mode interference; and check band,
geometry, budget, finite-input, cancellation and deterministic reuse boundaries.

## Finite load-resistor feedback on three-dimensional winding histories

`SpatialWindingThermalNetwork::step_material_candidate` binds every ordinary
native graph resistor to its existing `CellResistor` material body and explicit
electrical property basis. The same source-resolved `MaterialNetwork` supplies
endpoint `R(T)`, original conductor length/area and caloric state. A resistor
body cannot also claim a whole magnetic field volume; its heat route and source
geometry must describe the same separate finite receiver.

Before any trial, the accepted native resistance must equal its original source
evaluated at the accepted body temperature. Previously retained graph bindings
must also match the original material, query, source choice, basis and geometry.
Every endpoint refresh checks complete branch correspondence and original-source
temperature admission, then reuses `with_graph_resistors` without changing A,
winding current, capacitor charge or native flux history. Material providers
cannot silently substitute a new resistor while solving its temperature.

The common nonlinear spatial driver receives the already resolved thermal
coefficients, so graph loads and field materials share one endpoint source query.
Strand, bulk eddy and graph heat remain separate native receipts deposited once;
the original full system energy and magnetic entropy equations are unchanged.
`heat_assembly_candidate` supports these same complete material loads, cooling,
contacts and radiation through the existing heat assembly. Unrepresented electrical
contact work is still refused. The numeric-resistor API refuses omitted material
bindings on a circuit that already retains source-bound graph loads.

One additional regression in `tests/spatial_induction_winding_enthalpy.rs` uses
separate finite winding/core and load-resistor bodies over two accepted intervals,
including the full heat assembly. It checks independently reconstructed
`rho(T)*length/area`, each body's physical heat, changed current, retained source
metadata, and accepted-resistance/omitted-binding refusals. Pinned Rustfmt and
focused source/accounting review passed; native compilation/tests and the required
DSR gate remain unrun in this environment.

## Finite heat on nonlinear three-dimensional winding circuits

`circuit_enthalpy::spatial_induction::nonlinear::winding` couples the retained
`SpatialNonlinearWindingCircuit` to the existing finite-cell enthalpy residual.
`SpatialWindingThermalNetwork::step_candidate` and `heat_assembly_candidate`
return one `SpatialWindingEnthalpyCandidate` containing accepted circuit, A–phi
and caloric histories. Every trial begins at the same original accepted state;
source/domain refusal, failed local or global energy checks, exhausted work or
cancellation publishes none of those histories.

Original material cards are queried at each proposed caloric endpoint. Stranded
resistivity and bulk eddy conductivity therefore affect the actual native terminal
currents. Original source selection, electrical basis, geometry and conducting
support are retained. Every strand, bulk conductor and authored graph resistor
has one explicit finite heat destination. Intrinsic strand Joule heat, exact
tetrahedral P1 bulk `E dot J` heat and native graph heat are deposited once;
numerical magnetic and omitted-end dissipation remain a separate account.

The existing nonlinear spatial material/caloric owner now has a private generic
driver seam shared by prescribed-current and native-terminal histories. It keeps
the same source admission, caloric preparation, heat transport, Newton iteration,
budgets and acceptance equations. The winding candidate supplies its actual nested
field-trial count, clipped against the outer remaining work allowance. The complete
`HeatAssembly` entry preserves original boundary/contact/cooling/radiation payloads
and the existing closed-radiation balance check. Unrepresented electrical contacts
are refused; complete graph material bindings use the source-driven load extension.

Without explicit magnetic caloric laws, the original magnetic constitutive law
must remain temperature independent. Optional existing linear `mu(T)` or
source-separable `H(B,T)` bindings close magnetic parameter work and reversible
entropy exchange. The whole-system audit substitutes original magnetic internal
energy for magnetic free energy while retaining every native circuit and explicit
end-inductor store. Neither parameter work nor numerical damping is relabeled as
material heat. Fixed recoil temperature response still needs its own entropy law.

The finite-heat adapter is stationary and supplies zero velocity; thermal advection
is not silently omitted. Original fixed geometry, fully solid constant-density
caloric charts and first-order coupling remain authoritative. A thermally owned
stranded tetrahedron must be fully occupied by its one original material source;
multiple insulated bundles are allowed when their fractions sum to one within the
explicit dimensionless roundoff allowance, at most `1e-10`, without rescaling.
Partial-fill or mixed-material coil composites require their actual mixture caloric
volumes and refuse here. The field/circuit owner can still represent their explicit
electrical filling without making a finite-heat claim. `step_candidate` retains
explicitly authored numeric graph resistances with complete heat routing;
material-bound loads use `step_material_candidate` or the complete heat assembly.

Three G1/G3/G4 cases in `tests/spatial_induction_winding_enthalpy.rs` cover two
accepted intervals with original conductivity feedback, independently integrated
geometric winding resistance and each physical heat channel; complete heat assembly
with cooling and independent `mu(T)` field internal energy; and source/partial-fill/
actual-budget/late-cancellation refusal with unchanged histories and deterministic
retry. Pinned Rustfmt parsing, whitespace checks and source/physical-accounting
review passed. These newly authored native tests, compilation and the required
`dsr quality --tool frankensim` gate have not run in this environment.

## Tetrahedral material current closes finite enthalpy

`circuit_enthalpy::conduction::spatial::SpatialConductionThermalNetwork` closes
the original 3D Whitney current field inside the existing finite-cell
backward-Euler enthalpy/Newton owner. Each trial refreshes every original
electrical source role and heat-transport binding at its actual endpoint
temperature and liquid mass fraction. Scalar and fixed-orientation orthotropic
sources use the same field, thermal residual and source work account.

`SpatialConductionDrive` accepts the original voltage/current-density traces
or signed total current through two equipotential sets of exterior triangles.
The original `&SpatialConductionBoundary` calls remain accepted. Every nonlinear
trial rebuilds the native terminal response after refreshing all bulk and
electrical-contact sources at its actual temperature and liquid mass fraction.
The terminal voltage and heterogeneous local current therefore respond to
`R(T,f)` inside the same enthalpy residual; the accepted-state resistance is
never frozen. The same drive reaches the complete `HeatAssembly` entry point,
including sourced thermal contacts, layers, radiation and finite coolant.

The reference voltage is a finite positive numerical probe, not an additional
source of heat. Only the recovered actual-current bulk and contact powers enter
the cells; terminal reactions independently check the signed current and `I V`
work. Negative current reverses voltage/current while preserving Joule heat.
Zero current deposits no electrical heat but still requires an admitted finite
terminal resistance. The native terminal owner refuses unresolved, touching or
disconnected electrode sets. Current is held fixed over the authored interval;
drive changes require the existing interval/event splitting. Each native
reference solve counts against the original circuit-trial work budget, and all
source/cancellation/field refusals preserve the jointly accepted histories.

`SpatialConductionDrive::TerminalCurrents` adds separately authored balanced
signed currents on multiple equipotential electrodes. A zero-current terminal
is a floating electrode, not an insulating boundary: its potential is solved
while its local faces may exchange current. Every residual trial first refreshes
the original bulk/contact material state, then extracts the native conductance
response and recovers the combined local field before depositing bulk and
contact Joule heat once. All N-1 basis field solves are charged against the
remaining shared circuit-trial budget before extraction; changing current signs
does not create another heat source. The existing complete heat-assembly owner
consumes the same drive without a parallel thermal solver. A focused
three-terminal regression uses independent backward-Euler elimination to check
temperature-dependent voltage, cell heating, floating potential, polarity,
energy and the actual two-solves-per-trial budget.

`SpatialConductionDrive::MixedTerminals` prescribes total current on an arbitrary
proper subset and voltage on the remaining electrodes. Zero-current electrodes
float as bulk and contact resistance change with the actual thermal endpoint.
Each trial rebuilds the native mixed terminal block under its own storage/work
bounds and charges all N-1 basis solves to the original shared budget. Supply is
independently summed from authored voltage differences times native reactions
and authored currents times solved voltage differences; common bias is removed
before arithmetic. The original field, contact partitions and finite-cell heat
routes remain the only Joule deposition. A manufactured three-electrode fixture
in `tests/spatial_conduction_enthalpy/contacts/mixed.rs` checks simultaneous
bulk/contact temperature feedback, floating and injecting middle electrodes,
local heat, polarity, large common bias, a current-controlled reference,
source work, solve budgets and late cancellation with unchanged accepted state.

`tests/spatial_conduction_enthalpy/current.rs` constructs an exact two-body
backward-Euler endpoint for thermally linked, electrically parallel 3D conductors
with temperature-dependent resistance. It checks unequal current redistribution,
local cell work, terminal voltage, reversal, reference-amplitude independence,
zero drive and malformed-electrode refusal. The original latent scalar/tensor,
two-sided contact-temperature/phase and full contact/radiation regressions also
exercise total-current drive. These are authored synthetic numerical oracles,
not measured material/device validation.

Each original tetrahedron maps to exactly one thermal receiver; several tets
may cover one receiver only with matching full volume and original mass/density.
The existing planar/volume-independent receiving-inventory and heat-summing
owner is reused without constructing a surrogate planar mesh. Confined phase
charts require constant equal density, complete source support and an explicit
liquid-MASS-fraction coordinate. Latent heat stays inside the original enthalpy.

Native per-tet E.J powers enter their receiving cells once. Terminal work is
computed independently from authored current traces and solved voltage-source
reactions. Existing authored thermal sources, boundary transfer and material
conduction links remain in the same residual. The original circuit-trial budget
counts native field solves, including finite-difference and fresh final
acceptance. Cancellation, source refusal, budget exhaustion and unsupported
event windows publish neither thermal state nor endpoint electrical operator.

`PreparedHeatNetwork::step_spatial_conduction_candidate` now connects this SAME
3D field/enthalpy residual to source-defined thermal contacts, integrated
conductive layers and closed-enclosure radiation. It uses the original source
cards and receiving volumes; contact pressure, conductivity and emissivity data
changes therefore reach actual temperatures, current density and terminal work.
Every nonlinear trial refreshes the complete heat assembly. Its original
interface/emissivity samples are retained with the accepted tetrahedral field.
Graph electrical work cannot be silently omitted, and field-dependent thermal
boundaries refuse without an interface-field owner. Radiation's signed numerical
net is checked as internal exchange rather than external energy supply.

Matching tetrahedral electrical contacts use the original ordered material
interfaces on coincident exterior triangles with distinct vertices and opposite
adjacent volumes. The native P1 operator integrates the voltage jump with
`A/(12 R'') * [[2,1,1],[1,2,1],[1,1,2]]`; it does not replace the interface by a
mean jump or lumped resistor. Finite contacts join electrical gauge components,
and their faces cannot also carry external electrical traces. Bulk `E.J` and
surface work remain separate and close against terminal and assembled work.

`SpatialConductionThermalNetwork::step_contact_candidate` and
`PreparedHeatNetwork::step_spatial_conduction_contact_candidate` refresh the
original bulk and interface bundles at the same trial enthalpy endpoint. Each
interface keeps its two actual temperatures, explicitly owned liquid mass
fractions, fixed external pressure/history and selected source claims. Surface
work is deposited once into the original adjacent volume receivers using an
explicit physical heat fraction. Empty-binding entry points refuse unowned
contact heat. The planar and tetrahedral paths share source admission, caloric
coordinates and heat partitioning. Rejected/cancelled candidates preserve both
accepted histories; successful field and thermal endpoints publish jointly.

Synthetic G1/G3/G4 cases compare independent two-cell electrothermal elimination,
material swaps, voltage/current reversal, tetrahedral subdivision, continuation,
and source T/f surface feedback through an isothermal latent interval in both
conductivity/resistivity bases. Scalar and oriented material snapshots produce
the same admitted aligned limit. Final cancellation, work-budget, foreign-source
and receiving-volume refusals retain retryable original inputs.
`tests/spatial_conduction_enthalpy/heat.rs` adds independent contact/radiation
balance, source-only contact-pressure and emissivity substitutions, continuing
field/heat history, source-domain refusal and cancellation. These are authored
synthetic numerical cases, not measured-device qualification.
`tests/spatial_conduction_enthalpy/contacts.rs` adds independent two-body
current/heat balances, pressure/material substitutions, latent-fraction feedback,
continuation and transactional refusal. Native field tests independently
integrate nonuniform triangle voltage jumps and their work.

### Moving-conductor heat on the same material endpoint

`SpatialConductionThermalNetwork::step_moving_candidate` uses the original
native motional-current solve inside that SAME finite-enthalpy residual. Every
trial refreshes the retained scalar or fixed-orientation conductivity sources
at the original receiving cells' actual temperature/phase before resolving
`J = Sigma(T,f) (E + u cross B)`. Only the original tetrahedron integrals of
`(E + u cross B).J` enter those cells. The existing full-volume/mass/source
admission, local heat routing, thermal links, event rules and shared solve
budget remain authoritative. The accepted field arrays move into an owned
snapshot without cloning and publish with their exact refreshed operator and
thermal history.

`step_moving_contact_candidate` also admits explicitly bound resistive joints
whose adjacent original conductor cells are stationary. Moving bulk may remain
elsewhere in the same native domain. Each trial resolves the retained electrical
contact law at both original source-side material states and partitions its
surface loss once into the declared thermal receivers. Bulk cell heat excludes
that surface loss; the total Joule ledger includes both. A moving or sliding
resistive interface still refuses because this path owns no surface motion work.

`heat::step_spatial_moving_conduction_candidate` composes this same native solve
with any existing complete `HeatAssembly`, retaining its typed source sample and
errors. Thermal contacts, layers, reservoirs, radiation and finite coolant paths
therefore refresh at the same endpoint as current and electrical contact heat.
Thermal interface resistance remains distinct from electrical contact resistance;
this composition adds neither a second caloric inventory nor a second heat solve.
`tests/spatial_conduction_enthalpy/contacts/motion.rs` checks independent current-
driven two-cell temperatures, separate thermal and electrical contact laws,
subdivision, motoring/generating work, zero-motion parity and atomic refusals.

Signed boundary electrical work and Lorentz mechanical work are retained
separately. The shared driver receives electrical input minus mechanical
output, and independently closes that against Joule heat and finite thermal
storage. Short-circuit magnetic braking therefore heats from mechanical input
with zero electrical supply; generation can export electrical work while still
heating. Mechanical work is neither a thermal source nor numerical dissipation.
Voltage and distributed-current traces retain their native boundary semantics;
open components screen conservative motional fields through their actual
voltage solution.

This bounded path holds geometry, velocity and magnetic flux density over the
interval. The caller must preserve represented material/source occupancy and
finite thermal inventory; it does not add material advection, induction,
momentum evolution or magnet energy storage. Finite-resistance moving contacts
refuse without a surface-work owner. Every discarded, out-of-domain, over-budget
or cancelled trial leaves both accepted owners unchanged. First-order temporal
and spatial accuracy still require refinement independently of energy closure.

`tests/spatial_conduction_enthalpy/motion.rs` checks independent two-cell
backward-Euler temperatures, original-cell heat, changing resistivity and
continued current/material state in motoring and generating regimes. It also
checks short-circuit braking, open-circuit screening, zero-motion parity,
subdivision, source identity and transactional refusals. These remain synthetic
numerical checks, not measured moving-device validation.

This fixed-geometry, memoryless DC closure adds no thermal spatial-resolution
claim beyond its authored finite cells, phase-density motion, contact-layer
storage, opening/arcing, induction, dielectric storage,
unconditional nonlinear convergence or empirical material validation.

## Original latent fractions in constitutive graph feedback

`constitutive_graph::enthalpy::GraphLiquidMassFraction` binds an original finite
cell's borrowed caloric chart to a unique dimensionless graph input. Attachment
checks the actual caloric object, active cell, reference fraction and exclusive
input ownership. Each endpoint samples that chart at the cell's current
specific enthalpy, including changing fraction at constant temperature.
`GraphSpatialLiquidMassFraction` similarly samples the original nodal phase
values through explicit convex weights. Those weights define an observation;
they do not infer a material homogenization law or average unlike enthalpies.

Both existing thermal owners feed temperature, fraction and physical duration
to the same original graph solve. Their existing backward-Euler residual,
source/heat accounting, committed graph history and thermal inventory remain
authoritative. Additional heat assemblies still receive the identical endpoint;
latent storage is already in the caloric chart and is never deposited again.
Input/source/physical-domain, budget and cancellation refusals publish neither
history. Existing callers explicitly select empty fraction bindings.

The source-backed Ohmic consumer closes current, voltage and Joule heating on
these fractions. G1/G3/G4 tests independently eliminate the latent-cell BE
equation for sourced resistivity and conductivity, exercise data-only curve
and surface changes and current reversal, compare matching voltage/current
endpoints, continue the original phase history, preserve the constant-source
limit, and reject foreign caloric objects and late cancelled proposals. A
unit-tetra fixture checks original reference masses and P1 heat deposition
with unequal nodal fractions at equal temperature. Fixed-reference geometry,
explicit source support and the original first-order discretization remain
the model scope; no phase kinetics or numerical/experimental accuracy claim
is inferred from source inspection. Native execution remains unverified here.

## Finite heat on the shared hysteretic gap trajectory

`shared::motion::thermal::ThermalSharedHystereticGapMotion` uses the existing
sealed hysteretic thermal transaction for the complete moving native model.
`SharedHystereticGapHeatInput` retains the held external force and original
mechanical acceptance policy while the common owner clips the actual native
electrical trial allowance. Single/shared stationary callers keep their
original input type. No second thermal or mechanical integrator is introduced.

Each proposed finite temperature refreshes the same original core source query
and winding-resistance sources, then restarts accepted linkages, polarizations,
gap and momentum. One original core pinning heat and the distinct resistor
losses enter their finite caloric receivers exactly once. Complete native
magnetic, gap, leakage, circuit, kinetic and support energy stays in the motion
receipt. Magnetic force work is internal; external electrical work and held
external-force work supply the combined thermal/mechanical/electrical account.

The original `MaterialNetwork` and full `HeatAssembly` retain material heat
conduction, cooling, contacts, radiation and finite-coolant ownership. Source
refusal, trial exhaustion and cancellation preserve every accepted history.
Core fit coefficients stay constant over their explicitly source-admitted
temperature interval. The separately retained gap medium, body inertia and
frozen support sources do not silently follow the core temperature; live
thermoelastic support and unmeasured magnetic thermal fits remain outside this
owner. Native common-flux reversal cuts and first-order thermal refinement
remain required.

G0/G1/G3 fixtures in `hysteretic/tests/shared_motion/thermal.rs` independently
reconstruct complete magnetic/body/support storage, original pinning heat,
per-receiver temperature and R(T) feedback across accepted intervals. Full
heat-assembly conduction/cooling, source-domain refusal, bounded actual trial
work and late cancellation share the same production path. Native execution
remains a separate requirement where workspace dependencies are available.



## Native motion of a shared hysteretic series gap

`material_circuit::hysteretic::shared::motion::SharedHystereticGapMotion` binds
the original shared source/core/gap history to `ResolvedBody` inertia and an
optional existing frozen `AxialSupport`. One bounded vector closure advances
all total winding linkages, physical gap and native body momentum against the
existing circuit and descriptor steppers. Electrical trials consume the
material's exact path-mean winding currents; mechanics consumes its exact
path-mean attractive force. Every trial starts from the same accepted magnetic
memory, circuit stores and body momentum.

Independent terminal-voltage, displacement, electrical-work, mechanical-work
and total-energy gates precede publication. The total account retains physical
leakage, core and gap energy, unrelated native electrical stores, body/support
energy, external electrical/mechanical work, resistor heat and exactly one core
pinning heat. Shared native port preparation and interface Newton closure serve
stationary and moving owners while preserving original branch/winding order,
conserved initial stores and retained capacitor metadata.

Three native rate probes include `sum((r_i/D_i) V_i) - R' lambda_ref velocity`
and reject detectable common-flux reversals. These samples do not certify
unresolved turns; waveform/reversal cuts and interval refinement remain required.
Open gap limits, original source restrictions, incompatible native constraints,
work budgets and cancellation refuse without advancing any accepted history.
This fixed-medium, one-coordinate owner excludes thermal evolution, contact,
fringing, changing area and vector-field force qualification.

G0/G1/G3 fixtures in `hysteretic/tests/shared_motion.rs` compare native motion
with an independent full mutual-inductance ODE and refinement, fixed-linkage
reciprocal work, source pinning and body-density changes, original axial-support
mass/stiffness, and deterministic rollback/retry. Native execution requires the
workspace dependency repositories; authored sources are not device validation.



## Driven solid windings on the existing finite-enthalpy transaction

`WindingThermalNetwork::step_nonlinear_*` consumes driven solid regions through
the existing native circuit/field/enthalpy solve. Original conductivity sources
refresh at each proposed receiving-cell temperature, updating both the current-
constrained winding functional and DC drive. Every trial retains the same
accepted spatial potential, native circuit history and thermal state. Induction
uses the trial functional paired with the retained field increment; changing
the functional alone cannot create voltage in a stationary biased field.

The original field's total spatial E/J receipt supplies solid and floating-
conductor heat exactly once, including local drive/eddy cross terms. Distinct
stranded copper and graph-resistor heat retain their original routes. The
source-owned `HeatAssembly` and magnetic caloric owners remain authoritative;
source-supported linear permeability feedback retains its reversible magnetic
heat and internal-energy account. Numerical field loss is not deposited as heat.
The outer circuit-work limit now bounds all actual nested field trials, including
discarded interface/Jacobian evaluations and final replay, rather than counting
each entire nonlinear solve as one operation.

G0/G1/G3 fixtures in `tests/nonlinear_winding_enthalpy/solid.rs` independently
reconstruct region currents, P1 E/J/heat and induced voltage across successive
conductivity-profile refreshes. They check finite receiving temperatures,
original material conduction, affine permeability caloric exchange, exact
reported work allowances and deterministic retry after cancellation/refusal.
The model retains the original fixed geometry, source domains and first-order
coupling; it does not infer thermal evolution of nonlinear hysteresis fits.



## Source-dependent directional histories on sustained mixed trajectories

`finite_solid::contact::step::sustained` now shares ONE force/contact/drift and publication owner across circular, explicit directional, and original-source directional contact. The elliptic wrapper resolves each current original scene row at the actual post-force velocity, preserving pre-force Newton restitution targets. The source wrapper binds each stable feature to an original ordered interface, explicit static/kinetic coefficient role, and one selected material director. Full rigid inertia, consistent solid cross-point response, physical impulses, separate orbital/spin residuals, supports, work, geometry, clocks and noncancelling split-defect gates remain in that shared core.

`FiniteSolidSourceTrajectoryHistory` pins the original solid preparation, prepared scene namespace and accepted native clock. Immutable scene/index/source entries use copy-on-write checkpoints; free-flight retries clone Arc handles, and only active source entries clone their retained states. New features initialize through their original owner, while re-entry retains that feature's history. Changed source identity, material attachment, coefficient role, scene or stale clock refuses instead of silently rebinding. Histories do not transfer between unrelated features.

The explicitly selected R/F material director defines a convected ORTHONORMAL tangent-frame model. Native texture and slip components are retained exactly in those convected coordinates through the original identity/state admission, avoiding repeated trigonometric roundoff. The original minimum-normal transport supplies the physical director-twist receipt and preserves antipodal-normal refusal. Stretch does not rescale coefficients or accumulated slip; relative torsion is not treated as common material spin.

Original source coefficients are sampled once per geometrically active row before the one simultaneous contact solve. After native drift and whole-path admission, loaded nonopening zero-target contacts may advance source/native history by the declared first-order held post-contact slip h*v. Positive restitution targets remain rebound even below numerical velocity tolerance. Rebound/unloaded events accumulate no free-flight contact travel. Sustained history additionally requires the same endpoint feature and a directed whole-path nearness bound, separate from nonpenetration. An unresolved sustained path rejects the complete private attempt and uses the existing bounded half-step retry. Source failures, late physical audits and cancellation publish neither body nor source candidates. Explicit conditions supply temperature/pressure; neither an impulse nor an integration defect supplies heat or contact pressure.

The geometry owner follows convex endpoint barycentric witnesses through each moving triangle, using the quadratic Bernstein control hull. Its relative edge form removes common translation and avoids the artificial O(sliding-distance squared/radius) clearance increase of a fixed starting witness. Native post-contact point velocity includes actual rigid spin; no COM-only travel surrogate is used. The finite-slip quadrature and source-coefficient holding remain first-order discretizations, not exact curved-slip integration or a general Moreau-Jean/existence theorem.

Focused fixtures cover full-tensor circular equivalence, analytic anisotropic spin response, actual speed-dependent slip/path recurrence, force-budget refinement, late source refusal and atomic history rollback, scene/clock pins, native separation/re-entry, elastic and sub-tolerance rebounds, and convected rotating texture. Pinned rustfmt and scoped whitespace checks pass. Independent algebra checks on1200 deforming paths (49200 samples) confirm the Bernstein identity/bound, including one metre of steady planar sliding; these are not Rust-kernel tests. Native compilation and the new Rust regressions remain unexecuted under the established full-disk/dependency-memory and constellation-preflight block. Broad nonsmooth-contact and flagship beads remain open.

## Material-attached directional contact frames

`finite_solid::contact::texture` carries one explicitly selected participating material line with the original rigid orientation R or solid TOTAL deformation gradient F, then projects it into the actual contact plane. Wrong endpoint kinds and unresolved projections refuse. The original solid owner now exposes a sparse prepared deformation-gradient observation using exactly four cached node identities and the same translation-invariant arithmetic and positive-det-F admission as its force/tangent assembly. Coupling reads only twelve native displacement components, including prescribed DOFs, without allocating a full-mesh displacement buffer or evaluating a constitutive law.

This is an orthonormal single-director orientation model: stretching/shearing can turn the director but does not infer strain-dependent friction coefficients, average two surface textures, or transform an ellipse metric. Physical director evolution remains distinct from objective constitutive-history transport. The original mixed solver retains its full inertia, consistent solid mass, physical impulse/work, and angular-momentum gates.

Focused fixtures cover physical quarter-turn anisotropic recoil, full world-rotation covariance, explicit first/second endpoint attachment, affine sheared solid motion with prescribed displacement and large translation, sparse original-node sampling, inversion/degeneracy and cancellation. Pinned rustfmt and scoped whitespace checks passed; these new native fixtures have not been executed because the available environment has no free disk and its native dependency build is memory-blocked. No new native pass or broad contact-bead completion is claimed.

## Sustained mixed contact through original force and kinetic subflows

`finite_solid::split` supplies the original prepared hyperelastic, held dead-load,
explicit uniform-pressure and bound stationary-gravity forces to the existing
consistent-mass impulse owner. After a separately
admitted mixed contact event, the original generalized-alpha owner advances the
mass-only kinetic subflow and records its genuine solve/time/history receipt.
The original mass and force owners recompute smooth endpoint acceleration;
contact impulses are not disguised as smooth acceleration. An immutable private
receipt binds the original material, configuration, time, history, loads and
work, even if a caller edits public diagnostic fields.

`contact::step::sustained::sustained_sphere_trajectory_candidate` composes this
solid path with native full-inertia rigid force/torque impulses, ONE mixed contact
solve and original torque-free drift. Incoming pre-force velocities determine
Newton restitution targets, preventing a held load from creating a fresh bounce
at an otherwise resting contact. Original current-arm stationary-support moments
and world orbital/spin momentum are retained. Kick orbital error, kick spin
error, drift orbital increment and drift spin increment have separate absolute
bounds before accumulation; neither within-body nor cross-event cancellation
can hide lost angular transport.

Source, internal-force, support, kinetic and strain transfers remain distinct.
Source/internal impulse and work decomposition roundoff is explicitly measured.
Force-displacement work and body-torque kick work use the declared first-order
convention. Signed and absolute source/elastic/rigid splitting defects remain
visible; contact-map loss can contain splitting loss and is not friction heat.
The whole-call absolute defect budget is allocated by interval duration. Excess
splitting error, unresolved affine geometry, original trial-material refusal or
rotational nonconvergence can shorten and recompute the complete original map
from the accepted state. Admitted step sizes are retained, with conservative
work-margin growth, avoiding repeated rediscovery of a stable contact step.

Whole affine sphere-pair, deforming-triangle and element-orientation admission
is shared with the transient driver; every negative actual witness gap refuses.
Circular, explicit directional and source-dependent directional dispatches share
this owner. Uniform pressure enters through the explicit entries below. This is
a bounded FIRST-order split requiring refinement, not general Moreau--Jean,
continuous forced Euler-top quadrature, true-flow CCD, contact heat partition, or
guaranteed resolution of grazing/Zeno behavior.

Eight independent fixtures cover finite stretch, original mass/current-arm force
transfer, immutable continuation and atomic refusal, loaded rest with restitution
one, a loaded mixed chain, applied torque, and genuinely deforming face recoil.
The resting fixture checks that the declared splitting loss halves with the step;
the deforming fixture independently uses the consistent face-mode mass 0.6 kg
and fixed-apex reaction. Pinned formatting and source/math review passed. Native
compilation and execution remain unverified because the pinned sparse dependency
preflight and shared storage/memory constraints block the native test path.

## Held uniform pressure in sustained finite-solid contact

The explicit `force_kick_candidate_with_pressure`,
`sustained_sphere_trajectory_candidate_with_pressure`,
`sustained_elliptic_sphere_trajectory_candidate_with_pressure` and
`sustained_source_elliptic_sphere_trajectory_candidate_with_pressure` entries
extend the same original first-order force/contact/kinetic-drift owner. One
finite signed uniform pressure is held on the body's original closed material
boundary for the whole call. Legacy entries retain their explicit-pressure
refusal; explicit zero pressure takes the unchanged no-pressure numerical path.

The kick evaluates the original instantaneous pressure force at the SAME current
configuration used by the consistent-mass impulse and simultaneous contact maps.
Closed-pressure force and torque resultants pass independent impulse gates.
Pressure, elastic and dead-load impulses enter one native mass solve, retaining
actual stationary support impulses and their current-arm moments. A midpoint
discrete-gradient pressure force is deliberately not substituted into this
old-geometry split. The immutable continuation receipt retains the pressure
surface and value independently of editable diagnostic fields.

The original pressure owner supplies `-p Delta V` exactly once in displacement
source work; fixed-position impacts contribute no pressure work. Original endpoint
acceleration includes current pressure force. Pressure kick work, pressure
displacement work and their signed splitting defect are separately visible.
The existing refinement gate sums absolute dead-load, pressure and elastic
splitting defects separately, so those errors cannot cancel or become heat.
Circular, directional and original-source dispatches retain the same pressure
through every shortened attempt; body candidates, accepted clocks and private
source histories publish together after the existing final gates.

This is still a first-order split requiring refinement. It does not admit
nonuniform follower pressure, supply an exact pressure flow, partition contact
loss into heat, or certify the full contact trajectory. Original topology,
geometry, source-history and cancellation admission remain in their native owners.

Five new public regression cases cover pressure force/work/endpoint acceleration,
zero-pressure parity and refusal, circular/directional support balance, refinement
with late step-budget rollback, and late source-history refusal. All nine changed
Rust files passed pinned formatting and static native-owner review. The complete
public fs-couple/Cx target could not compile because the local constellation lacks
required siblings, including asupersync; these five public cases were not executed.

A separate selected native numerical lane executed the original pressure
force/volume/work arithmetic, native consistent-mass FGMRES, contact solve,
generalized-alpha kinetic drift and rigid momentum maps on a canonical unit
tetrahedron. Its pressure extraction omits Cx parameters/checkpoints and topology
construction, supplying that known topology explicitly; public callback, retry,
cancellation and topology admission are outside this execution claim. A free
recoiling pressure-loaded contact balances work within 5.6e-16 J and momentum
within 3.4e-25 in the fixture; halving the interval reduces the splitting defect
approximately fourfold. A supported contact over one fixed interval with 4/16/64
native steps has absolute defects 2.99017e-4/7.47559e-5/1.86891e-5 J, with work
closure below 1.2e-15 J and support-inclusive momentum residuals below 4.2e-17.
These analytic fixtures are refinement evidence, not experimental validation.

## Automatic complete mixed sphere-contact trajectories

`finite_solid::contact::scene::FiniteSolidSphereScene` binds one COM-centred
sphere per native rigid body and one explicit source-resolved law per canonical
sphere pair. Pair laws are never inferred by mixing each sphere's solid-interface
coefficient. Original closed outward AABB broad-phase bounds prune rigid pairs;
original current-triangle BVH witnesses supply the solid contacts. Both analytic
centre gaps and actual reconstructed material gaps must be nonnegative. Even
surface IDs and odd canonical pair IDs retain disjoint, input-order-independent
feature identity within the fixed prepared scene.

`mixed_sphere_trajectory_candidate` passes every detected rigid/rigid and
rigid/solid row into the same simultaneous response, then uses the original
objective energy/momentum solid and full-inertia rigid smooth owners. The
explicit `mixed_sphere_trajectory_candidate_with_pressure` entry preserves
uniform pressure on the original closed solid boundary throughout retries. Whole affine paths check every rigid pair,
swept triangle candidate and original tetrahedron orientation. An unresolved
possible encounter shortens and recomputes original native motion; endpoints are
never projected. The surface-only compatibility entry still refuses rigid pairs
because its API supplies no pair law. This is transient event refinement with
explicit positive activation distance, not sustained contact or true-flow CCD.

The original rigid sphere geometry now has one reusable paired-material witness.
The shared interval sphere-pair path retains both original radii in its exact
axis-separation fallback: real touching/receding pairs can be admitted while even
one representable step of penetration refuses. Native local residuals and whole
event finite-gap angular/work defects accumulate without cross-event cancellation.
Independent fixtures cover elastic pair trajectories, simultaneous solid-supported
chains, full rotated inertia matrices and original source-history routing. Native
compilation/execution is still required; source/format checks alone are not proof.

## Simultaneous mixed rigid and finite-solid impacts

`finite_solid::contact::mixed` assembles rigid/rigid and rigid/solid rows in one
bounded original Coulomb or directional-ellipse solve. Signed full-inertia rigid
incidences retain every shared-body interaction; original supported consistent
solid mobility retains every cross-point block. Rigid-only rows introduce no
solid mobility columns. Existing rigid/solid entry points adapt to this owner.

Each rigid body receives one native resultant linear and COM angular impulse.
Actual paired material points remain distinct: a finite-gap tangential couple
must pass the world angular-momentum gate, never be hidden by moving witnesses.
Original endpoint velocities, stable kinetic/relative impulse work, passivity,
solid support reactions, and noncancelling sums of per-body residuals gate the
entire candidate. Both spin and orbital angular momentum use their actual frames.

Mixed source-backed directional events preserve original input indices, sample
actual preimpact slip in the solved frame, and return history only with admitted
mechanics. No impact duration, pressure law, heat partition, penetration repair,
continuous geometry, or transported material texture is inferred. Independent
fixtures cover coupled free/supported chains, rigid-only chains, directional spin,
input-order invariance and atomic refusal. Native execution remains required;
format/source review is not a test-suite pass.

## Simultaneous deformable-solid contact assemblies

`finite_solid::contact::solids::contact_candidate` resolves caller-supplied
material-point contacts among one to eight original prepared finite solids at
one common physical time. Bodies may retain different source materials,
reference densities and stationary supports. Up to 21 contact rows use one
original simultaneous Coulomb/restitution solve; contact IDs determine stable
row order while returned bodies retain input order.

Signed material-point traces use each body's native supported consistent mass.
The assembled mobility retains every directed cross-contact entry and all four
signed blocks for distinct witnesses on the same body. Each body receives one
summed native point-impulse candidate, preserving shared-body coupling and
quadratic work. There is no lumped point mass, pairwise restitution sweep,
numerical symmetry repair or replacement mass/contact solver.

Acceptance checks actual endpoint point velocities, native quadratic kinetic
increments, relative impulse work, nonnegative loss and total linear/angular
momentum about the caller's origin, with stationary-support reactions accounted
separately. Absolute sums of per-body momentum/work residuals prevent
cancellation. The independent `absolute_contact_couple_kg_m2_s` gate sums
`|(x_first - x_second) cross J|` across pairs: opposing finite-gap friction
couples cannot hide in the global angular balance. Witnesses are never relocated
to repair that balance.

Candidates are returned together only after final admission and cancellation.
Sources, configuration, acceleration, time and accepted-step history are
unchanged by this instantaneous event. Contact discovery, trajectory integration,
moving supports, penetration repair and per-contact thermal allocation remain
separate concerns. Supplied witnesses must lie on the same normal line within
the declared geometric tolerance; the original contact law owns activation,
restitution, nonuniqueness and passivity admission.

Seven numerical cases passed with actual native CSR, complete consistent-mass
initialization/FGMRES, the original contact solver and the exact signed-scatter
helper. They cover unequal densities under a common boost, sliding loss, a
shared middle body, same-body touching arms, a supported counterpart and single
or globally cancelling finite-gap couples. The full Cx-bound public candidate
and five new public regression functions were statically reviewed but could not
be type-checked or executed without the missing material/runtime dependency
closure; these native numerical cases are not a full public-target Cargo pass.

## Spatial layer heating and native laminate deformation

`laminate_thermal::LayeredThermalShell` binds the original sourced thermal stack
to the original MITC shell through `LaminateThermalSource`. Mechanical and
thermal inputs must agree on each card, reference thickness, typed fixed context,
patch area and per-cell reference mass. Every thermal control volume becomes
one sublamina at its actual temperature, with current pinned elastic properties
and the original integrated expansion path. The native 128-lamina limit refuses
unsupported resolution before subdivision; temperatures and CTEs are never
averaged to fit that limit.

Current source-backed sections drive native load assembly, bubble condensation,
compatible static displacement, physical stress resultants and stored energy.
Mesh, layer orientation, reference geometry and all mass moments stay fixed.
One prepared laminate is shared across facets; each field supplies its exact
cached oriented section. The resulting native pencil is exposed for existing
vibration consumers. `LaminateShellSnapshot::membrane_prestress` supplies the
native bounded planar initial-stress pencil from that same compatible state;
vibration-history transfer remains a separate owner. Accepted strain polynomials
and the exact thermal load are retained for reusable work/observation queries.
Explicit translation/rotation limits and whole-element reduced elastic-strain
bounds enforce the caller's linear-model envelope, including the thermally
driven internal director field.

`advance_to` uses the original sparse implicit layered heat owner and publishes
heat plus mechanical endpoints only after both succeed. Source-domain exits,
mechanical refusal and cancellation retain the original accepted thermal state.
This is one-way quasistatic heat-to-mechanics coupling: the reported heat audit
excludes elastic entropy, stiffness-change work and mechanical feedback.
Mechanical energy is separate, with no claim of reciprocal energy conservation,
finite rotations, transverse-normal/interlaminar recovery or measured accuracy.
`tests/laminate_thermal.rs` exercises native two-cell heat exchange against an
independent CLT solution, five-cell gradient retention, compatible expansion,
reference mass and late atomic refusal.

### Reciprocal finite heat and layered-shell equilibrium

`laminate_thermal::reciprocal::ReciprocalLayeredThermalShell` adds mechanical
feedback to the SAME original caloric timeline and compatible native shell.
One declared temperature interval per actual thermal cell pins the original
plane-stress potential and integrated expansion. Native Bernstein controls and
both layer faces admit the complete total/elastic principal strain envelope.
The original piecewise-affine caloric capacity plus a conservative full-interval
in-plane contribution must remain positive. With varying elasticity, the
resulting entropy-derivative diagonal must also dominate the whole-temperature-
box absolute row sums of the nonlocal shear entropy Hessian. This sufficient
Gershgorin bound admits a positive coupled entropy-capacity matrix throughout
the declared envelope. It can conservatively refuse broad boxes; callers may
declare a physically justified narrower temperature/strain domain. Pointwise
capacity samples never replace that admission. Caloric/source population and work
constraints must describe the same unloaded material path; card identity alone
does not establish that experimental compatibility.

Original elastic coefficients and instantaneous CTE may vary affinely, with
every cell retaining its independent temperature and source-integrated thermal
strain. The source-derived Q(T) and corrected shear H(T) need not themselves be
affine. Each varying-material proposal rebuilds the native section, equilibrium
operator, thermal load and recovered internal director displacement. Reference
mass, geometry, support constraints and the original caloric owner stay fixed.
Exactly constant elastic sources keep their cached native factorization and
allocation-free quadrature shear path.

At the original triangle quadrature and two points through every layer, endpoint
average stress times compatible strain change gives mechanical work. The SAME
potential supplies `S=-Psi_T` and `U=Psi+T*S`; in-plane reversible heat is
mechanical work minus the change of that internal energy. The independent free-
energy identity retains both integrated-expansion work and polarized modulus
work `V*e0^T*(Q1-Q0)*e1/2`. Neither becomes dissipative heat by itself.

The corrected whole-stack shear potential supplies one internal energy
`U_s=Psi_s+sum(T_i*S_s,i)`. Analytic nonlocal gradients provide its entropy
conjugates. Symmetric spatial-order/reverse-order coordinate sweeps partition
`H1-H0` into `delta_H_i`; physical shear heat to original cell i is
`-A*g0^T*delta_H_i*g1/2-delta(T_i*S_s,i)`, integrated over the native quadrature.
No stored shear energy or solver residual is arbitrarily split among cells.
Section derivatives and coordinate increments are prepared once per endpoint
pair, outside facet quadrature. These finite-step paths require refinement and
do not establish exact continuum entropy transport or arbitrary layer-order
invariance. Native physical energy and nodal work independently cross-check
the complete in-plane and shear integration.

Numerical drilling storage and work remain separately reported. When its
stiffness changes, the independently polarized numerical parameter work
`(q0·fg1-q1·fg0)/2` enters only the numerical/external account. It never heats a
material cell. The combined gate closes caloric plus plane/shear internal
storage plus gauge storage against external heat, external mechanical work and
that explicit numerical parameter input.

Bounded outer iterations restart the original sparse implicit conduction solver
from accepted heat history. Interval-only reversible joules update through its
existing source port. Constant elasticity shares one scaled LDLT factor; varying
elasticity solves the actual source-resolved operator. Publication requires local heat consistency,
summed absolute heat defects, physical/gauge work closure and a combined energy
gate, followed by a final cancellation checkpoint. Failure leaves both accepted
histories unchanged; stale temperature/time pairing refuses. The report separates
external heat/work from reversible exchange and exposes each original cell.

This is quasistatic, small-strain endpoint coupling with explicit iteration and
energy budgets. It needs time/space refinement and makes no exact discrete
isentropy, dynamic-history, phase-flow, finite-rotation or moving-transport claim.
`tests/laminate_thermoelastic.rs` checks independent adiabatic-traction cooling,
restrained-heating capacity, the accepted thermal-prestress pencil, nonuniform
bending with nonzero transverse work, energy closure and atomic refusals.
The existing one-way fixtures share source/mesh helpers without changing their
CLT oracles. Targeted compiler metadata checks cover these entries; numerical
execution and experimental material validation remain separate requirements.
`tests/laminate_thermoelastic_softening.rs` adds independent restrained-heating
cubic and constant-traction references, accepted-temperature prestress, a
two-layer through-thickness shear oracle, source/capacity refusals and late
cancellation. Those cases exercise stiffness changes through the original
finite heat cells without replacing their mass or energy datum.

### Reciprocal finite heat with retained shell inertia

`ReciprocalLayeredThermalShell::begin_transient` transfers the same admitted
equilibrium, original source intervals and caloric history into
`laminate_thermal::reciprocal::dynamic::DynamicLayeredThermalShell`. Explicit
nodal and residual-director velocities initialize the full native consistent
mass. Both internal director coordinates per facet then evolve independently
in the invertible native `y=Tq+Br` basis. Observation uses their actual strain
polynomials and stress; no later heat trial resets them to a static thermal
lift. The full INITIAL coordinate basis, physical mass, reference geometry and
stationary supports stay fixed. Original affine elastic coefficients and CTE
may vary with each cell's actual temperature; the parent's complete-box capacity
and whole-field strain admission remain required. Updating material stiffness
retains all nodal/internal cross terms in that original basis, so no moving-
basis lift, momentum remap or inferred temperature-rate actuator is introduced.
The original numerical drilling stiffness also stays fixed throughout dynamics.

Every outer trial restarts heat and motion from the same accepted histories.
Constant elastic sources reuse one midpoint effective sparse factor across all
trials and accepted intervals with the same dt. Varying sources assemble and
factor their actual proposed material operator in the original basis; the
momentum equation uses both accepted and proposed endpoint stiffness forces. The
shared source-potential quadrature computes actual nodal AND internal strain
work; explicit internal-coordinate forces independently check that work.
In-plane reversible cell heat stays `W_stress-delta(Psi+T*S)`; the same nonlocal
shear entropy and symmetric coordinate-temperature work used by the quasistatic
owner supply the shear contribution. The combined gate includes
native kinetic storage, physical corrected shear and separately named numerical
drilling storage, without adding any of them twice or depositing gauge work as
heat. The current consistent mass contains physical translation/director
inertia, assembled directly from the original material velocity field. The
modal pencil's numerical relative-spin inertia is not counted as physical heat
or kinetic energy in this transient owner. Numerical drilling force/storage
come from the explicit fixed reference gauge, rather than a current condensed
static operator that would omit actual internal forces. Its parameter work is
exactly zero. The native diagnostic separately reports physical stiffness-
parameter work; it is not added again to the source-derived heat account.
Compatible initial gauge displacement/velocity and absence of effort on free
massless director-axis directions are required by the native dynamic owner.

Physical strain, complete director motion, source support, exact accepted
time/temperature pairing, native mechanical balance and local/global heat
defects all pass before the final atomic publication. Failure or cancellation
leaves accepted motion and caloric history unchanged. The report retains the
original heat audit, per-cell exchange, native midpoint identity, kinetic
change and internal-director work. No finite-rotation, prestressed dynamics,
moving transport, material damping or exact discrete-isentropy claim is made.
Time/space refinement remains required.
`tests/laminate_thermal_dynamic.rs` supplies independent scalar-oscillator,
adiabatic inertial-cooling and analytic internal-inertia oracles, nonuniform
heating and transactional refusal checks. These synthetic implementation
fixtures are not experimental material validation.

`tests/laminate_thermal_dynamic_softening.rs` adds an independent heated scalar
oscillator with simultaneous momentum and finite-heat equations across several
changing-load steps. A heterogeneous laminate fixture checks original reference
kinematics and exact physical bubble inertia while the current condensation map
changes, nonlocal shear heat, full storage balance and late rollback. The native
stiffness-parameter diagnostic is checked separately from the thermodynamic
exchange. These authored references require numerical execution and timestep/
space refinement before making accuracy or experimental claims.

## One shared hysteretic core in a native multiwinding circuit

`material_circuit::hysteretic::shared::SharedHystereticCircuit` connects the
source-owned `SharedHystereticCore` to multiple original winding branches.
Selected inductance placeholders become private native current ports; all
physical leakage and the one common polarization history remain in the material
owner. Original branch/node identities, external source controls and index gaps,
resistor heat routes, unrelated R/L/C stores and caller-supplied winding order
are retained. One bounded native Newton-Krylov solve closes every
`delta(lambda_i) = dt V_i` through the original `CircuitDae`.

Every trial starts from the same accepted total linkages and common memory.
The material owner's exact linkage-chord mean currents supply native electrical
work, admitting pinning and saturation crossings without inconsistent endpoint
magnetomotive balance. Actual endpoint currents are used for native consistency;
incompatible ideal-current cutsets refuse without publication. The step reports
terminal voltages, endpoint and mean currents, complete leakage/core storage,
and exactly one pinning-heat account. Private-port work cancels the complete
material work; actual resistor heat remains separate. External supplied work
closes independently against total stored change and physical heat. Numerical
residuals are never relabeled as heat.

A monotone common reference flux is required per interval; individual winding
linkages may move in different directions. Three native start/mid/end voltage
probes check the sign of `sum((r_i/D_i) V_i)`, whose sign follows the common
flux rate by the positive load-line tangent. Detectable interior turns return
`SharedMagneticReversal`; those samples do not certify unresolved oscillations.
Callers must retain waveform/reversal cuts and refinement. Original constant-fit
source queries may refresh without changing any history; changing coefficients,
moving geometry, eddy/excess-loss fits and empirical qualification are not inferred.
Actual thermal deposition and temperature publication require a thermal owner.

Five authored-source regression groups cover the full native mutual-L linear
limit with unrelated RLC and nonascending branch order, an independently
eliminated two-winding RL/pinning-crossing solution, one shared major-loop area
and remanence, physical polarity reversal, and source-domain/budget/cancellation/
turning-point refusals with unchanged accepted history. Native test execution
remains a separate requirement where workspace dependencies are available.

### Finite core and winding heat with retained shared magnetic memory

`shared::thermal::ThermalSharedHystereticCircuit` and its
`SharedHystereticEnthalpyCandidate` reuse the original hysteretic/enthalpy owner
through sealed native history and step interfaces. The public single-core types
retain their original default backend. Every thermal trial refreshes the same
original core source at its actual caloric temperature and restarts the same
accepted shared winding linkages and polarizations. No discarded thermal or
electrical trial advances the accepted magnetic memory.

Exactly one pinning-heat receipt enters the whole fixed solid core cell; each
original resistor's Joule heat enters its distinct source-bound material cell.
Leakage storage stays reversible. The original `MaterialNetwork` and complete
`HeatAssembly` cooling, thermal-contact, radiation and finite-coolant paths retain
their existing source, cancellation, budget and total-energy checks. Winding
resistance changes with its actual source temperature, feeding back into the
shared magnetic trajectory. Reports preserve all original winding states and
the complete shared magnetic work/heat receipt alongside the finite thermal state.

The magnetic fit remains constant over its explicit source-admitted temperature
interval. Geometry, density and whole-cell core volume stay fixed; common-flux
reversal cuts and the existing first-order thermal coupling remain authoritative.
No unmeasured temperature-dependent hysteresis fit or empirical loss map is added.
Three focused authored-source regressions check per-receiver finite heating,
original source and resistance feedback, retained memory across intervals,
native heat-conduction feedback, and unchanged histories with deterministic retry
after late source-temperature refusal or cancellation.



## Shared conservative gravity on deforming solids

`solid_gravity::bind_field` binds the original `fs-mbd::gravity::PointGravityField`
and an explicit inertial frame directly to the original hyper3 gravity callback.
The field owns specific potential, acceleration and its analytic spatial
derivative. The solid evaluates them at current quadrature positions and owns
reference material mass, consistent loads, external-force tangents, equilibrium
and support reactions. The adapter allocates nothing and introduces no second
gravity, material, mass or solver owner. Field-domain failures propagate and
cancellation retains its native error. World coordinates must agree explicitly;
moving-source work, rotating-frame forces and self-gravity are not inferred.

`tests/solid_gravity.rs` exercises the shared central field through current solid
force/potential/tangent differences and the shared uniform field through native
static reactions and stretch-invariant reference weight. These are authored
analytic/numerical checks, not measured material validation.

## Finite water mass, phase inventory and thermal cavitation

`cavitation::water_cell::EquilibriumWaterCell` evolves total water mass and
internal energy in a rigid volume. The original `IsochoricWater` owner recovers
temperature, pressure and finite liquid/vapor inventories from `V/M` and `U/M`
using the retained IF97 Region 1, 2 or 4 state. Conserved energy is never
replaced by an approximate source inverse. Phase reconstruction, source energy
and volume defects, and the accepted first-law residual have separate budgets.

Incoming and outgoing `WaterParcel` values retain original IF97 enthalpy and
entropy. Flow work is included once through enthalpy; rigid-boundary work is
zero and no additional latent-heat source is invented. `transfer_water` prepares
both finite cells, checks total mass/energy and cancellation, and publishes both
or neither. Withdrawal uses the initial well-mixed donor state, an explicit
first-order partition requiring interval refinement. Gross opposing flows count
against the exchange cap. Hydraulic flow direction and amount remain the
caller’s responsibility; the cell is not a momentum, valve or choking solver.

`water_cell::acoustics::FrozenWaterCellAcoustics` borrows one accepted snapshot
and reads its exact retained IF97 endpoints. Its existing `FluidAcoustics`
implementation feeds native ducts through single-phase response or the original
Wood compliance with mass-to-volume conversion. Acoustic perturbations freeze
phase masses and assume common pressure/velocity, adiabatic individual phases
and long wavelengths. They do not use an equilibrated boiling sound speed or
infer a resonant bubble radius from vapor volume.

`EquilibriumWaterCellAcoustics` provides the separate rapid-equilibration
limit. The reusable material `phase::water::acoustics::HomogeneousEquilibriumWater`
uses the original saturation pair, analytic Gibbs derivatives and Region 4
pressure slope. Its phase-fraction tangent satisfies original `du + p dv = 0`,
and its compressibility includes the resulting phase exchange. Original IF97
Gibbs/Clapeyron matching defects and the raw entropy tangent remain visible;
no source data are shifted or declared entropy production. Production evaluation
needs no inverse solve, finite differencing or workspace allocation.

Single-phase accepted water retains its exact original acoustic response.
Two-phase equilibrium observation requires both phase masses strictly positive;
an exact saturation endpoint refuses the ambiguous two-sided tangent. Caller
frequency and amplitude must admit rapid equilibration, negligible slip and an
interior coexistence perturbation. There is no guessed relaxation frequency or
interpolation between this response and frozen-phase acoustics. Native checks
compare the analytic tangent with independent reversible-work compression and
exercise distinct propagation and heat-driven updates from the same inventory.

`cavitation::water_cell::thermal::advance_water_wall` couples that finite
water inventory to an original `EquilibriumPhaseState` and enthalpy chart.
One bounded native heat-feedback solve debits Q from the finite wall and adds
the identical Q to conserved water U. The original `LayerStackConduction`
evaluates k(T), contact resistance and geometry at both proposed endpoint
temperatures. Both owners publish only after source recovery, energy and
temperature budgets, interface convergence and the final cancellation poll.
Conserved storage, IF97 source-inversion defects, conducting-layer continuity
and heat-root residuals remain separate. Entropy retains the native chart and
raw IF97 endpoint changes; the remainder is not certified physical production.
This backward-Euler closure requires uniform water/wall temperatures,
negligible wall pV work, separate quasi-steady conducting layers with negligible
storage, and negligible unresolved fluid-film resistance. It supplies no
invented boiling coefficient or extra latent-heat source. Four new native
regressions cover an independently manufactured finite-capacity endpoint,
nonlinear layers driving evaporation/condensation, repeated energy accounting,
and atomic late-cancellation/root/temperature refusals.

`cavitation::water_cell::volume` supplies sealed prescribed-volume motion on
the same finite water owner. `propose_volume` and `advance_volume` solve
`U1-U0=Q-(p0+p1)(V1-V0)/2` through the original IF97 fixed-v source and
existing bounded scalar solver. The exposed mean absolute pressure supplies
the identical signed work to an enclosing mechanical solve. No liquid EOS,
wall inertia, inferred evaporation rate or extra latent heat is introduced.
A dropped proposal or failed/cancelled solve retains M, U, V and the clock.
The rigid-volume endpoint is exactly the original heat update.

The original phase masses/volumes, conserved energy, source-image defects and
raw entropy remain separate. Explicit volume, temperature, time and nonlinear
work limits admit the uniform quasistatic equilibrium step; they are not
trajectory-error bounds. Trapezoidal pressure work requires path/phase-crossing
refinement, and no heat-entropy or wall-momentum closure is inferred.
Five native tests cover liquid expansion into cavitation and recompression,
isothermal coexistence with original latent enthalpy and pressure work,
single-phase adiabatic entropy convergence, the native equilibrium acoustic
tangent, rigid reduction, deterministic proposals and atomic refusal.

Eight native regressions exercise withdrawal-induced vapor, refill condensation,
conservative two-cell transfer, late cancellation, nonaccumulating inverse error,
original phase acoustics and actual phase-change feedback into duct propagation.
The thermal source domain remains 273.16..=623.15 K with explicit pressure and
iteration bounds. Finite-rate nucleation, metastability, noncondensable gas,
interface energy and resolved thermal gradients are outside this closure.
Heat-boundary temperature is unspecified, so the raw entropy remainder is not
certified entropy production.

## Resonant dilute bubbly-liquid propagation

The opt-in `cavitation::bubbly::DiluteBubblyLiquid` binds the original
`CavitationProperties` source to an identical native `KellerMiksisBubble`.
It consumes that owner's stable reference linearization, including polytropic
gas and capillary stiffness, viscous inertia, radiation damping and the
pressure-derivative forcing numerator. No independent gas EOS, damping fit or
second bubble equation is introduced.

For the exp(-i omega t) fractional-radius susceptibility `chi`, Foldy's
coherent dispersion is `k^2 = omega^2/c_l^2 - 3 alpha rho_l omega^2 chi`.
The native duct momentum relation uses `z = rho_l omega/k`; bubble mass is
neglected to leading dilute order. An admitted response reuses the original
cylinder/cone wave bases, smooth-profile Magnus integration, compliant-wall
application and stable chain composition through the explicit
`fs-duct::fluid::*_with_wavenumber` entry points. The pure-liquid endpoint
is exactly the original wave and chain result.

Admission bounds void fraction, bubble compactness, wavelength relative to
mean spacing, fractional-radius amplitude and the native radius/Mach/pressure
domain. It additionally evaluates neighboring scattered pressure and
multiple-scattering strength: small void fraction alone does not admit a
near-resonant bubble population. The supplied pressure envelope is the maximum
LOCAL amplitude anywhere in the bubbly region; reflected/amplified fields must
respect it, and inlet amplitude alone is insufficient.

This is a frozen, dilute coherent-wave model around stable
equilibrium. It excludes thermal gas relaxation, evolving vapor inventory,
spatial correlations, nonlinear cloud collapse and individual wall-near bubble
dynamics. Coherent attenuation includes scattered energy and cannot be
deposited wholesale as thermal heat. Physical references are
[Foldy (1945), equations 15 and 54](https://doi.org/10.1103/PhysRev.67.107)
and [Commander and Prosperetti (1989), section VI](https://pages.jh.edu/aprospe1/publications/PapersPublished/BubblyLiquids/Commander.pdf).

Focused regression sources check the exact IF97 pure-water endpoint across
cylinders, cones and profiles; independent static gas/capillary compliance;
native nonlinear-to-linear response including zero-viscosity radiation;
matched-pipe phase/attenuation; and actual amplitude, scattering,
homogenization, source-mismatch and cancellation refusals.

`PolydisperseBubblyLiquid` consumes an original `fs-flux::BubblePopulation`
without replacing its classes by a mean radius. Each occupied class retains
its native Keller--Miksis gas inventory and stable equilibrium linearization;
all classes share the exact liquid and the same ambient pressure to roundoff.
The coherent compliance is `3 sum_j alpha_j chi_j`, with
`alpha_j = (4 pi/3) n_j R_j^3`. Scattered power uses the separate positive sum
`sum_j n_j 4 pi |f_j|^2`, so opposite radial phases cannot cancel scattering
admission. Mean spacing comes from total number density; neighbor admission
uses the largest individual scattering amplitude. Splitting a class into
identical bins cannot relax these physical limits.

Mono- and polydisperse routes share one native susceptibility/admission kernel.
Each occupied radius must independently pass its nonlinear owner's radius,
Mach, pressure and amplitude domain. Fixed-order compensated accumulation
preserves number, volume, complex response and scattering moments; harmonic
evaluation allocates no workspace and polls cancellation per class. The caller
sets an explicit class budget. Empty populations recover original pure-liquid
propagation. This frozen distribution does not evolve nucleation, transport,
coalescence, breakup or bubble-cloud interactions. Six native population and
seven coupled regressions cover the true two-radius response, independent
scattering sums, one-bin and split-bin identities, static compliance, actual
duct phase/attenuation, and refusal of a rare but overdriven radius class.

## Event-resolved vapor cavities in liquid ducts

The opt-in `cavitation::duct::pockets::CavitatingDuctPockets` retains an original
material `DgDuct`, its complete native state and wall memories, the material
source, and one accepted clock. Up to eight declared endpoint/interior ports
may rupture into discrete isothermal pure-vapor pockets. Native closed,
pressure-release and matched boundaries remain original. `CavitatingDuct`
delegates its single-endpoint API to this same solver. An interior station uses
the original shared pressure and exact flow jump `Q = q_right - q_left`;
there is no second liquid inertia or invented hydraulic solver.

For original energy coordinates `x`, generator `L` and active inflow columns
`B_A`, every active pressure is constrained together:
`(B_A^T B_A) Q_A = -B_A^T L x`, `x' = L x + B_A Q_A`, `V_j' = Q_j`.
Filled ports have zero inflow. The existing diagonally scaled Cholesky solve
prepares the joint flow gain. Original mass-map observations and active gains
are cached; the RHS allocates no workspace. Original `fs-time` multiple dense
events locate pressure-driven rupture and volume-driven collapse on one
trajectory. Overlapping root brackets must admit one common contact state;
a bounded contact-complementarity solve chooses all zero-volume branches
together. Ambiguous flow and unresolved grazing refuse. Neither liquid state
nor volume is reset or projected. Small signed localization residuals remain
visible under explicit pressure, volume and event-time tolerances; finite event
scans do not certify unresolved tangencies or multiple crossings.

The original duct adopts independently integrated pressure impulse and actual
`p Q` work at every port, plus original generator dissipation. Its own trajectory
admission and combined storage `H_duct + sum_j (p_ref-p_sat) V_j` close before
publication. Refusal, cancellation and exhausted event/attempt budgets preserve
the complete accepted liquid/cavity history.

The bounded dense lane admits at most 128 complete coordinates and eight
pockets, with every native port explicitly assigned. It checks source/acoustic
identity, liquid strain, cavity size and available vapor-density data. Interior
Mach admission uses both original through-flow traces with their respective
bore areas, not only the net pocket flow. Local volume ceilings are bounded
against adjacent column spans, with a separate total ceiling. Cavitation at an
undeclared native liquid node refuses because the finite prescribed pocket
set does not resolve arbitrary distributed multiphase physics.
The fixed liquid geometry, negligible vapor mass and equilibrium interface are
explicit approximations. Trapped gas, capillarity, nucleation delay, moving
meshes, latent heat and strong-collapse shocks are not inferred.

Focused regression sources compare the filled branch to the original exact
DG propagator; recover analytical reservoir/column acceleration and collapse
at both endpoint orientations and an interior station; exercise coupled
end-pocket exchange, simultaneous interior collapse and joint empty-contact
growth; and check actual rupture, native loss, complete energy, both interface
Mach limits, and atomic cancellation/event exhaustion. All six single-pocket
and five multi-pocket tests execute against the same shared implementation.
The simultaneous-collapse fixture also preserves atomic refusal when the
subsequent dispersive collision cavitates at an undeclared liquid station.

### Prescribed flow driving the same liquid/pocket trajectory

`CavitatingDuctPockets::new_with_driven_port` assigns one original endpoint
port to a prescribed flow source, separately from one to eight pocket ports.
`advance_with_port_flow` holds a signed into-duct flow for one interval;
`advance_with_port_flow_and_gas_heat` adds the existing caloric heat callback
on that same trajectory. Flow may change between accepted intervals without
resetting waves, wall memory, gas inventory or the native clock. Existing
advance methods retain zero prescribed drive.

The original generator receives `b_drive Q_drive`. Active vapor/gas pressure
constraints and initial/event contact resolution all include
`B_A^T b_drive Q_drive`; the pure-vapor path caches its affine response.
Two additional original RK coordinates integrate drive pressure impulse and
acoustic work. Complete original-port receipts pass through
`candidate_from_trajectory` before the combined balance subtracts that
external work exactly once. `DuctDriveStep` exposes transported volume,
endpoint absolute pressure, pressure-perturbation impulse and acoustic work.
Reference-pressure transport `p_ref integral Q_drive dt` is reported
separately; it is neither acoustic storage nor generated heat.

Original material, liquid-strain, pocket-volume, interface-Mach and cancellation
admissions remain active. The drive's declared flow also passes its endpoint
Mach and represented volume/work bounds. Cavitation at the drive endpoint
refuses; this boundary does not silently become another pocket. This remains
the original fixed-geometry, small-strain liquid model, with no inferred pump
curve, valve law, moving piston mesh or strong-collapse shock treatment.
The focused driven-port regressions are authored; native execution is not
claimed by this checkpoint.

### Trapped noncondensable gas on the same pocket trajectory

`CavitatingDuctPockets::new_with_trapped_gas` attaches an original
`BubbleGas` to a declared pocket through `pockets::gas::TrappedGasPocket`.
Positive fixed gas inventory contributes p_g(V) to the saturation pressure.
Differentiating p_abs=p_sat+p_g(V) and V'=Q changes the joint constraint to
`(B_A^T B_A + diag(-dp_g/dV)) Q_A = -B_A^T L x`.
The existing liquid generator, ports, adaptive trajectory and simultaneous
vapor-contact solve remain shared. A bounded native LU workspace is allocated
once per interval and reused by gas-bearing RHS solves; the pure-vapor
constant-gain path is unchanged.

Gas-pocket mechanical storage adds -integral(V0..V) p_g dV to the existing
vapor displacement potential. The original duct pressure-work account and
this combined storage close before publication. Positive gas never silently
closes: its authored positive minimum volume is an atomic refusal boundary.
Exactly zero authored inventory reproduces the pure-vapor states and reports.
Nonlinear pressure constraints are checked at accepted endpoints; no liquid
pressure or pocket volume is projected. Original strain, interface-Mach and
local/total pocket-size admissions remain in force.

Four new coupling regressions cover an independently derived distributed
liquid-column/gas-cushion resonance and complete-cycle work, exact zero-gas
reduction, simultaneous gas support and neighboring vapor collapse, and atomic
cancellation/minimum-volume/initial-pressure refusals. All eleven existing
native duct/pocket regressions also pass. These short uniform-pressure pockets
retain fixed liquid geometry and negligible vapor mass. Gas transport,
dissolution, finite thermal exchange and finite-rate phase kinetics require
their own physical owners; the polytropic work is not labelled latent heat.

## Shared material cavitation, films and bubble acoustic history

`cavitation::MaterialCavitatingFilm` retains the original `CavitationProperties`
source and binds its density, viscosity and equilibrium vapor pressure to the
native `fs-tribo` cavitating Reynolds graph. Candidate states, pressure forces,
transport and liquid-mass diagnostics come from that original solver. The
surrounding mechanical owner must admit the pressure force in its own solve;
prescribed-configuration endpoint forces do not establish integrated work or
passivity. Material coefficients remain frozen within each advance. For a
`ResolvedCavitationState`, `advance_with_temperature_source` adopts a refreshed
temperature through the SAME native conservative solve. The new snapshot must
retain the original immutable card, every selected claim and property
requirement, the original temperature axis, and all nonthermal coordinates and
quantity conventions. Dynamic viscosity, density and saturation pressure refresh
together; kinematic-source conversion uses the new density. The original
dilute-vapor cap is checked again.

The returned `MaterialCavitationContinuation` contains the new source binding
and matching candidate inventory for joint acceptance. Its mass balance retains
`rho_old*V_old*theta_old` and uses the new density for candidate storage and face
transport; no existing cavity is reinitialized. Geometry, ambient boundary and
numerical limits remain original. Source, native domain/mass or late cancellation
refusals leave the accepted binding and liquid state untouched. This is
first-order endpoint constitutive splitting with caller-owned temperature,
energy and time refinement. It does not integrate heat, supply vapor/latent-heat
inventory, or certify time-integrated work from endpoint pressure forces.

Focused integration regressions cover source-backed drained heating/cooling
with dynamic and kinematic viscosity, cavity refill at refreshed saturation
pressure, constant-coefficient native identity, original card/claim/context and
vapor-density refusal, and late-cancellation retry. Native film continuation
passes all 136 tests in the actual film-module harness. These four source
continuation regressions now run with the nonlinear source lane below.

### Temperature-dependent nonlinear cavitating films

`MaterialCavitatingFilm::with_power_law` binds explicit reference channel laws
to the SAME native inventory, geometry, pressure supplies and mass controls.
`cavitation::power_law::MaterialCavitatingPowerLawFilm` retains the original
reference material values and caller-authored `K_ref`, `n` and core shear rate.
At an admitted source-temperature endpoint it uses the explicit separable
approximation `K(T) = K_ref*mu_source(T)/mu_source(T_ref)` with fixed `n` and
core rate. Reference viscosity is not substituted for consistency, and a
Newtonian viscosity source does not establish the independent rheological
shape. The caller must admit this separable approximation over its actual
temperature and shear range. Every refresh uses the original reference values,
so repeated thermal cycles cannot accumulate consistency rescaling.

The existing pinned-source continuation check retains the same card, claims,
requirements, temperature axis and all nonthermal query coordinates. Density,
dynamic viscosity and saturation pressure refresh together; kinematic-source
viscosity uses the new density before forming its ratio. The original dilute
vapor cap still applies. The new native conservative continuation carries
`rho_old*V_old*theta_old` into the endpoint pressure/content solve without
reinitializing cavities. Nonlinear moving-wall transport remains the original
combined Couette--Poiseuille profile. Independently prescribed reservoirs retain
their original pressures and contents and refuse if a new rupture pressure
overtakes them. A native Barus law cannot silently become a power-law film.

The returned `MaterialCavitatingPowerLawContinuation` publishes the matching
source binding and native candidate together. All source, constitutive, mass,
domain and cancellation refusals preserve the accepted pair. The endpoint
normal force remains distinct from integrated mechanical work. This first-order
endpoint material split supplies no caloric clock, phase heat, cavitated
tangential heat partition, pressure-dependent rheology or trapped-liquid
compressibility; temperature, energy and temporal refinement stay caller-owned.

Nine focused original-source tests pass: the four existing temperature/source
tests and five new G1/G3/G4 groups. New tests compare the nonlinear thermal
expansion pressure with an independent dimensional n=1/2 closed form; exercise
dynamic and kinematic source conversion, reference-consistency cycling, starved
sliding, retained supplies, exact native Newtonian source parity and atomic
lineage/cancellation refusal. The shared foreign-card fixture now uses the real
`MaterialCard::supersede` API. This scoped run retains the actual material-state
resolver, cavitation material functions, full matdb/thermochem dependencies and
native film solver. Unrelated adapters are excluded; missing `asupersync`
prevents a full workspace run, and no full-constellation result is claimed.

### Source-bound compressive storage and thermal expansion in films

`MaterialCavitatingFilm<ResolvedCavitationState>::with_isothermal_bulk_modulus`
binds the original source to the native `CompressibleCavitatingFilm`. The caller
names an explicit static semantic Pressure query coordinate, which must equal
the native ambient reference pressure. This establishes the meaning of source
`rho(T,p_ref)`; density at saturation or an unspecified pressure cannot silently
become the reference density. An independently authored positive ISOTHERMAL
bulk modulus supplies `rho(T,p)=rho(T,p_ref)*exp((p-p_ref)/B)`. Source sound speed
is not substituted for that modulus. The native geometry, viscosity, absolute
supplies and solve controls remain shared. The compressible owner has a distinct
inventory identity and must supply the initial state.

`MaterialCompressibleCavitatingFilm::advance_with_temperature_source` refreshes
the same pinned card/claims, with every nonthermal coordinate and quantity
convention fixed. The endpoint source supplies density, dynamic viscosity and
saturation pressure together; kinematic viscosity uses the new reference
density. The native solver reconstructs the OLD compressed mass before applying
endpoint coefficients. At constant sealed volume, full liquid therefore obeys
`p_new-p_old = B*ln(rho_ref,old/rho_ref,new)`. Existing cavities and explicit
absolute supplies are retained. The dilute-vapor cap is checked against native
liquid density AT CAVITATION, including its pressure correction, before either
initial binding or an endpoint pair is published.

The returned `MaterialCompressibleCavitationContinuation` contains the matching
source, native owner and candidate state. All refusals preserve the accepted
pair. This first-order endpoint material split retains fixed B and requires
caller-owned temperature, energy and temporal refinement. Reference-card
pressure validity does not certify the independent barotropic law across all
film pressures. No heat capacity, latent heat, compressive-energy receipt or
integrated mechanical work is inferred from endpoint pressure forces.

All thirteen focused original-source film tests pass: the nine retained thermal
and nonlinear groups plus four compressible groups covering already-compressed
sealed heating/cooling, dynamic and kinematic source viscosity, drained
expansion against an independent closed-form finite-B balance, retained supplies,
typed pressure/lineage/vapor-cap refusal, exact native coefficient parity and
late-cancellation retry. Validation uses the original material resolver and
cavitation code, full native fs-tribo, and selected original fs-couple film
modules; it is not a full workspace or distributed-build claim.

When saturated vapor density is supplied, its ratio to liquid density must fit
the explicitly declared dilute-vapor cap; a dense vapor state is refused.
Missing vapor density retains an explicit caller applicability assumption.

The opt-in `cavitation` feature adds `cavitation::bubble::MaterialBubble`. One
source and accepted physical clock drive the selected native Rayleigh--Plesset
or Keller--Miksis owner and
the original compact monopole observer. Density and sound speed must match
exactly. Both the complete integration interval and observer-domain admission
finish before state publication; cancellation or a downstream observer failure
leaves radius, velocity and time unchanged. Acoustic arrival time includes
distance/c. The source bandwidth is explicitly declared and must include emitted
harmonics. `new_compressible` requires the pressure and its analytic derivative
and uses the same accepted-state and observer path. Keller--Miksis owns its
first-order liquid compressibility and acoustic back-reaction; the passive
observer adds no second radiation loss. The Rayleigh--Plesset path remains
incompressible. Source radius, wall-Mach, pressure/bulk-modulus and observer
limits remain explicit and are never relaxed by the coupling wrapper.

`examples/cavitation.rs` uses original IF97 saturation and transport for a
separating/refilling film and a driven Keller--Miksis gas/vapor bubble. Surface tension remains
explicit input. Focused integration regressions cover two-temperature material
substitution, analytical mass/force, deterministic histories, exact medium
binding, the analytic initial compressibility response, exact reuse of the
native trajectory, and atomic cancellation/observer refusal. They do not claim finite
latent-heat feedback, nonspherical collapse, wall jets, erosion, fully coupled
multiphase CFD, or experimental validation.

