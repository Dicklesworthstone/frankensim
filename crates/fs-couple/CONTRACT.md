## Source-bound enthalpy endpoints with prescribed support motion

`SolidEnthalpyState::advance_transient_with_supports_candidate` and
`advance_transient_with_support_accelerations_candidate` pass source-resolved
thermal fields into the native thermal/support Newmark candidates. The existing
specimen, original property-source, reference-inventory, mesh-frame and exact
shared-clock checks apply identically to all three transient endpoint routes.

The returned state retains support displacement, velocity, acceleration,
refreshed thermal stress and reactions together. Complete stiffness stays fixed
on the support path; reference inertia stays frozen. Signed thermal work and
support work remain distinct, with no inferred caloric feedback. Endpoints
still require explicit admission of intervening thermal events and timestep
refinement. The focused `solid_enthalpy` regression advances cold/warm/hot
source states through both endpoint conventions and checks clock/source refusal
without altering accepted thermal or mechanical state.

## Complete heat assemblies on current and voltage motors

Current- and finite-inductance voltage-driven material motors accept
`MotorHeatNetwork<A: HeatAssembly>`, deriving their original material cells
from the exact prepared heat owner. Complete thermal-contact, integrated-layer,
closed-radiation and composed reservoir/finite-link source evaluations
participate in the same native motor/enthalpy endpoint residual. Winding and
bridge heat enter their authored cells exactly once, including an explicitly
homogenized common receiver.

Full associated source samples and typed refusals survive. Graph electrical
ports and prepared contact/emissivity laws with unowned field dependence
refuse, and the existing corrected combined-energy gate admits the closed-radiation numerical net.
The original cached rotor/current owners, material winding claims, source
limits and joint publication retain their existing responsibilities. Original
geometry, fixed mechanical/L/Kt coefficients, source domains, midpoint-current
loss quadrature and first-order thermal limitations remain.

Focused regression definitions cover independent contact/current/voltage
balances, original-entry recovery, and typed late failure/cancellation rollback.
A cross-crate regression uses the actual finite-coolant convection assembly:
changing liquid viscosity reaches winding resistance, source-limited current,
rotor motion and both finite caloric bodies through the same residual.
Native execution remains unverified because the execution workspace became
unavailable before the targeted tests could run.

## Translating supports with omitted-mode contact compliance

Fixed-map normal support translation composes with the original omitted-mode
compliance in `MultiFrictionModalSystem`. The accepted contact closure is
`N^T q - G R - s`; each candidate uses its proposed support offset and the
same original contact potential, secant reaction, retained modal response and
joint friction solve. Translation and compliance may bind in either order.
Changing attachment maps and nonlinear storage retain their existing refusals.

The native account retains prescribed support work `-R delta_s`, omitted-mode
storage `R^T G R / 2`, and the declared numerical loss
`delta_R^T G delta_R / 2` as separate terms. Numerical loss is not material
heat. Held offsets, original material-source admission, released preload and
joint publication keep their existing owners; a refused sample publishes no
motion, contact reaction, offset or clock change.

Focused G0/G1/G4/G5 regression definitions cover independent force/work/storage
relations, both binding orders, opposed supports, held offsets, zero-motion
parity and refusal/retry. Native execution remains unverified because the
execution workspace became unavailable before the targeted tests could run.
The fixed-map quasistatic omitted-mode approximation and time-resolution
limits remain unchanged.

## Complete force histories at moving acoustic receivers

`pcm_wav::observation::moving::force::MovingForceRenderer` now accepts the
original direct and reflected moving source/receiver owners through
`with_moving_receiver` and `with_moving_receiver_reflections`. The existing
accepted C1 volume/force histories, causal source clocks, receiver identity,
whole-window bandwidth admission and `PressureRenderer` output remain the
owners. Callers supply accepted histories and explicit zero force prehistory,
rather than implementing the past-force integral themselves.

For each arrival, every source/image evaluates its complete force past at that
image's CURRENT virtual microphone position. Only the immutable original
force/geometry knot partition is shared between arrivals; a spatial integral
from another microphone event is never reused. Fixed fluid-frame microphones
retain their cached-prefix fast path. Moving queries allocate no buffers and
charge every integrand evaluation, including failed attempts, to the original
finite lifetime work budget. Complete blocks and their receiver clock publish
atomically.

The original continuous source-speed bound resolves narrow history peaks and
checks separation from earlier source positions at the current microphone,
which is stronger than clearance of the current causal pair alone. Unresolved
past clearance, missing geometry, quadrature depth/work exhaustion and
cancellation refuse the complete block. Adaptive Simpson discrepancies remain
error estimates, not interval certificates or evidence of producer interpolation
accuracy. Both force and force derivative must vanish at the first strength
knot; unknown initial acoustic fields are not inferred.

History integrals, pressure and polar-vector velocity retain DD accumulation
through all sources and image transforms before total instantaneous `p*u`
and ideal first-order pickup. Receiver speed changes the sampled source clock,
not pressure gain or fluid perturbation velocity. The original fluid frame,
including admitted tangential uniform flow, is retained; pressure-work flux
excludes mean-flow advection. Focused regressions compare current-position
histories against analytic reflected arrays and independent dense quadrature,
retain weak force fields after volume-source cancellation, and check complete
past clearance, exact block replay and transactional work refusal.

## Tensile admission from original ribbon material histories

Fiber tensile admission composes sourced Kelvin–Voigt stress with each original
Maxwell arm, using actual endpoint velocities, nonlinear axial relaxed-strain
memory and exact linear bending memory. Signed per-layer coefficients combine
before conservative span envelopes. Original source arms, bending projections,
stiffnesses and relaxation times are checked when the native bank is bound.

Cold binding, restore and converged native candidates share this guard. Refusal
publishes no mechanical, history, energy, heat or clock change; prepared buffers
avoid new per-step allocations. Preload remains nonrelaxing. No rupture, stress
clipping or additional material evolution is inferred. Focused regressions
cover axial and bending history overload, retained-memory restore, signed
elastic/viscous/Maxwell cancellation and width partition invariance.

## Exact held-flow window-radiation damping

`WindowRadiationLoad::new_with_sampling` explicitly selects
`WindowRadiationSampling::ExactHeldFlow` for the existing cylinder-window
parallel RLC realization. The original constructor and omitted file selection
retain `BilinearMidpoint`. File-driven voices accept
`mouth_reaction window-cylinder exact-held-flow` or
`mouth_reaction window-cylinder-piston exact-held-flow`.

The continuous fit, static-mass subtraction and radiation geometry are unchanged.
Exact cell histories follow their continuous dissipative poles under the solved
held mouth flow. The implicit mouth equation uses continuous interval-mean
pressure, and the radiation ledger uses independently integrated resistor heat.
The reported sampled resistance and discrete band error use that same mean
pressure transfer; they do not reuse the bilinear frequency warp. Preparation
is shared immutably by candidate histories, with no per-sample allocation or
transition rebuilding. This forcing declaration resolves no variations inside
the held interval and adds no claim about unmodeled exterior geometry.

## Nonlinear electrostatic volume forces on native solids

`electrostatic::spatial_nonlinear_coulomb_volume_load` transfers the accepted
native tetrahedral nonlinear electric field into the existing solid P1 load
owner. It shares the original linear adapter's exact mesh/frame/connectivity
admission and integrates only the selected original free-volume-charge term
rho E. Each nodal contribution is Volume*rho*E/4; original cell order, force
resultant, first moment, affine-displacement work and fixed-DOF reactions are
retained without remeshing or a separate electric solve.

The source field retains its full material response and all sources of E.
Only free volume charge assigned mechanically to the selected material enters
this adapter. Electrode/interface/polarization forces and electrostriction are
separate contributions; an overlapping enclosing stress load must not count
the same charge again. This is an instantaneous load snapshot, with no implicit
charge transport or field/mechanical time integration. Cancellation and mesh
refusal leave the field reusable. Synthetic nonlinear charged-tetrahedron
regressions specify force, torque, affine work and transactional refusals.

## Finite-temperature feedback on the three-dimensional induction rotor

`SpatialMagneticThermalNetwork::step_rotor_candidate` advances the original
spatial field, native finite-inertia rotor and finite enthalpy together. The
linear material/entropy residual is shared with standalone induction and
physical electrode circuits. Every temperature proposal resolves original
conductivity and magnetic sources, then restarts the same accepted field and
mechanical histories for the bounded speed closure. Source-weighted operators
are prepared once per temperature proposal and shared across its field trials.

The joint driver includes magnetic free-energy change minus original parameter
work, native kinetic-energy change, impressed supply and external shaft work.
The original magnetic-caloric owner deposits reversible field exchange and
recovers magnetic internal energy. Actual conductor-frame Joule heat enters
the original finite body once; Lorentz work cancels internally against rotor
motion. Final admission independently closes magnetic internal energy, kinetic
energy and finite caloric storage against external electrical/shaft work and
external heat. Bearing reactions retain the native full impulse account.

All inner speed/field trials consume the original thermal work budget. Invalid
caller controls are admitted before remaining-budget restriction. Cancellation,
source refusal and exhausted solves publish no field, angle, speed or thermal
state. Complete heat-assembly entry points retain contact/layer/radiation and
finite-coolant source payloads and typed errors; closed-radiation roundoff is
not external supply, and unowned electrical contact/resistor work refuses.

The rotor remains a homogeneous isotropic cylinder/tube with uniform scalar
permeability throughout the field domain. A whole-body finite caloric cell
provides its uniform temperature; a nonuniform rotating material source refuses.
Endpoint thermal coefficients make the combined method generally first order.
Moving interfaces, transported magnetic memory, phase transitions and thermal
expansion remain outside this model. One focused source-bound regression
definition covers heating-dependent braking through two retained intervals,
independent current/field/kinetic/enthalpy balance and trial-budget rollback.
Native execution remains unverified while constellation dependencies are absent.

## Native three-dimensional induction with finite rotor inertia

`spatial_induction::SpatialInductionRotor` closes the original Whitney A–phi
field against one native centered, axisymmetric revolute-z body. A bounded
implicit iteration uses the same accepted field and momentum on every trial;
its converged speed supplies `v = omega cross (x - axis)` and its original
current supplies both conductor-frame Joule heat and the full Lorentz wrench.
Three-dimensional charge continuity and end effects remain in the native field.

The `induction_rotor` stage is shared with the original planar rotor. It calls
the original articulated forward kinematics and constrained inverse-inertia
map, retaining shaft input, full body force/torque and ideal-bearing reaction.
Actual load work, kinetic work, linear/angular impulse and combined field/rotor
energy must pass before the joint candidate publishes. Coupling exhaustion or
cancellation leaves the accepted cochains, source, angle, speed and time intact.
Geometry bindings and the mechanical model are shared across accepted steps;
a single velocity buffer serves every speed trial of an interval.

Admission requires homogeneous scalar conductivity, uniform scalar induced-
linear permeability throughout the field domain, and no impressed current in
the rotor. The original fixed tetrahedra approximate a complete circular
cylinder or coaxial tube whose analytic occupancy is rotation-invariant.
Constructor checks include connected oriented cells, closed manifold boundary,
full facet/cell radial bounds, end coverage and bounded chord error. They do
not certify arbitrary intersections or field convergence. Noncircular moving
interfaces, anisotropic rotating material, hysteretic-memory transport and
magnetization-interface loads require additional owners.

Three focused G0/G3/G4 regression definitions cover transverse braking in both
directions, independently integrated P1 force/torque/Joule and magnetic energy,
shaft drive followed by coasting, bearing reactions and transactional refusal.
The fixture uses real material-resolved spatial fields and native inertia;
synthetic coefficients do not establish experimental machine validation.
Pinned formatting/syntax and whitespace checks pass. Native compilation and
regression execution remain unverified while the seven constellation sibling
repositories are absent.

## Reciprocal heat and finite porous compaction

`solid_enthalpy::finite_plastic::porous::FiniteThermalPorousSolid` couples retained non-isochoric plastic history to the original `DepositionSimulation` caloric and transport owner. A native finite porous tetrahedral equilibrium is solved at each proposed endpoint temperature; original yield heat plus signed reversible caloric exchange enters the thermal timeline exactly once. Every mechanical and thermal proposal starts from the same accepted endpoint. The complete thermal, structural and observer candidate is accepted together, after the existing temperature, per-cell heat, thermal-energy and combined-energy gates. The original bounded Picard/IQN-ILS driver is shared with dense finite plasticity and sheet forming.

`FiniteThermalPorousMaterialBinding` resolves the complete source-bound thermal porous family and preserves the exact original caloric chart, typed context, pinned Cp/density receipts and same material-card identity. Common caloric admission is shared with the original J2 owner: an explicitly selected fully solid, stress-free zero-pressure reference chart, its complete temperature interval, and a positive conservative lower bound on caloric plus mechanical heat capacity. The bound is a constitutive admission check, not proof that the nonlinear coupled return has a unique solution.

The original reference tetrahedral volume and source density fix thermal inventory and optional consistent inertial mass. Actual density uses `rho0/det(F)`, with retained plastic volume available separately. Caloric density metadata at hot states does not replace original mass. Shared mechanical reports retain physical yield heat, signed reversible exchange, internal/free-energy changes, coefficient storage and signed endpoint-work defects. Kinetic change and backward-Euler loss are accounted independently; neither numerical defect nor kinetic loss is deposited as heat.

Thermal cells and conductances remain a supplied reference network, retaining the original observer sources, transport, events and histories. No moving-area Fourier operator, pore-fluid pressure, phase flash, damage, rate dependence or calibration is inferred. The native dry associated law, affine thermal envelope, P1 locking and physical interval/refinement limits remain explicit.

The four `finite_porous_enthalpy` regressions specify a loaded adiabatic crush against an independent hydrostatic/Voce/caloric oracle, two-cell physical heat redistribution, retained hot inertial history through late refusal/retry, and original source/chart/capacity/inventory admission. Native execution remains a separate verification requirement.

## Native electrode circuits and coupled finite-temperature feedback

`spatial_induction::terminals::circuit::SpatialTerminalCircuitThermalNetwork`
advances native circuits driving nonlinear or hysteretic solid electrodes with
the original finite-enthalpy state. The nonlinear and hysteretic child modules
reuse their existing field thermal residuals; no additional thermal integrator
or approximate conductor law is introduced. The joint candidate retains the
native circuit, complete field/magnetic memory and full thermal source receipt.

Each endpoint trial resolves bulk scalar/tensor conductivity from its original
material. Optional complete graph-resistor bindings resolve every authored R(T)
from a separate whole-cell finite body. Accepted source identity, electrical
basis, geometry, caloric temperature and original resistance must agree, and
trial refresh cannot replace them. Numeric authored resistors retain explicit
heat destinations. Graph and field histories restart together on every trial.

Actual local bulk Joule and irreversible pinning heat enter the original field
receivers; authored graph-resistor heat enters its explicit receivers once.
Canceled numerical reference heat and numerical magnetic dissipation never
enter caloric storage. Stable graph/magnetic storage, finite thermal storage
and numerical loss close against original external/impressed supply and external
heat. The nonlinear path preserves its optional source-separable magnetic
entropy and parameter-work account. Hysteretic fits remain constant throughout
the admitted caloric interval, while their vector memories remain retained.

Complete heat-assembly adapters preserve original contact, layer, radiation and
finite-coolant payloads, typed failures and closed-radiation accounting.
Unowned electrical contact work refuses. Nested circuit/field trials consume
the shared thermal budget, and cancellation or any failed source, solve or
balance leaves all accepted histories unchanged. These adapters describe
stationary, fixed-geometry, fully solid bodies within the existing MQS scope.

## Actual-temperature feedback into finite-strain membrane relaxation

`membrane_thermal::visco::temperature::advance_candidate` feeds each original
caloric receiver's actual temperature into native Maxwell relaxation times.
`ViscoMembraneKinetics` retains the original unshifted spectra and one explicit
WLF/Arrhenius calibration per triangle. Every preparation uses the shared
`TemperatureShift::shift_hyperelastic` owner; repeated intervals never rescale
an already shifted spectrum. The native owned-motion kinetic gate pins all
energy parameters, natural references, branch ordering and solver controls.

The original thermal timeline chooses the interval. Its accepted starting
temperatures set the kinetics; native physical branch loss returns to the same
enthalpy cells. Accepted steps and thermal events synchronize endpoint times
without changing history, physical storage, momentum or interval heat. The
frozen Kelvin--Voigt, frozen Maxwell and kinetic-feedback routes share receiver
admission, source/inventory checks, physical/numerical work accounting and atomic
candidate publication. Numerical kinetic damping never heats a cell.

Both kinetic calibration support and the explicitly admitted fully-solid
receiver interval are checked, including after final heating. Cancellation,
foreign energy/reference laws, out-of-domain endpoints and late observer
refusals preserve all accepted owners. Temperature-independent free energies
make this kinetic feedback compatible with the existing caloric account; no
reversible thermal work is inferred. This is first-order lagged-temperature
coupling, with explicit timestep refinement still required. Thermal expansion,
temperature-dependent elastic moduli and reciprocal entropy exchange remain
separate owners.

`tests/visco_membrane_thermal/temperature.rs` covers native shifted/reduced-time
agreement, self-heating at three timestep rungs, heterogeneous WLF/Arrhenius
regions, event/continuation equality and late whole-call refusals.

## Retained kinetic storage through full and partial rigid liquid release

Rigid liquid release compares original native momentum/inertia storage with
the independently integrated P1 field at a fixed COM velocity reference,
retaining mass and momentum defects. Its final gate adds the native transfer's
retained P1-to-cell increment. Selected-cell release independently checks the
actual native rigid partition, then composes that increment with the complete
selected-piece release. Successive releases can consume retained increments
without subtracting rounded totals. Original inertia, mass/momentum, caloric
state and independent remap budgets remain unchanged; numerical loss never
becomes heat or release work. Four executed production-leaf tests cover hidden
spin/orbital storage, actual mass/momentum defects and invalid inertia. Native
whole/partial regression assertions additionally cover analytic P1 spin,
thermal attachment and successive releases; full native execution is pending.
The source remains an already-liquid no-storage rigid constraint, not a
stressed solid whose constitutive history can be discarded.

## Native solid-electrode circuits with temperature feedback and finite heat

`SpatialMagneticThermalNetwork::step_terminal_circuit_candidate` closes the
original 3D Whitney A–phi electrodes, native RLC graph and finite caloric cells
in one accepted candidate. Every thermal proposal resolves scalar or oriented
tensor conductivity and induced-linear permeability from the original material
sources. Field-response probes retain the same accepted magnetic potential,
native charge/flux and thermal history; their actual count consumes the shared
thermal trial budget. Endpoint operators are prepared once per proposal and
the native field response is invalidated when material coefficients change.

Original local conductor-frame Joule heat enters whole-volume field receivers.
Every authored graph resistor additionally requires one explicit finite heat
destination, which may be a separate thermal body. These physical losses are
deposited once; artificial reference-resistor dissipation is excluded. Native
capacitor/inductor storage and external source work retain their own accounts,
while electrode interface power cancels internally. The existing magnetic
caloric owner preserves signed parameter work, field entropy and reversible
heat instead of treating every magnetic energy change as Joule loss.

The joint receipt checks `delta U_field + delta W_circuit + delta H_caloric
+ W_mechanical = W_impressed + W_external + Q_external`. Original sources,
clocks, volume correspondences, thermal laws, material and solve budgets remain
admitted before joint publication. Refused heat routing, source changes,
cancellation and failed field or enthalpy trials preserve all three histories.

The legacy fixed-resistance entry retains authored R/L/C coefficients,
geometry and conducting occupancy; prescribed motion represents external mechanical work. This is a linear
magnetic, fully solid MQS coupling with endpoint material/enthalpy feedback,
generally first order in the coupled temperature update. Independent Whitney
element elimination tests exercise both electrical property bases, signed
permeability slopes, separate resistor heating and retained capacitor charge.

## Complete heat transport and material resistor feedback for solid electrodes

`spatial_induction::circuit::heat_assembly_candidate` continues the same native
linear A–phi electrode circuit through any original `HeatAssembly`. Every
proposed caloric state resolves contact resistance, integrated conduction
layers, gray-diffuse enclosure radiation and provider-owned finite thermal
links or coolant together with conductor conductivity and magnetic entropy.
The original backward-Euler residual advances every finite receiver once;
there is no independent cooling step or extra temperature history. The full
associated source sample and provider error type remain in the accepted
candidate or typed refusal.

`SpatialTerminalCircuitThermalNetwork::step_linear_material_candidate`
also binds every authored graph resistor to its original electrical material,
basis and whole-volume finite body. The complete heat adapter consumes those
same bindings from its prepared electrical/material network. Accepted graph
resistance and retained source state must agree before stepping. Every trial
pins original claim selection, full fixed conditions, geometry and receiver,
then rebinds the original R(T) without resetting native charge, flux, time or
magnetic potential. A changed resistance invalidates the native response cache.
Source-backed graph resistors occupy distinct caloric bodies outside the field
receivers; their length times area must equal the corresponding original volume.
The earlier fixed-resistance entry retains its original API and numerical
resistances; it refuses to discard an already retained material binding.

Both entry points share the same linear field/circuit driver, graph-source
admission and authored-loss deposit owner used by the original coupled paths.
Conductor-frame Joule loss and actual graph resistor loss enter their explicit
finite receivers once. Artificial interface-reference dissipation is canceled
inside the native circuit and never becomes physical heat. Magnetic parameter
work, entropy and reversible heat retain their original separate accounts.
Finite-link exchange and closed-enclosure radiation are internal; both the
native and physical energy checks exclude the enclosure's numerical net from
external heat. The joint account remains
`delta U_field + delta W_circuit + delta H_caloric + W_mechanical
= W_impressed + W_external + Q_external`.

Fixed geometry, fully solid field receivers, linear magnetic source domains,
reactive coefficients and first-order thermal feedback remain the admitted
scope. Prescribed motion exchanges external work without advancing mechanics.
Thermal boundaries requiring an unowned E/J/B/H coordinate and electrical
interface contacts without their native work owner refuse. Cancellation,
changed source identity or geometry, invalid routing and exhausted budgets
return no partially advanced field, circuit or thermal history.

Focused regressions use independent Whitney-element elimination and two-body
caloric balances for both electrical source bases, continued drive reversal,
R(T), sigma(T), mu(T) and magnetic entropy. Complete-assembly cases exercise
source-controlled contact and emissivity changes, integrated layers, accepted
provider payloads, refreshed finite links, typed source refusal, cancellation
and unowned field-conditioned radiation.

## Physical nonlinear and hysteretic electrodes with finite heat

## Finite-strain Maxwell membrane relaxation into original enthalpy cells

`membrane_thermal::visco::advance_candidate` composes the native
`FiniteViscoMembraneMotion` with the existing finite enthalpy timeline. Each
actual transport, observer or event interval advances the original momentum
equation and Maxwell histories once. Physical per-triangle interval loss is
handed to the matching original caloric receiver before cumulative-ledger
rounding; variable heat capacity, original graded reference mass and volume,
source identity, fully solid phase and the declared temperature interval retain
the existing membrane thermal admission.

The Kelvin--Voigt and Maxwell entry points share the same private receiver,
timeline and heat-handoff implementation; each mechanical owner retains its
native physics. Backward-Euler kinetic loss is reported separately from heat.
Material endpoint defects, pressure path/endpoint work differences and inertial
or equilibrium residuals contribute to an absolute numerical allowance without
cancelling one another. Held supports must have their accepted displacement and
zero velocity throughout the coupled call. Dead loads and admitted current-area
pressure remain held.

The result publishes motion, thermal inventory and observer together. A late
observer refusal, source/inventory mismatch, work refusal or cancellation
preserves all accepted inputs. Mechanical coefficients remain frozen over the
explicitly admitted fully solid interval: this is one-way self-heating, not
temperature-driven coefficient or natural-reference evolution.
`tests/visco_membrane_thermal.rs` covers regional Maxwell heat into variable-Cp
enthalpy, exclusion of numerical kinetic loss, event/continuation history and
whole-call refusal behavior.

`SpatialMagneticThermalNetwork` supplies `step_nonlinear_voltage_candidate` and
`step_nonlinear_current_candidate`. With `magnetic-hysteresis`, the corresponding
`step_hysteretic_voltage_candidate` and `step_hysteretic_current_candidate`
advance actual solid-conductor electrodes, the original 3D A–phi field and
finite caloric inventory together. The shared `SpatialTerminalThermalCandidate`
retains the typed native electrode owner, finite thermal state and full coupled
receipt. Candidate failure or cancellation publishes neither history.

Private read-only history traits let the existing finite-enthalpy drivers read
accepted electrode geometry and sources without exposing isolated field
stepping. Every thermal trial starts from the same accepted magnetic history,
retains electrode charge topology, and resolves scalar/tensor conductivity at
the actual endpoint temperature. Nonlinear electrodes preserve the existing
optional source-separable magnetic entropy law and effective-capacity checks.
Hysteretic electrodes retain vector polarization and require fixed fitted
coefficients over the full admitted caloric interval.

Actual terminal V*I work joins impressed supply once. Native local Joule heat,
and irreversible pinning heat where present, enter their original whole-volume
solid receivers once; numerical magnetic dissipation remains outside deposited
heat. The joint work balance closes magnetic/caloric storage and numerical loss
against electrode work, impressed work and external heat. The nonlinear path
also retains the original reversible magnetic entropy/parameter-work account.
Thermal clocks, complete magnetic source identity, receiver volume, field and
finite-enthalpy budgets all retain their original admission.

The `nonlinear::terminals` and `hysteretic::terminals` modules additionally expose
`voltage_heat_assembly_candidate` and `current_heat_assembly_candidate` for the
complete original contact, layer, radiation and finite-coolant heat provider.
Provider source payloads, typed errors and closed-radiation accounting are
preserved. Electrical graph/contact work without an electrical owner refuses.
These thermal adapters describe stationary, fixed-geometry, fully solid bodies;
material transport, phase transitions and transported magnetic memory remain
outside this scope.

## Native 3D hysteretic winding circuits and finite enthalpy

With `magnetic-hysteresis`, the existing `SpatialWindingThermalNetwork` supplies
`step_hysteretic_candidate` and `step_hysteretic_material_candidate`. The native
3D graph/winding, A–phi/vector-memory field and finite caloric inventory advance
in one candidate. The material form additionally resolves every ordinary graph
resistor from its separate finite body. The complete
`spatial_induction::hysteretic::winding::heat_assembly_candidate` retains original
contact, layer, radiation and coolant source payloads and closed-heat accounting.

The original field thermal residual and the existing winding source/heat routes
serve both excitation forms. Each thermal endpoint resolves actual strand
resistivity, core scalar/tensor conductivity and sourced graph R(T). Magnetic fit
coefficients remain constant over the complete admitted caloric interval.
Neither thermal nor terminal trials reset or partially publish accepted circuit,
field, polarization or caloric history. The shared thermal budget counts nested
field trials, including discarded terminal proposals.

Core eddy, all-tetrahedron pinning, strand copper and graph-resistor heat are
deposited once into their original finite receivers. Insulating dissipative
cores require receivers too. Strand fills must cover their whole original
caloric volume; partial composites require a separate mixture owner. Ordinary
material graph bodies must be distinct from those whole magnetic volumes.
Numerical field/end loss remains outside caloric heating. The stable account is
`delta W_graph+field+end + delta H_caloric + D_numerical = W_external + Q_external`.

All source, original-volume, clock, terminal, finite-enthalpy and total-work gates
precede joint publication. Source changes, cancellation and exhausted budgets
leave every accepted history intact. These adapters remain stationary and fully
solid; they do not infer transported/rotated magnetic memory or a
temperature-dependent hysteresis free energy from a fixed fitted law.

## Three-dimensional hysteretic induction and finite heat

With `magnetic-hysteresis` (also selected by `hysteretic-eddy-heating`),
`SpatialMagneticThermalNetwork::step_hysteretic_candidate` advances the retained
three-dimensional Whitney A–phi field, all vector branch polarizations and
the original finite caloric inventory in one immutable candidate.
`spatial_induction::hysteretic_heat_assembly_candidate` supplies the same solve
with the complete original contact, layer, radiation and finite-coolant heat
assembly. No independent field or thermal history is restarted.

Every endpoint trial resolves scalar or oriented tensor conductivity from its
original material source at the actual proposed temperature. Planar and spatial
hysteretic consumers share the same fitted-source and full caloric-interval
admission. Magnetic branch coefficients remain constant over that supported
interval; changing temperature cannot silently replace the fit or erase its
memory. Complete magnetic source receipts retain the actual endpoint query.

Exact native Joule and irreversible pinning heat enter their assigned finite
cells once. Electrically insulating magnetic tetrahedra also need a real
thermal destination whenever they dissipate. Original whole-volume, fully
solid material correspondence, source identity and clocks are admitted before
thermal trials. A fixed region with unreceived pinning heat refuses.
Numerical endpoint dissipation is a separate energy export and never a heat
source. The stable work account closes
`delta W_magnetic + delta E_caloric + D_numerical = W_impressed + Q_external`.
Both physical heat channels remain individually available in the native step;
`deposited_heat_j()` reports their already-deposited total.

Complete heat providers retain their original typed errors and closed-radiation
balance, and refuse graph electrical work without an electrical owner. Native
field and finite-enthalpy budgets, source failures and cancellation publish
neither accepted history. This is stationary, fixed-geometry, first-order MQS
with source-admitted isotropic vector hysteresis and constant magnetic fits.
Transported/rotated magnetic memory, temperature-dependent fit entropy, phase
changes and empirical device qualification require their respective extensions.

## Fatigue on original finite-rotation elastic rod motion

With `rod-dynamics`, `rod_motion::fatigue::RodMotionFatiguePath` owns one native
rod trajectory and the shared fixed-material-point fatigue observer. Both
centered and eccentric reference inertia retain their original native motion
step, support reactions and energy/momentum accounts. World-frame nodal loads
and prepared material-point spatial loads use their existing kick/drift/kick
owners. The fatigue observer receives the accepted pose, so unforced vibration
still contributes actual signed stress cycles.

Admission observes the supplied mechanical state once, including its preload
and momentum, without reconstructing earlier damage or inserting a zero-stress
sample. The rod observer retains original per-layer axial/bending/torsional
stress, material-frame plane directions, typed source conditions, bounded S-N
rainflow and optional explicitly closed Findley windows. Its explicit
`new_with_transverse_shear` constructor includes the original material owner's
Saint-Venant transverse field for homogeneous isotropic circles and tubes.
It retains both free contours, source roles and whole-layer shear bounds.
Other constructors retain the solid observer's omitted-shear allowance.
No constructor treats a section-averaged force as a local stress or changes
the native shear correction to fit the observation.

A native candidate is accepted only after every gauge accepts its tensor.
Mechanical refusal, late gauge failure, sample/storage budget exhaustion and
cancellation leave pose, momenta, clock, accepted-step count and all histories
unchanged. Ordinary observation clones bounded rainflow residue and rolls back
appended Findley prefixes, without copying the retained tensor window each
step. Clones preserve continuation; provisional endpoint damage does not seal
the live path. Fatigue introduces no mechanical force or heat.

This consumer exposes immutable motion and uses the native conservative elastic
step at fixed source conditions. Viscous constitutive stress, source-temperature
refresh, impulses, contact and fracture are not silently omitted from a coupled
evolution; those operations require their own compatible fatigue continuation
and are not exposed by this path. Caller-selected time steps must resolve the
actual stress extrema and any rotating shear trajectory. A finite plane set
and empirical Miner or Findley estimate do not establish component service life.

## Retained kinetic changes through layer heat and cohesive motion

The rigid third-body shear-layer owner evaluates accepted storage through
`fs_mbd::RigidKineticEnergyChange`, reducing both original bodies before final
rounding. Small layer elasticity and irreversible work survive large common
translation. The finite layer's existing thermal feedback continues to consume
its independent constitutive loss and original caloric inventory; storage does
not define heat. The source-bound scalar, thermoelastic and tensor routes share
this same native mechanical owner.

Deformable cohesive heat, finite adhesive layers, prescribed source-temperature
transfer, opening/cooling, accepted intervals and sustained/event trajectories
consume the native `kinetic_energy_change_j`. They do not reconstruct it from
absolute endpoint totals. Existing bulk/cohesive storage, fracture work, explicit
heat-conversion fractions, source work and caloric loss accounts remain separate.
The original fixed-mass, fixed-reference, source and topology limitations apply.
G3 regression definitions carry actual rigid-layer heat/temperature feedback and
deformable bond-viscosity heat through large unchanged kinetic inventories.
Arithmetic-leaf tests alone do not qualify the complete coupled integrations.

## Native kinetic increments through spatial heat and thermal shells

The spatial viscoelastic thermal owner consumes the accepted native
`ViscoDynamicStep::kinetic_energy_change_j` and
`material_stored_energy_change_j`. These come from the original full or
supported consistent mass and original per-element constitutive increments.
Explicit support and eigenstrain work remain separate, and only independently
integrated Kelvin/Maxwell loss enters the original finite caloric cells.
Unchanged background translation cannot hide a small physical heat transfer
or enlarge the native interval energy allowance.

The native step's read-only `retained_energy_change_j` and `retained_work_j`
carry the combined storage and independently supplied work in double-double
form through the thermal consumer's complete balance. Opposing large kinetic
and material increments are reduced before diagnostic rounding; wrapping their
already-rounded public scalar fields in compensated arithmetic is insufficient.

The reciprocal dynamic laminate owner consumes the original enriched MITC
midpoint's retained kinetic increment directly. Its original director mass,
nodal/director cross terms, material stress work, finite heat cells and launch
coordinate basis remain shared across temperature trials. Reconstructed
absolute kinetic totals cannot reject or erase the original small exchange.
These remain fixed-reference, fixed-mass models with the existing material and
time-resolution obligations. Focused regression definitions cover supported
spatial strain heat and independent director motion beside dominant unchanged
kinetic storage; isolated arithmetic checks do not execute those integrations.

Finite-solid time steps, impulses, force kicks and kinetic drift retain actual
kinetic increments through the unchanged original consistent mass. Endpoint
difference and mean-velocity tails remain compensated through independent
storage and supplied-work checks, including pressure and gravity ports.
Background kinetic storage cannot erase a resolved local exchange. Numerical
integration, impulse-work and force-splitting defects retain their original
gates and never supply physical heat. `FiniteSolidStep::kinetic_change_j`
exposes the same increment used by the original acceptance calculation.

Dynamic source-bound tetrahedral thermoelasticity consumes the native Newmark
owner's retained kinetic increment, independent force/support work and explicit
Rayleigh loss. It retains those quantities through stress-work and combined
caloric acceptance instead of subtracting absolute kinetic diagnostics. Rayleigh
energy stays explicitly exported; local reversible constitutive heat alone
feeds its original material cells. The underlying consistent mass, stationary
support restriction, source ranges and time-refinement obligations are unchanged.
Its finite-gas cavity composition consumes the same retained mechanical storage
and work rather than rounding them again before the combined gas/wall balance.

The original rigid caloric partial-slip owner and its shared-body network now
retain native fixed-inertia kinetic increments across all participating bodies
and through their complete material/thermal balances. They recover relative
midpoint motion from original COM/spin/arm components before common translation
can swallow real slip. Contact free energy, independently computed constitutive
work, per-surface reversible entropy heat and physical microslip heat keep their
original owners; return-map numerical loss never enters a thermal cell. The
existing absolute energy tolerances do not grow with unchanged kinetic storage.
Focused regression definitions compare actual shear/spin impulses, histories,
physical losses and finite caloric endpoints with and without a dominant common
normal translation. Full coupled execution remains separate from arithmetic
leaf checks.

## Original specific heat in coupled capacity admission

Electrical, dielectric, magnetic, piezoelectric, thermoelastic, creep and
tensor/friction-layer thermal admission now use the original caloric interval's
minimum Cp, rather than its average enthalpy/temperature chord. Clipped
temperature envelopes use conservative intersecting-interval source bounds.
Creep retains its outward arithmetic bounds. Existing source identity,
reference mass/density, field-work and mechanical-energy accounts remain
authoritative; no second caloric or apparent field-capacity term is inserted.

String and concentric-layer thermoelastic admission evaluate the quadratic
TOTAL capacity, including the original affine Cp slope, at endpoints and its
interior stationary point. Local finite-solid and small-strain adiabatic wave
gates compare original Cp at the actual state. Continuous Cp remains
differentiable when its temperature slope changes; an actual capacity jump
still refuses a classical wave linearization.

Native regressions cover continuous source-Cp wave states and a two-knot
capacitor/inductor caloric source whose average capacity hides negative total
capacity elsewhere in its admitted envelope. Refusal leaves the original
thermal and circuit state unchanged. These whole-consumer tests were not
executed in this change because the real asupersync dependency build exceeded
the available memory; numerical-leaf results do not substitute for them.

## Source-dependent finite-solid heat capacity and adiabatic waves

`FiniteSolidMaterial::adiabatic_acoustic_state_with_coefficients` and the
corresponding `ConservativeThermalField` element view retain the original
Neo-Hookean E(T)/nu(T) source derivative owner alongside the existing caloric
chart, density, expansion path and manufacturing reference. The selected
model, complete fixed context, card and actual temperature must match. These
are the same prepared mechanical states used by the ordinary fixed-coefficient
path; the added source capability supplies their thermal derivatives.

The shared native thermoelastic potential combines coefficient curvature,
expansion derivatives and their mixed terms. Its `-T Psi_TT` enters total
fixed-F capacity once, together with the original source capacity at the
actual state and reference density. The full beta=dP/dT and positive total
capacity then feed the existing compressible adiabatic acoustic tensor,
signed modes and energy velocities. Original source receipts remain borrowed
by the result. Nonpositive capacity, foreign evidence, unsupported smooth
neighborhoods and unresolved states refuse.

The original stress-free zero-pressure caloric convention remains explicit.
This is a local reversible full-compressible response; constrained mixed or
incompressible projections, heat transport and attenuation require their
existing separate owners. Finite temperature evolution at held deformation is
available through the caloric transaction below.

## Finite heating and cooling at held deformation

`FiniteSolidMaterial::fixed_deformation_caloric_state` retains one immutable
original material, coefficient source, fully solid caloric curve, total F,
manufacturing reference map, reference mass and reference volume.
`FiniteSolidCaloricState` explicitly uses internal-energy storage with the
existing `CaloricReference::StressFreeZeroPressure` convention:
`U = m*h0(T) + V0*(Psi - T*Psi_T)` and
`dU/dT = m*Cp(T) - V0*T*Psi_TT`. Mechanical caloric coupling enters once.
Each temperature query evaluates the original E(T), nu(T) and natural stretch,
then returns the evolved native stress, tangent and adiabatic response.

`add_heat` accepts signed finite heat and reuses the shared safeguarded caloric
inverse. Original Cp and Cp/T finite-pair integrals cross retained source knots;
direct native mechanical increments avoid cancellation against large absolute
caloric or elastic baselines. A bounded sampled temperature bracket requires
positive combined capacity and increasing energy at every accepted sample.
The final result independently checks heat closure and endpoint Clausius
bounds. Explicit spacing, state-evaluation, root and entropy controls bound
work. Refusal or cancellation leaves the original state unchanged; zero heat
returns that state without claiming an unsampled branch.

This is a fixed-F, fixed-reference-volume transaction on the original fully
solid smooth source branch. Expansion must remain affine across the prepared
path; unsupported kinks, foreign sources, unbracketed heat and nonpositive
sampled capacity refuse. Numerical sampling does not establish global
monotonicity, stability between samples or a thermal CFL bound. Free mechanical
motion, phase change, finite-pressure Cp conversion and spatial heat transport
remain outside this transaction.

Six original consumer tests pass in the native selected finite-solid closure,
using rebuilt original material and thermochemical code with real dependencies.
They cover heating/cooling across a Cp knot, refreshed stress and waves, a
sub-ULP caloric datum, weak mechanical energy increments, constant-source
reduction and atomic source/budget/cancellation refusals. Unrelated FE/time
adapters and thermal-field wrappers were not included in that validation;
this is not a full workspace build.

### Timed reservoir exchange at held deformation

`FiniteSolidCaloricState::exchange_with_reservoir` advances that same original
inventory for an explicit duration against an authored
`fs_thermal::enthalpy::Reservoir`. The descriptor retains its covered area,
prescribed bath temperature and total conductance G in W/K; area does not
multiply G again. One backward-Euler inverse solves
`Delta U(T) + dt*G*(T - T_bath) = 0` with the actual tangent `c_F(T) + dt*G`.
It reuses the finite-heat source cache, bounded branch sampling and direct
energy/entropy increments. Temperature-dependent Cp, elastic coefficients and
expansion continue to affect the same evolved stress and acoustic response.

`FiniteThermoelasticReservoirBalance` retains the actual endpoint power and
`Q = dt*G*(T_bath - T_final)` separately from the independently evaluated
internal-energy change. It exposes the body's original heat balance, endpoint
conduction entropy `Q*(1/T_final - 1/T_bath)`, temporal entropy defect
`Delta S_body - Q/T_final` and combined body-plus-bath entropy change. The
time-discretization contribution is not relabeled as physical generated heat
or boundary dissipation. The existing energy/entropy allowances govern final
acceptance; late cancellation leaves the supplied state unchanged.

The prescribed reservoir has no represented finite heat capacity. Its
temperature may lie outside the material's source interval, while every body
state and the offered branch remain supported. Conductance and bath
temperature are held over the explicit interval; duration refinement remains
independent of nonlinear solver tolerances. Insulation, zero duration and an
equal-temperature bath retain the original state after input admission. This
transaction does not add spatial contact topology or infer a global thermal
slope bound from sampled material capacity.

Three additional original consumer tests cover the constant-Cp backward-Euler
closed form, refinement toward exponential cooling/heating with decreasing
temporal entropy defect, variable-source heating/cooling across a Cp knot,
updated mechanics, exact large-datum invariance, supported body states with an
external bath temperature, no-op behavior and atomic domain/budget/cancellation
refusal. The selected native consumer root passes all nine tests using the
same real dependency and original-source closure described above.

## Circuit-driven spatial hysteresis and finite heat

The `hysteretic-eddy-heating` feature is an alias of `magnetic-hysteresis`;
both select one shared field/heat owner, including complete heat assemblies.
Original fitted magnetic claims are retained at the actual endpoint temperature.
Static applicability and the full caloric source interval are admitted once per
candidate before thermal trials.

With the default-off `magnetic-hysteresis` feature,
`WindingThermalNetwork::step_hysteretic_candidate` and
`step_hysteretic_material_candidate` advance the original native circuit,
vector magnetic history and finite caloric state in one immutable candidate.
The material form also resolves original graph resistors at their receiving
cells' endpoint temperatures. Stranded copper and core conductivity are
resolved independently from their original material cards without replacing
the fitted magnetic sources or resetting remanence.

The existing terminal and enthalpy solvers retain ownership. Every discarded
thermal or terminal trial restarts the same accepted circuit, field and branch
polarizations; nested field work respects the original shared trial budget.
Original cell routes receive core Joule, stranded copper, graph-resistor and
all-triangle pinning heat once. Insulating magnetic cells require finite heat
receivers when pinning dissipates. End and field numerical losses remain
outside caloric deposition.

The combined work check uses the native stable graph/field/end storage
increment, external circuit work and actual finite enthalpy increment.
Source-bound graph resistors require exact original routes, distinct whole
thermal-cell volumes, material selections and physical conductor geometry.
A source refusal, work mismatch, exhausted budget or cancellation returns no
accepted circuit, magnetic or thermal history. Fixed geometry, constant
source-admitted magnetic coefficients and wholly solid caloric charts retain
the same first-order estimated MQS limitations as the field-driven adapter.

## Spatial hysteretic induction on the original finite thermal timeline

The default-off `magnetic-hysteresis` feature adds
`MagneticThermalNetwork::step_hysteretic_eddy_candidate` and
`circuit_enthalpy::eddy::thermal::heat::hysteretic_heat_assembly_candidate`.
They compose the original planar vector-history induction owner with the
existing finite-enthalpy residual and complete heat assembly. Every thermal
trial restarts the same accepted potential and all accepted vector branch
polarizations. Conductivity is resolved from the original material card at the
actual trial temperature; fitted magnetic coefficients retain their original
sources, values and explicit constant-temperature validity. Temperature cannot
silently change the magnetic fit or erase remanence.

Native per-triangle Joule heat and irreversible pinning heat are deposited once
into their actual caloric receivers. Pinning heat from electrically insulating
core cells is included. A dissipating fixed region without a finite thermal
destination refuses the candidate. Original magnetic/thermal cell identities,
whole represented solid volumes, source selections and clocks are checked.
Passive coolant cells can belong to the complete heat assembly without lying
inside the magnetic mesh; existing contact, conduction and radiation owners
continue to resolve at the same endpoint.

The work account uses the native stable magnetic energy increment:
`delta W_magnetic + delta E_thermal + D_numerical = W_impressed + Q_external`.
The report's existing `deposited_joule_j` field includes both physical loss
channels, while `electrical.joule_heat_j` and `electrical.hysteresis_heat_j`
retain their separate values. Backward-Euler numerical dissipation never heats
the material. Complete heat providers retain the original enclosure-balance
gate and refuse independent graph electrical work without its circuit owner.
All source, work-budget, numerical and cancellation refusals publish neither
accepted history. No residual correction or extra empirical core-loss term is
introduced.

This is first-order, fixed-geometry, isotropic estimated MQS coupling with
constant fitted magnetic coefficients. Source-supported electrical transport
can vary with temperature; temperature-dependent hysteresis parameters,
magnetic entropy, phase-changing geometry and measured device validation are
not implied. The focused `hysteretic_eddy_enthalpy` tests use manufactured
scalar-play/P1 balances, physical heat separation and transactional refusals.

## Nonlinear tetrahedral current and finite thermal storage

`circuit_enthalpy::conduction::spatial::response` closes the original
`SpatialConductionResponse` and finite caloric cells at one backward-Euler
endpoint. Each proposed temperature/phase refreshes the same complete pinned
bulk sources, then the native tetrahedral owner closes their actual E/J
coordinates. Boundary voltage/current density and mixed equipotential
voltage/current electrodes retain their native meanings. Immutable geometry,
source identity, query schema, caloric charts and the selected source-coordinate
relaxation policy persist through feedback and continuation.

`step_response_contact_candidate` refreshes both actual contact temperatures,
declared mean temperature and phase coordinates while retaining the accepted
electrical coordinates as native solve seeds. Explicit native nonlinear
contact bindings own any voltage-jump/current-density response. Their uniform
trace admission, source bounds and dimensional constitutive gates remain
mandatory. Each contact requires its original adjacent-volume heat receivers
and an explicit physical partition. Native volume E.J and exact surface work
enter those original finite enthalpies once; constitutive defects are only
diagnostics. Independent terminal source work closes against that deposited
heat before publication.

`step_spatial_conduction_response_candidate` composes any complete
`HeatAssembly`, retaining its source samples, thermal contacts, layers,
radiation, reservoirs and finite links at that same endpoint. Caloric and heat
transport claims must be independent of electrical response/applicability;
this path supplies no reversible field storage. Every native basis solve,
including discarded constitutive and thermal stencils and final fresh
admission, consumes the original common field-work allowance. Source-domain
refusal, work exhaustion, late cancellation and failed balances publish
neither history. Original fixed-geometry and confined-phase restrictions apply.

Planar `ConductionResponseCircuit` retains the complete original bulk/contact
source bundle and relaxation policy across frozen-operator transfers. Its
report retains actual-field contact source evaluations beside bulk source
defects and native circuit work.

## Relaxed nonlinear electrical source closure

`ConductivitySolve::with_relaxation(weight)` explicitly selects a finite weight
in `(0,1]` for the original circuit and reciprocal field/motion conductivity
iteration. Finite-enthalpy coupling retains that selection through each thermal
trial, its remaining electrical-work budget and final joint admission.
Only bound electrical query coordinates are relaxed. Every next
bulk/interface resistance comes from the same pinned source at that query;
temperature, pressure, source identity, geometry and accepted state stay fixed.
Plain controls and weight one preserve the original arithmetic. Constant laws
still finish in one native trial. Intermediate queries count in
`source_evaluations`, and no query is performed after the native-trial budget
is exhausted. Original actual-field resistance and constitutive-power gates
remain mandatory; native Joule heat is counted once. Source-domain failures,
nonconvergence and cancellation remain transactional. This can resolve passive
laws with oscillating Picard iterates but supplies no general convergence
guarantee or automatic relaxation choice. Focused circuit regressions use
strictly monotone terminal laws with independent quadratic midpoint solutions;
the coupled regression deposits their original Joule heat once in a finite
body with source-dependent heat capacity and checks the original enthalpy.

## Finite-solid material potential work

The finite-solid material binding forwards direct elastic energy increments
through the same held natural reference as its energy, stress and tangent,
including the original reference-volume scale. Midpoint candidates use the
prepared solid's single-pass material interval for endpoint storage and its
separately accumulated change. Force-kick/contact drift, mixed free flight and
sustained-contact accounting carry that direct `strain_change_j` into their
original balance equations instead of reconstructing it from rounded totals.
Mixed free-flight kinetic/source/defect sums remain compensated through their
final residual, retaining small work that rounded aggregate totals would lose.
Each original material and total deformation still admits both endpoints.
The change introduces no new constitutive potential, damping, heat or
time-integration scheme; temporal and nonlinear-solve errors remain measured
against the existing physical tolerances.

## Stationary gravitational potential work

`solid_gravity::bind_potential_increment` transports the original native point
field's optional analytic potential difference into the existing tetrahedral
reference-mass quadrature. Uniform and central sources retain small work that
independently rounded absolute potentials would erase. Custom fields without
this capability retain sampled endpoint differences; no additional precision
or gauge-invariance guarantee is inferred for that fallback.

Finite-solid midpoint accounting and force-kick/contact drift now obtain
gravity work from this shared increment. Mixed rigid/solid free flight consumes
the native rigid receipt's potential increment. Endpoint absolute potentials
remain diagnostics. Gravity is still counted once as conservative source work;
the measured difference from the integrator's force work remains a numerical
defect, never heat. The midpoint integrator is unchanged and nonlinear fields
still require temporal and spatial refinement. Endpoint diagnostics no longer
assemble unused Newton tangents.

Reciprocal finite-solid trajectories retain explicit native source snapshots
at both configurations and integrate their direct interaction-potential
increments over the original cellwise or nodal P1 reference-mass quadrature.
Both source recoil and quadrupole orientation participate before rounded
absolute binding-energy totals are formed. Central snapshots retain G M;
quadrupole snapshots retain the original mass, inertia, G, extent and frame.
The native force, recoil, torque, kick/drift stages and domain/energy gates
remain unchanged. The new reciprocal regression covers both mass inventories
and source selections; complete coupled execution still requires the original
constellation dependencies.

The new public-adapter and native free-flight regressions compare zero and
1e20 m uniform datums, including cellwise and graded mass inventories, one-joule
lifting work and 0.125-joule kinetic/potential exchange. Whole-workspace and
full solid/coupling execution remain unavailable with missing constellation
dependencies; these integration regressions are retained for that lane.

## Nonlinear membrane damping and source-bound finite heat storage

`membrane_thermal::advance_candidate` composes the existing nonlinear
`DampedMembraneMotion` with `DepositionSimulation::advance_observed`. One
unique active caloric receiver belongs to each original reference triangle.
Its `SheetCaloric` retains the identical complete original sheet source and
must be the exact caloric law already borrowed by that thermal cell. Volume is
original area times reference thickness; mass is the original P1 corner-density
integral, including nonuniform reference inventories and constrained coordinates.
An initial 512-epsilon relative allowance covers arithmetic only. Subsequent
cell/material/body/segment, geometry, mass and caloric-instance checks are exact.

The native thermal owner chooses every conductive/event/source interval, also
bounded by the mechanical step ceiling. Analytically integrated per-triangle
Kelvin--Voigt heat becomes average power on that same interval. Original source,
conductance, ambient/reservoir, pairwise-transfer, transition and observer hooks
remain active. No instantaneous-Cp temperature increment, second thermal solver,
numerical-defect heat or support-impulse heat is introduced. The source chart
performs the original nonlinear enthalpy-to-temperature conversion.

Mechanical state is borrowed until the first actual interval, then cloned once;
its prepared inertia and damping buffers are reused across substeps. The thermal
owner supplies whole-call transactional cloning. Motion, thermal inventory,
observer history and energy records publish together after final checks;
refusals and cancellation expose none of their half-steps. Compensated accounts
retain physical heat, held-load work, support impulses and signed numerical
defects. Absolute substep mechanical defects and per-cell heat-handoff errors
have whole-call budgets, preventing cancellation of opposite errors.

The combined physical residual is mechanical storage change plus native thermal
storage change minus incoming inventory and external work/heat. It remains an
uncorrected numerical residual; native mechanical defect is separately disclosed,
never deposited. Frozen elastic/viscous coefficients, natural reference and
geometry must stay within the caller's explicit fully solid temperature interval.
This is one-way finite self-heating, not evolving-temperature material feedback,
reversible thermoelastic work, phase-dependent mechanics or a stability/accuracy
certificate. Native transport order and nonlinear split error require refinement.



## Explicit rigid-impact loss heating in original finite caloric cells

`friction_heat::directional::impact::heat::apply_source_rigid_elliptic_impacts_with_heat`
composes the original simultaneous directional impact owner with
`DepositionSimulation::apply_heat_j`. A unique-cell-ID map explicitly assigns
nonnegative fractions of the one GLOBAL native impact dissipation. Fractions
must total at most one; compensated admission retains even sub-ulp overspend.
The remainder stays in `ImpactHeatReceipt::unallocated_energy_j`. This policy
does not infer a heat-conversion fraction or attribute total collision loss to
individual contacts, normal/restitution loss, friction, wear or radiation.

The original impact/source solve executes once using incoming conditions and
the complete native body/contact mobility. Its dissipated joules reach the
existing finite cells through their actual masses and original nonlinear
caloric laws. Cp/H, Cv/U and latent coordinates retain their source identity and
storage basis. An H receiver retains its pressure-work convention; the combined
account does not relabel H or declared sensible storage as internal energy.
No snapshot heat capacity, impulse-derived pressure, contact area or event
duration is introduced. Empty/zero allocation retains the original mechanical
endpoint and leaves thermal inventory unchanged.

The event occurs at the receiver's accepted physical time; the caller owns
mechanical event alignment and must process due thermal events first. Original
geometry, contact links, source watts and the thermal clock are preserved.
Native inventory limits and slice consistency are checked before staging.
Bodies, source histories and thermal state remain candidates until caloric,
partition, combined-work and final cancellation checks pass. External effects
inside caller callbacks remain the caller's responsibility.

The combined residual sums the original mechanical defect, original caloric
storage defect and explicit partition rounding defect:
`Delta K + Delta E_original + unallocated_loss - support_work`.
Each upstream owner retains its own bounds. No residual is credited as physical
heat, and represented cell-energy errors cannot cancel away inside the original
thermal audit. The partition helper allocates one map per event; no extra
thermal solver, physical time integration or constitutive model is added.

Three native tests execute the exact production partition/deposition helper
against the original thermal crate: distinct analytic Cp(T)/Cv(T) and latent
storage with fixed time/identity, fraction and representability limits, and
late cancellation/domain refusal with deterministic retry. Three full native
mechanical integration regressions are included for an independent 1.36-J
collision, empty/zero mechanical equivalence and joint body/history/caloric
rollback, but were not executed here because external workspace dependencies
are absent. The complete thermal production library compiles; its broader unit
test target separately remains blocked by two pre-existing advection observer
fixtures missing required `event`/`step` methods. No full-workspace verification
claim follows from the focused helper execution.

## Evolving native joint motion with finite-cell friction heat

`friction_heat::joints::trajectory::advance_joint_friction_heat_trajectory` composes the original constrained native half steps and `propose_joint_friction_heat` over a nonempty, strictly increasing, caller-authored time grid. Actual native midpoint states and physical time enter the temperature/material source callback. Every interval recomputes joint geometry and solves friction strength, native recoil, interface temperatures and original finite-cell heat together. Only the native rigid owner advances poses.

The original thermal clock determines each represented interval. Known drive/load switches belong at grid endpoints. Controls bound steps, interval duration, aggregate nonlinear trials, native/thermal inventories, absolute accumulated reaction work, signed total energy residual and the sum of absolute interval energy residuals. Only final body/thermal/source receipts and at most six cumulative joint/link work accounts are retained; memory does not grow with the number of steps. Late failures or cancellation publish none of the intermediate candidates.

Cumulative work accounts distinguish original wrench work, axial drive work, continuous friction loss, converted heat, explicit nonthermal exports, outside heat and measured bilateral residual work. Endpoint mechanical change uses paired native momentum/velocity increments with the original full inertia and gravity potential; caloric change uses each original cell's mass and stored coordinate. No balance residual is added to heat. After a phase with equal actual receiver temperatures, that caloric temperature seeds the next interface solve; it remains a Newton guess, avoiding spurious tiny energy corrections when friction stops.

The thermal routes remain explicitly reduced, fixed-geometry bearing/rail receivers; no moving footprint, material remap, surface recession or flash-temperature law is inferred. Backward-Euler heat makes this a first-order thermal composition despite symmetric mechanical splitting. Caller time refinement remains required: exact static balance belongs in joint drive ports, while separately split body wrenches can introduce artificial axial slip/heat at finite step size. A small energy residual is not a trajectory/heating accuracy certificate.

Four direct regressions cover event-aligned stopping/reversal poses with 1/2 J and 7/18 J loss, eight-step temperature-dependent friction against an independent quadratic recurrence, an authored drive switch with gravity/wrench/outside-heat and nonthermal-export accounts, and late source/cancellation or malformed-grid atomic refusal.

## Hardness-limited sliding, wear and finite caloric feedback

`friction_heat::plastic_wear::PlasticSlidingHeatContact` binds an accepted
`SourcePlasticSlidingWear` checkpoint to its original ordered material cards,
source-resolved Cp/density charts and an existing `FrictionHeatRoute`.
`step_candidate` solves original cell-energy increments and the massless
interface temperature using the existing `FrictionHeatStep`,
`FallibleFiniteDifference` and Newton--Krylov owners. No parallel thermal
balance or nonlinear solver is introduced.

Each residual evaluation, directional probe and final replay starts from the
same accepted source and complete indentation, shear and wear history. The
original caloric objects give separate body temperatures; the normal/wear
temperature factory refreshes original elasticity and each body's pressure
hardness at those temperatures, while original friction and Archard coefficients
follow the solved interface temperature. The selected normal boundary is held
approach, maintained force or a held finite-stiffness device. Native Coulomb
motion uses the resulting actual cell pressures, and only its solved slip
length drives each body's original Archard removal and subsequent normal
re-equilibration. Thermal coefficients are endpoint values; wear remains the
native beginning-pressure split, requiring physical time/motion refinement.

The original finite-cell geometry and endpoint conductivities determine signed
friction-interface heat into both bodies and a solved common interface temperature.
No prescribed friction-heat fraction or arithmetic temperature mean is inserted. The source callback
owns original conductivity claims and other transport, retaining its query and
selection evidence in the returned payload. The route excludes its old direct
link to avoid duplicate conduction. Only actual native Coulomb dissipation
enters the interface heat source. Signed temperature/recession parameter work
remains separate; its explicit external-work budget is the sum
of absolute parameter contributions. Native contact identities and the thermal
cell, interface and aggregate balances have independent acceptance tolerances.

`with_plastic_work_heat(PlasticWorkHeatPartition)` optionally converts declared
fractions `f_A, f_B >= 0`, `f_A + f_B <= 1` of the two native normal plastic-work
increments into heat in the original receiving bodies. The same coupled trial
uses `Q_A = f_A*delta_Wp` and `Q_B = f_B*delta_Wp`; its endpoint temperatures
therefore feed back into original hardness and pressure. These are direct bulk
sources, never extra friction-interface heat. Accepted cumulative plastic work
is not reheated. The callback must exclude both exported friction and configured
plastic heat, while other authored heat at the same cell is retained.
The default remains zero conversion. Fractions are explicit reduced-model
assumptions, not a universal Taylor--Quinney constant, a calibrated material law,
or an inferred split of relative indentation between the two solids. See
Knysh and Korkolis, *Mechanics of Materials* 86 (2015), 71--80,
<https://doi.org/10.1016/j.mechmat.2015.03.006>.

The report separates total incremental normal plastic work, per-body plastic
heat, its nonthermal remainder, and other external heat. The remainder is not
claimed to be recoverable storage. The lower-level thermal report counts direct
plastic sources as external to its thermal-only boundary; `other_heat_j`
removes those sources from the coupled report's external heat. The combined
local identity is `delta_Ucontact + delta_Uthermal + Wp_nonthermal = Wnormal
+ Wparameter_temperature + Wparameter_recession + Wfriction + Qother`.
This retains the externally supplied material-parameter work of a thermal
refresh; its conversion does not establish reciprocal thermoelastic closure.
Focused G0/G1/G3 regressions cover conservative partition admission, asymmetric
body heating, additive authored heat, independent two-cell balances and no
reheating of accepted indentation history.

Admission retains the original receiving body IDs, cell geometry, material IDs,
masses, caloric objects and clock. Original caloric and normal sources must share
ordered card identities, complete fixed conditions and quantity axes, and the
selected density claim. Accepted bulk source temperatures must equal the
thermal checkpoint; melted material refuses. The original active grid-cell
mask must have the thermal face's area and remain unchanged after both
temperature refresh and wear. Pressure may redistribute within that mask.

This is an explicit fixed-footprint, fixed-thermal-inventory approximation.
Cumulative removed volume is bounded against each ORIGINAL receiving cell
volume, never a renewed per-step denominator; zero permitted removal is valid.
The fraction bounds an omitted effect, not its thermal error or experimental
validity. The default heat step does not export material; the optional endpoint
operation below removes actual mass and carried enthalpy. Debris evolution
remains outside this model.
Reciprocal thermoelastic free energy, a resolved subsurface plastic heat field,
moving contact coverage and storage in caller-supplied tangential body response
are not closed by this owner. Its local contact identities and fixed-inventory
heat balance imply no complete body-energy, entropy or passivity theorem.

The candidate seals successor source, contact history and thermal state for
joint adoption. Bounded trials, Newton iterations and thermal inventory,
source-domain or footprint refusal, excessive cumulative removal or parameter
work, nonconvergence and cancellation publish no partial physical history.
The report retains both native mechanical operations, source callback evidence,
friction heat, normal/plastic/parameter accounts and cumulative removed volumes.

### Depth-dependent hardness in live caloric feedback

`friction_heat::plastic_wear::depth::DepthSlidingHeatContact` binds the
original `SourceDepthSlidingWear` history to the same finite-cell caloric
engine. Each trial resolves each body's original H(T,z) at its own temperature
and exposed material depth, with independent elastic properties. Original wear
and friction sources use the solved interface temperature and actual pressure.
Actual unsigned Coulomb slip drives each body's own depth-hardness integral;
signed shear travel, retained indentation and removed inventories survive.

Held approach, maintained force and a fixed finite-stiffness drive use the
native boundary continuations. Temperature refresh, recession and exposure of
the new hardness cap retain three separate normal-work receipts. All three
enter the work audit, sum-of-absolute-parameter-work limit and optional explicit
plastic-heat conversion; parameter work remains separate from friction heat.
The complete thermal, depth, shear and contact candidates are adopted together.
The ordinary and depth owners share the existing Newton/caloric engine and
source/footprint/inventory admission, with no second contact or thermal law.

Fixed and evolving active footprints retain their original geometry limits.
Contact resistance is the explicitly supplied fixed route resistance. This
owner does not add adaptive depth marching, lateral temperature fields or a
reciprocal thermoelastic free-energy model. Thermal mass and geometry remain
fixed within each interval, and cumulative removal uses the original receiving
volumes. The regression definitions cover an independent coupled quadratic,
actual reversal under force and spring control, all three work accounts and
transactional source, removal, iteration and cancellation refusal. Native Rust
compilation and test execution were unavailable during this implementation.

### Conservative export after depth-hardness sliding

`DepthSlidingHeatContact::export_wear_candidate` removes newly worn material
through the same native receiver-export path as ordinary hardness-limited
sliding. Original sourced density sets removed mass, and the thermal inventory
owner carries the material's current signed enthalpy out with that mass.
Remaining specific enthalpy and temperature stay unchanged, including for a
negative enthalpy datum. This material outflow is never deposited again as heat.

The export checkpoint excludes wear already present at binding and prevents a
second export of the same volume. Original cumulative removal limits survive
export. Explicit shortened prismatic geometry and affected contact patches
retain the ordered bodies, cross sections, contact area and tangent registration;
subsequent heat steps use the remaining masses and new conduction distances.
Hardness profiles, actual recession, indentation and shear history remain with
the original contact owner. Geometry, inventory or cancellation refusal returns
no partial state. Debris evolution and resolved remeshing remain outside this
bounded heat-then-export approximation. Three additional regression definitions
cover signed enthalpy, the next finite-capacity heat balance, repeated/prebinding
export, the cumulative envelope and transactional refusal; native execution
remains unverified because the runtime was unavailable.

## Static breakaway and kinetic sliding through joint heat and motion

Existing joint-friction thermal source callbacks now accept explicitly authored
static/kinetic native port variants. Every nonlinear temperature trial resolves
the selected native Coulomb/stiction law and supplies its continuous friction
work to the original finite thermal network. Holding generates no artificial
friction heat; released and already sliding joints use kinetic capacity. Static
capacity, kinetic capacity and conductivity may be resolved at the same
candidate interface temperature. Ambiguous mechanical continuation remains a
refusal, and accepted inputs remain borrowed until a complete candidate passes.

Three cross-module tests distinguish hold from slide under identical load,
solve temperature-dependent sliding friction against a finite two-cell thermal
response and an independent closed-form balance, and evolve an authored motion
grid through hold, breakaway and continued sliding below the static limit.
The last case checks actual poses, endpoint rates, drive work, kinetic change,
integrated friction loss and thermal storage together. No new capture speed,
implicit presliding history, inferred heat fraction or imposed thermostat is
introduced.


## Scenario pressure histories on native fluid reservoirs

`scenario_pressure` connects declared `IncompressibleFlow/PressureOutlet`
boundaries to the original specialized D3Q19 and generic D3Q19/D3Q27 open
boundary owners. `PressureReference` explicitly distinguishes signed gauge
pressure from positive absolute pressure with a supplied positive reference.
Both use `MaterialScaling::density_lattice_from_pressure_perturbation`; the
adapter does not infer an acoustic equation of state or a pressure-work model.

Uniform values and `TimeSignal` histories retain SI dimensions, the original
interpolation/clamping rules and cancellable validation/evaluation. Scalar
pressure is invariant under a checked frame transform. Spatial profiles require
an explicit geometry binding and refuse this uniform-time lane. Selected native
pressure profiles also refuse, including unchanged dormant uniform targets;
unselected controls and profiles remain active.

A six-face batch prepares every target before one native rebind. Sampling,
reference, units, density, native-kind or cancellation refusal preserves all
controls and fluid state. This updates held controls; native stepping retains
collision, force correction and reservoir transport. G0/G3 tests in
`tests/scenario_pressure.rs` exercise both native quadratures, actual boundary
density after steps, gauge/absolute equivalence, histories, profile precedence,
transactional refusal and cancellation. No new dependency or solver is added.

## Changing contact area with physical conduction and exposed cooling

`FrictionHeatStep::evaluate_with_area` and its entropy variant stage mechanical
area first, then use the native contact-only implicit endpoint ingress before
resolving material transport. The callback sees actual covered/exposed surface
areas and a suppressed direct link. Half-cell fluxes use original conductivity
times the new area divided by the original distance; conductivity itself is
never rescaled to imitate geometry. The fixed and changing-area APIs share
their physical transfer, energy and entropy assembly.

The accepted thermal state and `FrictionHeatAreaTrial` successor route retain
the same new coverage and source-derived series conductance. Actual ambient
cooling follows native coverage without a compensating source. Areas must be
positive and fit both original prismatic faces and aggregate surface coverage.
Complete separation refuses because a vanished massless interface has no
unique temperature. Cells, masses, caloric objects and geometry remain fixed.

This is a uniform-face endpoint discretization requiring time/spatial
refinement, not lateral heat spreading or a swept-contact model. Discrete
active-cell changes may prevent a nonlinear root from converging. Three
focused regressions specify independent shrink/growth energy and entropy
balances, callback-visible coverage and source-k rebinding, exact-target event
order, and transactional late failure including a caught contact-update error.

## Spatial thermal networks in impact pressure playback

`render::plate::impact::audio::network::ThermalNetworkSource` retains an original
heated or liquid-film impact source, its `DepositionSimulation`, thermal observer,
reservoir/cell bindings, step controls and contact-material resolver. Its
`ImpactSource` implementation calls the original `step_with_thermal_network`
once per mechanical sample. It does not first take an adiabatic step, reconstruct
caloric properties, or deposit mechanical dissipation a second time.

`ImpactPressureRenderer` therefore observes the accepted network-coupled
trajectory through its existing radiation and propagation histories. The same
source can pass through causal decimation, pressure mixing and WAV streaming.
`last_mechanical_frame` retains the original `NetworkImpactFrame`, including
thermal substeps, external heat and the combined mechanical/thermal account.
Source, thermal cells and observer are available read-only between callbacks.

The original coupled owner admits complete bindings, matching physical clocks,
caloric identity, source domains, thermal and mechanical budgets, and cancellation.
A rejected physical sample publishes none of those states; the pressure host
retains its existing failed-callback poisoning contract. External effects in a
captured resolver are not rolled back. Thermal splitting, fixed geometry,
unloaded-only contact refresh and one-way acoustic observation keep their
original limitations. Construction retains inputs without advancing or admitting
a physical step. The focused regressions compare direct/network-observed states
and audits, conduction-driven viscosity and pressure, and late-observer refusal.


## Native nodal heat and continuing tetrahedral dynamics

`solid_enthalpy::spatial::SpatialSolidEnthalpy::resolve_at` binds an explicitly
timed native nodal enthalpy field to the existing source-backed solid state.
`SpatialSolidState::advance_transient_candidate` passes that state to the
original tetrahedral Newmark thermal-update owner. Displacement, velocity,
acceleration and step history continue; the original operator owner retains
reference mass and rebuilds changed stiffness.

The previous mechanical material/eigenstrain state and clock must match the
accepted run. Both endpoints retain their reference mesh, frame, material
assignment, caloric chart, source laws, strain envelope and reference inventory.
Untimed equilibrium states, stale history and incompatible replacement states
refuse. Native support, load, cancellation and energy checks remain active;
failed proposals publish no heat or motion state.

This is prescribed one-way thermal loading on stationary reference geometry.
Signed thermal/stiffness work and Rayleigh loss retain the native mechanical
accounts; none is silently deposited into or removed from nodal enthalpy.
The caller accepts heat and motion together and resolves thermal events and
time refinement explicitly. It does not add reciprocal thermoelastic feedback
or a moving-mesh heat law.

`SpatialSolidEnthalpy::advance_transient` completes that handoff as one candidate:
it admits the caller's accepted nodal heat against the retained mechanical fields,
solves the original heterogeneous conduction/radiation owner, resolves its actual
thermal endpoint and continues native motion using exactly the same supplied
timestep. It never reconstructs that duration by subtracting rounded clocks.
All source-bound heat consumers share the original conductivity/inventory binding.

`SpatialSolidTransientCandidate` returns the original `EnthalpyHeatStep`, timed
source state, native `TetTransientCheckpoint3` and original mechanical work report
together. The checkpoint retains the trajectory and existing shared matrices
across caller-owned material lifetimes; rebinding its exact endpoint problem
continues the accepted clock and step budget. A late structural or cancellation
refusal publishes neither heat nor motion. Original radiative energy reports and
mechanical work remain separate. Regressions cover exact-duration clock rounding,
two-step continuation/replay and structural refusal after the heat solve.

## One scenario control transaction for fluid inlets and outlets

`scenario_open_boundary::rebind_open_faces3` composes the existing velocity and
pressure preparations, then calls the native boundary owner once. Its six-slot
`ScenarioBoundary3` selection admits one control family per face; all selected
histories share one physical time and material scaling. Unselected faces keep
their complete existing controls. Selected spatial profiles retain each
original adapter's explicit refusal rule.

A late pressure failure, velocity-speed refusal, native kind/collision refusal
or cancellation cannot advance only one side of the flow. Both specialized
D3Q19 and generic D3Q19/D3Q27 owners retain their own collision, reconstruction,
force correction and step accounting. The transaction changes held targets;
it does not step the fluid, infer a geometry mapping, enforce total mass-flow
balance or claim a pressure-work certificate. Fixed six-slot preparation
copies no histories or population fields. G0/G3 regressions in
`tests/scenario_open_boundary.rs` exercise mixed native velocity/density
response and all-or-none refusal, including the original collision Mach bound.

For the shared native D3Q19/D3Q27 owner,
`scenario_open_boundary::lattice3::{step_serial, step_pooled}` samples the same
typed histories and advances the original fluid operator as one transaction.
Candidate targets drive collision/stream reconstruction and measured reservoir
transport. Held controls, populations and completed-step count publish together
only after the final native checkpoint; a late local-Mach refusal or cancellation
retains the entire accepted prefix. `BoundaryStep3` returns the accepted physical
samples alongside the original native mass/momentum transport receipt.

All generic mixed and pressure-only adapters reject a selected installed pressure
profile, even when a sampled uniform target equals its dormant default. Unselected
profiles keep their original face fields and still drive the step. The caller
retains the physical sample time and material scaling; no pressure-work or thermal
energy certificate is inferred. Focused tests exercise changing histories on both
lattices, serial/pooled parity, late refusal, cancellation/retry and profile
ownership across the mixed and pressure-only entry points.



## Hardness-dependent changing contact area in finite heat feedback

`friction_heat::plastic_wear::PlasticSlidingHeatContact::new_evolving_area`
opts the existing source-backed sliding/Archard/thermal owner into changes of the
positive-pressure cells on its original finite-gap grid. The complete projected
grid area must fit both uniform thermal face areas. Every nonlinear trial starts
from the exact accepted plastic, signed tangential and removed-volume histories;
recruitment and opening never reset or heal them. Original material cards still
supply temperature-dependent hardness, dry friction and caloric state separately.

Frictional dissipation and native Archard increments use the refreshed beginning
pressure. The final worn pressure footprint determines actual endpoint contact
area, source-derived half-cell conductance and ambient exposure through the same
backward-Euler contact-patch ingress as all other thermal links. Reports retain
the previous, refreshed and final areas and the exact accepted conductivity.
Successful candidates own their successor route and active mask, so the next
interval validates the accepted coverage and physical conductance. The default
constructor continues to require its original active cells.

This resolves total area between two uniform finite caloric receivers; it does
not resolve lateral thermal fields, constriction resistance, contact residence
or swept patches. The native beginning-pressure wear/endpoint-area split requires
physical time and grid refinement, and an activation discontinuity may refuse
nonlinear convergence. Complete detachment refuses an undefined massless
interface temperature. This base adapter still holds thermal mass and geometry
fixed within its explicit original-inventory removal envelope. Source regressions
exercise thermal-softening recruitment followed by wear-induced opening, actual
friction/volume/heat accounts, original irreversible histories, incomplete face
support and late source refusal using independent two-cell elastic/heat references.

## Pressure-dependent resistance in hardness-coupled friction heat

`FrictionHeatRoute::from_link_with_resistance` admits two finite nonnegative
area-specific resistances around the existing massless heat-generation plane.
Each thermal leg is `A / (d_i/k_i + R_i)`; the accepted direct link retains the
same series sum `R_A + R_B`. The original direct path is suppressed once during
the native heat solve. Zero resistance preserves the original half-cell
arithmetic. The declared resistance split locates the generation plane; it does
not prescribe the resulting heat fraction or add an interface heat store.

`PlasticSlidingHeatContact::with_sourced_contact_resistance` pins the canonical
thermal-resistance claim from the original ordered interface, original fixed
conditions and explicitly typed static-pressure axis. Every thermal trial uses
the solved interface temperature and mean pressure over the actual positive-
pressure cells of the worn native endpoint. Equal native cell areas make this
mean exactly the resultant/loaded-area convention, without a nominal-area or
pressure/hardness area substitution. Sources calibrated to different pressure
or temperature conventions are outside this reduced model.

The same native cell-energy and massless-interface residuals own all transfers.
Accepted reports retain the actual resistance query and split; the successor
thermal patch and route retain its series conductance. Adaptive continuation and
actual material export preserve this resistance/source checkpoint. Claim support,
identity, geometry or nonlinear refusal publishes no partial physical state.
The model still has two uniform bulk receivers, beginning-pressure sliding/wear,
endpoint heat coefficients and no resolved lateral constriction or third-body
storage. Physical time/grid refinement remains necessary. Focused regressions
cover an independent asymmetric quadratic heat oracle, actual worn pressure,
export/continuation, exact zero-resistance replay and late source refusal.

## Runtime controls through pressure observation

`PressureRenderer` admits and applies batches of the existing `ControlDelta`
inputs at the current physical source boundary. Validation never advances time;
a refused batch leaves held inputs, pending events and control logs unchanged.
Producers without an explicit compatible input topology refuse nonempty batches.

Boxing and borrowing a producer preserve its input sink. `ScheduledRenderer`,
`EnsembleRenderer` and `DecimatedRenderer` forward controls to their original
retained mechanics. Decimation does not scale force, reindex source voices, move
source events, or reset the pressure filter. Both clocks remain unchanged by an
input update, and a poisoned producer refuses controls as well as rendering.

Live scheduled inputs apply immediately. A pending event on that same source
sample still runs in the schedule's existing stable order before the next sample,
so a scheduled assignment to the same input takes precedence. Live admission
reserves log capacity for the new batch and every remaining scheduled record;
the next pressure callback does not allocate merely because live logging consumed
the schedule's original capacity. This is a cold control operation, not a new
allocation-free control-thread guarantee. Schedule-only applied-event history
remains separate from the context's complete control log.

`ImpactPressureRenderer` accepts `SetModalForce` for its sole source slot zero.
The mode index, N/sqrt(kg) units and force ceiling are those of its original
`ImpactSource`, including a retained spatial thermal network. The complete batch
is checked before any held force changes; no mechanical, caloric, material,
radiation or propagation history is rebuilt. A force update uses the same held
input path as the concrete renderer's original `set_forces` operation.


## Exact reciprocal longitudinal piezoelectric storage

The longitudinal rod exposes the complete original quadratic storage block in
[extension, momentum, charge] order:
Q = [[k_D, 0, -theta/C], [0, 1/m, 0], [-theta/C, 0, 1/C]].
Its same source-bound energy and amplitude checks remain authoritative.
Native descriptor and prepared solvers can use analytic charge derivatives;
they do not need dimensionless finite-difference probes into microscopic
electrodes. The continuous storage tangent retains both reciprocal cross terms.
The quadratic declaration is exact for this fixed-state, small-signal rod and
is not a linearization of ferroelectric, hysteretic or nonlinear storage.

The G0 regression in tests/piezoelectric.rs compares every quadratic effort
with the original material response, independently reconstructs stored energy,
checks reciprocity and the analytic directional action, and verifies refusal
outside the original strain envelope. This adds no measured material
validation or hard real-time timing claim.


### Native circuit charge scales

The piezoelectric circuit now declares its complete reciprocal quadratic
differential block to native descriptor admission and stepping. Original
circuit storage, charge/strain cross terms, fixed-charge stiffness and inertia
are retained; algebraic multipliers carry zero storage. Analytic Jacobians
avoid generic finite differences of physical charge, including microscopic
electrodes. Material support and independent closed-work acceptance remain
authoritative. Regression definitions in tests/piezoelectric/circuit/charge_scale.rs
compare subfemtocoulomb response with independently eliminated loaded-rod
equations, both poling signs, zero coupling, domain refusal and retry.
These authored ceramic-scale fixtures are not measured material calibration.


### Prepared exact longitudinal current drive

`piezoelectric::prepared::PreparedPiezoelectricCurrent` prepares the original
three-coordinate [extension, momentum, free charge] system with the shared
`fs-phs::linear_exact::prepared::PreparedQuadraticFlow`. It retains every
reciprocal material entry and consistent rod inertia. Cold preparation has
explicit dimension, work and memory budgets; warm steps reuse fixed scratch
and the cached full flow and voltage integral. They perform no owner allocation,
Newton iteration, source clone or new matrix exponential.

Held force/current steps include the open-circuit case at zero current.
Mechanical supply is force times actual extension change; electrical supply
is current times the integrated trajectory voltage. The mean sensor voltage
is not an endpoint trapezoid. A computed conservative full-window energy bound
checks both original strain and field limits, including excursions between
equal endpoints. Independent charge/work acceptance and final cancellation
precede the sole physical state publication. The bound may refuse a physically
admissible step conservatively; it is not an interval-arithmetic certificate.

Synthetic G1/G4 definitions in tests/piezoelectric/prepared.rs compare an
independent held-current forced oscillator and its exact voltage integral,
repeated open-circuit frequency/charge, hidden interior excursions, and atomic
cancellation/retry. Exactness applies to the fixed linear held-input model
up to floating-point error; the existing spatial rod reduction, material
validity and absence of a hard real-time benchmark remain explicit limits.


### Moving-support spatial piezoelectric sensing and harvesting

`TetPiezoelectricTransient::step_with_support_motion_candidate` accepts one
endpoint velocity for every original constrained mechanical DOF. Support
coordinates retain their topology, while displacement follows the same midpoint
kinematics as the free coordinates. Full consistent-mass, constitutive and
reciprocal electrical coupling columns remain in the original solve. The
existing initializer still requires stationary supports initially.

Independent mechanical supply is applied-force work plus support-reaction
work. Returned midpoint support forces include inertia and intrinsic material
memory, with zero entries on free coordinates. Missing, duplicate, free,
nonfinite or unrepresentable support values refuse. Stationary APIs refuse
to silently stop a moving support; an explicit endpoint velocity accounts
for its deceleration and returned work.

The same motion is available through the resistive load owner's
`step_with_loaded_support_motion_candidate` and the reactive R/L/C owner's
`step_with_support_motion_candidate`. Original resistor heat, reactive storage,
field storage and source work participate in the enclosing balance.

The finite-temperature spatial owner exposes the same support motion through
`step_with_thermal_support_motion_candidate` and
`step_with_material_thermal_support_motion_candidate`. Every thermal proposal
starts from the accepted mechanical, electrical and material-memory history.
Actual support work enters joint supply once, while resistor and intrinsic
ceramic losses reach their explicitly bound finite caloric cells. The report
retains complete midpoint support reactions. Source-domain, inventory,
unchanged-coefficient, cancellation and physical-energy gates govern joint
publication. Without external resistors, nonlinear thermal trials reuse one
immutable native moving-support candidate.

The three synthetic journeys in
`tests/piezoelectric/spatial_transient_load/thermal/relaxation/support_motion.rs`
compare independent Maxwell/Debye, resistor-feedback, caloric and inertial
references through acceleration, continuation, reversal and stopping. Native
Rust compilation/execution remain unverified under the editing environment's
disk/memory limits; syntax/format and independent analytic checks were completed.
G1/G4 definitions in spatial_transient/support_motion.rs cover analytic rigid
inertia/reactions, independently condensed base-excited sensing, resistor and
capacitor harvesting, and transactional refusals. Small-strain assumptions,
source validity and temporal/spatial refinement obligations remain.

The spatial midpoint owner evaluates kinetic increments through its unchanged
full consistent mass and the shared `fs-math::quadratic` kernel, including all
prescribed-support and free-coordinate cross terms. Elastic/electric increments
use the original local `C^E` and `epsilon^S` tensors at the recovered endpoint
fields. Native Maxwell/Debye storage increments survive integration into the
cell account; their independently evaluated physical heat is unchanged.
Applied-force, support-reaction and electrode work retain endpoint low words.
The ordinary, resistive, reactive and finite-caloric gates consume those
independent increments before scalar diagnostic rounding. They neither subtract
large absolute kinetic inventories nor widen tolerances with background motion.
External capacitor/winding storage remains separate from conductor heat, and
the shared internal thermal-driver boundary retains storage minus source work
until the caloric balance is complete. All original field, charge, source,
amplitude, cancellation and publication checks remain authoritative.

Focused G0 definitions in the private `transient/energy.rs` exercise coupled
tensor cancellation, endpoint work and a small storage/heat exchange. G3
definitions in the existing support-motion suites carry a small transverse
kinetic change beside large common rigid motion through native sensing,
resistor/capacitor harvesting, intrinsic ceramic loss and actual finite caloric
cells. These synthetic dynamic-range fixtures do not establish material
calibration, relativistic physics, integration error bounds or performance.

## Annular no-slip supply passages

`FixedChannel::Annular` carries the original `AnnularChannel` into the shared
two-sided and DG-backed aperture supply connection. Both physical radii, all
radial velocities, the inner and outer no-slip conductances, and the actual
gas source survive the connection. Pressure condensation, stored fluid energy,
physical viscous heat, numerical loss and gas-intervention work use the same
native momentum owner as the standalone annulus. No equivalent circular
radius or hydraulic-diameter resistance replaces the annular solution.

The common channel API also retains independent axial inner/outer wall drive;
coupled supply constructors retain their existing stationary-wall semantics.
Geometry is fixed and concentric, the fluid is incompressible and Newtonian,
and axial flow is fully developed. Rotating cylinders, normal wall motion,
eccentricity, entrance effects and transition are outside this model.

The existing DG supply startup/reversal/passive-release battery includes the
annular geometry. `tests/annular_channel.rs` exercises momentum-preserving
gas changes through the common channel API alongside the original analytic
Poiseuille, logarithmic Couette and narrow-gap checks.

## Prescribed-temperature moisture histories

`MoistureSlab::at_prescribed_temperature` and
`MoistureCartesian::at_prescribed_temperature` propose a uniform-temperature
event at the accepted clock. They retain every cell's dry-basis water content,
dry skeleton, reference geometry, time, definition identity, and original
density/stress-free datum. They re-evaluate the originally pinned diffusivity
claims at the new absolute temperature and every occupied moisture state;
Cartesian fields retain all three explicitly selected principal sources. Any
unsupported cell, source-domain hole at an occupied endpoint, invalid
temperature, or inadmissible diffusivity refuses without publishing a state.
The transport anchor is an occupied cell, so the historical reference moisture
does not unnecessarily constrain a new temperature cut of the transport source.
Reference density is never queried to create or remove either inventory.

`temperature_kelvin` reports the current environment; `transport_context`
exposes its pinned source bundle, while `reference` remains the original datum.
Ordinary and SDIRK stages, finite-atmosphere/solution participants, paired and
adaptive reservoir controllers, layered sorption admission, and sorption-energy
admission all use the current context. Existing boundary laws and reservoirs
must independently qualify the same temperature. Kirchhoff paths are prepared
from the current source context on the next proposal; no old-temperature chart
is retained. Return temperature events use the same original claim pins.

Each transport interval remains isothermal: callers cut intervals at prescribed
temperature changes. This API supplies neither a solid heat equation nor the
heat needed for the imposed event, and does not interpolate an environmental
history across an unqualified source path. Native elasticity can resolve its
qualified current environment; the existing fixed-temperature hygroscopic
expansion integral still refuses a temperature change from its stress-free
datum. Thermal eigenstrain is never inferred from a moisture coefficient.
Authored G1 regressions exercise analytic hot-source flux, principal-axis
refresh, degenerate-front wetting, history invariance, late-cell refusal,
boundary/reservoir temperature agreement and the retained strain refusal.

## Second-order spatial moisture diffusion

`MoistureCartesian::advance_sdirk2` applies the existing Alexander two-stage
SDIRK tableau to the original three-dimensional finite-volume diffusion owner.
Both implicit stages resolve the pinned principal diffusivity laws and reuse
the native seven-point operator, checked CG, nonlinear controls and source
admission. Kirchhoff mode retains the original integrated D(M) face law,
including a degenerate dry state.

The first-stage contribution is evaluated from equal-and-opposite actual face
rates. The endpoint account contains (1-gamma)*dt times those rates plus
gamma*dt times the final rates through each of the six exterior faces.
The intermediate report is diagnostic and must not be added again. Dormant
surface transfer and impermeable directions cannot enlarge the admitted
second-stage moisture range. A nonphysical stage returns `StageRange` before
material evaluation; it never clips water or publishes a partial field.

G1 regressions specify anisotropic discrete-mode convergence, six distinct
boundary accounts, original nonlinear source refresh and analytic Kirchhoff
wetting. G0/G4 cases cover exact equilibrium, deterministic replay, dormant
reservoirs and cancellation. These are numerical checks with authored material
data, not experimental validation. Temperature and external faces are constant
within an interval. Ordinary `advance` remains backward Euler; no second-order
finite-atmosphere/solution coupling or global error bound is implied.

`MoistureCartesian::advance_sdirk2_adaptive_to` applies the same complete-duration
controller as the slab SDIRK path. The accepted two-half result uses the
second-order local-error factor (one third for equal halves), cubic step-growth
criterion and two-stage work charge per attempted endpoint. Range failures may
shorten the interval; source/domain failures remain terminal. Only accepted
RK-weighted six-face transfers contribute to the fixed whole-window water
budget. Rejected probes and accepted prefixes of a later failed call remain
private. Focused regressions specify the analytic error ratio, stiff-drying
recovery, complete-stage charging and late work-limit refusal.

## Simultaneous finite-wall, condensate and gas endpoints

`WallFilmThermalEnvironment::advance_implicit` couples original wall enthalpy,
volatile inventory and gas temperature in one shared Newton solve. The material
wall-film endpoint owns saturation, dry/wet complementarity, upstream vapor
enthalpy and finite-gas storage. The conduction endpoint owns wall calorics,
convection, radiation, source conductivity and the original Biot gate. Neither
an inner gas inversion nor an inner wall solve is performed during nonlinear
trials. After convergence each original caloric owner is reconstructed once,
and the complete phase/heat residual is checked again against both actual
returned states before the joint candidate is returned.

The finite wall receives condensation heat and supplies evaporation heat through
the same source-relative port as the material endpoint. Sensible gas-wall heat
has one reciprocal debit and credit; radiation remains an external enclosure.
Existing and newly formed film storage uses the liquid's original chart, never
the wall heat capacity. Exactly dry film has no thermal storage. At complete
evaporation, outgoing vapor still carries the wall-temperature-dependent
enthalpy. Wall melting/freezing plateaus remain primary-enthalpy states: latent
storage can change phase fraction at constant temperature without introducing
an apparent heat capacity or a duplicate heat source.

Numerical controls bound the shared Newton iterations and residual evaluations,
including the predictor, fresh converged trial and reconstruction rechecks.
Physical source, temperature, donor, mass and energy allowances remain those
of the original environment. Signed and absolute wall/film/gas balances retain
source-combination and boundary-port rounding; independent physical gates apply
after the numerical norm converges. Failure, source refusal, exhausted work or
cancellation returns no partial wall, film or gas state.

This opt-in endpoint path leaves the existing split advance and its continuation
semantics unchanged. It is first-order in physical time. The film exhaustion
time freezes the endpoint rate and requires timestep refinement for an accurate
continuous event time. Fully wetting, pure-liquid, fixed-area, fixed-gas-volume
and uniform-wall assumptions remain explicit. It does not add nucleation,
boiling, drainage, spatial gas gradients, displacement work or an entropy
theorem. Focused coupled regressions check actual finite-wall/gas feedback,
complete dry-out, retained latent wall storage and transactional refusal.

### Local discrete magnetic work in nonlinear voice coils

`speaker::magnetic` supplies a storage-owned symmetric magnetic chord for
`H_m = (lambda - psi_pm(x))^2 / (2 L(x))`. Endpoint-mean current and analytic
polynomial linkage/inductance secants close magnetic work in linkage and
excursion alone, including permanent-field and reluctance force. Independent
cone/panel inertia and reciprocal cavity energy keep their original quadratic
mean efforts and represented-midpoint energy remainder. The rule remains
well-defined at coincident displacements without dividing by their increment.

The original descriptor solver consumes the analytic endpoint derivative of
this complete rule. Constant `Bl`/`L` profiles retain the exact-quadratic
factorization. The existing thermal coupling observes the same discrete
current for coil-resistance heating. Source-support and arithmetic refusals
remain transactional; finite-step response still requires timestep refinement.
This displacement-dependent reversible law adds no current saturation,
hysteresis, skin/proximity loss or measured coefficient claims.

Focused regression definitions cover local and total work, independence from
unrelated mechanical momenta, reverse-step symmetry, endpoint derivatives,
coincident states, the constant-law path and invalid-support refusal.

### Physical momentum in reciprocal electric and magnetic gaps

`material_circuit::motion` confines nonlinear field work to the original
charge/linkage and gap coordinates. Its cubic uniform-gap energy uses an exact
average-vector-field chord with an analytic endpoint derivative. The mechanical
effort is the mean physical momentum divided by the original body/support mass;
unrelated circuit storage retains its own quadratic efforts. A full-state
Gonzalez correction no longer changes mechanical velocity to close field energy.
For a held current or voltage source, the actual force impulse matches the
analytic integral of the original quadratic Maxwell force, including reversal.
`motion::dielectric` uses the same average-vector-field rule on each original
piecewise-affine E(D) segment. Exact polynomial quadrature splits only at source
knots; kinetic effort and unrelated nonlinear capacitor secants retain their
physical coordinates. Its linear material limit matches the uniform-gap owner.

`motion::coupled_gap` uses symmetric fixed-gap/fixed-linkage work chords through
the same winding leakage and optional reversible-core source. Core saturation
work is projected only through its original scalar flux constraint. Fixed-linkage
gap chords remain attractive through field polarity reversals. The no-core
case has an analytic endpoint action; saturation and nonlinear dielectric paths
retain the native bounded Newton probes. All four shared-gap corner states must
remain inside original source support; an unsupported corner refuses the trial.
These symmetric, consistent time maps retain the original support, circuit
constraint and transactional refusal boundaries and require time refinement.
Focused regressions check physical displacement/momentum kinematics, reciprocal
energy and source work, polarity and source-knee crossings, and the existing
independent RK4 refinement reference. These are numerical-law checks, not new
measured material or actuator validation claims.

## Annular DG terminal passages

Fixed DG terminal passages accept concentric axial annuli through
`DgTerminalChannelKind::Annular` and `dg_terminal_annular`. The original annular
momentum owner retains the actual local gas, unequal radial cell masses, both
no-slip walls, viscous history and separate pressure, kinetic, viscous, numerical
and material-change work. Admission requires positive ordered radii,
representable area, radial gap <= length/10, bounded radial cells and a resolved
speed Mach ceiling <= 0.3. Attached-area and combined-state budgets still apply.
Harmonic coupling retains every radial state. Finite-gas inventory includes the
actual annular passage volume, with main and branch gas ownership kept distinct.
This is a fixed concentric, axially developed Newtonian reduction. Eccentricity,
radial moving-pad openings, entrance jets and transition prediction remain
outside its scope. Existing terminal-group, gas-refresh, caloric-heat and
harmonic-field regressions now include annular geometry and independent ring
mass and logarithmic-shear calculations.

## Continuous simultaneous wall-film-gas windows

`WallFilmThermalEnvironment::advance_implicit_window` advances the coupled
enthalpy/volatile/gas endpoint over the complete requested duration. It shares
the original window clock, cumulative source apportionment and accepted-step
account with the split continuation. Wall internal heat and gas external heat
are full-window joules, assigned once by represented interval endpoints.
Accepted wall storage, liquid-film storage, gas internal energy, phase transfer,
sensible heat and radiation remain separate physical accounts. The implicit
gas storage already contains its sensible supply and is never counted twice.

`WallFilmImplicitWindowControls` bounds attempts and the complete window's
Newton iterations, residual evaluations and gas inversion work. Coarse or
rejected proposals cannot renew these budgets. Gas finalization reserves its
remaining original inversion bound; success refunds unused iterations, while
a failed finalizer retains the reservation when no partial receipt is available.
The report distinguishes known iterations from conservatively charged work.
Only an original phase donor `StepSizeLimit` may shorten and retry an interval;
source, caloric, nonlinear, balance and cancellation refusals remain terminal.

Every accepted step contributes to compensated whole-window balances. Final
admission independently checks actual wall storage, total and volatile mass,
carrier preservation, source allocation and accumulated absolute energy and
mass defects. Rejected proposals contribute no heat, mass, events or physical
time. A late failure publishes none of the accepted prefix. Exactly dry film
may form, evaporate completely and form again across successive windows;
original wall phase plateaus retain their enthalpy and phase-fraction storage.

This fixed-grid, first-order continuation inherits the endpoint model's
uniform-wall, pure-liquid, fixed-area and fixed-volume assumptions. Endpoint
rate exhaustion estimates still require physical timestep refinement.
Four focused groups cover manual endpoint composition and source allocation,
dry/wet/dry/wet continuation, native wall latent storage and late atomic
refusal; existing window and single-endpoint regression groups also pass.

### Second-order Cartesian diffusion with finite vapor

`FiniteMoistureAtmosphere::advance_cartesian_sdirk2` and
`advance_shared_cartesian_sdirk2` solve the actual finite vapor inventory at
both SDIRK stages. Shared specimens use one RH per stage, stable-ID ordering,
their original anisotropic/nonlinear material owners and one compensated sum
of local water changes. Stage two retains the ORIGINAL gas and solid storage
references and adds the explicitly weighted first-stage physical face fluxes;
stage-one uptake is never debited twice. Its first-stage transport receipt is
diagnostic. Endpoint face transfers, vapor, carrier species and gas energy
changes describe only the original-to-final physical interval.

Sorption and optional surface emission refresh at each trial RH. The accepted
stage envelope includes original field extrema and only equilibria connected
by nonzero first-stage or current second-stage physical face coefficients.
Both pre-source and final admission use that envelope. Inadmissible stages
refuse without clipping or committing a partial atmosphere/specimen state.
Every scalar-bracket trial must remain inside admitted source support; an
inadmissible extreme trial may require a shorter timestep.

The original prescribed-temperature, fixed-volume dilute-vapor scope remains.
There is no condensation, sorption heat, moving interface, global temporal
error certificate or experimental-material validation implied here. G1
manufactured coupled relaxation, nonlinear flux/source and shared-inventory
fixtures check the method; G0/G4 fixtures cover equilibrium, cancellation and
stiff-stage refusal. Existing backward-Euler entry points retain their method.

The finite-vapor Cartesian SDIRK path also provides
`advance_cartesian_sdirk2_adaptive_to` for a complete requested duration.
The existing controller checks moisture, RH and vapor mass independently with
second-order local step-doubling scaling, charges both implicit stages per
proposal, and sums only accepted physical half-interval RK transfers. Final
water, carrier and total-mass checks admit both inventories together. Failed
stages, rejected probes and a late work-limit failure publish no partial
duration. G1/G4 fixtures independently check the coupled relaxation, all three
error coordinates, stiff-stage retries and complete-stage work reservation.

## Adaptive time accuracy through wall and film latent storage

`WallFilmThermalEnvironment::advance_implicit_adaptive` compares one simultaneous
endpoint with two consecutive half endpoints from the same accepted state.
The original thermal first-order controller supplies the local-error correction,
including unequal represented halves, and bounded step growth/reduction.
Only the fine pair enters the shared physical window account and retained
history; all three solves and rejected comparisons spend the same global
Newton, residual, gas-inversion and attempt budgets.

The componentwise error estimate observes original wall enthalpy, liquid phase
fraction, wall and gas temperatures, actual film inventory, every catalog gas
species, and liquid-film specific enthalpy when both endpoints are wet.
Relative allowances scale changes from the interval's accepted start, so an
arbitrary caloric datum cannot loosen accuracy. A latent wall plateau can
therefore require refinement even when its temperature is exactly constant.
A dry film has no invented caloric coordinate or positive mass seed.

Uniform body and gas source totals are apportioned at absolute cumulative
window endpoints. Coarse and rejected sources never reach accepted balances.
Only an original phase donor step-limit refusal permits admission subdivision;
other source, nonlinear, conservation and cancellation errors remain terminal.
Represented admission spans never exceed the source's maximum. A bounded
comparison/reduction count and minimum span prevent endless accuracy retries.

The report exposes accepted/discarded work, actual fine samples, component
differences, maximum accepted normalized error and an explicit next-span
suggestion. This remains a first-order physical discretization with a local
error estimate, not a global trajectory or continuous dry-out-time bound.
Known changes of forcing require a new authored window. The endpoint's
uniform-wall, pure-liquid, fixed-area and fixed-volume limits remain.

Focused regressions compare tightened controls against independent continuous
wall/film/gas balance equations, isolate phase-fraction refinement at constant
wall temperature, and verify fine-only source composition, complete discarded
work charging, bounded rejection and atomic refusal.

### Reciprocal acoustic guide loading of spatial piezoelectric transducers

`harmonic::acoustic::PlaneWaveAcousticPort` prepares an original exterior-face
pressure/flow map. The exact P1 outward normal weights b give U=b^T v and
traction -b p. One original `FluidAcoustics` owner supplies the coherent
density, adiabatic sound speed and bulk modulus at the evaluated state.
The retained uniform guide mode has Z=rho*c/A and p=Z U+2 p_in under the
original exp(-i omega t) convention. The guide section can include a declared
rigid area; moving-face geometry remains the original tetrahedral trace.

`solve_harmonic_with_acoustic_ports` adds -i omega Z b b^T and -2 b p_in to
the same sparse electroelastic pencil, including intrinsic relaxation and
external electrical shunts. Upper-triangle entries are mirrored exactly.
Pressure recovery retains complete fluid loading in support reactions.
Independent incoming/outgoing wave powers and signed p/U work remain distinct
from conductor and ceramic heat. Net outward acoustic power enters the same
physical work acceptance gate; incident sound may supply a real sensor/shunt.

Ports retain their exact model, unique consistently outward planar exterior
faces, complete guide area and explicit surface-coordinate work bound.
Duplicate/foreign ports, inconsistent geometry, invalid original fluid states,
out-of-band wavenumbers, excessive pressure/flow amplitudes and cancellation
refuse before publishing a response. Material evaluation occurs once per port
before sparse assembly, never in Krylov iterations.

This is a declared single uniform-mode, lossless semi-infinite guide model.
The caller owns its transverse-mode cutoff/accuracy domain and rigid remainder.
It is not free-space diffraction, a higher-mode/evanescent approximation,
viscothermal guide loss, a nonlinear acoustic model or empirical validation.
Four synthetic G1/G3/G4 regression definitions specify electrical actuation,
incoming-sound sensing with shunt heat, rotated pressure traces, rigid reflection
and refusal. Rust syntax/format and source/equation review passed; independent
scalar work calculations agreed. Native compile/tests remain unverified under
the editing environment's disk and memory limits.

## Source-bound saturating core, gap and leakage inductors

`material_circuit::inductance::MaterialReversibleInductor` maps an admitted
reversible H(B) source into the native circuit's bounded I(lambda) storage.
The declared uniform-flux reduction is Phi=A*B,
N*I=l*H(B)+R_gap*Phi, lambda=N*Phi+L_leak*I. Original source breakpoints are
transformed exactly, with both polarities of the isotropic zero-remanence law.
No sampled fit, constant apparent inductance or second time integrator replaces
that relation. Source snapshot, query context and fixed geometry stay inspectable.

Its stored energy is the same potential in both descriptions:
A*l*integral H dB + R_gap*Phi^2/2 + L_leak*I^2/2.
Operating points expose core B/H, current, differential inductance and separate
core/gap/leakage energies. These partition the native storage; they are never
additional energy debits. B is recovered from the retained transformed segments,
avoiding cancellation of lambda-L_leak*I when leakage dominates linkage.
The native right-linkage derivative owns source kinks, except at the final
endpoint where the left derivative applies.

`ReversibleInductorCircuit` binds these sources to original inductor branches,
sets their zero-field differential coefficients, and delegates complete
Kirchhoff assembly, initial consistency, analytic actions and audited stepping
to the existing `CircuitDae`. Nonlinear capacitors, reciprocal capacitance
blocks, other authored scalar inductors and disjoint linear mutual links keep
their native owners. No extra accepted history or time integrator is introduced.

The source must be a field-only reversible potential. Explicit thermal scaling
or additional response axes refuse, because freezing a thermal free-energy
snapshot cannot supply its magnetic entropy account. Source-domain exits,
duplicate ownership, scalar nonlinear/mutual overlap and unrepresentable
geometry refuse. Gap reluctance and leakage are explicit fixed inputs; this
reduction does not infer fringing, distributed field resolution, hysteresis,
eddy loss or temperature evolution. Graph resistance retains physical heat.

Focused regressions cover every source knot and polarity, independent
core/gap/leakage work, turns scaling, leakage-dominated observations, a native
saturating RL trajectory and physical/source refusals. Independent arithmetic
checks over 10,000 varied geometries agree with the source potential within
1.5e-15 relative error. Rustfmt parsed the implementation and tests; complete
source-to-native circuit tests remain unexecuted in the dependency environment.


### Prepared source-bound piezoelectric circuits

`material_circuit::piezoelectric::prepared::PreparedPiezoelectricCircuit` caches the original fixed-material index-one Kirchhoff system for a fixed positive clock. Native algebraic elimination retains every differential flux, charge and mechanical coordinate. Native exact quadratic propagation and independently prepared physical source observations preserve resistor feedthrough, separate electrical/mechanical supply and per-resistor heat. Original endpoint Kirchhoff rows, source-domain limits, material conductor power and total stored-energy balance are checked before transactional publication. Warm stepping performs no owner allocation, factorization, exponential or Newton solve. Singular algebraic blocks refuse explicitly without reducing conserved state; ordinary descriptor stepping remains available. This is the linear held-input model up to floating-point error, with a conservative numerical whole-window domain bound, not measured material validation or a hard real-time guarantee.


### Prepared fixed-voltage piezoelectric motion

`piezoelectric::prepared::PreparedPiezoelectricVoltage` owns the original rod and complete accepted extension, momentum and free charge under the terminal voltage inherited at preparation. Native exact quadratic propagation eliminates only the algebraic voltage constraint, retaining the short-circuit stiffness, consistent inertia and reciprocal force theta*V. Each candidate reconstructs charge from the original preparation anchor and independently checks the original full stored energy, force work, voltage-source work, induced charge and constitutive voltage before publication. Native whole-window mechanical norm admission bounds strain, and the fixed voltage with its explicit tolerance bounds electric field. Warm stepping uses fixed arrays and performs no owner allocation, factorization, exponential or Newton solve. Refusal and cancellation preserve the complete accepted body. This is the held-input linear model up to floating-point error with conservative numerical domain admission; changing to an incompatible ideal voltage requires another physical boundary model and cannot silently reset charge.

## Regional coexistence in the conserved finite water cell

`EquilibriumWaterCell::from_isochoric_state` admits an original prepared
`IsochoricWaterState` directly. Its mass, fixed specific volume, raw internal
energy, source mode and dense-root controls remain the source of the initial
M,U,V inventory. Initialization requires the exact bounded source controls
after applying the cell's absolute budgets per unit mass.

The original stream, finite-wall heat and sealed moving-volume paths now
preserve that source mode when rebinding V/M. Regional coexistence can therefore
evolve liquid/vapor masses, pressure and density above 623.15 K through the
same conservative ledgers and reciprocal pressure-work solve. Phase volumes
use the retained original regional endpoints; mass quality is never substituted
for vapor volume fraction. The constructor does not create a second EOS or
replace conserved energy with an approximate recovered source image.

Both native water-cell acoustic observers support the regional pair and use
exactly the retained thermal source density, with the selected original bulk
modulus. They expose the signed arithmetic density difference separately.
The equilibrium observer preserves the original Gibbs accessor and adds a
regional tangent accessor retaining its actual root/pressure defects. These
are explicit asymptotic acoustic closures; the observer does not add phase
kinetics or advance the thermal inventory.

Regional operation remains a bounded coexistence window below the critical
temperature, on one side of the 623.15 K source-fit seam. Losing either phase,
crossing that seam, or exceeding the original source/residual budgets refuses
the proposed interval. No noncondensable gas, slip, nucleation, moving interface
or thermal-gradient model is inferred. Stream transfer keeps its existing
held-donor approximation; moving-volume work keeps its existing trapezoidal
time accuracy. Coupled material and source errors remain distinct from
conserved mass/energy residuals.

Five added native integration tests cover 640 K heat/phase/acoustic updates,
two-cell source-retaining transfer, reciprocal isothermal phase expansion and
return using raw Delta U plus p Delta V, source-window/late-cancellation
rollback, and finite conducting-wall exchange on the same regional inventory.

All 23 exact native water-cell tests passed, including the five additions and
18 existing regressions. The source closure uses original feedback helpers and
genuine native material, duct and thermal dependencies; no whole-workspace
Cargo/DSR pass is implied.

## Finite-rate water compression in reciprocal cavity networks

`cavity_phs::relaxation::RelaxingWaterCavity` retains the original
`RelaxingRegionalWater` and `HelmholtzCavity` geometry. Native neck inertia uses
the exact thermal density and declared end corrections. Original equilibrium
bulk storage and one positive compression arm enter the existing
`RelaxationBranch` and `PreparedExactZoh` owners. Pressure and volume flow remain
reciprocal, including through native mechanical transformers.

Fixed-volume net inflow gives first-order excess mass and density. Pressure,
neck flow, storage and instantaneous compression/neck losses come from the
accepted native state. Relaxed or explicit frozen memory is supported. Native
prepared stepping retains integrated input work and dissipation with reusable
warm scratch.

This homogeneous small-signal lumped model retains caller-established
relaxation kinetics, explicit neck resistance and end geometry. It does not
evolve nonlinear phase inventory or temperature or infer exterior radiation
or viscosity. Four focused native regressions pass: original 100 kPa/20 MPa
density and inertia; independent refined mass/momentum/memory ODE plus
work/heat; actual mechanical-transformer power and total energy; native
admission and tiny density perturbations. Exact-source tests used genuine
material and pHS owners; no full workspace Cargo/DSR pass is claimed.

## Finite-body thermal feedback through spherical Hertz spin events

`friction_heat::sphere_plane::propose_sphere_plane_hertz_heat` composes the
original sphere Hertz support trajectory with the original finite enthalpy
network. Actual spin arrest, static retention, breakaway, motor reversal and
attitude remain with the native mechanical owner. Every coupled trial starts
from the same opaque accepted motion. Its solved common face temperature
refreshes the original pinned static/kinetic claims once; that snapshot retains
the actual native normal reaction and separately declared footprint radius.

Only integrated sustained physical friction loss enters the existing massless
face. Both finite bodies retain their own temperatures and conductivities;
signed body heat includes ordinary hot-to-cold exchange. Native motor work,
external thermal sources and the independently bounded static constraint-work
residual remain separate. The latter is never thermalized. The combined gate
reconstructs actual compensated native kinetic change plus caloric change,
less outside heat, mechanical supply and the explicit static residual. Motion,
thermal cells, original source receipt and material callback payload leave the
owner together after all original and combined gates pass.

This is a stationary circular contact-column, area-averaged thermal model.
The original thermal face must match the actual Hertz centre, normal and radius;
ordered thermal body IDs and the sphere-A/plane-B side geometry are checked.
The exact pure normal-axis spin subspace admits no COM or plane motion,
tangential drive, impact heat or displaced footprint. No velocity, force or
momentum is projected to meet these limits. Original material/caloric identity
and the applicability of the contact-column approximation remain with their
material owner. Radial flash gradients, material transport through a spinning
footprint, normal compliance, thermal expansion and whole-sphere thermal
geometry are not inferred. Endpoint coefficients are held over each native
window; backward-Euler thermal transport and thermal sampling require temporal
refinement independently of exact held-coefficient stopping.

Five focused G1/G3/G4 native tests pass in
`tests/sphere_plane_friction_heat.rs`: independent analytic spin/attitude and
two-body conduction, first-order convergence against an independently integrated
continuous thermal-fade ODE, heat-induced static breakaway with motor work,
opposed-motor reversal and signed-axis covariance, and late source/cancellation
rollback with exact retry. The 8/16/32-window maximum physical-state errors are
0.0531229, 0.0268324 and 0.0134851. These verify the bounded model with synthetic
source cards, not measured calibration, full-workspace Cargo/DSR or throughput.

## Motor torque in source-derived directional sphere contact

`friction_heat::directional::sphere_plane::advance_source_sphere_plane_with_torque`
threads a prescribed constant world couple through the original native sphere
trajectory and its original source-resolved support and impact laws. Actual
torque-driven spin reaches preimpact speed queries; supported motor rolling,
breakaway and sliding reach the motion-consistent speed/path source closure.
The full static and kinetic material claims, original selection receipts and
plane-owned directional frame remain authoritative. The existing entry
delegates zero couple.

The returned native motion receipt separates signed motor work and angular
impulse from centre-force work, moving-plane work, sustained slip heat and
unallocated impact loss. Query callbacks see the latest accepted private
mechanics/history. A late source refusal or cancellation discards the entire
candidate, including earlier motor motion and collision source updates.

Four G1/G4 regressions in
`tests/rigid_directional_heat/sphere_plane/impacts/torque.rs` cover sourced
static rolling with zero path/heat, motor breakaway with actual accumulated
slip, torque-dependent preimpact source coefficients and bounce work, and late
source/cancellation rollback. The fixture ledgers were independently derived
and all changed Rust sources parse with native rustfmt. These tests have not
been executed against the current fs-mbd/fs-couple dependency closure; no
Cargo/DSR gate is claimed. This remains an isotropic rigid point-contact model
with a held world couple, explicit temperature/pressure owners, and no inferred
contact area or impact-to-heat allocation.

## Adaptive second-order diffusion with shared finite vapor

`FiniteMoistureAtmosphere::advance_shared_cartesian_sdirk2_adaptive_to`
extends the same second-order controller to separate native Cartesian specimens
sharing one vapor inventory. Each field retains independent local error and
six-face water accounts, including simultaneous uptake/desorption with zero net
reservoir change; RH and vapor mass are checked separately. Only accepted
RK-weighted half intervals contribute transfers. Both implicit stages count per
coupled probe, stable IDs determine output order, and final whole-window water,
carrier and total-mass checks admit all endpoints together. The local error
estimate remains asymptotic, without a global-error guarantee.

G1/G4 regression definitions cover coupled mean/difference modes, stationary
vapor with opposing transfers, ID permutations and late stage-budget failure.

## Coupled orifice law on the first N open tone holes (`open_hole_orifice N`)

The optional reed-card record `open_hole_orifice N` (N = 1..=16, `fs_duct::exterior::MAX_ORIFICE_PORTS`; after `constriction_noise`, before `compile_limits`; once
only) solves the quasi-steady sharp-edged orifice law `delta_p_i = K_i |U_i| U_i` (`K = rho / (2 C_vc^2 A_h^2)`, the
first-hole shunt's own law and coefficient) of the first N open tone holes TOGETHER instead of the first alone
(`RegionalDuct::with_open_hole_orifice_ports`; `RegisterHoleOrificePorts` in `regional/outlet/orifice_ports.rs`). The
holes couple through the bore, so the per-sample solve is an N-vector one, a damped Newton on `U = u_th - A Phi(U)` with
`A = R^-1 (I - D)(I + D)^-1`, `D` the direct taps of the N x N reflection kernels `S = (I + R Y)^-1 (I - R Y)` and `R`
the holes' numerical port references (the ordering `R Y`, not `Y R`, matters as soon as the references differ). `Y` is
the admittance of the linear TMM chain seen at the N orifice terminals
(`fs_duct::exterior::input_impedance_holes_and_ports`: a nodal ladder of the blocks' own admittances, a pure-shunt
block such as a sealed hole without series length sitting on its node; every other hole stays linear), and it equals the
inverse of the one-port network impedance for N = 1. The reed receives the sum of the N `s12` kernels applied to the
accepted drops; radiated hole sound still comes from the linear network, as for the first-hole shunt.

An absent row or `N = 1` runs the unchanged first-hole shunt, bit for bit (`RegisterHoleShunt::Single`), and an N above
the number of open holes uses all of them. For N >= 2 the realized matrix is made passive by bounding the matrix
itself: the direct matrix is clamped to power-wave norm 0.999 and the whole tail scaled by the largest factor (found by
bisection, 2^-14 resolution) for which `max_w || D + lambda T(w) ||` stays at most 0.999 on the 2n grid. The
first-hole line's own rule, `|D| + max_w |T(w)| <= 0.999`, ignores the phase between `D` and `T`: it costs 0.1 % of the
tail for one port but 0.4 % for two and 3.3 % for three on the test bore; with it the coupled flows departed from the TMM
chain by up to 15 % at 100 Hz on a 4096 grid, with the bisected rule by at most 3.4 %. That legacy rule is kept (`PassivityRule::Triangle`) only so the N = 1 equivalence
test is exact; N = 1 in the card never reaches it. The row takes a constant-medium run or a `gas_field continuous`
binding (`GasFieldDuct::response_with_wall_holes_and_ports`: each orifice reads its own local gas, K and numerical
reference included). A molecular bore, a midpoint `gas_field` (N > 1) and a segment-gas duct with several regions
refuse by name. On a continuous field with `open_hole_interaction` there is ONE approximation, stated on the page too:
each orifice sees its own UNCOUPLED radiation load (`Y` is the uncoupled chain's) while the radiated sound uses the
coupled flows; the fully coupled orifice network is not built. Cost: `N^2 + 2N` kernels of the grid length instead of 3
(15 at N = 3, 143 at N = 11, 288 at N = 16), and one more frequency sweep at build; the passivity step holds one `NP x NP` complex matrix
per frequency bin (4.1 KB at the cap of 16, about 67 MB on a 8192-point grid, about 540 MB on the largest 65536 grid; 2.3 KB at the cap of 12 that
preceded it). The cap is a budget (stack arrays and that matrix), not a validity limit. It was 3 until the oboe D5 page (11 open holes), then 12 until the
clarinet G4 (14 open holes) and Gb4 (13); the solver residual test runs n up to the cap
(1200 random passive systems) and `compose::port_tests` checks an eight-hole and a fourteen-hole ladder (reciprocal, passive, leading block
independent of the port count). Measured on the oboe D5 card (11 open holes, `gas_field continuous`, tmmi 2, 96 kHz, 2.65 s):
`open_hole_orifice 3` renders the posted wav byte for byte on the build with the cap at 12, and `open_hole_orifice 11` (all
open holes) changes the steady tone by +0.8 c, the level by 0.00 dB and the harmonic levels by none to the digit, with the
grader count unchanged (6 of 14); user CPU 591 s at N = 3 and 1118 s at N = 11, so the cost is far below proportional to the
`N^2 + 2N` kernel count (15 to 143). Measured on the clarinet pages (24-hole Buffet bore, `gas_field continuous`, tmmi 2, 48 kHz, 2.65 s; 13 pages
with N from 2 to 12 equal to their own open holes, the controls byte-identical to the posted pages): user CPU 4.9 to 5.6 times the control's at every N
(3330 s at N = 3 to 4158 s at N = 12, controls about 700 s; Bb4 at N = 2 took 6609 s against a control of 1340 s), peak RSS 61 to 208 MB; the cost does not
scale with the audio length (A5 and Bb5 rendered for 0.9 s took 3551 to 3693 s at N = 2 to 6, the same as for 2.65 s) or with N, so it is a fixed cost of the render's build, not of the per-sample solve;
the steady tone moved by 0.5 to 4.4 c and the level by 0.2 to 5.5 dB on eleven pages, and on A5 (N = 10) and Bb5 (N = 11) the card played a low tone, 325 Hz 8 dB over A5's own line and 339 Hz with no line left at
Bb5 (the cause is not isolated; the family card records the same kind of tone for F6 without this row). Not measured: an N above 12 on a real card, and an N above 3 on any other family.

Evidence (G0/G1): `fs-duct` `compose::port_tests` (the ladder's admittance is the inverse of the one-port network
impedance to 1e-9, reciprocal, passive, independent of the port count, and the sign of its off-diagonal term fixed by a
coincident-pair test); `orifice_ports::tests` (the solver's residual over 1200 random passive systems, n = 1 up to the cap, and a stiff drive;
the realized matrix is passive and is the admittance scattering, not its transpose; the time-domain steady state with
linear drops equals the network built from the realized kernels to 1e-14 and the TMM chain to a few percent at a 4096
grid; the general solver run as one port follows the first-hole line through a strongly nonlinear drive to 1e-6);
card level `reed::orifice_tests` (an absent row and `open_hole_orifice 1` are bit-identical in audio and reed state, a
one-hole duct ignores the count, 2 and 3 differ from 1 and from each other, bad counts refuse, a fourteen-hole card takes fourteen ports); on a continuous field
`reed::orifice_tests::on_a_field` (an absent row and `open_hole_orifice 1` are bit-identical, 2 and 3 act and differ,
3 composes with `open_hole_interaction tmmi` and the printed lengths, a midpoint `gas_field` refuses by name) and
`fs-duct` `profile::gas_field::ports_tests` (uniform field = homogeneous ports bit for bit, one port = inverse of the
network impedance on a graded field, three ports reciprocal, passive and independent of the port count, the interaction
leaves `Y` alone and changes the radiated flows).

## Second-order diffusion coupled to finite solution activity

`FiniteMoistureSolution::advance_cartesian_sdirk2` and
`advance_shared_cartesian_sdirk2` solve native Cartesian moisture fields and
the same finite solvent inventory at both implicit stages. Actual water uptake
concentrates the solution; desorption dilutes it. Every trial refreshes solution
activity, original sorption and surface transfer, and nonlinear principal
diffusivity. Original-to-final compensated water increments close the reservoir
once. Stage-one fluxes enter the second-stage load exactly once; the retained
first-stage receipts are diagnostic, never additional solvent debits.

Dissolved ideal/nonideal solutions, neutral equilibrium solubility and ionic
equilibrium keep their original thermodynamic and total-component owners.
Second-order accuracy assumes a smooth constitutive/phase branch. Kinetic
solubility is explicitly refused because its endpoint-frozen phase update does
not supply the required coupled RK equations. Equal prescribed temperature and
pressure, negligible vapor storage and the original source domain remain in
force; no liquid-immersion permeability, ion transport, thermal evolution or
mechanical feedback is inferred. Shared participants retain unique stable IDs,
native geometry, source laws and aggregate cell admission. Borrowed original
fields and loads avoid cloning either for each scalar reservoir trial.

`advance_cartesian_sdirk2_adaptive_to` completes one requested duration with
the original temporal controller. Native content, solvent mass and activity
receive independent error scales; both stages spend work even in rejected
probes. Only accepted physical half-interval transfers enter the face ledger.
One fixed whole-duration water budget admits both endpoints together; a failed
duration returns no partial bath or field.

G1/G4 regression definitions cover independently integrated finite-bath
relaxation/refinement, equilibrium phase buffering and solute conservation,
opposing shared specimen transfers, kinetic/cancellation refusals, and adaptive
whole-duration accuracy and work limits.

## Condensate films on the original spatial thermal graph

`finite_gas_thermal::wall_film::deposition::advance_deposition_with_wall_film`
attaches one independent Dry/Wet film to an existing native thermal cell.
The original `BackwardEulerStep` supplies every cell's caloric storage, contact
flux, source and prescribed bath. The original material wall-film endpoint
supplies evaporation/condensation, zero-inventory complementarity, liquid storage,
vapor-carried energy and finite-gas feedback. A shared Newton solve uses native
cell-energy increments, one volatile inventory coordinate and gas temperature.

The selected patch occupies actual exposed wall area. Its gas convection is
reciprocal; the remaining area retains the original ambient cooling. Native
reservoir coverage and contacts remain authoritative. An ambiguous callback
replacement of that occupied ambient channel refuses. Other endpoint sources,
link conductances, internal transfers and bath replacements remain available,
and their callback sees the actual trial wall graph and gas.

The underlying wall keeps its mass, material identity and original phase chart.
Liquid-film h(T) comes from its separate source; melting/freezing storage in the
wall is never treated as film mass or apparent heat capacity. A conductive
neighbor may supply evaporation heat or absorb condensation heat while the
wall stays on its latent plateau.

No gas or wall caloric root is nested inside residual evaluations. After one
original gas reconstruction, a fresh native heat/phase evaluation checks the
actual joint endpoint. The independent external-heat audit runs only on this
final replay, reusing the same resolved coefficients without another callback.
Local cell defects, the complete signed/absolute energy account, source-port
roundoff and original species/mass balances all gate publication. The shared
implicit work meter also supports enclosing continuations.

This event-free first-order endpoint preserves native exact-target event
semantics and fixed geometry between events. Each cell remains a uniform
control volume requiring spatial/time refinement; the film is fully wetting,
pure liquid and thermally equilibrated with its wall patch. It adds no drainage,
moving contact line, boiling or liquid-mixture activity model. Focused tests
cover analytic conductive supply on a wall fusion plateau, dry condensation,
the uniform-wall limit, manufactured partial-area/source feedback and atomic
physical/source/event/work refusal.

## Native thermal events with a continuing condensate film

`deposition::timeline::advance_deposition_with_wall_film_timeline` advances
the spatial wall-film endpoint through the original thermal schedule. Births,
contacts, separations and boundary changes apply in authored equal-time order,
including events at the initial and requested final time. Every heat/phase
interval stops before a native event; the event pass permits zero thermal
substeps, so it cannot insert an independent sensible-only heat advance.

The initially active wall keeps the same physical wetting patch and separate
film inventory. New contacts must leave enough native exposed area for that
patch; incompatible coverage refuses the whole candidate instead of hiding,
draining or discarding condensate. Unselected cells may be born during the
window and subsequently exchange heat through original contacts. Birth mass
and energy retain each new cell's original material datum. Final wall storage
subtracts each cell's initial or birth enthalpy before multiplication and
summation, so incoming formation energy cannot masquerade as latent heat.

External gas heat is a full-window uniform source, apportioned by cumulative
physical endpoints. Wall sources retain their native watts and boundary events.
Only original phase step-limit refusals shorten a trial; other source, nonlinear
and conservation failures propagate. All proposals share the original implicit
work meter, the complete-call attempt limit and the native thermal substep
allowance. Accepted intervals alone contribute physical accounts. A final
source, coverage, conservation or cancellation failure returns no partial graph,
gas, film or event history.

Compensated totals retain accepted phase, sensible, storage and external-source
accounts. Final admission checks actual native wall storage against interval
sums, exact total gas-source allocation, film/volatile/whole-gas mass, carrier
preservation, and cumulative absolute mass/energy defects. The final payload
belongs to the last accepted heat endpoint before any final-time events.

This bounded continuation remains first-order with a fixed maximum interval;
native event clipping and donor admission do not estimate temporal accuracy.
The three focused groups compare manual event-cut composition, isolate a
large independent birth datum from actual conductive/latent heat, and exercise
late physical/source/work/cancellation refusals.

## Adaptive second-order diffusion with shared finite solutions

`FiniteMoistureSolution::advance_shared_cartesian_sdirk2_adaptive_to` completes
a requested duration through the existing shared temporal controller. Every
original specimen, solution activity and solvent inventory receives separate
local error control; opposing transfers cannot cancel specimen error estimates.
Complete-window accounts include only accepted RK-weighted transfers, preserve
the native nonvolatile stocks, and charge both implicit stages of every proposal.
All fields and the one finite reservoir are admitted together after water
conservation checks. Native kinetic-phase refusal remains; temporal estimates
assume smooth source and algebraic phase branches and certify neither global
error nor phase-crossing order.

G1/G4 regression definitions cover opposing fields with a stationary pure bath,
finite concentration/activity feedback against an independent continuous
reference, accepted signed transfers and late stage-work exhaustion.

## Adaptive exponential binary diffusion in the live material assembly

`MaterialExchangeState::adaptive_binary_diffusion_step_candidate` and
`advance_with_adaptive_binary_diffusion` connect the original finite binary
gas network to the existing reaction, thermal, electrical and piston owners.
The gas stage refreshes actual source states and partial enthalpies, resolves
local species/temperature error, and enforces separate trial/pair-work limits.
Every accepted gross transfer marks both chambers' closed two-port audits open,
including counterdiffusion with zero net mass.

The original assembly balance encloses the pre-diffusion inventories. Actual
diffused species feed the original reaction owner once; transported and reaction
enthalpy remain in their gas inventories and are not redeposited as external
heat. The common feedback and whole-window mass/energy gates must pass before
gas, thermal field, charge, motion and clock publish together. Source failures,
pair-work exhaustion, late feedback refusal and cancellation preserve the
entire accepted assembly.

The complete diffusion/reaction/heat-motion split remains first order, and the
gas-only local estimate does not certify coupled splitting error. This path
retains the original binary molar-frame and fixed-volume transport substage;
pressure/thermal diffusion, resolved nozzle flow and a simultaneous
moving-volume constitutive solve require their existing distinct owners.
G1/G4 regression definitions cover an independent stiff exponential limit,
original-owner composition with the full feedback solve, reaction consumption,
zero-net-mass open audits and transactional late failures.


### Continuous prepared piezoelectric force/current tones

`PreparedPiezoelectricTones` advances continuous fixed-frequency force/current programs through the original three-state source-bound rod and native extended quadratic-flow owner. One shared numerical cosine/sine carrier per frequency retains both ports' independent phase offsets without introducing physical storage. Native signed quadratic observations integrate mechanical and electrical work independently; native linear observations integrate current and terminal voltage. The actual accepted carrier norms enter a conservative whole-window strain/field bound, including numerical phase-norm drift. Body, phase, accepted clock and reports publish atomically after original source-domain, charge, energy and cancellation checks. Warm stepping allocates no owner buffers and performs no exponential, factorization or Newton solve. Scope is the fixed-material linear sinusoidal/DC source model up to floating-point error, without phase normalization, energy rescaling, empirical material validation or a real-time deadline claim. Six focused regression definitions cover resonance, weak coupling, multitone interference, independent source work, DC parity, cancellation and whole-window/budget refusals; native execution remains a separate verification obligation.


### Transient piezoelectric acoustic guide coupling

The original spatial midpoint solve accepts real incident/outgoing guide waves through the same `PlaneWaveAcousticPort` pressure/flow trace used by harmonic acoustics. The actual immutable `FluidAcoustics` owner supplies density and adiabatic sound speed. Reciprocal acoustic traction, consistent inertia, charge constraints, optional Maxwell/Debye ceramic memories and physical resistive harvesting share one native solve; prescribed support motion pays its complete reaction work. Acoustic outflow is audited against independently recovered incoming/outgoing wave energies and is counted once, separately from ceramic and resistor heat. Time-domain midpoint energies use no peak-phasor factor of one half. The accepted endpoint retains constitutive forces without spurious acoustic history. Foreign/overlapping traces, unrepresentable waves, amplitude/domain violations and cancellation refuse transactionally. This is a fixed-geometry uniform-guide model with caller-declared waveform/interface bandwidth; a timestep is not a spectral certificate. It does not imply free-field radiation, diffraction, a spatial acoustic volume, or time-discretization convergence. Three focused synthetic regression definitions cover sound-to-electric conversion with/without intrinsic memories, moving/rotated interfaces, reflection and domain/cancellation behavior; native execution remains a separate verification obligation.


### Finite-temperature full-tensor spatial piezoelectric dynamics

`spatial::transient::load::thermal::tensor::ThermalTensorPiezoelectricBody` couples the full native 3-D spatial fields and consistent-mass motion to one uniform finite ceramic temperature, with original per-cell tensor orientation and complete source-bound caloric inventory. Temperature proposals refresh the complete affine stress-charge C/e/epsilon source tensors and reassemble the same native operators under original budgets. The model may own refreshed tensors while borrowing unchanged geometry; identical-temperature trials reuse accepted operators. Each changed-temperature trial first solves electrostatics at fixed original displacement/free charge, restores physical momentum/history, and then runs the original shunt/current/support-motion midpoint solve. Signed reversible field/caloric transfer is derived from the full Helmholtz and entropy endpoints; actual resistor Joule heat is deposited only through its original receiver. Independent acceptance checks total field internal energy, kinetic energy, all finite caloric energy and actual external work/heat. Complete source/context/density/inventory checks, declared strain/electric-displacement bounds and a whole-interval caloric-plus-field capacity floor precede coupled Newton publication. This is fixed-reference geometry with one uniform ceramic temperature, fixed zero-voltage grounds and charge/current electrodes; intrinsic memory, nonzero ideal-voltage constraints, expansion, pyroelectric offsets, phase transitions and independent cell temperatures require additional accounted models. Endpoint-temperature coefficients make the thermal feedback generally first order in time. Four focused synthetic regression definitions cover moving/rotated supports, reversible caloric changes, heating-dependent capacitance at conserved charge across mesh refinements, separate finite shunt heating and transactional refusals; native execution remains a separate verification obligation.

## Moving stress histories on the pressure/PCM path

`pcm_wav::observation::quadrupole::moving::MovingQuadrupoleRenderer` reads
original medium-frame Q/Qdot/Qddot at each source's retarded time and sums
native moving quadrupole pressure before audio conversion. Its existing
`PressureRenderer` interface feeds the original decimator and WAV writer.
Complete-window admission covers source support, continuous closest approach,
maximum Doppler clock rate, finite work and scratch budgets. Received
bandwidth must cover Doppler and geometric amplitude modulation below Nyquist;
a declaration does not establish source interpolation or spectral accuracy.

Callbacks allocate nothing, preserve global integer sample clocks, use DD
source summation and publish output/cursor only after complete success.
Out-of-support queries, failed reads and cancellation preserve accepted state.
Scope includes uniform and admitted quintic source translations, a fixed
omnidirectional receiver and homogeneous quiescent lossless propagation on the original moving operator's
f64 clock domain. No source evolution, interpolation, particle velocity or
directional/reflected moving field is inferred.
`tests/moving_quadrupole_pcm.rs` covers history/bandwidth/clearance refusal,
coherent cancellation, block invariance, atomic failure and native WAV parity.

## Weak directional modal damping

`ModalPlateViscosity::from_psd_modal_matrix` retains positive dissipative
directions resolved by its compensated arithmetic even beside a much stronger
direction and an exact null. Its pivoted factor uses the existing double-double
owner for LDL Schur updates, with locally propagated arithmetic uncertainty.
Square roots and the arbitrary reference resistance enter only the realized
channels. Larger coupling ratios are formed before smaller cross terms to
avoid losing representable updates across widely separated scales. The former
`1e-10 * max_diagonal` positive-loss cutoff is absent; negative-pivot admission
and full reconstruction checks retain their existing boundaries. The ordinary
strictly positive-definite Cholesky path is unchanged.

This matrix-only API cannot identify positive source/projection error as an
exact physical null. Producers remain responsible for preserving their known
null spaces, and numerical Ritz vectors are checked by their actual physical
power rather than an assumed exact channel count. Focused native regressions
retain `1e-14` and `1e-100` loss ratios, correlated weak directions, dependent
rows, dense basis covariance, extreme-scale null motion and finite decay/heat
over a weak channel's physical relaxation time. These checks are ordinary
floating-point validation, not an interval or source-applicability certificate.

### Original stiffness in hereditary factorization

The shared regional relaxation factor supplies its original dimensional
stiffness entries to the compensated PSD fallback. Diagonal normalization is
confined to the established strictly positive-definite Cholesky path; a failed
Cholesky must not make rounded normalized entries the material source. This
preserves resolved weak strain energy and exact dependencies for regional
plate memory, shell substrate relaxation and coating relaxation, with the
existing full-form reconstruction and cancellation checks. Successful
positive-definite channel arithmetic remains unchanged.

## CST/DKT shell material damping

`ModalPlateViscosity::from_cst_dkt_shell` consumes the original CST/DKT shell
mesh, sparse mass pencil, reduced basis and one `ShellViscosity` per facet.
The spatial owner projects the passive membrane/bending resultant tensor,
including signed membrane-bending coupling, before the existing common owner
factors every retained modal cross term into dashpot channels. Physical strain
rates exclude numerical drilling stabilization. The supplied mass-normalized
basis is checked against the complete original mass matrix, including its
consistent and translation/director terms; no mode, frequency or normalization
is changed. Work admission covers projection, mass products, Gram checks and
factorization, with cancellation before publication.

The same channels enter `CoupledModalSystem` and its separate physical-loss
accounts. They introduce no spring, numerical damping, thermal destination or
new integrator. Source validity, frequency-band applicability, modal
completeness and use of the same original mesh remain caller obligations.
MITC and tied-shell kinematics retain their existing owners. The focused tests
compare sparse spatial power, the full modal quadratic form and channel power;
exercise lumped and consistent inertia, basis permutations, zero loss and
refusals; and compare coupled decay with an independently refined RK4 solution
of the original full modal equation while checking physical heat closure.

### Source-bound homogeneous shell viscosity

`plate_viscosity::shell::material::ShellMaterialViscosity` retains an immutable
`PlaneStressViscosity`, its original source state and plane, physical thickness
and orientation in the original facet. The existing material tensor owner
rotates the complete engineering-shear law and integrates homogeneous
resultants A = h V and D = h³ V / 12; B remains exactly zero. The preparation
refuses unrepresentable resultants and does not alter stiffness or inertia.
Material records are shared through `Arc` when multiple facets use one source.

`ModalShellMaterialViscosity::project` checks every retained eigenpair's
frequency against the intersection of all facet source bands, including zero
loss sources, then uses the original mass-normalized full CST/DKT projection.
The result retains facet sources, frequencies and full modal channels for the
existing physical-loss runtime. Cancellation and work refusal publish no
partial image. The source query point alone never establishes applicability;
bare mesh/model identity, homogeneous thickness correspondence, measured source
authority and modal completeness remain caller obligations. This adapter does
not infer layered, three-dimensional or hereditary material laws.

Focused tests compare rotated engineering-strain power independently, verify
linear/cubic thickness scaling, retain original records through modal
projection, and reject unsupported actual frequencies even for exact zero
viscosity. These are synthetic constitutive checks, not a material calibration.

## Material tensile admission on nonlinear acoustic strings

`NonlinearMaterialString::with_tensile_limit` explicitly binds one named static
strength claim per original homogeneous or perfectly bonded concentric layer.
The shared `AxialTensileCapacity` retains every source receipt, modulus and
annular geometry. No default strength, safety factor or allowable is inferred.
Missing, nonstatic or invalid strength refuses; measured whole-string and
helical winding sources require an explicit constituent stress model. Prony
sections and homogeneous hereditary bending sources refuse until their actual
layer memory stresses are included.

The existing Kirchhoff--Carrier storage supplies the actual common tension,
including geometric stretching from both transverse polarizations. Its exact
modal normalization supplies a continuous-span curvature envelope
`K = sum_n k_n^2 sqrt(2/(mu L)) hypot(q_xn, q_yn)`. Each layer is assessed using
`sigma_i <= E_i (T/EA + r_outer_i K)`, with the first maximum-utilization layer
governing. A single sine mode attains its curvature envelope; multimode
cancellation can make the envelope conservative. The two-plane norm keeps
this criterion independent of the transverse coordinate frame. The quantity
is a floating-point elastic model estimate, not an interval certificate.

Selected limits enter the original initial-state and candidate acceptance
checks. A crossed envelope refuses before physical state, material histories,
energy accounts or accepted time publish. Forces and stored energy remain the
original laws. Successful checks allocate nothing; preparation retains the
bounded modal factors and original capacity. Native reference, prepared and
pressure-rendering paths share the same check. Transversely moving supported
strings retain the guard through their original kinetic transformation: the
affine endpoint adds its actual stretching but has zero curvature. The
longitudinal moving-support extension below includes the original owner's
nonuniform axial stress. Transfer of a guarded potential
to the bow bank's geometric-only storage also refuses: that owner can add
independent axial material histories. Independently constructed bow banks do
not acquire this opt-in through source-property presence alone.

`tensile_assessment` exposes the selected elastic demand without enforcing its
limit, allowing diagnosis of overlimit candidates. This is a small-strain
normal-stress criterion at fixed material conditions and accepted temporal
endpoints. It supplies no viscous stress, torsional/shear failure, rate effects,
necking, fatigue, fracture energy, rupture topology or guarantee between time
samples. G0/G1/G3 regressions in `tests/material_string_nonlinear/tensile.rs`
specify axial-plus-bending demand, per-layer load sharing, transverse-frame
invariance, source refusals, unchanged admitted mechanics and native rollback,
including supported conversion. Native execution is a separate verification
obligation.

### Tensile admission on material fiber ribbons

`TensionedFiberRibbon::with_tensile_limit` binds every original fiber layer's
named static strength to the existing nonlinear `FiberStrain` owner. Its
retained longitudinal strain supplies the axial upper bound, including the
maintained preload and geometric mean. Original mass-normalized transverse
coordinates supply a continuous-span curvature envelope. Elastic layers use
`E_i (epsilon_max + r_outer_i K)` with their own physical fiber radii: refining
the ribbon width changes neither stress nor utilization.

Sourced Kelvin--Voigt layers retain independent original axial and bending
viscosities through `KelvinVoigtTensileCapacity`. The original longitudinal
owner supplies the complete strain rate, including `y_x*y_xt` and its mean.
For each layer, signed `E_i*epsilon_j + eta_ai*epsilon_dot_j` axial harmonics
and `E_i*kappa_n + eta_bi*kappa_dot_n` bending modes combine before their
span magnitudes are summed. Unloading can therefore cancel tensile stress;
no section-average viscosity redistributes it between layers. The maintained
preload has zero rate. `tensile_state_assessment` reports total stress, while
`tensile_assessment` explicitly reports the elastic contribution alone.
Both retain the original section/strength receipts and force per physical fiber.

Preload and the original coupled candidate gate refuse a crossed envelope;
state, histories, work accounts and time remain transactional. Forces and energy
are unchanged. Velocity-dependent limits use the actual original modal velocities
at cold binding, component-state restoration, accepted states, converged dynamic
candidates, and final staged static equilibria, before publication. Cold restores
also admit all displacement-only constitutive domains and retain their accepted
material histories. The held-force displacement increment
divided by dt is not substituted for an endpoint velocity. Mapped admission
buffers are prepared once and reused; displacement-only laws gather no extra
velocity data. The fiber assessment uses bounded stack workspaces and cached
per-layer coefficients without new heap allocation or runtime source lookup.

This remains a fixed-state normal-stress envelope at temporal endpoints, with
no guarantee between samples, rate-dependent strength, yielding, fatigue or
rupture topology. Hereditary combinations refuse in either builder order until
actual layer memories can participate. Tests in
`tests/material_multi_bow/fiber_ribbon/tensile.rs` and its `tensile/viscous.rs`
child cover analytical demand, width refinement, distinct layer viscosities,
signed cancellation, source/preload refusal, unchanged admitted mechanics,
and native cold/restore/static/dynamic admission and rollback.


## Actual worn material leaving hardness-coupled thermal receivers

`PlasticSlidingHeatContact::export_wear_candidate` transfers newly accumulated
original A/B Archard volumes to the existing
`DepositionSimulation::remove_material_before_pending_events` owner. Each
receiver supplies its original sourced caloric density at unchanged current
specific enthalpy. Native inventory checks require that density to agree with
stored mass and volume; temperature-dependent density is not silently replaced
with m/V and no thermal-expansion law is inferred. The exact caloric objects,
original material/body IDs, pressure/indentation/shear/wear histories, and clock
survive. Removed mass carries its current signed enthalpy, including negative
datums; remaining specific enthalpy and temperature do not change.

The caller authors complete remaining prism geometry, every incident contact
patch, and current original conductivities. Admission retains the same contact
link, ordered cells, actual loaded area, original prismatic face areas and axis
registration. Only normal shortening/translation is supported; rotation,
tangential translation and a silently enlarged receiver face refuse. The next
coupled heat interval uses actual reduced masses and new half-cell distances.
Original cumulative wear-fraction and homogeneous-depth envelopes are retained.

A per-body checkpoint records the cumulative requested wear already represented
by the supplied thermal binding or an accepted export. Pre-binding wear is never
removed retroactively; repeated export requests zero additional material, and
later sliding exports only its new wear. The native report separately retains
requested versus represented volume and mass discrepancies, outward enthalpy,
and stable inventory/energy balances. The contact and thermal successors must
be adopted together. Any source, geometry, density, patch, representability or
cancellation refusal leaves the inputs untouched; pending endpoint events stay
pending.

This is an explicit heat-then-export composition requiring time refinement,
with caller-supplied uniform thermal control volumes. It does not reconstruct
a spatial thermal mesh from local wear, model debris evolution or surface
creation energy, or account for mechanical work needed to move those volumes.
G0/G1/G3 regression definitions cover separate material inventories and
enthalpy, an independent two-capacity heat oracle after export, unchanged
mechanical history, no duplicate or retroactive removal, original cumulative
bounds, and transactional patch/registration/cancellation refusals.

### Adaptive prescribed sliding with hardness, wear and finite heat

`PlasticSlidingHeatContact::march_prescribed_velocity_candidate` refines a
constant actual per-cell velocity span with the original normal boundary held.
A full interval and two native half intervals start from the same accepted
history; only the half-step candidates and their original source payloads may
advance. Componentwise differences control thermal-cell energy, interface
temperature, pressure, approach, indentation, each body's recession, accumulated
slip and physical work. Relative scales use changes from the interval origin.
No extrapolated material, contact or thermal state is manufactured.

The adapter uses zero tangential body compliance. Subdividing a released body
displacement while retaining nonzero compliance would change the motion problem
and is not supported. Whole-span contact, thermal and combined work gates are
independent of per-step residual tolerances and substep count. Compensated sums
use native incremental storage and preserve signed parameter work, friction,
direct plastic heat and nonthermal plastic work separately. They do not supply
a full body thermodynamic or passivity theorem.

Attempt and accepted-receipt budgets bound work and storage. Source, numerical,
coverage, inventory or cancellation refusal publishes no partial span, including
after private accepted intervals. Native failures remain explicit; only excessive
coarse/fine defects trigger refinement. No thermal event or material export is
consumed, and cumulative removal limits/checkpoints survive. The defect is a
local estimate for the first-order coupled split, not a global error bound.

### Spatial magnetic hysteresis and finite enthalpy

The opt-in `hysteretic-eddy-heating` feature connects the original
`PlanarHystereticEddy` vector-history owner to `MagneticThermalNetwork` and
the existing finite-cell enthalpy Newton driver. Every thermal trial starts
from the same accepted potential and vector polarizations. Native per-triangle
pinning, intrinsic viscous switching and spatial Joule heat enter their actual
receiving cells exactly once. Insulating magnetic triangles can heat through
pinning and finite-rate switching without Ohmic conduction; any
dissipating triangle without a finite receiver refuses. The endpoint field's
numerical dissipation remains in the work account and never heats a cell.

Actual endpoint temperature refreshes the original conductivity source.
Every mapped magnetic fit must be source-constant over the complete original
caloric temperature interval, with the same material card, selected claims,
quantity schema and fixed external conditions. The report retains complete
fitted magnetic receipts separately from the electrical bundles. Constant
solid density and whole-cell volume correspondence remain mandatory. Stable
magnetic storage change, physical heat, numerical loss, thermal storage and
source work share one acceptance audit; failure/cancellation publishes neither
history. The shared report's `deposited_joule_j` counter denotes total deposited
physical heat in this driver; `deposited_heat_j()` names that sum explicitly,
and its native field report separates Joule, pinning and viscous heat.

An explicitly source-selected viscosity uses the same physical interval as the
thermal driver and every native field trial. The finite receiver deposits the
full constitutive heat, never the half-viscous variational objective. Its source
bundle includes all original viscosity claims, even an all-zero profile, so the
complete caloric interval must support the kinetic fit too. Native winding
circuits share this existing eddy and thermal route. Synthetic regressions
check finite-rate polarization and continued hold heating against an independent
scalar kinetic balance, exact deposited energy, retained source queries and
atomic refusal of insufficient viscosity temperature support.

This is fixed-geometry, fixed-fit, first-order magnetoquasistatic/thermal
coupling. It supplies no changing-fit thermomagnetic entropy/parameter work,
moving media, phase change, fitted core-loss overlay or experimental material
validation. Focused synthetic regressions compare insulating minor-loop
heating and conducting temperature feedback with independent scalar play and
P1 balances, distinguish numerical loss from heat, and exercise missing
receivers, source-temperature support, volume, cancellation and work budgets.

### Retained dielectric switching in native circuits

The separate opt-in `experimental-dielectric-hysteresis` feature exposes
`material_circuit::dielectric_hysteretic::DielectricHystereticCircuit`.
One original capacitor is replaced by a private native voltage port. Every
charge trial retains the same accepted source-bound polarization memory;
the exact monotone-segment mean voltage closes `delta q = dt I` through the
original Kirchhoff DAE. Native R/L/C coordinates, source ports and material
resistor bindings retain their owners. Private work cancels internally, while
resistor heat and dielectric switching heat remain distinct and enter the
combined energy audit once. Failed trials never publish charge or polarization.

Three instantaneous-current probes refuse detected charge reversals. They do
not certify that an arbitrarily large interval contains no hidden turn; the
caller supplies waveform/reversal cuts and refines threshold crossings. The
zero-branch limit is the native midpoint circuit; fixed active-set RC charging
has an independent analytic refinement regression. Geometry and coefficients
remain fixed. Reported switching heat is not automatically deposited into a
thermal cell, and this owner supplies no leakage, electromechanical strain,
spatial electric solution or empirical calibration. Six native integration
regressions cover linear RC/RLC limits, exact loop work, charging refinement,
atomic refusal, detected LC turning points and small transfers beneath large
unrelated energy biases. Both native and combined audits use the original
solver's stable energy increment.

### Finite thermal feedback for retained dielectric switching

With the same default-off `experimental-dielectric-hysteresis` feature,
`material_circuit::dielectric_hysteretic::thermal::ThermalDielectricHystereticCircuit`
couples accepted charge/polarization to the existing finite-cell enthalpy driver.
The original `PreparedMaterialNetwork` resolves caloric storage, resistance and
conduction at each proposed temperature. `HeatAssembly` providers retain that
bulk network and can supply thermal contacts, layers, closed radiation and
finite coolant. Every trial restarts the SAME accepted electrical history;
this adds no thermal, constitutive or nonlinear solver.

Switching work enters the actual whole fixed solid dielectric volume; native
resistor work enters its original material cell. An explicitly authored leakage
resistor may share that receiver: both physical losses accumulate once. Their
separate receipts survive, while `deposited_heat_j()` reports the total.
Numerical dissipation is never deposited as heat. Original circuit, finite-cell
thermal, total-energy and closed-radiation gates remain active.

Source-bound switching requires the exact original caloric card and accepted
query, plus explicit support for the entire admitted temperature interval in
every selected constant coefficient. Each trial refreshes that same fit; source
query, polarization, charge and thermal clock publish together. Both actual
cell volume and mass/density must reproduce capacitor area times gap, with a
wholly solid constant-density caloric chart. Cancellation, source domains and
combined nested circuit-work budgets refuse the complete proposal.

Endpoint resistance with backward-Euler heat is generally first order; charge
reversals and source/thermal events require interval cuts and refinement.
Electrical contact work and field-dependent thermal boundaries remain outside
this lane. No changing switching coefficients, electrocaloric heat, geometry
evolution, breakdown or empirical accuracy is inferred.
`tests/dielectric_hysteretic_thermal.rs` supplies independent two-cell
leakage/switching/R(T) balances through the native heat assembly, major-loop
finite heating, source/history continuation and atomic refusal regressions.

### Spatial dielectric hysteresis on the finite enthalpy timeline

The default-off `experimental-dielectric-hysteresis` feature also enables the
native field lane and exposes
`circuit_enthalpy::dielectric::hysteretic::ThermalHystereticElectrostatics`.
It composes the original planar Gauss/electrode history owner, source-bound
material network and finite-cell enthalpy driver; no constitutive, field or
thermal solver is duplicated. `DielectricFieldSolve` accepts the native single
endpoint controls or the existing adaptive affine-path controls through the
same thermal entrypoints. The complete original HeatAssembly may supply
conduction, thermal contacts, layers, closed radiation and finite coolant.

Every original triangle has an explicit `DielectricRegion::Cell` or `Fixed`
correspondence. Mapped triangle volumes, actual receiver volume and mass/density
must agree within the declared absolute geometry allowance; multiple triangles
may partition one whole receiver. Its source must be the same original caloric
card and complete accepted query. Fixed geometry requires a wholly solid,
constant-density caloric chart, and every selected constant switching
coefficient must explicitly support that chart's entire temperature interval.
Any nonzero physical switching heat in an unmapped Fixed region refuses.

Fixed fits and prescribed loading make the electrical endpoint/path independent
of temperature. It is evaluated once before the thermal residual, with the
native field/path work budgets retained. Each thermal trial rechecks source
queries and transport while borrowing the SAME immutable electrical account;
it neither reruns Gauss nor copies or advances polarization memory. Per-cell
physical switching heat is accumulated once into receiving-cell powers.
Numerical electrical dissipation is retained separately and is never heat.
The original field/source-work, finite thermal, total-energy and closed-radiation
gates precede publication of both histories, final source queries and clock.
Cancellation, work budgets and source-domain refusals preserve accepted inputs.

The optional path summary retains every coarse/half/rejected field solve and
operator application, plus accepted-only empirical refinement estimates. Both
the candidate getter and consuming parts preserve it. Accumulated electrical
path heat Q is deposited as Q/dt over ONE backward-Euler thermal interval;
electrical refinement does not resolve heat-release chronology or increase
thermal time order. Supply waveform/event cuts and refine the thermal interval.
There is no changing switching fit, electrocaloric transfer, field-dependent
thermal boundary law, leakage or geometry evolution in this lane.

`tests/hysteretic_dielectric_enthalpy.rs` supplies six focused G1/G3/G4 fixtures:
independent scalar-play heating in a rotated multi-triangle receiver; independent
two-cell temperature-dependent conduction through HeatAssembly; Fixed-region
zero-loss/refusal behavior; original source and whole-volume admission; late
cancellation, budgets, caloric support and replay; and cumulative refined-path
heat with separate numerical loss and retained global-work receipts.

### Resolving dielectric heat-release chronology

The same default-off spatial dielectric owner exposes
`advance_affine_candidate` and `advance_affine_candidate_with_heat_assembly`
through its `trajectory` module. One authored affine voltage/charge and volume
charge segment is now refined in physical time. Full and two-half trials use
the original coupled endpoint solver; only the two-half field and thermal
histories can be retained. Each substep deposits its own physical switching
heat while the original conduction, reservoirs and complete HeatAssembly act.
A hold can therefore refine solely for thermal cooling with unchanged electric
field, zero switching heat and unchanged polarization.

`HystereticThermalPathControls` combines the original P/Q/work refinement
policy with a separate L1 cell-energy defect and explicit global thermal-driver
evaluation cap. Thermal relative errors scale interval energy changes, never
absolute reference enthalpy. Both per-cell heat defects and thermal defects
use absolute sums so opposing spatial errors cannot cancel. Coarse, half and
rejected trials consume one global field-solve/operator/evaluation budget.
The original per-step caps remain active. Duration-weighted native energy
allowances and thermal Newton absolute targets precede whole-path physical
balance gates; no converged-state or small-update fallback bypasses them.

Only interval descriptors and a bounded set of trial states are retained. A
rejected first half is reused only for its identical left-child starting
history. Cumulative accepted-only electrical/cell and thermal energy accounts
use compensated sums; the original stable thermal residuals are preserved.
The report explicitly identifies its source sample, nonlinear reports and
charge residuals as final-substep data. Global attempted-work counts and
empirical error estimates remain in `HystereticThermalPathSummary`.
Cancellation, refinement/work refusal or source-domain failure returns neither
new history. These are first-order thermal trajectories with empirical local
estimates, not certified global time errors, automatic event detection,
temperature-dependent switching or an increase in the thermal integrator order.

The three `hysteretic_dielectric_enthalpy::trajectory` G1/G3/G4 tests compare
switching during cooling with an independent piecewise analytical solution,
check exponential cooling through HeatAssembly while every electrical error is
zero, and exercise global work caps, depth, late cancellation and exact retry.


## Accelerating stress histories through the original PCM renderer

`MovingQuadrupoleSourceWindow<M>` and `MovingQuadrupoleRenderer<F,M>` now
accept uniform motion, a borrowed `PreparedSourceTrajectory`, or the existing
mixed `CompactSourceMotion` ensemble. `QuadrupoleMotion` extends the sealed
original motion owners with complete stress-potential derivatives; the
renderer retains one source-history, admission, coherent summation and atomic
streaming loop. Its default remains `UniformSourceMotion`.

Curved paths use their original quintic kinematics and analytic jerk at each
retarded source event. Their continuous clearance and window-local Doppler
bounds come from the existing `CompactMotion` implementation. Observer, medium,
limits, source support and receiver coverage must match. Each original root
keeps its finite refinement cap; `maximum_retarded_iterations()` reports the
checked complete-window cap. Admission and constant-work endpoint/final checks
are additional bounded work. Runtime trajectory lookup and the prepared
18-coefficient pressure contraction allocate nothing in the renderer.

The receiver remains fixed and omnidirectional in a homogeneous quiescent
lossless medium. The declared receive bandwidth must account separately for
geometric modulation: the Doppler clock-rate bound does not certify a curved
source's full emitted spectrum, especially across C2 knots with jerk jumps.
History reconstruction, trajectory convergence and antialias accuracy remain
producer responsibilities. The existing decimator/WAV path consumes the actual
coherent pressure without a second pitch transform or receiver multiplier.

G3/G4 `tests/accelerating_quadrupole_pcm.rs` specifies mixed-owner retarded
clocks and chunk invariance, observer/history/Doppler/work refusal, late failure
and cancellation atomicity, and parity through the native decimator and WAV
writer. Native execution remains a separate validation obligation.

## Double-reed staple pressure recovery and jet contraction (`staple_recovery`, `contraction`)

Two optional reed-card records change how a massive reed's opening passes its flow. Both are off unless declared; absent, or at their identity value (`C_P` 0, `C_c` 1), the voice is the original one, bit for bit.

**`staple_recovery C_P S_C_M2`** (Almeida, Vergez and Caussé 2007, J. Acoust. Soc. Am. 121(1) 536-546, arXiv physics/0607011v1 pp. 8-9, eqs. 14-17), right after `slit_length`, or right after the `reed` row on a card without one, selects `ReedBoreVoice::set_staple_recovery`. A double reed's channel widens from the constriction of area `S_c` into the conical staple, and a fraction `C_P` of the dynamic head at the constriction comes back as pressure. The flow through the opening follows the slit law on the constriction drop `(dp)_c = p_m - p_c`: the viscous slit law with `slit_length`, the inviscid Bernoulli law without it (Almeida's own model neglects the viscous loss along the reed). The blade is driven by `(dp)_c` (eq. 16), not by the total mouth-to-bore drop (`JunctionPressure::blade_drop_pa` in `reed_bore::contact`), and the bore input is `p_r = p_m - (dp)_c + C_P (1/2) rho (q/S_c)^2` (eq. 17), `q` being the flow through the opening. `ViscousSlit::characteristic_junction_with_recovery` and `moving::characteristic_junction_with_recovery` solve the same-sample characteristic junction on the two laws (the inviscid one is the viscous one with a negligible slit, to 1e-9): for forward flow the quadratic term takes the factor `1 - C_P (S/S_c)^2`. Recovery acts on forward flow only (in reverse the staple is upstream of the jet) and on the JET flow: the moving blade's swept flow and any body flow are not carried through the staple. The law is singular where `C_P (S/S_c)^2 >= 1`, where it would cancel the whole jet loss; there it refuses (`STAPLE_RECOVERY_DOMAIN`), and `set_staple_recovery` refuses a reed whose rest opening already is there. The only printed numbers are one measured reed's, `C_P` 0.8 (Fig. 15) and `S_c` 4e-6 m^2 (Fig. 14): a ready-to-play natural-cane oboe reed, humid air, masses added, no lips. A card that uses them borrows them for another reed, and says so. Inside the domain the element dissipates `(dp)_total q`, never negative; because the blade feels `p_c` while its swept flow enters the bore at `p_r`, the step's energy account carries one exchange term, `-(C_P (1/2) rho (q/S_c)^2)` times the swept flow, of either sign.

**`contraction C_C`**, right after `staple_recovery` (or after `slit_length`, or the `reed` row), on a card without `slit_length` and with a massive reed, selects `ReedBoreVoice::set_jet_contraction`: the jet leaves the opening with the area `C_c w h`, so the inviscid flow is `q = C_c w h sqrt(2 (dp)_c / rho)` (the `rho U|U| / (2 (C_c w h)^2)` of `fs_phs::GapFlowLaw` without its channel terms) and the staple law takes `C_P (S/S_c)^2` on that jet area; the moving blade's swept flow is geometric and keeps the full width. `moving::characteristic_junction_with_recovery` takes the coefficient (1 is the original junction, bit for bit; outside `(0, 1]` it refuses). There is no default, and **no printed value for a double reed is held**: the clarinet's 0.6 (Hirschberg, in Silva, Kergomard and Vergez 2008) and the two-dimensional free-streamline bound pi/(pi+2) = 0.611 are borrowed numbers; Almeida's flow fit used no contraction, so an area taken from it may already hold whatever contraction there is, and a value given here is an upper bound on its effect. Hirschberg's downstream-neck resistance, which Vergez et al. (2003) carry to double reeds, moves the characteristic the other way from the measured oboe (Carral, Vergez and Nederveen 2011). The channel's unsteady inertance (`GapFlowLaw::inertance`) is not part of this record: no printed effective length of a double reed's channel is held.

Tests: the flow maximum of the printed characteristic sits where the closed form of eqs. 15-17 puts it; the domain refusal by name; the junction's power balance and passivity; the inviscid junction against the viscous one with a negligible slit; the contracted junction against the manufactured flows and against `GapFlowLaw::volume_flow` with the same contraction; the contraction moving the staple domain to the jet area; the same-sample flow closed on the contracted area; the blade momentum under `(dp)_c` and the step's energy account; bit-identical absence for both rows; callback invariance; the card rows' refusals.


## Source-bound plastic shell and spatial enthalpy feedback

`shell_plastic_enthalpy::ThermoplasticShell` binds each physical layer of every original CST/DKT facet to a unique active enthalpy cell, its original caloric object, reference volume and sourced density. The native plane-stress section retains independent complete material histories at all thickness points of all three surface points. Source families resolve actual thermal expansion, affine energetic modulus/hardening changes, yield softening and selected Voce/Perzyna/Chaboche kinetics. The native material account includes coefficient free-energy work and mechanical entropy; physical/reversible heat is integrated conservatively to the original cells and supplied exactly once to the existing conduction/observer timeline.

`advance` and `advance_backward_euler` share frozen-history equilibrium and thermal proposals. Every trial starts from the same accepted coordinates, complete plastic histories, enthalpies, observer state and, for dynamics, velocity. The actual interval owns both rate-dependent material response and native momentum. Only matching temperature/heat and independent force/support, thermal and combined first-law gates publish the shell, heat, motion and observer candidate together. Signed transverse residual work stays explicit. Material numerical loss, retained reference drilling storage/loss and physical/gauge kinetic integration loss never heat a cell. The original consistent mass is prepared only when requested, retained across intervals and checked against each original layer density, not merely total mass.

Physical sources and facet bindings are shared immutably; borrowed endpoint sections avoid material-law clones, and reference mass/gauge operators outlive temporary thermal laws. Source-bound conductivity selections lower to the existing material-network observer, retaining its events, reservoirs, conductances and interval hooks. The reduction remains fully solid, fixed geometry/density, small strain and rotation, with one temperature per physical layer/facet and explicit thickness quadrature. Finite forming, melting, transverse-shear plasticity, evolving transport geometry and delamination are not claimed. Physical interval, thermal substep, spatial cell and thickness refinement remain independent accuracy requirements.

Four focused synthetic regressions cover an independent heated uniaxial plastic membrane/permanent-set oracle, layer-gradient thermal moment reversal with original conduction, late-refusal checkpoint preservation, and prescribed plastic motion with an independent P1 kinetic integral, incorrect eccentric mass rejection and no heat from numerical stopping loss. They also exercise a cold reference binding whose free-coordinate budget is smaller than its total physical coordinate count. Formatting and source checks do not substitute for execution of these native numerical tests.

## Spatial thermoplastic indentation hardness

`solid_enthalpy::plastic::hardness` binds one independently selected original
indentation test to each tetrahedron of an accepted `ThermoplasticSolid`.
`bind_hardness_field` bounds point count and retains original mesh, ordered cell,
caloric source and test identity. `at_accepted` observes each actual caloric
temperature and native accumulated equivalent plastic strain, retaining complete
source receipts in original element order. It preserves full J2/Voce/viscous and
N-branch Chaboche histories, including reverse loading and thermal branch transport.
There is no averaged temperature, yield-to-hardness conversion, repeated return
map, repeated mechanical coefficient refresh or second mechanical history.

`check_candidate` admits an existing native thermoplastic proposal and returns
that unchanged proposal with its complete hardness field. Source coverage,
inventory, clock, original-mesh and cancellation refusals publish no partial
field or advanced physical owners; callers adopt both outputs together. These
are endpoint observations, with no assertion of coverage between physical
times. Consumers still own spatial projection, contact/indentation scale and
avoiding duplicate resolved and reduced indentation energy; no additional
pressure cap, wear law, annealing kinetics or heat term is introduced.

The focused synthetic `solid_plastic_enthalpy::hardness` regressions specify
heterogeneous element responses against an independent uniaxial plastic/heat
oracle, complete cyclic branch and native work preservation, and late endpoint
coverage refusal with accepted-owner preservation. They also exercise bounded
point count and original-mesh/clock admission. These are numerical source/history
checks, not measured material calibration.

## Ductile initiation on coupled thermal plastic solids and motion

`solid_enthalpy::plastic::ductile` adds `advance_with_ductility` to the original
`ThermoplasticSolid` and `ThermoplasticDynamics` owners. The existing native
proposal retains all thermal, equilibrium, source refresh, physical/reversible
heat, inertia and energy-balance behavior. Only after that complete proposal
converges does `DuctileHistory` assess the actual endpoint stress and accumulated
plastic alpha. Complete J2 and Chaboche histories remain the original states;
backstress is not projected away and reverse flow accumulates ductile usage.

Each point's physical temperature is evaluated from its retained cell enthalpy
through the SAME original caloric source. The shared physical duration is the
target time minus the accepted coupled time, including after restart. Material
kinetics, momentum and ductile plastic-strain rate use that interval. No ambient
temperature, iteration time, total strain or kinetic/numerical loss substitutes
for those quantities. The bridge allocates one bounded temperature vector per
endpoint and performs no extra material-source query or constitutive solve.

Before a proposal, the law is admitted against the immutable original thermal
source reference and accepted endpoint material. The converged endpoint source
is admitted again before acceptance. Card identity, named state, all fixed
context values and quantity types remain binding. A resolved fracture table at
one fixed temperature refuses a changed temperature, even when that change lies
inside its source validity box. Re-resolving a table at the latest temperature
does not replace original source authority. Explicit authored Johnson--Cook
coefficients may carry temperature/rate dependence across their declared domain;
this does not infer material calibration.

The returned `DuctileAccepted` contains the complete native coupled candidate
and corresponding initiation history. The caller adopts them together. Source,
temperature, history-reset, initiation, native solve and cancellation refusals
publish none of the proposed enthalpy, displacement, full plastic state, velocity
or observer state. Callers refine from the complete accepted pair to locate
initiation. Endpoint-held accumulation is first order and its crossing fraction
is plastic accumulation, not resolved event time. This gate supplies no
post-initiation degradation, fracture heat, topology change or fragments.

Five synthetic regressions in `tests/solid_plastic_enthalpy/ductile.rs` cover
independent Johnson--Cook temperature/rate utilization over successive physical
intervals, equality with the original native coupled solution and heat accounts,
initiation/refined retry/history-reset rollback, inertia and numerical-loss
separation, complete cyclic Chaboche backstress, and original-card/fixed-temperature
source refusal including attempted endpoint rebinding. All five passed against
the original native material, solid and thermal owners in the focused source
closure; this is not a full-workspace Cargo test result.

## Source-bound Timoshenko frame self-heating

`frame_vibration::thermal` and its complete exterior-fluid thermal consumer
retain native axial, symmetric bending/shear and antisymmetric bending ports.
The `fs-solid` owner supplies their common damping and original-layer work
fractions. `FrameLayerHeat` reports axial, curvature and transverse shear heat
separately; corrected shear viscosity follows its original layer area, while
curvature loss follows its original second moment. Their sum is the actual
native dissipated work, deposited once into the same finite caloric inventory.

Temperature refresh preserves original viscosity selectors, named/derived G
parents, frequency support and the authored shear correction. Geometry, density
and the actual source parents of E and G must be temperature independent in
this fixed-elasticity owner. Ordinary thermal feedback uses accepted starting
viscosities per interval; the fluid thermal owner resolves endpoint viscosity
and the complete fluid-memory trajectory jointly. Their existing work, entropy,
finite-cell, source-domain and atomic publication gates remain in force.

Live thermoelasticity retains each original Timoshenko interpolation and
consistent mass as a fixed reference Galerkin basis. The original
`ConcentricRodThermoelasticLaw` supplies affine E(T)/G(T), integrated expansion,
Helmholtz energy and entropy. Because each layer has uniform temperature and
the potential is quadratic, its RMS-curvature evaluation integrates the whole
member exactly; a separate endpoint-curvature and shear check enforces the
actual strain envelope. Named and derived G roles remain source-owned.

Each layer's symmetric spring combines bending and shear storage on the same
native coordinate. Updated viscosity, heat partition and stiffness use the
retained reference strains, while the existing complete coefficient update
checks source parameter work, Psi+TS, finite-cell heat and atomic publication.
The original complete structural basis and the same fluid/electrical/shared
assembly owners remain in force. This is a fixed-shape beam approximation:
member refinement is required when changed stiffness ratios alter the
unresolved physical fields. It does not update geometry or inertia, infer
nonlinear G derivatives, or claim exact continuum dynamics.

Focused synthetic regressions cover independent Timoshenko oscillator and
strain-work references, affine E/G thermoelastic energy/entropy exchange,
distinct layer temperatures, shear-only dissipation, full exterior-fluid
memory and warm-state refusal/retry. Native fixed-basis regressions separately
compare layer energy, loss and parameter work with integrated physical strains.


## Source-bound creep and reciprocal heat in layered shells

`shell_plastic_enthalpy::ThermoplasticShell` also accepts `CaloricCreepSource` through the shared `ThermalSectionSourceFamily`. Original standalone or joint Norton sources retain primary resistance, Arrhenius kinetics, optional sourced expansion and all native creep/plastic memories at every surface/thickness point. A standalone creep card requires no invented yield or plastic-hardening properties. Each original layer retains its own caloric cell, fixed reference mass, material identity and temperature; native conduction/observer, equilibrium, momentum and atomic feedback owners are reused.

`ShellThermoplasticEnergy::creep_heat_j` and `ShellLayerHeat::creep_j` report native creep work already included in total physical heat. Only total physical plus signed reciprocal heat enters the caloric timeline, exactly once. Numerical material/inertial loss and drilling energy remain excluded. Accepted thermal histories must match their actual enthalpy-derived temperatures, and all trials restart from the same accepted state and physical clock. Original thermal transport pinning and static-domain checks apply to creep sources as well as plasticity sources.

Two focused synthetic regressions define an independently reconstructed load/hold/unload/reload primary-creep membrane with nonzero expansion and thermal activation, late-refusal preservation, and incompatible caloric/history/card rejection. They are regression definitions, not executed results. The original fully solid, small-strain/rotation and fixed-geometry shell reduction remains bounded; measured material calibration, tertiary creep, damage, rupture and qualified service life remain separate work.

## Tensile envelopes with longitudinal waves and moving string ends

`SupportedMaterialStrings::in_stationary_gas_with_longitudinal_waves` carries
each source's selected tensile capacity into its original longitudinal
potential. The native `fs-nlmodal::longitudinal` strain polynomial supplies
`T_max = T0 + EA (epsilon_mean + sum_j |epsilon_j|)`. Actual oblique endpoint
extension, both transverse polarizations, their affine-chord cross terms and
the explicitly selected solid-inertia axial modes all contribute. Testing
only mean tension would miss a local tensile peak. Unretained axial channels
remain statically equilibrated as in that original reduction; this extension
does not supply their inertia or claim converged propagation.

The existing common tensile kernel combines that continuous-span axial bound
with each original layer's `E_i`, outer radius and selected static strength:
`sigma_i <= E_i (T_max/EA + r_outer_i K)`. Curvature factors are cached from
the actual transverse normalization, independently of solid axial inertia or
sparse axial mode indices. In this elastic lane the guard adds no successful-check
allocation, new force, stiffness, energy, mode, clock or solver. Native initial admission
and accepted-step gates enforce it after the original kinetic pullback;
generic copied bodies retain it and a failed trial publishes no partial state.

`SupportedStrings::tensile_assessment` diagnoses elastic candidates in the
original support coordinates without enforcing their selected strength limit.
`MaterialSupportRun::tensile_assessment` and
`ImpactSystem::support_tensile_assessment(body, wire)` read accepted native
state, including actual normalized kinetic velocities and source-bound axial
material histories when selected.
The diagnostic uses the original material-memory offsets even when other
bodies or contact coordinates precede the supported strings. Displacement-only
assessment refuses when the selected criterion needs velocity or hereditary
history. `StringTensileAssessment::tension_n` is the upper local tension bound:
elastic in the original lane, and complete equilibrium-plus-Maxwell-plus-selected
axial-dashpot tension when those source laws are selected.

`with_material_longitudinal_relaxation` retains every active arm's original
layer index, modulus increment, physical area and normalized memory address.
For each layer and retained cosine coefficient, it first sums the signed
equilibrium and actual Maxwell stresses. Only then does it take the
continuous-span harmonic envelope. Its mean includes the original maintained
preload `E_i * T0 / EA`; a branch contributes
`DeltaE_i,a * (epsilon_j - z_i,a,j / sqrt(c_i,a,j))`, where
`c_i,a,0 = DeltaE_i,a * A_i * L` and nonzero cosine channels use half that
value. Layer stresses can reinforce or cancel, and aggregate force divided
by equilibrium EA is not substituted for their load sharing.

Independently selected source axial viscosity contributes each original layer's
signed `eta_i * epsilon_dot_j` before that same harmonic envelope, with or
without Maxwell arms. Copied native bodies retain exact section, material-query
and viscosity receipts; cold binding caches each layer's viscosity. Original
normalized kinetic velocities cross the same reciprocal support-coordinate
map as displacement. The selected linear or finite-strain law supplies its
native strain rate, including the corresponding moving-end contributions.
The guarded local total-tension bound includes the same dashpot contribution.
No strain rate or layer stress is inferred from displacement increments or
aggregate damping forces.

The complete-loss source constructor binds viscosity before initial tensile
admission, so a physically relieving signed dashpot stress is retained. The
existing mechanical gate handles viscosity-only wires without creating any
relaxation state; actual Maxwell histories continue through their original
memory gate. A public empty-arm declaration remains a physical no-op and
cannot suppress initial or accepted-step tensile admission. Raw aggregate
viscosity without original layer authority still refuses for a selected guard.

The same common tensile kernel adds each layer's elastic outer-fibre bending
stress and compares the result with its selected static source limit. Guarded
axial memory is checked in the original memory acceptance gate, after actual
candidate histories are available, without an earlier elastic-only tensile
refusal. Initial and failed trial admission preserve the original transaction
semantics. The layer projection adds no physical force, stiffness, energy,
state, clock or integration method. Public aggregate-EA arms without original
layer authority remain refused for selected tensile criteria.

The envelope can be conservative when axial and bending maxima occur at
different positions. It remains a fixed-condition normal-stress estimate,
with no viscous bending stress, shear failure, fracture evolution or guarantee
between accepted temporal samples. Source bending-memory variants still
require their actual bending-history stresses. No rate-dependent strength
calibration or additional failure law is inferred from the selected viscosity.

G0/G1/G3 cases in `tests/material_string_nonlinear/tensile/supported.rs` and
its `hereditary.rs` child specify independent sparse-wave and endpoint-extension
formulas, two-plane rotation, bonded-layer load sharing, relaxed and unrelaxed
memory stresses, continuous-span containment, original body/history offsets,
source mismatch refusals, and native/prepared rollback on an actual Maxwell
tensile crossing that the elastic-only criterion would miss. The native memory
strength unit regression specifies signed coefficient cancellation before the
spatial bound. The `hereditary/viscous.rs` cases additionally cover signed
viscous cancellation and reinforcement, both selected loss laws with and without
Maxwell, initial relief, support velocity, preceding-body offsets, raw-source
refusal, and native/prepared rollback on a rate-driven tensile crossing. They
also verify that empty-arm attachment cannot remove the viscosity-only gate.


### Coupled sustained plastic/creep shell contact and finite interface heat

`shell_plastic_enthalpy::contact` advances the same source-bound layered shell, native sustained normal/Coulomb contact, optional consistent-mass momentum and original finite enthalpy network at one shared endpoint. The existing `ThermalSectionSourceFamily` preserves complete plastic, cyclic and creep histories; the common section audit and exact native resultants determine bulk physical/reversible heat, including creep only once within physical heat. The existing friction-heat network solves each massless common interface temperature and both finite-body heat transfers from original conductivity/geometry, with no prescribed heat share or added interface storage. Endpoint friction coefficients and transport refresh at the actual candidate temperatures; the callback must not hide unresolved pressure/rate/history dependence.

Every mechanical witness maps once to the original exposed layer and a distinct rigid-obstacle thermal cell. Reference skin/face positions, inward force normal, prism thickness, area and original inventory must agree; fixed footprints must remain closed along the admitted affine gap path. Pressure is native solved normal force divided by explicit area, never inferred from the shell's plane-stress normal stress. The whole-system gate closes material/caloric/kinetic storage, actual support/dead-load/obstacle work and outside heat. Material, drilling, inertial and normal-approach numerical losses remain separate from physical heat. Thermal roots stage all histories and receipts, retain the exact final mechanical response, and publish only a complete balanced candidate. One bounded cache reuses exact repeated layer enthalpies and Coulomb coefficients while source/transport receipts always refresh.

This bounded path supports fixed reference skin patches and prescribed obstacle motion; opening/changing footprints, evolving contact resistance, general impact events, wear and flash temperatures are outside it. Scheduled thermal events require interval splitting. Three authored native regressions cover independent finite-cell heat partition with temperature-dependent friction and bulk yielding, consistent-mass dynamics with numerical loss excluded from heat, and late source refusal preserving accepted states. Formatting and source/API/energy review passed; native execution remains pending the already observed shared dependency build resource blocker.

## Recovering primary resistance in solid and shell heat feedback

The existing generic `CreepSolid` and `ThermoplasticShell` owners accept
`RecoveringNortonSource`, including sourced expansion through
`ExpandingCreepSource` and `CaloricCreepSource`. Each material point retains
permanent signed creep, its nondecreasing accumulated path and a separate
recovering resistance. Temperature resolves both original Arrhenius laws over
the same physical interval. All thermal proposals restart from the complete
accepted histories; a source resolution or inspection cannot recover resistance
or advance creep by itself.

Physical creep work enters caloric storage once. Unloaded resistance recovery
alone supplies no heat; reciprocal thermoelastic heat and material/inertial
numerical losses retain their original partitions. Two manufactured native
load/hold/unload/dwell/reload regressions use independent scalar references to
check deformation, both histories, endpoint temperature, complete source
receipts and atomic refusal. The shell case includes nonzero sourced expansion.
These regression definitions have not yet run on a prescribed build runner;
their independently checked scalar trajectories are not production test results.

## Reciprocal caloric feedback through finite ductile softening

`solid_enthalpy::finite_ductile::FiniteDuctileThermalSolid` binds the original
finite ductile tetrahedral specimen to its original fully solid caloric cells.
`ReferenceCaloricMaterial` shares the existing finite thermoplastic source
admission: immutable same-card caloric object, complete typed fixed context,
original pinned Cp/density receipts, positive heat capacity and an explicit
stress-free zero-pressure caloric reference. Original reference volume and
sourced reference density must reproduce cell inventory. Optional native
consistent inertia uses that same density once; current temperature, damage
and deformation never replace the mass inventory.

The native nonpotential finite solve retains its original material, plastic
metric, initiation/damage history and velocity through every coupled proposal.
The existing finite thermoplastic Picard/IQN-ILS driver advances the original
thermal timeline, observer, events and transport. Actual caloric endpoint
temperatures and the SAME physical interval reach initiation and momentum.
Each proposal starts from the same complete accepted state. Coupling secants
belong only to the current interval and cannot advance physical history.
Source-backed initiation must match the original caloric material reference;
a calibration fixed at one temperature refuses an unsupported heated endpoint.
Authored numerical initiation retains its explicit domain without acquiring
material-card authority from the caloric source.

Compact element receipts come from the existing final constitutive audit, with
no repeated return mapping or subtraction of large cumulative dissipation.
Explicit fractions partition ALREADY-ACCOUNTED plastic dissipation and damage
release between caloric heat and two retained nonthermal irreversible energy
inventories. These are not Taylor--Quinney fractions of total plastic work;
native energetic hardening is already mechanical storage. Unheated energy is
neither discarded nor claimed to be measured fracture-surface energy, and no
later release law is inferred. Signed material endpoint defects and numerical
kinetic loss never become heat or nonthermal physical energy.

Independent plastic/damage partition checks precede the shared temperature,
sum-of-absolute-cell-heat-discrepancy, caloric-balance and combined-work gates.
The combined account retains mechanical storage, caloric storage, nonthermal
energy, kinetic change and both numerical accounts without charging physical
dissipation twice. Only the complete mechanical, thermal, observer and
nonthermal state publishes. Refusals preserve the original state for retry.
Initial binding refuses prior plastic/damage loss without its original coupled
history; subsequent intervals continue from the returned complete owner.

This extension uses the supplied fixed-coefficient finite J2 backbone, with
temperature-dependent initiation. It adds no thermal expansion, native
coefficient softening, mechanical entropy, evolving Fourier geometry, phase
change, crack topology, contact or fragment transfer. Reference-network
conductances and the heat fractions remain explicit model inputs. Coupling
convergence and first-law closure do not certify temporal or experimental
fracture accuracy. The focused definitions in
`tests/finite_plastic_enthalpy/ductile.rs` exercise independent scalar heated
initiation, subsequent release of stored energy, retained nonthermal energy,
held replay, simultaneous free-coordinate motion/heat at two physical rates,
original-source/domain/retry and cancellation refusals.

The five new coupled ductile regressions pass in the selected native execution
closure. Six existing thermoplastic regressions also pass; five older selected
cases fail identically when the original `d624b404` coupling/material/feedback
bodies are rebuilt against the same native dependencies (four native Newton
stalls and one caloric-density fixture assertion before coupling). The full
original finite ductile solid integration target passes all twelve tests.
These results do not claim a full-workspace Cargo or constellation pass: the
ordinary Cargo attempt was killed in the pinned upstream execution dependency,
and the constellation preflight reported pre-existing NumPy/pandas checkout
mismatches. No native physics, solvers or cancellation owners were substituted.

The shared source/feedback owner also preserves the concurrent finite porous
thermal adapter. Reconciled native material/solid/coupling production bodies
compile, and all five ductile thermal tests still pass. Three upstream porous
regressions pass; the fourth adiabatic-crush Newton stall reproduces with a
byte-identical residual history using the unmodified `23de8912b` upstream
feedback/material/porous bodies and identical native dependencies. No observed
new failure is attributable to this integration; the existing failures remain
explicitly outside a passing-workspace claim.

### Simultaneous recovering creep and plasticity in thermal solids

`RecoveringNortonJ2Source` composes with the existing `ExpandingCreepSource`
and native `CreepSolid` without a second material clock or projected history.
The complete fifteen-claim source and both permanent strain fields remain with
each point while primary resistance recovers independently. Physical heat is
exactly the native creep-plus-plastic dissipation; reciprocal expansion heat,
`T` times mechanical entropy storage and numerical loss stay separate. Neither
source refresh nor an unloaded recovering dwell invents heat or resets plastic
hardening. All candidate thermal roots reuse the complete accepted state.

A manufactured load/hold/unload/heated-dwell/reload regression independently
solves primary flow, plastic consistency and the thermal first law, checks all
three histories and retained receipts, and verifies proposal evaluation leaves
accepted cells, solid solution and observer unchanged. Native Rust execution
remains pending; the independent scalar reference is not a production test run.

## Source-owned body temperatures in native solid friction and heat

With `experimental-solid-contact-heat`,
`friction_heat::solid::propose_solid_contact_heat_with_body_temperatures`
maps each additional declared friction-source temperature axis to side A or B
of its original finite thermal route. Every joint trial queries the actual
endpoint caloric cell temperature while retaining the solved common interface
temperature, native pressure/speed/distance history and original source claims.
The one-result mechanical cache includes every queried interface/body
temperature; body heating or cooling cannot reuse an older friction response.

Missing, duplicate, foreign-node, primary-axis, nonthermal or undeclared
temperature ownership refuses. The original entry remains available with no
additional mappings. Fixed-footprint plane-strain mechanics, physical depth,
deposited heat, exported work and the original energy gates retain their
existing owners. No prescribed body heat split or interface storage is added.
The extended native regression covers independent two-temperature friction,
both side assignments, fixed-interface-temperature cache liveness and physical
work/heat balance. Five exact-source admission checks passed; the full native
regression awaits the pinned sibling dependencies and prescribed build runner.

## Anisotropic thermal plastic shells

The original `shell_plastic_enthalpy` family boundary accepts source-bound
thermal Hill laws with complete local-axis histories. Layer equilibrium,
reference drilling/inertia, source-resolved thermal feedback and conservative
cell heat remain owned by the existing native paths. No J2 yield/hardening
claims or state projection are required for an explicitly sourced Hill card.
Physical yield heat enters the original caloric graph once; storage changes,
mechanical entropy, reversible exchange and numerical losses stay separate.

`tests/solid_plastic_enthalpy/shell/hill.rs` exercises heated rolling and oblique
sheets under two loaded intervals and actual free unloading. An independent
uniaxial/caloric oracle checks directional flow, softened reference yield,
retained plastic memory, transverse contraction/shear, endpoint temperature,
physical heat and permanent dimensions. Source fixtures are authored numerical
cases. Scope remains fully solid, fixed material metric/frame and small strain
and rotation; the original source, geometry, step and energy refusals apply.

## Magnetized plasma absorption on finite electron/heavy gas energy

With `ionized-wave-heating`, `MagnetizedIonizedGasHeating` binds the original
magnetized half-space slab absorption to the exact original frozen gas parcel.
Its physical area/depth volume must equal the retained gas volume to roundoff;
no density, area or heat rescaling is inferred. The parcel is explicitly
uniformly mixed, while its original wave absorption is integrated exactly.

The native `fs-gas::chemistry::energy` owner receives the two signed carrier
powers once over a declared positive held-field interval. Original species
calorics supply energy differences and heat capacities; every present neutral
thermal target comes from the SAME retained electron momentum cross section.
Effective scalar electron drag alone is refused because it cannot identify
target-resolved energy transfer. The original Coulomb law supplies baseline
thermal exchange, while its signed drift correction is already in the wave
pool powers. Internal exchange cancels; coherent carrier kinetic storage,
background field storage and gyroscopic reactive power are never heat.

Native adaptive implicit exchange, source/domain and energy admission, and
physical versus numerical entropy accounts remain intact. Caller-declared
temperature increments bound the held-source step. Publication refreshes every
original carrier source and the original bound-charge card at actual Th/Te/p,
preserving immutable claims and source limits; failure or cancellation returns
no partial thermal/material result. The returned transport and background can
seed a new Maxwell solve with its own complete-field drift checks.

This is a fixed-volume, frozen-population, cycle-averaged interval with endpoint
source refresh. It does not assert continuous implicit field/thermal feedback,
source validity along an unexamined material path, a second film interface,
species evolution, spatial thermal conduction, flow or transient carrier
energy dynamics. Five authored regressions cover Maxwell-to-caloric closure
and a new field solve, signed Coulomb allocations and internal exchange,
zero-field thermal relaxation, original source/volume/domain binding,
temperature limits, cancellation and refusal of scalar-only electron drag.

## Two-axis string bridges in omitted-mode contact response

Each material string's original plain bridge plus one independently directed
bridge shares the full two-support omitted-mode Green matrix
`G_ab = (d_a dot d_b) G(s_a,s_b)`. Nonorthogonal and parallel axes, distinct
stations, and distinct receiving bodies are admitted. Static spring reduction
is `K_eff = (I + K G)^-1 K`; viscous reduction remains the first-order
low-frequency approximation `C_eff = T^T C T`, `T = (I + G K)^-1`. This does not
resolve omitted resonances or a full frequency-dependent viscoelastic tail.

Positive rank-one factors use at most three existing signed modal connection
columns, with the original native connection/setup budgets. There is no new
integrator. Rank-deficient compliance and zero spring/dashpot coefficients need
no inverse of G or stiffness/compliance floor. Shared double-double arithmetic
retains a weak conditional direction; only negative correlation roundoff within
128 machine epsilons is projected onto the passive boundary. The original
orthogonality admission tolerance (1e-12) and scalar/orthogonal arithmetic remain.

Fixed residual bow and stop maps use the same elastic factors in the complete
Schur complement.

Fixed residual contacts retain the affine displacement caused by the original
bridge rest extensions. The same elastic Schur factors supply
`offset_i = g_i^T K_eff r`, including scalar and orthogonal supports as well as
the full coupled two-support case. Opposite contact normals reverse this shift;
the complete reciprocal force columns and residual compliance are unchanged.

A normal stop lowers its native reference clearance to
`gap_effective = gap_original - offset_i`. This is the coordinate substitution
`U(N^T q + offset_i - G R - gap_original)`, not a new contact law. The original
material-bank stop retains its authored clearance, coefficients, quadrature
weight, physical loss and provenance. Native potential/secant evaluation,
penetration/material admission and released static preload use the same
lowered law. Source binding after preload checks the complete accepted closure,
including the omitted modes' displacement from the accepted reactions. It
cannot validate retained displacement alone or read disposable trial reactions.
Negative effective clearance represents initial interference and retains its
contact energy; it is never clamped away.

A fixed offset cancels analytically from bow/shear-film slip increments. It
adds no artificial initial shear strain, dissipation or actuator work. The
native network retains the full bridge/rest elastic energy, and omitted-tail
storage and numerical loss keep their existing accounting. Zero offset keeps
the original law and clearance bits. Rest extensions remain fixed source data;
time-varying support rest, migrating residual maps and existing incompatible
nonlinear/history lanes remain unsupported.

Resolved torsional termination springs/dashpots may coexist with scalar or
two-polarization bending residual supports, in either binding order. The actual
source-derived angular row has an exactly zero transverse prefix; rotational
lowering checks that invariant explicitly. Eliminating omitted bending modes
therefore leaves the original torque law, angular receiving-body map, rest angle
and resolved rotary inertia unchanged. The complete native modal operator still
exchanges twist/transverse energy through a shared body.

The original work-conjugate lowering remains `x = r theta`, `F = tau/r`,
`k_x = k_theta/r^2`, with the material section's actual outer radius. Native
spring/rest energy and viscous torque loss remain with the original network.
No omitted torsional compliance or extra torsional loss is inferred from the
bending Green function. This does not add a six-DOF or geometrically exact rod
termination, contact-foot separation, or dynamic unresolved poles. Bow residual
corrections leave the retained twist row unchanged; normal stop corrections
retain zero twist participation. Every torque law consumes its own native
connection column, alongside the coupled bending factorization.

More than one extra directed bridge per string, constant-loss residual bridges
and the existing mixed-contact exclusions remain outside this lane. Preparation
occurs once at binding, with no additional per-sample projection or integration
path.

Focused execution: 54 actual-source angular/torque/Schur/column cases at 4/16/64
retained bending modes, three support-axis relationships, two receiving bodies
and signed rest angles agreed with independent 100,000-mode bending tails and
the original physical angular work/energy. Full bank binding-order, torque-driven
release and bow/stop continuation regressions are authored in
`tests/material_string_torsion/polarization/rotational_residual.rs`; their native
integration execution remains unperformed in this resource-constrained session.

Support corrections retain every original receiving body. The generic friction
and normal-contact owners accept one sparse collection of signed attachments
per contact, with unique ascending original component indices; same-body rows
are combined by the material producer. Their sum enters the complete reciprocal
force/displacement column before any joint response or cross-compliance is
computed. Empty collections add no support. Total support-map work is bounded
by the original retained mode count and occurs only at binding. The actual
opposite contact surface remains independent and retains its reaction. The omitted bending
response remains quasi-static, and its numerical loss stays separate from
physical dissipation. Existing restrictions on moving contacts and incompatible
nonlinear/history mechanisms remain in force.
`tests/material_string_torsion/polarization/residual.rs` specifies a bowed X/Y
mirror with energy checks, an independent omitted-mode sine-series comparison
for oblique bow cross-compliance, active oblique-stop rotation covariance, and
equivalence when independent body coordinates are stored together or separately.
The existing material shear-film and coupled-stop residual tests independently
check each support's signed reaction, physical contact work and energy balance.

Three focused mathematical tests execute the actual two-support condensation,
including rank-one compliance at stiffness 1e30 and a weak nonparallel direction.
A bounded actual-source probe also compares 120 coupled cases against an
independent 100,000-term sine tail and matrix inverse: passive elastic/damping
operators, signed native columns, contact Schur maps and virtual work agree
within their recorded numerical tolerances. New complete-bank and contact
regressions cover nonorthogonal supports and native energy/storage behavior;
those material-admission and dynamic integration tests have not been executed.

## Nonlinear port harmonics and sidebands through physical BEM receivers

`phs_orbit::acoustic::PortAcousticSpectrum` connects the original nonlinear
port-spectrum and infinitesimal sideband observers to the existing prepared
BEM pressure/particle-velocity banks. Ordinary harmonic peak phasors retain
their negative-time convention. Positive-time sidebands fold at their actual
signed carrier-plus-orbit frequencies: positive bands conjugate, negative bands
do not, and coincident physical lines combine coherently before source mapping.
DC remains explicit in the original port units, without an invented static
acoustic solution. No amplitude threshold removes a retained harmonic.

The caller supplies each exact-frequency source-by-port map in declared SI
source units, with its physical sign and phase. All banks must share the given
receiver identities and geometry. The returned `PreparedHarmonicField<3>`
receivers use the existing finite-window interference/intensity and calibrated
PCM owners, retaining inharmonic beats and the original source clock. Preparation
and propagation admit complete count/work/storage budgets and cancellation;
they perform no new field solve. Native BEM identities and producer uncertainty
remain borrowed in the result. This is observation of the supplied periodic
regime or small probe, not radiation feedback, an attack transient, source-unit
inference, stability proof, or a new source/spatial accuracy certificate.

`tests/phs_harmonic_acoustic.rs` exercises an original cubic-storage output
with DC/even/odd harmonics through actual incident-plus-scattered BEM fields
and the PCM stream, conjugate sideband collisions, independently integrated
finite-window beats and vector intensity, static-only output, and physical
binding/budget/cancellation refusals. Forced and autonomous cycle-work audits
also refuse grids that cannot resolve the authored input or retained state
harmonics; nonlinear power still needs successive-grid convergence evidence.


## Exact compact radiation from linear string harmonics

The fixed-end linear string renderer and its material, thermal, concentric and
finite-gas runtimes evaluate compact-strip pressure from the accepted modal
state and authoritative damping rates. Under held force and coefficients,
`a = f - omega^2 q - d v` and `j = -omega^2 v - d a`. The original signed area
and first-moment weights retain odd-mode acceleration and even-mode jerk
radiation. The observer adds no allocation or sampled derivative history, so
even partials avoid backward differences' attenuation and half-sample phase lag.
Material and gas updates evaluate the same expression under their final held law.

Constrained moving-end interval-mean observations retain their clock. This
remains a compact strip observer without radiation loading or common
propagation delay; ideal force, release and parameter-jump impulses require a
resolved transition and are not
reconstructed. Analytic high-band amplitude/phase, two-grid public-runtime,
damped held-load and material/gas-update regressions cover the changed behavior.

## Complete nonlinear string endpoint radiation

Fixed-end Kirchhoff--Carrier strings now evaluate acceleration and jerk at the
proposed accepted endpoint using the original full pHS law. With constant
structure and held input, `xdot = (J-R) grad H + G u` and
`xddot = (J-R) D(grad H) xdot`; canonical mass-normalized momentum rows supply
the observer's modal acceleration and jerk. The storage-owned continuous
tangent retains shared cubic stretching, bending, elastic distributed contact
and all Maxwell memory coordinates. Four reusable state buffers and one zero
input buffer are prepared once; this observation and these admitted borrowed
storage actions allocate nothing while running. Pressure admission precedes
mechanical publication and preserves the existing texture rollback.

This differentiates the physical storage at the numerical endpoint, not the
numerical step map, and does not make nonlinear propagation exact. Existing
bow, tangential-friction and Hunt--Crossley input forces are held over each
step: their continuous feedback derivatives and ideal jump impulses are not
claimed. Alpha=1 contact at exactly zero penetration and nonzero normal
velocity has no classical two-sided effort derivative and refuses explicitly;
quiet contact, active/separated contact and alpha>1 entry retain their actual
derivatives. No finite-difference fallback changes that meaning. Joined
moving-end reactions remain interval-mean quantities. Compact-strip geometry,
sampling/aliasing limits and the existing polarization model are unchanged.

Independent force/jerk regressions cover bending dispersion, shared stretching,
modal damping, held-load changes, mixed active/inactive elastic contact, and
radiation caused by memory relaxation with zero instantaneous modal velocity
and acceleration. A renderer regression checks amplitude and phase of an even
partial at 80% of Nyquist; contact-corner coverage checks the named refusal.
These native regressions require execution through the real dependency graph.

### Coupled rod fatigue with bonded layered transverse shear

`RodMotionFatiguePath::new_with_layered_transverse_shear` selects the original
solid observation owner's bonded common-Poisson circular recovery for native
motion. Distinct original layer moduli and source identities reach the same
S-N and Findley tensor histories. The wrapper retains the existing mechanical
trajectory, accepted-pose observer, work and rollback owners; it supplies no
new integrator, averaged material law or additional mechanical energy.

### Native sustained finite-solid friction heat with source feedback

`friction_heat::finite_solid::propose_finite_solid_contact_heat` closes endpoint
body/common/declared-mean source temperatures, the original sourced sustained
sphere/finite-solid trajectory and the existing finite enthalpy network together.
Each original feature retains its source binding, ordered materials, caloric
objects and bounded reference thermal patch. Only unique loaded, nonopening,
passive midpoint tangential impulse work is eligible for declared heat
conversion. Exported tangential work, signed normal/contact remainder and native
integration defect stay separate in the combined first-law balance; no residual
or normal loss is deposited as friction heat.

Every trial restarts the accepted solid, rigid bodies and original source history.
Exact repeated full source queries reuse one bounded mechanical result while
transport refreshes; common interface temperature alone is never a cache key.
Source, geometry, authority, thermal, clock or balance refusal returns no partial
accepted state. Fixed witness neighborhoods and ordered face normals are checked
at native source evaluations and the final endpoint. This does not add a swept
thermal-footprint certificate. The joint candidate retains original synchronous
endpoint times and all mechanical/history/thermal reports.

This is fixed-reference-patch interface temperature feedback with frozen bulk
hyperelasticity. Pressure-dependent source coordinates, moving footprints, bulk
expansion and impact heating remain outside the model. The existing trajectory
and thermal interval still need time refinement; nonlinear convergence is not
time-accuracy evidence. The native regression definition in
`tests/finite_solid/contact/step/sustained/source/heat.rs` compares an independent
body-temperature-dependent friction/impulse/caloric solution, zero conversion,
source path, fixed-patch refusal and late source rollback. Isolated execution of
the production work-account function with metadata fixtures validates its heat,
export and remainder partition, not the complete native multiphysics solve.

## Cylindrical and spherical moisture diffusion

`moisture::radial::MoistureRadial` advances bound-water diffusion on nonuniform
cylindrical or spherical shells using their actual dry-skeleton storage volumes
and exact radial face resistances. Original source-integrated Kirchhoff laws
retain moisture-dependent and degenerate dry-front transport; finite sourced
surface exchange uses each boundary's curved area. Conservative backward Euler
uses a reusable native tridiagonal workspace, bounded nonlinear iteration and
independent endpoint residual/water-balance checks before publication. A solid
center is sealed. This is a fixed, homogeneous, isothermal radial reduction;
axial/angular gradients, arbitrary tensor orientation and swelling feedback
require their respective spatial/coupling owners.
Analytic cylindrical/spherical steady potentials, a spherical transient
eigenmode, nonlinear dry-front balance and finite surface flux define the focused
verification in `tests/moisture/radial.rs`; they are not measured validation or
a continuum-error certificate for an arbitrary mesh or time step.

### Prescribed temperature histories in radial diffusion

`MoistureRadial::at_prescribed_temperature` re-evaluates the original pinned
multiaxis diffusivity at every occupied shell content while retaining each
shell's dry mass, water, geometry, accepted time and original density/reference
identity. A supported enclosing temperature/moisture box cannot override the
actual source triangulation. Degenerate dry fronts retain the original
Kirchhoff path, and an unsupported cell or cancellation refuses the whole event.

Radial backward Euler, SDIRK2 and finite-atmosphere admission consume the current
transport temperature. Source-qualified sorption and finite surface laws must
match that temperature; the cold density datum does not silently drive a warmed
boundary. Finite solutions already consume the same current inventory context.
Transport intervals must end at externally prescribed temperature events.
This path supplies no heat balance, thermal expansion or temperature evolution;
the native mechanical owner still refuses an unqualified thermal strain.

### Source-partitioned coatings and curved material layers

`moisture::radial::layered::{MoistureRadialLayer, MoistureLayeredRadial}` joins
contiguous native cylindrical or spherical shells while retaining each original
dry inventory, nonlinear diffusivity and same-card monotone sorption source.
The shared layered activity/storage/Picard kernel serves both slabs and curved
stacks; each geometry supplies its actual face paths. Radial layers use exact
cylindrical/spherical resistance integrals, unequal shell masses and the actual
curved interface and exterior areas. Existing slab deformation-metric feedback
continues to resolve the proposed state at each nonlinear assembly.

Ideal contact equalizes activity while allowing material-specific dry-basis
moisture jumps. Explicit finite activity conductance permits an interface
activity jump with one equal/opposite water transfer. A zero-area center stays
sealed. Source/context, geometry, temperature, clock and support mismatches
refuse; native source/range/work/residual/water checks precede publication of
the entire stack. Geometry and temperature are held fixed; this rung does not
add hysteresis, moving interfaces, sorption heat or liquid capillary flow.
Five authored physical regressions define curved steady partition/transfer,
unequal-density closed transients, homogeneous-owner reduction, nonlinear dry
core/coating uptake and transactional geometry/source/admission refusals.

### Finite solution activity on radial geometry

`FiniteMoistureSolution::advance_radial` and `advance_shared_radial` connect
native cylindrical and spherical diffusion to the original finite aqueous
inventory. Each shell contributes its actual dry-mass-weighted water change;
one common activity root retains all nonvolatile solutes and any admitted
neutral/ionic or kinetic phase owner. Curved half-cell resistance, nonlinear
Kirchhoff transport and sourced surface exchange remain in the radial solver.
Stable specimen IDs fix signed transfer summation, and all fields and the
solution publish together after the original water and source checks.

The corresponding `advance_radial_sdirk2` and `advance_shared_radial_sdirk2`
solve activity at both physical stages. Both start from the original solution
inventory; the final report already contains the full weighted water transfer.
Kinetic crystallization explicitly refuses this second-order path because its
endpoint-frozen phase law supplies no RK stage equations. Equilibrium crystal
stocks remain finite and may buffer the activity. Smooth-branch temporal
accuracy requires refinement; phase events need their own resolution.

The coupling retains prescribed common temperature, equal pressure and an
equilibrated vapor path with negligible vapor storage. It does not reinterpret
vapor transfer coefficients as liquid-immersion permeability or move solutes
into the solid. Five focused authored regressions define analytic curved uptake
and release, concentration feedback, unequal-shell shared conservation/order,
second-order convergence, finite crystal buffering and transactional refusal.

### Second-order heterogeneous water storage and interface transfer

`MoistureLayeredSlab::advance_sdirk2` and
`MoistureLayeredRadial::advance_sdirk2` use Alexander SDIRK2 through the same
native layered storage/activity kernel. Both stages retain each original
nonlinear sorption and Kirchhoff diffusivity, actual dry-cell inventory and
explicit interface/exterior law. Stage two starts from the original water
storage and receives the actual first-stage face rates as physical water loads;
neither an activity increment nor a frozen capacity substitutes for stored water.

`MoistureLayeredSdirk2Report` separates the first diagnostic stage from the
complete endpoint. Exterior and material-interface water transfers in the
endpoint already include both RK weights; surface activities and source
receipts describe the final stage. Shared source preparation, native
tridiagonal iteration and admission remain in use for both spatial geometries.
The existing backward-Euler and deforming-slab paths retain their original
behavior. This SDIRK entry point holds geometry, temperature and reservoir
activities fixed. A nonmonotone stage returns typed `ActivityStageRange`,
which the existing duration controller may retry at a shorter interval.
No clipping, first-order replacement or finite-reservoir SDIRK claim is implied.

### Finite vapor exchange with coated radial bodies

`moisture::radial::layered::atmosphere` couples single or shared coated spheres
and cylinders to the original finite-vapor inventory. Each humidity trial
advances the complete original body with every material-specific sorption,
nonlinear diffusion and interface law intact. Compensated local moisture
increments, weighted by the actual unequal dry shell masses, supply one gas
water account. Canonical specimen IDs give deterministic shared summation and
result order. An interior material interface is never an atmospheric port.

Exterior source receipts resolve the pinned end-material sorption claim at the
accepted solid-side surface activity; explicit film resistance can separate it
from gas humidity. Authored exterior seals remain sealed while internal
redistribution continues. Source/context, common clock/temperature, zero-area
center and total cell budget admission precede the original water/carrier/total
mass checks and atomic publication. Temperature, free gas volume and reference
geometry stay fixed; this path supplies no condensation, carrier permeation,
hysteresis or thermal feedback.

### Finite atmosphere exchange on radial geometry

`FiniteMoistureAtmosphere::advance_radial`, `advance_shared_radial` and
`advance_paired_radial` couple the original radial fluxes to the existing finite
gas inventory owners. The actual implicit endpoint RH drives sourced sorption
and surface transfer; every probe restarts the original specimen. Compensated
local water increments use each shell's actual dry mass. Shared chambers debit
one gas inventory for all stable-ID-ordered specimens; paired chambers retain
separate signed inner/outer transfers even when total solid water is unchanged.
The original atmosphere owners resolve changed composition and pressure and
check carrier, water and total mass before all endpoints publish together.

Native solid adapters resolve compatible deformation and wet vibration inertia
after the accepted single or paired atmosphere candidate. A late constitutive
or mechanical refusal preserves every original gas and solid input. Gas
pressure is not automatically imposed as a mechanical traction. Temperature and
gas volumes are prescribed; condensation, carrier permeation and sorption heat
remain separate physics. Source-supported RH brackets and explicit work budgets
remain caller inputs. Seven radial gas regressions cover analytic exchange,
finite surface resistance, opposing chamber flow, unequal shared inventories,
dry-front wetting, sealed redistribution and atomic refusals; a native solid
regression adds finite-vapor shape/inertia and late mechanical refusal.

### Second-order radial transients

`MoistureRadial::advance_sdirk2` uses the same two-stage L-stable SDIRK method
as the existing slab and Cartesian owners. Both stages reuse the original radial
Kirchhoff law, unequal shell storage, curved surface exchange and bounded native
tridiagonal nonlinear solve. The first-stage explicit load is assembled from
actual conservative face rates, never rounded water differences. The endpoint
water account uses the RK-weighted sum of both stage boundary transfers.
`MoistureRadialSolid::advance_sdirk2` resolves compatible native mechanics after
the complete admitted endpoint, preserving atomic transport/mechanics publication.

L-stability does not guarantee monotonicity: a nonphysical stage returns the
existing typed `StageRange` refusal and requires a shorter time interval. It is
not clipped or replaced by backward Euler. Conditions are constant within the
call; callers split intervals at environmental events. Independent unequal-shell
linear modes, a nonlinear dry-front ODE, finite surface resistance and weighted
water balances define six focused native regressions. A converged stage is not
a global temporal-error estimate.

### Second-order and adaptive finite-vapor radial transport

`advance_radial_sdirk2` and `advance_shared_radial_sdirk2` on the original
finite-atmosphere owner solve the actual gas humidity at each implicit stage.
First-stage face rates supply the explicit RK load; the second stage balances
the original gas inventory and original solid storage against the full weighted
water transfer. Diagnostic intermediate gas/solid states are never debited twice.
The prepared Kirchhoff source interval retains the first stage's admitted range
for later humidity probes, while the algebraic reference shift and actual
active-face physical range checks remain unchanged.

`advance_radial_adaptive_to` and `advance_radial_sdirk2_adaptive_to` control
a complete finite-air duration using the existing temporal driver. Acceptance
checks every shell's moisture, noncancelling local water error, gas vapor mass
and RH. Accepted fine-step transfers alone enter the fixed whole-window water,
carrier and total-mass account. Source/bracket limits remain explicit; a failed
duration publishes neither inventory. Matching native-solid methods resolve the
compatible final shape and wet vibration pencil after these gates, with late
mechanical refusal preserving all inputs. Second-order endpoints also
support multiple specimens sharing a chamber. Adaptive finite-air durations
currently couple one specimen and one chamber; paired-chamber exchange retains
its backward-Euler path. Eight native regressions define
independent coupled exponential limits, nonlinear curved-film fluxes, stage
source support, unequal-shell/gas accuracy and whole-call rollback. Their error
estimates are local, not global continuum or experimental-validation claims.

### Accuracy-controlled radial durations

`MoistureRadial::advance_adaptive_to` and `advance_sdirk2_adaptive_to` reuse
the existing whole-duration temporal controller around the original backward
Euler and SDIRK2 endpoints. The shared borrowed inventory now explicitly owns
either uniform or per-cell dry masses. Generic water, caloric, reservoir and
solid accounts consume that physical measure, so nonuniform radial shells are
never interpreted as equal storage cells.

Each comparison resolves both componentwise moisture and the L1 sum of actual
shell water discrepancies; opposite internal errors cannot cancel in a closed
specimen. Relative scales use transported changes. For radial calls,
`reservoir_water_absolute_tolerance_kg` also controls this shell-water norm.
Only accepted fine-step transfers enter one fixed whole-duration water budget.
The existing controller charges rejected/stage work, bounds retries, handles
typed physical stage refusals and publishes no successful prefix after a later
failure. Matching native-solid convenience methods resolve the mechanical
endpoint after the complete admitted duration. The estimate remains local and
asymptotic; boundaries are held within the duration and environmental events
require explicit cuts. Six focused native regressions use unequal-shell modes,
nonlinear dry-front solutions, an independent mass-weighted error ratio and
work/cancellation rollback.

### Adaptive radial diffusion with a finite solution

`FiniteMoistureSolution::advance_radial_adaptive_to` and
`advance_radial_sdirk2_adaptive_to` use the original bounded duration controller
with complete native finite-solution backward-Euler or SDIRK2 endpoints.
Every trial solves actual solvent-volume/activity feedback with the original
curved transport paths. Error estimates compare each shell content, the
dry-mass-weighted L1 shell water discrepancy, finite solvent mass and solution
activity. Opposing shell errors cannot cancel, and pure-solvent uptake remains
controlled even when its activity stays one.

Only accepted fine endpoint exchanges enter the cumulative boundary account.
The original-to-final shell/face balance and combined solvent/shell water balance
must each meet one fixed whole-duration tolerance. Both SDIRK stages count
against the work budget, and any late failure discards the entire candidate.
Adaptive radial solid wrappers resolve shape, stress and wet inertia once from
the accepted final moisture, then publish the complete dependent result together.

The existing prescribed-temperature, equilibrated vapor bridge remains in use;
these methods do not add liquid immersion or osmotic traction. Native kinetic
phase evolution remains first order under BE and is refused by SDIRK2. Individual
phase-component errors and phase heat are not separately controlled. Local
step-doubling estimates do not promise a global accuracy bound.

### Radial sorption energy and finite-solution mechanics

`FiniteMoistureAtmosphere::advance_radial_with_sorption_energy` and its SDIRK2
counterpart integrate the original differential isosteric heat over each actual
shell's dry-mass-weighted loading change. The existing caloric owner closes
original gas plus bound-water internal energy against the heat supplied by the
external isothermal thermostat. Loading-dependent internal redistribution can
require heat even when no water leaves the specimen. The intermediate SDIRK
stage is never added as a separate energy debit. Source work and energy-balance
failure withhold the complete transport/energy candidate. This is global heat
duty, not a reconstructed local heat equation or temperature feedback.

`MoistureRadialSolid::advance_finite_solution` and its SDIRK2 counterpart resolve
compatible native deformation, stress and the wet prestressed vibration pencil
from the exact finite-solution diffusion endpoint. Late mechanical/source
refusal withholds both water inventories and their dependent solid state.
Osmotic pressure remains diagnostic; no immersion flow, inferred traction or
elastic-work feedback into chemical potential is added. Five authored tests
define unequal-shell heat integrals, sealed redistribution, late caloric
refusal, native affine swelling/wet inertia and late mechanical refusal.

### Radial diffusion driving native solid mechanics

`moisture::solid::MoistureRadialSolid` maps every original cylindrical or
annular water shell into a conforming three-dimensional tetrahedral mesh. Local
transported moisture selects the original elastic law and integrated directional
swelling source; the existing native solid owner resolves compatible displacement,
stress, reactions, deformed geometry and the signed prestressed K_e+K_g/M pencil.
A nonuniform water profile therefore produces compatible stress and local wet
inertia rather than independent shell expansions. Material axes are fixed in the
xyz mesh frame; radial/tangential grain rotation is not inferred.

Each shell's exact curved dry inventory is distributed over its actual polygonal
tetrahedral volumes. The user explicitly bounds the curved/polygonal volume
discrepancy; the corresponding inertia-density correction vanishes under angular
refinement. Axial, circumferential and radial mechanics still require ordinary
mesh convergence. The water metric remains on the original reference chart:
small deformation does not imply stress-assisted or spatially swelling diffusion.
This owner admits cylinders and annuli; spheres require a separate native mesh
owner. Source-domain, geometric, work, cancellation and late mechanical refusals
leave the original moisture state untouched. Geometry tests independently check
shell volumes, seam conformity and refinement; native integration tests define
affine free swelling, gradient-induced stress, exact conserved wet inertia and a
one-DOF generalized modal response using authored verification laws.

## Coating strain about a supplied static modal offset

`MaterialShellImage::prepare_coating_strain_with_offset` and the corresponding
assembly-image method admit `PreparedShellCoatingStrain::evaluate_with_static_offset`.
The supplied real modal offset and coherent complex modal amplitude use the same
original source-bound mode basis, facet geometry and coating sections. Each
reported film bound covers the full physical layer, complete triangle and
complete harmonic cycle about that offset; independently acceptable offset and
oscillation are not assumed acceptable when combined.

Preparation admits the extra real field, modal accumulation, storage and bounded
static-plus-harmonic maximization work before observation. A successful query
retains complete film coverage; invalid input, unsupported preparation and
cancellation do not consume the successful-query allowance. The reported
sufficient amplitude scale multiplies both supplied fields together. It is not
an oscillation-only allowance at fixed preload. The observer does not infer
thermal/cure eigenstrain or solve nonlinear prestress equilibrium.

Three integration regressions cover combined-limit crossing, phase invariance,
native snapshot containment and transactional admission/refusal. Their complete
source-material/modal fixture remains unexecuted in the available dependency
environment. The underlying native shell bound passed five focused numerical
tests, including complete-domain containment, degeneracies and extreme scaling.

## Original thermal hardening surfaces on reciprocal caloric solids

`FiniteThermoplasticMaterial::try_new_with_tabulated_hardening` binds an
explicit tabulated finite thermal model and its original strain/temperature
observations to the exact same-card caloric chart. The complete original
hardening source and every selected observation remain inspectable through the
resolved family. Original density at the stress-free reference temperature still
sets reference mass; expansion and density metadata cannot change that inventory.
Both scalar and tabulated constructors use the same caloric source/context and
positive-total-heat-capacity admission.

`FiniteThermoplasticSolid::try_new_with_tabulated_hardening` prepares these
families through the same cell mapping and native tetrahedral constructor.
Explicit claim selection and a per-source preparation budget are supplied by
the caller. Existing equilibrium, retained inertia, accepted-state transactions,
physical interval and reciprocal heat iteration remain the only step owners.
The prepared material can also enter the existing native finite thermal sheet
and `FiniteThermoplasticSheetCoupling`; no second heat solver is introduced.

Heating can cross a measured hardening knot and advance plastic strain before
melting. Its hardening-entropy change is retained in signed reversible caloric
exchange, independently from positive physical plastic heat. Numerical endpoint
work defects remain separate. The same accepted plastic metric, accumulations,
temperature and caloric inventory start every feedback proposal.

Two coupled regressions author an independent coaxial table/entropy oracle,
solve a free displacement with temperature, cross a hardening knot on heating,
and cool/unload to retained permanent shape. They also check original receipts,
foreign caloric rejection, bounded temperature/strain support, complete
nonvirgin state retention after refusal, and successful retry. Their independent
scalar energy identities and Rust syntax were checked; native regression
execution remains unverified in the resource-constrained session.

## Spatial pressure in moving partial-slip friction and heat

`friction_heat::moving::partial_slip::MovingPartialSlip::new_with_pressure_profile`
accepts `MovingPressureProfile` from the original normal-patch authority. Its
positive, piecewise-constant pressure cells exactly cover the unchanged strip
footprint and integrate to the original normal load, to arithmetic precision.
The original patch identity, source, authority, load and geometry must match the
retained checkpoint. Gaps, overlap, foreign receipts and load mismatch refuse;
there is no extrapolation, load renormalization or inferred pressure shape.

Pressure stays in footprint coordinates while the two native Jenkins-region
histories travel with surface A material. Exact intersections split parcels at
thermal/material and pressure boundaries without averaging strains or resetting
kinetic branches. Sorted interval cursors keep this splitting proportional to
the actual intersections. Each source callback receives its local pressure,
actual regional normal force and original caloric cell. The original return
map still owns local force, microslip, stored energy and physical heat; its
reference-normalized response is multiplied by the actual regional area share.

Stiffness, region decomposition and the storage-source identity remain fixed.
A lower local pressure may relax retained shear energy through the original
return map, even at zero additional imposed slip. Only its physical plastic
loss heats the finite cells. Numerical return-map loss and elastic energy
leaving the footprint remain separate accounts; no parameter work or normal
work is invented or thermalized. Raising pressure preserves loaded history.

Six new physical regressions cover uniform-profile equivalence, early edge
microslip with a sticking core and spatial heat, transport across pressure
boundaries with independent release/outflow balances, local pressure plus
endpoint-temperature source feedback, exact source/load/coverage admission,
and bounded transactional refusal/retry. They accompany the existing moving
brush tests for work, frame covariance, history and temporal refinement.
Pressure-space, thermal-space and time refinement remain independent. The
existing pure-shear, fixed-load, fixed-width and first-order transport/contact
composition is retained; zero-pressure gaps require a separate normal/contact
topology, and no resolved continuum traction or measured calibration is inferred.

### Evolving normal load and pressure without resetting moving shear

`MovingPartialSlip::step_candidate_with_pressure_profile` accepts the original
normal owner's next positive load and piecewise-constant pressure field on the
same footprint. The native partial-slip checkpoint admission retains exact
patch/model/source/authority and ordered-interface lineage. Semiaxes, pressure
moment, shear stiffness and its source stay fixed. New pressure intervals must
cover the unchanged physical footprint and integrate to the supplied new load;
no coefficient or pressure rescaling repairs an inconsistent receipt.

Transported material strains and kinetic branches enter the original regional
return map unchanged. Unloading can release stored shear at zero imposed slip;
only its native plastic dissipation heats the original finite thermal bodies.
Numerical return loss and elastic outflow remain separate. Increasing load
retains the surviving strains without inventing preload. Local source queries
receive the new regional pressures together with both actual caloric cells and
their solved interface temperatures. Pressure, all material histories and heat
publish atomically. `normal_patch()` exposes the accepted normal receipt, while
`reference()` retains the original construction checkpoint; ordinary subsequent
steps continue with the accepted pressure automatically.

This is a first-order endpoint-pressure split. Normal-path refinement is
independent of pressure-cell and thermal-mesh refinement. Normal indentation,
normal storage and actuator work belong to the supplied normal owner; this
brush neither infers nor thermalizes them. Separation, footprint growth,
pressure-dependent shear stiffness and a coupled normal equilibrium solve
remain outside this path.

Five G1/G3/G4 regressions verify analytic unloading/reloading, local pressure
redistribution at constant total load, convergence of physical heat under
16/32/64 load subdivisions, lineage/geometry/load refusal before source work,
and atomic late refusal, cancellation, budget exhaustion and retry. Together
with the preceding moving-memory, pressure and two-temperature regressions,
all 22 focused native tests pass. No measured calibration is inferred.

## CST/DKT hereditary material memory on the original modal basis

`plate_viscosity::shell::relaxation::ModalShellRelaxation` projects per-facet
passive excess A/B/D stiffness spectra through the original CST/DKT strain
owner. Its modes must be orthonormal in the original sparse mass, including
consistent director/translation terms; the instantaneous-viscosity mass
check is shared unchanged. The equilibrium pencil and retained coordinates
remain fixed. Source applicability, matching mesh/basis identity and the
interpretation of supplied tensors as excess above relaxed stiffness remain
explicit caller obligations.

Exactly equal relaxation times aggregate before the shared compensated PSD
factor sees each original dimensional modal Gram matrix. Different times,
resolved weak directions and all modal cross terms remain present. Each
resulting row c [1/s] carries the existing Maxwell storage `(c.q-z)^2/2`
and irreversible power `(c.q-z)^2/tau`; numerical drilling acquires no
constitutive state. Spatial work, complete factorization, branch embedding
and material-state counts have explicit checked budgets. No dense projection
or factorization occurs during playback.

`branches` embeds these rows into existing fs-phs interleaved mechanical
coordinates. `attach_to` composes them through the original
`CoupledModalSystem::with_additional_material_relaxation`: relaxed history
sets new z=c.q, unrelaxed history sets new z=0, and earlier material histories
remain unchanged. The target component must have zero intrinsic modal
damping; separately selected viscous, support, glue and radiation channels
can coexist. Existing cold-construction, instantaneous-frequency, energy and
augmented-state admissions remain binding. The existing exact propagator and
independent material-loss integral are the only time and heat owners.

Three focused native tests pass: full spatial/modal/channel storage for both
mass models and every retained time, branch embedding/refusal/zero cases,
and actual fs-phs exact propagation against independent Maxwell equations
with separately integrated heat and input work. Both initial histories pass
the 16/32/64-subdivision RK4 convergence checks and energy closure without
tolerance changes. The additional full coupled-attachment regression covers
component offsets, prior histories, initialization, loss counting and cold
binding; it is authored but unexecuted in the available dependency closure.

## Thermal porous follower pressure through enthalpy feedback

`finite_plastic::porous::FiniteThermalPorousSolid::advance_with_pressure` and
`advance_with_pressure_and_coupling` carry the original finite-solid follower
pressure through every native thermal proposal. The existing bounded Picard or
IQN-ILS iteration repeatedly solves from the same accepted deformation,
temperature, plastic metric and optional velocity. Original material/caloric
objects, reference volume, reference density and cell mass remain bound.

The pressure-aware candidate retains the original material report plus the
native `FinitePorousPressureEnergy`. Uniform pressure uses `-p delta(V)`;
per-face pressure retains the straight-displacement-segment work convention,
and spatial pressure fields retain their native potentials. The coupled energy
identity includes pressure path work once and the complete signed external-work
defect (material endpoint defect plus pressure endpoint correction). The
material's endpoint defect remains independently visible. Only original
physical yield heat and reciprocal thermal exchange are deposited into the
caloric timeline. Pressure correction and backward-Euler numerical kinetic loss
never become heat. Existing unpressured calls retain their report type and
arithmetic path; the shared feedback owner subtracts additional pressure work
only when it is nonzero. That pressure-aware identity uses compensated sums
for opposing boundary work and keeps kinetic change and numerical loss separate
inside the reduction.

Focused G1/G3 regressions extend `finite_porous_enthalpy`: a free displacement
and its temperature/plastic volume are recovered from an independent
hydrostatic pressure and caloric root with both coupling methods; prescribed
motion checks exact volume work, source reference mass, numerical-loss
separation, zero-pressure compatibility, late observer/iteration rollback and
successful retry. These are authored constitutive fixtures, not experimental
validation. The supplied thermal network still uses reference conductances; no
deformed heat-transfer area, wetting law, exterior-fluid evolution, pore-fluid
transport or global time-accuracy certificate is inferred.

Native regression execution remains unverified on the exhausted shared build host.

## Physical impedance back-reaction in nonlinear harmonic balance

`phs_orbit::impedance::PeriodicPortImpedance` admits the complete original
physical port impedance at DC and every supplied harmonic of one fixed clock.
The source matrices use `exp(-i omega t)`; preparation conjugates them for
fs-orbit's positive-frequency Fourier chart. DC is explicitly finite and real.
The shared fs-modal Hermitian PSD admission preserves all mutual and
nonreciprocal entries. It establishes pointwise numerical admission only;
physical source validity, unsampled bands and causality stay with the caller.

`PreparedPhsImpedance` augments the original pHS state with physical flow
unknowns, using `P(s)=[[sI,G Z(s)],[0,I]]` and the nonlinear island
`[(J-R) grad H+G f, G^T grad H]`. The load therefore changes the solved motion;
it is not a receiver observation or an effective damping substituted for the
source. Original nonlinear storage, physical port ordering and every Maxwell
memory coordinate remain intact. Complete native/caller storage tangents
differentiate the same original law; exact quadratic tangents are cached.

Fixed-order and adaptive scaled seeded solves use the existing fs-orbit
Newton/AFT/refinement owners. Seeds and scales include the appended flows;
their algebraic equations have flow units. The drive must fit the initial
order and the table must cover every potentially retained order before source
evaluation. No omitted matrix defaults to zero. The independent omitted-band
assessment requires no load beyond retained order because those state/flow
coefficients are exactly zero. No time-domain RHS, shooting capability,
infinite-tail bound or stability claim is inferred from the frequency table.

The returned state orbit retains the original coordinates for acoustic
observation. Its full-balance diagnostics also include the flow constraints.
The independent cycle audit reevaluates original storage and ports, separately
reporting authored source work, native resistance loss, absorbed boundary work
and the RMS difference between nonlinear physical flow and retained flow.
Boundary work is not inferred to be heat or radiation, and a sample ladder
remains necessary for nonlinear work convergence. Edited/foreign boundary
efforts, unsupported sampling, absent bands, refusing sources/tangents,
resource limits and cancellation refuse without an unloaded fallback.

Five G1/G3 regression definitions exercise the actual baffled-piston owner,
native nonlinear Duffing/Maxwell storage, nonlinear port harmonics, full mutual
nonreciprocal loads, DC work and power-conjugate unit changes, plus admission
and cancellation. Rust syntax/formatting were checked; native execution remains
unverified in the current dependency/resource environment.

## Exact distributed viscothermal memory

`WideTubeLoss::new_with_sampling` and
`with_regional_viscothermal_sections_options` explicitly select exact viscous
and thermal relaxation through `ViscothermalTimeOptions`. The original gas
transport coefficients, eight-arm Foster spectrum, wide-tube applicability and
spatial splitting remain unchanged. Both generated series impedance and dual
thermal admittance carry the same rule into the original acoustic network.
Exact viscous branches use held flow/mean pressure; exact thermal branches use
held pressure/mean flow. Their native owners integrate independent positive
resistor heat and retain the same wave-plus-material work balance.

The 5% complex/real-loss fit gate and 0.03 section-scattering gate use the actual
selected sampled port transfer, conjugated to the duct frequency convention.
They retain the same analogue target and can refuse an unresolved exact clock.
Integer and finite-band fractional propagation remain independent explicit
choices; original lengths, section gases, allpass phase and every delay enter
the cold scattering check. Existing constructors retain midpoint defaults.
`FlueVoice::from_network_with_time_options` carries the selection through the
original jet and optional sourced walls, while `TwoSidedAperture::with_time_options`
selects each original reed tract independently without replacing its mechanics.
Native synthetic tests compare generated load responses, actual integer and
fractional network transmission, independent storage/work, geometry, default
compatibility and existing refusals. This remains a first-order wide-boundary-
layer model, not all-shear Bessel physics or measured instrument qualification.


## Exact source-bound Maxwell wall memory

`WallRelaxation::with_network_sampling(ExactHeldFlow)` carries an explicit time
representation from the original generalized-Maxwell source into its acoustic
shunt. The shared material-file selector `maxwell-exact` takes the same modulus
and time property pairs as `maxwell`, and reaches both reed and flue network
walls. The selected material snapshot, claims, applicability, geometry, mass,
equilibrium stiffness and all positive source arms remain unchanged.

Exact branches use the solved held wall flow, continuous mean stress pressure
in the reciprocal port equation, and independently integrated physical heat.
The network may therefore admit poles faster than its step without midpoint
ringing. Its base inertial membrane still uses midpoint and must satisfy the
independent instantaneous-frequency limit. `MAX_DT_OVER_TAU` remains part of
the source declaration and still governs distributed FE midpoint memory; DG
retains its existing coupled propagator. The original `maxwell` selection is
unchanged. Pinned source refresh preserves the selected sampling, and a live
history transfer refuses a change of time rule. Focused source-card, clock and
full wave-network tests compare analytic branch decay, mean pressure, membrane
momentum, spring storage and resistor heat. No within-step flow resolution,
nonlocal shell coupling or full voice qualification is implied.

## Pair-resolved material and interface feedback in moving partial slip

The existing `friction_heat::moving::partial_slip` owner intersects each
transported A-material parcel with both original surface cell grids at the
proposed endpoint, as well as the original supplied pressure intervals.
`MovingSlipRegion::surface_b_cell` retains the opposing caloric identity;
`endpoint_interface_index` identifies that exact A/B pair in the callback's
swept patches and solved interface-temperature coordinates. A constitutive
source can therefore depend on either original body temperature and that
pair's massless interface temperature. One large A cell no longer forces the
same traction law across differently heated or heterogeneous B cells. No
averaged interface temperature, new caloric store or body heat fraction is
introduced.

Monotone cursors merge the actual material, A-cell, B-cell and pressure
boundaries rather than constructing a Cartesian product of the surfaces.
Original normal-load coverage is independently checked. Every resulting
region copies its surviving native elastic strains and kinetic flags exactly;
crossing an opposing cell boundary neither resets nor averages that memory.
Only the original law's static/kinetic coefficients may vary. The source,
storage law, pressure authority and generalized-work owner remain unchanged.
The original region limit includes all newly resolved intersections and
refuses before a constitutive source call if they do not fit.

Constitutive coefficients are endpoint source values over the existing
first-order transport/return-map interval. Regional physical plastic heat still
follows the complete swept footprint, including cell pairs visited before the
endpoint. Assigning a source to its endpoint pair does not reassign that
history of heat deposition. The original finite two-body network determines
signed body heat and each interface temperature together; numerical return
loss and elastic outflow retain their separate mechanical accounts. Time
refinement across heterogeneous material boundaries remains independent of
thermal-cell and pressure-field refinement. Transverse transport, torsion,
changing normal load and native body momentum are still outside this strip
owner.

Five focused regressions check independent unequal-traction and swept-heat
integrals, coupled opposing-body/interface feedback against analytically
condensed finite-caloric equations, exact strain/kinetic history through a
moving B boundary, constant-source mechanical invariance under opposing-grid
refinement, and bounded refusal with retry. Existing pressure tests retain
their independent load, storage, physical heat, numerical loss and outflow
oracles while allowing the actual opposing grid to subdivide each pressure
region. These synthetic constitutive tests establish no measured material
calibration or whole-system convergence claim.

## Accelerated and reversing moving shear with finite heat

`MovingPartialSlip::advance_quadratic_candidate` carries the original brush,
pressure field and two finite thermal bodies through one prescribed quadratic
placement path. It splits at every relative-slip reversal and every reversal
of footprint residence on either material, then at the explicit time limit.
Returning to the starting placement therefore retains shear hysteresis,
plastic heat and expelled elastic storage. Entering material remains unloaded;
a returning footprint does not resurrect an earlier departed checkpoint.

The existing quadratic geometry owner computes complete swept overlaps and
the actual mean force moment. Full-path coverage includes extrema between
the three supplied placements. The local endpoint source can read actual
incoming/outgoing slip velocities and both native cell temperatures. Each
regional integrated plastic heat is held as constant mean power over its
subinterval; the original first-order transport/return-map and thermal split
still needs time refinement. Exact overlap geometry does not certify an
unresolved plastic-power history. The current normal pressure and footprint
remain fixed over this path; the separate load-path entry retains its scope.

All subintervals share step and nonlinear-trial caps and a cumulative energy
budget. Per-interval work owners derive distinct IDs from the caller's parent
interval; every physical heat, numerical return loss, elastic outflow and
external thermal term has one account. No source payload or mechanical/thermal
prefix publishes if a later subinterval refuses or cancellation arrives.

Six G1/G3/G4 regressions cover a closed shear cycle with independently derived
0.03875 J plastic heat despite zero net displacement; decreasing numerical
loss without heating it; actual quadratic residence and convergence to the
continuous power/residence integral; reversible common motion and its mean
shear moment; non-resurrecting outflow memory; and hidden excursions, bounded
work, late source failure and exact retry. These are numerical synthetic-law
tests, not measured friction calibration or a native body-motion solver.

## Frequency sensitivity of an original forced periodic response

`PreparedPhsPortSidebands::forced_frequency_tangent` factors at zero carrier;
`FactoredPhsPortSidebands::forced_frequency_tangent` reuses that retained factor.
For fixed physical source parameters and authored input phasors, it solves
`(omega*d/dtheta - Df) z = -omega*x_theta` over every original state, including
explicit memory, with both signed Fourier indices. The physical output is
`G^T Hxx(x_base) z`, independently observed with the existing refinement gates.
The result reports `z = omega*dX/domega` and the corresponding output derivative,
retaining original physical units; division by the base frequency gives the
derivative per rad/s. Phase is held fixed, not physical time.

Original base balance, tangent sampling, retained/omitted equations, source
admission, cancellation and numerical budgets remain active. Autonomous sources,
nonzero carriers and omitted active base harmonics refuse. No second Jacobian,
finite-difference frequency solve, reduced resonance fit or fold regularization
is introduced. Delayed/dispersive forcing needs its original state realization;
fixed authored phasors do not supply that missing frequency dependence.

Four integration regressions specify analytic mass/inertia/Maxwell-memory
response through multiple harmonics; an original nonlinear storage response
compared with independently re-solved harmonic-balance orbits; source/resolution
and cancellation refusals; and an autonomous-clock refusal. They are authored
but unexecuted without the complete constellation dependency closure. Targeted
formatting and API review passed; neither full integration nor continuum,
branch-selection, Floquet or forward-error qualification is inferred.

## Independent MITC material spectra in the shared memory runtime

`ModalShellRelaxation::from_mitc_shell` carries arbitrary passive per-facet
A/B/D and transverse-shear increments into the same complete Maxwell memory
owner used by CST/DKT. The original sparse mass admits every input basis Gram
entry, including consistent or lumped reference-surface inertia and its first
mass moment. Original support frames and elastic bubble condensation remain
owned by the retained MITC model. Eigenvalues and mode shapes are not refitted.

Every exactly equal material time is factored from its complete projected
stiffness in original units. Its channel has storage `(c.q-z)^2/2` and physical
heat rate `(c.q-z)^2/tau`; all resolved tensor directions and modal cross terms
remain. The existing `branches` bridge and cold `attach_to` path supply the
original fs-phs/coupled runtime, including its relaxed/unrelaxed initial
history and prohibition on double-counted intrinsic scalar damping. No new
time integrator or material viscosity surrogate is introduced.

The cold MITC adapter admits up to 1024 original modes and 8192 material
coordinates subject to checked full projection, mass-Gram and factor work.
The selected runtime retains its own tighter state, work, frequency and
initial-energy ceilings. Source bands and the interpretation as increments
above relaxed equilibrium remain caller declarations. Tests compare full
source strain storage, original mass choices, support-frame and basis
covariance, both initial histories, and actual fs-phs motion, independently
integrated heat and input work against refined continuous Maxwell equations.

### Original causal rational boundary states

`phs_orbit::rational::PreparedPhsMimo` attaches the original admitted
`fs_vfit::mimo::load::PreparedMimoLoad` to an original nonlinear pHS through a
caller-authored state-by-boundary-port map `S`. The complete state is
`[x_host, z_source]`, in the original coordinates of both owners. Its equations
are `y=S^T grad H`, `z'=A z+B y`, `e_boundary=C z+D y`, and
`x'=(J-R) grad H+G f-S e_boundary`. The same signed map carries source flow
and opposing host reaction. External forcing ports remain independent of the
boundary map. Full mutual matrices, singular or zero direct terms, real poles,
conjugate-pair coordinates and source states invisible at the boundary survive.
There is no fit, pole condensation, direct-term inverse or new time integrator.
The original continuous owner supplies all A/B/C/D/E matrices. This explicit
ODE path requires E to be exactly zero; differentiating loads need a separate
descriptor/mass owner and are never silently truncated.

The existing fs-orbit owner supplies fixed/adaptive scaled harmonic balance
and variational shooting with the complete discrete Floquet spectrum. The
original native/caller storage Hessian differentiates the same nonlinear
effort law, including material-memory states. Exact quadratic joint tangents
are cached; an unavailable or refusing declared tangent cannot silently become
a numerical tangent for variational shooting. Source contractions reuse
prepared scratch. Native J/R/G and the boundary reaction share one original
pHS compensated row in both the RHS and host-state tangent; a rounded native
field is never subsequently corrected by a separately rounded reaction.
Original flow and effort entries remain f64 owner outputs, and the independent
work audit reports the source/host exchange residual. No arbitrary relative
accuracy under severe coupled cancellation is asserted.
Explicit state/port/Fourier/payload caps and per-query
contraction limits remain separate from the shared solver's numerical budgets;
callback internals remain caller-owned. Cancellation and refusal publish no
partial orbit or accepted-source mutation.

This is the continuous ODE of the supplied rational model, including its
original real DC dynamics, rather than a realization inferred from harmonic
samples. Every retained positive harmonic must lie in the unchanged source
band. Physical qualification of DC, transient content and unresolved higher
harmonics still requires source evidence; source-band or broadband acoustic
validity is not inferred from mathematical solvability. The original continuous
and sampled passivity reports remain exposed with their original classification.
One shooting grid is numerical evidence, not a temporal convergence certificate.

`initial_state` explicitly imports the original accepted real pole history;
steady harmonic solves never insert that history as an extra forcing.
`source_candidate` delegates endpoint conversion to the original history owner
and returns a separate candidate retaining its model, reports and physical
clock. The enclosing trajectory owns acceptance. Independent cycle auditing
reevaluates original host H/R/G, authored supply, signed source exchange and
host-side physical reaction work. Arbitrary rational coordinates acquire no
invented stored energy, heat or radiation classification. Nonlinear work still
requires a quadrature sample ladder.

Five real-owner integration regressions cover full mutual real/pair harmonic
response and history, loaded and dark Floquet modes, a manufactured nonlinear
Duffing response, original DC/singular-D/band/budget/cancellation behavior, and
joint host/reaction cancellation in the actual RHS and complete tangent.
The actual production contraction module was executed with the real fs-math
dependency: four regressions passed for an independent passive Debye storage
identity, complex pole/residue response including DC, cancellation-sensitive
work, and dark-state/history evolution. The exact native pHS row plus its
original accumulator also passed 12 executable regressions, including the
three new joint-force cases. Native adapter integration execution
remains unverified because the unchanged full dependency build is blocked by
the upstream asupersync memory limit.

### Original-state sidebands of the loaded nonlinear orbit

`phs_orbit::impedance::sideband` binds the same prepared physical source and
its complete `[x; flow]` base to the existing fs-orbit periodic-linearization,
signed-sideband LU and analytic carrier-sensitivity owners. A finite admitted
table must cover every `abs(carrier+k*omega)` for every declared carrier before
any storage evaluation. Exact shared base frequencies reuse the base load;
additional supplied values there must agree. There is no interpolation, missing
load default, constant damping substitution or source-band extension. The full
native Hessian and all material-memory coordinates remain original.

Nonnegative physical samples use `exp(-iwt)`. Positive-time values at `w>=0`
are conjugated; at `w<0` the original sample at `abs(w)` is used. Authored
frequency slopes follow the same convention with the additional derivative of
`abs(w)`: at negative frequencies the original slope changes sign. Finite DC
values must be real and a two-sided DC slope purely imaginary. Full mutual and
nonreciprocal entries survive unchanged in port order. Samples establish only
their admitted pointwise load; they supply no causal realization or Floquet
stability evidence. The base orbit and base clock stay fixed in the derivative.

One checked complete tangent is reused across declared carriers, and one
factorization across physical drives and their exact derivatives. Original
physical output is independently refined through `G^T Hxx(x_base) dx`, using
the existing output observer and cumulative primal/derivative budgets. Returned
auxiliary flow is checked against that observation in original physical units.
Boundary effort derivatives retain both `Z*dy'` and `Z'*dy`. Returned physical
state/output spectra feed the existing acoustic observers; their numerical
residuals still assess all augmented equations. Source identity, base balance,
coverage, derivative declarations, payload/work limits and cancellation can
refuse without an unloaded fallback or a partial response.

Five real-owner integration regression definitions cover analytic dispersive
loads on both sides of DC, same-LU carrier derivatives, full mutual nonreciprocal
loads under conjugate port-unit changes, native nonlinear Duffing/Maxwell
sideband mixing, and coverage/base/cancellation refusal. Three tests of the
actual shared signed-frequency production leaf passed against the real fs-math
dependency. Full native adapter execution remains unverified while the unchanged
dependency graph is blocked by the upstream asupersync memory limit.

## Thermal porous solids with sealed source cavities (2026-11-06)

The original `solid_enthalpy::finite_plastic::porous::FiniteThermalPorousSolid`
now exposes `advance_with_cavities` and
`advance_with_cavities_and_coupling`. Every Picard or IQN-ILS proposal uses the
same accepted material, motion and thermal/observer history. The native thermal
porous solid solves original sealed gas pressure, current cavity geometry,
thermal constitutive response and optional consistent inertia together. The
final actual-temperature proposal must pass the existing cellwise heat,
temperature, thermal transport and complete energy acceptance gates before any
endpoint is returned.

`FiniteThermalPorousCavitiesCandidate` retains the original caloric timeline,
source-bound solid, independent observer and `FiniteThermalPorousCavityEquilibrium`
containing both native solid equilibrium and caller-ordered gas endpoints.
`FiniteThermalPorousCavitiesEnergy` retains material heat/storage, optional held
pressure, each original gas/exterior/work-correction receipt, compensated gas
storage/exterior-work totals and the complete signed external-work defect.

The combined identity is material internal-energy change + caloric storage
change + gas storage change + kinetic change + numerical kinetic loss + complete
external-work defect = dead/support work + held pressure work + explicit cavity
exterior work + external/incoming thermal energy. Internal gas work cancels
between solid and gas. Only native signed material caloric heat enters the solid
cells; gas storage, numerical pressure quadrature corrections and kinetic losses
are never deposited again as heat. Default zero cavity accounts preserve the
existing feedback arithmetic for uncoupled owners.

Cavities remain resolved homogeneous isentropic sources. Original undeformed
volume bindings survive all increments; output deformed gas parcels cannot
replace them. This interface does not model gas-wall heat transfer, pore-fluid
transport or gas waves. Solid conductances retain the supplied reference thermal
network. Synthetic source coefficients establish mathematical regression
expectations, not measured calibration.

Focused regressions in `tests/finite_porous_enthalpy/cavity.rs` cover an actual
free inner-wall displacement and independent hydrostatic/caloric/gas predictions
through both coupling methods; direct gas energy with a large arbitrary datum;
prescribed-motion heat parity with and without gas; retained consistent inertia
and compaction; original source identity; and rollback after observer, coupling
iteration and replacement-reference refusals. The pinned formatter parsed the
new source. Native compilation and execution were not available in this shared
workspace because disk and memory were exhausted; these tests are authored,
not reported as passed.

## Source-bound thermal frame malleability with finite enthalpy transport

`frame_enthalpy::ThermoplasticFrame` connects retained tabulated thermal
fibers to the existing quasistatic or backward-Euler frame solver and the
original deposition/enthalpy transport owner. Binding requires complete,
nonduplicated fiber coverage and checks each original cell's reference volume,
mass, material, caloric-law identity and accepted temperature. Original fiber
quadrature volumes and member inertia remain authoritative. Bound inventories
and accepted enthalpies are checked again before each advance.

Each feedback proposal prepares one source-bound thermal interval per cell,
shared by its fibers, using the accepted old temperature, proposed endpoint
temperature and actual elapsed time. Native frame mechanics retains plastic,
kinematic and rate-dependent history. All proposals start from the same
accepted history; changing a feedback temperature never advances time or
accumulates plastic flow. Original source-family and caloric bindings survive
subsequent intervals.

The observer adapter deposits the frame's physical material heat plus signed
reversible thermal exchange exactly once, as a mean source over the physical
interval. It forwards the original observer's sources, conductances, reservoirs,
ambient conditions, transitions, time-step limits, events and step hooks.
Native transport still owns thermal substeps, incoming deposition energy and
external heat. Material numerical defects and backward-Euler inertial loss
remain separately audited and never become heat.

After transport returns, the full frame is solved again at the actual candidate
enthalpies. Publication requires both temperature convergence and agreement
between the material heat injected into transport and that of the actual
temperature frame. Per-cell heat defects, the native thermal balance and the
combined material, reciprocal thermal, kinetic and caloric balance must all
meet their configured tolerances. Compensated sums retain opposing heat/work
terms. The accepted frame therefore has the exact transported temperature for
the next interval. Frame, thermal and observer candidates remain private until
all gates pass; failure or cancellation leaves accepted owners unchanged.

Four focused regression definitions cover loaded heating/cooling/release with
retained permanent extension and an independent first-law oracle; unequal fiber
temperatures changing extension and curvature; successive nonlinear dynamic
steps for both native mass forms with numerical inertia excluded from heat;
and wrong mass, failed feedback or observer failure with successful retry.
New production and test files pass Rust formatting checks, and independent
scalar arithmetic checks the fixtures. Cargo compilation and native regression
execution remain unverified on the resource-constrained shared build host.

This is a first-order coupled time advance within the native frame's linear
geometry, small-strain, fixed-reference-mass and stationary-support domains.
Thermal cell/fiber resolution and reference transport coefficients remain
explicit caller choices. The extension does not establish finite-rotation
forming, deformed heat-transfer geometry, experimental material validation or
a global time-accuracy certificate.

## Source-bound shell viscosity at a prepared temperature

`ShellMaterialViscosity::at_temperature` uses the existing calibrated material
shift and retains the original homogeneous rotated A/D section and its passive
factor. Every snapshot shares the same immutable reference source and section;
a later temperature selection always starts from those references and cannot
compound an earlier shift. `ShellViscosity::scaled` scales the original
coefficients and cached strain channels without re-factoring rounded entries,
preserving weak correlated physical work. The actual `TemperatureShiftState`
and transported frequency support accompany the prepared facet. Existing
constructors keep their unshifted defaults and original source getter semantics.

`ModalShellMaterialViscosity::project` checks actual retained eigenfrequencies
against the effective facet bands, including zero tensors, before the original
complete spatial/modal projection. Thickness, orientation, original mass,
elasticity and modes are unchanged; existing physical channels own the heat.
Tests verify rotated A/D and full modal-matrix/power scaling, repeatable
temperature changes, reference-source identity and both expanded and contracted
frequency admission. This is fixed-temperature preparation without live heat
feedback or changes to elastic stiffness.

## Resolved NSF convection against finite spatial solid storage

`nsf_thermal::advance_nsf_solid` advances the original `NsfGrid3` and original
`DepositionSimulation` together for one lattice step. Each unique native fluid
face maps to an explicitly bound solid cell and its same registered caloric
object. The face area is the original dx²; mapped patches must fit the cell's
exposed surface. The adapter removes their area from authored ambient exchange
and refuses a competing ambient-conductance override on those cells.

The actual solid endpoint temperature drives an isothermal interface or the
original molecular gas half-cell in series with `d_s/k_s(h) + R_contact`.
The explicitly bound original conductivity law retains its enthalpy/phase
admission; zero conductivity insulates. Solid conductivity is held through its
declared half-layer, without inventing an integrated profile or phase boundary.
The finite cell still owns all storage. Native reciprocal solid links and
independent sources/reservoirs remain in the same backward-Euler residual.

Each fluid trial retains its complete original fields, material and flux
history. The SAME local-face trapezoid joules consumed by the native fluid
energy equation are debited once from their respective solid cells. Opposing
face heat is never apportioned from a net wall total. Newton/Krylov probes
start from the accepted states and discard their private candidates; the final
thermal payload retains the exact final fluid trial for joint publication.
Original local/total thermal checks, independent fluid and per-cell/interface
transfer checks, and a combined storage/external-supply check must all pass.
Cancellation, source refusal, exhausted work and failed energy gates preserve
both owners. Cell/contact limits and a shared endpoint-query cap bound every
native fluid trial, including derivative probes and final acceptance.

This is fixed-geometry, fixed-inventory, first-order coupled time sampling:
the solid endpoint boundary is held over one original fluid step. Both clocks
use the same representable lattice timestamp. Native thermal events still
require explicit splits/processing. Molecular/local SGS, held open controls and
the original periodic/channel dynamic SGS closures are admitted. Dynamic trials
retain accepted transport and flux history while rebuilding only deterministic
filter scratch under the original controls, geometry and retained byte bound.
Channel fits consume each proposed wall temperature; modeled heat remains on
interior shared faces and molecular conduction owns the wall exchange exactly
once. Manufactured fluid source callbacks still refuse because their external
state cannot participate in the transaction. Original exploratory SGS limits,
source calorics, fluid Mach/temperature/stability bounds and chosen spatial/lattice
resolution remain authoritative. No moving solid, boundary-layer wall function,
global time-error bound or experimental convection qualification is implied.

Focused native regression sources cover an independent uniform-fluid and finite
wall capacity balance with series resistance, separate heating/cooling patches,
solid conduction and remaining ambient area, conductivity/insulation response,
source identity and work limits, and transactional late source/cancellation
refusal. These regression sources have not been executed on the current
resource-constrained build host; parsing is not native numerical validation.

## Local resonance surveys of complete nonlinear forced responses

`phs_orbit::forced_resonance::survey_forced_port_resonances` connects the original
forced harmonic-balance solution, complete analytic frequency tangent, nonlinear
physical port observer and shared driven-peak search. Every query predicts all
original state coefficients, including material memory, from the same admitted
anchor and its logarithmic frequency derivative, then corrects through the
existing physically scaled HB solve. Per-state physical RMS correction limits
and a declared logarithmic frequency interval bound this local operation.
Query order does not choose the seed, and correction limits do not prove branch
uniqueness. Singular frequency graphs and source/resolution refusals remain
failures; pseudo-arclength continuation and stability analysis retain their
separate existing owners.

The output is an authored port/harmonic of `G^T grad(H)`. Its derivative uses
`G^T Hxx dX/dln(omega)`, the original coherent input phase, and both original
output-resolution gates. Peak phasors use `2 conj(c_k)`; derivatives are per
DRIVE fundamental rad/s even when observing harmonic k. The shared survey owns
query/peak caps and bracket refinement. Each synchronous HB correction is
bounded by its native caps; cancellation is polled between owners and inside
the existing cancellable tangent/observation/search operations. No orbit cache,
scalar fitted resonance, artificial loss, completeness or stability claim is
introduced.

Three focused integration regressions cover independent mass/spring/damper plus
Maxwell response and complex phase at harmonics one and three, finite-amplitude
nonlinear restoring-effort peaks with changing waveforms, and source-clock,
anchor, branch-distance, resolution and cancellation refusals. These tests are
authored but not claimed executed without the complete `fs-couple` dependency
closure. Formatting and bounded native-API/physics review passed.

## File-selected distributed viscothermal sampling in flues

The existing flue performance source accepts an optional final sampling token
on its `viscothermal` header: `bilinear-midpoint` or `exact-held-flow`. This
works after `none`, a tube's `CELLS MIN_HZ MAX_HZ`, or `sections COUNT` followed
by the original addressed `loss_section` records. Omission selects the same
midpoint arithmetic as legacy files; unsupported or additional tokens refuse.
The selection covers both viscous series and thermal shunt memories of every
selected physical section.

One `ViscothermalTimeOptions` value carries that choice and the existing
fractional transit policy into the original network constructor. The shared
baffled-environment constructor accepts the same options, so source-bound
material walls, outlet flanges and reaction/observer composition cannot reset
the selected distributed memory. Their own constitutive and radiation laws
remain separately selected. Regional gases, physical Foster spectra, source
bands, fit/scattering gates and interval-mean work ports retain their original
owners. Exact selection adds no acoustic filter or new fitted coefficients.

## Five-strain constrained tensor memory

`PlaneStressPoint::new_shell` retains engineering strains
`[epsilon_xx, epsilon_yy, gamma_xy, gamma_yz, gamma_zx]` and enforces only
`sigma_zz = 0`. The original `new` constructor and descriptor retain their
three-traction plane-stress meaning. Both constructors use the same complete
equilibrium/Maxwell energy reduction; all six coordinates of every original
branch and all in-plane/shear cross terms survive. The shell coordinate budget
includes one local thickness-normal strain. `normal_stress_pa` reports actual
`[sigma_zz, sigma_yz, sigma_zx]`; its transverse components need not vanish.
`thickness_normal_stress_pa` observes the one shell constraint explicitly.

The existing full-inertia body owns motion, work and independently integrated
heat. A resolved thickness-normal Kelvin viscosity retains one massless normal
history and requires velocity-aware traction observation. Original passivity
and exact-propagation admission still apply; arbitrary unresolved singular
Kelvin resistance is not newly admitted. `relaxation_image` continues to
refuse every nonzero Kelvin term, preserving the pure-Maxwell meaning.

A unit-volume identity five-strain `PlaneStressRelaxation` can be mapped with
`project_shell` into physical columns and positive quadrature volume. Its
rates and original-history transformation are copied exactly, while strain
rows scale by the square root of volume and equilibrium work by volume. The
original `new_shell` point remains the source observation owner; its complete
energy must agree with the mapped reference energy. Preparation and comparison
share an explicit work budget. No new eigensolve, approximate pole merge,
branchwise constraint, damping fit or time integrator is introduced.

## Oriented anisotropic porous-body momentum exchange

`porous_momentum::free_body::FreePorousBody::new_with_support` binds centred
spheres or body-aligned ellipsoids to the original D3Q19 fluid and finite rigid
body. `step_candidate_with_principal_resistance` carries three nonnegative
Darcy coefficients in the body's material axes and one isotropic Forchheimer
coefficient. Native lattice force density uses physical fluid/skeleton slip;
principal values are not accelerations or unscaled intrinsic permeabilities.
The existing `free_spherical` module and `FreeSphericalPorousBody/Step` names
remain aliases of this common owner, and the scalar sphere entry retains the
native isotropic representation.

Every six-impulse trial derives one geodesic midpoint attitude from the actual
Cayley endpoints. That frame rotates both the cell-centre ellipsoid mask and
the resistance tensor. The physical skeleton velocity remains the original
work-conjugate COM velocity plus angular velocity crossed with the actual arm.
Full dry inertia, its gyroscopic motion and the geometry/resistance axes remain
independent. The original native flow solve supplies the reaction impulse and
couple; its existing checks gate both participants' actual linear/angular
momentum, midpoint work and nonnegative relative dissipation before publication.
Incoming physical fluid momentum survives native force-frame rebinding when
the material mask moves.

Support dimensions, dry mass and inertia stay fixed. The enclosing sphere must
fit within the closed stationary box at both endpoints, conservatively covering
the straight COM segment and intervening rotations. The actual drag support is
the oriented ellipsoid, not that bounding sphere. No excluded volume, interior
fluid subtraction, porosity, thermal sink, moving impermeable wall, deformation,
or subcell overlap model is inferred. Cell-centre mask discontinuities can
require a smaller tick; failure or cancellation retains both original owners.
Material axes are prepared once per trial, and the existing scratch grid and
bounded native trial/cell work remain shared; no second integrator or per-cell
matrix factorization is added.

Seven focused native tests of the production support/frame code with original
MBD vector/quaternion arithmetic pass. Three original-owner integration
regressions specify anisotropic force/torque/work from independently reconstructed
population velocities, scalar/sphere equivalence, and transactional refusal.
Those complete coupled regressions remain unexecuted on the exhausted build host;
no full-workspace or throughput claim is made.

### Intact viscoelastic bonds on the cohesive load line

`CohesiveLoadLine::trial_relaxation` consumes the original source-bound
`fs-material::cohesive::relaxation` history and a prepared exact Maxwell ramp.
It solves `g = g_free - S A mean_traction(g)` through the same full 3x3
mechanical compliance, area and LU workspace. The hereditary mean map is affine,
so one direct solve replaces nonlinear iterations; endpoint traction is never
substituted for the force whose work advances the bodies. Original jump-increment,
force, material-work, compatibility and cancellation gates remain active. Refusal
leaves accepted internal strains and elapsed time unchanged. The caller retains
its native free-motion map and complete energy audit. This is an intact massless
thin bond, with independent normal/shear relaxation; no fracture, friction,
finite-rotation objectivity or exact continuous inertial trajectory is inferred.

The native regression advances a bonded mass using the physical mean force,
checks independent relaxed displacement limits and closes external work against
kinetic energy, stored bond energy and independently integrated physical heat.
A separate ordered-interface regression supplies the original source spectra and
checks missing-mode, units, external-temperature and elapsed-time refusals.

## Attached adhesive microslip on the complete native modal network

`contact::adhesive_dmt::mindlin::AdhesiveDmtMindlinModalSystem` connects the
original sourced low-Tabor DMT potential and changing-load Hertz--Mindlin
history to actual native body motion. Two complete work columns retain all
modes, connections and the full reciprocal normal/shear response. Each bounded
normal probe reuses the allocation-free original shear recoil solve; a
caller-declared attached indentation bracket closes normal compatibility.
Only replay on the actual native endpoint creates the successor history.

Normal adhesion storage, reversible tangential storage, physical microslip
heat, native physical losses and numerical loss remain separate. Contact-sized
force, coordinate, transpose-work and storage/heat gates pass independently of
unrelated bulk storage. The full balance composes the original native increment
account with the contact account, preserving small contact work beside a large
preload. Receiver refusal, domain/solve refusal or cancellation publishes no
body state, contact history, clock or last frame. The actual native network and
history can be recovered without resetting either.

The admitted model has nonnegative indentation endpoints, fixed source
conditions and work columns, affine local paths and one proportional shear
direction. The caller resolves loading reversals and supplies the root bracket;
no uniqueness, contact detachment/reattachment, JKR microslip or automatically
routed heat is inferred. Explicit step, travel, speed, force and evaluation
budgets remain active. G1/G3/G4 regression definitions cover independently
integrated Mindlin loading at zero/tensile net normal force, the analytical full
two-mass/oblique-spring response including its passive third coordinate,
loading/unloading and shear reversals, complete work, and atomic refusal/retry.
Changed Rust parses with the pinned toolchain; these native integration
regressions remain unexecuted in the full-storage workspace environment.

## Source-resolved full tensor memory on the primary shear shell

`plate_viscosity::shell::relaxation::tensor::TensorShellRelaxation` binds each
original MITC facet to a homogeneous `ResolvedTensorMaterial`, explicit layer
faces, proper material rotation and positive shear correction. Its original
seven area sites and two thickness Gauss points retain membrane/bending work,
both transverse shears and arbitrary Maxwell in-plane/shear couplings. The
complete material constraint eliminates only sigma_zz, through the original
`PlaneStressPoint` owner; branches are never reduced independently.

Each complete equilibrium A/B/D/shear section and all three corner reference
inertia moments must agree with the original source. Unsupported equilibrium
cross terms refuse. Original support expansion, static bubble condensation,
mass choice and mass-orthonormal modes remain authoritative. The explicit shear
correction acts once as a congruence on both retained shears and mixed work;
this equivalent-single-layer model does not recover a three-dimensional shear
profile. Numerical drilling never becomes material memory or physical heat.

The original compensated projection Phi^T K Phi must be nonnegative and
diagonal in the supplied basis. Its retained diagonal also checks the target
component's frequencies before attachment; equal mode count is insufficient.
Spatial mode signs and orientation in degenerate eigenspaces remain caller
obligations because the generic modal runtime carries no spatial mode shapes.

The unit-volume material coordinate transform is prepared once per facet.
Exact linear projection retains identical source-constrained rates; the
existing equal-time coordinate owner shares driven histories without fitting,
pole merging or rank truncation. Every original point history reconstructs
before strain/traction observation. Setup, retained states and observation
points have explicit capacities. Relaxed and unrelaxed initialization use the
original modal-memory owner; no new propagator or loss account is introduced.

The coupled attachment checks its complete linear frequency bound against
all source bands and installs accepted-endpoint quadrature strain checks,
including a conservative allowance for unscaled transverse strain. These are
point checks, not spatial or between-step peak certificates. Kelvin viscosity
is explicitly refused by this pure-Maxwell image; temperature is frozen and
no live thermal work, measured damping authority or nonlinear failure law is
inferred. Independent full-tensor complex constitutive contraction and the
existing exact pHS/RK4 heat comparison are the focused numerical test seams.

## Pressure-driven finite plasticity through original caloric transport

`solid_enthalpy::finite_plastic::FiniteThermoplasticSolid::advance_with_pressure`
and `advance_with_pressure_and_coupling` retain the original held pressure field
through every native thermal proposal. The existing finite material source,
complete plastic/creep metric and accumulated history, physical interval,
caloric binding and reference inertia remain authoritative. Both relaxed Picard
and IQN-ILS use the shared enthalpy feedback owner; its actual-transport-
temperature re-evaluation and per-cell heat/combined-energy gates decide
acceptance.

The pressure-aware candidate retains native material energy and the original
pressure path work, endpoint work and geometric correction. The complete
external-work defect includes the material and pressure terms once. Physical
yield/viscous/creep heat and signed reciprocal exchange enter the original
caloric transport exactly once. Pressure quadrature correction and numerical
kinetic loss never become heat. The native absolute-sum work budget prevents
opposing integration errors from hiding behind a small signed total. Existing
unpressured result types retain their default energy account and arithmetic.

Two focused regression definitions use the original sourced thermal hardening
surface with explicit Arrhenius Perzyna kinetics. Independent finite coaxial,
pressure-volume, caloric and consistent-mass equations cover loading, heating
under held pressure, release and successive dynamic intervals; both feedback
methods, exact accepted temperatures, original source retention, late refusal
and successful retry are exercised. Independent fixture arithmetic was checked.
No Rust formatting, compilation or native tests ran after execution-host loss.

This retains the fully solid, compressible P1 finite-strain and first-order
time domains, original reference mass and reference transport coefficients.
Deformed thermal conductances, contact tooling, exterior-fluid evolution and
experimental forming calibration are not inferred.

## Inertial dielectric heat on the shared enthalpy timeline

The existing `step_dielectric_field_candidate` thermal bridge accepts the native
planar, axisymmetric and tetrahedral Lorentz candidates through
`DielectricFieldState`. Its associated material exposes the complete original
source snapshot and exact coefficient roles through `DielectricMaterialState`.
Debye retains epsilon/tau/increment roles; Lorentz retains epsilon/increment,
angular frequency and damping roles. Every selected coefficient must remain
validity-wide constant, and the same original claims are requeried at each
actual receiving temperature. Temperature-dependent pole changes are refused;
this path supplies no missing parameter-work term.

Only independently computed physical pole friction enters the existing local
enthalpy residual. Instantaneous, restoring and kinetic polarization stores
remain in the field, and the combined energy window closes against original
electrode work. One immutable field candidate serves all thermal Newton trials;
no trial repeats its physical history or resets its polarization current.
Accepted field and caloric states publish together. Existing source identity,
volume correspondence, fully solid caloric charts, failure rollback, thermal
budgets and optional Ohmic accounting remain active.

Two added integration regressions specify damped/lossless original current
histories, source-bound thermal snapshots, combined work/storage/heat closure,
and refusal of varying damping or cancellation without partial publication.
They remain unexecuted without the full native `fs-couple` dependency closure.

### Original rational-load sidebands and signed carrier derivatives

`phs_orbit::rational::sideband` binds an existing `PreparedPhsMimo` and a complete original joint `HbOrbit`. `PhsMimoSidebandSource::prepare` delegates independent base-balance and exact periodic-tangent admission to `fs_orbit::sideband::PreparedPeriodicLinearization`; stored base residual fields are not accepted as evidence. The original host coordinates and every supplied continuous A/B/C/D source coordinate remain in that tangent, including zero-residue states and real/negative-imaginary conjugate-pair coordinates. This is the exact small-signal evolution of the already admitted proper rational model. The enclosing causal owner still refuses nonzero E; the adapter does not drop differentiator physics, refit, invert D or B, sample a replacement impedance, or infer a source storage/heat/passivity witness.

One prepared tangent can factor multiple signed probe carriers. Each `FactoredPhsMimoSidebands` factor can solve multiple original drive-port probes. External forcing maps through the original host G and appends exact zero source-coordinate forcing; source excitation occurs only through its actual B and the physical S attachment. `solve_with_carrier_derivative` reuses that same joint LU with exact `dP(i nu)/d nu = i I` on both host and source coordinates. Optional authored external forcing slopes are mapped through the same G. The base orbit, base clock and original source coefficients stay fixed; this operation does not differentiate the nonlinear operating point. Original and differentiated complete-equation backward errors and finite omitted-equation residuals are retained separately.

Returned `PhsMimoSidebandResponse` keeps the existing physical host response in its honest n-state, m-drive-port chart and carries q original source-coordinate coefficients in a separate sideband-major vector. Together these vectors are the complete joint state response; the host chart's retained `SidebandResponse` diagnostics continue to cover all joint equations. Boundary port count p and S may differ from m and G. The shared original-Hessian observation/refinement owner evaluates `G^T Hxx dx`, `S^T Hxx dx` and `D S^T Hxx dx` during the same real/imaginary Hessian pass. Exact retained `C dz` coefficients are then added to the sampled direct effort. Reported output order must cover state order, so this retained C contribution introduces no additional output tail. Separate positive drive-flow, boundary-flow and boundary-effort tolerances admit retained sampling change plus measured omitted-envelope norm; flow convergence alone never establishes effort convergence. Response and derivative return separate drive, boundary-flow and boundary-effort traces. Fixed-order finite observation/omitted-equation checks are not infinite-tail bounds or stability certificates. Original flow, direct-effort and Fourier interfaces remain admitted f64 outputs; this does not promise arbitrary relative accuracy at severe coupled cancellation.

All input, state, source and output coefficients use peak complex envelopes
`Re(exp(i nu t) sum_k X_k exp(i k omega t))`, with explicit signed sidebands, no conjugate folding and no factor-of-two conversion. Carrier derivatives follow this signed convention through zero. Every nonzero retained, reported and finite omitted-equation guard-band frequency must pass the unchanged original positive source-band admission. Exact zero uses the supplied model's real DC law, is disclosed by `probe_includes_model_dc`, and does not extend physical source qualification to DC or transients. The returned positive `source_band_rad_s` remains the original admitted band. Source reports and committed history are unchanged.

A source holds an exclusive constructor borrow; prepared tangents and factors retain that discipline. Native tangent, LU and derivative caps remain separate from the adapter cap. The adapter reserves its host-base copy, owned tolerances, original joint/source response capacities, input mappings, mapped output spectra, and split traces before queries. Shared physical actions/work apply cumulatively across the primal and derivative observations. The additional D contraction is charged on every observation sample; original source contraction caps apply during actual base/tangent evaluation. Cancellation, malformed charts, source-band gaps, insufficient full-state dimensions, missing exact storage Hessians, unresolved physical envelopes and budget refusals publish no partial response and never mutate source history.

`tests/phs_rational_sideband.rs` supplies focused real-owner regressions for singular D and retained dark source states; signed carriers, DC and authored forcing slopes; full mutual real-plus-pair realization with p different from m; manufactured nonlinear frequency conversion and exact source/derivative equations with a single shared Hessian pass; complete-base revalidation, original-band/guard-band limits, cumulative pair budgets and cancellation; and static zero-state realization plus unavailable-Hessian refusal.

Source/API review and exact GitHub blob readback completed for this increment. Native compilation, Rust formatting and test execution remain unverified while the shared workspace is unavailable.

## Heated cyclic frames on the shared enthalpy owner

`frame_enthalpy::FrameThermalSource` admits the original tabulated family and
one through eight native Chaboche branches through one source-bound cell
interface. A single frame may mix these families while retaining every original
fiber history, quadrature volume, member mass and caloric object. Each family
prepares one immutable physical interval per cell; fiber clones share it during
the existing equilibrium and backward-Euler iterations.

The complete thermal inventory admits family, branch ordering and accepted
temperature before any coupled proposal. Original nonlinear branch recovery
and viscous overstress are physical heat channels already included once in the
native material total. No branch aggregation, source substitution, recovery
reset or numerical-inertia heating is introduced. Existing per-cell heat,
actual-transport-temperature re-evaluation and combined first-law gates publish
the complete frame, thermal timeline and observer together.

Two G0/G1 regression definitions cover heating, reversal, cooling and unloading
with independent axial branch recurrences, recovery dissipation and retained
shape, and a mixed tabulated/rate-dependent cyclic member whose two original
thermal cells receive different heat. Original-source mismatch and late
observer failure leave accepted history unchanged and permit retry. Independent
scalar force/heat calculations checked the cyclic fixture. Rust formatting,
compilation and native test execution were unavailable after the execution
host failed; these are authored numerical regressions, not executed or
experimental validation.

The source family changes no geometry domain: native frames remain small-strain
and linear-geometry with fixed reference mass and stationary supports. Finite
plastic pressure forming uses the separate existing finite-solid owner.

### Finite solution reservoirs for original heterogeneous layers

`FiniteMoistureSolution::advance_layered`, `advance_shared_layered`,
`advance_layered_radial` and `advance_shared_layered_radial` now couple the
original slab or cylindrical/spherical layer stacks to one finite solution.
Each bounded solvent-volume trial advances the complete native layered body
from its accepted state, including each source-backed sorption relation, actual
dry mass, geometry, contact resistance and effective seal. No homogenized
material or replacement field is introduced. A private participant interface
shares the existing solution root and phase owner with homogeneous specimens;
shared IDs, clocks, temperatures and cell limits retain their admission gates.

Only the original exterior faces exchange with the finite bath. Their net
water exchange closes against actual solvent mass, while internal interface
transfers cancel between neighboring materials. Every candidate solution
retains the original nonvolatile solute totals and native equilibrium or kinetic
phase partition. Sealed bodies may redistribute water internally while the
original kinetic phases evolve. Exposed unsupported pressure/composition source
axes and zero-area radial-center exposure are refused.

This is the existing isothermal vapor-equilibrium water-activity bridge, not a
liquid-immersion, osmotic-traction or heat model. Native body steps and endpoint
volume/kinetic phase coupling remain first order. Five authored regressions
cover independent finite-bath balances, both curved geometries, opposing shared
exchange, sealed internal diffusion with kinetic dissolution, and source/budget
refusals. Rust parsing and formatting were checked; these new integration
regressions have not been executed in the full dependency build.

### Adaptive heterogeneous slab and radial diffusion

`MoistureLayeredSlab::advance_sdirk2_adaptive_to` and the backward-Euler and
SDIRK2 `MoistureLayeredRadial` duration methods reuse the original bounded
step-doubling controller. Every original layer contributes its local moisture
and source activity, and an additional L1 physical-water discrepancy weights
each local difference by that cell's actual dry mass. Unequal shell storage and
opposing layer exchanges cannot disappear through signed cancellation. Its
absolute scale is `reservoir_water_absolute_tolerance_kg`; relative scales use
transported water from the last accepted state.

The second-order controller charges both native implicit stages, including
failed trials. Activity-stage range refusals shorten the proposal without
clipping the state or changing the method. Only accepted fine transfers enter
the complete-duration ledger: each internal interface has one debit and one
credit, and the whole-body account contains only exterior exchange.
`FiniteMoistureAtmosphere::advance_layered_radial_adaptive_to` uses the original
first-order coupled endpoint and additionally compares actual vapor mass and
humidity; final checks retain both total and per-carrier gas mass.

Geometry, temperature and boundary conditions stay fixed within a duration.
Local error estimates do not certify global accuracy, and cancellation or a
late failure publishes no accepted prefix. Five authored regressions cover
independent continuous exchange, both time methods and curved geometries,
mass-weighted error with sealed exteriors, an inadmissible SDIRK stage followed
by a bounded retry, work exhaustion after an accepted interval, and finite-air
layer/gas accounts. Source and fixture review and Rust formatting passed; the
new native integration tests remain unrun in the full dependency build.


### Heat exchange through retained finite porous gas cavities

`solid_enthalpy::finite_plastic::porous::cavity::heat::ThermalPorousGasSystem`
retains original ideal or Peng--Robinson gas inventories, actual prepared closed cavities,
finite thermal porous solid history and explicit fixed gas-wall conductances.
This adds heat transfer to the native finite-geometry compression path. It
does not replace the source calorics, material return map, pressure tangent,
thermal cell network, optional inertia or retained plastic volume.

Construction admits the exact original mesh, gas volume at accepted deformed
geometry, independent nonoverlapping cavity faces, original structural cell
identities and unique heat-link IDs. It bounds the system to 64 cavities and
16,384 explicit links. Gas mass, source model and auxiliary conditions remain
with the original parcel; cell mass and source caloric identity remain with the
native thermal/solid binding. Conductance is supplied as a total W/K per link.
There is no inferred moving area, source conductivity, convection correlation,
pore-fluid transport or new gas equation of state.

Each interval uses an explicitly first-order heat-before-motion split. A
bounded Newton--Krylov solve determines one rigid heat increment per cavity.
The original gas energy inversion proposes the heated state at the actual
accepted volume. A native `CavityGasEpoch` binds that physical state to the
already-admitted geometry; no fictitious expansion to undeformed volume or
extrapolated thermodynamic state is constructed. Accepted compact anchors are
prepared once and shared through `with_gas_at_anchor` across caloric probes.
Every mechanical proposal still admits its actual trial geometry natively.

The heated gas temperature is held during the original explicit wall
substeps, while the existing finite thermal porous feedback closes caloric
temperature, permanent compaction, pressure and optional inertia. Every outer
probe starts from the same accepted solid, gas, thermal clock and independent
observer clone. Native heat transport and all original observer source,
conductance, reservoir, ambient, internal-transfer, interval, transition,
event and timestep hooks remain active. Multiple links can share a cell; their
combined nonnegative decay conductance enters the original mass/caloric-slope
timestep restriction. Native event and thermal work budgets remain authoritative.

Per-gas balances require rigid gas energy change plus integrated link heat
into the solid to equal prescribed external gas heat. Gas rigid heating and
native adiabatic energy increments are reported separately; their compensated
sum is the actual gas energy change. Internal film heat is removed from native
external solid heat exactly once. The complete first law is the native
solid/gas/mechanical residual plus rigid gas heat plus represented film heat,
minus prescribed external gas heat. Native plastic heat, gas pressure work,
exterior work, support/dead loads, endpoint defects and numerical kinetic
losses retain their original owners and are not reclassified as heat.

The wall observer retains per-link original heat integrals and a separate
integral of the represented cell source. It rejects vanished contributions
and bounds the sum of absolute cell source-representation errors over actual
accepted durations. The shared finite-plastic heat composer now also checks
both observer-plus-material and authored-boundary-plus-composed additions.
This closes a real conservation hole: rounding away a weak interface source
beside a large material source could otherwise be misclassified as an external
sink while aggregate balances passed. Compensated absolute errors cannot
cancel between cells or substeps, and the existing feedback energy gate
applies before every complete candidate can be published.

The returned `CavityHeatCandidate` contains the new thermal state, retained
gas/solid system, caller's original observer, complete native reports,
per-gas and per-link heat accounts, nonlinear report and actual trial count.
Gas/source refusal, late observer failure, insufficient trial work, native
feedback failure or cancellation publishes none of these endpoints. The final
physical evaluation consumes the same global trial budget as Newton and
Jacobian probes; every gas, interface and combined energy gate must pass.

Three focused integration regressions in
`tests/finite_porous_enthalpy/cavity/heat.rs` specify an independent finite-gas
and wall-capacity balance, shared-cell links on the actual native substep grid,
heated variable-Cv compression over retained plastic/inertial intervals,
source/work/late-observer rollback, and heat-induced motion of an actually
free wall coordinate under native cavity and exterior pressure. Two focused
shared-feedback arithmetic regressions retain IEEE rounding counterexamples
for swallowed heat and noncancelling absolute discrepancies. Existing native
epoch tests also exercise explicit constant-cost source replacement.

The execution host was unavailable during implementation. Source/API and
physical-account reviews, plus an independent scalar heat-balance calculation,
were completed; no new formatting, compilation or native test pass is claimed.
Energy closure is not a time-accuracy or material-calibration certificate.
Refine the complete physical interval for simultaneous heat/motion accuracy.
Gas remains well mixed and frozen in composition in this solid-wall owner;
reacting thermal compression has its separate original gas reactor owner.

### Adaptive finite-solution exchange across unlike materials

`FiniteMoistureSolution::advance_layered_adaptive_to` and
`advance_layered_radial_adaptive_to` now complete a finite-bath duration using
the original first-order step-doubling controller. Each full/half probe solves
the native solvent-volume root and every original layer, including contact
resistance, source sorption and unequal curved storage. Temporal comparisons
include actual solvent mass and activity, each local material activity/content,
and the L1 dry-mass-weighted local water discrepancy. The last comparison
remains effective for sealed internal redistribution with constant pure-solvent
activity and zero net uptake.

Accepted fine-step receipts alone form the per-material interface accounts and
the exterior whole-body account. The original bound-water-plus-solvent inventory
is checked over the complete duration against one fixed tolerance. Work
exhaustion after a successful private interval, source refusal or cancellation
returns no partially advanced bath or body. Kinetic phases retain their native
first-order endpoint-volume approximation; phase-component temporal error is
not separately estimated. Fixed geometry/temperature and the equal-pressure
vapor-equilibrium bridge retain their existing scope.

Four authored regression definitions compare finite-bath time refinement with
independent continuous references for slabs, spheres and cylinders, isolate
physical-water-driven refinement in a sealed body, and exercise late work
rollback. Rust formatting and analytical fixture review passed. These new
native integration tests remain unrun in the full dependency build.

### Original-coordinate tensor shell radiation

`tensor_shell_radiation::TensorShellRadiation` binds the immutable original `TensorPronyShell`, its complete supported physical coordinate map, and an explicitly supplied homogeneous fluid. The original oriented shell midsurface is the closed exterior acoustic boundary. Uniform planar quartering preserves the supplied geometry and P1 translational trace; it does not project to a guessed shape, add thickness offsets, cap an opening, flip winding, infer a cavity, or create an open-screen pressure-jump model. Every original free translation and physical rotation remains in the original material/inertia solve, including dark directions and zero direct rotational normal-velocity rows. Clamped faces remain as acoustic scatterers.

Each child panel receives the area mean of the original P1 translation dotted with its exact outward normal. The existing reciprocal P1/P0 velocity-port owner solves all original unit-coordinate columns with one boundary factor and forms the full self/mutual impedance using the transpose velocity map and actual panel areas. Identity generalized ports retain the original displacement/rotation and force/torque units; panel areas are not applied a second time. The existing `TensorPronyShell::prepare_harmonic_response_with_port_impedance` owns the complete dispersive material reduction and loaded factor. Repeated forces reuse both preparations; no new eigensolver, modal basis, radiation resistance, damping repair, material-loss fit, or pressure-density substitution is introduced.

One coherent field uses the solved original mechanical velocity to combine the same P1/P0 source densities. Original adaptive receiver rows evaluate coherent pressure and particle velocity, including the existing 1..=64 receiver and directional-microphone routes. Receiver geometry can be admitted before frequency source solves, and its aggregate quadrature budget remains with the native observer. The loaded preparation exposes the native material preparation for source-layer stress/heat observation of the same displacement.

Publication requires three explicit numerical work checks: projected physical-port work versus coherent surface work, coherent surface work versus independently integrated outgoing flux, and original input work versus intrinsic tensor-Prony/Kelvin heat plus outgoing flux. These quantities remain distinct; no absorbed-work field is silently relabelled as material heat. Signed roundoff is retained, normalization uses the supplied physical powers without a one-watt floor, and a failed check refuses the result. Native full-cycle strain and physical-row residual admission remain authoritative.

Explicit setup, original-coordinate/source-count, per-frequency boundary, native mechanical/port, coherent-field and outgoing-flux work limits admit complete requests without dropping directions. Native closed topology/orientation, wavelength, angular and slab sub-resonance restrictions remain unchanged. Acoustic planar refinement, shell spatial/material refinement and angular integration are independent accuracy axes. This forced-frequency bridge does not establish a causal time realization, loaded resonances/decay, continuum convergence, a steel spherical-shell benchmark, open two-sided radiation or closure of bead .106.

Three native regression definitions cover independent original P1 interpolation and a separate combined-velocity BEM solve; original loaded motion, material heat and coherent outgoing work; reusable complex-force scaling and pressure/velocity/directional arrays; dark rotations and clamped scattering faces; and original material/band/geometry/resource/cancellation refusals. Native Rust tests and formatting remain unexecuted because the workspace is unavailable. Exact GitHub blob readback verifies stored source integrity only.

## Material inertia on resolved axial pressure waves

With `experimental-shock-loading`, `cavity_phs::piston::material::resolved::MaterialResolvedPiston` retains the original `ResolvedBody` and a private native `fs-gas` `Piston`. Mass comes exclusively from the body's canonical spatial inertia; immutable views preserve original region geometry, density queries and source receipts without a competing mass parameter. The gas mesh remains the sole owner of wall area and position.

Frozen and evolving-mixture entries delegate the original free or suspended piston solves unchanged. Mixture contact-equilibrium variants retain caller-owned scratch and the native signed/absolute numerical-defect account. Source calorics, species references, midpoint pressure impulse/work, spring storage, damping loss, acceptance budgets and final publication retain their existing owners. Duration arguments are upper bounds; native reports supply the accepted interval. Cancellation and refusal preserve all accepted physical state.

The retained body is a guided rigid axial specimen: gas heating does not re-query density or silently heat the wall, and rotational inertia is not substituted for total translational mass. This binding does not supply shell deformation, changing area or multidimensional FSI. Three focused regression definitions cover source-density response and native frozen/suspended identity, all four mixture report/field paths, and warmed-state final cancellation/refusal/retry. They are authored but have not been compiled or executed in the unavailable session executor.


#### Native real-gas heat and retained compression


The retained gas-wall heat owner is parameterized by the original gas parcel:
`ThermalPorousGasSystem<'a, G = ClosedGasParcel>` and its cavity/candidate/report
types preserve the ideal default and also admit `ClosedRealGasParcel`.
The tiny sealed `ThermalCavityGas` seam adds only original temperature, conserved
mass and native fixed-volume `add_heat` to `AdiabaticCavityGas`. Both owners
return the same original `GasHeatStep` and `GasEnergyError`; zero-heat behavior,
caloric energy datums, cancellation and atomic refusal remain native. Dense
Peng--Robinson parcels retain their actual density-dependent caloric departures,
homogeneous-phase admission and acoustic `rho*c^2`; no ideal or constant-Cv
approximation enters the heat solve. The shared geometry anchors, nonlinear
heat closure, wall observer and complete energy accounts are unchanged.


Two focused coupled regressions in
`tests/finite_porous_enthalpy/cavity/heat/real_gas.rs` use independent pure-fluid
PR Helmholtz derivatives. A finite gas/wall heat balance checks density-dependent
caloric storage, the true acoustic bulk modulus and invariance under a 1e24 J/kg
formation-energy shift. A retained two-interval heated compression checks each
new physical entropy epoch, density-dependent work, original inventory and
source identities, geometric volume agreement and source-failure rollback.
The fixtures are explicitly synthetic, not measured methane calibration.

An independent scalar calculation gives the rigid fixture gas temperature
300.742136567 K and wall temperature 300.253455915 K; signed gas and wall energy
changes are -0.529154483 J and +0.549154483 J for +0.02 J external heat. Its
acoustic bulk modulus differs from gamma times pressure by about 8.6 percent,
so an ideal-gas stiffness substitution is detectable. These reference
calculations and source reviews do not replace native execution: the host
remained unavailable, including a final bounded compiler-availability check.
No new formatting, compilation or native regression pass is claimed.

## Source-bound two-chamber pressure and finite-wall heat

The same `MaterialResolvedPiston` can drive both resolved mixture chambers through the native shared-piston owner, with optional stationary suspension and one explicit finite thermal inventory. Both pressure fields use the original material-derived mass. Native gas geometry, species, caloric storage, spring energy, physical damping and accepted pressure work remain independently accounted. The suspended heat path deposits only the authored fraction of actual damper loss into the original wall store; the remainder stays an external loss. Gas, wall inventory and piston publish together.

The thermal store does not implicitly re-resolve the rigid specimen's frozen density or infer expansion; representing the same specimen requires a compatible authored caloric inventory. The original negligible-thickness axial partition and split heat/pressure/heat assumptions apply. Three focused regression definitions compare all four entries with native reports and states, independently reconstruct the combined physical energy balance, and exercise warmed final-cancellation/retry across both chambers and the finite wall. They remain uncompiled and unexecuted in the unavailable session executor.

## Material propagation on the complete modal mouth field

The circular modal mouth observer applies the original material gas law to each signed surface path. It shares the existing five-rung spatial quadrature and original gas filter/convolution owner; all rungs consume one material-work allowance, and convergence comparisons remove only differences in reported numerical padding. The final receiver keeps its actual source report and latency.

Native brass `modal_baffle` and `observer_gas` compose using the authored exterior gas while retaining bore/lip dynamics and all mouth-mode flows. Scope remains a fixed circular aperture in an infinite rigid baffle. Focused regression definitions cover signed-mode pressure, original-owner parity, padding alignment, aggregate work limits, the actual three-mode brass render and callback partition equality. Native compilation and execution remain pending because the session executor is unavailable.

### Per-component accuracy for finite-rate solution phases

The backward-Euler finite-solution duration APIs now accept
`SolutionMoistureTimeControls`. Existing `MoistureTimeControls` arguments
retain their original behavior; `with_kinetic_phase_accuracy` explicitly adds
a positive molar allowance and a change-relative allowance for each actual
dissolved and crystalline stock. This applies to slab, Cartesian, radial and
layered bodies, and the existing shared slab/Cartesian reservoir paths.

The existing full-versus-two-half controller combines this component norm with
its water/content/activity norms and applies the representable-half correction
once. Both phase stocks are read directly from the original finite-rate owner;
a small crystal tail is not reconstructed by subtracting two large amounts.
Component sources/layout remain fixed, and bounded cancellation checks cover
the phase comparison. Aggregate water activity or total component conservation
cannot suppress a separately requested phase discrepancy.

This changes temporal error control, not the native endpoint-volume kinetic
law or its first-order accuracy. Only accepted fine histories contribute to
the final candidate; the original complete-duration water and work gates remain.
Explicit phase controls require an actual kinetic reservoir and are checked
before proposals. Algebraic-phase SDIRK retains its existing kinetic refusal.

Four authored regressions cover the coupled continuous water/phase trajectory,
equivalent slab/Cartesian and shared paths, exact fixed-volume phase evolution,
invalid controls and late work refusal. An independent RK4 reference checked
the physical fixture and its mesh refinement; Rust parsing and formatting
passed. The new native coupling tests remain unexecuted because the complete
dependency build is unavailable.

### Sorption heat across unlike original material layers

`FiniteMoistureAtmosphere::advance_layered_with_sorption_energy` and
`advance_layered_radial_with_sorption_energy` now pair the original layered
backward-Euler endpoints with one closed isothermal gas/body energy account.
Each material supplies its own pinned differential isosteric heat law on the
same original card and current transport context. The common cell-integration
helper uses actual dry masses and loading-dependent source integrals for both
homogeneous and heterogeneous bodies.

Per-layer reports retain their source, water change and actual cell bound/binding
energy increments. The original vapor caloric datum is shared and the native
gas-energy change enters once. A transfer between unlike materials can therefore
require nonzero thermostat duty even with sealed exteriors and no net gas water
exchange. Internal interfaces create neither a new gas inventory nor a second
heat debit.

One aggregate source-segment cap and one final energy tolerance apply to the
complete body. Unsupported heat-source endpoints, a late layer refusal, work
exhaustion or cancellation publish neither transport nor energy. Temperature
remains prescribed: this is the physical thermostat duty, not a spatial heat
equation, a uniquely located per-cell heat source or thermal feedback.

Five authored tests cover the exact sealed two-material binding-energy
difference, finite-air uptake/desorption with loading-dependent heat, unequal
spherical/cylindrical shells, partition equilibrium with zero integration work,
and late source/budget/cancellation refusals. Rust formatting and independent
analytical fixture calculations passed; native integration tests remain unrun.

### Conducting resonant dielectric heat in the native receiver

The existing `step_dielectric_field_candidate` receiver now accepts Lorentz
fields constructed with independent source-bound DC conductivity through the
same `DielectricFieldState` interface. No new thermal integration or deposition
path is introduced. The original field candidate supplies separate pole-friction
and Ohmic Joule losses, their sum is deposited once into the original physical
cell volumes, and restoring/kinetic polarization storage remains in the field.
The total balance includes both displacement and conduction electrode work.

Every original epsilon, pole increment, angular frequency, damping and selected
conductivity role survives the complete source snapshot. Existing validity-wide
coefficient-constancy checks apply to conductivity even when its selected value
is zero. Thermal trial snapshots refresh their actual temperature query while
the retained field owner and its fixed coefficients remain unchanged. Variable
conductivity, source mismatch and cancellation refuse the joint candidate
without consuming field history or heat.

Two focused integration regressions are authored for all four zero/finite
conductivity and pole-damping combinations over retained-history steps, original
source-role/temperature recovery, independent physical loss accounting, and
variable-conductivity/cancellation refusal followed by replay. They are not
claimed executed without the full fs-couple dependency closure. Carrier
transport, temperature-dependent resonant/conductive coefficient work and
induction remain outside this fixed-source bridge.

## Combined Kelvin-Maxwell constrained shell image

`PlaneStressPoint::viscoelastic_image` preserves the complete source Kelvin
matrix together with every Maxwell history. The existing `relaxation_image`
keeps its pure-Maxwell refusal and zero-Kelvin arithmetic. At most one retained
viscous normal coordinate is admitted by the combined image; the full body and
descriptor remain the owners for multiple coupled viscous normal histories.

With `L = D_nn^-1 D_nu`, the coordinate `y = sqrt(V) (n + L q)` removes the
normal/velocity skew coupling. The original positive normal mobility and
Maxwell rates enter the same symmetric relaxation transform. Its channels
include thickness-normal Kelvin relaxation, even for a source with no Maxwell
branches. `feedthrough_viscosity` retains the complete remaining matrix
`D_uu - D_un D_nn^-1 D_nu`, evaluated by the original factor-congruence owner.
Its velocity power and the relaxation-channel power must both be included
exactly once; no diagonal decay substitution or extra loss fit is implied.

`original_memory_at(q, z)` reconstructs original Maxwell coordinates and the
physical normal strain, including the `-L q` contribution. The displacement-free
observation refuses when that contribution is required. Stress observations
accept actual velocity and retain physical transverse shear. Unit-volume
`project_shell` continues to reuse identical rates; Maxwell history-map rows
remain unchanged, while physical normal-strain rows scale by `1/sqrt(V)` and
their displacement offset follows the same strain-column map. Complete physical
energy reconstruction, passivity, source domain and work-capacity gates remain
binding. No new constitutive law, eigensolver or time integrator is introduced.


## Native sourced intact adhesive relaxation (G2)

`render::schedule::force::coupled::contact::cohesive::relaxation::ViscoelasticBondModalSystem`
retains the original ordered interface thickness, independent normal/shear
Prony spectra, and original native mechanical owner. The explicitly bonded
sample-zero configuration has zero separation and branch strain; initial body
velocities remain physical input. Complete work columns and their transpose
include every body mode and the full three-by-three native response, including
cross-axis coupling.

Cached exact ramp coefficients and the existing LU workspace solve the affine
mean-traction load line once per step. An independent constitutive replay uses
the actual propagated coordinates. Two material evaluations, original material
jump/rate/time bounds, native limits, force compatibility, transpose work and a
bond-sized increment energy gate must pass before body and branch history
publish together. Receiver refusal or late cancellation retains accepted motion,
elapsed time and the last frame. Relaxation heat, other physical losses and
native numerical dissipation are separate accounts.

This is a small-strain, massless, intact layer in a fixed local frame with a
straight constitutive jump ramp. It does not infer damage, thermal strain,
rebonding or a heat destination. Refine the physical timestep for trajectory
convergence. G2 regressions in `tests/cohesive_modal/relaxation.rs` exercise
original two-body motion, independent normal/shear memory, complete cross-axis
response, energy/momentum and rejection/replay.
