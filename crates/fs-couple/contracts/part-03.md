<!-- fs-couple contract, part 3 of 8; a higher number is newer. Index: ../CONTRACT.md. -->
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

### Original viscous and molecular heat in the complete finite-gas load

`DgPortLoad::step_with_finite_gas_viscous_molecular_heat_controlled` explicitly
returns physical heat belonging to the main finite gas inventory. One original
native trial solves the inlet, outlet and finite-branch junction flows. The
main bore's source-positive viscous and molecular observers consume those
EXACT held inputs and the complete original state. Independent budgets cover
their sequential preparations; no state is truncated to fit the 128-coordinate
source-integral limit. Absent bore sources require no preparation budget.

Original resolved slit/circular terminal passages belonging to the main gas
contribute their positive viscous loss. Branch passages retain their independent
gas ownership. Backward-Euler numerical loss, solid and thermal-wall loss,
independent-branch loss and signed exterior exchange are not assigned to this
gas. Channel pressure work also pays kinetic storage and numerical loss and
cannot be substituted for the viscous term.

The existing caloric/parameter-work solve deposits the internal heat once and
deducts it once from exported loss. `FiniteGasDgPortStep` keeps raw native loss,
main molecular heat, main viscous heat, terminal viscous heat and exported loss
separate. The original source, volume, material-memory and cancellation gates
remain atomic. A prepared complete trial is reusable without advancing an
accepted state or repeating the acoustic solve. Active bore heat on the
continuous-port backend still requires its own trajectory integral and refuses.

Focused regressions independently reconstruct slit and circular shear heat,
NASA caloric energy, actual held-port source integrals, branch ownership and
late-cancellation rollback. They are authored and parsed; native compilation
and execution remain unverified in this environment.

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
partition does not infer mean-pressure prestress or thermodynamic mass flow
from acoustic flow. Default callbacks retain explicit exported dissipation.
Focused synthetic regressions cover scalar caloric/wave-work feedback,
independent paired-blade motion, exact retry after source refusal or final
cancellation, and rejection of an unowned second inventory. Native execution
and measured physical validation remain separate gates.

### Physical gas heat from the original moving-aperture load

`ApertureDg::step_with_finite_gas_viscous_molecular_heat_controlled` returns
main-bore molecular/viscous and main-inventory terminal-passage shear heat
through the same caloric/material-work solve. Scalar, independent-pair and
material-modal mechanics retain their original nonlinear/ExactZoh solve.
One private marched moving owner supplies the actual held DG inlet and
implicit outlet/branch flows. Later thermal and material-work trials reuse
that proposal without repeating mechanics or acoustic source integration.

Raw coupled-frame losses remain intact; the gas enclosure removes the internal
transfer once. Mechanical damping, contact, gap/jet losses, independent-branch
loss, wall/solid and numerical losses retain their exported ownership. No
moving-gap, supply-passage or second-tract gas inventory is inferred. Original
complete-state heat budgets, source refresh, the fixed-volume partition and
atomic refusal remain in force.

Four focused regression groups cover independent slit shear/caloric balance,
actual moving-port molecular/viscous integrals with unequal blades, zero-main-
heat/independent-branch identity, and budget/late-cancellation rollback. Source
parsing and scoped whitespace checks completed before executor loss; recovered
manual source patches were reviewed and verified through GitHub. Native
coupled compilation and execution remain unverified.

## Finite regional gas transport on shared acoustic motion

`DoubleActingGasPiston::transported_acoustic_candidate` binds each additional
acoustic face to one distinct finite regional gas inventory. Additional
unbound parcels are finite reservoirs. The existing material transaction owns
nozzle, diffusion, porous and membrane transport over one original snapshot;
each mechanism retains its native species, partial-enthalpy and endpoint
driving-force gates. Acoustic displacement flow is never a material mass flow.

The original shared-wall candidate first resolves all acoustic reactions and
advances every complete load once in its initial medium. A private staged node
then delegates to `dg_port::finite_gas::prepare_after_step` inside the existing
transport transaction. Gas-work trials reuse the same transported inventories
and already-marched acoustic endpoints. Native operator rebinding supplies the
actual parameter work, which is debited from each finite gas. No piston,
acoustic branch, terminal channel or exterior history is replayed by that solve.

Admission requires unique parcel bindings, complete original gas-source
agreement, the actual fixed circular-bore inventory including main-gas terminal
passages, and the original common clock and work limits. Caller input ports
must be empty; only the solved displacement supplies the native inlet flow.
Independent branch gases retain their original separate sources. Every
endpoint remains private through final transport driving-force, extensive
mass/energy and cancellation checks.

`FiniteGasAcousticCandidate<N>` preserves the original native candidate,
regional parcels and complete loads keyed by the original parcel IDs. Its
report retains face order, native transport reports and actual interval force,
pressure and work. Refreshed endpoint pressure/storage and
`AcousticMotionStep::parameter_work_j` use the final gas; the interval
integrals retain the initial medium. The updated acoustic balance subtracts
native parameter work. The enclosing residual replaces the original frozen
acoustic residual with the native gas/transport/acoustic residual, so neither
propagation work, transported enthalpy nor operator work is counted twice.

This is the existing first-order frozen-medium partition on distinct,
fixed-volume, well-mixed acoustic regions. Their updated operators change the
next interval's shared-wall reaction. The regions are separate from the two
moving-volume piston chambers; uncompensated mean-pressure forces require an
authored external force or balanced reference. It does not implement a
continuous heterogeneous through-flow scattering law, changing acoustic
geometry, jet noise or implicit deposition of acoustic loss into gas heat.

Focused regressions exercise native transport/species conservation, independent
extensive energy, one-step clock/history continuation, thermal reaction on the
next interval, admission, late cancellation and exact retry. These tests are
authored but unexecuted: DSR/RCH and the pinned toolchain are unavailable,
and all seven required constellation checkouts are missing. The available
cached rustfmt parsed and formatted the changed sources. This is source-level
evidence, not compilation or numerical proof.

## Transported gas feedback through electrical and finite-heat motion

`transported_electroacoustic_candidate` and
`transported_electroacoustic_network_candidate` use the same finite regional
gas finisher after the existing reciprocal electrical/shared-motion solve.
The original candidate owns current, back emf, both piston gases and the only
mechanical trajectory. In the network path, the same original material-cell
solve owns temperature-dependent resistance, Joule heat and both gas contacts.
No electrical or thermal trial is repeated by regional gas-work feedback.

The current interval retains its original acoustic-medium pressure law.
Native transport, caloric heat and operator work produce the continuing
regional gas and load endpoints; those operators affect the next interval's
current, motion and thermal network. The complete candidate retains the
original circuit and source reports, and its final enclosing gate uses the
strictest applicable native, acoustic and transport energy tolerance.
The source-bound material-piston entry point retains canonical inertia,
specimen receipts and the full accepted motion clock.

In the ordinary electrical path resistor heat remains the original exported
loss. In the finite-network path it has already entered the authored cell
enthalpies and is not charged again. Regional supplied heat is separately
authored; neither the network cells nor either moving piston chamber may be
duplicated in the regional parcel graph. Native transported enthalpy and
gas-to-wave parameter work remain internal transfers.

Regressions cover unchanged current-interval native results, next-interval
gas-to-current feedback, independent complete energy/mass accounts, finite
network continuation and late cancellation/retry. They are authored but
unexecuted because native build prerequisites are unavailable. Cached rustfmt
parsed and formatted the changed sources; no numerical result is claimed.
Fixed geometry, frozen
linear transducer, balanced mean loads and the original first-order thermal
and acoustic-gas partitions remain explicit limits.

## Real-gas feedback through finite electrical and thermal networks

`DoubleActingGasPiston<G>::electromagnetic_network_candidate` and the
material-bound forwarding path admit the original ideal or PR gas owner.
Both chamber density, pressure and caloric laws remain in the original shared
impulse solve. Paired boundary callbacks retain complete `G::State` midpoints;
source-bound contacts consume their actual `FluidTransport` snapshots.

Every finite-cell trial resolves the same material resistance for circuit
current, reciprocal force/back emf and Joule heat. Both gas contact heats debit
their named finite cells once, including shared-cell contacts. Original circuit
storage, cell enthalpy, both gas internal-energy increments and one canonical
wall kinetic energy close against external supply. The source-resolved moving
body retains its density, geometry, inertia and complete accepted port history.

No circuit, heat or momentum solver is duplicated. Existing nested work limits,
source admission, thermal residuals and final atomic cancellation remain active.
Focused dense-mixture regressions check an independent whole-system storage
account, reciprocal electrical work, endpoint resistance and heat deposition,
continuing source-density substitutions and late contact/cancellation refusal.
The transducer remains the existing frozen linear reciprocal port, the solid
inertia is rigid, and each chamber is a homogeneous closed single-phase gas.

## Acoustic loading through finite electrothermal material networks

`DoubleActingGasPiston<G>::electroacoustic_network_candidate` closes the same
acoustic reaction through the original `MaterialNetwork`,
`circuit_enthalpy::driven_candidate`, electromagnetic velocity residual and
paired gas/thermal boundary owner. Each thermal trial uses its own resolved
material resistance for both native circuit current and Joule heat. Acoustic
back-pressure changes that trial's motion, back-emf, gas temperatures and
interface transfers before the finite-cell residual is admitted.

Both native and material-bound candidates retain the original ideal or PR gas
owner. Density, calorics and pressure remain in that owner's shared motion
solve; contacts use its actual `FluidTransport` midpoint. Dense-mixture
regressions independently account for both gas stores, capacitor energy, finite
cell enthalpy, one material-derived kinetic store and original DG storage,
losses and terminal work. They also check exact acoustic replay and continuing
source-bound clocks, including cancellation at the wrapper's final poll.

The additional affine resistance reaches every original heat/momentum trial.
DG histories remain immutable throughout the nested electrical and thermal
searches; complete acoustic loads advance once at the final displacement
flows. Existing total circuit-trial, boundary-solve, heat and gas proposal
budgets remain effective, including rejected trials and final admission.
There is no second circuit, gas, heat or acoustic integration engine.

The original circuit/enthalpy owner deposits resistor heat into its authored
finite cells once. Both gas-contact heats debit their actual counterpart
cells, including shared-cell accumulation. Native non-electromagnetic work
uses actual wall work after removing electrical transfer, allowing the
surrounding acoustic account to cancel only acoustic reaction. The enclosing
balance sums electrical, finite thermal, gas, kinetic and acoustic storage,
native acoustic losses and signed terminal work against external supplies.
It neither adds deposited Joule heat again nor manufactures heat from
unmodeled acoustic loss. The smaller native/acoustic energy tolerance gates
the complete candidate.

The `MaterialDoubleActingGasPiston` entry point retains the original
source-resolved specimen and canonical inertia. Its private binding preserves
the complete proposed dynamics owner and checks unchanged mass and face
geometry. The existing material electromagnetic network path uses this same
binding, fixing reconstruction that reset accepted port clocks and cumulative
history. Subsequent acoustic intervals therefore see the original continuing
mechanical clock, not a newly initialized piston.

All endpoints return together; late source, energy, global-work or cancellation
failure returns none. Zero-area faces reproduce the original finite network
while acoustic history continues freely. The moving specimen's source state
remains frozen; finite thermal cells are explicitly authored counterparts.
Fixed thermal/acoustic geometry, the native first-order thermal discretization,
held interval acoustic flow, uniform closed gas chambers and frozen linear
transducer assumptions remain. This is not a valve jet, a moving-mesh duct,
magnetic saturation or an inferred thermal-deformation model.

Four native G1/G3/G4 regression groups cover independent capacitor, caloric,
finite-cell, kinetic and DG energy over multiple steps; current/Joule/heat
feedback; exact zero-area reduction; source/budget/cancellation refusal and
retry; and two-step source-bound continuation with original port histories.
These tests are authored but unexecuted: DSR/RCH, the pinned toolchain and the
required constellation checkouts are unavailable. Parser/formatter checks
provide source-level evidence only.

## Reciprocal electrical and acoustic drive on one shared piston

`DoubleActingGasPiston::electroacoustic_candidate` composes the original
material circuit/transducer, both finite thermochemical gases, the single
moving body and independently sourced `DgPortLoad` faces. Acoustic pressure
changes the actual midpoint velocity used by the original back-emf residual;
the resulting circuit current changes force in that same solve. The circuit,
gas, momentum and acoustic owners remain the only integrators.

The original affine acoustic relation is prepared once. Its history force
and nonnegative velocity resistance enter each original electrical/momentum
trial, including the existing impulse tangent. Each complete acoustic load
advances once after the electrical residual has converged. No DG state is
marched inside electrical finite-difference stencils, and no propagator is
rebuilt. This is an operation-count property, not a measured speedup claim.

The native circuit retains source work, stored charge/flux energy and
resistor heat. Actual mechanical force work separates authored drive,
electromagnetic transfer and acoustic transfer. Enclosing the participants
cancels the two internal work pairs; original duct storage, physical loss and
signed terminal work remain explicit. The final gate uses the smaller
electrical/acoustic energy tolerance, without energy renormalization.
Resistor heat is an output of this base path, not an implicit finite-body
temperature update.

Every original owner stays immutable until the full candidate is returned.
Clock/face admission, source refusal, bounded electrical work or final
cancellation returns no partial endpoint. Zero-area faces retain the exact
original electrical/mechanical path while their free acoustic histories
continue. Additional acoustic domains are separate from the two finite
chambers; fixed acoustic geometry, held interval flow, frozen linear
transducer factor and the original circuit/gas assumptions remain.

Three focused G1/G3/G4 regression groups exercise independent midpoint RC
feedback, native force and complete energy accounts over twelve steps,
original DG replay, exact zero-area reduction, late refusal and repeatable
retry/continuation. They are authored but unexecuted in this session:
native build prerequisites are unavailable. Source parsing/format checks do
not establish numerical or experimental validation.

## Reciprocal acoustic loads on a shared gas-driven wall

`cavity_phs::piston::double_acting::acoustic` connects additional
signed-area faces of the original translating body to independently sourced
`DgPortLoad` owners. The finite thermodynamic chambers remain the original
two source-backed gas inventories. The acoustic domains are distinct physical
regions; their pressure storage is not a second copy of a chamber.
DG pressures are perturbations. Any uncompensated mean pressure on the added
faces belongs in the prescribed external force or a balanced reference state.

For each accepted load the original mean relation is
`p_i = H_i + Z_i * Q_i`, with `Q_i = a_i * dx/dt`. Its reciprocal reaction
is `F = -sum(a_i*H_i) - sum(a_i*a_i*Z_i)*dx/dt`. The constant history force and
nonnegative velocity resistance enter the original paired gas/impulse solve.
Its exact displacement derivative includes the resistance contribution;
neither a delayed pressure sample nor another mechanical integrator is added.

Each complete DG load advances once at the accepted displacement flow.
The existing propagator, finite branches, wall/diffusion histories, terminal
channels and rational exterior owner retain their own dynamics and admissions.
No acoustic exponential is rebuilt and no DG trial runs inside the nonlinear
gas solve. This removes repeated acoustic marching by construction; no measured
runtime-speedup claim is made.

The candidate retains all original gas/wall and acoustic endpoints until
common-clock, face-count, source, pressure, work, storage and cancellation
gates pass. Reusing the same live load as two faces is refused. Driven exterior-mouth
topologies remain explicitly unsupported by this face coupling. Actual native
inlet work cancels actual mechanical reaction work; original acoustic loss and
signed terminal exchange remain explicit. The enclosing residual includes the
gas and single-wall kinetic account without depositing acoustic losses as
unmodeled gas heat. Zero signed area preserves the original mechanical path
while the corresponding acoustic history still advances on the common clock.

This is a small-displacement fixed-geometry mean-port coupling. Acoustic gases
may differ from each other and from both chambers. It does not turn valve mass
flow into equal volume flow, introduce a moving-mesh bore, or model a
heterogeneous through-flow junction. A rational terminal's signed work is not
an invented storage certificate or an external microphone waveform.

Three focused native regression groups cover independent caloric/acoustic
energy and reciprocal work, reaction-dependent motion and retained terminal
history, exact decoupled reduction, refusal/cancellation and retry/suffix
continuation. They are authored but unexecuted in the unavailable native
workspace; parser/formatter checks are not numerical or experimental proof.

### Real-gas and source-bound material acoustic actuators

`AcousticPistonCandidate<G>` and `AcousticThermalPistonCandidate<S, G>` admit
the original ideal or PR chamber owner through the shared `PistonGas` seam.
Actual chamber density, caloric departures, pressure and acoustic stiffness
remain in the native momentum and thermal proposals. Each additional acoustic
domain retains its own sourced medium and complete history; it need not have
the same gas density or composition as either thermodynamic chamber.

The affine acoustic reaction is prepared once, participates in every native
heat/impulse trial, and each original DG load advances once after convergence.
No acoustic march is nested inside a heat search. Accepted mechanical work
cancels inlet acoustic work; the enclosing account retains both native gas
energy increments, one wall inertia, any finite-body enthalpy, original
acoustic storage, physical losses and signed terminal work.

`material::acoustic` forwards these same mechanical, thermal and source-layer
entries through `MaterialDoubleActingGasPiston<G>`. Its candidates retain the
original rigid-body density and geometry receipts, canonical inertia and
complete accepted port clocks. Binding verifies the native mass and face
geometry without rebuilding gas or acoustic snapshots, then performs a final
cancellation poll before returning any endpoint. The thermal body may be a
separately declared store; heating does not silently alter the rigid specimen.

Focused real-gas regressions cover continuing independent energy/impulse
balances, exact original DG replay, opposite face heats on actual PR
midpoints, zero-area reduction, source-density substitutions and rollback at
the final material-wrapper poll. This extends the closed-chamber acoustic and
finite-body paths; the metering, electrical and evolving acoustic-gas adapters
retain their existing constitutive scope. Fixed acoustic geometry and the
small-displacement mean-port approximation remain explicit.

## Acoustic back-pressure through metering and finite-body heat

`DoubleActingGasPiston::acoustic_valve_candidate` and
`acoustic_thermal_valve_candidate` use the same reciprocal acoustic-load
preparation and finalization as the original closed two-gas acoustic adapter.
They compose the existing metering and finite-enthalpy owners; no alternative
gas, nozzle, heat-transfer or acoustic solver is introduced.

`acoustic_thermal_candidate` connects the closed two-gas finite-body owner
straight to the same loads; it requires no artificial valve. Its
`acoustic_thermal_layers_candidate` variant resolves the existing source-bound
material conduction layers at the actual loaded gas/body midpoints.

Native species and donor total enthalpy are transferred once at the accepted
metering opening. Every subsequent gas/momentum or two-face thermal proposal
uses those same post-transfer inventories. The affine acoustic reaction enters
each original momentum solve, so both gas temperatures, pressure work and the
shared body heat laws respond to the acoustically loaded displacement.
Acoustic, branch, channel and terminal histories advance only once after the
original bounded thermal/momentum solve has converged.

The original valve finalization receives the actual accepted total wall force
and gas heats. Its source-stable caloric/species balance therefore remains
valid under velocity-dependent acoustic reaction. The shared body's actual
enthalpy change encloses both face heats exactly once. The final acoustic
account cancels actual mechanical reaction work and retains native acoustic
storage, loss and exterior exchange. No transported flow work, jet energy,
latent heat or emitted acoustic work is deposited twice.

Returned candidates contain the original valve or thermal-valve candidate
alongside all original acoustic loads and the measured enclosing energy
residual. Source refusal, nested work exhaustion, final pressure reversal,
acoustic refusal or cancellation returns no partially advanced participant.
Callers consume all returned endpoints together. Zero-area loading
reduces to the original native valve/thermal behavior.

The original first-order frozen-opening transport split, fixed finite-body
enthalpy geometry and separate acoustic-domain assumptions remain. There is no
new jet recoil, thermal deformation, moving-mesh duct or heterogeneous mass-flow
acoustic junction. Three focused native regression groups cover both donor
directions, independent species/caloric/acoustic work, sensible/latent body
feedback, exact native reductions and late refusal/retry/continuation. These
tests are authored but unexecuted; source/format checks do not establish
numerical or measured physical validation.

## One finite thermal body on the moving two-gas valve

`DoubleActingGasPiston<G>::thermal_valve_candidate` composes the metering valve
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

`DoubleActingGasPiston<G>::valve_candidate` and `step_with_valve` select the
donor from actual initial pressures, call the original compressible nozzle and
species/total-enthalpy transfer once, then call the original two-gas midpoint
pressure-work and single-wall momentum solve. Different gas compositions,
source caloric laws, temperatures and densities retain their original owners.
Native throat/choking, source-domain, gross donor-fraction and transfer
conservation gates remain active. Final pressure ordering is checked again
after the actual heat and motion.

`MeteredPistonGas` dispatches to the existing ideal or PR nozzle and transfer
owner, retaining its native control and report types. Compatible PR mixtures
transport actual species, with all new pair interactions explicitly sourced.
The subsequent closed motion measures energy from both post-transfer density
and temperature through the retained mixture law. It never applies an ideal
temperature-only caloric increment to a changing-density real gas.

The same generic thermal path supplies complete original gas midpoints to both
contacts. `MaterialDoubleActingGasPiston` forwards valve and thermal-valve
candidates without reconstructing the moving owner: original density claims,
geometry, canonical inertia and accepted port history remain attached. A
separate thermal body's heating does not silently change rigid specimen mass.
Additional dense-mixture regressions cover both native flow directions,
species and independent U+K, continuing position-controlled opening, reciprocal
finite-body heat, 1000/8000 kg/m3 source-inertia substitutions, and late atomic
refusal. These extend the original algorithms; no flow or motion solver is
duplicated.

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

## Finite transport networks drive original moving apertures

`bernoulli_aperture::dg::finite_gas::transport` admits single-coordinate,
independent-blade and material-modal sources through their original
`LoadedApertureDg` owners. They use the same finite-gas transaction as wall-DG
and passive complete loads. Typed source inputs retain prescribed pressure,
body flow, independent forces and original modal pressure tractions; callers
cannot also prescribe the acoustic flow computed by the moving source.
Reservoirs admit caloric controls but no source drive or acoustic ports.

Every material path consumes one original species/partial-enthalpy snapshot.
The original source and complete load advance once in their initial medium,
then the original caloric feedback pays actual bore/channel parameter work
from that transported inventory. Gas changes reach the same source density,
viscosity and exact acoustic impedance without resetting its motion, contact,
material memories, accepted clock, channel momentum or exterior pole history.
Immutable source potentials/reference geometry are shared across previews;
all evolving histories and work accounts are cloned independently.

The enclosing account contains gas plus full mechanical/acoustic storage,
native exported loss, actual heat, prescribed upstream/body/mechanical work
and signed exterior exchange. Valve-to-duct work is internal and cancels;
modal pressure work is already inside original mechanical work. The inherited
`acoustic_storage_change_j` and `acoustic_dissipation_j` fields therefore cover
the complete bound node, including mechanics. The existing final transport
driving-force, mass, energy and cancellation gates publish all parcels and
sources together. An immutable candidate exposes the same proposed endpoints.

Separate gas volumes retain explicit physical material-transport paths. This
first-order fixed-volume partition does not infer thermodynamic mass flow from
oscillatory valve flow, static pressure prestress, jet noise or loss-to-heat
partition. A two-sided acoustic source still needs its own joint inventory
model; it is refused here rather than assigning one parcel to both tracts.

Four public G1/G3/G4 regression groups in
`tests/dg_gas_rebind/finite_transport.rs` exercise species conservation, original
transport equivalence, independent caloric/mechanical/duct/channel energy,
changed subsequent source response, zero-link reduction, and complete rollback
with exact retry. They are authored but native compilation/execution remains
unverified in the current environment. Parser/format checks are not native
numerical evidence. This is an implementation increment under MR05/MR06/MR33/MR67.

## Mixed physical owners in one finite gas network

`dg_port::finite_gas::transport::assembly` admits original wall-DG nodes,
complete passive loads and moving apertures in one parcel-indexed assembly.
`FiniteGasAcousticNode` retains each complete native owner;
`FiniteGasAcousticInput` distinguishes held native ports from moving-source
drives. Plain reservoirs have caloric controls without an acoustic owner.
Every input kind must match its actual owner before material transport runs.

The enum dispatches the existing node adapters into the same transport kernel.
No constitutive, caloric, feedback, conservation or publication algorithm is
duplicated. A wall retains its complete physical port collection; a passive
load retains its one original inlet; an aperture computes its own acoustic
flow. Native gas/source/volume/history admissions and a common physical clock
remain mandatory. Original reports preserve those distinct physical accounts.

Gas and complete node storage, actual heat, exported loss and external work
close once over the entire heterogeneous enclosure. Valve-to-duct work remains
internal, while passive net port work already includes exterior exchange.
One material snapshot drives all owners; final transport driving forces and
energy/mass gates run after every endpoint. Immutable previews and accepted
publication retain every gas, mechanical, material, channel and terminal
history together. The existing fixed-volume first-order partition and explicit
material-path limits remain; it does not create an acoustic junction between
distinct chambers.

Two public regression groups in
`tests/dg_gas_rebind/finite_transport/assembly.rs` use a warm source-backed
elastic wall, a passive finite-channel/rational load, an independently moving
paired aperture and one finite reservoir. They compare each accepted endpoint
to its original owner using the same transported snapshot, check species and
independent caloric/external-work accounts, and exercise input-kind refusal,
immutable preview, final cancellation, exact retry and continued evolution.
These are authored synthetic implementation regressions; native execution and
measured physical validation remain separate gates.

### Complete-load viscous and molecular heat during material transport

`FiniteGasAcousticInput::ViscousMolecularPort` selects the complete-load heat
continuation described above for an original `FiniteGasAcousticNode::Port`.
Admission rejects other acoustic owners and unbound reservoirs before material
transport. Main-gas bore molecular and viscous losses and main-gas terminal
passage shear heat return through the same caloric, chemical and finite-body
transaction. Independent branch and numerical losses retain their original
exported ownership; no independent source is assigned to the main parcel.

Preparation retains one original marched load, its actual implicit port flows,
the native output and both physical heat accounts. Every later transported,
reacting, thermal-body or material-work trial borrows this same proposal; it
does not repeat the acoustic step or its dense source integrals. The original
kernel keeps the prepared owner beside its source heat until final atomic
publication. Plain port and existing wall/aperture inputs preserve their
previous explicit selection and zero-transfer behavior.

The common `viscous_heat_j` account includes bore and main-terminal viscosity,
while the native complete-load report retains those contributions separately.
Internal gas heat is included once in endpoint caloric energy and subtracted
once from exported acoustic loss. It is excluded from prescribed external
heat. Existing species, source, driving-force, parameter-work and cancellation
gates remain in force, including molecular source composition restrictions.

Focused assembly regressions independently integrate slit shear heating while
species transport runs, retain the original channel/exterior histories, and
check caloric/species balance, wrong-owner refusal and late-cancellation retry.
Source/parser checks and an independent bounded review passed; native coupled
compilation and execution remain unverified in this environment.

### Moving-source viscous heat during material transport

`FiniteGasAcousticInput::ViscousMolecularAperture` selects the same original
main-gas heat for an actual moving aperture. Its input retains the original
held drive, blade forces and material-modal pressures; it cannot prescribe an
acoustic flow. Other owners and unbound reservoirs refuse this input before
material transport. The original mechanical and acoustic trial solves that
flow once, and its private prepared owner is retained with its physical heat.

All caloric, chemical, finite-body and parameter-work endpoint trials borrow
the same prepared motion, waves, channel profiles and terminal histories.
The existing transport kernel includes main-bore and main-terminal gas heat
once and removes it once from exported loss. Prescribed external heat remains
separate. Mechanical, contact, jet and independently owned branch losses
retain their original exported accounts; no aperture-gap inventory is inferred.

The focused transport regressions reconstruct terminal shear heat independently,
check unchanged species transport, NASA caloric balance, original paired-blade
motion and exactly one accepted source/acoustic step, and exercise wrong-owner
and final-cancellation refusal with exact retry. Source/parser checks completed
before the execution service became unavailable; recovered source patches were
reviewed and verified through GitHub. Native coupled compilation and execution
remain unverified.

## Reacting finite gas in mixed acoustic assemblies

`dg_port::finite_gas::transport::assembly::reactive` connects the original
implicit `GasReaction` owner to the same mixed wall-DG, complete-load and
moving-aperture transport transaction. Sparse parcel-keyed bindings select
reacting acoustic chambers or finite reservoirs; unlisted parcels retain the
original caloric calls. The reaction catalog must preserve the exact sources
of incoming species and initially absent products. Gas inversion controls
must agree with the original parcel input.

Material transport runs once on the original network snapshot. Every native
acoustic work trial restarts chemistry from that same transported inventory.
The original boundary mechanism solves the endpoint at fixed volume with
`U_end = U_transported + Q_external - W_gas_to_wave`. Species formation
energy stays in gas storage and is never added again as reaction heat. Heat
and mechanical work retain their distinct entropy accounts.

Numerical work probes defer the native chemical entropy inequality; the
actual accepted background report must pass it after acoustic work closes.
Reacting reservoirs use the original boundary endpoint with zero mechanical
work. Final transport driving-force gates and the enclosing native energy
and mass gates still apply. An additional compensated species account checks
actual network amount changes against the retained stoichiometric reaction
changes. Candidates retain every endpoint privately, and final publication
changes all parcels and original physical owners together.

This is first-order snapshot transport followed by coupled endpoint chemistry,
heat and acoustic parameter work in well-mixed fixed volumes. Arriving species
may react during the interval; products enter outgoing transport on the next
interval. It supplies neither resolved flame dynamics nor chemical acoustic
relaxation, and does not invent pressure excitation or an acoustic junction
between chambers. Chemical entropy admission is stage-local, not a claim of
complete assembly entropy production.

Source parsing and a bounded API review passed for this implementation. Native
compilation and execution remain unverified because the available Cargo build
cannot load required sibling dependency sources. Experimental validation remains
a separate gate; no measured combustion or sound-accuracy claim is made.

## Finite material heat in reacting acoustic assemblies

`dg_port::finite_gas::transport::assembly::reactive::thermal` connects distinct
finite material bodies to the same mixed acoustic transport transaction.
Sparse thermal and reaction bindings are independent: a wall-DG volume, complete
passive load, moving aperture or plain reservoir may exchange heat with its own
uniform body, with or without chemistry. The accepted body-state map and each
binding's initial phase state must agree exactly. Original enthalpy and layer
conductivity must come from the same material card, and the contact area cannot
exceed the body's declared exposed area. The original graph budget applies
before any material source is evaluated.

The existing bounded scalar feedback owner solves
`Q = dt G(T_body,end - T_gas,end)`. Every heat and acoustic-work probe restarts
the same accepted body and the single transported gas snapshot. The original
phase chart gives `H_end - H_start = Q_body,external - Q`; the original caloric
or implicit chemistry endpoint receives `Q_gas,external + Q` and separate
gas-to-wave work. `MaterialConductionLayers` resolves its existing series
resistance law at both endpoints' temperatures, retaining the conductivity
claim selected before transport. Heat, formation energy and mechanical work
are never substituted for one another.

The gas/acoustic subaccount includes the actual heat across its gas boundary.
The enclosing account adds measured finite-body enthalpy changes and only the
original external heat sources: the interface Q cancels once. Native exterior
work, storage and exported losses remain unchanged. Local body and gas/body
energy gates, accepted heat direction, final transport driving forces, chemical
entropy, whole-network species and total energy admission all remain active.
Only the accepted work trial supplies the retained body/source report. Preview
does not mutate state; accepted publication changes gas, native acoustic owners
and body states together after a final cancellation check. Returned Q is already
applied and must not be deposited again.

These are distinct uniform thermal stores with fixed effective layer geometry
and negligible solid pV work. They do not automatically change an aperture's
mechanics, wall stiffness, dimensions, or general acoustic-loss partition. A body
shared by several gas volumes uses the separate `thermal::shared` owner below;
duplicate entries in this independent-body map remain unsupported. Phase follows the original chart without claiming melt
motion or structural recoupling. Uniform-body validity and material/density/
geometry consistency remain explicit caller obligations. Transport and frozen-
medium acoustic sampling retain the existing first-order split; source-root
convergence is not a temporal-error estimate or a uniqueness theorem.

Focused G1/G3/G4 regressions use the existing warm passive-load and paired-
aperture fixture, a source-bound finite material store, independent energy and
two-capacity backward-Euler references, exact empty-binding reduction, and
atomic refusal/cancellation with deterministic retry. These are synthetic
implementation checks, not experimental material or acoustic validation.

## One finite thermal body shared by reacting acoustic volumes

`dg_port::finite_gas::transport::assembly::reactive::thermal::shared` connects
several distinct gas parcels to one original uniform material body. One
canonical enthalpy endpoint supplies every contact temperature; the body's
mass, storage change and prescribed heat appear once in the enclosing account.
Each contact retains its own material layers, signed heat bracket and energy
allowances. All layers must use the body's original material card, and their
combined area cannot exceed its exposed surface. Disjoint patch placement and
uniform-body validity remain caller obligations.

The common transport owner now prepares one admitted immutable material
snapshot. Every body, contact-heat and acoustic-work trial reuses that snapshot
and the original accepted native histories. The existing scalar feedback owner
closes the body's energy-out coordinate against the independently evaluated
contact sum and any selected endpoint wall-material work.
Each gas contact closes `Q_i = dt G_i(T_body,end - T_gas,end)`, using the
original caloric or implicit chemistry endpoint and separate gas-to-wave work.
The conductivity claim selected before transport stays pinned at actual
endpoint conditions. Every signed contact is evaluated even when the net heat
is zero, so a hot gas can heat a cold gas through an unchanged-temperature body.

The body's original chart supplies
`H_end - H_start = Q_body,external + H_wall - Q_net - W_wall`. The accepted report retains
the applied heat coordinate, actual gas contact sum and signed wall work as
separate accounts; their difference is a bounded source residual, not an
additional heat source. The physical body-energy gate uses actual storage,
actual contact heat, original wall heat and native wall-material work. The
enclosing gas/acoustic/body account independently closes original external
heat/work, native storage and exported losses. Original transport driving-force,
chemical entropy and species gates remain active.

A whole-call endpoint evaluation budget bounds all discarded body, heat and
work trials, including parcels without a thermal contact. Native owners keep
their own internal limits. Missing sources, invalid charts, budget exhaustion,
unclosed residuals and cancellation leave the gas map, all acoustic histories
and the one body unchanged. Accepted publication changes those owners together;
reported contact heat has already been applied.

This remains a uniform-body, fixed-geometry, negligible-solid-pV-work model
with effective conduction layers and the existing first-order transport split.
It does not supply spatial solid gradients, thermal deformation, aperture
updates or general acoustic-loss partition. Explicit material-wall
feedback and selected molecular source heating use the original owners
described below. The G1/G4
regressions cover opposite nonzero contacts with zero net heat against an
independent constant-capacity solution, reacting transport with continuing warm
passive and paired-aperture owners, and atomic refusal/cancellation/retry. They
are formatter/source checked; native execution remains unverified.

## Finite-body temperature feedback and original acoustic-wall self-heating

Each shared-body thermal port may select `SharedFiniteGasAcousticWallSolve`
to refresh its original uniform or regional wall at the solved body temperature.
Isotropic elastic, Young-viscosity Kelvin--Voigt and source-pair Maxwell laws,
plus principal-aligned orthotropic elastic/proportional Kelvin--Voigt laws,
retain their complete original material bundles, selected claims, card identity
and every non-temperature condition. The original material-wall rebind retains
displacement, conjugate momentum and acoustic history. Each Maxwell arm keeps
its physical spring extension and viscous offset under the original source
correspondence; no arm is recreated, reset or matched by coefficient value.
A relaxation-time-only change affects subsequent loss without instantaneous
parameter work.

Each regional source keeps its original element, thickness, restraint, query
conditions and loss model. Inventory sums the native element length times its
own thickness and circumference; it never charges total duct length to every
region. All sections share the body's one caloric chart, card and endpoint
temperature. This is mechanically heterogeneous wall response without spatial
thermal gradients or an inferred multi-material heat-capacity mixture. Native
regional rebinding transfers the complete source list in one transaction.

Orthotropic refresh reuses the original principal hoop/axial axes and exact
normal-compliance reduction. An explicit proportional-viscosity time retains
its original property role; normal-plane tensor viscosity and generated axial
memory cannot silently become proportional damping. Source positivity and
fixed density are checked at every candidate temperature.

The optional `dissipation_heat` budget selects the original solid Kelvin/Maxwell
source for deposition into this same finite body. After original-owner and
transport admission, the native source integral is evaluated once from each
accepted initial state under its actual held ports. Every body/contact/work
trial reuses it. Gas boundary-layer, molecular, terminal and numerical losses
are excluded. With no selection, wall damping retains its exported meaning.

The outer energy-out coordinate remains `S = sum(Q_i) + W_wall`, while body
storage obeys `delta H = Q_body,external + H_wall - S`. Signed material work
and irreversible source heat have separate accounts. Positive wall work draws
from the body; negative work returns to it. Each endpoint resolves the original
wall at that same body temperature through the existing native transfer.

Reports retain each material receipt, native intervention and deposited wall
heat, together with total wall heat/work and independently measured final
acoustic storage. `SharedFiniteGasAcousticWallStep.materials` explicitly holds
`Uniform` or `Regional` receipts; regional entries follow original element
order. `uniform_material()` and `regional_materials()` expose the corresponding
layout. No first-element proxy, averaged material, or inferred regional heat/work
partition is reported. The transport report remains the preceding gas/acoustic
substage, which exports `D-M`. The enclosing report exports `D-M-H_wall`.
Each original owner's source heat is checked against its own native decrement
minus internal molecular heat and its explicit energy allowance. Remainders
are never clipped. The final enclosing energy gate uses actual post-rebind
storage and the reduced export, so the internal heat transfer cancels once.

Initial material/body temperature, all-solid phase, constant chart/material
density, fixed geometry/thickness/restraint and sufficient shared-body wall
inventory remain required. Each wall has an explicit parameter-work cap and
opted-in heat has an explicit complete-state integration budget. Source-arm
activation/removal, authored spectra without property identities, generated
normal-plane viscosity/memory, off-axis material frames, independent regional
temperatures, multi-card thermal mixtures, phase change and density changes refuse. Thermal expansion and a thermoelastic entropy model
remain outside this bounded path. The original frozen-sample, first-order
split and source domains still apply.

Failures and cancellation publish no partial body, gas, wall or history state.
Zero internal source heat preserves the previous body arithmetic. Authored
regressions cover signed work, self-heating, retained Maxwell memory, changing
subsequent response, source pins, regional inventories and receipts, directional
restraint/viscosity, opt-out and atomic refusal/retry. Formatter,
parser and bounded source checks are distinct from native runtime validation;
native compilation and execution remain unverified.

## Molecular acoustic heat in reacting thermal gas assemblies

`FiniteGasAcousticInput::MolecularWall` opts an original wall-DG gas volume
into native molecular relaxation heating within the existing heterogeneous
transport, chemistry and finite-body APIs. Uniform and regional material wall
owners retain their original histories. The selection requires the existing
`UniformMolecularDiffusion` refresh and an explicit native dense-exponential
budget. Binding it to a complete passive load, moving aperture or unbound gas
reservoir refuses, including when the expected heat would be zero.

The prepared transaction integrates the original molecular slip-square source
once from each selected accepted duct/state and its actual held ports, after
all owner admission and before material transport. Every subsequent body,
contact-heat and acoustic-work trial reuses this same source integral M. The
common caloric/chemical/thermal endpoint receives actual gas heat
`Q_external + M + Q_body`. Native gas-to-wave parameter work remains separate.
An internal callback seam carries M to its original native energy account;
there is no public override accepting an arbitrary deposited molecular heat.

Reports retain total native decrement D, internal molecular heat M and the
exported remainder D-M as distinct quantities. The gas/acoustic boundary heat
contains only external heat and body contact heat. Its balance therefore uses
`delta U_gas + delta E_acoustic + (D-M) - Q_boundary - W_external`.
The enclosing finite-body account adds body storage once and cancels each
body contact. Unselected wall, boundary, viscous and numerical losses keep their
existing exported meaning; an aggregate energy deficit is never converted
into molecular heat or clipped to a nonnegative remainder.

Ordinary inputs preserve their previous zero-transfer arithmetic. Missing
sources, unsupported owners, exhausted source budgets, invalid native refresh,
failed energy/species gates and cancellation publish no partial gas, acoustic
or body state. Molecular source observations are not repeated by nested roots.

The original molecular refresh still preserves species fractions and reservoir
support. Same-composition transport and changing thermodynamic conditions are
supported; a composition-changing reaction or inflow in the molecular source
duct refuses. Other owners in the same network may react, and the shared body
can transfer molecularly generated heat to them. This remains the original
fixed-geometry, frozen-medium sample with a first-order background split,
not a model for acoustic species perturbations or arbitrary loss partition.
Focused regressions exercise the native source/heat receipt, complete energy
enclosure, coupled body/chemical feedback and atomic refusal/retry. They are
authored and source/parser checked; native execution remains unverified.

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

`finite_solid::impulse` connects this same accepted solid to instantaneous nodal
and material-point impulses. `FiniteSolidPoint` uses original P1 barycentric
coordinates; its velocity trace and transpose impulse map are power-conjugate.
`point_directional_inverse_mass_per_kg` solves against the original consistent
reference mass, with stationary support elimination, so a contact consumer can
compute its effective mass without an unrelated lumped approximation.
`point_impulse_candidate` and `nodal_impulse_candidate` apply the resulting
velocity jump through the existing bounded native mass solver. Point queries
require no whole-mesh allocation; impulse solves retain sparse mass actions.

`momentum_about` integrates linear and angular momentum using consistent nodal
momenta and **current deformed** coordinates about an explicit world origin.
An impulse candidate retains all applied and support impulses, their angular
moments, actual representable momentum changes, and average-velocity work.
Absolute linear/angular/work gates refuse unresolved transfers; inputs remain
unchanged on late rejection or cancellation. Displacement, source/material
state, acceleration, clock and accepted-step history remain unchanged because
the event has zero duration and the admitted conservative force is independent
of velocity. No impact/friction law, event location or impulse adjoint is inferred.
The caller supplies a partner's equal/opposite impulse at the SAME physical
point. `tests/finite_solid/impulse.rs` checks an independent tetrahedral mass
inverse, off-center impulse/energy, fixed-support reactions, origin changes,
native continuation, rigid/solid sticking impact and atomic refusal.

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

The five original fixtures cover independent energy/stress/tangent formulas, free
expansion under rigid rotation, source-driven finite motion, heterogeneous
held temperatures, original mass, and source/domain/Jacobian refusals.
These authored numerical regression fixtures do not supply experimental
calibration.

`ExpansionThermalReference` additionally consumes the original
`fs-material::thermal::transport::expansion::PreparedFreeExpansion` path.
It retains the source engineering-length convention and exact datum
normalization `s_i=(1+epsilon_i(T))/(1+epsilon_i(T_ref))`. Isotropic sources
must be explicitly declared isotropic; directional sources require their
original principal-frame identity and a proper rotation into a named mesh
frame. Source and target frames with the same identity require the identity
map. Unsupported evolving source coordinates, contradictory temperature
coordinates, foreign material cards and extrapolation refuse.

The directional finite reference is `Fth=Q diag(s_i) Q^T`, `H=Fth^-1` and
`Jth=product(s_i)`. The native law sees `Fe=F H`; its response per original
reference volume is `W0=Jth We`, `P0=Jth Pe H^T`, and the exact tangent
`A0[iJ,kL]=Jth H[Ja] Ae[ia,kb] H[Lb]`. Reference tensors are prepared once,
the two tangent contractions are separated, and successful mechanical
evaluation performs no allocation, source lookup, interpolation or inverse.
Both existing mean and new engineering-curve references retain their scalar
fast path when explicitly isotropic. `FiniteSolidMaterial::resolve_thermal`
accepts the explicit `ThermalReference` source variant; the original
`resolve_mean_thermal` constructor uses that same binding owner. The
`thermal_reference` accessor now exposes the source variant, with shared
reference/current contexts, principal stretches and relaxed volume.

The new `tests/finite_solid/expansion.rs` fixtures exercise an independent
Ogden-invariant energy oracle, all nine Piola components and all 81 tangent
components, native tetrahedral assembly, rigid rotation of free expansion,
cooling and negative directional expansion, restraint stress, original mass,
nonzero source-datum rebasing, isotropic reduction and source/frame refusals.
They are numerical constitutive checks, not measured material validation.

This is thermal preparation followed by mechanics at held temperature per
element. It does not evolve a heat/entropy equation or account for work during
a temperature change. Fixed-principal thermal expansion is independent of
the native elastic-law symmetry; evolving principal axes, coupled evolving
finite thermomechanics, inelasticity and nearly-incompressible stabilization
remain open. MR58 and MR21 remain open.


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

`PreparedHeatNetwork::step_piezoelectric_candidate` uses this same ceramic
driver with the complete original thermal assembly. Source-integrated layers,
thermal contact resistance and closed gray-diffuse radiation are re-evaluated
at the same trial temperatures as the ceramic coefficients and resistors.
The accepted report retains the full `HeatNetworkSample`. Closed-radiation
net roundoff remains visible in `closed_radiation_net_j`, is excluded from
external heat, and enters the same physical total-energy acceptance gate.
Geometry, contact areas and view factors stay fixed during the small-signal
step. Electrical interfaces refuse because this rod/circuit driver does not
own their work and two-sided Joule heat. New G1/G3 regression definitions
exercise two finite bodies against independent thermal/mechanical balances
and data-only contact/emissivity changes. Native execution of this extension
remains unverified in the resource-constrained environment.

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
tetrahedron uses the source-backed affine thermoelastic potential. Stiffness
and expansion independently select isotropic or oriented orthotropic source
models, with independent affine engineering constants and instantaneous
principal expansion coefficients. The native compliance inverse and original
source-integrated free strain supply the full tensor thermal derivatives.
The same source/card/eigenstrain and delivered reference mass remain pinned. It admits
a fully solid constant-density interval within the original complete caloric
chart, preserving later latent/liquid segments. Any latent segment intersecting
the admitted interval, or any actual phase/source-domain exit, refuses. The
interval's actual caloric capacity plus the conservative fixed-total-strain
field capacity must stay positive.

The same tensor potential now serves the original static, Newmark inertial,
finite-heat-network and local adiabatic/prestress consumers. Elastic and
expansion axes remain fixed and may differ; temperature gradients retain their
individual physical cells. Whole-interval compliance and capacity bounds use
outward-rounded coefficient and matrix bounds and may conservatively refuse
otherwise valid wide intervals. No endpoint-only positivity claim, isotropic
replacement or new mechanical/thermal integrator is introduced. The original
isotropic path and endpoint thermal coupling convention remain intact.

Original independently affine bulk/shear cards now enter this same reciprocal
owner through `ThermalStiffnessModel::IsotropicBulkShear`. Their direct tensor
and exact first/second temperature derivatives feed the existing internal
energy, entropy, capacity, native stress/energy parity and reversible-heat
accounts. No derived E/nu source or reference-mass change is inserted. The same
binding serves equilibrium, inertial, acoustic and complete heat-network
continuations. Focused free/restrained heating and Newmark work/heat regressions
are in `tests/solid_enthalpy/thermoelastic/bulk_shear.rs`.
`ThermalStiffnessModel::IsotropicYoungShear` now enters the same owner using
the original affine E/G claims and exact rational stiffness derivatives,
including the second derivative needed by elastic heat capacity. Admission
requires positive E, G and `3G-E` over the full interval, the declared strain
envelope and positive combined heat capacity. The same continuations and
regressions cover both source schemas without synthesizing bulk or Poisson
claims; native regression definitions do not imply measured validation.

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

The same original momentum owner accepts prescribed **axial** no-slip wall
velocities through `step_with_wall_velocities`, in lower/upper order for a
slit and `[0, outer]` order for a circular passage. The axis is not a wall;
nonzero axis drive refuses. Boundary shear supplies the original implicit
finite-volume equations, without changing the cross-section, fitting a drag
coefficient, or replacing the transverse field. The factorization and trial
storage are reused, so the step adds no allocation. Stationary-wall stepping
retains its original arithmetic.

`SlitChannelStep` separates signed `wall_work_j`, both impulses **on the
fluid**, nonnegative relative-shear heat, pressure work, kinetic storage,
exit work and backward-Euler numerical dissipation. Mechanical walls receive
opposite impulses. A moving wall can extract work from faster fluid; that
work is not reclassified as heat. `pressure_load_with_wall_velocities` and
`FixedChannel::preview_at_flow_with_wall_velocities` condense the same drive
and history for a static-pressure terminal. Their offset includes the
specified wall source; any coupled acceptance must retain its separate work.
No boundary velocity or trial field publishes on refusal. Gas interventions
continue to preserve cell momentum through the existing material owner.

G1/G3 `channel_wall_motion` checks independent Couette profile/traction/heat,
adverse-pressure second-order spatial convergence, startup sine-series
first-order time convergence, momentum reaction, signed work extraction,
both geometries' pressure condensation and rollback, and circular viscous
startup towards uniform axial translation. This is fully developed laminar
Newtonian flow at fixed gap/radius. It does not supply normal wall motion,
rotating-cylinder flow, swept volume, entrance development, thermal evolution,
or transition prediction. Reference: NETL MFiX FLD05 steady Couette verification,
<https://mfix.netl.doe.gov/doc/vvuq-manual/main/html/fluid/fld-05.html>.

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
The material-modal DG grammar also accepts this same record, or
`supply slit LENGTH_M WIDTH_M HEIGHT_M CELLS`, immediately before
`slit_length`; omission or `supply none` retains its original behavior.
`dg::supply::inlet_state_with_supply` and
`LoadedApertureDg::new_loaded_with_supply` bind the passage to the original
single, independent-pair or full material-modal source, with an optional
upstream DG tract. The existing series elimination combines each exact DG mean
port with the channel's affine pressure law before the sole nonlinear valve
solve. Downstream body flow remains separate: only the actual valve through-flow
crosses the upstream tract and channel. Pressure, volume flow and work are
checked before all original candidates publish; a refusal retains every blade,
solid history, transverse velocity, acoustic state and terminal memory.
`DgApertureFrame::supply` retains the original channel account. Its kinetic
storage, wall heat and separately identified numerical dissipation enter total
accounts once. Source gas changes at construction reach both channel inertia
and no-slip diffusion. Live gas/finite-heat APIs refuse an attached supply
passage until its additional inventory and parameter work join those specific
transactions; they never silently drop its flow or energy.

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

### Prepared source coherence through the loaded exterior

`RadiatingNetwork::prepare_sources` explicitly supplies the complete source and
scalar-load partition of the original graph's non-opening ports. Prescribed
flow, finite scalar pressure generators and a full reciprocal source bank can
coexist. The returned `PreparedRadiatingNetwork` borrows the original geometry,
regional gas records and source/load definitions; source coordinates retain
their canonical labels and distinct open-circuit Pa versus graph-inflow m³/s
units. The original single-source driving-point partition is replaced explicitly.

One original batch BEM factorization supplies every unit-opening solution.
Those same retained solutions produce BOTH the complete mutual radiation load
and the existing receiver pressure/gradient rows; the harmonic observation owner
checks their exact source traces, surface, frequency and medium before reuse.
One `PreparedNetworkSources` boundary factorization supplies all source columns,
including the retained-wave route at an open-port pole. Cold composition forms
source-to-receiver p/u columns for coherent sweeps. No observation refactors a
boundary, substitutes internal pressure for microphone pressure, or adds a loss.

Factor-major source second moments satisfy C=B Bᴴ in the declared mixed SI
units. Each factor is mapped through the solved OUTWARD opening flows before
forming receiver products. Independent generators still create correlated
openings; the direct-column path preserves those correlations. Original owners
compute complete network p-p/U-U/p-U moments, signed source/interface work,
Hermitian internal-source dissipation, mean-square exterior pressure, velocity
cross moments, and active/reactive intensity. First-order microphones project
actual near-field p/u and retain complex cross-receiver spectra. The common
peak-phasor time average 1/2 occurs once; none of these are per-Hz densities.

Coherent and statistical buffers borrow their exact preparation and hide all
outputs after refusal or cancellation. Applications allocate nothing. Sources,
openings and factor rank are capped at 32; original duct and harmonic work
budgets remain explicit. Source-to-opening factor mapping adds at most 32³
complex products, polled between factors. Cold field storage and composition
work also fit the original exterior allowance. Exterior surface power is the
signed opening p/U work already present in network external-load power, never
a second debit or a claimed enclosing far-field power integral.

Focused regressions define independent assembled BEM field parity, an analytic
antisymmetric source mode with zero exterior radiation and positive internal
loss, mixed-unit dense covariance and microphone cross-spectrum references,
legacy flow parity, exact half-wave pole preparation, and transactional refusal.
This remains fixed-geometry, homogeneous-exterior, rigid-wall harmonic physics;
spatial refinement, measured validation and causal time-domain realization are
not established by these algebraic compositions or compiler checks.

### Loaded generator directivity and total outgoing power

`RadiatingNetwork::prepare_sources_with_radiation` additionally prepares the
original BEM angular bank while the same unit-opening solutions that produced
the complete mutual radiation load are still retained. The ordinary source
preparation has no angular cost. `PreparedRadiatingNetwork::far_field` returns
a view bound to that exact loaded network, exterior gas and source partition;
no repeated observation solves a network or boundary problem.

Canonical pressure/flow drives and coherence factors map through the original
loaded outward-opening transfer before complex angular fields are summed and
squared. Independent generators retain the induced cross-opening coherence;
diagonal opening powers cannot replace their factor columns. Peak exp(-iωt)
amplitudes F have units Pa m, with p ~ F exp(ikr)/r. Directional time/ensemble
power is E[|F|²]/(2ρc) in W/sr. Arbitrary observation directions remain separate
from the positive full-sphere quadrature ladder. The original owner requires
two successive grid refinements to agree for each actual drive's total power,
first directional moments and second directional moments. Quiet cancellation
is evaluated in the complex fields, not a rounded source-power Gram matrix.

Angular kernel work and storage have their own explicit cold-preparation
bounds. Repeated calls allocate nothing; their one harmonic work allowance
includes both the at-most-32³ source-to-opening products and the original
angular contraction. Cancellation remains bounded between mapped factors and
angular work tiles. After an error, the caller discards its output slice and
receives no converged integrated result. Zero coherence factors are silence.

Radiated power is an outgoing-field observation, never an additional network
loss or a repair of signed surface work. Its angular refinement estimates do
not establish continuum, panel-integration, port-basis or boundary-solve
accuracy. Focused G1/G3 regressions define a separately assembled physical
boundary solve and independent Simpson/azimuthal angular integral, mixed-unit phase
and amplitude scaling, a dark generator mode with positive internal loss,
partial coherence, and a case where dropping induced opening correlations
changes the predicted radiation. Budget, cancellation and silence refusals use
the same native owners.

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
existing source-pinned independently affine E(T), nu(T) and instantaneous
alpha(T) potential, finite caloric charts, stationary supports, fixed reference
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

The affine-coefficient extension uses this same quasistatic and inertial
assembly. Its focused `varying_coefficients` regressions compare free/clamped
finite heating and heterogeneous compatible stress/displacement with independent
Lamé-coordinate energy and thermal roots, then exercise an actual inertial step.
The full-interval combined-capacity gate includes nu-induced stiffness curvature
even at zero expansion. These added regressions passed source/parser and diff
checks; native execution of the extension remains unverified in this worker.

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

The single-patch and shared-network rigid consumers, and the Hertz/adhesive
island heat consumers, use the native impulse-based `kinetic_change_j` in their combined
energy audits. Rounded total kinetic energies remain diagnostics: subtracting
them can erase finite frictional work under a large common velocity. The
common-normal-velocity regression retains the same slip history, temperature
feedback, body heat and combined residual without enlarging any tolerance.

Native physical tests cover sticking with zero heat, partial slip with only
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

`tribothermal::normal::material::MaterialNormalModulus` binds the original `ThermalHertzPair` to `CaloricNormalContact`. Every thermal iterate re-queries each solid at its own temperature and uses the original nonlinear reduced modulus, gradient and Hessian in native normal response and reversible heat. Both ordered caloric charts must share the exact original material cards, external query context and density claims; the normal owner then retains those exact caloric objects. Density matching uses the solid projection's actual named density requirement and selected claim at every caloric knot, so an equal-valued alternate density source cannot silently replace it. The reference normal request must match the native reduced modulus and limiting yield strength. No representative bulk temperature, reset of normal storage, coefficient-work heat credit or changed reference mass is introduced.

Three focused native G1/G3 integrations, run alongside the seven original normal cases, verify source-backed loading against an independent rational-compliance caloric root, heating at held indentation with changing contact INTERNAL energy, and original-card/context/density/domain refusal with exact replay. In the authored loading reference, initially 300 K / 360 K bodies reach approximately 288.238949147 K / 360 K; mechanical work is 6.832367020 J, contact internal energy increases 8.008472106 J and physical damping heat is zero. The difference is reversible caloric cooling. These are synthetic implementation references. The selected per-solid Young/Poisson, Young/shear or bulk/shear coefficients follow the admitted affine source paths. Fixed yield, fixed geometry, finite source-domain and first-order coupling limits remain in force. Added native regression sources exercise bulk/shear material heating at held indentation with independently reconstructed temperature, force and internal-energy balance, named-density retention and alternate-density refusal. They retain the original accepted state until the complete coupled candidate succeeds.


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


### Native SDF solid friction and temperature feedback

The opt-in `experimental-solid-contact-heat` feature exposes
`friction_heat::solid::propose_solid_contact_heat`. It solves the existing native
`ContactTrajectory` and `FrictionHeatNetworkStep` at one common endpoint using
the existing bounded Newton owner. Every thermal trial starts from a clone of
the same accepted displacement, physical time, sliding history and work. Its
solved common face temperature refreshes each original `SourceThermalCoulomb`
primary temperature coordinate. Native contact still owns pressure, secant
speed, accumulated distance, normal load and force equilibrium. Additional
temperature coordinates, including semantic fixed-external temperatures and
non-primary legacy `T`, refuse until their thermal body/side ownership exists.

One bounded exact-temperature cache avoids repeating the finite-element solve
when thermal-cell perturbations leave all interface temperatures unchanged.
Accepted mechanics, load, motion and source bindings are immutable within the
call. The cache holds one native result and trajectory, with no Hessian or
thermal/source callback payload. Conductivity and other heat resolve on every
residual; actual mechanical solves and thermal trials are reported separately.

Each active node has one explicit original thermal route. Initial separation,
loss of the mapped contact, duplicate/missing coverage and pressure-width times
depth inconsistent with thermal face area refuse. Native nodal Huber work is
converted from J/m to J using a positive supplied physical depth. A declared
per-node heat fraction enters the massless thermal face; the remainder is
reported as exported nonthermal work, without inventing a damage or wear store.
Original material conductivities and finite caloric bodies determine the signed
two-body heat partition, including ordinary hot-to-cold exchange. No side split,
new heat capacity or second conduction model is introduced.

The exact final nonlinear payload retains the mechanical solution, source
queries and trajectory; trial callbacks never publish accepted history. Nodal
work closure, original cell/interface/thermal gates and the independent port
balance `W_body_friction + delta(H) + D_other - W_support - Q_external` all precede
publication. This is a friction/thermal port balance, not an integrated bulk
hyperelastic energy assertion. Source/domain, nonlinear-budget and cancellation
refusals leave both original owners unchanged.

This remains fixed-footprint, quasistatic plane strain with Huber creep and
endpoint work/path quadrature. Contact transitions and tangent reversals require
time refinement; no continuous-path contact certificate is added. Bulk thermal
strain/softening, moving thermal coverage, finite interface resistance and flash
temperature remain outside this entry. The native integration fixture in
`tests/friction_heat/solid.rs` compares heated friction and two finite bodies to
an independently eliminated scalar reference, checks final T/path queries and
late rollback, and exercises explicit nonthermal loss and ownership refusals.
The two original work-conversion kernel tests pass; the new integrated fixture
has not run under the current executor resource limit.

### Deformable damaged adhesive joints and finite-body heat

`friction_heat::cohesive::propose_deformable_cohesive_heat` advances the original
`CohesiveMotion` once and routes its accepted per-patch physical friction and
surviving-bond viscous losses into the original `FrictionHeatNetworkStep`.
Every source-bound mechanical patch has exactly one distinct, matching-area
thermal face. Original patch identity maps to canonical thermal link order;
several faces may share a cell without duplicating its capacity or heat.
The caller declares the fixed-reference face coincidence and original ordered
surface/material correspondence, which a reduced cohesive jump map cannot infer.

Both ends of each native straight normal-jump path must remain closed. The
existing thermal route retains original geometry, materials, masses and surface
coverage, checks its perfect-contact `k A/d` construction, and suppresses the
replaced direct link once. Common face temperatures and signed hot-to-cold
exchange are solved from endpoint conductivities and finite caloric laws. There
is no prescribed heat share or added interface capacity. Fracture work stays in
its original nonthermal account, and exact static friction contributes no heat.

Mechanical coefficients are frozen for this split increment. Thermal Newton
trials reuse the actual mechanical losses rather than rerunning the mechanics
or its damage update. Source-temperature transfer remains explicit between
increments; timestep refinement is required for thermal-feedback accuracy.
Sub-resolution numerical enthalpy increments use the original resolvable-iterate
projection while the physical residual still charges every source and flux.
The native mechanical, thermal, per-face and combined
`delta(K + U_bulk + U_bond + H) + D_fracture - W_applied - W_prescribed - Q_other`
gates precede publication of both candidates. Late source, opening, domain,
event, work-budget and cancellation refusals preserve both accepted owners.
Moving coverage, bond thermal resistance/storage and simultaneous mechanical
source-temperature feedback remain outside this fixed-reference split.

Four native regressions in `tests/friction_heat/cohesive.rs` compare scalar
momentum and independent finite-volume temperature balances, separate fracture
from friction during growing damage, permute two unequal patch losses sharing a
finite cell, and check static zero heat, continuation and late rollback.
These are manufactured analytical references, not specimen calibration.

`friction_heat::cohesive::interval::propose_deformable_cohesive_interval_heat`
completes an original bounded mechanical interval, including unaligned stopping,
static recapture and subsequent breakaway, before its accepted per-patch losses
enter one original finite-body heat solve. Constant reduced forces and affine
prescribed jumps retain their actual work. Initial and every accepted native
straight local normal path must remain closed; final reclosure cannot hide an
earlier opening. The native and thermal endpoints retain the exact requested
common clock, and rejected stopping trials supply neither heat nor path bounds.
Interval-average physical powers feed endpoint backward Euler, so refinement
is required to resolve heat timing and source feedback. Fracture remains outside
heat. An independent whole-interval energy gate, plus native and thermal audits,
precedes return of either candidate. Four regressions compare unaligned stopping
with analytical kinetic loss and finite-volume temperatures, static hold and
breakaway, two growing-damage patches sharing finite thermal capacity, late
thermal refusal/cancellation, and opening followed by reclosure.

`friction_heat::cohesive::cooling::propose_deformable_cohesive_cooling`
connects original native inertia/deformation to `FractureCooling` through actual
debonding, detached motion and reclosure. Every original patch must have the
same ID, exact reference area, immutable law and complete accepted point in
both owners; mechanical, thermal and last-synchronized clocks must agree.
The caller declares the fixed-reference face coincidence and ordered thermal
cell/body/material correspondence, which the reduced jump map cannot infer.
An exact target clock supplies the same represented interval to native mechanics
and timed cohesive loss replay. The original structural and thermal bounds apply.

The declared staggered order is one native mechanical proposal, native thermal
transport with OLD closure, then endpoint contact/history synchronization and
physical heat deposition. Opening/reclosure changes the NEXT thermal interval.
The original thermal owner retains caloric/latent storage, scheduled external
inputs and births, and independent bonded/reclosed conductances. Reclosure
cannot heal the failed bond. Every patch requires explicit, independent fracture
and bond-viscosity heat fractions and receiving-side splits; each unconverted
remainder remains in the report. Native mechanical and timed material replay
losses must agree, and only NEW losses are deposited once. Nonzero residual
friction refuses until its own receiving/heat owner is declared.

An independent compensated endpoint balance reduces original per-cell enthalpy
changes and `delta(K + U_bulk + U_bond)` plus nonthermal fracture/viscosity, minus
applied/prescribed work, incoming birth energy and external thermal input. It
must agree with thermal-transport plus loss-transfer residuals minus the native
input-oriented mechanical residual. Converted heat is internal to this combined
system and is not deducted again as external input. Every native and combined
gate, including final cancellation, precedes publication of both candidates.
The source coefficients remain frozen: thermal source refresh, event localization,
moving contact area/frames and strong thermal feedback remain explicit separate
operations. Refine coupling times for the stated first-order split.

Four focused regressions in `tests/friction_heat/cohesive/cooling.rs` cover an
analytical inertial rupture with independently integrated fracture/viscous work,
finite-cell heat, nonbinary synchronized clocks, continued detached motion and
reclosure with changed conductance; native latent storage and zero conversion;
foreign law/history/area/clock and unowned-friction refusal; and late caloric
domain/cancellation rollback with exact retry. These are manufactured analytical
references, not specimen calibration.

`friction_heat::cohesive::cooling::temperature` closes the explicit source
feedback between those two original owners. Preparation samples each patch's
actual caloric-cell temperature with an explicitly declared A/B weight and
refreshes its original immutable card through `at_temperature`. Each changed
patch performs one source query; an unchanged temperature validates and retains
the exact source clone. Inputs and returned sources follow original mechanical
patch order, checked by IDs. No new interface, body temperature average, claim
selection, nonthermal coordinate, model envelope or material limit is invented.

The caller prepares the existing native `CohesiveSystem` from those sources.
`propose_deformable_cohesive_temperature` then privately transfers the same
accepted fixed-jump history through BOTH native temperature owners. Current
thermal, mechanical and synchronized-history clocks must agree; source/card,
original area, law and complete old point must match before either transfer.
The original mass, bulk elastic operator, displacements, velocities, trajectory
count and mechanical controls remain with the native motion owner. The cooling
owner independently verifies the new source temperature against its actual
caloric cells and changes the same bonded/failed thermal-contact topology.
Complete resulting point histories and both native work/storage/fracture reports
must agree before either result is returned.

An independent compensated sum of original-area changes in represented bond
storage and cumulative fracture, minus signed parameter work, obeys the shared
joule gate and agrees with both native residuals. The two native reports are
observations of ONE physical transition, never two additive losses. Parameter
work needs the caller's declared physical supply/receiver, and new thermal
damage remains nonthermal; no heat, time, displacement work, thermal expansion
or strong implicit thermodynamic closure is inferred. Thermal weakening can
open a contact at fixed gap; cooling cannot heal it. The next ordinary coupled
cooling step consumes the changed traction, viscosity, damage and conductance.
Four focused native regressions cover an actual heating/weakening/cooling cycle
followed by independently calculated inertial motion, mixed/viscous source and
mode history, exact no-op preparation, stale/foreign/expired bindings, changed
mechanics, and late cancellation with a physical retry. These integration
regressions are authored; full execution remains blocked by the dependency-build
resource limits in this environment.


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

### Explicit attached Maugis branch on native mechanical coordinates

`BoundNormalContact` also accepts `MaugisAttachedPotential` through the original
`NormalLawSource` ingress. The original contacting Maugis equations supply
negative-indentation tensile forces and absolute elastic/surface storage.
Body masses, affine contact rows, inner mechanical storage and the native
`PortHamiltonian`/`StepWorkspace` remain the existing owners.

For this explicit branch, `BoundNormalStorage` supplies a complete native
`ConstitutiveDiscreteGradient`: the inner owner supplies its own rule or an
exact quadratic Hamiltonian, and each supported adhesive contact supplies its
original path-mean force. Each force is projected through the same contact row
and never used to correct another body's momentum coordinate. Unequal masses
therefore retain pairwise action/reaction. Attached composites require complete
inner work and original contact paths; incomplete inner storage refuses instead
of silently reverting to a global work correction. Native numerical coordinate
probes use each retained contact's physical displacement scale. The original
radius-coordinate path derivative also supplies the native Newton action,
including the explicitly admitted continuous noncontact branch.

Quadrature and original-source refusals produce nonfinite trial efforts in the
infallible native interface; the source APIs retain typed errors. A failed work
evaluation cannot become an optional-numerator fallback. The prepared solver
publishes no partial endpoint on refusal. No phase switch is inferred under
the default contacting-only policy. No snap impulse, damping,
heat, or evolving material law is inferred. The actual source branch
and envelope are checked at every endpoint, midpoint and quadrature sample;
mechanical timestep refinement and the local quasistatic contact assumption
remain necessary. Focused `adhesive_maugis_storage` tests cover independent
original-coordinate force work, source identity, tensile motion and reversal,
unequal-mass energy/momentum, prepared/reference agreement, tiny paths and
source/work-budget refusal without state mutation.

The potential's explicit `with_noncontact_continuation` policy also admits the
original stable deformed noncontact disk and continuous loss/recovery of solid
contact. Native bodies can move from positive adhesive gap through the
zero-area junction into compression and back. Original disk work is analytic;
only the contacting portion of a crossing path consumes quadrature. Gap states
retain attractive normal force but exactly zero solid friction area. The same
complete native discrete gradient preserves pairwise action/reaction and the
original absolute energy, without an impulse or loss at the continuous
junction. The original terminal disk fold and unsupported high-Tabor
continuation refuse without publishing a native state. Additional focused
regressions check independent disk integrals, path reversal/subdivision, tiny
paths at the fold, and inertial gap/contact/gap motion with unequal body masses.

### Explicit attached JKR branch on native mechanical coordinates

`NormalLawSource::JkrAttached` also drives the existing `BoundNormalStorage`
through the complete native contact-wise discrete gradient. Its original
high-Tabor source supplies signed tensile/compressive force, absolute
elastic/surface storage and an allocation-free polynomial path mean. The
original displacement scale sets numerical coordinate probes. Each effort
scatters only through its physical contact row; unequal-mass action/reaction
and inner kinetic/modal work are retained. Missing complete inner work refuses.

The original JKR source admits every endpoint and retains finite force/energy
at its fold without a smooth tangent. Fold crossing and source-envelope escape
refuse without publishing a prepared state. Local contact is quasistatic and
the caller owns body inertia; snap kinetics, branch selection, elastic waves,
evolving material state and jump heat require separate owners. Focused
`adhesive_jkr_storage` regressions cover original-source identity/work, tiny
paths, tensile acceleration and reversal, unequal-body energy/momentum,
prepared/reference agreement, cancellation/retry and atomic branch refusal.

JKR assemblies, including mixed Hertz/Hunt--Crossley, DMT and retained
finite-gap contacts, supply an allocation-free analytic derivative of this same
path effort when the complete inner native rule supplies its derivative, or
when that inner rule is the exact quadratic energy. Each original scalar
endpoint derivative projects through its contact row's
outer product. Valid near-fold motion therefore does not require central
numerical probes outside the attached domain. The coincident fold retains its
singular derivative and refuses; no finite tangent is invented. A native inner
rule without a derivative is never replaced by its quadratic view. Maugis uses
its original radius-coordinate analytic-slope quadrature and closed-form disk
derivative, retaining the source work budget, zero-area crossing and terminal
fold refusal. No endpoint-force subtraction or finite-difference contact
probe is required. Focused checks compare the
complete analytic action with directional work differences and advance a
strictly attached near-fold state through both native and prepared execution.

### Mixed adhesive and unilateral assemblies

The same complete native rule admits all existing `NormalLawSource` owners in
one assembly: JKR, Maugis, DMT, Hertz/Hunt--Crossley reversible storage and retained
finite-gap response curves. The inner mechanical owner still supplies its own
complete rule or exact quadratic energy. Each footprint scatters its original
scalar path force only through its declared physical row; no common power-law
fit or correction to unrelated momenta is introduced.

Work paths retain the raw signed affine approach, including the inactive
fraction before a unilateral patch closes. Hertz storage uses its original
stable polynomial mean; finite-gap curves integrate their retained linear
segments with compensated accumulation and unchanged estimate authority. The
ordinary energy/force observations retain their existing unilateral clipping.
Original source support, table continuity and bounded table size are checked;
failure cannot become an optional work fallback. Native probe scales also
cover coordinates connected only to Hertz, DMT or finite-gap footprints.

Hunt--Crossley loss remains on the original separate dissipative port, which
checks the actual rate and step. A conservative path mean supplies no physical
duration and cannot waive those checks. Composites without attached JKR/Maugis
retain their existing stepping behavior. The `adhesive_mixed_storage` regressions exercise
all four owners together on three unequal bodies, independent onset work,
retained BEM source/authority, tiny and reversing paths, native Hertz/table
closure and reopening, original rate/step refusals and atomic failed-step retry.

The `adhesive_mixed_jacobian` regressions independently check the complete
mixed-source Newton action and advance native/prepared three-body motion just
inside the JKR fold with the original energy and momentum accounts.

### Shared source normal work on the complete continuing modal network

`contact::source_normal::SourceNormalModalSystem` binds one original
`BoundNormalContact` to the complete `CoupledModalSystem`. Its unchanged
physical work column spans all participating bodies and modes. Original modal
stiffness and damping, multiple bilateral springs and dashpots, rest loads,
material-memory states and their separate stores/losses remain with the native
owner. The same native prediction and held-force response supply a cached scalar
compliance; the original network stages and publishes the resulting body,
material and admitted radiation states.

The scalar solve consumes `NormalLawSource::mean_force_between` and its original
analytic endpoint derivative, also used by `BoundNormalStorage`. Conservative
Hertz, explicit attached JKR, finite-zone Maugis, finite-range DMT and retained
finite-gap tables share this path without a fitted common law. Maugis retains
its original selected contact/noncontact domain and bounded work/derivative
budgets. Signed attraction and negative
adhesive potential remain physical. Unilateral laws retain the inactive
negative-approach fraction through contact onset and release. Positive
Hunt--Crossley dissipation refuses; native network damping does not replace
its contact loss.

Safeguarded analytic Newton refinement retains the positive-slope mechanical
root connected to the accepted JKR state and respects explicit source-domain,
iteration, force and coordinate bounds. Maugis requires the original positive
load-line-slope bracket; the JKR stable-root search does not infer a Maugis
post-fold branch. A coarse ambiguous DMT step, source fold/domain escape or
unresolved root refuses. Contact uses the straight
approach path between accepted native endpoints, with a held work-conjugate
reaction; continuous within-sample turning points, attachment and snap release
are not resolved by this lane.

The adapter reevaluates the original source on restored native geometry and
checks the force correction actually represented in every mode. Contact
potential change, source path work and native contact work close independently
at the contact energy scale; unrelated network storage cannot inflate that
tolerance. Authored external work, every original native loss and total storage
then pass the original network and combined energy checks. The network's own
frame retains its complete applied-force port; the outer frame splits authored
and contact work once. `step_with_admission` stages receiving state before final
native publication and cancellation. Refusals preserve accepted coordinates,
material memory, clocks, ports and previous observations.

Focused evidence uses the original ordinary complete-pencil propagator: four
stiff modes and two Kelvin--Voigt links under tensile JKR; twenty steps of native
Maxwell memory, storage and independently integrated loss; source-domain,
receiving-state and cancellation rollback; and Hertz/table onset and release,
DMT gap attraction and explicit lossy-contact refusal. Finite-zone Maugis cases
at lambda 0.2 and 5 compare the complete native modal/Maxwell propagation,
original source work, material loss and total energy closure. The retained-table test
does not rerun its BEM producer. Prepared radiation/material-ingress paths and
full typed-port/workspace audits are not newly claimed by the bounded
verification projection.

#### Prescribed-temperature continuation of attached JKR on the full network

`SourceNormalModalSystem::apply_jkr_source_refresh[_with_admission]` accepts
the original immutable held-indentation `SourceJkrThermalRefresh` at actual
accepted native approach. It validates the complete original JKR source and
canonical predecessor/successor responses, attached phase, geometry, force and
energy limits, and positive original mechanical load-line slope. The same
physical work row and gap bind the refreshed original source. Native
coordinates, velocities, sample clock, modal/link coefficients, Maxwell memory,
all native stores and accumulated losses remain unchanged; the cached native
compliance stays valid.

The original signed parameter supply is checked against actual contact
potential change independently of the network's energy datum. It is not a
caloric transfer or a thermal update of the native bulk models. A receiving
candidate and final cancellation check precede publication. Source/history
mismatch, lost equilibrium, work/force/domain refusal or cancellation preserve
the previous source and all accepted native history.

Heating can move the original attached fold inside the caller's old trial
interval. Subsequent scalar solves intersect the unchanged authored approach
limits with the current original JKR domain; accepted geometry is never clipped
and no caller bound is widened. This permits a still-stable heated contact to
continue while refusing motion through its new pull-off boundary. Focused tests
retain native Maxwell history through changed source forces and signed work,
exercise an outward step using this shifted lower bracket, and check no-op,
foreign/stale-source, receiving-refusal and late-cancellation behavior.

#### Prescribed-temperature continuation of finite-zone Maugis adhesion

`SourceNormalModalSystem::apply_maugis_source_refresh[_with_admission]`
consumes the original sealed held-indentation Maugis receipt through the
retained mechanical potential. The original contact/noncontact policy, source
work budgets, physical row and gap survive. A selected deformed noncontact
disk keeps its true zero solid area and original attractive force; losing
its physical support refuses rather than downgrading the policy.

The successor is admitted at the unchanged actual native approach, with
original force limits and positive mechanical load-line slope. Actual contact
potential change independently closes against the signed parameter work.
Native modal/link coefficients, bodies, velocities, sample clock, Maxwell
memory, all native stores and physical losses remain unchanged. The explicit
temperature driver owns parameter work; it is not inserted into native heat.
A receiving candidate and final cancellation check precede source publication.

Future scalar trials intersect authored limits with the refreshed original
Maugis contact domain or explicitly selected noncontact fold. This only narrows
trial bounds and never clips accepted geometry or widens caller controls.
The JKR-specific root isolation still makes no Maugis snap assumption.
Focused contact/noncontact cases preserve complete native state bit for bit,
continue outward motion after a shifted fold, retain exact no-op evolution,
and exercise stale/foreign/load-control, receiving and cancellation refusal.

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
not the signed net adhesive reaction. All shared dry pressure contributions
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
three-component path mean and its optional derivative at fixed accepted history
and interval. The Newton matrix is the full `I + S A D(t_mean)`; damage and
surviving-bond viscosity enter through the same original path. Current and
candidate derivatives are retained with their material trial, so smooth paths
need no displaced constitutive probes. Unsupported nonsmooth cases retain the
existing bounded numerical fallback. The native wrapper stages the original mechanics
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

### Prescribed temperature on the original cohesive modal network

`CohesiveModalSystem::apply_temperature_successor` continues the accepted
source-bound interface at fixed local jump. An enabled residual or
damaged-area friction mode requires the paired
`apply_temperature_and_friction_successors` operation; the original two
sources retain their card, selected claims, surface order and held coordinates
and must agree on the complete new query. The material owner's existing
temperature transfer runs once. It preserves irreversible history, accounts
for signed coefficient work and additional fracture release separately, and
cannot heal or rebond on cooling. Friction occupancy is retained exactly;
coefficient changes supply no slip displacement or frictional heat.

The returned report records the unchanged sample and jump, old/new damage,
interface and total storage, signed parameter work, additional damage loss,
and an independent double-double energy residual. Its tolerance scales only
with the interface energy account, so unrelated native energy cannot mask
contact defects. The original material tolerance and total native energy
ceiling also remain in force. Existing total-force limits retain their
actual-mechanical-audit semantics.

Only the original source, material point and enabled local friction binding
are replaced after all audits. Bodies, every modal coordinate and native
material/Maxwell history, work columns, full compliance,
clock and last accepted mechanical frame are retained. The next ordinary step
uses the new local coefficients. No mechanical time, viscosity loss,
temperature trajectory, bulk expansion or heat destination is inferred.
`apply_temperature_successors_with_admission` lets a receiving owner stage its
account before publication; receiving refusal or late cancellation leaves the
whole accepted mechanical/interface state unchanged.

G1/G3 regressions in `tests/cohesive_modal/temperature.rs` cover original
moving bodies, independently computed coefficient work and damage release,
cooling without healing, same-source continuation, paired friction occupancy,
next-step mechanical response and atomic refusal/retry.

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
implicitly selects coast. `GasRigidRocketState<I>` retains either native
fixed-COM inertia model, with `RigidRocketInertia` as the compatible default.
`TensorRocketInertia` keeps independently oriented dry and propellant tensors
in their original material axes as mass is consumed; its environment carries
prepared `SpecificInertia` for angular outflow without a diagonal projection.
The native sealed inertia owner supplies exact zero exhaust for closed valves.
All gas/source, inventory, axial-impulse, cancellation and native mechanical
budgets are shared between the two models. Neither model infers slosh, shifting
COM, outlet geometry, internal-flow momentum or a complete world exhaust-energy
audit. New full-tensor coupling regressions exercise original-owner parity,
noncommuting angular outflow, world spin balance, exact-gas coast and atomic
refusal; native execution of those added cases remains pending.

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

The initial and predicted midpoint drag must each resolve the complete outer
interval. For locally frozen density and coefficient, quadratic drag has
Jacobian norm `2|F|/|u|`. The point mobility is
`G = I/m - [r]x I_body^-1 [r]x`; its largest eigenvalue is bounded by `1/m`
plus the trace of its positive rotational part, evaluated in the original
principal frame. The outer condition is `dt * 2|F|/|u| * mobility_bound <= 1`.
`AtmosphericFlightError::DragResolution` identifies the failed sample and its
maximum duration without publishing a candidate. Refining native attitude
substeps alone cannot resolve a drag force held over the outer interval.
This prevents the demonstrated coarse-step reversed-predictor acceleration;
it is a local frozen-coefficient resolution condition, not a bound on arbitrary
coefficient derivatives, evolving wind or other forcing, and not a global
passivity/accuracy certificate. The condition uses air-relative slip and allows
wind to do positive inertial-frame work on the vehicle.

This smooth composition is second order, with physical outer timestep and
native substep controls; it does not locate impact, cutoff, staging or source-
domain events. The caller splits intervals at known mode changes and refines
changing flow regimes. The original inertia model remains typed through
`RocketFlightPropulsion<E>`, `RocketFlightForcing<E>`, the midpoint predictor,
accepted `RocketFlightStep<I>` and `AtmosphericFlightCandidate<I>`. Full-tensor
vehicles retain moving principal axes, original material coordinates and native
angular outflow. Each environment sample recovers point angular velocity and
drag mobility from that sample's actual remaining full inertia. Default type
parameters preserve the common-principal-axis API and arithmetic. Fixed common
COM, ideal frozen nozzle and source-domain limits remain in force. General
internal-flow/transverse outlet linear momentum, slosh, aerodynamic heating and
a full exhaust/world-energy audit are not inferred. Added tensor-flight cases
cover powered/native-owner equivalence, point-flow recovery, closed-valve
continuation and late source refusal; their native execution remains pending.

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

### Rotating aerodynamic attachments in native momentum evolution

`GasRigidRocketState::{advance_spatial,coast_spatial}` carry the original
`PreparedSpatialLoads` through the complete finite-gas/native-mechanics
transaction. Atmospheric flight prepares each sampled aerodynamic force as a
world-held force at its declared body attachment. Every native torque kick
reevaluates the rotating arm's moment; the sample's reported body torque
remains an observation at that pose. Additional free body couples retain their
own frame, and the force resultant enters translation once.

Both powered flight and closed-valve coast use this path, including the original
full-tensor inertia model. Propellant, nozzle impulse and rotational exhaust
retain their existing owners. The independent outer drag-timescale admission
and atmosphere/source refusals remain in force; resolving attachment rotation
does not resolve an otherwise frozen drag or environment law. Focused tests in
`rocket_blowdown::flight::tests` compare spinning coast against an independent
planar RK4 moment reference and verify the powered path against the original
native propellant/spatial-load composition (G0/G1).

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

The linear thermal path admits scalar or oriented positive-definite conductivity
and induced-linear globally affine source-backed mu(T), fixed geometry and
prescribed velocity. Linear, nonlinear and complete-heat candidates retain the
original conductor type, including mixed scalar/orthotropic cells. Each thermal
trial refreshes all original principal conductivity/resistivity keys from the
same complete material bundle and preserves the declared physical frame. No
scalar surrogate, averaged conductivity or independently prescribed heat is used.
The native charge solve, conductor-frame E dot sigma E heat and signed Lorentz
work remain the only electrical-transfer owners.

The source-refresh helpers are shared with existing scalar planar consumers;
source identity, phase, volume and caloric admissions are unchanged. Tests in
`spatial_induction_enthalpy/orthotropic.rs` compare rotated, temperature-dependent
source substitutions with an independently eliminated Whitney-cube/finite-heat
balance, check each refreshed principal property, and continue the same histories.
These are synthetic equation checks, not empirical material data.

`SpatialMagneticThermalNetwork::step_nonlinear_candidate` connects the original
three-dimensional saturating/fixed-recoil A–phi owner to this same finite-cell
residual, source preparation and whole-volume local heat routing. Every endpoint
trial re-evaluates the original conductivity or resistivity claim at the actual
local caloric temperature and continues the accepted magnetic history. Magnetic
H(B), linear and recoil laws remain temperature independent over their complete
source ranges, so field internal energy equals the existing magnetic storage.
The physical account includes the original signed impressed and mechanical work
and the positive backward-Euler Bregman loss. Only E_conductor dot sigma E_conductor
heats the finite stores; numerical magnetic loss is never credited as heat.

Focused nonlinear regressions compare independent charge-eliminated cube
solutions with linear and saturating material-data substitutions, both original
electrical source bases, and continued field/thermal history. They separately
check native Joule heat, finite-cell energy, numerical loss and full balance.
A prescribed-motion reference also separates actual motional-current heat,
signed Lorentz work and an independent external thermal source.
Changed source identity, missing caloric receivers, lost conducting support,
trial-budget exhaustion and late cancellation must leave both accepted owners
unchanged. These are synthetic equation checks, not empirical material data.

Temperature-dependent nonlinear/recoil magnetism, hysteresis, thermal expansion,
moving interfaces and mechanical integration require compatible additional
physical owners; they are not inferred by either thermal entry point.

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

## Nonlinear prestressed-film material memory

`MembranePotential::with_proportional_relaxation` adds a supplied uniform
proportional plane-stress Maxwell spectrum to the original statically
condensed film stretching. Ratios multiply the complete equilibrium membrane
tensor, retaining directional and Poisson cross terms. The same scalar
spectrum applies to every facet. Zero or fully relaxed initial histories lie
in the original equilibrated strain range; its common branch rates preserve
that range. This is the restriction that makes the existing cold in-plane
condensation valid. Heterogeneous spectra and independently prescribed
in-plane histories require another owner and are not inferred here.

Original fs-la Cholesky factors each positive local energy metric, then its
Householder QR transforms the complete weighted strain field. Every triangular
row survives, including zero and weak rows. R^T R is independently checked
against the original physical Gram; no inverse, eigenvalue/rank cutoff, modal
loss fit or separate time integrator is introduced. Cached quadratic maps need
n(n+1)/2 histories per positive source arm, independent of the facet count.
Full source-arm, setup and joint-state budgets remain binding.

A film retains its law when cloned into an ImpactSystem. Its original fs-phs
mapped-strain owner evolves all material memories together with mechanics,
contact and compatible supported-string histories. The native discrete effort
owns positive arm heat, counted once in total dissipation. Installed tension,
linear bending, additional physical resistance ports and exterior radiation
remain distinct. Intrinsic per-mode damping on a hereditary film refuses.
No material-card claim, automatic heat-cell deposition, live thermal evolution,
prestress relaxation or heterogeneous in-plane evolution is invented.

Explicit modal-displacement limits bound all quadratic strains and gradients;
accepted histories obey the corresponding strain envelopes. The complete
nonlinear tangent bound includes both strain-gradient energy and hereditary
stress curvature. Source frequency, midpoint clock, dt/tau<=2, original slope,
energy and cancellation gates remain active. Subsequent radiation attachment
rechecks the complete source bound; nonlinear gap loading requires a new
strain/history envelope and refuses through the existing composition gate.

The bounded native proof used actual new source projection, cached maps and
body-memory composition, original strain energy/tangents, fs-la factors and
the full fs-phs discrete-gradient owner. Six tests passed, including full-rank
and singular mixed-pair fields, independent Hessian differences, domain/budget
refusal and independent Duffing-Maxwell state/heat convergence. At 100/200/400
steps, normalized errors were 1.06214e-4, 2.65565e-5, 6.63930e-6 for unrelaxed
history and 1.06077e-4, 2.65222e-5, 6.63074e-6 for relaxed history. These are
manufactured numerical checks. The existing `membrane_impact` suite now also
contains physical disk, mixed-pair, acoustic-composition and refusal
regressions; the full integration suite was not run in this dependency-limited
session. Mesh/material/acoustic convergence is not established by this slice.

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

## Simultaneous impacts of retained cohesive bodies

`cohesive_contact::impact::normal_impacts_candidate` connects the original
`CohesiveMotion` states to the existing simultaneous normal-contact owner.
Caller-owned signed sparse traces sample relative velocity and scatter their
transpose impulses. For the actual retained reduced masses, checked CG forms
every entry of `W = sum B_b M_b^-1 B_b^T`, including cross-contact coupling.
The bounded contact matrix is passed unchanged to
`fs-contact::multicontact::solve_normal_impacts`. No diagonal approximation,
matrix averaging, successive pairwise impact law or dense structural inverse
replaces those owners.

Each body receives its complete summed generalized impulse once through
`CohesiveMotion::impulse_candidate`. Original damage, source-bound viscosity,
friction/contact bindings, displacement, clock and step budget survive. Exact
zero-impulse bodies retain their accepted state without another mass solve.
All bodies must share physical time; supplying the same accepted body twice
refuses. All candidates are staged before any result is returned, so one
body's velocity-limit, numerical or cancellation refusal rejects the event.

Returned velocities are re-sampled through the original traces. With prescribed
relative velocity `c`, external work on the retained bodies is `-sum c*p`.
Separate gates check the native velocity prediction, the sum of absolute body
work defects, `sum Delta K - prescribed_work - actual_relative_work`, and
`sum Delta K - prescribed_work + native_impact_loss`. Compensated accounts
retain opposing body-work contributions. The contact owner's positive-work
allowance remains visible in its original receipt; it is not subtracted from
the physical dissipation residual or relabeled as heat.

Mass solves, native SPD/PSD contact work, traces and numerical workspace have
explicit bounds; the workspace admission includes the original factorization
packing scratch. Redundant contacts retain the native velocity-determinacy
witness and representative-only impulse semantics. This supplies no individual
pressure or heat authority for an impulse gauge.

Contact discovery, nonpenetrating gaps, frozen event geometry, common frames
and system-level restitution are supplied explicitly. Arbitrary reduced traces
establish generalized work/momentum, not Cartesian angular momentum without
geometric point maps. Tangential impact, event-time search, persistent contact
and automatic heat/damage allocation are outside this normal-event adapter.

### Coupled normal and tangential impacts

`cohesive_contact::impact::frictional::impacts_candidate` extends the retained
body event to one normal and two tangential coordinates per contact. Explicit
orthonormal normal-first traces, Coulomb coefficients, native restitution or
target-speed laws, and prescribed relative velocities define the event.
The full contact-major `3m x 3m` mobility retains normal/tangent and cross-contact
coupling. The original `fs-contact::frictional_multicontact` owner resolves
normal complementarity and final-velocity tangential maximum dissipation
together, preserving its original determinacy, work and refusal results.

Normal and Coulomb adapters share original checked-CG mass response,
compensated transpose impulse distribution, whole-body state staging and
actual-successor work audits. Every body receives one complete summed impulse
through its existing cohesive owner. Damage, displacement, source-bound
material/contact histories, physical time and step budgets retain that owner's
semantics. Prescribed work is `-sum dot(c,p)`; both contact tangents contribute.
The final physical dissipation residual retains native positive-work allowances
rather than using them to erase an unexplained energy gain.

Admission bounds up to 21 contacts and accounts for native mass histories,
dense range/factorization storage, packed GEMM scratch, and Newton/FGMRES
workspace using the actual configured restart. Mass receipts identify body,
contact and normal/tangent component. Representative impulses for a redundant
event retain their native gauge semantics and confer no individual heat or
load authority.

Exact original-source numerical regressions cover sliding, sticking,
off-diagonal inertia, moving support work, unequal-body recoil, shared-body
multiple impacts, zero friction, cancellation, budgets and physical audit
refusals. The existing normal regressions exercise the shared implementation.
Full cohesive-owner regressions additionally specify damaged-state retention,
subsequent unloading and all-body refusal/replay. Frozen event geometry remains
caller-owned; event-time search, persistent-contact stepping and allocation of
instantaneous loss to heat or fracture require their own physical owners.

### Cohesive trajectories with affine contact geometry

`cohesive_contact::trajectory::advance_candidate` advances retained cohesive
bodies through a complete bounded interval with fixed normal-first contact
frames, sparse traces and affine prescribed surface motion. Constant reduced
loads and affine prescribed interface offsets enter each original native
midpoint/material-path solve. Native offset accessors reconstruct accepted
offsets through the original material map; target-time advancement gives every
body the same exact represented clock with its actual positive solve duration.

Contact gaps are reconstructed from initial gaps, actual displacement
increments and prescribed motion using compensated sums. Every accepted
single native displacement path is straight, so nonnegative start and endpoint
gaps admit the entire declared affine contact path. A penetrating proposal
shortens the duration and recomputes all native candidates from the unchanged
predecessor. Rejected material histories, interpolated physical states and
projected gaps never enter the accepted trajectory. All active contacts enter
the existing simultaneous Coulomb impact owner at the actual common event time.

Native residual/Newton work is reserved before each body trial; successful
trials refund unused work, and refused trials consume their reservation.
Accepted-segment, body-trial, refinement, event, workspace and minimum-progress
bounds apply. Whole-window publication is atomic. The report retains original
contact determinacy, event times and losses, physical storage/work channels,
and charged computation. One compensated account retains all accepted physical
increments; sums of absolute native smooth and impact defects prevent numerical
errors from cancelling into a claimed balance.

Positive activation distance is an explicit spatial approximation. Timestep
and activation-distance refinement govern trajectory accuracy. This admits the
declared discrete affine paths, not the exact continuous nonlinear flow,
undeclared geometric features or changing contact frames. Loaded resting
contact, unresolved chatter or exhausted time progress returns a refusal;
this adapter does not invent sustained reactions or claim a general
Moreau--Jean method. Instantaneous losses remain separate from material
fracture, viscosity and internal friction, with no automatic heat partition.

### Sustained contact on original cohesive midpoint evolution

`cohesive_contact::trajectory::sustained::step_candidate` supplies the finite
loaded-contact step that the transient impact trajectory deliberately refuses.
It solves simultaneous normal reactions and Coulomb forces for already active,
nonpenetrating, normal-stationary contacts. Fixed sparse traces in a declared
orthonormal frame remain power conjugate: each trial scatters `B^T lambda` into
the original reduced mean forces, then advances every original `CohesiveMotion`
to one common endpoint with its native midpoint/material-path owner. Complete
inertia, stiffness, native bonded damage, viscosity and internal friction enter
the response. No frozen mass-only mobility stands in for finite-step mechanics.

The shared `fs-contact::sustained` force law and original Newton/FGMRES solve
couple all declared contacts through that actual response. Numerical N/m
scales and force perturbations are explicit conditioning controls; forces and
physical states are never projected. Warm starts reuse force guesses only.
Every response starts from unchanged accepted bodies; obsolete trial histories
are released, and the last successful solve callback is a fresh native replay.
Force/path bit identity is checked before the replay's actual original body
steps become the candidate. Native body calls, residuals, Newton work, trace
storage, retained histories and force applicability have explicit budgets.

The independent audit checks actual nonnegative gaps and normal forces,
midpoint trace kinematics, loaded endpoint stationarity or separating release,
normal complementarity, the Coulomb disk and actual sliding direction. A
positive-capacity sliding reversal or false zero-slip reflection refuses and
requires interval refinement. The declared straight displacement path and
nonnegative endpoints admit the fixed affine gap throughout this numerical
interval; they do not certify exact continuous nonlinear motion.

With relative trace increment `s=B Delta u+h*c`, contact body work is
`lambda dot s-h*c dot lambda`, prescribed-surface work is `-h*c dot lambda`,
and dry-friction loss is `-lambda_t dot s_t`. Original authored force work,
prescribed jump work, kinetic/bulk/cohesive storage, fracture, viscosity and
internal friction are retained separately. The full first law and absolute
native/contact defects must pass. Signed and absolute normal constraint work
remain visible numerical defects; they are never renamed friction or heat.
No partial native state publishes on a refusal or cancellation.

This is local continuation on the original material branch for one bounded
interval, with fixed geometry and prescribed affine surface motion. Incoming
impacts and new-contact discovery retain their original owners. It does not
claim global frictional uniqueness, general Painleve resolution, automatic
stop-event discovery, a general Moreau--Jean method, or a thermal partition.
Generalized momentum/work follows the supplied traces; Cartesian force/torque
conservation additionally requires the caller's actual geometric point maps.
Focused core and audit tests cover cold support/stick/slide/breakaway solves,
full coupled inertia, elastic storage, shared contacts, moving surfaces and
physical refusals. Real damaged-cohesive-body integration cases are retained
for the complete native runtime build.

### Native sliding stops on the original contact response

`cohesive_contact::trajectory::sustained::segment_candidate` advances supported
motion to the requested endpoint or an admitted relative sliding stop. It
reuses the existing complete held-force solve and original cohesive midpoint
owner. An unpublished reflected trial supplies a signed time bracket; every
subsequent time query recomputes contact forces and native material evolution
from unchanged accepted bodies. Prescribed jump offsets are interpolated to
each actual represented endpoint from their original physical time path.

The allocation-free stopping locator alternates safeguarded secant and
bisection. An approximate stop requires both a narrow actual time bracket and
a small complete relative tangent-velocity vector. A perpendicular turn is not
a stop. An exact full-vector zero may identify a point bracket. The final query
freshly replays the original response; actual duration and velocity identity
are checked, then the unchanged all-contact and full-work audit must pass.
Another observed contact reversal triggers its earlier bracket before adoption.
Simultaneous stops are listed only from actual relative velocities; a moving
surface or second body can leave nonzero world velocities. No velocity, gap,
force or material state is projected.

Every failed or discarded probe spends bounded contact/native work but adds
no physical displacement, loss or irreversible history. Final body receipts
retain the work of the one actual admitted segment. Stopping brackets describe
the local sampled native endpoint response, not earliest exact continuous-flow
events or globally unique frictional dynamics. Minimum-time, transverse-turn,
nonrepresentability, source, work-budget and cancellation refusals stay explicit.
The original fixed-duration sustained API retains its strict behavior.

Direct locator tests include constant/viscous response, full coupled inertia,
moving surfaces, two-body momentum and false midpoint reflections. Native
damaged-body cases exercise stopping followed by rest or opposite breakaway,
moving-wall entrainment and coupled pair momentum through the original owners.

### Mixed free, impact and supported cohesive trajectories

`cohesive_contact::trajectory::evolution::advance_candidate` advances a whole
requested window through the existing transient loop and original physical
owners. Incoming contacts enter the simultaneous native Coulomb impact owner.
Already active, normal-stationary contacts enter the stop-aware sustained
segment; separating or open contacts resume native free motion. A completed
relative sliding stop ends a segment, and the next actual force solve determines
sticking or opposite breakaway under the continuing authored loads. The
existing transient-only API retains its original route.

All declared contacts remain in the actual endpoint geometry checks, including
those omitted from the current sustained force solve. A penetrating proposal
shortens the duration and recomputes mechanics from unchanged predecessors.
Original interface offsets follow their affine physical time path across both
ordinary substeps and shortened stopping intervals. Warm forces are numerical
guesses in stable original contact order; an impact clears them.

The wrapped report separates sustained Coulomb loss from native material
friction and instantaneous impact loss. Independently evaluated authored work
replaces native total force work in the shared window account, then prescribed
surface work and dry-friction loss enter exactly once. Original storage and
material losses retain their meanings. Sums of absolute native, contact,
friction, normal-constraint and segment energy defects must also pass the
whole-window gate. Geometrically rejected successful probes spend work but add
no motion or loss. Native body/residual/Newton, contact solve/response, retained
event and workspace budgets apply across the entire window; publication is
atomic for all original histories.

The geometry and numerical-path limits of both existing owners remain. A
positive activation-distance impact is not silently relabeled as zero-gap
support: its unchanged normal-kinematic/work gates can refuse that handoff.
Transverse directional crossings are not invented stopping events. This does
not certify exact-flow event discovery, arbitrary frictional uniqueness or a
general Moreau--Jean scheme. Native regression cases cover sliding to sustained
rest, stopping then genuine opposite breakaway, zero-gap plastic impact followed
by loaded support, and contact release followed by free flight.

### Normal arrival before cohesive impact and support

Mixed trajectories now use `EvolutionControls::arrival` to distinguish a
declared coarse activation band from the requested physical gap accuracy.
Impact selection, all metadata passed to the simultaneous impact solve, and
sustained-force selection use the smaller of those two distances. An initially
open body inside the coarse band therefore retains its original free mechanics
until the refined band is reached. A zero arrival gap tolerance requests an
actual numerical zero; it never authorizes changing a position or gap.

When an original free or supported proposal crosses an open contact, the
arrival locator reevaluates the original response from the same accepted
predecessors. Scaled secant proposals alternate with bisection; neither an
interpolated state nor a rejected material history becomes a successor.
Approximate acceptance requires both a time bracket within the requested bound
and an actual nonnegative lower-endpoint gap within the spatial bound. An
actual zero gap can give an exact numerical point bracket. Every result is
freshly replayed at its represented native clock, with unchanged gap and normal
velocity, before the retained original body candidate can be adopted.

Every declared contact still participates in geometry admission. If another
contact penetrates at the proposed arrival, its earlier bracket is resolved
from the same original predecessors. An earlier actual sliding stop is retained
as its own endpoint before normal arrival continues. Arrivals can therefore
occur during an existing sustained segment, while simultaneous impulse
coupling and the subsequent force solve remain with their original owners.
Starting at zero gap with outward motion still needs a positive-gap predecessor;
the existing bounded shortening and explicit no-progress refusals apply there.

`EvolutionReport::arrival_events` retains the selected contact, actual native
time, relative bracket, outward-rounded absolute bracket, actual gap, relative
normal velocity, and response evaluation count for each admitted arrival.
Retained receipts are bounded by the accepted-step limit. All discarded native
probes spend the shared body/residual/Newton budget; supported probes also spend
the shared contact-solve and mechanical-response budgets. They contribute no
physical work, loss, history, or warm-force acceptance. Original source failures,
cancellation, clock-resolution limits and exhausted budgets remain refusals.

The unchanged normal-kinematic, momentum and work audits remain decisive after
arrival and impact. A small positive gap is still positive and does not
guarantee a support handoff at every later step size. Brackets describe the
sampled original numerical response, not the earliest collision of an exact
continuous trajectory or arbitrary changing contact geometry. Native cases
exercise accelerated arrival, a translating wall with its actual source work,
an arrival during existing support, activation-band refinement, and atomic
budget/cancellation refusal.

### Angular resolution of sustained sliding

`SegmentControls::turning` bounds the actual change in a loaded contact's
relative tangential velocity direction. Its minimum duration, angle strictly
below pi/2, and refinement count are numerical resolution controls; they do not
change the Coulomb coefficient or declare a new physical stopping event.
Normalized full tangent vectors and a compensated determinant distinguish
resolved transverse motion from a collinear reversal without overflowing large
velocity products or manufacturing a transverse remainder.

A resolved turn exceeding the angle bound shortens the original interval and
re-solves all held forces and native mechanics from the same accepted bodies.
Every contact contributes to the final angular bound, and the complete original
gap, momentum, force and work audit remains mandatory. Original prescribed
offsets follow their actual shorter physical time. No endpoint velocity or
history is projected. Collinear reversal continues through the original full
relative-velocity stopping locator. Its intermediate native probes are also
checked for transverse turning, so an apparently collinear long interval can
reveal its change of direction before the scalar stop search exhausts.

The segment's optional `TurningReport` records the trigger contact, actual
duration, maximum relevant endpoint angle and refinement work.
`EvolutionReport::turning_segments` retains these bounded numerical receipts
separately from physical stop and impact reports. A segment can carry both a
turning receipt and a subsequent genuine stop at the same native endpoint;
that is one accepted segment. Arrival searches preserve an actual earlier
angular substep before continuing normal contact localization.

All force solves, mechanical responses, native residuals and Newton work from
discarded angular trials remain charged to the same segment/window budgets.
Minimum duration, unrepresentable shortening, source failure and exhausted
refinement/work budgets refuse without partial publication. The angle bound
describes the original discrete response endpoints, not a global integration
error estimate or an exact continuous-flow event-time theorem.

## Finite cohesive patch load line

`cohesive_contact::patch::CohesivePatchLoadLine` solves all opening and two shear coordinates of up to 64 original bonded quadrature sites together. It consumes the full caller mechanical cross-compliance, heterogeneous original `CohesivePoint` histories and reference areas. The original analytic mean traction supplies physical work. Finite lever arms transmit moments; local failure redistributes force through the retained off-diagonal response. No average damage or independent-site solution replaces the patch.

Prepared scratch and LU storage allocate only at construction. Each current/candidate material evaluation also retains its optional original path-mean tangent. The complete cross-site matrix multiplies each local three-by-three derivative by that site's reference area and every original compliance row. Only a site without an available tangent retains local finite-difference columns; those probes execute the original mean-only path. Smooth heterogeneous elastic/viscous fixtures require two material calls per site for one Newton correction, including the initial and candidate trials. This is bounded evaluation-count coverage, not a measured wall-clock speedup. Explicit total point-law evaluation, iteration, backtracking, per-site force/jump and compatibility budgets apply. Refusal/cancellation mutates scratch only; callers must independently audit actual mechanical motion before publishing all histories. Softening continuation does not certify root uniqueness, spatial resolution, finite-rotation objectivity or timestep accuracy. Tests cover heterogeneous full cross-compliance, force redistribution after local failure, exact finite-body force/torque/work balance with progressive edge peeling, the one-site limit, and atomic refusal.

## Native finite cohesive patches

`contact::cohesive::patch::CohesivePatchModalSystem` attaches heterogeneous `CohesiveModalSite` quadrature to the original `CoupledModalSystem`. Each site retains a borrowed material law or the original ordered source, reference area, and three complete modal work columns. Setup computes every cross-site/cross-axis compliance term through the original held-force response under an explicit projection-work budget. Translational, rotational and elastic participation therefore retain their reciprocal forces and moments. Construction explicitly bonds the current sample-zero configuration; it does not infer bonding from proximity.

The prepared patch solve stages original native mechanics with the common held force field. It independently reevaluates every material point on actual body motion and gates each jump and force component, then the complete native-plus-bond energy window. Material calls reserve one final audit per site. Only after all checks and the final cancellation poll do all body states, irreversible histories and per-site observations publish together. Original fracture, surviving-bond viscosity, native losses, external work and prescribed jump work remain separate. Contact scratch allocates at construction; the original propagator retains its existing allocation rules.

The explicit `patch::temperature` continuation accepts one original-source successor per site, including heterogeneous temperatures and quadrature areas. Each successor retains its own card, selected claims, ordered surfaces and non-temperature coordinates. Coefficients change at fixed accepted jumps before the original material owner admits irreversible damage release. Reports retain signed parameter work and new fracture release separately for every site and for the whole interface. Independent per-site gates precede compensated batch storage/work/release accounting; its tolerance scales only with interface accounts, so unrelated body energy or cancellation between sites cannot conceal a local defect.

Enabled damaged-area friction requires every paired successor at the exact corresponding query. The shared local continuation retains each canonical point, original site identity and branch occupancy without inventing slip velocity or frictional work. Its original per-site work gate and the tighter native/friction material tolerance apply. Only source-dependent local laws and cone scales change; full cross-compliance, geometry-based area scaling, prepared dense algebra and all native mechanical/Maxwell state remain in place. Candidate storage is linear in site count and lies outside mechanical stepping. Receiver admission and the final cancellation poll precede every publication, including source references and friction state. The native clock, work columns and last mechanical observations remain unchanged.

The native integration regressions cover quadrature-area partition, finite lever-arm force/torque and edge peeling, and multi-site refusal/replay. Temperature fixtures additionally compare heterogeneous-area coefficient work and damage release with independent common-mode expressions, heating/cooling without healing, natural load redistribution through complete off-diagonal compliance, continued friction with refreshed coefficients, and exact retry after source/receiver refusal or late cancellation. The full native executable requires the repository runtime build. Fixed frames, linearized geometry, declared spatial quadrature and straight local jump paths remain explicit; no finite-rotation contact search, remeshing/history transfer, spontaneous rebonding or automatic fracture heat partition is inferred.

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

## P1 reference inventory in full tensor-memory membranes

`TensorMembrane::new_with_reference_areal_density` retains three positive
reference areal densities per original facet and the existing `MembraneMass`
choice. Density is linear within each facet; incident values at a shared node
may differ. These are immutable kg/m² inventories on the original positive
facet area, not current-geometry samples. Thickness is not applied twice.
The explicit field replaces only scalar inertial density; authored tensor laws,
thickness, material frames and complete Maxwell memory retain their owners.
This entry does not infer constitutive properties or measured source authority.

The checked public `fs_solid::membrane::inertia` primitives
`p1_reference_row_masses` and `p1_reference_consistent_mass` are the same kernels
used by `MembranePatch`. They integrate the linear density and both native P1
velocity shapes, including all cubic kinetic moments. Tensor membrane assembly
restricts both consistent-mass axes by the original support map. Its full-mesh
row sums drive the existing uniform-acceleration loads without donating fixed
inventory to free vertices. Original scalar constructor arithmetic is retained.

The existing tensor-material inertia and exact propagator still own canonical
momentum, force inverse, physical kinetic storage, hereditary stress, source
work and dashpot heat. Assigning an initial inventory adds no heat or parameter
work. No live mass change, constitutive resampling, finite deformation, modal
truncation or acoustic loading is introduced. Original full-state and material
work caps remain; malformed inventory or cancellation returns no partial body.

Focused regression sources cover scalar parity, independent cubic velocity
quadrature, winding/corner permutation, equal-total-mass oscillators, native
hereditary work and heat, discontinuous density under native acceleration loads,
and large-contrast/invalid inventory admission. These are manufactured numerical
checks and do not assert experimental material validation.

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

`FrameThermoelasticSystem::step_heat_assembly` now accepts the original
`HeatAssembly` contract on this same frame/circuit/enthalpy residual. Finite
coolant links, prescribed reservoirs, thermal interfaces, integrated conduction
and enclosure radiation retain their complete associated source sample in
`FrameThermoelasticReport<S>::source`; concrete source errors survive refusal.
The existing prepared-network entry retains its original report shape. Circuit
Joule heat, reversible material-field heat and reversible frame heat remain
with their existing owners, and internal radiation remains separate from
external supply. Both cancellation authorities admit before joint publication.

The explicit `*_with_link_conductances` cold constructors admit a disjoint
partition of native thermal links alongside the original Fourier bindings.
They retain the original cell/link slices and do not cache external transport
coefficients. Ordinary stepping continues to refuse incomplete coverage; each
complete-assembly trial resolves the actual source laws afresh. The same
admission reaches standalone/shared frames, original contact and material-field
circuits, and the existing full-inertia fluid-loaded frame wrappers. Exterior
mechanical-fluid losses retain their separate account. Reference geometry,
source support, solid phase, capacity/strain bounds and first-order splitting
retain their existing limits.

The focused native fixture in
`fs-airflow/tests/material_cell_convection/film_source/heat_network/finite_links/frame.rs`
defines powered and unpowered finite-coolant cases, independent two-cell/RC
balances, source-only coolant substitutions, subsequent exact forced-frame
motion, complete source receipts, in-residual source-domain rollback and
separate-authority cancellation. The changed Rust files pass the pinned syntax
parser; native test execution remains unverified in the shared dependency and
resource environment.

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

## Native multi-electrode circuit storage on the spatial enthalpy timeline

`circuit_enthalpy::conduction::spatial::circuit::SpatialConductionCircuit`
binds original electrode face sets one-to-one to native circuit nodes and retains
branch-ordered capacitor charge and inductor flux.
`SpatialConductionThermalNetwork::step_circuit_candidate` refreshes complete
scalar/tensor and two-sided electrical-contact sources at each receiving
temperature/phase proposal, extracts the full native terminal matrix, and
advances the existing finite C/L descriptor from the same accepted coordinates.
Actual midpoint terminal voltages recover combined E/J and bulk/contact Joule
power, deposited once through the original caloric routes and explicit interface
partitions. Electrical storage, thermal history and source field publish jointly.

Independent gates compare circuit currents with electrode reactions, compensated
terminal sum I_i(V_i-V_ref) with the positive conductance factor, recovered local
heat and native dissipation, fixed-C/L coordinate energy increments with
descriptor storage, and authored source-port work with native supply. Each
residual evaluation charges N-1 field basis solves plus one native circuit step
before extraction. Original source support, volume/phase/contact ownership,
cancellation and joint acceptance remain in the existing owners. The shared
driver seam preserves the prescribed/moving consumer's accepted and trial
mechanical contact-pressure inputs.

The circuit consumer admits fixed geometry, fixed authored C/L, source ports
and ideal transformers. `new_with_couplings` additionally accepts fixed signed
mutual-inductance links by original graph branch index. Their full positive-
definite magnetic storage is preserved in every native field/thermal trial;
independent storage work uses mean current `L^-1 mean(flux)` times flux change,
including all cross terms. The same E/J solution supplies spatial heating, so
finite magnetizing/leakage dynamics transfer power into the original thermal
body without a separate transformer or heat approximation. Links and total
flux linkages persist through continuation. Magnetic parameter changes still
require an explicit reversible-work owner. The spatial circuit regression
independently reconstructs a loaded transformer, local heat, complete stored
energy and source work over three accepted steps and both winding orientations.

`new_with_material_resistors` and
`step_material_circuit_candidate` additionally admit ordinary graph resistors
with complete original material sources and explicit whole-cell caloric receivers
disjoint from spatial FE receiving volumes. Every nonlinear endpoint refreshes
the original resistivity/conductivity bundle at the receiver's actual T/phase,
retaining source coordinates and authored geometry. Accepted sources must agree
with the accepted caloric state before advancement.

The original resistor lowering and native midpoint branch-power accountant are
shared with `MaterialCircuit`. Spatial bulk, electrical-contact and lumped
branch heat remain separately reported; each enters the original caloric residual
once. Their sum independently closes native dissipation and fixed-C/L storage
plus authored source work. Missing or changed sources, duplicate or missing
branch routes, overlapping field/resistor volumes, unsupported domains, exhausted
shared solve budgets and cancellation refuse without publishing either history.
Source-dependent storage still requires an explicit reversible-work owner.
Native midpoint electrical work enters the existing backward-Euler thermal
residual; independent spatial and time refinement remain necessary.

Two additional integration regressions cover signed three-electrode RLC plus
temperature-dependent resistor feedback and capacitor discharge through spatial
bulk, two-sided contacts and a conductivity-based resistor retaining its original
lot coordinate. Independent scalar energy oracles, formatting and source/API
reviews passed. Native compilation and integration execution remain unverified
because the shared dependency build exceeded workspace memory/storage capacity.

Synthetic regressions cover signed two-interval three-electrode RLC/source
feedback and capacitor discharge through temperature-dependent two-sided
contacts, including source-support, work-budget and cancellation refusal.

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

The planar owner's `apply_source_refresh` accepts the original bridge's
immutable material-refresh candidate only against its exact current
predecessor. It retains native mode states, all connections, time/sample
count, mechanical gap and work column, rebuilding only the sourced load line.
An independent combined-energy limit and final cancellation check precede
publication. Returned signed surface-free-energy parameter work must be
retained by the caller separately from mechanical work and heat; the latest
rendered frame remains its historical sample. A native moving-body regression
checks source-force feedback without reset, complete energy accounting with
that parameter term, stale-candidate and cancellation refusal, and unchanged
body/observation history on both refusal and successful source adoption.

The axisymmetric owner's `apply_source_refresh` now accepts the original 3-D
Young--Laplace material endpoint through that same continuing network. Its
source can conserve closed-liquid mass while density and volume change. Exact
predecessor matching, the existing mechanical gap tolerance, native energy limit
and the original source's quadrature/root-work error estimate precede publication.
The true mechanical gap is retained independently of the newly converged profile;
neither numerical root residuals nor a material update move or reset the bodies.
Clock, coordinates, connections, work column and prior observation remain intact.
Returned surface parameter work stays separate from mechanical work and heat;
bulk thermal energy, evaporation and phase transitions require their own owners.
A native regression couples temperature-driven density/tension changes to actual
subsequent force and motion while checking mass, energy, stale-candidate refusal,
cancellation and preservation of all accepted mechanical/observation history.

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

The explicit `transfer_temperature` operation retains the accepted frame,
reference offset, actual world separation and original bonded area bit-for-bit.
It admits only a temperature successor of the same ordered source and model
family, then invokes the original common- or mixed-mode material history
transfer. Signed parameter work and fixed-separation damage release are scaled
by the original area and independently balanced against extensive storage.
The thermal stage has no mechanical displacement work or viscous rate.

`trial_with_temperature` composes that stage with the original three-component
finite-rotation physical trial. The original force pair and balancing frame
couple still supply spatial mechanical work. Its complete old-to-new balance
includes coefficient work, both thermal and mechanical fracture increments,
and the original timed viscosity; it never treats their numerical residual as
heat. Source, frame, separation and irreversible history remain unchanged if
either stage refuses. This is a prescribed-temperature staggered update;
continuous temperature paths, thermal transport, cure and fracture-heat
partition remain outside it. Source integration regressions exercise held
geometry, anisotropic work, finite-rotation covariance, rupture/reclosure and
source/model/late-mechanical refusal. Their full runtime execution remains
blocked in this worker by the pinned executor dependency's memory requirement.

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

The native `ResolvedThermoviscousDgInput::from_bytes_with_sources` forwards
the same original thermochemical catalogue and solid-material callback into
its embedded resolved source. Its cold composition moves the actual finite
supply inventory, NSF zone, primary DG load and unadvanced exterior receivers
into the complete field. The inventory is never reconstructed from pressure
and temperature; duct/mouth receivers retain their original geometry, signed
flow, material propagation kernels and output clock. Advanced histories,
upstream acoustic loads and reciprocal-mouth attachments still refuse because
their acoustic port basis needs an explicit complete-field binding.
`with_physical_renderer` exposes the original live participants and balances
within the same borrowed-field callback. `music_render` uses its existing
restricted attachment loader and emits both the complete-field and original
embedded-source information as valid JSON.

`resolved-thermoviscous-breath.performance` composes the existing manufactured
field geometry with a finite supply, volume/heat actuation and an observed
escape aperture. Native regression sources cover inventory/energy retention,
the independent original receiver-kernel sum, material caloric substitutions,
callback identity, source refusals and cancellation. These physics regressions
remain unexecuted here; only the extracted original JSON formatting and strict
JSON-reader check executed. The inserted thermoviscous field remains linear
about its reference gas; it does not transport mean breath mass or enthalpy.
Finite-source coupling retains the original held-state first-order partition
and still requires physical-clock and observation-clock refinement.

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

`AdhesiveDmtModalSystem::apply_source_refresh` continues the same native owner
at prescribed separate bulk-contact and interface temperatures. The exact
accepted indentation authenticates the original normal refresh; the original
selected shear claims follow its interface temperature. All native positions,
velocities, work columns, full cross-axis compliance, clock and last historical
frame remain unchanged. The new normal potential changes total native energy
by independently reported signed parameter work, with the original energy
limit and residual gate applied before publication. The next mechanical step
uses the new force and actual-area static/kinetic capacities. No heat is made
by changing junction strength at fixed coordinates. This operation changes
local contact parameters only; modal bulk stiffness, thermal transport,
expansion and surface entropy need their own physical owners. Two focused
native temperature groups exercise contact/gap/separated moving bodies,
independent potential and friction work, heating/cooling continuation, pinned
claims, source/control/cancellation refusal and exact retry.

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
autonomous source evolution, evolving frames, heat recipients and exact stopping events
remain outside this owner; a full Cargo/workspace check is separate.



### Fixed-jump temperature continuation of cohesive friction

`ReclosureLoadLine::transfer_temperature` and
`DamageFrictionLoadLine::transfer_temperature` stage the same original cohesive
point and residual Coulomb history at paired source-temperature successors.
Each original source must retain its card, selected claims, surface ordering,
non-temperature coordinates and constitutive bounds; the two new sources must
resolve the identical query. Area, full reciprocal compliance, numerical
controls and material/work allowances remain those of the original owner.
Admission is repeated at the new coefficients, including the reclosure solver's
normal/friction feedback screen.

The material temperature operation runs once at the held accepted jump. Signed
coefficient work, storage change and nonnegative additional fracture are
returned separately and checked with the original area-integrated work gate.
The operation preserves irreversible fracture history and the exact friction
occupancy, including an occupied kinetic branch; cooling cannot heal and
coefficient changes cannot rebond. Rigid Coulomb friction has no stored slip
energy, so its coefficient update supplies no friction work or heat. No
velocity or static reaction is inferred. Actual mechanical force limits remain
enforced by the original mechanical audit. The immutable candidate owner and
rebound state permit atomic publication by the receiving caller; every refusal
leaves the original owner and state usable. This operation supplies neither a
temperature trajectory, thermal expansion, heat transport nor heat partition.

G1/G3 tests in `tests/cohesive_friction_temperature.rs` cover independently
computed hot sliding forces with full cross-axis compliance, retained kinetic
occupancy versus static readoption, temperature-induced first damage without
repeated release, cooling without healing, exact same-query transfer, source
and state refusals, and the original reclosure feedback-admission bound.

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

## Regional reed pressure with independent mouth propagation

The physical `ReedPerformance` observer now binds an independent atmospheric
filter to every open or vented regional tone-hole mouth. It uses the original
regional site's effective radius, exterior gas and translated receiver position
`[receiver.x + L - x_hole, receiver.y, receiver.z]`. The bell and each mouth
therefore retain their own propagation range, attenuation, dispersive FIR delay
admission and omission policy. Each original Rayleigh receiver consumes that
mouth's filter delay allocation; pressures are filtered separately before
summation. No extra delay or pressure gain is applied after the sum.

The regional outlet owner accepts one positive finite propagation speed per
open/vented site in original bore order and validates the complete count before
constructing receivers. Closed sites have no receiver. The same original joint
transfer sweep, exterior density and loading, axial/chimney propagation, shared
replay gate and latency, numerical shunt port resistance and source priming
remain authoritative. Existing common-speed entry points retain their original
arithmetic. Source mechanics and reed/bore feedback do not consume atmospheric
observation filters.

`ReedPerformance::open_hole_observers` exposes each actual receiver and filter
receipt; `molecular_open_hole_observers` remains a compatibility alias. Molecular
and regional playback share the existing per-mouth filtering owner. The
lossless law preserves the original receiver arithmetic and pressure addition
order. The parallel-baffle geometry, independent mouth-load approximation and
existing finite-history limitations remain; no mutual exterior scattering or
new radiation power model is introduced.

Authored native regressions cover independent open/vented geometry and filter
receipts, per-mouth Rayleigh/FIR delay accounting, samplewise reconstruction
from independently filtered pressures, unchanged reed/bore mechanics, block
replay, common-speed/port-scale equivalence and invalid speed/count refusal.

## Native normal JKR attachment, capture and rebound

`contact::adhesive_normal_events::JkrNormalEventModalSystem` couples the
original source-owned JKR normal trajectory to original free modal
coordinates. A complete mass-normalized work column `B` defines the actual
relative approach and effective mass `1/(B dot B)`. Initial source history
anchors approach at native positions; subsequent motion uses position
differences, preserving common translation and uncontacted modes. Source,
ordered material identities, fixed normal and numerical budgets remain bound.

The free constructor retains its original scope.
`new_with_elastic_connection` admits exactly one original positive undamped
spring. `new_with_viscoelastic_connection` additionally admits that link's
finite nonnegative viscous coefficient. Both use the original column,
stiffness and rest extension without moving accepted bodies. This is one
Kelvin--Voigt mechanical connection, not an inferred rate-dependent adhesive
material law. Arbitrary elastic networks, component damping, memory, radiating
coordinates and harmonic-only preparation remain outside this adapter.

The source scalar discrete gradient resolves elastic load and passive
midpoint viscous work together with adhesion. Separate contact, spring and
dashpot impulses and first time moments restore the original bodies through
the existing event restoration owner. Branch changes preserve velocity and
release contact storage; spring storage is continuous and support dissipation
remains a distinct finite-time account. Original body limits and the combined
connection force bound over the complete trajectory still apply.

Native geometry independently supplies spring extension, storage, work and
the actual endpoint spring-plus-dashpot reaction. A separate loss-sized gate
reconstructs `sum(c * duration * midpoint_velocity^2)` from the original
coefficient and sealed source segments. The native connection retains its
spring storage, positive dissipation and internal port work
`-delta_spring_energy - support_loss`. Adding that storage change and loss
to the total reaction work isolates contact boundary work. Contact work,
coordinate/speed compatibility and whole mechanical energy retain separate
checks. No restitution, extra snap loss or heat fraction is invented.

Accepted bodies, native clock, source history and observations publish
together. `step_with_admission` stages a receiving candidate before the final
shared cancellation/publication boundary. Frames are cloneable and retain a
private optional shared immutable dissipation trace only for positive damping;
the original scalar step remains copyable. Full source/mass/support/policy,
time and step continuation prevents replay even after a closed contact cycle.

### Joint native mechanics and event-timed thermal transport

`adhesive_normal_events::thermal::JkrThermalModalSystem` owns the original
mechanics, source-bound heat history and native thermal timeline. Read-only
participant access prevents independent advancement. The original constructor
accepts free or undamped supported mechanics. Damped mechanics require
`new_with_support_heat` with an explicit separate `JkrViscousHeatRoute`;
omitting the support destination refuses construction.

Each interval transports snap heat at the original branch-event timestamps
and, when damping is present, adds the declared smooth support shares through
native source powers over sealed mechanical segments. Coincident scheduled
thermal events retain native ordering. Both channels have distinct work
ledgers, outside destinations and independent mechanical-to-thermal transfer
checks. A positive support loss never becomes adhesive snap heat.

The shared gate composes native relative-energy residuals and both transfer
defects with double-double accumulation. The frame's
`external_transport_heat_j` excludes both internal transfers; the underlying
`heat.transport.external_heat_j` still reports the original native source
account, including any support-generated heat. Boundary/source heat, incoming
material energy and both declared outside destinations remain explicit.
Large absolute caloric reference energies cannot replace small-work gates.
Any caloric, transport, mechanical or cancellation refusal preserves bodies,
spring/contact storage, both clocks, heat history and accepted observations.

Normal coefficients and both heat fractions remain frozen. This one-way
coupling does not infer temperature-dependent adhesion, contact conducting
area, topology changes, rotating normals or contact waves. Native explicit
thermal discretization, mechanical midpoint refinement and source applicability
limits remain binding. The same-time heat proposal remains available within
its original scope.

The original free, elastic and thermal regressions remain passing. Native
viscous tests compare underdamped, critically damped and overdamped oscillators
with independent continuum solutions, demonstrate damping-driven capture,
and check original ports, momentum, energy, force caps and retry. Joint
thermal tests independently integrate continuous damped-oscillator heating,
exercise capture/pull-off with distinct snap/support destinations, and verify
explicit-route refusal plus caloric/cancellation rollback.

### Prescribed source temperature with continuous native and caloric history

`JkrNormalEventModalSystem::continue_at_temperature` transfers the original
held-indentation JKR source at an accepted boundary. It retains all native
positions, velocities, sample count, physical support and complete normal-motion
history. The successor force, area and branch limits drive later motion.
Original contact storage is reevaluated at actual native geometry; its change
is independently checked against the source's signed parameter supply at the
contact energy scale. Native total-energy and force limits remain active.

`JkrSnapHeat::propose_temperature_continuation` transfers the exact sealed
motion predecessor and successor while retaining original ordered caloric
objects, receiver identities and geometry, accepted thermal time, heat
partitions, outside destinations and accumulated snap ledger. It increments
the original history version; stale stamps or substituted histories refuse.
No caloric operation, zero-time transport or duplicate release is performed.
`JkrThermalModalSystem::continue_at_temperature` stages this heat-history
successor before publishing the mechanical source, then publishes both
together. Any source, receiver, work or cancellation refusal leaves all
accepted owners and observations unchanged. Last interval frames remain
historical observations of those intervals.

The returned signed coefficient supply belongs to the explicitly prescribed
temperature driver. It changes contact potential once and enters neither snap
nor support heat partition. Coefficients stay frozen within each subsequent
mechanical interval. The three source temperatures remain independently
prescribed; finite-cell temperatures do not acquire an inferred interface
temperature, entropy or autonomous contact caloric closure.

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
and viscothermal duct and reuses it across samples. Nonlinear storage uses the
shared `PreparedNonlinearZoh` exponential discrete-gradient owner described
below. Both routes reuse candidate-state and mean-output buffers at the render
clock. The shared complete-state capacity remains; the separate Dirac route
retains its original nonlinear temporal admission.

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


