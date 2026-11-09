<!-- fs-couple contract, part 7 of 8; a higher number is newer. Index: ../CONTRACT.md. -->
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

### Imperfect moving thermal contact on the same friction residual

`MovingFrictionHeatStep::evaluate_with_resistance` resolves an explicit pair of
area-specific resistances from the massless heat-generation plane to the two
original bodies for every actual swept cell pair. Each leg uses
`A_mean / (d_i/k_i + R_i)`. The resistance split locates that plane; signed heat
partition and its temperature still follow the existing finite-cell equations.
The same source callback may depend on each actual trial body/interface
temperature. Missing, nonfinite or negative resistances refuse, one insulated
body is supported, and two insulated legs refuse. Exact zero resistance uses
the original arithmetic without allocating a new default resistance array.

`MovingPartialSlipSource::area_resistance_m2_k_w` forwards these original source
values through the existing linear, pressure-changing and quadratic-path
transactions. `None` explicitly selects perfect thermal contact; `Some` needs
one A/B pair per current swept patch. The source re-resolves the current cell
pairs as geometry changes, so a patch index cannot become a persistent material
identity. Local friction coefficients and resistance share the same thermal
trial; retained elastic/kinetic histories and physical-loss ownership stay with
the original partial-slip owner. Accepted patch reports retain their actual
resistances. No fixed link, finite layer store or prescribed body heat fraction
is introduced, and no numerical loss becomes heat.

Independent finite-capacity temperature equations check asymmetric moving
resistance with decreasing `mu(T)`, local heat, interface entropy and a successor
step with different patch indexing. Focused regressions also cover exact
perfect-contact replay, an insulating side, malformed-source refusal and late
cancellation. Conductance coefficients remain endpoint approximations over the
original swept residence; spatial/time refinement and source material validity
remain separate obligations.

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

## Physical frequency and waveform slopes of self-excited resonance

`phs_orbit::PreparedPhsAutonomousSensitivity` differentiates the frequency and
complete Fourier waveform of an original self-excited pHS cycle. It borrows the
unchanged source and feedback, uses the existing complete feedback tangent, and
supplies the original first-order ODE frequency column `i*n*X_n` to fs-orbit's
shared phase-bordered sensitivity owner. Storage-owned or explicitly supplied
complete Hessians retain every nonlinear and material-memory coordinate.

`PhsAutonomousParameterDirection` combines original J/R/G directions, a complete
fixed-state storage-effort partial and a required `PhsFeedbackParameterAction`.
The source-bound physical parameter owner computes the reciprocal direct flow
derivative `G^T*h_p+G_p^T*e`; feedback consumes it and supplies the TOTAL fixed-
state input derivative. The separate state chain is already owned by the cached
total feedback Jacobian. Signed resistance directions need not be positive;
original passive resistance and external active supply remain distinct.

One cached factor serves multiple physical directions without source perturbation
or nonlinear re-solving. Numerical and adapter payloads are admitted together.
Source callback work remains explicitly excluded from the numerical work cap,
with calls bounded by the declared grids. Returned coefficients and frequency
slopes are per physical parameter unit in the fixed authored phase chart.

Three native adapter tests compare simultaneous storage, interconnection,
damping, port and feedback changes with an independent oscillator containing
two explicit relaxation stores, check native/authored Hessian parity, and verify
refusal/cancellation followed by factor reuse. These and the ten forced-parameter
regressions pass together against the current complete fs-phs/fs-orbit owners
and actual coupled numerical bodies. The full Cargo graph remains blocked by
the missing sibling manifest. These local finite-order slopes establish no
stability, infinite harmonic resolution or constitutive-validity guarantee.

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

## Physical-parameter sensitivity of an original forced periodic response

`PreparedPhsPortSidebands::forced_parameter_tangent` and its factored
counterpart reuse the same zero-carrier periodic factor for
`D*x_p=(J-R)*h_p+(J_p-R_p)*e+G_p*u+G*u_p`, where e=grad H and
h_p=partial_p grad H, at fixed original physical coordinates and drive clock.
`PhsForcedParameterDirection` supplies complete skew/symmetric J/R directions,
the same reciprocal G direction for excitation and observation, an optional
complete storage derivative, and optional original-port Fourier input
derivatives. The source-bound fs-phs owner contracts every direct flow term
before rounding. Every nonlinear and material-memory coordinate remains.
`forced_storage_parameter_tangent` and its existing control/result names
remain compatibility wrappers/aliases of this same implementation.
SOS sources expose `SosModalStorage::gradient_parameter_derivative` for exact
original frequency, stress-coefficient and symmetric coupling directions.

The shared independently refined FFT observer first resolves the complete
parameter forcing, then observes `G^T*(Hxx*x_p+h_p)+G_p^T*e`. Direct terms
enter once on the real envelope; both signed Fourier coefficients remain.
Measured forcing-grid change and omitted forcing use the original equation
scales and share the native omitted-equation allowance with the solved response.
This is finite-grid evidence, not an infinite-tail or forward-error bound.
Both projections share work, payload and Hessian caps and one real parameter
sample cap, including any original-gradient/storage callbacks. Input directions
retain original clock/ports and must resolve their declared bandwidth.

Periodically driven feedback uses the SAME original total state Jacobian and
zero-carrier factor. Its required `PhsFeedbackParameterAction` receives original
state, `y=G^T*e` and the direct `y_p|x=G^T*h_p+G_p^T*e`. It supplies the COMPLETE
fixed-state feedback-input derivative, including explicit parameter dependence
and the flow chain exactly once. The imposed Fourier drive derivative remains
separate. `G_p*u` uses both original feedback and imposed drive; `G*u_p` combines
both direct derivatives before contraction. No state chain is added twice.
An explicit zero callback declares zero complete feedback derivative; omission
refuses even for a fixed feedback law. Feedback actions are inadmissible on an
ordinary forced source. Every direct feedback harmonic shares the same measured
forcing-grid/tail allowance and physical work/payload/parameter-sample caps.

Autonomous clocks, changed coordinates/clock and missing chain rules refuse.
Singular, unresolved, source-refused or cancelled operations publish no partial
tangent and leave the factor reusable. There is no finite parameter interval,
perturbed production solve, regularization or branch switch.

Ten full parameter regressions pass through the actual native production
owners, extending the previously executed coefficient/frequency/sideband lane.
These include analytic resonant
damping, reciprocal gain, complex drive and newly excited harmonics; simultaneous
nonlinear J/R/G/storage/input changes compared with independently re-solved
original models; Maxwell time/damping with original memory; and admission,
forcing tails, cancellation and retry. Added feedback oracles independently
differentiate mass, stiffness, damping, gain, feedback, drive and Maxwell time;
resolve nonlinear storage and cubic flow-feedback changes against separate
original perturbed HB solves; and check explicit zero/missing chains, unresolved
feedback harmonics, source refusal and reusable-factor cancellation. The lane uses the
complete fs-phs owner and original coupled FFT/periodic numerical bodies with
real numerical dependencies, without executor substitutes. This does not
establish the blocked full Cargo workspace, continuum or stability qualification.

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
`DepositionSimulation` together for one lattice step. The multirate
`advance_nsf_solid_interval` entry advances a declared positive number of native
fluid steps against one solid endpoint. Every candidate retains the original
fluid timestep and chronological flux histories; no large fluid step or heat
extrapolation replaces them. Each unique native fluid
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
native fluid trial, including derivative probes and final acceptance. An
independent whole-call cell-step budget reserves each full native interval
before computation, including solid cells in the fluid grid. Reports retain
accepted fluid steps, physical duration and all reserved probe work. Each
candidate copies the fluid state once and retains one face vector, independent
of the number of fluid steps; it does not retain a full substep trajectory.

This is fixed-geometry, fixed-inventory, first-order coupled time sampling:
the solid endpoint boundary is held over the chosen native fluid interval.
Both clocks use the same representable lattice timestamp; the thermal residual
uses the original step count times lattice duration, even when subtraction of
rounded absolute timestamps differs. Larger solid update intervals require
independent time refinement; nonlinear tolerances do not bound that sampling
error. Native thermal events still require explicit splits/processing.
Molecular/local SGS, held open controls and the original periodic/channel
dynamic SGS closures are admitted. Dynamic trials retain accepted transport
and flux history while reusing deterministic filter scratch under the original
controls, geometry and retained byte bound. Channel fits consume each proposed
wall temperature; modeled heat remains on interior shared faces and molecular
conduction owns the wall exchange exactly once. Manufactured fluid source
callbacks still refuse because their external state cannot participate in the
transaction. Original exploratory SGS limits, source calorics, fluid
Mach/temperature/stability bounds and chosen spatial/lattice
resolution remain authoritative. No moving solid, boundary-layer wall function,
global time-error bound or experimental convection qualification is implied.

Focused native regression sources cover an independent uniform-fluid and finite
wall capacity balance with series resistance, separate heating/cooling patches,
solid conduction and remaining ambient area, conductivity/insulation response,
source identity and work limits, and transactional late source/cancellation
refusal. Additional G1/G3/G4 definitions compare a multirate interval with
independently solved per-fluid-step heat equations and two finite capacities,
one-step equivalence, continuation, shared work exhaustion and late cancellation.
These regression sources have not been executed on the current
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

### Finite solution baths with compatible swelling feedback

`MoistureLayeredSolid::advance_spatial_diffusion_paired_solutions` couples two
independently finite aqueous inventories to the original nonlinear spatial
moisture/solid transaction. Each exterior face explicitly selects one bath or
is sealed; a sealed physical transfer law cannot be opened by an assignment.
The original layered constructor supplies common temperature and material
activity support. Bath temperature, solvent molar mass and exposed source
context are admitted before native transport, and each bath retains its own
solute totals, nonideal activity law and optional equilibrium/kinetic phases.

The original paired solvent roots now share one private callback transaction
with ordinary slab/Cartesian/radial paired transport. The first-bath inner and
second-bath outer order is unchanged. Every activity pair restarts the original
field and phase stocks. Native compatible mechanics evaluates source swelling,
elasticity, reference wet inertia and the signed prestressed vibration pencil.
The same candidate's actual inverse axial stretch squared enters bulk flux;
its actual current area enters finite surface/contact conductance. The dry
skeleton and water capacity stay on their original reference inventory. Neither
independent cell expansion nor a refitted diffusivity replaces these owners.

Only each bath's assigned signed face transfers enter its solvent balance.
Separate bath checks, field/face conservation, whole-system water increments
and a compensated direct sum of represented cell/bath endpoint/start inventories
must pass before either solution, the field or its actual solid state returns.
One explicit transport-proposal cap covers all nested bath trials; a separate
cumulative native-solid cap includes rejected geometric and bath proposals and
is reserved before each solve. Source, geometry, work and cancellation failures
publish no partial state.

This remains first-order isothermal exchange through the equal-pressure
vapor-equilibrium bridge, with current-frame scalar diffusivity and the existing
small-strain separable-prism reduction. Shear, bending or taper beyond the
specified geometry tolerance refuses. No immersed-liquid permeability, osmotic
traction, chemical-potential stress feedback, displaced bath-vessel volume or
thermal feedback is inferred; reservoir chemical work is outside the mechanical
energy account. Seven authored G1/G3/G4 groups specify independent length/area
flux balances, native source-driven shape/modal scaling, heterogeneous contact
through-flow with finite crystals, zero-swelling compatibility, retained wet/dry
histories, sealed kinetics and cumulative-work/source/cancellation refusals.
All eight changed Rust files passed the pinned formatter/parser check; the
integration groups remain unexecuted because dependencies block fs-couple.

### Adaptive compatible swelling with finite solution baths

`MoistureLayeredSolid::advance_spatial_diffusion_paired_solutions_adaptive_to`
uses the original first-order full/two-half controller around the complete
native geometric transport, mechanics and two-bath transaction. Every probe
retains the same original solid binding, reference mesh, dry skeleton and
swelling datum. Accepted final stress, compatible geometry, wet inertia and
the prestressed modal pencil come from the actual fine endpoint.

Temporal error includes each layer's moisture/activity, an L1 discrepancy
weighted by original dry cell mass, each bath's independent solvent/activity
and optional dissolved/crystal stocks, and every native nodal displacement.
Displacement has an explicit metre scale; relative scales use changes from the
last accepted state. The existing representable-half correction is applied
once to the combined comparison. This is a local asymptotic estimate, separate
from equilibrium/geometry tolerances, and does not certify global error.

Only accepted fine transfers enter layer/interface and assigned bath ledgers.
Whole-window original-to-final layer, independent bath, component/phase and
direct cell/bath inventory audits retain fixed budgets regardless of step count.
The field cap includes every nested root trial; the native-solid cap spans the
whole duration, including the initial comparison datum and discarded or failed
probes. Original native state snapshots share their sparse assemblies. A source,
geometry, work or cancellation refusal returns no accepted prefix.

The original isothermal vapor-equilibrium bridge, scalar current-frame
diffusivity and small-strain separable-prism scope apply. Three authored G1/G3/G4
groups cover displacement-only refinement against an independent moving-length
continuous solution, accepted-fine state/shape and opposing finite-bath accounts
with exact work replay, cumulative caps after rejected comparisons, and sealed
heterogeneous redistribution with two original kinetic phase histories.
All seven changed Rust files passed pinned parsing/formatting checks. These
integration groups remain unexecuted because dependencies block fs-couple.

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

## Elastic reflection coupled to local support mechanics

`elastic_reflection::OscillatorSupportedSurface` retains the original elastic
material state, two explicit `fs-phs::WallPin` supports and their common angular
frequency. The original wall impedance supplies normal and isotropic tangential
boundary response under `exp(-i omega t)`; `fs-material` owns P/SV/SH reflection,
mode conversion, coherent boundary traces, strain limits and positive support
work. Preparation reuses the native bounded solve and cancellation points.
Repeated coherent drives allocate nothing and perform no new factorization.

Each drive also returns support displacement and cycle-mean kinetic and spring
energy per unit surface area from the same physical velocity. The native
absorbed power remains the sole loss account; reactive storage is not heat.
Lossless resonance admits zero impedance without a resistance floor. The
reduced support is stationary, linear and locally reacting; its authored mass,
stiffness and damping do not assert measured material fidelity, resolved layers,
nonlocal shell response or a transient trajectory. G1/G3 regressions compare
independent normal reflection, resonant absorption, reactive work/storage,
coherent mode conversion, material contrast and support-mass changes.

`OrientedOscillatorSupport` retains three original `WallPin` laws and an admitted
`fs-ga::So3` principal-to-surface rotation. Its reciprocal impedance is
`R diag(Z_i) R^T`, assembled with double-double product/reduction arithmetic.
The original `R diag(sqrt(r_i))` dissipative channels are retained directly,
so refactoring rounded Cartesian coefficients cannot invent principal loss.
`OrientedOscillatorSupportedSurface` passes this tensor to the native material
boundary, allowing direction-dependent reflection and SH/P/SV conversion from
tilted supports. The prepared native velocity map is reused on every drive;
no reflected-field subtraction or new factorization reconstructs wall motion.

Principal velocities and displacements retain each pin's original kinetic and
spring potentials through the same storage evaluator as the two-pin owner.
Surface displacement remains in the material boundary frame. Native absorbed
power is the sole dissipative account, including rank-deficient and lossless
resonant supports; reactive storage is not heat. This is a common principal
frame for the authored mass, stiffness and resistance, not resolved layers or
independently oriented noncommuting mechanical tensors. Five focused regression
definitions cover the original diagonal limit, analytic rotated shear and
tilted longitudinal/shear conversion, coherent work and storage, and source
domains/cancellation. Native execution status is reported separately.

### Resolved elastic coatings on original oscillator supports

`elastic_reflection::OscillatorBackedElasticCoating` joins an original incident
solid and an original finite elastic coating to the same three principal
`WallPin` laws and checked `So3` rotation. `fs-material` owns the complete
coating propagation, P/SV/SH conversion, internal interference, critical and
evanescent fields, and both original strain domains. The coupling layer passes
the original support tensor to that native boundary, without fitting a scalar
coating impedance or discarding internal modes. Frequency is authored in Hz;
both the coating propagator and wall law use the same `TAU * frequency_hz`.

The actual rear support velocity feeds the same principal-coordinate motion
and storage evaluator used by the direct oriented surface. A front-face field
is never substituted for rear motion. Returned kinetic/spring energies belong
only to the original pins; coating bulk storage is not counted as pin storage.
The native positive support work is the sole dissipative account. Both original
solids remain lossless, isotropic and perfectly bonded, and no transient clock
or nonlocal support response is inferred. Exactly zero thickness recovers the
direct supported surface and has no interior coating field. Prepared drives
and native interior field queries allocate nothing and reuse their factorizations.

Five focused G1/G3/G4 regression definitions cover exact absence, homogeneous
propagation with polarization conversion, quarter-wave matching and rear
storage, a full rotated-support fixture against an independent 70-digit
continuum transfer exponential, and source-domain/cancellation refusals.
Native compilation and execution status are reported separately.

### Causal transient reflection from the original support motion

`elastic_reflection::PreparedTransientElasticReflection` retains all six
principal displacement/velocity coordinates of the original oriented pins.
Normal-incidence real P/SV/SH characteristics couple through the original
material impedances to those coordinates. The existing native
`PreparedQuadraticFlow` owns exact held-input evolution, exponential and
trajectory-integral preparation under explicit dimension/byte/work limits.
Warm steps use fixed-size buffers and perform no factorization or allocation.
Compatible prepared durations can continue the same original state.

Original pin heat, outgoing-wave energy and changes in kinetic/spring storage
remain separate accounts on one trajectory. Outgoing-wave work uses independent
quadratic and linear observations, never a residual energy balance. Real
instantaneous wave power has no phasor cycle-average factor. A pulse can leave
stored support motion that later radiates even when incident input is zero;
lossless pins emit this energy without creating fictitious heat.

The source-domain check bounds incoming/emitted strain throughout the interval
using the passive envelope `H(t) <= H(initial) + incoming wave energy`, rather
than checking only endpoints. Refusal or cancellation preserves accepted
motion and time. This remains a stationary local boundary at normal wavevector
incidence in a lossless solid; oblique history kernels, resolved coating time
propagation and previously emitted spatial fields belong to their own owners.
Five focused regression definitions cover the analytic pulse and causal tail,
rotated interval composition, lossless radiation, transactional/source/budget
refusals and interval-wide strain admission. Native execution is reported
separately from formula, source and formatter checks.

### Original anisotropic modes in harmonic and causal support reflection

`elastic_reflection::NormalOscillatorSupportedSurface` binds the material's
`NormalElasticCharacteristic` to the original three principal `WallPin` laws.
The native full-tensor Christoffel modes retain their ascending-speed order,
polarizations, strain domain and checked material/surface frame. No isotropic
projection or shear/longitudinal relabeling substitutes for that source. The
same polarization describes each mode in both propagation directions; the
legacy P/SV/SH API keeps its original direction-dependent convention.

Harmonic preparation reuses the native passive tensor boundary solve. Actual
surface motion feeds the existing principal-coordinate displacement and pin
storage evaluator, with support work remaining the sole dissipative account.
`PreparedTransientElasticReflection::prepare_normal_characteristic` connects
those same incoming characteristics to the existing exact held-input pHS flow.
Modal projections form the original radiation impedance and forcing directly;
the six physical support coordinates and fixed-size warm workspace are shared
with the isotropic path. Source identity, support values and surface frame bind
accepted motion; compatible physical states may pass between channel bases.

Anisotropic energy may walk laterally even when the wavevector is normal. This
boundary retains the original local normal work and material-frame strain,
without claiming finite-aperture transport or an oblique anisotropic solution.
Four focused G1/G3/G4 regression definitions compare a mixed-polarization SPD
tensor against independent 70-digit Cartesian harmonic and direct physical
transient references, recover the legacy isotropic fields and memory, and check
transactional frame binding and cancellation. Native test execution is reported
separately from these references and source/formatter checks.

## Compact outlet reaction from an impedance wall

`compact_wall_radiation::PreparedCompactWallRadiation` closes original finite-pressure or prescribed-flow network sources against the complete exterior impedance of one or more unflanged compact circular outlets. Each diagonal retains `fs-duct`'s native unflanged radiation impedance, every off-diagonal retains the original 3D monopole mutual impedance, and every pair, including the finite diagonal self return, receives the original spherical impedance-wall response. These use the same fluid, frequency, peak phasor and physical volume-flow units. Reciprocal pairs are evaluated once; all reflected preparations share one aggregate work allowance and retain their native error diagnostics. Cancellation or refusal publishes no partial matrix.

The complete free-space and wall-loaded matrices use the native full Hermitian passivity admission. The signed wall correction is never required to be passive by itself: a passive wall can reduce radiation resistance through interference. No artificial resistance, spectral clipping or fitted source impedance repairs a load. `drive_network` invokes the original duct source solver with that complete coupled termination, preserving finite generator impedances, source dissipation, mutual radiation work and the original network power account. `observe` maps the solved outward outlet flows (`Q=-U`, since network currents point inward) into the matching native total-field array; repeated observation needs no quadrature.

The caller retains and binds the actual outlet labels, radii and positions. Explicit positive `ka <= 0.5` and radius/distance bounds (at most 0.25) admit the compact approximation; smaller authored limits may be required. Wall clearance and pair separation are checked without invented distances. These bounds are not approximation-error certificates. The original wall remains infinite, stationary and locally reacting; outlet orientation, pipe-body scattering, finite-aperture multipoles, mean flow and nonlinear source changes remain outside this model. Network/exterior gas snapshots must agree.

Four focused regression definitions check independent rigid/release image impedances, signed reflected resistance with passive complete loads, single- and two-branch pressure-source feedback against independent transmission-line equations, complete work balance, actual outward-flow observation and shared quadrature/geometry admission. A bounded native probe passed these four plus seventeen existing acoustic/Bessel regressions using the actual wall, free-field, complete-load admission, source closure, LU and work-accounting bodies with real `fs_math`. Its lossless cylindrical graph assembly used an analytic nodal-admittance fixture, and gas construction supplied a fixed physical snapshot; full native graph/material admission and the full Cargo crate suite were not executed. The repository tests retain the original graph and gas constructors.

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

### Native second-order heterogeneous diffusion against finite solutions

`advance_layered_sdirk2`, `advance_shared_layered_sdirk2`,
`advance_layered_radial_sdirk2` and `advance_shared_layered_radial_sdirk2`
now solve the actual finite solution and original material layers at both SDIRK
stages. The first stage uses the existing native shared BE/root. Its converged
physical face loads are captured once without another nonlinear solve. The
second root starts from the ORIGINAL bath and body storage with those retained
weighted loads; its own activity drives the final native faces.

Both native geometries share the existing profile preparation, nonlinear
storage and stage solver. Each material retains its original sorption/source,
dry mass, contact resistance and exact slab or curved path. One private stage
participant supplies complete heterogeneous bodies to the original solution
root. Intermediate body copies are released after load capture while retaining
their stage diagnostics. Internal transfers have one adjacent-material debit
and credit; only exterior exchange debits finite solvent.

The endpoint report contains both RK-weighted contributions once. Diagnostic
first-stage water is not an additional debit, and root work includes both
solves. Canonical shared IDs and original source, clock, cell, phase and water
gates remain. Algebraic phase partitioning is native; kinetic crystal evolution
is refused because its endpoint-volume owner has no corresponding second-order
stage equations. Geometry and temperature remain held, and smooth-branch local
order does not certify behavior across a phase event.

Five authored regressions cover independent finite-bath stage equations and
weighted faces, second-order continuous-reference refinement for all geometries,
shared opposing exchange, algebraic crystal buffering, and late stage/source/
budget refusals. Rust formatting and independent reference calculations passed;
the new native integration tests remain unrun.


## Native relaxing-bond temperature continuation (G2)

The native intact adhesive owner accepts an explicit temperature successor only
through the original source gate. Ordered surfaces, card identity, pinned
coefficient claims, thickness, branch indexing, nonthermal coordinates and
runtime limits remain fixed. At held actual jump, the material owner retains
branch strains, past mechanical work/heat and elapsed time. Changed moduli
supply signed parameter work; changing relaxation times alone supplies none.
Prepared ramp coefficients are rebuilt for the existing native timestep.

`apply_temperature_successor_with_admission` stages that parameter-work account
before publishing the source, retained point and prepared coefficients together.
Native displacement, velocity, sample count and the last mechanical frame stay
at their accepted values. The separate temperature report retains constitutive
storage increments even beside large unchanged preload and admits no invented
heat or thermal expansion. The next native timestep uses the refreshed law and
the same branch memory. G2 regressions in
`tests/cohesive_modal/relaxation/temperature.rs` cover actual continued body
motion, source/receiver refusal and late cancellation with replay.

## Combined source Kelvin–Maxwell loss on the primary MITC shell

`TensorShellMaterial` accepts the original resolved tensor's explicitly sourced
Kelvin matrix alongside its Maxwell branches. `TensorShellRelaxation` uses the
same five-strain constrained owner, original MITC trial fields, positive area
and thickness quadrature, source rotation and shear-energy correction. A
viscous thickness-normal constraint contributes its actual additional history;
it is not condensed independently at each branch or replaced by a modal decay.
The original equilibrium A/B/D/shear and all inertia checks remain binding.

The complete instantaneous modal viscosity is summed with compensated work
arithmetic and factored by the shared original-unit passive matrix owner.
`kelvin_connections(component)` supplies its existing bilateral dashpot
channels for cold coupled construction. `attach_to_with_kelvin_connections`
then checks the exact contiguous connection range, every complete signed modal
column, zero spring/rest extension and original resistance. Constant-loss
substitutions refuse. The original `attach_to` refuses a nonzero feedthrough
unless the explicit Kelvin binding is used. Kelvin heat remains in the native
connection account; Maxwell and viscous-normal relaxation heat remain in the
native material account. Both use the same existing exact held-force evolution.

Source observations use actual displacement, velocity and complete reconstructed
history, including the displacement-dependent physical normal strain. Initial
and staged source checks retain the full original strain domain. A displaced
body with retained normal viscosity refuses the ambiguous zero-history
`Unrelaxed` shorthand; fully relaxed initialization and zero-displacement
unrelaxed initialization keep their exact physical meanings. The caller still
owns spatial modal signs, degenerate-basis identity and avoiding duplicate
physical material bindings. Source-band and runtime work/state/energy gates
remain unchanged; integration-point checks do not certify all spatial or
between-step peaks.

Focused native regressions compare every modal cross term to direct complex
six-dimensional source reduction and original spatial quadrature, compare
combined exact motion with independent continuous equations and separately
integrated Kelvin/relaxation heat, check full channel admission and retain
zero-Kelvin Maxwell identity. They exercise the actual material, geometry,
projection, passive factor and pHS propagation owners; they do not substitute
for execution of the complete coupled rendering wrapper.

## Physical forced-response branches through frequency folds

`phs_orbit::continue_forced_port_response_scaled_seeded` composes the existing
full-state pseudo-arclength, physical port-spectrum and optional variational
Floquet owners. Its analytic counterpart admits the complete continuous storage
Hessian, including every nonlinear coupling and explicit material-memory row.
The same declared Hessian enters the original continuation state Jacobian and
optional shooting tangent. Storage-owned exact quadratic metadata or continuous
tangents are discovered automatically; an unavailable or refusing tangent does
not become a quadratic surrogate or an undeclared state finite difference.
Both original forced-frequency adapters supply the exact original parameter
column `dR_n/domega = i*n*X_n`, with zero DC rows, for their fixed-phase
`PeriodicPortInput` and identity frequency map. Every explicit material-memory
coordinate remains. The shared corrector applies its original equation and
parameter scales and arclength row; a declared refusing, incomplete or nonfinite
derivative never falls back to a frequency probe. When both state and parameter
columns are analytic, the unchanged source samples need no restoration replay.
The cancellable branch polls before the analytic parameter column.

The forcing's original DC values, harmonic amplitudes and relative phases stay
fixed while its common fundamental changes. All Fourier coordinates, memory
states, normalized HB/AFT/arclength defects and physical cycle work are retained.
Returned points remain in continuation traversal order, including frequency
reversals; they are not sorted into a single-valued frequency graph. Each point
observes the original nonlinear `G^T gradient(H)` through the existing shifted
FFT-grid refinement. The output preserves real DC and negative-time peak phasors
`2 conj(c_k)`; output order is independent of retained state order.

Requested stability starts all three existing variational time grids from that
point's complete physical phase-zero state at its own authored clock. Every
original multiplier remains present. Physical seed-agreement, grid comparison
and unit-circle guards retain their existing semantics, including `Unresolved`.
The returned HB state is not replaced by a shooting correction. Omitted stability
returns `None`, never an implicit stable label. There is no global attraction,
unique-branch, exact bifurcation-location, infinite-tail or physical hysteresis
trajectory claim. Source laws remain smooth, deterministic and read-only.

The original continuation/Newton/AFT caps bound branch construction. At most
`max_steps+1` points receive an audit and bounded spectrum query, plus three
individually bounded variational queries when requested. A conservative bound
for all retained numerical-result arrays is checked before source evaluations;
it includes every optional monodromy, spectrum and refinement trace, excluding
allocator metadata and separately owned transient workspaces. The original ODE
adapter is wrapped only to poll cancellation at AFT source/tangent evaluations;
cycle audits poll per sample. Existing observation/shooting cancellation remains
active. Native dense solves and individual source calls remain atomic. A refusal
or cancellation publishes no partial branch and leaves the authored input intact.

Four G1/G3/G4 integration regressions cover independent original Maxwell phasors,
physical work and all three decay multipliers; two frequency folds with original
nonlinear restoring effort checked by independent Fourier convolution and a
complete memory Hessian; a manufactured nonlinear periodic response with stable
and parametrically unstable original Floquet maps; and tangent, payload,
cancellation and deterministic-replay boundaries. These tests are authored but
not claimed executed while the full `fs-couple` dependency closure is unavailable.
Formatting and the bounded native-API/physics review are separate from numerical
execution; no complete workspace or experimental-validation claim is made.

### Original conductivity and convection on retained porous gas walls

`ThermalPorousGasSystem<ClosedGasParcel>::with_source_contacts` adds original
`SolidGasContact` patches to the existing finite-gas, porous-compaction and
thermal-feedback owner. It reuses `PreparedContact` for two-layer conduction,
source k(T), forced/natural convection, and the actual solid-surface heat
balance. Original correlation domains and evidence remain authoritative.

The original material card, exact caloric curve and complete fixed query
context come from each retained mechanical family. Conductivity is selected
at its original stress-free reference and remains pinned. Integrated pure-phase
conductivity profiles are prepared once and shared across thermal probes and
accepted compression epochs. Each heated gas endpoint supplies one original
isobaric property law; native thermal substeps resolve the actual cell state.
No gas conductivity, heat capacity, density, pressure or material source is
replaced with a constant surrogate.

Source patches and fixed W/K links share one 16,384-entry cap and unique ID
space. Patch areas must fit the original exposed cell area, and source-contact
cells cannot also apply native ambient convection. Explicitly declared area,
layer distances and correlation geometry remain held reference-network data;
this does not infer a moving boundary layer or deformed Fourier operator.
The initial source-contact entry admits ideal parcels because that is the
original film-property owner's supported gas law. Existing Peng--Robinson
fixed-link heat and native dense-fluid pressure/caloric behavior remain intact.

All transfers enter the same compensated per-gas, per-cell and complete energy
accounts as fixed links. Source convection contributes its actual transferred
power divided by the gas/cell contrast to the native explicit decay budget;
its series diagnostic is not silently used as a nonlinear tangent. Reports
retain actual source queries, surface temperatures, source profiles and signed
heat integrals. Late source-domain failures, native solve refusal and
cancellation publish no gas, thermal, mechanical or observer history.

Four focused regression definitions compare an independent finite-capacity
series-conduction balance, integrated affine k(T) and the original sphere
Nu=2 limit across two heated compaction epochs, retained profile identity,
late conductivity-domain rollback/retry, and patch coverage/identity admission.
They retain the existing first-order heat-before-motion and thermal-substep
accuracy boundaries; energy closure is not a time-accuracy or calibration
claim. Validation status is reported with the implementation commit.


## Native signed adhesive release and recontact (G1/G3)

The existing DMT--Mindlin native owner accepts finite ordered signed indentation
brackets. The source's original affine contact-window integration owns release,
zero-load flight and reattachment; normal adhesive attraction may persist
while the Hertz contact radius and friction vanish. Native observers retain the
active contact fraction and normalized interval independently of endpoint
force. Total macroscopic slip includes flight, while physical microslip heat
and branch evolution belong only to actual contact.

All normal/shear cross compliance, actual-body force/work gates and atomic
history publication remain active. The whole-step mean force drives the
original bodies; active-window force is not applied for the whole timestep.
G1/G3 regressions in `tests/adhesive_dmt_modal/mindlin/gaps.rs` compare complete
release/flight/recontact and reversal trajectories with an independent coupled
two-mass solution, preserving all three momenta and complete work. Initially
separated contact also supports receiver refusal and late-cancellation replay.
Paths remain affine with one proportional shear direction; callers resolve
within-step reversals and supply the physical root bracket.

### Original spatial structural/acoustic resonance responses

`cavity_modal::interior::structure::StructuralInterior` now recovers physical
harmonic storage and supports analytic frequency response and driven receiver
surveys on its existing complete spatial pencil. Original structural mass,
stiffness, damping and free-coordinate order, heterogeneous fluid inverse-density
and compressibility forms, and the reciprocal interface map remain unchanged.
No structural or pressure mode is selected, normalized, truncated or refitted.

Every ordinary response carries separate structural kinetic/strain and fluid
kinetic/compression storage. Peak-phasor energies are `v^H Ms v/4`, `u^H Ks u/4`,
`psi^H Ka psi/4` and `p^H Ma p/4`. The acoustic coordinate is pressure potential,
not pressure: under the owner's `exp(-i omega t)` convention, `p=-i omega psi`.
Original signed interface force and displaced flow share B and B^T. Compensated
contractions preserve coherent cancellation; small signed energy roundoff stays
visible. Structural and acoustic source/exchange/heat balances are independently
observable, and no heat or storage is assigned from a work deficit.

`PreparedStructuralInterior::solve_with_frequency_derivative` reuses the existing
canonical harmonic factor and its original-coordinate residual gate. Returned
derivatives include displacement, velocity, pressure potential, physical pressure,
fluid force and displaced flow, all per rad/s. Both force and flow source slopes
are admitted explicitly; None fixes both phasors. The pressure derivative includes
`-i(psi+omega psi_omega)`. Independent power derivatives retain both source chain
rules, reciprocal exchange and the original symmetric structural/relative-wall
dissipation. Gyroscopic/interface skew work is excluded from heat.

The additional field-recovery budget is `256*n+256` bytes and work score
`128*n*n+256*n`, on top of the native harmonic-frequency budget; n includes every
structural and pressure coordinate. These bounds exclude caller arrays, the
already prepared owner and allocator metadata. Preparation is reused across
drives, and a response plus tangent requires one factorization per frequency.

`StructuralInterior::survey_resonances` uses the existing shared peak search with
fixed coherent source phasors and either an original pressure point or a real
weighted original structural displacement/velocity receiver. Observations retain
the original negative-time phasor convention; conjugation belongs to effort/flow
work, not the observation row. Queries remain in rad/s without a Hz round trip.
Refinement, peak/query caps, native singularity refusal and cancellation preserve
the shared transactional behavior. The survey reports grid-resolved driven peaks,
not a complete free-mode census, half-power width or decay Q. Material, geometry,
small-motion coupling and real memoryless wall-admittance assumptions remain those
of the original owner; no exterior radiation or new causal material law is inferred.

Four native tests executed the unchanged production response-contraction module
against the real `fs_math` owner and passed: physical cycle-averaged storage with
distinct original mass/compliance metrics; reciprocal signed interface and
coherent cancellation; original source/heat derivatives excluding gyroscopic work;
and invalid-data/cancellation/replay. Three public FE regressions additionally
specify complete differentiated source rows and work channels, distinct coupled
physical receiver peaks and phases, and source-tangent/resource-refusal behavior.
Those FE regressions are authored but unexecuted while the complete `fs-couple`
dependency build is unavailable. MR56's broader qualification remains open.

`StructuralInterior::prepare_resonance_modes` now prepares the same complete
original M/C/K arrays through `fs-modal::quadratic_modes::PreparedQuadraticModes`.
Repeated `refine` calls recover seeded local free modes without a bare-fluid
inverse, pressure-mode truncation, proportional-damping approximation or new
eigensolver. Exact binary equilibration and the full mass factor are reused.
Every structural displacement and pressure-potential coordinate survives,
including constant fluid potential; original-row residuals and mass
normalization are checked by the shared physical quadratic-mode owner.

The free-mode convention is `exp(s t)`: original velocity is `s*u` and physical
pressure is `s*psi`, retaining the real decay term. Returned pressure, fluid
force `B^T*p`, displaced flow `B*v`, oscillation frequency `abs(Im(s))` and signed
amplitude decay `-Re(s)` describe this coupled mode. Original structural and
fluid storage use Hermitian modal-envelope factor 1/2, not harmonic-cycle
factor 1/4. Structural loss uses only the symmetric part of the original
velocity operator; wall loss and signed interface work are independently
evaluated. The reported decay residual is `2*Re(s)*H + D_structure + D_wall`.
Interface and gyroscopic skew work create no heat. No linewidth is converted
to decay, and growth is never clipped to stability.

Shared `PoleOptions` additionally reserve `128*n+256` numerical payload bytes
and work score `128*n*n+256*n` for original seed, field and physical observation
recovery. The limits cover preparation and each refinement independently,
excluding the borrowed original model, caller arrays, earlier results and
allocator metadata. Cancellation or refusal publishes no partial mode and
leaves preparation reusable. Exact stiffness symmetry and original physical
storage checks remain in force. These are local numerical free modes of the
frozen spatial model, not a complete spectrum, continuum enclosure, nonlinear
resonance prediction or arbitrary frequency-dependent wall realization;
repeated or defective roots retain the native solver's refusal.

`structure/poles_tests.rs` specifies independent damped/overdamped scalar roots,
growth sign, a two-coordinate gyroscopic frequency split without dissipation,
exact affine-cavity compression loading, the same free trajectory under the
existing time propagator, complete source-row and separate envelope-work
identities, spatial wall losses, and budget/seed/cancellation/replay behavior.
These full-FE tests are authored but unexecuted in this lane while the complete
`fs-couple` dependency closure is unavailable. Their definitions alone do not
qualify MR56 or establish performance.

`PreparedStructuralInteriorModes::refine_with_sensitivities` refines one genuine
native coupled mode and prepares the shared quadratic sensitivity factor.
`PreparedStructuralInteriorSensitivities::parameter_derivatives` reuses that
factor for a complete nonempty batch of caller-declared physical directions.
Each direction supplies full signed derivatives of original structural
mass/damping/stiffness, acoustic mass/stiffness, symmetric wall damping and
pressure-row/structural-column coupling. One dB assembles both opposite
interface blocks; structural gyroscopic and nonproportional terms remain.
Constitutive and geometry derivatives, including trace/normal/measure changes,
must be supplied explicitly.

Results retain original-mass normalization and one fixed original phase
coordinate. They include complex ds, signed decay, all original field tangents,
dp=ds*psi+s*dpsi, dF=dB^T*p+B^T*dp and dQ=dB*v+B*dv. Both interface
product-rule terms accumulate together before rounding. Original-polynomial
and normalization residuals accompany each result. A simple real pole retains
its complex ds without inventing a derivative of abs(Im(s)) at zero. Repeated
or unresolved poles and invalid/over-budget/cancelled directions refuse;
no additional eigensolve, finite-difference query or damping floor is introduced.

For n=ns+na and k=options.numerical.max_vectors, the adapter reserves
24*k*n*n+128*(k+1)*n+512*(k+1) physical bytes and work score
128*(k+1)*n*n+256*(k+1)*n for assembly, the retained baseline and batch output
before passing the remaining physical ceilings to the shared owner. Native
numerical-factor ceilings remain separate. Caller directions, original
model/source preparation and earlier batches are excluded. A refused batch
leaves the prepared baseline reusable.

Seven original quadratic-sensitivity tests and an independent full five-DOF
P1 tetrahedral cavity fixture pass in a native lane using the actual shared
numerical owners. The cavity oracle checks structural mass/damping/stiffness,
gas bulk modulus, coupling strength, common coefficient scaling and zero
directions, including the complete pressure derivative. The new Cx-bound
adapter and its two full-FE integration regressions are authored but unexecuted
while the complete dependency graph is unavailable; these numerical checks do
not establish whole-adapter compilation or close the broader MR56 qualification.

### Full-field cavity free modes with original relaxation memory

`StructuralInterior::prepare_relaxation_modes` passes the same complete
original M/C/K arrays and explicit `MaxwellDampingBranch` channels to
`fs-modal::relaxation_modes::PreparedMaxwellModes`. Branch strains use original
[structural displacement, pressure potential] order. K must be relaxed;
branch stiffness must not already be included. Structural channels describe
Maxwell material relaxation; potential channels contribute the explicit
passive admittance k*tau/(1+s*tau) b*b^T. Mixed channels retain the caller's
declared law. No material fit or tensor-shell constitutive image is inferred.

Repeated `refine` calls retain every original field and branch stress through
the full rational pencil. Physical pressure is s*psi, including real decay.
The original B/B^T interface force and flow retain reciprocal exchange.
Relaxed structural/fluid envelope storage, structural viscous loss, memoryless
wall loss, and internal-memory storage/heat remain separate. Each branch
reports |sigma|²/(2*k) storage and |sigma|²/(k*tau) dissipation. Their totals
are subaccount sums, not additional debits. Binary-normalized DD arithmetic
avoids unrepresentable intermediate squares and coefficient products.
The independently evaluated decay residual is
2*Re(s)*(H_spatial+H_memory)+D_structure+D_wall+D_memory.
Native original mechanical and branch-history residuals accompany the result.

For n=ns+na and b original branches, the adapter reserves
160*n+64*b+512 payload bytes and work score 160*n*n+256*n+256*b
before native preparation. These bounds cover original seeds, output recovery
and physical observation; borrowed sources, earlier results and allocator
overhead are excluded. Native matrix/branch preparation retains its own
bounds. Cancellation and refusal publish no partial result and preserve reuse.
This is a seeded local simple-mode calculation. Native refusals at s=-1/tau,
pure material dark modes with zero displacement, unresolved multiplicity and
unrepresentable physical output remain explicit; no full-spectrum or continuum
accuracy claim is made.

Eight native tests pass using the actual shared Maxwell numerical owner,
the original response-contraction module and the exact new memory-energy
helper. Independent original P1 matrix fixtures cover all three analytic
fluid-loaded cubic poles, causal wall admittance, gyroscopic handedness with
mutual relaxation, the zero-memory limit and original energy balances.
Four complete FE adapter regressions are additionally authored but remain
unexecuted while the canonical dependency closure is unavailable. The native
checks do not establish whole-adapter compilation or broader MR56 qualification.

### Excitation residues and parameter slopes of full-field memory poles

`PreparedStructuralInteriorRelaxationModes::refine_residue` resolves one
original right/transpose pair through the existing shared Maxwell owner.
The transpose uses T(s)^T without conjugation: original reciprocal B blocks,
gyroscopic terms and nonproportional damping remain in the nonsymmetric
pencil. Independent original transpose residuals are evaluated at the right
pole, with both equilibrated derivative alignment and original cancellation
admitted against the caller's floor. The physical baseline is observed from
that genuine retained right mode without a second refinement.

`force_residues` reuses the pair for complete batches of original generalized
structural forces and inward pressure-test-function flow loads. Results are
Laurent coefficients of 1/(s-s0), including every original displacement,
velocity, pressure potential, pressure, reciprocal interface force/flow and
branch stress. Pressure residue is s0 times potential residue, including
real decay. Field names preserve the additional inverse-second residue unit.
No analytic background, initial history, other poles, complete frequency
response or full modal expansion is supplied by this local coefficient.

`pole_derivatives` contracts the same pair through complete signed original
M/C/K directions and both reciprocal dB blocks, sharing the existing physical
direction assembly with the quadratic-pole owner. Explicit branch derivatives
include stiffness, d ln(tau) and every original strain entry. Native factored
DD source contractions retain cancellation; no eigensolve, finite difference,
mode normalization or numerical equilibration is differentiated. Returned
slopes include complex ds, signed amplitude-decay change and d|Im(s)| only
for oscillatory poles. A real pole reports no selected absolute-frequency
derivative. Mode-shape and force-residue derivatives are outside this API.

For n original coordinates and b branches, separate physical preparation
ceilings admit 320*n+128*b+1024 bytes and work 256*n*n+512*n+512*b.
A k-load batch admits 256*n+128*b+1024+k*(80*n+32*b+512) bytes and work
k*(64*ns*na+512*n+128*b+512). A k-direction batch assembles one direction
at a time, admitting 24*n*n+256*n+128*b+1024+64*k bytes and work
k*(1024*n*n*(b+1)+128*n+256*b+512). These include the retained physical
baseline and complete current batch, excluding the original/native source,
caller arrays, earlier results and allocator overhead. Original PoleOptions
independently admit both native refinements. Cancellation, ill-conditioned
normalization, invalid directions and resource refusals return no partial
batch and preserve the immutable pair. Material poles, dark material modes,
repeated roots and mode crossings retain the shared owner's local limits.

The native validation lane exercises the actual shared Maxwell residue owner,
the shared original-block direction assembly and physical residue mapping
against the exact five-coordinate P1 fluid-loaded cubic. It covers force and
inward-flow excitation at all three roots, full pressure/force/flow/stress,
material and spatial pole slopes, gyroscopic direction retention, invalid
source dimensions, and cancellation. Three full FE/Cx integration tests are
also authored; the missing frankenscipy/fsci-special manifest prevents their
Cargo compilation in this environment. Native evidence is not a claim that
the complete adapter or workspace has compiled.

### Complete driven cavity response with original relaxation memory

`StructuralInterior::prepare_relaxation_response` binds original complete
M/C/K/B and explicit Maxwell channels through `LinearMechanicalSystem` and
`HarmonicModel`. The canonical state retains every [u,psi], momentum and
viscous-strain coordinate, including constant potential and dark memory.
Original M/K must be symmetric and K relaxed; native passive storage/resistance
admission applies without coefficient repair. Structural-only branches are
material channels; potential-only branches supply `k*tau/(1+s*tau)*b*b^T`
admittance; mixed channels retain the caller's declared law. The complete driven
solution includes off-resonance background as well as resonance.

Native assembly is reused across frequencies. One full LU per positive rad/s
frequency serves coherent original force/inward-flow drives and analytic
frequency tangents. Peak phasors use exp(-i*omega*t); velocity and pressure come
from the power-conjugate outputs. Separate componentwise checks enforce original
mechanical, kinematic and branch-history equations, including differentiated
equations and both supplied source slopes. Pressure's tangent includes the
derivative of -i*omega*psi. No initial-history forcing, mass lumping, modal
truncation, sampled damping, finite-difference solve or extra solver is introduced.

Responses retain original structural/fluid cycle-mean storage, each arm's
`|sigma|^2/(4*k)` storage and `|sigma|^2/(2*k*tau)` heat, structural viscous heat,
memoryless-wall heat, source work and signed reciprocal interface work. Each arm
separately reports work drawn from structure and acoustic field, including mixed
channels. Compensated contractions retain coherent cancellation. Aggregate arm
totals are subaccount sums; structural, acoustic, memory and total power residuals
remain separately observable. No heat is manufactured from a deficit.

Prepared `driven_resonance` and `survey_resonances` reuse the shared analytic
search for actual original pressure, displacement or velocity receivers. Local
queries can return selected squared-amplitude half-height crossings and observed
bandwidth Q, distinct from pole-decay Q. Surveys resolve the supplied grid and
do not establish a complete spectrum. Spatial/material validity and mesh/order
convergence remain caller obligations.

For n coordinates, b branches and d=2*n+b, the native complete-state cap remains.
Separate adapter ceilings admit `8*(32*d^2+16*n^2+8*b*n+64*d+64*n+32*b)` payload
bytes and work score `1024*d^3+128*n^2+128*b*n` before physical construction;
native harmonic ceilings separately admit preparation, factors and derivative
queries. Borrowed arrays, earlier results and allocator metadata are excluded.
Cancellation, true singularity, invalid sources, original-row failure or capacity
refusal returns no partial response and leaves preparation reusable.

Ten native tests pass with actual shared owners and unchanged production
recovery: six new cases cover the full five-coordinate P1 cubic response and
accounts, both source slopes, causal wall admittance, actual peak/bandwidth,
nonuniform mutual-memory/gyroscopic fields, and zero-memory/refusal/cancellation/
replay; four shared contraction cases also pass. One complete FE/Cx integration
test is authored but unexecuted because the missing frankenscipy/fsci-special
manifest blocks Cargo loading. This native evidence does not establish complete
outer-adapter compilation or broader MR56 qualification.

### Physical source columns on original structural/acoustic timelines

`HeterogeneousHelmholtzFe::prepare_source_loads` retains the original real
quadratic boundary-velocity and distributed volume-injection columns for a
complete pressure-potential pencil. Its exact heterogeneous material/geometry
binding, original source-count caps, byte ceiling and per-projection product
ceiling are admitted before column allocation. Signed and zero-integral
sources retain every original pressure coefficient. Cached row-wise compensated
sums require no geometric quadrature, material query or new spatial state.
This map supplies full-coordinate owners without pressure-release elimination.

`with_sources` binds those columns to the existing exact structural/interior,
explicit-Maxwell and sparse midpoint preparations. Ordinary and hereditary
`step_into` paths reuse the original native workspace plus one pressure-load
buffer and allocate nothing during stepping. Sparse forwarding adds only the
bounded source-column map and O(N_pressure) load scratch; it never constructs
dense spatial M/K/C forms. The original integrators, all pressure/structural
coordinates, reciprocal interface blocks and every retained Maxwell history
remain unchanged.

Boundary amplitudes are actual inward velocities [m/s]; volume amplitudes are
actual volume flows [m³/s]. Exact owners hold them over the declared interval;
the sparse owner receives their midpoint values. The pressure-potential RHS
receives the original projected flow without density, frequency or phasor
scaling. Existing acoustic input work already includes the physical pressure
trajectory integral dotted with the held load, or the native midpoint work,
exactly once. A source may remove stored compression energy. Material/wall heat
and every branch's original work, stored energy and relaxation heat remain
separate from that signed input.

Native state/workspace and branch-array bindings, complete physical balances
and final cancellation admission still control publication. Rejection leaves
the accepted fields, histories and clock unchanged. These are prescribed real
source amplitudes on the original frozen linear geometry/material laws, not a
conversion of a frequency-local impedance into a transient law. Source time
resolution, sparse midpoint convergence and spatial accuracy remain independent
obligations; source preparation does not establish those accuracies.

### Exact full-field cavity evolution with original relaxation memory

`StructuralInterior::prepare_relaxation_exact_time` passes the complete original
M/C/K and explicit Maxwell channels to `fs-phs::mechanical::LinearMechanicalSystem`,
then caches the native `PreparedQuadraticFlow` for held original generalized
forces and projected inward-flow loads. K is relaxed. Original [u,psi], velocity
[v,p], full consistent mass, gyroscopic/nonproportional damping, reciprocal
interface coupling and all native viscous-strain histories survive. Initial
stresses are supplied explicitly, or initialized at zero by `relaxed_state`;
a nonzero history that cannot survive the native coordinate conversion refuses.
Constant acoustic potential and pure dark material histories remain, including
the latter's independent relaxation heat.

Native mechanical stiffness admission uses its own diagonal-congruence scales
and the existing storage PSD tolerance. This retains original acoustic gauge
roundoff without kinetic or Maxwell scales concealing a resolved negative
elastic direction. Exact negative diagonals and coupled zero diagonals refuse;
the original stiffness entries and strict physical-heat admission are retained.

The same exact trajectory supplies original structural/fluid storage, separate
structural and memoryless-wall loss, signed B/B^T interface work, and each
branch's sigma^2/(2*k) storage, integral sigma^2/(k*tau) heat, and separately
integrated work drawn from structure and acoustic field. No heat is inferred
from a deficit. Original structural, acoustic, per-branch memory, aggregate
memory and complete-system balances gate publication. The original input
convention and pressure p=psi_dot are unchanged.

Preparation admits every 2*n+b state and the native doubled augmented bound
2*(3*n+b+1)<=256. Explicit adapter payload/work-score and native cache
byte/multiply budgets precede allocation. Warm `step_into` uses fixed scratch,
cached contractions and no allocation, mass solve or exponential. Source and
branch-array binding, numerical refusal and cancellation preserve the complete
accepted field, histories and clock; no truncation or positive-frequency floor
is introduced. Exactness concerns the frozen semidiscrete held-input model up
to numerical propagation error, not changing materials, continuum accuracy or
output antialiasing.

Seven actual production numerical-kernel regressions pass, covering all three
analytic affine-cavity cubic modes and independently integrated heat; forced
partial fractions; dark histories; mixed and wall channels against the original
physical ODE; zero-memory compression dynamics; late cancellation/invalid
inputs/clock rollback; and native passive/resource refusal. The real mechanical
and exact propagation owners are executed. The additional full FE/Cx zero-memory
and source-binding regression is authored but unexecuted while the canonical
dependency closure remains unavailable. Independent source review found no
blocker; the numerical evidence does not establish outer-adapter compilation.

### Adaptive physical time for retained porous gas/wall compression

`ThermalPorousGasSystem::advance_adaptive` resolves a complete physical stroke
through the existing heat-before-motion owner, including ideal or homogeneous
Peng--Robinson gas, permanent porous compaction, original wall calorics, native
transport and optional retained inertia. Prescribed coordinates move linearly
from their actual accepted positions to the requested endpoint. Other supplied
loads and pressure are held; external gas heat is supplied at constant power.
These loading assumptions are explicit. Birth/remap changes require a separate
admitted stroke boundary.

Each comparison evaluates one full interval and two sequential half intervals
from the same accepted histories. Only the two-half path becomes physical state.
No extrapolated plastic metric, gas state or observer history is published.
Local error includes gas temperature and volume, direct original caloric storage,
wall temperature and specific enthalpy, displacement, velocity, irreversible
porous metric/strain/volume, and native heat/work/numerical-loss channels. Relative
scales use interval changes and do not grow with a formation-enthalpy datum.
All coarse/rejected/fine native heat probes share one bounded global trial count.

Accepted energy increments use compensated accumulation, retain numerical
endpoint-work and kinetic-loss accounts separately, and exclude every rejected
receipt. Whole-window gates bound summed absolute energy residuals, interface
representation loss and direct original-to-final gas storage closure. The final
interval receives the requested external heat remainder. The report keeps compact
stroke totals and the complete last native candidate; storage does not grow with
the number of accepted meshes. A late source, observer, nonlinear, cancellation,
budget or accuracy refusal publishes no partial stroke.

Three focused regression definitions cover convergence toward an independent
continuous two-capacity heat solution, coarse/rejected work budgets and late
atomic refusal, and retained heated compaction with velocity, permanent volume
change and an independently recombined first law. Local full-versus-two-half
control resolves the declared first-order split; it is not a global error bound,
second-order integrator or calibrated moving-wall transport model.

### Adaptive second-order heterogeneous finite-solution diffusion

`FiniteMoistureSolution::advance_layered_sdirk2_adaptive_to` and `advance_layered_radial_sdirk2_adaptive_to` complete a slab, coated sphere, or coated cylinder duration with the original coupled SDIRK stages. Each full and half proposal solves the material storage, source diffusion, interface and exposed-face resistances, and finite solution activity together. The existing temporal controller uses its second-order correction once, charges both implicit stages (including rejected probes), and compares every local content/activity, actual solvent mass/activity, and the L1 water discrepancy weighted by each original cell's dry mass. Only the accepted fine path contributes exterior and intermaterial transfers; first-stage diagnostic reports are not additional physical transfers.

Each original layer retains its own water account, and one fixed original-to-final conservation tolerance applies to the whole requested duration. A work, source, balance, native-solve, or cancellation failure discards the complete private bath/body candidate, including previously accepted intervals. This is local temporal error control for held geometry, temperature, and exposed faces; smooth source and algebraic phase branches are required for the second-order estimate, and finite-rate phases retain the native SDIRK refusal. Three authored regressions cover refinement against independent continuous finite-bath equations for planar/spherical/cylindrical geometry, exact recovery of two native half steps without stage double counting, and work exhaustion after an accepted private prefix. These integration tests were source-reviewed and formatted but were not executed in this environment; independent fixture calculations do not substitute for their native execution.

## Original-layer Kelvin heat in the dynamic thermal shell

`DynamicLayeredThermalShell::with_layer_viscosity` explicitly binds one existing
`PlaneStressViscosity` to each original thermal control volume / mechanical
sublamina. Card identity, physical reference conditions and material axes 1,2
must match that same source layer. Frequency applicability coordinates do not
replace physical conditions; every viscosity must cover the caller-declared
motion bandwidth and the original temperature interval. Runtime temperature
refreshes retain the exact selected claims. The interval law is evaluated at
the mean accepted/trial temperature and frozen on both native mechanical
endpoint owners. No accepted state or source snapshot advances on a failed
thermal, mechanical, material or convergence trial.

Original layer faces, orientations and the laminate integration owner supply
the complete viscous A/B/D, including eccentric membrane/bending coupling.
The original enriched MITC midpoint owner retains physical mass, support
coordinates, independent director motion and all nodal/internal cross terms.
Unchanged operators retain its prepared factorization. The reduced source law
is viscous stress = V times total engineering strain rate. Transverse viscous
resultants are explicitly zero because this source does not characterize them;
elastic shear coefficients, correction factors and thermal eigenstrain rates
do not manufacture a different viscosity law.

Native midpoint strain-rate readback uses mean velocities directly in the
original kinematic evaluator. Each source's nonnegative work contraction is
integrated with the original positive triangle rule and exact two-point
through-thickness rule, then summed over facets into its own original thermal
cell. This physical heat must independently agree with native mechanical
dissipation at a heat-scaled roundoff allowance. It enters the existing heat
solve once, separately from signed reversible exchange, numerical drilling,
kinetic storage and external heating. Mechanical and combined energy accounts
include the same dissipation without treating a residual as heat.

The opt-in path retains the original zero-loss trajectory. Focused regressions
exercise independent scalar Kelvin decay, unequal original-layer heat,
temperature-dependent pinned midpoint coefficients, finite caloric heating
and rejected-trial preservation. This is reduced in-plane Kelvin damping on
fixed reference geometry, not full 3-D dissipation, hereditary memory, a
spectral certificate, or lateral thermal transport.


### Coupled finite forming with sustained tools and current heat paths

With `experimental-solid-contact-heat`, the finite thermoplastic owner advances
against explicitly translating rigid planes through the original sustained
contact solve. Every thermal trial starts at the same accepted plastic metric,
source family, temperature, deformation and optional reference-mass velocity.
The final mechanical state is re-solved at the actual returned caloric
temperature. Native force, complete affine-path, material and independent
contact-work gates remain authoritative; a failed proposal publishes no owner.

`FiniteThermoplasticContactEnergy` retains the native material account alongside
original dead-load, body-contact and prescribed-tool work. The coupled first law
counts tool work once, removes body-contact work from the dead-load port, and
retains the signed normal endpoint defect with the material defect. Only native
material caloric transfer enters thermal cells. Contact quadrature defects and
numerical kinetic losses are never deposited as heat. Picard and IQN-ILS reuse
this same proposal; numerical convergence is not a time-accuracy certificate.

`advance_with_deformed_conduction` prepares complete current normal transport
metrics from the SAME proposed finite solid before each native heat replay.
Actual-temperature mechanics separately supplies final metrics. Every used
area and normal distance must meet an explicit relative discrepancy allowance,
independently of temperature, per-cell heat and whole-energy gates. Original
conductivity/resistance laws, reference mass/geometry, caloric identity, clocks,
scheduled baths and source heat remain retained. The returned observer stores
the metrics actually used, with their discrepancy reported separately.

`advance_with_plane_contact_and_deformed_conduction` composes these capabilities
through one shared mechanical proposal and one shared transport helper. A moving
tool can permanently deform a heated specimen while its current areas and normal
heat paths alter the thermal response. Geometry callbacks must be pure and
justify the normal two-point reduction. Topology remains fixed; arbitrary
skewed-tetrahedron Fourier accuracy, inferred anisotropy, opening interfaces,
friction, curved/rotating tools and remeshing remain outside this composition.
Independent physical-interval, thermal-substep and spatial refinement is needed.

Four synthetic integration regression definitions cover independent diagonal
forming/heating/release and permanent metrics, original consistent inertia and
late observer refusal/retry, current-geometry conduction with an independent
four-step thermal recurrence and geometry-only refusal, and simultaneous moving
tooling/current-area heat. The final test compares actual reaction, temperature,
material heat and geometry through the public combined entry. Focused Rust
formatting/parsing and TOML checks passed; independent scalar calculations are
not Rust execution. Native compilation and these tests remain unverified because
Cargo and the configured remote runners are unavailable in this environment.

### Native second-order heterogeneous diffusion against finite vapor

`FiniteMoistureAtmosphere::{advance_layered_sdirk2, advance_shared_layered_sdirk2, advance_layered_radial_sdirk2, advance_shared_layered_radial_sdirk2}` couple both native layered SDIRK stages to the original finite atmosphere. Slabs, coated spheres and cylinders retain each material's original nonlinear sorption and diffusivity, exact native geometric storage/resistance, and explicit exterior/interface transfer laws. The first finite-gas root solves humidity at the implicit first-stage time; its converged physical loads are captured once. The final root starts from original gas and material storage with the retained weighted first load, rather than advancing or debiting the temporary first-stage atmosphere again.

One generic participant orchestration retains the original native slab/radial owners, canonical stable IDs, original clocks, aggregate cell budget and source/environment admission. Intermediate complete bodies are released before final-root candidates are allocated. The endpoint report carries the complete weighted exterior/interface transfers; the first report is diagnostic. Original-to-final vapor, carrier and total mass checks and one gas-energy increment remain in the finite-atmosphere owner. Source, stage-range, work and cancellation failure returns no accepted gas/body candidate. Geometry, gas volume and temperature remain held; this path does not infer condensation, thermal feedback or a temporal accuracy certificate, and second-order behavior assumes smooth admitted source branches.

Six authored native regressions cover independently assembled finite-gas storage graphs, both physical stages and their weighted flows, continuous second-order refinement, uptake/drying with continuation, unequal curved capacities and finite surface/interface resistance, original source surface receipts, opposing shared exchange in canonical order, original-to-final gas energy, and late-stage/work refusals. The native integration tests remain unexecuted under the existing dependency/resource block. Their separate Rust analytical reference helpers compiled and ran; the finite-vapor reference refinement ratios are approximately 3.846 and 4.015, which is fixture validation rather than native execution.

### Sorption energy on native layered finite-air SDIRK endpoints

`advance_layered_sdirk2_with_sorption_energy` and `advance_layered_radial_sdirk2_with_sorption_energy` compose the original two-stage finite-atmosphere diffusion with the same layered caloric owner. Each material's pinned loading-dependent isosteric heat is integrated over its actual original-to-final cell moisture and dry mass, including transfer between materials at zero net vapor uptake. The temporary first-stage transport report is diagnostic: it adds no second heat, gas-energy, or source-segment debit. One original vapor datum, one whole-interval gas-energy increment, one aggregate heat-source integration budget and one closed-system energy tolerance apply. Source, budget, energy or cancellation refusal returns neither native endpoint. The existing backward-Euler candidate type remains source-compatible through a default report parameter.

The result remains held-temperature external thermostat duty, with no inferred spatial heat deposition or thermal feedback. Three authored regressions cover the independent closed two-store SDIRK stability solution in both transfer directions, loading-dependent caloric primitives on actual planar and unequal radial cells, and late caloric source/work refusal after both native transport stages. These tests were source-reviewed and formatted, not executed here; the independent sealed-case fixture gives a signed transfer of approximately 8.95486053864e-4 kg and thermostat heat of -895.486053864 J in its forward direction.

### Adaptive second-order heterogeneous finite-air durations

`FiniteMoistureAtmosphere::advance_layered_sdirk2_adaptive_to` and `advance_layered_radial_sdirk2_adaptive_to` complete layered slab, coated-sphere and coated-cylinder durations using the original coupled native SDIRK stages and existing temporal controller. Local material content/activity, actual gas vapor mass/RH, and the L1 water discrepancy weighted by each original dry cell mass jointly govern acceptance. The second-order fine-path correction is applied once; both stages and all rejected probes consume the shared duration work allowance. First-stage diagnostics are not additional water transfers.

Only accepted fine endpoints accumulate exterior and paired intermaterial exchange. Each original layer and the complete gas/body system retain one original-to-final water tolerance, with the existing final carrier and total-mass admission. A late work, source, stage, conservation or cancellation refusal discards the full private duration. Temperature, free gas volume, geometry and exposed faces remain held; local error control does not add condensation, heat feedback or a global accuracy certificate. The generic atmosphere-report forwarding also preserves existing backward-Euler consumers.

Three authored regression groups check refinement against independent continuous finite-air storage equations in all three geometries, exact agreement with two native half steps and accepted-only transfer accounts, and work exhaustion after an accepted private prefix. All four changed Rust files passed the pinned formatter check. Native tests remain unexecuted under the known dependency/resource build block; independent graph/controller fixture checks are not native test execution.

### Relaxed finite-solid release with original reference density

The opt-in finite-solid liquid release binds the actual accepted plain hyperelastic
body and complete original thermal timeline at the same clock. Original reference
positions, per-element or retained P1 nodal densities, current motion and the identical
caloric object enter the existing reference-density/native thermal ingress. The
original consistent mass operator independently supplies M, first moment, P, L and
datum-centered K. Different reference densities and volume Jacobians can enter one
native liquid only when their current physical densities agree within the explicit
kg/m3 bound. Actual native occupied fill volume and mass remain separately audited.

Each thermal cell maps to one original tetrahedron with its retained reference
mass/volume and fully liquid caloric state at the original mechanical temperature.
Shared vertices require equal incident-cell h; no P0-to-P1 averaging is invented.
Original external conditions and density claims stay pinned. Detached-component
mass, momentum, angular momentum and enthalpy cannot cancel between pieces.

Native Cauchy stress, signed/absolute elastic storage and accepted M*a must fit
explicit near-zero gates. The signed mechanical discrepancy is liquid K minus
original consistent-mass K minus original elastic storage. Absolute source kinetic
identity error, per-component projection error and absolute elastic storage share
the projection budget; none becomes heat. Supports, pressure/gravity, original
nonzero loads, receiving gravity, thermal ports/links and future events refuse.

This zero-duration representation change does not supply melting kinetics,
stress relaxation or plastic/damage/viscoelastic/contact-history discharge. Both
borrowed sources remain unchanged on any outcome; accepting the native thermal
candidate replaces both source representations at their existing clock. Native
regression fixtures cover unequal reference densities/Jacobians, analytic cuboid
inventories, actual fill geometry, continuing thermal liquid motion and refusal/
retry for source identity, clock, heat inventory, retained owners and work budgets.

## Conservative mixing across original liquid density references

`liquid_transfer::transfer_interface_mass` accepts distinct temperature/source
bindings and physical density/dx/dt scalings when the exact original caloric
object, immutable material-card identity, complete optional property
requirements and selected claims agree. The named temperature value may
differ; typed axes, nonthermal coordinates, viscosity and surface conventions
remain fixed. Equal sampled coefficients or material names are insufficient.

`transfer_interface_mass_streams` withdraws one actual packet from each
bounded original donor and mixes all streams once at one coincident receiver
interface and world point. Their native withdrawals supply physical kg,
original specific enthalpy and Guo-corrected velocity; rotations retain the
common inertial-frame momentum account. The existing multistream receiver
retains velocity dispersion and deposits the physical complete-mixing loss
as heat exactly once. Collective capacity is checked before donor proposals.
Actual mass, world momentum and enthalpy-plus-kinetic defects pass their
bounds before a final cancellation poll and infallible joint publication.
No clock advances.

Incoming pure-liquid calorics need not have the receiver's reference density.
The accepted receiver must still satisfy its own original storage-density
admission. Balanced hot/cold feeds and source-density plateaus can therefore
mix different-density streams; an unmodeled endpoint expansion refuses
atomically. No population-density reset, density averaging or replacement
composition/EOS law repairs an incompatible result.

The reference-volume report divides actual withdrawn/received masses by
their own physical reference densities. It exposes their signed volume
difference separately from low-Mach pressure geometry. This instantaneous
coincident-port model adds no resolved connector flow, pressure work,
capillary energy or general total-energy theorem.


## Reciprocal liquid/rigid viscous momentum source

`liquid_momentum::apply_rigid_drag_native` and its material-bound adapter exchange momentum between a compact stencil of original 3-D free-surface liquid cells and the original finite `RigidBodyState`. The declared homogeneous inverse-time rate multiplies each cell’s actual retained liquid mass, including partial interface fill. The six-coordinate kinetic metric includes native mass, arbitrary principal-axis orientation and full rigid inertia. A bounded symmetric eigensolve integrates both rigid-response modes and remaining liquid slip modes exactly for the frozen geometry; large source durations relax slip without explicit-step reversal. Work and storage use SI units and the same actual momenta.

This is an operator-split viscous **source**: grid positions, rigid pose, source material/scaling, force, fill, collision/transport clocks and native free-rotation integration are held fixed. The source duration is not a drift/flow step. The authored rate is not inferred viscosity, hydrodynamic added mass, contact, excluded volume or an impermeable boundary. Reported nonnegative drag dissipation is mechanical kinetic loss; no heat owner is implicitly changed.

All selected physical velocities, true wet masses and body-point speeds are admitted before solving. Actual population readback and native mechanical impulse receipts supply linear momentum, moment about the fixed native COM, signed kinetic change and a relative-velocity work-loss account. Local readback defects, global conservation defects and independently integrated body-response defects cannot cancel. The complete fluid candidate and finite-body reaction are prepared before either publishes. Refusal or cancellation preserves both owners. Memory is proportional to selected cells, with one fixed six-dimensional solve and no dense domain impulse field or whole-grid clone.

The material adapter obtains dx, velocity and mass units from the original immutable liquid binding and enforces its scenario speed before staging. Source density and scaling cannot be independently replaced by the caller.

Fourteen focused native tests passed in this session: six exact-map cases use the production kinetic-metric solve with original Jacobi eigenvectors and double-double arithmetic; five coupling cases use production staging, original sparse population storage, the shared impulse kernel and original full-inertia MBD state/work/impulse methods; three compact-stencil cases check actual populations, readback and refusal. These cover analytical finite two-mass relaxation, independent vector-force/torque RK4, mixed retained masses, rotated inertia, semigroup and Galilean behavior, stiff relaxation, physical balances and late rollback. The selected-source harness does not execute unrelated collision/streaming or the complete material graph. Two original-material regressions are included but were not executed under the current full-Cargo resource limit.

## Sharp-edged tone-hole orifice on the DG route (`dg_terminal_orifice PORT_ID C_VC`)

`DgTerminalChannelKind::Orifice` is a terminal channel like `Closed`, `Slit`, `Circular` and `Annular`, declared in the
`dg_terminal_channels` block: `dg_terminal_orifice PORT_ID C_VC`. It puts the open tone hole's quasi-steady sharp-edged law
`p_duct - p_exterior = K |Q| Q`, `K = rho / (2 C_vc^2 A^2)`, in series with the chimney's exterior load. `A` is the DG chimney's actual
outlet area and `rho` the terminal's own gas (it follows a gas refresh); the law has the form of the live reed card's
`open_hole_orifice` (`fs_duct::open_hole_orifice_quadratic_coeff`). `C_VC` is the card's number, finite and in (0, 2]; the posted
oboe pages use `fs_duct::SHARP_EDGE_VENA_CONTRACTA` = 0.61 = pi/(pi + 2), which that constant's doc (fb6ff3ced) names the steady
thin-plate vena-contracta value (Kirchhoff free-streamline contraction of a plane jet), not a tone-hole measurement. The printed anchors
quoted there, in this same `K |Q| Q` form, are C_vc = 1 (Atig, Dalmont and Gilbert 2004, thin-walled unflanged open end), 1.18 (flanged)
and 1.096 (the side-hole term of Dalmont et al. 2002 as quoted by Szwarcberg et al. 2024): a `K` 2.7 to 3.2 times smaller than the 0.61
gives. A ruling on that constant therefore moves both routes: put the same value on the card.

The row needs a `dg_side_chimney` on its port (the axial outlet and flush side ports are refused) and no driven exterior mouth. It
takes no transverse cell. A port with no row, and a card with no row, are unchanged bit for bit.

Scheme. Each mechanical sample solves the linear terminal group with the series resistance `Z = K |Q_prev|`, `Q_prev` the previous
accepted sample's outward flow (zero before the first, so the first sample has no resistance). `Z >= 0` makes the group passive by
construction: the resistance returns no work. The coupling (`TerminalCoupling`) is rebuilt for the NEXT sample inside the fallible
phase of `preview_step` and published by `accept`, so a dropped or refused trial changes nothing and `port_response` and
`preview_step` always read the same accepted resistance. The channel frame reports the drop the sample was solved with
(`K |Q_prev| Q`), which is what the existing same-step pressure and work audits see; its work `Z Q^2 dt` is entered in the load's
ledger as port work and dissipation, with no storage.

Lag. The law is `K |Q| Q`; the solved drop uses `|Q_prev|`. `DgTerminalChannelLoad::orifice_report` keeps, over the accepted
samples, the largest `|K (|Q| - |Q_prev|) Q|` [Pa], `sum |K (|Q| - |Q_prev|) Q| |Q| dt` and `sum K |Q|^3 dt` [J] (their ratio is the share of
the law's dissipation the lag misattributes) and the work taken. The lag is `O(dt)`: tests settle a constant flow to the law
and show the share falling below 0.7 times its value when the step is halved on the same 1 kHz drive. The sidecar carries the law and
its coefficients (`sharp-edged-orifice`); the end-of-render lag account is the accessor, not a sidecar field.

Evidence (G0/G1): `dg_port::channel::orifice_tests` (the first sample is the plain group's bit for bit; every sample's solved drop
is `K |Q_prev| Q` and equals the observed duct-minus-exterior drop; the work is nonnegative and equals `Z Q^2 dt`; the ledger drift
stays under 1e-9; a settled flow has no lag; halving the step cuts the lag share; refusals; a dropped trial restores everything),
`dg_port::input::tests::terminal_channels` (the row, its refusals), and the existing terminal-channel tests unchanged.

## Modal aperture flow balance at a near-closed gap

`ModalAperture` certifies the flow balance `jet + swept + body - wave` of every step against the card's relative
tolerance times `flow_scale = |jet| + |swept| + |body| + |wave|`. `wave = (outgoing - incoming) / z` is a difference of
pressures over the bore impedance, so its arithmetic carries the rounding of those pressures, `eps max|P| / z` with `P` over
the upstream, no-drop, junction-drop, incoming and load pressures, however small the net flows are. As the gap closes the net
flows go to zero while that rounding does not, and the declared budget `relative x flow_scale` falls below it.

The flow gate therefore takes the default-API certificate: its rounding floor is taken over `flow_scale + max|P| / z` and,
where that floor is above the declared ceiling, the allowance is taken and recorded (`floor_raised`, `effective_budget` equal
to the floor). Where the floor is under the declared ceiling the certificate is exactly the explicit one (same effective
budget, `floor_raised` false), so every step that was admitted before is admitted with the same ceiling. The work gate is
unchanged and still an explicit declaration. A residual above the raised floor still refuses (`modal/refusal_tests.rs`:
1e-15 at the failing scale, six orders above the rounding). The gate admits or refuses; it changes no number of the step.
Measured on the refusal that motivated it (a 2.65 s D4 render on the typed DG reed, 1.047 s in): residual 4.2e-21 m^3/s on a
flow scale of 6.2e-11, `eps x 6178 Pa / 1.7e8 Pa s/m^3` = 8e-21 m^3/s.

## Positive exact blade damping heat

The shared `ExactBlade` used by scalar, paired and diagonal modal apertures
retains physical viscous and Maxwell loss through the existing
`fs-phs::linear_exact::PreparedSourceLoss` owner. Its source rows are
`sqrt(c) v` and `(c_j x-z_j)/sqrt(tau_j)` on the original held-force trajectory.
The unit-interval generator is the same `A dt` as the motion map, with source
rows multiplied by `sqrt(dt)`. No auxiliary inverse physical exponential or
expanded energy difference supplies the heat. Positive Gauss/QR factors retain
nonnegative loss, exact zero loss at resolved static equilibrium, and small
vibration heat about a large held-load displacement.

Preparation retains the existing 256-coordinate doubled-block capacity and
explicit byte/work bounds. Warm heat evaluation uses bounded stack scratch;
it allocates no vector, advances no history, and changes no motion or contact
solve. The material-source clock and applicability restrictions are unchanged.
The existing Maxwell ODE regression includes stiff relaxation rungs; independent
overdamped poles are checked through the live paired aperture, and an analytic
loaded-vibration regression distinguishes tiny physical heat from its much
larger equilibrium spring storage. These are synthetic G1/G3 checks of the
represented equations, not measured material or instrument validation.

## Per-opening quadrature of a circular opening array

An array of coplanar circular openings (`BaffledPressure::circular_openings*`) used to give every opening the receiver's one
declared (rings, angular points). The surface rule (cell phase `2 pi f (R/sqrt(rings) + 2 pi R/angles)/c <= pi/4` over the
declared band) is per opening, so one large opening could refuse a card the small ones resolve at the same grid. Now an opening
the declared grid resolves keeps it exactly (the local receiver is the one the unrefined build used, bit for bit, so every
downstream arithmetic step is untouched); an opening it does not resolve takes the grid with rings' in rings..=128 and
angles' in angles..=512 that resolves the band with the fewest points (ties to fewer rings); where none does the declared grid
stays and `circular_surface` refuses with its own text. The shared point budget counts each opening's own grid and refuses
above 120 000 with the same text (with equal grids it is the former `openings x rings x angles` rule, which is still checked on
the declared grid first). `BaffledPressure::resolved_circular_grids(openings, receiver, medium)` returns the grids in opening
order. The grid rule's phase test copies `circular_surface`'s expression and is checked against it on a sweep.
On the oboe DG card (the 19 mm bell and its three open chimneys, receiver on axis at 0.5 m, 20 kHz band, declared 32 x 128)
the bell refines to (128, 257), the minimum by exhaustive search, and the chimneys keep (32, 128): 45 184 points, with 6 287
distinct source-to-receiver paths for the gas-absorption observer against 43 202 for the shared 128 x 224 grid that resolves
the same band (computed from the surface rule; the card's `observer_gas` budget is 16 384 paths). The 20 kHz pilot that the
shared grid refused at the gesture phase renders on this rule (rc 0), and its wav equals the 9.6 kHz-observed pilot's byte for
byte (both quadratures are converged where the note has energy, and its 9.6 to 19 kHz band is 67 dB under the total, about
one 16-bit step).
This is a quadrature resolution of a prescribed uniform-velocity piston in an infinite rigid baffle; it adds no cross-mode
physics of the real mouth above the bell's first cutoff (about 5.3 kHz for 19 mm).

## Finite battery work and material self-heating

`circuit_enthalpy::battery::BatteryThermalNetwork` couples the original finite
battery's chemical charge, polarization RC pair and bus capacitor to the existing
material/enthalpy residual. Series and polarization heat have separate explicit
receivers; a common receiver is allowed as declared thermal homogenization.
Native discrete-gradient OCV work integrates across crossed SOC knots, and the
original two nonnegative branch heat integrals enter storage exactly once.

Resistances are explicitly fixed or selected positive scalar-ohm properties from
the receivers' original material cards. Selected claims refresh at actual trial
T/phase while retaining source identity and fixed external context. Each trial
starts from the same accepted charges. Battery and thermal clocks, finite SOC
support, thermal events, numerical budgets and cancellation admit both histories
together or publish neither. OCV and both capacitances remain fixed: no OCV(T),
reversible entropic heat, aging, chemistry fit or protection control is inferred.
Coupled temporal accuracy is generally first order and still requires refinement.
Focused `tests/battery_enthalpy.rs` cases cover crossed-knot work, independent
calorimetry, R(T) feedback, charging/discharging, shared receivers and refusals.

## Finite battery, converter, regulated load and heat

`circuit_enthalpy::converter::battery::BatteryFedConverterThermalNetwork` advances
the original six electrical coordinates and finite material enthalpy together.
Battery chemical/polarization/bus storage, converter flux/input/output storage,
the positive physical lead resistance, and optional current-limited/foldback
regulated demand share the native electrical solve at each thermal trial.
Modulation and demand may change between accepted intervals without resetting
any charge or flux. The same path supports reverse current and regeneration.

All eight physical heat channels have explicit receivers: the converter's five
conduction channels, battery series/polarization losses, and lead Joule heat.
Shared receivers sum those joules before one deposition. Exported load work
remains a separate electrical boundary and is never deposited as converter heat.
Independent battery, lead, converter and whole-system balances retain numerical
closure defects separately. Every trial resolves the original material sources
at its actual T/phase and restarts the same accepted electrical state; identical
resolved resistance tuples can reuse one bounded cached electrical result.

This is averaged conversion with fixed OCV/C/L and endpoint resistance feedback,
not an inferred PWM, entropic battery heat, aging or protection-control model.
Source support, exhaustion, thermal events, budgets and cancellation refuse
without publishing either state. `tests/converter_enthalpy/battery.rs` covers an
independently calculated steady power path, eight-channel calorimetry, nonlinear
load export, temperature feedback and transactional refusal/retry.

`step_with_load_heat_candidate` can explicitly bind that regulated load to a
finite material receiver when its physical model is wholly dissipative. The
accepted native load work, including voltage foldback and current limiting,
heats that cell once. Requested or curtailed power never supplies heat. The
same work cancels as an internal transfer in the combined electrical/caloric
balance; the original eight resistance channels and native electrical audits
remain separately visible. The accepted source sample records the load's cell
and energy. A separate passive load body or a receiver shared with a resistance
channel uses the same material inventory, endpoint residual and source-domain
admission, so thermal resistance feedback remains coupled to the electrical
trajectory. An omitted receiver preserves exported-load semantics. No useful
work, inferred efficiency, extra load storage or controller bandwidth is claimed.
`tests/converter_enthalpy/battery/load_heat.rs` checks independent steady-state
chemical depletion/calorimetry, foldback and current limits, shared-body thermal
feedback, disabled-load history, invalid mapping, caloric-domain refusal and
late-cancellation rollback/retry.

## Source-bound semiconductor commutation and material heat

`converter::commutation::CommutatedConverterThermalNetwork` joins the original
three-state finite-resistance commutation owner to the shared finite-cell
enthalpy residual. Four explicit device packages retain their original material
cards, selected claims, caloric bodies and seven typed source keys: channel ON
and OFF resistance, body-diode forward threshold and differential resistance,
maximum reverse voltage, maximum channel current and maximum diode current.
The original source owner refreshes the complete bundle at each proposed
junction temperature and declared phase coordinate. Unowned current, voltage,
gate and frequency dependent source laws are refused; no nominal parameter,
leakage, limit or device dataset is invented.

Three ordinary channels retain the existing scalar/bulk material resistance
adapters for the winding and both supply resistors. Actual winding, four
channel/leakage, four body-diode and two source-resistor heats remain separately
visible in the native cycle report. Each channel and its body diode heat the
declared package body once; shared receivers accumulate in double-double
arithmetic before deposition. Optional wholly dissipative load receivers absorb
only actual delivered work. Regeneration remains signed ideal-source work.
Independent device/converter, each source/load path and complete storage/work
balances precede internal routing.

Every proposal restarts the same accepted flux and charges, gate intervals and
source/load commands. An exact 31-word coefficient key includes all device
thresholds/slopes/limits and the three ordinary resistances; unchanged resistance
alone cannot reuse a trajectory when diode voltage or applicability changes.
Source and caloric admission still run on cache hits, and the full electrical
subdivision count is charged to the shared trial budget. The complete electrical
cycle must reach the same represented endpoint as the thermal solve. Refusal or
cancellation publishes neither history.

This is the native quasi-static conduction model with finite authored OFF
resistance and endpoint material feedback. Temporal refinement is required;
sampled device limits do not certify unobserved waveform peaks or a full
semiconductor safe operating area. Junction/output capacitance, gate dynamics,
reverse recovery, avalanche and transition losses are outside this lane.
`tests/commutated_converter_enthalpy.rs` defines independent steady KCL/KVL
calorimetry with regeneration and shared load heat, diode-threshold feedback
with every resistance held constant, and source-domain/quantity/limit,
clock/budget/cancellation refusal. Its material curves are manufactured fixtures,
not measured device qualification.

### Finite battery-to-shaft commutation with material feedback

`converter::commutation::battery_motor` connects the original nine-state
finite battery/lead/converter/motor commutation owner to the shared finite-cell
enthalpy driver. `CommutatedBatteryMotorThermalNetwork` retains five original
resistance bindings, four existing seven-role device bundles and an explicit
bearing receiver. `CommutatedBatteryMotorHeatNetwork` uses the same resolved
path through the complete original `HeatAssembly`, including finite coolant,
contact, radiation and reservoir transport.

Every thermal proposal refreshes actual resistor/device material coefficients,
then restarts the native gate program from the accepted chemical/polarization
charge, three bus charges, two fluxes and both absolute shaft momenta. All
fourteen native heats remain distinct in the report. Channel and diode heat
enter their declared package once; shared receivers accumulate in double-double
arithmetic. Optional wholly dissipative auxiliary receivers receive actual
delivered work only, with one internal-transfer cancellation in the combined
external account.

Independent battery/lead/converter/output/magnetic/shaft/native/full energy
gates precede deposition. A 33-value exact coefficient key includes five ordinary
resistances and all device coefficients/limits; fresh source admission still
runs on cache hits and every electrical subdivision charges the trial budget.
The represented electrical/thermal endpoint must agree. Source/operating
support, unowned electrical graph work, clock/budget or cancellation refusal
publishes no physical history.

Fixed OCV/storage, motor reciprocity and bearing drag, sampled device support
and endpoint-feedback temporal accuracy retain the native limitations. No
entropic battery heat, transition/recovery loss or unobserved-peak/safe-area
qualification is inferred. Four focused regressions define independent
motoring/regenerative fourteen-heat calorimetry, threshold-only feedback into
motion, finite-liquid-coolant exchange and transactional refusal. Sources are
manufactured fixtures, not measured device qualification.


### Dynamic device output charge with fixed material capacitance

`converter::commutation::battery_motor::capacitance` connects the original
eleven-state finite battery/converter/motor owner to the existing material and
HeatAssembly enthalpy path. Both switch-node charges remain continuous across
gate and thermal boundaries. Incoming bus/node coordinates are total electrode
charges. Binding requires four explicit authored capacitances and a bitwise
match to original positive farad claims on the corresponding semiconductor
packages; replacing capacitance at fixed incoming charge is refused.

`DeviceCapacitanceBinding` and `DeviceCapacitanceMaterial` retain the original
card, claim, property key, package and held external context. Every actual
caloric proposal refreshes source validity, even when the electrical trajectory
cache is reused. Only original zero-axis `ConstantWithinValidity` laws with
owned static applicability are admitted. A sampled temperature- or
voltage-dependent capacitance is refused. All storage coefficients remain fixed;
variable charge laws and their reversible parameter/entropy work require a
separate coupled model.

The actual eleven-state receipt retains four signed device-capacitor energy
changes and the original bus/inductor storage. Their independently summed
partition precedes the seven-constituent and fourteen-loss physical balances.
Only those fourteen dissipative terms and optionally routed actual auxiliary
work heat the finite bodies. Stored capacitor energy is never deposited as
package heat or reconstructed from a balance remainder.

The nine-state public network/candidate/sample interfaces remain available.
Both lanes share the same source-refresh and enthalpy loop, exact coefficient
cache, full subdivision budgets, represented endpoint clocks and original
HeatAssembly/radiation closure. The capacitor lane has explicit eleven-state
accepted histories, source bundles and native receipts; no synthetic nine-state
receipt or sequential electrical/thermal trajectory replaces them.

Four regression definitions cover steady motoring/regeneration calorimetry,
dynamic capacitor energy with hot-diode feedback and gate continuation,
source/charge/package/domain refusal, and two finite thermal inventories.
Fixture values are synthetic. These new coupled regressions were not compiled
or executed in this delivery; the dependency set was incomplete and the runtime
became unavailable. Gate-drive energy, reverse recovery, avalanche, nonlinear
Coss and unresolved-peak qualification remain outside this constant-C lane.

### Full device charge curves through the material and heat owners

`converter::commutation::battery_motor::charge` connects four original Q(V)
laws to the same eleven-state finite battery/converter/motor owner and shared
finite-cell/`HeatAssembly` thermal path. Curves follow input high/low then
output high/low drain-to-source signs. Binding compares every authored
charge/voltage knot and its full support against the material law before
accepting total electrode charges. Constant-capacitance coordinates and a
sampled differential capacitance cannot replace that law.

Each source declares a coulomb property, an owned typed volt axis and finite
support on both voltage polarities. The complete selected `PropertyLineSlice`
and claim remain retained. An exactly uncharged point inside an original
affine segment may be inserted without changing its law or integral.
Temperature/phase applicability refreshes at each actual caloric proposal;
only voltage may change constitutive charge. A temperature-dependent charge
law, changed support/claim/held context or changed curve refuses. Once the
fresh complete path matches the retained law, its original admitted slope
and energy tables are shared rather than rebuilt.

Native device energy is the original integral of V dQ. Its four signed
changes remain separate in the actual receipt and contribute to converter
storage, while the existing fourteen physical losses alone heat their declared
bodies. Optional auxiliary routing uses actual delivered work. The same
independent energy gates, coefficient cache, full electrical trial charge,
represented clock and original heat/radiation transport serve nonlinear and
constant laws. Exact native chord traversal has an explicit numerical budget
recorded with each proposal; exhaustion refuses instead of resampling or
weakening the law.

Four regression definitions exercise independent motoring/regenerative heat,
integral energy distinct from QV/2 with all eleven discrete balance equations
and hot-diode continuation, full-law/source/domain/temperature refusals and
affine zero-reference insertion, and original two-inventory heat transport.
The new coupled source and tests pass a syntax parse with the pinned Rust
formatter. Coupled typechecking and execution remain unverified because the
available native dependency set is incomplete. Endpoint feedback requires
time refinement. Synthetic curves do not qualify measured devices, gate
dynamics, recovery, hysteresis, avalanche, temperature-dependent reversible
field energy or unresolved operating peaks.

### Finite gate charge, Miller feedback and separately received driver heat

`converter::commutation::battery_motor::gate` connects the original
fifteen-state gate-controlled finite battery/converter/motor trajectory to the
same material and finite-cell/`HeatAssembly` enthalpy owners. Its nine converter
coordinates retain four gate charges in addition to flux, two bus charges and
two node charges grouped with the high-side gate electrodes. All accepted
coordinates continue through commands and thermal boundaries.

`GateMaterialBindings<B>` requires four complete intrinsic drain-source Q(V)
sources, four positive constant Cgs sources, four nonnegative constant Cgd
sources, four complete channel sources and eight original on/off
driver-resistance bindings. Its historical `channel_conductance` field defaults
to positive G(Vgs,T) lines and also accepts `GateCurrentSurfaceBinding` for
complete original channel-only Iref(Vds,Vgs) surfaces. Both material forms
implement the sealed `GateChannelBinding`/`GateChannelMaterial` interface and
share the same preparation, electrical solve, thermal iteration and accounting.
Charge, capacitance and channel response belong to the corresponding
semiconductor package. Each on/off resistor pair names one explicit
driver-heat receiver; shared receivers sum their actual heats. Gate-clamped
Coss is not silently reinterpreted as intrinsic Qds.

A current binding owns two distinct existing typed voltage coordinates, an
ampere property and a finite operating rectangle containing both drain signs.
It retains the complete original claim and `PreparedSurfaceResponse`, including
its original samples and triangles. The original full-rectangle envelope
query proves continuous coverage including the interior. No conductance-grid
conversion, I/V division, nominal drain voltage, resampling or extrapolation
is used. Original affine cells supply current and both analytical partials;
source-owned cutoff, saturation and passive negative differential slope are
admitted. Same-sign drain cells, zero current on the authored zero-drain mesh
and correctly signed vertex currents establish passivity over each original
cell. Cells crossing zero drain bias and isolated operating samples without
a tangent refuse. The channel excludes body-diode current, which remains with
its original separate owner.

An optional `GateCurrentThermalScale` supplies an explicitly authored positive
dimensionless s(T,phase), giving I=s*Iref and scaling both voltage partials.
Without it the admitted surface remains fixed. The scale may depend only on
owned thermal coordinates and permitted fixed external conditions; evolving
voltage response or voltage applicability cannot be frozen at a reference
point. Both the complete current claim and any scale must qualify the actual
thermal proposal. This separable law does not infer a three-coordinate
response from two-coordinate data.

Full typed voltage support, material claims, package, held context and caloric
applicability remain retained. Current/scale source receipts and full-rectangle
support refresh before numerical reuse, as do legacy conductance and driver
resistance sources. Changed source authority, unowned bias/frequency,
reinterpreted incoming charge and incomplete receivers refuse. Intrinsic Qds,
Cgs and Cgd remain unchanged along accepted charge history; temperature-dependent
reversible storage requires its own caloric model and is refused here.

Every proposal uses the same native gate-controlled solve and original
electrical/thermal loop. The numerical cache compares the complete current
response axes, samples, triangulation, voltage ranges and exact scale, or the
complete resolved legacy G(Vgs) curves, together with all driver coefficients
and the existing 33-word drive key. Retained immutable surface data is shared;
allocation identity alone never proves coefficient equality. Fresh source
admission and full subdivision-budget charging precede reuse.

Converter storage includes the separate intrinsic drain, Cgs and Cgd energy
changes. The same represented channel current supplies electrical flow and
Vds*I device heat. Four signed gate-supply works and four nonnegative driver
heats extend the original fourteen-loss physical balance; driver recovery
does not become negative heat or an inferred battery transfer. The original
finite coolant, contact and radiation owners transport the deposited heat.

Focused G0/G1 fixtures specify source-owned current saturation/cutoff and signed
reverse current, original-cell derivatives, explicit thermal scaling,
full-rectangle support and source/receiver/continuation refusals, in addition
to the existing finite-RC/temperature-feedback and Miller energy-account cases.
The coupled extension and fixtures passed pinned Rust syntax formatting and
source/API review. Full coupled typechecking and execution remain unverified
because the complete native dependency graph is unavailable. Endpoint thermal
feedback requires time refinement. Manufactured surfaces and thermal scales
do not qualify measured devices or unresolved waveform peaks; reverse
recovery, avalanche and finite-energy gate-rail dynamics remain separate gaps.

## Finite source impedance, regulated demand and material heat

`converter::voltage::power_load` adds `step_power_loaded_candidate` to the
existing seven-channel voltage-driven thermal network. Both held ideal voltage
sources, their physical series resistances, the original converter flux/charges
and both optional regulated loads share the native nonlinear electrical solve.
Every caloric trial resolves the same seven material resistance sources at its
actual endpoint; source sag and negative incremental load conductance therefore
respond to that same thermal state. Equal resistance tuples reuse the original
electrical result without bypassing material or caloric admission.

Each load may remain exported or explicitly heat one finite material receiver
as a wholly dissipative device. The actual accepted load work, including current
limiting and voltage foldback, is deposited once; curtailed/requested demand
never supplies heat. Shared receivers sum resistance and load energy before
deposition. The native converter, each source/resistor/load path and complete
electrical balance are checked before internal routing. Routed work cancels in
the combined external account; signed ideal-source regeneration remains visible.
The source sample separately records each load's receiver and actual joules.

These are ideal supplies with finite impedance, not finite-energy reservoirs.
The new path uses the original discrete-gradient effort and midpoint work,
with endpoint material feedback requiring temporal refinement. It does not
change the existing exact linear voltage-driven/PWM paths or infer efficiency,
source capacity, C/L parameter work, or semiconductor switching losses.
Receiver mapping, original material support, clocks, numerical budgets and
cancellation admit both histories together. Five focused regression definitions
cover independent two-bus calorimetry, shared/partial routing and turn-off,
source-resistance feedback, regenerative RC work, actual curtailed delivery,
unsupported caloric states and transactional cancellation/retry.

### Regulated demand through the original PWM and material timeline

`step_power_loaded_switching_candidate` advances every prescribed gate interval
with its own held source/load command through
`PreparedPowerLoadedSwitchCycle`. Explicit subdivisions resolve the same native
three-state nonlinear converter; they do not average gate/load correlation or
reset flux and charge at an edge. The complete authored cycle must reach the
same represented endpoint as the thermal interval, including at large absolute
timestamps. Its full native integration count is charged to the shared
circuit-trial budget on every proposal, including cache reuse.

The averaged and switched paths share one material, receiver and energy-account
adapter. Each proposal resolves the original seven R(T,phase) channels, then
restarts the electrical cycle from the accepted state. Actual substep
squared-current work supplies all seven resistor heats; squared mean current
cannot replace ripple heating. Each load's accepted work may enter its declared
finite receiver through the same compensated internal-transfer account.
Source/resistor/load, converter and whole-cycle balances remain independent.
The original native step or complete-cycle report is retained without
reconstructing physical heat from an energy residual.

All seven material coefficients use the proposed caloric endpoint over a cycle;
refining the thermal interval and electrical subdivisions is necessary to assess
time accuracy. The source programs, storage parameters and prescribed ideal
complementary gates stay fixed within each trial. This extension supplies no
semiconductor transition, gate-charge, dead-time or recovery law. Four regression
definitions cover an independent switched capacitor/load quadratic and exact
calorimetry, temperature-dependent ripple and shared load heat, bounded work and
transactional retry, and distinct represented endpoints at large clock epochs.

## Native structural vibration from harmonic induction forces

`spatial_induction_vibration::SpatialInductionVibration` transfers the original
harmonic induction response's consistent P1 conduction-current Lorentz forces
to an original `TetLinearElasticProblem`. Explicit one-to-one mechanical-node
and tetrahedron maps must match the original EM coordinates and connectivity;
permutations are admitted, while interpolation, changed geometry and duplicate
entities refuse. Only mapped conducting cells contribute force. Nonconducting
structural support cells may share the mesh. Original mechanical material
states, full consistent mass/stiffness and homogeneous stationary constraints
remain owned by `fs-solid`; constrained forces remain observable and do no work.

The existing `ExactHarmonicSolid` advances simultaneous DC and twice-frequency
forces, preserving mean load, peak phasor phase, startup transients and resonant
growth. With `exp(-i omega t)`, `Re(F2)` multiplies cosine and `Im(F2)` multiplies
positive sine. Native complete force-work and optional physical-viscosity heat
include all transient/DC/AC cross terms; they are not isolated cycle-average
estimates. Cached phase rotations are prepared before each native transaction
and published only after accepted motion. Refusal preserves phase, state, time
and accounts. Joule heat is not deposited into or added to structural damping.

This is an explicitly **one-way, fixed-geometry, linear forced-vibration**
response. Its mechanical source work is not a claimed closed debit from the
unchanged EM solve. It does not introduce motion feedback, sidebands,
magnetostriction, impressed-winding forces, magnetic interface traction,
moving-support reactions, deformation updates or acoustic radiation. The
existing two-tone exponential/work owner admits **at most 32 free mechanical
coordinates**, with declared payload/work/step budgets. No modes are truncated
and no sparse industrial-mesh scalability is claimed.

`tests/spatial_induction_vibration.rs` defines independent instantaneous
`J(t) cross B(t)` tetrahedral quadrature, damped scalar motion and work/heat
integrals, exact twice-frequency resonance and interval subdivision, original
mesh permutations, correspondence/support/resource refusals and atomic phase
rollback. These are synthetic implementation references, not experimental
validation of an electromagnetic device or its emitted sound.

## Source-bound rotor damping and finite material heat

`circuit_enthalpy::rotor` binds explicit pure-viscous interface channels to the
original rotor and caloric material network. Each `MaterialRotorDamper` retains
its original `InterfaceSystemCard`/`ResolvedLinearInterface`, ordered receiving
material cells, physical slip row, temperature side and authored two-body heat
partition. Exact source pins and fixed context survive endpoint T refresh.
Elastic interface coefficients, constant-loss anchors and unsupported frequency
bands refuse; no missing viscosity or heat partition is inferred.

Each thermal trial assembles `C_baseline + sum c(T) b b^T` and uses the native
`RotorBearingModel`/`RotorMidpoint` equations. Each channel's physical work is
independently `h*c(T)*(b.v_mid)^2`; it enters the two original caloric stores
once, with double-double aggregation for shared receivers. Arbitrary baseline C
and K_skew work are contracted directly as signed external exchanges. They are
never obtained by subtracting heat from total damping work or by projecting C.
Original symmetric-K storage is measured from actual endpoint configurations;
full M/K/G, spin and geometry stay fixed, so no hidden parameter work is added.

The existing `MaterialNetwork` and `cell_driver_candidate` retain caloric,
conduction, reservoir, event, nonlinear and local/global energy ownership.
Accepted motion, endpoint source bundles and thermal history publish together.
Immutable network binding prevents silent changes on continuation. Source,
receiver, clock, step/trial budget, numerical and cancellation failures publish
no half-state. A fresh native trial costs O(n³ + channels*n²), within the
existing 32-coordinate rotor bound and an explicit 256-channel cap. One
exact-coefficient native motion result is cached in O(channels+n) storage;
identical ordered c(T) tuples reuse it while source resolution, independent
work checks and trial accounting remain fresh.

This is the original midpoint mechanical discretization with endpoint material
and thermal feedback, generally first order as a coupled method. The caller
supplies a timestep ceiling and qualifies the physical response bandwidth;
tight energy residuals do not prove time accuracy, broadband material validity,
passivity of active baseline operators, thermal expansion or moving geometry.
`tests/rotor_enthalpy.rs` defines independent coupled c(T)/finite-heat solutions,
active baseline supply despite zero net damping, full mutual mass/gyro/K_skew
work, shared receivers and atomic source-domain/anchor/cancellation refusal.
These are synthetic constitutive fixtures; native execution is reported
separately from parsing and static source review.

## Finite battery-to-shaft drive with material heating

`circuit_enthalpy::converter::battery::motor::BatteryMotorThermalNetwork` advances
the original battery/converter/motor nine-state owner and finite material
enthalpy in one accepted interval. Chemical/polarization/bus charge, both winding
fluxes, converter bus charges, and finite rotor/stator momenta remain intact.
Reciprocal converter-to-motor voltage/current exchange supports motoring and
regeneration; both external shaft torques pair with their absolute speeds, while
optional regulated auxiliary demand exports electrical work.

Nine original resistance bindings refresh at the actual trial T/phase: converter
five, battery series/polarization, lead and motor winding. Fixed native bearing
drag in N m s/rad has a separate explicit receiver. All ten physical heat
integrals deposit once, including explicitly shared bodies; shaft work and
auxiliary export remain separate. Independent battery, lead, converter, motor and
whole-system balances accompany the existing local/global thermal residual.
Identical resistance tuples reuse bounded native trial results, and every trial
restarts the same accepted nine coordinates.

Original OCV/C/L, inertias, reciprocal coupling and bearing coefficient remain
fixed; no temperature-dependent storage, entropic heat, iron/PWM loss or
protection law is inferred. Endpoint resistance coupling is generally first
order. Material support, finite SOC, receiver identity, clocks, thermal events,
budgets and cancellation admit all histories together or publish none.
`tests/converter_enthalpy/battery/motor.rs` covers ten-channel calorimetry with
both shaft works and auxiliary export, finite-inertia regeneration, all-nine
R(T) feedback, shared receivers and atomic refusal/retry.


### Finite-battery PWM drive with material feedback

`BatteryMotorThermalNetwork::step_switch_candidate` and the existing
`BatteryMotorHeatNetwork` heat-assembly entry point advance the original
nine-state switched drive and all ten original caloric receivers together.
Every thermal proposal refreshes the same nine resistance source laws at its
actual material state, then restarts the native gate program from the accepted
charges, fluxes and shaft momenta. The ten original heat integrals enter the
shared finite-cell residual once; actual auxiliary demand and signed shaft work
remain external.

Averaged and switched entry points share coefficient resolution,
exact-coefficient trajectory cache, constituent work/heat gates, deposition
and publication. The generic candidate defaults to the original averaged
receipt; switched reports retain the complete cycle receipt. Each numerical
subdivision charges one circuit trial even on cache reuse; fresh source
resolution and physical-account admission still run. The represented
`start + period` must equal the requested endpoint, so local clock allowance
cannot admit a different end.

The original heat assembly continues to own contact, radiation, reservoir and
finite coolant transport, retaining its complete source sample and independent
closed-radiation check. No electrical or caloric solver is duplicated. Fixed
storage, linear motor reciprocity, authored bearing drag, endpoint-resistance
temporal accuracy and ideal complementary-gate limits remain. Every
electrical/mechanical/caloric history publishes atomically. Focused tests cover
ten-body calorimetry, actual R(T)-dependent motor motion with native replay,
finite-coolant exchange, clock/refusal/cancellation and full-subdivision trial
charging.

## Complete cooling for the finite battery-to-shaft drive

`circuit_enthalpy::converter::battery::motor::heat::BatteryMotorHeatNetwork` consumes the existing prepared `HeatAssembly` on the original nine-state battery/converter/motor and finite enthalpy timeline. Thermal contacts, integrated layers, closed radiation, prescribed reservoir exchange and finite coolant links resolve at each actual candidate state before the same nine resistance bindings and native drive solve. Cells, topology and full accepted source samples remain owned by that assembly; no additional electrical or thermal integrator is introduced. The ordinary material-network entry point shares those internals and retains its prior API through default sample parameters.

All ten drive-loss integrals deposit once. Finite coolant exchange remains internal to the combined caloric inventory; prescribed cooling remains signed external heat. Existing closed-radiation balance admission excludes its numerical net from external supply. Unowned graph resistors, electrical contacts, unsupported sources, numerical-budget refusals and cancellation publish no drive or caloric history. Fixed-storage, fixed-bearing and endpoint-resistance accuracy limitations remain unchanged.

`tests/converter_enthalpy/battery/motor/heat.rs` checks an independent eleven-cell caloric solution across continued drive steps, native-correlation coolant feedback into all nine resistances and solved motor current, separate reservoir heat, exact bare-network parity, and late typed-refusal/cancellation atomicity. The coolant provider is a manufactured fixture using the existing correlation owner and the same finite-link seam as downstream source-resolved airflow assemblies; it does not establish empirical motor or coolant validity.

## Tongue contact on the typed reed of the DG-valve file

The DG-valve file's typed reed (`reed HEIGHT WIDTH CLOSING_PA MASS STIFFNESS ZETA SOURCE_TOKEN`, one structural mode of the
`ModalAperture`) takes an optional `tongue_contact AREA EFFECTIVE_AREA FORCE RELEASE_START RELEASE` row directly after its
`contraction` row (before `constriction_noise`), the row the reed-performance file already reads in the same place. The
five numbers and their refusals are the reed-performance file's: a nonpositive area, effective area, force or release time;
an area above the reed blade area, which is the cantilever-tip factor's own domain; a negative release start; and, on the
DG route only, a release that starts after the last observed control tick. `render::schedule::reed::TongueContactInput`
admits the numbers once for both routes. The row is a lumped contact force with borrowed numbers on the oboe pages (no oboe
tongue area, force or contact time is printed); it is not a model of the tongue's own mechanics, and the lay's onset opening,
which the reed-performance `reed` row carries and the DG `contact` row does not, is not part of this change.

`TongueContactInput::force_series` compiles the closing-force history for both routes, so they cannot drift: the tongue
holds the reed against the blowing pressure's push (its net closing load is constant) until RELEASE_START, then lets go
linearly over RELEASE. The reed-performance file maps the series to `SetReedClosingForce` controls on its lumped reed, as
before (a test builds the controls of the F4 tongue card with the loop as it stood before the series was shared and compares
them sample for sample and bit for bit); the DG file gives the reed's one mode a force projection of -1.0 on blade 0, so a
closing force of F newtons is the positive force F on that port, and feeds the series to it as force events through the shared
source force scheduler (`LoadedDgPerformance::with_reed_closing_force`; a test compares those events with the series on the
same numbers). Only a card that states the row gives the reed's mode a force projection: without the row the projection stays
zero, no force schedule exists and the render is the unchanged one (the tests check the projection, the absence of events and
the absence of a force port, and a render of a card without the row was compared byte for byte, against the binary built
before this change, on one host before landing). The attachment is refused on a plate or solid source and after the first
accepted sample.

## Held start of the DG typed reed under a tongue

A DG-valve card that states `tongue_contact` starts its typed reed HELD: the reed's cold displacement is the static
equilibrium under the tongue's first load, minus the constant net closing load that the tongue holds (its force at the
first control tick plus the blowing pressure's push on the reed's face, `TongueContactInput::start_load_n`) over the
reed's stiffness, with velocity zero. The file passes it to the aperture constructor as `ReedHost::start_displacement_m`;
no parameter is added, and the card without the row starts at rest, bit for bit as before. Without this a tongue's first
load, applied at the first sample to a reed at its rest opening, presses the reed toward its lay and the reed's swept
volume radiates before the blow (read on the oboe D4 card: 0.114 Pa rms in the first 10 ms, nothing by 70 ms). A held
load that reaches the lay (the start displacement at or beyond minus the rest opening) is refused for its own reason,
not clamped: the contact would need its own law there.

## Broadband molecular absorption on finite gas storage

`molecular_heating::spectrum::MolecularGasSpectrumHeating` binds distinct original
prepared BEM absorption banks to one supplied thermochemical parcel and molecular
relaxation law. Each complete evaluated medium/reservoir snapshot must agree,
and each physical volume rule must represent the parcel volume. The caller owns
the source-model identity and sampled region; snapshot and volume agreement do
not establish provenance or spatial quadrature convergence.

Coherent phasors, independent sources and low-rank coherence factors retain their
original meaning within each frequency. Different frequencies contribute mean
powers in frequency order using double-double accumulation. Duplicate frequencies
refuse so same-frequency interference cannot silently disappear. Excitations
already include their spectral quadrature weights and use native peak units.
This is stationary ensemble or sufficiently long-time heating, not a finite-window
beat integral, transient acoustic storage or a spectral-convergence guarantee.

The sum converts to joules once, then uses the original nonlinear fixed-volume
`u(T)` inversion and relaxation-law refresh shared with the single-frequency
path. Alternatively it supplies exactly that energy and the ORIGINAL parcel to
`FiniteGasThermalEnvironment` for reciprocal finite-body exchange. Per-bin powers
are diagnostics, not extra deposits; classical shear/conduction and outgoing
radiation are not inferred as molecular heat. A subsequent feedback step must
re-solve every acoustic frequency at the accepted new state.

The bin cap bounds extra storage/work before allocation, each native readout
retains its own source/mode/product budgets, and repeated queries allocate no
field arrays or dense covariance. Borrowed loads freeze their readouts; refusal
or cancellation publishes no gas candidate and permits a complete retry.
`tests/molecular_heating.rs` includes independent nonlinear calorimetry, coherent
versus independent frequency powers, a finite-body quadratic energy reference,
source/volume and duplicate refusals, and cancellation/retry regressions.

## Complete thermal assemblies for induction motors

`circuit_enthalpy::motor::induction::heat::InductionHeatNetwork<A: HeatAssembly>`
uses the original assembly's material/caloric cells and thermal transport
providers in the same endpoint residual as the source-bound stator and rotor
resistances. Contact resistance, layers, radiation and finite-coolant
assemblies therefore reach the native flux/current/shaft evolution through
the same temperature proposal, without a separate heat timeline.

The bare-network and complete-assembly routes share the original seven-state
induction driver and one prepared machine workspace per thermal interval.
Stator copper, rotor copper and bearing heat retain their explicit receivers
and enter once. Full assembly source payloads and typed errors are preserved.
Graph resistor/electrical-contact work without a matching electrical owner
refuses. The existing closed-radiation balance treats its net roundoff as
internal transfer, retaining the original external-heat convention.
Cancellation, source refusal and exhausted work budgets publish neither
thermal nor machine history.

The original `InductionThermalNetwork` remains available through defaulted
sample/candidate source types. Both routes retain an optional original
magnetizing saturation curve through `bind_with_magnetizing_curve`, including
complete-assembly/bare-network equivalence over repeated thermal intervals.
This composition does not add cage-bar skin effect, changed geometry or
temperature-dependent magnetic storage. It retains the reduced machine's
generally first-order electrothermal feedback.
Regressions in `tests/induction_motor_enthalpy/heat.rs` compare an
independent two-coil/caloric/contact-resistance solution over two intervals,
check equivalence with the bare network when no extra paths are attached,
and exercise typed failure, cancellation, budgets and retry. Synthetic
material curves qualify the discrete coupling rather than an empirical motor.


## Material-bound interacting cavitation sources

With the existing cavitation feature, `cavitation::cluster::MaterialBubbleCluster`
retains the original liquid source and native `fs_flux::cavitation::BubbleCluster`.
Binding requires the complete same material snapshot, including viscosity and
surface tension. Each compact source uses that member's jointly solved volume
rate and acceleration, its restrained center, and the common accepted physical
clock. It never substitutes isolated-bubble accelerations or a summed source at
an arbitrary common distance.

The native cluster proposal solves the original joint DP54/work/energy kernel
without publishing its endpoint. Dropping a proposal changes no member state.
The material adapter admits every acoustic endpoint and checks cancellation
before committing all radii, velocities and the shared clock. No second solver,
cluster clone, extra integration or acoustic radiation damping is introduced.

Prepared monopoles retain separate distances and travel times in the original
liquid. Their endpoints generally have different arrival times. Coherent
pressure on one receiver clock uses the exported per-member histories and
stationary motions with the existing `MovingCompactRenderer`, which evaluates
each history at its own retarded time before summing signed pressures.
Source history resolution, nonlinear source bandwidth and interpolation
convergence remain explicit producer responsibilities. The cavitation example
includes this complete two-member source-to-receiver path.

Applicability remains restrained, separated spherical bubbles, frozen material,
incompressible long-wave mutual interaction and one-way linear compact
propagation. It does not supply finite-delay cloud dynamics, translational
forces, thermal phase change, shocks, nonspherical collapse or erosion.
The acoustic source relation is Brennen, Cavitation and Bubble Dynamics,
[chapter 3, section 3.8, equation 3.10](https://media.library.caltech.edu/CaltechBOOK:1995.001/chap3.htm).

Six regression tests cover the joint source and material identity, native
proposal rollback, independent per-path pressure and delay, late acoustic and
cancellation rollback with exact retry, mismatched observers, and receiver
history refinement against separately integrated retarded states. All six
cases passed with the original cluster solver, material coefficient definitions,
LU, math, monopole propagation and uniform-history renderer compiled through
focused dependency facades. The same run passed all 54 native cavitation tests.
This is module-scoped validation; the full Cargo workspace and IF97 example
executable were not exercised.

## Native adhesive rolling and torsion with finite-cell heat

`friction_heat::adhesive::moment::AdhesiveMomentHeat` advances the original
source-bound rigid rolling or bonded-torsion history and its actual finite
caloric receivers together. It accepts sealed native memory, including an
already loaded elastic spring. Internally it calls the original wrench-driven
mechanical preparation once; public diagnostic reports cannot supply heat or
successor histories. The source, native A/B pairing, frame and held geometry
retain their original mechanical checks. The real attached JKR area remains
conducting at zero or tensile normal resultant.

An explicit `AdhesiveMomentHeatPartition` supplies the physically justified
caloric fraction of irreversible moment loss. There is no default, including
no assumption that all adhesive hysteresis becomes heat. A nonzero remainder
requires a named noncaloric destination and persists in the accepted cumulative
loss account. The ledger begins when this coupling is bound; it does not infer
past loss from angle or relabel pre-existing spring storage. Applied COM-force
and body-torque work, recoverable bond storage, generated heat and retained
noncaloric loss remain independently visible. Reversible unloading releases
storage into native motion without generating heat.

The existing frozen-loss finite network solves one massless interface
temperature and the original cells' enthalpy increments. This uses original
areas, centre-to-face distances, optional contact resistances and endpoint
conductivities. The original direct link is suppressed only in the residual,
so conduction cannot be counted twice. Signed A/B heat follows the finite
transport solution and may be negative on a hot body. It is not the declared
caloric fraction divided arbitrarily between bodies. Other sources, latent
caloric capacity, reservoir exchange and optional caloric entropy audit keep
their original owners; the supplied callback must not reinsert mechanical
loss. The reported entropy production only includes the caloric fraction and
transport; no entropy model is invented for the noncaloric remainder.

The paired caloric curves must identify the original ordered source materials
and be the exact laws registered on the original receiving cells. Native
reaction markers, thermal face centre and ordered A-to-B normal must coincide.
Cell identity, mass, geometry, material, body assignment and fixed contact area
remain pinned. Changing topology, source, native geometry, or the caloric laws
cannot silently rebind the accepted history. Finite material-domain exhaustion,
failed source/transport solves and cancellation publish no owner. Private
candidates and the atomic `apply` path publish both native bodies, history,
cumulative loss and finite temperatures only after independent stage gates.
Unrelated native bodies remain unchanged.

The combined audit is Delta K + Delta U_bond + Delta H + noncaloric loss minus
actual applied body work minus other thermal supply. Separate mechanical,
thermal and partition residual magnitudes are also gated, so errors cannot
cancel between stages. The native exact monotone constitutive work and its
midpoint inertia discretization are unchanged. Thermal transport is backward
Euler at frozen mechanical coefficients, fixed area and held poses. This is
an explicitly split approximation: no temperature dependence, no-slip reaction,
normal motion, moving footprint, history rotation or normal/twist mixed-mode
evolution is inferred. Those require their original physical owners and
timestep/spatial refinement when coupled.

`tests/adhesive_moment_heat.rs` covers both source laws under tensile attachment:
driven yielding with zero/partial/full caloric conversion, independent finite
two-cell heat solutions with signed conduction and other supply, retained
noncaloric loss over subsequent steps, elastic storage and unloading, original
preload adoption, source/area/caloric mismatch, finite-domain refusal, final
thermal-trial cancellation and an identical retry. Its thresholds and caloric
data are manufactured discrete-physics tests, not material calibration claims.

## Shared finite-atmosphere sorption calorics

The moisture caloric owner now composes the original shared-atmosphere slab,
Cartesian backward-Euler and Cartesian SDIRK2 transports with one complete
gas-plus-specimens energy account. Each stable specimen ID binds exactly one
original isosteric-heat source on its own material card and fixed transport
context. The native common humidity solve still owns all vapor, carrier and
bound-water inventories; caloric evaluation never repeats or replaces it.

Every physical cell contributes its actual original-to-final integral of
the differential binding energy on the same original vapor datum. Opposing
specimen transfers and sealed internal redistribution remain visible even
when net atmospheric water exchange vanishes. One actual gas energy increment
enters the closure once, and one aggregate source-segment allowance covers all
specimens. Results use the native stable-ID order independently of caller or
source-binding enumeration. Source, identity, cell-budget, energy or
cancellation refusal publishes no partial transport or caloric endpoint.

Cartesian SDIRK calorics uses the complete final inventory; the first stage is
a diagnostic and is not charged again. A single-field Cartesian SDIRK entry
shares the same existing scalar-field caloric closure. Per-specimen reports
reuse the original per-material storage report, including cell contributions
and source identities. This is prescribed-temperature, rigid-system thermostat
duty. It does not infer spatial heat deposition, thermal feedback, hysteresis,
or consistency between independently authored sorption and caloric laws.

Focused regressions cover unlike specimens exchanging water at zero net vapor
change, canonical source/participant ordering, loading-dependent Cartesian
SDIRK storage, one-specimen equivalence, original first-order transport parity,
and aggregate work/source refusal followed by reuse.

## Settled three-dimensional AC induction heating with thermal feedback

`circuit_enthalpy::spatial_induction::harmonic::heat_assembly_candidate` refreshes
original scalar or oriented conductivity and re-solves stationary 3D A-phi at
every caloric Newton proposal and final acceptance. It reuses the native
source-temperature material trial, harmonic field solve and finite-cell
enthalpy residual. Actual cellwise cycle-mean Joule heat enters the assigned
caloric receivers exactly once, with the complete `HeatAssembly` evaluated at
that same proposed temperature. Contacts, layers, radiation and finite coolant
retain their source payloads and typed errors; closed-radiation net remains an
internal transfer. Unowned graph electrical work refuses.

Mapped permeability sources must be temperature independent, including when
sampled endpoint values happen to coincide. Original material domains, source
identity and whole-volume caloric correspondence are checked. The accepted
report retains original material snapshots, thermal-provider samples, cell
heat destinations, native field diagnostics and both steady AC and thermal
balances. Existing circuit-trial budgets count actual field evaluations;
cancellation, source refusal or failed nonlinear acceptance publishes no
thermal endpoint and leaves the template's transient cochain/clock unchanged.

This is a held-excitation, settled-carrier/slow-thermal approximation with
backward-Euler thermal rates. The caller must justify the time-scale separation;
no automatic settling or error certificate is supplied. Startup, envelope
magnetic-energy changes, finite-window phase/beat heat, motion, saturation and
hysteresis require their separate owners. Mean power is integrated over the
thermal interval without an additional peak-phasor factor of one half.

Three G1/G3/G4 regression definitions in
`tests/spatial_induction_enthalpy/harmonic.rs` compare independent AC diffusion
and caloric balances over consecutive intervals, distinguish conductivity
feedback from frozen-field heat, exercise temperature-dependent contact and
closed radiation transport, and cover source/fixed-permeability refusal,
transactional failure, budgets, cancellation and retry. Synthetic references
qualify the discrete model; they do not establish experimental device accuracy.

## Finite-temperature feedback on the loaded normal/twist JKR source

`friction_heat::adhesive::torsional_temperature::propose_jkr_torsional_temperature_heat`
connects the actual finite caloric cells to the original coupled I/III JKR
material queries. The input native endpoint is staged first; its original
coordinates, rates, bodies, full mobility, chronology and pending compatible
geometry remain with the native owner. The thermal phase advances from its
accepted time to exactly that native endpoint time. The incoming route must
already describe this staged mechanical endpoint's original area and coincident
ordered reaction plane. Stale normal approach, moving face centres or changed
mechanical coverage need their existing geometry/coverage owners first.

The bounded thermal Newton residual uses both actual receiving-cell Kelvin
values and a separately solved massless interface-temperature node. Those three
coordinates refresh their own original selected material/interface claims;
neither body is replaced by the interface temperature or an arithmetic average.
Each trial internally constructs the source-authenticated native continuation.
Original material-card identities, registered caloric objects, body assignments,
face geometry, source support and the complete loaded attached branch are
admitted again. No caller-authored energy report or reconstructed initial twist
can authorize a temperature successor.

Incoming bulk source T coordinates must agree with the accepted receiving-cell
temperatures, allowing only 64*epsilon relative representation roundoff. This
also applies to constant property claims: their original query must still own
T. Temperature-free source queries refuse. Thus an unrelated stale mechanical
temperature cannot silently be repaired by a nominally zero-supply heat step.
The massless interface node remains free to readjust under its physical flux
balance; it is not an extra finite thermal store.

The original source's relaxed contact area enters the existing variable-area
thermal owner during the same residual. Actual k*A/d, optional retained contact
resistance, ambient coverage and signed interbody heat therefore respond to
the source-supported temperature change. The accepted thermal state and route
retain the endpoint area and conductivity. This is backward-Euler endpoint
coverage, not a swept footprint, integrated mean area or guarantee of stability
on an unspecified intervening thermal path. Endpoint area change and signed
coefficient work each have explicit finite bounds.

JKR provides normal, torsional and surface mechanical storage and their signed
parameter work. It does not provide a contact entropy/internal-energy law.
Accordingly this adapter requires a named parameter-work source and keeps that
work in an explicit external port. It never subtracts coefficient work from
body enthalpy, invents reversible caloric exchange, or treats it as dissipation.
This reversible contact generates exactly zero heat; the callback supplies
original endpoint conductivity and OTHER thermal sources only. Optional entropy
audit covers the existing caloric/transport problem, with no mechanical entropy
or whole-system thermoelastic closure inferred.

The combined phase audit is Delta U_contact + Delta H minus parameter work
minus other thermal supply; native kinetic energy is unchanged. Independent
source and thermal residual magnitudes are gated so opposite-sign defects
cannot conceal each other. The private candidate owns the new source, sealed
native state, unchanged A/B bodies, finite thermal state, route and report.
Adopt them together and capture the next native pair from that source and the
retained state; do not initialize a new twist or clock. Original-source refusal,
thermal exhaustion, trial-budget failure or cancellation publishes no owner.

`tests/adhesive_torsional_temperature_heat.rs` independently solves the coupled
finite-cell heat equations and dimensional held-approach I/III radius equation,
checks full loaded coefficient work without caloric charging, continues the
original native preload through heating and cooling, and checks source/caloric,
body/area/clock/work-port refusals plus final-trial cancellation and retry.
The material and heat coefficients are manufactured qualification inputs.

## Original BEM harmonic microphone channels in the PCM stream

`pcm_wav::observation::harmonic::BemHarmonicPressureFrameRenderer` connects
`fs_bem::observation_transfers::PreparedHarmonicPressureSignal` to the existing
`PressureFrameRenderer`
and `render_pressure_frames_pcm16` owner. Original pressure/velocity microphone
projection, signed delayed taps, fixed source draws and relative channel phases
are retained. Native BEM angular frequencies never round-trip through Hz, and
the sink must use the declared sample clock and explicit pressure full-scale.
There is no additional WAV encoder, propagation, normalization or source onset.

Construction admits every original frequency and the complete callback work
before allocating one bounded scratch buffer. Whole-stream clock admission
precedes source/sink advancement. Successful callbacks allocate nothing; any
physical refusal or cancellation leaves the caller's entire channel block and
frame count unchanged, permitting retry. BEM refusals retain their original
`TransferError` in `RenderError::BemObservation`. Channel identities and patterns
remain available beside WAV bytes. Original producer accuracy remains unresolved.

## Forced nonlinear waveform parameter directions

`phs_orbit::PreparedPhsForcedSensitivity` prepares the complete original pHS
continuous tangent through the shared fs-orbit forced-sensitivity owner. Native
storage Hessians are discovered as in HB; an explicit complete authored Hessian
may also be supplied. One factor is reused across all physical port Fourier
partials and original storage-effort partials. The existing physical field owner
applies `G du/dp` or `(J-R) partial(gradient H)/partial p`, preserving nonlinear
modal coupling, damping and every explicit material-memory coordinate. Inputs
retain real DC and complex half-amplitudes, and unresolved nonzero forcing
bandwidth or inconsistent clocks refuse before direction evaluation.

These derivatives hold the clock and structure matrices fixed. A complete
source/port/clock residual declaration can use the shared prepared entry. All
coefficients are derivatives per physical parameter unit at fixed phase; they
are not finite-probe responses, peak-amplitude derivatives or stability evidence.
The same source must remain unchanged, and harmonic order needs independent
refinement. Tests cover phased multiharmonic force/stiffness mechanics and a
nonlinear original storage with an explicit relaxing strain, including its
material derivative, native/declared Hessian agreement and bandwidth refusals.

## Resolved thermoviscous interior with full aperture reaction

`cavity_modal::interior::thermoviscous::exterior::RadiatingThermoviscousInterior` connects the original source-bound P2/P1/P2 velocity/pressure/temperature owner to the native BEM self/mutual aperture matrix. One original loaded fluid LU and one small port factor serve coherent feed flow, original structural force, body force and heat-release excitation. Molecular Cp, shear/bulk viscosity, fluid and finite-solid conduction, structural damping and Maxwell memory remain on their original operators. No pressure-only replacement or added attenuation changes the fluid equations.

Radiators must pass the native unit, planar, exclusive traction-port admission. `plane_wave_port_triangles` exposes only those original physical triangles and outward normals; no second fluid mesh is assembled. The BEM applies the explicitly declared constant-velocity lift U/A on the same triangles, with U=-Q. This scalar work pairing does not establish a uniform resolved opening velocity: transverse opening motion outside that lift remains a port-model approximation. Signed/scaled/overlapping traction modes require a dual map and refuse. The caller supplies the independent homogeneous lossless exterior `Medium`; interior molecular loss does not imply an exterior absorption law. The exterior stays rigid outside selected openings, including moving structural walls.

The native outgoing BEM matrix is conjugated once for exp(+i omega t), with both geometric area factors retained. Original physical interface residuals are checked in power coordinates. Driven-port work, original body/heat/structural source work, exterior boundary work and native loss channels are accounted once. Conductive/caloric availability and numerical residuals are not deposited as additional heat. The three focused tests compare complete fields against independent native traction-admittance/BEM equations, vary numerical references, retain coherent Maxwell/body/heat source work, and exercise physical geometry, owner, budget, band and cancellation refusals. These are numerical regression definitions, not mesh-convergence or measured-validation claims.

## Full material damping in the resolved thermoviscous interior

`ThermoviscousInterior::with_structural_viscosity` binds the existing
`ModalPlateViscosity` passive channels to the original resolved moving-wall
coordinates. Each channel is scaled by the square root of its original
reference resistance, then admitted once by
`ThermoviscousFe::with_moving_wall_damping_channels`. The full positive
semidefinite material matrix, including signed cross-mode terms and loss null
spaces, reaches both harmonic impedance and the original causal energy/heat
operator. No additional constitutive fit, modal diagonalization, or per-step
material evaluation is introduced. Resolved fluid losses remain separate.

The existing modal damping rates must be zero, and a second material attachment
refuses. Channel mode count must match; physical basis identity remains the
caller's obligation. The same existing `ModalPlateViscosity` and FE channel
limits apply. The focused regression varies reference resistance over 240
orders of magnitude while holding the physical channels fixed, and checks
harmonic response, causal motion, heat, and extent refusal. It is authored but
not claimed executed: native Rust validation remains blocked by the bounded
shared dependency build.

## Axial Maxwell history on native finite-rotation rods

`rod_motion::PreparedRodAxialRelaxation` binds one original
`ConcentricPronyAxialSection` to each native reference segment. Equilibrium
Young stiffness stays with the existing elastic rod; each retained material arm
adds its own area-weighted excess axial rigidity, relaxation time and independent
history. Segment lengths, complete material queries, geometry, original elastic
selectors and source receipts must agree. A finite-frequency sample cannot
silently authorize this transient law: explicit broadband source support is
required. The scalar small-strain axial material model does not infer hereditary
shear, bending or torsion.

`RodMotion::with_relaxed_axial_relaxation` declares relaxed initial memory;
`with_axial_relaxation` accepts explicit dimensionless spring strains in original
segment/layer/declaration order. Both attach only at initialization. The native
motion state retains history, storage and cumulative physical heat through
clones, accepted continuation and instantaneous momentum changes. Rejected
candidates do not mutate that accepted state. Fluid refresh rebuilds the native
response with the same memory; changed solid source queries refuse. The existing
finite-enthalpy composer refuses this history until an explicit hereditary
caloric coupling can retain its physical heat.

The actual geodesic-midpoint axial rate and its full work-conjugate force/couple
pullback come from `fs-solid`'s original finite-rotation kinematics. The original
centroid/spin mobility and support-impulse owner are shared with axial viscosity.
Each frozen-pose arm forms a reciprocal two-state kinetic/spring system and uses
the existing exact quadratic pHS propagator. A zero-held reciprocal observation
port supplies the integrated force impulse directly; subtraction of nearly equal
endpoint velocities does not determine it. Under strong momentum cancellation,
retained native covector invariants and the exact responsive endpoint preserve
representable motion tails; weak updates and physical support accounts keep the
original integrated impulse. The same original exact owner
independently integrates nonnegative material heat. Symmetric original-arm sweeps
surround the native equilibrium-elastic mechanical step. Branch forces and strain
advection are not repeated in that mechanical step.

Returned `RodMotionStep::axial_relaxation_heat` separates original branch and
segment heat; the original energy gate includes branch storage and this physical
heat while retaining the signed numerical integration defect separately. Existing
Kelvin loss channels remain independently selectable. No residual is renamed as
heat and no heat is implicitly deposited in a thermal owner. Fixed supports
retain the complete time-integrated material reaction. Clamped arms use the
original exact held-strain Prony map, including positive heat and force integral;
zero-strength arms retain their authored history without force or dissipation.

The symmetric composition is second order for resolved smooth motion and still
requires elastic-wave and geometric-splitting refinement. The local reciprocal
flow treats fast Maxwell dynamics jointly: at fixed excess stiffness its force
and heat vanish as relaxation time tends to zero, while fixed excess stiffness
times relaxation time approaches the physical Kelvin resistance. No explicit
branch-force kick leaves artificial step-dependent drag. This does not establish
exact complete moving-rod trajectories, finite-strain material calibration or a
general hereditary tensor law. Existing pHS exponential resource and energy gates
remain authoritative and can refuse unrepresentable intervals. Immutable exact
plans reuse matching mobility and duration across candidates; one caller-owned
workspace serves all arms without warm per-arm allocations, and original
strain-covector scratch remains bounded on the stack.

## Physical-electrode settled AC heating

`circuit_enthalpy::spatial_induction::harmonic::terminals` adds
`voltage_heat_assembly_candidate` and `current_heat_assembly_candidate` to the
same material admission, complete `HeatAssembly`, cellwise heat routing and
enthalpy Newton residual. Every temperature proposal and final acceptance
refreshes the original electrode owner and solves impressed and electrode
excitation coherently in one native A-phi field. Total-current drive retains
the solved spatial current distribution and scalar electrode unknowns; it does
not freeze terminal impedance or impose uniform current density.

Actual original-order terminal reactions and separate signed impressed/terminal
real and reactive supplies remain in the accepted report. Only native conductor
Joule heat is deposited, exactly once. Fixed permeability, stationary geometry,
whole-volume solid caloric correspondence, original source identity and the
settled-carrier/slow-thermal approximation retain their existing admission.
Neither the original magnetic history nor its clock advances. Electrode contact
impedance, lead fields, loaded circuit heat, startup and envelope magnetic-energy
exchange require separate physical owners.

Two focused regression definitions in
`tests/spatial_induction_enthalpy/harmonic/terminals.rs` compare independent AC
voltage/current and two-body caloric references over consecutive accepted
intervals, conductivity feedback versus frozen heat, coherent cross heat,
separate source work, common-mode invariance, typed refusal and cancellation.

## Heterogeneous directional liquid/rigid drag

`liquid_momentum::heterogeneous::apply_rigid_drag_native` and
`liquid_momentum::material::apply_heterogeneous_rigid_drag` extend the same
compact native momentum transaction to one authored principal resistance
tensor per selected wet cell. Each tensor is a proper world-frame rotation
of three finite nonnegative inverse-time rates. Exact zero rates retain
free-slip directions. Partial cells use their original tracked liquid mass,
and the rigid reaction retains the complete original mass and inertia
tensor. The material entry derives every physical scale and the speed
domain from the original liquid source.

Each local cell/body exchange uses the symmetric kinetic metric of the
resistance factor and relative mobility; resistance and off-center mobility
need not commute. Principal-rate normalization and kinetic-mode normalization
avoid inverting zero damping or forming an overflowing inverse weak rate.
Exactly integrated local impulses supply equal and opposite body force/couple,
and a separate analytic integral supplies nonnegative drag dissipation.
A forward-half/reverse-half composition, with the central pair combined into
one full exchange, uses at most
`(2 * selected_cells - 1) * substeps` exchanges. The caller supplies the
positive substep count and an explicit work ceiling. Response factors are
prepared once per selected cell; storage is linear in stencil size.

The heterogeneous simultaneous source is second order under subdivision,
not globally exact. Stiff local modes decay without timestep-induced modal
sign reversal, but a coarse multi-cell composition still requires temporal
refinement. Pose, grid geometry, clocks, original material/force bindings and
thermal state remain fixed. This adds no impermeable boundary, excluded
volume, guessed viscosity, heat partition or full flow/drift integration order.
The existing uniform-rate six-coordinate exact source and its public report
remain available unchanged.

The same final native population readback and original rigid impulse receipt
independently gate linear momentum, angular momentum and kinetic work before
either owner publishes. Cumulative analytic local dissipation is also compared
with the actual aggregate endpoint kinetic loss; conservation or dissipation
errors cannot be hidden by assigning heat. Cancellation, invalid coefficients,
unresolved kinetic modes, unrepresentable positive exposure and exhausted
work budgets preserve both owners.

Focused regressions cover exact finite-mass recoil, noncommuting anisotropy
against an independent simultaneous exponential, second-order heterogeneous
refinement, full-inertia momentum/energy closure, stiff null directions,
weak-rate extreme scales, underflow refusal, bitwise idle behavior and
cancellation/work limits. Native liquid/body and original-material regressions
cover the actual partial-liquid inertia, physical unit bindings, reciprocal
accepted momentum and kinetic receipts, and state-preserving refusal.

## Original material-shell modes in resolved thermoviscous fluid

`cavity_modal::interior::thermoviscous::shell::ShellWallAssembly` transfers a
`MaterialShellImage` into the existing mixed pressure/velocity/temperature
owner. The source's complete signed, mass-normalized global P1 translations
become hierarchical P2 vertex coefficients; edge bubbles remain zero.
Rotations retain their original structural influence and are not mistaken for
midsurface translations. The native no-slip congruence supplies reciprocal
pressure and shear reaction, so zero compact volume participation does not
remove a mode's spatial fluid loading.

An explicit injective node map covers exactly the source's oriented cavity
facets. Coordinates must coincide, each mapped facet must be an actual exterior
no-slip face with the correct outward side, and an independently accumulated
signed volume-work witness must agree with the original source projection.
Incomplete maps, stale geometry, implicit welding, offset skins and nonmatching
meshes refuse. All original modes and frequencies survive; native silent-trace
and stationary-intersection refusals remain. Material passive channels and
original interface dashpots enter the full damping operator once. Interface
stiffness already belongs to the original eigenbasis and is not reapplied.

The assembly binds the existing gas, molecular-relaxation or generic material
fluid owner. The requested band must fit the original shell material and
interface support; sampling metadata is caller-declared without advancing a
source clock. The initial structural law envelope is instantaneous Kelvin--Voigt
material and linear interfaces. Frequency-local coating/constant-loss slices
and separate substrate/coating memories refuse because their original causal
storage/loss is not supplied by this attachment.

`ShellFluidInterior` retains the original image/map, offers allocation-free DD
projection of physical force-port amplitudes, and exposes the native harmonic,
causal, field and duct operations. Its consuming `into_interior` permits further
existing native source attachments. Explicit source, mapping/work, byte, mode,
channel and native FE budgets bound setup; cancellation publishes no partial
source or force projection. No new eigensolver, material fit or warm assembly
is introduced.

Four regression groups specify a genuinely reduced material shell's tangential
loading with zero compact volume, full material/interface damping, force-work
reciprocity, native harmonic/time parity, vertex relabeling, source binding and
geometry/law/resource refusals. Source/API and independent physical review,
pinned formatting and whitespace checks passed. Native Rust execution remains
unverified under the known real-dependency build memory limitation.

## Finite recovering creep on the shared caloric and contact timeline

`FiniteThermoplasticSolid::try_new_with_recovery` binds exact recovering finite
model cards to the original caloric chart, reference density and mesh inventory.
The existing material/solid/candidate types retain their old defaults and now
carry a complete typed material family. Total caloric plus mechanical capacity
must stay positive over the admitted domain. Physical creep/yield losses and
signed reversible exchange are deposited once; endpoint-work defects and
numerical kinetic loss remain outside caloric heat.

The original feedback solver repeats trials from the same accepted material,
thermal and observer state. Thermal substeps do not repeat the physical recovery
interval. Accepted candidates commit deformation, permanent histories, resistance,
temperature, optional momentum and observer history together. The heated-body
regression independently solves finite equilibrium, primary/recovery flow and
caloric balance through loaded heating, unloading, hot dwell and reload for both
standalone and simultaneous yielding families.

Existing pressure, rigid-tool, impact and deformation-dependent conduction
entrypoints carry the same generic family. Conductivity and admissible current
normal-flow geometry remain with the original observer/callback. Thermal and
isothermal Tet4 impact/surface/scene adapters preserve complete inelastic memory
through zero-duration events; Hex8 support is unchanged. These extensions infer
no impact-heat partition, distorted-tet Fourier law, remeshing, damage or rupture.

## Source-owned material-resistor heat from native AC phasors

`material_circuit::harmonic::account_harmonic_resistor_heat` retains actual peak
complex branch efforts and original resistor source/geometry. Mean heat comes
from the original conductor kernels applied to both quadratures with the peak
factor one half. Complete unique branch/source correspondence, source-owned
current, real/reactive work and summed native resistor loss are independently
checked. Original storage/source branches and canceled field reference
terminations create no resistor heat. Instantaneous E/J-, frequency- or
rate-dependent sources refuse, including applicability-only restrictions;
resolved static temperature-dependent resistance remains supported.

Field-aware source lowering skips numeric inductance checks only at explicitly
named, unique native field-owned inductor slots. It never stamps substitute
storage. Ordinary material circuits retain all physical L/C checks and existing
bounded graph capacities. Focused regression definitions cover source/phase
orientation, peak scaling, original branch ordering, placeholder ownership,
source-domain refusal, inconsistent receipts and cancellation.

## Reciprocal drag on thermal and non-Newtonian liquid owners

Both `liquid_momentum::material::apply_rigid_drag` and
`apply_heterogeneous_rigid_drag` accept the sealed
`fs_lbm::material::free_surface::MaterialFreeSurfaceMomentum3` interface.
The existing exact uniform source and bounded directional source now operate
directly on frozen material, source-rheology and enthalpy free-surface owners.

The physical unit map and wet-cell stage belong to the same original liquid.
Actual partial mass reacts against the finite body's original mass and complete
inertia. Existing actual linear/angular momentum, kinetic work, independent
dissipation and paired-publication gates are reused without an owner clone or
duplicate mechanical solver. Accepted viscosity/pressure samples, thermal
storage, source identity and chronology survive; reported source dissipation
still requires an explicit thermal allocation.

Integration regressions use previously evolved rheological and thermal/rheology
owners, compare recoil and dissipation with independent finite-mass solutions,
retain material and thermal histories, and advance the original flow afterwards.
A directional endpoint-speed violation occurs after the rigid candidate is
formed, checking refusal of both owners, retained enthalpy, cancellation and
exact retry.

### Native Hex8 plasticity with reciprocal spatial heat

`solid_enthalpy::hex::ThermoplasticHexSolid` binds all eight original Hex8
quadrature histories to distinct active caloric cells. Each cell must reproduce
its own native reference integration volume and source density; matching only
whole-element mass is insufficient on a warped brick. Original caloric chart
objects, source families, geometry and formulation remain bound. Virgin binding
uses the native geometry-only quadrature accessor and initial histories, without
a material increment or elapsed physical time; initially stressed sources refuse.

The native Standard/B-bar thermal Hex8 owner supplies equilibrium, or
consistent-reference-mass backward-Euler motion after explicit initial velocity.
The original Picard/IQN-ILS feedback driver now also accepts exact enthalpy
coordinates; existing temperature-only callers retain their original behavior.
Every trial reuses the same accepted displacement, complete material histories,
velocity, enthalpy and clock. The source law receives exactly one physical
duration. Physical plus signed reversible heat enters each original cell once
through the existing conduction/events/observer timeline. Source-pinned
conductivity can lower to the existing material-network observer.

Free-energy plus temperature-times-mechanical-entropy change forms mechanical
internal energy. Material and inertial numerical losses remain separate from
physical heat. Native work, per-cell heat mismatch, thermal balance and combined
balance must all close before the complete caloric/mechanical/observer candidate
returns. Late source, observer, convergence or cancellation refusal changes no
accepted state. A proposal evaluated at another rounded enthalpy cannot masquerade
as the returned exact caloric endpoint.

The supplied point network is not a generated geometric dual mesh. Thermal
geometry, conductances and spatial refinement remain explicit; this is fixed
reference, fully solid, small strain and rotation, without moving transport
geometry, finite forming, contact or phase change. Four authored regressions
cover force-held self-heating/softening/permanent set, timed relaxation,
nonuniform warped B-bar point heat/conduction and consistent inertia with atomic
late refusal. Native Cargo execution remains a separate verification requirement.

### Hardness on the actual accepted Hex8 mechanics/heat state

`ThermoplasticHexSolid::bind_hardness_field` pins one exact original
indentation selection per element and admits all eight native quadrature
points. `HexThermoplasticHardnessBinding::at_accepted` observes each point's
actual caloric temperature and accumulated plastic strain directly from the
private accepted coupled owner. Original mesh slices, formulation, ordered
shared material families, cell identities, clocks, reference geometry/mass,
caloric law objects and exact accepted enthalpies remain binding. A different
family allocation or an edited thermal inventory cannot substitute for these
owners.

`check_candidate` consumes the existing converged
`HexThermoplasticCandidate` and returns that unchanged candidate together
with the complete hardness field only after every endpoint succeeds. It does
not reconstruct a thermal step, repeat a return map, rerun transport, clone
native history/energy, or alter the original heat, work, inertia, observer and
feedback reports. A late unsupported hardness point returns no advanced owner.
The reusable binding shares the original family allocation across adopted
successors.

Coverage is checked at actual accepted endpoints; it does not establish
unobserved path coverage, uniform material hardness or empirical qualification.
Standard/B-bar heated loading and reversal regressions specify heterogeneous
point temperature/alpha, pinned test identity, unchanged native accounts and
combined energy balance. Additional cases place a support boundary before the
last point and refuse changed families, geometry/mass and stale thermal clocks
without partial adoption.

## Loaded physical-electrode AC circuits on the finite caloric residual

`circuit_enthalpy::spatial_induction::harmonic::circuit::heat_assembly_candidate`
closes the native affine electrode field, original RLC/transformer circuit and
complete `HeatAssembly` on the same proposed temperatures. Each proposal refreshes
the original conductor conductivity and every source-bound external resistor,
then recovers actual loaded fields and branch efforts. Spatial Joule heat and
original resistor heat enter their finite caloric receivers exactly once.
Graph-source plus impressed-field work is external supply; electrode transfer
remains an internal field diagnostic and is not added to total supply again.

Graph resistor destinations come from the same material bindings that resolve
their laws. Every resistor needs one complete whole-volume receiver separate
from resolved magnetic bodies. Source identity and geometry stay pinned across
trials. The accepted state is source-lowered before admission, so consecutive
thermal intervals reuse the original immutable graph template. Instantaneous
E/J-, frequency- or rate-dependent resistor laws refuse before field work;
electrical interface contacts retain their separate owner and are not discarded.

The shared existing field-material refresh and Joule recovery also serve the
isolated, prescribed-current and prescribed-voltage thermal entrypoints. No
second electromagnetic, circuit or thermal nonlinear solver is introduced.
Each loaded thermal proposal reserves and reports one native solve per port,
one additional solve for nonzero independent excitation, and one actual loaded
recovery. Discarded Newton stencils count too; descriptor and assembly work keep
their own native limits. Cancellation and typed provider refusal publish no
thermal endpoint and leave the borrowed field and graph templates unchanged.

The admitted regime is stationary, fixed-permeability, settled-carrier thermal
feedback. Reactive field/network energy remains diagnostic; carrier startup,
envelope-storage change, magnetic entropy, hysteresis and moving occupancy need
their physical transient owners. Complete thermal contacts, radiation, layers
and coolant retain the original heat provider and its same-temperature residual.
Independent regression definitions solve a driven and a passive loaded cube
against separate field/circuit and two-body caloric equations over consecutive
intervals. They check both temperature-dependent electrical laws, source work,
internal transfer, local heat, reactive energy, original source identity,
ownership refusals, exact work counts, typed failures, cancellation and retry.

## Source-shell relaxation in resolved fluid coupling

`ShellWallAssembly` now admits the MaterialShellImage's explicitly selected
coating and substrate relaxation owners through the original native moving-wall
memory equations. Every signed mass-normalized c row and tau is retained verbatim,
in coating-channel order followed by substrate-channel order. Equal-time and
zero rows remain distinct histories. Original modal frequencies already contain
relaxed stiffness; no unrelaxed stiffness, refit or second normalization is added.

The original instantaneous material channels and interface dashpots remain once
in the full damping matrix. Selected material memories already exclude their loss
from that projection. Causal preload is explicitly x=c dot q-z in
wall_relaxation_coordinates_sqrt_j with the same channel order; it is never inferred
from displacement or silently reset. Native harmonic/causal heat accounting and
all existing pressure, velocity, field and duct operations use the same interior.

ShellWallControls carries explicit complete relaxation channel, workspace and
projection limits in addition to existing mapping/native budgets. Binding admits
the requested drive band against original region, interface, coating and memory
fit support. Passive continuation does not widen that declared support. Unlowered
frequency-local coating or constant-loss slices still refuse instead of silently
losing material storage. Existing conforming geometry/source identity obligations
and complete-mode admission remain unchanged.

## Periodic quadratic Hamiltonian modulation and parametric excitation

`phs_orbit::parametric` admits a complete original `PortHamiltonian`, its
`PeriodicPortInput`, and an explicitly authored physical energy addition
`H(t,x) = H0(x) + x^T Q(t) x / 2`. This enables periodically varying stiffness,
inverse inertia and coupled energy terms without replacing the original
nonlinear storage, changing coordinates, eliminating memory states, fitting a
frequency shift, or introducing another numerical solver. The author remains
responsible for the constitutive meaning and support of that energy addition.

`PeriodicQuadraticModulation` stores complete row-major Fourier matrices in
original state units: `Q(t) = Q0 + 2 Re sum_{n>=1} Qn exp(+i n omega t)`.
A peak cosine matrix therefore supplies half its peak amplitude at the positive
harmonic. Every complex matrix is exactly **transpose symmetric**, including its
imaginary part; Hermitian symmetry is a different proposition and is refused.
DC must be real. Units are energy divided by the product of the two corresponding
state units. All supplied entries are finite, every matrix has the original
square shape, and the pump and port drive share the exact same positive finite
clock. The period and active harmonic frequencies must also be representable.

Source admission caps the original dimension at 256, positive harmonic index at
256, and total supplied full-matrix complex entries at 1,048,576. Preparation
keeps only exactly nonzero upper-triangle entry sequences and their actual
nonzero bandwidth; it applies no magnitude threshold, frequency clipping or
symmetrization. Reusable source scratch caches one physical-time matrix and,
when needed, its analytic time derivative. Positive definiteness is deliberately
not inferred or imposed: signed modulation, including an instantaneous
stiffness crossing zero, needs its author's physical justification. Symmetry
alone establishes neither passivity, bounded motion nor Floquet stability.

The physical field is evaluated as
`(J-R) (gradient H0(x) + Q(t)x) + G u(t)` through the **original compensated pHS
field owner in one contraction**. The complete source Hessian action is retained
and `Q(t)` is added exactly once before the same `J-R` map. Native storage
tangents survive boxing, nonlinear terms and explicit memory coordinates.
An authored Hessian override must cover all of H0 in the original physical
coordinate order. A declared tangent refusal remains a named solver refusal;
there is no numerical replacement for an authored action that fails.

`power` reports the original H0 plus the quadratic energy, explicit pump power
`x^T Qdot(t) x / 2`, original physical port power, resistance power evaluated
using the complete modulated effort, and total energy rate from that same effort
and actual original field plus pump power. Pump work is external supplied work;
it is not dissipative heat or unexplained energy creation. The reported energy
rate is evaluated independently of the port-minus-resistance expression, so it
does not manufacture closure by construction. These are instantaneous
diagnostics, not an integrated-work, stability or passivity certificate.

The prepared source delegates fixed-order and adaptive physically scaled HB to
the existing AFT/Newton owners, variational shooting to the existing original-ODE
RK4/monodromy owner, and physical-port sideband response to the existing periodic
linearization/complex-factor owner. No state reduction, new Newton iteration,
Fourier transform, time integrator, linear solver or eigensolver is added.
Initial HB/base/sideband orders must resolve every actually nonzero authored
drive and pump harmonic. Shooting additionally requires the complete pump and
drive bandwidth below its declared grid's Nyquist index; that admission alone
does not establish temporal accuracy. High, tiny pump coefficients cannot be
silently dropped to make these gates pass. Adaptive HB still measures actual
successive solves, full equations, omitted harmonics and independent quadrature
changes. Callers must supply their intended period chart and appropriate seeds;
a response driven at twice that chart's frequency may be subharmonic relative
to the pump, without implying primitive-period or branch selection.

Physical sideband probes use the shared owner's **peak complex amplitudes**
at signed `carrier + k*omega` frequencies. The original `G` maps every real and
imaginary port component into the complete state forcing. Both converted
positive and negative sidebands remain present. Adapter forcing scratch and a
declared mapping score `64*(2N+1)*dimension*max(ports,1)` are charged against the
supplied sideband payload/work limits before the shared owner runs. Already
prepared source storage and source callbacks' internal work are excluded, as
in the original shared-owner contracts. Cancellation and original source,
conditioning, residual, sampling and omitted-band refusals remain visible.
Finite-order responses and finite-grid multipliers do not imply a converged
continuous spectrum, stability, passivity, attraction or a quasiperiodic model.

Focused regressions compare a forced damped Mathieu response against an
independently written time-domain physical equation at two RK4 resolutions;
compare explicit pump work against independently integrated energy change;
retain a nonlinear cubic storage and its explicit relaxing material state under
actual adaptive HB; check mixed-coordinate energy/tangent terms; compare the
shared variational shooting state with HB and its monodromy determinant with
Liouville's identity; and check physical-port conversion into both signed
sidebands against weak-pump perturbation coefficients. Admission regressions
cover Fourier symmetry/DC, clock mismatch, tiny high-harmonic bandwidth,
explicit sideband payload limits and declared-Hessian refusal.

## Coherent three-dimensional induction heat in finite caloric cells

With `harmonic-winding-heating`,
`circuit_enthalpy::spatial_induction::harmonic::coherent::heat_assembly_candidate`
connects the original `SpatialHarmonicInductionExposure` bank to the existing
finite-cell enthalpy/Newton owner and complete `HeatAssembly`. One native
observation integrates the actual thermal interval on the authored source clock.
Spatial P1 Joule heat, carrier terms and inter-frequency beats survive; the
observation and cell routes are reused throughout thermal trials. Original
cycle means, magnetic storage and mechanical work are not deposited again.

Every conductor requires an original whole-volume solid caloric receiver.
Prescribed conductor motion refuses. Original scalar/oriented electrical
sources and induced-linear magnetic laws must be independent of the evolving
temperature, preserve their selected sources and resolve to the same prepared
coefficients at every proposal. Coincident endpoint values or zero heating do
not bypass the source-law admission. Material identity, nonthermal conditions,
fixed density/volume and caloric domain remain checked by the existing owners.
Source-temperature-dependent heat transport, contacts, layers, closed radiation
and finite coolant continue to resolve at each actual candidate temperature.

The public report retains the original exposure and accepted material payload,
plus deposited heat, other external heat and a corrected caloric residual.
It closes only thermal storage against these transfers; it does not infer
electrical supply, electromagnetic storage changes, entropy, startup or feedback.
The original backward-Euler thermal discretization still needs time refinement.
Refusal/cancellation returns no endpoint and leaves accepted field and caloric
histories reusable. Tests compare two consecutive beat-heated intervals against
independent physical-time integration and nonlinear contact/radiation balances,
and cover evolving-source refusal, late heat failure, cancellation and retry.

## Closed two-fluid bulk loss on the native thermal timeline

With `superfluid`, `advance_closed_friction_heat` also routes the original
channel's configured longitudinal shear/second-viscosity heat to its declared
`DepositionSimulation` cells. It uses the existing per-cell physical Joule
intake exactly once per accepted interval, preserving the native caloric law,
timeline and source accounts. Fourier conduction destroys quadratic wave
availability while redistributing thermal energy; its separate channel loss is
included in the channel balance and is never deposited as additional heat.

The composition remains a frozen-reference, first-order endpoint heat split;
it does not turn the linear perturbation temperature into absolute caloric
storage or relinearize material coefficients. Source-domain refusal and late
cancellation preserve both accepted owners and all heat/conduction histories.
`tests/superfluid_heat.rs` checks repeated physical bulk-heat intake, distinct
conductive availability loss, native-law retention and repeated thermal refusal.

## Transient porous absorption on native acoustic and thermal cells

`porous_heating::transient::TransientPorousHeating` binds each complete original
porous DG element to one active, adjoining, full-area native thermal prism. Its
optional local positive-source caches integrate Darcy drag and all original
viscous/thermal memory losses along the unchanged complete DG trajectory. One
total dense byte/work allowance is partitioned across these caches; the native
128-state source-integral cap and original constitutive source band still apply.
The held ports used for heat are exactly those passed to the original acoustic
propagator. Source observation and thermal preparation never alter that step.

Each original element's integrated Joules become interval-average watts on the
existing thermal event/CFL timeline, preserving cell masses, storage laws,
conduction, boundaries and authored sources. The acoustic interval must be
refined to resolve thermal source timing; refining thermal substeps alone does
not remove the interval-average approximation. This is frozen-material one-way
coupling, without inferred fluid/skeleton temperature partition or feedback.
Each thermal step receives its represented share of the original source Joules
through the native interval-source hook, using the exact duration its integrator
will apply. Whole and substep clock differences must resolve their original
durations within `sqrt(f64::EPSILON)` relative error; larger epoch-rounding
distortions refuse and require a local time origin.

Both acoustic and thermal candidates publish together. Physical mechanism/cell
heat and actual thermal credits are reported separately from the original port
work and native thermal audit, which already contains deposited heat once.
The signed remainder of acoustic energy plus work minus stored energy and
porous heat contains any other loss and numerical residual; it is explicitly
not an independent loss measurement and never becomes deposited heat. Clocks,
geometry, mass and thermal-law binding are checked, including final events.
Focused regressions exercise localized startup, native mechanics equality,
held ports, conduction and source events, excluded matched-boundary loss,
latent storage, shifted-clock energy preservation and underresolved-clock
refusal, cancellation and rollback on thermal/material refusal.


## Cavitation populations at changing operating pressure

`PolydisperseBubblyLiquid::at_equilibrium_pressure` continues the original native
number distribution through stable equilibria and binds the result to the SAME
material source. Every retained gas inventory, radius domain, original class
order and absolute number density stays unchanged. Updated radius, occupied
volume, restoring stiffness, viscous damping and Keller--Miksis compressible
coefficients enter coherent propagation together. The native root work occurs
once during preparation; frequency response retains its allocation-free,
fixed-order compensated accumulation.

`DiluteBubblyLiquid::from_population_class` consumes the current equilibrium
of a native class, preserving one-bin equivalence with the polydisperse owner.
It does not recreate a bubble at the gas reference radius. Source identity and
all existing dilution, amplitude, compactness, homogenization, neighboring
scattering and multiple-scattering limits remain in force. The complete new
population is re-admitted before publication, so decompression cannot silently
turn a previously dilute cloud into an admitted dense cloud. Cancellation or
failure leaves the previous immutable acoustic owner usable.

The same pressure-dependent coherent wavenumber supplies native segment-wave
momentum impedance and the existing native duct-chain entry point. Scattering
loss is still not an unrestricted thermal deposition source. These are
held-volume, frozen-liquid equilibrium queries with an explicitly authored
polytropic gas law, not closed-parcel compression, evolving temperature,
phase-change kinetics, gas diffusion or an interacting cloud model. Both
requested and returned equilibrium pressures must satisfy the native original-
reference pressure/bulk-modulus bound, independent of numerical root tolerance.

The `cavitation` example adds a fixed-number two-class pressure sweep reporting
radii, void fraction, phase speed, coherent attenuation and native root residuals.
Run `cargo run -p fs-couple --features cavitation --example cavitation` in a
complete constellation checkout. Focused regressions independently reconstruct
pressure-dependent complex susceptibility, Foldy dispersion, scattered power
and liquid-momentum impedance; they also check one-bin equivalence, original
gas/source retention and transactional dilution/cancellation refusals.

All 65 native cavitation tests and 12 focused bubbly-response tests passed using
the original numerical sources with restricted dependency facades. Three
existing full duct-chain/IF97 tests are excluded from that focused executable.
Full Cargo/DSR and the complete IF97 example remain unverified because the seven
constellation siblings and DSR/RCH are unavailable. No experimental-validation
or model-error certificate is claimed by these numerical checks.


## Material-driven thermal response of bubbly liquids

`cavitation::bubbly::thermal` prepares both `ThermalDiluteBubblyLiquid` and
`ThermalPolydisperseBubblyLiquid` by borrowing the accepted original acoustic
owner. Each occupied class receives its own coherent gas material snapshot at
its current gas partial pressure and liquid/interface temperature. Empty classes
require no fictitious gas source. Preparation retains the native isothermal
inventory, current radius, absolute number density and original liquid source.
A changed operating pressure requires re-evaluated gas snapshots from the
original material law; stale reference gas density/transport is rejected.

The native spherical heat equation supplies complex gas stiffness and thermal
absorption. Both thermal and original polytropic paths share the same native
radius, amplitude, wall-Mach, pressure-jump, compactness and scattering admission,
compensated Foldy accumulation, and existing duct-wave/duct-chain propagation.
The full complex wall pressure enters the native pressure-jump bound. Splitting
a histogram class cannot change its physical response or relax collective limits.
No per-frequency allocation, material query, root solve or state mutation occurs.
The population preparation has one bounded vector and cancellation per class.

`ThermalBubblyResponse` retains the admitted acoustic response and separate
thermal/viscous absorption coefficients. `mean_absorption_w_m3` evaluates local
cycle-mean power only for an actual peak pressure within the admitted envelope.
Gas thermal loss follows the native pressure-volume work, while Newtonian liquid
loss is the cycle average of 16 pi mu R Rdot^2. Outgoing radiation and coherent
scattering are not deposited again as heat. These local coefficients do not
invent a solved duct field, a spatial heat integral or a finite body's evolution.

The model remains a dilute, small-signal, ideal noncondensable gas with an
isothermal interface and negligible gas inertia and vapor thermal contribution.
Gas acoustic compactness is checked separately from liquid compactness. It does
not extend to nonlinear thermal collapse, vapor kinetics, liquid thermal
impedance or interacting clouds. The zero-population path is the exact original
liquid, and cancellation/failure publishes no partial preparation or response.

Four new coupled regression groups independently reconstruct spherical diffusion,
complex dispersion, pressure-volume thermal work, viscosity cycle quadrature and
native liquid-momentum impedance. They also cover mono/poly and split-bin parity,
current-pressure gas refresh, source/inventory retention and transactional limits.
All 16 focused coupled tests passed alongside all 71 native cavitation tests,
using original numerical and verbatim material sources with restricted dependency
facades. Three full duct-chain/IF97 tests remain excluded. The example's updated
pressure-sweep function, including a gas-diffusion spectrum and separate heating,
compiled and ran against explicit liquid coefficients. Full IF97 example and
Cargo/DSR remain unverified without the missing constellation siblings and tools.

## Adaptive finite gas and enthalpy-body windows

`FiniteGasThermalEnvironment::advance_adaptive` and
`advance_adaptive_with_convection` continue the original reciprocal thermal
transaction across a complete physical-time window. One native backward-Euler
interval is compared with two native half intervals; only the fine states and
their heat/entropy accounts are accepted. The existing first-order thermal
controller compares body enthalpy, liquid fraction, both temperatures and native
fixed-density gas internal energy. Formation references never scale time error.

Prescribed body and gas heats are total window energies delivered at constant
power. Cumulative source endpoints apportion them once through all accepted
intervals. Rejected comparisons supply no physical heat or elapsed time.
Source-domain failures propagate; only the existing exchange bracket/increment
admission can trigger subdivision. Every native interval spends the shared
attempt budget and retains all inner work limits. Final original body/gas
storage, radiation, source apportionment and summed absolute local energy
defects must pass a separate whole-window joule budget before publication.

Compact accepted history retains convection sources and original entropy
accounts without copying a gas catalog per sample. Fixed geometry, mass and
composition remain the native owners' assumptions. This is local time-error
control, not a global-error or unconditional-convergence certificate. The
`thermovisco_dynamics::finite_gas::window` regressions specify three-level
two-capacity decay, total-source ownership, external radiation, latent-only
accuracy control and atomic budget/source/cancellation refusal. They are
source/formatter checked; native execution remains unverified in this session.

## Multiple evaporating faces on one native liquid cell

The shared `liquid_gas::evaporation::multiple` transaction accepts distinct
measured Cartesian gas faces of one occupied VOF liquid cell. Every original
face passes native adjacency, exposure and per-face area admission; a repeated
liquid/gas face refuses. Canonical `(liquid, gas)` ordering makes input order
irrelevant. The input work ceiling counts faces, while physical inventory,
caloric inversion, donor allowance, mass selection and native receipts remain
one per liquid cell. Different cells retain the original common-atmosphere
species and caloric account.

Faces of a cell must retain the same original phase-law object and volatile
binding. Their coefficients remain separate in the original
`ParallelFilmEvaporation` owner, which sums actual per-face Stefan and finite
interface-kinetic fluxes at the temperature used by the phase solve. Averaging
mass-transfer coefficients before a nonlinear kinetic resistance is excluded.
The physical area is the sum of admitted face areas and may exceed one cell
face. A singleton retains its original law and arithmetic; zero-area groups
retain zero exposure.

Optional sensible heat uses the sum of face conductances against the one
liquid/gas temperature difference, with compatible original heat controls.
After that reciprocal heat stage, every phase film uses the actual heated
state. One summed native mass/enthalpy transaction refreshes material sources,
and one native flow step advances time. Existing momentum, conduction, reactive
and gas-plumbing compositions share this owner. No extra liquid inventory,
per-face gas reconstruction, phase heat or native timestep is introduced.

The original frozen-rate/time-splitting and fixed-volume common-gas assumptions
remain. This does not add spatial gas fields, overlapping sub-face geometry,
nucleation, interface exhaustion or overflow remapping. Four focused regressions
cover a two/three-face native corner, nonlinear heterogeneous kinetic rates in
both directions, actual post-heat rates, independent whole-system mass/energy,
one donor allowance, canonical ordering, and atomic geometry/source refusal.

## Shared finite gas with distinct enthalpy bodies

`finite_gas_thermal::bodies::FiniteGasBodiesEnvironment` advances one original
`IsochoricGasStorage` and a bounded set of independent native lumped bodies at
one shared gas endpoint. Each body retains its own caloric/phase curve, geometry,
convection law, internal heat and radiation enclosure. The scalar gas exchange
is matched to the sum of actual native body convective transfers; individual
transfers are never corrected or redistributed. Opposing transfers still solve
every body when the gas has exactly zero net heat. Canonical physical body IDs
make input ordering immaterial without imposing sequential reservoir updates.

`advance_with_wall_convection` also resolves each body's own unknown wall
temperature through the original single-body wall-feedback helper. Natural and
film-temperature-dependent material laws see that wall and the same shared gas
endpoint. Each body retains both queried and actual-endpoint source evaluations,
its temperature and boundary-work residuals, and its complete wall trial count.
The sum of absolute wall-law work mismatches is bounded by the shared source
tolerance, divided numerically between bodies without changing any heat flow.
No extra gas storage or sequential reservoir update is introduced. Focused
regressions cover independently solved nonlinear heat equations, the existing
single-body path, opposing latent transfers and late source/cancellation refusal.

Original fixed-density gas internal-energy and entropy increments are counted
once, including supported real-gas departures. Body storage uses each original
enthalpy and relative entropy law, including latent plateaus. The shared entropy
account sums each interface at its own body temperature and assigns the scalar
source mismatch at the shared gas temperature. One whole-account joule allowance
bounds the sum of absolute native body, gas and matching defects. Radiation
remains exchange with each explicitly prescribed external enclosure.

All endpoints are immutable proposals published together. Native source/phase
refusals, callback errors and cancellation retain body identity and publish no
partial update; external callback effects remain the caller's responsibility.
The body count and existing scalar/native iteration limits bound work. This is
a fixed-volume backward-Euler lumped operator, without direct body/body contact,
spatial gas transport or a replacement for existing spatial wall/piston networks.
Five native regression groups cover independent simultaneous-capacity solutions,
opposing latent transfers, gas-source/radiation ownership, original single-body
parity and transactional refusals, and dense-gas departure-energy conservation.
They are source/formatter checked; native execution remains unverified here.


## Simultaneous finite moving-gas molecular networks

`molecular_exchange::network::step_candidate` composes bounded ID-keyed moving
gas inventories and disjoint molecular apertures in one transaction. Every
aperture samples the original snapshot. Gross outgoing species are pooled per
gas and admitted against its original inventory before any incoming mass can
be credited. Each gas then receives one existing source-relative caloric and
momentum reconstruction. Link ordering is canonical by ID. A zero-area link
leaves its source and drift inactive; isolated gases retain their original state.

Apertures may have distinct prescribed velocities. Their molecular energy and
vector momentum are transformed analytically into an explicit common energy
frame. Actual species, total mass, momentum, reciprocal stream energy, rounded
pooled energy and final common/inertial-frame energy accounts must close before
all final endpoints are published. Source catalogs, final states and drift
support remain admitted by the original owners. Active apertures incident on
one gas must retain the same original catalog; a shared species identity cannot
substitute a different caloric source in another connected or disconnected part.

This is a frozen-rate, first-order fixed-volume interface operator. It does not
supply container motion, in-flight particle storage, other-wall forces, or a
closed-system passivity claim. The existing pair and pair-trajectory APIs remain
unchanged. Seven focused regression groups cover pair parity, three-gas balance,
canonical ordering, shared withdrawal, orientation, Galilean/caloric-reference
invariance, newly arriving species, topology refusal and terminal cancellation.
Native crate execution remains pending after the dependency build exhausted disk
space. Following executor recovery, recovered files match the original
implementation after Rustfmt; formatting and whitespace checks pass.

## Finite gas calorics in spherical cavitation sound

With the opt-in `caloric-cavitation` feature,
`cavitation::bubble::caloric::MaterialCaloricBubble` binds the original
liquid source to the native `fs-flux` caloric radial owner. The original finite
noncondensable gas supplies mass, composition, temperature, pressure and
temperature-dependent heat capacity. Explicit signed heat and pressure
histories drive the same joint gas/radius trajectory. The adapter does not
replace gas calorics with a fitted exponent or evolve a second state.

Native borrowing proposals permit admission of compact sources and observers
before publishing gas, radius, wall velocity and physical time together. A
failed observer, source export or final cancellation drops that proposal.
The endpoint exports native `Q = V'` and analytic `Q' = V''`, including the
heat-driven pressure derivative in Keller--Miksis. The existing observer uses
the original liquid rho/c, signed monopole pressure and causal arrival time;
an observer from a different medium refuses. Acoustic observation adds no
second radiation damping. Exported histories can use the existing coherent
receiver-clock renderer under its bandwidth and compactness requirements.

Liquid coefficients, liquid temperature and equilibrium vapor pressure remain
frozen. Gas heat is an explicit external port; it does not silently heat a
finite liquid or supply latent heat. This spatially uniform ideal-gas model
does not resolve gas thermal boundary layers, vapor mass kinetics, shock
collapse, erosion, translating bubbles or interacting clouds. The gas/liquid
density, vapor-pressure fraction and gas wall-Mach limits are the native
owner's additional admission conditions.

`tests/caloric_bubble.rs` specifies both native radial equations through the
adapter, material/query retention, caloric heat-dependent motion and sound,
the analytic heat-pressure-rate source, independent source caloric energy,
causal microphone pressure, and late observer/cancellation rollback with
deterministic retry. These are manufactured model and composition regressions,
not measured validation. All seven caloric adapter tests now pass in a
focused original-source build, including exact native/adapted trajectory
agreement and late rollback. The example below completed 100 samples in each
branch; its accumulated original-law balance residuals were 1.765e-18 J
(adiabatic) and 6.657e-18 J (heat exchange). These are synthetic model results,
not full workspace/DSR, experimental accuracy or performance measurements.

Running `cargo run -p fs-couple --features caloric-cavitation --example
caloric_cavitation -- crates/fs-material/examples/material_gas_spring/mixture-a.tsv
2500 100` compares the same finite caller-authored gas under adiabatic evolution and
explicit `h A (T_liquid-T_gas)` exchange during a smooth compression pulse.
The liquid properties come from the original IF97 water owner. Output retains
radius, gas temperature/pressure, causal microphone pressure, cumulative heat
and gas work, plus an independently source-integrated final energy residual.
The example's gas and heat-transfer coefficient are explicitly synthetic;
neither this example nor the lumped heat law supplies measured bubble accuracy.


## Conservative physical backbone sensitivities

`phs_orbit::PreparedPhsBackboneSensitivity` binds the original `PortHamiltonian`
and constant physical port bias to fs-orbit's conservative rectangular derivative
owner. `new` discovers the complete storage-owned Hessian; `new_analytic` accepts
an explicitly authored complete Hessian. Every original state and Fourier
coefficient remains, with an arbitrary original-coordinate fundamental-amplitude
anchor and a fixed complete-waveform phase reference. The original `s I` port
gives the exact frequency column. No loss, physical balance row, or material state
is deleted.

`solve_parameter` combines complete storage-gradient, skew interconnection,
symmetric resistance, reciprocal port-map and constant-bias directions through
fs-phs's compensated physical parameter owner:
`(J-R)e_p+(J_p-R_p)e+G_p u+G u_p`. The cached full state tangent supplies the
separate unknown state chain. A zero source direction and unit amplitude
derivative give the nonlinear backbone slope; physical parameters can instead
hold amplitude fixed. One QR factor serves all directions without perturbing or
refitting the source.

Passivity alone does not establish a conservative periodic family. The complete
base and differentiated original balances, phase and amplitude must pass their
independent two-grid numerical gates; a dissipative direction that destroys the
family refuses. Invalid or incomplete physical directions, unavailable/refused
tangents, resource overruns and cancellation publish no derivative and leave a
prepared factor reusable. Results are local fixed-order numerical derivatives,
not proofs of infinite-tail resolution, stability, or branch identity.

Three focused tests cover Duffing amplitude/material derivatives against
neighboring original nonlinear solves and an independently differentiated
continuous energy-period quadrature; a four-state elastically coupled nonlinear
system with bias-generated even harmonics and simultaneous storage/J/G/bias
changes; and tangent/resource admission plus incompatible loss, incomplete/invalid
direction, cancellation and reuse. They passed in the bounded native lane using
complete original fs-orbit and fs-phs, actual fs-la/FFT numerical modules, and
original PHS-adapter/tangent source. This is not a full fs-couple
Cargo/workspace/DSR gate.


## Coherent tensor-solid harmonic heating on the original thermal residual

`tensor_material::solid::thermal::harmonic::heat_assembly_candidate` observes
one `ViscoHarmonicMaterialSpectrumObserver` over the actual thermal interval
and deposits each original element's Kelvin/Maxwell heat once into its explicit
whole-element caloric receiver. Carrier, beat and equal-frequency interference
remain in the original quadrature account. The original `HeatAssembly` and
`cell_driver_candidate` own nonlinear transport, contacts, reservoirs,
radiation, source evaluation, thermal residuals and final acceptance. Repeated
thermal trials reuse the same material observation.

Each element retains its reference mass/volume, source identity and exact caloric
law. Complete nonthermal queries and selected source claims remain pinned.
Original tensor, relaxation, Kelvin and density source laws must be independent
of the evolving temperature: source axes are checked even when two sampled
values coincide. The explicit fully-solid temperature interval and original
source selections are re-admitted at actual proposed endpoints. Unsupported
source changes, incomplete/duplicate routing, strain excess, changed inventory,
late provider refusal and cancellation publish no candidate or partial state.

The report retains the original material field, element heat, native thermal
audit, source sample, heat-handoff residual, other external heat and separate
caloric and constitutive-plus-caloric balances. Stored spring energy, signed
stress work and numerical defects are never deposited as heat. The combined
account encloses constitutive and caloric storage only; it does not establish
kinetic-energy or external-force-work balance.

Exact interval heat supplies an interval-mean source to native backward Euler.
Time refinement is required for beat-driven thermal transport; fixed mechanical
coefficients, stationary reference geometry and explicit source-temperature
independence bound this model. Startup memory, thermal expansion and changing
mechanical coefficients require their existing transient/coupled owners.

Four integration tests in
`tests/tensor_material_source/thermal_relaxation/harmonic.rs` cover independent
real Maxwell/Kelvin heat-work-storage quadrature, two-cell Fourier conduction
across consecutive beat windows, temperature-dependent and flat-but-authored
source refusal, original strain/receiver/inventory admission, late provider
failure, cancellation and retry. Formatting and source/API checks passed;
these solid/thermal integration tests were not executed in the available
runtime/constellation environment, so no full coupled-runtime gate is claimed.


### Finite material-wall heating of a freely recoiling effusive vessel

`rocket_blowdown::effusion::thermal` composes the original implicit finite
wall/gas thermal owner with the original finite-mass molecular recoil owner.
`GasRocketState::advance_effusion_with_wall` and
`advance_effusion_with_wall_convection` apply both stages over the same requested
physical duration as a first-order heat-then-escape split. The actual heated gas
drives emitted species, molecular energy/momentum and subsequent recoil. The
callback variant evaluates the original state-dependent heat-transfer law;
the caller qualifies its molecular or continuum applicability.

Wall mass is already included in vehicle dry mass and is never added again.
Constant-density wall calorics and consistent fixed volume are required;
unmodeled expansion work refuses. Both material and vessel endpoints remain
tentative until all original source, Biot, caloric, withdrawal and mechanical
gates pass. The combined reference-free energy budget bounds the sum of absolute
stage defects plus native mass-remapping and gravity-history reserves. Opposite
defects cannot cancel through this gate. Original thermal and recoil reports,
external thermal heat, mass/impulse and mechanical work remain distinct.
Rotation, radiation momentum, structural integrity and wall-expansion mechanics
remain outside this adapter.

Six real-owner G0/G1/G3/G4 regressions specify analytic finite-capacity heating
and changed molecular rates, insulated parity, closed-port heat, radiation and
external loads, frame shifts, state-dependent h, atomic late/cancellation
refusals, and inconsistent mass/volume/expansion admission. Pinned Rustfmt and
whitespace checks pass. Native fs-couple integration is unexecuted: DSR/RCH are
unavailable and shared build storage cannot support the required dependencies.
Related broader porous/capillary transport scope remains open.

## Pressure-dependent hydraulic ports and passive reservoirs

`cavitation::duct::pockets::CavitatingDuctPockets` accepts a smooth deterministic
`Q(time_s, absolute_port_pressure_pa)` boundary through
`advance_with_port_flow_law[_and_gas_heat]`. Initial contact selection, every
native RK stage, localized rupture/collapse event and final observation sample
the actual current port pressure. The existing DG liquid, active-set pocket
constraints and original gas inventories keep one accepted trajectory; no
pressure clamp or separately advanced boundary state is introduced.

The declared endpoint is positive INTO the duct at either end. Actual signed
transported volume and acoustic pressure work are integrated as `integral Q dt`
and `integral p' Q dt`. Reports distinguish initial, endpoint and mean flow;
endpoint flow times duration is not used for a varying drive. Existing held
flow retains its original adaptive dimensions/tolerances and exact `Q*dt`
volume. Feedback adds one scalar volume coordinate, and a reservoir adds one
loss coordinate, without per-stage allocation. Gas caloric stages retain the
same conservative energy coordinate through the source-owned sub-quantum
sampler; original-source endpoint reconstruction remains strict.

`DuctReservoirResistance` provides the explicit quasi-steady law
`Q=(p_res-p_port)/R` for the same frozen liquid with positive supplied R.
Signed reverse flow is retained. Reservoir pressure must stay above the
original vapor pressure and within the retained liquid's strain bound.
The original liquid/cavity/Mach/domain limits remain active at every stage.

The same trajectory independently integrates `integral R Q^2 dt`. Before
publication, reservoir work, native acoustic work and this nonnegative loss
must close their balance. Reference-pressure transport work is reported
separately from the acoustic Hamiltonian. External resistance heat is a
reported port; it is not silently deposited in the isothermal liquid. Callback
failure, invalid physical state, work exhaustion or cancellation preserves
all accepted liquid, pocket and gas state. Callbacks must be side-effect-free
and valid at rejected trial states; their bandwidth needs resolved time steps,
and known discontinuities must be split at an accepted endpoint.

The IF97-water `cavitating_reservoir` example switches between withdrawal and
supply reservoirs and exports physical pressure/flow, cavity volume and events,
transported volume, reservoir/native work, resistance loss and both residuals.
Run `cargo run -p fs-couple --features cavitation --example cavitating_reservoir`.
Its geometry, initial field and resistance are demonstration inputs, not a
measured valve calibration. No valve inertia, finite reservoir thermal state,
vapor mass/latent heat, nucleation delay, pump curve or moving-mesh model is
claimed.

Focused regressions specify an independent acoustic polynomial with Q(t,p),
a distributed reservoir decay and lumped RC limit, signed work, cavity
rupture/collapse, variable-Cv gas heat plus feedback, and transactional failure.
The changed source passed rustfmt/whitespace checks and review against the
original DG/event APIs and balance identities. These DG tests and the new
example have NOT been executed in the constrained dependency environment;
no hydraulic runtime, full-workspace or performance result is claimed.

## Source-driven liquid thermal diffusion in bubbly acoustics

`ThermalDiluteBubblyLiquid::new_with_conducting_liquid` and
`ThermalPolydisperseBubblyLiquid::new_with_conducting_liquid` prepare the native
coupled spherical gas/liquid diffusion model from the retained liquid source's
`ThermalFluidTransport` view. In particular, `SaturatedWaterCavitation` delegates
that view to its original IF97 transport state. Temperature, density and
viscosity must match the original acoustic/cavitation snapshot, including for
empty populations. Each occupied class keeps its original gas, actual radius,
number density and KM tangent. A class split with unchanged total number
weight leaves dispersion and all loss channels invariant.

The same native bubble response drives susceptibility, coherent scattering,
Foldy dispersion, source amplitude admission and existing scalar/duct-wave
operations. Thermal damping is not fitted and acoustic radiation is not added
a second time. Prepared models are stored once in original class order;
queries allocate nothing, query no material sources and do no equilibrium
root solves. `thermal_bubbles()` now exposes a zero-allocation exact-size
iterator over the retained base gas preparations, rather than a borrowed slice.

`ConductingLiquidPopulationLimits` bounds liquid thermal penetration
`sqrt(2 D_l/omega)` relative to the original population's mean spacing and
bounds actual peak interface temperature relative to its equilibrium value.
The permitted ceilings are 0.1 and 0.05 respectively. These conservative model
admissions do not prove spatial separation for an arbitrary unseen cloud.
Frozen liquid mechanics, isolated thermal fields, small acoustic amplitude,
noncondensable gas, negligible latent heat and the original scattering bounds
remain required. A failed preparation or query leaves accepted owners usable.

Diagnostics report combined conductive loss and its independently computed
gas and surrounding-liquid contributions, separate from liquid viscosity.
`BubblyAbsorption` now includes `liquid_thermal_w_m3`; each channel uses the
actual local peak pressure squared. Coherent scattering attenuation is not
converted into heat. The response additionally exposes the largest liquid
diffusion parameter and admitted fractional interface-temperature amplitude.
The existing `cavitation` example compares fixed-temperature and conducting
interfaces using original dry-air and IF97-water material evaluations.

All seven tests in the actual thermal-population test file passed in a focused
original-source build: independent modal/interface reconstruction, original
scalar duct impedance, number/class weighting, heat-channel closure, IF97
source retention, and material/thermal-domain/cancellation refusals. The broad
DG/duct-chain dependency closure and full workspace/DSR were not executed.

## Frequency-qualified modal radiation reaction

`vibroacoustic::ModalRadiationLoading` replaces the unqualified exterior matrix
in `VibroacousticModel::try_new`. A direct sample carries its exact angular
frequency; `from_bem_bake` retains the native BEM/Rayleigh producer's complete
frequency grid, matrices and declaration identity. A loaded FRF selects its
original row by binary search. An absent frequency refuses on all full,
truncated and convergence solve paths; no frozen loading, interpolation or
extrapolation is inferred. A sweep reuses the original structural/cavity bases.

The shared `fs-modal` physical-port owner admits the COMPLETE Hermitian work
form at preparation, with its existing dense-work and cancellation bounds.
Positive diagonal resistance alone is insufficient. Nonfinite or active mutual
directions refuse without symmetrization, clipping or loss floors. Original
mutual reaction and exact loss nullspaces remain in the coupled equations.
Radiated power contracts the Hermitian part before multiplying fields, using
double-double accumulation; large reactive exchange cannot become fictitious
outgoing work. Truncation uses the same principal matrix and original stride,
without allocating another velocity or impedance matrix.

This is sampled passivity, not a causal realization or an inter-sample proof.
The producer still owns geometry, medium, modal normalization and spatial
convergence. Focused G0/G1/G3 regressions specify spherical added mass/radiation,
coherent mutual work and nullspaces, frequency refusal, active complex
directions and cancellation. The existing BEM casebook consumes one complete
loading bake across its physical sweep and checks the original loaded equations.


## Adaptive finite-wall effusion and vessel recoil

`GasRocketState::advance_effusion_with_wall_adaptive` and its convection-callback variant evolve the original implicit finite wall/gas heat step and original effusive recoil step over one physical window. Full versus two-half-step estimates include original gas temperature, pressure, source-qualified species masses, drift-free vessel motion, finite wall temperature/storage, and exported momentum/energy. Only the accepted two half steps enter state, history and totals. Thermal and recoil stages share one physical duration; the material wall remains part of the existing dry mass. The estimator controls the original first-order split's local temporal error, not constitutive error or a certified global bound.

Exported moments use one declared inertial frame. Energy comparisons remove fixed original species formation baselines analytically through their caloric integrals, then include the native kinetic transport and source correction once. Wall error uses the endpoint enthalpy difference before multiplication, with relative scaling from storage increments. The original phase owner still rejects changes below its represented-enthalpy resolution. Whole-window prescribed heat is interpreted as held power; actual accepted allocations are summed and any difference from authored totals consumes the absolute energy budget. The common temporal helper retains exact representable half steps, progress, final time and bounded work/history.

Fixed window budgets include actual species/total mass closure, accumulated absolute native defects, initial mass remapping, momentum/energy reserves and heat-allocation error. Opposite signed residuals cannot enlarge a budget, and no budget scales with step count. Only the original typed gross-withdrawal refusal can cause a physical retry; original source, caloric, mechanical, callback and cancellation failures remain typed refusals with no endpoint publication. Original convection callbacks may observe rejected trials and must be repeatable.

Six focused regression groups cover independent continuous monatomic heat/effusion/recoil equations, exact accepted-history replay and heat allocation, inertial-frame transformation of exports, emitted-energy-driven refinement and fixed budgets, represented wall-reference shifts, and retry/work/minimum-step/source/cancellation behavior. The fixed-budget tests include a stationary closed-outlet vessel with an admitted initial native mass mismatch. Pinned formatting and whitespace checks pass. These full-owner regressions are authored but not executed under the current dependency-build capacity limit.


## Relaxed solid bias in coherent thermal handoff

`tensor_material::solid::thermal::harmonic::heat_assembly_candidate_with_bias`
shares the original unbiased candidate and caloric solve. The complete represented
static displacement reaches original constitutive and compatible B-bar strain
admission before caloric trials. Relaxed bias changes equilibrium storage and
stress work while leaving independently integrated Kelvin/Maxwell heat unchanged;
only that original physical AC heat is deposited. Source-temperature independence,
reference mass/volume, exact caloric ownership, nonlinear transport, cancellation
and constitutive/caloric energy gates remain unchanged. Neither static force
balance, geometric prestress, thermal expansion nor startup memory is inferred.

## Reference-only glass optics on the accepted thermal timeline

`layered_wave_heating::SchottThermoOpticalSource` supplies the native
monochromatic and broadband thermo-optical heating owners with an immutable
original dielectric spectrum plus the shared SCHOTT absolute-index law.
Each original thermal cell's accepted temperature refreshes its separately
selected material/permeability requirements on the same immutable card. The
reference spectrum stays at its original T0, even when its claims have no hot
support. All supplied current-state requirements, source identities, coordinate
semantics, interpolation gaps and intersected wavelength/temperature bounds
remain binding. Reference data/law are Arc-shared across cells and endpoints;
no fitted spectrum or accumulated temperature increment replaces them.

The explicit lossless mode refuses original absorption. The separate reference-
absorption mode freezes original alpha(lambda) while the native full complex
Maxwell response changes refraction, standing waves and positive local volume
loss. Existing native owners sum physical spectral-line watts once, preserve
two-face coherence within each line, solve heat once, and admit all endpoint
optical states before publication. Original work bounds and cancellation apply.
This is the existing quasi-steady explicit optical/thermal split with fixed
geometry and caller-owned interval convergence, not a supplied absorption(T),
transient EM storage, expansion, thermal-emission or implicit feedback law.

Three focused regressions in `tests/layered_wave_heating/thermo_optical/schott.rs`
specify exact reference retention and source refusal, nonuniform native cell
heat/cool refresh with an independent two-cell caloric solution, and broadband
complex Airy absorption with endpoint-domain atomicity. Source/API, independent
scalar formula and pinned formatter checks were completed; these native
integration tests remain unexecuted in this environment because the dependency
build was killed before the physics crates were reached.

## Evaporation into distinct native gas regions

`liquid_gas::evaporation::multiple::regions` binds original native phase or
phase/sensible-heat faces to explicit finite gas parcel IDs. Each region uses
its actual pressure, temperature, vapor inventory and original gas calorics;
its liquid cells share the original chamber donor account. All regions prepare
against the same unchanged liquid. One combined native mass/enthalpy transaction
refreshes the original material, followed by exactly one native flow step.
The existing single-atmosphere and regional paths share source preparation and
native caloric accounting; momentum re-admission and flow callbacks remain with
their existing single-atmosphere owners.

All faces of one liquid cell and each exposed native gas coordinate require one
region. The native pure liquid retains one volatile identity and caloric datum
across all regions. Distinct faces of a cell keep their nonlinear parallel-film
rates and one inventory. Canonical region/cell ordering is independent of input
face enumeration. The total original face count bounds work. Per-region and
whole-interval absolute mass/energy gates prevent opposing source defects from
cancelling. Unused gas nodes remain unchanged; source refusal, native refusal
or cancellation publishes neither the liquid nor any regional gas.

The heat variant runs each region's original reciprocal sensible-heat then
phase stage before the one native batch. Latent/carried vapor energy and
sensible heat are each counted once. This is first-order splitting with
prescribed regional membership and fixed external gas volumes. It adds no
regional pressure traction, inter-region gas transport, moving gas momentum,
volume work, resolved gas CFD or exhausted-cell/topology remapping.

Four focused native regression groups cover different humidity at identical
liquid temperature, signed per-region species/energy accounts, reciprocal heat,
permutation invariance, one clock, inactive nodes, invalid regional mappings,
source refusal and retry. They are added but unexecuted under the existing
pinned dependency/resource build block; formatter checks do not establish a
native compile or runtime result.

## Continuous multi-tone forces on retained tensor solids

PreparedTensorBody::from_mass_matrix_tones and
PreparedTensorSolidBody::new_tones extend the existing complete tensor-memory
owner with continuous cosine/sine force records. The physical state, exact
source-time groups, original consistent mass, quadrature points and per-element
Kelvin/Maxwell heat remain unchanged. Source records carry no material energy.
All tone interference enters the same native PreparedQuadraticFlow trajectory
and the original positive PreparedSourceLoss channels.

Each interval supplies coefficients relative to its own start. The native
quadratic integral of F(t) dot v(t) is checked independently against cosine/sine
velocity moments from a one-way rotating observer. The solid adapter expands
those moments through its original displacement basis and checks full nodal
work before publishing the candidate. Fixed support forces do no work.
advance_force_phase rotates unchanged source coefficients without advancing
mechanical history. Source frequency order, including zero or repeated tones,
is retained; no damping fit, mode deletion or independent-tone heat sum occurs.

Each frequency must lie in every original TensorProny law's admitted band.
The original endpoint strain admission remains; neither full-interval material
applicability nor the transient spectrum's containment in that band is claimed.
Force records, moment observers and all heat caches enter the existing explicit
state, preparation byte/work and warm-work budgets. Cancellation, malformed
sources, exceeded budgets and failed work accounts preserve accepted history.

Focused regressions specify nonproportional two-coordinate tensor response
with full consistent mass against independent original-coordinate RK4, both
weighted velocity moments, signed source work and separate positive heat;
zero-frequency held-force parity; source-band, capacity and atomic refusals;
and original two-element nodal loading with full/half-clock agreement.
Pinned rustfmt parsing completed. Native Cargo compilation and these regressions
remain unexecuted in the current development environment.


## Inertial phase-changing water chamber

`cavitation::water_cell::piston::WaterPiston` couples the original finite IF97 water inventory to a constant-area inertial wall. Total liquid-plus-vapor volume is the sole position coordinate; original water mass and conserved internal energy remain in `EquilibriumWaterCell`. Each bounded mechanical trial calls its existing volume/heat proposal. The exact same trapezoidal endpoint pressure supplies water pressure work and the shared native piston momentum kernel, with midpoint wall kinematics. The existing gas piston delegates to this extracted mechanical kernel without changing its arithmetic or source law.

Held exterior absolute pressure, additional mechanical force and integrated water heat are explicit ports. Accepted steps expose vapor mass/volume changes, pressure/force work, wall kinetic energy, kinematic/impulse defects, and independently admitted conserved-total and original-source energy residuals. The ideal fixed support supplies reaction impulse with zero motion/work and can only attach to a wall already at rest. Proposals publish water, wall velocity, phase inventory and physical clock together after all checks and final cancellation; dropping a proposal preserves all accepted history.

This is a uniform, sealed, rapidly equilibrating chamber, not a resolved spherical bubble or water-hammer model. Nucleation delay, noncondensable gas, phase-front geometry, wall impacts and resolved thermal gradients remain outside its scope. Numerical interval, source bracket, temperature and volume bounds remain explicit. Conserved energy alone does not certify timing: focused tests additionally compare oscillation phase with the original homogeneous-equilibrium acoustic tangent and demonstrate convergence under interval refinement.

Fifteen focused original-source tests passed for the new mechanical and retained cell/volume paths. `cavitating_water_piston` executes a 1 g inventory under an 80 kPa then 140 kPa load, showing cavitation and subsequent complete condensation without prescribing its volume. The example exposes independent cumulative energy and accepts an integer time-refinement factor. The broader gas-piston graph was source-reviewed for the shared-kernel extraction; it was not executed in the focused dependency build.

## Molecular porous refraction on the accepted thermal timeline

`porous_heating::feedback::IsobaricRelaxingPorousGas` retains the original
isobaric gas specification, a complete `RelaxingMedium` composition and its
`RelaxationAbsorption` law. The reference equilibrium state must match that gas
law exactly. Both owners admit the declared thermal support and every later
native cell temperature; no extrapolation or independent absorption(T) fit is
introduced. Each changed cell refreshes density, sound speed, viscosity,
conductivity, molecular heat capacities and collision times together once.

The existing single-frequency, coherent spectral and covariance/independent
angular-mode feedback owners consume that same immutable state through
`PorousTemperatureFluid::prepare_layer`. The default preserves the original
equilibrium gas/liquid path. The molecular adapter uses the native finite
`PorousLayer` response, whose shared source ownership permits cache movement
and cloning without self-references or per-bin reservoir copies. Original
porous source support, conserved tangential wavenumber and full-area
pressure/volume-flow phase conventions remain binding.

Coupled pore-thermal/molecular compression and pore viscosity supply the
existing positive interior heat integrals. They are not independent added
attenuation channels. Every native explicit step deposits its accepted START
heat once; changed endpoint states prepare the next response and are admitted
before publication. Unchanged/latent-plateau cells reuse their original state.
Spectral bins share one source evaluation per changed cell; each frequency
retains its original coherent stack solve. The fixed exterior wave reference,
authored thermal mass/caloric law and separate passive backing-load account
remain unchanged. This is an isobaric, frozen-background-per-step quasisteady
model, not finite pore mass transport or transient molecular-memory evolution.

Four focused tests in `tests/porous_heating_feedback/molecular.rs` compare
oblique temperature-gradient fields and cell heat with the original material
and duct owners, selective accepted-state refresh with start-step energy,
shared spectral states and equivalent covariance heating, exact empty-reservoir
recovery, and source-domain/cancellation refusal without candidate publication.
They are authored and syntax checked; native execution remains blocked by the
dependency compilation resource limit.

## Finite-cell heat from adhesive rolling separation

`friction_heat::adhesive::snap::rolling::JkrRollingSnapHeat` owns an original
`SourceJkrRollingEvents` cycle and pins the original ordered material caloric
receivers. It regenerates bounded, explicitly ordered normal and angular
segments from sealed history. Original normal snap and released rolling spring
storage keep separate source reports and compensated transfer ledgers; their
declared surface shares enter actual finite caloric cells, with any remainder
assigned to a named outside destination. Prior source releases are never
re-deposited when a binding is created. Signed parameter work and accumulated
plastic rolling loss are not silently converted into snap heat.

Preparation stages the complete lifecycle, release ledgers and thermal cells.
Versioned source stamps prevent stale replay; caloric, source, budget or
cancellation refusal leaves every accepted owner unchanged. Detached rotation
creates no spring preload and recontact retains prior release/loss accounts.
The composition applies its releases at the current thermal time under the
original frozen-source approximation. It introduces no conducting contact link,
intervening thermal transport or temperature feedback; separate calls at actual
event times are required when transport between releases matters. Focused
cases in `tests/adhesive_snap_heat/rolling.rs` specify analytical spring release,
actual cell/outside energy closure, fresh recontact, reversible area work with
zero snap heat, and atomic caloric/cancellation refusal and retry.

## Regional evaporation coupled to original gas transport

`liquid_gas::evaporation::transport::regions` runs the existing complete
`GasTransportPaths` network once, then supplies its actual transported parcels
to the regional native phase or phase/sensible-heat owner. Different interfaces
therefore see their bound region's transported composition, temperature and
pressure. Arrivals can enter the following phase stage, while the original
transport snapshot still forbids funding outgoing paths with those arrivals.
All gas nodes retain their identities and fixed external volumes, including
nodes without liquid interfaces; the native liquid advances exactly once.

The actual final complete gas map must pass the original transport endpoint
driving-force and source gates. A phase/heat-induced pressure or concentration
reversal refuses every private gas node and native liquid together. The
single-atmosphere and regional compositions share the original residual gate:
absolute transport defects add to absolute native phase/heat/flow defects and
are bounded by the complete interval's explicit allowances. Internal partial
enthalpy, reciprocal sensible heat and latent subaccounts are not reapplied.
Transport entropy remains a statement about its original gas-only stage.

Both phase-only and sensible-heat APIs retain original source identities and
cancellation, first-order transport/heat/phase/flow splitting, one gas region per
liquid cell and prescribed regional membership. Refine the whole interval.
This composition adds no resolved gas CFD, regional gas traction, pressure
volume work, reaction stage, finite-wall/conduction stage or topology remap.
Three focused regressions exercise real nozzle arrivals and local phase/heat
feedback, independently reconstructed whole-system mass/species/energy, final
pressure reversal with atomic retry, empty-path equivalence and cancellation.
They are added but unexecuted under the existing dependency/resource build block.


## Vapor phase dynamics and displaced-liquid acoustic observation

`cavitation::bubble::vapor::MaterialVaporBubble` binds the original saturated-water source and native phase-changing radial owner to the existing compact acoustic observer. Evaporation and condensation make the geometric interface speed differ from liquid velocity. The adapter therefore transmits `4 pi R^2 U_l` and its native analytic derivative as the monopole source. It does not substitute geometric bubble-volume rate. Density, sound speed and causal propagation delay come from the same retained bulk-water source.

Native mass, temperature, radius, liquid velocity and their accepted clock remain a single history. Compact-source and microphone admission occur while the complete native advance is still a borrowing proposal. A wrong acoustic medium, late compactness/pressure refusal or final cancellation preserves the entire original physical state. The common compact-source formatter is shared by the existing gas-bubble adapters; observation adds neither radiation damping nor a second evolution of the vapor inventory.

Five new integration tests passed for identical direct/coupled phase trajectories, actual liquid displacement during evaporation, heat-driven mass/acoustic feedback, receiver-medium and late-observer rejection, and cancellation after source admission. Existing gas caloric and conducting-bubble integration tests remain passing. The complete focused source run passed 146 selected tests across material, native cavitation and coupled owners; it was not a full-workspace build.

The `vapor_cavitation` example compares zero heat with an explicitly authored smooth 1 mW pulse under the same pressure cycle, emitting 320 samples of interface/liquid velocities, vapor mass transfer, source flow, causal receiver pressure, heat and independently reconstructed cumulative balances. Both cases include evaporation and condensation. The heat pulse is an external thermal port, not a liquid-conduction correlation. Applicability remains the native incompressible, spherical, locally saturated pure-vapor model plus the original compact-observer limits. Uniform vapor temperature, continuum behavior and rapid phase equilibration remain caller assumptions; the mechanical Mach checks alone do not establish them.

## Hereditary material in cavity/exterior-loaded modal response

`vibroacoustic::ModalStructuralRelaxation` admits the original `fs-phs`
Maxwell arms through the existing BEM material-memory owner. Structural
frequencies describe relaxed equilibrium stiffness; every original signed
strain row adds its excess storage and loss exactly once. A source-owned plate
adapter preserves the original material/property selection and intersects all
regional frequency bands. Full, truncated and convergence solves enforce those
bands at the actual forcing frequency. A scalar structural loss and duplicate
attachment refuse; a separately identified Kelvin mechanism may coexist.

The native `fs-material` Prony spectral weights assemble the complete complex
modal stiffness, including mutual terms, without extra factorization dimensions,
fitted damping ratios or equal-time/rank merging. The original row prefix is
used under truncation. Branch channel squares independently provide positive
heat and cycle-mean spring storage. Heat is included once in total structural
loss while cavity, boundary and exterior-radiation work remain distinct.
The source owner must establish the actual relaxed basis and material validity;
equal dimensions do not establish identity. This harmonic composition adds no
time stepping, nonlinear constitutive or measured-audio claim.

Focused G1/G3 tests compare the coupled two-mode/cavity/radiation response with
independent spring-dashpot and cavity-elimination equations, including injected
volume flow, branch heat/storage, loss nullspace, principal-block truncation,
frequency support, cancellation and zero-arm parity. Native numerical-slice
verification compiled the original `fs-math` and exact new spectral contraction
bodies: ten tests passed, comprising six original Prony spectrum controls and
four contraction controls. Full `fs-couple` Cargo integration tests were not run
in the disk-constrained session; syntax/static checks and independent complex
work-balance calculations are distinct from that missing execution.

## Nonlinear multi-terminal electrothermal feedback

`step_response_terminal_candidate` and
`step_axisymmetric_response_terminal_candidate` place held voltage, conserved
current and mixed electrode controls inside the original nonlinear electrical
and finite enthalpy residual. Floating equipotential electrodes can transfer
current locally while their total current remains zero. Each thermal trial
refreshes the original bulk/contact temperature and declared phase context;
every constitutive trial solves the actual grouped-electrode operating state
once, retaining the original source laws, local fields and contact gates.

`ConductionResponseThermalStep::terminals` retains the accepted ordered terminal
potentials, reactions and independent source-work report. Current-source work
uses authored currents; voltage-source work uses native reactions and retained
physical voltage differences. Original swept/planar cell heat and explicitly
partitioned contact heat enter the same enthalpy account once. Original Fourier
links, fixed geometry and field-independent caloric/Fourier admission remain
binding. Every native field solve consumes the existing circuit budget,
including discarded constitutive trials and fresh final admission; failures
publish neither accepted owner.

Focused G1/G3/G4 regression definitions cover E/J/T-dependent materials in
planar and axisymmetric bodies, floating-face current transfer,
voltage/current/mixed controls, large common references, reordered continuation,
exact solve budgets, late cancellation/retry, nonlinear contact heat partitions
and latent-to-sensible phase continuation. Pinned syntax checks and independent
reference balances were performed; these new native integration regressions
remain unexecuted in the current dependency environment.

## Axisymmetric nonlinear electrothermal feedback

`step_axisymmetric_response_candidate` and
`step_axisymmetric_response_contact_candidate` advance the original finite
enthalpy history with `AxisymmetricConductionResponse`. Each thermal residual
trial refreshes the complete original temperature/phase sources, then closes
their electrical-field and current-density coordinates at the actual drive.
Voltage, total current, current-density traces and regulated power use the
native nonlinear electrical iteration; no frozen field is rescaled.

Explicit electrically varying contacts close their original voltage-jump or
normal-current coordinates under the native uniform-trace requirement. The
existing contact bindings retain both receiving bodies, their actual thermal
and phase states, ordered source identities, fixed external context and authored
heat fractions. Exact swept bulk and surface losses enter the original
enthalpy residual once, alongside the original Fourier links. Independent
terminal work must close with total deposited heat; constitutive defects are
not thermal sources.

Every native field solve, including discarded nonlinear iterates and fresh
final admission, consumes the existing circuit budget. Accepted continuation
returns both updated owners and their common report. Source/domain refusal,
budget exhaustion, compliance failure and cancellation publish neither owner.
Caloric response/applicability and Fourier conductivity dependent on electrical
coordinates remain refused; the scalar memoryless model retains fixed geometry
and the original constant-density caloric charts.

Focused regression definitions cover E/J-dependent bulk laws, all four drives,
sector scaling and unequal swept cell heat, nonlinear contact/temperature
feedback, two-sided surface heat, latent then sensible heating, source pins,
continuation, budget accounting and atomic retry. Pinned syntax checks and
independent body/contact/phase balance calculations were performed; the new
native coupled integration regressions remain unexecuted in the current
dependency environment.

## Axisymmetric multi-terminal electrothermal controls

`axisymmetric::terminals::AxisymmetricTerminalDrive` binds complete terminal
voltages, conserved terminal currents, or mixed voltage/current partitions to
`ConductionThermalNetwork::step_axisymmetric_terminal_candidate`. Floating
equipotential electrodes retain zero net current while allowing finite current
transfer between their physical faces. Each thermal trial refreshes the original
bulk and contact sources and rebuilds the native N-1 conductance basis; every
basis solve consumes the existing circuit budget, with preflight admission
before any basis solve starts.

The existing enthalpy residual, Fourier links and explicit two-body contact heat
partitions remain the sole thermal account. Source work uses authored currents
on current-controlled terminals and native reactions on voltage-controlled
terminals, each once. Retained relative/reference potential differences are
combined in double-double precision before rounding, so a large common voltage
reference does not erase a finite physical drop. Accepted continuation retains
the refreshed electrical owner and thermal history together; invalid controls,
budget exhaustion and cancellation leave both inputs reusable.

Focused regression definitions cover floating-disk transfer, three sweep angles,
temperature feedback, unequal caloric capacities, Fourier exchange, all three
control modes, large common references, terminal reordering, finite-contact
heating, exact basis budgets and deterministic retry. These new native coupled
integration regressions have not been executed in the current dependency
environment.

## Axisymmetric spatial electrothermal feedback

`circuit_enthalpy::conduction::axisymmetric` advances the original finite
enthalpy history with source-resolved meridional conduction at the same
thermal endpoint. Native swept-cell Joule heat enters that residual once;
the original material Fourier links remain active. Voltage, total-current,
regulated-power and distributed-current-density drives refresh the original
conductivity at each receiver's temperature and declared liquid fraction.
The returned candidate contains both accepted owners and their common
electrical-work/thermal-energy accounting.

Original cell routes must match the complete swept volume and constant-density
caloric chart. `step_axisymmetric_contact_candidate` adds finite electrical
interfaces through the existing `ConductionContactThermal` bindings. Each
original contact refreshes its source from both actual thermal/phase states;
fixed pressure/history, material order and selected sources remain pinned.
Its exact swept-surface loss enters the original adjacent receivers under the
authored two-sided heat fraction, separately from bulk cell heat and exactly
once in the common energy account. The volume-only entry continues to refuse
contact-bearing electrical owners.

No radial thermal mesh, expansion law, inferred heat split or second heat
integrator is introduced. Original Fourier links remain independently
authored; electrical contact resistance does not imply thermal resistance.
See [axisymmetric electrothermal semantics](docs/axisymmetric-electrothermal.md)
for source, geometry, budget and continuation details. Focused manufactured
regressions cover sector scaling, unequal radial volumes, Fourier transfer,
temperature/latent-phase feedback, swept contact work, asymmetric heat
partitions, source continuation, compliance, refusal and atomic retry.

## Native adhesive rolling with finite-cell snap heat

`friction_heat::adhesive::snap::rolling::native` stages the original unforced
native rolling event phase and finite-cell release together. Native and heat
owners must retain the same sealed lifecycle before preparation. The source
comes only from the heat binding; its original ordered normal/angular replay
must reproduce the native candidate's complete endpoint, normal report and
rolling report exactly. No constitutive or mechanical solver is duplicated.

The combined compensated audit uses actual native kinetic and mechanical
rolling storage changes, original plastic loss, actual caloric storage and
named outside work, with normal and rolling snap releases kept separate until
final rounding. Existing native and thermal gates remain active. Plastic loss
remains explicitly unconverted; it is never deposited a second time as snap
heat. Bodies, native lifecycle, heat lifecycle and caloric cells publish only
after every check and a final cancellation poll after the caller checkpoint.
The G1/G3/G4 cases in `tests/adhesive_snap_heat/rolling/native.rs` specify real
yielding recoil, loaded separation and fresh recontact, independent physical
energy closure, and atomic caloric, history, stale-replay and final-callback
cancellation refusal with deterministic retry. This remains a held-pose,
prescribed-normal, frozen-source composition: mechanical duration does not
advance the thermal clock or infer intervening transport or temperature feedback.


## Regional evaporation and original native gas-link pressure

`liquid_gas::evaporation::multiple::regions::pressure` adds regional pressure
loading after the original phase and optional reciprocal sensible-heat solves.
Every applied pressure comes from that finite parcel's actual post-exchange
EOS state. All regional sources still enter one native mass/enthalpy transaction
and exactly one original native flow step.

`RegionalGasPressure` supplies an explicit integer gas-coordinate membership
map and the original maximum lattice-density excursion. Every native D3Q19
gas-facing link must be covered, including diagonals and interfaces without
phase patches. Each authored phase face must agree with that map. Missing or
contradictory membership refuses the private candidate; there is no ambient
pressure fallback. The original fixed pressure reference, stationary solid
geometry, gas-volume, source, native conservation and cancellation gates remain
active. Receipts retain actual initial/applied pressures and native link counts,
which are streaming characteristics rather than interfacial areas.

The unpressurized APIs retain their original flow route. Regional pressure is
mechanical boundary loading with prescribed membership and fixed external gas
volumes. It does not infer gas momentum, moving gas-volume work, pressure heat,
topology remapping or closure of total mechanical plus thermal energy. Phase
and sensible heat enter the original caloric account once.

Focused regressions cover local pressure response, post-phase/post-heat EOS
loading, unpatched interfaces, complete native link coverage and atomic refusal.
Runtime execution remains unavailable under the session's dependency/resource
block and subsequent loss of its execution workspace.


## Gas transport followed by regional evaporation pressure

`liquid_gas::evaporation::transport::regions` now supplies explicit pressure
variants for both phase-only and sensible-heat/phase coupling. Original gas
paths first sample one network snapshot. The transported parcels then feed
their regional original phase/heat owners, and the resulting local EOS
pressures load the original native links during exactly one liquid flow step.
The existing final gas-map driving-force and aggregate absolute conservation
checks remain the publication gates; a late pressure or network refusal
discards gas transport and liquid changes together.

These routes reuse the same transport/coupling helper and pressure map as the
separate native regional pressure API. Pressure receipts begin at the
post-transport, pre-phase/heat gas state and end at the actually applied state.
Transported partial enthalpy, sensible heat and carried-vapor energy are each
accounted once. Fixed external gas volumes, prescribed regional membership,
first-order splitting and the absence of gas-volume work remain explicit.
Focused regressions specify original nozzle transport, independent final EOS
and native velocity, complete-system species/mass/caloric balances, empty-path
equivalence, late coverage refusal, retry and cancellation. These newly added
regressions have been reviewed statically but not executed.


## Finite material-wall heat with bidirectional molecular supply and ambient gas

`rocket_blowdown::effusion::thermal::exchange` composes the original implicit finite material-wall/gas heat transaction with the original bidirectional molecular particle transaction. Vacuum outlets, comoving reservoirs and fixed-inertial ambient Maxwellians sample one heated gas state, pool gross original-species withdrawals before incoming credit, and update the actual vessel once. Active apertures admit the unheated original gas first; zero-area sources remain inactive. Incoming molecular energy and tangential momentum remain in their particle streams and are not added as a second heat source.

Wall storage remains part of the existing dry vessel, with fixed geometry and constant-density calorics. Both stages share one physical interval and retain first-order splitting; source support, Biot, finite-inventory, drift, frame and mechanical gates remain authoritative. The shared source-relative energy account sums absolute thermal/particle defects and native inventory/gravity reserves under one fixed budget. Wall, gas and vehicle publish together; late source, withdrawal, callback, balance or cancellation refusal leaves inputs unchanged. Prescribed source inventories remain external and undepleted. Exterior wall-pressure forces, rotation, spatial thermal equilibration and empirical heat-transfer accuracy are not inferred.

Four focused regression groups cover independently calculated comoving supply/discharge and finite heat, mixed ambient/comoving/vacuum ports and inertial transforms, exact insulated-wall parity and inventory reserves, and atomic source/withdrawal/callback/cancellation refusals. Production formatting/source checks and an independent monatomic balance calculation completed before the executor outage; full owner tests were not compiled or run. Production sources were recovered byte-for-byte; the new test source was recovered from its complete pre-format patch.

## Causal material memory and caloric heat acceptance

`causal_heat_assembly_candidate` connects the original solid causal history
to the existing harmonic solid-to-thermal candidate solver. It deposits the
actual interval Maxwell/Kelvin heat, including energy released from retained
memory with zero AC forcing. The shared caloric solver retains its original
source/inventory admission, cell mapping, thermal residual and energy account.
The observation is evaluated once before its nonlinear thermal trials.

The returned candidate owns both unpublished endpoints. Only successful
material evolution, thermal solution and final cancellation admission return
the new material history and thermal simulation together; late failures leave
both accepted inputs unchanged. Interval work/heat/storage remain distinct
from cumulative material ledgers. The existing thermal owner has an f64 clock,
so its start and target must match the retained DD material interval exactly;
the standalone solid owner retains its wider DD-clock capability.

This is a fixed-law, prescribed-displacement material/heat handoff using the
existing caloric time discretization. Temperature-dependent material feedback,
thermal eigenstrain evolution and a mechanical equilibrium solve are not
inferred. Three G1/G3/G4 regression groups specify retained-memory relaxation
with zero AC forcing, consecutive coherent forced intervals against independent
real-time work/heat integrals, physical caloric endpoints, and late-failure/
clock/cancellation atomicity. They were authored and independently reviewed
but not native-compiled or executed after the execution environment failed.

## Sourced current-temperature absorption in native coherent heat

`SchottThermoOpticalSource::try_with_current_absorption` composes its original
T0 reference phase with a `SpectralDielectricAbsorptionStatePoint` carrying the
actual cell's original wavelength/alpha claims and same-state permeability.
The shared Maxwell response applies the catalog thermal phase correction from
the immutable reference and evaluates the supplied current loss once. The
source records separate original phase and absorption gap admissions; their
supported wavelength bands must intersect.

Both original monochromatic and broadband heat owners refresh that complete
actual source at each accepted cell temperature, retain its selected claims,
and reuse the existing complex field deposition and implicit thermal owner.
Current loss is never frozen by an ordinary state-only refresh, and reference
phase claims are never queried at a hot temperature. Applied and endpoint
optical receipts retain the full current absorption source. A late alpha-domain
or endpoint response refusal discards the candidate heat step atomically.

The lossless and explicitly frozen-reference-alpha constructors retain their
original behavior. This increment does not infer an absorption-temperature fit,
material expansion, stress, pressure feedback or additional heat. Authored
regressions specify original-source identity, hot-state alpha/mu refresh,
independent complex Airy monochromatic and broadband caloric response,
accepted-endpoint continuation and atomic current-source-bound refusal.
Native regression execution and final formatter verification were unavailable
after the session's executor outage; independent formula and source reviews
do not substitute for that pending execution.


### Physical forced-cycle objective adjoints

`PreparedPhsForcedSensitivity::mean_objective` exposes the shared two-grid value
and full-covector admission on unchanged physical state trajectories.
`mean_stored_energy` evaluates the original nonlinear Hamiltonian and complete
gradient, including explicit memory. `prepare_adjoint` and
`PreparedPhsForcedAdjoint` reuse the original coupled factor and exact source
direction owner for physical Fourier-port and storage/material gradients;
there is no source fitting or per-parameter waveform solve. Returned derivatives
are the direct objective partial minus the original residual-dual contraction.
Port count, physical clock and retained bandwidth remain mandatory; local-island
paths refuse prescribed clock changes. The author supplies direct objective
material/clock partials separately, including `partial_p H` for a storage
objective. Cancellation and refusal publish no result and preserve reuse.
Finite-order/sampling evidence does not prove harmonic-tail convergence or
stability.


### Physical autonomous-cycle objective adjoints

`PreparedPhsAutonomousSensitivity::mean_objective` and `mean_stored_energy`
expose the complete original cycle, with independently admitted value and
coefficient/frequency covector. Stored energy evaluates the source Hamiltonian
and every memory-coordinate partial and has zero explicit clock partial.
`PreparedPhsAutonomousAdjoint` retains the original phase border and solved
frequency. Its `gradient_parameter` shares the exact J/R/G/storage/feedback
direction owner with `solve_parameter`, including the direct conjugate-flow
chain once. The caller supplies the direct objective parameter partial at
fixed coefficients and frequency; the frequency partial belongs in adjoint
preparation, including changing-period effects for per-cycle energy. There is
no waveform-direction solve, source perturbation, refit or state reduction per
objective parameter query. Existing physical source refusals, cancellation,
budgets and two-grid direction admission remain. This is a local finite-harmonic
fixed-phase derivative, not a proof of stability or harmonic-tail convergence.


### Shared intrinsic creep in piano hammer and damper contacts

The grand-piano `felt_relaxation` adapter delegates Prony conversion, prepared
Kelvin histories and positive heat to
`fs_material::fiber::felt_rate::CreepSpectrum`. It supplies only the
work-conjugate area/thickness map; mechanical Memory stays in metres.
Hammer and finite-gap damper forces retain the original simultaneous contact
solve and WoolFelt discrete-gradient quadrature. Reference acceptance
propagates `try_advance` refusal through existing frame transactions; the
live path include retains its admitted-history compatibility entry point.

Zero-arm compliance, memory and heat stay exact zero. Resolvable weak arms
formerly discarded by the example cutoff and representable long recovery
tails are now retained. Rate-enabled audio is a declared numerical/physical
change; unchanged rendered audio and measured calibration are not claimed.

Seven focused adapter groups pass natively against the actual shared source:
original complex compliance, analytic SI creep/recovery, held-force work/heat,
area partition, weak loss/tails, zero arms and refusal/retry. The combined
source run passes all 24 shared-kernel, original felt/Prony and adapter groups.
Full reference rendering and engine-piano integration remain unrun. Generic
DrivenBody assembly and the remaining felt-adapter unification remain open.

## Temperature-dependent Lorentz refraction in coherent heating

The existing `layered_wave_heating::interference` and
`spatial_interference::feedback` owners admit original Lorentz polarization,
with optional independently sourced Ohmic leakage. Each accepted temperature
refreshes the original complete material selection, preserving the named
high-frequency permittivity, pole increments, angular frequencies, damping,
permeability and additional requirements. A combined Lorentz/Ohmic law resolves
that bundle once and shares the same refreshed snapshot between both laws.

Every frequency rebuilds the original common field stack from those local
cell states. Two-face illumination, polarization, transverse refraction,
continuous-footprint interference and real/proper source factors retain the
native electromagnetic owners. The original inertial reservoirs provide
positive damping heat `gamma Pdot²/(Delta_epsilon omega0²)`; their reversible
storage and independent Ohmic work remain separate. Native spectral products
retain frequency cross terms rather than substituting per-line mean losses.

The applied and endpoint source bundles travel with the original native
thermal transaction. The existing source-domain, frequency-response,
temperature-excursion, work/storage and cancellation gates remain authoritative;
a late refusal publishes neither heat history nor phase-clock advancement.
This remains an explicit quasi-steady split with settled periodic reservoirs
held over each interval. It does not transport transient polarization or
parameter work across a temperature refresh, infer causal oscillators from
optical tables, move geometry or resolve lateral heat transport.

Three focused G1/G3/G4 regressions specify independent oscillator-ODE heat
under mixed coherent illumination, nonuniform native-cell refraction and
caloric continuation, and atomic late Lorentz-domain refusal. Native Rust
execution is pending because the recovered environment has no Rust toolchain;
an independent scalar field/heat calculation and source checks are bounded
implementation evidence, not a runtime or measured-material validation claim.

## Finite-cell conversion of rolling plastic loss and junction release

`friction_heat::adhesive::snap::rolling::dissipation::JkrRollingDissipationHeat` preserves the original loaded rolling lifecycle and caloric bindings while adding an explicit partition and named outside destination for new irreversible rolling loss. Normal snap, remaining spring energy released at separation, and continuing plastic rolling loss retain distinct work ledgers. Earlier losses are not deposited again. Signed normal/area parameter work, elastic unloading and recoverable spring storage do not become heat.

A shared source-replay kernel serves both the unchanged snap-only route and the complete-loss route. It regenerates every original normal/angular segment, retains the sum of original source losses in double-double arithmetic, and applies all admitted heat to the original finite caloric cells in one proposal. Capacity/domain refusal, stale history, cancellation or energy imbalance leave every accepted owner unchanged. Explicit partitions model conversion/destination; they do not claim a solved interfacial conduction split.

The `dissipation::native` composition prepares actual body recoil, replays exactly the accepted source path, and checks the complete sealed successor, original normal report, rolling step and new plastic-loss increment before publication. Its combined body/contact/caloric balance counts converted plastic loss once, as caloric or named outside energy. Bodies, native lifecycle, conversion history and thermal state publish together after the final cancellation checks. The source remains frozen, normal motion prescribed, and heat a same-time impulse with no intervening transport or temperature feedback.

Six added source/native integration groups cover elastic/no-heat behavior, separate partitions across yielding/separation/recontact, retained prior loss, sub-ULP loss beside snap release, actual body/contact/cell energy, caloric refusal after a real recoil proposal, final-callback cancellation and retry/stale replay. Pinned formatting and recursive syntax/module resolution passed; these cross-crate groups were not type-checked or executed under the documented unavailable remote runners and shared build-resource limits.


## Finite and sealed acoustic source boundaries

`upstream_airway::AirwayChain` admits the same passive far-source law in impedance
or admittance coordinates. The admittance route uses
`U_source = Y_L (P_s - p_R)` in `exp(-i omega t)`; zero admittance is sealed,
including a compliant source's DC limit and finite-frequency antiresonance.
Scaled boundary-wave coefficients avoid both reciprocal overflow and growing
transfer matrices. Reflection, source-pressure transfer and input impedance
retain the original tube radius/count, Bessel viscothermal gas, yielding walls
and native two-port composition. An impedance pole may refuse while its
reflected wave remains finite.

`UpstreamAirway::realize_with_source_admittance` reuses the existing bounded
reflection and source-transfer FIR ladders, one shared material/source evaluation
per distinct frequency, their original passivity/convergence admission, and the
same implicit realtime valve junction. Exact zero DC conductance sets
`R(0)=1, G(0)=0`; no large finite resistance substitutes for a sealed boundary.
Receipt resistance/conductance fields are optional finite values: the authored
coordinate is exact, its reciprocal is present only when representable, and
`None` never means zero. Existing finite-impedance construction retains its
physical response.

This is a fixed-geometry, zero-state linear wave realization. It supplies no
independently initialized reservoir pressure, changing lung volume, local
stored-energy ledger, or causal/passivity certificate for an arbitrary supplied
frequency law. Four focused G0/G1/G3/G4/G5 groups cover loaded ABCD and closed
limits, compliant reflection phase and FIR DC bins, reciprocal-coordinate
realtime agreement, extreme scales, invalid inputs and callback cancellation.

## Driven spatial damping with finite caloric feedback

`rigid_damping_heat::driven` connects the original simultaneous driven-pair
motion owner to the existing native thermal event/transport transaction.
The input is exactly two original bodies and one [0, 1] connector; arbitrary
driven graphs are not inferred. Original attachment arms, full inertias,
passive spatial cross-coupling and material/world-held effort retain their
native interpretation. Every actual thermal interval performs the original
two endpoint source preparations and native finite motion.

The shared report now exposes signed `actuator_work_j` separately from held
external applied work, gravitational storage and positive damping heat.
Only the independently integrated native damping loss is routed to finite
caloric cells. Actuator recovery remains negative supplied work; it cannot
become heat. The same original caloric inversions, heat fractions, unallocated
outflow, source/transport hooks, event clipping and whole-window energy gates
remain active. A small private mechanical adapter shares this existing thermal
transaction; it does not introduce another time integrator or thermal store.

The temperature entry reuses `TemperatureDampingLaw` and
`CalibratedKelvinDamping` with the original reference form and an explicitly
selected source CellId. Accepted start temperatures select each interval's
viscous form; the original caloric endpoint must also admit its refreshed law
before publication. Held effort itself does not acquire a temperature model.
Source, observer, caloric, budget and cancellation refusals discard both physical
candidates. Moving geometry and lagged-temperature feedback still require
independent temporal refinement; feedback remains first order.

Six focused regression groups specify analytic terminal slip with positive
actuator supply and negative recovery, zero-law heat, native event/source clocks,
an independent coupled mechanical/caloric ODE ladder, endpoint-source and late
cancellation atomicity, and topology/coordinate admission. Pinned rustfmt passed.
An independent scalar calculation gives successive feedback error ratios
2.00107, 2.00103 and 2.00064, with balance residuals below 1.95e-15 J.
Native integration execution remains outstanding under the existing dependency
build resource constraint; these calculations are not full crate validation.


## Finite diffuse-channel recoil with material-wall feedback

`rocket_blowdown::effusion::channel` adds
`GasRocketState::advance_effusion_channels` for a translating finite gas vessel
with independently oriented molecular channel arrays and finite original
material walls. The existing material vacuum-channel transaction owns all
gross reflected and escaping species, one gas reconstruction, and
accommodation heat pooled by physical wall identity. Each qualified
`DiffuseChannelMomentumInterface` enriches that same accepted outlet receipt
with geometry-resolved axial momentum. The vessel encloses both source gas
and channel walls, so only escaping momentum enters the native forced exhaust
trajectory. Local source-gas and wall impulses remain separate receipts and
are not added again as thrust. Complete escaped molecular energy includes
accommodation before its inertial-frame transform.

Distinct wall inventories must fit within the vehicle's existing dry mass,
with shared wall IDs counted once. Gas mass remains the authoritative
propellant inventory. Native source/rate/remapping corrections and gravity
mass-history reserves retain the original physical budgets. The complete
energy gate also bounds the sum of absolute material and mechanical defects.
Gas, vehicle and all wall endpoints publish atomically. Geometry, axes, loads
and initial molecular rates are held for this first-order quasistatic interval;
repeated calls use actual evolving gas and wall states. Rotation, transit delay,
spatial rarefied-flow resolution and additional continuum pressure thrust are
outside this point-vessel model.

Five G1/G3/G4 test groups cover independent monatomic hot/cold gas and nonlinear
wall storage balances; local traction versus escaping recoil; shared-wall
opposed mixture outlets and canonical replay; loaded Galilean momentum/energy
transport; closed, elastic and zero-length limits; and dry-mass, frame, gross
withdrawal and late-cancellation refusal. Pinned formatting and scoped
whitespace checks pass. Focused source/API/physics review found no concrete
blocker. These native integration tests have not been compiled or executed
in the constrained workspace; geometry-only native tests do not prove this
cross-crate integration.

## Layer-resolved tensor damping on the original MITC shell

`TensorShellRelaxation::from_layered_mitc_shell` accepts one explicitly ordered,
contiguous stack of original `TensorShellMaterial` sources per physical facet.
Each layer retains its own orientation, signed faces, shear correction, material
band, complete Kelvin matrix and every Maxwell branch. The existing five-strain
shell constraint eliminates only thickness-normal traction on the whole source.
No independently condensed branch, proportional coating fit or numerical
drilling loss replaces that constitutive response.

The complete integrated eight-coordinate section must match the original
A/B/D/shear operator, including absent membrane/shear and bending/shear cross
blocks. The signed density moments must match every original corner's reference
inertia. A gap, overlap, changed reference surface or incompatible source refuses.
The existing support map, elastic bubble condensation, mass choice and modal
pencil remain authoritative. Geometry and in-plane strain samples prepare once
per facet; each layer contributes its original seven area and two thickness
points, counted against the shared point, history and work limits.

One facet-major material catalog identifies every point's actual source layer.
The original equal-time coordinate owner shares driven histories while retaining
reconstruction of every original point. All Kelvin channels and Maxwell memory
attach through the existing coupled runtime and its independent physical heat
accounts. The original single-layer entry uses this same implementation and
keeps its established homogeneous admission and arithmetic.

This is the native equivalent-single-layer displacement model with authored
work-conjugate shear corrections. It does not resolve interlaminar stress, slip,
delamination or warping. Material state is frozen; point strain checks do not
certify spatial or between-step maxima, and no experimental damping fit is
supplied.

Two isolated native tests execute the actual section-integration function with
real fs-math, checking independent original-strain quadrature for 129 fields and
identical-layer partition invariance, including signed inertia. Four additional
integration groups specify original 3D complex-work agreement, native coupled
heat closure, homogeneous/partition parity and source/geometry/budget refusals.
Pinned formatting passed; the full fs-couple integration groups remain unrun.

## Held expansion through causal material-to-thermal acceptance

The causal harmonic heat candidate accepts the solid owner's original frozen
thermal snapshot. It preserves that snapshot and the existing eigenstrain-work
ledger while depositing only interval Maxwell/Kelvin heat in the existing
caloric solver. Both endpoints still publish together after source, residual
and cancellation admission. Caloric temperature changes do not implicitly
change the held expansion or the constitutive law.

Two additional tests use actual thermal-strain preloading and the independent
mechanical-memory oracle. They specify repeated static heat release, coherent
work/heat, snapshot and nonzero prior eigenstrain-work preservation, and late
thermal refusal/retry/cancellation. Formatting and bounded review completed;
full original Cx-backed coupling execution remains unavailable. The expansion and
constitutive coefficients remain fixed.

## Applied-body work and rolling heat across contact events

`JkrRollingDissipationHeat::propose_native_with_wrenches` and `apply_native_with_wrenches` combine actual held COM forces/body-frame torques with the original attachment/separation/recontact lifecycle and all three finite-cell loss channels. Existing unforced methods delegate with zero loads. The combined balance uses the native owner's actual kinetic increment and external work in both attached and detached branches. Converted rolling loss is already represented by caloric or named outside energy and is never counted again as outgoing dissipation. Detached applied work changes kinetic energy without inventing contact heat; subsequent attachment starts an unstrained junction.

The same source replay, original caloric domains, independent work gates and all-owner cancellation/publication semantics remain active. Two further integration groups exercise sustained applied torque from rest through spring loading and plastic heat over eight steps, plus detached force/torque acceleration with zero heating followed by fresh recontact. Actual body/contact/cell energy is checked independently. Pinned formatting and recursive syntax/module resolution passed; the new integration groups remain uncompiled and unexecuted under the documented unavailable remote runners and shared resource limits. This is still prescribed-normal, held-geometry, frozen-source impulse coupling: mechanical duration does not advance the thermal clock, and conduction or temperature feedback must be supplied by their existing owners.


## Joint native MITC body preload

`fs_couple::body_assembly` solves thermal eigenstrain and explicit dead nodal loads jointly across the original `PreparedThermalShell` owners and native rigid/compliant junctions. Each reference must be explicitly unloaded. Homogeneous global support coordinates and the native flat-reference, small-displacement membrane initial-stress envelope apply. Prescribed slave supports that the native mortar elimination cannot compose are refused.

The original sparse elastic factor/refinement solves the tied equilibrium once. Admission checks the returned coordinates against both the original assembled force/moment residual and the original component forces plus glue forces, using a caller-declared moment length and absolute force tolerance. Recovery uses each same thermal/material owner at that shared displacement. Existing junction rows then assemble the signed native component tangents and add original glue stiffness exactly once; reference mass, supports, damping and maps are unchanged.

The accepted body and its native `fs-modal` windows remain bound by immutable borrows. Physical mode expansion and force projection use the same work-conjugate junction map. Negative eigenvalues remain signed, and a positive requested window alone does not establish global stability. Static displacements are reference offsets, with no fictitious time, heat or work credit. This slice does not provide general nonlinear curved-shell equilibrium, follower loads, changing contact, beam/solid prestress composition or post-buckling.

## Joint wall and opening radiation of modal cavities

`vibroacoustic::CavityOpeningRadiation` adds actual outward volume-flow
unknowns to a rigid-wall cavity pressure basis. One original sampled exterior
matrix acts on [full structural modal velocities, opening volume flows],
retaining wall/mouth and mouth/mouth reaction and interference. Its real cavity
map and transpose transfer outward modal flow and physical inner pressure.
Opening flows stay in the simultaneous solve; no cavity or neck resonance is
divided out. The complete external matrix uses the existing full Hermitian
passivity admission and exact source frequency grid.

The native Bessel duct owner supplies each internal cylinder's complex series
impedance at every retained frequency. No resonance-frozen resistance or
additional exterior end mass is inserted. A resolved cavity owns its entrance
inertia; an optional fs-phs inner-end correction is restricted to a single
uniform pressure mode, preventing double counting as higher pressure modes
resolve that field. The declared bound on |k| times the larger of radius and
effective neck length cannot exceed 0.25. These are short, straight, quiescent
necks with one volume flow, not distributed ducts, resolved entry-edge loss,
thermal compression storage, transverse modes or nonlinear jets.

The cavity's uniform zero mode is required. An already neck-loaded Helmholtz
mode or an independent structural exterior attachment cannot be combined with
this joint loading. Original structural viscosity/relaxation, cavity molecular
compressibility and disjoint local wall boundaries remain in the same equations.
Neck and cavity gas must agree; the external producer retains its own medium.

Independent accounts report signed cavity-to-opening work, positive internal
neck heat and complete coherent exterior work:
`input = structural + cavity + opening_necks + radiated` and
`interface + acoustic_input = cavity + cavity_openings`.
Opening transfer is not a second total-system debit. Returned original port
amplitudes can excite the same native BEM source solutions for field observation.
Truncated solves retain every opening and zero omitted structural coordinates
at their original exterior indices. Convergence additionally observes the full
complex opening-flow vector.

Four native numerical tests passed with actual fs-math and the original fs-la
complex LU, using the exact new operator, physical response and Hermitian work
bodies: independent Helmholtz transfer; full mutual work and pressure-basis
rescaling; coherent wall/opening cancellation; and retained-prefix port
indices/flow convergence. The isolated carrier narrows only the error type.
Three public integration tests additionally specify native Bessel geometry,
a genuine joint Rayleigh bake with hereditary structure, source-frequency
identity, double-counting refusals and compact/source/cancellation boundaries.
Their full Cargo execution remains unrun; pinned formatting and source/API
review do not substitute for that execution. This is a sampled harmonic modal
composition, without a causal time realization or continuum error certificate.

## Coupled resolved laminar momentum and native finite enthalpy

The `laminar-enthalpy` feature exposes `LaminarEnthalpy`: one original
fixed-mesh incompressible BDM continuum and native finite-cell enthalpy timeline.
Every implicit thermal trial resolves the original material BDF1 convection and
rheology equations from the SAME accepted velocity and held physical boundary
data. Each original constitutive source receives the actual trial caloric
temperature and native objective shear. Its thermal and momentum transport views
must agree and retain the original constant density; the exact source object and
accepted endpoint state are returned. No transport heat-capacity sample is
converted into a caloric law. Original identity-bearing enthalpy materials,
triangle-prism geometry and finite masses remain authoritative.

The actual sign-split affine BDM normal trace supplies two directed mass streams
per face when needed, including zero-net counterflow. Rates are original volume
flux times the declared extrusion depth and density. Native upwind enthalpy
transport uses each local donor's original enthalpy; every actually incoming
boundary portion needs an explicit external material and specific enthalpy.
One interior transfer is shared by its two inventories. Original fixed-mass
continuity and caloric-datum admission remain active, with no mass-flux repair.
Temperature-dependent physical surface/link conductances remain supplied by the
original heat callback at that same coupled endpoint.

Original physical bulk strain heating is integrated over each triangle prism.
Every Navier-slip edge needs an explicit original wall receiver and fluid heat
fraction; pressure-invariant native slip loss is partitioned once between the
two receivers. SIP, upwind, open-boundary stabilization and BDF numerical losses
are not deposited as material heat. Original authored heat and external thermal
exchange are measured independently before this composition adds dissipation
and advected enthalpy. Acceptance requires native momentum convergence,
local/global original cell-energy residuals, and the independently evaluated
identity `delta H = bulk heat + slip heat + net stream enthalpy + other heat`.
No residual is reclassified as heat.

The existing fallible finite-difference and Newton--Krylov owners perform bounded
thermal iteration; callers choose a directional difference scale in joules and
an inexact linear forcing floor appropriate to that numerical action (for
example `forcing_minimum = 1e-4`). This does not relax final physical thermal
tolerances. Source/domain, stream, work-budget, nonlinear, thermal-balance or
cancellation refusal returns neither replacement history. Both clocks and all
original inventories remain unchanged. Callbacks must be pure across nonlinear
probes; cancellation is checked at composition and source boundaries. Native
sparse-loop interruption is not added.

This is first-order endpoint coupling at fixed density and geometry, with native
mixed velocity, slip and explicit open-pressure conditions. It does not infer
thermal expansion, whole-system mechanical/thermal passivity, within-step flow
reversal times or a temporal-error certificate. Five authored native regression
groups specify bidirectional face transport; independently derived two-cell
throughflow heating under reversal and shifted enthalpy datums; original
thermal-Carreau viscosity, convection, strain heat and storage; physical slip
partitioning; and transactional source, density, exchange, inflow, budget and
cancellation failures. These full native FEM/thermal tests remain unexecuted in
this session. Rustfmt parsing, focused whitespace checks and bounded independent
source/API review passed; earlier isolated material/thermal checks do not
substitute for this pending integration execution.

## Coherent field heat over accepted thermal intervals

`SpatialInterferenceHeating::advance_to_with_field_heat` explicitly reconstructs
and integrates the original coherent reservoir fields for each native thermal
cell. Its per-cell `AngularSpectrumHeatWindowLimits` admit time projections,
original spatial integration and scratch before field reconstruction. The
original continuous footprint, positive depth samples, material cards, source
statistics and optical films remain authoritative. Source phase retains the
accepted double-double epoch and the actual interval selected by the native
thermal step. The EM interval query returns joules; no additional area, volume
or duration scaling is applied at the thermal boundary.

The field and cached-spectrum routes share cell order, the independent positive
sum of cell source-ripple bounds, finite nonnegative energy admission, and one
native heat debit. Original caloric laws, conduction, contacts, baths, authored
sources and nonlinear acceptance still belong to the thermal owner. Refusal or
cancellation publishes neither temperatures nor source-clock advancement.
Frozen optics, through-thickness thermal resolution and expected-heat statistics
retain the existing approximation boundaries; this query does not infer lateral
thermal resolution, reservoir switch-on dynamics or thermo-optical feedback.

Seven focused actual-source collector tests pass with warnings denied, including
fourth-order dark-fringe joules, native cell order, retained sub-ulp source
epochs, ripple refusal before field work, and cancellation. These exercise the
shared harmonic source and accepted-interval energy seam. Full native material,
field reconstruction and thermal-solver integration execution remains unrun.

## Adaptive finite-wall bidirectional and inertial ambient effusion

`rocket_blowdown::effusion::thermal::trajectory::exchange` adds
`GasRocketState::advance_effusion_exchange_with_wall_adaptive` and its
convection-source variant. Vacuum discharge, comoving supply and fixed-inertial
ambient capture share the existing finite-wall step-doubling controller.
Every full and half trial calls the original thermal/particle exchange owner
at the actual gas composition, gas temperature, finite-wall state and vessel
velocity. An inertial ambient velocity is never silently made comoving.

The existing vacuum API/report shape is preserved. One private generic driver
owns time, wall/species error control, accepted history and the fixed-window
audit; native report adapters supply the different molecular accounts. Incoming,
outgoing and authoritative signed-net momentum and energy participate
separately in temporal error control, so nearly cancelling streams do not hide
gross exchange errors. A species entering from zero joins the species error
axis without a fictitious seed inventory.

Molecular reference energies use original caloric integrals and accepted
crossing-energy corrections, including drift, before multiplication by actual
species mass. Incoming energy is not reconstructed from net plus outgoing
energy or from subtraction of large formation datums. Each species has one
fixed reference for the window: the original vessel temperature when present,
otherwise the first active external source in canonical port-ID order.
Original material models retain only positive-composition species; inactive
catalog entries do not acquire a baseline from an absent gas. Native boundary
admission still requires the actual heated pre-exchange gas temperature.

Accepted native reports retain complete inertial molecular energy, separate
external heat/radiation and force/work, and signed source-coordinate net totals.
The original fixed-window species, mass, momentum, energy, rounding/reserve and
prescribed-heat-prorating budgets remain binding and do not grow with step count.
Only typed gross-withdrawal errors permit physical retry; source errors,
cancellation, balance failures and exhausted work refuse atomically.
External reservoirs remain prescribed/undepleted and receive the opposite
original port accounts. Export tolerances are defined in the declared inertial
frame; a common observer boost need not produce identical adaptive step choices.
The estimates concern local temporal splitting, not certified global or
constitutive error.

Four G0/G1/G3/G4 test groups cover entering-species refinement and exact native
history replay; an independent continuous monatomic gas/wall/supply/discharge/
recoil balance; independent incoming/outgoing/net reference energies under
formation-datum shifts and common boosts; and fixed-window budget, work,
minimum-step, unsupported-source, closed-source and cancellation behavior.
Pinned formatter parsing/checks and scoped whitespace checks pass. Source/API
and physical-account review found no concrete blocker. Native compilation
and execution of these cross-crate integration tests were not performed in the
constrained workspace.

## Coherent interval heat through temperature-dependent optical feedback

`QuasiSteadySpatialInterferenceHeating::advance_to_with_field_heat` refreshes
every original native cell's complete Debye, Lorentz or Ohmic source bundle at
its accepted temperature, rebuilds the original common frequency stacks and
passes the original coherent fields to the bounded volume-heat interval query.
Per-cell `AngularSpectrumHeatWindowLimits` are explicit and independent of the
existing preparation limit. Original transverse carrier, finite footprint,
real/proper source statistics, native cell order and DD epoch are preserved.

Both advance routes share the same applied/endpoint transaction, source-ripple
gate, one native heat debit, per-cell temperature-excursion gate and complete
endpoint material/frequency admission. Endpoint snapshots add no second heat.
Late refusal or cancellation publishes neither the thermal proposal nor an
advanced source clock. This remains a quasi-steady split with settled optical
reservoirs held within each interval; it adds no transient polarization state,
parameter-work transfer, moving geometry or lateral thermal resolution.

Existing native fixtures now include three-tone mixed banks and three field
route cases for local Debye/Lorentz refresh, precise epoch continuation,
projection-budget refusal, caloric balance and late domain/cancellation refusal.
Formatting, source/API checks and a native DD clock probe passed. The complete
material/EM/thermal fixture graph was not compiled or executed in this session.

## Rate-dependent felt in original exponential impact mechanics

`FeltPad::from_material` binds the admitted intrinsic `RateFelt` envelope
and up to eight original Kelvin arms to explicit contact area and thickness.
The map `k=E*A/h`, `c=eta*A/h` preserves work, stored energy and dissipation.
It supplies no material defaults, fitted loss curve or separate held-stress
update. Existing raw pad descriptions and the zero-arm force/history arithmetic
remain available through the original owner.

`ExponentialImpactSystem::new_with_felt` retains all Kelvin storage coordinates
in the existing complete quadratic block and their original resistance in the
same pHS equation. The original nonlinear WoolFelt potential and contact force
remain in the implicit discrete-gradient remainder. Its compaction histories
are frozen during that solve and commit through the existing joint acceptance
gate with mechanics, creep, positive loss and time. Every material coordinate
counts against the unchanged dense propagation capacity; none is truncated.

The felt remainder reuses fs-dcontact's original unilateral power primitive,
splitting only at the original conditioning maximum. Compensated normalization
retains small endpoint work that would otherwise contaminate a cubic remainder.
Geometry, material domains and source parameters remain authoritative. The
opt-in exact linear block removes midpoint decay/phase error from original
linear modes and Kelvin recovery; nonlinear contact and conditioning still
require temporal refinement. This route does not yet admit nonlinear contact
drag, shell potentials, transverse-only support images, undeclared or
nonquadratic strain maps, or thermal feedback.
It does not certify a rendered waveform or measured felt calibration.

Seven focused native source tests passed with pinned nightly-2026-07-06:
three original felt checks plus four groups covering independent piecewise
energy/work, small-increment cubic remainders, all-eight-arm physical geometry,
and zero-rate bit preservation/refusals. The isolated harness used actual
WoolFelt, RateFelt construction/accessors, FeltPad and shared power-remainder
source with real fs-math/fs-evidence/fs-qty. Four full pHS integration groups
specify analytic recovery/heat, original nonlinear ODE refinement, conditioning
and rollback, and complete-state capacity. Those full fs-couple integration
groups remain unrun under the existing dependency build resource constraint.

## Rolling plastic heat over physical thermal time

`friction_heat::adhesive::snap::rolling::dissipation::timed` adds positive-duration proposals and applications to the complete original rolling heat owner. Accepted contact and thermal clocks must match exactly. Original scheduled thermal events at the start precede the prescribed normal/area update and its caloric release impulse. The original monotone angular increment then supplies only its irreversible rolling loss as explicit interval-average power into the pinned finite receiving cells.

Native thermal conduction, reservoir cooling, ordinary heat sources and scheduled changes evolve over the same actual interval with their original caloric bounds and bounded substep work. The rolling source integrates power over accepted native durations and independently checks each receiving-cell total against the declared partition. Clock or power rounding must satisfy the explicit transfer gate; no endpoint heat correction or guessed dissipation repairs it. Source memory, distinct cumulative loss/release ledgers, caloric cells and accepted end time publish together only after every source, transport, energy and cancellation gate admits.

The combined transfer audit adds original compensated relative native audits and independently checked contact-transfer defects. Its absolute gate keeps normal partition, normal caloric, transport and per-cell rolling defects separate. Large rounded aggregate formation enthalpies are not subtracted to invent a different residual. The report exposes initial events, normal impulse, distributed rolling loss, actual injected heat, ordinary transport heat and the advanced thermal interval.

This is a frozen-source split with prescribed normal loading at interval start and average-power reconstruction of a monotone rolling step. Temporal refinement is required for within-step transients. Contact heat partition and thermal topology remain explicit caller-selected physics; temperature feedback and reciprocal contact conductance are separate operations.

Four focused integration cases cover independent two-cell conduction/cooling eigenmodes, snap followed by a scheduled cooling change, budget/final-callback/stale-clock rollback, and finite caloric refusal. Pinned `rustfmt --check` passed. These integration cases were not compiled or executed under this session's resource ceiling.

## Initialized acoustic source storage and characteristic reflection

`upstream_airway::StoredPressureSource` retains the original
`fs_vfit::impedance::SeriesImpedance`, including explicitly initialized
compliance pressure and inertive flow. Its existing `UpstreamPort` law is
`p = P_s + p_load(-U_out)`: one same-sample affine solve followed by a
transactional storage commit, with no delayed compliance response. Initial
energy remains initial storage, not source work or irreversible loss.

The native driven-network path subtracts half the ideal series source from
the arriving wave, solves the original R-L-C boundary, and adds it back to
the departing wave. It accepts the source and complete original network
only after one native preview and the combined finite-energy/cancellation
check. It requires one additional driven endpoint and identical physical
clock and adjacent characteristic reference. No wave history is cloned;
source work is O(1), and the original network's work and memory bounds apply.

Original quadratic storage, resistor heat, ideal-source work and signed
source-to-network work remain separate. The combined account cancels that
internal exchange exactly once. A charged compliance can release pressure
with zero ideal-source drive and later recover reflected wave energy without
calling that exchange heat. Cancellation or either owner's refusal retains
both accepted histories, the original network clock and pending valve port.

Seven focused native tests execute the actual source, R-L-C, propagation and
driven-network owners. They cover a three-rung second-order compliance decay,
known transit, source/duct energy conservation under repeated reflections,
initialized driven R-L-C work/heat, original wave/flow parity and transactional
refusal. Five existing original reactive-impedance regressions also pass.
The isolated build omits unused exact-parallel-cell and MIMO module wiring;
full workspace and native valve/airway integration tests remain unrun.

This is the fixed-coefficient linear acoustic model, with explicit bilinear
frequency warping. It is not a finite thermochemical inventory, changing
reservoir volume, inferred anatomical law, validated player population, or
storage certificate for an independently fitted airway FIR. Such a FIR must
not receive a second copy of its already included termination.


### Complete physical forced directions and harmonic power objectives

`PreparedPhsForcedSensitivity::solve_parameter` and
`PreparedPhsForcedAdjoint::gradient_parameter` accept one complete
`PhsForcedParameterDirection` for simultaneous original J/R/G, storage and
Fourier-drive changes, with a prescribed clock derivative in the existing
`ForcedSensitivityDirection`. The shared original physical parameter owner
forms `(J-R)h_p + (J_p-R_p)e + G_p*u + G*u_p`. The shared fs-orbit split
projection adds the actual pHS `s*I` clock column `i*n*Omega_p*X_n` while
retaining every nonlinear and memory state. Fourier input coefficients are
held in the fixed-phase chart; adding their fixed-time carrier derivative
would double-count the clock. Returned coefficients are fixed-phase
sensitivities; fixed-time observations additionally differentiate phase.

Preparation retains an immutable copy of the exact original drive, with its
payload and physical direction scratch reserved before numerical preparation.
Existing constructor input lifetimes remain unchanged. Direction port count,
base clock and retained bandwidth must match, including arbitrarily small
nonzero coefficients. Signed structure directions retain original skew and
symmetry admission. No perturbed model, fitted replacement storage or
per-parameter objective waveform solve is constructed. Source/arithmetic
refusal and cancellation publish no result and preserve reusable factors.
Tests compare simultaneous matrix/material/source/clock changes with analytic
mechanics and independently solved nonlinear-memory neighbors, including
changed-period objectives and malformed/cancelled direction recovery.

Both physical forced and autonomous sensitivity preparations expose
`mean_supplied_power` and `mean_resistance_power`, with `_analytic` variants
for explicit complete original storage Hessians. Values use the original
compensated pHS supplied/resistive-power owners independently; equality is
never imposed through energy balance. State covectors use the same complete
StorageTangent owner as HB, including nonlinear and memory states. Forced
supply uses the exact prepared drive snapshot. Autonomous supply requires the
same TOTAL original feedback Jacobian and includes
`Hxx*G*u + (du_total/dx)^T*y`, with its direct flow chain once.

The mean and complete coefficient covector pass shared two-grid admission;
autonomous fixed-state clock partials are zero. Per-cycle energy additionally
requires the period chain. Direct source/material/port objective partials
remain explicit caller inputs. Signed raw R remains signed; no heat destination,
passivity, stability or harmonic-tail certificate is inferred. Observation
scratch is O(n+nports) plus original Hessian metadata, with O(n*nports) extra
feedback-Jacobian scratch for autonomous supply. Source callback work/metadata
retain the generic mean-objective budget exclusion; evaluations are bounded
by the original two AFT grids. Cancellation, unavailable/partial Hessians,
partial feedback derivatives and arithmetic refusal publish no objective.
Six physical power regression groups verify complete analytic covectors,
independent nonlinear-memory neighbors, changed-frequency feedback cycles and
per-cycle energy, original nonsymmetric/cancelling resistance, and refusal/retry.

## Driven native rolling with physical-time heat transport

`friction_heat::adhesive::snap::rolling::dissipation::timed::native` composes the original driven finite-body rolling lifecycle with the positive-duration heat owner. Held COM forces and body-frame torques drive the actual native motion. Only the internally accepted native angular increment supplies the source replay and irreversible rolling loss; public edited reports or caller-authored heat cannot replace it. The entire retained source state, original normal segment, rolling step and plastic loss must match exactly before acceptance.

Original scheduled thermal events and the prescribed normal/area stage occur at the interval start; normal snap and remaining spring release use their original explicit release partition there. Plastic rolling heat uses its independent partition as average power during native conduction, cooling and scheduled events over that same physical interval. Detached driven work remains kinetic and supplies no rolling heat. Recontact starts a fresh spring without converting historical loss again.

The combined audit cancels internal plastic transfer once, using actual native kinetic change, original incremental spring storage and body work together with the heat owner's compensated relative audits. Independent mechanical and thermal absolute defects cannot cancel. Both native momenta, native contact history, complete heat history, finite caloric cells and thermal clock publish together after the final caller checkpoint and cancellation check. Late transport-budget, caloric, geometry, history or energy refusal leaves every original owner unchanged.

The physical scope remains held poses/frame, prescribed normal loading at interval start and frozen material properties. Interval-average heat and the original inertial midpoint phase require temporal refinement. Finite orientation integration, normal body motion and temperature-dependent source feedback are separate operations.

Three focused integration cases cover eight driven intervals with independent kinetic/contact/cell/outside energy and background power, separation and fresh recontact while physical time advances, and post-snap transport failure, both final cancellation forms, stale time/history and deterministic retry. Pinned format checks and recursive Rust syntax/module-path parsing passed. These new integration cases were not compiled or run: DSR and RCH are unavailable and local memory/disk headroom is exhausted.

## Accepted-temperature water/silica phase with original current absorption

`AnalyticThermoOpticalSource::{try_iapws_water,try_corning_7980_silica}` binds
one explicitly selected original analytic phase family to the canonical
independent current-state wavelength/alpha spectrum. The complete original
snapshot also owns scalar permeability and, for water, bulk mass density.
Accepted native cell temperatures refresh all retained selected claims once;
the existing water/silica Maxwell adapter then evaluates phase and loss at
that SAME state and retains its original finite support and gap admission.

The existing monochromatic and mutually incoherent broadband thermal owners
consume this sealed source through their original coherent film fields,
positive local absorbed work, caloric capacity, conduction and transaction
boundaries. Applied/endpoint snapshots preserve the full original source;
subsequent intervals start from the accepted endpoint. The original source
card, absolute-temperature coordinate, alpha domain and analytic phase domain
remain binding, including the CHARMS fit's 300 K ceiling. Cancellation or an
endpoint source/phase refusal returns no partial optical/thermal update.

The source owns no extra Beer-law or Ohmic heating term. Current source
coefficients are held over each original quasi-steady thermal interval; it
supplies no causal field reservoir, transient constitutive parameter work,
flow/EOS closure, phase classification, geometry change or expansion work.
Fixed thermal mass/geometry remain explicit caller-owned approximations when
an independently supplied density path changes. This does not qualify an
empirical loss spectrum or infer a material identity from its name.

Two authored native regression groups exercise BOTH original phase families:
nonuniform two-cell coherent heating against a constant-capacity backward
Euler oracle, next-interval refresh and late cancellation; and two-line
broadband complex Airy heat/enthalpy, original alpha refresh and transactional
phase/absorption domain exits. Native tests are defined but unexecuted here;
independent scalar Airy arithmetic confirms the finite heating and domain-exit
fixtures, and source/API review checks the existing native integration.

## Simultaneous driven rigid graphs with per-connector finite heat

`rigid_damping_heat::driven::network` connects the native simultaneous driven
spatial graph to the same finite caloric transaction as the passive network and
driven pair. Every original link has an explicit held material/world effort,
heat route and optional temperature coordinate. Cycles, parallel links, shared
bodies and isolated-body drift retain the native graph and original geometry.
Its complete dense preparation limits are explicit; it does not use passive
pair sweeps or invent a source subdivision.

The shared report retains each original connector's signed actuator work as well
as its independently integrated positive heat. Opposite actuator supply and
recovery may cancel in the total work while physical damping remains positive.
Neither recovered work nor source/drift defects become thermal energy. Passive
links report exact zero actuator work. Connector order, compensated interval
accumulation, shared receivers, partial heat allocation and unallocated outflow
are preserved through the original thermal event and observer timeline.

Both graph entries move the native endpoint body bank and work/heat vectors into
the existing transaction without cloning them. Temperature feedback reuses the
original reference forms and `TemperatureDampingLaw`; source or endpoint refusal,
observer refusal, budget exhaustion and cancellation discard the complete
mechanical/thermal candidate. Frozen graph source exactness does not establish
moving-geometry or thermal accuracy; temperature feedback remains first order.

Four focused integration groups specify an analytic three-body driven cycle
with nonlinear calorics and partial shared heat allocation, parallel actuator
supply/recovery, refinement to the independent coupled mechanical/caloric ODE,
and native events with late refusal/retry. An independent scalar calculation
gives feedback error ratios 1.99746 and 1.99883; work/heat balance stays below
3.5e-15 J in that calculation. Pinned rustfmt passes. Full fs-couple integration
execution remains outstanding under the existing dependency resource limit.

## Original wall mechanics behind dissipative solid coatings

`elastic_reflection::OscillatorBackedViscoelasticLayers` joins the original
incident fluid and bonded viscoelastic sources to the original normal
`WallPin` mass, spring and resistance. The native material owner retains
every layer's bulk/shear Maxwell spectra, frequency/strain limits, resolved
fields and positive constitutive heat. The adapter passes the same physical
support impedance at the same angular frequency, with explicit conjugation
from the pin's `exp(-i omega t)` convention to the stack's `exp(+i omega t)`.

Actual rear velocity determines backing displacement `x = -i v / omega` and
the original cycle-mean kinetic/spring storage. Both per-area storage and total
storage over the drive's positive area are returned. The native material
response keeps total intrinsic bulk heat, intrinsic shear heat and backing
resistance work separate; pin storage is never counted as heat or coating
storage. An empty coating recovers the original direct fluid/pin boundary.
Prepared drives allocate nothing and query no constitutive source again.

This owner is a stationary harmonic normal-longitudinal model, with no
fluid-shear input, oblique conversion, thermal time advance or transient
realization inferred from a frequency snapshot. Four focused G1/G3/G4
regression definitions cover the direct analytic pin response and phasor
conversion, original one/two-layer material sources against an independent
80-digit physical transfer/positive-work reference, area and reactive limits,
and original domain/budget/cancellation admission. Native Rust execution and
formatting status are reported separately.

## Distributed opening ducts in the joint vibroacoustic response

`vibroacoustic::DistributedCavityOpeningRadiation` extends the joint wall/mouth
loading owner to finite cylindrical ducts. The original sealed-cavity pressure
basis still starts with its uniform zero mode; cavity volume/modal norms exclude
the physical duct volume. `with_distributed_cavity_openings` is exclusive with
compact opening or structural-only exterior loading. One complete original
exterior impedance retains every wall/mouth and mouth/mouth interference term,
its original source declaration, exact frequency grid, basis and physical units.
The original full structural indices remain fixed under modal truncation.

Each duct retains two distinct physical unknowns: flow from cavity into duct
and flow from duct into exterior. The native `fs_duct::compose::TwoPort`
scattering equations enter the same simultaneous structure/pressure/flow solve.
Only inner flow enters the cavity weak equation; only outer flow excites the
joint exterior. No division by transmission, duct admittance, cavity resonance
or pressure-node determinant occurs. Native decaying-wave construction retains
finite-length resonances, compressibility and attenuation without a growing
ABCD transfer matrix. Each native cylinder is prepared once and its original
response reused over the declared exterior samples.

The equilibrium constructor uses native Bessel viscous/thermal propagation.
The molecular constructor uses the original `MolecularDuct`, retaining the
complete `RelaxationState` and explicitly chosen wall-diffusion resolution.
It binds that same state to cavity compression or requires exact equality with
an existing law; either attachment order refuses equilibrium substitution or
replacement of the molecular reservoirs. `Some(wall)` retains the original
joint thermal/molecular law and its admitted frequency band; `None` explicitly
selects adiabatic slip walls. No parallel equilibrium-Cp thermal shunt or
empirical molecular attenuation is added. Original scalar cavity-loss exclusion
continues to apply to molecular cavities.

Controls independently declare the reference transverse domain
`omega * radius / sqrt(R_specific * T) <= maximum_isothermal_ka <= 1`, native
entrance-composition conditioning, and scattering-face/work residual tolerance.
There is no compact axial `k*length` ceiling. Optional native inner-end inertia
is still exclusive to ONE uniform pressure mode with reference radius ka <= 0.25;
a resolved multimode cavity
uses no added entrance correction. The exterior owns the outer end exactly
once. This is frozen, linear, quiescent, plane-pressure gas propagation through
straight mechanically rigid ducts. The reference transverse bound and native
sampled wall-error admission are approximation domains, not a continuum error
certificate. Transverse modes, mean flow, evolving gas, nonlinear edge jets,
and compliant or spatially varying duct walls are outside this attachment.

The response exposes both flow vectors and native two-sided `TwoPortPower`
observations. Native right-face inward flow is the negative of exterior flow.
Original scattering fields are checked against solved physical face pressures
and flows, and absorption is checked independently against signed
`Re(p_inner conj(Q_inner) - p_outer conj(Q_outer))/2`. Materially active or
inconsistent work refuses; signed numerical roundoff is retained. Internal
viscous/thermal/molecular duct absorption is debited once as `opening_necks`;
joint exterior work is debited once as `radiated`. Cavity outward opening work
is an internal transfer. Modal convergence compares inner and outer flow
vectors separately in addition to the existing modal response metric.

Focused verification: eight native physical tests passed with the pinned July
2026 compiler using real `fs-math`, original `fs-la` complex LU, native
`TwoPort` constitutive/composition/power methods, and the exact numeric source
leaves. Four existing compact tests still pass. New tests cover independent
quarter/half/full-wave continuum transfer under three scattering references
with and without the original series entrance,
opaque transmission with finite absorption, signed cavity-basis invariance and
complete mutual wall/cavity/duct work, and truncation at original exterior
indices while retaining both end flows. Three additional public integration
tests cover native Bessel/Rayleigh/structural relaxation composition, the same
molecular duct and complex cavity compliance, and physical/source admission.
Those public tests were source/API reviewed and syntax parsed, but full Cargo
execution remains unavailable in this session because sibling workspace
sources and writable dependency/build storage are unavailable. The native
source-slice result does not claim a whole-workspace build.

## Adaptive simultaneous moving-gas molecular networks

`molecular_exchange::network::trajectory::advance_candidate` evolves the
original bounded network transaction over one requested physical window.
Every full/half trial calls `network::step_candidate` on its actual finite
gas inventories, temperatures, compositions and velocities. All apertures
sample that same original snapshot before reconstruction; simultaneous
incoming material cannot fund outgoing mass. Only accepted two-half-step
histories survive. Native maximum duration caps trials, and only the original
typed gross-withdrawal refusal permits physical retry. Source, drift, caloric,
closure and cancellation failures refuse the complete candidate.

The original pair trajectory's gas-state control type is reused. Local error
includes temperature, pressure, each original catalog species' mass (including
entry from exact zero), and velocity increments rather than common drift.
Each link/species also controls both gross directions and the authoritative
signed net mass, vector momentum and reference energy separately. Opposing
streams cannot hide their individual integration errors. Stream moments use
the explicitly declared common translating energy frame; normal directions
and aperture velocities are held in their original inertial axes.

Each species' fixed caloric reference is selected from a real positive initial
inventory in canonical gas-ID order. Active gases must preserve one original
caloric source per species. The original NASA caloric integrals and retained
molecular crossing corrections remove arbitrary formation datums before
arithmetic. Fixed caloric reference selection uses a real initial inventory rather than an
absent species at a receiver's window-start temperature. The original per-step
enthalpy crosswalk and its source-domain admission remain binding. Signed net stream energy retains the native source-relative
coordinate, including small-drift net mass, instead of subtracting rounded
gross formation energies. Closed apertures leave their source, normal and
drift inactive; isolated gases remain exact.

Every accepted interval audits actual original-to-current total/species mass,
inertial vector momentum and source-relative internal-plus-bulk energy in both
common and laboratory frames. Kinetic changes reuse the original changing-mass
momentum identity. Fixed window budgets also debit absolute per-gas native
defects, pooled impulse/energy defects, gross/net mass uncertainty and actual
accepted-step closure. Opposite signed errors cannot refund these budgets,
and no allowance grows with step count. Canonical compensated link totals and
native chronological receipts retain all accepted transfer accounts. Attempts,
accepted intervals, native gas/link/species caps and checked history cardinality
bound work and storage; cancellation is polled through reporting as well as
physical trials. Refusal publishes no partial state or history.

This controls local first-order temporal error of the existing fixed-volume
aperture operator. It does not supply other-wall pressure/support forces,
container motion, molecular flight storage, spatial gas gradients, a closed
vessel equilibrium/passivity claim, or a certified global error bound.

Six focused G0/G1/G3/G4 groups cover entering-species propagation and exact
native-history replay; convergence to an independent continuous three-gas
Gaussian-half-space moment/RK4 reference; stream-only refinement with loose
endpoint tolerances; formation-reference/common-boost invariance; fixed window
budgets and late atomic cancellation; and native withdrawal, source, work,
minimum-step, zero-area and empty-network behavior. Before an executor outage,
pinned Rustfmt parsed and formatted the new production and test files, and
the parent-module whitespace check passed. Recovery recreated the three new
files from their complete original patches, including the final reporting
cancellation polls; their recovered formatting was not rechecked. These
cross-crate native regressions were authored and source/API reviewed but
were not compiled or executed; no repository-wide DSR verdict is claimed.

## Causal shell material heat into original layered caloric cells

`laminate_thermal::harmonic::CausalLayeredShell` constructs the original CST/DKT
shell from complete resolved tensor material sources and the existing
`BoundLayeredThermalBody`. Every thermal subdivision becomes an original
mechanical sublamina. Exact cards and fixed query contexts, validity-wide
temperature independence, patch area, cell volume and reference mass must
agree; coincident authored coefficients do not establish source identity.

`heat_candidate` advances the general causal shell observer once, integrates
its actual positive quadrature heat per original thermal cell, and supplies
that energy to the existing nonlinear
`LayeredThermalStep::with_cell_energy_j(...).solve(...)` owner. Thermal
iterations reuse this fixed heat input. Original memory and caloric endpoints
are returned together only after source, complete DD/f64 clock, handoff and
combined-energy gates. Refusal or cancellation leaves both accepted histories
unchanged, including a late thermal-solver failure.

The receipt separates original material stress work, direct storage change,
physical heat, deposited cell energy, other external thermal heat and the
combined constitutive/caloric residual. Generated heat is an internal transfer
and is removed once from the external-heat account. No Kelvin/Maxwell partition
is invented when the material controller provides combined dissipation.

The thermal model resolves its original through-thickness cells and homogenizes
lateral temperature over their patch area. Interval heat becomes mean power in
the original backward-Euler residual; temporal refinement is needed to resolve
beat-driven transport. Mechanical coefficients, reference geometry and thermal
strain remain frozen; reciprocal temperature feedback or force equilibrium is
not inferred. `tests/laminate_harmonic_heat.rs` contains independent startup
heat and three-cell thermal references plus retained-memory, source, late
failure, retry and cancellation cases. Source checks do not establish executed
runtime validation.


## Endpoint source-density natural convection on finite enthalpy cells

`LaminarEnthalpy::advance_boussinesq` adds natural convection to the original
joint BDF1 momentum/rheology and finite enthalpy endpoint. Each actual thermal
and objective-shear query retains its original coherent `ThermalFluidTransport`
owner. The shared native `BoussinesqMaterial::try_new` and `BoussinesqConfig<2>`
admit that source's density/Cp departures and preserve reference rho0/Cp with
the SAME source's local dynamic viscosity and conductivity. The unprojected
thermal source, original source object/receipts, and projected numerical view
remain separately readable. No second density fit or material adapter is
invented.

The original triangle geometry, reference fluid mass and identity-bearing
caloric law remain unchanged. That native caloric law still determines h(T)
and stored energy; projected transport Cp never manufactures or replaces it.
The caller associates the source with that caloric model and selects physically
small departure bounds. These bounds admit the declared Boussinesq
approximation; they are not a quantified modeling-error or validity certificate.

Only `(rho_source - rho0)*gravity` is added to momentum. The held additional
force must exclude buoyancy and the removed hydrostatic reference load.
Boundary traction and returned pressure/reactions use
`p_reduced = p_physical - rho0*gravity.(x-x_reference)`; physical absolute
pressure/stress reconstruction requires the caller's explicit hydrostatic
datum. The original native cell-force entry prepares viscosity and buoyancy
together for each candidate and again for the final nonlinear residual,
including original constrained reactions. No coordinate lookup or lagged
thermal force is substituted.

The existing sign-split mass/enthalpy transport, original bulk/slip heat,
bounded nonlinear work and independent local/global thermal acceptance are
shared with ordinary `advance`. Source gravity does not appear as a second
thermal source. `report.buoyancy` separately retains the accepted force of each
original cell, actual density/Cp departures, and the physical volume integral
of force dot accepted affine BDM velocity. Its reported work is that endpoint
power times the actual BDF1 interval. It is not gravitational-potential change,
viscous heat or a whole-system passivity certificate.

Source bounds, inconsistent thermal views, invalid/unrepresentable forces,
native momentum/thermal failure, budgets and cancellation preserve both input
histories. The ordinary fixed-density entry retains its public signature and
thermal gates. Its previously invalid material-error variant was corrected to
the original native `MaterialError::State` as a concrete compilation blocker.
This remains fixed-mesh/reference-inventory, first-order Boussinesq coupling,
without thermal expansion/remeshing, variable-mass transport or an inferred EOS.

Independent bounded source/API and energy-account reviews checked the original
projection, native force coupling and thermal accounting. Final compatibility
checking corrected the import to the existing PUBLIC
`fs_lbm::material::transport::boussinesq` path; no dependency or feature was
added. The coupled regression groups specify zero-departure compatibility,
endpoint heating with both signs of independent buoyancy work and unchanged
native caloric storage, and source/departure/cancellation atomicity.

After executor recovery, parser/formatter passes completed on the coupled
production files and new regression child, and focused whitespace checks
passed. The parent test module changed only by adding its child declaration.
The native hydrostatic/circulation/refusal tests and coupled FEM/thermal tests
remain UNEXECUTED. An earlier bounded genuine-source compiler run reached the
native force entry and exposed two API errors that were corrected, but did not
produce a complete native artifact. Recovery removed its temporary numerical
dependencies; no broad rebuild or passing integration run is claimed.

## Conducting finite channel walls in molecular vessel recoil

`GasRocketState::advance_effusion_channels_with_conduction` connects the
existing finite-channel recoil transaction to prepared material conduction
links. The original no-conduction API delegates with an empty link set. The
material owner samples all molecular and conductive rates from one original
snapshot, pools both heat contributions by physical wall ID, and performs one
original caloric inversion per affected wall. Walls without an outlet can
participate in conduction; every distinct inventory still fits within the
vessel's existing dry mass. Original source identities, aligned distances,
contact resistance, area, and conductivity support remain authoritative.

Conductive heat is an internal transfer. The native accepted gas/wall energy
account enters the existing complete vessel balance, while only escaping
molecular momentum drives recoil. No conduction force, second thermal energy,
or invented external work is added. Its physical effect on molecular escape
appears when subsequent accepted steps sample the evolving wall temperature.
All gas, wall, and vessel endpoints remain one atomic candidate.

The model retains fixed geometry, fixed wall composition and reference volume,
reversible isochoric sensible storage, and quasistatic collisionless transport.
Each individual bore remains isothermal at its associated wall inventory's
temperature; the conduction graph resolves heat between separate inventories,
not a temperature profile along a bore. This does not add wall expansion,
phase change, rotation, or particle residence time.

Three focused G1/G3/G4 integration tests specify nonlinear original wall
storage under conduction alone with zero additional force/work; pooled
molecular and conductive heat with increased subsequent escape energy and
momentum; and source-basis mismatch plus late cancellation with deterministic
retry. The workspace implementation passed pinned Rustfmt before executor
loss. Source recovered from its exact focused edits has received static review;
its formatting was not rerun after recovery. These material and mechanics
integration tests have not been compiled or executed in this environment;
no repository-wide Cargo/DSR result is claimed.

## Adaptive finite-channel recoil with conducting finite walls

`rocket_blowdown::effusion::channel::trajectory` adds
`GasRocketState::advance_effusion_channels_adaptive` and
`advance_effusion_channels_with_conduction_adaptive`. Both use one private
multiwall controller and the original finite-channel/material/mechanics step
for every full and half trial. Accepted halves supply the next actual gas,
wall and vessel states; no source flux, geometry or accommodation solve is
duplicated to construct their impulses. The original thermal/recoil tolerance
type is reused as `EffusionChannelAdaptiveSolve`.

Local error control includes original species inventories, gas temperature and
pressure, drift-relative vessel motion, every wall's temperature and storage
increment, each outlet's separate gas/wall/escaping impulse, accommodation heat
and source-relative escape energy, total inertial exported moments, separate
molecular/conductive wall heat, and each directional conduction-link heat.
Opposed channels or opposing physical heat paths cannot hide an inaccurate
individual transfer behind their vanishing net result. Conduction-only wall
inventories remain in the same state and error-control map.

Escaped reference energy uses actual accepted species mass, original caloric
integrals relative to the fixed window-start temperature, the molecular RT/2
crossing correction, and original wall accommodation heat. It never subtracts
large gross formation energies. The window thermal audit independently
reconstructs original-to-final gas and every material-wall storage; internal
wall conduction cancels over that complete inventory. Original variable-mass
mechanical defects and remapping/gravity reserves remain explicit. Actual
endpoint momentum is checked against accumulated original exported momentum
and external impulses with compensated arithmetic. Species, mass, momentum
and energy allowances stay fixed over the complete requested duration; sums
of absolute defects and reserves cannot buy acceptance through cancellation.

The existing dimensional time, species and drift-relative motion helpers own
admission and exact final-time handling. Attempted intervals and retained
accepted half steps are explicitly bounded. Only the reachable original
temporal refusals allow halving: gross withdrawal, candidate gas caloric
bracket/convergence failure, or finite-wall inversion/overshoot. Material source,
geometry, arithmetic and physical-balance errors remain hard refusals. Failed
windows publish no gas, vessel or wall state. Export tolerances are specific to
the declared inertial frame; an observer boost need not preserve the adaptive
partition. These are local first-order numerical estimates, not certified
global or physical-model error bounds. The existing collisionless, quasistatic
channel, finite isochoric wall and translating point-vessel qualifications
remain binding; no entropy is assigned to the absorbing vacuum.

Five G1/G3/G4 integration test groups cover independent continuous monatomic
partial-accommodation gas/wall/recoil equations; opposed shared-wall mixture
outlets, beam-driven refinement and exact accepted-history replay; nonlinear
conduction-only wall convergence and overshoot retry; formation-reference and
common-boost molecular exports; and fixed-window, work, minimum-step,
unsupported-source, inactive-wall and cancellation refusals. The completed
implementation passed pinned rustfmt parsing/checks and scoped whitespace
checks before the executor became unavailable. The durable recovered source
retains its pre-format layout; no final formatting pass is asserted for those
reconstructed blobs. Source/API/physics review is supplementary; these
cross-crate integration tests have not been compiled or executed, and no full
Cargo or DSR pass is claimed.

## Calibrated torsion release and plastic heat in finite receivers

`friction_heat::adhesive::snap::torsion::JkrTorsionHeat` carries the original calibrated torsion lifecycle into exact registered finite caloric receivers. Every ordered normal or twist segment is regenerated from sealed source history. Normal snap, remaining twist-spring release and irreversible torsional loss keep distinct ledgers. Release and plastic loss require independently declared surface/outside partitions and named noncaloric destinations; reversible storage and signed normal/stiffness work stay mechanical. The heat ledgers begin at binding without resetting the loaded spring or converting past source losses again.

One caloric proposal admits all requested heat shares together through the original thermal impulse operator. It preserves actual receiver identity, material objects, body, mass and geometry. Compensated source-channel sums and the original relative caloric audit supply the signed transfer residual; separate partition, accumulation and caloric defects supply its absolute gate. Final callback, budget, source, caloric or stale-history refusal leaves source history, all loss accounts and calorics unchanged. A marginal attached state is admitted with zero twist without inventing a snap.

This base path is a same-time frozen-source impulse split. It supplies no conducting area, temperature feedback, coupled mixed-mode stability or native body impulse. Normal and retained spring release follow the original prescribed mode-I fold and calibrated homogeneous twist model, with its existing mode-III applicability ceiling.

Four focused integration cases cover independently derived fold twist energy and normal release, fresh recontact and retained prior loss, elastic/area work versus plastic heat, finite caloric refusal, and final cancellation/stale history with deterministic retry. Pinned `rustfmt --check` passed. The new cases were not compiled or run: DSR and RCH are unavailable and local build memory/disk headroom is exhausted.

## Resolved liquid thermal cavitation and its causal acoustic source

`cavitation::bubble::vapor::conducting::MaterialConductingVaporBubble`
binds the original water source to the native resolved-liquid vapor owner and
the existing compact observer. The native owner retains the entire moving
liquid temperature profile with phase mass, vapor temperature and radial
motion on one adaptive trajectory. Every stage receives its conductive
interface heat from that profile; this composition supplies no replacement
bubble heat coefficient or separately advanced thermal field.

Compact sound uses the actual displaced-liquid flow `4*pi*R^2*U_liquid`
and its original analytic derivative. During phase exchange this differs
from the geometric bubble-volume rate. Receiver density and sound speed
come from the same original water source. Observation introduces no added
radiation damping or acoustic backreaction to the incompressible model.

The direct, compact-source and observed advance paths retain the same
complete native proposal. A mismatched medium, late observer refusal or
final cancellation leaves every liquid cell and the complete phase/radial
history unchanged. The explicitly authored outer bath carries the finite
grid's conductive heat and advected enthalpy. Native source-property,
pressure-heating, capacity, outflow-Peclet, phase-equilibrium, balance and
radial/Mach restrictions remain binding.

`conducting_vapor_cavitation` is a runnable example under the existing
`cavitation` feature. It follows a 100 micrometre vapor bubble in 373.15 K
water through a 40 microsecond sinusoidal pressure cycle of 1 kPa amplitude.
It retains a graded 40 micrometre liquid layer, samples a receiver at 0.1 m,
and compares 16, 32 and 64 liquid cells plus a halved maximum time step.
Its TSV output includes interface and first-cell temperatures, signed
phase transfer and interface heat, causal receiver pressure, viscous heat,
and independently reconstructed cumulative phase/liquid thermal and mass
balances. Heat and phase transfer both reverse direction during the cycle.

All four new coupling groups and five existing vapor-coupling groups passed
in focused builds of the actual sources. The new groups compare the entire
adapter result against the native thermal trajectory, verify finite liquid
heat changes phase inventory and sound, and exercise observer refusal and
final cancellation with all cells retained. The example completed all
640 samples. At 64 cells it gave a sampled peak receiver magnitude of 1.48027 Pa and
peak interface heat magnitude of 16.6006 mW; the maximum independently
reconstructed cumulative thermal and mass residuals were 1.80e-18 J and
3.25e-26 kg. Changing 32 to 64 cells changed those pressure and heat peaks
by 0.1413% and 0.1530%; halving the time-step ceiling changed them by less
than 1.5e-13 and 2.8e-9 relative. These are case-specific numerical
comparisons, not general accuracy bounds or experimental qualification.
The complete workspace was not executed.

## Physical-time transport of calibrated torsion heat

`friction_heat::adhesive::snap::torsion::timed` advances the complete original torsion heat owner and finite calorics over a positive physical interval. The accepted contact and thermal clocks must match exactly. Original scheduled events at interval start precede the prescribed normal/area stage and its normal and retained-spring release. The source-regenerated monotone twist step then supplies only its plastic loss as interval-average power, with the independent declared torsional partition.

The existing native thermal solver owns conduction, reservoir cooling, ordinary sources, scheduled changes, caloric domains and bounded transport work. A private surface-power observer now serves both rolling and torsion. It retains the original receiver admission and compensated per-cell integration over actual accepted durations; changing clocks or rounded average watts must pass the explicit transfer gate, without an endpoint correction. The rolling extraction preserves the original arithmetic and callback order.

Original relative thermal audits and independently checked contact-transfer defects close the combined account without subtracting large rounded formation enthalpies. Absolute source, caloric and transport defects remain separate. All three loss ledgers, complete retained contact state, finite cells and accepted end time publish together after the final checkpoint. A late budget, domain, clock/history or cancellation refusal discards the whole interval.

This is a prescribed-normal, frozen-source split with interval-average plastic power. The original homogeneous torsion law and mode-III applicability ceiling remain in force. Temporal refinement is required for within-step dynamics; reciprocal mixed-mode stability, resolved within-step twist and temperature feedback are separate physics.

Three focused integration cases compare against independent discrete conduction/cooling modes, verify original normal and fold spring release before a scheduled cooling change, and exercise transport/caloric/cancellation/stale-clock rollback with exact retry. Pinned formatting and recursive Rust syntax/module-path checks passed. Integration compilation and execution remain unavailable under this session's runner and memory/disk limits.

## Direct sparse tetrahedral electrodes in finite thermal feedback

`SpatialConductionDrive::Terminals` and
`SpatialConductionResponseDrive::Terminals` accept the shared explicit
`ConductionTerminalDrive`. Every trial refreshes original bulk/contact sources
and closes held voltage/total-current/floating controls with one native sparse
field solve per constitutive evaluation. Existing matrix drives retain their
explicit N-1 response/budget semantics. Sparse workspace limits remain separate
from the common field-trial budget.

The retained terminal report owns physical voltage differences and independently
authored source work, including under large common bias; work is never rebuilt
from rounded absolute potentials or inferred from heat. Native volume/contact
work is deposited once through original receivers and contact partitions.
Existing caloric/thermal field-independence gates, phase/history ownership, late
cancellation and transactional publication remain binding.

Direct electrodes also compose with prescribed motional EMF through
`step_moving_candidate` and `step_moving_contact_candidate`. Each trial closes
original charge continuity and electrode-total constraints together in one
native sparse solve. The retained report supplies signed authored electrical
input; separate Lorentz output is subtracted before the finite-heat balance.
Only conductor-frame Joule loss and admitted stationary interface loss enter
thermal cells. Open circuits may retain a nonzero motional voltage at zero
current/heat, and generation may return negative terminal work. Original
fixed-occupancy motion limits and refusal of moving/sliding resistive interfaces
remain binding.

Focused regressions extend the analytic current/thermal and nonlinear
finite-contact oracles and add bounded-current material support, floating
electrodes, large common bias, two-step continuation, exact field-solve
accounting and rejection paths. The moving-terminal oracle also covers both
current directions, changing R(T), open circuits, generation, huge bias and exact
work-budget/late-cancellation behavior. These are manufactured numerical checks,
not measured material validation or mesh/time convergence claims.

## Driven native torsion and heat on one accepted clock

`friction_heat::adhesive::snap::torsion::timed::native` advances the original finite-body torsion lifecycle and complete torsional heat owner over one positive physical interval. Actual held COM forces and body-frame torques determine native motion. The internally accepted twist increment supplies the source-owned plastic loss; authored work or edited reports cannot create heat. The complete retained successor, original normal segment, angular constitutive step and irreversible loss must match exactly between native and thermal proposals.

Normal and remaining twist-spring release enter at the prescribed initial normal stage, after original scheduled thermal events at that time. The independently partitioned plastic heat is transported during the interval by the original finite-cell solver. Detached body work remains kinetic, detached travel creates no future spring preload, and recontact preserves accumulated source accounts while starting a fresh junction.

The combined energy audit uses actual native kinetic change, original incremental twist storage and body work together with the heat owner's compensated relative audits. Internal plastic transfer cancels once; independent mechanical and thermal absolute defects cannot conceal each other. Final caller checkpoint and cancellation checks precede publication of both momenta, native source history, all heat ledgers, finite calorics and thermal time. Late transport, caloric, budget, source, history or energy refusal leaves every original owner unchanged.

The source retains its homogeneous torsion calibration and mode-III applicability ceiling. Poses/frame and material properties remain held; normal approach is prescribed at interval start. Inertial midpoint motion and average-power heat require temporal refinement. This composition supplies neither normal body integration nor reciprocal mixed-mode stability or temperature feedback.

Three integration cases cover eight driven intervals with independently reconstructed body/contact/cell/outside/background energy and changing normal area, analytic fold release with detached work and fresh recontact, and post-release transport refusal, both final cancellation forms, stale time/history and deterministic retry. Pinned format and recursive syntax/module-path checks passed. These new integration cases were not compiled or executed because the required remote runners are absent and local build headroom is exhausted.

## Electrical contact heat in thermally live actuator supports

`ElectricalContactNetwork::field_motion_support_candidate` connects the original
source-pinned series contacts to the same thermoelastic member, native field/body
descriptor and finite enthalpy residual used by `MaterialNetwork`. Actual trial
temperatures refresh interface resistance; declared A/B partitions deposit each
contact loss once. Bulk Joule heat, reversible support heat, field work and
external supplies remain distinct, and accepted reports retain the heat routes.
Source, work-budget or cancellation refusal publishes no partial history.

Interface area and heat partitions are fixed within a step. The original
constant-reference-density solid support, affine modulus, scalar expansion and
constant field-coefficient restrictions remain. This adds no opening, arc,
contact mechanics, interface storage or strain-dependent electrical geometry.
Endpoint thermal feedback remains generally first order in time. Focused G1/G4
regressions exercise per-cell and total energy, source-dependent support motion,
source continuity, budget refusal and final cancellation with exact retry.

These regression definitions were parsed but not compiled or executed in the
implementation environment; full coupled validation remains outstanding.

### Thermal-only nonlinear heat assembly

`contacts::heat::thermal::step_heat_assembly_candidate` resolves a thermal-only
`HeatAssembly` through the existing cell-energy backward-Euler/Newton--Krylov
owner. All nonlinear probes and fresh final admission re-evaluate the complete
assembly at actual trial enthalpies; reciprocal finite-link conductances and
prescribed source powers enter once. Electrical bindings and coefficients
refuse. Public controls bound thermal unknowns, assembly evaluations and
numerical work without requiring circuit state or controls.

Original phase, pending-event/time, cancellation, local/global thermal and
closed-enclosure energy gates remain authoritative. Reports retain the accepted
source sample, native thermal/nonlinear accounts and external heat with internal
radiative net removed once. The existing first-order numerical-Jacobian method
remains subject to independent temporal refinement.

