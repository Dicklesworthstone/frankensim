# fs-thermal contract

Layer L3. Shared conservative thermal transport and solid/vapor exchange.
Runtime dependencies: std, fs-la, fs-math. Spatial discretizations and material
constitutive data remain owned by their existing crates.

## Existing thermal timeline

The enthalpy, advection, contact, observation and event APIs retain their
documented control-volume semantics. A cell stores total energy as mass times
specific enthalpy. Every material chart in one balance uses a compatible energy
datum. The shared timeline is re-exported by fs-conduction; consumers must not
implement independent copies of its storage update.

`DepositionSimulation::apply_heat_j` adds signed endpoint heat to those same
cells without advancing time or re-integrating sources and contacts. The
original material law resolves `h_new = h_old + Q/m`, including latent storage.
Same-time scheduled events must be processed first. The method checks actual
cell-relative stored energy against the supplied joules, with separate local
and aggregate absolute-error tolerances; opposite lost increments cannot hide
behind zero net heat. Compensated sums preserve small changes beside large
caloric reference energies. Domain, arithmetic and cancellation failures return
no candidate. The energy audit identifies the input as external to the thermal
subsystem; a coupling owner must account for the corresponding upstream debit.
No mass, geometry, boundary, event cursor or accepted time changes, and no
mechanical-to-heat conversion fraction is inferred. Five focused regressions
cover latent/signed inputs, unchanged clocks/sources, large reference datums,
local/aggregate loss gates and transactional refusal.

### Enthalpy transported by balanced mass flow

The existing upwind graph evaluates transport relative to one canonical
start-state enthalpy per connected flow component. Fixed mass makes the common
reference cancel analytically; its numerical value cannot create heating or
erase a small source by subtraction of gross formation-energy powers. Component
identities are prepared once, and the same relative fluxes feed thermal and
optional entropy balances. Latent enthalpy still moves at equal temperature.
The gross absolute transported-enthalpy diagnostic retains its original datum.

The existing cell mass-balance admission remains 64 machine epsilons of the
incoming/outgoing rate scale. Inside that allowance the equations enforce zero
continuity; the time-integrated absolute cell mass residual and its
`|h_reference * mass_residual|` energy equivalent are reported explicitly as
roundoff diagnostics. Neither is a physical heat source, mass loss, or second
energy debit. G1/G3 regressions cover a flowing calorimeter at large positive
and negative reference offsets, small heating in a uniform circulating loop,
disconnected latent-heat loops, canonical edge orientation, admitted mass
roundoff, and matching energy/entropy reference invariance. This does not add a
pressure, velocity, variable-mass or chemical-reaction closure.

## Sublimation and deposition

`sublimation::SublimationLaw` supplies a bounded immutable saturation-pressure,
latent-enthalpy, density and molecular-mass chart. The material crate supplies
reference ice/CO2 correlations and an explicitly parameterized
Clausius–Clapeyron approximation. A numerical correlation is not a sourced
specimen certificate, and supplied density/accommodation/film coefficients
retain their caller-declared applicability.

Positive mass flux means solid becomes vapor; negative flux deposits vapor.
The input pressure is the volatile partial pressure. Outgoing molecules use
surface temperature and incoming molecules use gas temperature. One common
accommodation coefficient preserves detailed balance at equal temperature.
The optional dilute gas film is combined with interface kinetics by series
resistance. Zero accommodation or film coefficient seals the surface.

`surface_flux` returns signed mass and latent-energy fluxes, their analytic
surface-temperature derivatives, and outward-normal interface velocity.
`isothermal_mass_transfer` integrates an ideal, thermostated finite vapor
reservoir exponentially at fixed surface temperature. It uses the shared
deterministic expm1 implementation to preserve short timesteps and does not
overshoot its frozen-temperature equilibrium. Its caller supplies the thermal
work needed to hold the surface at that temperature.

`advance_surface` couples that transfer to the caller's original T(h) chart
with a bounded implicit endpoint solve. The endpoint vapor carries h_end + L.
The energy equation is:

    m_initial * (h_end - h_initial) + transferred_mass * L(T_end) = 0

This is algebraically the balance of remaining solid energy and transported
vapor enthalpy. The latent account is part of vapor energy, not another
independent energy debit. All mass limits are applied before deriving heat,
vapor inventory and recession. A shifted common caloric datum leaves the
physical temperature and mass solution unchanged.

The nonlinear solver requires a continuous nondecreasing caloric chart and a
bracketable endpoint in the declared constitutive domain. It refuses
unsupported temperatures and failed convergence rather than clamping them.
The default budget permits 96 evaluations; callers may select 1 through 4096.
Each evaluation polls the supplied checkpoint. Functions return candidate
values and mutate no caller state on failure.

## Numerical and physical limits

The coupled endpoint is a first-order reduced model. Area remains frozen over
one phase step; recession uses endpoint solid density. Finite gas temperature
is prescribed: its thermostat owns subsequent sensible-heat equilibration.
The dilute film does not model non-dilute Stefan flow, pore resistance, plumes,
curvature, nucleation, melt competition, radiation or unresolved gradients.
Those mechanisms require explicit compatible models, not hidden coefficients.

Returned zero solid mass denotes extinction. A spatial consumer must update
topology, contacts and conservative remapping or refuse that endpoint.
Reporting recession does not itself move a mesh, alter a level set or change a
flow/electrical geometry map. Uniform-cell use requires a justified internal
temperature approximation. Refining time cannot repair an invalid spatial
model.

The scalar flux and finite-reservoir kernels allocate no memory and perform
constant work. The implicit solver uses constant memory and at most its
declared evaluation budget. This is an algorithmic statement, not a measured
hardware throughput or cross-ISA reproducibility claim. Transcendentals use
the existing fs-math owner; current executed evidence is same-ISA.

## Focused verification

G0/G1/G3 tests cover signed deposition, equilibrium, unequal gas temperature,
analytic tangents, kinetic/film limits, finite-reservoir exact relaxation and
semigroup behavior, tiny timesteps, inventory limits, implicit cooling and
heating, mass/energy/recession consistency, caloric-datum shifts, first-order
refinement, cancellation and exhausted nonlinear budgets. Material tests
compare the actual correlations to published reference values and their
analytic derivatives to finite differences. These tests establish the
declared reduced model; they are not experimental validation of a complete
sublimating object or a general spatial Stefan solution.


## Phase exchange on the existing heat timeline

`enthalpy::sublimation` composes the surface operator with the original
`DepositionSimulation`. `sublimate_cell` atomically changes the active cell's
actual mass and specific enthalpy plus its vapor inventory, without advancing
the clock. `advance_sublimating_to` performs bounded first-order heat/phase
splitting with one total thermal work budget and one global mass/energy audit.
It honors birth and boundary/contact times exactly, uses the exposed area from
the start of each interval, and never sublimates a not-yet-born cell. The
existing caloric chart remains authoritative; rejected out-of-domain trial
points do not invalidate a successfully converged in-domain endpoint. All
failures and cancellation preserve the original timeline and reservoir.

The whole cell must remain a pure solid on its original caloric branch.
Temperature alone cannot establish solid fraction on a melting plateau;
mixed-phase inventories require a separately owned partition. The supplied
finite reservoir belongs to this one surface and remains thermostated. A
joint chamber shared by several surfaces needs a joint backpressure solve.
Reference geometry/contact areas remain fixed and recession is reported.
Complete cell exhaustion refuses until a consumer owns topology removal.
There is no implied moving mesh, gas-energy solve or material calibration.

Vapor energy includes the solid's original sensible-enthalpy datum plus latent
heat; the separately reported latent account is never an additional debit.
Ten new native integration regressions cover balances, frost growth, covered
surfaces, birth/contact timing, nonlinear chart trials, transactional refusal,
whole-advance budgets and first-order refinement. The combined exact-source
thermal/material suite passed 95 tests with one existing performance test
ignored. `cargo run -p fs-material --example sublimating_ice` demonstrates the
real ice law with explicit illustrative calorimetry in a maintained vacuum and
a finite chamber; both native example cases passed mass/energy assertions.


## Several surfaces sharing one vapor chamber

`sublimation::chamber::advance_chamber` solves one shared endpoint vapor
pressure together with the original conservative surface enthalpy operators.
Its backward-Euler mass balance admits warm-to-cold redistribution, including
condensation of vapor emitted during the same interval into an initially empty
chamber. Depositing surfaces never receive independent copies of the original
vapor inventory. One volatile species with constant molar mass, pure-solid
caloric branches and compatible enthalpy references are required. Gas volume
and temperature stay prescribed.

A frozen-temperature predictor and safeguarded scalar pressure interpolation
avoid a dense coupled matrix; each trial costs O(surfaces), with explicit
surface, local-iteration and total pressure-trial ceilings. Unsupported outer
trial bounds are backed away from. The same accepted transfers reconstruct
final vapor mass and drive every solid mass, latent/sensible-energy and
recession account. The nonlinear pressure/mass closure residual is reported
separately from actual published mass conservation. Refining time controls the
first-order integration error; the residual does not certify model accuracy.

Eight new native regressions check independent simultaneous algebra, transient
refinement, stationary-pressure circulation, initially empty chambers, shared
inventory, area partition/order, refusal and finite-capacity energy closure.
The combined thermal/material suite passed 113 tests with one existing benchmark
ignored. Thermally evolving gas, molecular heat accommodation, mixed phases and
mesh/topology changes remain separate physical owners.

### Relative entropy on the shared thermal timeline

ThermalLaw may provide its exact relative integral of dh/T at the original
fixed-pressure/composition convention. Missing entropy remains unavailable;
no substitute heat capacity or entropy reference is invented. The opt-in
advance_to_with_entropy and BackwardEulerStep::evaluate_with_entropy use the
original accepted thermal states, actual interface coefficients and the same
start/endpoint flux quadrature. EntropyAudit separates exact material storage,
actual ambient/support reservoir entropy, conductive production, unclassified
source energy divided by cell temperature, and the numerical balance defect.
The source flag survives cancellation of signed sources. Birth entropy is
excluded on both sides; subsequent storage starts from actual birth enthalpy.
A trial audit is accepted only with its trial. Ordinary energy-only APIs do no
entropy integration. G1/G3 tests cover isolated exchange, refinement, reservoirs,
latent storage, equilibrium, deferred birth, resolved implicit sources, missing
material entropy and cancellation without publication.


## Shared chambers on actual multi-cell heat histories

`enthalpy::sublimation::chamber` applies one joint phase exchange to canonical
cell identities using each cell's original caloric law and exposed area.
`sublimate_chamber` validates every candidate before committing any solid or
vapor state. `advance_chamber_sublimating_to` composes original heat transport
with joint exchange, honors start-of-interval exposure and exact scheduled
birth/contact boundaries, and retains one whole-advance thermal work budget.
The report checks combined solid/vapor mass and the original total energy
account, while keeping each surface's recession separate.

Whole cells must remain pure solid; gas temperature and reference geometry
remain fixed. Complete extinction requires a separate topology owner and
refuses atomically. Ten native integration regressions passed for inventory,
redistribution, original caloric charts and heat graph, canonical order,
event timing, work ceilings and late cancellation. The runnable fs-material
example `ice_cold_trap` uses the real ice law to compute frost growth on a cold
film supplied by a warmer film through an initially empty shared chamber.


## Molecular counterflux for nonisothermal interfaces

`sublimation::molecular` returns both nonnegative half-space mass fluxes and
one signed net flux from the original bounded constitutive and kinetic owners.
The equal-temperature path retains stable pressure-difference evaluation.
Different surface and local vapor temperatures remain distinct; volatile
partial pressure is never replaced by total carrier pressure. Scalar work is
constant and allocation-free. Four native regressions cover equilibrium,
unequal temperatures, pressure semantics and refusal boundaries.

An energy consumer must retain the opposing streams: zero net mass does not
imply zero heat exchange. For Maxwellian ideal vapor with equilibrium internal
modes, each stream carries h(T)-R_specific*T/2 in an explicitly aligned energy
reference. The material gas adapter owns that caloric evaluation and both
endpoint energy accounts. This interface assumes energy/angle-independent
equal sticking and specular nonsticking reflection; it supplies no gas-film,
diffuse thermal accommodation, Knudsen-layer or moving-geometry model.


## Original thermal storage for external solid/gas coupling

`enthalpy::sublimation::external` exposes the original pending events, work
budget and cell caloric evaluation, and transactionally applies prepared solid
mass/enthalpy endpoints to the existing thermal storage. Every full expected
cell, positive remaining inventory, caller-admitted pure-solid h interval,
original caloric endpoint and original-geometry exposure is checked before a
single write. Endpoint contact events may change current exposure; they cannot
retroactively alter the preceding interval's declared area. A surrounding
coupled owner mutates its private trial without another graph clone per phase
step. Joint gas/solid conservation and physical phase admission remain with
that owner; this boundary does not invent calorimetry or geometry evolution.

The material gas timeline now composes this boundary with the actual original
heat operator and one jointly evolving gas. Four native owner tests cover
state preservation, stale/bad batches, endpoint covering and late cancellation.
Both existing isothermal sublimation timelines also jump directly to the next
event/target while the WHOLE graph is empty, avoiding unbounded tiny pre-birth
intervals. Two regressions prove a long pre-birth gap with zero thermal work
allowance. The native sublimation suite passed all 51 selected tests.

## Heat-driven condensed-fuel pyrolysis

`enthalpy::pyrolysis::PyrolysisSimulation` extends the original activated
thermal cells with separate virgin and char inventories. It retains their
identity, geometry, contacts, boundary conditions and extensive storage, and
uses the same conductive/boundary heat operator through a private caloric
callback. No second conduction graph or alternate transport algebra is used.
`CondensedCaloric` adds source-owned h(T) to the original `ThermalLaw::T(h)`;
admission checks that the charts agree. Virgin and char may have distinct
nonlinear sensible calorics. `VolatileEnthalpy` supplies the escaping gas's
physical enthalpy on a compatible formation reference.

One first-order non-oxidative reaction converts virgin material into an explicit
char mass yield and outgoing volatiles. Its accepted endpoint temperature
determines the backward-Euler inventory, retained char and outgoing enthalpy:

```text
k(T) = A exp(-E/(R T))
m_virgin_new = m_virgin_old / (1 + dt k(T_new))
extent = m_virgin_old - m_virgin_new
mass_out = (1 - char_yield) extent
energy_out = mass_out h_volatile(T_new)
```

The local energy equation balances the change in remaining solid storage plus
that exact outgoing enthalpy against supplied heat. Reaction enthalpy is already
in the compatible material references; it is never a second source or debit.
The independent virgin/char caloric changes and reaction enthalpy differences
are evaluated before gross formation energies are combined. Per-cell releases
retain their integrated mass, enthalpy and consumed virgin mass for downstream
gas coupling, even when emission temperatures change within one advance.

Tiny rate-duration products use stable loss fractions. Published consumption
is reconciled to actual retained inventory, so an unrepresentable decrement
cannot emit fictitious mass or heat. The absolute kinetic mass-roundoff and
reaction-energy equivalents are reported separately and checked against the
caller budgets; they are numerical discrepancies, not physical transfer streams.
Local, cumulative absolute and whole-advance mass/energy residuals cannot be
hidden by cancellation between cells. Errors and cancellation publish no
candidate. Boundary and conductance changes use original transactional admission.

This initial lane requires endothermic single-step kinetics, monotone sensible
caloric charts, explicit temperature/source domains and fully activated fixed
topology. It refuses pending geometry events, exhausted work budgets and
unrepresentable depletion. It does not seed residual solid mass. Conduction
retains the original explicit stability bound using the current condensed mass
and mixture capacity; local reaction is implicit in its final temperature.
Near depletion this fixed conductive mesh can become stiff. Moving surfaces,
porosity/permeability, pore-gas transport, oxidation, latent-phase state,
multi-step wood chemistry and fire-spread closure remain separate capabilities.
The coupled method is first order, with no unconditional stability or entropy
theorem. Conductivity changes require an explicit caller-supplied law.

Eleven focused regressions cover independent isothermal and adiabatic
solutions, distinct char calorics, first-order temperature/mass refinement,
equivalence to original heat transport with reaction disabled, spatial heat-
driven mass release, caloric-datum shifts, small-rate representation, and
domain/cancellation/work-budget refusal. The native full thermal library suite
passes 164 tests with one existing performance test ignored. This is focused
native Cargo evidence, not a repository DSR gate or experimental validation.

`heated_pyrolysis_slab` exercises the actual 1D graph with a resolved front
half-cell heat balance, irradiation, emission and convection. It uses retained
NIST/MaCFP scalar PMMA inputs in an authored heating scenario; its optional
20% char sensitivity is explicitly not a wood material fit. `macfp_tga` uses
the same kinetic update against the retained MaCFP FDS material-temperature
history. Source data, units, calibration/verification distinctions and commands
are in `data/reference/compartment_fire/pyrolysis/README.md`.

### Implicit conduction and finite-rate decomposition

`pyrolysis::implicit` removes the explicit conduction stability limit while
retaining the original physical model. Its nonlinear unknown is absorbed cell
energy. Every trial calls the existing `solve_cell` for temperature, virgin/char
inventory and escaping volatile enthalpy, then the original
`heat_rates_with_calorics` evaluates endpoint conductive/Robin/source powers.
The actual residual is absorbed energy minus duration times endpoint power.
No alternate heat graph, heat capacity, kinetic law or heat-of-reaction source
is introduced.

`implicit_workspace` prepares original dense LU and numerical buffers for at
most 128 cells. The Jacobian combines the original link/boundary conductances
with local caloric-response probes, avoiding a complete graph reevaluation for
every column. Original bounded heat-rate maps and candidate/result objects
still allocate; this is not an allocation-free whole-advance claim.
`advance_implicit_to` returns the unchanged original release/audit format plus
actual nonlinear, local-solve, linear-solve and backtracking work counts. Local
and whole-advance conservation gates include original kinetic roundoff and the
separately reported graph heat-transfer summation discrepancy. Cancellation or
nonlinear/domain refusal publishes no candidate.

Seven focused native tests cover an independent stiff two-capacity endpoint,
simultaneous heating/conduction/reaction, char and caloric-datum invariance,
actual Robin coverage, temporal refinement, activated Arrhenius comparison
with the original fine explicit method, and refusal/cancellation/reuse. The
stiff case completes in one implicit interval while the explicit owner exhausts
its four-substep budget. This remains first-order backward Euler with caller
time-resolution controls; removing the conduction CFL is not an unconditional
nonlinear-convergence, temporal-accuracy or entropy theorem. Geometry and
authored coefficients stay fixed within each original advance.

### Conservative noncharring slab recession

`pyrolysis::slab::RecedingSlab::new` adopts the original homogeneous x-aligned cell
chain with caller-sourced constant density, conductivity and cross-section.
Actual retained mass sets thickness through `L = m / (rho A)` with a fixed
back face. This is the single-material density/length relation in the pinned
FDS Solid Chapter's shrinking/swelling section; fixed-back anchoring and
conservative rezoning are explicit choices of this implementation.

Each bounded interval uses the original implicit conduction/reaction owner.
Excessive actual per-cell mass loss retries the complete unpublished interval,
without accumulating rejected volatile releases. Accepted specific enthalpy is
reconstructed with limited linear slopes in retained-mass coordinate and
integrated over overlaps onto equal-mass cells. The original caloric inverse
recovers temperature; original cell geometry and series conductances are rebuilt.
Mass and extensive enthalpy are transferred, not averaged temperatures. Cell IDs
are spatial slots after rezoning rather than tracked material particles.

Original volatile mass/enthalpy integrals are unchanged and retain their
production slots and homogeneous source binding. The original aggregate audit
includes actual remap storage/discrepancy once; separate mass/energy residuals,
rounded-geometry mass equivalent and enthalpy redistribution are returned.
The redistribution diagnostic is not lost heat or a PDE error bound. No seed
mass, heat-capacity floor or release clipping extends the physical domain.
Geometry and coefficient evolution remain first-order split from each implicit
interval. The charring extension below supplies separate residue geometry.
Pores, mechanical work, thermal expansion, gas-domain volume displacement and
complete burnout require separate models.

Seven native groups verify independent isothermal inventory/recession,
nonlinear-caloric overlap conservation, monotone reconstruction/refinement,
formation-datum shifts, reaction-disabled identity, complete retries and
cancellation/replay/source/work refusals. The combined native thermal suite
passes 178 tests with one existing performance test ignored; DSR is unavailable
in this execution environment, so this is not a repository-level gate result.

`heated_pyrolysis_slab` retains its original five arguments and fixed default;
an optional sixth `receding` argument enables the new noncharring path with the
actual front half-cell resistance. Original NIST/MaCFP scalar inputs and the
effective volatile enthalpy remain unchanged. The 120-second 24/48-cell runs
complete, with emitted-mass spatial sensitivity of about 1.88%; their actual
accounts are retained in
`data/reference/compartment_fire/pyrolysis/receding_slab_example_results.json`.
These authored heating cases do not establish experimental PMMA validation.

### Conservative charring slab shrinkage and swelling

`RecedingSlab::with_char` adds retained residue to the same finite-rate
virgin-to-char-and-volatile material, implicit heat solver and release accounts.
`CharProperties` supplies independent char density, conductivity, source and
a declared `ConductivityMixing::Parallel` or `Series` arrangement. The
existing `new` constructor remains the noncharring path.

The declared condensed volume is `V = m_virgin/rho_virgin + m_char/rho_char`.
The back remains fixed; front motion may be positive recession or negative
recession (swelling). This additive-volume choice is explicit, not the general
conditional FDS shrinking/swelling algorithm. Parallel conductivity is the
volume-fraction arithmetic mixture used by the pinned FDS Solid Chapter;
series conductivity is its harmonic counterpart for the separately declared
unresolved arrangement. All rebuilt internal links and the public
`cell_conductivity_w_m_k` and `front_half_cell_conductance_w_k` queries use the
same current component mixture.

Equal-retained-mass rezoning conserves total mass, char mass and physical
enthalpy, with virgin mass recovered as the complementary component. Char
fraction and sensible enthalpy use a paired limited reconstruction that
preserves donor means, nonnegative component inventories and the common
caloric-domain bounds. Formation references travel with their component
masses. The original scalar energy owner recovers mixed temperature with
reaction disabled during inversion; rezoning creates no volatile stream or
additional reaction heat. Separate virgin/char remap residuals prevent total
mass cancellation from hiding composition drift.

Frozen-geometry intervals are bounded by both the existing mass-loss limit
and `maximum_cell_volume_change_fraction` (default 0.25). This also controls
swelling with little mass loss. An excessive change retries the entire
unpublished interval, retaining no rejected release. Reports include maximum
volume change, component and total mass discrepancy, physical-energy
discrepancy, geometry roundoff and original caloric evaluation work.
Cancellation or failed budgets return no candidate.

Five additional native groups establish independent component inventory,
volume, fixed-back motion, outgoing enthalpy and both conductivity limits;
nonlinear-caloric remap and series-link consistency; better-than-threefold
remap refinement with independent phase energy-datum shifts; swelling retries;
and source refusal, cancellation and deterministic reuse. All 12 slab tests
pass, and the focused native thermal library passes 183 tests with one existing
performance test ignored. This is not a full-workspace DSR result.

Geometry/coefficient evolution remains split and first order in time, and
remap refinement does not certify global PDE accuracy or an entropy theorem.
No gas-volume displacement or pressure work is inferred. Pore transport,
permeability/porosity evolution, char oxidation, cracking, multi-step wood
chemistry, latent phases, complete burnout transitions and an experimentally
validated wood material remain outside this closure.


## Geometric conduction with evolving material conductivity

`enthalpy::conduction::ConductionObserver` prepares existing thermal links with
explicit centre-to-face distances, contact resistances and immutable scalar
`ConductivityLaw` bindings. It queries each participating cell once per actual
explicit substep and shares that value across incident faces. Internal links
reuse the original series resistance owner. A surface reservoir adds the
half-cell material resistance and an explicit areal film/contact resistance;
zero external resistance prescribes temperature at the face, not at the cell
centre. Zero conductivity explicitly insulates. Geometry, reference mass,
material identities and face bindings remain fixed and are checked on reuse.
Material selection and mixed-phase microstructure remain the constitutive
owner's responsibility.

The original observer protocol now accepts step-local reservoir conductances.
Both pairwise links and reservoirs enter the same original heat-rate and CFL
calculation. Their authored values are restored before accepted-step hooks and
events. No reservoir heat is reclassified as an applied source. Composed
observers reject conflicting reservoir ownership, and the existing advection
wrapper forwards the hook. `advance_observed_with_entropy` forwards each
physical hook once and retains the exact applied sources, links and reservoirs
for the original entropy account, including actual bath temperatures. General
internal transfers and applied sources retain explicitly unclassified entropy
contributions; temporal entropy defects remain distinct from physical
conduction production. Failure or cancellation publishes no partial physical
state.

The prepared graph uses O(cells + faces) storage and
O((cells + faces) log(cells)) work per substep with ordered maps. This is an
algorithmic bound, not a measured throughput claim. The existing explicit
coefficient-freezing and caloric CFL applies; it is not an implicit solver or
a time/spatial error certificate. Focused G0/G1/G3 tests cover series/Robin
resistance, dynamic boundary CFL and entropy, exact event splitting, independent
nonlinear transient and Stefan front/heat refinement, composition, insulation,
domain/work limits and transactional cancellation.


## Implicit enthalpy transport on geometric conduction chains

`enthalpy::backward_euler::BackwardEulerStep` additionally accepts optional
endpoint reservoir-conductance replacements through
`evaluate_with_reservoirs` and `evaluate_with_entropy_and_reservoirs`.
Original bath temperatures, areas, heat rates and entropy remain authoritative;
authored coefficients remain unchanged through evaluation and acceptance. Existing
coefficient callbacks retain their original behavior. Unknown cells, absent
baths, negative/nonfinite conductances and cancellation refuse transactionally.

`enthalpy::conduction::implicit::ImplicitConductionChain` solves fixed-mass,
nearest-neighbor chains using that same residual, not a second heat-rate owner.
Every caloric and conductivity law is bound to the original borrowed object.
A fused constitutive query supplies T(h), k(h), their analytic local derivatives
and an explicit admitted enthalpy interval. Physical latent plateaus retain
zero dT/dh while their conductivity may continue changing with phase inventory.
The exact endpoint Jacobian differentiates both the temperature difference and
the half-cell/contact/surface conductance, in cell ENERGY coordinates.

Only three matrix diagonals are assembled. The reusable `fs-la` adjacent-pivot
workspace keeps the linear solve O(cells) in work and storage; original
map-based graph/constitutive bookkeeping remains. Newton candidates whose
nonzero increment cannot change the stored enthalpy are represented by exactly
zero increment before evaluation. This does not bypass or relax the original
underflow policy, heat rates, per-cell residual or total-energy acceptance.
Every resulting candidate is reevaluated by the original owner.

Declared controls bound timestep, accepted intervals, Newton iterations,
backtracks, interval halvings and total residual evaluations, including
rejected trials and final entropy queries. Bounded backtracking respects the
material domains; nonconvergence or a singular linear system may trigger a
whole-interval retry. Scheduled events split intervals at their exact times.
Only fully accepted intervals contribute heat or entropy. A complete call
also checks independently summed original stored energy against all supplied
heat. Failure may alter scratch but publishes no physical history. The report
retains original heat/entropy, maximum accepted local residual and work counts.
Every admitted caloric law must provide the original entropy integral.

This path removes the explicit diffusion CFL restriction for an accepted
implicit interval. Backward Euler remains first order; interval halving here
is nonlinear globalization, not a time-error estimator. Global uniqueness,
unconditional nonlinear convergence, a general graph solve, evolving geometry
and constitutive accuracy beyond the supplied laws are not claimed.

Five focused chain regressions cover stiff exact solutions, latent traversal,
original-residual Jacobian differences, nonlinear time refinement, event
splitting, cancellation/budget rollback, singular falling-conductivity retries,
and a 96-cell diffusion tail with sub-ULP energy updates. Two added reservoir
regressions check original heat and actual-bath entropy plus refusal/rollback.
The full pinned-main native thermal source passes 191 tests, with one existing
ignored test; required workspace/DSR gates remain a separate unavailable check.

## Additional material conductances on the native explicit timeline

`DepositionSimulation::limit_timestep_for_additional_conductances` composes
nonnegative extra W/K with an already admitted explicit step using the original
cell masses and registered global caloric dT/de bounds. It adds the largest
extra decay rate harmonically to the current rate allowance, rather than taking
two independent CFL minima. An empty or zero map leaves the proposed step
bit-for-bit unchanged; zero caloric slope adds no restriction. Unknown cells,
invalid rates, unrepresentable steps and checkpoint failures refuse.

This helper applies no heat and changes no thermal state or authored boundary.
Call it before event clipping through the existing observer timestep hook.
It can be conservative when different cells or a caller maximum step dominate.
It bounds frozen heat-transfer coefficients, not their nonlinear derivatives
or time-discretization error. Focused native regressions cover combined
conductive/cooling limits, sequential composition, zero slopes and refusal.
The airflow material observer uses the exact hA from its start-state source
evaluation, so changing fluid data changes both cooling and native subcycling.

## Shared endpoint ambient and reservoir conductances

`backward_euler::{BoundaryConductances, BoundaryStepCoefficients}` pass endpoint
total ambient/reservoir G directly into the original heat-rate and entropy
owners. The boundary-aware residual and external-power APIs retain prescribed
temperatures, reference coverage, authored h/G and cell watts. Replacement G
is applied once; no boundary input is mutated or reclassified as generated heat.
The same G participates in signed physical watts, outgoing-conductance/CFL
accounting and boundary entropy Q/T. Existing ordinary callbacks use empty maps;
the reservoir-only APIs delegate their maps through this same general ingress.
Finite thermal stores remain original internal links, not prescribed baths.

Unknown cells, absent reservoirs, negative/nonfinite G and nonzero Robin G
without exposed reference coverage refuse. Three added native regressions cover
independent finite-cell energy/entropy/CFL, unscaled applied watts and empty-map
parity, invalid coverage and final cancellation/replay. These tests are compile
checked in this change; their Cargo runtime build was killed while compiling
the existing asupersync dependency. The earlier native test counts above are
prior validation, not runtime evidence for this extension. Geometry-model
validity and nonlinear/time accuracy remain the composing caller's responsibility.

## Declared thermal storage and closed-volume admission

`ThermalLaw::storage_basis` declares whether the original caloric law stores
specific enthalpy or fixed-volume specific internal energy. Existing laws retain
the enthalpy default. The inherited `specific_h_j_kg` field carries the selected
coordinate; the timeline applies the same conservative heat exchange and
`integral de/T` entropy account without converting one thermodynamic potential
into another. Each law can admit its own mass/geometry relationship through
`validate_cell`, called before a birth changes inventory.

`material_storage_basis` exposes the original registered law's declaration.
Nonzero prescribed mass flow and inlet calorics require enthalpy storage. The
single-surface, joint-chamber and external-endpoint sublimation paths also require
enthalpy before constructing or publishing open-mass exchange. Batch admission
finishes before any endpoint mutation. Zero-flow connections remain available
for heat-only internal-energy storage.

Focused regressions exercise two rigid bodies exchanging heat with conserved
internal energy and increasing thermal entropy; invalid mass/volume births;
forward/reverse circulation and boundary-flow refusal; and all three sublimation
entry paths. Mixed enthalpy/internal-energy batches preserve the original solids,
vapor, mass, stored energy and clock on refusal. This adds no open-mass internal-
energy model, pressure-work conversion or equation of state. Such conversions
remain the explicit responsibility of the physical material/equipment owner.

## Adaptive implicit thermal accuracy

`enthalpy::conduction::implicit::adaptive` adds first-order step-doubling
control around the original `ImplicitConductionChain`. One full backward-Euler
interval is compared with two actual half intervals; only the conservative
half-step trajectory is accepted. Comparison trials never contribute physical
heat, entropy or accepted time, and no Richardson-extrapolated state is formed.
Unequal representable half spans receive their first-order error correction.

`ImplicitConductionLaw::time_observables` supplies original temperature and an
optional original phase mass fraction. The per-cell estimator checks stored
specific energy, temperature and, when requested, phase inventory. Relative
scales use changes from the accepted state, rather than absolute enthalpy or
Kelvin values. Phase accuracy refuses when an original observable is unavailable
rather than inferring phase from temperature. The h-labelled coordinate retains
the declared thermal storage basis.

Controls retain the original local/whole-call joule gates and bound comparison
spans, comparison count, accepted half steps, consecutive reductions, Newton
work and the shared residual-evaluation budget. Scheduled events are hit exactly.
Source refusal, exhausted work, cancellation and unresolved accuracy return no
partial trajectory. Reports separate local time-error estimates from algebraic
energy residuals, count all attempted nonlinear work, and expose the next proposed
comparison span for continuation.

Four native regressions cover analytic cooling and caloric-datum invariance;
phase-driven refinement on an exactly isothermal latent plateau; non-dyadic bath
events and fine-only heat/entropy; and atomic refusal for missing observables,
minimum-step limits, trial/residual budgets and late cancellation. A local
step-doubling indicator is not a certified global time error or spatial error
bound. The existing fixed-chain geometry and material-domain restrictions remain.
