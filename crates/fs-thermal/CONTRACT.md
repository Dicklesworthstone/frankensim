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
