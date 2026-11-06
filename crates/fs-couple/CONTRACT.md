

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
heat into both bodies and a solved common interface temperature. No prescribed
heat fraction or arithmetic temperature mean is inserted. The source callback
owns original conductivity claims and other transport, retaining its query and
selection evidence in the returned payload. The route excludes its old direct
link to avoid duplicate conduction. Only actual native Coulomb dissipation
enters this heat source. Normal plastic work and signed temperature/recession
parameter work remain separate; their explicit external-work budget is the sum
of absolute parameter contributions. Native contact identities and the thermal
cell, interface and aggregate balances have independent acceptance tolerances.

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
validity. Removed material has no modeled mass, enthalpy or debris export.
Reciprocal thermoelastic free energy, conversion of normal plastic work to heat,
moving contact coverage and storage in caller-supplied tangential body response
are not closed by this owner. Its local contact identities and fixed-inventory
heat balance imply no complete body-energy, entropy or passivity theorem.

The candidate seals successor source, contact history and thermal state for
joint adoption. Bounded trials, Newton iterations and thermal inventory,
source-domain or footprint refusal, excessive cumulative removal or parameter
work, nonconvergence and cancellation publish no partial physical history.
The report retains both native mechanical operations, source callback evidence,
friction heat, normal/plastic/parameter accounts and cumulative removed volumes.
