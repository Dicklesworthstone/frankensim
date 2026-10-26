# Live liquid material feedback on the enthalpy timeline

The source-bound `LiquidFlowEnthalpyInput` now has a `feedback` entry for a
temperature-dependent liquid viscosity. It uses the same pinned
`LiquidTemperatureLaw`, sourced heat-capacity chart, pressure/composition
context and actual borrowed caloric law as the existing frozen-flow binding.

```rust,ignore
// binding: LiquidFlowEnthalpyInput, already naming the original liquid/caloric
// sources, actual voxel-to-thermal cells and incoming reservoir enthalpies.
let observer = binding.feedback(
    &accepted_thermal, &voxel_domain, &simple_config,
    existing_thermal_observer, &cancel_gate, max_momentum_solves,
)?;
let (next_thermal, next_observer, energy) = accepted_thermal.advance_observed(
    &observer, target_time_s, thermal_config, checkpoint,
)?;
```

Every actual thermal substep follows one physical sequence: obtain T from that
cell's current h and original caloric law; resolve the pinned liquid material;
assemble its local molecular viscosity in the [existing full-stress momentum
operator](MATERIAL_FLOW_FEEDBACK.md); pass the actual conservative face fluxes
to the original `Advected` owner. That owner supplies the same-step mass and h
fluxes, boundary-energy account and combined conductive/advective CFL. The
wrapped source/contact/history observer remains on those same substeps. This
adds no alternate heat integrator or fabricated Cp-based temperature update.

`applied()` retains the last solved flow and its exact starting epoch and source
samples. `endpoint()` retains material evaluations at the accepted thermal
endpoint. `audit()` sums the original advective accounting over the advance.
Boundary enthalpy is already included in the thermal audit and must not be
added a second time. `inner()` preserves the original physical history.

This source-bound entry admits fully liquid, validity-wide constant density.
Viscosity can vary locally and between substeps, and source-domain failures
refuse instead of switching claims or extrapolating. The existing Cp/rho claim
pins, complete non-temperature context and actual caloric-object identity are
still checked. Both the initial and proposed endpoint must be admitted. Any
material, momentum, budget, history or cancellation failure discards the whole
tentative thermal/flow/history candidate.

Momentum is quasisteady: its equilibration time must be small relative to the
physical thermal evolution time. The returned flow belongs to the explicitly
reported starting temperature, not an implicit endpoint solve. Time and space
refinement remain necessary. Density-change continuity, phase mobility, moving
surfaces, inertia history, pressure-work heat and viscous dissipation are not
supplied. Unassigned thermal/contact links retain their original owner.
Existing full fluid-face links can opt into the same-sample conduction below.

`tests/liquid_flow_enthalpy/feedback.rs` exercises data-only viscosity curves,
native Darcy flow and donor-cell enthalpy, fixed-viscosity parity, pinned source
support, borrowed law/context consistency and atomic refusal. The independent
Darcy-series energy recurrence is a reference check, not native Rust test
execution or experimental material qualification.

## Molecular conduction from the same material epoch

Call `observer.with_fluid_conduction(&accepted_thermal, &fluid_face_link_ids)?`
to assign existing full fluid-face thermal links. Each actual substep then
uses the **same resolved transport sample** for molecular viscosity and thermal
conductivity. Changing conductivity data changes the next temperature field,
which changes local viscosity and the next momentum solve.

For neighboring cells with normal widths `dx_a, dx_b`, the applied link law is

```
G = A / (dx_a / (2 k_a) + dx_b / (2 k_b))
```

The shared `fs_thermal::enthalpy::series_conductance_w_k` also supplies the
existing material circuit network. The native graded voxel geometry supplies
area and both half-cell distances; no material query, enthalpy law or interface
resistance is added. The original thermal link owner applies opposite heat
rates and uses these same conductances in its combined conductive/advective
timestep bound. Authored links and both covered areas remain unchanged.

This opt-in assignment requires an existing link across adjacent mapped fluid
voxels, their native total surface areas, and full face coverage on **both**
sides of the authoritative contact patch. Partial/solid interfaces, changed
endpoints, duplicate or parallel links on an owned face, and a wrapped observer
supplying that same link's conductance refuse. Unassigned contacts retain their
original owner. Source time and every applied cell enthalpy must still match
the actual thermal state before the stored transport bundle can supply k.
The assignment is bounded by the existing connection budget.

The added native tests exercise a conductivity-only card swap affecting heat
spread, outlet energy, subsequent viscosity-driven flow and the actual CFL
substeps; a graded-cell series-resistance/mass balance; and atomic refusal for
partial coverage, parallel/competing ownership and stale source epochs. These
remain numerical coupling fixtures, not experimentally qualified material laws.

## Conducting liquid with a moving free surface

`MaterialFreeSurfaceEnthalpy3::step_with_conduction`, enabled by the existing
`fs-lbm/enthalpy-coupling` feature, now connects source `k(T)` to the native
free-surface momentum and extensive-enthalpy transport. The original
`step`/`step_with_heat` path remains available. This entry first advances the
actual positive VOF cell masses and the original same-card caloric chart on
`fs_thermal::enthalpy::DepositionSimulation`. The existing thermal owner applies
opposite face heat, chooses its explicit CFL substeps and audits storage. Each
substep re-evaluates the original pinned liquid claims at its actual temperature.
The accepted thermal field then supplies local viscosity to the original native
collision and mass-carried enthalpy step. Both stages cover one flow interval.

Each `FreeSurfaceConductionFace` supplies two adjacent occupied cells, their
actual wetted contact area and both liquid-centroid distances to that face.
Conductance uses the shared series-resistance law and both source conductivities.
The area must fit the Cartesian face; duplicate faces, nonpositive distances
and endpoints without positive native liquid mass refuse. Undeclared faces are
explicitly insulated. No averaged fill fraction invents a wetted contact, gas
heat transfer or ambient surface. Supply updated physical geometry on each call.

The report distinguishes the last applied source samples from the thermal
endpoint before VOF advection. Internal conduction is absent from external heat;
explicit `FreeSurfaceCellHeat` inputs are delivered uniformly over this interval
and counted once. The combined audit retains the original pre-thermal enthalpy.
No heat, flow, conversion, source history or transport state commits unless both
stages and their endpoint sources pass, including work and cancellation checks.

This first-order thermal/flow split needs time refinement for nonlinear feedback.
It retains constant density, one admitted liquid phase and neglected capillarity.
Phase-volume change, surface energy and pressure work remain separate physical
owners. The optional viscous-heating composition below supplies molecular heat. The focused `material_free_surface_enthalpy` integration
target adds an independent two-cell nonlinear-conductivity balance, a
conductivity-only card swap affecting actual native momentum, external-source
accounting and atomic failure/retry regressions under `enthalpy-coupling`.

## Source-density buoyancy in a moving liquid

`MaterialFreeSurface3::with_boussinesq_enthalpy` connects the same pinned liquid
density law to the original sparse free-surface momentum solver. The existing
`BoussinesqMaterial` owner bounds source density and heat-capacity departures.
The original VOF reference mass/volume and caloric h chart remain fixed; actual
source `mu(T)/rho0` supplies native relaxation and `k(T)` supplies conduction.
The unprojected material samples and their receipts remain observable. No
ideal-gas expansion coefficient is inferred for a liquid.

The supplied gravity must equal the existing domain gravity. Its reference
load `rho0*g` stays present at the free surface, with `(rho(T)-rho0)*g` added as
a spatial perturbation. Native collision, pressure-boundary reconstruction,
new-interface initialization and velocity sampling carry that actual force.
The existing uniform-force boundary approximation is translated through the
perturbation source frame, including both anti-bounce-back populations, so its
zero-density-departure limit is continuous and the exact limit retains the
original path. Newly created interface populations use the same full Guo
half-source transform as force rebinding, including its second moment.

Attachment preserves physical velocity and force-corrected stress even for
an evolved flow; raw populations change only to express the new force frame.
Each step holds the post-thermal source force through native momentum. The
accepted advected h then refreshes source support and rebinds the endpoint
force frame without adding an extra half-force impulse. The same-force fast
path preserves population bits. Thermal/conduction, native force, VOF routes,
conversion counters and all enthalpy/source history are accepted together.

This is a first-order held-force thermal/flow split with a bounded Boussinesq
momentum approximation. It does not supply a thermodynamic volume source,
phase conversion, capillary energy, pressure work, gravitational total-energy
closure or physical viscous heating. Existing native low-Mach and boundary
accuracy limitations remain. `material_free_surface_enthalpy` adds six focused
regressions: exact reference-density behavior; independent held-force impulse
through changed heat; source-only density replacement with and without gravity;
a partial-fill density-refinement limit under gravity; transactional source,
cancellation and late-flow failures; and common conduction/buoyancy temperature.
Compilation and source review are not execution or experimental validation.

## Molecular viscous work heats the moving liquid

`MaterialFreeSurfaceEnthalpy3::step_with_viscous_heating` converts the
accepted molecular stress power into finite liquid enthalpy before the same
native momentum and VOF transport step. `step_with_conduction_and_viscous_heating`
also carries that power through the existing thermal substeps, alongside
source-dependent conductivity and authored heat. Each entry advances flow
once and composes with the existing source-density buoyancy binding.

The fixed-domain thermal-flow and free-surface consumers share the native
D3Q19 force-corrected non-equilibrium strain estimator. The physical rate is
`sqrt(2 D:D)`, and the original liquid source supplies dynamic viscosity.
Applied molecular power is `mu * rate^2`; each located cell receives that
power times its reference liquid volume and the physical timestep. Reference
volume is its actual VOF carrier mass divided by retained reference density.
Zero-mass interface ghosts and numerical conversion carry have no resolved
cell stress and receive no invented viscous heat. Reframing the Guo force
must preserve this physical stress, including when buoyancy changes.

`FreeSurfaceEnthalpyReport` keeps `external_heat_j` and
`viscous_heat_j` separate. Its enclosing balance subtracts both from the
extensive enthalpy change. The conduction subsystem's native thermal audit
sees their sum crossing its boundary; internal conduction is still internal.
Returned per-cell samples expose the original material query, strain rate,
reference volume, power and heat actually used.

This is an explicit first-order source split: stress and molecular viscosity
are sampled at the accepted initial state, heat is applied, then the original
owner refreshes collision viscosity and optional buoyancy before flow.
Strong feedback requires timestep refinement. No kinetic-energy discrepancy,
numerical damping, bulk-viscosity loss, pressure work or surface energy is
converted into heat. The existing single-liquid, stationary-wall,
reference-density and neglected-capillarity limits remain. This closes the
constitutive heat source and extensive enthalpy account; it does not claim
exact discrete conservation of combined kinetic and thermal energy.
