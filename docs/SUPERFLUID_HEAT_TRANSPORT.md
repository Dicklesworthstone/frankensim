# Superfluid heat transport, viscous walls and vortex friction

The `superfluid` feature connects the finite-temperature Landau two-fluid
material, the original two-sound channel, native rigid pistons, and the existing
finite thermal timeline. A heater can now launch counterflow and second sound
while its own finite caloric storage changes the heat it supplies. Two finite
contacts can connect a heater and a thermal sink through the same channel,
including authored normal-fluid wall drag or unsteady transverse-profile memory.
Distributed irreversible heat can enter explicit native caloric receivers on
the same timeline as both finite contacts.

## Normal-fluid wall viscosity and fourth sound

`NormalViscousWall` supplies an explicit permeability per fluid volume. The
stationary-wall force is `F_n = -eta_n*v_n/k`; only the normal component feels
this viscous force. The superfluid retains slip. The normal inertia is
`rho_n`, and the walls receive an external impulse opposite to the fluid's.

```rust,ignore
use fs_material::superfluid::NormalViscousWall;
use fs_duct::superfluid::TwoFluidChannel;

let wall = NormalViscousWall::circular_capillary(radius_m)?; // k = R^2 / 8
let channel = TwoFluidChannel::try_new_with_normal_wall(
    material, length_m, area_m2, cells, limits, wall,
)?;
```

`parallel_plates(full_gap_m)` supplies `k = gap^2/12`; `try_new(k)` accepts
authored permeability. The geometry helpers are the steady laminar Poiseuille
resistances. Their use requires a continuum no-slip normal fluid, a developed
profile, small normal Reynolds number and low Womersley number. A material
name or channel area alone does not establish those conditions. Porous-medium
permeability must be converted to the stated fluid-volume basis; porosity and
tortuosity are not inferred.

The channel integrates normal-wall relaxation together with its midpoint
pressure/temperature gradients and body forces. Cell balances use the exact
interval-mean velocities. This avoids a force-splitting error in the stiff
equilibrium `v_n = F_n*k/eta_n`, and lets a narrow passage approach the
clamped-normal fourth-sound limit without an explicit viscous stability step.
The existing wave Courant policy still governs temporal resolution. With zero
viscosity, the original first/second-sound numerical path is retained exactly.

The dissipative operator remains two scalar solves by the original native
tridiagonal owner, preceded by a two-by-two rotation. Its work and storage are
O(cells); there is no new dense solver, per-step heap allocation or per-face
transcendental evaluation. Both mechanical and thermal endpoint responses
use the same viscous operator and retain reciprocity.

The independent clamped-normal limit is
`c4^2 = f_s * [c_ad^2*(1-alpha*s*T0/cp)^2 + s^2*T0/cp]`.
At zero thermal expansion this becomes `f_s*c1^2 + f_n*c2^2`. Tests retain the
corresponding frozen entropy-per-volume polarization at finite expansion.
Fourth sound emerges from the original two-fluid equations; a replacement
sound speed is not fitted to the channel.

`ChannelStep.normal_viscous_heat_j` and `interior_wall_impulse_n_s` keep wall
heat and support impulse distinct from mutual friction and authored forces.
`dissipated_heat_j()` sums both heat mechanisms. Existing mutual-friction-only
accessors preserve their meaning; the channel additionally exposes normal-wall
and combined heat histories. The native closed-channel and piston thermal
compositions receive both losses through their existing mapped caloric cells,
once per accepted step. Refusal or cancellation publishes neither loss.
End-half-cell inertia and viscous forces remain the boundary owner's separate
responsibility; the channel's distributed drag covers resolved interior dual
volumes. The current piston composition includes their inertia and reports the
interior wall reaction; it does not infer side-wall drag on those endpoint
half volumes.

This scalar Darcy model does not resolve transverse startup or high-frequency
Stokes boundary layers. The spectrum below supplies that unsteady response.
Neither wall closure resolves turbulence, longitudinal or second viscosities. The
frozen two-fluid properties and small-signal thermal limits remain explicit.
The wall-force convention and its frequency-dependent limitations follow
[Kekutia and Chkhaidze, section 3](https://arxiv.org/abs/cond-mat/0604193).

## Unsteady transverse normal-fluid profiles

`NormalViscousSpectrum` retains the viscous memory of a circular capillary or
wide parallel-plate slit. It is bound to the same frozen `TwoFluidState` as the
channel, including normal density and viscosity. Each mode obeys
`dU_m/dt + lambda_m*U_m = a_n`; its positive weight determines the mean flow,
kinetic storage and irreversible heat:

| Quantity per fluid volume | Modal expression |
| --- | --- |
| Section-mean normal velocity | `sum(w_m*U_m)` |
| Normal kinetic energy | `rho_n/2 * sum(w_m*U_m^2)` |
| Normal viscous heating rate | `rho_n * sum(w_m*lambda_m*U_m^2)` |

Circular modes use J0 zeros with `w_m=4/j_m^2` and
`lambda_m=(eta_n/rho_n)*j_m^2/R^2`. Slit modes use the full gap H,
weights `8/((2m+1)*pi)^2` and rates `(eta_n/rho_n)*((2m+1)*pi/H)^2`.
These are classical transverse Stokes eigenfunctions, not an empirical
frequency fit. See [Szymanski's tube solution](https://www.numdam.org/item/JMPA_1932_9_11__67_0.pdf),
[Womersley's oscillatory-flow analysis](https://doi.org/10.1113/jphysiol.1955.sp005276),
and [NIST's Bessel-root conventions](https://dlmf.nist.gov/10.21).

```rust,ignore
use fs_material::superfluid::NormalViscousSpectrum;

let spectrum = NormalViscousSpectrum::circular_capillary(material, radius_m, 8)?;
let mut channel = TwoFluidChannel::try_new_with_normal_profile(
    material, length_m, area_m2, cells, limits, spectrum,
)?;
channel.initialize(&pressure, &temperature, &normal, &superfluid, || false)?;
// Optional: interpret the authored mean flow as already fully developed.
channel.initialize_normal_profile_steady(|| false)?;
```

Choose 1 through 32 resolved modes. One additional positive mode represents
the unresolved infinite tail by matching its inertia and static mobility.
The complete approximation retains unit normal inertia and the Poiseuille
resistance. It does not discard or renormalize omitted inertia.
`harmonic_mobility_s(omega)` returns velocity/acceleration for the
`exp(+i*omega*t)` convention. `harmonic_mobility_error_bound_s(omega)` bounds
only the analytical error from the aggregate tail; it does not cover
floating-point, longitudinal discretization, timestep or physical-model errors.
Increase the mode count to meet the frequency band's accuracy needs.

Ordinary initialization gives each mode the authored mean velocity, a uniform
initial profile with zero excess energy. The optional steady initializer seeds
inverse-rate modal velocities while preserving the mean. It admits the
profile's positive kinetic energy before time advances. Zero viscosity has
inviscid modal dynamics but no selected steady profile. The public modal
values are spectral amplitudes, not point samples across the capillary.
Admission bounds the profile RMS speed as well as the original section means.

The channel stores accepted and candidate profile memory at interior faces.
`normal_profile_excess_energy_j()` reports storage above section-mean kinetic
energy; this is already included in `wave_energy_j()` and every step report.
It must not be added again by a coupler. Existing closed thermal and native
rigid-piston compositions therefore retain the complete storage. Exact positive
heat integrals and wall impulses use the material owner's shared factors.
Every accepted loss is exported once through the original per-cell heat arrays.
Rejected response probes, caloric intake or cancellation publish no profile
memory, clock, heat history or field.

The longitudinal solve remains two original scalar tridiagonals. Profile work
and memory are O(cells * modes), with at most 33 modes. Exponentials are cached
once per distinct timestep outside the face loops; the channel step performs
no heap allocation. Both endpoint response matrices use the same operator.
The wall integration has no explicit viscous stability restriction; the wave
Courant policy and timestep refinement still govern temporal accuracy.

This resolves fully developed linear transverse momentum diffusion, including
startup and oscillatory phase lag, without a low-Womersley assumption. It
requires continuum no-slip normal flow, fixed geometry and frozen properties.
The section-mean mutual-friction closure accelerates all normal modes equally;
it does not resolve local radial vortex forces. Entrance regions, moving side
walls, normal turbulence, vortex evolution and transverse thermodynamics remain
outside the model. Endpoint half volumes retain their separate existing owner;
interior profile memory does not infer an additional end-wall force.

## Run the viscous heater-to-sink example

```bash
cargo run -p fs-couple --features superfluid --example superfluid_heat_link
```

This connects two finite native caloric cells through the sourced **1.8 K SVP
helium-II** reference in an authored one-micrometre-radius, one-metre capillary.
Eight resolved circular Stokes modes and one matched tail retain the normal
profile's memory. The example starts from an independently derived steady
counterflow/fountain-pressure profile and seeds its developed modal shape,
then advances both finite temperatures and both contacts together.
Its geometry, caloric capacities and conductances are demonstration inputs,
not measured apparatus calibration. It prints normal relaxation rate, both
temperature histories, heat power at each end, each endpoint's caloric balance,
distributed viscous heat, profile excess energy and wave availability.
One explicitly authored finite caloric receiver per channel cell collects the
local irreversible heat through the composition below. The example checks each
receiver's intake against the original per-cell loss and the summed stored
heat against cumulative channel dissipation. These receiver capacities are
demonstration inputs and do not update the helium's frozen coefficients.

## Run the finite heater example

```bash
cargo run -p fs-couple --features superfluid --example superfluid_heater
```

This uses the existing sourced helium-II point at **1.8 K and saturated vapor
pressure**. It prints heater temperature, channel temperature at two locations,
heat flow, and independent reference-heat balances. The heater capacitance,
interface conductance, geometry, source power and frozen vortex population are
explicit demonstration inputs, not fitted measurements of an apparatus.

The example propagates both sound modes, retains thermal expansion, and applies
isotropic vortex friction on the original channel. Exported friction heat and
wave availability are printed separately from heat transported by counterflow.

## Finite caloric contact

`fs_couple::superfluid::counterflow::advance_counterflow_heat` connects a selected
channel end to an active `fs_thermal::enthalpy::CellId`. The opposite end is
closed. The contact has an authored **total conductance in W/K**, including area.
Both endpoints carry zero net mass current; heat flow is positive from the
channel into the finite thermal cell.

For each interval it:

1. Prepares the original channel's reciprocal midpoint temperature response.
2. Advances the original thermal sources, links and scheduled events once.
3. Solves the interface law `Q = G * (T_channel_mid - T_cell_mid)` against the
   original cell's caloric chart, with bounded scalar probes. The channel's
   wave impedance and the finite cell's heat capacity both feed back into Q.
4. Applies exactly `Q * dt` through native `apply_heat_j`, preserving original
   material domains and any latent plateau.
5. Publishes both histories only after physical limits, transfer checks and
   cancellation admission succeed.

The native thermal network followed by the contact is a **first-order split**;
refine the timestep for coupling accuracy. The contact itself uses midpoint
temperatures. Probes do not repeatedly solve the wave equation or clone the
thermal network: they query its original caloric chart at candidate energies.
An out-of-chart trial is shortened while bracketing; it never extrapolates.

The source and sink must remain in the caller's declared small-signal range
about the channel's frozen reference temperature. No nonlinear helium equation
of state, phase crossing, evolving material coefficients, or uncalibrated
Kapitza conductance is inferred.

## Two finite thermal contacts

`fs_couple::superfluid::counterflow_pair::advance_counterflow_heat_pair` accepts
two `CounterflowThermalContact` values, one for each channel end. Contacts may
be supplied in either order; reports use `[left, right]`. The cells must be
distinct and belong to the same original thermal simulation, whose complete
network and sources are advanced once.

The two heat powers satisfy the contact laws together with the channel's full
reciprocal 2x2 response. Each finite cell supplies its original monotone caloric
chart, including nonlinear heat capacity or latent plateaux. The bounded
coordinate solve makes only caloric probes; it does not rerun the wave or
thermal network at each trial. Intermediate power or chart projections do not
prematurely reject an admissible joint solution. Both final contact residuals
must pass before one combined native heat intake and the physical wave step
publish atomically.

The existing `CounterflowHeatControls` apply to both contacts, with
`max_iterations` bounding probes per endpoint. Each conductance may be zero;
one active contact reduces to the single-contact equations. As in the original
single-contact composition, thermal-network coupling is a first-order split
and reference-state heat is distinct from quadratic wave availability.

## Distributed heating during finite-contact transport

`counterflow_pair::dissipation::advance_counterflow_heat_pair_with_dissipation`
adds an explicit `destinations: &[CellId]` to the paired-contact composition.
Exactly one active, unique receiving cell is required per channel cell. Both
endpoint contact cells must be distinct from those receivers: adding heat to a
contact after solving its midpoint temperature would invalidate the contact
equation. A route is never inferred from geometry or a material name.

The original native network advances once. Both finite-contact powers solve
against its original caloric laws; their combined heat intake is prepared.
The actual wave step supplies its current per-cell mutual-friction plus
normal-viscous heat to a second native caloric intake for distinct distributed
receivers. Only after all laws, balances and cancellation checks pass do the
wave, profile and complete thermal candidate publish together. This uses one
physical wave step and one native network march. It does not repeat damping,
estimate loss from kinetic-energy subtraction or replay previously accepted heat.

`CounterflowDissipationStep.contacts` preserves the original reference-heat,
entropy, wave-port and contact-intake accounts. Its separate
`dissipation_intake` and `dissipation_transfer_residual_j` account for the
irreversible distributed heat. Each receiver's original caloric chart remains
authoritative, including phase plateaux. The thermal-network split remains
first order and the channel coefficients remain frozen. Choose and justify
the physical heat route for the intended apparatus.

## Three distinct thermal quantities

For reference temperature `T0` and channel temperature perturbation `theta`:

| Quantity | Outward rate | Meaning |
| --- | --- | --- |
| Reference-state counterflow heat | `Q` [W] | Heat exchanged with the caloric owner |
| Reference entropy flow | `Q / T0` [W/K] | Linearized entropy transport |
| Wave-availability power | `theta * Q / T0` [W] | Boundary power in the quadratic wave model |

These rates are not interchangeable. The report separately retains native
caloric energy, reference heat and entropy residuals, the wave port residual,
and the contact's linearized availability loss. That availability loss is not
extra heat to deposit. The coupled reference-heat account is not a nonlinear
total-internal-energy theorem for helium or an arbitrary enthalpy-based solid.

Optional mutual-friction and normal-wall losses remain in the channel's
per-cell heat fields and report. This boundary contact does not redirect all
distributed heating into the heater. Existing closed-channel and piston thermal
compositions remain available for explicitly mapped dissipative-heat destinations.
The paired distributed-heating composition above provides such destinations
while both finite boundary contacts operate.

## Thermal boundary response for other couplings

`TwoFluidChannel::counterflow_boundary_velocities([Q_left, Q_right])` converts
authored total watts into opposite normal and superfluid velocities.
`counterflow_ends_response` returns the full reciprocal 2x2 midpoint-temperature
response, retaining propagation between both ends, both sound modes, body forces,
thermal expansion and the original friction law. It reuses the cached native
tridiagonal factors and endpoint columns in O(cells) work and storage.

The response is an algebraic proposal. It never publishes a field, clock or heat
history, and does not replace the final step's physical amplitude admission.

## Vortex-dependent mutual friction

`fs_material::superfluid::VortexMutualFriction` requires explicit circulation
magnitude, vortex-line density, dissipative Hall–Vinen B and reactive B-prime.
`TwoFluidState::relax_hall_vinen` analytically advances homogeneous flow through
a frozen aligned straight-vortex population. It leaves parallel counterflow
unchanged, damps and rotates transverse counterflow, preserves barycentric
momentum, and exports heat from dissipative transverse motion only. Pure
reactive rotation produces no heat.

The scalar channel supports `MutualFriction::IsotropicVortexTangle`. Its force
average is `gamma = (rho_n*rho_s/rho) * B*kappa*L/3`; it uses the existing exact
linear-drag step and heat destinations. A one-dimensional channel does not
resolve oriented transverse vortex motion. Neither API evolves or nucleates
vortices, models tension/reconnections, or predicts a critical velocity.

The coefficient conventions and angular average follow
[Jou, Mongiovi and Sciacca, equations 2.5–2.9](https://arxiv.org/abs/0806.1889).
The helium reference and B/alpha convention follow
[Donnelly and Barenghi (1998), section 11](https://doi.org/10.1063/1.556028).
The latter source documents frequency and velocity dependence; the API therefore
does not silently assign a universal B value to a material name.

## Focused regression targets

```bash
cargo test -p fs-material --test superfluid
cargo test -p fs-material --test superfluid_viscosity
cargo test -p fs-material --test superfluid_normal_spectrum
cargo test -p fs-duct --features superfluid --test superfluid_channel --test superfluid_ports
cargo test -p fs-duct --features superfluid --test superfluid_viscosity
cargo test -p fs-duct --features superfluid --test superfluid_normal_profile
cargo test -p fs-couple --features superfluid --test superfluid_counterflow --test superfluid_heat
cargo test -p fs-couple --features superfluid --test superfluid_counterflow_pair
cargo test -p fs-couple --features superfluid --test superfluid_transport_dissipation
```

The contact tests include independent finite-capacity/two-cell equations,
heating and cooling at either end, nonlinear caloric feedback, original sources,
latent storage, reference heat and linear availability balances, zero-contact
identity, and atomic refusal/cancellation/retry. Vortex regressions cover
direction, frame covariance, opposite impulses, decay, reactive rotation,
weak/strong damping, the angular force average and channel equivalence.
Viscous regressions add exact normal relaxation, stiff force balance,
fourth-sound speed and refinement, zero-viscosity identity, reciprocal ports,
external wall momentum, both heat mechanisms on native caloric cells, native
piston reactions and all-owner refusal/cancellation.
Pair-contact regressions add independent coupled finite-capacity equations,
nonlinear heater-to-sink transport with viscous walls, simultaneous latent
plateaux, zero/one-contact limits, temporary power/chart projection, and
complete thermal/wave/viscous-history refusal and cancellation.

Profile regressions add independent circular/slit startup and harmonic
references, exact inertia and steady resistance, positive tail bounds,
developed-profile storage, full energy/impulse/heat balances, flow reversal,
native piston coupling and transactional mode memory. Distributed transport
checks finite contacts and local native calorics together, original source
intake, mapping and chart refusals, and early/final cancellation with accepted
nonzero profile history.
