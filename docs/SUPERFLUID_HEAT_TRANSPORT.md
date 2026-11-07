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

## Evolving vortex tangles with represented line energy

The optional `fs_material::superfluid::vinen::VinenTangle` advances an
unpolarized, unresolved local population instead of holding its density fixed:

```text
dL/dt = a |v_n-v_s| L^(3/2) - b L^2
d(v_n-v_s)/dt = -(B kappa/3) L (v_n-v_s)
```

Here L is line length per total volume [1/m²], a is dimensionless, and b has
units m²/s and already includes circulation. All coefficients are authored.
The homogeneous law requires a remanent seed: L=0 remains exactly zero. It
does not supply a critical velocity or an automatic nucleation law.

A constant effective line energy e [J/m] is REQUIRED. Represented storage is
`mu |v_n-v_s|²/2 + e L`, where `mu=rho_n*rho_s/rho`. Growth reserves its
line-energy increase from the lost relative kinetic energy; decay releases
line energy. Only the remaining nonnegative heat reaches a caloric owner.
Admission requires `e a² <= 4 mu (B kappa/3) b`, the condition making the
closure's heat-production quadratic nonnegative at every state. Incompatible
coefficients refuse rather than clipping production or inventing energy.

The numerical source uses the exposure `tau=integral(L dt)`. In that
coordinate the relative speed and sqrt(L) satisfy a constant triangular linear
system, advanced by a positive exponential map. Physical-clock integration is
second order, with explicit fractional-rate and substep budgets. Production-
free decay uses its exact physical-time solution, including drag integrated
over the changing L. Weak impulses, coincident decay rates and representable
strongly attenuated tails are retained. Source updates allocate no memory.

Configure one line density on each INTERIOR face of an existing channel before
the first step. This example uses demonstration coefficients, not measured
vortex-growth or line-energy calibration:

```rust
use fs_material::superfluid::vinen::{
    VinenParameters, VinenStepControls, VinenTangle,
};
use fs_duct::superfluid::MutualFriction;

let tangle = VinenTangle::try_new(
    *channel.material(),
    VinenParameters {
        circulation_quantum_m2_s: 9.97e-8,
        dissipative_b: 1.0,
        production_coefficient: 0.05,
        decay_coefficient_m2_s: 1e-8,
        line_energy_j_m: 1e-12,
    },
    VinenStepControls {
        maximum_fractional_step: 0.05,
        maximum_substeps: 512,
    },
)?;
let seeds = vec![1e6; channel.cells() - 1];
channel.configure_vinen_tangle(tangle, &seeds, || false)?;
let friction = MutualFriction::EvolvingVinenTangle;
```

Pass this selector to the existing channel or `fs-couple` stepping operation.
Once configured, other friction selectors refuse: they would silently freeze
the tangle or apply another drag law. Density is accepted/pending state alongside
the original wave and transverse normal-profile memory. Endpoint response
queries never publish it. A source-budget refusal, failed heat acceptance or
cancellation leaves the accepted density, velocity, heat and clock unchanged.

`vortex_line_density_per_m2(face)` observes interior faces 1..N-1.
`vortex_energy_j()` returns their line store; `represented_energy_j()`
adds it to the original `wave_energy_j()`. Reports keep wave-only fields
and expose separate `vortex_energy_before_j` / `vortex_energy_after_j`.
When configured, `wave_balance_defect_j` includes the line-storage change.
The original thermal callbacks receive heat after line growth or decay has
been accounted for, through the existing per-cell route.

The local source is a frozen-reference, constant-effective-line-energy closure.
It does not resolve line geometry, tension, reconnections, polarization, radial
variations in vortex drag, or the core/spacing dependence of a vortex equation
of state. Optional line transport is described below. Heat deposition does not relinearize helium properties or
supply phase-transition feedback. The represented balance is not a theorem for
the full nonlinear internal energy of helium.

The Vinen convention is equation (1.3) of
[Sciacca, Mongiovi and Jou](https://cab.unime.it/mus/501/1/C1S0801017.pdf).
The need to couple evolving lines to thermodynamics is discussed by
[Sciacca](https://arxiv.org/abs/1008.4390); the restricted effective-energy
closure and coefficient inequality above specify this implementation.

## Stiff linear mutual friction and driven second sound

With no active normal-wall or configured transverse-profile operator, positive
`MutualFriction::Linear` and its frozen `IsotropicVortexTangle` equivalent now
share the original exponentially fitted wave solve. A uniform counterflow
driven by opposite component forces retains `gamma*w = F_n` even when the
friction time is much shorter than the wave step. Its endpoint velocity,
interval-mean transport, body-force work and positive friction heat come from
the same exact held-gradient evolution. The caller's wave Courant and amplitude
limits remain in force; there is no additional friction substep loop or spatial
allocation. Factors and endpoint responses are invalidated when the friction
coefficient changes, including changes at the same timestep.

This resolves the strong-drag failure of separately applying friction half
steps around an undamped force/wave step. For the homogeneous forced problem
that split has the stationary-flow ratio `x/(2*sinh(x/2))`, where
`x=gamma*dt/(rho*f_s*f_n)`, which tends to zero instead of one as drag stiffens.
The fitted path retains the physical stationary solution and the correct
transient relaxation and dissipated work.

Nonlinear Gorter-Mellink/evolving-Vinen friction and friction combined with
wall/profile dynamics still use the original symmetric split. Their competing
source and wave errors require timestep refinement; the exact linear-path
claim does not extend to those mixed or nonlinear cases. The existing heat,
counterflow and piston compositions automatically use the corrected native
path whenever its conditions apply.


## Connected channels and ideal superleaks

`fs_duct::superfluid::TwoFluidConnection` joins the right end of one native
channel to the left end of another. Both channels retain their original wave,
wall/profile, body-force and vortex/friction owners. They must share exactly the
same immutable reference material and accepted time; their lengths, areas and
cell counts may differ.

The connection stores the two endpoint half volumes' inertia. Simply equating
adjacent-cell pressures and temperatures would omit that physical storage and
make a split channel differ from the original unsplit discretization. With
`D = dx_left/(2*A_left) + dx_right/(2*A_right)`, the connection stores

```text
K_junction = rho*D*U^2/2 + rho*f_s*f_n*D*W^2/2,
U = A*v,  W = A*(v_n-v_s).
```

Flows point from the first channel to the second. A small inertia-normalized
reciprocal solve couples both native endpoint responses to this storage.
Authored body forces act on the joined half volumes as well as channel
interiors, and their work remains an external source. No dense channel solve,
finite-difference response probe or new numerical dependency is introduced.

`ConnectionKind::Open` conserves both normal and superfluid volume currents.
`ConnectionKind::Superleak` constrains normal flow to zero while retaining
superfluid inertia. Writing superfluid volume flow as `Q_s`, its constraint is
`(U,W)=(f_s*Q_s,-Q_s)`; its effective inertia is `rho_s*D`. The driving effort
is proportional to `delta p-rho*s*delta T`. It therefore preserves the
small-signal fountain equilibrium `p_hot-p_cold = rho*s*(T_hot-T_cold)`, and an
initially equal-pressure pair drives superfluid toward the warmer side.
Admission rejects an existing nonzero normal current instead of silently
removing its kinetic energy.

Run the sourced 1.8 K helium comparison with:

```sh
cargo run -p fs-duct --features superfluid --example superfluid_superleak
```

Create the connection with `TwoFluidConnection::try_new`, then call `step`
using `ConnectionDrive` for the two outer boundary velocities, per-channel
forces and friction selectors. `ConnectionStep` reports junction energy and
work, combined represented-energy residual, mass and entropy residuals, and
both unchanged native channel reports. The connection's `flow_m3_s()` is the
instantaneous accepted flow; joined channel endpoint faces contain the last
interval's mean flow, as required by the original boundary convention. Both
instantaneous component speeds and the native staged channels must satisfy
their original physical limits.

### Full thermomechanical endpoint response

`thermomechanical_ends_response` exposes the complete reciprocal 4-by-4
response in end-major order:
`[left bulk, left relative, right bulk, right relative]`. At each end the
outward volume flows are `(U,W)` and the conjugate efforts are
`(p',rho_s*s*T')`. Their product is wave-availability power. It retains
thermal expansion, communication between ends and the same fitted or split
dissipation operator used by actual steps. Cached endpoint columns come from
the native scalar tridiagonal factors.

This general response does not identify wave work with transported heat or
entropy. Full outward entropy current is `rho*s*(U+f_s*W)`; relative
counterflow heat is `rho_s*s*T0*W`. An ideal superleak can transport mass with
zero normal entropy current. The response query publishes no accepted state;
the final physical proposal independently checks every original limit.

### Connection scope

The joined half volumes use an inviscid plug-flow closure. Interior wall,
profile and mutual-friction models remain active, but no extra half-volume
drag, entrance loss, porous critical velocity or junction vortex transport is
inferred. Lossless equal-grid splitting is compared against the original
unsplit solver; that exact segmentation claim does not extend to a missing
junction dissipation model. Area-change and superleak support reactions do not
carry a whole-system axial-momentum claim. Vortex reservoirs and authored
vortex drift remain separate ports. Material coefficients remain frozen: this
adds neither nonlinear advection nor finite-temperature or lambda-transition
feedback.

Both channels and the connection publish together. A refused physical limit,
heat callback or final cancellation preserves all three accepted histories.
The focused regressions cover split/unsplit waves, unequal-area component
currents, fountain equilibrium and flow direction, half-volume forcing,
profile-memory rollback, and the instantaneous junction speed gate.


## Connected-channel dissipation on the native thermal timeline

`fs_couple::superfluid::connection::advance_connected_friction_heat` advances
two connected channels, their finite-inertia junction and one original
`DepositionSimulation` as a single accepted interval. Supply one explicit heat
destination per channel cell, with distinct active destinations across both
routes. Both outer ends must be closed; the internal connection may be open or
an ideal superleak.

The first thermal candidate advances authored sources, contacts and events to
the common endpoint. The second route accepts its actual per-cell dissipative
joules at that same endpoint, so external thermal sources advance once. Neither
channel's cumulative heat is replayed. Junction kinetic energy, boundary wave
work and body-force work remain in the original mechanical account; only
dissipative heat crosses into the caloric receivers.

The report retains both native heat-intake/march reports and the connection
report. Its combined residual adds thermal and transfer residuals to the
connection account without counting either channel wave residual twice.
Rejection by the second caloric law, an energy gate or the final cancellation
poll preserves both channels, connection and thermal simulation. Retrying
therefore cannot deposit an earlier rejected proposal.

This remains an endpoint heat split with the original caloric law's storage
convention and frozen wave reference. It does not feed the receiving cells'
finite temperatures back into the linear two-fluid coefficients. The
`fs-couple` `superfluid_connection` target exercises repeated acceptance with
native sources, second-route caloric refusal, final cancellation and retry,
overlapping-route rejection and refusal of unhandled outer transport.

## Discrete sourced helium reference temperatures

`HeliumSvpReference::at_temperature_k` now admits exactly 1.4, 1.6, 1.8 and
2.0 K at saturated vapor pressure. `SUPPORTED_TEMPERATURES_K` lists the
available points, `raw()` retains the selected published row, and
`at_1_8_k()` retains the original constructor. The existing channel examples
continue to use their declared 1.8 K point.

Density, superfluid density, entropy, heat capacity, expansion, sound and
viscosity come from [Donnelly and Barenghi](https://srd.nist.gov/jpcrdreprint/1.556028.pdf).
The original conversion from saturation-path derivatives to the isobaric
two-fluid inputs remains in use. The independently fitted tables are not
exactly mutually consistent: the predicted second-sound discrepancy is about
2.8% at 1.4 K, below 0.8% at 1.6/1.8 K, and below 0.1% at 2.0 K. These
discrepancies are retained and tested; entropy or density is not retuned to
hide them. Predicted fourth sound agrees within 0.5% at all retained points.

Intermediate temperatures and extrapolation refuse. Selecting another
equilibrium reference does not create a temperature-dependent EOS or
conservatively reinitialize a running channel.

## Authored counterflow onset

The optional `VinenTangle::with_counterflow_onset_m_s(W0)` changes production
to `a * max(|v_n-v_s|-W0,0) * L^(3/2)`. Decay and mutual friction remain
active below the threshold. Zero preserves the original source arithmetic;
an empty tangle stays empty at every speed. W0 is an authored or measured
application parameter, not a geometry-derived critical velocity or a
spontaneous nucleation law.

At an externally maintained counterflow speed, the nonzero source equilibrium
is `sqrt(L)=(a/b)*max(W-W0,0)`; the helper
`imposed_counterflow_equilibrium_per_m2` evaluates it for calibration checks.
The freely decaying velocity/tangle system does not maintain that speed.
The offset form is motivated by equation (8) in
[Baggaley and Laizet](https://arxiv.org/html/1310.5890v1); a steady fit does not
uniquely determine a transient law, so the chosen transient closure is explicit.

The active source freezes a reduced production coefficient at a predicted
midpoint and uses the same positive, dissipative exponential map. Every frozen
coefficient lies in [0,a], so the existing global heat condition remains
sufficient. The smooth active regime is second order; onset crossing has a
continuous vector field with a derivative kink and requires time refinement.
Once the source is subcritical, its remaining decay and drag are integrated
exactly. The one original substep budget still bounds the whole source call.

## Conservative vortex-line transport

A configured tangle can also carry an independently authored drift field:

```rust
use fs_duct::superfluid::{
    VortexTransportBoundaries, VortexTransportControls,
};

let n = channel.cells();
let mut drift_m_s = vec![0.02; n];
drift_m_s[0] = 0.0;
drift_m_s[n - 1] = 0.0;
channel.configure_vortex_line_transport(
    &drift_m_s,
    VortexTransportBoundaries::CLOSED,
    VortexTransportControls {
        maximum_outgoing_courant: 0.4,
        maximum_substeps: 256,
    },
    || false,
)?;
```

There are N-1 interior dual cells and N transport interfaces. Their transport
domain runs from dx/2 to length-dx/2; endpoint half-cell vortex inventories
remain the boundary model's responsibility. Drift is NOT silently equated to
normal, superfluid or barycentric velocity. It can be reconfigured atomically
between accepted steps.

The law `L_t + (u_L L)_x = 0` uses limited linear upwinding and SSPRK2.
The explicit budget limits the sum of outgoing Courant numbers to at most
one half per substep. Each internal interface transfers the same line length
out of one cell and into its neighbour; positivity requires no clipping.
Preparation allocates reusable buffers, and stepping uses O(N) work per
bounded transport substep with no transport allocation.

Reflecting ends require zero endpoint drift. A
`VortexTransportBoundary::Reservoir { incoming_line_density_per_m2 }`
uses that density only on inflow; outflow uses the channel's adjacent density.
Reports expose signed inward `vortex_transport_boundary_line_m` and
`vortex_transport_boundary_energy_j`, a line-conservation defect, and the
transport substep count. Each energy port is the constant effective line energy
times its line-length port. Transport itself generates no heat.

The channel orders transport/2, local source/2, wave, source/2, transport/2.
Its represented energy balance includes the external line-energy ports.
Existing heat and piston compositions retain those ports in their ChannelStep
and compose the corrected residual; a mechanically closed channel may still
have an explicitly open vortex reservoir. Transport shares accepted/pending
publication with the original owners. Queries, cancellation, work-budget
refusals and failed caloric acceptance never publish transported density.

This supplies unresolved line advection, including external seeds arriving
through an authored reservoir. It does not supply filament stretching,
orientation, tension, reconnections, or a dynamical equation for the drift.

## Continuous calorics on an authored isobar

`fs_material::superfluid::caloric::IsobaricSuperfluidLaw` supplies a continuous
finite-temperature material family on ONE externally maintained pressure.
Its knots author Cp, superfluid fraction, thermodynamic adiabatic sound speed,
expansion and normal viscosity. Preparation integrates the piecewise linear
inputs consistently:

```text
h(T)   = h0 + integral Cp(T) dT
s(T)   = s0 + integral Cp(T)/T dT
rho(T) = rho0 exp(-integral alpha(T) dT)
```

Entropy and density are therefore derived from the same capacity/expansion
data used by the local response. The entropy datum is physical because it
enters the two-fluid sound equations; the enthalpy datum is explicitly supplied.
The law identity binds all constitutive data, pressure, domain and datums.

`state_at_temperature_k` and `state_at_specific_enthalpy` return a compatible
`TwoFluidState` and caloric storage. Enthalpy inversion starts with a scaled
analytic quadratic estimate and has a 64-step bracketed fallback. Ordinary
queries locate a segment in O(log knots), integrate analytically and allocate
nothing. Entropy exchange integrates crossed intervals and retains the
original exchanged enthalpy even below absolute-temperature resolution.

The law directly implements `fs_thermal::enthalpy::ThermalLaw`, declaring
enthalpy storage, its caloric identity, a global temperature-slope bound and
entropy increments. It can be supplied directly to the original thermal
timeline and existing superfluid heat/contact compositions. Native coupling
regressions exercise temperature-dependent heating and cooling, combined
vortex/thermal energy, and complete rollback at the finite caloric-domain limit.

This is an integrable authored isobar, not a calibrated helium EOS. The four
measured helium references above lie on different saturation pressures and
must not be relabelled as one isobaric curve. The API does not infer arbitrary
pressure response, extrapolate, cross the lambda transition, or admit a pure
normal/superfluid endpoint. Thermal geometry is reference heat-exchange
geometry under the declared pressure load, not a rigid-vessel constraint.
A running wave channel retains its selected linearization; caloric queries do
not silently remap its wave, kinetic or vortex-line storage.

## Focused regression targets

```bash
cargo test -p fs-material --test superfluid
cargo test -p fs-material --test superfluid_vinen --test superfluid_vinen_onset
cargo test -p fs-material --lib superfluid::caloric::tests
cargo test -p fs-material --test superfluid_viscosity
cargo test -p fs-material --test superfluid_normal_spectrum
cargo test -p fs-duct --features superfluid --test superfluid_channel --test superfluid_ports
cargo test -p fs-duct --features superfluid --test superfluid_vinen --test superfluid_vortex_transport
cargo test -p fs-duct --features superfluid --test superfluid_viscosity
cargo test -p fs-duct --features superfluid --test superfluid_normal_profile
cargo test -p fs-couple --features superfluid --test superfluid_counterflow --test superfluid_heat
cargo test -p fs-couple --features superfluid --test superfluid_counterflow_pair
cargo test -p fs-couple --features superfluid --test superfluid_transport_dissipation
cargo test -p fs-couple --features superfluid --test superfluid_caloric
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
