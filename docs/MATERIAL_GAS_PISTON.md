# Material-driven gas and piston dynamics

`fs_couple::cavity_phs::piston::GasPiston` couples an existing
`ClosedGasParcel` to a constant-area translating inertial wall. This is a
nonlinear closed-cavity component alongside the existing linear vented
plate/Helmholtz cavity; it does not replace that model or pretend to add a neck.

Supply the gas parcel, piston mass [kg], area [m2] and velocity [m/s] to
`GasPiston::new`. Each `step` receives explicit duration, integrated gas heat,
constant exterior pressure and additional applied force, plus a final-volume
bracket, gas solve controls, mechanical impulse tolerance, combined energy
tolerance and iteration cap. Positive motion increases gas volume.

The momentum solve calls the existing
`ClosedGasParcel::midpoint_heat_volume` material update. That owner uses
midpoint **internal energy and volume**, not averaged temperature or fixed
gamma. Its actual midpoint pressure and discrete volume derivative supply the
load and Newton tangent. The mechanical equations are

`V1-V0 = A*dt*(v0+v1)/2`

`M*(v1-v0) = dt*(A*(p_mid-p_exterior)+F_external)`.

The accepted states independently close

`delta_U_gas + delta_K_piston + p_exterior*delta_V - F_external*delta_x = Q`.

Gas work is not applied twice, and exterior-pressure work is not lost. Every
material trial restarts from the accepted gas; no trial changes persistent
state. Source refusal, unsupported volume/temperature, exhausted work, failed
energy admission and cancellation leave both the gas and velocity unchanged.
Volume is the single position owner; the effective wall coordinate is V/A.
The returned report retains the complete midpoint source state for thermal,
transport and acoustic consumers as well as actual mechanical/thermal balances.

Changing the supplied gas composition or caloric material data changes the
actual trajectory and its pressure/temperature/sound-speed response, not only
diagnostics. No gas name selects a solver branch. Thermodynamic/calibration
support stays with fs-material/fs-thermochem, while the moving wall and coupled
momentum/energy step live in fs-couple.

The ordinary `step` is one spatially uniform, closed ideal-gas cavity and one
translational inertial wall. `step_with_nozzle`, described below, adds conservative
exchange with a finite reservoir. Neither models combustion, valve actuation,
friction, piston-wall contact, a crankshaft, radiation feedback, shocks or exact
isentropic integration. Entropy is reported by the gas owner, not certified nondecreasing.
A material source model is not experimental qualification of a real engine.

Six focused Rust regressions cover stationary equilibrium, heat-driven motion,
source/composition swaps, whole-trajectory energy with external loading, the
shared acoustic/mechanical small-signal period, and refusal/cancellation rollback.
Run:

```sh
cargo test -p fs-couple --lib cavity_phs::piston::
```

The synthetic source fixture uses fs-thermochem and fs-ga as test-only
dependencies. The material-bound adapter consumes the existing fs-mbd-material
production owner rather than adding another density-to-inertia implementation. Native tests and Cargo lockfile regeneration
must be run with the workspace toolchain before claiming native verification.

## Solid-material and composite piston binding

`piston::material::MaterialGasPiston::new(gas, body, area, velocity)` takes the
existing `fs_mbd_material::ResolvedBody` produced by `resolve_body`. It uses the
canonical spatial inertia's total mass and retains every region's geometry,
source density, conditions and usage receipt. There is no separate mass knob
that can silently disagree with the source. The body is immutable during
evolution; `step` delegates to the exact same gas/piston dynamics above.

Swap the rigid specimen's density cards or regional material assignments and
recompile that body to change actual piston inertia, acceleration, boundary
work and the gas's resulting thermal/acoustic state. Composite regions use the
existing body's canonical reduction, not another mass calculation. Gas and
solid source conditions remain independent; heating the gas does not imply
that the rigid piston heats or expands. Instantaneously changing a live mass
is not modeled, because that needs an explicit momentum/energy transfer law.

The exposed area is an independent geometry quantity, not the specimen volume
or gas cavity volume. Translation is constrained: rotational inertia and
centre-of-mass offsets are retained on the body but do not generate rotation.
Three additional Rust tests exercise density-only trajectory changes, composite
source retention and exact agreement with the numerical mass path, and atomic
failure/cancellation including immutable solid-source history.

## Conservative nozzle exchange and moving walls

`GasPiston::step_with_nozzle` accepts a mutable finite `ClosedGasParcel` reservoir,
the existing `GasNozzle` area/discharge coefficient and `GasNozzleStepSolve`, the
ordinary piston drive/controls, and a whole-pair energy tolerance. The
`MaterialGasPiston` adapter exposes the same operation with its source-resolved
solid inertia. Both gases and the piston publish together only after all source,
inventory, clock, momentum and combined energy checks succeed.

The upstream side is selected from the actual initial pressures. The original
material nozzle computes subsonic or choked flow using the source's enthalpy,
entropy and sound speed. Its paired transfer moves actual species and upstream
total enthalpy at fixed volumes. The existing midpoint piston step then uses
that newly resolved chamber composition, density, temperature, pressure,
transport and acoustic state. Incoming enthalpy is already present; supply only
independent external heat in `PistonInput::heat_j`.

This is an explicitly first-order transfer-then-motion split over one physical
time interval. Refine that interval to assess splitting error. The nozzle's
donor-fraction and pressure-crossing checks remain active; no transfer is clipped
to fit them. The reservoir is a well-mixed fixed volume and incoming jet energy
thermalizes there or in the chamber. No resolved jet thrust, shock, valve motion
or simultaneous compressible-flow solution is implied.

`PistonNozzleStep` preserves both original stage reports and independently checks
`delta_U_chamber + delta_U_reservoir + delta_K + exterior_work - force_work - Q`.
Total gas mass and the transfer owner's per-species balances remain explicit.
The transferred enthalpy cancels inside this two-gas accounting boundary.

After nonzero mass exchange, the piston still exposes its real wall and heat
transfers, but its existing two-port `PortedModule` window audit is unavailable:
that schema lacks the stream port needed for an open chamber. The returned
paired report supplies the complete energy/mass account for this composition.
A closed nozzle retains the ordinary closed-step audit and arithmetic. Failed
transfer, failed subsequent motion and late cancellation preserve both parcels,
the wall state and the prior accepted port account.

## Nozzle exchange with a finite material wall

`GasPiston::nozzle_thermal_candidate` joins the same nozzle transfer to the
existing reciprocal thermal/piston solver. It accepts an immutable reservoir,
the ordinary nozzle and piston controls, a `UniformEnthalpyStepInput` for the
finite wall or lining, `HeatExchangeSolve`, and the existing conductance callback.
The returned `NozzleThermalPistonCandidate` contains the piston, reservoir and
finite thermal endpoint. Publish all three returned endpoints together.

The nozzle stage runs once from the accepted initial pressures. Every trial of
the subsequent coupled heat/momentum solve starts from that same transferred
chamber. Incoming species therefore affect the gas midpoint's caloric, transport,
pressure and acoustic properties used by the wall-heating and motion equations.
Thermal iteration cannot withdraw or deposit nozzle enthalpy a second time.

`nozzle_thermal_layers_candidate` supplies the existing
`MaterialConductionLayers` callback directly. The wall's conductivity comes from
the same source card as its enthalpy/phase chart; gas conductivity comes from the
actual transferred mixture at its solved energy midpoint. Changing supported
conductivity or caloric data, inlet composition or the layer geometry changes
heat exchange, finite-wall temperature/phase, piston motion and the resulting
gas state through those original owners. `MaterialGasPiston` exposes both methods
and retains its original source-resolved rigid-body inertia in the candidate.

The full accepted balance is

`delta_U_chamber + delta_U_reservoir + delta_H_wall + delta_K + exterior_work - force_work - Q_gas_external - Q_wall_internal`.

Both nozzle enthalpy and reciprocal wall heat cancel inside this boundary.
`HeatExchangeSolve::energy_tolerance_j` admits this complete balance as well as
the original thermal stage. The report retains original nozzle/species and
thermal/momentum results plus independently measured endpoint storage and mass.
Nonzero flow keeps the same explicit missing-stream port-audit boundary described
above; it never emits a falsely complete closed-chamber window.

The transfer-then-coupled-heat/motion composition is first order over one physical
interval and requires timestep refinement. Its finite thermal store uses the
existing fixed-volume, negligible-pV-work enthalpy model and can retain latent
heat. That store may be a lining separate from the rigid moving specimen; its
temperature or phase does not silently update rigid-body shape or mass. Missing
material support, a failed thermal trial, failed energy admission or late
cancellation returns none of the tentative endpoints. A closed nozzle recovers
the original thermal solve, and an insulated wall recovers the original nozzle
and piston step.

The focused regression module is
`cavity_phs::piston::exchange::thermal::tests`. It exercises conductivity and inlet
composition changes, both limiting cases, finite latent storage with independent
whole-trajectory energy, post-transfer refusal/cancellation/retry and refinement.
