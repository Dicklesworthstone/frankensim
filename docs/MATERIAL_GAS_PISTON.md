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

This is one spatially uniform, closed ideal-gas cavity and one translational
inertial wall. It is not combustion, valves, gas exchange, friction, piston-wall
contact, a crankshaft, radiation feedback, a shock solver or an exact-isentropic
integrator. Entropy is reported by the gas owner, not certified nondecreasing.
A material source model is not experimental qualification of a real engine.

Six focused Rust regressions cover stationary equilibrium, heat-driven motion,
source/composition swaps, whole-trajectory energy with external loading, the
shared acoustic/mechanical small-signal period, and refusal/cancellation rollback.
Run:

```sh
cargo test -p fs-couple --lib cavity_phs::piston::tests
```

The synthetic source fixture uses fs-thermochem as a test-only dependency; no
production dependency was added. Native tests and Cargo lockfile regeneration
must be run with the workspace toolchain before claiming native verification.
