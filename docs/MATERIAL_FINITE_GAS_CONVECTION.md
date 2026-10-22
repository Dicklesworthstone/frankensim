# Finite material gas/body convection

`fs_airflow::graph::thermal::medium::convection::finite_gas::ForcedGasSurface`
connects the existing material convection correlation to the original finite-body
enthalpy and gas/piston heat-work solve. Neither participant has an externally
fixed temperature. Call `surface.finite_body_candidate(&piston, piston_input,
body_input, piston_controls, heat_controls, cancelled)` and publish its two
endpoints together. Do not deposit the returned internal heat again.

The original correlation receives density, heat capacity, viscosity and thermal
conductivity from ONE source-bearing gas energy midpoint. h=Nu*k/L, G=h*area;
the original residual solves Q=dt*G*(T_body_mid-T_gas_mid) along with the original
piston pressure-work and momentum equations. One Q leaves the body and enters
the gas. Midpoints are in gas internal energy and body enthalpy, not averages of
endpoint temperatures. Latent heat stays in the existing enthalpy chart.

The public `resolve` function has the original conductance-callback signature;
reactive, powered or multi-boundary callers can use it without capturing a stale
initial atmosphere. Their coverage, whole-state publication and splitting rules
remain in force. This interface does not add separate chemistry or thermal steppers.
Actual hot/cold ordering selects Dittus-Boelter direction. A zero-area patch is
insulated and does not query an inapplicable flow correlation. The complete gas
sources, body phase state, transport snapshot and original coefficient evaluation
are retained in the accepted source sample.

## Limits

This interface explicitly uses BULK gas midpoint properties. Select a correlation
and temperature-contrast regime admitting that reference. Film-temperature
sampling, arbitrary large property variations, resolved boundary layers and
stagnant-gas fallback are not inferred. Flow speed is prescribed relative-flow
data, not piston speed or a solved pump. No jet thrust, fan energy, viscous
heating or turbulence field is added. The uniform body uses the original
negligible-pV-work enthalpy approximation. Melting does not silently change rigid
piston inertia or provide deformable/free-surface mechanics. Patch area must fit
the declared body surface and not be charged to another boundary law.

The original midpoint solver can overshoot equilibrium at large timesteps.
Refine time and separately justify uniform-body and correlation assumptions.
Conservation is necessary, not experimental validation. The heat-transfer
sign/area convention follows the documented control-volume wall balance:
https://www.cantera.org/stable/reference/reactors/interactions.html#walls
Cantera is a formulation reference, not a runtime dependency or the source of
this project's inertial-wall or correlation implementations.

## Focused checks

    cargo test -p fs-airflow --test finite_gas_convection

Seven authored tests use canonical manufactured species/caloric/transport packs:
pressure/composition sensitivity; reciprocal heat/motion energy; melting latent
heat; zero-area equivalence; geometry/regime/late-cancellation failure; repeated
material feedback; and actual direction selection. There is no new dependency.
Native Rust compilation/tests remain UNVERIFIED: Cargo exited 127, cargo/rustc
are absent. The working environment reset after the upload block; source was
recovered from session records, and the original parent/source hashes matched.

Independent Python equations checked 256 random exchanges, 64 repeated
steps and a melting case. Maximum step energy discrepancy is about 1.15e-13 J.
Changing only A mole fraction from 0.75 to 0.25 changes final 32-step wall velocity
from 0.0001314 to 0.0001995 m/s. These are manufactured equation checks, NOT Rust
execution or named-material measurements. Patch-application and source-byte checks also passed.
