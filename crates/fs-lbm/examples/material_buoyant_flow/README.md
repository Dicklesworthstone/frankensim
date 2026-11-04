# Material-driven buoyant flow with physical heating

This runnable example connects one material card to the existing local density,
viscosity, conductivity, caloric, momentum and enthalpy owners. Change the input
file to change the material response; no solver code changes are needed.

```bash
cargo run -p fs-lbm --features enthalpy-coupling --example material_buoyant_flow -- \
  crates/fs-lbm/examples/material_buoyant_flow/expanding.txt 40 400

cargo run -p fs-lbm --features enthalpy-coupling --example material_buoyant_flow -- \
  crates/fs-lbm/examples/material_buoyant_flow/weaker_expansion.txt 40 400
```

Arguments are `MATERIAL.txt [STEPS=40] [HOT_K=400]`. Each reported interval is
0.05 seconds; the original thermal owner may shorten internal steps for its
conductive and advective CFL limits. The accepted momentum and material history
survive every interval.

## Material input

The first non-comment row is the constant specific heat capacity in J/(kg K).
It is followed by 2–64 rows with four SI columns:

| Column | Meaning |
| --- | --- |
| T | Absolute temperature, K; strictly increasing |
| rho | Density at that temperature, kg/m³ |
| mu | Dynamic viscosity, Pa s |
| k | Thermal conductivity, W/(m K) |

All values must be finite and positive. Input is bounded to 64 KiB. The example
compiles these values into a `MaterialCard`; the existing material evaluator
owns bounded interpolation and source selection. `LiquidTemperatureLaw` and
`ResolvedHeatCapacityCurve` read that same card at a declared pressure of
100,000 Pa. The example does not implement another property interpolator.

Both supplied files are **synthetic numerical inputs**, not measured data or
calibration for a named real liquid. Their Cp, viscosity and conductivity agree;
only the density curve changes. The weaker density response produces a weaker
buoyant drive under the same geometry, temperatures and gravity.

## Physical scenario and energy ownership

The fixed geometry is a deliberately coarse 2×2×1 grid of unit-metre cells. The
left column starts at the material's first temperature, and the right column at
`HOT_K`. The four side walls are stationary and impermeable. The two depth faces
are symmetry boundaries. Gravity is `[0, -9.81, 0]` m/s². There is no fan, heater,
external body actuator, or moving wall.

At every accepted source epoch:

1. The original caloric chart converts cell enthalpy to temperature.
2. The same material evaluation supplies local density, viscosity and k.
3. Boussinesq buoyancy uses `(rho(T) - rho_reference) * gravity`; cell mass and
   inertia retain the reference density. Density departures are limited to 5%.
4. The original native Stokes solver advances retained momentum and computes
   physical molecular loss using its own stress quadrature.
5. The original thermal solver deposits that physical loss once, conducts heat
   across the four interior faces using source k, and advects enthalpy with the
   accepted old velocity.

The CSV separates physical heat, native kinetic energy, backward-Euler numerical
loss and signed buoyancy work. Numerical loss is **not** deposited as heat.
Buoyancy work is already part of total body work, so the energy balance never
adds that subset twice. The final defect compares thermal change plus kinetic
storage plus numerical loss with the signed mechanical supply and bounded
pressure-constraint work. A defect above 1e-6 J refuses the run.

## Model limits and checks

This is a closed, single-phase, incompressible Boussinesq/unsteady-Stokes model
with a cell Reynolds ceiling of 0.1. Momentum convection is omitted. It uses the
existing first-order split and requires temporal and spatial refinement for an
accuracy claim. Buoyancy work is an explicit exchange with the declared forcing
model; a gravitational-potential reservoir is not resolved. The example does
not model melting, a free surface, turbulent mixing, compressible work, or a
calibrated real substance.

The example tests exercise actual source evaluation and both material files
through this same simulation entry. They check changed circulation and physical
heating, source-domain refusal and the combined energy account. The separate
library regressions independently reduce the native 2×2 mechanical equations.
Metadata compilation alone does not establish executed numerical-test results.
