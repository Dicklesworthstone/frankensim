# Prescribed regional gas in valve performances

The same nonlinear valve, contact and solid-material state can now drive a duct
whose sections have different prescribed temperatures and humidity. Each section
uses its own gas-derived density and sound speed for characteristic impedance
and flight time, and its own viscosity, conductivity, heat-capacity ratio and
Prandtl number when viscothermal loss is selected. It is not a pitch adjustment
or an independently filtered voice. This advances the regional-gas part of MR33.

```bash
cargo run --release -p fs-couple --bin music_render -- \
  wind crates/fs-couple/examples/plate-valve-regional-gas.performance \
  /tmp/regional-valve.wav --decimate --block 37
```

The example retains the earlier illustrative plate, spatial lay, Maxwell solid
memory and pressure phrase. The inlet section is dry air at 293.15 K. The outlet
section is at 313.15 K and 50% relative humidity; a side branch and its enclosed
chamber are dry at 283.15 K. The common static pressure is 101325 Pa. These are
prescribed synthetic conditions, not a prediction of warmed-breath transport.
Both main sections select six existing viscothermal cells over 100..1000 Hz;
the narrow side branch explicitly does not select that wide-tube approximation.
The outlet uses the existing 1 kHz compact baffled radiation load and a receiver
0.2 m above the mouth. Mechanics run at 96 kHz, observed PCM at 48 kHz, with an
explicit 0.01 Pa full scale. There is no normalization or pressure retuning.

## Input

After a `duct_section` geometry, place an optional `gas TEMPERATURE_K RH`
before any optional viscothermal suffix:

```text
duct_section 1 2 0.125 0.007 0.002 gas 313.15 0.5 viscothermal 100 1000 6 8
```

All values are finite. RH is a fraction in [0,1], and the existing moist-air
model's further validity limits apply. Static pressure is inherited from the
file's `ambient` record: no unsupported mean pressure gradient is introduced.
Omitting `gas` uses that original ambient state. Invalid or repeated suffixes
refuse rather than selecting ambient. A single uniform `tube` still uses its
ambient record; split it into explicit network sections to supply a gradient.

The inlet section's gas also supplies the aperture's Bernoulli density and
characteristic impedance. There is no additional separately specified upstream
reservoir gas or density-mixing model. An enclosed cavity and radiating terminal
inherit the gas of their unique adjacent section. A wall's supplied solid law
is not implicitly changed by a nearby gas temperature.

Subdivision retains a section's full gas on every generated interval. Original
node addresses and each original integer transit are preserved; losses receive
represented local lengths, not an averaged sound speed. The ordinary uniform
construction keeps its original arithmetic. The same input works with
`ensemble --valve`: per-source scales remain encoding declarations, not gains.
Existing cancellation, finite-window and failed-callback behavior is retained.

## Library and observation

`ApertureNetwork::with_section_gases(valve, graph, gases)` accepts exactly one
`GasState` per section, in section order. `graph.inlet_impedance_with_gases()`
provides the physical load for constructing the valve. The graph speed field
must equal its inlet gas speed, and the valve density must equal inlet density;
no legacy input is silently ignored. States remain immutable through playback.
`section_gases()` and `section_medium()` expose the actual retained media.

`with_regional_viscothermal_sections` expands the same physical loss cells as
the uniform owner and returns the complete new gas mapping with its graph and
reports. Bind that result with the regional constructor; the scalar uniform
constructor is not a substitute. The existing finite-band, shear, temporal,
spatial, fit and memory gates are unchanged.

The existing Rayleigh receiver uses the **observed outlet's** medium, not the
inlet's. This assumes a homogeneous exterior continuation of that outlet gas.
No independent exterior-temperature field, interface refraction or room mixing
is modeled. Radiation-load and receiver geometry assumptions remain unchanged.
The sidecar records actual section gas values and the actual outlet medium.

## Physical scope

This is linear acoustics about a frozen piecewise-uniform, common-pressure gas
field. The existing pressure/volume-flow junction performs the impedance-step
scattering and work accounting. See J. O. Smith, *Physical Audio Signal Processing*,
“Lossless Scattering,” for that parallel acoustic junction relation.

No mean flow, advection, thermal evolution, diffusion, phase change, turbulent
jet mixing or moving thermodynamic interface is computed. GasState transport
and saturation approximations retain their original limitations. The wide-tube
loss and compact radiation declarations are finite-band sampled checks, not
bounds on all nonlinear transient content or measured-instrument calibration.
The plate is still one linear mode. The example demonstrates causal parameter
wiring, not a complete thermofluid instrument reconstruction.

The native regression targets remain `plate_aperture` and `music_render_wind`.
They compare raw network and direct physical constructions, preserve uniform
parity, exercise delayed regional feedback, check the local loss/load/receiver
medium, and cover cancellation, malformed inputs and actual CLI PCM.
