# Finite material thermal storage coupled to gas and piston motion

`fs_couple::cavity_phs::piston::GasPiston::thermal_candidate` connects the
existing gas/momentum step to a finite thermal body using the canonical
`fs_material::phase::UniformEnthalpyStepInput` and `UniformEnthalpyStep`.
The body is an enthalpy/phase chart, not a constant-Cp temperature increment.
A signed interface heat Q is added to the gas and subtracted from the body.

The nonlinear solve evaluates BOTH materials at their energy midpoints:

- gas: the existing internal-energy/volume midpoint and pressure-work solve;
- thermal body: `h_mid = h_old + (Q_body_internal - Q)/(2*m_body)`;
- interface: `Q = dt * G(gas_mid, body_mid) * (T_body_mid - T_gas_mid)`.

An isothermal latent plateau changes phase fraction rather than inventing a
large heat capacity or forcing temperature to rise. Independent gas and body
heat sources remain separate from the internal exchange. Each trial reuses the
same old gas/momentum state. Absolute source, body-energy and combined-energy
residuals must pass; no bracket-width-only or frozen-property fallback is used.

## Source-bound conductive layers

`thermal::layers::MaterialConductionLayers` supplies a concrete interface:

`G = area / (gas_distance / k_gas + body_distance / k_body)`.

The gas coefficient comes from its actual midpoint. The body coefficient comes
from the SAME material card as the body's enthalpy chart, queried at actual
midpoint temperature and, when explicitly declared, liquid mass fraction. All
other source conditions, quantity descriptors, exact property keys and claim
selection are preserved. Scalar and declared T/fraction curves are supported;
frequency-dependent coefficients cannot masquerade as static conduction.
Missing data, wrong card identities or unsupported states refuse, even with
zero area. Explicit zero body conductivity or zero area means insulation.

```rust,ignore
use fs_couple::cavity_phs::piston::thermal::layers::MaterialConductionLayers;

// card, context, conductivity_key, body_input, input and both solve controls
// are caller-owned existing source/geometry/thermal types. body_input carries
// the original phase chart, accepted state, conserved mass and physical clock.
let path = MaterialConductionLayers {
    material: &card, context: &context, temperature_axis: "T",
    liquid_fraction_axis: None, conductivity: &conductivity_key,
    selection, area_m2: 0.02, gas_distance_m: 0.001, body_distance_m: 0.02,
};
let proposal = piston.thermal_candidate(input, body_input, piston_controls,
    exchange_controls, |gas, body| path.resolve(gas, body), || false)?;
let (next_piston, next_thermal, report) = proposal.into_parts();
// Publish both endpoints in the application's existing physical state owner.
// next_thermal.report retains the exact conductivity source bundle.
// Do NOT also apply next_thermal.external_heat_j to next_thermal.state:
// that signed boundary heat is already included in its accepted enthalpy.
```

Use the generic callback instead for an existing contact, convection or other
interface owner; return its concrete source record in `ConductanceSample`.
Callbacks must be repeatable and must not debit a live external reservoir.
No input state is mutated. Publish both candidate endpoints together only after
all other coupled owners have accepted. Cancellation may use an existing gate.

## Scope and verification

This is a uniform thermal body with fixed geometry and negligible pV/deformation
work, a closed frozen-composition gas, and one translating inertial wall. The
thermal body may be a separate chamber wall or store; it is not implicitly the
piston. Its phase change does not change piston mass, shape or strength. The
layer lengths and spatially uniform coefficients are explicit approximations,
not resolved conduction profiles, convection, contacts or moving phase fronts.
Midpoint is not unconditionally monotone or a nonnegative entropy-production
certificate; timestep refinement and material support remain necessary.

Six thermal-coupling and four source-layer Rust tests are included. Run:

```sh
cargo test -p fs-couple --lib cavity_phs::piston:: -- --test-threads=1
```

The existing instrument-material workflow runs this focused battery. Local
Rust execution was unavailable during authoring; independent Python equations
checked finite heat capacity, latent storage, two heat sources, data-only
conductivity changes and actual coupled energy/motion. Those synthetic equation
checks are not native test passes or experimental validation.
