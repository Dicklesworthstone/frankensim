# Independent exterior gas for material plate-valve performances

The optional `exterior` record belongs immediately after `ambient`, before
`tube` or `network`. Both records use the SAME material-owned parser:

```text
ambient mixture 293.15 101325 mole 2 N2 0.7 O2 0.3
exterior mixture 293.15 101325 mole 3 N2 0.4 O2 0.2 CO2 0.4
```

`mass` selects mass fractions. The legacy `exterior T_K P_PA RH` form supplies
humid air, with explicit absolute pressure and fractional relative humidity.
Existing named-species, composition, saturation and temperature/pressure limits
apply. Unsupported states never fall back to air. All internal regions and the
exterior must share one static pressure; this is not a mean-flow interface.

## One declaration, two physical consumers

Every explicitly declared `baffled-low-ka BAND_HZ` endpoint uses its ACTUAL
adjacent radius and the exterior gas in the existing `BaffledRadiationLoad`.
The selected baffled receiver uses the same gas in the existing Rayleigh
integral and its propagation delays. The exterior is not a gain or a replacement
for the gas inside the bore. Its load reacts on the actual traveling waves and
moving plate, including the original spatial contact and material memory.

The aperture and each duct section retain their internal density, sound speed,
viscosity, heat capacities and conductivity. Enclosed cavities keep their own
adjacent gas. Wall laws, explicit lumped resistances, solid material, geometry,
initial state, pressure gestures and mechanical forces are unchanged. No new
constitutive law, solver, radiation integral or PCM encoder is introduced.

Omitting `exterior` takes the original local-gas path verbatim. An explicit
matched exterior preserves the regional reference sound. Internal-pressure
observations may also select exterior-derived reaction without a microphone.
The library's `AperturePerformance::new_with_exterior_gas` is observation-only;
it does not silently replace an arbitrary already-authored load. The file
adapter is what binds the common exterior to its declared radiation loads.

An exterior record requires at least one declared baffled radiation endpoint.
Numeric tube reflections are not replaced. A baffled observation with an
explicit exterior must name a declared radiation endpoint, not an enclosed
cavity, junction or arbitrary impedance. Duplicate exterior records refuse.
All input, geometry, material and radiation checks complete before CLI files
are created. Cancellation returns no partially running performance.

## Controlled comparison

`plate-valve-exterior.performance` is the existing regional-gas material-valve
fixture with ONE added exterior record. The humid hot outlet, colder cavity,
viscothermal sections, plate and material memory are the original inputs.
The 0.01 Pa PCM scale is identical, not separately normalized.

```bash
cargo run --release -p fs-couple --bin music_render -- \
  wind crates/fs-couple/examples/plate-valve-regional-gas.performance \
  /tmp/plate-inside.wav --decimate --block 37
cargo run --release -p fs-couple --bin music_render -- \
  wind crates/fs-couple/examples/plate-valve-exterior.performance \
  /tmp/plate-outside.wav --decimate --block 37
```

`ensemble --valve` uses the same finite pressure producer. Its source PCM scales
are not gains. Both CLI routes report the resolved exterior separately from
ambient and section states, and report the ACTUAL receiver density/sound speed.
Physical propagation is not removed by output resampling or ensemble alignment.

Focused checks use the existing test seams:

```bash
cargo test -p fs-couple --lib bernoulli_aperture::performance::exterior_tests
cargo test -p fs-couple --lib performance::file::duct::exterior::tests
cargo test -p fs-couple --test music_render_wind gas_mixture::exterior
```

These tests are authored assertions, not evidence of a passing native run in
the implementation session. They cover direct owner equality, unchanged internal
coefficients, exterior feedback into plate motion, matched-medium/callback parity,
all declared endpoints, actual-exterior radiation admission, CLI/ensemble PCM,
pressure-only comparisons at fixed mole fractions and refusal before output.

The model remains stationary prescribed atmospheres, independent compact
baffled loads and ONE observed outlet. Existing low-ka and numerical band limits
remain active. There is no mutual radiation reaction, room/body scattering,
thermal/composition transport, live gas replacement, or measured-timbre claim.
