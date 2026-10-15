# Source-driven layered sections: electricity, heat and sound

`fs_material::section_temperature::ConcentricSectionTemperatureLaw` retains
one original material card and source-claim selection per layer. A temperature
update re-evaluates every requirement in the original bundle, not only Young's
modulus. Mechanical, thermal, electrical and auxiliary properties keep the same
source selection and non-temperature coordinates. Unsupported states refuse.

`section_temperature::heating::HeatedConcentricSection` connects those updates
to the existing `fs-thermal` enthalpy timeline. Each annulus retains its original
mass and its source-backed Cp/density chart. Every thermal substep uses one
material evaluation for per-layer Joule heat and link conductance, then refreshes
the endpoint mechanical section. A failed advance returns no candidate. Generated
heat does not become a persistent authored source, and the thermal audit already
includes it: never add the electrical-energy report to that audit again.

## Run the actual acoustic consumer

From the configured FrankenSim workspace:

```sh
cargo run -p fs-couple --example heated_section_sound -- \
  --demo 2.0 0.2 1.0 0.01 comparison
```

The arguments are sleeve electrical-resistivity scale, terminal voltage [V],
heating duration [s], PCM full scale [Pa], and output prefix. The core is unchanged.
The demo is synthetic data, not a claim about gold, silver, steel or a measured
specimen. Repeat with another sleeve scale to change data without editing physics.

For real input, use two existing canonical `NormalizedMaterialCardPack` files:

```sh
cargo run -p fs-couple --example heated_section_sound -- \
  --packs core.pack sleeve.pack 0.2 1.0 0.01 comparison_from_packs
```

The example's declared source context is `T` in kelvin over [300,400] K; it
requires positive density, Young's modulus, specific heat, thermal conductivity
and electrical resistivity. Density must be constant on that interval. Missing,
ambiguous, differently conditioned or unsupported data refuse; no nominal-metal
fallback is used. The library API additionally supports full caller-supplied
nonthermal context and per-layer absolute-temperature axis names.

The same section geometry, fixed 10 N tension, pluck, ambient gas, listener and
acoustic renderer produce the two signals. `comparison.csv` contains physical
pressure samples for cold and heated states. `comparison.wav` is a stereo
comparison, **cold on the left and heated on the right**, using the existing
PCM encoder with one shared explicit full scale. No separate peak normalization
hides material effects. Clipping is reported. Existing output files are not
overwritten. The example also reports layer temperatures, axial stiffness,
integrated electrical heat, stored-energy change and the thermal residual.

These are two frozen-state 50 ms sound observations, not continuous sound during
heating. The model is a perfectly bonded concentric slender member, not a helical
winding model. Each annulus has one temperature; interface area and half-path
thermal resistance geometry are explicit modelling inputs. This rung does not
include thermal expansion/work feedback, density-driven geometry changes, phase
transitions, skin effects, or a resolved intralayer temperature field. The existing
conduction CFL does not certify nonlinear electrical-source stability; independently
refine the requested maximum step. Source/model uncertainty remains upstream.

Focused native checks (not executed in the connector-only implementation session):

```sh
cargo test -p fs-material section_temperature
cargo test -p fs-couple --example heated_section_sound
```
