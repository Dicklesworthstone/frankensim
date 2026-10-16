# Named gas mixtures in wind performances

`frankensim-plate-valve-performance-v1` admits named frozen-composition mixtures
at both the ambient/default and individual duct-section positions. The same
input works with `music_render wind` and `ensemble --valve`. It changes the
physical medium, not a pitch shifter, post-filter or instrument-specific voice.

## Input and compatibility

The existing forms are unchanged:

```text
ambient 293.15 101325 0
duct_section 1 2 0.125 0.007 0.002 gas 313.15 0.5 viscothermal 100 1000 6 8
```

Ambient is temperature in kelvin, total static pressure in pascals, then relative
humidity as a fraction. The legacy regional `gas T RH` explicitly means moist
air at the ambient pressure, even when another region uses a named mixture.
Zero humidity still takes the original dry-air arithmetic path.

The additive mixture form is:

```text
ambient mixture TEMPERATURE_K PRESSURE_PA BASIS COUNT SPECIES FRACTION ...
duct_section ... gas mixture TEMPERATURE_K PRESSURE_PA BASIS COUNT SPECIES FRACTION ...
```

`BASIS` is exactly `mole` or `mass`. Supply one positive fraction for each of
1..5 distinct species; the shared composition owner requires a unit sum. These
are fractions, not percentages. There is no implicit balance gas, normalization,
or volume-fraction reinterpretation. Omit zero components. Supported symbols and
aliases are `N2`/`nitrogen`, `O2`/`oxygen`, `Ar`/`argon`,
`CO2`/`carbon-dioxide`, and `H2O`/`water-vapor`; names are case-sensitive.

```text
ambient mixture 293.15 100000 mole 2 N2 0.7 O2 0.3
duct_section 1 2 0.125 0.007 0.002 gas mixture 303.15 100000 mole 5 N2 0.70 O2 0.20 Ar 0.05 CO2 0.04 H2O 0.01 viscothermal 100 1000 6 8
```

A mass-basis alternative, in the ambient position, is:

```text
ambient mixture 300 100000 mass 3 nitrogen 0.75 oxygen 0.20 argon 0.05
```

The gas declaration precedes the optional existing viscothermal suffix. Omission
copies the complete ambient state. A uniform `tube` uses that same ambient state;
use explicit network sections to vary composition or temperature spatially.
Repeated or malformed declarations never fall back to air.

## Material ownership and refusal boundaries

`fs_material::gas_species::{GasSpecies, resolve_gas_mixture}` owns the reusable
named data selection and domain admission. It delegates constitutive evaluation
to the existing `GasState::try_new_mixture`, with composition semantics and
mass/mole conversion owned by `fs_qty::semantic::Composition`. The input parser
does not compute density, heat capacity, sound speed or transport coefficients.

The named convenience model admits **273.15..=323.15 K** and
**10000..=200000 Pa**. These are conservative declared model-use limits, NOT a
measured accuracy certificate. N2/O2/CO2 heat capacities are frozen at their
NIST-JANAF 298.15 K values; argon is a calorically perfect monatomic gas.
Sutherland coefficients and their primary references are recorded in
`fs-material/src/gas_species.rs`. H2O reuses the existing water-vapor model;
conductivity uses the existing Eucken approximation, and mixtures use the
existing Wilke/Wassiljewa-Mason-Saxena owners. This does not create measured
material-card or transport-accuracy authority. The programmatic supplied-GasSpec
mixture API and legacy humid-air constructor retain their own documented domains.

Water partial pressure is checked with the existing saturation-pressure law.
Supersaturation is refused because condensation is not modeled. A physically
admitted pure-vapor state is not rejected merely because its water mole fraction
is one. Unknown species, duplicates, missing/extra fields, invalid fractions,
nonfinite values and out-of-domain states produce explicit source-line errors;
material failures retain their original cause.

**All connected sections must share the ambient static pressure**, including
when every section supplies an explicit override. This is a stationary acoustic
network, not a mean-flow solver. Change the common pressure for the whole render;
do not prescribe an unsupported steady pressure jump across a junction.

## One resolved state through the instrument

The inlet state supplies Bernoulli density and characteristic impedance. Each
duct section uses its own sound speed and impedance for flight time and junction
scattering. Viscothermal subdivisions retain the complete local state and use
its viscosity, conductivity, heat-capacity ratio and Prandtl number. A cavity or
radiation boundary inherits its unique adjacent section's gas. The existing
Rayleigh receiver uses the observed outlet's density and sound speed, assuming
a homogeneous exterior continuation of that gas. No gas law is reimplemented
in any of these source adapters.

Authored solid material, wall laws and explicit lumped resistance values are not
silently retuned from nearby gas conditions. Air-only ISO 9613 absorption is not
made valid for an arbitrary mixture by this input extension.

`PlateValvePerformance::ambient_gas()` exposes the admitted default.
`ApertureNetwork::section_gases()` exposes actual regional states, including
loss-cell mappings. WAV sidecars now retain an `ambient_gas` snapshot on both
uniform-tube and graph paths; regional snapshots include heat capacity, specific
gas constant and characteristic impedance as well as the previous properties.
The existing input hash binds the exact species, basis and fractions.

## Controlled sound comparisons

Three complete fixtures use the identical plate, closure, solid memory, cavity,
loss cells, mouth-pressure phrase, receiver, clocks and 0.01 Pa PCM full scale:

| Fixture | Mole fractions | Static pressure |
| --- | --- | --- |
| `plate-valve-mixture.performance` | N2 0.7, O2 0.3 | 100000 Pa |
| `plate-valve-mixture-co2.performance` | N2 0.4, O2 0.2, CO2 0.4 | 100000 Pa |
| `plate-valve-mixture-low-pressure.performance` | N2 0.7, O2 0.3 | 50000 Pa |

All three use 293.15 K. Render from the repository root:

```bash
for case in mixture mixture-co2 mixture-low-pressure; do
  cargo run --release -p fs-couple --bin music_render -- \
    wind "crates/fs-couple/examples/plate-valve-${case}.performance" \
    "/tmp/plate-valve-${case}.wav" --decimate --block 37
done
```

Composition changes density, sound speed, heat capacity and transport. At fixed
composition and temperature, the ideal-gas model makes sound speed independent
of static pressure: reducing pressure changes density, impedance, gas loading
and boundary-layer response, not an invented pitch-control parameter. No output
normalization or gesture retuning is applied. These are synthetic engineering
comparisons, not reconstructed or experimentally validated instrument timbres.

Focused native tests:

```bash
cargo test -p fs-material --lib gas_species:: -- --test-threads=1
cargo test -p fs-couple --lib performance::file::gas::tests -- --nocapture --test-threads=1
cargo test -p fs-couple --test music_render_wind gas_mixture:: -- --test-threads=1
```

The coupled tests compare every sample and aperture state against independent
assembly through the existing physical owners, check discrete energy balance,
legacy parity, mass/mole equivalence, local loss/cavity/radiation/observation
states, changed rendered pressure, callback partitioning and PCM cancellation/
resume. The CLI tests compare actual WAV bytes with the library PCM owner,
exercise all five species in a regional outlet and two-mixture ensembles, and
require invalid inputs to fail before output creation. `--nocapture` prints the
physical parameters and actual raw-pressure RMS from the coupled comparison.

## Frozen fields, not evolving gas transport

Temperature, composition and static pressure are fixed after admission. This
does **not** implement time-varying gas replacement during one playback,
advection, diffusion, reactions, real-gas thermodynamics, condensation, gas/solid
thermal evolution or an independent exterior-composition field. The existing
finite-band loss, compact-radiation, receiver, geometry and time-resolution
refusals remain active; valid material data do not override them.
