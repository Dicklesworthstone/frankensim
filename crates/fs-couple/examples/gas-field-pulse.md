# A pressure pulse through a source-backed gas field

[gas_field_pulse.rs](gas_field_pulse.rs) runs the existing DG time-domain duct
with the same `AxialGasField` used by the frequency-domain and reed paths.
The supplied [gas-field-pulse.pulse](gas-field-pulse.pulse) has a 0.57 m cylinder,
a 1 Pa pressure pulse initially at rest, a closed inlet, and a matched outlet.
Nitrogen is pure at the inlet; carbon dioxide reaches mole fraction 0.5 at the
outlet. Temperature is 300 K, 320 K and 310 K at the three stations.

The material owner reevaluates density, sound speed and transport from the
explicit species, NASA-9 and pure-transport packs. DG uses local density and
bulk modulus in its mass matrices and the exact outlet gas for the matched
boundary. This example adds no gas equation or species-name constants.

## Run

From the repository root, build the six nitrogen/carbon-dioxide packs using
**Build the supplied source packs and render** in
[reed-axial-gas-field.md](reed-axial-gas-field.md#build-the-supplied-source-packs-and-render).
The same `target/reed-field-packs` directory serves this example:

```sh
cargo run -p fs-couple --example gas_field_pulse -- \
  crates/fs-couple/examples/gas-field-pulse.pulse target/gas-field-pulse.csv \
  --gas-component target/reed-field-packs/nitrogen.pack \
    target/reed-field-packs/nitrogen-nasa9-mcbride2002.pack \
    target/reed-field-packs/nitrogen-transport-lemmon2004.pack \
  --gas-component target/reed-field-packs/carbon-dioxide.pack \
    target/reed-field-packs/carbon-dioxide-nasa9-mcbride2002.pack \
    target/reed-field-packs/carbon-dioxide-transport-laesecke2017-huber2016.pack \
  --gas-water none --gas-phase frozen
```

Choose a new output path for each run. Parsing, material evaluation, propagator
admission and simulation finish before the CSV is created. An existing path is
never overwritten. Source errors and unresolved propagation fail without
producing a numerical output file.

## Input

These whitespace-delimited records occur in order; the version is literal.
Each station supplies one mole fraction per declared species in that order.

```text
gas_field_pulse_v1
cylinder RADIUS_M LENGTH_M
(or: cone INLET_RADIUS_M OUTLET_RADIUS_M LENGTH_M)
pressure STATIC_PRESSURE_PA
mesh ELEMENTS POLYNOMIAL_ORDER
time DT_S STEPS
pulse AMPLITUDE_PA CENTER_M HALF_WIDTH_M
species COUNT CANONICAL_ID...
stations COUNT
station X_M T_K MOLE_FRACTION...
```

After the last station, an optional record enables physical wall diffusion:

```text
viscothermal POLES MINIMUM_HZ MAXIMUM_HZ MAXIMUM_RELATIVE_ERROR
```

The supplied [gas-field-pulse-viscothermal.pulse](gas-field-pulse-viscothermal.pulse)
uses a 0.5 mm radius bore and a uniform 70% nitrogen / 30% carbon-dioxide
atmosphere. Substitute this input and a new CSV name in the command above;
the same source packs serve both examples. Density, viscosity, heat-capacity
ratio and Prandtl number select the positive wall-memory coefficients. Change
the static pressure, temperature or mixture at **both** stations to compare
the resulting attenuation and dispersion. No duct code or loss gain changes.

The `viscothermal` record uses the exact-mass uniform-cylinder path. For
nonuniform temperature/composition or a cone, use the same wall model with
explicit spatial quadrature:

```text
spatial_viscothermal POLES MINIMUM_HZ MAXIMUM_HZ MAXIMUM_RELATIVE_ERROR GAUSS_POINTS
```

[gas-field-pulse-spatial-viscothermal.pulse](gas-field-pulse-spatial-viscothermal.pulse)
supplies a 0.5-to-0.7 mm radius cone with three different nitrogen/carbon-dioxide
and temperature stations. It uses the same source packs and command. Gas
properties are reevaluated at each positive Gauss point; the source radius
sets the local viscous and thermal diffusion times. Pressure and every wall
memory advance in the original joint exponential. Changing the station data
therefore changes actual dissipation, dispersion, stored wall energy and the
pressure trace. This is a frozen gas field with rigid, isothermal walls and
the local-cylinder Webster approximation, not a full 3-D cone calculation.

GAUSS_POINTS must be at least ORDER+1 and at most 64. Refine it separately from
the mesh and pole count; its error is not included in the constitutive band
check. The example's 256-DOF cap includes every spatial wall state (the supplied
input has 216). The final diagnostic reports wall stored energy and the physical
instantaneous dissipation rate separately from the cumulative deficit.

The pole count
is a resolution axis (1 through 64), with complex wall response and real loss
checked independently against the exact Bessel owner over the requested band.
An unresolved count/band refuses. The complete acoustic and memory state uses
the existing DG exponential. Refine space separately: the declared wall-model
error does not establish spatial accuracy, and the initial pulse has broadband
content outside any finite admitted band. There is no automatic spectral filter.

Initial pressure is the DG projection of
`AMPLITUDE_PA * (1 - ((x-CENTER_M)/HALF_WIDTH_M)^2)^4` inside the support and
zero elsewhere. Initial volume flow is zero: the pulse splits in both
directions, and the left-going part reflects at the closed inlet. Its support
must lie strictly inside the duct. Choose an acoustic perturbation small
relative to static pressure, consistent with the linear acoustic equations.

Stations increase strictly from exactly zero to exactly the cylinder length.
Temperature and mole fractions interpolate; derived coefficients do not.
Authored fractions must sum to one under the material owner's canonical order
and are never repaired. Species IDs must match loaded pack associations.
Reordering IDs together with fraction columns preserves the field. Exact zero
fractions are retained, including a species absent at the inlet or everywhere.

The example bounds input to 64 KiB, 16 species, 64 stations, 64 elements,
order 1 through 8, 256 total DG degrees of freedom **including wall memory**,
10000 steps and 4 MiB of CSV.
The mesh must provide at least one element per material span. Source packs
retain the shared loader's 16 MiB per-pack and 64 MiB aggregate limits.
Increase element count or polynomial order to refine space. The certified
exponential advances the frozen operator; the time step sets output sampling.

## Read the CSV

There are `STEPS + 1` rows, including the projected initial state.
Values use 17 digits after the decimal point:

| Column | Meaning |
| --- | --- |
| `time_s` | Time since initialization |
| `inlet_pressure_pa` | Physical acoustic pressure at the inlet endpoint |
| `outlet_pressure_pa` | Physical acoustic pressure at the outlet endpoint |
| `stored_energy_j` | Current acoustic energy in the field-weighted DG state |
| `inferred_matched_loss_j` | Initial energy minus current energy |

With viscothermal memory the last column is named `inferred_total_loss_j`:
it includes both outlet absorption and wall dissipation. `stored_energy_j`
then includes the retained viscous and thermal memory as well as the acoustic
field. The two losses are not independently separated in this ledger.

Pressure uses LGL endpoint traces after physical-state reconstruction, with no
nearest-node choice or normalization. Matched loss is inferred from the energy
deficit. It is **not an independent conservation check**: the existing ledger's
`balance_checked` is false for this boundary.

For a uniform comparison, set all stations to 300 K and fractions `1 0`, retaining
both source IDs and the same mesh and pulse. Compare arrival times, reflections
and energy against the supplied gradient. At fixed temperature and composition,
this ideal-gas source model's sound speed does not change with static pressure;
density, characteristic impedance and the energy of a fixed pressure pulse do.

This is frozen, no-mean-flow plane-wave propagation, lossless in the interior
unless the optional viscothermal record is supplied. The matched outlet
absorbs energy; it is not a radiation receiver. Evolution of the equilibrium
temperature/species and molecular relaxation are outside this DG example. The source
manifests linked from the reed example retain references and validity domains.
