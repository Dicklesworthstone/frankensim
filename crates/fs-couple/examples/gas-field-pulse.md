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

After the last station, an optional `steady_wall` block can first compute the
gas temperatures from an original wall material source, as described below.
After that block, or directly after the stations, an optional record enables
physical acoustic wall diffusion:

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
The mesh must provide at least one element per material span. Gas source packs
retain the shared loader's 16 MiB per-pack and 64 MiB aggregate limits.
Increase element count or polynomial order to refine space. The certified
exponential advances the frozen operator; the time step sets output sampling.

## Generate the field from wall material data

[gas-field-pulse-material-wall.pulse](gas-field-pulse-material-wall.pulse)
uses the same pulse renderer after the existing steady gas owner has solved
enthalpy transport and heat exchange with a material wall. The wall is the
repository's first-heating EN AW 2024-T3 source: its density, heat capacity and
conductivity curves cover 298.15–373.15 K. The conductivity and warm density
are derived quantities and Cp is a published CALPHAD calculation; these are
not new measured data or a general qualification of aluminum tubing. The
original notes and citations are retained in
[`aluminum-2024-t3-first-heating-km2024`](../../../data/matdb/seed-v1/aluminum-2024-t3-first-heating-km2024/manifest.tsv).

Compile the existing source into a material-card pack and copy the pulse input
beside it. Pack paths in a `wall_material` record are relative to that input;
the existing bounded native loader accepts material-card packs up to 8 MiB,
and the material record refuses paths escaping the input directory.

```sh
mkdir -p target/gas-wall-pulse
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/aluminum-2024-t3-first-heating-km2024/manifest.tsv \
  --out target/gas-wall-pulse/aluminum-wall.fsmatpk \
  --material-state EN-AW-2024-T3 solid first-heating 0
cp crates/fs-couple/examples/gas-field-pulse-material-wall.pulse \
  target/gas-wall-pulse/aluminum.pulse
```

Use the same six gas packs prepared above:

```sh
cargo run -p fs-couple --example gas_field_pulse -- \
  target/gas-wall-pulse/aluminum.pulse target/gas-wall-pulse/aluminum.csv \
  --gas-component target/reed-field-packs/nitrogen.pack \
    target/reed-field-packs/nitrogen-nasa9-mcbride2002.pack \
    target/reed-field-packs/nitrogen-transport-lemmon2004.pack \
  --gas-component target/reed-field-packs/carbon-dioxide.pack \
    target/reed-field-packs/carbon-dioxide-nasa9-mcbride2002.pack \
    target/reed-field-packs/carbon-dioxide-transport-laesecke2017-huber2016.pack \
  --gas-water none --gas-phase frozen
```

The diagnostic output records the material pack hash, state and fixed source
coordinates, the conservative steady power and energy residuals, and each
emitted station's temperature, density and sound speed. The CSV contains the
original DG pressure and energy observations for that generated field.

### Swap the material through data

The existing NASA generic 316 source supplies density, Cp and conductivity
over 273–400 K. Its source identifies neither a UNS designation nor a heat
treatment or product form. Those unknowns remain explicit in
[`stainless-316-nasa-tp216435`](../../../data/matdb/seed-v1/stainless-316-nasa-tp216435/manifest.tsv)
and in the required coordinates below.

```sh
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/stainless-316-nasa-tp216435/manifest.tsv \
  --out target/gas-wall-pulse/stainless-wall.fsmatpk \
  --material-state generic-316-stainless solid nasa-table23 0
cp target/gas-wall-pulse/aluminum.pulse target/gas-wall-pulse/stainless.pulse
```

In `stainless.pulse`, replace only the `wall_material` record with:

```text
wall_material stainless-wall.fsmatpk nasa-table23 5 source-generic-316-identity=1 source-uns-designation-known=0 source-heat-treatment-known=0 source-product-form-known=0 source-pressure-known=0
```

Run the same rendering command using `stainless.pulse` and a new
`stainless.csv` output path. Geometry, inlet state, flow, outer reservoir,
film and numerical controls stay identical. The original source law
determines the difference. These two conductive metals may produce a small
effect when the gas-side film dominates the resistance; no audible difference
or validated instrument response is asserted.

### Optional block grammar and physical scope

The block appears after the station rows and before any viscothermal record:

```text
steady_wall
wall_material PACK STATE AXIS_COUNT AXIS=VALUE...
wall_source TEMPERATURE_AXIS solid|liquid DENSITY_KEY CP_KEY CONDUCTIVITY_KEY MIN_K MAX_K CALORIC_INTERPOLATION_TOLERANCE_J_KG
wall_geometry THICKNESS_M FILM_AREAL_RESISTANCE_M2_K_W
wall_reservoir EXTERNAL_T_K
wall_flow MASS_FLOW_KG_S
gas_solve MIN_T MAX_T H_TOL MAX_ITER DT MAX_RELAX MAX_SUBSTEPS POWER_TOL ENERGY_TOL MAX_MACH
wall_solve MAX_ITER ABSOLUTE_HEAT_TOL_W RELATIVE_HEAT_TOL
```

Property keys and coordinates are the **normalized names in the compiled
source**, including hyphens. Their registered quantity semantics are retained;
the original caloric and conductivity owners reject incompatible units, kinds
or forms. The material state and pure phase must match the compiled card. The
temperature span must have continuous support for all three wall properties.
The source temperature coordinate cannot also occur among the fixed axes.

This mode requires a cylindrical bore. Station zero supplies the inlet; every
later station ends one well-mixed gas cell and supplies its initial temperature
for steady relaxation. All stations must have the same composition because
this example has one inlet and no gas injections. Downstream temperatures
are outputs of the enthalpy solve. They are not fitted to the initial guesses.

Each cell uses the exact annular-cylinder conductivity integral between its
inner material face and the prescribed outer reservoir. The explicit film
resistance connects bulk gas to the inner face; zero selects perfect contact.
Dividing the areal film resistance by the actual inner area gives its total
K/W resistance, so changing the cell grid does not change that physical film.
The example does not estimate this film from flow, infer isotropy from a scalar
source claim, or model wall heat storage, phase change or axial wall conduction.
It explicitly assigns the source scalar conductivity to the radial direction.
Refine the thermal cells separately from the DG mesh and acoustic wall memory.

The accepted field feeds the existing frozen, no-mean-flow linear acoustic
model. Prescribed through-flow belongs to steady field preparation; acoustic
convection and pressure drop are not included. The reported mean Mach limit
is an admission check, not an error bound for those omissions. The optional
viscothermal model retains its existing rigid, isothermal acoustic wall
assumption; this block changes the equilibrium gas state, not the wall's
frequency-dependent acoustic thermal impedance.

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
temperature/species during the pulse and molecular relaxation are outside this
DG example. The optional steady-wall block prepares the equilibrium field before
the pulse; it does not add a transient gas or wall solve to DG. The source
manifests linked from the reed example retain references and validity domains.
