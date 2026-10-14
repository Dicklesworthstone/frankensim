# Heat-driven material rod

The `heated_material_rod` executable connects existing owners rather than
introducing a heat or mechanics solver:

`material card -> Cp/enthalpy and conductivity -> thermal cell history ->
RodMaterialSource -> current E/G/rho and integrated expansion -> rod equilibrium`.

A heater supplies the first segment for the first half of the experiment, then
switches off. Heat continues to move between segments and escape through the
external film. Temperatures are solved states, not an animation or prescribed
ramp. The thermal timeline and mechanical snapshot are published together only
after every source and solver accepts the candidate.

## Run

From the workspace root, these two commands differ ONLY in synthetic material
heat-capacity data (J/(kg K)); all geometry, supports, power and timing match:

```sh
cargo run --release -p fs-couple --example heated_material_rod -- --demo 500 free 1 4 20
cargo run --release -p fs-couple --example heated_material_rod -- --demo 1000 free 1 4 20
```

Use `fixed` instead of `free` to restrain axial extension at the far end and
observe thermal compression. Both modes use explicit straight axial guides:
transverse positions and frames are fixed at all nodes. This is not an
unsupported assertion that a freely bending member cannot buckle.

A canonical `NormalizedMaterialCardPack` uses the existing decoder:

```sh
cargo run --release -p fs-couple --example heated_material_rod -- --pack material.fsmat fixed 1 4 20
```

The final arguments are support, heater watts, duration seconds, and number of
coupling/observation intervals. This bounded example fixes length at 0.04 m,
radius at 0.002 m, four segments, external film coefficient at 10 W/(m2 K),
reference/initial/ambient temperature at 300 K and the caloric span at 300..500 K.
It prints these choices and the source identity. It admits T-only source
context under the legacy `T` axis, not unprovided pressure/moisture/frequency
coordinates. Missing conditions refuse; they are never populated by a preset.

The card needs density, specific_heat_capacity, thermal_conductivity,
young_modulus, shear_modulus and linear_thermal_expansion_coefficient. The
mechanical and caloric properties use their canonical dimensional schemas;
conductivity also accepts the existing static thermal-conductivity kind.
Cp/density must support the whole caloric span. Mechanical properties and
conductivity are checked at each used state. Claims are pinned instead of
switching datasets as temperature changes. The demo is explicitly synthetic,
not a named metal or measured material.

CSV reports every segment temperature, extension, axial force, elastic energy,
stored heat change, external heat, heat-balance defect, and the first segment's
applied conductivity. The original caloric source and applied conductivity
receipts are retained. Cell masses come from the mechanical reference source;
heating does not recompute mass from current density on frozen geometry.

## Reuse and limits

`fs_solid::rod::material_field::thermal::source::RodMaterialSource` is the
reusable library entry point. It supports different cards and temperatures in
every concentric layer and spatial segment, retaining the reference source and
all current section data. `resolve_at` takes temperatures flattened in
segment/layer order. The executable is its smaller single-material example.

Heat uses the existing conservative explicit enthalpy network with CFL-limited
substeps (maximum 0.01 s). Adjacent conductance is the series resistance of two
half cells. k(T) is refreshed at coupling boundaries and held through that
interval, with endpoint source checks. This is a lagged coefficient scheme;
refine coupling intervals AND thermal timesteps before claiming accuracy.
The external-film Biot check is necessary, not a full spatial-error certificate.
Uniform cell/cross-section temperature, fixed thermal geometry and explicit
external-film data remain model assumptions.

Mechanics is quasistatic and one-way. Elastic energy is NOT counted as heat,
but its change and constraint work are not fed back into the thermal balance.
There is no total thermomechanical conservation, dynamics, branch tracking,
melting, yielding, creep, thermal radial strain, remeshing or experimental
qualification claim. A changed mechanical domain can reject an otherwise
valid heat step without altering the accepted timeline or pose.

Five library and four executable regressions are authored:

```sh
rch exec -- cargo test -p fs-solid --test rod_material_source
rch exec -- cargo test -p fs-couple --example heated_material_rod
```

They were not compiled or executed in the implementation environment, which
has no Rust/RCH toolchain. Independent four-cell heat equations and series
axial compliance calculations are reference checks, not native test results.
