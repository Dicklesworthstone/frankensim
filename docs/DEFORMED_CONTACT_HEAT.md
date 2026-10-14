# Deformed facing surfaces to conservative thermal contact

The explicit deposition lane can now keep each material's own covered surface
measure while exchanging one equal-and-opposite heat flow. This connects the
facing-surface query to heat transfer without pretending that differently
stretched roads have the same reference contact area.

## Facing geometry

`fs_solid::bridge::contact::pair` supplies `Query`, `Budget`, `Row` and
`Footprint`. Call `lower.neighbour_contact(&upper, query, budget, cancelled)`
on two accepted `Bridge` states at the same time. It merges the two world-x
partitions, inverts each surface's own axial/rotation map, and checks the
quadratic relative gap, including interior extrema and capture roots.

Rows retain both source-cell IDs, their reference-coordinate intervals,
minimum/maximum gaps, two captured material-reference areas and a common
captured projected area. Nonmatching meshes and different thermal strains
are allowed. No overlap is distinct from overlap with zero captured area.
Excessive penetration or exhausted work refuses the whole query, not a
truncated footprint. Shared x/z frame directions, translated origins, centred
widths and single-valued facing top/bottom surfaces are required. These are
nominal floating-point geometry results, not interval clearance certificates.

## Heat exchange and side coverage

`fs_conduction::transient::deposition::ContactPatch` consists of a `ThermalLink`
and `covered_area_m2: [area_on_a, area_on_b]`. Each area must use the same
measure as that incident cell's `surface_area_m2` and ambient coefficient.
The two reference measures can differ even when they describe the same world
footprint. The link's nominal area can be that footprint's projected area;
it does not silently replace either covered area or rescale conductance.

For an explicitly supplied conductance per projected area, a footprint row
maps to:

```rust
ContactPatch {
    link: ThermalLink {
        id,
        a: lower_thermal_cell,
        b: upper_thermal_cell,
        area_m2: row.projected_area_m2,
        conductance_w_k: coefficient_w_m2_k * row.projected_area_m2,
    },
    covered_area_m2: [row.lower_reference_area_m2, row.upper_reference_area_m2],
}
```

A caller using another conductance law must supply its total W/K value
explicitly. Reference-area boundary coefficients are not interchangeable with
current-area coefficients without a corresponding Jacobian. Width changes,
curved-surface area and deformation-dependent convection are not supplied here.

`simulation.contact_patch(id)` reads the complete current interface.
`replace_contact_patches(&updates, checkpoint)` admits a batch of
`ContactPatchUpdate { id, expected, replacement }`. `None` means separated.
Remove all changed interfaces before adding replacements when the common
partition changes. The method does this internally, checks expected link AND
side areas, and returns a candidate only after the complete batch passes.
Due scheduled events must be processed first. This operation changes neither
mass, enthalpy, time nor boundary coefficients.

Each side's covered area enters its actual exposed-surface heat loss and CFL
bound. The same total conductance supplies opposite internal powers, including
when the two areas differ. Insulating patches still cover both surfaces.
Separation clears coverage; a subsequent ordinary `Event::Contact` is symmetric
again. Existing symmetric events and `replace_contacts` keep their original
meaning. A symmetric update to an asymmetric current patch refuses rather than
silently dropping its second measure, even when the exchange link is unchanged.
The separate implicit parcel lane is not modified by this API.

## Coupled example

```sh
cargo test -p fs-solid bridge::contact::pair
cargo test -p fs-conduction transient::deposition
cargo test -p fs-couple --example deformed_contact_cooling
cargo run -p fs-couple --example deformed_contact_cooling
```

The example uses three lower cells and five upper cells, each with actual
thermomechanical birth references. Both roads have free axial contraction and
zero transverse load. They start flat and touching at different temperatures;
heat exchange expands one road and contracts the other. Each mechanical endpoint
recomputes both solved facing surfaces and replaces the thermal common partition
for the next interval. The two side areas stay distinct, and source-cell identity,
original mass and enthalpy survive every remap. No fake surface adjustment,
reservoir heat, or temperature reset compensates for a mismatched area.

The conducting case is compared with a geometrically identical insulating-law
case. The numerical fixture is synthetic. Flat faces remain in geometric contact
without a solved pressure or adhesion model; capture alone does not establish a
bond. The example does not generalize that assumption to arbitrary deformable
contact. It uses staggered endpoint coupling: thermal and coupling timesteps
need independent refinement, and a remap is not within-step event localization.

Six facing-geometry tests, seven thermal-coverage tests and two example tests
are authored. Native Rust compilation and test execution remain unverified in
the authoring environment: DSR, its documented absolute fallback, RCH and Cargo
are unavailable. Executed Python/SciPy references exercise geometry, incidence
heat equations and the synthetic coupled limit, not the compiled Rust program.
No calibrated printer, pressure, healing, friction, impact or whole-print
prediction follows from this implementation alone.
