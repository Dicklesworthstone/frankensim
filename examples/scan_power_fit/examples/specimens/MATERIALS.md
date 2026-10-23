# Heterogeneous material inputs for CSV scan fitting

`spatial_measurement_fit --materials FILE` uses distinct element-assigned
conductivities, reference densities and equilibrium phase charts, with explicit
initial enthalpy for each material and optional nodal overrides. It uses the
original heterogeneous enthalpy producer, its source normalization and its
adjoints; the input loader does not implement another heat solver.

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example spatial_measurement_fit -- \
  examples/scan_power_fit/examples/specimens/two_materials_measurements.csv \
  --specimen examples/scan_power_fit/examples/specimens/two_materials.scan \
  --materials examples/scan_power_fit/examples/specimens/two_materials.materials \
  --heaters examples/scan_power_fit/examples/specimens/two_materials.heaters \
  --substeps 2
```

The manufactured fixture contains two *disconnected* tetrahedral bodies of
volume 1/6 m3 each. Their masses are 1 kg and 2 kg, and their latent intervals
are 100..200 J/kg at 310 K and 100..300 J/kg at 330 K. Separate heaters fit
2/5 W and 6/10 W, respectively. These are analytic fixture expectations, not
experimental data. The CSV includes separate region enthalpies, exposure and
lagged readings, and liquid mass fractions evaluated through each body's own
phase chart. A volume-averaged mass fraction is not a bulk-mass-weighted fraction.

Append `--population 12 5000 --seed 17 --objective-target 1e-8 --workers 4`
to use the existing bounded value-only population engine. Its target is not
guaranteed within the requested budget. Hard requirements, robust loss, exact
shutdowns, shared evaluation accounting and final physical replay are unchanged.
No native runtime or speedup result is claimed by this document.

## Format and units

```text
FS_SCAN_MATERIALS_V1
material 10 2 6 150
phase_knot 10 -200 280 0
phase_knot 10 100 310 0
phase_knot 10 200 310 1
phase_knot 10 900 380 1
element 0 10
initial 1 160
```

`material ID K RHO INITIAL_H` declares constant isotropic conductivity in
W/(m K), frozen reference density in kg/m3, and initial specific enthalpy in
J/kg. Every material needs its own ordered `phase_knot ID H T FRACTION` rows:
J/kg, absolute kelvin, and liquid mass fraction in [0,1], not percent.
`element TET_INDEX ID` uses the specimen's zero-based original tetrahedron order.
Every tetrahedron must be assigned once; unknown, duplicate and unused material
IDs refuse. IDs need not be contiguous. Declare a material before its knots.

Optional `initial VERTEX_INDEX H` overrides that node's initial specific
enthalpy, within its assigned chart. Other nodes retain their material's
explicit initial value. A thermal history is never reset to a single global
initial value when candidates are rebuilt. No initial temperature is inverted
on a latent plateau and no phase composition is inferred.

The file replaces **all** homogeneous conductivity, density, caloric and initial
enthalpy values from the specimen for this fit. The existing specimen format
still requires those fields for its standalone meaning; its mesh, uniform
exterior boundary, observation horizon and source defaults remain operative.
Without `--materials`, existing commands retain their original homogeneous
physics. Supplying a bad material file never falls back to that physics.

The executable admits at most 128 KiB, 64 materials and 4096 total phase knots.
Its existing mesh limits bound assignment and initial-field sizes. Original
material, solver, energy, source and cancellation refusals remain authoritative.
Output includes the actual material count and total reference mass, not the
specimen's superseded homogeneous density times volume.

## Interface and modeling boundaries

Distinct caloric materials **cannot share a vertex**. The native producer
requires separate enthalpy traces, even when two declared curves have equal
numbers. This command currently has no contact-file input: separated bodies
are supported, but touching split traces requiring contact are refused by the
native no-undeclared-contact check. No adiabatic gap, perfect contact, material
mixture or weld is inferred. Native library callers can combine the immutable
material map with separately declared native contact operators.

This is fixed-geometry, fixed-reference-mass thermal modeling, with declared
constant isotropic conductivity per material. It does not calibrate properties,
model expansion/motion, create a conforming mesh, supply material-card provenance,
or prove continuum accuracy or identifiability. Multiple separated bodies in
one fit are not a coupled composite contact model.
