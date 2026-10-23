# Fit measurements on an explicitly supplied physical specimen

`spatial_measurement_fit --specimen FILE` uses supplied tetrahedral geometry,
reference mass density, caloric curve, initial enthalpy, conductivity, exterior
heat loss, source positions/footprint and switching times. It does not select a
unit cube or infer material/boundary data from the measured values. Without the
flag, the historical manufactured reference specimen is unchanged.

Run from the repository root:

```sh
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example spatial_measurement_fit -- \
  examples/scan_power_fit/examples/specimens/rectangular_measurements.csv 0 0 \
  --specimen examples/scan_power_fit/examples/specimens/rectangular.scan --substeps 4
```

This manufactured fixture has volume 2 m^3, reference density 2 kg/m^3 and mass
4 kg. Starting at 150 J/kg, absorbed powers 8 W over 0..1.5 s and 20 W over
1.5..4 s give 153 and 165.5 J/kg. It also includes a 1.5..4 s exposure mean,
a lagged sensor reading, and a liquid-mass-fraction measurement. These numbers
are analytical fixture balances, not a claim that the Rust example was run.

## Versioned, explicit SI input

The file starts with `FS_SCAN_SPECIMEN_V1`. Lines are whitespace-separated;
`#` starts a comment. Every scalar/tuple declaration below is required exactly
once. Unsupported/duplicate declarations refuse instead of acquiring defaults.

```text
conductivity_w_m_k 2
reference_density_kg_m3 2
initial_specific_enthalpy_j_kg 150
switch_time_s 1.5
end_time_s 4
first_center_m 1 0.5 0.5
second_center_m 1 0.5 0.5
half_width_m 3 2 2
exterior adiabatic
```

Alternatively use `exterior robin HTC_W_m2_K AMBIENT_K`, with both values
positive. This is applied to the entire extracted exterior. For example,
`exterior robin 0.1 300` removes 10 W from this fixture while every node remains
at its 310 K melting plateau. The same mean measurements then require 18 W and
30 W, not 8 W and 20 W. Treating that specimen as adiabatic would infer the wrong
powers even though it could fit the observations.

Provide two or more `phase_knot H_J_KG T_K LIQUID_MASS_FRACTION` rows in increasing
enthalpy order. Fraction is dimensionless, not percent. All knots use the one
explicit reference density. The original equilibrium-phase constructor checks
its monotonicity/domain rules; the loader checks that initial enthalpy is inside
the supplied range. Conductivity is a declared positive isotropic constant,
not a database-calibrated material or inferred temperature-dependent law.

Supply `vertex_m X Y Z` rows, followed by or interleaved with `tet I J K L`
connectivity. Indices are zero-based vertex-row numbers. Tet rows define the
indices used by `tet` and `volume` readouts in the measurement CSV; vertex and
weighted readouts refer to the same supplied vertices. No remeshing, nearest-node
snapping or unit conversion occurs. The native mesh constructor owns geometry
and boundary extraction. Repeated vertices in a tet, duplicate connectivity,
orphan vertices, out-of-range indices, nonmanifold faces and degenerate elements
refuse. Geometric intersection of otherwise unrelated tets is not checked; the
input must be a valid conforming mesh.

The example admits at most 2 MiB, 4,096 vertices, 16,384 tetrahedra and 256 phase
knots. Existing measurement limits remain 128 KiB, 256 rows, 4,096 refined
endpoints and 65,536 temporal contributions. These are input/retained-shape
bounds, not a total-memory guarantee. The declared clock must satisfy
`0 < switch_time_s < end_time_s`; measurement times must be positive and not
exceed shutdown. Exact source-switch/shutdown and measurement knots are inserted
before `--substeps 1..128` refinement. Exposure and sensor windows still start
at accepted endpoints, not at the initial state.

## What is fitted and what is retained

Exactly two nonnegative, held, domain-absorbed powers are fitted. Their positions
are the two supplied centres, their widths are the supplied half widths, and
the final power at `end_time_s` remains exactly zero. These are not incident
optical or electrical powers. Position/power jumps remain explicit held signals;
this mode makes no travel-speed or slew-rate claim. There is no cooling-only
continuation beyond shutdown in this consumer.

All candidate evaluations and the final physical replay borrow the same mesh,
material, chart and boundary, using the existing P1 enthalpy solver, source
normalization, spatial/temporal measurement operators, Huber/quadratic losses,
checkpointed adjoints and constrained SQP. Solver and cumulative evaluation
budgets remain bounded. Temperature/phase corners and failed physics remain
errors rather than invented gradients or penalties. Output reports actual
mesh counts, reference mass, source times, predictions and raw residuals.

This is a bounded standalone example input format, not a replacement for native
FrankenScript. This rung is homogeneous, fixed-geometry and fixed-density, with
one heater and a uniform adiabatic/Robin exterior. It does not add heterogeneous
material regions, interfaces, radiation, deformation, material inference,
general multi-heater scheduling, identifiability or experimental validation.
A converged fit is not a guarantee that its supplied physical model is correct.

## Focused checks

```sh
cargo test --manifest-path examples/scan_power_fit/Cargo.toml --lib specimen
cargo test --manifest-path examples/scan_power_fit/Cargo.toml --example spatial_measurement_fit
```

Three loader regressions and four executable regressions were authored for this
increment. They cover input admission, exact shutdown, different geometry/counts,
mass/clock recovery, Robin-loss recovery, temporal readouts and legacy defaults.
Native compilation/execution was unavailable in the authoring environment; no
native test pass or measured performance is claimed.
