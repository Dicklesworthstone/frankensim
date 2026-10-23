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

Without `--heaters`, exactly two nonnegative, held, domain-absorbed powers are fitted. Their positions
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
the declared heaters and a uniform adiabatic/Robin exterior. It does not add heterogeneous
material regions, interfaces, radiation, deformation, material inference,
identifiability or experimental validation.
A converged fit is not a guarantee that its supplied physical model is correct.

## Multiple heaters and arbitrary selected power stages

`--heaters FILE` replaces the legacy heater with a bounded program of native
sources. Each source has its own held or linear position/power table and fixed
footprint. Any nonterminal power knot can be selected for fitting, including an
initial zero. Other powers remain bit-for-bit fixed. A completely fixed heater
still contributes to every physical solve; it is not dropped from the model.

Run the three-source manufactured fixture from the repository root:

```sh
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example spatial_measurement_fit -- \
  examples/scan_power_fit/examples/specimens/multiple_measurements.csv \
  --specimen examples/scan_power_fit/examples/specimens/rectangular.scan \
  --heaters examples/scan_power_fit/examples/specimens/multiple.heaters --substeps 4
```

The file starts with `FS_SCAN_HEATERS_V1`, followed by these records:

```text
heater ID hold|linear HALF_WIDTH_X_M HALF_WIDTH_Y_M HALF_WIDTH_Z_M
knot ID TIME_S X_M Y_M Z_M POWER_W fixed
knot ID TIME_S X_M Y_M Z_M INITIAL_POWER_W fit REFERENCE_W
```

Heater IDs are contiguous from zero in declaration order. Knots may interleave
different sources but must be strictly time-ordered within each source. Every
heater starts at time zero, has at least two knots, and ends with an explicitly
fixed zero-power value no later than the specimen's `end_time_s`. Thereafter its
clamped power is zero. A heater that starts later must explicitly declare zero
before starting; pre-first-knot behavior is never inferred. `hold` and `linear`
apply to both positions and powers. Widths/positions use the specimen's world
metre frame. Powers are absorbed in the domain, not incident optical watts.

`REFERENCE_W` is a positive optimizer-coordinate scale: `power = reference*x`,
with `x >= 0`. It is not a maximum power. Fixed knots have no reference scale.
Selected controls are ordered by heater ID and then knot index, independently
of the order in which source records interleave. The two positional starting
powers cannot be combined with `--heaters`; the file supplies all starting
values. Invalid or duplicate options refuse rather than silently taking priority.

The V1 specimen's legacy `switch_time_s`, centre and width fields are still
required and validated by that format, but are unused with `--heaters`. Its
`end_time_s` is now the observation horizon and may be later than every heater's
shutdown, retaining a cooling-only tail. All declared heater event times and
measurement/window boundaries enter the physical grid before refinement. Close
events are never merged by a tolerance; unrepresentable substeps refuse.

The executable admits 16 heaters, 512 total knots, 32 selected powers, 4,096
authored position/power/width values and 128 KiB of heater input. It derives the
source-storage and dense SQP dimension allowances from the full program. The
existing maximum 127 fitting evaluations plus one final replay is unchanged.
Failure to converge is an unsuccessful exit, not a reason to omit a heater.

Output with `--heaters` lists every source, knot, interpolation, time, fitted
power and selection flag, including fixed background heaters and shutdowns.
Export checks the original clocks, positions, widths and unselected powers;
there is no clipping or snapping to zero. Legacy commands keep their previous
two-column power output. Spatial, exposure and sensor-lag measurements and
Huber/quadratic scoring use their original operators and adjoints.

For the supplied 4 kg insulated specimen, heater 0 fits 8 W over 0..1 s and
16 W over 1..2 s, heater 1 fits a triangular pulse peaking at 24 W at 3 s, and
heater 2 supplies a fixed 2 W until 4 s. Energy balance gives the declared means
152.5, 157 and 164 J/kg at 1, 2 and 4 s. Falsely omitting the fixed heater still
fits these measurements, but with 10, 18 and 28 W: fitting quality alone cannot
establish model correctness. These are analytical fixture values, not reported
native Rust results. Colocated heaters with indistinguishable temporal inputs
can be unidentifiable; this workflow does not assert uniqueness.

Paths, clocks, interpolation and footprints remain fixed design inputs. This
flag does not add speed, slew, thermal-state or energy-cap declarations to the
CSV executable; those remain separate existing library interfaces. In
particular, a held path is allowed to jump and carries no finite-speed claim.
Linear ramps still use the original time-discrete thermal solver; fitting on a
coarse grid is not a continuum-accuracy guarantee.

## Focused checks

```sh
cargo test --manifest-path examples/scan_power_fit/Cargo.toml --lib specimen
cargo test --manifest-path examples/scan_power_fit/Cargo.toml --lib heaters
cargo test --manifest-path examples/scan_power_fit/Cargo.toml --example spatial_measurement_fit
```

Three loader regressions and four executable regressions were authored for this
increment. They cover input admission, exact shutdown, different geometry/counts,
mass/clock recovery, Robin-loss recovery, temporal readouts and legacy defaults.
Native compilation/execution was unavailable in the authoring environment; no
native test pass or measured performance is claimed.

The multi-heater extension adds four loader tests and three executable tests
for asynchronous clocks, fixed/background power, full exports, input caps,
three-control inverse fitting, event/window alignment and option admission.
These tests were authored but could not be compiled or executed in the
authoring environment; DSR, RCH, Cargo and rustc were unavailable.
