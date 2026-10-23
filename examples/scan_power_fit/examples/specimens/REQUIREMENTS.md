# Physical constraints in CSV scan fitting

The optional `--requirements FILE` supplies hard inequalities to the same native
SQP problem that fits the CSV readings. These are not objective penalties or
post-fit clipping. The original thermal model, source normalization, measurement
loss, exposure/sensor interpretation and checkpointed adjoints remain in use.
Commands without the option retain their previous nonnegative-power constraints.

From the repository root:

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example spatial_measurement_fit -- \
  examples/scan_power_fit/examples/specimens/multiple_measurements.csv \
  --specimen examples/scan_power_fit/examples/specimens/rectangular.scan \
  --heaters examples/scan_power_fit/examples/specimens/multiple.heaters \
  --requirements examples/scan_power_fit/examples/specimens/multiple.requirements \
  --substeps 4
```

This uses the existing 4 kg rectangular manufactured specimen, two fitted held
stages, a fitted triangular pulse, and a fixed 2 W background heater. The desired
8/16/24 W schedule absorbs 56 J in total. The supplied 48 J cap includes the fixed
heater's 8 J; the analytical constrained reference is 8/16/16 W, with final
enthalpy 162 rather than 164 J/kg and liquid fraction 0.62 rather than 0.64. The
quadratic measurement objective is 4, not zero. These are manufactured discrete
balance values, not experimental validation or a native execution record.

## File grammar and units

The first noncomment line is `FS_SCAN_REQUIREMENTS_V1`. Blank lines and `#`
comments are ignored. Each subsequent line has exactly one form:

```text
power SOURCE KNOT LOWER_W_OR_* UPPER_W_OR_* SCALE_W
energy START_S END_S MAXIMUM_J SCALE_J
state TIME_S VERTEX OBSERVABLE LOWER_OR_* UPPER_OR_* SCALE
peak START_S END_S OBSERVABLE MAXIMUM SCALE SMOOTHING VERTICES
slew SOURCE MAXIMUM_RISE_W_S MAXIMUM_FALL_W_S SCALE_W_S
speed SOURCE MAXIMUM_M_S SCALE_M_S
```

Sources and knots use the zero-based heater-file identities, not optimizer
coordinate numbers. `power` requires a selected/fitted knot; a fixed value does
not become adjustable just because a requirement names it. The two sides must
be strictly ordered when both are given. Use `*` for an omitted side; at least
one side is required. Scales must be finite and positive.

`energy` integrates the authored absorbed-power schedules of **all** sources,
including fixed heaters and unselected values, over `0 <= START_S < END_S <=`
the observation horizon. It is absorbed joules, not electrical input energy,
net storage change or boundary heat. The existing exact schedule integral is
used; thermal stepping retains its original discrete quadrature.

`state` observes an actual mesh vertex at one exact positive physical time.
`OBSERVABLE` is `temperature_k`, `specific_enthalpy_j_kg` or
`liquid_mass_fraction`; limits and scale use that coordinate's units. Fraction
state limits lie in [0,1], not percentages. These constraints are not applied to
a time mean or lagged sensor output, even when the fitted readings are.

`peak` applies the existing smooth upper envelope to every selected vertex at
every physical-grid endpoint in the inclusive window. `VERTICES` is `all` or a
strictly increasing semicolon-separated list such as `0;3;7`. A single-time
window is permitted. Its start must be positive. `SMOOTHING` uses the observable
units and must be positive. The report separately exposes the sampled maximum
and positive smoothing margin. The envelope maximum is not clamped to [0,1]
for fractions: smoothing can legitimately raise it above the physical maximum.
This is not a continuous-time or unsampled-space peak certificate.

Constraint state/peak times enter the source/measurement event union **before**
substep refinement and exposure/sensor weight construction. Every refined peak
sample is counted. Adding times may change the underlying discrete solution;
it does not merely add labels to an already computed trajectory. Energy windows
need no extra thermal stop because their source integral is analytic.

`slew` covers every continuous power segment, including a ramp to a fixed
shutdown. `speed` covers the simultaneous 3D velocity on the complete declared
path, including heater-off travel and return motion. Native constant/ramp/linear
schedules are supported; held paths/powers are refused, not interpreted through
finite-difference chords. The heater-file loader uses one interpolation policy
per source, so the example applies these two declarations only to source 1.
Paths are fixed in this executable: it cannot repair an overspeed path by
changing heater powers. Neither speed nor slew establishes acceleration, jerk,
collision avoidance or hardware feasibility beyond the specified inequality.

## Results, tolerances and work

The executable prints raw measurement predictions/residuals independently from
physical constraint values. Within each requirement family, report order is
file declaration order; slew segments also follow each source's original
segments. `native_constraint_row` reports signed scaled inequalities (<= 0 is
feasible) in the native order: nonnegative powers, lower/upper selected powers,
energy, slew rise/fall, squared speed, lower/upper states, peak envelopes.

At the fixed feasibility tolerance 1e-8, an ordinary row permits at most
`1e-8 * SCALE` violation in its physical units. Speed instead uses
`(|v|^2 - maximum^2) / scale^2`, so its tolerance is in squared-speed units.
Numerical solver error and peak smoothing are distinct from this tolerance.
The final physical replay must reproduce every accepted inequality bit-for-bit;
a changed or missing constraint row is an error.

A successful exit requires convergence and feasibility, **not** zero measurement
loss. No convergence/global-optimum or experimental-validity claim follows from
merely parsing the file or exporting the sources. Invalid declarations, absent
IDs, unsupported signals, physical solver/derivative failures and cancellation
propagate; they never silently disable requirements or become penalty scores.

The command admits at most 64 KiB of requirement text, 64 declarations, dense
SQP dimension 256 (controls plus inequalities), 65,536 retained peak samples,
1,048,576 response-gradient entries and 1,048,576 objective/response transpose
solves per evaluation. These are explicit shape/work limits, not total memory
or wall-time guarantees. Original solver and cumulative evaluation limits remain.

## Focused checks

```bash
cargo test --manifest-path examples/scan_power_fit/Cargo.toml --lib requirements
cargo test --manifest-path examples/scan_power_fit/Cargo.toml \
  --example spatial_measurement_fit requirement_tests
```

The first checks lowering, exact clocks, identity, budgets and refusal. The
second exercises the real physical fit under binding energy, melt-fraction and
shutdown-rate limits; insertion of requirement times alongside exposure/sensor
readings; exact fixed-heater export; and unsuccessful input handling. These tests
were authored but not executed in the delivery environment: DSR, RCH, Cargo,
rustc and rustfmt could not start because their executables were absent.
