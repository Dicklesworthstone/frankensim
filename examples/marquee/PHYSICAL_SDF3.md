# Physical-size 3-D topology studies

The feature-gated `sdf3-study` command accepts `physical-curved-height-sdf`
domains in SI coordinates. Unlike the retained `curved-height-sdf` unit-cube
example, this explicitly selected domain uses the declared lower and upper
bounds in geometry integration, elasticity, load integration, density filters,
refinement and exported node positions. It does not rescale loads or modulus.

For lower corner `(x0,y0,z0)` and upper corner `(x1,y1,z1)`, material occupies

```text
z < z0 + height-m + curvature-per-m * (x - x0) * (x1 - x)
```

`height-m` is thickness above the lower z plane, not a global elevation.
`curvature-per-m` has inverse-length units. A translated design changes its
physical coordinates without changing its dimensions. Clamp choices are
`left` (x0), `right` (x1), `front` (y0), `back` (y1), and `bottom` (z0).
Each is a homogeneous three-component displacement constraint on that box
plane. `top` is not admitted because the graph is strictly below z1.

For example, replace the domain/physics/scenario values in
`bracket-3d-adaptive.fsim` with:

```lisp
(domain
  :type physical-curved-height-sdf
  :bounds ((0.0 0.0 0.0) (0.1 0.05 0.04))
  :height-m 0.028
  :curvature-per-m 1.0)
```

Declare the actual modulus in pascals, actual reference body-force densities
in N/m³, the chosen clamp, and a filter radius such as `0.01` metres. The
example above defines a 100 x 50 mm footprint, not a one-metre bracket with
scaled display labels. Each new mesh reintegrates the same physical load law.

```bash
cargo run -p fs-cli --features sdf3-study -- --json study \
  YOUR_PHYSICAL_STUDY.fsim /tmp/physical-3d.db
```

Admission currently bounds each box span to `[1e-6,1000]` metres, aspect ratio
to 64, `height/Lz` to `[0.1,0.8]`, and `curvature*Lx²/Lz` to `[0,0.4]`.
The filter radius lies between `0.001*min(span)` and `max(span)`. Coordinate
magnitudes must not exceed `1e8` times their corresponding span, so an offset
cannot collapse distinct finest-grid nodes. Existing geometry, work and
memory caps still apply. These are numerical admission limits, not a
conditioning, peak-memory, manufacturing or physical-validation guarantee.

The legacy domain retains its original unit-cube restrictions and field
arithmetic. The physical model participates in the canonical source identity
and existing stage checkpoint, resume and retained report/package paths.
Same-executable replay remains required. Neither producer moves the implicit
boundary or certifies continuum accuracy.

Native regressions in `sdf3/geometry_tests.rs` exercise graph enclosures,
translated clamps, volume scaling, the `s^5` compliance and `s^2` displacement
scaling under constant body-force density, and independent modulus/load
changes. They must be run natively; these identities do not substitute for
execution of the Rust tests.
