# Nonlinear sheet data through the existing mechanical solvers

A characterized stress/strain curve is a constitutive law, not an equivalent
constant modulus. This integration carries the explicitly piecewise-linear
`tabulated_green_strain/1` image from FrankenTextiles into FrankenSim material
claims, source-bound membrane equilibrium and prestressed vibration.

The coordinates are `[E11, E22, 2 E12]` Green strain and work-conjugate
second-Piola resultants `[S11, S22, S12]` in N/m. Nominal force/extension curves,
bulk stress in Pa and fiber names cannot be substituted. Energy is the analytic
integral of the three stress curves plus the explicit constant warp/weft
coupling term. Stress and tangent use those same curves. No scalar membrane
stiffness is added to the table response.

## Import an explicit nonlinear record

Choose an output directory that does not exist:

```sh
python3 scripts/import_franken_textiles.py \
  data/matdb/textile-examples/synthetic-hardening-sheet.json \
  synthetic-hardening-sheet /tmp/nonlinear-sheet \
  --nonlinear-membrane --allow-estimates \
  --citation 'Synthetic nonlinear sheet integration fixture' \
  --license 'Synthetic data; repository license applies' \
  --redistribution-terms 'Authored synthetic numerical data'

cargo run -p xtask -- matdb-pack \
  --manifest /tmp/nonlinear-sheet/manifest.tsv \
  --out /tmp/nonlinear-sheet/material.pack

cargo run -p fs-couple --example textile_sheet -- \
  /tmp/nonlinear-sheet/material.pack --membrane 300 1.02 0.05 9.81 0.0 0.1
```

The fixture is synthetic, not empirical textile calibration. The Python import
command is tested; native pack compilation and executable execution were not
available in the authoring host. The existing `textile_cavity` example also
resolves its walls through `resolve_sheet`, so this material selection reaches
its equilibrium/modal construction without a new cavity solver. Its time
evolution remains linearized around the nonlinear equilibrium, not a
finite-amplitude nonlinear cloth simulation.

The flag requires explicit `interpolation: "piecewise_linear"`. Missing/default
or `monotone_cubic` interpolation refuses rather than changing the source law.
Records with relaxation history refuse. Scalar membrane properties are excluded
by explicit table selection and remain in the retained original JSON;
explicitly requesting both profiles is an error. The original scalar/thermal
import behavior is unchanged without the flag.

## Shared material and numerical ownership

`ResolvedSheet::membrane_law` returns the complete selected law. The old
`membrane()` quadratic API refuses table-backed sources; it cannot return an
instantaneous tangent as a constant modulus. `try_from_material_with_card`
binds an existing complete selection to its exact card, preserving additional
thermal/contact properties and original claim receipts. Curves are compiled
once and shared immutably across triangles.

The original `MembranePatch` supplies reference geometry, grain, conserved areal
mass, Green strain and geometric stiffness. The original Newton and modal
solvers remain the numerical owners. Source strain intervals intersect knot
support, and neither trial evaluation nor accepted states extrapolate.

At table knots Newton uses the right-hand tangent, except at the upper endpoint
where it uses the left-hand tangent. Energy is C1, not C2 at a slope change.
Small-vibration extraction separately requires a unique two-sided local tangent:
it refuses slope-changing knots and source-domain boundaries. An eight-epsilon
strain neighborhood conservatively treats unresolved knot proximity as a
refusal. This numerical screen is not a certified material-accuracy bound.

The manufactured regression at 53.12 N expects 20 mm centre deflection and
7136 N/m tangent. Changing only the hardening data to a linear response raises
deflection to about 22.8804 mm at the same mass, geometry and supports. These
are mathematical references, not experimental data or Rust execution results.

## Focused checks and limits

```sh
python3 -m unittest scripts.tests.test_import_franken_textiles \
  scripts.tests.test_import_nonlinear_textile -v
cargo test -p fs-material --test sheet_tabulated
cargo test -p fs-solid --test membrane_tabulated
cargo test -p fs-solid --test membrane
```

All 15 Python importer tests passed locally, including nine unchanged upstream
regressions. Nine new Rust regression tests are authored but not compiled or
executed because Cargo/rustc are unavailable. Independent SciPy quadrature and
loaded-pyramid force/energy/tangent equations were checked separately; they do
not establish a native test pass.

Scope is fixed-state, reversible, separable nonlinear sheet elasticity with
constant cross coupling and strictly positive segment tangents. Coupling near
singularity is conservatively rejected with an eight-epsilon relative margin.
No monotone
cubic law, hysteresis, relaxation, plasticity, wrinkle relaxation,
bending/contact, heat/moisture evolution, full garment convergence or measured
material fidelity is implied. No runtime dependency is added.
