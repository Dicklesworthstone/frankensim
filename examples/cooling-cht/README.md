# cooling-cht: voxel conjugate heat transfer from a JSON scene

`frankensim [--json] cooling-cht <scene.json>` solves steady laminar airflow
by finite-volume SIMPLEC on a voxel grid (forced convection, or natural and
mixed convection with the Boussinesq force when `gravity_m_s2` is declared)
and one conservative energy equation over fluid and solid cells. It reports
per-material and per-source temperatures, the hottest solid voxel, flow and
energy residuals, and the energy balance.

Scene (`frankensim.cooling-cht.v1`):

- `size_m`, `voxel_m`: the box domain; each size must be a whole number of
  voxels. Voxel centres decide occupancy: a box `[min_m, max_m)` takes the
  voxels whose centres it contains, and an edge lying exactly on a voxel
  centre takes the voxel at its min edge and leaves the one at its max edge
  (deterministically, whatever the decimal rounding).
- `fluid`: `"dry-air-300k"` (default) or explicit `density_kg_m3`,
  `specific_heat_j_kg_k`, `conductivity_w_m_k`, `kinematic_viscosity_m2_s`.
- `materials`: `name`, `conductivity_w_m_k` (a number, or `[k_x, k_y, k_z]`
  for grid-aligned orthotropic solids such as PCB laminates), optional
  `emissivity` (0 to 1; any non-zero value enables radiation) and
  `volumetric_heat_capacity_j_m3_k` (transients).
- Radiation: exposed faces of emissive solids exchange gray diffuse
  radiation with each other and with the surroundings seen through
  openings, inlets and fans (at their temperatures); domain walls and
  non-emitting solids reflect perfectly, so a sealed box radiates from its
  hot parts to its emissive walls (model enclosure walls as emissive solid
  boxes). Exchange factors come from deterministic Monte Carlo rays on face
  patches (optional `radiation`: `rays_per_face`, default 256, `seed`,
  `patch_size`, default 4 faces, and `surface_exchange`, default true;
  false keeps only the escape to the surroundings). The result reports
  `radiation.radiated_w` (to the surroundings) and the energy balance's
  `sink_outflow_w`. Transients refuse radiation.
- `contacts` (optional): `between` (two material names) and
  `resistance_m2_k_w`, a per-area interface resistance (thermal interface
  material, bonded or pressed joint) on every face the two materials share.
- `solids`: `material` plus either a box (`min_m`, `max_m`) or a closed STL
  mesh (`stl`, a path relative to the scene file; optional `scale` and
  `offset_m`, placing `world = scale * stl + offset`). Mesh occupancy is the
  robust generalized winding number above one half (exact solid-angle sum,
  or the dipole octree above 4096 triangles). Later solids override
  earlier ones.
- `sources`: `name`, `power_w`, `min_m`, `max_m`; the power is spread
  uniformly over the solid voxels whose centres the box covers.
- `faces` (`x-`, `x+`, `y-`, `y+`, `z-`, `z+`; missing faces are adiabatic
  walls): `inlet` (`velocity_m_s`, `temperature_k`), `opening` (`ambient_k`;
  pressure zero, flow in either direction), `symmetry`, or `wall`
  (adiabatic, or one of `temperature_k`, `heat_flux_w_m2` into the domain,
  `htc_w_m2_k` with `ambient_k`).
- `fan` faces: `curve` `[[flow_m3_s, pressure_pa], ...]` (2 to 8 points,
  pressure non-increasing) and `temperature_k`; the delivered flow is the
  operating point where the curve meets the system resistance, reported as
  `fan_flow_m3_s` / `fan_pressure_pa`.
- Optional `transient`: `time_step_s`, `steps`, `power_schedule`
  `[[time_s, scale], ...]` (piecewise linear, constant beyond the ends;
  default 1), `initial_temperature_k` (default: the inlet temperature).
  Backward Euler marches the energy equation over the converged steady
  forced flow; materials then need `volumetric_heat_capacity_j_m3_k`. The
  result adds per-step peak solid temperatures (up to ~200 records) and the
  worst per-step energy closure. With `"flow": "unsteady"` the flow marches
  with the energy from rest instead (implicit SIMPLEC, `"scheme"`: `"bdf2"`,
  the default, or `"backward-euler"`; `inner_iterations`,
  `inner_tolerance`; `inlet_schedule` `[[time_s, scale], ...]` scales the
  inlet velocities, e.g. a fan start-up ramp), buoyant when gravity is
  declared, so flows that never settle (shedding wakes, plumes) are marched
  instead of refused. The unsteady result reports the march: steps and
  sweeps, final kinetic energy and flows, per-step records, the peak and
  final solid temperatures, and per source the final and time-averaged
  (second half) maxima. With `"energy": "steady-on-mean-flow"` the march
  carries the flow only, and the steady energy equation (with any
  radiation) runs on its time-averaged fluxes (second half of the march):
  the practical answer for a wake that never settles while the solids'
  thermal time constants are minutes. That neglects the unsteady
  correlation `<u' T'>` and needs a forced flow (buoyant scenes refuse).
  Frozen-flow transients still refuse buoyant scenes.
- Optional `gravity_m_s2`, `expansion_per_k` (default `1 / T_ref`),
  `reference_temperature_k` (default: the first inlet or opening
  temperature), `solver.tolerance`, `solver.max_iterations`,
  `solver.turbulence` (`"laminar"`, the default, or `"lvel"`: the LVEL
  algebraic eddy viscosity and its turbulent conductivity, for
  transitional or turbulent fan-driven flow; the result then reports
  `flow.max_eddy_viscosity_ratio`), `limits.wall_seconds`.

- Optional `internal_fans`: `name`, `axis` (`"x"`, `"y"`, `"z"`), `at_m` (a
  voxel face plane strictly inside the domain), `direction` (`"+"` or
  `"-"`), `min_m`/`max_m` (the transverse extent; the entries along `axis`
  are ignored), and `curve`. The pressure rises across the plane by the
  curve's value at the flow through it; the result lists each fan's
  operating point under `flow.internal_fans`.
- Optional `resistances`: `{"type": "grille", axis, at_m, min_m, max_m,
  loss_coefficient | free_area_ratio}` (pressure drop `1/2 rho K |u| u`;
  a free-area ratio uses Idelchik's thin perforated plate) or
  `{"type": "porous", min_m, max_m, permeability_m2, inertial_per_m}`
  (Darcy-Forchheimer, scalar or per axis; a missing permeability means no
  viscous term). Grilles on a vent sit one voxel inside the open face.
- A solid may be a parametric plate-fin heatsink instead of a box or STL:
  `{"material": ..., "heatsink": {"base_min_m", "base_size_m", "fin_count",
  "fin_thickness_m", "fin_height_m", "fins_along": "x" | "y"}}` (fins stand
  on the base top, the outer ones flush with its edges). A fin that covers
  no voxel centre, or a gap that keeps no fluid voxel, refuses: the grid
  cannot represent that design.
- Optional `study`: `{"parameters": [{"name", "path": [keys and indices
  into this scene], "values": [...]}, ...], "objective": {"minimize":
  quantity}, "constraints": [{"quantity", "min", "max"}]}` evaluates every
  combination of the values (at most 64) as an ordinary scene under its
  own wall budget and ranks the completed, feasible variants. Quantities:
  `max_solid_temperature_k`, `source:<name>`, `component:<name>`,
  `internal_fan:<name>`, `fan_flow_m3_s`, `inflow_m3_s`. Refused variants
  are listed with their refusal and never ranked; the result
  (`frankensim.cooling-cht.study.v1`) lists every evaluation and the best.
  Variants are solved on `parallelism` threads (default: the available
  cores); the report is identical for any thread count. A grid search: no
  optimality claim between grid points.
- Instead of `size_m` + `voxel_m`, a graded grid: `"grid": {"x": [{"to_m",
  "voxel_m"}, ...], "y": [...], "z": [...]}`, each zone a whole number of
  uniform cells from the previous zone's end (0 first), so fine voxels go
  only where features and boundary layers need them. Box, plane and
  source positions use the actual cell centres and faces; the result
  reports `graded` and `voxel_m` as the smallest width.
- Optional `components`: JEDEC two-resistor compact models (`name`,
  `min_m`/`max_m`, `board_side`, `power_w`, `junction_to_case_k_w`,
  `junction_to_board_k_w`). The box blocks flow; the junction reaches the
  case top and the board only through the two resistors (sides adiabatic),
  and the result reports `components[].junction_temperature_k`, `case_w`
  and `board_w`. Steady scenes only.

Examples:

- `heatsink-duct.json`: the ducted plate-fin heatsink of
  `crates/fs-lbm/examples/heatsink_cht.rs` at 1 mm voxels (2 W chip,
  0.25 m/s air at 300 K).
- `vented-heatsink-natural.json`: a vertical plate-fin heatsink (three 2 mm
  fins) in a 30 x 20 x 60 mm column open at the bottom and top, cooled by
  natural convection only, at 2 mm voxels. Measured (debug build, 315 s):
  181 energy couplings, 905 SIMPLEC sweeps, induced draft 2.23e-5 m^3/s
  (peak 0.15 m/s), 0.5 W chip at 332.21 K (64 K/W), all heat leaving by
  advection through the top opening (balance 2e-12). This scene declares
  no emissivity; adding one enables radiation through the openings.
- `fan-enclosure.json`: a 100 x 60 x 30 mm electronics enclosure at 2.5 mm
  voxels: a 50 % perforated vent grille behind the open x- face, an axial
  fan (40 Pa shut-off, 6 l/s free delivery) in an ABS baffle, an orthotropic
  FR4 board with a 3 W and a 2 W package, a porous card array, an open
  exhaust, and LVEL turbulence. Measured (debug build, 529 s): 420 SIMPLEC
  iterations, fan operating point 3.65e-3 m^3/s at 19.6 Pa (on its curve),
  peak eddy viscosity 86 x molecular, packages at 394.3 K and 372.5 K,
  energy balance 9e-14.
- `stl-heatsink-duct.json`: the Journey A body `../heatsink-fan/heatsink.stl`
  (80 x 60 mm base, four 6 mm fins) in a 100 mm duct along its fin channels
  at 2 mm voxels, with a 3 W, 20 x 20 mm die under the middle fins, at
  0.1 m/s (duct Reynolds number about 250). At 0.5 m/s (about 1250) the
  flow behind the blunt body does not settle (measured at 4 mm voxels: the
  SIMPLEC residuals cycle around 1e-2) and the run refuses as not steady or
  diverged: a steady laminar answer is not claimed there.
  Measured at 0.1 m/s (debug build, 262 s): 82 SIMPLEC iterations, junction
  323.43 K for the 3 W die, energy balance 6e-11. At 2 mm voxels the
  voxelized body holds 4.80e-5 m^3 of the STL's 5.28e-5 m^3 (its 5 mm base
  top lies on a voxel centre and falls outside under the tie rule): refine
  `voxel_m` before reading temperatures to better than that geometry.

Results are Estimated numerical evidence at one resolution: turbulence only
through the algebraic LVEL closure (its friction runs 13-16 % above
turbulent channel correlations; see `crates/fs-lbm/CONTRACT.md`), radiation
only to the surroundings, no temperature-dependent properties; staircase
geometry; not a ledger-backed `.fsim` run. Refine `voxel_m` to measure resolution
sensitivity.
