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
  for grid-aligned orthotropic solids such as PCB laminates).
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
  worst per-step energy closure. Buoyant scenes refuse a transient: their
  flow depends on temperature.
- Optional `gravity_m_s2`, `expansion_per_k` (default `1 / T_ref`),
  `reference_temperature_k` (default: the first inlet or opening
  temperature), `solver.tolerance`, `solver.max_iterations`,
  `limits.wall_seconds`.

Examples:

- `heatsink-duct.json`: the ducted plate-fin heatsink of
  `crates/fs-lbm/examples/heatsink_cht.rs` at 1 mm voxels (2 W chip,
  0.25 m/s air at 300 K).
- `vented-heatsink-natural.json`: a vertical plate-fin heatsink (three 2 mm
  fins) in a 30 x 20 x 60 mm column open at the bottom and top, cooled by
  natural convection only, at 2 mm voxels. Measured (debug build, 315 s):
  181 energy couplings, 905 SIMPLEC sweeps, induced draft 2.23e-5 m^3/s
  (peak 0.15 m/s), 0.5 W chip at 332.21 K (64 K/W), all heat leaving by
  advection through the top opening (balance 2e-12). No radiation is
  modelled, which matters for natural convection.
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

Results are Estimated numerical evidence at one resolution: no turbulence
model, radiation, or temperature-dependent properties; staircase geometry;
not a ledger-backed `.fsim` run. Refine `voxel_m` to measure resolution
sensitivity.
