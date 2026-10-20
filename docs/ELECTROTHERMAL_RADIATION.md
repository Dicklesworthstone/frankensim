# Source-backed radiation in electrothermal scenes

The existing `electrothermal scene.json` and all-solid
`electrothermal scene.json --transient time.json` paths accept an optional
`radiation` object. Emissivity, electrical resistivity, thermal conductivity
and caloric storage use the original material cards. No second radiation law,
heat solver or temperature history is introduced.

Add this fragment to a complete `frankensim.electrothermal.v1` scene:

```json
{
  "radiation": {
    "max_iterations": 250,
    "temperature_tolerance_k": 1e-9,
    "power_tolerance_w": 1e-9,
    "relative_tolerance": 1e-10,
    "relaxation": 0.5,
    "surfaces": [{
      "name": "exterior",
      "material_index": 0,
      "faces": [[1, 3]],
      "surroundings_k": 300,
      "initial_surface_temperature_k": 450,
      "convection": {"h_w_m2_k": 2, "ambient_k": 280}
    }]
  }
}
```

Each face is `[linear_cell_index, local_face]`, with x varying fastest and
local faces ordered x-, x+, y-, y+, z-, z+. It must be an exterior solid face
of the declared material, unique across all radiation patches. Its existing
`thermal_faces` rule must be `adiabatic`: that delegates boundary exchange to
this explicit model, rather than installing a second owner on a prescribed
flux, temperature or convection face. Omit `convection` for radiation alone.
Convection here is an explicitly supplied film law, not inferred gas dynamics.

The same material pack must supply exactly one static, dimensionless
`hemispherical-total-emissivity` claim in `(0,1]`. This retains the original
surface owner's range and source requirements. Queries use the surface
temperature coordinate named by the material binding and retain all its
pressure, process and quantity descriptors. A source is not selected by metal
name, and a visible-color sample does not become total emissivity. The initial
surface temperature is a numerical guess, not an independent stored body;
it and every trial/final surface state must have source support.

## Surface temperature and one conductive path

The existing `AmbientRadiationPatch` supplies the stable gray-radiation secant

```text
h_rad = epsilon(Ts) * sigma * (Ts + Tr) * (Ts^2 + Tr^2)
```

If convection is present, `h = h_rad + h_conv` and the effective reservoir is
`(h_rad*Tr + h_conv*Tair)/h`. The shared `PreparedBoundaryExchange` lowers this
through the solid half-cell resistance, `G = A / (w/(2*k) + 1/h)`. It recovers
surface temperature from the conductive drop, not by identifying Ts with the
cell-center temperature. Conductivity remains the resolved piecewise-constant
cell value; the adapter does not add subcell nonlinear diffusion.

At every thermal trial, material emissivity is re-resolved at the recovered Ts.
Both raw surface-temperature change and applied-versus-nonlinear heat must pass
their declared tolerances. A sum of absolute face defects also prevents opposite
local errors from cancelling. The existing electrical feedback and thermal
contact/conductivity/storage checks remain active. The inner solve budget is
the smaller of `max_thermal_passes` and radiation `max_iterations`.

Radiation and convection share ONE conductive half cell, but their signed heat
rates remain separate. Positive heat leaves the solid; a hotter reservoir can
heat it. Raw energy reports already contain the generated boundary sinks.
Transient accounting adds actual nonlinear surface outflow alongside original
physical sinks, never alongside those same generated sinks a second time.
The result contains per-face Ts, emissivity, radiative/convective/applied watts,
heat mismatch and the exact final material-state identity. Transient records
also report radiated watts and physical outflow joules. Source or convergence
refusal publishes no partial field/trajectory and does not change input files.

## Scope and focused checks

This is fixed-geometry, opaque gray diffuse radiation to a black isothermal
reservoir occupying view factor one for each declared face. It does not compute
occlusion, enclosure exchange, participating-gas radiation, spectral optics,
coating evolution or phase change. Prescribed-fluid steady scenes may radiate
from their exterior SOLID faces; the fluid-transient restriction is unchanged.
Source nominal values, algebraic convergence and energy closure do not establish
experimental accuracy or time/space discretization error.

Five shared-owner and four command regressions cover signed half-cell exchange,
storage, source conditions, material-only emissivity swaps, nonlinear emissivity
with a distinct convective reservoir, surface-versus-bulk source support,
duplicate/occupied faces, replay and transient cooling. Native Rust compilation
and tests were not run in the authoring environment (Cargo/rustc unavailable).

Independent NumPy/SciPy equations were executed: at 450 W, unit area and unit
conductivity, epsilon 0.2 gives Ts=467.53201285 K and Tcell=692.53201285 K;
epsilon 0.8 gives Ts=366.38575097 K and Tcell=591.38575097 K. Changing only the
emissivity data changes equilibrium, not the authored electrical drive. A
three-cell cooling reference loses 41.71330905 J versus 66.68712641 J in five
0.1-s steps; those match its caloric storage decreases, with worst per-step
reference defects below 6e-11 J. These are synthetic equation checks, not
outputs from an executed Rust binary or measured materials.

```sh
cargo test -p fs-lbm --lib conjugate::boundary_exchange::tests
cargo test -p fs-conduction --lib radiation::equilibrium::state::condition_tests
cargo test -p fs-cli --bin frankensim electrothermal_command::radiation_tests
```
