# Material-driven electrothermal time evolution

The existing `electrothermal` command can advance a material-bound scene in time:

```sh
cargo run -p fs-cli --bin frankensim -- --json electrothermal scene.json \
  --transient time.json > trajectory.json
```

Without `--transient`, the existing steady command and result schema remain in
use. `scene.json` retains the schema described in `MATERIAL_ELECTROTHERMAL.md`:
canonical material packs, geometry, material assignments, complete source state,
thermal faces, electrode voltage/current conditions, contacts and solve budgets.
No new material catalog or material-name presets are introduced.

## Time request

An illustrative `time.json` is below. The declared temperature interval must be
inside the actual density and heat-capacity source support for every material;
these bounds do not authorize extrapolation or establish phase stability.

```json
{
  "schema": "frankensim.electrothermal.transient.v1",
  "time_step_s": 0.1,
  "steps": 8,
  "temperature_min_k": 290.0,
  "temperature_max_k": 350.0,
  "enthalpy_interpolation_error_j_kg": 0.00001,
  "energy_absolute_tolerance_j": 0.0000001,
  "energy_relative_tolerance": 0.000000001,
  "schedule": [[1, 1.0, 1.0], [5, 0.0, 1.0]]
}
```

A schedule row is `[first_step, electrode_drive_scale, authored_power_scale]`.
Step numbers start at **1**, are strictly increasing, and cannot exceed `steps`.
The first row must start at 1. Each row applies to the complete named step and
subsequent steps up to the next row. Thus this example powers steps 1--4 and
switches off the electrical drives for steps 5--8: four intervals of 0.1 s are
powered. There is no interpolation, hidden mid-step switching or state reset.

Drive scale multiplies every declared electrode voltage and current, preserving
the drive type. Authored-power scale multiplies only the scene's `power_w`; it
does not also multiply Joule heating. Joule heat comes from the electrical
solution. Thermal boundary conditions remain fixed. Pressure, composition and
other source coordinates retain the scene's explicitly declared values.

The request is limited to 10,000 steps, 16,777,216 cell-steps and 256 KiB. These
are resource bounds, not guarantees of runtime. Existing inner thermal, outer
electrical/thermal and Krylov iteration budgets still apply at each step.

## Shared material and solver path

Each material must additionally supply `density` [kg/m3] and
`specific_heat_capacity` [J/(kg K)] under the existing dimensional property
schemas. The existing `ResolvedHeatCapacityCurve` resolves the complete source
context over the requested interval, retains the selected Cp/density receipts,
and constructs its bounded single-phase enthalpy chart. The explicit scene role
chooses solid or fluid storage; it does not discover phase stability. Reference
enthalpy is zero at the lower endpoint; only changes enter the energy equation.
This is not a mixture with an independent absolute-enthalpy source.

Cell masses are resolved once from the **actual initial-state density query**
and voxel volume. They do not change as density data vary with temperature.
Re-evaluating `rho(T)*V` after each step on a fixed grid would create or destroy
mass; this command does not do that. Source density curves remain available to
the chart but do not silently move geometry or generate mechanical work.

For each time step the code solves the backward-Euler balance

```text
m_c [h_c(T_new) - h_c(T_old)] / dt
    + sum(outward conductive and admitted advective heat fluxes)
    = authored_power_c + joule_power_c(T_new).
```

The thermal iteration uses the chart's conservative positive secant
`C_c = m_c * [h(T_iter)-h(T_old)]/[T_iter-T_old]`, with a declared one-sided
chord slope when the temperatures coincide. `T_old` and mass are immutable
throughout every nonlinear iteration. Signed enthalpy increments integrate
crossed chart slopes, avoiding cancellation against large reference enthalpies.

`EnergyConfig::step_with_conductivity` lowers `C/dt*(T_new-T_old)` into the
existing cell-diagonal exchange assembly and calls the same conservative
conductivity-field solver. Geometry, material IDs, physical thermal contacts,
spatial conductivity and boundary assembly are unchanged. Insulated domains
are anchored by storage, not by an invented fixed-temperature cell.

An optional fluid role uses this same path with its existing prescribed
volume-conservative face flow. Its source density and Cp must be constant over
the **complete** requested interval and match the original reference claims;
all source breakpoints are checked, not just endpoints. Inlet and backflow
temperatures must lie inside that interval. This keeps the native `rho*Cp*Q*T`
flux consistent with `m*delta(h)` storage, including the cancellation of the
arbitrary enthalpy reference. Neighboring solids may retain variable Cp.
See `MATERIAL_ELECTROTHERMAL_FLUID.md` for admission and physical boundaries.

The steady and transient command paths share the material/electrical/thermal
coupling driver. The transient path additionally checks actual chart storage
against the thermal linearization using the **sum of absolute cell defects**,
so opposing local errors cannot cancel. Every proposed endpoint must still be
supported by the material sources. The outer loop checks temperature and local
Joule-source agreement; the complete accepted step audits actual `m*delta(h)`
against endpoint electrical heating and physical heat outflow.

## Output and accounting

The transient result uses `frankensim.electrothermal.transient.result.v1` and
contains final temperatures, potentials, local Joule heat, final resolved
conductivities, invariant cell masses, work count, material/chart identities,
initial mass-source identities and final local material-state identities.
Fluid scenes also retain the original `fluid_reference` transport state.
`records` contains each accepted endpoint time, peak temperature, signed stored
energy change, deposited source energy, physical outward heat and energy defect.
A failed run emits a diagnostic and nonzero exit, with no partial trajectory on
stdout; neither scene nor material files are modified. Memory stores the current
spatial state and per-step summaries, not a full spatial field at every step.

The low-level step API returns both an `EnergySolution` and a `TransientRecord`.
The solution's balance describes the **storage-augmented linear system**, so its
sink total includes the storage diagonal. The record separately reports physical
outflow; storage must not be counted again as heat lost to an environment. The
command independently recomputes chart enthalpy changes and endpoint Joule heat
rather than promoting that linearization record to exact nonlinear authority.

The per-step final energy tolerance is `absolute_j + relative * scale`, with
scale equal to the sum of absolute cell enthalpy changes plus absolute source
and physical-outflow energies. This is a numerical consistency tolerance, not a
bound on temporal discretization, source error or experimental prediction error.

## Verification status and boundaries

Twelve new Rust regressions were authored across the storage interface (4),
material-chart operations (3), and file-driven transient path (5). They cover
insulated heating, cooling/time refinement, heterogeneous contact/capacity,
signed and reference-independent enthalpy, drive switch-off, data-only density
swaps, variable Cp, source limits, cancellation and unchanged-input replay.

These Rust tests have **not** been compiled or executed in the implementation
environment: Cargo/rustc are unavailable. Focused selectors, under the normal
repository execution policy, are:

```sh
cargo test -p fs-lbm --lib conjugate::energy::material_field::tests
cargo test -p fs-material --lib heat_capacity_state::tests
cargo test -p fs-cli --bin frankensim electrothermal_command::
```

Independent Python matrix/nonlinear-network calculations were executed. The
synthetic insulated three-cell fixture deposits 1.8 J over four powered steps.
After four further unpowered steps its reference centre temperatures are
300.2458631128 K for 2-kg cells and 300.1352040740 K for 4-kg cells. With the
fixture's variable Cp, conductivity and resistivity, the independently computed
centre is 300.2440891879 K and supplied heat is 1.79060906094 J. These are
reference calculations, **not outputs from the unexecuted Rust program**.
The standalone cooling check gave refinement error ratios 1.95 and 1.97.

The model is fixed-geometry, fixed cell mass, isotropic ohmic conduction with
optional prescribed constant-density fluid heat transport. Cp is used in the
negligible-volume-change solid and incompressible-fluid thermal approximations;
there is no thermal-expansion work, stress evolution, phase change, latent heat,
chemistry, fluid momentum/pressure evolution, variable-density or variable-Cp
fluid advection, AC current, uncertainty propagation or experimental qualification.
Exterior gray radiation uses its existing explicitly admitted source-backed
owner. The declared single phase is a caller assumption,
not a phase-boundary discovery. Generated-chart interpolation and source error
remain separate, and tighter chart tolerance can exhaust the existing knot cap.
Finite time-step accuracy requires its own refinement study. No arbitrary-law
Picard convergence, automatic timestep reduction or checkpoint/resume is claimed.
