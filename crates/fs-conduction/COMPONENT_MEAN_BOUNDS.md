# Bound a component's mean without disconnecting it from the assembly

The `thermal-verification` feature now exposes
`verification::region::{solve_with_region_mean_bound, bound_temperature_region_mean}`.
Both use the original complete `ConductionProblem` and all original matching
`InterfaceSurface` declarations. Supply an explicit nonempty list of original
mesh-element indices for the component or region. An empty surface list means
no contacts are declared, not permission to ignore existing coincident faces.

```rust,ignore
use fs_conduction::verification::{FluxBudget, region::{
    GoalResidualLimits, RegionMeanConfig, bound_temperature_region_mean,
}};

// Existing cx, original problem, complete contact declarations, and temperature.
let config = RegionMeanConfig {
    dual: original_linear_config,
    residual_limits: GoalResidualLimits {
        max_rows: maximum_free_dofs,
        max_nonzeros: maximum_matrix_nonzeros,
    },
    flux: FluxBudget { max_cells: maximum_domain_cells, max_iterations: 128 },
};
let result = bound_temperature_region_mean(
    cx, problem, &contacts, &temperature, &component_cells, config,
)?;
let temperature_interval_k = result.bound.enclosure;
```

The supplied-field API runs **no primal solve**. `solve_with_region_mean_bound`
takes the same inputs except temperature, and an additional explicit primal
`SolveConfig`; it returns the original physical primal plus this analysis.
Neither function changes the input model, source, geometry or temperature array.

The regional integral dual is solved on the **whole domain**. Its source is one
on selected cells and zero elsewhere, assembled as `volume/4` at each selected
cell's four nodes. It is not a nodal indicator: that would smear the functional
into neighboring cells sharing a vertex. Matching contact resistances and
separate traces are retained. No insulated boundary is invented at the edge of
the selection. Other components can carry heat and dual response even when they
are not included in the temperature average.

The verifier uses the exact cell indicator, full-domain primal/dual majorants
and residual correction, and divides by the outward selected volume **last**.
Rounded assembly only proposes a dual; its error relative to that exact
functional remains in the bound. Finite inexact duals remain usable. The returned
`dual_analysis` is the original stored-system report in *integral* units, not a
replacement continuum certificate. No extra stability-proposal solves are run.

For two unit slabs of conductivity 1 with 400 K / 300 K exterior temperatures
and area-specific contact resistance 2, the analytical heat flux is 25 W/m².
The hot slab mean is **387.5 K**, the cold slab mean **312.5 K**, while their
combined mean is **350 K**. A whole-assembly average must not be relabelled as
the component average. A second fixture heats only the insulated left slab
with source `6*x`; it drains through the unselected right slab and has regional
mean **309.75 K**, not the combined **305.625 K**. These are analytical fixture
values, not claims that the Rust regressions passed in the authoring environment.

Invalid or duplicated cells, unselected invalid material, missing/conflicting
contacts, changed prescribed temperatures, budget exhaustion and cancellation
remain refusals. Selection order is canonicalized and the original cell IDs
are retained in the result. Selection is not a mesh identity certificate; the
caller must use the same original element numbering.

This is a nominal **linear**, tensor-conductivity, affine-source,
matching-contact, fixed-polyhedral-domain *volume-mean* bound. It is not a point
maximum, nonlinear/nonmatching-contact, CAD, material-uncertainty or physical
validation certificate. It does not add a `.fsim` report output or close the
broader q61wp.11 product-maximum requirement. All existing whole-domain APIs
and numerical kernels remain unchanged.

Focused checks registered in the existing thermal workflow:

```sh
cargo test -p fs-verify --features certified-speculation --test region_mean
cargo test -p fs-conduction --features thermal-verification --test region_mean_bound
```

Local native execution remains unverified: the authoring runtime has no
`rch`, `cargo` or `rustc`. Independent exact-rational and NumPy/SciPy fixture
checks are not Rust execution or interval-rounding verification.
