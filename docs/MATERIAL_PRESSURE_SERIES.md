# Pressure-derived material flow and convection

`fs_airflow::graph::thermal::medium::series::hydraulic::pressure::SeriesPressureFlow`
solves one fixed-state series path from its differential pressure and the
existing `MaterialFlowResistance` sections. Each section retains its original
typed material owner. The solver samples each owner once and uses its density
and viscosity in the original loss law; no material name selects a model.

For conserved mass flow `m`, the native section law reduces to
`dp_total = A*m + B*m*abs(m)`, with `A = sum(a*mu/rho)` and `B = sum(b/rho)`.
The inverse uses scaled single-term solutions instead of subtracting similar
quadratic roots or squaring dimensional coefficients. It then re-evaluates
EVERY original section to check Reynolds support, pressure closure, tangent
and dissipated power. Local volume flows may differ when densities differ.
Total hydraulic heat is `sum(dp_i*Q_i)`, not total pressure times an unrelated
reference volume flow. Zero and reverse drive retain original admission.

```rust,ignore
use fs_airflow::graph::thermal::medium::series::hydraulic::pressure::{
    PressureSection, SeriesPressureBudget, SeriesPressureFlow,
};
use fs_airflow::graph::thermal::medium::convection::ForcedConvectionSettings;

let flow = SeriesPressureFlow::solve(
    vec![PressureSection { resistance: original_geometry, medium: original_fluid }],
    differential_pressure_pa,
    SeriesPressureBudget {
        maximum_sections: 16,
        absolute_pressure_tolerance_pa: 1e-8,
        relative_pressure_tolerance: 1e-12,
    },
    &gate,
)?;
let transfer = flow.forced_convection(
    0, selected_correlation, characteristic_length,
    ForcedConvectionSettings::default(),
)?;
let boundary = transfer.convection().robin_boundary(reference_temperature)?;
```

`forced_convection` obtains speed from the solved local volume flow and the
ORIGINAL hydraulic flow area. It uses the exact retained hydraulic transport
snapshot for Reynolds, Prandtl and conductivity through the existing convection
owner. It performs no second query of the original material. The result borrows
the complete pressure owner and exposes its original medium beside the native
convection result; a detached numeric tuple gains no material authority.
The original correlation, evidence and Robin-boundary conversion remain intact.
Heat geometry and correlation choice are explicit, with their original support
checks. No zero-speed regularizer or substitute correlation is introduced.

This closes a fixed-state pressure -> flow -> convection connection. It does
NOT implement temperature-feedback circulation, prescribe a new film state,
deposit hydraulic work, or refresh any existing simulation history. The caller
must account for the reported hydraulic dissipation and justify the frozen
section/material and heat-geometry assumptions. Pump curves, shaft efficiency,
transient momentum, compressible pressure work, cavitation and changing phase
are outside this model. Differing densities do not by themselves qualify a
compressible or variable-volume trajectory. The tangent is at fixed materials,
not the derivative of a temperature-coupled operating point.

Work and storage scale linearly with section count; there is no root-search
loop or additional global matrix. This is an algorithmic statement, not a
measured hardware speedup. Explicit section and pressure-closure budgets and
cancellation checks remain active. Thermal evaluation has the original
constant-size correlation cost and borrows source data.

## Focused verification

```sh
rch exec -- cargo test -p fs-airflow --lib graph::thermal::medium::series::hydraulic::pressure
```

Seven native regression groups are authored: independent heterogeneous pressure
and power equations, source-only substitutions, split-section and reverse-flow
invariance, one original query, extreme inverse scaling, regime/budget refusal,
and shared pressure/convection data with independent flat-plate arithmetic.
Compilation, native execution and rustfmt are UNEXECUTED in this editor because
Rust, Cargo, rustfmt, RCH and DSR are unavailable.

An independent 100-digit Decimal reference checked 300 scaled inversions,
with maximum relative discrepancy 3.44e-16. Separate manufactured convection
calculations verify pressure/velocity/Re/Pr/h propagation and signed-flow
symmetry. These are equation checks, not execution of Rust, experimental
material qualification, or a full-workspace or real-time result.
