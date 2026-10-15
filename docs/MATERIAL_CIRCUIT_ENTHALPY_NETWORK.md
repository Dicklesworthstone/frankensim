# Finite circuits and material-bound thermal networks

`fs_couple::circuit_enthalpy` connects the existing `MaterialCircuit` to the
shared `fs_thermal::enthalpy::DepositionSimulation`. Each resistor has an
explicit receiving thermal cell; capacitors and inductors retain their original
charge and flux histories. Thermal cells may conduct to one another, exchange
heat with authored reservoirs, contain independent sources, and cross a supplied
latent-heat plateau. Passive thermal cells need no electrical branch.

## One physical interval, existing owners

The unknowns are changes in **cell energy in joules**, ordered by `CellId`, not
circuit branch order. Every trial starts from the same accepted circuit and
thermal histories. The material callback sees the complete tentative thermal
endpoint, resolves every resistor and selected thermal-link coefficient, and
returns original material records. `MaterialCircuit::step_candidate` advances
the original Kirchhoff descriptor equations and supplies branch-local heat.
The thermal owner's `BackwardEulerStep` supplies the local residual:

`actual cell energy change - circuit branch heat - dt * other endpoint heat rates`.

Neither a trial circuit evaluation nor a numerical Jacobian probe consumes the
accepted capacitor charge or deposits heat into a live cell. Only an admitted
joint candidate returns both histories. Every resistor needs exactly one active
heat destination; missing, duplicate or foreign branch bindings refuse. The
generic driver can map several branches to one receiving cell when the caller
justifies that homogenization. The concrete material adapter instead requires
one whole-cell uniform conductor per resistor, with distinct receiving cells.

The final fresh evaluation must pass the maximum **local** cell residual,
global thermal residual, original circuit storage/supply and resistor-power
checks, deposited-heat check, and combined circuit/thermal energy balance.
Internal Joule heat is subtracted from the thermal external-source account
before the combined audit. Independent authored or callback heat is not lost
or counted as electrical supply. All evaluations, including the final one,
consume the explicit circuit-trial budget. Cancellation, unsupported material
queries, invalid finite-difference stencils and nonconvergence return neither
physical endpoint.

## Data-only material bindings

`circuit_enthalpy::material::MaterialNetwork` connects supplied caloric charts,
material cards, electrical geometry and thermal-link geometry. No handwritten
source-resolution callback is needed for this uniform-cell model:

```rust,ignore
use fs_couple::circuit_enthalpy::{CircuitEnthalpySolve, ResistorCell};
use fs_couple::circuit_enthalpy::material::{
    MaterialCell, CellResistor, ConductionLink, MaterialNetwork,
};
use fs_couple::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
use fs_thermal::enthalpy::{CellId, MaterialId};

// card_a/card_b, curve_a/curve_b, contexts, exact conductivity keys and claim
// selections are caller-supplied material data. The thermal timeline must have
// been initialized using THESE SAME borrowed curves and material-table IDs.
let cells = [
    MaterialCell {
        cell: CellId(30), material_id: MaterialId(0), curve: &curve_a,
        point: ThermalMaterialPoint {
            material: &card_a, context: &context_a, temperature_axis: "T",
            liquid_fraction_axis: Some("melt"),
        },
        conductivity: Some((&conductivity_a, thermal_selection_a)),
    },
    MaterialCell {
        cell: CellId(10), material_id: MaterialId(1), curve: &curve_b,
        point: ThermalMaterialPoint {
            material: &card_b, context: &context_b, temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: Some((&conductivity_b, thermal_selection_b)),
    },
];
let resistors = [
    CellResistor {
        destination: ResistorCell { branch: 1, cell: CellId(30) },
        selection: electrical_selection_a, length_m: length_a, area_m2: area_a,
    },
    CellResistor {
        destination: ResistorCell { branch: 2, cell: CellId(10) },
        selection: electrical_selection_b, length_m: length_b, area_m2: area_b,
    },
];
let links = [ConductionLink {
    link: 7, distance_a_m: half_distance_a, distance_b_m: half_distance_b,
}];
let binding = MaterialNetwork { cells: &cells, resistors: &resistors, links: &links };
// Empty external ports means a source-free circuit; energy then comes from
// its finite stored charge/flux. Target is an absolute time in seconds.
let proposal = binding.step_candidate(
    &thermal, target_s, &circuit, &[], controls, &initial_guess_j, || false,
)?;
let (next_thermal, next_circuit, report) = proposal.into_parts();
// Publish BOTH endpoints together. All returned branch heat is ALREADY included.
```

The binding covers every active cell and thermal link. The original timeline's
borrowed caloric object is checked, not merely a matching temperature sample or
numeric material ID. This in-process address check survives timeline cloning;
it is not a portable source identity or scientific validation of the chart.
The phase chart and electrical/thermal property queries must refer to the same
material card. Their source receipts, claim selections, dimensional/static
quantity keys and all undeclared physical coordinates remain intact.

Each conductor uses the existing `MaterialOhmicBody` source and geometry checks:
`length * current_area` must match its whole receiving cell volume. Thermal
transfer area comes from the existing link, not the electrical cross-section.
For the uniform two-layer approximation, the link conductance is

`G = link.area / (distance_a / k_a + distance_b / k_b)`.

Both conductivities are resolved at their actual endpoint T/phase. Explicit
zero conductivity or area gives isolation; missing linked conductivity refuses,
even if the other side is insulating. Only scalar static conductivity and
source curves on declared temperature/fraction coordinates are admitted. The
caloric chart may change phase fraction at fixed temperature, and a declared
fraction-dependent resistivity or conductivity responds without inventing a
Kelvin increment. No implicit alloy mixing rule or generic metal preset exists.

For additional repeatable source laws, call the generic `step_candidate` with a
callback that starts from `binding.resolve(...)` and adds only independent watts
to `sources_w`. Original authored boundary sources remain in place. Do not add
Joule power yourself or debit external physical state inside that callback.

## Shared numerical bridge and model limits

The Newton--Krylov owner still supplies iteration, Krylov actions and
globalization. `fs_solver::finite_difference::FallibleFiniteDifference` bridges
fallible physical residuals to its numerical interface, retaining the first
original source error and poisoning further callbacks instead of returning a
false zero residual. The spatial LBM implicit electrothermal path now uses this
same bridge; its current, contact and thermal residual equations are unchanged.
The difference scale is explicitly in joules, and unsupported probes refuse.

Thermal rates use **backward Euler**. Circuit currents use the existing midpoint
scheme with **endpoint material resistances**. The combined method is generally
**first order**, despite the circuit owner's isolated midpoint order. Conservation
and algebraic convergence do not establish time accuracy, uniqueness, entropy
production, unconditional stability, or measured-material fidelity. Refine the
time grid and numerical difference scale independently.

Geometry, circuit topology, capacitance/inductance and caloric charts stay fixed
during an interval. Process start-time thermal events first, split at events
strictly inside the interval, and explicitly process target-time events after
joint acceptance. Source discontinuities also need interval boundaries. This
is a bounded dense Kirchhoff circuit and a uniform-cell thermal network, not a
spatial AC/induction solver, moving melt, evolving C/L material law, switches/arcs,
convection calculation or simultaneous mechanical deformation integrator.
The separate single-body circuit/gas/piston adapter remains available; this
network path does not implicitly attach gas or mechanical state to each cell.

## Focused native checks

```sh
cargo test -p fs-solver --lib finite_difference -- --test-threads=1
cargo test --no-fail-fast -p fs-couple --test circuit_enthalpy --test circuit_material_network -- --test-threads=1
cargo test -p fs-lbm --features implicit-enthalpy-coupling --test joule_enthalpy_implicit -- --test-threads=1
```

The existing shared-thermal workflow selects these tests. Local Rust execution
was unavailable during authoring. Independent Python equations checked finite
supply depletion, exact two-body heat partition/conduction, endpoint feedback,
conductivity-only substitutions, latent-energy changes at identical Kelvin,
and first-order refinement. Those references are not native Rust test passes.
