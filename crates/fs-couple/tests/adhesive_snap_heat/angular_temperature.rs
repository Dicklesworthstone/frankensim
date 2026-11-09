//! Original temperature-dependent material cards and finite sensible receivers.
#[path = "../../../fs-contact/tests/adhesive_temperature/common.rs"]
mod material;

use fs_contact::normal_patch::adhesive::{JkrEnvelope, SourceJkrSphereContact};
use fs_couple::friction_heat::adhesive::snap::angular_temperature::JkrAngularTemperatureHeatControls;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_thermal::enthalpy::{CellId, DepositionSimulation};
pub(super) use material::Fixture;
use std::collections::BTreeMap;

const SENSIBLE_CAPACITY: f64 = 1e-3;

pub(super) fn material_fixture() -> Fixture {
    Fixture::with_options(10.0, 0.05, false, 320.0)
}

pub(super) fn normal(fixture: &Fixture) -> SourceJkrSphereContact {
    SourceJkrSphereContact::new(
        fixture.normal.clone(),
        0.01,
        1e-9,
        JkrEnvelope {
            maximum_indentation_m: 1e-5,
            ..fixture.jkr().envelope()
        },
    )
    .unwrap()
}

pub(super) fn curves(normal: &SourceJkrSphereContact) -> [EquilibriumEnthalpyPhaseCurve; 2] {
    let elastic = normal.normal_state().elastic();
    [elastic.surface_a_state(), elastic.surface_b_state()].map(|state| {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            state.resolved().card_identity(),
            SolidLiquidPhase::Solid,
            [280.0, 320.0]
                .map(|temperature_k| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: SENSIBLE_CAPACITY * (temperature_k - 300.0),
                    temperature_k,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 1000.0,
                })
                .to_vec(),
        )
        .unwrap()
    })
}

pub(super) fn warm_cells(thermal: &mut DepositionSimulation<'_>) {
    // Independently supplied finite heat changes the actual caloric cells at
    // the accepted time. The contact refresh may only observe this endpoint.
    let heat = BTreeMap::from([(CellId(0), 310.0), (CellId(1), 305.0)].map(|(id, t)| {
        let cell = &thermal.cells()[&id];
        (
            id,
            cell.mass_kg * (SENSIBLE_CAPACITY * (t - 300.0) - cell.specific_h_j_kg),
        )
    }));
    *thermal = thermal
        .apply_heat_j(&heat, 1e-22, 2e-22, || Ok(()))
        .unwrap()
        .0;
}

pub(super) fn controls() -> JkrAngularTemperatureHeatControls<'static> {
    JkrAngularTemperatureHeatControls {
        parameter_work_source: "prescribed material-temperature driver",
        maximum_absolute_parameter_work_j: 1e-6,
    }
}
