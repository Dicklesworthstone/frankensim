//! Original finite-gas transactions with independently calculated layer circuits.
use super::*;
use crate::sublimation::gas::{
    chamber::{GasChamberSolve, coating::GasSupportedChamberSolve},
    tests as fixture,
};
use fs_thermal::{
    enthalpy::DepositionError,
    sublimation::{SolidSurfaceState, SolverConfig, SublimationError},
};
use std::cell::Cell;

const DENSITY: f64 = 917.0;
const AREA: f64 = 1e-4;
const K: f64 = 1e-3;
const CONTACT: f64 = 2e-4;

struct Conductivity;
impl ConductivityLaw for Conductivity {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(K)
    }
}
fn controls() -> GasGrowingCoatingSolve {
    GasGrowingCoatingSolve {
        advance: GasResistiveAdvanceSolve {
            step: GasSupportedChamberSolve {
                chamber: GasChamberSolve {
                    surface: fixture::controls(),
                    maximum_surfaces: 8,
                },
                storage: SolverConfig {
                    relative_tolerance: 1e-12,
                    absolute_energy_tolerance_j: 1e-11,
                    max_iterations: 128,
                },
            },
            // Admit one coarse/fine comparison so an independent pair of
            // original physical intervals is an exact history reference.
            absolute_mass_tolerance_kg: 1.0,
            absolute_temperature_tolerance_k: 1000.0,
            absolute_specific_enthalpy_tolerance_j_kg: 1e7,
            absolute_heat_tolerance_j: 1e4,
            relative_tolerance: 1e-3,
            minimum_dt_s: 1e-10,
            maximum_dt_s: 0.02,
            maximum_trials: 10000,
        },
        maximum_deposit_volume_fraction: 0.01,
        absolute_extent_tolerance_m: 1.0,
        absolute_conductance_tolerance_w_k: 1e4,
    }
}
fn support(h: f64) -> Result<f64, SublimationError> {
    Ok(400.0 + h / 500.0)
}
fn state(chart: (f64, f64), mass: f64) -> ResistiveSurfaceState {
    ResistiveSurfaceState {
        substrate_mass_kg: 0.02,
        substrate_specific_enthalpy_j_kg: -50.0 * 500.0,
        deposit: SolidSurfaceState {
            mass_kg: mass,
            specific_enthalpy_j_kg: chart.0 - 50.0 * chart.1,
            // This is intentionally stale. The mass-owned geometry must bind
            // actual area before any molecular or reflected stream is taken.
            area_m2: 7.0,
        },
    }
}
fn planar_bound(mut row: GasResistiveSurface<'_>) -> GasResistiveSurface<'_> {
    row.state.deposit.area_m2 = AREA;
    let thickness = row.state.deposit.mass_kg / DENSITY / AREA;
    row.conductance_w_k = 1.0 / (thickness / (K * AREA) + CONTACT / AREA);
    row
}

#[test]
fn two_original_half_intervals_match_growth_refreshed_planar_circuit() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |h| fixture::temperature(h, chart);
    let initial = GasResistiveSurface {
        id: 4,
        law: &law,
        binding: &binding,
        state: state(chart, 2e-6),
        conductance_w_k: 1e9,
        substrate_temperature: &support,
        deposit_temperature: &deposit,
    };
    let geometry = ConductiveCoatingGeometry::planar(DENSITY, AREA).unwrap();
    let layers = BTreeMap::from([(
        4,
        GasGrowingCoating {
            geometry,
            conductivity: &Conductivity,
            contact_resistance_m2_k_w: CONTACT,
        },
    )]);
    let result = advance_growing_resistive_chamber(
        &[initial],
        &layers,
        None,
        None,
        &gas,
        &catalog,
        0.02,
        controls(),
        || false,
    )
    .unwrap();

    // Derive thickness and G directly from the mass balance/Fourier circuit,
    // then invoke only the original heat/molecular owner for each half step.
    let first = super::super::advance_resistive_chamber(
        &[planar_bound(initial)],
        None,
        &gas,
        &catalog,
        0.01,
        controls().advance.step,
        || false,
    )
    .unwrap();
    let mut middle = initial;
    middle.state = first.surfaces()[0].state;
    let second = super::super::advance_resistive_chamber(
        &[planar_bound(middle)],
        None,
        first.gas(),
        &catalog,
        0.01,
        controls().advance.step,
        || false,
    )
    .unwrap();
    let final_state = result.history().endpoint().surfaces()[0].state;
    let reference = second.surfaces()[0].state;
    fixture::near(
        final_state.deposit.mass_kg,
        reference.deposit.mass_kg,
        1e-16,
    );
    fixture::near(
        final_state.deposit.specific_enthalpy_j_kg,
        reference.deposit.specific_enthalpy_j_kg,
        1e-6,
    );
    fixture::near(
        final_state.substrate_specific_enthalpy_j_kg,
        reference.substrate_specific_enthalpy_j_kg,
        1e-8,
    );
    fixture::near(
        result.gas().state().state().temperature,
        second.gas().state().state().temperature,
        1e-9,
    );
    assert_eq!(result.history().audit().chamber.accepted_intervals, 1);
    let geometry = result.geometry()[0];
    assert!(geometry.final_geometry.front.mass_kg > geometry.initial.front.mass_kg);
    assert!(geometry.final_geometry.front.extent_m > geometry.initial.front.extent_m);
    assert!(geometry.final_geometry.conductance_w_k < geometry.initial.conductance_w_k);
    fixture::near(
        geometry.final_geometry.front.extent_m,
        final_state.deposit.mass_kg / (DENSITY * AREA),
        1e-18,
    );
    fixture::near(
        geometry.normal_recession_m,
        (initial.state.deposit.mass_kg - final_state.deposit.mass_kg) / (DENSITY * AREA),
        1e-18,
    );
    fixture::near(
        result.history().audit().thermal[0].interface_heat_j,
        first.surfaces()[0].interface_heat_j + second.surfaces()[0].interface_heat_j,
        1e-9,
    );
    fixture::near(
        result.history().audit().chamber.energy_residual_j,
        0.0,
        1e-7,
    );
    fixture::near(
        result.history().audit().chamber.mass_residual_kg,
        0.0,
        1e-14,
    );
}

#[test]
fn curved_growth_updates_outer_area_and_preserves_inner_contact_area() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |h| fixture::temperature(h, chart);
    let core: f64 = 0.001;
    let pi = core::f64::consts::PI;
    for (id, geometry) in [
        (
            3,
            ConductiveCoatingGeometry::spherical(DENSITY, core).unwrap(),
        ),
        (
            8,
            ConductiveCoatingGeometry::cylindrical(DENSITY, core, 0.004).unwrap(),
        ),
    ] {
        let initial = GasResistiveSurface {
            id,
            law: &law,
            binding: &binding,
            state: state(chart, 1e-6),
            conductance_w_k: 0.0,
            substrate_temperature: &support,
            deposit_temperature: &deposit,
        };
        let layers = BTreeMap::from([(
            id,
            GasGrowingCoating {
                geometry,
                conductivity: &Conductivity,
                contact_resistance_m2_k_w: CONTACT,
            },
        )]);
        let result = advance_growing_resistive_chamber(
            &[initial],
            &layers,
            None,
            None,
            &gas,
            &catalog,
            0.02,
            controls(),
            || false,
        )
        .unwrap();
        let history = result.geometry()[0];
        let mass = result.history().endpoint().surfaces()[0]
            .state
            .deposit
            .mass_kg;
        let (outer, area, contact) = if id == 3 {
            let outer = (core.powi(3) + 3.0 * mass / (4.0 * pi * DENSITY)).cbrt();
            (outer, 4.0 * pi * outer * outer, 4.0 * pi * core * core)
        } else {
            let outer = (core * core + mass / (pi * 0.004 * DENSITY)).sqrt();
            (outer, 2.0 * pi * outer * 0.004, 2.0 * pi * core * 0.004)
        };
        assert!(mass > initial.state.deposit.mass_kg);
        assert!(history.final_geometry.front.area_m2 > history.initial.front.area_m2);
        assert_eq!(
            history.initial.contact_area_m2,
            history.final_geometry.contact_area_m2
        );
        fixture::near(history.final_geometry.front.area_m2, area, 1e-18);
        fixture::near(history.final_geometry.front.extent_m, outer - core, 1e-17);
        fixture::near(history.final_geometry.contact_area_m2, contact, 1e-18);
        fixture::near(
            history.final_geometry.contact_resistance_k_w,
            CONTACT / contact,
            1e-10,
        );
        assert_eq!(
            result.history().endpoint().surfaces()[0]
                .exchange
                .solid
                .area_m2,
            history.final_geometry.front.area_m2
        );
    }
}

#[test]
fn original_conductivity_domain_failure_remains_typed_and_transactional() {
    struct LimitedConductivity {
        maximum_h: f64,
        queries: Cell<usize>,
    }
    impl ConductivityLaw for LimitedConductivity {
        fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
            self.queries.set(self.queries.get() + 1);
            if h > self.maximum_h {
                Err(DepositionError::MaterialDomain)
            } else {
                Ok(K)
            }
        }
    }
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let before_gas = gas.clone();
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |h| fixture::temperature(h, chart);
    let initial_state = state(chart, 2e-6);
    let conductivity = LimitedConductivity {
        maximum_h: initial_state.deposit.specific_enthalpy_j_kg + 1.0,
        queries: Cell::new(0),
    };
    let initial = GasResistiveSurface {
        id: 4,
        law: &law,
        binding: &binding,
        state: initial_state,
        conductance_w_k: 0.0,
        substrate_temperature: &support,
        deposit_temperature: &deposit,
    };
    let layers = BTreeMap::from([(
        4,
        GasGrowingCoating {
            geometry: ConductiveCoatingGeometry::planar(DENSITY, AREA).unwrap(),
            conductivity: &conductivity,
            contact_resistance_m2_k_w: CONTACT,
        },
    )]);
    let error = advance_growing_resistive_chamber(
        &[initial],
        &layers,
        None,
        None,
        &gas,
        &catalog,
        0.02,
        controls(),
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        GasResistiveAdvanceError::Thermal(DepositionError::MaterialDomain)
    ));
    assert!(conductivity.queries.get() >= 3);
    assert_eq!(initial.state, initial_state);
    assert_eq!(gas, before_gas);
}

#[test]
fn incompatible_density_volume_domain_and_cancellation_publish_no_state() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |h| fixture::temperature(h, chart);
    let initial = GasResistiveSurface {
        id: 4,
        law: &law,
        binding: &binding,
        state: state(chart, 2e-6),
        conductance_w_k: 0.0,
        substrate_temperature: &support,
        deposit_temperature: &deposit,
    };
    let model = |density| GasGrowingCoating {
        geometry: ConductiveCoatingGeometry::planar(density, AREA).unwrap(),
        conductivity: &Conductivity,
        contact_resistance_m2_k_w: CONTACT,
    };
    let mismatched = BTreeMap::from([(4, model(1000.0))]);
    assert!(matches!(
        advance_growing_resistive_chamber(
            &[initial],
            &mismatched,
            None,
            None,
            &gas,
            &catalog,
            0.0,
            controls(),
            || false
        ),
        Err(GasResistiveAdvanceError::Surface(GasSurfaceError::Surface(
            _
        )))
    ));
    let layers = BTreeMap::from([(4, model(DENSITY))]);
    let mut limited = controls();
    limited.maximum_deposit_volume_fraction = 1e-10;
    assert!(matches!(
        advance_growing_resistive_chamber(
            &[initial],
            &layers,
            None,
            None,
            &gas,
            &catalog,
            0.0,
            limited,
            || false
        ),
        Err(GasResistiveAdvanceError::Invalid(_))
    ));
    assert!(matches!(
        advance_growing_resistive_chamber(
            &[initial],
            &layers,
            None,
            None,
            &gas,
            &catalog,
            0.02,
            controls(),
            || true
        ),
        Err(GasResistiveAdvanceError::Surface(GasSurfaceError::Surface(
            SublimationError::Cancelled
        )))
    ));
}
