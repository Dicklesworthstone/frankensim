//! Original gas fixtures exercise source-integrated growing-layer transactions.
use super::*;
use crate::sublimation::gas::{
    chamber::{GasChamberSolve, coating::GasSupportedChamberSolve},
    tests as fixture,
};
use fs_thermal::sublimation::{SolidSurfaceState, SolverConfig, SublimationError};

const DENSITY: f64 = 917.0;
const AREA: f64 = 1e-4;
const CONTACT: f64 = 2e-4;
const CONSTANT_K: f64 = 1e-3;

struct ConstantConductivity;
impl ConductivityLaw for ConstantConductivity {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(CONSTANT_K)
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
            // The constant-k comparison admits the same pair of physical
            // half intervals in both independently selected interface paths.
            absolute_mass_tolerance_kg: 1.0,
            absolute_temperature_tolerance_k: 1000.0,
            absolute_specific_enthalpy_tolerance_j_kg: 1e7,
            absolute_heat_tolerance_j: 1e4,
            relative_tolerance: 1e-3,
            minimum_dt_s: 1e-10,
            maximum_dt_s: 0.02,
            maximum_trials: 20000,
        },
        maximum_deposit_volume_fraction: 0.01,
        absolute_extent_tolerance_m: 1.0,
        absolute_conductance_tolerance_w_k: 1e4,
    }
}

fn support(h: f64) -> Result<f64, SublimationError> {
    Ok(400.0 + h / 500.0)
}

fn state(ts: f64, td: f64, chart: (f64, f64), mass: f64) -> ResistiveSurfaceState {
    ResistiveSurfaceState {
        substrate_mass_kg: 0.02,
        substrate_specific_enthalpy_j_kg: 500.0 * (ts - 400.0),
        deposit: SolidSurfaceState {
            mass_kg: mass,
            specific_enthalpy_j_kg: chart.0 + chart.1 * (td - 400.0),
            // Actual mass-owned area must replace this stale input.
            area_m2: 7.0,
        },
    }
}

fn integrated_layer(curve: &ConductivityCurve, area: f64) -> GasIntegratedGrowingCoating<'_> {
    GasIntegratedGrowingCoating {
        geometry: ConductiveCoatingGeometry::planar(DENSITY, area).unwrap(),
        conductivity: curve,
        contact_resistance_m2_k_w: CONTACT,
        interface_config: LayerConductionConfig {
            max_iterations: 96,
            absolute_heat_tolerance_w: 1e-11,
            relative_heat_tolerance: 1e-12,
        },
    }
}

#[test]
fn constant_original_curve_matches_frozen_growth_as_deposition_changes_resistance() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let original_gas = gas.clone();
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |h| fixture::temperature(h, chart);
    let initial = GasResistiveSurface {
        id: 4,
        law: &law,
        binding: &binding,
        state: state(350.0, 350.0, chart, 2e-6),
        conductance_w_k: 1e9,
        substrate_temperature: &support,
        deposit_temperature: &deposit,
    };
    let curve = ConductivityCurve::new(vec![[200.0, CONSTANT_K], [900.0, CONSTANT_K]]).unwrap();
    let layer = integrated_layer(&curve, AREA);
    let integrated = BTreeMap::from([(4, layer)]);
    let frozen = BTreeMap::from([(
        4,
        GasGrowingCoating {
            geometry: layer.geometry,
            conductivity: &ConstantConductivity,
            contact_resistance_m2_k_w: CONTACT,
        },
    )]);
    let nonlinear = advance_growing_integrated_chamber(
        &[initial],
        &integrated,
        None,
        None,
        None,
        &gas,
        &catalog,
        0.02,
        controls(),
        || false,
    )
    .unwrap();
    let reference = advance_growing_resistive_chamber(
        &[initial],
        &frozen,
        None,
        None,
        &gas,
        &catalog,
        0.02,
        controls(),
        || false,
    )
    .unwrap();
    let end = nonlinear.history().endpoint().surfaces()[0].state;
    let expected = reference.history().endpoint().surfaces()[0].state;
    assert!(end.deposit.mass_kg > initial.state.deposit.mass_kg);
    fixture::near(end.deposit.mass_kg, expected.deposit.mass_kg, 1e-15);
    fixture::near(
        end.deposit.specific_enthalpy_j_kg,
        expected.deposit.specific_enthalpy_j_kg,
        1e-4,
    );
    fixture::near(
        end.substrate_specific_enthalpy_j_kg,
        expected.substrate_specific_enthalpy_j_kg,
        1e-7,
    );
    fixture::near(
        nonlinear.gas().state().state().temperature,
        reference.gas().state().state().temperature,
        1e-8,
    );
    fixture::near(
        nonlinear
            .gas()
            .species_mass_kg(&fixture::species())
            .unwrap(),
        reference
            .gas()
            .species_mass_kg(&fixture::species())
            .unwrap(),
        1e-14,
    );
    fixture::near(
        nonlinear.history().audit().thermal[0].interface_heat_j,
        reference.history().audit().thermal[0].interface_heat_j,
        1e-8,
    );
    assert_eq!(nonlinear.history().audit().chamber.accepted_intervals, 1);
    assert_eq!(reference.history().audit().chamber.accepted_intervals, 1);
    let geometry = nonlinear.geometry()[0];
    let thickness = end.deposit.mass_kg / (DENSITY * AREA);
    let expected_g = 1.0 / (thickness / (CONSTANT_K * AREA) + CONTACT / AREA);
    fixture::near(geometry.final_geometry.front.extent_m, thickness, 1e-18);
    fixture::near(geometry.final_geometry.conductance_w_k, expected_g, 1e-12);
    assert!(geometry.final_geometry.conductance_w_k < geometry.initial.conductance_w_k);
    assert_eq!(end.deposit.area_m2, AREA);
    assert_eq!(nonlinear.history().audit().chamber.supplied_heat_j, 0.0);
    fixture::near(
        nonlinear.history().audit().chamber.energy_residual_j,
        0.0,
        1e-7,
    );
    fixture::near(
        nonlinear.history().audit().chamber.mass_residual_kg,
        0.0,
        1e-14,
    );
    assert_eq!(gas, original_gas);
}

#[test]
fn nonlinear_source_and_signed_powers_preserve_separate_heat_and_energy_histories() {
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
        state: state(450.0, 350.0, chart, 0.002),
        conductance_w_k: 0.0,
        substrate_temperature: &support,
        deposit_temperature: &deposit,
    };
    // Original k(T)=1e-4*T. The interval crosses its 350 K source knot.
    let slope = 1e-4;
    let curve = ConductivityCurve::new(vec![[200.0, 0.02], [350.0, 0.035], [900.0, 0.09]]).unwrap();
    let area = 0.002;
    let layers = BTreeMap::from([(4, integrated_layer(&curve, area))]);
    let films = BTreeMap::from([(4, 0.0)]);
    let powers = BTreeMap::from([(
        4,
        GasResistivePower {
            substrate_power_w: -20.0,
            deposit_power_w: 3.0,
            interface: None,
        },
    )]);
    let duration = 2.0;
    let mut solve = controls();
    solve.advance.absolute_mass_tolerance_kg = 1e-9;
    solve.advance.absolute_temperature_tolerance_k = 0.01;
    solve.advance.absolute_specific_enthalpy_tolerance_j_kg = 2.0;
    solve.advance.absolute_heat_tolerance_j = 0.001;
    solve.advance.relative_tolerance = 1e-6;
    solve.advance.maximum_dt_s = duration;
    let result = advance_growing_integrated_chamber(
        &[initial],
        &layers,
        None,
        Some(&films),
        Some(&powers),
        &gas,
        &catalog,
        duration,
        solve,
        || false,
    )
    .unwrap();
    let audit = result.history().audit();
    let thermal = audit.thermal[0];
    let end = result.history().endpoint().surfaces()[0].state;
    let support_heat = -20.0 * duration;
    let deposit_heat = 3.0 * duration;
    assert!(audit.chamber.rejected_trials > 0);
    assert!(audit.chamber.accepted_intervals > 1);
    assert_eq!(audit.chamber.elapsed_s, duration);
    fixture::near(thermal.supplied_substrate_heat_j, support_heat, 1e-10);
    fixture::near(thermal.supplied_deposit_heat_j, deposit_heat, 1e-10);
    fixture::near(
        audit.chamber.supplied_heat_j,
        support_heat + deposit_heat,
        1e-10,
    );
    fixture::near(
        audit.chamber.surfaces[0].supplied_heat_j,
        support_heat + deposit_heat,
        1e-10,
    );

    // Independent endpoint energy differences, not the reported residuals.
    // The sealed mass film leaves both inventories and the actual gas unchanged.
    assert_eq!(end.substrate_mass_kg, initial.state.substrate_mass_kg);
    assert_eq!(end.deposit.mass_kg, initial.state.deposit.mass_kg);
    assert_eq!(result.gas(), &gas);
    assert_eq!(audit.chamber.surfaces[0].transferred_mass_kg, 0.0);
    assert_eq!(audit.chamber.surfaces[0].outgoing_mass_kg, 0.0);
    assert_eq!(audit.chamber.surfaces[0].incoming_mass_kg, 0.0);
    let support_storage = end.substrate_mass_kg
        * (end.substrate_specific_enthalpy_j_kg - initial.state.substrate_specific_enthalpy_j_kg);
    let deposit_storage = end.deposit.mass_kg
        * (end.deposit.specific_enthalpy_j_kg - initial.state.deposit.specific_enthalpy_j_kg);
    fixture::near(
        support_storage + deposit_storage,
        support_heat + deposit_heat,
        1e-7,
    );
    fixture::near(
        support_storage + thermal.interface_heat_j,
        support_heat,
        1e-7,
    );
    fixture::near(
        deposit_storage - thermal.interface_heat_j,
        deposit_heat,
        1e-7,
    );
    fixture::near(thermal.substrate_energy_change_j, support_storage, 1e-7);
    fixture::near(thermal.deposit_energy_change_j, deposit_storage, 1e-7);
    fixture::near(audit.chamber.energy_residual_j, 0.0, 1e-7);
    fixture::near(thermal.substrate_energy_residual_j, 0.0, 1e-7);
    fixture::near(thermal.deposit_energy_residual_j, 0.0, 1e-7);

    // An independent quadratic Fourier/contact solution checks final metadata.
    // b*T_face^2 + T_face = T_support + b*T_deposit^2.
    let ts = support(end.substrate_specific_enthalpy_j_kg).unwrap();
    let td = deposit(end.deposit.specific_enthalpy_j_kg).unwrap();
    let thickness = end.deposit.mass_kg / (DENSITY * area);
    let multiplier = area / thickness;
    let contact = CONTACT / area;
    let b = 0.5 * contact * multiplier * slope;
    let rhs = ts + b * td * td;
    let face = 2.0 * rhs / (1.0 + (1.0 + 4.0 * b * rhs).sqrt());
    let average_k = 0.5 * slope * (face + td);
    let expected_g = 1.0 / (1.0 / (multiplier * average_k) + contact);
    let frozen_g = 1.0 / (1.0 / (multiplier * slope * td) + contact);
    fixture::near(
        result.geometry()[0].final_geometry.conductance_w_k,
        expected_g,
        1e-10,
    );
    assert!(expected_g - frozen_g > 0.003);
    assert_eq!(result.geometry()[0].normal_recession_m, 0.0);
}

#[test]
fn unsupported_actual_coating_profile_refuses_without_source_extrapolation() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let original_gas = gas.clone();
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |h| fixture::temperature(h, chart);
    let original = state(450.0, 350.0, chart, 0.002);
    let rows = [GasResistiveSurface {
        id: 4,
        law: &law,
        binding: &binding,
        state: original,
        conductance_w_k: 0.0,
        substrate_temperature: &support,
        deposit_temperature: &deposit,
    }];
    let curve = ConductivityCurve::new(vec![[300.0, 0.03], [400.0, 0.04]]).unwrap();
    let layers = BTreeMap::from([(4, integrated_layer(&curve, 0.002))]);
    // The support's own caloric state is valid, but this small contact jump
    // cannot keep the ACTUAL coating face below its 400 K source ceiling.
    let result = advance_growing_integrated_chamber(
        &rows,
        &layers,
        None,
        Some(&BTreeMap::from([(4, 0.0)])),
        None,
        &gas,
        &catalog,
        0.02,
        controls(),
        || false,
    );
    assert!(matches!(
        result,
        Err(GasResistiveAdvanceError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(rows[0].state, original);
    assert_eq!(gas, original_gas);
}

#[test]
fn competing_interfaces_invalid_rosters_and_late_cancellation_publish_no_history() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let original_gas = gas.clone();
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |h| fixture::temperature(h, chart);
    let original = state(450.0, 350.0, chart, 0.002);
    let rows = [GasResistiveSurface {
        id: 4,
        law: &law,
        binding: &binding,
        state: original,
        conductance_w_k: 0.0,
        substrate_temperature: &support,
        deposit_temperature: &deposit,
    }];
    let curve = ConductivityCurve::new(vec![[200.0, 0.02], [900.0, 0.09]]).unwrap();
    let layer = integrated_layer(&curve, 0.002);
    let layers = BTreeMap::from([(4, layer)]);
    let films = BTreeMap::from([(4, 0.0)]);
    let valid_power = GasResistivePower {
        substrate_power_w: -20.0,
        deposit_power_w: 3.0,
        interface: None,
    };
    let competing = |_: f64, _: f64| Ok(0.0);
    for (id, power) in [
        (9, valid_power),
        (
            4,
            GasResistivePower {
                interface: Some(&competing),
                ..valid_power
            },
        ),
    ] {
        assert!(matches!(
            advance_growing_integrated_chamber(
                &rows,
                &layers,
                None,
                Some(&films),
                Some(&BTreeMap::from([(id, power)])),
                &gas,
                &catalog,
                0.02,
                controls(),
                || false,
            ),
            Err(GasResistiveAdvanceError::Invalid(_))
        ));
    }
    let wrong_layers = BTreeMap::from([(9, layer)]);
    assert!(matches!(
        advance_growing_integrated_chamber(
            &rows,
            &wrong_layers,
            None,
            Some(&films),
            None,
            &gas,
            &catalog,
            0.02,
            controls(),
            || false,
        ),
        Err(GasResistiveAdvanceError::Invalid(_))
    ));
    let powers = BTreeMap::from([(4, valid_power)]);
    let mut polls = 0;
    advance_growing_integrated_chamber(
        &rows,
        &layers,
        None,
        Some(&films),
        Some(&powers),
        &gas,
        &catalog,
        0.02,
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(polls > 16);
    let cancel_at = polls - 1;
    let mut replay_polls = 0;
    let result = advance_growing_integrated_chamber(
        &rows,
        &layers,
        None,
        Some(&films),
        Some(&powers),
        &gas,
        &catalog,
        0.02,
        controls(),
        || {
            replay_polls += 1;
            replay_polls == cancel_at
        },
    );
    // Source metadata and outer thermal polls retain their original typed
    // cancellation paths; neither can publish the completed candidate.
    assert!(matches!(
        result,
        Err(GasResistiveAdvanceError::Thermal(
            DepositionError::Cancelled
        )) | Err(GasResistiveAdvanceError::Surface(GasSurfaceError::Surface(
            SublimationError::Cancelled
        )))
    ));
    assert_eq!(replay_polls, cancel_at);
    assert_eq!(rows[0].state, original);
    assert_eq!(gas, original_gas);
}
