//! Declared analytic species/geometry, not measurements of a named instrument.
use super::*;
use crate::thermochemical_gas::tests::{model, point};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

const LENGTH: f64 = 0.57;
const RADIUS: f64 = 0.0073;
const FLOW: f64 = 0.00005;
const WALL_G: f64 = 0.03;

fn state(t: f64, fractions: Vec<f64>) -> PrimitiveGasState {
    PrimitiveGasState { temperature_k: t,
        composition: Composition::new(CompositionBasis::MoleFraction, fractions).unwrap() }
}
fn controls() -> SteadyGasSolve {
    SteadyGasSolve {
        caloric: GasEnergySolve { minimum_temperature_k: 275.0,
            maximum_temperature_k: 350.0, energy_tolerance_j_kg: 1e-8, maximum_iterations: 64 },
        maximum_dt_s: 0.05, maximum_relaxations: 2000, maximum_substeps: 100000,
        steady_power_tolerance_w: 1e-8, energy_tolerance_j: 1e-8, maximum_mean_mach: 0.05,
    }
}
fn fixture(count: usize) -> SteadyAxialGas {
    let dx = LENGTH / count as f64;
    SteadyAxialGas {
        source: model(vec![0.75, 0.25], 0.0), source_conditions: QueryPoint::new(),
        static_pressure_pa: 101325.0, phase: ThermochemicalPhase::FrozenGas,
        inlet: state(310.0, vec![0.75, 0.25]), inlet_mass_flow_kg_s: FLOW,
        cells: (0..count).map(|index| AxialGasCell {
            end_m: LENGTH * (index + 1) as f64 / count as f64,
            volume_m3: core::f64::consts::PI * RADIUS.powi(2) * dx,
            wall_area_m2: core::f64::consts::TAU * RADIUS * dx,
            wall_conductance_w_k: WALL_G / count as f64,
            wall_temperature_k: 290.0, initial_temperature_k: 290.0,
        }).collect(), injections: vec![],
    }
}

#[test]
fn wall_cooling_has_conservative_donor_cell_balance_and_observed_first_order() {
    let cp = R * 3.25 / (0.75 * 0.028 + 0.25 * 0.004);
    let ntu = WALL_G / (FLOW * cp);
    let continuum = 290.0 + 20.0 * (-ntu).exp();
    let mut errors = Vec::new();
    for count in [4, 8, 16, 32] {
        let problem = fixture(count);
        let solved = problem.solve(controls(), || false).unwrap();
        let outlet = solved.field().state_at(LENGTH).unwrap();
        let expected_discrete = 290.0 + 20.0 / (1.0 + ntu / count as f64).powi(count as i32);
        assert!((outlet.state().temperature - expected_discrete).abs() < 2e-4);
        let error = (outlet.state().temperature - continuum).abs();
        errors.push(error);
        assert!(solved.cell_power_residuals_w().iter()
            .all(|r| r.abs() <= controls().steady_power_tolerance_w));
        assert!(solved.maximum_energy_residual_j() <= controls().energy_tolerance_j);
        assert!(solved.substeps() > 0 && solved.substeps() <= controls().maximum_substeps);
        assert!(solved.maximum_mean_mach() <= controls().maximum_mean_mach);
        assert!((solved.field().state_at(0.0).unwrap().state().temperature - 310.0).abs() < 1e-9);
        outlet.state().validate_coherence().unwrap();
    }
    let n = errors.len();
    let order = (errors[n - 2] / errors[n - 1]).log2();
    assert!((0.90..1.05).contains(&order), "order={order}, errors={errors:?}");
}

#[test]
fn changing_only_species_data_changes_cooling_and_the_complete_acoustic_medium() {
    let mut problem = fixture(8);
    problem.inlet = state(310.0, vec![1.0, 0.0]);
    let a = problem.solve(controls(), || false).unwrap();
    problem.inlet = state(310.0, vec![0.0, 1.0]);
    let b = problem.solve(controls(), || false).unwrap();
    let a_out = a.field().state_at(LENGTH).unwrap();
    let b_out = b.field().state_at(LENGTH).unwrap();
    // The larger mass-specific cp of species B reduces cooling at fixed kg/s.
    assert!(b_out.state().temperature > a_out.state().temperature + 3.0);
    assert!(b_out.state().sound_speed > 2.0 * a_out.state().sound_speed);
    assert!(b_out.state().density < a_out.state().density / 4.0);
    assert_ne!(b_out.state().dynamic_viscosity, a_out.state().dynamic_viscosity);
    assert_ne!(b_out.state().thermal_conductivity, a_out.state().thermal_conductivity);
    assert_eq!(a.species().cells()[&CellId(7)], vec![1.0, 0.0]);
    assert_eq!(b.species().cells()[&CellId(7)], vec![0.0, 1.0]);
}

#[test]
fn downstream_species_injection_preserves_mass_and_species_specific_enthalpy() {
    let mut problem = fixture(4);
    problem.inlet = state(330.0, vec![1.0, 0.0]);
    for cell in &mut problem.cells {
        cell.wall_conductance_w_k = 0.0;
        cell.initial_temperature_k = 310.0;
    }
    problem.injections.push(GasInjection { cell: 2, mass_flow_kg_s: 3.0 * FLOW,
        state: state(290.0, vec![0.0, 1.0]) });
    let solution = problem.solve(controls(), || false).unwrap();
    assert_eq!(solution.species().cells()[&CellId(0)], vec![1.0, 0.0]);
    let fractions = &solution.species().cells()[&CellId(3)];
    assert!((fractions[0] - 0.25).abs() < 1e-15);
    assert!((fractions[1] - 0.75).abs() < 1e-15);
    assert_eq!(solution.mass_flow_kg_s(), &[FLOW, FLOW, 4.0 * FLOW, 4.0 * FLOW]);
    let expected = (0.25 * (3.5 / 0.028) * 330.0 + 0.75 * (2.5 / 0.004) * 290.0)
        / (0.25 * (3.5 / 0.028) + 0.75 * (2.5 / 0.004));
    let outlet = solution.field().state_at(LENGTH).unwrap();
    assert!((outlet.state().temperature - expected).abs() < 2e-4);
    assert!((expected - (0.25 * 330.0 + 0.75 * 290.0)).abs() > 5.0);
    // Every emitted local primitive evaluates through the same source/EOS owner.
    outlet.state().validate_coherence().unwrap();
}

#[test]
fn pressure_changes_density_and_flow_admission_without_reinventing_heat_capacity() {
    let low = fixture(4).solve(controls(), || false).unwrap();
    let mut high_problem = fixture(4);
    high_problem.static_pressure_pa *= 2.0;
    let high = high_problem.solve(controls(), || false).unwrap();
    let a = low.field().state_at(LENGTH).unwrap();
    let b = high.field().state_at(LENGTH).unwrap();
    assert!((b.state().temperature - a.state().temperature).abs() < 2e-4);
    assert!((b.state().density / a.state().density - 2.0).abs() < 2e-6);
    assert!(high.maximum_mean_mach() < 0.51 * low.maximum_mean_mach());
    assert_eq!(b.state().pressure, 202650.0);
}

#[test]
fn supersaturation_is_refused_instead_of_silently_condensing_or_clamping_species() {
    let mut problem = fixture(4);
    let original = model(vec![0.98, 0.02], 0.0);
    problem.source = ThermochemicalGasModel::new(original.thermodynamics().clone(),
        original.transport_materials().to_vec(), Some(&SpeciesId::new("B").unwrap())).unwrap();
    problem.phase = ThermochemicalPhase::LiquidWaterLimit;
    // 320 K: inside the Buck-1996 liquid-water window [253.15, 323.15] K.
    problem.inlet = state(320.0, vec![0.98, 0.02]);
    for cell in &mut problem.cells {
        cell.wall_temperature_k = 280.0;
        cell.initial_temperature_k = 280.0;
        cell.wall_conductance_w_k = 0.1;
    }
    assert!(problem.source.evaluate_with_liquid_water_limit(&point(320.0, 101325.0)).is_ok());
    let error = problem.solve(controls(), || false).unwrap_err();
    assert!(matches!(error, SteadyGasError::Material(_)));
    assert!(error.to_string().contains("supersaturated"));
}

#[test]
fn work_source_range_cancellation_and_unsupported_flow_are_explicit_refusals() {
    let problem = fixture(8);
    assert!(matches!(problem.solve(controls(), || true),
        Err(SteadyGasError::Thermal(DepositionError::Cancelled))));
    assert!(matches!(problem.solve(SteadyGasSolve {
        maximum_relaxations: 1, ..controls()
    }, || false), Err(SteadyGasError::NotConverged { .. })));
    let mut fast = fixture(1);
    fast.inlet_mass_flow_kg_s = 0.01;
    assert!(matches!(fast.solve(controls(), || false), Err(SteadyGasError::MeanFlow { .. })));
    let mut invalid = fixture(1);
    invalid.cells[0].initial_temperature_k = 360.0;
    assert!(invalid.solve(controls(), || false).is_err());
    invalid = fixture(1);
    invalid.cells[0].volume_m3 = f64::MIN_POSITIVE;
    assert!(invalid.solve(controls(), || false).is_err());
}

#[test]
fn loose_caloric_roots_cannot_cancel_a_declared_inlet_enthalpy_mismatch() {
    let mut problem = fixture(1);
    problem.inlet.temperature_k = 300.0;
    problem.cells[0].initial_temperature_k = 300.0;
    problem.cells[0].wall_conductance_w_k = 0.0;
    // Stored h is already stationary. Both coarse inverse roots would accept
    // 275 K and their identical errors would cancel in an emitted-grid balance.
    // The conservative closure-flux bound must still preserve the 300 K supply.
    let loose = SteadyGasSolve { caloric: GasEnergySolve {
        energy_tolerance_j_kg: 100000.0, ..controls().caloric
    }, ..controls() };
    assert!(matches!(problem.solve(loose, || false),
        Err(SteadyGasError::PublishedBalance { cell: 0, .. })));
    let tight = problem.solve(controls(), || false).unwrap();
    assert!((tight.field().state_at(LENGTH).unwrap().state().temperature - 300.0).abs() < 1e-9);
    assert!(tight.maximum_power_residual_bound_w() <= controls().steady_power_tolerance_w);
}
