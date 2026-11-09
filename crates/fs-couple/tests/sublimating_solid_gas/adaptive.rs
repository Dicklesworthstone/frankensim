use super::fixture::*;
use fs_couple::liquid_gas::evaporation::momentum::MovingGasParcel;
use fs_couple::solid_gas::{
    kinetic::{KineticSublimationControls, MovingSolidSurface},
    trajectory::{
        SublimationAdaptiveControls, SublimationTrajectoryCandidate, SublimationTrajectoryError,
        advance_sublimating_trajectory,
    },
};
use fs_material::sublimation::gas::{GasSurfaceBinding, chamber::GasChamberSurface};
use fs_material::thermochemical_gas::ThermochemicalGasModel;

fn accuracy() -> SublimationAdaptiveControls {
    SublimationAdaptiveControls {
        initial_step_s: 0.05,
        minimum_step_s: 1e-8,
        maximum_attempts: 4096,
        maximum_accepted_intervals: 2048,
        relative_tolerance: 1e-4,
        mass_tolerance_kg: 1e-11,
        temperature_tolerance_k: 1e-5,
        pressure_tolerance_pa: 0.01,
        specific_enthalpy_tolerance_j_kg: 0.02,
        velocity_tolerance_m_s: 1e-4,
    }
}

fn with_area(body: MovingSolidSurface, area_m2: f64) -> MovingSolidSurface {
    let mut state = body.solid();
    state.area_m2 = area_m2;
    MovingSolidSurface::from_momentum(state, body.momentum_kg_m_s()).unwrap()
}

#[allow(clippy::too_many_arguments)]
fn run(
    body: &MovingSolidSurface,
    reservoir: &MovingGasParcel,
    model: &ThermochemicalGasModel,
    binding: &GasSurfaceBinding,
    caloric: (f64, f64),
    duration: f64,
    watts: f64,
    physical: KineticSublimationControls,
    adaptive: SublimationAdaptiveControls,
    cancelled: impl FnMut() -> bool,
) -> Result<SublimationTrajectoryCandidate, SublimationTrajectoryError> {
    let law = phase();
    let temperature = |h| temperature(h, caloric);
    let source = GasChamberSurface {
        id: 7,
        law: &law,
        binding,
        solid: body.solid(),
        temperature_from_enthalpy: &temperature,
    };
    advance_sublimating_trajectory(&source, body, reservoir, model, frame(), enthalpy_range(caloric),
        None, duration, watts, physical, adaptive, cancelled)
}

#[test]
fn constant_heat_and_short_final_remainder_match_the_original_caloric_solution() {
    let model = catalog(0.0);
    let binding = binding(0.0);
    let caloric = chart(&model, &binding);
    let body = with_area(solid(400.0, caloric, [2.0, -3.0, 4.0]), 0.0);
    let reservoir = gas(&model, 400.0, 0.0, [0.0; 3]);
    let mut adaptive = accuracy();
    adaptive.minimum_step_s = 0.02;
    let result = run(&body, &reservoir, &model, &binding, caloric, 0.11, 0.2,
        controls(), adaptive, || false).unwrap();
    near(temperature(result.solid().solid().specific_enthalpy_j_kg, caloric).unwrap(),
        400.0 + 0.2 * 0.11 / MASS / caloric.1, 1e-11);
    assert_eq!(result.gas(), &reservoir);
    assert_eq!(result.solid().momentum_kg_m_s(), body.momentum_kg_m_s());
    near(result.step().supplied_heat_j, 0.022, 1e-15);
    near(result.step().accepted_steps.iter().map(|step| step.duration_s).sum(), 0.11, 1e-15);
    assert_eq!(result.step().accepted_intervals, 3);
    assert_eq!(result.step().accepted_steps.len(), 6);
    assert_eq!(result.step().source_evaluations, 1 + 3 * result.step().attempted_intervals);
    assert_eq!(result.step().outgoing_mass_kg, 0.0);
    assert_eq!(result.step().mass_residual_kg, 0.0);
    assert_eq!(result.step().momentum_residual_n_s, [0.0; 3]);
    near(total_energy(result.solid(), result.gas(), 0.0) - total_energy(&body, &reservoir, 0.0),
        0.022, 1e-8);
}

#[test]
fn admissible_minimum_is_tried_before_a_conservative_retry_refuses() {
    let model = catalog(0.0);
    let binding = binding(0.0);
    let caloric = chart(&model, &binding);
    let body = with_area(solid(400.0, caloric, [0.0; 3]), 0.0);
    let reservoir = gas(&model, 400.0, 0.0, [0.0; 3]);
    let mut adaptive = accuracy();
    adaptive.minimum_step_s = 0.04;
    let mut physical = controls();
    physical.maximum_temperature_change_k = 0.00045;
    let watts = MASS * caloric.1 * 0.01;
    let result = run(&body, &reservoir, &model, &binding, caloric, 0.08, watts,
        physical, adaptive, || false).unwrap();
    assert_eq!(result.step().rejected_intervals, 1);
    assert_eq!(result.step().accepted_intervals, 2);
    near(temperature(result.solid().solid().specific_enthalpy_j_kg, caloric).unwrap(),
        400.0008, 1e-10);
    near(result.step().supplied_heat_j, watts * 0.08, 1e-14);
}

#[test]
fn tighter_error_control_resolves_evolving_backpressure_and_recoil() {
    let model = catalog(0.0);
    let binding = binding(0.0);
    let caloric = chart(&model, &binding);
    let body = with_area(solid(430.0, caloric, [0.0; 3]), 100.0 * AREA);
    let reservoir = gas(&model, 350.0, 400.0, [5.0, -20.0, 0.0]);
    let mut loose_controls = accuracy();
    loose_controls.relative_tolerance = 1e-2;
    loose_controls.velocity_tolerance_m_s = 0.05;
    loose_controls.temperature_tolerance_k = 0.05;
    loose_controls.pressure_tolerance_pa = 10.0;
    loose_controls.specific_enthalpy_tolerance_j_kg = 100.0;
    loose_controls.mass_tolerance_kg = 1e-7;
    let mut tight_controls = accuracy();
    tight_controls.relative_tolerance = 1e-5;
    tight_controls.velocity_tolerance_m_s = 1e-5;
    let execute = |control| run(&body, &reservoir, &model, &binding, caloric, 0.1, 0.0,
        controls(), control, || false).unwrap();
    let loose = execute(loose_controls);
    let tight = execute(tight_controls);
    assert!(tight.step().accepted_intervals > loose.step().accepted_intervals);
    assert!(tight.step().rejected_intervals > 0);
    assert!(tight.step().maximum_accepted_error_ratio <= 1.0);
    let (mut reference_body, mut reference_gas) = (body, reservoir.clone());
    // Independently prescribed uniform refinement, without the adaptive controller.
    for _ in 0..1024 {
        let next = advance(&reference_body, &reference_gas, &model, &binding, caloric, None,
            0.1 / 1024.0, 0.0, controls()).unwrap();
        (reference_body, reference_gas, _) = next.into_parts();
    }
    let distance = |candidate: &SublimationTrajectoryCandidate| {
        let relative = candidate.gas().velocity_m_s()[0] - candidate.solid().velocity_m_s()[0];
        let reference = reference_gas.velocity_m_s()[0] - reference_body.velocity_m_s()[0];
        (relative - reference).abs()
    };
    assert!(distance(&tight) < 0.35 * distance(&loose));
    near(total_energy(tight.solid(), tight.gas(), 0.0) - total_energy(&body, &reservoir, 0.0),
        0.0, 1e-5);
    near(tight.step().mass_residual_kg, 0.0, 1e-12);
    for value in tight.step().momentum_residual_n_s { near(value, 0.0, 1e-10); }
    assert_eq!(tight, execute(tight_controls));
}

#[test]
fn donor_limits_retry_without_leaking_rejected_heat_mass_or_impulse() {
    let model = catalog(0.0);
    let binding = binding(0.0);
    let caloric = chart(&model, &binding);
    let body = with_area(solid(400.0, caloric, [0.0; 3]), 100.0 * AREA);
    let reservoir = gas(&model, 400.0, 1000.0, [0.0; 3]);
    let mut physical = controls();
    physical.maximum_gross_inventory_fraction = 1e-4;
    let result = run(&body, &reservoir, &model, &binding, caloric, 0.05, 0.3,
        physical, accuracy(), || false).unwrap();
    assert!(result.step().rejected_intervals > 0);
    assert_eq!(result.step().accepted_steps.len(), 2 * result.step().accepted_intervals);
    near(result.step().supplied_heat_j, 0.015, 1e-14);
    near(result.step().accepted_steps.iter().map(|s| s.duration_s).sum(), 0.05, 1e-14);
    near(body.solid().mass_kg - result.solid().solid().mass_kg,
        result.step().transferred_mass_kg, 1e-13);
    near(total_energy(result.solid(), result.gas(), 0.0) - total_energy(&body, &reservoir, 0.0),
        0.015, 1e-5);
    assert!(result.step().source_evaluations <= 1 + 3 * result.step().attempted_intervals);
    let carrier = fs_thermochem::SpeciesId::new("A").unwrap();
    near(result.gas().parcel().species_mass_kg(&carrier).unwrap(),
        reservoir.parcel().species_mass_kg(&carrier).unwrap(), 1e-12);
}

#[test]
fn budgets_minimum_steps_and_late_cancellation_preserve_the_original_pair() {
    let model = catalog(0.0);
    let binding = binding(0.0);
    let caloric = chart(&model, &binding);
    let body = with_area(solid(400.0, caloric, [0.0; 3]), 100.0 * AREA);
    let reservoir = gas(&model, 400.0, 1000.0, [0.0; 3]);
    let before = (body, reservoir.clone());
    let mut tight = accuracy();
    tight.relative_tolerance = 0.0;
    tight.velocity_tolerance_m_s = 1e-12;
    tight.maximum_attempts = 1;
    assert!(matches!(run(&body, &reservoir, &model, &binding, caloric, 0.1, 0.0,
        controls(), tight, || false), Err(SublimationTrajectoryError::WorkBudget { .. })));
    tight.maximum_attempts = 4096;
    tight.minimum_step_s = tight.initial_step_s;
    assert!(matches!(run(&body, &reservoir, &model, &binding, caloric, 0.1, 0.0,
        controls(), tight, || false), Err(SublimationTrajectoryError::MinimumStep { .. })));
    let mut calls = 0;
    run(&body, &reservoir, &model, &binding, caloric, 0.1, 0.0, controls(), accuracy(),
        || { calls += 1; false }).unwrap();
    let mut repeated = 0;
    assert!(matches!(run(&body, &reservoir, &model, &binding, caloric, 0.1, 0.0,
        controls(), accuracy(), || { repeated += 1; repeated == calls - 1 }),
        Err(SublimationTrajectoryError::Cancelled { .. })));
    assert_eq!((body, reservoir), before);
}

#[test]
fn zero_time_admits_sources_and_has_no_fictitious_heat_or_history() {
    let model = catalog(0.0);
    let binding = binding(0.0);
    let caloric = chart(&model, &binding);
    let body = solid(400.0, caloric, [0.0; 3]);
    let reservoir = gas(&model, 400.0, 0.0, [0.0; 3]);
    let result = run(&body, &reservoir, &model, &binding, caloric, 0.0, 1e9,
        controls(), accuracy(), || false).unwrap();
    assert_eq!(result.solid(), &body);
    assert_eq!(result.gas(), &reservoir);
    assert_eq!(result.step().supplied_heat_j, 0.0);
    assert!(result.step().accepted_steps.is_empty());
    assert_eq!(result.step().source_evaluations, 1);
    let mut bad = binding.clone();
    bad.solid_enthalpy_offset_j_kg += LATENT;
    assert!(matches!(run(&body, &reservoir, &model, &bad, caloric, 0.0, 0.0,
        controls(), accuracy(), || false), Err(SublimationTrajectoryError::Step { .. })));
}
