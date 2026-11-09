//! Native history adjoints against independent P1/contact/radiation equations.
use super::*;
use fs_conduction::bc::history::ThermalRegionHistory;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::ScanSensitivityError;
use fs_conduction::bc::history::march::enthalpy::sdirk::sensitivity::{
    Sdirk2HistoryScheduleGradient, Sdirk2SensitivityBudget,
};
use fs_conduction::transient::enthalpy::adjoint::EnthalpyAdjointError;
use fs_conduction::transient::enthalpy::adjoint::sdirk::{
    EnthalpySdirk2ParameterBudget, EnthalpySdirk2RadiationGradient,
};
use fs_scenario::bc::{BcKind, BcValue, BoundaryCondition, Physics};

const FIRST: [f64; 8] = [0.2, -0.1, 0.4, -0.3, 0.1, 0.3, -0.2, 0.5];
const LAST: [f64; 8] = [0.7, -0.4, 0.3, 0.2, -0.6, 0.9, -0.2, 0.5];
const HISTORY_DIRECTION: [f64; 8] = [0.8, -0.2, 0.4, -0.3, 0.1, 0.6, -0.5, 0.7];

fn budget() -> Sdirk2SensitivityBudget {
    Sdirk2SensitivityBudget { max_stage_source_vertices: 16, max_radiation_entries: 4096 }
}
fn linear() -> LinearConfig {
    LinearConfig { tolerance: 1e-12, max_iterations: 64, restart: 8 }
}
fn parameter_budget() -> EnthalpySdirk2ParameterBudget {
    EnthalpySdirk2ParameterBudget { max_contact_face_pairs_per_stage: 10_000,
        max_radiation_entries: 4096 }
}

// The no-contact paths need separated bodies, preserving both original volumes.
fn disconnected_tetrahedra() -> ConductionMesh {
    let mut positions = tetrahedra().positions().to_vec();
    for point in &mut positions[4..] { point[0] += 4.0; }
    ConductionMesh::new(TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
        positions).unwrap()
}

fn same_schedule_transport(a: &Sdirk2HistoryScheduleGradient, b: &Sdirk2HistoryScheduleGradient) {
    assert_eq!(a.sources, b.sources);
    assert_eq!(a.prescribed_temperatures, b.prescribed_temperatures);
    assert_eq!(a.physical.sources, b.physical.sources);
    let (a, b) = (&a.physical.transport, &b.physical.transport);
    assert_eq!(a.previous_specific_enthalpy, b.previous_specific_enthalpy);
    assert_eq!(a.source_density, b.source_density);
    assert_eq!(a.prescribed_temperature, b.prescribed_temperature);
    assert_eq!(a.stage_adjoint, b.stage_adjoint);
    assert_eq!(a.relative_residual, b.relative_residual);
    assert_eq!(a.iterations, b.iterations);
}

// Parameters are the original z-ramp endpoints, power-ramp endpoints and z width.
fn controlled(p: [f64; 5]) -> MovingHeatSource {
    let signal = |from, to, dims| TimeSignal::Ramp { t_start: 0.0, t_end: 0.1,
        from: QtyAny::new(from, dims), to: QtyAny::new(to, dims) };
    MovingHeatSource::new(
        [constant(0.2, LENGTH), constant(0.2, LENGTH), signal(p[0], p[1], LENGTH)],
        [2.0, 2.0, p[4]], signal(p[2], p[3], COMPONENT_POWER_DIMS), 0.0,
        MovingHeatBudget { max_vertices: 8, max_elements: 2,
            min_active_vertices: 8, power_relative_tolerance: 1e-12 },
    ).unwrap()
}

// Source normalization and P1 loads are independent of the production projection.
fn controlled_density(mesh: &ConductionMesh, p: [f64; 5], time: f64) -> [f64; 8] {
    let alpha = (time / 0.1).clamp(0.0, 1.0);
    let z = (1.0 - alpha) * p[0] + alpha * p[1];
    let watts = (1.0 - alpha) * p[2] + alpha * p[3];
    let weights: [f64; 8] = std::array::from_fn(|i| {
        let point = mesh.positions()[i];
        let r2 = ((point[0] - 0.2).powi(2) + (point[1] - 0.2).powi(2)) / 4.0
            + ((point[2] - z) / p[4]).powi(2);
        (1.0 - r2).max(0.0).powi(2)
    });
    let volume = weights.iter().sum::<f64>() / 24.0;
    weights.map(|weight| watts * weight / volume)
}

fn oracle(mesh: &ConductionMesh, p: [f64; 5], old: &[f64; 8], dt: f64,
    baths: Option<[f64; 2]>) -> [[f64; 8]; 2]
{
    let q = [G * dt, dt].map(|time| controlled_density(mesh, p, time));
    if let Some(baths) = baths {
        let first = radiative_stage(old, G * dt, &q[0], baths[0]);
        let carry = std::array::from_fn(|i| old[i] + (1.0 - G) / G * (first[i] - old[i]));
        [first, radiative_stage(&carry, G * dt, &q[1], baths[1])]
    } else { contact_stages(old, dt, &q) }
}
fn objective(stages: &[[f64; 8]; 2]) -> f64 {
    dot(&stages[0], &FIRST) + dot(&stages[1], &LAST)
}
fn derivative_close(actual: f64, expected: f64) {
    near(actual, expected, 3e-7 * (1.0 + actual.abs().max(expected.abs())));
}

#[test]
fn moving_sdirk_history_gradients_keep_both_stage_controls_and_original_history() {
    with_cx(|cx| {
        let mesh = tetrahedra(); let boundary = history(cx, &mesh); let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")]; let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0);
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let p = [-0.2, 0.3, 1.0, 2.0, 2.0]; let dt = 0.04;
        let create = || ContactEnthalpyHistoryRun::new_newton(cx,
            base(cx, &boundary, &conductivity, &curves, &OLD, zero(), 2),
            contacts(cx, &mesh, &initial, &cards, &curves, &source, &context,
                sides(&mesh, 4, 2, 0.0), 0, false), config()).unwrap().into_sdirk2().unwrap()
            .with_moving_sources(cx, vec![controlled(p)], 1).unwrap();
        let mut run = create();
        let binding = run.advance_linearized_to(cx, dt, config().step, budget()).unwrap();
        let primal = create().advance_to(cx, dt, config().step).unwrap();
        assert_eq!(binding.primal().first_stage.specific_enthalpy_j_kg, primal.first_stage.specific_enthalpy_j_kg);
        assert_eq!(binding.primal().second_stage.specific_enthalpy_j_kg, primal.second_stage.specific_enthalpy_j_kg);
        assert_eq!(binding.stage_times_s(), [G * dt, dt]);
        assert!(!binding.end_before());
        let expected = oracle(&mesh, p, &OLD, dt, None);
        admit(&binding.primal().first_stage, &expected[0], 1.0 + 10.0 * G * dt);
        admit(&binding.primal().second_stage, &expected[1], 1.0 + 10.0 * dt);
        // Nine authored coordinates are shared by the two stages, not charged twice.
        assert!(binding.pullback_schedule_stages(cx, [&FIRST, &LAST], linear(), 8).is_err());
        let gradient = binding.pullback_schedule_stages(cx, [&FIRST, &LAST], linear(), 9).unwrap();
        let full = binding.pullback_schedule_stages_with_parameters(cx, [&FIRST, &LAST],
            linear(), parameter_budget(), 9).unwrap();
        same_schedule_transport(&full.physical, &gradient);
        for stage in &full.stages {
            assert_eq!(stage.resistance_multipliers.len(), 1);
            let contact = &stage.resistance_multipliers[0];
            assert_eq!(contact.interface, "joint");
            assert_eq!(contact.card_identity, source.content_hash());
            assert_eq!(contact.face_pairs, 1);
            assert!(contact.derivative.is_finite());
            assert!(stage.radiation.is_none());
        }
        assert!(gradient.prescribed_temperatures.is_empty());
        let authored = &gradient.sources[0];
        let actual = [authored.center_values_m[2][0], authored.center_values_m[2][1],
            authored.absorbed_power_values_w[0], authored.absorbed_power_values_w[1], authored.half_width_m[2]];
        for coordinate in 0..5 {
            let (mut plus, mut minus) = (p, p); let epsilon = 2e-5;
            plus[coordinate] += epsilon; minus[coordinate] -= epsilon;
            let fd = (objective(&oracle(&mesh, plus, &OLD, dt, None))
                - objective(&oracle(&mesh, minus, &OLD, dt, None))) / (2.0 * epsilon);
            derivative_close(actual[coordinate], fd);
        }
        let epsilon = 1e-4;
        let plus = std::array::from_fn(|i| OLD[i] + epsilon * HISTORY_DIRECTION[i]);
        let minus = std::array::from_fn(|i| OLD[i] - epsilon * HISTORY_DIRECTION[i]);
        let fd = (objective(&oracle(&mesh, p, &plus, dt, None))
            - objective(&oracle(&mesh, p, &minus, dt, None))) / (2.0 * epsilon);
        derivative_close(dot(&gradient.physical.transport.previous_specific_enthalpy, &HISTORY_DIRECTION), fd);
        for stage in 0..2 {
            assert!(gradient.physical.sources[stage][0].absorbed_power_w.abs() > 1e-6);
            assert!(gradient.physical.transport.relative_residual[stage] <= linear().tolerance);
        }
        assert_eq!(binding.temperature_pullback(cx, [&FIRST, &LAST]).unwrap(), [vec![0.0; 8], vec![0.0; 8]]);
        let fraction = binding.liquid_mass_fraction_pullback(cx, [&FIRST, &LAST]).unwrap();
        for (stage, seed) in [&FIRST, &LAST].into_iter().enumerate() {
            for i in 0..8 { near(fraction[stage][i], seed[i] / 100.0, 1e-15); }
        }
        assert_eq!(run.conduction().accepted_steps(), 1);
    });
}

#[test]
fn moving_sdirk_held_schedule_gradients_use_original_stage_sides_without_extra_weights() {
    with_cx(|cx| {
        let mesh = tetrahedra(); let boundary = history(cx, &mesh); let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")]; let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0);
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let mass: [f64; 8] = std::array::from_fn(|i| if i < 4 { 1.0 / 24.0 } else { 1.0 / 12.0 });
        for left in [false, true] {
            let heater = moving(constant(-0.2, LENGTH), TimeSignal::Table {
                times: vec![0.0, 0.01, 0.02], values: vec![1.0, 0.0, 2.0],
                dims: COMPONENT_POWER_DIMS, interp: Interp::Hold });
            let mut run = ContactEnthalpyHistoryRun::new_newton(cx,
                base(cx, &boundary, &conductivity, &curves, &OLD, zero(), 2),
                contacts(cx, &mesh, &initial, &cards, &curves, &source, &context,
                    sides(&mesh, 4, 2, 0.0), 0, false), config()).unwrap().into_sdirk2().unwrap()
                .with_moving_sources(cx, vec![heater], 1).unwrap();
            let binding = if left { run.advance_linearized_event_step(cx, 0.03, 0.03, config().step, budget()).unwrap() }
                else { run.advance_linearized_to(cx, 0.01, config().step, budget()).unwrap() };
            assert_eq!(binding.end_time_s(), 0.01);
            assert_eq!(binding.end_before(), left);
            near(binding.primal().first_stage.source_w, 1.0, 1e-12);
            near(binding.primal().second_stage.source_w, if left { 1.0 } else { 0.0 }, 1e-12);
            let gradient = binding.pullback_schedule_stages(cx, [&[0.0; 8], &mass], linear(), 9).unwrap();
            let power = &gradient.sources[0].absorbed_power_values_w;
            near(power[0], if left { 0.01 } else { (1.0 - G) * 0.01 }, 1e-11);
            near(power[1], if left { 0.0 } else { G * 0.01 }, 1e-11);
            assert_eq!(power[2], 0.0);
            near(gradient.physical.sources[0][0].absorbed_power_w, (1.0 - G) * 0.01, 1e-11);
            near(gradient.physical.sources[1][0].absorbed_power_w, G * 0.01, 1e-11);
            for i in 0..8 { near(gradient.physical.transport.previous_specific_enthalpy[i], mass[i], 1e-11); }
        }
    });
}

#[test]
fn sdirk_prescribed_bath_controls_combine_original_stage_sides_without_reweighting() {
    with_cx(|cx| {
        let mesh = disconnected_tetrahedra();
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
            .remainder("bath", ThermalBc::dirichlet(350.0).unwrap()).unwrap(),
            &Environment::earth_lab(), vec![ThermalRegionHistory {
                row: BoundaryCondition { region: "bath".into(), physics: Physics::Thermal,
                    kind: BcKind::Dirichlet, value: Some(BcValue::Signal(TimeSignal::Table {
                        times: vec![0.0, 0.01, 0.02], values: vec![350.0, 360.0, 370.0],
                        dims: TEMPERATURE_DIMS, interp: Interp::Hold })),
                    compatibility: None, frame: 0 }, reference_temperature: None,
            }], 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")];
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 400.0)];
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        for left in [false, true] {
            let mut run = base(cx, &boundary, &conductivity, &curves, &[-50.0; 8], zero(), 1).into_sdirk2();
            let binding = if left { run.advance_linearized_event_step(cx, 0.03, 0.03, config().step, budget()).unwrap() }
                else { run.advance_linearized_to(cx, 0.01, config().step, budget()).unwrap() };
            assert_eq!(binding.end_time_s(), 0.01);
            assert_eq!(binding.end_before(), left);
            assert_eq!(binding.primal().first_stage.temperature, [350.0; 8]);
            assert_eq!(binding.primal().second_stage.temperature, [if left { 350.0 } else { 360.0 }; 8]);
            let seeds = binding.temperature_pullback(cx, [&FIRST, &LAST]).unwrap();
            assert!(binding.pullback_schedule_stages(cx, [&seeds[0], &seeds[1]], linear(), 2).is_err());
            let gradient = binding.pullback_schedule_stages(cx, [&seeds[0], &seeds[1]], linear(), 3).unwrap();
            assert!(gradient.sources.is_empty());
            assert_eq!(gradient.prescribed_temperatures.len(), 1);
            assert_eq!(gradient.prescribed_temperatures[0].region, "bath");
            let first = FIRST.iter().sum::<f64>(); let last = LAST.iter().sum::<f64>();
            let expected = if left { [first + last, 0.0, 0.0] } else { [first, last, 0.0] };
            for (actual, expected) in gradient.prescribed_temperatures[0].values_k.iter().zip(expected) {
                near(*actual, expected, 1e-12);
            }
            for stage in 0..2 {
                let expected = if stage == 0 { FIRST } else { LAST };
                for i in 0..8 { near(gradient.physical.transport.prescribed_temperature[stage][i], expected[i], 1e-12); }
            }
            for value in gradient.physical.transport.previous_specific_enthalpy { near(value, 0.0, 1e-12); }
        }
    });
}

#[test]
fn moving_sdirk_history_derivative_refusals_leave_all_accepted_state_unchanged() {
    with_cx(|cx| {
        let mesh = disconnected_tetrahedra(); let boundary = history(cx, &mesh);
        let cards = [material_card("cold"), material_card("hot")];
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let create = |old: &[f64; 8], power: f64| {
            let heater = MovingHeatSource::new(
                [constant(2.0, LENGTH), constant(0.2, LENGTH), constant(0.0, LENGTH)],
                [10.0; 3], constant(power, COMPONENT_POWER_DIMS), 0.0,
                MovingHeatBudget { max_vertices: 8, max_elements: 2,
                    min_active_vertices: 8, power_relative_tolerance: 1e-12 }).unwrap();
            base(cx, &boundary, &conductivity, &curves, old, zero(), 2)
                .into_sdirk2().with_moving_sources(cx, vec![heater], 1).unwrap()
        };
        let mut run = create(&OLD, 1.0);
        let initial_temperature = run.conduction().temperature().to_vec();
        let initial_fraction = run.conduction().liquid_mass_fraction().to_vec();
        support::with_cancelled_cx(|cancelled| assert!(run.advance_linearized_to(cancelled, 0.01,
            config().step, budget()).is_err()));
        assert!(run.advance_linearized_to(cx, 0.01, config().step,
            Sdirk2SensitivityBudget { max_stage_source_vertices: 15, ..budget() }).is_err());
        assert_eq!(run.conduction().specific_enthalpy_j_kg(), &OLD);
        assert_eq!(run.conduction().temperature(), initial_temperature);
        assert_eq!(run.conduction().liquid_mass_fraction(), initial_fraction);
        assert_eq!(run.conduction().time_s(), 0.0);
        assert_eq!(run.conduction().accepted_steps(), 0);
        let binding = run.advance_linearized_to(cx, 0.01, config().step, budget()).unwrap();
        let fresh = create(&OLD, 1.0).advance_to(cx, 0.01, config().step).unwrap();
        assert_eq!(binding.primal().second_stage.specific_enthalpy_j_kg, fresh.second_stage.specific_enthalpy_j_kg);
        let mut kink = create(&[0.0; 8], 0.0);
        let kink_temperature = kink.conduction().temperature().to_vec();
        let kink_fraction = kink.conduction().liquid_mass_fraction().to_vec();
        assert!(matches!(kink.advance_linearized_to(cx, 0.01, config().step, budget()),
            Err(ScanSensitivityError::Adjoint(EnthalpyAdjointError::ChartKink { .. }))));
        assert_eq!(kink.conduction().specific_enthalpy_j_kg(), &[0.0; 8]);
        assert_eq!(kink.conduction().temperature(), kink_temperature);
        assert_eq!(kink.conduction().liquid_mass_fraction(), kink_fraction);
        assert_eq!(kink.conduction().time_s(), 0.0);
        assert_eq!(kink.conduction().accepted_steps(), 0);
        kink.advance_to(cx, 0.01, config().step).unwrap();
    });
}

#[test]
fn moving_sdirk_history_adjoint_retains_material_reflection_and_phase_contact_feedback() {
    with_cx(|cx| {
        let mesh = tetrahedra(); let boundary = history(cx, &mesh); let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")]; let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0);
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let enclosure = open(cx, &mesh, &cards, &context);
        let radiation = PointwiseEnclosureRadiation::new_open(cx, &mesh, &initial, &enclosure,
            350.0, PointwiseEnclosureTolerance::default()).unwrap();
        let p = [-0.2, 0.3, 10.0, 510.0, 2.0]; let dt = 0.0004;
        let mut run = ContactEnthalpyHistoryRun::new_newton(cx,
            base(cx, &boundary, &conductivity, &curves, &OLD, zero(), 2),
            contacts(cx, &mesh, &initial, &cards, &curves, &source, &context,
                sides(&mesh, 4, 2, 0.0), 0, false), config()).unwrap()
            .with_material_enclosure(cx, &radiation, patches(&cards, &curves, &context), Some(bath(350.0))).unwrap()
            .into_sdirk2().unwrap().with_moving_sources(cx, vec![controlled(p)], 1).unwrap();
        assert!(run.advance_linearized_to(cx, dt, config().step,
            Sdirk2SensitivityBudget { max_radiation_entries: 0, ..budget() }).is_err());
        assert_eq!(run.conduction().specific_enthalpy_j_kg(), &OLD);
        assert_eq!(run.conduction().accepted_steps(), 0);
        let binding = run.advance_linearized_to(cx, dt, config().step, budget()).unwrap();
        let gradient = binding.pullback_schedule_stages(cx, [&FIRST, &LAST], linear(), 9).unwrap();
        let full = binding.pullback_schedule_stages_with_parameters(cx, [&FIRST, &LAST],
            linear(), parameter_budget(), 9).unwrap();
        same_schedule_transport(&full.physical, &gradient);
        for stage in &full.stages {
            assert_eq!(stage.resistance_multipliers.len(), 1);
            assert_eq!(stage.resistance_multipliers[0].interface, "joint");
            assert_eq!(stage.resistance_multipliers[0].card_identity, source.content_hash());
            let Some(EnthalpySdirk2RadiationGradient::MaterialEnclosure(radiation)) = &stage.radiation
                else { panic!("native material enclosure controls must retain their original owner"); };
            assert_eq!(radiation.patches.len(), 2);
            assert!(radiation.reservoir_temperature.is_some_and(f64::is_finite));
            for (body, patch) in radiation.patches.iter().enumerate() {
                assert_eq!(patch.region, REGIONS[body]);
                assert_eq!(patch.material_identity, cards[body].content_hash());
            }
        }
        let baths = Some([350.0; 2]);
        let expected = oracle(&mesh, p, &OLD, dt, baths);
        admit(&binding.primal().first_stage, &expected[0], 10.0 + 5000.0 * G * dt);
        admit(&binding.primal().second_stage, &expected[1], 10.0 + 5000.0 * dt);
        assert!(binding.primal().first_stage.pointwise_enclosure_radiation.is_some());
        assert!(binding.primal().second_stage.pointwise_enclosure_radiation.is_some());
        for (coordinate, actual) in [(2, gradient.sources[0].absorbed_power_values_w[0]),
            (4, gradient.sources[0].half_width_m[2])]
        {
            let (mut plus, mut minus) = (p, p); let epsilon = 1e-3;
            plus[coordinate] += epsilon; minus[coordinate] -= epsilon;
            let fd = (objective(&oracle(&mesh, plus, &OLD, dt, baths))
                - objective(&oracle(&mesh, minus, &OLD, dt, baths))) / (2.0 * epsilon);
            derivative_close(actual, fd);
        }
        let epsilon = 1e-4;
        let plus = std::array::from_fn(|i| OLD[i] + epsilon * HISTORY_DIRECTION[i]);
        let minus = std::array::from_fn(|i| OLD[i] - epsilon * HISTORY_DIRECTION[i]);
        let fd = (objective(&oracle(&mesh, p, &plus, dt, baths))
            - objective(&oracle(&mesh, p, &minus, dt, baths))) / (2.0 * epsilon);
        derivative_close(dot(&gradient.physical.transport.previous_specific_enthalpy, &HISTORY_DIRECTION), fd);
    });
}
