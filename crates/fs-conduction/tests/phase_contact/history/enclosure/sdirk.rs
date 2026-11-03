//! Stage-aware contact/reflection histories compared with independent heat equations.
use super::*;
use fs_conduction::bc::history::march::enthalpy::radiation::enclosure::MaterialEnclosureEnthalpyHistoryRun;
use fs_conduction::transient::enthalpy::sdirk::{EnthalpySdirk2Stage, EnthalpySdirk2Step, SDIRK2_GAMMA as G};

fn loads(h: &[f64; 8], mutual: f64, bath: f64, contact: bool) -> ([f64; 8], [f64; 2]) {
    let (mut outward, radiation) = flux(h, mutual, bath);
    if !contact {
        let fa = h[..3].iter().sum::<f64>() / 300.0;
        let fb = h[4..7].iter().sum::<f64>() / 300.0;
        let q = AREA * (400.0 - 500.0) / (1.0 + 2.0 * fa + fb);
        for i in 0..3 { outward[i] -= q / 3.0; outward[4 + i] += q / 3.0; }
    }
    (outward, radiation)
}

fn stage_rhs(old: &[f64; 8], h: &[f64; 8], dt: f64, mutual: f64, bath: f64,
    source: f64, contact: bool) -> [f64; 8]
{
    let (outward, _) = loads(h, mutual, bath, contact);
    std::array::from_fn(|i| old[i] + dt * (source / 24.0 - outward[i])
        * if i < 4 { 24.0 } else { 12.0 })
}

fn stages(old: &[f64; 8], dt: f64, mutual: f64, bath: f64, source: f64, contact: bool)
    -> [[f64; 8]; 2]
{
    let mut first = *old;
    for _ in 0..128 { first = stage_rhs(old, &first, G * dt, mutual, bath, source, contact); }
    // The RHS is algebra only: no source/chart is queried at this extrapolation.
    let rhs = stage_rhs(old, &first, (1.0 - G) * dt, mutual, bath, source, contact);
    let mut second = first;
    for _ in 0..128 { second = stage_rhs(&rhs, &second, G * dt, mutual, bath, source, contact); }
    [first, second]
}

fn admit_stage(stage: &EnthalpySdirk2Stage, expected: &[f64; 8], mutual: f64, bath: f64, contact: bool) {
    for i in 0..8 {
        near(stage.specific_enthalpy_j_kg[i], expected[i], 4e-8);
        near(stage.liquid_mass_fraction[i], expected[i] / 100.0, 4e-10);
    }
    assert_eq!(stage.temperature, [400.0, 400.0, 400.0, 400.0, 500.0, 500.0, 500.0, 500.0]);
    let radiation = stage.pointwise_enclosure_radiation.as_ref().unwrap();
    let (_, heat) = loads(expected, mutual, bath, contact);
    near(stage.radiation_out_w, heat.iter().sum(), 3e-7);
    near(radiation.system_energy_closure_w, 0.0, 1e-8);
    for (body, patch) in radiation.patches.iter().enumerate() {
        near(patch.outward_heat_w, heat[body], 3e-7);
        for state in &patch.material_states {
            near(state.emissivity.value(), 0.2 + 0.006 * expected[state.vertex], 4e-10);
        }
    }
    if let Some(bath_report) = &radiation.ambient {
        assert_eq!(bath_report.temperature_k, bath);
        near(radiation.outward_heat_w + bath_report.outward_heat_w, 0.0, 1e-8);
    } else { near(radiation.outward_heat_w, 0.0, 1e-8); }
    assert_eq!(stage.physical_contacts.is_some(), contact);
    if let Some(contacts) = &stage.physical_contacts {
        for trace in contacts.material_trace_states() {
            for body in 0..2 { near(trace.liquid_mass_fractions[body],
                expected[4 * body..4 * body + 3].iter().sum::<f64>() / 300.0, 4e-10); }
        }
    }
    assert!(stage.stage_energy_residual_j.abs() <= config().step.energy_tolerance_j);
    assert!(stage.stage_residual_norm_j <= config().residual_tolerance_j);
}

fn admit(step: &EnthalpySdirk2Step, old: &[f64; 8], expected: &[[f64; 8]; 2],
    dt: f64, mutual: f64, bath: f64, source: f64, contact: bool)
{
    admit_stage(&step.first_stage, &expected[0], mutual, bath, contact);
    admit_stage(&step.second_stage, &expected[1], mutual, bath, contact);
    let (_, q1) = loads(&expected[0], mutual, bath, contact);
    let (_, q2) = loads(&expected[1], mutual, bath, contact);
    let radiated = dt * ((1.0 - G) * q1.iter().sum::<f64>() + G * q2.iter().sum::<f64>());
    let stored: f64 = expected[1].iter().zip(old).enumerate()
        .map(|(i, (new, old))| (new - old) * if i < 4 { 1.0 / 24.0 } else { 1.0 / 12.0 }).sum();
    near(step.source_j, dt * source / 3.0, 2e-9);
    near(step.radiation_out_j, radiated, 2e-9);
    near(step.stored_energy_change_j, stored, 2e-9);
    near(stored + radiated - dt * source / 3.0, 0.0, 2e-9);
    assert!(step.energy_residual_j.abs() <= config().step.energy_tolerance_j);
}

#[test]
fn sdirk_contact_enclosure_events_evaluate_both_phase_laws_and_integrate_heat_once() {
    with_cx(|cx| {
        let mesh = tetrahedra(); let boundary = history(cx, &mesh); let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")]; let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0);
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let closed = closed(&mesh, &cards, &context); let open = open(cx, &mesh, &cards, &context);
        for mutual in [1.0, 0.4] {
            let radiation = if mutual == 1.0 { PointwiseEnclosureRadiation::new(cx, &mesh, &initial, &closed,
                PointwiseEnclosureTolerance::default()).unwrap() }
            else { PointwiseEnclosureRadiation::new_open(cx, &mesh, &initial, &open,
                350.0, PointwiseEnclosureTolerance::default()).unwrap() };
            let forcing = HeatSourceHistory::Uniform(TimeSignal::Table { times: vec![0.0, 0.0004],
                values: vec![0.0, 60.0], dims: VOLUMETRIC_SOURCE_DIMS, interp: Interp::Hold });
            let ambient = (mutual != 1.0).then_some(TimeSignal::Table { times: vec![0.0, 0.0002, 0.0004],
                values: vec![350.0, 600.0, 450.0], dims: TEMPERATURE_DIMS, interp: Interp::Hold });
            let mut run = ContactEnthalpyHistoryRun::new_newton(cx,
                base(cx, &boundary, &conductivity, &curves, &OLD, forcing, 3),
                contacts(cx, &mesh, &initial, &cards, &curves, &source, &context,
                    sides(&mesh, 4, 2, 0.0), 0, false), config()).unwrap()
                .with_material_enclosure(cx, &radiation, patches(&cards, &curves, &context), ambient).unwrap()
                .into_sdirk2().unwrap();
            assert!(run.contacts().is_some());
            assert_eq!(run.material_enclosure_patches().unwrap().len(), 2);
            let intervals = if mutual == 1.0 { vec![(0.0004, 350.0, 0.0), (0.0007, 350.0, 60.0)] }
                else { vec![(0.0002, 350.0, 0.0), (0.0004, 600.0, 0.0), (0.0007, 450.0, 60.0)] };
            let mut old = OLD; let mut from = 0.0;
            for (to, bath, q) in intervals {
                let dt = to - from; let expected = stages(&old, dt, mutual, bath, q, true);
                let actual = run.advance_event_step(cx, 0.0007, 1.0, config().step).unwrap();
                assert_eq!(actual.start_time_s, from); assert_eq!(actual.end_time_s, to);
                admit(&actual.report, &old, &expected, dt, mutual, bath, q, true);
                old = expected[1]; from = to;
            }
            assert_eq!(run.conduction().time_s(), 0.0007);
            assert_eq!(run.conduction().accepted_steps(), if mutual == 1.0 { 2 } else { 3 });
        }
    });
}

#[test]
fn sdirk_material_enclosure_adapts_on_latent_error_and_retains_bath_after_failed_proposals() {
    with_cx(|cx| {
        // Radiation-only bodies have disjoint physical boundaries. Unlike the
        // contact fixture, no coincident face may silently lose its interface.
        let mesh = ConductionMesh::new(TetComplex::from_tets(8,
            vec![[0, 1, 2, 3], [4, 5, 6, 7]]), vec![
            [0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., -1.],
            [2., 0., 0.], [3., 0., 0.], [2., 1., 0.], [2., 0., 1.],
        ]).unwrap();
        let boundary = history(cx, &mesh); let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")]; let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let enclosure = open(cx, &mesh, &cards, &context);
        let radiation = PointwiseEnclosureRadiation::new_open(cx, &mesh, &initial, &enclosure,
            350.0, PointwiseEnclosureTolerance::default()).unwrap();
        let create = |bath| MaterialEnclosureEnthalpyHistoryRun::new(cx,
            base(cx, &boundary, &conductivity, &curves, &OLD, zero(), 2), &radiation,
            patches(&cards, &curves, &context), Some(bath)).unwrap().into_sdirk2();
        let mut run = create(bath(350.0));
        assert!(run.contacts().is_none());
        let accepted = run.advance_adaptive_step(cx, 0.001, 0.001, config().step,
            AdaptiveHistoryConfig { absolute_temperature_k: 1e-10, relative_temperature: 0.0,
                min_step_s: 1e-8, max_step_s: 0.001, max_attempts: 16 },
            PhaseErrorTolerance { absolute_enthalpy_j_kg: 100.0, relative_enthalpy: 0.0,
                absolute_liquid_fraction: 1e-9 }).unwrap();
        assert!(accepted.rejected_attempts > 0 && accepted.end_time_s < 0.001);
        let first_dt = accepted.midpoint_time_s;
        let second_dt = accepted.end_time_s - accepted.midpoint_time_s;
        let first = stages(&OLD, first_dt, 0.4, 350.0, 0.0, false);
        let second = stages(&first[1], second_dt, 0.4, 350.0, 0.0, false);
        admit(&accepted.first_half, &OLD, &first, first_dt, 0.4, 350.0, 0.0, false);
        admit(&accepted.second_half, &first[1], &second, second_dt, 0.4, 350.0, 0.0, false);
        assert_eq!(run.conduction().accepted_steps(), 2);
        assert_eq!(run.conduction().specific_enthalpy_j_kg(), accepted.second_half.second_stage.specific_enthalpy_j_kg);

        let mut invalid = create(TimeSignal::Table { times: vec![0.0, 0.0004], values: vec![350.0, -10.0],
            dims: TEMPERATURE_DIMS, interp: Interp::Hold });
        invalid.advance_to(cx, 0.0002, config().step).unwrap();
        let saved = invalid.conduction().specific_enthalpy_j_kg().to_vec();
        assert!(invalid.advance_to(cx, 0.0005, config().step).is_err());
        support::with_cancelled_cx(|cancelled|
            assert!(invalid.advance_to(cancelled, 0.0003, config().step).is_err()));
        assert_eq!(invalid.conduction().specific_enthalpy_j_kg(), saved);
        assert_eq!(invalid.conduction().time_s(), 0.0002);
        assert_eq!(invalid.conduction().accepted_steps(), 1);
        let resumed = invalid.advance_to(cx, 0.0003, config().step).unwrap();
        let mut fresh = create(bath(350.0)); fresh.advance_to(cx, 0.0002, config().step).unwrap();
        let replay = fresh.advance_to(cx, 0.0003, config().step).unwrap();
        assert_eq!(resumed.second_stage.specific_enthalpy_j_kg, replay.second_stage.specific_enthalpy_j_kg);
    });
}

#[test]
fn sdirk_adjoint_retains_heterogeneous_latent_contact_and_material_reflection_feedback() {
    use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;
    use fs_conduction::transient::enthalpy::sdirk::EnthalpySdirk2Input;
    with_cx(|cx| {
        let mesh = tetrahedra(); let history = history(cx, &mesh); let boundary = history.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")]; let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0);
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let contact = contacts(cx, &mesh, &boundary, &cards, &curves, &source, &context,
            sides(&mesh, 4, 2, 0.0), 0, false);
        let enclosure = open(cx, &mesh, &cards, &context);
        let radiation = PointwiseEnclosureRadiation::new_open(cx, &mesh, &boundary, &enclosure,
            350.0, PointwiseEnclosureTolerance::default()).unwrap();
        let declarations = patches(&cards, &curves, &context);
        let material = MaterialPointwiseEnclosureRadiation::new(cx, &radiation, &declarations).unwrap();
        let engine = engine(cx, &mesh, &curves, 4);
        let zero = ScalarField::Uniform(0.0);
        let input = |q| EnthalpySdirk2Input {
            problem: ConductionProblem { mesh: &mesh, boundary: &boundary, material: &conductivity,
                element_materials: None, source: q }, interfaces: None,
            radiation: EnthalpyRadiation::MaterialEnclosure(&material), contacts: Some((&contact, config())),
        };
        let dt = 0.0005; let allowance = 2 * material.adjoint_workspace_entries().unwrap();
        assert!(engine.linearize_sdirk2(cx, [input(&zero); 2], &OLD, dt, config().step, allowance - 1).is_err());
        let response = engine.linearize_sdirk2(cx, [input(&zero); 2], &OLD, dt, config().step, allowance).unwrap();
        let expected = stages(&OLD, dt, 0.4, 350.0, 0.0, true);
        admit(response.primal(), &OLD, &expected, dt, 0.4, 350.0, 0.0, true);
        let first = [0.03, -0.01, 0.04, -0.02, 0.02, -0.03, 0.01, 0.05];
        let last = [0.7, -0.4, 0.3, 0.2, -0.6, 0.9, -0.2, 0.5];
        let direction = [0.8, -0.2, 0.4, -0.3, 0.1, 0.6, -0.5, 0.7];
        let gradient = response.pullback_stages(cx, [&first, &last], LinearConfig {
            tolerance: 1e-12, max_iterations: 64, restart: 8 }).unwrap();
        assert_eq!(response.temperature_pullback(cx, [&first, &last]).unwrap(), [vec![0.0; 8], vec![0.0; 8]]);
        // The oracle solves the original analytic two-body heat equations,
        // including R(fa,fb), correlated P1 epsilon and reflected irradiation.
        let objective = |h: [f64; 8]| {
            let stages = stages(&h, dt, 0.4, 350.0, 0.0, true);
            dot(&stages[0], &first) + dot(&stages[1], &last)
        };
        let epsilon = 0.001;
        let plus = objective(std::array::from_fn(|i| OLD[i] + epsilon * direction[i]));
        let minus = objective(std::array::from_fn(|i| OLD[i] - epsilon * direction[i]));
        near(dot(&gradient.previous_specific_enthalpy, &direction), (plus - minus) / (2.0 * epsilon), 2e-8);
        assert!(gradient.previous_specific_enthalpy.iter().enumerate()
            .any(|(i, g)| (g - first[i] - last[i]).abs() > 1e-3));
        let v = [direction, last]; let w = [last, first];
        let action = response.apply_jacobian(cx, [&v[0], &v[1]]).unwrap();
        let transpose = response.apply_jacobian_transpose(cx, [&w[0], &w[1]]).unwrap();
        near(dot(&action[0], &w[0]) + dot(&action[1], &w[1]),
            dot(&v[0], &transpose[0]) + dot(&v[1], &transpose[1]), 1e-13);
        for stage in 0..2 {
            let evaluate = |delta| {
                let q = ScalarField::nodal("perturbed stage deposition", 8,
                    direction.map(|d| delta * d).to_vec()).unwrap();
                let inputs = [0, 1].map(|sample| EnthalpySdirk2Input {
                    problem: ConductionProblem { mesh: &mesh, boundary: &boundary,
                        material: &conductivity, element_materials: None,
                        source: if sample == stage { &q } else { &zero } }, interfaces: None,
                    radiation: EnthalpyRadiation::MaterialEnclosure(&material),
                    contacts: Some((&contact, config())),
                });
                let step = engine.advance_sdirk2(cx, inputs, &OLD, dt, config().step).unwrap();
                dot(&step.first_stage.specific_enthalpy_j_kg, &first)
                    + dot(&step.second_stage.specific_enthalpy_j_kg, &last)
            };
            let epsilon = 0.05;
            let fd = (-evaluate(2.0 * epsilon) + 8.0 * evaluate(epsilon)
                - 8.0 * evaluate(-epsilon) + evaluate(-2.0 * epsilon)) / (12.0 * epsilon);
            near(dot(&gradient.source_density[stage], &direction), fd, 2e-8);
        }
    });
}
