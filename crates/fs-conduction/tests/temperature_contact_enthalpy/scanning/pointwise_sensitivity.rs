//! Local radiation and live contacts retain their original physical scan dual.
use super::*;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::{
    ScanObjectiveBudget, ScanObjectiveSampling, ScanObservable, ScanObservation,
    checkpoint::{ScanCheckpointBudget, evaluate_scan_objective_checkpointed_with_peaks},
    responses::{ScanResponse, ScanResponseBudget,
        peaks::{ScanPeakResponse, evaluate_scan_objective_with_peaks}},
    value::{ScanValueBudget,
        responses::{ScanValueResponseBudget, evaluate_scan_values_with_responses}},
};

const P: [f64; 7] = [0.15, 0.2, 0.3, 1.6, 1.5, 1.2, 3.];

fn patches(bath: TimeSignal, epsilon: f64) -> Vec<RadiationRegionHistory> {
    let dims = fs_conduction::EMISSIVITY_DIMS;
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(fs_conduction::SURFACE_EMISSIVITY_PROPERTY, dims),
        value: PropertyValue::Scalar { value: epsilon, dims },
        validity: ValidityDomain::unconstrained().with("T", 280., 380.),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance { source: "synthetic pointwise scan controls".into(),
            license: "internal-test-use".into(), artifact: None },
    }).unwrap();
    let card = fs_matdb::MaterialCard::assemble(MaterialStateId {
        chemistry: "pointwise scan surface".into(), phase: "solid".into(),
        process: "fixture".into(), revision: 0,
    }, claims, vec![]).unwrap();
    vec![RadiationRegionHistory { region: "exterior".into(),
        emissivity: SurfaceEmissivity::from_card("exterior", &card, 300.,
            SelectionPolicy::SingleClaimOnly).unwrap(), ambient_temperature: bath }]
}

#[allow(clippy::too_many_arguments)]
fn scan<'h, 'm, 'c>(cx: &Cx<'_>, history: &'h ThermalBoundaryHistory<'m>,
    material: &'h ConductivityModel, curves: &'c [EquilibriumEnthalpyPhaseCurve; 2],
    card: Option<&'c InterfaceSystemCard>, p: [f64; 7], bath: TimeSignal,
    epsilon: f64) -> ScanningEnthalpyRun<'h, 'm, 'c>
{
    let source = heater_with_geometry(
        [p[0], p[1], p[2]].map(|v| constant(v, LENGTH)),
        [p[3], p[4], p[5]], constant(p[6], COMPONENT_POWER_DIMS));
    let run = ScanningEnthalpyRun::new(cx, base(cx, history, material, curves, 2),
        vec![source], 1, None).unwrap()
        .with_pointwise_ambient_radiation(cx, patches(bath, epsilon)).unwrap();
    if let Some(card) = card {
        let initial = history.at(cx, 0.).unwrap();
        let contacts = TemperatureContacts::new(cx, history.mesh(), &initial,
            vec![], vec![law(history.mesh(), card)], 1).unwrap();
        run.with_contacts_newton(cx, contacts, coupled_newton()).unwrap()
    } else { run }
}

#[test]
fn g3_pointwise_scan_sources_and_radiation_controls_share_the_native_contact_dual() {
    context(&CancelGate::new_clock_free(), |cx| {
        let card = card();
        let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(280., 1.6), (380., 3.6)]).unwrap());
        let seed = [0.3, -0.1, 0.7, 0.5, -0.2, 0.4, 0.1, -0.3];
        for coupled in [false, true] {
            let mesh = if coupled { mesh() } else { separated_mesh() };
            let history = control_history(cx, &mesh);
            let mut run = scan(cx, &history, &material, &curves, coupled.then_some(&card),
                P, constant(300., TEMPERATURE_DIMS), 0.2);
            let endpoint = run.advance_linearized_to(cx, 0.025, coupled_newton().step,
                control_budget()).unwrap();
            assert!(endpoint.report().solution().pointwise_radiation.is_some());
            let gradient = endpoint.pullback_temperature(cx, &seed, control_linear()).unwrap();
            let radiation = gradient.radiation.as_ref().unwrap();
            assert_eq!(radiation.radiation_regions, vec!["exterior".to_string()]);
            assert_eq!(radiation.reservoir_temperatures.len(), 1);
            assert_eq!(radiation.emissivities.len(), 1);
            assert_eq!(gradient.contact_resistance_multipliers.len(), usize::from(coupled));
            let evaluate = |p: [f64; 7], bath: f64, epsilon: f64, resistance: f64| {
                let trial_card = card_with_curve(vec![(280., 0.1 * resistance), (380., 0.7 * resistance)]);
                let mut trial = scan(cx, &history, &material, &curves,
                    coupled.then_some(&trial_card), p, constant(bath, TEMPERATURE_DIMS), epsilon);
                let report = trial.advance_to(cx, 0.025, coupled_newton().step).unwrap();
                objective(&seed, &report.solution().temperature)
            };
            close(evaluate(P, 300., 0.2, 1.),
                objective(&seed, &endpoint.report().solution().temperature), 1e-12);
            let source = heater_controls(&gradient.sources[0]);
            for coordinate in 0..7 {
                let (mut plus, mut minus) = (P, P);
                plus[coordinate] += 1e-4; minus[coordinate] -= 1e-4;
                close_control(source[coordinate],
                    (evaluate(plus, 300., 0.2, 1.) - evaluate(minus, 300., 0.2, 1.)) / 2e-4);
            }
            close_control(radiation.reservoir_temperatures[0],
                (evaluate(P, 300.001, 0.2, 1.) - evaluate(P, 299.999, 0.2, 1.)) / 0.002);
            close_control(radiation.emissivities[0],
                (evaluate(P, 300., 0.20001, 1.) - evaluate(P, 300., 0.19999, 1.)) / 2e-5);
            if coupled {
                let resistance = &gradient.contact_resistance_multipliers[0];
                assert_eq!(resistance.interface, "joint");
                assert_eq!(resistance.card_identity, card.content_hash());
                close_control(resistance.derivative,
                    (evaluate(P, 300., 0.2, 1.0001) - evaluate(P, 300., 0.2, 0.9999)) / 2e-4);
            }
        }
    });
}

#[test]
fn g3_pointwise_reservoir_histories_retain_event_sides_and_checkpointed_contact_feedback() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let history = control_history(cx, &mesh); let card = card();
        let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(280., 1.6), (380., 3.6)]).unwrap());
        let times = [0.03, 0.08]; let bath = [298., 302., 320.];
        let make = |values: &[f64], epsilon| scan(cx, &history, &material, &curves,
            Some(&card), P, TimeSignal::Table { times: vec![0., 0.03, 0.08],
                values: values.to_vec(), dims: TEMPERATURE_DIMS, interp: Interp::Hold }, epsilon);
        let observations = [ScanObservation { endpoint: 1, vertex: 1,
            observable: ScanObservable::SpecificEnthalpy, target: 110., scale: 8. }];
        let points = [ScanResponse { endpoint: 0, vertex: 0, observable: ScanObservable::Temperature }];
        let peaks = [ScanPeakResponse { first_endpoint: 0, last_endpoint: 1,
            vertices: vec![0, 1], observable: ScanObservable::SpecificEnthalpy, smoothing_scale: 2. }];
        // Seven heater, three bath, one emissivity and one whole-contact-law control.
        let budget = ScanObjectiveBudget { max_steps: 2, max_observations: 1,
            max_retained_vertices: 16, max_retained_source_vertices: 16,
            max_schedule_controls: 12, max_retained_radiation_entries: 2048 };
        let responses = ScanResponseBudget { max_responses: 2, max_gradient_entries: 40,
            max_transpose_solves: 5 };
        let (result, rows, envelopes) = evaluate_scan_objective_with_peaks(cx,
            make(&bath, 0.2), &times, &observations, ScanObjectiveSampling::EventLeft,
            coupled_newton().step, control_linear(), budget, &points, &peaks, responses, 4).unwrap();
        let replay = evaluate_scan_objective_checkpointed_with_peaks(cx,
            make(&bath, 0.2), &times, &observations, ScanObjectiveSampling::EventLeft,
            coupled_newton().step, control_linear(), ScanObjectiveBudget {
                max_retained_vertices: 8, max_retained_source_vertices: 8,
                max_retained_radiation_entries: 1024, ..budget }, &points, &peaks, responses, 4,
            ScanCheckpointBudget { segment_steps: 1, max_checkpoint_vertices: 24, max_forward_solves: 4 }).unwrap();
        assert_eq!(result.radiation, replay.objective.radiation);
        assert_eq!(rows[0].radiation, replay.responses[0].radiation);
        assert_eq!(envelopes[0].radiation, replay.peaks[0].radiation);
        assert_eq!(result.contact_resistance_multipliers, replay.objective.contact_resistance_multipliers);
        assert_eq!(result.sources, replay.objective.sources);
        assert_eq!(result.initial_specific_enthalpy, replay.objective.initial_specific_enthalpy);
        assert_eq!(result.transpose_iterations, replay.objective.transpose_iterations);
        assert_eq!(replay.work.recomputed_forward_solves, 2);
        let gradients = [&result.radiation[0], &rows[0].radiation[0], &envelopes[0].radiation[0]];
        for gradient in gradients { assert_eq!(gradient.reservoir_values_k[2], 0.); }
        assert_eq!(rows[0].radiation[0].reservoir_values_k[1], 0.);
        let evaluate = |values: &[f64], epsilon| {
            let (objective, rows, peaks) = evaluate_scan_values_with_responses(cx,
                make(values, epsilon), &times, &observations, ScanObjectiveSampling::EventLeft,
                coupled_newton().step, ScanValueBudget { max_steps: 2, max_observations: 1,
                    max_schedule_controls: 7 }, &points, &peaks,
                ScanValueResponseBudget { max_responses: 2, max_peak_samples: 4 }).unwrap();
            [objective.value, rows[0].value, peaks[0].value]
        };
        assert_eq!(evaluate(&bath, 0.2), [result.value, rows[0].value, envelopes[0].value]);
        for coordinate in 0..4 {
            let (mut plus, mut minus) = (bath, bath);
            let (ep, em, delta) = if coordinate == 3 { (0.20001, 0.19999, 1e-5) }
                else { plus[coordinate] += 1e-3; minus[coordinate] -= 1e-3; (0.2, 0.2, 1e-3) };
            let (a, b) = (evaluate(&plus, ep), evaluate(&minus, em));
            for row in 0..3 {
                let actual = if coordinate == 3 { gradients[row].emissivity }
                    else { gradients[row].reservoir_values_k[coordinate] };
                close_control(actual, (a[row] - b[row]) / (2. * delta));
            }
        }
    });
}

#[test]
fn g4_pointwise_binding_keeps_native_relative_target_and_fixed_row_control_masks() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = separated_mesh(); let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(280., 1.6), (380., 3.6)]).unwrap());
        let history = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
            .region("bath", |face| face.centroid[2] > 0. && face.centroid[1] == 0.,
                ThermalBc::dirichlet(340.).unwrap()).unwrap()
            .region("exterior", |face| face.centroid[2] < 0. && face.centroid[1] == 0.,
                ThermalBc::neumann(0.).unwrap()).unwrap().adiabatic_remainder(),
            &Environment::earth_lab(), Vec::new(), 0.).unwrap();
        let mut run = scan(cx, &history, &material, &curves, None, P,
            constant(300., TEMPERATURE_DIMS), 0.2);
        let mut step = config();
        step.newton.absolute_tolerance = 0.;
        step.newton.relative_tolerance = 1e-10;
        // Binding may prepare a physical candidate, but no state is published
        // when its explicit radiation tangent allowance is insufficient.
        assert!(run.advance_linearized_to(cx, 0.025, step, ScanSensitivityBudget {
            max_radiation_feedback_entries: 0, ..control_budget() }).is_err());
        assert_eq!(run.conduction().specific_enthalpy_j_kg(), &OLD);
        assert_eq!(run.conduction().time_s(), 0.);
        assert_eq!(run.conduction().accepted_steps(), 0);
        let endpoint = run.advance_linearized_to(cx, 0.025, step, control_budget()).unwrap();
        let mut fixed_seed = [0.; 8]; fixed_seed[4] = 1.;
        let gradient = endpoint.pullback_temperature(cx, &fixed_seed, control_linear()).unwrap();
        close(gradient.storage.transport.prescribed_temperature[4], 1., 1e-10);
        for value in heater_controls(&gradient.sources[0]) { close(value, 0., 1e-12); }
        for value in &gradient.storage.transport.previous_specific_enthalpy { close(*value, 0., 1e-12); }
        let radiation = gradient.radiation.unwrap();
        close(radiation.reservoir_temperatures[0], 0., 1e-12);
        close(radiation.emissivities[0], 0., 1e-12);
        assert_eq!(run.conduction().accepted_steps(), 1);
    });
}
