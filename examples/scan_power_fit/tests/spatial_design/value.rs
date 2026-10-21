use super::*;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::loss::{
    ScanObservationWithLoss, ScanResidualLoss};
use fs_ascent::scan_design::value::ScanValueResponseBudget;

mod feasibility;
mod batch;

#[test]
fn g3_value_only_matches_physical_scores_and_preserves_shared_work_accounting() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|_, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .remainder("surface", ThermalBc::adiabatic()).unwrap(),
            &Environment::earth_lab(), vec![], 0.).unwrap();
        for loss in [ScanResidualLoss::Quadratic, ScanResidualLoss::Huber { threshold: 0.5 }] {
            let rows: Vec<_> = observations(cx, &m).into_iter()
                .map(|o| ScanObservationWithLoss::new(o, loss).unwrap()).collect();
            let mut study = ScanPowerStudy::new(cx, sources([2., 3., 0.]), controls(),
                objective(rows, false), factory(&boundary, &k, &c), 2).unwrap();
            let point = study.initial_coordinates().to_vec();
            let valued = study.evaluate_value(cx, &point).unwrap();
            let differentiated = study.evaluate(cx, &point).unwrap();
            assert_eq!(valued.physical.value, differentiated.physical.value);
            assert_eq!(valued.physical.predictions, differentiated.physical.predictions);
            assert_eq!(valued.physical.standardized_residuals, differentiated.physical.standardized_residuals);
            assert_eq!(valued.physical.final_step.solution().specific_enthalpy_j_kg,
                differentiated.physical.final_step.solution().specific_enthalpy_j_kg);
            assert_eq!(valued.powers_w, differentiated.powers_w);
            assert_eq!(valued.physical.endpoint_count, 2);
            assert_eq!(study.evaluations(), 2);
            assert!(study.evaluate_value(cx, &point).is_err());
            assert_eq!(study.evaluations(), 2);
            study.extend_evaluation_budget(3).unwrap();
            assert_eq!(study.evaluate_value(cx, &point).unwrap().physical.value, valued.physical.value);
        }
    });
}

#[test]
fn g1_value_only_accepts_a_physical_phase_corner_without_fabricating_a_gradient() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|_, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .remainder("surface", ThermalBc::adiabatic()).unwrap(),
            &Environment::earth_lab(), vec![], 0.).unwrap();
        let corner_factory = |cx: &Cx<'_>, sources| {
            let run = EnthalpyHistoryRun::new(cx, ThermalHistoryProblem {
                boundary: &boundary, material: &k, element_materials: None, interfaces: None,
                source: HeatSourceHistory::Fixed(ScalarField::Uniform(0.)) },
                EnthalpyHistoryMaterials { table: &[ReferenceEnthalpyMaterial {
                    curve: &c, reference_density_kg_m3: 1. }], of_element: &[0; 6], phase_conductivity: None },
                &[100.; 8], EnthalpyHistoryBudget {
                    spatial: EnthalpyBudget { max_vertices: 8, max_elements: 6 }, max_steps: 2 })?;
            Ok(ScanningEnthalpyRun::new(cx, run, sources, 1, None)?)
        };
        let rows = vec![ScanObservation { endpoint: 1, vertex: 0,
            observable: ScanObservable::LiquidMassFraction, target: 0.5, scale: 0.1 }];
        let mut study = ScanPowerStudy::new_joint_with_power_coordinates(cx, sources([0.; 3]),
            controls(), vec![], ScanPowerCoordinates::NonnegativeLinear,
            objective(rows, false), corner_factory, 3).unwrap();
        let point = study.initial_coordinates().to_vec();
        let valued = study.evaluate_value(cx, &point).unwrap();
        assert_eq!(valued.physical.predictions, vec![0.]);
        assert_eq!(valued.physical.value, 12.5);
        assert!(valued.physical.final_step.solution().temperature.iter().all(|&t| t == 310.));
        assert!(study.evaluate(cx, &point).is_err());
        assert_eq!(study.evaluations(), 2);
        assert_eq!(study.evaluate_value(cx, &point).unwrap().physical.value, 12.5);
    });
}

#[test]
fn g0_value_only_needs_no_tangent_budget_but_preserves_input_and_cancel_refusals() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|gate, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .remainder("surface", ThermalBc::adiabatic()).unwrap(),
            &Environment::earth_lab(), vec![], 0.).unwrap();
        let mut spec = objective(observations(cx, &m), false);
        spec.budget.max_retained_vertices = 0;
        spec.budget.max_retained_source_vertices = 0;
        let mut study = ScanPowerStudy::new(cx, sources([2., 3., 0.]), controls(), spec,
            factory(&boundary, &k, &c), 3).unwrap();
        let point = study.initial_coordinates().to_vec();
        assert!(study.evaluate_value(cx, &point).unwrap().physical.value.is_finite());
        assert!(study.evaluate(cx, &point).is_err());
        assert!(study.evaluate_value(cx, &[f64::NAN, 0.]).is_err());
        assert_eq!(study.evaluations(), 2);
        gate.request();
        assert!(study.evaluate_value(cx, &point).is_err());
        assert_eq!(study.evaluations(), 2);
    });
}

#[test]
fn g3_value_responses_match_all_three_physical_adjoints_and_do_not_follow_loss_weights() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|_, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .remainder("surface", ThermalBc::adiabatic()).unwrap(),
            &Environment::earth_lab(), vec![], 0.).unwrap();
        let quantities = [ScanObservable::Temperature, ScanObservable::SpecificEnthalpy,
            ScanObservable::LiquidMassFraction];
        let points: Vec<_> = quantities.iter().enumerate().map(|(i, &observable)|
            ScanResponse { endpoint: usize::from(i != 0), vertex: i, observable }).collect();
        let peaks: Vec<_> = quantities.iter().map(|&observable| ScanPeakResponse {
            first_endpoint: 0, last_endpoint: 1, vertices: (0..8).collect(), observable,
            smoothing_scale: if observable == ScanObservable::LiquidMassFraction { 0.01 } else { 0.2 },
        }).collect();
        let mut compared = Vec::new();
        for loss in [ScanResidualLoss::Quadratic, ScanResidualLoss::Huber { threshold: 0.1 }] {
            let rows: Vec<_> = observations(cx, &m).into_iter()
                .map(|row| ScanObservationWithLoss::new(row, loss).unwrap()).collect();
            let direct = objective(rows.clone(), false);
            let mut study = ScanPowerStudy::new(cx, sources([2., 3., 0.]), controls(),
                objective(rows.clone(), false), factory(&boundary, &k, &c), 1).unwrap();
            let point = study.initial_coordinates().to_vec();
            let controlled = study.sources_at(cx, &point).unwrap();
            let (value, actual_points, actual_peaks) = study.evaluate_value_with_responses(cx, &point,
                &points, &peaks, ScanValueResponseBudget { max_responses: 6, max_peak_samples: 48 }).unwrap();
            let (gradient, expected_points, expected_peaks, _) = evaluate_scan_observations(cx,
                factory(&boundary, &k, &c)(cx, controlled).unwrap(), &direct.endpoints_s, &rows,
                direct.sampling, direct.step, direct.transpose, direct.budget, &points, &peaks,
                ScanResponseBudget { max_responses: 6, max_gradient_entries: 256, max_transpose_solves: 13 },
                48, None).unwrap();
            assert_eq!(value.physical.value, gradient.value);
            assert_eq!(value.physical.predictions, gradient.predictions);
            assert_eq!(value.physical.standardized_residuals, gradient.standardized_residuals);
            assert_eq!(value.physical.endpoint_count, 2);
            assert_eq!(study.evaluations(), 1);
            for (actual, expected) in actual_points.iter().zip(expected_points) {
                assert_eq!(actual.response, expected.response);
                assert_eq!(actual.value, expected.value);
            }
            for (actual, expected) in actual_peaks.iter().zip(expected_peaks) {
                assert_eq!(actual.peak, expected.peak);
                assert_eq!(actual.value, expected.value);
                assert_eq!(actual.sampled_maximum, expected.sampled_maximum);
                assert_eq!(actual.smoothing_gap, expected.smoothing_gap);
                assert_eq!(actual.maximum_smoothing_gap, expected.maximum_smoothing_gap);
                assert_eq!(actual.sample_count, expected.sample_count);
            }
            compared.push((actual_points, actual_peaks));
        }
        assert_eq!(compared[0], compared[1]);
    });
}

#[test]
fn g1_value_responses_at_solidus_need_no_tangent_storage_or_phase_derivative() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|_, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .adiabatic_remainder(), &Environment::earth_lab(), vec![], 0.).unwrap();
        let corner_factory = |cx: &Cx<'_>, sources| {
            let run = EnthalpyHistoryRun::new(cx, ThermalHistoryProblem {
                boundary: &boundary, material: &k, element_materials: None, interfaces: None,
                source: HeatSourceHistory::Fixed(ScalarField::Uniform(0.)) },
                EnthalpyHistoryMaterials { table: &[ReferenceEnthalpyMaterial {
                    curve: &c, reference_density_kg_m3: 1. }], of_element: &[0; 6], phase_conductivity: None },
                &[100.; 8], EnthalpyHistoryBudget { spatial: EnthalpyBudget {
                    max_vertices: 8, max_elements: 6 }, max_steps: 2 })?;
            Ok(ScanningEnthalpyRun::new(cx, run, sources, 1, None)?)
        };
        let mut spec = objective(observations(cx, &m), false);
        spec.budget.max_retained_vertices = 0;
        spec.budget.max_retained_source_vertices = 0;
        let mut study = ScanPowerStudy::new_joint_with_power_coordinates(cx, sources([0.; 3]),
            controls(), vec![], ScanPowerCoordinates::NonnegativeLinear, spec, corner_factory, 1).unwrap();
        let points = [ScanResponse { endpoint: 0, vertex: 0, observable: ScanObservable::SpecificEnthalpy },
            ScanResponse { endpoint: 1, vertex: 7, observable: ScanObservable::LiquidMassFraction }];
        let peak = ScanPeakResponse { first_endpoint: 0, last_endpoint: 1, vertices: (0..8).collect(),
            observable: ScanObservable::LiquidMassFraction, smoothing_scale: 0.01 };
        let (_, values, peaks) = study.evaluate_value_with_responses(cx, &[0., 0.], &points, &[peak],
            ScanValueResponseBudget { max_responses: 3, max_peak_samples: 16 }).unwrap();
        assert_eq!(values[0].value, 100.);
        assert_eq!(values[1].value, 0.);
        assert_eq!(peaks[0].sampled_maximum, 0.);
        assert!((peaks[0].value - 0.01 * 16_f64.ln()).abs() < 1e-12);
        assert_eq!(peaks[0].value, peaks[0].smoothing_gap);
    });
}

#[test]
fn g0_value_responses_reject_bad_selection_caps_and_cancellation_without_a_partial_report() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|gate, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .adiabatic_remainder(), &Environment::earth_lab(), vec![], 0.).unwrap();
        let mut study = ScanPowerStudy::new(cx, sources([2., 3., 0.]), controls(),
            objective(observations(cx, &m), false), factory(&boundary, &k, &c), 5).unwrap();
        let point = study.initial_coordinates().to_vec();
        for (endpoint, vertex) in [(2, 0), (0, 8)] {
            let response = ScanResponse { endpoint, vertex, observable: ScanObservable::Temperature };
            assert!(study.evaluate_value_with_responses(cx, &point, &[response], &[],
                ScanValueResponseBudget { max_responses: 1, max_peak_samples: 0 }).is_err());
        }
        let peak = ScanPeakResponse { first_endpoint: 0, last_endpoint: 1, vertices: (0..8).collect(),
            observable: ScanObservable::SpecificEnthalpy, smoothing_scale: 0.1 };
        for budget in [ScanValueResponseBudget { max_responses: 0, max_peak_samples: 16 },
            ScanValueResponseBudget { max_responses: 1, max_peak_samples: 15 }] {
            assert!(study.evaluate_value_with_responses(cx, &point, &[], &[peak.clone()], budget).is_err());
        }
        assert_eq!(study.evaluations(), 4);
        gate.request();
        assert!(study.evaluate_value_with_responses(cx, &point, &[], &[peak],
            ScanValueResponseBudget { max_responses: 1, max_peak_samples: 16 }).is_err());
        assert_eq!(study.evaluations(), 4);
    });
}
