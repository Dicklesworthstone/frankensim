use super::*;
use fs_ascent::scan_design::constraints::{ScanControlBound, ScanEnergyCap};
use fs_ascent::scan_design::constraints::value::{ValueConstrainedScanStudy, ScanValueConstraintBudget};

#[test]
fn g3_value_constraints_match_sqp_units_order_and_fixed_heater_energy() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|_, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .adiabatic_remainder(), &Environment::earth_lab(), vec![], 0.).unwrap();
        for coordinates in [ScanPowerCoordinates::PositiveLog, ScanPowerCoordinates::NonnegativeLinear] {
            let make = || ScanPowerStudy::new_joint_with_power_coordinates(cx, sources([2., 3., 0.]),
                vec![controls()[0]], vec![ScanGeometryControl::HalfWidth { source: 0, axis: 0, reference_m: 1. }],
                coordinates, objective(observations(cx, &m), false), factory(&boundary, &k, &c), 3).unwrap();
            let bounds = vec![ScanControlBound { coordinate: 0, lower: None, upper: Some(1.5), scale: 2. },
                ScanControlBound { coordinate: 1, lower: Some(1.2), upper: Some(1.6), scale: 0.2 }];
            let caps = vec![ScanEnergyCap { start_time_s: 0., end_time_s: 2., maximum_j: 4., scale_j: 2. }];
            let states = vec![ScanStateBound { response: ScanResponse { endpoint: 0, vertex: 0,
                    observable: ScanObservable::SpecificEnthalpy }, lower: Some(150.), upper: Some(152.), scale: 2. },
                ScanStateBound { response: ScanResponse { endpoint: 1, vertex: 7,
                    observable: ScanObservable::Temperature }, lower: Some(309.), upper: Some(311.), scale: 1. }];
            let peaks = vec![ScanPeakBound { response: ScanPeakResponse { first_endpoint: 0, last_endpoint: 1,
                vertices: (0..8).collect(), observable: ScanObservable::LiquidMassFraction, smoothing_scale: 0.01 },
                maximum: 0.51, scale: 0.1 }];
            let mut valued = ValueConstrainedScanStudy::new(cx, make(), bounds.clone(), caps.clone(),
                states.clone(), peaks.clone(), ScanValueConstraintBudget { max_rows: 16, max_energy_weights: 1,
                    responses: ScanValueResponseBudget { max_responses: 3, max_peak_samples: 16 } }).unwrap();
            assert_eq!(valued.evaluations(), 0);
            let point = valued.initial_coordinates().to_vec();
            assert!(valued.constraints_at(cx, &point).is_err());
            assert_eq!(valued.evaluations(), 0);
            let (physical, values) = valued.evaluate_at(cx, &point).unwrap();
            let mut sqp = ConstrainedScanStudy::new_with_peak_bounds(cx, make(), bounds, caps, states, peaks,
                ScanResponseBudget { max_responses: 3, max_gradient_entries: 128, max_transpose_solves: 7 },
                16, 32).unwrap();
            let (differentiated, rows) = sqp.evaluate_at(cx, &point).unwrap();
            assert_eq!(physical.physical.value, differentiated.physical.value);
            assert_eq!(physical.physical.predictions, differentiated.physical.predictions);
            assert_eq!(values.nonnegative_power_values_w, rows.nonnegative_power_values_w);
            assert_eq!(values.bound_values, rows.bound_values);
            assert_eq!(values.absorbed_energy_j, rows.absorbed_energy_j);
            assert_eq!(values.state_values, rows.state_values);
            assert_eq!(values.peak_values, rows.peak_values);
            assert_eq!(values.sampled_peak_values, rows.sampled_peak_values);
            assert_eq!(values.peak_smoothing_gaps, rows.peak_smoothing_gaps);
            assert_eq!(values.residuals, rows.residuals);
            close(values.absorbed_energy_j[0], 5.); // Includes the unselected second heater value.
            assert!(!values.is_feasible(0.).unwrap());
            assert!(values.is_feasible(values.max_scaled_violation()).unwrap());
            for tolerance in [-1., f64::NAN, f64::INFINITY] { assert!(values.is_feasible(tolerance).is_err()); }
            assert_eq!(valued.evaluations(), 1);
        }
    });
}

#[test]
fn g1_phase_corner_feasibility_keeps_smoothing_margin_and_loss_separate() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|_, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .adiabatic_remainder(), &Environment::earth_lab(), vec![], 0.).unwrap();
        for (maximum, expected) in [(0.04, true), (0.02, false)] {
            let calls = std::cell::Cell::new(0usize);
            let create = |cx: &Cx<'_>, sources| {
                calls.set(calls.get() + 1);
                let run = EnthalpyHistoryRun::new(cx, ThermalHistoryProblem {
                    boundary: &boundary, material: &k, element_materials: None, interfaces: None,
                    source: HeatSourceHistory::Fixed(ScalarField::Uniform(0.)) },
                    EnthalpyHistoryMaterials { table: &[ReferenceEnthalpyMaterial {
                        curve: &c, reference_density_kg_m3: 1. }], of_element: &[0; 6], phase_conductivity: None },
                    &[100.; 8], EnthalpyHistoryBudget { spatial: EnthalpyBudget {
                        max_vertices: 8, max_elements: 6 }, max_steps: 2 })?;
                Ok(ScanningEnthalpyRun::new(cx, run, sources, 1, None)?)
            };
            let mut spec = objective(vec![ScanObservation { endpoint: 1, vertex: 0,
                observable: ScanObservable::LiquidMassFraction, target: 0.5, scale: 0.1 }], false);
            spec.budget.max_retained_vertices = 0;
            spec.budget.max_retained_source_vertices = 0;
            let study = ScanPowerStudy::new_joint_with_power_coordinates(cx, sources([0.; 3]), controls(),
                vec![], ScanPowerCoordinates::NonnegativeLinear, spec, create, 1).unwrap();
            let mut study = ValueConstrainedScanStudy::new(cx, study, vec![], vec![], vec![],
                vec![ScanPeakBound { response: ScanPeakResponse { first_endpoint: 0, last_endpoint: 1,
                    vertices: (0..8).collect(), observable: ScanObservable::LiquidMassFraction,
                    smoothing_scale: 0.01 }, maximum, scale: 0.1 }],
                ScanValueConstraintBudget { max_rows: 3, max_energy_weights: 0,
                    responses: ScanValueResponseBudget { max_responses: 1, max_peak_samples: 16 } }).unwrap();
            assert_eq!(calls.get(), 0);
            let (physical, constraints) = study.evaluate_at(cx, &[0., 0.]).unwrap();
            assert_eq!(calls.get(), 1);
            assert_eq!(physical.physical.value, 12.5); // Feasibility is NOT target recovery.
            assert_eq!(constraints.sampled_peak_values, vec![0.]);
            assert!((constraints.peak_values[0] - 0.01 * 16_f64.ln()).abs() < 1e-12);
            assert_eq!(constraints.is_feasible(0.).unwrap(), expected);
            assert!(study.evaluate_at(cx, &[0., 0.]).is_err());
            assert_eq!(calls.get(), 1);
            study.extend_evaluation_budget(2).unwrap();
            assert_eq!(study.evaluate_at(cx, &[0., 0.]).unwrap().1, constraints);
            assert_eq!(calls.get(), 2);
        }
    });
}

#[test]
fn g0_value_limit_admission_and_input_only_checks_do_not_spend_physical_evaluations() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|gate, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .adiabatic_remainder(), &Environment::earth_lab(), vec![], 0.).unwrap();
        let make = || ScanPowerStudy::new_joint_with_power_coordinates(cx, sources([2., 3., 0.]), controls(),
            vec![], ScanPowerCoordinates::NonnegativeLinear, objective(observations(cx, &m), false),
            factory(&boundary, &k, &c), 1).unwrap();
        let cap = ScanEnergyCap { start_time_s: 0., end_time_s: 2., maximum_j: 5., scale_j: 1. };
        let budget = ScanValueConstraintBudget { max_rows: 3, max_energy_weights: 2,
            responses: ScanValueResponseBudget { max_responses: 0, max_peak_samples: 0 } };
        assert!(ValueConstrainedScanStudy::new(cx, make(), vec![], vec![cap], vec![], vec![],
            ScanValueConstraintBudget { max_rows: 2, ..budget }).is_err());
        assert!(ValueConstrainedScanStudy::new(cx, make(), vec![], vec![cap], vec![], vec![],
            ScanValueConstraintBudget { max_energy_weights: 1, ..budget }).is_err());
        let bad = ScanControlBound { coordinate: 0, lower: Some(2.), upper: Some(2.), scale: 1. };
        assert!(ValueConstrainedScanStudy::new(cx, make(), vec![bad], vec![], vec![], vec![], budget).is_err());
        let mut study = ValueConstrainedScanStudy::new(cx, make(), vec![], vec![cap], vec![], vec![], budget).unwrap();
        let point = study.initial_coordinates().to_vec();
        let cheap = study.constraints_at(cx, &point).unwrap();
        assert_eq!(cheap.absorbed_energy_j, vec![5.]);
        assert!(cheap.is_feasible(0.).unwrap());
        assert_eq!(study.evaluations(), 0);
        assert!(study.evaluate_at(cx, &[-1., 1.]).is_err());
        assert_eq!(study.evaluations(), 0);
        assert_eq!(study.evaluate_at(cx, &point).unwrap().1, cheap);
        let exported = study.sources_at(cx, &point).unwrap();
        assert_eq!(exported[0].absorbed_power_values(16).unwrap()[2], 0.);
        study.extend_evaluation_budget(2).unwrap();
        gate.request();
        assert!(study.evaluate_at(cx, &point).is_err());
        assert_eq!(study.evaluations(), 1);
    });
}
