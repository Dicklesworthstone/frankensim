use super::*;
use fs_ascent::scan_design::constraints::{ScanPowerSlewBound, ScanEnergyCap};
use fs_ascent::scan_design::constraints::value::{ScanValueConstraintBudget, ValueConstrainedScanStudy};
use fs_dfo::differential_evolution::DifferentialEvolutionConfig;
use scan_power_fit::population::PopulationScanStudy;

fn linear_sources(values: [f64; 3]) -> Vec<MovingHeatSource> {
    let length = Dims([1, 0, 0, 0, 0, 0]);
    vec![MovingHeatSource::new(std::array::from_fn(|_|
        TimeSignal::Constant(QtyAny::new(0.5, length))), [1.4, 1.2, 1.1],
        TimeSignal::Table { times: vec![0.0, 1.0, 2.0], values: values.to_vec(),
            dims: COMPONENT_POWER_DIMS, interp: Interp::Linear }, 0.0,
        MovingHeatBudget { max_vertices: 8, max_elements: 6, min_active_vertices: 1,
            power_relative_tolerance: 1e-11 }).unwrap()]
}
fn bound() -> ScanPowerSlewBound {
    ScanPowerSlewBound { source: 0, maximum_rise_w_s: 2.0,
        maximum_fall_w_s: 2.0, scale_w_s: 1.0 }
}
fn budget() -> ScanValueConstraintBudget {
    ScanValueConstraintBudget { max_rows: 16, max_energy_weights: 4,
        responses: ScanValueResponseBudget { max_responses: 0, max_peak_samples: 0 } }
}
fn no_responses() -> ScanResponseBudget {
    ScanResponseBudget { max_responses: 0, max_gradient_entries: 0, max_transpose_solves: 0 }
}

#[test]
fn g3_slew_rows_and_coordinate_gradients_match_serial_batch_and_sqp() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.0).unwrap(); let c = curve();
    with_cx(|_, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .remainder("surface", ThermalBc::adiabatic()).unwrap(),
            &Environment::earth_lab(), vec![], 0.0).unwrap();
        for map in [ScanPowerCoordinates::PositiveLog, ScanPowerCoordinates::NonnegativeLinear] {
            let make = || ScanPowerStudy::new_joint_with_power_coordinates(cx, linear_sources([2.0, 3.0, 0.0]),
                controls(), vec![ScanGeometryControl::HalfWidth { source: 0, axis: 0, reference_m: 1.0 }],
                map, objective(observations(cx, &m), false), factory(&boundary, &k, &c), 20).unwrap();
            let cap = ScanEnergyCap { start_time_s: 0.0, end_time_s: 2.0, maximum_j: 100.0, scale_j: 10.0 };
            let mut values = ValueConstrainedScanStudy::new_with_power_slew_bounds(cx, make(),
                vec![], vec![cap], vec![], vec![], vec![bound()], 2, budget()).unwrap();
            let point = values.initial_coordinates().to_vec();
            assert_eq!(values.evaluations(), 0);
            let mut sqp = ConstrainedScanStudy::new_with_power_slew_bounds(cx, make(),
                vec![], vec![cap], vec![], vec![], no_responses(), 0, vec![bound()], 2, 32).unwrap();
            let (physical, rows) = values.evaluate_at(cx, &point).unwrap();
            let (derivative_physical, differentiated) = sqp.evaluate_at(cx, &point).unwrap();
            assert_eq!(rows.residuals, differentiated.residuals);
            assert_eq!(rows.power_slew_values_w_s, vec![1.0, -3.0]);
            assert_eq!(rows.power_slew_values_w_s, differentiated.power_slew_values_w_s);
            assert_eq!(physical.physical.predictions, derivative_physical.physical.predictions);
            assert!(!rows.is_feasible(0.0).unwrap());
            assert_eq!(rows.max_scaled_violation(), 1.0);
            assert_eq!(rows.absorbed_energy_j, vec![4.0]);
            for j in 0..point.len() {
                let (mut plus, mut minus) = (point.clone(), point.clone());
                plus[j] += 1e-5; minus[j] -= 1e-5;
                let p = sqp.constraints_at(cx, &plus).unwrap();
                let q = sqp.constraints_at(cx, &minus).unwrap();
                for row in 0..rows.residuals.len() {
                    close(differentiated.jacobian[row * point.len() + j],
                        (p.residuals[row] - q.residuals[row]) / 2e-5);
                }
            }
            let spent = values.evaluations();
            {
                // A FnMut-only fixture cannot execute a shared batch; the same
                // immutable rate rows are nevertheless exercised by the worker
                // path in the dedicated Fn factory test below.
                let input_rows = values.constraints_at(cx, &point).unwrap();
                assert_eq!(rows, input_rows);
            }
            assert_eq!(values.evaluations(), spent);
            assert_eq!(values.sources_at(cx, &point).unwrap()[0].absorbed_power_values(3).unwrap()[2], 0.0);
        }
    });
}

#[test]
fn g3_fixed_heaters_keep_their_violations_and_shared_batches_keep_all_rate_rows() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.0).unwrap(); let c = curve();
    with_cx(|_, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .remainder("surface", ThermalBc::adiabatic()).unwrap(),
            &Environment::earth_lab(), vec![], 0.0).unwrap();
        let mut templates = linear_sources([2.0, 3.0, 0.0]);
        templates.extend(linear_sources([5.0, 9.0, 0.0]));
        let mut spec = objective(observations(cx, &m), false);
        spec.budget.max_schedule_controls = 32;
        let make = |cx: &Cx<'_>, sources| -> Result<_, ScanSensitivityError> {
            let run = EnthalpyHistoryRun::new(cx, ThermalHistoryProblem {
                boundary: &boundary, material: &k, element_materials: None, interfaces: None,
                source: HeatSourceHistory::Fixed(ScalarField::Uniform(0.0)) },
                EnthalpyHistoryMaterials { table: &[ReferenceEnthalpyMaterial {
                    curve: &c, reference_density_kg_m3: 1.0 }], of_element: &[0; 6], phase_conductivity: None },
                &[150.0; 8], EnthalpyHistoryBudget {
                    spatial: EnthalpyBudget { max_vertices: 8, max_elements: 6 }, max_steps: 2 })?;
            Ok(ScanningEnthalpyRun::new(cx, run, sources, 2, None)?)
        };
        let study = ScanPowerStudy::new_joint_with_power_coordinates(cx, templates, controls(), vec![],
            ScanPowerCoordinates::NonnegativeLinear, spec, make, 8).unwrap();
        let mut limits = ValueConstrainedScanStudy::new_with_power_slew_bounds(cx, study,
            vec![], vec![], vec![], vec![], vec![ScanPowerSlewBound { source: 1, ..bound() }],
            2, budget()).unwrap();
        let point = limits.initial_coordinates().to_vec();
        let (_, original) = limits.evaluate_at(cx, &point).unwrap();
        assert_eq!(original.power_slew_values_w_s, vec![4.0, -9.0]);
        assert_eq!(original.max_scaled_violation(), 7.0);
        {
            let batch = limits.value_batch(cx, 2).unwrap();
            let (_, same) = batch.evaluate_at(cx, &point).unwrap();
            assert_eq!(same, original);
            let (_, changed) = batch.evaluate_at(cx, &[2.0, 2.0]).unwrap();
            assert_eq!(changed.power_slew_values_w_s, original.power_slew_values_w_s);
            assert_eq!(batch.started_evaluations(), 2);
        }
        assert_eq!(limits.evaluations(), 3);
    });
}

#[test]
fn g0_rate_admission_refuses_held_signals_bad_limits_and_excess_rows_before_physics() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.0).unwrap(); let c = curve();
    with_cx(|gate, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .remainder("surface", ThermalBc::adiabatic()).unwrap(),
            &Environment::earth_lab(), vec![], 0.0).unwrap();
        let calls = std::cell::Cell::new(0);
        let make = |held| ScanPowerStudy::new_joint_with_power_coordinates(cx,
            if held { sources([2.0, 3.0, 0.0]) } else { linear_sources([2.0, 3.0, 0.0]) },
            controls(), vec![], ScanPowerCoordinates::NonnegativeLinear,
            objective(observations(cx, &m), false), |cx: &Cx<'_>, sources| {
                calls.set(calls.get() + 1); factory(&boundary, &k, &c)(cx, sources)
            }, 10).unwrap();
        for (held, declarations, segments, max_rows) in [
            (true, vec![bound()], 2, 16),
            (false, vec![bound(), bound()], 4, 16),
            (false, vec![ScanPowerSlewBound { source: 1, ..bound() }], 2, 16),
            (false, vec![ScanPowerSlewBound { maximum_rise_w_s: -1.0, ..bound() }], 2, 16),
            (false, vec![ScanPowerSlewBound { scale_w_s: f64::NAN, ..bound() }], 2, 16),
            (false, vec![bound()], 1, 16),
            (false, vec![bound()], 2, 5),
        ] {
            assert!(ValueConstrainedScanStudy::new_with_power_slew_bounds(cx, make(held),
                vec![], vec![], vec![], vec![], declarations, segments,
                ScanValueConstraintBudget { max_rows, ..budget() }).is_err());
        }
        assert_eq!(calls.get(), 0);
        let study = make(false);
        gate.request();
        assert!(ValueConstrainedScanStudy::new_with_power_slew_bounds(cx, study,
            vec![], vec![], vec![], vec![], vec![bound()], 2, budget()).is_err());
        assert_eq!(calls.get(), 0);
    });
}

#[test]
fn g1_population_fit_enforces_the_fixed_shutdown_rate_not_just_power_or_energy() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.0).unwrap(); let c = curve();
    with_cx(|_, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .remainder("surface", ThermalBc::adiabatic()).unwrap(),
            &Environment::earth_lab(), vec![], 0.0).unwrap();
        let observations = (0..2).map(|endpoint| ScanObservation { endpoint, vertex: 0,
            observable: ScanObservable::SpecificEnthalpy, target: 160.0, scale: 10.0 }).collect();
        let study = ScanPowerStudy::new_joint_with_power_coordinates(cx, linear_sources([0.0, 1.0, 0.0]),
            vec![ScanPowerControl { source: 0, value: 1, reference_w: 1.0 }], vec![],
            ScanPowerCoordinates::NonnegativeLinear, objective(observations, false),
            factory(&boundary, &k, &c), 514).unwrap();
        let limits = ValueConstrainedScanStudy::new_with_power_slew_bounds(cx, study,
            vec![], vec![], vec![], vec![], vec![ScanPowerSlewBound {
                maximum_rise_w_s: 3.0, maximum_fall_w_s: 2.0, ..bound() }], 2, budget()).unwrap();
        let mut search = PopulationScanStudy::new(cx, limits, &[[0.0, 20.0]],
            DifferentialEvolutionConfig { population_size: 8, differential_weight: 0.8,
                crossover_probability: 0.9, seed: 17, max_evaluations: 512,
                max_population_entries: 8, target: Some(1e-12) }, 1e-10).unwrap();
        search.advance(cx, 512).unwrap();
        let replay = search.replay_best(cx).unwrap();
        assert!(replay.score.is_feasible());
        assert!(!replay.target_recovered());
        assert!((replay.physical.powers_w[0] - 2.0).abs() < 1e-3);
        assert!((replay.physical.physical.value - 0.64).abs() < 2e-4);
        assert!(replay.constraints.power_slew_values_w_s[1] >= -2.0 - 1e-10);
        assert_eq!(replay.sources[0].absorbed_power_values(3).unwrap()[2], 0.0);
    });
}
