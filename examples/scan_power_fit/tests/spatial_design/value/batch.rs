use super::*;
use fs_ascent::scan_design::constraints::ScanEnergyCap;
use fs_ascent::scan_design::constraints::value::{
    ScanConstraintValues, ScanValueConstraintBudget, ValueConstrainedScanStudy,
};
use fs_ascent::scan_design::value::ScanPowerValueEvaluation;
use std::sync::atomic::{AtomicBool, Ordering};

fn study<'a>(cx: &Cx<'_>, mesh: &'a ConductionMesh,
    boundary: &'a ThermalBoundaryHistory<'a>, material: &'a ConductivityModel,
    chart: &'a EquilibriumEnthalpyPhaseCurve, fail: &'a AtomicBool, limit: usize)
    -> ValueConstrainedScanStudy<impl Fn(&Cx<'_>, Vec<MovingHeatSource>)
        -> Result<ScanningEnthalpyRun<'a, 'a, 'a>, ScanSensitivityError> + Sync + 'a,
        ScanObservationWithLoss<ScanSpatialObservation<'a>>>
{
    let factory = move |worker: &Cx<'_>, sources| {
        if fail.swap(false, Ordering::Relaxed) {
            return Err(fs_conduction::ConductionError::Config {
                parameter: "injected batch factory failure", what: "retry the same declaration".into(),
            }.into());
        }
        super::factory(boundary, material, chart)(worker, sources)
    };
    let rows = observations(cx, mesh).into_iter().map(|row|
        ScanObservationWithLoss::new(row, ScanResidualLoss::Huber { threshold: 0.1 }).unwrap()).collect();
    let mut spec = objective(rows, false);
    spec.budget.max_retained_vertices = 0;
    spec.budget.max_retained_source_vertices = 0;
    let study = ScanPowerStudy::new_joint_with_power_coordinates(cx, sources([0.; 3]), controls(),
        vec![], ScanPowerCoordinates::NonnegativeLinear, spec, factory, limit).unwrap();
    ValueConstrainedScanStudy::new(cx, study, vec![],
        vec![ScanEnergyCap { start_time_s: 0., end_time_s: 2., maximum_j: 10., scale_j: 1. }],
        vec![ScanStateBound { response: ScanResponse { endpoint: 0, vertex: 0,
            observable: ScanObservable::SpecificEnthalpy }, lower: None, upper: Some(158.), scale: 10. }],
        vec![ScanPeakBound { response: ScanPeakResponse { first_endpoint: 0, last_endpoint: 1,
            vertices: (0..8).collect(), observable: ScanObservable::LiquidMassFraction,
            smoothing_scale: 0.01 }, maximum: 0.8, scale: 0.1 }],
        ScanValueConstraintBudget { max_rows: 5, max_energy_weights: 2,
            responses: ScanValueResponseBudget { max_responses: 2, max_peak_samples: 16 } }).unwrap()
}

fn signature((physical, constraints): (ScanPowerValueEvaluation, ScanConstraintValues))
    -> (u64, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, ScanConstraintValues)
{
    let solution = physical.physical.final_step.solution();
    (physical.physical.value.to_bits(), physical.physical.predictions.clone(),
        solution.temperature.clone(), solution.specific_enthalpy_j_kg.clone(),
        solution.liquid_mass_fraction.clone(), constraints)
}

#[test]
#[cfg(not(target_arch = "wasm32"))]
fn g3_shared_constrained_batches_match_serial_fields_rows_and_robust_loss_bitwise() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|_, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .adiabatic_remainder(), &Environment::earth_lab(), vec![], 0.).unwrap();
        let fail = AtomicBool::new(false);
        let mut owner = study(cx, &m, &boundary, &k, &c, &fail, 9);
        let points = [[0., 0.], [1., 1.], [2., 0.5], [0.5, 2.]];
        let expected: Vec<_> = points.iter().map(|p| signature(owner.evaluate_at(cx, p).unwrap())).collect();
        {
            let batch = owner.value_batch(cx, 4).unwrap();
            let actual = std::thread::scope(|scope| {
                let handles: Vec<_> = points.into_iter().map(|point| {
                    let batch = &batch;
                    scope.spawn(move || {
                        let mut value = None;
                        with_cx(|_, worker| value = Some(signature(batch.evaluate_at(worker, &point).unwrap())));
                        value.unwrap()
                    })
                }).collect();
                handles.into_iter().map(|h| h.join().unwrap()).collect::<Vec<_>>()
            });
            assert_eq!(actual, expected);
            assert_eq!(batch.started_evaluations(), 4);
            assert_eq!(batch.evaluations(), 8);
            assert!(batch.evaluate_at(cx, &[0., 0.]).is_err());
        }
        assert_eq!(owner.evaluations(), 8);
        assert_eq!(signature(owner.evaluate_at(cx, &points[1]).unwrap()), expected[1]);
        assert_eq!(owner.evaluations(), 9);
        assert!(owner.value_batch(cx, 1).is_err());
    });
}

#[test]
fn g4_batch_errors_cancellation_and_unwinding_preserve_started_work_not_unused_slots() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|gate, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .adiabatic_remainder(), &Environment::earth_lab(), vec![], 0.).unwrap();
        let fail = AtomicBool::new(true);
        let mut owner = study(cx, &m, &boundary, &k, &c, &fail, 8);
        assert!(owner.value_batch(cx, usize::MAX).is_err());
        assert_eq!(owner.evaluations(), 0);
        {
            let batch = owner.value_batch(cx, 0).unwrap();
            assert!(batch.evaluate_at(cx, &[1., 1.]).is_err());
            assert_eq!(batch.started_evaluations(), 0);
        }
        {
            let batch = owner.value_batch(cx, 4).unwrap();
            assert!(batch.evaluate_at(cx, &[f64::NAN, 1.]).is_err());
            assert_eq!(batch.started_evaluations(), 0);
            let failure = batch.evaluate_at(cx, &[1., 1.]).unwrap_err();
            assert!(matches!(failure, ScanSensitivityError::Source(
                fs_conduction::ConductionError::Config { parameter: "injected batch factory failure", .. })));
            assert!(batch.evaluate_at(cx, &[1., 1.]).is_ok());
            assert_eq!(batch.started_evaluations(), 2);
        }
        assert_eq!(owner.evaluations(), 2);
        let unwind = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let batch = owner.value_batch(cx, 3).unwrap();
            batch.evaluate_at(cx, &[1., 1.]).unwrap();
            panic!("unwind after one completed physical evaluation");
        }));
        assert!(unwind.is_err());
        assert_eq!(owner.evaluations(), 3);
        {
            let batch = owner.value_batch(cx, 3).unwrap();
            gate.request();
            assert!(batch.evaluate_at(cx, &[1., 1.]).is_err());
            assert_eq!(batch.started_evaluations(), 0);
        }
        assert_eq!(owner.evaluations(), 3);
        with_cx(|_, fresh| {
            assert!(owner.evaluate_at(fresh, &[1., 1.]).is_ok());
            assert_eq!(owner.evaluations(), 4);
            let batch = owner.value_batch(fresh, 4).unwrap();
            std::mem::forget(batch);
            assert_eq!(owner.evaluations(), 8);
            assert!(owner.evaluate_at(fresh, &[1., 1.]).is_err());
        });
    });
}

#[test]
#[cfg(not(target_arch = "wasm32"))]
fn g4_racing_workers_cannot_exceed_one_shared_physical_allowance() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|_, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .adiabatic_remainder(), &Environment::earth_lab(), vec![], 0.).unwrap();
        let fail = AtomicBool::new(false);
        let mut owner = study(cx, &m, &boundary, &k, &c, &fail, 4);
        {
            let batch = owner.value_batch(cx, 3).unwrap();
            let successes = std::thread::scope(|scope| {
                let handles: Vec<_> = (0..12).map(|_| {
                    let batch = &batch;
                    scope.spawn(move || {
                        let mut succeeded = false;
                        with_cx(|_, worker| succeeded = batch.evaluate_at(worker, &[1., 1.]).is_ok());
                        usize::from(succeeded)
                    })
                }).collect();
                handles.into_iter().map(|h| h.join().unwrap()).sum::<usize>()
            });
            assert_eq!(successes, 3);
            assert_eq!(batch.started_evaluations(), 3);
        }
        assert_eq!(owner.evaluations(), 3);
        assert!(owner.evaluate_at(cx, &[1., 1.]).is_ok());
        assert_eq!(owner.evaluations(), 4);
    });
}

#[test]
#[cfg(not(target_arch = "wasm32"))]
fn g4_population_batch_retains_successful_siblings_and_retries_only_the_failed_physics() {
    use fs_dfo::differential_evolution::DifferentialEvolutionConfig;
    use scan_power_fit::population::{PopulationScanStudy, PopulationScanParallelError};
    use fs_dfo::differential_evolution::parallel::DifferentialEvolutionParallelError;
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|_, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .adiabatic_remainder(), &Environment::earth_lab(), vec![], 0.).unwrap();
        let config = |limit| DifferentialEvolutionConfig { population_size:8,
            differential_weight:0.8, crossover_probability:0.9, seed:17,
            max_evaluations:limit, max_population_entries:16, target:None };
        let clean_flag = AtomicBool::new(false);
        let failure_flag = AtomicBool::new(true);
        let mut clean = PopulationScanStudy::new(cx,
            study(cx, &m, &boundary, &k, &c, &clean_flag, 17), &[[0.,4.];2], config(16), 0.).unwrap();
        clean.advance(cx, 16).unwrap();
        let mut retry = PopulationScanStudy::new(cx,
            study(cx, &m, &boundary, &k, &c, &failure_flag, 18), &[[0.,4.];2], config(17), 0.).unwrap();
        let pool = fs_exec::TilePool::for_host(4, 17);
        let lease = fs_alloc::OperationMemoryLease::bounded(64 << 20);
        let failure = retry.advance_parallel(cx, &pool, &lease, 4, 4).unwrap_err();
        assert!(matches!(failure, PopulationScanParallelError::Search(
            DifferentialEvolutionParallelError::Evaluation { source: ScanSensitivityError::Source(
                fs_conduction::ConductionError::Config { parameter: "injected batch factory failure", .. }), .. })));
        assert_eq!(retry.physical_evaluations(), retry.search_evaluations());
        assert!(retry.search_evaluations() >= 1 && retry.search_evaluations() <= 4);
        retry.advance(cx, 17).unwrap();
        assert_eq!(retry.search_evaluations(), 17);
        assert_eq!(retry.physical_evaluations(), 17);
        assert_eq!(retry.population(), clean.population());
        assert_eq!(retry.best(), clean.best());
        let replayed = retry.replay_best(cx).unwrap();
        let expected = clean.replay_best(cx).unwrap();
        assert_eq!(replayed.score, expected.score);
        assert_eq!(replayed.constraints, expected.constraints);
        assert_eq!(replayed.physical.physical.predictions, expected.physical.physical.predictions);
        assert_eq!(replayed.physical_evaluations, 18);
        let exhausted = retry.advance_parallel(cx, &pool, &lease, 20, 8).unwrap();
        assert_eq!(exhausted.physical_evaluations, 18);
        retry.extend_evaluation_budgets(25, 27).unwrap();
        retry.advance_parallel(cx, &pool, &lease, 8, 4).unwrap();
        assert_eq!(retry.physical_evaluations(), 26);
        assert_eq!(retry.replay_best(cx).unwrap().physical_evaluations, 27);
    });
}

#[test]
#[cfg(not(target_arch = "wasm32"))]
fn g0_population_dispatch_refusal_and_pre_cancel_return_unused_reservations() {
    use fs_dfo::differential_evolution::DifferentialEvolutionConfig;
    use scan_power_fit::population::PopulationScanStudy;
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|gate, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .adiabatic_remainder(), &Environment::earth_lab(), vec![], 0.).unwrap();
        let fail = AtomicBool::new(false);
        let mut owner = PopulationScanStudy::new(cx,
            study(cx, &m, &boundary, &k, &c, &fail, 17), &[[0.,4.];2],
            DifferentialEvolutionConfig { population_size:8, differential_weight:0.8,
                crossover_probability:0.9, seed:17, max_evaluations:16,
                max_population_entries:16, target:None }, 0.).unwrap();
        let pool = fs_exec::TilePool::for_host(2, 17);
        let lease = fs_alloc::OperationMemoryLease::bounded(64 << 20);
        assert!(owner.advance_parallel(cx, &pool, &lease, 6, 0).is_err());
        assert_eq!(owner.physical_evaluations(), 0);
        assert_eq!(owner.search_evaluations(), 0);
        assert_eq!(owner.advance_parallel(cx, &pool, &lease, 0, 8).unwrap().physical_evaluations, 0);
        gate.request();
        assert!(owner.advance_parallel(cx, &pool, &lease, 8, 8).is_err());
        assert_eq!(owner.physical_evaluations(), 0);
        assert_eq!(owner.search_evaluations(), 0);
        with_cx(|_, fresh| {
            owner.advance_parallel(fresh, &pool, &lease, 8, 8).unwrap();
            assert_eq!(owner.search_evaluations(), 8);
            assert_eq!(owner.physical_evaluations(), 8);
        });
    });
}
