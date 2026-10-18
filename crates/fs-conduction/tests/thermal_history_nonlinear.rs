//! Nonlinear scenario histories through the production endpoint solver.
mod support;

use fs_conduction::bc::history::{ThermalBoundaryHistory, ThermalRegionHistory};
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryMethod,
    ThermalHistoryProblem, ThermalHistoryRun, ThermalHistoryStep};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    ConductivityTable, LinearConfig, ScalarField, ThermalBc, ThermalBoundaryBuilder, TEMPERATURE_DIMS};
use fs_exec::Cx;
use fs_qty::QtyAny;
use fs_scenario::bc::{BcKind, BcValue, BoundaryCondition, Physics};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::TimeSignal;
use support::{with_cancelled_cx, with_cx};

fn mesh() -> ConductionMesh {
    let (complex, positions) = box_grid([2, 2, 2], [0.01, 0.01, 0.01]);
    ConductionMesh::new(complex, positions).unwrap()
}
fn capacity() -> VolumetricHeatCapacity { VolumetricHeatCapacity::declared(2e6).unwrap() }
fn step_config() -> StepConfig {
    StepConfig { linear: LinearConfig { tolerance: 1e-11, max_iterations: 10_000, restart: 32 },
        energy_tolerance_j: 1e-6 }
}
fn nonlinear() -> NonlinearStepConfig {
    NonlinearStepConfig { residual_rtol: 1e-11, residual_atol_j: 1e-11,
        ..NonlinearStepConfig::default() }
}
fn curved() -> ConductivityModel {
    ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(250.0, 2.0), (450.0, 42.0)]).unwrap())
}
fn history<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh, final_k: f64) -> ThermalBoundaryHistory<'m> {
    let builder = ThermalBoundaryBuilder::new(mesh)
        .region("hot", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::dirichlet(300.0).unwrap()).unwrap()
        .region("air", |f| on_box_face(f.centroid[0], 0.01), ThermalBc::robin(20.0, 300.0).unwrap()).unwrap()
        .adiabatic_remainder();
    let row = BoundaryCondition { region: "hot".into(), physics: Physics::Thermal,
        kind: BcKind::Dirichlet, compatibility: None, frame: 0,
        value: Some(BcValue::Signal(TimeSignal::Ramp {
            t_start: 0.0, t_end: 2.0, from: QtyAny::new(300.0, TEMPERATURE_DIMS),
            to: QtyAny::new(final_k, TEMPERATURE_DIMS) })) };
    ThermalBoundaryHistory::new(cx, builder, &Environment::earth_lab(),
        vec![ThermalRegionHistory { row, reference_temperature: None }], 0.0).unwrap()
}
fn problem<'h, 'm>(boundary: &'h ThermalBoundaryHistory<'m>, material: &'h ConductivityModel)
    -> ThermalHistoryProblem<'h, 'm>
{
    ThermalHistoryProblem { boundary, material, element_materials: None, interfaces: None,
        source: HeatSourceHistory::Fixed(ScalarField::Uniform(0.0)) }
}

#[test]
fn g1_nonlinear_histories_match_direct_endpoints_and_do_not_freeze_conductivity() {
    let mesh = mesh();
    let material = curved();
    let frozen = ConductivityModel::isotropic_declared(12.0).unwrap(); // k(300 K).
    with_cx(|cx| {
        let history = history(cx, &mesh, 400.0);
        let mut old = vec![300.0; mesh.vertex_count()];
        assert!(ThermalHistoryRun::new(cx, problem(&history, &material), capacity(), &old, 4).is_err());
        let method = ThermalHistoryMethod::Nonlinear(nonlinear());
        let mut run = ThermalHistoryRun::new_with_method(cx, problem(&history, &material),
            capacity(), &old, 4, method).unwrap();
        let mut frozen_run = ThermalHistoryRun::new(cx, problem(&history, &frozen), capacity(), &old, 4).unwrap();
        assert_eq!(run.method(), method);
        let engine = BackwardEuler::uniform(cx, &mesh, capacity()).unwrap();
        let mut previous = 0.0;
        let mut updates = 0;
        for time in [0.25, 0.5, 1.0, 2.0] {
            let endpoint = history.at(cx, time).unwrap();
            let source = ScalarField::Uniform(0.0);
            let expected = engine.advance_nonlinear_prescribed(cx, ConductionProblem {
                mesh: &mesh, boundary: &endpoint, material: &material, source: &source,
                element_materials: None }, None, &old, time - previous, step_config(), nonlinear()).unwrap();
            let ThermalHistoryStep::Nonlinear(actual) = run.advance_to_report(cx, time, step_config()).unwrap()
                else { panic!("selected nonlinear solve must retain its evidence") };
            assert_eq!(actual.step.temperature, expected.step.temperature);
            assert_eq!(actual.step.dirichlet_in_w, expected.step.dirichlet_in_w);
            assert_eq!(actual.step.energy_residual_j, expected.step.energy_residual_j);
            assert_eq!(actual.nonlinear_iterations, expected.nonlinear_iterations);
            assert_eq!(actual.residual_j, expected.residual_j);
            assert!(actual.residual_j <= actual.threshold_j);
            assert!(actual.step.energy_residual_j.abs() <= step_config().energy_tolerance_j);
            updates += actual.nonlinear_iterations;
            frozen_run.advance_to(cx, time, step_config()).unwrap();
            old = expected.step.temperature;
            previous = time;
        }
        assert!(updates > 0);
        let difference = run.temperature().iter().zip(frozen_run.temperature())
            .map(|(a, b)| (a - b).abs()).fold(0.0_f64, f64::max);
        assert!(difference > 0.01, "freezing k at the initial temperature must change this trajectory: {difference}");
        assert_eq!(run.accepted_steps(), 4);
        assert_eq!(run.time_s(), 2.0);
    });
}

#[test]
fn g1_constant_material_nonlinear_limit_agrees_with_linear_history() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(12.0).unwrap();
    with_cx(|cx| {
        let history = history(cx, &mesh, 400.0);
        let initial = vec![300.0; mesh.vertex_count()];
        let mut linear = ThermalHistoryRun::new(cx, problem(&history, &material), capacity(), &initial, 3).unwrap();
        let mut nonlinear_run = ThermalHistoryRun::new_with_method(cx, problem(&history, &material),
            capacity(), &initial, 3, ThermalHistoryMethod::Nonlinear(nonlinear())).unwrap();
        for time in [0.25, 1.0, 2.0] {
            let expected = linear.advance_to_report(cx, time, step_config()).unwrap();
            assert!(matches!(&expected, ThermalHistoryStep::Linear(_)));
            let actual = nonlinear_run.advance_to_report(cx, time, step_config()).unwrap();
            for (a, b) in actual.solution().temperature.iter().zip(&expected.solution().temperature) {
                assert!((a - b).abs() < 1e-7);
            }
        }
    });
}

#[test]
fn g0_nonlinear_cancellation_budget_and_material_refusals_preserve_resume() {
    let mesh = mesh();
    let material = curved();
    let history = with_cx(|cx| history(cx, &mesh, 2000.0));
    let initial = vec![300.0; mesh.vertex_count()];
    let mut run = with_cx(|cx| ThermalHistoryRun::new_with_method(cx, problem(&history, &material),
        capacity(), &initial, 1, ThermalHistoryMethod::Nonlinear(nonlinear())).unwrap());
    with_cancelled_cx(|cx| assert!(matches!(run.advance_to_report(cx, 0.01, step_config()),
        Err(ConductionError::Cancelled { .. }))));
    assert_eq!(run.temperature(), initial.as_slice());
    assert_eq!(run.time_s(), 0.0);
    with_cx(|cx| {
        run.advance_to(cx, 0.01, step_config()).unwrap();
        let accepted = run.temperature().to_vec();
        assert!(run.advance_to(cx, 0.02, step_config()).is_err());
        run.extend_step_budget(2).unwrap();
        assert!(run.advance_to(cx, 2.0, step_config()).is_err());
        assert_eq!(run.temperature(), accepted.as_slice());
        assert_eq!(run.time_s(), 0.01);
        assert_eq!(run.accepted_steps(), 1);
        run.advance_to(cx, 0.02, step_config()).unwrap();
        let mut direct = ThermalHistoryRun::new_with_method(cx, problem(&history, &material),
            capacity(), &initial, 2, ThermalHistoryMethod::Nonlinear(nonlinear())).unwrap();
        direct.advance_to(cx, 0.01, step_config()).unwrap();
        direct.advance_to(cx, 0.02, step_config()).unwrap();
        assert_eq!(run.temperature(), direct.temperature());
    });
}

#[test]
fn g0_invalid_nonlinear_policy_and_inner_budget_do_not_publish() {
    let mesh = mesh();
    let material = curved();
    with_cx(|cx| {
        let history = history(cx, &mesh, 400.0);
        let initial = vec![300.0; mesh.vertex_count()];
        let policy = NonlinearStepConfig { max_iterations: 0, ..nonlinear() };
        assert!(ThermalHistoryRun::new_with_method(cx, problem(&history, &material), capacity(),
            &initial, 1, ThermalHistoryMethod::Nonlinear(policy)).is_err());
        let mut run = ThermalHistoryRun::new_with_method(cx, problem(&history, &material), capacity(),
            &initial, 1, ThermalHistoryMethod::Nonlinear(nonlinear())).unwrap();
        let mut invalid = step_config();
        invalid.linear.restart = 0;
        assert!(run.advance_to(cx, 0.25, invalid).is_err());
        assert_eq!(run.temperature(), initial.as_slice());
        assert_eq!(run.time_s(), 0.0);
        assert_eq!(run.accepted_steps(), 0);
        run.advance_to(cx, 0.25, step_config()).unwrap();
    });
}
