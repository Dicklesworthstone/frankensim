//! End-to-end thermal history stepping and transactional resume regressions.
mod support;

use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::bc::history::{ThermalBoundaryHistory, ThermalRegionHistory};
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem, ThermalHistoryRun};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    LinearConfig, ScalarField, HEAT_FLUX_DIMS, HTC_DIMS, TEMPERATURE_DIMS, VOLUMETRIC_SOURCE_DIMS};
use fs_qty::{Dims, QtyAny};
use fs_scenario::bc::{BcKind, BcValue, BoundaryCondition, Physics};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::TimeSignal;
use support::{with_cancelled_cx, with_cx};

fn mesh(n: usize) -> ConductionMesh {
    let (complex, positions) = box_grid([n, n, n], [1.0, 1.0, 1.0]);
    ConductionMesh::new(complex, positions).unwrap()
}
fn capacity() -> VolumetricHeatCapacity { VolumetricHeatCapacity::declared(2e6).unwrap() }
fn config() -> StepConfig {
    StepConfig { linear: LinearConfig { tolerance: 1e-11, max_iterations: 60_000, restart: 60 },
        energy_tolerance_j: 1e-5 }
}
fn ramp(from: f64, to: f64, dims: Dims) -> TimeSignal {
    TimeSignal::Ramp { t_start: 0.0, t_end: 2.0,
        from: QtyAny::new(from, dims), to: QtyAny::new(to, dims) }
}
fn history_row(name: &str, kind: BcKind, from: f64, to: f64, dims: Dims) -> ThermalRegionHistory {
    ThermalRegionHistory { row: BoundaryCondition { region: name.into(), physics: Physics::Thermal,
        kind, value: Some(BcValue::Signal(ramp(from, to, dims))), compatibility: None, frame: 0 },
        reference_temperature: None }
}
fn problem<'h, 'm>(history: &'h ThermalBoundaryHistory<'m>, material: &'h ConductivityModel,
    source: HeatSourceHistory) -> ThermalHistoryProblem<'h, 'm>
{
    ThermalHistoryProblem { boundary: history, material, element_materials: None, interfaces: None, source }
}

#[test]
fn g1_prescribed_temperature_ramp_accounts_for_boundary_node_storage() {
    let mesh = mesh(1); // Every node lies on the prescribed exterior.
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    with_cx(|cx| {
        let history = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
            .remainder("surface", ThermalBc::dirichlet(300.0).unwrap()).unwrap(),
            &Environment::earth_lab(), vec![history_row("surface", BcKind::Dirichlet,
                300.0, 400.0, TEMPERATURE_DIMS)], 0.0).unwrap();
        let mut run = ThermalHistoryRun::new(cx, problem(&history, &material,
            HeatSourceHistory::Fixed(ScalarField::Uniform(0.0))), capacity(),
            &vec![300.0; mesh.vertex_count()], 3).unwrap();
        let mut previous = 0.0;
        for time in [0.5, 1.0, 2.0] {
            let step = run.advance_to(cx, time, config()).unwrap();
            assert!(step.temperature.iter().all(|&t| t == 300.0 + 50.0 * time));
            let expected_j = 2e6 * 50.0 * (time - previous);
            assert!((step.stored_energy_change_j - expected_j).abs() < 1e-5);
            assert!((step.dirichlet_in_w * (time - previous) - expected_j).abs() < 1e-5);
            assert!(step.energy_residual_j.abs() <= config().energy_tolerance_j);
            previous = time;
        }
        assert_eq!(run.accepted_steps(), 3);
        assert_eq!(run.time_s(), 2.0);
    });
}

#[test]
fn g1_mixed_boundary_and_source_histories_match_manually_assembled_endpoints() {
    let mesh = mesh(2);
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let builder = |time: f64| ThermalBoundaryBuilder::new(&mesh)
        .region("hot", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::dirichlet(300.0 + 50.0*time).unwrap()).unwrap()
        .region("air", |f| on_box_face(f.centroid[0], 1.0), ThermalBc::robin(10.0 + 10.0*time, 290.0 + 10.0*time).unwrap()).unwrap()
        .region("flux", |f| on_box_face(f.centroid[1], 0.0), ThermalBc::neumann(-100.0 - 100.0*time).unwrap()).unwrap()
        .adiabatic_remainder();
    with_cx(|cx| {
        let mut air = history_row("air", BcKind::Robin, 10.0, 30.0, HTC_DIMS);
        air.reference_temperature = Some(ramp(290.0, 310.0, TEMPERATURE_DIMS));
        let history = ThermalBoundaryHistory::new(cx, builder(0.0), &Environment::earth_lab(),
            vec![history_row("hot", BcKind::Dirichlet, 300.0, 400.0, TEMPERATURE_DIMS), air,
                history_row("flux", BcKind::Neumann, -100.0, -300.0, HEAT_FLUX_DIMS)], 0.0).unwrap();
        let mut old = vec![300.0; mesh.vertex_count()];
        let mut run = ThermalHistoryRun::new(cx, problem(&history, &material,
            HeatSourceHistory::Uniform(ramp(1e5, 3e5, VOLUMETRIC_SOURCE_DIMS))),
            capacity(), &old, 4).unwrap();
        let manual = BackwardEuler::uniform(cx, &mesh, capacity()).unwrap();
        let mut previous = 0.0;
        for time in [0.25, 0.5, 1.0, 1.5] {
            let boundary = builder(time).finish().unwrap();
            let source = ScalarField::Uniform(1e5 + 1e5*time);
            let expected = manual.advance_prescribed(cx, ConductionProblem {
                mesh: &mesh, boundary: &boundary, material: &material, source: &source,
                element_materials: None }, None, &old, time-previous, config()).unwrap();
            let actual = run.advance_to(cx, time, config()).unwrap();
            assert_eq!(actual.temperature, expected.temperature);
            assert_eq!(actual.dirichlet_in_w, expected.dirichlet_in_w);
            assert_eq!(actual.neumann_out_w, expected.neumann_out_w);
            assert_eq!(actual.robin_out_w, expected.robin_out_w);
            assert_eq!(actual.source_w, expected.source_w);
            assert_eq!(actual.energy_residual_j, expected.energy_residual_j);
            old = expected.temperature;
            previous = time;
        }
        assert_eq!(run.temperature(), old.as_slice());
    });
}

#[test]
fn g1_uniform_and_scaled_nodal_workloads_follow_endpoint_energy_quadrature() {
    let mesh = mesh(2);
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    with_cx(|cx| {
        let history = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder(), &Environment::earth_lab(), Vec::new(), 0.0).unwrap();
        let sources = [HeatSourceHistory::Uniform(ramp(4e5, 8e5, VOLUMETRIC_SOURCE_DIMS)),
            HeatSourceHistory::Scaled { field: ScalarField::Nodal(vec![4e5; mesh.vertex_count()]),
                scale: ramp(1.0, 2.0, Dims([0; 6])) }];
        for source in sources {
            let mut run = ThermalHistoryRun::new(cx, problem(&history, &material, source),
                capacity(), &vec![300.0; mesh.vertex_count()], 4).unwrap();
            let (mut previous, mut expected) = (0.0, 300.0);
            for time in [0.25, 0.5, 1.0, 2.0] {
                let q = 4e5 + 2e5*time;
                expected += q * (time-previous) / 2e6;
                let step = run.advance_to(cx, time, config()).unwrap();
                assert!(step.temperature.iter().all(|&t| (t-expected).abs() < 1e-8));
                assert!((step.source_w-q).abs() < 1e-7*q);
                assert!(step.energy_residual_j.abs() <= config().energy_tolerance_j);
                previous = time;
            }
        }
    });
}

#[test]
fn g0_cancelled_invalid_and_exhausted_attempts_preserve_resumable_state() {
    let mesh = mesh(2);
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let history = with_cx(|cx| ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder(), &Environment::earth_lab(), Vec::new(), 0.0).unwrap());
    let mut run = with_cx(|cx| ThermalHistoryRun::new(cx, problem(&history, &material,
        HeatSourceHistory::Uniform(ramp(4e5, 8e5, VOLUMETRIC_SOURCE_DIMS))), capacity(),
        &vec![300.0; mesh.vertex_count()], 1).unwrap());
    with_cancelled_cx(|cx| assert!(matches!(run.advance_to(cx, 0.25, config()), Err(ConductionError::Cancelled { .. }))));
    assert_eq!(run.accepted_steps(), 0);
    with_cx(|cx| { run.advance_to(cx, 0.25, config()).unwrap(); });
    let accepted = run.temperature().to_vec();
    with_cx(|cx| assert!(run.advance_to(cx, 0.5, config()).is_err()));
    assert_eq!(run.temperature(), accepted.as_slice());
    assert_eq!(run.time_s(), 0.25);
    run.extend_step_budget(2).unwrap();
    with_cx(|cx| {
        for time in [f64::NAN, f64::INFINITY, 0.25, -1.0] {
            assert!(run.advance_to(cx, time, config()).is_err());
            assert_eq!(run.temperature(), accepted.as_slice());
            assert_eq!(run.accepted_steps(), 1);
        }
        run.advance_to(cx, 0.5, config()).unwrap();
        let mut uninterrupted = ThermalHistoryRun::new(cx, problem(&history, &material,
            HeatSourceHistory::Uniform(ramp(4e5, 8e5, VOLUMETRIC_SOURCE_DIMS))), capacity(),
            &vec![300.0; mesh.vertex_count()], 2).unwrap();
        uninterrupted.advance_to(cx, 0.25, config()).unwrap();
        uninterrupted.advance_to(cx, 0.5, config()).unwrap();
        assert_eq!(run.temperature(), uninterrupted.temperature());
        assert_eq!(run.time_s(), uninterrupted.time_s());
    });
}

#[test]
fn g0_bad_future_robin_endpoint_does_not_commit_a_partial_thermal_step() {
    let mesh = mesh(2);
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    with_cx(|cx| {
        let history = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
            .remainder("air", ThermalBc::robin(10.0, 293.15).unwrap()).unwrap(),
            &Environment::earth_lab(), vec![history_row("air", BcKind::Robin, 10.0, 0.0, HTC_DIMS)], 0.0).unwrap();
        let mut run = ThermalHistoryRun::new(cx, problem(&history, &material,
            HeatSourceHistory::Fixed(ScalarField::Uniform(0.0))), capacity(),
            &vec![300.0; mesh.vertex_count()], 3).unwrap();
        run.advance_to(cx, 0.25, config()).unwrap();
        let accepted = run.temperature().to_vec();
        assert!(run.advance_to(cx, 2.0, config()).is_err());
        assert_eq!(run.temperature(), accepted.as_slice());
        assert_eq!(run.time_s(), 0.25);
        assert_eq!(run.accepted_steps(), 1);
        run.advance_to(cx, 0.5, config()).unwrap();
    });
}

#[test]
fn g0_initial_boundary_mismatch_and_wrong_source_units_are_not_silently_repaired() {
    let mesh = mesh(1);
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    with_cx(|cx| {
        let history = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
            .remainder("surface", ThermalBc::dirichlet(300.0).unwrap()).unwrap(),
            &Environment::earth_lab(), Vec::new(), 0.0).unwrap();
        let source = HeatSourceHistory::Fixed(ScalarField::Uniform(0.0));
        assert!(ThermalHistoryRun::new(cx, problem(&history, &material, source), capacity(),
            &vec![310.0; mesh.vertex_count()], 1).is_err());
        let source = HeatSourceHistory::Uniform(ramp(1.0, 2.0, HEAT_FLUX_DIMS));
        assert!(ThermalHistoryRun::new(cx, problem(&history, &material, source), capacity(),
            &vec![300.0; mesh.vertex_count()], 1).is_err());
    });
}
