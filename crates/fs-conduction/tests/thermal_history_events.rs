//! Event-resolved physical work, manually sampled boundary equivalence and resume.
mod support;
use fs_conduction::bc::history::{ThermalBoundaryHistory, ThermalRegionHistory};
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem, ThermalHistoryRun};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    LinearConfig, ScalarField, ThermalBc, ThermalBoundaryBuilder, HEAT_FLUX_DIMS,
    HTC_DIMS, TEMPERATURE_DIMS, VOLUMETRIC_SOURCE_DIMS};
use fs_exec::Cx;
use fs_qty::Dims;
use fs_scenario::bc::{BcKind, BcValue, BoundaryCondition, Physics};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{Interp, TimeSignal};
use support::{with_cancelled_cx, with_cx};

fn mesh() -> ConductionMesh {
    let (complex, positions) = box_grid([2; 3], [1.0; 3]);
    ConductionMesh::new(complex, positions).unwrap()
}
fn config() -> StepConfig {
    StepConfig { linear: LinearConfig { tolerance: 1e-12, ..LinearConfig::default() },
        energy_tolerance_j: 1e-7 }
}
fn held(times: Vec<f64>, values: Vec<f64>, dims: Dims) -> TimeSignal {
    TimeSignal::Table { times, values, dims, interp: Interp::Hold }
}
fn pulse() -> HeatSourceHistory {
    HeatSourceHistory::Uniform(held(vec![0.0, 0.125, 0.25], vec![0.0, 120.0, 0.0], VOLUMETRIC_SOURCE_DIMS))
}
fn run<'h, 'm>(cx: &Cx<'_>, history: &'h ThermalBoundaryHistory<'m>,
    material: &'h ConductivityModel, source: HeatSourceHistory, budget: usize) -> ThermalHistoryRun<'h, 'm>
{
    ThermalHistoryRun::new(cx, ThermalHistoryProblem { boundary: history, material,
        element_materials: None, interfaces: None, source }, VolumetricHeatCapacity::declared(30.0).unwrap(),
        &vec![300.0; history.mesh().vertex_count()], budget).unwrap()
}
fn row(name: &str, kind: BcKind, at: f64, old: f64, new: f64, dims: Dims) -> ThermalRegionHistory {
    ThermalRegionHistory { row: BoundaryCondition { region: name.into(), physics: Physics::Thermal,
        kind, value: Some(BcValue::Signal(held(vec![0.0, at], vec![old, new], dims))),
        compatibility: None, frame: 0 }, reference_temperature: None }
}

#[test]
fn g1_short_pulse_is_not_lost_or_backdated_at_its_switches() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    with_cx(|cx| {
        let history = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder(), &Environment::earth_lab(), vec![], 0.0).unwrap();
        let sources = [pulse(), HeatSourceHistory::Scaled {
            field: ScalarField::Nodal(vec![120.0; mesh.vertex_count()]),
            scale: held(vec![0.0, 0.125, 0.25], vec![0.0, 1.0, 0.0], Dims([0; 6])) }];
        for source in sources {
            let mut run = run(cx, &history, &material, source, 3);
            let mut work = 0.0;
            for (end, power) in [(0.125, 0.0), (0.25, 120.0), (1.0, 0.0)] {
                let step = run.advance_event_step(cx, 1.0, 1.0, config()).unwrap();
                assert_eq!(step.end_time_s, end);
                let solved = step.report.solution();
                assert!((solved.source_w - power).abs() < 1e-9);
                work += solved.source_w * (step.end_time_s - step.start_time_s);
            }
            assert!((work - 15.0).abs() < 1e-9);
            assert!(run.temperature().iter().all(|&t| (t - 300.5).abs() < 1e-9));
            assert_eq!(run.time_s(), 1.0);
        }
        // Negative control: the old endpoint API deliberately samples q(1)=0.
        let mut unaligned = run(cx, &history, &material, pulse(), 1);
        unaligned.advance_to(cx, 1.0, config()).unwrap();
        assert!(unaligned.temperature().iter().all(|&t| (t - 300.0).abs() < 1e-9));
    });
}

#[test]
fn g1_boundary_switches_use_pre_jump_data_and_the_original_physical_solver() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let builder = |g, h, reference, flux| ThermalBoundaryBuilder::new(&mesh)
        .region("left", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::dirichlet(g).unwrap()).unwrap()
        .region("right", |f| on_box_face(f.centroid[0], 1.0), ThermalBc::robin(h, reference).unwrap()).unwrap()
        .region("bottom", |f| on_box_face(f.centroid[1], 0.0), ThermalBc::neumann(flux).unwrap()).unwrap()
        .adiabatic_remainder();
    with_cx(|cx| {
        let mut air = row("right", BcKind::Robin, 0.25, 5.0, 10.0, HTC_DIMS);
        air.reference_temperature = Some(held(vec![0.0, 0.375], vec![300.0, 290.0], TEMPERATURE_DIMS));
        let history = ThermalBoundaryHistory::new(cx, builder(300.0, 5.0, 300.0, 0.0),
            &Environment::earth_lab(), vec![row("left", BcKind::Dirichlet, 0.125, 300.0, 330.0, TEMPERATURE_DIMS),
                air, row("bottom", BcKind::Neumann, 0.5, 0.0, -20.0, HEAT_FLUX_DIMS)], 0.0).unwrap();
        let mut run = run(cx, &history, &material, HeatSourceHistory::Fixed(ScalarField::Uniform(0.0)), 5);
        let manual = BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(30.0).unwrap()).unwrap();
        let mut old = run.temperature().to_vec();
        let mut start = 0.0;
        for (end, g, h, reference, flux) in [(0.125, 300.0, 5.0, 300.0, 0.0),
            (0.25, 330.0, 5.0, 300.0, 0.0), (0.375, 330.0, 10.0, 300.0, 0.0),
            (0.5, 330.0, 10.0, 290.0, 0.0), (0.75, 330.0, 10.0, 290.0, -20.0)]
        {
            let boundary = builder(g, h, reference, flux).finish().unwrap();
            let source = ScalarField::Uniform(0.0);
            let expected = manual.advance_prescribed(cx, ConductionProblem { mesh: &mesh,
                boundary: &boundary, material: &material, source: &source, element_materials: None },
                None, &old, end-start, config()).unwrap();
            let actual = run.advance_event_step(cx, 0.75, 1.0, config()).unwrap();
            assert_eq!(actual.end_time_s, end);
            assert_eq!(actual.report.solution().temperature, expected.temperature);
            assert_eq!(actual.report.solution().dirichlet_in_w, expected.dirichlet_in_w);
            assert_eq!(actual.report.solution().energy_residual_j, expected.energy_residual_j);
            old = expected.temperature;
            start = end;
        }
    });
}

#[test]
fn g0_cancellation_step_limits_and_budget_exhaustion_do_not_lose_schedule_position() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let history = with_cx(|cx| ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder(), &Environment::earth_lab(), vec![], 0.0).unwrap());
    let mut driver = with_cx(|cx| run(cx, &history, &material, pulse(), 1));
    with_cancelled_cx(|cx| assert!(matches!(driver.advance_event_step(cx, 1.0, 1.0, config()),
        Err(ConductionError::Cancelled { .. }))));
    assert_eq!(driver.time_s(), 0.0);
    with_cx(|cx| {
        driver.advance_event_step(cx, 1.0, 1.0, config()).unwrap();
        let old = driver.temperature().to_vec();
        assert!(driver.advance_event_step(cx, 1.0, 1.0, config()).is_err());
        driver.extend_step_budget(4).unwrap();
        for max_step in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MIN_POSITIVE] {
            assert!(driver.advance_event_step(cx, 1.0, max_step, config()).is_err());
            assert_eq!(driver.temperature(), old.as_slice());
            assert_eq!(driver.time_s(), 0.125);
        }
        // The cap inserts an actual extra physical step before the next event.
        let a = driver.advance_event_step(cx, 1.0, 0.0625, config()).unwrap();
        assert_eq!(a.end_time_s, 0.1875);
        driver.advance_event_step(cx, 1.0, 1.0, config()).unwrap();
        driver.advance_event_step(cx, 1.0, 1.0, config()).unwrap();
        assert!(driver.temperature().iter().all(|&t| (t - 300.5).abs() < 1e-9));
        assert_eq!(driver.accepted_steps(), 4);
    });
}
