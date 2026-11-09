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

fn typed_history(
    mut history: ThermalRegionHistory,
    outside: fs_scenario::payload::OutsideDomainPolicy,
    frame: fs_scenario::FrameId,
) -> ThermalRegionHistory {
    use fs_scenario::payload::{OrientationParity, Payload, PayloadId, PayloadMeta,
        QuantityContract, ReferenceSemantics, SampleSource, ScalarPayload, TableInterpolation};
    let Some(BcValue::Signal(TimeSignal::Table { times, values, dims, interp })) =
        history.row.value.take() else { panic!("test requires an explicit table"); };
    let interpolation = match interp {
        Interp::Hold => TableInterpolation::StepLeft,
        Interp::Linear => TableInterpolation::Linear,
    };
    history.row.frame = frame.0;
    history.row.value = Some(BcValue::Typed(Payload::Scalar(ScalarPayload::new(
        PayloadMeta::new(QuantityContract::Dimensions(dims),
            PayloadId::new("basis/scalar").unwrap(), frame, OrientationParity::Even,
            ReferenceSemantics::Continuous).unwrap(),
        SampleSource::table(
            times.into_iter().map(|t| fs_qty::QtyAny::new(t, Dims([0, 0, 1, 0, 0, 0]))).collect(),
            values.into_iter().map(|v| fs_qty::QtyAny::new(v, dims)).collect(),
            interpolation, outside,
        ).unwrap(),
    ).unwrap())));
    history
}

// The fixture supplies the local direction explicitly; frame rotation and the
// bottom face's -Y normal determine the physical inward/outward exchange.
fn vector_flux_history(mut history: ThermalRegionHistory, direction: [f64; 3]) -> ThermalRegionHistory {
    use fs_scenario::payload::{OrientationParity, Payload, PayloadId, PayloadMeta,
        QuantityContract, ReferenceSemantics, SampleSource, VectorPayload};
    let Some(BcValue::Typed(Payload::Scalar(scalar))) = history.row.value.take() else {
        panic!("fixture requires a typed scalar table");
    };
    let SampleSource::Table { times, values, interpolation, outside } = scalar.source() else {
        panic!("fixture requires a table");
    };
    let vectors = values.iter().map(|q| direction.map(|d| -q.value * d)
        .map(|v| fs_qty::QtyAny::new(v, HEAT_FLUX_DIMS)).to_vec()).collect();
    history.row.value = Some(BcValue::Typed(Payload::Vector(VectorPayload::new(
        PayloadMeta::new(QuantityContract::Dimensions(HEAT_FLUX_DIMS),
            PayloadId::new("basis/frame-cartesian").unwrap(), scalar.meta().frame(),
            OrientationParity::Even, ReferenceSemantics::Continuous).unwrap(),
        SampleSource::table(times.clone(), vectors, *interpolation, *outside).unwrap(),
    ).unwrap())));
    history
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
      for use_typed in [false, true] {
        let mut air = row("right", BcKind::Robin, 0.25, 5.0, 10.0, HTC_DIMS);
        air.reference_temperature = Some(held(vec![0.0, 0.375], vec![300.0, 290.0], TEMPERATURE_DIMS));
        let rows = vec![row("left", BcKind::Dirichlet, 0.125, 300.0, 330.0, TEMPERATURE_DIMS),
            air, row("bottom", BcKind::Neumann, 0.5, 0.0, -20.0, HEAT_FLUX_DIMS)];
        let rows = if use_typed {
            rows.into_iter().map(|row| typed_history(row,
                fs_scenario::payload::OutsideDomainPolicy::Clamp, fs_scenario::WORLD)).collect()
        } else { rows };
        let history = ThermalBoundaryHistory::new(cx, builder(300.0, 5.0, 300.0, 0.0),
            &Environment::earth_lab(), rows, 0.0).unwrap();
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

#[test]
fn g1_typed_surface_pulse_delivers_actual_face_heat_and_visits_ancestor_frame_knots() {
    use fs_scenario::{Frame, FrameId, FrameMotion, FrameTree, WORLD};
    use fs_scenario::payload::OutsideDomainPolicy;
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let mut frames = FrameTree::new();
    frames.add(Frame { id: FrameId(1), name: "source-axes".into(), parent: WORLD,
        motion: FrameMotion::Translating { direction: [1.0, 0.0, 0.0],
            displacement: held(vec![0.0, 0.0625], vec![0.0, 3.0],
                Dims([1, 0, 0, 0, 0, 0])) } });
    frames.add(Frame { id: FrameId(2), name: "turning-flux".into(), parent: FrameId(1),
        motion: FrameMotion::Tilt { axis: [0.0, 0.0, 1.0], center: Default::default(),
            angle: TimeSignal::Ramp { t_start: 0.0, t_end: 0.25,
                from: fs_qty::QtyAny::new(0.0, Dims([0; 6])),
                to: fs_qty::QtyAny::new(std::f64::consts::FRAC_PI_2, Dims([0; 6])) } } });
    for use_vector in [false, true] {
    let mut flux = row("bottom", BcKind::Neumann, 0.125, 0.0, -120.0, HEAT_FLUX_DIMS);
    flux.row.value = Some(BcValue::Signal(held(
        vec![0.0, 0.125, 0.25], vec![0.0, -120.0, 0.0], HEAT_FLUX_DIMS)));
    let flux = typed_history(flux, OutsideDomainPolicy::Clamp,
        if use_vector { FrameId(2) } else { FrameId(1) });
    // At the pulse's right knot, the held local +X vector rotates to +Y at
    // physical time 0.25. Sampling the new zero or freezing the previous frame
    // would both miss the original BE endpoint's 15 J inward contribution.
    let flux = if use_vector { vector_flux_history(flux, [1.0, 0.0, 0.0]) } else { flux };
    let history = with_cx(|cx| ThermalBoundaryHistory::new_in_frames(
        cx, ThermalBoundaryBuilder::new(&mesh)
            .region("bottom", |f| on_box_face(f.centroid[1], 0.0), ThermalBc::adiabatic()).unwrap()
            .adiabatic_remainder(),
        &Environment::earth_lab(), &frames,
        vec![flux], 0.0,
    ).unwrap());
    let mut driver = with_cx(|cx| run(cx, &history, &material,
        HeatSourceHistory::Fixed(ScalarField::Uniform(0.0)), 4));
    with_cancelled_cx(|cx| assert!(matches!(
        driver.advance_event_step(cx, 1.0, 1.0, config()), Err(ConductionError::Cancelled { .. }))));
    assert_eq!(driver.time_s(), 0.0);
    with_cx(|cx| {
        let mut heat_in = 0.0;
        let mut storage = 0.0;
        for (end, outward_w) in [(0.0625, 0.0), (0.125, 0.0), (0.25, -120.0), (1.0, 0.0)] {
            let step = driver.advance_event_step(cx, 1.0, 1.0, config()).unwrap();
            assert_eq!(step.end_time_s, end);
            let solved = step.report.solution();
            // This original exterior patch has area exactly 1 m^2.
            assert!((solved.neumann_out_w - outward_w).abs() < 1e-9);
            heat_in -= solved.neumann_out_w * (step.end_time_s - step.start_time_s);
            storage += solved.stored_energy_change_j;
        }
        assert!((heat_in - 15.0).abs() < 1e-9);
        assert!((storage - 15.0).abs() < 1e-7);
        assert_eq!(driver.accepted_steps(), 4);
        let old = driver.temperature().to_vec();
        assert!(driver.advance_event_step(cx, 2.0, 1.0, config()).is_err());
        assert_eq!(driver.time_s(), 1.0);
        assert_eq!(driver.temperature(), old.as_slice());
    });
    }
}

#[test]
fn g0_typed_tail_policies_refuse_without_losing_accepted_thermal_state() {
    use fs_scenario::payload::OutsideDomainPolicy;
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    with_cx(|cx| {
      for use_vector in [false, true] {
        for outside in [OutsideDomainPolicy::Refuse, OutsideDomainPolicy::Periodic] {
            let mut flux = row("bottom", BcKind::Neumann, 0.125, 0.0, -120.0, HEAT_FLUX_DIMS);
            flux.row.value = Some(BcValue::Signal(held(
                vec![0.0, 0.125, 0.25], vec![0.0, -120.0, 0.0], HEAT_FLUX_DIMS)));
            let flux = typed_history(flux, outside, fs_scenario::WORLD);
            let flux = if use_vector { vector_flux_history(flux, [0.0, 1.0, 0.0]) } else { flux };
            let history = ThermalBoundaryHistory::new(
                cx, ThermalBoundaryBuilder::new(&mesh)
                    .region("bottom", |f| on_box_face(f.centroid[1], 0.0), ThermalBc::adiabatic()).unwrap()
                    .adiabatic_remainder(),
                &Environment::earth_lab(),
                vec![flux], 0.0,
            ).unwrap();
            let mut driver = run(cx, &history, &material,
                HeatSourceHistory::Fixed(ScalarField::Uniform(0.0)), 4);
            if outside == OutsideDomainPolicy::Refuse {
                assert_eq!(driver.advance_event_step(cx, 1.0, 1.0, config()).unwrap().end_time_s, 0.125);
                assert_eq!(driver.advance_event_step(cx, 1.0, 1.0, config()).unwrap().end_time_s, 0.25);
            } else {
                // Endpoint sampling preserves periodic behavior, while finite
                // event stepping requires the author to unroll repeated knots.
                assert_eq!(history.at(cx, 0.375).unwrap().conditions()[0],
                    if use_vector { ThermalBc::neumann_vector([0.0, 120.0, 0.0]).unwrap() }
                    else { ThermalBc::neumann(-120.0).unwrap() });
            }
            let time = driver.time_s();
            let old = driver.temperature().to_vec();
            let steps = driver.accepted_steps();
            assert!(driver.advance_event_step(cx, 1.0, 1.0, config()).is_err());
            assert_eq!(driver.time_s(), time);
            assert_eq!(driver.temperature(), old.as_slice());
            assert_eq!(driver.accepted_steps(), steps);
        }
      }
    });
}
