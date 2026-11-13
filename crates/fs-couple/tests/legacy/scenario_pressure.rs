//! G0/G3: physical pressure histories drive the original native reservoirs.

use fs_couple::scenario_pressure::{
    PressureBoundary3, PressureReference, ScenarioPressureError, lattice3,
    rebind_pressure_face3, rebind_pressure_faces3, sample_pressure3,
};
use fs_exec::{CancelGate, Cx};
use fs_ga::{Quat, Vec3};
use fs_lbm::d3q19::{BoundaryGrid3, BoundarySpec3, CollisionModel3, Face3, FaceBoundary3};
use fs_lbm::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27};
use fs_lbm::lattice3::open_boundaries::OpenBoundaryGrid3;
use fs_lbm::lattice3::{D3q19, EquilibriumLattice3, Grid3, GridBudget3};
use fs_lbm::material::plan_scaling_from_material;
use fs_material::fluid_transport::FluidTransportState;
use fs_qty::{Dims, QtyAny};
use fs_scenario::{
    BcKind, BcValue, BoundaryCondition, Frame, FrameId, FrameMotion, FrameTree, Interp,
    Physics, TimeSignal, WORLD,
};
use std::ops::ControlFlow;

const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);

fn fluid() -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 4200.0, 0.6, 0.001).unwrap()
}

fn boundary(value: BcValue) -> BoundaryCondition {
    BoundaryCondition {
        region: "pressure_reservoir".into(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::PressureOutlet,
        value: Some(value),
        compatibility: None,
        frame: WORLD.0,
    }
}

fn uniform(value: f64) -> BoundaryCondition {
    boundary(BcValue::Uniform(QtyAny::new(value, PRESSURE)))
}

fn gauge(boundary: &BoundaryCondition) -> PressureBoundary3<'_> {
    PressureBoundary3 { boundary, reference: PressureReference::Gauge }
}

fn absolute(boundary: &BoundaryCondition) -> PressureBoundary3<'_> {
    PressureBoundary3 {
        boundary,
        reference: PressureReference::Absolute { reference_pa: 101_325.0 },
    }
}

fn faces() -> BoundarySpec3 {
    BoundarySpec3::new([
        FaceBoundary3::Periodic,
        FaceBoundary3::Periodic,
        FaceBoundary3::Periodic,
        FaceBoundary3::Periodic,
        FaceBoundary3::Pressure { density: 1.0 },
        FaceBoundary3::Pressure { density: 1.0 },
    ])
}

fn grid(tau: f64) -> BoundaryGrid3 {
    BoundaryGrid3::new(8, 8, 8, tau, [1e-5, -2e-5, 3e-5], faces())
}

fn density(pressure_pa: f64) -> f64 {
    // Independent SI mapping: rho_ref=1000 kg/m^3 and dx/dt=20 m/s.
    1.0 + pressure_pa / ((1.0 / 3.0) * 1000.0 * 20.0 * 20.0)
}

fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 5e-14, "{actual} != {expected}");
}

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            fs_exec::StreamKey { seed: 0xCE11, kernel_id: 1, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

#[test]
fn gauge_ramp_and_absolute_table_reach_original_d3q19_density_boundaries() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut flow = grid(scaling.flow().tau);
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    let inlet = boundary(BcValue::Signal(TimeSignal::Ramp {
        t_start: 0.0, t_end: 1.0,
        from: QtyAny::new(200.0, PRESSURE),
        to: QtyAny::new(-200.0, PRESSURE),
    }));
    let outlet = boundary(BcValue::Signal(TimeSignal::Table {
        times: vec![0.0, 0.5, 1.0],
        values: vec![101_205.0, 101_405.0, 101_285.0],
        dims: PRESSURE,
        interp: Interp::Linear,
    }));
    for (time, expected) in [
        (0.0, [200.0, -120.0]),
        (0.25, [100.0, -20.0]),
        (0.5, [0.0, 80.0]),
        (1.0, [-200.0, -40.0]),
    ] {
        let before = (flow.populations(3, 3, 0), flow.populations(3, 3, 7));
        let force = flow.force_density(3, 3, 0);
        let samples = rebind_pressure_faces3(
            &mut flow, [None, None, None, None, Some(gauge(&inlet)), Some(absolute(&outlet))],
            &frames, time, &scaling, &gate,
        ).unwrap();
        assert!(samples[..4].iter().all(Option::is_none));
        for (sample, pa) in samples[4..].iter().zip(expected) {
            let sample = sample.unwrap();
            assert_eq!(sample.time_s, time);
            assert_eq!(sample.perturbation_pa, pa);
            near(sample.density_lattice, density(pa));
        }
        assert_eq!((flow.populations(3, 3, 0), flow.populations(3, 3, 7)), before);
        assert_eq!(flow.force_density(3, 3, 0), force);
        flow.step();
        for (z, pa) in [0, 7].into_iter().zip(expected) {
            for y in 0..8 {
                for x in 0..8 {
                    near(flow.density(x, y, z), density(pa));
                }
            }
        }
    }
}

#[test]
fn scalar_pressure_keeps_its_reference_under_rigid_frame_changes_and_hold_sampling() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut frames = FrameTree::new();
    frames.add(Frame {
        id: FrameId(1), name: "reservoir_axes".into(), parent: WORLD,
        motion: FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), 0.7),
            translation: Vec3::new(1e6, -2e6, 3e6),
        },
    });
    let original = uniform(-80.0);
    let mut moved = uniform(101_245.0);
    moved.frame = 1;
    let mut checkpoint = || ControlFlow::Continue(());
    let first = sample_pressure3(gauge(&original), &frames, 0.3, &scaling, &mut checkpoint).unwrap();
    let second = sample_pressure3(absolute(&moved), &frames, 0.3, &scaling, &mut checkpoint).unwrap();
    assert_eq!(first.perturbation_pa, second.perturbation_pa);
    assert_eq!(first.density_lattice, second.density_lattice);
    assert_eq!(second.supplied_pa, 101_245.0);

    let held = boundary(BcValue::Signal(TimeSignal::Table {
        times: vec![0.0, 0.5, 1.0],
        values: vec![40.0, -20.0, 80.0], dims: PRESSURE, interp: Interp::Hold,
    }));
    for (time, pa) in [(-1.0, 40.0), (0.49, 40.0), (0.5, -20.0), (2.0, 80.0)] {
        let sample = sample_pressure3(gauge(&held), &frames, time, &scaling, &mut checkpoint).unwrap();
        assert_eq!(sample.perturbation_pa, pa);
    }
}

#[test]
fn a_later_invalid_history_or_native_face_keeps_the_entire_pressure_batch() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut flow = grid(scaling.flow().tau);
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    let inlet = uniform(100.0);
    let invalid = boundary(BcValue::Signal(TimeSignal::Table {
        times: vec![0.0, 1.0], values: vec![0.0, f64::NAN],
        dims: PRESSURE, interp: Interp::Linear,
    }));
    let before = flow.boundaries();
    let populations = (flow.populations(3, 3, 0), flow.populations(3, 3, 7));
    assert!(matches!(
        rebind_pressure_faces3(
            &mut flow, [None, None, None, None, Some(gauge(&inlet)), Some(gauge(&invalid))],
            &frames, 0.0, &scaling, &gate,
        ),
        Err(ScenarioPressureError::Scenario(_))
    ));
    // Every source is valid, but a periodic face is not a pressure outlet.
    assert!(matches!(
        rebind_pressure_faces3(
            &mut flow, [Some(gauge(&inlet)), None, None, None, Some(gauge(&inlet)), None],
            &frames, 0.0, &scaling, &gate,
        ),
        Err(ScenarioPressureError::Boundary(_))
    ));
    assert_eq!(flow.boundaries(), before);
    assert_eq!((flow.populations(3, 3, 0), flow.populations(3, 3, 7)), populations);
}

#[test]
fn selected_pressure_profile_refuses_even_an_unchanged_dormant_target() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut flow = grid(scaling.flow().tau);
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    flow.set_pressure_profile(Face3::ZMax, &gate, |_| Ok(1.002)).unwrap();
    let profile = flow.pressure_profile(Face3::ZMax).unwrap().to_vec();
    let before = flow.boundaries();
    let inlet = uniform(40.0);
    let dormant = uniform(0.0);
    assert!(matches!(
        rebind_pressure_faces3(
            &mut flow, [None, None, None, None, Some(gauge(&inlet)), Some(gauge(&dormant))],
            &frames, 0.0, &scaling, &gate,
        ),
        Err(ScenarioPressureError::Boundary(_))
    ));
    assert_eq!(flow.boundaries(), before);
    assert_eq!(flow.pressure_profile(Face3::ZMax).unwrap(), profile);
    rebind_pressure_face3(
        &mut flow, Face3::ZMin, gauge(&inlet), &frames, 0.0, &scaling, &gate,
    ).unwrap();
    assert_eq!(flow.pressure_profile(Face3::ZMax).unwrap(), profile);
    flow.step();
    near(flow.density(3, 3, 0), density(40.0));
    near(flow.density(3, 3, 7), 1.002);
}

#[test]
fn invalid_reference_units_frames_and_geometry_profiles_publish_no_target() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut flow = grid(scaling.flow().tau);
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    let source = uniform(101_325.0);
    let before = flow.boundaries();
    let populations = flow.populations(3, 3, 7);
    for reference_pa in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(rebind_pressure_face3(
            &mut flow, Face3::ZMax,
            PressureBoundary3 { boundary: &source, reference: PressureReference::Absolute { reference_pa } },
            &frames, 0.0, &scaling, &gate,
        ).is_err());
    }
    let mut bad_units = uniform(1.0);
    bad_units.value = Some(BcValue::Uniform(QtyAny::new(1.0, Dims([0, 0, 0, 1, 0, 0]))));
    let mut missing_frame = uniform(1.0);
    missing_frame.frame = 99;
    let mut wrong_role = uniform(1.0);
    wrong_role.kind = BcKind::Traction;
    for invalid in [&bad_units, &missing_frame, &wrong_role, &uniform(-1e6)] {
        assert!(rebind_pressure_face3(
            &mut flow, Face3::ZMax, gauge(invalid), &frames, 0.0, &scaling, &gate,
        ).is_err());
    }
    assert!(rebind_pressure_face3(
        &mut flow, Face3::ZMax, absolute(&uniform(0.0)), &frames, 0.0, &scaling, &gate,
    ).is_err());

    // The original decoder constructs an actual spatial Chebyshev profile;
    // the time-only bridge must not silently evaluate it at physical time.
    let parsed = fs_scenario::ir::parse_ir(r#"(scenario :version 2 "pressure-profile" 7
      (environment (qty 0 1 0 -2 0 0 0) (qty 0 1 0 -2 0 0 0)
        (qty 0 1 0 -2 0 0 0) (qty 293 0 0 0 1 0 0) (qty 101325 -1 1 -2 0 0 0))
      (frames)
      (bcs (bc "outlet" incompressible-flow pressure-outlet 0
        (profile (dims -1 1 -2 0 0 0) 0 1 (coeffs 20 5)) none))
      (cases) (combos) (ensembles) (contacts))"#).unwrap();
    assert!(matches!(
        rebind_pressure_face3(
            &mut flow, Face3::ZMax, gauge(&parsed.scenario().base_bcs[0]),
            &frames, 0.5, &scaling, &gate,
        ),
        Err(ScenarioPressureError::Invalid(_))
    ));
    assert_eq!(flow.boundaries(), before);
    assert_eq!(flow.populations(3, 3, 7), populations);
}

#[test]
fn cancellation_and_empty_batches_preserve_controls_and_signal_sampling_can_stop() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut flow = grid(scaling.flow().tau);
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    let before = flow.boundaries();
    let populations = flow.populations(3, 3, 7);
    assert_eq!(
        rebind_pressure_faces3(&mut flow, [None; 6], &frames, 0.0, &scaling, &gate).unwrap(),
        [None; 6]
    );
    assert!(rebind_pressure_faces3(
        &mut flow, [None; 6], &frames, f64::NAN, &scaling, &gate,
    ).is_err());
    gate.request();
    assert!(matches!(
        rebind_pressure_face3(
            &mut flow, Face3::ZMax, gauge(&uniform(20.0)), &frames, 0.0, &scaling, &gate,
        ),
        Err(ScenarioPressureError::Cancelled)
    ));
    let long = boundary(BcValue::Signal(TimeSignal::Table {
        times: (0..64).map(f64::from).collect(),
        values: vec![40.0; 64], dims: PRESSURE, interp: Interp::Linear,
    }));
    let mut polls = 0;
    assert!(matches!(
        sample_pressure3(gauge(&long), &frames, 0.0, &scaling, &mut || {
            polls += 1;
            if polls == 8 { ControlFlow::Break(()) } else { ControlFlow::Continue(()) }
        }),
        Err(ScenarioPressureError::Cancelled)
    ));
    assert_eq!(polls, 8);
    assert_eq!(flow.boundaries(), before);
    assert_eq!(flow.populations(3, 3, 7), populations);
}

fn check_generic<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let native = Grid3::<L, Q>::try_new_equilibrium(
            [4, 4, 5], collision, [2e-5, -1e-5, 3e-5], [true, true, false],
            GridBudget3 { maximum_cells: 80, maximum_mach: 0.2 },
            |_| false, |_| (1.0, [0.0; 3]),
        ).unwrap();
        let mut flow = OpenBoundaryGrid3::try_new(native, cx, faces()).unwrap();
        let inlet = boundary(BcValue::Signal(TimeSignal::Table {
            times: vec![0.0, 1.0], values: vec![80.0, 160.0],
            dims: PRESSURE, interp: Interp::Linear,
        }));
        let outlet = uniform(101_245.0);
        let before = flow.populations().to_vec();
        let force = flow.force_density();
        let samples = lattice3::rebind_pressure_faces3(
            &mut flow, [None, None, None, None, Some(gauge(&inlet)), Some(absolute(&outlet))],
            &frames, 0.5, &scaling, &gate,
        ).unwrap();
        assert_eq!(flow.populations(), before);
        assert_eq!(flow.completed_steps(), 0);
        assert_eq!(flow.force_density(), force);
        near(samples[4].unwrap().density_lattice, density(120.0));
        near(samples[5].unwrap().density_lattice, density(-80.0));
        flow.step_serial(cx, None, [0.0; 3]).unwrap();
        assert_eq!(flow.completed_steps(), 1);
        for y in 0..4 {
            for x in 0..4 {
                near(flow.macroscopic(y * 4 + x).unwrap().0, density(120.0));
                near(flow.macroscopic((4 * 4 + y) * 4 + x).unwrap().0, density(-80.0));
            }
        }
        let before = flow.boundaries();
        let populations = flow.populations().to_vec();
        let changed = uniform(60.0);
        let impossible = uniform(-1e6);
        assert!(matches!(
            lattice3::rebind_pressure_faces3(
                &mut flow, [None, None, None, None, Some(gauge(&changed)), Some(gauge(&impossible))],
                &frames, 1.0, &scaling, &gate,
            ),
            Err(lattice3::Error::Scenario(ScenarioPressureError::Scaling(_)))
        ));
        assert!(matches!(
            lattice3::rebind_pressure_face3(
                &mut flow, Face3::XMin, gauge(&changed), &frames, 1.0, &scaling, &gate,
            ),
            Err(lattice3::Error::Boundary(_))
        ));
        assert_eq!(flow.boundaries(), before);
        assert_eq!(flow.populations(), populations);
        assert_eq!(flow.completed_steps(), 1);
        assert_eq!(flow.force_density(), force);

        // This is the outlet's already-held uniform target. A native rebind
        // alone accepts it, even though the spatial profile still wins.
        let dormant = uniform(-80.0);
        flow.set_pressure_profile(Face3::ZMax, &gate, |point| {
            Ok(1.002 + 0.0001 * point[0] as f64)
        }).unwrap();
        let profile = flow.pressure_profile(Face3::ZMax).unwrap().to_vec();
        assert!(matches!(
            lattice3::rebind_pressure_faces3(
                &mut flow,
                [None, None, None, None, Some(gauge(&changed)), Some(gauge(&dormant))],
                &frames, 1.0, &scaling, &gate,
            ),
            Err(lattice3::Error::Scenario(ScenarioPressureError::Boundary(_)))
        ));
        assert_eq!(flow.boundaries(), before);
        assert_eq!(flow.populations(), populations);
        assert_eq!(flow.completed_steps(), 1);
        assert_eq!(flow.force_density(), force);
        assert_eq!(flow.pressure_profile(Face3::ZMax).unwrap(), profile);

        // An unselected profile remains a valid physical boundary while the
        // other face advances its scenario history.
        let samples = lattice3::rebind_pressure_faces3(
            &mut flow, [None, None, None, None, Some(gauge(&changed)), None],
            &frames, 1.0, &scaling, &gate,
        ).unwrap();
        assert!(samples[5].is_none());
        assert_eq!(flow.pressure_profile(Face3::ZMax).unwrap(), profile);
        flow.step_serial(cx, None, [0.0; 3]).unwrap();
        for y in 0..4 {
            for x in 0..4 {
                near(flow.macroscopic(y * 4 + x).unwrap().0, density(60.0));
                near(
                    flow.macroscopic((4 * 4 + y) * 4 + x).unwrap().0,
                    1.002 + 0.0001 * x as f64,
                );
            }
        }

        // Explicitly relinquishing the profile makes the uniform scenario
        // authoritative on the actual next native step.
        flow.clear_pressure_profile(Face3::ZMax);
        lattice3::rebind_pressure_face3(
            &mut flow, Face3::ZMax, gauge(&dormant), &frames, 1.5, &scaling, &gate,
        ).unwrap();
        flow.step_serial(cx, None, [0.0; 3]).unwrap();
        assert_eq!(flow.completed_steps(), 3);
        for y in 0..4 {
            for x in 0..4 {
                near(flow.macroscopic((4 * 4 + y) * 4 + x).unwrap().0, density(-80.0));
            }
        }
    });
}

#[test]
fn original_generic_d3q19_and_d3q27_share_pressure_history_and_atomic_refusal() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let tau = scaling.flow().tau;
    check_generic::<D3q19, 19>(CollisionModel3::Bgk { tau });
    check_generic::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk { tau, equilibrium: Equilibrium27::Factorized },
        maximum_mach: 0.2,
    });
}


fn typed_pressure(outside: fs_scenario::payload::OutsideDomainPolicy) -> BoundaryCondition {
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    use fs_scenario::payload::{OrientationParity, Payload, PayloadId, PayloadMeta,
        QuantityContract, ReferenceSemantics, SampleSource, ScalarPayload, TableInterpolation};
    boundary(BcValue::Typed(Payload::Scalar(ScalarPayload::new(
        PayloadMeta::new(
            QuantityContract::Semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Instantaneous)),
            PayloadId::new("basis/scalar").unwrap(), WORLD, OrientationParity::Even,
            ReferenceSemantics::Continuous,
        ).unwrap(),
        SampleSource::table(
            [0.0, 0.002, 0.004].map(|t| QtyAny::new(t, Dims([0, 0, 1, 0, 0, 0]))).to_vec(),
            [80.0, -80.0, 40.0].map(|p| QtyAny::new(p, PRESSURE)).to_vec(),
            TableInterpolation::Linear, outside,
        ).unwrap(),
    ).unwrap())))
}

fn check_typed_native<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    use fs_couple::scenario_open_boundary::{ScenarioBoundary3, lattice3::step_serial};
    use fs_scenario::payload::OutsideDomainPolicy;
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    let source = typed_pressure(OutsideDomainPolicy::Refuse);
    let legacy = boundary(BcValue::Signal(TimeSignal::Table {
        times: vec![0.0, 0.002, 0.004], values: vec![80.0, -80.0, 40.0],
        dims: PRESSURE, interp: Interp::Linear,
    }));
    with_cx(&gate, |cx| {
        let make_native = || Grid3::<L, Q>::try_new_equilibrium(
            [4, 4, 5], collision, [0.0; 3], [true, true, false],
            GridBudget3 { maximum_cells: 80, maximum_mach: 0.2 },
            |_| false, |_| (1.0, [0.0; 3]),
        ).unwrap();
        let mut typed = OpenBoundaryGrid3::try_new(make_native(), cx, faces()).unwrap();
        let mut control = OpenBoundaryGrid3::try_new(make_native(), cx, faces()).unwrap();
        for (time, pressure) in [(0.001, 0.0), (0.002, -80.0), (0.003, -20.0), (0.004, 40.0)] {
            for (grid, row) in [(&mut typed, &source), (&mut control, &legacy)] {
                step_serial(grid,
                    [None, None, None, None, Some(ScenarioBoundary3::Pressure(gauge(row))), None],
                    &frames, time, &scaling, cx, None, [0.0; 3]).unwrap();
            }
            assert_eq!(typed.populations(), control.populations());
            near(typed.macroscopic(5).unwrap().0, density(pressure));
        }
        let before = typed.populations().to_vec();
        let targets = typed.boundaries();
        let steps = typed.completed_steps();
        assert!(step_serial(&mut typed,
            [None, None, None, None, Some(ScenarioBoundary3::Pressure(gauge(&source))), None],
            &frames, 0.005, &scaling, cx, None, [0.0; 3]).is_err());
        assert_eq!(typed.populations(), before);
        assert_eq!(typed.boundaries(), targets);
        assert_eq!(typed.completed_steps(), steps);
        // Refusal did not consume a source sample or advance either native clock.
        for (grid, row) in [(&mut typed, &source), (&mut control, &legacy)] {
            step_serial(grid,
                [None, None, None, None, Some(ScenarioBoundary3::Pressure(gauge(row))), None],
                &frames, 0.004, &scaling, cx, None, [0.0; 3]).unwrap();
        }
        assert_eq!(typed.populations(), control.populations());
    });
}

#[test]
fn g3_typed_pressure_advances_original_d3q19_and_d3q27_with_exact_legacy_equivalence() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let tau = scaling.flow().tau;
    check_typed_native::<D3q19, 19>(CollisionModel3::Bgk { tau });
    check_typed_native::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk { tau, equilibrium: Equilibrium27::Factorized },
        maximum_mach: 0.2,
    });
}
