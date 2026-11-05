//! G0/G3: typed SI histories drive the native force-corrected velocity inlet.

#[path = "scenario_velocity/lattice3.rs"]
mod lattice3;

use fs_couple::scenario_velocity::{
    ScenarioVelocityError, rebind_velocity_face3, rebind_velocity_faces3, sample_velocity3,
};
use fs_exec::CancelGate;
use fs_ga::{Quat, Vec3};
use fs_lbm::d3q19::{BoundaryGrid3, BoundarySpec3, Face3, FaceBoundary3, InletFlux3};
use fs_lbm::material::plan_scaling_from_material;
use fs_material::fluid_transport::FluidTransportState;
use fs_qty::{Dims, QtyAny};
use fs_scenario::payload::{
    OrientationParity, OutsideDomainPolicy, Payload, PayloadId, PayloadMeta, QuantityContract,
    ReferenceSemantics, SampleSource, TableInterpolation, VectorPayload,
};
use fs_scenario::vector::VectorEvaluationError;
use fs_scenario::{
    BcKind, BcValue, BoundaryCondition, Frame, FrameId, FrameMotion, FrameTree, Physics, WORLD,
};
use std::ops::ControlFlow;

const VELOCITY: Dims = Dims([1, 0, -1, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);

fn values(vector: [f64; 3]) -> Vec<QtyAny> {
    vector.map(|v| QtyAny::new(v, VELOCITY)).to_vec()
}

fn boundary(first: [f64; 3], last: [f64; 3], frame: FrameId) -> BoundaryCondition {
    let meta = PayloadMeta::new(
        QuantityContract::Dimensions(VELOCITY),
        PayloadId::new("basis/frame-cartesian").unwrap(),
        frame,
        OrientationParity::Even,
        ReferenceSemantics::Continuous,
    )
    .unwrap();
    let source = SampleSource::table(
        vec![QtyAny::new(0.0, TIME), QtyAny::new(1.0, TIME)],
        vec![values(first), values(last)],
        TableInterpolation::Linear,
        OutsideDomainPolicy::Refuse,
    )
    .unwrap();
    BoundaryCondition {
        region: "inlet".into(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::Dirichlet,
        value: Some(BcValue::Typed(Payload::Vector(
            VectorPayload::new(meta, source).unwrap(),
        ))),
        compatibility: None,
        frame: frame.0,
    }
}

fn fluid() -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 4200.0, 0.6, 0.001).unwrap()
}

fn grid(tau: f64) -> BoundaryGrid3 {
    BoundaryGrid3::new(
        8,
        8,
        8,
        tau,
        [1e-5, -2e-5, 3e-5],
        BoundarySpec3::new([
            FaceBoundary3::Periodic,
            FaceBoundary3::Periodic,
            FaceBoundary3::Periodic,
            FaceBoundary3::Periodic,
            FaceBoundary3::Velocity { velocity: [0.0; 3] },
            FaceBoundary3::Pressure { density: 1.0 },
        ]),
    )
}

fn paired_grid(tau: f64) -> BoundaryGrid3 {
    BoundaryGrid3::new(
        8,
        8,
        8,
        tau,
        [1e-5, -2e-5, 3e-5],
        BoundarySpec3::new([
            FaceBoundary3::Periodic,
            FaceBoundary3::Periodic,
            FaceBoundary3::Periodic,
            FaceBoundary3::Periodic,
            FaceBoundary3::Velocity { velocity: [0.0; 3] },
            FaceBoundary3::Velocity { velocity: [0.0; 3] },
        ]),
    )
}

fn near(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    for axis in 0..3 {
        assert!(
            (actual[axis] - expected[axis]).abs() <= tolerance,
            "{actual:?} != {expected:?}"
        );
    }
}

#[test]
fn rotated_signed_history_reaches_native_inlet_in_physical_force_coordinates() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut flow = grid(scaling.flow().tau);
    let mut frames = FrameTree::new();
    frames.add(Frame {
        id: FrameId(1),
        name: "inlet_axes".into(),
        parent: WORLD,
        motion: FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(
                Vec3::new(0.0, 0.0, 1.0),
                std::f64::consts::FRAC_PI_2,
            ),
            translation: Vec3::new(1e6, -2e6, 3e6),
        },
    });
    let inlet = boundary([0.4, -0.2, 0.6], [-0.2, 0.4, -0.6], FrameId(1));
    let gate = CancelGate::new_clock_free();
    for (time_s, expected_si) in [
        (0.0, [0.2, 0.4, 0.6]),
        (0.5, [-0.1, 0.1, 0.0]),
        (1.0, [-0.4, -0.2, -0.6]),
    ] {
        let before = flow.populations(3, 3, 0);
        let sample = rebind_velocity_face3(
            &mut flow,
            Face3::ZMin,
            &inlet,
            &frames,
            time_s,
            &scaling,
            &gate,
        )
        .unwrap();
        near(sample.world_m_s, expected_si, 1e-14);
        near(sample.lattice, expected_si.map(|u| u / 20.0), 1e-15);
        // A control update publishes a target, not a hidden physics step.
        assert_eq!(flow.populations(3, 3, 0), before);
        flow.step();
        for y in 0..8 {
            for x in 0..8 {
                near(flow.velocity(x, y, 0), sample.lattice, 2e-14);
            }
        }
        assert_eq!(
            flow.boundaries().face(Face3::ZMax),
            FaceBoundary3::Pressure { density: 1.0 }
        );
    }
}

#[test]
fn full_vector_speed_refusal_and_native_kind_refusal_publish_no_target() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut flow = grid(scaling.flow().tau);
    let before = flow.boundaries();
    let populations = flow.populations(3, 3, 0);
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    // Every component is below 1 m/s; their vector speed is above it.
    let too_fast = boundary([0.7; 3], [0.7; 3], WORLD);
    assert!(matches!(
        rebind_velocity_face3(
            &mut flow,
            Face3::ZMin,
            &too_fast,
            &frames,
            0.0,
            &scaling,
            &gate,
        ),
        Err(ScenarioVelocityError::SpeedBound { .. })
    ));
    let inlet = boundary([0.0, 0.0, 0.5], [0.0, 0.0, 0.5], WORLD);
    assert!(matches!(
        rebind_velocity_face3(
            &mut flow,
            Face3::ZMax,
            &inlet,
            &frames,
            0.0,
            &scaling,
            &gate,
        ),
        Err(ScenarioVelocityError::Boundary(_))
    ));
    assert_eq!(flow.boundaries(), before);
    assert_eq!(flow.populations(3, 3, 0), populations);
}

#[test]
fn installed_profile_refuses_even_an_unchanged_uniform_scenario_target() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut flow = grid(scaling.flow().tau);
    let gate = CancelGate::new_clock_free();
    let prescribed = [0.01, -0.005, 0.02];
    flow.set_velocity_profile(Face3::ZMin, InletFlux3::Unconstrained, &gate, |_| {
        Ok(prescribed)
    })
    .unwrap();
    let profile = flow.velocity_profile(Face3::ZMin).unwrap().to_vec();
    let before = flow.boundaries();
    let populations = flow.populations(3, 3, 0);
    // This equals the stored uniform target. The native rebind alone would
    // take its unchanged-target shortcut while the profile stayed active.
    let inlet = boundary([0.0; 3], [0.0; 3], WORLD);
    assert!(matches!(
        rebind_velocity_face3(
            &mut flow,
            Face3::ZMin,
            &inlet,
            &FrameTree::new(),
            0.0,
            &scaling,
            &gate,
        ),
        Err(ScenarioVelocityError::Boundary(_))
    ));
    assert_eq!(flow.boundaries(), before);
    assert_eq!(flow.populations(3, 3, 0), populations);
    assert_eq!(flow.velocity_profile(Face3::ZMin).unwrap(), profile);
    flow.step();
    near(flow.velocity(3, 3, 0), prescribed, 2e-14);
}

#[test]
fn unsupported_scalar_time_and_cancellation_do_not_change_native_state() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut flow = grid(scaling.flow().tau);
    let frames = FrameTree::new();
    let mut inlet = boundary([0.0, 0.0, 0.5], [0.0, 0.0, 0.5], WORLD);
    let before = flow.boundaries();
    let populations = flow.populations(3, 3, 0);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        rebind_velocity_face3(
            &mut flow,
            Face3::ZMin,
            &inlet,
            &frames,
            0.0,
            &scaling,
            &cancelled,
        ),
        Err(ScenarioVelocityError::Scenario(
            VectorEvaluationError::Cancelled
        ))
    ));
    for time_s in [-1.0, f64::NAN] {
        assert!(
            rebind_velocity_face3(
                &mut flow,
                Face3::ZMin,
                &inlet,
                &frames,
                time_s,
                &scaling,
                &CancelGate::new_clock_free(),
            )
            .is_err()
        );
    }
    inlet.value = Some(BcValue::Uniform(QtyAny::new(0.5, VELOCITY)));
    assert!(matches!(
        sample_velocity3(&inlet, &frames, 0.0, &scaling, &mut || {
            ControlFlow::Continue(())
        }),
        Err(ScenarioVelocityError::Scenario(
            VectorEvaluationError::Unsupported(_)
        ))
    ));
    assert_eq!(flow.boundaries(), before);
    assert_eq!(flow.populations(3, 3, 0), populations);
}

#[test]
fn paired_histories_publish_one_time_on_both_native_velocity_faces() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut flow = paired_grid(scaling.flow().tau);
    let mut frames = FrameTree::new();
    frames.add(Frame {
        id: FrameId(1),
        name: "outlet_axes".into(),
        parent: WORLD,
        motion: FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(
                Vec3::new(0.0, 0.0, 1.0),
                std::f64::consts::FRAC_PI_2,
            ),
            translation: Vec3::new(-1e6, 2e6, -3e6),
        },
    });
    let inlet = boundary([0.4, 0.1, 0.6], [-0.2, 0.3, -0.6], WORLD);
    let mut outlet = boundary([0.1, 0.2, 0.4], [0.3, -0.2, -0.4], FrameId(1));
    outlet.region = "outlet".into();
    let gate = CancelGate::new_clock_free();
    for (time_s, expected) in [
        (0.0, [[0.4, 0.1, 0.6], [-0.2, 0.1, 0.4]]),
        (0.5, [[0.1, 0.2, 0.0], [0.0, 0.2, 0.0]]),
        (1.0, [[-0.2, 0.3, -0.6], [0.2, 0.3, -0.4]]),
    ] {
        let before = (flow.populations(3, 3, 0), flow.populations(3, 3, 7));
        let samples = rebind_velocity_faces3(
            &mut flow,
            [None, None, None, None, Some(&inlet), Some(&outlet)],
            &frames,
            time_s,
            &scaling,
            &gate,
        )
        .unwrap();
        assert!(samples[..4].iter().all(Option::is_none));
        for (sample, expected_si) in samples[4..].iter().zip(expected) {
            let sample = sample.unwrap();
            assert_eq!(sample.time_s, time_s);
            near(sample.world_m_s, expected_si, 1e-14);
            near(sample.lattice, expected_si.map(|u| u / 20.0), 1e-15);
        }
        assert_eq!(
            (flow.populations(3, 3, 0), flow.populations(3, 3, 7)),
            before
        );
        flow.step();
        for (z, expected_si) in [0, 7].into_iter().zip(expected) {
            for y in 0..8 {
                for x in 0..8 {
                    near(flow.velocity(x, y, z), expected_si.map(|u| u / 20.0), 3e-14);
                }
            }
        }
    }
}

#[test]
fn later_history_or_native_kind_refusal_keeps_all_face_targets() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut flow = paired_grid(scaling.flow().tau);
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    let inlet = boundary([0.0, 0.0, 0.5], [0.0, 0.0, 0.5], WORLD);
    let too_fast = boundary([0.7; 3], [0.7; 3], WORLD);
    let mut scalar = inlet.clone();
    scalar.value = Some(BcValue::Uniform(QtyAny::new(0.5, VELOCITY)));
    let missing_frame = boundary([0.0; 3], [0.0; 3], FrameId(99));
    let before = flow.boundaries();
    let populations = (flow.populations(3, 3, 0), flow.populations(3, 3, 7));
    for invalid in [&too_fast, &scalar, &missing_frame] {
        assert!(
            rebind_velocity_faces3(
                &mut flow,
                [None, None, None, None, Some(&inlet), Some(invalid)],
                &frames,
                0.5,
                &scaling,
                &gate,
            )
            .is_err()
        );
        assert_eq!(flow.boundaries(), before);
        assert_eq!(
            (flow.populations(3, 3, 0), flow.populations(3, 3, 7)),
            populations
        );
    }

    // Both histories sample successfully, but the second face is a pressure
    // boundary. Native admission must reject the whole prepared update.
    let mut pressure_outlet = grid(scaling.flow().tau);
    let before = pressure_outlet.boundaries();
    assert!(matches!(
        rebind_velocity_faces3(
            &mut pressure_outlet,
            [None, None, None, None, Some(&inlet), Some(&inlet)],
            &frames,
            0.5,
            &scaling,
            &gate,
        ),
        Err(ScenarioVelocityError::Boundary(_))
    ));
    assert_eq!(pressure_outlet.boundaries(), before);
}

#[test]
fn selected_profile_refuses_the_batch_and_unselected_profile_remains_active() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut flow = paired_grid(scaling.flow().tau);
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    let prescribed = [0.01, -0.005, 0.02];
    flow.set_velocity_profile(Face3::ZMax, InletFlux3::Unconstrained, &gate, |_| {
        Ok(prescribed)
    })
    .unwrap();
    let profile = flow.velocity_profile(Face3::ZMax).unwrap().to_vec();
    let before = flow.boundaries();
    let populations = (flow.populations(3, 3, 0), flow.populations(3, 3, 7));
    let inlet = boundary([0.1, 0.2, 0.3], [0.1, 0.2, 0.3], WORLD);
    let dormant_uniform = boundary([0.0; 3], [0.0; 3], WORLD);
    assert!(matches!(
        rebind_velocity_faces3(
            &mut flow,
            [None, None, None, None, Some(&inlet), Some(&dormant_uniform)],
            &frames,
            0.5,
            &scaling,
            &gate,
        ),
        Err(ScenarioVelocityError::Boundary(_))
    ));
    assert_eq!(flow.boundaries(), before);
    assert_eq!(
        (flow.populations(3, 3, 0), flow.populations(3, 3, 7)),
        populations
    );
    assert_eq!(flow.velocity_profile(Face3::ZMax).unwrap(), profile);

    let samples = rebind_velocity_faces3(
        &mut flow,
        [None, None, None, None, Some(&inlet), None],
        &frames,
        0.5,
        &scaling,
        &gate,
    )
    .unwrap();
    assert!(samples[5].is_none());
    assert_eq!(
        flow.boundaries().face(Face3::ZMax),
        before.face(Face3::ZMax)
    );
    assert_eq!(flow.velocity_profile(Face3::ZMax).unwrap(), profile);
    flow.step();
    near(flow.velocity(3, 3, 0), [0.005, 0.01, 0.015], 3e-14);
    near(flow.velocity(3, 3, 7), prescribed, 3e-14);
}

#[test]
fn empty_or_cancelled_velocity_batches_keep_native_state() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut flow = paired_grid(scaling.flow().tau);
    let frames = FrameTree::new();
    let before = flow.boundaries();
    let populations = (flow.populations(3, 3, 0), flow.populations(3, 3, 7));
    assert_eq!(
        rebind_velocity_faces3(
            &mut flow,
            [None; 6],
            &frames,
            0.5,
            &scaling,
            &CancelGate::new_clock_free(),
        )
        .unwrap(),
        [None; 6]
    );
    assert!(matches!(
        rebind_velocity_faces3(
            &mut flow,
            [None; 6],
            &frames,
            f64::NAN,
            &scaling,
            &CancelGate::new_clock_free(),
        ),
        Err(ScenarioVelocityError::Scenario(
            VectorEvaluationError::InvalidTime
        ))
    ));
    let inlet = boundary([0.0, 0.0, 0.5], [0.0, 0.0, 0.5], WORLD);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        rebind_velocity_faces3(
            &mut flow,
            [None, None, None, None, Some(&inlet), Some(&inlet)],
            &frames,
            0.5,
            &scaling,
            &cancelled,
        ),
        Err(ScenarioVelocityError::Scenario(
            VectorEvaluationError::Cancelled
        ))
    ));
    assert_eq!(flow.boundaries(), before);
    assert_eq!(
        (flow.populations(3, 3, 0), flow.populations(3, 3, 7)),
        populations
    );
}
