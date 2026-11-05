//! G0/G3: execute the existing typed payload, not an independent sample schema.
use fs_ga::{Quat, Vec3};
use fs_qty::{Dims, QtyAny};
use fs_scenario::{BcKind, BcValue, BoundaryCondition, Environment, Frame, FrameId,
    FrameMotion, FrameTree, Physics, Scenario, WORLD};
use fs_scenario::ir::{parse_ir, write_ir};
use fs_scenario::payload::{
    DistributionFamily, OrientationParity, OutsideDomainPolicy, Payload, PayloadId, PayloadMeta,
    QuantityContract, ReferenceSemantics, SampleSource, TableInterpolation, VectorPayload,
};
use fs_scenario::vector::VectorEvaluationError;
use std::ops::ControlFlow;

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
fn meta(frame: FrameId, basis: &str, reference: ReferenceSemantics) -> PayloadMeta {
    PayloadMeta::new(QuantityContract::Dimensions(LENGTH), PayloadId::new(basis).unwrap(),
        frame, OrientationParity::Even, reference).unwrap()
}
fn values(v: [f64; 3]) -> Vec<QtyAny> { v.map(|x| QtyAny::new(x, LENGTH)).to_vec() }
fn fixed(v: [f64; 3], frame: FrameId) -> VectorPayload {
    VectorPayload::new(meta(frame, "basis/frame-cartesian", ReferenceSemantics::Continuous),
        SampleSource::fixed(values(v))).unwrap()
}
fn table(interpolation: TableInterpolation, outside: OutsideDomainPolicy) -> VectorPayload {
    VectorPayload::new(meta(WORLD, "basis/world-cartesian", ReferenceSemantics::Continuous),
        SampleSource::table(vec![QtyAny::new(1.0, TIME), QtyAny::new(2.0, TIME), QtyAny::new(3.0, TIME)],
            vec![values([0.0, 2.0, 4.0]), values([2.0, 4.0, 6.0]), values([4.0, 6.0, 8.0])],
            interpolation, outside).unwrap()).unwrap()
}
fn sample(payload: &VectorPayload, frames: &FrameTree, time: f64) -> [f64; 3] {
    payload.world_vector_at(frames, time, &mut || ControlFlow::Continue(())).unwrap().values
}
fn near(a: [f64; 3], b: [f64; 3]) {
    for i in 0..3 { assert!((a[i] - b[i]).abs() < 1e-10, "{a:?} != {b:?}"); }
}

#[test]
fn linear_step_and_exact_knots_obey_declared_semantics() {
    let frames = FrameTree::new();
    let linear = table(TableInterpolation::Linear, OutsideDomainPolicy::Refuse);
    assert_eq!(sample(&linear, &frames, 1.5), [1.0, 3.0, 5.0]);
    assert_eq!(sample(&linear, &frames, 2.0), [2.0, 4.0, 6.0]);
    assert_eq!(sample(&linear, &frames, 3.0), [4.0, 6.0, 8.0]);
    let held = table(TableInterpolation::StepLeft, OutsideDomainPolicy::Refuse);
    assert_eq!(sample(&held, &frames, 1.999), [0.0, 2.0, 4.0]);
    assert_eq!(sample(&held, &frames, 2.0), [2.0, 4.0, 6.0]);
    assert_eq!(sample(&held, &frames, 3.0), [4.0, 6.0, 8.0]);
    let result = linear.world_vector_at(&frames, 1.5, &mut || ControlFlow::Continue(())).unwrap();
    assert_eq!(result.contract, QuantityContract::Dimensions(LENGTH));
    assert_eq!(result.orientation, OrientationParity::Even);
}

#[test]
fn outside_refuse_clamp_periodic_and_singleton_are_distinct() {
    let frames = FrameTree::new();
    let refused = table(TableInterpolation::Linear, OutsideDomainPolicy::Refuse);
    for t in [0.99, 3.01] {
        assert!(matches!(refused.world_vector_at(&frames, t, &mut || ControlFlow::Continue(())),
            Err(VectorEvaluationError::OutsideDomain { .. })));
    }
    let clamp = table(TableInterpolation::Linear, OutsideDomainPolicy::Clamp);
    assert_eq!(sample(&clamp, &frames, -99.0), [0.0, 2.0, 4.0]);
    assert_eq!(sample(&clamp, &frames, 99.0), [4.0, 6.0, 8.0]);
    let repeat = table(TableInterpolation::Linear, OutsideDomainPolicy::Periodic);
    assert_eq!(sample(&repeat, &frames, 3.0), [4.0, 6.0, 8.0]);
    assert_eq!(sample(&repeat, &frames, 5.0), [0.0, 2.0, 4.0]);
    assert_eq!(sample(&repeat, &frames, 3.5), [1.0, 3.0, 5.0]);
    assert_eq!(sample(&repeat, &frames, -0.5), [1.0, 3.0, 5.0]);
    let singleton = VectorPayload::new(meta(WORLD, "basis/world-cartesian", ReferenceSemantics::Continuous),
        SampleSource::table(vec![QtyAny::new(2.0, TIME)], vec![values([1.0; 3])],
            TableInterpolation::Linear, OutsideDomainPolicy::Clamp).unwrap()).unwrap();
    assert_eq!(sample(&singleton, &frames, -10.0), [1.0; 3]);
}

#[test]
fn nested_frames_rotate_without_adding_or_subtracting_translation() {
    let mut frames = FrameTree::new();
    frames.add(Frame { id: FrameId(1), name: "parent".into(), parent: WORLD,
        motion: FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), std::f64::consts::FRAC_PI_2),
            translation: Vec3::new(1e12, -1e12, 1e12),
        } });
    frames.add(Frame { id: FrameId(2), name: "child".into(), parent: FrameId(1),
        motion: FrameMotion::Fixed { orientation: Quat::identity(), translation: Vec3::new(10.0, 20.0, 30.0) } });
    near(sample(&fixed([2.0, 0.0, 0.0], FrameId(2)), &frames, 0.0), [0.0, 2.0, 0.0]);
    assert!(fixed([1.0; 3], FrameId(99)).world_vector_at(&frames, 0.0,
        &mut || ControlFlow::Continue(())).is_err());
    frames.frames[0].parent = FrameId(2);
    assert!(fixed([1.0; 3], FrameId(2)).world_vector_at(&frames, 0.0,
        &mut || ControlFlow::Continue(())).is_err());
}

#[test]
fn clamped_value_keeps_rotating_at_requested_physical_time() {
    let mut frames = FrameTree::new();
    frames.add(Frame { id: FrameId(1), name: "rotating".into(), parent: WORLD,
        motion: FrameMotion::Rotating { axis: [0.0, 0.0, 1.0], center: Vec3::new(10.0, 20.0, 30.0),
            rate: QtyAny::new(std::f64::consts::FRAC_PI_2, Dims([0, 0, -1, 0, 0, 0])) } });
    let payload = VectorPayload::new(meta(FrameId(1), "basis/frame-cartesian", ReferenceSemantics::Continuous),
        SampleSource::table(vec![QtyAny::new(0.0, TIME), QtyAny::new(1.0, TIME)],
            vec![values([1.0, 0.0, 0.0]); 2], TableInterpolation::Linear, OutsideDomainPolicy::Clamp).unwrap()).unwrap();
    near(sample(&payload, &frames, 1.0), [0.0, 1.0, 0.0]);
    near(sample(&payload, &frames, 2.0), [-1.0, 0.0, 0.0]);
}

#[test]
fn unsupported_basis_reference_shape_and_distribution_never_become_vectors() {
    let frames = FrameTree::new();
    let cases = [
        VectorPayload::new(meta(WORLD, "basis/material", ReferenceSemantics::Continuous), SampleSource::fixed(values([1.0; 3]))).unwrap(),
        VectorPayload::new(meta(FrameId(1), "basis/world-cartesian", ReferenceSemantics::Continuous), SampleSource::fixed(values([1.0; 3]))).unwrap(),
        VectorPayload::new(meta(WORLD, "basis/world-cartesian", ReferenceSemantics::ResetAtEvent(PayloadId::new("impact").unwrap())), SampleSource::fixed(values([1.0; 3]))).unwrap(),
        VectorPayload::new(meta(WORLD, "basis/world-cartesian", ReferenceSemantics::Continuous), SampleSource::fixed(vec![QtyAny::new(1.0, LENGTH); 2])).unwrap(),
        VectorPayload::new(meta(WORLD, "basis/world-cartesian", ReferenceSemantics::Continuous), SampleSource::distribution(DistributionFamily::Empirical, vec![values([1.0; 3])]).unwrap()).unwrap(),
    ];
    for payload in cases {
        assert!(matches!(payload.world_vector_at(&frames, 0.0, &mut || ControlFlow::Continue(())),
            Err(VectorEvaluationError::Unsupported(_))));
    }
}

#[test]
fn convex_interpolation_avoids_opposite_sign_sample_difference_overflow() {
    let payload = VectorPayload::new(meta(WORLD, "basis/world-cartesian", ReferenceSemantics::Continuous),
        SampleSource::table(vec![QtyAny::new(0.0, TIME), QtyAny::new(2.0, TIME)],
            vec![values([-1e308, 0.0, 0.0]), values([1e308, 0.0, 0.0])],
            TableInterpolation::Linear, OutsideDomainPolicy::Refuse).unwrap()).unwrap();
    assert_eq!(sample(&payload, &FrameTree::new(), 1.0), [0.0; 3]);
}

#[test]
fn cancellation_at_each_publication_boundary_and_nonfinite_times_refuse() {
    let payload = fixed([1.0; 3], WORLD);
    let frames = FrameTree::new();
    for stop in 1..=4 {
        let mut count = 0;
        let result = payload.world_vector_at(&frames, 0.0, &mut || {
            count += 1;
            if count == stop { ControlFlow::Break(()) } else { ControlFlow::Continue(()) }
        });
        assert_eq!(result, Err(VectorEvaluationError::Cancelled));
    }
    for time in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(payload.world_vector_at(&frames, time, &mut || ControlFlow::Continue(())),
            Err(VectorEvaluationError::InvalidTime));
    }
    assert_eq!(sample(&payload, &frames, 0.0), [1.0; 3]);
}

const VELOCITY: Dims = Dims([1, 0, -1, 0, 0, 0]);
fn velocity_values(value: [f64; 3]) -> Vec<QtyAny> {
    value.map(|x| QtyAny::new(x, VELOCITY)).to_vec()
}
fn velocity_condition(source: SampleSource<Vec<QtyAny>>, frame: FrameId) -> BoundaryCondition {
    let meta = PayloadMeta::new(QuantityContract::Dimensions(VELOCITY),
        PayloadId::new("basis/frame-cartesian").unwrap(), frame,
        OrientationParity::Even, ReferenceSemantics::Continuous).unwrap();
    BoundaryCondition { region: "inlet".into(), physics: Physics::IncompressibleFlow,
        kind: BcKind::Dirichlet,
        value: Some(BcValue::Typed(Payload::Vector(VectorPayload::new(meta, source).unwrap()))),
        compatibility: None, frame: frame.0 }
}
fn velocity_sample(row: &BoundaryCondition, frames: &FrameTree, time: f64) -> [f64; 3] {
    row.velocity_world_at(frames, time, &mut || ControlFlow::Continue(())).unwrap()
}

#[test]
fn g3_velocity_boundary_ir_preserves_history_and_physical_component_frame_time() {
    let mut scenario = Scenario::new("vector-velocity", 7, Environment::earth_lab());
    scenario.frames.add(Frame { id: FrameId(1), name: "rotating-axes".into(), parent: WORLD,
        motion: FrameMotion::Rotating { axis: [0.0, 0.0, 1.0],
            center: Vec3::new(10.0, 20.0, 30.0),
            rate: QtyAny::new(std::f64::consts::FRAC_PI_2, Dims([0, 0, -1, 0, 0, 0])) } });
    let source = SampleSource::table(vec![QtyAny::new(0.0, TIME), QtyAny::new(1.0, TIME)],
        vec![velocity_values([1.0, 0.0, 0.0]), velocity_values([2.0, 0.0, 0.0])],
        TableInterpolation::Linear, OutsideDomainPolicy::Clamp).unwrap();
    scenario.base_bcs.push(velocity_condition(source, FrameId(1)));
    assert!(scenario.validate().is_empty());
    let decoded = parse_ir(&write_ir(&scenario)).unwrap();
    let restored = decoded.scenario();
    assert_eq!(restored, &scenario);
    let row = &restored.base_bcs[0];
    near(velocity_sample(row, &restored.frames, 0.0), [1.0, 0.0, 0.0]);
    let diagonal = 1.5 * std::f64::consts::FRAC_1_SQRT_2;
    near(velocity_sample(row, &restored.frames, 0.5), [diagonal, diagonal, 0.0]);
    near(velocity_sample(row, &restored.frames, 1.0), [0.0, 2.0, 0.0]);
    near(velocity_sample(row, &restored.frames, 2.0), [-2.0, 0.0, 0.0]);
    // A component frame is not a material surface: zero absolute velocity
    // remains zero even with rotating axes about an offset origin.
    let stationary = velocity_condition(SampleSource::fixed(velocity_values([0.0; 3])), FrameId(1));
    near(velocity_sample(&stationary, &restored.frames, 2.0), [0.0; 3]);
}

#[test]
fn g0_velocity_boundary_does_not_infer_direction_or_reinterpret_a_different_role() {
    let frames = FrameTree::new();
    let original = velocity_condition(SampleSource::fixed(velocity_values([1.0, -2.0, 3.0])), WORLD);
    assert_eq!(velocity_sample(&original, &frames, 0.0), [1.0, -2.0, 3.0]);
    let mut bad = original.clone();
    bad.frame = 1;
    assert!(bad.velocity_world_at(&frames, 0.0, &mut || ControlFlow::Continue(())).is_err());
    bad = original.clone();
    bad.value = Some(BcValue::Uniform(QtyAny::new(2.0, VELOCITY)));
    assert!(matches!(bad.velocity_world_at(&frames, 0.0, &mut || ControlFlow::Continue(())),
        Err(VectorEvaluationError::Unsupported(_))));
    bad = original.clone();
    bad.kind = BcKind::Neumann;
    assert!(bad.velocity_world_at(&frames, 0.0, &mut || ControlFlow::Continue(())).is_err());
    bad = original.clone();
    bad.physics = Physics::Elasticity;
    assert!(bad.velocity_world_at(&frames, 0.0, &mut || ControlFlow::Continue(())).is_err());
    bad = original.clone();
    bad.value = Some(BcValue::Typed(Payload::Vector(fixed([1.0; 3], WORLD))));
    assert!(bad.velocity_world_at(&frames, 0.0, &mut || ControlFlow::Continue(())).is_err());
    let odd = PayloadMeta::new(QuantityContract::Dimensions(VELOCITY),
        PayloadId::new("basis/world-cartesian").unwrap(), WORLD,
        OrientationParity::Odd, ReferenceSemantics::Continuous).unwrap();
    bad.value = Some(BcValue::Typed(Payload::Vector(VectorPayload::new(odd,
        SampleSource::fixed(velocity_values([1.0; 3]))).unwrap())));
    assert!(bad.velocity_world_at(&frames, 0.0, &mut || ControlFlow::Continue(())).is_err());
    bad = velocity_condition(SampleSource::fixed(vec![QtyAny::new(1.0, VELOCITY); 2]), WORLD);
    assert!(bad.velocity_world_at(&frames, 0.0, &mut || ControlFlow::Continue(())).is_err());
    assert_eq!(original.velocity_world_at(&frames, 0.0, &mut || ControlFlow::Break(())),
        Err(VectorEvaluationError::Cancelled));
    assert_eq!(original.velocity_world_at(&frames, f64::NAN, &mut || ControlFlow::Continue(())),
        Err(VectorEvaluationError::InvalidTime));
    assert_eq!(velocity_sample(&original, &frames, 0.0), [1.0, -2.0, 3.0]);
}
