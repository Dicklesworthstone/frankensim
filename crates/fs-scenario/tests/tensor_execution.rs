//! G0/G3: existing tensor histories retain both Cartesian transformation indices.
use fs_ga::{Quat, Vec3};
use fs_qty::{Dims, QtyAny};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_scenario::{Frame, FrameId, FrameMotion, FrameTree, WORLD};
use fs_scenario::payload::{
    OrientationParity, OutsideDomainPolicy, Payload, PayloadId, PayloadMeta,
    QuantityContract, ReferenceSemantics, SampleSource, TableInterpolation,
    TensorPayload, canonical_payload_bytes, decode_payload,
};
use fs_scenario::tensor::TensorEvaluationError;
use std::ops::ControlFlow;

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
fn meta(frame: FrameId) -> PayloadMeta {
    PayloadMeta::new(
        QuantityContract::Semantic(SemanticType::new(QuantityKind::Stress, ValueForm::Instantaneous)),
        PayloadId::new("basis/frame-cartesian").unwrap(), frame,
        OrientationParity::Even, ReferenceSemantics::Continuous,
    ).unwrap()
}
fn entries(values: [[f64; 3]; 3]) -> Vec<QtyAny> {
    values.into_iter().flatten().map(|value| QtyAny::new(value, PA)).collect()
}
fn fixed(values: [[f64; 3]; 3], frame: FrameId) -> TensorPayload {
    TensorPayload::new(meta(frame), 3, 3, SampleSource::fixed(entries(values))).unwrap()
}
fn sample(value: &TensorPayload, frames: &FrameTree, time: f64) -> [[f64; 3]; 3] {
    value.world_tensor_at(frames, time, &mut || ControlFlow::Continue(())).unwrap().values
}
fn close(actual: [[f64; 3]; 3], expected: [[f64; 3]; 3]) {
    for (a, b) in actual.into_iter().flatten().zip(expected.into_iter().flatten()) {
        assert!((a-b).abs() < 1e-12, "{actual:?} != {expected:?}");
    }
}

#[test]
fn g3_two_index_rotation_preserves_symmetry_and_does_not_add_translation() {
    let mut frames = FrameTree::new();
    frames.add(Frame { id: FrameId(1), name: "rotated".into(), parent: WORLD,
        motion: FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), std::f64::consts::FRAC_PI_2),
            translation: Vec3::new(1e12, -1e12, 1e12),
        } });
    let value = fixed([[2.0, 1.0, 3.0], [1.0, 4.0, 5.0], [3.0, 5.0, 6.0]], FrameId(1));
    let result = value.world_tensor_at(&frames, 0.0, &mut || ControlFlow::Continue(())).unwrap();
    close(result.values, [[4.0, -1.0, -5.0], [-1.0, 2.0, 3.0], [-5.0, 3.0, 6.0]]);
    assert_eq!(result.contract, value.meta().contract());
    for i in 0..3 { for j in 0..3 { assert_eq!(result.values[i][j], result.values[j][i]); } }
    let unsymmetric = fixed([[0.0, 1.0, 0.0], [2.0, 0.0, 0.0], [0.0; 3]], FrameId(1));
    close(sample(&unsymmetric, &frames, 0.0), [[0.0, -2.0, 0.0], [-1.0, 0.0, 0.0], [0.0; 3]]);
}

#[test]
fn g3_canonical_tensor_history_reuses_sampling_and_keeps_physical_frame_time() {
    let mut frames = FrameTree::new();
    frames.add(Frame { id: FrameId(1), name: "rotation".into(), parent: WORLD,
        motion: FrameMotion::Rotating { axis: [0.0, 0.0, 1.0], center: Vec3::new(10.0, 20.0, 30.0),
            rate: QtyAny::new(std::f64::consts::FRAC_PI_2, Dims([0, 0, -1, 0, 0, 0])) } });
    let start = [[2.0, 0.0, 0.0], [0.0; 3], [0.0; 3]];
    let end = [[4.0, 0.0, 0.0], [0.0; 3], [0.0; 3]];
    for outside in [OutsideDomainPolicy::Clamp, OutsideDomainPolicy::Periodic, OutsideDomainPolicy::Refuse] {
        let original = Payload::Tensor(TensorPayload::new(meta(FrameId(1)), 3, 3,
            SampleSource::table(vec![QtyAny::new(0.0, TIME), QtyAny::new(1.0, TIME)],
                vec![entries(start), entries(end)], TableInterpolation::Linear, outside).unwrap()).unwrap());
        let decoded = decode_payload(&canonical_payload_bytes(&original)).unwrap();
        assert_eq!(decoded, original);
        let Payload::Tensor(value) = decoded else { panic!("tensor codec changed shape") };
        close(sample(&value, &frames, 0.5), [[1.5, 1.5, 0.0], [1.5, 1.5, 0.0], [0.0; 3]]);
        close(sample(&value, &frames, 1.0), [[0.0; 3], [0.0, 4.0, 0.0], [0.0; 3]]);
        match outside {
            OutsideDomainPolicy::Clamp => close(sample(&value, &frames, 2.0), end),
            OutsideDomainPolicy::Periodic => close(sample(&value, &frames, 2.0), start),
            OutsideDomainPolicy::Refuse => assert!(matches!(value.world_tensor_at(&frames, 2.0,
                &mut || ControlFlow::Continue(())), Err(TensorEvaluationError::OutsideDomain { .. }))),
        }
    }
}

#[test]
fn g0_shape_basis_reference_and_g4_cancel_refuse_without_samples() {
    let frames = FrameTree::new();
    let good = fixed([[1.0; 3]; 3], WORLD);
    for stop in 1..=4 {
        let mut count = 0;
        assert_eq!(good.world_tensor_at(&frames, 0.0, &mut || {
            count += 1;
            if count == stop { ControlFlow::Break(()) } else { ControlFlow::Continue(()) }
        }), Err(TensorEvaluationError::Cancelled));
    }
    let wrong_shape = TensorPayload::new(meta(WORLD), 1, 9, good.source().clone()).unwrap();
    assert!(wrong_shape.world_tensor_at(&frames, 0.0, &mut || ControlFlow::Continue(())).is_err());
    for (basis, reference) in [
        ("basis/material", ReferenceSemantics::Continuous),
        ("basis/frame-cartesian", ReferenceSemantics::ResetAtEvent(PayloadId::new("impact").unwrap())),
    ] {
        let metadata = PayloadMeta::new(good.meta().contract(), PayloadId::new(basis).unwrap(), WORLD,
            OrientationParity::Even, reference).unwrap();
        let bad = TensorPayload::new(metadata, 3, 3, good.source().clone()).unwrap();
        assert!(bad.world_tensor_at(&frames, 0.0, &mut || ControlFlow::Continue(())).is_err());
    }
    assert_eq!(good.world_tensor_at(&frames, f64::NAN, &mut || ControlFlow::Continue(())),
        Err(TensorEvaluationError::InvalidTime));
    assert_eq!(sample(&good, &frames, 0.0), [[1.0; 3]; 3]);
}
