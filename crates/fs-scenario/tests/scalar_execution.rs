//! G0/G3: scalar histories share time semantics with existing vector consumers.
use fs_ga::{Quat, Vec3};
use fs_qty::{Dims, QtyAny};
use fs_scenario::payload::{
    DistributionFamily, OrientationParity, OutsideDomainPolicy, PayloadId, PayloadMeta,
    QuantityContract, ReferenceSemantics, SampleSource, ScalarPayload, TableInterpolation,
    VectorPayload,
};
use fs_scenario::scalar::ScalarEvaluationError;
use fs_scenario::{Frame, FrameId, FrameMotion, FrameTree, WORLD};
use std::ops::ControlFlow;

const VOLTAGE: Dims = Dims([2, 1, -3, 0, -1, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
fn meta(frame: FrameId, parity: OrientationParity, reference: ReferenceSemantics) -> PayloadMeta {
    PayloadMeta::new(
        QuantityContract::Dimensions(VOLTAGE),
        PayloadId::new("basis/frame-cartesian").unwrap(),
        frame,
        parity,
        reference,
    )
    .unwrap()
}
fn scalar(source: SampleSource<QtyAny>) -> ScalarPayload {
    ScalarPayload::new(
        meta(
            WORLD,
            OrientationParity::Even,
            ReferenceSemantics::Continuous,
        ),
        source,
    )
    .unwrap()
}
fn sample(payload: &ScalarPayload, time: f64) -> f64 {
    let value = payload
        .scalar_at(&FrameTree::new(), time, &mut || ControlFlow::Continue(()))
        .unwrap();
    assert_eq!(value.contract, QuantityContract::Dimensions(VOLTAGE));
    value.value
}

#[test]
fn scalar_and_vector_tables_share_knots_interpolation_and_outside_policies() {
    for interpolation in [TableInterpolation::StepLeft, TableInterpolation::Linear] {
        for outside in [
            OutsideDomainPolicy::Refuse,
            OutsideDomainPolicy::Clamp,
            OutsideDomainPolicy::Periodic,
        ] {
            let times = [1.0, 2.0, 3.0].map(|v| QtyAny::new(v, TIME)).to_vec();
            let values = [0.0, 4.0, -2.0].map(|v| QtyAny::new(v, VOLTAGE)).to_vec();
            let scalar = scalar(
                SampleSource::table(times.clone(), values.clone(), interpolation, outside).unwrap(),
            );
            let vector = VectorPayload::new(
                meta(
                    WORLD,
                    OrientationParity::Even,
                    ReferenceSemantics::Continuous,
                ),
                SampleSource::table(
                    times,
                    values.into_iter().map(|v| vec![v; 3]).collect(),
                    interpolation,
                    outside,
                )
                .unwrap(),
            )
            .unwrap();
            for time in [-1.5, 1.0, 1.5, 2.0, 2.75, 3.0, 3.5, 5.0] {
                let a =
                    scalar.scalar_at(&FrameTree::new(), time, &mut || ControlFlow::Continue(()));
                let b = vector
                    .world_vector_at(&FrameTree::new(), time, &mut || ControlFlow::Continue(()));
                match (a, b) {
                    (Ok(a), Ok(b)) => assert_eq!(a.value.to_bits(), b.values[0].to_bits()),
                    (Err(a), Err(b)) => assert_eq!(a, b),
                    other => panic!("scalar/vector policy mismatch: {other:?}"),
                }
            }
            assert_eq!(sample(&scalar, 2.0), 4.0);
            assert_eq!(sample(&scalar, 3.0), -2.0);
            assert_eq!(
                sample(&scalar, 1.5),
                if interpolation == TableInterpolation::Linear {
                    2.0
                } else {
                    0.0
                }
            );
        }
    }
    let huge = scalar(
        SampleSource::table(
            vec![QtyAny::new(0.0, TIME), QtyAny::new(2.0, TIME)],
            vec![QtyAny::new(-1e308, VOLTAGE), QtyAny::new(1e308, VOLTAGE)],
            TableInterpolation::Linear,
            OutsideDomainPolicy::Refuse,
        )
        .unwrap(),
    );
    assert_eq!(sample(&huge, 1.0), 0.0);
}

#[test]
fn scalar_frames_do_not_add_voltage_but_missing_frames_and_unowned_semantics_refuse() {
    let mut frames = FrameTree::new();
    frames.add(Frame {
        id: FrameId(1),
        name: "rotated".into(),
        parent: WORLD,
        motion: FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 1.0),
            translation: Vec3::new(1e12, -1e12, 1e12),
        },
    });
    let payload = ScalarPayload::new(
        meta(
            FrameId(1),
            OrientationParity::Even,
            ReferenceSemantics::Continuous,
        ),
        SampleSource::fixed(QtyAny::new(-3.0, VOLTAGE)),
    )
    .unwrap();
    assert_eq!(
        payload
            .scalar_at(&frames, 20.0, &mut || ControlFlow::Continue(()))
            .unwrap()
            .value,
        -3.0
    );
    assert!(matches!(
        payload.scalar_at(&FrameTree::new(), 20.0, &mut || ControlFlow::Continue(())),
        Err(ScalarEvaluationError::Frame(_))
    ));
    for (parity, reference) in [
        (OrientationParity::Odd, ReferenceSemantics::Continuous),
        (
            OrientationParity::Even,
            ReferenceSemantics::ResetAtEvent(PayloadId::new("switch").unwrap()),
        ),
    ] {
        let bad = ScalarPayload::new(
            meta(WORLD, parity, reference),
            SampleSource::fixed(QtyAny::new(1.0, VOLTAGE)),
        )
        .unwrap();
        assert!(matches!(
            bad.scalar_at(&frames, 0.0, &mut || ControlFlow::Continue(())),
            Err(ScalarEvaluationError::Unsupported(_))
        ));
    }
    let distribution = scalar(
        SampleSource::distribution(
            DistributionFamily::Empirical,
            vec![QtyAny::new(1.0, VOLTAGE)],
        )
        .unwrap(),
    );
    assert!(matches!(
        distribution.scalar_at(&frames, 0.0, &mut || ControlFlow::Continue(())),
        Err(ScalarEvaluationError::Unsupported(_))
    ));
    for stop in 1..=3 {
        let mut calls = 0;
        assert_eq!(
            payload.scalar_at(&frames, 0.0, &mut || {
                calls += 1;
                if calls == stop {
                    ControlFlow::Break(())
                } else {
                    ControlFlow::Continue(())
                }
            }),
            Err(ScalarEvaluationError::Cancelled)
        );
    }
    assert_eq!(
        payload.scalar_at(&frames, f64::NAN, &mut || ControlFlow::Continue(())),
        Err(ScalarEvaluationError::InvalidTime)
    );
}

#[test]
fn oriented_scalar_preserves_parity_without_weakening_ordinary_scalar_admission() {
    let dims = Dims([0, 1, -2, 0, -1, 0]);
    let mut frames = FrameTree::new();
    frames.add(Frame {
        id: FrameId(1),
        name: "proper-rotation".into(),
        parent: WORLD,
        motion: FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 1.0),
            translation: Vec3::new(10.0, 20.0, 30.0),
        },
    });
    let payload = ScalarPayload::new(
        PayloadMeta::new(
            QuantityContract::Dimensions(dims),
            PayloadId::new("basis/scalar").unwrap(),
            FrameId(1),
            OrientationParity::Odd,
            ReferenceSemantics::Continuous,
        )
        .unwrap(),
        SampleSource::table(
            [0.0, 1.0].map(|v| QtyAny::new(v, TIME)).to_vec(),
            [-2.0, 4.0].map(|v| QtyAny::new(v, dims)).to_vec(),
            TableInterpolation::Linear,
            OutsideDomainPolicy::Refuse,
        )
        .unwrap(),
    )
    .unwrap();
    let sampled = payload
        .oriented_scalar_at(&frames, 0.5, &mut || ControlFlow::Continue(()))
        .unwrap();
    assert_eq!(sampled.value, 1.0);
    assert_eq!(sampled.contract, QuantityContract::Dimensions(dims));
    assert_eq!(sampled.orientation, OrientationParity::Odd);
    assert!(matches!(
        payload.scalar_at(&frames, 0.5, &mut || ControlFlow::Continue(())),
        Err(ScalarEvaluationError::Unsupported(_))
    ));
    assert!(matches!(
        payload.oriented_scalar_at(&FrameTree::new(), 0.5, &mut || ControlFlow::Continue(())),
        Err(ScalarEvaluationError::Frame(_))
    ));
    assert_eq!(
        payload.oriented_scalar_at(&frames, 0.5, &mut || ControlFlow::Break(())),
        Err(ScalarEvaluationError::Cancelled)
    );
}


fn pressure_row(contract: QuantityContract, outside: OutsideDomainPolicy) -> fs_scenario::BoundaryCondition {
    use fs_scenario::{BcKind, BcValue, BoundaryCondition, Physics};
    use fs_scenario::payload::Payload;
    let dims = Dims([-1, 1, -2, 0, 0, 0]);
    BoundaryCondition {
        region: "reservoir".into(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::PressureOutlet,
        frame: WORLD.0,
        compatibility: None,
        value: Some(BcValue::Typed(Payload::Scalar(ScalarPayload::new(
            PayloadMeta::new(contract, PayloadId::new("basis/scalar").unwrap(), WORLD,
                OrientationParity::Even, ReferenceSemantics::Continuous).unwrap(),
            SampleSource::table(
                [0.0, 1.0, 2.0].map(|t| QtyAny::new(t, TIME)).to_vec(),
                [40.0, 80.0, 20.0].map(|p| QtyAny::new(p, dims)).to_vec(),
                TableInterpolation::Linear, outside,
            ).unwrap(),
        ).unwrap()))),
    }
}

#[test]
fn g3_pressure_history_round_trips_and_keeps_its_native_scalar_time_policy() {
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    use fs_scenario::{Environment, Scenario};
    let contract = QuantityContract::Semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Instantaneous));
    let row = pressure_row(contract, OutsideDomainPolicy::Periodic);
    let mut scenario = Scenario::new("typed-reservoir", 9, Environment::earth_lab());
    scenario.base_bcs.push(row);
    assert!(scenario.validate().is_empty());
    let ir = fs_scenario::ir::write_ir(&scenario);
    let parsed = fs_scenario::ir::parse_ir(&ir).unwrap();
    let decoded = &parsed.scenario().base_bcs[0];
    assert_eq!(decoded, &scenario.base_bcs[0]);
    for (t, expected) in [(0.5, 60.0), (2.0, 20.0), (2.5, 60.0), (-1.5, 60.0)] {
        let sample = decoded.pressure_at(&scenario.frames, t, &mut || ControlFlow::Continue(())).unwrap();
        assert_eq!(sample.value, expected);
        assert_eq!(sample.contract, contract);
    }
    let refused = pressure_row(contract, OutsideDomainPolicy::Refuse);
    assert!(refused.pressure_at(&scenario.frames, 2.1, &mut || ControlFlow::Continue(())).is_err());
}

#[test]
fn g0_pressure_admission_and_execution_do_not_confuse_stress_or_amplitudes_with_pressure() {
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    for (kind, form) in [
        (QuantityKind::Stress, ValueForm::Static),
        (QuantityKind::AcousticPressure, ValueForm::Instantaneous),
        (QuantityKind::Pressure, ValueForm::Peak),
        (QuantityKind::Pressure, ValueForm::Rms),
    ] {
        let row = pressure_row(QuantityContract::Semantic(SemanticType::new(kind, form)), OutsideDomainPolicy::Clamp);
        let mut findings = Vec::new();
        row.check(&mut findings);
        assert!(findings.iter().any(|v| v.code == "bc-pressure-semantics"));
        assert!(matches!(row.pressure_at(&FrameTree::new(), 0.5, &mut || ControlFlow::Continue(())),
            Err(ScalarEvaluationError::Unsupported(_))));
    }
    let mut row = pressure_row(QuantityContract::Dimensions(Dims([-1, 1, -2, 0, 0, 0])), OutsideDomainPolicy::Clamp);
    assert_eq!(row.pressure_at(&FrameTree::new(), 0.5, &mut || ControlFlow::Break(())), Err(ScalarEvaluationError::Cancelled));
    row.frame = 7;
    let mut findings = Vec::new();
    row.check(&mut findings);
    assert!(findings.iter().any(|v| v.code == "bc-payload-frame"));
    assert!(matches!(row.pressure_at(&FrameTree::new(), 0.5, &mut || ControlFlow::Continue(())),
        Err(ScalarEvaluationError::Unsupported(_))));
}
