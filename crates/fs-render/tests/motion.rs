//! G0/G3 tests for deterministic shutter-time semantics.

use fs_render::motion::{
    MotionTimeError, NormalizedShutterTime, ShotTimeBounds, ShutterConvention, ShutterDistribution,
    ShutterInterval, TimedRay,
};

fn shot() -> ShotTimeBounds {
    ShotTimeBounds::try_new(0.0, 10.0).unwrap()
}

#[test]
fn centered_front_and_back_conventions_resolve_exactly() {
    let distribution = ShutterDistribution::UniformCounterV1;
    let centered =
        ShutterInterval::resolve(5.0, 2.0, ShutterConvention::Centered, distribution, shot())
            .unwrap();
    let front = ShutterInterval::resolve(
        5.0,
        2.0,
        ShutterConvention::FrontLoaded,
        distribution,
        shot(),
    )
    .unwrap();
    let back = ShutterInterval::resolve(
        5.0,
        2.0,
        ShutterConvention::BackLoaded,
        distribution,
        shot(),
    )
    .unwrap();
    assert_eq!((centered.open_s(), centered.close_s()), (4.0, 6.0));
    assert_eq!((front.open_s(), front.close_s()), (5.0, 7.0));
    assert_eq!((back.open_s(), back.close_s()), (3.0, 5.0));
}

#[test]
fn zero_width_shutter_reduces_every_ray_to_static_time() {
    let shutter = ShutterInterval::resolve(
        3.0,
        0.0,
        ShutterConvention::Centered,
        ShutterDistribution::UniformCounterV1,
        shot(),
    )
    .unwrap();
    for sample in 0..128 {
        let ray = TimedRay::from_sample([1_u8, 2, 3], shutter, 17, sample);
        assert_eq!(ray.spatial(), &[1, 2, 3]);
        assert_eq!(ray.absolute_time_s().to_bits(), 3.0_f64.to_bits());
    }
}

#[test]
fn explicit_endpoints_and_counter_sampling_are_deterministic_and_bounded() {
    let shutter = ShutterInterval::resolve(
        2.0,
        4.0,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::UniformCounterV1,
        shot(),
    )
    .unwrap();
    let open = TimedRay::at_normalized((), shutter, NormalizedShutterTime::try_new(0.0).unwrap());
    let close = TimedRay::at_normalized((), shutter, NormalizedShutterTime::try_new(1.0).unwrap());
    assert_eq!(open.absolute_time_s(), 2.0);
    assert_eq!(close.absolute_time_s(), 6.0);
    for sample in 0..128 {
        let first = TimedRay::from_sample((), shutter, 99, sample);
        let replay = TimedRay::from_sample((), shutter, 99, sample);
        assert_eq!(first, replay);
        assert!((2.0..6.0).contains(&first.absolute_time_s()));
    }
}

#[test]
fn malformed_and_out_of_shot_intervals_refuse() {
    assert_eq!(
        ShotTimeBounds::try_new(2.0, 1.0),
        Err(MotionTimeError::InvalidShotBounds)
    );
    assert_eq!(
        ShutterInterval::resolve(
            0.5,
            -0.1,
            ShutterConvention::Centered,
            ShutterDistribution::UniformCounterV1,
            shot(),
        ),
        Err(MotionTimeError::InvalidExposure)
    );
    assert_eq!(
        ShutterInterval::resolve(
            0.5,
            2.0,
            ShutterConvention::Centered,
            ShutterDistribution::UniformCounterV1,
            shot(),
        ),
        Err(MotionTimeError::ExposureOutsideShot)
    );
    for invalid in [f64::NAN, f64::INFINITY, -0.1, 1.1] {
        assert_eq!(
            NormalizedShutterTime::try_new(invalid),
            Err(MotionTimeError::InvalidNormalizedTime)
        );
    }
}
