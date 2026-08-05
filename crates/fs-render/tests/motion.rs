//! G0/G3 tests for deterministic shutter-time semantics.

use fs_render::halton;
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
fn zero_width_shutter_canonicalizes_signed_zero_frame_time() {
    let negative = ShutterInterval::resolve(
        -0.0,
        0.0,
        ShutterConvention::Centered,
        ShutterDistribution::UniformCounterV1,
        ShotTimeBounds::try_new(-1.0, 1.0).unwrap(),
    )
    .unwrap();
    let positive = ShutterInterval::resolve(
        0.0,
        0.0,
        ShutterConvention::Centered,
        ShutterDistribution::UniformCounterV1,
        ShotTimeBounds::try_new(-1.0, 1.0).unwrap(),
    )
    .unwrap();

    assert_eq!(negative, positive);
    assert_eq!(negative.open_s().to_bits(), 0.0_f64.to_bits());
    assert_eq!(negative.close_s().to_bits(), 0.0_f64.to_bits());
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
fn uniform_counter_v1_compatibility_stream_has_locked_output_bits() {
    let shutter = ShutterInterval::resolve(
        2.0,
        4.0,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::UniformCounterV1,
        shot(),
    )
    .unwrap();
    for ((pixel, sample), expected_bits) in [
        ((0, 0), 0x3fe8_175e_9a90_db40),
        ((1, 0), 0x3fec_c6ce_35a6_021a),
        ((99, 17), 0x3fda_007a_4dd6_49f0),
        ((u64::MAX, u64::MAX), 0x3fc6_7881_f485_56bc),
    ] {
        assert_eq!(
            shutter.sample(pixel, sample).value().to_bits(),
            expected_bits
        );
        assert_eq!(
            shutter
                .sample_for_stream(0, pixel, sample)
                .value()
                .to_bits(),
            expected_bits
        );
    }
}

#[test]
fn explicit_endpoints_survive_extreme_dynamic_range_exactly() {
    let close_s = -1.0e-300;
    let shutter = ShutterInterval::resolve(
        close_s,
        1.0e10,
        ShutterConvention::BackLoaded,
        ShutterDistribution::UniformCounterV1,
        ShotTimeBounds::try_new(-1.0e10, close_s).unwrap(),
    )
    .unwrap();

    assert_eq!(
        shutter
            .time_at(NormalizedShutterTime::try_new(0.0).unwrap())
            .to_bits(),
        shutter.open_s().to_bits()
    );
    assert_eq!(
        shutter
            .time_at(NormalizedShutterTime::try_new(1.0).unwrap())
            .to_bits(),
        close_s.to_bits()
    );
}

#[test]
fn full_shot_and_stratified_distributions_are_admitted_and_cover_every_stratum() {
    let full_shot = ShutterInterval::resolve(
        0.0,
        10.0,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::StratifiedCounterV1 { strata: 8 },
        shot(),
    )
    .unwrap();
    assert_eq!((full_shot.open_s(), full_shot.close_s()), (0.0, 10.0));

    let mut visited = [false; 8];
    for sample_identity in 0..8 {
        let coordinate = full_shot.sample(41, sample_identity).value();
        let stratum = (coordinate * 8.0).floor() as usize;
        assert!(stratum < visited.len(), "sample escaped admitted shutter");
        visited[stratum] = true;
    }
    assert!(visited.into_iter().all(core::convert::identity));

    // The stratum permutation repeats cyclically, so an arbitrarily aligned
    // progressive window also contains each stratum exactly once.
    let mut crossing_cycle = [false; 8];
    for sample_identity in 5..13 {
        let coordinate = full_shot.sample(41, sample_identity).value();
        let stratum = (coordinate * 8.0).floor() as usize;
        assert!(!crossing_cycle[stratum], "stratum {stratum} repeated");
        crossing_cycle[stratum] = true;
    }
    assert!(crossing_cycle.into_iter().all(core::convert::identity));
}

#[test]
fn partial_stratified_batches_are_not_biased_to_shutter_open() {
    let shutter = ShutterInterval::resolve(
        0.0,
        10.0,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::StratifiedCounterV1 { strata: 4_096 },
        shot(),
    )
    .unwrap();
    for stream_identity in 0..16 {
        let coordinates: Vec<_> = (0..64)
            .map(|sample_identity| {
                shutter
                    .sample_for_stream(stream_identity, 73 + stream_identity, sample_identity)
                    .value()
            })
            .collect();
        let minimum = coordinates.iter().copied().fold(f64::INFINITY, f64::min);
        let maximum = coordinates
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let mean = coordinates.iter().sum::<f64>() / coordinates.len() as f64;
        assert!(
            minimum < 0.2,
            "stream {stream_identity} missed shutter opening: {minimum}"
        );
        assert!(
            maximum > 0.8,
            "stream {stream_identity} missed shutter closing: {maximum}"
        );
        assert!(
            (0.3..0.7).contains(&mean),
            "stream {stream_identity} remained temporally biased: mean={mean}"
        );
    }
}

#[test]
fn logical_sample_ids_replay_identically_across_progressive_partitions() {
    let shutter = ShutterInterval::resolve(
        5.0,
        2.0,
        ShutterConvention::Centered,
        ShutterDistribution::StratifiedCounterV1 { strata: 16 },
        shot(),
    )
    .unwrap();
    let full: Vec<_> = (0..256).map(|sample| shutter.sample(7, sample)).collect();
    let mut partitioned = Vec::new();
    for range in [0..31, 31..113, 113..256] {
        partitioned.extend(range.map(|sample| shutter.sample(7, sample)));
    }
    assert_eq!(full, partitioned);
}

#[test]
fn explicit_stream_identity_replays_and_reseeds_temporal_samples() {
    let shutter = ShutterInterval::resolve(
        5.0,
        2.0,
        ShutterConvention::Centered,
        ShutterDistribution::StratifiedCounterV1 { strata: 64 },
        shot(),
    )
    .unwrap();
    let first: Vec<_> = (0..32)
        .map(|sample| shutter.sample_for_stream(0xaaaa, 7, sample))
        .collect();
    let replay: Vec<_> = (0..32)
        .map(|sample| shutter.sample_for_stream(0xaaaa, 7, sample))
        .collect();
    let reseeded: Vec<_> = (0..32)
        .map(|sample| shutter.sample_for_stream(0xbbbb, 7, sample))
        .collect();

    assert_eq!(first, replay);
    assert_ne!(first, reseeded);
}

#[test]
fn adding_shutter_time_does_not_perturb_existing_spatial_sample_dimensions() {
    let shutter = ShutterInterval::resolve(
        5.0,
        1.0,
        ShutterConvention::Centered,
        ShutterDistribution::UniformCounterV1,
        shot(),
    )
    .unwrap();
    for sample_identity in 1..128 {
        let spatial_dimensions = [
            halton(0, sample_identity),
            halton(1, sample_identity),
            halton(2, sample_identity),
            halton(3, sample_identity),
        ];
        let timed = TimedRay::from_sample(spatial_dimensions, shutter, 91, sample_identity);
        assert_eq!(timed.spatial(), &spatial_dimensions);
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
            1.0e20,
            1.0,
            ShutterConvention::FrontLoaded,
            ShutterDistribution::UniformCounterV1,
            ShotTimeBounds::try_new(0.0, 2.0e20).unwrap(),
        ),
        Err(MotionTimeError::CollapsedExposure)
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
    assert_eq!(
        ShutterInterval::resolve(
            0.5,
            0.1,
            ShutterConvention::Centered,
            ShutterDistribution::StratifiedCounterV1 { strata: 0 },
            shot(),
        ),
        Err(MotionTimeError::InvalidDistribution)
    );
    for invalid in [f64::NAN, f64::INFINITY, -0.1, 1.1] {
        assert_eq!(
            NormalizedShutterTime::try_new(invalid),
            Err(MotionTimeError::InvalidNormalizedTime)
        );
    }
}
