use super::*;
use crate::mann::MANN_TARGET_V1;

const MODES: usize = 8;
const BYTES: usize = 88 * MODES;

fn field(dt_s: f64, tau_s: f64, origin_time_s: f64) -> FixedStepMannField {
    let model = MannSpectralModel::build(
        17,
        MANN_TARGET_V1,
        MannQuadrature {
            n_modes: MODES,
            min_kl: 0.1,
            max_kl: 12.0,
        },
        [4.0, -0.3, 0.2],
    )
    .unwrap();
    let snapshot = HomogeneousMannField::at_anchor(model, 73, vec![tau_s; MODES]).unwrap();
    let mut checkpoint = FixedStepMannField::from_snapshot(snapshot, 193, dt_s, BYTES)
        .unwrap()
        .checkpoint();
    checkpoint.origin_time_s = origin_time_s;
    FixedStepMannField::resume(checkpoint, BYTES).unwrap()
}

fn bridge(source: &FixedStepMannField) -> MannOuBridge {
    MannOuBridge::prepare(
        source,
        0xabcdef,
        MAX_MANN_OU_BRIDGE_DEPTH,
        BYTES,
        4 * MODES,
        || false,
    )
    .unwrap()
}

#[test]
fn endpoints_preserve_the_original_owner_and_all_synthesis_bits() {
    let source = field(0.0079, 0.04, -32.0);
    let before = source.checkpoint();
    let bridge = bridge(&source);
    let mut frame = bridge.frame(BYTES, || false).unwrap();
    let mut expected = source.clone();
    let filter = MannBoxFilter::prepare(
        source.model(),
        [[0.1, 0.04, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.3]],
        8 * MODES,
        || false,
    )
    .unwrap();
    for step in 0..=1 {
        expected.advance_to(step, 4 * MODES, || false).unwrap();
        let report = bridge
            .sample_frame(expected.time_s(), &mut frame, 0, || false)
            .unwrap();
        assert_eq!(report.normal_draws, 0);
        assert_eq!(report.tree_depth, 0);
        assert_eq!(frame.amplitudes, expected.amplitudes);
        assert_eq!(frame.time_s().to_bits(), expected.time_s().to_bits());
        for point in [[0.0; 3], [0.2, -0.4, 0.7]] {
            assert_eq!(frame.sample(point), expected.sample(point));
            assert_eq!(
                frame.sample_filtered(point, &filter),
                expected.sample_filtered(point, &filter)
            );
        }
    }
    assert_eq!(source.checkpoint(), before);
    assert_eq!(bridge.workspace_bytes(), BYTES);
    assert_eq!(frame.workspace_bytes(), BYTES);
    assert_eq!(frame.correlation_times_s(), source.correlation_times_s());
}

#[test]
fn arbitrary_stage_times_replay_across_query_orders_retries_and_intervals() {
    let mut source = field(0.1, 0.3, 0.0);
    let first = bridge(&source);
    let mut forward = first.frame(BYTES, || false).unwrap();
    let mut reverse = first.frame(BYTES, || false).unwrap();
    let times = [
        0.013,
        0.05,
        f64::from_bits(0.1f64.to_bits() - 1),
        f64::from_bits(0.05f64.to_bits() + 1),
    ];
    let expected: Vec<_> = times
        .iter()
        .map(|&time| {
            let count = first.required_normal_draws(time).unwrap();
            let report = first
                .sample_frame(time, &mut forward, count, || false)
                .unwrap();
            assert_eq!(report.normal_draws, count);
            assert_eq!(report.sampled_time_s.to_bits(), time.to_bits());
            assert_eq!(forward.time_s().to_bits(), time.to_bits());
            assert_eq!(forward.workspace_bytes(), BYTES);
            forward.amplitudes.clone()
        })
        .collect();
    for i in (0..times.len()).rev().chain(0..times.len()) {
        first
            .sample_frame(times[i], &mut reverse, usize::MAX, || false)
            .unwrap();
        assert_eq!(reverse.amplitudes, expected[i]);
    }
    source.advance_to(1, 4 * MODES, || false).unwrap();
    let next = bridge(&source);
    next.sample_frame(next.start_time_s(), &mut reverse, 0, || false)
        .unwrap();
    assert_eq!(reverse.amplitudes, source.amplitudes);
    next.sample_frame(0.123, &mut reverse, usize::MAX, || false)
        .unwrap();
    assert_eq!(reverse.workspace_bytes(), BYTES);
    assert_eq!(reverse.time_s(), 0.123);
}

#[test]
fn conditional_joint_statistics_match_an_independent_gaussian_schur_complement() {
    // Conditioning the stationary covariance exp(-|t-s|/tau) on both
    // endpoints is independent of the production midpoint conditioning.
    let tau = 0.7;
    let duration = 1.0;
    let mut bridge = bridge(&field(duration, tau, 0.0));
    bridge.left.iter_mut().for_each(|a| a.fill(-0.4));
    bridge.right.iter_mut().for_each(|a| a.fill(0.8));
    let mut frame = bridge.frame(BYTES, || false).unwrap();
    let times: [f64; 2] = [0.25, 0.75];
    let rho = (-duration / tau).exp();
    let endpoint_inverse = 1.0 / (1.0 - rho * rho);
    let rows = times.map(|time| [(-time / tau).exp(), (-(duration - time) / tau).exp()]);
    let weights = rows.map(|r| {
        [
            (r[0] - rho * r[1]) * endpoint_inverse,
            (r[1] - rho * r[0]) * endpoint_inverse,
        ]
    });
    let mean = weights.map(|w| -0.4 * w[0] + 0.8 * w[1]);
    let mut covariance = [[0.0; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            covariance[i][j] = (-(times[i] - times[j]).abs() / tau).exp()
                - weights[i][0] * rows[j][0]
                - weights[i][1] * rows[j][1];
        }
    }
    let trials = 2048;
    let count = (trials * 4 * MODES) as f64;
    let mut sampled_mean = [0.0; 2];
    let mut sampled_covariance = [[0.0; 2]; 2];
    for seed in 0..trials {
        bridge.bridge_seed = seed as u64;
        frame.bridge_seed = seed as u64;
        bridge
            .sample_frame(times[0], &mut frame, 8 * MODES, || false)
            .unwrap();
        let first = frame.amplitudes.clone();
        bridge
            .sample_frame(times[1], &mut frame, 8 * MODES, || false)
            .unwrap();
        for path in 0..4 {
            for mode in 0..MODES {
                let values = [first[path][mode], frame.amplitudes[path][mode]];
                for i in 0..2 {
                    sampled_mean[i] += values[i] / count;
                    for j in 0..2 {
                        sampled_covariance[i][j] +=
                            (values[i] - mean[i]) * (values[j] - mean[j]) / count;
                    }
                }
            }
        }
    }
    for i in 0..2 {
        assert!((sampled_mean[i] - mean[i]).abs() < 0.015 * covariance[i][i].sqrt());
        for j in 0..2 {
            let scale = (covariance[i][i] * covariance[j][j]).sqrt();
            assert!((sampled_covariance[i][j] - covariance[i][j]).abs() < 0.025 * scale);
        }
    }
}

#[test]
fn conditioning_stays_representable_at_slow_fast_and_rounded_clock_limits() {
    let slow = coefficients(1e-120, 3e-120, 1e120).unwrap();
    assert!((slow[0] - 0.75).abs() < 2e-15);
    assert!((slow[1] - 0.25).abs() < 2e-15);
    assert!((slow[2] / 1.5e-240f64.sqrt() - 1.0).abs() < 2e-15);
    assert_eq!(coefficients(1e200, 2e200, 1e-200).unwrap(), [0.0, 0.0, 1.0]);

    // The represented physical gap differs from nominal dt. All geometry
    // nodes must use the same affine OU clock, while phase retains real time.
    let bridge = bridge(&field(0.003, 0.02, 1e12));
    assert_ne!(bridge.end_time_s() - bridge.start_time_s(), bridge.dt_s());
    let time = bridge.start_time_s() + 0.001;
    let path = bridge.query_path(time).unwrap();
    let root = path.nodes[0];
    assert!(((root.left_duration_s + root.right_duration_s) / bridge.dt_s() - 1.0).abs() < 3e-16);
    let mut frame = bridge.frame(BYTES, || false).unwrap();
    bridge
        .sample_frame(time, &mut frame, usize::MAX, || false)
        .unwrap();
    assert_eq!(frame.time_s().to_bits(), time.to_bits());

    let negative = super::tests::bridge(&field(0.2, 0.5, -0.1));
    let mut frame = negative.frame(BYTES, || false).unwrap();
    for time in [-0.073, -0.0, 0.0, 0.013] {
        negative
            .sample_frame(time, &mut frame, usize::MAX, || false)
            .unwrap();
        assert!(frame.amplitudes.iter().flatten().all(|x| x.is_finite()));
        assert_eq!(frame.time_s().to_bits(), time.to_bits());
    }
    assert_eq!(negative.normal(-0.0, 0, 0), negative.normal(0.0, 0, 0));
    assert_ne!(negative.normal(-0.05, 0, 0), negative.normal(0.05, 0, 0));
}

#[test]
fn work_depth_identity_and_cancellation_refusals_are_atomic() {
    let source = field(0.1, 0.3, 0.0);
    let bridge = bridge(&source);
    let mut frame = bridge.frame(BYTES, || false).unwrap();
    let initial = frame.amplitudes.clone();
    let initial_time = frame.time_s();
    let time = 0.013;
    let count = bridge.required_normal_draws(time).unwrap();
    assert!(
        bridge
            .sample_frame(time, &mut frame, count - 1, || false)
            .is_err()
    );
    assert_eq!(frame.amplitudes, initial);
    assert_eq!(frame.time_s(), initial_time);
    assert!(MannOuBridge::prepare(&source, 0, 128, BYTES - 1, usize::MAX, || false).is_err());
    assert!(MannOuBridge::prepare(&source, 0, 128, BYTES, 4 * MODES - 1, || false).is_err());
    assert!(bridge.frame(BYTES - 1, || false).is_err());
    let shallow = MannOuBridge::prepare(&source, 0xabcdef, 1, BYTES, 4 * MODES, || false).unwrap();
    assert!(
        shallow
            .sample_frame(time, &mut frame, usize::MAX, || false)
            .is_err()
    );
    assert!(bridge.required_normal_draws(f64::from_bits(1)).is_err());
    assert!(bridge.required_normal_draws(f64::NAN).is_err());
    assert!(bridge.required_normal_draws(-1.0).is_err());
    let mut incompatible = bridge.frame(BYTES, || false).unwrap();
    incompatible.bridge_seed ^= 1;
    assert!(
        bridge
            .sample_frame(0.0, &mut incompatible, 0, || false)
            .is_err()
    );

    let mut polls = 0;
    bridge
        .sample_frame(time, &mut frame, count, || {
            polls += 1;
            false
        })
        .unwrap();
    let expected = frame.amplitudes.clone();
    for stop_at in [1, polls / 2, polls] {
        bridge.sample_frame(0.0, &mut frame, 0, || false).unwrap();
        let mut calls = 0;
        assert!(
            bridge
                .sample_frame(time, &mut frame, count, || {
                    calls += 1;
                    calls == stop_at
                })
                .is_err()
        );
        assert_eq!(frame.amplitudes, initial);
        assert_eq!(frame.time_s(), initial_time);
        bridge
            .sample_frame(time, &mut frame, count, || false)
            .unwrap();
        assert_eq!(frame.amplitudes, expected);
    }
}
