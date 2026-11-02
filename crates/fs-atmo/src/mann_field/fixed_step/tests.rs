use super::*;
use crate::mann::MANN_TARGET_V1;
use crate::ou::STATIONARY_ANCHOR_TICK;

fn model() -> MannSpectralModel {
    MannSpectralModel::build(
        17,
        MANN_TARGET_V1,
        MannQuadrature {
            n_modes: 8,
            min_kl: 0.1,
            max_kl: 12.0,
        },
        [4.0, -0.3, 0.2],
    )
    .unwrap()
}
fn snapshot(seed: u64) -> HomogeneousMannField {
    HomogeneousMannField::at_anchor(model(), seed, vec![0.04; 8]).unwrap()
}
fn field(seed: u64) -> FixedStepMannField {
    FixedStepMannField::from_snapshot(snapshot(seed), seed ^ 0x9876, 0.0079, 88 * 8).unwrap()
}

#[test]
fn g0_snapshot_samples_filters_and_covariances_retain_the_original_bits() {
    let mut source = snapshot(29);
    source.advance_to(STATIONARY_ANCHOR_TICK + 11).unwrap();
    let filter = MannBoxFilter::prepare(
        source.model(),
        [[0.3, 0.2, 0.0], [-0.1, 0.4, 0.0], [0.0, 0.0, 0.2]],
        64,
        || false,
    )
    .unwrap();
    let fixed = FixedStepMannField::from_snapshot(source.clone(), 123, 0.0013, 88 * 8).unwrap();
    assert_eq!(fixed.step(), 0);
    assert_eq!(fixed.elapsed_time_s(), 0.0);
    assert_eq!(fixed.time_s(), source.tick() as f64 / TICK_HZ);
    assert_eq!(fixed.workspace_bytes(), 88 * 8);
    for p in [[0.0; 3], [0.1, -0.7, 0.9], [1.2, 2.3, -0.4]] {
        assert_eq!(fixed.sample(p), source.sample(p));
        assert_eq!(
            fixed.sample_filtered(p, &filter),
            source.sample_filtered(p, &filter)
        );
        for lag in [0.0, 0.01, -0.03] {
            assert_eq!(fixed.covariance(p, lag), source.covariance(p, lag));
            assert_eq!(
                fixed.covariance_filtered(p, lag, &filter),
                source.covariance_filtered(p, lag, &filter)
            );
        }
    }
    assert!(Arc::ptr_eq(&fixed.definition, &fixed.clone().definition));
}

#[test]
fn g0_fixed_step_matches_the_exact_two_step_gaussian_transition() {
    let mut f = field(42);
    let initial = f.amplitudes.clone();
    let h = f.dt_s() / 0.04;
    let rho = (-h).exp();
    let sigma = (-(-2.0 * h).exp_m1()).sqrt();
    assert!((rho * rho + sigma * sigma - 1.0).abs() < 2e-16);
    f.advance_to(2, 4 * 8 * 2, || false).unwrap();
    for path in 0..4 {
        for mode in 0..8 {
            let key = StreamKey {
                seed: (42 ^ 0x9876) ^ PATH_KEYS[path],
                kernel: FIXED_STEP_MANN_KERNEL,
                tile: mode as u32,
            };
            let x1 = Stream::resume(StreamCheckpoint::current(key, 2))
                .unwrap()
                .next_normal();
            let x2 = Stream::resume(StreamCheckpoint::current(key, 4))
                .unwrap()
                .next_normal();
            let expected = (-2.0 * h).exp() * initial[path][mode] + sigma * (rho * x1 + x2);
            assert!((f.amplitudes[path][mode] - expected).abs() < 3e-15);
        }
    }
    assert_eq!(f.elapsed_time_s(), 2.0 * 0.0079);
    let p = [0.2, 0.4, -0.3];
    let frozen_phase = sample_modes(
        f.model(),
        f.origin_time_s(),
        core::array::from_fn(|i| f.amplitudes[i].as_slice()),
        core::array::from_fn(|i| p[i] - f.advection_mps()[i] * f.elapsed_time_s()),
        None,
    )
    .unwrap();
    let actual = f.sample(p).unwrap();
    for (a, b) in actual.u.iter().zip(frozen_phase.u) {
        assert!((a - b).abs() < 2e-14);
    }
    assert!((actual.grad[0][0] + actual.grad[1][1] + actual.grad[2][2]).abs() < 2e-14);
}

#[test]
fn g1_cfd_rate_paths_have_the_declared_stationary_tensor_and_lag_covariance() {
    let m = model();
    let filter = MannBoxFilter::prepare(
        &m,
        [[0.3, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.1]],
        64,
        || false,
    )
    .unwrap();
    let point = [0.2, -0.4, 0.7];
    let separation = [0.4, -0.1, 0.3];
    let lag = 4.0 * 0.0079;
    let prototype = field(1);
    let zero = prototype
        .covariance_filtered([0.0; 3], 0.0, &filter)
        .unwrap();
    // Independent tensor sum uses the original spectral tensor, without the
    // sampled polarizations, production covariance helper or temporal cache.
    let gains = filter.gains(&m).unwrap().unwrap();
    let mut target = [[0.0; 3]; 3];
    for (mode, gain) in m.modes.iter().zip(gains) {
        let tensor = mann_tensor(&m.parameters, mode.k).unwrap();
        let phase: f64 = (0..3)
            .map(|a| mode.k[a] * (separation[a] - m.advection_mps[a] * lag))
            .sum();
        let weight = mode.weight * gain * gain * phase.cos() * (-lag / 0.04).exp();
        for i in 0..3 {
            for j in 0..3 {
                target[i][j] += weight * tensor[i][j];
            }
        }
    }
    let represented = prototype
        .covariance_filtered(separation, lag, &filter)
        .unwrap();
    for i in 0..3 {
        for j in 0..3 {
            assert!((represented[i][j] - target[i][j]).abs() < 2e-12);
        }
    }
    let mut first = [[0.0; 3]; 3];
    let mut last = [[0.0; 3]; 3];
    let mut cross = [[0.0; 3]; 3];
    let count = 4096;
    for seed in 0..count {
        let source = HomogeneousMannField::at_anchor(m.clone(), seed, vec![0.04; 8]).unwrap();
        let mut f =
            FixedStepMannField::from_snapshot(source, seed ^ 0x9876, 0.0079, 88 * 8).unwrap();
        let a = f.sample_filtered(point, &filter).unwrap().u;
        f.advance_to(4, 4 * 8 * 4, || false).unwrap();
        let b = f
            .sample_filtered(core::array::from_fn(|i| point[i] + separation[i]), &filter)
            .unwrap()
            .u;
        for i in 0..3 {
            for j in 0..3 {
                first[i][j] += a[i] * a[j] / count as f64;
                last[i][j] += b[i] * b[j] / count as f64;
                cross[i][j] += a[i] * b[j] / count as f64;
            }
        }
    }
    for i in 0..3 {
        for j in 0..3 {
            let scale = (zero[i][i] * zero[j][j]).sqrt();
            assert!((first[i][j] - zero[i][j]).abs() < 0.08 * scale);
            assert!((last[i][j] - zero[i][j]).abs() < 0.08 * scale);
            assert!((cross[i][j] - target[i][j]).abs() < 0.08 * scale);
        }
    }
}

#[test]
fn g3_checkpoint_and_storage_reuse_preserve_the_complete_future_path() {
    let mut original = field(17);
    original.advance_to(3, 4 * 8 * 3, || false).unwrap();
    let mut restored = FixedStepMannField::resume(original.checkpoint(), 88 * 8).unwrap();
    let mut staged = field(17);
    staged.copy_temporal_state_from(&original).unwrap();
    assert_eq!(restored, original);
    assert_eq!(staged, original);
    for target in [5, 11] {
        original.advance_to(target, 4 * 8 * 6, || false).unwrap();
        restored.advance_to(target, 4 * 8 * 6, || false).unwrap();
        staged.advance_to(target, 4 * 8 * 6, || false).unwrap();
        assert_eq!(restored.checkpoint(), original.checkpoint());
        assert_eq!(
            staged.sample([0.1, 0.2, 0.3]),
            original.sample([0.1, 0.2, 0.3])
        );
        assert_eq!(original.workspace_bytes(), 88 * 8);
    }
    for which in 0..5 {
        let mut cp = original.checkpoint();
        match which {
            0 => cp.dt_s *= 2.0,
            1 => cp.origin_time_s += 1.0,
            2 => cp.evolution_seed ^= 1,
            3 => cp.origin_evolution_seed ^= 1,
            _ => cp.correlation_times_s[0] *= 2.0,
        }
        let mismatched = FixedStepMannField::resume(cp, 88 * 8).unwrap();
        let before = original.checkpoint();
        assert!(original.copy_temporal_state_from(&mismatched).is_err());
        assert_eq!(original.checkpoint(), before);
    }
}

#[test]
fn g4_work_clock_storage_and_late_cancel_refuse_before_publication() {
    let initial = field(5);
    let mut expected = initial.clone();
    let mut polls = 0;
    expected
        .advance_to(3, 4 * 8 * 3, || {
            polls += 1;
            false
        })
        .unwrap();
    for at in [1, polls / 2, polls] {
        let mut cancelled = initial.clone();
        let mut calls = 0;
        assert_eq!(
            cancelled
                .advance_to(3, 4 * 8 * 3, || {
                    calls += 1;
                    calls == at
                })
                .unwrap_err()
                .code,
            "mann-fixed-step-cancelled"
        );
        assert_eq!(cancelled, initial);
        assert_eq!(cancelled.checkpoint(), initial.checkpoint());
        cancelled.advance_to(3, 4 * 8 * 3, || false).unwrap();
        assert_eq!(cancelled, expected);
    }
    let before = expected.checkpoint();
    for (step, budget) in [(2, usize::MAX), (4, 31), (MAX_FIXED_STEP + 1, usize::MAX)] {
        assert!(expected.advance_to(step, budget, || false).is_err());
        assert_eq!(expected.checkpoint(), before);
    }
    expected.advance_to(3, 0, || false).unwrap();
    assert_eq!(expected.checkpoint(), before);
    for dt in [0.0, -1.0, f64::INFINITY, f64::NAN, 1e-300] {
        assert!(FixedStepMannField::from_snapshot(snapshot(7), 8, dt, 88 * 8).is_err());
    }
    assert!(FixedStepMannField::from_snapshot(snapshot(7), 8, 0.001, 88 * 8 - 1).is_err());
    let mut cp = initial.checkpoint();
    cp.amplitudes[0].reserve_exact(64);
    assert!(FixedStepMannField::resume(cp, 88 * 8).is_err());
    let mut tiny = initial.checkpoint();
    tiny.origin_time_s = 0.0;
    tiny.dt_s = 1e-20;
    assert_eq!(
        FixedStepMannField::resume(tiny, usize::MAX)
            .unwrap_err()
            .code,
        "mann-fixed-step-clock"
    );
    for which in 0..4 {
        let mut cp = initial.checkpoint();
        match which {
            0 => cp.version += 1,
            1 => cp.amplitudes[0][0] = f64::NAN,
            2 => {
                cp.amplitudes[0].pop();
            }
            _ => cp.stream_semantics_version += 1,
        }
        assert!(FixedStepMannField::resume(cp, usize::MAX).is_err());
    }
}
