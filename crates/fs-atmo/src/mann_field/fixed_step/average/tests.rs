use super::*;
use crate::mann::MANN_TARGET_V1;

fn field(seed: u64) -> FixedStepMannField {
    let model = MannSpectralModel::build(
        17,
        MANN_TARGET_V1,
        MannQuadrature {
            n_modes: 8,
            min_kl: 0.1,
            max_kl: 12.0,
        },
        [4.0, -0.3, 0.2],
    )
    .unwrap();
    let snapshot = HomogeneousMannField::at_anchor(model, seed, vec![0.04; 8]).unwrap();
    FixedStepMannField::from_snapshot(snapshot, seed ^ 0x9876, 0.0079, 88 * 8).unwrap()
}

// Independent time-domain Simpson integration of the OU stochastic
// convolution. No production F, divided differences or variance identity.
fn convolution_oracle(l: f64, theta: f64) -> (Complex, Complex, f64) {
    const N: usize = 16384;
    let mut mean = [0.0; 2];
    let mut cross = [0.0; 2];
    let mut variance = 0.0;
    let norm = l * l + theta * theta;
    for i in 0..=N {
        let r = i as f64 / N as f64;
        let weight = if i == 0 || i == N {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let e = (-l * r).exp();
        let phase = theta * r;
        let real = -(-l * r).exp_m1() - 2.0 * (phase / 2.0).sin().powi(2);
        let imaginary = phase.sin();
        let kernel = [
            (l * real + theta * imaginary) / norm,
            (l * imaginary - theta * real) / norm,
        ];
        mean[0] += weight * e * (theta * (1.0 - r)).cos();
        mean[1] += weight * e * (theta * (1.0 - r)).sin();
        cross[0] += weight * 2.0 * l * e * kernel[0];
        cross[1] += weight * 2.0 * l * e * kernel[1];
        variance += weight * 2.0 * l * (kernel[0] * kernel[0] + kernel[1] * kernel[1]);
    }
    let denominator = 3.0 * N as f64;
    (
        mean.map(|x| x / denominator),
        cross.map(|x| x / denominator),
        variance / denominator,
    )
}

#[test]
fn coefficients_match_the_joint_continuous_ou_integral_including_small_and_rotating_limits() {
    for (l, theta) in [
        (1e-16, 0.0),
        (1e-12, 1e-12),
        (1e-12, 0.49),
        (0.001, core::f64::consts::TAU),
        (0.2, 0.7),
        (0.5, 0.5),
        (0.500001, 0.500001),
        (2.0, -3.0),
        (20.0, 25.0),
        (50.0, 0.0),
        (1e-12, 100.0),
    ] {
        let sigma = det::sqrt(-det::expm1(-2.0 * l));
        let c = Coefficients::prepare(l, theta, sigma).unwrap();
        let (mean, cross, noise_variance) = convolution_oracle(l, theta);
        for i in 0..2 {
            assert!(
                (c.mean[i] - mean[i]).abs() < 3e-10,
                "mean: l={l} theta={theta}"
            );
            assert!(
                (c.endpoint[i] * sigma - cross[i]).abs() < 3e-10 * l.min(1.0),
                "endpoint: l={l} theta={theta} actual={:?} expected={cross:?}",
                c.endpoint
            );
        }
        let expected_bridge =
            noise_variance - (cross[0] * cross[0] + cross[1] * cross[1]) / (sigma * sigma);
        assert!(
            (c.bridge * c.bridge - expected_bridge).abs() < 1e-8 * expected_bridge,
            "bridge: l={l} theta={theta}, actual={} expected={expected_bridge}",
            c.bridge * c.bridge
        );
    }
    let l = 1e-16;
    let c = Coefficients::prepare(l, 0.0, det::sqrt(-det::expm1(-2.0 * l))).unwrap();
    assert!((c.endpoint[0] * c.endpoint[0] / l - 0.5).abs() < 5e-16);
    assert!((c.bridge * c.bridge / l - 1.0 / 6.0).abs() < 2e-16);
    // Squared phase and bridge variance overflow/underflow here, but their
    // scaled coefficients are representable and must not be zeroed out.
    for (l, theta) in [(1e308, 1e308), (1e-16, 1e308), (1e308, -1e308)] {
        let c = Coefficients::prepare(l, theta, det::sqrt(-det::expm1(-2.0 * l))).unwrap();
        assert!(c.bridge.is_finite() && c.bridge > 0.0);
        assert!(c.mean.iter().chain(&c.endpoint).all(|x| x.is_finite()));
    }
}

#[test]
fn endpoint_bits_and_future_means_survive_ordinary_advances_clones_and_checkpoint_resume() {
    let mut ordinary = field(39);
    let mut joint = ordinary.clone();
    let mut average = MannStepAverage::prepare(&joint, 104 * 8, || false).unwrap();
    assert_eq!(average.workspace_bytes(), 104 * 8);
    assert_eq!(average.interval(), None);
    assert!(average.sample([0.0; 3]).is_err());
    let buffers: Vec<_> = average
        .amplitudes
        .iter()
        .chain(&average.trial)
        .map(|v| v.as_ptr())
        .collect();
    for step in 1..=7 {
        ordinary.advance_to(step, 32, || false).unwrap();
        joint
            .advance_one_with_average(&mut average, 64, || false)
            .unwrap();
        assert_eq!(joint.checkpoint(), ordinary.checkpoint());
        assert_eq!(average.end_step(), Some(step));
        assert_eq!(
            average.interval().unwrap(),
            [
                joint.origin_time_s() + (step - 1) as f64 * joint.dt_s(),
                joint.time_s()
            ]
        );
        assert!(
            average
                .amplitudes
                .iter()
                .chain(&average.trial)
                .all(|v| buffers.contains(&v.as_ptr()))
        );
    }
    let mut resumed = FixedStepMannField::resume(joint.checkpoint(), 88 * 8).unwrap();
    let mut restored_average = MannStepAverage::prepare(&resumed, 104 * 8, || false).unwrap();
    joint
        .advance_one_with_average(&mut average, 64, || false)
        .unwrap();
    resumed
        .advance_one_with_average(&mut restored_average, 64, || false)
        .unwrap();
    assert_eq!(average, restored_average);
    assert_eq!(joint.checkpoint(), resumed.checkpoint());
    // An ordinary catch-up does not need or consume previous bridge draws.
    resumed.advance_to(12, 4 * 32, || false).unwrap();
    restored_average.validate_source(&resumed).unwrap();
    resumed
        .advance_one_with_average(&mut restored_average, 64, || false)
        .unwrap();
    joint.advance_to(12, 4 * 32, || false).unwrap();
    joint
        .advance_one_with_average(&mut average, 64, || false)
        .unwrap();
    assert_eq!(average, restored_average);
}

#[test]
fn every_spatial_observation_retains_pure_samples_solenoidality_and_analytic_gradients() {
    let mut source = field(17);
    let mut average = MannStepAverage::prepare(&source, 104 * 8, || false).unwrap();
    source
        .advance_one_with_average(&mut average, 64, || false)
        .unwrap();
    let wall = MannPlaneWall::new(-0.2).unwrap();
    let edges = [[0.2, 0.1, 0.0], [-0.1, 0.15, 0.0], [0.0, 0.0, 0.1]];
    let filter = MannBoxFilter::prepare(source.model(), edges, 64, || false).unwrap();
    let wall_filter =
        MannWallBoxFilter::prepare(source.model(), wall, edges, 256, || false).unwrap();
    let zero = MannBoxFilter::prepare(source.model(), [[0.0; 3]; 3], 0, || false).unwrap();
    let point = [0.2, -0.4, 0.7];
    assert_eq!(average.sample(point), average.sample_filtered(point, &zero));
    let before = average.clone();
    for variant in 0..4 {
        let observe = |p| {
            match variant {
                0 => average.sample(p),
                1 => average.sample_filtered(p, &filter),
                2 => average.sample_wall_blocked(p, wall),
                _ => average.sample_wall_blocked_filtered(p, &wall_filter),
            }
            .unwrap()
        };
        let sample = observe(point);
        assert!((sample.grad[0][0] + sample.grad[1][1] + sample.grad[2][2]).abs() < 2e-13);
        for axis in 0..3 {
            let mut plus = point;
            let mut minus = point;
            plus[axis] += 2e-6;
            minus[axis] -= 2e-6;
            let up = observe(plus).u;
            let down = observe(minus).u;
            for component in 0..3 {
                let numerical = (up[component] - down[component]) / 4e-6;
                assert!(
                    (numerical - sample.grad[component][axis]).abs()
                        < 2e-8 * (1.0 + sample.grad[component][axis].abs())
                );
            }
        }
        assert_eq!(sample, observe(point));
    }
    assert_eq!(
        average
            .sample_wall_blocked([0.2, -0.4, -0.2], wall)
            .unwrap()
            .u[2],
        0.0
    );
    assert_eq!(average, before);
}

#[test]
fn averaged_realizations_have_the_continuous_integral_variance_and_both_endpoint_correlations() {
    let quadrature = MannQuadrature {
        n_modes: 1,
        min_kl: 0.1,
        max_kl: 12.0,
    };
    let geometry = MannSpectralModel::build(17, MANN_TARGET_V1, quadrature, [0.0; 3]).unwrap();
    let k = geometry.modes[0].k;
    let dt = 0.04;
    let theta = 3.0;
    let velocity = k.map(|ki| ki * theta / (dot(k, k) * dt));
    let model = MannSpectralModel::build(17, MANN_TARGET_V1, quadrature, velocity).unwrap();
    let l = 2.0;
    let (m, _, variance_noise) = convolution_oracle(l, theta);
    let (end, _, _) = convolution_oracle(l, -theta);
    // The unconditional mean variance includes initial stationary variance.
    let target_variance = m[0] * m[0] + m[1] * m[1] + variance_noise;
    // Cov(B,A1) is integral exp((-l+i theta)r)dr = conj(F).
    let endpoint_target = multiply([theta.cos(), theta.sin()], end);
    let mut variance = 0.0;
    let mut initial_cross = [0.0; 2];
    let mut final_cross = [0.0; 2];
    const SEEDS: u64 = 4096;
    for seed in 0..SEEDS {
        let snapshot = HomogeneousMannField::at_anchor(model.clone(), seed, vec![dt / l]).unwrap();
        let mut source =
            FixedStepMannField::from_snapshot(snapshot, seed ^ 0x9876, dt, 88).unwrap();
        let initial: [f64; 4] = core::array::from_fn(|p| source.amplitudes[p][0]);
        let mut average = MannStepAverage::prepare(&source, 104, || false).unwrap();
        source
            .advance_one_with_average(&mut average, 8, || false)
            .unwrap();
        for polarization in 0..2 {
            let a = 2 * polarization;
            let b = a + 1;
            let x = average.amplitudes[a][0];
            let y = average.amplitudes[b][0];
            variance += x * x + y * y;
            initial_cross[0] += x * initial[a] + y * initial[b];
            initial_cross[1] += x * initial[b] - y * initial[a];
            final_cross[0] += x * source.amplitudes[a][0] + y * source.amplitudes[b][0];
            final_cross[1] += x * source.amplitudes[b][0] - y * source.amplitudes[a][0];
        }
    }
    let count = (4 * SEEDS) as f64;
    assert!((variance / count / target_variance - 1.0).abs() < 0.05);
    for i in 0..2 {
        assert!((initial_cross[i] / count - m[i]).abs() < 0.025);
        assert!((final_cross[i] / count - endpoint_target[i]).abs() < 0.025);
    }
}

#[test]
fn budgets_mismatch_and_every_cancellation_boundary_preserve_both_accepted_states() {
    let mut source = field(43);
    assert!(MannStepAverage::prepare(&source, 104 * 8 - 1, || false).is_err());
    assert!(MannStepAverage::required_workspace_bytes(usize::MAX).is_err());
    let mut average = MannStepAverage::prepare(&source, 104 * 8, || false).unwrap();
    source
        .advance_one_with_average(&mut average, 64, || false)
        .unwrap();
    let before_source = source.clone();
    let before_average = average.clone();
    assert!(
        source
            .advance_one_with_average(&mut average, 63, || false)
            .is_err()
    );
    let mut wrong = MannStepAverage::prepare(&field(44), 104 * 8, || false).unwrap();
    assert!(
        source
            .advance_one_with_average(&mut wrong, 64, || false)
            .is_err()
    );
    assert_eq!(source, before_source);
    assert_eq!(average, before_average);
    // Initial admission, one eight-mode tile, and final publication poll.
    for stop in 1..=3 {
        let mut calls = 0;
        let refusal = source
            .advance_one_with_average(&mut average, 64, || {
                calls += 1;
                calls == stop
            })
            .unwrap_err();
        assert_eq!(refusal.code, "mann-fixed-step-cancelled");
        assert_eq!(source, before_source);
        assert_eq!(average, before_average);
    }
    let mut reference = before_source.clone();
    let mut reference_average = before_average.clone();
    reference
        .advance_one_with_average(&mut reference_average, 64, || false)
        .unwrap();
    source
        .advance_one_with_average(&mut average, 64, || false)
        .unwrap();
    assert_eq!(source, reference);
    assert_eq!(average, reference_average);
}
