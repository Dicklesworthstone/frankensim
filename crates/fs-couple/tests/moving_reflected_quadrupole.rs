//! G1/G3 original stress tensors at moving microphones and ideal image paths.
use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    IdealReflectionObserver, MovingReceiverBudget, MovingReceiverConfig, MovingReflectionObserver,
    PreparedKinematicTrajectory, PreparedMovingReceiverSource,
    PreparedReflectedMovingReceiverSource, PreparedReflectedSourceTrajectory,
    PreparedSourceTrajectory, ReflectedTrajectoryBudget, TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{MovingCompactLimits, MovingCompactObserver, SourceKinematics};
use fs_aeroac::curle3d::quadrupole::transient::TransientQuadrupoleSource;
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};
use fs_couple::pcm_wav::observation::moving::MovingPressureConfig;
use fs_couple::pcm_wav::observation::quadrupole::QuadrupolePressureBudget;
use fs_couple::pcm_wav::observation::quadrupole::moving::{
    MovingQuadrupoleRenderer, MovingQuadrupoleSourceWindow, QuadrupoleMotion,
};
use std::cell::{Cell, RefCell};

const C: f64 = 300.0;
const TENSORS: [[f64; 6]; 2] = [
    [0.7, -0.2, 0.4, 0.3, -0.1, 0.23],
    [-0.1, 0.8, 0.2, -0.25, 0.4, -0.15],
];
const ENTRIES: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
fn medium() -> Fluid {
    Fluid::try_new(1.2, C).unwrap()
}
type Fluid = fs_material::fluid_acoustics::FluidAcousticState;
fn bounds() -> TrajectoryBudget {
    TrajectoryBudget {
        maximum_segments: 16,
        maximum_payload_bytes: 131_072,
        maximum_bound_nodes: 16_384,
        maximum_subdivision_depth: 20,
        maximum_root_iterations: 96,
    }
}
fn image_bounds() -> ReflectedTrajectoryBudget {
    ReflectedTrajectoryBudget {
        maximum_payload_bytes: 1_000_000,
        maximum_bound_nodes: 16_384,
        maximum_subdivision_depth: 20,
    }
}
fn pair_bounds() -> MovingReceiverBudget {
    MovingReceiverBudget {
        maximum_segment_pairs: 512,
        maximum_bound_nodes: 32_768,
        maximum_subdivision_depth: 20,
        maximum_root_iterations: 96,
    }
}
fn pair_config() -> MovingReceiverConfig {
    MovingReceiverConfig {
        receiver_time_interval_s: [0.08, 0.18],
        minimum_distance_m: 0.001,
        maximum_time_residual_s: 1e-12,
    }
}
fn config() -> MovingPressureConfig {
    MovingPressureConfig {
        first_receiver_time_s: 0.10,
        sample_rate_hz: 1000,
        samples: 17,
        max_block: 17,
        receiver_bandwidth_hz: 20.0,
        maximum_time_residual_s: 1e-12,
    }
}
fn budget() -> QuadrupolePressureBudget {
    QuadrupolePressureBudget {
        maximum_sources: 8,
        maximum_output_samples: 4096,
        maximum_source_path_samples: 65_536,
        maximum_scratch_bytes: 1_000_000,
    }
}
fn knots(motion: impl Fn(f64) -> SourceKinematics) -> Vec<TimedSourceKinematics> {
    [0.0, 0.04, 0.11, 0.19, 0.3]
        .into_iter()
        .map(|source_time_s| TimedSourceKinematics {
            source_time_s,
            kinematics: motion(source_time_s),
        })
        .collect()
}
fn path(motion: impl Fn(f64) -> SourceKinematics) -> PreparedKinematicTrajectory {
    PreparedKinematicTrajectory::new(medium(), &knots(motion), 0.3, bounds(), &mut || false)
        .unwrap()
}
fn stationary(position_m: [f64; 3]) -> SourceKinematics {
    SourceKinematics {
        position_m,
        velocity_m_s: [0.0; 3],
        acceleration_m_s2: [0.0; 3],
    }
}
fn source_motion(index: usize, t: f64) -> SourceKinematics {
    if index == 0 {
        SourceKinematics {
            position_m: [
                2.0 + 12.0 * t + 15.0 * t * t + 6.0 * t * t * t,
                1.2 + 2.0 * t - 4.0 * t * t + 2.0 * t * t * t,
                0.5 + 0.3 * t * t - 0.2 * t * t * t,
            ],
            velocity_m_s: [
                12.0 + 30.0 * t + 18.0 * t * t,
                2.0 - 8.0 * t + 6.0 * t * t,
                0.6 * t - 0.6 * t * t,
            ],
            acceleration_m_s2: [30.0 + 36.0 * t, -8.0 + 12.0 * t, 0.6 - 1.2 * t],
        }
    } else {
        SourceKinematics {
            position_m: [
                1.1 + 5.0 * t + 4.0 * t * t + 3.0 * t * t * t,
                3.0 + 0.6 * t + 2.0 * t * t - 0.5 * t * t * t,
                1.0 + 0.2 * t + 0.1 * t * t * t,
            ],
            velocity_m_s: [
                5.0 + 8.0 * t + 9.0 * t * t,
                0.6 + 4.0 * t - 1.5 * t * t,
                0.2 + 0.3 * t * t,
            ],
            acceleration_m_s2: [8.0 + 18.0 * t, 4.0 - 3.0 * t, 0.6 * t],
        }
    }
}
fn receiver_motion(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [
            10.0 - 4.0 * t + 2.0 * t * t,
            6.0 + 3.0 * t - 2.0 * t * t,
            3.0 + 0.4 * t * t,
        ],
        velocity_m_s: [-4.0 + 4.0 * t, 3.0 - 4.0 * t, 0.8 * t],
        acceleration_m_s2: [4.0, -4.0, 0.8],
    }
}
fn stress(t: f64, tensor: [f64; 6], omega: f64, phase: f64) -> TransientQuadrupoleSource {
    let (s, c) = (omega * t + phase).sin_cos();
    TransientQuadrupoleSource {
        stress_moment_n_m: tensor.map(|q| q * c),
        stress_rate_n_m_s: tensor.map(|q| -q * omega * s),
        stress_acceleration_n_m_s2: tensor.map(|q| -q * omega * omega * c),
        ..Default::default()
    }
}
fn read(index: usize, t: f64) -> Result<TransientQuadrupoleSource, AeroacError> {
    Ok(stress(
        t,
        TENSORS[index],
        core::f64::consts::TAU * [7.0, 11.0][index],
        [0.2, -0.6][index],
    ))
}
fn window<M>(motion: M) -> MovingQuadrupoleSourceWindow<M> {
    MovingQuadrupoleSourceWindow {
        motion,
        source_time_interval_s: [0.0, 0.3],
        source_bandwidth_hz: 11.0,
    }
}
fn planes() -> [IdealAcousticPlane; 2] {
    [
        IdealAcousticPlane::new([0.0; 3], [1.0, 0.0, 0.0], IdealPlaneBoundary::Rigid).unwrap(),
        IdealAcousticPlane::new(
            [0.0; 3],
            [0.0, 1.0, 0.0],
            IdealPlaneBoundary::PressureRelease,
        )
        .unwrap(),
    ]
}
fn render<
    M: QuadrupoleMotion,
    F: Fn(usize, f64) -> Result<TransientQuadrupoleSource, AeroacError>,
>(
    observer: M::Observer,
    sources: &[MovingQuadrupoleSourceWindow<M>],
    reader: F,
) -> Vec<f64> {
    let mut renderer =
        MovingQuadrupoleRenderer::new(observer, sources, reader, config(), budget(), &mut || false)
            .unwrap();
    let mut output = vec![0.0; config().samples as usize];
    renderer.block(&mut output).unwrap();
    output
}
fn norm(x: [f64; 3]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}
fn close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:e}, expected={expected:e}"
    );
}

// Independent image-SOURCE construction, bisection and spatial differences.
// Production evaluates original tensors at inverse-image receivers instead.
fn potential(x: [f64; 3], t: f64, index: usize, image: [f64; 3]) -> f64 {
    let (mut lo, mut hi) = (0.0_f64, 0.3_f64);
    for _ in 0..80 {
        let tau = lo.midpoint(hi);
        let source = source_motion(index, tau);
        let r = norm(core::array::from_fn(|i| {
            x[i] - image[i] * source.position_m[i]
        }));
        if tau + r / C < t {
            lo = tau;
        } else {
            hi = tau;
        }
    }
    let tau = lo.midpoint(hi);
    let source = source_motion(index, tau);
    let r: [f64; 3] = core::array::from_fn(|i| x[i] - image[i] * source.position_m[i]);
    let distance = norm(r);
    let radial = (0..3)
        .map(|i| r[i] * image[i] * source.velocity_m_s[i])
        .sum::<f64>()
        / distance;
    let omega = core::f64::consts::TAU * [7.0, 11.0][index];
    (omega * tau + [0.2, -0.6][index]).cos()
        / (4.0 * core::f64::consts::PI * distance * (1.0 - radial / C))
}
fn hessian(x: [f64; 3], t: f64, index: usize, image: [f64; 3], h: f64) -> f64 {
    let mut total = 0.0;
    for (entry, (i, j)) in ENTRIES.into_iter().enumerate() {
        let derivative = if i == j {
            let mut plus = x;
            plus[i] += h;
            let mut minus = x;
            minus[i] -= h;
            (potential(plus, t, index, image) - 2.0 * potential(x, t, index, image)
                + potential(minus, t, index, image))
                / (h * h)
        } else {
            let mut pp = x;
            pp[i] += h;
            pp[j] += h;
            let mut pm = x;
            pm[i] += h;
            pm[j] -= h;
            let mut mp = x;
            mp[i] -= h;
            mp[j] += h;
            let mut mm = x;
            mm[i] -= h;
            mm[j] -= h;
            2.0 * (potential(pp, t, index, image)
                - potential(pm, t, index, image)
                - potential(mp, t, index, image)
                + potential(mm, t, index, image))
                / (4.0 * h * h)
        };
        // Full rank-two image, including trace and physical shear entries.
        total += TENSORS[index][entry] * image[i] * image[j] * derivative;
    }
    total
}

#[test]
fn accelerating_source_arrays_at_moving_microphones_match_independent_image_hessians() {
    let sources = [path(|t| source_motion(0, t)), path(|t| source_motion(1, t))];
    let receiver = path(receiver_motion);
    let direct: Vec<_> = sources
        .iter()
        .map(|source| {
            PreparedMovingReceiverSource::new(
                source,
                &receiver,
                pair_config(),
                pair_bounds(),
                &mut || false,
            )
            .unwrap()
        })
        .collect();
    let walls =
        MovingReflectionObserver::new(&receiver, &planes(), 4, image_bounds(), &mut || false)
            .unwrap();
    let reflected: Vec<_> = sources
        .iter()
        .map(|source| {
            PreparedReflectedMovingReceiverSource::new(
                source,
                &walls,
                pair_config(),
                image_bounds(),
                pair_bounds(),
                &mut || false,
            )
            .unwrap()
        })
        .collect();
    let direct_windows: Vec<_> = direct.iter().map(window).collect();
    let reflected_windows: Vec<_> = reflected.iter().map(window).collect();
    for (output, images) in [
        (render(&receiver, &direct_windows, read), 1usize),
        (render(&walls, &reflected_windows, read), 4usize),
    ] {
        let mut errors = [0.0_f64; 3];
        let mut signal = 0.0_f64;
        for sample in [0usize, 8, 16] {
            let t = config().first_receiver_time_s + sample as f64 / 1000.0;
            let x = receiver_motion(t).position_m;
            signal = signal.max(output[sample].abs());
            for (rung, h) in [0.08, 0.04, 0.02].into_iter().enumerate() {
                let mut expected = 0.0;
                for index in 0..2 {
                    for mask in 0..images {
                        let image = [
                            if mask & 1 == 0 { 1.0 } else { -1.0 },
                            if mask & 2 == 0 { 1.0 } else { -1.0 },
                            1.0,
                        ];
                        expected += image[1] * hessian(x, t, index, image, h);
                    }
                }
                errors[rung] = errors[rung].max((expected - output[sample]).abs());
            }
        }
        assert!(signal > 1e-8);
        assert!(
            errors[1] < 0.4 * errors[0] && errors[2] < 0.4 * errors[1],
            "{images} paths: {errors:?}"
        );
        assert!(
            errors[2] < 3e-5 * signal,
            "signal={signal:e}, errors={errors:?}"
        );
    }
    let t = 0.113;
    for motion in &reflected {
        for image in 0..motion.path_count() {
            let sample = motion
                .prepare_quadrupole_path(image, t, &mut || false)
                .unwrap();
            let original = motion.prepare_path(image, t, &mut || false).unwrap();
            close(
                sample
                    .observed_source_time_rate(receiver_motion(t).velocity_m_s)
                    .unwrap(),
                original.observed_source_time_rate,
                2e-14,
                2e-14,
            );
        }
    }
}

type Matrix = [[f64; 3]; 3];
fn transform(h: Matrix, x: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| (0..3).map(|j| h[i][j] * x[j]).sum())
}
fn tensor_image(h: Matrix, q: [f64; 6]) -> [f64; 6] {
    let q = [[q[0], q[3], q[4]], [q[3], q[1], q[5]], [q[4], q[5], q[2]]];
    ENTRIES.map(|(i, j)| {
        (0..3)
            .flat_map(|a| (0..3).map(move |b| h[i][a] * q[a][b] * h[j][b]))
            .sum()
    })
}
fn stationary_pressure(x: [f64; 3], source: [f64; 3], t: f64, tensor: [f64; 6]) -> f64 {
    let r = core::array::from_fn(|i| x[i] - source[i]);
    let distance = norm(r);
    let n = r.map(|v| v / distance);
    let q = stress(t - distance / C, tensor, core::f64::consts::TAU * 7.0, 0.2);
    let contraction = |tensor: [f64; 6]| {
        ENTRIES
            .iter()
            .enumerate()
            .map(|(entry, &(i, j))| tensor[entry] * n[i] * n[j] * if i == j { 1.0 } else { 2.0 })
            .sum::<f64>()
    };
    let a = 3.0 * contraction(q.stress_moment_n_m) - q.stress_moment_n_m[..3].iter().sum::<f64>();
    let b = 3.0 * contraction(q.stress_rate_n_m_s) - q.stress_rate_n_m_s[..3].iter().sum::<f64>();
    (a / distance.powi(3)
        + b / (C * distance * distance)
        + contraction(q.stress_acceleration_n_m_s2) / (C * C * distance))
        / (4.0 * core::f64::consts::PI)
}

#[test]
fn fixed_and_moving_wedges_transform_the_whole_tensor_once() {
    let source_position = [2.0, 0.6, 0.5];
    let microphone = |t: f64| SourceKinematics {
        position_m: [5.0 + t, 3.5 - 0.3 * t, 1.2 + 0.4 * t],
        velocity_m_s: [1.0, -0.3, 0.4],
        acceleration_m_s2: [0.0; 3],
    };
    let source = path(|_| stationary(source_position));
    let receiver = path(microphone);
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let theta = core::f64::consts::PI / 3.0;
        let planes = [
            IdealAcousticPlane::new([0.0; 3], [0.0, 1.0, 0.0], boundary).unwrap(),
            IdealAcousticPlane::new([0.0; 3], [theta.sin(), -theta.cos(), 0.0], boundary).unwrap(),
        ];
        let walls = MovingReflectionObserver::new_wedge(
            &receiver,
            planes,
            3,
            6,
            image_bounds(),
            &mut || false,
        )
        .unwrap();
        let moving = PreparedReflectedMovingReceiverSource::new(
            &source,
            &walls,
            pair_config(),
            image_bounds(),
            pair_bounds(),
            &mut || false,
        )
        .unwrap();
        let output = render(&walls, &[window(&moving)], |_, t| read(0, t));
        let fixed_position = microphone(0.0).position_m;
        let observer = MovingCompactObserver::new(
            medium(),
            fixed_position,
            MovingCompactLimits {
                maximum_mach: 0.3,
                minimum_distance_m: 0.001,
            },
        )
        .unwrap();
        let fixed_source = PreparedSourceTrajectory::new(
            observer,
            &knots(|_| stationary(source_position)),
            bounds(),
            &mut || false,
        )
        .unwrap();
        let fixed_walls = IdealReflectionObserver::new_wedge(observer, planes, 3, 6).unwrap();
        let fixed = PreparedReflectedSourceTrajectory::new(
            &fixed_source,
            &fixed_walls,
            image_bounds(),
            &mut || false,
        )
        .unwrap();
        let fixed_output = render(&fixed_walls, &[window(&fixed)], |_, t| read(0, t));
        for (sample, &actual) in output.iter().enumerate() {
            let t = config().first_receiver_time_s + sample as f64 / 1000.0;
            let mut expected = [0.0; 2];
            for rotation in 0..3 {
                let angle = 2.0 * theta * rotation as f64;
                let (s, c) = angle.sin_cos();
                for mirror in [false, true] {
                    let sign = if mirror { -1.0 } else { 1.0 };
                    let h = [[c, -s * sign, 0.0], [s, c * sign, 0.0], [0.0, 0.0, 1.0]];
                    let parity = if mirror && boundary == IdealPlaneBoundary::PressureRelease {
                        -1.0
                    } else {
                        1.0
                    };
                    let image_source = transform(h, source_position);
                    let image_tensor = tensor_image(h, TENSORS[0]);
                    for (i, x) in [microphone(t).position_m, fixed_position]
                        .into_iter()
                        .enumerate()
                    {
                        expected[i] +=
                            parity * stationary_pressure(x, image_source, t, image_tensor);
                    }
                }
            }
            close(actual, expected[0], 1e-14, 2e-9);
            close(fixed_output[sample], expected[1], 1e-14, 2e-9);
        }
        for image in 0..6 {
            let t = 0.111;
            let prepared = moving
                .prepare_quadrupole_path(image, t, &mut || false)
                .unwrap();
            let original = moving.prepare_path(image, t, &mut || false).unwrap();
            close(
                prepared
                    .observed_source_time_rate(microphone(t).velocity_m_s)
                    .unwrap(),
                original.observed_source_time_rate,
                2e-14,
                2e-14,
            );
        }
    }
}

#[test]
fn weak_shear_survives_strong_source_cancellation_before_output_rounding() {
    let source = path(|_| stationary([2.0, 1.0, 0.5]));
    let receiver = path(receiver_motion);
    let walls =
        MovingReflectionObserver::new(&receiver, &planes(), 4, image_bounds(), &mut || false)
            .unwrap();
    let motion = PreparedReflectedMovingReceiverSource::new(
        &source,
        &walls,
        pair_config(),
        image_bounds(),
        pair_bounds(),
        &mut || false,
    )
    .unwrap();
    let weak = [0.0, 0.0, 0.0, 2.0_f64.powi(-60), 0.0, 0.0];
    let mixed = [1.0, 0.0, 0.0, weak[3], 0.0, 0.0];
    let opposite = [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let expected = render(&walls, &[window(&motion)], |_, t| {
        Ok(stress(t, weak, 40.0, 0.2))
    });
    let actual = render(&walls, &[window(&motion), window(&motion)], |i, t| {
        Ok(stress(t, if i == 0 { mixed } else { opposite }, 40.0, 0.2))
    });
    assert!(expected.iter().any(|v| v.abs() > 1e-25));
    for (a, e) in actual.into_iter().zip(expected) {
        close(a, e, 1e-34, 2e-8);
    }
    let local = motion
        .prepare_quadrupole_path(0, 0.1, &mut || false)
        .unwrap();
    let t = local.source_time_s();
    let rounded = local.pressure_pa(stress(t, mixed, 40.0, 0.2)).unwrap()
        + local.pressure_pa(stress(t, opposite, 40.0, 0.2)).unwrap();
    let residual = local.pressure_pa(stress(t, weak, 40.0, 0.2)).unwrap();
    assert!(residual != 0.0 && (rounded - residual).abs() > 0.5 * residual.abs());
}

#[test]
fn complete_image_support_work_and_receiver_bandwidth_are_admitted() {
    let source = path(|t| SourceKinematics {
        position_m: [8.0 - 20.0 * t, 1.0, 0.5],
        velocity_m_s: [-20.0, 0.0, 0.0],
        acceleration_m_s2: [0.0; 3],
    });
    let receiver = path(|_| stationary([12.0, 2.0, 1.0]));
    let wall =
        IdealAcousticPlane::new([0.0; 3], [1.0, 0.0, 0.0], IdealPlaneBoundary::Rigid).unwrap();
    let walls = MovingReflectionObserver::new(&receiver, &[wall], 2, image_bounds(), &mut || false)
        .unwrap();
    let motion = PreparedReflectedMovingReceiverSource::new(
        &source,
        &walls,
        pair_config(),
        image_bounds(),
        pair_bounds(),
        &mut || false,
    )
    .unwrap();
    let good = window(&motion);
    let first = config().first_receiver_time_s;
    let direct = motion
        .prepare_quadrupole_path(0, first, &mut || false)
        .unwrap()
        .source_time_s();
    let image = motion
        .prepare_quadrupole_path(1, first, &mut || false)
        .unwrap()
        .source_time_s();
    assert!(image < direct);
    let refuse = |source, cfg, bounds| {
        assert!(
            MovingQuadrupoleRenderer::new(&walls, &[source], read, cfg, bounds, &mut || false,)
                .is_err()
        )
    };
    refuse(
        MovingQuadrupoleSourceWindow {
            source_time_interval_s: [image.midpoint(direct), 0.3],
            ..good
        },
        config(),
        budget(),
    );
    refuse(
        good,
        config(),
        QuadrupolePressureBudget {
            maximum_source_path_samples: 2 * config().samples - 1,
            ..budget()
        },
    );
    let reads = Cell::new(0);
    let accepted_sources = [good];
    let mut accepted = MovingQuadrupoleRenderer::new(
        &walls,
        &accepted_sources,
        |i, t| {
            reads.set(reads.get() + 1);
            read(i, t)
        },
        config(),
        budget(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(reads.get(), 0);
    assert_eq!(accepted.source_path_count(), 2);
    assert_eq!(
        accepted.maximum_retarded_iterations(),
        2 * 96 * config().samples
    );
    accepted.block(&mut [0.0; 1]).unwrap();
    assert_eq!(reads.get(), 2);
    // The direct path redshifts, but its image approaches the virtual receiver.
    let direct_pair = PreparedMovingReceiverSource::new(
        &source,
        &receiver,
        pair_config(),
        pair_bounds(),
        &mut || false,
    )
    .unwrap();
    let narrow = MovingPressureConfig {
        receiver_bandwidth_hz: 1.0,
        ..config()
    };
    let direct_window = MovingQuadrupoleSourceWindow {
        source_bandwidth_hz: 1.0,
        ..window(&direct_pair)
    };
    assert!(
        MovingQuadrupoleRenderer::new(
            &receiver,
            &[direct_window],
            read,
            narrow,
            budget(),
            &mut || false
        )
        .is_ok()
    );
    refuse(
        MovingQuadrupoleSourceWindow {
            source_bandwidth_hz: 1.0,
            ..good
        },
        narrow,
        budget(),
    );

    // A microphone's own admitted speed envelope need not equal its source's.
    let slow = PreparedKinematicTrajectory::new(
        medium(),
        &knots(|_| stationary([2.0, 1.0, 0.5])),
        0.01,
        bounds(),
        &mut || false,
    )
    .unwrap();
    let fast = path(|t| SourceKinematics {
        position_m: [10.0, 6.0 + 80.0 * t, 3.0],
        velocity_m_s: [0.0, 80.0, 0.0],
        acceleration_m_s2: [0.0; 3],
    });
    let pair =
        PreparedMovingReceiverSource::new(&slow, &fast, pair_config(), pair_bounds(), &mut || {
            false
        })
        .unwrap();
    let prepared = pair
        .prepare_quadrupole_receiver(0.12, &mut || false)
        .unwrap();
    let original = pair.prepare_receiver(0.12, &mut || false).unwrap();
    close(
        prepared
            .observed_source_time_rate([0.0, 80.0, 0.0])
            .unwrap(),
        original.observed_source_time_rate,
        2e-14,
        2e-14,
    );
}

#[test]
fn last_image_failure_keeps_all_output_and_original_source_clocks_transactional() {
    let sources = [path(|t| source_motion(0, t)), path(|t| source_motion(1, t))];
    let receiver = path(receiver_motion);
    let walls =
        MovingReflectionObserver::new(&receiver, &planes(), 4, image_bounds(), &mut || false)
            .unwrap();
    let motions: Vec<_> = sources
        .iter()
        .map(|source| {
            PreparedReflectedMovingReceiverSource::new(
                source,
                &walls,
                pair_config(),
                image_bounds(),
                pair_bounds(),
                &mut || false,
            )
            .unwrap()
        })
        .collect();
    let windows: Vec<_> = motions.iter().map(window).collect();
    let expected = render(&walls, &windows, read);
    let fail = Cell::new(true);
    let calls = RefCell::new(Vec::new());
    let reader = |i, t| {
        calls.borrow_mut().push((i, t));
        if fail.get() && calls.borrow().len() == 8 {
            Err(AeroacError::InvalidParameter {
                what: "last image stress unavailable",
            })
        } else {
            read(i, t)
        }
    };
    let mut renderer =
        MovingQuadrupoleRenderer::new(&walls, &windows, reader, config(), budget(), &mut || false)
            .unwrap();
    let mut output = [123.0; 17];
    assert!(renderer.block(&mut output).is_err());
    assert_eq!(output, [123.0; 17]);
    assert_eq!(renderer.samples_rendered(), 0);
    assert_eq!(
        calls.borrow().iter().map(|v| v.0).collect::<Vec<_>>(),
        vec![0, 0, 0, 0, 1, 1, 1, 1]
    );
    for (source_index, motion) in motions.iter().enumerate() {
        for image in 0..4 {
            let expected_time = motion
                .prepare_quadrupole_path(image, 0.1, &mut || false)
                .unwrap()
                .source_time_s();
            assert_eq!(calls.borrow()[source_index * 4 + image].1, expected_time);
        }
    }
    fail.set(false);
    let mut polls = 0;
    assert!(
        renderer
            .block_controlled(&mut output, &mut || {
                polls += 1;
                polls == 25
            })
            .is_err()
    );
    assert_eq!(output, [123.0; 17]);
    assert_eq!(renderer.samples_rendered(), 0);
    for chunk in output.chunks_mut(5) {
        renderer.block(chunk).unwrap();
    }
    assert_eq!(output.as_slice(), expected.as_slice());
    assert_eq!(renderer.samples_rendered(), 17);
}
