use super::*;
use fs_render::camera::{
    CoherentDiffractionInterferometer, CoherentJonesArm, CoherentJonesPupilLayer,
    MAX_COHERENT_JONES_ARM_LAYERS, PolarizedCoherentInterferometerSample,
};
use fs_render::polarization::{Mueller, Stokes};

fn receiver(phase: f64) -> CoherentDiffractionInterferometer {
    CoherentDiffractionInterferometer::try_new([[0.0; 2]; 2], 0.03, [0.5; 2], 0.0, phase).unwrap()
}

fn layer(power: [f64; 2], phase: [f64; 2], angle: f64) -> CoherentJonesPupilLayer {
    CoherentJonesPupilLayer::try_new(power, phase, angle, [0.0; 2]).unwrap()
}

fn arm(layers: &[CoherentJonesPupilLayer]) -> CoherentJonesArm {
    CoherentJonesArm::try_new(layers).unwrap()
}

fn observe(
    camera: &PhysicalCamera,
    mean: &[CoherentJonesPupilWave],
    factors: &[JonesPupilCoherenceFactor<'_>],
    receiver: CoherentDiffractionInterferometer,
    arms: [CoherentJonesArm; 2],
) -> PolarizedCoherentInterferometerSample {
    let query = camera
        .prepare_polarized_coherent_diffraction(
            LAMBDA_NM,
            CoherentJonesPupilIllumination { mean, factors },
            LIMITS,
            || false,
        )
        .unwrap();
    query
        .observe_interferometer(receiver, arms, 2 * query.wave_count(), || false)
        .unwrap()
}

fn components(s: Stokes) -> [f64; 4] {
    [s.i, s.q, s.u, s.v]
}

fn stokes_close(actual: Stokes, expected: Stokes, tolerance: f64) {
    for (a, e) in components(actual).into_iter().zip(components(expected)) {
        assert!(
            (a - e).abs() <= tolerance * expected.i.max(1e-280),
            "{actual:?} != {expected:?}"
        );
    }
}

fn closure(sample: PolarizedCoherentInterferometerSample) {
    let outputs = sample.output_coherency_density_w.map(|c| c.stokes().i);
    let scale = outputs
        .into_iter()
        .chain(sample.arm_power_density_w)
        .fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return;
    }
    let incoming = sample
        .arm_power_density_w
        .map(|p| p / scale)
        .iter()
        .sum::<f64>();
    let outgoing = outputs.map(|p| p / scale).iter().sum::<f64>();
    assert!(
        (incoming - outgoing).abs() < 3e-12,
        "{incoming} != {outgoing}"
    );
}

#[test]
fn g1_orthogonal_arms_keep_phase_in_polarization_and_analyzer_restores_lossy_fringe() {
    let camera = focused();
    let peak = intensity(&camera, &[wave([0.0; 2], C64::ONE)], &[], [0.0; 2]);
    let mean = [jones_wave([0.0; 2], [C64::ONE; 2])];
    let x = layer([1.0, 0.0], [0.0; 2], 0.0);
    let y = layer([0.0, 1.0], [0.0; 2], 0.0);
    let analyzer = layer([1.0, 0.0], [0.0; 2], core::f64::consts::FRAC_PI_4);
    for phase in [0.0, core::f64::consts::FRAC_PI_2, core::f64::consts::PI] {
        let bare = observe(&camera, &mean, &[], receiver(phase), [arm(&[x]), arm(&[y])]);
        assert_eq!(bare.mutual_coherence_density_w, C64::ZERO);
        assert_eq!(bare.fringe_visibility, 0.0);
        for port in 0..2 {
            let sign = if port == 0 { 1.0 } else { -1.0 };
            stokes_close(
                bare.output_coherency_density_w[port].stokes(),
                Stokes {
                    i: 0.5 * peak,
                    q: 0.0,
                    u: sign * 0.5 * peak * phase.cos(),
                    v: -sign * 0.5 * peak * phase.sin(),
                },
                3e-13,
            );
            close(bare.arm_power_density_w[port], 0.5 * peak, 3e-13);
        }
        closure(bare);
        let analysed = observe(
            &camera,
            &mean,
            &[],
            receiver(phase),
            [arm(&[x, analyzer]), arm(&[y, analyzer])],
        );
        close(analysed.fringe_visibility, 1.0, 1e-14);
        close(analysed.arm_power_density_w.iter().sum(), 0.5 * peak, 3e-13);
        for port in 0..2 {
            let expected = Mueller::linear_polarizer(core::f64::consts::FRAC_PI_4)
                .apply(bare.output_coherency_density_w[port].stokes());
            let actual = analysed.output_coherency_density_w[port].stokes();
            for (a, e) in components(actual).into_iter().zip(components(expected)) {
                assert!((a - e).abs() < 3e-13 * peak, "{actual:?} != {expected:?}");
            }
        }
        if phase == 0.0 {
            assert!(analysed.output_coherency_density_w[1].stokes().i < peak * 1e-28);
        }
        closure(analysed);
    }
    // Projection restores an interference fringe only for a source bank that
    // supplies both arms. Independent x/y populations retain no mutual phase.
    let x_source = [jones_wave([0.0; 2], [C64::ONE, C64::ZERO])];
    let y_source = [jones_wave([0.0; 2], [C64::ZERO, C64::ONE])];
    let factors = [
        JonesPupilCoherenceFactor { waves: &x_source },
        JonesPupilCoherenceFactor { waves: &y_source },
    ];
    for phase in [0.0, core::f64::consts::PI] {
        let sample = observe(
            &camera,
            &[],
            &factors,
            receiver(phase),
            [arm(&[x, analyzer]), arm(&[y, analyzer])],
        );
        assert_eq!(sample.fringe_visibility, 0.0);
        assert_eq!(sample.mutual_coherence_density_w, C64::ZERO);
        for output in sample.output_coherency_density_w {
            close(output.stokes().i, 0.25 * peak, 3e-13);
        }
        closure(sample);
    }
}

#[test]
fn g1_retardance_circular_sign_and_common_arm_phase_are_retained() {
    let camera = focused();
    let peak = intensity(&camera, &[wave([0.0; 2], C64::ONE)], &[], [0.0; 2]);
    let mean = [jones_wave([0.0; 2], [C64::ONE; 2])];
    for phase in [core::f64::consts::FRAC_PI_2, -core::f64::consts::FRAC_PI_2] {
        let retarder = arm(&[layer([1.0; 2], [0.0, phase], 0.0)]);
        let sample = observe(&camera, &mean, &[], receiver(0.0), [retarder; 2]);
        let expected = Mueller::linear_retarder(0.0, phase).apply(Stokes {
            i: 2.0 * peak,
            q: 0.0,
            u: 2.0 * peak,
            v: 0.0,
        });
        stokes_close(
            sample.output_coherency_density_w[0].stokes(),
            expected,
            3e-13,
        );
        assert_eq!(sample.output_coherency_density_w[1].stokes().i, 0.0);
        closure(sample);
    }
    let quarter = arm(&[layer([1.0; 2], [0.0, core::f64::consts::FRAC_PI_2], 0.0)]);
    let sample = observe(
        &camera,
        &mean,
        &[],
        receiver(0.0),
        [CoherentJonesArm::IDENTITY, quarter],
    );
    close(
        sample.fringe_visibility,
        core::f64::consts::FRAC_1_SQRT_2,
        3e-13,
    );
    close(sample.mutual_coherence_density_w.re, 0.5 * peak, 3e-13);
    close(sample.mutual_coherence_density_w.im, -0.5 * peak, 3e-13);
    stokes_close(
        sample.output_coherency_density_w[0].stokes(),
        Stokes {
            i: 1.5 * peak,
            q: 0.5 * peak,
            u: peak,
            v: -peak,
        },
        3e-13,
    );
    stokes_close(
        sample.output_coherency_density_w[1].stokes(),
        Stokes {
            i: 0.5 * peak,
            q: -0.5 * peak,
            u: 0.0,
            v: 0.0,
        },
        3e-13,
    );
    closure(sample);

    let phases = [0.3, 0.7].map(|p| arm(&[layer([1.0; 2], [p; 2], 0.0)]));
    let sample = observe(&camera, &mean, &[], receiver(0.0), phases);
    close(
        sample.output_coherency_density_w[0].stokes().i,
        2.0 * peak * 0.2_f64.cos().powi(2),
        3e-13,
    );
    close(
        sample.output_coherency_density_w[1].stokes().i,
        2.0 * peak * 0.2_f64.sin().powi(2),
        3e-13,
    );
    closure(sample);
}

#[test]
fn g3_independent_spatial_covariance_and_its_basis_preserve_both_ports() {
    let camera = defocused(1.7);
    let left = [jones_wave([-1.0, 0.2], [C64::ONE, C64::new(0.0, 0.4)])];
    let right = [jones_wave(
        [0.8, -0.3],
        [C64::new(0.3, -0.2), C64::new(0.7, 0.0)],
    )];
    let factors = [
        JonesPupilCoherenceFactor { waves: &left },
        JonesPupilCoherenceFactor { waves: &right },
    ];
    let receiver = CoherentDiffractionInterferometer::try_new(
        [[-2.0 * SCALE, 0.0], [1.4 * SCALE, 0.3 * SCALE]],
        0.03,
        [0.4, 0.6],
        0.17 * LAMBDA,
        0.2,
    )
    .unwrap();
    let arms = [
        arm(&[layer([0.8, 0.2], [0.2, -0.3], 0.31)]),
        arm(&[layer([1.0; 2], [0.0, 0.8], -0.2)]),
    ];
    let left_sample = observe(&camera, &left, &[], receiver, arms);
    let right_sample = observe(&camera, &right, &[], receiver, arms);
    let total = observe(&camera, &[], &factors, receiver, arms);
    for port in 0..2 {
        let a = left_sample.output_coherency_density_w[port].stokes();
        let b = right_sample.output_coherency_density_w[port].stokes();
        stokes_close(
            total.output_coherency_density_w[port].stokes(),
            Stokes {
                i: a.i + b.i,
                q: a.q + b.q,
                u: a.u + b.u,
                v: a.v + b.v,
            },
            3e-12,
        );
        close(
            total.arm_power_density_w[port],
            left_sample.arm_power_density_w[port] + right_sample.arm_power_density_w[port],
            3e-12,
        );
    }
    let trace = left_sample.mutual_coherence_density_w + right_sample.mutual_coherence_density_w;
    close(total.mutual_coherence_density_w.re, trace.re, 3e-12);
    close(total.mutual_coherence_density_w.im, trace.im, 3e-12);
    close(
        total.fringe_visibility,
        2.0 * trace.abs() / total.arm_power_density_w.iter().sum::<f64>(),
        3e-12,
    );
    let scaled = |mut w: CoherentJonesPupilWave, scale: f64| {
        w.amplitude_sqrt_w = w.amplitude_sqrt_w.map(|a| a.scale(scale));
        w
    };
    let mixed_a = [scaled(left[0], 0.6), scaled(right[0], 0.8)];
    let mixed_b = [scaled(left[0], -0.8), scaled(right[0], 0.6)];
    let mixed_factors = [
        JonesPupilCoherenceFactor { waves: &mixed_a },
        JonesPupilCoherenceFactor { waves: &mixed_b },
    ];
    let mixed = observe(&camera, &[], &mixed_factors, receiver, arms);
    for port in 0..2 {
        stokes_close(
            mixed.output_coherency_density_w[port].stokes(),
            total.output_coherency_density_w[port].stokes(),
            3e-12,
        );
    }
    close(mixed.fringe_visibility, total.fringe_visibility, 3e-12);
    closure(total);
    closure(mixed);
}

#[test]
fn g1_scalar_limit_uses_original_spatial_chirp_and_whole_cycle_path_owner() {
    let camera = defocused(2.4);
    let scalar = [
        wave([-0.7, 0.1], C64::new(0.8, -0.3)),
        wave([0.6, -0.2], C64::new(-0.2, 0.4)),
    ];
    let mean = scalar.map(|w| CoherentJonesPupilWave {
        tangent_direction: w.tangent_direction,
        amplitude_sqrt_w: [w.amplitude_sqrt_w, C64::ZERO],
    });
    let receiver = CoherentDiffractionInterferometer::try_new(
        [[0.0; 2], [0.003, -0.0001]],
        0.03,
        [0.5; 2],
        0.13 * LAMBDA,
        0.27,
    )
    .unwrap();
    let original = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &scalar,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap()
        .observe_interferometer(receiver, 4, || false)
        .unwrap();
    let actual = observe(
        &camera,
        &mean,
        &[],
        receiver,
        [CoherentJonesArm::IDENTITY; 2],
    );
    for port in 0..2 {
        close(
            actual.output_coherency_density_w[port].axis_powers()[0],
            original.output_power_density_w[port],
            3e-13,
        );
        assert_eq!(
            actual.output_coherency_density_w[port].axis_powers()[1],
            0.0
        );
        close(
            actual.arm_power_density_w[port],
            original.arm_power_density_w[port],
            3e-13,
        );
    }
    close(actual.fringe_visibility, original.fringe_visibility, 3e-13);
    close(
        actual.mutual_coherence_density_w.re,
        original.mutual_coherence_density_w.re,
        3e-13,
    );
    close(
        actual.mutual_coherence_density_w.im,
        original.mutual_coherence_density_w.im,
        3e-13,
    );
    closure(actual);

    let mean = [jones_wave([0.0; 2], [C64::ONE, C64::ZERO])];
    let query = camera
        .prepare_polarized_coherent_diffraction(
            512.0,
            CoherentJonesPupilIllumination {
                mean: &mean,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    for path in [0.0, 1.0, 2.0_f64.powi(30)] {
        let receiver =
            CoherentDiffractionInterferometer::try_new([[0.0; 2]; 2], 0.03, [0.5; 2], path, 0.0)
                .unwrap();
        let sample = query
            .observe_interferometer(receiver, [CoherentJonesArm::IDENTITY; 2], 2, || false)
            .unwrap();
        assert_eq!(sample.output_coherency_density_w[1].stokes().i, 0.0);
        closure(sample);
    }
}

#[test]
fn g1_weak_principal_populations_and_dark_port_survive_before_stokes_rounding() {
    let camera = focused();
    let peak = intensity(&camera, &[wave([0.0; 2], C64::ONE)], &[], [0.0; 2]);
    let x = [jones_wave([0.0; 2], [C64::ONE, C64::ZERO])];
    let y = [jones_wave([0.0; 2], [C64::ZERO, C64::ONE])];
    let factors = [
        JonesPupilCoherenceFactor { waves: &x },
        JonesPupilCoherenceFactor { waves: &y },
    ];
    for leakage in [1e-20_f64, 1e-80, 1e-200] {
        let first = layer([1.0, leakage], [0.0; 2], 0.0);
        let second = layer([leakage, 1.0], [0.0; 2], 0.0);
        let selected = observe(&camera, &[], &factors, receiver(0.0), [arm(&[first]); 2]);
        let s = selected.output_coherency_density_w[0].stokes();
        assert_eq!(s.i, s.q);
        close(
            selected.output_coherency_density_w[0].axis_powers()[1] / (peak * leakage),
            1.0,
            3e-12,
        );
        let crossed = observe(
            &camera,
            &[],
            &factors,
            receiver(0.0),
            [arm(&[first, second]); 2],
        );
        for power in crossed.output_coherency_density_w[0].axis_powers() {
            close(power / (peak * leakage), 1.0, 3e-12);
        }
        assert_eq!(crossed.output_coherency_density_w[1].stokes().i, 0.0);
        closure(crossed);
    }
    let leakage = 1e-80_f64;
    let phase = 1e-20;
    let mean = [jones_wave(
        [0.0; 2],
        [C64::ONE, C64::new(leakage.sqrt(), 0.0)],
    )];
    let sample = observe(
        &camera,
        &mean,
        &[],
        receiver(phase),
        [CoherentJonesArm::IDENTITY; 2],
    );
    let dark = sample.output_coherency_density_w[1];
    for (power, fraction) in dark.axis_powers().into_iter().zip([1.0, leakage]) {
        assert!(power > 0.0);
        close(power / (0.25 * peak * phase * phase * fraction), 1.0, 3e-12);
    }
    assert_eq!(dark.stokes().i, dark.stokes().q);
    closure(sample);
}

#[test]
fn g1_finite_bright_ports_survive_combiner_and_total_arm_overflow() {
    let camera = focused();
    let peak = intensity(&camera, &[wave([0.0; 2], C64::ONE)], &[], [0.0; 2]);
    let target = 0.75 * f64::MAX;
    let amplitude = C64::new((target / peak).sqrt(), 0.0);
    let mean = [jones_wave([0.0; 2], [amplitude, C64::ZERO])];
    let sample = observe(
        &camera,
        &mean,
        &[],
        receiver(0.0),
        [CoherentJonesArm::IDENTITY; 2],
    );
    close(
        sample.output_coherency_density_w[0].stokes().i / target,
        1.0,
        3e-13,
    );
    assert_eq!(sample.output_coherency_density_w[1].stokes().i, 0.0);
    closure(sample);

    let mean = [jones_wave([0.0; 2], [amplitude; 2])];
    let half_wave = arm(&[layer([1.0; 2], [0.0, core::f64::consts::PI], 0.0)]);
    let sample = observe(
        &camera,
        &mean,
        &[],
        receiver(0.0),
        [CoherentJonesArm::IDENTITY, half_wave],
    );
    assert!(!sample.arm_power_density_w.iter().sum::<f64>().is_finite());
    for output in sample.output_coherency_density_w {
        close(output.stokes().i / target, 1.0, 3e-13);
    }
    assert!(sample.fringe_visibility < 2e-15);
    closure(sample);
}

#[test]
fn g0_g4_arm_admission_work_underflow_and_cancellation_are_fail_closed() {
    assert!(matches!(
        CoherentJonesPupilLayer::try_new([1.01, 1.0], [0.0; 2], 0.0, [0.0; 2]),
        Err(CoherentDiffractionError::InvalidJonesPupil)
    ));
    let patterned =
        CoherentJonesPupilLayer::try_new([1.0; 2], [0.0; 2], 0.0, [0.001, 0.0]).unwrap();
    assert_eq!(
        CoherentJonesArm::try_new(&[patterned]),
        Err(CoherentDiffractionError::InvalidJonesPupil)
    );
    let identity = layer([1.0; 2], [0.0; 2], 0.0);
    assert!(matches!(
        CoherentJonesArm::try_new(&[identity; MAX_COHERENT_JONES_ARM_LAYERS + 1]),
        Err(CoherentDiffractionError::Limit {
            quantity: "Jones arm layers",
            ..
        })
    ));
    let camera = focused();
    let mean = [jones_wave([0.0; 2], [C64::ONE, C64::ZERO])];
    let query = camera
        .prepare_polarized_coherent_diffraction(
            LAMBDA_NM,
            CoherentJonesPupilIllumination {
                mean: &mean,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        query.observe_interferometer(receiver(0.0), [CoherentJonesArm::IDENTITY; 2], 1, || {
            polls += 1;
            false
        }),
        Err(CoherentDiffractionError::Limit {
            requested: 2,
            maximum: 1,
            ..
        })
    ));
    assert_eq!(polls, 1);
    let mut polls = 0;
    assert_eq!(
        query.observe_interferometer(receiver(0.0), [CoherentJonesArm::IDENTITY; 2], 2, || {
            polls += 1;
            polls == 3
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(polls, 3);
    query
        .observe_interferometer(receiver(0.0), [CoherentJonesArm::IDENTITY; 2], 2, || false)
        .unwrap();
    let dark = CoherentDiffractionInterferometer::try_new([[0.0; 2]; 2], 0.03, [0.0; 2], 0.0, 0.0)
        .unwrap();
    let dark = query
        .observe_interferometer(dark, [CoherentJonesArm::IDENTITY; 2], 0, || false)
        .unwrap();
    assert_eq!(dark.pupil_evaluations, 0);
    assert_eq!(
        dark.output_coherency_density_w[0].stokes(),
        Stokes::unpolarized(0.0)
    );
    assert_eq!(dark.fringe_visibility, 0.0);

    let mean = [jones_wave([0.0; 2], [C64::new(1e-300, 0.0), C64::ZERO])];
    let query = camera
        .prepare_polarized_coherent_diffraction(
            LAMBDA_NM,
            CoherentJonesPupilIllumination {
                mean: &mean,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let weak = arm(&[layer([1e-80, 1.0], [0.0; 2], 0.0)]);
    assert_eq!(
        query.observe_interferometer(receiver(0.0), [weak; 2], 2, || false),
        Err(CoherentDiffractionError::Unrepresentable)
    );
}
