use super::*;
use fs_render::camera::CoherentDiffractionInterferometer;

fn receiver(
    coordinates: [[f64; 2]; 2],
    powers: [f64; 2],
    path_difference: f64,
    phase: f64,
) -> CoherentDiffractionInterferometer {
    CoherentDiffractionInterferometer::try_new(coordinates, 0.03, powers, path_difference, phase)
        .unwrap()
}

#[test]
fn g1_path_phase_passive_output_closure_and_weak_dark_port() {
    let camera = focused();
    let mean = [wave([0.0; 2], C64::ONE)];
    let query = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &mean,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let one = core::f64::consts::PI * RADIUS * RADIUS / (LAMBDA * LAMBDA);
    for (path, phase, plus_fraction) in [
        (0.0, 0.0, 1.0),
        (0.25 * LAMBDA, 0.0, 0.5),
        (0.5 * LAMBDA, 0.0, 0.0),
        (0.0, core::f64::consts::PI, 0.0),
    ] {
        let sample = query
            .observe_interferometer(receiver([[0.0; 2]; 2], [0.5; 2], path, phase), 2, || false)
            .unwrap();
        close(sample.arm_power_density_w[0], 0.5 * one, 3e-12);
        close(sample.arm_power_density_w[1], 0.5 * one, 3e-12);
        close(sample.output_power_density_w.iter().sum(), one, 3e-12);
        assert!((sample.output_power_density_w[0] / one - plus_fraction).abs() < 3e-12);
        close(sample.fringe_visibility, 1.0, 1e-14);
        assert_eq!(sample.pupil_evaluations, 2);
        if path == 0.0 && phase == 0.0 {
            assert_eq!(sample.output_power_density_w[1], 0.0);
        }
    }
    let unequal = query
        .observe_interferometer(receiver([[0.0; 2]; 2], [0.25, 0.75], 0.0, 0.0), 2, || false)
        .unwrap();
    close(unequal.fringe_visibility, 3.0_f64.sqrt() / 2.0, 1e-14);
    close(unequal.output_power_density_w.iter().sum(), one, 3e-12);

    // The phase is far below one ulp of the bright intensity. A rounded-power
    // subtraction would erase this physically positive dark-port signal.
    let phase = 1e-20;
    let weak = query
        .observe_interferometer(receiver([[0.0; 2]; 2], [0.5; 2], 0.0, phase), 2, || false)
        .unwrap();
    assert!(weak.output_power_density_w[1] > 0.0);
    close(
        weak.output_power_density_w[1] / (0.25 * one * phase * phase),
        1.0,
        3e-12,
    );
}

#[test]
fn g3_whole_wavelength_paths_preserve_dark_ports_and_tiny_authored_phase() {
    let camera = focused();
    let mean = [wave([0.0; 2], C64::ONE)];
    // 10^9/512 = 1,953,125 exactly: every integer metre is a whole number
    // of these authored vacuum wavelengths, without a rounded SI conversion.
    let query = camera
        .prepare_coherent_diffraction(
            512.0,
            CoherentPupilIllumination {
                mean: &mean,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let one = query.observe([0.0; 2], || false).unwrap().power_density_w;
    for power in [0, 20, 40] {
        for sign in [-1.0, 1.0] {
            let path = sign * 2.0_f64.powi(power);
            for phase in [0.0, 1e-20] {
                let sample = query
                    .observe_interferometer(
                        receiver([[0.0; 2]; 2], [0.5; 2], path, phase),
                        2,
                        || false,
                    )
                    .unwrap();
                close(sample.output_power_density_w.iter().sum(), one, 3e-12);
                close(sample.fringe_visibility, 1.0, 1e-14);
                if phase == 0.0 {
                    assert_eq!(sample.output_power_density_w[1], 0.0);
                    assert_eq!(sample.mutual_coherence_density_w.im, 0.0);
                } else {
                    assert!(sample.output_power_density_w[1] > 0.0);
                    close(
                        sample.output_power_density_w[1] / (0.25 * one * phase * phase),
                        1.0,
                        3e-12,
                    );
                    assert!(sample.mutual_coherence_density_w.im < 0.0);
                }
            }
        }
    }
}

#[test]
fn g3_whole_wavelength_sensor_chirp_preserves_the_spatial_cross_phase() {
    let camera = focused();
    let mean = [wave([0.0; 2], C64::ONE)];
    let query = camera
        .prepare_coherent_diffraction(
            512.0,
            CoherentPupilIllumination {
                mean: &mean,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let center = [0.0; 2];
    let off_axis = [2.0_f64.powi(-10), 0.0];
    let a = query.observe(center, || false).unwrap().mean_field;
    let b = query.observe(off_axis, || false).unwrap().mean_field;
    assert_ne!(a.re * b.re, 0.0);
    for power in [21, 41, 61] {
        // d*theta^2/2 is exactly 1, 2^20 or 2^40 metres. Swapping the
        // coordinates reverses the chirp; both signs remain whole cycles.
        for coordinates in [[center, off_axis], [off_axis, center]] {
            let receiver = CoherentDiffractionInterferometer::try_new(
                coordinates,
                2.0_f64.powi(power),
                [0.5; 2],
                0.0,
                0.0,
            )
            .unwrap();
            let sample = query.observe_interferometer(receiver, 2, || false).unwrap();
            assert_eq!(sample.mutual_coherence_density_w.im, 0.0);
            close(
                sample.mutual_coherence_density_w.re,
                0.5 * a.re * b.re,
                2e-12,
            );
            close(
                sample.output_power_density_w[0],
                0.25 * (a + b).norm_sq(),
                2e-12,
            );
            close(
                sample.output_power_density_w[1],
                0.25 * (a - b).norm_sq(),
                2e-12,
            );
        }
    }
}

#[test]
fn g1_mutual_spatial_coherence_controls_visibility() {
    let camera = focused();
    let left = [wave([-1.0, 0.0], C64::ONE)];
    let right = [wave([1.0, 0.0], C64::ONE)];
    let factors = [
        PupilCoherenceFactor { waves: &left },
        PupilCoherenceFactor { waves: &right },
    ];
    let query = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &[],
                factors: &factors,
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let coordinates = [[-2.0 * SCALE, 0.0], [2.0 * SCALE, 0.0]];
    let receiver = receiver(coordinates, [0.5; 2], 0.0, 0.0);
    let sample = query.observe_interferometer(receiver, 4, || false).unwrap();
    let single = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &left,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let a = single
        .observe(coordinates[0], || false)
        .unwrap()
        .mean_field
        .re;
    let b = single
        .observe(coordinates[1], || false)
        .unwrap()
        .mean_field
        .re;
    close(
        sample.fringe_visibility,
        2.0 * (a * b).abs() / (a * a + b * b),
        2e-12,
    );
    assert!(sample.fringe_visibility < 0.8);
    close(
        sample.output_power_density_w[0],
        0.5 * (a + b).powi(2),
        2e-12,
    );
    close(
        sample.output_power_density_w[1],
        0.5 * (a - b).powi(2),
        2e-12,
    );
    close(sample.mutual_coherence_density_w.re, a * b, 2e-12);
    assert_eq!(sample.mutual_coherence_density_w.im, 0.0);

    let coherent = [left[0], right[0]];
    let query = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &coherent,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    close(
        query
            .observe_interferometer(receiver, 4, || false)
            .unwrap()
            .fringe_visibility,
        1.0,
        1e-14,
    );
}

#[test]
fn g1_sensor_chirp_and_authored_optical_path_retain_complex_pupil_phase() {
    let camera = defocused(2.4);
    let mean = [wave([-1.0, 0.0], C64::new(0.8, -0.3))];
    let query = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &mean,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let coordinates = [[0.0; 2], [0.003, 0.0]];
    let path = 0.13 * LAMBDA;
    let phase = 0.27;
    let sample = query
        .observe_interferometer(receiver(coordinates, [0.5; 2], path, phase), 2, || false)
        .unwrap();
    let first = query.observe(coordinates[0], || false).unwrap().mean_field;
    let second = query.observe(coordinates[1], || false).unwrap().mean_field;
    let relative =
        core::f64::consts::TAU / LAMBDA * (0.03 * 0.003_f64.powi(2) / 2.0 + path) + phase;
    let rotated = second * C64::new(relative.cos(), relative.sin());
    close(
        sample.output_power_density_w[0],
        0.25 * (first + rotated).norm_sq(),
        2e-12,
    );
    close(
        sample.output_power_density_w[1],
        0.25 * (first - rotated).norm_sq(),
        2e-12,
    );
    let mutual = (first * rotated.conj()).scale(0.5);
    close(sample.mutual_coherence_density_w.re, mutual.re, 2e-12);
    close(sample.mutual_coherence_density_w.im, mutual.im, 2e-12);
    assert_ne!(sample.mutual_coherence_density_w.im, 0.0);
}

#[test]
fn g0_g4_passivity_sensor_work_and_cancellation_refusals() {
    for powers in [
        [0.6, 0.6],
        [-0.1, 0.5],
        [f64::NAN, 0.0],
        [1.0, f64::EPSILON],
    ] {
        assert_eq!(
            CoherentDiffractionInterferometer::try_new([[0.0; 2]; 2], 0.03, powers, 0.0, 0.0),
            Err(CoherentDiffractionError::InvalidInterferometer),
        );
    }
    let camera = focused();
    let mean = [wave([0.0; 2], C64::ONE)];
    let query = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &mean,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let balanced = receiver([[0.0; 2]; 2], [0.5; 2], 0.0, 0.0);
    let mut polls = 0;
    assert!(matches!(
        query.observe_interferometer(balanced, 1, || {
            polls += 1;
            false
        }),
        Err(CoherentDiffractionError::Limit {
            quantity: "interferometer pupil evaluations",
            requested: 2,
            maximum: 1
        }),
    ));
    assert_eq!(polls, 1);
    assert_eq!(
        query.observe_interferometer(balanced, 2, || true),
        Err(CoherentDiffractionError::Cancelled)
    );
    close(
        query
            .observe_interferometer(balanced, 2, || false)
            .unwrap()
            .fringe_visibility,
        1.0,
        1e-14,
    );
    let blocked = query
        .observe_interferometer(receiver([[0.0; 2]; 2], [1.0, 0.0], 0.0, 0.0), 1, || false)
        .unwrap();
    assert_eq!(blocked.pupil_evaluations, 1);
    assert_eq!(blocked.fringe_visibility, 0.0);
    assert_eq!(
        blocked.output_power_density_w[0],
        blocked.output_power_density_w[1]
    );

    let physical = PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, -1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.024).unwrap(),
        2.0,
        Aperture::try_circular(RADIUS).unwrap(),
    )
    .unwrap()
    .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
    .unwrap();
    let query = physical
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &mean,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    assert_eq!(
        query.observe_interferometer(balanced, 2, || false),
        Err(CoherentDiffractionError::InvalidSensorPlane)
    );
    let conjugate = CoherentDiffractionInterferometer::try_new(
        [[0.0; 2]; 2],
        0.05 / (1.0 - 0.05 / 2.0),
        [0.5; 2],
        0.0,
        0.0,
    )
    .unwrap();
    query
        .observe_interferometer(conjugate, 2, || false)
        .unwrap();
}
