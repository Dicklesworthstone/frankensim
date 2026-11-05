use super::super::harmonic::{Field, Options as HarmonicOptions, Prepared};
use super::*;
use fs_math::c64::C64;

fn medium() -> RayleighMedium {
    RayleighMedium {
        density: 1.2,
        sound_speed: 320.0,
    }
}
fn source(x: f64, shape: f64) -> SurfaceSample {
    SurfaceSample {
        position_m: [x, 0.0, 0.0],
        area_m2: 0.01,
        mode_shape: vec![shape],
    }
}
fn phase(x: f64) -> C64 {
    C64::new(x.cos(), x.sin())
}
fn close(a: C64, b: C64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * b.abs().max(1e-200),
        "{a:?} versus {b:?}"
    );
}
fn harmonic(surface: &[SurfaceSample], aperture: &Aperture, hz: f64) -> Prepared {
    Prepared::new_aperture(
        surface,
        1,
        aperture,
        hz,
        medium(),
        HarmonicOptions::default(),
        || false,
    )
    .unwrap()
}

#[test]
fn g0_physical_area_normalization_and_rotated_disk_geometry() {
    let center = [0.1, -0.2, 0.4];
    let normal = [1.0, 2.0, -3.0];
    let aperture = Aperture::circular(center, normal, 0.03, 4, 16).unwrap();
    assert!((aperture.area_m2() - PI * 0.03_f64.powi(2)).abs() < 1e-18);
    for (a, b) in aperture.center_m().iter().zip(center) {
        assert!((a - b).abs() < 1e-15);
    }
    let mut weights = Dd::ZERO;
    for (&weight, sample) in aperture.weights().iter().zip(aperture.samples()) {
        weights = weights + weight;
        let delta = std::array::from_fn::<_, 3, _>(|i| sample.position_m[i] - center[i]);
        let dot: f64 = delta.iter().zip(normal).map(|(x, n)| x * n).sum();
        assert!(dot.abs() < 2e-16);
        assert!(delta.iter().map(|x| x * x).sum::<f64>() < 0.03_f64.powi(2));
    }
    assert!((weights - Dd::ONE).to_f64().abs() < 1e-30);
    // A change of physical area units scales numerator and denominator equally.
    let scaled = Aperture::new(
        aperture
            .samples()
            .iter()
            .map(|s| Sample {
                position_m: s.position_m,
                area_m2: s.area_m2 * 1e6,
            })
            .collect(),
    )
    .unwrap();
    assert_eq!(aperture.weights(), scaled.weights());
}

#[test]
fn g3_single_point_delegates_bit_for_bit_for_every_field_bank_and_clock() {
    let surface = [source(-0.2, 1.0), source(0.12, -0.3)];
    let position = [0.03, -0.1, 0.6];
    let aperture = Aperture::new(vec![Sample {
        position_m: position,
        area_m2: 0.001,
    }])
    .unwrap();
    for ratio in [1, 3, 4] {
        for alpha in [0.0, 0.5, 1.0] {
            let microphone = FirstOrder::new(alpha, [0.0, 0.0, -1.0]).unwrap();
            let mut point = BaffledPressure::new_multirate_first_order(
                &surface,
                1,
                48_000 * ratio,
                48_000,
                position,
                medium(),
                microphone,
            )
            .unwrap();
            let mut area = BaffledPressure::new_multirate_aperture_first_order(
                &surface,
                1,
                48_000 * ratio,
                48_000,
                &aperture,
                medium(),
                microphone,
                Options::default(),
                || false,
            )
            .unwrap();
            assert_eq!(point.kernels, area.kernels);
            assert_eq!(point.filter_delay_samples(), area.filter_delay_samples());
            assert_eq!(point.delay_samples, area.delay_samples);
            for frame in 0..300 {
                let trace: Vec<_> = (0..ratio)
                    .map(|s| ((frame * ratio + s) as f64 * 0.13).sin())
                    .collect();
                assert_eq!(
                    point.step_trace(&trace).unwrap(),
                    area.step_trace(&trace).unwrap()
                );
            }
        }
    }
}

#[test]
fn g1_circular_capsule_pressure_and_velocity_converge_to_exact_area_integrals() {
    // Independent disk integrals of a baffled point volume-velocity source:
    // p_bar=rho*c*Q/(pi*a²)(exp(ikz)-exp(ikR));
    // uz_bar=Q/(pi*a²)(exp(ikz)-z/R*exp(ikR)), R=sqrt(z²+a²).
    let (radius, height, hz) = (0.03_f64, 0.08_f64, 6000.0);
    let k = TAU * hz / medium().sound_speed;
    let edge = height.hypot(radius);
    let q = 0.01;
    let expected_p = (phase(k * height) - phase(k * edge))
        .scale(medium().density * medium().sound_speed * q / (PI * radius * radius));
    let expected_u = (phase(k * height) - phase(k * edge).scale(height / edge))
        .scale(q / (PI * radius * radius));
    let mut previous = f64::INFINITY;
    for rings in [2, 8, 32] {
        let aperture =
            Aperture::circular([0.0, 0.0, height], [0.0, 0.0, 1.0], radius, rings, 16).unwrap();
        let prepared = harmonic(&[source(0.0, 1.0)], &aperture, hz);
        assert_eq!(prepared.receiver_count(), 1);
        let field = prepared.transfer(0, 0).unwrap();
        let error = ((field.pressure_pa - expected_p).abs() / expected_p.abs())
            .max((field.particle_velocity_m_s[2] - expected_u).abs() / expected_u.abs());
        assert!(
            error < 0.08 * previous,
            "rings={rings}, error={error}, previous={previous}"
        );
        assert!(field.particle_velocity_m_s[0].abs() < 1e-14);
        assert!(field.particle_velocity_m_s[1].abs() < 1e-14);
        close(
            field.pressure_gradient_pa_m[2],
            field.particle_velocity_m_s[2] * C64::new(0.0, TAU * hz * medium().density),
            1e-15,
        );
        previous = error;
    }
    assert!(previous < 2e-5, "finest disk error {previous}");
}

#[test]
fn g3_area_pressure_null_retains_particle_velocity_and_a_weak_area_residual() {
    let surface = [source(-0.1, 1.0), source(0.1, -1.0)];
    let left = [-0.07, 0.0, 0.12];
    let right = [0.07, 0.0, 0.12];
    let aperture = Aperture::new(vec![
        Sample {
            position_m: left,
            area_m2: 1.0,
        },
        Sample {
            position_m: right,
            area_m2: 1.0,
        },
    ])
    .unwrap();
    let prepared = harmonic(&surface, &aperture, 1300.0);
    let field = prepared.transfer(0, 0).unwrap();
    assert_eq!(field.pressure_pa, C64::ZERO);
    assert!(field.particle_velocity_m_s[0].abs() > 0.1);
    let pattern = FirstOrder::new(0.0, [1.0, 0.0, 0.0]).unwrap();
    let mut output = [C64::ZERO];
    prepared
        .observe_into(&[C64::ONE], &[pattern], &mut output)
        .unwrap();
    close(
        output[0],
        field.particle_velocity_m_s[0].scale(-medium().density * medium().sound_speed),
        1e-15,
    );
    let mut silent = BaffledPressure::new_aperture_first_order(
        &surface,
        1,
        48_000,
        &aperture,
        medium(),
        FirstOrder::default(),
        Options::default(),
        || false,
    )
    .unwrap();
    let point = BaffledPressure::new_first_order(
        &surface,
        1,
        48_000,
        left,
        medium(),
        FirstOrder::default(),
    )
    .unwrap();
    let kernel_scale = point
        .kernels
        .iter()
        .map(|(_, w)| w[0].abs())
        .fold(0.0_f64, f64::max);
    let residual = silent
        .kernels
        .iter()
        .map(|(_, w)| w[0].abs())
        .fold(0.0_f64, f64::max);
    // Four source/capsule products can cancel after an intermediate DD sum
    // has rounded. Bound that numerical residue, never add a noise floor or
    // replace a small physical coefficient by zero.
    assert!(
        residual < 1e-29 * kernel_scale,
        "relative dark residual {}",
        residual / kernel_scale
    );
    for frame in 0..400 {
        let pressure = silent.step(&[(frame as f64 * 0.17).sin()]).unwrap();
        assert!(pressure.abs() < 1e-29 * kernel_scale * 48_000.0 * silent.kernels.len() as f64);
    }

    // The middle area is below f64's unit-relative precision. Rounding the
    // large weighted fields before cancellation would lose this physical patch.
    let weak = 2.0_f64.powi(-60);
    let aperture = Aperture::new(vec![
        Sample {
            position_m: left,
            area_m2: 1.0,
        },
        Sample {
            position_m: left,
            area_m2: weak,
        },
        Sample {
            position_m: right,
            area_m2: 1.0,
        },
    ])
    .unwrap();
    let point = Prepared::new(
        &surface,
        1,
        &[left],
        1300.0,
        medium(),
        HarmonicOptions::default(),
        || false,
    )
    .unwrap();
    let expected = point.transfer(0, 0).unwrap().pressure_pa.scale(weak / 2.0);
    close(
        harmonic(&surface, &aperture, 1300.0)
            .transfer(0, 0)
            .unwrap()
            .pressure_pa,
        expected,
        2e-13,
    );
    let area = BaffledPressure::new_aperture_first_order(
        &surface,
        1,
        48_000,
        &aperture,
        medium(),
        FirstOrder::default(),
        Options::default(),
        || false,
    )
    .unwrap();
    let point = BaffledPressure::new_first_order(
        &surface,
        1,
        48_000,
        left,
        medium(),
        FirstOrder::default(),
    )
    .unwrap();
    assert_eq!(area.kernels.len(), point.kernels.len());
    let scale = point
        .kernels
        .iter()
        .map(|(_, w)| w[0].abs())
        .fold(0.0_f64, f64::max)
        * weak
        / 2.0;
    for ((da, a), (dp, p)) in area.kernels.iter().zip(&point.kernels) {
        assert_eq!(da, dp);
        assert!((a[0] - p[0] * weak / 2.0).abs() < 3e-13 * scale);
    }
}

#[test]
fn g2_streamed_multirate_capsule_matches_its_continuous_vector_field() {
    let surface = [source(-0.04, 1.0), source(0.07, -0.2)];
    let aperture = Aperture::circular([0.03, 0.0, 0.12], [0.2, 0.1, 1.0], 0.018, 4, 16).unwrap();
    let pattern = FirstOrder::new(0.35, [0.0, 0.0, -1.0]).unwrap();
    let hz = 6000.0;
    let mut exact = [C64::ZERO];
    harmonic(&surface, &aperture, hz)
        .observe_into(&[C64::from_re(0.001)], &[pattern], &mut exact)
        .unwrap();
    for ratio in [1, 3, 4] {
        let mut receiver = BaffledPressure::new_multirate_aperture_first_order(
            &surface,
            1,
            48_000 * ratio,
            48_000,
            &aperture,
            medium(),
            pattern,
            Options::default(),
            || false,
        )
        .unwrap();
        let latency = receiver.filter_delay_samples();
        let settle = receiver.delay_samples.1 + (2.0 * latency).ceil() as usize + 16;
        let storage = (receiver.history.as_ptr(), receiver.kernels.as_ptr());
        let omega = TAU * hz / 48_000.0;
        let (mut error, mut power) = (0.0, 0.0);
        for frame in 0..settle + 2048 {
            let trace: Vec<_> = (0..ratio)
                .map(|sub| 0.001 * (omega * (frame as f64 + sub as f64 / ratio as f64)).cos())
                .collect();
            let p = receiver.step_trace(&trace).unwrap();
            if frame >= settle {
                let angle = omega * (frame as f64 + (ratio - 1) as f64 / ratio as f64 - latency);
                let expected = exact[0].re * angle.cos() + exact[0].im * angle.sin();
                error += (p - expected).powi(2);
                power += expected.powi(2);
            }
        }
        assert!(
            (error / power).sqrt() < 5e-4,
            "ratio {ratio}: {}",
            (error / power).sqrt()
        );
        assert_eq!(
            storage,
            (receiver.history.as_ptr(), receiver.kernels.as_ptr())
        );
        assert!(
            receiver
                .directional_velocity
                .as_ref()
                .unwrap()
                .directional_velocity
                .is_none()
        );
        let mut clone = receiver.clone();
        assert!(
            receiver
                .step_trace(&vec![f64::NAN; ratio as usize])
                .is_err()
        );
        for _ in 0..10 {
            assert_eq!(
                receiver.step_trace(&vec![0.0; ratio as usize]).unwrap(),
                clone.step_trace(&vec![0.0; ratio as usize]).unwrap()
            );
        }
    }
}

#[test]
fn g4_complete_work_clearance_and_cancellation_refuse_without_a_partial_receiver() {
    assert!(Aperture::new(vec![]).is_err());
    for area in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            Aperture::new(vec![Sample {
                position_m: [0.0, 0.0, 0.1],
                area_m2: area
            }])
            .is_err()
        );
    }
    assert!(Aperture::circular([0.0; 3], [0.0; 3], 0.01, 4, 8).is_err());
    assert!(Aperture::circular([0.0; 3], [0.0, 0.0, 1.0], 0.01, 512, 1024).is_err());
    let surface = [source(0.0, 1.0)];
    let aperture = Aperture::circular([0.0, 0.0, 0.1], [0.0, 0.0, 1.0], 0.01, 2, 8).unwrap();
    let pattern = FirstOrder::default();
    let required = (aperture.samples().len() + 1) * (path::TAPS + 2);
    assert!(
        BaffledPressure::new_aperture_first_order(
            &surface,
            1,
            48_000,
            &aperture,
            medium(),
            pattern,
            Options {
                maximum_kernel_terms: required - 1
            },
            || false
        )
        .is_err()
    );
    assert!(
        BaffledPressure::new_aperture_first_order(
            &surface,
            1,
            48_000,
            &aperture,
            medium(),
            pattern,
            Options {
                maximum_kernel_terms: required
            },
            || false
        )
        .is_ok()
    );
    assert!(
        Prepared::new_aperture(
            &surface,
            1,
            &aperture,
            1000.0,
            medium(),
            HarmonicOptions {
                maximum_kernel_terms: aperture.samples().len() - 1,
                ..HarmonicOptions::default()
            },
            || false
        )
        .is_err()
    );
    for after in [0, 2, 5] {
        let mut calls = 0;
        assert!(
            BaffledPressure::new_aperture_first_order(
                &surface,
                1,
                48_000,
                &aperture,
                medium(),
                pattern,
                Options::default(),
                || {
                    calls += 1;
                    calls > after
                }
            )
            .is_err()
        );
        let mut calls = 0;
        assert!(
            Prepared::new_aperture(
                &surface,
                1,
                &aperture,
                1000.0,
                medium(),
                HarmonicOptions::default(),
                || {
                    calls += 1;
                    calls > after
                }
            )
            .is_err()
        );
    }
    let low = Aperture::new(vec![
        Sample {
            position_m: [0.0, 0.0, 0.1],
            area_m2: 1.0,
        },
        Sample {
            position_m: [0.0, 0.0, 0.049],
            area_m2: 1.0,
        },
    ])
    .unwrap();
    assert!(
        BaffledPressure::new_aperture_first_order(
            &surface,
            1,
            48_000,
            &low,
            medium(),
            pattern,
            Options::default(),
            || false
        )
        .is_err()
    );
    assert!(
        Prepared::new_aperture(
            &surface,
            1,
            &low,
            1000.0,
            medium(),
            HarmonicOptions::default(),
            || false
        )
        .is_err()
    );
}
