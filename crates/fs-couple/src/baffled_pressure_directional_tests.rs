use super::super::harmonic::{Options, Prepared};
use super::*;
use fs_math::c64::C64;

fn medium() -> RayleighMedium {
    RayleighMedium {
        density: 1.2,
        sound_speed: 320.0,
    }
}
fn patch(x: f64, shapes: &[f64]) -> SurfaceSample {
    SurfaceSample {
        position_m: [x, 0.0, 0.0],
        area_m2: 0.01,
        mode_shape: shapes.to_vec(),
    }
}
fn phasor(phase: f64) -> C64 {
    C64::new(phase.cos(), phase.sin())
}

// An actual streamed waveform against a continuous-frequency reference, after
// finite startup support. The last mechanics substep determines the sample time.
fn waveform_error(mut mic: BaffledPressure, source: &[C64], hz: f64, exact: C64) -> f64 {
    let ratio = mic.decimator.input_frames();
    let omega = TAU * hz / mic.rate;
    let latency = mic.filter_delay_samples();
    let settle = mic.delay_samples.1 + (2.0 * latency).ceil() as usize + 16;
    let mut trace = vec![0.0; ratio * source.len()];
    let (mut error, mut signal) = (0.0, 0.0);
    for frame in 0..settle + 2048 {
        for substep in 0..ratio {
            let phase = omega * (frame as f64 + substep as f64 / ratio as f64);
            for (i, value) in source.iter().enumerate() {
                trace[substep * source.len() + i] = value.re * phase.cos() + value.im * phase.sin();
            }
        }
        let pressure = mic.step_trace(&trace).unwrap();
        if frame >= settle {
            let phase = omega * (frame as f64 + (ratio - 1) as f64 / ratio as f64 - latency);
            let reference = exact.re * phase.cos() + exact.im * phase.sin();
            error += (pressure - reference).powi(2);
            signal += reference.powi(2);
        }
    }
    assert!(signal > 0.0);
    (error / signal).sqrt()
}

#[test]
fn g2_streamed_cardioid_proximity_and_pattern_follow_the_point_field() {
    let (rate, hz, source) = (32_000, 50.0, C64::from_re(0.001));
    let gas = medium();
    let surface = [patch(0.0, &[1.0])];
    for kr in [0.05, 0.2, 1.0, 10.0, 50.0] {
        let radius = kr * gas.sound_speed / (TAU * hz);
        let pressure =
            C64::new(0.0, -gas.density * hz * surface[0].area_m2 / radius) * phasor(kr) * source;
        for (alpha, axis) in [
            (1.0, [0.0, 0.0, -1.0]),
            (0.5, [0.0, 0.0, -1.0]),
            (0.0, [0.0, 0.0, -1.0]),
            (0.5, [0.0, 0.0, 1.0]),
            (0.5, [1.0, 0.0, 0.0]),
        ] {
            // Euler momentum for exp(-iwt) gives u_r=p/(rho*c)*(1+i/(kR)).
            // The source-facing axis is -Rhat. No distance-gain filter is added.
            let exact = pressure
                * (C64::from_re(alpha) - C64::new(1.0, 1.0 / kr).scale((1.0 - alpha) * axis[2]));
            let microphone = FirstOrder::new(alpha, axis).unwrap();
            let receiver = BaffledPressure::new_first_order(
                &surface,
                1,
                rate,
                [0.0, 0.0, radius],
                gas,
                microphone,
            )
            .unwrap();
            let error = waveform_error(receiver, &[source], hz, exact);
            assert!(
                error < 5e-4,
                "kR={kr}, alpha={alpha}, axis={axis:?}: {error}"
            );
        }
    }
}

#[test]
fn g3_multirate_surface_projection_matches_the_harmonic_owner() {
    let surface = [patch(-0.2, &[1.0, -0.4]), patch(0.13, &[0.3, 0.7])];
    let position = [0.1, -0.2, 0.7];
    let source = [C64::new(0.0004, 0.0003), C64::new(-0.0002, 0.0001)];
    let microphone = FirstOrder::new(0.3, [0.0, 0.0, -1.0]).unwrap();
    for hz in [200.0, 6000.0, 21_600.0] {
        let field = Prepared::new(
            &surface,
            2,
            &[position],
            hz,
            medium(),
            Options::default(),
            || false,
        )
        .unwrap();
        let mut exact = [C64::ZERO];
        field
            .observe_into(&source, &[microphone], &mut exact)
            .unwrap();
        // Ratio 3 exercises fractional output-frame alignment; 4 exercises the
        // sparse half-band cascade; 1 exercises the 48-frame differentiator.
        for ratio in [1, 3, 4] {
            let receiver = BaffledPressure::new_multirate_first_order(
                &surface,
                2,
                48_000 * ratio,
                48_000,
                position,
                medium(),
                microphone,
            )
            .unwrap();
            let reference = BaffledPressure::new_multirate_differentiating(
                &surface,
                2,
                48_000 * ratio,
                48_000,
                position,
                medium(),
            )
            .unwrap();
            assert_eq!(
                receiver.filter_delay_samples(),
                reference.filter_delay_samples()
            );
            let error = waveform_error(receiver, &source, hz, exact[0]);
            assert!(error < 1e-3, "ratio={ratio}, hz={hz}: {error}");
        }
    }
}

#[test]
fn g3_directional_signal_survives_a_scalar_pressure_null() {
    let surface = [patch(-0.1, &[1.0]), patch(0.1, &[-1.0])];
    let position = [0.0, 0.0, 0.2];
    let microphone = FirstOrder::new(0.5, [1.0, 0.0, 0.0]).unwrap();
    let mut omni =
        BaffledPressure::new_differentiating(&surface, 1, 48_000, position, medium()).unwrap();
    let field = Prepared::new(
        &surface,
        1,
        &[position],
        800.0,
        medium(),
        Options::default(),
        || false,
    )
    .unwrap();
    let source = [C64::from_re(0.001)];
    let mut exact = [C64::ZERO];
    field
        .observe_into(&source, &[microphone], &mut exact)
        .unwrap();
    assert!(exact[0].abs() > 1e-3);
    for frame in 0..1000 {
        assert_eq!(
            omni.step(&[(TAU * 800.0 * frame as f64 / 48_000.0).cos() * 0.001])
                .unwrap(),
            0.0
        );
    }
    let directional =
        BaffledPressure::new_first_order(&surface, 1, 48_000, position, medium(), microphone)
            .unwrap();
    assert!(waveform_error(directional, &source, 800.0, exact[0]) < 5e-4);
}

#[test]
fn g1_circular_cardioid_converges_to_the_analytic_piston() {
    let (radius, height, hz) = (0.08_f64, 0.05_f64, 1000.0);
    let gas = medium();
    let source = C64::from_re(1e-5); // Actual outward volume flow, not surface velocity.
    let velocity = source.scale(1.0 / (std::f64::consts::PI * radius * radius));
    let edge_range = radius.hypot(height);
    let k = TAU * hz / gas.sound_speed;
    // On axis, direct integration gives p=rho*c*v*(exp(ikz)-exp(ikL))
    // and u_z=v*(exp(ikz)-(z/L)*exp(ikL)); project the analytic vector.
    let exact = (phasor(k * height)
        - phasor(k * edge_range).scale(0.5 + 0.5 * height / edge_range))
        * velocity.scale(gas.density * gas.sound_speed);
    let microphone = FirstOrder::new(0.5, [0.0, 0.0, -1.0]).unwrap();
    let errors: Vec<_> = [8, 32, 128]
        .into_iter()
        .map(|radial_rings| {
            let spec = CircularOutletReceiver {
                position_m: [0.0, 0.0, height],
                radial_rings,
                angular_points: 64,
                maximum_frequency_hz: hz,
            };
            let receiver = BaffledPressure::circular_outlet_first_order(
                radius, 48_000, 48_000, spec, gas, microphone,
            )
            .unwrap();
            waveform_error(receiver, &[source], hz, exact)
        })
        .collect();
    assert!(errors[1] < 0.4 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.7 * errors[1], "{errors:?}");
    assert!(errors[2] < 3e-4, "{errors:?}");
}

#[test]
fn g5_omni_keeps_the_original_arithmetic_and_single_bank() {
    let surface = [patch(0.0, &[1.0, -0.3])];
    for ratio in [1, 3, 4] {
        let mut original = BaffledPressure::new_multirate_differentiating(
            &surface,
            2,
            48_000 * ratio,
            48_000,
            [0.0, 0.0, 0.4],
            medium(),
        )
        .unwrap();
        let mut omni = BaffledPressure::new_multirate_first_order(
            &surface,
            2,
            48_000 * ratio,
            48_000,
            [0.0, 0.0, 0.4],
            medium(),
            FirstOrder::default(),
        )
        .unwrap();
        assert!(omni.directional_velocity.is_none());
        assert_eq!(
            omni.multiply_adds_per_sample(),
            original.multiply_adds_per_sample()
        );
        for frame in 0..400 {
            let trace: Vec<_> = (0..2 * ratio)
                .map(|i| (frame as f64 * 0.1 + f64::from(i) * 0.3).sin())
                .collect();
            assert_eq!(
                omni.step_trace(&trace).unwrap().to_bits(),
                original.step_trace(&trace).unwrap().to_bits()
            );
        }
    }
}

#[test]
fn g4_shared_histories_clone_reject_and_decompose_without_advancing() {
    let surface = [patch(-0.1, &[1.0, -0.2]), patch(0.1, &[0.7, 0.3])];
    let microphone = FirstOrder::new(0.5, [0.0, 0.0, -1.0]).unwrap();
    let mut receiver = BaffledPressure::new_multirate_first_order(
        &surface,
        2,
        144_000,
        48_000,
        [0.1, 0.0, 0.4],
        medium(),
        microphone,
    )
    .unwrap();
    let addresses = (
        receiver.history.as_ptr(),
        receiver.taps.as_ptr(),
        receiver
            .directional_velocity
            .as_ref()
            .unwrap()
            .history
            .as_ptr(),
        receiver
            .directional_velocity
            .as_ref()
            .unwrap()
            .taps
            .as_ptr(),
    );
    let mut control = receiver.clone();
    let mut components = [0.0; 2];
    for frame in 0..600 {
        assert!(receiver.step_trace(&[f64::NAN; 6]).is_err());
        assert!(receiver.step_trace(&[0.0; 5]).is_err());
        assert!(receiver.step_trace(&[f64::INFINITY; 6]).is_err());
        let trace: [f64; 6] =
            std::array::from_fn(|i| (frame as f64 * 0.07 + i as f64 * 0.31).sin());
        let value = receiver.step_trace(&trace).unwrap();
        assert_eq!(
            value.to_bits(),
            control.step_trace(&trace).unwrap().to_bits()
        );
        receiver.mode_pressures(&mut components).unwrap();
        assert!((components.iter().sum::<f64>() - value).abs() < 1e-9);
        if frame < receiver.delay_samples.0 {
            assert_eq!(value, 0.0);
        }
        if frame == 300 {
            control = receiver.clone();
        }
    }
    assert_eq!(
        addresses,
        (
            receiver.history.as_ptr(),
            receiver.taps.as_ptr(),
            receiver
                .directional_velocity
                .as_ref()
                .unwrap()
                .history
                .as_ptr(),
            receiver
                .directional_velocity
                .as_ref()
                .unwrap()
                .taps
                .as_ptr()
        )
    );
}

#[test]
fn g0_alignment_and_circular_band_respect_original_budgets() {
    let surface = [patch(0.0, &[1.0])];
    let position = [
        0.0,
        0.0,
        (MAX_HISTORY - TAPS - 20) as f64 * 320.0 / 32_000.0,
    ];
    assert!(BaffledPressure::new_differentiating(&surface, 1, 32_000, position, medium()).is_ok());
    assert!(
        BaffledPressure::new_first_order(
            &surface,
            1,
            32_000,
            position,
            medium(),
            FirstOrder::new(0.5, [0.0, 0.0, -1.0]).unwrap()
        )
        .is_err()
    );
    let spec = CircularOutletReceiver {
        position_m: [0.0, 0.0, 0.2],
        radial_rings: 8,
        angular_points: 32,
        maximum_frequency_hz: 21_600.0,
    };
    assert!(
        BaffledPressure::circular_outlet_first_order(
            0.001,
            144_000,
            48_000,
            spec,
            medium(),
            FirstOrder::default()
        )
        .is_ok()
    );
    assert!(
        BaffledPressure::circular_outlet_first_order(
            0.001,
            144_000,
            48_000,
            CircularOutletReceiver {
                maximum_frequency_hz: 21_601.0,
                ..spec
            },
            medium(),
            FirstOrder::default()
        )
        .is_err()
    );
}
