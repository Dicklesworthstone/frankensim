//! G1/G3: original liquid coefficients reach every signed physical sound path.
use super::super::path;
use super::*;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
use fs_material::fluid_transport::FluidTransportState;
use std::f64::consts::TAU;

fn propagation(bulk: f64) -> FluidPropagation {
    // Authored viscous-liquid fixture, no assertion of measured properties.
    let state = ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(1000.0, 1200.0).unwrap(),
        FluidTransportState::try_new(1000.0, 4000.0, 0.6, 5.0).unwrap(),
        3900.0,
    )
    .unwrap();
    FluidPropagation {
        absorption: ClassicalBulkAbsorption::from_material(&state, bulk).unwrap(),
        maximum_truncated_energy: 1e-12,
        limits: FluidPropagationLimits::surface_reference(),
    }
}

fn surface() -> Vec<SurfaceSample> {
    vec![
        SurfaceSample {
            position_m: [-0.3, 0.0, 0.0],
            area_m2: 0.2,
            mode_shape: vec![1.0, 0.2],
        },
        SurfaceSample {
            position_m: [0.3, 0.0, 0.0],
            area_m2: 0.1,
            mode_shape: vec![-0.7, 0.4],
        },
        SurfaceSample {
            position_m: [0.8, 0.4, 0.0],
            area_m2: 0.1,
            mode_shape: vec![0.3, -0.5],
        },
    ]
}

#[test]
fn g1_signed_fluid_paths_keep_distinct_attenuation_and_relative_phase() {
    let propagation = propagation(10.0);
    let surface = surface();
    let position = [0.0, 0.0, 2.0];
    let rate = 48_000;
    let (receiver, report) = BaffledPressure::new_multirate_fluid_differentiating(
        &surface,
        2,
        rate,
        rate,
        position,
        &propagation,
        &mut || false,
    )
    .unwrap();
    assert_eq!(report.surface_filter_indices, [0, 0, 1]);
    assert_ne!(report.filters[0].range_m, report.filters[1].range_m);
    for frequency in [100.0, 4000.0, 12_000.0, 21_000.0] {
        let omega_dt = TAU * frequency / f64::from(rate);
        for input in 0..2 {
            let (mut expected_re, mut expected_im, mut scale) = (0.0, 0.0, 0.0);
            for (sample, &index) in surface.iter().zip(&report.surface_filter_indices) {
                let receipt = report.filters[index];
                assert_eq!(receipt.absorption, propagation.absorption);
                let filter = FluidObserverFilter::new(
                    propagation.absorption,
                    receipt.range_m,
                    1.0 / f64::from(rate),
                    propagation.maximum_truncated_energy,
                )
                .unwrap();
                let (re, im) =
                    filter
                        .taps()
                        .iter()
                        .enumerate()
                        .fold((0.0, 0.0), |(re, im), (n, h)| {
                            (
                                re + h * (omega_dt * n as f64).cos(),
                                im - h * (omega_dt * n as f64).sin(),
                            )
                        });
                let phase = -omega_dt
                    * (receipt.receiver_delay_s * f64::from(rate)
                        + report.additional_latency_samples as f64
                        + path::HALF as f64);
                let gain =
                    1000.0 * sample.area_m2 * sample.mode_shape[input] / (TAU * receipt.range_m);
                expected_re += gain * (re * phase.cos() - im * phase.sin());
                expected_im += gain * (re * phase.sin() + im * phase.cos());
                scale += gain.abs() * re.hypot(im);
                assert!(
                    (receipt.receiver_delay_s + receipt.low_frequency_delay_s
                        - receipt.range_m / 1200.0)
                        .abs()
                        < 1e-17
                );
            }
            let (re, im) = receiver
                .kernels
                .iter()
                .fold((0.0, 0.0), |(re, im), (n, weights)| {
                    (
                        re + weights[input] * (omega_dt * *n as f64).cos(),
                        im - weights[input] * (omega_dt * *n as f64).sin(),
                    )
                });
            assert!((re - expected_re).hypot(im - expected_im) < 5e-4 * scale);
        }
    }
}

#[test]
fn g3_actual_liquid_pressure_has_the_independent_bulk_viscosity_attenuation_ratio() {
    let patch = [SurfaceSample {
        position_m: [0.0; 3],
        area_m2: 0.2,
        mode_shape: vec![1.0],
    }];
    let mut receivers: Vec<_> = [0.0, 10.0]
        .into_iter()
        .map(|bulk| {
            BaffledPressure::new_multirate_fluid_differentiating(
                &patch,
                1,
                48_000,
                48_000,
                [0.0, 0.0, 2.0],
                &propagation(bulk),
                &mut || false,
            )
            .unwrap()
            .0
        })
        .collect();
    let mut energy = [0.0; 2];
    for frame in 0..4096 {
        let velocity = [(TAU * 12_000.0 * frame as f64 / 48_000.0).sin()];
        for (i, receiver) in receivers.iter_mut().enumerate() {
            let pressure = receiver.step(&velocity).unwrap();
            if frame >= 2048 {
                energy[i] += pressure * pressure;
            }
        }
    }
    let omega = TAU * 12_000.0;
    let expected = (-omega * omega * 10.0 * 2.0 / (2.0 * 1000.0 * 1200.0_f64.powi(3))).exp();
    assert!(energy[0] > 1.0 && energy[1] < 0.95 * energy[0]);
    assert!(((energy[1] / energy[0]).sqrt() - expected).abs() < 2e-3);
    let mut replay = receivers[1].clone();
    for value in [0.2, -0.4, 0.0, 0.1] {
        assert_eq!(
            receivers[1].step(&[value]).unwrap().to_bits(),
            replay.step(&[value]).unwrap().to_bits()
        );
    }
}

#[test]
fn g0_fluid_surface_keeps_original_path_work_cancellation_and_applicability_limits() {
    let original = propagation(10.0);
    let surface = surface();
    let prepare = |p: &FluidPropagation, cancelled: &mut dyn FnMut() -> bool| {
        BaffledPressure::new_multirate_fluid_differentiating(
            &surface,
            2,
            48_000,
            48_000,
            [0.0, 0.0, 2.0],
            p,
            &mut || cancelled(),
        )
    };
    assert!(prepare(&original, &mut || true).is_err());
    let mut limited = original;
    limited.limits.max_unique_paths = 1;
    assert!(prepare(&limited, &mut || false).is_err());
    limited = original;
    limited.limits.max_filter_taps = 1;
    assert!(prepare(&limited, &mut || false).is_err());
    limited = original;
    limited.limits.max_kernel_terms = 1;
    assert!(prepare(&limited, &mut || false).is_err());
    limited = original;
    limited.maximum_truncated_energy = f64::NAN;
    assert!(prepare(&limited, &mut || false).is_err());
    assert!(prepare(&propagation(1e6), &mut || false).is_err());
    assert!(prepare(&original, &mut || false).is_ok());
}

fn kernel(receiver: &BaffledPressure, input: usize, frequency_hz: f64) -> [f64; 2] {
    let omega = TAU * frequency_hz / receiver.rate;
    receiver
        .kernels
        .iter()
        .fold([0.0; 2], |[re, im], (n, weights)| {
            [
                re + weights[input] * (omega * *n as f64).cos(),
                im - weights[input] * (omega * *n as f64).sin(),
            ]
        })
}

#[test]
fn g3_fluid_circular_flow_keeps_the_independent_piston_directivity() {
    use super::super::{CircularOutletReceiver, openings::CircularOpening};
    // Explicit nearly lossless liquid fixture. The comparison is the analytic
    // far-field uniform-piston ratio 2 J1(ka sin(theta))/(ka sin(theta)), not
    // another evaluation of the surface quadrature. Bulk viscosity is supplied.
    let state = ThermoviscousFluidState::try_new(
        FluidAcousticState::try_new(1000.0, 1200.0).unwrap(),
        FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.001).unwrap(),
        3900.0,
    )
    .unwrap();
    let selection = FluidPropagation {
        absorption: ClassicalBulkAbsorption::from_material(&state, 0.0).unwrap(),
        ..propagation(0.0)
    };
    let (radius, range, frequency, theta) = (0.015, 2.0, 18_000.0, 1.0_f64);
    let opening = [CircularOpening {
        center_m: [0.0; 2],
        radius_m: radius,
    }];
    let prepare = |position_m| {
        BaffledPressure::circular_openings_multirate_fluid_differentiating(
            &opening,
            48_000,
            48_000,
            CircularOutletReceiver {
                position_m,
                radial_rings: 8,
                angular_points: 32,
                maximum_frequency_hz: frequency,
            },
            &selection,
            &mut || false,
        )
        .unwrap()
        .0
    };
    let axial = kernel(&prepare([0.0, 0.0, range]), 0, frequency);
    let oblique = kernel(
        &prepare([range * theta.sin(), 0.0, range * theta.cos()]),
        0,
        frequency,
    );
    let argument = TAU * frequency * radius * theta.sin() / 1200.0;
    let expected = 2.0 * fs_math::special::bessel_j1(argument) / argument;
    let actual = oblique[0].hypot(oblique[1]) / axial[0].hypot(axial[1]);
    assert!((actual - expected).abs() < 3e-3, "{actual} != {expected}");
    assert!(
        actual < 0.9,
        "the physical aperture must not become a point source"
    );
}

#[test]
fn g1_aligned_fluid_openings_reconstruct_one_signed_surface_response() {
    use super::super::{CircularOutletReceiver, openings::CircularOpening};
    let openings = [
        CircularOpening {
            center_m: [0.0; 2],
            radius_m: 0.001,
        },
        CircularOpening {
            center_m: [2.0, 0.0],
            radius_m: 0.001,
        },
    ];
    let spec = CircularOutletReceiver {
        position_m: [0.0, 0.0, 0.05],
        radial_rings: 1,
        angular_points: 8,
        maximum_frequency_hz: 21_000.0,
    };
    let selection = propagation(10.0);
    let prepare = |openings: &[CircularOpening]| {
        BaffledPressure::circular_openings_multirate_fluid_differentiating(
            openings,
            48_000,
            48_000,
            spec,
            &selection,
            &mut || false,
        )
        .unwrap()
    };
    let (mut joint, mut joint_report) = prepare(&openings);
    let mut split: Vec<_> = openings
        .iter()
        .map(|opening| prepare(core::slice::from_ref(opening)))
        .collect();
    let target = split
        .iter()
        .map(|(_, report)| report.additional_latency_samples)
        .chain([joint_report.additional_latency_samples])
        .max()
        .unwrap()
        + 7;
    for (receiver, report) in &mut split {
        let before = receiver.filter_delay_samples();
        let extra = target - report.additional_latency_samples;
        let filters = report.filters.clone();
        receiver.align_fluid_latency(report, target).unwrap();
        assert_eq!(receiver.filter_delay_samples(), before + extra as f64);
        assert_eq!(report.additional_latency_samples, target);
        assert_eq!(
            report.filters, filters,
            "physical path receipts must not change"
        );
    }
    joint
        .align_fluid_latency(&mut joint_report, target)
        .unwrap();
    assert_eq!(joint_report.surface_filter_indices.len(), 16);
    assert_eq!(
        split[0].0.filter_delay_samples(),
        split[1].0.filter_delay_samples()
    );
    for frequency in [100.0, 4000.0, 12_000.0, 21_000.0] {
        let a = kernel(&split[0].0, 0, frequency);
        let b = kernel(&split[1].0, 0, frequency);
        let x = kernel(&joint, 0, frequency);
        let y = kernel(&joint, 1, frequency);
        let scale = a[0].hypot(a[1]) + 0.4 * b[0].hypot(b[1]);
        let error =
            (x[0] - 0.4 * y[0] - a[0] + 0.4 * b[0]).hypot(x[1] - 0.4 * y[1] - a[1] + 0.4 * b[1]);
        assert!(
            error < 3e-12 * scale,
            "signed coherent pressure at {frequency} Hz: {error}"
        );
    }
    // A native volume-flow history reaches the same pressure, not just the
    // static transfer: opposite signed flows remain coherent after alignment.
    for frame in 0..1024 {
        let flow = 1e-6 * (TAU * 1500.0 * frame as f64 / 48_000.0).sin();
        let expected = split[0].0.step(&[flow]).unwrap() - 0.4 * split[1].0.step(&[flow]).unwrap();
        let actual = joint.step(&[flow, -0.4 * flow]).unwrap();
        assert!((actual - expected).abs() < 2e-10 * expected.abs().max(1.0));
    }
}

#[test]
fn g0_fluid_opening_and_latency_refusals_preserve_the_original_clock_and_history() {
    use super::super::{CircularOutletReceiver, openings::CircularOpening};
    let opening = [CircularOpening {
        center_m: [0.0; 2],
        radius_m: 0.001,
    }];
    let spec = CircularOutletReceiver {
        position_m: [0.0, 0.0, 0.05],
        radial_rings: 1,
        angular_points: 8,
        maximum_frequency_hz: 21_000.0,
    };
    let selection = propagation(10.0);
    let prepare =
        |openings: &[CircularOpening], spec, rate, cancelled: &mut dyn FnMut() -> bool| {
            BaffledPressure::circular_openings_multirate_fluid_differentiating(
                openings,
                rate,
                48_000,
                spec,
                &selection,
                &mut || cancelled(),
            )
        };
    assert!(prepare(&[], spec, 192_000, &mut || false).is_err());
    assert!(prepare(&[opening[0], opening[0]], spec, 192_000, &mut || false).is_err());
    assert!(prepare(&opening, spec, 48_001, &mut || false).is_err());
    assert!(
        prepare(
            &opening,
            CircularOutletReceiver {
                maximum_frequency_hz: 21_601.0,
                ..spec
            },
            192_000,
            &mut || false
        )
        .is_err()
    );
    assert!(prepare(&opening, spec, 192_000, &mut || true).is_err());
    let (mut receiver, mut report) = prepare(&opening, spec, 192_000, &mut || false).unwrap();
    let initial = report.additional_latency_samples;
    let initial_kernels = receiver.kernels.clone();
    let initial_delay = receiver.delay_samples;
    let mut stale = report.clone();
    stale.additional_latency_samples += 1;
    assert!(
        receiver
            .align_fluid_latency(&mut stale, initial + 2)
            .is_err()
    );
    assert_eq!(stale.additional_latency_samples, initial + 1);
    assert!(
        receiver
            .align_fluid_latency(&mut report, usize::MAX)
            .is_err()
    );
    assert_eq!(receiver.kernels, initial_kernels);
    assert_eq!(receiver.delay_samples, initial_delay);
    assert_eq!(report.additional_latency_samples, initial);
    receiver
        .align_fluid_latency(&mut report, initial + 7)
        .unwrap();
    let shifted = receiver.kernels.clone();
    assert!(
        receiver
            .align_fluid_latency(&mut report, initial + 6)
            .is_err()
    );
    assert_eq!(receiver.kernels, shifted);
    assert_eq!(report.additional_latency_samples, initial + 7);

    // Rejected mechanical packets do not advance any accepted state. Even a
    // valid all-zero packet must close cold configuration before sound arrives.
    assert!(receiver.step_trace(&[0.0; 3]).is_err());
    assert!(receiver.step_trace(&[0.0, f64::NAN, 0.0, 0.0]).is_err());
    receiver
        .align_fluid_latency(&mut report, initial + 7)
        .unwrap();
    let mut silent = receiver.clone();
    let mut silent_report = report.clone();
    assert_eq!(silent.step_trace(&[0.0; 4]).unwrap(), 0.0);
    assert!(
        silent
            .align_fluid_latency(&mut silent_report, initial + 8)
            .is_err()
    );
    assert_eq!(receiver.step_trace(&[1e-6, 0.0, 0.0, 0.0]).unwrap(), 0.0);
    let mut replay = receiver.clone();
    let mut replay_report = report.clone();
    assert!(
        replay
            .align_fluid_latency(&mut replay_report, initial + 8)
            .is_err()
    );
    assert!(
        receiver
            .align_fluid_latency(&mut report, initial + 8)
            .is_err()
    );
    assert_eq!(report.additional_latency_samples, initial + 7);
    for frame in 0..2 * receiver.history_frames {
        let trace = [if frame % 5 == 0 { 1e-6 } else { 0.0 }; 4];
        assert_eq!(
            receiver.step_trace(&trace).unwrap().to_bits(),
            replay.step_trace(&trace).unwrap().to_bits()
        );
    }
    assert!(
        receiver
            .align_fluid_latency(&mut report, initial + 7)
            .is_err()
    );

    let (mut single_rate, mut single_report) =
        prepare(&opening, spec, 48_000, &mut || false).unwrap();
    single_rate.step(&[1e-6]).unwrap();
    let mut replay = single_rate.clone();
    assert!(single_rate.differentiate_pressure().is_err());
    assert!(
        single_rate
            .align_fluid_latency(&mut single_report, initial + 8)
            .is_err()
    );
    for frame in 0..2 * single_rate.history_frames + 128 {
        let value = [0.0, 2e-6, -1e-6, 0.0][frame % 4];
        assert_eq!(
            single_rate.step(&[value]).unwrap().to_bits(),
            replay.step(&[value]).unwrap().to_bits()
        );
    }

    // A first-order receiver has another physical velocity bank. Padding only
    // its pressure half is unsupported, including when both are still cold.
    let (mut directional, mut directional_report) =
        prepare(&opening, spec, 48_000, &mut || false).unwrap();
    directional.directional_velocity = Some(Box::new(directional.clone()));
    let before = directional.kernels.clone();
    let padding = directional_report.additional_latency_samples;
    assert!(
        directional
            .align_fluid_latency(&mut directional_report, padding + 1)
            .is_err()
    );
    assert_eq!(directional.kernels, before);
    assert_eq!(directional_report.additional_latency_samples, padding);
}
