//! Independent dynamics checks for conservative motion under the held preload.
use fs_nlmodal::{
    PreloadedModalStorage, PreloadedSpectrum, PreloadedSpectrumBudget, PreloadedSpectrumError,
    SosModalStorage, StressChannel, duffing_kinematic_spectrum_exact,
};
use fs_phs::Storage;

fn scalar(omega: f64, beta: f64, reference: f64) -> PreloadedModalStorage {
    PreloadedModalStorage::new(
        SosModalStorage {
            omegas: vec![omega],
            channels: vec![StressChannel {
                coefficient: beta,
                coupling: vec![1.0],
            }],
        },
        &[reference],
    )
    .unwrap()
}

fn budget() -> PreloadedSpectrumBudget {
    PreloadedSpectrumBudget {
        max_intervals: 32768,
        tolerance: 2e-9,
    }
}

fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= scale * tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, scale={scale:.17e}"
    );
}

fn reconstruct(s: &PreloadedSpectrum, phase: f64) -> (f64, f64, f64) {
    let mut q = s.displacement_dc;
    let mut v = 0.0;
    let mut a = 0.0;
    for j in 0..s.cosine_amplitudes.len() {
        let theta = (j + 1) as f64 * phase;
        q += s.cosine_amplitudes[j] * theta.cos();
        v += s.velocity_sine_amplitudes[j] * theta.sin();
        a += s.acceleration_cosine_amplitudes[j] * theta.cos();
    }
    (q, v, a)
}

#[test]
fn biased_orbit_preserves_original_force_energy_and_unequal_turning_points() {
    let storage = scalar(3.0, 20.0, 0.4);
    let spectrum = storage
        .invariant_axis_spectrum(0, 0.6, 24, budget())
        .unwrap();
    assert!(spectrum.displacement_dc < -0.02);
    assert!(spectrum.cosine_amplitudes[1] > 0.01);
    assert!(spectrum.negative_turning_displacement < -0.6);
    assert!(spectrum.intervals >= 4 * 512);
    assert!(spectrum.relative_frequency_change <= budget().tolerance);
    assert!(spectrum.normalized_kinematic_change <= budget().tolerance);
    let energy = storage.hamiltonian(&[0.6, 0.0]);
    for j in 0..=512 {
        let phase = std::f64::consts::TAU * j as f64 / 512.0;
        let (q, v, a) = reconstruct(&spectrum, phase);
        let mut effort = [0.0; 2];
        storage.gradient(&[q, v], &mut effort);
        close(a, -effort[0], 30.0, 3e-8);
        close(storage.hamiltonian(&[q, v]), energy, energy, 3e-8);
    }
    close(reconstruct(&spectrum, 0.0).0, 0.6, 0.6, 2e-9);
    close(
        reconstruct(&spectrum, std::f64::consts::PI).0,
        spectrum.negative_turning_displacement,
        1.0,
        2e-9,
    );
    close(
        storage.hamiltonian(&[spectrum.negative_turning_displacement, 0.0]),
        energy,
        energy,
        3e-15,
    );
    assert!(spectrum.guard.acceleration_l1 < 1e-7);
}

// Independent fixed-step RK4 integrates the ACTUAL original Storage::gradient;
// it does not use the spectrum's scalar reduction or energy quadrature.
fn step(storage: &PreloadedModalStorage, state: [f64; 2], dt: f64) -> [f64; 2] {
    let rhs = |x: [f64; 2]| {
        let mut effort = [0.0; 2];
        storage.gradient(&x, &mut effort);
        [effort[1], -effort[0]]
    };
    let a = rhs(state);
    let b = rhs([state[0] + 0.5 * dt * a[0], state[1] + 0.5 * dt * a[1]]);
    let c = rhs([state[0] + 0.5 * dt * b[0], state[1] + 0.5 * dt * b[1]]);
    let d = rhs([state[0] + dt * c[0], state[1] + dt * c[1]]);
    [
        state[0] + dt / 6.0 * (a[0] + 2.0 * b[0] + 2.0 * c[0] + d[0]),
        state[1] + dt / 6.0 * (a[1] + 2.0 * b[1] + 2.0 * c[1] + d[1]),
    ]
}

#[test]
fn period_dc_even_and_odd_coefficients_match_original_ode_trajectory() {
    let storage = scalar(3.0, 20.0, 0.4);
    let s = storage
        .invariant_axis_spectrum(0, 0.6, 16, budget())
        .unwrap();
    let samples = 16384;
    let dt = std::f64::consts::TAU / (s.angular_frequency * samples as f64);
    let mut state = [0.6, 0.0];
    let mut dc = 0.0;
    let mut q = vec![0.0; 16];
    let mut v = vec![0.0; 16];
    for j in 0..samples {
        let phase = std::f64::consts::TAU * j as f64 / samples as f64;
        dc += state[0] / samples as f64;
        for h in 1..=16 {
            q[h - 1] += 2.0 * state[0] * (h as f64 * phase).cos() / samples as f64;
            v[h - 1] += 2.0 * state[1] * (h as f64 * phase).sin() / samples as f64;
        }
        if j % 64 == 0 {
            let (spectral_q, spectral_v, _) = reconstruct(&s, phase);
            close(spectral_q, state[0], 1.0, 2e-9);
            close(spectral_v, state[1], s.angular_frequency, 3e-9);
        }
        state = step(&storage, state, dt);
    }
    close(state[0], 0.6, 0.6, 2e-10);
    close(state[1], 0.0, s.angular_frequency, 2e-10);
    close(dc, s.displacement_dc, 0.6, 2e-9);
    for j in 0..16 {
        close(q[j], s.cosine_amplitudes[j], 0.6, 2e-9);
        close(
            v[j],
            s.velocity_sine_amplitudes[j],
            s.angular_frequency,
            3e-9,
        );
    }
}

#[test]
fn reversing_preload_reverses_dc_and_even_orders_at_equal_energy() {
    let a = scalar(3.0, 20.0, 0.4)
        .invariant_axis_spectrum(0, 0.6, 20, budget())
        .unwrap();
    let b = scalar(3.0, 20.0, -0.4)
        .invariant_axis_spectrum(0, -a.negative_turning_displacement, 20, budget())
        .unwrap();
    close(
        a.angular_frequency,
        b.angular_frequency,
        a.angular_frequency,
        2e-10,
    );
    close(a.displacement_dc, -b.displacement_dc, 1.0, 2e-9);
    close(b.negative_turning_displacement, -0.6, 1.0, 3e-15);
    for j in 0..20 {
        // q_b(phi)=-q_a(phi+pi): odd terms agree, even terms reverse.
        let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
        close(
            b.cosine_amplitudes[j],
            sign * a.cosine_amplitudes[j],
            1.0,
            2e-9,
        );
        close(
            b.acceleration_cosine_amplitudes[j],
            sign * a.acceleration_cosine_amplitudes[j],
            a.angular_frequency * a.angular_frequency,
            3e-8,
        );
    }
}

#[test]
fn symmetric_and_zero_amplitude_paths_keep_exact_jacobi_owner() {
    let storage = scalar(3.0, 20.0, 0.0);
    let s = storage
        .invariant_axis_spectrum(0, 0.6, 16, budget())
        .unwrap();
    let exact = duffing_kinematic_spectrum_exact(3.0, 20.0, 0.6, 16).unwrap();
    assert_eq!(s.angular_frequency, exact.displacement.angular_frequency);
    assert_eq!(s.intervals, 0);
    assert_eq!(s.displacement_dc, 0.0);
    for j in 0..16 {
        if j % 2 == 0 {
            assert_eq!(
                s.cosine_amplitudes[j],
                exact.displacement.cosine_amplitudes[j / 2]
            );
            assert_eq!(
                s.velocity_sine_amplitudes[j],
                exact.velocity_sine_amplitudes[j / 2]
            );
            assert_eq!(
                s.acceleration_cosine_amplitudes[j],
                exact.acceleration_cosine_amplitudes[j / 2]
            );
        } else {
            assert_eq!(s.cosine_amplitudes[j], 0.0);
            assert_eq!(s.velocity_sine_amplitudes[j], 0.0);
            assert_eq!(s.acceleration_cosine_amplitudes[j], 0.0);
        }
    }
    let zero = scalar(3.0, 20.0, 0.4)
        .invariant_axis_spectrum(0, 0.0, 16, budget())
        .unwrap();
    close(
        zero.angular_frequency,
        18.6_f64.sqrt(),
        18.6_f64.sqrt(),
        3e-15,
    );
    assert!(
        zero.cosine_amplitudes
            .iter()
            .chain(&zero.velocity_sine_amplitudes)
            .chain(&zero.acceleration_cosine_amplitudes)
            .all(|x| *x == 0.0)
    );
    let linear = scalar(3.0, 0.0, 1.0)
        .invariant_axis_spectrum(0, 0.6, 16, budget())
        .unwrap();
    assert_eq!(linear.cosine_amplitudes[0], 0.6);
    assert_eq!(linear.displacement_dc, 0.0);
    assert_eq!(linear.guard.displacement_l1, 0.0);
}

#[test]
fn physical_units_survive_extreme_length_and_frequency_scales() {
    let original = scalar(3.0, 20.0, 0.4)
        .invariant_axis_spectrum(0, 0.6, 16, budget())
        .unwrap();
    for (length, rate) in [
        (2.0_f64.powi(-400), 1.0),
        (2.0_f64.powi(400), 1.0),
        (1.0, 2.0_f64.powi(-400)),
        (1.0, 2.0_f64.powi(508)),
    ] {
        let storage = scalar(
            3.0 * rate,
            20.0 * (rate / length) * (rate / length),
            0.4 * length,
        );
        let scaled = storage
            .invariant_axis_spectrum(0, 0.6 * length, 16, budget())
            .unwrap();
        close(
            scaled.angular_frequency / rate,
            original.angular_frequency,
            original.angular_frequency,
            5e-14,
        );
        close(
            scaled.displacement_dc / length,
            original.displacement_dc,
            0.6,
            5e-14,
        );
        for j in 0..16 {
            close(
                scaled.cosine_amplitudes[j] / length,
                original.cosine_amplitudes[j],
                1.0,
                5e-14,
            );
            close(
                scaled.velocity_sine_amplitudes[j] / length / rate,
                original.velocity_sine_amplitudes[j],
                10.0,
                5e-14,
            );
            close(
                scaled.acceleration_cosine_amplitudes[j] / length / rate / rate,
                original.acceleration_cosine_amplitudes[j],
                100.0,
                5e-14,
            );
        }
    }
}

#[test]
fn original_cross_mode_and_preload_coupling_refuse_without_cutoffs() {
    let source = |coupling| SosModalStorage {
        omegas: vec![3.0, 4.0],
        channels: vec![StressChannel {
            coefficient: 20.0,
            coupling,
        }],
    };
    let diagonal =
        PreloadedModalStorage::new(source(vec![1.0, 0.0, 0.0, 2.0]), &[0.4, 0.0]).unwrap();
    let scalar = scalar(3.0, 20.0, 0.4);
    assert_eq!(
        diagonal
            .invariant_axis_spectrum(0, 0.6, 16, budget())
            .unwrap(),
        scalar
            .invariant_axis_spectrum(0, 0.6, 16, budget())
            .unwrap()
    );
    let biased_other =
        PreloadedModalStorage::new(source(vec![1.0, 0.0, 0.0, 2.0]), &[0.4, 0.1]).unwrap();
    assert_eq!(
        biased_other.invariant_axis_spectrum(0, 0.6, 16, budget()),
        Err(PreloadedSpectrumError::NonInvariantAxis)
    );
    let tiny =
        PreloadedModalStorage::new(source(vec![1.0, 1e-280, 1e-280, 2.0]), &[0.0, 0.0]).unwrap();
    assert_eq!(
        tiny.invariant_axis_spectrum(0, 0.6, 16, budget()),
        Err(PreloadedSpectrumError::NonInvariantAxis)
    );
    let mixed = PreloadedModalStorage::new(source(vec![0.0, 1.0, 1.0, 0.0]), &[0.0, 0.3]).unwrap();
    assert_eq!(
        mixed.invariant_axis_spectrum(0, 0.6, 16, budget()),
        Err(PreloadedSpectrumError::NonInvariantAxis)
    );
}

#[test]
fn refinement_reports_guard_content_and_refuses_insufficient_work() {
    let storage = scalar(3.0, 20.0, 0.4);
    let low = storage
        .invariant_axis_spectrum(0, 0.6, 2, budget())
        .unwrap();
    let high = storage
        .invariant_axis_spectrum(0, 0.6, 16, budget())
        .unwrap();
    assert!(low.guard.acceleration_l1 > 1e-3);
    assert!(high.guard.acceleration_l1 < low.guard.acceleration_l1 * 1e-5);
    close(
        low.cosine_amplitudes[0],
        high.cosine_amplitudes[0],
        1.0,
        2e-9,
    );
    close(
        low.cosine_amplitudes[1],
        high.cosine_amplitudes[1],
        1.0,
        2e-9,
    );
    assert_eq!(low.guard.last_harmonic, 4);
    assert_eq!(high.guard.last_harmonic, 32);
    let too_small = PreloadedSpectrumBudget {
        max_intervals: 64,
        tolerance: 1e-12,
    };
    assert!(matches!(
        storage.invariant_axis_spectrum(0, 0.6, 2, too_small),
        Err(PreloadedSpectrumError::QuadratureUnresolved { intervals: 64, .. })
    ));
    for (mode, amplitude, order) in [(1, 0.6, 2), (0, -1.0, 2), (0, 0.6, 0), (0, 0.6, 129)] {
        assert_eq!(
            storage.invariant_axis_spectrum(mode, amplitude, order, budget()),
            Err(PreloadedSpectrumError::InvalidInput)
        );
    }
}
