//! Original mixed source records retain coherent physical sound power.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::array::recorded::SourceSampleClock;
use fs_aeroac::curle3d::array::recorded::joint::{
    JointRadiatedPowerBandBudget, JointSourceRecordingBudget, JointSourceWelchStream,
    RecordedJointSourceSpectrum,
};
use fs_aeroac::curle3d::power::{RadiatedPowerBudget, RadiatedPowerModeBudget};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchDetrend, WelchWindow};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.23, 343.0).unwrap()
}

fn spatial(n: usize) -> RadiatedPowerBudget {
    RadiatedPowerBudget {
        max_sources: n,
        max_source_pairs: n * n,
    }
}

fn modes(n: usize, factors: usize) -> RadiatedPowerModeBudget {
    RadiatedPowerModeBudget {
        max_modes: factors,
        max_source_pair_mode_products: n * n * factors,
    }
}

fn band_budget(n: usize, factors: usize, bins: usize) -> JointRadiatedPowerBandBudget {
    JointRadiatedPowerBandBudget {
        spatial: spatial(n),
        modes: modes(n, factors),
        max_bins: bins,
        max_source_pair_mode_bin_products: n * n * factors * bins,
    }
}

#[allow(clippy::too_many_arguments)]
fn record(
    compact_positions: &[[f64; 3]],
    stress_positions: &[[f64; 3]],
    n: usize,
    samples: usize,
    hop: usize,
    fs: f64,
    window: WelchWindow,
    sample: impl Fn(usize) -> (Vec<[f64; 4]>, Vec<[f64; 6]>),
) -> RecordedJointSourceSpectrum {
    let clock = SourceSampleClock::new(0.125, 1.0 / fs, 0.0).unwrap();
    let channels = compact_positions.len() * 4 + stress_positions.len() * 6;
    let mut stream = JointSourceWelchStream::new(
        medium(),
        compact_positions,
        stress_positions,
        clock,
        samples,
        WelchConfig {
            segment_samples: n,
            hop_samples: hop,
            window,
            detrend: WelchDetrend::None,
        },
        JointSourceRecordingBudget {
            maximum_compact_sources: compact_positions.len(),
            maximum_quadrupole_sources: stress_positions.len(),
            welch: WelchBudget {
                maximum_channels: channels,
                maximum_samples_per_channel: samples,
                maximum_transforms: channels * ((samples - n) / hop + 1),
                maximum_workspace_bytes: 4_000_000,
            },
        },
        || false,
    )
    .unwrap();
    for index in 0..samples {
        let (compact, stress) = sample(index);
        stream
            .push_sample(
                compact_positions,
                stress_positions,
                clock.time_at(index).unwrap(),
                &compact,
                &stress,
                || false,
            )
            .unwrap();
    }
    stream.finish().unwrap()
}

fn wave(value: C64, angle: f64) -> f64 {
    value.re * angle.cos() + value.im * angle.sin()
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-30),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn monopole_power(frequency: f64, q: C64) -> f64 {
    medium().density_kg_m3() * (TAU * frequency).powi(2) * (q.re * q.re + q.im * q.im)
        / (8.0 * PI * medium().sound_speed_m_s())
}

#[test]
fn colocated_compact_and_isotropic_stress_cancel_or_reinforce_before_power() {
    let (n, bin, fs) = (64, 7, 8192.0);
    let frequency = bin as f64 * fs / n as f64;
    let q = C64::new(2e-5, -3e-5);
    let k = TAU * frequency / medium().sound_speed_m_s();
    // rho*c*Q - i*k*T = 0 for this actual isotropic stress, with e^-iwt.
    let stress = C64::new(0.0, -medium().characteristic_impedance_pa_s_m() / k) * q;
    let single = monopole_power(frequency, q);
    for sign in [1.0, -1.0, 0.0] {
        let record = record(
            &[[0.2, 0.1, -0.3]],
            &[[0.2, 0.1, -0.3]],
            n,
            3 * n,
            n,
            fs,
            WelchWindow::Rectangular,
            |index| {
                let theta = TAU * bin as f64 * index as f64 / n as f64;
                let t = sign * wave(stress, theta);
                (
                    vec![[wave(q, theta), 0.0, 0.0, 0.0]],
                    vec![[t, t, t, 0.0, 0.0, 0.0]],
                )
            },
        );
        let mut prepared = record
            .prepare_radiated_power(bin, spatial(2), || false)
            .unwrap();
        assert!(core::ptr::eq(prepared.record(), &record));
        assert_eq!(prepared.record().clock(), record.clock());
        assert_eq!(prepared.power_operator().source_count(), 2);
        assert_eq!(
            prepared.power_operator().compact_power().position_count(),
            1
        );
        let row = prepared.evaluate(modes(2, 3), || false).unwrap();
        if sign == 1.0 {
            assert!(
                row.power_w < single * 2e-24,
                "coherent null lost: {}",
                row.power_w
            );
        } else {
            close(row.power_w, single * (1.0 - sign).powi(2), 2e-12);
        }
        close(
            row.power_density_w_hz * row.bin_width_hz,
            row.power_w,
            3e-15,
        );
        assert_eq!(row, prepared.evaluate(modes(2, 3), || false).unwrap());
    }
}

fn legendre(order: usize, x: f64) -> (f64, f64) {
    let (mut previous, mut value) = (0.0, 1.0);
    for degree in 1..=order {
        let next =
            ((2 * degree - 1) as f64 * x * value - (degree - 1) as f64 * previous) / degree as f64;
        previous = value;
        value = next;
    }
    (value, previous)
}

fn sphere_integral(order: usize, field: impl Fn([f64; 3]) -> f64) -> f64 {
    let mut integral = 0.0;
    for i in 0..order {
        let mut z = (PI * (i as f64 + 0.75) / (order as f64 + 0.5)).cos();
        for _ in 0..12 {
            let (p, previous) = legendre(order, z);
            let derivative = order as f64 * (z * p - previous) / (z * z - 1.0);
            z -= p / derivative;
        }
        let (p, previous) = legendre(order, z);
        let derivative = order as f64 * (z * p - previous) / (z * z - 1.0);
        let weight = 2.0 / ((1.0 - z * z) * derivative * derivative);
        for j in 0..2 * order {
            let phi = TAU * j as f64 / (2 * order) as f64;
            let radius = (1.0 - z * z).sqrt();
            integral += weight * TAU / (2 * order) as f64
                * field([radius * phi.cos(), radius * phi.sin(), z]);
        }
    }
    integral
}

#[test]
fn separated_mixed_record_matches_independent_radiative_sphere_flux() {
    let compact_positions = [[0.04, -0.1, 0.03], [-0.08, 0.02, -0.05]];
    let stress_positions = [[0.11, 0.05, -0.02]];
    let compact = [
        [
            C64::new(2e-5, -1e-5),
            C64::new(0.012, -0.006),
            C64::new(-0.003, 0.008),
            C64::new(0.002, 0.005),
        ],
        [
            C64::new(-1e-5, 3e-5),
            C64::new(0.004, 0.007),
            C64::new(0.009, -0.002),
            C64::new(-0.006, 0.003),
        ],
    ];
    let stress = [
        C64::new(0.002, -0.001),
        C64::new(-0.001, 0.0007),
        C64::new(0.0008, 0.0003),
        C64::new(0.0015, -0.0004),
        C64::new(-0.0006, 0.0009),
        C64::new(0.0005, 0.0002),
    ];
    let phases = [
        (C64::ONE, C64::ONE),
        (C64::new(0.0, 1.0), C64::new(0.0, -1.0)),
        (C64::new(-0.4, 0.3), C64::new(0.6, -0.2)),
    ];
    let (n, bin, fs) = (64, 5, 8192.0);
    let record = record(
        &compact_positions,
        &stress_positions,
        n,
        3 * n,
        n,
        fs,
        WelchWindow::Rectangular,
        |index| {
            let theta = TAU * bin as f64 * (index % n) as f64 / n as f64;
            let phase = phases[index / n];
            (
                vec![
                    compact[0].map(|c| wave(c, theta)),
                    compact[1].map(|c| wave(c * phase.0, theta)),
                ],
                vec![stress.map(|c| wave(c * phase.1, theta))],
            )
        },
    );
    let row = record
        .prepare_radiated_power(bin, spatial(3), || false)
        .unwrap()
        .evaluate(modes(3, 3), || false)
        .unwrap();
    let k = TAU * row.frequency_hz / medium().sound_speed_m_s();
    // Independent far-zone outgoing p/u flux: integrate the explicit physical
    // amplitude over solid angle, without production Bessel moments or tensors.
    let expected = phases
        .iter()
        .map(|&(compact_phase, stress_phase)| {
            sphere_integral(36, |direction| {
                let phase_at = |position: [f64; 3]| {
                    let angle = -k * (0..3).map(|a| direction[a] * position[a]).sum::<f64>();
                    C64::new(angle.cos(), angle.sin())
                };
                let mut amplitude = C64::ZERO;
                for source in 0..2 {
                    let mut a =
                        compact[source][0].scale(medium().characteristic_impedance_pa_s_m());
                    for axis in 0..3 {
                        a = a + compact[source][axis + 1].scale(direction[axis]);
                    }
                    if source == 1 {
                        a = a * compact_phase;
                    }
                    amplitude = amplitude + a * phase_at(compact_positions[source]);
                }
                let [x, y, z] = direction;
                let weights = [x * x, y * y, z * z, 2.0 * x * y, 2.0 * x * z, 2.0 * y * z];
                let tensor = (0..6).fold(C64::ZERO, |sum, a| sum + stress[a].scale(weights[a]));
                amplitude = amplitude
                    + C64::new(0.0, -k) * tensor * stress_phase * phase_at(stress_positions[0]);
                k * k * (amplitude.re * amplitude.re + amplitude.im * amplitude.im)
                    / (32.0 * PI * PI * medium().characteristic_impedance_pa_s_m())
            })
        })
        .sum::<f64>()
        / phases.len() as f64;
    close(row.power_w, expected, 3e-12);
}

#[test]
fn complete_bin_band_obeys_windowed_derivative_parseval_with_overlap() {
    let (n, fs) = (128, 8192.0);
    let tones = [(5.0, C64::new(3e-5, -2e-5)), (19.0, C64::new(-2e-5, -1e-5))];
    for window in [WelchWindow::Rectangular, WelchWindow::PeriodicHann] {
        let record = record(
            &[[0.0; 3]],
            &[[0.0; 3]],
            n,
            3 * n,
            n / 2,
            fs,
            window,
            |index| {
                let q = tones
                    .iter()
                    .map(|&(b, a)| wave(a, TAU * b * index as f64 / n as f64))
                    .sum();
                (vec![[q, 0.0, 0.0, 0.0]], vec![[0.0; 6]])
            },
        );
        assert_eq!(record.spectrum().segment_count(), 5);
        let band = record
            .radiated_power_band(1..n / 2, band_budget(2, 5, n / 2 - 1), || false)
            .unwrap();
        // Parseval applied to the actual windowed source, differentiated
        // analytically. Window energy, including the window derivative's
        // leakage contribution, is independent of the production FFT factors.
        let omega = TAU * fs / n as f64;
        let (mut square, mut energy) = (0.0, 0.0);
        for j in 0..n {
            let theta = TAU * j as f64 / n as f64;
            let q = tones.iter().map(|&(b, a)| wave(a, b * theta)).sum::<f64>();
            let dq = tones
                .iter()
                .map(|&(b, a)| b * omega * (-a.re * (b * theta).sin() + a.im * (b * theta).cos()))
                .sum::<f64>();
            let (w, dw) = if window == WelchWindow::PeriodicHann {
                (0.5 * (1.0 - theta.cos()), 0.5 * omega * theta.sin())
            } else {
                (1.0, 0.0)
            };
            square += (w * dq + dw * q).powi(2);
            energy += w * w;
        }
        let expected =
            medium().density_kg_m3() * square / (4.0 * PI * medium().sound_speed_m_s() * energy);
        close(band.power_w, expected, 4e-12);
        close(
            band.bins
                .iter()
                .map(|row| row.power_density_w_hz * row.bin_width_hz)
                .sum(),
            expected,
            4e-12,
        );
        assert_eq!(
            band.frequency_edges_hz,
            [fs / (2 * n) as f64, (n as f64 / 2.0 - 0.5) * fs / n as f64]
        );
        let lower = record
            .radiated_power_band(1..12, band_budget(2, 5, 11), || false)
            .unwrap();
        let upper = record
            .radiated_power_band(12..n / 2, band_budget(2, 5, n / 2 - 12), || false)
            .unwrap();
        close(lower.power_w + upper.power_w, band.power_w, 3e-15);
    }
}

#[test]
fn odd_last_bin_is_physical_and_even_nyquist_is_refused() {
    let fs = 9600.0;
    let q = C64::new(1e-5, -3e-5);
    for n in [15, 16] {
        let bin = n / 2;
        let record = record(
            &[[0.0; 3]],
            &[[0.1, 0.0, 0.0]],
            n,
            n,
            n,
            fs,
            WelchWindow::Rectangular,
            |index| {
                (
                    vec![[
                        wave(q, TAU * bin as f64 * index as f64 / n as f64),
                        0.0,
                        0.0,
                        0.0,
                    ]],
                    vec![[0.0; 6]],
                )
            },
        );
        assert!(
            record
                .prepare_radiated_power(0, spatial(2), || false)
                .is_err()
        );
        if n % 2 == 0 {
            assert!(
                record
                    .prepare_radiated_power(bin, spatial(2), || false)
                    .is_err()
            );
            assert!(
                record
                    .radiated_power_band(1..bin + 1, band_budget(2, 1, bin), || false)
                    .is_err()
            );
        } else {
            let row = record
                .prepare_radiated_power(bin, spatial(2), || false)
                .unwrap()
                .evaluate(modes(2, 1), || false)
                .unwrap();
            close(
                row.power_w,
                monopole_power(bin as f64 * fs / n as f64, q),
                4e-12,
            );
            let band = record
                .radiated_power_band(bin..bin + 1, band_budget(2, 1, 1), || false)
                .unwrap();
            assert_eq!(band.frequency_edges_hz[1], fs / 2.0);
        }
    }
}

#[test]
fn combined_work_budgets_and_cancellation_leave_original_record_reusable() {
    let record = record(
        &[[0.0; 3]],
        &[[0.1, 0.0, 0.0]],
        32,
        64,
        32,
        8192.0,
        WelchWindow::Rectangular,
        |index| {
            let angle = TAU * 3.0 * index as f64 / 32.0;
            (
                vec![[angle.cos() * 1e-5, angle.sin() * 0.002, 0.0, 0.0]],
                vec![[angle.sin() * 0.003; 6]],
            )
        },
    );
    let mut short_spatial = spatial(2);
    short_spatial.max_source_pairs -= 1;
    assert!(
        record
            .prepare_radiated_power(3, short_spatial, || false)
            .is_err()
    );
    assert!(matches!(
        record.prepare_radiated_power(3, spatial(2), || true),
        Err(AeroacError::Cancelled { .. })
    ));
    let mut prepared = record
        .prepare_radiated_power(3, spatial(2), || false)
        .unwrap();
    let expected = prepared.evaluate(modes(2, 2), || false).unwrap();
    assert!(prepared.evaluate(modes(2, 1), || false).is_err());
    let mut calls = 0;
    assert!(matches!(
        prepared.evaluate(modes(2, 2), || {
            calls += 1;
            calls == 7
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    assert_eq!(prepared.evaluate(modes(2, 2), || false).unwrap(), expected);
    let mut short_band = band_budget(2, 2, 4);
    short_band.max_source_pair_mode_bin_products -= 1;
    let mut polls = 0;
    assert!(
        record
            .radiated_power_band(1..5, short_band, || {
                polls += 1;
                false
            })
            .is_err()
    );
    assert_eq!(polls, 1, "total work must refuse before preparing any bin");
    let mut calls = 0;
    assert!(matches!(
        record.radiated_power_band(1..5, band_budget(2, 2, 4), || {
            calls += 1;
            calls == 35
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    let band = record
        .radiated_power_band(1..5, band_budget(2, 2, 4), || false)
        .unwrap();
    assert_eq!(band.bins[2], expected);
    assert_eq!(record.spectrum().used_sample_count(), 64);
    assert!(
        record
            .radiated_power_band(3..3, band_budget(2, 2, 4), || false)
            .is_err()
    );
}
