use super::*;
use core::f64::consts::{PI, TAU};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchDetrend, WelchWindow};

const N: usize = 32;
const BIN: usize = 3;
const DT: f64 = 1.0 / 4096.0;
const MEDIUM: Medium = Medium {
    density: 1.2,
    sound_speed: 340.0,
};
const DIRECTIONS: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]];
const WORK: HarmonicWorkBudget = HarmonicWorkBudget {
    maximum_products: MAX_PRODUCTS,
};

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * scale,
        "{actual:e} != {expected:e}"
    );
}

fn spectrum(segments: &[[C64; 2]], window: WelchWindow) -> WelchCrossSpectrum {
    let rows: Vec<[f64; 2]> = segments
        .iter()
        .flat_map(|a| {
            (0..N).map(move |j| {
                let theta = TAU * BIN as f64 * j as f64 / N as f64;
                a.map(|z| z.re * theta.cos() + z.im * theta.sin())
            })
        })
        .collect();
    WelchCrossSpectrum::estimate_interleaved(
        &rows,
        DT,
        WelchConfig {
            segment_samples: N,
            hop_samples: N,
            window,
            detrend: WelchDetrend::None,
        },
        WelchBudget {
            maximum_channels: 2,
            maximum_samples_per_channel: 128,
            maximum_transforms: 8,
            maximum_workspace_bytes: 100_000,
        },
        || false,
    )
    .unwrap()
}

fn bank(
    spectrum: &WelchCrossSpectrum,
    refinements: usize,
    field: impl Fn(usize, [f64; 3]) -> C64,
) -> PreparedRadiationObservation {
    PreparedRadiationObservation::prepare::<TransferError>(
        TAU * spectrum.frequency_hz(BIN).unwrap(),
        MEDIUM,
        &[
            SourceId::PointVolumeFlow("Q".into()),
            SourceId::PointForce("Fz".into()),
        ],
        ContentHash([19; 32]),
        &DIRECTIONS,
        PreparedRadiationOptions {
            angular: helmholtz::DirectivityProjectionOptions::default(),
            initial_latitudes: 8,
            refinements,
            maximum_kernel_evaluations: 1_000_000,
            maximum_source_entries: MAX_OBSERVATION_ENTRIES,
        },
        2,
        |directions| {
            let mut columns = Vec::new();
            for source in 0..2 {
                columns.extend(directions.iter().map(|&d| field(source, d)));
            }
            Ok(columns)
        },
        &|| false,
    )
    .unwrap()
}

#[test]
fn measured_mixed_source_power_and_signed_directional_moments_match_sphere_integrals() {
    let segments = [
        [C64::new(0.001, 0.002), C64::new(0.8, -0.3)],
        [C64::new(-0.0007, 0.001), C64::new(-0.4, 0.9)],
        [C64::new(0.0013, -0.0006), C64::new(0.6, 0.2)],
    ];
    for (window, fraction) in [
        (WelchWindow::Rectangular, 1.0),
        (WelchWindow::PeriodicHann, 2.0 / 3.0),
    ] {
        let spectrum = spectrum(&segments, window);
        let omega = TAU * spectrum.frequency_hz(BIN).unwrap();
        let k = omega / MEDIUM.sound_speed;
        let monopole = C64::new(0.0, -omega * MEDIUM.density / (4.0 * PI));
        let dipole = C64::new(0.0, -k / (4.0 * PI));
        let bank = bank(&spectrum, 4, |source, d| {
            if source == 0 {
                monopole
            } else {
                dipole.scale(d[2])
            }
        });
        let mut density = [0.0; 3];
        let actual = bank
            .evaluate_welch_bin(&spectrum, BIN, bank.sources(), WORK, &mut density, || false)
            .unwrap();
        let z = MEDIUM.density * MEDIUM.sound_speed;
        let mut total = 0.0;
        let mut first_z = 0.0;
        let mut second = [0.0; 3];
        let mut directional = [0.0; 3];
        for source in segments {
            let a = monopole * source[0];
            let b = dipole * source[1];
            let q_power = 2.0 * PI * a.norm_sq() / z;
            let f_power = 2.0 * PI * b.norm_sq() / (3.0 * z);
            let weight = fraction / segments.len() as f64;
            total += weight * (q_power + f_power);
            first_z += weight * 4.0 * PI * (a * b.conj()).re / (3.0 * z);
            second[0] += weight * (q_power / 3.0 + f_power / 5.0);
            second[1] = second[0];
            second[2] += weight * (q_power / 3.0 + 3.0 * f_power / 5.0);
            for d in 0..3 {
                directional[d] += weight * (a + b.scale(DIRECTIONS[d][2])).norm_sq() / (2.0 * z);
            }
        }
        close(actual.total_power_w, total, total);
        close(actual.first_directional_moment_w[2], first_z, total);
        for i in 0..2 {
            close(actual.first_directional_moment_w[i], 0.0, total);
        }
        for i in 0..3 {
            close(density[i], directional[i], total);
            for j in 0..3 {
                close(
                    actual.second_directional_moment_w[i][j],
                    if i == j { second[i] } else { 0.0 },
                    total,
                );
            }
        }
        assert!(first_z.abs() > 1e-4 * total);
        let mut direct_only = [0.0; 3];
        bank.observe_welch_bin(
            &spectrum,
            BIN,
            bank.sources(),
            WORK,
            &mut direct_only,
            || false,
        )
        .unwrap();
        assert_eq!(direct_only, density);
    }
}

#[test]
fn measured_radiated_power_retains_exact_nulls_and_partial_coherence() {
    let q = C64::new(0.001, 0.0002);
    for (segments, ratio) in [
        ([[q, q], [q, q]], 4.0),
        ([[q, -q], [q, -q]], 0.0),
        ([[q, q], [q, -q]], 2.0),
    ] {
        let spectrum = spectrum(&segments, WelchWindow::Rectangular);
        // Two identical original angular columns, with their measured
        // coordinate units explicitly chosen to give the same field per unit.
        let bank = bank(&spectrum, 4, |_, _| C64::new(0.7, -0.2));
        let mut density = [0.0; 3];
        let power = bank
            .evaluate_welch_bin(&spectrum, BIN, bank.sources(), WORK, &mut density, || false)
            .unwrap();
        let single =
            (C64::new(0.7, -0.2) * q).norm_sq() / (2.0 * MEDIUM.density * MEDIUM.sound_speed);
        close(power.total_power_w, 4.0 * PI * single * ratio, single);
        for value in density {
            close(value, single * ratio, single);
        }
        if ratio == 0.0 {
            assert_eq!(power.total_power_w, 0.0);
            assert_eq!(power.first_directional_moment_w, [0.0; 3]);
            assert_eq!(power.second_directional_moment_w, [[0.0; 3]; 3]);
        }
    }
}

#[test]
fn recorded_covariance_must_pass_its_own_angular_convergence_gate() {
    let simple = spectrum(&[[C64::ONE, C64::ZERO]], WelchWindow::Rectangular);
    let intricate = spectrum(&[[C64::ZERO, C64::ONE]], WelchWindow::Rectangular);
    let bank = bank(&simple, 3, |source, d| {
        if source == 0 {
            return C64::ONE;
        }
        let (mut a, mut b) = (1.0, d[2]);
        for l in 2..=20 {
            let c = ((2 * l - 1) as f64 * d[2] * b - (l - 1) as f64 * a) / l as f64;
            a = b;
            b = c;
        }
        C64::from_re(b)
    });
    bank.evaluate_welch_bin(&simple, BIN, bank.sources(), WORK, &mut [0.0; 3], || false)
        .unwrap();
    assert!(matches!(
        bank.evaluate_welch_bin(&intricate, BIN, bank.sources(), WORK, &mut [0.0; 3], || {
            false
        }),
        Err(TransferError::AngularPowerNotConverged { .. })
    ));
    // Requested directions are well-defined even when the full-sphere bank
    // is too coarse to claim its integrated outgoing power.
    bank.observe_welch_bin(&intricate, BIN, bank.sources(), WORK, &mut [0.0; 3], || {
        false
    })
    .unwrap();
}

#[test]
fn recorded_radiation_refuses_identity_work_and_cancellation_then_replays() {
    let spectrum = spectrum(&[[C64::ONE, C64::new(0.3, 0.7)]], WelchWindow::Rectangular);
    let bank = bank(&spectrum, 4, |source, d| {
        if source == 0 {
            C64::ONE
        } else {
            C64::new(0.0, d[2])
        }
    });
    let mut density = [0.0; 3];
    for bin in [0, N / 2, BIN + 1, usize::MAX] {
        assert!(
            bank.evaluate_welch_bin(&spectrum, bin, bank.sources(), WORK, &mut density, || false)
                .is_err()
        );
    }
    let reversed = [bank.sources()[1].clone(), bank.sources()[0].clone()];
    assert!(
        bank.evaluate_welch_bin(&spectrum, BIN, &reversed, WORK, &mut density, || false)
            .is_err()
    );
    assert!(
        bank.evaluate_welch_bin(
            &spectrum,
            BIN,
            bank.sources(),
            HarmonicWorkBudget {
                maximum_products: 1
            },
            &mut density,
            || false
        )
        .is_err()
    );
    assert!(
        bank.evaluate_welch_bin(
            &spectrum,
            BIN,
            bank.sources(),
            WORK,
            &mut density[..1],
            || false
        )
        .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        bank.evaluate_welch_bin(&spectrum, BIN, bank.sources(), WORK, &mut density, || {
            polls += 1;
            polls == 10
        }),
        Err(TransferError::Cancelled)
    ));
    let first = bank
        .evaluate_welch_bin(&spectrum, BIN, bank.sources(), WORK, &mut density, || false)
        .unwrap();
    let original_density = density;
    let replay = bank
        .evaluate_welch_bin(&spectrum, BIN, bank.sources(), WORK, &mut density, || false)
        .unwrap();
    assert_eq!(first, replay);
    assert_eq!(original_density, density);
}
