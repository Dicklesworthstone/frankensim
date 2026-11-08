//! Original mixed-source clocks, spectral normalization and passive wall dispersion.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::QuadrupolePhasor;
use fs_aeroac::curle3d::quadrupole::recorded::{
    QuadrupoleRecordingBudget, QuadrupoleSampleClock, QuadrupoleWelchStream,
};
use fs_aeroac::curle3d::refraction::full_wave::array::FullWaveCompactSourceLimits;
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::FullWaveQuadrupoleLimits;
use fs_aeroac::curle3d::refraction::full_wave::quadrupole::array::FullWaveQuadrupoleSourceLimits;
use fs_aeroac::curle3d::refraction::full_wave::wall::recorded::{
    ImpedanceWallWelchBandLimits, observe_welch_band,
};
use fs_aeroac::curle3d::refraction::full_wave::wall::{
    ImpedanceWallFieldKind as Kind, PlanarImpedanceWall,
    PreparedImpedanceWallCompactSourceObserver as Compact,
    PreparedImpedanceWallQuadrupoleObserver as Stress,
};
use fs_aeroac::curle3d::refraction::full_wave::{FullWaveForceLimits, FullWaveMonopoleLimits};
use fs_aeroac::curle3d::{
    CoherenceBudget, CoherentCompactField, CompactFieldStatistics, CompactSourcePhasor,
};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchDetrend, WelchWindow};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

const N: usize = 64;
const DT: f64 = 1.0 / 8192.0;
const RHO: f64 = 1.2;
const C: f64 = 343.0;
const COMPACT: [[f64; 3]; 2] = [[0.01, 0.0, -0.09], [0.05, -0.02, -0.12]];
const STRESS: [[f64; 3]; 1] = [[-0.01, 0.02, -0.07]];
const RECEIVER: [f64; 3] = [0.16, 0.04, -0.1];
const PAIRS: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

fn wall(z: Option<C64>) -> PlanarImpedanceWall {
    PlanarImpedanceWall::new(
        FluidAcousticState::try_new(RHO, C).unwrap(),
        [0.0; 3],
        [0.0, 0.0, 1.0],
        z,
    )
    .unwrap()
}
fn limits() -> ImpedanceWallWelchBandLimits {
    ImpedanceWallWelchBandLimits {
        compact: FullWaveCompactSourceLimits {
            max_sources: 2,
            maximum_total_evaluations: 524_288,
            monopole: FullWaveMonopoleLimits {
                absolute_pressure_tolerance_pa_per_m3_s: 1e-6,
                absolute_velocity_tolerance_per_m2: 1e-9,
                relative_tolerance: 1e-9,
                maximum_evaluations: 131_072,
                maximum_depth: 40,
            },
            force: FullWaveForceLimits {
                absolute_pressure_tolerance_pa_per_n: 1e-8,
                absolute_velocity_tolerance_m_s_per_n: 1e-11,
                relative_tolerance: 1e-9,
                maximum_evaluations: 131_072,
                maximum_depth: 40,
            },
        },
        quadrupole: FullWaveQuadrupoleSourceLimits {
            max_sources: 2,
            maximum_total_evaluations: 524_288,
            quadrupole: FullWaveQuadrupoleLimits {
                absolute_pressure_tolerance_pa_per_n_m: 1e-7,
                absolute_velocity_tolerance_m_s_per_n_m: 1e-10,
                relative_tolerance: 1e-9,
                maximum_evaluations: 131_072,
                maximum_depth: 40,
            },
        },
        maximum_bins: 6,
        maximum_total_evaluations: 8_388_608,
        coherence: CoherenceBudget {
            max_modes: 12,
            max_source_mode_pairs: 36,
        },
    }
}
fn config(window: WelchWindow) -> WelchConfig {
    WelchConfig {
        segment_samples: N,
        hop_samples: N,
        window,
        detrend: WelchDetrend::None,
    }
}
fn fft_budget() -> WelchBudget {
    WelchBudget {
        maximum_channels: 14,
        maximum_samples_per_channel: 2 * N,
        maximum_transforms: 28,
        maximum_workspace_bytes: 1_000_000,
    }
}
fn estimate<const D: usize>(samples: &[[f64; D]], window: WelchWindow) -> WelchCrossSpectrum {
    WelchCrossSpectrum::estimate_interleaved(samples, DT, config(window), fft_budget(), || false)
        .unwrap()
}
fn record<const D: usize>(tones: &[(usize, [[C64; D]; 2])]) -> Vec<[f64; D]> {
    (0..2 * N)
        .map(|j| {
            core::array::from_fn(|channel| {
                tones
                    .iter()
                    .map(|(bin, amplitudes)| {
                        let a = amplitudes[j / N][channel];
                        let angle = TAU * *bin as f64 * (j % N) as f64 / N as f64;
                        a.re * angle.cos() + a.im * angle.sin()
                    })
                    .sum()
            })
        })
        .collect()
}
fn amplitudes(segment: usize) -> [C64; 14] {
    core::array::from_fn(|i| {
        let scale = if i == 0 || i == 4 {
            1e-5
        } else if i < 8 {
            0.002
        } else {
            0.0001
        };
        let phase = if segment == 0 {
            C64::ONE
        } else if i < 8 {
            C64::new(0.2, 0.8)
        } else {
            C64::new(-0.7, 0.3)
        };
        C64::new((i + 1) as f64 * scale, (3.0 - i as f64) * 0.4 * scale) * phase
    })
}
fn sources(a: [C64; 14]) -> ([CompactSourcePhasor; 2], QuadrupolePhasor) {
    (
        core::array::from_fn(|i| CompactSourcePhasor {
            volume_flow_m3_s: a[4 * i],
            force_on_fluid_n: [a[4 * i + 1], a[4 * i + 2], a[4 * i + 3]],
        }),
        QuadrupolePhasor {
            stress_moment_n_m: core::array::from_fn(|i| a[8 + i]),
        },
    )
}
fn values(field: CoherentCompactField) -> [C64; 4] {
    [
        field.pressure_pa,
        field.particle_velocity_m_s[0],
        field.particle_velocity_m_s[1],
        field.particle_velocity_m_s[2],
    ]
}
fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 1e-17 + 2e-8 * expected.abs(),
        "{actual:?} != {expected:?}"
    );
}
fn moments(actual: &CompactFieldStatistics, factors: &[[C64; 4]]) {
    let entry = |i: usize, j: usize| {
        factors
            .iter()
            .fold(C64::ZERO, |sum, f| sum + f[i].scale(0.5) * f[j].conj())
    };
    close(C64::from_re(actual.mean_square_pressure_pa2), entry(0, 0));
    for i in 0..3 {
        for j in 0..3 {
            close(actual.velocity_covariance_m2_s2[i][j], entry(i + 1, j + 1));
        }
        close(
            C64::new(
                actual.active_intensity_w_m2[i],
                actual.reactive_intensity_w_m2[i],
            ),
            entry(0, i + 1),
        );
    }
}
fn observers(wall: PlanarImpedanceWall, f: f64, kind: Kind) -> (Compact, Stress) {
    let compact = match kind {
        Kind::Reflected => {
            Compact::new_reflected(wall, f, RECEIVER, &COMPACT, limits().compact, &mut || false)
        }
        Kind::Total => {
            Compact::new_total(wall, f, RECEIVER, &COMPACT, limits().compact, &mut || false)
        }
    }
    .unwrap();
    let stress = match kind {
        Kind::Reflected => {
            Stress::new_reflected(wall, f, RECEIVER, &STRESS, limits().quadrupole, &mut || {
                false
            })
        }
        Kind::Total => {
            Stress::new_total(wall, f, RECEIVER, &STRESS, limits().quadrupole, &mut || {
                false
            })
        }
    }
    .unwrap();
    (compact, stress)
}

// Independent radial derivatives of G=exp(i k r)/(4 pi r):
// p=-i rho omega Q G-F.grad(G)+M:Hess(G), u=grad(p)/(i rho omega).
fn green(
    f: f64,
    position: [f64; 3],
    qf: CompactSourcePhasor,
    moment: QuadrupolePhasor,
) -> [C64; 4] {
    let d: [f64; 3] = core::array::from_fn(|i| RECEIVER[i] - position[i]);
    let r = d.iter().map(|x| x * x).sum::<f64>().sqrt();
    let n = d.map(|x| x / r);
    let omega = TAU * f;
    let k = omega / C;
    let g = C64::new((k * r).cos(), (k * r).sin()).scale(1.0 / (4.0 * PI * r));
    let mut m = [[C64::ZERO; 3]; 3];
    for ((i, j), value) in PAIRS.into_iter().zip(moment.stress_moment_n_m) {
        m[i][j] = value;
        m[j][i] = value;
    }
    let mn: [C64; 3] =
        core::array::from_fn(|i| (0..3).fold(C64::ZERO, |s, j| s + m[i][j].scale(n[j])));
    let trace = m[0][0] + m[1][1] + m[2][2];
    let radial_m = (0..3).fold(C64::ZERO, |s, i| s + mn[i].scale(n[i]));
    let radial_f = (0..3).fold(C64::ZERO, |s, i| s + qf.force_on_fluid_n[i].scale(n[i]));
    let a = C64::new(-1.0 / r.powi(2), k / r);
    let b = C64::new(3.0 / r.powi(2) - k * k, -3.0 * k / r);
    let third = C64::new(
        6.0 * k * k / r - 15.0 / r.powi(3),
        -k.powi(3) + 15.0 * k / r.powi(2),
    );
    let radial = g * C64::new(1.0 / r, -k);
    let derivative = g * C64::new(k * k - 2.0 / r.powi(2), 2.0 * k / r);
    let u: [C64; 3] = core::array::from_fn(|i| {
        radial * qf.volume_flow_m3_s.scale(n[i])
            + (derivative * radial_f.scale(n[i])
                + radial.scale(1.0 / r) * (qf.force_on_fluid_n[i] - radial_f.scale(n[i]))
                + g * (b.scale(1.0 / r) * (mn[i].scale(2.0) + trace.scale(n[i]))
                    + third * radial_m.scale(n[i])))
                * C64::new(0.0, -1.0 / (RHO * omega))
    });
    [
        g * (qf.volume_flow_m3_s * C64::new(0.0, -RHO * omega) + a * trace + b * radial_m)
            + radial * radial_f,
        u[0],
        u[1],
        u[2],
    ]
}
fn ideal_field(f: f64, a: [C64; 14], sign: f64, kind: Kind) -> [C64; 4] {
    let (qf, moment) = sources(a);
    let mut result = [C64::ZERO; 4];
    for (position, compact, stress) in [
        (COMPACT[0], qf[0], QuadrupolePhasor::default()),
        (COMPACT[1], qf[1], QuadrupolePhasor::default()),
        (STRESS[0], CompactSourcePhasor::default(), moment),
    ] {
        let mut image = position;
        image[2] = -image[2];
        let mut image_qf = compact;
        image_qf.force_on_fluid_n[2] = -image_qf.force_on_fluid_n[2];
        let mut image_stress = stress;
        image_stress.stress_moment_n_m[4] = -image_stress.stress_moment_n_m[4];
        image_stress.stress_moment_n_m[5] = -image_stress.stress_moment_n_m[5];
        let direct = green(f, position, compact, stress);
        let reflected = green(f, image, image_qf, image_stress);
        for i in 0..4 {
            result[i] = result[i]
                + reflected[i].scale(sign)
                + if kind == Kind::Total {
                    direct[i]
                } else {
                    C64::ZERO
                };
        }
    }
    result
}

#[test]
fn original_joint_factors_match_independent_spherical_images_and_all_source_views() {
    let data = record(&[(7, [amplitudes(0), amplitudes(1)])]);
    let spectrum = estimate(&data, WelchWindow::Rectangular);
    let f = spectrum.frequency_hz(7).unwrap();
    for (z, sign) in [(None, 1.0), (Some(C64::ZERO), -1.0)] {
        for kind in [Kind::Reflected, Kind::Total] {
            let (compact, stress) = observers(wall(z), f, kind);
            let expected: Vec<_> = (0..2)
                .map(|s| {
                    ideal_field(f, amplitudes(s), sign, kind).map(|v| v.scale(1.0 / 2.0_f64.sqrt()))
                })
                .collect();
            let actual = compact
                .evaluate_welch_bin_with_quadrupoles(
                    &stress,
                    &spectrum,
                    7,
                    limits().coherence,
                    || false,
                )
                .unwrap();
            moments(&actual, &expected);
            let compact_record: Vec<[f64; 8]> = data
                .iter()
                .map(|row| core::array::from_fn(|i| row[i]))
                .collect();
            let stress_record: Vec<[f64; 6]> = data
                .iter()
                .map(|row| core::array::from_fn(|i| row[8 + i]))
                .collect();
            for only_stress in [false, true] {
                let expected: Vec<_> = (0..2)
                    .map(|s| {
                        let a = amplitudes(s);
                        ideal_field(
                            f,
                            core::array::from_fn(|i| {
                                if (i >= 8) == only_stress {
                                    a[i]
                                } else {
                                    C64::ZERO
                                }
                            }),
                            sign,
                            kind,
                        )
                        .map(|z| z.scale(1.0 / 2.0_f64.sqrt()))
                    })
                    .collect();
                let actual = if only_stress {
                    stress.evaluate_welch_bin(
                        &estimate(&stress_record, WelchWindow::Rectangular),
                        7,
                        limits().coherence,
                        || false,
                    )
                } else {
                    compact.evaluate_welch_bin(
                        &estimate(&compact_record, WelchWindow::Rectangular),
                        7,
                        limits().coherence,
                        || false,
                    )
                }
                .unwrap();
                moments(&actual, &expected);
            }
        }
    }
}

// Independent O(N²) windowed DFT, converted directly to peak-bin factors.
// For R segments: conj(DFT(w*x)) * 2/sqrt(N R sum(w²)). No FFT/PSD helper.
fn dft_factor(data: &[[f64; 14]], bin: usize, segment: usize) -> [C64; 14] {
    let weights: [f64; N] =
        core::array::from_fn(|j| 0.5 * (1.0 - (TAU * j as f64 / N as f64).cos()));
    let scale = 2.0 / (N as f64 * 2.0 * weights.iter().map(|w| w * w).sum::<f64>()).sqrt();
    core::array::from_fn(|channel| {
        (0..N).fold(C64::ZERO, |sum, j| {
            let angle = TAU * bin as f64 * j as f64 / N as f64;
            sum + C64::new(angle.cos(), angle.sin())
                .scale(weights[j] * data[segment * N + j][channel] * scale)
        })
    })
}
fn dispersive_wall(f: f64) -> Result<PlanarImpedanceWall, AeroacError> {
    // Passive series mass/spring/resistance under exp(-i omega t).
    let omega = TAU * f;
    Ok(wall(Some(C64::new(80.0, 3.0e6 / omega - 0.12 * omega))))
}

#[test]
fn broadband_wall_uses_each_native_frequency_and_retains_original_hann_factors() {
    let data = record(&[
        (5, [amplitudes(0), amplitudes(1)]),
        (
            8,
            [
                amplitudes(1).map(|z| z.scale(-0.7)),
                amplitudes(0).map(|z| z * C64::new(0.3, 0.4)),
            ],
        ),
    ]);
    let spectrum = estimate(&data, WelchWindow::PeriodicHann);
    let mut calls = Vec::new();
    let band = observe_welch_band(
        &spectrum,
        4..10,
        &COMPACT,
        &STRESS,
        RECEIVER,
        Kind::Total,
        |f| {
            calls.push(f);
            dispersive_wall(f)
        },
        limits(),
        || false,
    )
    .unwrap();
    let mut expected = Vec::new();
    for bin in 4..10 {
        let f = spectrum.frequency_hz(bin).unwrap();
        assert_eq!(calls[bin - 4].to_bits(), f.to_bits());
        let (compact, stress) = observers(dispersive_wall(f).unwrap(), f, Kind::Total);
        for segment in 0..2 {
            let (qf, moment) = sources(dft_factor(&data, bin, segment));
            let a = values(compact.evaluate(&qf).unwrap());
            let b = values(stress.evaluate(&[moment]).unwrap());
            expected.push(core::array::from_fn(|i| a[i] + b[i]));
        }
    }
    moments(&band.statistics, &expected);
    assert_eq!(calls.len(), 6);
    assert_eq!(band.bins, 4..10);
    assert_eq!(band.bin_width_hz, 128.0);
    assert_eq!(band.center_frequency_range_hz, [512.0, 1152.0]);
    assert!(band.evaluations > 0);
    let tight = ImpedanceWallWelchBandLimits {
        maximum_total_evaluations: band.evaluations - 1,
        ..limits()
    };
    assert!(
        observe_welch_band(
            &spectrum,
            4..10,
            &COMPACT,
            &STRESS,
            RECEIVER,
            Kind::Total,
            dispersive_wall,
            tight,
            || false
        )
        .is_err()
    );
}

#[test]
fn recorded_cross_mechanism_cancellation_preserves_pressure_null_velocity() {
    let frequency = 7.0 / (N as f64 * DT);
    let (compact, stress) = observers(dispersive_wall(frequency).unwrap(), frequency, Kind::Total);
    let qf = [
        CompactSourcePhasor::monopole(C64::from_re(2e-5)),
        CompactSourcePhasor::default(),
    ];
    let mut unit = QuadrupolePhasor::default();
    unit.stress_moment_n_m[0] = C64::ONE;
    let a = values(compact.evaluate(&qf).unwrap());
    let b = values(stress.evaluate(&[unit]).unwrap());
    let cancel = -a[0] / b[0];
    let channels = |stress_sign: f64| {
        let mut result = [C64::ZERO; 14];
        result[0] = qf[0].volume_flow_m3_s;
        result[8] = cancel.scale(stress_sign);
        result
    };
    let observe = |signs: [f64; 2]| {
        let spectrum = estimate(
            &record(&[(7, signs.map(channels))]),
            WelchWindow::Rectangular,
        );
        compact
            .evaluate_welch_bin_with_quadrupoles(&stress, &spectrum, 7, limits().coherence, || {
                false
            })
            .unwrap()
    };
    let destructive = observe([1.0, 1.0]);
    let constructive = observe([-1.0, -1.0]);
    let decorrelated = observe([1.0, -1.0]);
    let combined = core::array::from_fn(|i| a[i] + b[i] * cancel);
    moments(&destructive, &[combined]);
    assert!(destructive.mean_square_pressure_pa2 < 1e-20 * constructive.mean_square_pressure_pa2);
    assert!(
        destructive
            .velocity_covariance_m2_s2
            .iter()
            .enumerate()
            .map(|(i, row)| row[i].re)
            .sum::<f64>()
            > 1e-12
    );
    close(
        C64::from_re(decorrelated.mean_square_pressure_pa2),
        C64::from_re(0.5 * constructive.mean_square_pressure_pa2),
    );
}

#[test]
fn acquired_stress_clock_and_geometry_survive_wall_preparation_and_stress_only_bands() {
    let positions = [STRESS[0], COMPACT[1]];
    let moment = sources(amplitudes(0)).1;
    let source = [
        moment,
        QuadrupolePhasor {
            stress_moment_n_m: moment.stress_moment_n_m.map(|z| z * C64::new(-0.5, 0.7)),
        },
    ];
    let clock = QuadrupoleSampleClock::new(0.125, DT, 0.0).unwrap();
    let mut stream = QuadrupoleWelchStream::new(
        &positions,
        clock,
        N,
        config(WelchWindow::Rectangular),
        QuadrupoleRecordingBudget {
            maximum_sources: 2,
            welch: fft_budget(),
        },
        || false,
    )
    .unwrap();
    for j in 0..N {
        let angle = TAU * 7.0 * j as f64 / N as f64;
        let sample = source.map(|s| {
            s.stress_moment_n_m
                .map(|a| a.re * angle.cos() + a.im * angle.sin())
        });
        stream
            .push_sample(&positions, clock.time_at(j).unwrap(), &sample, || false)
            .unwrap();
    }
    let recorded = stream.finish().unwrap();
    assert_eq!(recorded.clock(), clock);
    let prepared = recorded
        .prepare_impedance_wall_total_observer(wall(None), RECEIVER, 7, limits().quadrupole, || {
            false
        })
        .unwrap();
    let mut analytic = [C64::ZERO; 4];
    for (position, stress) in positions.into_iter().zip(source) {
        let mut image = position;
        image[2] = -image[2];
        let mut image_stress = stress;
        image_stress.stress_moment_n_m[4] = -image_stress.stress_moment_n_m[4];
        image_stress.stress_moment_n_m[5] = -image_stress.stress_moment_n_m[5];
        let direct = green(
            prepared.frequency_hz(),
            position,
            CompactSourcePhasor::default(),
            stress,
        );
        let returned = green(
            prepared.frequency_hz(),
            image,
            CompactSourcePhasor::default(),
            image_stress,
        );
        for i in 0..4 {
            analytic[i] = analytic[i] + direct[i] + returned[i];
        }
    }
    moments(
        &recorded
            .evaluate_impedance_wall_bin(&prepared, 7, limits().coherence, || false)
            .unwrap(),
        &[analytic],
    );
    let band = observe_welch_band(
        recorded.spectrum(),
        7..8,
        &[],
        &positions,
        RECEIVER,
        Kind::Total,
        |_| Ok(wall(None)),
        ImpedanceWallWelchBandLimits {
            maximum_total_evaluations: 1,
            ..limits()
        },
        || false,
    )
    .unwrap();
    assert_eq!(band.evaluations, 0);
    moments(&band.statistics, &[analytic]);
    let reordered = Stress::new_total(
        wall(None),
        prepared.frequency_hz(),
        RECEIVER,
        &[positions[1], positions[0]],
        limits().quadrupole,
        &mut || false,
    )
    .unwrap();
    assert!(
        recorded
            .evaluate_impedance_wall_bin(&reordered, 7, limits().coherence, || false)
            .is_err()
    );
    assert!(
        recorded
            .evaluate_impedance_wall_bin(&prepared, 6, limits().coherence, || false)
            .is_err()
    );
}

#[test]
fn mixed_banks_require_the_same_wall_receiver_frequency_field_and_joint_record() {
    let data = record(&[(7, [amplitudes(0), amplitudes(1)])]);
    let spectrum = estimate(&data, WelchWindow::Rectangular);
    let f = spectrum.frequency_hz(7).unwrap();
    let (compact, stress) = observers(wall(None), f, Kind::Total);
    let shifted =
        PlanarImpedanceWall::new(wall(None).fluid(), [0.0, 0.0, 0.01], [0.0, 0.0, 1.0], None)
            .unwrap();
    let other_fluid = PlanarImpedanceWall::new(
        FluidAcousticState::try_new(1.3, C).unwrap(),
        [0.0; 3],
        [0.0, 0.0, 1.0],
        None,
    )
    .unwrap();
    for (plane, frequency, receiver, kind) in [
        (wall(Some(C64::ZERO)), f, RECEIVER, Kind::Total),
        (shifted, f, RECEIVER, Kind::Total),
        (other_fluid, f, RECEIVER, Kind::Total),
        (wall(None), f.next_up(), RECEIVER, Kind::Total),
        (wall(None), f, [0.17, 0.04, -0.1], Kind::Total),
        (wall(None), f, RECEIVER, Kind::Reflected),
    ] {
        let changed = match kind {
            Kind::Reflected => Stress::new_reflected(
                plane,
                frequency,
                receiver,
                &STRESS,
                limits().quadrupole,
                &mut || false,
            ),
            Kind::Total => Stress::new_total(
                plane,
                frequency,
                receiver,
                &STRESS,
                limits().quadrupole,
                &mut || false,
            ),
        }
        .unwrap();
        assert!(
            compact
                .evaluate_welch_bin_with_quadrupoles(
                    &changed,
                    &spectrum,
                    7,
                    limits().coherence,
                    || false
                )
                .is_err()
        );
    }
    assert!(
        compact
            .evaluate_welch_bin(&spectrum, 7, limits().coherence, || false)
            .is_err()
    );
    assert!(
        stress
            .evaluate_welch_bin(&spectrum, 7, limits().coherence, || false)
            .is_err()
    );
    for budget in [
        CoherenceBudget {
            max_modes: 1,
            max_source_mode_pairs: 6,
        },
        CoherenceBudget {
            max_modes: 2,
            max_source_mode_pairs: 5,
        },
    ] {
        assert!(
            compact
                .evaluate_welch_bin_with_quadrupoles(&stress, &spectrum, 7, budget, || false)
                .is_err()
        );
    }
    assert!(matches!(
        compact.evaluate_welch_bin_with_quadrupoles(
            &stress,
            &spectrum,
            7,
            limits().coherence,
            || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
}

#[test]
fn bands_admit_aggregate_work_before_material_calls_and_refuse_geometry_changes_atomically() {
    let spectrum = estimate(
        &record(&[(7, [amplitudes(0), amplitudes(1)])]),
        WelchWindow::Rectangular,
    );
    let calls = std::cell::Cell::new(0);
    for bins in [0..1, 32..33, 7..7, 5..12] {
        assert!(
            observe_welch_band(
                &spectrum,
                bins,
                &COMPACT,
                &STRESS,
                RECEIVER,
                Kind::Total,
                |_| {
                    calls.set(calls.get() + 1);
                    Ok(wall(None))
                },
                limits(),
                || false
            )
            .is_err()
        );
    }
    for budget in [
        CoherenceBudget {
            max_modes: 3,
            max_source_mode_pairs: 12,
        },
        CoherenceBudget {
            max_modes: 4,
            max_source_mode_pairs: 11,
        },
    ] {
        assert!(
            observe_welch_band(
                &spectrum,
                7..9,
                &COMPACT,
                &STRESS,
                RECEIVER,
                Kind::Total,
                |_| {
                    calls.set(calls.get() + 1);
                    Ok(wall(None))
                },
                ImpedanceWallWelchBandLimits {
                    coherence: budget,
                    ..limits()
                },
                || false
            )
            .is_err()
        );
    }
    assert_eq!(calls.get(), 0);
    let changed = observe_welch_band(
        &spectrum,
        7..9,
        &COMPACT,
        &STRESS,
        RECEIVER,
        Kind::Total,
        |_| {
            calls.set(calls.get() + 1);
            PlanarImpedanceWall::new(
                wall(None).fluid(),
                [0.0, 0.0, if calls.get() == 1 { 0.0 } else { 0.01 }],
                [0.0, 0.0, 1.0],
                None,
            )
        },
        limits(),
        || false,
    );
    assert!(changed.is_err());
    assert_eq!(calls.get(), 2);
    calls.set(0);
    let cancelled = observe_welch_band(
        &spectrum,
        7..9,
        &COMPACT,
        &STRESS,
        RECEIVER,
        Kind::Total,
        |_| {
            calls.set(calls.get() + 1);
            Ok(wall(None))
        },
        limits(),
        || calls.get() == 2,
    );
    assert!(matches!(cancelled, Err(AeroacError::Cancelled { .. })));
    assert_eq!(calls.get(), 2);
}
