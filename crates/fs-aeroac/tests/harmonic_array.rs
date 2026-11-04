//! Spatial and processed acoustic observations on the original common clock.
use fs_aeroac::AeroacError;
use fs_aeroac::observation::{
    HarmonicEnsembleBudget, HarmonicFactorStatistics, HarmonicFieldFactor, HarmonicFieldLine,
    PreparedHarmonicEnsemble,
    array::{
        AcousticCrossMoments, HarmonicArrayBudget, HarmonicArrayReceiver, HarmonicFactorCoordinate,
        HarmonicReadoutBudget, HarmonicReadoutTap, PreparedHarmonicArray,
    },
};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::TAU;

const COORDINATES: [HarmonicFactorCoordinate<'static>; 2] = [
    HarmonicFactorCoordinate {
        id: "shared amplitude fluctuation",
        statistics: HarmonicFactorStatistics::Real,
    },
    HarmonicFactorCoordinate {
        id: "shared proper fluctuation",
        statistics: HarmonicFactorStatistics::ProperComplex,
    },
];
const ARRAY_WORK: HarmonicArrayBudget = HarmonicArrayBudget {
    max_receivers: 3,
    max_factors: 2,
    max_total_lines: 24,
    max_line_pairs: 144,
};
const READOUT_WORK: HarmonicReadoutBudget = HarmonicReadoutBudget {
    max_taps: 8,
    max_total_lines: 100,
    max_line_pairs: 4000,
};

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.2, 340.0).unwrap()
}

fn line(frequency_hz: f64, pressure: C64, velocity: [C64; 3]) -> HarmonicFieldLine<3> {
    HarmonicFieldLine {
        frequency_hz,
        pressure_pa: pressure,
        particle_velocity_m_s: velocity,
    }
}

fn ensemble(
    mean: &[HarmonicFieldLine<3>],
    factors: &[HarmonicFieldFactor<'_, 3>],
) -> PreparedHarmonicEnsemble<3> {
    PreparedHarmonicEnsemble::new(
        mean,
        factors,
        HarmonicEnsembleBudget {
            max_factors: 2,
            max_total_lines: 12,
            max_line_pairs: 144,
        },
    )
    .unwrap()
}

fn receiver<'a>(x: f64, ensemble: &'a PreparedHarmonicEnsemble<3>) -> HarmonicArrayReceiver<'a, 3> {
    HarmonicArrayReceiver {
        position_m: [x, 0.0, 0.0],
        acoustic: medium(),
        ensemble,
    }
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * scale.max(1e-18),
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

fn matrix(value: AcousticCrossMoments<3>) -> [[f64; 4]; 4] {
    let mut result = [[0.0; 4]; 4];
    result[0][0] = value.pressure_pressure_pa2;
    for i in 0..3 {
        result[0][i + 1] = value.pressure_velocity_w_m2[i];
        result[i + 1][0] = value.velocity_pressure_w_m2[i];
        for j in 0..3 {
            result[i + 1][j + 1] = value.velocity_velocity_m2_s2[i][j];
        }
    }
    result
}

fn sample(lines: &[HarmonicFieldLine<3>], t: f64, coefficient: C64) -> [f64; 4] {
    let mut out = [0.0; 4];
    for line in lines {
        let phase = C64::new(
            (TAU * line.frequency_hz * t).cos(),
            -(TAU * line.frequency_hz * t).sin(),
        );
        for (i, amplitude) in core::iter::once(line.pressure_pa)
            .chain(line.particle_velocity_m_s)
            .enumerate()
        {
            out[i] += (amplitude * coefficient * phase).re;
        }
    }
    out
}

#[test]
fn spatial_real_and_proper_moments_match_independent_ensemble_time_quadrature() {
    // Distinct axis amplitudes/phases make the off-diagonal tensor independent
    // checks rather than repetitions of one velocity component.
    let axes = |re, im| {
        [
            C64::new(re, im),
            C64::new(-0.4 * im, 1.6 * re),
            C64::new(0.7 * re - 0.3 * im, -1.1 * im),
        ]
    };
    let mean = [
        [
            line(173.0, C64::new(0.8, -0.2), axes(0.002, 0.001)),
            line(177.0, C64::new(-0.3, 0.1), axes(-0.001, 0.002)),
        ],
        [
            line(173.0, C64::new(-0.2, 0.6), axes(0.001, -0.002)),
            line(177.0, C64::new(0.4, -0.5), axes(0.002, 0.003)),
        ],
    ];
    let real = [
        [line(173.0, C64::new(0.3, 0.2), axes(-0.002, 0.001))],
        [line(181.0, C64::new(-0.1, 0.4), axes(0.003, -0.001))],
    ];
    let proper = [
        [line(177.0, C64::new(0.2, -0.3), axes(0.001, 0.002))],
        [line(181.0, C64::new(0.4, 0.1), axes(-0.002, 0.001))],
    ];
    let ensembles: [_; 2] = core::array::from_fn(|r| {
        ensemble(
            &mean[r],
            &[
                HarmonicFieldFactor {
                    lines: &real[r],
                    statistics: HarmonicFactorStatistics::Real,
                },
                HarmonicFieldFactor {
                    lines: &proper[r],
                    statistics: HarmonicFactorStatistics::ProperComplex,
                },
            ],
        )
    });
    let array = PreparedHarmonicArray::new(
        "accepted physical clock",
        &COORDINATES,
        &[receiver(0.0, &ensembles[0]), receiver(0.2, &ensembles[1])],
        ARRAY_WORK,
    )
    .unwrap();
    let start = 0.00013;
    let duration = 0.0037;
    let window = HarmonicWindow::new(start, duration).unwrap();
    let actual = array.observe_pair([0, 1], window).unwrap();
    let reversed = array.observe_pair([1, 0], window).unwrap();
    let mut products = [[0.0; 4]; 4];
    let mut averaged = [[0.0; 4]; 4];
    let mut means = [[0.0; 4]; 2];
    let samples = 4096;
    for sign in [-1.0, 1.0] {
        for z in [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)] {
            let mut averages = [[0.0; 4]; 2];
            for j in 0..=samples {
                let t = start + duration * j as f64 / samples as f64;
                let weight = if j == 0 || j == samples {
                    1.0
                } else if j % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                let weight = weight / (3.0 * samples as f64);
                let fields: [[f64; 4]; 2] = core::array::from_fn(|r| {
                    let a = sample(&mean[r], t, C64::ONE);
                    let b = sample(&real[r], t, C64::from_re(sign));
                    let c = sample(&proper[r], t, z);
                    core::array::from_fn(|i| a[i] + b[i] + c[i])
                });
                for i in 0..4 {
                    for r in 0..2 {
                        averages[r][i] += weight * fields[r][i];
                    }
                    for k in 0..4 {
                        products[i][k] += weight * fields[0][i] * fields[1][k] / 8.0;
                    }
                }
            }
            for i in 0..4 {
                for r in 0..2 {
                    means[r][i] += averages[r][i] / 8.0;
                }
                for j in 0..4 {
                    averaged[i][j] += averages[0][i] * averages[1][j] / 8.0;
                }
            }
        }
    }
    let second = matrix(actual.second_moments);
    let covariance = matrix(actual.averaged_covariance);
    let reverse = matrix(reversed.second_moments);
    for i in 0..4 {
        for j in 0..4 {
            let scale = if i == 0 && j == 0 {
                1.0
            } else if i == 0 || j == 0 {
                0.01
            } else {
                1e-4
            };
            close(second[i][j], products[i][j], scale);
            close(
                covariance[i][j],
                averaged[i][j] - means[0][i] * means[1][j],
                scale,
            );
            close(reverse[j][i], second[i][j], scale);
        }
    }
    let local = array.observe_receiver(0, window, || false).unwrap();
    let self_pair = array.observe_pair([0, 0], window).unwrap().second_moments;
    close(
        self_pair.pressure_pressure_pa2,
        local.mean_square_pressure_pa2,
        1.0,
    );
    for i in 0..3 {
        close(
            self_pair.pressure_velocity_w_m2[i],
            local.mean_intensity_w_m2[i],
            0.01,
        );
    }
    assert_ne!(
        actual.second_moments.pressure_pressure_pa2,
        actual.averaged_covariance.pressure_pressure_pa2
    );
}

#[test]
fn causal_delays_and_directional_microphones_use_actual_velocity_and_original_material() {
    let f = 223.0;
    let omega = TAU * f;
    let distance = 0.31;
    let delay = distance / medium().sound_speed_m_s();
    let impedance = medium().characteristic_impedance_pa_s_m();
    let phase = C64::new((omega * delay).cos(), (omega * delay).sin());
    // A coherent wave travelling toward -x; the front axis points toward +x.
    let near = [line(
        f,
        C64::ONE,
        [C64::from_re(-1.0 / impedance), C64::ZERO, C64::ZERO],
    )];
    let far = [line(
        f,
        phase,
        [-phase.scale(1.0 / impedance), C64::ZERO, C64::ZERO],
    )];
    let fields = [ensemble(&near, &[]), ensemble(&far, &[])];
    let array = PreparedHarmonicArray::new(
        "wave clock",
        &[],
        &[receiver(0.0, &fields[0]), receiver(-distance, &fields[1])],
        ARRAY_WORK,
    )
    .unwrap();
    let cancel = array
        .prepare_readout(
            &[
                HarmonicReadoutTap::pressure(0, 1.0, delay),
                HarmonicReadoutTap::pressure(1, -1.0, 0.0),
            ],
            READOUT_WORK,
        )
        .unwrap();
    let window = HarmonicWindow::new(0.00031, 0.0023).unwrap();
    assert!(
        cancel
            .observe_window(window)
            .unwrap()
            .mean_square_pressure_pa2
            < 1e-29
    );
    let front = HarmonicReadoutTap {
        receiver: 0,
        gain: 1.0,
        delay_s: 0.0,
        pressure_fraction: 0.5,
        front_axis: [2.0, 0.0, 0.0],
    };
    let front = array.prepare_readout(&[front], READOUT_WORK).unwrap();
    let back = array
        .prepare_readout(
            &[HarmonicReadoutTap {
                receiver: 0,
                gain: 1.0,
                delay_s: 0.0,
                pressure_fraction: 0.5,
                front_axis: [-1.0, 0.0, 0.0],
            }],
            READOUT_WORK,
        )
        .unwrap();
    close(
        front
            .observe_window(window)
            .unwrap()
            .mean_square_pressure_pa2,
        fields[0]
            .observe_window(window)
            .unwrap()
            .mean_square_pressure_pa2,
        1.0,
    );
    assert!(
        back.observe_window(window)
            .unwrap()
            .mean_square_pressure_pa2
            < 1e-30
    );

    // At a pressure null a velocity-sensitive microphone still responds.
    let pressure_null = ensemble(
        &[line(
            f,
            C64::ZERO,
            [C64::new(0.0, 0.004), C64::ZERO, C64::ZERO],
        )],
        &[],
    );
    let array = PreparedHarmonicArray::new(
        "near field clock",
        &[],
        &[receiver(0.0, &pressure_null)],
        ARRAY_WORK,
    )
    .unwrap();
    let figure_eight = array
        .prepare_readout(
            &[HarmonicReadoutTap {
                receiver: 0,
                gain: 1.0,
                delay_s: 0.0,
                pressure_fraction: 0.0,
                front_axis: [1.0, 0.0, 0.0],
            }],
            READOUT_WORK,
        )
        .unwrap();
    assert!(
        figure_eight
            .observe_window(window)
            .unwrap()
            .mean_square_pressure_pa2
            > 0.1
    );
}

#[test]
fn readout_grouping_preserves_dd_cancellation_and_weak_fluctuation_covariance() {
    let f = 7.0;
    let mean_a = [
        line(f, C64::from_re(1e16), [C64::ZERO; 3]),
        line(f, C64::ONE, [C64::ZERO; 3]),
    ];
    let mean_b = [line(f, C64::from_re(-1e16), [C64::ZERO; 3])];
    let noise = [line(f, C64::from_re(3e-4), [C64::ZERO; 3])];
    let factor = [HarmonicFieldFactor {
        lines: &noise,
        statistics: HarmonicFactorStatistics::Real,
    }];
    let fields = [ensemble(&mean_a, &factor), ensemble(&mean_b, &factor)];
    let coordinates = [COORDINATES[0]];
    let array = PreparedHarmonicArray::new(
        "shared quiet clock",
        &coordinates,
        &[receiver(0.0, &fields[0]), receiver(0.1, &fields[1])],
        ARRAY_WORK,
    )
    .unwrap();
    let channel = array
        .prepare_readout(
            &[
                HarmonicReadoutTap::pressure(0, 1.0, 0.0),
                HarmonicReadoutTap::pressure(1, 1.0, 0.0),
            ],
            READOUT_WORK,
        )
        .unwrap();
    let window = HarmonicWindow::new(0.005, 0.013).unwrap();
    let reference = ensemble(
        &[line(f, C64::ONE, [C64::ZERO; 3])],
        &[HarmonicFieldFactor {
            lines: &[line(f, C64::from_re(6e-4), [C64::ZERO; 3])],
            statistics: HarmonicFactorStatistics::Real,
        }],
    )
    .observe_window(window)
    .unwrap();
    let actual = channel.observe_window(window).unwrap();
    close(actual.mean_pressure_pa, reference.mean_pressure_pa, 1.0);
    close(
        actual.mean_square_pressure_pa2,
        reference.mean_square_pressure_pa2,
        1.0,
    );
    let mean_cosine = ((TAU * f * 0.018).sin() - (TAU * f * 0.005).sin()) / (TAU * f * 0.013);
    let pair = array.observe_pair([0, 1], window).unwrap();
    close(
        pair.averaged_covariance.pressure_pressure_pa2,
        (3e-4 * mean_cosine).powi(2),
        1e-7,
    );
    let processed = channel.observe_pair(&channel, window, 4, || false).unwrap();
    close(
        processed.averaged_covariance_pa2,
        (6e-4 * mean_cosine).powi(2),
        1e-6,
    );
}

#[test]
fn readout_cross_covariance_preserves_shared_proper_factors_and_empty_channels() {
    let left_lines = [line(113.0, C64::new(0.8, -0.4), [C64::ZERO; 3])];
    let right_lines = [line(127.0, C64::new(-0.2, 0.9), [C64::ZERO; 3])];
    let left = ensemble(
        &[],
        &[HarmonicFieldFactor {
            lines: &left_lines,
            statistics: HarmonicFactorStatistics::ProperComplex,
        }],
    );
    let right = ensemble(
        &[],
        &[HarmonicFieldFactor {
            lines: &right_lines,
            statistics: HarmonicFactorStatistics::ProperComplex,
        }],
    );
    let coordinates = [COORDINATES[1]];
    let array = PreparedHarmonicArray::new(
        "joint phase clock",
        &coordinates,
        &[receiver(0.0, &left), receiver(0.2, &right)],
        ARRAY_WORK,
    )
    .unwrap();
    let a = array
        .prepare_readout(&[HarmonicReadoutTap::pressure(0, 1.0, 0.0)], READOUT_WORK)
        .unwrap();
    let b = array
        .prepare_readout(&[HarmonicReadoutTap::pressure(1, 1.0, 0.0)], READOUT_WORK)
        .unwrap();
    let dark = array.prepare_readout(&[], READOUT_WORK).unwrap();
    let window = HarmonicWindow::new(0.001, 0.004).unwrap();
    let physical = array.observe_pair([0, 1], window).unwrap();
    let processed = a.observe_pair(&b, window, 1, || false).unwrap();
    close(
        processed.mean_product_pa2,
        physical.second_moments.pressure_pressure_pa2,
        1.0,
    );
    close(
        processed.averaged_covariance_pa2,
        physical.averaged_covariance.pressure_pressure_pa2,
        1.0,
    );
    assert_eq!(
        dark.observe_window(window)
            .unwrap()
            .mean_square_pressure_pa2,
        0.0
    );
    assert_eq!(
        dark.observe_pair(&a, window, 0, || false)
            .unwrap()
            .averaged_covariance_pa2,
        0.0
    );
}

#[test]
fn array_factor_clock_work_and_physical_readout_refusals_are_atomic() {
    let tone = [line(137.0, C64::ONE, [C64::ZERO; 3])];
    let e = ensemble(&tone, &[]);
    let receivers = [receiver(0.0, &e), receiver(0.1, &e)];
    assert!(PreparedHarmonicArray::new("", &[], &receivers, ARRAY_WORK).is_err());
    assert!(PreparedHarmonicArray::new("clock", &COORDINATES, &receivers, ARRAY_WORK).is_err());
    assert!(
        PreparedHarmonicArray::new(
            "clock",
            &[],
            &receivers,
            HarmonicArrayBudget {
                max_line_pairs: 3,
                ..ARRAY_WORK
            }
        )
        .is_err()
    );
    let array = PreparedHarmonicArray::new("clock", &[], &receivers, ARRAY_WORK).unwrap();
    let window = HarmonicWindow::new(0.0, 0.001).unwrap();
    assert!(array.observe_pair([0, 2], window).is_err());
    assert!(matches!(
        array.observe_pair_cancellable([0, 1], window, || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        array
            .prepare_readout(&[HarmonicReadoutTap::pressure(0, 1.0, -0.1)], READOUT_WORK)
            .is_err()
    );
    assert!(
        array
            .prepare_readout(
                &[HarmonicReadoutTap {
                    pressure_fraction: 0.5,
                    ..HarmonicReadoutTap::pressure(0, 1.0, 0.0)
                }],
                READOUT_WORK
            )
            .is_err()
    );
    assert!(
        array
            .prepare_readout(
                &[HarmonicReadoutTap::pressure(0, 1.0, 0.0); 2],
                HarmonicReadoutBudget {
                    max_total_lines: 1,
                    ..READOUT_WORK
                }
            )
            .is_err()
    );
    let taps = [
        HarmonicReadoutTap::pressure(0, 1.0, 0.0),
        HarmonicReadoutTap::pressure(1, 0.5, 0.001),
    ];
    let mut polls = 0;
    assert!(matches!(
        array.prepare_readout_cancellable(&taps, READOUT_WORK, || {
            polls += 1;
            polls == 7
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    let a = array.prepare_readout(&taps, READOUT_WORK).unwrap();
    let different_clock =
        PreparedHarmonicArray::new("different clock", &[], &receivers, ARRAY_WORK).unwrap();
    let b = different_clock
        .prepare_readout(&taps, READOUT_WORK)
        .unwrap();
    assert!(a.observe_pair(&b, window, 10, || false).is_err());
    assert!(a.observe_pair(&a, window, 0, || false).is_err());
    let first = a.observe_window(window).unwrap();
    let second = a.observe_window(window).unwrap();
    assert_eq!(
        first.mean_square_pressure_pa2,
        second.mean_square_pressure_pa2
    );
}
