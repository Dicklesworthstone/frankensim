use super::*;
use crate::network::field::{
    NetworkFieldExposureLimits, NetworkFieldFactorStatistics, NetworkFieldFrequency,
    NetworkFieldSpectralFactor, NetworkFieldSpectralLine, NetworkFieldSpectrumBudget,
    NetworkFieldSpectrumRequest, PreparedNetworkFieldSpectrum,
};
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn spectrum(
    net: &DuctNetwork,
    gas: &GasState,
    omegas: &[f64],
    sources: &[PortSource],
    loads: &[PortTermination],
    points: &[NetworkFieldProbe],
) -> PreparedNetworkFieldSpectrum {
    let frequencies: Vec<_> = omegas
        .iter()
        .map(|&omega| NetworkFieldFrequency {
            omega,
            loss: LossModel::Lossless,
            sources,
            source_bank: None,
            loads,
            coupled_load: None,
        })
        .collect();
    PreparedNetworkFieldSpectrum::new(
        net,
        gas,
        &NetworkFieldSpectrumRequest {
            source_clock: "same-network-emission-clock",
            probes: points,
            frequencies: &frequencies,
            field_budget: NetworkFieldBudget::default(),
            spectrum_budget: NetworkFieldSpectrumBudget {
                maximum_frequencies: 8,
                maximum_retained_transfer_entries: 4096,
            },
        },
        || false,
    )
    .unwrap()
}

fn exposure_limits() -> NetworkFieldExposureLimits {
    NetworkFieldExposureLimits {
        maximum_lines: 64,
        maximum_factors: 8,
        maximum_transfer_products: 16384,
        maximum_probe_frequency_pairs: 16384,
        maximum_workspace_bytes: 4_000_000,
    }
}

fn spectral_line<'a>(
    frequency_index: usize,
    excitation: &'a [C64],
) -> NetworkFieldSpectralLine<'a> {
    NetworkFieldSpectralLine {
        frequency_index,
        excitation,
    }
}

fn raw(point: NetworkFieldMoments) -> [f64; 5] {
    [
        point.mean_pressure_pa,
        point.mean_volume_flow_m3_s,
        point.mean_square_pressure_pa2,
        point.mean_square_volume_flow_m6_s2,
        point.mean_power_w,
    ]
}

// Direct pressure-Dirichlet cylinder solution; no production propagation,
// loaded network, factor contraction or temporal moment is used here.
fn cylinder_field(gas: &GasState, length: f64, x: f64, omega: f64, ends: [C64; 2]) -> [C64; 2] {
    let k = omega / gas.sound_speed;
    let z = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.008_f64.powi(2));
    let p = (ends[0].scale((k * (length - x)).sin()) + ends[1].scale((k * x).sin()))
        .scale(1.0 / (k * length).sin());
    let q = (ends[0].scale((k * (length - x)).cos()) - ends[1].scale((k * x).cos()))
        * C64::new(0.0, 1.0 / z / (k * length).sin());
    [p, q]
}

/// G1: native fields from TWO physical frequencies, a deterministic mean,
/// phase-locked real amplitude noise and a jointly proper complex factor.
#[test]
fn spectral_fields_match_direct_cylinder_and_eight_state_time_integration() {
    let gas = gas();
    let length = 0.37;
    let net = cylinder_net(length);
    let points = probes(length, &[0.0, 0.31, 0.73, 1.0]);
    let sources = [
        pressure_source("right", C64::ZERO),
        pressure_source("left", C64::ZERO),
    ];
    let omegas = [900.0, 1700.0];
    let owner = spectrum(&net, &gas, &omegas, &sources, &[], &points);
    assert_eq!(owner.frequency_count(), 2);
    assert_eq!(owner.source_coordinates()[0].label, "left");
    assert_eq!(owner.source_coordinates()[1].label, "right");
    let amplitudes = [
        [
            [C64::new(0.7, 0.3), C64::new(-0.2, 0.4)],
            [C64::new(0.2, -0.1), C64::new(0.6, 0.2)],
        ],
        [
            [C64::new(0.3, -0.4), C64::new(0.5, 0.1)],
            [C64::new(-0.1, 0.8), C64::new(0.2, -0.3)],
        ],
        [
            [C64::new(-0.2, 0.9), C64::new(0.6, 0.4)],
            [C64::new(0.4, 0.3), C64::new(-0.7, 0.2)],
        ],
    ];
    let banks: [[NetworkFieldSpectralLine<'_>; 2]; 3] = core::array::from_fn(|b| {
        [
            spectral_line(0, &amplitudes[b][0]),
            spectral_line(1, &amplitudes[b][1]),
        ]
    });
    let factors = [
        NetworkFieldSpectralFactor {
            lines: &banks[1],
            statistics: NetworkFieldFactorStatistics::Real,
        },
        NetworkFieldSpectralFactor {
            lines: &banks[2],
            statistics: NetworkFieldFactorStatistics::ProperComplex,
        },
    ];
    let prepared = owner
        .prepare_observation(&banks[0], &factors, exposure_limits(), || false)
        .unwrap();
    let phases = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)];
    let impedance = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.008_f64.powi(2));
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.00017, 0.0031, profile).unwrap();
        let output = prepared.observe_window(window, || false).unwrap();
        assert_eq!(output.source_clock(), "same-network-emission-clock");
        assert_eq!(output.locations(), owner.locations());
        assert_eq!(output.window(), window);
        for (r, point) in output.moments().iter().copied().enumerate() {
            let x = owner.locations()[r].axial_position_m;
            let fields: [[[C64; 2]; 2]; 3] = core::array::from_fn(|b| {
                core::array::from_fn(|f| {
                    cylinder_field(&gas, length, x, omegas[f], amplitudes[b][f])
                })
            });
            let mut expected = [0.0; 5];
            let n = 4096;
            for state in 0..8 {
                let latent = [
                    C64::ONE,
                    C64::from_re(if state < 4 { -1.0 } else { 1.0 }),
                    phases[state % 4],
                ];
                for step in 0..=n {
                    let u = step as f64 / n as f64;
                    let time = window.start_s() + window.duration_s() * u;
                    let simpson = if step == 0 || step == n {
                        1.0
                    } else if step % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    let shutter = if profile == HarmonicWindowProfile::Uniform {
                        1.0
                    } else {
                        2.0 * (1.0 - (2.0 * u - 1.0).abs())
                    };
                    let mut value = [0.0; 2];
                    for b in 0..3 {
                        for f in 0..2 {
                            let (s, c) = (omegas[f] * time).sin_cos();
                            for component in 0..2 {
                                let a = fields[b][f][component] * latent[b];
                                value[component] += a.re * c + a.im * s;
                            }
                        }
                    }
                    let [p, q] = value;
                    for (sum, term) in expected.iter_mut().zip([p, q, p * p, q * q, p * q]) {
                        *sum += simpson * shutter * term / (24.0 * n as f64);
                    }
                }
            }
            let scale = [
                4.0,
                4.0 / impedance,
                8.0,
                8.0 / impedance.powi(2),
                8.0 / impedance,
            ];
            for ((got, want), scale) in raw(point).into_iter().zip(expected).zip(scale) {
                assert!(
                    (got - want).abs() < 3e-9 * scale,
                    "probe {r}: {got:e} != {want:e}"
                );
            }
            let area = owner.locations()[r].area_m2;
            assert_eq!(point.mean_intensity_w_m2, point.mean_power_w / area);
            assert_eq!(
                point.mean_particle_velocity_m_s,
                point.mean_volume_flow_m3_s / area
            );
        }
    }
}

/// G1: retaining two close tones changes the measured interior beat. Adding
/// separately averaged line powers would lose this term.
#[test]
fn spectral_real_and_proper_factors_keep_different_subperiod_physics() {
    let gas = gas();
    let net = cylinder_net(0.2);
    let points = probes(0.2, &[0.0, 0.4]);
    let sources = [
        pressure_source("left", C64::ZERO),
        pressure_source("right", C64::ZERO),
    ];
    let owner = spectrum(&net, &gas, &[900.0, 901.0], &sources, &[], &points);
    let a = [C64::new(0.0, 1.0), C64::ZERO];
    let one = [spectral_line(0, &a)];
    let real = [NetworkFieldSpectralFactor {
        lines: &one,
        statistics: NetworkFieldFactorStatistics::Real,
    }];
    let proper = [NetworkFieldSpectralFactor {
        lines: &one,
        statistics: NetworkFieldFactorStatistics::ProperComplex,
    }];
    let window = HarmonicWindow::new(0.0, 1e-8).unwrap();
    let r = owner
        .prepare_observation(&[], &real, exposure_limits(), || false)
        .unwrap();
    let p = owner
        .prepare_observation(&[], &proper, exposure_limits(), || false)
        .unwrap();
    let r = r.observe_window(window, || false).unwrap();
    let p = p.observe_window(window, || false).unwrap();
    assert_eq!(r.moments()[0].mean_pressure_pa, 0.0);
    assert_eq!(p.moments()[0].mean_pressure_pa, 0.0);
    let expected = (900.0 * 1e-8_f64).powi(2) / 3.0;
    assert!((r.moments()[0].mean_square_pressure_pa2 / expected - 1.0).abs() < 1e-6);
    assert!((p.moments()[0].mean_square_pressure_pa2 - 0.5).abs() < 1e-10);

    let positive = [C64::ONE, C64::ZERO];
    let negative = [-C64::ONE, C64::ZERO];
    let beat = [spectral_line(0, &positive), spectral_line(1, &negative)];
    let correlated = [NetworkFieldSpectralFactor {
        lines: &beat,
        statistics: NetworkFieldFactorStatistics::ProperComplex,
    }];
    let bank = owner
        .prepare_observation(&[], &correlated, exposure_limits(), || false)
        .unwrap();
    let w = HarmonicWindow::new(-0.05, 0.1).unwrap();
    let observed = bank.observe_window(w, || false).unwrap();
    let expected = 1.0 - 0.05_f64.sin() / 0.05;
    assert!((observed.moments()[0].mean_square_pressure_pa2 - expected).abs() < 5e-10);
    assert!(observed.moments()[0].mean_square_pressure_pa2 < 0.001);
}

/// G3: exact-frequency coalescing retains a weak survivor; translated source
/// phase and observation clock preserve the physical mean and power.
#[test]
fn spectral_grouping_and_clock_translation_preserve_interference() {
    let gas = gas();
    let net = cylinder_net(0.27);
    let points = probes(0.27, &[0.0, 0.35, 0.8]);
    let sources = [
        pressure_source("left", C64::ZERO),
        pressure_source("right", C64::ZERO),
    ];
    let omega = [700.0, 1100.0];
    let owner = spectrum(&net, &gas, &omega, &sources, &[], &points);
    let large = [C64::from_re(2.0_f64.powi(54)), C64::ZERO];
    let small = [C64::new(0.25, 0.125), C64::ZERO];
    let opposite = [-large[0], C64::ZERO];
    let supplied = [
        spectral_line(0, &large),
        spectral_line(0, &small),
        spectral_line(0, &opposite),
    ];
    let simple = [spectral_line(0, &small)];
    let grouped = owner
        .prepare_observation(&supplied, &[], exposure_limits(), || false)
        .unwrap();
    let plain = owner
        .prepare_observation(&simple, &[], exposure_limits(), || false)
        .unwrap();
    assert_eq!(grouped.probe_frequency_pairs(), points.len());
    let window = HarmonicWindow::new(0.0003, 0.0021).unwrap();
    let a = grouped.observe_window(window, || false).unwrap();
    let b = plain.observe_window(window, || false).unwrap();
    for (x, y) in a.moments().iter().zip(b.moments()) {
        for (x, y) in raw(*x).into_iter().zip(raw(*y)) {
            assert!(
                (x - y).abs() <= 1e-12 * y.abs().max(1e-18),
                "{x:e} != {y:e}"
            );
        }
    }

    let amplitudes = [
        [C64::new(0.7, 0.4), C64::new(-0.2, 0.1)],
        [C64::new(-0.3, 0.8), C64::new(0.1, 0.5)],
    ];
    let shift = 0.0037;
    let shifted: [[C64; 2]; 2] = core::array::from_fn(|f| {
        let (s, c) = (omega[f] * shift).sin_cos();
        amplitudes[f].map(|a| a * C64::new(c, s))
    });
    let original = [
        spectral_line(0, &amplitudes[0]),
        spectral_line(1, &amplitudes[1]),
    ];
    let translated = [spectral_line(0, &shifted[0]), spectral_line(1, &shifted[1])];
    let a = owner
        .prepare_observation(&original, &[], exposure_limits(), || false)
        .unwrap();
    let b = owner
        .prepare_observation(&translated, &[], exposure_limits(), || false)
        .unwrap();
    let a = a.observe_window(window, || false).unwrap();
    let b = b
        .observe_window(
            HarmonicWindow::new(window.start_s() + shift, window.duration_s()).unwrap(),
            || false,
        )
        .unwrap();
    for (x, y) in a.moments().iter().zip(b.moments()) {
        for (x, y) in raw(*x).into_iter().zip(raw(*y)) {
            assert!((x - y).abs() < 1e-10 * x.abs().max(y.abs()).max(1e-17));
        }
    }
}

/// G4: source-coordinate changes and ambiguous repeated model frequencies
/// cannot silently become one cross-frequency physical source.
#[test]
fn spectral_preparation_refuses_inconsistent_sources_frequencies_and_capacity() {
    let gas = gas();
    let net = cylinder_net(0.2);
    let points = probes(0.2, &[0.0, 0.5]);
    let both = [
        pressure_source("left", C64::ZERO),
        pressure_source("right", C64::ZERO),
    ];
    let left = [pressure_source("left", C64::ZERO)];
    let closed = [PortTermination {
        label: "right".into(),
        load: NetworkLoad::Closed,
    }];
    let frequencies = [
        NetworkFieldFrequency {
            omega: 700.0,
            loss: LossModel::Lossless,
            sources: &both,
            source_bank: None,
            loads: &[],
            coupled_load: None,
        },
        NetworkFieldFrequency {
            omega: 900.0,
            loss: LossModel::Lossless,
            sources: &left,
            source_bank: None,
            loads: &closed,
            coupled_load: None,
        },
    ];
    let budget = NetworkFieldSpectrumBudget {
        maximum_frequencies: 8,
        maximum_retained_transfer_entries: 1024,
    };
    let request = |f| NetworkFieldSpectrumRequest {
        source_clock: "shared",
        probes: &points,
        frequencies: f,
        field_budget: NetworkFieldBudget::default(),
        spectrum_budget: budget,
    };
    assert!(
        PreparedNetworkFieldSpectrum::new(&net, &gas, &request(&frequencies), || false).is_err()
    );
    let duplicate = [frequencies[0], frequencies[0]];
    assert!(PreparedNetworkFieldSpectrum::new(&net, &gas, &request(&duplicate), || false).is_err());
    let mut valid = duplicate;
    valid[1].omega = 900.0;
    let mut query = request(&valid);
    query.spectrum_budget.maximum_frequencies = 1;
    assert!(PreparedNetworkFieldSpectrum::new(&net, &gas, &query, || false).is_err());
    query.spectrum_budget.maximum_frequencies = 2;
    query.spectrum_budget.maximum_retained_transfer_entries = 15;
    assert!(PreparedNetworkFieldSpectrum::new(&net, &gas, &query, || false).is_err());
    query.spectrum_budget.maximum_retained_transfer_entries = 16;
    assert!(PreparedNetworkFieldSpectrum::new(&net, &gas, &query, || false).is_ok());
    assert!(matches!(
        PreparedNetworkFieldSpectrum::new(&net, &gas, &query, || true),
        Err(DuctError::Cancelled { .. })
    ));
}

/// G4/G5: bounded atomic observation and replay, including explicitly dark
/// fields. Indexed lines never accept independently prepared foreign owners.
#[test]
fn spectral_observation_budget_cancellation_darkness_and_replay() {
    let gas = gas();
    let net = cylinder_net(0.2);
    let points = probes(0.2, &[0.0, 0.6]);
    let sources = [
        pressure_source("left", C64::ZERO),
        pressure_source("right", C64::ZERO),
    ];
    let owner = spectrum(&net, &gas, &[700.0, 1000.0], &sources, &[], &points);
    let amplitudes = [C64::ONE, C64::ZERO];
    let line = [spectral_line(0, &amplitudes), spectral_line(1, &amplitudes)];
    let wrong = [spectral_line(2, &amplitudes)];
    assert!(
        owner
            .prepare_observation(&wrong, &[], exposure_limits(), || false)
            .is_err()
    );
    let bad_shape = [spectral_line(0, &[C64::ONE])];
    assert!(
        owner
            .prepare_observation(&bad_shape, &[], exposure_limits(), || false)
            .is_err()
    );
    for index in 0..4 {
        let mut limit = exposure_limits();
        match index {
            0 => limit.maximum_lines = 1,
            1 => limit.maximum_transfer_products = 0,
            2 => limit.maximum_probe_frequency_pairs = 0,
            _ => limit.maximum_workspace_bytes = 0,
        }
        assert!(
            owner
                .prepare_observation(&line, &[], limit, || false)
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        owner.prepare_observation(&line, &[], exposure_limits(), || {
            polls += 1;
            polls == 5
        }),
        Err(DuctError::Cancelled { .. })
    ));
    let ready = owner
        .prepare_observation(&line, &[], exposure_limits(), || false)
        .unwrap();
    let window = HarmonicWindow::new(0.0001, 0.004).unwrap();
    assert!(matches!(
        ready.observe_window(window, || true),
        Err(DuctError::Cancelled { .. })
    ));
    assert_eq!(
        ready.observe_window(window, || false).unwrap().moments(),
        ready.observe_window(window, || false).unwrap().moments()
    );
    let empty = owner
        .prepare_observation(&[], &[], exposure_limits(), || false)
        .unwrap();
    let dark = empty.observe_window(window, || false).unwrap();
    assert_eq!(
        dark.moments(),
        vec![NetworkFieldMoments::default(); points.len()]
    );
    assert_eq!(dark.locations(), owner.locations());
}
