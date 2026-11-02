//! G1/G3/G4 spatial and temporal oracles for the actual contraction path.
use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;
use fs_material::field_storage::dielectric::DebyePole;

const OMEGA: [f64; 2] = [2.0, 2.3];
const SIGMA: f64 = 0.3;
const GRID: AngularSpectrumGrid = AngularSpectrumGrid {
    width: 2,
    height: 2,
    spacing_m: [0.25, 0.4],
};
const DEPTH: [AngularSpectrumHeatDepthSample; 3] = [
    AngularSpectrumHeatDepthSample {
        depth_m: 0.0,
        weight_m: 0.1,
    },
    AngularSpectrumHeatDepthSample {
        depth_m: 0.3,
        weight_m: 0.4,
    },
    AngularSpectrumHeatDepthSample {
        depth_m: 0.6,
        weight_m: 0.1,
    },
];
const AREAS: [f64; 4] = [0.1, 0.0, 0.05, 0.1];

fn poles() -> [DebyePole; 2] {
    [(0.7, 0.4), (0.2, 0.9)].map(|(delta, tau)| DebyePole {
        permittivity_increment_f_m: delta,
        relaxation_time_s: tau,
        capacitance_f: 1.0,
        resistance_ohm: tau,
    })
}
fn quadrature() -> AngularSpectrumHeatVolumeQuadrature<'static> {
    AngularSpectrumHeatVolumeQuadrature {
        depth_span_m: [0.0, 0.6],
        depth_samples: &DEPTH,
        transverse: AngularSpectrumHeatAreaQuadrature::CellAreasM2(&AREAS),
    }
}
fn electric(bank: usize, group: usize, point: usize, depth: f64) -> [C64; 3] {
    let base = [
        [C64::new(0.8, 0.2), C64::new(-0.1, 0.3), C64::new(0.2, -0.4)],
        [C64::new(-0.3, 0.6), C64::new(0.4, 0.1), C64::new(-0.2, 0.1)],
    ];
    let bank_scale = [
        [C64::ONE, C64::ONE],
        [C64::new(0.3, -0.2), C64::new(-0.1, 0.5)],
        [C64::new(-0.4, 0.1), C64::new(0.6, -0.2)],
    ];
    let phase = 0.4 * point as f64 + (0.7 + 0.3 * group as f64) * depth;
    let (s, c) = phase.sin_cos();
    let forward = C64::new(c, -s);
    let backward = C64::new(c, s);
    core::array::from_fn(|axis| {
        (base[group][axis] * forward
            + base[1 - group][axis] * backward.scale(0.2 + 0.1 * axis as f64))
            * bank_scale[bank][group]
    })
}
fn prepared(
    quadrature: AngularSpectrumHeatVolumeQuadrature<'_>,
    channel: AngularSpectrumExposureHeatChannel,
) -> (PreparedHarmonicProductSpectrum, usize, usize) {
    let volume = admit_volume(GRID, 0.6, quadrature, &mut || false).unwrap();
    let cache: Vec<Vec<PoleTransfer>> = OMEGA
        .iter()
        .map(|&omega| {
            poles()
                .into_iter()
                .map(|pole| PoleTransfer::new(pole, omega).unwrap())
                .collect()
        })
        .collect();
    let (first, count) = reservoir_range(channel, SIGMA, 2).unwrap();
    let mut terms = Vec::new();
    let mut field_calls = 0;
    for bank in 0..3 {
        let proper = bank == 2;
        contract_bank(
            &OMEGA,
            proper,
            quadrature,
            volume,
            |i, j, cancelled| {
                reservoir_pair(SIGMA, &cache[i], &cache[j], first, count, proper, cancelled)
            },
            |group, depth, _cancelled| {
                field_calls += 1;
                Ok((0..4)
                    .map(|point| {
                        let e = electric(bank, group, point, depth);
                        [e[0], e[1], e[2], C64::ZERO, C64::ZERO, C64::ZERO]
                    })
                    .collect())
            },
            &mut terms,
            &mut || false,
        )
        .unwrap();
    }
    let count = terms.len();
    let spectrum = PreparedHarmonicProductSpectrum::from_angular_terms(
        terms,
        HarmonicSpectrumTermLimits {
            maximum_terms: count,
            maximum_workspace_bytes: count * core::mem::size_of::<HarmonicProductTone>(),
        },
    )
    .unwrap();
    (spectrum, field_calls, count)
}

// Independent positive-time field/current calculation at each ORIGINAL
// quadrature sample, followed by weighted physical-volume integration.
fn reference(quadrature: AngularSpectrumHeatVolumeQuadrature<'_>, time: f64) -> [f64; 4] {
    let mut total = [0.0; 4];
    for sample in quadrature.depth_samples {
        for point in 0..4 {
            let area = match quadrature.transverse {
                AngularSpectrumHeatAreaQuadrature::PeriodicCell => 0.25 * 0.4,
                AngularSpectrumHeatAreaQuadrature::CellAreasM2(areas) => areas[point],
            };
            for sign in [-1.0, 1.0] {
                for z in [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)] {
                    let coefficients = [C64::ONE, C64::from_re(sign), z];
                    let mut e = [0.0; 3];
                    let mut pd = [[0.0; 3]; 2];
                    for (bank, coefficient) in coefficients.into_iter().enumerate() {
                        for group in 0..2 {
                            let amplitude = electric(bank, group, point, sample.depth_m);
                            let omega = OMEGA[group];
                            let (s, c) = (omega * time).sin_cos();
                            let real = |a: C64| a.re * c - a.im * s;
                            for axis in 0..3 {
                                let a = amplitude[axis] * coefficient;
                                e[axis] += real(a);
                                for (index, pole) in poles().iter().enumerate() {
                                    let p = a.scale(pole.permittivity_increment_f_m)
                                        / C64::new(1.0, omega * pole.relaxation_time_s);
                                    pd[index][axis] += real(C64::new(-p.im, p.re).scale(omega));
                                }
                            }
                        }
                    }
                    let ohmic = SIGMA * e.iter().map(|x| x * x).sum::<f64>();
                    let polarization: [f64; 2] = core::array::from_fn(|p| {
                        let pole = poles()[p];
                        (pole.relaxation_time_s / pole.permittivity_increment_f_m)
                            * pd[p].iter().map(|x| x * x).sum::<f64>()
                    });
                    let values = [
                        ohmic,
                        polarization[0],
                        polarization[1],
                        ohmic + polarization[0] + polarization[1],
                    ];
                    for (sum, value) in total.iter_mut().zip(values) {
                        *sum += value * area * sample.weight_m / 8.0;
                    }
                }
            }
        }
    }
    total
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 4e-12 * (1.0 + b.abs()), "{a:e} != {b:e}");
}

#[test]
fn g1_volume_heat_matches_independent_real_fields_and_exact_mixed_states() {
    use AngularSpectrumExposureHeatChannel as Channel;
    let q = quadrature();
    let spectra = [
        Channel::Ohmic,
        Channel::DebyePole { index: 0 },
        Channel::DebyePole { index: 1 },
        Channel::Total,
    ]
    .map(|channel| prepared(q, channel));
    for (_, calls, terms) in &spectra {
        assert_eq!(*calls, 3 * 2 * DEPTH.len());
        assert_eq!(*terms, 15);
    }
    for time in [-0.4, 0.0, 0.17, 0.83] {
        let expected = reference(q, time);
        for ((spectrum, _, _), expected) in spectra.iter().zip(expected) {
            close(spectrum.evaluate(time).unwrap(), expected);
        }
    }
    let polarization = prepared(q, Channel::Polarization).0;
    close(
        polarization.evaluate(0.17).unwrap(),
        reference(q, 0.17)[1] + reference(q, 0.17)[2],
    );
    // Duplicating nodes with half weights changes reconstruction work but
    // never the physical integral or retained frequency-term count.
    let doubled: Vec<_> = DEPTH
        .iter()
        .flat_map(|s| {
            [AngularSpectrumHeatDepthSample {
                weight_m: 0.5 * s.weight_m,
                ..*s
            }; 2]
        })
        .collect();
    let q2 = AngularSpectrumHeatVolumeQuadrature {
        depth_samples: &doubled,
        ..q
    };
    let (same, calls, terms) = prepared(q2, Channel::Total);
    assert_eq!(calls, 2 * spectra[3].1);
    assert_eq!(terms, spectra[3].2);
    for time in [-0.4, 0.0, 0.17, 0.83] {
        close(
            same.evaluate(time).unwrap(),
            spectra[3].0.evaluate(time).unwrap(),
        );
    }
}

#[test]
fn g1_integrated_heat_means_match_both_physical_shutters_and_full_cell_support() {
    for transverse in [
        AngularSpectrumHeatAreaQuadrature::PeriodicCell,
        AngularSpectrumHeatAreaQuadrature::CellAreasM2(&AREAS),
    ] {
        let q = AngularSpectrumHeatVolumeQuadrature {
            transverse,
            ..quadrature()
        };
        let (spectrum, _, _) = prepared(q, AngularSpectrumExposureHeatChannel::Total);
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let window = HarmonicWindow::with_profile(0.17, 0.43, profile).unwrap();
            let n = 512;
            let mut expected = 0.0;
            for k in 0..=n {
                let u = k as f64 / n as f64;
                let weight = if k == 0 || k == n {
                    1.0
                } else if k % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                let shutter = if profile == HarmonicWindowProfile::Uniform {
                    1.0
                } else {
                    2.0 * (1.0 - (2.0 * u - 1.0).abs())
                };
                expected += weight * shutter * reference(q, 0.17 + 0.43 * u)[3] / (3.0 * n as f64);
            }
            close(spectrum.mean(window).unwrap(), expected);
        }
    }
}

#[test]
fn g4_physical_volume_admission_preserves_weights_and_refuses_invalid_support() {
    let q = quadrature();
    let volume = admit_volume(GRID, 0.6, q, &mut || false).unwrap();
    close(volume.area_m2.to_f64(), 0.25);
    close(volume.volume_m3.to_f64(), 0.15);
    let mut adjusted = DEPTH;
    adjusted[1].weight_m += 2e-15;
    let accepted = admit_volume(
        GRID,
        0.6,
        AngularSpectrumHeatVolumeQuadrature {
            depth_samples: &adjusted,
            ..q
        },
        &mut || false,
    )
    .unwrap();
    assert!(accepted.depth_weight_m.to_f64() > 0.6);
    assert_eq!(
        accepted.depth_weight_m,
        adjusted
            .iter()
            .fold(Dd::ZERO, |sum, s| sum + Dd::from_f64(s.weight_m))
    );
    for invalid in [
        AngularSpectrumHeatDepthSample {
            depth_m: -0.1,
            weight_m: 0.6,
        },
        AngularSpectrumHeatDepthSample {
            depth_m: 0.7,
            weight_m: 0.6,
        },
        AngularSpectrumHeatDepthSample {
            depth_m: 0.3,
            weight_m: 0.0,
        },
        AngularSpectrumHeatDepthSample {
            depth_m: 0.3,
            weight_m: 0.5,
        },
        AngularSpectrumHeatDepthSample {
            depth_m: f64::NAN,
            weight_m: 0.6,
        },
    ] {
        assert!(
            admit_volume(
                GRID,
                0.6,
                AngularSpectrumHeatVolumeQuadrature {
                    depth_samples: &[invalid],
                    ..q
                },
                &mut || false
            )
            .is_err()
        );
    }
    for areas in [
        [0.1, 0.1, 0.1, 0.1001],
        [-0.1, 0.1, 0.1, 0.1],
        [f64::NAN, 0.0, 0.0, 0.0],
    ] {
        assert!(
            admit_volume(
                GRID,
                0.6,
                AngularSpectrumHeatVolumeQuadrature {
                    transverse: AngularSpectrumHeatAreaQuadrature::CellAreasM2(&areas),
                    ..q
                },
                &mut || false
            )
            .is_err()
        );
    }
    let dark = admit_volume(
        GRID,
        0.6,
        AngularSpectrumHeatVolumeQuadrature {
            transverse: AngularSpectrumHeatAreaQuadrature::CellAreasM2(&[0.0; 4]),
            ..q
        },
        &mut || false,
    )
    .unwrap();
    assert_eq!(dark.volume_m3.to_f64(), 0.0);
    assert!(matches!(
        admit_volume(GRID, 0.6, q, &mut || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    let mut terms = Vec::new();
    assert!(matches!(
        contract_bank(
            &OMEGA,
            false,
            q,
            volume,
            |_, _, _| Ok(PairIntegral::default()),
            |_, _, _| Ok(vec![[C64::ZERO; 6]; 4]),
            &mut terms,
            &mut || true
        ),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(terms.is_empty());
    assert!(
        contract_bank(
            &OMEGA,
            false,
            q,
            volume,
            |_, _, _| Ok(PairIntegral::default()),
            |_, _, _| Ok(vec![[C64::ZERO; 6]; 3]),
            &mut terms,
            &mut || false
        )
        .is_err()
    );
    assert!(terms.is_empty());
}

// This calls the actual spatial/reservoir contraction with independently
// authored physical units. The reference filters/rescales before multiplication
// using ordinary finite f64 values, as the existing local heat owner does.
fn single_balanced_heat(
    omega: f64,
    electric: C64,
    sigma: f64,
    pole: Option<DebyePole>,
    proper: bool,
    area_m2: f64,
) -> Result<PreparedHarmonicProductSpectrum, AngularSpectrumError> {
    let grid = AngularSpectrumGrid {
        width: 1,
        height: 1,
        spacing_m: [1.0, 1.0],
    };
    let depth = [AngularSpectrumHeatDepthSample {
        depth_m: 0.5,
        weight_m: 1.0,
    }];
    let areas = [area_m2];
    let quadrature = AngularSpectrumHeatVolumeQuadrature {
        depth_span_m: [0.0, 1.0],
        depth_samples: &depth,
        transverse: AngularSpectrumHeatAreaQuadrature::CellAreasM2(&areas),
    };
    let volume = admit_volume(grid, 1.0, quadrature, &mut || false)?;
    let response: Vec<_> = pole
        .into_iter()
        .map(|p| PoleTransfer::new(p, omega))
        .collect::<Result<_, _>>()?;
    let first = usize::from(!response.is_empty());
    let mut terms = Vec::new();
    contract_bank(
        &[omega],
        proper,
        quadrature,
        volume,
        |_, _, cancelled| reservoir_pair(sigma, &response, &response, first, 1, proper, cancelled),
        |_, _, _| {
            Ok(vec![[
                electric,
                C64::ZERO,
                C64::ZERO,
                C64::ZERO,
                C64::ZERO,
                C64::ZERO,
            ]])
        },
        &mut terms,
        &mut || false,
    )?;
    PreparedHarmonicProductSpectrum::from_angular_terms(
        terms,
        HarmonicSpectrumTermLimits {
            maximum_terms: if proper { 1 } else { 2 },
            maximum_workspace_bytes: 2 * core::mem::size_of::<HarmonicProductTone>(),
        },
    )
    .map_err(spectrum_error)
}

fn relative_heat(actual: f64, expected: f64) {
    assert!(
        expected.is_finite() && expected > 0.0,
        "bad reference {expected:e}"
    );
    assert!(
        actual.is_finite() && actual > 0.0,
        "lost finite heat {actual:e}, expected {expected:e}"
    );
    assert!(
        (actual / expected - 1.0).abs() < 3e-12,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_exponent_balancing_preserves_finite_ohmic_debye_and_volume_heat() {
    let omega = 2.0;
    for proper in [false, true] {
        // E² underflows in an unscaled dot, although sigma*E² is representable.
        let electric = C64::new(1e-200, -0.25e-200);
        let sigma = 1e200;
        let spectrum = single_balanced_heat(omega, electric, sigma, None, proper, 1.0).unwrap();
        for phase in [0.0_f64, 0.17, 0.83] {
            let (sin, cos) = phase.sin_cos();
            let e = electric.re * cos - electric.im * sin;
            let expected = if proper {
                0.5 * ((sigma * electric.re) * electric.re + (sigma * electric.im) * electric.im)
            } else {
                (sigma * e) * e
            };
            relative_heat(spectrum.evaluate(phase / omega).unwrap(), expected);
        }
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let window = HarmonicWindow::with_profile(0.13, 0.37, profile).unwrap();
            let left = [
                Dd::from_f64(sigma) * Dd::from_f64(electric.re),
                -(Dd::from_f64(sigma) * Dd::from_f64(electric.im)),
            ];
            let right = [Dd::from_f64(electric.re), -Dd::from_f64(electric.im)];
            let expected = if proper {
                window
                    .proper_product_kernel(omega, omega)
                    .unwrap()
                    .mean_real_product_parts_dd(left, right)
                    .unwrap()
            } else {
                window
                    .product_kernel(omega, omega)
                    .unwrap()
                    .mean_real_product_parts_dd(left, right)
                    .unwrap()
            };
            relative_heat(spectrum.mean(window).unwrap(), expected.to_f64());
        }

        // The inverse case: an unweighted spatial E² overflows, while the
        // positive physical cell area makes the final load finite.
        let large = C64::new(1e200, -0.25e200);
        let area = 1e-200;
        let spectrum = single_balanced_heat(omega, large, 1.0, None, proper, area).unwrap();
        for phase in [0.0_f64, 0.17, 0.83] {
            let (sin, cos) = phase.sin_cos();
            let e = large.re * cos - large.im * sin;
            let expected = if proper {
                0.5 * ((area * large.re) * large.re + (area * large.im) * large.im)
            } else {
                (area * e) * e
            };
            relative_heat(spectrum.evaluate(phase / omega).unwrap(), expected);
        }

        // Original PoleTransfer values are finite in both directions. Their
        // squares are not, but the complete filtered physical fields are.
        for (delta, tau, omega, electric) in [
            (1e-300, 1e-100, 2.0, C64::new(1e150, 0.3e150)),
            (1e300, 1e-100, 1e200, C64::new(1e-200, -0.25e-200)),
        ] {
            let pole = DebyePole {
                permittivity_increment_f_m: delta,
                relaxation_time_s: tau,
                capacitance_f: 1.0,
                resistance_ohm: tau,
            };
            let filtered = electric * PoleTransfer::new(pole, omega).unwrap().heat;
            assert!(filtered.re.is_finite() && filtered.im.is_finite());
            let spectrum =
                single_balanced_heat(omega, electric, 0.0, Some(pole), proper, 1.0).unwrap();
            for phase in [0.0_f64, 0.17, 0.83] {
                let (sin, cos) = phase.sin_cos();
                let response = filtered.re * cos - filtered.im * sin;
                let expected = if proper {
                    0.5 * filtered.norm_sq()
                } else {
                    response * response
                };
                relative_heat(spectrum.evaluate(phase / omega).unwrap(), expected);
            }
        }
    }
}

#[test]
fn g4_scaled_contraction_refuses_unrepresentable_final_watts_without_darkening() {
    for proper in [false, true] {
        assert!(single_balanced_heat(2.0, C64::from_re(1e-200), 1.0, None, proper, 1.0).is_err());
        assert!(single_balanced_heat(2.0, C64::from_re(1e200), 1.0, None, proper, 1.0).is_err());
        let zero = single_balanced_heat(2.0, C64::ZERO, 1.0, None, proper, 1.0).unwrap();
        assert_eq!(zero.evaluate(0.2).unwrap(), 0.0);
    }
}
