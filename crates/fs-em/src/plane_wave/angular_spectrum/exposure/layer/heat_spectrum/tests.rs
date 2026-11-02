//! G1/G3 real-time reservoir oracles; no sampled random-phase approximation.
use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;
use fs_material::field_storage::dielectric::DebyePole;

type Electric = [[C64; 3]; 2];
const OMEGA: [f64; 2] = [2.0, 2.3];
const SIGMA: f64 = 0.3;

fn poles() -> [DebyePole; 2] {
    [(0.7, 0.4), (0.2, 0.9)].map(|(delta, tau)| DebyePole {
        permittivity_increment_f_m: delta,
        relaxation_time_s: tau,
        capacitance_f: 1.0,
        resistance_ohm: tau,
    })
}
fn fixture() -> [Electric; 3] {
    let mean = [
        [C64::new(0.8, 0.2), C64::new(-0.1, 0.3), C64::new(0.2, -0.4)],
        [C64::new(-0.3, 0.6), C64::new(0.4, 0.1), C64::new(-0.2, 0.1)],
    ];
    let scaled =
        |coefficients: [C64; 2]| core::array::from_fn(|i| mean[i].map(|e| e * coefficients[i]));
    [
        mean,
        scaled([C64::new(0.3, -0.2), C64::new(-0.1, 0.5)]),
        scaled([C64::new(-0.4, 0.1), C64::new(0.6, -0.2)]),
    ]
}
fn limits() -> HarmonicProductSpectrumLimits {
    HarmonicProductSpectrumLimits {
        maximum_products: 512,
        maximum_generated_terms: 1024,
        maximum_workspace_bytes: 1024 * core::mem::size_of::<HarmonicProductTone>(),
    }
}
fn spectrum(
    input: &[(HarmonicProductLaw, Electric)],
    channel: AngularSpectrumExposureHeatChannel,
) -> PreparedHarmonicProductSpectrum {
    let (first, count) = reservoir_range(channel, SIGMA, poles().len()).unwrap();
    let mut fields = Vec::new();
    let mut banks = Vec::new();
    for &(law, electric) in input {
        for reservoir in first..first + count {
            let start = fields.len();
            let root = Dd::from_f64(det::sqrt(SIGMA));
            let left_scale = if reservoir == 0 {
                Dd::from_f64(SIGMA) / (root * root)
            } else {
                Dd::ONE
            };
            for i in 0..2 {
                let transfer = if reservoir == 0 {
                    C64::from_re(root.to_f64())
                } else {
                    PoleTransfer::new(poles()[reservoir - 1], OMEGA[i])
                        .unwrap()
                        .heat
                };
                fields.push(HeatLine {
                    omega: OMEGA[i],
                    field: electric[i].map(|e| filtered_negative_time_parts(e, transfer)),
                });
            }
            banks.push(HeatBank {
                start,
                end: fields.len(),
                law,
                left_scale,
            });
        }
    }
    prepare_products(&banks, &fields, limits(), &mut || false).unwrap()
}

// Direct original ODE, independent of the normalized reservoir filter:
// P=Delta*E/(1+i omega tau), q= sigma E_real² + sum tau/Delta Pdot_real².
fn real_heat(input: &[Electric; 3], coefficients: [C64; 3], time: f64) -> [f64; 4] {
    let mut e = [0.0; 3];
    let mut pd = [[0.0; 3]; 2];
    for (bank, z) in input.iter().zip(coefficients) {
        for i in 0..2 {
            let (s, c) = (OMEGA[i] * time).sin_cos();
            let real = |a: C64| a.re * c - a.im * s;
            for axis in 0..3 {
                let amplitude = bank[i][axis] * z;
                e[axis] += real(amplitude);
                for (p, pole) in poles().iter().enumerate() {
                    let polarization = amplitude.scale(pole.permittivity_increment_f_m)
                        / C64::new(1.0, OMEGA[i] * pole.relaxation_time_s);
                    let derivative = C64::new(-polarization.im, polarization.re).scale(OMEGA[i]);
                    pd[p][axis] += real(derivative);
                }
            }
        }
    }
    let ohmic = SIGMA * e.iter().map(|v| v * v).sum::<f64>();
    let q: [f64; 2] = core::array::from_fn(|p| {
        let pole = poles()[p];
        (pole.relaxation_time_s / pole.permittivity_increment_f_m)
            * pd[p].iter().map(|v| v * v).sum::<f64>()
    });
    [ohmic, q[0], q[1], ohmic + q[0] + q[1]]
}
fn expected(input: &[Electric; 3], time: f64) -> [f64; 4] {
    let mut sum = [0.0; 4];
    for sign in [-1.0, 1.0] {
        for z in [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)] {
            for (sum, value) in
                sum.iter_mut()
                    .zip(real_heat(input, [C64::ONE, C64::from_re(sign), z], time))
            {
                *sum += value / 8.0;
            }
        }
    }
    sum
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 3e-12 * (1.0 + b.abs()), "{a:e} != {b:e}");
}

#[test]
fn g1_mixed_reservoir_spectra_match_exact_states_and_both_window_profiles() {
    use AngularSpectrumExposureHeatChannel as Channel;
    let input = fixture();
    let banks = [
        (HarmonicProductLaw::Real, input[0]),
        (HarmonicProductLaw::Real, input[1]),
        (HarmonicProductLaw::ProperComplex, input[2]),
    ];
    let spectra = [
        Channel::Ohmic,
        Channel::DebyePole { index: 0 },
        Channel::DebyePole { index: 1 },
        Channel::Total,
    ]
    .map(|channel| spectrum(&banks, channel));
    let polarization = spectrum(&banks, Channel::Polarization);
    for time in [-0.3, 0.0, 0.17, 0.49, 1.7] {
        let reference = expected(&input, time);
        for (actual, expected) in spectra.iter().zip(reference) {
            close(actual.evaluate(time).unwrap(), expected);
        }
        close(
            polarization.evaluate(time).unwrap(),
            reference[1] + reference[2],
        );
    }
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.43, profile).unwrap();
        let n = 2048;
        let mut quadrature = [0.0; 4];
        for k in 0..=n {
            let u = k as f64 / n as f64;
            let quadrature_weight = if k == 0 || k == n {
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
            for (sum, value) in quadrature.iter_mut().zip(expected(&input, 0.17 + 0.43 * u)) {
                *sum += quadrature_weight * shutter * value / (3.0 * n as f64);
            }
        }
        for (spectrum, expected) in spectra.iter().zip(quadrature) {
            close(spectrum.mean(window).unwrap(), expected);
        }
    }
}

#[test]
fn g3_proper_real_quadratures_native_carrier_beats_and_refusals() {
    use AngularSpectrumExposureHeatChannel as Channel;
    let field = fixture()[2];
    let proper = spectrum(
        &[(HarmonicProductLaw::ProperComplex, field)],
        Channel::Total,
    );
    let quadratures = [
        C64::from_re(core::f64::consts::FRAC_1_SQRT_2),
        C64::new(0.0, core::f64::consts::FRAC_1_SQRT_2),
    ]
    .map(|z| {
        (
            HarmonicProductLaw::Real,
            field.map(|line| line.map(|e| e * z)),
        )
    });
    let real = spectrum(&quadratures, Channel::Total);
    for time in [-0.17, 0.0, 0.23, 0.9] {
        close(proper.evaluate(time).unwrap(), real.evaluate(time).unwrap());
    }
    let carrier = 2.0_f64.powi(50);
    let fields = [
        (carrier + 2.0, C64::ONE),
        (carrier + 2.25, C64::new(0.0, 1.0)),
    ]
    .map(|(omega, e)| HeatLine {
        omega,
        field: [
            filtered_negative_time_parts(e, C64::ONE),
            [Dd::ZERO; 2],
            [Dd::ZERO; 2],
        ],
    });
    let banks = [HeatBank {
        start: 0,
        end: 2,
        law: HarmonicProductLaw::ProperComplex,
        left_scale: Dd::ONE,
    }];
    let beat = prepare_products(&banks, &fields, limits(), &mut || false).unwrap();
    for time in [-0.2, 0.0, 0.17, 3.0] {
        close(beat.evaluate(time).unwrap(), 1.0 - (0.25_f64 * time).sin());
    }
    assert!(reservoir_range(Channel::DebyePole { index: 2 }, SIGMA, 2).is_err());
    assert_eq!(reservoir_range(Channel::Ohmic, 0.0, 2).unwrap(), (0, 0));
    assert!(
        prepare_products(
            &banks,
            &fields,
            HarmonicProductSpectrumLimits {
                maximum_products: 0,
                ..limits()
            },
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        prepare_products(&banks, &fields, limits(), &mut || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    close(beat.evaluate(0.17).unwrap(), 1.0 - (0.25_f64 * 0.17).sin());
}
