//! Cross-receiver coherency through the same positive temporal factors.

use fs_aeroac::AeroacError;
use fs_aeroac::observation::{
    HarmonicEnsembleBudget, HarmonicFactorStatistics, HarmonicFieldFactor, HarmonicFieldLine,
    PreparedHarmonicEnsemble,
    array::{
        AcousticCrossMoments, HarmonicArrayBudget, HarmonicArrayReceiver, HarmonicFactorCoordinate,
        PreparedHarmonicArray,
    },
};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;

type Line = HarmonicFieldLine<2>;
const BUDGET: HarmonicArrayBudget = HarmonicArrayBudget {
    max_receivers: 2,
    max_factors: 2,
    max_total_lines: 18,
    max_line_pairs: 108,
};
const MIXED: [HarmonicFactorCoordinate<'static>; 2] = [
    HarmonicFactorCoordinate {
        id: "real",
        statistics: HarmonicFactorStatistics::Real,
    },
    HarmonicFactorCoordinate {
        id: "proper",
        statistics: HarmonicFactorStatistics::ProperComplex,
    },
];
const PROPER: [HarmonicFactorCoordinate<'static>; 1] = [HarmonicFactorCoordinate {
    id: "shared",
    statistics: HarmonicFactorStatistics::ProperComplex,
}];

fn ensemble(mean: &[Line], factors: &[HarmonicFieldFactor<'_, 2>]) -> PreparedHarmonicEnsemble<2> {
    PreparedHarmonicEnsemble::new(
        mean,
        factors,
        HarmonicEnsembleBudget {
            max_factors: 2,
            max_total_lines: 9,
            max_line_pairs: 27,
        },
    )
    .unwrap()
}

fn receiver(x: f64, ensemble: &PreparedHarmonicEnsemble<2>) -> HarmonicArrayReceiver<'_, 2> {
    HarmonicArrayReceiver {
        position_m: [x, 0.0],
        acoustic: FluidAcousticState::try_new(1.2, 343.0).unwrap(),
        ensemble,
    }
}

fn lines(offset: f64, phase: C64, velocity: [f64; 2]) -> [Line; 3] {
    // Authored order is deliberately not sorted.
    let frequencies = [1002.0, 1000.0, 1001.0];
    let amplitudes = [1.0, 1.0, -2.0];
    std::array::from_fn(|i| {
        let pressure = phase.scale(amplitudes[i]);
        Line {
            frequency_hz: frequencies[i] + offset,
            pressure_pa: pressure,
            particle_velocity_m_s: velocity.map(|v| pressure.scale(v)),
        }
    })
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite() && expected != 0.0);
    assert!(
        (actual / expected - 1.0).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

fn matrix(value: AcousticCrossMoments<2>) -> [[f64; 3]; 3] {
    [
        [
            value.pressure_pressure_pa2,
            value.pressure_velocity_w_m2[0],
            value.pressure_velocity_w_m2[1],
        ],
        [
            value.velocity_pressure_w_m2[0],
            value.velocity_velocity_m2_s2[0][0],
            value.velocity_velocity_m2_s2[0][1],
        ],
        [
            value.velocity_pressure_w_m2[1],
            value.velocity_velocity_m2_s2[1][0],
            value.velocity_velocity_m2_s2[1][1],
        ],
    ]
}

#[test]
fn short_band_self_pairs_equal_the_original_local_full_moments() {
    let source = lines(0.0, C64::ONE, [0.25, -0.5]);
    let quadrature = lines(0.0, C64::new(0.0, 0.5), [0.25, -0.5]);
    let field = ensemble(
        &source,
        &[
            HarmonicFieldFactor {
                lines: &quadrature,
                statistics: HarmonicFactorStatistics::Real,
            },
            HarmonicFieldFactor {
                lines: &source,
                statistics: HarmonicFactorStatistics::ProperComplex,
            },
        ],
    );
    let array = PreparedHarmonicArray::new(
        "physical source clock",
        &MIXED,
        &[receiver(0.0, &field)],
        BUDGET,
    )
    .unwrap();
    for duration in [1e-6, 1e-9, 1e-12] {
        let window = HarmonicWindow::new(-0.5 * duration, duration).unwrap();
        let local = array.observe_receiver(0, window, || false).unwrap();
        let pair = array.observe_pair([0, 0], window).unwrap();
        assert!(local.mean_square_pressure_pa2 > 0.0);
        assert_eq!(
            pair.second_moments.pressure_pressure_pa2,
            local.mean_square_pressure_pa2
        );
        assert_eq!(
            pair.second_moments.pressure_velocity_w_m2,
            local.mean_intensity_w_m2
        );
        assert_eq!(
            pair.second_moments.velocity_pressure_w_m2,
            local.mean_intensity_w_m2
        );
        assert_eq!(
            pair.second_moments.velocity_velocity_m2_s2,
            local.velocity_second_moment_m2_s2
        );
        assert!(pair.averaged_covariance.pressure_pressure_pa2 > 0.0);
        assert!(pair.averaged_covariance.pressure_pressure_pa2 <= local.mean_square_pressure_pa2);
    }
    let window = HarmonicWindow::new(-0.5e-9, 1e-9).unwrap();
    let mut polls = 0;
    assert!(matches!(
        array.observe_pair_cancellable([0, 0], window, || {
            polls += 1;
            polls == 8
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        array
            .observe_pair([0, 0], window)
            .unwrap()
            .second_moments
            .pressure_pressure_pa2
            > 0.0
    );
}

#[test]
fn different_spectral_supports_retain_complex_source_phase_and_cross_reciprocity() {
    let a_scales = [0.25, -0.5];
    let b_scales = [-0.75, 0.125];
    let source_a = lines(0.0, C64::ONE, a_scales);
    // B(t)=(-1/2+i/4)*exp(-i*2pi*t)*A(t). Different supports need
    // a common four-frequency basis, with the missing fields exactly zero.
    let source_b = lines(1.0, C64::new(-0.5, 0.25), b_scales);
    let a = ensemble(
        &[],
        &[HarmonicFieldFactor {
            lines: &source_a,
            statistics: HarmonicFactorStatistics::ProperComplex,
        }],
    );
    let b = ensemble(
        &[],
        &[HarmonicFieldFactor {
            lines: &source_b,
            statistics: HarmonicFactorStatistics::ProperComplex,
        }],
    );
    let array = PreparedHarmonicArray::new(
        "physical source clock",
        &PROPER,
        &[receiver(0.0, &a), receiver(1.0, &b)],
        BUDGET,
    )
    .unwrap();
    let duration = 1e-9;
    let window = HarmonicWindow::new(-0.5 * duration, duration).unwrap();
    let power = (core::f64::consts::TAU * duration).powi(4) / 160.0;
    let cross = -0.5 * power;
    let actual = array.observe_pair([0, 1], window).unwrap();
    let reversed = array.observe_pair([1, 0], window).unwrap();
    relative(actual.second_moments.pressure_pressure_pa2, cross, 3e-9);
    for i in 0..2 {
        relative(
            actual.second_moments.pressure_velocity_w_m2[i],
            b_scales[i] * cross,
            3e-9,
        );
        relative(
            actual.second_moments.velocity_pressure_w_m2[i],
            a_scales[i] * cross,
            3e-9,
        );
        for j in 0..2 {
            relative(
                actual.second_moments.velocity_velocity_m2_s2[i][j],
                a_scales[i] * b_scales[j] * cross,
                3e-9,
            );
        }
    }
    for (ab, ba) in [
        (
            matrix(actual.second_moments),
            matrix(reversed.second_moments),
        ),
        (
            matrix(actual.averaged_covariance),
            matrix(reversed.averaged_covariance),
        ),
    ] {
        for i in 0..3 {
            for j in 0..3 {
                relative(ab[i][j], ba[j][i], 3e-9);
            }
        }
    }
    let pa = array
        .observe_receiver(0, window, || false)
        .unwrap()
        .mean_square_pressure_pa2;
    let pb = array
        .observe_receiver(1, window, || false)
        .unwrap()
        .mean_square_pressure_pa2;
    relative(pa, power, 3e-9);
    relative(pb, 0.3125 * power, 3e-9);
    assert!(actual.second_moments.pressure_pressure_pa2.abs() <= (pa * pb).sqrt());

    // The same fields with independent source coordinates have no cross moment.
    let coordinates = [
        HarmonicFactorCoordinate {
            id: "A",
            statistics: HarmonicFactorStatistics::ProperComplex,
        },
        HarmonicFactorCoordinate {
            id: "B",
            statistics: HarmonicFactorStatistics::ProperComplex,
        },
    ];
    let independent_a = ensemble(
        &[],
        &[
            HarmonicFieldFactor {
                lines: &source_a,
                statistics: HarmonicFactorStatistics::ProperComplex,
            },
            HarmonicFieldFactor {
                lines: &[],
                statistics: HarmonicFactorStatistics::ProperComplex,
            },
        ],
    );
    let independent_b = ensemble(
        &[],
        &[
            HarmonicFieldFactor {
                lines: &[],
                statistics: HarmonicFactorStatistics::ProperComplex,
            },
            HarmonicFieldFactor {
                lines: &source_b,
                statistics: HarmonicFactorStatistics::ProperComplex,
            },
        ],
    );
    let independent = PreparedHarmonicArray::new(
        "physical source clock",
        &coordinates,
        &[receiver(0.0, &independent_a), receiver(1.0, &independent_b)],
        BUDGET,
    )
    .unwrap()
    .observe_pair([0, 1], window)
    .unwrap();
    assert_eq!(matrix(independent.second_moments), [[0.0; 3]; 3]);
    assert_eq!(matrix(independent.averaged_covariance), [[0.0; 3]; 3]);
}

#[test]
fn averaged_proper_covariance_uses_the_physical_carrier_and_finite_peak_normalization() {
    let source = [1024.0, 1025.0, 1026.0].map(|frequency_hz| Line {
        frequency_hz,
        pressure_pa: if frequency_hz == 1024.0 {
            C64::from_re(1.5e154)
        } else {
            C64::ZERO
        },
        particle_velocity_m_s: [C64::ZERO; 2],
    });
    let field = ensemble(
        &[],
        &[HarmonicFieldFactor {
            lines: &source,
            statistics: HarmonicFactorStatistics::ProperComplex,
        }],
    );
    let array = PreparedHarmonicArray::new(
        "physical source clock",
        &PROPER,
        &[receiver(0.0, &field)],
        BUDGET,
    )
    .unwrap();
    let short = array
        .observe_pair([0, 0], HarmonicWindow::new(0.0, 1e-9).unwrap())
        .unwrap();
    relative(short.second_moments.pressure_pressure_pa2, 1.125e308, 2e-14);
    relative(
        short.averaged_covariance.pressure_pressure_pa2,
        1.125e308,
        2e-11,
    );
    // A whole original-carrier cycle averages this pure harmonic to zero.
    // The analytic band's relative carrier must not become the microphone clock.
    let full = array
        .observe_pair([0, 0], HarmonicWindow::new(0.0, 1.0 / 1024.0).unwrap())
        .unwrap();
    relative(full.second_moments.pressure_pressure_pa2, 1.125e308, 2e-14);
    assert!(full.averaged_covariance.pressure_pressure_pa2 < 1e-50 * 1.125e308);
}
