use super::*;
use fs_fft::harmonic_window::{AngularHarmonicSpectrumTerm, HarmonicSpectrumTermLimits};

fn spectrum(dc: f64, cosine: f64, sine: f64, omega: f64) -> PreparedHarmonicProductSpectrum {
    PreparedHarmonicProductSpectrum::from_angular_terms(
        [
            AngularHarmonicSpectrumTerm {
                angular_frequency_rad_s: Dd::ZERO,
                peak_amplitude: [Dd::from_f64(dc), Dd::ZERO],
            },
            AngularHarmonicSpectrumTerm {
                angular_frequency_rad_s: Dd::from_f64(omega),
                peak_amplitude: [Dd::from_f64(cosine), Dd::from_f64(sine)],
            },
        ],
        HarmonicSpectrumTermLimits {
            maximum_terms: 2,
            maximum_workspace_bytes: 4096,
        },
    )
    .unwrap()
}
fn controls() -> SpatialInterferenceHeatingControls {
    SpatialInterferenceHeatingControls {
        maximum_thermal_step_s: 1.0,
        maximum_source_ripple_j: 100.0,
    }
}
fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 3e-14, "{actual} != {expected}");
}

#[test]
fn cell_joules_follow_the_physical_interval_without_volume_scaling() {
    let start = 0.3_f64;
    let dt = 0.6;
    let sources = [spectrum(2.0, 0.7, 0.0, 2.3), spectrum(3.0, 0.0, 0.4, 1.7)];
    let value =
        integrate_cells(&sources, Dd::from_f64(start), dt, controls(), &mut || false).unwrap();
    let first = 2.0 * dt + 0.7 / 2.3 * ((2.3 * (start + dt)).sin() - (2.3 * start).sin());
    let second = 3.0 * dt + 0.4 / 1.7 * ((1.7 * start).cos() - (1.7 * (start + dt)).cos());
    near(value.cell_energy_j[0], first);
    near(value.cell_energy_j[1], second);
    near(value.total_j.to_f64(), first + second);
    near(value.end_s.to_f64(), start + dt);
}

#[test]
fn continuation_preserves_sub_ulp_source_epochs_and_partitioned_energy() {
    let start = Dd::from_pair(9_007_199_254_740_992.0, 0.125);
    let sources = [spectrum(2.0, 0.8, 0.3, 0.7)];
    let first = integrate_cells(&sources, start, 0.25, controls(), &mut || false).unwrap();
    let second = integrate_cells(&sources, first.end_s, 0.25, controls(), &mut || false).unwrap();
    let whole = integrate_cells(&sources, start, 0.5, controls(), &mut || false).unwrap();
    assert_eq!(first.end_s.to_f64(), start.to_f64());
    near((first.end_s - start).to_f64(), 0.25);
    near((second.end_s - start).to_f64(), 0.5);
    near(
        (first.total_j + second.total_j).to_f64(),
        whole.total_j.to_f64(),
    );
    let reset = integrate_cells(&sources, start, 0.25, controls(), &mut || false).unwrap();
    assert!((second.total_j - reset.total_j).to_f64().abs() > 1e-3);
}

#[test]
fn opposite_cell_beats_cannot_hide_local_source_ripple() {
    let sources = [spectrum(2.0, 1.0, 0.0, 3.0), spectrum(2.0, -1.0, 0.0, 3.0)];
    let accepted = integrate_cells(&sources, Dd::ZERO, 0.7, controls(), &mut || false).unwrap();
    near(accepted.total_j.to_f64(), 2.8);
    assert!(accepted.ripple_j > 0.0);
    let refused = integrate_cells(
        &sources,
        Dd::ZERO,
        0.7,
        SpatialInterferenceHeatingControls {
            maximum_source_ripple_j: 0.0,
            ..controls()
        },
        &mut || false,
    );
    assert!(matches!(
        refused,
        Err(SpatialInterferenceHeatingError::SourceRipple { .. })
    ));
}

#[test]
fn darkness_stays_exact_and_negative_heat_is_not_clipped() {
    let dark = [spectrum(0.0, 0.0, 0.0, 1.0)];
    let accepted = integrate_cells(
        &dark,
        Dd::ZERO,
        0.4,
        SpatialInterferenceHeatingControls {
            maximum_source_ripple_j: 0.0,
            ..controls()
        },
        &mut || false,
    )
    .unwrap();
    assert_eq!(accepted.total_j.to_f64(), 0.0);
    assert_eq!(accepted.cell_energy_j, [0.0]);
    assert_eq!(accepted.ripple_j, 0.0);
    assert!(
        integrate_cells(
            &[spectrum(-1.0, 0.0, 0.0, 1.0)],
            Dd::ZERO,
            0.4,
            controls(),
            &mut || false
        )
        .is_err()
    );
}

#[test]
fn refusal_is_bounded_and_measures_are_never_renormalized() {
    let sources = [spectrum(2.0, 1.0, 0.0, 3.0), spectrum(2.0, -1.0, 0.0, 3.0)];
    for dt in [0.0, -1.0, 2.0, f64::NAN] {
        assert!(integrate_cells(&sources, Dd::ZERO, dt, controls(), &mut || false).is_err());
    }
    assert!(
        integrate_cells(
            &sources,
            Dd::from_f64(f64::INFINITY),
            0.4,
            controls(),
            &mut || false
        )
        .is_err()
    );
    let mut calls = 0;
    let value = integrate_cells(&sources, Dd::ZERO, 0.4, controls(), &mut || {
        calls += 1;
        calls >= 7
    });
    assert!(value.is_err());
    assert_eq!(calls, 7);
    assert!(match_measure(Dd::from_f64(0.5), Dd::ONE).is_err());
    assert!(match_measure(Dd::ONE, Dd::ONE).is_ok());
}

#[test]
fn field_energy_retains_coherent_joules_in_native_cell_order() {
    use fs_fft::harmonic_window::{
        AngularHarmonicProductInput, HarmonicProductLaw, HarmonicProductSpectrumLimits,
        HarmonicWindow,
    };
    let frequencies = [
        562_949_953_421_312.0,
        562_949_953_421_313.0,
        562_949_953_421_314.0,
    ];
    let amplitudes = [1.0, -2.0, 1.0];
    let source = |gain: f64| {
        PreparedHarmonicProductSpectrum::from_angular_frequencies(
            (0..3).flat_map(|i| {
                (0..3).map(move |j| AngularHarmonicProductInput {
                    left_angular_frequency_rad_s: frequencies[i],
                    right_angular_frequency_rad_s: frequencies[j],
                    left_amplitude: [Dd::from_f64(gain * amplitudes[i]), Dd::ZERO],
                    right_amplitude: [Dd::from_f64(gain * amplitudes[j]), Dd::ZERO],
                    law: HarmonicProductLaw::ProperComplex,
                })
            }),
            HarmonicProductSpectrumLimits {
                maximum_products: 9,
                maximum_generated_terms: 9,
                maximum_workspace_bytes: 4096,
            },
        )
        .unwrap()
    };
    let spectra = [source(1.0), source(2.0)];
    let duration = 1e-12_f64;
    let start = Dd::from_f64(-0.5 * duration);
    let mut order = Vec::new();
    let actual = integrate_cells_with(
        &spectra,
        start,
        duration,
        controls(),
        &mut || false,
        |index, _, epoch, dt, cancelled| {
            poll(cancelled)?;
            order.push(index);
            assert_eq!(epoch, start);
            assert_eq!(dt, duration);
            // The volume owner supplies joules from the original coherent
            // fields. This native temporal factor isolates the cell-energy
            // seam from unrelated material and geometric preparation.
            let window = HarmonicWindow::new(epoch.to_f64(), dt).unwrap();
            let band = window
                .prepare_uniform_band_rad_s(&frequencies)
                .unwrap()
                .unwrap();
            let fields =
                amplitudes.map(|value| [Dd::from_f64(value * (index + 1) as f64), Dd::ZERO]);
            Ok(band
                .mean_square_parts_dd_scaled(&fields, Dd::from_f64(0.5))
                .unwrap()
                * Dd::from_f64(dt))
        },
    )
    .unwrap();
    let expected = duration.powi(5) / 160.0;
    assert_eq!(order, [0, 1]);
    assert!((actual.cell_energy_j[0] / expected - 1.0).abs() < 2e-10);
    assert!((actual.cell_energy_j[1] / (4.0 * expected) - 1.0).abs() < 2e-10);
    assert!((actual.total_j.to_f64() / (5.0 * expected) - 1.0).abs() < 2e-10);
    assert_eq!(actual.end_s, start + Dd::from_f64(duration));
    let cached = integrate_cells(&spectra, start, duration, controls(), &mut || false);
    assert!(cached.is_err() || (cached.unwrap().cell_energy_j[0] / expected - 1.0).abs() > 0.5);
}

#[test]
fn field_energy_preserves_dd_epoch_and_common_ripple_refusal() {
    let start = Dd::from_pair(9_007_199_254_740_992.0, 0.125);
    let duration = 0.25;
    let sources = [spectrum(2.0, 0.0, 0.0, 1.0)];
    let result = integrate_cells_with(
        &sources,
        start,
        duration,
        controls(),
        &mut || false,
        |_, _, epoch, dt, cancelled| {
            poll(cancelled)?;
            assert_eq!(epoch, start);
            assert_eq!(dt, duration);
            Ok(Dd::from_f64(0.3))
        },
    )
    .unwrap();
    assert_eq!(result.cell_energy_j, [0.3]);
    assert_eq!(result.end_s, start + Dd::from_f64(duration));
    assert_eq!(result.end_s.to_f64(), start.to_f64());

    let varying = [spectrum(2.0, 1.0, 0.0, 3.0)];
    let mut queries = 0;
    let refused = integrate_cells_with(
        &varying,
        Dd::ZERO,
        0.7,
        SpatialInterferenceHeatingControls {
            maximum_source_ripple_j: 0.0,
            ..controls()
        },
        &mut || false,
        |_, _, _, _, _| {
            queries += 1;
            Ok(Dd::ONE)
        },
    );
    assert!(matches!(
        refused,
        Err(SpatialInterferenceHeatingError::SourceRipple { .. })
    ));
    assert_eq!(queries, 0);
    assert!(
        integrate_cells_with(
            &sources,
            start,
            duration,
            controls(),
            &mut || false,
            |_, _, _, _, _| Ok(-Dd::ONE),
        )
        .is_err()
    );
    assert!(matches!(
        integrate_cells_with(
            &sources,
            start,
            duration,
            controls(),
            &mut || false,
            |_, _, _, _, _| Err(SpatialInterferenceHeatingError::Cancelled),
        ),
        Err(SpatialInterferenceHeatingError::Cancelled)
    ));
}
