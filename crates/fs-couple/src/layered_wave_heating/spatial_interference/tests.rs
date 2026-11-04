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
