//! G0/G1/G3: independent force-spectrum and known-frequency projection checks.

use super::*;
use fs_fft::harmonic::{HarmonicFitBudget, HarmonicFitError};
use std::f64::consts::TAU;

fn budget() -> HarmonicFitBudget {
    HarmonicFitBudget {
        max_workspace_bytes: 4 * 1024 * 1024,
        min_reciprocal_condition: 1e-10,
    }
}

fn close(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "actual={actual:e}, expected={expected:e}, tolerance={absolute:e}"
    );
}

#[test]
fn g1_prepared_force_spectrum_matches_independent_hann_dft() {
    let n = 128;
    let dt = 0.0001;
    let forces: Vec<_> = (0..n)
        .map(|i| {
            let theta = TAU * i as f64 / n as f64;
            [
                2.0 + 1.3 * (4.0 * theta).cos() - 0.7 * (9.0 * theta).sin(),
                -3.0 + 0.4 * (4.0 * theta).sin() + 0.03 * (31.0 * theta).cos(),
            ]
        })
        .collect();
    let spectrum = ForceSpectrum::new(&forces, dt, n).unwrap();
    for bin in [1, 4, 9, 31, 63] {
        let (actual, frequency) = spectrum.line(bin).unwrap();
        close(frequency, bin as f64 / (n as f64 * dt), 0.0);
        for axis in 0..2 {
            let mean = forces.iter().map(|f| f[axis]).sum::<f64>() / n as f64;
            let mut expected = C64::ZERO;
            for (i, force) in forces.iter().enumerate() {
                let theta = TAU * bin as f64 * i as f64 / n as f64;
                let hann = 0.5 - 0.5 * (TAU * i as f64 / (n - 1) as f64).cos();
                expected = expected
                    + C64::new(theta.cos(), theta.sin())
                        .scale((force[axis] - mean) * hann * 4.0 / (n - 1) as f64);
            }
            assert!((actual[axis] - expected).abs() < 3e-14);
        }
        assert_eq!(spectrum.line(bin).unwrap().0, actual);
    }
}

fn harmonic_record(n: usize, dt: f64, fundamental: f64) -> Vec<[f64; 2]> {
    (0..n)
        .map(|i| {
            let theta = TAU * fundamental * i as f64 * dt;
            [
                0.4 + 1.3 * theta.cos() - 0.2 * theta.sin()
                    + 0.00001 * (2.0 * theta).cos()
                    + 0.00002 * (2.0 * theta).sin()
                    - 0.15 * (3.0 * theta).cos(),
                -0.8 + 0.3 * theta.sin() - 0.07 * (2.0 * theta).cos() + 0.25 * (3.0 * theta).sin(),
            ]
        })
        .collect()
}

#[test]
fn g1_joint_force_fit_resolves_weak_off_bin_even_harmonics_and_radiation() {
    let n = 512;
    let dt = 1.0 / 8000.0;
    let f0 = 132.7;
    let plan = HarmonicFitPlan::new(n, dt, f0, &[1, 2, 3], budget()).unwrap();
    let record = harmonic_record(n, dt, f0);
    let (force, summary) = fit_force_harmonics(&record, dt, &plan).unwrap();
    let expected = [
        [C64::new(1.3, -0.2), C64::new(0.0, 0.3)],
        [C64::new(0.00001, 0.00002), C64::new(-0.07, 0.0)],
        [C64::new(-0.15, 0.0), C64::new(0.0, 0.25)],
    ];
    for (actual, expected) in force.iter().zip(expected) {
        for axis in 0..2 {
            assert!((actual[axis] - expected[axis]).abs() < 3e-12);
        }
    }
    close(summary[0].dc, 0.4, 3e-12);
    close(summary[1].dc, -0.8, 3e-12);
    for fit in &summary {
        assert!(fit.residual_rms < 3e-12);
    }
    let acoustic = FluidAcousticState::try_new(1.2, 340.0).unwrap();
    let fitted = JetForceHarmonics {
        acoustic,
        harmonics: plan.harmonics().to_vec(),
        fundamental_hz: f0,
        force_n_m: force,
        dc_force_n_m: [summary[0].dc, summary[1].dc],
        residual_rms_n_m: [summary[0].residual_rms, summary[1].residual_rms],
    };
    let observed = fitted.observe_harmonic(2, [0.4, 0.3], [0.0, 0.0]).unwrap();
    assert_eq!(observed.frequency_hz, 2.0 * f0);
    let reference =
        dipole_pressure(expected[1], TAU * 2.0 * f0 / 340.0, [0.4, 0.3], [0.0, 0.0]).unwrap();
    assert!((observed.field.pressure - reference.pressure).abs() < 1e-11);
    assert_eq!(observed.field.scope, SCOPE_STATEMENT);
    assert!(fitted.observe_harmonic(4, [0.4, 0.3], [0.0, 0.0]).is_err());
}

#[test]
fn g0_force_spectrum_budget_bins_and_original_fit_clock_refuse() {
    let forces = harmonic_record(64, 0.001, 31.7);
    assert!(ForceSpectrum::new(&forces, 0.001, 63).is_err());
    assert!(ForceSpectrum::new(&forces[..63], 0.001, 64).is_err());
    assert!(ForceSpectrum::new(&forces, 0.0, 64).is_err());
    let spectrum = ForceSpectrum::new(&forces, 0.001, 64).unwrap();
    for bin in [0, 32, usize::MAX] {
        assert!(spectrum.line(bin).is_err());
    }
    let plan = HarmonicFitPlan::new(64, 0.001, 31.7, &[1, 2, 3], budget()).unwrap();
    assert!(fit_force_harmonics(&forces[..63], 0.001, &plan).is_err());
    assert!(fit_force_harmonics(&forces, 0.002, &plan).is_err());
    let mut invalid = forces;
    invalid[5][1] = f64::NAN;
    assert!(matches!(
        fit_force_harmonics(&invalid, 0.001, &plan),
        Err(AeroacError::HarmonicFit(HarmonicFitError::NonFinite { .. }))
    ));
}

#[test]
fn g3_joint_force_fit_preserves_phase_and_physical_force_scaling() {
    let n = 256;
    let dt = 1.0 / 8000.0;
    let plan = HarmonicFitPlan::new(n, dt, 132.7, &[1, 2, 3], budget()).unwrap();
    let record = harmonic_record(n, dt, 132.7);
    let scale = 0.037;
    let scaled: Vec<_> = record.iter().map(|f| f.map(|x| x * scale)).collect();
    let (original, baseline) = fit_force_harmonics(&record, dt, &plan).unwrap();
    let (changed, summary) = fit_force_harmonics(&scaled, dt, &plan).unwrap();
    for (actual, expected) in changed.iter().zip(original) {
        for axis in 0..2 {
            assert!((actual[axis] - expected[axis].scale(scale)).abs() < 3e-14);
        }
    }
    for axis in 0..2 {
        close(summary[axis].dc, baseline[axis].dc * scale, 3e-14);
    }
}
