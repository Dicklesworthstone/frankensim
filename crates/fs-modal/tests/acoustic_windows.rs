//! Wide piano/soundboard bands must not turn small inverse-eigenvalue units
//! into a false Lanczos breakdown. Analytic pencils, unchanged work budgets.
use fs_modal::{ModalError, SliceOptions, shift_invert_modes_mfree, slice_window};
use fs_sparse::{Coo, Csr};
use std::f64::consts::TAU;

fn diagonal(values: &[f64]) -> Csr {
    let mut c = Coo::new(values.len(), values.len());
    for (i, value) in values.iter().enumerate() { c.push(i, i, *value); }
    c.assemble()
}

#[test]
fn matrix_free_modes_survive_inverse_operator_rescaling_without_extra_iterations() {
    for scale in [1e-6, 1.0, 1e10, 1e20] {
        let n = 12;
        let eigenvalues: Vec<f64> = (1..=n).map(|i| scale * i as f64).collect();
        let sigma = 6.5 * scale;
        let mut solves = 0;
        let got = shift_invert_modes_mfree(n, sigma, n, n,
            |x, y| y.copy_from_slice(x),
            |rhs, out| {
                solves += 1;
                for i in 0..n { out[i] = rhs[i] / (eigenvalues[i] - sigma); }
            }).expect("unit changes cannot discard modes");
        assert_eq!(got.len(), n);
        assert!(solves <= n, "the correction must not raise the work budget");
        for (i, (value, phi)) in got.iter().enumerate() {
            assert!((value / eigenvalues[i] - 1.0).abs() < 1e-9);
            let residual = phi.iter().enumerate().map(|(j, p)|
                ((eigenvalues[j] - value) / eigenvalues[i] * p).powi(2)).sum::<f64>().sqrt();
            assert!(residual < 1e-8);
            assert!((phi.iter().map(|p| p * p).sum::<f64>() - 1.0).abs() < 1e-10);
            for (_, other) in &got[..i] {
                assert!(phi.iter().zip(other).map(|(a, b)| a * b).sum::<f64>().abs() < 1e-9);
            }
        }
    }
}

#[test]
fn a_wide_audio_band_recovers_all_36_known_resonances_with_mass_certification() {
    let n = 36;
    // All modes lie at 50..1500 Hz, but a user requests a complete 20 kHz
    // band. At its midpoint shift, inverse eigenvalues cluster near -1e-10;
    // useful residual directions are smaller than the old absolute 1e-14.
    let values: Vec<f64> = (0..n).map(|i|
        (TAU * 50.0 * 30.0_f64.powf(i as f64 / (n - 1) as f64)).powi(2)).collect();
    let mass: Vec<f64> = (0..n).map(|i| [0.02, 0.2, 1.0, 5.0][i % 4]).collect();
    let m = diagonal(&mass);
    let k = diagonal(&values.iter().zip(&mass).map(|(l, m)| l * m).collect::<Vec<_>>());
    for equilibrated in [false, true] {
        let opts = SliceOptions { max_lanczos: n, mass_diagonal_equilibration: equilibrated,
            ..SliceOptions::default() };
        let report = slice_window(&k, &m, (0.0, (TAU * 20_000.0).powi(2)), &opts)
            .expect("complete wide acoustic slice");
        assert_eq!((report.below_low, report.below_high, report.expected), (0, n, n));
        assert_eq!(report.modes.len(), n);
        assert!(report.stats.lanczos_iters <= n);
        for (i, mode) in report.modes.iter().enumerate() {
            assert!((mode.lambda / values[i] - 1.0).abs() < 1e-7);
            assert!(mode.residual / values[i] < 1e-7);
            assert!(mode.interval.0 <= values[i] + 1e-7 * values[i]);
            assert!(mode.interval.1 >= values[i] - 1e-7 * values[i]);
            let norm: f64 = mode.phi.iter().zip(&mass).map(|(p, m)| p * p * m).sum();
            assert!((norm - 1.0).abs() < 1e-10);
            for other in &report.modes[..i] {
                let product: f64 = mode.phi.iter().zip(&other.phi).zip(&mass)
                    .map(|((a, b), m)| a * b * m).sum();
                assert!(product.abs() < 1e-9);
            }
        }
    }
}

#[test]
fn actual_degeneracy_still_restarts_at_large_physical_eigenvalues() {
    for scale in [1.0, 1e20] {
        let k = diagonal(&[scale; 2]); let m = diagonal(&[1.0; 2]);
        let report = slice_window(&k, &m, (0.5 * scale, 1.5 * scale),
            &SliceOptions::default()).unwrap();
        assert_eq!(report.expected, 2);
        assert!(report.stats.restarts >= 1);
        for mode in &report.modes { assert!((mode.lambda / scale - 1.0).abs() < 1e-10); }
        let product: f64 = report.modes[0].phi.iter().zip(&report.modes[1].phi)
            .map(|(a, b)| a * b).sum();
        assert!(product.abs() < 1e-10);
    }
}

#[test]
fn an_insufficient_budget_still_refuses_the_entire_wide_slice() {
    let k = diagonal(&(1..=12).map(|i| 1e20 * i as f64).collect::<Vec<_>>());
    let m = diagonal(&[1.0; 12]);
    let opts = SliceOptions { max_lanczos: 2, max_restarts: 1, ..SliceOptions::default() };
    assert!(matches!(slice_window(&k, &m, (0.5e20, 12.5e20), &opts),
        Err(ModalError::WindowUnresolved { expected: 12, converged, .. }) if converged < 12));
}
