//! A signed observation is not an impedance: recover its actual direct terms.
use fs_math::c64::C64;
use fs_vfit::{FitOptions, PoleTerm, RationalModel, vector_fit, vector_fit_transfer};

fn grid() -> Vec<f64> {
    (0..120)
        .map(|i| 0.1 * (1000.0_f64).powf(i as f64 / 119.0))
        .collect()
}
fn truth(d: f64, e: f64) -> RationalModel {
    RationalModel {
        terms: vec![PoleTerm::Real {
            pole: -6.0,
            residue: 2.6,
        }],
        d,
        e,
    }
}
fn samples(model: &RationalModel, w: &[f64]) -> Vec<C64> {
    w.iter().map(|&w| model.eval_iw(w)).collect()
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 3e-10 * (1.0 + b.abs()),
        "{a:.17e} versus {b:.17e}"
    );
}

#[test]
fn signed_pressure_transfer_and_mode_reversal_recover_direct_term_and_off_grid_response() {
    let w = grid();
    let reference = truth(-0.65, 0.0);
    let h = samples(&reference, &w);
    let opts = FitOptions {
        fit_e: false,
        ..FitOptions::new(1)
    };
    let fitted = vector_fit_transfer(&w, &h, &opts).unwrap();
    assert!(fitted.model.d < 0.0);
    assert_eq!(fitted.model.e, 0.0);
    close(fitted.model.d, reference.d);
    assert!(fitted.model.is_stable());
    assert!(fitted.report.weighted_rms < 1e-11);
    // Signed modal normalization reverses the complete transfer, including D.
    let reversed: Vec<_> = h.iter().map(|v| v.scale(-7.5)).collect();
    let flipped = vector_fit_transfer(&w, &reversed, &opts).unwrap();
    close(flipped.model.d, -7.5 * reference.d);
    for w in [0.13, 0.7, 4.2, 21.3, 63.0] {
        let expected = reference.eval_iw(w);
        assert!((fitted.model.eval_iw(w) - expected).abs() < 2e-10 * expected.abs());
        assert!((flipped.model.eval_iw(w) - expected.scale(-7.5)).abs() < 2e-10 * expected.abs());
    }
    let repeated = vector_fit_transfer(&w, &h, &opts).unwrap();
    assert_eq!(repeated.model, fitted.model);
}

#[test]
fn signed_improper_term_is_recovered_and_disabled_direct_terms_remain_zero() {
    let w = grid();
    let reference = truth(-0.4, -0.015);
    let fitted = vector_fit_transfer(&w, &samples(&reference, &w), &FitOptions::new(1)).unwrap();
    close(fitted.model.d, reference.d);
    close(fitted.model.e, reference.e);
    assert!(fitted.model.d < 0.0 && fitted.model.e < 0.0);
    for w in [0.31, 3.7, 17.9, 79.1] {
        let expected = reference.eval_iw(w);
        assert!((fitted.model.eval_iw(w) - expected).abs() < 2e-10 * expected.abs());
    }
    let proper = truth(0.0, 0.0);
    let pinned = vector_fit_transfer(
        &w,
        &samples(&proper, &w),
        &FitOptions {
            fit_d: false,
            fit_e: false,
            ..FitOptions::new(1)
        },
    )
    .unwrap();
    assert_eq!((pinned.model.d, pinned.model.e), (0.0, 0.0));
    assert!(pinned.report.weighted_rms < 1e-11);
}

#[test]
fn original_impedance_policy_and_positive_data_arithmetic_are_unchanged() {
    let w = grid();
    let signed = truth(-0.65, -0.015);
    let h = samples(&signed, &w);
    let opts = FitOptions::new(1);
    let constrained = vector_fit(&w, &h, &opts).unwrap();
    assert!(constrained.model.d >= 0.0 && constrained.model.e >= 0.0);
    assert!(constrained.report.weighted_rms > 1e-3);
    let fixed = fs_vfit::vf::residue_fit_at_poles(&w, &h, &signed.terms, &opts).unwrap();
    assert!(fixed.model.d >= 0.0 && fixed.model.e >= 0.0);
    let positive = truth(0.65, 0.015);
    let h = samples(&positive, &w);
    let impedance = vector_fit(&w, &h, &opts).unwrap();
    let transfer = vector_fit_transfer(&w, &h, &opts).unwrap();
    assert_eq!(impedance.model, transfer.model);
    assert_eq!(impedance.report.weighted_rms, transfer.report.weighted_rms);
    assert_eq!(
        impedance.report.iterations_run,
        transfer.report.iterations_run
    );
}
