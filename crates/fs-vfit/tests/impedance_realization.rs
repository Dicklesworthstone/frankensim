//! Admission regressions for computed impedance to passive runtime filters.

use fs_math::c64::C64;
use fs_vfit::FitOptions;
use fs_vfit::discretize::{
    DiscretizeError, RealizeError, realize_tabulated_impedance,
    realize_tabulated_impedance_with_report,
};
use fs_vfit::passivity::{GlobalPassivityBudget, GlobalPassivityVerdict};

/// Positive-real operator 0.4 + 0.03*s + (4*s+16)/(s²+6*s+25).
/// The proper term's real numerator is 400 + 8*omega², strictly positive.
fn impedance(omega: f64) -> C64 {
    C64::new(0.4, 0.03 * omega)
        + C64::new(16.0, 4.0 * omega) * C64::new(25.0 - omega * omega, 6.0 * omega).recip()
}

fn fixture() -> (Vec<f64>, Vec<C64>, FitOptions) {
    let mut omega: Vec<_> = (0..96)
        .map(|index| 0.5 + 39.5 * index as f64 / 95.0)
        .collect();
    // The passivity/accuracy band belongs to all samples, not their endpoints.
    omega.rotate_left(37);
    let data = omega
        .iter()
        .map(|&frequency| impedance(frequency))
        .collect();
    (omega, data, FitOptions::new(2))
}

#[test]
fn reported_impedance_retains_physical_error_and_certifies_both_operators() {
    let (omega, data, options) = fixture();
    let (dt, prewarp) = (0.02, 8.0);
    let result = realize_tabulated_impedance_with_report(
        &omega,
        &data,
        dt,
        &options,
        prewarp,
        GlobalPassivityBudget::default(),
    )
    .expect("strict passive realization");
    assert_eq!(
        result.repair.certificate.verdict,
        GlobalPassivityVerdict::Certified
    );
    assert_eq!(
        result.digital_certificate.verdict,
        GlobalPassivityVerdict::Certified
    );
    assert!(result.filter.is_stable());
    assert!(
        result.sample_errors.continuous_max_abs < 2.0e-6,
        "{:?}",
        result.sample_errors
    );
    assert!(
        result.sample_errors.digital_max_abs > 1.0e-3,
        "frequency warping must not be hidden behind the continuous fit error"
    );
    assert_eq!(result.repair.rounds, 0, "a passive fit must not be changed");
    assert_eq!(result.repair.relative_perturbation, 0.0);
    // Independent second-order polynomial evaluates the known physical law;
    // it does not use the fitted poles/residues or the bilinear implementation.
    let map_scale = prewarp / (0.5 * prewarp * dt).tan();
    for frequency in [0.0, 0.7, 3.4, 8.0, 17.3, 39.0, 80.0] {
        let mapped = map_scale * (0.5 * frequency * dt).tan();
        let expected = impedance(mapped);
        let actual = result.filter.eval(frequency);
        assert!(
            (actual - expected).abs() < 4.0e-6 * expected.abs().max(1.0),
            "frequency={frequency}, actual={actual:?}, expected={expected:?}"
        );
    }
    let simple = realize_tabulated_impedance(&omega, &data, dt, &options, prewarp)
        .expect("compatibility entry uses the same admitted preparation");
    assert_eq!(simple, result.filter);
}

#[test]
fn exhausted_global_admission_is_not_silently_replaced_with_the_raw_fit() {
    let (omega, data, options) = fixture();
    let budget = GlobalPassivityBudget {
        max_order: 0,
        ..GlobalPassivityBudget::default()
    };
    assert!(matches!(
        realize_tabulated_impedance_with_report(&omega, &data, 0.02, &options, 8.0, budget,),
        Err(RealizeError::Passivity(_))
    ));
}

#[test]
fn sampled_impedance_refuses_an_unrepresentable_sample_clock_or_aliased_band() {
    let (omega, data, options) = fixture();
    assert!(matches!(
        realize_tabulated_impedance(&omega, &data, 0.0, &options, 8.0),
        Err(RealizeError::Discretize(DiscretizeError::BadSampleInterval))
    ));
    // Nyquist is about 31.4 rad/s; the requested band reaches 40 rad/s.
    assert!(matches!(
        realize_tabulated_impedance(&omega, &data, 0.1, &options, 8.0),
        Err(RealizeError::Discretize(
            DiscretizeError::BeyondNyquist { .. }
        ))
    ));
}
