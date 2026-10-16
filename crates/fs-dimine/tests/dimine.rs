//! Battery for dimensional mining (addendum Proposal 9's knowledge apex).
//! Covers power-law recovery (1D + multi-π), estimated-color, the fit-
//! significance gate (noise yields no significant law), extrapolation refusal
//! (the π-space envelope = convex hull in 1D), and the error paths
//! (too-few-samples, non-positive, singular, dim-mismatch), plus determinism.

use fs_dimine::{Color, MineError, Sample, fit_power_law};

fn one_d(coeff: f64, exp: f64, pis: &[f64]) -> Vec<Sample> {
    pis.iter()
        .map(|&p| Sample::new(vec![p], coeff * p.powf(exp)))
        .collect()
}

#[test]
fn recovers_a_one_dimensional_power_law() {
    // y = 2 * π^3 sampled exactly.
    let corpus = one_d(
        2.0,
        3.0,
        &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
    );
    let law = fit_power_law(&corpus).unwrap();
    assert!(
        (law.coefficient - 2.0).abs() < 1e-9,
        "C = {}",
        law.coefficient
    );
    assert!(
        (law.exponents[0] - 3.0).abs() < 1e-9,
        "a = {}",
        law.exponents[0]
    );
    assert!((law.r_squared - 1.0).abs() < 1e-9, "r2 = {}", law.r_squared);
    assert!(law.is_significant(0.99));
    // predict inside the trained range.
    assert!((law.predict(&[5.0]).unwrap() - 250.0).abs() < 1e-6);
    // it is estimated-color (a conjecture, never a certified bound).
    assert!(matches!(law.color, Color::Estimated { .. }));
}

#[test]
fn recovers_a_multi_pi_power_law() {
    // y = 1.5 * π1^2 * π2^{-1} over a 3x3 grid.
    let mut corpus = Vec::new();
    for &p1 in &[1.0, 2.0, 3.0] {
        for &p2 in &[1.0, 2.0, 3.0] {
            corpus.push(Sample::new(vec![p1, p2], 1.5 * p1.powi(2) * p2.powi(-1)));
        }
    }
    let law = fit_power_law(&corpus).unwrap();
    assert!((law.coefficient - 1.5).abs() < 1e-9);
    assert!((law.exponents[0] - 2.0).abs() < 1e-9);
    assert!((law.exponents[1] + 1.0).abs() < 1e-9);
    assert!(law.is_significant(0.99));
}

#[test]
fn noise_yields_no_significant_law() {
    // qoi unrelated to π (oscillating) -> a power-law fit is poor.
    let qois = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0, 5.0, 3.0, 5.0, 8.0];
    let corpus: Vec<Sample> = qois
        .iter()
        .enumerate()
        .map(|(i, &q)| Sample::new(vec![(i + 1) as f64], q))
        .collect();
    let law = fit_power_law(&corpus).unwrap();
    assert!(
        !law.is_significant(0.9),
        "noise must not fit: r2 = {}",
        law.r_squared
    );
}

#[test]
fn refuses_to_extrapolate_beyond_the_pi_envelope() {
    // trained on π in [1, 10]; the envelope IS the convex hull in 1D.
    let corpus = one_d(
        1.0,
        1.0,
        &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
    );
    let law = fit_power_law(&corpus).unwrap();
    assert_eq!(law.envelope, vec![(1.0, 10.0)]);
    // just inside is served; outside (either side) is refused.
    assert!(law.predict(&[5.0]).is_ok());
    assert!(law.predict(&[1.0]).is_ok()); // boundary inclusive
    assert!(matches!(
        law.predict(&[0.5]),
        Err(MineError::Extrapolation { .. })
    ));
    assert!(matches!(
        law.predict(&[10.5]),
        Err(MineError::Extrapolation { .. })
    ));
}

#[test]
fn too_few_samples_is_rejected() {
    // M=1 needs >= 3 samples.
    let corpus = one_d(2.0, 1.0, &[1.0, 2.0]);
    assert!(matches!(
        fit_power_law(&corpus),
        Err(MineError::TooFewSamples { have: 2, need: 3 })
    ));
}

#[test]
fn non_positive_values_are_rejected() {
    // a non-positive QoI (log undefined) fails closed.
    let mut corpus = one_d(2.0, 1.0, &[1.0, 2.0, 3.0, 4.0]);
    corpus[1].qoi = 0.0;
    assert!(matches!(
        fit_power_law(&corpus),
        Err(MineError::NonPositive { what: "qoi", .. })
    ));
    // a non-positive π too.
    let mut corpus2 = one_d(2.0, 1.0, &[1.0, 2.0, 3.0, 4.0]);
    corpus2[2].pi[0] = -1.0;
    assert!(matches!(
        fit_power_law(&corpus2),
        Err(MineError::NonPositive {
            what: "pi coordinate",
            ..
        })
    ));
}

#[test]
fn a_collinear_design_is_singular() {
    // all π identical -> the ln(π) column is constant -> rank-deficient.
    let corpus: Vec<Sample> = (0..4)
        .map(|i| Sample::new(vec![5.0], 2.0 + f64::from(i)))
        .collect();
    assert!(matches!(fit_power_law(&corpus), Err(MineError::Singular)));
}

#[test]
fn predict_rejects_a_dimension_mismatch() {
    let corpus = one_d(1.0, 1.0, &[1.0, 2.0, 3.0, 4.0]);
    let law = fit_power_law(&corpus).unwrap();
    assert!(matches!(
        law.predict(&[1.0, 2.0]),
        Err(MineError::DimMismatch {
            expected: 1,
            found: 2
        })
    ));
}

#[test]
fn fitting_is_deterministic() {
    let corpus = one_d(2.0, 3.0, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    let a = fit_power_law(&corpus).unwrap();
    let b = fit_power_law(&corpus).unwrap();
    assert_eq!(a.coefficient.to_bits(), b.coefficient.to_bits());
    assert_eq!(
        a.exponents.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        b.exponents.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
    );
}
