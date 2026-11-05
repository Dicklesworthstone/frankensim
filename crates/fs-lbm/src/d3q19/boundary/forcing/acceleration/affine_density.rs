//! Scalar closure for a regularized boundary whose retained mass depends on
//! its own density-proportional source. Population rows have the form
//! `f_q(rho) = base_q + response_q * rho`; unchanged wall rows have zero response.

/// Solve the actual zeroth moment of a density-dependent source reconstruction.
/// This direct rank-one closure also works when fixed-point iteration would
/// diverge. A singular or nonpositive density response is not an admitted root.
pub(super) fn solve<const Q: usize>(base: &[f64; Q], response: &[f64; Q]) -> Option<f64> {
    if !base.iter().chain(response).all(|value| value.is_finite()) {
        return None;
    }
    let mass = base.iter().sum::<f64>();
    let denominator = 1.0 - response.iter().sum::<f64>();
    let density = mass / denominator;
    (mass.is_finite()
        && mass > 0.0
        && denominator.is_finite()
        && denominator > 0.0
        && density.is_finite()
        && density > 0.0)
        .then_some(density)
}

/// Recheck the result after the native source transform, including its own
/// rest-population roundoff correction. Scale with population mass, without a
/// unit-sized absolute tolerance that could admit a wrong small-density root.
pub(super) fn matches<const Q: usize>(populations: &[f64; Q], expected: f64) -> bool {
    let density = populations.iter().sum::<f64>();
    let scale = populations
        .iter()
        .map(|value| value.abs())
        .sum::<f64>()
        .max(expected.abs());
    expected.is_finite()
        && expected > 0.0
        && density.is_finite()
        && density > 0.0
        && scale.is_finite()
        && (density - expected).abs() <= 64.0 * f64::EPSILON * scale
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reservoir_and_retained_wall_rows_close_actual_mass() {
        // Rows 0 and 3 are retained wall data. Rows 1 and 2 are reservoir-owned;
        // their combined unit-density response is -1/8, hence rho = 8/9.
        let base = [0.5, 0.2, 0.1, 0.2];
        let response = [0.0, -0.2, 0.075, 0.0];
        let rho = solve(&base, &response).unwrap();
        assert!((rho - 8.0 / 9.0).abs() < 2e-16);
        let populations = core::array::from_fn::<_, 4, _>(|q| base[q] + rho * response[q]);
        assert!(matches(&populations, rho));
        assert_eq!(populations[0].to_bits(), base[0].to_bits());
        assert_eq!(populations[3].to_bits(), base[3].to_bits());
        // Reconstructing with the nominal template mass (1) misses the closure.
        let lagged = core::array::from_fn::<_, 4, _>(|q| base[q] + response[q]);
        assert!(!matches(&lagged, rho));
    }

    #[test]
    fn direct_root_agrees_with_independent_converged_substitution() {
        let base = [0.35, 0.2, 0.2, 0.15, 0.1];
        let response = [0.0, 0.14, -0.03, 0.08, 0.0];
        let mut iterated = 1.0;
        for _ in 0..100 {
            iterated = base
                .iter()
                .zip(response)
                .map(|(base, response)| base + response * iterated)
                .sum();
        }
        let exact = solve(&base, &response).unwrap();
        assert!((iterated - exact).abs() < 4e-16);
        let populations = core::array::from_fn::<_, 5, _>(|q| base[q] + exact * response[q]);
        assert!(matches(&populations, exact));
    }

    #[test]
    fn direct_root_does_not_require_a_contractive_fixed_point() {
        let base = [0.6, 0.3, 0.1];
        let response = [-0.5, -1.0, -0.5];
        let rho = solve(&base, &response).unwrap();
        assert!((rho - 1.0 / 3.0).abs() < 1e-16);
        let populations = core::array::from_fn::<_, 3, _>(|q| base[q] + rho * response[q]);
        assert!(matches(&populations, rho));
    }

    #[test]
    fn closure_and_admission_scale_with_density_units() {
        for scale in [1e-200, 1.0, 1e200] {
            let base = [0.4 * scale, 0.6 * scale];
            let response = [0.125, 0.125];
            let rho = solve(&base, &response).unwrap();
            assert!((rho / scale - 4.0 / 3.0).abs() < 4e-16);
            assert!(matches(
                &[base[0] + rho * response[0], base[1] + rho * response[1]],
                rho
            ));
            assert!(!matches(&[rho * 0.5, rho * (0.5 + 1e-10)], rho));
        }
    }

    #[test]
    fn singular_negative_and_unrepresentable_roots_refuse() {
        for response in [
            [0.5, 0.5],
            [1.0, 0.5],
            [f64::NAN, 0.0],
            [f64::INFINITY, 0.0],
        ] {
            assert!(solve(&[0.5, 0.5], &response).is_none());
        }
        for base in [
            [0.0, 0.0],
            [-1.0, 0.0],
            [f64::MAX, f64::MAX],
            [f64::NAN, 1.0],
        ] {
            assert!(solve(&base, &[0.0, 0.0]).is_none());
        }
        assert!(solve(&[f64::MAX, 0.0], &[0.5, 0.0]).is_none());
        assert!(!matches(&[f64::INFINITY], 1.0));
        assert!(!matches(&[1.0], f64::NAN));
        assert!(!matches(&[0.0], 0.0));
    }
}
