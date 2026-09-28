//! Noisy expected improvement for expensive, stochastic objectives.
//!
//! Unlike EI against the smallest noisy observation, q-NEI integrates the
//! uncertain incumbent: E[(min f(baseline) - min f(candidates))_+]. Both sets
//! are drawn from ONE latent GP posterior. Observation noise belongs in the
//! GP fit, not in these latent draws. All values remain model-based Monte
//! Carlo estimates, not confidence bounds or certificates.

use crate::gp::Gp;

/// Fixed common-random-number bank with arbitrary positive width.
///
/// The leading columns use the existing scrambled-Sobol normal bank. Beyond
/// the embedded Sobol dimension ceiling, independently keyed Philox uniforms
/// pass through the same deterministic normal quantile. This tail is MC, NOT
/// additional Sobol dimensions. Its stream identity depends on sample index,
/// not on bank width, so appending columns preserves previous tail draws.
///
/// # Panics
/// Rejects empty shapes, sample counts outside the Sobol index space, and
/// row-major length overflow before allocating.
#[must_use]
pub fn joint_normal_bank(samples: usize, width: usize, seed: u64) -> Vec<f64> {
    assert!(samples > 0, "joint normal bank needs samples");
    assert!(width > 0, "joint normal bank needs columns");
    assert!(u32::try_from(samples).is_ok(), "too many joint normal samples");
    let len = samples.checked_mul(width).expect("joint normal bank length overflow");
    let leading = width.min(fs_rand::qmc::MAX_SOBOL_DIM);
    let sobol = crate::acq::normal_bank(samples, leading, seed);
    if leading == width {
        return sobol;
    }
    let mut bank = vec![0.0; len];
    for (s, row) in bank.chunks_exact_mut(width).enumerate() {
        row[..leading].copy_from_slice(&sobol[s * leading..(s + 1) * leading]);
        let mut stream = fs_rand::StreamKey {
            seed,
            kernel: 0x514E_4549,
            tile: u32::try_from(s).expect("validated sample count"),
        }
        .stream();
        for value in &mut row[leading..] {
            *value = crate::acq::phi_inv(stream.next_f64().clamp(1e-12, 1.0 - 1e-12));
        }
    }
    bank
}

/// q-noisy expected improvement for MINIMIZATION over a fixed normal bank.
///
/// `baseline` is the caller's incumbent set, normally all evaluated inputs.
/// `bank` is row-major with `baseline.len() + candidates.len()` columns:
/// baseline columns first, then candidates. Reuse the same bank while
/// optimizing an acquisition; take row-wise prefixes when growing a batch.
///
/// Equal coordinates are the SAME latent random variable, including signed
/// zero. Only their first occurrence consumes a bank column. In particular,
/// a candidate already in the baseline has exactly zero improvement, even
/// when its observation is noisy. This is improvement of the latent design,
/// not the information value of taking a replicate measurement.
///
/// Uses the existing `Gp::predict_joint` numerical model, including its
/// adaptive jitter and documented degenerate-covariance fallback. Neither
/// exact integration, cross-ISA replay nor global acquisition optimality is
/// asserted by this function.
///
/// # Panics
/// Rejects empty sets, wrong dimensions, non-finite coordinates, empty/ragged
/// or non-finite banks, and non-finite posterior/sample arithmetic.
#[must_use]
pub fn q_noisy_expected_improvement(
    gp: &Gp,
    baseline: &[Vec<f64>],
    candidates: &[Vec<f64>],
    bank: &[f64],
) -> f64 {
    assert!(!baseline.is_empty(), "q-NEI needs a nonempty baseline");
    assert!(!candidates.is_empty(), "q-NEI needs candidates");
    let width = baseline.len().checked_add(candidates.len()).expect("q-NEI width overflow");
    let dim = gp.kernel.lengthscales.len();
    assert!(dim > 0, "q-NEI needs a positive input dimension");
    for point in baseline.iter().chain(candidates) {
        assert_eq!(point.len(), dim, "q-NEI input dimension mismatch");
        assert!(point.iter().all(|v| v.is_finite()), "q-NEI coordinates must be finite");
    }
    assert!(!bank.is_empty(), "q-NEI bank must not be empty");
    assert_eq!(bank.len() % width, 0, "q-NEI bank must be rectangular");
    assert!(bank.iter().all(|v| v.is_finite()), "q-NEI bank must be finite");

    // Stable first-occurrence ordering preserves the baseline prefix and
    // does not turn Cholesky jitter at duplicate locations into fake gain.
    let mut unique: Vec<Vec<f64>> = Vec::new();
    let mut columns = Vec::new();
    let mut baseline_len = 0;
    for (column, point) in baseline.iter().chain(candidates).enumerate() {
        if !unique.iter().any(|previous| previous == point) {
            unique.push(point.clone());
            columns.push(column);
        }
        if column + 1 == baseline.len() {
            baseline_len = unique.len();
        }
    }
    if baseline_len == unique.len() {
        return 0.0;
    }
    let (mean, lower) = gp.predict_joint(&unique);
    assert!(mean.iter().chain(&lower).all(|v| v.is_finite()), "q-NEI posterior must be finite");
    let n = unique.len();
    let mut average = 0.0;
    for (sample, row) in bank.chunks_exact(width).enumerate() {
        let mut incumbent = f64::INFINITY;
        let mut candidate = f64::INFINITY;
        for i in 0..n {
            let mut value = mean[i];
            for j in 0..=i {
                value = lower[i * n + j].mul_add(row[columns[j]], value);
            }
            assert!(value.is_finite(), "q-NEI posterior draw overflow");
            if i < baseline_len {
                incumbent = incumbent.min(value);
            } else {
                candidate = candidate.min(value);
            }
        }
        let improvement = (incumbent - candidate).max(0.0);
        assert!(improvement.is_finite(), "q-NEI improvement overflow");
        // Nonnegative running mean avoids overflowing a representable average
        // merely because the unnormalized sum exceeds f64's range.
        average += (improvement - average) / (sample + 1) as f64;
    }
    average
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::{Kernel, Matern};
    use std::panic::{AssertUnwindSafe, catch_unwind};

    fn gp(y: f64, noise: f64) -> Gp {
        Gp::fit(
            &[vec![0.0]],
            &[y],
            Kernel { family: Matern::FiveHalves, signal: 1.0, lengthscales: vec![1.0] },
            noise,
        )
    }

    #[test]
    fn correlated_two_point_gaussian_matches_closed_form() {
        let gp = gp(0.0, 0.3);
        let points = vec![vec![0.0], vec![0.7]];
        let (_, l) = gp.predict_joint(&points);
        // f0-f1 is a centered Gaussian with variance (L00-L10)^2+L11^2.
        let sd = fs_math::det::sqrt((l[0] - l[2]).powi(2) + l[3].powi(2));
        let expected = sd / fs_math::det::sqrt(2.0 * core::f64::consts::PI);
        let bank = joint_normal_bank(32768, 2, 71);
        let value = q_noisy_expected_improvement(&gp, &points[..1], &points[1..], &bank);
        assert!(value > 0.05, "non-vacuous acquisition");
        assert!((value - expected).abs() < 0.005 * expected + 0.001, "{value} vs {expected}");
    }

    #[test]
    fn noisy_outlier_is_not_a_fixed_incumbent() {
        let gp = gp(-10.0, 100.0);
        let bank = joint_normal_bank(8192, 2, 19);
        let noisy = q_noisy_expected_improvement(&gp, &[vec![0.0]], &[vec![3.0]], &bank);
        let plug_in = crate::acq::expected_improvement(&gp, &[3.0], -10.0, 0.0);
        assert!(noisy > 0.1, "uncertain incumbent must not suppress exploration");
        assert!(plug_in < 1e-12, "fixture must expose the noisy-incumbent failure");
    }

    #[test]
    fn near_coincident_points_keep_their_cross_covariance() {
        let gp = gp(0.0, 0.3);
        let bank = joint_normal_bank(8192, 2, 27);
        let value = q_noisy_expected_improvement(&gp, &[vec![0.0]], &[vec![0.001]], &bank);
        let independent_sd = fs_math::det::sqrt(gp.predict(&[0.0]).1 + gp.predict(&[0.001]).1);
        assert!(independent_sd > 0.5);
        assert!(value < 0.01, "independent incumbent/candidate draws would fail: {value}");
    }

    #[test]
    fn baseline_duplicates_have_exactly_zero_improvement() {
        let gp = gp(-3.0, 2.0);
        let bank = joint_normal_bank(128, 3, 4);
        assert_eq!(q_noisy_expected_improvement(&gp, &[vec![0.0]], &[vec![-0.0], vec![0.0]], &bank), 0.0);
    }

    #[test]
    fn repeated_candidates_reuse_the_first_latent_draw() {
        let gp = gp(0.0, 0.2);
        let bank = joint_normal_bank(256, 2, 6);
        let mut repeated = Vec::new();
        for row in bank.chunks_exact(2) {
            repeated.extend_from_slice(&[row[0], row[1], 1e6]);
        }
        let one = q_noisy_expected_improvement(&gp, &[vec![0.0]], &[vec![1.0]], &bank);
        let two = q_noisy_expected_improvement(&gp, &[vec![0.0]], &[vec![1.0], vec![1.0]], &repeated);
        assert_eq!(one.to_bits(), two.to_bits());
        let mut repeated_baseline = Vec::new();
        for row in bank.chunks_exact(2) {
            repeated_baseline.extend_from_slice(&[row[0], -1e6, row[1]]);
        }
        let duplicate = q_noisy_expected_improvement(&gp, &[vec![0.0], vec![0.0]], &[vec![1.0]], &repeated_baseline);
        assert_eq!(one.to_bits(), duplicate.to_bits());
    }

    #[test]
    fn noiseless_incumbent_recovers_expected_improvement() {
        let gp = gp(0.0, 0.0);
        let bank = joint_normal_bank(32768, 2, 73);
        let noisy = q_noisy_expected_improvement(&gp, &[vec![0.0]], &[vec![1.0]], &bank);
        let ei = crate::acq::expected_improvement(&gp, &[1.0], 0.0, 0.0);
        assert!((noisy - ei).abs() < 0.005 * ei + 0.001);
    }

    #[test]
    fn wide_banks_are_replayable_with_stable_prefixes() {
        let width = fs_rand::qmc::MAX_SOBOL_DIM + 2;
        let small = joint_normal_bank(16, width, 83);
        let wide = joint_normal_bank(16, width + 3, 83);
        assert_eq!(small, joint_normal_bank(16, width, 83));
        assert!(wide.iter().all(|v| v.is_finite()));
        for (a, b) in small.chunks_exact(width).zip(wide.chunks_exact(width + 3)) {
            assert_eq!(a, &b[..width]);
        }
        assert_ne!(small, joint_normal_bank(16, width, 84));
    }

    #[test]
    fn invalid_banks_and_points_are_rejected() {
        let gp = gp(0.0, 0.2);
        let baseline = vec![vec![0.0]];
        let candidate = vec![vec![1.0]];
        for bank in [vec![], vec![0.0], vec![0.0, f64::NAN], vec![f64::INFINITY, 0.0]] {
            assert!(catch_unwind(AssertUnwindSafe(|| q_noisy_expected_improvement(&gp, &baseline, &candidate, &bank))).is_err());
        }
        for point in [vec![], vec![1.0, 2.0], vec![f64::NAN]] {
            assert!(catch_unwind(AssertUnwindSafe(|| q_noisy_expected_improvement(&gp, &baseline, &[point], &[0.0, 0.0]))).is_err());
        }
        assert!(catch_unwind(|| joint_normal_bank(2, usize::MAX, 0)).is_err());
        assert!(catch_unwind(|| joint_normal_bank(0, 1, 0)).is_err());
        assert!(catch_unwind(|| joint_normal_bank(1, 0, 0)).is_err());
    }
}
