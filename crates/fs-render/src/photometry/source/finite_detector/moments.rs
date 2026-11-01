//! Fixed-size moments for sparsely hit sensor pixels.
//!
//! Missing sample indices are zero receiving power, not missing observations.
//! Consecutive zeros are merged analytically in constant work per component,
//! so a sparse detector does not require touching every pixel on every launch.
//! Each component is normalized by its largest observed magnitude before the
//! variance update; dimensionful powers are never squared. These are empirical
//! sample errors, with independence and complete-spectrum grouping supplied by
//! the transport owner, not certified confidence or optical error bounds.

use super::{FiniteDetectorError, MAX_FINITE_DETECTOR_SAMPLES};
use fs_math::det;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Moments<const N: usize = 4> {
    scale: [f64; N],
    mean: [f64; N],
    m2: [f64; N],
    count: u32,
}

impl<const N: usize> Default for Moments<N> {
    fn default() -> Self {
        Self {
            scale: [0.0; N],
            mean: [0.0; N],
            m2: [0.0; N],
            count: 0,
        }
    }
}

impl<const N: usize> Moments<N> {
    /// Add one complete source sample at its original zero-based launch index.
    /// Earlier unrecorded indices contribute zero. A repeated/out-of-order
    /// index or refused value leaves every component and count unchanged.
    pub(super) fn add_at(
        &mut self,
        values: [f64; N],
        sample_index: u32,
    ) -> Result<(), FiniteDetectorError> {
        if sample_index < self.count || sample_index >= MAX_FINITE_DETECTOR_SAMPLES {
            return Err(FiniteDetectorError::InvalidRequest);
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(FiniteDetectorError::NumericalRange);
        }
        let mut next = *self;
        next.pad_to(sample_index)?;
        next.count += 1;
        for (axis, value) in values.into_iter().enumerate() {
            let magnitude = value.abs();
            if magnitude > next.scale[axis] {
                let ratio = next.scale[axis] / magnitude;
                next.mean[axis] *= ratio;
                next.m2[axis] *= ratio * ratio;
                next.scale[axis] = magnitude;
            }
            if next.scale[axis] == 0.0 {
                continue;
            }
            let scaled = value / next.scale[axis];
            let delta = scaled - next.mean[axis];
            next.mean[axis] += delta / f64::from(next.count);
            next.m2[axis] += delta * (scaled - next.mean[axis]);
        }
        next.check_range()?;
        *self = next;
        Ok(())
    }

    /// Account for all remaining launches through `total_samples`, including
    /// pixels never hit at all. Padding to the current count is a no-op.
    pub(super) fn pad_to(&mut self, total_samples: u32) -> Result<(), FiniteDetectorError> {
        if total_samples < self.count || total_samples > MAX_FINITE_DETECTOR_SAMPLES {
            return Err(FiniteDetectorError::InvalidRequest);
        }
        if total_samples == self.count {
            return Ok(());
        }
        let mut next = *self;
        if next.count != 0 {
            let retained = f64::from(next.count) / f64::from(total_samples);
            let zeros = f64::from(total_samples - next.count);
            // Parallel-variance merge with a block whose mean and M2 are zero:
            // M2' = M2 + mean² * n*z/(n+z), mean' = mean*n/(n+z).
            let cross_weight = f64::from(next.count) * (zeros / f64::from(total_samples));
            for axis in 0..N {
                let mean = next.mean[axis];
                next.m2[axis] += mean * mean * cross_weight;
                next.mean[axis] *= retained;
            }
        }
        next.count = total_samples;
        next.check_range()?;
        *self = next;
        Ok(())
    }

    pub(super) const fn count(&self) -> u32 {
        self.count
    }

    /// Mean component values. The receiving owner admits the resulting Stokes
    /// or color state; signed components are intentionally permitted here.
    pub(super) fn mean(&self) -> [f64; N] {
        core::array::from_fn(|axis| self.scale[axis] * self.mean[axis])
    }

    pub(super) fn standard_error(&self) -> Result<[f64; N], FiniteDetectorError> {
        if self.count < 2 {
            return Err(FiniteDetectorError::InvalidRequest);
        }
        let denominator = f64::from(self.count) * f64::from(self.count - 1);
        let errors =
            core::array::from_fn(|axis| self.scale[axis] * det::sqrt(self.m2[axis] / denominator));
        if errors.iter().all(|value| value.is_finite()) {
            Ok(errors)
        } else {
            Err(FiniteDetectorError::NumericalRange)
        }
    }

    fn check_range(&self) -> Result<(), FiniteDetectorError> {
        if (0..N).any(|axis| {
            !(self.scale[axis] * self.mean[axis]).is_finite()
                || !self.m2[axis].is_finite()
                || self.m2[axis] < 0.0
        }) {
            Err(FiniteDetectorError::NumericalRange)
        } else {
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
        );
    }

    #[test]
    fn g1_sparse_gaps_match_explicit_zero_samples_and_two_pass_statistics() {
        let samples = [
            [2.0, -4.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [8.0, 3.0],
            [0.0, 0.0],
            [-4.0, 1.0],
            [0.0, 0.0],
        ];
        let mut sparse = Moments::<2>::default();
        let mut dense = Moments::<2>::default();
        for (index, value) in samples.into_iter().enumerate() {
            dense.add_at(value, index as u32).unwrap();
            if value != [0.0; 2] {
                sparse.add_at(value, index as u32).unwrap();
            }
        }
        sparse.pad_to(samples.len() as u32).unwrap();
        assert_eq!(sparse.count(), 7);
        for axis in 0..2 {
            let mean = samples.iter().map(|value| value[axis]).sum::<f64>() / 7.0;
            let error = (samples
                .iter()
                .map(|value| (value[axis] - mean).powi(2))
                .sum::<f64>()
                / (7.0 * 6.0))
                .sqrt();
            close(sparse.mean()[axis], mean, 2e-15);
            close(sparse.standard_error().unwrap()[axis], error, 2e-15);
            close(sparse.mean()[axis], dense.mean()[axis], 2e-15);
            close(
                sparse.standard_error().unwrap()[axis],
                dense.standard_error().unwrap()[axis],
                2e-15,
            );
        }
    }

    #[test]
    fn g1_single_rare_hit_keeps_all_launches_without_squaring_large_power() {
        let mut sparse = Moments::<2>::default();
        sparse.add_at([1e200, -2e200], 97_531).unwrap();
        sparse.pad_to(MAX_FINITE_DETECTOR_SAMPLES).unwrap();
        for (axis, power) in [1e200_f64, -2e200].into_iter().enumerate() {
            // One nonzero observation among n launches has mean a/n and
            // standard error |a|/n, even though a² cannot fit in binary64.
            let mean = power / f64::from(MAX_FINITE_DETECTOR_SAMPLES);
            close(sparse.mean()[axis], mean, mean.abs() * 4e-15);
            close(
                sparse.standard_error().unwrap()[axis],
                mean.abs(),
                mean.abs() * 4e-15,
            );
        }
        let mut dark = Moments::<4>::default();
        dark.pad_to(MAX_FINITE_DETECTOR_SAMPLES).unwrap();
        assert_eq!(dark.count(), MAX_FINITE_DETECTOR_SAMPLES);
        assert_eq!(dark.mean(), [0.0; 4]);
        assert_eq!(dark.standard_error().unwrap(), [0.0; 4]);
    }

    #[test]
    fn g3_signed_large_stokes_components_retain_finite_standard_errors() {
        let mut moments = Moments::<2>::default();
        moments.add_at([2e200, 2e200], 0).unwrap();
        moments.add_at([2e200, -2e200], 2).unwrap();
        close(moments.mean()[0], (4.0 / 3.0) * 1e200, 1e185);
        close(moments.mean()[1], 0.0, 1e185);
        close(
            moments.standard_error().unwrap()[0],
            (2.0 / 3.0) * 1e200,
            1e185,
        );
        close(
            moments.standard_error().unwrap()[1],
            (2.0 / 3.0_f64.sqrt()) * 1e200,
            1e185,
        );
    }

    #[test]
    fn g4_invalid_order_values_and_budget_leave_accumulated_pixels_unchanged() {
        let mut moments = Moments::<2>::default();
        assert_eq!(
            moments.standard_error(),
            Err(FiniteDetectorError::InvalidRequest)
        );
        moments.add_at([1.0, -1.0], 0).unwrap();
        assert!(moments.standard_error().is_err());
        let original = moments;
        for (values, index) in [
            ([2.0, 3.0], 0),
            ([2.0, f64::NAN], 5),
            ([f64::INFINITY, 3.0], 10),
            ([2.0, 3.0], MAX_FINITE_DETECTOR_SAMPLES),
        ] {
            assert!(moments.add_at(values, index).is_err());
            assert_eq!(moments, original);
        }
        for count in [0, MAX_FINITE_DETECTOR_SAMPLES + 1] {
            assert!(moments.pad_to(count).is_err());
            assert_eq!(moments, original);
        }
        moments.pad_to(1).unwrap();
        assert_eq!(moments, original);
        moments
            .add_at([0.0; 2], MAX_FINITE_DETECTOR_SAMPLES - 1)
            .unwrap();
        assert_eq!(moments.count(), MAX_FINITE_DETECTOR_SAMPLES);
    }
}
