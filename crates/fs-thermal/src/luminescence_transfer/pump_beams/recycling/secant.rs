//! Bounded inverse-Broyden history for the normalized emitted-source closure.
//!
//! The inverse is I + sum(u v^T). Applying or updating it takes O(n m)
//! arithmetic and storage for n radiating states and at most m retained pairs.
//! A full history restarts BEFORE the next correction is formed, so the newest
//! secant equation still holds. This changes acceleration, never admission.

use super::{SingletTripletError, Sum, admit_storage, poll};

struct Correction {
    left: Vec<f64>,
    right: Vec<f64>,
}

pub(super) struct SecantInverse {
    dimension: usize,
    maximum_history: usize,
    history: Vec<Correction>,
    /// Conservative element bound; avoids constructing a dense matrix merely
    /// to reject a nonfinite update.
    element_bound: f64,
}

impl SecantInverse {
    pub(super) fn try_new(
        dimension: usize,
        maximum_history: usize,
    ) -> Result<Self, SingletTripletError> {
        let elements = dimension
            .checked_mul(maximum_history)
            .and_then(|count| count.checked_mul(2))
            .ok_or(SingletTripletError::Unrepresentable(
                "photon secant history size",
            ))?;
        admit_storage::<f64>(elements, "photon secant history storage")?;
        admit_storage::<f64>(dimension, "photon secant vector storage")?;
        admit_storage::<Correction>(maximum_history, "photon secant history entries")?;
        Ok(Self {
            dimension,
            maximum_history,
            history: Vec::new(),
            element_bound: 1.0,
        })
    }

    pub(super) fn reset(&mut self) {
        self.history.clear();
        self.element_bound = 1.0;
    }

    pub(super) fn apply(
        &self,
        input: &[f64],
        output: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), SingletTripletError> {
        debug_assert_eq!(input.len(), self.dimension);
        debug_assert_eq!(output.len(), self.dimension);
        poll(cancelled)?;
        let mut projections = Vec::with_capacity(self.history.len());
        for correction in &self.history {
            let mut product = Sum::default();
            for (index, (&right, &value)) in correction.right.iter().zip(input).enumerate() {
                if index % 256 == 0 {
                    poll(cancelled)?;
                }
                product.add(right * value)?;
            }
            projections.push(product.value()?);
        }
        for row in 0..self.dimension {
            if row % 256 == 0 {
                poll(cancelled)?;
            }
            let mut sum = Sum::default();
            sum.add(input[row])?;
            for (correction, &projection) in self.history.iter().zip(&projections) {
                sum.add(correction.left[row] * projection)?;
            }
            output[row] = sum.value()?;
        }
        Ok(())
    }

    pub(super) fn update(
        &mut self,
        displacement: &[f64],
        change: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), SingletTripletError> {
        debug_assert_eq!(displacement.len(), self.dimension);
        debug_assert_eq!(change.len(), self.dimension);
        poll(cancelled)?;
        if self.maximum_history == 0 {
            return Ok(());
        }
        let mut denominator = 0.0;
        for (index, &value) in change.iter().enumerate() {
            if index % 256 == 0 {
                poll(cancelled)?;
            }
            denominator += value * value;
        }
        if !denominator.is_finite() || denominator <= f64::MIN_POSITIVE {
            return Ok(());
        }
        if self.history.len() == self.maximum_history {
            self.reset();
        }
        let mut left = vec![0.0; self.dimension];
        self.apply(change, &mut left, cancelled)?;
        let mut right = Vec::with_capacity(self.dimension);
        // (s - H y) y^T / (y^T y), with balanced factors. Normalizing both
        // factors by ||y|| avoids the needlessly large y/(y^T y) vector.
        let norm = denominator.sqrt();
        let mut maximum_left = 0.0_f64;
        let mut maximum_right = 0.0_f64;
        for index in 0..self.dimension {
            if index % 256 == 0 {
                poll(cancelled)?;
            }
            left[index] = (displacement[index] - left[index]) / norm;
            let value = change[index] / norm;
            if !left[index].is_finite() || !value.is_finite() {
                self.reset();
                return Ok(());
            }
            maximum_left = maximum_left.max(left[index].abs());
            maximum_right = maximum_right.max(value.abs());
            right.push(value);
        }
        let bound = self.element_bound + maximum_left * maximum_right;
        if !bound.is_finite() {
            self.reset();
            return Ok(());
        }
        self.element_bound = bound;
        self.history.push(Correction { left, right });
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pump_recycling_secant_matches_dense_updates_and_bounded_restarts() {
        const N: usize = 4;
        const HISTORY: usize = 3;
        let mut inverse = SecantInverse::try_new(N, HISTORY).unwrap();
        let mut dense = [[0.0; N]; N];
        let samples = [
            ([0.7, 0.9, -0.1, 0.5], [1.0, 0.5, -0.3, 0.8]),
            ([0.4, -0.2, 0.6, 1.1], [0.2, -0.7, 0.8, 0.5]),
            ([-0.6, 0.4, 1.2, 0.3], [-0.5, 0.3, 0.9, -0.2]),
            ([0.5, 0.2, -0.4, 0.6], [0.9, -0.1, -0.6, 0.4]),
            ([-0.3, 0.8, 0.2, -0.5], [-0.4, 0.7, 0.1, -0.6]),
            ([0.2, -0.4, 0.3, 0.9], [0.1, -0.5, 0.6, 0.7]),
        ];
        for (ordinal, (displacement, change)) in samples.into_iter().enumerate() {
            if ordinal % HISTORY == 0 {
                dense = std::array::from_fn(|row| {
                    std::array::from_fn(|column| if row == column { 1.0 } else { 0.0 })
                });
            }
            let denominator: f64 = change.iter().map(|value| value * value).sum();
            let correction: [f64; N] = std::array::from_fn(|row| {
                displacement[row]
                    - dense[row]
                        .iter()
                        .zip(change)
                        .map(|(&coefficient, value)| coefficient * value)
                        .sum::<f64>()
            });
            for row in 0..N {
                for column in 0..N {
                    dense[row][column] += correction[row] * change[column] / denominator;
                }
            }
            inverse
                .update(&displacement, &change, &mut || false)
                .unwrap();
            assert_eq!(inverse.history.len(), ordinal % HISTORY + 1);
            for input in [[1.0, -0.5, 0.2, 0.8], [0.0, 1.0, 0.0, 0.0], change] {
                let mut actual = [0.0; N];
                inverse.apply(&input, &mut actual, &mut || false).unwrap();
                for row in 0..N {
                    let expected: f64 = dense[row]
                        .iter()
                        .zip(input)
                        .map(|(&coefficient, value)| coefficient * value)
                        .sum();
                    assert!((actual[row] - expected).abs() < 1e-12);
                    if input == change {
                        assert!((actual[row] - displacement[row]).abs() < 1e-12);
                    }
                }
            }
        }

        let input = [0.3, 0.7, -0.2, 0.5];
        let mut actual = [0.0; N];
        let mut plain = SecantInverse::try_new(N, 0).unwrap();
        plain
            .update(&samples[0].0, &samples[0].1, &mut || false)
            .unwrap();
        plain.apply(&input, &mut actual, &mut || false).unwrap();
        assert_eq!(actual, input);
        assert!(plain.history.is_empty());

        inverse
            .update(&[f64::MAX; N], &[1e-100; N], &mut || false)
            .unwrap();
        assert!(inverse.history.is_empty());
        inverse.apply(&input, &mut actual, &mut || false).unwrap();
        assert_eq!(actual, input);
        assert!(inverse.apply(&input, &mut actual, &mut || true).is_err());
        assert!(SecantInverse::try_new(usize::MAX, HISTORY).is_err());
    }
}
