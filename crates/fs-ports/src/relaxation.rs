//! Bounded vector Aitken relaxation for dimensioned coupling coordinates.
//!
//! For r_k = solved_k - applied_k, use
//! omega_k = -omega_(k-1) (r_(k-1) . (r_k-r_(k-1))) / |r_k-r_(k-1)|^2.
//! See Küttler & Wall, Computational Mechanics 43 (2008), 61–72,
//! doi:10.1007/s00466-008-0255-5. This selects an iteration weight;
//! it supplies neither a constitutive law nor a physical acceptance test.

/// Invalid relaxation controls or a changed/nonfinite coordinate vector.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AitkenRelaxationError {
    /// Require finite 0 < minimum <= initial <= maximum <= 1.
    InvalidWeights,
    /// Every applied and solved coordinate must be finite.
    NonFiniteCoordinate,
    /// One history must retain the same coordinate count and ordering.
    CoordinateCountChanged,
}

impl core::fmt::Display for AitkenRelaxationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "bounded Aitken relaxation: {self:?}")
    }
}

impl std::error::Error for AitkenRelaxationError {}

/// Candidate-local vector relaxation with independently scaled coordinates.
///
/// The same finite, ordered (applied, solved) coordinates must be supplied on
/// every trial. Each scale is the largest magnitude seen for that coordinate.
/// Old and new residuals use the same current metric, so changing independent
/// coordinate units does not choose the relaxation weight. Updates stay convex.
///
/// Storage and synchronous work are linear in the caller's coordinate count.
/// Physical owners retain their own work/cancellation limits and must validate
/// the actual constitutive/energy residuals. Convergence is not guaranteed.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundedAitkenRelaxation {
    weight: f64,
    minimum: f64,
    maximum: f64,
    dimension: Option<usize>,
    scales: Vec<f64>,
    previous: Vec<f64>,
    residuals: Vec<f64>,
}

impl BoundedAitkenRelaxation {
    /// Admit positive, finite, ordered weights bounded above by one.
    pub fn new(
        initial_weight: f64,
        minimum_weight: f64,
        maximum_weight: f64,
    ) -> Result<Self, AitkenRelaxationError> {
        if ![initial_weight, minimum_weight, maximum_weight].iter().all(|x| x.is_finite())
            || minimum_weight <= 0.0
            || minimum_weight > initial_weight
            || initial_weight > maximum_weight
            || maximum_weight > 1.0
        {
            return Err(AitkenRelaxationError::InvalidWeights);
        }
        Ok(Self {
            weight: initial_weight,
            minimum: minimum_weight,
            maximum: maximum_weight,
            dimension: None,
            scales: Vec::new(),
            previous: Vec::new(),
            residuals: Vec::new(),
        })
    }

    /// Select the next weight without changing any caller-owned coordinates.
    ///
    /// The first call returns the initial weight. Stagnant or numerically
    /// unresolved secants retain the preceding weight. Invalid inputs leave
    /// the complete history unchanged. Zero and signed coordinates are valid.
    pub fn next_weight(
        &mut self,
        coordinates: &[(f64, f64)],
    ) -> Result<f64, AitkenRelaxationError> {
        if coordinates.iter().any(|(a, b)| !a.is_finite() || !b.is_finite()) {
            return Err(AitkenRelaxationError::NonFiniteCoordinate);
        }
        if self.dimension.is_some_and(|n| n != coordinates.len()) {
            return Err(AitkenRelaxationError::CoordinateCountChanged);
        }
        let first = self.dimension.is_none();
        if first {
            self.dimension = Some(coordinates.len());
            self.scales.resize(coordinates.len(), 0.0);
            self.previous.resize(coordinates.len(), 0.0);
        }
        self.residuals.clear();
        let mut magnitude = 0.0_f64;
        for (index, &(applied, solved)) in coordinates.iter().enumerate() {
            let old_scale = self.scales[index];
            let scale = old_scale.max(applied.abs()).max(solved.abs());
            let residual = if scale == 0.0 {
                0.0
            } else {
                self.previous[index] *= old_scale / scale;
                // Preserve nearby differences; normalize first only when
                // opposite signed finite coordinates overflow subtraction.
                let difference = solved - applied;
                if difference.is_finite() {
                    difference / scale
                } else {
                    solved / scale - applied / scale
                }
            };
            self.scales[index] = scale;
            magnitude = magnitude.max(residual.abs()).max(self.previous[index].abs());
            self.residuals.push(residual);
        }
        if !first && magnitude > 0.0 {
            let mut numerator = 0.0;
            let mut denominator = 0.0;
            for (&previous, &residual) in self.previous.iter().zip(&self.residuals) {
                // Common scaling avoids underflow of squared small residuals.
                let previous = previous / magnitude;
                let difference = residual / magnitude - previous;
                numerator += previous * difference;
                denominator += difference * difference;
            }
            if denominator > 0.0 {
                let proposed = -self.weight * (numerator / denominator);
                if proposed.is_finite() {
                    self.weight = proposed.clamp(self.minimum, self.maximum);
                }
            }
        }
        std::mem::swap(&mut self.previous, &mut self.residuals);
        Ok(self.weight)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_affine_oscillation_recovers_the_exact_relaxation() {
        // F(x)=1-9x has x*=0.1; plain fixed point diverges.
        let mut acceleration = BoundedAitkenRelaxation::new(1.0, 1e-4, 1.0).unwrap();
        assert_eq!(acceleration.next_weight(&[(0.0, 1.0)]).unwrap(), 1.0);
        let weight = acceleration.next_weight(&[(1.0, -8.0)]).unwrap();
        assert!((weight - 0.1).abs() < 1e-15);
        assert!((1.0 + weight * (-8.0 - 1.0) - 0.1).abs() < 1e-15);
    }

    #[test]
    fn g1_passive_resistor_map_reaches_its_independent_operating_point() {
        // rho(J)=1, then 40|J|-35, then 9, with breaks at 0.9 and 1.1.
        // At E=+/-5 the unique passive solution is J=+/-1, rho=5.
        for sign in [-1.0, 1.0] {
            let mut acceleration = BoundedAitkenRelaxation::new(1.0, 1e-4, 1.0).unwrap();
            let mut current = 0.0_f64;
            let mut converged = false;
            for _ in 0..48 {
                let rho = if current.abs() <= 0.9 { 1.0 }
                    else if current.abs() >= 1.1 { 9.0 }
                    else { 40.0 * current.abs() - 35.0 };
                let solved = sign * 5.0 / rho;
                if (solved - current).abs() < 1e-12 {
                    converged = true;
                    break;
                }
                let weight = acceleration.next_weight(&[(current, solved)]).unwrap();
                current += weight * (solved - current);
            }
            assert!(converged);
            assert!((current - sign).abs() < 1e-12);
        }
    }

    #[test]
    fn g3_coordinate_units_and_polarity_do_not_choose_the_weight() {
        let pairs = [
            [(0.0, 1.0), (2.0, 0.0)],
            [(1.0, -8.0), (0.0, 3.0)],
            [(0.1, 0.2), (0.3, 0.4)],
        ];
        let mut reference = BoundedAitkenRelaxation::new(0.8, 1e-5, 1.0).unwrap();
        let expected: Vec<_> = pairs.iter()
            .map(|pair| reference.next_weight(pair).unwrap()).collect();
        for scales in [[1e-200, 1e200], [-1e200, -1e-200]] {
            let mut acceleration = BoundedAitkenRelaxation::new(0.8, 1e-5, 1.0).unwrap();
            for (pair, &expected) in pairs.iter().zip(&expected) {
                let scaled: Vec<_> = pair.iter().zip(scales)
                    .map(|(&(a, b), scale)| (a * scale, b * scale)).collect();
                assert!((acceleration.next_weight(&scaled).unwrap() - expected).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn g0_zero_stagnant_and_extreme_coordinates_keep_bounded_weights() {
        let mut acceleration = BoundedAitkenRelaxation::new(0.5, 0.01, 0.8).unwrap();
        assert_eq!(acceleration.next_weight(&[(0.0, 0.0)]).unwrap(), 0.5);
        assert!(acceleration.next_weight(&[(0.0, 1.0)]).unwrap() >= 0.01);
        for pair in [(f64::MAX, -f64::MAX), (-f64::MAX, f64::MAX), (1.0, 1.0)] {
            assert!((0.01..=0.8).contains(&acceleration.next_weight(&[pair]).unwrap()));
        }
        let mut stagnant = BoundedAitkenRelaxation::new(0.4, 0.01, 1.0).unwrap();
        for _ in 0..4 {
            assert_eq!(stagnant.next_weight(&[(0.0, 1.0)]).unwrap(), 0.4);
        }
        let mut empty = BoundedAitkenRelaxation::new(0.5, 0.01, 1.0).unwrap();
        for _ in 0..2 { assert_eq!(empty.next_weight(&[]).unwrap(), 0.5); }
    }

    #[test]
    fn g0_invalid_controls_and_changed_coordinates_refuse_atomically() {
        for weights in [(0.0, 0.0, 1.0), (0.1, 0.2, 1.0), (0.5, 0.1, 0.4),
            (0.5, 0.1, 1.1), (f64::NAN, 0.1, 1.0), (0.5, 0.1, f64::INFINITY)] {
            assert_eq!(BoundedAitkenRelaxation::new(weights.0, weights.1, weights.2),
                Err(AitkenRelaxationError::InvalidWeights));
        }
        let mut acceleration = BoundedAitkenRelaxation::new(0.5, 0.01, 1.0).unwrap();
        acceleration.next_weight(&[(0.0, 1.0)]).unwrap();
        let before = acceleration.clone();
        assert_eq!(acceleration.next_weight(&[(0.0, f64::NAN)]),
            Err(AitkenRelaxationError::NonFiniteCoordinate));
        assert_eq!(acceleration, before);
        assert_eq!(acceleration.next_weight(&[]),
            Err(AitkenRelaxationError::CoordinateCountChanged));
        assert_eq!(acceleration, before);
    }
}
