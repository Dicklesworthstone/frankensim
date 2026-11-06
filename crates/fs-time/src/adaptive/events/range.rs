//! Floating-point Bernstein hulls of the original accepted dense polynomial.

use super::{DenseStep, EventError};
use fs_math::dd::Dd;

impl DenseStep {
    /// Bound one component of the numerical interpolant from its start through
    /// `through_time`, including the actual accepted endpoint when applicable.
    ///
    /// Restricting the quartic to this prefix and expressing it in Bernstein
    /// form gives a convex hull of five coefficients. It therefore includes
    /// internal extrema without sampling or another root search. It may be
    /// wider than the actual range; callers can refine the original time step.
    ///
    /// The coefficients are accumulated with double-double arithmetic and
    /// rounded outward when their retained tail is nonzero. This is a numerical
    /// hull, not an interval certificate for arithmetic or for the exact ODE
    /// solution. The dense extension retains its original accuracy limitations.
    /// No extrapolation or allocation is performed.
    pub fn component_range(
        &self,
        component: usize,
        through_time: f64,
    ) -> Result<[f64; 2], EventError> {
        if component >= self.initial.len()
            || !through_time.is_finite()
            || through_time < self.start
            || through_time > self.end
        {
            return Err(EventError::InvalidSample);
        }
        let initial = self.initial[component];
        if through_time == self.start {
            return Ok([initial; 2]);
        }
        let theta = Dd::from_f64((through_time - self.start) / self.h);
        let mut power = Dd::from_f64(self.h);
        let mut c = [Dd::ZERO; 4];
        for (i, coefficient) in self.coefficients[component].iter().enumerate() {
            power = power * theta;
            c[i] = power * Dd::from_f64(*coefficient);
        }
        let y = Dd::from_f64(initial);
        let quarter = Dd::from_f64(0.25);
        let half = Dd::from_f64(0.5);
        let b = [
            y,
            y + c[0] * quarter,
            y + c[0] * half + c[1] / Dd::from_f64(6.0),
            y + c[0] * Dd::from_f64(0.75) + c[1] * half + c[2] * quarter,
            y + c[0] + c[1] + c[2] + c[3],
        ];
        let mut range = [initial; 2];
        for value in b {
            if !value.hi.is_finite() || !value.lo.is_finite() {
                return Err(EventError::NonFiniteInterpolation);
            }
            let rounded = value.to_f64();
            let bounds = if value.lo == 0.0 {
                [rounded; 2]
            } else {
                [fs_math::next_down(rounded), fs_math::next_up(rounded)]
            };
            if bounds.iter().any(|value| !value.is_finite()) {
                return Err(EventError::NonFiniteInterpolation);
            }
            range[0] = range[0].min(bounds[0]);
            range[1] = range[1].max(bounds[1]);
        }
        if through_time == self.end {
            range[0] = range[0].min(self.final_state[component]);
            range[1] = range[1].max(self.final_state[component]);
        }
        Ok(range)
    }
}
