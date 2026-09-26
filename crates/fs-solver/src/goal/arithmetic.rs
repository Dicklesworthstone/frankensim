//! Directed binary64 operations for the sparse residual evaluator. The
//! rounding vocabulary is fs-math's audited facade, not an ambient rounding
//! mode. FMA performs one correctly rounded operation before outward nudging.

use super::{GoalResidualError, ScalarEnclosure};
use fs_math::{next_down, next_up};

pub(super) fn finite(value: f64) -> Result<f64, GoalResidualError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(GoalResidualError::ArithmeticRange)
    }
}

pub(super) fn up(value: f64) -> Result<f64, GoalResidualError> {
    finite(next_up(finite(value)?))
}

pub(super) fn down(value: f64) -> Result<f64, GoalResidualError> {
    finite(next_down(finite(value)?))
}

pub(super) fn add_up(a: f64, b: f64) -> Result<f64, GoalResidualError> {
    if a == 0.0 {
        return Ok(b);
    }
    if b == 0.0 {
        return Ok(a);
    }
    up(a + b)
}

pub(super) fn mul_up(a: f64, b: f64) -> Result<f64, GoalResidualError> {
    if a == 0.0 || b == 0.0 {
        return Ok(0.0);
    }
    up(a * b)
}

impl ScalarEnclosure {
    pub(super) const fn point(value: f64) -> Self {
        Self {
            lower: value,
            upper: value,
        }
    }

    pub(super) fn add_scaled(
        self,
        value: Self,
        scale: f64,
    ) -> Result<Self, GoalResidualError> {
        if scale == 0.0 || (value.lower == 0.0 && value.upper == 0.0) {
            return Ok(self);
        }
        let (lo, hi) = if scale > 0.0 {
            (value.lower, value.upper)
        } else {
            (value.upper, value.lower)
        };
        Ok(Self {
            lower: down(scale.mul_add(lo, self.lower))?,
            upper: up(scale.mul_add(hi, self.upper))?,
        })
    }

    pub(super) fn widen(self, radius: f64) -> Result<Self, GoalResidualError> {
        if radius == 0.0 {
            return Ok(self);
        }
        Ok(Self {
            lower: down(self.lower - radius)?,
            upper: up(self.upper + radius)?,
        })
    }

    pub(super) fn radius_about(self, value: f64) -> Result<f64, GoalResidualError> {
        if self.lower == value && self.upper == value {
            return Ok(0.0);
        }
        up((self.lower - value).abs().max((self.upper - value).abs()))
    }
}
