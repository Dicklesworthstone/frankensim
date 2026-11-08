//! Explicit post-initiation attenuation in the original equivalent-work coordinate.
//!
//! `x = h/Gp * integral R(alpha) d(alpha)` is dimensionless and nondecreasing.
//! Both profiles have `integral_0^infinity g(x) dx = 1`, so Gp retains its
//! meaning as DEGRADED equivalent plastic-work capacity. It is not total
//! thermodynamic fracture dissipation, current crack area, or a deletion rule.
//! Selecting a shape is an explicit model assumption, not an inference from Gp.

use fs_math::det;

/// Authored post-initiation shape with the same original equivalent-work capacity.
/// Existing constructors select `Exponential`; a different profile must be
/// selected before an immutable material owner lends any accepted histories.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum DuctileWorkEvolution {
    /// `g=exp(-x)`; positive at every finite representable supported coordinate.
    /// Unrepresentable integrity still refuses instead of imposing a cutoff.
    #[default]
    Exponential,
    /// `g=max(1-x/3,0)^2`. Integrity and its first derivative reach zero
    /// continuously at x=3, after exactly Gp of integrated degraded equivalent
    /// work. This is an authored compact-support polynomial, not a fitted alloy
    /// law or permission to erase material, mass, caloric energy, or histories.
    FiniteSupportQuadratic,
}

/// Stable increments in damage and in the dimensionless primitive of g.
/// These are distinct for the finite-support shape.
pub(crate) struct WorkIncrement {
    pub damage: f64,
    pub degraded_work: f64,
}

impl DuctileWorkEvolution {
    // Native owners admit finite nonnegative coordinates before these kernels.
    pub(crate) fn integrity(self, x: f64) -> f64 {
        match self {
            Self::Exponential => det::exp(-x),
            Self::FiniteSupportQuadratic => {
                let remaining = (3.0 - x).max(0.0) / 3.0;
                remaining * remaining
            }
        }
    }

    pub(crate) fn damage(self, x: f64) -> f64 {
        match self {
            Self::Exponential => -det::expm1(-x),
            Self::FiniteSupportQuadratic if x >= 3.0 => 1.0,
            // x*(6-x)/9 avoids cancellation and division of a tiny x first.
            Self::FiniteSupportQuadratic => x * ((6.0 - x) / 9.0),
        }
    }

    /// Nonnegative -dg/dx, including the exact zero branch at finite failure.
    pub(crate) fn degradation_slope(self, x: f64) -> f64 {
        match self {
            Self::Exponential => self.integrity(x),
            Self::FiniteSupportQuadratic => (2.0 / 3.0) * ((3.0 - x).max(0.0) / 3.0),
        }
    }

    /// Phi(x)=integral_0^x g(s) ds, evaluated without subtracting nearby ones.
    pub(crate) fn degraded_work(self, x: f64) -> f64 {
        match self {
            Self::Exponential => self.damage(x),
            Self::FiniteSupportQuadratic if x >= 3.0 => 1.0,
            Self::FiniteSupportQuadratic if x <= 1.5 => x * (1.0 - x / 3.0 + x * x / 27.0),
            Self::FiniteSupportQuadratic => {
                let remaining = (3.0 - x) / 3.0;
                1.0 - remaining * remaining * remaining
            }
        }
    }

    /// Integrate from x through x+increment, including a terminal crossing.
    /// The original coordinate is never capped: only the specified profile's
    /// zero tail stops contributing damage or degraded equivalent work.
    pub(crate) fn increment(self, x: f64, increment: f64) -> WorkIncrement {
        match self {
            Self::Exponential => {
                let damage = self.integrity(x) * -det::expm1(-increment);
                WorkIncrement {
                    damage,
                    degraded_work: damage,
                }
            }
            Self::FiniteSupportQuadratic => {
                let remaining = (3.0 - x).max(0.0);
                let before = remaining / 3.0;
                if increment >= remaining {
                    // Consume the exact surviving tail, even if the authored
                    // interval crosses far beyond the finite support.
                    WorkIncrement {
                        damage: before * before,
                        degraded_work: before * before * before,
                    }
                } else {
                    let after = (remaining - increment) / 3.0;
                    WorkIncrement {
                        damage: increment * ((before + after) / 3.0),
                        degraded_work: increment
                            * ((before * before + before * after + after * after) / 3.0),
                    }
                }
            }
        }
    }

    pub(crate) fn is_fully_failed(self, x: f64) -> bool {
        matches!(self, Self::FiniteSupportQuadratic) && x >= 3.0
    }
}
