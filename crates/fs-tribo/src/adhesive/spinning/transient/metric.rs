//! Small original-matrix metric and conservative finite-extinction bounds.

use super::*;
use crate::sliding_patch::contact::{Factor, factor_circular_response};

pub(super) struct Metric {
    factor: Factor,
    lower: [[f64; 3]; 3],
    root_scale: f64,
    root_lower: f64,
    root_upper: f64,
    capacity: f64,
    radial_error: f64,
    row_bounds: [f64; 3],
    pub minimum_decay: f64,
    pub maximum_speed: f64,
}
impl Metric {
    pub fn new(
        matrix: [[f64; 3]; 3],
        capacity: f64,
        power_norm_bounds: [f64; 2],
        config: AdhesiveWrenchConfig,
    ) -> Result<Self, Error> {
        if (0..3).any(|i| (0..i).any(|j| matrix[i][j] != matrix[j][i])) {
            return Err(Error::Invalid(
                "transient physical mobility must be symmetric",
            ));
        }
        let factor = factor_circular_response(matrix, config)?;
        let (lower, scale) = factor.normalized_lower();
        let normalized = matrix.map(|row| row.map(|v| v / scale));
        let mut inverse = [[0.0; 3]; 3];
        for j in 0..3 {
            let mut unit = [0.0; 3];
            unit[j] = 1.0;
            let column = factor
                .solve(unit)
                .ok_or(Error::Invalid("unrepresentable mobility inverse"))?;
            for i in 0..3 {
                inverse[i][j] = finite(column[i] * scale)?;
            }
        }
        let upper = row_norm(normalized);
        let inverse_norm = row_norm(inverse);
        let inverse_residual = row_norm(core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                (if i == j { 1.0 } else { 0.0 })
                    - (0..3)
                        .map(|k| normalized[i][k] * inverse[k][j])
                        .sum::<f64>()
            })
        })) + 64.0 * f64::EPSILON * upper * inverse_norm;
        if !inverse_residual.is_finite() || inverse_residual >= 0.25 {
            return Err(Error::Invalid(
                "mobility inverse conditioning prevents a bounded tail",
            ));
        }
        // ||W^-1|| <= ||B||/(1-||I-WB||), without assuming Gershgorin
        // positivity of W itself. Its upper norm bound gives lambda_min.
        let lower_eigenvalue = (1.0 - inverse_residual) / inverse_norm;
        let factor_residual = row_norm(core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                normalized[i][j] - (0..3).map(|k| lower[i][k] * lower[j][k]).sum::<f64>()
            })
        })) + 64.0 * f64::EPSILON * upper;
        let radial_error = factor_residual / lower_eigenvalue;
        if !radial_error.is_finite() || radial_error >= 0.25 || lower_eigenvalue <= 0.0 {
            return Err(Error::Invalid(
                "mobility factor conditioning prevents a bounded tail",
            ));
        }
        let root_scale = det::sqrt(scale);
        let root_lower = finite(root_scale * det::sqrt(lower_eigenvalue))?;
        let root_upper = finite(root_scale * det::sqrt(upper * (1.0 + 64.0 * f64::EPSILON)))?;
        let minimum_decay = finite(power_norm_bounds[0] * capacity * root_lower)?;
        let maximum_speed = finite(power_norm_bounds[1] * capacity * root_upper)?;
        if minimum_decay < 0.0
            || (power_norm_bounds[0] > 0.0 && minimum_decay == 0.0)
            || maximum_speed <= 0.0
        {
            return Err(Error::Invalid("unrepresentable positive extinction bounds"));
        }
        Ok(Self {
            factor,
            lower,
            root_scale,
            root_lower,
            root_upper,
            capacity,
            radial_error,
            row_bounds: core::array::from_fn(|i| {
                det::sqrt(matrix[i][i]) * (1.0 + 64.0 * f64::EPSILON)
            }),
            minimum_decay,
            maximum_speed,
        })
    }
    pub fn radius(&self, q: [f64; 3]) -> Result<(f64, f64), Error> {
        let radius = norm(
            self.factor
                .mass_coordinates(q)
                .ok_or(Error::Invalid("unrepresentable mass coordinates"))?,
        );
        Ok((
            finite(radius * det::sqrt(1.0 - self.radial_error))?,
            finite(radius * det::sqrt(1.0 + self.radial_error))?,
        ))
    }
    pub fn tail(
        &self,
        q: [f64; 3],
        conditional: Option<CircularTransientExtinction>,
    ) -> Result<Option<([f64; 2], [f64; 3])>, Error> {
        let (lower, upper) = self.radius(q)?;
        let minimum_decay = if self.minimum_decay > 0.0 {
            self.minimum_decay
        } else if let Some(bound) = conditional
            && upper <= bound.maximum_mass_radius
        {
            let decay = finite(bound.minimum_power_norm * self.capacity * self.root_lower)?;
            if decay <= 0.0 {
                return Err(Error::Invalid(
                    "unrepresentable conditional extinction bound",
                ));
            }
            decay
        } else {
            return Ok(None);
        };
        let time = [
            finite(lower / self.maximum_speed)?,
            finite(upper / minimum_decay)?,
        ];
        let integrated_radius = finite(0.5 * upper * time[1])?;
        let travel = self.row_bounds.map(|v| v * integrated_radius);
        for value in travel {
            finite(value)?;
        }
        Ok(Some((time, travel)))
    }
    pub fn safe_step(&self, q: [f64; 3], linear_bound: Option<f64>) -> Result<f64, Error> {
        let (lower, upper) = self.radius(q)?;
        let mut speed = self.maximum_speed;
        if let Some(bound) = linear_bound {
            // |w|<=b|q| implies |y'|<=lambda_max(W)*b*|y| in mass
            // coordinates. Every positive stage stays in the 1.25-radius
            // ball under this step; unlike a constant force bound, the local
            // rate bound shrinks with genuine asymptotic viscous motion.
            let local = (1.25 * self.root_upper * bound) * (self.root_upper * upper);
            if local > 0.0 {
                speed = speed.min(local);
            }
        }
        Ok(0.25 * lower / speed)
    }
    pub fn stop_impulse(&self, q: [f64; 3]) -> Result<[f64; 3], Error> {
        self.factor
            .solve(q.map(|v| -v))
            .ok_or(Error::Invalid("unrepresentable stopping impulse"))
    }
    pub fn rate_error(&self, impulse: [f64; 3]) -> Result<f64, Error> {
        let mass_change = core::array::from_fn(|j| {
            self.root_scale * (0..3).map(|i| self.lower[i][j] * impulse[i]).sum::<f64>()
        });
        // The computed factor obeys L L^T >= (1-eta) W. Inflate its
        // impulse-induced mass change before bounding the original W rate.
        finite(self.root_upper * norm(mass_change) / det::sqrt(1.0 - self.radial_error))
    }
}
fn row_norm(matrix: [[f64; 3]; 3]) -> f64 {
    matrix
        .map(|row| row.into_iter().map(f64::abs).sum::<f64>())
        .into_iter()
        .fold(0.0, f64::max)
}
