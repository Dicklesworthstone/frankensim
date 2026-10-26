//! Reciprocal pressure/wall coordinates for an already assembled gas mass.
//!
//! G is the exact geometric polynomial Gram matrix, independent of the gas.
//! With wall mass P*sigma*G, C = sqrt(P/sigma) Mp^(-1/2) G^(1/2).
//! C is generally NOT symmetric: pressure receives -C, velocity receives +C^T.
//! No nodal gas sample or mean bulk modulus enters this reduction.

use fs_math::det;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum CoordinateError {
    Invalid,
    Cancelled,
}

#[derive(Debug)]
pub(super) struct WallCoordinates {
    local: usize,
    coupling: Vec<f64>,
    pub(super) velocity_map: Vec<f64>,
}

impl WallCoordinates {
    pub(super) fn new(
        local: usize,
        pressure_inv_sqrt: &[f64],
        gram_sqrt: &[f64],
        gram_inv_sqrt: &[f64],
        perimeter: f64,
        surface_density: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, CoordinateError> {
        if cancelled() {
            return Err(CoordinateError::Cancelled);
        }
        let entries = local.checked_mul(local).ok_or(CoordinateError::Invalid)?;
        if local == 0
            || [pressure_inv_sqrt, gram_sqrt, gram_inv_sqrt]
                .iter()
                .any(|a| a.len() != entries || a.iter().any(|v| !v.is_finite()))
            || !perimeter.is_finite()
            || perimeter <= 0.0
            || !surface_density.is_finite()
            || surface_density <= 0.0
        {
            return Err(CoordinateError::Invalid);
        }
        let coupling_scale = det::sqrt(perimeter / surface_density);
        let velocity_scale = 1.0 / det::sqrt(perimeter * surface_density);
        if [coupling_scale, velocity_scale]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(CoordinateError::Invalid);
        }
        let mut coupling = vec![0.0; entries];
        let mut velocity_map = vec![0.0; entries];
        for i in 0..local {
            if cancelled() {
                return Err(CoordinateError::Cancelled);
            }
            for j in 0..local {
                let mut value = 0.0;
                for k in 0..local {
                    value += pressure_inv_sqrt[i * local + k] * gram_sqrt[k * local + j];
                }
                coupling[i * local + j] = coupling_scale * value;
                velocity_map[i * local + j] = velocity_scale * gram_inv_sqrt[i * local + j];
            }
            if coupling[i * local..(i + 1) * local]
                .iter()
                .all(|v| *v == 0.0)
                || velocity_map[i * local + i] <= 0.0
            {
                return Err(CoordinateError::Invalid);
            }
        }
        if coupling.iter().chain(&velocity_map).any(|v| !v.is_finite()) {
            return Err(CoordinateError::Invalid);
        }
        Ok(Self {
            local,
            coupling,
            velocity_map,
        })
    }

    /// Stamp only an unpublished candidate. Cancellation discards that candidate
    /// in the owning duct constructor, so no accepted state/operator is changed.
    pub(super) fn stamp(
        &self,
        operator: &mut [f64],
        total: usize,
        pressure_start: usize,
        velocity_start: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), CoordinateError> {
        if total.checked_mul(total) != Some(operator.len())
            || pressure_start
                .checked_add(self.local)
                .is_none_or(|end| end > total)
            || velocity_start
                .checked_add(self.local)
                .is_none_or(|end| end > total)
            || !(pressure_start + self.local <= velocity_start
                || velocity_start + self.local <= pressure_start)
        {
            return Err(CoordinateError::Invalid);
        }
        for i in 0..self.local {
            if cancelled() {
                return Err(CoordinateError::Cancelled);
            }
            for j in 0..self.local {
                let pressure = pressure_start + i;
                let velocity = velocity_start + j;
                let value = self.coupling[i * self.local + j];
                operator[pressure * total + velocity] = -value;
                operator[velocity * total + pressure] = value;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "coordinates/tests.rs"]
mod tests;
