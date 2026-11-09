//! Sensor-protected dynamic WALE with the existing Favre--Germano filters.
//!
//! For `G_ij = du_i/dx_j`, let `S = sym(G)` and
//! `D = dev(sym(G G))`. The WALE rate is
//! `W = (D:D)^(3/2) / [(S:S)^(5/2) + (D:D)^(5/4)]` and
//! `nu_t = C_w^2 Delta^2 W`. Only the momentum basis changes from
//! `rho |S| S_dev` to `rho W S_dev`; scalar transport retains its independent
//! strain-based Germano identity and its original arithmetic.
//!
//! An unprotected dynamic WALE fit is ill-conditioned near a no-slip wall.
//! Baya Toda, Truffin & Nicoud, *Is the Dynamic Procedure Appropriate for All
//! SGS Models?*, ECCOMAS CFD 2010, §§5.1--5.2, equations (24)--(26), instead
//! use the shear/vortex sensor `SVS = (D:D)^(3/2) / [(D:D)^(3/2) + (S:S)^3]`.
//! Above its threshold the Germano fit supplies `C_w^2`; at or below it the
//! authored fixed coefficient supplies `C_w^2`. Their threshold is 0.09 and
//! fixed `C_w` is 0.5. Both are explicit here; there is no unprotected default.
//! <https://imag.umontpellier.fr/~nicoud/PDF/ECCOMAS_BayaToda.pdf>
//!
//! The same sensor applies to periodic flow and to the existing tangentially
//! filtered channel. It does not add a wall function, shock model, dynamic
//! scalar WALE model, or claim the paper's channel calibration for this grid.

use super::*;
use crate::sgs3::{SgsModel, eddy_viscosity};

/// Explicit parameters for the published sensor-protected dynamic WALE model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DynamicWaleConfig {
    /// Dimensionless SVS threshold, strictly between zero and one.
    pub sensor_threshold: f64,
    /// Fixed `C_w`, not its square, used at or below the sensor threshold.
    /// Zero disables only this shear-dominated branch, not the dynamic fit.
    pub near_wall_coefficient: f64,
}

impl DynamicWaleConfig {
    pub fn validate(self) -> Result<(), DynamicSgsError> {
        let squared = self.near_wall_coefficient * self.near_wall_coefficient;
        if !self.sensor_threshold.is_finite()
            || self.sensor_threshold <= 0.0
            || self.sensor_threshold >= 1.0
            || !self.near_wall_coefficient.is_finite()
            || self.near_wall_coefficient < 0.0
            || !squared.is_finite()
            || (squared == 0.0 && self.near_wall_coefficient != 0.0)
        {
            return Err(DynamicSgsError::Invalid(
                "WALE requires 0 < sensor threshold < 1 and representable nonnegative C_w squared",
            ));
        }
        Ok(())
    }

    pub(crate) fn validate_bound(self, maximum: f64) -> Result<(), DynamicSgsError> {
        self.validate()?;
        if self.near_wall_coefficient * self.near_wall_coefficient > maximum {
            return Err(DynamicSgsError::Invalid(
                "fixed WALE coefficient squared exceeds the momentum admission bound",
            ));
        }
        Ok(())
    }
}

pub(super) struct WaleWorkspace {
    pub(super) config: DynamicWaleConfig,
    pub(super) fixed_squared: f64,
    // [W(Delta G / velocity_scale), SVS]. Neither replaces |S| in the heat fit.
    pub(super) response: Vec<[f64; 2]>,
}

impl PeriodicDynamicSgs {
    /// Nominal numerical payload for dynamic WALE: 152 bytes per cell.
    /// The additional two scalars retain its rate and dimensionless sensor.
    pub fn required_wale_bytes(cells: usize) -> Result<usize, DynamicSgsError> {
        cells
            .checked_mul(19 * core::mem::size_of::<f64>())
            .ok_or(DynamicSgsError::Invalid("WALE workspace size overflow"))
    }

    /// Select sensor-protected dynamic WALE on this filter geometry.
    /// The byte limit covers all existing and added numerical capacities.
    /// Selection consumes the workspace; prepare it again before using outputs.
    /// The original Smagorinsky constructor and its byte budget are unchanged.
    pub fn with_wale(
        mut self,
        config: DynamicWaleConfig,
        maximum_bytes: usize,
    ) -> Result<Self, DynamicSgsError> {
        config.validate()?;
        if self.wale.is_some() {
            return Err(DynamicSgsError::Invalid("dynamic WALE is already selected"));
        }
        let extra = self
            .viscosity
            .len()
            .checked_mul(core::mem::size_of::<[f64; 2]>())
            .ok_or(DynamicSgsError::Invalid("WALE workspace size overflow"))?;
        let required = self
            .workspace_bytes()
            .checked_add(extra)
            .ok_or(DynamicSgsError::Invalid("WALE workspace size overflow"))?;
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            });
        }
        self.wale = Some(WaleWorkspace {
            config,
            fixed_squared: config.near_wall_coefficient * config.near_wall_coefficient,
            response: allocated(self.viscosity.len(), [0.0; 2])?,
        });
        let required = self.workspace_bytes();
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            });
        }
        Ok(self)
    }

    /// None identifies the original strain-based dynamic Smagorinsky closure.
    pub fn wale_config(&self) -> Option<DynamicWaleConfig> {
        self.wale.as_ref().map(|wale| wale.config)
    }
}

// Reuse the original WALE kernel with C_w = Delta = 1: its returned viscosity
// is exactly the rate W in the gradient's units. No second primal formula.
pub(super) fn rate(gradient: [[f64; 3]; 3], cell: usize) -> Result<f64, DynamicSgsError> {
    eddy_viscosity(SgsModel::Wale { cw: 1.0 }, 1.0, gradient).map_err(|_| DynamicSgsError::Cell {
        cell,
        reason: "unrepresentable normalized WALE rate",
    })
}

// Normalize before the sixth-degree sensor contractions. The common scale
// cancels exactly; no dimensionful floor or wall distance enters the model.
pub(super) fn sensor(gradient: [[f64; 3]; 3]) -> f64 {
    let scale = gradient
        .iter()
        .flatten()
        .fold(0.0_f64, |a, b| a.max(b.abs()));
    if scale == 0.0 {
        return 0.0;
    }
    let g = gradient.map(|row| row.map(|value| value / scale));
    let square: [[f64; 3]; 3] =
        core::array::from_fn(|i| core::array::from_fn(|j| (0..3).map(|k| g[i][k] * g[k][j]).sum()));
    let trace = (square[0][0] + square[1][1] + square[2][2]) / 3.0;
    let mut strain2 = 0.0;
    let mut deviator2 = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            let s = 0.5 * (g[i][j] + g[j][i]);
            let d = 0.5 * (square[i][j] + square[j][i]) - if i == j { trace } else { 0.0 };
            strain2 += s * s;
            deviator2 += d * d;
        }
    }
    let numerator = deviator2 * det::sqrt(deviator2);
    let denominator = numerator + strain2 * strain2 * strain2;
    if denominator == 0.0 {
        0.0
    } else {
        numerator / denominator
    }
}

#[cfg(test)]
#[path = "wale/tests.rs"]
mod tests;
