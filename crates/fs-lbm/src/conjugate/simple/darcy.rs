//! Symmetric cross-coordinate Darcy action on the original MAC velocities.
//!
//! The diagonal uses the existing half-cell quadrature without alteration.
//! With P the average of each cell's two normal faces, the cross term is
//! P^T V (B-diag(B)) P. Its transpose distribution is the native drag
//! action, and u^T times that SAME action supplies the original-cell heat.
//! The complete cell quadratic is
//! V [ubar^T B ubar + sum_a B_aa (u_a+ - u_a-)^2 / 4] >= 0 for SPD B.
//! This retains diagonal face coercivity and avoids a checkerboard null mode.

use super::*;

/// Symmetric xy, xz, yz entries, compiled once per material source epoch.
#[derive(Clone, Copy, Default)]
pub(super) struct DarcyCross([f64; 3]);

impl DarcyCross {
    pub(super) fn from_tensor(matrix: &[[f64; 3]; 3]) -> Result<Self, ChtError> {
        let refuse = || {
            ChtError::InvalidInput {
            field: "simple.darcy_resistance_pa_s_m2",
            reason: "Darcy resistance must be finite symmetric positive definite, or a nonnegative diagonal".into(),
        }
        };
        for a in 0..3 {
            for b in 0..3 {
                if !matrix[a][b].is_finite() || matrix[a][b] != matrix[b][a] {
                    return Err(refuse());
                }
            }
            if matrix[a][a] < 0.0 {
                return Err(refuse());
            }
        }
        let cross = Self([matrix[0][1], matrix[0][2], matrix[1][2]]);
        if cross.is_zero() {
            return Ok(cross);
        }
        // Scale before Cholesky: the criterion depends on conditioning, not
        // whether a permeability was stated in a large or small SI magnitude.
        let scale = (0..3).map(|a| matrix[a][a]).fold(0.0_f64, f64::max);
        if scale <= 0.0 {
            return Err(refuse());
        }
        let a = matrix[0][0] / scale;
        if a <= 0.0 {
            return Err(refuse());
        }
        let l00 = fs_math::det::sqrt(a);
        let l10 = (matrix[1][0] / scale) / l00;
        let l20 = (matrix[2][0] / scale) / l00;
        let d1 = matrix[1][1] / scale - l10 * l10;
        if !d1.is_finite() || d1 <= 0.0 {
            return Err(refuse());
        }
        let l21 = (matrix[2][1] / scale - l20 * l10) / fs_math::det::sqrt(d1);
        let d2 = matrix[2][2] / scale - l20 * l20 - l21 * l21;
        if !d2.is_finite() || d2 <= 0.0 {
            return Err(refuse());
        }
        Ok(cross)
    }

    pub(super) fn is_zero(self) -> bool {
        self.0 == [0.0; 3]
    }

    fn action(self, u: [f64; 3]) -> [f64; 3] {
        let [xy, xz, yz] = self.0;
        [
            xy * u[1] + xz * u[2],
            xy * u[0] + yz * u[2],
            xz * u[0] + yz * u[1],
        ]
    }

    pub(super) fn rescale_viscosity(&mut self, old: f64, new: f64) -> Result<(), ChtError> {
        for value in &mut self.0 {
            *value = (*value / old) * new;
            finite("simple.material_darcy_resistance", *value)?;
        }
        Ok(())
    }
}

impl Solver<'_> {
    /// Half of each physical cell contributes to this normal face. A coupled
    /// porous cell cannot touch a mirrored pressure/inlet momentum boundary;
    /// those need a distinct open-port work reduction, not altered weights.
    pub(super) fn darcy_cross_force(
        &self,
        axis: usize,
        minus: [usize; 3],
        plus: [usize; 3],
    ) -> f64 {
        if self.porous_cross.is_empty() {
            return 0.0;
        }
        [minus, plus]
            .into_iter()
            .map(|at| {
                let c = cell_of(self.n, at);
                0.5 * self.domain.volume(c)
                    * self.porous_cross[c].action(self.cell_velocity(c))[axis]
            })
            .sum()
    }

    /// Signed cross contribution; add to the original nonnegative diagonal
    /// half-face quadrature, then admit the complete cell loss. The cross term
    /// alone is not a dissipation channel and may legitimately be negative.
    pub(super) fn darcy_cross_power(&self, cell: usize) -> f64 {
        if self.porous_cross.is_empty() {
            return 0.0;
        }
        let u = self.cell_velocity(cell);
        let force = self.porous_cross[cell].action(u);
        self.domain.volume(cell) * (u[0] * force[0] + u[1] * force[1] + u[2] * force[2])
    }

    pub(super) fn admit_darcy_opening(
        &self,
        cell: usize,
        cross: DarcyCross,
    ) -> Result<(), ChtError> {
        if !cross.is_zero()
            && (0..6).any(|side| {
                self.domain.neighbor(cell, side).is_none()
                    && matches!(
                        self.faces[side],
                        FvBoundary::Inlet { .. } | FvBoundary::Outlet | FvBoundary::Pressure { .. }
                    )
            })
        {
            return Err(ChtError::InvalidInput {
                field: "simple.darcy_resistance_pa_s_m2",
                reason: "rotated Darcy cells need closed normal boundaries or a clear-fluid layer before pressure/inlet openings".into(),
            });
        }
        Ok(())
    }

    /// Reassemble all three residuals on ONE frozen current velocity/pressure
    /// iterate. Sequential component solves alone can hide stale cross forces
    /// or a nonlinear Forchheimer coefficient evaluated at an earlier speed.
    /// This probe never changes velocities, pressure or correction diagonals.
    pub(super) fn coupled_darcy_residual(&mut self, gate: &CancelGate) -> Result<f64, ChtError> {
        let mut rms = [0.0_f64; 3];
        for (axis, value) in rms.iter_mut().enumerate() {
            *value = self.momentum(axis, false, gate)?.rms;
        }
        let scale = self.momentum_rms.iter().copied().fold(0.0_f64, f64::max);
        let defect = rms.into_iter().fold(0.0_f64, f64::max);
        let residual = if scale > 0.0 { defect / scale } else { defect };
        finite("simple.coupled_darcy_residual", residual)?;
        Ok(residual)
    }

    /// The final tight pressure projection must not invalidate momentum.
    pub(super) fn admit_projected_darcy(&mut self, gate: &CancelGate) -> Result<f64, ChtError> {
        if self.porous_cross.is_empty() && !self.porous_nonlinear && self.pressure_work_axis.is_none() {
            return Ok(0.0);
        }
        let residual = self.coupled_darcy_residual(gate)?;
        if residual > self.config.tolerance {
            return Err(ChtError::FlowNotSteady {
                steps: self.sweeps,
                last_change: residual,
                tolerance: self.config.tolerance,
            });
        }
        Ok(residual)
    }
}
