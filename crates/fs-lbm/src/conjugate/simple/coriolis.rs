//! Constant local-frame rotation on the original staggered momentum owner.
//!
//! With P the average of each cell's two normal faces and V its physical
//! volume, the left-hand-side action is C = P^T V [2 rho Omega]x P.
//! C is skew even on graded cells: u^T C u = 0. Closed normal boundaries
//! eliminate fixed-face work; moving tangential walls retain their ordinary
//! viscous actuator account. Open momentum rows mirror boundary cells and
//! need a distinct adjoint boundary treatment, so are not admitted here.
//!
//! This source is evaluated at the current coupled velocity, including the
//! native BE/BDF2 endpoint. It never enters material force refresh, Darcy
//! resistance or physical dissipation. The final frozen momentum residual
//! includes this action after pressure projection.

use super::*;

impl SimpleConfig {
    pub(super) fn admit_coriolis(&self, density_kg_m3: f64) -> Result<(), ChtError> {
        if self.coriolis_angular_velocity_rad_s == [0.0; 3] {
            return Ok(());
        }
        for omega in self.coriolis_angular_velocity_rad_s {
            finite("simple.coriolis_angular_velocity_rad_s", omega)?;
            finite("simple.coriolis_mass_rate", 2.0 * (density_kg_m3 * omega))?;
        }
        if self
            .faces
            .iter()
            .any(|face| !matches!(face, FvBoundary::Wall { .. } | FvBoundary::Symmetry))
        {
            return Err(ChtError::InvalidInput {
                field: "simple.coriolis_angular_velocity_rad_s",
                reason: "local Coriolis requires closed normal wall/symmetry boundaries".into(),
            });
        }
        Ok(())
    }
}

impl Solver<'_> {
    pub(in crate::conjugate) fn has_coriolis(&self) -> bool {
        self.config.coriolis_angular_velocity_rad_s != [0.0; 3]
    }

    /// The separate actual-mass viscosity substep creates a zero-spin config.
    /// Refuse any future combination until rotation shares that inertia metric.
    pub(super) fn admit_coriolis_inertia(&self) -> Result<(), ChtError> {
        if self.has_coriolis()
            && self
                .time
                .as_ref()
                .is_some_and(|time| time.face_mass_kg.is_some())
        {
            return Err(ChtError::InvalidInput {
                field: "simple.coriolis_inertia",
                reason: "local Coriolis requires the native constant-density face inertia".into(),
            });
        }
        Ok(())
    }

    /// Integrated positive left-hand-side force [N]; momentum subtracts it
    /// from the right-hand side to implement -2*rho*Omega cross u.
    pub(super) fn coriolis_face_force(
        &self,
        axis: usize,
        minus: [usize; 3],
        plus: [usize; 3],
    ) -> f64 {
        let j = (axis + 1) % 3;
        let k = (axis + 2) % 3;
        [minus, plus]
            .into_iter()
            .map(|at| {
                let c = cell_of(self.n, at);
                let u = self.cell_velocity(c);
                let cross = self.coriolis_mass_rate[j] * u[k]
                    - self.coriolis_mass_rate[k] * u[j];
                0.5 * self.domain.volume(c) * cross
            })
            .sum()
    }
}
