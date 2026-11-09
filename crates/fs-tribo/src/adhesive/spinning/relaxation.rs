//! Continuous uniform-disk relaxation on an invariant slip/spin ray.
//!
//! In conjugate coordinates `q = [u_x,u_y,a*omega]`, the real disk supplies
//! `w(q) = [F_x,F_y,M/a]`. Its gross-sliding wrench is homogeneous of degree
//! zero. If the **complete** original mobility satisfies `W*w(q) = -q/t_stop`,
//! that same wrench acts until the exact finite event `t_stop`, with
//! `q(t) = (1-t/t_stop)*q(0)`. The remaining no-drive interval is stationary:
//! zero reaction belongs to the original static set and adds no impulse or heat.
//!
//! This includes pure translation and spin when their recoil stays in that
//! subspace, and combined motion for appropriate full body mobilities. A curved
//! velocity path is refused at machine-roundoff ray admission; it is never
//! replaced by an ellipse, an independent force/moment clamp, or clipped rates.
//! Geometry, real contact area, strengths and body response are held fixed.
//! Pose/normal evolution and changing sources belong to the surrounding owner.

use super::AdhesiveSlidingWrench;
use super::contact::{AdhesiveWrenchConfig, AdhesiveWrenchContact, AdhesiveWrenchError};
use crate::sliding_patch::contact::relaxation::{
    CircularRelaxationError, CircularRelaxationInput, CircularTraction, relax_circular,
};

/// Actual accepted rates and the original complete finite-body response.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveRelaxationInput {
    /// Initial relative velocity at the disk centre [m/s].
    pub velocity_m_s: [f64; 2],
    /// Initial relative rotation rate about the disk normal [rad/s].
    pub spin_rate_rad_s: f64,
    /// Original `S*J*M^-1*J^T*S`, `S = diag(1,1,a)`, in [1/kg].
    /// Positive definite; no artificial mass, regularization or diagonalization.
    pub scaled_mobility_per_kg: [[f64; 3]; 3],
    /// Complete requested physical duration, including any stationary remainder.
    pub dt_s: f64,
}

/// Integrated original traction along the admitted continuous trajectory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveRelaxationReport {
    /// Original kinetic disk wrench at the initial physical motion.
    pub initial_wrench: AdhesiveSlidingWrench,
    /// Analytic endpoint velocity; a located event gives exactly zero.
    /// Native bodies must independently check this target after applying impulses.
    pub velocity_after_m_s: [f64; 2],
    /// Analytic endpoint relative spin, zero only after an admitted event.
    pub spin_rate_after_rad_s: f64,
    /// Acting force integrated only over the moving phase [N s].
    pub force_impulse_ns: [f64; 2],
    /// Acting normal-axis moment integrated over the moving phase [N m s].
    pub spin_angular_impulse_nms: f64,
    /// Integrated actual centre slip, excluding stationary dwell [m].
    pub slip_increment_m: [f64; 2],
    /// Integrated relative spin, excluding stationary dwell [rad].
    pub spin_increment_rad: f64,
    /// Original local junction loss integrated over the moving phase [J].
    pub dissipated_energy_j: f64,
    /// Positive absorbed work from the original full mobility and total impulse.
    pub contact_work_j: f64,
    /// Loss minus actual absorbed impulse work; never added to heat [J].
    pub work_residual_j: f64,
    /// Time occupied by the moving branch, or zero for initially stationary data.
    pub moving_duration_s: f64,
    /// Located stop within this interval, including `Some(0)` for initial rest.
    pub kinetic_phase_end_s: Option<f64>,
    /// Duration held at the located stop with zero external drive [s].
    pub stationary_duration_s: f64,
    /// Original-matrix endpoint minus the analytic ray endpoint [m/s].
    /// The third coordinate is radius times the spin-rate discrepancy.
    pub scaled_velocity_residual_m_s: [f64; 3],
    /// Relative acceleration misalignment admitted only at arithmetic roundoff.
    pub ray_residual_relative: f64,
}

/// Invalid mechanics, unsupported curved motion, and bounded refusals stay distinct.
#[derive(Clone, Debug, PartialEq)]
pub enum AdhesiveRelaxationError {
    /// Original material, numerical-control or matrix admission refusal.
    Contact(AdhesiveWrenchError),
    /// Nonfinite, nonpositive or unrepresentable physical input/result.
    Invalid(&'static str),
    /// The actual mobility rotates the velocity; another trajectory owner is needed.
    CurvedTrajectory {
        /// Maximum dimensionless acceleration departure from the initial ray.
        ray_residual_relative: f64,
    },
    /// The actual initial force or spin moment exceeds its declared budget.
    EffortBudget,
    /// Original work or displacement gates reject the roundoff-level candidate.
    Unresolved {
        /// Integrated original response minus the analytic path [m].
        displacement_residual_m: [f64; 3],
        /// Integrated native loss minus actual absorbed work [J].
        work_residual_j: f64,
    },
    /// No candidate is accepted after a cooperative cancellation request.
    Cancelled,
}
impl From<AdhesiveWrenchError> for AdhesiveRelaxationError {
    fn from(error: AdhesiveWrenchError) -> Self {
        Self::Contact(error)
    }
}
impl core::fmt::Display for AdhesiveRelaxationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adhesive fixed-ray relaxation: {self:?}")
    }
}
impl std::error::Error for AdhesiveRelaxationError {}

impl AdhesiveWrenchContact {
    /// Resolve one passive fixed-ray trajectory through its finite stopping event.
    ///
    /// Unlike a held-effort midpoint sample, time after the actual stop receives
    /// no kinetic wrench. The larger static strength never captures a moving
    /// junction early. The support evaluation and three-coordinate algebra are
    /// bounded and allocation-free; no time stepping or root iteration is needed.
    pub fn relax_velocity(
        &self,
        input: &AdhesiveRelaxationInput,
        config: AdhesiveWrenchConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AdhesiveRelaxationReport, AdhesiveRelaxationError> {
        let result = relax_circular(
            CircularRelaxationInput {
                velocity_m_s: input.velocity_m_s,
                spin_rate_rad_s: input.spin_rate_rad_s,
                scaled_mobility_per_kg: input.scaled_mobility_per_kg,
                dt_s: input.dt_s,
            },
            self.parameters().contact_radius_m,
            config,
            || {
                let native = self
                    .kinetic_patch()
                    .evaluate(input.velocity_m_s, input.spin_rate_rad_s)
                    .map_err(AdhesiveWrenchError::Law)?;
                Ok(CircularTraction {
                    original: native,
                    force_n: native.force_n,
                    spin_moment_nm: native.spin_moment_nm,
                    dissipated_power_w: native.dissipated_power_w,
                })
            },
            cancelled,
        )?;
        Ok(AdhesiveRelaxationReport {
            initial_wrench: result.initial_wrench,
            velocity_after_m_s: result.velocity_after_m_s,
            spin_rate_after_rad_s: result.spin_rate_after_rad_s,
            force_impulse_ns: result.force_impulse_ns,
            spin_angular_impulse_nms: result.spin_angular_impulse_nms,
            slip_increment_m: result.slip_increment_m,
            spin_increment_rad: result.spin_increment_rad,
            dissipated_energy_j: result.dissipated_energy_j,
            contact_work_j: result.contact_work_j,
            work_residual_j: result.work_residual_j,
            moving_duration_s: result.moving_duration_s,
            kinetic_phase_end_s: result.kinetic_phase_end_s,
            stationary_duration_s: result.stationary_duration_s,
            scaled_velocity_residual_m_s: result.scaled_velocity_residual_m_s,
            ray_residual_relative: result.ray_residual_relative,
        })
    }
}

impl From<CircularRelaxationError> for AdhesiveRelaxationError {
    fn from(error: CircularRelaxationError) -> Self {
        match error {
            CircularRelaxationError::Contact(error) => Self::Contact(error),
            CircularRelaxationError::Invalid(field) => Self::Invalid(field),
            CircularRelaxationError::NonFinite => {
                Self::Invalid("unrepresentable adhesive trajectory")
            }
            CircularRelaxationError::CurvedTrajectory {
                ray_residual_relative,
            } => Self::CurvedTrajectory {
                ray_residual_relative,
            },
            CircularRelaxationError::EffortBudget => Self::EffortBudget,
            CircularRelaxationError::Unresolved {
                displacement_residual_m,
                work_residual_j,
            } => Self::Unresolved {
                displacement_residual_m,
                work_residual_j,
            },
            CircularRelaxationError::Cancelled => Self::Cancelled,
        }
    }
}
