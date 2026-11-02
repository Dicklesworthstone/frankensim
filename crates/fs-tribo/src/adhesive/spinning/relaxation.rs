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
use crate::sliding_patch::contact::validate_circular_response;
use fs_math::det;

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
    #[allow(clippy::too_many_lines)]
    pub fn relax_velocity(
        &self,
        input: &AdhesiveRelaxationInput,
        config: AdhesiveWrenchConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdhesiveRelaxationReport, AdhesiveRelaxationError> {
        poll(&mut cancelled)?;
        config.validate()?;
        validate_circular_response(input.scaled_mobility_per_kg, config)?;
        if !input.dt_s.is_finite()
            || input.dt_s <= 0.0
            || input.velocity_m_s.into_iter().any(|x| !x.is_finite())
            || !input.spin_rate_rad_s.is_finite()
            || input.scaled_mobility_per_kg == [[0.0; 3]; 3]
        {
            return Err(AdhesiveRelaxationError::Invalid(
                "finite rates, positive duration and finite-body mobility required",
            ));
        }
        let radius = self.parameters().contact_radius_m;
        let initial = [
            input.velocity_m_s[0],
            input.velocity_m_s[1],
            finite(radius * input.spin_rate_rad_s)?,
        ];
        if input.spin_rate_rad_s != 0.0 && initial[2] == 0.0 {
            return Err(AdhesiveRelaxationError::Invalid(
                "unrepresentable nonzero radius-scaled spin rate",
            ));
        }
        let native = self
            .kinetic_patch()
            .evaluate(input.velocity_m_s, input.spin_rate_rad_s)
            .map_err(AdhesiveWrenchError::Law)?;
        if det::hypot(native.force_n[0], native.force_n[1]) > config.maximum_force_n
            || native.spin_moment_nm.abs() > config.maximum_spin_moment_nm
        {
            return Err(AdhesiveRelaxationError::EffortBudget);
        }
        let wrench = [
            native.force_n[0],
            native.force_n[1],
            finite(native.spin_moment_nm / radius)?,
        ];
        let acceleration = multiply(input.scaled_mobility_per_kg, wrench);
        for value in acceleration {
            finite(value)?;
        }
        poll(&mut cancelled)?;
        let at_rest = initial == [0.0; 3];
        let frictionless = wrench == [0.0; 3];
        let (moving, stop, fraction, ray_residual) = if at_rest {
            (0.0, Some(0.0), 1.0, 0.0)
        } else if frictionless {
            (input.dt_s, None, 0.0, 0.0)
        } else {
            // Pivot on the largest physical rate, not a small component. Test
            // direction before constructing time, so real cross-coupled paths
            // cannot be accepted merely because the requested interval is tiny.
            let pivot = (1..3).fold(0, |p, i| {
                if initial[i].abs() > initial[p].abs() {
                    i
                } else {
                    p
                }
            });
            let scale = acceleration.into_iter().map(f64::abs).fold(0.0, f64::max);
            if scale == 0.0 {
                return Err(AdhesiveRelaxationError::Invalid(
                    "positive kinetic resistance needs a representable finite stop",
                ));
            }
            let slope = -acceleration[pivot] / scale;
            let residual = (0..3)
                .map(|i| (acceleration[i] / scale + slope * (initial[i] / initial[pivot])).abs())
                .fold(0.0, f64::max);
            if residual > 256.0 * f64::EPSILON {
                return Err(AdhesiveRelaxationError::CurvedTrajectory {
                    ray_residual_relative: residual,
                });
            }
            let stop_time = finite(-initial[pivot] / acceleration[pivot])?;
            if stop_time <= 0.0 {
                return Err(AdhesiveRelaxationError::Invalid(
                    "original mobility must give a positive physical stopping time",
                ));
            }
            let moving = input.dt_s.min(stop_time);
            (
                moving,
                (stop_time <= input.dt_s).then_some(stop_time),
                moving / stop_time,
                residual,
            )
        };
        let remainder = if stop.is_some() {
            input.dt_s - moving
        } else {
            0.0
        };
        let factor = moving * (1.0 - 0.5 * fraction);
        let endpoint = initial.map(|value| value * (1.0 - fraction));
        let impulse = wrench.map(|value| value * moving);
        let change = multiply(input.scaled_mobility_per_kg, impulse);
        let residual = core::array::from_fn(|i| initial[i] + change[i] - endpoint[i]);
        let displacement_residual = residual.map(|value| 0.5 * moving * value);
        let absorbed = -(0..3)
            .map(|i| impulse[i] * (initial[i] + 0.5 * change[i]))
            .sum::<f64>();
        let heat = finite(native.dissipated_power_w * factor)?;
        let work_residual = finite(heat - absorbed)?;
        if displacement_residual
            .into_iter()
            .any(|v| !v.is_finite() || v.abs() > config.displacement_tolerance_m)
            || work_residual.abs() > config.energy_tolerance_j
            || !absorbed.is_finite()
            || absorbed < -config.energy_tolerance_j
            || heat < 0.0
        {
            return Err(AdhesiveRelaxationError::Unresolved {
                displacement_residual_m: displacement_residual,
                work_residual_j: work_residual,
            });
        }
        let result = AdhesiveRelaxationReport {
            initial_wrench: native,
            velocity_after_m_s: [endpoint[0], endpoint[1]],
            spin_rate_after_rad_s: finite(input.spin_rate_rad_s * (1.0 - fraction))?,
            force_impulse_ns: [finite(impulse[0])?, finite(impulse[1])?],
            spin_angular_impulse_nms: finite(native.spin_moment_nm * moving)?,
            slip_increment_m: [
                finite(input.velocity_m_s[0] * factor)?,
                finite(input.velocity_m_s[1] * factor)?,
            ],
            spin_increment_rad: finite(input.spin_rate_rad_s * factor)?,
            dissipated_energy_j: heat,
            contact_work_j: absorbed,
            work_residual_j: work_residual,
            moving_duration_s: moving,
            kinetic_phase_end_s: stop,
            stationary_duration_s: remainder,
            scaled_velocity_residual_m_s: residual,
            ray_residual_relative: ray_residual,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }
}

fn multiply(matrix: [[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    matrix.map(|row| row[0] * vector[0] + row[1] * vector[1] + row[2] * vector[2])
}
fn finite(value: f64) -> Result<f64, AdhesiveRelaxationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(AdhesiveRelaxationError::Invalid(
            "unrepresentable adhesive trajectory",
        ))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AdhesiveRelaxationError> {
    if cancelled() {
        Err(AdhesiveRelaxationError::Cancelled)
    } else {
        Ok(())
    }
}
