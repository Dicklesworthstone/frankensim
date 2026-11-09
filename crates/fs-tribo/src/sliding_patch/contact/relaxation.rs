//! Continuous Hertz sliding/spinning through an actual finite stopping event.
//!
//! In radius-scaled coordinates `q = [u_x,u_y,a*omega]`, the original Hertz
//! wrench is homogeneous of degree zero. If the complete native mobility maps
//! that wrench to `-q/t_stop`, the direction stays fixed until exact arrest.
//! The stationary remainder has zero impulse and zero heat. Pure translation,
//! pure spin and compatible combined motion use the same full-matrix path.
//!
//! This shares the invariant-direction integrator with uniform adhesive disks;
//! each supplies its own original pressure/traction integral exactly once.
//! Curved velocity paths need another trajectory owner. Geometry, normal load,
//! material coefficients and inertia remain held, with no pose integration,
//! external drive, velocity capture band or independent force/moment clipping.

use super::{HertzWrenchConfig, HertzWrenchContact, HertzWrenchError, validate_circular_response};
use crate::sliding_patch::HertzSlidingWrench;
use fs_math::det;

/// Actual accepted rates and complete finite-body response, held for the phase.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzRelaxationInput {
    pub velocity_m_s: [f64; 2],
    pub spin_rate_rad_s: f64,
    /// `S*J*M^-1*J^T*S`, `S = diag(1,1,a)`, in [1/kg].
    pub scaled_mobility_per_kg: [[f64; 3]; 3],
    /// Requested duration, including any stationary remainder [s].
    pub dt_s: f64,
}

/// Integrated original Hertz traction and actual continuous travel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HertzRelaxationReport {
    pub initial_wrench: HertzSlidingWrench,
    /// Analytic event endpoint; native bodies independently check their recoil.
    pub velocity_after_m_s: [f64; 2],
    pub spin_rate_after_rad_s: f64,
    pub force_impulse_ns: [f64; 2],
    pub spin_angular_impulse_nms: f64,
    pub slip_increment_m: [f64; 2],
    pub spin_increment_rad: f64,
    /// Physical pressure-resolved slip loss, without a heat partition [J].
    pub dissipated_energy_j: f64,
    /// Absorbed impulse work against the original complete mobility [J].
    pub contact_work_j: f64,
    /// Physical loss minus impulse work; never added to heat [J].
    pub work_residual_j: f64,
    pub moving_duration_s: f64,
    /// Located stop in the interval; `Some(0)` denotes initial rest.
    pub kinetic_phase_end_s: Option<f64>,
    pub stationary_duration_s: f64,
    /// Original impulse endpoint minus the analytic event endpoint [m/s].
    pub scaled_velocity_residual_m_s: [f64; 3],
    /// Directional mismatch admitted only at arithmetic roundoff.
    pub ray_residual_relative: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum HertzRelaxationError {
    Contact(HertzWrenchError),
    Invalid(&'static str),
    CurvedTrajectory {
        ray_residual_relative: f64,
    },
    EffortBudget,
    Unresolved {
        displacement_residual_m: [f64; 3],
        work_residual_j: f64,
    },
    Cancelled,
}
impl From<HertzWrenchError> for HertzRelaxationError {
    fn from(error: HertzWrenchError) -> Self {
        Self::Contact(error)
    }
}
impl core::fmt::Display for HertzRelaxationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Hertz fixed-direction relaxation: {self:?}")
    }
}
impl std::error::Error for HertzRelaxationError {}

impl HertzWrenchContact {
    /// Integrate original kinetic traction up to arrest, then stationary dwell.
    ///
    /// The static coefficient never captures a moving phase early. One native
    /// Hertz evaluation and fixed-size algebra resolve the entire duration;
    /// no time stepping or root search is needed for an admitted direction.
    pub fn relax_velocity(
        &self,
        input: &HertzRelaxationInput,
        config: HertzWrenchConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HertzRelaxationReport, HertzRelaxationError> {
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
                    .kinetic_patch
                    .evaluate(input.velocity_m_s, input.spin_rate_rad_s)
                    .map_err(HertzWrenchError::Law)?;
                Ok(CircularTraction {
                    original: native,
                    force_n: native.force_n,
                    spin_moment_nm: native.spin_moment_nm,
                    dissipated_power_w: native.dissipated_power_w,
                })
            },
            cancelled,
        )?;
        Ok(HertzRelaxationReport {
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

// Shared mechanics contain no Hertz or adhesive constitutive approximation.
// Public wrappers retain their original law-specific types and diagnostics.
pub(crate) struct CircularRelaxationInput {
    pub velocity_m_s: [f64; 2],
    pub spin_rate_rad_s: f64,
    pub scaled_mobility_per_kg: [[f64; 3]; 3],
    pub dt_s: f64,
}
pub(crate) struct CircularTraction<W> {
    pub original: W,
    pub force_n: [f64; 2],
    pub spin_moment_nm: f64,
    pub dissipated_power_w: f64,
}
pub(crate) struct CircularRelaxationReport<W> {
    pub initial_wrench: W,
    pub velocity_after_m_s: [f64; 2],
    pub spin_rate_after_rad_s: f64,
    pub force_impulse_ns: [f64; 2],
    pub spin_angular_impulse_nms: f64,
    pub slip_increment_m: [f64; 2],
    pub spin_increment_rad: f64,
    pub dissipated_energy_j: f64,
    pub contact_work_j: f64,
    pub work_residual_j: f64,
    pub moving_duration_s: f64,
    pub kinetic_phase_end_s: Option<f64>,
    pub stationary_duration_s: f64,
    pub scaled_velocity_residual_m_s: [f64; 3],
    pub ray_residual_relative: f64,
}
pub(crate) enum CircularRelaxationError {
    Contact(HertzWrenchError),
    Invalid(&'static str),
    NonFinite,
    CurvedTrajectory {
        ray_residual_relative: f64,
    },
    EffortBudget,
    Unresolved {
        displacement_residual_m: [f64; 3],
        work_residual_j: f64,
    },
    Cancelled,
}
impl From<HertzWrenchError> for CircularRelaxationError {
    fn from(error: HertzWrenchError) -> Self {
        Self::Contact(error)
    }
}
impl From<CircularRelaxationError> for HertzRelaxationError {
    fn from(error: CircularRelaxationError) -> Self {
        match error {
            CircularRelaxationError::Contact(error) => Self::Contact(error),
            CircularRelaxationError::Invalid(field) => Self::Invalid(field),
            CircularRelaxationError::NonFinite => Self::Invalid("unrepresentable Hertz trajectory"),
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

/// Exact degree-zero-traction flow when the full response preserves its ray.
#[allow(clippy::too_many_lines)]
pub(crate) fn relax_circular<W>(
    input: CircularRelaxationInput,
    radius: f64,
    config: HertzWrenchConfig,
    evaluate: impl FnOnce() -> Result<CircularTraction<W>, HertzWrenchError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<CircularRelaxationReport<W>, CircularRelaxationError> {
    poll(&mut cancelled)?;
    config.validate()?;
    validate_circular_response(input.scaled_mobility_per_kg, config)?;
    if !input.dt_s.is_finite()
        || input.dt_s <= 0.0
        || input.velocity_m_s.into_iter().any(|x| !x.is_finite())
        || !input.spin_rate_rad_s.is_finite()
        || input.scaled_mobility_per_kg == [[0.0; 3]; 3]
    {
        return Err(CircularRelaxationError::Invalid(
            "finite rates, positive duration and finite-body mobility required",
        ));
    }
    let initial = [
        input.velocity_m_s[0],
        input.velocity_m_s[1],
        finite(radius * input.spin_rate_rad_s)?,
    ];
    if input.spin_rate_rad_s != 0.0 && initial[2] == 0.0 {
        return Err(CircularRelaxationError::Invalid(
            "unrepresentable nonzero radius-scaled spin rate",
        ));
    }
    let native = evaluate()?;
    if det::hypot(native.force_n[0], native.force_n[1]) > config.maximum_force_n
        || native.spin_moment_nm.abs() > config.maximum_spin_moment_nm
    {
        return Err(CircularRelaxationError::EffortBudget);
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
        // Pivot on the largest physical rate. A tiny duration cannot hide
        // genuine transverse acceleration in the original full mobility.
        let pivot = (1..3).fold(0, |p, i| {
            if initial[i].abs() > initial[p].abs() {
                i
            } else {
                p
            }
        });
        let scale = acceleration.into_iter().map(f64::abs).fold(0.0, f64::max);
        if scale == 0.0 {
            return Err(CircularRelaxationError::Invalid(
                "positive kinetic resistance needs a representable finite stop",
            ));
        }
        let slope = -acceleration[pivot] / scale;
        let residual = (0..3)
            .map(|i| (acceleration[i] / scale + slope * (initial[i] / initial[pivot])).abs())
            .fold(0.0, f64::max);
        if residual > 256.0 * f64::EPSILON {
            return Err(CircularRelaxationError::CurvedTrajectory {
                ray_residual_relative: residual,
            });
        }
        let stop_time = finite(-initial[pivot] / acceleration[pivot])?;
        if stop_time <= 0.0 {
            return Err(CircularRelaxationError::Invalid(
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
        return Err(CircularRelaxationError::Unresolved {
            displacement_residual_m: displacement_residual,
            work_residual_j: work_residual,
        });
    }
    let result = CircularRelaxationReport {
        initial_wrench: native.original,
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

fn multiply(matrix: [[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    matrix.map(|row| row[0] * vector[0] + row[1] * vector[1] + row[2] * vector[2])
}
fn finite(value: f64) -> Result<f64, CircularRelaxationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CircularRelaxationError::NonFinite)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CircularRelaxationError> {
    if cancelled() {
        Err(CircularRelaxationError::Cancelled)
    } else {
        Ok(())
    }
}
