//! A stationary linear spring and viscous damper in the native pressure solve.

use super::{PistonControl, PistonError, PistonLoad, PistonStep, StepReport};

/// Passive axial suspension anchored in the gas mesh's world coordinate.
/// Stiffness and damping are explicitly supplied physical coefficients; neither
/// gas pressure nor acoustic impedance is inferred from them. Rest position is
/// the piston coordinate at zero spring force, not the initial gas length.
/// Coefficients and anchor remain fixed during each reported step. Changing
/// stiffness or rest position between calls changes stored energy outside that
/// step; the caller must account for that constitutive or anchor work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PistonSuspension {
    stiffness_n_m: f64,
    damping_n_s_m: f64,
    rest_position_m: f64,
}

impl PistonSuspension {
    /// Admit nonnegative finite stiffness [N/m], damping [N s/m] and a finite
    /// zero-force position [m]. Zero coefficients retain the free-piston limit.
    pub fn new(
        stiffness_n_m: f64,
        damping_n_s_m: f64,
        rest_position_m: f64,
    ) -> Result<Self, PistonError> {
        if !stiffness_n_m.is_finite()
            || stiffness_n_m < 0.0
            || !damping_n_s_m.is_finite()
            || damping_n_s_m < 0.0
            || !rest_position_m.is_finite()
        {
            return Err(PistonError::Invalid(
                "piston suspension requires nonnegative finite stiffness/damping and finite rest position",
            ));
        }
        Ok(Self {
            stiffness_n_m,
            damping_n_s_m,
            rest_position_m,
        })
    }

    /// Authored spring stiffness [N/m].
    #[must_use]
    pub const fn stiffness_n_m(self) -> f64 {
        self.stiffness_n_m
    }

    /// Authored viscous coefficient [N s/m].
    #[must_use]
    pub const fn damping_n_s_m(self) -> f64 {
        self.damping_n_s_m
    }

    /// Zero-force world position [m].
    #[must_use]
    pub const fn rest_position_m(self) -> f64 {
        self.rest_position_m
    }

    /// Spring storage at one actual mesh endpoint [J], with zero at rest.
    pub fn stored_energy_j(self, position_m: f64) -> Result<f64, PistonError> {
        Ok(self.bind(position_m)?.energy_before_j)
    }

    pub(crate) fn bind(self, position_m: f64) -> Result<BoundSuspension, PistonError> {
        if !position_m.is_finite() {
            return Err(PistonError::Invalid("nonfinite piston suspension position"));
        }
        let extension = if self.stiffness_n_m == 0.0 {
            0.0
        } else {
            position_m - self.rest_position_m
        };
        let energy_before_j = self.energy_at_extension(extension)?;
        Ok(BoundSuspension {
            law: self,
            extension,
            energy_before_j,
        })
    }

    fn energy_at_extension(self, extension: f64) -> Result<f64, PistonError> {
        let energy = (0.5 * (self.stiffness_n_m * extension)) * extension;
        if !extension.is_finite() || !energy.is_finite() || energy < 0.0 {
            return Err(PistonError::Invalid(
                "piston spring storage is unrepresentable",
            ));
        }
        Ok(energy)
    }
}

/// Accepted source gas/piston exchange with its distinct mechanical storage/loss.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SuspendedPistonStep<R = StepReport> {
    /// Original pressure-field and inertia report. Its additional-force work
    /// includes the spring and damper; `applied_force_work_j` below isolates
    /// the authored drive. Original gas/species/work observations are retained.
    pub piston: PistonStep<R>,
    /// Initial/final spring stored energy [J].
    pub spring_energy_j: [f64; 2],
    /// Exact quadratic-potential increment, evaluated before cancellation [J].
    pub spring_energy_change_j: f64,
    /// Nonnegative midpoint viscous loss `c*w²*dt` [J], not gas heating.
    pub damping_energy_j: f64,
    /// Spring-plus-damper impulse ON the piston [N s]; the anchor receives
    /// its opposite. The stationary anchor does no mechanical work.
    pub suspension_impulse_n_s: f64,
    /// Work from the separately held `PistonLoad::external_force_n` [J].
    pub applied_force_work_j: f64,
    /// Native mechanical balance including spring storage and viscous loss [J].
    pub mechanical_energy_residual_j: f64,
    /// Gas energy + piston kinetic + spring storage + damping loss + exterior
    /// pressure work - independently authored drive work [J].
    pub total_energy_residual_j: f64,
}

#[derive(Clone, Copy)]
pub(crate) struct BoundSuspension {
    law: PistonSuspension,
    extension: f64,
    energy_before_j: f64,
}

impl BoundSuspension {
    fn forces(self, speed: f64, displacement: f64) -> Result<(f64, f64), PistonError> {
        let spring = -self.law.stiffness_n_m * (self.extension + 0.5 * displacement);
        let damper = -self.law.damping_n_s_m * speed;
        if !spring.is_finite() || !damper.is_finite() {
            return Err(PistonError::Invalid(
                "piston suspension force is unrepresentable",
            ));
        }
        Ok((spring, damper))
    }

    pub(crate) fn load(
        self,
        drive: PistonLoad,
        speed: f64,
        displacement: f64,
    ) -> Result<PistonLoad, PistonError> {
        if self.law.stiffness_n_m == 0.0 && self.law.damping_n_s_m == 0.0 {
            return Ok(drive);
        }
        let (spring, damper) = self.forces(speed, displacement)?;
        let external_force_n = drive.external_force_n + spring + damper;
        if !external_force_n.is_finite() {
            return Err(PistonError::Invalid(
                "combined piston drive and suspension overflow",
            ));
        }
        Ok(PistonLoad {
            external_force_n,
            ..drive
        })
    }
}

pub(crate) fn finish<R>(
    suspension: Option<BoundSuspension>,
    piston: PistonStep<R>,
    drive: PistonLoad,
    duration_s: f64,
    controls: PistonControl,
) -> Result<SuspendedPistonStep<R>, PistonError> {
    let mut step = SuspendedPistonStep {
        spring_energy_j: [0.0; 2],
        spring_energy_change_j: 0.0,
        damping_energy_j: 0.0,
        suspension_impulse_n_s: 0.0,
        applied_force_work_j: piston.external_force_work_j,
        mechanical_energy_residual_j: piston.mechanical_energy_residual_j,
        total_energy_residual_j: piston.total_energy_residual_j,
        piston,
    };
    let Some(suspension) = suspension else {
        return Ok(step);
    };
    let dx = step.piston.displacement_m;
    let speed = step.piston.midpoint_velocity_m_s;
    let (spring, damper) = suspension.forces(speed, dx)?;
    step.spring_energy_j = [
        suspension.energy_before_j,
        suspension
            .law
            .energy_at_extension(suspension.extension + dx)?,
    ];
    // The original quadratic potential difference, retaining small increments
    // beside a large preload. Damping is independently integrated in time,
    // so represented dx != dt*w cannot silently manufacture work closure.
    step.spring_energy_change_j = -spring * dx;
    step.damping_energy_j = (-damper * speed) * duration_s;
    step.suspension_impulse_n_s = (spring + damper) * duration_s;
    step.applied_force_work_j = drive.external_force_n * dx;
    let correction = (step.piston.external_force_work_j - step.applied_force_work_j)
        + step.spring_energy_change_j
        + step.damping_energy_j;
    step.mechanical_energy_residual_j += correction;
    step.total_energy_residual_j += correction;
    if ![
        step.spring_energy_change_j,
        step.damping_energy_j,
        step.suspension_impulse_n_s,
        step.applied_force_work_j,
        step.mechanical_energy_residual_j,
        step.total_energy_residual_j,
    ]
    .iter()
    .all(|value| value.is_finite())
        || step.damping_energy_j < 0.0
    {
        return Err(PistonError::Invalid(
            "piston suspension impulse or work is unrepresentable",
        ));
    }
    let energy_residual_j = step
        .mechanical_energy_residual_j
        .abs()
        .max(step.total_energy_residual_j.abs());
    if energy_residual_j > controls.energy_tolerance_j {
        return Err(PistonError::Balance {
            momentum_residual_n_s: step
                .piston
                .momentum_residual_n_s
                .abs()
                .max(step.piston.combined_momentum_residual_n_s.abs()),
            energy_residual_j,
        });
    }
    Ok(step)
}
