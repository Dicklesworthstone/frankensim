//! Simultaneous native shell impacts against explicitly prescribed planes.

use super::{
    MovingLoadError, MovingShellImpulse, MovingShellModel, MovingShellMotion, ShellLoadPoint,
    ShellPointImpulse, finite, linear, poll,
};
use fs_contact::multicontact::{
    MultiImpactConfig, MultiImpactError, MultiImpactResult, NormalContactImpact,
    solve_normal_impacts,
};
use fs_math::dd::Dd;

/// One fixed material attachment and its declared separating normal law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellNormalImpactSite {
    pub point: ShellLoadPoint,
    /// Unit Cartesian normal into the allowed half-space; never normalized here.
    pub normal: [f64; 3],
    /// Nonnegative proximity allowance [m], never a penetration repair.
    pub activation_tolerance_m: f64,
    /// Native Newton restitution in [0, 1]; zero is perfectly inelastic.
    pub restitution: f64,
}

/// Instantaneous prescribed plane geometry and local surface velocity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellObstacleMotion {
    /// A point on the plane at the shell's accepted time, in world metres.
    pub point_m: [f64; 3],
    /// Surface velocity at this contact, held throughout the event [m/s].
    /// Tangential motion supplies no work to this frictionless normal law.
    pub velocity_m_s: [f64; 3],
}

/// Native law controls and separate measured-state acceptance allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellNormalImpactControls {
    pub contact: MultiImpactConfig,
    pub speed_tolerance_m_s: f64,
    pub linear_impulse_tolerance_ns: f64,
    pub angular_impulse_tolerance_nms: f64,
    pub work_tolerance_j: f64,
}
impl Default for ShellNormalImpactControls {
    fn default() -> Self {
        Self {
            contact: MultiImpactConfig::default(),
            speed_tolerance_m_s: 1e-10,
            linear_impulse_tolerance_ns: 1e-10,
            angular_impulse_tolerance_nms: 1e-10,
            work_tolerance_j: 1e-10,
        }
    }
}

/// Every refusal leaves the accepted shell motion and held inputs unchanged.
#[derive(Debug)]
pub enum ShellNormalImpactError {
    Invalid(&'static str),
    Shell(MovingLoadError),
    Contact(MultiImpactError),
    Audit {
        quantity: &'static str,
        error: f64,
        tolerance: f64,
    },
    Cancelled,
}
impl core::fmt::Display for ShellNormalImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native shell normal impact refused: {self:?}")
    }
}
impl std::error::Error for ShellNormalImpactError {}
impl From<MovingLoadError> for ShellNormalImpactError {
    fn from(error: MovingLoadError) -> Self {
        if matches!(error, MovingLoadError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Shell(error)
        }
    }
}
impl From<MultiImpactError> for ShellNormalImpactError {
    fn from(error: MultiImpactError) -> Self {
        if matches!(error, MultiImpactError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Contact(error)
        }
    }
}

/// One accepted simultaneous event, with native determinacy retained.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellNormalImpactReport {
    pub time_s: f64,
    /// Redundant impulses are only a checked representative. They do not
    /// determine individual pressures or separate obstacle/support work shares.
    pub contact: MultiImpactResult,
    pub relative_normal_speeds_after_m_s: Vec<f64>,
    pub maximum_speed_error_m_s: f64,
    pub maximum_complementarity_residual_m_s: f64,
    /// Sum p*n dot prescribed surface velocity, ON the shell [J].
    /// On a redundant lane this and native support work describe the selected
    /// representative; their combined energy balance is independently checked.
    pub obstacle_work_j: f64,
    /// Actual midpoint relative velocity paired with normal impulse [J].
    pub relative_work_j: f64,
    /// Lost relative mechanical energy; no thermal or damage partition.
    pub dissipated_energy_j: f64,
    pub passivity_residual_j: f64,
    /// Actual kinetic change minus support, obstacle and relative work [J].
    pub work_balance_residual_j: f64,
    /// Includes actual full-mass momentum, support reaction and separate work.
    /// None means every contact impulse was zero and native state stayed exact.
    pub impulse: Option<MovingShellImpulse>,
}

/// Prepared complete normal mobility, bound to the original shell model.
/// Points/normals and reference inertia stay fixed; obstacle samples may move.
#[derive(Debug)]
pub struct PreparedShellNormalImpacts<'a> {
    model: &'a MovingShellModel,
    sites: Vec<ShellNormalImpactSite>,
    delassus: Vec<f64>,
    controls: ShellNormalImpactControls,
}
impl PreparedShellNormalImpacts<'_> {
    pub fn sites(&self) -> &[ShellNormalImpactSite] {
        &self.sites
    }
    /// Complete row-major J_f M_ff^-1 J_f^T [1/kg], including cross contacts.
    pub fn inverse_effective_mass(&self) -> &[f64] {
        &self.delassus
    }
}

type Result<T> = core::result::Result<T, ShellNormalImpactError>;

impl<'a> MovingShellMotion<'a> {
    /// Prepare one native mass solve per normal and retain all cross response.
    ///
    /// This is the normal contraction of the complete point inverse inertia.
    /// Existing point projection, material kinematics and retained mass LU own
    /// the map; no diagonal effective masses or new factorization are used.
    /// The native contact cap (64), original point cap and checked original
    /// matrix-entry/byte bounds apply before allocation. Preparation changes
    /// only solver scratch and can be reused across this model's trajectories.
    pub fn prepare_normal_impacts(
        &mut self,
        sites: &[ShellNormalImpactSite],
        controls: ShellNormalImpactControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedShellNormalImpacts<'a>> {
        poll(&mut cancelled)?;
        let count = sites.len();
        if count > 64
            || count > controls.contact.solver.max_variables
            || count > self.model.max_load_points
        {
            return Err(ShellNormalImpactError::Invalid(
                "shell normal contact count budget",
            ));
        }
        for value in [
            controls.speed_tolerance_m_s,
            controls.linear_impulse_tolerance_ns,
            controls.angular_impulse_tolerance_nms,
            controls.work_tolerance_j,
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(ShellNormalImpactError::Invalid(
                    "finite nonnegative physical audit allowances required",
                ));
            }
        }
        // The existing owner validates its own numerical configuration.
        solve_normal_impacts(&[], &[], controls.contact, &mut cancelled)?;
        let entries = count
            .checked_mul(count)
            .ok_or(ShellNormalImpactError::Invalid(
                "normal mobility size overflow",
            ))?;
        let bytes = entries.checked_mul(core::mem::size_of::<f64>()).ok_or(
            ShellNormalImpactError::Invalid("normal mobility byte size overflow"),
        )?;
        if entries > self.model.plate.budget.max_matrix_entries || bytes > isize::MAX as usize {
            return Err(ShellNormalImpactError::Invalid(
                "normal mobility matrix budget",
            ));
        }
        let mut nodes = Vec::with_capacity(count);
        for site in sites {
            poll(&mut cancelled)?;
            finite(&site.normal)?;
            if (rounded(dot(site.normal, site.normal))? - 1.0).abs() > 64.0 * f64::EPSILON
                || !site.activation_tolerance_m.is_finite()
                || site.activation_tolerance_m < 0.0
                || !site.restitution.is_finite()
                || !(0.0..=1.0).contains(&site.restitution)
            {
                return Err(ShellNormalImpactError::Invalid(
                    "unit normal, nonnegative activation and restitution in [0,1] required",
                ));
            }
            nodes.push(self.model.point(site.point)?.0);
        }
        let mut delassus = vec![0.0; entries];
        let mut response = vec![0.0; self.state.velocity.len()];
        for (column, site) in sites.iter().enumerate() {
            poll(&mut cancelled)?;
            let load = self.model.project_components(
                1,
                core::iter::once((site.point, site.normal, [0.0; 3])),
                &mut cancelled,
            )?;
            self.mass_solver
                .solve_factored_into(&load.reduced, &mut response)
                .map_err(linear)?;
            finite(&response)?;
            for (row, (site, &nodes)) in sites.iter().zip(&nodes).enumerate() {
                let value = self
                    .model
                    .material_components(site.point, nodes, &response, None);
                finite(&value)?;
                delassus[row * count + column] =
                    rounded(dot(site.normal, [value[0], value[1], value[2]]))?;
            }
        }
        poll(&mut cancelled)?;
        Ok(PreparedShellNormalImpacts {
            model: self.model,
            sites: sites.to_vec(),
            delassus,
            controls,
        })
    }

    /// Resolve every supplied plane contact in one native simultaneous event.
    ///
    /// Gaps and relative speeds come from actual material-point observations.
    /// Every geometry-active row participates, including initially separating
    /// ones. The original contact owner supplies restitution, complementarity,
    /// passivity and SPD/PSD determinacy; no pairwise sequence or new law is used.
    ///
    /// One native aggregate impulse trial retains payload and constrained mass
    /// coupling. Its actually rounded velocity, endpoint contact law, full
    /// linear/reference-angular momentum and obstacle/support/relative work
    /// must pass before publication. Held forces, support motion, displacement,
    /// accepted time and step count stay fixed. A no-impulse event changes none
    /// of the native state, including acceleration and regular reactions.
    ///
    /// The caller owns event timing and completeness of the prescribed plane
    /// set. No CCD, penetration repair, sustained force, finite-rotation shell
    /// or contact-loss-to-heat allocation is inferred.
    pub fn apply_normal_impacts(
        &mut self,
        prepared: &PreparedShellNormalImpacts<'_>,
        obstacles: &[ShellObstacleMotion],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ShellNormalImpactReport> {
        poll(&mut cancelled)?;
        if !core::ptr::eq(self.model, prepared.model) || obstacles.len() != prepared.sites.len() {
            return Err(ShellNormalImpactError::Invalid(
                "original model identity and complete obstacle samples required",
            ));
        }
        let controls = prepared.controls;
        let count = prepared.sites.len();
        let mut contacts = Vec::with_capacity(count);
        let mut surface_speeds = Vec::with_capacity(count);
        for (site, obstacle) in prepared.sites.iter().zip(obstacles) {
            poll(&mut cancelled)?;
            finite(&obstacle.point_m)?;
            finite(&obstacle.velocity_m_s)?;
            let point = self.observe(site.point)?;
            let gap = (0..3).fold(Dd::ZERO, |sum, axis| {
                sum + (Dd::from_f64(point.position_m[axis]) - Dd::from_f64(obstacle.point_m[axis]))
                    * Dd::from_f64(site.normal[axis])
            });
            let surface = dot(site.normal, obstacle.velocity_m_s);
            rounded(surface)?;
            contacts.push(NormalContactImpact {
                gap: rounded(gap)?,
                activation_tolerance: site.activation_tolerance_m,
                relative_normal_speed: rounded(dot(site.normal, point.velocity_m_s) - surface)?,
                restitution: site.restitution,
            });
            surface_speeds.push(surface);
        }
        let law = solve_normal_impacts(
            &contacts,
            &prepared.delassus,
            controls.contact,
            &mut cancelled,
        )?;
        let trial = if law.normal_impulses.iter().any(|p| *p != 0.0) {
            let impulses: Vec<_> = prepared
                .sites
                .iter()
                .zip(&law.normal_impulses)
                .map(|(site, &p)| ShellPointImpulse {
                    point: site.point,
                    linear_impulse_ns: site.normal.map(|n| n * p),
                    angular_impulse_nms: [0.0; 3],
                })
                .collect();
            Some(self.trial_impulses(&impulses, &mut cancelled)?)
        } else {
            None
        };
        let state = trial.as_ref().map_or(&self.state, |trial| &trial.state);
        let impulse = trial.as_ref().map(|trial| trial.report);
        let mut after = Vec::with_capacity(count);
        let (mut maximum_speed_error_m_s, mut maximum_complementarity_residual_m_s) =
            (0.0_f64, 0.0_f64);
        let (mut obstacle_work, mut relative_work) = (Dd::ZERO, Dd::ZERO);
        let contact_tolerance = rounded(
            Dd::from_f64(controls.speed_tolerance_m_s)
                + Dd::from_f64(controls.contact.solver.gradient_tolerance),
        )?;
        for (i, site) in prepared.sites.iter().enumerate() {
            poll(&mut cancelled)?;
            let point = self.model.observe_with_supports(
                site.point,
                &state.displacement,
                &state.velocity,
                self.supports.as_ref(),
            )?;
            let speed = rounded(dot(site.normal, point.velocity_m_s) - surface_speeds[i])?;
            let defect =
                rounded(Dd::from_f64(speed) - Dd::from_f64(law.relative_normal_speeds_after[i]))?;
            audit("actual contact speed", defect, controls.speed_tolerance_m_s)?;
            maximum_speed_error_m_s = maximum_speed_error_m_s.max(defect.abs());
            if contacts[i].gap <= contacts[i].activation_tolerance {
                let gradient =
                    rounded(Dd::from_f64(speed) - Dd::from_f64(law.restitution_targets[i]))?;
                let residual = if law.normal_impulses[i] > 0.0 {
                    gradient.abs()
                } else {
                    (-gradient).max(0.0)
                }
                .max((-speed).max(0.0));
                audit(
                    "actual contact complementarity",
                    residual,
                    contact_tolerance,
                )?;
                maximum_complementarity_residual_m_s =
                    maximum_complementarity_residual_m_s.max(residual);
            }
            let p = Dd::from_f64(law.normal_impulses[i]);
            obstacle_work = obstacle_work + p * surface_speeds[i];
            relative_work = relative_work
                + p * (Dd::from_f64(0.5 * contacts[i].relative_normal_speed)
                    + Dd::from_f64(0.5 * speed));
            after.push(speed);
        }
        let obstacle_work_j = rounded(obstacle_work)?;
        let relative_work_j = rounded(relative_work)?;
        audit(
            "contact relative work",
            rounded(relative_work - Dd::from_f64(law.relative_work))?,
            controls.work_tolerance_j,
        )?;
        audit(
            "positive relative work",
            relative_work_j.max(0.0),
            controls.contact.energy_tolerance,
        )?;
        let applied_work = impulse.map_or(0.0, |p| p.applied_work_j);
        let support_work = impulse.map_or(0.0, |p| p.support_work_j);
        let kinetic_change = impulse.map_or(0.0, |p| p.kinetic_energy_change_j);
        audit(
            "point impulse work",
            rounded(Dd::from_f64(applied_work) - obstacle_work - relative_work)?,
            controls.work_tolerance_j,
        )?;
        let work_balance_residual_j = rounded(
            Dd::from_f64(kinetic_change)
                - Dd::from_f64(support_work)
                - obstacle_work
                - relative_work,
        )?;
        audit(
            "kinetic/source/support balance",
            work_balance_residual_j,
            controls.work_tolerance_j,
        )?;
        if let Some(impulse) = impulse {
            for error in impulse.momentum_residual.linear {
                audit(
                    "linear impulse",
                    error,
                    controls.linear_impulse_tolerance_ns,
                )?;
            }
            for error in impulse.momentum_residual.angular {
                audit(
                    "reference angular impulse",
                    error,
                    controls.angular_impulse_tolerance_nms,
                )?;
            }
            audit(
                "native impulse energy",
                impulse.energy_residual_j,
                controls.work_tolerance_j,
            )?;
        }
        poll(&mut cancelled)?;
        if let Some(trial) = trial {
            self.accept_impulses(trial);
        }
        Ok(ShellNormalImpactReport {
            time_s: self.state.time_s,
            contact: law,
            relative_normal_speeds_after_m_s: after,
            maximum_speed_error_m_s,
            maximum_complementarity_residual_m_s,
            obstacle_work_j,
            relative_work_j,
            dissipated_energy_j: (-relative_work_j).max(0.0),
            passivity_residual_j: relative_work_j.max(0.0),
            work_balance_residual_j,
            impulse,
        })
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> Dd {
    a.into_iter().zip(b).fold(Dd::ZERO, |sum, (a, b)| {
        sum + Dd::from_f64(a) * Dd::from_f64(b)
    })
}
fn rounded(value: Dd) -> Result<f64> {
    let result = value.to_f64();
    finite(&[value.hi, value.lo, result])?;
    Ok(result)
}
fn audit(quantity: &'static str, error: f64, tolerance: f64) -> Result<()> {
    finite(&[error])?;
    if error.abs() > tolerance {
        Err(ShellNormalImpactError::Audit {
            quantity,
            error,
            tolerance,
        })
    } else {
        Ok(())
    }
}
