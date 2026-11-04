//! Simultaneous native normal/Coulomb shell impacts on prescribed planes.

use super::{
    MovingLoadError, MovingShellImpulse, MovingShellModel, MovingShellMotion, ShellLoadPoint,
    ShellObstacleMotion, ShellPointImpulse, finite, linear, poll,
};
pub use fs_contact::frictional_multicontact::NormalVelocityLaw;
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactError, FrictionalContactResiduals,
    FrictionalContactResult, solve_frictional_contacts,
};
use fs_math::{dd::Dd, det};

/// A fixed material attachment with an explicitly authored Coulomb coefficient.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellFrictionalImpactSite {
    /// Unique stable identity; native solution and aggregate scatter use this order.
    pub id: u64,
    pub point: ShellLoadPoint,
    /// Orthonormal Cartesian basis [normal, tangent 1, tangent 2].
    /// The normal points into the allowed half-space. No frame is inferred or
    /// normalized; either orthonormal tangent orientation describes the same disk.
    pub frame: [[f64; 3]; 3],
    /// Nonnegative proximity allowance [m], never a penetration repair.
    pub activation_tolerance_m: f64,
    pub friction_coefficient: f64,
    /// Existing native Newton or explicit normal endpoint target [m/s].
    pub normal_law: NormalVelocityLaw,
}

/// Native numerical/law controls and separate actual-state audit allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellFrictionalImpactControls {
    pub contact: FrictionalContactConfig,
    pub speed_tolerance_m_s: f64,
    pub linear_impulse_tolerance_ns: f64,
    pub angular_impulse_tolerance_nms: f64,
    pub work_tolerance_j: f64,
}
impl Default for ShellFrictionalImpactControls {
    fn default() -> Self {
        Self {
            contact: FrictionalContactConfig::default(),
            speed_tolerance_m_s: 1e-10,
            linear_impulse_tolerance_ns: 1e-10,
            angular_impulse_tolerance_nms: 1e-10,
            work_tolerance_j: 1e-10,
        }
    }
}

/// A refusal publishes no change to accepted motion, supports or held inputs.
#[derive(Debug)]
pub enum ShellFrictionalImpactError {
    Invalid(&'static str),
    Shell(MovingLoadError),
    Contact(FrictionalContactError),
    Audit {
        quantity: &'static str,
        error: f64,
        tolerance: f64,
    },
    Cancelled,
}
impl core::fmt::Display for ShellFrictionalImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native shell Coulomb impact refused: {self:?}")
    }
}
impl std::error::Error for ShellFrictionalImpactError {}
impl From<MovingLoadError> for ShellFrictionalImpactError {
    fn from(error: MovingLoadError) -> Self {
        if matches!(error, MovingLoadError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Shell(error)
        }
    }
}
impl From<FrictionalContactError> for ShellFrictionalImpactError {
    fn from(error: FrictionalContactError) -> Self {
        if matches!(error, FrictionalContactError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Contact(error)
        }
    }
}

/// Accepted native law and its independently observed physical endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellFrictionalImpactReport {
    pub time_s: f64,
    /// Original input order; native determinacy/witnesses are retained. On a
    /// redundant lane, individual impulses and obstacle/support work shares
    /// describe a representative, not uniquely determined loads or heat.
    pub contact: FrictionalContactResult,
    /// Actual stored motion in each site's [normal, tangent 1, tangent 2] frame.
    pub relative_velocities_after_m_s: Vec<[f64; 3]>,
    pub maximum_speed_error_m_s: f64,
    /// Physical residuals reconstructed from the actually stored endpoint.
    pub measured_residuals: FrictionalContactResiduals,
    /// Sum of all three impulse components paired with prescribed velocity [J].
    pub obstacle_work_j: f64,
    pub relative_work_j: f64,
    /// Lost relative mechanical energy, with no thermal/damage allocation.
    pub dissipated_energy_j: f64,
    pub passivity_residual_j: f64,
    /// Native energy allowance, including its work-scale roundoff term [J].
    pub allowed_relative_work_j: f64,
    /// Actual kinetic change minus support, obstacle and relative work [J].
    pub work_balance_residual_j: f64,
    /// Actual full-mass momentum, support impulse and separate support work.
    /// None means all contact impulses vanished and native state stayed exact.
    pub impulse: Option<MovingShellImpulse>,
}

/// Complete contact mobility bound to one original shell model and fixed frames.
/// Reference inertia, payloads, supports and material points retain the original
/// small-rotation scope. Reuse preparation while these inputs stay unchanged.
#[derive(Debug)]
pub struct PreparedShellFrictionalImpacts<'a> {
    model: &'a MovingShellModel,
    sites: Vec<ShellFrictionalImpactSite>,
    order: Vec<usize>,
    delassus: Vec<f64>,
    controls: ShellFrictionalImpactControls,
}
impl PreparedShellFrictionalImpacts<'_> {
    pub fn sites(&self) -> &[ShellFrictionalImpactSite] {
        &self.sites
    }
    /// Complete row-major J_f M_ff^-1 J_f^T [1/kg], contact-major and normal-first.
    /// All normal/tangent and intercontact blocks remain, including null directions.
    pub fn inverse_effective_mass(&self) -> &[f64] {
        &self.delassus
    }
}
type Result<T> = core::result::Result<T, ShellFrictionalImpactError>;

impl<'a> MovingShellMotion<'a> {
    /// Three native mass solves per contact, without another factorization.
    /// Original exact point kinematics/transpose and retained mass LU include
    /// off-center rotation, full payload inertia and prescribed-pin constraints.
    /// The native 21-contact cap and existing point/matrix/byte budgets apply.
    /// Preparation changes only solver scratch and is reusable across events.
    pub fn prepare_frictional_impacts(
        &mut self,
        sites: &[ShellFrictionalImpactSite],
        controls: ShellFrictionalImpactControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedShellFrictionalImpacts<'a>> {
        poll(&mut cancelled)?;
        let count = sites.len();
        if count > 21
            || count > controls.contact.maximum_contacts
            || count > self.model.max_load_points
        {
            return Err(ShellFrictionalImpactError::Invalid(
                "shell Coulomb contact count budget",
            ));
        }
        for value in [
            controls.speed_tolerance_m_s,
            controls.linear_impulse_tolerance_ns,
            controls.angular_impulse_tolerance_nms,
            controls.work_tolerance_j,
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(ShellFrictionalImpactError::Invalid(
                    "finite nonnegative physical audit allowances required",
                ));
            }
        }
        solve_frictional_contacts(&[], &[], controls.contact, &mut cancelled)?;
        let dimension = count
            .checked_mul(3)
            .ok_or(ShellFrictionalImpactError::Invalid(
                "mobility dimension overflow",
            ))?;
        let entries =
            dimension
                .checked_mul(dimension)
                .ok_or(ShellFrictionalImpactError::Invalid(
                    "mobility size overflow",
                ))?;
        let bytes = entries.checked_mul(core::mem::size_of::<f64>()).ok_or(
            ShellFrictionalImpactError::Invalid("mobility byte size overflow"),
        )?;
        if entries > self.model.plate.budget.max_matrix_entries || bytes > isize::MAX as usize {
            return Err(ShellFrictionalImpactError::Invalid(
                "Coulomb mobility matrix budget",
            ));
        }
        let mut order: Vec<_> = (0..count).collect();
        order.sort_unstable_by_key(|&i| sites[i].id);
        if order.windows(2).any(|w| sites[w[0]].id == sites[w[1]].id) {
            return Err(ShellFrictionalImpactError::Invalid("duplicate contact ID"));
        }
        let mut nodes = Vec::with_capacity(count);
        for site in sites {
            poll(&mut cancelled)?;
            for (i, axis) in site.frame.iter().enumerate() {
                finite(axis)?;
                for (j, other) in site.frame.iter().enumerate() {
                    let target = if i == j { 1.0 } else { 0.0 };
                    if (rounded(dot(*axis, *other))? - target).abs() > 64.0 * f64::EPSILON {
                        return Err(ShellFrictionalImpactError::Invalid(
                            "orthonormal normal-first contact frame required",
                        ));
                    }
                }
            }
            let valid_law = match site.normal_law {
                NormalVelocityLaw::NewtonRestitution { coefficient } => {
                    coefficient.is_finite() && (0.0..=1.0).contains(&coefficient)
                }
                NormalVelocityLaw::Target { speed } => speed.is_finite() && speed >= 0.0,
            };
            if !site.activation_tolerance_m.is_finite()
                || site.activation_tolerance_m < 0.0
                || !site.friction_coefficient.is_finite()
                || site.friction_coefficient < 0.0
                || !valid_law
            {
                return Err(ShellFrictionalImpactError::Invalid(
                    "finite nonnegative activation/friction and native normal law required",
                ));
            }
            nodes.push(self.model.point(site.point)?.0);
        }
        let mut delassus = vec![0.0; entries];
        let mut response = vec![0.0; self.state.velocity.len()];
        for &column in &order {
            let site = &sites[column];
            for (b, &direction) in site.frame.iter().enumerate() {
                poll(&mut cancelled)?;
                let load = self.model.project_components(
                    1,
                    core::iter::once((site.point, direction, [0.0; 3])),
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
                    for (a, &direction) in site.frame.iter().enumerate() {
                        delassus[(3 * row + a) * dimension + 3 * column + b] =
                            rounded(dot(direction, [value[0], value[1], value[2]]))?;
                    }
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(PreparedShellFrictionalImpacts {
            model: self.model,
            sites: sites.to_vec(),
            order,
            delassus,
            controls,
        })
    }

    /// Solve simultaneous normal and two-axis non-associated Coulomb impacts.
    ///
    /// Actual point positions define plane gaps; all geometry-active contacts
    /// enter one native solve with full cross mobility, even if initially
    /// separating. Native restitution/target, stick/slip and SPD/PSD determinacy
    /// remain authoritative; no sequential friction clipping or compliance.
    ///
    /// One aggregate native impulse trial is checked against actual rounded
    /// velocity, the physical contact law, full linear/reference-angular momentum
    /// and separate obstacle/support work before acceptance. Time, displacement,
    /// support motion, held force and step count stay fixed. The force owner
    /// refreshes acceleration/reactions only when a nonzero event is accepted.
    ///
    /// Callers own event timing, complete plane geometry and coefficient evidence.
    /// This is not finite-rotation contact, penetration repair, sustained friction,
    /// accumulated slip history, a partial-slip material model or a heat law.
    pub fn apply_frictional_impacts(
        &mut self,
        prepared: &PreparedShellFrictionalImpacts<'_>,
        obstacles: &[ShellObstacleMotion],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ShellFrictionalImpactReport> {
        poll(&mut cancelled)?;
        if !core::ptr::eq(self.model, prepared.model) || obstacles.len() != prepared.sites.len() {
            return Err(ShellFrictionalImpactError::Invalid(
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
                    * Dd::from_f64(site.frame[0][axis])
            });
            let mut surface = [Dd::ZERO; 3];
            let mut relative = [0.0; 3];
            for a in 0..3 {
                surface[a] = dot(site.frame[a], obstacle.velocity_m_s);
                rounded(surface[a])?;
                relative[a] = rounded(dot(site.frame[a], point.velocity_m_s) - surface[a])?;
            }
            contacts.push(FrictionalContact {
                id: site.id,
                gap: rounded(gap)?,
                activation_tolerance: site.activation_tolerance_m,
                relative_velocity: relative,
                friction_coefficient: site.friction_coefficient,
                normal_law: site.normal_law,
            });
            surface_speeds.push(surface);
        }
        let law = solve_frictional_contacts(
            &contacts,
            &prepared.delassus,
            controls.contact,
            &mut cancelled,
        )?;
        let trial = if law.impulses.iter().flatten().any(|p| *p != 0.0) {
            let mut impulses = Vec::with_capacity(count);
            for &i in &prepared.order {
                poll(&mut cancelled)?;
                let site = &prepared.sites[i];
                let mut linear_impulse_ns = [0.0; 3];
                for (axis, value) in linear_impulse_ns.iter_mut().enumerate() {
                    *value = rounded((0..3).fold(Dd::ZERO, |sum, a| {
                        sum + Dd::from_f64(law.impulses[i][a]) * Dd::from_f64(site.frame[a][axis])
                    }))?;
                }
                impulses.push(ShellPointImpulse {
                    point: site.point,
                    linear_impulse_ns,
                    angular_impulse_nms: [0.0; 3],
                });
            }
            Some(self.trial_impulses(&impulses, &mut cancelled)?)
        } else {
            None
        };
        let state = trial.as_ref().map_or(&self.state, |trial| &trial.state);
        let impulse = trial.as_ref().map(|trial| trial.report);
        let mut after = vec![[0.0; 3]; count];
        let mut maximum_speed_error_m_s = 0.0_f64;
        let (mut obstacle_work, mut relative_work, mut work_scale) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
        for &i in &prepared.order {
            poll(&mut cancelled)?;
            let site = &prepared.sites[i];
            let point = self.model.observe_with_supports(
                site.point,
                &state.displacement,
                &state.velocity,
                self.supports.as_ref(),
            )?;
            for a in 0..3 {
                let speed = rounded(dot(site.frame[a], point.velocity_m_s) - surface_speeds[i][a])?;
                after[i][a] = speed;
                let defect = rounded(
                    Dd::from_f64(speed) - Dd::from_f64(law.relative_velocities_after[i][a]),
                )?;
                audit("actual contact speed", defect, controls.speed_tolerance_m_s)?;
                maximum_speed_error_m_s = maximum_speed_error_m_s.max(defect.abs());
                let p = Dd::from_f64(law.impulses[i][a]);
                obstacle_work = obstacle_work + p * surface_speeds[i][a];
                let term = p
                    * (Dd::from_f64(0.5 * contacts[i].relative_velocity[a])
                        + Dd::from_f64(0.5 * speed));
                relative_work = relative_work + term;
                work_scale = work_scale + term.abs();
            }
        }
        let measured_residuals =
            measured_law(&contacts, &law, &after, controls.contact, &mut cancelled)?;
        let obstacle_work_j = rounded(obstacle_work)?;
        let relative_work_j = rounded(relative_work)?;
        audit(
            "contact relative work",
            rounded(relative_work - Dd::from_f64(law.relative_work))?,
            controls.work_tolerance_j,
        )?;
        let allowed_relative_work_j = rounded(
            Dd::from_f64(controls.contact.energy_tolerance)
                + Dd::from_f64(64.0 * f64::EPSILON) * work_scale,
        )?;
        audit(
            "positive relative work",
            relative_work_j.max(0.0),
            allowed_relative_work_j,
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
        Ok(ShellFrictionalImpactReport {
            time_s: self.state.time_s,
            contact: law,
            relative_velocities_after_m_s: after,
            maximum_speed_error_m_s,
            measured_residuals,
            obstacle_work_j,
            relative_work_j,
            dissipated_energy_j: (-relative_work_j).max(0.0),
            passivity_residual_j: relative_work_j.max(0.0),
            allowed_relative_work_j,
            work_balance_residual_j,
            impulse,
        })
    }
}

fn measured_law(
    contacts: &[FrictionalContact],
    law: &FrictionalContactResult,
    after: &[[f64; 3]],
    config: FrictionalContactConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<FrictionalContactResiduals> {
    let mut r = FrictionalContactResiduals::default();
    for &i in &law.active_indices {
        poll(cancelled)?;
        let (p, v) = (law.impulses[i], after[i]);
        let g = v[0] - law.normal_targets[i];
        let radius = contacts[i].friction_coefficient * p[0];
        let speed = det::hypot(v[1], v[2]);
        let normal_work = p[0] * g;
        let cone = det::hypot(p[1], p[2]) - radius;
        finite(&[g, radius, speed, normal_work, cone])?;
        r.normal_impulse_violation = r.normal_impulse_violation.max((-p[0]).max(0.0));
        r.normal_velocity_violation = r.normal_velocity_violation.max((-g).max(0.0));
        r.normal_complementarity = r.normal_complementarity.max(normal_work.abs());
        r.cone_violation = r.cone_violation.max(cone.max(0.0));
        let duality = rounded(
            Dd::from_f64(radius) * Dd::from_f64(speed)
                + Dd::from_f64(p[1]) * Dd::from_f64(v[1])
                + Dd::from_f64(p[2]) * Dd::from_f64(v[2]),
        )?;
        r.tangential_duality_gap = r.tangential_duality_gap.max(duality.abs());
        if speed > config.velocity_tolerance {
            let residual = det::hypot(
                p[1] + radius * (v[1] / speed),
                p[2] + radius * (v[2] / speed),
            );
            finite(&[residual])?;
            r.sliding_impulse_residual = r.sliding_impulse_residual.max(residual);
        }
    }
    if r.normal_impulse_violation > config.impulse_tolerance
        || r.normal_velocity_violation > config.velocity_tolerance
        || r.normal_complementarity > config.complementarity_tolerance
        || r.cone_violation > config.impulse_tolerance
        || r.sliding_impulse_residual > config.impulse_tolerance
        || r.tangential_duality_gap > config.complementarity_tolerance
    {
        return Err(FrictionalContactError::ResidualNotSatisfied(r).into());
    }
    Ok(r)
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
    finite(&[error, tolerance])?;
    if error.abs() > tolerance {
        Err(ShellFrictionalImpactError::Audit {
            quantity,
            error,
            tolerance,
        })
    } else {
        Ok(())
    }
}
