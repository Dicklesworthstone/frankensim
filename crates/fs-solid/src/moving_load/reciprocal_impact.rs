//! Reciprocal Coulomb events between two continuing native shell trajectories.

use super::frictional_impact::{audit, dot, measured_law, rounded};
use super::{
    MovingShellImpulse, MovingShellModel, MovingShellMotion, ShellFrictionalImpactControls,
    ShellFrictionalImpactError, ShellFrictionalImpactSite, ShellLoadPoint, ShellPointImpulse,
    ShellSpatialResultant, poll,
};
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactResiduals, FrictionalContactResult,
    solve_frictional_contacts,
};
use fs_math::dd::Dd;

/// Original material witnesses in the same reference-world contact frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellReciprocalImpactSite {
    /// Receives +p; its normal points from the second shell toward the first.
    /// ID, frame, activation, friction and normal law retain their native meaning.
    pub first: ShellFrictionalImpactSite,
    /// Receives -p, without an added free couple or relocated attachment.
    pub second: ShellLoadPoint,
}

/// Physical geometry admission plus original endpoint/law audit controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellReciprocalImpactControls {
    pub impact: ShellFrictionalImpactControls,
    /// Maximum Cartesian separation of the original reference attachments [m].
    /// Any remaining impulse-weighted reference couple is independently gated.
    pub reference_point_tolerance_m: f64,
    /// Maximum separation along either contact tangent at the actual event [m].
    /// Normal separation uses the original activation/penetration law.
    pub tangential_position_tolerance_m: f64,
}
impl Default for ShellReciprocalImpactControls {
    fn default() -> Self {
        Self {
            impact: ShellFrictionalImpactControls::default(),
            reference_point_tolerance_m: 1e-10,
            tangential_position_tolerance_m: 1e-10,
        }
    }
}

/// One complete relative mobility bound to both original shell models.
#[derive(Debug)]
pub struct PreparedShellReciprocalImpacts<'a, 'b> {
    first_model: &'a MovingShellModel,
    second_model: &'b MovingShellModel,
    sites: Vec<ShellReciprocalImpactSite>,
    order: Vec<usize>,
    delassus: Vec<f64>,
    body_delassus: [Vec<f64>; 2],
    controls: ShellReciprocalImpactControls,
}
impl PreparedShellReciprocalImpacts<'_, '_> {
    pub fn sites(&self) -> &[ShellReciprocalImpactSite] {
        &self.sites
    }
    /// Complete row-major (J_A M_A,ff^-1 J_A^T + J_B M_B,ff^-1 J_B^T),
    /// contact-major [normal, tangent 1, tangent 2], in 1/kg.
    pub fn inverse_effective_mass(&self) -> &[f64] {
        &self.delassus
    }
}

/// Two accepted native impulses with one simultaneous contact law.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellReciprocalImpactReport {
    pub time_s: f64,
    /// Native determinacy and total-work witnesses remain authoritative.
    /// Redundant individual loads/support work shares remain representative.
    pub contact: FrictionalContactResult,
    pub relative_velocities_after_m_s: Vec<[f64; 3]>,
    pub maximum_speed_error_m_s: f64,
    /// Each body's actual point increment against its own original mobility.
    pub maximum_body_speed_error_m_s: [f64; 2],
    pub measured_residuals: FrictionalContactResiduals,
    /// First and second native receipts, or None for an exactly zero event.
    pub impulses: Option<[MovingShellImpulse; 2]>,
    /// Actual net commanded contact couple about the reference-world origin.
    /// Geometry tolerance alone does not authorize a finite internal couple.
    pub contact_reference_couple_nms: [f64; 3],
    /// Combined actual momentum change minus both external support impulses.
    /// Angular components are the original reference-frame inventory.
    pub momentum_residual: ShellSpatialResultant,
    pub kinetic_energy_change_j: f64,
    pub support_work_j: f64,
    pub relative_work_j: f64,
    pub dissipated_energy_j: f64,
    pub passivity_residual_j: f64,
    pub allowed_relative_work_j: f64,
    /// Kinetic change minus support and actual relative impulse work [J].
    pub work_balance_residual_j: f64,
}
type Result<T> = core::result::Result<T, ShellFrictionalImpactError>;

impl<'a> MovingShellMotion<'a> {
    /// Prepare reciprocal finite-inertia response using each original mass LU.
    ///
    /// Reference attachment agreement gives this linearized shell interface a
    /// common moment arm. It does not require zero displacement: both continuing
    /// shells may be deformed and translated. Original material witnesses stay
    /// fixed; finite-rotation spatial angular conservation is not asserted.
    /// Each model's point/matrix budgets and the native 21-contact cap apply.
    /// Both body matrices and their sum fit checked aggregate entry/byte bounds;
    /// retained individual responses prevent cancellation of endpoint errors.
    pub fn prepare_reciprocal_impacts<'b>(
        &mut self,
        second: &mut MovingShellMotion<'b>,
        sites: &[ShellReciprocalImpactSite],
        controls: ShellReciprocalImpactControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedShellReciprocalImpacts<'a, 'b>> {
        poll(&mut cancelled)?;
        for value in [
            controls.reference_point_tolerance_m,
            controls.tangential_position_tolerance_m,
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(ShellFrictionalImpactError::Invalid(
                    "finite nonnegative reciprocal geometry allowances required",
                ));
            }
        }
        // Bound all site-sized allocations before deriving either point list.
        if sites.len() > 21
            || sites.len() > controls.impact.contact.maximum_contacts
            || sites.len() > self.model.max_load_points
            || sites.len() > second.model.max_load_points
        {
            return Err(ShellFrictionalImpactError::Invalid(
                "reciprocal shell contact count budget",
            ));
        }
        let dimension = sites
            .len()
            .checked_mul(3)
            .ok_or(ShellFrictionalImpactError::Invalid(
                "reciprocal dimension overflow",
            ))?;
        let entries =
            dimension
                .checked_mul(dimension)
                .ok_or(ShellFrictionalImpactError::Invalid(
                    "reciprocal matrix overflow",
                ))?;
        let total_entries = entries
            .checked_mul(3)
            .ok_or(ShellFrictionalImpactError::Invalid(
                "reciprocal matrix count overflow",
            ))?;
        let bytes = total_entries
            .checked_mul(core::mem::size_of::<f64>())
            .ok_or(ShellFrictionalImpactError::Invalid(
                "reciprocal byte count overflow",
            ))?;
        if total_entries > self.model.plate.budget.max_matrix_entries
            || total_entries > second.model.plate.budget.max_matrix_entries
            || bytes > isize::MAX as usize
        {
            return Err(ShellFrictionalImpactError::Invalid(
                "reciprocal total matrix budget",
            ));
        }
        for site in sites {
            poll(&mut cancelled)?;
            let first_position = self.model.point(site.first.point)?.1;
            let second_position = second.model.point(site.second)?.1;
            for a in 0..3 {
                audit(
                    "shared reference contact point",
                    rounded(Dd::from_f64(first_position[a]) - Dd::from_f64(second_position[a]))?,
                    controls.reference_point_tolerance_m,
                )?;
            }
        }
        let first_sites: Vec<_> = sites.iter().map(|site| site.first).collect();
        let second_sites: Vec<_> = sites
            .iter()
            .map(|site| ShellFrictionalImpactSite {
                point: site.second,
                ..site.first
            })
            .collect();
        let first =
            self.prepare_frictional_impacts(&first_sites, controls.impact, &mut cancelled)?;
        let other =
            second.prepare_frictional_impacts(&second_sites, controls.impact, &mut cancelled)?;
        let mut delassus = Vec::with_capacity(entries);
        for (a, b) in first.delassus.iter().zip(&other.delassus) {
            delassus.push(rounded(Dd::from_f64(*a) + Dd::from_f64(*b))?);
        }
        poll(&mut cancelled)?;
        Ok(PreparedShellReciprocalImpacts {
            first_model: self.model,
            second_model: second.model,
            sites: sites.to_vec(),
            order: first.order,
            delassus,
            body_delassus: [first.delassus, other.delassus],
            controls,
        })
    }

    /// Resolve all paired contacts once and atomically accept both trajectories.
    ///
    /// Accepted times must match. Actual deformed witnesses must lie on the same
    /// normal line; their separation supplies the native gap, without position
    /// repair. The fixed orthonormal reference frame is retained throughout the
    /// instantaneous event. Reprepare when material witnesses or frames change.
    ///
    /// Equal/opposite point impulses drive two full-mass trials. Each native
    /// endpoint error is checked separately before combined momentum/work checks,
    /// so opposing numerical errors cannot hide one failed body. Original support
    /// impulses/work are external; there is no prescribed obstacle energy supply.
    /// Any refusal, including the second trial or final cancellation, publishes
    /// neither body. Held forces, displacement, time and support motion stay fixed.
    ///
    /// This shared-reference small-rotation interface does not add finite rigid
    /// motion, contact discovery, sustained friction, or per-contact heat authority.
    pub fn apply_reciprocal_impacts(
        &mut self,
        second: &mut MovingShellMotion<'_>,
        prepared: &PreparedShellReciprocalImpacts<'_, '_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ShellReciprocalImpactReport> {
        poll(&mut cancelled)?;
        if !core::ptr::eq(self.model, prepared.first_model)
            || !core::ptr::eq(second.model, prepared.second_model)
            || self.state.time_s != second.state.time_s
        {
            return Err(ShellFrictionalImpactError::Invalid(
                "both original models and synchronized accepted times required",
            ));
        }
        let controls = prepared.controls.impact;
        let count = prepared.sites.len();
        let mut contacts = Vec::with_capacity(count);
        let mut before_body_velocities = Vec::with_capacity(count);
        for site in &prepared.sites {
            poll(&mut cancelled)?;
            let a = self.observe(site.first.point)?;
            let b = second.observe(site.second)?;
            let mut separation = [0.0; 3];
            let mut relative_velocity = [0.0; 3];
            for k in 0..3 {
                separation[k] = difference(site.first.frame[k], a.position_m, b.position_m)?;
                relative_velocity[k] =
                    difference(site.first.frame[k], a.velocity_m_s, b.velocity_m_s)?;
                if k > 0 {
                    audit(
                        "current tangent witness",
                        separation[k],
                        prepared.controls.tangential_position_tolerance_m,
                    )?;
                }
            }
            before_body_velocities.push([a.velocity_m_s, b.velocity_m_s]);
            contacts.push(FrictionalContact {
                id: site.first.id,
                gap: separation[0],
                activation_tolerance: site.first.activation_tolerance_m,
                relative_velocity,
                friction_coefficient: site.first.friction_coefficient,
                normal_law: site.first.normal_law,
            });
        }
        let law = solve_frictional_contacts(
            &contacts,
            &prepared.delassus,
            controls.contact,
            &mut cancelled,
        )?;
        let trials = if law.impulses.iter().flatten().any(|p| *p != 0.0) {
            let mut first_impulses = Vec::with_capacity(count);
            let mut second_impulses = Vec::with_capacity(count);
            for &i in &prepared.order {
                poll(&mut cancelled)?;
                let site = prepared.sites[i];
                let mut world = [0.0; 3];
                for (axis, value) in world.iter_mut().enumerate() {
                    *value = rounded((0..3).fold(Dd::ZERO, |sum, k| {
                        sum + Dd::from_f64(law.impulses[i][k])
                            * Dd::from_f64(site.first.frame[k][axis])
                    }))?;
                }
                first_impulses.push(ShellPointImpulse {
                    point: site.first.point,
                    linear_impulse_ns: world,
                    angular_impulse_nms: [0.0; 3],
                });
                second_impulses.push(ShellPointImpulse {
                    point: site.second,
                    linear_impulse_ns: world.map(|p| -p),
                    angular_impulse_nms: [0.0; 3],
                });
            }
            Some((
                self.trial_impulses(&first_impulses, &mut cancelled)?,
                second.trial_impulses(&second_impulses, &mut cancelled)?,
            ))
        } else {
            None
        };
        let first_state = trials
            .as_ref()
            .map_or(&self.state, |trials| &trials.0.state);
        let second_state = trials
            .as_ref()
            .map_or(&second.state, |trials| &trials.1.state);
        let impulses = trials
            .as_ref()
            .map(|trials| [trials.0.report, trials.1.report]);
        let mut after = vec![[0.0; 3]; count];
        let mut maximum_speed_error_m_s = 0.0_f64;
        let mut maximum_body_speed_error_m_s = [0.0_f64; 2];
        let (mut relative_work, mut work_scale) = (Dd::ZERO, Dd::ZERO);
        for &i in &prepared.order {
            poll(&mut cancelled)?;
            let site = prepared.sites[i];
            let a = self.model.observe_with_supports(
                site.first.point,
                &first_state.displacement,
                &first_state.velocity,
                self.supports.as_ref(),
            )?;
            let b = second.model.observe_with_supports(
                site.second,
                &second_state.displacement,
                &second_state.velocity,
                second.supports.as_ref(),
            )?;
            for k in 0..3 {
                for (body, velocity) in [a.velocity_m_s, b.velocity_m_s].into_iter().enumerate() {
                    let mut predicted = Dd::ZERO;
                    for &j in &prepared.order {
                        for axis in 0..3 {
                            let p = if body == 0 {
                                law.impulses[j][axis]
                            } else {
                                -law.impulses[j][axis]
                            };
                            predicted = predicted
                                + Dd::from_f64(
                                    prepared.body_delassus[body]
                                        [(3 * i + k) * (3 * count) + 3 * j + axis],
                                ) * Dd::from_f64(p);
                        }
                    }
                    let actual = dot(site.first.frame[k], velocity)
                        - dot(site.first.frame[k], before_body_velocities[i][body]);
                    let defect = rounded(actual - predicted)?;
                    audit(
                        "individual contact velocity increment",
                        defect,
                        controls.speed_tolerance_m_s,
                    )?;
                    maximum_body_speed_error_m_s[body] =
                        maximum_body_speed_error_m_s[body].max(defect.abs());
                }
                let speed = difference(site.first.frame[k], a.velocity_m_s, b.velocity_m_s)?;
                after[i][k] = speed;
                let defect = rounded(
                    Dd::from_f64(speed) - Dd::from_f64(law.relative_velocities_after[i][k]),
                )?;
                audit(
                    "actual reciprocal contact speed",
                    defect,
                    controls.speed_tolerance_m_s,
                )?;
                maximum_speed_error_m_s = maximum_speed_error_m_s.max(defect.abs());
                let work = Dd::from_f64(law.impulses[i][k])
                    * (Dd::from_f64(0.5 * contacts[i].relative_velocity[k])
                        + Dd::from_f64(0.5 * speed));
                relative_work = relative_work + work;
                work_scale = work_scale + work.abs();
            }
        }
        let measured_residuals =
            measured_law(&contacts, &law, &after, controls.contact, &mut cancelled)?;
        let relative_work_j = rounded(relative_work)?;
        audit(
            "reciprocal relative work",
            rounded(relative_work - Dd::from_f64(law.relative_work))?,
            controls.work_tolerance_j,
        )?;
        let allowed_relative_work_j = rounded(
            Dd::from_f64(controls.contact.energy_tolerance)
                + Dd::from_f64(64.0 * f64::EPSILON) * work_scale,
        )?;
        audit(
            "positive reciprocal work",
            relative_work_j.max(0.0),
            allowed_relative_work_j,
        )?;
        let mut momentum_residual = ShellSpatialResultant::default();
        let mut contact_reference_couple_nms = [0.0; 3];
        let (mut kinetic_change, mut support_work, mut applied_work) =
            (Dd::ZERO, Dd::ZERO, Dd::ZERO);
        if let Some(reports) = &impulses {
            for report in reports {
                poll(&mut cancelled)?;
                for error in report.momentum_residual.linear {
                    audit(
                        "individual linear impulse",
                        error,
                        controls.linear_impulse_tolerance_ns,
                    )?;
                }
                for error in report.momentum_residual.angular {
                    audit(
                        "individual reference angular impulse",
                        error,
                        controls.angular_impulse_tolerance_nms,
                    )?;
                }
                audit(
                    "individual native impulse energy",
                    report.energy_residual_j,
                    controls.work_tolerance_j,
                )?;
                kinetic_change = kinetic_change + Dd::from_f64(report.kinetic_energy_change_j);
                support_work = support_work + Dd::from_f64(report.support_work_j);
                applied_work = applied_work + Dd::from_f64(report.applied_work_j);
            }
            for axis in 0..3 {
                let mut linear = Dd::ZERO;
                let mut angular = Dd::ZERO;
                let mut applied_linear = Dd::ZERO;
                let mut couple = Dd::ZERO;
                for report in reports {
                    linear = linear + Dd::from_f64(report.momentum_change.linear[axis])
                        - Dd::from_f64(report.support_impulse.linear[axis]);
                    angular = angular + Dd::from_f64(report.momentum_change.angular[axis])
                        - Dd::from_f64(report.support_impulse.angular[axis]);
                    applied_linear =
                        applied_linear + Dd::from_f64(report.applied_impulse.linear[axis]);
                    couple = couple + Dd::from_f64(report.applied_impulse.angular[axis]);
                }
                momentum_residual.linear[axis] = rounded(linear)?;
                momentum_residual.angular[axis] = rounded(angular)?;
                contact_reference_couple_nms[axis] = rounded(couple)?;
                audit(
                    "reciprocal commanded linear impulse",
                    rounded(applied_linear)?,
                    controls.linear_impulse_tolerance_ns,
                )?;
                audit(
                    "reference offset contact couple",
                    contact_reference_couple_nms[axis],
                    controls.angular_impulse_tolerance_nms,
                )?;
                audit(
                    "combined linear momentum",
                    momentum_residual.linear[axis],
                    controls.linear_impulse_tolerance_ns,
                )?;
                audit(
                    "combined reference angular momentum",
                    momentum_residual.angular[axis],
                    controls.angular_impulse_tolerance_nms,
                )?;
            }
        }
        audit(
            "reciprocal point impulse work",
            rounded(applied_work - relative_work)?,
            controls.work_tolerance_j,
        )?;
        let kinetic_energy_change_j = rounded(kinetic_change)?;
        let support_work_j = rounded(support_work)?;
        let work_balance_residual_j = rounded(kinetic_change - support_work - relative_work)?;
        audit(
            "reciprocal kinetic/support balance",
            work_balance_residual_j,
            controls.work_tolerance_j,
        )?;
        poll(&mut cancelled)?;
        if let Some((first, other)) = trials {
            self.accept_impulses(first);
            second.accept_impulses(other);
        }
        Ok(ShellReciprocalImpactReport {
            time_s: self.state.time_s,
            contact: law,
            relative_velocities_after_m_s: after,
            maximum_speed_error_m_s,
            maximum_body_speed_error_m_s,
            measured_residuals,
            impulses,
            contact_reference_couple_nms,
            momentum_residual,
            kinetic_energy_change_j,
            support_work_j,
            relative_work_j,
            dissipated_energy_j: (-relative_work_j).max(0.0),
            passivity_residual_j: relative_work_j.max(0.0),
            allowed_relative_work_j,
            work_balance_residual_j,
        })
    }
}
fn difference(axis: [f64; 3], first: [f64; 3], second: [f64; 3]) -> Result<f64> {
    rounded(dot(axis, first) - dot(axis, second))
}
