//! Native cohesive trajectories through fixed affine contact geometry.
//!
//! Every free segment is the original implicit midpoint/material-path step.
//! Its displacement path and each supplied contact gap are affine, so checking
//! both actual endpoints admits that entire declared numerical path. A refused
//! path shortens the time proposal and RECOMPUTES native mechanics from the
//! unchanged predecessor. Rejected damage or interpolated states never publish.
//!
//! Contact frames, sparse traces and prescribed contact velocities stay fixed
//! over this interval. Positive activation distances are explicit spatial
//! approximations. Maximum-step and activation-distance refinement are needed
//! for trajectory accuracy; this is not exact-flow collision detection or a
//! general sustained-contact integrator. Loaded rest and unresolved progress
//! refuse, without projection, artificial support or repeated zero-time kicks.
//! Generalized work follows the supplied traces; Cartesian force/torque claims
//! additionally require the caller's actual geometric point maps.

use fs_contact::frictional_multicontact::FrictionalContactResult;
use fs_exec::Cx;
use fs_solid::cohesive::Error as CohesiveError;
use fs_solid::cohesive::dynamics::{CohesiveMotion, DynamicsError, DynamicsStep};
use fs_sparse::precond::Precond;

use super::impact::frictional::{
    FrictionalContact, FrictionalImpactControls, impacts_candidate_controlled,
};
use super::impact::{self, ImpactError};

mod geometry;
pub use geometry::GeometryError;

/// Constant reduced force and the original prescribed jump offsets at the
/// requested interval endpoint. Offsets evolve affinely from accepted native
/// offsets; a shorter trial receives its own correct physical endpoint.
#[derive(Clone, Copy, Debug)]
pub struct BodyLoad<'l> {
    pub mean_force_n: &'l [f64],
    pub endpoint_prescribed_jumps_m: &'l [[f64; 3]],
}

#[derive(Clone, Copy, Debug)]
pub struct TrajectoryControls {
    pub impact: FrictionalImpactControls,
    pub maximum_step_s: f64,
    pub minimum_step_s: f64,
    pub maximum_accepted_steps: usize,
    pub maximum_body_trials: usize,
    pub maximum_refinements_per_step: usize,
    pub maximum_impacts: usize,
    /// Native residual allowance reserved before each original body trial.
    /// Successful trials refund unused work; refused trials consume the reserve.
    pub maximum_residual_evaluations: usize,
    pub maximum_newton_iterations: usize,
    /// Adapter numeric buffers and retained native contact reports. Original
    /// material/operator storage and each owner's transient solver workspace
    /// remain governed by that owner's controls. Each retained event reserves
    /// the contact-response owner's complete numerical workspace conservatively.
    pub maximum_workspace_entries: usize,
    pub closing_velocity_tolerance_m_s: f64,
    /// Absolute whole-window work gate, also applied to sums of ABSOLUTE native
    /// step/event defects so opposing local errors cannot cancel into acceptance.
    pub energy_tolerance_j: f64,
}

/// Original native determinacy and loss receipt at the actual common clock.
/// Impulses remain in the caller's original contact order. Representative-only
/// native impulses do not acquire individual heat or load authority here.
#[derive(Debug)]
pub struct TimedImpactReport {
    pub time_s: f64,
    pub contact: FrictionalContactResult,
    pub kinetic_change_j: f64,
    pub prescribed_work_j: f64,
    pub work_energy_residual_j: f64,
    pub dissipation_residual_j: f64,
    pub absolute_body_work_residual_j: f64,
}

#[derive(Debug)]
pub struct TrajectoryReport {
    pub elapsed_s: f64,
    pub accepted_steps: usize,
    pub body_trials: usize,
    pub rejected_trials: usize,
    /// Charged native work; a refused trial consumes its reserved maximum.
    pub residual_evaluations: usize,
    /// Charged native work, with the same refused-trial reservation policy.
    pub newton_iterations: usize,
    pub contact_ids: Vec<u64>,
    pub events: Vec<TimedImpactReport>,
    pub initial_gaps_m: Vec<f64>,
    pub final_gaps_m: Vec<f64>,
    pub minimum_gap_m: f64,
    pub external_work_j: f64,
    pub prescribed_jump_work_j: f64,
    pub contact_prescribed_work_j: f64,
    pub kinetic_change_j: f64,
    pub bulk_stored_change_j: f64,
    pub cohesive_stored_change_j: f64,
    pub fracture_loss_j: f64,
    pub viscous_loss_j: f64,
    pub internal_friction_loss_j: f64,
    /// Native instantaneous impact loss, with no inferred thermal partition.
    pub impact_loss_j: f64,
    /// All supplied work minus actual storage changes and physical losses.
    pub energy_residual_j: f64,
    pub absolute_smooth_energy_defect_j: f64,
    pub absolute_impact_energy_defect_j: f64,
}

pub struct TrajectoryCandidate<'m, 'a> {
    pub motions: Vec<CohesiveMotion<'m, 'a>>,
    pub report: TrajectoryReport,
}

#[derive(Debug)]
pub enum TrajectoryError {
    Invalid(&'static str),
    Budget(&'static str),
    NonFinite(&'static str),
    Cancelled,
    Geometry(GeometryError),
    Impact(ImpactError),
    Body {
        body: usize,
        error: DynamicsError,
    },
    /// The original free/contact split cannot supply an admitted positive path.
    /// This includes unresolved loaded rest; no sustained reaction is invented.
    NoProgress {
        time_s: f64,
        trial_step_s: f64,
    },
    Refinement {
        time_s: f64,
        attempts: usize,
    },
    Energy {
        quantity: &'static str,
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl core::fmt::Display for TrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cohesive contact trajectory refused: {self:?}")
    }
}
impl std::error::Error for TrajectoryError {}
impl From<GeometryError> for TrajectoryError {
    fn from(error: GeometryError) -> Self {
        Self::Geometry(error)
    }
}
impl From<ImpactError> for TrajectoryError {
    fn from(error: ImpactError) -> Self {
        Self::Impact(error)
    }
}

/// Advance one complete requested interval, preserving all original owners.
/// No successor is returned until every segment, impact and final work audit
/// passes. `gap_m` describes the input state; prescribed contact velocity and
/// each sparse contact frame remain constant through this interval.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    contacts: &[FrictionalContact<'_>],
    loads: &[BodyLoad<'_>],
    preconditioners: &[&dyn Precond],
    duration_s: f64,
    controls: TrajectoryControls,
    cx: &Cx<'_>,
) -> Result<TrajectoryCandidate<'m, 'a>, TrajectoryError> {
    advance_candidate_controlled(
        motions,
        contacts,
        loads,
        preconditioners,
        duration_s,
        controls,
        &mut || cx.checkpoint().is_err(),
    )
}

#[allow(clippy::too_many_arguments)]
pub fn advance_candidate_controlled<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    contacts: &[FrictionalContact<'_>],
    loads: &[BodyLoad<'_>],
    preconditioners: &[&dyn Precond],
    duration_s: f64,
    controls: TrajectoryControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<TrajectoryCandidate<'m, 'a>, TrajectoryError> {
    admit(
        motions,
        contacts,
        loads,
        preconditioners,
        duration_s,
        controls,
        cancelled,
    )?;
    let start_s = motions[0].time_s();
    let target_s = finite(start_s + duration_s)?;
    let actual_duration_s = finite(target_s - start_s)?;
    if target_s <= start_s {
        return Err(TrajectoryError::NoProgress {
            time_s: start_s,
            trial_step_s: duration_s,
        });
    }
    let reference_positions: Vec<_> = motions
        .iter()
        .map(|m| m.state().displacements_m())
        .collect();
    geometry::admit(
        contacts,
        &reference_positions,
        controls.impact.response.maximum_trace_terms,
        controls.maximum_workspace_entries,
        cancelled,
    )?;
    let mut initial_offsets = Vec::with_capacity(motions.len());
    let mut maximum_step_s = controls.maximum_step_s;
    for (body, motion) in motions.iter().enumerate() {
        poll(cancelled)?;
        maximum_step_s = maximum_step_s.min(motion.dynamics_config().maximum_step_s);
        initial_offsets.push(
            motion
                .prescribed_jumps_m_controlled(cancelled)
                .map_err(|error| TrajectoryError::Body { body, error })?,
        );
    }
    if maximum_step_s < controls.minimum_step_s {
        return Err(TrajectoryError::Invalid(
            "native maximum step is below minimum step",
        ));
    }
    let mut current: Vec<_> = motions.iter().map(|m| (**m).clone()).collect();
    let mut sampled = sample(contacts, &reference_positions, &current, 0.0, cancelled)?;
    let initial_gaps_m = sampled.gaps_m.clone();
    let mut minimum_gap_m = sampled.minimum_gap_m;
    let mut work = Work::default();
    let mut balance = Balance::default();
    let mut events = Vec::new();
    let mut accepted_steps = 0;
    let mut rejected_trials = 0;
    let mut last_impact_time = None;

    loop {
        poll(cancelled)?;
        let time_s = current[0].time_s();
        if sampled.gaps_m.iter().any(|gap| *gap < 0.0) {
            return Err(TrajectoryError::Invalid(
                "accepted contact geometry penetrates",
            ));
        }
        if sampled
            .closing_contact(contacts, controls.closing_velocity_tolerance_m_s)
            .is_some()
        {
            if last_impact_time == Some(time_s) {
                return Err(TrajectoryError::NoProgress {
                    time_s,
                    trial_step_s: 0.0,
                });
            }
            if events.len() >= controls.maximum_impacts {
                return Err(TrajectoryError::Budget("impact events"));
            }
            let active: Vec<_> = contacts
                .iter()
                .zip(&sampled.gaps_m)
                .map(|(contact, &gap_m)| FrictionalContact { gap_m, ..*contact })
                .collect();
            let borrowed: Vec<_> = current.iter().collect();
            let event = impacts_candidate_controlled(
                &borrowed,
                &active,
                preconditioners,
                controls.impact,
                cancelled,
            )?;
            let r = event.report.response;
            balance.add_impact(&r)?;
            events.push(TimedImpactReport {
                time_s,
                kinetic_change_j: r.kinetic_change_j,
                prescribed_work_j: r.prescribed_work_j,
                work_energy_residual_j: r.work_energy_residual_j,
                dissipation_residual_j: r.dissipation_residual_j,
                absolute_body_work_residual_j: r.absolute_body_work_residual_j,
                contact: r.contact,
            });
            current = event.motions;
            sampled = sample(
                contacts,
                &reference_positions,
                &current,
                time_s - start_s,
                cancelled,
            )?;
            last_impact_time = Some(time_s);
            continue;
        }
        if time_s == target_s {
            break;
        }
        if time_s > target_s {
            return Err(TrajectoryError::Invalid(
                "native clock exceeds requested endpoint",
            ));
        }
        if accepted_steps >= controls.maximum_accepted_steps {
            return Err(TrajectoryError::Budget("accepted native segments"));
        }
        let remaining_s = finite(target_s - time_s)?;
        let mut h = remaining_s.min(maximum_step_s);
        let mut refinements = 0;
        loop {
            poll(cancelled)?;
            if h < controls.minimum_step_s && h != remaining_s {
                return Err(TrajectoryError::NoProgress {
                    time_s,
                    trial_step_s: h,
                });
            }
            let mut trial_target_s = if h == remaining_s {
                target_s
            } else {
                finite(time_s + h)?
            };
            if trial_target_s - time_s > maximum_step_s {
                trial_target_s = trial_target_s.next_down();
            }
            let dt_s = finite(trial_target_s - time_s)?;
            if dt_s <= 0.0 || (dt_s < controls.minimum_step_s && trial_target_s != target_s) {
                return Err(TrajectoryError::NoProgress {
                    time_s,
                    trial_step_s: dt_s,
                });
            }
            let fraction = finite((trial_target_s - start_s) / actual_duration_s)?;
            let mut stages = Vec::with_capacity(current.len());
            let mut reverse = false;
            for (body, motion) in current.iter().enumerate() {
                poll(cancelled)?;
                let offsets = offsets_at(
                    &initial_offsets[body],
                    loads[body].endpoint_prescribed_jumps_m,
                    fraction,
                    cancelled,
                )?;
                let reserve = work.reserve(motion, controls)?;
                let trial = motion.advance_to_time_controlled(
                    loads[body].mean_force_n,
                    &offsets,
                    trial_target_s,
                    cancelled,
                );
                match trial {
                    Ok(step) => {
                        work.refund(reserve, &step)?;
                        stages.push(step);
                    }
                    Err(DynamicsError::Model(CohesiveError::UnresolvedFrictionReversal {
                        ..
                    })) => {
                        reverse = true;
                        break;
                    }
                    Err(error) => return Err(TrajectoryError::Body { body, error }),
                }
            }
            if reverse {
                rejected_trials += 1;
                refinements += 1;
                if refinements > controls.maximum_refinements_per_step {
                    return Err(TrajectoryError::Refinement {
                        time_s,
                        attempts: refinements,
                    });
                }
                h = finite(0.5 * dt_s)?;
                continue;
            }
            let positions: Vec<_> = stages
                .iter()
                .map(|s| s.motion.state().displacements_m())
                .collect();
            let velocities: Vec<_> = stages.iter().map(|s| s.motion.velocity_m_s()).collect();
            let endpoint = geometry::sample(
                contacts,
                &reference_positions,
                &positions,
                &velocities,
                trial_target_s - start_s,
                cancelled,
            )?;
            if endpoint.gaps_m.iter().any(|gap| *gap < 0.0) {
                rejected_trials += 1;
                refinements += 1;
                if refinements > controls.maximum_refinements_per_step {
                    return Err(TrajectoryError::Refinement {
                        time_s,
                        attempts: refinements,
                    });
                }
                h = geometry::shorter_step(contacts, &sampled.gaps_m, &endpoint.gaps_m, dt_s)?;
                continue;
            }
            // Both endpoints are actual native states. Every declared contact
            // gap is affine on this one native straight displacement path.
            for step in &stages {
                balance.add_smooth(step)?;
            }
            minimum_gap_m = minimum_gap_m.min(endpoint.minimum_gap_m);
            current = stages.into_iter().map(|step| step.motion).collect();
            sampled = endpoint;
            accepted_steps += 1;
            break;
        }
    }
    let report = balance.finish(
        TrajectoryReport {
            elapsed_s: actual_duration_s,
            accepted_steps,
            body_trials: work.body_trials,
            rejected_trials,
            residual_evaluations: work.residual_evaluations,
            newton_iterations: work.newton_iterations,
            contact_ids: contacts.iter().map(|c| c.id).collect(),
            events,
            initial_gaps_m,
            final_gaps_m: sampled.gaps_m,
            minimum_gap_m,
            external_work_j: 0.0,
            prescribed_jump_work_j: 0.0,
            contact_prescribed_work_j: 0.0,
            kinetic_change_j: 0.0,
            bulk_stored_change_j: 0.0,
            cohesive_stored_change_j: 0.0,
            fracture_loss_j: 0.0,
            viscous_loss_j: 0.0,
            internal_friction_loss_j: 0.0,
            impact_loss_j: 0.0,
            energy_residual_j: 0.0,
            absolute_smooth_energy_defect_j: 0.0,
            absolute_impact_energy_defect_j: 0.0,
        },
        controls.energy_tolerance_j,
    )?;
    poll(cancelled)?;
    Ok(TrajectoryCandidate {
        motions: current,
        report,
    })
}

fn sample(
    contacts: &[FrictionalContact<'_>],
    reference_positions: &[&[f64]],
    motions: &[CohesiveMotion<'_, '_>],
    elapsed_s: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<geometry::Sample, TrajectoryError> {
    let positions: Vec<_> = motions
        .iter()
        .map(|m| m.state().displacements_m())
        .collect();
    let velocities: Vec<_> = motions.iter().map(CohesiveMotion::velocity_m_s).collect();
    Ok(geometry::sample(
        contacts,
        reference_positions,
        &positions,
        &velocities,
        elapsed_s,
        cancelled,
    )?)
}

#[allow(clippy::too_many_arguments)]
fn admit(
    motions: &[&CohesiveMotion<'_, '_>],
    contacts: &[FrictionalContact<'_>],
    loads: &[BodyLoad<'_>],
    preconditioners: &[&dyn Precond],
    duration_s: f64,
    controls: TrajectoryControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), TrajectoryError> {
    poll(cancelled)?;
    impact::admit_motions(
        motions,
        preconditioners,
        controls.impact.response.maximum_bodies,
        controls.maximum_workspace_entries,
        cancelled,
    )?;
    if loads.len() != motions.len()
        || contacts.len() > controls.impact.response.contact_solve.maximum_contacts
        || !positive(duration_s)
        || !positive(controls.maximum_step_s)
        || !positive(controls.minimum_step_s)
        || controls.minimum_step_s > controls.maximum_step_s
        || !positive(controls.energy_tolerance_j)
        || !positive(controls.closing_velocity_tolerance_m_s)
        || controls.closing_velocity_tolerance_m_s < controls.impact.response.velocity_tolerance_m_s
        || controls.closing_velocity_tolerance_m_s
            < controls.impact.response.contact_solve.velocity_tolerance
        || controls.maximum_accepted_steps == 0
        || controls.maximum_body_trials == 0
        || controls.maximum_refinements_per_step == 0
        || controls.maximum_residual_evaluations == 0
        || controls.maximum_newton_iterations == 0
    {
        return Err(TrajectoryError::Invalid(
            "finite bounded trajectory/load controls required",
        ));
    }
    let mut workspace = checked_mul(motions.len(), 32)?;
    workspace = checked_add(workspace, checked_mul(contacts.len(), 32)?)?;
    workspace = checked_add(
        workspace,
        checked_mul(
            controls.maximum_impacts,
            controls.impact.response.maximum_workspace_entries,
        )?,
    )?;
    for (motion, load) in motions.iter().zip(loads) {
        poll(cancelled)?;
        let n = motion.velocity_m_s().len();
        let patches = motion.state().points().len();
        if load.mean_force_n.len() != n || load.endpoint_prescribed_jumps_m.len() != patches {
            return Err(TrajectoryError::Invalid(
                "original body load and prescribed-jump dimensions",
            ));
        }
        workspace = checked_add(workspace, checked_mul(n, 8)?)?;
        workspace = checked_add(workspace, checked_mul(patches, 64)?)?;
        if workspace > controls.maximum_workspace_entries {
            return Err(TrajectoryError::Budget(
                "trajectory buffers and retained contact reports",
            ));
        }
        for chunk in load.mean_force_n.chunks(256) {
            poll(cancelled)?;
            for &value in chunk {
                finite(value)?;
            }
        }
        for values in load.endpoint_prescribed_jumps_m {
            poll(cancelled)?;
            for &value in values {
                finite(value)?;
            }
        }
    }
    Ok(())
}

fn offsets_at(
    initial: &[[f64; 3]],
    target: &[[f64; 3]],
    fraction: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<[f64; 3]>, TrajectoryError> {
    let mut result = Vec::with_capacity(initial.len());
    for (a, b) in initial.iter().zip(target) {
        poll(cancelled)?;
        let mut value = [0.0; 3];
        for d in 0..3 {
            value[d] = if fraction == 1.0 {
                b[d]
            } else {
                finite(a[d] + fraction * finite(b[d] - a[d])?)?
            };
        }
        result.push(value);
    }
    Ok(result)
}

#[derive(Default)]
struct Work {
    body_trials: usize,
    residual_evaluations: usize,
    newton_iterations: usize,
}
impl Work {
    fn reserve(
        &mut self,
        motion: &CohesiveMotion<'_, '_>,
        controls: TrajectoryControls,
    ) -> Result<(usize, usize), TrajectoryError> {
        let config = motion.dynamics_config().solve;
        let body_trials = checked_add(self.body_trials, 1)?;
        let residuals = checked_add(self.residual_evaluations, config.max_operator_evaluations)?;
        let iterations = checked_add(self.newton_iterations, config.max_newton_iterations)?;
        if body_trials > controls.maximum_body_trials
            || residuals > controls.maximum_residual_evaluations
            || iterations > controls.maximum_newton_iterations
        {
            return Err(TrajectoryError::Budget("original body trial work"));
        }
        self.body_trials = body_trials;
        self.residual_evaluations = residuals;
        self.newton_iterations = iterations;
        Ok((
            config.max_operator_evaluations,
            config.max_newton_iterations,
        ))
    }
    fn refund(
        &mut self,
        reserve: (usize, usize),
        step: &DynamicsStep<'_, '_>,
    ) -> Result<(), TrajectoryError> {
        if step.residual_evaluations > reserve.0 || step.nonlinear.iterations > reserve.1 {
            return Err(TrajectoryError::Invalid(
                "native work exceeded reserved allowance",
            ));
        }
        self.residual_evaluations -= reserve.0 - step.residual_evaluations;
        self.newton_iterations -= reserve.1 - step.nonlinear.iterations;
        Ok(())
    }
}

// Per-quantity compensated sums retain actual accepted increments and keep
// rejected-trial losses entirely outside the physical trajectory account.
#[derive(Default)]
struct Balance {
    external: Sum,
    prescribed: Sum,
    contact_work: Sum,
    kinetic: Sum,
    bulk: Sum,
    cohesive: Sum,
    fracture: Sum,
    viscous: Sum,
    friction: Sum,
    impact: Sum,
    smooth_defect: Sum,
    impact_defect: Sum,
    closure: Sum,
}
impl Balance {
    fn add_smooth(&mut self, step: &DynamicsStep<'_, '_>) -> Result<(), TrajectoryError> {
        for value in [
            step.external_work_j,
            step.prescribed_jump_work_j,
            -(step.kinetic_after_j - step.kinetic_before_j),
            -(step.bulk_stored_after_j - step.bulk_stored_before_j),
            -step.cohesive_energy_j.stored_change,
            -step.cohesive_energy_j.dissipated,
            -step.viscous_dissipated_j,
            -step.friction_dissipated_j,
        ] {
            self.closure.add(value)?;
        }
        self.external.add(step.external_work_j)?;
        self.prescribed.add(step.prescribed_jump_work_j)?;
        self.kinetic
            .add(finite(step.kinetic_after_j - step.kinetic_before_j)?)?;
        self.bulk.add(finite(
            step.bulk_stored_after_j - step.bulk_stored_before_j,
        )?)?;
        self.cohesive.add(step.cohesive_energy_j.stored_change)?;
        self.fracture.add(step.cohesive_energy_j.dissipated)?;
        self.viscous.add(step.viscous_dissipated_j)?;
        self.friction.add(step.friction_dissipated_j)?;
        self.smooth_defect.add(step.energy_residual_j.abs())
    }
    fn add_impact(
        &mut self,
        r: &super::impact::frictional::FrictionalResponseReport,
    ) -> Result<(), TrajectoryError> {
        for value in [
            r.prescribed_work_j,
            -r.kinetic_change_j,
            -r.contact.dissipated_energy,
        ] {
            self.closure.add(value)?;
        }
        self.kinetic.add(r.kinetic_change_j)?;
        self.contact_work.add(r.prescribed_work_j)?;
        self.impact.add(r.contact.dissipated_energy)?;
        self.impact_defect.add(r.work_energy_residual_j.abs())?;
        self.impact_defect.add(r.dissipation_residual_j.abs())?;
        self.impact_defect.add(r.absolute_body_work_residual_j)
    }
    fn finish(
        self,
        mut r: TrajectoryReport,
        tolerance_j: f64,
    ) -> Result<TrajectoryReport, TrajectoryError> {
        r.external_work_j = self.external.total()?;
        r.prescribed_jump_work_j = self.prescribed.total()?;
        r.contact_prescribed_work_j = self.contact_work.total()?;
        r.kinetic_change_j = self.kinetic.total()?;
        r.bulk_stored_change_j = self.bulk.total()?;
        r.cohesive_stored_change_j = self.cohesive.total()?;
        r.fracture_loss_j = self.fracture.total()?;
        r.viscous_loss_j = self.viscous.total()?;
        r.internal_friction_loss_j = self.friction.total()?;
        r.impact_loss_j = self.impact.total()?;
        r.absolute_smooth_energy_defect_j = self.smooth_defect.total()?;
        r.absolute_impact_energy_defect_j = self.impact_defect.total()?;
        // Retain every original increment in one compensated account. Summing
        // rounded channel totals here could discard a small physical defect
        // beside large opposing work/storage contributions.
        r.energy_residual_j = self.closure.total()?;
        for (quantity, residual_j) in [
            ("full trajectory first law", r.energy_residual_j.abs()),
            (
                "absolute native smooth work defects",
                r.absolute_smooth_energy_defect_j,
            ),
            (
                "absolute native impact work defects",
                r.absolute_impact_energy_defect_j,
            ),
        ] {
            if residual_j > tolerance_j {
                return Err(TrajectoryError::Energy {
                    quantity,
                    residual_j,
                    tolerance_j,
                });
            }
        }
        Ok(r)
    }
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), TrajectoryError> {
        finite(value)?;
        let next = finite(self.value + value)?;
        let residual = if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.correction = finite(self.correction + residual)?;
        self.value = next;
        Ok(())
    }
    fn total(self) -> Result<f64, TrajectoryError> {
        finite(self.value + self.correction)
    }
}
fn finite(value: f64) -> Result<f64, TrajectoryError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(TrajectoryError::NonFinite("derived trajectory value"))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn checked_add(a: usize, b: usize) -> Result<usize, TrajectoryError> {
    a.checked_add(b)
        .ok_or(TrajectoryError::Budget("count/workspace overflow"))
}
fn checked_mul(a: usize, b: usize) -> Result<usize, TrajectoryError> {
    a.checked_mul(b)
        .ok_or(TrajectoryError::Budget("count/workspace overflow"))
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), TrajectoryError> {
    if cancelled() {
        Err(TrajectoryError::Cancelled)
    } else {
        Ok(())
    }
}
