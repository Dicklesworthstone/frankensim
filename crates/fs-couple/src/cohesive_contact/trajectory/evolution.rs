//! Complete fixed-frame trajectories with free motion, impacts and support.
//!
//! This owner shares the original transient loop. Closed stationary contacts
//! enter the original sustained force solve; a bracketed sliding stop ends one
//! native segment and the next force solve decides sticking or breakaway. No
//! body velocity or gap is projected. Open contacts remain in geometry checks.
//! Arrivals narrow activation to a separately requested gap accuracy and retain
//! the original response at the accepted endpoint of the arrival bracket.
//!
//! Geometry and straight native displacement paths have the same approximation
//! as the transient owner. A positive impact activation gap is not silently
//! treated as zero for support: the unchanged physical work/kinematic audit may
//! refuse that handoff. This does not certify exact-flow event times, arbitrary
//! nonlinear contact uniqueness, or a general Moreau--Jean discretization.

pub use super::arrival::{ArrivalControls, ArrivalReport};
use super::sustained::{
    NaturalMapScales, SegmentCandidate, SegmentControls, StoppingReport, SustainedContact,
    segment_candidate_controlled,
};
use super::{
    BodyLoad, FrictionalContact, Sum, TrajectoryControls, TrajectoryError, TrajectoryReport, Work,
    advance_shared, checked_add, checked_mul, finite, geometry, offsets_at, poll, positive,
};
use fs_exec::Cx;
use fs_solid::cohesive::dynamics::CohesiveMotion;
use fs_sparse::precond::Precond;

/// Global trajectory limits and the original per-segment solver controls.
#[derive(Clone, Copy, Debug)]
pub struct EvolutionControls {
    /// Native free/impact limits, including the whole-trajectory energy gate.
    /// Workspace includes the supported adapter and contact solver reservations.
    pub trajectory: TrajectoryControls,
    pub sustained: SegmentControls,
    /// Arrival search accuracy. Impacts and support activate only inside the
    /// smaller of this gap tolerance and each authored activation band. Zero
    /// requests an actual zero gap; it never permits position projection.
    pub arrival: ArrivalControls,
    /// Across all accepted and geometrically rejected sustained segments.
    pub maximum_contact_solves: usize,
    pub maximum_response_evaluations: usize,
    /// Retained stop events, each with actual simultaneous stopped contacts.
    pub maximum_stopping_events: usize,
}

/// One actual native endpoint with vanishing relative tangential velocity.
#[derive(Debug)]
pub struct TimedStoppingReport {
    pub time_s: f64,
    pub segment_start_time_s: f64,
    pub contact_ids: Vec<u64>,
    /// Bracket coordinates are durations from `segment_start_time_s`.
    pub stopping: StoppingReport,
}

/// One original native endpoint accepted at a contact arrival. The absolute
/// bracket includes outward rounding of the relative bracket's clock offset.
#[derive(Clone, Copy, Debug)]
pub struct TimedArrivalReport {
    pub contact_id: u64,
    pub time_s: f64,
    pub segment_start_time_s: f64,
    pub bracket_time_s: [f64; 2],
    pub arrival: ArrivalReport,
}

/// The complete physical window. External Coulomb loss is separate from the
/// original material friction and instantaneous impact loss in `trajectory`.
#[derive(Debug)]
pub struct EvolutionReport {
    /// Its external work contains only authored body loads; contact-prescribed
    /// work includes impact and sustained moving-surface work. Its energy
    /// residual also includes this report's `contact_friction_loss_j`.
    pub trajectory: TrajectoryReport,
    pub sustained_steps: usize,
    pub stopping_events: Vec<TimedStoppingReport>,
    /// At most one arrival is retained per accepted positive native segment.
    /// All unpublished arrival probes still charge `trajectory`'s body work.
    pub arrival_events: Vec<TimedArrivalReport>,
    pub contact_friction_loss_j: f64,
    pub absolute_sustained_energy_defect_j: f64,
    pub absolute_sustained_native_energy_defect_j: f64,
    pub absolute_contact_work_defect_j: f64,
    pub absolute_friction_work_defect_j: f64,
    pub absolute_normal_constraint_work_j: f64,
    pub contact_solves: usize,
    pub response_evaluations: usize,
    /// Original contact order; numerical guesses for a subsequent solve, not
    /// a claim of instantaneous force at a final impact or open contact.
    pub warm_start_forces_n: Vec<[f64; 3]>,
}

pub struct EvolutionCandidate<'m, 'a> {
    pub motions: Vec<CohesiveMotion<'m, 'a>>,
    pub report: EvolutionReport,
}

/// Advance the entire requested interval through original native owners.
/// Scales and warm-force guesses follow the complete original contact order.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    contacts: &[FrictionalContact<'_>],
    loads: &[BodyLoad<'_>],
    preconditioners: &[&dyn Precond],
    scales: &[NaturalMapScales],
    initial_forces_n: &[[f64; 3]],
    duration_s: f64,
    controls: EvolutionControls,
    cx: &Cx<'_>,
) -> Result<EvolutionCandidate<'m, 'a>, TrajectoryError> {
    advance_candidate_controlled(
        motions,
        contacts,
        loads,
        preconditioners,
        scales,
        initial_forces_n,
        duration_s,
        controls,
        &mut || cx.checkpoint().is_err(),
    )
}

/// Predicate counterpart. Every candidate remains private until all native
/// states, geometry, contact laws and whole-window energy gates pass.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate_controlled<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    contacts: &[FrictionalContact<'_>],
    loads: &[BodyLoad<'_>],
    preconditioners: &[&dyn Precond],
    scales: &[NaturalMapScales],
    initial_forces_n: &[[f64; 3]],
    duration_s: f64,
    controls: EvolutionControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<EvolutionCandidate<'m, 'a>, TrajectoryError> {
    poll(cancelled)?;
    let mut state = State::new(motions, contacts, scales, initial_forces_n, controls)?;
    let candidate = advance_shared(
        motions,
        contacts,
        loads,
        preconditioners,
        duration_s,
        state.base_controls,
        Some(&mut state),
        cancelled,
    )?;
    let report = state.finish(candidate.report)?;
    poll(cancelled)?;
    Ok(EvolutionCandidate {
        motions: candidate.motions,
        report,
    })
}

pub(super) struct State<'s> {
    controls: EvolutionControls,
    base_controls: TrajectoryControls,
    scales: &'s [NaturalMapScales],
    warm_forces: Vec<[f64; 3]>,
    sustained_steps: usize,
    stopping_events: Vec<TimedStoppingReport>,
    arrival_events: Vec<TimedArrivalReport>,
    contact_solves: usize,
    response_evaluations: usize,
    friction_loss: Sum,
    energy_defect: Sum,
    native_defect: Sum,
    contact_defect: Sum,
    friction_defect: Sum,
    normal_defect: Sum,
}

impl<'s> State<'s> {
    fn new(
        motions: &[&CohesiveMotion<'_, '_>],
        contacts: &[FrictionalContact<'_>],
        scales: &'s [NaturalMapScales],
        initial_forces_n: &[[f64; 3]],
        controls: EvolutionControls,
    ) -> Result<Self, TrajectoryError> {
        if contacts.len() > 21
            || scales.len() != contacts.len()
            || initial_forces_n.len() != contacts.len()
            || controls.maximum_contact_solves == 0
            || controls.maximum_response_evaluations == 0
            || !positive(controls.sustained.step.audit.velocity_tolerance_m_s)
            || controls.trajectory.closing_velocity_tolerance_m_s
                > controls.sustained.step.audit.velocity_tolerance_m_s
            || controls.sustained.stopping.minimum_step_s < controls.trajectory.minimum_step_s
            || !positive(controls.arrival.minimum_step_s)
            || !positive(controls.arrival.time_tolerance_s)
            || !controls.arrival.gap_tolerance_m.is_finite()
            || controls.arrival.gap_tolerance_m < 0.0
            || controls.arrival.maximum_trials == 0
            || controls.arrival.minimum_step_s < controls.trajectory.minimum_step_s
        {
            return Err(TrajectoryError::Invalid(
                "complete supported controls and compatible normal/progress tolerances required",
            ));
        }
        for scale in scales {
            if !positive(scale.normal_n_per_m) || !positive(scale.tangent_n_per_m) {
                return Err(TrajectoryError::Invalid(
                    "positive finite contact scales required",
                ));
            }
        }
        if initial_forces_n
            .iter()
            .flatten()
            .any(|f| !f.is_finite() || f.abs() > controls.sustained.step.maximum_contact_force_n)
        {
            return Err(TrajectoryError::Invalid(
                "finite bounded force warm starts required",
            ));
        }
        // The original transient adapter keeps its own workspace allowance.
        // Reserve the simultaneously live sustained buffers and contact-solver
        // scratch separately; original body/material scratch stays owner-bound.
        // Includes the refined impact-contact view kept by the shared loop.
        let mut extra = checked_add(256, checked_mul(contacts.len(), 48)?)?;
        extra = checked_add(extra, checked_mul(motions.len(), 16)?)?;
        for motion in motions {
            extra = checked_add(extra, checked_mul(motion.state().points().len(), 3)?)?;
        }
        let event_header =
            core::mem::size_of::<TimedStoppingReport>().div_ceil(core::mem::size_of::<f64>());
        let event_entries = checked_add(event_header, checked_mul(contacts.len(), 2)?.max(4))?;
        extra = checked_add(
            extra,
            checked_mul(controls.maximum_stopping_events, event_entries)?,
        )?;
        let arrival_entries =
            core::mem::size_of::<TimedArrivalReport>().div_ceil(core::mem::size_of::<f64>());
        extra = checked_add(
            extra,
            checked_mul(controls.trajectory.maximum_accepted_steps, arrival_entries)?,
        )?;
        if !contacts.is_empty() {
            extra = checked_add(extra, controls.sustained.step.maximum_workspace_entries)?;
            extra = checked_add(
                extra,
                controls.sustained.step.solve.maximum_workspace_entries,
            )?;
        }
        let mut base_controls = controls.trajectory;
        base_controls.maximum_workspace_entries = base_controls
            .maximum_workspace_entries
            .checked_sub(extra)
            .ok_or(TrajectoryError::Budget("supported trajectory workspace"))?;
        Ok(Self {
            controls,
            base_controls,
            scales,
            warm_forces: initial_forces_n.to_vec(),
            sustained_steps: 0,
            stopping_events: Vec::with_capacity(controls.maximum_stopping_events),
            arrival_events: Vec::with_capacity(controls.trajectory.maximum_accepted_steps),
            contact_solves: 0,
            response_evaluations: 0,
            friction_loss: Sum::default(),
            energy_defect: Sum::default(),
            native_defect: Sum::default(),
            contact_defect: Sum::default(),
            friction_defect: Sum::default(),
            normal_defect: Sum::default(),
        })
    }

    pub(super) fn after_impact(&mut self) {
        // The event changed velocity; a previously sustained force is only a
        // guess, so restart its continuation from zero after this mode change.
        self.warm_forces.fill([0.0; 3]);
    }

    pub(super) fn arrival_controls(&self) -> ArrivalControls {
        self.controls.arrival
    }

    /// Called once after a complete positive native segment and every contact's
    /// geometry have passed. Scalar probes are already charged through the
    /// shared trajectory Work; retaining metadata does not adopt another state.
    pub(super) fn accept_arrival(
        &mut self,
        contact_id: u64,
        segment_start_time_s: f64,
        arrival: ArrivalReport,
    ) -> Result<(), TrajectoryError> {
        if self.arrival_events.len() >= self.controls.trajectory.maximum_accepted_steps {
            return Err(TrajectoryError::Budget("retained contact-arrival events"));
        }
        if !segment_start_time_s.is_finite()
            || segment_start_time_s < 0.0
            || !positive(arrival.elapsed_s)
            || !arrival.bracket_s[0].is_finite()
            || !arrival.bracket_s[1].is_finite()
            || arrival.bracket_s[0] < 0.0
            || arrival.bracket_s[0] > arrival.elapsed_s
            || arrival.elapsed_s > arrival.bracket_s[1]
        {
            return Err(TrajectoryError::Invalid(
                "complete finite arrival bracket required",
            ));
        }
        let time_s = finite(segment_start_time_s + arrival.elapsed_s)?;
        if time_s <= segment_start_time_s {
            return Err(TrajectoryError::NoProgress {
                time_s: segment_start_time_s,
                trial_step_s: arrival.elapsed_s,
            });
        }
        let bracket_time_s = if arrival.bracket_s[0] == arrival.bracket_s[1] {
            // An exact numerical endpoint is a point in the represented clock.
            [time_s; 2]
        } else {
            let low = if arrival.bracket_s[0] == 0.0 {
                segment_start_time_s
            } else {
                finite(finite(segment_start_time_s + arrival.bracket_s[0])?.next_down())?
            };
            let high = finite(finite(segment_start_time_s + arrival.bracket_s[1])?.next_up())?;
            [low, high]
        };
        self.arrival_events.push(TimedArrivalReport {
            contact_id,
            time_s,
            segment_start_time_s,
            bracket_time_s,
            arrival,
        });
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn trial<'m, 'a>(
        &mut self,
        motions: &[CohesiveMotion<'m, 'a>],
        contacts: &[FrictionalContact<'_>],
        sampled: &geometry::Sample,
        loads: &[BodyLoad<'_>],
        initial_offsets: &[Vec<[f64; 3]>],
        start_s: f64,
        requested_duration_s: f64,
        target_s: f64,
        trajectory: TrajectoryControls,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<SegmentCandidate<'m, 'a>>, TrajectoryError> {
        poll(cancelled)?;
        let mut active = Vec::with_capacity(contacts.len());
        let mut scales = Vec::with_capacity(contacts.len());
        let mut initial = Vec::with_capacity(contacts.len());
        for (i, contact) in contacts.iter().enumerate() {
            poll(cancelled)?;
            let activation_tolerance_m = contact
                .activation_tolerance_m
                .min(self.controls.arrival.gap_tolerance_m);
            if sampled.gaps_m[i] <= activation_tolerance_m
                && sampled.normal_velocities_m_s[i].abs()
                    <= trajectory.closing_velocity_tolerance_m_s
            {
                active.push(SustainedContact {
                    id: contact.id,
                    gap_m: sampled.gaps_m[i],
                    activation_tolerance_m,
                    friction_coefficient: contact.friction_coefficient,
                    prescribed_velocity_m_s: contact.prescribed_velocity_m_s,
                    trace: contact.trace,
                });
                scales.push(self.scales[i]);
                initial.push(self.warm_forces[i]);
            }
        }
        if active.is_empty() {
            self.warm_forces.fill([0.0; 3]);
            return Ok(None);
        }
        let fraction = finite((target_s - start_s) / requested_duration_s)?;
        let mut offsets = Vec::with_capacity(motions.len());
        for (body, load) in loads.iter().enumerate() {
            offsets.push(offsets_at(
                &initial_offsets[body],
                load.endpoint_prescribed_jumps_m,
                fraction,
                cancelled,
            )?);
        }
        let local_loads: Vec<_> = loads
            .iter()
            .zip(&offsets)
            .map(|(load, offsets)| BodyLoad {
                mean_force_n: load.mean_force_n,
                endpoint_prescribed_jumps_m: offsets,
            })
            .collect();
        let borrowed: Vec<_> = motions.iter().collect();
        let mut controls = self.controls.sustained;
        controls.step.maximum_body_trials = controls
            .step
            .maximum_body_trials
            .min(trajectory.maximum_body_trials - work.body_trials);
        controls.step.maximum_residual_evaluations = controls
            .step
            .maximum_residual_evaluations
            .min(trajectory.maximum_residual_evaluations - work.residual_evaluations);
        controls.step.maximum_native_newton_iterations = controls
            .step
            .maximum_native_newton_iterations
            .min(trajectory.maximum_newton_iterations - work.newton_iterations);
        controls.maximum_contact_solves = controls
            .maximum_contact_solves
            .min(self.controls.maximum_contact_solves - self.contact_solves);
        controls.maximum_response_evaluations = controls
            .maximum_response_evaluations
            .min(self.controls.maximum_response_evaluations - self.response_evaluations);
        if controls.step.maximum_body_trials == 0
            || controls.step.maximum_residual_evaluations == 0
            || controls.step.maximum_native_newton_iterations == 0
            || controls.maximum_contact_solves == 0
            || controls.maximum_response_evaluations == 0
        {
            return Err(TrajectoryError::Budget("whole-trajectory supported work"));
        }
        // A refused nested solve cannot publish a partial work report. Reserve
        // its complete bounded allowance before entering it; successful solves
        // refund only the unused part, including all stopping-time probes.
        work.body_trials = checked_add(work.body_trials, controls.step.maximum_body_trials)?;
        work.residual_evaluations = checked_add(
            work.residual_evaluations,
            controls.step.maximum_residual_evaluations,
        )?;
        work.newton_iterations = checked_add(
            work.newton_iterations,
            controls.step.maximum_native_newton_iterations,
        )?;
        self.contact_solves = checked_add(self.contact_solves, controls.maximum_contact_solves)?;
        self.response_evaluations = checked_add(
            self.response_evaluations,
            controls.maximum_response_evaluations,
        )?;
        let candidate = segment_candidate_controlled(
            &borrowed,
            &active,
            &local_loads,
            &scales,
            &initial,
            finite(target_s - motions[0].time_s())?,
            controls,
            cancelled,
        )
        .map_err(TrajectoryError::Sustained)?;
        let step = &candidate.step.report;
        if step.body_trials > controls.step.maximum_body_trials
            || step.residual_evaluations > controls.step.maximum_residual_evaluations
            || step.native_newton_iterations > controls.step.maximum_native_newton_iterations
            || candidate.report.contact_solves > controls.maximum_contact_solves
            || candidate.report.response_evaluations > controls.maximum_response_evaluations
        {
            return Err(TrajectoryError::Invalid(
                "supported work exceeded reserved allowance",
            ));
        }
        work.body_trials -= controls.step.maximum_body_trials - step.body_trials;
        work.residual_evaluations -=
            controls.step.maximum_residual_evaluations - step.residual_evaluations;
        work.newton_iterations -=
            controls.step.maximum_native_newton_iterations - step.native_newton_iterations;
        self.contact_solves -= controls.maximum_contact_solves - candidate.report.contact_solves;
        self.response_evaluations -=
            controls.maximum_response_evaluations - candidate.report.response_evaluations;
        Ok(Some(candidate))
    }

    pub(super) fn accept(
        &mut self,
        start_time_s: f64,
        candidate: &mut SegmentCandidate<'_, '_>,
        contacts: &[FrictionalContact<'_>],
    ) -> Result<(), TrajectoryError> {
        let step = &candidate.step.report;
        if step.contact_ids.len() != step.contact.forces_n.len() {
            return Err(TrajectoryError::Invalid(
                "complete admitted contact forces required",
            ));
        }
        self.warm_forces.fill([0.0; 3]);
        for (&id, &force) in step.contact_ids.iter().zip(&step.contact.forces_n) {
            let index = contacts
                .iter()
                .position(|contact| contact.id == id)
                .ok_or(TrajectoryError::Invalid("unknown admitted contact ID"))?;
            self.warm_forces[index] = force;
        }
        if let Some(stopping) = candidate.report.stopping.take() {
            if self.stopping_events.len() >= self.controls.maximum_stopping_events {
                return Err(TrajectoryError::Budget("retained sliding-stop events"));
            }
            if candidate.report.stopped_contact_ids.capacity()
                > checked_mul(contacts.len(), 2)?.max(4)
            {
                return Err(TrajectoryError::Budget("retained stopped-contact IDs"));
            }
            self.stopping_events.push(TimedStoppingReport {
                time_s: step.target_time_s,
                segment_start_time_s: start_time_s,
                contact_ids: core::mem::take(&mut candidate.report.stopped_contact_ids),
                stopping,
            });
        }
        self.sustained_steps = checked_add(self.sustained_steps, 1)?;
        self.friction_loss.add(step.audit.friction_loss_j)?;
        self.energy_defect.add(step.audit.energy_residual_j.abs())?;
        self.native_defect
            .add(step.audit.absolute_native_energy_defect_j)?;
        self.contact_defect
            .add(step.audit.contact_work_residual_j.abs())?;
        self.friction_defect
            .add(step.audit.absolute_friction_work_defect_j)?;
        self.normal_defect
            .add(step.audit.absolute_normal_constraint_work_j)
    }

    fn finish(self, trajectory: TrajectoryReport) -> Result<EvolutionReport, TrajectoryError> {
        let report = EvolutionReport {
            trajectory,
            sustained_steps: self.sustained_steps,
            stopping_events: self.stopping_events,
            arrival_events: self.arrival_events,
            contact_friction_loss_j: self.friction_loss.total()?,
            absolute_sustained_energy_defect_j: self.energy_defect.total()?,
            absolute_sustained_native_energy_defect_j: self.native_defect.total()?,
            absolute_contact_work_defect_j: self.contact_defect.total()?,
            absolute_friction_work_defect_j: self.friction_defect.total()?,
            absolute_normal_constraint_work_j: self.normal_defect.total()?,
            contact_solves: self.contact_solves,
            response_evaluations: self.response_evaluations,
            warm_start_forces_n: self.warm_forces,
        };
        for (quantity, residual_j) in [
            (
                "absolute sustained first-law defects",
                report.absolute_sustained_energy_defect_j,
            ),
            (
                "absolute sustained native defects",
                report.absolute_sustained_native_energy_defect_j,
            ),
            (
                "absolute sustained contact-work defects",
                report.absolute_contact_work_defect_j,
            ),
            (
                "absolute sustained friction-work defects",
                report.absolute_friction_work_defect_j,
            ),
            (
                "absolute sustained normal constraint work",
                report.absolute_normal_constraint_work_j,
            ),
        ] {
            if residual_j > self.controls.trajectory.energy_tolerance_j {
                return Err(TrajectoryError::Energy {
                    quantity,
                    residual_j,
                    tolerance_j: self.controls.trajectory.energy_tolerance_j,
                });
            }
        }
        Ok(report)
    }
}

#[cfg(test)]
#[path = "evolution/tests.rs"]
mod tests;
