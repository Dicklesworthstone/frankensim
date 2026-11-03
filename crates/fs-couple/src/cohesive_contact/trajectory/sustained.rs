//! Held contact forces coupled to original cohesive midpoint evolution.
//!
//! Each contact-force trial advances every original body under F + B^T lambda.
//! Original mass, stiffness, material path means, damage, viscosity and internal
//! friction therefore participate in the response. Only the contact law is
//! linearized analytically; its original Newton owner differences the actual
//! mechanical response. No frozen mass mobility replaces finite-step physics.
//!
//! This is a local continuation candidate for already active, normal-stationary
//! contacts with fixed affine frames. It inherits the native nonlinear material
//! branch and does not certify global contact uniqueness or arbitrary Painleve
//! regimes. Incoming impacts, new contact geometry and unresolved slip reversals
//! require their own event/refinement owner. No accepted state is projected.

pub use fs_contact::sustained::NaturalMapScales;
use fs_contact::sustained::solve::{SolveControls, SolveError, SolveReport, solve};
use fs_exec::Cx;
use fs_solid::cohesive::dynamics::{CohesiveMotion, DynamicsError, DynamicsStep};

use super::super::impact::{self, ImpactError, TraceTerm};
use super::BodyLoad;
use super::geometry::{self, GeometryError, Sum};

mod audit;
#[cfg(test)]
#[path = "sustained/tests.rs"]
mod tests;
pub use audit::{AuditControls, AuditError, AuditReport, BodyWork, ContactAuditInput};

/// A fixed orthonormal contact frame in original reduced body coordinates.
/// Positive normal displacement separates the sides. The initial gap belongs
/// to the supplied accepted states; prescribed relative velocity is constant.
#[derive(Clone, Copy, Debug)]
pub struct SustainedContact<'t> {
    pub id: u64,
    pub gap_m: f64,
    pub activation_tolerance_m: f64,
    pub friction_coefficient: f64,
    pub prescribed_velocity_m_s: [f64; 3],
    pub trace: [&'t [TraceTerm]; 3],
}

#[derive(Clone, Copy, Debug)]
pub struct SustainedControls {
    /// Original contact Newton/FGMRES and physical-response difference controls.
    pub solve: SolveControls,
    pub audit: AuditControls,
    pub maximum_bodies: usize,
    pub maximum_trace_terms: usize,
    pub maximum_body_trials: usize,
    /// Each native body trial reserves its entire residual allowance. Success
    /// refunds unused work; failure retains its conservative charge.
    pub maximum_residual_evaluations: usize,
    pub maximum_native_newton_iterations: usize,
    /// Absolute componentwise limit on both trial and accepted contact forces.
    pub maximum_contact_force_n: f64,
    /// Adapter buffers and one retained set of native body receipts. Native
    /// material/operator storage and solver scratch retain their own controls.
    pub maximum_workspace_entries: usize,
}

pub struct SustainedCandidate<'m, 'a> {
    /// One actual native step per original body, including its accepted motion,
    /// material losses, local tractions and nonlinear receipt. Its external
    /// force work includes the solved contact force; `report.body_work` splits
    /// authored work for the whole coupled physical audit.
    pub bodies: Vec<DynamicsStep<'m, 'a>>,
    pub report: SustainedReport,
}

#[derive(Debug)]
pub struct SustainedReport {
    pub elapsed_s: f64,
    pub target_time_s: f64,
    pub contact_ids: Vec<u64>,
    /// Original force-law solve, with complete force and path ordering.
    pub contact: SolveReport,
    pub audit: AuditReport,
    pub body_work: Vec<BodyWork>,
    pub body_trials: usize,
    /// Charged work includes full reservations for any refused native trials.
    pub residual_evaluations: usize,
    pub native_newton_iterations: usize,
}

/// Concrete refusal from the original mechanical response or its admission.
#[derive(Debug)]
pub enum ResponseError {
    Invalid(&'static str),
    Budget(&'static str),
    NonFinite(&'static str),
    Cancelled,
    Geometry(GeometryError),
    State(ImpactError),
    Body { body: usize, error: DynamicsError },
}
impl core::fmt::Display for ResponseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cohesive held-force response refused: {self:?}")
    }
}
impl std::error::Error for ResponseError {}
impl From<GeometryError> for ResponseError {
    fn from(error: GeometryError) -> Self {
        Self::Geometry(error)
    }
}
impl From<ImpactError> for ResponseError {
    fn from(error: ImpactError) -> Self {
        Self::State(error)
    }
}

#[derive(Debug)]
pub enum SustainedError {
    Response(ResponseError),
    Solve(SolveError<ResponseError>),
    Audit(AuditError),
}
impl core::fmt::Display for SustainedError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cohesive sustained contact refused: {self:?}")
    }
}
impl std::error::Error for SustainedError {}
impl From<ResponseError> for SustainedError {
    fn from(error: ResponseError) -> Self {
        Self::Response(error)
    }
}
impl From<GeometryError> for SustainedError {
    fn from(error: GeometryError) -> Self {
        Self::Response(ResponseError::Geometry(error))
    }
}

impl From<AuditError> for SustainedError {
    fn from(error: AuditError) -> Self {
        Self::Audit(error)
    }
}

/// Solve one complete finite interval of already active, normally stationary
/// contact. Forces are in newtons. Reusing the previous admitted force is a
/// warm start, not a replacement for solving the new actual body response.
/// The numerical N/m scales affect conditioning, never material coefficients.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    contacts: &[SustainedContact<'_>],
    loads: &[BodyLoad<'_>],
    scales: &[NaturalMapScales],
    initial_forces_n: &[[f64; 3]],
    duration_s: f64,
    controls: SustainedControls,
    cx: &Cx<'_>,
) -> Result<SustainedCandidate<'m, 'a>, SustainedError> {
    step_candidate_controlled(
        motions,
        contacts,
        loads,
        scales,
        initial_forces_n,
        duration_s,
        controls,
        &mut || cx.checkpoint().is_err(),
    )
}

#[allow(clippy::too_many_arguments)]
pub fn step_candidate_controlled<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    contacts: &[SustainedContact<'_>],
    loads: &[BodyLoad<'_>],
    scales: &[NaturalMapScales],
    initial_forces_n: &[[f64; 3]],
    duration_s: f64,
    controls: SustainedControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SustainedCandidate<'m, 'a>, SustainedError> {
    admit(
        motions,
        contacts,
        loads,
        scales,
        initial_forces_n,
        duration_s,
        controls,
        cancelled,
    )?;
    let start_s = motions[0].time_s();
    let target_s = finite(start_s + duration_s)?;
    let dt_s = finite(target_s - start_s)?;
    if dt_s <= 0.0 {
        return Err(ResponseError::Invalid("unrepresentable positive interval").into());
    }
    let positions: Vec<_> = motions
        .iter()
        .map(|m| m.state().displacements_m())
        .collect();
    let before_velocities: Vec<_> = motions.iter().map(|m| m.velocity_m_s()).collect();
    let before = sample_paths(
        contacts,
        &positions,
        &positions,
        &before_velocities,
        0.0,
        cancelled,
    )?;
    // An incoming normal impact must pass through the existing impact owner;
    // a finite position constraint must not silently turn it into restitution.
    for velocity in &before.velocities_m_s {
        if velocity[0].abs() > controls.audit.velocity_tolerance_m_s {
            return Err(ResponseError::Invalid(
                "sustained contact requires stationary relative normals",
            )
            .into());
        }
    }
    let mut work = Work::default();
    let mut latest: Option<Trial<'m, 'a>> = None;
    let mu: Vec<_> = contacts.iter().map(|c| c.friction_coefficient).collect();
    let mut response = |forces: &[[f64; 3]], poller: &mut dyn FnMut() -> bool| {
        // Only the final replay may publish; release the obsolete native
        // candidate before allocating the next complete material response.
        drop(latest.take());
        let trial = evaluate(
            motions,
            contacts,
            loads,
            &positions,
            forces,
            target_s,
            dt_s,
            controls,
            &mut work,
            &mut || poller(),
        )?;
        let paths = trial.sample.paths_m.clone();
        latest = Some(trial);
        Ok::<_, ResponseError>(paths)
    };
    let contact = solve(
        &mu,
        scales,
        initial_forces_n,
        controls.solve,
        &mut response,
        cancelled,
    )
    .map_err(SustainedError::Solve)?;
    drop(response);
    poll(cancelled)?;
    let trial = latest.ok_or(ResponseError::Invalid("missing final native response"))?;
    if !same_triples(&trial.forces_n, &contact.forces_n)
        || !same_triples(&trial.sample.paths_m, &contact.paths_m)
    {
        return Err(
            ResponseError::Invalid("final native replay differs from contact result").into(),
        );
    }
    let mut inputs = Vec::with_capacity(contacts.len());
    for (index, c) in contacts.iter().enumerate() {
        poll(cancelled)?;
        inputs.push(ContactAuditInput {
            id: c.id,
            gap_before_m: c.gap_m,
            gap_after_m: trial.sample.paths_m[index][0],
            slip_increment_m: trial.sample.increments_m[index],
            relative_velocity_before_m_s: before.velocities_m_s[index],
            relative_velocity_after_m_s: trial.sample.velocities_m_s[index],
            force_n: trial.forces_n[index],
            friction_coefficient: c.friction_coefficient,
            prescribed_velocity_m_s: c.prescribed_velocity_m_s,
            activation_tolerance_m: c.activation_tolerance_m,
        });
    }
    let mut body_work = Vec::with_capacity(motions.len());
    for (body, step) in trial.steps.iter().enumerate() {
        poll(cancelled)?;
        let mut authored = Sum::default();
        for (i, ((&force, &after), &before)) in loads[body]
            .mean_force_n
            .iter()
            .zip(step.motion.state().displacements_m())
            .zip(positions[body])
            .enumerate()
        {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            authored.add(force * finite(after - before)?)?;
        }
        let mut material = Sum::default();
        material.add(step.cohesive_energy_j.dissipated)?;
        material.add(step.viscous_dissipated_j)?;
        material.add(step.friction_dissipated_j)?;
        body_work.push(BodyWork {
            kinetic_change_j: finite(step.kinetic_after_j - step.kinetic_before_j)?,
            bulk_change_j: finite(step.bulk_stored_after_j - step.bulk_stored_before_j)?,
            cohesive_change_j: step.cohesive_energy_j.stored_change,
            material_loss_j: material.total()?,
            authored_work_j: authored.total()?,
            prescribed_jump_work_j: step.prescribed_jump_work_j,
            native_total_external_work_j: step.external_work_j,
            native_energy_defect_j: step.energy_residual_j,
        });
    }
    let accepted = audit::audit(&inputs, &body_work, dt_s, controls.audit, cancelled)?;
    poll(cancelled)?;
    Ok(SustainedCandidate {
        bodies: trial.steps,
        report: SustainedReport {
            elapsed_s: dt_s,
            target_time_s: target_s,
            contact_ids: contacts.iter().map(|c| c.id).collect(),
            contact,
            audit: accepted,
            body_work,
            body_trials: work.body_trials,
            residual_evaluations: work.residual_evaluations,
            native_newton_iterations: work.newton_iterations,
        },
    })
}

struct Trial<'m, 'a> {
    steps: Vec<DynamicsStep<'m, 'a>>,
    forces_n: Vec<[f64; 3]>,
    sample: Paths,
}
struct Paths {
    paths_m: Vec<[f64; 3]>,
    increments_m: Vec<[f64; 3]>,
    velocities_m_s: Vec<[f64; 3]>,
}

#[allow(clippy::too_many_arguments)]
fn evaluate<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    contacts: &[SustainedContact<'_>],
    loads: &[BodyLoad<'_>],
    reference: &[&[f64]],
    forces: &[[f64; 3]],
    target_s: f64,
    dt_s: f64,
    controls: SustainedControls,
    work: &mut Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Trial<'m, 'a>, ResponseError> {
    poll(cancelled)?;
    if forces.len() != contacts.len() {
        return Err(ResponseError::Invalid(
            "one complete force per original contact",
        ));
    }
    let mut assembled = Vec::with_capacity(motions.len());
    for load in loads {
        poll(cancelled)?;
        let mut body = Vec::with_capacity(load.mean_force_n.len());
        for (i, &value) in load.mean_force_n.iter().enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            let mut sum = Sum::default();
            sum.add(value)?;
            body.push(sum);
        }
        assembled.push(body);
    }
    for (contact, force) in contacts.iter().zip(forces) {
        poll(cancelled)?;
        for axis in 0..3 {
            finite(force[axis])?;
            if force[axis].abs() > controls.maximum_contact_force_n {
                return Err(ResponseError::Budget("contact force applicability"));
            }
            for (i, term) in contact.trace[axis].iter().enumerate() {
                if i % 256 == 0 {
                    poll(cancelled)?;
                }
                assembled[term.body][term.dof].add(term.coefficient * force[axis])?;
            }
        }
    }
    let mut steps = Vec::with_capacity(motions.len());
    for (body, (motion, sums)) in motions.iter().zip(assembled).enumerate() {
        poll(cancelled)?;
        let mut total_force = Vec::with_capacity(sums.len());
        for (i, sum) in sums.into_iter().enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            total_force.push(sum.total()?);
        }
        let reserve = work.reserve(motion, controls)?;
        let step = motion
            .advance_to_time_controlled(
                &total_force,
                loads[body].endpoint_prescribed_jumps_m,
                target_s,
                cancelled,
            )
            .map_err(|error| ResponseError::Body { body, error })?;
        work.refund(reserve, &step)?;
        if step.nonlinear.history.capacity() > history_capacity(motion)? {
            return Err(ResponseError::Budget("retained native Newton history"));
        }
        steps.push(step);
    }
    let positions: Vec<_> = steps
        .iter()
        .map(|s| s.motion.state().displacements_m())
        .collect();
    let velocities: Vec<_> = steps.iter().map(|s| s.motion.velocity_m_s()).collect();
    let sample = sample_paths(
        contacts,
        reference,
        &positions,
        &velocities,
        dt_s,
        cancelled,
    )?;
    poll(cancelled)?;
    Ok(Trial {
        steps,
        forces_n: forces.to_vec(),
        sample,
    })
}

fn sample_paths(
    contacts: &[SustainedContact<'_>],
    reference: &[&[f64]],
    positions: &[&[f64]],
    velocities: &[&[f64]],
    dt_s: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Paths, ResponseError> {
    geometry::validate_vectors(positions, reference, cancelled)?;
    geometry::validate_vectors(velocities, reference, cancelled)?;
    let mut paths_m = Vec::with_capacity(contacts.len());
    let mut increments_m = Vec::with_capacity(contacts.len());
    let mut velocities_m_s = Vec::with_capacity(contacts.len());
    for contact in contacts {
        poll(cancelled)?;
        let mut increment = [0.0; 3];
        let mut velocity = [0.0; 3];
        let mut normal_gap = Sum::default();
        normal_gap.add(contact.gap_m)?;
        normal_gap.add(dt_s * contact.prescribed_velocity_m_s[0])?;
        for axis in 0..3 {
            let mut delta = Sum::default();
            let mut speed = Sum::default();
            delta.add(dt_s * contact.prescribed_velocity_m_s[axis])?;
            speed.add(contact.prescribed_velocity_m_s[axis])?;
            for (i, term) in contact.trace[axis].iter().enumerate() {
                if i % 256 == 0 {
                    poll(cancelled)?;
                }
                let displacement = finite(
                    geometry::coordinate(positions, term)? - geometry::coordinate(reference, term)?,
                )?;
                let contribution = finite(term.coefficient * displacement)?;
                delta.add(contribution)?;
                if axis == 0 {
                    normal_gap.add(contribution)?;
                }
                speed.add(term.coefficient * geometry::coordinate(velocities, term)?)?;
            }
            increment[axis] = delta.total()?;
            velocity[axis] = speed.total()?;
        }
        paths_m.push([normal_gap.total()?, increment[1], increment[2]]);
        increments_m.push(increment);
        velocities_m_s.push(velocity);
    }
    poll(cancelled)?;
    Ok(Paths {
        paths_m,
        increments_m,
        velocities_m_s,
    })
}

#[allow(clippy::too_many_arguments)]
fn admit(
    motions: &[&CohesiveMotion<'_, '_>],
    contacts: &[SustainedContact<'_>],
    loads: &[BodyLoad<'_>],
    scales: &[NaturalMapScales],
    initial: &[[f64; 3]],
    duration_s: f64,
    controls: SustainedControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ResponseError> {
    poll(cancelled)?;
    impact::admit_states(
        motions,
        controls.maximum_bodies,
        controls.maximum_workspace_entries,
        cancelled,
    )?;
    if contacts.is_empty()
        || contacts.len() > 21
        || contacts.len() > controls.solve.maximum_contacts
        || scales.len() != contacts.len()
        || initial.len() != contacts.len()
        || loads.len() != motions.len()
        || !positive(duration_s)
        || !positive(controls.maximum_contact_force_n)
        || controls.maximum_body_trials == 0
        || controls.maximum_residual_evaluations == 0
        || controls.maximum_native_newton_iterations == 0
        || !positive(controls.audit.force_tolerance_n)
        || !controls.audit.gap_tolerance_m.is_finite()
        || controls.audit.gap_tolerance_m < 0.0
        || !positive(controls.audit.slip_tolerance_m)
        || !positive(controls.audit.velocity_tolerance_m_s)
        || !positive(controls.audit.energy_tolerance_j)
    {
        return Err(ResponseError::Invalid(
            "bounded complete sustained-contact inputs required",
        ));
    }
    let reference: Vec<_> = motions
        .iter()
        .map(|m| m.state().displacements_m())
        .collect();
    geometry::validate_vectors(&reference, &reference, cancelled)?;
    let mut workspace = mul(motions.len(), 48)?;
    workspace = add(workspace, mul(contacts.len(), 64)?)?;
    for (motion, load) in motions.iter().zip(loads) {
        poll(cancelled)?;
        let n = motion.velocity_m_s().len();
        let patches = motion.state().points().len();
        if load.mean_force_n.len() != n || load.endpoint_prescribed_jumps_m.len() != patches {
            return Err(ResponseError::Invalid(
                "original load and prescribed-jump dimensions",
            ));
        }
        // Retained native reports include configuration and full Newton
        // history in addition to the original physical vectors and patches.
        workspace = add(
            workspace,
            core::mem::size_of::<DynamicsStep<'_, '_>>().div_ceil(core::mem::size_of::<f64>()),
        )?;
        workspace = add(
            workspace,
            mul(
                history_capacity(motion)?,
                core::mem::size_of::<fs_solver::NewtonIteration>()
                    .div_ceil(core::mem::size_of::<f64>()),
            )?,
        )?;
        workspace = add(workspace, mul(n, 12)?)?;
        workspace = add(workspace, mul(patches, 96)?)?;
        if workspace > controls.maximum_workspace_entries {
            return Err(ResponseError::Budget("held-force adapter workspace"));
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
    let mut terms = 0usize;
    for (i, contact) in contacts.iter().enumerate() {
        poll(cancelled)?;
        if contacts[..i].iter().any(|c| c.id == contact.id)
            || !contact.gap_m.is_finite()
            || contact.gap_m < 0.0
            || !contact.activation_tolerance_m.is_finite()
            || contact.activation_tolerance_m < 0.0
            || contact.gap_m > contact.activation_tolerance_m
            || !contact.friction_coefficient.is_finite()
            || contact.friction_coefficient < 0.0
            || contact
                .prescribed_velocity_m_s
                .iter()
                .any(|x| !x.is_finite())
            || initial[i]
                .iter()
                .any(|x| !x.is_finite() || x.abs() > controls.maximum_contact_force_n)
        {
            return Err(ResponseError::Invalid(
                "unique already-active nonpenetrating contacts and finite law data required",
            ));
        }
        for row in contact.trace {
            terms = add(terms, row.len())?;
            if terms > controls.maximum_trace_terms {
                return Err(ResponseError::Budget("sustained contact trace terms"));
            }
            geometry::admit_trace(row, &reference, cancelled)?;
        }
    }
    poll(cancelled)
}

// Native vector growth is admitted conservatively, then the actual retained
// capacity is checked before a body receipt enters the candidate. Transient
// allocations inside the original solver retain that owner's controls.
fn history_capacity(motion: &CohesiveMotion<'_, '_>) -> Result<usize, ResponseError> {
    Ok(mul(motion.dynamics_config().solve.max_newton_iterations, 2)?.max(4))
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
        c: SustainedControls,
    ) -> Result<(usize, usize), ResponseError> {
        let native = motion.dynamics_config().solve;
        let trials = add(self.body_trials, 1)?;
        let residuals = add(self.residual_evaluations, native.max_operator_evaluations)?;
        let iterations = add(self.newton_iterations, native.max_newton_iterations)?;
        if trials > c.maximum_body_trials
            || residuals > c.maximum_residual_evaluations
            || iterations > c.maximum_native_newton_iterations
        {
            return Err(ResponseError::Budget("original held-force body work"));
        }
        self.body_trials = trials;
        self.residual_evaluations = residuals;
        self.newton_iterations = iterations;
        Ok((
            native.max_operator_evaluations,
            native.max_newton_iterations,
        ))
    }
    fn refund(
        &mut self,
        reserve: (usize, usize),
        step: &DynamicsStep<'_, '_>,
    ) -> Result<(), ResponseError> {
        if step.residual_evaluations > reserve.0 || step.nonlinear.iterations > reserve.1 {
            return Err(ResponseError::Invalid(
                "native response exceeded reserved work",
            ));
        }
        self.residual_evaluations -= reserve.0 - step.residual_evaluations;
        self.newton_iterations -= reserve.1 - step.nonlinear.iterations;
        Ok(())
    }
}
fn same_triples(a: &[[f64; 3]], b: &[[f64; 3]]) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b)
            .all(|(a, b)| a.iter().zip(b).all(|(a, b)| a.to_bits() == b.to_bits()))
}
fn finite(value: f64) -> Result<f64, ResponseError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ResponseError::NonFinite("derived held-force response"))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn add(a: usize, b: usize) -> Result<usize, ResponseError> {
    a.checked_add(b)
        .ok_or(ResponseError::Budget("work/count overflow"))
}
fn mul(a: usize, b: usize) -> Result<usize, ResponseError> {
    a.checked_mul(b)
        .ok_or(ResponseError::Budget("workspace overflow"))
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ResponseError> {
    if cancelled() {
        Err(ResponseError::Cancelled)
    } else {
        Ok(())
    }
}
