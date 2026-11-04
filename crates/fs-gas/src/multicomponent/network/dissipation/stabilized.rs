//! RKL2 on the original branch/finite-plenum spatial owner.
//!
//! The axial owner supplies the Legendre coefficients and branch recurrence.
//! Plenum inventory and every signed endpoint transfer use those same weights;
//! there is no independently subcycled plenum, approximate port or flux repair.

use super::evolve::{self, Input, Progress, compensated_add, forward, spatial};
use super::{
    Cell, DissipationError, MixtureNetwork1d, NetworkDissipationOperator, NetworkDissipationReport,
    Time, TransportLaw, allocated, poll, positive,
};
use crate::multicomponent::dissipation::StabilizedDissipationControl;
use crate::multicomponent::dissipation::evolve::Stage;
use crate::multicomponent::dissipation::stabilized::{
    b, check_rate, extent, recurrence, reserve_evaluation,
};

/// Complete network balances and measured RKL2 spatial work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedNetworkDissipationReport {
    /// Original species, energy, branch and finite-port accounting. Accepted
    /// steps are RKL2 macrosteps; the reported parabolic number is unscaled.
    pub transport: NetworkDissipationReport,
    /// Actual RHS calls, including refused trials and endpoint admission.
    pub operator_evaluations: usize,
    /// Largest accepted stage count; zero for a zero-duration interval.
    pub maximum_stages_used: usize,
    /// Largest accepted `dt * rate / ((s*s+s-2)/4)`, bounded by the base CFL.
    pub maximum_scaled_parabolic_number: f64,
}

impl<L: TransportLaw> NetworkDissipationOperator<L> {
    /// Advance branches and finite plenums together with bounded RKL2.
    ///
    /// This opt-in path retains actual finite storage, all incident port
    /// conductances, and the original Maxwell--Stefan/heat/stress face kernel
    /// at every stage. Port receipts and exterior transfers use the same RKL2
    /// recurrence as the field, including its negative weights. The existing
    /// three state arrays are reused; extra face, plenum-RHS and endpoint
    /// buffers are allocated lazily once, never once per stage.
    ///
    /// `controls.transport.maximum_step_s` is the caller's accuracy cap on a
    /// physical macrostep. There is no temporal error estimator or unconditional
    /// nonlinear positivity guarantee. Every stage, including the endpoint,
    /// must pass the original source, caloric and residual checks and the
    /// globally estimated rate bound. Refused trials are halved without clipping.
    /// Complete-duration failure preserves fields, clocks and accepted reports.
    pub fn advance_stabilized_for(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: StabilizedDissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedNetworkDissipationReport, DissipationError> {
        let report = self.advance_stabilized_for_trial(
            network,
            duration,
            controls,
            None,
            &mut 0,
            &mut cancelled,
        )?;
        self.publish_trial_report();
        Ok(report)
    }

    /// The same stabilized interval with stationary prescribed exterior baths.
    ///
    /// Bath ordering, canonical energy gauges and finite-port exclusion follow
    /// [`Self::advance_with_reservoirs_for`]. All actual half-cell distances,
    /// bath enthalpies and exterior balances remain in the shared face owner.
    pub fn advance_stabilized_with_reservoirs_for(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: StabilizedDissipationControl,
        reservoirs: &[[Option<Cell<'_>>; 2]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedNetworkDissipationReport, DissipationError> {
        let report = self.advance_stabilized_for_trial(
            network,
            duration,
            controls,
            Some(reservoirs),
            &mut 0,
            &mut cancelled,
        )?;
        self.publish_trial_report();
        Ok(report)
    }

    // A composing owner carries the same actual RHS counter through both
    // molecular halves and discarded complete coupled trials.
    pub(in crate::multicomponent::network) fn advance_stabilized_for_trial(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: StabilizedDissipationControl,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StabilizedNetworkDissipationReport, DissipationError> {
        self.admit_plenum_volumes(network, cancelled)?;
        self.advance_stabilized_for_moving_trial(
            network,
            duration,
            controls,
            reservoirs,
            evaluations,
            cancelled,
        )
    }

    // The original composing owner must already have admitted geometry motion.
    // Geometry and reports remain unpublished until its whole-step commit.
    pub(in crate::multicomponent::network) fn advance_stabilized_for_moving_trial(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: StabilizedDissipationControl,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StabilizedNetworkDissipationReport, DissipationError> {
        self.trial_report = None;
        let report = advance(
            self,
            network,
            duration,
            controls,
            reservoirs,
            evaluations,
            cancelled,
        )?;
        self.trial_report = Some(report.transport);
        Ok(report)
    }
}

fn advance<L: TransportLaw>(
    owner: &mut NetworkDissipationOperator<L>,
    network: &mut MixtureNetwork1d,
    duration: Time,
    controls: StabilizedDissipationControl,
    reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StabilizedNetworkDissipationReport, DissipationError> {
    controls.validate()?;
    let base = controls.transport;
    let (before, active_reservoirs) =
        evolve::prepare(owner, network, duration, base, reservoirs, cancelled)?;
    let reservoirs = if active_reservoirs { reservoirs } else { None };
    let initial_evaluations = *evaluations;
    let mut elapsed = 0.0;
    let mut progress = Progress::default();
    let mut maximum_stages_used = 0;
    let mut maximum_scaled_parabolic = 0.0_f64;
    if duration.value() > 0.0 {
        prepare_work(owner)?;
    }
    while elapsed < duration.value() {
        poll(cancelled)?;
        if progress.accepted_steps + progress.rejected_steps >= base.maximum_steps {
            return Err(DissipationError::WorkLimit {
                elapsed_time_s: elapsed,
                accepted_steps: progress.accepted_steps,
                rejected_steps: progress.rejected_steps,
            });
        }
        reserve_evaluation(evaluations, controls.maximum_operator_evaluations, elapsed)?;
        let a = spatial(owner, Input::Current, base, reservoirs, cancelled)?;
        cache_anchor(owner);
        let remaining = duration.value() - elapsed;
        let mut step = remaining.min(base.maximum_step_s);
        if a.rate > 0.0 {
            step = step.min(base.parabolic_cfl * extent(controls.maximum_stages) / a.rate);
        }
        let mut accepted = false;
        let mut last = DissipationError::Invalid("no admissible stabilized molecular network step");
        for _ in 0..=base.maximum_halvings {
            poll(cancelled)?;
            if progress.accepted_steps + progress.rejected_steps >= base.maximum_steps {
                return Err(DissipationError::WorkLimit {
                    elapsed_time_s: elapsed,
                    accepted_steps: progress.accepted_steps,
                    rejected_steps: progress.rejected_steps,
                });
            }
            if !positive(step)
                || elapsed + step == elapsed
                || (step < base.minimum_step_s && remaining > base.minimum_step_s)
            {
                return Err(DissipationError::StepFailure {
                    elapsed_time_s: elapsed,
                    step_s: step,
                    last: Box::new(last),
                });
            }
            let stages = (2..=controls.maximum_stages)
                .find(|s| step * a.rate <= base.parabolic_cfl * extent(*s))
                .unwrap_or(controls.maximum_stages);
            let extension = extent(stages);
            let attempt = macrostep(
                owner,
                Step {
                    step_s: step,
                    stages,
                    anchor: a,
                    elapsed,
                },
                controls,
                reservoirs,
                evaluations,
                cancelled,
            );
            match attempt {
                Ok(peak) => {
                    for k in 0..owner.work.endpoint_changes.len() {
                        if k % 1024 == 0 {
                            poll(cancelled)?;
                        }
                        compensated_add(
                            &mut owner.work.endpoint_changes[k],
                            &mut owner.work.endpoint_corrections[k],
                            owner.work.stabilized_endpoint_changes[k],
                        )?;
                    }
                    maximum_stages_used = maximum_stages_used.max(stages);
                    progress.maximum_parabolic = progress.maximum_parabolic.max(step * peak.rate);
                    maximum_scaled_parabolic =
                        maximum_scaled_parabolic.max(step * peak.rate / extension);
                    progress.maximum_gradient =
                        progress.maximum_gradient.max(peak.gradient_residual);
                    progress.maximum_mass_frame =
                        progress.maximum_mass_frame.max(peak.mass_frame_residual);
                    core::mem::swap(&mut owner.work.current, &mut owner.work.stage);
                    elapsed = if step == remaining {
                        duration.value()
                    } else {
                        elapsed + step
                    };
                    progress.accepted_steps += 1;
                    accepted = true;
                    break;
                }
                Err(
                    error @ (DissipationError::Cancelled
                    | DissipationError::OperatorWorkLimit { .. }),
                ) => {
                    return Err(error);
                }
                Err(error) => {
                    progress.rejected_steps += 1;
                    last = error;
                    step *= 0.5;
                }
            }
        }
        if !accepted {
            return Err(DissipationError::StepFailure {
                elapsed_time_s: elapsed,
                step_s: step,
                last: Box::new(last),
            });
        }
    }
    let transport = evolve::finish(owner, network, duration, before, progress, cancelled)?;
    Ok(StabilizedNetworkDissipationReport {
        transport,
        operator_evaluations: *evaluations - initial_evaluations,
        maximum_stages_used,
        maximum_scaled_parabolic_number: maximum_scaled_parabolic,
    })
}

pub(super) fn prepare_work<L: TransportLaw>(
    owner: &mut NetworkDissipationOperator<L>,
) -> Result<(), DissipationError> {
    let work = &mut owner.work;
    if work.stabilized_faces.len() != work.faces.len() {
        work.stabilized_faces = allocated(work.faces.len(), 0.0)?;
    }
    if work.stabilized_plenum_rhs.len() != work.plenum_rhs.len() {
        work.stabilized_plenum_rhs = allocated(work.plenum_rhs.len(), 0.0)?;
    }
    if work.stabilized_endpoint_changes.len() != work.endpoints_a.len() {
        work.stabilized_endpoint_changes = allocated(work.endpoints_a.len(), 0.0)?;
    }
    Ok(())
}

pub(super) fn cache_anchor<L: TransportLaw>(owner: &mut NetworkDissipationOperator<L>) {
    owner
        .work
        .stabilized_faces
        .copy_from_slice(&owner.work.faces);
    owner
        .work
        .stabilized_plenum_rhs
        .copy_from_slice(&owner.work.plenum_rhs);
}

pub(super) struct Step {
    pub step_s: f64,
    pub stages: usize,
    pub anchor: Stage,
    pub elapsed: f64,
}

// One original RKL2 recurrence. `current` and the cached anchor RHS are input;
// `stage`, endpoint transfers and the fully admitted endpoint RHS are output.
// Keeping this one kernel lets accuracy probes share the constitutive physics
// and the exact same fixed-stage integration formula with the original path.
pub(super) fn macrostep<L: TransportLaw>(
    owner: &mut NetworkDissipationOperator<L>,
    trial: Step,
    controls: StabilizedDissipationControl,
    reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Stage, DissipationError> {
    let Step {
        step_s: step,
        stages,
        anchor,
        elapsed,
    } = trial;
    let base = controls.transport;
    let m = owner.gas.species_count() + 2;
    let extension = extent(stages);
    let limit = base.parabolic_cfl * extension;
    let w1 = 1.0 / extension;
    let first_weight = b(1) * w1;
    let mut peak = anchor;
    check_rate(step, anchor.rate, limit)?;
    // A refused trial overwrote the spatial scratch. Restore its cached RHS.
    owner
        .work
        .faces
        .copy_from_slice(&owner.work.stabilized_faces);
    owner
        .work
        .plenum_rhs
        .copy_from_slice(&owner.work.stabilized_plenum_rhs);
    forward(owner, Input::Current, step * first_weight, cancelled)?;
    owner.work.candidate.copy_from_slice(&owner.work.current);
    owner.work.endpoints_a.fill(0.0);
    for (index, shape) in owner.branches.iter().enumerate() {
        poll(cancelled)?;
        for end in 0..2 {
            let face = shape.face_start + if end == 0 { 0 } else { shape.mesh.cell_count() };
            let sign = if end == 0 { 1.0 } else { -1.0 };
            for k in 0..m {
                owner.work.stabilized_endpoint_changes[(2 * index + end) * m + k] = step
                    * first_weight
                    * shape.mesh.face_areas()[0]
                    * sign
                    * owner.work.stabilized_faces[face * m + k];
            }
        }
    }
    for j in 2..=stages + 1 {
        poll(cancelled)?;
        reserve_evaluation(evaluations, controls.maximum_operator_evaluations, elapsed)?;
        let evaluated = spatial(owner, Input::Stage, base, reservoirs, cancelled)?;
        check_rate(step, evaluated.rate, limit)?;
        peak.rate = peak.rate.max(evaluated.rate);
        peak.gradient_residual = peak.gradient_residual.max(evaluated.gradient_residual);
        peak.mass_frame_residual = peak.mass_frame_residual.max(evaluated.mass_frame_residual);
        if j > stages {
            break; // endpoint passed the full original spatial admission
        }
        let jf = j as f64;
        let mu = ((2.0 * jf - 1.0) / jf) * b(j) / b(j - 1);
        let nu = -((jf - 1.0) / jf) * b(j) / b(j - 2);
        let mu_tilde = mu * w1;
        let gamma_tilde = -(1.0 - b(j - 1)) * mu_tilde;
        stage_recurrence(owner, step, [mu, nu, mu_tilde, gamma_tilde], cancelled)?;
        core::mem::swap(&mut owner.work.stage, &mut owner.work.candidate);
    }
    Ok(peak)
}

fn stage_recurrence<L: TransportLaw>(
    owner: &mut NetworkDissipationOperator<L>,
    step: f64,
    coefficients: [f64; 4],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DissipationError> {
    let [mu, nu, mu_tilde, gamma_tilde] = coefficients;
    let NetworkDissipationOperator {
        gas,
        branches,
        plenum_volumes,
        plenum_cell_start,
        work,
        ..
    } = owner;
    let n = gas.species_count();
    let m = n + 2;
    for (index, shape) in branches.iter().enumerate() {
        let count = shape.mesh.cell_count();
        let cs = shape.cell_start * m;
        let ce = cs + count * m;
        let fs = shape.face_start * m;
        let fe = fs + (count + 1) * m;
        recurrence(
            gas,
            &shape.mesh,
            &work.current[cs..ce],
            &work.stage[cs..ce],
            &mut work.candidate[cs..ce],
            &work.stabilized_faces[fs..fe],
            &work.faces[fs..fe],
            step,
            coefficients,
            cancelled,
        )?;
        for end in 0..2 {
            poll(cancelled)?;
            let face = shape.face_start + if end == 0 { 0 } else { count };
            let sign = if end == 0 { 1.0 } else { -1.0 };
            for k in 0..m {
                let target = (2 * index + end) * m + k;
                let previous = work.stabilized_endpoint_changes[target];
                let next = mu * previous
                    + nu * work.endpoints_a[target]
                    + step
                        * shape.mesh.face_areas()[0]
                        * sign
                        * (mu_tilde * work.faces[face * m + k]
                            + gamma_tilde * work.stabilized_faces[face * m + k]);
                if !next.is_finite() {
                    return Err(DissipationError::Invalid(
                        "stabilized molecular endpoint overflow",
                    ));
                }
                work.endpoints_a[target] = previous;
                work.stabilized_endpoint_changes[target] = next;
            }
        }
    }
    for (index, &volume) in plenum_volumes.iter().enumerate() {
        poll(cancelled)?;
        let start = (*plenum_cell_start + index) * m;
        let destination = &mut work.candidate[start..start + m];
        for (k, value) in destination.iter_mut().enumerate() {
            let anchor = work.current[start + k];
            *value = anchor
                + mu * (work.stage[start + k] - anchor)
                + nu * (*value - anchor)
                + (step / volume)
                    * (mu_tilde * work.plenum_rhs[index * m + k]
                        + gamma_tilde * work.stabilized_plenum_rhs[index * m + k]);
        }
        destination[n] = 0.0; // fixed support receives each port's axial impulse
        if !gas.admissible(destination) {
            return Err(DissipationError::Invalid(
                "stabilized molecular plenum stage left species/caloric domain",
            ));
        }
    }
    Ok(())
}
