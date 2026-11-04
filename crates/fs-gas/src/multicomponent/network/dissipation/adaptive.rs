//! Physical local-error control on the original conservative network RKL2.

use super::evolve::{self, Input, Progress, compensated_add, spatial};
use super::stabilized::{self, Step, macrostep};
use super::{
    Cell, DissipationError, MixtureNetwork1d, NetworkDissipationOperator,
    StabilizedNetworkDissipationReport, Time, TransportLaw, allocated, poll, positive,
};
use crate::multicomponent::dissipation::StabilizedDissipationControl;
use crate::multicomponent::dissipation::stabilized::{extent, reserve_evaluation};

/// Accuracy and bounded work for full/two-half RKL2 molecular intervals.
///
/// Every species density, temperature and axial velocity is checked separately
/// in every branch cell and finite plenum. Each signed port transfer is also
/// checked, including external baths and the impulse delivered to fixed supports.
/// Relative scales use the actual quantity, not total mixture mass or absolute
/// formation energy. These are estimated LOCAL errors of the molecular operator;
/// spatial, splitting and global-error certification remain separate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveNetworkDissipationControl {
    /// Original stage, constitutive, interval-attempt and complete-RHS budgets.
    /// The maximum/minimum steps bound an accuracy-controlled interval; its
    /// accepted solution contains two half-steps. `maximum_halvings` also bounds
    /// adaptive retries, whose shrink factor lies in [0.2,0.8].
    pub stabilized: StabilizedDissipationControl,
    /// Finite relative accuracy in [0,1e-3], applied component by component.
    pub relative_tolerance: f64,
    /// Positive absolute accuracy for each original species density [kg/m³].
    pub absolute_partial_density_tolerance_kg_m3: f64,
    /// Positive absolute temperature accuracy [K].
    pub absolute_temperature_tolerance_k: f64,
    /// Positive absolute axial-velocity accuracy [m/s].
    pub absolute_velocity_tolerance_m_s: f64,
}

impl Default for AdaptiveNetworkDissipationControl {
    fn default() -> Self {
        Self {
            stabilized: StabilizedDissipationControl::default(),
            relative_tolerance: 1e-6,
            absolute_partial_density_tolerance_kg_m3: 1e-10,
            absolute_temperature_tolerance_k: 1e-6,
            absolute_velocity_tolerance_m_s: 1e-7,
        }
    }
}

impl AdaptiveNetworkDissipationControl {
    /// Admit local accuracy independently of the molecular constitutive budget.
    pub fn validate(self) -> Result<(), DissipationError> {
        self.stabilized.validate()?;
        if !self.relative_tolerance.is_finite()
            || !(0.0..=1e-3).contains(&self.relative_tolerance)
            || !positive(self.absolute_partial_density_tolerance_kg_m3)
            || !positive(self.absolute_temperature_tolerance_k)
            || !positive(self.absolute_velocity_tolerance_m_s)
        {
            return Err(DissipationError::Invalid(
                "adaptive molecular accuracy controls",
            ));
        }
        Ok(())
    }
}

/// Accepted physical interval, complete probe work and measured local accuracy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveNetworkDissipationReport {
    /// The original balances and actual RHS calls, including all full-step,
    /// half-step and refused probes. Accepted steps count paired half-steps.
    pub stabilized: StabilizedNetworkDissipationReport,
    /// Largest accepted maximum component/port error ratio, at most one.
    pub maximum_estimated_local_error: f64,
    /// Complete probes rejected specifically by the local accuracy test.
    pub temporal_rejections: usize,
}

#[derive(Default)]
pub(super) struct Workspace {
    anchor: Vec<f64>,
    coarse: Vec<f64>,
    coarse_endpoints: Vec<f64>,
    first_endpoints: Vec<f64>,
}

impl<L: TransportLaw> NetworkDissipationOperator<L> {
    /// Evolve a complete rigid molecular interval to declared local accuracy.
    ///
    /// Full and half probes use the SAME RKL2 stage count. Thus changing the
    /// stability-selected stage count cannot invalidate the second-order
    /// difference/3 estimate within a probe. The admitted two-half field and its
    /// original conservative transfers are published; no Richardson-extrapolated
    /// state, species clipping or flux correction is used. Two additional fields
    /// and two endpoint arrays are retained lazily, with no per-probe allocation.
    pub fn advance_adaptive_for(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: AdaptiveNetworkDissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptiveNetworkDissipationReport, DissipationError> {
        let report = self.advance_adaptive_for_trial(
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

    /// The same accuracy-controlled interval with stationary exterior baths.
    pub fn advance_adaptive_with_reservoirs_for(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: AdaptiveNetworkDissipationControl,
        reservoirs: &[[Option<Cell<'_>>; 2]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptiveNetworkDissipationReport, DissipationError> {
        let report = self.advance_adaptive_for_trial(
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

    pub(in crate::multicomponent::network) fn advance_adaptive_for_trial(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: AdaptiveNetworkDissipationControl,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AdaptiveNetworkDissipationReport, DissipationError> {
        self.admit_plenum_volumes(network, cancelled)?;
        self.advance_adaptive_for_moving_trial(
            network,
            duration,
            controls,
            reservoirs,
            evaluations,
            cancelled,
        )
    }

    // A composing owner admits the actual geometry before each molecular half.
    pub(in crate::multicomponent::network) fn advance_adaptive_for_moving_trial(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: AdaptiveNetworkDissipationControl,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AdaptiveNetworkDissipationReport, DissipationError> {
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
        self.trial_report = Some(report.stabilized.transport);
        Ok(report)
    }
}

fn advance<L: TransportLaw>(
    owner: &mut NetworkDissipationOperator<L>,
    network: &mut MixtureNetwork1d,
    duration: Time,
    controls: AdaptiveNetworkDissipationControl,
    reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AdaptiveNetworkDissipationReport, DissipationError> {
    controls.validate()?;
    let stabilized = controls.stabilized;
    let base = stabilized.transport;
    let (before, active) = evolve::prepare(owner, network, duration, base, reservoirs, cancelled)?;
    let reservoirs = if active { reservoirs } else { None };
    let initial_evaluations = *evaluations;
    let mut progress = Progress::default();
    let (mut elapsed, mut next_step) = (0.0, base.maximum_step_s);
    let (mut maximum_stages_used, mut temporal_rejections) = (0, 0);
    let mut maximum_scaled_parabolic = 0.0_f64;
    let mut maximum_error = 0.0_f64;
    if duration.value() > 0.0 {
        prepare_work(owner)?;
    }
    while elapsed < duration.value() {
        let remaining = duration.value() - elapsed;
        let mut step = remaining.min(next_step);
        copy(
            &mut owner.work.adaptive.anchor,
            &owner.work.current,
            cancelled,
        )?;
        let mut accepted = false;
        let mut last = DissipationError::Invalid("no accurate molecular network step");
        for _ in 0..=base.maximum_halvings {
            poll(cancelled)?;
            if progress.accepted_steps + progress.rejected_steps >= base.maximum_steps {
                return Err(DissipationError::WorkLimit {
                    elapsed_time_s: elapsed,
                    accepted_steps: progress.accepted_steps,
                    rejected_steps: progress.rejected_steps,
                });
            }
            reserve_evaluation(
                evaluations,
                stabilized.maximum_operator_evaluations,
                elapsed,
            )?;
            let anchor = spatial(owner, Input::Current, base, reservoirs, cancelled)?;
            if anchor.rate > 0.0 {
                step =
                    step.min(base.parabolic_cfl * extent(stabilized.maximum_stages) / anchor.rate);
            }
            if !positive(0.5 * step)
                || elapsed + step == elapsed
                || (step < base.minimum_step_s && remaining > base.minimum_step_s)
            {
                return Err(DissipationError::StepFailure {
                    elapsed_time_s: elapsed,
                    step_s: step,
                    last: Box::new(last),
                });
            }
            let stages = (2..=stabilized.maximum_stages)
                .find(|s| step * anchor.rate <= base.parabolic_cfl * extent(*s))
                .unwrap_or(stabilized.maximum_stages);
            let half = 0.5 * step;
            stabilized::cache_anchor(owner);
            let attempt = (|| {
                let coarse = macrostep(
                    owner,
                    Step {
                        step_s: step,
                        stages,
                        anchor,
                        elapsed,
                    },
                    stabilized,
                    reservoirs,
                    evaluations,
                    cancelled,
                )?;
                copy(
                    &mut owner.work.adaptive.coarse,
                    &owner.work.stage,
                    cancelled,
                )?;
                copy(
                    &mut owner.work.adaptive.coarse_endpoints,
                    &owner.work.stabilized_endpoint_changes,
                    cancelled,
                )?;
                // The unchanged cached original RHS serves the full and first
                // half probe. The admitted midpoint RHS serves the second half.
                let first = macrostep(
                    owner,
                    Step {
                        step_s: half,
                        stages,
                        anchor,
                        elapsed,
                    },
                    stabilized,
                    reservoirs,
                    evaluations,
                    cancelled,
                )?;
                copy(
                    &mut owner.work.adaptive.first_endpoints,
                    &owner.work.stabilized_endpoint_changes,
                    cancelled,
                )?;
                core::mem::swap(&mut owner.work.current, &mut owner.work.stage);
                stabilized::cache_anchor(owner);
                let second = macrostep(
                    owner,
                    Step {
                        step_s: half,
                        stages,
                        anchor: first,
                        elapsed,
                    },
                    stabilized,
                    reservoirs,
                    evaluations,
                    cancelled,
                )?;
                let error = estimate(owner, controls, cancelled)?;
                if error > 1.0 {
                    return Err(DissipationError::EstimatedTemporalError {
                        observed: error,
                        maximum: 1.0,
                    });
                }
                Ok((coarse, first, second, error))
            })();
            match attempt {
                Ok((coarse, first, second, error)) => {
                    for k in 0..owner.work.endpoint_changes.len() {
                        if k % 1024 == 0 {
                            poll(cancelled)?;
                        }
                        compensated_add(
                            &mut owner.work.endpoint_changes[k],
                            &mut owner.work.endpoint_corrections[k],
                            owner.work.adaptive.first_endpoints[k],
                        )?;
                        compensated_add(
                            &mut owner.work.endpoint_changes[k],
                            &mut owner.work.endpoint_corrections[k],
                            owner.work.stabilized_endpoint_changes[k],
                        )?;
                    }
                    let number = (step * coarse.rate).max(half * first.rate.max(second.rate));
                    progress.maximum_parabolic = progress.maximum_parabolic.max(number);
                    maximum_scaled_parabolic =
                        maximum_scaled_parabolic.max(number / extent(stages));
                    progress.maximum_gradient = progress.maximum_gradient.max(
                        coarse
                            .gradient_residual
                            .max(first.gradient_residual)
                            .max(second.gradient_residual),
                    );
                    progress.maximum_mass_frame = progress.maximum_mass_frame.max(
                        coarse
                            .mass_frame_residual
                            .max(first.mass_frame_residual)
                            .max(second.mass_frame_residual),
                    );
                    maximum_stages_used = maximum_stages_used.max(stages);
                    maximum_error = maximum_error.max(error);
                    core::mem::swap(&mut owner.work.current, &mut owner.work.stage);
                    elapsed = if step == remaining {
                        duration.value()
                    } else {
                        elapsed + step
                    };
                    progress.accepted_steps += 1;
                    next_step = (step * factor(error, true)).min(base.maximum_step_s);
                    accepted = true;
                    break;
                }
                Err(error) => {
                    copy(
                        &mut owner.work.current,
                        &owner.work.adaptive.anchor,
                        cancelled,
                    )?;
                    if matches!(
                        error,
                        DissipationError::Cancelled | DissipationError::OperatorWorkLimit { .. }
                    ) {
                        return Err(error);
                    }
                    let shrink =
                        if let DissipationError::EstimatedTemporalError { observed, .. } = error {
                            temporal_rejections += 1;
                            factor(observed, false)
                        } else {
                            0.5
                        };
                    progress.rejected_steps += 1;
                    last = error;
                    step *= shrink;
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
    Ok(AdaptiveNetworkDissipationReport {
        stabilized: StabilizedNetworkDissipationReport {
            transport,
            operator_evaluations: *evaluations - initial_evaluations,
            maximum_stages_used,
            maximum_scaled_parabolic_number: maximum_scaled_parabolic,
        },
        maximum_estimated_local_error: maximum_error,
        temporal_rejections,
    })
}

fn prepare_work<L: TransportLaw>(
    owner: &mut NetworkDissipationOperator<L>,
) -> Result<(), DissipationError> {
    stabilized::prepare_work(owner)?;
    let work = &mut owner.work;
    if work.adaptive.anchor.len() != work.current.len() {
        work.adaptive.anchor = allocated(work.current.len(), 0.0)?;
    }
    if work.adaptive.coarse.len() != work.current.len() {
        work.adaptive.coarse = allocated(work.current.len(), 0.0)?;
    }
    if work.adaptive.coarse_endpoints.len() != work.endpoints_a.len() {
        work.adaptive.coarse_endpoints = allocated(work.endpoints_a.len(), 0.0)?;
    }
    if work.adaptive.first_endpoints.len() != work.endpoints_a.len() {
        work.adaptive.first_endpoints = allocated(work.endpoints_a.len(), 0.0)?;
    }
    Ok(())
}

fn copy(
    destination: &mut [f64],
    source: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DissipationError> {
    for (to, from) in destination.chunks_mut(1024).zip(source.chunks(1024)) {
        poll(cancelled)?;
        to.copy_from_slice(from);
    }
    Ok(())
}

// Second-order local truncation scales with h³. Growth is deliberately bounded
// because state-dependent molecular rates and material domains can change.
fn factor(error: f64, accepted: bool) -> f64 {
    if error == 0.0 {
        2.0
    } else {
        (0.9 / error.cbrt()).clamp(0.2, if accepted { 2.0 } else { 0.8 })
    }
}

fn ratio(coarse: f64, fine: f64, absolute: f64, relative: f64) -> Result<f64, DissipationError> {
    let tolerance = absolute + relative * coarse.abs().max(fine.abs());
    if !tolerance.is_finite() || tolerance < 0.0 || !coarse.is_finite() || !fine.is_finite() {
        return Err(DissipationError::Invalid(
            "non-finite adaptive molecular error scale",
        ));
    }
    let difference = (fine - coarse).abs();
    Ok(if difference == 0.0 {
        0.0
    } else {
        (difference / tolerance) / 3.0
    })
}

fn estimate<L: TransportLaw>(
    owner: &NetworkDissipationOperator<L>,
    controls: AdaptiveNetworkDissipationControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, DissipationError> {
    let work = &owner.work;
    let n = owner.gas.species_count();
    let m = n + 2;
    let relative = controls.relative_tolerance;
    let mut error = 0.0_f64;
    for (index, (coarse, fine)) in work
        .adaptive
        .coarse
        .chunks_exact(m)
        .zip(work.stage.chunks_exact(m))
        .enumerate()
    {
        poll(cancelled)?;
        for k in 0..n {
            error = error.max(ratio(
                coarse[k],
                fine[k],
                controls.absolute_partial_density_tolerance_kg_m3,
                relative,
            )?);
        }
        // The admitted fine endpoint already has original-source primitives.
        // Reconstruct the coarse state through that SAME caloric owner; total
        // formation energy never enters the temperature/velocity error norm.
        let fine_primitive = work.primitive[index];
        let coarse_primitive =
            owner
                .gas
                .decode(coarse, Some(fine_primitive.temperature), cancelled)?;
        error = error.max(ratio(
            coarse_primitive.temperature,
            fine_primitive.temperature,
            controls.absolute_temperature_tolerance_k,
            relative,
        )?);
        error = error.max(ratio(
            coarse_primitive.velocity,
            fine_primitive.velocity,
            controls.absolute_velocity_tolerance_m_s,
            relative,
        )?);
    }
    for (index, shape) in owner.branches.iter().enumerate() {
        for end in 0..2 {
            poll(cancelled)?;
            let local_cell = if end == 0 {
                0
            } else {
                shape.mesh.cell_count() - 1
            };
            let cell = shape.cell_start + local_cell;
            let branch_volume = shape.mesh.volumes()[local_cell];
            let mut volume = branch_volume;
            let momentum_capacity = branch_volume * work.primitive[cell].density;
            let mut heat_capacity = branch_volume * work.heat_capacity_v[cell];
            if let Some(connection) = owner.endpoint_connections[index][end] {
                let plenum = owner.connections[connection].plenum;
                volume = volume.min(owner.plenum_volumes[plenum]);
                heat_capacity = heat_capacity.min(
                    owner.plenum_volumes[plenum]
                        * work.heat_capacity_v[owner.plenum_cell_start + plenum],
                );
            }
            // Convert physical field accuracy to inventory accuracy using the
            // actual receiving storage. Relative port scales use the transferred
            // quantity, so opposite ports cannot conceal quadrature error.
            for k in 0..m {
                let entry = (2 * index + end) * m + k;
                let coarse = work.adaptive.coarse_endpoints[entry];
                let fine =
                    work.adaptive.first_endpoints[entry] + work.stabilized_endpoint_changes[entry];
                let absolute = if k < n {
                    volume * controls.absolute_partial_density_tolerance_kg_m3
                } else if k == n {
                    momentum_capacity * controls.absolute_velocity_tolerance_m_s
                } else {
                    heat_capacity * controls.absolute_temperature_tolerance_k
                };
                error = error.max(ratio(coarse, fine, absolute, relative)?);
            }
        }
    }
    Ok(error)
}
