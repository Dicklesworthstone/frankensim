//! Joint temporal accuracy around the original coupled BDF1 continuum owner.
//!
//! A coarse interval is compared with two native fine intervals. Only the fine
//! pair contributes physical heat/work and continuation history; no field or
//! material is Richardson-extrapolated. Both ordinary and source-density
//! Boussinesq entry points use the same controller and physical acceptance.
//! These are LOCAL error estimates, not a global trajectory-error certificate.
//!
//! Existing held-input entries keep their original force/boundary semantics.
//! Source-driven entries share this controller with exact mechanical-event
//! stops. Callbacks must be pure, smooth within authored segments and valid on
//! original nonlinear probes. Original thermal inventory events inside either
//! window still refuse; target events remain pending for explicit handling.

use super::drive::LaminarDriveQuery;
use super::*;
use fs_thermal::enthalpy::conduction::implicit::adaptive::ConductionTimeOrder;

pub(super) mod estimate;
use estimate::{Account, compare};

/// First-order comparison spans and shared work/physical accuracy ceilings.
#[derive(Clone, Copy, Debug)]
pub struct LaminarAdaptiveControls {
    pub step: LaminarEnthalpyControls,
    pub initial_step_s: f64,
    /// Bounds full comparison spans. Accepted fine intervals are shorter.
    /// A final target-clipped comparison may also be shorter than this floor.
    pub minimum_step_s: f64,
    /// A one-ULP endpoint gap is rounded to the exact target. All work and
    /// error estimates use those actual spans, including cap-rounding effects.
    pub maximum_step_s: f64,
    /// All attempted native intervals, including coarse/rejected work.
    pub maximum_native_steps: usize,
    /// Retained native fine reports; every accepted comparison adds two.
    /// The original thermal inventory's max_substeps ceiling also applies.
    pub maximum_accepted_steps: usize,
    /// One shared ceiling, including coarse and rejected thermal evaluations.
    /// Each native call receives at most the remaining allowance.
    pub maximum_endpoint_evaluations: usize,
    /// Original per-cell stored-energy error [J], including latent storage.
    pub energy_absolute_tolerance_j: f64,
    pub temperature_absolute_tolerance_k: f64,
    /// Euclidean physical BDM velocity error at every cell vertex [m/s].
    pub velocity_absolute_tolerance_m_s: f64,
    /// Applied separately to integrated bulk/slip/stream/other heat, buoyancy
    /// work, and reported body/prescribed-velocity/open-traction/backflow/
    /// carried-kinetic exchanges. None is deposited as heat by this controller.
    /// Signed and resolved-magnitude defects are checked independently, so
    /// opposing accounts or reversals cannot hide behind net heat/work.
    pub transfer_absolute_tolerance_j: f64,
    /// Scales changes from each accepted start, never formation h, Kelvin or
    /// absolute translation velocity. Must be in [0,1).
    pub relative_tolerance: f64,
    /// Independent whole-window thermal storage/source acceptance [J].
    pub window_tolerance_j: f64,
}

/// Numerical work, separate from retained physical heat and work.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct LaminarAdaptiveWork {
    pub native_steps: usize,
    /// Completed native results, including discarded coarse/rejected results.
    /// On a terminal native error, its unreported work is bounded separately.
    pub endpoint_evaluations: usize,
    pub rejected_comparisons: usize,
}

/// Actual coarse/fine discrepancy converted to the retained fine local error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminarTimeEstimate {
    pub start_s: f64,
    pub end_s: f64,
    pub energy_error_ratio: f64,
    pub temperature_error_ratio: f64,
    pub velocity_error_ratio: f64,
    pub transfer_error_ratio: f64,
    pub maximum_energy_difference_j: f64,
    pub maximum_temperature_difference_k: f64,
    pub maximum_velocity_difference_m_s: f64,
    pub maximum_transfer_difference_j: f64,
    /// Present only when the same controller advances native RCR/inertance.
    pub hydraulic: Option<windkessel::adaptive::LaminarWindkesselTimeEstimate>,
}
impl LaminarTimeEstimate {
    pub fn maximum_error_ratio(&self) -> f64 {
        self.energy_error_ratio
            .max(self.temperature_error_ratio)
            .max(self.velocity_error_ratio)
            .max(self.transfer_error_ratio)
            .max(self.hydraulic.map_or(0.0, |h| h.maximum_error_ratio))
    }
}

/// An original accepted fine interval; the callback payload is never fabricated.
pub struct LaminarAdaptiveStep<T> {
    pub start_s: f64,
    pub end_s: f64,
    pub report: LaminarEnthalpyReport<T>,
}

/// Accepted-only observations plus all numerical work. Buoyancy stays mechanical.
pub struct LaminarAdaptiveReport<T> {
    pub steps: Vec<LaminarAdaptiveStep<T>>,
    pub estimates: Vec<LaminarTimeEstimate>,
    pub work: LaminarAdaptiveWork,
    /// Authored mechanical events reached, including one exactly at the target.
    pub events_reached: usize,
    /// Accepted-only selected mechanical exchanges; numerical loss stays separate.
    pub mechanical: Option<drive::LaminarMechanicalWork>,
    /// Controller proposal for a subsequent window's initial comparison span.
    /// It remains subject to new admission and is not an error guarantee.
    pub suggested_step_s: f64,
    pub stored_energy_change_j: f64,
    pub bulk_viscous_heat_j: f64,
    pub slip_heat_j: f64,
    pub external_stream_enthalpy_j: f64,
    pub other_external_heat_j: f64,
    pub buoyancy_work_j: f64,
    pub thermal_balance_residual_j: f64,
    /// Accepted-only native external-circuit work/storage, never caloric heat.
    pub hydraulic: Option<Vec<windkessel::adaptive::LaminarWindkesselWindowBalance>>,
}

/// Final original states, published only after the WHOLE window succeeds.
/// Intermediate fine reports are retained without cloning full field histories.
pub struct LaminarAdaptiveCandidate<'m, G, T> {
    thermal: DepositionSimulation<'m>,
    flow: MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>>,
    pub(super) report: LaminarAdaptiveReport<T>,
}
impl<'m, G, T> LaminarAdaptiveCandidate<'m, G, T> {
    pub fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub fn flow(&self) -> &MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>> {
        &self.flow
    }
    pub fn report(&self) -> &LaminarAdaptiveReport<T> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'m>,
        MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>>,
        LaminarAdaptiveReport<T>,
    ) {
        (self.thermal, self.flow, self.report)
    }
}

pub type LaminarBoussinesqAdaptiveCandidate<'m, G, T> =
    LaminarAdaptiveCandidate<'m, boussinesq::LaminarBoussinesqProjection<G>, T>;

#[derive(Debug)]
pub enum LaminarAdaptiveFailure<E = DepositionError> {
    Invalid(&'static str),
    /// Source, momentum, nonlinear, thermal and cancellation failures are
    /// terminal, not invitations to hide invalid physics under a smaller dt.
    Native {
        error: LaminarEnthalpyError<E>,
        /// Maximum additional thermal endpoint evaluations spent by this
        /// failed call; the original owner does not return a partial counter.
        endpoint_allowance: usize,
    },
    WorkBudget,
    AccuracyFloor(LaminarTimeEstimate),
    UnrepresentableTime,
}

/// No accepted prefix is returned on failure. The caller's histories stay intact.
#[derive(Debug)]
pub struct LaminarAdaptiveError<E = DepositionError> {
    pub failure: LaminarAdaptiveFailure<E>,
    pub work: LaminarAdaptiveWork,
}
impl<E: core::fmt::Debug> core::fmt::Display for LaminarAdaptiveError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive laminar enthalpy: {:?}", self.failure)
    }
}
impl<E: core::fmt::Debug> std::error::Error for LaminarAdaptiveError<E> {}

impl<'s, 'mesh, 'source, F: FluidTransport + ?Sized> LaminarEnthalpy<'s, 'mesh, 'source, F> {
    /// Advance both original histories through a bounded adaptive window.
    /// Each comparison and each finite-difference probe resolves actual source
    /// viscosity and the original momentum residual. Algebraic tolerances in
    /// controls.step remain independent from temporal error tolerances.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive<'m, G: ThermalFluidTransport, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        input: LaminarEnthalpyInput<'_>,
        controls: LaminarAdaptiveControls,
        resolve_material: impl Fn(MaterialFluxRheologyPoint, f64, &F) -> Result<G, MaterialError>,
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &MaterialNavierStokesRheologyStep<LaminarThermalMaterial<G>>,
        ) -> Result<LaminarHeatExchange<T>, E>,
    ) -> Result<LaminarAdaptiveCandidate<'m, G, T>, LaminarAdaptiveError<E>> {
        self.advance_adaptive_native(
            cx,
            thermal,
            target_s,
            input.history,
            &[],
            controls,
            |old, query, history, step, guess| {
                self.advance(
                    cx,
                    old,
                    query.time_s,
                    LaminarEnthalpyInput { history, ..input },
                    step,
                    guess,
                    &resolve_material,
                    &exchange,
                )
            },
        )
    }

    /// The same joint controller with SAME-endpoint original density buoyancy.
    /// Reference inventory/projection bounds and reduced-pressure conventions
    /// are exactly those of advance_boussinesq; no new EOS is introduced.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_boussinesq_adaptive<'m, G: ThermalFluidTransport, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        input: LaminarEnthalpyInput<'_>,
        config: boussinesq::BoussinesqConfig<2>,
        controls: LaminarAdaptiveControls,
        resolve_material: impl Fn(MaterialFluxRheologyPoint, f64, &F) -> Result<G, MaterialError>,
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &MaterialNavierStokesRheologyStep<
                LaminarThermalMaterial<boussinesq::LaminarBoussinesqProjection<G>>,
            >,
        ) -> Result<LaminarHeatExchange<T>, E>,
    ) -> Result<LaminarBoussinesqAdaptiveCandidate<'m, G, T>, LaminarAdaptiveError<E>> {
        self.advance_adaptive_native(
            cx,
            thermal,
            target_s,
            input.history,
            &[],
            controls,
            |old, query, history, step, guess| {
                self.advance_boussinesq(
                    cx,
                    old,
                    query.time_s,
                    LaminarEnthalpyInput { history, ..input },
                    config,
                    step,
                    guess,
                    &resolve_material,
                    &exchange,
                )
            },
        )
    }

    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub(super) fn advance_adaptive_native<'m, G, T, E>(
        &self,
        cx: &Cx<'_>,
        initial: &DepositionSimulation<'m>,
        target_s: f64,
        history: LaminarEnthalpyHistory<'_>,
        event_times_s: &[f64],
        controls: LaminarAdaptiveControls,
        native: impl Fn(
            &DepositionSimulation<'m>,
            LaminarDriveQuery,
            LaminarEnthalpyHistory<'_>,
            LaminarEnthalpyControls,
            &[f64],
        ) -> Result<LaminarEnthalpyCandidate<'m, G, T>, LaminarEnthalpyError<E>>,
    ) -> Result<LaminarAdaptiveCandidate<'m, G, T>, LaminarAdaptiveError<E>> {
        self.advance_adaptive_with_port_history(
            cx,
            initial,
            target_s,
            history,
            event_times_s,
            controls,
            |old, query, history, _, step, guess| native(old, query, history, step, guess),
            |_, _, _, _, _| Ok(None),
        )
    }

    // The same native controller optionally carries an original external
    // history through the accepted report and checks its additional physics.
    // It owns no separate speculative circuit state and no second time loop.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub(super) fn advance_adaptive_with_port_history<'m, G, T, E>(
        &self,
        cx: &Cx<'_>,
        initial: &DepositionSimulation<'m>,
        target_s: f64,
        history: LaminarEnthalpyHistory<'_>,
        event_times_s: &[f64],
        controls: LaminarAdaptiveControls,
        native: impl Fn(
            &DepositionSimulation<'m>,
            LaminarDriveQuery,
            LaminarEnthalpyHistory<'_>,
            Option<&LaminarEnthalpyReport<T>>,
            LaminarEnthalpyControls,
            &[f64],
        ) -> Result<LaminarEnthalpyCandidate<'m, G, T>, LaminarEnthalpyError<E>>,
        compare_ports: impl Fn(
            Option<&LaminarEnthalpyReport<T>>,
            &LaminarEnthalpyCandidate<'m, G, T>,
            &LaminarEnthalpyCandidate<'m, G, T>,
            &LaminarEnthalpyCandidate<'m, G, T>,
            f64,
        ) -> Result<
            Option<windkessel::adaptive::LaminarWindkesselTimeEstimate>,
            LaminarEnthalpyError<E>,
        >,
    ) -> Result<LaminarAdaptiveCandidate<'m, G, T>, LaminarAdaptiveError<E>> {
        let mut work = LaminarAdaptiveWork::default();
        let wrap = |error| failed_native(error, 0, work);
        checkpoint(cx).map_err(wrap)?;
        if !valid_controls(controls)
            || !target_s.is_finite()
            || !(target_s > initial.time_s())
            || !(target_s - initial.time_s()).is_finite()
            || initial.cells().len() > controls.step.maximum_active_cells
            || history.time_s != initial.time_s()
            || history.velocity_m_s.len() != self.system.n_u
            || history.velocity_m_s.iter().any(|x| !x.is_finite())
        {
            return Err(failed(
                LaminarAdaptiveFailure::Invalid(
                    "finite bounded controls and matching complete histories required",
                ),
                work,
            ));
        }
        for (index, &event) in event_times_s.iter().enumerate() {
            if index % 128 == 0 {
                checkpoint(cx).map_err(wrap)?;
            }
            if !event.is_finite() || event < 0.0 || (index > 0 && event <= event_times_s[index - 1])
            {
                return Err(failed(
                    LaminarAdaptiveFailure::Invalid(
                        "mechanical events must be finite, nonnegative and strictly increasing",
                    ),
                    work,
                ));
            }
        }
        self.binding.check(cx, self.system, initial).map_err(wrap)?;
        // Original inventory/event admission for the ENTIRE window. This is
        // not a coarse solve and does not consume or activate any event.
        BackwardEulerStep::new(initial, target_s, || thermal_checkpoint(cx))
            .map_err(|error| failed_native(error.into(), 0, work))?;
        let guess = vec![0.0; initial.cells().len()];
        let mut current = initial.clone();
        let mut velocity = history.velocity_m_s.to_vec();
        let mut final_flow = None;
        let accepted_limit = controls
            .maximum_accepted_steps
            .min(initial.budget().max_substeps);
        let mut steps: Vec<LaminarAdaptiveStep<T>> = Vec::new();
        let mut estimates = Vec::new();
        let mut account = Account::default();
        let mut proposal = controls.initial_step_s;
        let mut segment = event_times_s.partition_point(|&time| time <= initial.time_s());
        let mut events_reached = 0usize;
        let order = ConductionTimeOrder::First;
        while current.time_s() < target_s {
            checkpoint(cx).map_err(|error| failed_native(error, 0, work))?;
            if controls.maximum_native_steps - work.native_steps < 3
                || accepted_limit.saturating_sub(steps.len()) < 2
            {
                return Err(failed(LaminarAdaptiveFailure::WorkBudget, work));
            }
            let start_s = current.time_s();
            let event = event_times_s.get(segment).copied();
            let stop = event.map_or(target_s, |time| time.min(target_s));
            let remaining = stop - start_s;
            let unclipped_proposal = proposal;
            let mut end_s = if proposal >= remaining {
                stop
            } else {
                (start_s + proposal).min(stop)
            };
            // Preserve the exact event/target, including a one-ULP decimal
            // tail. Every original solve still uses its actual time span.
            if end_s < stop && end_s.next_up() == stop {
                end_s = stop;
            }
            let middle_s = f64::midpoint(start_s, end_s);
            if !(start_s < middle_s && middle_s < end_s) {
                return Err(failed(LaminarAdaptiveFailure::UnrepresentableTime, work));
            }
            let correction = order
                .fine_error_correction(middle_s - start_s, end_s - middle_s)
                .ok_or_else(|| failed(LaminarAdaptiveFailure::UnrepresentableTime, work))?;
            let old = LaminarEnthalpyHistory {
                velocity_m_s: &velocity,
                time_s: start_s,
            };
            let end_query = LaminarDriveQuery {
                time_s: end_s,
                segment,
            };
            let middle_query = LaminarDriveQuery {
                time_s: middle_s,
                segment,
            };
            let previous_report = steps.last().map(|step| &step.report);
            let coarse = native_step(
                &native,
                &current,
                end_query,
                old,
                previous_report,
                controls,
                &guess,
                &mut work,
            )?;
            let first = native_step(
                &native,
                &current,
                middle_query,
                old,
                previous_report,
                controls,
                &guess,
                &mut work,
            )?;
            let fine = native_step(
                &native,
                first.thermal(),
                end_query,
                LaminarEnthalpyHistory {
                    velocity_m_s: first.flow().step.solution.velocity_dofs(),
                    time_s: middle_s,
                },
                Some(first.report()),
                controls,
                &guess,
                &mut work,
            )?;
            let mut estimate = compare(
                cx,
                self.system,
                &current,
                &velocity,
                &coarse,
                &first,
                &fine,
                correction,
                controls,
            )
            .map_err(|error| failed_native(error, 0, work))?;
            estimate.hydraulic = compare_ports(previous_report, &coarse, &first, &fine, correction)
                .map_err(|error| failed_native(error, 0, work))?;
            let error = estimate.maximum_error_ratio();
            if error > 1.0 {
                work.rejected_comparisons += 1;
                let floor = controls.minimum_step_s.min(remaining);
                let span = end_s - start_s;
                if span <= floor {
                    return Err(failed(
                        LaminarAdaptiveFailure::AccuracyFloor(estimate),
                        work,
                    ));
                }
                let reduction = order
                    .rejected_reduction_factor(error)
                    .ok_or_else(|| failed(LaminarAdaptiveFailure::UnrepresentableTime, work))?;
                let reduced = (span * reduction).max(floor);
                if !(reduced < span && start_s + reduced < end_s) {
                    return Err(failed(
                        LaminarAdaptiveFailure::AccuracyFloor(estimate),
                        work,
                    ));
                }
                proposal = reduced;
                continue;
            }
            let (_, _, first_report) = first.into_parts();
            let (next, flow, fine_report) = fine.into_parts();
            account
                .add(&first_report)
                .map_err(|error| failed_native(error, 0, work))?;
            account
                .add(&fine_report)
                .map_err(|error| failed_native(error, 0, work))?;
            steps.push(LaminarAdaptiveStep {
                start_s,
                end_s: middle_s,
                report: first_report,
            });
            steps.push(LaminarAdaptiveStep {
                start_s: middle_s,
                end_s,
                report: fine_report,
            });
            estimates.push(estimate);
            velocity = flow.step.solution.velocity_dofs().to_vec();
            current = next;
            final_flow = Some(flow);
            let growth = order
                .accepted_growth_factor(error)
                .ok_or_else(|| failed(LaminarAdaptiveFailure::UnrepresentableTime, work))?;
            proposal = ((end_s - start_s).min(controls.maximum_step_s / growth) * growth)
                .clamp(controls.minimum_step_s, controls.maximum_step_s);
            if event == Some(end_s) {
                segment += 1;
                events_reached += 1;
                // A short exact event interval must not shrink the nominal
                // step for the following authored continuous segment.
                proposal = proposal
                    .max(unclipped_proposal)
                    .min(controls.maximum_step_s);
            }
        }
        let (values, mechanical) = account
            .finish(cx, initial, &current, controls.window_tolerance_j)
            .map_err(|error| failed_native(error, 0, work))?;
        checkpoint(cx).map_err(|error| failed_native(error, 0, work))?;
        let flow = final_flow.ok_or_else(|| {
            failed(
                LaminarAdaptiveFailure::Invalid("empty adaptive trajectory"),
                work,
            )
        })?;
        Ok(LaminarAdaptiveCandidate {
            thermal: current,
            flow,
            report: LaminarAdaptiveReport {
                steps,
                estimates,
                work,
                events_reached,
                mechanical,
                suggested_step_s: proposal,
                stored_energy_change_j: values[0],
                bulk_viscous_heat_j: values[1],
                slip_heat_j: values[2],
                external_stream_enthalpy_j: values[3],
                other_external_heat_j: values[4],
                buoyancy_work_j: values[5],
                thermal_balance_residual_j: values[6],
                hydraulic: None,
            },
        })
    }
}

fn valid_controls(c: LaminarAdaptiveControls) -> bool {
    [
        c.initial_step_s,
        c.minimum_step_s,
        c.maximum_step_s,
        c.energy_absolute_tolerance_j,
        c.temperature_absolute_tolerance_k,
        c.velocity_absolute_tolerance_m_s,
        c.transfer_absolute_tolerance_j,
        c.window_tolerance_j,
    ]
    .iter()
    .all(|x| x.is_finite() && *x > 0.0)
        && c.minimum_step_s <= c.initial_step_s
        && c.initial_step_s <= c.maximum_step_s
        && c.relative_tolerance.is_finite()
        && (0.0..1.0).contains(&c.relative_tolerance)
        && c.maximum_native_steps >= 3
        && c.maximum_accepted_steps >= 2
        && c.maximum_endpoint_evaluations > 0
        && c.step.maximum_active_cells > 0
        && c.step.maximum_endpoint_evaluations > 0
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn native_step<'m, G, T, E>(
    native: &impl Fn(
        &DepositionSimulation<'m>,
        LaminarDriveQuery,
        LaminarEnthalpyHistory<'_>,
        Option<&LaminarEnthalpyReport<T>>,
        LaminarEnthalpyControls,
        &[f64],
    ) -> Result<LaminarEnthalpyCandidate<'m, G, T>, LaminarEnthalpyError<E>>,
    thermal: &DepositionSimulation<'m>,
    query: LaminarDriveQuery,
    history: LaminarEnthalpyHistory<'_>,
    previous_report: Option<&LaminarEnthalpyReport<T>>,
    controls: LaminarAdaptiveControls,
    guess: &[f64],
    work: &mut LaminarAdaptiveWork,
) -> Result<LaminarEnthalpyCandidate<'m, G, T>, LaminarAdaptiveError<E>> {
    let remaining = controls.maximum_endpoint_evaluations - work.endpoint_evaluations;
    if work.native_steps == controls.maximum_native_steps || remaining == 0 {
        return Err(failed(LaminarAdaptiveFailure::WorkBudget, *work));
    }
    let mut step = controls.step;
    step.maximum_endpoint_evaluations = step.maximum_endpoint_evaluations.min(remaining);
    work.native_steps += 1;
    let solved = native(thermal, query, history, previous_report, step, guess)
        .map_err(|error| failed_native(error, step.maximum_endpoint_evaluations, *work))?;
    let used = solved.report().endpoint_evaluations;
    if used > step.maximum_endpoint_evaluations {
        return Err(failed(
            LaminarAdaptiveFailure::Invalid("native endpoint accounting exceeded its allowance"),
            *work,
        ));
    }
    work.endpoint_evaluations += used;
    Ok(solved)
}

fn failed<E>(
    failure: LaminarAdaptiveFailure<E>,
    work: LaminarAdaptiveWork,
) -> LaminarAdaptiveError<E> {
    LaminarAdaptiveError { failure, work }
}
fn failed_native<E>(
    error: LaminarEnthalpyError<E>,
    endpoint_allowance: usize,
    work: LaminarAdaptiveWork,
) -> LaminarAdaptiveError<E> {
    failed(
        LaminarAdaptiveFailure::Native {
            error,
            endpoint_allowance,
        },
        work,
    )
}
