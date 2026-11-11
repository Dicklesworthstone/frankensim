//! Prescribed/lagged sources and transport on the original implicit heat solver.
use super::*;

/// Native implicit refusal or an unchanged source-owner refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum ImplicitConductionObservedError<E> {
    Implicit(ImplicitConductionError),
    Observer(E),
}
impl<E> From<ImplicitConductionError> for ImplicitConductionObservedError<E> {
    fn from(error: ImplicitConductionError) -> Self {
        Self::Implicit(error)
    }
}
impl<E> From<DepositionError> for ImplicitConductionObservedError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Implicit(error.into())
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for ImplicitConductionObservedError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "observed implicit conduction: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for ImplicitConductionObservedError<E> {}

// Reuse the original event timeline without restarting the source work budget
// at every accepted time. This wrapper is only advanced to its current time.
#[derive(Clone)]
struct Events<'a, O> {
    conduction: ConductionObserver<'a>,
    source: O,
}
impl<O: DepositionObserver> DepositionObserver for Events<'_, O> {
    type Error = ImplicitConductionObservedError<O::Error>;

    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.conduction.event(event, thermal)?;
        self.source
            .event(event, thermal)
            .map_err(ImplicitConductionObservedError::Observer)
    }

    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        Err(ImplicitConductionError::Invalid("event-only implicit observer cannot step").into())
    }
}

fn restore_sources(thermal: &mut DepositionSimulation<'_>, authored: Vec<(CellId, f64)>) {
    for (id, watts) in authored {
        thermal
            .cells
            .get_mut(&id)
            .expect("validated active source")
            .boundary
            .source_w = watts;
    }
}

#[derive(Clone, Copy)]
enum ObservedMethod {
    BackwardEuler,
    Sdirk2,
}

impl ImplicitConductionChain<'_> {
    /// Advance native conduction/cooling implicitly with prescribed or lagged
    /// external heat and conservative internal transfers. The original source
    /// hooks see the accepted START state; their watts are held fixed during
    /// each endpoint Newton solve. Fixed
    /// absorption responses therefore need no optical Jacobian. A state-dependent
    /// optical source is explicit in this IMEX discretization and still needs
    /// its own timestep restriction and independent temporal refinement.
    ///
    /// Source transitions and native events cut intervals exactly. Interval
    /// means are recomputed after nonlinear halving, using the actual native
    /// duration. `begin` runs once. Preparation hooks may replace pending state
    /// and caches but must publish physical history only in `step`; repeated
    /// preparations at the same accepted time must retain cumulative work costs.
    /// A rejected interval never calls `step`, and any final refusal returns no
    /// thermal or observer candidate. Clone must retain the original observer's
    /// transactional, non-shared mutable history contract.
    ///
    /// Balanced mass-carried enthalpy uses the original `Advected` observer,
    /// including its datum, storage-basis, flow and timestep admission. Flow is
    /// prepared once per accepted start state and retained through nonlinear
    /// retries. This is explicit advection with implicit conduction; the flow
    /// restriction and independent temporal refinement still apply.
    ///
    /// Nonempty conductance, ambient-conductance or reservoir-conductance hooks
    /// refuse explicitly; all conduction remains owned by this solver. Native
    /// audits count external heat once and cancel internal transfers. Transfer
    /// entropy remains unclassified unless its owner supplies its own account.
    /// This is not the adaptive step-doubling path or a moving-mass model.
    /// Independent start/interval contributions remain compensated alongside
    /// authored and endpoint sources through the original heat, residual and
    /// entropy owners. The native local and whole-call energy gates account
    /// for the actual stored heat without a scalar source handoff.
    pub fn advance_observed_to<'m, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'m>,
        observer: &O,
        target_s: f64,
        controls: ImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (DepositionSimulation<'m>, O, ImplicitConductionReport),
        ImplicitConductionObservedError<O::Error>,
    > {
        self.advance_observed_with_method(
            thermal, observer, target_s, controls, workspace, ObservedMethod::BackwardEuler, checkpoint,
        )
    }

    /// Advance native SDIRK2 conduction with transactional external heat owners.
    ///
    /// The source/transfer hooks see the accepted start state. One preparation
    /// supplies the actual whole-interval mean heat to both native stages; their
    /// positive quadrature weights integrate that heat exactly once. `step` runs
    /// once only after both stages pass the original energy/domain gates. A
    /// retry recomputes interval heat without resetting the source work budget.
    /// No observer history is advanced or published at the internal stage time.
    ///
    /// Source timestep restrictions, exact transitions, authored boundary events
    /// and unsupported conductance hooks follow [`Self::advance_observed_to`].
    /// The complete window must satisfy [`Self::advance_sdirk2_to`]'s fixed-cell
    /// event admission before any observer work. The same native stage solve,
    /// entropy accounting and sparse/direct linear owner remain authoritative.
    ///
    /// SDIRK2 is second order for its admitted thermal endpoint laws. Holding
    /// state-dependent observed forcing or explicit transfers at the start can
    /// reduce the coupled method to first order; retain the source owner's CFL
    /// and independently refine the coupled time grid. An interval heat mean
    /// does not supply a source tangent or imply strong coupling.
    pub fn advance_sdirk2_observed_to<'m, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'m>,
        observer: &O,
        target_s: f64,
        controls: ImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (DepositionSimulation<'m>, O, ImplicitConductionReport),
        ImplicitConductionObservedError<O::Error>,
    > {
        self.advance_observed_with_method(
            thermal, observer, target_s, controls, workspace, ObservedMethod::Sdirk2, checkpoint,
        )
    }

    fn advance_observed_with_method<'m, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'m>,
        observer: &O,
        target_s: f64,
        controls: ImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        method: ObservedMethod,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (DepositionSimulation<'m>, O, ImplicitConductionReport),
        ImplicitConductionObservedError<O::Error>,
    > {
        use ImplicitConductionObservedError::Observer;
        checkpoint()?;
        validate_advance(
            thermal.time_s,
            target_s,
            self.ids.len(),
            controls,
            workspace.cells(),
        )?;
        self.admit(thermal, &mut checkpoint)?;
        if matches!(method, ObservedMethod::Sdirk2) {
            self.admit_sdirk2_window(thermal, target_s, &mut checkpoint)?;
        }
        let mut next = thermal.clone();
        let mut observers = Events {
            conduction: self.conduction.clone(),
            source: observer.clone(),
        };
        observers.source.begin().map_err(Observer)?;
        let mut report = empty_report();
        let mut supplied = Sum::default();
        // Residual magnitudes of the accepted intervals, for their allowances.
        let mut spent = Sum::default();
        loop {
            checkpoint()?;
            let (with_events, observed, events) = next.advance_observed(
                &observers,
                next.time_s,
                StepConfig {
                    max_dt_s: controls.maximum_dt_s,
                    energy_tolerance_j: controls.total_tolerance_j,
                },
                &mut checkpoint,
            )?;
            next = with_events;
            observers = observed;
            report.energy.events += events.events;
            self.admit(&next, &mut checkpoint)?;
            if next.time_s == target_s {
                break;
            }
            if report.energy.substeps >= controls.maximum_steps.min(next.budget.max_substeps) {
                return Err(DepositionError::Budget("implicit thermal intervals").into());
            }
            let start = next.time_s;
            let source = &mut observers.source;
            let mut stop = next
                .schedule
                .get(next.cursor)
                .map_or(target_s, |e| e.time_s.min(target_s));
            if let Some(transition) = self.next_forcing_transition_time_s(start) {
                if !transition.is_finite() || transition <= start {
                    return Err(ImplicitConductionError::Invalid(
                        "source transition must be finite and strictly future",
                    ).into());
                }
                stop = stop.min(transition);
            }
            if let Some(transition) = source
                .next_transition_time_s(&next, &mut checkpoint)
                .map_err(Observer)?
            {
                if !transition.is_finite() || transition <= start {
                    return Err(DepositionError::InvalidInput(
                        "observer transition must be finite and strictly future",
                    )
                    .into());
                }
                stop = stop.min(transition);
            }
            let sources = crate::enthalpy::SourcePowers::capture(
                &next, source, &mut checkpoint,
            ).map_err(Observer)?;
            if !source
                .conductances(&next, &mut checkpoint)
                .map_err(Observer)?
                .is_empty()
                || !source
                    .ambient_conductances(&next, &mut checkpoint)
                    .map_err(Observer)?
                    .is_empty()
                || !source
                    .reservoir_conductances(&next, &mut checkpoint)
                    .map_err(Observer)?
                    .is_empty()
            {
                return Err(ImplicitConductionError::Invalid(
                    "observed implicit conduction owns all conductive coefficients",
                )
                .into());
            }
            // Preserve the native source -> transfer -> timestep hook order.
            // Advected consumes its source epoch here and retains the pending
            // mass/enthalpy account until exactly one accepted step callback.
            let transfers = source.heat_transfers(&next, &mut checkpoint).map_err(Observer)?;
            let authored = next.add_step_source_powers(&sources, &mut checkpoint)?;
            let limit = source.limit_timestep_s(&next, controls.maximum_dt_s, &mut checkpoint);
            restore_sources(&mut next, authored);
            let limit = limit.map_err(Observer)?;
            if !positive(limit) || limit > controls.maximum_dt_s {
                return Err(
                    DepositionError::InvalidInput("invalid observer timestep bound").into(),
                );
            }
            let mut old_temperatures = BTreeMap::new();
            for &id in &self.ids {
                checkpoint()?;
                old_temperatures.insert(id, next.temperature_k(id)?);
            }
            let mut dt = limit.min(stop - start);
            let mut halvings = 0;
            loop {
                checkpoint()?;
                // Backward Euler has no stage; an SDIRK2 step absorbs a remainder
                // that could not hold one (sdirk::sdirk2_step_end).
                let end = match method {
                    ObservedMethod::BackwardEuler if dt == stop - start => stop,
                    ObservedMethod::BackwardEuler => start + dt,
                    ObservedMethod::Sdirk2 => super::sdirk::sdirk2_step_end(start, dt, stop)?,
                };
                if !end.is_finite() || end <= start {
                    return Err(DepositionError::Arithmetic.into());
                }
                let duration = finite(end - start)?;
                let authored = next.add_step_source_powers(&sources, &mut checkpoint)?;
                let interval = crate::enthalpy::SourcePowers::capture_interval(
                    &next, source, start, duration, end, &mut checkpoint,
                );
                restore_sources(&mut next, authored);
                let interval = interval.map_err(Observer)?;
                let mut combined = sources.clone();
                let mut visited = 0usize;
                interval.visit(&mut |id, watts| {
                    if visited % 256 == 0 { checkpoint()?; }
                    visited += 1;
                    combined.add(&next, id, watts)
                })?;
                let allowance =
                    controls.interval_allowance(thermal.time_s, target_s, end, spent.total()?)?;
                let trial = match method {
                    ObservedMethod::BackwardEuler => self.solve_step_with_sources(
                        &next, end, &combined, &transfers, allowance, workspace,
                        &mut report, &mut checkpoint,
                    ),
                    ObservedMethod::Sdirk2 => self.solve_sdirk2_step_with_sources(
                        &next, end, &combined, &transfers, allowance, workspace,
                        &mut report, &mut checkpoint,
                    ),
                };
                match trial {
                    Ok((accepted, energy, entropy, maximum)) => {
                        spent.add(energy.residual_j.abs())?;
                        source
                            .step(
                                start,
                                duration,
                                &old_temperatures,
                                &accepted,
                                &mut checkpoint,
                            )
                            .map_err(Observer)?;
                        supplied.add(energy.external_heat_j)?;
                        accumulate_entropy(&mut report.entropy, entropy)?;
                        report.maximum_cell_residual_j =
                            report.maximum_cell_residual_j.max(maximum);
                        report.energy.substeps += 1;
                        next = accepted;
                        break;
                    }
                    Err(
                        error @ (ImplicitConductionError::NotConverged { .. }
                        | ImplicitConductionError::EntropyDecrease { .. }
                        | ImplicitConductionError::KrylovNotConverged { .. }
                        | ImplicitConductionError::Linear(TridiagonalError::Singular {
                            ..
                        })),
                    ) => {
                        report.rejected_intervals += 1;
                        if halvings >= controls.maximum_halvings || dt * 0.5 < controls.minimum_dt_s
                        {
                            return Err(error.into());
                        }
                        halvings += 1;
                        dt *= 0.5;
                    }
                    Err(error) => return Err(error.into()),
                }
            }
        }
        let mut storage = Sum::default();
        for &id in &self.ids {
            checkpoint()?;
            storage.add(
                thermal.cells[&id].mass_kg
                    * (next.cells[&id].specific_h_j_kg - thermal.cells[&id].specific_h_j_kg),
            )?;
        }
        report.energy.external_heat_j = supplied.total()?;
        report.energy.stored_energy_change_j = storage.total()?;
        report.energy.residual_j =
            finite(report.energy.stored_energy_change_j - report.energy.external_heat_j)?;
        if report.energy.residual_j.abs() > controls.total_tolerance_j {
            return Err(DepositionError::EnergyBalance(report.energy.residual_j).into());
        }
        checkpoint()?;
        Ok((next, observers.source, report))
    }
}
