//! Two native implicit stages on the same fixed-mass enthalpy network.
//!
//! With gamma = 1 - 1/sqrt(2), the stage equations in cell joules are
//!   H1 - Hn = gamma*dt*Q(H1),
//!   H2 - Hn = (1-gamma)*dt*Q(H1) + gamma*dt*Q(H2).
//! Q is the original conservative heat-rate owner. The first-stage rate is
//! retained directly, never inferred from a rounded stored-energy difference.
//! H2 is the accepted endpoint; H1 does not become the physical old state of
//! the second solve. Both stages share the original material/domain, Newton,
//! sparse/direct linear solve, residual, entropy and cancellation gates.
use super::*;

const GAMMA: f64 = 1.0 - core::f64::consts::FRAC_1_SQRT_2;

// The first-stage weight GAMMA*dt and time of a step from `start_s` to `end_s`,
// or None when f64 has no stage time strictly inside the step. The step solve
// refuses such a step; the fixed-step drivers use the same test to merge an
// unusable remainder into the step before it.
fn sdirk2_stage(start_s: f64, end_s: f64) -> Result<Option<(f64, f64)>, DepositionError> {
    let implicit_dt_s = finite(GAMMA * finite(end_s - start_s)?)?;
    let stage_time_s = finite(start_s + implicit_dt_s)?;
    let inside = positive(implicit_dt_s) && start_s < stage_time_s && stage_time_s < end_s;
    Ok(inside.then_some((implicit_dt_s, stage_time_s)))
}

// The end of a fixed SDIRK2 step of at most `dt_s` from `start_s` toward `stop_s`.
// When the remainder after `start_s + dt_s` cannot hold a stage (sdirk2_stage),
// the step ends at `stop_s` instead. That remainder r is below
// U / (GAMMA * (1 - 2^-53)), U the larger ulp of the remainder's two ends: at
// that length GAMMA*r rounds to at least one ulp past its start and (1 - GAMMA)*r
// to at least one ulp before its end, so a stage time exists. The last step thus
// exceeds `dt_s` by less than that bound plus the rounding of `start_s + dt_s`.
pub(super) fn sdirk2_step_end(
    start_s: f64,
    dt_s: f64,
    stop_s: f64,
) -> Result<f64, DepositionError> {
    if dt_s == stop_s - start_s {
        return Ok(stop_s);
    }
    let end = start_s + dt_s;
    if end.is_finite() && end < stop_s && sdirk2_stage(end, stop_s)?.is_none() {
        return Ok(stop_s);
    }
    Ok(end)
}

#[path = "sdirk/adaptive.rs"]
mod time_control;

#[path = "sdirk/entropy.rs"]
mod entropy;

impl ImplicitConductionChain<'_> {
    // Check the complete requested window BEFORE an adaptive controller can
    // process events. Source/bath/native-Robin jumps preserve prepared bindings;
    // a later geometry or law-owned coefficient change must not admit a prefix.
    pub(super) fn admit_sdirk2_window(
        &self,
        thermal: &DepositionSimulation<'_>,
        target_s: f64,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitConductionError> {
        // A prescribed flow is sampled at the actual time of EACH stage.
        // Its original endpoint owner admits continuity/calorics at EACH stage,
        // and StageHeat retains its paired transfers and external stream power
        // with the same quadrature as conduction. No flow history is advanced.
        for (index, event) in thermal.schedule[thermal.cursor..]
            .iter()
            .take_while(|event| event.time_s <= target_s)
            .enumerate()
        {
            // The caller has already polled at entry. Keep early unsupported
            // events a preflight refusal while bounding long admitted schedules.
            if index != 0 && index % 256 == 0 {
                checkpoint()?;
            }
            let Event::Boundary(id, boundary) = &event.event else {
                return Err(ImplicitConductionError::Invalid(
                    "SDIRK2 admits only scheduled boundary changes on fixed cells",
                ));
            };
            let cell = thermal
                .cells
                .get(id)
                .ok_or(ImplicitConductionError::Invalid(
                    "SDIRK2 boundary event requires an active cell",
                ))?;
            crate::enthalpy::validate_boundary(*boundary)?;
            let mut expected = cell.boundary;
            expected.source_w = boundary.source_w;
            expected.ambient_k = boundary.ambient_k;
            // An unbound ambient channel uses area*htc in both the native rate
            // and tangent. Bound channels own their complete conductance instead;
            // admitting an authored htc change there would silently ignore it.
            if !self.conduction.ambient_faces.contains_key(id) {
                expected.htc_w_m2_k = boundary.htc_w_m2_k;
            }
            if let (Some(expected), Some(reservoir)) =
                (expected.reservoir.as_mut(), boundary.reservoir)
            {
                expected.temperature_k = reservoir.temperature_k;
            }
            if *boundary != expected {
                return Err(ImplicitConductionError::Invalid(
                    "SDIRK2 boundary events must preserve geometry and law-owned coefficients",
                ));
            }
        }
        Ok(())
    }

    fn solve_sdirk2_step<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: ImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        report: &mut ImplicitConductionReport,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, EnergyAudit, EntropyAudit, f64), ImplicitConductionError>
    {
        self.solve_sdirk2_step_with_sources(
            thermal, target_s, &SourcePowers::default(), &[], controls, workspace, report, checkpoint,
        )
    }

    // Observer-owned start powers, whole-interval means and transfers are held
    // through both native stages. Only the enclosing observed driver prepares
    // or accepts physical source history; an internal stage is not a substep.
    pub(super) fn solve_sdirk2_step_with_sources<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        sources: &SourcePowers,
        transfers: &[HeatTransfer],
        controls: ImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        report: &mut ImplicitConductionReport,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, EnergyAudit, EntropyAudit, f64), ImplicitConductionError> {
        let dt = finite(target_s - thermal.time_s)?;
        let Some((implicit_dt_s, stage_time_s)) = sdirk2_stage(thermal.time_s, target_s)? else {
            return Err(DepositionError::Arithmetic.into());
        };
        // Generated contributions add to each cell's authored power.
        // Native StageHeat carries them and any attached endpoint advection with
        // the same quadrature as conduction, including external stream energy
        // and unclassified at-cell entropy. Reconstructed faces and their
        // tangents are evaluated at each trial, never frozen at the first stage.
        let first = self.solve_stage(
            thermal,
            stage_time_s,
            sources,
            transfers,
            Some(StageQuadrature {
                implicit_dt_s,
                earlier: None,
                retain_heat: true,
            }),
            controls,
            workspace,
            report,
            checkpoint,
        )?;
        let first_maximum = first.maximum_cell_residual_j();
        // The passive Clausius check has complete thermal-bath semantics only
        // in this scope. General two-bath callbacks do not assert passivity;
        // source/work and advective entropy require their own closures.
        let passive = self.advection.is_none()
            && self.heat_sources.is_empty()
            && self.coupled_heat_source.is_none()
            && self.reservoir_exchanges.is_empty()
            && self.ambient_exchanges.is_empty()
            && first
                .entropy_audit()
                .is_some_and(|audit| !audit.has_unclassified_sources);
        let (_, first_endpoint, _, heat) = first.accept_with_stage_heat(
            controls.cell_tolerance_j,
            controls.total_tolerance_j,
            &mut *checkpoint,
        )?;
        // Retain already evaluated coefficients only in the admitted scope.
        // They are used for an uncertainty reduction only if entropy decreases;
        // no additional constitutive solve or caloric integral is introduced.
        let first_endpoint = passive.then_some(first_endpoint);
        let carried = heat.scaled((1.0 - GAMMA) / GAMMA, &mut *checkpoint)?;
        let second = self.solve_stage(
            thermal,
            target_s,
            sources,
            transfers,
            Some(StageQuadrature {
                implicit_dt_s,
                earlier: Some(&carried),
                retain_heat: false,
            }),
            controls,
            workspace,
            report,
            checkpoint,
        )?;
        if let Some(first_endpoint) = &first_endpoint {
            self.admit_passive_sdirk2_entropy(
                thermal,
                first_endpoint,
                &second,
                dt,
                controls.cell_tolerance_j,
                checkpoint,
            )?;
        }
        let (next, energy, entropy, maximum) = accept_trial(second, controls, checkpoint)?;
        Ok((next, energy, entropy, first_maximum.max(maximum)))
    }

    /// Advance fixed-mass convection/conduction with second-order, L-stable SDIRK2.
    ///
    /// All geometric, nonlinear material, contact and integrated boundary laws
    /// of this prepared solver remain active. Authored cell source/sink powers,
    /// ambient/reservoir temperatures and unbound ambient Robin htc may change
    /// through boundary events. Every interval ends at the next event and uses
    /// its original left-limit boundary values in both stages. Events at the
    /// target update the returned boundary without heat under its new values.
    /// Bound ambient htc and reservoir conductance declarations must stay fixed:
    /// their geometric/complete-law owners supply conductance instead. Reservoir
    /// presence/area and all prepared geometry/material bindings stay fixed.
    /// Both stages integrate powers in the original energy residual. Attached
    /// local source laws receive each stage's actual physical time and retain
    /// their full enthalpy feedback; their declared continuous slope changes
    /// cap every interval exactly. Source entropy retains its original
    /// unclassified at-cell convention.
    /// Attached prescribed implicit advection remains active, including limited
    /// reconstruction on the graph path. Both stages transport their actual
    /// endpoint enthalpy with the original flow, continuity and caloric admission;
    /// rates/inlet enthalpies are fixed unless a continuous advection schedule
    /// is attached; scheduled values use each stage's actual time.
    /// External observer histories use [`Self::advance_sdirk2_observed_to`].
    /// Moving mass and all other scheduled events through the target refuse
    /// before stage work. Boundary events retain
    /// authored order and the same whole-call nonlinear and accepted-step budgets.
    ///
    /// Both stages use the SAME original enthalpy storage. Internal first-stage
    /// heat and advective pairs cancel from the external account. Bath heat,
    /// inlet/outlet enthalpy, bath entropy, conductive production and transport
    /// entropy use the two positive Runge--Kutta weights. The optional
    /// `with_implicit_advection_entropy` binding supplies original upwind boundary
    /// entropy and spatial mixing; otherwise transport Q/T stays unclassified.
    /// Exact caloric storage and the signed temporal defect remain separate;
    /// no chemical entropy or constitutive fluid dissipation is inferred.
    ///
    /// Passive networks without advection, source laws, unclassified power or
    /// generalized two-bath callbacks additionally admit the exact cell entropy
    /// change minus the actual bath entropy input. A decrease beyond the local
    /// joule-residual and summation allowance retries the whole interval. This
    /// includes passive single-bath laws in their existing thermal-bath
    /// convention; it does not infer photon entropy or change their closure.
    ///
    /// Each accepted physical interval counts once in `energy.substeps`; all
    /// stage/retry work consumes the existing residual and linear budgets.
    /// Failed intervals may be halved for convergence or passive entropy
    /// admission, without publishing either stage. This is not adaptive
    /// time-error control.
    /// Second-order accuracy requires a smooth admitted solution and independent
    /// time refinement within each event interval; constitutive/phase/limiter
    /// corners can reduce the observed order.
    /// No explicit flow CFL is imposed; spatial transport error still needs
    /// independent refinement and the supplied flow's history remains unchanged.
    /// L-stability does not guarantee positivity or a maximum principle for
    /// arbitrarily large intervals. Original constitutive domains are enforced
    /// at every stage; temperatures and heat are never clipped.
    pub fn advance_sdirk2_to<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: ImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, ImplicitConductionReport), ImplicitConductionError> {
        checkpoint()?;
        validate_advance(
            thermal.time_s,
            target_s,
            self.ids.len(),
            controls,
            workspace.cells(),
        )?;
        self.admit(thermal, &mut checkpoint)?;
        self.admit_sdirk2_window(thermal, target_s, &mut checkpoint)?;
        let mut next = thermal.clone();
        let mut observer = None;
        let mut report = empty_report();
        let mut supplied = Sum::default();
        // Residual magnitudes of the accepted steps, for the step allowances.
        let mut spent = Sum::default();
        loop {
            checkpoint()?;
            // Reuse the native event owner only when events are due. The common
            // constant-boundary path does not clone/process a timeline per step.
            if next
                .schedule
                .get(next.cursor)
                .is_some_and(|event| event.time_s <= next.time_s)
            {
                let current_observer = observer.get_or_insert_with(|| self.conduction.clone());
                let (with_events, observed, events) = next.advance_observed(
                    current_observer,
                    next.time_s,
                    StepConfig {
                        max_dt_s: controls.maximum_dt_s,
                        energy_tolerance_j: controls.total_tolerance_j,
                    },
                    &mut checkpoint,
                )?;
                next = with_events;
                observer = Some(observed);
                report.energy.events += events.events;
                self.admit(&next, &mut checkpoint)?;
            }
            if next.time_s == target_s {
                break;
            }
            if report.energy.substeps >= controls.maximum_steps.min(next.budget.max_substeps) {
                return Err(DepositionError::Budget("implicit thermal intervals").into());
            }
            let start = next.time_s;
            let mut stop = next
                .schedule
                .get(next.cursor)
                .map_or(target_s, |event| event.time_s.min(target_s));
            if let Some(transition) = self.next_forcing_transition_time_s(start) {
                if !transition.is_finite() || transition <= start {
                    return Err(ImplicitConductionError::Invalid(
                        "source transition must be finite and strictly future",
                    ));
                }
                stop = stop.min(transition);
            }
            let remaining = finite(stop - start)?;
            let mut dt = controls.maximum_dt_s.min(remaining);
            let mut halvings = 0;
            loop {
                checkpoint()?;
                let end = sdirk2_step_end(start, dt, stop)?;
                if !end.is_finite() || end <= start {
                    return Err(DepositionError::Arithmetic.into());
                }
                match self.solve_sdirk2_step(
                    &next,
                    end,
                    controls.interval_allowance(thermal.time_s, target_s, end, spent.total()?)?,
                    workspace,
                    &mut report,
                    &mut checkpoint,
                ) {
                    Ok((accepted, energy, entropy, maximum)) => {
                        supplied.add(energy.external_heat_j)?;
                        spent.add(energy.residual_j.abs())?;
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
                            return Err(error);
                        }
                        halvings += 1;
                        dt *= 0.5;
                    }
                    Err(error) => return Err(error),
                }
            }
        }
        // Independent whole-call storage reduction: accepted stage residuals
        // cannot hide a cumulative energy discrepancy or internal-heat credit.
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
        Ok((next, report))
    }
}

#[cfg(test)]
#[path = "sdirk/tests.rs"]
mod tests;
