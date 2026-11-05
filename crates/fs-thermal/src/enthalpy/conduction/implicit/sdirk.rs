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

#[path = "sdirk/adaptive.rs"]
mod time_control;

impl ImplicitConductionChain<'_> {
    // Check the complete requested window BEFORE an adaptive controller can
    // process events. Source-only jumps preserve the prepared geometry/material
    // bindings; a later bath or topology change must not admit a solved prefix.
    fn admit_sdirk2_window(
        &self,
        thermal: &DepositionSimulation<'_>,
        target_s: f64,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitConductionError> {
        // A prepared implicit flow has fixed mass rates and inlet enthalpies.
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
                    "SDIRK2 admits only scheduled source changes on fixed cells",
                ));
            };
            let cell = thermal
                .cells
                .get(id)
                .ok_or(ImplicitConductionError::Invalid(
                    "SDIRK2 source event requires an active cell",
                ))?;
            crate::enthalpy::validate_boundary(*boundary)?;
            let mut expected = cell.boundary;
            expected.source_w = boundary.source_w;
            if *boundary != expected {
                return Err(ImplicitConductionError::Invalid(
                    "SDIRK2 source events must preserve ambient, film and reservoir boundaries",
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
        let dt = finite(target_s - thermal.time_s)?;
        let implicit_dt_s = finite(GAMMA * dt)?;
        let stage_time_s = finite(thermal.time_s + implicit_dt_s)?;
        if !positive(implicit_dt_s) || stage_time_s <= thermal.time_s || stage_time_s >= target_s {
            return Err(DepositionError::Arithmetic.into());
        }
        // An empty override map preserves each cell's authored, constant power.
        // Native StageHeat carries it and any attached endpoint advection with
        // the same quadrature as conduction, including external stream energy
        // and unclassified at-cell entropy. Reconstructed faces and their
        // tangents are evaluated at each trial, never frozen at the first stage.
        let sources = BTreeMap::new();
        let first = self.solve_stage(
            thermal,
            stage_time_s,
            &sources,
            &[],
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
        let (_, _, _, heat) = first.accept_with_stage_heat(
            controls.cell_tolerance_j,
            controls.total_tolerance_j,
            &mut *checkpoint,
        )?;
        let carried = heat.scaled((1.0 - GAMMA) / GAMMA, &mut *checkpoint)?;
        let second = self.solve_stage(
            thermal,
            target_s,
            &sources,
            &[],
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
        let (next, energy, entropy, maximum) = accept_trial(second, controls, checkpoint)?;
        Ok((next, energy, entropy, first_maximum.max(maximum)))
    }

    /// Advance fixed-mass convection/conduction with second-order, L-stable SDIRK2.
    ///
    /// All geometric, nonlinear material, contact and integrated boundary laws
    /// of this prepared solver remain active. Ambient and reservoir temperatures
    /// must remain constant throughout the call. Authored cell source/sink powers
    /// may change through source-only boundary events: every interval ends at
    /// the next event and uses its original left-limit power in both stages.
    /// Events at the target are applied without aging under their new power.
    /// Both stages integrate powers in the original energy residual;
    /// source entropy retains its original unclassified at-cell convention.
    /// Attached prescribed implicit advection remains active, including limited
    /// reconstruction on the graph path. Both stages transport their actual
    /// endpoint enthalpy with the original flow, continuity and caloric admission;
    /// mass rates and inlet enthalpies remain fixed throughout the call.
    /// External observers, moving mass and all other scheduled events through
    /// the target refuse before stage work. Source events retain
    /// authored order and the same whole-call nonlinear and accepted-step budgets.
    ///
    /// Both stages use the SAME original enthalpy storage. Internal first-stage
    /// heat and advective pairs cancel from the external account. Bath heat,
    /// inlet/outlet enthalpy, bath entropy, conductive production and signed
    /// unclassified transport Q/T use the two positive Runge--Kutta weights.
    /// Exact caloric entropy storage and the signed temporal defect remain
    /// separate; no stream entropy or mixing-production law is inferred.
    ///
    /// Each accepted physical interval counts once in `energy.substeps`; all
    /// stage/retry work consumes the existing residual and linear budgets.
    /// Failed intervals may be halved for nonlinear convergence, without
    /// publishing either stage. This is not adaptive time-error control.
    /// Second-order accuracy requires a smooth admitted solution and independent
    /// time refinement; phase/source/limiter corners can reduce the observed order.
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
        loop {
            checkpoint()?;
            // Reuse the native event owner only when events are due. The common
            // constant-source path does not clone/process a timeline per step.
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
            let stop = next
                .schedule
                .get(next.cursor)
                .map_or(target_s, |event| event.time_s.min(target_s));
            let remaining = finite(stop - start)?;
            let mut dt = controls.maximum_dt_s.min(remaining);
            let mut halvings = 0;
            loop {
                checkpoint()?;
                let end = if dt == remaining { stop } else { start + dt };
                if !end.is_finite() || end <= start {
                    return Err(DepositionError::Arithmetic.into());
                }
                match self.solve_sdirk2_step(
                    &next,
                    end,
                    controls,
                    workspace,
                    &mut report,
                    &mut checkpoint,
                ) {
                    Ok((accepted, energy, entropy, maximum)) => {
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
