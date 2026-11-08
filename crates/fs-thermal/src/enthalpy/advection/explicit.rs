//! SSPRK2 on the original explicit enthalpy and prescribed-flow equations.
//!
//! For cell energy E=m*h, stages are E1=E0+dt*Q(t0,E0) and
//! Enew=E0+dt*(Q(t0,E0)+Q(t0+dt,E1))/2. Positive weights preserve the
//! original paired internal transfers and use the same quadrature for every
//! external heat/mass account. Temperature and phase still come from T(h).
//! Both Euler stage states must satisfy the native combined conductive and
//! advective CFL. A failed second-stage bound halves the interval and retries.
//!
//! Only pure Advected<()> and ReconstructedAdvected<()> owners are admitted.
//! Coefficient callbacks are repeatable queries, not physical history advances.
//! Euler entropy histories, arbitrary observer.step histories, momentum,
//! interval-integrated sources and kinetics cannot be averaged by this method.

use super::{Advected, AdvectionAudit, Pending, ReconstructedAdvected, finite, invalid};
use crate::enthalpy::backward_euler::{BoundaryStepCoefficients, StepCoefficients};
use crate::enthalpy::observe::EnergySum;
use crate::enthalpy::{
    DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, HeatRates, StepConfig,
};
use std::collections::BTreeMap;

#[path = "explicit/adaptive.rs"]
pub mod adaptive;
pub use adaptive::{
    AdaptiveSsprk2Controls, AdaptiveSsprk2Error, AdaptiveSsprk2Report, Ssprk2TimeError,
};

/// Physical interval/energy controls plus a shared cap on stage work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ssprk2Controls {
    /// The existing maximum dt and absolute whole-advance joule allowance.
    /// Also choose dt for accuracy and any nonlinear source stiffness; the
    /// native CFL bounds conduction/advection, not an arbitrary source law.
    pub step: StepConfig,
    /// Counts every coefficient/transport stage query, including retries.
    /// Each accepted interval needs two queries; a retry reuses its first stage.
    /// Also caps rejected intervals, including caloric-domain failures before
    /// a second query. Neither kind of speculative work can run unbounded.
    pub maximum_stage_evaluations: usize,
}

/// Accepted physical energy and actual bounded computational work.
#[derive(Clone, Debug, PartialEq)]
pub struct Ssprk2Report {
    /// Original incoming inventory, external heat and stored enthalpy account.
    /// Each accepted physical interval counts once in `substeps`.
    pub energy: EnergyAudit,
    pub stage_evaluations: usize,
    /// Intervals halved for a second-stage CFL or intermediate caloric domain.
    /// Rejected trials contribute no heat, mass throughput or flow history.
    pub rejected_intervals: usize,
}

trait PureFlow: DepositionObserver<Error = DepositionError> {
    fn flow(&self) -> &Advected<()>;
    fn flow_mut(&mut self) -> &mut Advected<()>;
}
impl PureFlow for Advected<()> {
    fn flow(&self) -> &Advected<()> {
        self
    }
    fn flow_mut(&mut self) -> &mut Advected<()> {
        self
    }
}
impl PureFlow for ReconstructedAdvected<()> {
    fn flow(&self) -> &Advected<()> {
        self.advection()
    }
    fn flow_mut(&mut self) -> &mut Advected<()> {
        self.advection_mut()
    }
}

fn empty_energy() -> EnergyAudit {
    EnergyAudit {
        incoming_mass_kg: 0.0,
        incoming_energy_j: 0.0,
        external_heat_j: 0.0,
        stored_energy_change_j: 0.0,
        residual_j: 0.0,
        substeps: 0,
        events: 0,
    }
}

fn unchanged(
    _: &DepositionSimulation<'_>,
) -> Result<BoundaryStepCoefficients<()>, DepositionError> {
    Ok(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    }
    .into())
}

impl Advected<()> {
    /// Advance prescribed donor-cell flow with second-order SSP time sampling.
    ///
    /// Uses native heat/flow equations and the original caloric chart, including
    /// exact latent plateaus. Authored births, contacts and boundary events are
    /// processed only at accepted endpoints. No observer history is advanced
    /// at speculative stages. An attached Euler entropy history refuses.
    /// Smooth solutions have second-order time accuracy; phase/limiter corners
    /// still require refinement. This is CFL-controlled, not time-error adaptive.
    pub fn advance_ssprk2_to<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: Ssprk2Controls,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, Self, Ssprk2Report), DepositionError> {
        self.advance_ssprk2_with_coefficients(thermal, target_s, controls, unchanged, checkpoint)
    }

    /// The same SSP stages with additional source watts and link/bath
    /// conductances supplied through the existing thermal coefficient ingress.
    /// The callback sees each stage's actual time and enthalpy, with authored
    /// boundary values. It must be repeatable and must not advance physical
    /// history: rejected stage queries can be discarded. Place discontinuities
    /// on native thermal events; this callback supplies no hidden event clock.
    pub fn advance_ssprk2_with_coefficients<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: Ssprk2Controls,
        coefficients: impl FnMut(
            &DepositionSimulation<'m>,
        ) -> Result<BoundaryStepCoefficients<()>, DepositionError>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, Self, Ssprk2Report), DepositionError> {
        advance(self, thermal, target_s, controls, coefficients, checkpoint)
    }
}

impl ReconstructedAdvected<()> {
    /// SSPRK2 with the original limited spatial fluxes and optional continuous
    /// pump/inlet schedule. Each stage samples its own physical time; schedule
    /// knots and native events cut intervals exactly. Both stages use the
    /// existing reconstruction-aware combined CFL, including pump startup.
    ///
    /// The semidiscrete thermal-entropy face limiter remains active when bound;
    /// no finite-step entropy claim is inferred from it. Euler entropy-history
    /// bindings refuse. No nonlinear source or history is integrated twice.
    pub fn advance_ssprk2_to<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: Ssprk2Controls,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, Self, Ssprk2Report), DepositionError> {
        self.advance_ssprk2_with_coefficients(thermal, target_s, controls, unchanged, checkpoint)
    }

    /// Stage-local thermal sources/conductances with the same repeatability and
    /// history restrictions as [`Advected::advance_ssprk2_with_coefficients`].
    pub fn advance_ssprk2_with_coefficients<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: Ssprk2Controls,
        coefficients: impl FnMut(
            &DepositionSimulation<'m>,
        ) -> Result<BoundaryStepCoefficients<()>, DepositionError>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, Self, Ssprk2Report), DepositionError> {
        advance(self, thermal, target_s, controls, coefficients, checkpoint)
    }
}

struct Stage<O> {
    rates: HeatRates,
    owner: O,
}

fn evaluate<'m, O: PureFlow>(
    thermal: &mut DepositionSimulation<'m>,
    owner: &O,
    controls: Ssprk2Controls,
    report: &mut Ssprk2Report,
    coefficients: &mut impl FnMut(
        &DepositionSimulation<'m>,
    ) -> Result<BoundaryStepCoefficients<()>, DepositionError>,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<Stage<O>, DepositionError> {
    checkpoint()?;
    if report.stage_evaluations >= controls.maximum_stage_evaluations {
        return Err(DepositionError::Budget("SSPRK2 stage evaluations"));
    }
    report.stage_evaluations += 1;
    let coefficients = coefficients(thermal)?;
    let original_links =
        thermal.set_step_conductances(&coefficients.bulk.conductances_w_k, checkpoint)?;
    let result = (|| {
        let mut owner = owner.clone();
        let flow_sources = owner.sources(thermal, checkpoint)?;
        let transfers = owner.heat_transfers(thermal, checkpoint)?;
        let mut rates = thermal.heat_rates_with_boundaries(
            controls.step.max_dt_s,
            &coefficients.boundaries,
            checkpoint,
        )?;
        thermal.add_rate_sources(&mut rates, &coefficients.bulk.sources_w, checkpoint)?;
        thermal.add_rate_sources(&mut rates, &flow_sources, checkpoint)?;
        thermal.add_step_transfers(&mut rates, &transfers, checkpoint)?;
        let dt = owner.limit_timestep_s(thermal, rates.explicit_dt_s, checkpoint)?;
        if !dt.is_finite() || dt <= 0.0 || dt > rates.explicit_dt_s {
            return Err(invalid("invalid SSPRK2 transport timestep bound"));
        }
        rates.explicit_dt_s = dt;
        Ok(Stage { rates, owner })
    })();
    // Trials keep original authored coefficients, including on a retry.
    for (id, conductance) in original_links {
        thermal
            .links
            .get_mut(&id)
            .expect("validated SSP link")
            .conductance_w_k = conductance;
    }
    result
}

fn average_rates(
    first: &HeatRates,
    second: &HeatRates,
    dt: f64,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<HeatRates, DepositionError> {
    let mut watts = BTreeMap::new();
    for (index, (&id, first)) in first.watts.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        let second = second.watts.get(&id).ok_or(DepositionError::Arithmetic)?;
        let mut power = EnergySum::default();
        power.add_scaled(first, 0.5)?;
        power.add_scaled(second, 0.5)?;
        watts.insert(id, power);
    }
    if watts.len() != second.watts.len() {
        return Err(DepositionError::Arithmetic);
    }
    let mut external_w = EnergySum::default();
    external_w.add_scaled(&first.external_w, 0.5)?;
    external_w.add_scaled(&second.external_w, 0.5)?;
    Ok(HeatRates {
        temperatures: BTreeMap::new(),
        watts,
        external_w,
        explicit_dt_s: dt,
    })
}

fn accumulate_flow(
    audit: &mut AdvectionAudit,
    first: &Pending,
    second: &Pending,
    dt: f64,
) -> Result<(), DepositionError> {
    let integrate = |old: f64, a: f64, b: f64| -> Result<f64, DepositionError> {
        let mut sum = EnergySum::default();
        sum.add(old)?;
        sum.add(finite((0.5 * dt) * a)?)?;
        sum.add(finite((0.5 * dt) * b)?)?;
        sum.value()
    };
    audit.transported_mass_kg =
        integrate(audit.transported_mass_kg, first.mass_kg_s, second.mass_kg_s)?;
    audit.absolute_enthalpy_transfer_j = integrate(
        audit.absolute_enthalpy_transfer_j,
        first.absolute_power_w,
        second.absolute_power_w,
    )?;
    audit.inlet_mass_kg = integrate(audit.inlet_mass_kg, first.inlet_kg_s, second.inlet_kg_s)?;
    audit.outlet_mass_kg = integrate(audit.outlet_mass_kg, first.outlet_kg_s, second.outlet_kg_s)?;
    audit.external_enthalpy_j = integrate(
        audit.external_enthalpy_j,
        first.external_power_w,
        second.external_power_w,
    )?;
    audit.absolute_mass_balance_roundoff_kg = integrate(
        audit.absolute_mass_balance_roundoff_kg,
        first.mass_balance_roundoff_kg_s,
        second.mass_balance_roundoff_kg_s,
    )?;
    audit.absolute_datum_energy_roundoff_j = integrate(
        audit.absolute_datum_energy_roundoff_j,
        first.datum_energy_roundoff_w,
        second.datum_energy_roundoff_w,
    )?;
    audit.maximum_courant = audit
        .maximum_courant
        .max(finite(dt * first.decay_per_s)?)
        .max(finite(dt * second.decay_per_s)?);
    audit.substeps = audit
        .substeps
        .checked_add(1)
        .ok_or(DepositionError::Arithmetic)?;
    Ok(())
}

fn halve(dt: &mut f64, report: &mut Ssprk2Report, maximum: usize) -> Result<(), DepositionError> {
    if report.rejected_intervals >= maximum {
        return Err(DepositionError::Budget("SSPRK2 rejected intervals"));
    }
    let shorter = *dt * 0.5;
    if shorter <= 0.0 || shorter >= *dt {
        return Err(DepositionError::Arithmetic);
    }
    *dt = shorter;
    report.rejected_intervals = report
        .rejected_intervals
        .checked_add(1)
        .ok_or(DepositionError::Arithmetic)?;
    Ok(())
}

struct Interval<'m, O> {
    thermal: DepositionSimulation<'m>,
    owner: O,
    external_heat: EnergySum,
    unrepresented_heat_j: f64,
}

/// Exactly one SSPRK2 interval, without hidden subcycling or event execution.
/// A speculative caloric-domain or second-stage CFL failure asks the caller
/// to shorten the entire interval; all queried stages still consume its budget.
fn try_interval<'m, O: PureFlow>(
    initial: &DepositionSimulation<'m>,
    first: &mut Stage<O>,
    end: f64,
    dt: f64,
    unrepresented_heat_j: f64,
    controls: Ssprk2Controls,
    report: &mut Ssprk2Report,
    coefficients: &mut impl FnMut(
        &DepositionSimulation<'m>,
    ) -> Result<BoundaryStepCoefficients<()>, DepositionError>,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<Option<Interval<'m, O>>, DepositionError> {
    let result = (|| {
        first.rates.explicit_dt_s = dt;
        let mut predictor = initial.clone();
        let mut trial_roundoff = 0.0;
        predictor.thermal_step_with_rates(
            &first.rates,
            &mut empty_energy(),
            &mut trial_roundoff,
            None,
            checkpoint,
        )?;
        if trial_roundoff > controls.step.energy_tolerance_j {
            return Err(DepositionError::Arithmetic);
        }
        predictor.time_s = end;
        let mut second = evaluate(
            &mut predictor,
            &first.owner,
            controls,
            report,
            coefficients,
            checkpoint,
        )?;
        if dt > second.rates.explicit_dt_s {
            return Ok(None);
        }
        // Admit both Euler endpoints before mixing. A valid average must not
        // conceal an out-of-chart SSP stage; reuse the predictor allocation.
        second.rates.explicit_dt_s = dt;
        predictor.thermal_step_with_rates(
            &second.rates,
            &mut empty_energy(),
            &mut trial_roundoff,
            None,
            checkpoint,
        )?;
        if trial_roundoff > controls.step.energy_tolerance_j {
            return Err(DepositionError::Arithmetic);
        }
        let rates = average_rates(&first.rates, &second.rates, dt, checkpoint)?;
        let mut accepted = initial.clone();
        let mut accepted_roundoff = unrepresented_heat_j;
        accepted.thermal_step_with_rates(
            &rates,
            &mut empty_energy(),
            &mut accepted_roundoff,
            None,
            checkpoint,
        )?;
        if accepted_roundoff > controls.step.energy_tolerance_j {
            return Err(DepositionError::Arithmetic);
        }
        let mut flow_audit = first.owner.flow().audit;
        accumulate_flow(
            &mut flow_audit,
            first
                .owner
                .flow()
                .pending
                .as_ref()
                .ok_or(DepositionError::Arithmetic)?,
            second
                .owner
                .flow()
                .pending
                .as_ref()
                .ok_or(DepositionError::Arithmetic)?,
            dt,
        )?;
        let mut external_heat = EnergySum::default();
        external_heat.add_scaled(&rates.external_w, dt)?;
        accepted.time_s = end;
        let mut owner = second.owner;
        let flow = owner.flow_mut();
        flow.audit = flow_audit;
        flow.pending = None;
        flow.boundary_pending = None;
        Ok(Some(Interval {
            thermal: accepted,
            owner,
            external_heat,
            unrepresented_heat_j: accepted_roundoff,
        }))
    })();
    match result {
        Err(DepositionError::MaterialDomain) => Ok(None),
        other => other,
    }
}

fn advance<'m, O: PureFlow>(
    owner: &O,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    controls: Ssprk2Controls,
    mut coefficients: impl FnMut(
        &DepositionSimulation<'m>,
    ) -> Result<BoundaryStepCoefficients<()>, DepositionError>,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<(DepositionSimulation<'m>, O, Ssprk2Report), DepositionError> {
    checkpoint()?;
    if !target_s.is_finite()
        || target_s < thermal.time_s()
        || !controls.step.max_dt_s.is_finite()
        || controls.step.max_dt_s <= 0.0
        || !controls.step.energy_tolerance_j.is_finite()
        || controls.step.energy_tolerance_j < 0.0
    {
        return Err(invalid(
            "invalid SSPRK2 target, timestep or energy allowance",
        ));
    }
    if owner.flow().entropy.is_some() {
        return Err(invalid(
            "SSPRK2 requires its own entropy history; Euler flow entropy is unsupported",
        ));
    }
    finite(thermal.energy_j())?;
    let mut next = thermal.clone();
    let mut owner = owner.clone();
    owner.begin()?;
    let mut report = Ssprk2Report {
        energy: empty_energy(),
        stage_evaluations: 0,
        rejected_intervals: 0,
    };
    let mut external_heat = EnergySum::default();
    let mut unrepresented_heat = 0.0;
    loop {
        checkpoint()?;
        while next.cursor < next.schedule.len() && next.schedule[next.cursor].time_s == next.time_s
        {
            checkpoint()?;
            let event = next.schedule[next.cursor].event.clone();
            next.apply(event.clone(), &mut report.energy)?;
            next.cursor += 1;
            report.energy.events += 1;
            owner.event(&event, &next)?;
        }
        if next.time_s == target_s {
            break;
        }
        let mut stop = next
            .schedule
            .get(next.cursor)
            .map_or(target_s, |event| event.time_s.min(target_s));
        if let Some(transition) = owner.next_transition_time_s(&next, &mut checkpoint)? {
            if !transition.is_finite() || transition <= next.time_s {
                return Err(invalid("invalid SSPRK2 source transition"));
            }
            stop = stop.min(transition);
        }
        if next.cells.is_empty() {
            next.time_s = stop;
            continue;
        }
        if report.energy.substeps >= next.budget.max_substeps {
            return Err(DepositionError::Budget("SSPRK2 accepted intervals"));
        }
        let mut first = evaluate(
            &mut next,
            &owner,
            controls,
            &mut report,
            &mut coefficients,
            &mut checkpoint,
        )?;
        let start = next.time_s;
        let mut dt = first.rates.explicit_dt_s.min(stop - start);
        loop {
            checkpoint()?;
            let end = if dt == stop - start {
                stop
            } else {
                finite(start + dt)?
            };
            if end <= start {
                return Err(DepositionError::Arithmetic);
            }
            let Some(accepted) = try_interval(
                &next,
                &mut first,
                end,
                dt,
                unrepresented_heat,
                controls,
                &mut report,
                &mut coefficients,
                &mut checkpoint,
            )?
            else {
                halve(&mut dt, &mut report, controls.maximum_stage_evaluations)?;
                continue;
            };
            external_heat.add_scaled(&accepted.external_heat, 1.0)?;
            next = accepted.thermal;
            owner = accepted.owner;
            unrepresented_heat = accepted.unrepresented_heat_j;
            report.energy.substeps += 1;
            break;
        }
    }
    thermal.admit_advance_energy(
        &next,
        &external_heat,
        &mut report.energy,
        controls.step.energy_tolerance_j,
        &mut checkpoint,
    )?;
    Ok((next, owner, report))
}

#[cfg(test)]
#[path = "explicit/tests.rs"]
mod tests;
