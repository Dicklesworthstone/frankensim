//! Transactional substep/event observation of the existing explicit thermal lane.
use super::*;

/// Instantaneous thermal balance from the same source/flux owner used to step.
/// No state is advanced and no boundary energy is counted twice.
#[derive(Clone, Debug, PartialEq)]
pub struct PowerBalance {
    /// Actual caloric-law temperatures at this accepted state [K].
    pub temperatures_k: BTreeMap<CellId, f64>,
    /// Net inward watts for every active cell; a steady solution makes these zero.
    pub net_power_w: BTreeMap<CellId, f64>,
    /// Sum of authored and observed external sources/boundaries [W].
    pub external_power_w: f64,
}

/// History coupling on the actual thermal grid, without storing whole fields.
/// `Clone` must deeply copy mutable physical state. Shared immutable laws are
/// fine; external effects or shared mutable state cannot be rolled back here.
/// An observer is called after each individual event in supplied order, and
/// after each CFL-limited step BEFORE any events at its ending time. An optional
/// source hook and then a conductance hook run before that step, after all
/// events at its starting time. Both see the unchanged start-state history.
pub trait DepositionObserver: Clone {
    type Error: From<DepositionError>;
    /// Start one candidate advance; useful for resetting a per-call work counter.
    fn begin(&mut self) -> Result<(), Self::Error> { Ok(()) }
    /// Additional signed heat rates [W] at active cells, evaluated at the start
    /// of EACH thermal substep, after same-time events. Omitted cells get zero.
    /// These rates are added to authored boundary sources for this step only;
    /// they are never persisted in boundaries or re-added on the next step.
    /// Unknown cells and nonfinite rates refuse before thermal integration.
    ///
    /// The existing conduction CFL does not bound a nonlinear source's stiffness.
    /// Select `max_dt_s` for that source and perform time refinement. A source
    /// discontinuity must coincide with a scheduled event or an advance endpoint.
    /// `step` receives the actual CFL-limited duration for integrating a source's
    /// own work/energy account. Generated heat is already in `external_heat_j`;
    /// do not add it again to that thermal audit. On failure, the source owner's
    /// cloned state and the entire proposed thermal advance are discarded.
    fn sources(&mut self, _thermal: &DepositionSimulation<'_>,
        _checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error> { Ok(BTreeMap::new()) }
    /// Step-local thermal-link conductances [W/K], evaluated after `sources`
    /// against the SAME start state. A source owner may retain one prepared
    /// material evaluation for both hooks. Omitted links keep their authored
    /// conductance; zero explicitly insulates without removing covered area.
    /// Only active link IDs and finite nonnegative values are admitted.
    ///
    /// The original thermal kernel uses these values for BOTH pairwise heat
    /// exchange and its timestep bound. This frozen-coefficient bound is not
    /// a general nonlinear-stability or time-accuracy guarantee. Refine max_dt.
    /// Geometry, per-side coverage, material laws and reference mass never change.
    /// Authored conductances are restored before `step` and scheduled events;
    /// retain applied values in the observer when observations need them.
    fn conductances(&mut self, _thermal: &DepositionSimulation<'_>,
        _checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<u64, f64>, Self::Error> { Ok(BTreeMap::new()) }
    /// Internal pairwise watts at the same start state, after `conductances`.
    /// Unlike `sources`, these are conservative exchanges, not external heat.
    /// At most `DepositionBudget::max_links` transfers may be submitted per step.
    /// They do not cover surfaces or create conductive contacts. The supplying
    /// owner must also restrict the timestep when its exchange law requires it.
    fn heat_transfers(&mut self, _thermal: &DepositionSimulation<'_>,
        _checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<Vec<HeatTransfer>, Self::Error> { Ok(Vec::new()) }
    /// Tighten the conductive/source step bound; return finite 0 < dt <= proposed.
    /// Called after all rate hooks. This sees step-local source/conductance
    /// overrides but unchanged enthalpy. Compose coupled decay rates, not just
    /// independent CFL minima. This hook does not itself prove stability.
    ///
    /// The proposal excludes distance to the next event/end time. The kernel
    /// clips the result to that distance AFTER this hook, so repeated tightening
    /// cannot asymptotically approach an event without ever reaching it.
    fn limit_timestep_s(&mut self, _thermal: &DepositionSimulation<'_>, proposed_s: f64,
        _checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<f64, Self::Error> { Ok(proposed_s) }
    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>)
        -> Result<(), Self::Error>;
    /// Old temperatures and new enthalpy-backed temperatures on the same step.
    /// `duration_s` is the exact duration used for thermal integration, not a
    /// separately chosen resampling interval. No event occurs inside this step.
    fn step(&mut self, start_s: f64, duration_s: f64,
        old_temperatures_k: &BTreeMap<CellId, f64>, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<(), Self::Error>;
}

impl DepositionObserver for () {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(),Self::Error> { Ok(()) }
    fn step(&mut self, _: f64, _: f64, _: &BTreeMap<CellId,f64>, _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<(),Self::Error> { Ok(()) }
}

impl<'a> DepositionSimulation<'a> {
    // Validate the complete generated field before changing any boundary. The
    // ordinary no-source path does not even add +0, preserving its arithmetic.
    pub(super) fn add_step_sources(&mut self, sources: &BTreeMap<CellId, f64>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>)
        -> Result<Vec<(CellId, f64)>, DepositionError> {
        if sources.len() > self.cells.len() {
            return Err(DepositionError::InvalidInput("too many generated heat sources"));
        }
        let mut changes = Vec::with_capacity(sources.len());
        for (index, (&id, &watts)) in sources.iter().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            let cell = self.cells.get(&id).ok_or(DepositionError::InvalidInput("inactive generated-source cell"))?;
            finite(watts)?;
            if watts != 0.0 {
                changes.push((id, cell.boundary.source_w, finite(cell.boundary.source_w + watts)?));
            }
        }
        checkpoint()?;
        let mut authored = Vec::with_capacity(changes.len());
        for (id, old, combined) in changes {
            self.cells.get_mut(&id).expect("validated active source").boundary.source_w = combined;
            authored.push((id, old));
        }
        Ok(authored)
    }

    pub(super) fn set_step_conductances(&mut self, values: &BTreeMap<u64, f64>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>)
        -> Result<Vec<(u64, f64)>, DepositionError> {
        if values.len() > self.links.len() {
            return Err(DepositionError::InvalidInput("too many step-local conductances"));
        }
        let mut authored = Vec::with_capacity(values.len());
        for (index, (&id, &g)) in values.iter().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            let link = self.links.get(&id)
                .ok_or(DepositionError::InvalidInput("inactive step-local thermal link"))?;
            if !nonnegative(g) {
                return Err(DepositionError::InvalidInput("step-local conductance must be finite and nonnegative"));
            }
            authored.push((id, link.conductance_w_k));
        }
        checkpoint()?;
        for (&id, &g) in values {
            self.links.get_mut(&id).expect("validated active link").conductance_w_k = g;
        }
        Ok(authored)
    }

    fn add_step_transfers(&self, rates: &mut HeatRates, transfers: &[HeatTransfer],
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>)
        -> Result<(), DepositionError> {
        if transfers.len() > self.budget.max_links {
            return Err(DepositionError::Budget("internal heat transfers"));
        }
        for (index, transfer) in transfers.iter().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            if transfer.a == transfer.b || !self.cells.contains_key(&transfer.a)
                || !self.cells.contains_key(&transfer.b) {
                return Err(DepositionError::InvalidInput("invalid internal heat-transfer endpoints"));
            }
            finite(transfer.power_w)?;
        }
        for (index, transfer) in transfers.iter().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            if transfer.power_w == 0.0 { continue; }
            let a = rates.watts.get_mut(&transfer.a).expect("validated transfer endpoint");
            *a = finite(*a - transfer.power_w)?;
            let b = rates.watts.get_mut(&transfer.b).expect("validated transfer endpoint");
            *b = finite(*b + transfer.power_w)?;
        }
        // Deliberately do not modify rates.external_w.
        checkpoint()
    }

    /// Evaluate the actual instantaneous heat balance without advancing time.
    /// Uses the same observer source/conductance/transfer hooks and heat-rate
    /// owner as advance_observed. All mutable preparation is confined to clones.
    /// This is a residual for steady solving, not a finite-difference estimate
    /// from two advances or a second implementation of the transport equations.
    pub fn power_balance_observed<O: DepositionObserver>(
        &self, observer: &O,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<PowerBalance, O::Error> {
        checkpoint()?;
        let mut candidate = self.clone();
        let mut observer = observer.clone();
        observer.begin()?;
        let sources = observer.sources(&candidate, &mut checkpoint)?;
        let conductances = observer.conductances(&candidate, &mut checkpoint)?;
        let transfers = observer.heat_transfers(&candidate, &mut checkpoint)?;
        candidate.set_step_conductances(&conductances, &mut checkpoint)?;
        candidate.add_step_sources(&sources, &mut checkpoint)?;
        let mut rates = candidate.heat_rates(f64::MAX, &mut checkpoint)?;
        candidate.add_step_transfers(&mut rates, &transfers, &mut checkpoint)?;
        checkpoint()?;
        Ok(PowerBalance { temperatures_k: rates.temperatures, net_power_w: rates.watts,
            external_power_w: rates.external_w })
    }

    /// Advance thermal state and a deeply cloned observer as one candidate.
    /// Any thermal, observer, conservation, budget or cancellation error drops
    /// both. The ordinary `advance_to` delegates here with a no-op observer;
    /// default hooks preserve thermal stepping, event ordering and arithmetic.
    pub fn advance_observed<O: DepositionObserver>(
        &self, observer: &O, target_s: f64, config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, O, EnergyAudit), O::Error> {
        self.advance_observed_limited(observer, target_s, config, self.budget.max_substeps, checkpoint)
    }

    /// The same thermal/observer operation with an additional per-invocation
    /// substep allowance. Event-splitting owners pass their remaining shared
    /// allowance here instead of renewing the whole thermal budget at each cut.
    /// This can only tighten the stored ceiling, never alter it. Zero still
    /// admits same-time events or movement through an empty interval, but not
    /// a thermal step. Work refusal occurs BEFORE the next thermal evaluation.
    pub fn advance_observed_limited<O: DepositionObserver>(
        &self, observer: &O, target_s: f64, config: StepConfig, max_substeps: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, O, EnergyAudit), O::Error> {
        checkpoint()?;
        if !target_s.is_finite() || target_s < self.time_s || !positive(config.max_dt_s)
            || !nonnegative(config.energy_tolerance_j) {
            return Err(DepositionError::InvalidInput("invalid target, timestep or joule tolerance").into());
        }
        let initial_energy = finite(self.energy_j())?;
        let mut next = self.clone();
        let mut observer = observer.clone();
        observer.begin()?;
        let mut audit = EnergyAudit { incoming_mass_kg: 0.0, incoming_energy_j: 0.0,
            external_heat_j: 0.0, stored_energy_change_j: 0.0, residual_j: 0.0,
            substeps: 0, events: 0 };
        loop {
            checkpoint()?;
            while next.cursor < next.schedule.len()
                && next.schedule[next.cursor].time_s == next.time_s {
                checkpoint()?;
                let event = next.schedule[next.cursor].event.clone();
                next.apply(event.clone(), &mut audit)?;
                next.cursor += 1;
                audit.events += 1;
                observer.event(&event, &next)?;
                checkpoint()?;
            }
            if next.time_s == target_s { break; }
            let stop = next.schedule.get(next.cursor).map_or(target_s, |e| e.time_s.min(target_s));
            if next.cells.is_empty() { next.time_s = stop; continue; }
            if audit.substeps >= next.budget.max_substeps.min(max_substeps) {
                return Err(DepositionError::Budget("thermal substeps").into());
            }
            let start_s = next.time_s;
            let sources = observer.sources(&next, &mut checkpoint)?;
            checkpoint()?;
            let conductances = observer.conductances(&next, &mut checkpoint)?;
            checkpoint()?;
            let transfers = observer.heat_transfers(&next, &mut checkpoint)?;
            checkpoint()?;
            let original_links = next.set_step_conductances(&conductances, &mut checkpoint)?;
            let authored = next.add_step_sources(&sources, &mut checkpoint)?;
            let step: Result<_, O::Error> = (|| {
                let mut rates = next.heat_rates(config.max_dt_s, &mut checkpoint)?;
                next.add_step_transfers(&mut rates, &transfers, &mut checkpoint)?;
                let dt = observer.limit_timestep_s(&next, rates.explicit_dt_s, &mut checkpoint)?;
                if !positive(dt) || dt > rates.explicit_dt_s {
                    return Err(DepositionError::InvalidInput("invalid observer timestep bound").into());
                }
                rates.explicit_dt_s = dt.min(stop-start_s);
                next.thermal_step(rates, &mut audit, &mut checkpoint).map_err(Into::into)
            })();
            // Restore even on thermal refusal. Observations/events see authored
            // boundaries, never a generated source masquerading as a new input.
            for (id, watts) in authored {
                next.cells.get_mut(&id).expect("validated active source").boundary.source_w = watts;
            }
            for (id, g) in original_links {
                next.links.get_mut(&id).expect("validated active link").conductance_w_k = g;
            }
            let (dt, old_temperatures) = step?;
            let time = finite(start_s+dt)?;
            if time <= start_s { return Err(DepositionError::Arithmetic.into()); }
            next.time_s = if dt == stop-start_s { stop } else { time };
            audit.substeps += 1;
            observer.step(start_s, dt, &old_temperatures, &next, &mut checkpoint)?;
            checkpoint()?;
        }
        audit.stored_energy_change_j = finite(next.energy_j()-initial_energy)?;
        finite(next.mass_kg())?;
        audit.residual_j = finite(audit.stored_energy_change_j
            - audit.incoming_energy_j - audit.external_heat_j)?;
        if audit.residual_j.abs() > config.energy_tolerance_j {
            return Err(DepositionError::EnergyBalance(audit.residual_j).into());
        }
        checkpoint()?;
        Ok((next, observer, audit))
    }
}

/// Checked replacement of one current thermal contact at the accepted time.
/// `None` means absent/open. The expected value prevents a stale mechanical
/// proposal from overwriting a changed area, conductance or endpoint binding.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactUpdate {
    pub id: u64,
    pub expected: Option<ThermalLink>,
    pub replacement: Option<ThermalLink>,
}

impl DepositionSimulation<'_> {
    /// Apply a complete same-time contact change without advancing heat or
    /// consuming scheduled events. Process events due now with `advance_to`
    /// first. All expected contacts are checked against the original state;
    /// duplicate IDs refuse. Only the returned candidate changes interfaces.
    ///
    /// Removing all changed contacts before adding replacements allows a
    /// conservative area repartition without transient double coverage. Existing
    /// contact/cell/area admission is reused. Mass, enthalpy and time are untouched;
    /// no mechanical work is silently injected as heat. Cost is bounded by the
    /// event and active-link budgets, with existing O(links) area checks per add.
    /// This symmetric API refuses an asymmetric current patch, including a
    /// no-op, rather than losing its side areas. Use `replace_contact_patches`
    /// with the complete expected patch for such an interface.
    pub fn replace_contacts(
        &self,
        updates: &[ContactUpdate],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if updates.len() > self.budget.max_events {
            return Err(DepositionError::Budget("contact updates"));
        }
        let mut patches = Vec::with_capacity(updates.len());
        for update in updates {
            checkpoint()?;
            patches.push(ContactPatchUpdate {
                id: update.id,
                expected: update.expected.map(ContactPatch::symmetric),
                replacement: update.replacement.map(ContactPatch::symmetric),
            });
        }
        self.replace_contact_patches(&patches, checkpoint)
    }
}
