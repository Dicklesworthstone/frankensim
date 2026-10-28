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
/// After all timestep restrictions and event clipping, `step_sources` may
/// prepare duration-dependent energy exchange on that same candidate.
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
    /// Called after the start-state rate hooks. This sees step-local source/conductance
    /// overrides but unchanged enthalpy. Compose coupled decay rates, not just
    /// independent CFL minima. This hook does not itself prove stability.
    ///
    /// The proposal excludes distance to the next event/end time. The kernel
    /// clips the result to that distance AFTER this hook, so repeated tightening
    /// cannot asymptotically approach an event without ever reaching it.
    fn limit_timestep_s(&mut self, _thermal: &DepositionSimulation<'_>, proposed_s: f64,
        _checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<f64, Self::Error> { Ok(proposed_s) }
    /// Whether this observer needs an accepted interval to report all its heat.
    /// Interval-dependent owners must return true, and wrappers must propagate
    /// it. An instantaneous power balance refuses such an owner instead of
    /// silently omitting its duration-dependent contribution.
    fn requires_step_interval(&self) -> bool { false }
    /// Additional signed average heat rates [W] for the FINAL accepted interval.
    /// Called exactly once after all CFL restrictions and event/end clipping,
    /// before changing enthalpy. `start_s` is the current thermal time; `end_s`
    /// is the exact endpoint that `step` will see (including event-time rounding).
    /// `dt_s` is the duration used for every energy integral in this step.
    ///
    /// Use this for an owned transient solve whose physical loss depends on dt.
    /// Retain its tentative history in the cloned observer and accept it in
    /// `step`. Refusal or cancellation discards that history and the thermal
    /// candidate together. The unchanged thermal state still exposes start-state
    /// source/conductance overrides, but these rates are added directly to the
    /// already prepared heat balance and never persisted as boundary data.
    /// All returned heat is counted once in `external_heat_j`.
    ///
    /// This hook cannot alter the chosen timestep, topology, material law or
    /// conductances. Any needed stability bound belongs in `limit_timestep_s`;
    /// a new nonlinear heat source still requires time refinement. Owners that
    /// need this interval must also override `requires_step_interval`. Generic
    /// observer wrappers must forward both hooks.
    fn step_sources(&mut self, _thermal_start: &DepositionSimulation<'_>,
        _start_s: f64, _dt_s: f64, _end_s: f64,
        _checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error> { Ok(BTreeMap::new()) }
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

    // Duration-dependent losses enter the SAME rate field and external-energy
    // audit after dt is fixed. Never evaluate boundaries/conduction a second time.
    fn add_interval_sources(&self, rates: &mut HeatRates, sources: &BTreeMap<CellId, f64>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>)
        -> Result<(), DepositionError> {
        if sources.is_empty() { return Ok(()); }
        if sources.len() > self.cells.len() {
            return Err(DepositionError::InvalidInput("too many interval heat sources"));
        }
        for (index, (&id, &watts)) in sources.iter().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            if !self.cells.contains_key(&id) {
                return Err(DepositionError::InvalidInput("inactive interval-source cell"));
            }
            finite(watts)?;
        }
        for (index, (&id, &watts)) in sources.iter().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            if watts == 0.0 { continue; }
            let total = rates.watts.get_mut(&id).expect("validated interval-source cell");
            *total = finite(*total + watts)?;
            rates.external_w = finite(rates.external_w + watts)?;
        }
        checkpoint()
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
    /// An observer requiring accepted-interval heat refuses: there is no duration
    /// here, so returning its other rates would silently report an incomplete balance.
    pub fn power_balance_observed<O: DepositionObserver>(
        &self, observer: &O,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<PowerBalance, O::Error> {
        checkpoint()?;
        let mut candidate = self.clone();
        let mut observer = observer.clone();
        observer.begin()?;
        if observer.requires_step_interval() {
            return Err(DepositionError::InvalidInput(
                "instantaneous power balance cannot evaluate interval-dependent heat",
            ).into());
        }
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
        finite(self.energy_j())?;
        let mut next = self.clone();
        let mut observer = observer.clone();
        observer.begin()?;
        let mut unrepresented_heat_j = 0.0;
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
                let time = finite(start_s + rates.explicit_dt_s)?;
                if time <= start_s { return Err(DepositionError::Arithmetic.into()); }
                let end_s = if rates.explicit_dt_s == stop-start_s { stop } else { time };
                let interval_sources = observer.step_sources(
                    &next, start_s, rates.explicit_dt_s, end_s, &mut checkpoint,
                )?;
                next.add_interval_sources(&mut rates, &interval_sources, &mut checkpoint)?;
                let (dt, old_temperatures) = next.thermal_step(
                    rates, &mut audit, &mut unrepresented_heat_j, &mut checkpoint,
                )?;
                if unrepresented_heat_j > config.energy_tolerance_j {
                    return Err(DepositionError::Arithmetic.into());
                }
                Ok((dt, end_s, old_temperatures))
            })();
            // Restore even on thermal refusal. Observations/events see authored
            // boundaries, never a generated source masquerading as a new input.
            for (id, watts) in authored {
                next.cells.get_mut(&id).expect("validated active source").boundary.source_w = watts;
            }
            for (id, g) in original_links {
                next.links.get_mut(&id).expect("validated active link").conductance_w_k = g;
            }
            let (dt, end_s, old_temperatures) = step?;
            next.time_s = end_s;
            audit.substeps += 1;
            observer.step(start_s, dt, &old_temperatures, &next, &mut checkpoint)?;
            checkpoint()?;
        }
        // Subtract each cell's original enthalpy before reducing. Two large
        // absolute totals can lose a representable heat increment; subtracting
        // incoming birth energy afterwards repeats that cancellation.
        let mut birth_enthalpies = BTreeMap::new();
        let mut incoming = EnergySum::default();
        for (index, event) in next.schedule[self.cursor..next.cursor].iter().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            if let Event::Birth(cell) = &event.event {
                birth_enthalpies.insert(cell.id, cell.specific_h_j_kg);
                incoming.add(finite(cell.mass_kg * cell.specific_h_j_kg)?)?;
            }
        }
        let mut thermal_change = EnergySum::default();
        let mut stored_change = EnergySum::default();
        let mut originals = self.cells.iter().peekable();
        for (index, (&id, cell)) in next.cells.iter().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            let original = if originals.peek().is_some_and(|(old_id, _)| **old_id == id) {
                originals.next().map(|(_, cell)| cell)
            } else {
                None
            };
            let baseline_h = match original {
                Some(old) => old.specific_h_j_kg,
                None => *birth_enthalpies.get(&id).ok_or(DepositionError::Arithmetic)?,
            };
            let change = finite(cell.mass_kg * finite(cell.specific_h_j_kg - baseline_h)?)?;
            thermal_change.add(change)?;
            stored_change.add(if original.is_some() {
                change
            } else {
                finite(cell.mass_kg * cell.specific_h_j_kg)?
            })?;
        }
        if originals.next().is_some() { return Err(DepositionError::Arithmetic.into()); }
        finite(next.energy_j())?;
        audit.incoming_energy_j = incoming.value()?;
        audit.stored_energy_change_j = stored_change.value()?;
        finite(next.mass_kg())?;
        audit.residual_j = finite(thermal_change.value()? - audit.external_heat_j)?;
        if audit.residual_j.abs() > config.energy_tolerance_j {
            return Err(DepositionError::EnergyBalance(audit.residual_j).into());
        }
        checkpoint()?;
        Ok((next, observer, audit))
    }
}

// Fixed cell order and Neumaier correction retain small terms of either sign.
// The sum contains physical changes, not arbitrary formation-energy datums.
#[derive(Default)]
struct EnergySum {
    sum: f64,
    correction: f64,
}
impl EnergySum {
    fn add(&mut self, value: f64) -> Result<(), DepositionError> {
        finite(value)?;
        let next = finite(self.sum + value)?;
        let lost = if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.correction = finite(self.correction + lost)?;
        self.sum = next;
        Ok(())
    }
    fn value(self) -> Result<f64, DepositionError> {
        finite(self.sum + self.correction)
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
    /// Apply signed endpoint heat [J] through the original enthalpy owner.
    ///
    /// This operation does not advance time, re-integrate sources or contacts,
    /// consume events, or change mass and geometry. Process due events first.
    /// The original caloric law resolves temperature, including latent storage.
    /// Both the local error and the sum of absolute cell-storage errors must
    /// fit the explicit joule allowances; opposite lost increments cannot cancel.
    /// The returned heat is external to this thermal subsystem. Its coupling
    /// owner must account for the corresponding upstream physical energy debit.
    /// All changes remain on a candidate until validation and cancellation pass.
    pub fn apply_heat_j(
        &self,
        heat_by_cell_j: &BTreeMap<CellId, f64>,
        cell_energy_tolerance_j: f64,
        energy_tolerance_j: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, EnergyAudit), DepositionError> {
        checkpoint()?;
        if !nonnegative(cell_energy_tolerance_j) || !nonnegative(energy_tolerance_j) {
            return Err(DepositionError::InvalidInput("finite nonnegative heat-impulse tolerances required"));
        }
        if heat_by_cell_j.len() > self.cells.len() {
            return Err(DepositionError::InvalidInput("too many heat-impulse cells"));
        }
        if self.schedule.get(self.cursor).is_some_and(|e| e.time_s <= self.time_s) {
            return Err(DepositionError::InvalidInput("process due scheduled events before heat impulses"));
        }
        let mut next = self.clone();
        let mut supplied = EnergySum::default();
        let mut stored = EnergySum::default();
        let mut residual = EnergySum::default();
        let mut absolute_residual = EnergySum::default();
        for (index, (&id, &heat)) in heat_by_cell_j.iter().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            finite(heat)?;
            let cell = next.cells.get_mut(&id)
                .ok_or(DepositionError::InvalidInput("inactive heat-impulse cell"))?;
            if heat == 0.0 { continue; }
            let old_h = cell.specific_h_j_kg;
            let new_h = finite(old_h + finite(heat / cell.mass_kg)?)?;
            let temperature = next.materials[&cell.material_id].temperature_k(new_h)?;
            if !positive(temperature) { return Err(DepositionError::MaterialDomain); }
            finite(cell.mass_kg * new_h)?;
            let change = finite(cell.mass_kg * finite(new_h - old_h)?)?;
            let defect = finite(change - heat)?;
            if defect.abs() > cell_energy_tolerance_j {
                return Err(DepositionError::EnergyBalance(defect));
            }
            supplied.add(heat)?;
            stored.add(change)?;
            residual.add(defect)?;
            absolute_residual.add(defect.abs())?;
            cell.specific_h_j_kg = new_h;
        }
        let absolute_residual = absolute_residual.value()?;
        if absolute_residual > energy_tolerance_j {
            return Err(DepositionError::EnergyBalance(absolute_residual));
        }
        finite(next.energy_j())?;
        let audit = EnergyAudit {
            incoming_mass_kg: 0.0,
            incoming_energy_j: 0.0,
            external_heat_j: supplied.value()?,
            stored_energy_change_j: stored.value()?,
            residual_j: residual.value()?,
            substeps: 0,
            events: 0,
        };
        checkpoint()?;
        Ok((next, audit))
    }

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

#[cfg(test)]
#[path = "observe/energy_tests.rs"]
mod energy_tests;
