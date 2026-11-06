//! Entropy transported by the original fixed-mass upwind graph.
//!
//! Every flowing component uses one immutable caloric law and convention.
//! Only NET boundary entropy is reported: a common relative entropy reference
//! cancels because the existing owner balances mass. Chemical mixing, unequal
//! pressure conventions and unlike-law transport require another entropy owner.
//! With a stationary solid skeleton, only the moving liquid's enthalpy and
//! entropy cross faces; original mixture entropy remains in cell storage.

use super::liquid::LiquidEnthalpyCarrier;
use super::{Advected, MassBoundary, MassFlow, finite, invalid};
use crate::enthalpy::{CellId, DepositionError, DepositionObserver, DepositionSimulation,
    EnergyAudit, MaterialId, StepConfig};
use crate::enthalpy::backward_euler::TransportStepCoefficients;
use crate::enthalpy::entropy::{AdvectiveEntropyContribution, EntropyAudit, EntropyObserver};
use std::collections::BTreeMap;

/// Entropy balance over the last successful advance, in J/K.
///
/// Storage equals boundary input + upwind mixing + temporal defect + other
/// energy/T + advective balance residue, up to floating-point reduction.
/// This is not a physical entropy-production certificate for a coupled flow.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AdvectionEntropyAudit {
    /// Exact sum of m*integral(dh/T), starting from each actual substep state.
    /// Excludes entropy carried by scheduled Birth events between substeps.
    pub thermal_storage_change_j_k: f64,
    /// Net inlet donor entropy minus CURRENT outlet donor entropy, integrated
    /// on the same starting-state quadrature as the original enthalpy flux.
    /// This is independent of a common absolute entropy datum for balanced mass.
    /// Liquid-only carriage uses the pure-liquid side at coexistence, even
    /// when a cell retains stationary solid and its bulk entropy is smaller.
    pub advective_boundary_entropy_input_j_k: f64,
    /// Sum dm*[(h_donor-h_receiver)/T_receiver - (s_donor-s_receiver)].
    /// This is spatial upwind/control-volume mixing, not constitutive fluid
    /// dissipation. It is nonnegative for monotone T(h) in exact arithmetic;
    /// small signed roundoff remains visible and is never clamped into heat.
    /// Donor and receiver h/s are those of the actual carrier. An isothermal
    /// liquid loop has zero mixing regardless of stationary solid fractions.
    pub upwind_mixing_j_k: f64,
    /// Sum m*[integral(dh/T) - delta_h/T_start] for the actual accepted update.
    /// Separates explicit temporal entropy error from spatial upwind mixing.
    pub temporal_defect_j_k: f64,
    /// Actual stored energy minus this observer's advective energy, divided
    /// by each starting cell temperature. Includes other heat/work and their
    /// roundoff; it is NOT an assigned external source entropy or production.
    pub other_energy_at_cell_temperature_j_k: f64,
    /// Advective energy/T minus boundary entropy and upwind mixing. Keeps
    /// mass-balance tolerance, flux arithmetic and reduction error separate.
    pub advective_balance_residual_j_k: f64,
}

/// Instantaneous entropy terms from the original upwind carrier endpoint, W/K.
/// No physical interval, accepted history or extra heat source is implied.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AdvectionEntropyRates {
    pub boundary_entropy_input_w_k: f64,
    pub upwind_mixing_w_k: f64,
    /// Original transport power/T minus boundary entropy and upwind mixing.
    pub balance_residual_w_k: f64,
}

pub(in crate::enthalpy) struct EndpointAdvectionEntropy {
    pub(in crate::enthalpy) rates: AdvectionEntropyRates,
    pub(in crate::enthalpy) temperatures: BTreeMap<CellId, f64>,
}

impl AdvectionEntropyRates {
    pub(in crate::enthalpy) fn integrated(
        self,
        duration_s: f64,
    ) -> Result<AdvectiveEntropyContribution, DepositionError> {
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(invalid("advective entropy requires positive quadrature"));
        }
        Ok(AdvectiveEntropyContribution {
            boundary_entropy_input_j_k: finite(duration_s * self.boundary_entropy_input_w_k)?,
            upwind_mixing_j_k: finite(duration_s * self.upwind_mixing_w_k)?,
            balance_residual_j_k: finite(duration_s * self.balance_residual_w_k)?,
        })
    }
}

impl Advected<()> {
    /// Query native upwind carried enthalpy and entropy at the SAME endpoint.
    ///
    /// The original transport owner admits fixed mass, continuity, caloric
    /// storage and source conventions before entropy is evaluated. Inlet/outlet
    /// entropy and control-volume mixing use the same immutable law as their
    /// donor/receiver enthalpy, with no arbitrary absolute entropy datum.
    /// Temperature and caloric entropy must be supported by the original laws.
    ///
    /// Ordinary bulk flow carries the stored h. An admitted stationary-skeleton
    /// liquid carrier instead transports its original pure-liquid h and entropy,
    /// including on the latent plateau. Mixture h and entropy remain in storage;
    /// local liquid mass balance cancels the stationary-phase offset. No second
    /// latent heat or entropy source is introduced.
    ///
    /// This stateless query does not advance an observer or impose an explicit
    /// flow CFL. Stationary species still require a chemical-mixing entropy
    /// closure even when their energy transport is implicit. Selective outlets
    /// and attached flow histories also refuse. Reconstructed faces require their
    /// own entropy closure and are not substituted here.
    pub fn endpoint_transport_with_entropy(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(TransportStepCoefficients<()>, AdvectionEntropyRates), DepositionError> {
        self.endpoint_transport_entropy_state(thermal, checkpoint)
            .map(|(transport, entropy)| (transport, entropy.rates))
    }

    pub(in crate::enthalpy) fn endpoint_transport_entropy_state(
        &self,
        thermal: &DepositionSimulation<'_>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(TransportStepCoefficients<()>, EndpointAdvectionEntropy), DepositionError> {
        checkpoint()?;
        self.admit_implicit_transport()?;
        if self.stationary_species.is_some() {
            return Err(invalid("thermal advection entropy does not cover stationary species mixing"));
        }
        let transport = self.endpoint_transport(thermal, &mut checkpoint)?;
        let mut enthalpies = BTreeMap::new();
        let mut temperatures = BTreeMap::new();
        for (i, (&id, cell)) in thermal.cells().iter().enumerate() {
            if i % 256 == 0 { checkpoint()?; }
            enthalpies.insert(id, cell.specific_h_j_kg);
            temperatures.insert(id, thermal.temperature_k(id)?);
        }
        let (boundary, mixing) = flux_entropy_rates(
            thermal, &self.flows, &self.boundaries, &self.reference_cells,
            self.liquid_carrier, &enthalpies, &temperatures, None, &mut checkpoint,
        )?;
        // Use the exact powers returned by the physical transport owner, not
        // a reconstructed stored-energy difference or an independently evolved
        // flow. Continuity/flux roundoff stays in its own balance residue.
        let mut at_cell = Sum::default();
        for (&id, &watts) in &transport.heat.bulk.sources_w {
            checkpoint()?;
            at_cell.add(finite(watts / temperature(&temperatures, id)?)?)?;
        }
        for transfer in &transport.transfers {
            checkpoint()?;
            at_cell.add(finite(-transfer.power_w / temperature(&temperatures, transfer.a)?)?)?;
            at_cell.add(finite(transfer.power_w / temperature(&temperatures, transfer.b)?)?)?;
        }
        let rates = AdvectionEntropyRates {
            boundary_entropy_input_w_k: boundary,
            upwind_mixing_w_k: mixing,
            balance_residual_w_k: finite(at_cell.total()? - boundary - mixing)?,
        };
        checkpoint()?;
        Ok((transport, EndpointAdvectionEntropy { rates, temperatures }))
    }
}

impl<O> Advected<O> {
    // Preserve accepted history and every transport binding while composing
    // the two existing audit owners. No topology or physical owner is rebuilt.
    fn map_inner<N>(self, map: impl FnOnce(O) -> N) -> Advected<N> {
        Advected {
            inner: map(self.inner), flows: self.flows, boundaries: self.boundaries,
            reference_cells: self.reference_cells, selective_outlets: self.selective_outlets,
            stationary_species: self.stationary_species, liquid_carrier: self.liquid_carrier,
            boundary_pending: self.boundary_pending, pending: self.pending,
            audit: self.audit, entropy: self.entropy,
        }
    }

    /// Audit entropy on the SAME accepted upwind trajectory. Every active law
    /// must provide its exact dh/T integral. Flowing cells and inlets must use
    /// the same immutable MaterialId; a shared enthalpy datum alone is not
    /// enough. Stationary composition mixing deliberately refuses this audit.
    ///
    /// For an admitted liquid carrier, fluxes use `h_liquid=max(h,h_liquidus)`
    /// and the original law's entropy at that liquid state. Cell storage still
    /// uses the actual bulk h. At pure-substance coexistence,
    /// `s_liquid-s_bulk=(h_liquid-h_bulk)/T_m`; local liquid mass balance
    /// cancels this stationary-phase offset. No additional melting entropy or
    /// latent energy source is introduced. This does not cover alloy segregation,
    /// phase kinetics, unequal phase densities or a moving solid skeleton.
    ///
    /// Extra work is sparse and only occurs when enabled. The original energy,
    /// mass, CFL, event and transactional refusal semantics remain unchanged.
    /// Nonadvective exchange remains explicitly unclassified until its owner
    /// supplies actual reservoir temperatures and heat/work entropy semantics.
    #[must_use]
    pub fn with_entropy_audit(mut self) -> Self {
        self.entropy = Some(EntropyState::default());
        self
    }

    /// Accepted entropy account; None when not requested. Reset by begin() on
    /// the tentative observer, together with the ordinary advection account.
    #[must_use]
    pub fn entropy_audit(&self) -> Option<AdvectionEntropyAudit> {
        self.entropy.as_ref().map(|state| state.audit)
    }
}

impl DepositionSimulation<'_> {
    /// Advance native fixed-mass flow with one classified thermal entropy audit.
    ///
    /// The original `Advected` owner supplies net stream entropy, upwind mixing
    /// and its balance residue on the SAME starting-state quadrature as its
    /// enthalpy flux. Conduction retains actual bath temperatures and applied
    /// coefficients; other sources/transfers remain explicitly unclassified.
    /// No advective energy/T is counted again as an unclassified heat source.
    /// The signed temporal defect includes caloric/time/roundoff error, not
    /// constitutive fluid dissipation or chemical entropy production.
    ///
    /// Every energy update, CFL restriction, event and original observer hook
    /// remains on the existing trajectory. The returned flow retains its inner
    /// history and has `entropy_audit()` enabled; per-call accounts reset through
    /// the usual `begin()` hook. The input thermal and observer states are
    /// unchanged on failure, including missing entropy or late cancellation.
    /// Same-law balanced bulk flow and the existing liquid carrier are supported;
    /// stationary species and selective outlets still require another closure.
    pub fn advance_advected_with_entropy<O: DepositionObserver>(
        &self, observer: &Advected<O>, target_s: f64, config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, Advected<O>, EnergyAudit, EntropyAudit), O::Error> {
        // Put thermal observation INSIDE advection: it retains only genuine
        // nonadvective sources, including their flags when signed terms cancel.
        // Subtracting accumulated transport from a combined Q/T would lose that
        // distinction and magnify cancellation at disparate source scales.
        let wrapped = observer.clone().with_entropy_audit().map_inner(EntropyObserver::new);
        let (next, wrapped, energy) = self.advance_observed(&wrapped, target_s, config, checkpoint)?;
        let flow = wrapped.entropy_audit().ok_or(DepositionError::Arithmetic)?;
        let mut entropy = wrapped.inner.audit;
        entropy.advection = Some(AdvectiveEntropyContribution {
            boundary_entropy_input_j_k: flow.advective_boundary_entropy_input_j_k,
            upwind_mixing_j_k: flow.upwind_mixing_j_k,
            balance_residual_j_k: flow.advective_balance_residual_j_k,
        });
        entropy.refresh_temporal_defect()?;
        Ok((next, wrapped.map_inner(|observer| observer.inner), energy, entropy))
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), DepositionError> {
        let term = finite(term)?;
        let next = finite(self.value + term)?;
        let lost = if self.value.abs() >= term.abs() {
            (self.value - next) + term
        } else {
            (term - next) + self.value
        };
        self.correction = finite(self.correction + lost)?;
        self.value = next;
        Ok(())
    }
    fn total(&self) -> Result<f64, DepositionError> {
        finite(self.value + self.correction)
    }
}

#[derive(Clone, Default)]
pub(super) struct EntropyState {
    initial_h: BTreeMap<CellId, f64>,
    sums: [Sum; 6],
    audit: AdvectionEntropyAudit,
}

fn difference(
    thermal: &DepositionSimulation<'_>,
    material: MaterialId,
    from: f64,
    to: f64,
) -> Result<f64, DepositionError> {
    let law = thermal
        .materials
        .get(&material)
        .ok_or(invalid("unregistered advective entropy law"))?;
    finite(
        law.specific_entropy_change_j_kg_k(from, to)?
            .ok_or(invalid(
                "thermal law does not provide advective entropy change",
            ))?,
    )
}

fn temperature(temperatures: &BTreeMap<CellId, f64>, id: CellId) -> Result<f64, DepositionError> {
    let t = *temperatures.get(&id).ok_or(DepositionError::Arithmetic)?;
    if t.is_finite() && t > 0.0 {
        Ok(t)
    } else {
        Err(DepositionError::MaterialDomain)
    }
}

fn add_power(
    powers: &mut BTreeMap<CellId, f64>,
    id: CellId,
    watts: f64,
) -> Result<(), DepositionError> {
    let power = powers.get_mut(&id).ok_or(DepositionError::Arithmetic)?;
    *power = finite(*power + watts)?;
    Ok(())
}

// One original entropy-flux law for explicit start states and implicit
// endpoint stages. Only explicit history needs the per-cell power reconstruction;
// implicit callers retain the native transport owner's actual returned powers.
#[allow(clippy::too_many_arguments)]
fn flux_entropy_rates(
    thermal: &DepositionSimulation<'_>,
    flows: &[MassFlow],
    boundaries: &[MassBoundary],
    energy_reference_cells: &BTreeMap<CellId, CellId>,
    liquid_carrier: Option<LiquidEnthalpyCarrier>,
    enthalpies: &BTreeMap<CellId, f64>,
    temperatures: &BTreeMap<CellId, f64>,
    mut powers: Option<&mut BTreeMap<CellId, f64>>,
    check: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(f64, f64), DepositionError> {
    // At coexistence the carrier uses original pure-liquid h/s; stationary
    // solid remains in cell storage. Explicit and implicit callers supply the
    // same original descriptor; ordinary bulk flow passes no carrier.
    let carried = |id| -> Result<f64, DepositionError> {
        let h = *enthalpies.get(&id).ok_or(DepositionError::Arithmetic)?;
        match liquid_carrier {
            Some(carrier) => carrier.carried_h_at(h),
            None => Ok(h),
        }
    };
    let mut references = BTreeMap::new();
    let mut boundary = Sum::default();
    let mut mixing = Sum::default();
    for (i, entry) in boundaries.iter().enumerate() {
        if i % 256 == 0 { check()?; }
        let (_, id, rate) = entry.identity();
        let cell = thermal.cells().get(&id).ok_or(DepositionError::Arithmetic)?;
        let h = carried(id)?;
        let reference = *references.entry(cell.material_id).or_insert(h);
        let t = temperature(temperatures, id)?;
        let (boundary_s, carried_h, sign) = match *entry {
            MassBoundary::Inlet { material_id, specific_h_j_kg: donor_h, .. } => {
                if material_id != cell.material_id {
                    return Err(invalid("advective entropy requires one immutable caloric law"));
                }
                let incoming_s = difference(thermal, material_id, reference, donor_h)?;
                let ds = difference(thermal, material_id, h, donor_h)?;
                mixing.add(finite(rate * finite((donor_h - h) / t - ds)?)?)?;
                (incoming_s, donor_h, 1.0)
            }
            MassBoundary::Outlet { .. } => {
                (difference(thermal, cell.material_id, reference, h)?, h, -1.0)
            }
        };
        boundary.add(sign * finite(rate * boundary_s)?)?;
        if let Some(powers) = powers.as_deref_mut() {
            let energy_reference = carried(*energy_reference_cells.get(&id).unwrap_or(&id))?;
            add_power(powers, id, sign * finite(rate * finite(carried_h - energy_reference)?)?)?;
        }
    }
    for (i, flow) in flows.iter().enumerate() {
        if i % 256 == 0 { check()?; }
        let Some((up, down, rate)) = flow.directed() else { continue; };
        let donor = thermal.cells().get(&up).ok_or(DepositionError::Arithmetic)?;
        let receiver = thermal.cells().get(&down).ok_or(DepositionError::Arithmetic)?;
        if donor.material_id != receiver.material_id {
            return Err(invalid("advective entropy requires one immutable caloric law"));
        }
        let hu = carried(up)?;
        let hd = carried(down)?;
        let td = temperature(temperatures, down)?;
        let ds = difference(thermal, donor.material_id, hd, hu)?;
        mixing.add(finite(rate * finite((hu - hd) / td - ds)?)?)?;
        if let Some(powers) = powers.as_deref_mut() {
            let energy_reference = carried(*energy_reference_cells.get(&up).unwrap_or(&up))?;
            let power = finite(rate * finite(hu - energy_reference)?)?;
            add_power(powers, up, -power)?;
            add_power(powers, down, power)?;
        }
    }
    Ok((boundary.total()?, mixing.total()?))
}

impl EntropyState {
    pub(super) fn reset(&mut self) {
        *self = Self::default();
    }

    pub(super) fn capture(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        self.initial_h.clear();
        for (i, (&id, cell)) in thermal.cells().iter().enumerate() {
            if i % 256 == 0 {
                check()?;
            }
            self.initial_h.insert(id, cell.specific_h_j_kg);
        }
        Ok(())
    }

    pub(super) fn record(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        flows: &[MassFlow],
        boundaries: &[MassBoundary],
        energy_reference_cells: &BTreeMap<CellId, CellId>,
        liquid_carrier: Option<LiquidEnthalpyCarrier>,
        temperatures: &BTreeMap<CellId, f64>,
        dt: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        check()?;
        if self.initial_h.len() != thermal.cells().len()
            || temperatures.len() != thermal.cells().len()
            || !dt.is_finite()
            || dt <= 0.0
        {
            return Err(invalid("advective entropy must retain its active cells"));
        }
        let initial = &self.initial_h;
        let mut powers: BTreeMap<_, _> = initial.keys().map(|&id| (id, 0.0)).collect();
        let (boundary, mixing) = flux_entropy_rates(
            thermal, flows, boundaries, energy_reference_cells, liquid_carrier,
            initial, temperatures, Some(&mut powers), check,
        )?;
        let mut storage = Sum::default();
        let mut temporal = Sum::default();
        let mut other = Sum::default();
        let mut advective = Sum::default();
        for (i, (&id, cell)) in thermal.cells().iter().enumerate() {
            if i % 256 == 0 {
                check()?;
            }
            let h0 = *initial.get(&id).ok_or(DepositionError::Arithmetic)?;
            let t0 = temperature(temperatures, id)?;
            let dh = finite(cell.specific_h_j_kg - h0)?;
            let ds = difference(thermal, cell.material_id, h0, cell.specific_h_j_kg)?;
            let energy = finite(cell.mass_kg * dh)?;
            let transported = finite(dt * powers[&id])?;
            storage.add(finite(cell.mass_kg * ds)?)?;
            temporal.add(finite(cell.mass_kg * finite(ds - dh / t0)?)?)?;
            other.add(finite((energy - transported) / t0)?)?;
            advective.add(finite(transported / t0)?)?;
        }
        let boundary = finite(dt * boundary)?;
        let mixing = finite(dt * mixing)?;
        let values = [
            storage.total()?,
            boundary,
            mixing,
            temporal.total()?,
            other.total()?,
            finite(advective.total()? - boundary - mixing)?,
        ];
        for (sum, value) in self.sums.iter_mut().zip(values) {
            sum.add(value)?;
        }
        self.audit = AdvectionEntropyAudit {
            thermal_storage_change_j_k: self.sums[0].total()?,
            advective_boundary_entropy_input_j_k: self.sums[1].total()?,
            upwind_mixing_j_k: self.sums[2].total()?,
            temporal_defect_j_k: self.sums[3].total()?,
            other_energy_at_cell_temperature_j_k: self.sums[4].total()?,
            advective_balance_residual_j_k: self.sums[5].total()?,
        };
        check()?;
        Ok(())
    }
}

#[cfg(test)]
#[path = "entropy/tests.rs"]
mod tests;
