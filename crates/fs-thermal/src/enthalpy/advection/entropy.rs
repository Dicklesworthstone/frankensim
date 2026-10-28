//! Entropy transported by the original fixed-mass upwind graph.
//!
//! Every flowing component uses one immutable caloric law and convention.
//! Only NET boundary entropy is reported: a common relative entropy reference
//! cancels because the existing owner balances mass. Chemical mixing, unequal
//! pressure conventions and unlike-law transport require another entropy owner.

use super::{Advected, MassBoundary, MassFlow, finite, invalid};
use crate::enthalpy::{CellId, DepositionError, DepositionSimulation, MaterialId};
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
    pub advective_boundary_entropy_input_j_k: f64,
    /// Sum dm*[(h_donor-h_receiver)/T_receiver - (s_donor-s_receiver)].
    /// This is spatial upwind/control-volume mixing, not constitutive fluid
    /// dissipation. It is nonnegative for monotone T(h) in exact arithmetic;
    /// small signed roundoff remains visible and is never clamped into heat.
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

impl<O> Advected<O> {
    /// Audit entropy on the SAME accepted upwind trajectory. Every active law
    /// must provide its exact dh/T integral. Flowing cells and inlets must use
    /// the same immutable MaterialId; a shared enthalpy datum alone is not
    /// enough. Stationary composition mixing deliberately refuses this audit.
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
        // A single relative reference per immutable material law. Separate
        // disconnected components with this same law can share the reference.
        let mut references = BTreeMap::new();
        let mut boundary = Sum::default();
        let mut mixing = Sum::default();
        for (i, entry) in boundaries.iter().enumerate() {
            if i % 256 == 0 {
                check()?;
            }
            let (_, id, rate) = entry.identity();
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(DepositionError::Arithmetic)?;
            let h = *initial.get(&id).ok_or(DepositionError::Arithmetic)?;
            let reference = *references.entry(cell.material_id).or_insert(h);
            let t = temperature(temperatures, id)?;
            match *entry {
                MassBoundary::Inlet {
                    material_id,
                    specific_h_j_kg: donor_h,
                    ..
                } => {
                    if material_id != cell.material_id {
                        return Err(invalid(
                            "advective entropy requires one immutable caloric law",
                        ));
                    }
                    let incoming_s = difference(thermal, material_id, reference, donor_h)?;
                    boundary.add(finite(rate * incoming_s)?)?;
                    let ds = difference(thermal, material_id, h, donor_h)?;
                    mixing.add(finite(rate * finite((donor_h - h) / t - ds)?)?)?;
                    add_power(&mut powers, id, finite(rate * donor_h)?)?;
                }
                MassBoundary::Outlet { .. } => {
                    let outgoing_s = difference(thermal, cell.material_id, reference, h)?;
                    boundary.add(-finite(rate * outgoing_s)?)?;
                    add_power(&mut powers, id, -finite(rate * h)?)?;
                }
            }
        }
        for (i, flow) in flows.iter().enumerate() {
            if i % 256 == 0 {
                check()?;
            }
            let Some((up, down, rate)) = flow.directed() else {
                continue;
            };
            let donor = thermal
                .cells()
                .get(&up)
                .ok_or(DepositionError::Arithmetic)?;
            let receiver = thermal
                .cells()
                .get(&down)
                .ok_or(DepositionError::Arithmetic)?;
            if donor.material_id != receiver.material_id {
                return Err(invalid(
                    "advective entropy requires one immutable caloric law",
                ));
            }
            let hu = *initial.get(&up).ok_or(DepositionError::Arithmetic)?;
            let hd = *initial.get(&down).ok_or(DepositionError::Arithmetic)?;
            let td = temperature(temperatures, down)?;
            let ds = difference(thermal, donor.material_id, hd, hu)?;
            mixing.add(finite(rate * finite((hu - hd) / td - ds)?)?)?;
            let power = finite(rate * hu)?;
            add_power(&mut powers, up, -power)?;
            add_power(&mut powers, down, power)?;
        }
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
        let boundary = finite(dt * boundary.total()?)?;
        let mixing = finite(dt * mixing.total()?)?;
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
