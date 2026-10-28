//! A bounded graph of conservative nozzle transfers on existing gas parcels.
//!
//! Each physical edge uses the original nozzle and species/enthalpy transfer
//! owner once. Edges are Lie-split in increasing stable-ID order: each acts for
//! the SAME interval on the result of preceding edges. This is first order,
//! not a simultaneous pressure solve. Input vector order has no effect, but
//! changing edge IDs can change finite-step splitting error; refine time.
/// Unsplit alternative: evaluate all physical links from one accepted material snapshot.
pub mod snapshot;

use super::{GasNozzle, GasNozzleFlow, GasNozzleStepSolve};
use super::super::{inventory, ClosedGasParcel, GasEnergyError, GasMassTransferStep};
use std::collections::BTreeMap;

/// One physical restriction. Direction follows its current staged pressures.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasNetworkLink {
    pub id: u64,
    /// Two distinct existing parcel IDs; either may be upstream.
    pub parcels: [u64; 2],
    pub nozzle: GasNozzle,
    /// Original throat and pair-transfer admission, including donor fraction.
    pub solve: GasNozzleStepSolve,
}
/// Global admission in addition to every original pair's physical checks.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasNetworkSolve {
    pub maximum_parcels: usize,
    pub maximum_links: usize,
    /// Total GROSS mass withdrawn from any parcel / its initial mass, in (0,1).
    /// Incoming mass cannot conceal excessive withdrawal through parallel edges.
    pub maximum_withdrawn_fraction: f64,
    /// Absolute net mass AND every global species-mass residual [kg].
    pub mass_tolerance_kg: f64,
    /// Absolute sum of original-caloric link residuals after species datums cancel [J].
    pub energy_tolerance_j: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub enum GasNetworkError {
    Invalid(&'static str),
    Link { id: u64, error: GasEnergyError },
    Parcel { id: u64, error: GasEnergyError },
    Withdrawal { id: u64, fraction: f64, maximum: f64 },
    Cancelled,
}
impl core::fmt::Display for GasNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "gas transfer network: {self:?}")
    }
}
impl std::error::Error for GasNetworkError {}
#[derive(Clone, Debug, PartialEq)]
pub struct GasNetworkLinkStep {
    pub id: u64,
    /// Actual upwind order. At zero flow this is only the evaluated orientation.
    pub donor: u64,
    pub receiver: u64,
    pub flow: GasNozzleFlow,
    /// Original source-bearing pair account; never an extra heat input.
    pub transfer: GasMassTransferStep,
}
#[derive(Clone, Debug, PartialEq)]
pub struct GasNetworkParcelStep {
    pub id: u64,
    pub mass_change_kg: f64,
    /// Sum of this parcel's accepted species-reference and sensible changes [J].
    /// Inventory rounding is separately bounded by the actual mass/species audit.
    pub energy_change_j: f64,
    pub gross_withdrawn_mass_kg: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct GasNetworkStep {
    pub duration_s: f64,
    pub links: Vec<GasNetworkLinkStep>,
    pub parcels: Vec<GasNetworkParcelStep>,
    pub total_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    /// Telescoping sum of accepted link energy residuals, after datum cancellation.
    /// Subtracting the rounded absolute parcel reports need not resolve this value.
    pub total_energy_residual_j: f64,
}
/// All endpoints, or none. No mutation of accepted parcels occurs in a trial.
#[derive(Clone, Debug, PartialEq)]
pub struct GasNetworkCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    step: GasNetworkStep,
}
impl GasNetworkCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> { &self.parcels }
    pub const fn step(&self) -> &GasNetworkStep { &self.step }
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, GasNetworkStep) {
        (self.parcels, self.step)
    }
}
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), GasNetworkError> {
    if cancelled() { Err(GasNetworkError::Cancelled) } else { Ok(()) }
}
fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }

#[derive(Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), GasNetworkError> {
        let adjusted = term - self.correction;
        let value = self.value + adjusted;
        let correction = (value - self.value) - adjusted;
        if ![term, adjusted, value, correction].iter().all(|x| x.is_finite()) {
            return Err(GasNetworkError::Invalid("gas network energy accumulation is unrepresentable"));
        }
        self.value = value; self.correction = correction;
        Ok(())
    }
}

/// Compose arbitrary, including branching/cyclic, fixed nozzle graphs.
/// Each pair's source consistency, donor limit and pressure-crossing refusal
/// remains active. Shared nodes retain their current staged composition, not
/// an averaged medium. Final mass and every species are independently checked
/// from actual endpoints. Energy telescopes the original-caloric link residuals:
/// each link has already checked its actual sensible endpoint storage against
/// the conserved species enthalpy flux. Re-subtracting large absolute formation
/// energies would lose that physical balance; inventory rounding has its own
/// mass allowance. This is a sequential balance, not a time-discretization bound.
///
/// No clock advances during edge composition; surrounding heat, reaction and
/// motion operators act over this same interval. Network storage is the sum of
/// original parcel storage; carried enthalpy and thermalized jet energy are
/// internal transfers. No extra pressure work, heat, pipe storage or jet thrust
/// is inferred. A closed/no-flow edge preserves the original pair behavior.
/// Callers supply physical initial states and fixed source/geometry bindings.
/// Refusal calls for a shorter/admitted interval, never clipped mass or energy.
pub fn step_candidate(
    parcels: &BTreeMap<u64, ClosedGasParcel>, links: &[GasNetworkLink], duration_s: f64,
    controls: GasNetworkSolve, mut cancelled: impl FnMut() -> bool,
) -> Result<GasNetworkCandidate, GasNetworkError> {
    use GasNetworkError as Error;
    check(&mut cancelled)?;
    if !positive(duration_s) || controls.maximum_parcels == 0
        || parcels.len() > controls.maximum_parcels || links.len() > controls.maximum_links
        || !positive(controls.maximum_withdrawn_fraction) || controls.maximum_withdrawn_fraction >= 1.0
        || !positive(controls.mass_tolerance_kg) || !positive(controls.energy_tolerance_j) {
        return Err(Error::Invalid("bounded gas graph, positive time/tolerances and withdrawal fraction in (0,1) required"));
    }
    let mut ordered = BTreeMap::new();
    for link in links {
        check(&mut cancelled)?;
        if link.parcels[0] == link.parcels[1] || link.parcels.iter().any(|id| !parcels.contains_key(id))
            || ordered.insert(link.id, link).is_some() {
            return Err(Error::Invalid("unique link IDs and two distinct existing parcels required"));
        }
    }
    let mut next = BTreeMap::new();
    for (&id, parcel) in parcels {
        check(&mut cancelled)?;
        next.insert(id, parcel.clone());
    }
    let mut withdrawn = BTreeMap::<u64, f64>::new();
    let mut steps = Vec::with_capacity(links.len());
    let mut storage_changes = BTreeMap::<u64, Sum>::new();
    let mut energy = Sum::default();
    for (&id, link) in &ordered {
        check(&mut cancelled)?;
        let [a, b] = link.parcels;
        let (donor, receiver) = if next[&a].state().state().pressure >= next[&b].state().state().pressure {
            (a, b)
        } else { (b, a) };
        let (candidate, flow) = next[&donor].through_nozzle_to(&next[&receiver], link.nozzle,
            duration_s, link.solve, &mut cancelled).map_err(|error| {
                if error == GasEnergyError::Cancelled { Error::Cancelled } else { Error::Link { id, error } }
            })?;
        let gross = withdrawn.get(&donor).copied().unwrap_or(0.0) + candidate.step().transferred_mass_kg;
        let fraction = gross / parcels[&donor].mass_kg();
        if !fraction.is_finite() || fraction > controls.maximum_withdrawn_fraction {
            return Err(Error::Withdrawal { id: donor, fraction, maximum: controls.maximum_withdrawn_fraction });
        }
        let (upstream, downstream, transfer) = candidate.into_parts();
        storage_changes.entry(donor).or_default().add(transfer.donor_energy_change_j)?;
        storage_changes.entry(receiver).or_default().add(transfer.receiver_energy_change_j)?;
        // The staged endpoint of one link is exactly the initial state of the
        // next link. Internal enthalpy fluxes cancel before these small residuals
        // are accumulated; no new storage or heat is introduced by the network.
        energy.add(transfer.energy_residual_j)?;
        withdrawn.insert(donor, gross);
        next.insert(donor, upstream);
        next.insert(receiver, downstream);
        steps.push(GasNetworkLinkStep { id, donor, receiver, flow, transfer });
    }
    let mut mass = 0.0;
    let mut species = BTreeMap::new();
    let mut accounts = Vec::with_capacity(parcels.len());
    for (&id, old) in parcels {
        check(&mut cancelled)?;
        let new = &next[&id];
        let du = storage_changes.get(&id).map_or(0.0, |sum| sum.value);
        let dm = new.mass_kg() - old.mass_kg();
        mass += dm;
        if !du.is_finite() || !dm.is_finite() || !mass.is_finite() {
            return Err(Error::Invalid("gas network extensive balance is unrepresentable"));
        }
        let map_error = |error| Error::Parcel { id, error };
        let before = inventory(old.model(), old.mass_kg()).map_err(map_error)?;
        let mut after = inventory(new.model(), new.mass_kg()).map_err(map_error)?;
        // Subtract per cell before accumulating globally, avoiding cancellation
        // between large nearly equal whole-network inventories where possible.
        for (name, value) in before {
            check(&mut cancelled)?;
            let change = after.remove(&name).unwrap_or(0.0) - value;
            *species.entry(name).or_insert(0.0) += change;
        }
        for (name, value) in after { *species.entry(name).or_insert(0.0) += value; }
        accounts.push(GasNetworkParcelStep { id, mass_change_kg: dm, energy_change_j: du,
            gross_withdrawn_mass_kg: withdrawn.get(&id).copied().unwrap_or(0.0) });
    }
    let mut maximum_species: f64 = 0.0;
    for value in species.values() {
        if !value.is_finite() { return Err(Error::Invalid("gas network species balance is unrepresentable")); }
        maximum_species = maximum_species.max(value.abs());
    }
    if mass.abs() > controls.mass_tolerance_kg || maximum_species > controls.mass_tolerance_kg
        || energy.value.abs() > controls.energy_tolerance_j {
        return Err(Error::Invalid("gas network endpoints do not conserve total mass, each species and internal energy"));
    }
    check(&mut cancelled)?;
    Ok(GasNetworkCandidate { parcels: next, step: GasNetworkStep { duration_s, links: steps,
        parcels: accounts, total_mass_residual_kg: mass,
        maximum_species_mass_residual_kg: maximum_species, total_energy_residual_j: energy.value } })
}

#[cfg(test)]
mod tests;
