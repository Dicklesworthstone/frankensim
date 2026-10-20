//! Unsplit extensive transfer from one accepted material-network snapshot.
//!
//! All outlets carry their ORIGINAL donor composition and enthalpy. Combine
//! all inlets/outlets before reconstructing each affected volume exactly once.
//! Newly arrived species therefore cannot leave a downstream outlet in this
//! same explicit interval. No source law, caloric inversion or transport mixing
//! rule is duplicated: these remain the pair-transfer owner's implementations.
//!
//! Fixed volumes, well-mixed ideal gases and first-order explicit transport.
//! Prescribed transfers do not establish a pressure-flow or pumping-work law.
//! No reactions, vacuum, spatial jet momentum or entropy-production guarantee.
//! Balance convention: Cantera's Control Volume Reactor / Reactor Interactions.

use std::collections::BTreeMap;
use super::{ClosedGasParcel, GasEnergyError, GasEnergySolve, GasSpeciesTransfer,
    SpeciesId, endpoint, inventory, mixed_model_from_portions, positive, stored_change};

/// Hard bound on material volumes per proposal.
pub const MAX_SNAPSHOT_PARCELS: usize = 256;
/// Hard bound on submitted transfer records, including zero/parallel records.
pub const MAX_SNAPSHOT_TRANSFERS: usize = 4096;

/// Prescribed integrated mass, not a rate. Opposite directions are separate rows.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasSnapshotTransfer {
    /// Stable original upstream parcel id.
    pub donor: u64,
    /// Stable original receiving parcel id.
    pub receiver: u64,
    /// Nonnegative whole-interval upstream mass [kg].
    pub mass_kg: f64,
}
/// Absolute local AND global physical acceptance, not relative to stored energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasSnapshotSolve {
    /// Gross outgoing mass / ORIGINAL mass, in (0,1). Inflow cannot replenish it.
    pub maximum_withdrawn_fraction: f64,
    /// Total and individual species mass tolerance [kg].
    pub mass_tolerance_kg: f64,
    /// Extensive internal-energy tolerance [J].
    pub energy_tolerance_j: f64,
}
/// Any refusal returns no endpoints; accepted inputs are immutable.
#[derive(Clone, Debug, PartialEq)]
pub enum GasSnapshotError {
    /// Invalid inputs or unresolved extensive arithmetic.
    Invalid(&'static str),
    /// An original material/energy owner refused this parcel.
    Parcel { id: u64, error: GasEnergyError },
    /// Aggregate withdrawal exceeded the initial donor inventory allowance.
    Withdrawal { id: u64, fraction: f64, maximum: f64 },
    /// Cancelled before publishing the complete proposal.
    Cancelled,
}
impl core::fmt::Display for GasSnapshotError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "gas snapshot transfer: {self:?}")
    }
}
impl std::error::Error for GasSnapshotError {}

/// Parallel transfers combined in canonical (donor,receiver) order.
#[derive(Clone, Debug, PartialEq)]
pub struct GasSnapshotFlux {
    /// Original donor, receiver and aggregate mass.
    pub transfer: GasSnapshotTransfer,
    /// Energy carried once, with the original upstream enthalpy reference [J].
    pub enthalpy_j: f64,
    /// Original upstream species inventory, in canonical species order.
    pub species: Vec<GasSpeciesTransfer>,
}
/// Actual endpoint response to all incoming and outgoing fluxes.
#[derive(Clone, Debug, PartialEq)]
pub struct GasSnapshotParcelStep {
    /// Original parcel id.
    pub id: u64,
    /// Gross integrated inflow [kg].
    pub incoming_mass_kg: f64,
    /// Gross integrated outflow [kg].
    pub outgoing_mass_kg: f64,
    /// Actual change in extensive internal energy [J].
    pub energy_change_j: f64,
    /// Actual energy change minus the net upstream enthalpy flux [J].
    pub energy_residual_j: f64,
    /// Actual mass change minus the net prescribed mass flux [kg].
    pub mass_residual_kg: f64,
    /// Maximum local individual-species discrepancy [kg].
    pub maximum_species_mass_residual_kg: f64,
    /// Interior iterations of the single endpoint energy inversion (zero if idle).
    pub iterations: usize,
}
/// Global checks come from actual endpoints, not just cancelling flux records.
#[derive(Clone, Debug, PartialEq)]
pub struct GasSnapshotStep {
    /// Canonical upstream-state transfers.
    pub fluxes: Vec<GasSnapshotFlux>,
    /// Local acceptance accounts in parcel-id order.
    pub parcels: Vec<GasSnapshotParcelStep>,
    /// Sum of actual mass changes [kg].
    pub total_mass_residual_kg: f64,
    /// Largest global individual-species mass change [kg].
    pub maximum_species_mass_residual_kg: f64,
    /// Sum of actual internal-energy changes [J].
    pub total_energy_residual_j: f64,
}
/// Publish every endpoint together only after the surrounding coupling accepts.
#[derive(Clone, Debug, PartialEq)]
pub struct GasSnapshotCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    step: GasSnapshotStep,
}
impl GasSnapshotCandidate {
    /// Updated complete media retaining original source/pack bindings.
    #[must_use]
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> { &self.parcels }
    /// Checked local/global balances.
    #[must_use]
    pub const fn step(&self) -> &GasSnapshotStep { &self.step }
    /// Move all endpoints without repeating any material query.
    #[must_use]
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, GasSnapshotStep) {
        (self.parcels, self.step)
    }
}

use GasSnapshotError as Error;
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}
fn material(id: u64, error: GasEnergyError) -> Error {
    if error == GasEnergyError::Cancelled { Error::Cancelled } else { Error::Parcel { id, error } }
}
#[derive(Clone, Copy, Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), Error> {
        let adjusted = term - self.correction;
        let next = self.value + adjusted;
        let correction = (next - self.value) - adjusted;
        if !term.is_finite() || !adjusted.is_finite() || !next.is_finite() || !correction.is_finite() {
            return Err(Error::Invalid("snapshot extensive sum is unrepresentable"));
        }
        self.value = next;
        self.correction = correction;
        Ok(())
    }
}
#[derive(Default)]
struct Budget {
    incoming: Vec<(u64, f64)>,
    mass_in: Sum,
    mass_out: Sum,
    energy: Sum,
    // Sum dm*(h_upstream-u_old) avoids cancellation of formation-energy totals.
    relative_energy: Sum,
    species: BTreeMap<SpeciesId, Sum>,
}

/// Transfer all prescribed inlet/outlet masses from one immutable snapshot.
///
/// `endpoints` must provide exactly one original energy-inversion control per
/// parcel, including idle parcels. Parallel transfers are summed in canonical
/// numerical order before species/energy assembly. Reordering transfer records
/// cannot change results. No edge id participates in physical arithmetic.
///
/// Each active node is reconstructed once from its retained initial species
/// plus every original incoming species inventory. Sources must agree for
/// shared species. Gross donor bounds, local mass/species/energy balances and
/// global balances are checked before ANY endpoint is returned. Idle nodes
/// remain bitwise unchanged, and closed transfers do not join incompatible data.
/// No mass, energy, species or temperature is clipped on failure.
pub fn transfer_candidate(
    parcels: &BTreeMap<u64, ClosedGasParcel>, transfers: &[GasSnapshotTransfer],
    endpoints: &BTreeMap<u64, GasEnergySolve>, controls: GasSnapshotSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasSnapshotCandidate, Error> {
    poll(&mut cancelled)?;
    if parcels.is_empty() || parcels.len() > MAX_SNAPSHOT_PARCELS
        || transfers.len() > MAX_SNAPSHOT_TRANSFERS
        || endpoints.len() != parcels.len() || !parcels.keys().eq(endpoints.keys())
        || !positive(controls.maximum_withdrawn_fraction) || controls.maximum_withdrawn_fraction >= 1.0
        || !positive(controls.mass_tolerance_kg) || !positive(controls.energy_tolerance_j) {
        return Err(Error::Invalid("bounded nonempty gas graph, exact endpoint controls and positive physical allowances required"));
    }
    for (&id, solve) in endpoints {
        poll(&mut cancelled)?;
        solve.validate().map_err(|e| material(id, e))?;
    }
    for row in transfers {
        poll(&mut cancelled)?;
        if row.donor == row.receiver || !parcels.contains_key(&row.donor)
            || !parcels.contains_key(&row.receiver) || !row.mass_kg.is_finite() || row.mass_kg < 0.0 {
            return Err(Error::Invalid("snapshot transfer needs distinct existing parcels and finite nonnegative mass"));
        }
    }
    let mut rows = transfers.to_vec();
    rows.sort_by(|a, b| (a.donor, a.receiver).cmp(&(b.donor, b.receiver))
        .then_with(|| a.mass_kg.total_cmp(&b.mass_kg)));
    let mut grouped = BTreeMap::<(u64, u64), Sum>::new();
    for row in rows {
        if row.mass_kg > 0.0 { grouped.entry((row.donor, row.receiver)).or_default().add(row.mass_kg)?; }
    }
    let mut budgets: BTreeMap<u64, Budget> = parcels.keys().map(|&id| (id, Budget::default())).collect();
    let mut fluxes = Vec::with_capacity(grouped.len());
    for ((donor, receiver), mass) in grouped {
        poll(&mut cancelled)?;
        let dm = mass.value;
        if dm <= controls.mass_tolerance_kg {
            return Err(Error::Invalid("positive aggregate transfer is below absolute mass resolution"));
        }
        let source = &parcels[&donor];
        let h = source.state().specific_enthalpy_j_kg();
        let carried = dm * h;
        let species = inventory(source.model(), dm).map_err(|e| material(donor, e))?;
        for (id, sign) in [(donor, -1.0), (receiver, 1.0)] {
            let budget = budgets.get_mut(&id).expect("validated parcel key");
            if id == donor { budget.mass_out.add(dm)?; } else {
                budget.mass_in.add(dm)?;
                budget.incoming.push((donor, dm));
            }
            budget.energy.add(sign * carried)?;
            budget.relative_energy.add(sign * dm * (h - parcels[&id].state().specific_internal_energy_j_kg()))?;
            for (name, amount) in &species {
                budget.species.entry(name.clone()).or_default().add(sign * amount)?;
            }
        }
        fluxes.push(GasSnapshotFlux { transfer: GasSnapshotTransfer { donor, receiver, mass_kg: dm },
            enthalpy_j: carried, species: species.into_iter()
                .map(|(species, mass_kg)| GasSpeciesTransfer { species, mass_kg }).collect() });
    }
    // Reject aggregate depletion before resolving any new material endpoint.
    for (&id, budget) in &budgets {
        let fraction = budget.mass_out.value / parcels[&id].mass_kg();
        if !fraction.is_finite() || fraction > controls.maximum_withdrawn_fraction {
            return Err(Error::Withdrawal { id, fraction, maximum: controls.maximum_withdrawn_fraction });
        }
    }
    let mut next = BTreeMap::new();
    let mut accounts = Vec::with_capacity(parcels.len());
    let (mut total_mass, mut total_energy) = (Sum::default(), Sum::default());
    let mut global_species = BTreeMap::<SpeciesId, Sum>::new();
    for (&id, old) in parcels {
        poll(&mut cancelled)?;
        let budget = &budgets[&id];
        let (inflow, outflow) = (budget.mass_in.value, budget.mass_out.value);
        let retained = old.mass_kg() - outflow;
        let mass = retained + inflow;
        let mass_residual = (mass - old.mass_kg()) - (inflow - outflow);
        if !positive(retained) || !positive(mass) || !mass_residual.is_finite()
            || mass_residual.abs() > controls.mass_tolerance_kg
            || (outflow > 0.0 && retained == old.mass_kg()) || (inflow > 0.0 && mass == retained) {
            return Err(Error::Invalid("snapshot endpoint mass change is unresolved"));
        }
        let (new, iterations) = if inflow == 0.0 && outflow == 0.0 {
            (old.clone(), 0)
        } else {
            let model = if budget.incoming.is_empty() { old.model().clone() } else {
                let mut portions = vec![(old.model(), retained)];
                for &(donor, amount) in &budget.incoming { portions.push((parcels[&donor].model(), amount)); }
                mixed_model_from_portions(&portions, mass).map_err(|e| material(id, e))?
            };
            endpoint(old, model, mass, budget.relative_energy.value / mass,
                endpoints[&id], &mut cancelled).map_err(|e| material(id, e))?
        };
        let change = stored_change(old, &new);
        let energy_residual = change - budget.energy.value;
        if !change.is_finite() || !energy_residual.is_finite() || energy_residual.abs() > controls.energy_tolerance_j {
            return Err(Error::Invalid("snapshot endpoint failed local energy balance"));
        }
        let before = inventory(old.model(), old.mass_kg()).map_err(|e| material(id, e))?;
        let after = inventory(new.model(), new.mass_kg()).map_err(|e| material(id, e))?;
        let mut keys = std::collections::BTreeSet::new();
        keys.extend(before.keys().cloned()); keys.extend(after.keys().cloned());
        keys.extend(budget.species.keys().cloned());
        let mut maximum: f64 = 0.0;
        for name in keys {
            poll(&mut cancelled)?;
            let actual = after.get(&name).copied().unwrap_or(0.0) - before.get(&name).copied().unwrap_or(0.0);
            let expected = budget.species.get(&name).map_or(0.0, |s| s.value);
            let residual = actual - expected;
            if !residual.is_finite() || residual.abs() > controls.mass_tolerance_kg {
                return Err(Error::Invalid("snapshot endpoint failed individual species balance"));
            }
            maximum = maximum.max(residual.abs());
            global_species.entry(name).or_default().add(actual)?;
        }
        total_mass.add(new.mass_kg() - old.mass_kg())?;
        total_energy.add(change)?;
        accounts.push(GasSnapshotParcelStep { id, incoming_mass_kg: inflow, outgoing_mass_kg: outflow,
            energy_change_j: change, energy_residual_j: energy_residual, mass_residual_kg: mass_residual,
            maximum_species_mass_residual_kg: maximum, iterations });
        next.insert(id, new);
    }
    let maximum_species = global_species.values().map(|sum| sum.value.abs()).fold(0.0_f64, f64::max);
    if total_mass.value.abs() > controls.mass_tolerance_kg || maximum_species > controls.mass_tolerance_kg
        || total_energy.value.abs() > controls.energy_tolerance_j {
        return Err(Error::Invalid("snapshot endpoints failed global mass, species or energy balance"));
    }
    poll(&mut cancelled)?;
    Ok(GasSnapshotCandidate { parcels: next, step: GasSnapshotStep { fluxes, parcels: accounts,
        total_mass_residual_kg: total_mass.value, maximum_species_mass_residual_kg: maximum_species,
        total_energy_residual_j: total_energy.value } })
}

#[cfg(test)]
mod tests;
