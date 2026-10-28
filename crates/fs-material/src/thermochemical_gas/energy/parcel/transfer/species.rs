//! Selective species exchange between two original, fixed-volume gas states.
//!
//! Each signed species flux carries its own upstream partial specific enthalpy,
//! not the donor MIXTURE enthalpy. Opposite species can cross simultaneously.
//! Reconstruct each endpoint once with the original caloric/EOS/transport owners.
//! No rate law, membrane work, diffusion model or entropy guarantee is implied by
//! prescribed fluxes. The physical transport owner must provide those separately.

/// Conservative prescribed streams between the original gas and an external phase.
pub mod boundary;
/// Source-backed multicomponent molecular mixing and conservative pair advancement.
pub mod diffusion;
/// Reciprocal free-molecular transport through qualified thin apertures.
pub mod effusion;
/// Shared-snapshot selective transport on branching and cyclic volume networks.
pub mod network;
mod entropy;

use super::{ClosedGasParcel, Composition, CompositionBasis, GasEnergyError, GasEnergySolve,
    GasMassTransferSolve, SpeciesId, ThermochemicalGasModel, MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1,
    endpoint, invalid, inventory, mixed_model_from_portions, poll, positive, source_error, stored_change};
use std::collections::BTreeMap;
use fs_qty::Temperature;

/// Integrated species mass, positive left -> right; negative right -> left.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesMassTransfer {
    pub species: SpeciesId,
    pub mass_left_to_right_kg: f64,
}
/// Original endpoint controls plus bounds on actual individual donor inventories.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasSpeciesExchangeSolve {
    /// Existing [left, right] caloric reconstruction controls.
    pub gas: [GasEnergySolve; 2],
    /// Largest withdrawal / initial mass of ANY consumed species, in (0,1).
    pub maximum_species_fraction: f64,
    /// Absolute local/global total and individual-species mass tolerance [kg].
    pub mass_tolerance_kg: f64,
    /// Absolute local/global extensive internal-energy tolerance [J].
    pub energy_tolerance_j: f64,
}
impl GasSpeciesExchangeSolve {
    fn validate(self) -> Result<(), GasEnergyError> {
        GasMassTransferSolve { donor: self.gas[0], receiver: self.gas[1],
            mass_tolerance_kg: self.mass_tolerance_kg, energy_tolerance_j: self.energy_tolerance_j }.validate()?;
        if !positive(self.maximum_species_fraction) || self.maximum_species_fraction >= 1.0 {
            return Err(invalid("species withdrawal fraction must lie in (0,1)"));
        }
        Ok(())
    }
}
/// Original signed flux and the actual source partial enthalpy it transported.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesEnergyTransfer {
    pub transfer: SpeciesMassTransfer,
    pub source_enthalpy_j_kg: f64,
    pub enthalpy_left_to_right_j: f64,
}
/// Checked endpoint accounts; chemical formation energy is not an extra heat.
#[derive(Clone, Debug, PartialEq)]
pub struct GasSpeciesExchangeStep {
    pub species: Vec<SpeciesEnergyTransfer>,
    pub net_mass_left_to_right_kg: f64,
    pub enthalpy_left_to_right_j: f64,
    pub energy_change_j: [f64; 2],
    pub mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub energy_residual_j: f64,
    pub iterations: usize,
}
/// Both proposed states, or neither. No accepted input is changed by a trial.
#[derive(Clone, Debug, PartialEq)]
pub struct GasSpeciesExchangeCandidate {
    left: ClosedGasParcel,
    right: ClosedGasParcel,
    step: GasSpeciesExchangeStep,
}
impl GasSpeciesExchangeCandidate {
    pub const fn left(&self) -> &ClosedGasParcel { &self.left }
    pub const fn right(&self) -> &ClosedGasParcel { &self.right }
    pub const fn step(&self) -> &GasSpeciesExchangeStep { &self.step }
    pub fn into_parts(self) -> (ClosedGasParcel, ClosedGasParcel, GasSpeciesExchangeStep) {
        (self.left, self.right, self.step)
    }
}

#[derive(Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), GasEnergyError> {
        let adjusted = term - self.correction;
        let value = self.value + adjusted;
        let correction = (value - self.value) - adjusted;
        if ![term, adjusted, value, correction].iter().all(|x| x.is_finite()) {
            return Err(invalid("selective gas exchange arithmetic is unrepresentable"));
        }
        self.value = value; self.correction = correction;
        Ok(())
    }
}

// The existing source join owns consistency, transport, water and pack identity.
// Its composition is a catalog join only; no fictitious gas state is evaluated.
fn catalog(left: &ClosedGasParcel, right: &ClosedGasParcel) -> Result<ThermochemicalGasModel, GasEnergyError> {
    mixed_model_from_portions(&[(left.model(), left.mass_kg()), (right.model(), right.mass_kg())],
        left.mass_kg() + right.mass_kg())
}
fn model_for_inventory(
    catalog: &ThermochemicalGasModel, target: &BTreeMap<SpeciesId, f64>,
) -> Result<ThermochemicalGasModel, GasEnergyError> {
    let mut amounts = Vec::with_capacity(catalog.thermodynamics().components().len());
    for species in catalog.thermodynamics().components() {
        let mass = target.get(species.species()).copied().unwrap_or(0.0);
        let amount = mass / species.molar_mass().value();
        if !mass.is_finite() || mass < 0.0 || !amount.is_finite() || (mass > 0.0 && amount <= 0.0) {
            return Err(invalid("exchanged species inventory is unrepresentable"));
        }
        amounts.push(amount);
    }
    let total: f64 = amounts.iter().sum();
    if !positive(total) { return Err(invalid("selective exchange cannot produce a vacuum")); }
    let mut fractions: Vec<_> = amounts.iter().map(|n| n / total).collect();
    let sum: f64 = fractions.iter().sum();
    if !sum.is_finite() || (sum - 1.0).abs() > 256.0 * f64::EPSILON {
        return Err(invalid("selective composition exceeds normalization roundoff"));
    }
    if sum != 1.0 {
        let last = amounts.iter().rposition(|n| *n > 0.0).ok_or(invalid("empty exchanged composition"))?;
        fractions[last] = 1.0 - fractions[..last].iter().sum::<f64>();
    }
    if amounts.iter().zip(&fractions).any(|(&n, &x)| n > 0.0 && !positive(x)) {
        return Err(invalid("selective exchange would discard a positive trace species"));
    }
    let composition = Composition::new(CompositionBasis::MoleFraction, fractions).map_err(source_error)?;
    catalog.with_composition(composition).map_err(GasEnergyError::from)
}

// Shared by pair and network exchange: selective transport carries partial h_i,
// not the mixture-average enthalpy. The original NASA-9 owner supplies the law.
fn source_enthalpy(source: &ClosedGasParcel, id: &SpeciesId) -> Result<f64, GasEnergyError> {
    let species = source.model().thermodynamics().components().iter()
        .find(|s| s.species() == id).ok_or(invalid("species absent from its declared upstream"))?;
    let properties = species.evaluate_properties(Temperature::new(source.state().state().temperature)).map_err(source_error)?;
    let h = properties.h().value() / species.molar_mass().value();
    if !h.is_finite() { return Err(invalid("source species enthalpy is unrepresentable")); }
    Ok(h)
}

impl ClosedGasParcel {
    /// Apply all prescribed selective species transfers to the SAME pair snapshot.
    /// Duplicate/unknown species, vanished increments, depleted donors, conflicting
    /// source laws, failed energy balances and cancellation refuse the whole pair.
    /// Rows are canonicalized by original SpeciesId; caller row order is immaterial.
    ///
    /// A species carries h_i(T_upstream) including flow work/formation references.
    /// No jet energy or chemical heat may be added a second time. Selective
    /// pumping or membrane separation may require additional work not modeled here.
    /// Endpoint volumes and non-temperature/pressure material axes stay unchanged.
    pub fn exchange_species_with(
        &self, right: &Self, transfers: &[SpeciesMassTransfer], controls: GasSpeciesExchangeSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasSpeciesExchangeCandidate, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if transfers.len() > MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1 {
            return Err(invalid("selective exchange exceeds the shared species bound"));
        }
        let mut ordered = BTreeMap::new();
        for row in transfers {
            poll(&mut cancelled)?;
            if !row.mass_left_to_right_kg.is_finite()
                || ordered.insert(row.species.clone(), row.mass_left_to_right_kg).is_some()
                || ![self, right].iter().any(|p| p.model().thermodynamics().components()
                    .iter().any(|s| s.species() == &row.species)) {
                return Err(invalid("finite unique transfers naming an existing species required"));
            }
        }
        if ordered.values().all(|dm| *dm == 0.0) {
            poll(&mut cancelled)?;
            return Ok(GasSpeciesExchangeCandidate { left: self.clone(), right: right.clone(),
                step: GasSpeciesExchangeStep { species: Vec::new(), net_mass_left_to_right_kg: 0.0,
                    enthalpy_left_to_right_j: 0.0, energy_change_j: [0.0; 2], mass_residual_kg: 0.0,
                    maximum_species_mass_residual_kg: 0.0, energy_residual_j: 0.0, iterations: 0 } });
        }
        if std::ptr::eq(self, right) {
            return Err(invalid("nonzero species exchange needs two distinct physical gas owners"));
        }
        let catalog = catalog(self, right)?;
        let before = [inventory(self.model(), self.mass_kg())?, inventory(right.model(), right.mass_kg())?];
        let mut targets = before.clone();
        let mut mass = Sum::default(); let mut energy = Sum::default();
        let mut relative = [Sum::default(), Sum::default()];
        let mut rows = Vec::with_capacity(transfers.len());
        for (id, &dm) in &ordered {
            poll(&mut cancelled)?;
            if dm == 0.0 { continue; }
            let donor = usize::from(dm < 0.0);
            let present = before[donor].get(id).copied().unwrap_or(0.0);
            if !positive(present) || dm.abs() / present > controls.maximum_species_fraction {
                return Err(invalid("selective transfer exceeds an original species inventory allowance"));
            }
            for (side, sign) in [(0, -1.0), (1, 1.0)] {
                let old = before[side].get(id).copied().unwrap_or(0.0);
                let next = old + sign * dm;
                if !positive(next) || next == old {
                    return Err(invalid("selective species increment vanished or exhausted its inventory"));
                }
                targets[side].insert(id.clone(), next);
            }
            let source = [self, right][donor];
            let h = source_enthalpy(source, id)?;
            let carried = dm * h;
            mass.add(dm)?; energy.add(carried)?;
            relative[0].add(-dm * (h - self.state().specific_internal_energy_j_kg()))?;
            relative[1].add(dm * (h - right.state().specific_internal_energy_j_kg()))?;
            rows.push(SpeciesEnergyTransfer { transfer: SpeciesMassTransfer {
                species: id.clone(), mass_left_to_right_kg: dm }, source_enthalpy_j_kg: h,
                enthalpy_left_to_right_j: carried });
        }
        let masses = [self.mass_kg() - mass.value, right.mass_kg() + mass.value];
        if masses.iter().any(|m| !positive(*m)) {
            return Err(invalid("exchanged total gas mass is unrepresentable"));
        }
        let (left, li) = endpoint(self, model_for_inventory(&catalog, &targets[0])?, masses[0],
            relative[0].value / masses[0], controls.gas[0], &mut cancelled)?;
        let (next_right, ri) = endpoint(right, model_for_inventory(&catalog, &targets[1])?, masses[1],
            relative[1].value / masses[1], controls.gas[1], &mut cancelled)?;
        let du = [stored_change(self, &left), stored_change(right, &next_right)];
        let er = du[0] + du[1];
        let dm = [left.mass_kg() - self.mass_kg(), next_right.mass_kg() - right.mass_kg()];
        let mr = dm[0] + dm[1];
        if ![du[0], du[1], er, mr].iter().all(|v| v.is_finite())
            || (du[0] + energy.value).abs() > controls.energy_tolerance_j
            || (du[1] - energy.value).abs() > controls.energy_tolerance_j || er.abs() > controls.energy_tolerance_j
            || (dm[0] + mass.value).abs() > controls.mass_tolerance_kg
            || (dm[1] - mass.value).abs() > controls.mass_tolerance_kg || mr.abs() > controls.mass_tolerance_kg {
            return Err(invalid("selective gas exchange failed local/global mass or energy balance"));
        }
        let actual = [inventory(left.model(), left.mass_kg())?, inventory(next_right.model(), next_right.mass_kg())?];
        let mut maximum: f64 = 0.0;
        for species in catalog.thermodynamics().components() {
            poll(&mut cancelled)?;
            let id = species.species();
            let mut global = 0.0;
            for side in 0..2 {
                let after = actual[side].get(id).copied().unwrap_or(0.0);
                let expected = targets[side].get(id).copied().unwrap_or(0.0);
                maximum = maximum.max((after - expected).abs());
                global += after - before[side].get(id).copied().unwrap_or(0.0);
            }
            if !global.is_finite() { return Err(invalid("selective species balance overflow")); }
            maximum = maximum.max(global.abs());
        }
        if maximum > controls.mass_tolerance_kg {
            return Err(invalid("selective gas endpoints did not retain individual species inventories"));
        }
        poll(&mut cancelled)?;
        Ok(GasSpeciesExchangeCandidate { left, right: next_right,
            step: GasSpeciesExchangeStep { species: rows, net_mass_left_to_right_kg: mass.value,
                enthalpy_left_to_right_j: energy.value, energy_change_j: du, mass_residual_kg: mr,
                maximum_species_mass_residual_kg: maximum, energy_residual_j: er, iterations: li + ri } })
    }
}

#[cfg(test)]
mod tests;
