//! Original fixed-composition calorics for one native Fourier interval.
//!
//! Each cell owns one source model for this interval. Trial queries reuse the
//! original safeguarded NASA inverse and its Cp sample; they never clone source
//! cards or resolve transport. Full source/card storage owned by each cell law
//! is additional to the compact native field arrays and source-query counters.

use super::super::{GasSpeciesEnthalpySource, Sum, caloric_error, invalid, poll, poll_every};
use crate::conjugate::{ChtError, TensorCaloricLaw};
use fs_conduction::enthalpy::advection::species::{CaloricComposition, CaloricEntropyAnchor};
use fs_conduction::enthalpy::{DepositedCell, DepositionError, ThermalLaw, ThermalStorageBasis};
use fs_exec::CancelGate;
use fs_material::thermochemical_gas::energy::GasEnergyError;
use fs_material::thermochemical_gas::energy::enthalpy::{GasEnthalpyCaloricState, GasEnthalpyLaw};
use fs_qty::Amount;
use std::cell::Cell;

/// One shared budget for all cells, nonlinear trials and native stages.
///
/// Counts direct forward/inverse caloric queries, including the two source
/// bracket endpoints and prescribed cell-temperature evaluations. Original
/// source-branch admission and opening/endpoint full gas preparations remain
/// separate bounded owners. Each inverse retains its source iteration limit.
pub(super) struct GasConductionQueryBudget<'g> {
    gate: &'g CancelGate,
    maximum_caloric_queries: usize,
    queries: Cell<usize>,
}

impl<'g> GasConductionQueryBudget<'g> {
    pub(super) fn new(
        maximum_caloric_queries: usize,
        gate: &'g CancelGate,
    ) -> Result<Self, ChtError> {
        poll(gate)?;
        if maximum_caloric_queries == 0 {
            return Err(invalid(
                "gas conduction needs a positive caloric-query budget",
            ));
        }
        Ok(Self {
            gate,
            maximum_caloric_queries,
            queries: Cell::new(0),
        })
    }

    pub(super) fn queries(&self) -> usize {
        self.queries.get()
    }

    fn poll(&self) -> Result<(), DepositionError> {
        poll(self.gate).map_err(|_| DepositionError::Cancelled)
    }

    fn consume(&self, count: usize) -> Result<(), DepositionError> {
        self.poll()?;
        let next = self
            .queries
            .get()
            .checked_add(count)
            .filter(|&next| next <= self.maximum_caloric_queries)
            .ok_or(DepositionError::InvalidInput(
                "gas conduction caloric-query budget exhausted",
            ))?;
        self.queries.set(next);
        Ok(())
    }
}

fn native_error(error: DepositionError) -> ChtError {
    match error {
        DepositionError::Cancelled => ChtError::Cancelled,
        other => invalid(format!("original gas conduction calorics refused: {other}")),
    }
}

fn thermal_error(error: GasEnergyError) -> DepositionError {
    match error {
        GasEnergyError::Cancelled => DepositionError::Cancelled,
        _ => DepositionError::MaterialDomain,
    }
}

/// Original mixture chart plus its exact accepted carrier binding.
///
/// The full canonical bank identity stays in GasEnthalpyLaw's composition and
/// enthalpy-reference descriptors even when its active model omits an absent
/// species. This is an isobaric fixed-carrier heat split, not sealed-volume Cv
/// heating. The native Fourier owner remains the only heat integrator.
pub(super) struct GasConductionCaloric<'q, 'g> {
    law: GasEnthalpyLaw,
    bounds_j_kg: [f64; 2],
    reference_mass_kg: f64,
    budget: &'q GasConductionQueryBudget<'g>,
    maximum_query_caloric_defect_j: Cell<f64>,
}

impl<'q, 'g> GasConductionCaloric<'q, 'g> {
    pub(super) fn new<const N: usize>(
        source: &GasSpeciesEnthalpySource<'_, N>,
        species_mass_kg: [f64; N],
        reference_mass_kg: f64,
        budget: &'q GasConductionQueryBudget<'g>,
    ) -> Result<Self, ChtError> {
        poll(budget.gate)?;
        if !reference_mass_kg.is_finite()
            || reference_mass_kg <= 0.0
            || species_mass_kg.iter().any(|m| !m.is_finite() || *m < 0.0)
        {
            return Err(invalid(
                "gas conduction needs finite positive carrier and nonnegative species",
            ));
        }
        let mut partition = Sum::default();
        partition.add(-reference_mass_kg);
        let mut amounts = [Amount::new(0.0); N];
        for (k, model) in source
            .source
            .thermodynamics()
            .components()
            .iter()
            .enumerate()
        {
            poll_every(budget.gate, k)?;
            partition.add(species_mass_kg[k]);
            let amount = species_mass_kg[k] / model.molar_mass().value();
            if !amount.is_finite() || (species_mass_kg[k] > 0.0 && amount == 0.0) {
                return Err(invalid(
                    "gas conduction loses a positive constituent molar amount",
                ));
            }
            amounts[k] = Amount::new(amount);
        }
        let relative_partition = partition.value.abs() / reference_mass_kg;
        if !relative_partition.is_finite() || relative_partition > 64.0 * N as f64 * f64::EPSILON {
            return Err(invalid(
                "gas conduction species must partition the carrier to machine roundoff",
            ));
        }
        budget.consume(2).map_err(native_error)?;
        let law = GasEnthalpyLaw::from_mole_amounts(
            source.source,
            &amounts,
            source.context.clone(),
            source.controls,
        )
        .map_err(caloric_error)?;
        let bounds_j_kg = law.enthalpy_bounds_j_kg().map_err(caloric_error)?;
        poll(budget.gate)?;
        Ok(Self {
            law,
            bounds_j_kg,
            reference_mass_kg,
            budget,
            maximum_query_caloric_defect_j: Cell::new(0.0),
        })
    }

    /// Largest actual native storage-inverse residual times this cell's mass.
    /// Rejected Newton probes consume work; their residual is a source-query
    /// diagnostic, never an external heat deposit or accepted physical history.
    pub(super) fn maximum_query_caloric_defect_j(&self) -> f64 {
        self.maximum_query_caloric_defect_j.get()
    }

    /// A prescribed Kelvin cell state has one original h on this increasing
    /// branch. Its source query is budgeted separately from inverse queries.
    pub(super) fn specific_enthalpy_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<f64, ChtError> {
        self.budget.consume(1).map_err(native_error)?;
        let h = self
            .law
            .enthalpy_at_temperature(temperature_k)
            .map_err(caloric_error)?;
        poll(self.budget.gate)?;
        Ok(h)
    }

    fn admit_energy(&self, energy: f64) -> Result<(), DepositionError> {
        if !energy.is_finite() || !(self.bounds_j_kg[0]..=self.bounds_j_kg[1]).contains(&energy) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(())
    }

    fn caloric(&self, energy: f64) -> Result<GasEnthalpyCaloricState, DepositionError> {
        self.budget.consume(1)?;
        self.admit_energy(energy)?;
        let state = self
            .law
            .caloric_at_enthalpy(energy, || poll(self.budget.gate).is_err())
            .map_err(thermal_error)?;
        let defect = self.reference_mass_kg * state.residual_j_kg.abs();
        if !defect.is_finite() {
            return Err(DepositionError::MaterialDomain);
        }
        self.maximum_query_caloric_defect_j
            .set(self.maximum_query_caloric_defect_j.get().max(defect));
        Ok(state)
    }
}

impl ThermalLaw for GasConductionCaloric<'_, '_> {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.law.storage_basis()
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        self.law.caloric_identity()
    }

    fn validate_cell(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
        self.law.validate_cell(cell)?;
        self.validate_reference_cell(cell.mass_kg, cell.geometry.volume_m3)
    }

    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        Ok(self.caloric(energy)?.temperature_k)
    }

    fn maximum_temperature_slope(&self) -> f64 {
        self.law.maximum_temperature_slope()
    }

    fn bulk_density_kg_m3(&self, energy: f64) -> Result<Option<f64>, DepositionError> {
        self.budget.consume(1)?;
        self.admit_energy(energy)?;
        let density = self.law.bulk_density_kg_m3(energy)?;
        self.budget.poll()?;
        Ok(density)
    }

    fn specific_entropy_change_j_kg_k(
        &self,
        initial: f64,
        final_energy: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.budget
            .consume(if initial == final_energy { 1 } else { 2 })?;
        self.admit_energy(initial)?;
        self.admit_energy(final_energy)?;
        let entropy = self
            .law
            .specific_entropy_change_j_kg_k(initial, final_energy)?;
        self.budget.poll()?;
        Ok(entropy)
    }

    fn advected_composition(&self) -> Option<CaloricComposition<'_>> {
        self.law.advected_composition()
    }

    fn advected_entropy_anchor(&self) -> Option<CaloricEntropyAnchor<'_>> {
        self.law.advected_entropy_anchor()
    }
}

impl TensorCaloricLaw for GasConductionCaloric<'_, '_> {
    fn specific_energy_bounds_j_kg(&self) -> [f64; 2] {
        self.bounds_j_kg
    }

    fn temperature_slope_kg_k_j(&self, energy: f64) -> Result<f64, DepositionError> {
        Ok(self.temperature_and_slope_at_specific_energy(energy)?[1])
    }

    fn temperature_and_slope_at_specific_energy(
        &self,
        energy: f64,
    ) -> Result<[f64; 2], DepositionError> {
        let state = self.caloric(energy)?;
        let slope = state.heat_capacity_j_kg_k.recip();
        if !slope.is_finite() || slope <= 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok([state.temperature_k, slope])
    }

    fn validate_reference_cell(
        &self,
        reference_mass_kg: f64,
        reference_volume_m3: f64,
    ) -> Result<(), DepositionError> {
        self.budget.poll()?;
        if reference_mass_kg != self.reference_mass_kg
            || !reference_volume_m3.is_finite()
            || reference_volume_m3 <= 0.0
        {
            return Err(DepositionError::InvalidInput(
                "gas conduction cell differs from its accepted carrier inventory",
            ));
        }
        Ok(())
    }
}
