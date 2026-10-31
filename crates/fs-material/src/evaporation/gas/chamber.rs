//! Simultaneous finite-rate exchange between retained liquid films and one gas.
//!
//! All surfaces use the SAME old gas for the original Stefan-film rate and
//! upstream vapor enthalpy. Each retains its original liquid enthalpy chart and
//! explicit caloric-datum binding. Gross condensation of each species shares
//! the old vapor inventory: evaporation elsewhere cannot finance an oversized
//! withdrawal in this interval. The original gas owner reconstructs one joint
//! composition, temperature, pressure and transport endpoint.
//!
//! Stable surface IDs give deterministic summation and mass-selection order.
//! Opposing streams retain source-relative sensible energy even when their net
//! mass is exactly zero. Every film checks equilibrium overshoot against the
//! common endpoint. Failed proposals and cancellation publish no partial state.
//!
//! This is the existing first-order frozen-rate, positive wet-film reduction,
//! with fixed gas volume and caller-owned liquid geometry. Original thermal
//! transport may cool/warm each retained film before an interval. This module
//! adds no sensible gas-film heat, nucleation, film exhaustion or hydrodynamics.

use std::collections::BTreeMap;

use super::{
    ClosedGasParcel, EvaporationError, EvaporationLaw, GasEvaporationError, GasEvaporationSolve,
    GasLiquidBinding, GasSpeciesBoundaryCandidate, LiquidSurfaceState, PreparedGasLiquidSurface,
    SpeciesBoundaryRelativeTransfer, SpeciesId, ThermochemicalGasModel, check_surface_endpoint,
    poll, prepare_surface, validate_controls,
};

/// One existing positive liquid inventory and its explicit original sources.
/// Identity describes a physical surface; duplicate IDs refuse before selection.
#[derive(Clone, Copy)]
pub struct GasLiquidChamberSurface<'a> {
    pub id: u64,
    pub law: &'a dyn EvaporationLaw,
    pub binding: &'a GasLiquidBinding,
    pub liquid: LiquidSurfaceState,
    pub temperature_from_enthalpy: &'a dyn Fn(f64) -> Result<f64, EvaporationError>,
    pub mass_transfer_coefficient_m_s: f64,
}

/// Original physical/caloric controls plus an explicit allocation/work ceiling.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasLiquidChamberSolve {
    pub surface: GasEvaporationSolve,
    pub maximum_surfaces: usize,
}

/// An actual retained liquid endpoint in canonical surface-ID order.
#[derive(Clone, Debug, PartialEq)]
pub struct GasLiquidChamberSurfaceStep {
    pub id: u64,
    pub liquid: LiquidSurfaceState,
    pub requested_transferred_mass_kg: f64,
    /// Positive evaporation; negative condensation.
    pub transferred_mass_kg: f64,
    pub surface_temperature_k: [f64; 2],
    pub vapor_partial_pressure_pa: [f64; 2],
    /// Signed diagnostic subaccount, already present in carried vapor energy.
    pub latent_energy_j: f64,
    pub energy_into_gas_j: f64,
    pub liquid_enthalpy_change_j: f64,
    /// Local original-liquid caloric defect. The common gas defect appears
    /// exactly once in the chamber report.
    pub energy_residual_j: f64,
    /// Actual local liquid mass change plus its selected signed transfer.
    pub mass_residual_kg: f64,
    pub recession_m: f64,
}

/// Shared original species inventory, including zero-net opposing streams.
#[derive(Clone, Debug, PartialEq)]
pub struct GasLiquidChamberSpeciesStep {
    pub species: SpeciesId,
    pub evaporated_mass_kg: f64,
    pub condensed_mass_kg: f64,
    pub transferred_mass_kg: f64,
    pub vapor_mass_kg: [f64; 2],
    pub mass_residual_kg: f64,
}

/// One whole-interval account, not a per-film multiplication of tolerances.
#[derive(Clone, Debug, PartialEq)]
pub struct GasLiquidChamberStep {
    pub duration_s: f64,
    pub transferred_mass_kg: f64,
    pub evaporated_mass_kg: f64,
    pub condensed_mass_kg: f64,
    pub latent_energy_j: f64,
    pub energy_into_gas_j: f64,
    pub liquid_enthalpy_change_j: f64,
    pub gas_internal_energy_change_j: f64,
    /// Stable signed sum of local liquid and the sole gas caloric defect.
    pub energy_residual_j: f64,
    /// Sum of absolute local liquid defects plus absolute gas defect.
    /// This, rather than cancellation of opposing defects, gates acceptance.
    pub sum_absolute_energy_residual_j: f64,
    pub mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub gas_temperature_k: [f64; 2],
    pub species: Vec<GasLiquidChamberSpeciesStep>,
}

/// All retained films and the original finite gas are accepted together.
#[derive(Clone, Debug, PartialEq)]
pub struct GasLiquidChamberCandidate<'catalog> {
    gas: GasSpeciesBoundaryCandidate<'catalog>,
    surfaces: Vec<GasLiquidChamberSurfaceStep>,
    step: GasLiquidChamberStep,
}

impl<'catalog> GasLiquidChamberCandidate<'catalog> {
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.gas.gas()
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.gas.catalog()
    }
    pub const fn gas_exchange(&self) -> &GasSpeciesBoundaryCandidate<'catalog> {
        &self.gas
    }
    pub fn surfaces(&self) -> &[GasLiquidChamberSurfaceStep] {
        &self.surfaces
    }
    pub const fn step(&self) -> &GasLiquidChamberStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        Vec<GasLiquidChamberSurfaceStep>,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasLiquidChamberStep,
    ) {
        let (gas, catalog, _) = self.gas.into_parts();
        (self.surfaces, gas, catalog, self.step)
    }
}

/// Deterministic compensated accumulation through the shared arithmetic owner.
#[derive(Default)]
struct Total {
    sum: f64,
    correction: f64,
}
impl Total {
    fn add(&mut self, value: f64) -> Result<(), GasEvaporationError> {
        let (sum, error) = fs_math::eft::two_sum(self.sum, value);
        let correction = self.correction + error;
        if ![value, sum, error, correction]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid("joint liquid/gas accumulation is unrepresentable"));
        }
        self.sum = sum;
        self.correction = correction;
        Ok(())
    }
    fn value(&self) -> Result<f64, GasEvaporationError> {
        finite(self.sum + self.correction)
    }
}

struct SpeciesAccount {
    old_mass_kg: f64,
    reference_h_j_kg: f64,
    evaporated: Total,
    condensed: Total,
    gross_withdrawal: Total,
    transferred: Total,
    relative_energy: Total,
    liquid_mass_change: Total,
}
impl SpeciesAccount {
    fn new(old_mass_kg: f64, reference_h_j_kg: f64) -> Self {
        Self {
            old_mass_kg,
            reference_h_j_kg,
            evaporated: Total::default(),
            condensed: Total::default(),
            gross_withdrawal: Total::default(),
            transferred: Total::default(),
            relative_energy: Total::default(),
            liquid_mass_change: Total::default(),
        }
    }
    fn add(&mut self, prepared: &PreparedGasLiquidSurface) -> Result<(), GasEvaporationError> {
        self.transferred.add(prepared.dm)?;
        self.relative_energy
            .add(prepared.relative_energy_into_gas_j)?;
        self.liquid_mass_change
            .add(prepared.liquid.mass_kg - prepared.original.mass_kg)?;
        if prepared.dm < 0.0 {
            self.condensed.add(-prepared.dm)?;
            self.gross_withdrawal
                .add(prepared.dm.abs().max(prepared.requested_dm.abs()))?;
        } else {
            self.evaporated.add(prepared.dm)?;
        }
        Ok(())
    }
}

fn invalid(message: &'static str) -> GasEvaporationError {
    GasEvaporationError::Invalid(message)
}
fn finite(value: f64) -> Result<f64, GasEvaporationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("joint liquid/gas account is unrepresentable"))
    }
}

/// Prepare every existing film against the same original gas and reconstruct
/// that gas once. Empty active rosters admit/preserve the original gas, allowing
/// a thermal timeline to wait for an authored film birth without inventing one.
pub fn advance_chamber_with_gas<'catalog>(
    surfaces: &[GasLiquidChamberSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasLiquidChamberSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasLiquidChamberCandidate<'catalog>, GasEvaporationError> {
    advance_chamber_with_gas_mass_selection(
        surfaces,
        gas,
        catalog,
        duration_s,
        controls,
        |_, mass| Ok(mass),
        cancelled,
    )
}

/// Same joint physical step, with one original native-mass selection per unique
/// surface in canonical ID order. The selector may choose only a nearby
/// representable increment under the original single-surface mass allowance.
/// Gross vapor admission considers BOTH the requested and selected withdrawals.
/// Callbacks should stage isolated state; their own external effects cannot be
/// rolled back if a later surface, joint caloric solve or cancellation refuses.
#[allow(clippy::too_many_arguments)]
pub fn advance_chamber_with_gas_mass_selection<'catalog>(
    surfaces: &[GasLiquidChamberSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasLiquidChamberSolve,
    mut select_mass: impl FnMut(u64, f64) -> Result<f64, GasEvaporationError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasLiquidChamberCandidate<'catalog>, GasEvaporationError> {
    poll(&mut cancelled)?;
    validate_controls(controls.surface, duration_s, 0.0)?;
    if controls.maximum_surfaces == 0 || surfaces.len() > controls.maximum_surfaces {
        return Err(invalid(
            "joint liquid surface count exceeds its positive work ceiling",
        ));
    }
    let mut ordered: Vec<_> = surfaces.iter().collect();
    ordered.sort_unstable_by_key(|surface| surface.id);
    if ordered.windows(2).any(|pair| pair[0].id == pair[1].id) {
        return Err(invalid("joint liquid surfaces require unique stable IDs"));
    }
    poll(&mut cancelled)?;
    let initial_temperature = gas.state().state().temperature;
    let mut accounts = BTreeMap::<SpeciesId, SpeciesAccount>::new();
    let mut prepared = Vec::with_capacity(ordered.len());
    for surface in &ordered {
        poll(&mut cancelled)?;
        let proposal = prepare_surface(
            surface.law,
            surface.binding,
            surface.liquid,
            gas,
            catalog,
            surface.temperature_from_enthalpy,
            surface.mass_transfer_coefficient_m_s,
            duration_s,
            controls.surface,
            |mass| select_mass(surface.id, mass),
            &mut cancelled,
        )?;
        let species = &surface.binding.species;
        if !accounts.contains_key(species) {
            accounts.insert(
                species.clone(),
                SpeciesAccount::new(
                    proposal.old_vapor_mass_kg,
                    catalog.species_specific_enthalpy_j_kg(species, initial_temperature)?,
                ),
            );
        }
        accounts
            .get_mut(species)
            .expect("inserted species")
            .add(&proposal)?;
        prepared.push(proposal);
    }
    let mut transfers = Vec::with_capacity(accounts.len());
    for (species, account) in &accounts {
        poll(&mut cancelled)?;
        let withdrawal = account.gross_withdrawal.value()?;
        if withdrawal > 0.0
            && (account.old_mass_kg <= 0.0
                || withdrawal / account.old_mass_kg > controls.surface.maximum_inventory_fraction)
        {
            let maximum_duration = duration_s
                * ((controls.surface.maximum_inventory_fraction * account.old_mass_kg)
                    / withdrawal);
            return Err(GasEvaporationError::StepSizeLimit {
                maximum_duration_s: maximum_duration.next_down().max(0.0),
                reason: "shared initial vapor inventory withdrawal limit",
            });
        }
        // Keep the row when net mass vanishes: opposing hot/cold streams can
        // still transport energy through this original species.
        transfers.push(SpeciesBoundaryRelativeTransfer {
            species: species.clone(),
            mass_into_gas_kg: account.transferred.value()?,
            relative_energy_into_gas_j: account.relative_energy.value()?,
        });
    }
    poll(&mut cancelled)?;
    let endpoint = gas.exchange_species_boundary_relative(
        catalog,
        &transfers,
        controls.surface.gas,
        &mut cancelled,
    )?;
    let mut reports = Vec::with_capacity(prepared.len());
    let mut liquid_mass_change = Total::default();
    let mut absolute_mass_defect = Total::default();
    let mut liquid_energy = Total::default();
    let mut latent_energy = Total::default();
    let mut energy_defect = Total::default();
    let mut absolute_energy_defect = Total::default();
    absolute_mass_defect.add(endpoint.step().mass_residual_kg.abs())?;
    energy_defect.add(endpoint.step().energy_residual_j)?;
    absolute_energy_defect.add(endpoint.step().energy_residual_j.abs())?;
    for (surface, proposal) in ordered.iter().zip(prepared) {
        poll(&mut cancelled)?;
        let final_pressure = check_surface_endpoint(
            &proposal,
            surface.law,
            surface.binding,
            endpoint.gas(),
            surface.mass_transfer_coefficient_m_s,
            duration_s,
            controls.surface,
        )?;
        let account = &accounts[&surface.binding.species];
        let carried_energy =
            finite(proposal.dm * account.reference_h_j_kg + proposal.relative_energy_into_gas_j)?;
        liquid_mass_change.add(proposal.liquid.mass_kg - proposal.original.mass_kg)?;
        absolute_mass_defect.add(proposal.liquid_mass_residual_kg.abs())?;
        liquid_energy.add(proposal.liquid_enthalpy_change_j)?;
        latent_energy.add(proposal.latent_energy_j)?;
        energy_defect.add(proposal.liquid_energy_residual_j)?;
        absolute_energy_defect.add(proposal.liquid_energy_residual_j.abs())?;
        reports.push(GasLiquidChamberSurfaceStep {
            id: surface.id,
            liquid: proposal.liquid,
            requested_transferred_mass_kg: proposal.requested_dm,
            transferred_mass_kg: proposal.dm,
            surface_temperature_k: proposal.surface_temperature_k,
            vapor_partial_pressure_pa: [proposal.vapor_partial_pressure_pa, final_pressure],
            latent_energy_j: proposal.latent_energy_j,
            energy_into_gas_j: carried_energy,
            liquid_enthalpy_change_j: proposal.liquid_enthalpy_change_j,
            energy_residual_j: proposal.liquid_energy_residual_j,
            mass_residual_kg: proposal.liquid_mass_residual_kg,
            recession_m: proposal.recession_m,
        });
    }
    let mut evaporated = Total::default();
    let mut condensed = Total::default();
    let mut species_reports = Vec::with_capacity(accounts.len());
    let mut maximum_species_mass_residual = endpoint.step().maximum_species_mass_residual_kg;
    for (species, account) in accounts {
        poll(&mut cancelled)?;
        let final_mass = endpoint.gas().species_mass_kg(&species)?;
        let mass_residual =
            finite(account.liquid_mass_change.value()? + (final_mass - account.old_mass_kg))?;
        maximum_species_mass_residual = maximum_species_mass_residual.max(mass_residual.abs());
        let evaporation = account.evaporated.value()?;
        let condensation = account.condensed.value()?;
        evaporated.add(evaporation)?;
        condensed.add(condensation)?;
        species_reports.push(GasLiquidChamberSpeciesStep {
            species,
            evaporated_mass_kg: evaporation,
            condensed_mass_kg: condensation,
            transferred_mass_kg: account.transferred.value()?,
            vapor_mass_kg: [account.old_mass_kg, final_mass],
            mass_residual_kg: mass_residual,
        });
    }
    let mass_residual =
        finite(liquid_mass_change.value()? + (endpoint.gas().mass_kg() - gas.mass_kg()))?;
    if mass_residual.abs() > controls.surface.gas.mass_tolerance_kg
        || maximum_species_mass_residual > controls.surface.gas.mass_tolerance_kg
        || absolute_mass_defect.value()? > controls.surface.gas.mass_tolerance_kg
    {
        return Err(invalid(
            "joint liquid/gas inventory fails its whole-interval mass allowance",
        ));
    }
    let sum_absolute_energy_residual = absolute_energy_defect.value()?;
    if sum_absolute_energy_residual > controls.surface.gas.energy_tolerance_j {
        return Err(GasEvaporationError::EnergyBalance {
            residual_j: sum_absolute_energy_residual,
            tolerance_j: controls.surface.gas.energy_tolerance_j,
        });
    }
    let step = GasLiquidChamberStep {
        duration_s,
        transferred_mass_kg: endpoint.step().mass_into_gas_kg,
        evaporated_mass_kg: evaporated.value()?,
        condensed_mass_kg: condensed.value()?,
        latent_energy_j: latent_energy.value()?,
        energy_into_gas_j: endpoint.step().energy_into_gas_j,
        liquid_enthalpy_change_j: liquid_energy.value()?,
        gas_internal_energy_change_j: endpoint.step().stored_energy_change_j,
        energy_residual_j: energy_defect.value()?,
        sum_absolute_energy_residual_j: sum_absolute_energy_residual,
        mass_residual_kg: mass_residual,
        maximum_species_mass_residual_kg: maximum_species_mass_residual,
        gas_temperature_k: [
            initial_temperature,
            endpoint.gas().state().state().temperature,
        ],
        species: species_reports,
    };
    poll(&mut cancelled)?;
    Ok(GasLiquidChamberCandidate {
        gas: endpoint,
        surfaces: reports,
        step,
    })
}

#[cfg(test)]
mod tests;
