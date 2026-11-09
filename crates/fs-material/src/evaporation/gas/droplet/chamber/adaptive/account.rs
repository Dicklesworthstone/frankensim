//! Accepted fine intervals and independent original-to-final inventories.
//!
//! Original caloric defects are accumulated in their source-relative coordinates.
//! Absolute formation-energy totals are never subtracted to manufacture a
//! window energy residual. Coarse and rejected candidates never enter Account.

use std::collections::BTreeMap;

use super::{AdaptiveGasDropletControls, GasDropletSeed};
use crate::evaporation::gas::{
    GasEvaporationError, chamber::Total, droplet::chamber::GasDropletChamberCandidate, poll,
};
use crate::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_thermochem::SpeciesId;

/// Accepted physical transfers for one original droplet, in stable-ID order.
#[derive(Clone, Debug, PartialEq)]
pub struct GasDropletEvolutionSurfaceAudit {
    pub id: u64,
    pub mass_kg: [f64; 2],
    /// Positive evaporation; negative condensation.
    pub transferred_mass_kg: f64,
    pub evaporated_mass_kg: f64,
    pub condensed_mass_kg: f64,
    /// Diagnostic subaccount already included in vapor-carried energy.
    pub latent_energy_j: f64,
    /// Vapor stream only, excluding reciprocal sensible heat.
    pub energy_into_gas_j: f64,
    pub sensible_heat_into_liquid_j: f64,
    /// Original-to-final change in the liquid's admitted gas-compatible datum.
    pub liquid_enthalpy_change_j: f64,
    /// Actual final-minus-initial liquid mass plus the accepted signed transfer.
    pub mass_residual_kg: f64,
    pub accumulated_mass_residual_kg: f64,
    /// Local phase caloric defect. Joint heat defects remain in the window sum.
    pub phase_energy_residual_j: f64,
}

/// Every original catalog species, including untransferred carrier inventories.
#[derive(Clone, Debug, PartialEq)]
pub struct GasDropletEvolutionSpeciesAudit {
    pub species: SpeciesId,
    pub gas_mass_kg: [f64; 2],
    pub transferred_mass_kg: f64,
    pub evaporated_mass_kg: f64,
    pub condensed_mass_kg: f64,
    /// Actual gas inventory change minus accepted transfer into that species.
    pub mass_residual_kg: f64,
    /// Actual gas change plus changes of all liquids bound to this species.
    pub closed_mass_residual_kg: f64,
}

/// One physical accounting window, independent of its adaptive partition.
#[derive(Clone, Debug, PartialEq)]
pub struct GasDropletEvolutionAudit {
    pub droplets: Vec<GasDropletEvolutionSurfaceAudit>,
    pub species: Vec<GasDropletEvolutionSpeciesAudit>,
    pub transferred_mass_kg: f64,
    pub evaporated_mass_kg: f64,
    pub condensed_mass_kg: f64,
    pub latent_energy_j: f64,
    /// Vapor-carried energy; latent energy is already part of it.
    pub energy_into_gas_j: f64,
    pub sensible_heat_into_liquid_j: f64,
    pub liquid_enthalpy_change_j: f64,
    pub gas_internal_energy_change_j: f64,
    /// Independent actual final gas/liquid mass change.
    pub mass_residual_kg: f64,
    /// Independent actual gas mass change minus the accepted transfer.
    pub gas_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    /// Signed sum of the original accepted interval defects.
    pub accumulated_mass_residual_kg: f64,
    /// Larger of accumulated absolute interval defects and independent final
    /// local inventory defects. One allowance bounds the whole window.
    pub sum_absolute_mass_residual_kg: f64,
    /// Signed original source-relative caloric defects, including joint heat.
    pub energy_residual_j: f64,
    pub sum_absolute_energy_residual_j: f64,
}

struct SurfaceTotals {
    id: u64,
    species: SpeciesId,
    initial_mass_kg: f64,
    transferred: Total,
    evaporated: Total,
    condensed: Total,
    latent: Total,
    vapor: Total,
    sensible: Total,
    liquid_energy: Total,
    mass_defect: Total,
    phase_energy_defect: Total,
}

#[derive(Default)]
struct SpeciesTotals {
    initial_mass_kg: f64,
    transferred: Total,
    evaporated: Total,
    condensed: Total,
    liquid_change: Total,
}

pub(super) struct Account {
    surfaces: Vec<SurfaceTotals>,
    species: BTreeMap<SpeciesId, SpeciesTotals>,
    initial_gas_mass_kg: f64,
    gas_volume_m3: f64,
    transferred: Total,
    evaporated: Total,
    condensed: Total,
    latent: Total,
    vapor: Total,
    sensible: Total,
    liquid_energy: Total,
    gas_energy: Total,
    mass_defect: Total,
    absolute_mass_defect: Total,
    energy_defect: Total,
    absolute_energy_defect: Total,
}

impl Account {
    /// Called after the original zero-duration physical/source admission.
    pub(super) fn new(
        seeds: &[GasDropletSeed<'_>],
        gas: &ClosedGasParcel,
        catalog: &ThermochemicalGasModel,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, GasEvaporationError> {
        poll(cancelled)?;
        if seeds.windows(2).any(|pair| pair[0].id >= pair[1].id) {
            return Err(invalid(
                "droplet accounting needs canonical unique identities",
            ));
        }
        let initial = gas.species_inventory_kg()?;
        let mut species = BTreeMap::new();
        for component in catalog.thermodynamics().components() {
            poll(cancelled)?;
            let id = component.species();
            species.insert(
                id.clone(),
                SpeciesTotals {
                    initial_mass_kg: initial.get(id).copied().unwrap_or(0.0),
                    ..SpeciesTotals::default()
                },
            );
        }
        if initial.keys().any(|id| !species.contains_key(id)) {
            return Err(invalid(
                "initial droplet atmosphere leaves its original catalog",
            ));
        }
        let mut surfaces = Vec::with_capacity(seeds.len());
        for seed in seeds {
            poll(cancelled)?;
            if !species.contains_key(&seed.binding.species) {
                return Err(invalid(
                    "droplet volatile is absent from the original catalog",
                ));
            }
            surfaces.push(SurfaceTotals {
                id: seed.id,
                species: seed.binding.species.clone(),
                initial_mass_kg: seed.droplet.mass_kg,
                transferred: Total::default(),
                evaporated: Total::default(),
                condensed: Total::default(),
                latent: Total::default(),
                vapor: Total::default(),
                sensible: Total::default(),
                liquid_energy: Total::default(),
                mass_defect: Total::default(),
                phase_energy_defect: Total::default(),
            });
        }
        Ok(Self {
            surfaces,
            species,
            initial_gas_mass_kg: gas.mass_kg(),
            gas_volume_m3: gas.volume_m3(),
            transferred: Total::default(),
            evaporated: Total::default(),
            condensed: Total::default(),
            latent: Total::default(),
            vapor: Total::default(),
            sensible: Total::default(),
            liquid_energy: Total::default(),
            gas_energy: Total::default(),
            mass_defect: Total::default(),
            absolute_mass_defect: Total::default(),
            energy_defect: Total::default(),
            absolute_energy_defect: Total::default(),
        })
    }

    pub(super) fn add(
        &mut self,
        endpoint: &GasDropletChamberCandidate<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), GasEvaporationError> {
        poll(cancelled)?;
        if endpoint.droplets().len() != self.surfaces.len()
            || endpoint.gas().volume_m3().to_bits() != self.gas_volume_m3.to_bits()
        {
            return Err(invalid(
                "accepted droplet interval changed its roster or gas volume",
            ));
        }
        let report = endpoint.report();
        let part = &report.exchange;
        if !part.duration_s.is_finite() || part.duration_s <= 0.0 {
            return Err(invalid(
                "only positive accepted intervals enter droplet accounting",
            ));
        }
        for (total, row) in self.surfaces.iter_mut().zip(endpoint.droplets()) {
            poll(cancelled)?;
            if total.id != row.id {
                return Err(invalid(
                    "accepted droplet interval changed a stable identity",
                ));
            }
            let phase = &row.exchange;
            total.transferred.add(phase.transferred_mass_kg)?;
            total.evaporated.add(phase.transferred_mass_kg.max(0.0))?;
            total.condensed.add((-phase.transferred_mass_kg).max(0.0))?;
            total.latent.add(phase.latent_energy_j)?;
            total.vapor.add(phase.energy_into_gas_j)?;
            total
                .sensible
                .add(row.heat.map_or(0.0, |heat| heat.heat_into_liquid_j))?;
            total.liquid_energy.add(phase.liquid_enthalpy_change_j)?;
            total.mass_defect.add(phase.mass_residual_kg)?;
            total.phase_energy_defect.add(row.phase.energy_residual_j)?;
        }
        for row in &part.species {
            poll(cancelled)?;
            let total = self.species.get_mut(&row.species).ok_or_else(|| {
                invalid("accepted droplet interval introduced an unknown gas species")
            })?;
            total.transferred.add(row.transferred_mass_kg)?;
            total.evaporated.add(row.evaporated_mass_kg)?;
            total.condensed.add(row.condensed_mass_kg)?;
        }
        self.transferred.add(part.transferred_mass_kg)?;
        self.evaporated.add(part.evaporated_mass_kg)?;
        self.condensed.add(part.condensed_mass_kg)?;
        self.latent.add(part.latent_energy_j)?;
        self.vapor.add(part.energy_into_gas_j)?;
        self.sensible
            .add(report.heat.map_or(0.0, |heat| heat.heat_into_liquid_j))?;
        self.liquid_energy.add(part.liquid_enthalpy_change_j)?;
        self.gas_energy.add(part.gas_internal_energy_change_j)?;
        self.mass_defect.add(part.mass_residual_kg)?;
        self.absolute_mass_defect
            .add(report.sum_absolute_mass_residual_kg)?;
        self.energy_defect.add(part.energy_residual_j)?;
        self.absolute_energy_defect
            .add(part.sum_absolute_energy_residual_j)?;
        Ok(())
    }

    pub(super) fn gate(
        &self,
        controls: AdaptiveGasDropletControls,
    ) -> Result<(), GasEvaporationError> {
        admit_mass(
            self.absolute_mass_defect.value()?,
            controls.mass_conservation_tolerance_kg,
        )?;
        let residual_j = self.absolute_energy_defect.value()?;
        if residual_j > controls.energy_conservation_tolerance_j {
            return Err(GasEvaporationError::EnergyBalance {
                residual_j,
                tolerance_j: controls.energy_conservation_tolerance_j,
            });
        }
        Ok(())
    }

    pub(super) fn finish(
        mut self,
        endpoint: &GasDropletChamberCandidate<'_>,
        controls: AdaptiveGasDropletControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<GasDropletEvolutionAudit, GasEvaporationError> {
        poll(cancelled)?;
        self.gate(controls)?;
        if endpoint.droplets().len() != self.surfaces.len()
            || endpoint.gas().volume_m3().to_bits() != self.gas_volume_m3.to_bits()
        {
            return Err(invalid(
                "final droplet accounting lost its roster or fixed volume",
            ));
        }
        let tolerance = controls.mass_conservation_tolerance_kg;
        let mut closed_mass = Total::default();
        let mut final_absolute_mass = Total::default();
        let mut droplets = Vec::with_capacity(self.surfaces.len());
        for (total, row) in self.surfaces.into_iter().zip(endpoint.droplets()) {
            poll(cancelled)?;
            if total.id != row.id {
                return Err(invalid(
                    "final droplet accounting changed a stable identity",
                ));
            }
            let change = finite(row.droplet.mass_kg - total.initial_mass_kg)?;
            let transferred_mass_kg = total.transferred.value()?;
            let mass_residual_kg = finite(change + transferred_mass_kg)?;
            admit_mass(mass_residual_kg, tolerance)?;
            final_absolute_mass.add(mass_residual_kg.abs())?;
            closed_mass.add(change)?;
            self.species
                .get_mut(&total.species)
                .ok_or_else(|| invalid("final droplet lost its bound gas species"))?
                .liquid_change
                .add(change)?;
            droplets.push(GasDropletEvolutionSurfaceAudit {
                id: total.id,
                mass_kg: [total.initial_mass_kg, row.droplet.mass_kg],
                transferred_mass_kg,
                evaporated_mass_kg: total.evaporated.value()?,
                condensed_mass_kg: total.condensed.value()?,
                latent_energy_j: total.latent.value()?,
                energy_into_gas_j: total.vapor.value()?,
                sensible_heat_into_liquid_j: total.sensible.value()?,
                liquid_enthalpy_change_j: total.liquid_energy.value()?,
                mass_residual_kg,
                accumulated_mass_residual_kg: total.mass_defect.value()?,
                phase_energy_residual_j: total.phase_energy_defect.value()?,
            });
        }
        let final_species = endpoint.gas().species_inventory_kg()?;
        if final_species
            .keys()
            .any(|id| !self.species.contains_key(id))
        {
            return Err(invalid(
                "final droplet atmosphere leaves its original catalog",
            ));
        }
        let mut species = Vec::with_capacity(self.species.len());
        let mut maximum_species_mass_residual_kg = 0.0_f64;
        for (id, total) in self.species {
            poll(cancelled)?;
            let mass = final_species.get(&id).copied().unwrap_or(0.0);
            let change = finite(mass - total.initial_mass_kg)?;
            let transferred_mass_kg = total.transferred.value()?;
            let mass_residual_kg = finite(change - transferred_mass_kg)?;
            let closed_mass_residual_kg = finite(change + total.liquid_change.value()?)?;
            admit_mass(mass_residual_kg, tolerance)?;
            admit_mass(closed_mass_residual_kg, tolerance)?;
            maximum_species_mass_residual_kg = maximum_species_mass_residual_kg
                .max(mass_residual_kg.abs())
                .max(closed_mass_residual_kg.abs());
            species.push(GasDropletEvolutionSpeciesAudit {
                species: id,
                gas_mass_kg: [total.initial_mass_kg, mass],
                transferred_mass_kg,
                evaporated_mass_kg: total.evaporated.value()?,
                condensed_mass_kg: total.condensed.value()?,
                mass_residual_kg,
                closed_mass_residual_kg,
            });
        }
        let gas_change = finite(endpoint.gas().mass_kg() - self.initial_gas_mass_kg)?;
        let transferred_mass_kg = self.transferred.value()?;
        let gas_mass_residual_kg = finite(gas_change - transferred_mass_kg)?;
        admit_mass(gas_mass_residual_kg, tolerance)?;
        final_absolute_mass.add(gas_mass_residual_kg.abs())?;
        closed_mass.add(gas_change)?;
        let mass_residual_kg = closed_mass.value()?;
        admit_mass(mass_residual_kg, tolerance)?;
        let sum_absolute_mass_residual_kg = self
            .absolute_mass_defect
            .value()?
            .max(final_absolute_mass.value()?)
            .max(mass_residual_kg.abs())
            .max(maximum_species_mass_residual_kg);
        admit_mass(sum_absolute_mass_residual_kg, tolerance)?;
        poll(cancelled)?;
        Ok(GasDropletEvolutionAudit {
            droplets,
            species,
            transferred_mass_kg,
            evaporated_mass_kg: self.evaporated.value()?,
            condensed_mass_kg: self.condensed.value()?,
            latent_energy_j: self.latent.value()?,
            energy_into_gas_j: self.vapor.value()?,
            sensible_heat_into_liquid_j: self.sensible.value()?,
            liquid_enthalpy_change_j: self.liquid_energy.value()?,
            gas_internal_energy_change_j: self.gas_energy.value()?,
            mass_residual_kg,
            gas_mass_residual_kg,
            maximum_species_mass_residual_kg,
            accumulated_mass_residual_kg: self.mass_defect.value()?,
            sum_absolute_mass_residual_kg,
            energy_residual_j: self.energy_defect.value()?,
            sum_absolute_energy_residual_j: self.absolute_energy_defect.value()?,
        })
    }
}

fn invalid(message: &'static str) -> GasEvaporationError {
    GasEvaporationError::Invalid(message)
}
fn finite(value: f64) -> Result<f64, GasEvaporationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "adaptive droplet inventory account is unrepresentable",
        ))
    }
}
fn admit_mass(residual_kg: f64, tolerance_kg: f64) -> Result<(), GasEvaporationError> {
    if !residual_kg.is_finite() || residual_kg.abs() > tolerance_kg {
        Err(invalid(
            "adaptive droplets exceed their whole-window mass allowance",
        ))
    } else {
        Ok(())
    }
}
