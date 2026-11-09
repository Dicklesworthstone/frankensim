//! Accepted interval accounts for one droplet cloud and its complete gas network.
//! No trial history or gas solver is retained. Original caloric residuals close
//! energy; independent original-to-final inventories close mass and species.

use std::collections::{BTreeMap, btree_map::Entry};

use super::{GasDropletNetworkLimits, admit_work, invalid};
use super::super::{
    AdaptiveGasDropletControls, ClosedGasParcel, GasDropletEvolutionAudit,
    GasDropletEvolutionSpeciesAudit, GasDropletEvolutionSurfaceAudit, GasDropletSeed,
    GasDropletSurfaceStep, GasEvaporationError, poll,
};
use crate::evaporation::gas::{
    chamber::Total,
    droplet::chamber::GasDropletChamberReport,
};
use crate::thermochemical_gas::energy::parcel::transfer::species::network::transport::GasTransportStep;
use fs_thermochem::SpeciesId;

/// Original-to-final accounting for one species in one finite gas parcel.
#[derive(Clone, Debug, PartialEq)]
pub struct GasDropletNetworkParcelSpeciesAudit {
    pub species: SpeciesId,
    pub mass_kg: [f64; 2],
    /// Signed accepted network transport into this parcel.
    pub transported_mass_kg: f64,
    /// Signed evaporation into this parcel, excluding gas/gas transport.
    pub phase_transferred_mass_kg: f64,
    pub mass_residual_kg: f64,
}

/// Every original gas parcel participates, including idle and droplet-free ones.
#[derive(Clone, Debug, PartialEq)]
pub struct GasDropletNetworkParcelAudit {
    pub id: u64,
    pub mass_kg: [f64; 2],
    pub net_transported_mass_kg: f64,
    pub phase_transferred_mass_kg: f64,
    /// Sum of original source-relative storage reports. Their independently
    /// bounded inventory rounding is not multiplied by a formation reference.
    pub gas_internal_energy_change_j: f64,
    pub mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub species: Vec<GasDropletNetworkParcelSpeciesAudit>,
}

/// Physical transfers cover only accepted original intervals in their order.
#[derive(Clone, Debug, PartialEq)]
pub struct GasDropletNetworkAudit {
    /// Compensated sum of the actual accepted original interval durations.
    pub duration_s: f64,
    pub accepted_intervals: usize,
    /// Existing liquid/phase accounts; gas totals include the COMPLETE network.
    /// Original transport and phase caloric defects share one window allowance.
    pub exchange: GasDropletEvolutionAudit,
    pub parcels: Vec<GasDropletNetworkParcelAudit>,
}

#[derive(Default)]
struct TransferTotals {
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
struct SurfaceTotals {
    id: u64,
    species: SpeciesId,
    initial_mass: f64,
    last_droplet: super::super::SphericalDropletState,
    totals: TransferTotals,
}
#[derive(Default)]
struct LocalSpecies {
    initial_mass: f64,
    transported: Total,
    phase: Total,
}
struct ParcelTotals {
    initial_mass: f64,
    volume: f64,
    transported: Total,
    phase: Total,
    energy: Total,
    species: BTreeMap<SpeciesId, LocalSpecies>,
}
#[derive(Default)]
struct SpeciesTotals {
    initial_mass: Total,
    final_mass: Total,
    gas_change: Total,
    transferred: Total,
    evaporated: Total,
    condensed: Total,
    liquid_change: Total,
}

/// Online physical account for consecutive original transport/chamber intervals.
///
/// Call record only for accepted intervals whose input is the preceding accepted
/// endpoint. The original interval owners admit sources, rates and calorics;
/// this account does not independently solve or certify their physical laws.
/// Record consumes the account so a refused partial update cannot be reused.
/// Storage follows the droplet roster and actual per-parcel species; there is
/// no list of past trajectories or source samples.
///
/// The gas network is closed except for the retained liquid phase. External
/// gas streams or heat need their own opposite finite storage/boundary owner
/// and are refused here. Internal molecular energy corrections remain original.
pub struct GasDropletNetworkWindow {
    gas_region: u64,
    controls: AdaptiveGasDropletControls,
    limits: GasDropletNetworkLimits,
    surfaces: Vec<SurfaceTotals>,
    parcels: BTreeMap<u64, ParcelTotals>,
    species: BTreeMap<SpeciesId, SpeciesTotals>,
    species_entries: usize,
    intervals: usize,
    duration: Total,
    phase: TransferTotals,
    gas_energy: Total,
    mass_defect: Total,
    absolute_mass_defect: Total,
    energy_defect: Total,
    absolute_energy_defect: Total,
}

impl GasDropletNetworkWindow {
    /// Read existing inventory metadata only, never a liquid or film source.
    pub fn new(
        initial: &[GasDropletSeed<'_>],
        gas_region: u64,
        before: &BTreeMap<u64, ClosedGasParcel>,
        controls: AdaptiveGasDropletControls,
        limits: GasDropletNetworkLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, GasEvaporationError> {
        poll(cancelled)?;
        controls.validate(0.0, initial.len())?;
        admit_work(before, limits, cancelled)?;
        if !before.contains_key(&gas_region)
            || initial.windows(2).any(|rows| rows[0].id >= rows[1].id)
        {
            return Err(invalid("network window needs original region and canonical droplet IDs"));
        }
        let mut species = BTreeMap::<SpeciesId, SpeciesTotals>::new();
        let mut parcels = BTreeMap::new();
        let mut species_entries = 0usize;
        for (&id, gas) in before {
            poll(cancelled)?;
            let inventory = gas.species_inventory_kg()?;
            species_entries = species_entries
                .checked_add(inventory.len())
                .ok_or_else(|| invalid("network window species count overflow"))?;
            if species_entries > limits.maximum_species_entries {
                return Err(invalid("network window exceeds its species work ceiling"));
            }
            let mut local = BTreeMap::new();
            for (id, mass) in inventory {
                poll(cancelled)?;
                species.entry(id.clone()).or_default().initial_mass.add(mass)?;
                local.insert(id, LocalSpecies {
                    initial_mass: mass,
                    ..LocalSpecies::default()
                });
            }
            parcels.insert(id, ParcelTotals {
                initial_mass: gas.mass_kg(),
                volume: gas.volume_m3(),
                transported: Total::default(),
                phase: Total::default(),
                energy: Total::default(),
                species: local,
            });
        }
        let mut surfaces = Vec::with_capacity(initial.len());
        for seed in initial {
            poll(cancelled)?;
            if !seed.droplet.mass_kg.is_finite() || seed.droplet.mass_kg <= 0.0 {
                return Err(invalid("network window needs positive original liquid inventories"));
            }
            species.entry(seed.binding.species.clone()).or_default();
            surfaces.push(SurfaceTotals {
                id: seed.id,
                species: seed.binding.species.clone(),
                initial_mass: seed.droplet.mass_kg,
                last_droplet: seed.droplet,
                totals: TransferTotals::default(),
            });
        }
        Ok(Self {
            gas_region,
            controls,
            limits,
            surfaces,
            parcels,
            species,
            species_entries,
            intervals: 0,
            duration: Total::default(),
            phase: TransferTotals::default(),
            gas_energy: Total::default(),
            mass_defect: Total::default(),
            absolute_mass_defect: Total::default(),
            energy_defect: Total::default(),
            absolute_energy_defect: Total::default(),
        })
    }

    /// Add one accepted original transport-then-chamber interval.
    /// Coarse/rejected trials must never be passed to this method.
    pub fn record(
        mut self,
        transport: &GasTransportStep,
        coupled: &GasDropletChamberReport,
        droplets: &[GasDropletSurfaceStep],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, GasEvaporationError> {
        poll(cancelled)?;
        let part = &coupled.exchange;
        if !transport.duration_s.is_finite()
            || transport.duration_s <= 0.0
            || transport.duration_s.to_bits() != part.duration_s.to_bits()
            || transport.duration_s.to_bits() != coupled.phase.duration_s.to_bits()
            || droplets.len() != self.surfaces.len()
            || !self.parcels.keys().copied().eq(
                transport.exchange.parcels.iter().map(|row| row.id),
            )
            || !transport.exchange.boundaries.is_empty()
            || !transport.exchange.external_heat.is_empty()
        {
            return Err(invalid("network window needs one closed original transport/phase interval"));
        }
        let mut local_mass_defect = Total::default();
        let mut local_energy_defect = Total::default();
        for row in &transport.exchange.parcels {
            poll(cancelled)?;
            local_mass_defect.add(nonnegative(row.mass_residual_kg.abs())?
                .max(nonnegative(row.maximum_species_mass_residual_kg)?))?;
            local_energy_defect.add(row.energy_residual_j.abs())?;
            let parcel = self.parcels.get_mut(&row.id)
                .ok_or_else(|| invalid("accepted gas transport changed parcel identity"))?;
            parcel.energy.add(row.energy_change_j)?;
            self.gas_energy.add(row.energy_change_j)?;
        }
        for flux in &transport.exchange.fluxes {
            poll(cancelled)?;
            let flow = &flux.transfer;
            nonnegative(flow.mass_kg)?;
            if !self.species.contains_key(&flow.species) {
                return Err(invalid("accepted transport introduced an unknown source species"));
            }
            for (id, mass) in [(flow.donor, -flow.mass_kg), (flow.receiver, flow.mass_kg)] {
                let parcel = self.parcels.get_mut(&id)
                    .ok_or_else(|| invalid("accepted transport left its original gas roster"))?;
                parcel.transported.add(mass)?;
                local_species(
                    parcel, &flow.species, &mut self.species_entries,
                    self.limits.maximum_species_entries,
                )?.transported.add(mass)?;
            }
        }
        for (total, row) in self.surfaces.iter_mut().zip(droplets) {
            poll(cancelled)?;
            let geometry = super::super::SphericalDropletGeometry::from_mass(
                total.last_droplet.mass_kg,
                row.geometry[0].liquid_density_kg_m3(),
            )?;
            if total.id != row.id || geometry != row.geometry[0] {
                return Err(invalid("accepted droplet interval is not contiguous in identity or mass"));
            }
            total.last_droplet = row.droplet;
            let phase = &row.exchange;
            total.totals.transferred.add(phase.transferred_mass_kg)?;
            total.totals.evaporated.add(phase.transferred_mass_kg.max(0.0))?;
            total.totals.condensed.add((-phase.transferred_mass_kg).max(0.0))?;
            total.totals.latent.add(phase.latent_energy_j)?;
            total.totals.vapor.add(phase.energy_into_gas_j)?;
            total.totals.sensible.add(row.heat.map_or(0.0, |heat| heat.heat_into_liquid_j))?;
            total.totals.liquid_energy.add(phase.liquid_enthalpy_change_j)?;
            total.totals.mass_defect.add(phase.mass_residual_kg)?;
            total.totals.phase_energy_defect.add(row.phase.energy_residual_j)?;
        }
        let gas = self.parcels.get_mut(&self.gas_region)
            .ok_or_else(|| invalid("accepted phase lost its named gas region"))?;
        gas.phase.add(part.transferred_mass_kg)?;
        gas.energy.add(part.gas_internal_energy_change_j)?;
        self.gas_energy.add(part.gas_internal_energy_change_j)?;
        for row in &part.species {
            poll(cancelled)?;
            let total = self.species.get_mut(&row.species)
                .ok_or_else(|| invalid("accepted phase introduced an unknown bound species"))?;
            total.transferred.add(row.transferred_mass_kg)?;
            total.evaporated.add(row.evaporated_mass_kg)?;
            total.condensed.add(row.condensed_mass_kg)?;
            local_species(
                gas, &row.species, &mut self.species_entries,
                self.limits.maximum_species_entries,
            )?.phase.add(row.transferred_mass_kg)?;
        }
        self.phase.transferred.add(part.transferred_mass_kg)?;
        self.phase.evaporated.add(part.evaporated_mass_kg)?;
        self.phase.condensed.add(part.condensed_mass_kg)?;
        self.phase.latent.add(part.latent_energy_j)?;
        self.phase.vapor.add(part.energy_into_gas_j)?;
        self.phase.sensible.add(coupled.heat.map_or(0.0, |heat| heat.heat_into_liquid_j))?;
        self.phase.liquid_energy.add(part.liquid_enthalpy_change_j)?;
        self.mass_defect.add(transport.exchange.total_mass_residual_kg)?;
        self.mass_defect.add(part.mass_residual_kg)?;
        self.absolute_mass_defect.add(local_mass_defect.value()?
            .max(nonnegative(transport.exchange.total_mass_residual_kg.abs())?)
            .max(nonnegative(transport.exchange.maximum_species_mass_residual_kg)?))?;
        self.absolute_mass_defect.add(nonnegative(coupled.sum_absolute_mass_residual_kg)?)?;
        self.energy_defect.add(transport.exchange.total_energy_residual_j)?;
        self.energy_defect.add(part.energy_residual_j)?;
        self.absolute_energy_defect.add(local_energy_defect.value()?
            .max(nonnegative(transport.exchange.total_energy_residual_j.abs())?))?;
        self.absolute_energy_defect.add(nonnegative(part.sum_absolute_energy_residual_j)?)?;
        self.duration.add(transport.duration_s)?;
        self.intervals = self.intervals.checked_add(1)
            .ok_or_else(|| invalid("accepted network interval count overflow"))?;
        self.gate()?;
        poll(cancelled)?;
        Ok(self)
    }

    fn gate(&self) -> Result<(), GasEvaporationError> {
        admit_mass(
            self.absolute_mass_defect.value()?.max(self.mass_defect.value()?.abs()),
            self.controls.mass_conservation_tolerance_kg,
        )?;
        let residual_j = self.absolute_energy_defect.value()?
            .max(self.energy_defect.value()?.abs());
        if residual_j > self.controls.energy_conservation_tolerance_j {
            return Err(GasEvaporationError::EnergyBalance {
                residual_j,
                tolerance_j: self.controls.energy_conservation_tolerance_j,
            });
        }
        Ok(())
    }

    /// Finish against actual final inventories, independently of rounded reports.
    pub fn finish(
        mut self,
        after: &BTreeMap<u64, ClosedGasParcel>,
        droplets: &[GasDropletSurfaceStep],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<GasDropletNetworkAudit, GasEvaporationError> {
        poll(cancelled)?;
        self.gate()?;
        admit_work(after, self.limits, cancelled)?;
        if self.intervals == 0
            || !self.parcels.keys().eq(after.keys())
            || droplets.len() != self.surfaces.len()
        {
            return Err(invalid("network window requires its accepted final physical roster"));
        }
        let tolerance = self.controls.mass_conservation_tolerance_kg;
        let mut closed_mass = Total::default();
        let mut gas_change = Total::default();
        let mut final_absolute_mass = Total::default();
        let mut liquid_audits = Vec::with_capacity(self.surfaces.len());
        for (total, row) in self.surfaces.into_iter().zip(droplets) {
            poll(cancelled)?;
            if total.id != row.id || total.last_droplet != row.droplet {
                return Err(invalid("network window final liquid differs from its accepted endpoint"));
            }
            let change = finite(row.droplet.mass_kg - total.initial_mass)?;
            let transferred = total.totals.transferred.value()?;
            let residual = difference(change, -transferred, 0.0)?;
            admit_mass(residual, tolerance)?;
            final_absolute_mass.add(residual.abs())?;
            closed_mass.add(change)?;
            self.species.get_mut(&total.species)
                .ok_or_else(|| invalid("network window lost a liquid species binding"))?
                .liquid_change.add(change)?;
            liquid_audits.push(GasDropletEvolutionSurfaceAudit {
                id: total.id,
                mass_kg: [total.initial_mass, row.droplet.mass_kg],
                transferred_mass_kg: transferred,
                evaporated_mass_kg: total.totals.evaporated.value()?,
                condensed_mass_kg: total.totals.condensed.value()?,
                latent_energy_j: total.totals.latent.value()?,
                energy_into_gas_j: total.totals.vapor.value()?,
                sensible_heat_into_liquid_j: total.totals.sensible.value()?,
                liquid_enthalpy_change_j: total.totals.liquid_energy.value()?,
                mass_residual_kg: residual,
                accumulated_mass_residual_kg: total.totals.mass_defect.value()?,
                phase_energy_residual_j: total.totals.phase_energy_defect.value()?,
            });
        }
        let mut parcel_audits = Vec::with_capacity(self.parcels.len());
        let mut maximum_species_residual = 0.0_f64;
        for (id, mut total) in self.parcels {
            poll(cancelled)?;
            let gas = &after[&id];
            if gas.volume_m3().to_bits() != total.volume.to_bits() {
                return Err(invalid("network window changed an original fixed gas volume"));
            }
            let inventory = gas.species_inventory_kg()?;
            for species in inventory.keys() {
                poll(cancelled)?;
                if !self.species.contains_key(species) {
                    return Err(invalid("network window final gas has an unknown source species"));
                }
                local_species(
                    &mut total, species, &mut self.species_entries,
                    self.limits.maximum_species_entries,
                )?;
            }
            let mut local_audits = Vec::with_capacity(total.species.len());
            let mut local_maximum = 0.0_f64;
            for (species, local) in total.species {
                poll(cancelled)?;
                let mass = inventory.get(&species).copied().unwrap_or(0.0);
                let change = finite(mass - local.initial_mass)?;
                let transported = local.transported.value()?;
                let phase = local.phase.value()?;
                let residual = difference(change, transported, phase)?;
                admit_mass(residual, tolerance)?;
                local_maximum = local_maximum.max(residual.abs());
                let global = self.species.get_mut(&species)
                    .ok_or_else(|| invalid("network window final species lost its source"))?;
                global.final_mass.add(mass)?;
                global.gas_change.add(change)?;
                local_audits.push(GasDropletNetworkParcelSpeciesAudit {
                    species,
                    mass_kg: [local.initial_mass, mass],
                    transported_mass_kg: transported,
                    phase_transferred_mass_kg: phase,
                    mass_residual_kg: residual,
                });
            }
            let change = finite(gas.mass_kg() - total.initial_mass)?;
            gas_change.add(change)?;
            closed_mass.add(change)?;
            let transported = total.transported.value()?;
            let phase = total.phase.value()?;
            let residual = difference(change, transported, phase)?;
            admit_mass(residual, tolerance)?;
            final_absolute_mass.add(residual.abs().max(local_maximum))?;
            maximum_species_residual = maximum_species_residual.max(local_maximum);
            parcel_audits.push(GasDropletNetworkParcelAudit {
                id,
                mass_kg: [total.initial_mass, gas.mass_kg()],
                net_transported_mass_kg: transported,
                phase_transferred_mass_kg: phase,
                gas_internal_energy_change_j: total.energy.value()?,
                mass_residual_kg: residual,
                maximum_species_mass_residual_kg: local_maximum,
                species: local_audits,
            });
        }
        let mut species_audits = Vec::with_capacity(self.species.len());
        for (species, total) in self.species {
            poll(cancelled)?;
            let change = total.gas_change.value()?;
            let transferred = total.transferred.value()?;
            let residual = difference(change, transferred, 0.0)?;
            let closed_residual = difference(change, -total.liquid_change.value()?, 0.0)?;
            admit_mass(residual, tolerance)?;
            admit_mass(closed_residual, tolerance)?;
            maximum_species_residual = maximum_species_residual
                .max(residual.abs()).max(closed_residual.abs());
            species_audits.push(GasDropletEvolutionSpeciesAudit {
                species,
                gas_mass_kg: [total.initial_mass.value()?, total.final_mass.value()?],
                transferred_mass_kg: transferred,
                evaporated_mass_kg: total.evaporated.value()?,
                condensed_mass_kg: total.condensed.value()?,
                mass_residual_kg: residual,
                closed_mass_residual_kg: closed_residual,
            });
        }
        let transferred = self.phase.transferred.value()?;
        let gas_residual = difference(gas_change.value()?, transferred, 0.0)?;
        let mass_residual = closed_mass.value()?;
        let absolute_mass = self.absolute_mass_defect.value()?
            .max(final_absolute_mass.value()?)
            .max(gas_residual.abs())
            .max(mass_residual.abs())
            .max(maximum_species_residual);
        admit_mass(absolute_mass, tolerance)?;
        let energy_residual = self.energy_defect.value()?;
        poll(cancelled)?;
        Ok(GasDropletNetworkAudit {
            duration_s: self.duration.value()?,
            accepted_intervals: self.intervals,
            exchange: GasDropletEvolutionAudit {
                droplets: liquid_audits,
                species: species_audits,
                transferred_mass_kg: transferred,
                evaporated_mass_kg: self.phase.evaporated.value()?,
                condensed_mass_kg: self.phase.condensed.value()?,
                latent_energy_j: self.phase.latent.value()?,
                energy_into_gas_j: self.phase.vapor.value()?,
                sensible_heat_into_liquid_j: self.phase.sensible.value()?,
                liquid_enthalpy_change_j: self.phase.liquid_energy.value()?,
                gas_internal_energy_change_j: self.gas_energy.value()?,
                mass_residual_kg: mass_residual,
                gas_mass_residual_kg: gas_residual,
                maximum_species_mass_residual_kg: maximum_species_residual,
                accumulated_mass_residual_kg: self.mass_defect.value()?,
                sum_absolute_mass_residual_kg: absolute_mass,
                energy_residual_j: energy_residual,
                sum_absolute_energy_residual_j: self.absolute_energy_defect.value()?
                    .max(energy_residual.abs()),
            },
            parcels: parcel_audits,
        })
    }
}

fn local_species<'a>(
    parcel: &'a mut ParcelTotals,
    species: &SpeciesId,
    entries: &mut usize,
    maximum: usize,
) -> Result<&'a mut LocalSpecies, GasEvaporationError> {
    match parcel.species.entry(species.clone()) {
        Entry::Occupied(entry) => Ok(entry.into_mut()),
        Entry::Vacant(entry) => {
            let next = entries.checked_add(1)
                .ok_or_else(|| invalid("network window species count overflow"))?;
            if next > maximum {
                return Err(invalid("network window exceeds its retained species work ceiling"));
            }
            *entries = next;
            Ok(entry.insert(LocalSpecies::default()))
        }
    }
}

fn finite(value: f64) -> Result<f64, GasEvaporationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("network window arithmetic is unrepresentable"))
    }
}
fn nonnegative(value: f64) -> Result<f64, GasEvaporationError> {
    let value = finite(value)?;
    if value < 0.0 {
        Err(invalid("network window requires nonnegative absolute defects and gross mass"))
    } else {
        Ok(value)
    }
}
fn difference(value: f64, first: f64, second: f64) -> Result<f64, GasEvaporationError> {
    let mut total = Total::default();
    total.add(value)?;
    total.add(-first)?;
    total.add(-second)?;
    total.value()
}
fn admit_mass(residual: f64, tolerance: f64) -> Result<(), GasEvaporationError> {
    if !residual.is_finite() || residual.abs() > tolerance {
        Err(invalid("droplet network exceeded its whole-window mass allowance"))
    } else {
        Ok(())
    }
}
