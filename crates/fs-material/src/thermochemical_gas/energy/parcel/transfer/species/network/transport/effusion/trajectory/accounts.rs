//! Accepted original streams and independent whole-window inventory balances.

use super::{
    ClosedGasParcel, EffusionAdvanceControls, Error, GasTransportSolve, MAX_EFFUSION_RETAINED_ROWS,
    ReservoirEffusionCandidate, ReservoirEffusionStep, SpeciesId, inventory, poll, ratio,
};
use fs_math::dd::Dd;
use fs_qty::Temperature;
use fs_thermochem::{Nasa9StandardStateModelV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use std::collections::BTreeMap;

fn add(total: &mut Dd, value: f64) -> Result<(), Error> {
    add_dd(total, Dd::from_f64(value))
}
fn add_dd(total: &mut Dd, value: Dd) -> Result<(), Error> {
    *total = *total + value;
    if !value.hi.is_finite()
        || !value.lo.is_finite()
        || !total.hi.is_finite()
        || !total.lo.is_finite()
    {
        return Err(Error::Invalid(
            "unrepresentable integrated effusion account",
        ));
    }
    Ok(())
}

#[derive(Clone, Debug, PartialEq)]
pub struct EffusionIntegratedBoundary {
    pub id: u64,
    pub parcel: u64,
    pub species: SpeciesId,
    pub incoming_mass_kg: f64,
    pub outgoing_mass_kg: f64,
    /// Independently summed signed stream, including resolved tiny drift.
    pub net_mass_into_gas_kg: f64,
    /// Complete original energy for reporting. An arbitrary formation datum
    /// can dominate this number; it never scales a conservation/error budget.
    pub energy_into_gas_j: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffusionWindowBalance {
    pub total_mass_residual_kg: f64,
    pub maximum_parcel_species_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    /// Actual endpoint internal-energy change minus accepted boundary energy,
    /// using a fixed internal-energy reference for each exact original source.
    /// Species representation rounding therefore enters this independent audit.
    pub energy_residual_j: f64,
    /// Compensated sum of accepted ORIGINAL caloric closure residuals; retained
    /// separately because requested and represented species increments can differ.
    pub accumulated_energy_residual_j: f64,
}

struct CaloricReference {
    source: Nasa9StandardStateModelV1,
    temperature_k: f64,
}
// Disconnected parcels may use different source records for the same name.
// Only exact original source equality shares an energy reference.
type CaloricReferences = BTreeMap<SpeciesId, Vec<CaloricReference>>;

fn source<'a>(
    gas: &'a ClosedGasParcel,
    species: &SpeciesId,
) -> Result<&'a Nasa9StandardStateModelV1, Error> {
    let components = gas.model().thermodynamics().components();
    components
        .binary_search_by(|s| s.species().cmp(species))
        .map(|i| &components[i])
        .map_err(|_| Error::Invalid("energy audit lost its original species source"))
}
fn reference(
    references: &CaloricReferences,
    source: &Nasa9StandardStateModelV1,
) -> Result<f64, Error> {
    references
        .get(source.species())
        .and_then(|rows| {
            rows.iter()
                .find(|r| r.source == *source)
                .map(|r| r.temperature_k)
        })
        .ok_or(Error::Invalid(
            "energy audit lost its exact original caloric reference",
        ))
}
fn internal_change(source: &Nasa9StandardStateModelV1, from: f64, to: f64) -> Result<Dd, Error> {
    let molar = source
        .molar_internal_energy_change(Temperature::new(from), Temperature::new(to))
        .map_err(|_| Error::Invalid("energy audit left its original caloric support"))?
        .value();
    Ok(Dd::from_f64(molar) / Dd::from_f64(source.molar_mass().value()))
}

#[derive(Default)]
pub(super) struct Accepted {
    inventory_change: BTreeMap<(u64, SpeciesId), Dd>,
    boundaries: BTreeMap<(u64, u64, SpeciesId), [Dd; 4]>,
    energy_residual: Dd,
    calorics: CaloricReferences,
    source_count: usize,
    reduced_boundary_energy: Dd,
}
impl Accepted {
    pub(super) fn record(
        &mut self,
        before: &BTreeMap<u64, ClosedGasParcel>,
        result: &ReservoirEffusionCandidate,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        // Register old sources before newly arriving ones, at an actual supported
        // temperature. Reference temperatures stay fixed for the whole window.
        for gas in before.values().chain(result.parcels().values()) {
            for source in gas.model().thermodynamics().components() {
                poll(cancelled)?;
                let rows = self.calorics.entry(source.species().clone()).or_default();
                if !rows.iter().any(|r| r.source == *source) {
                    if self.source_count == MAX_EFFUSION_RETAINED_ROWS {
                        return Err(Error::Invalid(
                            "energy source table exceeds retained work bound",
                        ));
                    }
                    rows.push(CaloricReference {
                        source: source.clone(),
                        temperature_k: gas.state().state().temperature,
                    });
                    self.source_count += 1;
                }
            }
        }
        let step = result.step();
        for flux in &step.exchange.fluxes {
            poll(cancelled)?;
            let row = &flux.transfer;
            for (parcel, mass) in [(row.donor, -row.mass_kg), (row.receiver, row.mass_kg)] {
                add(
                    self.inventory_change
                        .entry((parcel, row.species.clone()))
                        .or_default(),
                    mass,
                )?;
            }
        }
        for flux in &step.exchange.boundaries {
            poll(cancelled)?;
            let row = &flux.transfer;
            let original = source(&result.parcels()[&flux.parcel], &row.species)?;
            let t0 = before[&flux.parcel].state().state().temperature;
            let tref = reference(&self.calorics, original)?;
            // E - dm*u(Tref) = relative_E + dm*(Delta u(Tref,T0)+R_i*T0).
            // The source owns Delta u; no absolute formation energies subtract.
            let specific = internal_change(original, tref, t0)?
                + Dd::from_f64(R) * Dd::from_f64(t0) / Dd::from_f64(original.molar_mass().value());
            add_dd(
                &mut self.reduced_boundary_energy,
                Dd::from_f64(row.relative_energy_into_gas_j)
                    + Dd::from_f64(row.net_mass_into_gas_kg) * specific,
            )?;
            add(
                self.inventory_change
                    .entry((flux.parcel, row.species.clone()))
                    .or_default(),
                row.net_mass_into_gas_kg,
            )?;
            let totals = self
                .boundaries
                .entry((flux.id, flux.parcel, row.species.clone()))
                .or_default();
            for (sum, value) in totals.iter_mut().zip([
                row.incoming_mass_kg,
                row.outgoing_mass_kg,
                row.net_mass_into_gas_kg,
                flux.energy_into_gas_j,
            ]) {
                add(sum, value)?;
            }
        }
        add(
            &mut self.energy_residual,
            step.exchange.total_energy_residual_j,
        )
    }

    pub(super) fn finish(
        self,
        before: &BTreeMap<u64, ClosedGasParcel>,
        after: &BTreeMap<u64, ClosedGasParcel>,
        controls: GasTransportSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<EffusionIntegratedBoundary>, EffusionWindowBalance), Error> {
        let mut residuals: BTreeMap<_, _> = self
            .inventory_change
            .into_iter()
            .map(|(key, value)| (key, Dd::ZERO - value))
            .collect();
        let mut mass = Dd::ZERO;
        let mut energy = Dd::ZERO - self.reduced_boundary_energy;
        for (&id, old) in before {
            poll(cancelled)?;
            let next = &after[&id];
            add(&mut mass, next.mass_kg())?;
            add(&mut mass, -old.mass_kg())?;
            let (old_inventory, next_inventory) = (inventory(old)?, inventory(next)?);
            for (inventory, sign) in [(&old_inventory, -1.0), (&next_inventory, 1.0)] {
                for (species, &value) in inventory {
                    poll(cancelled)?;
                    add(
                        residuals.entry((id, species.clone())).or_default(),
                        sign * value,
                    )?;
                }
            }
            for species in old_inventory.keys().chain(
                next_inventory
                    .keys()
                    .filter(|s| !old_inventory.contains_key(*s)),
            ) {
                poll(cancelled)?;
                let m0 = old_inventory.get(species).copied().unwrap_or(0.0);
                let m1 = next_inventory.get(species).copied().unwrap_or(0.0);
                let original = source(if m1 > 0.0 { next } else { old }, species)?;
                let tref = reference(&self.calorics, original)?;
                let (t0, t1) = (
                    old.state().state().temperature,
                    next.state().state().temperature,
                );
                let change = if m0 == 0.0 {
                    Dd::from_f64(m1) * internal_change(original, tref, t1)?
                } else if m1 == 0.0 {
                    Dd::from_f64(-m0) * internal_change(original, tref, t0)?
                } else {
                    if source(old, species)? != original {
                        return Err(Error::Invalid(
                            "one parcel changed its original caloric source",
                        ));
                    }
                    // Stable endpoint difference, preserving a small thermal
                    // change beside a large fixed reference-temperature offset.
                    Dd::from_f64(m1) * internal_change(original, t0, t1)?
                        + (Dd::from_f64(m1) - Dd::from_f64(m0))
                            * internal_change(original, tref, t0)?
                };
                add_dd(&mut energy, change)?;
            }
        }
        let mut boundaries = Vec::with_capacity(self.boundaries.len());
        for ((id, parcel, species), values) in self.boundaries {
            poll(cancelled)?;
            mass = mass - values[2];
            boundaries.push(EffusionIntegratedBoundary {
                id,
                parcel,
                species,
                incoming_mass_kg: values[0].to_f64(),
                outgoing_mass_kg: values[1].to_f64(),
                net_mass_into_gas_kg: values[2].to_f64(),
                energy_into_gas_j: values[3].to_f64(),
            });
        }
        let mut maximum_local = 0.0_f64;
        let mut species_residuals: BTreeMap<SpeciesId, Dd> = BTreeMap::new();
        for ((_, species), value) in residuals {
            poll(cancelled)?;
            maximum_local = maximum_local.max(value.to_f64().abs());
            let sum = species_residuals.entry(species).or_default();
            *sum = *sum + value;
        }
        let maximum_species = species_residuals
            .values()
            .map(|v| v.to_f64().abs())
            .fold(0.0_f64, f64::max);
        let balance = EffusionWindowBalance {
            total_mass_residual_kg: mass.to_f64(),
            maximum_parcel_species_mass_residual_kg: maximum_local,
            maximum_species_mass_residual_kg: maximum_species,
            energy_residual_j: energy.to_f64(),
            accumulated_energy_residual_j: self.energy_residual.to_f64(),
        };
        if !balance.total_mass_residual_kg.is_finite()
            || balance.total_mass_residual_kg.abs() > controls.exchange.mass_tolerance_kg
            || maximum_local > controls.exchange.mass_tolerance_kg
            || maximum_species > controls.exchange.mass_tolerance_kg
            || !balance.energy_residual_j.is_finite()
            || balance.energy_residual_j.abs() > controls.exchange.energy_tolerance_j
            || !balance.accumulated_energy_residual_j.is_finite()
            || balance.accumulated_energy_residual_j.abs() > controls.exchange.energy_tolerance_j
        {
            return Err(Error::Invalid(
                "complete effusion window exceeded its fixed species/mass/energy budget",
            ));
        }
        Ok((boundaries, balance))
    }
}

// All original reports are retained, so count inactive faces too. The hard
// row bound prevents many large snapshots from bypassing the interval budget.
pub(super) fn row_count(step: &ReservoirEffusionStep) -> usize {
    let rates = &step.transport;
    1 + step.exchange.fluxes.len()
        + step.exchange.heat.len()
        + step.exchange.external_heat.len()
        + step.exchange.boundaries.len()
        + step.exchange.parcels.len()
        + rates.nozzles.len()
        + rates.diffusion.len()
        + rates.porous.len()
        + rates.membranes.len()
        + rates.dusty_gas.len()
        + rates
            .effusion
            .iter()
            .map(|r| 1 + r.flux.species.len())
            .sum::<usize>()
        + step
            .vacuum_outlets
            .iter()
            .map(|r| 1 + r.flux.species.len())
            .sum::<usize>()
        + step
            .reservoirs
            .iter()
            .map(|r| 1 + r.flux.species.len())
            .sum::<usize>()
        + step
            .drifting_reservoirs
            .iter()
            .map(|r| 1 + r.incoming.species.len() + r.outgoing.species.len())
            .sum::<usize>()
}

type StreamKey = (u8, u64, u64, SpeciesId);
#[derive(Default)]
struct Streams {
    /// Internal keys preserve donor/receiver direction. Boundary keys preserve
    /// independent gross in/out and signed mass plus a common caloric chart.
    species: BTreeMap<StreamKey, [Dd; 4]>,
    heat: BTreeMap<[u64; 2], [Dd; 2]>,
}
impl Streams {
    fn record(
        &mut self,
        reference: &BTreeMap<u64, ClosedGasParcel>,
        original: &BTreeMap<u64, ClosedGasParcel>,
        result: &ReservoirEffusionCandidate,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        for flux in &result.step().exchange.fluxes {
            poll(cancelled)?;
            let row = &flux.transfer;
            let sum = &mut self
                .species
                .entry((0, row.donor, row.receiver, row.species.clone()))
                .or_default()[0];
            add(sum, row.mass_kg)?;
        }
        for flux in &result.step().exchange.boundaries {
            poll(cancelled)?;
            let row = &flux.transfer;
            let values = self
                .species
                .entry((1, flux.id, flux.parcel, row.species.clone()))
                .or_default();
            for (sum, mass) in values.iter_mut().zip([
                row.incoming_mass_kg,
                row.outgoing_mass_kg,
                row.net_mass_into_gas_kg,
            ]) {
                add(sum, mass)?;
            }
            add(&mut values[3], row.relative_energy_into_gas_j)?;
            let from = reference[&flux.parcel].state().state().temperature;
            let to = original[&flux.parcel].state().state().temperature;
            if from != to && row.net_mass_into_gas_kg != 0.0 {
                // A second half-step's reference is different. Move it to the
                // trial's ORIGINAL partial enthalpy using source-owned changes;
                // never subtract two rounded absolute formation enthalpies.
                let components = result.parcels()[&flux.parcel]
                    .model()
                    .thermodynamics()
                    .components();
                let index = components
                    .binary_search_by(|s| s.species().cmp(&row.species))
                    .map_err(|_| Error::Invalid("molecular boundary lost its original species"))?;
                let source = &components[index];
                let change = source
                    .molar_enthalpy_change(Temperature::new(from), Temperature::new(to))
                    .map_err(|_| {
                        Error::Invalid(
                            "boundary energy comparison left its original caloric support",
                        )
                    })?
                    .value()
                    / source.molar_mass().value();
                add(&mut values[3], row.net_mass_into_gas_kg * change)?;
            }
        }
        for row in &result.step().exchange.heat {
            poll(cancelled)?;
            add(
                &mut self.heat.entry(row.parcels).or_default()
                    [usize::from(row.heat_left_to_right_j < 0.0)],
                row.heat_left_to_right_j.abs(),
            )?;
        }
        Ok(())
    }
}

pub(super) fn stream_error(
    before: &BTreeMap<u64, ClosedGasParcel>,
    coarse: &ReservoirEffusionCandidate,
    first: &ReservoirEffusionCandidate,
    fine: &ReservoirEffusionCandidate,
    controls: EffusionAdvanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let (mut a, mut b) = (Streams::default(), Streams::default());
    a.record(before, before, coarse, cancelled)?;
    b.record(before, before, first, cancelled)?;
    b.record(before, first.parcels(), fine, cancelled)?;
    let mut error = 0.0_f64;
    for key in a.species.keys().chain(b.species.keys()) {
        poll(cancelled)?;
        for i in 0..4 {
            let x = a.species.get(key).map_or(0.0, |s| s[i].to_f64());
            let y = b.species.get(key).map_or(0.0, |s| s[i].to_f64());
            let absolute = if i == 3 {
                controls.stream_energy_tolerance_j
            } else {
                controls.species_mass_tolerance_kg
            };
            error = error.max(ratio(
                x,
                y,
                x.abs().max(y.abs()),
                absolute,
                controls.relative_tolerance,
            )?);
        }
    }
    for key in a.heat.keys().chain(b.heat.keys()) {
        poll(cancelled)?;
        for direction in 0..2 {
            let x = a.heat.get(key).map_or(0.0, |s| s[direction].to_f64());
            let y = b.heat.get(key).map_or(0.0, |s| s[direction].to_f64());
            error = error.max(ratio(
                x,
                y,
                x.max(y),
                controls.stream_energy_tolerance_j,
                controls.relative_tolerance,
            )?);
        }
    }
    Ok(error)
}
