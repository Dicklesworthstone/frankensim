//! Directional error estimates and fixed-budget gas/thermostat window accounts.

use super::super::super::trajectory::ratio;
use super::{
    BTreeMap, ClosedGasParcel, DiffuseChannelNetworkStep, EffusionAdvanceControls, Error,
    GasTransportError, GasTransportSolve, MAX_EFFUSION_RETAINED_ROWS, physical, poll,
};
use crate::thermochemical_gas::energy::parcel::transfer::species::{
    diffusion::network::DiffusionNetworkError, entropy,
};
use fs_math::dd::Dd;
use fs_thermochem::SpeciesId;

fn add(total: &mut Dd, value: f64) -> Result<(), Error> {
    *total = *total + Dd::from_f64(value);
    if !value.is_finite() || !total.hi.is_finite() || !total.lo.is_finite() {
        return Err(Error::Invalid(
            "integrated channel account is unrepresentable",
        ));
    }
    Ok(())
}

/// Actual integrated thermostat exchange, positive into the named gas or wall.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffuseChannelIntegratedWall {
    pub id: u64,
    pub parcels: [u64; 2],
    pub wall_temperature_k: f64,
    pub heat_into_gas_j: [f64; 2],
    pub heat_into_wall_j: f64,
    pub entropy_change_j_k: f64,
}

/// Original-to-final storage and accepted original streams, with fixed budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiffuseChannelWindowBalance {
    pub total_mass_residual_kg: f64,
    pub maximum_parcel_species_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub energy_residual_j: f64,
    pub accumulated_energy_residual_j: f64,
    pub heat_into_walls_j: f64,
    pub gas_entropy_change_j_k: f64,
    pub wall_entropy_change_j_k: f64,
    pub total_entropy_change_j_k: f64,
}

struct WallTotal {
    parcels: [u64; 2],
    temperature: f64,
    heat: [Dd; 3],
    entropy: Dd,
}
#[derive(Default)]
pub(super) struct Accepted {
    inventory_change: BTreeMap<(u64, SpeciesId), Dd>,
    walls: BTreeMap<u64, WallTotal>,
    native_residual: Dd,
}
impl Accepted {
    pub(super) fn record(
        &mut self,
        step: &DiffuseChannelNetworkStep,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
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
        for boundary in &step.exchange.boundaries {
            poll(cancelled)?;
            if boundary.transfer.net_mass_into_gas_kg != 0.0
                || boundary.transfer.relative_energy_into_gas_j != 0.0
            {
                return Err(Error::Invalid(
                    "held-wall channel history contains an external material stream",
                ));
            }
        }
        for channel in &step.channels {
            poll(cancelled)?;
            let total = self.walls.entry(channel.id).or_insert_with(|| WallTotal {
                parcels: channel.parcels,
                temperature: channel.wall_temperature_k,
                heat: [Dd::ZERO; 3],
                entropy: Dd::ZERO,
            });
            // Closed channels may deliberately carry an inactive NaN wall.
            // Identity is bitwise; no dormant thermal coordinate is queried.
            if total.parcels != channel.parcels
                || total.temperature.to_bits() != channel.wall_temperature_k.to_bits()
            {
                return Err(Error::Invalid(
                    "channel window changed its original wall identity",
                ));
            }
            for (sum, value) in total.heat.iter_mut().zip([
                channel.wall.heat_into_gas_j[0],
                channel.wall.heat_into_gas_j[1],
                channel.wall.heat_into_wall_j,
            ]) {
                add(sum, value)?;
            }
            if channel.wall.heat_into_wall_j != 0.0 {
                add(
                    &mut total.entropy,
                    channel.wall.heat_into_wall_j / channel.wall_temperature_k,
                )?;
            }
        }
        add(
            &mut self.native_residual,
            step.exchange.total_energy_residual_j,
        )?;
        if self.inventory_change.len() + self.walls.len() > MAX_EFFUSION_RETAINED_ROWS {
            return Err(Error::Invalid(
                "channel accounts exceed retained work bounds",
            ));
        }
        poll(cancelled)
    }

    pub(super) fn finish(
        self,
        before: &BTreeMap<u64, ClosedGasParcel>,
        after: &BTreeMap<u64, ClosedGasParcel>,
        duration_s: f64,
        controls: GasTransportSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<
        (
            Vec<DiffuseChannelIntegratedWall>,
            DiffuseChannelWindowBalance,
        ),
        Error,
    > {
        if !before.keys().eq(after.keys()) {
            return Err(Error::Invalid(
                "channel window lost original parcel identities",
            ));
        }
        let mut residuals: BTreeMap<_, _> = self
            .inventory_change
            .into_iter()
            .map(|(key, mass)| (key, -mass))
            .collect();
        for (&id, old) in before {
            poll(cancelled)?;
            for (gas, sign) in [(old, -1.0), (&after[&id], 1.0)] {
                for (species, mass) in gas
                    .species_inventory_kg()
                    .map_err(|_| Error::Invalid("channel window inventory is unrepresentable"))?
                {
                    poll(cancelled)?;
                    add(residuals.entry((id, species)).or_default(), sign * mass)?;
                }
            }
            if residuals.len() > MAX_EFFUSION_RETAINED_ROWS {
                return Err(Error::Invalid(
                    "channel inventory audit exceeds retained work bounds",
                ));
            }
        }
        let maximum_parcel_species_mass_residual_kg = residuals
            .values()
            .map(|value| value.to_f64().abs())
            .fold(0.0_f64, f64::max);
        let change = entropy::closed_inventory_change(
            before.iter().map(|(&id, old)| (id, old, &after[&id])),
            controls.exchange.mass_tolerance_kg,
            &mut *cancelled,
        )
        .map_err(|error| {
            physical(
                match error {
                    entropy::Error::Invalid(reason) => GasTransportError::Invalid(reason),
                    entropy::Error::Parcel { id, error } => {
                        GasTransportError::Diffusion(DiffusionNetworkError::Parcel { id, error })
                    }
                    entropy::Error::Cancelled => GasTransportError::Cancelled,
                },
                duration_s,
                duration_s,
            )
        })?;

        let (mut heat, mut wall_entropy) = (Dd::ZERO, Dd::ZERO);
        let mut walls = Vec::with_capacity(self.walls.len());
        for (id, total) in self.walls {
            poll(cancelled)?;
            // Preserve compensated per-wall sums through the global balance.
            heat = heat + total.heat[2];
            wall_entropy = wall_entropy + total.entropy;
            walls.push(DiffuseChannelIntegratedWall {
                id,
                parcels: total.parcels,
                wall_temperature_k: total.temperature,
                heat_into_gas_j: [total.heat[0].to_f64(), total.heat[1].to_f64()],
                heat_into_wall_j: total.heat[2].to_f64(),
                entropy_change_j_k: total.entropy.to_f64(),
            });
        }
        let energy_residual_j = (Dd::from_f64(change.energy_change_j) + heat).to_f64();
        let accumulated_energy_residual_j = self.native_residual.to_f64();
        let total_entropy_change_j_k =
            (Dd::from_f64(change.entropy_change_j_k) + wall_entropy).to_f64();
        if ![
            energy_residual_j,
            accumulated_energy_residual_j,
            heat.to_f64(),
            wall_entropy.to_f64(),
            total_entropy_change_j_k,
            maximum_parcel_species_mass_residual_kg,
        ]
        .into_iter()
        .all(f64::is_finite)
            || change.total_mass_residual_kg.abs() > controls.exchange.mass_tolerance_kg
            || change.maximum_species_mass_residual_kg > controls.exchange.mass_tolerance_kg
            || maximum_parcel_species_mass_residual_kg > controls.exchange.mass_tolerance_kg
            || energy_residual_j.abs() > controls.exchange.energy_tolerance_j
            || accumulated_energy_residual_j.abs() > controls.exchange.energy_tolerance_j
        {
            return Err(Error::Invalid(
                "channel window exceeded its fixed mass or energy budget",
            ));
        }
        if total_entropy_change_j_k < -controls.entropy_tolerance_j_k {
            return Err(physical(
                GasTransportError::Diffusion(DiffusionNetworkError::EntropyDecrease {
                    change_j_k: total_entropy_change_j_k,
                    tolerance_j_k: controls.entropy_tolerance_j_k,
                }),
                duration_s,
                duration_s,
            ));
        }
        poll(cancelled)?;
        Ok((
            walls,
            DiffuseChannelWindowBalance {
                total_mass_residual_kg: change.total_mass_residual_kg,
                maximum_parcel_species_mass_residual_kg,
                maximum_species_mass_residual_kg: change.maximum_species_mass_residual_kg,
                energy_residual_j,
                accumulated_energy_residual_j,
                heat_into_walls_j: heat.to_f64(),
                gas_entropy_change_j_k: change.entropy_change_j_k,
                wall_entropy_change_j_k: wall_entropy.to_f64(),
                total_entropy_change_j_k,
            },
        ))
    }
}

pub(super) fn row_count(step: &DiffuseChannelNetworkStep) -> usize {
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
        + rates.layered_dusty_gas.len()
        + rates
            .effusion
            .iter()
            .map(|row| 1 + row.flux.species.len())
            .sum::<usize>()
        + step
            .channels
            .iter()
            .map(|row| 1 + row.transmitted.species.len() + row.wall.reflected_species.len())
            .sum::<usize>()
}

// Namespaces distinguish internal transfers, reflected withdrawals, internal
// heat, external accommodation and thermostat heat without merging directions.
type SpeciesKey = (u8, u64, u64, SpeciesId);
type HeatKey = (u8, u64, u64);
#[derive(Default)]
struct Streams {
    species: BTreeMap<SpeciesKey, [Dd; 2]>,
    heat: BTreeMap<HeatKey, [Dd; 2]>,
}
impl Streams {
    fn record_heat(&mut self, key: HeatKey, value: f64) -> Result<(), Error> {
        let direction = usize::from(value < 0.0);
        let sums = self.heat.entry(key).or_insert([Dd::ZERO; 2]);
        add(&mut sums[direction], value.abs())
    }

    fn record(
        &mut self,
        step: &DiffuseChannelNetworkStep,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        for flux in &step.exchange.fluxes {
            poll(cancelled)?;
            let row = &flux.transfer;
            let sums = self
                .species
                .entry((0, row.donor, row.receiver, row.species.clone()))
                .or_insert([Dd::ZERO; 2]);
            add(&mut sums[0], row.mass_kg)?;
        }
        for flux in &step.exchange.boundaries {
            poll(cancelled)?;
            let row = &flux.transfer;
            if row.net_mass_into_gas_kg != 0.0 || row.relative_energy_into_gas_j != 0.0 {
                return Err(Error::Invalid(
                    "channel error estimate contains an external material stream",
                ));
            }
            let sums = self
                .species
                .entry((1, flux.id, flux.parcel, row.species.clone()))
                .or_insert([Dd::ZERO; 2]);
            add(&mut sums[0], row.incoming_mass_kg)?;
            add(&mut sums[1], row.outgoing_mass_kg)?;
        }
        for row in &step.exchange.heat {
            poll(cancelled)?;
            self.record_heat(
                (0, row.parcels[0], row.parcels[1]),
                row.heat_left_to_right_j,
            )?;
        }
        for row in &step.exchange.external_heat {
            poll(cancelled)?;
            self.record_heat((1, row.id, row.parcel), row.heat_into_gas_j)?;
        }
        for row in &step.channels {
            poll(cancelled)?;
            self.record_heat((2, row.id, 0), row.wall.heat_into_wall_j)?;
        }
        if self.species.len() + self.heat.len() > MAX_EFFUSION_RETAINED_ROWS {
            return Err(Error::Invalid(
                "channel error estimate exceeds its work bound",
            ));
        }
        poll(cancelled)
    }
}

pub(super) fn stream_error(
    coarse: &DiffuseChannelNetworkStep,
    first: &DiffuseChannelNetworkStep,
    second: &DiffuseChannelNetworkStep,
    controls: EffusionAdvanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let (mut a, mut b) = (Streams::default(), Streams::default());
    a.record(coarse, cancelled)?;
    b.record(first, cancelled)?;
    b.record(second, cancelled)?;
    let mut error = 0.0_f64;
    for key in a.species.keys().chain(b.species.keys()) {
        poll(cancelled)?;
        let x = a.species.get(key).copied().unwrap_or([Dd::ZERO; 2]);
        let y = b.species.get(key).copied().unwrap_or([Dd::ZERO; 2]);
        for (x, y) in x.into_iter().zip(y) {
            let (x, y) = (x.to_f64(), y.to_f64());
            error = error.max(ratio(
                x,
                y,
                x.abs().max(y.abs()),
                controls.species_mass_tolerance_kg,
                controls.relative_tolerance,
            )?);
        }
    }
    for key in a.heat.keys().chain(b.heat.keys()) {
        poll(cancelled)?;
        let x = a.heat.get(key).copied().unwrap_or([Dd::ZERO; 2]);
        let y = b.heat.get(key).copied().unwrap_or([Dd::ZERO; 2]);
        for (x, y) in x.into_iter().zip(y) {
            let (x, y) = (x.to_f64(), y.to_f64());
            error = error.max(ratio(
                x,
                y,
                x.max(y),
                controls.stream_energy_tolerance_j,
                controls.relative_tolerance,
            )?);
        }
    }
    poll(cancelled)?;
    Ok(error)
}
