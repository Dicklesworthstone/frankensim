//! Filling an exactly evacuated receiver through a finite thermal channel wall.
//!
//! The original finite-wall vacuum transaction supplies the first zero-return
//! beam. Its accepted escaped species and wall-adjusted energy are captured by
//! the original evacuated-volume caloric owner. That energy determines the new
//! receiver temperature; there is no seed gas, pressure floor or vacuum state.
//! The next call uses the ordinary two-gas finite-wall transaction, including
//! the receiver's actual return beam. Gas and every finite wall endpoint are atomic.
//!
//! This is a first-order, quasistatic stage with well-mixed ideal gases, rigid
//! fixed volumes, a stationary diffuse channel and negligible flight time. The
//! fixed vessels absorb molecular momentum without boundary work. Their motion,
//! angular momentum and spatial filling transient are absent. Original sourced
//! conduction may connect the channel wall to other finite material inventories.

/// Error-controlled startup and reciprocal evolution with the same finite wall.
pub mod trajectory;

pub use super::super::filling::EffusionReceiver;
use super::thermal::{
    self, FiniteChannelError, FiniteChannelNetworkLink, FiniteChannelNetworkStep,
    FixedVolumeChannelWall, FiniteChannelWallStep,
    conduction::ChannelWallConduction,
    vacuum::{self, FiniteVacuumChannelLink, FiniteVacuumChannelStep},
};
use crate::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve,
    parcel::{
        ClosedGasParcel,
        evacuated::{EvacuatedGasFillStep, EvacuatedGasVolume},
        transfer::species::{
            effusion::{channel::DiffuseMolecularChannel, positive_product_ratio},
            entropy,
            network::transport::{GasTransportPaths, GasTransportSolve},
        },
    },
};
use fs_math::dd::Dd;
use fs_thermal::enthalpy::ThermalLaw;
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use std::collections::BTreeMap;

/// Stable original passage, gas-volume and material-wall identities.
#[derive(Clone, Copy, Debug)]
pub struct FiniteChannelFillLink<'a> {
    pub id: u64,
    pub donor: u64,
    pub receiver: u64,
    pub wall: u64,
    pub channel: &'a DiffuseMolecularChannel,
    pub area_m2: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteChannelFirstFill {
    /// Original vacuum-stage account, captured INSIDE this finite receiver.
    /// Its escaped species/energy are not losses from the complete system.
    pub transport: FiniteVacuumChannelStep,
    pub fill: Option<EvacuatedGasFillStep>,
    pub captured_species_kg: BTreeMap<SpeciesId, f64>,
    pub captured_mass_kg: f64,
    /// Caloric chart anchor, never an assigned vacuum temperature.
    pub reference_temperature_k: f64,
    /// Captured U minus original species U at the chart anchor, including the
    /// original thermal channel's heat into the transmitted molecules [J].
    pub reference_relative_energy_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FiniteChannelFillTransfer {
    Initial(FiniteChannelFirstFill),
    Exchange(FiniteChannelNetworkStep),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteChannelFillBalance {
    pub total_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub gas_energy_change_j: f64,
    pub wall_energy_change_j: f64,
    pub energy_residual_j: f64,
    pub gas_entropy_change_j_k: f64,
    pub wall_entropy_change_j_k: f64,
    /// Complete closed gas-pair plus all finite-wall storage entropy change.
    pub total_entropy_change_j_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteChannelFillStep {
    pub duration_s: f64,
    pub transfer: FiniteChannelFillTransfer,
    pub balance: FiniteChannelFillBalance,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FiniteChannelFillError {
    Invalid(&'static str),
    Physical(FiniteChannelError),
    Receiver(GasEnergyError),
    Endpoint { id: u64, error: GasEnergyError },
    Conservation(FiniteChannelFillBalance),
    EntropyDecrease { change_j_k: f64, tolerance_j_k: f64 },
    Cancelled,
}
impl core::fmt::Display for FiniteChannelFillError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite thermal channel filling: {self:?}")
    }
}
impl std::error::Error for FiniteChannelFillError {}
type Error = FiniteChannelFillError;

#[derive(Clone, Debug)]
pub struct FiniteChannelFillCandidate<'l, 'p> {
    donor: ClosedGasParcel,
    receiver: EffusionReceiver,
    wall: FixedVolumeChannelWall<'l, 'p>,
    step: FiniteChannelFillStep,
}
impl<'l, 'p> FiniteChannelFillCandidate<'l, 'p> {
    pub const fn donor(&self) -> &ClosedGasParcel {
        &self.donor
    }
    pub const fn receiver(&self) -> &EffusionReceiver {
        &self.receiver
    }
    pub const fn wall(&self) -> &FixedVolumeChannelWall<'l, 'p> {
        &self.wall
    }
    pub const fn step(&self) -> &FiniteChannelFillStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ClosedGasParcel,
        EffusionReceiver,
        FixedVolumeChannelWall<'l, 'p>,
        FiniteChannelFillStep,
    ) {
        (self.donor, self.receiver, self.wall, self.step)
    }
}

/// Atomic finite receiver and every original material wall in its assembly.
#[derive(Clone, Debug)]
pub struct FiniteChannelFillAssemblyCandidate<'l, 'p> {
    donor: ClosedGasParcel,
    receiver: EffusionReceiver,
    walls: BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    step: FiniteChannelFillStep,
}
impl<'l, 'p> FiniteChannelFillAssemblyCandidate<'l, 'p> {
    pub const fn donor(&self) -> &ClosedGasParcel { &self.donor }
    pub const fn receiver(&self) -> &EffusionReceiver { &self.receiver }
    pub const fn walls(&self) -> &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>> { &self.walls }
    pub const fn step(&self) -> &FiniteChannelFillStep { &self.step }
    pub fn into_parts(self) -> (
        ClosedGasParcel, EffusionReceiver,
        BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>, FiniteChannelFillStep,
    ) {
        (self.donor, self.receiver, self.walls, self.step)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn physical(error: FiniteChannelError) -> Error {
    if error == FiniteChannelError::Cancelled {
        Error::Cancelled
    } else {
        Error::Physical(error)
    }
}
fn receiver_error(error: GasEnergyError) -> Error {
    if error == GasEnergyError::Cancelled {
        Error::Cancelled
    } else {
        Error::Receiver(error)
    }
}
fn endpoint(id: u64, error: GasEnergyError) -> Error {
    if error == GasEnergyError::Cancelled {
        Error::Cancelled
    } else {
        Error::Endpoint { id, error }
    }
}

#[allow(clippy::too_many_arguments)]
fn first_fill<'l, 'p>(
    donor: &ClosedGasParcel,
    receiver: &EvacuatedGasVolume,
    walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    link: FiniteChannelFillLink<'_>,
    conduction: &[ChannelWallConduction<'_, '_>],
    duration_s: f64,
    donor_solve: GasEnergySolve,
    receiver_solve: GasEnergySolve,
    controls: GasTransportSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<
    (
        ClosedGasParcel,
        EffusionReceiver,
        BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
        FiniteChannelFirstFill,
    ),
    Error,
> {
    let (mut gases, walls, transport) = vacuum::step_candidate_with_conduction(
        &BTreeMap::from([(link.donor, donor.clone())]),
        walls,
        GasTransportPaths::default(),
        &[FiniteVacuumChannelLink {
            id: link.id,
            parcel: link.donor,
            channel: link.channel,
            area_m2: link.area_m2,
            wall: link.wall,
        }],
        conduction,
        duration_s,
        &BTreeMap::from([(link.donor, donor_solve)]),
        controls,
        &mut *cancelled,
    )
    .map_err(physical)?
    .into_parts();
    let outlet = transport.gas.outlets.first().ok_or(Error::Invalid(
        "original thermal channel lost its captured outlet",
    ))?;
    let mut captured = BTreeMap::new();
    let mut mass = Dd::ZERO;
    let mut relative = Dd::from_f64(outlet.wall.heat_into_escaped_j);
    let reference_temperature_k = donor.state().state().temperature;
    for row in &outlet.escaped_species {
        poll(cancelled)?;
        if captured.insert(row.species.clone(), row.mass_kg).is_some() {
            return Err(Error::Invalid(
                "duplicate originally captured molecular species",
            ));
        }
        let source = link
            .channel
            .mouth()
            .gas_catalog()
            .thermodynamics()
            .components()
            .iter()
            .find(|source| source.species() == &row.species)
            .ok_or(Error::Invalid(
                "captured molecule lost its original caloric source",
            ))?;
        mass = mass + Dd::from_f64(row.mass_kg);
        // Original aperture e_i(T)=u_i(T)+R*T/(2*M_i), plus ONLY the
        // transmitted molecules' wall heat. Returned heat stays with donor.
        relative = relative
            + Dd::from_f64(
                positive_product_ratio(
                    &[row.mass_kg, 0.5 * R, reference_temperature_k],
                    &[source.molar_mass().value()],
                )
                .map_err(receiver_error)?,
            );
    }
    if !mass.to_f64().is_finite() || !relative.to_f64().is_finite() {
        return Err(Error::Invalid(
            "captured molecular energy or mass is unrepresentable",
        ));
    }
    let (next_receiver, fill) = if captured.is_empty() {
        (EffusionReceiver::Evacuated(receiver.clone()), None)
    } else {
        let (gas, fill) = receiver
            .fill_candidate(
                donor.model(),
                &captured,
                reference_temperature_k,
                relative.to_f64(),
                receiver_solve,
                controls.exchange.mass_tolerance_kg,
                controls.exchange.energy_tolerance_j,
                &mut *cancelled,
            )
            .map_err(receiver_error)?
            .into_parts();
        (EffusionReceiver::Filled(gas), Some(fill))
    };
    let next_donor = gases
        .remove(&link.donor)
        .ok_or(Error::Invalid("original thermal channel lost its donor"))?;
    Ok((
        next_donor,
        next_receiver,
        walls,
        FiniteChannelFirstFill {
            transport,
            fill,
            captured_species_kg: captured,
            captured_mass_kg: mass.to_f64(),
            reference_temperature_k,
            reference_relative_energy_j: relative.to_f64(),
        },
    ))
}

/// Fill an empty finite receiver, then continue through the original reciprocal
/// gas/wall owner on subsequent calls. Every stage samples one original state;
/// the caller refines time. Species capture preserves accepted masses exactly.
/// Source-relative gas storage and actual finite material storage close under
/// one fixed pair budget, including the newly filled receiver's true entropy.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'l, 'p>(
    donor: &ClosedGasParcel,
    receiver: &EffusionReceiver,
    wall: &FixedVolumeChannelWall<'l, 'p>,
    link: FiniteChannelFillLink<'_>,
    duration_s: f64,
    donor_solve: GasEnergySolve,
    receiver_solve: GasEnergySolve,
    controls: GasTransportSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<FiniteChannelFillCandidate<'l, 'p>, Error> {
    let (donor, receiver, mut walls, step) = step_candidate_with_conduction(
        donor, receiver, &BTreeMap::from([(link.wall, wall.clone())]), link, &[],
        duration_s, donor_solve, receiver_solve, controls, cancelled,
    )?.into_parts();
    let wall = walls.remove(&link.wall)
        .ok_or(Error::Invalid("finite filling lost its original single wall"))?;
    Ok(FiniteChannelFillCandidate { donor, receiver, wall, step })
}

/// Fill or continue the actual receiver while its finite material assembly conducts.
///
/// The original vacuum/reciprocal transactions pool all molecular and conductive
/// heat before one material inversion. Captured species and energy retain the
/// original evacuated-volume caloric owner; auxiliary wall heat is never added
/// to the receiver's molecular energy. Every original wall publishes atomically.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate_with_conduction<'l, 'p>(
    donor: &ClosedGasParcel,
    receiver: &EffusionReceiver,
    walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    link: FiniteChannelFillLink<'_>,
    conduction: &[ChannelWallConduction<'_, '_>],
    duration_s: f64,
    donor_solve: GasEnergySolve,
    receiver_solve: GasEnergySolve,
    controls: GasTransportSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<FiniteChannelFillAssemblyCandidate<'l, 'p>, Error> {
    poll(&mut cancelled)?;
    if link.donor == link.receiver
        || !duration_s.is_finite()
        || duration_s <= 0.0
        || [
            controls.exchange.mass_tolerance_kg,
            controls.exchange.energy_tolerance_j,
            controls.entropy_tolerance_j_k,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v <= 0.0)
    {
        return Err(Error::Invalid(
            "distinct gas identities, positive time and fixed physical budgets required",
        ));
    }
    receiver_solve.validate().map_err(receiver_error)?;
    if !walls.contains_key(&link.wall) {
        return Err(Error::Invalid("finite filling references an absent channel wall"));
    }
    let (next_donor, next_receiver, next_walls, transfer) = match receiver {
        EffusionReceiver::Evacuated(empty) => {
            let (gas, receiver, walls, first) = first_fill(
                donor,
                empty,
                walls,
                link,
                conduction,
                duration_s,
                donor_solve,
                receiver_solve,
                controls,
                &mut cancelled,
            )?;
            (
                gas,
                receiver,
                walls,
                FiniteChannelFillTransfer::Initial(first),
            )
        }
        EffusionReceiver::Filled(gas) => {
            let (mut gases, walls, report) = thermal::step_candidate_with_conduction(
                &BTreeMap::from([(link.donor, donor.clone()), (link.receiver, gas.clone())]),
                walls,
                GasTransportPaths::default(),
                &[FiniteChannelNetworkLink {
                    id: link.id,
                    parcels: [link.donor, link.receiver],
                    channel: link.channel,
                    area_m2: link.area_m2,
                    wall: link.wall,
                }],
                conduction,
                duration_s,
                &BTreeMap::from([(link.donor, donor_solve), (link.receiver, receiver_solve)]),
                controls,
                &mut cancelled,
            )
            .map_err(physical)?
            .into_parts();
            (
                gases
                    .remove(&link.donor)
                    .ok_or(Error::Invalid("finite channel lost donor"))?,
                EffusionReceiver::Filled(
                    gases
                        .remove(&link.receiver)
                        .ok_or(Error::Invalid("finite channel lost receiver"))?,
                ),
                walls,
                FiniteChannelFillTransfer::Exchange(report),
            )
        }
    };
    // The initial vacuum owner cannot admit an as-yet nonexistent receiving
    // gas. Check that actual new gas against the SAME channel and final wall.
    if link.area_m2 > 0.0 {
        if let Some(gas) = next_receiver.gas() {
            link.channel
                .admit_wall_for_parcel(gas, next_walls[&link.wall].temperature_k(), &mut cancelled)
                .map_err(|error| endpoint(link.receiver, error))?;
        }
    }
    let balance = audit_assembly(
        donor,
        receiver,
        walls,
        &next_donor,
        &next_receiver,
        &next_walls,
        link,
        controls,
        &mut cancelled,
    )?;
    poll(&mut cancelled)?;
    Ok(FiniteChannelFillAssemblyCandidate {
        donor: next_donor,
        receiver: next_receiver,
        walls: next_walls,
        step: FiniteChannelFillStep {
            duration_s,
            transfer,
            balance,
        },
    })
}

// Shared by one physical stage and the complete accepted trajectory. The
// original empty volume contributes no invented thermodynamic state.
#[allow(clippy::too_many_arguments)]
fn audit_assembly(
    donor: &ClosedGasParcel,
    receiver: &EffusionReceiver,
    walls: &BTreeMap<u64, FixedVolumeChannelWall<'_, '_>>,
    next_donor: &ClosedGasParcel,
    next_receiver: &EffusionReceiver,
    next_walls: &BTreeMap<u64, FixedVolumeChannelWall<'_, '_>>,
    link: FiniteChannelFillLink<'_>,
    controls: GasTransportSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<FiniteChannelFillBalance, Error> {
    poll(cancelled)?;
    if !walls.keys().eq(next_walls.keys()) {
        return Err(Error::Invalid("finite filling changed its original wall identities"));
    }
    let gas = match (receiver.gas(), next_receiver.gas()) {
        (Some(before), Some(after)) => entropy::closed_inventory_change(
            [
                (link.donor, donor, next_donor),
                (link.receiver, before, after),
            ],
            controls.exchange.mass_tolerance_kg,
            &mut *cancelled,
        ),
        (None, Some(after)) => entropy::newly_filled_pair_change(
            donor,
            (link.donor, next_donor),
            (link.receiver, after),
            controls.exchange.mass_tolerance_kg,
            &mut *cancelled,
        ),
        (None, None) => entropy::closed_inventory_change(
            [(link.donor, donor, next_donor)],
            controls.exchange.mass_tolerance_kg,
            &mut *cancelled,
        ),
        (Some(_), None) => return Err(Error::Invalid("finite receiver was silently evacuated")),
    }
    .map_err(|error| match error {
        entropy::Error::Invalid(why) => Error::Invalid(why),
        entropy::Error::Parcel { id, error } => endpoint(id, error),
        entropy::Error::Cancelled => Error::Cancelled,
    })?;
    let mut wall_energy = Dd::ZERO;
    let mut wall_entropy = Dd::ZERO;
    for (&id, wall) in walls {
        poll(cancelled)?;
        let next_wall = &next_walls[&id];
        if wall.mass_kg().to_bits() != next_wall.mass_kg().to_bits()
            || wall.volume_m3().to_bits() != next_wall.volume_m3().to_bits()
            || wall.law().identity() != next_wall.law().identity()
        {
            return Err(Error::Invalid("finite filling changed an original material inventory"));
        }
        wall_energy = wall_energy
            + (Dd::from_f64(next_wall.specific_energy_j_kg())
                - Dd::from_f64(wall.specific_energy_j_kg())) * Dd::from_f64(wall.mass_kg());
        let ds = wall.law().specific_entropy_change_j_kg_k(
            wall.specific_energy_j_kg(), next_wall.specific_energy_j_kg(),
        ).map_err(|error| physical(FiniteChannelError::WallDomain { id, error }))?
            .ok_or(Error::Invalid("finite wall entropy law became unavailable"))?;
        wall_entropy = wall_entropy + Dd::from_f64(ds) * Dd::from_f64(wall.mass_kg());
    }
    let balance = FiniteChannelFillBalance {
        total_mass_residual_kg: gas.total_mass_residual_kg,
        maximum_species_mass_residual_kg: gas.maximum_species_mass_residual_kg,
        gas_energy_change_j: gas.energy_change_j,
        wall_energy_change_j: wall_energy.to_f64(),
        energy_residual_j: (Dd::from_f64(gas.energy_change_j) + wall_energy).to_f64(),
        gas_entropy_change_j_k: gas.entropy_change_j_k,
        wall_entropy_change_j_k: wall_entropy.to_f64(),
        total_entropy_change_j_k: (Dd::from_f64(gas.entropy_change_j_k) + wall_entropy).to_f64(),
    };
    if ![
        balance.wall_energy_change_j,
        balance.energy_residual_j,
        balance.wall_entropy_change_j_k,
        balance.total_entropy_change_j_k,
    ]
    .into_iter()
    .all(f64::is_finite)
        || balance.total_mass_residual_kg.abs() > controls.exchange.mass_tolerance_kg
        || balance.maximum_species_mass_residual_kg > controls.exchange.mass_tolerance_kg
        || balance.energy_residual_j.abs() > controls.exchange.energy_tolerance_j
    {
        return Err(Error::Conservation(balance));
    }
    if balance.total_entropy_change_j_k < -controls.entropy_tolerance_j_k {
        return Err(Error::EntropyDecrease {
            change_j_k: balance.total_entropy_change_j_k,
            tolerance_j_k: controls.entropy_tolerance_j_k,
        });
    }
    poll(cancelled)?;
    Ok(balance)
}

/// One original-to-final storage allowance per wall, independent of step count.
fn audit_wall_history<'r>(
    walls: &BTreeMap<u64, FixedVolumeChannelWall<'_, '_>>,
    next_walls: &BTreeMap<u64, FixedVolumeChannelWall<'_, '_>>,
    reports: impl Iterator<Item = &'r FiniteChannelWallStep>,
    tolerance_j: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    let mut heat = BTreeMap::<u64, Dd>::new();
    for row in reports {
        poll(cancelled)?;
        if !walls.contains_key(&row.id) {
            return Err(Error::Invalid("finite filling history lost an original wall"));
        }
        let sum = heat.entry(row.id).or_insert(Dd::ZERO);
        *sum = *sum + Dd::from_f64(row.heat_into_wall_j);
    }
    for (&id, wall) in walls {
        poll(cancelled)?;
        let next = next_walls.get(&id)
            .ok_or(Error::Invalid("finite filling history lost a wall endpoint"))?;
        let stored = (Dd::from_f64(next.specific_energy_j_kg())
            - Dd::from_f64(wall.specific_energy_j_kg())) * Dd::from_f64(wall.mass_kg());
        let residual_j = (stored - heat.get(&id).copied().unwrap_or(Dd::ZERO)).to_f64();
        if !residual_j.is_finite() || residual_j.abs() > tolerance_j {
            return Err(physical(FiniteChannelError::WallEnergyBalance { id, residual_j }));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
