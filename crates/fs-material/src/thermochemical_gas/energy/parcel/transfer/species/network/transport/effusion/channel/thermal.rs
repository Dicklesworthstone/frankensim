//! Finite fixed-volume material walls for diffuse molecular channels.
//!
//! Channels sample one original gas/wall snapshot. Every channel sharing a wall
//! contributes heat to the same pool before that wall is reconstructed. The
//! original material law owns capacity, inversion and sensible entropy; no
//! constant heat capacity or thermostat entropy is substituted.
//!
//! This is an explicit, first-order stage with fixed composition and geometry.
//! The caller must select volumetric capacity appropriate to reversible,
//! isochoric sensible storage. The fixed-volume source adapter does not convert
//! an isobaric Cp or an apparent capacity measurement into Cv. Expansion work,
//! phase change, wall conduction, adsorption and wall mass exchange are absent.
//! An unsupported or overlarge step refuses atomically; refine the duration.

/// Adaptive gas-plus-material-wall trajectories with fixed window audits.
pub mod trajectory;

use super::{DiffuseChannelLinkStep, DiffuseChannelNetworkLink};
use crate::thermal::transport::{
    ThermalTransportError, volumetric::timeline::FixedVolumeThermalLaw,
};
use crate::thermochemical_gas::energy::{
    GasEnergySolve,
    parcel::{
        ClosedGasParcel,
        transfer::{
            snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
            species::{
                effusion::channel::DiffuseMolecularChannel,
                entropy,
                network::{
                    SpeciesNetworkStep,
                    transport::{
                        GasTransportError, GasTransportPaths, GasTransportRates, GasTransportSolve,
                    },
                },
            },
        },
    },
};
use fs_blake3::ContentHash;
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};
use std::collections::BTreeMap;

/// Immutable wall inventory on one original fixed-volume caloric law.
///
/// The law, reference volume and mass remain bound across candidates. The
/// stored J/kg coordinate is interpreted as internal energy under the explicit
/// isochoric contract of `try_new_isochoric`; its zero is the law's own datum.
#[derive(Clone, Debug)]
pub struct FixedVolumeChannelWall<'l, 'p> {
    law: &'l FixedVolumeThermalLaw<'p>,
    volume_m3: f64,
    mass_kg: f64,
    specific_energy_j_kg: f64,
    temperature_k: f64,
}
impl<'l, 'p> FixedVolumeChannelWall<'l, 'p> {
    /// Bind an isochoric sensible material inventory without a new source query.
    ///
    /// Calling this constructor declares that the original capacity is
    /// appropriate to the rigid wall's fixed composition and work constraint.
    /// Select `with_sensible_entropy()` on the original law first. An isobaric
    /// chart cannot be made isochoric merely by selecting this interpretation.
    pub fn try_new_isochoric(
        law: &'l FixedVolumeThermalLaw<'p>,
        volume_m3: f64,
        temperature_k: f64,
    ) -> Result<Self, FiniteChannelError> {
        if !law.sensible_entropy_enabled() {
            return Err(FiniteChannelError::Invalid(
                "finite channel walls require the original sensible entropy law",
            ));
        }
        let mass_kg = law.reference_mass_kg(volume_m3)?;
        let specific_energy_j_kg = law.specific_storage_j_kg(temperature_k)?;
        let actual_temperature = law.temperature_k(specific_energy_j_kg)?;
        if !actual_temperature.is_finite() || actual_temperature <= 0.0 {
            return Err(FiniteChannelError::Invalid(
                "invalid material wall temperature",
            ));
        }
        Ok(Self {
            law,
            volume_m3,
            mass_kg,
            specific_energy_j_kg,
            temperature_k: actual_temperature,
        })
    }

    #[must_use]
    pub const fn law(&self) -> &'l FixedVolumeThermalLaw<'p> {
        self.law
    }
    #[must_use]
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    #[must_use]
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    #[must_use]
    pub const fn specific_energy_j_kg(&self) -> f64 {
        self.specific_energy_j_kg
    }
    #[must_use]
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
}

/// Each physical channel refers to one stable finite wall identity.
#[derive(Clone, Copy, Debug)]
pub struct FiniteChannelNetworkLink<'a> {
    pub id: u64,
    pub parcels: [u64; 2],
    pub channel: &'a DiffuseMolecularChannel,
    pub area_m2: f64,
    pub wall: u64,
}

/// Actual original-to-final material storage and entropy, in wall-ID order.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteChannelWallStep {
    pub id: u64,
    pub caloric_identity: ContentHash,
    pub volume_m3: f64,
    pub mass_kg: f64,
    pub temperatures_k: [f64; 2],
    pub specific_energies_j_kg: [f64; 2],
    /// Pooled molecular accommodation heat, positive into this wall.
    pub heat_into_wall_j: f64,
    /// Actual stored change, which is independently checked against pooled heat.
    pub energy_change_j: f64,
    /// Original material integral of du/T, multiplied by the fixed mass.
    pub entropy_change_j_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteChannelNetworkStep {
    pub duration_s: f64,
    pub transport: GasTransportRates,
    /// Wall temperatures in these kinetic samples are the initial temperatures.
    pub channels: Vec<DiffuseChannelLinkStep>,
    /// Channel ID to its original wall ID; stable under caller input reordering.
    pub channel_walls: BTreeMap<u64, u64>,
    pub exchange: SpeciesNetworkStep,
    pub walls: Vec<FiniteChannelWallStep>,
    pub heat_into_walls_j: f64,
    pub wall_energy_change_j: f64,
    /// Actual total gas plus material-wall storage change in the closed system.
    pub energy_residual_j: f64,
    pub gas_entropy_change_j_k: f64,
    pub wall_entropy_change_j_k: f64,
    /// Actual gas entropy plus finite material entropy; contains no Q/T bath.
    pub total_entropy_change_j_k: f64,
}

/// Gas and wall endpoints are published together only after all admissions.
#[derive(Clone, Debug)]
pub struct FiniteChannelNetworkCandidate<'l, 'p> {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    walls: BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    step: FiniteChannelNetworkStep,
}
impl<'l, 'p> FiniteChannelNetworkCandidate<'l, 'p> {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub fn walls(&self) -> &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>> {
        &self.walls
    }
    pub const fn step(&self) -> &FiniteChannelNetworkStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, ClosedGasParcel>,
        BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
        FiniteChannelNetworkStep,
    ) {
        (self.parcels, self.walls, self.step)
    }
}

/// A refusal never changes the caller's gas, wall, or material sources.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteChannelError {
    Invalid(&'static str),
    Transport(GasTransportError),
    Material(ThermalTransportError),
    Thermal(DepositionError),
    WallDomain {
        id: u64,
        error: DepositionError,
    },
    /// The explicit wall endpoint left the original connected temperature range.
    WallOvershoot {
        id: u64,
        temperature_k: f64,
        range_k: [f64; 2],
    },
    WallEnergyBalance {
        id: u64,
        residual_j: f64,
    },
    TotalEnergyBalance {
        residual_j: f64,
    },
    EntropyDecrease {
        change_j_k: f64,
        tolerance_j_k: f64,
    },
    Cancelled,
}
impl std::fmt::Display for FiniteChannelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "finite molecular channel refusal: {self:?}")
    }
}
impl std::error::Error for FiniteChannelError {}
impl From<GasTransportError> for FiniteChannelError {
    fn from(error: GasTransportError) -> Self {
        if matches!(error, GasTransportError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Transport(error)
        }
    }
}
impl From<ThermalTransportError> for FiniteChannelError {
    fn from(error: ThermalTransportError) -> Self {
        Self::Material(error)
    }
}
impl From<DepositionError> for FiniteChannelError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}

/// Couple all native gas paths and molecular channels to finite material walls.
///
/// Reuses the original channel stage, including gross incident withdrawals and
/// one caloric solve per gas. Its initial-temperature thermostat account is only
/// an intermediate necessary admission. The returned receipt replaces it with
/// independently evaluated actual gas and finite-wall storage and entropy.
/// Shared walls are reconstructed once after all incident heat is pooled; two
/// individually unsupported heats may therefore cancel before any inversion.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'l, 'p>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    paths: GasTransportPaths<'_>,
    links: &[FiniteChannelNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: GasTransportSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<FiniteChannelNetworkCandidate<'l, 'p>, FiniteChannelError> {
    check(&mut cancelled)?;
    if walls.len() > MAX_SNAPSHOT_PARCELS || links.len() > MAX_SNAPSHOT_TRANSFERS {
        return Err(FiniteChannelError::Invalid(
            "finite wall/channel work budget exceeded",
        ));
    }
    let mut ordered: Vec<_> = links.iter().collect();
    ordered.sort_by_key(|link| link.id);
    let mut channel_walls = BTreeMap::new();
    let mut pools = BTreeMap::<u64, Dd>::new();
    let mut ranges = BTreeMap::<u64, [f64; 2]>::new();
    let mut held = Vec::with_capacity(ordered.len());
    for link in &ordered {
        check(&mut cancelled)?;
        if channel_walls.insert(link.id, link.wall).is_some() {
            return Err(FiniteChannelError::Invalid(
                "duplicate finite channel identity",
            ));
        }
        let wall = walls.get(&link.wall).ok_or(FiniteChannelError::Invalid(
            "finite channel references an absent material wall",
        ))?;
        let range = ranges.entry(link.wall).or_insert([wall.temperature_k; 2]);
        if link.area_m2 > 0.0 && link.channel.geometry().length_m() > 0.0 {
            for id in link.parcels {
                let gas = parcels.get(&id).ok_or(FiniteChannelError::Invalid(
                    "finite channel references an absent gas",
                ))?;
                let temperature = gas.state().state().temperature;
                range[0] = range[0].min(temperature);
                range[1] = range[1].max(temperature);
            }
        }
        pools.entry(link.wall).or_insert(Dd::from_f64(0.0));
        held.push(DiffuseChannelNetworkLink {
            id: link.id,
            parcels: link.parcels,
            channel: link.channel,
            area_m2: link.area_m2,
            wall_temperature_k: wall.temperature_k,
        });
    }
    let gas_candidate = super::step_candidate(
        parcels,
        paths,
        &held,
        duration_s,
        endpoints,
        controls,
        &mut cancelled,
    )?;
    for channel in &gas_candidate.step().channels {
        check(&mut cancelled)?;
        let id = channel_walls[&channel.id];
        let pooled = pools.get_mut(&id).expect("admitted original wall identity");
        *pooled = *pooled + Dd::from_f64(channel.wall.heat_into_wall_j);
    }

    let mut next_walls = walls.clone();
    let mut reports = Vec::with_capacity(pools.len());
    let mut wall_energy = Dd::from_f64(0.0);
    let mut wall_entropy = Dd::from_f64(0.0);
    let mut wall_heat = Dd::from_f64(0.0);
    for (&id, pooled) in &pools {
        check(&mut cancelled)?;
        let old = &walls[&id];
        let heat = pooled.to_f64();
        let new_energy =
            (Dd::from_f64(old.specific_energy_j_kg) + *pooled / Dd::from_f64(old.mass_kg)).to_f64();
        if !heat.is_finite()
            || !new_energy.is_finite()
            || (heat != 0.0 && new_energy == old.specific_energy_j_kg)
        {
            return Err(FiniteChannelError::Invalid(
                "finite wall heat is not representable",
            ));
        }
        let temperature = old
            .law
            .temperature_k(new_energy)
            .map_err(|error| FiniteChannelError::WallDomain { id, error })?;
        let range = ranges[&id];
        if temperature < range[0] || temperature > range[1] {
            return Err(FiniteChannelError::WallOvershoot {
                id,
                temperature_k: temperature,
                range_k: range,
            });
        }
        let energy = (Dd::from_f64(new_energy) - Dd::from_f64(old.specific_energy_j_kg))
            * Dd::from_f64(old.mass_kg);
        let residual = (energy - *pooled).to_f64();
        if !residual.is_finite() || residual.abs() > controls.exchange.energy_tolerance_j {
            return Err(FiniteChannelError::WallEnergyBalance {
                id,
                residual_j: residual,
            });
        }
        let specific_entropy = old
            .law
            .specific_entropy_change_j_kg_k(old.specific_energy_j_kg, new_energy)
            .map_err(|error| FiniteChannelError::WallDomain { id, error })?
            .ok_or(FiniteChannelError::Invalid(
                "original wall entropy became unavailable",
            ))?;
        let entropy_change = Dd::from_f64(specific_entropy) * Dd::from_f64(old.mass_kg);
        if !entropy_change.to_f64().is_finite() {
            return Err(FiniteChannelError::Invalid(
                "finite wall entropy is not representable",
            ));
        }
        wall_heat = wall_heat + *pooled;
        wall_energy = wall_energy + energy;
        wall_entropy = wall_entropy + entropy_change;
        let mut next = old.clone();
        next.specific_energy_j_kg = new_energy;
        next.temperature_k = temperature;
        next_walls.insert(id, next);
        reports.push(FiniteChannelWallStep {
            id,
            caloric_identity: old.law.identity(),
            volume_m3: old.volume_m3,
            mass_kg: old.mass_kg,
            temperatures_k: [old.temperature_k, temperature],
            specific_energies_j_kg: [old.specific_energy_j_kg, new_energy],
            heat_into_wall_j: heat,
            energy_change_j: energy.to_f64(),
            entropy_change_j_k: entropy_change.to_f64(),
        });
    }
    for link in &ordered {
        check(&mut cancelled)?;
        if link.area_m2 == 0.0 || link.channel.geometry().length_m() == 0.0 {
            continue;
        }
        let temperature = next_walls[&link.wall].temperature_k;
        for id in link.parcels {
            link.channel
                .admit_wall_for_parcel(&gas_candidate.parcels()[&id], temperature, &mut cancelled)
                .map_err(|error| {
                    GasTransportError::from(super::super::face_error(link.id, error))
                })?;
        }
    }
    let gas_change = entropy::closed_inventory_change(
        parcels.iter().map(|(&id, old)| (id, old, &gas_candidate.parcels()[&id])),
        controls.exchange.mass_tolerance_kg, &mut cancelled,
    ).map_err(|error| match error {
        entropy::Error::Invalid(reason) => FiniteChannelError::Invalid(reason),
        entropy::Error::Parcel { id, error } => FiniteChannelError::Transport(
            GasTransportError::Diffusion(
                crate::thermochemical_gas::energy::parcel::transfer::species::diffusion::network::DiffusionNetworkError::Parcel { id, error }
            )
        ),
        entropy::Error::Cancelled => FiniteChannelError::Cancelled,
    })?;
    let energy_residual_j = (Dd::from_f64(gas_change.energy_change_j) + wall_energy).to_f64();
    let total_entropy_change_j_k =
        (Dd::from_f64(gas_change.entropy_change_j_k) + wall_entropy).to_f64();
    if !energy_residual_j.is_finite()
        || energy_residual_j.abs() > controls.exchange.energy_tolerance_j
    {
        return Err(FiniteChannelError::TotalEnergyBalance {
            residual_j: energy_residual_j,
        });
    }
    if !total_entropy_change_j_k.is_finite()
        || total_entropy_change_j_k < -controls.entropy_tolerance_j_k
    {
        return Err(FiniteChannelError::EntropyDecrease {
            change_j_k: total_entropy_change_j_k,
            tolerance_j_k: controls.entropy_tolerance_j_k,
        });
    }
    check(&mut cancelled)?;
    let (parcels, stage) = gas_candidate.into_parts();
    Ok(FiniteChannelNetworkCandidate {
        parcels,
        walls: next_walls,
        step: FiniteChannelNetworkStep {
            duration_s,
            transport: stage.transport,
            channels: stage.channels,
            channel_walls,
            exchange: stage.exchange,
            walls: reports,
            heat_into_walls_j: wall_heat.to_f64(),
            wall_energy_change_j: wall_energy.to_f64(),
            energy_residual_j,
            gas_entropy_change_j_k: gas_change.entropy_change_j_k,
            wall_entropy_change_j_k: wall_entropy.to_f64(),
            total_entropy_change_j_k,
        },
    })
}

fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), FiniteChannelError> {
    if cancelled() {
        Err(FiniteChannelError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
