//! True-vacuum channel blowdown coupled to finite original material walls.
//!
//! The native vacuum-channel transaction owns every molecular beam, reflected
//! inventory, escaped species and gas reconstruction. Its accommodation heat is
//! pooled by physical wall identity before the same material storage inversion
//! used by closed channels. Actual wall storage replaces the prescribed-bath
//! account in the final source-relative energy balance.
//!
//! This is a first-order, fixed-geometry stage; callers refine time and resample
//! the evolving gas and wall temperatures. The wall has reversible isochoric
//! sensible storage, not an invented constant capacity. Vacuum has no assigned
//! temperature or entropy. Reported wall entropy is a material storage change,
//! not total entropy production of the unmodeled escaped gas.

/// Adaptive gas and finite-wall feedback with complete-window conservation.
pub mod trajectory;

use super::super::vacuum::{self as native, VacuumChannelNetworkLink, VacuumChannelNetworkStep};
use super::{
    FiniteChannelError, FiniteChannelWallStep, FixedVolumeChannelWall, PooledWallHeat, check,
    conduction::{self, ChannelWallConduction, ChannelWallConductionStep},
};
use crate::thermochemical_gas::energy::{
    GasEnergySolve,
    parcel::{
        ClosedGasParcel,
        transfer::{
            snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
            species::{
                effusion::channel::DiffuseMolecularChannel,
                network::transport::{GasTransportError, GasTransportPaths, GasTransportSolve},
            },
        },
    },
};
use fs_math::dd::Dd;
use std::collections::BTreeMap;

/// One physical channel array exhausting through a finite material wall.
#[derive(Clone, Copy, Debug)]
pub struct FiniteVacuumChannelLink<'a> {
    pub id: u64,
    pub parcel: u64,
    pub channel: &'a DiffuseMolecularChannel,
    /// Total exposed mouth area [m²]; zero closes the outlet.
    pub area_m2: f64,
    /// Stable finite wall identity. Several outlets may share the same wall.
    pub wall: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteVacuumChannelStep {
    /// Original molecular transaction sampled at the initial wall temperatures.
    /// Its `audit.heat_into_walls_j` is replaced by ACTUAL finite storage below.
    pub gas: VacuumChannelNetworkStep,
    pub channel_walls: BTreeMap<u64, u64>,
    pub walls: Vec<FiniteChannelWallStep>,
    /// Conducted heat is internal to these same original finite wall inventories.
    pub conduction: Vec<ChannelWallConductionStep>,
    pub wall_energy_change_j: f64,
    pub wall_entropy_change_j_k: f64,
    /// Delta U_gas + E_escaped + Delta U_wall, in the native source-relative
    /// chart. Species/mass residuals remain in `gas.audit` and use their own
    /// fixed physical budgets. No entropy production is assigned to vacuum.
    pub energy_residual_j: f64,
}

/// Both finite gas and material endpoints, published only after every check.
#[derive(Clone, Debug)]
pub struct FiniteVacuumChannelCandidate<'l, 'p> {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    walls: BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    step: FiniteVacuumChannelStep,
}
impl<'l, 'p> FiniteVacuumChannelCandidate<'l, 'p> {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub fn walls(&self) -> &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>> {
        &self.walls
    }
    pub const fn step(&self) -> &FiniteVacuumChannelStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, ClosedGasParcel>,
        BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
        FiniteVacuumChannelStep,
    ) {
        (self.parcels, self.walls, self.step)
    }
}

/// Evolve native internal gas paths and finite-wall vacuum outlets together.
///
/// All outgoing and returning molecules spend the original shared inventories.
/// The source-relative native audit and each actual wall storage increment are
/// checked before either endpoint is returned. An oversize interval, source
/// refusal or cancellation changes no caller state. The accepted material wall
/// temperatures become the next call's kinetic temperatures, closing feedback.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'l, 'p>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    paths: GasTransportPaths<'_>,
    links: &[FiniteVacuumChannelLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: GasTransportSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<FiniteVacuumChannelCandidate<'l, 'p>, FiniteChannelError> {
    step_candidate_with_conduction(
        parcels, walls, paths, links, &[], duration_s, endpoints, controls, cancelled,
    )
}

/// Advance vacuum channel escape and sourced finite-wall conduction together.
///
/// All molecular and conductive heat uses the same initial wall snapshot and
/// pools before the one original material inversion per wall. Conductive
/// transfers cancel internally; escaped molecular energy retains the native
/// source-relative account. Vacuum still has no modeled entropy production.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate_with_conduction<'l, 'p>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    paths: GasTransportPaths<'_>,
    links: &[FiniteVacuumChannelLink<'_>],
    conduction: &[ChannelWallConduction<'_, '_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: GasTransportSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<FiniteVacuumChannelCandidate<'l, 'p>, FiniteChannelError> {
    check(&mut cancelled)?;
    if walls.len() > MAX_SNAPSHOT_PARCELS || links.len() > MAX_SNAPSHOT_TRANSFERS {
        return Err(FiniteChannelError::Invalid(
            "finite vacuum wall/channel work budget exceeded",
        ));
    }
    let mut ordered: Vec<_> = links.iter().collect();
    ordered.sort_by_key(|link| link.id);
    let mut channel_walls = BTreeMap::new();
    let mut pools = BTreeMap::<u64, PooledWallHeat>::new();
    let mut ranges = BTreeMap::<u64, [f64; 2]>::new();
    let mut held = Vec::with_capacity(ordered.len());
    for link in &ordered {
        check(&mut cancelled)?;
        if channel_walls.insert(link.id, link.wall).is_some() {
            return Err(FiniteChannelError::Invalid(
                "duplicate finite vacuum channel identity",
            ));
        }
        let wall = walls.get(&link.wall).ok_or(FiniteChannelError::Invalid(
            "finite vacuum channel references an absent material wall",
        ))?;
        let range = ranges.entry(link.wall).or_insert([wall.temperature_k(); 2]);
        if link.area_m2 > 0.0 && link.channel.geometry().length_m() > 0.0 {
            let gas = parcels
                .get(&link.parcel)
                .ok_or(FiniteChannelError::Invalid(
                    "finite vacuum channel references an absent gas",
                ))?;
            let temperature = gas.state().state().temperature;
            range[0] = range[0].min(temperature);
            range[1] = range[1].max(temperature);
        }
        pools.entry(link.wall).or_default();
        held.push(VacuumChannelNetworkLink {
            id: link.id,
            parcel: link.parcel,
            channel: link.channel,
            area_m2: link.area_m2,
            wall_temperature_k: wall.temperature_k(),
        });
    }
    let conduction_steps = conduction::pool(
        walls, conduction, duration_s, &mut pools, &mut ranges, &mut cancelled,
    )?;
    let gas = native::step_candidate(
        parcels,
        paths,
        &held,
        duration_s,
        endpoints,
        controls,
        &mut cancelled,
    )?;
    for outlet in &gas.step().outlets {
        check(&mut cancelled)?;
        let pool = pools
            .get_mut(&channel_walls[&outlet.id])
            .expect("admitted original finite wall identity");
        pool.molecular = pool.molecular + Dd::from_f64(outlet.wall.heat_into_wall_j);
    }
    let mut next_walls = walls.clone();
    let mut reports = Vec::with_capacity(pools.len());
    let mut wall_energy = Dd::ZERO;
    let mut wall_entropy = Dd::ZERO;
    // Subtract each ORIGINAL sampled heat before summation so a large
    // inter-wall exchange cannot hide a small storage mismatch by cancellation.
    let mut storage_correction = Dd::ZERO;
    let mut conductive_balance = Dd::ZERO;
    for (&id, &pool) in &pools {
        check(&mut cancelled)?;
        let (next, report, energy, entropy) = walls[&id].heated_candidate(
            id,
            pool,
            ranges[&id],
            controls.exchange.energy_tolerance_j,
        )?;
        wall_energy = wall_energy + energy;
        wall_entropy = wall_entropy + entropy;
        storage_correction = storage_correction + (energy - pool.total());
        conductive_balance = conductive_balance + pool.conduction;
        next_walls.insert(id, next);
        reports.push(report);
    }
    for link in &ordered {
        check(&mut cancelled)?;
        if link.area_m2 == 0.0 || link.channel.geometry().length_m() == 0.0 {
            continue;
        }
        link.channel
            .admit_wall_for_parcel(
                &gas.parcels()[&link.parcel],
                next_walls[&link.wall].temperature_k(),
                &mut cancelled,
            )
            .map_err(|error| {
                GasTransportError::from(super::super::super::face_error(link.id, error))
            })?;
    }
    conduction::admit_endpoints(&next_walls, conduction, &mut cancelled)?;
    let energy_residual_j =
        (Dd::from_f64(gas.step().audit.energy_residual_j) + storage_correction
            + conductive_balance).to_f64();
    if !energy_residual_j.is_finite()
        || energy_residual_j.abs() > controls.exchange.energy_tolerance_j
    {
        return Err(FiniteChannelError::TotalEnergyBalance {
            residual_j: energy_residual_j,
        });
    }
    if !wall_energy.to_f64().is_finite() || !wall_entropy.to_f64().is_finite() {
        return Err(FiniteChannelError::Invalid(
            "finite vacuum wall account is unrepresentable",
        ));
    }
    check(&mut cancelled)?;
    let (parcels, gas) = gas.into_parts();
    Ok(FiniteVacuumChannelCandidate {
        parcels,
        walls: next_walls,
        step: FiniteVacuumChannelStep {
            gas,
            channel_walls,
            walls: reports,
            conduction: conduction_steps,
            wall_energy_change_j: wall_energy.to_f64(),
            wall_entropy_change_j_k: wall_entropy.to_f64(),
            energy_residual_j,
        },
    })
}

#[cfg(test)]
mod tests;
