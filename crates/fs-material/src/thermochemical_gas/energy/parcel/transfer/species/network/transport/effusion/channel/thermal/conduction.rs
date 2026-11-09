//! Sourced conduction between the original finite molecular-channel walls.
//!
//! A prepared native interface binds endpoint IDs, original caloric/source laws
//! and fixed transfer geometry. The channel owner supplies the actual matching
//! wall temperatures; the original conduction kernel supplies its heat flow.
//! Every opposite pair enters the same wall heat pools as molecular scattering
//! before ONE original caloric inversion. No duplicate capacity, thermostat or
//! artificial deposition geometry is created to adapt the state owners.
//!
//! Cellwise links use each endpoint's k(T) and an explicit series/contact
//! resistance. Integrated links retain the original source-free, one-dimensional
//! layer solve and both original conductivity charts. Each finite wall remains
//! uniform in temperature; this does not resolve an axial temperature profile
//! inside an individual bore. Link IDs identify distinct physical heat paths.

use super::{FiniteChannelError, FixedVolumeChannelWall, PooledWallHeat, check};
use crate::thermal::{
    conduction::{LayerConductionConfig, LayerConductionResponse},
    transport::volumetric::timeline::{FixedVolumeConductionLink, FixedVolumeLayerConduction},
};
use crate::thermochemical_gas::energy::parcel::transfer::snapshot::MAX_SNAPSHOT_TRANSFERS;
use fs_blake3::ContentHash;
use fs_math::dd::Dd;
use fs_thermal::enthalpy::DepositionError;
use std::collections::BTreeMap;

/// Original material interfaces evaluated on the matching finite wall states.
///
/// Construct the underlying link with wall IDs as `CellId`s and the SAME wall
/// laws. The original constructor declares aligned scalar conductivity and
/// retains distances, area and area-specific contact resistance. Converting to
/// integrated layers prepares their source curves once, outside time stepping.
#[derive(Clone, Debug)]
pub enum ChannelWallConduction<'l, 'p> {
    Cellwise(FixedVolumeConductionLink<'l, 'p>),
    Integrated {
        link: FixedVolumeLayerConduction<'l, 'p>,
        controls: LayerConductionConfig,
    },
}

impl<'l, 'p> ChannelWallConduction<'l, 'p> {
    pub const fn binding(&self) -> &FixedVolumeConductionLink<'l, 'p> {
        match self {
            Self::Cellwise(link) => link,
            Self::Integrated { link, .. } => link.binding(),
        }
    }
    pub const fn id(&self) -> u64 {
        self.binding().id()
    }
}

/// Accepted explicit heat transfer; positive power/heat leaves a and enters b.
#[derive(Clone, Debug, PartialEq)]
pub struct ChannelWallConductionStep {
    pub id: u64,
    pub walls: [u64; 2],
    pub caloric_identities: [ContentHash; 2],
    pub temperatures_k: [f64; 2],
    pub conductance_w_k: f64,
    pub power_a_to_b_w: f64,
    pub heat_a_to_b_j: f64,
    /// Native integrated-layer interfaces, derivatives and residuals when used.
    /// A cellwise resistance has no inferred internal temperature profile.
    pub layers: Option<LayerConductionResponse>,
}

/// One original snapshot and one deterministic, opposite heat pair per link.
pub(super) fn pool(
    walls: &BTreeMap<u64, FixedVolumeChannelWall<'_, '_>>,
    links: &[ChannelWallConduction<'_, '_>],
    duration_s: f64,
    pools: &mut BTreeMap<u64, PooledWallHeat>,
    ranges: &mut BTreeMap<u64, [f64; 2]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<ChannelWallConductionStep>, FiniteChannelError> {
    if links.len() > MAX_SNAPSHOT_TRANSFERS {
        return Err(FiniteChannelError::Invalid(
            "finite wall conduction work budget exceeded",
        ));
    }
    if !links.is_empty() && (!duration_s.is_finite() || duration_s <= 0.0) {
        return Err(FiniteChannelError::Invalid(
            "finite wall conduction duration must be finite and positive",
        ));
    }
    let mut ordered: Vec<_> = links.iter().collect();
    ordered.sort_by_key(|link| link.id());
    if ordered.windows(2).any(|pair| pair[0].id() == pair[1].id()) {
        return Err(FiniteChannelError::Invalid(
            "duplicate finite wall conduction identity",
        ));
    }
    let mut steps = Vec::with_capacity(ordered.len());
    for link in ordered {
        check(cancelled)?;
        let binding = link.binding();
        let ids = binding.cell_ids().map(|id| id.0);
        let laws = binding.endpoint_laws();
        let caloric_identities = laws.map(|law| law.identity());
        let mut temperatures = [0.0; 2];
        for i in 0..2 {
            let wall = walls.get(&ids[i]).ok_or(FiniteChannelError::Invalid(
                "sourced conduction references an absent finite wall",
            ))?;
            if wall.law().identity() != caloric_identities[i] {
                return Err(FiniteChannelError::Invalid(
                    "sourced conduction requires the original finite wall caloric law and mass basis",
                ));
            }
            temperatures[i] = wall.temperature_k();
        }
        // A closed cellwise path keeps its binding but queries no inactive k(T).
        let (conductance, power, layers) = if binding.area_m2() == 0.0 {
            (0.0, 0.0, None)
        } else {
            match link {
                ChannelWallConduction::Cellwise(native) => {
                    let conductance = native
                        .conductance_at_temperatures(temperatures)
                        .map_err(|error| failure(binding.id(), error))?;
                    let power = (Dd::from_f64(conductance)
                        * (Dd::from_f64(temperatures[0]) - Dd::from_f64(temperatures[1])))
                    .to_f64();
                    (conductance, power, None)
                }
                ChannelWallConduction::Integrated { link, controls } => {
                    let response = link
                        .evaluate_at_temperatures(temperatures, *controls, || {
                            if cancelled() { Err(DepositionError::Cancelled) } else { Ok(()) }
                        })
                        .map_err(|error| failure(binding.id(), error))?;
                    (response.conductance_w_k, response.power_a_to_b_w, Some(response))
                }
            }
        };
        let heat = (Dd::from_f64(power) * Dd::from_f64(duration_s)).to_f64();
        if !conductance.is_finite() || conductance < 0.0
            || !power.is_finite() || !heat.is_finite()
            || (conductance > 0.0 && temperatures[0] != temperatures[1] && power == 0.0)
            || (power != 0.0 && heat == 0.0)
        {
            return Err(FiniteChannelError::Invalid(
                "finite wall conduction heat is not representable",
            ));
        }
        if binding.area_m2() > 0.0 {
            for i in 0..2 {
                let range = ranges.entry(ids[i]).or_insert([temperatures[i]; 2]);
                range[0] = range[0].min(temperatures[1 - i]);
                range[1] = range[1].max(temperatures[1 - i]);
                let pool = pools.entry(ids[i]).or_default();
                let signed_heat = if i == 0 { -heat } else { heat };
                pool.conduction = pool.conduction + Dd::from_f64(signed_heat);
            }
        }
        steps.push(ChannelWallConductionStep {
            id: binding.id(),
            walls: ids,
            caloric_identities,
            temperatures_k: temperatures,
            conductance_w_k: conductance,
            power_a_to_b_w: power,
            heat_a_to_b_j: heat,
            layers,
        });
    }
    Ok(steps)
}

fn failure(id: u64, error: DepositionError) -> FiniteChannelError {
    if matches!(error, DepositionError::Cancelled) {
        FiniteChannelError::Cancelled
    } else {
        FiniteChannelError::Conduction { id, error }
    }
}

/// Integrated layers require each original source across the whole final span,
/// even when molecular heating remains inside each wall's own caloric domain.
/// This admission does not run another nonlinear layer solve.
pub(super) fn admit_endpoints(
    walls: &BTreeMap<u64, FixedVolumeChannelWall<'_, '_>>,
    links: &[ChannelWallConduction<'_, '_>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), FiniteChannelError> {
    let mut ordered: Vec<_> = links.iter().collect();
    ordered.sort_by_key(|link| link.id());
    for link in ordered {
        check(cancelled)?;
        if let ChannelWallConduction::Integrated { link, .. } = link {
            let temperatures = link.binding().cell_ids().map(|id| walls[&id.0].temperature_k());
            let range = [temperatures[0].min(temperatures[1]), temperatures[0].max(temperatures[1])];
            for source in link.source_paths() {
                let span = source.temperature_span_k();
                if range[0] < span[0] || range[1] > span[1] {
                    return Err(failure(link.id(), DepositionError::MaterialDomain));
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
