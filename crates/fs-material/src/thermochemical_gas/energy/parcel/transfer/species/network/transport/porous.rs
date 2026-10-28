//! Pressure-driven porous leakage on the original finite-gas snapshot.
//!
//! The material barrier owns the compressible Darcy/Klinkenberg closure. This
//! adapter supplies the actual upstream viscosity, EOS and composition, then
//! lowers its mass flux to the SAME conservative inventory owner as the nozzle
//! and selective transport paths. There is no second energy or species update.
//!
//! Pores are quasi-steady, dry, rigid and have negligible gas storage. The
//! isothermal slab uses the original upstream temperature and frozen mixture;
//! unequal endpoint compositions therefore use an explicit upwind advective
//! approximation. Molecular mixing needs its own authored diffusion path.
//! This is neither a resolved pore-filling transient nor a nonisothermal porous
//! heat exchanger. Refine the physical timestep and respect the barrier's
//! calibrated pressure, temperature and superficial-speed support.

use std::collections::BTreeMap;

use super::{ClosedGasParcel, MAX_SNAPSHOT_TRANSFERS};
use crate::MaterialError;
use crate::gas::GasState;
use crate::porous::permeability::{DarcyGasBarrier, DarcyGasFlow};
use crate::thermochemical_gas::energy::parcel::transfer::snapshot::GasSnapshotTransfer;

/// An explicitly oriented, source-bearing porous connection.
#[derive(Clone, Copy, Debug)]
pub struct PorousNetworkLink<'a> {
    /// Reporting identity; never selects an inventory-update order.
    pub id: u64,
    /// Left/right finite volumes. Flow direction follows their original pressure.
    pub parcels: [u64; 2],
    /// Compiled intrinsic permeability, layer geometry and material support.
    pub barrier: &'a DarcyGasBarrier,
    /// Exposed area [m²]; zero closes the connection exactly.
    pub area_m2: f64,
    /// Explicit isothermal-face admission [K]; no temperature is averaged away.
    pub maximum_temperature_difference_k: f64,
}

/// Original material flow and its integrated mass, before any mixing.
#[derive(Clone, Debug, PartialEq)]
pub struct PorousNetworkLinkStep {
    pub id: u64,
    pub parcels: [u64; 2],
    pub donor: u64,
    pub receiver: u64,
    /// Coherent viscosity, density and gas constant actually used by the law.
    pub upstream_state: GasState,
    /// Signed rates retain the connection's left/right orientation.
    pub flow: DarcyGasFlow,
    /// Nonnegative mass withdrawn from the original upstream inventory [kg].
    pub transferred_mass_kg: f64,
}

/// A refusal never mutates any accepted parcel.
#[derive(Clone, Debug, PartialEq)]
pub enum PorousNetworkError {
    Invalid(&'static str),
    Material { id: u64, error: MaterialError },
    /// The combined endpoint overshot an initially active pressure ordering.
    /// Shorten the interval; do not clip an individual edge's mass.
    EquilibriumCrossed { id: u64 },
    Cancelled,
}

impl core::fmt::Display for PorousNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material porous gas network: {self:?}")
    }
}
impl std::error::Error for PorousNetworkError {}

use PorousNetworkError as Error;

fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}

/// Sample constitutive rates only; the enclosing transport owner spends mass.
pub(super) fn sample_transfers(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    links: &[PorousNetworkLink<'_>],
    duration_s: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(Vec<GasSnapshotTransfer>, Vec<PorousNetworkLinkStep>), Error> {
    check(&mut cancelled)?;
    if !duration_s.is_finite() || duration_s <= 0.0 || links.len() > MAX_SNAPSHOT_TRANSFERS {
        return Err(Error::Invalid("positive duration and bounded porous connection count required"));
    }
    let mut ordered = BTreeMap::new();
    for link in links {
        check(&mut cancelled)?;
        if link.parcels[0] == link.parcels[1]
            || link.parcels.iter().any(|id| !parcels.contains_key(id))
            || !link.area_m2.is_finite() || link.area_m2 < 0.0
            || !link.maximum_temperature_difference_k.is_finite()
            || link.maximum_temperature_difference_k < 0.0
            || ordered.insert(link.id, link).is_some()
        {
            return Err(Error::Invalid(
                "unique porous IDs, distinct existing endpoints, nonnegative area and temperature allowance required",
            ));
        }
    }
    let mut transfers = Vec::with_capacity(links.len());
    let mut reports = Vec::with_capacity(links.len());
    for (&id, link) in &ordered {
        check(&mut cancelled)?;
        let [left, right] = link.parcels;
        let a = parcels[&left].state().state();
        let b = parcels[&right].state().state();
        let active = link.area_m2 > 0.0
            && link.barrier.geometric_resistance_per_m().is_some()
            && a.pressure != b.pressure;
        if active && (a.temperature - b.temperature).abs() > link.maximum_temperature_difference_k {
            return Err(Error::Invalid("porous gas endpoints exceed the declared isothermal-face allowance"));
        }
        let (donor, receiver) = if a.pressure > b.pressure
            || (a.pressure == b.pressure && left < right)
        {
            (left, right)
        } else {
            (right, left)
        };
        let upstream = &parcels[&donor];
        let flow = link.barrier.mass_flow_for_model(
            upstream.state().state(), upstream.model(), a.pressure, b.pressure, link.area_m2,
        ).map_err(|error| Error::Material { id, error })?;
        let mass_kg = flow.mass_rate_left_to_right_kg_s.abs() * duration_s;
        if !mass_kg.is_finite()
            || (flow.mass_rate_left_to_right_kg_s != 0.0 && mass_kg == 0.0)
        {
            return Err(Error::Invalid("integrated porous gas mass is unrepresentable"));
        }
        transfers.push(GasSnapshotTransfer { donor, receiver, mass_kg });
        reports.push(PorousNetworkLinkStep {
            id, parcels: link.parcels, donor, receiver,
            upstream_state: *upstream.state().state(), flow, transferred_mass_kg: mass_kg,
        });
    }
    check(&mut cancelled)?;
    Ok((transfers, reports))
}

/// Check pressure ordering on the FINAL shared transport/heat/work endpoint.
pub(super) fn admit_pressures(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    reports: &[PorousNetworkLinkStep],
    mut cancelled: impl FnMut() -> bool,
) -> Result<(), Error> {
    for report in reports {
        check(&mut cancelled)?;
        if report.transferred_mass_kg > 0.0
            && parcels[&report.donor].state().state().pressure
                < parcels[&report.receiver].state().state().pressure
        {
            return Err(Error::EquilibriumCrossed { id: report.id });
        }
    }
    check(&mut cancelled)
}
