//! Open molecular escape shares the original finite-gas conservation owner.
//!
//! Internal paths and vacuum outlets sample the SAME original gases. Every
//! outgoing species spends the same gross allowance, so arrivals cannot fund
//! a leak in the same step. The source-relative boundary coordinate debits
//! molecular energy rather than bulk enthalpy without subtracting rounded
//! formation energies. Each finite gas is reconstructed exactly once.
//!
//! Escaped mass and energy leave the modeled system. No receiver inventory,
//! vacuum temperature, closed-system entropy condition or pump work is invented.
//! Fixed-volume well-mixed gases and the original aperture qualifications apply.
//! Time integration is an explicit first-order candidate; callers refine time
//! and accept a complete candidate atomically. Complete gas depletion refuses.

use super::{check, face_error};
use crate::thermochemical_gas::energy::parcel::transfer::{
    snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
    species::{
        effusion::{EffusiveInterface, positive_product_ratio, vacuum::VacuumEffusiveFlux},
        network::{
            SpeciesNetworkStep,
            boundary::{BoundarySpeciesTransfer, SpeciesNetworkBoundary},
            exchange_candidate_with_boundaries,
            transport::{
                self, GasTransportError, GasTransportPaths, GasTransportRates, GasTransportSolve,
            },
        },
    },
};
use crate::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use std::collections::BTreeMap;

/// One separately exposed physical outlet to an ideal absorbing vacuum.
#[derive(Clone, Copy, Debug)]
pub struct VacuumEffusionLink<'a> {
    /// Unique outlet identity, independent of internal path IDs.
    pub id: u64,
    pub parcel: u64,
    pub interface: &'a EffusiveInterface,
    /// Total exposed area of qualified thin apertures [m2]; zero closes the outlet.
    pub area_m2: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct VacuumEffusionOutletStep {
    pub id: u64,
    pub parcel: u64,
    pub flux: VacuumEffusiveFlux,
}

/// Original rates and actual open-system species/energy accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct VacuumEffusionStep {
    pub duration_s: f64,
    /// Original internal path rates; this open subsystem has no closed entropy claim.
    pub transport: GasTransportRates,
    pub outlets: Vec<VacuumEffusionOutletStep>,
    /// `boundaries` gives actual escaped species masses and complete energies
    /// with negative sign into gas. Source-relative residuals remain independent
    /// of arbitrary formation-energy reference magnitudes.
    pub exchange: SpeciesNetworkStep,
}

impl VacuumEffusionStep {
    /// Re-admit original active path sources after enclosing heat/work coupling.
    /// This samples no additional flow and spends no inventory a second time.
    pub fn admit_endpoint_driving_forces(
        &self,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), GasTransportError> {
        self.transport
            .admit_endpoint_driving_forces(parcels, &mut cancelled)?;
        for row in &self.outlets {
            check(&mut cancelled)?;
            let gas = parcels.get(&row.parcel).ok_or(GasTransportError::Invalid(
                "vacuum endpoint admission lost its original gas identity",
            ))?;
            if !row.flux.species.is_empty() {
                row.flux
                    .interface
                    .admit_parcel(gas, &mut cancelled)
                    .map_err(|error| face_error(row.id, error))?;
            }
        }
        check(&mut cancelled)?;
        Ok(())
    }
}

/// All finite endpoints or none. Inputs are never mutated by an attempted step.
#[derive(Clone, Debug, PartialEq)]
pub struct VacuumEffusionCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    step: VacuumEffusionStep,
}
impl VacuumEffusionCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn step(&self) -> &VacuumEffusionStep {
        &self.step
    }
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, VacuumEffusionStep) {
        (self.parcels, self.step)
    }
}

/// Advance native internal paths and vacuum outlets in one finite transaction.
///
/// Outgoing species carry `h_i(T)-R*T/(2*M_i)` rather than nozzle enthalpy.
/// In the boundary's original-partial-enthalpy chart, this is `dm_into=-dm_out`
/// and `Q_relative=dm_out*R*T/(2*M_i)`. The shared caloric owner then obtains
/// the cooling target `-dm_out*R*T/(2*M_i)` after inventory-reference changes
/// cancel. No constant heat capacity or second gas equation of state is used.
///
/// `controls.exchange` governs the shared gross withdrawals and absolute mass
/// and energy budgets. `controls.entropy_tolerance_j_k` is retained for native
/// path input admission but cannot constrain the gas-only entropy of an open
/// vacuum system. If no outlet actually emits, the original closed-network
/// entropy admission remains in force. Outlets and internal paths must describe
/// distinct physical apertures. Refusals/cancellation publish no partial state.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    paths: GasTransportPaths<'_>,
    vacuum: &[VacuumEffusionLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: GasTransportSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<VacuumEffusionCandidate, GasTransportError> {
    check(&mut cancelled)?;
    let connection_count = [
        paths.nozzles.len(),
        paths.diffusion.len(),
        paths.porous.len(),
        paths.membranes.len(),
        paths.dusty_gas.len(),
        paths.effusion.len(),
        vacuum.len(),
    ]
    .into_iter()
    .try_fold(0usize, |total, n| total.checked_add(n));
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || parcels.is_empty()
        || parcels.len() > MAX_SNAPSHOT_PARCELS
        || connection_count.is_none_or(|n| n > MAX_SNAPSHOT_TRANSFERS)
    {
        return Err(GasTransportError::Invalid(
            "vacuum transport needs a bounded nonempty graph and positive duration",
        ));
    }
    let mut ordered = BTreeMap::new();
    for link in vacuum {
        check(&mut cancelled)?;
        if !parcels.contains_key(&link.parcel)
            || !link.area_m2.is_finite()
            || link.area_m2 < 0.0
            || ordered.insert(link.id, link).is_some()
        {
            return Err(GasTransportError::Invalid(
                "vacuum outlets need unique IDs, existing gases and nonnegative areas",
            ));
        }
    }
    let sampled = transport::sample_transfers(
        parcels,
        paths,
        duration_s,
        endpoints,
        controls,
        &mut cancelled,
    )?;
    let mut work = sampled.transfers.len() + sampled.heat.len();
    let mut outlets = Vec::with_capacity(vacuum.len());
    let mut species_rows = Vec::with_capacity(vacuum.len());
    for (id, link) in ordered {
        check(&mut cancelled)?;
        let flux = parcels[&link.parcel]
            .effusive_flux_to_vacuum(link.interface, link.area_m2, &mut cancelled)
            .map_err(|error| face_error(id, error))?;
        let mut rows = Vec::with_capacity(flux.species.len());
        for species in &flux.species {
            check(&mut cancelled)?;
            let beam = &species.outward;
            if beam.mass_rate_kg_s == 0.0 {
                continue;
            }
            // The correction uses the actual integrated mass spent below,
            // never a separately rounded energy-rate-times-duration product.
            let mass = positive_product_ratio(&[beam.mass_rate_kg_s, duration_s], &[])
                .map_err(|error| face_error(id, error))?;
            let specific =
                beam.additional_specific_energy_j_kg
                    .ok_or(GasTransportError::Invalid(
                        "vacuum beam lost its molecular energy convention",
                    ))?;
            let relative = positive_product_ratio(&[mass, -specific], &[])
                .map_err(|error| face_error(id, error))?;
            work += 1 + usize::from(rows.is_empty());
            if work > MAX_SNAPSHOT_TRANSFERS {
                return Err(GasTransportError::Invalid(
                    "vacuum and internal species streams exceed the shared work bound",
                ));
            }
            rows.push(BoundarySpeciesTransfer {
                species: species.species.clone(),
                incoming_mass_kg: 0.0,
                outgoing_mass_kg: mass,
                net_mass_into_gas_kg: -mass,
                relative_energy_into_gas_j: relative,
            });
        }
        outlets.push(VacuumEffusionOutletStep {
            id,
            parcel: link.parcel,
            flux,
        });
        species_rows.push(rows);
    }
    let boundaries: Vec<_> = outlets
        .iter()
        .zip(&species_rows)
        .filter(|(_, rows)| !rows.is_empty())
        .map(|(outlet, rows)| SpeciesNetworkBoundary {
            id: outlet.id,
            parcel: outlet.parcel,
            // Do not introduce inactive interface catalog species into gas.
            catalog: parcels[&outlet.parcel].model(),
            species: rows,
        })
        .collect();
    let candidate = exchange_candidate_with_boundaries(
        parcels,
        &sampled.transfers,
        &sampled.heat,
        &boundaries,
        endpoints,
        controls.exchange,
        &mut cancelled,
    )?;
    if boundaries.is_empty() {
        transport::audit_closed_window(parcels, candidate.parcels(), controls, &mut cancelled)?;
    }
    let (parcels, exchange) = candidate.into_parts();
    let step = VacuumEffusionStep {
        duration_s,
        transport: sampled.rates,
        outlets,
        exchange,
    };
    step.admit_endpoint_driving_forces(&parcels, &mut cancelled)?;
    check(&mut cancelled)?;
    Ok(VacuumEffusionCandidate { parcels, step })
}

#[cfg(test)]
mod tests;
