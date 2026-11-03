//! Prescribed molecular sources share the original finite-gas boundary solve.
//!
//! Both beams remain distinct even at zero net mass. The external reservoir is
//! a prescribed thermochemical source, not a finite endpoint silently depleted
//! here. Its owner accepts the opposite exported species and energy accounts.
//! The scalar molecular law uses a common rest frame; moving-source kinematics
//! and exterior-wall forces belong to the composing mechanical owner.

use super::{check, face_error, vacuum};
use crate::thermochemical_gas::energy::parcel::transfer::{
    snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
    species::{
        Sum,
        effusion::{EffusiveFlux, EffusiveInterface, positive_product_ratio},
        network::{
            SpeciesNetworkStep,
            boundary::{BoundarySpeciesTransfer, SpeciesNetworkBoundary},
            exchange_candidate_with_boundaries,
            transport::{
                self, GasTransportError, GasTransportPaths, GasTransportRates, GasTransportSolve,
            },
        },
        source_error,
    },
};
use crate::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_qty::Temperature;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use std::collections::{BTreeMap, BTreeSet};
use vacuum::{VacuumEffusionLink, VacuumEffusionOutletStep};

/// A finite gas (left) exchanging with a prescribed molecular reservoir (right).
#[derive(Clone, Copy, Debug)]
pub struct EffusionReservoirLink<'a> {
    /// Unique across reservoir and vacuum ports in this open transaction.
    pub id: u64,
    pub parcel: u64,
    /// Prescribed source state. Its volume does not impose a finite donor budget.
    pub reservoir: &'a ClosedGasParcel,
    pub interface: &'a EffusiveInterface,
    /// Separately exposed qualified area [m2]; zero is inactive.
    pub area_m2: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EffusionReservoirStep {
    pub id: u64,
    pub parcel: u64,
    /// Original source-bearing two-beam sample; left is the modeled finite gas.
    pub flux: EffusiveFlux,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ReservoirEffusionStep {
    pub duration_s: f64,
    pub transport: GasTransportRates,
    pub vacuum_outlets: Vec<VacuumEffusionOutletStep>,
    pub reservoirs: Vec<EffusionReservoirStep>,
    /// Actual gross species transfers and stable net caloric accounts.
    pub exchange: SpeciesNetworkStep,
}
impl ReservoirEffusionStep {
    /// Re-admit actual endpoints after enclosing work/heat without resampling flow.
    pub fn admit_endpoint_driving_forces(
        &self,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), GasTransportError> {
        self.transport
            .admit_endpoint_driving_forces(parcels, &mut cancelled)?;
        vacuum::admit_outlets(&self.vacuum_outlets, parcels, &mut cancelled)?;
        for row in &self.reservoirs {
            check(&mut cancelled)?;
            if row.flux.species.is_empty() {
                continue;
            }
            let gas = parcels.get(&row.parcel).ok_or(GasTransportError::Invalid(
                "reservoir endpoint lost its original finite gas",
            ))?;
            let interface = &row.flux.interface;
            interface
                .admit_parcel(gas, &mut cancelled)
                .map_err(|e| face_error(row.id, e))?;
            // A selected zero transmission still qualifies both states, even
            // when the species exists only on the prescribed source side.
            for species in &row.flux.species {
                check(&mut cancelled)?;
                let index = interface
                    .transmissions()
                    .binary_search_by(|t| t.species().cmp(&species.species))
                    .map_err(|_| {
                        GasTransportError::Invalid("reservoir lost original species transmission")
                    })?;
                let state = gas.state().state();
                interface.transmissions()[index]
                    .admit(state.temperature, state.pressure)
                    .map_err(|e| face_error(row.id, e))?;
            }
        }
        check(&mut cancelled)?;
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ReservoirEffusionCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    step: ReservoirEffusionStep,
}
impl ReservoirEffusionCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn step(&self) -> &ReservoirEffusionStep {
        &self.step
    }
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, ReservoirEffusionStep) {
        (self.parcels, self.step)
    }
}

fn signed_product(id: u64, a: f64, b: f64) -> Result<f64, GasTransportError> {
    if !a.is_finite() || !b.is_finite() {
        return Err(GasTransportError::Invalid(
            "nonfinite reservoir energy product",
        ));
    }
    if a == 0.0 || b == 0.0 {
        return Ok(0.0);
    }
    let magnitude =
        positive_product_ratio(&[a.abs(), b.abs()], &[]).map_err(|e| face_error(id, e))?;
    Ok(magnitude.copysign(a.signum() * b.signum()))
}

/// Advance internal paths, vacuum escape and prescribed reservoirs simultaneously.
///
/// Gross outgoing species spend the ORIGINAL inventory, including when an inlet
/// replaces the same species. All boundary rows enter one original caloric
/// reconstruction. Incoming molecular energy is already complete; add no second
/// thermalization term. A source-relative enthalpy difference preserves heat
/// beside arbitrary formation datums and at zero net mass. Open-system entropy
/// is not constrained as though the external reservoir were a closed endpoint.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    paths: GasTransportPaths<'_>,
    vacuum_links: &[VacuumEffusionLink<'_>],
    reservoir_links: &[EffusionReservoirLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: GasTransportSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ReservoirEffusionCandidate, GasTransportError> {
    check(&mut cancelled)?;
    let count = [
        paths.nozzles.len(),
        paths.diffusion.len(),
        paths.porous.len(),
        paths.membranes.len(),
        paths.dusty_gas.len(),
        paths.effusion.len(),
        vacuum_links.len(),
        reservoir_links.len(),
    ]
    .into_iter()
    .try_fold(0usize, |n, k| n.checked_add(k));
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || parcels.is_empty()
        || parcels.len() > MAX_SNAPSHOT_PARCELS
        || count.is_none_or(|n| n > MAX_SNAPSHOT_TRANSFERS)
    {
        return Err(GasTransportError::Invalid(
            "reservoir transport requires bounded original gases, ports and positive time",
        ));
    }
    let mut ids = BTreeSet::new();
    for link in vacuum_links {
        check(&mut cancelled)?;
        if !ids.insert(link.id) {
            return Err(GasTransportError::Invalid(
                "duplicate external molecular port",
            ));
        }
    }
    let mut ordered = BTreeMap::new();
    for link in reservoir_links {
        check(&mut cancelled)?;
        if !parcels.contains_key(&link.parcel)
            || !link.area_m2.is_finite()
            || link.area_m2 < 0.0
            || !ids.insert(link.id)
        {
            return Err(GasTransportError::Invalid(
                "reservoir ports need unique IDs, existing gases and nonnegative areas",
            ));
        }
        ordered.insert(link.id, link);
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
    let (vacuum_outlets, vacuum_rows) =
        vacuum::sample_outlets(parcels, vacuum_links, duration_s, &mut work, &mut cancelled)?;
    let mut reservoirs = Vec::with_capacity(reservoir_links.len());
    let mut reservoir_rows = Vec::with_capacity(reservoir_links.len());
    for (id, link) in ordered {
        check(&mut cancelled)?;
        let gas = &parcels[&link.parcel];
        let flux = gas
            .effusive_flux_to(link.reservoir, link.interface, link.area_m2, &mut cancelled)
            .map_err(|e| face_error(id, e))?;
        let (tl, tr) = (
            gas.state().state().temperature,
            link.reservoir.state().state().temperature,
        );
        let mut rows = Vec::with_capacity(flux.species.len());
        for species in &flux.species {
            check(&mut cancelled)?;
            let (incoming, outgoing) = (species.right_to_left, species.left_to_right);
            if incoming.mass_rate_kg_s == 0.0 && outgoing.mass_rate_kg_s == 0.0 {
                continue;
            }
            let mass = |rate| {
                if rate == 0.0 {
                    return Ok(0.0);
                }
                positive_product_ratio(&[rate, duration_s], &[]).map_err(|e| face_error(id, e))
            };
            let mi = mass(incoming.mass_rate_kg_s)?;
            let mo = mass(outgoing.mass_rate_kg_s)?;
            let net = mi - mo;
            let original = link
                .interface
                .gas_catalog()
                .thermodynamics()
                .components()
                .iter()
                .find(|row| row.species() == &species.species)
                .ok_or(GasTransportError::Invalid(
                    "reservoir beam lost its original caloric source",
                ))?;
            let molar_mass = original.molar_mass().value();
            let a_left = -positive_product_ratio(&[0.5, R, tl], &[molar_mass])
                .map_err(|e| face_error(id, e))?;
            let mut relative = Sum::default();
            relative
                .add(signed_product(id, net, a_left)?)
                .map_err(|e| face_error(id, e))?;
            if mi > 0.0 && tl != tr {
                let dh = original
                    .molar_enthalpy_change(Temperature::new(tl), Temperature::new(tr))
                    .map_err(|e| face_error(id, source_error(e)))?
                    .value()
                    / molar_mass;
                let correction = positive_product_ratio(&[0.5, R, (tr - tl).abs()], &[molar_mass])
                    .map_err(|e| face_error(id, e))?
                    .copysign(tr - tl);
                relative
                    .add(signed_product(id, mi, dh - correction)?)
                    .map_err(|e| face_error(id, e))?;
            }
            work += 1 + usize::from(rows.is_empty());
            if work > MAX_SNAPSHOT_TRANSFERS {
                return Err(GasTransportError::Invalid(
                    "combined molecular species streams exceed shared work bound",
                ));
            }
            rows.push(BoundarySpeciesTransfer {
                species: species.species.clone(),
                incoming_mass_kg: mi,
                outgoing_mass_kg: mo,
                net_mass_into_gas_kg: net,
                relative_energy_into_gas_j: relative.value,
            });
        }
        reservoirs.push(EffusionReservoirStep {
            id,
            parcel: link.parcel,
            flux,
        });
        reservoir_rows.push(rows);
    }
    let boundaries: Vec<_> = vacuum_outlets
        .iter()
        .zip(&vacuum_rows)
        .filter(|(_, rows)| !rows.is_empty())
        .map(|(row, rows)| SpeciesNetworkBoundary {
            id: row.id,
            parcel: row.parcel,
            catalog: parcels[&row.parcel].model(),
            species: rows,
        })
        .chain(
            reservoirs
                .iter()
                .zip(&reservoir_rows)
                .filter(|(_, rows)| !rows.is_empty())
                .map(|(row, rows)| SpeciesNetworkBoundary {
                    id: row.id,
                    parcel: row.parcel,
                    catalog: row.flux.interface.gas_catalog(),
                    species: rows,
                }),
        )
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
    let step = ReservoirEffusionStep {
        duration_s,
        transport: sampled.rates,
        vacuum_outlets,
        reservoirs,
        exchange,
    };
    step.admit_endpoint_driving_forces(&parcels, &mut cancelled)?;
    check(&mut cancelled)?;
    Ok(ReservoirEffusionCandidate { parcels, step })
}
