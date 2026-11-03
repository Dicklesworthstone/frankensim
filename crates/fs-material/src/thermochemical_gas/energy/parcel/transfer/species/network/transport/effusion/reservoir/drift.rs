//! Qualified incoming drift and outgoing effusion on the same original budget.

use super::{
    BoundarySpeciesTransfer, ClosedGasParcel, GasTransportError, MAX_SNAPSHOT_TRANSFERS, Sum,
    check, face_error, relative_energy, signed_product,
};
use crate::thermochemical_gas::energy::parcel::transfer::species::effusion::{
    drift::{DriftingEffusiveFlux, EffusiveDriftInterface},
    vacuum::VacuumEffusiveFlux,
};
use std::collections::BTreeMap;

/// One prescribed Maxwellian source moving relative to a finite gas aperture.
#[derive(Clone, Copy, Debug)]
pub struct EffusionDriftingReservoirLink<'a> {
    /// Unique across every external molecular port in this transaction.
    pub id: u64,
    pub parcel: u64,
    pub reservoir: &'a ClosedGasParcel,
    pub interface: &'a EffusiveDriftInterface,
    pub area_m2: f64,
    /// Source bulk velocity MINUS aperture velocity, sampled in one common basis.
    pub relative_bulk_velocity_m_s: [f64; 3],
    /// Unit normal from the finite gas into the prescribed source.
    pub outward_unit_normal: [f64; 3],
}

/// Both original samples used by one simultaneous finite transaction.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionDriftingReservoirStep {
    pub id: u64,
    pub parcel: u64,
    /// The finite gas is locally at rest relative to the aperture.
    pub outgoing: VacuumEffusiveFlux,
    /// Distinct shifted-Maxwellian incidence; never a rest-beam momentum law.
    pub incoming: DriftingEffusiveFlux,
}

pub(super) fn sample(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    links: &[EffusionDriftingReservoirLink<'_>],
    duration_s: f64,
    work: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<
    (
        Vec<EffusionDriftingReservoirStep>,
        Vec<Vec<BoundarySpeciesTransfer>>,
    ),
    GasTransportError,
> {
    // The enclosing owner has already admitted identities and the graph bound.
    let ordered: BTreeMap<_, _> = links.iter().map(|link| (link.id, link)).collect();
    let mut reports = Vec::with_capacity(links.len());
    let mut boundaries = Vec::with_capacity(links.len());
    for (id, link) in ordered {
        check(cancelled)?;
        let gas = &parcels[&link.parcel];
        let interface = link.interface.momentum_interface().interface();
        let outgoing = gas
            .effusive_flux_to_vacuum(interface, link.area_m2, &mut *cancelled)
            .map_err(|e| face_error(id, e))?;
        let incoming = link
            .reservoir
            .drifting_effusive_flux(
                link.interface,
                link.area_m2,
                link.relative_bulk_velocity_m_s,
                link.outward_unit_normal,
                &mut *cancelled,
            )
            .map_err(|e| face_error(id, e))?;
        let (tl, tr) = (
            gas.state().state().temperature,
            link.reservoir.state().state().temperature,
        );
        let mut rows = Vec::new();
        for original in interface.gas_catalog().thermodynamics().components() {
            check(cancelled)?;
            let out = outgoing
                .species
                .binary_search_by(|s| s.species.cmp(original.species()))
                .ok()
                .map(|i| &outgoing.species[i].outward);
            let inc = incoming
                .species
                .binary_search_by(|s| s.species.cmp(original.species()))
                .ok()
                .map(|i| &incoming.species[i]);
            let out_rate = out.map_or(0.0, |b| b.mass_rate_kg_s);
            let in_rate = inc.map_or(0.0, |b| b.incoming.mass_rate_kg_s);
            if in_rate == 0.0 && out_rate == 0.0 {
                continue;
            }
            let mi = signed_product(id, in_rate, duration_s)?;
            let mo = signed_product(id, out_rate, duration_s)?;
            let mut net = mi - mo;
            let mut kinetic_change = 0.0;
            if let Some(inc) = inc.filter(|b| b.incoming.mass_rate_kg_s > 0.0) {
                let baseline = inc.baseline_incoming.mass_rate_kg_s;
                let change = inc.incoming_mass_rate_change_from_zero_drift_kg_s;
                // Near zero drift, keep the small change even if the two gross
                // streams round equally. In the far outward tail, use the
                // DIRECT incoming mass: baseline + change would erase it.
                if change.abs() <= 0.5 * baseline {
                    let mut sum = Sum::default();
                    sum.add(signed_product(id, baseline, duration_s)? - mo)
                        .map_err(|e| face_error(id, e))?;
                    sum.add(signed_product(id, change, duration_s)?)
                        .map_err(|e| face_error(id, e))?;
                    net = sum.value;
                }
                let correction = |beam: &crate::thermochemical_gas::energy::parcel::transfer::species::effusion::EffusiveBeam| {
                    beam.additional_specific_energy_j_kg.ok_or(GasTransportError::Invalid(
                        "drifting reservoir lost its original molecular energy correction",
                    ))
                };
                kinetic_change = correction(&inc.incoming)? - correction(&inc.baseline_incoming)?;
            }
            let relative = relative_energy(id, original, tl, tr, net, mi, kinetic_change)?;
            *work += 1 + usize::from(rows.is_empty());
            if *work > MAX_SNAPSHOT_TRANSFERS {
                return Err(GasTransportError::Invalid(
                    "combined drifting molecular species exceed shared work bound",
                ));
            }
            rows.push(BoundarySpeciesTransfer {
                species: original.species().clone(),
                incoming_mass_kg: mi,
                outgoing_mass_kg: mo,
                net_mass_into_gas_kg: net,
                relative_energy_into_gas_j: relative,
            });
        }
        reports.push(EffusionDriftingReservoirStep {
            id,
            parcel: link.parcel,
            outgoing,
            incoming,
        });
        boundaries.push(rows);
    }
    Ok((reports, boundaries))
}

pub(super) fn admit_endpoints(
    rows: &[EffusionDriftingReservoirStep],
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), GasTransportError> {
    for row in rows {
        check(cancelled)?;
        if row.outgoing.species.is_empty() && row.incoming.species.is_empty() {
            continue;
        }
        let gas = parcels.get(&row.parcel).ok_or(GasTransportError::Invalid(
            "drifting reservoir endpoint lost its finite gas",
        ))?;
        let interface = row.incoming.interface.momentum_interface().interface();
        interface
            .admit_parcel(gas, &mut *cancelled)
            .map_err(|e| face_error(row.id, e))?;
        // Calibrated zero transmissions still retain their authored support,
        // including species absent from the other side of the original port.
        for species in row
            .outgoing
            .species
            .iter()
            .map(|s| &s.species)
            .chain(row.incoming.species.iter().map(|s| &s.species))
        {
            check(cancelled)?;
            let index = interface
                .transmissions()
                .binary_search_by(|t| t.species().cmp(species))
                .map_err(|_| {
                    GasTransportError::Invalid("drifting source lost original transmission")
                })?;
            let state = gas.state().state();
            interface.transmissions()[index]
                .admit(state.temperature, state.pressure)
                .map_err(|e| face_error(row.id, e))?;
        }
    }
    check(cancelled)?;
    Ok(())
}
