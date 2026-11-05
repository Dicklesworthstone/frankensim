//! Pressure-driven and selective real-gas streams on one immutable snapshot.
//!
//! Every nozzle reads the original upstream PR state and resolves its original
//! isentrope, enthalpy drop and sonic/back-pressure condition. All bulk streams
//! are lowered to original species amounts before the common partial-enthalpy
//! network update reconstructs each affected volume exactly once.
//!
//! This is a first-order explicit, quasi-steady interval. Active-path pressure
//! reversal and aggregate inventory exhaustion refuse instead of clipping flux.
//! The caller refines the interval and resamples the rates. No incoming stream
//! can immediately fund another outflow. The receiving volume thermalizes the
//! original stagnation enthalpy; jet kinetic energy is not deposited again.

use super::super::species::network::{
    self as species, RealGasSpeciesNetworkStep, SpeciesNetworkError, SpeciesNetworkExternalHeat,
    SpeciesNetworkHeatTransfer, SpeciesNetworkSolve, SpeciesNetworkTransfer,
};
use super::{GasNozzle, IsochoricRealGasParcel, RealGasNozzleFlow, RealGasNozzleStepSolve};
use crate::thermochemical_gas::energy::parcel::transfer::{
    inventory,
    snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
};
use crate::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use crate::thermochemical_gas::positive;
use std::collections::BTreeMap;

/// One distinct physical restriction between two finite real-gas owners.
/// Array order does not prescribe flow direction; actual upstream pressure does.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasNozzleNetworkLink {
    pub id: u64,
    pub parcels: [u64; 2],
    pub nozzle: GasNozzle,
    /// Original throat/isentropic controls. The supplied node controls own the
    /// final collective endpoint solves; this flow query retains the original
    /// control admission, including its transfer-control validation.
    /// Each active path's donor fraction limits ALL withdrawal from that owner.
    pub solve: RealGasNozzleStepSolve,
}

/// Original source-resolved flow and actual upstream identity for one path.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasNozzleNetworkLinkStep {
    pub id: u64,
    pub donor: u64,
    pub receiver: u64,
    pub transferred_mass_kg: f64,
    pub flow: RealGasNozzleFlow,
}

/// All physical paths and one conservative species/heat inventory account.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasNozzleNetworkStep {
    pub duration_s: f64,
    /// Canonical physical path-ID order, including closed/equal-pressure paths.
    pub links: Vec<RealGasNozzleNetworkLinkStep>,
    pub exchange: RealGasSpeciesNetworkStep,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RealGasNozzleNetworkCandidate {
    parcels: BTreeMap<u64, IsochoricRealGasParcel>,
    step: RealGasNozzleNetworkStep,
}
impl RealGasNozzleNetworkCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, IsochoricRealGasParcel> {
        &self.parcels
    }
    pub const fn step(&self) -> &RealGasNozzleNetworkStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, IsochoricRealGasParcel>,
        RealGasNozzleNetworkStep,
    ) {
        (self.parcels, self.step)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum RealGasNozzleNetworkError {
    Invalid(&'static str),
    Link {
        id: u64,
        error: GasEnergyError,
    },
    Withdrawal {
        donor: u64,
        fraction: f64,
        maximum: f64,
    },
    /// A smaller interval must resolve the first change in flow direction.
    EquilibriumCrossing {
        link: u64,
        donor: u64,
        receiver: u64,
    },
    Exchange(SpeciesNetworkError),
    Cancelled,
}
impl core::fmt::Display for RealGasNozzleNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "real-gas nozzle network: {self:?}")
    }
}
impl std::error::Error for RealGasNozzleNetworkError {}
impl From<SpeciesNetworkError> for RealGasNozzleNetworkError {
    fn from(error: SpeciesNetworkError) -> Self {
        match error {
            SpeciesNetworkError::Cancelled => Self::Cancelled,
            error => Self::Exchange(error),
        }
    }
}
use RealGasNozzleNetworkError as Error;

fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn link_error(id: u64, error: GasEnergyError) -> Error {
    if error == GasEnergyError::Cancelled {
        Error::Cancelled
    } else {
        Error::Link { id, error }
    }
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), Error> {
        let adjusted = term - self.correction;
        let value = self.value + adjusted;
        let correction = (value - self.value) - adjusted;
        if ![term, adjusted, value, correction]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(Error::Invalid(
                "real-gas nozzle network withdrawal is unrepresentable",
            ));
        }
        self.value = value;
        self.correction = correction;
        Ok(())
    }
}

/// Advance all pressure-driven paths over one explicit physical interval.
/// Distinct IDs denote distinct physical flow paths, including parallel paths.
/// All streams use the original donor state and its full stagnation enthalpy.
pub fn advance(
    parcels: &BTreeMap<u64, IsochoricRealGasParcel>,
    links: &[RealGasNozzleNetworkLink],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<RealGasNozzleNetworkCandidate, Error> {
    advance_with_exchange(
        parcels,
        links,
        duration_s,
        &[],
        &[],
        &[],
        endpoints,
        controls,
        cancelled,
    )
}

/// Combine nozzle outflow, authored selective streams and additional heat in
/// the SAME reconstruction. Additional transfers are integrated species masses;
/// heat records are integrated joules over duration_s, not rates. Their physical
/// owner supplies any membrane work, diffusion law or opposite thermal storage.
/// No entropy claim is inferred for arbitrary prescribed species/heat records.
///
/// Every original species allowance applies to gross outgoing amounts before
/// receipts. Additionally, the strictest active nozzle donor-fraction limit
/// applies to the sum of all nozzle AND selective outflow from that donor.
/// The 4096-record snapshot limit includes all generated species and heat rows.
/// Cancelled, unsupported, overdrawn or pressure-crossing trials publish nothing.
#[allow(clippy::too_many_arguments)]
pub fn advance_with_exchange(
    parcels: &BTreeMap<u64, IsochoricRealGasParcel>,
    links: &[RealGasNozzleNetworkLink],
    duration_s: f64,
    transfers: &[SpeciesNetworkTransfer],
    heat: &[SpeciesNetworkHeatTransfer],
    external_heat: &[SpeciesNetworkExternalHeat],
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RealGasNozzleNetworkCandidate, Error> {
    check(&mut cancelled)?;
    let initial_rows = transfers
        .len()
        .checked_add(heat.len())
        .and_then(|n| n.checked_add(external_heat.len()))
        .ok_or(Error::Invalid(
            "real-gas nozzle network record count overflow",
        ))?;
    if !positive(duration_s)
        || parcels.is_empty()
        || parcels.len() > MAX_SNAPSHOT_PARCELS
        || links.len() > MAX_SNAPSHOT_TRANSFERS
        || initial_rows > MAX_SNAPSHOT_TRANSFERS
    {
        return Err(Error::Invalid(
            "real-gas nozzle network needs positive duration and bounded original owners/paths",
        ));
    }
    let mut ordered = BTreeMap::new();
    for link in links {
        check(&mut cancelled)?;
        if link.parcels[0] == link.parcels[1]
            || link.parcels.iter().any(|id| !parcels.contains_key(id))
            || ordered.insert(link.id, link).is_some()
        {
            return Err(Error::Invalid(
                "real-gas nozzle paths require unique IDs and distinct existing owners",
            ));
        }
    }
    // Canonical accumulation also keeps user selective-row order immaterial.
    let mut extra: Vec<_> = transfers.iter().collect();
    extra.sort_by(|a, b| {
        a.donor
            .cmp(&b.donor)
            .then(a.receiver.cmp(&b.receiver))
            .then(a.species.cmp(&b.species))
            .then(a.mass_kg.total_cmp(&b.mass_kg))
    });
    let mut gross: BTreeMap<u64, Sum> = BTreeMap::new();
    for row in extra {
        check(&mut cancelled)?;
        if !row.mass_kg.is_finite() || row.mass_kg < 0.0 || !parcels.contains_key(&row.donor) {
            return Err(Error::Invalid(
                "real-gas nozzle network needs admitted selective donor masses",
            ));
        }
        gross.entry(row.donor).or_default().add(row.mass_kg)?;
    }
    let mut rows = transfers.to_vec();
    let mut steps = Vec::with_capacity(ordered.len());
    let mut limits: BTreeMap<u64, f64> = BTreeMap::new();
    for (id, link) in ordered {
        check(&mut cancelled)?;
        let [a, b] = link.parcels;
        let pressure_a = parcels[&a].state().pressure_pa();
        let pressure_b = parcels[&b].state().pressure_pa();
        let (donor, receiver) = if pressure_a > pressure_b || (pressure_a == pressure_b && a < b) {
            (a, b)
        } else {
            (b, a)
        };
        let upstream = &parcels[&donor];
        let flow = upstream
            .nozzle_flow_to(
                parcels[&receiver].state().pressure_pa(),
                link.nozzle,
                link.solve,
                &mut cancelled,
            )
            .map_err(|e| link_error(id, e))?;
        let mass = flow.mass_flow_kg_s * duration_s;
        if !mass.is_finite() || (flow.mass_flow_kg_s > 0.0 && mass == 0.0) {
            return Err(Error::Invalid(
                "real-gas nozzle network integrated mass is unrepresentable",
            ));
        }
        if mass > 0.0 {
            limits
                .entry(donor)
                .and_modify(|maximum| *maximum = maximum.min(link.solve.maximum_donor_fraction))
                .or_insert(link.solve.maximum_donor_fraction);
            gross.entry(donor).or_default().add(mass)?;
            let count = upstream
                .model()
                .source()
                .thermodynamics()
                .components()
                .len();
            let used = rows
                .len()
                .checked_add(heat.len())
                .and_then(|n| n.checked_add(external_heat.len()))
                .and_then(|n| n.checked_add(count))
                .ok_or(Error::Invalid(
                    "real-gas nozzle network species expansion overflow",
                ))?;
            if used > MAX_SNAPSHOT_TRANSFERS {
                return Err(Error::Invalid(
                    "real-gas nozzle species expansion exceeds the shared snapshot bound",
                ));
            }
            for (species, mass_kg) in
                inventory(upstream.model().source(), mass).map_err(|e| link_error(id, e))?
            {
                rows.push(SpeciesNetworkTransfer {
                    donor,
                    receiver,
                    species,
                    mass_kg,
                });
            }
        }
        steps.push(RealGasNozzleNetworkLinkStep {
            id,
            donor,
            receiver,
            transferred_mass_kg: mass,
            flow,
        });
    }
    for (donor, maximum) in limits {
        check(&mut cancelled)?;
        let fraction = gross[&donor].value / parcels[&donor].mass_kg();
        if !fraction.is_finite() || fraction > maximum {
            return Err(Error::Withdrawal {
                donor,
                fraction,
                maximum,
            });
        }
    }
    let candidate = species::exchange_candidate_with_external_heat(
        parcels,
        &rows,
        heat,
        external_heat,
        endpoints,
        controls,
        &mut cancelled,
    )?;
    for step in &steps {
        check(&mut cancelled)?;
        if step.transferred_mass_kg > 0.0
            && candidate.parcels()[&step.donor].state().pressure_pa()
                < candidate.parcels()[&step.receiver].state().pressure_pa()
        {
            return Err(Error::EquilibriumCrossing {
                link: step.id,
                donor: step.donor,
                receiver: step.receiver,
            });
        }
    }
    check(&mut cancelled)?;
    let (parcels, exchange) = candidate.into_parts();
    Ok(RealGasNozzleNetworkCandidate {
        parcels,
        step: RealGasNozzleNetworkStep {
            duration_s,
            links: steps,
            exchange,
        },
    })
}

#[cfg(test)]
mod tests;
