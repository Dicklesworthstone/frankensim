//! Finite gas blowdown through fully accommodating circular molecular channels.
//!
//! Native paths and every incident channel beam spend the same ORIGINAL species
//! inventory. Returned molecules count as gross withdrawals and arrivals; only
//! escaping molecules change mass. Their return heat enters the SAME caloric
//! reconstruction as the original molecular escape energy.
//!
//! The vacuum is an absorbing boundary with no invented gas state or entropy.
//! Escaping molecules and the prescribed wall thermostat have explicit energy
//! accounts. A source-relative original-to-final audit cancels formation datums
//! before arithmetic, including across several accepted steps.

use super::super::{
    check, face_error,
    vacuum::{self as thin_vacuum, VacuumEffusionLink},
};
use crate::thermochemical_gas::energy::parcel::transfer::{
    snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
    species::{
        effusion::{
            channel::{DiffuseMolecularChannel, vacuum::VacuumChannelWallExchange},
            vacuum::VacuumEffusiveFlux,
        },
        network::{
            SpeciesNetworkError, SpeciesNetworkExternalHeat, SpeciesNetworkStep,
            boundary::SpeciesNetworkBoundary,
            exchange_candidate_with_external_heat,
            transport::{
                self, GasTransportError, GasTransportPaths, GasTransportRates, GasTransportSolve,
            },
        },
        source_error,
    },
};
use crate::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_math::dd::Dd;
use fs_qty::Temperature;
use fs_thermochem::{
    Nasa9StandardStateModelV1, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};
use std::collections::{BTreeMap, BTreeSet};

/// One physical array of circular channels emptying into an absorbing vacuum.
#[derive(Clone, Copy, Debug)]
pub struct VacuumChannelNetworkLink<'a> {
    pub id: u64,
    pub parcel: u64,
    pub channel: &'a DiffuseMolecularChannel,
    /// Sum of physical mouth areas [m2]; zero closes the entire outlet.
    pub area_m2: f64,
    pub wall_temperature_k: f64,
}

/// Actual integrated escaping mass, separate from the reflected gross inventory.
#[derive(Clone, Debug, PartialEq)]
pub struct VacuumChannelEscapedSpecies {
    pub species: SpeciesId,
    pub mass_kg: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct VacuumChannelOutletStep {
    pub id: u64,
    pub parcel: u64,
    pub channel: DiffuseMolecularChannel,
    /// Original gas temperature retained for independent complete-window audits.
    pub source_temperature_k: f64,
    pub wall_temperature_k: f64,
    /// Original source-energy beam after geometric transmission. Accommodation
    /// of escaped molecules is additional and appears in `wall`.
    pub transmitted: VacuumEffusiveFlux,
    pub wall: VacuumChannelWallExchange,
    pub escaped_species: Vec<VacuumChannelEscapedSpecies>,
    /// Complete escaped energy including accommodation [J]. This diagnostic
    /// depends on the original species energy datums; balance residuals do not.
    pub escaped_energy_j: f64,
}

/// Actual finite gas storage plus escaped molecules and the wall thermostats.
/// No entropy is assigned to vacuum or to an unmodeled outgoing gas inventory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VacuumChannelWindowAudit {
    pub escaped_mass_kg: f64,
    pub escaped_energy_j: f64,
    pub heat_into_walls_j: f64,
    pub total_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    /// Delta U_gas + E_escaped + Q_wall, with every species datum removed
    /// BEFORE arithmetic. Actual mass rounding has its own independent budget.
    pub energy_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct VacuumChannelNetworkStep {
    pub duration_s: f64,
    pub transport: GasTransportRates,
    pub outlets: Vec<VacuumChannelOutletStep>,
    pub exchange: SpeciesNetworkStep,
    pub audit: VacuumChannelWindowAudit,
}

impl VacuumChannelNetworkStep {
    /// Recheck gas, collision and active-species wall calorics after coupling.
    /// Species supplied by native paths must support this held wall even when
    /// they were absent from the original outlet beam. This spends no additional
    /// inventory and samples no additional stream.
    pub fn admit_endpoint_driving_forces(
        &self,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), GasTransportError> {
        self.transport
            .admit_endpoint_driving_forces(parcels, &mut cancelled)?;
        for row in &self.outlets {
            check(&mut cancelled)?;
            if row.transmitted.species.is_empty() {
                continue;
            }
            let gas = parcels.get(&row.parcel).ok_or(GasTransportError::Invalid(
                "vacuum channel endpoint lost its original gas identity",
            ))?;
            row.channel
                .admit_wall_for_parcel(gas, row.wall_temperature_k, &mut cancelled)
                .map_err(|error| face_error(row.id, error))?;
        }
        check(&mut cancelled)?;
        Ok(())
    }
}

/// All reconstructed finite gases or none; caller inventories remain immutable.
#[derive(Clone, Debug, PartialEq)]
pub struct VacuumChannelNetworkCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    step: VacuumChannelNetworkStep,
}
impl VacuumChannelNetworkCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn step(&self) -> &VacuumChannelNetworkStep {
        &self.step
    }
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, VacuumChannelNetworkStep) {
        (self.parcels, self.step)
    }
}

/// Advance native internal paths and thermal vacuum channels simultaneously.
///
/// Geometric transmission multiplies the exposed area; the original vacuum
/// boundary still carries h_i(T)-R*T/(2*M_i). Reflected molecules are added to
/// BOTH gross directions without changing that authoritative net escape. Only
/// their accommodation heat enters the finite gas. Accommodation of transmitted
/// molecules instead changes escaped energy, with opposite thermostat heat.
///
/// Every path describes a separately exposed physical passage. This is an
/// explicit first-order, quasistatic stage; callers refine time. Complete
/// depletion, source extrapolation, overspent gross inventories, failed energy
/// budgets and cancellation refuse atomically. A closed-network entropy check
/// applies only when no outlet emits; no gas-only entropy gate is imposed on
/// an open vacuum system.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    paths: GasTransportPaths<'_>,
    links: &[VacuumChannelNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: GasTransportSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<VacuumChannelNetworkCandidate, GasTransportError> {
    check(&mut cancelled)?;
    let count = [
        paths.nozzles.len(),
        paths.diffusion.len(),
        paths.porous.len(),
        paths.membranes.len(),
        paths.dusty_gas.len(),
        paths.layered_dusty_gas.len(),
        paths.effusion.len(),
        links.len(),
    ]
    .into_iter()
    .try_fold(0usize, |total, n| total.checked_add(n));
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || parcels.is_empty()
        || parcels.len() > MAX_SNAPSHOT_PARCELS
        || count.is_none_or(|n| n > MAX_SNAPSHOT_TRANSFERS)
    {
        return Err(GasTransportError::Invalid(
            "vacuum channels need a bounded nonempty graph and positive duration",
        ));
    }
    let mut ordered = BTreeMap::new();
    for link in links {
        check(&mut cancelled)?;
        if !parcels.contains_key(&link.parcel) || ordered.insert(link.id, link).is_some() {
            return Err(GasTransportError::Invalid(
                "vacuum channels need unique IDs and existing finite gases",
            ));
        }
    }
    let mut molecular = Vec::with_capacity(links.len());
    for link in ordered.values() {
        check(&mut cancelled)?;
        let area_m2 = link
            .channel
            .transmitted_area_m2(link.area_m2)
            .map_err(|error| face_error(link.id, error))?;
        if area_m2 > 0.0 {
            link.channel
                .admit_parcel(&parcels[&link.parcel], &mut cancelled)
                .map_err(|error| face_error(link.id, error))?;
        }
        molecular.push(VacuumEffusionLink {
            id: link.id,
            parcel: link.parcel,
            interface: link.channel.mouth(),
            area_m2,
        });
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
    let (reports, mut species_rows) =
        thin_vacuum::sample_outlets(parcels, &molecular, duration_s, &mut work, &mut cancelled)?;
    let mut outlets = Vec::with_capacity(links.len());
    let mut external_heat = Vec::with_capacity(links.len());
    for ((link, report), rows) in ordered.values().zip(reports).zip(&mut species_rows) {
        check(&mut cancelled)?;
        let gas = &parcels[&link.parcel];
        let wall = link
            .channel
            .integrate_vacuum_wall_exchange(
                gas,
                &report.flux,
                link.wall_temperature_k,
                duration_s,
                &mut cancelled,
            )
            .map_err(|error| face_error(link.id, error))?;
        let mut escaped_species = Vec::with_capacity(rows.len());
        let mut escaped_energy = Dd::from_f64(wall.heat_into_escaped_j);
        for row in rows.iter() {
            check(&mut cancelled)?;
            let beam = report
                .flux
                .species
                .iter()
                .find(|entry| entry.species == row.species)
                .ok_or(GasTransportError::Invalid(
                    "vacuum channel lost an escaped species beam",
                ))?;
            let h = beam
                .outward
                .source_enthalpy_j_kg
                .ok_or(GasTransportError::Invalid(
                    "vacuum channel beam lost its original partial enthalpy",
                ))?;
            // Use the ACTUAL mass and relative energy selected by the shared
            // boundary owner, not a separately rounded energy_rate * duration.
            escaped_energy = escaped_energy + Dd::from_f64(row.outgoing_mass_kg) * Dd::from_f64(h)
                - Dd::from_f64(row.relative_energy_into_gas_j);
            escaped_species.push(VacuumChannelEscapedSpecies {
                species: row.species.clone(),
                mass_kg: row.outgoing_mass_kg,
            });
        }
        for reflected in &wall.reflected_species {
            check(&mut cancelled)?;
            if reflected.mass_kg == 0.0 {
                continue;
            }
            let row = rows
                .iter_mut()
                .find(|entry| entry.species == reflected.species)
                .ok_or(GasTransportError::Invalid(
                    "vacuum channel returned a species absent from its incident beam",
                ))?;
            row.incoming_mass_kg = reflected.mass_kg;
            row.outgoing_mass_kg =
                (Dd::from_f64(row.outgoing_mass_kg) + Dd::from_f64(reflected.mass_kg)).to_f64();
            // The explicit net mass and molecular energy correction remain
            // the original escape values, even for nearly complete reflection.
        }
        if wall.heat_into_gas_j != 0.0 {
            work = work.checked_add(1).ok_or(GasTransportError::Invalid(
                "vacuum channel heat count overflowed",
            ))?;
            external_heat.push(SpeciesNetworkExternalHeat {
                id: link.id,
                parcel: link.parcel,
                heat_into_gas_j: wall.heat_into_gas_j,
            });
        }
        if work > MAX_SNAPSHOT_TRANSFERS || !escaped_energy.to_f64().is_finite() {
            return Err(GasTransportError::Invalid(
                "vacuum channel heat/energy exceeded its representable work bound",
            ));
        }
        outlets.push(VacuumChannelOutletStep {
            id: link.id,
            parcel: link.parcel,
            channel: link.channel.clone(),
            source_temperature_k: gas.state().state().temperature,
            wall_temperature_k: link.wall_temperature_k,
            transmitted: report.flux,
            wall,
            escaped_species,
            escaped_energy_j: escaped_energy.to_f64(),
        });
    }
    let boundaries: Vec<_> = outlets
        .iter()
        .zip(&species_rows)
        .filter(|(_, rows)| !rows.is_empty())
        .map(|(outlet, rows)| SpeciesNetworkBoundary {
            id: outlet.id,
            parcel: outlet.parcel,
            catalog: parcels[&outlet.parcel].model(),
            species: rows,
        })
        .collect();
    let candidate = exchange_candidate_with_external_heat(
        parcels,
        &sampled.transfers,
        &sampled.heat,
        &boundaries,
        &external_heat,
        endpoints,
        controls.exchange,
        &mut cancelled,
    )?;
    if boundaries.is_empty() {
        transport::audit_closed_window(parcels, candidate.parcels(), controls, &mut cancelled)?;
    }
    let audit = audit_window(
        parcels,
        candidate.parcels(),
        &outlets,
        controls,
        &mut cancelled,
    )?;
    let (parcels, exchange) = candidate.into_parts();
    let step = VacuumChannelNetworkStep {
        duration_s,
        transport: sampled.rates,
        outlets,
        exchange,
        audit,
    };
    step.admit_endpoint_driving_forces(&parcels, &mut cancelled)?;
    check(&mut cancelled)?;
    Ok(VacuumChannelNetworkCandidate { parcels, step })
}

struct CaloricReference<'a> {
    source: &'a Nasa9StandardStateModelV1,
    temperature_k: f64,
    mass_residual: Dd,
}

fn same_source(a: &Nasa9StandardStateModelV1, b: &Nasa9StandardStateModelV1) -> bool {
    a.species() == b.species()
        && a.molar_mass().value().to_bits() == b.molar_mass().value().to_bits()
        && a.convention() == b.convention()
        && a.caloric_reference() == b.caloric_reference()
        && a.regions()
            .iter()
            .map(|r| r.card_identity())
            .eq(b.regions().iter().map(|r| r.card_identity()))
}

fn amount(gas: &ClosedGasParcel, index: usize) -> Result<f64, GasTransportError> {
    let fraction = gas.model().thermodynamics().mole_composition().fractions()[index];
    let m = gas
        .state()
        .thermodynamics()
        .properties()
        .molar_mass()
        .value();
    let n = (gas.mass_kg() / m) * fraction;
    if n.is_finite() && n >= 0.0 && (fraction == 0.0 || n > 0.0) {
        Ok(n)
    } else {
        Err(GasTransportError::Invalid(
            "vacuum channel audit species amount vanished or overflowed",
        ))
    }
}

fn du(
    source: &Nasa9StandardStateModelV1,
    from: f64,
    to: f64,
    parcel: u64,
) -> Result<f64, GasTransportError> {
    source
        .molar_internal_energy_change(Temperature::new(from), Temperature::new(to))
        .map(|value| value.value())
        .map_err(|error| {
            GasTransportError::Exchange(SpeciesNetworkError::Parcel {
                id: parcel,
                error: source_error(error),
            })
        })
}

/// Audit actual original/final gas storage against accepted channel outlets.
///
/// A window may flatten the outlet records of several accepted steps. Each
/// record retains its own source temperature and actual escaped species masses;
/// absolute tolerances are never multiplied by the number of steps. Native
/// internal paths cancel within the gas inventory. No additional external mass,
/// heat, work, reactions or moving volumes may be omitted from this balance.
///
/// The common original internal-energy reference for each source-qualified
/// species cancels analytically against its escaped mass. The remaining escape
/// term is dm/M * [Delta u(T_ref,T_source)+R*T_source/2], plus accommodation.
/// Large formation constants therefore cannot erase cooling or wall heat.
/// This audit makes no claim about vacuum entropy or time-discretization error.
pub fn audit_window(
    before: &BTreeMap<u64, ClosedGasParcel>,
    after: &BTreeMap<u64, ClosedGasParcel>,
    outlets: &[VacuumChannelOutletStep],
    controls: GasTransportSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<VacuumChannelWindowAudit, GasTransportError> {
    audit_window_borrowed(before, after, outlets.iter(), controls, cancelled)
}

/// Share the native source-relative audit without cloning source catalogs in
/// each outlet of an accepted trajectory.
pub(super) fn audit_window_borrowed<'a>(
    before: &BTreeMap<u64, ClosedGasParcel>,
    after: &BTreeMap<u64, ClosedGasParcel>,
    outlets: impl Iterator<Item = &'a VacuumChannelOutletStep> + Clone,
    controls: GasTransportSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<VacuumChannelWindowAudit, GasTransportError> {
    check(&mut cancelled)?;
    let work = outlets.clone().try_fold(0usize, |sum, row| {
        sum.checked_add(1 + row.escaped_species.len())
    });
    if before.is_empty()
        || before.len() > MAX_SNAPSHOT_PARCELS
        || !before.keys().eq(after.keys())
        || work.is_none_or(|n| n > MAX_SNAPSHOT_TRANSFERS)
        || [
            controls.exchange.mass_tolerance_kg,
            controls.exchange.energy_tolerance_j,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v <= 0.0)
    {
        return Err(GasTransportError::Invalid(
            "vacuum channel audit needs bounded original identities and fixed positive budgets",
        ));
    }
    let mut references: BTreeMap<SpeciesId, Vec<CaloricReference<'_>>> = BTreeMap::new();
    let mut total_mass = Dd::ZERO;
    for (&id, old) in before {
        check(&mut cancelled)?;
        let next = &after[&id];
        if old.volume_m3() != next.volume_m3() {
            return Err(GasTransportError::Invalid(
                "vacuum channel audit requires unchanged finite gas volumes",
            ));
        }
        total_mass = total_mass + Dd::from_f64(next.mass_kg()) - Dd::from_f64(old.mass_kg());
        for (index, source) in old.model().thermodynamics().components().iter().enumerate() {
            check(&mut cancelled)?;
            if amount(old, index)? == 0.0 {
                continue;
            }
            let group = references.entry(source.species().clone()).or_default();
            if !group.iter().any(|r| same_source(r.source, source)) {
                group.push(CaloricReference {
                    source,
                    temperature_k: old.state().state().temperature,
                    mass_residual: Dd::ZERO,
                });
            }
        }
    }
    let mut energy = Dd::ZERO;
    for (&id, old) in before {
        check(&mut cancelled)?;
        let next = &after[&id];
        let original = old.model().thermodynamics().components();
        let final_sources = next.model().thermodynamics().components();
        let (mut i, mut j) = (0, 0);
        while i < original.len() || j < final_sources.len() {
            check(&mut cancelled)?;
            let order = match (original.get(i), final_sources.get(j)) {
                (Some(a), Some(b)) => a.species().cmp(b.species()),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => unreachable!(),
            };
            let (source, n0, n1) = match order {
                std::cmp::Ordering::Less => {
                    let row = (&original[i], amount(old, i)?, 0.0);
                    i += 1;
                    row
                }
                std::cmp::Ordering::Greater => {
                    let row = (&final_sources[j], 0.0, amount(next, j)?);
                    j += 1;
                    row
                }
                std::cmp::Ordering::Equal => {
                    if !same_source(&original[i], &final_sources[j]) {
                        return Err(GasTransportError::Invalid(
                            "vacuum channel audit endpoint changed an original caloric source",
                        ));
                    }
                    let row = (&original[i], amount(old, i)?, amount(next, j)?);
                    i += 1;
                    j += 1;
                    row
                }
            };
            if n0 == 0.0 && n1 == 0.0 {
                continue;
            }
            let reference = references
                .get_mut(source.species())
                .and_then(|group| group.iter_mut().find(|r| same_source(r.source, source)))
                .ok_or(GasTransportError::Invalid(
                    "vacuum channel audit introduced an unqualified species source",
                ))?;
            let t0 = old.state().state().temperature;
            let t1 = next.state().state().temperature;
            if n0 == 0.0 {
                energy = energy
                    + Dd::from_f64(n1) * Dd::from_f64(du(source, reference.temperature_k, t1, id)?);
            } else if n1 == 0.0 {
                energy = energy
                    - Dd::from_f64(n0) * Dd::from_f64(du(source, reference.temperature_k, t0, id)?);
            } else {
                energy = energy
                    + Dd::from_f64(n1) * Dd::from_f64(du(source, t0, t1, id)?)
                    + (Dd::from_f64(n1) - Dd::from_f64(n0))
                        * Dd::from_f64(du(source, reference.temperature_k, t0, id)?);
            }
            reference.mass_residual = reference.mass_residual
                + (Dd::from_f64(n1) - Dd::from_f64(n0)) * Dd::from_f64(source.molar_mass().value());
        }
    }
    let mut escaped_mass = Dd::ZERO;
    let mut escaped_energy = Dd::ZERO;
    let mut wall_heat = Dd::ZERO;
    for outlet in outlets {
        check(&mut cancelled)?;
        if !before.contains_key(&outlet.parcel)
            || ![
                outlet.wall.heat_into_gas_j,
                outlet.wall.heat_into_escaped_j,
                outlet.wall.heat_into_wall_j,
                outlet.escaped_energy_j,
            ]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(GasTransportError::Invalid(
                "vacuum channel audit lost an original gas or finite energy account",
            ));
        }
        let wall_balance = (Dd::from_f64(outlet.wall.heat_into_gas_j)
            + Dd::from_f64(outlet.wall.heat_into_escaped_j)
            + Dd::from_f64(outlet.wall.heat_into_wall_j))
        .to_f64();
        if !wall_balance.is_finite() || wall_balance.abs() > controls.exchange.energy_tolerance_j {
            return Err(GasTransportError::Invalid(
                "vacuum channel wall account did not retain its opposite heat",
            ));
        }
        let mut unique = BTreeSet::new();
        for row in &outlet.escaped_species {
            check(&mut cancelled)?;
            if !row.mass_kg.is_finite() || row.mass_kg <= 0.0 || !unique.insert(&row.species) {
                return Err(GasTransportError::Invalid(
                    "vacuum channel audit needs unique positive escaped species masses",
                ));
            }
            let source = outlet
                .channel
                .mouth()
                .gas_catalog()
                .thermodynamics()
                .components()
                .iter()
                .find(|model| model.species() == &row.species)
                .ok_or(GasTransportError::Invalid(
                    "vacuum channel audit escape lost its original source catalog",
                ))?;
            let reference = references
                .get_mut(&row.species)
                .and_then(|group| group.iter_mut().find(|r| same_source(r.source, source)))
                .ok_or(GasTransportError::Invalid(
                    "vacuum channel escaped species does not belong to the original inventory",
                ))?;
            let n = row.mass_kg / source.molar_mass().value();
            if !n.is_finite() || n <= 0.0 {
                return Err(GasTransportError::Invalid(
                    "vacuum channel escaped amount is unrepresentable",
                ));
            }
            let sensible = du(
                source,
                reference.temperature_k,
                outlet.source_temperature_k,
                outlet.parcel,
            )?;
            energy = energy
                + Dd::from_f64(n)
                    * (Dd::from_f64(sensible)
                        + Dd::from_f64(0.5 * R) * Dd::from_f64(outlet.source_temperature_k));
            reference.mass_residual = reference.mass_residual + Dd::from_f64(row.mass_kg);
            escaped_mass = escaped_mass + Dd::from_f64(row.mass_kg);
        }
        energy = energy
            + Dd::from_f64(outlet.wall.heat_into_escaped_j)
            + Dd::from_f64(outlet.wall.heat_into_wall_j);
        escaped_energy = escaped_energy + Dd::from_f64(outlet.escaped_energy_j);
        wall_heat = wall_heat + Dd::from_f64(outlet.wall.heat_into_wall_j);
    }
    let mut maximum_species_mass_residual_kg = 0.0_f64;
    for group in references.into_values() {
        let mut species_residual = Dd::ZERO;
        for reference in group {
            check(&mut cancelled)?;
            let residual = reference.mass_residual.to_f64();
            if !residual.is_finite() {
                return Err(GasTransportError::Invalid(
                    "vacuum channel species mass residual is unrepresentable",
                ));
            }
            maximum_species_mass_residual_kg = maximum_species_mass_residual_kg.max(residual.abs());
            species_residual = species_residual + reference.mass_residual;
        }
        maximum_species_mass_residual_kg =
            maximum_species_mass_residual_kg.max(species_residual.to_f64().abs());
    }
    let result = VacuumChannelWindowAudit {
        escaped_mass_kg: escaped_mass.to_f64(),
        escaped_energy_j: escaped_energy.to_f64(),
        heat_into_walls_j: wall_heat.to_f64(),
        total_mass_residual_kg: (total_mass + escaped_mass).to_f64(),
        maximum_species_mass_residual_kg,
        energy_residual_j: energy.to_f64(),
    };
    if ![
        result.escaped_mass_kg,
        result.escaped_energy_j,
        result.heat_into_walls_j,
        result.total_mass_residual_kg,
        result.maximum_species_mass_residual_kg,
        result.energy_residual_j,
    ]
    .into_iter()
    .all(f64::is_finite)
        || result.total_mass_residual_kg.abs() > controls.exchange.mass_tolerance_kg
        || result.maximum_species_mass_residual_kg > controls.exchange.mass_tolerance_kg
        || result.energy_residual_j.abs() > controls.exchange.energy_tolerance_j
    {
        return Err(GasTransportError::Invalid(
            "vacuum channel window exceeded its fixed mass or internal-energy budget",
        ));
    }
    check(&mut cancelled)?;
    Ok(result)
}

#[cfg(test)]
mod tests;
