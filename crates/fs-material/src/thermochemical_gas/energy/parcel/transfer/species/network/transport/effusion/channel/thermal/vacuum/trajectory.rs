//! Adaptive molecular blowdown with simultaneous finite gas and wall feedback.
//!
//! The same finite-wall time driver used for closed channels compares one full
//! step with two resampled half steps. Only accepted halves enter the output.
//! Error control includes gas and wall states, each outlet's escaped/returned
//! species, source-relative escaped energy and all accommodation heat channels.
//! One original-to-final mass/energy audit covers the entire requested window;
//! its tolerances are not multiplied by the number of accepted steps.
//!
//! This is local first-order temporal error estimation, not a rigorous global
//! bound. Geometry, exposed areas and native internal transport paths remain
//! fixed. Vacuum has no assigned entropy, temperature or storage law.

use super::super::super::vacuum::{
    self as native, VacuumChannelOutletStep, VacuumChannelWindowAudit,
};
use super::super::trajectory::{advance_history, ratio};
use super::super::{FiniteChannelError, FixedVolumeChannelWall};
use super::{FiniteVacuumChannelLink, FiniteVacuumChannelStep, step_candidate};
use crate::thermochemical_gas::energy::{
    GasEnergySolve,
    parcel::{
        ClosedGasParcel,
        transfer::{
            snapshot::MAX_SNAPSHOT_TRANSFERS,
            species::network::transport::{GasTransportPaths, GasTransportSolve},
        },
    },
};
use fs_math::dd::Dd;
use fs_qty::Temperature;
use fs_thermal::enthalpy::ThermalLaw;
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use std::collections::BTreeMap;

pub use super::super::trajectory::{EffusionAdvanceControls, FiniteChannelAdvanceError};
type Error = FiniteChannelAdvanceError;

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteVacuumChannelWindowBalance {
    /// Native source-relative gas/escape balance with sampled wall heat.
    pub gas: VacuumChannelWindowAudit,
    pub wall_energy_change_j: f64,
    /// Material storage entropy only; no total open-system entropy assertion.
    pub wall_entropy_change_j_k: f64,
    /// Native gas/escape balance with ACTUAL finite material wall storage.
    pub energy_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteVacuumChannelAdvanceStep {
    pub duration_s: f64,
    pub attempted_steps: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub maximum_accepted_error: f64,
    pub accepted_steps: Vec<FiniteVacuumChannelStep>,
    pub balance: FiniteVacuumChannelWindowBalance,
}

#[derive(Clone, Debug)]
pub struct FiniteVacuumChannelAdvanceCandidate<'l, 'p> {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    walls: BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    step: FiniteVacuumChannelAdvanceStep,
}
impl<'l, 'p> FiniteVacuumChannelAdvanceCandidate<'l, 'p> {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub fn walls(&self) -> &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>> {
        &self.walls
    }
    pub const fn step(&self) -> &FiniteVacuumChannelAdvanceStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, ClosedGasParcel>,
        BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
        FiniteVacuumChannelAdvanceStep,
    ) {
        (self.parcels, self.walls, self.step)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

fn mass_error<'a>(
    coarse: impl Iterator<Item = (&'a SpeciesId, f64)>,
    fine: impl Iterator<Item = (&'a SpeciesId, f64)>,
    controls: EffusionAdvanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let mut values = BTreeMap::<&SpeciesId, [Dd; 2]>::new();
    for (side, (species, mass)) in coarse.map(|row| (0, row)).chain(fine.map(|row| (1, row))) {
        poll(cancelled)?;
        let value = values.entry(species).or_insert([Dd::ZERO; 2]);
        value[side] = value[side] + Dd::from_f64(mass);
    }
    let mut error = 0.0_f64;
    for [a, b] in values.into_values() {
        poll(cancelled)?;
        let (a, b) = (a.to_f64(), b.to_f64());
        error = error.max(ratio(
            a,
            b,
            a.abs().max(b.abs()),
            controls.species_mass_tolerance_kg,
            controls.relative_tolerance,
        )?);
    }
    Ok(error)
}

/// Remove every formation datum BEFORE summing escaped energy. One reference
/// temperature per outlet is shared by the coarse step and both fine steps.
fn relative_escape(
    outlet: &VacuumChannelOutletStep,
    reference: &VacuumChannelOutletStep,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Dd, Error> {
    let mut energy = Dd::from_f64(outlet.wall.heat_into_escaped_j);
    for row in &outlet.escaped_species {
        poll(cancelled)?;
        let source = reference
            .channel
            .mouth()
            .gas_catalog()
            .thermodynamics()
            .components()
            .iter()
            .find(|source| source.species() == &row.species)
            .ok_or(Error::Invalid(
                "adaptive vacuum escape lost its original species source",
            ))?;
        let du = source
            .molar_internal_energy_change(
                Temperature::new(reference.source_temperature_k),
                Temperature::new(outlet.source_temperature_k),
            )
            .map_err(|_| {
                Error::Invalid("adaptive vacuum escape left its original caloric support")
            })?;
        energy = energy
            + Dd::from_f64(row.mass_kg / source.molar_mass().value())
                * (Dd::from_f64(du.value())
                    + Dd::from_f64(0.5 * R) * Dd::from_f64(outlet.source_temperature_k));
    }
    Ok(energy)
}

fn stream_error(
    coarse: &FiniteVacuumChannelStep,
    first: &FiniteVacuumChannelStep,
    second: &FiniteVacuumChannelStep,
    controls: EffusionAdvanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    if coarse.gas.outlets.len() != first.gas.outlets.len()
        || coarse.gas.outlets.len() != second.gas.outlets.len()
    {
        return Err(Error::Invalid("adaptive vacuum outlet identities changed"));
    }
    let mut error = 0.0_f64;
    for ((a, b), c) in coarse
        .gas
        .outlets
        .iter()
        .zip(&first.gas.outlets)
        .zip(&second.gas.outlets)
    {
        poll(cancelled)?;
        if a.id != b.id
            || a.id != c.id
            || a.parcel != b.parcel
            || a.parcel != c.parcel
            || a.channel != b.channel
            || a.channel != c.channel
        {
            return Err(Error::Invalid(
                "adaptive vacuum changed its original channel",
            ));
        }
        error = error.max(mass_error(
            a.escaped_species
                .iter()
                .map(|row| (&row.species, row.mass_kg)),
            b.escaped_species
                .iter()
                .chain(&c.escaped_species)
                .map(|row| (&row.species, row.mass_kg)),
            controls,
            cancelled,
        )?);
        error = error.max(mass_error(
            a.wall
                .reflected_species
                .iter()
                .map(|row| (&row.species, row.mass_kg)),
            b.wall
                .reflected_species
                .iter()
                .chain(&c.wall.reflected_species)
                .map(|row| (&row.species, row.mass_kg)),
            controls,
            cancelled,
        )?);
        let coarse_energy = relative_escape(a, a, cancelled)?.to_f64();
        let fine_energy =
            (relative_escape(b, a, cancelled)? + relative_escape(c, a, cancelled)?).to_f64();
        for (x, y) in [
            (coarse_energy, fine_energy),
            (
                a.wall.heat_into_gas_j,
                (Dd::from_f64(b.wall.heat_into_gas_j) + Dd::from_f64(c.wall.heat_into_gas_j))
                    .to_f64(),
            ),
            (
                a.wall.heat_into_escaped_j,
                (Dd::from_f64(b.wall.heat_into_escaped_j)
                    + Dd::from_f64(c.wall.heat_into_escaped_j))
                .to_f64(),
            ),
            (
                a.wall.heat_into_wall_j,
                (Dd::from_f64(b.wall.heat_into_wall_j) + Dd::from_f64(c.wall.heat_into_wall_j))
                    .to_f64(),
            ),
        ] {
            error = error.max(ratio(
                x,
                y,
                x.abs().max(y.abs()),
                controls.stream_energy_tolerance_j,
                controls.relative_tolerance,
            )?);
        }
    }
    Ok(error)
}

fn audit_window<'l, 'p>(
    before: &BTreeMap<u64, ClosedGasParcel>,
    after: &BTreeMap<u64, ClosedGasParcel>,
    walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    next_walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    steps: &[FiniteVacuumChannelStep],
    duration_s: f64,
    controls: GasTransportSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<FiniteVacuumChannelWindowBalance, Error> {
    let gas = native::audit_window_borrowed(before, after,
        steps.iter().flat_map(|step| step.gas.outlets.iter()), controls, &mut *cancelled,
    ).map_err(|error| {
        if error == crate::thermochemical_gas::energy::parcel::transfer::species::network::transport::GasTransportError::Cancelled {
            Error::Cancelled
        } else {
            Error::Step { elapsed_s: 0.0, duration_s, error: FiniteChannelError::Transport(error) }
        }
    })?;
    if !walls.keys().eq(next_walls.keys()) {
        return Err(Error::Invalid(
            "adaptive vacuum window lost finite wall identities",
        ));
    }
    let mut sampled = BTreeMap::<u64, Dd>::new();
    for step in steps {
        for outlet in &step.gas.outlets {
            poll(cancelled)?;
            let wall = step
                .channel_walls
                .get(&outlet.id)
                .ok_or(Error::Invalid("adaptive vacuum window lost an outlet wall"))?;
            let sum = sampled.entry(*wall).or_insert(Dd::ZERO);
            *sum = *sum + Dd::from_f64(outlet.wall.heat_into_wall_j);
        }
    }
    let (mut stored, mut entropy, mut correction) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
    for (&id, wall) in walls {
        poll(cancelled)?;
        let next = &next_walls[&id];
        if wall.mass_kg() != next.mass_kg()
            || wall.volume_m3() != next.volume_m3()
            || wall.law().identity() != next.law().identity()
        {
            return Err(Error::Invalid(
                "adaptive vacuum changed its original finite material",
            ));
        }
        let energy = (Dd::from_f64(next.specific_energy_j_kg())
            - Dd::from_f64(wall.specific_energy_j_kg()))
            * Dd::from_f64(wall.mass_kg());
        let defect = energy - sampled.get(&id).copied().unwrap_or(Dd::ZERO);
        if !defect.to_f64().is_finite()
            || defect.to_f64().abs() > controls.exchange.energy_tolerance_j
        {
            return Err(Error::Invalid(
                "adaptive vacuum wall exceeded its fixed storage budget",
            ));
        }
        let ds = wall
            .law()
            .specific_entropy_change_j_kg_k(
                wall.specific_energy_j_kg(),
                next.specific_energy_j_kg(),
            )
            .map_err(|_| Error::Invalid("adaptive vacuum wall entropy left its original support"))?
            .ok_or(Error::Invalid(
                "adaptive vacuum wall entropy law unavailable",
            ))?;
        stored = stored + energy;
        entropy = entropy + Dd::from_f64(ds) * Dd::from_f64(wall.mass_kg());
        correction = correction + defect;
    }
    let energy_residual_j = (Dd::from_f64(gas.energy_residual_j) + correction).to_f64();
    if !stored.to_f64().is_finite()
        || !entropy.to_f64().is_finite()
        || !energy_residual_j.is_finite()
        || energy_residual_j.abs() > controls.exchange.energy_tolerance_j
    {
        return Err(Error::Invalid(
            "adaptive vacuum window exceeded its fixed total energy budget",
        ));
    }
    poll(cancelled)?;
    Ok(FiniteVacuumChannelWindowBalance {
        gas,
        wall_energy_change_j: stored.to_f64(),
        wall_entropy_change_j_k: entropy.to_f64(),
        energy_residual_j,
    })
}

/// Advance one complete gas/finite-wall blowdown window or publish no endpoint.
/// Native species withdrawals, wall-domain refusals and source qualification
/// remain authoritative at every resampled half step. Retained output is capped
/// by the native window audit's row budget, in addition to caller step budgets.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate<'l, 'p, C: FnMut() -> bool>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    paths: GasTransportPaths<'_>,
    links: &[FiniteVacuumChannelLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    step_controls: GasTransportSolve,
    controls: EffusionAdvanceControls,
    mut cancelled: C,
) -> Result<FiniteVacuumChannelAdvanceCandidate<'l, 'p>, Error> {
    let history = advance_history(
        parcels,
        walls,
        duration_s,
        controls,
        MAX_SNAPSHOT_TRANSFERS,
        |gas, wall, dt, cancelled| {
            step_candidate(
                gas,
                wall,
                paths,
                links,
                dt,
                endpoints,
                step_controls,
                cancelled,
            )
            .map(|candidate| candidate.into_parts())
        },
        |coarse, first, second, cancelled| stream_error(coarse, first, second, controls, cancelled),
        |step| {
            step.walls.len()
                + step.gas.exchange.fluxes.len()
                + step.gas.exchange.heat.len()
                + step.gas.exchange.external_heat.len()
                + step.gas.exchange.boundaries.len()
                + step.gas.exchange.parcels.len()
                + step
                    .gas
                    .outlets
                    .iter()
                    .map(|outlet| {
                        1 + outlet.escaped_species.len()
                            + outlet.transmitted.species.len()
                            + outlet.wall.reflected_species.len()
                    })
                    .sum::<usize>()
        },
        &mut cancelled,
    )?;
    let balance = audit_window(
        parcels,
        &history.parcels,
        walls,
        &history.walls,
        &history.steps,
        duration_s,
        step_controls,
        &mut cancelled,
    )?;
    poll(&mut cancelled)?;
    Ok(FiniteVacuumChannelAdvanceCandidate {
        parcels: history.parcels,
        walls: history.walls,
        step: FiniteVacuumChannelAdvanceStep {
            duration_s,
            attempted_steps: history.attempted_steps,
            accepted_intervals: history.accepted_intervals,
            rejected_intervals: history.rejected_intervals,
            maximum_accepted_error: history.maximum_accepted_error,
            accepted_steps: history.steps,
            balance,
        },
    })
}

#[cfg(test)]
mod tests;
