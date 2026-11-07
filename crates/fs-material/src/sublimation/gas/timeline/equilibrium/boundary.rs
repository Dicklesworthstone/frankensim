//! Conservative complete melting on the original native heat timeline.
//!
//! One explicitly exposed solid interface shares the original gas, molecular
//! streams, reflected heat, and all native thermal neighbors. A clipped phase
//! event retries the complete interval from its accepted start. Discarded
//! network heat never survives with a shorter molecular clock. The first fully
//! liquid native substep also bounds the search, even on a covered interface.
//!
//! This is event localization for the original first-order split, not another
//! heat or phase solver. Physical timestep refinement is still required. The
//! resulting liquid body needs an explicitly selected next interface owner.

use super::super::{
    GasTimelineCandidate, GasTimelineEquilibriumSurface, GasTimelineError,
    GasTimelineObservedError, ThermalWork, advance_deposition_observed_mode, normalize_surfaces,
};
use crate::sublimation::gas::{
    GasSurfaceBinding, GasSurfaceError, GasSurfaceSolve, SolidSurfaceState, SublimationLaw,
    accommodation::{ReflectedHeatStep, ThermalAccommodation},
    admit_surface_inventory,
    chamber::{
        GasChamberCandidate, GasChamberSolve, GasChamberSpeciesStep, GasChamberStep,
        GasChamberSurfaceStep,
        mixed::{CondensedInventory, SharpMeltingCaloric},
    },
    equilibrium::boundary::EquilibriumMeltingCandidate,
    invalid, molecular_stream_energy, phase_at, poll, snapshot_surface_at_temperature,
    species_snapshot, validate_surface,
};
use crate::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, StepConfig};
use std::collections::BTreeMap;

/// One physical duration consumed by native heat, molecules and reflection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeltingTimelineBoundary {
    pub elapsed_duration_s: f64,
    pub remaining_duration_s: f64,
    pub melting_completed: bool,
}

/// Original network and finite gas at the same accepted physical time.
#[derive(Clone)]
pub struct EquilibriumMeltingTimelineCandidate<'thermal, 'catalog> {
    exchange: GasTimelineCandidate<'thermal, 'catalog>,
    boundary: MeltingTimelineBoundary,
}
impl<'thermal, 'catalog> EquilibriumMeltingTimelineCandidate<'thermal, 'catalog> {
    #[must_use]
    pub const fn exchange(&self) -> &GasTimelineCandidate<'thermal, 'catalog> {
        &self.exchange
    }
    #[must_use]
    pub const fn boundary(&self) -> &MeltingTimelineBoundary {
        &self.boundary
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        GasTimelineCandidate<'thermal, 'catalog>,
        MeltingTimelineBoundary,
    ) {
        (self.exchange, self.boundary)
    }
}

/// Advance at most to `target_s`, returning early at the original liquid knot.
///
/// The selected active body and its scheduled birth require actual solid. All
/// native neighbors remain on the same accepted clock. Source identity, donor,
/// event ordering, cancellation and whole-window balance gates remain original.
/// Every discarded event probe spends the original native thermal work budget.
/// Thermal-only exact-knot selection spends the original thermal joule allowance
/// and appears in stored energy and residual, never as additional applied heat.
/// Active molecular exchange must reach the knot over the full native interval.
#[allow(clippy::too_many_arguments)]
pub fn advance_deposition_to_melting<'thermal, 'catalog>(
    thermal: &DepositionSimulation<'thermal>,
    surface: GasTimelineEquilibriumSurface<'_>,
    accommodation: Option<&ThermalAccommodation>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasChamberSolve,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<EquilibriumMeltingTimelineCandidate<'thermal, 'catalog>, GasTimelineError> {
    let surfaces = normalize_surfaces(thermal, &[surface], phase_config, &mut checkpoint)?;
    let accommodation = accommodation.map(|law| BTreeMap::from([(surface.cell_id, law.clone())]));
    let work = ThermalWork::new(thermal.work_budget().max_substeps);
    let (exchange, ()) = advance_deposition_observed_mode(
        thermal,
        &surfaces,
        accommodation.as_ref(),
        &(),
        gas,
        catalog,
        target_s,
        thermal_config,
        phase_config,
        Some(&work),
        true,
        true,
        checkpoint,
    )
    .map_err(|error| match error {
        GasTimelineObservedError::Timeline(error) => error,
        GasTimelineObservedError::Observer(error) => error.into(),
    })?;
    let boundary = MeltingTimelineBoundary {
        elapsed_duration_s: exchange.thermal().time_s() - thermal.time_s(),
        remaining_duration_s: target_s - exchange.thermal().time_s(),
        melting_completed: exchange
            .audit()
            .condensed_phases
            .get(&surface.cell_id)
            .is_some_and(|phase| phase.final_inventory.solid_mass_kg == 0.0),
    };
    Ok(EquilibriumMeltingTimelineCandidate { exchange, boundary })
}

// Normalize already accepted reports only. Neither conversion recomputes
// transport nor reconstructs gas; the legacy chamber keeps its exact reports.
pub(in crate::sublimation::gas::timeline) struct PhaseInterval {
    pub surfaces: Vec<GasChamberSurfaceStep>,
    pub gas: ClosedGasParcel,
    pub report: GasChamberStep,
    pub reflection: BTreeMap<u64, ReflectedHeatStep>,
    pub gas_residual_j: f64,
    pub gas_iterations: usize,
}
impl PhaseInterval {
    pub(in crate::sublimation::gas::timeline) fn from_chamber(
        candidate: GasChamberCandidate<'_>,
        reflection: BTreeMap<u64, ReflectedHeatStep>,
    ) -> Self {
        let gas_residual_j = candidate.gas_exchange().step().energy_residual_j;
        let gas_iterations = candidate.gas_exchange().step().iterations;
        let (surfaces, gas, _, report) = candidate.into_parts();
        Self {
            surfaces,
            gas,
            report,
            reflection,
            gas_residual_j,
            gas_iterations,
        }
    }

    pub(in crate::sublimation::gas::timeline) fn from_scalar(
        id: u64,
        binding: &GasSurfaceBinding,
        initial: SolidSurfaceState,
        original_gas: &ClosedGasParcel,
        candidate: EquilibriumMeltingCandidate<'_, '_>,
    ) -> Result<Self, GasSurfaceError> {
        let (endpoint, boundary) = candidate.into_parts();
        let (exchange, _, reflection) = endpoint.into_parts();
        let gas_step = exchange.gas_exchange().step();
        let gas_residual_j = gas_step.energy_residual_j;
        let gas_iterations = gas_step.iterations;
        // These are the ORIGINAL accepted species energy transfers, including
        // carrier reflection. Retain the volatile account even when stationary.
        let mut species_energy: BTreeMap<_, _> = gas_step
            .species
            .iter()
            .map(|row| (row.species.clone(), row.energy_into_gas_j))
            .collect();
        species_energy.entry(binding.species.clone()).or_insert(0.0);
        let (solid, gas, _, step) = exchange.into_parts();
        let mut species = Vec::with_capacity(species_energy.len());
        for (id, energy_into_gas_j) in species_energy {
            let volatile = id == binding.species;
            let old_mass = original_gas.species_mass_kg(&id)?;
            let new_mass = gas.species_mass_kg(&id)?;
            species.push(GasChamberSpeciesStep {
                species: id,
                outgoing_mass_kg: if volatile { step.outgoing_mass_kg } else { 0.0 },
                incoming_mass_kg: if volatile { step.incoming_mass_kg } else { 0.0 },
                transferred_mass_kg: if volatile {
                    step.transferred_mass_kg
                } else {
                    0.0
                },
                energy_into_gas_j,
                vapor_mass_kg: [old_mass, new_mass],
                mass_residual_kg: (new_mass - old_mass)
                    + if volatile {
                        solid.mass_kg - initial.mass_kg
                    } else {
                        0.0
                    },
            });
        }
        let surfaces = vec![GasChamberSurfaceStep {
            id,
            solid,
            outgoing_mass_kg: step.outgoing_mass_kg,
            incoming_mass_kg: step.incoming_mass_kg,
            transferred_mass_kg: step.transferred_mass_kg,
            surface_temperature_k: step.surface_temperature_k,
            outgoing_specific_energy_j_kg: step.outgoing_specific_energy_j_kg,
            incoming_specific_energy_j_kg: step.incoming_specific_energy_j_kg,
            energy_into_gas_j: step.energy_into_gas_j,
            supplied_heat_j: 0.0,
            solid_energy_change_j: step.solid_energy_change_j,
            energy_residual_j: boundary.condensed_energy_residual_j,
            recession_m: step.recession_m,
        }];
        let report = GasChamberStep {
            outgoing_mass_kg: step.outgoing_mass_kg,
            incoming_mass_kg: step.incoming_mass_kg,
            transferred_mass_kg: step.transferred_mass_kg,
            energy_into_gas_j: step.energy_into_gas_j,
            supplied_heat_j: 0.0,
            solid_energy_change_j: step.solid_energy_change_j,
            gas_energy_change_j: step.gas_energy_change_j,
            energy_residual_j: step.energy_residual_j,
            mass_residual_kg: step.mass_residual_kg,
            gas_temperature_k: step.gas_temperature_k,
            species,
        };
        let reflection = reflection.map_or_else(BTreeMap::new, |row| BTreeMap::from([(id, row)]));
        Ok(Self {
            surfaces,
            gas,
            report,
            reflection,
            gas_residual_j,
            gas_iterations,
        })
    }
}

// Query the limiting ORIGINAL solid constituent at its fusion temperature.
// The retained donor is the actual solid at the start of the WHOLE interval,
// never a fabricated post-heat epsilon seed. This separate terminal operation
// is admitted only when net mass AND stream/reflected energy are exactly zero:
// the resulting identity phase update commutes with native heat at the fusion
// endpoint. It is not an ordinary post-heat positive-solid-seed phase step.
// A unit exposure first classifies rates without a tiny physical area hiding
// nonzero transfer. The accepted physical exposure retains gross turnover.
#[allow(clippy::too_many_arguments)]
pub(in crate::sublimation::gas::timeline) fn inactive_at_melting(
    id: u64,
    source: SharpMeltingCaloric<'_>,
    law: &dyn SublimationLaw,
    binding: &GasSurfaceBinding,
    original: SolidSurfaceState,
    gas: &ClosedGasParcel,
    catalog: &ThermochemicalGasModel,
    accommodation: Option<&ThermalAccommodation>,
    duration_s: f64,
    controls: GasSurfaceSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Option<PhaseInterval>, GasSurfaceError> {
    poll(&mut cancelled)?;
    validate_surface(binding, original)?;
    let initial = CondensedInventory::new(source, original)?;
    if initial.solid_mass_kg <= 0.0 {
        return Err(invalid("melting event requires original solid"));
    }
    let terminal_h = source.liquid_at_melting_h_j_kg();
    let terminal = source.resolve(terminal_h)?;
    let context = species_snapshot(gas, catalog, &binding.species)?;
    let properties = phase_at(law, terminal.temperature_k(), context.molar_mass)?;
    source.admit_density(terminal, properties.solid_density_kg_m3)?;
    let mut evaluate = |area, dt| {
        poll(&mut cancelled)?;
        let snapshot = snapshot_surface_at_temperature(
            law,
            binding,
            SolidSurfaceState {
                mass_kg: initial.solid_mass_kg,
                specific_enthalpy_j_kg: source.solid_enthalpy(terminal),
                area_m2: area,
            },
            context,
            catalog,
            terminal.temperature_k(),
            dt,
            controls,
        )?;
        let reflected = match accommodation {
            Some(law) => Some(law.evaluate(
                binding,
                gas,
                terminal.temperature_k(),
                area,
                dt,
                &mut cancelled,
            )?),
            None => None,
        };
        let energy = molecular_stream_energy(
            snapshot,
            reflected.as_ref().map_or(0.0, |row| row.heat_into_gas_j),
        )?;
        Ok::<_, GasSurfaceError>((snapshot, reflected, energy))
    };
    let (probe, _, energy) = evaluate(if original.area_m2 == 0.0 { 0.0 } else { 1.0 }, 1.0)?;
    if probe.transferred_mass_kg != 0.0 || energy.relative_energy_into_gas_j != 0.0 {
        return Ok(None);
    }
    let (snapshot, reflection, energy) = evaluate(original.area_m2, duration_s)?;
    if snapshot.transferred_mass_kg != 0.0 || energy.relative_energy_into_gas_j != 0.0 {
        return Err(invalid(
            "inactive melting endpoint changed its physical stream classification",
        ));
    }
    // The original gross solid/OLD vapor gates and their exact stationary-
    // equilibrium exception still own this turnover; no incoming stream funds
    // outgoing material and no new gas parcel is manufactured.
    admit_surface_inventory(snapshot, duration_s, controls)?;
    let surface = GasChamberSurfaceStep {
        id,
        solid: SolidSurfaceState {
            specific_enthalpy_j_kg: terminal_h,
            ..original
        },
        outgoing_mass_kg: snapshot.outgoing_mass_kg,
        incoming_mass_kg: snapshot.incoming_mass_kg,
        transferred_mass_kg: 0.0,
        surface_temperature_k: [terminal.temperature_k(); 2],
        outgoing_specific_energy_j_kg: snapshot.outgoing_specific_energy_j_kg,
        incoming_specific_energy_j_kg: snapshot.incoming_specific_energy_j_kg,
        energy_into_gas_j: 0.0,
        supplied_heat_j: 0.0,
        solid_energy_change_j: 0.0,
        energy_residual_j: 0.0,
        recession_m: 0.0,
    };
    let report = GasChamberStep {
        outgoing_mass_kg: snapshot.outgoing_mass_kg,
        incoming_mass_kg: snapshot.incoming_mass_kg,
        transferred_mass_kg: 0.0,
        energy_into_gas_j: 0.0,
        supplied_heat_j: 0.0,
        solid_energy_change_j: 0.0,
        gas_energy_change_j: 0.0,
        energy_residual_j: 0.0,
        mass_residual_kg: 0.0,
        gas_temperature_k: [context.temperature_k; 2],
        species: vec![GasChamberSpeciesStep {
            species: binding.species.clone(),
            outgoing_mass_kg: snapshot.outgoing_mass_kg,
            incoming_mass_kg: snapshot.incoming_mass_kg,
            transferred_mass_kg: 0.0,
            energy_into_gas_j: 0.0,
            vapor_mass_kg: [context.vapor_mass_kg; 2],
            mass_residual_kg: 0.0,
        }],
    };
    let reflection = reflection.map_or_else(BTreeMap::new, |row| BTreeMap::from([(id, row)]));
    poll(&mut cancelled)?;
    Ok(Some(PhaseInterval {
        surfaces: vec![surface],
        gas: gas.clone(),
        report,
        reflection,
        gas_residual_j: 0.0,
        gas_iterations: 0,
    }))
}

// A bounded root bracket for the existing physical owners. Residuals are the
// original source-relative extensive headroom to complete melting.
pub(in crate::sublimation::gas::timeline) struct MeltingSearch {
    lower: (f64, f64),
    upper: Option<(f64, f64)>,
    trials: usize,
}
impl MeltingSearch {
    pub(in crate::sublimation::gas::timeline) fn new(start_s: f64, headroom_j: f64) -> Self {
        Self {
            lower: (start_s, headroom_j),
            upper: None,
            trials: 0,
        }
    }
    pub(in crate::sublimation::gas::timeline) const fn is_bracketed(&self) -> bool {
        self.upper.is_some()
    }
    pub(in crate::sublimation::gas::timeline) fn observe(
        &mut self,
        time_s: f64,
        residual_j: f64,
    ) -> Result<Option<f64>, DepositionError> {
        if !time_s.is_finite() || !residual_j.is_finite() || self.lower.1 <= 0.0 {
            return Err(DepositionError::Arithmetic);
        }
        if residual_j <= 0.0 {
            self.upper = Some((time_s, residual_j));
        } else if self.upper.is_some() {
            self.lower = (time_s, residual_j);
        } else {
            return Ok(None);
        }
        self.trials += 1;
        if self.trials > 64 {
            return Err(DepositionError::Budget("melting event localization trials"));
        }
        let (hi, r_hi) = self.upper.ok_or(DepositionError::Arithmetic)?;
        let (lo, r_lo) = self.lower;
        if hi <= lo {
            return Err(DepositionError::Arithmetic);
        }
        let scale = r_lo.max(r_hi.abs());
        let fraction = ((r_lo / scale) / (r_lo / scale + r_hi.abs() / scale)).clamp(0.05, 0.95);
        let mut next = lo + (hi - lo) * fraction;
        if next <= lo || next >= hi {
            next = lo + (hi - lo) * 0.5;
        }
        if !next.is_finite() || next <= lo || next >= hi {
            return Err(DepositionError::InvalidInput(
                "melting event has no representable common heat/phase interval",
            ));
        }
        Ok(Some(next))
    }
}

#[cfg(test)]
#[path = "boundary/tests.rs"]
mod tests;
