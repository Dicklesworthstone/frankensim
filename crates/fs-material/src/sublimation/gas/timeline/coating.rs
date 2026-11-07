//! Supported deposition, original heat transport and one evolving gas clock.
//!
//! The coating owner retains substrate mass, its original caloric chart and
//! reference geometry. Deposited mass has its own h(T), may start at zero and
//! participates in every native heat step. The original gas parcel owns all
//! vapor mass and internal energy; coating reservoirs are prescribed diagnostic
//! observations, never copies of that finite inventory.
//!
//! Heat precedes one joint molecular solve per event-clipped interval. Births at
//! its endpoint cannot collect molecules from the preceding interval, and end
//! contact changes cannot rewrite its exposure. This is first-order splitting:
//! refine max_dt_s for accuracy and for explicit gross-inventory refusals.
//! Thin coatings share the substrate temperature and authored reference area.
//! Nucleation, internal film resistance, moving geometry and liquid exchange are
//! not inferred from a bare surface or from the depositing species' name.

/// Error-controlled complete heat/coating/gas histories, including bare surfaces.
#[path = "coating/adaptive.rs"]
pub mod adaptive;

use std::cell::RefCell;
use std::collections::BTreeMap;

use super::super::accommodation::ThermalAccommodation;
use super::super::chamber::coating::{
    GasSupportedChamberSolve, GasSupportedSurface, GasSupportedSurfaceStep,
    advance_supported_chamber, advance_supported_chamber_with_films,
};
use super::super::{GasSurfaceBinding, GasSurfaceError, compatibility, phase_at, species_snapshot};
use super::{GasCellAudit, GasTimelineAudit, HeatTotals, Sum, ThermalWork, add, finite};
use crate::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_thermal::enthalpy::sublimation::{
    ThermalSublimationError,
    coating::{CoatingSimulation, external::CoatingCellExchange},
};
use fs_thermal::enthalpy::{CellId, DepositionError, EnergyAudit, Event, StepConfig};
use fs_thermal::sublimation::{SublimationError, VaporReservoir, coating::SupportedSurfaceState};

/// Bind the original registered coating to its original gas species/calorics.
/// Every configured coating must occur exactly once, including future births.
#[derive(Clone, Copy)]
pub struct GasCoatingTimelineSurface<'a> {
    pub cell_id: CellId,
    pub binding: &'a GasSurfaceBinding,
}

#[derive(Clone, Debug, PartialEq)]
pub enum GasCoatingTimelineError {
    Thermal(ThermalSublimationError),
    Exchange(GasSurfaceError),
    MassBalance { residual_kg: f64, tolerance_kg: f64 },
}
impl From<DepositionError> for GasCoatingTimelineError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error.into())
    }
}
impl From<ThermalSublimationError> for GasCoatingTimelineError {
    fn from(error: ThermalSublimationError) -> Self {
        Self::Thermal(error)
    }
}
impl From<GasSurfaceError> for GasCoatingTimelineError {
    fn from(error: GasSurfaceError) -> Self {
        Self::Exchange(error)
    }
}
impl core::fmt::Display for GasCoatingTimelineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "supported-coating/evolving-gas advance refused: {self:?}"
        )
    }
}
impl std::error::Error for GasCoatingTimelineError {}

/// Accept the original substrate/coating and gas endpoints together.
#[derive(Clone)]
pub struct GasCoatingTimelineCandidate<'thermal, 'catalog> {
    coating: CoatingSimulation<'thermal>,
    gas: ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    audit: GasTimelineAudit,
    // Source-relative defects, retained separately so adaptive histories do not
    // subtract large formation/reference energies or hide opposing residuals.
    solid_stream_energy_residual_j: f64,
    gas_stream_energy_residual_j: f64,
}
impl<'thermal, 'catalog> GasCoatingTimelineCandidate<'thermal, 'catalog> {
    #[must_use]
    pub const fn coating(&self) -> &CoatingSimulation<'thermal> {
        &self.coating
    }
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    #[must_use]
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    /// `solid_energy_change_j` includes BOTH substrate and coating storage.
    /// Fixed gas-datum conversions apply only to deposit mass and deposit births.
    /// `condensed_phases` is empty: no solid/liquid fraction is manufactured.
    #[must_use]
    pub const fn audit(&self) -> &GasTimelineAudit {
        &self.audit
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        CoatingSimulation<'thermal>,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasTimelineAudit,
    ) {
        (self.coating, self.gas, self.catalog, self.audit)
    }
}

// Admit all active endpoints in one original-owner batch, including unchanged
// newborns whose diagnostic observation must follow the accepted shared gas.
fn accept_endpoints(
    coating: &mut CoatingSimulation<'_>,
    bindings: &BTreeMap<CellId, &GasSurfaceBinding>,
    endpoints: &[GasSupportedSurfaceStep],
    gas: &ClosedGasParcel,
    catalog: &ThermochemicalGasModel,
    latent_tolerance_j_kg: f64,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<(), GasCoatingTimelineError> {
    let mut observations = BTreeMap::new();
    let mut changes = Vec::with_capacity(coating.coatings().len());
    let mut endpoints = endpoints.iter().peekable();
    for (&id, state) in coating.coatings() {
        checkpoint()?;
        let binding = bindings[&id];
        let source = match observations.get(&binding.species) {
            Some(value) => *value,
            None => {
                let source = species_snapshot(gas, catalog, &binding.species)?;
                observations.insert(binding.species.clone(), source);
                source
            }
        };
        let substrate = &coating.cells()[&id];
        let replacement = if endpoints.peek().is_some_and(|step| step.id == id.0) {
            endpoints.next().unwrap().state
        } else {
            // Newborn and zero-window coatings did not enter a molecular
            // solve. Admit their original phase/species/datum here without
            // inventing a prebirth stream or another gas reconstruction.
            let model = &coating.coating_surfaces()[&id];
            let temperature =
                coating.temperature_at_specific_enthalpy_k(id, substrate.specific_h_j_kg)?;
            let phase = phase_at(model.law, temperature, source.molar_mass)?;
            let deposit_h = (model.deposit_enthalpy)(temperature).map_err(GasSurfaceError::from)?;
            let vapor_h = catalog
                .species_specific_enthalpy_j_kg(&binding.species, temperature)
                .map_err(GasSurfaceError::from)?;
            compatibility(
                vapor_h,
                deposit_h,
                binding.solid_enthalpy_offset_j_kg,
                phase.latent_heat_j_kg,
                temperature,
                latent_tolerance_j_kg,
            )?;
            SupportedSurfaceState {
                substrate_mass_kg: substrate.mass_kg,
                substrate_specific_enthalpy_j_kg: substrate.specific_h_j_kg,
                deposit_mass_kg: state.mass_kg,
                area_m2: 0.0,
            }
        };
        changes.push(CoatingCellExchange {
            expected_substrate: substrate.clone(),
            expected_coating: *state,
            interval_exposed_area_m2: replacement.area_m2,
            replacement,
            replacement_reservoir: VaporReservoir::Prescribed {
                partial_pressure_pa: source.partial_pressure_pa,
                temperature_k: source.temperature_k,
            },
        });
    }
    if endpoints.next().is_some() {
        return Err(DepositionError::InvalidInput("supported endpoint identity changed").into());
    }
    coating.apply_coating_exchanges(&changes, checkpoint)?;
    Ok(())
}

/// Advance supported deposition without replacing either physical owner.
///
/// The selected binding coefficient must equal the registered coating's
/// molecular coefficient. Its optional prescribed film-transfer resistance is
/// retained in the same gross molecular streams, so finite films limit both
/// mass and crossing-molecule energy consistently. Optional reflected-carrier
/// heat uses the same gas and supported endpoint.
/// One original thermal substep allowance covers the complete requested window.
/// Refusal/cancellation discards all candidate heat, births, phase and gas work.
#[allow(clippy::too_many_arguments)]
pub fn advance_coatings_with_gas<'thermal, 'catalog>(
    coating: &CoatingSimulation<'thermal>,
    surfaces: &[GasCoatingTimelineSurface<'_>],
    accommodation: Option<&BTreeMap<CellId, ThermalAccommodation>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    controls: GasSupportedChamberSolve,
    checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<GasCoatingTimelineCandidate<'thermal, 'catalog>, GasCoatingTimelineError> {
    advance_coatings_metered(coating, surfaces, accommodation, gas, catalog,
        target_s, thermal_config, controls, None, checkpoint)
}

// Share the original physical march with adaptive coarse/fine probes. Completed
// heat work is charged BEFORE a molecular inventory refusal can request retry.
#[allow(clippy::too_many_arguments)]
fn advance_coatings_metered<'thermal, 'catalog>(
    coating: &CoatingSimulation<'thermal>,
    surfaces: &[GasCoatingTimelineSurface<'_>],
    accommodation: Option<&BTreeMap<CellId, ThermalAccommodation>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    controls: GasSupportedChamberSolve,
    work: Option<&ThermalWork>,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<GasCoatingTimelineCandidate<'thermal, 'catalog>, GasCoatingTimelineError> {
    checkpoint()?;
    if !target_s.is_finite()
        || target_s < coating.time_s()
        || !thermal_config.max_dt_s.is_finite()
        || thermal_config.max_dt_s <= 0.0
        || !thermal_config.energy_tolerance_j.is_finite()
        || thermal_config.energy_tolerance_j < 0.0
        || surfaces.is_empty()
        || surfaces.len() != coating.coating_surfaces().len()
    {
        return Err(DepositionError::InvalidInput(
            "finite forward coating window and complete surface selection required",
        )
        .into());
    }
    let budget = coating.work_budget();
    if surfaces.len() > budget.max_cells || surfaces.len() > controls.chamber.maximum_surfaces {
        return Err(DepositionError::Budget("evolving-gas coating selection").into());
    }
    let mut bindings = BTreeMap::new();
    for surface in surfaces {
        checkpoint()?;
        let model = coating
            .coating_surfaces()
            .get(&surface.cell_id)
            .ok_or(DepositionError::InvalidInput("unregistered gas coating"))?;
        if bindings.insert(surface.cell_id, surface.binding).is_some()
            || !surface.binding.solid_enthalpy_offset_j_kg.is_finite()
            || model.kinetics.accommodation_coefficient != surface.binding.accommodation_coefficient
        {
            return Err(DepositionError::InvalidInput(
                "duplicate coating, invalid datum or incompatible molecular kinetics",
            )
            .into());
        }
        species_snapshot(gas, catalog, &surface.binding.species)?;
    }
    if let Some(laws) = accommodation {
        if laws.len() != bindings.len() || bindings.keys().any(|id| !laws.contains_key(id)) {
            return Err(
                DepositionError::InvalidInput("coating accommodation roster mismatch").into(),
            );
        }
    }
    // Preserve original birth inventories for independent whole-window checks.
    let mut births = BTreeMap::new();
    for (index, event) in coating.pending_events().iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        if event.time_s <= target_s {
            if let Event::Birth(cell) = &event.event {
                births.insert(cell.id, cell);
            }
        }
    }
    // Admit original gas/catalog and controls even in a zero-length window or
    // when every selected coating is scheduled strictly after the endpoint.
    let mut refusal = None;
    let admission =
        advance_supported_chamber(
            &[],
            None,
            gas,
            catalog,
            0.0,
            controls,
            || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    refusal = Some(error);
                    true
                }
            },
        );
    if let Some(error) = refusal {
        return Err(error.into());
    }
    admission?;
    let (mut next, first_heat) = coating.advance_heat_to_limited(
        coating.time_s(),
        thermal_config,
        controls.storage,
        work.map_or(budget.max_substeps, ThermalWork::remaining),
        &mut checkpoint,
    )?;
    if let Some(work) = work { work.charge_steps(first_heat.substeps)?; }
    let mut heat = HeatTotals::default();
    heat.add(first_heat)?;
    let mut next_gas = gas.clone();
    accept_endpoints(
        &mut next,
        &bindings,
        &[],
        &next_gas,
        catalog,
        controls.chamber.surface.latent_compatibility_tolerance_j_kg,
        &mut checkpoint,
    )?;
    let mut accounts: BTreeMap<_, _> = bindings
        .keys()
        .map(|&id| (id, GasCellAudit::default()))
        .collect();
    let mut active_accommodation = BTreeMap::new();
    let mut outgoing = Sum::default();
    let mut incoming = Sum::default();
    let mut transferred = Sum::default();
    let mut stream_energy = Sum::default();
    let mut reflected_energy = Sum::default();
    let mut phase_storage = Sum::default();
    let mut gas_storage = Sum::default();
    let mut gas_defect = Sum::default();
    let mut phase_defect = Sum::default();
    let mut phase_steps = 0usize;
    let mut gas_iterations = 0usize;
    let mut last_step = None;
    while next.time_s() < target_s {
        checkpoint()?;
        let start = next.time_s();
        let mut active = Vec::with_capacity(bindings.len());
        for &id in bindings.keys() {
            checkpoint()?;
            if next.coatings().contains_key(&id) {
                active.push(id);
                if let Some(laws) = accommodation {
                    active_accommodation
                        .entry(id.0)
                        .or_insert_with(|| laws[&id].clone());
                }
            }
        }
        let areas = next.exposed_areas_m2(&active, &mut checkpoint)?;
        let mut stop = if next.cells().is_empty() {
            target_s
        } else {
            (start + thermal_config.max_dt_s).min(target_s)
        };
        if let Some(event) = next.pending_events().first() {
            stop = stop.min(event.time_s);
        }
        if !stop.is_finite() || stop <= start {
            return Err(DepositionError::Arithmetic.into());
        }
        let remaining =
            budget
                .max_substeps
                .checked_sub(heat.substeps)
                .ok_or(DepositionError::Budget(
                    "evolving-gas coating thermal substeps",
                ))?;
        let (heated, report) = next.advance_heat_to_limited(
            stop,
            thermal_config,
            controls.storage,
            work.map_or(remaining, |meter| remaining.min(meter.remaining())),
            &mut checkpoint,
        )?;
        if let Some(work) = work { work.charge_steps(report.substeps)?; }
        next = heated;
        heat.add(report)?;
        if active.is_empty() {
            accept_endpoints(
                &mut next,
                &bindings,
                &[],
                &next_gas,
                catalog,
                controls.chamber.surface.latent_compatibility_tolerance_j_kg,
                &mut checkpoint,
            )?;
            continue;
        }
        let chart_error = RefCell::new(None);
        let calorics: Vec<_> = active
            .iter()
            .map(|&id| {
                let owner = &next;
                let errors = &chart_error;
                move |h| {
                    owner
                        .temperature_at_specific_enthalpy_k(id, h)
                        .map_err(|error| {
                            *errors.borrow_mut() = Some(error);
                            SublimationError::Invalid(
                                "original substrate chart refused molecular endpoint",
                            )
                        })
                }
            })
            .collect();
        let inputs: Vec<_> = active
            .iter()
            .zip(&calorics)
            .map(|(&id, caloric)| {
                let substrate = &next.cells()[&id];
                let model = &next.coating_surfaces()[&id];
                GasSupportedSurface {
                    id: id.0,
                    law: model.law,
                    binding: bindings[&id],
                    state: SupportedSurfaceState {
                        substrate_mass_kg: substrate.mass_kg,
                        substrate_specific_enthalpy_j_kg: substrate.specific_h_j_kg,
                        deposit_mass_kg: next.coatings()[&id].mass_kg,
                        area_m2: areas[&id],
                    },
                    substrate_temperature: caloric,
                    deposit_enthalpy: model.deposit_enthalpy,
                }
            })
            .collect();
        let film_coefficients: BTreeMap<_, _> = active
            .iter()
            .filter_map(|id| {
                next.coating_surfaces()[id]
                    .kinetics
                    .mass_transfer_coefficient_m_s
                    .map(|coefficient| (id.0, coefficient))
            })
            .collect();
        let mut refusal = None;
        let exchange = advance_supported_chamber_with_films(
            &inputs,
            accommodation.map(|_| &active_accommodation),
            &film_coefficients,
            &next_gas,
            catalog,
            stop - start,
            controls,
            || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    refusal = Some(error);
                    true
                }
            },
        );
        drop(inputs);
        drop(calorics);
        if let Some(error) = refusal {
            return Err(error.into());
        }
        let candidate = match exchange {
            Ok(value) => value,
            Err(error) => {
                return Err(match chart_error.into_inner() {
                    Some(error) => error.into(),
                    None => error.into(),
                });
            }
        };
        if candidate.surfaces().len() != active.len()
            || candidate
                .surfaces()
                .iter()
                .zip(&active)
                .any(|(step, id)| step.id != id.0)
        {
            return Err(DepositionError::InvalidInput("supported endpoint roster changed").into());
        }
        let iterations = candidate.gas_exchange().step().iterations;
        gas_defect.add(candidate.gas_exchange().step().energy_residual_j)?;
        accept_endpoints(
            &mut next,
            &bindings,
            candidate.surfaces(),
            candidate.gas(),
            catalog,
            controls.chamber.surface.latent_compatibility_tolerance_j_kg,
            &mut checkpoint,
        )?;
        let (endpoints, accepted_gas, _, report) = candidate.into_parts();
        next_gas = accepted_gas;
        outgoing.add(report.outgoing_mass_kg)?;
        incoming.add(report.incoming_mass_kg)?;
        transferred.add(report.transferred_mass_kg)?;
        stream_energy.add(report.energy_into_gas_j)?;
        phase_storage.add(report.solid_energy_change_j)?;
        gas_storage.add(report.gas_energy_change_j)?;
        phase_steps = phase_steps
            .checked_add(1)
            .ok_or(DepositionError::Arithmetic)?;
        gas_iterations = gas_iterations
            .checked_add(iterations)
            .ok_or(DepositionError::Arithmetic)?;
        for endpoint in endpoints {
            checkpoint()?;
            let exchange = endpoint.exchange;
            phase_defect.add(exchange.energy_residual_j)?;
            reflected_energy.add(endpoint.reflected_heat_into_gas_j)?;
            let account = accounts.get_mut(&CellId(endpoint.id)).unwrap();
            account.outgoing_mass_kg = add(account.outgoing_mass_kg, exchange.outgoing_mass_kg)?;
            account.incoming_mass_kg = add(account.incoming_mass_kg, exchange.incoming_mass_kg)?;
            account.transferred_mass_kg =
                add(account.transferred_mass_kg, exchange.transferred_mass_kg)?;
            account.energy_into_gas_j = add(account.energy_into_gas_j, exchange.energy_into_gas_j)?;
            account.reflected_heat_into_gas_j = add(
                account.reflected_heat_into_gas_j,
                endpoint.reflected_heat_into_gas_j,
            )?;
            account.recession_m = add(account.recession_m, exchange.recession_m)?;
            account.phase_steps = account
                .phase_steps
                .checked_add(1)
                .ok_or(DepositionError::Arithmetic)?;
        }
        last_step = Some(report);
    }
    let thermal_audit: EnergyAudit = heat.audit()?;
    let mut deposit_change = Sum::default();
    let mut species_change: BTreeMap<_, Sum> = BTreeMap::new();
    let mut birth_offsets = Sum::default();
    for (&id, substrate) in next.cells() {
        checkpoint()?;
        let original = coating
            .cells()
            .get(&id)
            .or_else(|| births.get(&id).copied())
            .ok_or(DepositionError::InvalidInput(
                "substrate endpoint has no original or birth",
            ))?;
        if original.mass_kg != substrate.mass_kg {
            return Err(DepositionError::InvalidInput(
                "supported deposition changed substrate mass",
            )
            .into());
        }
        if let Some(deposit) = next.coatings().get(&id) {
            let binding = bindings[&id];
            let old_mass = match coating.coatings().get(&id) {
                Some(old) => old.mass_kg,
                None => {
                    let mass = coating.coating_surfaces()[&id].initial_mass_kg;
                    birth_offsets.add(mass * binding.solid_enthalpy_offset_j_kg)?;
                    mass
                }
            };
            let dm = finite(deposit.mass_kg - old_mass)?;
            deposit_change.add(dm)?;
            species_change
                .entry(binding.species.clone())
                .or_default()
                .add(dm)?;
        }
    }
    let mass_residual = add(deposit_change.value()?, next_gas.mass_kg() - gas.mass_kg())?;
    let tolerance = controls.chamber.surface.gas.mass_tolerance_kg;
    let mut largest = mass_residual.abs();
    for source in catalog.thermodynamics().components() {
        checkpoint()?;
        let id = source.species();
        let condensed = species_change.get(id).map_or(Ok(0.0), Sum::value)?;
        let gas_change = next_gas
            .species_mass_kg(id)
            .map_err(GasSurfaceError::from)?
            - gas.species_mass_kg(id).map_err(GasSurfaceError::from)?;
        largest = largest.max(finite(condensed + gas_change)?.abs());
    }
    if largest > tolerance {
        return Err(GasCoatingTimelineError::MassBalance {
            residual_kg: largest,
            tolerance_kg: tolerance,
        });
    }
    // Original local caloric increments cancel independent substrate datums
    // and transferred formation energy before whole-window accumulation.
    let supported_defect = add(thermal_audit.residual_j, phase_defect.value()?)?;
    let gas_defect = gas_defect.value()?;
    let energy_residual = add(supported_defect, gas_defect)?;
    let energy_tolerance = thermal_config
        .energy_tolerance_j
        .max(controls.chamber.surface.gas.energy_tolerance_j);
    for residual in [supported_defect, gas_defect, energy_residual] {
        if residual.abs() > energy_tolerance {
            return Err(GasSurfaceError::EnergyBalance {
                residual_j: residual,
                tolerance_j: energy_tolerance,
            }
            .into());
        }
    }
    let incoming_energy = add(thermal_audit.incoming_energy_j, birth_offsets.value()?)?;
    let mut storage = Sum::default();
    storage.add(thermal_audit.stored_energy_change_j)?;
    storage.add(phase_storage.value()?)?;
    storage.add(birth_offsets.value()?)?;
    let audit = GasTimelineAudit {
        thermal: thermal_audit,
        cells: accounts,
        condensed_phases: BTreeMap::new(),
        outgoing_mass_kg: outgoing.value()?,
        incoming_mass_kg: incoming.value()?,
        transferred_mass_kg: transferred.value()?,
        energy_into_gas_j: stream_energy.value()?,
        reflected_heat_into_gas_j: reflected_energy.value()?,
        solid_energy_change_j: storage.value()?,
        gas_energy_change_j: gas_storage.value()?,
        incoming_energy_j: incoming_energy,
        mass_residual_kg: mass_residual,
        energy_residual_j: energy_residual,
        phase_steps,
        gas_iterations,
        last_step,
    };
    checkpoint()?;
    Ok(GasCoatingTimelineCandidate {
        coating: next,
        gas: next_gas,
        catalog,
        audit,
        solid_stream_energy_residual_j: supported_defect,
        gas_stream_energy_residual_j: gas_defect,
    })
}

#[cfg(test)]
#[path = "coating/tests.rs"]
mod tests;
