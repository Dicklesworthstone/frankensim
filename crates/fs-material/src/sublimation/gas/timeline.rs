//! Original thermal transport and a jointly evolving gas on one timeline.
//!
//! Every selected cell is a caller-admitted pure solid, on its original caloric
//! chart and a declared enthalpy interval. The existing heat observer checks that
//! branch after each accepted thermal substep and birth. One joint molecular
//! step follows heat over each event-clipped interval; refine max_dt_s for this
//! first-order splitting and the explicit molecular rates. Gross-inventory
//! refusals require a shorter physical step, never independently clipped streams.
//!
//! Optional per-cell thermal accommodation reaches the SAME joint molecular
//! endpoint. Both sticking and reflected streams use the preceding interval's
//! exposure; neither can act before birth or use endpoint contact coverage
//! retroactively. Coefficients are activated once per cell and reused thereafter.
//!
//! Gas temperature, species inventories and internal energy come from the
//! original gas parcel. Solid datum offsets are fixed reference conversions,
//! included in birth and storage accounts; latent energy is never added twice.
//! Reflected heat is included once in stream energy and separately reported as
//! a subaccount. Reference geometry is retained, with recession per surface.

use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};

use super::accommodation::ThermalAccommodation;
use super::chamber::{
    GasChamberSolve, GasChamberStep, GasChamberSurface, advance_chamber_with_accommodation,
    advance_chamber_with_gas,
};
use super::{GasSurfaceBinding, GasSurfaceError};
use crate::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit,
    Event, StepConfig, sublimation::external::SolidCellExchange,
};
use fs_thermal::sublimation::{SolidSurfaceState, SublimationError, SublimationLaw};

/// Original cell identity and explicitly admitted pure-solid phase binding.
#[derive(Clone, Copy)]
pub struct GasTimelineSurface<'a> {
    pub cell_id: CellId,
    pub law: &'a dyn SublimationLaw,
    pub binding: &'a GasSurfaceBinding,
    /// Inclusive pure-solid branch in the ORIGINAL thermal law's h datum.
    pub pure_solid_enthalpy_range_j_kg: [f64; 2],
}

#[derive(Clone, Debug, PartialEq)]
pub enum GasTimelineError {
    Thermal(DepositionError),
    Exchange(GasSurfaceError),
    MassBalance { residual_kg: f64, tolerance_kg: f64 },
}
impl From<DepositionError> for GasTimelineError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<GasSurfaceError> for GasTimelineError {
    fn from(error: GasSurfaceError) -> Self {
        Self::Exchange(error)
    }
}
impl core::fmt::Display for GasTimelineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal/evolving-gas advance refused: {self:?}")
    }
}
impl std::error::Error for GasTimelineError {}

/// Local accepted molecular accounts; negative recession denotes growth.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct GasCellAudit {
    pub outgoing_mass_kg: f64,
    pub incoming_mass_kg: f64,
    pub transferred_mass_kg: f64,
    pub energy_into_gas_j: f64,
    /// Reflected heat INCLUDED in energy_into_gas_j; zero on the specular path.
    pub reflected_heat_into_gas_j: f64,
    pub recession_m: f64,
    pub phase_steps: usize,
}

/// Actual whole-window conservation in the original gas energy datum.
#[derive(Clone, Debug, PartialEq)]
pub struct GasTimelineAudit {
    pub thermal: EnergyAudit,
    pub cells: BTreeMap<CellId, GasCellAudit>,
    pub outgoing_mass_kg: f64,
    pub incoming_mass_kg: f64,
    pub transferred_mass_kg: f64,
    /// Complete stream energy, already including the phase energy difference.
    pub energy_into_gas_j: f64,
    /// Reflected heat INCLUDED in complete stream energy, not external heat.
    pub reflected_heat_into_gas_j: f64,
    /// All solid storage changes, including births and the fixed datum offsets.
    pub solid_energy_change_j: f64,
    pub gas_energy_change_j: f64,
    /// Original birth energy plus selected cells' fixed datum conversions.
    pub incoming_energy_j: f64,
    /// Actual final solid/gas mass change minus incoming births.
    pub mass_residual_kg: f64,
    /// Actual joint storage change minus births and external thermal heat.
    pub energy_residual_j: f64,
    pub phase_steps: usize,
    pub gas_iterations: usize,
    pub last_step: Option<GasChamberStep>,
}

/// The original thermal and gas inputs stay untouched until the caller accepts
/// these endpoints together. The complete catalog retains depleted species.
#[derive(Clone)]
pub struct GasTimelineCandidate<'thermal, 'catalog> {
    thermal: DepositionSimulation<'thermal>,
    gas: ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    audit: GasTimelineAudit,
}
impl<'thermal, 'catalog> GasTimelineCandidate<'thermal, 'catalog> {
    #[must_use]
    pub const fn thermal(&self) -> &DepositionSimulation<'thermal> {
        &self.thermal
    }
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    #[must_use]
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    #[must_use]
    pub const fn audit(&self) -> &GasTimelineAudit {
        &self.audit
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'thermal>,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasTimelineAudit,
    ) {
        (self.thermal, self.gas, self.catalog, self.audit)
    }
}

fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}
fn add(a: f64, b: f64) -> Result<f64, DepositionError> {
    finite(a + b)
}

#[derive(Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, x: f64) -> Result<(), DepositionError> {
        finite(x)?;
        let next = finite(self.sum + x)?;
        let correction = if self.sum.abs() >= x.abs() {
            (self.sum - next) + x
        } else {
            (x - next) + self.sum
        };
        self.correction = finite(self.correction + correction)?;
        self.sum = next;
        Ok(())
    }
    fn value(&self) -> Result<f64, DepositionError> {
        finite(self.sum + self.correction)
    }
}
#[derive(Default)]
struct HeatTotals {
    mass: Sum,
    incoming: Sum,
    external: Sum,
    stored: Sum,
    residual: Sum,
    substeps: usize,
    events: usize,
}
impl HeatTotals {
    fn add(&mut self, part: EnergyAudit) -> Result<(), DepositionError> {
        self.mass.add(part.incoming_mass_kg)?;
        self.incoming.add(part.incoming_energy_j)?;
        self.external.add(part.external_heat_j)?;
        self.stored.add(part.stored_energy_change_j)?;
        self.residual.add(part.residual_j)?;
        self.substeps = self
            .substeps
            .checked_add(part.substeps)
            .ok_or(DepositionError::Arithmetic)?;
        self.events = self
            .events
            .checked_add(part.events)
            .ok_or(DepositionError::Arithmetic)?;
        Ok(())
    }
    fn audit(&self) -> Result<EnergyAudit, DepositionError> {
        Ok(EnergyAudit {
            incoming_mass_kg: self.mass.value()?,
            incoming_energy_j: self.incoming.value()?,
            external_heat_j: self.external.value()?,
            stored_energy_change_j: self.stored.value()?,
            residual_j: self.residual.value()?,
            substeps: self.substeps,
            events: self.events,
        })
    }
}

#[derive(Clone)]
struct SolidBranchGuard {
    ranges: BTreeMap<CellId, [f64; 2]>,
}
impl SolidBranchGuard {
    fn check(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
        if let Some(&[lo, hi]) = self.ranges.get(&cell.id) {
            if !(lo..=hi).contains(&cell.specific_h_j_kg) {
                return Err(DepositionError::MaterialDomain);
            }
        }
        Ok(())
    }
    fn check_all(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        for (index, id) in self.ranges.keys().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if let Some(cell) = thermal.cells().get(id) {
                self.check(cell)?;
            }
        }
        Ok(())
    }
}
impl DepositionObserver for SolidBranchGuard {
    type Error = DepositionError;
    fn event(&mut self, event: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        if let Event::Birth(cell) = event {
            self.check(cell)?;
        }
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.check_all(thermal, checkpoint)
    }
}

/// Advance original heat, scheduled geometry exposure and ONE evolving gas.
///
/// A birth at an interval endpoint first exchanges molecules in the following
/// interval. An endpoint contact change cannot rewrite the preceding exposure.
/// The original thermal substep budget is shared across the COMPLETE request.
/// No phase-specific timestep controller is added; reduce max_dt_s when the
/// explicit molecular owner refuses a gross inventory or caloric endpoint.
///
/// Whole-window mass tolerance is the gas owner's absolute mass tolerance.
/// Whole-window energy tolerance is the larger of the two owners' absolute
/// energy tolerances; neither tolerance is renewed or multiplied by step count.
/// A late failure or cancellation discards all heat, phase, gas and event work.
#[allow(clippy::too_many_arguments)]
pub fn advance_deposition_with_gas<'thermal, 'catalog>(
    thermal: &DepositionSimulation<'thermal>,
    surfaces: &[GasTimelineSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasChamberSolve,
    checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<GasTimelineCandidate<'thermal, 'catalog>, GasTimelineError> {
    advance_deposition(thermal, surfaces, None, gas, catalog, target_s,
        thermal_config, phase_config, checkpoint)
}

/// Original heat timeline with jointly evolving phase exchange and reflected heat.
///
/// The accommodation map must cover exactly ALL selected CellIds, including
/// cells scheduled beyond the requested endpoint. Inactive cells contribute no
/// flux. Each coefficient set is cloned once when its cell first becomes active;
/// the active map is then reused, with no material/source-card cloning here.
/// Contact coverage at the START of each interval applies to both mechanisms.
///
/// Reflection is an internal solid/gas energy transfer, not extra boundary heat.
/// It reaches each original cell's mass/enthalpy endpoint together with phase
/// exchange, before either owner is published. Per-cell and total reflected
/// heat fields are subaccounts of complete stream energy, not additional heat.
///
/// This uses the existing event-clipped, first-order heat/phase splitting and
/// fixed maximum physical timestep. Refine max_dt_s for both coupling and
/// molecular accuracy; the standalone adaptive chamber does not estimate this
/// timeline's splitting error. All original branch, work and rollback gates hold.
#[allow(clippy::too_many_arguments)]
pub fn advance_deposition_with_accommodation<'thermal, 'catalog>(
    thermal: &DepositionSimulation<'thermal>,
    surfaces: &[GasTimelineSurface<'_>],
    accommodation: &BTreeMap<CellId, ThermalAccommodation>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasChamberSolve,
    checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<GasTimelineCandidate<'thermal, 'catalog>, GasTimelineError> {
    advance_deposition(thermal, surfaces, Some(accommodation), gas, catalog, target_s,
        thermal_config, phase_config, checkpoint)
}

#[allow(clippy::too_many_arguments)]
fn advance_deposition<'thermal, 'catalog>(
    thermal: &DepositionSimulation<'thermal>,
    surfaces: &[GasTimelineSurface<'_>],
    accommodation: Option<&BTreeMap<CellId, ThermalAccommodation>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasChamberSolve,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<GasTimelineCandidate<'thermal, 'catalog>, GasTimelineError> {
    checkpoint()?;
    if !target_s.is_finite()
        || target_s < thermal.time_s()
        || !thermal_config.max_dt_s.is_finite()
        || thermal_config.max_dt_s <= 0.0
        || !thermal_config.energy_tolerance_j.is_finite()
        || thermal_config.energy_tolerance_j < 0.0
        || surfaces.is_empty()
    {
        return Err(DepositionError::InvalidInput(
            "finite forward target, positive timestep and selected solid cells required",
        )
        .into());
    }
    let budget = thermal.work_budget();
    if surfaces.len() > budget.max_cells || surfaces.len() > phase_config.maximum_surfaces {
        return Err(DepositionError::Budget("evolving-gas surface selection").into());
    }
    let mut known = BTreeSet::new();
    for (index, id) in thermal.cells().keys().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        known.insert(*id);
    }
    let mut births = BTreeMap::new();
    for (index, event) in thermal.pending_events().iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        if let Event::Birth(cell) = &event.event {
            known.insert(cell.id);
            if event.time_s <= target_s {
                births.entry(cell.id).or_insert(cell);
            }
        }
    }
    let mut ordered: Vec<_> = surfaces.iter().collect();
    ordered.sort_unstable_by_key(|surface| surface.cell_id);
    let mut ranges = BTreeMap::new();
    let mut offsets = BTreeMap::new();
    for (index, surface) in ordered.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        let [lo, hi] = surface.pure_solid_enthalpy_range_j_kg;
        if !known.contains(&surface.cell_id)
            || !lo.is_finite()
            || !hi.is_finite()
            || lo > hi
            || !surface.binding.solid_enthalpy_offset_j_kg.is_finite()
            || ranges.insert(surface.cell_id, [lo, hi]).is_some()
        {
            return Err(DepositionError::InvalidInput(
                "unknown/duplicate cell or invalid pure-solid branch/datum",
            )
            .into());
        }
        offsets.insert(surface.cell_id, surface.binding.solid_enthalpy_offset_j_kg);
    }
    if let Some(laws) = accommodation {
        if laws.len() != ordered.len() {
            return Err(DepositionError::InvalidInput("accommodation cell roster mismatch").into());
        }
        for surface in &ordered {
            checkpoint()?;
            if !laws.contains_key(&surface.cell_id) {
                return Err(DepositionError::InvalidInput("missing selected cell accommodation").into());
            }
        }
    }
    let guard = SolidBranchGuard { ranges };
    guard.check_all(thermal, &mut checkpoint)?;
    // Validate original gas/catalog/controls even for a zero-duration window or
    // when all selected cells are scheduled strictly after the target.
    let mut checkpoint_error = None;
    let validation = advance_chamber_with_gas(&[], gas, catalog, 0.0, phase_config, || {
        match checkpoint() {
            Ok(()) => false,
            Err(error) => {
                checkpoint_error = Some(error);
                true
            }
        }
    });
    if let Some(error) = checkpoint_error {
        return Err(error.into());
    }
    validation?;

    let (mut next, _, initial_heat) = thermal.advance_observed_limited(
        &guard,
        thermal.time_s(),
        thermal_config,
        budget.max_substeps,
        &mut checkpoint,
    )?;
    let mut heat = HeatTotals::default();
    heat.add(initial_heat)?;
    let mut next_gas = gas.clone();
    let mut accounts: BTreeMap<_, _> = ordered
        .iter()
        .map(|surface| (surface.cell_id, GasCellAudit::default()))
        .collect();
    // Cell extinction/removal is refused by the original whole-cell model.
    // Each selected cell therefore enters this map at most once per advance.
    let mut active_accommodation = BTreeMap::new();
    let mut outgoing = Sum::default();
    let mut incoming = Sum::default();
    let mut transferred = Sum::default();
    let mut stream_energy = Sum::default();
    let mut reflected_energy = Sum::default();
    let mut phase_steps = 0usize;
    let mut gas_iterations = 0usize;
    let mut last_step = None;
    while next.time_s() < target_s {
        checkpoint()?;
        let start_s = next.time_s();
        let mut exposed = Vec::with_capacity(ordered.len());
        for (index, surface) in ordered.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if next.cells().contains_key(&surface.cell_id) {
                exposed.push((*surface, next.exposed_area_m2(surface.cell_id)?));
                if let Some(laws) = accommodation {
                    active_accommodation.entry(surface.cell_id.0)
                        .or_insert_with(|| laws[&surface.cell_id].clone());
                }
            }
        }
        // Empty pre-birth time carries no thermal or phase work. Jump to the
        // next event/target instead of spending unbudgeted max-dt iterations.
        let mut stop_s = if next.cells().is_empty() {
            target_s
        } else {
            (start_s + thermal_config.max_dt_s).min(target_s)
        };
        if let Some(event) = next.pending_events().first() {
            stop_s = stop_s.min(event.time_s);
        }
        if !stop_s.is_finite() || stop_s <= start_s {
            return Err(DepositionError::Arithmetic.into());
        }
        let remaining = budget
            .max_substeps
            .checked_sub(heat.substeps)
            .ok_or(DepositionError::Budget("evolving-gas thermal substeps"))?;
        let (heated, _, report) = next.advance_observed_limited(
            &guard,
            stop_s,
            thermal_config,
            remaining,
            &mut checkpoint,
        )?;
        next = heated;
        heat.add(report)?;
        if exposed.is_empty() {
            continue;
        }

        let chart_error = RefCell::new(None);
        let calorics: Vec<_> = exposed
            .iter()
            .map(|(surface, _)| {
                let state = &next;
                let errors = &chart_error;
                let id = surface.cell_id;
                let [lo, hi] = surface.pure_solid_enthalpy_range_j_kg;
                move |h| {
                    let result = if (lo..=hi).contains(&h) {
                        state.temperature_at_specific_enthalpy_k(id, h)
                    } else {
                        Err(DepositionError::MaterialDomain)
                    };
                    result.map_err(|error| {
                        *errors.borrow_mut() = Some(error);
                        SublimationError::Invalid(
                            "original pure-solid caloric chart refused endpoint",
                        )
                    })
                }
            })
            .collect();
        let inputs: Vec<_> = exposed
            .iter()
            .zip(&calorics)
            .map(|((surface, area), caloric)| {
                let cell = &next.cells()[&surface.cell_id];
                GasChamberSurface {
                    id: surface.cell_id.0,
                    law: surface.law,
                    binding: surface.binding,
                    solid: SolidSurfaceState {
                        mass_kg: cell.mass_kg,
                        specific_enthalpy_j_kg: cell.specific_h_j_kg,
                        area_m2: *area,
                    },
                    temperature_from_enthalpy: caloric,
                }
            })
            .collect();
        let mut checkpoint_error = None;
        let result = {
            let mut cancelled = || match checkpoint() {
                Ok(()) => false,
                Err(error) => {
                    checkpoint_error = Some(error);
                    true
                }
            };
            if accommodation.is_some() {
                advance_chamber_with_accommodation(&inputs, &active_accommodation,
                    &next_gas, catalog, stop_s - start_s, phase_config, &mut cancelled)
            } else {
                advance_chamber_with_gas(&inputs, &next_gas, catalog,
                    stop_s - start_s, phase_config, &mut cancelled)
                    .map(|candidate| (candidate, BTreeMap::new()))
            }
        };
        drop(inputs);
        drop(calorics);
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        let (candidate, reflection) = match result {
            Ok(value) => value,
            Err(error) => {
                return Err(match chart_error.into_inner() {
                    Some(error) => error.into(),
                    None => error.into(),
                });
            }
        };
        if candidate.surfaces().len() != exposed.len() {
            return Err(DepositionError::InvalidInput("joint phase endpoint count changed").into());
        }
        let mut updates = Vec::with_capacity(exposed.len());
        for (index, ((surface, area), endpoint)) in
            exposed.iter().zip(candidate.surfaces()).enumerate()
        {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if endpoint.id != surface.cell_id.0 {
                return Err(
                    DepositionError::InvalidInput("joint phase endpoint identity changed").into(),
                );
            }
            updates.push(SolidCellExchange {
                expected: next.cells()[&surface.cell_id].clone(),
                interval_exposed_area_m2: *area,
                replacement: endpoint.solid,
                pure_solid_enthalpy_range_j_kg: surface.pure_solid_enthalpy_range_j_kg,
            });
        }
        let iterations = candidate.gas_exchange().step().iterations;
        next.apply_solid_exchanges(&updates, &mut checkpoint)?;
        let (endpoints, accepted_gas, _, report) = candidate.into_parts();
        next_gas = accepted_gas;
        outgoing.add(report.outgoing_mass_kg)?;
        incoming.add(report.incoming_mass_kg)?;
        transferred.add(report.transferred_mass_kg)?;
        stream_energy.add(report.energy_into_gas_j)?;
        phase_steps = phase_steps
            .checked_add(1)
            .ok_or(DepositionError::Arithmetic)?;
        gas_iterations = gas_iterations
            .checked_add(iterations)
            .ok_or(DepositionError::Arithmetic)?;
        for endpoint in endpoints {
            let account = accounts
                .get_mut(&CellId(endpoint.id))
                .expect("selected gas surface");
            account.outgoing_mass_kg = add(account.outgoing_mass_kg, endpoint.outgoing_mass_kg)?;
            account.incoming_mass_kg = add(account.incoming_mass_kg, endpoint.incoming_mass_kg)?;
            account.transferred_mass_kg =
                add(account.transferred_mass_kg, endpoint.transferred_mass_kg)?;
            account.energy_into_gas_j = add(account.energy_into_gas_j, endpoint.energy_into_gas_j)?;
            if let Some(reflected) = reflection.get(&endpoint.id) {
                account.reflected_heat_into_gas_j =
                    add(account.reflected_heat_into_gas_j, reflected.heat_into_gas_j)?;
                reflected_energy.add(reflected.heat_into_gas_j)?;
            }
            account.recession_m = add(account.recession_m, endpoint.recession_m)?;
            account.phase_steps = account
                .phase_steps
                .checked_add(1)
                .ok_or(DepositionError::Arithmetic)?;
        }
        last_step = Some(report);
    }
    let thermal_audit = heat.audit()?;
    let mut solid_mass_without_births = Sum::default();
    let mut solid_energy_without_births = Sum::default();
    let mut birth_offsets = Sum::default();
    // Compare each actual endpoint to its original state or incoming birth.
    // This retains small heat against large formation-enthalpy references.
    for (index, (&id, cell)) in next.cells().iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        let offset = offsets.get(&id).copied().unwrap_or(0.0);
        let original = if let Some(original) = thermal.cells().get(&id) {
            original
        } else {
            let original = births.get(&id).ok_or(DepositionError::InvalidInput(
                "thermal endpoint has no original cell or incoming birth",
            ))?;
            birth_offsets.add(original.mass_kg * offset)?;
            original
        };
        let dm = finite(cell.mass_kg - original.mass_kg)?;
        solid_mass_without_births.add(dm)?;
        solid_energy_without_births.add(
            cell.mass_kg * (cell.specific_h_j_kg - original.specific_h_j_kg)
                + dm * (original.specific_h_j_kg + offset),
        )?;
    }
    let gas_mass_change = finite(next_gas.mass_kg() - gas.mass_kg())?;
    let initial_u = gas.state().specific_internal_energy_j_kg();
    let gas_energy_change = finite(
        next_gas.mass_kg() * (next_gas.state().specific_internal_energy_j_kg() - initial_u)
            + gas_mass_change * initial_u,
    )?;
    let mass_residual = add(solid_mass_without_births.value()?, gas_mass_change)?;
    let mass_tolerance = phase_config.surface.gas.mass_tolerance_kg;
    if mass_residual.abs() > mass_tolerance {
        return Err(GasTimelineError::MassBalance {
            residual_kg: mass_residual,
            tolerance_kg: mass_tolerance,
        });
    }
    let energy_residual = finite(
        solid_energy_without_births.value()? + gas_energy_change - thermal_audit.external_heat_j,
    )?;
    let energy_tolerance = thermal_config
        .energy_tolerance_j
        .max(phase_config.surface.gas.energy_tolerance_j);
    if energy_residual.abs() > energy_tolerance {
        return Err(GasSurfaceError::EnergyBalance {
            residual_j: energy_residual,
            tolerance_j: energy_tolerance,
        }
        .into());
    }
    let incoming_energy = add(thermal_audit.incoming_energy_j, birth_offsets.value()?)?;
    let audit = GasTimelineAudit {
        thermal: thermal_audit,
        cells: accounts,
        outgoing_mass_kg: outgoing.value()?,
        incoming_mass_kg: incoming.value()?,
        transferred_mass_kg: transferred.value()?,
        energy_into_gas_j: stream_energy.value()?,
        reflected_heat_into_gas_j: reflected_energy.value()?,
        solid_energy_change_j: add(solid_energy_without_births.value()?, incoming_energy)?,
        gas_energy_change_j: gas_energy_change,
        incoming_energy_j: incoming_energy,
        mass_residual_kg: mass_residual,
        energy_residual_j: energy_residual,
        phase_steps,
        gas_iterations,
        last_step,
    };
    checkpoint()?;
    Ok(GasTimelineCandidate {
        thermal: next,
        gas: next_gas,
        catalog,
        audit,
    })
}

#[cfg(test)]
#[path = "timeline/tests.rs"]
mod tests;

#[cfg(test)]
#[path = "timeline/accommodation_tests.rs"]
mod accommodation_tests;
