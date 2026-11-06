//! Native thermal events and finite-gas film exchange on one physical clock.
//!
//! Only the original implicit endpoint advances heat and phase. Native event
//! application is explicitly limited to zero thermal substeps. Film inventory
//! survives boundary/contact changes on its fixed exposed patch; a newly
//! covering contact refuses rather than silently removing liquid or wall area.

use super::super::{
    implicit::WorkMeter,
    window::{Sum, WallFilmWindowSolve, increment},
};
use super::{
    DepositionWallFilmBinding, DepositionWallFilmControls, DepositionWallFilmReport,
    advance_deposition_with_wall_film_metered, work_error,
};
use fs_material::{
    evaporation::gas::{
        GasEvaporationError, implicit::ImplicitGasError, wall_film::GasWallFilmInventory,
    },
    thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel},
};
use fs_thermal::enthalpy::{
    DepositionError, DepositionSimulation, StepConfig,
    backward_euler::BoundaryExchangeStepCoefficients,
};

mod account;
use account::Account;

/// Fixed nominal intervals, native event cuts and one complete-request budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DepositionWallFilmTimelineControls {
    /// Numerical iteration/evaluation limits span the complete timeline call.
    pub step: DepositionWallFilmControls,
    /// Whole-call physical balances, attempts and nominal interval bounds.
    /// Event cuts and the final remainder may be below the retry minimum.
    pub window: WallFilmWindowSolve,
    /// Includes original gas inversions on every attempted endpoint.
    pub maximum_gas_iterations: usize,
}
impl DepositionWallFilmTimelineControls {
    pub fn new(step: DepositionWallFilmControls, window: WallFilmWindowSolve) -> Self {
        Self {
            step,
            window,
            maximum_gas_iterations: 100_000,
        }
    }
}

/// Accepted physical accounts; work includes refused proposals as well.
#[derive(Clone, Debug, PartialEq)]
pub struct DepositionWallFilmTimelineReport {
    pub initial_time_s: f64,
    pub final_time_s: f64,
    pub accepted_intervals: usize,
    pub attempted_intervals: usize,
    pub rejected_intervals: usize,
    pub events: usize,
    pub incoming_mass_kg: f64,
    /// Native birth storage in each ORIGINAL wall material's declared datum.
    pub incoming_energy_j: f64,
    /// Native sources and prescribed baths, excluding internal gas/film heat.
    pub external_thermal_heat_j: f64,
    pub external_gas_heat_j: f64,
    /// Original-cell and birth-relative storage increments, excluding birth h.
    pub wall_storage_change_j: f64,
    pub film_enthalpy_change_j: f64,
    pub gas_internal_energy_change_j: f64,
    pub evaporated_mass_kg: f64,
    pub condensed_mass_kg: f64,
    /// Diagnostic phase component already included in the full energy account.
    pub latent_energy_j: f64,
    pub vapor_energy_into_gas_j: f64,
    pub sensible_heat_into_gas_j: f64,
    pub wall_heat_into_phase_system_j: f64,
    pub gas_heat_apportionment_residual_j: f64,
    pub thermal_telescoping_residual_j: f64,
    pub combined_energy_residual_j: f64,
    pub sum_absolute_energy_residual_j: f64,
    pub mass_residual_kg: f64,
    pub volatile_mass_residual_kg: f64,
    pub maximum_carrier_mass_residual_kg: f64,
    pub sum_absolute_mass_residual_kg: f64,
    pub newton_iterations: usize,
    pub residual_evaluations: usize,
    pub gas_iterations: usize,
    pub charged_gas_iterations: usize,
    /// Last physical endpoint, before any native events at the final time.
    pub last_step: Option<DepositionWallFilmReport>,
}

/// Native geometry/history, gas and the separate film are adopted together.
pub struct DepositionWallFilmTimelineCandidate<'m, 'catalog, T> {
    thermal: DepositionSimulation<'m>,
    film: GasWallFilmInventory,
    gas: ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    payload: Option<T>,
    report: DepositionWallFilmTimelineReport,
}
impl<'m, 'catalog, T> DepositionWallFilmTimelineCandidate<'m, 'catalog, T> {
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn film(&self) -> GasWallFilmInventory {
        self.film
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub fn payload(&self) -> Option<&T> {
        self.payload.as_ref()
    }
    pub const fn report(&self) -> &DepositionWallFilmTimelineReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'m>,
        GasWallFilmInventory,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        Option<T>,
        DepositionWallFilmTimelineReport,
    ) {
        (
            self.thermal,
            self.film,
            self.gas,
            self.catalog,
            self.payload,
            self.report,
        )
    }
}

/// Advance the same initially active wall patch through native authored events.
///
/// Native births, contacts, separations and boundary changes apply in their
/// original equal-time order, including events at the initial and final times.
/// The heat callback sees the left-limit state at an event endpoint. Its last
/// retained payload precedes final-time events; it is not recomputed as heat.
///
/// `external_gas_heat_j` is TOTAL uniform-source heat over the requested window.
/// Native wall sources keep their original watts and event schedule. Only an
/// original phase StepSizeLimit permits a shorter retry. Every proposal spends
/// shared work, while only accepted intervals contribute heat, mass or time.
///
/// The physical patch area is fixed, and must fit the actual exposed wall after
/// every event. New contacts do not displace, drain or hide film inventory.
/// The original wall cell must already exist. Unselected native cells may be
/// born later, carrying their own original mass and caloric datum. This remains
/// first-order fixed-grid integration; refine time and spatial cells separately.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_deposition_with_wall_film_timeline<'m, 'catalog, T>(
    thermal: &DepositionSimulation<'m>,
    binding: DepositionWallFilmBinding<'_>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    external_gas_heat_j: f64,
    controls: DepositionWallFilmTimelineControls,
    mut resolve_heat: impl FnMut(
        &DepositionSimulation<'m>,
        &ClosedGasParcel,
    ) -> Result<BoundaryExchangeStepCoefficients<T>, ImplicitGasError>,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<DepositionWallFilmTimelineCandidate<'m, 'catalog, T>, ImplicitGasError> {
    checkpoint()?;
    let origin = thermal.time_s();
    let duration = finite(target_s - origin)?;
    let window = controls.window;
    if !target_s.is_finite()
        || duration <= 0.0
        || !external_gas_heat_j.is_finite()
        || !positive(window.maximum_step_s)
        || !positive(window.minimum_step_s)
        || window.minimum_step_s > window.maximum_step_s
        || window.maximum_attempts == 0
        || !positive(window.energy_tolerance_j)
        || !positive(window.mass_tolerance_kg)
    {
        return Err(ImplicitGasError::Invalid(
            "invalid wall-film thermal timeline controls",
        ));
    }
    admit_patch(thermal, binding)?;
    let meter = WorkMeter::new(controls.step.numerics(), controls.maximum_gas_iterations)
        .map_err(work_error)?;
    let mut next = thermal.clone();
    let mut next_gas = gas.clone();
    let mut film = binding.film;
    let mut account = Account::default();
    let mut accepted = 0usize;
    let mut attempted = 0usize;
    let mut rejected = 0usize;
    let mut payload = None;
    let mut last_step = None;
    loop {
        checkpoint()?;
        let now = next.time_s();
        if next
            .pending_events()
            .first()
            .is_some_and(|event| event.time_s == now)
        {
            // The native event owner cannot insert an explicit heat step here.
            let (advanced, (), audit) = next.advance_observed_limited(
                &(),
                now,
                StepConfig {
                    max_dt_s: window.maximum_step_s,
                    energy_tolerance_j: window.energy_tolerance_j,
                },
                0,
                &mut checkpoint,
            )?;
            if advanced.time_s() != now || audit.substeps != 0 {
                return Err(ImplicitGasError::Invalid(
                    "native event pass advanced thermal time",
                ));
            }
            account.add_events(next.pending_events(), &audit, &mut checkpoint)?;
            admit_patch(&advanced, DepositionWallFilmBinding { film, ..binding })?;
            account.gate(window)?;
            next = advanced;
        }
        if next.time_s() == target_s {
            break;
        }
        if accepted >= thermal.work_budget().max_substeps {
            return Err(DepositionError::Budget("thermal substeps").into());
        }
        let start = next.time_s();
        let mut end = if window.maximum_step_s >= target_s - start {
            target_s
        } else {
            let mut nominal = finite(start + window.maximum_step_s)?;
            if nominal - start > window.maximum_step_s {
                nominal = nominal.next_down();
            }
            nominal
        };
        if let Some(event) = next.pending_events().first() {
            end = end.min(event.time_s);
        }
        if !(start < end && end <= target_s) {
            return Err(ImplicitGasError::Invalid(
                "wall-film timeline cannot represent a forward interval",
            ));
        }
        loop {
            checkpoint()?;
            if attempted >= window.maximum_attempts {
                return Err(DepositionError::Budget("wall-film timeline attempts").into());
            }
            let part_heat = increment(external_gas_heat_j, start - origin, end - origin, duration)
                .map_err(|_| DepositionError::Arithmetic)?;
            attempted = count(attempted, 1)?;
            let candidate = advance_deposition_with_wall_film_metered(
                &next,
                DepositionWallFilmBinding { film, ..binding },
                &next_gas,
                catalog,
                end,
                part_heat,
                controls.step,
                &meter,
                &mut resolve_heat,
                &mut checkpoint,
            );
            match candidate {
                Ok(candidate) => {
                    let (advanced, next_film, advanced_gas, _, source, step) =
                        candidate.into_parts();
                    if advanced.time_s() != end || step.phase.duration_s != end - start {
                        return Err(ImplicitGasError::Invalid(
                            "wall-film endpoint changed its shared clock",
                        ));
                    }
                    account.add(&step)?;
                    account.gate(window)?;
                    accepted = count(accepted, 1)?;
                    next = advanced;
                    film = next_film;
                    next_gas = advanced_gas;
                    payload = Some(source);
                    last_step = Some(step);
                    break;
                }
                Err(
                    error @ ImplicitGasError::Phase(GasEvaporationError::StepSizeLimit {
                        maximum_duration_s,
                        ..
                    }),
                ) => {
                    rejected = count(rejected, 1)?;
                    let dt = end - start;
                    let wanted = maximum_duration_s.min(0.5 * dt);
                    if !positive(maximum_duration_s)
                        || !positive(wanted)
                        || wanted < window.minimum_step_s
                    {
                        return Err(error);
                    }
                    let mut shorter = start + wanted;
                    if shorter - start > wanted {
                        shorter = shorter.next_down();
                    }
                    if !(start < shorter && shorter < end)
                        || shorter - start < window.minimum_step_s
                        || shorter - start > wanted
                    {
                        return Err(error);
                    }
                    end = shorter;
                }
                Err(error) => return Err(error),
            }
        }
    }
    let report = account.finish(
        thermal,
        &next,
        binding,
        gas,
        &next_gas,
        film,
        external_gas_heat_j,
        window,
        accepted,
        attempted,
        rejected,
        &meter,
        last_step,
        &mut checkpoint,
    )?;
    checkpoint()?;
    Ok(DepositionWallFilmTimelineCandidate {
        thermal: next,
        film,
        gas: next_gas,
        catalog,
        payload,
        report,
    })
}

fn admit_patch(
    thermal: &DepositionSimulation<'_>,
    binding: DepositionWallFilmBinding<'_>,
) -> Result<(), ImplicitGasError> {
    if !binding.boundary.area_m2.is_finite()
        || binding.boundary.area_m2 < 0.0
        || binding.boundary.area_m2 > thermal.exposed_area_m2(binding.cell_id)?
    {
        return Err(ImplicitGasError::Invalid(
            "wall-film patch must fit its active exposed native wall",
        ));
    }
    Ok(())
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn finite(value: f64) -> Result<f64, ImplicitGasError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}
fn count(value: usize, increment: usize) -> Result<usize, ImplicitGasError> {
    value
        .checked_add(increment)
        .ok_or_else(|| DepositionError::Arithmetic.into())
}
fn add(sum: &mut Sum, value: f64) -> Result<(), ImplicitGasError> {
    sum.add(value)
        .map_err(|_| DepositionError::Arithmetic.into())
}
fn difference(mut sum: Sum, reference: f64) -> Result<f64, ImplicitGasError> {
    add(&mut sum, -reference)?;
    finite(sum.value())
}
