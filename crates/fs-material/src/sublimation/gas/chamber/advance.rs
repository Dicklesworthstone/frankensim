//! Bounded adaptive histories for a closed nonisothermal sublimation chamber.
//!
//! Each trial compares one original explicit interval with two half intervals.
//! The accepted state is the two-half-step conservative endpoint; there is no
//! extrapolation, clipping, or replacement caloric law. Temperature and EVERY
//! species inventory participate in the local error estimate, including trace
//! vapor in a carrier. Fixed canonical surface order makes retry decisions
//! independent of the input order on the same ISA.
//!
//! Optional per-surface thermal accommodation uses the SAME chamber kernel,
//! error controller and accepted-state transaction. Reflected heat participates
//! in both temperatures even when mass transfer is exactly zero. Only accepted
//! half steps contribute to its history; coarse and rejected trials contribute
//! neither state nor heat. Coefficient maps are borrowed, not cloned per trial.
//!
//! The moving-front entry point refreshes mass-owned geometry on every trial
//! and half step and also controls radius/thickness error. Existing entry points
//! retain fixed areas and their original physical update. The partitioned entry
//! point additionally controls BOTH condensed phase masses, resolving latent
//! plateaus even when total mass and temperature do not change. Its source
//! charts are borrowed across trials, not rebuilt or replaced by apparent Cp.
//!
//! The estimate controls local temporal error, not a rigorous global bound or
//! experimental accuracy. It supplies no external heat, topology change,
//! gas-flow model or thermal event scheduling. Mixed inventory and moving
//! pure-body geometry are separate model choices, not automatically composable.
//! Inventory-limit refusals alone are retryable. Invalid constitutive data,
//! unsupported caloric endpoints and cancellation remain immediate failures.

/// Mass-owned changing area on this same adaptive, conservative history owner.
pub mod moving;
/// Actual solid/liquid inventory error control on this same adaptive owner.
pub mod partitioned;

use super::super::{ClosedGasParcel, GasSurfaceError, SublimationError, ThermochemicalGasModel};
use super::super::accommodation::{ReflectedHeatStep, ThermalAccommodation};
use super::mixed::{self as phase_inventory, SharpMeltingCaloric};
use super::{
    GasChamberCandidate, GasChamberSolve, GasChamberSurface, GasChamberSurfaceStep,
    advance_chamber_with_accommodation, advance_chamber_with_gas,
};
use fs_thermal::sublimation::SolidSurfaceState;
use std::collections::BTreeMap;

/// Local physical-state tolerances and a strict work ceiling for the whole call.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasChamberAdvanceSolve {
    pub step: GasChamberSolve,
    /// Absolute error scale for each condensed body and EACH gas species [kg].
    /// The partitioned lane also applies this to each SOLID and LIQUID mass.
    pub absolute_mass_tolerance_kg: f64,
    /// Absolute error scale for each surface and the gas [K].
    pub absolute_temperature_tolerance_k: f64,
    /// Added to each absolute scale as relative_tolerance * max(|before|,|after|).
    pub relative_tolerance: f64,
    pub minimum_dt_s: f64,
    pub maximum_dt_s: f64,
    /// Includes rejected trials; each trial uses at most three kernel calls.
    pub maximum_trials: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum GasChamberAdvanceError {
    Invalid(&'static str),
    Surface(GasSurfaceError),
    TrialBudget { maximum_trials: usize },
    MinimumStep { requested_dt_s: f64, minimum_dt_s: f64 },
}
impl core::fmt::Display for GasChamberAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive molecular chamber: {self:?}")
    }
}
impl std::error::Error for GasChamberAdvanceError {}
impl From<GasSurfaceError> for GasChamberAdvanceError {
    fn from(value: GasSurfaceError) -> Self { Self::Surface(value) }
}

/// Integrated streams along the accepted history, in canonical surface order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasSurfaceHistory {
    pub id: u64,
    pub outgoing_mass_kg: f64,
    pub incoming_mass_kg: f64,
    pub transferred_mass_kg: f64,
    pub energy_into_gas_j: f64,
    /// Signed reflected heat INCLUDED in energy_into_gas_j, never an addend.
    /// This sums only accepted half steps and is zero for the specular API.
    pub reflected_heat_into_gas_j: f64,
    /// Sum of accepted normal retreat. Fixed-area paths use dm/(rho*A);
    /// the moving path uses the geometry's actual shell displacement.
    /// Mixed-body output is vapor-transfer-equivalent thickness, not geometry
    /// change from fusion or a reconstructed wetting/contact boundary.
    pub recession_m: f64,
    pub surface_temperature_k: [f64; 2],
}

#[derive(Clone, Debug, PartialEq)]
pub struct GasChamberAdvanceAudit {
    pub elapsed_s: f64,
    pub trials: usize,
    pub rejected_trials: usize,
    /// Each accepted interval contains two physical half steps.
    pub accepted_intervals: usize,
    pub kernel_evaluations: usize,
    pub maximum_accepted_error_ratio: f64,
    pub surfaces: Vec<GasSurfaceHistory>,
    pub energy_into_gas_j: f64,
    /// Signed reflected heat INCLUDED in the complete stream energy above.
    pub reflected_heat_into_gas_j: f64,
    /// Sum of original-reference storage increments on accepted half steps.
    /// Rounded diagnostic totals need not reproduce the small residuals below.
    pub solid_energy_change_j: f64,
    pub gas_energy_change_j: f64,
    /// Telescoping caloric defects, with conserved species datums canceled
    /// inside each accepted physical interval before accumulation.
    pub energy_residual_j: f64,
    pub gas_stream_energy_residual_j: f64,
    pub solid_stream_energy_residual_j: f64,
    pub mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
}

/// A complete proposed history. Caller-owned states remain unchanged on failure.
#[derive(Clone, Debug, PartialEq)]
pub struct GasChamberAdvanceCandidate<'catalog> {
    endpoint: GasChamberCandidate<'catalog>,
    audit: GasChamberAdvanceAudit,
}
impl<'catalog> GasChamberAdvanceCandidate<'catalog> {
    /// Endpoint physical states; the endpoint's step report is ONLY its last half step.
    #[must_use]
    pub const fn endpoint(&self) -> &GasChamberCandidate<'catalog> { &self.endpoint }
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel { self.endpoint.gas() }
    #[must_use]
    pub const fn audit(&self) -> &GasChamberAdvanceAudit { &self.audit }
    #[must_use]
    pub fn into_parts(self) -> (GasChamberCandidate<'catalog>, GasChamberAdvanceAudit) {
        (self.endpoint, self.audit)
    }
}

#[derive(Clone, Copy, Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), GasChamberAdvanceError> {
        let adjusted = value - self.correction;
        let next = self.value + adjusted;
        self.correction = (next - self.value) - adjusted;
        self.value = next;
        if !self.value.is_finite() || !self.correction.is_finite() {
            return Err(GasChamberAdvanceError::Invalid("history account overflow"));
        }
        Ok(())
    }
}
#[derive(Default)]
struct SurfaceSums {
    outgoing: Sum,
    incoming: Sum,
    transferred: Sum,
    energy: Sum,
    stored: Sum,
    residual: Sum,
    reflected: Sum,
    recession: Sum,
    initial_temperature: Option<f64>,
}
impl SurfaceSums {
    fn add(&mut self, row: &GasChamberSurfaceStep, reflection: Option<&ReflectedHeatStep>)
        -> Result<(), GasChamberAdvanceError> {
        self.initial_temperature.get_or_insert(row.surface_temperature_k[0]);
        self.outgoing.add(row.outgoing_mass_kg)?;
        self.incoming.add(row.incoming_mass_kg)?;
        self.transferred.add(row.transferred_mass_kg)?;
        self.energy.add(row.energy_into_gas_j)?;
        self.stored.add(row.solid_energy_change_j)?;
        self.residual.add(row.energy_residual_j)?;
        if let Some(reflection) = reflection { self.reflected.add(reflection.heat_into_gas_j)?; }
        self.recession.add(row.recession_m)
    }
}

fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GasChamberAdvanceError> {
    if cancelled() { Err(GasSurfaceError::Surface(SublimationError::Cancelled).into()) } else { Ok(()) }
}
fn with_states<'a>(templates: &[GasChamberSurface<'a>], states: &[SolidSurfaceState]) -> Vec<GasChamberSurface<'a>> {
    templates.iter().zip(states).map(|(template, &solid)| GasChamberSurface { solid, ..*template }).collect()
}
fn error_ratio(a: f64, b: f64, initial: f64, absolute: f64, relative: f64) -> Result<f64, GasChamberAdvanceError> {
    let scale = absolute + relative * initial.abs().max(a.abs()).max(b.abs());
    let ratio = (a - b).abs() / scale;
    if !positive(scale) || !ratio.is_finite() {
        Err(GasChamberAdvanceError::Invalid("unrepresentable local error scale"))
    } else { Ok(ratio) }
}
fn compare(
    before: &[GasChamberSurface<'_>], gas: &ClosedGasParcel, catalog: &ThermochemicalGasModel,
    coarse: &GasChamberCandidate<'_>, fine: &GasChamberCandidate<'_>, controls: GasChamberAdvanceSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, GasChamberAdvanceError> {
    let mut maximum = error_ratio(coarse.gas().state().state().temperature,
        fine.gas().state().state().temperature, gas.state().state().temperature,
        controls.absolute_temperature_tolerance_k, controls.relative_tolerance)?;
    for ((original, a), b) in before.iter().zip(coarse.surfaces()).zip(fine.surfaces()) {
        poll(cancelled)?;
        maximum = maximum.max(error_ratio(a.solid.mass_kg, b.solid.mass_kg, original.solid.mass_kg,
            controls.absolute_mass_tolerance_kg, controls.relative_tolerance)?);
        maximum = maximum.max(error_ratio(a.surface_temperature_k[1], b.surface_temperature_k[1],
            a.surface_temperature_k[0], controls.absolute_temperature_tolerance_k, controls.relative_tolerance)?);
    }
    for species in catalog.thermodynamics().components() {
        poll(cancelled)?;
        let id = species.species();
        maximum = maximum.max(error_ratio(coarse.gas().species_mass_kg(id).map_err(GasSurfaceError::from)?,
            fine.gas().species_mass_kg(id).map_err(GasSurfaceError::from)?,
            gas.species_mass_kg(id).map_err(GasSurfaceError::from)?,
            controls.absolute_mass_tolerance_kg, controls.relative_tolerance)?);
    }
    Ok(maximum)
}
fn smaller_step(requested: f64, old: f64, minimum: f64) -> Result<f64, GasChamberAdvanceError> {
    let dt = requested.min(0.5 * old);
    if !positive(dt) || dt < minimum {
        Err(GasChamberAdvanceError::MinimumStep { requested_dt_s: dt, minimum_dt_s: minimum })
    } else { Ok(dt) }
}

// One physical kernel dispatch, shared by trial, half-step and zero-time paths.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn evaluate_interval<'catalog>(
    surfaces: &[GasChamberSurface<'_>], accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    geometry: Option<&moving::MovingChamberGeometry>,
    calorics: Option<&BTreeMap<u64, SharpMeltingCaloric<'_>>>, gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel, duration_s: f64, controls: GasChamberSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(GasChamberCandidate<'catalog>, BTreeMap<u64, ReflectedHeatStep>), GasSurfaceError> {
    if let Some(geometry) = geometry { geometry.admit_interval(surfaces, gas, &mut cancelled)?; }
    let mut result = if let Some(calorics) = calorics {
        // Reuse the joint snapshot/aggregation owner; do not construct a second
        // per-trial phase report or reconstruct gas once per condensed body.
        super::advance_chamber(surfaces, accommodation, Some(calorics), gas, catalog,
            duration_s, controls, &mut cancelled)
    } else {
        match accommodation {
            Some(laws) => advance_chamber_with_accommodation(surfaces, laws, gas, catalog, duration_s, controls, &mut cancelled),
            None => advance_chamber_with_gas(surfaces, gas, catalog, duration_s, controls, &mut cancelled)
                .map(|endpoint| (endpoint, BTreeMap::new())),
        }
    }?;
    if let Some(geometry) = geometry { geometry.refresh_endpoint(surfaces, &mut result.0, &mut cancelled)?; }
    Ok(result)
}

/// Advance all surfaces and the gas through a complete physical duration.
///
/// Temperatures and species masses set the error estimate; formation-energy
/// datums do not directly enter the error norm. The finer conservative
/// history is accepted directly, without Richardson extrapolation. Every
/// rejected trial is discarded. Total actual solid/gas mass and energy must
/// meet controls.step.surface.gas tolerances over the WHOLE duration.
/// A final remainder shorter than minimum_dt_s may be attempted once; any
/// needed further reduction refuses. All models, bindings and exposed areas
/// stay fixed during the call. Use the thermal timeline adapter for events.
pub fn advance_chamber_adaptive<'catalog>(
    surfaces: &[GasChamberSurface<'_>], gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel, duration_s: f64, controls: GasChamberAdvanceSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasChamberAdvanceCandidate<'catalog>, GasChamberAdvanceError> {
    advance_history(surfaces, None, None, gas, catalog, duration_s, controls, cancelled)
}

/// Adaptive phase exchange AND reflected heat using one original chamber owner.
///
/// Coefficients must cover exactly the surface identities, including explicit
/// zero coefficients for specular surfaces. They remain fixed through the call;
/// original gas and surface temperatures are reevaluated on EVERY half step.
/// No coefficient or source catalog is cloned by this driver. Reflected heat
/// enters the existing temperature error estimates, even at zero net mass.
/// The audit's reflection fields are subaccounts of its complete stream energy.
///
/// The same trial/minimum-step limits, gross inventory restrictions, source
/// domains, cancellation and whole-window conservation gates apply. Source
/// refusal is not reclassified as timestep error; an unsupported coarse trial
/// still refuses. Choose max_dt_s within the admitted physical regime.
#[allow(clippy::too_many_arguments)]
pub fn advance_chamber_adaptive_with_accommodation<'catalog>(
    surfaces: &[GasChamberSurface<'_>], accommodation: &BTreeMap<u64, ThermalAccommodation>,
    gas: &ClosedGasParcel, catalog: &'catalog ThermochemicalGasModel, duration_s: f64,
    controls: GasChamberAdvanceSolve, cancelled: impl FnMut() -> bool,
) -> Result<GasChamberAdvanceCandidate<'catalog>, GasChamberAdvanceError> {
    advance_history(surfaces, Some(accommodation), None, gas, catalog, duration_s, controls, cancelled)
}

// Preserve the original entry used by pure-solid and moving-body callers.
#[allow(clippy::too_many_arguments)]
fn advance_history<'catalog>(
    surfaces: &[GasChamberSurface<'_>], accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    geometry: Option<&moving::MovingChamberGeometry>, gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel, duration_s: f64, controls: GasChamberAdvanceSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasChamberAdvanceCandidate<'catalog>, GasChamberAdvanceError> {
    advance_history_with_calorics(surfaces, accommodation, geometry, None,
        gas, catalog, duration_s, controls, cancelled)
}

#[allow(clippy::too_many_arguments)]
fn advance_history_with_calorics<'catalog>(
    surfaces: &[GasChamberSurface<'_>], accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    geometry: Option<&moving::MovingChamberGeometry>,
    calorics: Option<&BTreeMap<u64, SharpMeltingCaloric<'_>>>, gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel, duration_s: f64, controls: GasChamberAdvanceSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasChamberAdvanceCandidate<'catalog>, GasChamberAdvanceError> {
    poll(&mut cancelled)?;
    if !duration_s.is_finite() || duration_s < 0.0
        || !positive(controls.absolute_mass_tolerance_kg) || !positive(controls.absolute_temperature_tolerance_k)
        || !controls.relative_tolerance.is_finite() || controls.relative_tolerance < 0.0
        || !positive(controls.minimum_dt_s) || !positive(controls.maximum_dt_s)
        || controls.minimum_dt_s > controls.maximum_dt_s || controls.maximum_trials == 0
        || surfaces.len() > controls.step.maximum_surfaces {
        return Err(GasChamberAdvanceError::Invalid("finite nonnegative duration and bounded positive error/work controls required"));
    }
    let mut templates = surfaces.to_vec();
    templates.sort_by_key(|surface| surface.id);
    if templates.windows(2).any(|pair| pair[0].id == pair[1].id) {
        return Err(GasChamberAdvanceError::Invalid("duplicate surface identity"));
    }
    if let Some(laws) = accommodation {
        if laws.len() != templates.len() { return Err(GasChamberAdvanceError::Invalid("accommodation surface roster mismatch")); }
        for surface in &templates {
            poll(&mut cancelled)?;
            if !laws.contains_key(&surface.id) { return Err(GasChamberAdvanceError::Invalid("missing surface accommodation")); }
        }
    }
    if let Some(calorics) = calorics {
        phase_inventory::validate_roster(&templates, calorics)?;
        if geometry.is_some() {
            return Err(GasChamberAdvanceError::Invalid("mixed inventories cannot infer pure-body moving geometry"));
        }
    }
    if let Some(geometry) = geometry { geometry.validate(&templates)?; }
    let mut states: Vec<_> = templates.iter().map(|s| s.solid).collect();
    let mut sums: Vec<_> = (0..templates.len()).map(|_| SurfaceSums::default()).collect();
    let mut gas_energy = Sum::default();
    let mut gas_stream_residual = Sum::default();
    let mut endpoint: Option<GasChamberCandidate<'catalog>> = None;
    let mut elapsed = 0.0;
    // Sum accepted physical durations rather than accumulating clock roundoff
    // into an artificial final interval whose mass transfer cannot be represented.
    let mut elapsed_sum = Sum::default();
    let mut dt = controls.maximum_dt_s.min(duration_s);
    let mut trials = 0;
    let mut rejected = 0;
    let mut accepted = 0;
    let mut evaluations = 0;
    let mut maximum_error: f64 = 0.0;
    while elapsed < duration_s {
        poll(&mut cancelled)?;
        if trials == controls.maximum_trials {
            return Err(GasChamberAdvanceError::TrialBudget { maximum_trials: controls.maximum_trials });
        }
        trials += 1;
        let remaining = duration_s - elapsed;
        dt = dt.min(remaining);
        if !positive(dt) || elapsed + dt == elapsed || dt * 0.5 == 0.0 {
            return Err(GasChamberAdvanceError::MinimumStep { requested_dt_s: dt, minimum_dt_s: controls.minimum_dt_s });
        }
        let before = with_states(&templates, &states);
        let current_gas = endpoint.as_ref().map_or(gas, |value| value.gas());
        let mut evaluate = |rows: &[GasChamberSurface<'_>], gas: &ClosedGasParcel, duration| {
            evaluations += 1;
            evaluate_interval(rows, accommodation, geometry, calorics, gas, catalog, duration, controls.step, &mut cancelled)
        };
        let trial = (|| {
            // A coarse trial is an error probe, never an accepted heat history.
            let (coarse, _) = evaluate(&before, current_gas, dt)?;
            let (first, first_heat) = evaluate(&before, current_gas, dt * 0.5)?;
            let middle: Vec<_> = first.surfaces().iter().map(|row| row.solid).collect();
            let (second, second_heat) = evaluate(&with_states(&templates, &middle), first.gas(), dt * 0.5)?;
            Ok::<_, GasSurfaceError>((coarse, first, first_heat, second, second_heat))
        })();
        let (coarse, first, first_heat, second, second_heat) = match trial {
            Ok(value) => value,
            Err(GasSurfaceError::InventoryLimit { maximum_duration_s }) => {
                rejected += 1;
                dt = smaller_step(0.8 * maximum_duration_s, dt, controls.minimum_dt_s)?;
                continue;
            }
            Err(error) => return Err(error.into()),
        };
        let mut ratio = compare(&before, current_gas, catalog, &coarse, &second, controls, &mut cancelled)?;
        if let Some(calorics) = calorics {
            ratio = ratio.max(partitioned::compare(&before, calorics, &coarse, &second, controls, &mut cancelled)?);
        }
        if let Some(geometry) = geometry {
            ratio = ratio.max(geometry.compare(&before, &coarse, &second, controls.relative_tolerance, &mut cancelled)?);
        }
        if ratio > 1.0 {
            rejected += 1;
            dt = smaller_step(dt * 0.8 / fs_math::det::sqrt(ratio), dt, controls.minimum_dt_s)?;
            continue;
        }
        for ((sum, a), b) in sums.iter_mut().zip(first.surfaces()).zip(second.surfaces()) {
            poll(&mut cancelled)?;
            sum.add(a, first_heat.get(&a.id))?;
            sum.add(b, second_heat.get(&b.id))?;
        }
        // Only the two accepted physical endpoints belong to this history.
        // Their source-relative caloric defects telescope without subtracting
        // rounded absolute formation energies at the window endpoints.
        for part in [&first, &second] {
            gas_energy.add(part.step().gas_energy_change_j)?;
            gas_stream_residual.add(part.gas_exchange().step().energy_residual_j)?;
        }
        accepted += 1;
        maximum_error = maximum_error.max(ratio);
        for (state, row) in states.iter_mut().zip(second.surfaces()) { *state = row.solid; }
        endpoint = Some(second);
        elapsed_sum.add(dt)?;
        elapsed = if dt == remaining { duration_s } else { elapsed_sum.value };
        let factor = if ratio == 0.0 { 2.0 } else { (0.8 / fs_math::det::sqrt(ratio)).clamp(0.5, 2.0) };
        dt = (dt * factor).min(controls.maximum_dt_s).max(controls.minimum_dt_s);
    }
    let endpoint = match endpoint {
        Some(endpoint) => endpoint,
        None => {
            evaluations += 1;
            evaluate_interval(&templates, accommodation, geometry, calorics, gas, catalog, 0.0, controls.step, &mut cancelled)?.0
        }
    };
    let mut solid_energy = Sum::default();
    let mut solid_stream_residual = Sum::default();
    let mut solid_mass = Sum::default();
    let mut stream_energy = Sum::default();
    let mut reflected_energy = Sum::default();
    let mut species_transfers = BTreeMap::new();
    let mut histories = Vec::with_capacity(templates.len());
    for ((original, final_row), sum) in templates.iter().zip(endpoint.surfaces()).zip(sums) {
        poll(&mut cancelled)?;
        let start = original.solid;
        let end = final_row.solid;
        solid_energy.add(sum.stored.value)?;
        solid_stream_residual.add(sum.residual.value)?;
        solid_mass.add(end.mass_kg - start.mass_kg)?;
        stream_energy.add(sum.energy.value)?;
        reflected_energy.add(sum.reflected.value)?;
        species_transfers.entry(original.binding.species.clone()).or_insert_with(Sum::default).add(sum.transferred.value)?;
        histories.push(GasSurfaceHistory {
            id: original.id, outgoing_mass_kg: sum.outgoing.value, incoming_mass_kg: sum.incoming.value,
            transferred_mass_kg: sum.transferred.value, energy_into_gas_j: sum.energy.value,
            reflected_heat_into_gas_j: sum.reflected.value, recession_m: sum.recession.value,
            surface_temperature_k: [sum.initial_temperature.unwrap_or(final_row.surface_temperature_k[0]), final_row.surface_temperature_k[1]],
        });
    }
    let gas_energy = gas_energy.value;
    let gas_stream_residual = gas_stream_residual.value;
    let solid_stream_residual = solid_stream_residual.value;
    let energy_residual = solid_stream_residual + gas_stream_residual;
    let mass_residual = solid_mass.value + (endpoint.gas().mass_kg() - gas.mass_kg());
    let mut maximum_species_mass_residual: f64 = 0.0;
    for species in catalog.thermodynamics().components() {
        poll(&mut cancelled)?;
        let id = species.species();
        let change = endpoint.gas().species_mass_kg(id).map_err(GasSurfaceError::from)?
            - gas.species_mass_kg(id).map_err(GasSurfaceError::from)?;
        let transported = species_transfers.get(id).map_or(0.0, |sum| sum.value);
        let residual = change - transported;
        if !residual.is_finite() { return Err(GasChamberAdvanceError::Invalid("whole-history species balance overflow")); }
        maximum_species_mass_residual = maximum_species_mass_residual.max(residual.abs());
    }
    let energy_tolerance = controls.step.surface.gas.energy_tolerance_j;
    if ![gas_energy, energy_residual, gas_stream_residual, solid_stream_residual].iter().all(|x| x.is_finite())
        || energy_residual.abs().max(gas_stream_residual.abs()).max(solid_stream_residual.abs()) > energy_tolerance {
        return Err(GasSurfaceError::EnergyBalance {
            residual_j: energy_residual.abs().max(gas_stream_residual.abs()).max(solid_stream_residual.abs()),
            tolerance_j: energy_tolerance,
        }.into());
    }
    if !mass_residual.is_finite() || mass_residual.abs() > controls.step.surface.gas.mass_tolerance_kg
        || maximum_species_mass_residual > controls.step.surface.gas.mass_tolerance_kg {
        return Err(GasChamberAdvanceError::Invalid("whole-history solid/gas mass balance failed"));
    }
    poll(&mut cancelled)?;
    Ok(GasChamberAdvanceCandidate { endpoint, audit: GasChamberAdvanceAudit {
        elapsed_s: elapsed, trials, rejected_trials: rejected, accepted_intervals: accepted,
        kernel_evaluations: evaluations, maximum_accepted_error_ratio: maximum_error, surfaces: histories,
        energy_into_gas_j: stream_energy.value, reflected_heat_into_gas_j: reflected_energy.value,
        solid_energy_change_j: solid_energy.value, gas_energy_change_j: gas_energy,
        energy_residual_j: energy_residual, gas_stream_energy_residual_j: gas_stream_residual,
        solid_stream_energy_residual_j: solid_stream_residual, mass_residual_kg: mass_residual,
        maximum_species_mass_residual_kg: maximum_species_mass_residual,
    } })
}

#[cfg(test)]
mod tests;

#[cfg(test)]
#[path = "advance/accommodation_tests.rs"]
mod accommodation_tests;
