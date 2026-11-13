//! Adaptive sublimation/deposition and recoil through the original source step.
//!
//! Compare one full interval with two half intervals, accepting only the latter.
//! Every trial reevaluates the original vapor pressure, molecular distributions,
//! caloric sources and both momenta. Conservation tolerances are separate from
//! this first-order LOCAL truncation estimate; neither is a global error bound.
//! Supplied heat is constant power over the requested window, counted once on
//! accepted steps. No rejected trial contributes mass, heat, recoil or history.

use super::kinetic::{
    KineticSublimationCandidate, KineticSublimationControls, KineticSublimationError,
    KineticSublimationStep, MovingSolidSurface, advance_solid_with_moving_gas,
};
use crate::liquid_gas::evaporation::momentum::MovingGasParcel;
use fs_material::sublimation::gas::{
    chamber::GasChamberSurface,
    momentum::{MaxwellReflection, MolecularSurfaceFrame},
};
use fs_material::thermochemical_gas::{ThermochemicalGasModel, energy::GasEnergyError};
use fs_math::dd::Dd;

/// Physical error scales and bounded work. The native controls retain the
/// maximum physical step, donor-inventory and temperature-increment limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SublimationAdaptiveControls {
    pub initial_step_s: f64,
    /// Minimum RETRY interval. A shorter final remainder may be tried once.
    pub minimum_step_s: f64,
    /// Every attempted full/two-half comparison spends one attempt, even when
    /// an original inventory limit refuses before all three calls complete.
    pub maximum_attempts: usize,
    /// Each accepted interval retains exactly two original source reports.
    pub maximum_accepted_intervals: usize,
    /// Scales changes during the trial, never a formation datum or common drift.
    pub relative_tolerance: f64,
    pub mass_tolerance_kg: f64,
    pub temperature_tolerance_k: f64,
    pub pressure_tolerance_pa: f64,
    pub specific_enthalpy_tolerance_j_kg: f64,
    /// Per-component error in gas velocity relative to the solid.
    pub velocity_tolerance_m_s: f64,
}

#[derive(Debug)]
pub enum SublimationTrajectoryError {
    Invalid(&'static str),
    Step { elapsed_s: f64, error: KineticSublimationError },
    Inventory { elapsed_s: f64, error: GasEnergyError },
    MinimumStep { elapsed_s: f64, requested_s: f64, minimum_s: f64 },
    WorkBudget { elapsed_s: f64, attempts: usize, accepted_intervals: usize },
    Cancelled { elapsed_s: f64 },
}
impl core::fmt::Display for SublimationTrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive solid/gas recoil refused: {self:?}")
    }
}
impl std::error::Error for SublimationTrajectoryError {}
type Error = SublimationTrajectoryError;

/// Only accepted native half steps enter this chronological physical history.
#[derive(Clone, Debug, PartialEq)]
pub struct SublimationTrajectoryStep {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub rejected_intervals: usize,
    pub accepted_intervals: usize,
    /// Includes the initial zero-time admission and rejected/coarse work.
    pub source_evaluations: usize,
    pub maximum_accepted_error_ratio: f64,
    pub accepted_steps: Vec<KineticSublimationStep>,
    pub supplied_heat_j: f64,
    pub outgoing_mass_kg: f64,
    pub incoming_mass_kg: f64,
    pub transferred_mass_kg: f64,
    pub recession_m: f64,
    /// Independent actual endpoint solid-plus-gas mass change.
    pub mass_residual_kg: f64,
    /// Independent actual endpoint solid-plus-gas momentum change.
    pub momentum_residual_n_s: [f64; 3],
    /// Compensated sum of accepted original LABORATORY-frame residuals.
    /// Local translating-frame residuals cannot be added without a crosswalk.
    pub accumulated_energy_residual_j: f64,
}

/// Both complete-window endpoints publish atomically. No partial trajectory
/// is returned when a later source, accuracy, budget or cancellation check fails.
#[derive(Clone, Debug, PartialEq)]
pub struct SublimationTrajectoryCandidate {
    solid: MovingSolidSurface,
    gas: MovingGasParcel,
    step: SublimationTrajectoryStep,
}
impl SublimationTrajectoryCandidate {
    pub const fn solid(&self) -> &MovingSolidSurface { &self.solid }
    pub const fn gas(&self) -> &MovingGasParcel { &self.gas }
    pub const fn step(&self) -> &SublimationTrajectoryStep { &self.step }
    pub fn into_parts(self) -> (MovingSolidSurface, MovingGasParcel, SublimationTrajectoryStep) {
        (self.solid, self.gas, self.step)
    }
}

fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }
fn finite(value: Dd) -> Result<f64, Error> {
    let result = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && result.is_finite() {
        Ok(result)
    } else {
        Err(Error::Invalid("unrepresentable adaptive state or account"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool, elapsed_s: f64) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled { elapsed_s }) } else { Ok(()) }
}
fn physical(error: KineticSublimationError, elapsed_s: f64) -> Error {
    if matches!(error, KineticSublimationError::Cancelled) {
        Error::Cancelled { elapsed_s }
    } else {
        Error::Step { elapsed_s, error }
    }
}

// Increment scaling preserves accuracy near equilibrium without granting a
// large allowance merely because a carrier inventory or energy datum is large.
fn ratio(old: f64, full: f64, fine: f64, absolute: f64, relative: f64) -> Result<f64, Error> {
    let full_change = finite(Dd::from_f64(full) - Dd::from_f64(old))?;
    let fine_change = finite(Dd::from_f64(fine) - Dd::from_f64(old))?;
    let scale = absolute + relative * full_change.abs().max(fine_change.abs());
    let difference = finite(Dd::from_f64(fine) - Dd::from_f64(full))?;
    let value = difference.abs() / scale;
    if !positive(scale) || !value.is_finite() {
        Err(Error::Invalid("unrepresentable adaptive error allowance"))
    } else { Ok(value) }
}

fn relative_velocity(solid: &MovingSolidSurface, gas: &MovingGasParcel) -> Result<[f64; 3], Error> {
    let mut result = [0.0; 3];
    for (d, value) in result.iter_mut().enumerate() {
        *value = finite(
            Dd::from_f64(gas.momentum_kg_m_s()[d]) / Dd::from_f64(gas.parcel().mass_kg())
                - Dd::from_f64(solid.momentum_kg_m_s()[d]) / Dd::from_f64(solid.solid().mass_kg),
        )?;
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn compare(
    before_solid: &MovingSolidSurface,
    before_gas: &MovingGasParcel,
    full: &KineticSublimationCandidate,
    fine: &KineticSublimationCandidate,
    catalog: &ThermochemicalGasModel,
    controls: SublimationAdaptiveControls,
    elapsed_s: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let mut error = 0.0_f64;
    let old = before_solid.solid();
    let a = full.solid().solid();
    let b = fine.solid().solid();
    for (old, full, fine, absolute) in [
        (old.mass_kg, a.mass_kg, b.mass_kg, controls.mass_tolerance_kg),
        (old.specific_enthalpy_j_kg, a.specific_enthalpy_j_kg, b.specific_enthalpy_j_kg,
            controls.specific_enthalpy_tolerance_j_kg),
        (full.step().surface_temperature_k[0], full.step().surface_temperature_k[1],
            fine.step().surface_temperature_k[1], controls.temperature_tolerance_k),
    ] {
        error = error.max(ratio(old, full, fine, absolute, controls.relative_tolerance)?);
    }
    let [old, a, b] = [before_gas, full.gas(), fine.gas()].map(|gas| gas.parcel().state().state());
    for (old, full, fine, absolute) in [
        (old.temperature, a.temperature, b.temperature, controls.temperature_tolerance_k),
        (old.pressure, a.pressure, b.pressure, controls.pressure_tolerance_pa),
    ] {
        error = error.max(ratio(old, full, fine, absolute, controls.relative_tolerance)?);
    }
    let old = relative_velocity(before_solid, before_gas)?;
    let a = relative_velocity(full.solid(), full.gas())?;
    let b = relative_velocity(fine.solid(), fine.gas())?;
    for d in 0..3 {
        error = error.max(ratio(old[d], a[d], b[d], controls.velocity_tolerance_m_s,
            controls.relative_tolerance)?);
    }
    // Query each original composition once, not once per retained species.
    let inventory = |gas: &MovingGasParcel| gas.parcel().species_inventory_kg()
        .map_err(|error| Error::Inventory { elapsed_s, error });
    let old = inventory(before_gas)?;
    let a = inventory(full.gas())?;
    let b = inventory(fine.gas())?;
    for source in catalog.thermodynamics().components() {
        poll(cancelled, elapsed_s)?;
        let id = source.species();
        error = error.max(ratio(*old.get(id).unwrap_or(&0.0), *a.get(id).unwrap_or(&0.0),
            *b.get(id).unwrap_or(&0.0), controls.mass_tolerance_kg, controls.relative_tolerance)?);
    }
    Ok(error)
}

fn retry_limit(error: &KineticSublimationError) -> Option<f64> {
    match error {
        KineticSublimationError::InventoryLimit { maximum_duration_s }
        | KineticSublimationError::Duration { maximum_duration_s }
        | KineticSublimationError::TemperatureIncrement { maximum_duration_s } => Some(*maximum_duration_s),
        // Source domains, caloric compatibility, arithmetic and conservation
        // failures are not recast as temporal error or silently retried.
        _ => None,
    }
}
fn smaller(old: f64, requested: f64, minimum_s: f64, elapsed_s: f64) -> Result<f64, Error> {
    // Try the admissible minimum once even when a conservative controller
    // proposal undershoots it. The native source and error check still decide
    // whether that interval is physically and numerically acceptable.
    let bounded = requested.max(minimum_s);
    if !requested.is_finite() || requested < 0.0 || old <= minimum_s || bounded >= old {
        Err(Error::MinimumStep { elapsed_s, requested_s: requested, minimum_s })
    } else { Ok(bounded) }
}

/// Evolve a finite pure solid and the original gas through one complete window.
///
/// `supplied_heat_w` is constant signed power into the SOLID. Its integral is
/// applied only through the original accepted source transactions. Source and
/// geometry definitions remain fixed, while both temperatures, backpressure,
/// composition and bulk velocities are reevaluated after every native half step.
/// Inventory and temperature-increment limits request a smaller physical step;
/// constitutive, conservation and cancellation failures remain named refusals.
///
/// This controls the temporal error of the uniform fixed-reference source
/// model. It does not advance pose or a spatial gas field. The local ratio uses
/// both temperatures, pressure, every retained species, solid h and gas/solid
/// relative drift. Reflected carrier loads therefore participate even when net
/// phase mass is zero. No Richardson extrapolation can overshoot inventories.
#[allow(clippy::too_many_arguments)]
pub fn advance_sublimating_trajectory(
    surface: &GasChamberSurface<'_>,
    solid: &MovingSolidSurface,
    gas: &MovingGasParcel,
    catalog: &ThermochemicalGasModel,
    frame: MolecularSurfaceFrame,
    pure_solid_enthalpy_range_j_kg: [f64; 2],
    scattering: Option<&MaxwellReflection>,
    duration_s: f64,
    supplied_heat_w: f64,
    native: KineticSublimationControls,
    adaptive: SublimationAdaptiveControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SublimationTrajectoryCandidate, Error> {
    poll(&mut cancelled, 0.0)?;
    if !duration_s.is_finite() || duration_s < 0.0 || !supplied_heat_w.is_finite()
        || !positive(adaptive.initial_step_s) || !positive(adaptive.minimum_step_s)
        || adaptive.initial_step_s < adaptive.minimum_step_s
        || !positive(native.maximum_duration_s) || native.maximum_duration_s < adaptive.minimum_step_s
        || adaptive.maximum_attempts == 0 || adaptive.maximum_accepted_intervals == 0
        || adaptive.maximum_accepted_intervals.checked_mul(2).is_none()
        || !adaptive.relative_tolerance.is_finite() || !(0.0..1.0).contains(&adaptive.relative_tolerance)
        || ![adaptive.mass_tolerance_kg, adaptive.temperature_tolerance_k, adaptive.pressure_tolerance_pa,
            adaptive.specific_enthalpy_tolerance_j_kg, adaptive.velocity_tolerance_m_s]
            .into_iter().all(positive)
    {
        return Err(Error::Invalid("invalid sublimation trajectory controls, duration or heat power"));
    }
    // A zero-time request still admits the original source, state and controls.
    let admission = advance_solid_with_moving_gas(surface, solid, gas, catalog, frame,
        pure_solid_enthalpy_range_j_kg, scattering, 0.0, 0.0, native, &mut cancelled)
        .map_err(|error| physical(error, 0.0))?;
    let (mut current_solid, mut current_gas, _) = admission.into_parts();
    let mut elapsed = 0.0;
    let mut next_dt = adaptive.initial_step_s.min(native.maximum_duration_s);
    let mut report = SublimationTrajectoryStep {
        duration_s, attempted_intervals: 0, rejected_intervals: 0, accepted_intervals: 0,
        source_evaluations: 1, maximum_accepted_error_ratio: 0.0, accepted_steps: Vec::new(),
        supplied_heat_j: 0.0, outgoing_mass_kg: 0.0, incoming_mass_kg: 0.0,
        transferred_mass_kg: 0.0, recession_m: 0.0, mass_residual_kg: 0.0,
        momentum_residual_n_s: [0.0; 3], accumulated_energy_residual_j: 0.0,
    };
    while elapsed < duration_s {
        poll(&mut cancelled, elapsed)?;
        if report.attempted_intervals >= adaptive.maximum_attempts
            || report.accepted_intervals >= adaptive.maximum_accepted_intervals
        {
            return Err(Error::WorkBudget { elapsed_s: elapsed, attempts: report.attempted_intervals,
                accepted_intervals: report.accepted_intervals });
        }
        let remaining = duration_s - elapsed;
        let proposed = next_dt.min(remaining).min(native.maximum_duration_s);
        let mut end = if proposed == remaining { duration_s } else { elapsed + proposed };
        // Rounding the absolute clock upward must not exceed a native duration
        // cap by one ulp. Keep the represented interval within its proposal.
        if end - elapsed > proposed { end = end.next_down(); }
        let dt = end - elapsed;
        if !positive(dt) || dt > native.maximum_duration_s || 0.5 * dt == 0.0 {
            return Err(Error::Invalid("adaptive physical clock cannot represent its next interval"));
        }
        report.attempted_intervals += 1;
        let mut run = |solid: &MovingSolidSurface, gas: &MovingGasParcel, dt: f64| {
            report.source_evaluations = report.source_evaluations.checked_add(1)
                .ok_or(KineticSublimationError::Invalid("adaptive source work counter overflow"))?;
            let heat = supplied_heat_w * dt;
            if !heat.is_finite() || (supplied_heat_w != 0.0 && dt != 0.0 && heat == 0.0) {
                return Err(KineticSublimationError::Invalid("unrepresentable supplied heat increment"));
            }
            let source = GasChamberSurface { id: surface.id, law: surface.law, binding: surface.binding,
                solid: solid.solid(), temperature_from_enthalpy: surface.temperature_from_enthalpy };
            advance_solid_with_moving_gas(&source, solid, gas, catalog, frame,
                pure_solid_enthalpy_range_j_kg, scattering, dt, heat, native, &mut cancelled)
        };
        let trial = (|| {
            let full = run(&current_solid, &current_gas, dt)?;
            let first = run(&current_solid, &current_gas, 0.5 * dt)?;
            let second = run(first.solid(), first.gas(), dt - 0.5 * dt)?;
            Ok::<_, KineticSublimationError>((full, first, second))
        })();
        let (full, first, second) = match trial {
            Ok(trial) => trial,
            Err(error) => {
                if let Some(limit) = retry_limit(&error) {
                    report.rejected_intervals += 1;
                    next_dt = smaller(dt, (0.5 * dt).min(0.9 * limit), adaptive.minimum_step_s, elapsed)?;
                    continue;
                }
                return Err(physical(error, elapsed));
            }
        };
        let error = compare(&current_solid, &current_gas, &full, &second, catalog, adaptive,
            elapsed, &mut cancelled)?;
        if error > 1.0 {
            report.rejected_intervals += 1;
            let factor = (0.9 / fs_math::det::sqrt(error)).clamp(0.1, 0.5);
            next_dt = smaller(dt, dt * factor, adaptive.minimum_step_s, elapsed)?;
            continue;
        }
        poll(&mut cancelled, elapsed)?;
        report.accepted_steps.try_reserve(2)
            .map_err(|_| Error::Invalid("accepted sublimation history allocation refused"))?;
        let (_, _, first_step) = first.into_parts();
        let (solid, gas, second_step) = second.into_parts();
        report.accepted_steps.extend([first_step, second_step]);
        current_solid = solid;
        current_gas = gas;
        report.accepted_intervals += 1;
        report.maximum_accepted_error_ratio = report.maximum_accepted_error_ratio.max(error);
        elapsed = end;
        // A bounded binary growth rule is deterministic and avoids large
        // controller jumps immediately after an inventory or accuracy rejection.
        next_dt = if error < 0.125 { (2.0 * dt).min(native.maximum_duration_s) } else { dt };
    }
    let mut sums = [Dd::ZERO; 6];
    for step in &report.accepted_steps {
        poll(&mut cancelled, elapsed)?;
        for (sum, value) in sums.iter_mut().zip([step.supplied_heat_j, step.outgoing_mass_kg,
            step.incoming_mass_kg, step.transferred_mass_kg, step.recession_m,
            step.laboratory_energy_residual_j]) {
            *sum = *sum + Dd::from_f64(value);
        }
    }
    report.supplied_heat_j = finite(sums[0])?;
    report.outgoing_mass_kg = finite(sums[1])?;
    report.incoming_mass_kg = finite(sums[2])?;
    report.transferred_mass_kg = finite(sums[3])?;
    report.recession_m = finite(sums[4])?;
    report.accumulated_energy_residual_j = finite(sums[5])?;
    report.mass_residual_kg = finite(
        (Dd::from_f64(current_solid.solid().mass_kg) - Dd::from_f64(solid.solid().mass_kg))
            + (Dd::from_f64(current_gas.parcel().mass_kg()) - Dd::from_f64(gas.parcel().mass_kg())),
    )?;
    for d in 0..3 {
        report.momentum_residual_n_s[d] = finite(
            (Dd::from_f64(current_solid.momentum_kg_m_s()[d]) - Dd::from_f64(solid.momentum_kg_m_s()[d]))
                + (Dd::from_f64(current_gas.momentum_kg_m_s()[d]) - Dd::from_f64(gas.momentum_kg_m_s()[d])),
        )?;
    }
    poll(&mut cancelled, elapsed)?;
    Ok(SublimationTrajectoryCandidate { solid: current_solid, gas: current_gas, step: report })
}
