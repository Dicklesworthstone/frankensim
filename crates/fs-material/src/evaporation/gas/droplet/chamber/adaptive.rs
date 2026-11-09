//! Adaptive physical time for original collective spherical evaporation.
//!
//! One original full interval is compared with two original fine intervals.
//! Only the fine pair is accepted. Each interval reconstructs spherical area
//! and transfer from its actual input inventories; source adapters can refresh
//! molecular and momentum properties through the same controller. Coarse and
//! rejected work never contributes physical transfers or final source payload.
//!
//! This estimates local first-order time error in each droplet's mass, caloric
//! coordinate, temperature and radius, and in the shared gas temperature and
//! every original species inventory. Relative scales use changes from the
//! interval start, never absolute formation enthalpy or Kelvin temperature.
//! The existing method-order arithmetic accounts for unequal representable
//! half spans. The conservative fine endpoint is not Richardson-extrapolated.
//!
//! A complete window publishes atomically under explicit comparison, reduction
//! and conservation budgets. Positive droplet seeds and the original phase
//! admissibility remain mandatory. Complete dry-out, nucleation, trajectories,
//! shielding and moving-volume work are not inferred; a step floor or exhausted
//! work allowance returns refusal, never a clipped physical state. The local
//! estimator is not a certified bound on accumulated trajectory error.

/// Complete gas-network accuracy and accepted physical-window accounting.
pub mod network;

mod account;
pub use account::{
    GasDropletEvolutionAudit, GasDropletEvolutionSpeciesAudit, GasDropletEvolutionSurfaceAudit,
};

use super::{
    ClosedGasParcel, EvaporationError, EvaporationLaw, GasDropletChamberCandidate,
    GasDropletHeatSurface, GasDropletSurface, GasDropletSurfaceStep, GasEvaporationError,
    GasLiquidBinding, GasLiquidChamberHeatSolve, GasLiquidChamberSolve, SphericalDropletGeometry,
    SphericalDropletState, ThermochemicalGasModel, advance_droplets_with_gas,
    advance_droplets_with_gas_and_heat, poll,
};
use fs_thermal::enthalpy::conduction::implicit::adaptive::ConductionTimeOrder;

/// Physical identities and original calorics, independent of film selection.
/// A source-backed stepper supplies its own D/Sh rather than inventing numerical
/// placeholders merely to participate in the common time controller.
#[derive(Clone, Copy)]
pub struct GasDropletSeed<'a> {
    pub id: u64,
    pub law: &'a dyn EvaporationLaw,
    pub binding: &'a GasLiquidBinding,
    pub droplet: SphericalDropletState,
    pub temperature_from_enthalpy: &'a dyn Fn(f64) -> Result<f64, EvaporationError>,
}
impl<'a> From<GasDropletSurface<'a>> for GasDropletSeed<'a> {
    fn from(surface: GasDropletSurface<'a>) -> Self {
        Self {
            id: surface.id,
            law: surface.law,
            binding: surface.binding,
            droplet: surface.droplet,
            temperature_from_enthalpy: surface.temperature_from_enthalpy,
        }
    }
}

/// Local time accuracy, original per-interval physics and one window budget.
/// Step bounds apply to the full comparison span; its fine intervals are shorter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveGasDropletControls {
    pub solve: GasLiquidChamberSolve,
    pub initial_dt_s: f64,
    pub minimum_dt_s: f64,
    pub maximum_dt_s: f64,
    /// Includes coarse/fine accuracy rejections and original step-size refusals.
    /// Each comparison invokes at most three bounded original chamber solves;
    /// one additional zero-duration solve admits the original physical state.
    pub maximum_comparisons: usize,
    pub maximum_reductions: usize,
    pub enthalpy_absolute_tolerance_j_kg: f64,
    pub temperature_absolute_tolerance_k: f64,
    pub mass_absolute_tolerance_kg: f64,
    pub radius_absolute_tolerance_m: f64,
    /// Applies to each coordinate's CHANGE from the accepted comparison start.
    pub relative_tolerance: f64,
    /// Noncancelling accepted defects share these complete-window allowances.
    pub mass_conservation_tolerance_kg: f64,
    pub energy_conservation_tolerance_j: f64,
}
impl AdaptiveGasDropletControls {
    /// Numerical defaults only; all original physical inputs remain explicit.
    pub fn new(solve: GasLiquidChamberSolve) -> Self {
        Self {
            solve,
            initial_dt_s: 0.1,
            minimum_dt_s: 1e-10,
            maximum_dt_s: 1.0,
            maximum_comparisons: 4096,
            maximum_reductions: 32,
            enthalpy_absolute_tolerance_j_kg: 1.0,
            temperature_absolute_tolerance_k: 1e-3,
            mass_absolute_tolerance_kg: 1e-12,
            radius_absolute_tolerance_m: 1e-9,
            relative_tolerance: 1e-3,
            mass_conservation_tolerance_kg: solve.surface.gas.mass_tolerance_kg,
            energy_conservation_tolerance_j: solve.surface.gas.energy_tolerance_j,
        }
    }

    /// Admit numerical time, accuracy and roster controls without querying laws.
    /// Original interval owners still admit their physical solver/source inputs.
    pub fn validate(self, duration_s: f64, count: usize) -> Result<(), GasEvaporationError> {
        if !duration_s.is_finite()
            || duration_s < 0.0
            || self.solve.maximum_surfaces == 0
            || count > self.solve.maximum_surfaces
            || self.maximum_comparisons == 0
            || [
                self.initial_dt_s,
                self.minimum_dt_s,
                self.maximum_dt_s,
                self.enthalpy_absolute_tolerance_j_kg,
                self.temperature_absolute_tolerance_k,
                self.mass_absolute_tolerance_kg,
                self.radius_absolute_tolerance_m,
            ]
            .into_iter()
            .any(|x| !x.is_finite() || x <= 0.0)
            || self.minimum_dt_s > self.maximum_dt_s
            || self.initial_dt_s < self.minimum_dt_s
            || !self.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_tolerance)
            || [
                self.mass_conservation_tolerance_kg,
                self.energy_conservation_tolerance_j,
            ]
            .into_iter()
            .any(|x| !x.is_finite() || x < 0.0)
        {
            return Err(GasEvaporationError::Invalid(
                "invalid adaptive spherical time, roster, accuracy or conservation controls",
            ));
        }
        Ok(())
    }
}

/// Actual coarse/fine differences and their first-order local-error estimate.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct GasDropletTimeError {
    pub maximum_enthalpy_difference_j_kg: f64,
    pub maximum_temperature_difference_k: f64,
    pub maximum_liquid_mass_difference_kg: f64,
    pub maximum_radius_difference_m: f64,
    pub gas_temperature_difference_k: f64,
    pub maximum_species_mass_difference_kg: f64,
    pub normalized: f64,
}

/// Physical totals describe only accepted fine intervals; work counts every try.
#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveGasDropletReport {
    pub duration_s: f64,
    pub exchange: GasDropletEvolutionAudit,
    pub attempted_steps: usize,
    pub comparisons: usize,
    pub accepted_comparisons: usize,
    pub accepted_intervals: usize,
    pub accuracy_rejections: usize,
    pub physical_rejections: usize,
    pub maximum_accepted_normalized_error: f64,
    pub last_error: Option<GasDropletTimeError>,
    /// Explicit next-call controller state; no hidden PI history is retained.
    pub suggested_dt_s: f64,
}

/// Original step failures keep their owning source error and attempted clock.
#[derive(Debug)]
pub enum AdaptiveGasDropletError<E = GasEvaporationError> {
    Step {
        time_s: f64,
        duration_s: f64,
        error: E,
    },
    Material(GasEvaporationError),
    Invalid(&'static str),
    ComparisonBudget,
    TimeAccuracy {
        time_s: f64,
        span_s: f64,
        error: GasDropletTimeError,
    },
}
impl<E> From<GasEvaporationError> for AdaptiveGasDropletError<E> {
    fn from(error: GasEvaporationError) -> Self {
        Self::Material(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for AdaptiveGasDropletError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive collective spherical evaporation: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for AdaptiveGasDropletError<E> {}

/// Complete fine-trajectory endpoint, cumulative physical audit and last payload.
#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveGasDropletCandidate<'catalog, T = ()> {
    endpoint: GasDropletChamberCandidate<'catalog>,
    payload: T,
    report: AdaptiveGasDropletReport,
}
impl<'catalog, T> AdaptiveGasDropletCandidate<'catalog, T> {
    /// Last accepted original fine interval, or the zero-duration admission.
    /// Its local receipts are not the complete-window transfer totals.
    pub const fn endpoint(&self) -> &GasDropletChamberCandidate<'catalog> {
        &self.endpoint
    }
    /// Accepted final droplets. Geometry/receipts here describe the last fine
    /// interval; the cumulative report retains the full original-to-final audit.
    pub fn droplets(&self) -> &[GasDropletSurfaceStep] {
        self.endpoint.droplets()
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.endpoint.gas()
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.endpoint.catalog()
    }
    pub const fn payload(&self) -> &T {
        &self.payload
    }
    pub const fn report(&self) -> &AdaptiveGasDropletReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        GasDropletChamberCandidate<'catalog>,
        T,
        AdaptiveGasDropletReport,
    ) {
        (self.endpoint, self.payload, self.report)
    }
}

/// Evolve positive droplets while refreshing geometry at every original step.
/// Molecular D and explicit Sh remain those authored in each input row.
pub fn advance_droplets_with_gas_adaptive<'catalog>(
    surfaces: &[GasDropletSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: AdaptiveGasDropletControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AdaptiveGasDropletCandidate<'catalog>, AdaptiveGasDropletError> {
    poll(&mut cancelled)?;
    controls.validate(duration_s, surfaces.len())?;
    let mut rows = surfaces.to_vec();
    rows.sort_unstable_by_key(|surface| surface.id);
    let seeds: Vec<_> = rows.iter().copied().map(GasDropletSeed::from).collect();
    advance_droplets_with_gas_adaptive_using(
        &seeds,
        gas,
        catalog,
        duration_s,
        controls,
        |current, gas, dt, check| {
            for (row, current) in rows.iter_mut().zip(current) {
                row.droplet = current.droplet;
            }
            advance_droplets_with_gas(&rows, gas, catalog, dt, controls.solve, check)
                .map(|candidate| (candidate, ()))
        },
        material_step_limit,
        cancelled,
    )
}

/// Refine the complete original heat-then-phase split with evolving geometry.
#[allow(clippy::too_many_arguments)]
pub fn advance_droplets_with_gas_and_heat_adaptive<'catalog>(
    surfaces: &[GasDropletHeatSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: AdaptiveGasDropletControls,
    heat_controls: GasLiquidChamberHeatSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AdaptiveGasDropletCandidate<'catalog>, AdaptiveGasDropletError> {
    poll(&mut cancelled)?;
    controls.validate(duration_s, surfaces.len())?;
    let mut rows = surfaces.to_vec();
    rows.sort_unstable_by_key(|row| row.surface.id);
    let seeds: Vec<_> = rows
        .iter()
        .map(|row| GasDropletSeed::from(row.surface))
        .collect();
    advance_droplets_with_gas_adaptive_using(
        &seeds,
        gas,
        catalog,
        duration_s,
        controls,
        |current, gas, dt, check| {
            for (row, current) in rows.iter_mut().zip(current) {
                row.surface.droplet = current.droplet;
            }
            advance_droplets_with_gas_and_heat(
                &rows,
                gas,
                catalog,
                dt,
                controls.solve,
                heat_controls,
                check,
            )
            .map(|candidate| (candidate, ()))
        },
        material_step_limit,
        cancelled,
    )
}

/// Common conservative controller for original material/source interval owners.
///
/// The callback must evaluate exactly one original chamber interval from the
/// supplied physical states, without retaining mutable physical history or
/// internally advancing extra time. Each invocation may resample its own film
/// sources; only the last ACCEPTED fine payload is retained. Fixed physical
/// identities, original calorics/catalog and external gas volume cannot change.
///
/// `retry_step_size` may identify only an original phase `StepSizeLimit` and
/// return that refusal's positive maximum duration. Source support, caloric,
/// conservation, cancellation and other errors remain terminal. The initial
/// zero-duration callback precedes caloric queries, preserving source owners'
/// identity-before-query admission. It also supplies the zero-window payload.
/// Empty rosters retain that admission with no fictitious physical intervals.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_droplets_with_gas_adaptive_using<'surface, 'catalog, T, E>(
    seeds: &[GasDropletSeed<'surface>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: AdaptiveGasDropletControls,
    mut step: impl FnMut(
        &[GasDropletSeed<'surface>],
        &ClosedGasParcel,
        f64,
        &mut dyn FnMut() -> bool,
    ) -> Result<(GasDropletChamberCandidate<'catalog>, T), E>,
    retry_step_size: impl Fn(&E) -> Option<f64>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AdaptiveGasDropletCandidate<'catalog, T>, AdaptiveGasDropletError<E>> {
    poll(&mut cancelled)?;
    controls.validate(duration_s, seeds.len())?;
    let mut current = seeds.to_vec();
    current.sort_unstable_by_key(|seed| seed.id);
    if current.windows(2).any(|pair| pair[0].id == pair[1].id) {
        return Err(AdaptiveGasDropletError::Invalid(
            "adaptive spherical droplets require unique stable IDs",
        ));
    }
    let mut attempted_steps = 0;
    let (mut endpoint, mut payload) = run_step(
        &current,
        gas,
        catalog,
        0.0,
        0.0,
        &mut attempted_steps,
        &mut step,
        &mut cancelled,
    )?;
    let mut account = account::Account::new(&current, gas, catalog, &mut cancelled)?;
    let mut half_seeds = current.clone();
    let mut time = 0.0;
    let mut span = controls.initial_dt_s.min(controls.maximum_dt_s);
    let mut comparisons = 0;
    let mut accepted_comparisons = 0;
    let mut accepted_intervals = 0usize;
    let mut accuracy_rejections = 0;
    let mut physical_rejections = 0;
    let mut reductions = 0;
    let mut maximum_accepted_normalized_error = 0.0_f64;
    let mut last_error = None;
    let order = ConductionTimeOrder::First;
    while time < duration_s && !current.is_empty() {
        poll(&mut cancelled)?;
        if comparisons >= controls.maximum_comparisons {
            return Err(AdaptiveGasDropletError::ComparisonBudget);
        }
        let end = if span >= duration_s - time {
            duration_s
        } else {
            finite(time + span)?
        };
        let middle = f64::midpoint(time, end);
        if !(time < middle && middle < end && end <= duration_s) {
            return Err(AdaptiveGasDropletError::Invalid(
                "no representable adaptive spherical midpoint",
            ));
        }
        let dt = finite(end - time)?;
        comparisons = increment(comparisons)?;
        let comparison = (|| -> Result<_, AdaptiveGasDropletError<E>> {
            let (coarse, _) = run_step(
                &current,
                endpoint.gas(),
                catalog,
                time,
                dt,
                &mut attempted_steps,
                &mut step,
                &mut cancelled,
            )?;
            let (first, _) = run_step(
                &current,
                endpoint.gas(),
                catalog,
                time,
                middle - time,
                &mut attempted_steps,
                &mut step,
                &mut cancelled,
            )?;
            for (seed, updated) in half_seeds.iter_mut().zip(first.droplets()) {
                seed.droplet = updated.droplet;
            }
            let (fine, last_payload) = run_step(
                &half_seeds,
                first.gas(),
                catalog,
                middle,
                end - middle,
                &mut attempted_steps,
                &mut step,
                &mut cancelled,
            )?;
            let correction = order
                .fine_error_correction(middle - time, end - middle)
                .ok_or(AdaptiveGasDropletError::Invalid(
                    "adaptive spherical error correction is unrepresentable",
                ))?;
            let error = time_error(
                &endpoint,
                &coarse,
                &fine,
                controls,
                correction,
                &mut cancelled,
            )?;
            Ok((first, fine, last_payload, error))
        })();
        let (first, fine, last_payload, error) = match comparison {
            Ok(result) => result,
            Err(cause) => {
                let factor = match &cause {
                    AdaptiveGasDropletError::Step {
                        duration_s, error, ..
                    } => retry_step_size(error)
                        .filter(|hint| hint.is_finite() && *hint > 0.0)
                        .map(|hint| (hint / *duration_s).min(0.5)),
                    _ => None,
                };
                let Some(factor) = factor else {
                    return Err(cause);
                };
                physical_rejections = increment(physical_rejections)?;
                if reductions >= controls.maximum_reductions || dt <= controls.minimum_dt_s {
                    return Err(cause);
                }
                let Some(smaller) = smaller_span(time, end, dt * factor, controls.minimum_dt_s)
                else {
                    return Err(cause);
                };
                span = smaller;
                reductions = increment(reductions)?;
                continue;
            }
        };
        if error.normalized > 1.0 {
            accuracy_rejections = increment(accuracy_rejections)?;
            let cause = AdaptiveGasDropletError::TimeAccuracy {
                time_s: time,
                span_s: dt,
                error,
            };
            if reductions >= controls.maximum_reductions || dt <= controls.minimum_dt_s {
                return Err(cause);
            }
            let factor = order.rejected_reduction_factor(error.normalized).ok_or(
                AdaptiveGasDropletError::Invalid("invalid adaptive spherical time error"),
            )?;
            let Some(smaller) = smaller_span(time, end, dt * factor, controls.minimum_dt_s) else {
                return Err(cause);
            };
            span = smaller;
            reductions = increment(reductions)?;
            continue;
        }
        account.add(&first, &mut cancelled)?;
        account.add(&fine, &mut cancelled)?;
        account.gate(controls)?;
        accepted_comparisons = increment(accepted_comparisons)?;
        accepted_intervals = accepted_intervals
            .checked_add(2)
            .ok_or(AdaptiveGasDropletError::ComparisonBudget)?;
        maximum_accepted_normalized_error = maximum_accepted_normalized_error.max(error.normalized);
        last_error = Some(error);
        let growth = order.accepted_growth_factor(error.normalized).ok_or(
            AdaptiveGasDropletError::Invalid("invalid accepted spherical time error"),
        )?;
        let suggested = if growth >= 1.0 {
            dt.min(controls.maximum_dt_s / growth) * growth
        } else {
            dt * growth
        };
        span = suggested.clamp(controls.minimum_dt_s, controls.maximum_dt_s);
        for (seed, updated) in current.iter_mut().zip(fine.droplets()) {
            seed.droplet = updated.droplet;
        }
        endpoint = fine;
        payload = last_payload;
        time = end;
        reductions = 0;
    }
    let exchange = account.finish(&endpoint, controls, &mut cancelled)?;
    poll(&mut cancelled)?;
    Ok(AdaptiveGasDropletCandidate {
        endpoint,
        payload,
        report: AdaptiveGasDropletReport {
            duration_s,
            exchange,
            attempted_steps,
            comparisons,
            accepted_comparisons,
            accepted_intervals,
            accuracy_rejections,
            physical_rejections,
            maximum_accepted_normalized_error,
            last_error,
            suggested_dt_s: span,
        },
    })
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn run_step<'surface, 'catalog, T, E>(
    seeds: &[GasDropletSeed<'surface>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    time_s: f64,
    duration_s: f64,
    attempted_steps: &mut usize,
    step: &mut impl FnMut(
        &[GasDropletSeed<'surface>],
        &ClosedGasParcel,
        f64,
        &mut dyn FnMut() -> bool,
    ) -> Result<(GasDropletChamberCandidate<'catalog>, T), E>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(GasDropletChamberCandidate<'catalog>, T), AdaptiveGasDropletError<E>> {
    poll(cancelled)?;
    *attempted_steps = increment(*attempted_steps)?;
    let (next, payload) =
        step(seeds, gas, duration_s, cancelled).map_err(|error| AdaptiveGasDropletError::Step {
            time_s,
            duration_s,
            error,
        })?;
    if !std::ptr::eq(next.catalog(), catalog)
        || next.gas().volume_m3().to_bits() != gas.volume_m3().to_bits()
        || next.report().exchange.duration_s.to_bits() != duration_s.to_bits()
        || next.report().phase.duration_s.to_bits() != duration_s.to_bits()
        || next.droplets().len() != seeds.len()
    {
        return Err(AdaptiveGasDropletError::Invalid(
            "original adaptive droplet step changed its catalog, volume, interval or roster",
        ));
    }
    for (seed, updated) in seeds.iter().zip(next.droplets()) {
        poll(cancelled)?;
        let initial_geometry = SphericalDropletGeometry::from_mass(
            seed.droplet.mass_kg,
            updated.geometry[0].liquid_density_kg_m3(),
        )
        .map_err(GasEvaporationError::from)?;
        if seed.id != updated.id || initial_geometry != updated.geometry[0] {
            return Err(AdaptiveGasDropletError::Invalid(
                "original adaptive droplet step changed an input identity or mass",
            ));
        }
        if duration_s == 0.0 && seed.droplet != updated.droplet {
            return Err(AdaptiveGasDropletError::Invalid(
                "zero-duration spherical admission changed a liquid inventory",
            ));
        }
    }
    if duration_s == 0.0 && next.gas() != gas {
        return Err(AdaptiveGasDropletError::Invalid(
            "zero-duration spherical admission changed the original atmosphere",
        ));
    }
    poll(cancelled)?;
    Ok((next, payload))
}

fn material_step_limit(error: &GasEvaporationError) -> Option<f64> {
    match error {
        GasEvaporationError::StepSizeLimit {
            maximum_duration_s, ..
        } => Some(*maximum_duration_s),
        _ => None,
    }
}

fn smaller_span(start: f64, end: f64, wanted: f64, floor: f64) -> Option<f64> {
    let mut next = start + wanted.max(floor);
    if next - start < floor {
        next = next.next_up();
    }
    (start < next && next < end).then_some(next - start)
}

fn increment<E>(value: usize) -> Result<usize, AdaptiveGasDropletError<E>> {
    value
        .checked_add(1)
        .ok_or(AdaptiveGasDropletError::ComparisonBudget)
}

fn finite(value: f64) -> Result<f64, GasEvaporationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(GasEvaporationError::Invalid(
            "adaptive spherical arithmetic is unrepresentable",
        ))
    }
}

fn component_error(
    old: f64,
    coarse: f64,
    fine: f64,
    absolute: f64,
    relative: f64,
    correction: f64,
) -> Result<(f64, f64), GasEvaporationError> {
    let difference = finite(fine - coarse)?.abs();
    let change = finite(coarse - old)?.abs().max(finite(fine - old)?.abs());
    let allowance = finite(absolute + finite(relative * change)?)?;
    let ratio = (difference / allowance) * correction;
    // An infinite ratio remains a rejection, not a fabricated finite estimate.
    if ratio.is_nan() || ratio < 0.0 || allowance <= 0.0 {
        return Err(GasEvaporationError::Invalid(
            "invalid adaptive spherical error scale",
        ));
    }
    Ok((difference, ratio))
}

fn time_error(
    old: &GasDropletChamberCandidate<'_>,
    coarse: &GasDropletChamberCandidate<'_>,
    fine: &GasDropletChamberCandidate<'_>,
    controls: AdaptiveGasDropletControls,
    correction: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<GasDropletTimeError, GasEvaporationError> {
    let mut error = GasDropletTimeError::default();
    for ((old, coarse), fine) in old
        .droplets()
        .iter()
        .zip(coarse.droplets())
        .zip(fine.droplets())
    {
        poll(cancelled)?;
        let (difference, ratio) = component_error(
            old.droplet.specific_enthalpy_j_kg,
            coarse.droplet.specific_enthalpy_j_kg,
            fine.droplet.specific_enthalpy_j_kg,
            controls.enthalpy_absolute_tolerance_j_kg,
            controls.relative_tolerance,
            correction,
        )?;
        error.maximum_enthalpy_difference_j_kg =
            error.maximum_enthalpy_difference_j_kg.max(difference);
        error.normalized = error.normalized.max(ratio);
        let (difference, ratio) = component_error(
            old.exchange.surface_temperature_k[1],
            coarse.exchange.surface_temperature_k[1],
            fine.exchange.surface_temperature_k[1],
            controls.temperature_absolute_tolerance_k,
            controls.relative_tolerance,
            correction,
        )?;
        error.maximum_temperature_difference_k =
            error.maximum_temperature_difference_k.max(difference);
        error.normalized = error.normalized.max(ratio);
        let (difference, ratio) = component_error(
            old.droplet.mass_kg,
            coarse.droplet.mass_kg,
            fine.droplet.mass_kg,
            controls.mass_absolute_tolerance_kg,
            controls.relative_tolerance,
            correction,
        )?;
        error.maximum_liquid_mass_difference_kg =
            error.maximum_liquid_mass_difference_kg.max(difference);
        error.normalized = error.normalized.max(ratio);
        let (difference, ratio) = component_error(
            old.geometry[1].radius_m(),
            coarse.geometry[1].radius_m(),
            fine.geometry[1].radius_m(),
            controls.radius_absolute_tolerance_m,
            controls.relative_tolerance,
            correction,
        )?;
        error.maximum_radius_difference_m = error.maximum_radius_difference_m.max(difference);
        error.normalized = error.normalized.max(ratio);
    }
    let (difference, ratio) = component_error(
        old.gas().state().state().temperature,
        coarse.gas().state().state().temperature,
        fine.gas().state().state().temperature,
        controls.temperature_absolute_tolerance_k,
        controls.relative_tolerance,
        correction,
    )?;
    error.gas_temperature_difference_k = difference;
    error.normalized = error.normalized.max(ratio);
    let initial_inventory = old.gas().species_inventory_kg()?;
    let coarse_inventory = coarse.gas().species_inventory_kg()?;
    let fine_inventory = fine.gas().species_inventory_kg()?;
    for component in old.catalog().thermodynamics().components() {
        poll(cancelled)?;
        let species = component.species();
        let (difference, ratio) = component_error(
            initial_inventory.get(species).copied().unwrap_or(0.0),
            coarse_inventory.get(species).copied().unwrap_or(0.0),
            fine_inventory.get(species).copied().unwrap_or(0.0),
            controls.mass_absolute_tolerance_kg,
            controls.relative_tolerance,
            correction,
        )?;
        error.maximum_species_mass_difference_kg =
            error.maximum_species_mass_difference_kg.max(difference);
        error.normalized = error.normalized.max(ratio);
    }
    poll(cancelled)?;
    Ok(error)
}

#[cfg(test)]
mod tests;
