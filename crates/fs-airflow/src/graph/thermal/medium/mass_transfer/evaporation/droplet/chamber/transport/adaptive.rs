//! Adaptive evaporation in a complete original finite-gas network.
//!
//! Each coarse or fine trial runs the original transport-then-droplet owner
//! from its own immutable gas map and liquid inventories. Every gas node's
//! temperature, total mass and species masses enter the local time estimate,
//! including nodes without droplets and species newly arriving at a receiver.
//! Only the two accepted fine intervals enter the physical window account;
//! only their final candidate retains endpoint source and transport receipts.
//!
//! Transport samples before droplet films, so source admission always sees the
//! actual transported atmosphere. No initial zero-time source query or hidden
//! mutable gas history changes this sequence. Positive duration is required by
//! the original network owner, and an empty droplet roster still evolves gas.
//!
//! Typed phase step limits, original species withdrawal and pressure/permeation
//! overshoot may refine physical time. Source, caloric, conservation and
//! cancellation refusals retain their original errors and remain terminal.
//! Fixed volumes, positive liquid seeds, prescribed relative speeds and heat
//! coefficients retain their original meaning. The estimator controls local
//! first-order splitting error, not a certified global trajectory bound.

use std::collections::BTreeMap;

use super::{
    ClosedGasParcel, DropletGasTransportError, GasDropletSurfaceStep,
    GasDropletTransportCandidate, GasDropletTransportSolve, GasEnergySolve,
    GasLiquidChamberHeatSolve, GasTransportError, GasTransportPaths, GasTransportSolve,
    SourceDropletFilmStep, SourceGasDropletHeatSurface, SourceGasDropletSurface,
    SphericalDropletFilmError, ThermochemicalGasModel,
    advance_droplets_with_gas_transport as advance_interval,
    advance_droplets_with_gas_transport_and_heat as advance_interval_with_heat,
};
use super::super::GasFilmHeatSolve;
use fs_material::evaporation::gas::{
    GasEvaporationError,
    droplet::chamber::adaptive::{
        AdaptiveGasDropletError, GasDropletSeed, GasDropletTimeError,
        network::{self, GasDropletNetworkAudit, GasDropletNetworkWindow},
    },
};
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::network::snapshot::SnapshotNozzleError,
    species::network::{
        SpeciesNetworkError,
        transport::{membrane::MembraneNetworkError, porous::PorousNetworkError},
    },
};
use fs_thermal::enthalpy::conduction::implicit::adaptive::ConductionTimeOrder;

pub use fs_material::evaporation::gas::droplet::chamber::adaptive::{
    AdaptiveGasDropletControls, network::GasDropletNetworkLimits,
};

/// One original phase/transport allowance and one complete-window time budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveGasDropletTransportControls {
    pub droplets: AdaptiveGasDropletControls,
    pub transport: GasTransportSolve,
    pub network: GasDropletNetworkLimits,
}
impl AdaptiveGasDropletTransportControls {
    pub fn new(solve: GasDropletTransportSolve) -> Self {
        let mut droplets = AdaptiveGasDropletControls::new(solve.phase);
        droplets.mass_conservation_tolerance_kg = solve.mass_tolerance_kg;
        droplets.energy_conservation_tolerance_j = solve.energy_tolerance_j;
        Self {
            droplets,
            transport: solve.transport,
            network: GasDropletNetworkLimits::default(),
        }
    }

    fn interval(self) -> GasDropletTransportSolve {
        GasDropletTransportSolve {
            transport: self.transport,
            phase: self.droplets.solve,
            mass_tolerance_kg: self.droplets.mass_conservation_tolerance_kg,
            energy_tolerance_j: self.droplets.energy_conservation_tolerance_j,
        }
    }
}

/// Original coupled refusals retain the exact attempted time and duration.
pub type AdaptiveGasDropletTransportError = AdaptiveGasDropletError<DropletGasTransportError>;

/// Whole-window accepted physics plus bounded local time-control work.
#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveGasDropletTransportReport {
    pub duration_s: f64,
    pub exchange: GasDropletNetworkAudit,
    pub attempted_steps: usize,
    pub comparisons: usize,
    pub accepted_comparisons: usize,
    pub accepted_intervals: usize,
    pub accuracy_rejections: usize,
    pub physical_rejections: usize,
    pub maximum_accepted_normalized_error: f64,
    pub last_error: Option<GasDropletTimeError>,
    pub suggested_dt_s: f64,
}

/// Complete final gas map, droplets and final original source evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveGasDropletTransportCandidate<'catalog> {
    endpoint: GasDropletTransportCandidate<'catalog>,
    report: AdaptiveGasDropletTransportReport,
}
impl<'catalog> AdaptiveGasDropletTransportCandidate<'catalog> {
    /// Last accepted fine interval; its local receipts are not window totals.
    pub const fn endpoint(&self) -> &GasDropletTransportCandidate<'catalog> {
        &self.endpoint
    }
    pub const fn report(&self) -> &AdaptiveGasDropletTransportReport {
        &self.report
    }
    pub fn droplets(&self) -> &[GasDropletSurfaceStep] {
        self.endpoint.droplets()
    }
    pub const fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        self.endpoint.parcels()
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.endpoint.catalog()
    }
    pub fn films(&self) -> &[SourceDropletFilmStep] {
        self.endpoint.films()
    }
    pub fn into_parts(
        self,
    ) -> (
        GasDropletTransportCandidate<'catalog>,
        AdaptiveGasDropletTransportReport,
    ) {
        (self.endpoint, self.report)
    }
}

/// Evolve original gas transport and sourced phase exchange for positive time.
#[allow(clippy::too_many_arguments)]
pub fn advance_droplets_with_gas_transport_adaptive<'catalog>(
    surfaces: &[SourceGasDropletSurface<'_, '_>],
    gas_region: u64,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    paths: GasTransportPaths<'_>,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    duration_s: f64,
    controls: AdaptiveGasDropletTransportControls,
    cancelled: impl FnMut() -> bool,
) -> Result<AdaptiveGasDropletTransportCandidate<'catalog>, AdaptiveGasDropletTransportError> {
    advance(
        surfaces.iter().copied().map(|surface| (surface, None)),
        gas_region,
        parcels,
        catalog,
        paths,
        endpoints,
        duration_s,
        controls,
        None,
        cancelled,
    )
}

/// The same complete trajectory with original collective reciprocal heat.
#[allow(clippy::too_many_arguments)]
pub fn advance_droplets_with_gas_transport_and_heat_adaptive<'catalog>(
    surfaces: &[SourceGasDropletHeatSurface<'_, '_>],
    gas_region: u64,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    paths: GasTransportPaths<'_>,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    duration_s: f64,
    controls: AdaptiveGasDropletTransportControls,
    heat_controls: GasLiquidChamberHeatSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<AdaptiveGasDropletTransportCandidate<'catalog>, AdaptiveGasDropletTransportError> {
    advance(
        surfaces.iter().map(|row| {
            (
                row.surface,
                Some((row.heat_transfer_coefficient_w_m2_k, row.heat_controls)),
            )
        }),
        gas_region,
        parcels,
        catalog,
        paths,
        endpoints,
        duration_s,
        controls,
        Some(heat_controls),
        cancelled,
    )
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn advance<'surface, 'm: 'surface, 'catalog>(
    surfaces: impl ExactSizeIterator<
        Item = (
            SourceGasDropletSurface<'surface, 'm>,
            Option<(f64, GasFilmHeatSolve)>,
        ),
    >,
    gas_region: u64,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    paths: GasTransportPaths<'_>,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    duration_s: f64,
    controls: AdaptiveGasDropletTransportControls,
    heat_controls: Option<GasLiquidChamberHeatSolve>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AdaptiveGasDropletTransportCandidate<'catalog>, AdaptiveGasDropletTransportError> {
    poll(&mut cancelled)?;
    if !duration_s.is_finite() || duration_s <= 0.0 {
        return Err(AdaptiveGasDropletError::Invalid(
            "adaptive droplet gas transport requires positive finite duration",
        ));
    }
    let numerical = controls.droplets;
    numerical.validate(duration_s, surfaces.len())?;
    network::admit_work(parcels, controls.network, &mut cancelled)?;
    let mut ordered: Vec<_> = surfaces.collect();
    ordered.sort_unstable_by_key(|(row, _)| row.id);
    if ordered.windows(2).any(|rows| rows[0].0.id == rows[1].0.id) {
        return Err(AdaptiveGasDropletError::Invalid(
            "adaptive transported droplets require unique stable IDs",
        ));
    }
    let mut current: Vec<_> = ordered
        .iter()
        .map(|(row, _)| GasDropletSeed {
            id: row.id,
            law: row.law,
            binding: row.binding,
            droplet: row.droplet,
            temperature_from_enthalpy: row.temperature_from_enthalpy,
        })
        .collect();
    let mut account = GasDropletNetworkWindow::new(
        &current,
        gas_region,
        parcels,
        numerical,
        controls.network,
        &mut cancelled,
    )?;
    // Exactly one source-row buffer is allocated. Trials overwrite only the
    // physical droplet inventory; the original law and film selections persist.
    let mut phase_rows = Vec::new();
    let mut heat_rows = Vec::new();
    if heat_controls.is_some() {
        heat_rows = ordered
            .into_iter()
            .map(|(surface, heat)| {
                let (coefficient, heat_controls) = heat.ok_or(
                    AdaptiveGasDropletError::Invalid("adaptive transported heat row is missing"),
                )?;
                Ok(SourceGasDropletHeatSurface {
                    surface,
                    heat_transfer_coefficient_w_m2_k: coefficient,
                    heat_controls,
                })
            })
            .collect::<Result<Vec<_>, AdaptiveGasDropletTransportError>>()?;
    } else {
        phase_rows = ordered.into_iter().map(|(surface, _)| surface).collect();
    }
    let mut evaluate = |seeds: &[GasDropletSeed<'surface>],
                        gases: &BTreeMap<u64, ClosedGasParcel>,
                        duration: f64,
                        check: &mut dyn FnMut() -> bool| {
        if let Some(heat_controls) = heat_controls {
            for (row, seed) in heat_rows.iter_mut().zip(seeds) {
                if check() {
                    return Err(DropletGasTransportError::Cancelled);
                }
                row.surface.droplet = seed.droplet;
            }
            advance_interval_with_heat(
                &heat_rows,
                gas_region,
                gases,
                catalog,
                paths,
                endpoints,
                duration,
                controls.interval(),
                heat_controls,
                check,
            )
        } else {
            for (row, seed) in phase_rows.iter_mut().zip(seeds) {
                if check() {
                    return Err(DropletGasTransportError::Cancelled);
                }
                row.droplet = seed.droplet;
            }
            advance_interval(
                &phase_rows,
                gas_region,
                gases,
                catalog,
                paths,
                endpoints,
                duration,
                controls.interval(),
                check,
            )
        }
    };
    let mut half_seeds = current.clone();
    let mut endpoint: Option<GasDropletTransportCandidate<'catalog>> = None;
    let mut time = 0.0;
    let mut span = numerical.initial_dt_s.min(numerical.maximum_dt_s);
    let mut attempted_steps = 0;
    let mut comparisons = 0;
    let mut accepted_comparisons = 0;
    let mut accepted_intervals = 0usize;
    let mut accuracy_rejections = 0;
    let mut physical_rejections = 0;
    let mut reductions = 0;
    let mut maximum_accepted_normalized_error = 0.0_f64;
    let mut last_error = None;
    let order = ConductionTimeOrder::First;
    while time < duration_s {
        poll(&mut cancelled)?;
        if comparisons >= numerical.maximum_comparisons {
            return Err(AdaptiveGasDropletError::ComparisonBudget);
        }
        let end = trial_endpoint(time, span, duration_s)?;
        let middle = f64::midpoint(time, end);
        if !(time < middle && middle < end && end <= duration_s) {
            return Err(AdaptiveGasDropletError::Invalid(
                "no representable adaptive transported droplet midpoint",
            ));
        }
        let dt = finite(end - time)?;
        comparisons = increment(comparisons)?;
        let before = endpoint.as_ref().map_or(parcels, |last| last.parcels());
        let comparison = (|| -> Result<_, AdaptiveGasDropletTransportError> {
            let coarse = run_step(
                &current, before, catalog, gas_region, time, dt,
                controls.network, &mut attempted_steps, &mut evaluate, &mut cancelled,
            )?;
            let first = run_step(
                &current, before, catalog, gas_region, time, middle - time,
                controls.network, &mut attempted_steps, &mut evaluate, &mut cancelled,
            )?;
            for (seed, row) in half_seeds.iter_mut().zip(first.droplets()) {
                seed.droplet = row.droplet;
            }
            let fine = run_step(
                &half_seeds, first.parcels(), catalog, gas_region, middle, end - middle,
                controls.network, &mut attempted_steps, &mut evaluate, &mut cancelled,
            )?;
            let correction = order
                .fine_error_correction(middle - time, end - middle)
                .ok_or(AdaptiveGasDropletError::Invalid(
                    "adaptive transported droplet time correction is unrepresentable",
                ))?;
            let error = network::time_error(
                &current,
                before,
                (coarse.droplets(), coarse.parcels()),
                (fine.droplets(), fine.parcels()),
                numerical,
                correction,
                controls.network,
                &mut cancelled,
            )?;
            Ok((first, fine, error))
        })();
        let (first, fine, error) = match comparison {
            Ok(result) => result,
            Err(cause) => {
                let factor = match &cause {
                    AdaptiveGasDropletError::Step { duration_s, error, .. } => {
                        physical_reduction(error, *duration_s)
                    }
                    _ => None,
                };
                let Some(factor) = factor else {
                    return Err(cause);
                };
                physical_rejections = increment(physical_rejections)?;
                if reductions >= numerical.maximum_reductions || dt <= numerical.minimum_dt_s {
                    return Err(cause);
                }
                let Some(smaller) =
                    smaller_span(time, end, dt * factor, numerical.minimum_dt_s)
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
            if reductions >= numerical.maximum_reductions || dt <= numerical.minimum_dt_s {
                return Err(cause);
            }
            let factor = order.rejected_reduction_factor(error.normalized).ok_or(
                AdaptiveGasDropletError::Invalid("invalid transported droplet time error"),
            )?;
            let Some(smaller) =
                smaller_span(time, end, dt * factor, numerical.minimum_dt_s)
            else {
                return Err(cause);
            };
            span = smaller;
            reductions = increment(reductions)?;
            continue;
        }
        account = account.record(
            &first.report().transport,
            &first.report().coupled,
            first.droplets(),
            &mut cancelled,
        )?;
        account = account.record(
            &fine.report().transport,
            &fine.report().coupled,
            fine.droplets(),
            &mut cancelled,
        )?;
        accepted_comparisons = increment(accepted_comparisons)?;
        accepted_intervals = accepted_intervals
            .checked_add(2)
            .ok_or(AdaptiveGasDropletError::ComparisonBudget)?;
        maximum_accepted_normalized_error = maximum_accepted_normalized_error.max(error.normalized);
        last_error = Some(error);
        let growth = order.accepted_growth_factor(error.normalized).ok_or(
            AdaptiveGasDropletError::Invalid("invalid accepted transported droplet time error"),
        )?;
        let suggested = if growth >= 1.0 {
            dt.min(numerical.maximum_dt_s / growth) * growth
        } else {
            dt * growth
        };
        span = suggested.clamp(numerical.minimum_dt_s, numerical.maximum_dt_s);
        for (seed, row) in current.iter_mut().zip(fine.droplets()) {
            seed.droplet = row.droplet;
        }
        endpoint = Some(fine);
        time = end;
        reductions = 0;
    }
    let endpoint = endpoint.ok_or(AdaptiveGasDropletError::Invalid(
        "positive adaptive transport produced no accepted physical interval",
    ))?;
    let exchange = account.finish(endpoint.parcels(), endpoint.droplets(), &mut cancelled)?;
    poll(&mut cancelled)?;
    Ok(AdaptiveGasDropletTransportCandidate {
        endpoint,
        report: AdaptiveGasDropletTransportReport {
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
fn run_step<'surface, 'catalog>(
    seeds: &[GasDropletSeed<'surface>],
    before: &BTreeMap<u64, ClosedGasParcel>,
    catalog: &'catalog ThermochemicalGasModel,
    gas_region: u64,
    time_s: f64,
    duration_s: f64,
    limits: GasDropletNetworkLimits,
    attempted_steps: &mut usize,
    evaluate: &mut impl FnMut(
        &[GasDropletSeed<'surface>],
        &BTreeMap<u64, ClosedGasParcel>,
        f64,
        &mut dyn FnMut() -> bool,
    ) -> Result<GasDropletTransportCandidate<'catalog>, DropletGasTransportError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<GasDropletTransportCandidate<'catalog>, AdaptiveGasDropletTransportError> {
    poll(cancelled)?;
    *attempted_steps = increment(*attempted_steps)?;
    let next = evaluate(seeds, before, duration_s, cancelled).map_err(|error| {
        AdaptiveGasDropletError::Step {
            time_s,
            duration_s,
            error,
        }
    })?;
    let report = next.report();
    if !std::ptr::eq(next.catalog(), catalog)
        || report.gas_region != gas_region
        || report.duration_s.to_bits() != duration_s.to_bits()
        || report.transport.duration_s.to_bits() != duration_s.to_bits()
        || report.coupled.exchange.duration_s.to_bits() != duration_s.to_bits()
        || !before.keys().eq(next.parcels().keys())
        || seeds.len() != next.droplets().len()
        || seeds.iter().zip(next.droplets()).any(|(seed, row)| seed.id != row.id)
    {
        return Err(AdaptiveGasDropletError::Invalid(
            "original transported droplet trial changed its roster, catalog or interval",
        ));
    }
    network::admit_work(next.parcels(), limits, cancelled)?;
    for (&id, old) in before {
        poll(cancelled)?;
        if old.volume_m3().to_bits() != next.parcels()[&id].volume_m3().to_bits() {
            return Err(AdaptiveGasDropletError::Invalid(
                "original transported droplet trial changed a fixed gas volume",
            ));
        }
    }
    poll(cancelled)?;
    Ok(next)
}

fn physical_reduction(error: &DropletGasTransportError, duration_s: f64) -> Option<f64> {
    match error {
        DropletGasTransportError::Droplets(SphericalDropletFilmError::Phase(
            GasEvaporationError::StepSizeLimit { maximum_duration_s, .. },
        )) if maximum_duration_s.is_finite() && *maximum_duration_s > 0.0 => {
            Some((maximum_duration_s / duration_s).min(0.5))
        }
        DropletGasTransportError::Transport(
            GasTransportError::Exchange(SpeciesNetworkError::Withdrawal { .. })
            | GasTransportError::Nozzle(SnapshotNozzleError::EquilibriumCrossed { .. })
            | GasTransportError::Porous(PorousNetworkError::EquilibriumCrossed { .. })
            | GasTransportError::Membrane(MembraneNetworkError::EquilibriumCrossed { .. }),
        ) => Some(0.5),
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

/// Endpoint of the next step-doubling trial: `elapsed_s + step_s`, or `target_s` once the step
/// reaches it. Each trial is compared with two half steps, so the interval it leaves behind must
/// itself be divisible. When the midpoint of the leftover `[endpoint, target_s]` is not strictly
/// inside it (a remainder of one ulp, left by rounded sums of steps), no later trial could take
/// it, and this trial runs on to `target_s` instead. Every other endpoint is the one the step
/// alone gives. This is the twin of `crate::miscibility::adaptive::trial_endpoint` in fs-flux
/// (e38444e640), which fixed the same loop there; one shared copy is a later job (bead .7).
fn trial_endpoint(
    elapsed_s: f64,
    step_s: f64,
    target_s: f64,
) -> Result<f64, AdaptiveGasDropletTransportError> {
    let endpoint = if step_s >= target_s - elapsed_s {
        target_s
    } else {
        finite(elapsed_s + step_s)?
    };
    let leftover_midpoint = endpoint + 0.5 * (target_s - endpoint);
    let divisible = endpoint < leftover_midpoint && leftover_midpoint < target_s;
    Ok(if endpoint < target_s && !divisible {
        target_s
    } else {
        endpoint
    })
}

#[cfg(test)]
mod trial_endpoint_tests {
    use super::trial_endpoint;

    #[test]
    fn a_trial_that_would_strand_an_indivisible_remainder_runs_to_the_target() {
        // The droplet fixture's loop reached 9.999999999999998e-3 and took a 1e-2 step: the rounded
        // sum is one ulp below the 2e-2 target, and the leftover has no representable midpoint.
        let (elapsed, step, target) = (9.999999999999998e-3_f64, 1e-2_f64, 2e-2_f64);
        assert_eq!(elapsed + step, 1.9999999999999997e-2);
        assert_eq!(1.9999999999999997e-2_f64.next_up(), target);
        assert_eq!(trial_endpoint(elapsed, step, target).unwrap(), target);
        // From one ulp below the target any step runs on to it.
        assert_eq!(trial_endpoint(1.9999999999999997e-2, 1e-2, target).unwrap(), target);
    }

    #[test]
    fn an_ordinary_trial_endpoint_is_the_one_the_step_alone_gives() {
        assert_eq!(trial_endpoint(0.0, 5e-3, 2e-2).unwrap(), 5e-3);
        assert_eq!(trial_endpoint(5e-3, 1e-2, 2e-2).unwrap(), 5e-3 + 1e-2);
        // A step that reaches the target ends on it.
        assert_eq!(trial_endpoint(1.5e-2, 2e-2, 2e-2).unwrap(), 2e-2);
        // A leftover of two ulps still has a midpoint strictly inside it: the endpoint is kept.
        let target = 2e-2_f64;
        let endpoint = target.next_down().next_down();
        assert_eq!(trial_endpoint(0.0, endpoint, target).unwrap(), endpoint);
        // A step that is not a number is refused as before.
        assert!(trial_endpoint(0.0, f64::NAN, target).is_err());
    }
}

fn increment(value: usize) -> Result<usize, AdaptiveGasDropletTransportError> {
    value.checked_add(1).ok_or(AdaptiveGasDropletError::ComparisonBudget)
}

fn finite(value: f64) -> Result<f64, AdaptiveGasDropletTransportError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(AdaptiveGasDropletError::Invalid(
            "adaptive transported droplet arithmetic is unrepresentable",
        ))
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AdaptiveGasDropletTransportError> {
    if cancelled() {
        Err(AdaptiveGasDropletError::Material(GasEvaporationError::Cancelled))
    } else {
        Ok(())
    }
}
