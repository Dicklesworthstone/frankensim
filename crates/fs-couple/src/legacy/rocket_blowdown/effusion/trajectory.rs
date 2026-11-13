//! Bounded adaptive time advancement of the original finite-gas recoil owner.
//!
//! Each trial compares one complete step with two half steps, then accepts only
//! the two half steps. Their underlying gas update remains first order. The
//! difference estimates LOCAL truncation error; it is not a certified or global
//! error bound. Original source, inventory and physical conservation gates stay
//! active in every trial, with no acceptance of a failed coarse physical solve.

/// Adaptive supply and discharge using the same bounded controller.
pub mod exchange;

use super::{
    EffusionOutlet, EffusionRocketEnvironment, EffusionRocketError, EffusionRocketSolve,
    EffusionRocketStep, GasRocketState, Sum,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::{
    SpeciesNetworkError, transport::GasTransportError,
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_mbd::Vec3;
use fs_qty::MolarMass;
use fs_thermochem::Nasa9StandardStateModelV1;
use fs_thermochem::mixture::MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1;
use std::collections::BTreeMap;

/// Temporal accuracy and independent, finite work bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffusionAdaptiveSolve {
    pub initial_step_s: f64,
    /// Minimum trial interval, except a shorter final remainder. Each accepted
    /// interval advances through two half steps; a failed shorter tail refuses.
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    pub maximum_attempts: usize,
    /// At most twice this many native half-step reports are retained.
    pub maximum_accepted_intervals: usize,
    pub relative_tolerance: f64,
    pub temperature_tolerance_k: f64,
    pub pressure_tolerance_pa: f64,
    pub species_mass_tolerance_kg: f64,
    /// Relative scaling removes both the coordinate origin and initial drift.
    pub position_tolerance_m: f64,
    /// Relative scaling uses the velocity change, not the common inertial boost.
    pub velocity_tolerance_m_s: f64,
}

/// Only accepted native reports contribute to these compensated totals.
/// Physical `EffusionRocketSolve` budgets apply to each accepted half step;
/// these totals report the whole interval without claiming a global error budget.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionTrajectoryStep {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub rejected_intervals: usize,
    pub accepted_intervals: usize,
    pub accepted_steps: Vec<EffusionRocketStep>,
    pub maximum_accepted_error: f64,
    pub total_escaped_mass_kg: f64,
    pub total_escaped_momentum_world_kg_m_s: Vec3,
    pub total_escaped_energy_world_j: f64,
    pub total_external_impulse_n_s: Vec3,
    pub total_external_work_j: f64,
    pub total_mass_residual_kg: f64,
    pub total_momentum_residual_kg_m_s: Vec3,
    pub total_energy_residual_j: f64,
    pub total_gravity_impulse_uncertainty_n_s: Vec3,
    pub total_gravity_work_uncertainty_j: f64,
}

/// A failed trajectory exposes no partially advanced gas or vessel.
#[derive(Clone, Debug, PartialEq)]
pub enum EffusionTrajectoryError {
    Invalid(&'static str),
    Unrepresentable(&'static str),
    Step {
        elapsed_s: f64,
        duration_s: f64,
        error: EffusionRocketError,
    },
    MinimumStep {
        duration_s: f64,
        normalized_error: f64,
    },
    WorkBudget {
        attempted_intervals: usize,
        accepted_intervals: usize,
    },
    Accounting(EffusionRocketError),
    Cancelled,
}
impl core::fmt::Display for EffusionTrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive effusive trajectory: {self:?}")
    }
}
impl std::error::Error for EffusionTrajectoryError {}

#[derive(Clone, Debug, PartialEq)]
pub struct EffusionTrajectoryCandidate {
    endpoint: GasRocketState,
    step: EffusionTrajectoryStep,
}
impl EffusionTrajectoryCandidate {
    #[must_use]
    pub const fn endpoint(&self) -> &GasRocketState {
        &self.endpoint
    }
    #[must_use]
    pub const fn step(&self) -> &EffusionTrajectoryStep {
        &self.step
    }
    #[must_use]
    pub fn into_parts(self) -> (GasRocketState, EffusionTrajectoryStep) {
        (self.endpoint, self.step)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), EffusionTrajectoryError> {
    if cancelled() {
        Err(EffusionTrajectoryError::Cancelled)
    } else {
        Ok(())
    }
}
pub(super) fn ratio(
    coarse: f64,
    fine: f64,
    scale: f64,
    absolute: f64,
    relative: f64,
) -> Result<f64, EffusionTrajectoryError> {
    let denominator = relative.mul_add(scale, absolute);
    let value = (fine - coarse).abs() / denominator;
    if !denominator.is_finite() || denominator <= 0.0 || !value.is_finite() {
        return Err(EffusionTrajectoryError::Unrepresentable(
            "local error or dimensional tolerance scale",
        ));
    }
    Ok(value)
}
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

// One original-source axis for the entire requested interval. Inflow may add a
// previously absent species, so its zero original mass must still participate
// in every local error comparison. Sources remain borrowed, without rebuilding
// a composition or copying material cards for each trial.
pub(super) struct SpeciesAxis<'a> {
    pub(super) sources: Vec<&'a Nasa9StandardStateModelV1>,
}
impl<'a> SpeciesAxis<'a> {
    pub(super) fn new(
        gases: impl IntoIterator<Item = &'a ClosedGasParcel>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, EffusionTrajectoryError> {
        let mut sources = BTreeMap::new();
        for gas in gases {
            poll(cancelled)?;
            for source in gas.model().thermodynamics().components() {
                poll(cancelled)?;
                if sources
                    .insert(source.species(), source)
                    .is_some_and(|old| old != source)
                {
                    return Err(EffusionTrajectoryError::Invalid(
                        "adaptive species union substitutes an original caloric source",
                    ));
                }
                if sources.len() > MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1 {
                    return Err(EffusionTrajectoryError::Invalid(
                        "adaptive species union exceeds the original gas work bound",
                    ));
                }
            }
        }
        Ok(Self {
            sources: sources.into_values().collect(),
        })
    }

    pub(super) fn masses(
        &self,
        state: &GasRocketState,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1], EffusionTrajectoryError> {
        poll(cancelled)?;
        let mixture = state.reservoir().model().thermodynamics();
        if mixture.components().len() > MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1 {
            return Err(EffusionTrajectoryError::Invalid(
                "adaptive endpoint species bound",
            ));
        }
        let mut molar = [MolarMass::new(0.0); MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1];
        let mut indices = [0; MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1];
        for (i, source) in mixture.components().iter().enumerate() {
            poll(cancelled)?;
            let index = self
                .sources
                .binary_search_by(|s| s.species().cmp(source.species()))
                .map_err(|_| {
                    EffusionTrajectoryError::Invalid(
                        "adaptive endpoint introduced an unqualified species",
                    )
                })?;
            if self.sources[index] != source {
                return Err(EffusionTrajectoryError::Invalid(
                    "adaptive endpoint changed its original caloric source",
                ));
            }
            indices[i] = index;
            molar[i] = self.sources[index].molar_mass();
        }
        let fractions = mixture
            .mole_composition()
            .to_mass_fractions(&molar[..mixture.components().len()])
            .map_err(|_| {
                EffusionTrajectoryError::Unrepresentable("original species mass fractions")
            })?;
        let mut masses = [0.0; MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1];
        for (i, fraction) in fractions.fractions().iter().enumerate() {
            masses[indices[i]] = state.reservoir().mass_kg() * fraction;
        }
        Ok(masses)
    }
}

pub(super) fn normalized_error(
    before: &GasRocketState,
    coarse: &GasRocketState,
    fine: &GasRocketState,
    species: &SpeciesAxis<'_>,
    h: f64,
    controls: EffusionAdaptiveSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, EffusionTrajectoryError> {
    let mut error = 0.0_f64;
    let [a, b, c] = [before, coarse, fine].map(|s| s.reservoir().state().state());
    for (old, full, half, tolerance) in [
        (
            a.temperature,
            b.temperature,
            c.temperature,
            controls.temperature_tolerance_k,
        ),
        (
            a.pressure,
            b.pressure,
            c.pressure,
            controls.pressure_tolerance_pa,
        ),
    ] {
        error = error.max(ratio(
            full,
            half,
            old.abs().max(full.abs()).max(half.abs()),
            tolerance,
            controls.relative_tolerance,
        )?);
    }
    let masses = [
        species.masses(before, cancelled)?,
        species.masses(coarse, cancelled)?,
        species.masses(fine, cancelled)?,
    ];
    for i in 0..species.sources.len() {
        poll(cancelled)?;
        let old = masses[0][i];
        let full = masses[1][i];
        let half = masses[2][i];
        error = error.max(ratio(
            full,
            half,
            old.abs().max(full.abs()).max(half.abs()),
            controls.species_mass_tolerance_kg,
            controls.relative_tolerance,
        )?);
    }
    let v0 = before.vehicle().velocity_m_s;
    let drift = v0.scale(h);
    for (full, half, tolerance) in [
        (
            coarse.vehicle().velocity_m_s.sub(v0),
            fine.vehicle().velocity_m_s.sub(v0),
            controls.velocity_tolerance_m_s,
        ),
        (
            coarse
                .vehicle()
                .position_m
                .sub(before.vehicle().position_m)
                .sub(drift),
            fine.vehicle()
                .position_m
                .sub(before.vehicle().position_m)
                .sub(drift),
            controls.position_tolerance_m,
        ),
    ] {
        for (a, b) in components(full).into_iter().zip(components(half)) {
            error = error.max(ratio(
                a,
                b,
                a.abs().max(b.abs()),
                tolerance,
                controls.relative_tolerance,
            )?);
        }
    }
    poll(cancelled)?;
    Ok(error)
}

fn withdrawal(error: &EffusionTrajectoryError) -> bool {
    matches!(
        error,
        EffusionTrajectoryError::Step {
            error: EffusionRocketError::Gas(GasTransportError::Exchange(
                SpeciesNetworkError::Withdrawal { .. }
            )),
            ..
        }
    )
}
fn values(s: &EffusionRocketStep) -> [f64; 18] {
    let p = s.exhaust_momentum_world_kg_m_s;
    let impulse = s.exhaust.external_impulse_n_s;
    let residual = s.momentum_residual_kg_m_s;
    let bound = s.gravity_impulse_uncertainty_n_s;
    [
        s.expelled_mass_kg,
        p.x,
        p.y,
        p.z,
        s.exhaust_energy_world_j,
        impulse.x,
        impulse.y,
        impulse.z,
        s.exhaust.external_work_j,
        s.gas.exchange.total_mass_residual_kg,
        residual.x,
        residual.y,
        residual.z,
        s.energy_residual_j,
        bound.x,
        bound.y,
        bound.z,
        s.gravity_work_uncertainty_j,
    ]
}

// One controller owns physical time, speculative trials, acceptance and all
// history bounds. Wrappers supply only the original physical step and receipt
// coordinates; rejected reports never reach an accumulator.
struct ControlledTrajectory<S, const N: usize> {
    endpoint: GasRocketState,
    duration_s: f64,
    attempted_intervals: usize,
    rejected_intervals: usize,
    accepted_intervals: usize,
    accepted_steps: Vec<S>,
    maximum_accepted_error: f64,
    totals: [f64; N],
}

pub(super) fn admit_adaptive(
    duration_s: f64,
    adaptive: EffusionAdaptiveSolve,
) -> Result<(), EffusionTrajectoryError> {
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || ![
            adaptive.initial_step_s,
            adaptive.minimum_step_s,
            adaptive.maximum_step_s,
            adaptive.temperature_tolerance_k,
            adaptive.pressure_tolerance_pa,
            adaptive.species_mass_tolerance_kg,
            adaptive.position_tolerance_m,
            adaptive.velocity_tolerance_m_s,
        ]
        .into_iter()
        .all(|v| v.is_finite() && v > 0.0)
        || adaptive.minimum_step_s > adaptive.initial_step_s
        || adaptive.initial_step_s > adaptive.maximum_step_s
        || !adaptive.relative_tolerance.is_finite()
        || adaptive.relative_tolerance < 0.0
        || adaptive.maximum_attempts == 0
        || adaptive.maximum_accepted_intervals == 0
        || adaptive.maximum_accepted_intervals.checked_mul(2).is_none()
    {
        return Err(EffusionTrajectoryError::Invalid(
            "positive bounded time, dimensional tolerances and work limits required",
        ));
    }
    Ok(())
}

// Shared by thermal/recoil evolution: neither a split thermal stage nor a
// rounded final tail may create or discard elapsed physical time.
pub(super) fn trial_times(
    elapsed: f64,
    next_step: f64,
    duration_s: f64,
    minimum_step_s: f64,
) -> Result<(f64, f64, f64, f64), EffusionTrajectoryError> {
    let remaining = duration_s - elapsed;
    let mut h = next_step.min(remaining);
    let tail = remaining - h;
    if tail > 0.0 && tail < minimum_step_s && 0.5 * remaining >= minimum_step_s {
        h = 0.5 * remaining;
    }
    let half = 0.5 * h;
    let end = elapsed + h;
    let midpoint = elapsed + half;
    if half <= 0.0
        || half + half != h
        || end <= elapsed
        || end > duration_s
        || midpoint <= elapsed
        || midpoint >= end
        || (h == remaining && end != duration_s)
    {
        return Err(EffusionTrajectoryError::Unrepresentable(
            "trial halves or final physical time cannot advance exactly",
        ));
    }
    Ok((h, half, midpoint, end))
}

fn integrate<'a, S, const N: usize, C: FnMut() -> bool>(
    original: &'a GasRocketState,
    source_gases: impl IntoIterator<Item = &'a ClosedGasParcel>,
    duration_s: f64,
    adaptive: EffusionAdaptiveSolve,
    mut advance: impl FnMut(
        &GasRocketState,
        f64,
        &mut C,
    ) -> Result<(GasRocketState, S), EffusionRocketError>,
    receipt_values: impl Fn(&S) -> [f64; N],
    cancelled: &mut C,
) -> Result<ControlledTrajectory<S, N>, EffusionTrajectoryError> {
    poll(cancelled)?;
    admit_adaptive(duration_s, adaptive)?;
    let species = SpeciesAxis::new(
        std::iter::once(original.reservoir()).chain(source_gases),
        cancelled,
    )?;
    let mut current = original.clone();
    let (mut elapsed, mut next_step, mut maximum_error) = (0.0, adaptive.initial_step_s, 0.0_f64);
    let (mut attempts, mut rejected, mut accepted) = (0, 0, 0);
    let mut accepted_steps = Vec::new();
    let mut totals: [Sum; N] = core::array::from_fn(|_| Sum::default());
    while elapsed < duration_s {
        poll(cancelled)?;
        if attempts == adaptive.maximum_attempts || accepted == adaptive.maximum_accepted_intervals
        {
            return Err(EffusionTrajectoryError::WorkBudget {
                attempted_intervals: attempts,
                accepted_intervals: accepted,
            });
        }
        let (h, half, midpoint_time, end) =
            trial_times(elapsed, next_step, duration_s, adaptive.minimum_step_s)?;
        attempts += 1;
        let trial = (|| {
            let mut perform = |state: &GasRocketState, dt: f64, at: f64, cancel: &mut C| {
                advance(state, dt, cancel).map_err(|error| EffusionTrajectoryError::Step {
                    elapsed_s: at,
                    duration_s: dt,
                    error,
                })
            };
            let (coarse, _) = perform(&current, h, elapsed, cancelled)?;
            let (middle, first) = perform(&current, half, elapsed, cancelled)?;
            let (fine, second) = perform(&middle, half, midpoint_time, cancelled)?;
            let error =
                normalized_error(&current, &coarse, &fine, &species, h, adaptive, cancelled)?;
            Ok::<_, EffusionTrajectoryError>((fine, [first, second], error))
        })();
        let (fine, steps, error) = match trial {
            Ok(result) => result,
            Err(error) if withdrawal(&error) && half >= adaptive.minimum_step_s => {
                rejected += 1;
                next_step = half;
                continue;
            }
            Err(error) => return Err(error),
        };
        if error > 1.0 {
            rejected += 1;
            if half < adaptive.minimum_step_s {
                return Err(EffusionTrajectoryError::MinimumStep {
                    duration_s: h,
                    normalized_error: error,
                });
            }
            next_step = half;
            continue;
        }
        for step in steps {
            poll(cancelled)?;
            for (total, value) in totals.iter_mut().zip(receipt_values(&step)) {
                total
                    .add(value)
                    .map_err(EffusionTrajectoryError::Accounting)?;
            }
            accepted_steps.push(step);
        }
        current = fine;
        elapsed = end;
        accepted += 1;
        maximum_error = maximum_error.max(error);
        next_step = if error <= 0.25 {
            (2.0 * h).min(adaptive.maximum_step_s)
        } else {
            h
        };
    }
    if elapsed != duration_s {
        return Err(EffusionTrajectoryError::Unrepresentable(
            "accepted interval did not reach requested time",
        ));
    }
    poll(cancelled)?;
    Ok(ControlledTrajectory {
        endpoint: current,
        duration_s: elapsed,
        attempted_intervals: attempts,
        rejected_intervals: rejected,
        accepted_intervals: accepted,
        accepted_steps,
        maximum_accepted_error: maximum_error,
        totals: totals.map(|sum| sum.total()),
    })
}

impl GasRocketState {
    /// Advance a complete positive interval through accepted native half steps.
    ///
    /// All outlets, normals and loads are held until the end. Source and physical
    /// refusals propagate unchanged, except the original typed gross-withdrawal
    /// limit can trigger halving. Ordinary error rejection also halves; a small
    /// accepted error (<=1/4) permits doubling. No rejected state or report enters
    /// the accepted history. Cancellation, work exhaustion or unresolved final
    /// time returns neither endpoint. A shorter final remainder is permitted
    /// once, but cannot be refined below the declared minimum trial interval.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_adaptive(
        &self,
        environment: EffusionRocketEnvironment,
        outlets: &[EffusionOutlet<'_>],
        duration_s: f64,
        endpoint: GasEnergySolve,
        solve: EffusionRocketSolve,
        adaptive: EffusionAdaptiveSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EffusionTrajectoryCandidate, EffusionTrajectoryError> {
        let window = integrate(
            self,
            std::iter::empty(),
            duration_s,
            adaptive,
            |state, dt, cancel| {
                state
                    .advance_effusion_with_environment(
                        environment,
                        outlets,
                        dt,
                        endpoint,
                        solve,
                        cancel,
                    )
                    .map(|candidate| candidate.into_parts())
            },
            values,
            &mut cancelled,
        )?;
        let t = window.totals;
        Ok(EffusionTrajectoryCandidate {
            endpoint: window.endpoint,
            step: EffusionTrajectoryStep {
                duration_s: window.duration_s,
                attempted_intervals: window.attempted_intervals,
                rejected_intervals: window.rejected_intervals,
                accepted_intervals: window.accepted_intervals,
                accepted_steps: window.accepted_steps,
                maximum_accepted_error: window.maximum_accepted_error,
                total_escaped_mass_kg: t[0],
                total_escaped_momentum_world_kg_m_s: Vec3::new(t[1], t[2], t[3]),
                total_escaped_energy_world_j: t[4],
                total_external_impulse_n_s: Vec3::new(t[5], t[6], t[7]),
                total_external_work_j: t[8],
                total_mass_residual_kg: t[9],
                total_momentum_residual_kg_m_s: Vec3::new(t[10], t[11], t[12]),
                total_energy_residual_j: t[13],
                total_gravity_impulse_uncertainty_n_s: Vec3::new(t[14], t[15], t[16]),
                total_gravity_work_uncertainty_j: t[17],
            },
        })
    }
}
