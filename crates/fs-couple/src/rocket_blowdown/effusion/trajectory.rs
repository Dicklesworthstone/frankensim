//! Bounded adaptive time advancement of the original finite-gas recoil owner.
//!
//! Each trial compares one complete step with two half steps, then accepts only
//! the two half steps. Their underlying gas update remains first order. The
//! difference estimates LOCAL truncation error; it is not a certified or global
//! error bound. Original source, inventory and physical conservation gates stay
//! active in every trial, with no acceptance of a failed coarse physical solve.

use super::{
    EffusionOutlet, EffusionRocketEnvironment, EffusionRocketError, EffusionRocketSolve,
    EffusionRocketStep, GasRocketState, Sum,
};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::{
    SpeciesNetworkError, transport::GasTransportError,
};
use fs_mbd::Vec3;
use fs_qty::MolarMass;

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
fn ratio(
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

#[allow(clippy::too_many_arguments)]
fn normalized_error(
    before: &GasRocketState,
    coarse: &GasRocketState,
    fine: &GasRocketState,
    original: &GasRocketState,
    molar_masses: &[MolarMass],
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
    let catalog = original.reservoir().model().thermodynamics().components();
    let mut mass_fractions = |state: &GasRocketState| {
        poll(cancelled)?;
        let mixture = state.reservoir().model().thermodynamics();
        if !mixture
            .components()
            .iter()
            .map(|s| s.species())
            .eq(catalog.iter().map(|s| s.species()))
        {
            return Err(EffusionTrajectoryError::Invalid(
                "effusion changed its original species axis",
            ));
        }
        mixture
            .mole_composition()
            .to_mass_fractions(molar_masses)
            .map_err(|_| {
                EffusionTrajectoryError::Unrepresentable("original species mass fractions")
            })
    };
    let fractions = [
        mass_fractions(before)?,
        mass_fractions(coarse)?,
        mass_fractions(fine)?,
    ];
    for i in 0..molar_masses.len() {
        poll(cancelled)?;
        let old = before.reservoir().mass_kg() * fractions[0].fractions()[i];
        let full = coarse.reservoir().mass_kg() * fractions[1].fractions()[i];
        let half = fine.reservoir().mass_kg() * fractions[2].fractions()[i];
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
        poll(&mut cancelled)?;
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
        let molar_masses: Vec<_> = self
            .reservoir()
            .model()
            .thermodynamics()
            .components()
            .iter()
            .map(|s| s.molar_mass())
            .collect();
        let mut current = self.clone();
        let (mut elapsed, mut next_step, mut maximum_error) =
            (0.0, adaptive.initial_step_s, 0.0_f64);
        let (mut attempts, mut rejected, mut accepted) = (0, 0, 0);
        let mut accepted_steps = Vec::new();
        let mut totals: [Sum; 18] = core::array::from_fn(|_| Sum::default());
        while elapsed < duration_s {
            poll(&mut cancelled)?;
            if attempts == adaptive.maximum_attempts
                || accepted == adaptive.maximum_accepted_intervals
            {
                return Err(EffusionTrajectoryError::WorkBudget {
                    attempted_intervals: attempts,
                    accepted_intervals: accepted,
                });
            }
            let remaining = duration_s - elapsed;
            let mut h = next_step.min(remaining);
            // A rounded clock can put the last requested endpoint just beyond
            // the held step. Split the remaining pair instead of manufacturing
            // an unrepresentably small final gas withdrawal or skipping time.
            let tail = remaining - h;
            if tail > 0.0
                && tail < adaptive.minimum_step_s
                && 0.5 * remaining >= adaptive.minimum_step_s
            {
                h = 0.5 * remaining;
            }
            let half = 0.5 * h;
            let end = elapsed + h;
            let midpoint_time = elapsed + half;
            if half <= 0.0
                || half + half != h
                || end <= elapsed
                || end > duration_s
                || midpoint_time <= elapsed
                || midpoint_time >= end
                || (h == remaining && end != duration_s)
            {
                return Err(EffusionTrajectoryError::Unrepresentable(
                    "trial halves or final physical time cannot advance exactly",
                ));
            }
            attempts += 1;
            let trial = (|| {
                let advance = |state: &GasRocketState, dt: f64, at: f64, cancel: &mut _| {
                    state
                        .advance_effusion_with_environment(
                            environment,
                            outlets,
                            dt,
                            endpoint,
                            solve,
                            cancel,
                        )
                        .map_err(|error| EffusionTrajectoryError::Step {
                            elapsed_s: at,
                            duration_s: dt,
                            error,
                        })
                };
                let (coarse, _) = advance(&current, h, elapsed, &mut cancelled)?.into_parts();
                let (middle, first) =
                    advance(&current, half, elapsed, &mut cancelled)?.into_parts();
                let (fine, second) =
                    advance(&middle, half, midpoint_time, &mut cancelled)?.into_parts();
                let error = normalized_error(
                    &current,
                    &coarse,
                    &fine,
                    self,
                    &molar_masses,
                    h,
                    adaptive,
                    &mut cancelled,
                )?;
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
                poll(&mut cancelled)?;
                for (total, value) in totals.iter_mut().zip(values(&step)) {
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
        poll(&mut cancelled)?;
        let t = totals.map(|sum| sum.total());
        Ok(EffusionTrajectoryCandidate {
            endpoint: current,
            step: EffusionTrajectoryStep {
                duration_s: elapsed,
                attempted_intervals: attempts,
                rejected_intervals: rejected,
                accepted_intervals: accepted,
                accepted_steps,
                maximum_accepted_error: maximum_error,
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
