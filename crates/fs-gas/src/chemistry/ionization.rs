//! Conservative two-temperature evolution of one ground-state atomic channel.
//!
//! The L1 source owns `e + A <=> A+ + 2e`, including its Maxwellian event
//! cross section and original-source reverse detailed balance. This L3 owner
//! evolves the reaction extent in a rigid parcel, without projecting to LTE.
//! For the admitted monatomic channel, one heavy particle remains one heavy
//! particle and both have Cv=3R/2. Chemistry therefore preserves heavy sensible
//! energy while electrons supply the original ionization energy AND the
//! sensible energy of the new electron. Recombination reverses that transfer.
//!
//! With extent xi [mol], the exact electron budget is
//! `3R/2 * (ne0+xi) * Te = 3R/2 * ne0 * Te0 + Pe*t - chi*xi - Q`.
//! The optional shared collision sources supply the equal-and-opposite
//! electron-to-heavy transfer Q. Reaction and heat transfer are solved
//! simultaneously, retaining the changing electron heat capacity.
//! Heavy sensible heating uses the same source-difference caloric owner as
//! frozen electron/heavy exchange. Separate, constant externally calculated
//! powers may act on either pool. A power sampled from an electrical solve is
//! held over this explicit interval; it is not a continuous voltage solve.
//!
//! Adaptive backward Euler with step doubling solves only ONE scalar extent.
//! Analytic constrained derivatives, positivity and source-domain brackets,
//! actual two-pool entropy, and source energy balances gate publication. No
//! nonlinear trial allocates. Collision coefficients are refreshed at each
//! substep's actual populations and temperatures, then held in its implicit
//! solve; step doubling controls this first-order approximation. Radiative
//! recombination, excited states, secondary-electron spectra and spatial
//! transport remain outside this channel. Estimated local errors are not
//! empirical or global error certificates.

use super::energy::{ElectronHeavyCollisionSources, ExchangeConductance, PreparedExchange};
use fs_qty::{Amount, Energy, Power, Temperature, Time};
use fs_thermochem::equilibrium::composition::{
    FrozenIonizedGasStateErrorV1, FrozenIonizedGasStateV1,
};
use fs_thermochem::kinetics::atomic_ionization::{AtomicIonizationErrorV1, AtomicIonizationV1};
use fs_thermochem::kinetics::{KineticsErrorV1, admit_caloric_interval};
use fs_thermochem::{SpeciesId, ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

mod solve;
use solve::{Point, Problem};

/// Source-domain bounds and explicit numerical work/error controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AtomicIonizationSolve {
    /// Minimum source-supported heavy temperature [K].
    pub minimum_heavy_temperature_k: f64,
    /// Maximum source-supported heavy temperature [K].
    pub maximum_heavy_temperature_k: f64,
    /// Minimum source-supported electron temperature [K].
    pub minimum_electron_temperature_k: f64,
    /// Maximum source-supported electron temperature [K].
    pub maximum_electron_temperature_k: f64,
    /// Estimated local relative amount/temperature error tolerance.
    pub relative_tolerance: f64,
    /// Estimated local absolute amount error scale [mol].
    pub absolute_amount_tolerance_mol: f64,
    /// Estimated local absolute temperature error scale [K].
    pub temperature_tolerance_k: f64,
    /// Original-source pool and total energy closure tolerance [J].
    pub energy_tolerance_j: f64,
    /// Numerical allowance in each step's physical entropy inequality [J/K].
    pub entropy_tolerance_j_k: f64,
    /// Carrier-scaled charge admission at publication, not diluted by neutrals.
    pub relative_charge_tolerance: f64,
    /// First proposed step [s].
    pub initial_step_s: f64,
    /// Smallest retry [s]; a shorter final remainder is permitted.
    pub minimum_step_s: f64,
    /// Largest proposed step [s].
    pub maximum_step_s: f64,
    /// Maximum accepted plus rejected attempts.
    pub maximum_steps: usize,
    /// Maximum safeguarded scalar or caloric iterations per solve.
    pub maximum_iterations: usize,
}
impl AtomicIonizationSolve {
    /// Set a common source interval; callers may then select separate bounds.
    /// Numerical defaults express tolerances, not measured physical accuracy.
    pub fn new(minimum: Temperature, maximum: Temperature) -> Self {
        Self {
            minimum_heavy_temperature_k: minimum.value(),
            maximum_heavy_temperature_k: maximum.value(),
            minimum_electron_temperature_k: minimum.value(),
            maximum_electron_temperature_k: maximum.value(),
            relative_tolerance: 1e-5,
            absolute_amount_tolerance_mol: 1e-14,
            temperature_tolerance_k: 1e-6,
            energy_tolerance_j: 1e-8,
            entropy_tolerance_j_k: 1e-9,
            relative_charge_tolerance: 1e-10,
            initial_step_s: 1e-10,
            minimum_step_s: 1e-24,
            maximum_step_s: 1e-3,
            maximum_steps: 100_000,
            maximum_iterations: 96,
        }
    }
    fn validate(self) -> Result<(), AtomicIonizationError> {
        if [
            self.minimum_heavy_temperature_k,
            self.maximum_heavy_temperature_k,
            self.minimum_electron_temperature_k,
            self.maximum_electron_temperature_k,
            self.relative_tolerance,
            self.absolute_amount_tolerance_mol,
            self.temperature_tolerance_k,
            self.energy_tolerance_j,
            self.relative_charge_tolerance,
            self.initial_step_s,
            self.minimum_step_s,
            self.maximum_step_s,
        ]
        .iter()
        .any(|x| !positive(*x))
            || self.minimum_heavy_temperature_k >= self.maximum_heavy_temperature_k
            || self.minimum_electron_temperature_k >= self.maximum_electron_temperature_k
            || self.relative_tolerance > 0.1
            || self.relative_charge_tolerance >= 1.0
            || !self.entropy_tolerance_j_k.is_finite()
            || self.entropy_tolerance_j_k < 0.0
            || self.minimum_step_s > self.initial_step_s
            || self.initial_step_s > self.maximum_step_s
            || !(1..=1_000_000).contains(&self.maximum_steps)
            || !(8..=256).contains(&self.maximum_iterations)
        {
            return Err(AtomicIonizationError::Invalid(
                "atomic integration controls",
            ));
        }
        Ok(())
    }
}

/// Atomic refusal: caller-owned source and state are never mutated.
#[derive(Clone, Debug, PartialEq)]
pub enum AtomicIonizationError {
    /// Incompatible state/source binding or invalid input/control.
    Invalid(&'static str),
    /// Original channel or caloric admission refused.
    Kinetics(KineticsErrorV1),
    /// Atomic level statistics, detailed balance, or impact source refused.
    Channel(AtomicIonizationErrorV1),
    /// Original source species refused a physical temperature.
    Species {
        species: SpeciesId,
        source: ThermochemErrorV1,
    },
    /// Shared heavy-pool caloric evaluator refused.
    Caloric(super::energy::ElectronHeavyEnergyError),
    /// A physical source/positivity/energy domain has been exhausted.
    OutsideDomain,
    /// Bounded nonlinear solve failed.
    Nonconvergence,
    /// Explicit step or minimum-step budget exhausted.
    StepBudget,
    /// A needed arithmetic value was unrepresentable.
    Arithmetic(&'static str),
    /// Source energy or entropy inequality failed its numerical allowance.
    Balance(&'static str),
    /// Shared final physical-state publication refused.
    State(FrozenIonizedGasStateErrorV1),
    /// Cancellation before publication.
    Cancelled,
}
impl core::fmt::Display for AtomicIonizationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "two-temperature atomic ionization: {self:?}")
    }
}
impl std::error::Error for AtomicIonizationError {}

/// Actual source balances and bounded adaptive work, not measurement accuracy.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AtomicIonizationBalance {
    /// Integrated external electron power [J].
    pub electron_deposited_energy_j: f64,
    /// Integrated external heavy power [J].
    pub heavy_deposited_energy_j: f64,
    /// Original reaction formation energy stored, chi*xi [J].
    pub chemical_storage_change_j: f64,
    /// Signed cumulative electron-to-heavy thermal transfer [J].
    pub electron_to_heavy_energy_j: f64,
    /// Coulomb part of the cumulative electron-to-heavy transfer [J].
    pub electron_to_heavy_coulomb_energy_j: f64,
    /// Neutral-collision part of the cumulative electron-to-heavy transfer [J].
    pub electron_to_heavy_neutral_energy_j: f64,
    /// Electron sensible + chemical storage + transfer - electron work [J].
    pub electron_sensible_energy_residual_j: f64,
    /// Heavy sensible change - transfer - heavy work [J].
    pub heavy_sensible_energy_residual_j: f64,
    /// Original-source total energy increment - all external work [J].
    pub energy_residual_j: f64,
    /// Actual ideal-mixture entropy difference, including changing populations [J/K].
    pub entropy_change_j_k: f64,
    /// Smallest accepted entropy increment above implicit power/temperature sum [J/K].
    pub minimum_step_entropy_excess_j_k: f64,
    /// Largest accepted full/two-half normalized local error estimate.
    pub maximum_normalized_local_error: f64,
    /// Accepted complete steps; each publishes its two-half-step solution.
    pub accepted_steps: usize,
    /// Rejected attempts.
    pub rejected_steps: usize,
}

/// Source-retained endpoint, reusable directly by transport and electrical owners.
#[derive(Clone, Debug, PartialEq)]
pub struct AtomicIonizationResult {
    initial: FrozenIonizedGasStateV1,
    state: FrozenIonizedGasStateV1,
    channel: AtomicIonizationV1,
    collisions: Option<ElectronHeavyCollisionSources>,
    controls: AtomicIonizationSolve,
    elapsed: Time,
    extent: Amount,
    balance: AtomicIonizationBalance,
}
impl AtomicIonizationResult {
    /// Actual final populations, separate temperatures and original thermodynamics.
    pub const fn state(&self) -> &FrozenIonizedGasStateV1 {
        &self.state
    }
    /// Original source-bound initial state.
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.initial
    }
    /// Original atomic channel, cross section, and statistical weights.
    pub const fn channel(&self) -> &AtomicIonizationV1 {
        &self.channel
    }
    /// Original collision laws and initial Coulomb population receipt. Changed
    /// endpoint populations require explicit rebinding for material transport.
    pub const fn collision_sources(&self) -> Option<&ElectronHeavyCollisionSources> {
        self.collisions.as_ref()
    }
    /// Numerical controls and source intervals used.
    pub const fn controls(&self) -> AtomicIonizationSolve {
        self.controls
    }
    /// Completed physical duration [s].
    pub const fn elapsed_time(&self) -> Time {
        self.elapsed
    }
    /// Signed extent [mol]; negative means net three-body recombination.
    pub const fn reaction_extent(&self) -> Amount {
        self.extent
    }
    /// Actual source balances and integration work.
    pub const fn balance(&self) -> AtomicIonizationBalance {
        self.balance
    }
    /// Original formation energy change [J].
    pub const fn chemical_storage_change(&self) -> Energy {
        Energy::new(self.balance.chemical_storage_change_j)
    }
}

/// Evolve one atomic collisional channel at fixed volume with held pool powers.
///
/// Explicit zero electrons stay zero; this operator never manufactures a seed.
/// Electron power without an electron population refuses. The original bank
/// must exactly match the channel. Inert species keep their amounts, contribute
/// their original heavy calorics, and need not be evaluated at electron T.
pub fn advance_atomic_ionization(
    initial: &FrozenIonizedGasStateV1,
    channel: &AtomicIonizationV1,
    electron_power: Power,
    heavy_power: Power,
    duration: Time,
    controls: AtomicIonizationSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<AtomicIonizationResult, AtomicIonizationError> {
    advance(
        initial,
        channel,
        None,
        electron_power,
        heavy_power,
        duration,
        controls,
        cancelled,
    )
}

/// Evolve atomic chemistry and electron/heavy thermal relaxation in one
/// conservative implicit solve. Sources for neutrals created by recombination
/// are required even when their initial abundance is zero. The original
/// Coulomb receipt must match the initial population; subsequent moments are
/// explicitly refreshed at actual populations, density and temperatures.
///
/// Held powers may already include field-derived friction and signed drift
/// corrections. Do not add total E.J or zero-drift exchange to them again.
pub fn advance_atomic_ionization_with_collisions(
    initial: &FrozenIonizedGasStateV1,
    channel: &AtomicIonizationV1,
    collisions: &ElectronHeavyCollisionSources,
    electron_power: Power,
    heavy_power: Power,
    duration: Time,
    controls: AtomicIonizationSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<AtomicIonizationResult, AtomicIonizationError> {
    advance(
        initial,
        channel,
        Some(collisions),
        electron_power,
        heavy_power,
        duration,
        controls,
        cancelled,
    )
}

fn advance(
    initial: &FrozenIonizedGasStateV1,
    channel: &AtomicIonizationV1,
    collisions: Option<&ElectronHeavyCollisionSources>,
    electron_power: Power,
    heavy_power: Power,
    duration: Time,
    controls: AtomicIonizationSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<AtomicIonizationResult, AtomicIonizationError> {
    controls.validate()?;
    poll(&mut cancelled)?;
    if !duration.value().is_finite()
        || duration.value() < 0.0
        || !electron_power.value().is_finite()
        || !heavy_power.value().is_finite()
    {
        return Err(AtomicIonizationError::Invalid(
            "finite powers and nonnegative duration",
        ));
    }
    let problem = Problem::new(
        initial,
        channel,
        collisions,
        [heavy_power.value(), electron_power.value()],
        controls,
        &mut cancelled,
    )?;
    let mut point = Point {
        time: 0.0,
        extent: 0.0,
        heavy_temperature: initial.heavy_temperature().value(),
        electron_temperature: initial.electron_temperature().value(),
        exchange: 0.0,
        coulomb_exchange: 0.0,
        electron_temperature_extent_derivative: 0.0,
    };
    let start = point;
    let mut balance = AtomicIonizationBalance::default();
    let mut step = controls.initial_step_s;
    // The accepted endpoint's source evaluation is exactly the next step's
    // starting coefficient. Reuse it across accepted steps and rejected
    // retries instead of repeating Maxwellian quadrature and Coulomb refresh.
    let mut coupling = problem.coupling(&point, &mut cancelled)?;
    while point.time < duration.value() {
        poll(&mut cancelled)?;
        if balance.accepted_steps + balance.rejected_steps >= controls.maximum_steps {
            return Err(AtomicIonizationError::StepBudget);
        }
        let remaining = duration.value() - point.time;
        let dt = step.min(remaining);
        let end = if dt == remaining {
            duration.value()
        } else {
            point.time + dt
        };
        let middle = point.time + 0.5 * (end - point.time);
        if middle <= point.time || middle >= end {
            return Err(AtomicIonizationError::StepBudget);
        }
        let attempt = (|| {
            let full = problem.step(&point, end, coupling, &mut cancelled)?;
            let half = problem.step(&point, middle, coupling, &mut cancelled)?;
            let middle_coupling = problem.coupling(&half, &mut cancelled)?;
            let final_point = problem.step(&half, end, middle_coupling, &mut cancelled)?;
            // Admission uses the endpoint's real abundances, not the lagged
            // state that supplied this first-order substep's conductance.
            let endpoint_coupling = problem.coupling(&final_point, &mut cancelled)?;
            let error = problem.error(&full, &final_point);
            let entropy = problem.entropy_change(&point, &final_point, &mut cancelled)?;
            // For backward Euler, concavity gives Delta S >= sum(h*P/T_end)
            // plus nonnegative chemical production. Audit both actual halves.
            let supplied = problem.entropy_supply(middle - point.time, &half)
                + problem.entropy_supply(end - middle, &final_point);
            let first_excess = problem.entropy_change(&point, &half, &mut cancelled)?
                - problem.entropy_supply(middle - point.time, &half);
            let second_excess = problem.entropy_change(&half, &final_point, &mut cancelled)?
                - problem.entropy_supply(end - middle, &final_point);
            if !error.is_finite() || !entropy.is_finite() || !supplied.is_finite() {
                return Err(AtomicIonizationError::Arithmetic("step error or entropy"));
            }
            if first_excess < -controls.entropy_tolerance_j_k
                || second_excess < -controls.entropy_tolerance_j_k
                || entropy - supplied < -controls.entropy_tolerance_j_k
            {
                return Err(AtomicIonizationError::Balance(
                    "chemical entropy inequality",
                ));
            }
            Ok((
                final_point,
                error,
                (entropy - supplied).min(first_excess).min(second_excess),
                endpoint_coupling,
            ))
        })();
        match attempt {
            Ok((next, error, excess, next_coupling)) if error <= 1.0 => {
                balance.maximum_normalized_local_error =
                    balance.maximum_normalized_local_error.max(error);
                balance.minimum_step_entropy_excess_j_k = if balance.accepted_steps == 0 {
                    excess
                } else {
                    balance.minimum_step_entropy_excess_j_k.min(excess)
                };
                balance.accepted_steps += 1;
                point = next;
                coupling = next_coupling;
                let factor = if error == 0.0 {
                    2.0
                } else {
                    (0.85 / fs_math::det::sqrt(error)).clamp(0.25, 2.0)
                };
                step = (dt * factor)
                    .min(controls.maximum_step_s)
                    .max(controls.minimum_step_s);
            }
            Ok((_, error, _, _)) => {
                balance.rejected_steps += 1;
                step = dt * (0.8 / fs_math::det::sqrt(error)).clamp(0.1, 0.5);
                if step < controls.minimum_step_s {
                    return Err(AtomicIonizationError::StepBudget);
                }
            }
            Err(
                error @ (AtomicIonizationError::OutsideDomain
                | AtomicIonizationError::Nonconvergence
                | AtomicIonizationError::Balance(_)),
            ) => {
                balance.rejected_steps += 1;
                step = dt * 0.5;
                if step < controls.minimum_step_s {
                    return Err(error);
                }
            }
            Err(error) => return Err(error),
        }
    }
    let state = problem.publish(&point, &mut balance, &mut cancelled)?;
    balance.entropy_change_j_k = problem.entropy_change(&start, &point, &mut cancelled)?;
    if electron_power.value() == 0.0
        && heavy_power.value() == 0.0
        && balance.entropy_change_j_k < -controls.entropy_tolerance_j_k
    {
        return Err(AtomicIonizationError::Balance(
            "complete isolated chemical entropy",
        ));
    }
    poll(&mut cancelled)?;
    Ok(AtomicIonizationResult {
        initial: initial.clone(),
        state,
        channel: channel.clone(),
        collisions: collisions.cloned(),
        controls,
        elapsed: duration,
        extent: Amount::new(point.extent),
        balance,
    })
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AtomicIonizationError> {
    if cancelled() {
        Err(AtomicIonizationError::Cancelled)
    } else {
        Ok(())
    }
}
fn kinetic_error(error: KineticsErrorV1) -> AtomicIonizationError {
    match error {
        KineticsErrorV1::Cancelled => AtomicIonizationError::Cancelled,
        error => AtomicIonizationError::Kinetics(error),
    }
}
fn channel_error(error: AtomicIonizationErrorV1) -> AtomicIonizationError {
    match error {
        AtomicIonizationErrorV1::Cancelled => AtomicIonizationError::Cancelled,
        error => AtomicIonizationError::Channel(error),
    }
}

#[cfg(test)]
mod tests;
