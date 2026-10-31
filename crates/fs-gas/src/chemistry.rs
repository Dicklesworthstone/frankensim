//! Adaptive finite-rate chemistry in a sealed, rigid ideal-gas volume.
//!
//! The immutable L1 mechanism owns reaction laws and thermodynamics. This L3
//! owner advances reaction extents with damped backward Euler or the optional
//! two-stage, second-order L-stable SDIRK method, both with step doubling.
//! Each trial inverts the same total internal energy, including formation energy
//! exactly once. Analytic constrained Jacobians handle stiffness; no explicit
//! positivity clipping, trace seeding or equilibrium projection is used.
//!
//! Step error is an estimated local error, not a global error certificate.
//! Accepted endpoints are nonnegative, conserve the original element and charge
//! inventories and U within the declared tolerances, and pass the actual ideal-mixture
//! entropy check. No ignition, flame, transport or empirical mechanism claim is
//! implied by a successful zero-dimensional homogeneous reactor integration.

use fs_qty::{Amount, Energy, Mass, Pressure, Temperature, Time, Volume};
use fs_thermochem::kinetics::{
    IdealGasKineticsV1, KineticThermodynamicsV1, KineticWorkspaceV1, KineticsErrorV1,
};
use fs_thermochem::{Nasa9EvaluationReceiptV1, SpeciesId};

mod prepared;
mod solve;

pub use prepared::{IsochoricKineticsView, IsochoricReactor};

/// Time method for the homogeneous constrained reactor, not its outer coupling.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum IsochoricKineticsMethod {
    /// First-order backward Euler; retained as the default.
    #[default]
    BackwardEuler,
    /// Alexander's two-stage, second-order L-stable SDIRK method, with
    /// gamma = 1 - 1/sqrt(2). Both stages use the same original fixed U,V.
    ///
    /// Reference: R. Alexander, SIAM J. Numer. Anal. 14(6), 1006–1021 (1977),
    /// <https://doi.org/10.1137/0714068>. The order is for the smooth reactor
    /// ODE; it does not establish global order for split transport or fire
    /// coupling. Positivity and entropy remain actual trial acceptance gates,
    /// not unconditional properties of this Runge–Kutta formula.
    Sdirk2,
}

impl IsochoricKineticsMethod {
    fn error_divisor(self) -> f64 {
        match self {
            Self::BackwardEuler => 1.0,
            Self::Sdirk2 => 3.0,
        }
    }

    fn error_root(self, error: f64) -> f64 {
        match self {
            Self::BackwardEuler => fs_math::det::sqrt(error),
            Self::Sdirk2 => fs_math::det::pow(error, 1.0 / 3.0),
        }
    }
}

/// Explicit physical tolerances and bounded adaptive integration work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsochoricKineticsSolve {
    /// Reactor integration method; backward Euler unless explicitly selected.
    pub method: IsochoricKineticsMethod,
    /// Common lower source-supported temperature, K.
    pub minimum_temperature_k: f64,
    /// Common upper source-supported temperature, K.
    pub maximum_temperature_k: f64,
    /// Absolute tolerance on the complete parcel internal energy, J.
    pub energy_tolerance_j: f64,
    /// Estimated local relative species and temperature tolerance.
    pub relative_tolerance: f64,
    /// Estimated local absolute amount tolerance per species, mol.
    pub absolute_amount_tolerance_mol: f64,
    /// Maximum numerical entropy decrease across any accepted step and the
    /// complete integration, J/K. This is not source-data uncertainty.
    pub entropy_tolerance_j_k: f64,
    /// Maximum relative element, charge-carrier and molecular-mass discrepancy.
    pub relative_conservation_tolerance: f64,
    /// Initial proposed step, s; the adaptive controller may reduce it.
    pub initial_step_s: f64,
    /// Maximum accepted step, s.
    pub maximum_step_s: f64,
    /// Minimum proposed step after a rejection, s.
    pub minimum_step_s: f64,
    /// Maximum total step attempts, including rejections; at most100000.
    pub maximum_steps: usize,
    /// Maximum Newton iterations for each implicit stage; at most64.
    pub maximum_newton_iterations: usize,
    /// Maximum positivity/residual line-search trials; at most64.
    pub maximum_backtracks: usize,
    /// Maximum temperature iterations per constrained trial; at most256.
    pub maximum_temperature_iterations: usize,
}
impl Default for IsochoricKineticsSolve {
    fn default() -> Self {
        Self {
            method: IsochoricKineticsMethod::BackwardEuler,
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 6000.0,
            energy_tolerance_j: 1e-7,
            relative_tolerance: 1e-6,
            absolute_amount_tolerance_mol: 1e-12,
            entropy_tolerance_j_k: 1e-9,
            relative_conservation_tolerance: 1e-10,
            initial_step_s: 1e-6,
            maximum_step_s: 1.0,
            minimum_step_s: 1e-16,
            maximum_steps: 10000,
            maximum_newton_iterations: 24,
            maximum_backtracks: 24,
            maximum_temperature_iterations: 64,
        }
    }
}
impl IsochoricKineticsSolve {
    fn validate(self) -> Result<(), IsochoricKineticsError> {
        if !positive(self.minimum_temperature_k)
            || !positive(self.maximum_temperature_k)
            || self.maximum_temperature_k <= self.minimum_temperature_k
            || !positive(self.energy_tolerance_j)
            || !positive(self.relative_tolerance)
            || self.relative_tolerance > 0.1
            || !positive(self.absolute_amount_tolerance_mol)
            || !self.entropy_tolerance_j_k.is_finite()
            || self.entropy_tolerance_j_k < 0.0
            || !positive(self.relative_conservation_tolerance)
            || self.relative_conservation_tolerance > 1e-6
            || !positive(self.initial_step_s)
            || !positive(self.maximum_step_s)
            || !positive(self.minimum_step_s)
            || self.minimum_step_s > self.initial_step_s
            || self.initial_step_s > self.maximum_step_s
            || !(1..=100000).contains(&self.maximum_steps)
            || !(1..=64).contains(&self.maximum_newton_iterations)
            || !(1..=64).contains(&self.maximum_backtracks)
            || !(1..=256).contains(&self.maximum_temperature_iterations)
        {
            return Err(IsochoricKineticsError::Invalid(
                "invalid kinetic integration controls",
            ));
        }
        Ok(())
    }
}

/// Actual final balances and bounded work, not uncertainty or global error bounds.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IsochoricKineticsBalance {
    /// Final-minus-initial original source molecular mass, kg.
    pub mass_residual_kg: f64,
    /// Maximum relative residual of an original elemental inventory row.
    pub maximum_relative_element_residual: f64,
    /// Largest charge-inventory drift across accepted steps and publication,
    /// divided by the larger initial/current sum of absolute charged amounts.
    /// Neutral species never dilute this residual. The original net charge is
    /// preserved, including a nonzero authored feed; no electrostatic field or
    /// physical admissibility of a net-charged homogeneous parcel is inferred.
    pub maximum_relative_charge_residual: f64,
    /// Final actual internal energy minus the one conserved target, J.
    pub energy_residual_j: f64,
    /// Actual final-minus-initial ideal-mixture entropy, J/K.
    pub entropy_change_j_k: f64,
    /// Smallest accepted step entropy increment, J/K; zero for zero duration.
    pub minimum_step_entropy_change_j_k: f64,
    /// Maximum admitted dimensionless local species/temperature error estimate.
    pub maximum_estimated_local_error: f64,
    /// Accepted complete steps (each publishes the two-half-step result).
    pub accepted_steps: usize,
    /// Rejected proposed steps.
    pub rejected_steps: usize,
    /// Actual nonlinear iterations across all accepted and rejected stages.
    pub nonlinear_iterations: usize,
    /// Actual caloric trial evaluations across all stages.
    pub temperature_evaluations: usize,
    /// Newton corrections recovered by the prepared full reaction solve after
    /// the reduced solve missed the unchanged linear residual tolerance.
    /// Includes corrections in stages that are subsequently rejected.
    pub direct_fallbacks: usize,
}

/// Fully accepted finite-duration endpoint retaining the exact mechanism.
#[derive(Clone, Debug, PartialEq)]
pub struct IsochoricKineticsResult {
    mechanism: IdealGasKineticsV1,
    amounts: Vec<Amount>,
    temperature: Temperature,
    pressure: Pressure,
    internal_energy: Energy,
    target_internal_energy: Energy,
    mass: Mass,
    volume: Volume,
    elapsed_time: Time,
    balance: IsochoricKineticsBalance,
    sources: Vec<Nasa9EvaluationReceiptV1>,
}
impl IsochoricKineticsResult {
    /// Original immutable mechanism, kinetic laws and source bank.
    pub const fn mechanism(&self) -> &IdealGasKineticsV1 {
        &self.mechanism
    }
    /// Final nonnegative amounts in the mechanism's canonical species order.
    pub fn amounts(&self) -> &[Amount] {
        &self.amounts
    }
    /// Accepted final source temperature.
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    /// Final pressure from the actual molecular count and temperature.
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    /// Actual final total internal energy, J.
    pub const fn internal_energy(&self) -> Energy {
        self.internal_energy
    }
    /// Original single energy target, J.
    pub const fn target_internal_energy(&self) -> Energy {
        self.target_internal_energy
    }
    /// Final source molecular mass, kg.
    pub const fn mass(&self) -> Mass {
        self.mass
    }
    /// Fixed original volume, m3.
    pub const fn volume(&self) -> Volume {
        self.volume
    }
    /// Complete requested duration; partial integrations never publish a result.
    pub const fn elapsed_time(&self) -> Time {
        self.elapsed_time
    }
    /// Actual balances and solver work.
    pub const fn balance(&self) -> IsochoricKineticsBalance {
        self.balance
    }
    /// Accepted original NASA-9 evaluation receipts, including possible products.
    pub fn sources(&self) -> &[Nasa9EvaluationReceiptV1] {
        &self.sources
    }
}

/// Refusal without mutating caller state or publishing a partial integration.
#[derive(Clone, Debug, PartialEq)]
pub enum IsochoricKineticsError {
    /// Invalid input identity, extensive state or numerical controls.
    Invalid(&'static str),
    /// Original source or constitutive kinetics refusal.
    Chemistry(KineticsErrorV1),
    /// The declared source bracket cannot represent the prescribed total U.
    EnergyNotBracketed {
        lower_energy_j: f64,
        upper_energy_j: f64,
    },
    /// Temperature root did not close actual energy within its work budget.
    TemperatureNoConvergence { iterations: usize, residual_j: f64 },
    /// No admitted Newton step or adaptive step remains.
    StepFailure {
        elapsed_time_s: f64,
        step_s: f64,
        reason: &'static str,
    },
    /// Requested duration was not completed within the explicit attempt budget.
    WorkLimit {
        elapsed_time_s: f64,
        accepted_steps: usize,
        rejected_steps: usize,
    },
    /// Actual conservation or entropy balances exceeded the declared tolerances.
    Balance(IsochoricKineticsBalance),
    /// Bounded workspace/result reservation failed.
    AllocationRefused,
    /// Cancellation before final publication.
    Cancelled,
}
impl core::fmt::Display for IsochoricKineticsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "isochoric kinetic reactor: {self:?}")
    }
}
impl std::error::Error for IsochoricKineticsError {}
impl From<KineticsErrorV1> for IsochoricKineticsError {
    fn from(error: KineticsErrorV1) -> Self {
        match error {
            KineticsErrorV1::Cancelled => Self::Cancelled,
            other => Self::Chemistry(other),
        }
    }
}

/// Advance a caller-authored mechanism for an actual duration at fixed U,V.
/// The initial source temperature supplies the original total internal energy.
#[allow(clippy::too_many_arguments)]
pub fn advance_adiabatic(
    mechanism: &IdealGasKineticsV1,
    feed: &[(SpeciesId, Amount)],
    volume: Volume,
    initial_temperature: Temperature,
    duration: Time,
    controls: IsochoricKineticsSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<IsochoricKineticsResult, IsochoricKineticsError> {
    advance(
        mechanism,
        feed,
        None,
        volume,
        initial_temperature,
        duration,
        controls,
        cancelled,
    )
}

/// Same fixed-volume evolution from one explicit original-source energy target.
/// The supplied temperature is a bounded root hint; initial composition is
/// retained and resolved to the target before kinetic time begins.
#[allow(clippy::too_many_arguments)]
pub fn advance_adiabatic_energy(
    mechanism: &IdealGasKineticsV1,
    feed: &[(SpeciesId, Amount)],
    target: Energy,
    volume: Volume,
    initial_temperature: Temperature,
    duration: Time,
    controls: IsochoricKineticsSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<IsochoricKineticsResult, IsochoricKineticsError> {
    advance(
        mechanism,
        feed,
        Some(target),
        volume,
        initial_temperature,
        duration,
        controls,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn advance(
    mechanism: &IdealGasKineticsV1,
    feed: &[(SpeciesId, Amount)],
    target: Option<Energy>,
    volume: Volume,
    initial_temperature: Temperature,
    duration: Time,
    controls: IsochoricKineticsSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<IsochoricKineticsResult, IsochoricKineticsError> {
    poll(&mut cancelled)?;
    let amounts = feed_amounts(mechanism, feed)?;
    let mut prepared = IsochoricReactor::new(mechanism.clone(), controls, &mut cancelled)?;
    advance_prepared(
        &mut prepared,
        &amounts,
        target,
        volume,
        initial_temperature,
        duration,
        &mut cancelled,
    )?
    .to_owned(cancelled)
}

#[allow(clippy::too_many_arguments)]
fn advance_prepared<'a>(
    prepared: &'a mut IsochoricReactor,
    feed: &[Amount],
    target: Option<Energy>,
    volume: Volume,
    initial_temperature: Temperature,
    duration: Time,
    mut cancelled: impl FnMut() -> bool,
) -> Result<IsochoricKineticsView<'a>, IsochoricKineticsError> {
    poll(&mut cancelled)?;
    let controls = prepared.controls();
    if !positive(volume.value())
        || !duration.value().is_finite()
        || duration.value() < 0.0
        || !initial_temperature.value().is_finite()
        || !(controls.minimum_temperature_k..=controls.maximum_temperature_k)
            .contains(&initial_temperature.value())
        || target.is_some_and(|u| !u.value().is_finite())
    {
        return Err(IsochoricKineticsError::Invalid(
            "invalid reactor volume, energy, temperature or duration",
        ));
    }
    let IsochoricReactor {
        reactor,
        current,
        full,
        half,
        fine,
        stage_forcing,
        current_concentrations,
        full_concentrations,
        amounts,
    } = prepared;
    let target = reactor.initialize(feed, volume, initial_temperature, target, &mut cancelled)?;
    current.fill(0.0);
    full.fill(0.0);
    half.fill(0.0);
    fine.fill(0.0);
    current_concentrations.copy_from_slice(&reactor.initial);
    let initial = reactor.resolve(current, initial_temperature.value(), &mut cancelled)?;
    let mut accepted = initial;
    let mut elapsed = 0.0;
    let mut step = controls.initial_step_s.min(duration.value());
    let mut balance = IsochoricKineticsBalance::default();
    while elapsed < duration.value() {
        poll(&mut cancelled)?;
        if balance.accepted_steps + balance.rejected_steps >= controls.maximum_steps {
            return Err(IsochoricKineticsError::WorkLimit {
                elapsed_time_s: elapsed,
                accepted_steps: balance.accepted_steps,
                rejected_steps: balance.rejected_steps,
            });
        }
        let remaining = duration.value() - elapsed;
        step = step.min(remaining).min(controls.maximum_step_s);
        if !positive(step) || elapsed + step == elapsed {
            return Err(IsochoricKineticsError::StepFailure {
                elapsed_time_s: elapsed,
                step_s: step,
                reason: "time increment is not representable",
            });
        }
        let attempt = (|| {
            let full_state = reactor.implicit_step(
                current,
                step,
                full,
                stage_forcing,
                accepted,
                &mut cancelled,
            )?;
            full_concentrations.copy_from_slice(reactor.concentrations());
            let half_state = reactor.implicit_step(
                current,
                step * 0.5,
                half,
                stage_forcing,
                accepted,
                &mut cancelled,
            )?;
            let fine_state = reactor.implicit_step(
                half,
                step * 0.5,
                fine,
                stage_forcing,
                half_state,
                &mut cancelled,
            )?;
            let mut error = ((fine_state.temperature.value() - full_state.temperature.value())
                / (controls.relative_tolerance
                    * fine_state
                        .temperature
                        .value()
                        .max(accepted.temperature.value())))
            .abs();
            for (i, &c) in reactor.concentrations().iter().enumerate() {
                let scale = controls.absolute_amount_tolerance_mol / volume.value()
                    + controls.relative_tolerance * c.abs().max(current_concentrations[i]);
                error = error.max((c - full_concentrations[i]).abs() / scale);
            }
            // The two-half-step endpoint has estimated local error
            // (fine - full)/(2^p - 1); no Richardson extrapolation is applied
            // to the conserved state or its nonnegative species inventories.
            error /= controls.method.error_divisor();
            let entropy_change = (fine_state.entropy_density_j_m3_k
                - accepted.entropy_density_j_m3_k)
                * volume.value();
            if !error.is_finite()
                || !entropy_change.is_finite()
                || entropy_change < -controls.entropy_tolerance_j_k
            {
                return Err(IsochoricKineticsError::StepFailure {
                    elapsed_time_s: elapsed,
                    step_s: step,
                    reason: "trial entropy decreases or error estimate is unrepresentable",
                });
            }
            Ok((fine_state, error, entropy_change))
        })();
        match attempt {
            Err(IsochoricKineticsError::Cancelled) => {
                return Err(IsochoricKineticsError::Cancelled);
            }
            Ok((state, error, entropy_change)) if error <= 1.0 => {
                audit_charge(
                    &reactor.mechanism,
                    |i| reactor.initial[i],
                    |i| reactor.concentrations()[i],
                    controls.relative_conservation_tolerance,
                    &mut balance,
                )?;
                current_concentrations.copy_from_slice(reactor.concentrations());
                reactor.rebase(current_concentrations);
                current.fill(0.0);
                accepted = state;
                elapsed = if step == remaining {
                    duration.value()
                } else {
                    elapsed + step
                };
                if balance.accepted_steps == 0 {
                    balance.minimum_step_entropy_change_j_k = entropy_change;
                } else {
                    balance.minimum_step_entropy_change_j_k =
                        balance.minimum_step_entropy_change_j_k.min(entropy_change);
                }
                balance.accepted_steps += 1;
                balance.maximum_estimated_local_error =
                    balance.maximum_estimated_local_error.max(error);
                let factor = if error == 0.0 {
                    2.0
                } else {
                    (0.85 / controls.method.error_root(error)).clamp(0.5, 2.0)
                };
                step = (step * factor).min(controls.maximum_step_s);
            }
            other => {
                balance.rejected_steps += 1;
                let factor = match other {
                    Ok((_, error, _)) => (0.8 / controls.method.error_root(error)).clamp(0.1, 0.5),
                    _ => 0.5,
                };
                step *= factor;
                if step < controls.minimum_step_s || elapsed + step == elapsed {
                    return Err(IsochoricKineticsError::StepFailure {
                        elapsed_time_s: elapsed,
                        step_s: step,
                        reason: "no admissible step above the declared minimum",
                    });
                }
            }
        }
    }
    // Resolve the exact accepted composition after the last coordinate rebase.
    // Each local update used stoichiometric extents without clipping; the final
    // mass/element/charge audit compares against the untouched original feed.
    accepted = reactor.resolve(current, accepted.temperature.value(), &mut cancelled)?;
    for (amount, concentration) in amounts.iter_mut().zip(reactor.concentrations()) {
        *amount = Amount::new(*concentration * volume.value());
    }
    balance.mass_residual_kg =
        (accepted.mass_density_kg_m3 - initial.mass_density_kg_m3) * volume.value();
    balance.energy_residual_j =
        accepted.internal_energy_density_j_m3 * volume.value() - target.value();
    balance.entropy_change_j_k =
        (accepted.entropy_density_j_m3_k - initial.entropy_density_j_m3_k) * volume.value();
    balance.nonlinear_iterations = reactor.stats.nonlinear_iterations;
    balance.temperature_evaluations = reactor.stats.temperature_evaluations;
    balance.direct_fallbacks = reactor.stats.direct_fallbacks;
    balance.maximum_relative_element_residual = element_residual(
        &reactor.mechanism,
        &reactor.initial,
        reactor.concentrations(),
    );
    audit_charge(
        &reactor.mechanism,
        |i| feed[i].value(),
        |i| amounts[i].value(),
        controls.relative_conservation_tolerance,
        &mut balance,
    )?;
    let final_mass = accepted.mass_density_kg_m3 * volume.value();
    let final_energy = accepted.internal_energy_density_j_m3 * volume.value();
    if !positive(final_mass)
        || !final_energy.is_finite()
        || !balance.mass_residual_kg.is_finite()
        || !balance.maximum_relative_element_residual.is_finite()
        || !balance.energy_residual_j.is_finite()
        || !balance.entropy_change_j_k.is_finite()
        || !balance.minimum_step_entropy_change_j_k.is_finite()
        || !balance.maximum_estimated_local_error.is_finite()
        || balance.maximum_relative_element_residual > controls.relative_conservation_tolerance
        || balance.mass_residual_kg.abs()
            > controls.relative_conservation_tolerance * initial.mass_density_kg_m3 * volume.value()
        || balance.energy_residual_j.abs() > controls.energy_tolerance_j
        || balance.entropy_change_j_k < -controls.entropy_tolerance_j_k
        || amounts
            .iter()
            .any(|n| !n.value().is_finite() || n.value() < 0.0)
    {
        return Err(IsochoricKineticsError::Balance(balance));
    }
    let result = IsochoricKineticsView {
        mechanism: &reactor.mechanism,
        amounts,
        temperature: accepted.temperature,
        pressure: accepted.pressure,
        internal_energy: Energy::new(final_energy),
        target_internal_energy: target,
        mass: Mass::new(final_mass),
        volume,
        elapsed_time: duration,
        balance,
    };
    poll(&mut cancelled)?;
    Ok(result)
}

fn feed_amounts(
    mechanism: &IdealGasKineticsV1,
    feed: &[(SpeciesId, Amount)],
) -> Result<Vec<Amount>, IsochoricKineticsError> {
    let species = mechanism.stoichiometry().species();
    if feed.is_empty() || feed.len() > species.len() {
        return Err(IsochoricKineticsError::Invalid("invalid kinetic feed size"));
    }
    let mut amounts = Vec::new();
    amounts
        .try_reserve_exact(species.len())
        .map_err(|_| IsochoricKineticsError::AllocationRefused)?;
    amounts.resize(species.len(), Amount::new(0.0));
    let mut present = vec![false; species.len()];
    for (id, amount) in feed {
        let i = species
            .binary_search(id)
            .map_err(|_| IsochoricKineticsError::Invalid("unknown kinetic feed species"))?;
        if present[i] || !amount.value().is_finite() || amount.value() < 0.0 {
            return Err(IsochoricKineticsError::Invalid(
                "duplicate species or invalid feed amount",
            ));
        }
        present[i] = true;
        amounts[i] = *amount;
    }
    Ok(amounts)
}
// Fixed-volume concentrations and amounts have the same relative charge
// residual. Scale both charges and charged concentrations before multiplying:
// no large neutral inventory, overflowing product, or underflowed absolute
// carrier total can erase a physically relevant charge discrepancy.
fn audit_charge(
    mechanism: &IdealGasKineticsV1,
    before: impl Fn(usize) -> f64,
    after: impl Fn(usize) -> f64,
    tolerance: f64,
    balance: &mut IsochoricKineticsBalance,
) -> Result<(), IsochoricKineticsError> {
    let charges = mechanism.thermodynamics().charge_vector();
    let mut concentration_scale = 0.0_f64;
    let mut charge_scale = 0.0_f64;
    for i in 0..charges.species().len() {
        let z = charges.get(i).expect("admitted charge shape") as f64;
        if z == 0.0 {
            continue;
        }
        let (initial, current) = (before(i), after(i));
        if !initial.is_finite() || initial < 0.0 || !current.is_finite() || current < 0.0 {
            balance.maximum_relative_charge_residual = f64::INFINITY;
            return Err(IsochoricKineticsError::Balance(*balance));
        }
        if initial > 0.0 || current > 0.0 {
            concentration_scale = concentration_scale.max(initial).max(current);
            charge_scale = charge_scale.max(z.abs());
        }
    }
    if concentration_scale == 0.0 {
        return Ok(());
    }
    let mut initial_carriers = Sum::default();
    let mut final_carriers = Sum::default();
    let mut change = Sum::default();
    for i in 0..charges.species().len() {
        let z = charges.get(i).expect("admitted charge shape") as f64 / charge_scale;
        if z == 0.0 {
            continue;
        }
        let (initial, current) = (before(i), after(i));
        initial_carriers.add(z.abs() * (initial / concentration_scale));
        final_carriers.add(z.abs() * (current / concentration_scale));
        change.add(z * ((current - initial) / concentration_scale));
    }
    let carrier_sum = initial_carriers.value().max(final_carriers.value());
    let residual = change.value().abs() / carrier_sum;
    if !residual.is_finite() || residual > tolerance {
        balance.maximum_relative_charge_residual = if residual.is_finite() {
            balance.maximum_relative_charge_residual.max(residual)
        } else {
            f64::INFINITY
        };
        return Err(IsochoricKineticsError::Balance(*balance));
    }
    balance.maximum_relative_charge_residual =
        balance.maximum_relative_charge_residual.max(residual);
    Ok(())
}

fn element_residual(mechanism: &IdealGasKineticsV1, before: &[f64], after: &[f64]) -> f64 {
    let matrix = mechanism.thermodynamics().elemental_matrix();
    let mut largest = 0.0_f64;
    for row in 0..matrix.row_count() {
        let mut initial = Sum::default();
        let mut change = Sum::default();
        for i in 0..before.len() {
            let a = matrix.get(row, i).expect("admitted element shape") as f64;
            initial.add(a * before[i]);
            change.add(a * (after[i] - before[i]));
        }
        if !initial.value().is_finite() || !change.value().is_finite() {
            return f64::INFINITY;
        }
        let residual = if initial.value() == 0.0 {
            if change.value() == 0.0 {
                0.0
            } else {
                f64::INFINITY
            }
        } else {
            change.value().abs() / initial.value()
        };
        if !residual.is_finite() {
            return f64::INFINITY;
        }
        largest = largest.max(residual);
    }
    largest
}
fn zeros(n: usize) -> Result<Vec<f64>, IsochoricKineticsError> {
    let mut v = Vec::new();
    v.try_reserve_exact(n)
        .map_err(|_| IsochoricKineticsError::AllocationRefused)?;
    v.resize(n, 0.0);
    Ok(v)
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), IsochoricKineticsError> {
    if cancelled() {
        Err(IsochoricKineticsError::Cancelled)
    } else {
        Ok(())
    }
}
#[derive(Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, x: f64) {
        let n = self.sum + x;
        self.correction += if self.sum.abs() >= x.abs() {
            (self.sum - n) + x
        } else {
            (x - n) + self.sum
        };
        self.sum = n;
    }
    fn value(&self) -> f64 {
        self.sum + self.correction
    }
}

#[cfg(test)]
mod tests;
