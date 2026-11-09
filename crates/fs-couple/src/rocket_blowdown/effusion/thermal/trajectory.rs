//! Adaptive physical windows of the original finite-wall/recoiling-gas step.
//!
//! Full versus two-half-step comparisons include the wall, gas composition,
//! vessel motion and emitted moments. Only accepted half steps contribute to
//! history and fixed whole-window budgets. Wall/gas heat is never counted as
//! extra elapsed time or extra vessel mass. This controls a local estimate of
//! first-order splitting error, not a certified global or constitutive error.

/// Molecular supply and ambient capture through the same adaptive wall owner.
pub mod exchange;

use super::{
    ConvectionCoefficient, EffusionOutlet, EffusionRocketEnvironment, EffusionRocketError,
    EffusionWallError, EffusionWallSolve, EffusionWallStep, GasRocketState,
    LumpedThermalEnvironment, ResolvedThermochemicalGasState, Sum, UniformEnthalpyStepInput,
};
use crate::rocket_blowdown::effusion::trajectory::{
    EffusionAdaptiveSolve, EffusionTrajectoryError, SpeciesAxis, admit_adaptive, normalized_error,
    ratio, trial_times,
};
use fs_exec::Cx;
use fs_material::{
    phase::EquilibriumPhaseState,
    thermochemical_gas::energy::parcel::{
        ClosedGasParcel,
        transfer::species::network::{
            SpeciesNetworkError, SpeciesNetworkStep, transport::GasTransportError,
        },
    },
};
use fs_mbd::Vec3;
use fs_qty::Temperature;
use fs_thermochem::{ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

/// Existing gas/motion control plus wall/export estimates and fixed budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffusionWallAdaptiveSolve {
    pub motion: EffusionAdaptiveSolve,
    pub wall_temperature_tolerance_k: f64,
    /// Wall storage error, independent of the enthalpy reference [J].
    pub wall_energy_tolerance_j: f64,
    /// Componentwise molecular impulse in the declared inertial frame [N s].
    /// Exchange histories apply this separately to incoming, outgoing and net.
    pub emitted_momentum_tolerance_n_s: f64,
    /// Emitted energy in that frame, with original species formation baselines
    /// removed before arithmetic. All trials use the same species baselines;
    /// exchange histories control incoming, outgoing and signed net separately.
    pub emitted_energy_tolerance_j: f64,
    pub window_mass_tolerance_kg: f64,
    pub window_momentum_tolerance_n_s: f64,
    pub window_energy_tolerance_j: f64,
}

/// Signed closure or a nonnegative componentwise debit, in coherent SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffusionWallWindowBalance {
    pub mass_kg: f64,
    pub momentum_n_s: Vec3,
    pub energy_j: f64,
}

/// Every refusal leaves the original wall, gas and vehicle unchanged.
#[derive(Debug)]
pub enum EffusionWallTrajectoryError<E = core::convert::Infallible> {
    /// Original temporal admission, minimum-step and work-budget refusals.
    Control(EffusionTrajectoryError),
    Step {
        elapsed_s: f64,
        duration_s: f64,
        error: EffusionWallError<E>,
    },
    Caloric(ThermochemErrorV1),
    Accounting(EffusionRocketError),
    WindowBalance {
        residual: EffusionWallWindowBalance,
        debit: EffusionWallWindowBalance,
        maximum_species_mass_residual_kg: f64,
    },
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for EffusionWallTrajectoryError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive finite-wall effusive recoil: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for EffusionWallTrajectoryError<E> {}
impl<E> From<EffusionTrajectoryError> for EffusionWallTrajectoryError<E> {
    fn from(error: EffusionTrajectoryError) -> Self {
        match error {
            EffusionTrajectoryError::Cancelled => Self::Cancelled,
            error => Self::Control(error),
        }
    }
}
impl<E> From<EffusionRocketError> for EffusionWallTrajectoryError<E> {
    fn from(error: EffusionRocketError) -> Self {
        Self::Accounting(error)
    }
}

/// Accepted native reports and compensated totals over exactly one window.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionWallTrajectoryStep<S = ()> {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub rejected_intervals: usize,
    pub accepted_intervals: usize,
    /// Exactly two original thermal/recoil steps per accepted interval.
    pub accepted_steps: Vec<EffusionWallStep<S>>,
    pub maximum_accepted_error: f64,
    pub escaped_mass_kg: f64,
    pub escaped_momentum_world_kg_m_s: Vec3,
    pub escaped_energy_world_j: f64,
    /// Same inertial frame, but zero species internal energies at window start T.
    /// This diagnostic removes arbitrary formation datums before evaluation.
    pub escaped_reference_energy_world_j: f64,
    pub external_impulse_n_s: Vec3,
    pub external_work_j: f64,
    /// Authored body/gas heat plus accepted enclosing radiation, never wall Q.
    pub external_thermal_heat_j: f64,
    pub wall_storage_change_j: f64,
    /// Actual allocated body/gas heat minus each authored window total [J].
    pub prescribed_heat_residual_j: [f64; 2],
    pub residual: EffusionWallWindowBalance,
    /// Sum of absolute step defects and native reserves, including initial
    /// inventory remapping and heat prorating. Signed cancellation buys no budget.
    pub debit: EffusionWallWindowBalance,
    pub maximum_species_mass_residual_kg: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EffusionWallTrajectoryCandidate<S = ()> {
    endpoint: GasRocketState,
    wall: EquilibriumPhaseState,
    step: EffusionWallTrajectoryStep<S>,
}
impl<S> EffusionWallTrajectoryCandidate<S> {
    #[must_use]
    pub const fn endpoint(&self) -> &GasRocketState {
        &self.endpoint
    }
    #[must_use]
    pub const fn wall(&self) -> &EquilibriumPhaseState {
        &self.wall
    }
    #[must_use]
    pub const fn step(&self) -> &EffusionWallTrajectoryStep<S> {
        &self.step
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        GasRocketState,
        EquilibriumPhaseState,
        EffusionWallTrajectoryStep<S>,
    ) {
        (self.endpoint, self.wall, self.step)
    }
}

fn check<E>(cx: &Cx<'_>) -> Result<(), EffusionWallTrajectoryError<E>> {
    cx.checkpoint()
        .map_err(|_| EffusionWallTrajectoryError::Cancelled)
}
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn vector(s: &[Sum; 3]) -> Vec3 {
    Vec3::new(s[0].total(), s[1].total(), s[2].total())
}
fn add_vector(s: &mut [Sum; 3], v: Vec3) -> Result<(), EffusionRocketError> {
    for (sum, value) in s.iter_mut().zip(components(v)) {
        sum.add(value)?;
    }
    Ok(())
}
fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64, EffusionRocketError> {
    let mut total = Sum::default();
    for value in values {
        total.add(value)?;
    }
    Ok(total.total())
}

// Physical energy in one common inertial frame, with u_i(T_window_start)
// canceled before arithmetic. Effusive molecules carry u_i(T) + R*T/(2*M_i).
// This uses accepted species masses and the original source's caloric integral,
// not a resampled flux or subtraction of two huge formation energies.
fn emitted_reference_energy<S, E>(
    cx: &Cx<'_>,
    step: &EffusionWallStep<S>,
    species: &SpeciesAxis<'_>,
    reference_t: f64,
) -> Result<f64, EffusionWallTrajectoryError<E>> {
    let temperature = step.thermal.report.next_gas.state().state().temperature;
    let mut energy = Sum::default();
    for row in &step.recoil.gas.exchange.boundaries {
        check(cx)?;
        let index = species
            .sources
            .binary_search_by(|s| s.species().cmp(&row.transfer.species))
            .map_err(|_| {
                EffusionTrajectoryError::Invalid("escaped species lost its original source")
            })?;
        let source = species.sources[index];
        let delta_u = source
            .molar_internal_energy_change(
                Temperature::new(reference_t),
                Temperature::new(temperature),
            )
            .map_err(EffusionWallTrajectoryError::Caloric)?
            .value();
        energy.add(
            row.transfer.outgoing_mass_kg
                * ((delta_u + 0.5 * R * temperature) / source.molar_mass().value()),
        )?;
    }
    energy.add(step.recoil.exhaust.relative_kinetic_transport_energy_j)?;
    energy.add(step.recoil.source_kinetic_correction_j)?;
    Ok(energy.total())
}

fn retry_withdrawal<E>(error: &EffusionWallTrajectoryError<E>) -> bool {
    matches!(
        error,
        EffusionWallTrajectoryError::Step {
            error: EffusionWallError::Effusion(EffusionRocketError::Gas(
                GasTransportError::Exchange(SpeciesNetworkError::Withdrawal { .. })
            )),
            ..
        }
    )
}

#[derive(Default)]
struct Totals {
    mass: Sum,
    momentum: [Sum; 3],
    energy: Sum,
    reference_energy: Sum,
    impulse: [Sum; 3],
    work: Sum,
    radiation: Sum,
    prescribed: [Sum; 2],
    mass_debit: Sum,
    momentum_residual: [Sum; 3],
    momentum_debit: [Sum; 3],
    energy_residual: Sum,
    energy_debit: Sum,
}

// Only native receipt projection varies between vacuum and reservoir histories.
// Time, wall/species error control and fixed-window accounting have one owner.
trait HistoryStep {
    fn exchange(&self) -> &SpeciesNetworkStep;
    fn momenta(&self) -> [Vec3; 3];
}
impl<S> HistoryStep for EffusionWallStep<S> {
    fn exchange(&self) -> &SpeciesNetworkStep {
        &self.recoil.gas.exchange
    }
    fn momenta(&self) -> [Vec3; 3] {
        let outgoing = self.recoil.exhaust_momentum_world_kg_m_s;
        [Vec3::ZERO, outgoing, outgoing.scale(-1.0)]
    }
}

struct BalanceReceipt<'a> {
    exchange: &'a SpeciesNetworkStep,
    initial_mass_residual_kg: f64,
    momentum_residual: Vec3,
    initial_momentum_correction: Vec3,
    gravity_impulse_uncertainty: Vec3,
    energy_residual_j: f64,
    absolute_energy_residual_sum_j: f64,
    energy_uncertainty_bound_j: f64,
}

struct History<T> {
    endpoint: GasRocketState,
    wall: EquilibriumPhaseState,
    duration_s: f64,
    attempted_intervals: usize,
    rejected_intervals: usize,
    accepted_intervals: usize,
    accepted_steps: Vec<T>,
    maximum_accepted_error: f64,
    totals: Totals,
    external_thermal_heat_j: f64,
    wall_storage_change_j: f64,
    prescribed_heat_residual_j: [f64; 2],
    residual: EffusionWallWindowBalance,
    debit: EffusionWallWindowBalance,
    maximum_species_mass_residual_kg: f64,
}

impl Totals {
    fn add<S>(
        &mut self,
        step: &EffusionWallStep<S>,
        energy: f64,
        heat: [f64; 2],
    ) -> Result<(), EffusionRocketError> {
        let recoil = &step.recoil;
        self.mass.add(recoil.expelled_mass_kg)?;
        add_vector(&mut self.momentum, recoil.exhaust_momentum_world_kg_m_s)?;
        self.energy.add(recoil.exhaust_energy_world_j)?;
        self.reference_energy.add(energy)?;
        add_vector(&mut self.impulse, recoil.exhaust.external_impulse_n_s)?;
        self.work.add(recoil.exhaust.external_work_j)?;
        self.radiation.add(step.thermal.report.radiation_heat_j)?;
        for (sum, value) in self.prescribed.iter_mut().zip(heat) {
            sum.add(value)?;
        }
        self.add_balance(BalanceReceipt {
            exchange: &recoil.gas.exchange,
            initial_mass_residual_kg: recoil.initial_mass_residual_kg,
            momentum_residual: recoil.momentum_residual_kg_m_s,
            initial_momentum_correction: recoil.initial_inventory_momentum_correction_n_s,
            gravity_impulse_uncertainty: recoil.gravity_impulse_uncertainty_n_s,
            energy_residual_j: step.combined_energy_residual_j,
            absolute_energy_residual_sum_j: step.absolute_energy_residual_sum_j,
            energy_uncertainty_bound_j: step.energy_uncertainty_bound_j,
        })
    }

    fn add_balance(&mut self, receipt: BalanceReceipt<'_>) -> Result<(), EffusionRocketError> {
        self.mass_debit.add(
            receipt
                .exchange
                .total_mass_residual_kg
                .abs()
                .max(receipt.exchange.maximum_species_mass_residual_kg),
        )?;
        self.mass_debit
            .add(receipt.initial_mass_residual_kg.abs())?;
        add_vector(&mut self.momentum_residual, receipt.momentum_residual)?;
        for ((total, residual), (initial, gravity)) in self
            .momentum_debit
            .iter_mut()
            .zip(components(receipt.momentum_residual))
            .zip(
                components(receipt.initial_momentum_correction)
                    .into_iter()
                    .zip(components(receipt.gravity_impulse_uncertainty)),
            )
        {
            total.add(residual.abs())?;
            total.add(initial.abs())?;
            total.add(gravity)?;
        }
        self.energy_residual.add(receipt.energy_residual_j)?;
        self.energy_debit
            .add(receipt.absolute_energy_residual_sum_j)?;
        self.energy_debit.add(receipt.energy_uncertainty_bound_j)?;
        Ok(())
    }
}

impl GasRocketState {
    /// Evolve finite wall heat replenishment, effusive cooling and actual recoil.
    ///
    /// `wall.duration_s` is the entire window. `wall.internal_heat_j` and
    /// `gas_heat_j` are whole-window totals interpreted as held powers, prorated
    /// for every physical trial. Actual accepted sums are compared with those
    /// authored totals and their rounding consumes the fixed energy budget.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_with_wall_adaptive(
        &self,
        cx: &Cx<'_>,
        environment: EffusionRocketEnvironment,
        outlets: &[EffusionOutlet<'_>],
        wall: UniformEnthalpyStepInput<'_>,
        body: &LumpedThermalEnvironment,
        gas_heat_j: f64,
        solve: EffusionWallSolve,
        adaptive: EffusionWallAdaptiveSolve,
    ) -> Result<EffusionWallTrajectoryCandidate, EffusionWallTrajectoryError> {
        self.advance_effusion_with_wall_convection_adaptive(
            cx,
            environment,
            outlets,
            wall,
            body,
            gas_heat_j,
            solve,
            adaptive,
            |_| {
                Ok(ConvectionCoefficient {
                    coefficient_w_per_m2_k: body.convection_w_per_m2_k,
                    source: (),
                })
            },
        )
    }

    /// Resolve the repeatable original convection callback at every thermal trial.
    ///
    /// Coarse and fine steps all pass original physical/source gates. Only a
    /// typed gross-withdrawal refusal may trigger physical retry. A temporal
    /// estimate above one halves; an accepted estimate at most one quarter may
    /// double. A failed shorter final remainder refuses instead of skipping time.
    /// Neither rejected states nor reports enter accepted totals; callbacks may
    /// still observe rejected evaluations, and their side effects are not undone.
    ///
    /// Temperature/species and drift-free motion use the existing recoil error
    /// owner. Wall energy uses increments. Emitted moments share the declared
    /// inertial frame; emitted energy additionally cancels fixed species
    /// formation baselines before evaluation. Inertial export tolerances are
    /// frame-specific, not asserted to be invariant under a change of observer.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_with_wall_convection_adaptive<S, E>(
        &self,
        cx: &Cx<'_>,
        environment: EffusionRocketEnvironment,
        outlets: &[EffusionOutlet<'_>],
        wall: UniformEnthalpyStepInput<'_>,
        body: &LumpedThermalEnvironment,
        gas_heat_j: f64,
        solve: EffusionWallSolve,
        adaptive: EffusionWallAdaptiveSolve,
        mut resolve: impl FnMut(&ResolvedThermochemicalGasState) -> Result<ConvectionCoefficient<S>, E>,
    ) -> Result<EffusionWallTrajectoryCandidate<S>, EffusionWallTrajectoryError<E>> {
        let history = integrate(
            self,
            std::iter::empty(),
            cx,
            wall,
            gas_heat_j,
            adaptive,
            |state, input, heat| {
                state
                    .advance_effusion_with_wall_convection(
                        cx,
                        environment,
                        outlets,
                        input,
                        body,
                        heat,
                        solve,
                        &mut resolve,
                    )
                    .map(|candidate| candidate.into_parts())
            },
            |step, species| {
                let energy = emitted_reference_energy(
                    cx,
                    step,
                    species,
                    self.reservoir().state().state().temperature,
                )?;
                Ok([0.0, energy, -energy])
            },
            |totals, step, energies, heat| totals.add(step, energies[1], heat),
        )?;
        let totals = history.totals;
        Ok(EffusionWallTrajectoryCandidate {
            endpoint: history.endpoint,
            wall: history.wall,
            step: EffusionWallTrajectoryStep {
                duration_s: history.duration_s,
                attempted_intervals: history.attempted_intervals,
                rejected_intervals: history.rejected_intervals,
                accepted_intervals: history.accepted_intervals,
                accepted_steps: history.accepted_steps,
                maximum_accepted_error: history.maximum_accepted_error,
                escaped_mass_kg: totals.mass.total(),
                escaped_momentum_world_kg_m_s: vector(&totals.momentum),
                escaped_energy_world_j: totals.energy.total(),
                escaped_reference_energy_world_j: totals.reference_energy.total(),
                external_impulse_n_s: vector(&totals.impulse),
                external_work_j: totals.work.total(),
                external_thermal_heat_j: history.external_thermal_heat_j,
                wall_storage_change_j: history.wall_storage_change_j,
                prescribed_heat_residual_j: history.prescribed_heat_residual_j,
                residual: history.residual,
                debit: history.debit,
                maximum_species_mass_residual_kg: history.maximum_species_mass_residual_kg,
            },
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn integrate<'a, T: HistoryStep, E>(
    initial: &'a GasRocketState,
    sources: impl IntoIterator<Item = &'a ClosedGasParcel>,
    cx: &Cx<'_>,
    wall: UniformEnthalpyStepInput<'_>,
    gas_heat_j: f64,
    adaptive: EffusionWallAdaptiveSolve,
    mut perform: impl FnMut(
        &GasRocketState,
        UniformEnthalpyStepInput<'_>,
        f64,
    )
        -> Result<(GasRocketState, EquilibriumPhaseState, T), EffusionWallError<E>>,
    mut energies: impl FnMut(&T, &SpeciesAxis<'_>) -> Result<[f64; 3], EffusionWallTrajectoryError<E>>,
    mut account: impl FnMut(&mut Totals, &T, [f64; 3], [f64; 2]) -> Result<(), EffusionRocketError>,
) -> Result<History<T>, EffusionWallTrajectoryError<E>> {
    check(cx)?;
    let duration = wall.duration_s;
    let motion = adaptive.motion;
    admit_adaptive(duration, motion)?;
    if ![
        adaptive.wall_temperature_tolerance_k,
        adaptive.wall_energy_tolerance_j,
        adaptive.emitted_momentum_tolerance_n_s,
        adaptive.emitted_energy_tolerance_j,
    ]
    .into_iter()
    .all(|v| v.is_finite() && v > 0.0)
        || ![
            adaptive.window_mass_tolerance_kg,
            adaptive.window_momentum_tolerance_n_s,
            adaptive.window_energy_tolerance_j,
        ]
        .into_iter()
        .all(|v| v.is_finite() && v >= 0.0)
        || !wall.internal_heat_j.is_finite()
        || !gas_heat_j.is_finite()
    {
        return Err(EffusionTrajectoryError::Invalid(
                "finite heat, positive wall/export error scales and nonnegative fixed window budgets required",
            ).into());
    }
    let mut cancelled = || cx.checkpoint().is_err();
    let species = SpeciesAxis::new(
        std::iter::once(initial.reservoir()).chain(sources),
        &mut cancelled,
    )?;
    let initial_masses = species.masses(initial, &mut cancelled)?;
    let mut escaped_species: Vec<Sum> =
        (0..species.sources.len()).map(|_| Sum::default()).collect();
    let mut current = initial.clone();
    let mut current_wall = wall.initial;
    let (mut elapsed, mut next_step, mut maximum_error) = (0.0, motion.initial_step_s, 0.0_f64);
    let (mut attempts, mut rejected, mut accepted) = (0, 0, 0);
    let mut accepted_steps = Vec::new();
    let mut totals = Totals::default();
    while elapsed < duration {
        check(cx)?;
        if attempts == motion.maximum_attempts || accepted == motion.maximum_accepted_intervals {
            return Err(EffusionTrajectoryError::WorkBudget {
                attempted_intervals: attempts,
                accepted_intervals: accepted,
            }
            .into());
        }
        let (h, half, midpoint, end) =
            trial_times(elapsed, next_step, duration, motion.minimum_step_s)?;
        attempts += 1;
        let heat = |dt: f64| {
            [
                wall.internal_heat_j * (dt / duration),
                gas_heat_j * (dt / duration),
            ]
        };
        let trial = (|| {
            let mut trial_step = |state: &GasRocketState, initial, dt, at| {
                let q = heat(dt);
                let input = UniformEnthalpyStepInput {
                    initial,
                    internal_heat_j: q[0],
                    duration_s: dt,
                    ..wall
                };
                perform(state, input, q[1]).map_err(|error| {
                    if cx.checkpoint().is_err() {
                        EffusionWallTrajectoryError::Cancelled
                    } else {
                        EffusionWallTrajectoryError::Step {
                            elapsed_s: at,
                            duration_s: dt,
                            error,
                        }
                    }
                })
            };
            let (coarse, coarse_wall, coarse_step) =
                trial_step(&current, current_wall, h, elapsed)?;
            let (first, first_wall, first_step) =
                trial_step(&current, current_wall, half, elapsed)?;
            let (fine, fine_wall, second_step) = trial_step(&first, first_wall, half, midpoint)?;
            let mut error = normalized_error(
                &current,
                &coarse,
                &fine,
                &species,
                h,
                motion,
                &mut cancelled,
            )?;
            let coarse_storage = wall.mass_kg
                * (coarse_wall.specific_enthalpy_j_kg() - current_wall.specific_enthalpy_j_kg());
            let fine_storage = wall.mass_kg
                * (fine_wall.specific_enthalpy_j_kg() - current_wall.specific_enthalpy_j_kg());
            // Use the endpoint difference before multiplying: separately
            // rounded storage products must not erase a small discrepancy.
            // The original thermal owner still refuses unresolved changes
            // in a wall chart whose enthalpy datum exhausts f64 resolution.
            let wall_energy_error = wall.mass_kg
                * (fine_wall.specific_enthalpy_j_kg() - coarse_wall.specific_enthalpy_j_kg());
            error = error.max(ratio(
                0.0,
                wall_energy_error,
                coarse_storage.abs().max(fine_storage.abs()),
                adaptive.wall_energy_tolerance_j,
                motion.relative_tolerance,
            )?);
            error = error.max(ratio(
                coarse_wall.temperature_k(),
                fine_wall.temperature_k(),
                current_wall
                    .temperature_k()
                    .max(coarse_wall.temperature_k())
                    .max(fine_wall.temperature_k()),
                adaptive.wall_temperature_tolerance_k,
                motion.relative_tolerance,
            )?);
            for (full, (a, b)) in coarse_step
                .momenta()
                .into_iter()
                .zip(first_step.momenta().into_iter().zip(second_step.momenta()))
            {
                for (full, (a, b)) in components(full)
                    .into_iter()
                    .zip(components(a).into_iter().zip(components(b)))
                {
                    let fine = sum([a, b])?;
                    error = error.max(ratio(
                        full,
                        fine,
                        full.abs().max(fine.abs()),
                        adaptive.emitted_momentum_tolerance_n_s,
                        motion.relative_tolerance,
                    )?);
                }
            }
            let full_energies = energies(&coarse_step, &species)?;
            let half_energies = [
                energies(&first_step, &species)?,
                energies(&second_step, &species)?,
            ];
            for (i, full) in full_energies.into_iter().enumerate() {
                let fine = sum([half_energies[0][i], half_energies[1][i]])?;
                error = error.max(ratio(
                    full,
                    fine,
                    full.abs().max(fine.abs()),
                    adaptive.emitted_energy_tolerance_j,
                    motion.relative_tolerance,
                )?);
            }
            Ok::<_, EffusionWallTrajectoryError<E>>((
                fine,
                fine_wall,
                [first_step, second_step],
                half_energies,
                error,
            ))
        })();
        let (fine, next_wall, steps, energies, error) = match trial {
            Ok(result) => result,
            Err(error) if retry_withdrawal(&error) && half >= motion.minimum_step_s => {
                rejected += 1;
                next_step = half;
                continue;
            }
            Err(error) => return Err(error),
        };
        if error > 1.0 {
            rejected += 1;
            if half < motion.minimum_step_s {
                return Err(EffusionTrajectoryError::MinimumStep {
                    duration_s: h,
                    normalized_error: error,
                }
                .into());
            }
            next_step = half;
            continue;
        }
        for (step, energy) in steps.into_iter().zip(energies) {
            check(cx)?;
            account(&mut totals, &step, energy, heat(half))?;
            for row in &step.exchange().boundaries {
                check(cx)?;
                let index = species
                    .sources
                    .binary_search_by(|s| s.species().cmp(&row.transfer.species))
                    .map_err(|_| {
                        EffusionTrajectoryError::Invalid("accepted escape changed species axis")
                    })?;
                escaped_species[index].add(-row.transfer.net_mass_into_gas_kg)?;
            }
            accepted_steps.push(step);
        }
        current = fine;
        current_wall = next_wall;
        elapsed = end;
        accepted += 1;
        maximum_error = maximum_error.max(error);
        next_step = if error <= 0.25 {
            (2.0 * h).min(motion.maximum_step_s)
        } else {
            h
        };
    }
    let final_masses = species.masses(&current, &mut cancelled)?;
    let mut maximum_species_mass_residual_kg = 0.0_f64;
    for (i, escaped) in escaped_species.iter().enumerate() {
        maximum_species_mass_residual_kg = maximum_species_mass_residual_kg
            .max(sum([final_masses[i], -initial_masses[i], escaped.total()])?.abs());
    }
    let prescribed_heat_residual_j = [
        sum([totals.prescribed[0].total(), -wall.internal_heat_j])?,
        sum([totals.prescribed[1].total(), -gas_heat_j])?,
    ];
    let residual = EffusionWallWindowBalance {
        mass_kg: sum([
            current.reservoir().mass_kg(),
            -initial.reservoir().mass_kg(),
            totals.mass.total(),
        ])?,
        momentum_n_s: vector(&totals.momentum_residual),
        energy_j: sum([
            totals.energy_residual.total(),
            prescribed_heat_residual_j[0],
            prescribed_heat_residual_j[1],
        ])?,
    };
    let debit = EffusionWallWindowBalance {
        mass_kg: totals
            .mass_debit
            .total()
            .max(residual.mass_kg.abs())
            .max(maximum_species_mass_residual_kg),
        momentum_n_s: vector(&totals.momentum_debit),
        energy_j: sum([
            totals.energy_debit.total(),
            prescribed_heat_residual_j[0].abs(),
            prescribed_heat_residual_j[1].abs(),
        ])?,
    };
    if debit.mass_kg > adaptive.window_mass_tolerance_kg
        || components(debit.momentum_n_s)
            .into_iter()
            .any(|v| v > adaptive.window_momentum_tolerance_n_s)
        || debit.energy_j > adaptive.window_energy_tolerance_j
    {
        return Err(EffusionWallTrajectoryError::WindowBalance {
            residual,
            debit,
            maximum_species_mass_residual_kg,
        });
    }
    let external_thermal_heat_j =
        sum([wall.internal_heat_j, gas_heat_j, totals.radiation.total()])?;
    let wall_storage_change_j = wall.mass_kg
        * (current_wall.specific_enthalpy_j_kg() - wall.initial.specific_enthalpy_j_kg());
    if !wall_storage_change_j.is_finite() {
        return Err(EffusionTrajectoryError::Unrepresentable("whole-window wall storage").into());
    }
    check(cx)?;
    Ok(History {
        endpoint: current,
        wall: current_wall,
        duration_s: elapsed,
        attempted_intervals: attempts,
        rejected_intervals: rejected,
        accepted_intervals: accepted,
        accepted_steps,
        maximum_accepted_error: maximum_error,
        totals,
        external_thermal_heat_j,
        wall_storage_change_j,
        prescribed_heat_residual_j,
        residual,
        debit,
        maximum_species_mass_residual_kg,
    })
}
