//! One simultaneous original-wall, gas and film solve in energy coordinates.
//!
//! No wall or gas caloric root runs inside a numerical residual. After Newton
//! convergence each original owner reconstructs once, then the phase and heat
//! ports are evaluated against BOTH actual returned states before acceptance.

use super::*;
use fs_conduction::lumped::endpoint::{LumpedEnthalpyEndpoint, LumpedEnthalpyEndpointSample};
use fs_material::evaporation::gas::wall_film::endpoint::{
    GasWallFilmEndpoint, GasWallFilmEndpointTrial,
};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use std::cell::{Cell, RefCell};

/// Numerical work only; original wall, phase, source and balance controls remain
/// on `WallFilmThermalEnvironment::solve` and its original thermal body.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallFilmImplicitControls {
    pub newton: NewtonKrylovConfig,
    pub difference_step_j: f64,
    pub mass_coordinate_scale_j_kg: f64,
    pub maximum_iterations: usize,
    /// Includes the initial predictor, fresh convergence check, gas inversion's
    /// phase recheck and the final returned-wall/returned-gas phase recheck.
    pub maximum_evaluations: usize,
}
impl Default for WallFilmImplicitControls {
    fn default() -> Self {
        Self {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-11,
                relative_tolerance: 0.0,
                ..NewtonKrylovConfig::default()
            },
            difference_step_j: 1e-3,
            mass_coordinate_scale_j_kg: 1e6,
            maximum_iterations: 64,
            maximum_evaluations: 10_000,
        }
    }
}

#[derive(Debug)]
pub enum WallFilmImplicitError {
    Physical(WallFilmThermalError),
    Numerical(NewtonError),
    NonConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    EvaluationBudget,
    IterationBudget,
    GasIterationBudget,
}
impl core::fmt::Display for WallFilmImplicitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "implicit finite wall-film exchange: {self:?}")
    }
}
impl std::error::Error for WallFilmImplicitError {}
impl From<WallFilmThermalError> for WallFilmImplicitError {
    fn from(value: WallFilmThermalError) -> Self {
        Self::Physical(value)
    }
}
impl From<GasEvaporationError> for WallFilmImplicitError {
    fn from(value: GasEvaporationError) -> Self {
        Self::Physical(value.into())
    }
}
impl From<ConductionError> for WallFilmImplicitError {
    fn from(value: ConductionError) -> Self {
        Self::Physical(super::super::Error::Interface(value).into())
    }
}
impl From<NewtonError> for WallFilmImplicitError {
    fn from(value: NewtonError) -> Self {
        Self::Numerical(value)
    }
}
type Result<T> = core::result::Result<T, WallFilmImplicitError>;

/// Shared numerical allowance for complete physical proposals. Failed trials
/// keep their evaluation/Newton charges. A gas reconstruction reserves its
/// bounded loop only when called; success refunds unused iterations. A failed
/// finalizer exposes no partial receipt, so its full reservation stays charged.
pub(super) struct WorkMeter {
    evaluations: Cell<usize>,
    iterations: Cell<usize>,
    gas_iterations: Cell<usize>,
    gas_charged: Cell<usize>,
    maximum_evaluations: usize,
    maximum_iterations: usize,
    maximum_gas_iterations: usize,
}
impl WorkMeter {
    pub(super) fn new(
        controls: WallFilmImplicitControls,
        maximum_gas_iterations: usize,
    ) -> Result<Self> {
        if controls.maximum_evaluations == 0
            || controls.maximum_iterations == 0
            || maximum_gas_iterations == 0
        {
            return Err(Error::Invalid("positive implicit work allowances required").into());
        }
        Ok(Self {
            evaluations: Cell::new(0),
            iterations: Cell::new(0),
            gas_iterations: Cell::new(0),
            gas_charged: Cell::new(0),
            maximum_evaluations: controls.maximum_evaluations,
            maximum_iterations: controls.maximum_iterations,
            maximum_gas_iterations,
        })
    }
    pub(super) fn newton_iterations(&self) -> usize {
        self.iterations.get()
    }
    pub(super) fn residual_evaluations(&self) -> usize {
        self.evaluations.get()
    }
    pub(super) fn gas_iterations(&self) -> usize {
        self.gas_iterations.get()
    }
    pub(super) fn charged_gas_iterations(&self) -> usize {
        self.gas_charged.get()
    }
    fn remaining_iterations(&self) -> Result<usize> {
        let remaining = self.maximum_iterations - self.iterations.get();
        if remaining == 0 {
            Err(WallFilmImplicitError::IterationBudget)
        } else {
            Ok(remaining)
        }
    }
    fn gas_limit(&self, original_limit: usize) -> Result<usize> {
        let remaining = self.maximum_gas_iterations - self.gas_charged.get();
        if remaining == 0 {
            Err(WallFilmImplicitError::GasIterationBudget)
        } else {
            Ok(original_limit.min(remaining))
        }
    }
    fn refund_gas(&self, reserved: usize, actual: usize) -> Result<()> {
        if actual > reserved {
            return Err(Error::Invalid("original gas work exceeded its reserved loop").into());
        }
        self.gas_charged
            .set(self.gas_charged.get() - (reserved - actual));
        self.gas_iterations.set(self.gas_iterations.get() + actual);
        Ok(())
    }
}

/// Actual original endpoints and source-relative physical receipts. Latent,
/// vapor sensible, liquid sensible and nonadvective sensible heat stay separate.
#[derive(Clone, Debug, PartialEq)]
pub struct WallFilmImplicitReport {
    pub next_gas: ClosedGasParcel,
    pub next_film: GasWallFilmInventory,
    pub phase: GasWallFilmStep,
    pub body: LumpedEnthalpyMarch,
    pub solid_to_gas_heat_j: f64,
    pub wall_heat_into_phase_system_j: f64,
    pub external_gas_heat_j: f64,
    pub radiation_heat_j: f64,
    pub wall_energy_residual_j: f64,
    pub gas_energy_residual_j: f64,
    pub combined_energy_residual_j: f64,
    /// Fresh wall/film/gas defects, source-combination and returned-heat
    /// roundoff, plus the original gas inversion defect independently.
    pub sum_absolute_energy_residual_j: f64,
    pub mass_residual_kg: f64,
    pub maximum_mass_residual_kg: f64,
    /// Fresh phase inventory defect plus original gas reconstruction defects.
    pub sum_absolute_mass_residual_kg: f64,
    pub wall_trial_energy_residual_j: f64,
    pub wall_temperature_residual_k: f64,
    pub source_residual_j: f64,
    pub gas_source_combination_residual_j: f64,
    pub wall_source_combination_residual_j: f64,
    pub external_heat_roundoff_j: f64,
    pub newton_iterations: usize,
    pub residual_evaluations: usize,
    pub gas_iterations: usize,
    pub maximum_accepted_residual_norm_j: f64,
}

struct Prepared<'e, 'a, 'w> {
    environment: &'e WallFilmThermalEnvironment<'a>,
    input: UniformEnthalpyStepInput<'w>,
    wall: LumpedEnthalpyEndpoint<'w>,
    phase: GasWallFilmEndpoint<'a, 'a>,
    gas_capacity_j_k: f64,
    initial_vapor_mass_kg: f64,
    mass_scale_j_kg: f64,
}
struct Trial {
    phase: GasWallFilmEndpointTrial,
    wall: LumpedEnthalpyEndpointSample,
    q: f64,
}
impl Trial {
    fn residual(&self, mass_scale: f64) -> Result<Vec<f64>> {
        Ok(vec![
            self.wall.residual_j,
            finite(self.phase.gas_mass_residual_kg() * mass_scale)?,
            self.phase.gas_energy_residual_j(),
        ])
    }
}

impl<'e, 'a, 'w> Prepared<'e, 'a, 'w> {
    fn new(
        environment: &'e WallFilmThermalEnvironment<'a>,
        input: UniformEnthalpyStepInput<'w>,
        controls: WallFilmImplicitControls,
    ) -> Result<Self> {
        controls.newton.validate()?;
        if controls.maximum_iterations == 0
            || controls.maximum_evaluations == 0
            || !positive(controls.difference_step_j)
            || !positive(controls.mass_coordinate_scale_j_kg)
            || !environment.gas_heat_j.is_finite()
            || !environment.boundary.area_m2.is_finite()
            || environment.boundary.area_m2 < 0.0
            || environment.boundary.area_m2 > input.surface_area_m2
            || !positive(environment.solve.temperature_tolerance_k)
            || !positive(environment.solve.wall.source_tolerance_j)
            || !positive(environment.solve.wall.energy_tolerance_j)
            || !positive(environment.solve.sensible.source_tolerance_j)
            || !positive(environment.solve.sensible.energy_tolerance_j)
        {
            return Err(Error::Invalid(
                "invalid implicit wall-film source, geometry or work controls",
            )
            .into());
        }
        let state = environment.gas.state().state();
        let gas_capacity_j_k =
            finite(environment.gas.mass_kg() * state.specific_heat_cp / state.gamma)?;
        if !positive(gas_capacity_j_k) {
            return Err(Error::Invalid("invalid original gas caloric scaling").into());
        }
        let mut body = environment.body.clone();
        body.temperature_k = state.temperature;
        let wall = body.prepare_endpoint(input)?;
        let exposure =
            finite(body.convection_w_per_m2_k * input.surface_area_m2 * input.duration_s)?;
        if body.convection_w_per_m2_k > 0.0 && exposure == 0.0 {
            return Err(Error::Invalid("nonzero wall convection exposure underflows").into());
        }
        let liquid_h = liquid_h(environment, input.initial.temperature_k())?;
        let phase = GasWallFilmEndpoint::new(
            GasWallFilmSurface {
                law: environment.boundary.law,
                binding: environment.boundary.binding,
                film: environment.film,
                area_m2: environment.boundary.area_m2,
                interface_liquid_enthalpy_j_kg: liquid_h,
                temperature_from_enthalpy: environment.boundary.liquid_temperature_from_enthalpy,
                mass_transfer_coefficient_m_s: environment.boundary.mass_transfer_coefficient_m_s,
            },
            environment.gas,
            environment.catalog,
            input.duration_s,
            environment.solve.phase,
        )?;
        let initial_vapor_mass_kg = environment
            .gas
            .species_mass_kg(&environment.boundary.binding.species)
            .map_err(GasEvaporationError::from)?;
        Ok(Self {
            environment,
            input,
            wall,
            phase,
            gas_capacity_j_k,
            initial_vapor_mass_kg,
            mass_scale_j_kg: controls.mass_coordinate_scale_j_kg,
        })
    }
    fn project(&self, coordinates: &[f64]) -> Vec<f64> {
        let origins = [
            self.input.initial.specific_enthalpy_j_kg(),
            self.initial_vapor_mass_kg,
            self.environment.gas.state().state().temperature,
        ];
        let scales = [
            self.input.mass_kg,
            self.mass_scale_j_kg,
            self.gas_capacity_j_k,
        ];
        coordinates
            .iter()
            .zip(origins)
            .zip(scales)
            .map(|((&x, origin), scale)| {
                if x.is_finite() && origin + x / scale == origin {
                    0.0
                } else {
                    x
                }
            })
            .collect()
    }
    fn evaluate(&self, coordinates: &[f64], cx: &Cx<'_>) -> Result<Trial> {
        poll(cx)?;
        let gas_t = finite(
            self.environment.gas.state().state().temperature
                + coordinates[2] / self.gas_capacity_j_k,
        )?;
        let wall = self
            .wall
            .evaluate(coordinates[0], gas_t, self.input.internal_heat_j)?;
        let h = liquid_h(self.environment, wall.temperature_k)?;
        let q = -wall.convection_heat_j;
        let (gas_heat, _) = combine(self.environment.gas_heat_j, q)?;
        let phase = self.phase.evaluate(
            h,
            coordinates[1] / self.mass_scale_j_kg,
            gas_t,
            gas_heat,
            || cx.checkpoint().is_err(),
        )?;
        let (wall_heat, _) = combine(
            self.input.internal_heat_j,
            -phase.step().wall_heat_into_system_j,
        )?;
        let wall = self.wall.evaluate(coordinates[0], gas_t, wall_heat)?;
        Ok(Trial { phase, wall, q })
    }
}

impl WallFilmThermalEnvironment<'_> {
    /// Advance the full requested interval with endpoint wall, film and gas
    /// feedback in ONE original Newton--Krylov solve. Wall enthalpy, including
    /// native latent plateaus, is independent of the original liquid chart.
    ///
    /// The original Stefan rate sees final composition, pressure and gas/wall
    /// temperatures. A dry patch can condense without a seed; evaporation stops
    /// at exactly zero film. Reported dry-out time freezes this endpoint rate
    /// and is an estimate requiring physical timestep refinement.
    ///
    /// Source, donor, Biot, cancellation and final balance failures propagate.
    /// Legacy scalar brackets/iteration counts are unused; their source/energy
    /// tolerances and all original phase/body controls remain enforced. This
    /// fully-wetting closure introduces no boiling, nucleation or moving film
    /// geometry. Phase heat remains homogenized over the original uniform wall.
    pub fn advance_implicit(
        &self,
        cx: &Cx<'_>,
        input: UniformEnthalpyStepInput<'_>,
        controls: WallFilmImplicitControls,
    ) -> Result<UniformEnthalpyStep<WallFilmImplicitReport>> {
        let meter = WorkMeter::new(controls, self.solve.phase.gas.gas.maximum_iterations)?;
        self.advance_implicit_metered(cx, input, controls, &meter)
    }

    pub(super) fn advance_implicit_metered(
        &self,
        cx: &Cx<'_>,
        input: UniformEnthalpyStepInput<'_>,
        controls: WallFilmImplicitControls,
        meter: &WorkMeter,
    ) -> Result<UniformEnthalpyStep<WallFilmImplicitReport>> {
        poll(cx)?;
        let remaining_iterations = meter
            .remaining_iterations()?
            .min(controls.maximum_iterations);
        let gas_limit = meter.gas_limit(self.solve.phase.gas.gas.maximum_iterations)?;
        let mut bounded = self.clone();
        bounded.solve.phase.gas.gas.maximum_iterations = gas_limit;
        let prepared = Prepared::new(&bounded, input, controls)?;
        let evaluations = &meter.evaluations;
        let initial_evaluations = evaluations.get();
        spend(evaluations, meter.maximum_evaluations)?;
        let zero = prepared.project(&[0.0; 3]);
        let predictor = prepared.evaluate(&zero, cx)?;
        let guess = prepared.project(&[
            0.0,
            finite(
                predictor.phase.step().transferred_mass_kg * controls.mass_coordinate_scale_j_kg,
            )?,
            0.0,
        ]);
        let cache = RefCell::new(vec![(
            zero,
            predictor.residual(controls.mass_coordinate_scale_j_kg)?,
        )]);
        drop(predictor);
        let problem = FallibleFiniteDifference::new(
            3,
            controls.difference_step_j,
            |coordinates: &[f64]| {
                poll(cx)?;
                let point = prepared.project(coordinates);
                {
                    let mut entries = cache.borrow_mut();
                    if let Some(index) = entries.iter().position(|(x, _)| {
                        x.iter()
                            .zip(&point)
                            .all(|(a, b)| a.to_bits() == b.to_bits())
                    }) {
                        let entry = entries.remove(index);
                        let residual = entry.1.clone();
                        entries.push(entry);
                        return Ok(residual);
                    }
                }
                spend(evaluations, meter.maximum_evaluations)?;
                let residual = prepared
                    .evaluate(&point, cx)?
                    .residual(controls.mass_coordinate_scale_j_kg)?;
                let mut entries = cache.borrow_mut();
                if entries.len() == 2 {
                    entries.remove(0);
                }
                entries.push((point, residual.clone()));
                Ok::<_, WallFilmImplicitError>(residual)
            },
        )
        .map_err(difference_error)?;
        let initialization = NewtonKrylovState::new(&problem, guess, controls.newton);
        if let Some(error) = problem.take_failure() {
            return Err(difference_error(error));
        }
        let mut solver = initialization?;
        let mut cancelled = false;
        let run = solver.run_controlled(&problem, remaining_iterations, &mut || {
            if problem.has_failed() {
                return true;
            }
            cancelled = cx.checkpoint().is_err();
            cancelled
        });
        meter
            .iterations
            .set(meter.iterations.get() + run.report.iterations);
        if let Some(error) = problem.take_failure() {
            return Err(difference_error(error));
        }
        if cancelled {
            return Err(Error::Cancelled.into());
        }
        if run.progress != SolverRunProgress::Complete {
            return Err(WallFilmImplicitError::Interrupted(run.progress));
        }
        if !run.report.converged {
            return Err(WallFilmImplicitError::NonConverged(run.report));
        }
        spend(evaluations, meter.maximum_evaluations)?;
        let trial = prepared.evaluate(&prepared.project(&solver.x), cx)?;
        prepared.phase.check_balance(&trial.phase)?;
        gate(trial.wall.residual_j, self.solve.wall.source_tolerance_j)?;
        // The material finalizer performs exactly one fresh physical recheck.
        spend(evaluations, meter.maximum_evaluations)?;
        meter.gas_charged.set(meter.gas_charged.get() + gas_limit);
        let candidate = prepared
            .phase
            .finalize(&trial.phase, || cx.checkpoint().is_err())?;
        meter.refund_gas(gas_limit, candidate.gas_exchange().step().iterations)?;
        let mut body_environment = self.body.clone();
        let gas_t = candidate.gas().state().state().temperature;
        body_environment.temperature_k = gas_t;
        let (internal_heat_j, _) = combine(
            input.internal_heat_j,
            -candidate.step().wall_heat_into_system_j,
        )?;
        let body = body_environment.advance(
            cx,
            UniformEnthalpyStepInput {
                internal_heat_j,
                ..input
            },
        )?;
        let storage = finite(
            input.mass_kg
                * (body.state.specific_enthalpy_j_kg() - input.initial.specific_enthalpy_j_kg()),
        )?;
        let wall = prepared
            .wall
            .evaluate(storage, gas_t, input.internal_heat_j)?;
        let q = -wall.convection_heat_j;
        let (gas_heat, gas_roundoff) = combine(self.gas_heat_j, q)?;
        spend(evaluations, meter.maximum_evaluations)?;
        let phase = prepared.phase.evaluate_at_gas(
            liquid_h(self, body.state.temperature_k())?,
            candidate.gas(),
            gas_heat,
            || cx.checkpoint().is_err(),
        )?;
        prepared.phase.check_balance(&phase)?;
        // The returned wall can change the last phase rate. Reapply the SAME
        // material aggregate allowance to this fresh trial and original gas
        // reconstruction, rather than granting each a new full allowance.
        gate(
            sum(&[
                phase.step().sum_absolute_energy_residual_j,
                candidate.gas_exchange().step().energy_residual_j.abs(),
            ])?,
            self.solve.phase.gas.energy_tolerance_j,
        )?;
        let phase_heat = phase.step().wall_heat_into_system_j;
        let (wall_heat, wall_roundoff) = combine(input.internal_heat_j, -phase_heat)?;
        let wall = prepared.wall.evaluate(storage, gas_t, wall_heat)?;
        gate(wall.residual_j, self.solve.wall.source_tolerance_j)?;
        let wall_trial_energy_residual_j = finite(storage - trial.wall.stored_energy_change_j)?;
        let wall_temperature_residual_k =
            finite(body.state.temperature_k() - trial.wall.temperature_k)?;
        gate(
            wall_trial_energy_residual_j,
            self.solve.wall.source_tolerance_j,
        )?;
        if wall_temperature_residual_k.abs() > self.solve.temperature_tolerance_k {
            return Err(Error::Invalid(
                "original reconstructed wall left its converged phase temperature",
            )
            .into());
        }
        if (body.state.temperature_k() - input.initial.temperature_k()).abs()
            > self.solve.phase.maximum_temperature_increment_k
        {
            return Err(GasEvaporationError::StepSizeLimit {
                maximum_duration_s: (input.duration_s * 0.5).next_down().max(0.0),
                reason: "complete finite wall-film temperature departure",
            }
            .into());
        }
        let sample = body
            .report
            .samples()
            .last()
            .ok_or(Error::Invalid("original wall returned no endpoint"))?;
        let source_residual_j = finite(q + input.duration_s * sample.convection_into_body_w)?;
        gate(source_residual_j, self.solve.sensible.source_tolerance_j)?;
        gate(q - trial.q, self.solve.sensible.source_tolerance_j)?;
        let gas_source_combination_residual_j = -gas_roundoff;
        let wall_source_combination_residual_j = -wall_roundoff;
        let wall_energy_residual_j = sum(&[
            storage,
            -input.internal_heat_j,
            -wall.radiation_heat_j,
            q,
            phase_heat,
        ])?;
        let combined_energy_residual_j = sum(&[
            wall_energy_residual_j,
            phase.step().energy_residual_j,
            gas_source_combination_residual_j,
        ])?;
        let (external, roundoff_a) = fs_math::eft::two_sum(wall.radiation_heat_j, -q);
        let (external_heat_j, roundoff_b) = fs_math::eft::two_sum(external, -phase_heat);
        let external_heat_j = finite(external_heat_j)?;
        let external_heat_roundoff_j = sum(&[roundoff_a, roundoff_b])?;
        let sum_absolute_energy_residual_j = sum(&[
            wall_energy_residual_j
                .abs()
                .max(sample.step_energy_residual_j.abs()),
            phase.step().sum_absolute_energy_residual_j,
            candidate.gas_exchange().step().energy_residual_j.abs(),
            gas_roundoff.abs(),
            wall_roundoff.abs(),
            external_heat_roundoff_j.abs(),
        ])?;
        let mass_residual_kg = sum(&[
            phase.film().mass_kg() - self.film.mass_kg(),
            candidate.gas().mass_kg() - self.gas.mass_kg(),
        ])?;
        let maximum_mass_residual_kg = mass_residual_kg
            .abs()
            .max(phase.maximum_mass_residual_kg())
            .max(
                candidate
                    .gas_exchange()
                    .step()
                    .maximum_species_mass_residual_kg,
            )
            .max(candidate.gas_exchange().step().mass_residual_kg.abs());
        let sum_absolute_mass_residual_kg = sum(&[
            phase.maximum_mass_residual_kg(),
            candidate.gas_exchange().step().mass_residual_kg.abs(),
            candidate
                .gas_exchange()
                .step()
                .maximum_species_mass_residual_kg,
        ])?;
        let worst_mass = maximum_mass_residual_kg.max(sum_absolute_mass_residual_kg);
        if worst_mass > self.solve.phase.gas.mass_tolerance_kg {
            return Err(Error::MassBalance {
                residual_kg: worst_mass,
                tolerance_kg: self.solve.phase.gas.mass_tolerance_kg,
            }
            .into());
        }
        gate(
            candidate.gas_exchange().step().energy_residual_j,
            self.solve.sensible.energy_tolerance_j,
        )?;
        gate(
            sum_absolute_energy_residual_j.max(combined_energy_residual_j.abs()),
            self.solve.wall.energy_tolerance_j,
        )?;
        poll(cx)?;
        Ok(UniformEnthalpyStep {
            state: body.state,
            external_heat_j,
            energy_residual_tolerance_j: self.solve.wall.energy_tolerance_j,
            report: WallFilmImplicitReport {
                next_gas: candidate.gas().clone(),
                next_film: phase.film(),
                phase: phase.step().clone(),
                body: body.report,
                solid_to_gas_heat_j: q,
                wall_heat_into_phase_system_j: phase_heat,
                external_gas_heat_j: self.gas_heat_j,
                radiation_heat_j: wall.radiation_heat_j,
                wall_energy_residual_j,
                gas_energy_residual_j: phase.gas_energy_residual_j(),
                combined_energy_residual_j,
                sum_absolute_energy_residual_j,
                mass_residual_kg,
                maximum_mass_residual_kg,
                sum_absolute_mass_residual_kg,
                wall_trial_energy_residual_j,
                wall_temperature_residual_k,
                source_residual_j,
                gas_source_combination_residual_j,
                wall_source_combination_residual_j,
                external_heat_roundoff_j,
                newton_iterations: run.report.iterations,
                residual_evaluations: evaluations.get() - initial_evaluations,
                gas_iterations: candidate.gas_exchange().step().iterations,
                maximum_accepted_residual_norm_j: run.report.residual_norm,
            },
        })
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint().map_err(|_| Error::Cancelled.into())
}
fn liquid_h(environment: &WallFilmThermalEnvironment<'_>, temperature: f64) -> Result<f64> {
    let h = (environment.boundary.liquid_enthalpy_at_temperature)(temperature)
        .map_err(GasEvaporationError::from)?;
    let t = (environment.boundary.liquid_temperature_from_enthalpy)(h)
        .map_err(GasEvaporationError::from)?;
    if !h.is_finite()
        || !t.is_finite()
        || (t - temperature).abs() > environment.solve.temperature_tolerance_k
    {
        return Err(Error::Invalid("original liquid h(T) and T(h) disagree with the wall").into());
    }
    Ok(h)
}
fn gate(value: f64, tolerance: f64) -> Result<()> {
    if !value.is_finite() || value.abs() > tolerance {
        Err(Error::EnergyBalance {
            residual_j: value,
            tolerance_j: tolerance,
        }
        .into())
    } else {
        Ok(())
    }
}
fn spend(count: &Cell<usize>, limit: usize) -> Result<()> {
    if count.get() >= limit {
        return Err(WallFilmImplicitError::EvaluationBudget);
    }
    count.set(count.get() + 1);
    Ok(())
}
fn difference_error(error: FiniteDifferenceError<WallFilmImplicitError>) -> WallFilmImplicitError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => Error::Invalid(message).into(),
    }
}
fn sum(values: &[f64]) -> Result<f64> {
    let (mut total, mut correction) = (0.0, 0.0);
    for &value in values {
        let (next, remainder) = fs_math::eft::two_sum(total, value);
        total = finite(next)?;
        correction = finite(correction + remainder)?;
    }
    Ok(finite(total + correction)?)
}
