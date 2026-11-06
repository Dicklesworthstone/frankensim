//! Original bounded Newton--Krylov and one final original gas reconstruction.
use super::*;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonKrylovState, SolverRunProgress};

/// Solve one event-free native thermal-network/film/gas interval.
///
/// The selected cell keeps its original mass, solid/liquid phase chart and
/// graph contacts. The patch occupies only currently exposed ambient area;
/// callback replacements of that cell's ambient channel refuse as ambiguous.
/// Other native source, contact, reservoir and boundary channels retain their
/// original endpoint resolution. The callback sees the actual trial gas and
/// thermal coordinates; only its final accepted payload is returned.
///
/// No wall or gas roots occur in residual evaluations. Original gas calorics
/// reconstruct once after convergence, followed by a fresh native heat/phase
/// check against that actual gas. Events at target remain pending in the
/// native left-limit state. This is first-order in physical time; source and
/// model-domain refusals propagate without extrapolation or artificial films.
#[allow(clippy::too_many_arguments)]
pub fn advance_deposition_with_wall_film<'m, 'catalog, T>(
    thermal: &DepositionSimulation<'m>,
    binding: DepositionWallFilmBinding<'_>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    external_gas_heat_j: f64,
    controls: DepositionWallFilmControls,
    resolve_heat: impl FnMut(
        &DepositionSimulation<'m>,
        &ClosedGasParcel,
    ) -> Result<BoundaryExchangeStepCoefficients<T>, ImplicitGasError>,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<DepositionWallFilmCandidate<'m, 'catalog, T>, ImplicitGasError> {
    checkpoint()?;
    let meter = WorkMeter::new(
        controls.numerics(),
        controls.implicit.phase.gas.gas.maximum_iterations,
    )
    .map_err(work_error)?;
    advance_deposition_with_wall_film_metered(
        thermal,
        binding,
        gas,
        catalog,
        target_s,
        external_gas_heat_j,
        controls,
        &meter,
        resolve_heat,
        checkpoint,
    )
}

#[allow(clippy::too_many_arguments)]
pub(in super::super) fn advance_deposition_with_wall_film_metered<'m, 'catalog, T>(
    thermal: &DepositionSimulation<'m>,
    binding: DepositionWallFilmBinding<'_>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    external_gas_heat_j: f64,
    mut controls: DepositionWallFilmControls,
    meter: &WorkMeter,
    mut resolve_heat: impl FnMut(
        &DepositionSimulation<'m>,
        &ClosedGasParcel,
    ) -> Result<BoundaryExchangeStepCoefficients<T>, ImplicitGasError>,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<DepositionWallFilmCandidate<'m, 'catalog, T>, ImplicitGasError> {
    checkpoint()?;
    let iterations = meter
        .remaining_iterations()
        .map_err(work_error)?
        .min(controls.implicit.thermal.maximum_iterations);
    let gas_limit = meter
        .gas_limit(controls.implicit.phase.gas.gas.maximum_iterations)
        .map_err(work_error)?;
    controls.implicit.phase.gas.gas.maximum_iterations = gas_limit;
    let prepared = Prepared::new(
        thermal,
        binding,
        gas,
        catalog,
        target_s,
        external_gas_heat_j,
        controls,
        &mut checkpoint,
    )?;
    let poll = RefCell::new(&mut checkpoint);
    let source = RefCell::new(&mut resolve_heat);
    let initial_evaluations = meter.residual_evaluations();
    meter.spend_evaluation().map_err(work_error)?;
    let zero = prepared.project(&vec![0.0; prepared.dimension()]);
    let predictor = prepared.evaluate(&zero, None, &mut *source.borrow_mut(), &mut || {
        (poll.borrow_mut())()
    })?;
    let mut guess = zero.clone();
    guess[prepared.dimension() - 2] = finite(
        predictor.heat.payload().phase.step().transferred_mass_kg
            * controls.implicit.mass_coordinate_scale_j_kg,
    )?;
    let cache = RefCell::new(vec![(zero, predictor.residual_j)]);
    drop(predictor.heat);
    let problem = FallibleFiniteDifference::new(
        prepared.dimension(),
        controls.implicit.thermal.difference_step_j,
        |point: &[f64]| {
            (poll.borrow_mut())()?;
            let point = prepared.project(point);
            {
                let mut entries = cache.borrow_mut();
                if let Some(index) = entries.iter().position(|(old, _)| {
                    old.iter()
                        .zip(&point)
                        .all(|(a, b)| a.to_bits() == b.to_bits())
                }) {
                    let entry = entries.remove(index);
                    let residual = entry.1.clone();
                    entries.push(entry);
                    return Ok(residual);
                }
            }
            meter.spend_evaluation().map_err(work_error)?;
            let trial = prepared.evaluate(&point, None, &mut *source.borrow_mut(), &mut || {
                (poll.borrow_mut())()
            })?;
            let mut entries = cache.borrow_mut();
            if entries.len() == 2 {
                entries.remove(0);
            }
            entries.push((point, trial.residual_j.clone()));
            Ok::<_, ImplicitGasError>(trial.residual_j)
        },
    )
    .map_err(difference_error)?;
    let initialization = NewtonKrylovState::new(
        &problem,
        prepared.project(&guess),
        controls.implicit.thermal.newton,
    );
    if let Some(error) = problem.take_failure() {
        return Err(difference_error(error));
    }
    let mut solver = initialization?;
    let mut interrupted = None;
    let run = solver.run_controlled(&problem, iterations, &mut || {
        if problem.has_failed() {
            return true;
        }
        match (poll.borrow_mut())() {
            Ok(()) => false,
            Err(error) => {
                interrupted = Some(error);
                true
            }
        }
    });
    meter
        .record_iterations(run.report.iterations)
        .map_err(work_error)?;
    if let Some(error) = problem.take_failure() {
        return Err(difference_error(error));
    }
    if let Some(error) = interrupted {
        return Err(error.into());
    }
    if run.progress != SolverRunProgress::Complete {
        return Err(ImplicitGasError::Interrupted(run.progress));
    }
    if !run.report.converged {
        return Err(ImplicitGasError::NonConverged(run.report));
    }
    let point = prepared.project(&solver.x);
    meter.spend_evaluation().map_err(work_error)?;
    let trial = prepared.evaluate(&point, None, &mut *source.borrow_mut(), &mut || {
        (poll.borrow_mut())()
    })?;
    prepared.check(&trial, &mut || (poll.borrow_mut())())?;
    meter.spend_evaluation().map_err(work_error)?;
    let reserved = meter.reserve_gas(gas_limit).map_err(work_error)?;
    let candidate = phase_call(&mut || (poll.borrow_mut())(), |cancel| {
        prepared.phase.finalize(&trial.heat.payload().phase, cancel)
    })?;
    meter
        .refund_gas(reserved, candidate.gas_exchange().step().iterations)
        .map_err(work_error)?;
    drop(trial);
    meter.spend_evaluation().map_err(work_error)?;
    let fresh = prepared.evaluate(
        &point,
        Some(candidate.gas()),
        &mut *source.borrow_mut(),
        &mut || (poll.borrow_mut())(),
    )?;
    prepared.accept(
        fresh,
        candidate.gas_exchange().clone(),
        run.report.iterations,
        meter.residual_evaluations() - initial_evaluations,
        run.report.residual_norm,
        &mut || (poll.borrow_mut())(),
    )
}

fn difference_error(error: FiniteDifferenceError<ImplicitGasError>) -> ImplicitGasError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => ImplicitGasError::Invalid(message),
    }
}
