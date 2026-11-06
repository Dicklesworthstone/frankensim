//! One original Newton--Krylov solve, followed by one gas reconstruction and
//! a fresh coupled physical check. Cached values are residuals, never accepted
//! source payloads or gas/thermal endpoints.

use super::endpoint::Prepared;
use super::*;
use fs_solver::NewtonKrylovState;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use std::cell::{Cell, RefCell};

/// Solve one event-free thermal/finite-gas interval on the original owners.
///
/// Every actual residual evaluation, including the physical initial predictor,
/// fresh converged endpoint and reconstructed-gas recheck, spends the same
/// explicit evaluation budget. Only the original gas boundary owner reconstructs
/// its accepted caloric state, once after nonlinear convergence. Both endpoints
/// are published together only after rechecking their phase and heat streams.
/// Source and cancellation failures propagate without a substitute atmosphere,
/// caloric law or numerical retry that resets the work allowance.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn solve_deposition_with_gas<'m, 'catalog, T>(
    thermal: &DepositionSimulation<'m>,
    surfaces: &[ImplicitGasSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    controls: ImplicitGasControls,
    mut resolve_heat: impl FnMut(
        &DepositionSimulation<'m>,
        &BTreeMap<CellId, LiquidSurfaceState>,
        &ClosedGasParcel,
    ) -> Result<BoundaryExchangeStepCoefficients<T>, ImplicitGasError>,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<ImplicitGasCandidate<'m, 'catalog, T>, ImplicitGasError> {
    checkpoint()?;
    let prepared = Prepared::new(
        thermal,
        surfaces,
        gas,
        catalog,
        target_s,
        controls,
        &mut checkpoint,
    )?;
    let poll = RefCell::new(&mut checkpoint);
    let source = RefCell::new(&mut resolve_heat);
    let evaluations = Cell::new(0usize);

    // The physical zero-coordinate endpoint supplies the represented initial
    // phase transfer. This is a bounded starting guess, not manufactured trace
    // vapor or an independently advanced heat/gas state.
    spend_evaluation(&evaluations, controls.thermal.maximum_evaluations)?;
    let zero = prepared.project(&vec![0.0; prepared.dimension()]);
    let predictor = prepared.evaluate(&zero, &mut *source.borrow_mut(), &mut || {
        (poll.borrow_mut())()
    })?;
    let guess = prepared.project(&prepared.initial_guess(&predictor));
    let zero_residual = predictor.residual_j();
    drop(predictor);
    let mut initial_cache = Vec::with_capacity(2);
    initial_cache.push((zero, zero_residual));
    let cache = RefCell::new(initial_cache);

    let problem = FallibleFiniteDifference::new(
        prepared.dimension(),
        controls.thermal.difference_step_j,
        |coordinates: &[f64]| {
            (poll.borrow_mut())()?;
            let coordinates = prepared.project(coordinates);
            {
                let mut entries = cache.borrow_mut();
                if let Some(index) = entries
                    .iter()
                    .position(|(point, _)| same_point(point, &coordinates))
                {
                    let entry = entries.remove(index);
                    let residual = entry.1.clone();
                    entries.push(entry);
                    return Ok(residual);
                }
            }
            spend_evaluation(&evaluations, controls.thermal.maximum_evaluations)?;
            let trial = prepared.evaluate(&coordinates, &mut *source.borrow_mut(), &mut || {
                (poll.borrow_mut())()
            })?;
            let residual = trial.residual_j();
            let mut entries = cache.borrow_mut();
            if entries.len() == 2 {
                entries.remove(0);
            }
            entries.push((coordinates, residual.clone()));
            Ok::<_, ImplicitGasError>(residual)
        },
    )
    .map_err(difference_error)?;
    let initialization = NewtonKrylovState::new(&problem, guess, controls.thermal.newton);
    if let Some(error) = problem.take_failure() {
        return Err(difference_error(error));
    }
    let mut solver = initialization?;
    let mut checkpoint_error = None;
    let run = solver.run_controlled(&problem, controls.thermal.maximum_iterations, &mut || {
        if problem.has_failed() {
            return true;
        }
        match (poll.borrow_mut())() {
            Ok(()) => false,
            Err(error) => {
                checkpoint_error = Some(error);
                true
            }
        }
    });
    if let Some(error) = problem.take_failure() {
        return Err(difference_error(error));
    }
    if let Some(error) = checkpoint_error {
        return Err(error.into());
    }
    if run.progress != SolverRunProgress::Complete {
        return Err(ImplicitGasError::Interrupted(run.progress));
    }
    if !run.report.converged {
        return Err(ImplicitGasError::NonConverged(run.report));
    }

    let coordinates = prepared.project(&solver.x);
    spend_evaluation(&evaluations, controls.thermal.maximum_evaluations)?;
    let trial = prepared.evaluate(&coordinates, &mut *source.borrow_mut(), &mut || {
        (poll.borrow_mut())()
    })?;
    prepared.check_balance(&trial)?;
    let gas = prepared.reconstruct(&trial, &mut || (poll.borrow_mut())())?;
    drop(trial);

    // Reconstructing the actual gas can change endpoint films and streams by
    // its admitted numerical tolerance. Re-evaluate against that real state;
    // never install a thermal endpoint whose opposite gas account is stale.
    spend_evaluation(&evaluations, controls.thermal.maximum_evaluations)?;
    let trial = prepared.evaluate_at_gas(
        &coordinates[..prepared.dimension() - 2],
        gas.gas(),
        &mut *source.borrow_mut(),
        &mut || (poll.borrow_mut())(),
    )?;
    prepared.accept(
        trial,
        gas,
        run.report.iterations,
        evaluations.get(),
        run.report.residual_norm,
        &mut || (poll.borrow_mut())(),
    )
}

fn spend_evaluation(count: &Cell<usize>, limit: usize) -> Result<(), ImplicitGasError> {
    if count.get() >= limit {
        return Err(DepositionError::Budget("implicit finite-gas residual evaluations").into());
    }
    count.set(count.get() + 1);
    Ok(())
}

fn same_point(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(a, b)| a.to_bits() == b.to_bits())
}

fn difference_error(error: FiniteDifferenceError<ImplicitGasError>) -> ImplicitGasError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => ImplicitGasError::Invalid(message),
    }
}
