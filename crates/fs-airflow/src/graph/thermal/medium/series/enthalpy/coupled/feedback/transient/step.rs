//! Joint residual adapter; native thermal, momentum and Newton owners do the work.
use super::*;
use fs_conduction::transient::deposition::{
    advection::{Advected, MassFlow},
    backward_euler::{
        BackwardEulerStep, BoundaryStepCoefficients, StepCoefficients, TransportStepCoefficients,
    },
};
use fs_solver::NewtonKrylovState;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use std::cell::Cell;

struct Prepared<T> {
    hydraulics: FiniteLoopHydraulics,
    heat: T,
}

impl FiniteLoopControls {
    fn admit<E>(self) -> Result<(), FiniteLoopError<E>> {
        self.solver.validate().map_err(FiniteLoopError::Solver)?;
        if !(1..=1024).contains(&self.maximum_newton_iterations)
            || !(1..=1_000_000).contains(&self.maximum_residual_evaluations)
            || ![
                self.energy_scale_j,
                self.flow_scale_m3_s,
                self.finite_difference_step,
            ]
            .into_iter()
            .all(positive)
            || ![
                self.cell_tolerance_j,
                self.total_thermal_tolerance_j,
                self.momentum_tolerance_pa_s,
                self.mechanical_energy_tolerance_j,
                self.maximum_neglected_power_w,
                self.maximum_neglected_energy_j,
            ]
            .into_iter()
            .all(nonnegative)
            || !positive(self.maximum_drag_step)
            || self.maximum_drag_step > 2.0
        {
            return Err(FiniteLoopError::Invalid(
                "invalid finite-loop solve scales, tolerances or budgets",
            ));
        }
        Ok(())
    }
}

fn difference_error<E>(error: FiniteDifferenceError<FiniteLoopError<E>>) -> FiniteLoopError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => FiniteLoopError::Invalid(message),
    }
}

impl FiniteLiquidLoop<'_, '_> {
    /// Advance retained pipe momentum and ALL original thermal cells together.
    ///
    /// The pure heat callback supplies the native source/link/boundary
    /// coefficients at the proposed endpoint and actual native average flow.
    /// It may be called repeatedly, including numerical derivative probes.
    /// Its accepted payload is returned verbatim. It must not mutate external
    /// physical histories; callback, source, budget and cancellation failures
    /// leave this owner and the thermal timeline untouched.
    ///
    /// Energy guesses use native stable CellId order (thermal.cells().keys()).
    /// Signed upwind advection uses rho0 times the PHS step-average flow. Thus
    /// an interval ending with reversed velocity can still transport net heat
    /// in the original direction. At exact rest, the same original conductive
    /// and prescribed heat sources continue to change finite pipe enthalpy.
    ///
    /// This is one algebraically coupled discrete step, not an adaptive event
    /// detector for the instant of zero crossing. Refine time and cell spacing.
    /// Pending thermal events obey BackwardEulerStep's original split rules.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'m, T, E>(
        &mut self,
        cx: &Cx<'_>,
        thermal: &mut DepositionSimulation<'m>,
        target_s: f64,
        driving_pressure_pa: f64,
        heat: impl Fn(
            &DepositionSimulation<'m>,
            &FiniteLoopHydraulics,
        ) -> Result<BoundaryStepCoefficients<T>, E>,
        controls: FiniteLoopControls,
        energy_guess_j: &[f64],
    ) -> Result<FiniteLoopStep<T>, FiniteLoopError<E>> {
        poll(cx)?;
        controls.admit()?;
        if thermal.time_s() != self.accepted_time_s || !driving_pressure_pa.is_finite() {
            return Err(FiniteLoopError::Invalid(
                "momentum and original thermal time must agree; finite pressure required",
            ));
        }
        let (next, prepared, thermal_audit, maximum_cell_residual, evaluations, nonlinear) = {
            let step = BackwardEulerStep::new(thermal, target_s, || thermal_poll(cx))?;
            let count = step.cell_ids().len();
            if count > self.limits.maximum_cells
                || energy_guess_j.len() != count
                || energy_guess_j.iter().any(|value| !value.is_finite())
            {
                return Err(FiniteLoopError::Invalid(
                    "complete bounded thermal energy guess required",
                ));
            }
            let evaluations = Cell::new(0usize);
            let evaluate = |x: &[f64]| {
                poll(cx)?;
                if evaluations.get() >= controls.maximum_residual_evaluations {
                    return Err(FiniteLoopError::WorkBudget {
                        performed: evaluations.get(),
                        maximum: controls.maximum_residual_evaluations,
                    });
                }
                evaluations.set(evaluations.get() + 1);
                if x.len() != count + 1 || x.iter().any(|value| !value.is_finite()) {
                    return Err(FiniteLoopError::Invalid(
                        "complete finite loop residual coordinates required",
                    ));
                }
                let energies: Vec<_> = x[..count]
                    .iter()
                    .map(|value| value * controls.energy_scale_j)
                    .collect();
                let candidate_average = x[count] * controls.flow_scale_m3_s;
                if !candidate_average.is_finite() {
                    return Err(FiniteLoopError::Invalid(
                        "candidate average flow is unrepresentable",
                    ));
                }
                step.evaluate_with_transport(
                    &energies,
                    |trial| {
                        let hydraulics = self.hydraulics(
                            cx,
                            trial,
                            candidate_average,
                            driving_pressure_pa,
                            step.duration_s(),
                            controls,
                        )?;
                        let mass = self.reference_density * hydraulics.average_volume_flow_m3_s;
                        let flows = self
                            .cells
                            .iter()
                            .enumerate()
                            .map(|(i, cell)| MassFlow {
                                id: i as u64,
                                a: cell.section.cell,
                                b: self.cells[(i + 1) % self.cells.len()].section.cell,
                                kg_s: mass,
                            })
                            .collect();
                        let advection = Advected::new((), flows, self.cells.len())?
                            .endpoint_transport(trial, || thermal_poll(cx))?;
                        // A closed balanced graph cannot introduce an external stream.
                        if !advection.heat.bulk.sources_w.is_empty() {
                            return Err(FiniteLoopError::Invalid(
                                "closed loop generated external advective heat",
                            ));
                        }
                        let coefficients =
                            heat(trial, &hydraulics).map_err(FiniteLoopError::Heat)?;
                        poll(cx)?;
                        let BoundaryStepCoefficients { bulk, boundaries } = coefficients;
                        let StepCoefficients {
                            sources_w,
                            conductances_w_k,
                            payload,
                        } = bulk;
                        Ok(TransportStepCoefficients {
                            heat: BoundaryStepCoefficients {
                                bulk: StepCoefficients {
                                    sources_w,
                                    conductances_w_k,
                                    payload: Prepared {
                                        hydraulics,
                                        heat: payload,
                                    },
                                },
                                boundaries,
                            },
                            transfers: advection.transfers,
                        })
                    },
                    || thermal_poll(cx),
                )
            };
            let problem = FallibleFiniteDifference::new(
                count + 1,
                controls.finite_difference_step,
                |x: &[f64]| {
                    let trial = evaluate(x)?;
                    let mut residual: Vec<_> = trial
                        .residual_j()
                        .iter()
                        .map(|value| value / controls.energy_scale_j)
                        .collect();
                    residual.push(
                        x[count]
                            - trial.payload().hydraulics.average_volume_flow_m3_s
                                / controls.flow_scale_m3_s,
                    );
                    Ok(residual)
                },
            )
            .map_err(difference_error)?;
            let mut guess: Vec<_> = energy_guess_j
                .iter()
                .map(|value| value / controls.energy_scale_j)
                .collect();
            guess.push(self.volume_flow_m3_s() / controls.flow_scale_m3_s);
            let initialized = NewtonKrylovState::new(&problem, guess, controls.solver);
            if let Some(error) = problem.take_failure() {
                return Err(difference_error(error));
            }
            let mut solver = initialized.map_err(FiniteLoopError::Solver)?;
            let run =
                solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                    problem.has_failed() || cx.checkpoint().is_err()
                });
            if let Some(error) = problem.take_failure() {
                return Err(difference_error(error));
            }
            poll(cx)?;
            if !run.report.converged {
                return Err(FiniteLoopError::NotConverged(run.report));
            }
            let final_trial = evaluate(&solver.x)?;
            let hydraulic = &final_trial.payload().hydraulics;
            if hydraulic.momentum_residual_pa_s.abs() > controls.momentum_tolerance_pa_s
                || hydraulic.mechanical_balance_residual_j.abs()
                    > controls.mechanical_energy_tolerance_j
            {
                return Err(FiniteLoopError::MomentumResidual {
                    impulse_pa_s: hydraulic.momentum_residual_pa_s,
                    energy_j: hydraulic.mechanical_balance_residual_j,
                });
            }
            let maximum_cell_residual = final_trial.maximum_cell_residual_j();
            let (next, prepared, audit) = final_trial.accept(
                controls.cell_tolerance_j,
                controls.total_thermal_tolerance_j,
                || thermal_poll(cx),
            )?;
            (
                next,
                prepared,
                audit,
                maximum_cell_residual,
                evaluations.get(),
                run.report,
            )
        };
        // Last cancellation boundary precedes both publications; nothing after
        // this point can fail or advance an independently owned physical clock.
        poll(cx)?;
        self.momentum_pa_s = prepared.hydraulics.momentum.x[0];
        self.accepted_time_s = target_s;
        *thermal = next;
        Ok(FiniteLoopStep {
            hydraulics: prepared.hydraulics,
            heat: prepared.heat,
            thermal: thermal_audit,
            maximum_cell_residual_j: maximum_cell_residual,
            residual_evaluations: evaluations,
            nonlinear,
        })
    }
}
