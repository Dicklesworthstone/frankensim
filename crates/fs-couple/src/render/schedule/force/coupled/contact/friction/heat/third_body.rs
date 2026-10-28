//! Original modal motion and a finite third-body enthalpy cell accept together.
//! The layer's existing caloric law owns temperature; there is no massless
//! temperature coordinate or added heat capacity. The same material endpoint
//! supplies friction, all three conductivities and both contact resistances.

use super::*;
use crate::friction_heat::third_body::{
    ThirdBodyHeatInput, ThirdBodyHeatReport, ThirdBodyHeatRoute, ThirdBodyHeatStep,
    ThirdBodyHeatTolerances,
};

pub struct ModalThirdBodyHeatSource<P> {
    pub law: FrictionLaw,
    /// Original endpoint conductivities [body A, body B, layer], W/(m K).
    pub conductivity_w_m_k: [f64; 3],
    /// Original source contact resistances [A/layer, B/layer], m² K/W.
    pub area_resistance_m2_k_w: [f64; 2],
    /// Other heat only. Physical friction work comes from the motion owner.
    pub transport: StepCoefficients<P>,
}

#[derive(Clone, Copy, Debug)]
pub struct ModalThirdBodyHeatSolve {
    pub maximum_cells: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    pub heat: ThirdBodyHeatTolerances,
    pub combined_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

pub struct ThirdBodyModalFrictionFrame<P> {
    pub mechanical: FrictionModalFrame,
    pub heat: ThirdBodyHeatReport,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub combined_energy_residual_j: f64,
}

impl FrictionModalSystem {
    /// Resolve motion, finite layer storage and conductive body heating in one
    /// transaction. The callback sees the SAME trial layer temperature used by
    /// the original thermal network. Its candidate source/material history must
    /// remain unpublished until this function returns success.
    ///
    /// The route's body cells must represent this port's contacting bodies.
    /// Retained dry-interface cards preserve their source and applicability
    /// through the parent owner's temperature continuation. Caller laws keep
    /// caller authority. This method does not turn a dry law into an elastic
    /// third-body law; retained shear layers use their existing constitutive
    /// owner and the general `ThirdBodyHeatStep` or original layer thermal path.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_with_third_body_heat<'m, P, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'m>,
        route: &ThirdBodyHeatRoute,
        external: &[f64],
        drive: FrictionDrive,
        controls: ModalThirdBodyHeatSolve,
        energy_guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<ModalThirdBodyHeatSource<P>, E>,
        gate: &CancelGate,
    ) -> Result<ThirdBodyModalFrictionFrame<P>, ModalFrictionHeatError<E>> {
        if self.adhesive.is_some() {
            return Err(ModalFrictionHeatError::Invalid(
                "dry-law thermal refresh cannot replace a frozen adhesive normal/shear source",
            ));
        }
        let checkpoint = || {
            if gate.is_requested() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        controls.newton.validate()?;
        let sample = self.samples_rendered();
        let next_sample = sample
            .checked_add(1)
            .ok_or(ModalFrictionHeatError::Invalid(
                "modal sample clock overflow",
            ))?;
        let start = sample as f64 * self.sample_period_s();
        let target = next_sample as f64 * self.sample_period_s();
        if thermal.time_s() != start
            || !target.is_finite()
            || target <= start
            || controls.maximum_cells == 0
            || thermal.cells().len() > controls.maximum_cells
            || controls.maximum_trials == 0
            || controls.maximum_newton_iterations == 0
            || !controls.difference_step_j.is_finite()
            || controls.difference_step_j <= 0.0
            || !controls.combined_energy_tolerance_j.is_finite()
            || controls.combined_energy_tolerance_j < 0.0
            || energy_guess_j.len() != thermal.cells().len()
            || energy_guess_j.iter().any(|x| !x.is_finite())
        {
            return Err(ModalFrictionHeatError::Invalid(
                "synchronized clocks and finite bounded third-body coordinates required",
            ));
        }
        let step = ThirdBodyHeatStep::new(thermal, target, route, checkpoint)?;
        let cells = step.cell_ids().len();
        let trials = Cell::new(0usize);
        let owner = RefCell::new(&mut *self);
        let evaluate = |x: &[f64]| {
            checkpoint()?;
            if trials.get() >= controls.maximum_trials {
                return Err(ModalFrictionHeatError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            if x.len() != cells {
                return Err(ModalFrictionHeatError::Invalid(
                    "incomplete third-body heat iterate",
                ));
            }
            step.evaluate(
                x,
                |endpoint, layer_temperature| {
                    let ModalThirdBodyHeatSource {
                        law,
                        conductivity_w_m_k,
                        area_resistance_m2_k_w,
                        transport,
                    } = resolve(endpoint, layer_temperature)
                        .map_err(ModalFrictionHeatError::Source)?;
                    let StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload,
                    } = transport;
                    let staged = owner.borrow_mut().stage_at_temperature(
                        law,
                        layer_temperature,
                        external,
                        drive,
                        payload,
                        gate,
                    )?;
                    Ok(ThirdBodyHeatInput {
                        dissipated_j: staged.mechanical.friction_dissipation_j,
                        conductivity_w_m_k,
                        area_resistance_m2_k_w,
                        transport: StepCoefficients {
                            sources_w,
                            conductances_w_k,
                            payload: staged,
                        },
                    })
                },
                checkpoint,
            )
        };
        let problem =
            FallibleFiniteDifference::new(cells, controls.difference_step_j, |x: &[f64]| {
                evaluate(x).map(|trial| trial.residual_j().to_vec())
            })
            .map_err(numerical)?;
        let construction =
            NewtonKrylovState::new(&problem, energy_guess_j.to_vec(), controls.newton);
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        let mut solver = construction?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                gate.is_requested() || problem.has_failed()
            });
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        checkpoint()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(ModalFrictionHeatError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(ModalFrictionHeatError::NotConverged(result.report));
        }
        let trial = evaluate(&solver.x)?;
        let combined =
            trial.payload().mechanical.energy_residual_j + trial.report().total_energy_residual_j;
        if !combined.is_finite() || combined.abs() > controls.combined_energy_tolerance_j {
            return Err(ModalFrictionHeatError::Invalid(
                "joint modal/third-body work does not balance",
            ));
        }
        let (next_thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
        checkpoint()?;
        drop(problem);
        let owner = owner.into_inner();
        // No fallible work or external callback follows this common boundary.
        owner.port.law = staged.law;
        owner.interface = staged.binding;
        owner.load_line = staged.load_line;
        owner.staged_branch = staged.branch;
        owner.publish_candidate(staged.mechanical, staged.nonlinear);
        *thermal = next_thermal;
        Ok(ThirdBodyModalFrictionFrame {
            mechanical: staged.mechanical,
            heat,
            source: staged.source,
            nonlinear: result.report,
            trials: trials.get(),
            combined_energy_residual_j: combined,
        })
    }
}
