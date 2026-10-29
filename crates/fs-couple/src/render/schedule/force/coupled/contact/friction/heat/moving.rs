//! Distributed Coulomb friction and moving finite-network heat in one sample.
//!
//! Each swept patch resolves its own normal load and temperature-dependent
//! law. All patches share the original scalar port's actual relative slip;
//! their static and kinetic capacities add. Local heat follows local traction
//! times that same slip. Neither uniform pressure nor uniform heating is
//! inferred. The prescribed geometric sweep is independent input: applications
//! must derive it from their original motion/geometry owner and refine in time.

use super::*;
use crate::friction_heat::moving::{
    MovingContactPath, MovingFrictionHeatReport, MovingFrictionHeatStep, MovingFrictionInput,
    MovingFrictionRoute, SweptContactPatch,
};
use fs_thermal::enthalpy::CellId;
use std::collections::BTreeMap;

/// Original local normal load and friction law at one swept patch temperature.
pub struct MovingModalPatchSource {
    /// Mean compressive force [N] over the interval, from the normal owner.
    pub normal_force_n: f64,
    /// Coulomb static/kinetic coefficients resolved at this patch's temperature.
    pub law: FrictionLaw,
}

/// Endpoint sources in the exact canonical order returned by the moving step.
pub struct MovingModalHeatSource<P> {
    pub patches: Vec<MovingModalPatchSource>,
    pub conductivity_w_m_k: BTreeMap<CellId, f64>,
    /// Other thermal transport only; mechanics supplies all friction heat.
    pub transport: StepCoefficients<P>,
}

pub struct HeatedMovingModalFrictionFrame<P> {
    pub mechanical: FrictionModalFrame,
    pub heat: MovingFrictionHeatReport,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub combined_energy_residual_j: f64,
}

impl FrictionModalSystem {
    /// Advance the original mechanical port and distributed moving heat together.
    ///
    /// Every patch uses its own solved interface temperature and source law.
    /// Local normal forces must sum to the original drive's normal force. This
    /// is the common-slip Coulomb rung; non-Coulomb laws and a retained single
    /// interface binding refuse rather than discard their source/history rules.
    /// Finite-patch partial slip belongs to the tangential-island owner.
    ///
    /// The supplied path must describe the actual material/footprint sweep over
    /// this sample. It is prescribed geometry, not inferred from a guessed slip.
    /// A full moving thermo-mechanical trajectory needs temporal refinement.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_with_moving_friction_heat<'m, P, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'m>,
        route: &MovingFrictionRoute,
        path: MovingContactPath,
        external: &[f64],
        drive: FrictionDrive,
        controls: ModalFrictionHeatSolve,
        energy_guess_j: &[f64],
        temperature_guess_k: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'m>,
            &[SweptContactPatch],
            &[f64],
        ) -> Result<MovingModalHeatSource<P>, E>,
        gate: &CancelGate,
    ) -> Result<HeatedMovingModalFrictionFrame<P>, ModalFrictionHeatError<E>> {
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
        let next = sample
            .checked_add(1)
            .ok_or(ModalFrictionHeatError::Invalid(
                "modal sample clock overflow",
            ))?;
        let start = sample as f64 * self.sample_period_s();
        let target = next as f64 * self.sample_period_s();
        if self.interface.is_some()
            || thermal.time_s() != start
            || !target.is_finite()
            || target <= start
            || controls.maximum_cells == 0
            || thermal.cells().len() > controls.maximum_cells
            || controls.maximum_trials == 0
            || controls.maximum_newton_iterations == 0
            || !controls.difference_step_j.is_finite()
            || controls.difference_step_j <= 0.0
            || !controls.temperature_scale_j_per_k.is_finite()
            || controls.temperature_scale_j_per_k <= 0.0
            || !controls.combined_energy_tolerance_j.is_finite()
            || controls.combined_energy_tolerance_j < 0.0
            || !drive.normal_force_n.is_finite()
            || drive.normal_force_n < 0.0
            || energy_guess_j.len() != thermal.cells().len()
            || energy_guess_j.iter().any(|x| !x.is_finite())
        {
            return Err(ModalFrictionHeatError::Invalid(
                "moving heat needs an unbound common-slip port, synchronized clocks and finite bounded controls",
            ));
        }
        let step = MovingFrictionHeatStep::new(thermal, target, route, path, checkpoint)?;
        let cells = step.cell_ids().len();
        let patch_count = step.patches().len();
        if patch_count == 0
            || temperature_guess_k.len() != patch_count
            || temperature_guess_k
                .iter()
                .any(|t| !t.is_finite() || *t <= 0.0)
        {
            return Err(ModalFrictionHeatError::Invalid(
                "every swept patch needs a positive interface-temperature guess",
            ));
        }
        let dimension = cells
            .checked_add(patch_count)
            .ok_or(ModalFrictionHeatError::Invalid(
                "moving heat coordinate count overflow",
            ))?;
        let normal_tolerance = self.config.force_absolute_tolerance_n
            + self.config.force_relative_tolerance * drive.normal_force_n;
        let trials = Cell::new(0usize);
        let owner = RefCell::new(&mut *self);
        let evaluate = |x: &[f64]| {
            checkpoint()?;
            if trials.get() >= controls.maximum_trials {
                return Err(ModalFrictionHeatError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            if x.len() != dimension {
                return Err(ModalFrictionHeatError::Invalid(
                    "incomplete moving thermal iterate",
                ));
            }
            let temperatures: Vec<f64> = temperature_guess_k
                .iter()
                .zip(&x[cells..])
                .map(|(t, increment)| t + increment / controls.temperature_scale_j_per_k)
                .collect();
            step.evaluate(&x[..cells], &temperatures, |endpoint, patches, local_t| {
                let MovingModalHeatSource { patches: local, conductivity_w_m_k, transport } =
                    resolve(endpoint, patches, local_t).map_err(ModalFrictionHeatError::Source)?;
                if local.len() != patch_count { return Err(ModalFrictionHeatError::Invalid("moving mechanical source omitted a patch")); }
                let mut normal_sum = 0.0;
                let mut static_mu = 0.0;
                let mut kinetic_mu = 0.0;
                let mut kinetic_weight = Vec::with_capacity(patch_count);
                for patch in &local {
                    if !patch.normal_force_n.is_finite() || patch.normal_force_n < 0.0 {
                        return Err(ModalFrictionHeatError::Invalid("moving patch normal force must be finite and nonnegative"));
                    }
                    patch.law.static_coefficient().map_err(ModalFrictionError::from)?;
                    let FrictionLaw::Coulomb { static_mu: ms, kinetic_mu: mk } = patch.law else {
                        return Err(ModalFrictionHeatError::Invalid("moving common-slip composition currently requires local Coulomb laws"));
                    };
                    normal_sum += patch.normal_force_n;
                    let weight = if drive.normal_force_n > 0.0 { patch.normal_force_n / drive.normal_force_n } else { 0.0 };
                    static_mu += weight * ms;
                    kinetic_mu += weight * mk;
                    kinetic_weight.push(weight * mk);
                }
                if !normal_sum.is_finite() || (normal_sum - drive.normal_force_n).abs() > normal_tolerance {
                    return Err(ModalFrictionHeatError::Invalid("local patch normal forces do not equal the original normal load"));
                }
                let law = FrictionLaw::Coulomb { static_mu, kinetic_mu };
                let StepCoefficients { sources_w, conductances_w_k, payload } = transport;
                // This port has no single retained interface binding. Its local
                // laws above already consumed ALL temperatures; the temperature
                // argument below is used only by that excluded binding path.
                let staged = owner.borrow_mut().stage_at_temperature(law, local_t[0], external, drive, payload, gate)?;
                let work = staged.mechanical.friction_dissipation_j;
                if work > 0.0 && kinetic_mu <= 0.0 {
                    return Err(ModalFrictionHeatError::Invalid("positive friction work has no local kinetic capacity"));
                }
                let dissipated_j = kinetic_weight.into_iter().map(|weight|
                    if work > 0.0 { work * (weight / kinetic_mu) } else { 0.0 }).collect();
                Ok(MovingFrictionInput { dissipated_j, conductivity_w_m_k,
                    transport: StepCoefficients { sources_w, conductances_w_k, payload: staged } })
            }, checkpoint)
        };
        let problem =
            FallibleFiniteDifference::new(dimension, controls.difference_step_j, |x: &[f64]| {
                evaluate(x).map(|trial| trial.residual_j().to_vec())
            })
            .map_err(numerical)?;
        let mut initial = energy_guess_j.to_vec();
        initial.resize(dimension, 0.0);
        let construction = NewtonKrylovState::new(&problem, initial, controls.newton);
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
                "moving mechanical work and original thermal storage do not balance",
            ));
        }
        let (next_thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
        checkpoint()?;
        drop(problem);
        let owner = owner.into_inner();
        owner.network.validate_staged_material_points(Some(gate))
            .map_err(ModalFrictionError::from)?;
        owner.port.law = staged.law;
        owner.interface = staged.binding;
        owner.load_line = staged.load_line;
        owner.staged_branch = staged.branch;
        owner.publish_candidate(staged.mechanical, staged.nonlinear);
        *thermal = next_thermal;
        Ok(HeatedMovingModalFrictionFrame {
            mechanical: staged.mechanical,
            heat,
            source: staged.source,
            nonlinear: result.report,
            trials: trials.get(),
            combined_energy_residual_j: combined,
        })
    }
}
