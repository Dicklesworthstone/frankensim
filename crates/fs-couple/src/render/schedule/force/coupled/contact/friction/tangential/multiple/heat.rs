//! Simultaneous modal shear and shared finite-cell heat in one transaction.
//!
//! Each original contact has a distinct massless interface temperature. All
//! mechanical trials retain the complete contact-to-contact modal response;
//! all thermal trials retain the same original finite cells and conductances.
//! Shared body motion and shared caloric storage are each counted once.

use super::{
    MultiTangentialModalCandidate, MultiTangentialModalError, MultiTangentialModalFrame,
    MultiTangentialModalSystem, PlanarTangentialModalDrive, TangentialLaw,
};
use crate::friction_heat::FrictionHeatRoute;
use crate::friction_heat::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use fs_exec::CancelGate;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovState, NewtonReport, SolverRunProgress};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::{Cell, RefCell};

pub use super::super::heat::ModalFrictionHeatSolve;

/// Explicit one-to-one binding of mechanical history to its original heat route.
#[derive(Clone, Copy)]
pub struct ModalTangentialHeatRoute<'a> {
    /// Index in the mechanical owner's original contact construction order.
    pub contact_index: usize,
    /// This contact's actual two-body thermal link, in its original A/B order.
    pub route: &'a FrictionHeatRoute,
    /// Numerical initial guess only; the common temperature is solved here.
    pub temperature_guess_k: f64,
}

/// Bounded joint solve using the existing cell, work and nonlinear controls.
#[derive(Clone, Copy, Debug)]
pub struct MultiTangentialHeatSolve {
    pub solve: ModalFrictionHeatSolve,
    pub maximum_interfaces: usize,
    /// Bound original cells plus common interface temperatures before allocation.
    pub maximum_unknowns: usize,
    /// Audit original caloric entropy at the final root without another solve.
    /// All active material laws must supply their exact entropy integral.
    pub audit_entropy: bool,
}

/// Original conductivities for one named physical interface.
#[derive(Clone, Copy, Debug)]
pub struct ModalContactHeatSource {
    pub link_id: u64,
    /// Ordered like the bound route's original A/B cells.
    pub conductivity_w_m_k: [f64; 2],
}

/// One complete endpoint material/transport evaluation for the entire network.
pub struct MultiTangentialHeatSource<P> {
    /// Exactly the callback's canonical ascending original link ID order.
    pub interfaces: Vec<ModalContactHeatSource>,
    /// Other heat only; neither contact heat nor bound links may be reinserted.
    pub transport: StepCoefficients<P>,
}

/// Original source, simultaneous mechanical, thermal or nonlinear refusal.
#[derive(Debug)]
pub enum MultiTangentialHeatError<E> {
    Invalid(&'static str),
    Mechanical(MultiTangentialModalError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}

impl<E> From<MultiTangentialModalError> for MultiTangentialHeatError<E> {
    fn from(error: MultiTangentialModalError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<DepositionError> for MultiTangentialHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<NewtonError> for MultiTangentialHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for MultiTangentialHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "simultaneous modal friction/heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for MultiTangentialHeatError<E> {}

/// One accepted mechanical/contact/thermal interval for the whole contact set.
pub struct HeatedMultiTangentialModalFrame<P> {
    /// Contacts retain their original mechanical construction order.
    pub mechanical: MultiTangentialModalFrame,
    /// Interfaces follow canonical ascending original thermal link ID order.
    pub heat: FrictionHeatNetworkReport,
    /// `contact_indices[i]` names the mechanical contact for `heat.interfaces[i]`.
    pub contact_indices: Vec<usize>,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    /// Complete mechanical plus caloric work, with friction an internal transfer.
    pub combined_energy_residual_j: f64,
}

struct Staged<S, P> {
    mechanical: MultiTangentialModalCandidate<S>,
    source: P,
}

fn numerical<E>(
    error: FiniteDifferenceError<MultiTangentialHeatError<E>>,
) -> MultiTangentialHeatError<E> {
    match error {
        FiniteDifferenceError::Source(source) => source,
        FiniteDifferenceError::Invalid(message) => MultiTangentialHeatError::Invalid(message),
    }
}

impl<L: TangentialLaw> MultiTangentialModalSystem<L> {
    /// Solve all retained contacts and their shared original finite cells together.
    ///
    /// Bind each mechanical contact exactly once to a distinct original thermal
    /// link. Binding enumeration is arbitrary; callbacks and heat reports use
    /// ascending link IDs. Drives and mechanical history retain construction
    /// order. The explicit index map relates those orders without reordering or
    /// reconstructing constitutive states. The caller owns the physical mapping
    /// between modal bodies and route cells; material identity is not inferred.
    ///
    /// Every drive temperature must be absent because these temperatures are
    /// simultaneous unknowns. The callback receives the whole proposed thermal
    /// endpoint and all common temperatures once per joint trial. It must not
    /// publish external state, replace laws or reinsert their friction heat.
    /// Source IDs must exactly match the supplied canonical link IDs.
    ///
    /// The original laws own their accepted memory and thermal applicability.
    /// Only irreversible contact work heats the cells; reversible shear and
    /// original network damping retain their mechanical accounts. A frozen
    /// elastic source must refuse unsupported temperature changes: variable
    /// elastic free energy needs its own reversible heat/entropy port, which
    /// this adapter does not manufacture by replacing parameters or history.
    ///
    /// Cell-energy guesses use sorted CellId order. Normal loads and frames are
    /// fixed per sample. Held-effort mechanics and backward-Euler heat require
    /// refinement; nonlinear convergence alone does not establish uniqueness.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_with_friction_heat<'m, P, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'m>,
        bindings: &[ModalTangentialHeatRoute<'_>],
        external: &[f64],
        drives: &[PlanarTangentialModalDrive],
        controls: MultiTangentialHeatSolve,
        energy_guess_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'m>,
            &[u64],
            &[f64],
        ) -> Result<MultiTangentialHeatSource<P>, E>,
        gate: &CancelGate,
    ) -> Result<HeatedMultiTangentialModalFrame<P>, MultiTangentialHeatError<E>> {
        let checkpoint = || {
            if gate.is_requested() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        let solve = controls.solve;
        solve.newton.validate()?;
        let sample = self.network.samples_rendered();
        let next_sample = sample
            .checked_add(1)
            .ok_or(MultiTangentialHeatError::Invalid(
                "modal sample clock overflow",
            ))?;
        let dt = self.network.sample_period_s();
        let start = sample as f64 * dt;
        let target = next_sample as f64 * dt;
        let contacts = self.laws.len();
        let cells = thermal.cells().len();
        let dimension = cells
            .checked_add(contacts)
            .ok_or(MultiTangentialHeatError::Invalid(
                "joint thermal coordinate count overflow",
            ))?;
        if thermal.time_s() != start
            || !target.is_finite()
            || target <= start
            || bindings.len() != contacts
            || drives.len() != contacts
            || contacts == 0
            || contacts > controls.maximum_interfaces
            || dimension > controls.maximum_unknowns
            || solve.maximum_cells == 0
            || cells > solve.maximum_cells
            || solve.maximum_trials == 0
            || solve.maximum_newton_iterations == 0
            || !solve.difference_step_j.is_finite()
            || solve.difference_step_j <= 0.0
            || !solve.temperature_scale_j_per_k.is_finite()
            || solve.temperature_scale_j_per_k <= 0.0
            || !solve.combined_energy_tolerance_j.is_finite()
            || solve.combined_energy_tolerance_j < 0.0
            || energy_guess_j.len() != cells
            || energy_guess_j.iter().any(|x| !x.is_finite())
            || drives
                .iter()
                .any(|drive| drive.contact_temperature_k.is_some())
        {
            return Err(MultiTangentialHeatError::Invalid(
                "synchronized bounded cells/interfaces and solved contact temperatures required",
            ));
        }
        let mut seen = vec![false; contacts];
        for binding in bindings {
            checkpoint()?;
            if binding.contact_index >= contacts
                || seen[binding.contact_index]
                || !binding.temperature_guess_k.is_finite()
                || binding.temperature_guess_k <= 0.0
            {
                return Err(MultiTangentialHeatError::Invalid(
                    "each mechanical contact requires one heat route and positive temperature guess",
                ));
            }
            seen[binding.contact_index] = true;
        }
        let mut ordered: Vec<_> = bindings.iter().collect();
        ordered.sort_by_key(|binding| binding.route.link_id());
        if ordered
            .windows(2)
            .any(|pair| pair[0].route.link_id() == pair[1].route.link_id())
        {
            return Err(MultiTangentialHeatError::Invalid(
                "each physical thermal link requires exactly one mechanical heat owner",
            ));
        }
        let routes: Vec<_> = ordered.iter().map(|binding| binding.route).collect();
        let contact_indices: Vec<_> = ordered
            .iter()
            .map(|binding| binding.contact_index)
            .collect();
        let step = FrictionHeatNetworkStep::new(
            thermal,
            target,
            &routes,
            controls.maximum_interfaces,
            checkpoint,
        )?;
        let trials = Cell::new(0usize);
        let owner = RefCell::new(&mut *self);
        let evaluate = |x: &[f64], with_entropy: bool| {
            checkpoint()?;
            if trials.get() >= solve.maximum_trials {
                return Err(MultiTangentialHeatError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            if x.len() != dimension || x.iter().any(|x| !x.is_finite()) {
                return Err(MultiTangentialHeatError::Invalid(
                    "complete finite simultaneous friction/heat coordinates required",
                ));
            }
            let temperatures: Vec<_> = ordered
                .iter()
                .enumerate()
                .map(|(i, binding)| {
                    binding.temperature_guess_k + x[cells + i] / solve.temperature_scale_j_per_k
                })
                .collect();
            let endpoint = |thermal: &DepositionSimulation<'m>, temperatures: &[f64]| {
                let MultiTangentialHeatSource {
                    interfaces,
                    transport,
                } = resolve(thermal, step.link_ids(), temperatures)
                    .map_err(MultiTangentialHeatError::Source)?;
                if interfaces.len() != contacts
                    || interfaces
                        .iter()
                        .zip(step.link_ids())
                        .any(|(source, id)| source.link_id != *id)
                {
                    return Err(MultiTangentialHeatError::Invalid(
                        "material conductivities must retain every canonical original link identity",
                    ));
                }
                let mut candidate_drives = drives.to_vec();
                for (i, &contact) in contact_indices.iter().enumerate() {
                    candidate_drives[contact].contact_temperature_k = Some(temperatures[i]);
                }
                let mechanical =
                    owner
                        .borrow_mut()
                        .stage_candidate(external, &candidate_drives, Some(gate))?;
                let interfaces = interfaces
                    .iter()
                    .enumerate()
                    .map(|(i, source)| InterfaceHeatInput {
                        dissipated_j: mechanical.frame.contacts[contact_indices[i]].dissipated_j,
                        conductivity_w_m_k: source.conductivity_w_m_k,
                    })
                    .collect();
                let StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload,
                } = transport;
                Ok(FrictionHeatNetworkInput {
                    interfaces,
                    transport: StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: Staged {
                            mechanical,
                            source: payload,
                        },
                    },
                })
            };
            if with_entropy {
                step.evaluate_with_entropy(&x[..cells], &temperatures, endpoint, checkpoint)
            } else {
                step.evaluate(&x[..cells], &temperatures, endpoint, checkpoint)
            }
        };
        let problem =
            FallibleFiniteDifference::new(dimension, solve.difference_step_j, |x: &[f64]| {
                evaluate(x, false).map(|trial| trial.residual_j().to_vec())
            })
            .map_err(numerical)?;
        let mut initial = energy_guess_j.to_vec();
        initial.resize(dimension, 0.0);
        let construction = NewtonKrylovState::new(&problem, initial, solve.newton);
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        let mut solver = construction?;
        let result = solver.run_controlled(&problem, solve.maximum_newton_iterations, &mut || {
            gate.is_requested() || problem.has_failed()
        });
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        checkpoint()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(MultiTangentialHeatError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(MultiTangentialHeatError::NotConverged(result.report));
        }
        // Final probes may leave body scratch at another temperature vector.
        // Restage all contacts and their summed body force at the accepted root.
        let trial = evaluate(&solver.x, controls.audit_entropy)?;
        let combined = trial.payload().mechanical.frame.energy_residual_j
            + trial.report().total_energy_residual_j;
        if !combined.is_finite() || combined.abs() > solve.combined_energy_tolerance_j {
            return Err(MultiTangentialHeatError::Invalid(
                "joint multi-contact mechanical and caloric work does not balance",
            ));
        }
        let (next_thermal, staged, heat) = trial.accept(solve.heat, checkpoint)?;
        checkpoint()?;
        drop(problem);
        let owner = owner.into_inner();
        owner.network.validate_staged_material_points(Some(gate))
            .map_err(MultiTangentialModalError::from)?;
        let mechanical = staged.mechanical.frame.clone();
        // No callback, root solve or fallible physical check follows publication.
        owner.publish_candidate(staged.mechanical);
        *thermal = next_thermal;
        Ok(HeatedMultiTangentialModalFrame {
            mechanical,
            heat,
            contact_indices,
            source: staged.source,
            nonlinear: result.report,
            trials: trials.get(),
            combined_energy_residual_j: combined,
        })
    }
}
