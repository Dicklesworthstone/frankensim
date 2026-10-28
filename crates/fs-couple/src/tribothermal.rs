//! Closed dry-friction/enthalpy feedback at a prescribed sliding interface.
//!
//! The caller prescribes a positive normal load, a fixed nominal area and the
//! relative displacement of two surfaces in one unchanged tangent frame. There
//! is no normal displacement or normal work. The original fs-tribo law computes
//! tangential effort and physical heat; the original fs-thermal backward-Euler
//! residual advances finite caloric stores, conduction and authored boundaries.
//! Friction claims are re-evaluated at the SAME endpoint's two temperatures,
//! pressure, speed and optional accumulated sliding distance. Every Newton
//! trial starts from the same accepted slip history; only a balanced candidate
//! can publish both histories. Heat exported elsewhere remains in the balance.
//!
//! This entry is a fixed-footprint, isotropic, prescribed-motion reduction;
//! [`dynamics`] closes the same contact with finite-body tangent inertia.
//! No normal contact solve, inferred interface temperature, flash estimate,
//! wear geometry, lubricant, adhesion or friction entropy law is introduced.
//! Dry contact refuses liquid material states. Source curves remain their
//! original one-dimensional laws and validity domains; no multivariable law is
//! fabricated. Thermal discretization and source data require their own physical
//! justification. Endpoint coupling is first order; convergence is not a time
//! accuracy, uniqueness or experimental-validation certificate.

pub mod material;
/// Force-driven finite-body motion using the same material/thermal contact.
pub mod dynamics;
/// Native electric/magnetic field motion with material guide friction and heat.
pub mod field_motion;
pub use material::{EffusivityPartitionSample, MaterialSlidingContact, SemiInfiniteHeatPartition,
    SlidingHeatPartition, SlidingInterface, SlidingMaterialSample};

use fs_contact::interface_binding::InterfaceBindingError;
use fs_material::phase::PhaseStateError;
use fs_material::state_point::MaterialStatePointError;
use fs_material::thermal::source::ThermalMaterialError;
use fs_math::det;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, EnergyAudit};
use fs_tribo::{
    DissipationStep, DryTangentialState, TangentialLaw, TangentialSolveConfig, TangentialStep,
    TangentialStepInput, TriboError,
};
use std::cell::{Cell, RefCell};
use std::convert::Infallible;

/// Held load and actual A-minus-B displacement in the bound tangent frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlidingDrive {
    pub normal_load_n: f64,
    pub relative_increment_m: [f64; 2],
}

/// Explicit nonlinear work, effort and independent local/global energy gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlidingSolve {
    pub maximum_cells: usize,
    pub maximum_friction_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Numerical Jacobian perturbation in cell ENERGY, not temperature [J].
    pub difference_step_j: f64,
    pub cell_energy_tolerance_j: f64,
    pub thermal_energy_tolerance_j: f64,
    pub total_energy_tolerance_j: f64,
    pub tangential: TangentialSolveConfig,
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum SlidingError<E = Infallible> {
    Invalid(&'static str),
    Binding(InterfaceBindingError),
    Material(ThermalMaterialError),
    Phase(PhaseStateError),
    Tribology(TriboError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    FrictionBudget,
    MotionBudget,
    Cancelled,
}
impl<E> From<DepositionError> for SlidingError<E> {
    fn from(e: DepositionError) -> Self {
        if e == DepositionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Thermal(e)
        }
    }
}
impl<E> From<TriboError> for SlidingError<E> {
    fn from(e: TriboError) -> Self {
        Self::Tribology(e)
    }
}
impl<E> From<ThermalMaterialError> for SlidingError<E> {
    fn from(e: ThermalMaterialError) -> Self {
        Self::Material(e)
    }
}
impl<E> From<MaterialStatePointError> for SlidingError<E> {
    fn from(e: MaterialStatePointError) -> Self {
        Self::Material(ThermalMaterialError::Material(e))
    }
}
impl<E> From<PhaseStateError> for SlidingError<E> {
    fn from(e: PhaseStateError) -> Self {
        Self::Phase(e)
    }
}
impl<E> From<NewtonError> for SlidingError<E> {
    fn from(e: NewtonError) -> Self {
        Self::Newton(e)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SlidingError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "dry sliding/enthalpy coupling: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for SlidingError<E> {}

fn number<E>(value: f64) -> Result<f64, SlidingError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SlidingError::Invalid(
            "unrepresentable sliding coordinate, power or energy",
        ))
    }
}
fn numerical<E>(error: FiniteDifferenceError<SlidingError<E>>) -> SlidingError<E> {
    match error {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(s) => SlidingError::Invalid(s),
    }
}

/// Exact endpoint source, original force/work and finite-network accounts.
pub struct SlidingReport<S> {
    /// Effort into the interface. Friction acts on A with -traction and B with
    /// +traction; prescribed drives supply the opposite forces.
    pub friction: TangentialStep<DryTangentialState>,
    pub material: SlidingMaterialSample,
    pub partition: DissipationStep,
    /// Original thermal account, INCLUDING the internal friction deposition.
    pub thermal: EnergyAudit,
    pub source: S,
    pub nonlinear: NewtonReport,
    pub friction_trials: usize,
    pub maximum_cell_residual_j: f64,
    /// Independently reconstructed effort dotted with actual solved relative motion.
    pub drive_work_j: f64,
    pub deposited_heat_j: f64,
    /// Authored/callback thermal supply, independently evaluated WITHOUT friction.
    pub external_heat_j: f64,
    /// Thermal storage + interface storage + outside share - drive - other heat.
    pub total_energy_residual_j: f64,
}

pub struct SlidingCandidate<'a, 't, S> {
    thermal: DepositionSimulation<'t>,
    contact: MaterialSlidingContact<'a>,
    report: SlidingReport<S>,
}
impl<'a, 't, S> SlidingCandidate<'a, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn contact(&self) -> &MaterialSlidingContact<'a> {
        &self.contact
    }
    pub const fn report(&self) -> &SlidingReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialSlidingContact<'a>,
        SlidingReport<S>,
    ) {
        (self.thermal, self.contact, self.report)
    }
}

impl<'a> MaterialSlidingContact<'a> {
    /// Solve one interval against the original finite thermal network.
    ///
    /// The callback supplies only OTHER cell sources and existing-link material
    /// conductances at the proposed endpoint. It may use prepared material
    /// transport; it must not advance history or include friction heat. Its
    /// original payload and errors survive the nonlinear solve. Energy guesses
    /// use sorted CellId order. Cut at source/topology events first; target events
    /// remain pending in the thermal owner's accepted left-limit state.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: SlidingDrive,
        controls: SlidingSolve,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SlidingCandidate<'a, 't, S>, SlidingError<E>> {
        self.step_candidate_with_motion(thermal, target_s, drive, 0.0,
            drive.relative_increment_m, controls, guess_j, |state, _| resolve(state), cancelled)
    }

    /// Internal coupled proposal through the SAME friction and thermal owners.
    /// `drive` supplies released motion; compliance belongs to the mechanical
    /// owner. `material_increment_m` is its source-query iterate. Work, heat and
    /// accepted history always use the ORIGINAL return map's solved slip.
    /// A compliant caller MUST gate source-query versus solved motion before
    /// publishing this internal candidate. The callback sees that same native
    /// response and supplies only other heat/conduction and its own payload.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn step_candidate_with_motion<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: SlidingDrive,
        compliance_m_per_n: f64,
        material_increment_m: [f64; 2],
        controls: SlidingSolve,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>, &TangentialStep<DryTangentialState>)
            -> Result<StepCoefficients<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SlidingCandidate<'a, 't, S>, SlidingError<E>> {
        let cancel = RefCell::new(cancelled);
        let checkpoint = || {
            if (*cancel.borrow_mut())() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        controls.newton.validate()?;
        if self.time_s != thermal.time_s()
            || controls.maximum_cells == 0
            || thermal.cells().len() > controls.maximum_cells
            || guess_j.len() != thermal.cells().len()
            || controls.maximum_friction_trials == 0
            || controls.maximum_newton_iterations == 0
            || !compliance_m_per_n.is_finite() || compliance_m_per_n < 0.0
            || [
                controls.cell_energy_tolerance_j,
                controls.thermal_energy_tolerance_j,
                controls.total_energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(SlidingError::Invalid(
                "incompatible clocks, thermal unknowns or work/energy controls",
            ));
        }
        let step = BackwardEulerStep::new(thermal, target_s, checkpoint)?;
        let dt = step.duration_s();
        for value in drive.relative_increment_m.into_iter().chain(material_increment_m) {
            number::<E>(value)?;
        }
        let distance = number(det::hypot(
            material_increment_m[0],
            material_increment_m[1],
        ))?;
        let speed = number(distance / dt)?;
        let accumulated = number(self.distance_m + distance)?;
        if distance > 0.0 && (speed == 0.0 || accumulated == self.distance_m) {
            return Err(SlidingError::Invalid(
                "sliding motion vanished in speed or history conversion",
            ));
        }
        let trials = Cell::new(0usize);
        let evaluate = |x: &[f64]| {
            step.evaluate(
                x,
                |state| {
                    checkpoint()?;
                    if trials.get() >= controls.maximum_friction_trials {
                        return Err(SlidingError::FrictionBudget);
                    }
                    trials.set(trials.get() + 1);
                    let material = self.sample(state, drive.normal_load_n, speed, accumulated)?;
                    let friction = material.law.advance(
                        &self.history,
                        &TangentialStepInput {
                            normal_n: drive.normal_load_n,
                            slip_increment_m: drive.relative_increment_m,
                            dt_s: dt,
                            compliance_m_per_n,
                            contact_temperature_k: self
                                .data
                                .point
                                .mean_temperature_axis
                                .and_then(|axis| material.resolved.state_coordinate(axis)),
                            nominal_pressure_pa: Some(material.pressure_pa),
                            solve: controls.tangential,
                        },
                    )?;
                    let drive_work_j = number(
                        number(friction.traction_n[0] * friction.slip_increment_m[0])?
                            + number(friction.traction_n[1] * friction.slip_increment_m[1])?,
                    )?;
                    if compliance_m_per_n == 0.0
                        && friction.slip_increment_m != drive.relative_increment_m
                        || number(friction.stored_change_j + friction.dissipated_j - drive_work_j)?
                            .abs()
                            > controls.total_energy_tolerance_j
                    {
                        return Err(SlidingError::Invalid(
                            "original friction work does not close actual interface motion",
                        ));
                    }
                    let partition = friction.dissipation_step(material.partition)?;
                    let StepCoefficients {
                        mut sources_w,
                        conductances_w_k,
                        payload,
                    } = resolve(state, &friction).map_err(SlidingError::Source)?;
                    // Independent external supply: do not reconstruct it by subtracting
                    // internal heat from the final thermal balance.
                    let external_heat_j = number(
                        dt * state.external_step_power_w(
                            &sources_w,
                            &conductances_w_k,
                            checkpoint,
                        )?,
                    )?;
                    let heat = [partition.surface_a_heat_j(), partition.surface_b_heat_j()];
                    let mut deposited_heat_j = 0.0;
                    for (id, heat_j) in self.data.cells.into_iter().zip(heat) {
                        checkpoint()?;
                        let power = number(heat_j / dt)?;
                        if heat_j != 0.0 && power == 0.0 {
                            return Err(SlidingError::Invalid(
                                "friction heat vanished during power conversion",
                            ));
                        }
                        deposited_heat_j = number(deposited_heat_j + number(dt * power)?)?;
                        let total = sources_w.entry(id).or_insert(0.0);
                        *total = number(*total + power)?;
                    }
                    if number(
                        deposited_heat_j + partition.other_work_j() - partition.total_work_j(),
                    )?
                    .abs()
                        > controls.total_energy_tolerance_j
                    {
                        return Err(SlidingError::Invalid(
                            "friction power destinations do not close dissipative work",
                        ));
                    }
                    Ok::<_, SlidingError<E>>(StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: (
                            material,
                            friction,
                            partition,
                            payload,
                            drive_work_j,
                            deposited_heat_j,
                            external_heat_j,
                        ),
                    })
                },
                checkpoint,
            )
        };
        let problem = FallibleFiniteDifference::new(
            step.cell_ids().len(),
            controls.difference_step_j,
            |x: &[f64]| evaluate(x).map(|t| t.residual_j().to_vec()),
        )
        .map_err(numerical)?;
        let construction = NewtonKrylovState::new(&problem, guess_j.to_vec(), controls.newton);
        if let Some(e) = problem.take_failure() {
            return Err(numerical(e));
        }
        let mut solver = construction?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                (*cancel.borrow_mut())() || problem.has_failed()
            });
        if let Some(e) = problem.take_failure() {
            return Err(numerical(e));
        }
        checkpoint()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(SlidingError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(SlidingError::NotConverged(result.report));
        }
        let trial = evaluate(&solver.x)?;
        let maximum_cell_residual_j = trial.maximum_cell_residual_j();
        let (
            next,
            (
                material,
                friction,
                partition,
                source,
                drive_work_j,
                deposited_heat_j,
                external_heat_j,
            ),
            thermal_audit,
        ) = trial.accept(
            controls.cell_energy_tolerance_j,
            controls.thermal_energy_tolerance_j,
            checkpoint,
        )?;
        let total = number(
            thermal_audit.stored_energy_change_j
                + friction.stored_change_j
                + partition.other_work_j()
                - drive_work_j
                - external_heat_j,
        )?;
        if total.abs() > controls.total_energy_tolerance_j {
            return Err(SlidingError::Invalid(
                "finite heat, outside share and drive work do not balance",
            ));
        }
        // Publish one work increment and one original slip update only after
        // all endpoint residuals and the whole-interval account are admitted.
        let mut contact = self.clone();
        contact.work.record(partition)?;
        contact.history = friction.state.clone();
        let actual_distance = number(det::hypot(friction.slip_increment_m[0], friction.slip_increment_m[1]))?;
        contact.distance_m = number(self.distance_m + actual_distance)?;
        if actual_distance > 0.0 && contact.distance_m == self.distance_m {
            return Err(SlidingError::Invalid("actual sliding history increment is unrepresentable"));
        }
        contact.time_s = target_s;
        if let Some(sample) = &material.effusivity {
            contact.partition_maximum_diffusivity_m2_s = sample.maximum_diffusivity_m2_s;
        }
        checkpoint()?;
        Ok(SlidingCandidate {
            thermal: next,
            contact,
            report: SlidingReport {
                friction,
                material,
                partition,
                thermal: thermal_audit,
                source,
                nonlinear: result.report,
                friction_trials: trials.get(),
                maximum_cell_residual_j,
                drive_work_j,
                deposited_heat_j,
                external_heat_j,
                total_energy_residual_j: total,
            },
        })
    }
}

/// Finite caloric and entropy feedback for loaded native partial-slip history.
pub mod partial_slip;

/// Reciprocal finite-caloric normal contact with evolving indentation and patch.
pub mod normal;
