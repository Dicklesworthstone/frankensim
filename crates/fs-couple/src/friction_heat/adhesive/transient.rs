//! Original adhesive trajectories and finite caloric storage on one time step.
//!
//! Every thermal trial refreshes the same pinned junction-strength claims at
//! its solved face temperature and integrates the original native trajectory
//! through curved slip/spin and finite arrest. Only physical trajectory work
//! enters the original massless-face heat balance. Native body geometry and
//! mobility are prepared once, outside the thermal nonlinear iterations.
//!
//! Strengths are sampled at the thermal endpoint and held during that phase.
//! Temporal subdivision is required to resolve changing strength and heat-release
//! timing. Normal coefficients, solid area, poses and inertia remain fixed;
//! normal thermal storage, expansion and attachment changes have their own owners.

use super::validate_geometry;
use crate::friction_heat::{
    FrictionHeatInput, FrictionHeatReport, FrictionHeatRoute, FrictionHeatStep,
    FrictionHeatTolerances,
};
use fs_contact::interface_binding::adhesive_shear::spinning::temperature::{
    SourceThermalAdhesivePatch, SourceThermalAdhesivePatchError,
};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::adhesive::transient::{
    AdhesiveTransientControls, PreparedRigidAdhesiveTransient, RigidAdhesiveTransientCandidate,
    RigidAdhesiveTransientError, RigidAdhesiveTransientReport,
};
use fs_mbd::tangential_friction::adhesive::{RigidAdhesiveContact, RigidAdhesiveControls};
use fs_mbd::{MassProperties, RigidBodyState};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::Cell;

/// One original solid-contact disk and its existing thermal face, in A/B order.
#[derive(Clone, Copy)]
pub struct RigidAdhesiveTransientHeatContact<'a> {
    pub contact: RigidAdhesiveContact<'a>,
    pub temperature_axis: &'a str,
    pub route: &'a FrictionHeatRoute,
    /// Positive nonlinear seed; the physical face temperature is solved.
    pub temperature_guess_k: f64,
}

/// Bounded mechanics, thermal unknowns and original dimensional energy gates.
#[derive(Clone, Copy, Debug)]
pub struct RigidAdhesiveTransientHeatControls {
    pub mechanical: RigidAdhesiveControls,
    pub transient: AdhesiveTransientControls,
    pub maximum_cells: usize,
    pub maximum_unknowns: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Numerical coordinate scaling [J/K]; no interface heat capacity is added.
    pub temperature_scale_j_per_k: f64,
    pub heat: FrictionHeatTolerances,
    pub combined_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
    /// Audit exact original caloric entropy on the final root only.
    pub audit_entropy: bool,
}

/// Original endpoint bulk transport and other heat; friction has one owner.
pub struct RigidAdhesiveTransientHeatSource<P> {
    pub conductivity_w_m_k: [f64; 2],
    pub transport: StepCoefficients<P>,
}

/// The same accepted trajectory, source query and finite thermal endpoint.
pub struct HeatedRigidAdhesiveTransientReport<'a, P> {
    pub mechanical: RigidAdhesiveTransientReport,
    pub heat: FrictionHeatReport,
    pub friction_source: SourceThermalAdhesivePatch<'a>,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub target_time_s: f64,
    /// Actual delta(K + H_cells) minus external heat [J].
    pub combined_energy_residual_j: f64,
    /// Independently gated native mechanical and thermal residuals summed [J].
    pub summed_energy_residual_j: f64,
}

/// Immutable proposal; application publishes both original bodies and heat once.
pub struct HeatedRigidAdhesiveTransientCandidate<'a, 'm, P> {
    pub bodies_after: [RigidBodyState; 2],
    pub thermal: DepositionSimulation<'m>,
    pub report: HeatedRigidAdhesiveTransientReport<'a, P>,
}

#[derive(Debug)]
pub enum RigidAdhesiveTransientHeatError<E> {
    Invalid(&'static str),
    Mechanical(RigidAdhesiveTransientError),
    FrictionSource(SourceThermalAdhesivePatchError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<RigidAdhesiveTransientError> for RigidAdhesiveTransientHeatError<E> {
    fn from(error: RigidAdhesiveTransientError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<SourceThermalAdhesivePatchError> for RigidAdhesiveTransientHeatError<E> {
    fn from(error: SourceThermalAdhesivePatchError) -> Self {
        Self::FrictionSource(error)
    }
}
impl<E> From<DepositionError> for RigidAdhesiveTransientHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<NewtonError> for RigidAdhesiveTransientHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidAdhesiveTransientHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adhesive trajectory/heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for RigidAdhesiveTransientHeatError<E> {}

struct Staged<'a, P> {
    mechanical: RigidAdhesiveTransientCandidate,
    friction_source: SourceThermalAdhesivePatch<'a>,
    source: P,
}

/// Solve original finite-body motion through arrest with finite-cell heating.
///
/// Each nonlinear trial starts from the same accepted mechanical and caloric
/// state. The original source follows the solved common face temperature while
/// each bulk conductivity is resolved at that body's own endpoint state. The
/// callback returns other transport only and must stage external payloads.
///
/// A zero/tensile resultant normal force retains the actual adhesive solid area.
/// Attractive noncontact states cannot own this solid conducting route. The
/// existing thermal link is replaced inside the private heat trial, retaining
/// physical surface coverage and signed hot-to-cold exchange exactly once.
///
/// This endpoint source/heat discretization requires temporal refinement. The
/// mechanical stop bracket and travel bounds retain their native meaning; they
/// do not bound the thermal sampling error or resolve pointwise varying strength.
/// Geometry and normal storage remain with their original fixed owners.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_rigid_adhesive_transient_heat<'a, 'm, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidAdhesiveTransientHeatContact<'a>,
    thermal: &DepositionSimulation<'m>,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidAdhesiveTransientHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidAdhesiveTransientHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidAdhesiveTransientCandidate<'a, 'm, P>, RigidAdhesiveTransientHeatError<E>> {
    let checkpoint = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    checkpoint()?;
    controls.newton.validate()?;
    let cells = thermal.cells().len();
    let dimension = cells
        .checked_add(1)
        .ok_or(RigidAdhesiveTransientHeatError::Invalid(
            "thermal dimension overflow",
        ))?;
    let target = accepted_time_s + dt_s;
    if !accepted_time_s.is_finite()
        || accepted_time_s != thermal.time_s()
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !target.is_finite()
        || target <= accepted_time_s
        || cells == 0
        || cells > controls.maximum_cells
        || dimension > controls.maximum_unknowns
        || controls.maximum_trials == 0
        || controls.maximum_newton_iterations == 0
        || !controls.difference_step_j.is_finite()
        || controls.difference_step_j <= 0.0
        || !controls.temperature_scale_j_per_k.is_finite()
        || controls.temperature_scale_j_per_k <= 0.0
        || energy_guess_j.len() != cells
        || energy_guess_j.iter().any(|v| !v.is_finite())
        || !contact.temperature_guess_k.is_finite()
        || contact.temperature_guess_k <= 0.0
        || contact.contact.source.contact_area_m2() <= 0.0
        || [
            controls.heat.cell_j,
            controls.heat.thermal_j,
            controls.heat.interface_j,
            controls.heat.total_j,
            controls.combined_energy_tolerance_j,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v < 0.0)
    {
        return Err(RigidAdhesiveTransientHeatError::Invalid(
            "bounded attached disk and synchronized finite thermal coordinates required",
        ));
    }
    let prepared = PreparedRigidAdhesiveTransient::new(
        bodies,
        properties,
        contact.contact,
        target - accepted_time_s,
        controls.mechanical,
        controls.transient,
        gate,
    )?;
    let point = prepared.contact_points_world()[0];
    validate_geometry(
        [point.x, point.y, point.z],
        contact.contact,
        contact.route,
        controls.mechanical.geometry_tolerance_m,
    )
    .map_err(RigidAdhesiveTransientHeatError::Invalid)?;
    let step = FrictionHeatStep::new(thermal, target, contact.route, checkpoint)?;
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64], with_entropy: bool| {
        checkpoint()?;
        if trials.get() >= controls.maximum_trials {
            return Err(RigidAdhesiveTransientHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
            return Err(RigidAdhesiveTransientHeatError::Invalid(
                "complete finite thermal iterate required",
            ));
        }
        let temperature =
            contact.temperature_guess_k + x[cells] / controls.temperature_scale_j_per_k;
        let endpoint = |endpoint: &DepositionSimulation<'m>, temperature: f64| {
            checkpoint()?;
            let friction_source = contact
                .contact
                .source
                .at_temperature(contact.temperature_axis, temperature)?;
            let mechanical = prepared.advance_thermal(&friction_source, gate)?;
            let RigidAdhesiveTransientHeatSource {
                conductivity_w_m_k,
                transport,
            } = resolve(endpoint, temperature).map_err(RigidAdhesiveTransientHeatError::Source)?;
            let StepCoefficients {
                sources_w,
                conductances_w_k,
                payload,
            } = transport;
            Ok(FrictionHeatInput {
                dissipated_j: mechanical.report.dissipated_energy_j,
                conductivity_w_m_k,
                transport: StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload: Staged {
                        mechanical,
                        friction_source,
                        source: payload,
                    },
                },
            })
        };
        if with_entropy {
            step.evaluate_with_entropy(&x[..cells], temperature, endpoint, checkpoint)
        } else {
            step.evaluate(&x[..cells], temperature, endpoint, checkpoint)
        }
    };
    let problem =
        FallibleFiniteDifference::new(dimension, controls.difference_step_j, |x: &[f64]| {
            evaluate(x, false).map(|trial| trial.residual_j().to_vec())
        })
        .map_err(numerical)?;
    let mut initial = energy_guess_j.to_vec();
    initial.push(0.0);
    let construction = NewtonKrylovState::new(&problem, initial, controls.newton);
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    let mut solver = construction?;
    let result = solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
        gate.is_requested() || problem.has_failed()
    });
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    checkpoint()?;
    if result.progress != SolverRunProgress::Complete {
        return Err(RigidAdhesiveTransientHeatError::Interrupted(
            result.progress,
        ));
    }
    if !result.report.converged {
        return Err(RigidAdhesiveTransientHeatError::NotConverged(result.report));
    }
    let trial = evaluate(&solver.x, controls.audit_entropy)?;
    let native = &trial.payload().mechanical.report;
    let heat = trial.report();
    let combined =
        native.kinetic_energy_change_j + heat.thermal.stored_energy_change_j - heat.external_heat_j;
    let native_audit = native
        .audit
        .ok_or(RigidAdhesiveTransientHeatError::Invalid(
            "attached trajectory requires its native energy audit",
        ))?;
    let summed = native_audit.total_energy_residual_j + heat.total_energy_residual_j;
    if [combined, summed]
        .into_iter()
        .any(|v| !v.is_finite() || v.abs() > controls.combined_energy_tolerance_j)
    {
        return Err(RigidAdhesiveTransientHeatError::Invalid(
            "original trajectory kinetic and finite-cell energy does not balance",
        ));
    }
    let (thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
    checkpoint()?;
    Ok(HeatedRigidAdhesiveTransientCandidate {
        bodies_after: staged.mechanical.bodies_after,
        thermal,
        report: HeatedRigidAdhesiveTransientReport {
            mechanical: staged.mechanical.report,
            heat,
            friction_source: staged.friction_source,
            source: staged.source,
            nonlinear: result.report,
            trials: trials.get(),
            target_time_s: target,
            combined_energy_residual_j: combined,
            summed_energy_residual_j: summed,
        },
    })
}

/// Publish both finite bodies and the original thermal network together.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_adhesive_transient_heat<'a, 'm, P, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidAdhesiveTransientHeatContact<'a>,
    thermal: &mut DepositionSimulation<'m>,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidAdhesiveTransientHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidAdhesiveTransientHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidAdhesiveTransientReport<'a, P>, RigidAdhesiveTransientHeatError<E>> {
    let candidate = propose_rigid_adhesive_transient_heat(
        bodies,
        properties,
        contact,
        thermal,
        accepted_time_s,
        dt_s,
        controls,
        energy_guess_j,
        resolve,
        gate,
    )?;
    if gate.is_requested() {
        return Err(DepositionError::Cancelled.into());
    }
    bodies[contact.contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.contact.binding.body_b] = candidate.bodies_after[1];
    *thermal = candidate.thermal;
    Ok(candidate.report)
}

fn numerical<E>(
    error: FiniteDifferenceError<RigidAdhesiveTransientHeatError<E>>,
) -> RigidAdhesiveTransientHeatError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(field) => RigidAdhesiveTransientHeatError::Invalid(field),
    }
}
