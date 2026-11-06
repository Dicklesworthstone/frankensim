//! Continuous native joint friction and the original finite thermal network.
//!
//! Every nonlinear trial starts from the accepted body momenta and caloric
//! state. The original constrained flow resolves arrest, holding and reversal;
//! its integrated friction work, not an endpoint impulse pairing, supplies the
//! heat. Conductivity and friction effort are resolved at the same candidate
//! temperatures. Shared cells store energy once and retain ordinary conduction.

/// Actual native pose evolution and thermal feedback on an authored time grid.
#[path = "joints/trajectory.rs"]
pub mod trajectory;

use super::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use super::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_mbd::rigid_rattle::distance::graph::ArticulatedJoint;
use fs_mbd::rigid_rattle::distance::graph::joint_friction::MAX_JOINT_FRICTION_PORTS;
use fs_mbd::rigid_rattle::distance::graph::joint_friction::relaxation::{
    JointFrictionFlowControls, JointFrictionFlowError, JointFrictionFlowPort,
    JointFrictionFlowReport, advance_joint_friction,
};
use fs_mbd::{MassProperties, RigidBodyState};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::Cell;

/// Source-authored correspondence between one axial loss port and thermal face.
///
/// A bearing or rail can have distributed load and heating. This binding is the
/// caller's reduced physical model, not an inferred contact point or bearing
/// load. One distinct original link belongs to each selected native joint.
#[derive(Clone, Copy)]
pub struct JointFrictionHeatRoute<'a> {
    pub joint_index: usize,
    pub route: &'a FrictionHeatRoute,
    /// Initial nonlinear coordinate only, never a thermostat.
    pub temperature_guess_k: f64,
}

/// Original source data at one candidate common interface temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JointFrictionHeatInterfaceSource {
    pub link_id: u64,
    /// Held physical effort and common static/kinetic bound for this phase.
    pub friction: JointFrictionFlowPort,
    /// Original materials in the route's A/B order.
    pub conductivity_w_m_k: [f64; 2],
    /// Fraction of irreversible mechanical loss converted to heat, in [0, 1].
    /// The remainder leaves this subsystem through a separately reported
    /// nonthermal work port. It is not an invented heat store or hidden loss.
    pub heat_fraction: f64,
}

pub struct JointFrictionHeatSource<P> {
    /// Exactly the callback's ascending original thermal-link order.
    pub interfaces: Vec<JointFrictionHeatInterfaceSource>,
    /// Other heat only. Friction heat and bound link conductances have this owner.
    pub transport: StepCoefficients<P>,
}

/// Independent native, nonlinear, resource and physical energy acceptance gates.
#[derive(Clone, Copy, Debug)]
pub struct JointFrictionHeatControls {
    pub mechanical: JointFrictionFlowControls,
    pub maximum_cells: usize,
    /// Cell increments plus one common temperature for each selected joint.
    pub maximum_unknowns: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Newton coordinate scale [J/K]; no physical capacitance is introduced.
    pub temperature_scale_j_per_k: f64,
    pub heat: FrictionHeatTolerances,
    /// Applied separately to loss allocation and to the combined energy balance.
    pub combined_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
    /// Original caloric entropy audit of the final candidate only. Nonthermal
    /// exports are outside its scope; this is not a whole-system entropy claim.
    pub audit_entropy: bool,
}

#[derive(Debug)]
pub enum JointFrictionHeatError<E> {
    Invalid(&'static str),
    Mechanical(JointFrictionFlowError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<JointFrictionFlowError> for JointFrictionHeatError<E> {
    fn from(error: JointFrictionFlowError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<DepositionError> for JointFrictionHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<NewtonError> for JointFrictionHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for JointFrictionHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native joint friction/heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for JointFrictionHeatError<E> {}

/// One source's physical loss, with its thermal and nonthermal destinations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JointFrictionHeatInterfaceReport {
    pub source: JointFrictionHeatInterfaceSource,
    /// Continuous friction path loss, including motion before a reversal [J].
    pub mechanical_loss_j: f64,
    pub heat_j: f64,
    /// Signed work convention: positive energy exported out of this subsystem.
    pub nonthermal_export_j: f64,
    /// Heat + exported nonthermal work - physical loss; retained, never heated.
    pub partition_residual_j: f64,
}

/// Complete private candidate. Publish the native bodies, thermal timeline and
/// source-owned nonthermal exports together, after any surrounding owners accept.
#[must_use]
pub struct HeatedJointFrictionCandidate<'m, P> {
    pub mechanical: JointFrictionFlowReport,
    pub thermal: DepositionSimulation<'m>,
    /// Original heat network; its dissipated_j is the heat share only.
    pub heat: FrictionHeatNetworkReport,
    /// Ascending original link order, independently of native joint order.
    pub interfaces: Vec<JointFrictionHeatInterfaceReport>,
    pub nonthermal_export_j: f64,
    /// Sum of all heat/exports minus the native owner's total physical loss.
    pub loss_partition_residual_j: f64,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub target_time_s: f64,
    /// Kinetic + caloric changes + nonthermal export - axial drive work -
    /// external heat - bilateral residual work. Numerical residual work is
    /// separately visible in mechanics and is never added to the heat sources.
    pub combined_energy_residual_j: f64,
}

struct Staged<P> {
    mechanical: JointFrictionFlowReport,
    interfaces: Vec<JointFrictionHeatInterfaceReport>,
    nonthermal_export_j: f64,
    loss_partition_residual_j: f64,
    source: P,
}

/// Propose one native joint-friction phase with simultaneous thermal feedback.
///
/// Poses, bilateral geometry and source efforts are held during each trial.
/// The common static/kinetic Coulomb law retains its exact arrest/reversal
/// segments; material/thermal continuation uses backward-Euler endpoint data.
/// Time and spatial refinement remain necessary. No bearing preload, pressure,
/// distinct static coefficient, elastic joint storage or pose advance is inferred.
/// The surrounding native integrator owns geometry updates between phases.
///
/// `accepted_time_s` must match the original thermal timeline. Both owners use
/// its represented target-minus-accepted duration, so clock rounding cannot
/// give the mechanical and thermal phases different elapsed times. Energy
/// guesses are in ascending CellId order. Bindings may be arbitrarily ordered;
/// the callback receives their ascending link IDs, corresponding original
/// joint indices, and solved candidate interface temperatures in the same order.
///
/// `resolve` must retain those identities and stage its source receipt without
/// publishing state. It owns material applicability, capacities, drives and
/// heat fractions. The nonthermal work port may feed source-owned wear/plastic
/// or other irreversible processes; their state is not synthesized here. Only
/// the reported heat share enters the original caloric cells. The existing
/// network solves the split between bodies from conductivity and temperature.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_joint_friction_heat<'m, P, E>(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    thermal: &DepositionSimulation<'m>,
    accepted_time_s: f64,
    duration_s: f64,
    bindings: &[JointFrictionHeatRoute<'_>],
    controls: JointFrictionHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        &[u64],
        &[usize],
        &[f64],
    ) -> Result<JointFrictionHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedJointFrictionCandidate<'m, P>, JointFrictionHeatError<E>> {
    let checkpoint = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    checkpoint()?;
    controls.newton.validate()?;
    let ports = bindings.len();
    let cells = thermal.cells().len();
    let dimension = cells
        .checked_add(ports)
        .ok_or(JointFrictionHeatError::Invalid(
            "thermal dimension overflow",
        ))?;
    let target_time_s = accepted_time_s + duration_s;
    if !accepted_time_s.is_finite()
        || accepted_time_s != thermal.time_s()
        || !duration_s.is_finite()
        || duration_s <= 0.0
        || !target_time_s.is_finite()
        || target_time_s <= accepted_time_s
        || states.len() != properties.len()
        || states.len() > controls.mechanical.graph.maximum_bodies
        || ports == 0
        || ports > MAX_JOINT_FRICTION_PORTS
        || cells == 0
        || cells > controls.maximum_cells
        || dimension > controls.maximum_unknowns
        || controls.maximum_trials == 0
        || controls.maximum_newton_iterations == 0
        || !controls.difference_step_j.is_finite()
        || controls.difference_step_j <= 0.0
        || !controls.temperature_scale_j_per_k.is_finite()
        || controls.temperature_scale_j_per_k <= 0.0
        || [
            controls.heat.cell_j,
            controls.heat.thermal_j,
            controls.heat.interface_j,
            controls.heat.total_j,
            controls.combined_energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        || energy_guess_j.len() != cells
        || energy_guess_j.iter().any(|x| !x.is_finite())
    {
        return Err(JointFrictionHeatError::Invalid(
            "synchronized bounded native/thermal phase and finite controls required",
        ));
    }
    for (i, binding) in bindings.iter().enumerate() {
        checkpoint()?;
        if binding.joint_index >= joints.len()
            || bindings[..i]
                .iter()
                .any(|previous| previous.joint_index == binding.joint_index)
            || !binding.temperature_guess_k.is_finite()
            || binding.temperature_guess_k <= 0.0
        {
            return Err(JointFrictionHeatError::Invalid(
                "distinct original joints and positive interface guesses required",
            ));
        }
    }
    let mut ordered: Vec<_> = bindings.iter().collect();
    ordered.sort_by_key(|binding| binding.route.link_id());
    if ordered
        .windows(2)
        .any(|pair| pair[0].route.link_id() == pair[1].route.link_id())
    {
        return Err(JointFrictionHeatError::Invalid(
            "each original thermal link needs one joint heat owner",
        ));
    }
    let routes: Vec<_> = ordered.iter().map(|binding| binding.route).collect();
    let joint_indices: Vec<_> = ordered.iter().map(|binding| binding.joint_index).collect();
    let step = FrictionHeatNetworkStep::new(
        thermal,
        target_time_s,
        &routes,
        MAX_JOINT_FRICTION_PORTS,
        checkpoint,
    )?;
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64], with_entropy: bool| {
        checkpoint()?;
        if trials.get() >= controls.maximum_trials {
            return Err(JointFrictionHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|x| !x.is_finite()) {
            return Err(JointFrictionHeatError::Invalid(
                "complete finite native joint/thermal iterate required",
            ));
        }
        let temperatures: Vec<_> = ordered
            .iter()
            .enumerate()
            .map(|(i, binding)| {
                binding.temperature_guess_k + x[cells + i] / controls.temperature_scale_j_per_k
            })
            .collect();
        let endpoint = |thermal: &DepositionSimulation<'m>, temperatures: &[f64]| {
            let JointFrictionHeatSource {
                interfaces,
                transport,
            } = resolve(thermal, step.link_ids(), &joint_indices, temperatures)
                .map_err(JointFrictionHeatError::Source)?;
            if interfaces.len() != ports
                || interfaces.iter().enumerate().any(|(i, source)| {
                    source.link_id != step.link_ids()[i]
                        || source.friction.joint_index() != joint_indices[i]
                        || !source.heat_fraction.is_finite()
                        || !(0.0..=1.0).contains(&source.heat_fraction)
                })
            {
                return Err(JointFrictionHeatError::Invalid(
                    "source must retain link/joint identities and a heat fraction in [0,1]",
                ));
            }
            let efforts: Vec<_> = interfaces.iter().map(|source| source.friction).collect();
            let mechanical = advance_joint_friction(
                states,
                properties,
                joints,
                &efforts,
                step.duration_s(),
                controls.mechanical,
                || gate.is_requested(),
            )?;
            let mut heat_inputs = Vec::with_capacity(ports);
            let mut reports = Vec::with_capacity(ports);
            let mut exported = Dd::ZERO;
            let mut converted = Dd::ZERO;
            for source in interfaces {
                checkpoint()?;
                let port = mechanical
                    .ports
                    .iter()
                    .find(|port| port.joint_index == source.friction.joint_index())
                    .ok_or(JointFrictionHeatError::Invalid("missing native joint loss"))?;
                let loss = -port.friction_work_j;
                let heat = loss * source.heat_fraction;
                let nonthermal = loss - heat;
                if !loss.is_finite()
                    || loss < 0.0
                    || !heat.is_finite()
                    || heat < 0.0
                    || !nonthermal.is_finite()
                    || nonthermal < 0.0
                    || (loss > 0.0 && source.heat_fraction > 0.0 && heat == 0.0)
                    || (loss > 0.0 && source.heat_fraction < 1.0 && nonthermal == 0.0)
                {
                    return Err(JointFrictionHeatError::Invalid(
                        "representable nonnegative physical loss allocation required",
                    ));
                }
                let residual =
                    (Dd::from_f64(heat) + Dd::from_f64(nonthermal) - Dd::from_f64(loss)).to_f64();
                if !residual.is_finite() || residual.abs() > controls.combined_energy_tolerance_j {
                    return Err(JointFrictionHeatError::Invalid(
                        "physical joint loss allocation does not balance",
                    ));
                }
                exported = exported + Dd::from_f64(nonthermal);
                converted = converted + Dd::from_f64(heat);
                heat_inputs.push(InterfaceHeatInput {
                    dissipated_j: heat,
                    conductivity_w_m_k: source.conductivity_w_m_k,
                });
                reports.push(JointFrictionHeatInterfaceReport {
                    source,
                    mechanical_loss_j: loss,
                    heat_j: heat,
                    nonthermal_export_j: nonthermal,
                    partition_residual_j: residual,
                });
            }
            let nonthermal_export_j = exported.to_f64();
            let loss_partition_residual_j =
                (converted + exported - Dd::from_f64(mechanical.dissipated_j)).to_f64();
            if !nonthermal_export_j.is_finite()
                || !loss_partition_residual_j.is_finite()
                || loss_partition_residual_j.abs() > controls.combined_energy_tolerance_j
            {
                return Err(JointFrictionHeatError::Invalid(
                    "total physical joint loss allocation does not balance",
                ));
            }
            let StepCoefficients {
                sources_w,
                conductances_w_k,
                payload,
            } = transport;
            Ok(FrictionHeatNetworkInput {
                interfaces: heat_inputs,
                transport: StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload: Staged {
                        mechanical,
                        interfaces: reports,
                        nonthermal_export_j,
                        loss_partition_residual_j,
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
        FallibleFiniteDifference::new(dimension, controls.difference_step_j, |x: &[f64]| {
            evaluate(x, false).map(|trial| trial.residual_j().to_vec())
        })
        .map_err(numerical)?;
    let mut initial = energy_guess_j.to_vec();
    initial.resize(dimension, 0.0);
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
        return Err(JointFrictionHeatError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(JointFrictionHeatError::NotConverged(result.report));
    }
    let trial = evaluate(&solver.x, controls.audit_entropy)?;
    let staged = trial.payload();
    let mechanical = &staged.mechanical;
    let heat = trial.report();
    let combined = (Dd::from_f64(mechanical.kinetic_energy_change_j)
        + Dd::from_f64(heat.thermal.stored_energy_change_j)
        + Dd::from_f64(staged.nonthermal_export_j)
        - Dd::from_f64(mechanical.drive_work_j)
        - Dd::from_f64(heat.external_heat_j)
        - Dd::from_f64(mechanical.constraint_work_j))
    .to_f64();
    if !combined.is_finite() || combined.abs() > controls.combined_energy_tolerance_j {
        return Err(JointFrictionHeatError::Invalid(
            "joint kinetic/caloric/exported work does not balance",
        ));
    }
    let (thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
    checkpoint()?;
    drop(problem);
    Ok(HeatedJointFrictionCandidate {
        mechanical: staged.mechanical,
        thermal,
        heat,
        interfaces: staged.interfaces,
        nonthermal_export_j: staged.nonthermal_export_j,
        loss_partition_residual_j: staged.loss_partition_residual_j,
        source: staged.source,
        nonlinear: result.report,
        trials: trials.get(),
        target_time_s,
        combined_energy_residual_j: combined,
    })
}

fn numerical<E>(
    error: FiniteDifferenceError<JointFrictionHeatError<E>>,
) -> JointFrictionHeatError<E> {
    match error {
        FiniteDifferenceError::Source(source) => source,
        FiniteDifferenceError::Invalid(message) => JointFrictionHeatError::Invalid(message),
    }
}
