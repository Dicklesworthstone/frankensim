//! Native articulated motion and finite-cell joint heating on an authored grid.
//!
//! Each interval advances the original constrained rigid owner for half a step,
//! solves the existing joint-friction/heat phase at that actual geometry, then
//! advances the rigid owner for the remaining half. Only native rigid steps
//! update poses. Thermal cells retain their source-authored reduced geometry.

use super::{
    JointFrictionHeatControls, JointFrictionHeatError, JointFrictionHeatRoute,
    JointFrictionHeatSource, propose_joint_friction_heat,
};
use crate::friction_heat::network::FrictionHeatNetworkReport;
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_mbd::rigid_rattle::distance::graph::joint_friction::MAX_JOINT_FRICTION_PORTS;
use fs_mbd::rigid_rattle::distance::graph::{
    ArticulatedJoint, DistanceGraphControls, DistanceGraphError, advance_articulated_graph,
};
use fs_mbd::{DynamicsError, MassProperties, RigidBodyIntegrator, RigidBodyState, Vec3, Wrench};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};

/// Native geometry gates are `heat.mechanical.graph` for both rigid half steps
/// and the original friction phase. Thermal/Newton gates remain independent.
#[derive(Clone, Copy, Debug)]
pub struct JointFrictionHeatTrajectoryControls {
    pub heat: JointFrictionHeatControls,
    pub maximum_step_s: f64,
    pub maximum_steps: usize,
    /// Total nonlinear trials across all accepted phases and the current trial.
    pub maximum_thermal_trials: usize,
    pub energy_residual_tolerance_j: f64,
    /// Sum of absolute complete-interval energy residuals, preventing cancelled
    /// time defects from passing only a signed final balance.
    pub absolute_energy_residual_tolerance_j: f64,
    /// Sum of absolute native half-step and friction-phase reaction work.
    pub constraint_work_tolerance_j: f64,
}

/// Integrated physical work at one selected source-defined bearing/rail route.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JointFrictionHeatTrajectoryPort {
    pub joint_index: usize,
    pub link_id: u64,
    pub drive_work_j: f64,
    pub mechanical_loss_j: f64,
    pub heat_j: f64,
    pub nonthermal_export_j: f64,
}

/// Only the final original owners and bounded work accounts are retained.
/// Publish these states and the reported nonthermal exports together.
#[must_use]
pub struct HeatedJointFrictionTrajectory<'m, P> {
    pub states_after: Vec<RigidBodyState>,
    pub thermal: DepositionSimulation<'m>,
    /// Ascending original link order, with cumulative work over the whole grid.
    pub ports: Vec<JointFrictionHeatTrajectoryPort>,
    /// Source and network receipts for the final interval, not cumulative heat.
    pub last_source: P,
    pub last_heat: FrictionHeatNetworkReport,
    pub start_time_s: f64,
    pub target_time_s: f64,
    pub steps: usize,
    pub thermal_trials: usize,
    /// Actual endpoint kinetic plus original uniform-gravity potential change.
    pub mechanical_energy_change_j: f64,
    /// Original finite-cell endpoint energy change, with no additional store.
    pub thermal_energy_change_j: f64,
    pub wrench_work_j: f64,
    pub drive_work_j: f64,
    pub mechanical_loss_j: f64,
    pub heat_j: f64,
    pub nonthermal_export_j: f64,
    pub external_heat_j: f64,
    pub constraint_work_j: f64,
    pub absolute_constraint_work_j: f64,
    /// Mechanical + caloric change + export - wrench/drive work - outside heat
    /// - measured bilateral reaction work. No residual is converted into heat.
    pub energy_residual_j: f64,
    pub absolute_step_energy_residual_j: f64,
    /// Maximum individual scalar-row closure and rate residuals from both
    /// original native half steps. Original vector-group gates also apply.
    pub maximum_constraint_gap_m: f64,
    pub maximum_constraint_rate_m_s: f64,
}

#[derive(Debug)]
pub enum JointFrictionHeatTrajectoryError<E> {
    Invalid(&'static str),
    Native(DistanceGraphError),
    Heat(JointFrictionHeatError<E>),
    Thermal(DepositionError),
    Limit(&'static str),
    Unresolved(&'static str),
    Allocation,
    Cancelled,
}
impl<E> From<DistanceGraphError> for JointFrictionHeatTrajectoryError<E> {
    fn from(error: DistanceGraphError) -> Self {
        Self::Native(error)
    }
}
impl<E> From<DynamicsError> for JointFrictionHeatTrajectoryError<E> {
    fn from(error: DynamicsError) -> Self {
        Self::Native(error.into())
    }
}
impl<E> From<JointFrictionHeatError<E>> for JointFrictionHeatTrajectoryError<E> {
    fn from(error: JointFrictionHeatError<E>) -> Self {
        Self::Heat(error)
    }
}
impl<E> From<DepositionError> for JointFrictionHeatTrajectoryError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for JointFrictionHeatTrajectoryError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native joint friction/heat trajectory: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for JointFrictionHeatTrajectoryError<E> {}

fn checkpoint<E>(gate: &CancelGate) -> Result<(), JointFrictionHeatTrajectoryError<E>> {
    if gate.is_requested() {
        Err(JointFrictionHeatTrajectoryError::Cancelled)
    } else {
        Ok(())
    }
}
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn dot(a: Vec3, b: Vec3) -> Dd {
    components(a)
        .into_iter()
        .zip(components(b))
        .fold(Dd::ZERO, |sum, (a, b)| {
            sum + Dd::from_f64(a) * Dd::from_f64(b)
        })
}

struct Smooth {
    states: Vec<RigidBodyState>,
    wrench_work: Dd,
    constraint_work: Dd,
    gap: f64,
    rate: f64,
}
#[allow(clippy::too_many_arguments)]
fn smooth<E>(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    integrator: RigidBodyIntegrator,
    wrenches: &[Wrench],
    duration_s: f64,
    controls: DistanceGraphControls,
    gate: &CancelGate,
) -> Result<Smooth, JointFrictionHeatTrajectoryError<E>> {
    checkpoint(gate)?;
    let out = advance_articulated_graph(
        states,
        properties,
        joints,
        integrator,
        wrenches,
        duration_s,
        controls,
        || gate.is_requested(),
    )?;
    let mut wrench_work = Dd::ZERO;
    for ((step, property), wrench) in out.step.free_steps.iter().zip(properties).zip(wrenches) {
        checkpoint(gate)?;
        let displacement_work = components(step.state_after.pose().position_world())
            .into_iter()
            .zip(components(step.state_before.pose().position_world()))
            .zip(components(wrench.force_world))
            .fold(Dd::ZERO, |sum, ((a, b), force)| {
                sum + (Dd::from_f64(a) - Dd::from_f64(b)) * Dd::from_f64(force)
            });
        let before =
            property.angular_velocity_body_checked(step.state_before.angular_momentum_body())?;
        let after =
            property.angular_velocity_body_checked(step.state_after.angular_momentum_body())?;
        let torque_work = (dot(wrench.torque_body, before) + dot(wrench.torque_body, after))
            * Dd::from_f64(0.5 * step.duration_seconds);
        wrench_work = wrench_work + displacement_work + torque_work;
    }
    if !wrench_work.to_f64().is_finite() {
        return Err(JointFrictionHeatTrajectoryError::Unresolved(
            "finite native wrench work",
        ));
    }
    let gap = out
        .step
        .residuals
        .iter()
        .fold(0.0_f64, |x, r| x.max(r.gap_m.abs()));
    let rate = out
        .step
        .residuals
        .iter()
        .fold(0.0_f64, |x, r| x.max(r.range_rate_mps.abs()));
    Ok(Smooth {
        states: out.step.states_after,
        wrench_work,
        constraint_work: Dd::from_f64(out.constraint_work_j),
        gap,
        rate,
    })
}

/// Paired native kinetic increments avoid subtracting large common kinetic
/// energies. Both angular velocities use the same original full inertia owner.
fn mechanical_change<E>(
    before: &[RigidBodyState],
    after: &[RigidBodyState],
    properties: &[MassProperties],
    integrator: RigidBodyIntegrator,
    gate: &CancelGate,
) -> Result<Dd, JointFrictionHeatTrajectoryError<E>> {
    let mut change = Dd::ZERO;
    for ((&a, &b), &p) in before.iter().zip(after).zip(properties) {
        checkpoint(gate)?;
        integrator.diagnostics(a, p)?;
        integrator.diagnostics(b, p)?;
        let velocities = [
            a.center_of_mass_velocity_world(p)?,
            b.center_of_mass_velocity_world(p)?,
        ];
        let spins = [
            p.angular_velocity_body_checked(a.angular_momentum_body())?,
            p.angular_velocity_body_checked(b.angular_momentum_body())?,
        ];
        for (rates, old, new) in [
            (
                velocities,
                a.linear_momentum_world(),
                b.linear_momentum_world(),
            ),
            (spins, a.angular_momentum_body(), b.angular_momentum_body()),
        ] {
            for (((old_rate, new_rate), old), new) in components(rates[0])
                .into_iter()
                .zip(components(rates[1]))
                .zip(components(old))
                .zip(components(new))
            {
                change = change
                    + (Dd::from_f64(0.5 * old_rate) + Dd::from_f64(0.5 * new_rate))
                        * (Dd::from_f64(new) - Dd::from_f64(old));
            }
        }
        for ((old, new), gravity) in components(a.pose().position_world())
            .into_iter()
            .zip(components(b.pose().position_world()))
            .zip(components(integrator.gravity().acceleration_world()))
        {
            change = change
                - Dd::from_f64(p.mass())
                    * Dd::from_f64(gravity)
                    * (Dd::from_f64(new) - Dd::from_f64(old));
        }
    }
    if !change.to_f64().is_finite() {
        return Err(JointFrictionHeatTrajectoryError::Unresolved(
            "finite native mechanical change",
        ));
    }
    Ok(change)
}
fn thermal_change<E>(
    before: &DepositionSimulation<'_>,
    after: &DepositionSimulation<'_>,
    gate: &CancelGate,
) -> Result<Dd, JointFrictionHeatTrajectoryError<E>> {
    if before.cells().len() != after.cells().len() {
        return Err(JointFrictionHeatTrajectoryError::Invalid(
            "thermal inventory changed",
        ));
    }
    let mut change = Dd::ZERO;
    for ((id, a), (next_id, b)) in before.cells().iter().zip(after.cells()) {
        checkpoint(gate)?;
        if id != next_id || a.mass_kg != b.mass_kg || a.material_id != b.material_id {
            return Err(JointFrictionHeatTrajectoryError::Invalid(
                "thermal inventory changed",
            ));
        }
        change = change
            + Dd::from_f64(a.mass_kg)
                * (Dd::from_f64(b.specific_h_j_kg) - Dd::from_f64(a.specific_h_j_kg));
    }
    if !change.to_f64().is_finite() {
        return Err(JointFrictionHeatTrajectoryError::Unresolved(
            "finite caloric change",
        ));
    }
    Ok(change)
}

/// Advance actual body poses, momenta and finite temperatures on a bounded grid.
///
/// `step_end_times_s` must be nonempty and strictly increase from the accepted
/// thermal time. Author known source/load switches as grid endpoints. The source
/// callback sees the physical midpoint time and the actual native midpoint
/// bodies, then each candidate thermal endpoint. It must stage data without
/// publishing external state, as in [`propose_joint_friction_heat`].
///
/// The two native half steps use the supplied fixed world forces/body torques
/// and original gravity. Supply axial loads that must balance static friction
/// through the source's joint `drive_*` port. Separate body-wrench kicks can
/// create O(h) artificial axial slip and heat in a nominally stuck balance.
/// Changing geometry and other split loads also require time refinement.
///
/// Heat uses backward Euler, so symmetric mechanical splitting does not make
/// this coupled method second order. A successful energy gate bounds numerical
/// work defects; it does not certify trajectory or heating accuracy on a coarse
/// grid. Within each phase the original friction owner resolves its complete
/// piecewise velocity path, but native pose travel converges to that path only
/// under refinement or when the grid resolves its slope changes.
///
/// Thermal routes remain the source-defined reduced bearing/rail interfaces.
/// Their original cell geometry, mass, materials and coverage stay fixed; no
/// moving footprint, material remapping or flash-temperature model is inferred.
/// Only the final states/receipts and at most six cumulative port accounts are
/// retained. Memory is bounded independently of grid length by original native
/// and thermal inventory limits, excluding callback-owned source payloads.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_joint_friction_heat_trajectory<'m, P, E>(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    integrator: RigidBodyIntegrator,
    wrenches: &[Wrench],
    thermal: &DepositionSimulation<'m>,
    step_end_times_s: &[f64],
    bindings: &[JointFrictionHeatRoute<'_>],
    controls: JointFrictionHeatTrajectoryControls,
    resolve: impl Fn(
        f64,
        &[RigidBodyState],
        &DepositionSimulation<'m>,
        &[u64],
        &[usize],
        &[f64],
    ) -> Result<JointFrictionHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedJointFrictionTrajectory<'m, P>, JointFrictionHeatTrajectoryError<E>> {
    checkpoint(gate)?;
    if states.len() != properties.len()
        || states.len() != wrenches.len()
        || states.is_empty()
        || states.len() > controls.heat.mechanical.graph.maximum_bodies
        || thermal.cells().len() > controls.heat.maximum_cells
        || bindings.is_empty()
        || bindings.len() > MAX_JOINT_FRICTION_PORTS
        || step_end_times_s.is_empty()
        || step_end_times_s.len() > controls.maximum_steps
        || controls.maximum_thermal_trials == 0
        || !controls.maximum_step_s.is_finite()
        || controls.maximum_step_s <= 0.0
        || [
            controls.energy_residual_tolerance_j,
            controls.absolute_energy_residual_tolerance_j,
            controls.constraint_work_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(JointFrictionHeatTrajectoryError::Invalid(
            "bounded nonempty joint/heat trajectory required",
        ));
    }
    let start_time_s = thermal.time_s();
    let mut time = start_time_s;
    for &end in step_end_times_s {
        checkpoint(gate)?;
        let dt = end - time;
        let middle = time + 0.5 * dt;
        if !end.is_finite()
            || !dt.is_finite()
            || dt <= 0.0
            || dt > controls.maximum_step_s
            || !(middle > time && middle < end)
        {
            return Err(JointFrictionHeatTrajectoryError::Invalid(
                "strictly increasing representable bounded time grid required",
            ));
        }
        time = end;
    }
    let target_time_s = time;
    let mut current_states = Vec::new();
    current_states
        .try_reserve_exact(states.len())
        .map_err(|_| JointFrictionHeatTrajectoryError::Allocation)?;
    current_states.extend_from_slice(states);
    let mut current_thermal = thermal.clone();
    let mut current_bindings = bindings.to_vec();
    current_bindings.sort_by_key(|binding| binding.route.link_id());
    let energy_guess = vec![0.0; thermal.cells().len()];
    // Per-port: drive, mechanical loss, converted heat, exported nonthermal work.
    let mut port_work = [[Dd::ZERO; 4]; MAX_JOINT_FRICTION_PORTS];
    let mut wrench = Dd::ZERO;
    let mut external_heat = Dd::ZERO;
    let mut constraint = Dd::ZERO;
    let mut absolute_constraint = Dd::ZERO;
    let mut absolute_step_residual = Dd::ZERO;
    let mut gap = 0.0_f64;
    let mut rate = 0.0_f64;
    let mut thermal_trials = 0usize;
    let mut last = None;
    time = start_time_s;
    for &end in step_end_times_s {
        checkpoint(gate)?;
        let remaining = controls.maximum_thermal_trials - thermal_trials;
        if remaining == 0 {
            return Err(JointFrictionHeatTrajectoryError::Limit(
                "total thermal trials",
            ));
        }
        let dt = end - time;
        let middle = time + 0.5 * dt;
        let first = smooth(
            &current_states,
            properties,
            joints,
            integrator,
            wrenches,
            middle - time,
            controls.heat.mechanical.graph,
            gate,
        )?;
        let mut phase_controls = controls.heat;
        phase_controls.maximum_trials = phase_controls.maximum_trials.min(remaining);
        let phase = propose_joint_friction_heat(
            &first.states,
            properties,
            joints,
            &current_thermal,
            time,
            dt,
            &current_bindings,
            phase_controls,
            &energy_guess,
            |state, links, indices, temperatures| {
                resolve(middle, &first.states, state, links, indices, temperatures)
            },
            gate,
        )?;
        if phase.target_time_s != end {
            return Err(JointFrictionHeatTrajectoryError::Invalid(
                "thermal clock differs from authored endpoint",
            ));
        }
        thermal_trials += phase.trials;
        let second = smooth(
            &phase.mechanical.states_after,
            properties,
            joints,
            integrator,
            wrenches,
            end - middle,
            controls.heat.mechanical.graph,
            gate,
        )?;
        let interval_wrench = first.wrench_work + second.wrench_work;
        let interval_constraint = first.constraint_work
            + second.constraint_work
            + Dd::from_f64(phase.mechanical.constraint_work_j);
        let step_residual = mechanical_change(
            &current_states,
            &second.states,
            properties,
            integrator,
            gate,
        )? + thermal_change(&current_thermal, &phase.thermal, gate)?
            + Dd::from_f64(phase.nonthermal_export_j)
            - interval_wrench
            - Dd::from_f64(phase.mechanical.drive_work_j)
            - Dd::from_f64(phase.heat.external_heat_j)
            - interval_constraint;
        absolute_step_residual =
            absolute_step_residual + Dd::from_f64(step_residual.to_f64().abs());
        absolute_constraint = absolute_constraint
            + Dd::from_f64(first.constraint_work.to_f64().abs())
            + Dd::from_f64(second.constraint_work.to_f64().abs())
            + Dd::from_f64(phase.mechanical.constraint_work_j.abs());
        if !absolute_step_residual.to_f64().is_finite()
            || absolute_step_residual.to_f64() > controls.absolute_energy_residual_tolerance_j
            || !absolute_constraint.to_f64().is_finite()
            || absolute_constraint.to_f64() > controls.constraint_work_tolerance_j
        {
            return Err(JointFrictionHeatTrajectoryError::Unresolved(
                "cumulative native work defect",
            ));
        }
        for (i, report) in phase.interfaces.iter().enumerate() {
            checkpoint(gate)?;
            let port = phase
                .mechanical
                .ports
                .iter()
                .find(|port| port.joint_index == report.source.friction.joint_index())
                .ok_or(JointFrictionHeatTrajectoryError::Invalid(
                    "missing native work port",
                ))?;
            for (sum, value) in port_work[i].iter_mut().zip([
                port.drive_work_j,
                report.mechanical_loss_j,
                report.heat_j,
                report.nonthermal_export_j,
            ]) {
                *sum = *sum + Dd::from_f64(value);
            }
            // A predictor for Newton only; the next temperature remains solved.
            // Equal receiving temperatures provide the exact zero-power
            // interface coordinate after arrest. Carrying the previous heated
            // interface peak instead can force Newton to manufacture a tiny,
            // unrepresentable caloric correction in a truly unchanged state.
            let cells = current_bindings[i].route.cells();
            let temperatures = [
                phase.thermal.temperature_k(cells[0])?,
                phase.thermal.temperature_k(cells[1])?,
            ];
            current_bindings[i].temperature_guess_k = if temperatures[0] == temperatures[1] {
                temperatures[0]
            } else {
                phase.heat.interfaces[i].interface_temperature_k
            };
        }
        wrench = wrench + interval_wrench;
        external_heat = external_heat + Dd::from_f64(phase.heat.external_heat_j);
        constraint = constraint + interval_constraint;
        gap = gap.max(first.gap).max(second.gap);
        rate = rate.max(first.rate).max(second.rate);
        current_states = second.states;
        current_thermal = phase.thermal;
        last = Some((phase.source, phase.heat));
        time = end;
    }
    let mechanical_change =
        mechanical_change(states, &current_states, properties, integrator, gate)?;
    let caloric_change = thermal_change(thermal, &current_thermal, gate)?;
    let mut total_work = [Dd::ZERO; 4];
    let ports: Vec<_> = current_bindings
        .iter()
        .enumerate()
        .map(|(i, binding)| {
            for (sum, value) in total_work.iter_mut().zip(port_work[i]) {
                *sum = *sum + value;
            }
            JointFrictionHeatTrajectoryPort {
                joint_index: binding.joint_index,
                link_id: binding.route.link_id(),
                drive_work_j: port_work[i][0].to_f64(),
                mechanical_loss_j: port_work[i][1].to_f64(),
                heat_j: port_work[i][2].to_f64(),
                nonthermal_export_j: port_work[i][3].to_f64(),
            }
        })
        .collect();
    let residual = (mechanical_change + caloric_change + total_work[3]
        - wrench
        - total_work[0]
        - external_heat
        - constraint)
        .to_f64();
    if !residual.is_finite()
        || residual.abs() > controls.energy_residual_tolerance_j
        || total_work.into_iter().any(|x| !x.to_f64().is_finite())
        || !wrench.to_f64().is_finite()
        || !constraint.to_f64().is_finite()
        || !external_heat.to_f64().is_finite()
    {
        return Err(JointFrictionHeatTrajectoryError::Unresolved(
            "complete native/caloric/exported work balance",
        ));
    }
    checkpoint(gate)?;
    let (last_source, last_heat) = last.ok_or(JointFrictionHeatTrajectoryError::Invalid(
        "empty trajectory",
    ))?;
    Ok(HeatedJointFrictionTrajectory {
        states_after: current_states,
        thermal: current_thermal,
        ports,
        last_source,
        last_heat,
        start_time_s,
        target_time_s,
        steps: step_end_times_s.len(),
        thermal_trials,
        mechanical_energy_change_j: mechanical_change.to_f64(),
        thermal_energy_change_j: caloric_change.to_f64(),
        wrench_work_j: wrench.to_f64(),
        drive_work_j: total_work[0].to_f64(),
        mechanical_loss_j: total_work[1].to_f64(),
        heat_j: total_work[2].to_f64(),
        nonthermal_export_j: total_work[3].to_f64(),
        external_heat_j: external_heat.to_f64(),
        constraint_work_j: constraint.to_f64(),
        absolute_constraint_work_j: absolute_constraint.to_f64(),
        energy_residual_j: residual,
        absolute_step_energy_residual_j: absolute_step_residual.to_f64(),
        maximum_constraint_gap_m: gap,
        maximum_constraint_rate_m_s: rate,
    })
}
