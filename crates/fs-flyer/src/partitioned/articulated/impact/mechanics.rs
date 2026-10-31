//! Shared point-impact mechanics for independently assembled fluid-bearing trees.

use super::super::{ArticulatedState, MAX_ARTICULATED_LINKS};
use super::{finite, native, poll, refusal};
use crate::Refusal;
use crate::addedmass::joint::native_index;
use crate::addedmass::{AeroGeneralizedLoads, MAX_CONTROL_COORDS, PreparedArticulatedMass};
use crate::partitioned::MAX_CERT_STATES;
use fs_ga::{Se3, Twist, Vec3, Wrench};
use fs_mbd::articulated::contact_geometry::{NormalContactGeometry, PointDirection};
use fs_mbd::articulated::generalized::{GeneralizedKinematics, free_floating_terms};
use fs_mbd::articulated::{ArticulatedModel, BaseState, Kinematics, forward_kinematics};

/// One independently assembled dry/fluid tree participating in a common event.
#[derive(Clone, Copy, Debug)]
pub struct ArticulatedImpactSystem<'a> {
    pub model: &'a ArticulatedModel,
    pub state: &'a ArticulatedState,
}

/// Actual endpoint mechanics shared by normal and frictional contact laws.
#[derive(Clone, Debug, PartialEq)]
pub struct ImpactMechanics {
    /// Accepted states in original system order; configuration/memory unchanged.
    pub states: Vec<ArticulatedState>,
    /// Flyer `[linear, angular, joints]` rates, one vector per original tree.
    pub velocities_after: Vec<Vec<f64>>,
    /// Native material-point speeds in caller row order.
    pub relative_speeds_after: Vec<f64>,
    pub kinetic_before: f64,
    pub kinetic_after: f64,
    /// Stable quadratic polarization of the complete kinetic change [J].
    pub kinetic_change: f64,
    pub support_work: f64,
    pub energy_balance_residual: f64,
    pub linear_momentum_residual: Vec3,
    pub angular_momentum_residual: Vec3,
    pub endpoint_velocity_residual: f64,
}

/// A direction already prepared by the native material-point geometry owner.
#[derive(Clone, Copy)]
pub(crate) struct PointImpactRow {
    pub first: PointDirection,
    pub second: Option<PointDirection>,
    pub axis_world: Vec3,
    pub world_point: Vec3,
    pub plane_speed: f64,
    pub relative_speed_before: f64,
}
impl From<NormalContactGeometry> for PointImpactRow {
    fn from(row: NormalContactGeometry) -> Self {
        Self {
            first: row.first,
            second: row.second,
            axis_world: row.axis_world,
            world_point: row.world_point,
            plane_speed: row.plane_speed,
            relative_speed_before: row.relative_normal_speed,
        }
    }
}

struct SystemWorkspace<'a> {
    system: ArticulatedImpactSystem<'a>,
    mass: PreparedArticulatedMass,
    map: GeneralizedKinematics,
    velocity: Vec<f64>,
}
pub(crate) struct ImpactWorkspace<'a> {
    systems: Vec<SystemWorkspace<'a>>,
}
struct RowPart {
    system: usize,
    jacobian: Vec<f64>,
    response: Vec<f64>,
}
struct PreparedRow {
    geometry: PointImpactRow,
    parts: Vec<RowPart>,
}
pub(crate) struct PreparedRows {
    rows: Vec<PreparedRow>,
    pub delassus: Vec<f64>,
    pub relative_speeds_before: Vec<f64>,
}

fn vector(value: Vec3) -> Result<Vec3, Refusal> {
    finite(value.x)?;
    finite(value.y)?;
    finite(value.z)?;
    Ok(value)
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn velocity(state: &ArticulatedState) -> Vec<f64> {
    let native = state.motion.base.twist_body.to_array();
    let mut result: Vec<_> = (0..6).map(|i| native[native_index(i)]).collect();
    result.extend_from_slice(&state.motion.joint_velocity);
    result
}
fn xyz(value: &[f64]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}

fn spatial_scale(wrench: Wrench, pose: &Se3) -> Result<(f64, f64), Refusal> {
    let rotate = |value| {
        pose.rotation()
            .rotate(value)
            .map_err(|e| refusal("articulated-impact-invalid", e))
    };
    let force = rotate(wrench.force)?;
    let spin = rotate(wrench.torque)?;
    let orbital = pose.translation().cross(force);
    // A point at the world origin can have cancelling local and orbital
    // torque. Bound roundoff before that cancellation, not by its residual.
    Ok((finite(force.norm())?, finite(spin.norm() + orbital.norm())?))
}

impl<'a> ImpactWorkspace<'a> {
    pub(crate) fn new(
        systems: &[ArticulatedImpactSystem<'a>],
        mut evaluate_fluid: impl FnMut(
            usize,
            &GeneralizedKinematics,
        ) -> Result<AeroGeneralizedLoads, Refusal>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Refusal> {
        poll(cancelled)?;
        if systems.len() > 128 {
            return Err(refusal(
                "articulated-impact-invalid",
                "more than 128 participating trees",
            ));
        }
        let mut prepared = Vec::with_capacity(systems.len());
        for (index, &system) in systems.iter().enumerate() {
            poll(cancelled)?;
            if system.model.link_count() > MAX_ARTICULATED_LINKS
                || system.model.dof_count() > MAX_CONTROL_COORDS
                || system.state.memory.len() > MAX_CERT_STATES
            {
                return Err(refusal(
                    "articulated-impact-invalid",
                    "participating tree exceeds its size cap",
                ));
            }
            let motion = &system.state.motion;
            let dry = free_floating_terms(
                system.model,
                motion.base,
                &motion.joint_position,
                &motion.joint_velocity,
                &mut *cancelled,
            )
            .map_err(native)?;
            let fluid = evaluate_fluid(index, &dry.kinematics)?;
            poll(cancelled)?;
            if fluid
                .q_rigid_nonaccel
                .iter()
                .chain(&fluid.q_control_nonaccel)
                .any(|x| *x != 0.0)
            {
                return Err(refusal(
                    "coupled-load-ownership",
                    "instantaneous fluid inertia must not contain applied force loads",
                ));
            }
            let mass = PreparedArticulatedMass::new(&dry, &fluid)?;
            prepared.push(SystemWorkspace {
                system,
                mass,
                map: dry.kinematics,
                velocity: velocity(system.state),
            });
        }
        Ok(Self { systems: prepared })
    }

    pub(crate) fn kinematics(&self) -> Vec<&Kinematics> {
        self.systems
            .iter()
            .map(|s| s.map.native_kinematics())
            .collect()
    }

    pub(crate) fn project_rows(
        &self,
        geometries: &[PointImpactRow],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedRows, Refusal> {
        if geometries.len() > 64 {
            return Err(refusal(
                "articulated-impact-invalid",
                "more than 64 point-direction rows",
            ));
        }
        let mut rows = Vec::with_capacity(geometries.len());
        let mut speeds = Vec::with_capacity(geometries.len());
        for &geometry in geometries {
            poll(cancelled)?;
            finite(geometry.plane_speed)?;
            finite(geometry.relative_speed_before)?;
            vector(geometry.axis_world)?;
            vector(geometry.world_point)?;
            let mut parts: Vec<RowPart> = Vec::with_capacity(2);
            for term in core::iter::once(geometry.first).chain(geometry.second) {
                let system = self.systems.get(term.point.system_index).ok_or_else(|| {
                    refusal(
                        "articulated-impact-invalid",
                        "point system index out of range",
                    )
                })?;
                if term.point.link_index >= system.map.link_count() {
                    return Err(refusal(
                        "articulated-impact-invalid",
                        "point link index out of range",
                    ));
                }
                let part = if let Some(index) = parts
                    .iter()
                    .position(|p| p.system == term.point.system_index)
                {
                    &mut parts[index]
                } else {
                    parts.push(RowPart {
                        system: term.point.system_index,
                        jacobian: vec![0.0; system.mass.coordinate_count()],
                        response: Vec::new(),
                    });
                    parts.last_mut().expect("one inserted point direction")
                };
                let columns = system.map.body_jacobian(term.point.link_index);
                for (i, value) in part.jacobian.iter_mut().enumerate() {
                    *value = finite(*value + term.direction.pairing(columns[native_index(i)]))?;
                }
            }
            let mut speed = -geometry.plane_speed;
            for part in &mut parts {
                let system = &self.systems[part.system];
                speed = finite(speed + dot(&part.jacobian, &system.velocity))?;
                part.response = system.mass.solve(&part.jacobian)?;
            }
            speeds.push(speed);
            rows.push(PreparedRow { geometry, parts });
        }
        let count = rows.len();
        let mut delassus = vec![0.0; count * count];
        for i in 0..count {
            poll(cancelled)?;
            for j in 0..=i {
                let mut value = 0.0;
                // Only trees incident to BOTH rows contribute. No global
                // body-coordinate matrix or unrelated inverse action is built.
                for a in &rows[i].parts {
                    if let Some(b) = rows[j].parts.iter().find(|b| b.system == a.system) {
                        value = finite(
                            value
                                + 0.5 * dot(&a.jacobian, &b.response)
                                + 0.5 * dot(&b.jacobian, &a.response),
                        )?;
                    }
                }
                delassus[i * count + j] = value;
                delassus[j * count + i] = value;
            }
        }
        Ok(PreparedRows {
            rows,
            delassus,
            relative_speeds_before: speeds,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn apply(
        &self,
        rows: &PreparedRows,
        impulses: &[f64],
        predicted_speeds: &[f64],
        dissipated_energy: f64,
        velocity_tolerance: f64,
        energy_tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ImpactMechanics, Refusal> {
        poll(cancelled)?;
        if impulses.len() != rows.rows.len()
            || predicted_speeds.len() != rows.rows.len()
            || !impulses
                .iter()
                .chain(predicted_speeds)
                .all(|x| x.is_finite())
            || ![dissipated_energy, velocity_tolerance, energy_tolerance]
                .iter()
                .all(|x| x.is_finite() && *x >= 0.0)
        {
            return Err(refusal(
                "articulated-impact-invalid",
                "invalid impulse response or endpoint budgets",
            ));
        }
        let mut generalized: Vec<Vec<f64>> = self
            .systems
            .iter()
            .map(|s| vec![0.0; s.mass.coordinate_count()])
            .collect();
        let mut support_work = 0.0;
        let mut external_linear = Vec3::default();
        let mut external_angular = Vec3::default();
        let (mut linear_scale, mut angular_scale) = (0.0, 0.0);
        for (row, &impulse) in rows.rows.iter().zip(impulses) {
            poll(cancelled)?;
            for part in &row.parts {
                for (total, coefficient) in generalized[part.system].iter_mut().zip(&part.jacobian)
                {
                    *total = finite(*total + impulse * coefficient)?;
                }
            }
            let g = row.geometry;
            support_work = finite(support_work + impulse * g.plane_speed)?;
            for term in core::iter::once(g.first).chain(g.second) {
                let system = &self.systems[term.point.system_index];
                let (linear, angular) = spatial_scale(
                    term.direction,
                    &system.map.world_from_link(term.point.link_index),
                )?;
                linear_scale = finite(linear_scale + impulse.abs() * linear)?;
                angular_scale = finite(angular_scale + impulse.abs() * angular)?;
            }
            if g.second.is_none() {
                let force = vector(g.axis_world.scale(impulse))?;
                external_linear = vector(external_linear + force)?;
                external_angular = vector(external_angular + g.world_point.cross(force))?;
            }
        }
        let mut states = Vec::with_capacity(self.systems.len());
        let mut after_maps = Vec::with_capacity(self.systems.len());
        let mut velocities_after = Vec::with_capacity(self.systems.len());
        let (mut kinetic_before, mut kinetic_after, mut kinetic_change, mut work_scale) =
            (0.0, 0.0, 0.0, 0.0);
        let mut linear_change = Vec3::default();
        let mut angular_change = Vec3::default();
        for (system, impulse) in self.systems.iter().zip(&generalized) {
            poll(cancelled)?;
            let delta = system.mass.solve(impulse)?;
            let after: Vec<_> = system
                .velocity
                .iter()
                .zip(delta)
                .map(|(u, du)| finite(u + du))
                .collect::<Result<_, _>>()?;
            let mut state = system.system.state.clone();
            state.motion.base.twist_body = Twist::new(xyz(&after[3..6]), xyz(&after[..3]));
            state.motion.joint_velocity.copy_from_slice(&after[6..]);
            after_maps.push(
                forward_kinematics(
                    system.system.model,
                    BaseState::prescribed(
                        state.motion.base.world_from_base,
                        state.motion.base.twist_body,
                        Twist::zero(),
                    ),
                    &state.motion.joint_position,
                    &state.motion.joint_velocity,
                )
                .map_err(|e| refusal("articulated-impact-invalid", e))?,
            );
            kinetic_before =
                finite(kinetic_before + system.mass.kinetic_energy(&system.velocity)?)?;
            kinetic_after = finite(kinetic_after + system.mass.kinetic_energy(&after)?)?;
            let midpoint: Vec<_> = system
                .velocity
                .iter()
                .zip(&after)
                .map(|(a, b)| f64::midpoint(*a, *b))
                .collect();
            let momentum_midpoint = system.mass.momentum(&midpoint)?;
            let actual_delta: Vec<_> = after
                .iter()
                .zip(&system.velocity)
                .map(|(a, b)| a - b)
                .collect();
            kinetic_change = finite(kinetic_change + dot(&actual_delta, &momentum_midpoint))?;
            for (du, momentum) in actual_delta.iter().zip(&momentum_midpoint) {
                work_scale = finite(work_scale + (du * momentum).abs())?;
            }
            let dp = system.mass.momentum(&actual_delta)?;
            let local = Wrench::new(xyz(&dp[3..6]), xyz(&dp[..3]));
            let world = local.transform_by(&state.motion.base.world_from_base);
            linear_change = vector(linear_change + world.force)?;
            angular_change = vector(angular_change + world.torque)?;
            let (linear, angular) = spatial_scale(local, &state.motion.base.world_from_base)?;
            linear_scale = finite(linear_scale + linear)?;
            angular_scale = finite(angular_scale + angular)?;
            states.push(state);
            velocities_after.push(after);
        }
        let mut relative_speeds_after = Vec::with_capacity(rows.rows.len());
        let mut endpoint_velocity_residual: f64 = 0.0;
        for (row, &predicted) in rows.rows.iter().zip(predicted_speeds) {
            poll(cancelled)?;
            let mut speed = -row.geometry.plane_speed;
            let mut scale =
                row.geometry.plane_speed.abs() + row.geometry.relative_speed_before.abs();
            for term in core::iter::once(row.geometry.first).chain(row.geometry.second) {
                let twist = after_maps[term.point.system_index].body_twist[term.point.link_index];
                speed = finite(speed + term.direction.pairing(twist))?;
                let wrench = term.direction.to_array();
                let twist = twist.to_array();
                for (a, b) in wrench.iter().zip(twist) {
                    scale = finite(scale + (a * b).abs())?;
                }
            }
            let residual = finite((speed - predicted).abs())?;
            endpoint_velocity_residual = endpoint_velocity_residual.max(residual);
            if residual > finite(velocity_tolerance + 256.0 * f64::EPSILON * scale)? {
                return Err(refusal(
                    "articulated-impact-unresolved",
                    "native endpoint point speed does not match global contact law",
                ));
            }
            relative_speeds_after.push(speed);
        }
        let energy_balance_residual = finite(kinetic_change - support_work + dissipated_energy)?;
        if energy_balance_residual.abs()
            > finite(
                energy_tolerance
                    + 512.0 * f64::EPSILON * (work_scale + support_work.abs() + dissipated_energy),
            )?
        {
            return Err(refusal(
                "articulated-impact-unresolved",
                "complete dry/fluid kinetic work does not close",
            ));
        }
        let linear_momentum_residual = vector(linear_change - external_linear)?;
        let angular_momentum_residual = vector(angular_change - external_angular)?;
        if finite(linear_momentum_residual.norm())? > 2048.0 * f64::EPSILON * linear_scale
            || finite(angular_momentum_residual.norm())? > 2048.0 * f64::EPSILON * angular_scale
        {
            return Err(refusal(
                "articulated-impact-unresolved",
                "complete dry/fluid spatial impulse does not close",
            ));
        }
        poll(cancelled)?;
        Ok(ImpactMechanics {
            states,
            velocities_after,
            relative_speeds_after,
            kinetic_before,
            kinetic_after,
            kinetic_change,
            support_work,
            energy_balance_residual,
            linear_momentum_residual,
            angular_momentum_residual,
            endpoint_velocity_residual,
        })
    }
}
