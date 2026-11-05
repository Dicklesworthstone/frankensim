//! Simultaneous axial friction on finite revolute and prismatic mechanisms.
//!
//! Let `J` contain the existing bilateral rows and `P` the selected friction
//! ports. A port impulse `p` needs the simultaneous original reactions
//! `r = -(J M^-1 J^T)^-1 J M^-1 P^T p`. The friction mobility is therefore
//! `P M^-1 P^T - (J M^-1 P^T)^T (J M^-1 J^T)^-1 (J M^-1 P^T)`.
//! Independent joint clamps lose those cross responses and generally leave a
//! supposedly arrested neighbour moving again.
//!
//! These are fixed-pose impulse events. Frozen capacities are caller-declared
//! physical inputs, not inferred static/kinetic phases, normal loads or a
//! timestep. Signed work at individual ports can exchange energy between
//! coupled modes; only the total mechanical loss is a passivity statement.

use super::revolute::friction::account;
use super::{
    ArticulatedJoint, DistanceGraphControls, DistanceGraphError, DistanceStepError, Geometry,
    JointReactionImpulse, MAX_DISTANCE_ROWS, Mobility, NativeReaction, PointPlaneJoint, Problem,
    Row, Rows, limit, planned_bytes_for, poll, reserved, zeros,
};
use crate::{
    DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, Vec3, Wrench,
};
use fs_contact::tangential::impact::multiple::{
    CoupledFrictionImpulseControls, CoupledFrictionImpulseError, CoupledFrictionImpulseInput,
    CoupledFrictionImpulseReport, solve_coupled_friction_impulses,
};
use fs_contact::tangential::impact::{
    CoulombImpulseError, CoulombImpulseInput, FrictionLaw, FrictionRegime,
    solve_isotropic_coulomb_impulse,
};
use fs_math::dd::Dd;

/// Each distinct supported joint consumes five of the original bilateral rows.
pub const MAX_JOINT_FRICTION_PORTS: usize = MAX_DISTANCE_ROWS / 5;

/// One existing joint and an explicitly supplied frozen impulse capacity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum JointFrictionPort {
    /// Positive torque impulse acts on the second body along the admitted mean
    /// hinge axis. The first body receives its negative.
    Revolute {
        /// Index of an original `ArticulatedJoint::Revolute`.
        joint_index: usize,
        /// Nonnegative angular impulse capacity [N m s].
        capacity_nms: f64,
    },
    /// Positive force impulse acts on the carriage (first body), along the
    /// carrier's positive joint-frame Z axis. The carrier receives its negative.
    Prismatic {
        /// Index of an original `ArticulatedJoint::Prismatic`.
        joint_index: usize,
        /// Nonnegative linear impulse capacity [N s].
        capacity_ns: f64,
    },
}
impl JointFrictionPort {
    /// Index in the original authored graph. A joint may be selected only once.
    #[must_use]
    pub const fn joint_index(self) -> usize {
        match self {
            Self::Revolute { joint_index, .. } | Self::Prismatic { joint_index, .. } => joint_index,
        }
    }
    fn capacity(self) -> f64 {
        match self {
            Self::Revolute { capacity_nms, .. } => capacity_nms,
            Self::Prismatic { capacity_ns, .. } => capacity_ns,
        }
    }
}

/// Explicit single-slider phase, in physical force-impulse units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PrismaticFrictionInput {
    /// A moving slider uses only its kinetic capacity, including when it stops.
    Sliding {
        /// Nonnegative kinetic capacity [N s].
        kinetic_capacity_ns: f64,
    },
    /// Requires exactly zero supplied axial carriage speed. Drive and friction
    /// act together; a held drive has exactly zero resultant impulse and work.
    RestingDrive {
        /// Signed drive on the carriage along the carrier's positive Z axis [N s].
        drive_impulse_ns: f64,
        /// Nonnegative static capacity, at least the kinetic capacity [N s].
        static_capacity_ns: f64,
        /// Nonnegative kinetic capacity after breakaway [N s].
        kinetic_capacity_ns: f64,
    },
}

/// Original graph admission plus impulse, work and actual native momentum gates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JointFrictionControls {
    /// Original geometry, bilateral-reaction, pivot, body/row and memory bounds.
    /// Time-integration iteration and travel fields are unused by this event.
    pub graph: DistanceGraphControls,
    /// Shared contact-kernel controls in work-conjugate Ns and m/s coordinates.
    /// A hinge uses `p = angular_impulse/L` and `v = L*relative_spin`, where
    /// `L` is that original hinge's authored angular length. These same physical
    /// velocity, impulse and support-gap gates audit the actual native endpoint.
    /// For the analytical single-slider phase, the iteration fields are unused.
    pub solve: CoupledFrictionImpulseControls,
    /// Absolute native work closure and bilateral residual-work allowance [J].
    pub work_tolerance_j: f64,
    /// Actual complete internal linear-momentum increment allowance [N s].
    pub linear_impulse_tolerance_ns: f64,
    /// Actual spin-plus-orbital angular-momentum increment allowance [N m s].
    pub angular_impulse_tolerance_nms: f64,
}

/// Meaning of the explicitly selected event; no numerical speed capture band.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JointFrictionPhase {
    /// Several original ports solved together under their declared fixed sets.
    FrozenCapacities,
    /// A single initially resting slider held its drive without motion or work.
    Held,
    /// A single slider retained motion or broke away under its declared drive.
    Sliding,
    /// Kinetic capacity arrested a single slider. Its actual residual is retained.
    Arrested,
}

/// Actual native motion and physical impulses at one joint port.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum JointFrictionMotion {
    /// Axial relative spin, second body minus first body.
    Revolute {
        /// Original authored length used to scale the port [m].
        angular_length_m: f64,
        /// Actual initial relative spin [rad/s].
        rate_before_rad_s: f64,
        /// Actual final relative spin [rad/s], without zero snapping.
        rate_after_rad_s: f64,
        /// Signed friction impulse about the world axis on body two [N m s].
        friction_impulse_nms: f64,
    },
    /// Carriage speed along the rotating rail and force impulses on the carriage.
    Prismatic {
        /// Initial derivative of the actual moving-rail axial coordinate [m/s].
        velocity_before_m_s: f64,
        /// Final derivative of the same coordinate [m/s], without zero snapping.
        velocity_after_m_s: f64,
        /// Signed friction impulse on the carriage [N s].
        friction_impulse_ns: f64,
        /// Authored drive impulse, zero for the frozen-capacity operation [N s].
        drive_impulse_ns: f64,
    },
}

/// One port's native action and its simultaneous, work-conjugate accounting.
#[derive(Debug, Clone, PartialEq)]
pub struct JointFrictionPortReport {
    /// Original authored joint index. Reports are in increasing index order.
    pub joint_index: usize,
    /// Unit hinge axis or the carrier's positive joint-frame Z axis.
    pub axis_world: Vec3,
    /// Native physical rate and force/couple impulses with explicit units.
    pub motion: JointFrictionMotion,
    /// Original ordered point-pair or pure-couple momentum application. Its
    /// stage-local work allocation differs from the simultaneous port work below.
    pub native_axial_impulse: JointReactionImpulse,
    /// Actual scaled native endpoint minus the complete mobility prediction [m/s].
    pub scaled_velocity_residual_m_s: f64,
    /// Signed drive work against the original/final native midpoint rate [J].
    pub drive_work_j: f64,
    /// Signed midpoint friction work [J]. An individual coupled port can be
    /// positive; only the aggregate event has a mechanical-loss guarantee.
    pub friction_work_j: f64,
}

/// Atomic original-body candidate, preserving every native work and momentum gate.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct JointFrictionReport {
    /// Original body order, with poses unchanged and momenta updated natively.
    pub states_after: Vec<RigidBodyState>,
    /// Explicit frozen-set event or the selected single-slider phase.
    pub phase: JointFrictionPhase,
    /// Canonical increasing authored joint-index order, independent of port order.
    pub ports: Vec<JointFrictionPortReport>,
    /// Complete constrained mobility in the returned port order [1/kg]. Only
    /// the first `ports.len()` rows/columns are active; the remainder is zero.
    pub scaled_mobility_per_kg: [[f64; MAX_JOINT_FRICTION_PORTS]; MAX_JOINT_FRICTION_PORTS],
    /// Shared box projection for frozen capacities; the exact scalar phase has
    /// no iterative projection receipt and returns `None`.
    pub projection: Option<CoupledFrictionImpulseReport>,
    /// Bilateral scalar impulses in original authored graph-row order [N s].
    pub constraint_impulses_ns: Vec<f64>,
    /// Original native point/couple receipts for every bilateral scalar row.
    pub native_constraint_impulses: Vec<JointReactionImpulse>,
    /// Sum of signed drive work [J].
    pub drive_work_j: f64,
    /// Nonpositive aggregate midpoint friction work [J]; not a resolved heat share.
    pub friction_work_j: f64,
    /// Bilateral residual work [J], separately bounded and never thermalized.
    pub constraint_work_j: f64,
    /// Reconstructed from actual native momentum differences [J].
    pub kinetic_energy_change_j: f64,
    /// Native kinetic change minus drive, friction and bilateral work [J].
    pub work_energy_residual_j: f64,
    /// Box support gap evaluated at actual native final velocities [J].
    pub native_duality_gap_j: f64,
    /// Actual complete linear-momentum increment [N s].
    pub linear_momentum_change_ns: Vec3,
    /// Actual complete angular-momentum increment about the world origin [N m s].
    pub angular_momentum_change_nms: Vec3,
    /// Admitted logical heap payload, excluding the fixed bounded stack banks.
    pub planned_live_bytes: usize,
}

/// Every refusal leaves all borrowed native states and graph inputs unchanged.
#[derive(Debug, Clone, PartialEq)]
pub enum JointFrictionError {
    /// Original native graph, geometry, state, resource or cancellation refusal.
    Native(DistanceGraphError),
    /// Shared coupled friction projection refusal.
    Coupled(CoupledFrictionImpulseError),
    /// Original scalar Coulomb projection refusal for the single-slider phase.
    Scalar(CoulombImpulseError),
    /// Invalid physical input or control.
    Invalid(&'static str),
    /// Static holding cannot be used on an already moving rail coordinate.
    NotResting {
        /// Selected original joint index.
        joint_index: usize,
        /// Actual nonzero input axial speed [m/s].
        velocity_m_s: f64,
    },
    /// The actual native system failed a named residual or physical gate.
    Unresolved(&'static str),
}
impl From<DistanceGraphError> for JointFrictionError {
    fn from(error: DistanceGraphError) -> Self {
        Self::Native(error)
    }
}
impl From<DistanceStepError> for JointFrictionError {
    fn from(error: DistanceStepError) -> Self {
        Self::Native(error.into())
    }
}
impl From<DynamicsError> for JointFrictionError {
    fn from(error: DynamicsError) -> Self {
        Self::Native(error.into())
    }
}
impl From<CoupledFrictionImpulseError> for JointFrictionError {
    fn from(error: CoupledFrictionImpulseError) -> Self {
        Self::Coupled(error)
    }
}
impl core::fmt::Display for JointFrictionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native joint friction: {self:?}")
    }
}
impl std::error::Error for JointFrictionError {}

/// Apply all selected friction ports simultaneously through the original graph.
///
/// The complete cross mobility is solved once. Native impulses are then applied
/// in canonical joint-index order and all original bilateral rows react. Port
/// capacities are explicit frozen physical sets; static/kinetic evolution,
/// continuous trajectories, normal-load inference and heat partition are not
/// supplied by this fixed-pose operation. Duplicate ports, locked/dependent
/// friction directions and unresolved physical responses refuse atomically.
///
/// The graph retains its original bounded dense factorization. At most six
/// selected hinges/sliders fit its 32 bilateral rows; all port matrices and
/// cross-response banks use fixed-size stack storage. Empty port sets retain
/// the original admitted states without a mobility factorization.
pub fn apply_joint_friction_impulses(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    ports: &[JointFrictionPort],
    controls: JointFrictionControls,
    cancelled: impl FnMut() -> bool,
) -> Result<JointFrictionReport, JointFrictionError> {
    run(states, properties, joints, ports, controls, None, cancelled)
}

/// Apply an explicitly selected kinetic or resting-drive phase to one slider.
///
/// This uses the same full graph response, actual rail-coordinate derivative,
/// native application and accounting as the coupled operation. Static holding
/// requires an exactly resting coordinate; a kinetic impulse never borrows its
/// larger static capacity to stop moving bodies. All other joints keep their
/// supplied ideal bilateral laws and participate in the finite recoil.
#[allow(clippy::too_many_arguments)]
pub fn apply_prismatic_friction_impulse(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    joint_index: usize,
    input: PrismaticFrictionInput,
    controls: JointFrictionControls,
    cancelled: impl FnMut() -> bool,
) -> Result<JointFrictionReport, JointFrictionError> {
    let capacity = match input {
        PrismaticFrictionInput::Sliding {
            kinetic_capacity_ns,
        } => kinetic_capacity_ns,
        PrismaticFrictionInput::RestingDrive {
            drive_impulse_ns,
            static_capacity_ns,
            kinetic_capacity_ns,
        } => {
            if !drive_impulse_ns.is_finite()
                || !static_capacity_ns.is_finite()
                || !kinetic_capacity_ns.is_finite()
                || static_capacity_ns < kinetic_capacity_ns
                || kinetic_capacity_ns < 0.0
            {
                return Err(JointFrictionError::Invalid("single-slider resting drive"));
            }
            if drive_impulse_ns.abs() <= static_capacity_ns {
                static_capacity_ns
            } else {
                kinetic_capacity_ns
            }
        }
    };
    run(
        states,
        properties,
        joints,
        &[JointFrictionPort::Prismatic {
            joint_index,
            capacity_ns: capacity,
        }],
        controls,
        Some(input),
        cancelled,
    )
}

#[derive(Clone, Copy)]
struct Port {
    source: JointFrictionPort,
    bodies: [usize; 2],
    axis: Vec3,
    length: f64,
    geometry: Geometry,
    axial_plane: Option<Row>,
    rate_before: f64,
    capacity: f64,
    linear_response: [Vec3; 2],
    spin_response: [Vec3; 2],
}
impl Port {
    const EMPTY: Self = Self {
        source: JointFrictionPort::Revolute {
            joint_index: 0,
            capacity_nms: 0.0,
        },
        bodies: [0; 2],
        axis: Vec3::ZERO,
        length: 1.0,
        geometry: Geometry {
            normal: Vec3::ZERO,
            arms_body_m: [Vec3::ZERO; 2],
            gap: 0.0,
            rate: 0.0,
        },
        axial_plane: None,
        rate_before: 0.0,
        capacity: 0.0,
        linear_response: [Vec3::ZERO; 2],
        spin_response: [Vec3::ZERO; 2],
    };

    fn angular(self) -> bool {
        matches!(self.source, JointFrictionPort::Revolute { .. })
    }

    fn prepare(
        source: JointFrictionPort,
        joints: &[ArticulatedJoint],
        problem: &Problem<'_, &[Wrench]>,
    ) -> Result<Self, JointFrictionError> {
        let states = problem.states;
        let properties = problem.properties;
        let mut out = Self {
            source,
            ..Self::EMPTY
        };
        match (source, joints.get(source.joint_index()).copied()) {
            (JointFrictionPort::Revolute { .. }, Some(ArticulatedJoint::Revolute(joint))) => {
                out.bodies = joint.bodies();
                let axes = core::array::from_fn::<_, 2, _>(|side| {
                    states[out.bodies[side]]
                        .pose()
                        .orientation()
                        .rotate_body_to_world(joint.axes_body()[side])
                });
                out.axis = axes[0]
                    .add(axes[1])
                    .normalized("joint-friction hinge axis")?;
                out.length = joint.angular_length_m();
                out.geometry.normal = out.axis.scale(out.length);
                let (rate, scaled) = out.velocity(states, properties, problem)?;
                out.rate_before = rate;
                out.geometry.rate = scaled;
                out.capacity = scaled_capacity(source.capacity(), out.length)?;
            }
            (JointFrictionPort::Prismatic { .. }, Some(ArticulatedJoint::Prismatic(joint))) => {
                out.bodies = joint.bodies();
                let points = joint.points_body_m();
                let axis_body =
                    joint.frames_body()[1].rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
                // The free coordinate is the moving rail's signed axial
                // point/plane separation. Reuse its exact rotational Jacobian
                // and reciprocal current-point force, not a frozen anchor arm.
                let row =
                    PointPlaneJoint::new(out.bodies, points[0], points[1], axis_body)?.row(0, 1);
                out.geometry = problem.geometry(row, out.bodies.map(|b| states[b]), None)?;
                out.axis = out.geometry.normal.scale(-1.0);
                out.axial_plane = Some(row);
                out.rate_before = out.geometry.rate;
                out.capacity = source.capacity();
            }
            _ => return Err(JointFrictionError::Invalid("selected joint kind and index")),
        }
        if !out.geometry.normal.is_finite() || !out.geometry.rate.is_finite() {
            return Err(JointFrictionError::Invalid("finite axial port"));
        }
        for side in 0..2 {
            let body = out.bodies[side];
            let sign = if side == 0 { -1.0 } else { 1.0 };
            let orientation = states[body].pose().orientation();
            let direction = orientation.rotate_world_to_body(out.geometry.normal.scale(sign));
            let couple = if out.angular() {
                direction
            } else {
                out.geometry.arms_body_m[side].cross(direction)
            };
            out.spin_response[side] = orientation
                .rotate_body_to_world(properties[body].angular_velocity_body_checked(couple)?);
            if !out.angular() {
                out.linear_response[side] =
                    out.geometry.normal.scale(sign / properties[body].mass());
            }
            if !out.spin_response[side].is_finite() || !out.linear_response[side].is_finite() {
                return Err(JointFrictionError::Invalid("finite axial response"));
            }
        }
        Ok(out)
    }

    // Physical rate and its exact native work-conjugate scaled coordinate.
    fn velocity(
        self,
        states: &[RigidBodyState],
        properties: &[MassProperties],
        problem: &Problem<'_, &[Wrench]>,
    ) -> Result<(f64, f64), JointFrictionError> {
        if let Some(row) = self.axial_plane {
            let g = problem.geometry(row, self.bodies.map(|b| states[b]), None)?;
            return Ok((g.rate, g.rate));
        }
        let spin = [
            states[self.bodies[0]]
                .point_kinematics(properties[self.bodies[0]], Vec3::ZERO)?
                .angular_velocity_world,
            states[self.bodies[1]]
                .point_kinematics(properties[self.bodies[1]], Vec3::ZERO)?
                .angular_velocity_world,
        ];
        let relative = difference(spin[1], spin[0]);
        let rate = project(relative, self.axis).to_f64();
        let scaled = project(relative, self.geometry.normal).to_f64();
        if !rate.is_finite() || !scaled.is_finite() || (rate != 0.0 && scaled == 0.0) {
            return Err(JointFrictionError::Invalid(
                "representable scaled axial rate",
            ));
        }
        Ok((rate, scaled))
    }

    fn action(
        self,
        states: &[RigidBodyState],
        properties: &[MassProperties],
        impulse: f64,
    ) -> Result<JointReactionImpulse, JointFrictionError> {
        let world = self.geometry.normal.scale(impulse);
        if self.angular() {
            Ok(JointReactionImpulse::couple(
                self.bodies.map(|b| states[b]),
                self.bodies.map(|b| properties[b]),
                world,
            )?)
        } else {
            Ok(JointReactionImpulse::Point(
                crate::apply_equal_and_opposite_impulse_at_body_points(
                    states[self.bodies[0]],
                    properties[self.bodies[0]],
                    self.geometry.arms_body_m[0],
                    states[self.bodies[1]],
                    properties[self.bodies[1]],
                    self.geometry.arms_body_m[1],
                    world,
                )?,
            ))
        }
    }

    fn motion(self, rate_after: f64, friction: f64, drive: f64) -> JointFrictionMotion {
        if self.angular() {
            JointFrictionMotion::Revolute {
                angular_length_m: self.length,
                rate_before_rad_s: self.rate_before,
                rate_after_rad_s: rate_after,
                friction_impulse_nms: friction * self.length,
            }
        } else {
            JointFrictionMotion::Prismatic {
                velocity_before_m_s: self.rate_before,
                velocity_after_m_s: rate_after,
                friction_impulse_ns: friction,
                drive_impulse_ns: drive,
            }
        }
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn run(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    sources: &[JointFrictionPort],
    controls: JointFrictionControls,
    single: Option<PrismaticFrictionInput>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<JointFrictionReport, JointFrictionError> {
    poll(&mut cancelled)?;
    validate_controls(controls)?;
    let graph = controls.graph;
    let n = sources.len();
    limit(
        "friction ports",
        n as u128,
        MAX_JOINT_FRICTION_PORTS as u128,
    )?;
    if states.len() != properties.len() {
        return Err(JointFrictionError::Invalid("body/property lengths"));
    }
    limit("bodies", states.len() as u128, graph.maximum_bodies as u128)?;
    let maximum_rows = graph.maximum_constraints.min(MAX_DISTANCE_ROWS);
    limit("joints", joints.len() as u128, maximum_rows as u128)?;
    let mut order = [0; MAX_JOINT_FRICTION_PORTS];
    for (i, source) in sources.iter().enumerate() {
        if !source.capacity().is_finite() || source.capacity() < 0.0 {
            return Err(JointFrictionError::Invalid(
                "finite nonnegative axial capacity",
            ));
        }
        order[i] = i;
    }
    order[..n].sort_unstable_by_key(|&i| sources[i].joint_index());
    if order[..n]
        .windows(2)
        .any(|p| sources[p[0]].joint_index() == sources[p[1]].joint_index())
    {
        return Err(JointFrictionError::Invalid("duplicate joint friction port"));
    }
    let mut bank = [Row::EMPTY; MAX_DISTANCE_ROWS];
    let mut count = 0;
    for &joint in joints {
        poll(&mut cancelled)?;
        let width = joint.row_count();
        limit("constraints", (count + width) as u128, maximum_rows as u128)?;
        joint.write_rows(&mut bank[count..count + width]);
        count += width;
    }
    let rows = Rows::Mixed(&bank[..count]);
    let bytes = planned_bytes_for::<JointReactionImpulse>(states.len(), count)
        + n as u128 * core::mem::size_of::<JointFrictionPortReport>() as u128;
    limit("live bytes", bytes, graph.maximum_live_bytes as u128)?;
    for (&state, &property) in states.iter().zip(properties) {
        poll(&mut cancelled)?;
        state.validate()?;
        property.validate()?;
    }
    let problem = Problem {
        states,
        properties,
        links: rows,
        integrator: RigidBodyIntegrator::new(Gravity::ZERO),
        flow: &[] as &[Wrench],
        duration_s: 0.0,
        controls: graph,
        targets: None,
    };
    let mut before: Vec<Geometry> = reserved(count)?;
    for (i, row) in rows.iter().enumerate() {
        poll(&mut cancelled)?;
        if row.bodies[0] == row.bodies[1] || row.bodies.iter().any(|&b| b >= states.len()) {
            return Err(JointFrictionError::Invalid(
                "distinct existing joint bodies",
            ));
        }
        let g = problem.geometry(row, row.bodies.map(|b| states[b]), None)?;
        if !g.rate.is_finite()
            || g.gap.abs() > graph.position_tolerance_m
            || g.rate.abs() > graph.velocity_tolerance_mps
        {
            return Err(DistanceGraphError::InitialConstraint { constraint: i }.into());
        }
        before.push(g);
    }
    if let Some(constraint) = rows.grouped_violation(
        &before,
        graph.position_tolerance_m,
        graph.velocity_tolerance_mps,
    ) {
        return Err(DistanceGraphError::InitialConstraint { constraint }.into());
    }
    let mut ports = [Port::EMPTY; MAX_JOINT_FRICTION_PORTS];
    let mut velocity = [0.0; MAX_JOINT_FRICTION_PORTS];
    let mut capacity = [0.0; MAX_JOINT_FRICTION_PORTS];
    for i in 0..n {
        poll(&mut cancelled)?;
        ports[i] = Port::prepare(sources[order[i]], joints, &problem)?;
        velocity[i] = ports[i].geometry.rate;
        capacity[i] = ports[i].capacity;
    }
    if matches!(single, Some(PrismaticFrictionInput::RestingDrive { .. })) && velocity[0] != 0.0 {
        return Err(JointFrictionError::NotResting {
            joint_index: ports[0].source.joint_index(),
            velocity_m_s: velocity[0],
        });
    }
    let mut coupling = [[0.0; MAX_DISTANCE_ROWS]; MAX_JOINT_FRICTION_PORTS];
    let mut response = [[0.0; MAX_DISTANCE_ROWS]; MAX_JOINT_FRICTION_PORTS];
    if n != 0 {
        let factor = Mobility::new(&problem, states, &before, &mut cancelled)?;
        for i in 0..n {
            for (k, (row, geometry)) in rows.iter().zip(&before).enumerate() {
                poll(&mut cancelled)?;
                let b = response_rate(row.bodies, row.is_angular(), *geometry, ports[i], states)
                    .to_f64();
                coupling[i][k] = b;
                response[i][k] = -b;
            }
            factor.solve(&mut response[i][..count])?;
        }
    }
    let mut matrix = [0.0; MAX_JOINT_FRICTION_PORTS * MAX_JOINT_FRICTION_PORTS];
    let mut published_matrix = [[0.0; MAX_JOINT_FRICTION_PORTS]; MAX_JOINT_FRICTION_PORTS];
    for i in 0..n {
        poll(&mut cancelled)?;
        for j in 0..n {
            let raw = response_rate(
                ports[i].bodies,
                ports[i].angular(),
                ports[i].geometry,
                ports[j],
                states,
            );
            let value = (0..count)
                .fold(raw, |sum, k| {
                    sum + Dd::from_f64(coupling[i][k]) * Dd::from_f64(response[j][k])
                })
                .to_f64();
            if !value.is_finite()
                || (i == j
                    && (!raw.to_f64().is_finite()
                        || raw.to_f64() <= 0.0
                        || value <= 0.0
                        || value / raw.to_f64() < graph.minimum_scaled_pivot))
            {
                return Err(JointFrictionError::Unresolved(
                    "locked or ill-conditioned friction direction",
                ));
            }
            matrix[i * n + j] = value;
            published_matrix[i][j] = value;
        }
    }
    let mut friction = [0.0; MAX_JOINT_FRICTION_PORTS];
    let mut drive = [0.0; MAX_JOINT_FRICTION_PORTS];
    let mut phase = JointFrictionPhase::FrozenCapacities;
    let projection = if let Some(input) = single {
        match input {
            PrismaticFrictionInput::Sliding { .. } => {
                let law = FrictionLaw::Coulomb {
                    static_mu: 1.0,
                    kinetic_mu: 1.0,
                };
                let projected = solve_isotropic_coulomb_impulse(
                    CoulombImpulseInput {
                        normal_impulse: capacity[0],
                        relative_tangent_velocity: [velocity[0], 0.0],
                        inverse_tangent_mass: matrix[0],
                    },
                    &law,
                )
                .map_err(JointFrictionError::Scalar)?;
                friction[0] = projected.tangent_impulse[0];
                phase = if projected.regime == FrictionRegime::Sticking {
                    JointFrictionPhase::Arrested
                } else {
                    JointFrictionPhase::Sliding
                };
            }
            PrismaticFrictionInput::RestingDrive {
                drive_impulse_ns,
                static_capacity_ns,
                kinetic_capacity_ns,
            } => {
                drive[0] = drive_impulse_ns;
                if drive_impulse_ns.abs() <= static_capacity_ns {
                    friction[0] = -drive_impulse_ns;
                    phase = JointFrictionPhase::Held;
                } else {
                    friction[0] = -drive_impulse_ns.signum() * kinetic_capacity_ns;
                    phase = JointFrictionPhase::Sliding;
                }
            }
        }
        None
    } else {
        let projected = solve_coupled_friction_impulses(
            &CoupledFrictionImpulseInput {
                velocity_before_m_s: &velocity[..n],
                mobility_per_kg: &matrix[..n * n],
                capacity_ns: &capacity[..n],
            },
            controls.solve,
            &mut cancelled,
        )?;
        friction[..n].copy_from_slice(&projected.impulse_ns[..n]);
        Some(projected)
    };
    let mut total = [0.0; MAX_JOINT_FRICTION_PORTS];
    for i in 0..n {
        total[i] = (Dd::from_f64(friction[i]) + Dd::from_f64(drive[i])).to_f64();
        let physical = friction[i] * ports[i].length;
        if [friction[i], drive[i], total[i]]
            .iter()
            .any(|x| !x.is_finite() || x.abs() > controls.solve.maximum_impulse_ns)
            || friction[i].abs() > capacity[i]
            || !physical.is_finite()
            || physical.abs() > ports[i].source.capacity()
            || (friction[i] != 0.0 && physical == 0.0)
        {
            return Err(JointFrictionError::Unresolved(
                "physical axial impulse bounds",
            ));
        }
    }
    let mut reactions = zeros(count)?;
    for (k, reaction) in reactions.iter_mut().enumerate() {
        *reaction = (0..n)
            .fold(Dd::ZERO, |sum, i| {
                sum + Dd::from_f64(response[i][k]) * Dd::from_f64(total[i])
            })
            .to_f64();
    }
    if !rows.admits_impulses(&reactions, &before, graph.step.maximum_impulse_ns) {
        return Err(JointFrictionError::Unresolved("bilateral reaction bound"));
    }
    let mut after_states = reserved(states.len())?;
    after_states.extend_from_slice(states);
    let mut reports: Vec<JointFrictionPortReport> = reserved(n)?;
    for i in 0..n {
        poll(&mut cancelled)?;
        let event = ports[i].action(&after_states, properties, total[i])?;
        let after = event.states_after();
        for side in 0..2 {
            after_states[ports[i].bodies[side]] = after[side];
        }
        reports.push(JointFrictionPortReport {
            joint_index: ports[i].source.joint_index(),
            axis_world: ports[i].axis,
            motion: ports[i].motion(ports[i].rate_before, friction[i], drive[i]),
            native_axial_impulse: event,
            scaled_velocity_residual_m_s: 0.0,
            drive_work_j: 0.0,
            friction_work_j: 0.0,
        });
    }
    let mut native_constraints = reserved(count)?;
    problem.apply(
        &mut after_states,
        &reactions,
        &before,
        &mut native_constraints,
        &mut cancelled,
    )?;
    let mut drive_work = Dd::ZERO;
    let mut friction_work = Dd::ZERO;
    let mut gap = Dd::ZERO;
    for i in 0..n {
        poll(&mut cancelled)?;
        let (physical_rate, actual) = ports[i].velocity(&after_states, properties, &problem)?;
        let target = (0..n).fold(Dd::from_f64(velocity[i]), |sum, j| {
            sum + Dd::from_f64(matrix[i * n + j]) * Dd::from_f64(total[j])
        });
        let residual = (Dd::from_f64(actual) - target).to_f64();
        let stationarity = if capacity[i] == 0.0 {
            0.0
        } else if friction[i] == -capacity[i] {
            (-actual).max(0.0)
        } else if friction[i] == capacity[i] {
            actual.max(0.0)
        } else {
            actual.abs()
        };
        if !residual.is_finite()
            || !actual.is_finite()
            || residual.abs() > controls.solve.velocity_tolerance_m_s
            || stationarity > controls.solve.velocity_tolerance_m_s
        {
            return Err(JointFrictionError::Unresolved(
                "actual native axial endpoint",
            ));
        }
        gap = gap
            + Dd::from_f64(capacity[i]) * Dd::from_f64(actual.abs())
            + Dd::from_f64(friction[i]) * Dd::from_f64(actual);
        let midpoint = (Dd::from_f64(velocity[i]) + Dd::from_f64(actual)) * Dd::from_f64(0.5);
        let dw = Dd::from_f64(drive[i]) * midpoint;
        let fw = Dd::from_f64(friction[i]) * midpoint;
        drive_work = drive_work + dw;
        friction_work = friction_work + fw;
        reports[i].motion = ports[i].motion(physical_rate, friction[i], drive[i]);
        reports[i].scaled_velocity_residual_m_s = residual;
        reports[i].drive_work_j = dw.to_f64();
        reports[i].friction_work_j = fw.to_f64();
    }
    let mut constraint_work = Dd::ZERO;
    let mut after = reserved(count)?;
    for (i, row) in rows.iter().enumerate() {
        poll(&mut cancelled)?;
        let g = problem.geometry(row, row.bodies.map(|b| after_states[b]), None)?;
        constraint_work = constraint_work
            + Dd::from_f64(reactions[i])
                * (Dd::from_f64(before[i].rate) + Dd::from_f64(g.rate))
                * Dd::from_f64(0.5);
        after.push(g);
    }
    if rows
        .grouped_violation(
            &after,
            graph.position_tolerance_m,
            graph.velocity_tolerance_mps,
        )
        .is_some()
        || after
            .iter()
            .any(|g| !g.rate.is_finite() || g.rate.abs() > graph.velocity_tolerance_mps)
    {
        return Err(JointFrictionError::Unresolved(
            "actual native bilateral endpoint",
        ));
    }
    let (energy, linear, angular) = account(states, &after_states, properties, &mut cancelled)?;
    let work_residual = (energy - drive_work - friction_work - constraint_work).to_f64();
    if [
        energy.to_f64(),
        drive_work.to_f64(),
        friction_work.to_f64(),
        constraint_work.to_f64(),
        gap.to_f64(),
        work_residual,
    ]
    .iter()
    .any(|x| !x.is_finite())
        || !linear.is_finite()
        || !angular.is_finite()
        || friction_work.to_f64() > 0.0
        || constraint_work.to_f64().abs() > controls.work_tolerance_j
        || work_residual.abs() > controls.work_tolerance_j
        || gap.to_f64().abs() > controls.solve.energy_tolerance_j
    {
        return Err(JointFrictionError::Unresolved(
            "native work, support and passivity",
        ));
    }
    if norm(linear) > controls.linear_impulse_tolerance_ns
        || norm(angular) > controls.angular_impulse_tolerance_nms
    {
        return Err(JointFrictionError::Unresolved("native total momentum"));
    }
    poll(&mut cancelled)?;
    Ok(JointFrictionReport {
        states_after: after_states,
        phase,
        ports: reports,
        scaled_mobility_per_kg: published_matrix,
        projection,
        constraint_impulses_ns: reactions,
        native_constraint_impulses: native_constraints,
        drive_work_j: drive_work.to_f64(),
        friction_work_j: friction_work.to_f64(),
        constraint_work_j: constraint_work.to_f64(),
        kinetic_energy_change_j: energy.to_f64(),
        work_energy_residual_j: work_residual,
        native_duality_gap_j: gap.to_f64(),
        linear_momentum_change_ns: linear,
        angular_momentum_change_nms: angular,
        planned_live_bytes: bytes as usize,
    })
}

fn validate_controls(controls: JointFrictionControls) -> Result<(), JointFrictionError> {
    controls.solve.validate()?;
    let graph = controls.graph;
    if [
        graph.position_tolerance_m,
        graph.velocity_tolerance_mps,
        graph.step.maximum_impulse_ns,
        controls.work_tolerance_j,
        controls.linear_impulse_tolerance_ns,
        controls.angular_impulse_tolerance_nms,
    ]
    .iter()
    .any(|x| !x.is_finite() || *x <= 0.0)
        || !graph.minimum_scaled_pivot.is_finite()
        || !(4096.0 * f64::EPSILON..1.0).contains(&graph.minimum_scaled_pivot)
    {
        return Err(JointFrictionError::Invalid(
            "finite positive native controls",
        ));
    }
    Ok(())
}

// Conversion is conservative by at most one representable step. In particular,
// changing an authored length cannot round an angular impulse above its original
// physical capacity. Unrepresentable nonzero ports refuse instead of disappearing.
fn scaled_capacity(capacity: f64, length: f64) -> Result<f64, JointFrictionError> {
    let mut scaled = capacity / length;
    if !scaled.is_finite() || (capacity != 0.0 && scaled == 0.0) {
        return Err(JointFrictionError::Invalid("representable scaled capacity"));
    }
    if (Dd::from_f64(scaled) * Dd::from_f64(length) - Dd::from_f64(capacity)).to_f64() > 0.0 {
        scaled = scaled.next_down();
    }
    if scaled < 0.0 || (capacity != 0.0 && scaled == 0.0) {
        return Err(JointFrictionError::Invalid("representable scaled capacity"));
    }
    Ok(scaled)
}

fn response_rate(
    bodies: [usize; 2],
    angular: bool,
    geometry: Geometry,
    source: Port,
    states: &[RigidBodyState],
) -> Dd {
    let mut value = Dd::ZERO;
    for side in 0..2 {
        for other in 0..2 {
            if bodies[side] == source.bodies[other] {
                let response = if angular {
                    dot(source.spin_response[other], geometry.normal)
                } else {
                    let arm = states[bodies[side]]
                        .pose()
                        .orientation()
                        .rotate_body_to_world(geometry.arms_body_m[side]);
                    dot(source.linear_response[other], geometry.normal)
                        + dot(source.spin_response[other].cross(arm), geometry.normal)
                };
                value = value + response * Dd::from_f64(if side == 0 { -1.0 } else { 1.0 });
            }
        }
    }
    value
}

fn difference(a: Vec3, b: Vec3) -> [Dd; 3] {
    let a = [a.x, a.y, a.z];
    let b = [b.x, b.y, b.z];
    core::array::from_fn(|i| Dd::from_f64(a[i]) - Dd::from_f64(b[i]))
}
fn dot(a: Vec3, b: Vec3) -> Dd {
    project([a.x, a.y, a.z].map(Dd::from_f64), b)
}
fn project(a: [Dd; 3], b: Vec3) -> Dd {
    a.into_iter()
        .zip([b.x, b.y, b.z])
        .fold(Dd::ZERO, |sum, (a, b)| sum + a * Dd::from_f64(b))
}
fn norm(a: Vec3) -> f64 {
    a.x.hypot(a.y).hypot(a.z)
}
