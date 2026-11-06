//! Finite-time Coulomb motion on the original constrained hinge/slider ports.
//!
//! The geometry, inertia and bilateral Jacobian stay fixed during this phase.
//! Constant physical drives and force capacities generate a piecewise-linear
//! velocity trajectory with isolated arrest events and set-valued holding.
//! Its integrated actions update the original body momenta once; elapsed time,
//! travel and physical path work are retained separately from the frozen poses.

/// Adaptive native pose evolution with changing authored axial loads.
#[path = "relaxation/trajectory.rs"]
pub mod trajectory;

use super::{
    ArticulatedJoint, Dd, DistanceGraphControls, DistanceGraphError, DistanceStepError,
    JointFrictionControls, JointFrictionError, JointFrictionPort, JointReactionImpulse,
    MAX_JOINT_FRICTION_PORTS, MassProperties, RigidBodyState, Vec3, apply_native, limit, poll,
    prepare, reserved, scaled_capacity,
};
use fs_contact::tangential::impact::multiple::relaxation::stiction::{
    CoupledStictionFlowControls, CoupledStictionFlowDiagnostics, CoupledStictionFlowError,
    CoupledStictionFlowInput, MAX_STICTION_MODES_PER_PHASE, advance_coupled_stiction,
};
use fs_contact::tangential::impact::multiple::relaxation::{
    CoupledFrictionFlowControls, CoupledFrictionFlowError, CoupledFrictionFlowInput,
    CoupledFrictionFlowReport, advance_coupled_friction,
};

/// Complete six-port static/release mode inventory per isolated phase.
/// The whole native flow admits at most this value times its segment limit.
pub const MAX_JOINT_STICTION_MODES_PER_PHASE: usize = MAX_STICTION_MODES_PER_PHASE;

/// A held physical drive and source-declared friction capacities at one joint.
///
/// The original variants retain a common static/kinetic bound. The explicit
/// stiction variants use a larger static bound and lower moving bound, with the
/// shared contact owner's maximal-sticking continuation and ambiguity refusal.
/// These are declared physical inputs, not inferred bearing loads or materials.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum JointFrictionFlowPort {
    /// Positive torque acts on the second body along the original hinge axis.
    Revolute {
        /// Index of an original revolute joint.
        joint_index: usize,
        /// Nonnegative Coulomb torque capacity [N m].
        capacity_nm: f64,
        /// Signed held drive torque [N m].
        drive_nm: f64,
    },
    /// Positive force acts on the carriage along the carrier's joint-frame Z.
    Prismatic {
        /// Index of an original prismatic joint.
        joint_index: usize,
        /// Nonnegative Coulomb force capacity [N].
        capacity_n: f64,
        /// Signed held drive force [N].
        drive_n: f64,
    },
    /// Axial hinge friction with an independently declared breakaway capacity.
    RevoluteStiction {
        /// Index of an original revolute joint.
        joint_index: usize,
        /// Nonnegative static torque capacity, at least the kinetic value [N m].
        static_capacity_nm: f64,
        /// Nonnegative sliding torque capacity [N m].
        kinetic_capacity_nm: f64,
        /// Signed held torque on the second body along the hinge axis [N m].
        drive_nm: f64,
    },
    /// Axial carriage friction with an independently declared breakaway force.
    PrismaticStiction {
        /// Index of an original prismatic joint.
        joint_index: usize,
        /// Nonnegative static force capacity, at least the kinetic value [N].
        static_capacity_n: f64,
        /// Nonnegative sliding force capacity [N].
        kinetic_capacity_n: f64,
        /// Signed held force on the carriage along the carrier's joint Z [N].
        drive_n: f64,
    },
}
impl JointFrictionFlowPort {
    /// Original joint index; caller ordering has no physical effect.
    #[must_use]
    pub const fn joint_index(self) -> usize {
        match self {
            Self::Revolute { joint_index, .. }
            | Self::Prismatic { joint_index, .. }
            | Self::RevoluteStiction { joint_index, .. }
            | Self::PrismaticStiction { joint_index, .. } => joint_index,
        }
    }
    fn capacity(self) -> f64 {
        match self {
            Self::Revolute { capacity_nm, .. } => capacity_nm,
            Self::Prismatic { capacity_n, .. } => capacity_n,
            Self::RevoluteStiction {
                static_capacity_nm, ..
            } => static_capacity_nm,
            Self::PrismaticStiction {
                static_capacity_n, ..
            } => static_capacity_n,
        }
    }
    fn kinetic_capacity(self) -> f64 {
        match self {
            Self::Revolute { capacity_nm, .. } => capacity_nm,
            Self::Prismatic { capacity_n, .. } => capacity_n,
            Self::RevoluteStiction {
                kinetic_capacity_nm,
                ..
            } => kinetic_capacity_nm,
            Self::PrismaticStiction {
                kinetic_capacity_n, ..
            } => kinetic_capacity_n,
        }
    }
    fn drive(self) -> f64 {
        match self {
            Self::Revolute { drive_nm, .. } | Self::RevoluteStiction { drive_nm, .. } => drive_nm,
            Self::Prismatic { drive_n, .. } | Self::PrismaticStiction { drive_n, .. } => drive_n,
        }
    }
    fn impulse_bound(self, duration_s: f64) -> Result<JointFrictionPort, JointFrictionFlowError> {
        let capacity = self.capacity();
        let kinetic = self.kinetic_capacity();
        let bound = capacity * duration_s;
        if !capacity.is_finite()
            || capacity < 0.0
            || !kinetic.is_finite()
            || kinetic < 0.0
            || kinetic > capacity
            || !self.drive().is_finite()
            || !bound.is_finite()
            || (capacity != 0.0 && bound == 0.0)
        {
            return Err(JointFrictionFlowError::Invalid(
                "representable held joint efforts",
            ));
        }
        Ok(match self {
            Self::Revolute { joint_index, .. } | Self::RevoluteStiction { joint_index, .. } => {
                JointFrictionPort::Revolute {
                    joint_index,
                    capacity_nms: bound,
                }
            }
            Self::Prismatic { joint_index, .. } | Self::PrismaticStiction { joint_index, .. } => {
                JointFrictionPort::Prismatic {
                    joint_index,
                    capacity_ns: bound,
                }
            }
        })
    }
}

/// Physical flow controls and the original native graph/momentum gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JointFrictionFlowControls {
    /// Original geometry, constraint, reaction, factor and memory bounds.
    /// Free-flight timestep/rotation iteration fields are unused by this phase.
    pub graph: DistanceGraphControls,
    /// Continuous force-law controls in work-conjugate linear coordinates.
    /// A hinge uses velocity `L*spin` and effort `torque/L`.
    pub flow: CoupledFrictionFlowControls,
    /// Native path-work closure and bounded bilateral residual work [J].
    pub work_tolerance_j: f64,
    /// Total internal linear-momentum increment allowance [N s].
    pub linear_impulse_tolerance_ns: f64,
    /// Total spin-plus-orbital angular-momentum increment allowance [N m s].
    pub angular_impulse_tolerance_nms: f64,
}
impl JointFrictionFlowControls {
    fn native(self) -> JointFrictionControls {
        JointFrictionControls {
            graph: self.graph,
            solve: self.flow.projection,
            work_tolerance_j: self.work_tolerance_j,
            linear_impulse_tolerance_ns: self.linear_impulse_tolerance_ns,
            angular_impulse_tolerance_nms: self.angular_impulse_tolerance_nms,
        }
    }
}

/// Original physical rates, integrated actions and path travel at one joint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum JointFrictionFlowMotion {
    /// Axial spin of the second body relative to the first body.
    Revolute {
        /// Original scale used in the constrained solve [m].
        angular_length_m: f64,
        /// Actual initial relative angular velocity [rad/s].
        rate_before_rad_s: f64,
        /// Actual final relative angular velocity, retaining roundoff [rad/s].
        rate_after_rad_s: f64,
        /// Integrated signed relative rotation along the fixed-pose phase [rad].
        rotation_rad: f64,
        /// Integrated absolute relative angular travel [rad].
        angular_distance_rad: f64,
        /// Integrated signed friction torque on body two [N m s].
        friction_impulse_nms: f64,
        /// Integrated signed drive torque on body two [N m s].
        drive_impulse_nms: f64,
    },
    /// Carriage motion along the actual moving-rail axial Jacobian.
    Prismatic {
        /// Actual initial axial speed [m/s].
        velocity_before_m_s: f64,
        /// Actual final axial speed, retaining roundoff [m/s].
        velocity_after_m_s: f64,
        /// Integrated signed axial travel during the frozen-geometry phase [m].
        displacement_m: f64,
        /// Integrated absolute axial travel [m].
        distance_m: f64,
        /// Integrated signed friction force on the carriage [N s].
        friction_impulse_ns: f64,
        /// Integrated signed drive force on the carriage [N s].
        drive_impulse_ns: f64,
    },
}

/// One original joint's actual native update and continuous path accounting.
#[derive(Clone, Debug, PartialEq)]
pub struct JointFrictionFlowPortReport {
    /// Original joint index, in canonical increasing order.
    pub joint_index: usize,
    /// Actual unit axis used for the original point-force or couple port.
    pub axis_world: Vec3,
    /// Physical native endpoint, impulses and time-integrated coordinate travel.
    pub motion: JointFrictionFlowMotion,
    /// Minimum nonnegative speed over the integrated phase, also enclosing the
    /// measured native start/end speeds [m/s]. A hinge uses `L*|omega|`, where
    /// L is its reported angular length, not an inferred physical rubbing radius.
    pub minimum_scaled_speed_m_s: f64,
    /// Maximum nonnegative speed over that same complete phase [m/s]. These
    /// computed extrema are not a validated trajectory tube.
    pub maximum_scaled_speed_m_s: f64,
    /// Native momentum application of the summed drive and friction action.
    /// Its stage-local work is distinct from the path work below.
    pub native_axial_impulse: JointReactionImpulse,
    /// Actual native endpoint minus the original mobility prediction [m/s].
    pub scaled_velocity_residual_m_s: f64,
    /// Actual native endpoint minus the piecewise physical flow endpoint [m/s].
    pub flow_velocity_residual_m_s: f64,
    /// Actual integrated capacity excess, retained as a numerical residual [N s].
    /// Actions are never clipped; the shared impulse tolerance bounds this.
    pub capacity_residual_ns: f64,
    /// Scaled actual drive impulse minus original held drive times duration [N s].
    pub drive_impulse_residual_ns: f64,
    /// Signed integral of held drive times actual coordinate rate [J].
    pub drive_work_j: f64,
    /// Nonpositive integral of friction times coordinate rate [J].
    /// This is a mechanical loss, without an inferred heat partition.
    pub friction_work_j: f64,
}

/// Atomic native candidate plus the complete finite-time friction trajectory.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct JointFrictionFlowReport {
    /// Original body order, poses unchanged, momenta updated by integrated loads.
    pub states_after: Vec<RigidBodyState>,
    /// Physical elapsed duration [s].
    pub duration_s: f64,
    /// Canonically ordered original joint reports.
    pub ports: Vec<JointFrictionFlowPortReport>,
    /// Complete constrained mobility in canonical port order [1/kg].
    pub scaled_mobility_per_kg: [[f64; MAX_JOINT_FRICTION_PORTS]; MAX_JOINT_FRICTION_PORTS],
    /// Shared flow in the same canonical radius/length-scaled coordinates.
    pub flow: CoupledFrictionFlowReport,
    /// Distinct static/kinetic continuation diagnostics in canonical port order.
    /// Absent when all bounds are equal and the original flow path is used.
    pub stiction: Option<CoupledStictionFlowDiagnostics>,
    /// Summed original bilateral impulses in authored row order [N s].
    pub constraint_impulses_ns: Vec<f64>,
    /// Original native point/couple receipts for the bilateral reactions.
    pub native_constraint_impulses: Vec<JointReactionImpulse>,
    /// Total signed physical drive work along the complete trajectory [J].
    pub drive_work_j: f64,
    /// Total physical mechanical loss from the continuous friction law [J].
    pub dissipated_j: f64,
    /// Bilateral residual work, separately bounded and never thermalized [J].
    pub constraint_work_j: f64,
    /// Measured from actual original-body momentum differences [J].
    pub kinetic_energy_change_j: f64,
    /// Summed integrated actions paired with original/final midpoint rates [J].
    /// Only their sum, not its drive/friction split, equals physical path work.
    pub native_impulse_work_j: f64,
    /// Native kinetic change minus drive work plus loss minus constraint work [J].
    pub work_energy_residual_j: f64,
    /// Actual total internal linear-momentum increment [N s].
    pub linear_momentum_change_ns: Vec3,
    /// Actual total internal angular-momentum increment [N m s].
    pub angular_momentum_change_nms: Vec3,
    /// Admitted logical heap payload, excluding bounded stack banks.
    pub planned_live_bytes: usize,
}

/// Failed geometry, flow, resource or native audits preserve every borrowed input.
#[derive(Clone, Debug, PartialEq)]
pub enum JointFrictionFlowError {
    /// Original graph mechanics or native-state admission failure.
    Native(JointFrictionError),
    /// Original continuous Coulomb owner refusal.
    Flow(CoupledFrictionFlowError),
    /// Larger-static continuation refused, including ambiguous coupled modes.
    Stiction(CoupledStictionFlowError),
    /// Invalid or unrepresentable physical input.
    Invalid(&'static str),
    /// Actual native response failed a named physical residual.
    Unresolved(&'static str),
}
impl From<JointFrictionError> for JointFrictionFlowError {
    fn from(error: JointFrictionError) -> Self {
        Self::Native(error)
    }
}
impl From<DistanceGraphError> for JointFrictionFlowError {
    fn from(error: DistanceGraphError) -> Self {
        Self::Native(error.into())
    }
}
impl From<DistanceStepError> for JointFrictionFlowError {
    fn from(error: DistanceStepError) -> Self {
        Self::Native(error.into())
    }
}
impl From<CoupledFrictionFlowError> for JointFrictionFlowError {
    fn from(error: CoupledFrictionFlowError) -> Self {
        Self::Flow(error)
    }
}
impl From<CoupledStictionFlowError> for JointFrictionFlowError {
    fn from(error: CoupledStictionFlowError) -> Self {
        match error {
            // Existing trajectory owners may refine an interval after a flow
            // impulse/segment budget refusal. Preserve that shared error path.
            CoupledStictionFlowError::Flow(flow) => Self::Flow(flow),
            mode => Self::Stiction(mode),
        }
    }
}
impl core::fmt::Display for JointFrictionFlowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native joint friction flow: {self:?}")
    }
}
impl std::error::Error for JointFrictionFlowError {}

/// Advance coupled axial Coulomb motion under held physical efforts.
///
/// The shared force law resolves real arrest times and zero-speed holding.
/// Original equal-bound ports retain the original algorithm. Explicit stiction
/// ports use their larger static and lower kinetic bounds, preserving maximal
/// admissible sticking sets; multiple undominated continuations refuse instead
/// of selecting a mode by joint order. Each phase has at most 3^6 mode candidates,
/// and the segment limit bounds the complete continuation work. The original
/// full constrained mobility is shared across every selected joint.
///
/// No numerical speed capture, inferred bearing pressure, pose advance or thermal
/// partition is supplied. Rest is an exact input velocity or a resolved zero
/// event within this phase. Native endpoint roundoff remains visible; it is not
/// silently replaced by a remembered stick state on the next call. A persistent
/// geometrically locked trajectory needs its surrounding constraint/history
/// owner. Successive calls can change declared loads/drives, and moved geometry
/// must be re-admitted by the original native graph.
///
/// Unlike an impact projection, continuous work integrates the whole velocity
/// path. This function therefore audits actual native endpoint and path energy,
/// without requiring the integrated friction impulse to maximize dissipation at
/// the final velocity. Native states are returned only after every physical gate
/// and a final cancellation poll succeeds.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_joint_friction(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    sources: &[JointFrictionFlowPort],
    duration_s: f64,
    controls: JointFrictionFlowControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<JointFrictionFlowReport, JointFrictionFlowError> {
    poll(&mut cancelled)?;
    if !duration_s.is_finite() || duration_s <= 0.0 {
        return Err(JointFrictionFlowError::Invalid(
            "positive finite flow duration",
        ));
    }
    let n = sources.len();
    limit(
        "friction ports",
        n as u128,
        MAX_JOINT_FRICTION_PORTS as u128,
    )?;
    let mut bounds = [JointFrictionPort::Prismatic {
        joint_index: 0,
        capacity_ns: 0.0,
    }; MAX_JOINT_FRICTION_PORTS];
    for i in 0..n {
        bounds[i] = sources[i].impulse_bound(duration_s)?;
    }
    let has_stiction = sources.iter().any(|s| s.capacity() > s.kinetic_capacity());
    let native = controls.native();
    let prepared = prepare(
        states,
        properties,
        joints,
        &bounds[..n],
        native,
        core::mem::size_of::<JointFrictionFlowPortReport>(),
        &mut cancelled,
    )?;
    let mut capacities = [0.0; MAX_JOINT_FRICTION_PORTS];
    let mut kinetic = [0.0; MAX_JOINT_FRICTION_PORTS];
    let mut drives = [0.0; MAX_JOINT_FRICTION_PORTS];
    for i in 0..n {
        let source = sources[prepared.order[i]];
        let length = prepared.ports[i].length;
        if has_stiction {
            // The static threshold and drive must use the same conversion.
            // Conservatively rounding only Cs down would cause an artificial
            // finite breakaway when the authored drive is exactly Cs. Native
            // integrated impulse bounds remain independently audited below.
            capacities[i] = source.capacity() / length;
            kinetic[i] = source.kinetic_capacity() / length;
            if !capacities[i].is_finite()
                || !kinetic[i].is_finite()
                || (source.capacity() != 0.0 && capacities[i] == 0.0)
                || (source.kinetic_capacity() != 0.0 && kinetic[i] == 0.0)
            {
                return Err(JointFrictionFlowError::Invalid(
                    "representable scaled static/kinetic capacities",
                ));
            }
        } else {
            capacities[i] = scaled_capacity(source.capacity(), length)?;
            kinetic[i] = capacities[i];
        }
        if source.capacity() > source.kinetic_capacity() && capacities[i] == kinetic[i] {
            return Err(JointFrictionFlowError::Invalid(
                "representable static/kinetic capacity gap",
            ));
        }
        drives[i] = source.drive() / length;
        if !drives[i].is_finite() || (source.drive() != 0.0 && drives[i] == 0.0) {
            return Err(JointFrictionFlowError::Invalid(
                "representable scaled drive",
            ));
        }
    }
    let (flow, stiction) = if has_stiction {
        let maximum_mode_evaluations = MAX_JOINT_STICTION_MODES_PER_PHASE
            .checked_mul(controls.flow.maximum_segments)
            .ok_or(JointFrictionFlowError::Invalid("stiction mode work bound"))?;
        let resolved = advance_coupled_stiction(
            &CoupledStictionFlowInput {
                velocity_before_m_s: &prepared.velocity[..n],
                mobility_per_kg: &prepared.matrix[..n * n],
                static_capacity_n: &capacities[..n],
                kinetic_capacity_n: &kinetic[..n],
                drive_n: &drives[..n],
                duration_s,
            },
            CoupledStictionFlowControls {
                flow: controls.flow,
                maximum_mode_evaluations,
            },
            &mut cancelled,
        )?;
        (resolved.flow, Some(resolved.diagnostics))
    } else {
        let flow = advance_coupled_friction(
            &CoupledFrictionFlowInput {
                velocity_before_m_s: &prepared.velocity[..n],
                mobility_per_kg: &prepared.matrix[..n * n],
                capacity_n: &capacities[..n],
                drive_n: &drives[..n],
                duration_s,
            },
            controls.flow,
            &mut cancelled,
        )?;
        (flow, None)
    };
    let mut friction = [0.0; MAX_JOINT_FRICTION_PORTS];
    let mut drive = [0.0; MAX_JOINT_FRICTION_PORTS];
    friction[..n].copy_from_slice(&flow.friction_impulse_ns[..n]);
    drive[..n].copy_from_slice(&flow.drive_impulse_ns[..n]);
    let applied = apply_native(
        &prepared,
        states,
        properties,
        &friction,
        &drive,
        native,
        controls.flow.projection.impulse_tolerance_ns,
        &mut cancelled,
    )?;
    let mut ports = reserved(n)?;
    let mut drive_work = Dd::ZERO;
    let mut friction_work = Dd::ZERO;
    for i in 0..n {
        poll(&mut cancelled)?;
        let port = prepared.ports[i];
        let source = sources[prepared.order[i]];
        let endpoint = applied.velocity_after[i];
        let flow_residual = endpoint - flow.flow_velocity_after_m_s[i];
        let impulse_residual = endpoint - flow.velocity_after_m_s[i];
        let drive_residual = (Dd::from_f64(drive[i]) * Dd::from_f64(port.length)
            - Dd::from_f64(source.drive()) * Dd::from_f64(duration_s))
        .to_f64()
            / port.length;
        if [flow_residual, impulse_residual, drive_residual]
            .iter()
            .any(|x| !x.is_finite())
            || flow_residual.abs() > controls.flow.projection.velocity_tolerance_m_s
            || impulse_residual.abs() > controls.flow.projection.velocity_tolerance_m_s
            || drive_residual.abs() > controls.flow.projection.impulse_tolerance_ns
        {
            return Err(JointFrictionFlowError::Unresolved(
                "original native flow endpoint",
            ));
        }
        let travel = flow.displacement_m[i] / port.length;
        let distance = flow.distance_m[i] / port.length;
        let physical_friction = friction[i] * port.length;
        let physical_drive = drive[i] * port.length;
        if [travel, distance, physical_friction, physical_drive]
            .iter()
            .any(|x| !x.is_finite())
            || (flow.displacement_m[i] != 0.0 && travel == 0.0)
            || (flow.distance_m[i] != 0.0 && distance == 0.0)
            || (drive[i] != 0.0 && physical_drive == 0.0)
            || !flow.friction_work_j[i].is_finite()
            || flow.friction_work_j[i] > 0.0
            || !flow.drive_work_j[i].is_finite()
        {
            return Err(JointFrictionFlowError::Unresolved(
                "representable native path quantities",
            ));
        }
        let motion = if port.angular() {
            JointFrictionFlowMotion::Revolute {
                angular_length_m: port.length,
                rate_before_rad_s: port.rate_before,
                rate_after_rad_s: applied.physical_rate_after[i],
                rotation_rad: travel,
                angular_distance_rad: distance,
                friction_impulse_nms: physical_friction,
                drive_impulse_nms: physical_drive,
            }
        } else {
            JointFrictionFlowMotion::Prismatic {
                velocity_before_m_s: port.rate_before,
                velocity_after_m_s: applied.physical_rate_after[i],
                displacement_m: travel,
                distance_m: distance,
                friction_impulse_ns: physical_friction,
                drive_impulse_ns: physical_drive,
            }
        };
        drive_work = drive_work + Dd::from_f64(flow.drive_work_j[i]);
        friction_work = friction_work + Dd::from_f64(flow.friction_work_j[i]);
        ports.push(JointFrictionFlowPortReport {
            joint_index: source.joint_index(),
            axis_world: port.axis,
            motion,
            minimum_scaled_speed_m_s: flow.minimum_speed_m_s[i]
                .min(prepared.velocity[i].abs())
                .min(endpoint.abs()),
            maximum_scaled_speed_m_s: flow.maximum_speed_m_s[i]
                .max(prepared.velocity[i].abs())
                .max(endpoint.abs()),
            native_axial_impulse: applied.axial_events[i].ok_or(
                JointFrictionFlowError::Unresolved("missing native axial update"),
            )?,
            scaled_velocity_residual_m_s: applied.velocity_residual[i],
            flow_velocity_residual_m_s: flow_residual,
            capacity_residual_ns: applied.capacity_residual[i],
            drive_impulse_residual_ns: drive_residual,
            drive_work_j: flow.drive_work_j[i],
            friction_work_j: flow.friction_work_j[i],
        });
    }
    let residual = (applied.energy - drive_work - friction_work - applied.constraint_work).to_f64();
    let path_residual = (applied.impulse_work - drive_work - friction_work).to_f64();
    let energy_residual =
        (applied.energy - applied.constraint_work - Dd::from_f64(flow.kinetic_energy_change_j))
            .to_f64();
    if [
        residual,
        path_residual,
        energy_residual,
        drive_work.to_f64(),
        friction_work.to_f64(),
    ]
    .iter()
    .any(|x| !x.is_finite())
        || friction_work.to_f64() > 0.0
        || residual.abs() > controls.work_tolerance_j
        || path_residual.abs() > controls.work_tolerance_j
        || energy_residual.abs() > controls.work_tolerance_j
        || (drive_work.to_f64() - flow.total_drive_work_j).abs() > controls.work_tolerance_j
        || (friction_work.to_f64() + flow.dissipated_j).abs() > controls.work_tolerance_j
    {
        return Err(JointFrictionFlowError::Unresolved(
            "native path work and dissipation",
        ));
    }
    poll(&mut cancelled)?;
    Ok(JointFrictionFlowReport {
        states_after: applied.states_after,
        duration_s,
        ports,
        scaled_mobility_per_kg: prepared.published_matrix,
        flow,
        stiction,
        constraint_impulses_ns: applied.reactions,
        native_constraint_impulses: applied.native_constraints,
        drive_work_j: drive_work.to_f64(),
        dissipated_j: -friction_work.to_f64(),
        constraint_work_j: applied.constraint_work.to_f64(),
        kinetic_energy_change_j: applied.energy.to_f64(),
        native_impulse_work_j: applied.impulse_work.to_f64(),
        work_energy_residual_j: residual,
        linear_momentum_change_ns: applied.linear,
        angular_momentum_change_nms: applied.angular,
        planned_live_bytes: prepared.bytes,
    })
}
