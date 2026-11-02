//! Native central Hertz contact with simultaneous retained Mindlin shear.
//!
//! The normal Hertz force and tangential history are evaluated along the SAME
//! affine indentation/shear segment. Actual masses, full inertias and contact
//! arms determine midpoint recoil. The source-owned conservative normal mean
//! drives both normal momentum and work; the source-owned Mindlin mean drives
//! shear momentum, storage and heat. No caller effective mass is accepted.
//! Held world forces and body torques can maintain compression and drive the
//! same retained shear path. Their actual native work remains an external port.
//!
//! This is a frozen-geometry momentum phase, like the adjacent rigid-friction
//! owners. The caller owns initial indentation, circular surface geometry,
//! material-point association, pose integration and trajectory refinement.
//! Central force lines decouple normal and shear mobility. Noncentral coupling,
//! transverse shear, torsional slip and intervals reaching separation refuse;
//! no unrepresented coordinate is projected away. This is not an impact law,
//! contact locator, rolling-memory transport or a finite-rotation integrator.

use crate::{
    DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3, Wrench,
    apply_equal_and_opposite_impulse_at_body_points,
};
use fs_contact::interface_binding::mindlin::loading::{
    SourceHertzMindlinContact, SourceHertzMindlinState,
    path::{SourceMindlinPathError, SourceMindlinPathInput, SourceMindlinPathStep},
};
use fs_contact::normal_patch::{
    NormalPatchError,
    source::{PointContactGeometry, SourcePointContact, path::SourceElasticPathAverage},
};
use fs_tribo::{TangentialLaw, TangentialSolveConfig, partial_slip::TangentFrame};

#[path = "mindlin/driven.rs"]
mod driven;
#[path = "mindlin/multiple.rs"]
pub mod multiple;
#[path = "mindlin/normal.rs"]
mod normal;

use driven::HeldDrive;
use normal::solve_normal;

/// Original A/B source ordering and physical attachment, held over one phase.
#[derive(Clone, Copy, Debug)]
pub struct RigidMindlinContact<'a> {
    pub law: &'a SourceHertzMindlinContact,
    /// COM-relative points in each body's declared frame. The world points
    /// must coincide, and both arms must lie on the declared normal line.
    pub arms_body: [Vec3; 2],
    /// Normal points from A toward B; the longitudinal tangent is scalar shear.
    /// Positive A-minus-B motion increases indentation and interface shear.
    pub frame: TangentFrame,
}

/// Independent dimensional gates and the original constitutive solve budget.
#[derive(Clone, Copy, Debug)]
pub struct RigidMindlinControls {
    pub solve: TangentialSolveConfig,
    /// Bounded monotone normal solve, in 1..=128 iterations.
    pub maximum_normal_iterations: usize,
    pub maximum_step_s: f64,
    /// Bounds both actual native endpoint speeds and their midpoint [m/s].
    pub maximum_slip_speed_m_s: f64,
    pub geometry_tolerance_m: f64,
    pub displacement_tolerance_m: f64,
    pub velocity_tolerance_m_s: f64,
    pub angular_velocity_tolerance_rad_s: f64,
    pub linear_momentum_tolerance_ns: f64,
    pub angular_momentum_tolerance_nms: f64,
    pub energy_tolerance_j: f64,
}

/// Actual native motion and assembled contact/body balance for one phase.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidMindlinReport {
    /// Normal, then longitudinal shear, with no surrogate effective mass.
    pub mobility_kg_inverse: [f64; 2],
    pub normal_average_force_n: f64,
    pub normal_iterations: usize,
    pub relative_velocity_before_m_s: [f64; 2],
    pub relative_velocity_after_m_s: [f64; 2],
    /// Actual normal/shear contact acceleration from the held external wrenches.
    pub driving_acceleration_m_s2: [f64; 2],
    /// Actual midpoint approach and shear from the two native endpoints [m].
    pub relative_displacement_m: [f64; 2],
    /// COM midpoint travel, exposed for the surrounding pose integrator [m].
    pub body_displacement_world_m: [Vec3; 2],
    /// Contact impulse on B, excluding the external wrench on either body.
    pub impulse_on_b_world_ns: Vec3,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Work of the contact impulses alone, evaluated on actual native motion.
    pub midpoint_impulse_work_j: f64,
    pub external_wrench_work_j: f64,
    /// Delta(K) minus actual contact and external work [J].
    pub work_energy_residual_j: f64,
    pub displacement_residual_m: [f64; 2],
    /// Delta(K + normal storage + shear storage) + heat - external work [J].
    pub total_energy_residual_j: f64,
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    pub external_linear_impulse_ns: Vec3,
    pub external_angular_impulse_nms: Vec3,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
}

/// Publish both original body states and source-bound contact history together.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidMindlinCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub contact: SourceMindlinPathStep,
    pub report: RigidMindlinReport,
}

/// No refusal changes either input body or accepted source/history.
#[derive(Debug)]
pub enum RigidMindlinError {
    Invalid(&'static str),
    Dynamics(DynamicsError),
    Normal(NormalPatchError),
    Contact(SourceMindlinPathError),
    /// The requested interval reaches or passes opening; locate/refine it.
    Separation,
    Unresolved(&'static str),
    NormalSolve {
        residual_m: f64,
        iterations: usize,
    },
    Cancelled,
}
impl From<DynamicsError> for RigidMindlinError {
    fn from(value: DynamicsError) -> Self {
        Self::Dynamics(value)
    }
}
impl From<NormalPatchError> for RigidMindlinError {
    fn from(value: NormalPatchError) -> Self {
        Self::Normal(value)
    }
}
impl From<SourceMindlinPathError> for RigidMindlinError {
    fn from(value: SourceMindlinPathError) -> Self {
        Self::Contact(value)
    }
}
impl core::fmt::Display for RigidMindlinError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native Hertz--Mindlin phase: {self:?}")
    }
}
impl std::error::Error for RigidMindlinError {}

/// Resolve normal loading and retained shear, then apply their combined impulse.
///
/// No external force is inferred; the enclosing composition owns driving and
/// other forces. Both average contact efforts act oppositely on the same native
/// material points. Initial/preloaded storage is preserved, and elastic shear
/// can return energy to the bodies rather than being mislabeled friction heat.
pub fn advance_rigid_mindlin(
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    contact: RigidMindlinContact<'_>,
    history: &SourceHertzMindlinState,
    dt_s: f64,
    controls: RigidMindlinControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RigidMindlinCandidate, RigidMindlinError> {
    advance_rigid_mindlin_driven(
        bodies,
        properties,
        contact,
        history,
        [Wrench::ZERO; 2],
        dt_s,
        controls,
        cancelled,
    )
}

/// Resolve held external drive, Hertz normal loading and retained Mindlin shear.
///
/// External COM forces are world-frame vectors; free torques are body-frame
/// vectors, in original A/B order. Native acceleration enters both released
/// displacements before the original source solves. The resulting external and
/// contact impulses are applied together once per body, without a pre-kick.
/// Full-body midpoint external work is separate from reversible contact storage
/// and physical microslip heat, including when a preloaded contact returns work.
///
/// Geometry, arms, orientations and wrenches are held for this momentum phase.
/// The source paths are affine in indentation and solved shear. Pose transport,
/// gyroscopic evolution, curved trajectories and separation/recontact retain
/// their surrounding owners; unsupported motion and source domains still refuse.
#[allow(clippy::too_many_arguments)]
pub fn advance_rigid_mindlin_driven(
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    contact: RigidMindlinContact<'_>,
    history: &SourceHertzMindlinState,
    external: [Wrench; 2],
    dt_s: f64,
    controls: RigidMindlinControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RigidMindlinCandidate, RigidMindlinError> {
    poll(cancelled)?;
    admit(controls, dt_s)?;
    let stored_before = finite(contact.law.stored_energy_j(history))?;
    if stored_before < 0.0 || history.indentation_m() <= 0.0 {
        return Err(RigidMindlinError::Invalid(
            "original loaded source history required",
        ));
    }
    let source = contact.law.fixed_source().normal_source();
    if !matches!(source.geometry(), PointContactGeometry::SpherePlane { .. }) {
        return Err(RigidMindlinError::Invalid(
            "original circular Hertz geometry required",
        ));
    }
    source.validate_kinematics(history.indentation_m(), 0.0, dt_s)?;
    let [normal, tangent, lateral] = [
        contact.frame.normal(),
        contact.frame.longitudinal(),
        contact.frame.lateral(),
    ]
    .map(|v| Vec3::new(v[0], v[1], v[2]));
    let points = [
        bodies[0].point_kinematics(properties[0], contact.arms_body[0])?,
        bodies[1].point_kinematics(properties[1], contact.arms_body[1])?,
    ];
    if norm(points[0].point_world.sub(points[1].point_world))? > controls.geometry_tolerance_m {
        return Err(RigidMindlinError::Invalid(
            "common point and central normal force lines required",
        ));
    }
    for point in points {
        if norm(point.arm_world.cross(normal))? > controls.geometry_tolerance_m {
            return Err(RigidMindlinError::Invalid(
                "central normal force lines required",
            ));
        }
    }
    if points[0].arm_world.dot(normal) < 0.0
        || points[1].arm_world.dot(normal) > 0.0
        || bodies[1]
            .pose()
            .position_world()
            .sub(bodies[0].pose().position_world())
            .dot(normal)
            <= 0.0
    {
        return Err(RigidMindlinError::Invalid(
            "opposed central attachments and A-to-B normal required",
        ));
    }
    let before_relative = points[0]
        .point_velocity_world
        .sub(points[1].point_velocity_world);
    admit_unmodeled(
        before_relative,
        points[0]
            .angular_velocity_world
            .sub(points[1].angular_velocity_world),
        normal,
        lateral,
        controls,
    )?;
    let before = [
        finite(before_relative.dot(normal))?,
        finite(before_relative.dot(tangent))?,
    ];
    source.validate_kinematics(history.indentation_m(), before[0], dt_s)?;
    if before[1].abs() > controls.maximum_slip_speed_m_s {
        return Err(RigidMindlinError::Unresolved(
            "initial native shear speed exceeds admitted phase",
        ));
    }
    let axes = [normal, tangent];
    let mut mobility = [[0.0; 2]; 2];
    for (column, direction) in axes.into_iter().enumerate() {
        let mut response = Vec3::ZERO;
        for i in 0..2 {
            let orientation = bodies[i].pose().orientation();
            let angular = properties[i].angular_velocity_body_checked(
                contact.arms_body[i].cross(orientation.rotate_world_to_body(direction)),
            )?;
            response = response
                .add(direction.scale(properties[i].mass().recip()))
                .add(
                    orientation
                        .rotate_body_to_world(angular)
                        .cross(points[i].arm_world),
                );
        }
        for (row, axis) in axes.into_iter().enumerate() {
            mobility[row][column] = finite(axis.dot(response))?;
        }
    }
    if mobility[0][0] <= 0.0 || mobility[1][1] <= 0.0 {
        return Err(RigidMindlinError::Invalid(
            "positive native normal and shear mobility required",
        ));
    }
    // Centrality makes these entries mathematically zero. Admit only reduction
    // roundoff, never discard a physically coupled normal/shear response.
    let cross_roundoff = 256.0 * f64::EPSILON * mobility[0][0].sqrt() * mobility[1][1].sqrt();
    if mobility[0][1].abs() > cross_roundoff || mobility[1][0].abs() > cross_roundoff {
        return Err(RigidMindlinError::Invalid(
            "normal/shear mechanical coupling requires its own solver",
        ));
    }
    let compliance = [
        finite(0.5 * dt_s * dt_s * mobility[0][0])?,
        finite(0.5 * dt_s * dt_s * mobility[1][1])?,
    ];
    let drive = HeldDrive::new(bodies, properties, points, axes, external)?;
    let mut free = [finite(dt_s * before[0])?, finite(dt_s * before[1])?];
    if external != [Wrench::ZERO; 2] {
        let half_dt_squared = finite(0.5 * dt_s * dt_s)?;
        for axis in 0..2 {
            free[axis] = finite(free[axis] + half_dt_squared * drive.acceleration[axis])?;
        }
    }
    let (indentation, normal_path, normal_iterations) = solve_normal(
        source,
        history.indentation_m(),
        free[0],
        compliance[0],
        dt_s,
        controls,
        cancelled,
    )?;
    let path = contact.law.advance_coupled_path(
        history,
        &SourceMindlinPathInput {
            indentation_m: indentation,
            free_slip_increment_m: free[1],
            dt_s,
            compliance_m_per_n: compliance[1],
            solve: controls.solve,
        },
        cancelled,
    )?;
    poll(cancelled)?;
    let force = path.tangential.traction_n[0];
    let impulse = normal
        .scale(finite(dt_s * normal_path.mean_force_n)?)
        .add(tangent.scale(finite(dt_s * force)?));
    let native = drive.apply(
        bodies,
        properties,
        contact.arms_body,
        points,
        impulse,
        dt_s,
        cancelled,
    )?;
    let after_bodies = native.bodies;
    let after_points = native.points;
    let after_relative = after_points[0]
        .point_velocity_world
        .sub(after_points[1].point_velocity_world);
    admit_unmodeled(
        after_relative,
        after_points[0]
            .angular_velocity_world
            .sub(after_points[1].angular_velocity_world),
        normal,
        lateral,
        controls,
    )?;
    let after = [
        finite(after_relative.dot(normal))?,
        finite(after_relative.dot(tangent))?,
    ];
    source.validate_kinematics(indentation, after[0], dt_s)?;
    if after[1].abs() > controls.maximum_slip_speed_m_s {
        return Err(RigidMindlinError::Unresolved(
            "final native shear speed exceeds admitted phase",
        ));
    }
    let travel = [
        finite(dt_s * (0.5 * before[0] + 0.5 * after[0]))?,
        finite(dt_s * (0.5 * before[1] + 0.5 * after[1]))?,
    ];
    let displacement_residual = [
        travel[0] - (indentation - history.indentation_m()),
        travel[1] - path.tangential.slip_increment_m[0],
    ];
    for axis in 0..2 {
        audit(
            displacement_residual[axis],
            controls.displacement_tolerance_m,
            "actual native midpoint displacement",
        )?;
        let mean = if axis == 0 {
            normal_path.mean_force_n
        } else {
            force
        };
        audit(
            after[axis]
                - (before[axis] - dt_s * mobility[axis][axis] * mean
                    + dt_s * drive.acceleration[axis]),
            controls.velocity_tolerance_m_s,
            "actual native endpoint velocity",
        )?;
    }
    if finite(travel[1] / dt_s)?.abs() > controls.maximum_slip_speed_m_s {
        return Err(RigidMindlinError::Unresolved(
            "native slip speed exceeds admitted phase",
        ));
    }
    let kinetic_before = native.kinetic_before_j;
    let kinetic_after = native.kinetic_after_j;
    let work_residual =
        finite(native.kinetic_energy_change_j - native.contact_work_j - native.external_work_j)?;
    let shear_change = finite(contact.law.stored_energy_j(&path.tangential.state) - stored_before)?;
    let energy_residual = finite(
        native.kinetic_energy_change_j
            + path.normal_work_j
            + shear_change
            + path.tangential.dissipated_j
            - native.external_work_j,
    )?;
    for (residual, what) in [
        (native.work_energy_residual_j, "native impulse work"),
        (work_residual, "native contact and external work"),
        (
            normal_path.work_j - path.normal_work_j,
            "original normal potential work",
        ),
        (
            normal_path.mean_force_n * travel[0] - path.normal_work_j,
            "actual normal work",
        ),
        (
            shear_change + path.tangential.dissipated_j - force * travel[1],
            "actual shear work",
        ),
        (
            shear_change - path.tangential.stored_change_j,
            "retained shear storage",
        ),
        (energy_residual, "assembled native contact/body energy"),
    ] {
        audit(residual, controls.energy_tolerance_j, what)?;
    }
    let linear = native.linear_change_ns;
    let angular = angular_momentum(after_bodies).sub(angular_momentum(bodies));
    let linear_residual = linear.sub(native.external_linear_ns);
    let angular_residual = angular.sub(native.external_angular_nms);
    audit(
        norm(linear_residual)?,
        controls.linear_momentum_tolerance_ns,
        "native pair linear momentum",
    )?;
    audit(
        norm(angular_residual)?,
        controls.angular_momentum_tolerance_nms,
        "native pair angular momentum",
    )?;
    let mut body_displacement = [Vec3::ZERO; 2];
    for i in 0..2 {
        body_displacement[i] = points[i]
            .center_of_mass_velocity_world
            .scale(0.5 * dt_s)
            .add(
                after_points[i]
                    .center_of_mass_velocity_world
                    .scale(0.5 * dt_s),
            );
        norm(body_displacement[i])?;
    }
    poll(cancelled)?;
    Ok(RigidMindlinCandidate {
        bodies_after: after_bodies,
        contact: path,
        report: RigidMindlinReport {
            mobility_kg_inverse: [mobility[0][0], mobility[1][1]],
            normal_average_force_n: normal_path.mean_force_n,
            normal_iterations,
            relative_velocity_before_m_s: before,
            relative_velocity_after_m_s: after,
            driving_acceleration_m_s2: drive.acceleration,
            relative_displacement_m: travel,
            body_displacement_world_m: body_displacement,
            impulse_on_b_world_ns: impulse,
            kinetic_before_j: kinetic_before,
            kinetic_after_j: kinetic_after,
            midpoint_impulse_work_j: native.contact_work_j,
            external_wrench_work_j: native.external_work_j,
            work_energy_residual_j: work_residual,
            displacement_residual_m: displacement_residual,
            total_energy_residual_j: energy_residual,
            linear_momentum_change_ns: linear,
            angular_momentum_change_nms: angular,
            external_linear_impulse_ns: native.external_linear_ns,
            external_angular_impulse_nms: native.external_angular_nms,
            linear_momentum_residual_ns: linear_residual,
            angular_momentum_residual_nms: angular_residual,
        },
    })
}

fn admit(controls: RigidMindlinControls, dt: f64) -> Result<(), RigidMindlinError> {
    if !(1..=128).contains(&controls.maximum_normal_iterations)
        || [
            dt,
            controls.maximum_step_s,
            controls.maximum_slip_speed_m_s,
            controls.geometry_tolerance_m,
            controls.displacement_tolerance_m,
            controls.velocity_tolerance_m_s,
            controls.angular_velocity_tolerance_rad_s,
            controls.linear_momentum_tolerance_ns,
            controls.angular_momentum_tolerance_nms,
            controls.energy_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
        || dt > controls.maximum_step_s
    {
        return Err(RigidMindlinError::Invalid(
            "finite positive native phase gates and bounded normal solve required",
        ));
    }
    Ok(())
}
fn admit_unmodeled(
    relative: Vec3,
    spin: Vec3,
    normal: Vec3,
    lateral: Vec3,
    controls: RigidMindlinControls,
) -> Result<(), RigidMindlinError> {
    if finite(relative.dot(lateral))?.abs() > controls.velocity_tolerance_m_s
        || finite(spin.dot(normal))?.abs() > controls.angular_velocity_tolerance_rad_s
    {
        return Err(RigidMindlinError::Invalid(
            "scalar Mindlin refuses transverse shear or relative torsional slip",
        ));
    }
    Ok(())
}
fn angular_momentum(bodies: [RigidBodyState; 2]) -> Vec3 {
    bodies.into_iter().fold(Vec3::ZERO, |sum, body| {
        sum.add(
            body.pose()
                .position_world()
                .cross(body.linear_momentum_world()),
        )
        .add(
            body.pose()
                .orientation()
                .rotate_body_to_world(body.angular_momentum_body()),
        )
    })
}
fn finite(value: f64) -> Result<f64, RigidMindlinError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RigidMindlinError::Invalid(
            "unrepresentable native contact arithmetic",
        ))
    }
}
fn norm(value: Vec3) -> Result<f64, RigidMindlinError> {
    Ok(value.stable_norm("native Mindlin vector")?)
}
fn audit(value: f64, allowance: f64, what: &'static str) -> Result<(), RigidMindlinError> {
    if finite(value)?.abs() <= allowance {
        Ok(())
    } else {
        Err(RigidMindlinError::Unresolved(what))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RigidMindlinError> {
    if cancelled() {
        Err(RigidMindlinError::Cancelled)
    } else {
        Ok(())
    }
}
