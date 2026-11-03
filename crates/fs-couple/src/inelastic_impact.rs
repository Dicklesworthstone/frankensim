//! One frictionless rigid/material-point impact with complete solid history.
//!
//! Current witnesses and velocities determine the impulse through the original
//! fs-contact Newton law. The original consistent solid mass and full rigid
//! inertia supply its effective mass. Both native states are staged before any
//! result is returned; material deformation, memory and physical time stay fixed.
//!
//! The scalar functions below resolve one frictionless contact. [`multiple`]
//! resolves a complete simultaneous frictional manifold through the existing
//! coupled contact owner. Detection, event time, persistent contact and heat
//! partition remain with their owners. No penetration repair occurs.

/// Coupled frictional events retaining plastic, Maxwell and thermal histories.
pub mod multiple;

use fs_contact::impact::{
    NormalImpactError, NormalImpactInput, NormalImpactResult, solve_normal_impact,
};
use fs_exec::Cx;
use fs_material::plastic::finite::FiniteJ2Material;
use fs_mbd::{
    DynamicsError, MassProperties, MomentumImpulseReceipt, PointKinematics, RigidBodyState, Vec3,
};
use fs_solid::hyper3::plastic::{FinitePlasticTetError, dynamics::FinitePlasticTetMotion};
use fs_solid::hyper3::visco::dynamics::{
    SolidImpulseControls, SolidImpulseReport, SolidPointKinematics, SolidPointMobility,
    TetImpulsePoint,
};
use fs_solid::hyper3::visco::{FiniteViscoTetError, dynamics::FiniteViscoTetMotion};
use fs_solver::SolveReport;

/// Original material witnesses; the normal points from solid toward rigid.
#[derive(Clone, Copy, Debug)]
pub struct InelasticRigidContact {
    pub solid_point: TetImpulsePoint,
    /// Rigid COM-relative material arm [m].
    pub rigid_arm_body_m: Vec3,
    /// Nonzero world direction; the native rigid owner normalizes it once.
    pub normal_world: Vec3,
    pub activation_tolerance_m: f64,
    /// Authored system-level Newton coefficient in [0,1], not a bulk property.
    pub restitution: f64,
}

/// Separate native solve and observable physical acceptance limits.
#[derive(Clone, Copy, Debug)]
pub struct InelasticImpactControls {
    pub solid_impulse: SolidImpulseControls,
    pub geometry_tolerance_m: f64,
    pub velocity_tolerance_m_s: f64,
    pub linear_tolerance_kg_m_s: f64,
    pub angular_tolerance_kg_m2_s: f64,
    pub energy_tolerance_j: f64,
}

#[derive(Clone, Debug)]
pub struct InelasticImpactReport {
    pub contact: NormalImpactResult,
    pub gap_m: f64,
    pub inverse_effective_mass_per_kg: f64,
    pub solid_mobility_solve: SolveReport,
    pub solid_impulse: SolidImpulseReport,
    pub rigid_impulse: MomentumImpulseReceipt,
    pub relative_normal_speed_before_m_s: f64,
    pub relative_normal_speed_after_m_s: f64,
    pub endpoint_velocity_residual_m_s: f64,
    /// Actual native kinetic increments, without subtracting stored energies.
    pub kinetic_change_j: f64,
    /// Signed external work at held solid velocities; it need not vanish.
    pub support_work_j: f64,
    pub relative_impulse_work_j: f64,
    /// Kinetic change minus support work minus actual contact-point work.
    pub work_energy_residual_j: f64,
    /// Kinetic change minus support work plus the original contact-law loss.
    pub dissipation_residual_j: f64,
    /// Native total change minus the solid support reaction.
    pub linear_momentum_residual_kg_m_s: Vec3,
    pub angular_momentum_residual_kg_m2_s: Vec3,
    /// Body-local defects are summed in magnitude so partners cannot hide them.
    pub absolute_linear_residual_kg_m_s: f64,
    pub absolute_angular_residual_kg_m2_s: f64,
    pub absolute_work_residual_j: f64,
}

#[derive(Debug)]
pub struct InelasticImpactCandidate<M> {
    pub solid: M,
    pub rigid: RigidBodyState,
    pub report: InelasticImpactReport,
}

/// Original solid/material failures retain their concrete source type.
#[derive(Debug)]
pub enum InelasticImpactError<E> {
    Invalid(&'static str),
    Solid(E),
    Rigid(DynamicsError),
    Contact(NormalImpactError),
    Cancelled,
    Audit(Box<InelasticImpactReport>),
}
impl<E: core::fmt::Debug> core::fmt::Display for InelasticImpactError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "inelastic rigid impact: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for InelasticImpactError<E> {}

/// Stage one rigid/plastic event using original consistent mass and history.
/// Held DOFs retain their CURRENT velocity, including nonzero driven supports.
/// Dropping a candidate or any refusal leaves both accepted states intact.
#[allow(clippy::too_many_arguments)]
pub fn plastic_rigid_impact_candidate<'a, L: FiniteJ2Material>(
    solid: &FinitePlasticTetMotion<'a, L>,
    held_velocity_dofs: &[usize],
    rigid: RigidBodyState,
    properties: MassProperties,
    contact: InelasticRigidContact,
    origin_m: [f64; 3],
    controls: InelasticImpactControls,
    cx: &Cx<'_>,
) -> Result<
    InelasticImpactCandidate<FinitePlasticTetMotion<'a, L>>,
    InelasticImpactError<FinitePlasticTetError>,
> {
    let prepared = prepare(rigid, properties, contact, origin_m, controls, cx)?;
    let before = solid
        .point_kinematics(contact.solid_point, cx)
        .map_err(InelasticImpactError::Solid)?;
    let mobility = solid
        .point_directional_inverse_mass_per_kg(
            held_velocity_dofs,
            contact.solid_point,
            array(prepared.normal),
            controls.solid_impulse,
            cx,
        )
        .map_err(InelasticImpactError::Solid)?;
    finish(
        prepared,
        before,
        mobility,
        |impulse| {
            let event = solid.point_impulse_candidate(
                held_velocity_dofs,
                contact.solid_point,
                impulse,
                origin_m,
                controls.solid_impulse,
                cx,
            )?;
            let after = event.motion.point_kinematics(contact.solid_point, cx)?;
            Ok((event.motion, event.impulse, after))
        },
        cx,
    )
}

/// Same contact law and coupled audit for every retained Maxwell branch.
/// No constitutive step or elapsed material time is added by the impact.
#[allow(clippy::too_many_arguments)]
pub fn visco_rigid_impact_candidate<'a>(
    solid: &FiniteViscoTetMotion<'a>,
    held_velocity_dofs: &[usize],
    rigid: RigidBodyState,
    properties: MassProperties,
    contact: InelasticRigidContact,
    origin_m: [f64; 3],
    controls: InelasticImpactControls,
    cx: &Cx<'_>,
) -> Result<
    InelasticImpactCandidate<FiniteViscoTetMotion<'a>>,
    InelasticImpactError<FiniteViscoTetError>,
> {
    let prepared = prepare(rigid, properties, contact, origin_m, controls, cx)?;
    let before = solid
        .point_kinematics(contact.solid_point, cx)
        .map_err(InelasticImpactError::Solid)?;
    let mobility = solid
        .point_directional_inverse_mass_per_kg(
            held_velocity_dofs,
            contact.solid_point,
            array(prepared.normal),
            controls.solid_impulse,
            cx,
        )
        .map_err(InelasticImpactError::Solid)?;
    finish(
        prepared,
        before,
        mobility,
        |impulse| {
            let event = solid.point_impulse_candidate(
                held_velocity_dofs,
                contact.solid_point,
                impulse,
                origin_m,
                controls.solid_impulse,
                cx,
            )?;
            let after = event.motion.point_kinematics(contact.solid_point, cx)?;
            Ok((event.motion, event.impulse, after))
        },
        cx,
    )
}

struct Prepared {
    rigid: RigidBodyState,
    properties: MassProperties,
    contact: InelasticRigidContact,
    origin: Vec3,
    controls: InelasticImpactControls,
    normal: Vec3,
    rigid_point: PointKinematics,
    rigid_inverse_mass: f64,
}
fn prepare<E>(
    rigid: RigidBodyState,
    properties: MassProperties,
    contact: InelasticRigidContact,
    origin: [f64; 3],
    controls: InelasticImpactControls,
    cx: &Cx<'_>,
) -> Result<Prepared, InelasticImpactError<E>> {
    poll(cx)?;
    if !vector(origin).is_finite()
        || [
            controls.geometry_tolerance_m,
            controls.velocity_tolerance_m_s,
            controls.linear_tolerance_kg_m_s,
            controls.angular_tolerance_kg_m2_s,
            controls.energy_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(InelasticImpactError::Invalid(
            "finite origin and positive physical tolerances required",
        ));
    }
    let rigid_point = rigid
        .point_kinematics(properties, contact.rigid_arm_body_m)
        .map_err(InelasticImpactError::Rigid)?;
    let response = rigid
        .directional_effective_mass(properties, contact.rigid_arm_body_m, contact.normal_world)
        .map_err(InelasticImpactError::Rigid)?;
    Ok(Prepared {
        rigid,
        properties,
        contact,
        origin: vector(origin),
        controls,
        normal: response.unit_direction_world,
        rigid_point,
        rigid_inverse_mass: response.inverse_mass_kg_inverse,
    })
}

fn finish<M, E>(
    p: Prepared,
    before: SolidPointKinematics,
    mobility: SolidPointMobility,
    stage: impl FnOnce([f64; 3]) -> Result<(M, SolidImpulseReport, SolidPointKinematics), E>,
    cx: &Cx<'_>,
) -> Result<InelasticImpactCandidate<M>, InelasticImpactError<E>> {
    poll(cx)?;
    let separation = p.rigid_point.point_world.sub(vector(before.position_m));
    let gap = separation.dot(p.normal);
    if magnitude(separation.sub(p.normal.scale(gap))) > p.controls.geometry_tolerance_m {
        return Err(InelasticImpactError::Invalid(
            "original contact witnesses must share the normal line",
        ));
    }
    let relative_before = p
        .rigid_point
        .point_velocity_world
        .sub(vector(before.velocity_mps))
        .dot(p.normal);
    let inverse_mass = p.rigid_inverse_mass + mobility.inverse_mass_per_kg;
    let contact = solve_normal_impact(NormalImpactInput {
        gap,
        activation_tolerance: p.contact.activation_tolerance_m,
        relative_normal_speed: relative_before,
        inverse_effective_mass: inverse_mass,
        restitution: p.contact.restitution,
    })
    .map_err(InelasticImpactError::Contact)?;
    let impulse = p.normal.scale(contact.normal_impulse);
    let (solid, solid_impulse, after) =
        stage(array(impulse.scale(-1.0))).map_err(InelasticImpactError::Solid)?;
    let impulse_body = p.rigid.pose().orientation().rotate_world_to_body(impulse);
    let rigid_impulse = p
        .rigid
        .apply_momentum_impulse(
            p.properties,
            impulse,
            p.contact.rigid_arm_body_m.cross(impulse_body),
        )
        .map_err(InelasticImpactError::Rigid)?;
    let rigid_after = rigid_impulse
        .state_after
        .point_kinematics(p.properties, p.contact.rigid_arm_body_m)
        .map_err(InelasticImpactError::Rigid)?;
    let relative_after = rigid_after
        .point_velocity_world
        .sub(vector(after.velocity_mps))
        .dot(p.normal);
    let relative_work = (0.5 * relative_before + 0.5 * relative_after) * contact.normal_impulse;
    let change = solid_impulse.kinetic_change_j + rigid_impulse.work.kinetic_energy_change_j;
    let support_work = solid_impulse.support_work_j;
    let arm = p.rigid.pose().position_world().sub(p.origin);
    let linear_residual = vector(solid_impulse.applied.linear_kg_m_s)
        .add(vector(solid_impulse.linear_residual_kg_m_s))
        .add(rigid_impulse.linear_increment_world);
    let angular_residual = vector(solid_impulse.applied.angular_kg_m2_s)
        .add(vector(solid_impulse.angular_residual_kg_m2_s))
        .add(arm.cross(rigid_impulse.linear_increment_world))
        .add(
            p.rigid
                .pose()
                .orientation()
                .rotate_body_to_world(rigid_impulse.angular_increment_body),
        );
    let absolute_linear = magnitude(vector(solid_impulse.linear_residual_kg_m_s))
        + magnitude(rigid_impulse.linear_residual_world);
    let absolute_angular = magnitude(vector(solid_impulse.angular_residual_kg_m2_s))
        + magnitude(arm.cross(rigid_impulse.linear_residual_world))
        + magnitude(rigid_impulse.angular_residual_body);
    let absolute_work =
        solid_impulse.work_residual_j.abs() + rigid_impulse.work.work_energy_residual_j.abs();
    let report = InelasticImpactReport {
        contact,
        gap_m: gap,
        inverse_effective_mass_per_kg: inverse_mass,
        solid_mobility_solve: mobility.solve,
        solid_impulse,
        rigid_impulse,
        relative_normal_speed_before_m_s: relative_before,
        relative_normal_speed_after_m_s: relative_after,
        endpoint_velocity_residual_m_s: relative_after - contact.relative_normal_speed_after,
        kinetic_change_j: change,
        support_work_j: support_work,
        relative_impulse_work_j: relative_work,
        work_energy_residual_j: change - support_work - relative_work,
        dissipation_residual_j: change - support_work + contact.dissipated_energy,
        linear_momentum_residual_kg_m_s: linear_residual,
        angular_momentum_residual_kg_m2_s: angular_residual,
        absolute_linear_residual_kg_m_s: absolute_linear,
        absolute_angular_residual_kg_m2_s: absolute_angular,
        absolute_work_residual_j: absolute_work,
    };
    if [
        change,
        support_work,
        relative_work,
        report.endpoint_velocity_residual_m_s,
        report.work_energy_residual_j,
        report.dissipation_residual_j,
        magnitude(linear_residual),
        magnitude(angular_residual),
        absolute_linear,
        absolute_angular,
        absolute_work,
    ]
    .iter()
    .any(|v| !v.is_finite())
        || report.endpoint_velocity_residual_m_s.abs() > p.controls.velocity_tolerance_m_s
        || magnitude(linear_residual) > p.controls.linear_tolerance_kg_m_s
        || absolute_linear > p.controls.linear_tolerance_kg_m_s
        || magnitude(angular_residual) > p.controls.angular_tolerance_kg_m2_s
        || absolute_angular > p.controls.angular_tolerance_kg_m2_s
        || absolute_work > p.controls.energy_tolerance_j
        || report.work_energy_residual_j.abs() > p.controls.energy_tolerance_j
        || report.dissipation_residual_j.abs() > p.controls.energy_tolerance_j
    {
        return Err(InelasticImpactError::Audit(Box::new(report)));
    }
    poll(cx)?;
    Ok(InelasticImpactCandidate {
        solid,
        rigid: rigid_impulse.state_after,
        report,
    })
}
fn poll<E>(cx: &Cx<'_>) -> Result<(), InelasticImpactError<E>> {
    cx.checkpoint().map_err(|_| InelasticImpactError::Cancelled)
}
fn vector(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn magnitude(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
