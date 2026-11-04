//! Rigid/nodal impacts without restarting a tensor-Prony solid's trajectory.
//!
//! The admitted solid supplies its own reference nodes, current displacement,
//! supported consistent-mass response and complete constitutive history. The
//! rigid owner supplies its full inertia. One original Newton contact law
//! determines equal/opposite impulses; both native endpoints and their actual
//! linear/angular momentum and work are checked before either is published.
//!
//! This is a single frictionless event at an original solid node. Detection,
//! event time and persistent contact remain with their owners. Contact loss is
//! reported separately; it is not deposited into Maxwell/Kelvin material heat.

use super::{InelasticRigidContact, array, magnitude, vector};
use fs_contact::impact::{
    NormalImpactError, NormalImpactInput, NormalImpactResult, solve_normal_impact,
};
use fs_exec::Cx;
use fs_math::dd::Dd;
use fs_mbd::{DynamicsError, MassProperties, MomentumImpulseReceipt, RigidBodyState, Vec3};
use fs_solid::linear3::visco::{ViscoDynamicImpulse, ViscoDynamicState, ViscoElementError};

/// Independent coupled acceptance limits. The solid retains its prepared CG
/// and energy controls; this adapter cannot replace or relax them.
#[derive(Clone, Copy, Debug)]
pub struct LinearViscoImpactControls {
    pub geometry_tolerance_m: f64,
    pub velocity_tolerance_m_s: f64,
    pub linear_tolerance_kg_m_s: f64,
    pub angular_tolerance_kg_m2_s: f64,
    pub energy_tolerance_j: f64,
}

#[derive(Clone, Debug)]
pub struct LinearViscoImpactReport {
    pub contact: NormalImpactResult,
    pub gap_m: f64,
    pub inverse_effective_mass_per_kg: f64,
    pub solid_impulse: ViscoDynamicImpulse,
    pub rigid_impulse: MomentumImpulseReceipt,
    pub relative_normal_speed_before_m_s: f64,
    pub relative_normal_speed_after_m_s: f64,
    pub endpoint_velocity_residual_m_s: f64,
    /// Actual total kinetic change, less support work, plus contact loss [J].
    pub dissipation_residual_j: f64,
    /// Actual total kinetic change, less support and contact-point work [J].
    pub work_energy_residual_j: f64,
    /// Actual two-owner momentum change minus external support impulse.
    pub linear_momentum_residual_kg_m_s: Vec3,
    pub angular_momentum_residual_kg_m2_s: Vec3,
    /// Sum of native residual magnitudes: errors cannot cancel across owners.
    pub absolute_linear_residual_kg_m_s: f64,
    pub absolute_angular_residual_kg_m2_s: f64,
    pub absolute_work_residual_j: f64,
}

#[derive(Debug)]
pub enum LinearViscoImpactError {
    Invalid(&'static str),
    Solid(ViscoElementError),
    Rigid(DynamicsError),
    Contact(NormalImpactError),
    Cancelled,
    Audit(Box<LinearViscoImpactReport>),
}
impl core::fmt::Display for LinearViscoImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "linear viscoelastic rigid impact: {self:?}")
    }
}
impl std::error::Error for LinearViscoImpactError {}

/// Apply one normal impact using the original nodal trace and supported mass.
///
/// `contact.solid_point` is an original node index. Its physical location comes
/// from the solid's own admitted geometry and accepted displacement. The normal
/// points from solid to rigid; restitution is an explicit system contact law.
/// Supports retain their actual velocities, including driven nonzero rates.
/// All nodes and rigid coordinates use the same inertial Cartesian frame.
///
/// The event preserves the solid's complete material/thermal/source history,
/// time and interval budget. It changes only native momentum and its work
/// accounts. A native exclusive proposal stages only velocity and work, without
/// copying the solid's constitutive history. Any refusal or cancellation leaves
/// both input owners unchanged; their publication performs no further checks.
#[allow(clippy::too_many_arguments)]
pub fn apply_linear_visco_rigid_impact<'p, 'a>(
    solid: &mut ViscoDynamicState<'p, 'a>,
    rigid: &mut RigidBodyState,
    properties: MassProperties,
    contact: InelasticRigidContact<usize>,
    origin_m: [f64; 3],
    controls: LinearViscoImpactControls,
    cx: &Cx<'_>,
) -> Result<LinearViscoImpactReport, LinearViscoImpactError> {
    poll(cx)?;
    let origin = vector(origin_m);
    if !origin.is_finite()
        || [
            controls.geometry_tolerance_m,
            controls.velocity_tolerance_m_s,
            controls.linear_tolerance_kg_m_s,
            controls.angular_tolerance_kg_m2_s,
            controls.energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || contact.solid_point >= solid.reference_nodes_m().len()
    {
        return Err(LinearViscoImpactError::Invalid(
            "original node, finite origin and positive physical tolerances required",
        ));
    }
    let node = contact.solid_point;
    let point = position(solid, node)?;
    let before_velocity = nodal(solid.velocity_mps(), node);
    let rigid_point = rigid
        .point_kinematics(properties, contact.rigid_arm_body_m)
        .map_err(LinearViscoImpactError::Rigid)?;
    let mobility = rigid
        .directional_effective_mass(properties, contact.rigid_arm_body_m, contact.normal_world)
        .map_err(LinearViscoImpactError::Rigid)?;
    let normal = mobility.unit_direction_world;
    let separation = rigid_point.point_world.sub(point);
    let gap_m = separation.dot(normal);
    if !gap_m.is_finite()
        || magnitude(separation.sub(normal.scale(gap_m))) > controls.geometry_tolerance_m
    {
        return Err(LinearViscoImpactError::Invalid(
            "native contact witnesses must share the normal line",
        ));
    }
    let mut nodal_impulse = vec![0.0; solid.velocity_mps().len()];
    nodal_impulse[3 * node..3 * node + 3].copy_from_slice(&array(normal));
    let response = solid
        .momentum_response(&nodal_impulse, cx)
        .map_err(LinearViscoImpactError::Solid)?;
    let inverse_solid_mass = normal.dot(nodal(&response, node));
    if !inverse_solid_mass.is_finite() || inverse_solid_mass < 0.0 {
        return Err(LinearViscoImpactError::Invalid(
            "native supported contact mobility must be finite and nonnegative",
        ));
    }
    let relative_before = rigid_point
        .point_velocity_world
        .sub(before_velocity)
        .dot(normal);
    let inverse_effective_mass = inverse_solid_mass + mobility.inverse_mass_kg_inverse;
    let event = solve_normal_impact(NormalImpactInput {
        gap: gap_m,
        activation_tolerance: contact.activation_tolerance_m,
        relative_normal_speed: relative_before,
        inverse_effective_mass,
        restitution: contact.restitution,
    })
    .map_err(LinearViscoImpactError::Contact)?;
    let impulse = normal.scale(event.normal_impulse);
    nodal_impulse[3 * node..3 * node + 3].copy_from_slice(&array(impulse.scale(-1.0)));
    let proposal = solid
        .prepare_momentum_impulse(&nodal_impulse, controls.linear_tolerance_kg_m_s, cx)
        .map_err(LinearViscoImpactError::Solid)?;
    let solid_impulse = proposal.report();
    let solid = proposal.original();
    let rigid_impulse = rigid
        .apply_momentum_impulse(
            properties,
            impulse,
            contact
                .rigid_arm_body_m
                .cross(rigid.pose().orientation().rotate_world_to_body(impulse)),
        )
        .map_err(LinearViscoImpactError::Rigid)?;
    let rigid_after = rigid_impulse
        .state_after
        .point_kinematics(properties, contact.rigid_arm_body_m)
        .map_err(LinearViscoImpactError::Rigid)?;
    let relative_after = rigid_after
        .point_velocity_world
        .sub(nodal(proposal.velocity_mps(), node))
        .dot(normal);

    // Reduce actual full nodal momentum changes, removing only the original
    // support reactions. The complete mass coupling is included even though
    // the applied point impulse has only three nonzero entries.
    let mut linear = [Dd::ZERO; 3];
    let mut angular = [Dd::ZERO; 3];
    let mut absolute_linear = 0.0;
    let mut absolute_angular = 0.0;
    for i in 0..solid.reference_nodes_m().len() {
        if i % 256 == 0 {
            poll(cx)?;
        }
        let arm = position(solid, i)?.sub(origin);
        let r = array(arm).map(Dd::from_f64);
        let p: [Dd; 3] = core::array::from_fn(|d| {
            Dd::from_f64(solid_impulse.momentum_change_kg_m_s[3 * i + d])
                - Dd::from_f64(solid_impulse.support_impulse_kg_m_s[3 * i + d])
        });
        for d in 0..3 {
            linear[d] = linear[d] + p[d];
        }
        angular[0] = angular[0] + r[1] * p[2] - r[2] * p[1];
        angular[1] = angular[1] + r[2] * p[0] - r[0] * p[2];
        angular[2] = angular[2] + r[0] * p[1] - r[1] * p[0];
        let residual = nodal(&solid_impulse.impulse_residual_kg_m_s, i);
        absolute_linear += magnitude(residual);
        absolute_angular += magnitude(arm.cross(residual));
    }
    let rigid_arm = rigid.pose().position_world().sub(origin);
    // Keep orbital products and spin separate until the compensated world
    // reduction; large opposing terms must not round away a contact moment.
    let r = array(rigid_arm).map(Dd::from_f64);
    let p = array(rigid_impulse.linear_increment_world).map(Dd::from_f64);
    let rigid_orbital = [
        r[1] * p[2] - r[2] * p[1],
        r[2] * p[0] - r[0] * p[2],
        r[0] * p[1] - r[1] * p[0],
    ];
    let rigid_spin = array(
        rigid
            .pose()
            .orientation()
            .rotate_body_to_world(rigid_impulse.angular_increment_body),
    )
    .map(Dd::from_f64);
    let linear_residual = vector(core::array::from_fn(|d| {
        (linear[d] + Dd::from_f64(array(rigid_impulse.linear_increment_world)[d])).to_f64()
    }));
    let angular_residual = vector(core::array::from_fn(|d| {
        (angular[d] + rigid_orbital[d] + rigid_spin[d]).to_f64()
    }));
    absolute_linear += magnitude(rigid_impulse.linear_residual_world);
    absolute_angular += magnitude(rigid_arm.cross(rigid_impulse.linear_residual_world))
        + magnitude(rigid_impulse.angular_residual_body);
    let change = Dd::from_f64(solid_impulse.kinetic_change_j)
        + Dd::from_f64(rigid_impulse.work.kinetic_energy_change_j)
        - Dd::from_f64(solid_impulse.support_work_j);
    let relative_work = (Dd::from_f64(0.5 * relative_before) + Dd::from_f64(0.5 * relative_after))
        * Dd::from_f64(event.normal_impulse);
    let absolute_work =
        solid_impulse.balance_residual_j.abs() + rigid_impulse.work.work_energy_residual_j.abs();
    let endpoint_velocity_residual = relative_after - event.relative_normal_speed_after;
    let dissipation_residual = (change + Dd::from_f64(event.dissipated_energy)).to_f64();
    let work_energy_residual = (change - relative_work).to_f64();
    let report = |solid_impulse| LinearViscoImpactReport {
        contact: event,
        gap_m,
        inverse_effective_mass_per_kg: inverse_effective_mass,
        solid_impulse,
        rigid_impulse,
        relative_normal_speed_before_m_s: relative_before,
        relative_normal_speed_after_m_s: relative_after,
        endpoint_velocity_residual_m_s: endpoint_velocity_residual,
        dissipation_residual_j: dissipation_residual,
        work_energy_residual_j: work_energy_residual,
        linear_momentum_residual_kg_m_s: linear_residual,
        angular_momentum_residual_kg_m2_s: angular_residual,
        absolute_linear_residual_kg_m_s: absolute_linear,
        absolute_angular_residual_kg_m2_s: absolute_angular,
        absolute_work_residual_j: absolute_work,
    };
    if [
        (
            endpoint_velocity_residual.abs(),
            controls.velocity_tolerance_m_s,
        ),
        (magnitude(linear_residual), controls.linear_tolerance_kg_m_s),
        (
            magnitude(angular_residual),
            controls.angular_tolerance_kg_m2_s,
        ),
        (absolute_linear, controls.linear_tolerance_kg_m_s),
        (absolute_angular, controls.angular_tolerance_kg_m2_s),
        (absolute_work, controls.energy_tolerance_j),
        (dissipation_residual.abs(), controls.energy_tolerance_j),
        (work_energy_residual.abs(), controls.energy_tolerance_j),
    ]
    .iter()
    .any(|&(value, limit)| !value.is_finite() || value > limit)
    {
        return Err(LinearViscoImpactError::Audit(Box::new(report(
            solid_impulse.clone(),
        ))));
    }
    poll(cx)?;
    let report = report(proposal.publish());
    *rigid = rigid_impulse.state_after;
    Ok(report)
}

fn nodal(values: &[f64], node: usize) -> Vec3 {
    Vec3::new(values[3 * node], values[3 * node + 1], values[3 * node + 2])
}
fn position(
    solid: &ViscoDynamicState<'_, '_>,
    node: usize,
) -> Result<Vec3, LinearViscoImpactError> {
    let point = vector(solid.reference_nodes_m()[node])
        .add(nodal(solid.material_history().displacement_m(), node));
    if !point.is_finite() {
        return Err(LinearViscoImpactError::Invalid(
            "native contact position is not representable",
        ));
    }
    Ok(point)
}
fn poll(cx: &Cx<'_>) -> Result<(), LinearViscoImpactError> {
    cx.checkpoint()
        .map_err(|_| LinearViscoImpactError::Cancelled)
}
