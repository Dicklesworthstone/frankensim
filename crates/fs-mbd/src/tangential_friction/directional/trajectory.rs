//! Original finite-body recoil from continuous directional friction relaxation.
//!
//! This impulse phase retains the original frozen poses, arms and normal load.
//! It integrates the native associated kinetic law through curved sliding and
//! finite arrest, then applies its integrated opposing impulses to the bodies.
//! Returned travel and heat belong to the complete curve, not a midpoint chord.
//! Normal constraints, pose transport, changing sources and thermal storage
//! remain with their surrounding owners.

use super::*;
use fs_tribo::tangential::elliptic::relaxation::{
    EllipticRelaxationControls, EllipticRelaxationError, EllipticRelaxationInput,
    EllipticRelaxationReport,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidDirectionalTrajectoryAudit {
    pub body_linear_momentum_residual_ns: [Vec3; 2],
    pub body_angular_momentum_residual_nms: [Vec3; 2],
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    pub endpoint_velocity_residual_m_s: f64,
    pub work_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
}

#[derive(Clone, Debug)]
pub struct RigidDirectionalTrajectoryReport {
    pub relaxation: EllipticRelaxationReport,
    pub mobility_per_kg: [[f64; 2]; 2],
    pub slip_velocity_before_m_s: [f64; 2],
    pub slip_velocity_after_m_s: [f64; 2],
    pub normal_velocity_before_m_s: f64,
    pub normal_velocity_after_m_s: f64,
    pub force_impulse_on_a_world_ns: Vec3,
    pub kinetic_energy_change_j: f64,
    pub impulse_work_j: f64,
    pub audit: RigidDirectionalTrajectoryAudit,
}

#[derive(Debug)]
pub enum RigidDirectionalTrajectoryError<E = EllipticRelaxationError> {
    Mechanics(RigidDirectionalError),
    Relaxation(E),
    Audit(RigidDirectionalTrajectoryAudit),
    Cancelled,
}
impl<E> From<DynamicsError> for RigidDirectionalTrajectoryError<E> {
    fn from(e: DynamicsError) -> Self {
        Self::Mechanics(RigidDirectionalError::Dynamics(e))
    }
}
impl<E> From<RigidDirectionalError> for RigidDirectionalTrajectoryError<E> {
    fn from(e: RigidDirectionalError) -> Self {
        Self::Mechanics(e)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidDirectionalTrajectoryError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid directional trajectory: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for RigidDirectionalTrajectoryError<E> {}

/// Advance a constant-law, held-load two-body friction phase through arrest.
/// Exactly two bodies and the original history publish only after independent
/// impulse, endpoint-velocity and work/heat checks and a final cancellation poll.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_directional_trajectory(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidDirectionalContact,
    law: &EllipticCoulomb,
    history: &mut EllipticCoulombState,
    dt_s: f64,
    controls: RigidDirectionalControls,
    relaxation: EllipticRelaxationControls,
    gate: &CancelGate,
) -> Result<RigidDirectionalTrajectoryReport, RigidDirectionalTrajectoryError> {
    apply_rigid_directional_trajectory_with(
        bodies,
        properties,
        contact,
        history,
        dt_s,
        controls,
        |state, input| {
            let report =
                law.relax_velocity(state, input, relaxation, &mut || gate.is_requested())?;
            let next = report.state.clone();
            Ok((report, next))
        },
        gate,
    )
}

/// Apply a source owner's continuous trajectory to the original body pair.
///
/// The callback receives the actual full body mobility and initial velocity,
/// and returns its native trajectory plus its own source/history candidate.
/// It may resolve changing coefficient queries without owning another body
/// response or applying impulses. This owner independently audits endpoint
/// velocity, opposing impulses and physical heat before publishing exactly two
/// bodies and the callback history together. Callback refusals are preserved.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn apply_rigid_directional_trajectory_with<S, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidDirectionalContact,
    history: &mut S,
    dt_s: f64,
    controls: RigidDirectionalControls,
    mut resolve: impl FnMut(&S, &EllipticRelaxationInput) -> Result<(EllipticRelaxationReport, S), E>,
    gate: &CancelGate,
) -> Result<RigidDirectionalTrajectoryReport, RigidDirectionalTrajectoryError<E>> {
    if gate.is_requested() {
        return Err(RigidDirectionalTrajectoryError::Cancelled);
    }
    let binding = contact.binding;
    if bodies.len() != properties.len()
        || bodies.len() > controls.maximum_bodies
        || binding.body_a >= bodies.len()
        || binding.body_b >= bodies.len()
        || binding.body_a == binding.body_b
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || [
            controls.geometry_tolerance_m,
            controls.velocity_tolerance_m_s,
            controls.linear_momentum_tolerance_ns,
            controls.angular_momentum_tolerance_nms,
            controls.energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(RigidDirectionalError::Invalid(
            "bounded original pair and finite physical gates required",
        )
        .into());
    }
    let before = [bodies[binding.body_a], bodies[binding.body_b]];
    let props = [properties[binding.body_a], properties[binding.body_b]];
    let arms = [binding.arm_a_body, binding.arm_b_body];
    let axes = [
        vector(contact.frame.longitudinal()),
        vector(contact.frame.lateral()),
        vector(contact.frame.normal()),
    ];
    let old = points(before, props, arms)?;
    if magnitude(old[0].point_world.sub(old[1].point_world)) > controls.geometry_tolerance_m {
        return Err(
            RigidDirectionalError::Invalid("original contact points do not coincide").into(),
        );
    }
    let mobility = response::<TriboError>(before, props, arms, axes)?;
    let initial = old[0].point_velocity_world.sub(old[1].point_velocity_world);
    let velocity = [initial.dot(axes[0]), initial.dot(axes[1])];
    let (resolved, next) = resolve(
        history,
        &EllipticRelaxationInput {
            velocity_m_s: velocity,
            mobility_per_kg: mobility,
            normal_n: contact.normal_load_n,
            dt_s,
            maximum_force_n: controls.contact.maximum_force_n,
            contact_temperature_k: contact.contact_temperature_k,
            nominal_pressure_pa: contact.nominal_pressure_pa,
        },
    )
    .map_err(RigidDirectionalTrajectoryError::Relaxation)?;
    if !resolved.distance_m.is_finite()
        || resolved.distance_m < 0.0
        || !resolved.dissipated_j.is_finite()
        || resolved.dissipated_j < 0.0
        || resolved.slip_increment_m.iter().any(|x| !x.is_finite())
    {
        return Err(RigidDirectionalError::Invalid(
            "finite nonnegative native trajectory distance and heat required",
        )
        .into());
    }
    let impulse = axes[0]
        .scale(-resolved.traction_impulse_ns[0])
        .add(axes[1].scale(-resolved.traction_impulse_ns[1]));
    let mut after = before;
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        after[i] = before[i]
            .apply_impulse_at_body_point(props[i], arms[i], impulse.scale(sign))?
            .state_after;
    }
    let new = points(after, props, arms)?;
    let final_rate = new[0].point_velocity_world.sub(new[1].point_velocity_world);
    let end = [final_rate.dot(axes[0]), final_rate.dot(axes[1])];
    let mut kinetic = 0.0;
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    let mut body_linear = [Vec3::ZERO; 2];
    let mut body_angular = [Vec3::ZERO; 2];
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let dp = after[i]
            .linear_momentum_world()
            .sub(before[i].linear_momentum_world());
        let dl = after[i]
            .angular_momentum_body()
            .sub(before[i].angular_momentum_body());
        let mean_v = old[i]
            .center_of_mass_velocity_world
            .scale(0.5)
            .add(new[i].center_of_mass_velocity_world.scale(0.5));
        let mean_w = old[i]
            .angular_velocity_body
            .scale(0.5)
            .add(new[i].angular_velocity_body.scale(0.5));
        kinetic += dp.dot(mean_v) + dl.dot(mean_w);
        let expected = impulse.scale(sign);
        body_linear[i] = dp.sub(expected);
        body_angular[i] = dl.sub(
            arms[i].cross(
                before[i]
                    .pose()
                    .orientation()
                    .rotate_world_to_body(expected),
            ),
        );
        linear = linear.add(dp);
        angular = angular
            .add(before[i].pose().position_world().cross(dp))
            .add(before[i].pose().orientation().rotate_body_to_world(dl));
    }
    // For a constant symmetric body mobility the integrated impulse's endpoint
    // pairing equals kinetic change even when instantaneous force turns. It is
    // independent of the net displacement integral and does not replace it.
    let impulse_work = impulse.dot(initial.scale(0.5).add(final_rate.scale(0.5)));
    let audit = RigidDirectionalTrajectoryAudit {
        body_linear_momentum_residual_ns: body_linear,
        body_angular_momentum_residual_nms: body_angular,
        linear_momentum_change_ns: linear,
        angular_momentum_change_nms: angular,
        endpoint_velocity_residual_m_s: det::hypot(
            end[0] - resolved.velocity_after_m_s[0],
            end[1] - resolved.velocity_after_m_s[1],
        ),
        work_energy_residual_j: kinetic - impulse_work,
        total_energy_residual_j: kinetic + resolved.dissipated_j,
    };
    if !audit.endpoint_velocity_residual_m_s.is_finite()
        || audit.endpoint_velocity_residual_m_s > controls.velocity_tolerance_m_s
        || [audit.work_energy_residual_j, audit.total_energy_residual_j]
            .into_iter()
            .any(|x| !x.is_finite() || x.abs() > controls.energy_tolerance_j)
        || body_linear
            .into_iter()
            .chain([linear])
            .any(|x| !x.is_finite() || magnitude(x) > controls.linear_momentum_tolerance_ns)
        || body_angular
            .into_iter()
            .chain([angular])
            .any(|x| !x.is_finite() || magnitude(x) > controls.angular_momentum_tolerance_nms)
    {
        return Err(RigidDirectionalTrajectoryError::Audit(audit));
    }
    let report = RigidDirectionalTrajectoryReport {
        relaxation: resolved,
        mobility_per_kg: mobility,
        slip_velocity_before_m_s: velocity,
        slip_velocity_after_m_s: end,
        normal_velocity_before_m_s: initial.dot(axes[2]),
        normal_velocity_after_m_s: final_rate.dot(axes[2]),
        force_impulse_on_a_world_ns: impulse,
        kinetic_energy_change_j: kinetic,
        impulse_work_j: impulse_work,
        audit,
    };
    if gate.is_requested() {
        return Err(RigidDirectionalTrajectoryError::Cancelled);
    }
    bodies[binding.body_a] = after[0];
    bodies[binding.body_b] = after[1];
    *history = next;
    Ok(report)
}
