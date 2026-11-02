//! Source-owned uniform-junction relaxation through a real stopping event.
//!
//! Actual masses, full inertia tensors, arms and tangent gauge determine the
//! original three-coordinate mobility. The junction owner admits an invariant
//! slip/spin ray and integrates its original traction to the exact finite stop;
//! any remaining no-drive interval has zero reaction and zero heat. Native
//! bodies receive only those integrated impulses, with no velocity assignment
//! or endpoint clipping. Their actual recoil and kinetic loss are checked again.
//!
//! This is a held-geometry impulse trajectory, not pose or normal integration.
//! Material/source changes, external driving and curved velocity trajectories
//! require their surrounding owners and are not approximated by this branch.

use super::*;
use fs_tribo::adhesive::spinning::relaxation::{
    AdhesiveRelaxationError, AdhesiveRelaxationInput, AdhesiveRelaxationReport,
};

/// Independent checks of original body impulses, endpoint motion and heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidAdhesiveTrajectoryAudit {
    /// Applied minus requested linear impulse for each original body [N s].
    pub body_linear_momentum_residual_ns: [Vec3; 2],
    /// Applied minus requested body-frame angular impulse for each body [N m s].
    pub body_angular_momentum_residual_nms: [Vec3; 2],
    /// Total pair linear-momentum change [N s].
    pub linear_momentum_change_ns: Vec3,
    /// Orbital plus spin angular-momentum change [N m s].
    pub angular_momentum_change_nms: Vec3,
    /// Actual endpoint velocity minus the analytic event target [m/s].
    pub endpoint_velocity_residual_m_s: f64,
    /// Actual endpoint spin minus the analytic event target [rad/s].
    pub endpoint_spin_residual_rad_s: f64,
    /// Largest initial/final normal rate, held at zero on an attached disk [m/s].
    pub normal_velocity_residual_m_s: f64,
    /// Actual midpoint native force minus the integrated constant phase force [N].
    pub force_residual_n: f64,
    /// Actual midpoint native moment minus the constant phase moment [N m].
    pub spin_moment_residual_nm: f64,
    /// Actual body kinetic change minus native impulse work [J].
    pub work_energy_residual_j: f64,
    /// Actual body kinetic change plus integrated original junction loss [J].
    pub total_energy_residual_j: f64,
    /// Original trajectory loss minus independently reevaluated native loss [J].
    pub heat_residual_j: f64,
    /// Positive assisting endpoint force/spin power times moving duration [J].
    pub maximum_assisting_work_j: f64,
}
impl RigidAdhesiveTrajectoryAudit {
    const ZERO: Self = Self {
        body_linear_momentum_residual_ns: [Vec3::ZERO; 2],
        body_angular_momentum_residual_nms: [Vec3::ZERO; 2],
        linear_momentum_change_ns: Vec3::ZERO,
        angular_momentum_change_nms: Vec3::ZERO,
        endpoint_velocity_residual_m_s: 0.0,
        endpoint_spin_residual_rad_s: 0.0,
        normal_velocity_residual_m_s: 0.0,
        force_residual_n: 0.0,
        spin_moment_residual_nm: 0.0,
        work_energy_residual_j: 0.0,
        total_energy_residual_j: 0.0,
        heat_residual_j: 0.0,
        maximum_assisting_work_j: 0.0,
    };
}

/// Physical trajectory and independently verified finite support recoil.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveTrajectoryReport {
    /// Original disk trajectory; absent for a detached zero-area source.
    pub relaxation: Option<AdhesiveRelaxationReport>,
    /// Complete actual body response in `[u_x,u_y,a*spin]`; absent if detached.
    pub scaled_mobility_per_kg: Option<[[f64; 3]; 3]>,
    /// Actual relative point velocity before the impulse phase [m/s].
    pub slip_velocity_before_m_s: [f64; 2],
    /// Actual relative point velocity after the impulse phase [m/s].
    pub slip_velocity_after_m_s: [f64; 2],
    /// Actual initial relative normal-axis angular velocity [rad/s].
    pub spin_before_rad_s: f64,
    /// Actual final relative normal-axis angular velocity [rad/s].
    pub spin_after_rad_s: f64,
    /// Original normal-owner rate before the impulse phase [m/s].
    pub normal_velocity_before_m_s: f64,
    /// Original normal-owner rate after the impulse phase [m/s].
    pub normal_velocity_after_m_s: f64,
    /// Actual equal/opposite point impulse acting on side A [N s].
    pub force_impulse_on_a_world_ns: Vec3,
    /// Actual equal/opposite free angular impulse acting on side A [N m s].
    pub spin_angular_impulse_on_a_world_nms: Vec3,
    /// Actual finite-pair kinetic change, computed without subtracting totals [J].
    pub kinetic_energy_change_j: f64,
    /// Work of the actual complete point and free-couple impulse [J].
    pub impulse_work_j: f64,
    /// Integrated local junction loss, unpartitioned and allocated at most once [J].
    pub dissipated_energy_j: f64,
    /// Independent original-body, actual-motion and original-source checks.
    pub audit: RigidAdhesiveTrajectoryAudit,
}

/// Complete staged pair; a surrounding owner can compose its own atomic publish.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveTrajectoryCandidate {
    /// Actual native momentum endpoints in the original binding's A/B order.
    pub bodies_after: [RigidBodyState; 2],
    /// Physical event, travel, impulse, heat and native-body evidence.
    pub report: RigidAdhesiveTrajectoryReport,
}

/// Source trajectory, actual-body mechanics and late audit refusals are preserved.
#[derive(Debug)]
pub enum RigidAdhesiveTrajectoryError {
    /// Invalid geometry, duration, body binding or dimensional controls.
    Invalid(&'static str),
    /// Original native kinematics, impulse or cancellation refusal.
    Mechanics(RigidHertzError),
    /// Original disk trajectory refused its physical or numerical admission.
    Relaxation(AdhesiveRelaxationError),
    /// The candidate actual bodies did not pass the original dimensional gates.
    Audit(RigidAdhesiveTrajectoryAudit),
}
impl From<RigidHertzError> for RigidAdhesiveTrajectoryError {
    fn from(error: RigidHertzError) -> Self {
        Self::Mechanics(error)
    }
}
impl From<DynamicsError> for RigidAdhesiveTrajectoryError {
    fn from(error: DynamicsError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl From<super::super::RigidFrictionError> for RigidAdhesiveTrajectoryError {
    fn from(error: super::super::RigidFrictionError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl From<AdhesiveRelaxationError> for RigidAdhesiveTrajectoryError {
    fn from(error: AdhesiveRelaxationError) -> Self {
        Self::Relaxation(error)
    }
}
impl From<AdhesiveWrenchError> for RigidAdhesiveTrajectoryError {
    fn from(error: AdhesiveWrenchError) -> Self {
        Self::Relaxation(error.into())
    }
}
impl core::fmt::Display for RigidAdhesiveTrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid adhesive trajectory: {self:?}")
    }
}
impl std::error::Error for RigidAdhesiveTrajectoryError {}

/// Resolve finite arrest and publish both original momenta after every audit.
pub fn apply_rigid_adhesive_trajectory(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidAdhesiveContact<'_>,
    dt_s: f64,
    controls: RigidAdhesiveControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveTrajectoryReport, RigidAdhesiveTrajectoryError> {
    let candidate =
        prepare_rigid_adhesive_trajectory(bodies, properties, contact, dt_s, controls, gate)?;
    poll(gate)?;
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    Ok(candidate.report)
}

/// Stage a continuous fixed-ray impulse phase from the unchanged original source.
/// Poses and contact area remain frozen; no body or source is changed on refusal.
#[allow(clippy::too_many_lines)]
pub fn prepare_rigid_adhesive_trajectory(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidAdhesiveContact<'_>,
    dt_s: f64,
    controls: RigidAdhesiveControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveTrajectoryCandidate, RigidAdhesiveTrajectoryError> {
    poll(gate)?;
    controls.contact.validate()?;
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
            controls.angular_velocity_tolerance_rad_s,
            controls.linear_momentum_tolerance_ns,
            controls.angular_momentum_tolerance_nms,
            controls.energy_tolerance_j,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x < 0.0)
    {
        return Err(RigidAdhesiveTrajectoryError::Invalid(
            "bounded distinct finite bodies, positive time and finite controls required",
        ));
    }
    let before = [bodies[binding.body_a], bodies[binding.body_b]];
    let properties = [properties[binding.body_a], properties[binding.body_b]];
    let arms = [binding.arm_a_body, binding.arm_b_body];
    let axes = [
        vector(contact.frame.longitudinal()),
        vector(contact.frame.lateral()),
        vector(contact.frame.normal()),
    ];
    let old = points(before, properties, arms)?;
    let initial = relative_point(old);
    let angular_initial = relative_spin(old);
    let velocity = [initial.dot(axes[0]), initial.dot(axes[1])];
    let spin = angular_initial.dot(axes[2]);
    let radius = contact.source.patch().parameters().contact_radius_m;
    if radius == 0.0 {
        poll(gate)?;
        return Ok(RigidAdhesiveTrajectoryCandidate {
            bodies_after: before,
            report: RigidAdhesiveTrajectoryReport {
                relaxation: None,
                scaled_mobility_per_kg: None,
                slip_velocity_before_m_s: velocity,
                slip_velocity_after_m_s: velocity,
                spin_before_rad_s: spin,
                spin_after_rad_s: spin,
                normal_velocity_before_m_s: initial.dot(axes[2]),
                normal_velocity_after_m_s: initial.dot(axes[2]),
                force_impulse_on_a_world_ns: Vec3::ZERO,
                spin_angular_impulse_on_a_world_nms: Vec3::ZERO,
                kinetic_energy_change_j: 0.0,
                impulse_work_j: 0.0,
                dissipated_energy_j: 0.0,
                audit: RigidAdhesiveTrajectoryAudit::ZERO,
            },
        });
    }
    if magnitude(old[0].point_world.sub(old[1].point_world)) > controls.geometry_tolerance_m
        || initial.dot(axes[2]).abs() > controls.velocity_tolerance_m_s
    {
        return Err(RigidAdhesiveTrajectoryError::Invalid(
            "attached material points must coincide with held normal motion",
        ));
    }
    let law = AdhesiveWrenchContact::new(AdhesiveWrenchParameters {
        contact_radius_m: radius,
        static_shear_strength_pa: contact.source.static_patch().parameters().shear_strength_pa,
        kinetic_shear_strength_pa: contact.source.patch().parameters().shear_strength_pa,
    })?;
    let mobility = scaled_response(before, properties, arms, axes, [0.0, 0.0, radius])?;
    let resolved = law.relax_velocity(
        &AdhesiveRelaxationInput {
            velocity_m_s: velocity,
            spin_rate_rad_s: spin,
            scaled_mobility_per_kg: mobility,
            dt_s,
        },
        controls.contact,
        || gate.is_requested(),
    )?;
    let impulse = axes[0]
        .scale(resolved.force_impulse_ns[0])
        .add(axes[1].scale(resolved.force_impulse_ns[1]));
    let angular_impulse = axes[2].scale(resolved.spin_angular_impulse_nms);
    let after = apply_pair_wrench(before, properties, arms, impulse, angular_impulse, gate)?;
    let new = points(after, properties, arms)?;
    let endpoint = relative_point(new);
    let angular_endpoint = relative_spin(new);
    let midpoint = initial.scale(0.5).add(endpoint.scale(0.5));
    let angular_midpoint = angular_initial.scale(0.5).add(angular_endpoint.scale(0.5));
    let native = contact
        .source
        .patch()
        .evaluate(
            [midpoint.dot(axes[0]), midpoint.dot(axes[1])],
            angular_midpoint.dot(axes[2]),
        )
        .map_err(AdhesiveWrenchError::Law)?;
    let native_heat = native.dissipated_power_w * resolved.moving_duration_s;
    let force = axes[0]
        .scale(resolved.initial_wrench.force_n[0])
        .add(axes[1].scale(resolved.initial_wrench.force_n[1]));
    let couple = axes[2].scale(resolved.initial_wrench.spin_moment_nm);
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
        let orientation = before[i].pose().orientation();
        body_linear[i] = dp.sub(expected);
        body_angular[i] = dl.sub(
            arms[i]
                .cross(orientation.rotate_world_to_body(expected))
                .add(orientation.rotate_world_to_body(angular_impulse.scale(sign))),
        );
        linear = linear.add(dp);
        angular = angular
            .add(before[i].pose().position_world().cross(dp))
            .add(orientation.rotate_body_to_world(dl));
    }
    // Impulse paired with endpoint mean rates is the exact native kinetic
    // change for frozen symmetric mobility. The complete requested dt is not
    // multiplied by a moving-phase power after the stop.
    let work = impulse.dot(midpoint) + angular_impulse.dot(angular_midpoint);
    let audit = RigidAdhesiveTrajectoryAudit {
        body_linear_momentum_residual_ns: body_linear,
        body_angular_momentum_residual_nms: body_angular,
        linear_momentum_change_ns: linear,
        angular_momentum_change_nms: angular,
        endpoint_velocity_residual_m_s: (endpoint.dot(axes[0]) - resolved.velocity_after_m_s[0])
            .hypot(endpoint.dot(axes[1]) - resolved.velocity_after_m_s[1]),
        endpoint_spin_residual_rad_s: (angular_endpoint.dot(axes[2])
            - resolved.spin_rate_after_rad_s)
            .abs(),
        normal_velocity_residual_m_s: initial.dot(axes[2]).abs().max(endpoint.dot(axes[2]).abs()),
        force_residual_n: (native.force_n[0] - resolved.initial_wrench.force_n[0])
            .hypot(native.force_n[1] - resolved.initial_wrench.force_n[1]),
        spin_moment_residual_nm: (native.spin_moment_nm - resolved.initial_wrench.spin_moment_nm)
            .abs(),
        work_energy_residual_j: kinetic - work,
        total_energy_residual_j: kinetic + resolved.dissipated_energy_j,
        heat_residual_j: resolved.dissipated_energy_j - native_heat,
        maximum_assisting_work_j: [
            force.dot(initial),
            force.dot(endpoint),
            couple.dot(angular_initial),
            couple.dot(angular_endpoint),
        ]
        .into_iter()
        .fold(0.0, f64::max)
            * resolved.moving_duration_s,
    };
    if [
        (
            audit.endpoint_velocity_residual_m_s,
            controls.velocity_tolerance_m_s,
        ),
        (
            audit.endpoint_spin_residual_rad_s,
            controls.angular_velocity_tolerance_rad_s,
        ),
        (
            audit.normal_velocity_residual_m_s,
            controls.velocity_tolerance_m_s,
        ),
        (audit.force_residual_n, controls.contact.force_tolerance_n),
        (
            audit.spin_moment_residual_nm,
            controls.contact.moment_tolerance_nm,
        ),
    ]
    .into_iter()
    .any(|(value, tolerance)| !value.is_finite() || value > tolerance)
        || [
            audit.work_energy_residual_j,
            audit.total_energy_residual_j,
            audit.heat_residual_j,
            audit.maximum_assisting_work_j,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v.abs() > controls.energy_tolerance_j)
        || body_linear
            .into_iter()
            .chain([linear])
            .any(|v| !v.is_finite() || magnitude(v) > controls.linear_momentum_tolerance_ns)
        || body_angular
            .into_iter()
            .chain([angular])
            .any(|v| !v.is_finite() || magnitude(v) > controls.angular_momentum_tolerance_nms)
    {
        return Err(RigidAdhesiveTrajectoryError::Audit(audit));
    }
    let report = RigidAdhesiveTrajectoryReport {
        relaxation: Some(resolved),
        scaled_mobility_per_kg: Some(mobility),
        slip_velocity_before_m_s: velocity,
        slip_velocity_after_m_s: [endpoint.dot(axes[0]), endpoint.dot(axes[1])],
        spin_before_rad_s: spin,
        spin_after_rad_s: angular_endpoint.dot(axes[2]),
        normal_velocity_before_m_s: initial.dot(axes[2]),
        normal_velocity_after_m_s: endpoint.dot(axes[2]),
        force_impulse_on_a_world_ns: impulse,
        spin_angular_impulse_on_a_world_nms: angular_impulse,
        kinetic_energy_change_j: kinetic,
        impulse_work_j: work,
        dissipated_energy_j: resolved.dissipated_energy_j,
        audit,
    };
    poll(gate)?;
    Ok(RigidAdhesiveTrajectoryCandidate {
        bodies_after: after,
        report,
    })
}
