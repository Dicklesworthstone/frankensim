//! Finite-duration Hertz slip/spin with actual arrest and finite-body recoil.
//!
//! The original Hertz law integrates an invariant velocity direction through
//! its stopping event using the complete native mobility. Only the moving
//! interval contributes impulse, travel and heat. Both bodies then receive
//! those integrated impulses once, without assigning or clipping velocities.
//! Actual native endpoints, momentum and work independently check the result.
//!
//! Poses, contact arms, normal load, radius and inertia stay fixed. This is an
//! unforced contact-momentum trajectory; normal/pose evolution and curved
//! velocity paths remain with their surrounding owners.

use super::*;
use fs_tribo::sliding_patch::contact::relaxation::{
    HertzRelaxationError, HertzRelaxationInput, HertzRelaxationReport,
};

#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzTrajectoryReport {
    /// Original physical trajectory; its dwell generates neither impulse nor heat.
    pub relaxation: HertzRelaxationReport,
    pub scaled_mobility_per_kg: [[f64; 3]; 3],
    pub slip_velocity_before_m_s: [f64; 2],
    pub slip_velocity_after_m_s: [f64; 2],
    pub spin_before_rad_s: f64,
    pub spin_after_rad_s: f64,
    /// Exposed to the original normal owner; no normal constraint is imposed here.
    pub normal_velocity_before_m_s: f64,
    pub normal_velocity_after_m_s: f64,
    pub force_impulse_world_ns: Vec3,
    pub spin_angular_impulse_world_nms: Vec3,
    /// Actual momentum increments contracted with native mean rates [J].
    pub kinetic_change_j: f64,
    /// Actual point-force and free-couple impulse work [J].
    pub midpoint_impulse_work_j: f64,
    pub work_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub contact_work_residual_j: f64,
    pub endpoint_velocity_residual_m_s: f64,
    pub endpoint_spin_residual_rad_s: f64,
    /// Moving-phase native mean rate times its duration, minus original travel.
    pub scaled_displacement_residual_m: [f64; 3],
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    /// Positive native endpoint power times moving duration; dwell is excluded.
    pub maximum_assisting_work_j: f64,
}

/// Staged original pair in binding A/B order, for an enclosing atomic owner.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzTrajectoryCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub report: RigidHertzTrajectoryReport,
}

/// Original finite pair and full mobility captured once for coefficient trials.
/// Every advance starts at these same momenta. Normal load, radius, geometry,
/// duration and inertia stay fixed; a coupled source may refresh only friction
/// coefficients. No midpoint-duration compliance is constructed.
#[derive(Clone, Debug)]
pub struct PreparedRigidHertzTrajectory {
    pub(super) before: [RigidBodyState; 2],
    pub(super) properties: [MassProperties; 2],
    pub(super) arms: [Vec3; 2],
    pub(super) axes: [Vec3; 3],
    pub(super) points_before: [PointKinematics; 2],
    pub(super) mobility: [[f64; 3]; 3],
    pub(super) radius: f64,
    pub(super) dt_s: f64,
    pub(super) controls: RigidHertzControls,
    normal_force_n: f64,
}

#[derive(Debug)]
pub enum RigidHertzTrajectoryError {
    Mechanics(RigidHertzError),
    Relaxation(HertzRelaxationError),
    Audit(Box<RigidHertzTrajectoryReport>),
}
impl From<RigidHertzError> for RigidHertzTrajectoryError {
    fn from(error: RigidHertzError) -> Self {
        Self::Mechanics(error)
    }
}
impl From<DynamicsError> for RigidHertzTrajectoryError {
    fn from(error: DynamicsError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl From<RigidFrictionError> for RigidHertzTrajectoryError {
    fn from(error: RigidFrictionError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl From<HertzWrenchError> for RigidHertzTrajectoryError {
    fn from(error: HertzWrenchError) -> Self {
        Self::Relaxation(error.into())
    }
}
impl From<HertzRelaxationError> for RigidHertzTrajectoryError {
    fn from(error: HertzRelaxationError) -> Self {
        Self::Relaxation(error)
    }
}
impl core::fmt::Display for RigidHertzTrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid Hertz trajectory: {self:?}")
    }
}
impl std::error::Error for RigidHertzTrajectoryError {}

/// Advance the complete requested interval, including time after actual arrest.
/// All original native body changes publish together after independent checks.
pub fn apply_rigid_hertz_trajectory(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzContact<'_>,
    dt_s: f64,
    controls: RigidHertzControls,
    gate: &CancelGate,
) -> Result<RigidHertzTrajectoryReport, RigidHertzTrajectoryError> {
    let candidate =
        prepare_rigid_hertz_trajectory(bodies, properties, contact, dt_s, controls, gate)?;
    poll(gate)?;
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    Ok(candidate.report)
}

/// Stage exact fixed-direction Hertz relaxation on the unchanged original pair.
///
/// One original kinetic traction evaluation supplies the entire moving phase.
/// Homogeneity and the full-matrix direction check determine the physical path;
/// actual native impulse work checks its integrated heat independently. A source
/// owner can retry coefficient trials from the same accepted bodies. Geometry,
/// laws and body momenta remain untouched on any refusal or cancellation.
pub fn prepare_rigid_hertz_trajectory(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzContact<'_>,
    dt_s: f64,
    controls: RigidHertzControls,
    gate: &CancelGate,
) -> Result<RigidHertzTrajectoryCandidate, RigidHertzTrajectoryError> {
    PreparedRigidHertzTrajectory::new(bodies, properties, contact, dt_s, controls, gate)?
        .advance(contact.law, gate)
}

impl PreparedRigidHertzTrajectory {
    /// Capture the actual original pair without changing any accepted state.
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contact: RigidHertzContact<'_>,
        dt_s: f64,
        controls: RigidHertzControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidHertzError> {
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
            .any(|value| !value.is_finite() || value < 0.0)
        {
            return Err(RigidHertzError::Invalid(
                "matching bounded native bodies, distinct pair, positive duration and finite gates required",
            ));
        }
        let before = [bodies[binding.body_a], bodies[binding.body_b]];
        let pair_properties = [properties[binding.body_a], properties[binding.body_b]];
        let arms = [binding.arm_a_body, binding.arm_b_body];
        let axes = [
            vector(contact.frame.longitudinal()),
            vector(contact.frame.lateral()),
            vector(contact.frame.normal()),
        ];
        let old = points(before, pair_properties, arms)?;
        if magnitude(old[0].point_world.sub(old[1].point_world)) > controls.geometry_tolerance_m {
            return Err(RigidHertzError::Invalid(
                "original contact points do not coincide",
            ));
        }
        let radius = contact.law.parameters().contact_radius_m;
        let mobility = response(before, pair_properties, arms, axes, radius)?;
        poll(gate)?;
        Ok(Self {
            before,
            properties: pair_properties,
            arms,
            axes,
            points_before: old,
            mobility,
            radius,
            dt_s,
            controls,
            normal_force_n: contact.law.parameters().normal_force_n,
        })
    }

    /// Captured A/B material contact points in the original world frame.
    #[must_use]
    pub fn contact_points_world(&self) -> [Vec3; 2] {
        self.points_before.map(|point| point.point_world)
    }

    pub(super) fn admit_law(&self, law: &HertzWrenchContact) -> Result<(), RigidHertzError> {
        let parameters = law.parameters();
        if parameters.normal_force_n != self.normal_force_n
            || parameters.contact_radius_m != self.radius
        {
            return Err(RigidHertzError::Invalid(
                "trajectory coefficient trials must preserve original normal load and radius",
            ));
        }
        Ok(())
    }

    /// Resolve one exact-ray coefficient trial from the unchanged original pair.
    #[allow(clippy::too_many_lines)]
    pub fn advance(
        &self,
        law: &HertzWrenchContact,
        gate: &CancelGate,
    ) -> Result<RigidHertzTrajectoryCandidate, RigidHertzTrajectoryError> {
        poll(gate)?;
        self.admit_law(law)?;
        let before = self.before;
        let pair_properties = self.properties;
        let arms = self.arms;
        let axes = self.axes;
        let old = self.points_before;
        let radius = self.radius;
        let mobility = self.mobility;
        let dt_s = self.dt_s;
        let controls = self.controls;
        let initial = relative_point(old)?;
        let angular_initial = relative_spin(old);
        let velocity = [initial.dot(axes[0]), initial.dot(axes[1])];
        let spin = angular_initial.dot(axes[2]);
        let resolved = law.relax_velocity(
            &HertzRelaxationInput {
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
        let native = prepare_pair_wrench(
            before,
            pair_properties,
            arms,
            impulse,
            angular_impulse,
            gate,
        )?;
        let new = points(native.bodies_after, pair_properties, arms)?;
        let endpoint = relative_point(new)?;
        let angular_endpoint = relative_spin(new);
        let midpoint = initial.scale(0.5).add(endpoint.scale(0.5));
        let angular_midpoint = angular_initial.scale(0.5).add(angular_endpoint.scale(0.5));
        let work = impulse.dot(midpoint) + angular_impulse.dot(angular_midpoint);
        let moving = resolved.moving_duration_s;
        let force = axes[0]
            .scale(resolved.initial_wrench.force_n[0])
            .add(axes[1].scale(resolved.initial_wrench.force_n[1]));
        let couple = axes[2].scale(resolved.initial_wrench.spin_moment_nm);
        let report = RigidHertzTrajectoryReport {
            relaxation: resolved,
            scaled_mobility_per_kg: mobility,
            slip_velocity_before_m_s: velocity,
            slip_velocity_after_m_s: [endpoint.dot(axes[0]), endpoint.dot(axes[1])],
            spin_before_rad_s: spin,
            spin_after_rad_s: angular_endpoint.dot(axes[2]),
            normal_velocity_before_m_s: initial.dot(axes[2]),
            normal_velocity_after_m_s: endpoint.dot(axes[2]),
            force_impulse_world_ns: impulse,
            spin_angular_impulse_world_nms: angular_impulse,
            kinetic_change_j: native.kinetic_change_j,
            midpoint_impulse_work_j: work,
            work_energy_residual_j: native.kinetic_change_j - work,
            total_energy_residual_j: native.kinetic_change_j + resolved.dissipated_energy_j,
            contact_work_residual_j: work + resolved.contact_work_j,
            endpoint_velocity_residual_m_s: (endpoint.dot(axes[0])
                - resolved.velocity_after_m_s[0])
                .hypot(endpoint.dot(axes[1]) - resolved.velocity_after_m_s[1]),
            endpoint_spin_residual_rad_s: (angular_endpoint.dot(axes[2])
                - resolved.spin_rate_after_rad_s)
                .abs(),
            scaled_displacement_residual_m: [
                moving * midpoint.dot(axes[0]) - resolved.slip_increment_m[0],
                moving * midpoint.dot(axes[1]) - resolved.slip_increment_m[1],
                radius * (moving * angular_midpoint.dot(axes[2]) - resolved.spin_increment_rad),
            ],
            linear_momentum_change_ns: native.linear_change_ns,
            angular_momentum_change_nms: native.angular_change_nms,
            maximum_assisting_work_j: [
                force.dot(initial),
                force.dot(endpoint),
                couple.dot(angular_initial),
                couple.dot(angular_endpoint),
            ]
            .into_iter()
            .fold(0.0, f64::max)
                * moving,
        };
        if [
            (
                report.endpoint_velocity_residual_m_s,
                controls.velocity_tolerance_m_s,
            ),
            (
                report.endpoint_spin_residual_rad_s,
                controls.angular_velocity_tolerance_rad_s,
            ),
            (
                magnitude(report.linear_momentum_change_ns),
                controls.linear_momentum_tolerance_ns,
            ),
            (
                magnitude(report.angular_momentum_change_nms),
                controls.angular_momentum_tolerance_nms,
            ),
        ]
        .into_iter()
        .any(|(value, tolerance)| !value.is_finite() || value > tolerance)
            || [
                report.work_energy_residual_j,
                report.total_energy_residual_j,
                report.contact_work_residual_j,
                report.maximum_assisting_work_j,
            ]
            .into_iter()
            .any(|value| !value.is_finite() || value.abs() > controls.energy_tolerance_j)
            || report
                .scaled_displacement_residual_m
                .into_iter()
                .any(|value| {
                    !value.is_finite() || value.abs() > controls.contact.displacement_tolerance_m
                })
        {
            return Err(RigidHertzTrajectoryError::Audit(Box::new(report)));
        }
        poll(gate)?;
        Ok(RigidHertzTrajectoryCandidate {
            bodies_after: native.bodies_after,
            report,
        })
    }
}
