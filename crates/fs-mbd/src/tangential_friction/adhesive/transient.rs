//! Original finite-body recoil from curved uniform-junction trajectories.
//!
//! The original source, true contact radius and actual complete body mobility
//! feed the bounded transient law. Both bodies receive only its integrated
//! point/free-couple impulses. Native endpoint kinetic work is independent of
//! the constitutive physical-power quadrature and analytic extinction tail.
//! There is no average-wrench endpoint friction test: a curved path dissipates
//! through its instantaneous original traction, not a fictitious held force.
//!
//! This remains a frozen-geometry, no-drive impulse phase. A conservative bound
//! on the actual normal-response row keeps normal motion closed throughout the
//! entire path; unsupported normal/shear coupling refuses without a corrective
//! normal impulse. Poses, normal evolution and changing sources have other owners.

use super::*;
use crate::PointKinematics;
use fs_contact::interface_binding::adhesive_shear::spinning::temperature::SourceThermalAdhesivePatch;
use fs_tribo::adhesive::spinning::relaxation::AdhesiveRelaxationInput;
use fs_tribo::adhesive::spinning::transient::{AdhesiveTransientError, AdhesiveTransientReport};

/// Original numerical trajectory accuracy and work budgets.
pub use fs_tribo::adhesive::spinning::transient::AdhesiveTransientControls;

/// Actual original-body endpoint and full trajectory work checks.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidAdhesiveTransientAudit {
    /// Applied minus integrated requested impulse for each body [N s].
    pub body_linear_momentum_residual_ns: [Vec3; 2],
    /// Applied minus requested arm/free-couple impulse for each body [N m s].
    pub body_angular_momentum_residual_nms: [Vec3; 2],
    /// Pair linear momentum closure [N s].
    pub linear_momentum_change_ns: Vec3,
    /// Orbital plus spin angular momentum closure [N m s].
    pub angular_momentum_change_nms: Vec3,
    /// Actual centre slip minus transient target [m/s].
    pub endpoint_velocity_residual_m_s: f64,
    /// Actual relative spin minus transient target [rad/s].
    pub endpoint_spin_residual_rad_s: f64,
    /// Conservative absolute normal-rate bound over the entire moving path [m/s].
    pub normal_velocity_bound_m_s: f64,
    /// Native body kinetic change minus actual endpoint-paired impulse work [J].
    pub work_energy_residual_j: f64,
    /// Native body kinetic change plus original physical trajectory loss [J].
    pub total_energy_residual_j: f64,
}

/// Source-owned trajectory and actual finite support recoil.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveTransientReport {
    /// Original bounded transient, absent for detached zero-area sources.
    pub transient: Option<AdhesiveTransientReport>,
    /// Actual complete radius-scaled body response, absent after detachment.
    pub scaled_mobility_per_kg: Option<[[f64; 3]; 3]>,
    /// Actual initial and final relative centre velocities [m/s].
    pub slip_velocity_m_s: [[f64; 2]; 2],
    /// Actual initial and final relative normal-axis spin [rad/s].
    pub spin_rad_s: [f64; 2],
    /// Equal/opposite point impulse acting on the original side A [N s].
    pub force_impulse_on_a_world_ns: Vec3,
    /// Equal/opposite free angular impulse acting on original side A [N m s].
    pub spin_angular_impulse_on_a_world_nms: Vec3,
    /// Actual body kinetic change without subtracting large absolute totals [J].
    pub kinetic_energy_change_j: f64,
    /// Actual original-body work of the complete integrated impulse [J].
    pub impulse_work_j: f64,
    /// Original path quadrature plus physical extinction-tail heat [J].
    pub dissipated_energy_j: f64,
    /// Independent native checks; absent for a detached exact no-op.
    pub audit: Option<RigidAdhesiveTransientAudit>,
}

/// Private staged endpoints; callers can compose their own atomic publication.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveTransientCandidate {
    /// Original binding's A/B order; poses are unchanged.
    pub bodies_after: [RigidBodyState; 2],
    /// Native and original-law evidence from the same candidate.
    pub report: RigidAdhesiveTransientReport,
}

/// Original source, frozen finite bodies and complete admitted contact mobility.
///
/// Construction captures the actual pair and computes its geometry, point
/// kinematics, radius-scaled mobility and normal-response row once. Every
/// advance starts from these same endpoints; no trial accumulates a previous
/// trial's impulses. This is a reusable no-drive impulse phase, not a mutable
/// accepted state or a thermal/normal-motion integrator.
#[derive(Clone, Debug)]
pub struct PreparedRigidAdhesiveTransient<'a> {
    contact: RigidAdhesiveContact<'a>,
    before: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    arms: [Vec3; 2],
    axes: [Vec3; 3],
    points_before: [PointKinematics; 2],
    attached: Option<AttachedTransient>,
    dt_s: f64,
    controls: RigidAdhesiveControls,
    transient: AdhesiveTransientControls,
}

#[derive(Clone, Debug)]
struct AttachedTransient {
    law: AdhesiveWrenchContact,
    mobility: [[f64; 3]; 3],
    normal_response_row: Vec3,
}

/// Source trajectory, native mechanics and whole-pair refusals remain distinct.
#[derive(Debug)]
pub enum RigidAdhesiveTransientError {
    /// Invalid original binding, controls or held-normal geometry.
    Invalid(&'static str),
    /// Original finite-body mechanics or cancellation refusal.
    Mechanics(RigidHertzError),
    /// Original bounded traction integration refusal.
    Transient(AdhesiveTransientError),
    /// Actual native body candidate failed its dimensional gates.
    Audit(RigidAdhesiveTransientAudit),
}
impl From<RigidHertzError> for RigidAdhesiveTransientError {
    fn from(e: RigidHertzError) -> Self {
        Self::Mechanics(e)
    }
}
impl From<DynamicsError> for RigidAdhesiveTransientError {
    fn from(e: DynamicsError) -> Self {
        Self::Mechanics(e.into())
    }
}
impl From<super::super::RigidFrictionError> for RigidAdhesiveTransientError {
    fn from(e: super::super::RigidFrictionError) -> Self {
        Self::Mechanics(e.into())
    }
}
impl From<AdhesiveWrenchError> for RigidAdhesiveTransientError {
    fn from(e: AdhesiveWrenchError) -> Self {
        Self::Transient(e.into())
    }
}
impl From<AdhesiveTransientError> for RigidAdhesiveTransientError {
    fn from(e: AdhesiveTransientError) -> Self {
        Self::Transient(e)
    }
}
impl core::fmt::Display for RigidAdhesiveTransientError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid adhesive transient: {self:?}")
    }
}
impl std::error::Error for RigidAdhesiveTransientError {}

/// Publish both original momenta only after the entire curved path is admitted.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_adhesive_transient(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidAdhesiveContact<'_>,
    dt_s: f64,
    controls: RigidAdhesiveControls,
    transient: AdhesiveTransientControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveTransientReport, RigidAdhesiveTransientError> {
    let candidate = prepare_rigid_adhesive_transient(
        bodies, properties, contact, dt_s, controls, transient, gate,
    )?;
    poll(gate)?;
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    Ok(candidate.report)
}

/// Integrate the original source law against actual finite-body mechanics.
#[allow(clippy::too_many_arguments)]
pub fn prepare_rigid_adhesive_transient(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidAdhesiveContact<'_>,
    dt_s: f64,
    controls: RigidAdhesiveControls,
    transient: AdhesiveTransientControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveTransientCandidate, RigidAdhesiveTransientError> {
    PreparedRigidAdhesiveTransient::new(
        bodies, properties, contact, dt_s, controls, transient, gate,
    )?
    .advance(gate)
}

impl<'a> PreparedRigidAdhesiveTransient<'a> {
    /// Admit the original source-owned disk and capture its actual body response.
    /// Detached sources retain the original exact no-op behavior and do not
    /// acquire contact mobility, a solid heat face or geometric constraints.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contact: RigidAdhesiveContact<'a>,
        dt_s: f64,
        controls: RigidAdhesiveControls,
        transient: AdhesiveTransientControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidAdhesiveTransientError> {
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
            .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(RigidAdhesiveTransientError::Invalid(
                "bounded original bodies, positive time and finite gates required",
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
        let radius = contact.source.patch().parameters().contact_radius_m;
        let attached = if radius == 0.0 {
            None
        } else {
            if magnitude(old[0].point_world.sub(old[1].point_world)) > controls.geometry_tolerance_m
                || initial.dot(axes[2]).abs() > controls.velocity_tolerance_m_s
            {
                return Err(RigidAdhesiveTransientError::Invalid(
                    "attached points must coincide with held normal motion",
                ));
            }
            Some(AttachedTransient {
                law: AdhesiveWrenchContact::new(AdhesiveWrenchParameters {
                    contact_radius_m: radius,
                    static_shear_strength_pa: contact
                        .source
                        .static_patch()
                        .parameters()
                        .shear_strength_pa,
                    kinetic_shear_strength_pa: contact
                        .source
                        .patch()
                        .parameters()
                        .shear_strength_pa,
                })?,
                mobility: scaled_response(before, properties, arms, axes, [0.0, 0.0, radius])?,
                normal_response_row: normal_response_row(before, properties, arms, axes, radius)?,
            })
        };
        poll(gate)?;
        Ok(Self {
            contact,
            before,
            properties,
            arms,
            axes,
            points_before: old,
            attached,
            dt_s,
            controls,
            transient,
        })
    }

    /// Captured world positions of the original A/B material contact points.
    /// Attached pairs have already passed the original coincidence gate.
    #[must_use]
    pub fn contact_points_world(&self) -> [Vec3; 2] {
        self.points_before.map(|point| point.point_world)
    }

    /// Integrate the original source law from the unchanged captured endpoints.
    pub fn advance(
        &self,
        gate: &CancelGate,
    ) -> Result<RigidAdhesiveTransientCandidate, RigidAdhesiveTransientError> {
        poll(gate)?;
        if let Some(attached) = &self.attached {
            return self.advance_law(attached, &attached.law, gate);
        }
        let initial = relative_point(self.points_before);
        let velocity = [initial.dot(self.axes[0]), initial.dot(self.axes[1])];
        let spin = relative_spin(self.points_before).dot(self.axes[2]);
        poll(gate)?;
        Ok(RigidAdhesiveTransientCandidate {
            bodies_after: self.before,
            report: RigidAdhesiveTransientReport {
                transient: None,
                scaled_mobility_per_kg: None,
                slip_velocity_m_s: [velocity; 2],
                spin_rad_s: [spin; 2],
                force_impulse_on_a_world_ns: Vec3::ZERO,
                spin_angular_impulse_on_a_world_nms: Vec3::ZERO,
                kinetic_energy_change_j: 0.0,
                impulse_work_j: 0.0,
                dissipated_energy_j: 0.0,
                audit: None,
            },
        })
    }

    /// Integrate the original selected junction claims at a proposed face
    /// temperature using the same captured geometry and full finite-body recoil.
    /// The COMPLETE original patch must match, including both source bundles,
    /// original strengths, radius and area. Detached faces cannot enter this
    /// thermal operation. This changes neither the normal law nor body poses.
    /// The enclosing owner must close the temperature and finite heat balance.
    pub fn advance_thermal(
        &self,
        source: &SourceThermalAdhesivePatch<'_>,
        gate: &CancelGate,
    ) -> Result<RigidAdhesiveTransientCandidate, RigidAdhesiveTransientError> {
        poll(gate)?;
        if source.original() != self.contact.source {
            return Err(RigidAdhesiveTransientError::Invalid(
                "thermal transient requires its complete original source-owned patch",
            ));
        }
        let attached = self
            .attached
            .as_ref()
            .ok_or(RigidAdhesiveTransientError::Invalid(
                "thermal transient requires an attached original solid disk",
            ))?;
        self.advance_law(attached, source.law(), gate)
    }

    #[allow(clippy::too_many_lines)]
    fn advance_law(
        &self,
        attached: &AttachedTransient,
        law: &AdhesiveWrenchContact,
        gate: &CancelGate,
    ) -> Result<RigidAdhesiveTransientCandidate, RigidAdhesiveTransientError> {
        poll(gate)?;
        let before = self.before;
        let properties = self.properties;
        let arms = self.arms;
        let axes = self.axes;
        let old = self.points_before;
        let initial = relative_point(old);
        let initial_spin = relative_spin(old);
        let velocity = [initial.dot(axes[0]), initial.dot(axes[1])];
        let spin = initial_spin.dot(axes[2]);
        let mobility = attached.mobility;
        let dt_s = self.dt_s;
        let controls = self.controls;
        let path = law.advance_transient(
            &AdhesiveRelaxationInput {
                velocity_m_s: velocity,
                spin_rate_rad_s: spin,
                scaled_mobility_per_kg: mobility,
                dt_s,
            },
            controls.contact,
            self.transient,
            || gate.is_requested(),
        )?;
        let moving_upper = path.stopping_time_bracket_s.map_or(dt_s, |t| t[1]);
        let normal_bound = initial.dot(axes[2]).abs()
            + magnitude(attached.normal_response_row)
                * law.kinetic_patch().sliding_capacity_n()
                * (13.0_f64.sqrt() / 3.0)
                * moving_upper;
        if !normal_bound.is_finite() || normal_bound > controls.velocity_tolerance_m_s {
            return Err(RigidAdhesiveTransientError::Invalid(
                "curved traction cannot preserve the held normal-motion bound",
            ));
        }
        let impulse = axes[0]
            .scale(path.force_impulse_ns[0])
            .add(axes[1].scale(path.force_impulse_ns[1]));
        let couple = axes[2].scale(path.spin_angular_impulse_nms);
        let after = apply_pair_wrench(before, properties, arms, impulse, couple, gate)?;
        let new = points(after, properties, arms)?;
        let endpoint = relative_point(new);
        let endpoint_spin = relative_spin(new);
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
            let orientation = before[i].pose().orientation();
            let expected = impulse.scale(sign);
            body_linear[i] = dp.sub(expected);
            body_angular[i] = dl.sub(
                arms[i]
                    .cross(orientation.rotate_world_to_body(expected))
                    .add(orientation.rotate_world_to_body(couple.scale(sign))),
            );
            linear = linear.add(dp);
            angular = angular
                .add(before[i].pose().position_world().cross(dp))
                .add(orientation.rotate_body_to_world(dl));
        }
        let work = impulse.dot(initial.scale(0.5).add(endpoint.scale(0.5)))
            + couple.dot(initial_spin.scale(0.5).add(endpoint_spin.scale(0.5)));
        let audit = RigidAdhesiveTransientAudit {
            body_linear_momentum_residual_ns: body_linear,
            body_angular_momentum_residual_nms: body_angular,
            linear_momentum_change_ns: linear,
            angular_momentum_change_nms: angular,
            endpoint_velocity_residual_m_s: (endpoint.dot(axes[0]) - path.velocity_after_m_s[0])
                .hypot(endpoint.dot(axes[1]) - path.velocity_after_m_s[1]),
            endpoint_spin_residual_rad_s: (endpoint_spin.dot(axes[2]) - path.spin_rate_after_rad_s)
                .abs(),
            normal_velocity_bound_m_s: normal_bound.max(endpoint.dot(axes[2]).abs()),
            work_energy_residual_j: kinetic - work,
            total_energy_residual_j: kinetic + path.dissipated_energy_j,
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
                audit.normal_velocity_bound_m_s,
                controls.velocity_tolerance_m_s,
            ),
        ]
        .into_iter()
        .any(|(v, t)| !v.is_finite() || v > t)
            || [
                audit.work_energy_residual_j,
                audit.total_energy_residual_j,
                work + path.contact_work_j,
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
            return Err(RigidAdhesiveTransientError::Audit(audit));
        }
        let report = RigidAdhesiveTransientReport {
            scaled_mobility_per_kg: Some(mobility),
            slip_velocity_m_s: [velocity, [endpoint.dot(axes[0]), endpoint.dot(axes[1])]],
            spin_rad_s: [spin, endpoint_spin.dot(axes[2])],
            force_impulse_on_a_world_ns: impulse,
            spin_angular_impulse_on_a_world_nms: couple,
            kinetic_energy_change_j: kinetic,
            impulse_work_j: work,
            dissipated_energy_j: path.dissipated_energy_j,
            audit: Some(audit),
            transient: Some(path),
        };
        poll(gate)?;
        Ok(RigidAdhesiveTransientCandidate {
            bodies_after: after,
            report,
        })
    }
}

/// Actual normal point-rate response to each of [Fx,Fy,M/a], including all arms
/// and inertia terms. Its norm bounds every possible time-dependent disk wrench.
fn normal_response_row(
    before: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    arms: [Vec3; 2],
    axes: [Vec3; 3],
    radius: f64,
) -> Result<Vec3, RigidAdhesiveTransientError> {
    let mut row = [0.0; 3];
    for i in 0..2 {
        let orientation = before[i].pose().orientation();
        for k in 0..3 {
            let force = if k < 2 { axes[k] } else { Vec3::ZERO };
            let couple = if k == 2 {
                axes[2].scale(radius)
            } else {
                Vec3::ZERO
            };
            let torque = arms[i]
                .cross(orientation.rotate_world_to_body(force))
                .add(orientation.rotate_world_to_body(couple));
            let omega = orientation
                .rotate_body_to_world(properties[i].angular_velocity_body_checked(torque)?);
            let rate = force
                .scale(properties[i].mass().recip())
                .add(omega.cross(orientation.rotate_body_to_world(arms[i])));
            row[k] += rate.dot(axes[2]);
        }
    }
    if row.into_iter().any(|v| !v.is_finite()) {
        return Err(RigidAdhesiveTransientError::Invalid(
            "unrepresentable normal response",
        ));
    }
    Ok(Vec3::new(row[0], row[1], row[2]))
}
