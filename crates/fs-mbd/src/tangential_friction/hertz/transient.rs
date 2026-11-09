//! Curved Hertz shear/spin trajectories on the original finite-body pair.
//!
//! The original pressure-resolved law supplies instantaneous traction and
//! physical-power quadrature. Its integrated point/couple impulses update both
//! native momenta once. Actual endpoint-paired impulse work independently checks
//! the heat; numerical work defects are never deposited as frictional heat.
//!
//! Preparation captures geometry and full mobility once for coefficient trials.
//! A conservative normal-response bound admits only a held closed contact over
//! the whole path. Poses, normal-load evolution and external drive remain with
//! their original owners. Tail time/travel uncertainty stays in the law report.

use super::trajectory::PreparedRigidHertzTrajectory;
use super::*;
use fs_tribo::sliding_patch::contact::relaxation::HertzRelaxationInput;
use fs_tribo::sliding_patch::contact::transient::{HertzTransientError, HertzTransientReport};

/// Physical trajectory accuracy and bounded numerical work.
pub use fs_tribo::sliding_patch::contact::transient::HertzTransientControls;

/// Independent checks against the actual native endpoints and original moments.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidHertzTransientAudit {
    pub body_linear_momentum_residual_ns: [Vec3; 2],
    pub body_angular_momentum_residual_nms: [Vec3; 2],
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    pub endpoint_velocity_residual_m_s: f64,
    pub endpoint_spin_residual_rad_s: f64,
    /// Bound over the entire moving path, including any finite-extinction tail.
    pub normal_velocity_bound_m_s: f64,
    pub work_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub contact_work_residual_j: f64,
}

/// Original physical path and finite-body response from one admitted candidate.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzTransientReport {
    pub transient: HertzTransientReport,
    pub scaled_mobility_per_kg: [[f64; 3]; 3],
    pub slip_velocity_m_s: [[f64; 2]; 2],
    pub spin_rad_s: [f64; 2],
    pub normal_velocity_m_s: [f64; 2],
    /// Actual requested integrated resultants on A; B receives their opposites.
    pub force_impulse_world_ns: Vec3,
    pub spin_angular_impulse_world_nms: Vec3,
    /// Native momentum increments contracted with native mean rates [J].
    pub kinetic_change_j: f64,
    /// Integrated impulse paired with the actual native endpoint mean [J].
    /// This is a work identity, not a held-wrench approximation to the path.
    pub impulse_work_j: f64,
    /// Original positive path quadrature plus physical extinction-tail loss [J].
    pub dissipated_energy_j: f64,
    pub audit: RigidHertzTransientAudit,
}

/// Staged original A/B endpoints; the caller owns joint thermal publication.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzTransientCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub report: RigidHertzTransientReport,
}

/// Full original response prepared once, including the normal-coupling row.
#[derive(Clone, Debug)]
pub struct PreparedRigidHertzTransient {
    pub(super) pair: PreparedRigidHertzTrajectory,
    normal_response_row: Vec3,
    pub(super) transient: HertzTransientControls,
}

#[derive(Debug)]
pub enum RigidHertzTransientError {
    Invalid(&'static str),
    Mechanics(RigidHertzError),
    Transient(HertzTransientError),
    Audit(RigidHertzTransientAudit),
}
impl From<RigidHertzError> for RigidHertzTransientError {
    fn from(error: RigidHertzError) -> Self {
        Self::Mechanics(error)
    }
}
impl From<DynamicsError> for RigidHertzTransientError {
    fn from(error: DynamicsError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl From<RigidFrictionError> for RigidHertzTransientError {
    fn from(error: RigidFrictionError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl From<HertzTransientError> for RigidHertzTransientError {
    fn from(error: HertzTransientError) -> Self {
        Self::Transient(error)
    }
}
impl core::fmt::Display for RigidHertzTransientError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid Hertz transient: {self:?}")
    }
}
impl std::error::Error for RigidHertzTransientError {}

/// Publish the pair only after original-law and independent native work checks.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_hertz_transient(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzContact<'_>,
    dt_s: f64,
    controls: RigidHertzControls,
    transient: HertzTransientControls,
    gate: &CancelGate,
) -> Result<RigidHertzTransientReport, RigidHertzTransientError> {
    let candidate = prepare_rigid_hertz_transient(
        bodies, properties, contact, dt_s, controls, transient, gate,
    )?;
    poll(gate)?;
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    Ok(candidate.report)
}

/// Stage curved Hertz relaxation on the unchanged original finite bodies.
#[allow(clippy::too_many_arguments)]
pub fn prepare_rigid_hertz_transient(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzContact<'_>,
    dt_s: f64,
    controls: RigidHertzControls,
    transient: HertzTransientControls,
    gate: &CancelGate,
) -> Result<RigidHertzTransientCandidate, RigidHertzTransientError> {
    PreparedRigidHertzTransient::new(bodies, properties, contact, dt_s, controls, transient, gate)?
        .advance(contact.law, gate)
}

impl PreparedRigidHertzTransient {
    /// Capture both original bodies and the same geometry used by the exact ray.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contact: RigidHertzContact<'_>,
        dt_s: f64,
        controls: RigidHertzControls,
        transient: HertzTransientControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidHertzTransientError> {
        let pair =
            PreparedRigidHertzTrajectory::new(bodies, properties, contact, dt_s, controls, gate)?;
        Self::from_pair(pair, transient, gate)
    }

    /// Shared actual geometry and normal response for circular pressure laws.
    pub(super) fn from_pair(
        pair: PreparedRigidHertzTrajectory,
        transient: HertzTransientControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidHertzTransientError> {
        if relative_point(pair.points_before)?.dot(pair.axes[2]).abs()
            > pair.controls.velocity_tolerance_m_s
        {
            return Err(RigidHertzTransientError::Invalid(
                "curved Hertz trajectory requires held closed normal motion",
            ));
        }
        let normal_response_row = normal_response_row(&pair)?;
        poll(gate)?;
        Ok(Self {
            pair,
            normal_response_row,
            transient,
        })
    }

    /// Original A/B material contact points, already checked for coincidence.
    #[must_use]
    pub fn contact_points_world(&self) -> [Vec3; 2] {
        self.pair.contact_points_world()
    }

    /// Retry coefficients from the same original state, retaining load and radius.
    #[allow(clippy::too_many_lines)]
    pub fn advance(
        &self,
        law: &HertzWrenchContact,
        gate: &CancelGate,
    ) -> Result<RigidHertzTransientCandidate, RigidHertzTransientError> {
        poll(gate)?;
        let pair = &self.pair;
        pair.admit_law(law)?;
        let axes = pair.axes;
        let controls = pair.controls;
        let initial = relative_point(pair.points_before)?;
        let initial_spin = relative_spin(pair.points_before);
        let velocity = [initial.dot(axes[0]), initial.dot(axes[1])];
        let spin = initial_spin.dot(axes[2]);
        let path = law.advance_transient(
            &HertzRelaxationInput {
                velocity_m_s: velocity,
                spin_rate_rad_s: spin,
                scaled_mobility_per_kg: pair.mobility,
                dt_s: pair.dt_s,
            },
            controls.contact,
            self.transient,
            || gate.is_requested(),
        )?;
        let capacity = law.parameters().kinetic_mu * law.parameters().normal_force_n;
        self.finish_path(
            path,
            capacity.hypot((3.0 * core::f64::consts::PI / 16.0) * capacity),
            gate,
        )
    }

    /// Apply an original constitutive path exactly once and audit it against
    /// the original native bodies. The supplied bound covers every actual
    /// radius-scaled wrench, including any terminal tail.
    #[allow(clippy::too_many_lines)]
    pub(super) fn finish_path(
        &self,
        path: HertzTransientReport,
        scaled_wrench_bound_n: f64,
        gate: &CancelGate,
    ) -> Result<RigidHertzTransientCandidate, RigidHertzTransientError> {
        poll(gate)?;
        if !scaled_wrench_bound_n.is_finite() || scaled_wrench_bound_n < 0.0 {
            return Err(RigidHertzTransientError::Invalid(
                "finite nonnegative original traction bound required",
            ));
        }
        let pair = &self.pair;
        let axes = pair.axes;
        let controls = pair.controls;
        let initial = relative_point(pair.points_before)?;
        let initial_spin = relative_spin(pair.points_before);
        let velocity = [initial.dot(axes[0]), initial.dot(axes[1])];
        let spin = initial_spin.dot(axes[2]);
        // |F| <= mu*N and |M/a| <= (3*pi/16)*mu*N for the same Hertz patch.
        // Cauchy--Schwarz bounds every intermediate normal rate, including the
        // unresolved extinction tail; checking endpoints alone is insufficient.
        let moving_upper = path.stopping_time_bracket_s.map_or(pair.dt_s, |t| t[1]);
        let normal_change = if scaled_wrench_bound_n == 0.0 || moving_upper == 0.0 {
            0.0
        } else {
            magnitude(self.normal_response_row) * (moving_upper * scaled_wrench_bound_n)
        };
        let normal_bound = initial.dot(axes[2]).abs() + normal_change;
        if !normal_bound.is_finite() || normal_bound > controls.velocity_tolerance_m_s {
            return Err(RigidHertzTransientError::Invalid(
                "curved traction cannot preserve the held normal-motion bound",
            ));
        }
        let impulse = axes[0]
            .scale(path.force_impulse_ns[0])
            .add(axes[1].scale(path.force_impulse_ns[1]));
        let couple = axes[2].scale(path.spin_angular_impulse_nms);
        let native = prepare_pair_wrench(
            pair.before,
            pair.properties,
            pair.arms,
            impulse,
            couple,
            gate,
        )?;
        let after = native.bodies_after;
        let new = points(after, pair.properties, pair.arms)?;
        let endpoint = relative_point(new)?;
        let endpoint_spin = relative_spin(new);
        let work = impulse.dot(initial.scale(0.5).add(endpoint.scale(0.5)))
            + couple.dot(initial_spin.scale(0.5).add(endpoint_spin.scale(0.5)));
        let mut body_linear = [Vec3::ZERO; 2];
        let mut body_angular = [Vec3::ZERO; 2];
        for (i, sign) in [(0, 1.0), (1, -1.0)] {
            let orientation = pair.before[i].pose().orientation();
            let expected = impulse.scale(sign);
            body_linear[i] = after[i]
                .linear_momentum_world()
                .sub(pair.before[i].linear_momentum_world())
                .sub(expected);
            body_angular[i] = after[i]
                .angular_momentum_body()
                .sub(pair.before[i].angular_momentum_body())
                .sub(
                    pair.arms[i]
                        .cross(orientation.rotate_world_to_body(expected))
                        .add(orientation.rotate_world_to_body(couple.scale(sign))),
                );
        }
        let audit = RigidHertzTransientAudit {
            body_linear_momentum_residual_ns: body_linear,
            body_angular_momentum_residual_nms: body_angular,
            linear_momentum_change_ns: native.linear_change_ns,
            angular_momentum_change_nms: native.angular_change_nms,
            endpoint_velocity_residual_m_s: (endpoint.dot(axes[0]) - path.velocity_after_m_s[0])
                .hypot(endpoint.dot(axes[1]) - path.velocity_after_m_s[1]),
            endpoint_spin_residual_rad_s: (endpoint_spin.dot(axes[2]) - path.spin_rate_after_rad_s)
                .abs(),
            normal_velocity_bound_m_s: normal_bound.max(endpoint.dot(axes[2]).abs()),
            work_energy_residual_j: native.kinetic_change_j - work,
            total_energy_residual_j: native.kinetic_change_j + path.dissipated_energy_j,
            contact_work_residual_j: work + path.contact_work_j,
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
                audit.contact_work_residual_j,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v.abs() > controls.energy_tolerance_j)
            || body_linear
                .into_iter()
                .chain([native.linear_change_ns])
                .any(|v| !v.is_finite() || magnitude(v) > controls.linear_momentum_tolerance_ns)
            || body_angular
                .into_iter()
                .chain([native.angular_change_nms])
                .any(|v| !v.is_finite() || magnitude(v) > controls.angular_momentum_tolerance_nms)
        {
            return Err(RigidHertzTransientError::Audit(audit));
        }
        let report = RigidHertzTransientReport {
            scaled_mobility_per_kg: pair.mobility,
            slip_velocity_m_s: [velocity, [endpoint.dot(axes[0]), endpoint.dot(axes[1])]],
            spin_rad_s: [spin, endpoint_spin.dot(axes[2])],
            normal_velocity_m_s: [initial.dot(axes[2]), endpoint.dot(axes[2])],
            force_impulse_world_ns: impulse,
            spin_angular_impulse_world_nms: couple,
            kinetic_change_j: native.kinetic_change_j,
            impulse_work_j: work,
            dissipated_energy_j: path.dissipated_energy_j,
            audit,
            transient: path,
        };
        poll(gate)?;
        Ok(RigidHertzTransientCandidate {
            bodies_after: after,
            report,
        })
    }
}

/// Actual normal point-rate response to [Fx,Fy,M/a] with original full inertia.
fn normal_response_row(pair: &PreparedRigidHertzTrajectory) -> Result<Vec3, RigidHertzError> {
    let mut row = [0.0; 3];
    for i in 0..2 {
        let orientation = pair.before[i].pose().orientation();
        for k in 0..3 {
            let force = if k < 2 { pair.axes[k] } else { Vec3::ZERO };
            let couple = if k == 2 {
                pair.axes[2].scale(pair.radius)
            } else {
                Vec3::ZERO
            };
            let torque = pair.arms[i]
                .cross(orientation.rotate_world_to_body(force))
                .add(orientation.rotate_world_to_body(couple));
            let omega = orientation
                .rotate_body_to_world(pair.properties[i].angular_velocity_body_checked(torque)?);
            let rate = force
                .scale(pair.properties[i].mass().recip())
                .add(omega.cross(orientation.rotate_body_to_world(pair.arms[i])));
            row[k] = finite(row[k] + rate.dot(pair.axes[2]))?;
        }
    }
    Ok(Vec3::new(row[0], row[1], row[2]))
}
