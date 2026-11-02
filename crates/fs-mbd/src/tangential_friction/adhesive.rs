//! Uniform adhesive junction shear and spin on two actual finite bodies.
//!
//! The original JKR, DMT or Maugis response supplies the real attached disk;
//! its original ordered interface supplies both shear strengths. The same disk
//! supports translation and twisting, including at zero or tensile net load.
//! Body masses, full inertia tensors and contact arms determine all nine recoil
//! entries. The existing circular-support solver stages the joint wrench, and
//! native point/angular velocities independently check its work and loss.
//!
//! This is a held-geometry impulse phase. Normal evolution, poses, mixed-mode
//! peeling, partial-slip elasticity and the destination of frictional heat
//! remain with their existing owners. Endpoint reversals require a refined
//! physical interval; a zero midpoint is never silently turned into arrest.

/// Simultaneous source-owned junctions sharing actual native body recoil.
pub mod multiple;

use fs_contact::interface_binding::adhesive_shear::spinning::SourceAdhesiveSlidingPatch;
use fs_exec::CancelGate;
use fs_tribo::adhesive::spinning::contact::{
    AdhesiveWrenchContact, AdhesiveWrenchError, AdhesiveWrenchInput, AdhesiveWrenchMode,
    AdhesiveWrenchParameters, AdhesiveWrenchStep,
};
use fs_tribo::partial_slip::TangentFrame;

use super::hertz::{
    RigidHertzError, accepted_phase, apply_pair_wrench, points, relative_point, relative_spin,
    scaled_response,
};
use super::{RigidContactBinding, magnitude, poll, totals, vector};
use crate::{DynamicsError, MassProperties, RigidBodyState, Vec3};

/// Common dimensional native-body and circular-wrench acceptance controls.
pub use super::hertz::RigidHertzControls as RigidAdhesiveControls;

/// Actual source-owned junction and its original coincident body-frame points.
#[derive(Clone, Copy, Debug)]
pub struct RigidAdhesiveContact<'a> {
    /// Physical pair and COM-relative material arms.
    pub binding: RigidContactBinding,
    /// Original normal and tangent gauge, held over this impulse phase.
    pub frame: TangentFrame,
    /// Borrowed normal and shear sources plus their actual attached geometry.
    pub source: SourceAdhesiveSlidingPatch<'a>,
}

/// Native recoil and independent work checks for one accepted physical sample.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveReport {
    /// Original coupled constitutive result; absent for a detached zero-area disk.
    pub contact: Option<AdhesiveWrenchStep>,
    /// Complete `J M^-1 J^T` in `[u_x,u_y,a*spin]`, absent after detachment.
    pub scaled_mobility_per_kg: Option<[[f64; 3]; 3]>,
    /// Equal/opposite physical point impulse acting on body A.
    pub force_impulse_world_ns: Vec3,
    /// Equal/opposite normal-axis free angular impulse acting on body A.
    pub spin_angular_impulse_world_nms: Vec3,
    /// Actual initial relative point slip.
    pub slip_velocity_before_m_s: [f64; 2],
    /// Actual final relative point slip.
    pub slip_velocity_after_m_s: [f64; 2],
    /// Actual initial relative normal-axis spin.
    pub spin_before_rad_s: f64,
    /// Actual final relative normal-axis spin.
    pub spin_after_rad_s: f64,
    /// Native pair kinetic energy before this sample.
    pub kinetic_before_j: f64,
    /// Native pair kinetic energy after this sample.
    pub kinetic_after_j: f64,
    /// Signed work of actual point and free-couple impulses.
    pub midpoint_impulse_work_j: f64,
    /// Original uniform-junction loss at actual midpoint motion, unpartitioned.
    pub dissipated_energy_j: f64,
    /// Kinetic-energy change plus the original junction loss.
    pub energy_residual_j: f64,
    /// Actual force minus the original native law at actual body motion.
    pub force_residual_n: f64,
    /// Actual moment minus the original native law at actual body motion.
    pub spin_moment_residual_nm: f64,
    /// Equal/opposite linear-impulse closure.
    pub linear_momentum_change_ns: Vec3,
    /// Total orbital plus spin angular-momentum closure.
    pub angular_momentum_change_nms: Vec3,
}

/// Complete private proposal; no accepted body is mutated by preparation.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveCandidate {
    /// Native endpoints in the original binding's A/B order.
    pub bodies_after: [RigidBodyState; 2],
    /// Physical response and independent work accounting.
    pub report: RigidAdhesiveReport,
}

/// Source, native constitutive and finite-body refusals remain distinguishable.
#[derive(Debug)]
pub enum RigidAdhesiveError {
    /// Invalid physical controls, geometry or nonfinite derived data.
    Invalid(&'static str),
    /// Original circular-junction projection refused the sample.
    Contact(AdhesiveWrenchError),
    /// Shared native kinematics, impulse or momentum owner refused.
    Mechanics(RigidHertzError),
    /// The held wrench would assist actual endpoint slip or spin.
    UnresolvedSlipEvent {
        /// Positive assisting endpoint work [J].
        assisting_work_j: f64,
    },
    /// The actual endpoints did not pass the original dimensional balance gates.
    Audit(Box<RigidAdhesiveReport>),
}
impl From<RigidHertzError> for RigidAdhesiveError {
    fn from(error: RigidHertzError) -> Self {
        Self::Mechanics(error)
    }
}
impl From<DynamicsError> for RigidAdhesiveError {
    fn from(error: DynamicsError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl From<super::RigidFrictionError> for RigidAdhesiveError {
    fn from(error: super::RigidFrictionError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl From<AdhesiveWrenchError> for RigidAdhesiveError {
    fn from(error: AdhesiveWrenchError) -> Self {
        Self::Contact(error)
    }
}
impl core::fmt::Display for RigidAdhesiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid adhesive junction: {self:?}")
    }
}
impl std::error::Error for RigidAdhesiveError {}

/// Publish both original body momenta together only after every native check.
pub fn apply_rigid_adhesive(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidAdhesiveContact<'_>,
    dt_s: f64,
    controls: RigidAdhesiveControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveReport, RigidAdhesiveError> {
    let candidate = prepare_rigid_adhesive(bodies, properties, contact, dt_s, controls, gate)?;
    poll(gate)?;
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    Ok(candidate.report)
}

/// Stage one source-preserving sample for composition with an external normal
/// or thermal owner. Retrying starts from the same accepted physical endpoints.
#[allow(clippy::too_many_lines)]
pub fn prepare_rigid_adhesive(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidAdhesiveContact<'_>,
    dt_s: f64,
    controls: RigidAdhesiveControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveCandidate, RigidAdhesiveError> {
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
        return Err(RigidAdhesiveError::Invalid(
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
    let initial = points(before, properties, arms)?;
    let relative_before = relative_point(initial);
    let angular_before = relative_spin(initial);
    let slip_before = [relative_before.dot(axes[0]), relative_before.dot(axes[1])];
    let spin_before = angular_before.dot(axes[2]);
    let radius = contact.source.patch().parameters().contact_radius_m;
    let (energy_before, linear_before, momentum_before) = totals(&before, &properties, gate)?;
    if radius == 0.0 {
        // An attractive DMT gap has no solid junction, so arbitrary free motion
        // and separated material points cannot produce shear, torque or heat.
        poll(gate)?;
        return Ok(RigidAdhesiveCandidate {
            bodies_after: before,
            report: RigidAdhesiveReport {
                contact: None,
                scaled_mobility_per_kg: None,
                force_impulse_world_ns: Vec3::ZERO,
                spin_angular_impulse_world_nms: Vec3::ZERO,
                slip_velocity_before_m_s: slip_before,
                slip_velocity_after_m_s: slip_before,
                spin_before_rad_s: spin_before,
                spin_after_rad_s: spin_before,
                kinetic_before_j: energy_before,
                kinetic_after_j: energy_before,
                midpoint_impulse_work_j: 0.0,
                dissipated_energy_j: 0.0,
                energy_residual_j: 0.0,
                force_residual_n: 0.0,
                spin_moment_residual_nm: 0.0,
                linear_momentum_change_ns: Vec3::ZERO,
                angular_momentum_change_nms: Vec3::ZERO,
            },
        });
    }
    if magnitude(initial[0].point_world.sub(initial[1].point_world)) > controls.geometry_tolerance_m
    {
        return Err(RigidAdhesiveError::Invalid(
            "attached material points must coincide",
        ));
    }
    let law = AdhesiveWrenchContact::new(AdhesiveWrenchParameters {
        contact_radius_m: radius,
        static_shear_strength_pa: contact.source.static_patch().parameters().shear_strength_pa,
        kinetic_shear_strength_pa: contact.source.patch().parameters().shear_strength_pa,
    })?;
    let mobility = scaled_response(before, properties, arms, axes, [0.0, 0.0, radius])?;
    let half_dt_squared = 0.5 * dt_s * dt_s;
    if !half_dt_squared.is_finite() || half_dt_squared == 0.0 {
        return Err(RigidAdhesiveError::Invalid(
            "unrepresentable midpoint duration",
        ));
    }
    let input = AdhesiveWrenchInput {
        free_displacement_m: slip_before.map(|v| dt_s * v),
        free_spin_angle_rad: dt_s * spin_before,
        scaled_compliance_m_per_n: mobility.map(|row| row.map(|v| half_dt_squared * v)),
    };
    let resolved = law.advance_phase(
        &input,
        accepted_phase(&[slip_before[0], slip_before[1], spin_before]),
        controls.contact,
        || gate.is_requested(),
    )?;
    let force = axes[0]
        .scale(resolved.force_n[0])
        .add(axes[1].scale(resolved.force_n[1]));
    let couple = axes[2].scale(resolved.spin_moment_nm);
    let impulse = force.scale(dt_s);
    let angular_impulse = couple.scale(dt_s);
    let after = apply_pair_wrench(before, properties, arms, impulse, angular_impulse, gate)?;
    let endpoint = points(after, properties, arms)?;
    let relative_after = relative_point(endpoint);
    let angular_after = relative_spin(endpoint);
    let midpoint = relative_before.add(relative_after).scale(0.5);
    let angular_midpoint = angular_before.add(angular_after).scale(0.5);
    let slip_midpoint = [midpoint.dot(axes[0]), midpoint.dot(axes[1])];
    let spin_midpoint = angular_midpoint.dot(axes[2]);
    let work = impulse.dot(midpoint) + angular_impulse.dot(angular_midpoint);
    let (heat, force_residual, moment_residual) = if resolved.mode == AdhesiveWrenchMode::Sticking {
        (0.0, 0.0, 0.0)
    } else {
        let native = contact
            .source
            .patch()
            .evaluate(slip_midpoint, spin_midpoint)
            .map_err(AdhesiveWrenchError::Law)?;
        (
            dt_s * native.dissipated_power_w,
            (native.force_n[0] - resolved.force_n[0])
                .hypot(native.force_n[1] - resolved.force_n[1]),
            (native.spin_moment_nm - resolved.spin_moment_nm).abs(),
        )
    };
    let (energy_after, linear_after, momentum_after) = totals(&after, &properties, gate)?;
    let report = RigidAdhesiveReport {
        contact: Some(resolved),
        scaled_mobility_per_kg: Some(mobility),
        force_impulse_world_ns: impulse,
        spin_angular_impulse_world_nms: angular_impulse,
        slip_velocity_before_m_s: slip_before,
        slip_velocity_after_m_s: [relative_after.dot(axes[0]), relative_after.dot(axes[1])],
        spin_before_rad_s: spin_before,
        spin_after_rad_s: angular_after.dot(axes[2]),
        kinetic_before_j: energy_before,
        kinetic_after_j: energy_after,
        midpoint_impulse_work_j: work,
        dissipated_energy_j: heat,
        energy_residual_j: energy_after - energy_before + heat,
        force_residual_n: force_residual,
        spin_moment_residual_nm: moment_residual,
        linear_momentum_change_ns: linear_after.sub(linear_before),
        angular_momentum_change_nms: momentum_after.sub(momentum_before),
    };
    let velocity_residual = (slip_midpoint[0] - resolved.slip_displacement_m[0] / dt_s)
        .abs()
        .max((slip_midpoint[1] - resolved.slip_displacement_m[1] / dt_s).abs());
    let spin_residual = (spin_midpoint - resolved.spin_angle_rad / dt_s).abs();
    let balances = [
        report.energy_residual_j,
        energy_after - energy_before - work,
        work + resolved.contact_work_j,
        heat - resolved.dissipated_energy_j,
    ];
    if balances
        .into_iter()
        .any(|x| !x.is_finite() || x.abs() > controls.energy_tolerance_j)
        || !force_residual.is_finite()
        || force_residual > controls.contact.force_tolerance_n
        || !moment_residual.is_finite()
        || moment_residual > controls.contact.moment_tolerance_nm
        || !velocity_residual.is_finite()
        || velocity_residual > controls.velocity_tolerance_m_s
        || !spin_residual.is_finite()
        || spin_residual > controls.angular_velocity_tolerance_rad_s
        || magnitude(report.linear_momentum_change_ns) > controls.linear_momentum_tolerance_ns
        || magnitude(report.angular_momentum_change_nms) > controls.angular_momentum_tolerance_nms
    {
        return Err(RigidAdhesiveError::Audit(Box::new(report)));
    }
    let assisting_work = [
        force.dot(relative_before),
        force.dot(relative_after),
        couple.dot(angular_before),
        couple.dot(angular_after),
    ]
    .into_iter()
    .fold(0.0, f64::max)
        * dt_s;
    if !assisting_work.is_finite() || assisting_work > controls.energy_tolerance_j {
        return Err(RigidAdhesiveError::UnresolvedSlipEvent {
            assisting_work_j: assisting_work,
        });
    }
    poll(gate)?;
    Ok(RigidAdhesiveCandidate {
        bodies_after: after,
        report,
    })
}
