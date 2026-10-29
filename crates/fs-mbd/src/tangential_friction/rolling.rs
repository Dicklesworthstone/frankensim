//! Two-axis viscoelastic sphere rolling on original finite rigid-body momenta.
//!
//! The original material law supplies `M = -c omega`, `c = A N R`. Isotropy
//! extends its scalar response along the actual tangent relative angular
//! velocity. Two workless point-traction reactions preserve rolling without
//! slip; their magnitudes must fit the original static friction capacity.
//!
//! With the native contact mobility partitioned into shear and rolling blocks,
//! `F = -W_tt^-1 W_tr M` and `omega_dot = K M`, where
//! `K = W_rr - W_rt W_tt^-1 W_tr`. The small matrices come from the original
//! masses, orientations, principal inertias and contact arms. No effective
//! inertia, rolling speed, damping coefficient or thermal store is authored.
//!
//! This is a frozen-configuration supported-contact phase, not a pose or normal
//! integrator. Initially slipping/spinning contacts, induced normal separation
//! or relative spin, and insufficient static traction refuse. Another contact
//! lane must resolve those modes. There is no added shear/spin dissipation.
//! The normal patch, temperature and frequency stay held over this phase.

use super::{RigidContactBinding, RigidFrictionError, components, magnitude, poll, totals, vector};
use crate::{DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3};
use fs_exec::CancelGate;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::rolling_loss::viscoelastic::DeformingSphereSurface;
use fs_tribo::rolling_loss::{
    RollingKinematics, RollingLossCheckpoint, RollingLossError, RollingLossLaw, RollingLossState,
    RollingLossStep, RollingPatchReceipt, RollingWorkOwnership,
};
use fs_tribo::{FrictionLaw, InterfaceSystemRef, TriboError};

/// One original sphere/rigid-plane patch with both bodies free to recoil.
/// The frame normal points from body B toward body A. The sphere's COM-relative
/// arm must be minus that normal times its actual radius on A, or plus on B.
#[derive(Debug, Clone, Copy)]
pub struct RigidRollingContact<'a> {
    pub binding: RigidContactBinding,
    pub frame: TangentFrame,
    pub patch: &'a RollingPatchReceipt,
    pub interface: &'a InterfaceSystemRef,
    /// Must be the original material-derived `ViscoelasticSphere` variant.
    pub law: &'a RollingLossLaw,
    /// Only its native static capacity is used; no sliding law is added here.
    pub static_friction: &'a FrictionLaw,
    /// Actual sphere material temperature, held over the phase [K].
    pub temperature_k: f64,
    /// Original material-card frequency coordinate [Hz], never inferred from dt.
    pub excitation_frequency_hz: f64,
    /// One separate rolling-vector work channel, excluding shear and spin work.
    pub work: &'a RollingWorkOwnership,
}

/// Original rolling history and its original source-bound checkpoint.
/// No mechanical or caloric state is copied into this owner.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct RigidRollingHistory {
    state: RollingLossState,
    checkpoint: Option<RollingLossCheckpoint>,
}
impl RigidRollingHistory {
    /// A fresh, zero-heat history; its first accepted phase binds the source.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
    /// Restore the original checkpoint without losing its cumulative physical heat.
    pub fn from_checkpoint(
        contact: RigidRollingContact<'_>,
        checkpoint: RollingLossCheckpoint,
    ) -> Result<Self, RigidRollingError> {
        let state = contact.law.restore_checkpoint(
            contact.patch,
            contact.interface,
            contact.work,
            &checkpoint,
        )?;
        Ok(Self {
            state,
            checkpoint: Some(checkpoint),
        })
    }
    #[must_use]
    pub const fn state(&self) -> &RollingLossState {
        &self.state
    }
    #[must_use]
    pub const fn checkpoint(&self) -> Option<&RollingLossCheckpoint> {
        self.checkpoint.as_ref()
    }
}

/// Explicit accuracy/work limits and independent native-mechanics tolerances.
#[derive(Debug, Clone, Copy)]
pub struct RigidRollingControls {
    pub maximum_bodies: usize,
    /// At most one million midpoint substeps; exhaustion refuses before mutation.
    pub maximum_substeps: usize,
    /// User time-resolution ceiling [s]. A stricter monotone-decay bound is
    /// derived from the actual rolling mobility and original material law.
    pub maximum_step_s: f64,
    pub geometry_tolerance_m: f64,
    pub velocity_tolerance_m_s: f64,
    pub angular_velocity_tolerance_rad_s: f64,
    pub linear_momentum_tolerance_ns: f64,
    pub angular_momentum_tolerance_nms: f64,
    pub energy_tolerance_j: f64,
}

/// Physical impulse/work observables of the accepted supported rolling phase.
#[derive(Debug, Clone, PartialEq)]
pub struct RigidRollingReport {
    /// The only heat-generating body for this bulk viscoelastic model.
    pub deforming_body: usize,
    /// Original `[shear_x, shear_y, roll_x, roll_y]` response, in mixed SI units.
    pub mobility: [[f64; 4]; 4],
    /// Derived constrained rolling response [1/(kg m²)], retaining cross terms.
    pub rolling_mobility: [[f64; 2]; 2],
    pub rolling_rate_before_rad_s: [f64; 2],
    pub rolling_rate_after_rad_s: [f64; 2],
    /// Integrated force on A at its actual contact point [N s]; B receives minus.
    pub traction_impulse_world_ns: Vec3,
    /// Integrated rolling couple on A [N m s]; B receives minus.
    pub rolling_angular_impulse_world_nms: Vec3,
    pub static_capacity_n: f64,
    pub maximum_static_force_n: f64,
    pub substeps: usize,
    pub dissipated_heat_j: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub midpoint_impulse_work_j: f64,
    /// Work of the no-slip reactions; independently required to vanish.
    pub static_reaction_work_j: f64,
    pub work_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    pub maximum_point_speed_m_s: f64,
    pub maximum_relative_spin_rad_s: f64,
    pub endpoint_rolling_residual_rad_s: f64,
    /// Original response for the final substep only. Total heat is above;
    /// its native next state/checkpoint retain all accepted substep heat.
    pub final_step: RollingLossStep,
}

/// All errors leave the original bodies and rolling history unchanged.
#[derive(Debug)]
pub enum RigidRollingError {
    Invalid(&'static str),
    Native(RollingLossError),
    Friction(TriboError),
    Rigid(RigidFrictionError),
    Dynamics(DynamicsError),
    /// The requested time resolution/monotonicity exceeds the explicit work cap.
    SubstepBudget {
        required: f64,
        maximum: usize,
    },
    /// No-slip reaction cannot be supported by the original static friction law.
    StaticTraction {
        required_n: f64,
        capacity_n: f64,
    },
    Audit(Box<RigidRollingReport>),
}
impl From<RollingLossError> for RigidRollingError {
    fn from(e: RollingLossError) -> Self {
        Self::Native(e)
    }
}
impl From<TriboError> for RigidRollingError {
    fn from(e: TriboError) -> Self {
        Self::Friction(e)
    }
}
impl From<RigidFrictionError> for RigidRollingError {
    fn from(e: RigidFrictionError) -> Self {
        Self::Rigid(e)
    }
}
impl From<DynamicsError> for RigidRollingError {
    fn from(e: DynamicsError) -> Self {
        Self::Dynamics(e)
    }
}
impl core::fmt::Display for RigidRollingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid viscoelastic rolling: {self:?}")
    }
}
impl std::error::Error for RigidRollingError {}

/// Advance both bodies and the original rolling-loss state atomically.
///
/// The original Hertz owner verifies patch load/area/radius/material agreement.
/// Native midpoint moment work is exactly the exported physical heat, up to the
/// independently checked arithmetic gates. No numerical damping is deposited.
/// Every eigenmode of the constrained viscous response has a nonnegative
/// midpoint decay factor: automatic substeps enforce `c lambda_max h <= 1`.
/// Viscous arrest is asymptotic; zero speed stays zero without a cutoff or snap.
/// Accuracy at frozen geometry is second order under timestep refinement.
///
/// A surrounding trajectory owner must advance the poses and normal contact.
/// This pure-rolling phase cannot be superposed with a second shear reaction
/// or microslip model on the same phase; those modes need a joint contact solve.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn apply_rigid_rolling(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    history: &mut RigidRollingHistory,
    contact: RigidRollingContact<'_>,
    dt_s: f64,
    controls: RigidRollingControls,
    gate: &CancelGate,
) -> Result<RigidRollingReport, RigidRollingError> {
    poll(gate)?;
    validate(bodies, properties, contact, dt_s, controls)?;
    let RollingLossLaw::ViscoelasticSphere(card) = contact.law else {
        return Err(invalid(
            "this consumer requires the original viscoelastic sphere law",
        ));
    };
    if let Some(checkpoint) = &history.checkpoint {
        let restored = contact.law.restore_checkpoint(
            contact.patch,
            contact.interface,
            contact.work,
            checkpoint,
        )?;
        if restored != history.state {
            return Err(invalid(
                "rolling history differs from its original checkpoint",
            ));
        }
    } else if history.state != RollingLossState::zero() {
        return Err(invalid(
            "nonzero rolling history needs its original checkpoint",
        ));
    }
    let axes = [
        vector(contact.frame.longitudinal()),
        vector(contact.frame.lateral()),
    ];
    let normal = vector(contact.frame.normal());
    let before = points(bodies, properties, contact.binding)?;
    let sphere = card.parameters().sphere_surface;
    let (sphere_slot, sphere_sign, deforming_body) = match sphere {
        DeformingSphereSurface::SurfaceA => (0, -1.0, contact.binding.body_a),
        DeformingSphereSurface::SurfaceB => (1, 1.0, contact.binding.body_b),
    };
    let radius = card.parameters().radius_m;
    if magnitude(before[0].point_world.sub(before[1].point_world)) > controls.geometry_tolerance_m
        || magnitude(
            before[sphere_slot]
                .arm_world
                .sub(normal.scale(sphere_sign * radius)),
        ) > controls.geometry_tolerance_m
    {
        return Err(invalid(
            "original coincident points and actual sphere radius/normal required",
        ));
    }
    let loaded = contact.patch.normal_load_n() > 0.0;
    let mut maximum_point_speed = magnitude(relative_point(before));
    let mut maximum_spin = relative_spin(before).dot(normal).abs();
    if loaded
        && (maximum_point_speed > controls.velocity_tolerance_m_s
            || maximum_spin > controls.angular_velocity_tolerance_rad_s)
    {
        return Err(invalid(
            "supported rolling requires no initial slip, separation or relative spin",
        ));
    }
    let q_before = tangent(relative_spin(before), axes);
    let kinematics = |rate: f64, interval: f64| {
        RollingKinematics::new(
            radius * rate,
            rate,
            0.0,
            contact.temperature_k,
            contact.excitation_frequency_hz,
            interval,
        )
    };
    // Check the fastest accepted state against the original Hertz/material
    // derivation; accepting only the smaller midpoint rate would widen its box.
    let initial_response = card.response(
        contact.patch,
        contact.interface,
        kinematics(norm(q_before), dt_s)?,
    )?;
    let capacity =
        finite(contact.static_friction.static_coefficient()? * contact.patch.normal_load_n())?;
    let (mobility, reaction, rolling_mobility) =
        response(bodies, properties, contact.binding, axes)?;
    let c = finite((card.rolling_time_s() * radius) * contact.patch.normal_load_n())?;
    if c < 0.0 || (c == 0.0 && card.rolling_time_s() > 0.0 && loaded) {
        return Err(invalid("unrepresentable material rolling damping"));
    }
    let largest = finite(
        0.5 * rolling_mobility[0][0]
            + 0.5 * rolling_mobility[1][1]
            + (0.5 * (rolling_mobility[0][0] - rolling_mobility[1][1]))
                .hypot(rolling_mobility[0][1]),
    )?;
    let active = c > 0.0 && norm(q_before) > 0.0;
    let decay_rate = if active { finite(c * largest)? } else { 0.0 };
    let limit = if decay_rate > 0.0 {
        controls.maximum_step_s.min(decay_rate.recip())
    } else {
        controls.maximum_step_s
    };
    let required = if active {
        (dt_s / limit).ceil().max(1.0)
    } else {
        1.0
    };
    if !required.is_finite() || required > controls.maximum_substeps as f64 {
        return Err(RigidRollingError::SubstepBudget {
            required,
            maximum: controls.maximum_substeps,
        });
    }
    let substeps = required as usize;
    let h = dt_s / substeps as f64;
    if !h.is_finite() || h <= 0.0 {
        return Err(invalid("unrepresentable rolling substep duration"));
    }
    let half_damping = if active { 0.5 * h * c } else { 0.0 };
    let midpoint_matrix = [
        [
            1.0 + half_damping * rolling_mobility[0][0],
            half_damping * rolling_mobility[0][1],
        ],
        [
            half_damping * rolling_mobility[1][0],
            1.0 + half_damping * rolling_mobility[1][1],
        ],
    ];
    let (energy_before, linear_before, angular_before) = totals(bodies, properties, gate)?;
    let mut next = bodies.to_vec();
    let mut next_history = history.clone();
    let mut final_step = None;
    let mut heat = 0.0;
    let mut midpoint_work = 0.0;
    let mut static_work = 0.0;
    let mut maximum_force = require_capacity(
        reaction,
        q_before,
        initial_response.rolling_moment_nm,
        capacity,
    )?;
    let mut endpoint_residual = 0.0_f64;
    let mut traction_impulse = Vec3::ZERO;
    let mut rolling_impulse = Vec3::ZERO;
    for _ in 0..substeps {
        poll(gate)?;
        let old_points = points(&next, properties, contact.binding)?;
        let q = tangent(relative_spin(old_points), axes);
        let q_mid = solve_spd(midpoint_matrix, q)?;
        let rate = norm(q_mid);
        let step = contact.law.advance(
            contact.patch,
            contact.interface,
            kinematics(rate, h)?,
            contact.work,
            &next_history.state,
        )?;
        let moment = along(q_mid, step.wrench.rolling_moment_nm);
        let force = mul(reaction, moment);
        let demand = norm(force);
        if !demand.is_finite() || demand > capacity {
            return Err(RigidRollingError::StaticTraction {
                required_n: demand,
                capacity_n: capacity,
            });
        }
        maximum_force = maximum_force.max(demand);
        let impulse = world(force, axes).scale(h);
        let couple = world(moment, axes).scale(h);
        for (body, arm, sign) in [
            (contact.binding.body_a, contact.binding.arm_a_body, 1.0),
            (contact.binding.body_b, contact.binding.arm_b_body, -1.0),
        ] {
            let after_force = next[body]
                .apply_impulse_at_body_point(properties[body], arm, impulse.scale(sign))?
                .state_after;
            let angular = after_force
                .pose()
                .orientation()
                .rotate_world_to_body(couple.scale(sign));
            next[body] = RigidBodyState::new(
                after_force.pose(),
                after_force.linear_momentum_world(),
                after_force.angular_momentum_body().add(angular),
            )?;
        }
        let new_points = points(&next, properties, contact.binding)?;
        let q_after = tangent(relative_spin(new_points), axes);
        let predicted = mul(rolling_mobility, moment);
        for axis in 0..2 {
            endpoint_residual =
                endpoint_residual.max((q_after[axis] - q[axis] - h * predicted[axis]).abs());
        }
        maximum_point_speed = maximum_point_speed.max(magnitude(relative_point(new_points)));
        maximum_spin = maximum_spin.max(relative_spin(new_points).dot(normal).abs());
        if loaded
            && (maximum_point_speed > controls.velocity_tolerance_m_s
                || maximum_spin > controls.angular_velocity_tolerance_rad_s)
        {
            return Err(invalid(
                "rolling recoil induced unresolved slip, separation or relative spin",
            ));
        }
        // The material approximation is checked at actual endpoints as well as
        // the actual midpoint constitutive state, including zero-loss cases.
        let response = card.response(
            contact.patch,
            contact.interface,
            kinematics(norm(q_after), h)?,
        )?;
        maximum_force = maximum_force.max(require_capacity(
            reaction,
            q_after,
            response.rolling_moment_nm,
            capacity,
        )?);
        let force_work = impulse.dot(
            relative_point(old_points)
                .add(relative_point(new_points))
                .scale(0.5),
        );
        let moment_work = couple.dot(
            relative_spin(old_points)
                .add(relative_spin(new_points))
                .scale(0.5),
        );
        if !force_work.is_finite()
            || !moment_work.is_finite()
            || (moment_work + step.dissipation.total_heat_j).abs() > controls.energy_tolerance_j
            || step.wrench.contour_force_n != 0.0
            || step.wrench.spin_moment_nm != 0.0
            || step.stored_energy_j != 0.0
            || step.storage_change_j != 0.0
        {
            return Err(invalid(
                "original rolling moment work and physical heat do not agree",
            ));
        }
        heat = finite(heat + step.dissipation.total_heat_j)?;
        static_work = finite(static_work + force_work)?;
        midpoint_work = finite(midpoint_work + force_work + moment_work)?;
        traction_impulse = traction_impulse.add(impulse);
        rolling_impulse = rolling_impulse.add(couple);
        next_history.state = step.next_state.clone();
        final_step = Some(step);
    }
    let (energy_after, linear_after, angular_after) = totals(&next, properties, gate)?;
    let final_step = final_step.ok_or_else(|| invalid("missing native rolling substep"))?;
    next_history.checkpoint = Some(final_step.checkpoint.clone());
    let report = RigidRollingReport {
        deforming_body,
        mobility,
        rolling_mobility,
        rolling_rate_before_rad_s: q_before,
        rolling_rate_after_rad_s: tangent(
            relative_spin(points(&next, properties, contact.binding)?),
            axes,
        ),
        traction_impulse_world_ns: traction_impulse,
        rolling_angular_impulse_world_nms: rolling_impulse,
        static_capacity_n: capacity,
        maximum_static_force_n: maximum_force,
        substeps,
        dissipated_heat_j: heat,
        kinetic_before_j: energy_before,
        kinetic_after_j: energy_after,
        midpoint_impulse_work_j: midpoint_work,
        static_reaction_work_j: static_work,
        work_energy_residual_j: energy_after - energy_before - midpoint_work,
        total_energy_residual_j: energy_after - energy_before + heat,
        linear_momentum_change_ns: linear_after.sub(linear_before),
        angular_momentum_change_nms: angular_after.sub(angular_before),
        maximum_point_speed_m_s: maximum_point_speed,
        maximum_relative_spin_rad_s: maximum_spin,
        endpoint_rolling_residual_rad_s: endpoint_residual,
        final_step,
    };
    let history_heat = next_history.state.cumulative_heat_j() - history.state.cumulative_heat_j();
    if [
        energy_before,
        energy_after,
        heat,
        midpoint_work,
        static_work,
        history_heat,
        report.work_energy_residual_j,
        report.total_energy_residual_j,
        endpoint_residual,
        magnitude(report.linear_momentum_change_ns),
        magnitude(report.angular_momentum_change_nms),
        magnitude(traction_impulse),
        magnitude(rolling_impulse),
    ]
    .iter()
    .any(|v| !v.is_finite())
        || heat < 0.0
        || (history_heat - heat).abs() > controls.energy_tolerance_j
        || static_work.abs() > controls.energy_tolerance_j
        || report.work_energy_residual_j.abs() > controls.energy_tolerance_j
        || report.total_energy_residual_j.abs() > controls.energy_tolerance_j
        || endpoint_residual > controls.angular_velocity_tolerance_rad_s
        || magnitude(report.linear_momentum_change_ns) > controls.linear_momentum_tolerance_ns
        || magnitude(report.angular_momentum_change_nms) > controls.angular_momentum_tolerance_nms
    {
        return Err(RigidRollingError::Audit(Box::new(report)));
    }
    poll(gate)?;
    bodies.copy_from_slice(&next);
    *history = next_history;
    Ok(report)
}

fn validate(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidRollingContact<'_>,
    dt: f64,
    c: RigidRollingControls,
) -> Result<(), RigidRollingError> {
    if bodies.len() != properties.len()
        || bodies.len() > c.maximum_bodies
        || c.maximum_bodies == 0
        || !(1..=1_000_000).contains(&c.maximum_substeps)
        || contact.binding.body_a >= bodies.len()
        || contact.binding.body_b >= bodies.len()
        || contact.binding.body_a == contact.binding.body_b
        || [
            dt,
            c.maximum_step_s,
            c.geometry_tolerance_m,
            c.velocity_tolerance_m_s,
            c.angular_velocity_tolerance_rad_s,
            c.linear_momentum_tolerance_ns,
            c.angular_momentum_tolerance_nms,
            c.energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(invalid(
            "bounded original bodies, positive duration and finite gates required",
        ));
    }
    Ok(())
}

fn points(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    binding: RigidContactBinding,
) -> Result<[PointKinematics; 2], RigidRollingError> {
    Ok([
        bodies[binding.body_a].point_kinematics(properties[binding.body_a], binding.arm_a_body)?,
        bodies[binding.body_b].point_kinematics(properties[binding.body_b], binding.arm_b_body)?,
    ])
}
fn relative_point(p: [PointKinematics; 2]) -> Vec3 {
    p[0].point_velocity_world.sub(p[1].point_velocity_world)
}
fn relative_spin(p: [PointKinematics; 2]) -> Vec3 {
    p[0].angular_velocity_world.sub(p[1].angular_velocity_world)
}
fn tangent(v: Vec3, axes: [Vec3; 2]) -> [f64; 2] {
    [v.dot(axes[0]), v.dot(axes[1])]
}
fn world(v: [f64; 2], axes: [Vec3; 2]) -> Vec3 {
    axes[0].scale(v[0]).add(axes[1].scale(v[1]))
}
fn norm(v: [f64; 2]) -> f64 {
    v[0].hypot(v[1])
}
fn along(direction: [f64; 2], magnitude: f64) -> [f64; 2] {
    let length = norm(direction);
    if length == 0.0 {
        [0.0; 2]
    } else {
        direction.map(|v| magnitude * (v / length))
    }
}
fn require_capacity(
    reaction: [[f64; 2]; 2],
    rate: [f64; 2],
    moment: f64,
    capacity: f64,
) -> Result<f64, RigidRollingError> {
    let demand = norm(mul(reaction, along(rate, moment)));
    if !demand.is_finite() || demand > capacity {
        return Err(RigidRollingError::StaticTraction {
            required_n: demand,
            capacity_n: capacity,
        });
    }
    Ok(demand)
}
fn mul(m: [[f64; 2]; 2], v: [f64; 2]) -> [f64; 2] {
    [
        m[0][0] * v[0] + m[0][1] * v[1],
        m[1][0] * v[0] + m[1][1] * v[1],
    ]
}
fn finite(v: f64) -> Result<f64, RigidRollingError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(invalid("unrepresentable rolling response"))
    }
}
fn invalid(message: &'static str) -> RigidRollingError {
    RigidRollingError::Invalid(message)
}

/// Scaled 2x2 Cholesky; no friction or body dynamics is approximated here.
fn solve_spd(m: [[f64; 2]; 2], rhs: [f64; 2]) -> Result<[f64; 2], RigidRollingError> {
    let scale = m[0][0].max(m[1][1]);
    if !scale.is_finite() || scale <= 0.0 {
        return Err(invalid(
            "positive definite representable rolling response required",
        ));
    }
    let a = (m[0][0] / scale).sqrt();
    let b = (m[1][0] / scale) / a;
    let d2 = m[1][1] / scale - b * b;
    if !a.is_finite() || a <= 0.0 || !d2.is_finite() || d2 <= 0.0 {
        return Err(invalid("singular rolling or no-slip response"));
    }
    let d = d2.sqrt();
    let y0 = (rhs[0] / scale) / a;
    let y1 = (rhs[1] / scale - b * y0) / d;
    let x1 = y1 / d;
    Ok([finite((y0 - b * x1) / a)?, finite(x1)?])
}

/// Build J M^-1 J^T from the original two bodies, then eliminate the two
/// workless shear reactions. Projected Gram rows avoid subtracting nearly
/// equal positive inertias when forming the Schur complement.
fn response(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    binding: RigidContactBinding,
    axes: [Vec3; 2],
) -> Result<([[f64; 4]; 4], [[f64; 2]; 2], [[f64; 2]; 2]), RigidRollingError> {
    let mut factor = [[0.0; 12]; 4];
    for (slot, body, arm, sign) in [
        (0, binding.body_a, binding.arm_a_body, 1.0),
        (1, binding.body_b, binding.arm_b_body, -1.0),
    ] {
        let inertia = components(properties[body].principal_inertia_body());
        let inverse_roots = inertia.map(|v| v.sqrt().recip());
        let inverse_mass_root = properties[body].mass().sqrt().recip();
        let orientation = bodies[body].pose().orientation();
        for axis in 0..2 {
            let direction = orientation.rotate_world_to_body(axes[axis]);
            let translation = components(axes[axis].scale(sign * inverse_mass_root));
            let force_rotation = components(arm.cross(direction).scale(sign));
            let moment_rotation = components(direction.scale(sign));
            for j in 0..3 {
                factor[axis][6 * slot + j] = translation[j];
                factor[axis][6 * slot + 3 + j] = force_rotation[j] * inverse_roots[j];
                factor[axis + 2][6 * slot + 3 + j] = moment_rotation[j] * inverse_roots[j];
            }
        }
    }
    let dot = |a: &[f64; 12], b: &[f64; 12]| -> Result<f64, RigidRollingError> {
        finite(a.iter().zip(b).map(|(x, y)| x * y).sum())
    };
    let mut mobility = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..=i {
            mobility[i][j] = dot(&factor[i], &factor[j])?;
            mobility[j][i] = mobility[i][j];
        }
    }
    let shear = [
        [mobility[0][0], mobility[0][1]],
        [mobility[1][0], mobility[1][1]],
    ];
    let mut reaction = [[0.0; 2]; 2];
    let mut projected = [[0.0; 12]; 2];
    for axis in 0..2 {
        let solution = solve_spd(shear, [-mobility[0][axis + 2], -mobility[1][axis + 2]])?;
        for j in 0..2 {
            reaction[j][axis] = solution[j];
        }
        for j in 0..12 {
            projected[axis][j] =
                factor[axis + 2][j] + solution[0] * factor[0][j] + solution[1] * factor[1][j];
        }
    }
    let off = dot(&projected[0], &projected[1])?;
    let rolling = [
        [dot(&projected[0], &projected[0])?, off],
        [off, dot(&projected[1], &projected[1])?],
    ];
    solve_spd(rolling, [0.0; 2])?;
    Ok((mobility, reaction, rolling))
}
