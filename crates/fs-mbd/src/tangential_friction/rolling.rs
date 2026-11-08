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

use super::{
    RigidContactBinding, RigidFrictionError, components, kinetic_change, magnitude, poll, totals,
    vector,
};
use crate::{DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3};
use fs_exec::CancelGate;
use fs_math::{dd::Dd, det};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::rolling_loss::viscoelastic::DeformingSphereSurface;
use fs_tribo::rolling_loss::viscoelastic::path::{ExponentialPathError, ExponentialRatePath};
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

    /// Internal composition retains the original admitted source and history;
    /// it cannot synthesize a checkpoint from separate numerical observables.
    pub(super) fn from_native_step(step: &RollingLossStep) -> Self {
        Self {
            state: step.next_state.clone(),
            checkpoint: Some(step.checkpoint.clone()),
        }
    }
}

/// Explicit accuracy/work limits and independent native-mechanics tolerances.
#[derive(Debug, Clone, Copy)]
pub struct RigidRollingControls {
    pub maximum_bodies: usize,
    /// At most one million constitutive substeps; exhaustion refuses before mutation.
    pub maximum_substeps: usize,
    /// User time-resolution ceiling [s]. Each held-geometry substep integrates
    /// both rolling eigenmodes exactly, with no damping-stiffness step limit.
    pub maximum_step_s: f64,
    pub geometry_tolerance_m: f64,
    pub velocity_tolerance_m_s: f64,
    pub angular_velocity_tolerance_rad_s: f64,
    pub linear_momentum_tolerance_ns: f64,
    pub angular_momentum_tolerance_nms: f64,
    pub energy_tolerance_j: f64,
}

/// Physical impulse/work observables of the accepted supported rolling phase.
/// Kinetic totals and momentum changes cover only the contacted pair; its
/// changes are also the whole-system changes because other bodies are untouched.
#[derive(Debug, Clone, PartialEq)]
pub struct RigidRollingReport {
    /// The only heat-generating body's original input-slice index.
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
    /// Contacted-pair kinetic energy before the phase [J].
    pub kinetic_before_j: f64,
    /// Contacted-pair kinetic energy after the phase [J].
    pub kinetic_after_j: f64,
    /// Actual momentum-based change; retains losses beside common motion [J].
    pub kinetic_energy_change_j: f64,
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
    Path(ExponentialPathError),
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
impl From<ExponentialPathError> for RigidRollingError {
    fn from(e: ExponentialPathError) -> Self {
        Self::Path(e)
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
/// Native impulse work is exactly the exported physical heat, up to the
/// independently checked arithmetic gates. The two orthogonal modes of the
/// original constrained mobility use exact material-owned exponential paths,
/// including signed moment impulses and squared-rate heat integrals. Static
/// traction and unsupported normal/spin motion are checked over each complete
/// path, including interior extrema. No numerical damping is deposited.
/// Viscous arrest is asymptotic; zero speed stays zero without a cutoff or snap.
/// Both isotropic and anisotropic decay are exact at frozen geometry.
/// Only the bound pair is staged and audited, so unrelated kinetic energy
/// cannot hide its energy change. Every other original body stays untouched.
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
    let mut next = [
        bodies[contact.binding.body_a],
        bodies[contact.binding.body_b],
    ];
    let pair_properties = [
        properties[contact.binding.body_a],
        properties[contact.binding.body_b],
    ];
    let pair_binding = RigidContactBinding {
        body_a: 0,
        body_b: 1,
        ..contact.binding
    };
    let before = points(&next, &pair_properties, pair_binding)?;
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
    let mut maximum_point_speed = magnitude(relative_point(before)?);
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
    let ContactResponse {
        mobility,
        reaction,
        rolling_mobility,
        unresolved_mobility,
    } = response(&next, &pair_properties, pair_binding, axes, normal)?;
    let c = finite((card.rolling_time_s() * radius) * contact.patch.normal_load_n())?;
    if c < 0.0 || (c == 0.0 && card.rolling_time_s() > 0.0 && loaded) {
        return Err(invalid("unrepresentable material rolling damping"));
    }
    let active = c > 0.0 && norm(q_before) > 0.0;
    let modes = if active {
        let (values, basis) = rolling_modes(rolling_mobility)?;
        Some((basis, [finite(c * values[0])?, finite(c * values[1])?]))
    } else {
        None
    };
    let required = if active {
        (dt_s / controls.maximum_step_s).ceil().max(1.0)
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
    let before = next;
    let (energy_before, linear_before, angular_before) = totals(&before, &pair_properties, gate)?;
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
        let old_points = points(&next, &pair_properties, pair_binding)?;
        let q = tangent(relative_spin(old_points), axes);
        let (step, moment_impulse, predicted_endpoint) = if let Some((basis, decay)) = modes {
            let modal_rate = [dot2(basis[0], q), dot2(basis[1], q)];
            let paths = [
                ExponentialRatePath::new(modal_rate[0], 0.0, decay[0], h, || gate.is_requested())?,
                ExponentialRatePath::new(modal_rate[1], 0.0, decay[1], h, || gate.is_requested())?,
            ];
            let native = contact.law.advance_orthogonal_exponential_sphere_paths(
                contact.patch,
                contact.interface,
                &paths,
                0.0,
                contact.temperature_k,
                contact.excitation_frequency_hz,
                contact.work,
                &next_history.state,
                || gate.is_requested(),
            )?;
            let modal_moment = [-c * modal_rate[0], -c * modal_rate[1]];
            let force = [
                mul(reaction, basis[0].map(|v| v * modal_moment[0])),
                mul(reaction, basis[1].map(|v| v * modal_moment[1])),
            ];
            let demand = decaying_norm_max(force, decay, h)?;
            if demand > capacity {
                return Err(RigidRollingError::StaticTraction {
                    required_n: demand,
                    capacity_n: capacity,
                });
            }
            maximum_force = maximum_force.max(demand);
            // Shear stays constrained. The omitted normal and spin coordinates
            // must also remain admissible between endpoints; two distinct
            // eigenmodes can otherwise cancel only at the accepted endpoint.
            let point = relative_point(old_points)?;
            let mut unsupported = [0.0; 2];
            for (i, initial) in [point.dot(normal), relative_spin(old_points).dot(normal)]
                .into_iter()
                .enumerate()
            {
                let rates = [
                    finite(dot2(unresolved_mobility[i], basis[0]) * modal_moment[0])?,
                    finite(dot2(unresolved_mobility[i], basis[1]) * modal_moment[1])?,
                ];
                unsupported[i] = integrated_rate_max(initial, rates, decay, h, gate)?;
            }
            maximum_point_speed =
                maximum_point_speed.max(norm(tangent(point, axes)).hypot(unsupported[0]));
            maximum_spin = maximum_spin.max(unsupported[1]);
            let impulse = combine(native.rolling_moment_impulse_nms, basis);
            let endpoint = combine(paths.map(|path| path.endpoint_rate()), basis);
            (native.rolling, impulse, endpoint)
        } else {
            let step = contact.law.advance(
                contact.patch,
                contact.interface,
                kinematics(norm(q), h)?,
                contact.work,
                &next_history.state,
            )?;
            (step, [0.0; 2], q)
        };
        let impulse = world(mul(reaction, moment_impulse), axes);
        let couple = world(moment_impulse, axes);
        for (body, arm, sign) in [
            (pair_binding.body_a, pair_binding.arm_a_body, 1.0),
            (pair_binding.body_b, pair_binding.arm_b_body, -1.0),
        ] {
            let after_force = next[body]
                .apply_impulse_at_body_point(pair_properties[body], arm, impulse.scale(sign))?
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
        let new_points = points(&next, &pair_properties, pair_binding)?;
        let q_after = tangent(relative_spin(new_points), axes);
        for axis in 0..2 {
            endpoint_residual =
                endpoint_residual.max((q_after[axis] - predicted_endpoint[axis]).abs());
        }
        maximum_point_speed = maximum_point_speed.max(magnitude(relative_point(new_points)?));
        maximum_spin = maximum_spin.max(relative_spin(new_points).dot(normal).abs());
        if loaded
            && (maximum_point_speed > controls.velocity_tolerance_m_s
                || maximum_spin > controls.angular_velocity_tolerance_rad_s)
        {
            return Err(invalid(
                "rolling recoil induced unresolved slip, separation or relative spin",
            ));
        }
        // Check the material approximation at the actual momentum endpoint as
        // well as the complete analytic path, including zero-loss cases.
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
            relative_point(old_points)?
                .add(relative_point(new_points)?)
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
    let (energy_after, linear_after, angular_after) = totals(&next, &pair_properties, gate)?;
    let kinetic_energy_change_j = kinetic_change(&before, &next, &pair_properties, gate)?;
    let final_step = final_step.ok_or_else(|| invalid("missing native rolling substep"))?;
    next_history.checkpoint = Some(final_step.checkpoint.clone());
    let report = RigidRollingReport {
        deforming_body,
        mobility,
        rolling_mobility,
        rolling_rate_before_rad_s: q_before,
        rolling_rate_after_rad_s: tangent(
            relative_spin(points(&next, &pair_properties, pair_binding)?),
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
        kinetic_energy_change_j,
        midpoint_impulse_work_j: midpoint_work,
        static_reaction_work_j: static_work,
        work_energy_residual_j: kinetic_energy_change_j - midpoint_work,
        total_energy_residual_j: kinetic_energy_change_j + heat,
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
    bodies[contact.binding.body_a] = next[0];
    bodies[contact.binding.body_b] = next[1];
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
fn relative_point(p: [PointKinematics; 2]) -> Result<Vec3, DynamicsError> {
    p[0].relative_velocity_world(p[1])
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
fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0].mul_add(b[0], a[1] * b[1])
}
fn combine(v: [f64; 2], basis: [[f64; 2]; 2]) -> [f64; 2] {
    [
        basis[0][0].mul_add(v[0], basis[1][0] * v[1]),
        basis[0][1].mul_add(v[0], basis[1][1] * v[1]),
    ]
}

/// One scaled Jacobi rotation, with no eigenvalue-gap division or tolerance
/// that would discard a small cross term. The positive Cholesky determinant
/// retains the smaller eigenvalue without subtracting two large eigenvalues.
fn rolling_modes(m: [[f64; 2]; 2]) -> Result<([f64; 2], [[f64; 2]; 2]), RigidRollingError> {
    if m[0][1] == 0.0 {
        return Ok(([m[0][0], m[1][1]], [[1.0, 0.0], [0.0, 1.0]]));
    }
    let scale = m[0][0].max(m[1][1]);
    let (a, b, d) = (m[0][0] / scale, m[0][1] / scale, m[1][1] / scale);
    if b == 0.0 {
        return Err(invalid("unrepresentable rolling eigenvector"));
    }
    let delta = 0.5 * a - 0.5 * d;
    let t = b / (delta.abs() + delta.hypot(b));
    let cosine = (1.0 + t * t).sqrt().recip();
    let sine = t * cosine;
    let (high, first) = if delta >= 0.0 {
        (a + b * t, [cosine, sine])
    } else {
        (d + b * t, [sine, cosine])
    };
    let cholesky = b / a.sqrt();
    let low = finite(m[0][0] * ((d - cholesky * cholesky) / high))?;
    let high = finite(scale * high)?;
    if low <= 0.0 || high <= 0.0 {
        return Err(invalid("unresolved positive rolling eigenvalues"));
    }
    Ok(([high, low], [first, [-first[1], first[0]]]))
}

/// Maximum norm of two decaying vectors over the whole held interval.
/// With x=exp(-(a-b)t), stationary squared norms are the positive roots of
/// a*|u|²*x² + (a+b)*(u.v)*x + b*|v|². Scaling and double-double products
/// retain nearly cancelling force directions; repeated rates need no roots.
fn decaying_norm_max(
    mut vectors: [[f64; 2]; 2],
    mut decay: [f64; 2],
    duration: f64,
) -> Result<f64, RigidRollingError> {
    if decay[0] < decay[1] {
        decay.swap(0, 1);
        vectors.swap(0, 1);
    }
    let at = |t: f64| -> Result<f64, RigidRollingError> {
        let relative_exponent = (decay[0] - decay[1]) * t;
        if relative_exponent <= 0.5 {
            // Separate exponentials can round to the same value although
            // releasing their cancellation produces a resolved physical force.
            let shift = Dd::from_f64(det::expm1(-relative_exponent));
            let slow = Dd::from_f64(det::exp(-decay[1] * t));
            let components = [0, 1].map(|i| {
                ((Dd::from_f64(vectors[0][i])
                    + Dd::from_f64(vectors[1][i])
                    + Dd::from_f64(vectors[0][i]) * shift)
                    * slow)
                    .to_f64()
            });
            return finite(norm(components));
        }
        let factors = decay.map(|rate| det::exp(-rate * t));
        finite(norm(combine(factors, vectors)))
    };
    let mut maximum = at(0.0)?.max(at(duration)?);
    if decay[0] == decay[1] {
        return Ok(maximum);
    }
    let scale = finite(norm(vectors[0]).max(norm(vectors[1])))?;
    if scale == 0.0 {
        return Ok(0.0);
    }
    let [u, v] = vectors.map(|row| row.map(|x| x / scale));
    let dot = |x: [f64; 2], y: [f64; 2]| {
        Dd::from_f64(x[0]) * Dd::from_f64(y[0]) + Dd::from_f64(x[1]) * Dd::from_f64(y[1])
    };
    let one = Dd::from_f64(1.0);
    let zero = Dd::from_f64(0.0);
    let ratio = Dd::from_f64(decay[1]) / Dd::from_f64(decay[0]);
    let a = dot(u, u);
    let uv = dot(u, v);
    let b = uv * (one + ratio);
    let c = dot(v, v) * ratio;
    if b.to_f64() >= 0.0 {
        return Ok(maximum);
    }
    let roots = if a.to_f64() == 0.0 {
        [-c / b, zero]
    } else {
        // The Gram identity avoids subtracting O(1) terms to recover the
        // squared rate gap for almost collinear, almost equally decaying forces.
        let cross =
            Dd::from_f64(u[0]) * Dd::from_f64(v[1]) - Dd::from_f64(u[1]) * Dd::from_f64(v[0]);
        let gap = (Dd::from_f64(decay[0]) - Dd::from_f64(decay[1])) / Dd::from_f64(decay[0]);
        let discriminant = gap * gap * uv * uv - ratio * cross * cross * Dd::from_f64(4.0);
        if discriminant.to_f64() <= 0.0 {
            return Ok(maximum);
        }
        let numerator = (-b + discriminant.sqrt()) * Dd::from_f64(0.5);
        [numerator / a, c / numerator]
    };
    for root in roots {
        if zero.lt(root) && root.lt(one) {
            let offset = (root - one).to_f64();
            let logarithm = if offset.abs() <= 0.0625 {
                log1p_near_zero(offset)
            } else {
                det::ln(root.to_f64())
            };
            let t = -logarithm / (decay[0] - decay[1]);
            if t > 0.0 && t < duration {
                maximum = maximum.max(at(t)?);
            }
        }
    }
    Ok(maximum)
}

// ln(1+x)=2 atanh(x/(2+x)); eight terms suffice on |x|<=1/16.
// Preserve the small increment instead of rounding 1+x to one first.
fn log1p_near_zero(x: f64) -> f64 {
    if x.abs() <= f64::EPSILON {
        return x;
    }
    let s = x / (2.0 + x);
    let square = s * s;
    let mut term = s;
    let mut sum = s;
    for k in 1..8 {
        term *= square;
        sum += term / (2 * k + 1) as f64;
    }
    2.0 * sum
}

/// An omitted scalar velocity is an initial value plus two integrated
/// exponentials. Its derivative has at most one interior zero. This checks
/// normal separation and relative spin without sampling or damping substeps.
fn integrated_rate_max(
    initial: f64,
    mut rates: [f64; 2],
    mut decay: [f64; 2],
    duration: f64,
    gate: &CancelGate,
) -> Result<f64, RigidRollingError> {
    if decay[0] < decay[1] {
        decay.swap(0, 1);
        rates.swap(0, 1);
    }
    let at = |t| -> Result<f64, RigidRollingError> {
        poll(gate)?;
        let change = if decay[1] > 0.0 {
            // Integrate the pair before rounding: independently rounded path
            // integrals can erase a small released cancellation completely.
            let p0 = Dd::from_f64(rates[0]) / Dd::from_f64(decay[0]);
            let p1 = Dd::from_f64(rates[1]) / Dd::from_f64(decay[1]);
            finite(p0.to_f64())?;
            finite(p1.to_f64())?;
            let slow_release = Dd::from_f64(-det::expm1(-decay[1] * t));
            let slow_remaining = Dd::from_f64(det::exp(-decay[1] * t));
            let relative_release = Dd::from_f64(-det::expm1(-(decay[0] - decay[1]) * t));
            (p0 + p1) * slow_release + p0 * slow_remaining * relative_release
        } else if decay[0] == 0.0 {
            (Dd::from_f64(rates[0]) + Dd::from_f64(rates[1])) * Dd::from_f64(t)
        } else {
            let path =
                ExponentialRatePath::new(rates[0], 0.0, decay[0], t, || gate.is_requested())?;
            if decay[0] * t <= 0.5 {
                // integral(q)-q(0)*t = -decay*integrated_displacement.
                (Dd::from_f64(rates[0]) + Dd::from_f64(rates[1])) * Dd::from_f64(t)
                    - Dd::from_f64(decay[0]) * Dd::from_f64(path.integrated_displacement())
            } else {
                Dd::from_f64(path.rate_integral()) + Dd::from_f64(rates[1]) * Dd::from_f64(t)
            }
        };
        let value = Dd::from_f64(initial) + change;
        Ok(finite(value.to_f64())?.abs())
    };
    let mut maximum = initial.abs().max(at(duration)?);
    if decay[0] != decay[1]
        && rates[0] != 0.0
        && rates[1] != 0.0
        && rates[0].is_sign_positive() != rates[1].is_sign_positive()
    {
        let (a, b) = (rates[0].abs(), rates[1].abs());
        let offset = (a - b) / b;
        let ratio = a / b;
        let logarithm = if offset.abs() <= 0.0625 {
            log1p_near_zero(offset)
        } else if ratio.is_finite() && ratio > 0.0 {
            det::ln(ratio)
        } else {
            det::ln(a) - det::ln(b)
        };
        let t = logarithm / (decay[0] - decay[1]);
        if t > 0.0 && t < duration {
            maximum = maximum.max(at(t)?);
        }
    }
    Ok(maximum)
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
struct ContactResponse {
    mobility: [[f64; 4]; 4],
    reaction: [[f64; 2]; 2],
    rolling_mobility: [[f64; 2]; 2],
    /// Normal point velocity and relative normal-axis spin per rolling impulse.
    unresolved_mobility: [[f64; 2]; 2],
}

fn response(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    binding: RigidContactBinding,
    axes: [Vec3; 2],
    normal: Vec3,
) -> Result<ContactResponse, RigidRollingError> {
    let mut factor = [[0.0; 12]; 6];
    for (slot, body, arm, sign) in [
        (0, binding.body_a, binding.arm_a_body, 1.0),
        (1, binding.body_b, binding.arm_b_body, -1.0),
    ] {
        let inverse_mass_root = properties[body].mass().sqrt().recip();
        let orientation = bodies[body].pose().orientation();
        for (axis, world_direction) in [axes[0], axes[1], normal].into_iter().enumerate() {
            let direction = orientation.rotate_world_to_body(world_direction);
            let translation = components(world_direction.scale(sign * inverse_mass_root));
            let force_rotation = components(
                properties[body].angular_mobility_factor(arm.cross(direction).scale(sign)),
            );
            let moment_rotation =
                components(properties[body].angular_mobility_factor(direction.scale(sign)));
            let force_row = if axis == 2 { 4 } else { axis };
            let moment_row = if axis == 2 { 5 } else { axis + 2 };
            for j in 0..3 {
                factor[force_row][6 * slot + j] = translation[j];
                factor[force_row][6 * slot + 3 + j] = force_rotation[j];
                factor[moment_row][6 * slot + 3 + j] = moment_rotation[j];
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
    let unresolved_mobility = [
        [
            dot(&factor[4], &projected[0])?,
            dot(&factor[4], &projected[1])?,
        ],
        [
            dot(&factor[5], &projected[0])?,
            dot(&factor[5], &projected[1])?,
        ],
    ];
    Ok(ContactResponse {
        mobility,
        reaction,
        rolling_mobility: rolling,
        unresolved_mobility,
    })
}

#[cfg(test)]
mod continuous_path_tests {
    use super::*;

    #[test]
    fn cancelling_decay_modes_reveal_interior_unsupported_motion() {
        // Integral(2 exp(-2t)-exp(-t)) = exp(-t)-exp(-2t).
        // Its initial and long-time limits vanish, but it reaches 1/4 at ln 2.
        let maximum =
            integrated_rate_max(0.0, [2.0, -1.0], [2.0, 1.0], 20.0, &CancelGate::new()).unwrap();
        assert!((maximum - 0.25).abs() < 2e-15);
        let force = decaying_norm_max([[1.0, -2.0], [-1.0, 2.0]], [2.0, 1.0], 20.0).unwrap();
        assert!((force - 5.0_f64.sqrt() / 4.0).abs() < 2e-15);
        let slow = 2.0_f64.next_down();
        let nearly_repeated =
            decaying_norm_max([[1e16, 0.0], [-1e16, 0.0]], [2.0, slow], 20.0).unwrap();
        // The relative-rate gap is one ULP, but its released cancellation is
        // macroscopically resolvable. The first-order peak has O(epsilon)
        // relative error here, below the stated floating comparison tolerance.
        let expected = 1e16 * ((2.0 - slow) / 2.0) * (-1.0_f64).exp();
        assert!((nearly_repeated - expected).abs() < 2e-15);
        let omitted =
            integrated_rate_max(0.0, [2.0, -slow], [2.0, slow], 20.0, &CancelGate::new()).unwrap();
        let expected = ((2.0 - slow) / 2.0) * (-1.0_f64).exp();
        assert!((omitted - expected).abs() < 2e-31);
        assert_eq!(
            integrated_rate_max(0.0, [2.0, -2.0], [2.0; 2], 20.0, &CancelGate::new()).unwrap(),
            0.0
        );
        let zero_mode =
            integrated_rate_max(0.0, [1.0, -1.0], [1e-20, 0.0], 1.0, &CancelGate::new()).unwrap();
        assert!((zero_mode - 0.5e-20).abs() < 1e-35);
    }
}
