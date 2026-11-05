//! Original JKR attachment and pull-off on two finite rigid bodies.
//!
//! Actual masses, full inertia tensors and reaction arms supply the scalar
//! normal mobility. The existing JKR trajectory integrates its original force,
//! impulse and first time moment. These restore every native momentum and the
//! fixed-Jacobian COM/marker travel, including trajectories with a reversal or
//! snap. A snap changes only contact storage; its release remains an explicit
//! outgoing energy account, never an invented impulse or automatically heat.
//!
//! This frictionless normal model supplies no tangential force or free spin
//! couple; the full inertial cross response remains unconstrained.
//! Geometry markers authenticate approach against an immutable reference datum;
//! forces act at their common midpoint reaction plane. Poses and orientations
//! stay held during this phase. The surrounding pose owner must consume travel
//! and refresh geometry before continuation. No finite-rotation trajectory,
//! contact discovery, tangential law, support or additional event law is inferred.

use crate::{DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3};
use fs_contact::normal_patch::adhesive::{
    JkrState, SourceJkrSphereContact,
    motion::{JkrNormalMotion, JkrNormalMotionBudget, JkrNormalMotionError, JkrNormalMotionStep},
};
use fs_exec::CancelGate;
use fs_tribo::partial_slip::TangentFrame;

use super::super::hertz::{
    RigidHertzError, points, prepare_pair_wrench, relative_point, scaled_response,
};
use super::super::{RigidContactBinding, RigidFrictionError, magnitude, poll, vector};

/// Original normal source and physical pair with COM-relative geometry markers.
#[derive(Clone, Copy, Debug)]
pub struct RigidJkrContact<'a> {
    /// Original body indices and body-frame geometry-marker arms.
    pub binding: RigidContactBinding,
    /// Normal points from A toward B; positive A-minus-B motion closes the gap.
    pub frame: TangentFrame,
    /// Approach at coincident markers [m], fixed by geometry, not by history.
    pub reference_approach_m: f64,
    /// Original source-bound large-Tabor frictionless contact.
    pub normal: &'a SourceJkrSphereContact,
}

/// Original motion budget plus independent dimensional native-body checks.
#[derive(Clone, Copy, Debug)]
pub struct RigidJkrControls {
    /// Original scalar motion, event and constitutive-domain budget.
    pub motion: JkrNormalMotionBudget,
    /// Maximum size of the original body array.
    pub maximum_bodies: usize,
    /// Initial source/marker approach and marker collinearity gate [m].
    pub geometry_tolerance_m: f64,
    /// Native first-moment travel minus source trajectory travel [m].
    pub displacement_tolerance_m: f64,
    /// Native normal-rate minus source trajectory rate [m/s].
    pub velocity_tolerance_m_s: f64,
    /// Individual and total pair linear-impulse gate [N s].
    pub linear_momentum_tolerance_ns: f64,
    /// Individual and total pair angular-impulse gate [N m s].
    pub angular_momentum_tolerance_nms: f64,
    /// Native work, original normal work and total energy gates [J].
    pub energy_tolerance_j: f64,
}

/// Independent checks on actual body endpoints and the original normal source.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrAudit {
    /// Native reaction-point travel minus original source travel [m].
    pub normal_displacement_residual_m: f64,
    /// Geometry-marker travel minus reaction-point normal travel [m].
    pub marker_displacement_residual_m: f64,
    /// Actual endpoint normal rate minus the original trajectory endpoint [m/s].
    pub normal_velocity_residual_m_s: f64,
    /// Largest marker-versus-reaction normal-rate difference [m/s].
    pub marker_velocity_residual_m_s: f64,
    /// Native kinetic change minus native impulse work [J].
    pub work_energy_residual_j: f64,
    /// Native kinetic change plus original integrated contact work [J].
    pub normal_work_residual_j: f64,
    /// Native kinetic + contact storage change + explicit snap release [J].
    pub total_energy_residual_j: f64,
    /// Actual minus commanded impulse, in original A/B order [N s].
    pub body_linear_momentum_residual_ns: [Vec3; 2],
    /// Actual minus commanded body-frame angular impulse [N m s].
    pub body_angular_momentum_residual_nms: [Vec3; 2],
    /// Complete pair linear-momentum change [N s].
    pub linear_momentum_change_ns: Vec3,
    /// Complete pair orbital-plus-spin change at the held geometry [N m s].
    pub angular_momentum_change_nms: Vec3,
}

/// Original events and physical travel for the surrounding native pose owner.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrReport {
    /// Source trajectory, including branch events and unpartitioned snap release.
    /// Its clock starts at zero for this captured native phase.
    pub motion: JkrNormalMotionStep,
    /// Actual normal `J M^-1 J^T`, retaining both full inertia tensors [1/kg].
    pub mobility_per_kg: f64,
    /// Equal/opposite physical point impulse on side A [N s].
    pub force_impulse_on_a_world_ns: Vec3,
    /// Integral t times the force on A, with time local to this phase [N s²].
    pub force_time_moment_on_a_world_ns2: Vec3,
    /// Actual initial/final A-minus-B normal rates [m/s].
    pub relative_normal_velocity_m_s: [f64; 2],
    /// Actual normal travel reconstructed from the impulse first moment [m].
    pub relative_normal_displacement_m: f64,
    /// Integrated COM velocity in original A/B order [m].
    pub body_displacement_world_m: [Vec3; 2],
    /// Integrated angular velocity in the held body frame [rad], not a finite
    /// rotation or an instruction to overwrite the actual body orientation.
    pub angular_displacement_body_rad: [Vec3; 2],
    /// Integrated geometry-marker velocity in the held Jacobian [m].
    pub marker_displacement_world_m: [Vec3; 2],
    /// Integrated common-reaction-point velocity in the held Jacobian [m].
    pub reaction_displacement_world_m: [Vec3; 2],
    /// Compensated kinetic change of the original two finite bodies [J].
    pub kinetic_energy_change_j: f64,
    /// Actual impulse contracted with actual midpoint point velocity [J].
    pub impulse_work_j: f64,
    /// Original endpoint elastic-plus-surface potential change [J].
    pub normal_storage_change_j: f64,
    /// Independent dimensional checks, without source-authority promotion.
    pub audit: RigidJkrAudit,
}

/// Staged endpoints; failed preparation never changes bodies or contact history.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrCandidate {
    /// Original A/B order; poses remain unchanged.
    pub bodies_after: [RigidBodyState; 2],
    /// Original source-authenticated successor branch and approach.
    pub state_after: JkrState,
    /// Original events, native travel and physical balances.
    pub report: RigidJkrReport,
}

/// Input, source, native mechanics or independently audited balance refusal.
#[derive(Debug)]
pub enum RigidJkrError {
    /// Invalid or unrepresentable native input.
    Invalid(&'static str),
    /// The original scalar trajectory refused its source, domain or budget.
    Motion(JkrNormalMotionError),
    /// Actual kinematics or momentum application refused.
    Mechanics(RigidHertzError),
    /// Accepted source approach disagrees with actual geometry-marker approach.
    GeometryMismatch {
        /// Source minus geometry approach [m].
        residual_m: f64,
    },
    /// Actual native endpoints failed one or more original dimensional gates.
    Audit(Box<RigidJkrAudit>),
}
impl From<JkrNormalMotionError> for RigidJkrError {
    fn from(value: JkrNormalMotionError) -> Self {
        Self::Motion(value)
    }
}
impl From<RigidHertzError> for RigidJkrError {
    fn from(value: RigidHertzError) -> Self {
        Self::Mechanics(value)
    }
}
impl From<DynamicsError> for RigidJkrError {
    fn from(value: DynamicsError) -> Self {
        Self::Mechanics(value.into())
    }
}
impl From<RigidFrictionError> for RigidJkrError {
    fn from(value: RigidFrictionError) -> Self {
        Self::Mechanics(value.into())
    }
}
impl core::fmt::Display for RigidJkrError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native JKR normal trajectory: {self:?}")
    }
}
impl std::error::Error for RigidJkrError {}

/// Captured native geometry and original constitutive trajectory. Repeated
/// trials always start at these same physical body endpoints and never mutate
/// accepted history. No caller effective mass replaces the actual mobility.
pub struct PreparedRigidJkrPair<'a> {
    contact: RigidJkrContact<'a>,
    before: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    reaction_arms: [Vec3; 2],
    normal: Vec3,
    markers_before: [PointKinematics; 2],
    reaction_before: [PointKinematics; 2],
    approach_m: f64,
    mobility: f64,
    motion: JkrNormalMotion<'a>,
    controls: RigidJkrControls,
}

impl<'a> PreparedRigidJkrPair<'a> {
    /// Capture actual marker geometry, full native normal response and bounds.
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contact: RigidJkrContact<'a>,
        controls: RigidJkrControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidJkrError> {
        poll(gate)?;
        let binding = contact.binding;
        if bodies.len() != properties.len()
            || bodies.len() > controls.maximum_bodies
            || binding.body_a >= bodies.len()
            || binding.body_b >= bodies.len()
            || binding.body_a == binding.body_b
            || !contact.reference_approach_m.is_finite()
            || [
                controls.geometry_tolerance_m,
                controls.displacement_tolerance_m,
                controls.velocity_tolerance_m_s,
                controls.linear_momentum_tolerance_ns,
                controls.angular_momentum_tolerance_nms,
                controls.energy_tolerance_j,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(RigidJkrError::Invalid(
                "bounded distinct bodies and finite nonnegative controls required",
            ));
        }
        let before = [bodies[binding.body_a], bodies[binding.body_b]];
        let properties = [properties[binding.body_a], properties[binding.body_b]];
        let marker_arms = [binding.arm_a_body, binding.arm_b_body];
        let normal = vector(contact.frame.normal());
        let markers_before = points(before, properties, marker_arms)?;
        let separation = markers_before[1]
            .point_world
            .sub(markers_before[0].point_world);
        let separation_normal = finite(separation.dot(normal))?;
        if magnitude(separation.sub(normal.scale(separation_normal)))
            > controls.geometry_tolerance_m
        {
            return Err(RigidJkrError::Invalid(
                "geometry markers must lie on the declared normal line",
            ));
        }
        let approach_m = finite(contact.reference_approach_m - separation_normal)?;
        let reaction_arms = core::array::from_fn(|i| {
            let sign = if i == 0 { 0.5 } else { -0.5 };
            marker_arms[i].add(
                before[i]
                    .pose()
                    .orientation()
                    .rotate_world_to_body(separation.scale(sign)),
            )
        });
        let reaction_before = points(before, properties, reaction_arms)?;
        if magnitude(
            reaction_before[0]
                .point_world
                .sub(reaction_before[1].point_world),
        ) > controls.geometry_tolerance_m
        {
            return Err(RigidJkrError::Invalid(
                "unrepresentable common reaction plane",
            ));
        }
        // For N=1 this existing Gram helper has one point-force row, no couple.
        let mobility = scaled_response(before, properties, reaction_arms, [normal], [0.0])?[0][0];
        let motion = JkrNormalMotion::new(contact.normal, mobility.recip(), controls.motion)?;
        poll(gate)?;
        Ok(Self {
            contact,
            before,
            properties,
            reaction_arms,
            normal,
            markers_before,
            reaction_before,
            approach_m,
            mobility,
            motion,
            controls,
        })
    }

    /// Original source, marker binding and approach datum.
    pub const fn contact(&self) -> RigidJkrContact<'a> {
        self.contact
    }
    /// Actual signed approach from captured geometry [m].
    pub const fn approach_m(&self) -> f64 {
        self.approach_m
    }
    /// Actual scalar normal mobility, without a caller effective mass [1/kg].
    pub const fn mobility_per_kg(&self) -> f64 {
        self.mobility
    }
    /// COM-relative arms of the common reaction plane, in A/B body frames.
    pub const fn reaction_arms_body(&self) -> [Vec3; 2] {
        self.reaction_arms
    }

    /// Advance the explicit original branch at the ACTUAL initial body rate.
    /// Source tokens retain attachment history; body kinematics, not a cached
    /// effective mass or external scalar velocity, initialize every phase.
    #[allow(clippy::too_many_lines)]
    pub fn advance(
        &self,
        previous: JkrState,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<RigidJkrCandidate, RigidJkrError> {
        poll(gate)?;
        let geometry_residual = finite(previous.indentation_m() - self.approach_m)?;
        if geometry_residual.abs() > self.controls.geometry_tolerance_m {
            return Err(RigidJkrError::GeometryMismatch {
                residual_m: geometry_residual,
            });
        }
        let initial = relative_point(self.reaction_before);
        let velocity_before = finite(initial.dot(self.normal))?;
        let state = self.motion.initialize(previous, velocity_before)?;
        let before_response = self.motion.response(state)?;
        let step = self.motion.advance(state, 0.0, dt_s, gate)?;
        let impulse = self.normal.scale(-step.normal_impulse_n_s);
        let time_moment = self.normal.scale(-step.normal_time_moment_n_s2);
        let phase = prepare_pair_wrench(
            self.before,
            self.properties,
            self.reaction_arms,
            impulse,
            Vec3::ZERO,
            gate,
        )?;
        let after = phase.bodies_after;
        let reaction_after = points(after, self.properties, self.reaction_arms)?;
        let marker_arms = [
            self.contact.binding.arm_a_body,
            self.contact.binding.arm_b_body,
        ];
        let markers_after = points(after, self.properties, marker_arms)?;
        let endpoint = relative_point(reaction_after);
        let velocity_after = finite(endpoint.dot(self.normal))?;
        // Integral (h-t) F(t) dt, not h/2 times the total impulse: branch
        // events and varying forces require the original first time moment.
        let displacement_impulse = impulse.scale(dt_s).sub(time_moment);
        let mut body_travel = [Vec3::ZERO; 2];
        let mut angular_travel = [Vec3::ZERO; 2];
        let mut marker_travel = [Vec3::ZERO; 2];
        let mut reaction_travel = [Vec3::ZERO; 2];
        let mut linear_residual = [Vec3::ZERO; 2];
        let mut angular_residual = [Vec3::ZERO; 2];
        for (i, sign) in [(0, 1.0), (1, -1.0)] {
            poll(gate)?;
            let orientation = self.before[i].pose().orientation();
            let expected = impulse.scale(sign);
            linear_residual[i] = after[i]
                .linear_momentum_world()
                .sub(self.before[i].linear_momentum_world())
                .sub(expected);
            angular_residual[i] = after[i]
                .angular_momentum_body()
                .sub(self.before[i].angular_momentum_body())
                .sub(self.reaction_arms[i].cross(orientation.rotate_world_to_body(expected)));
            let integrated_force = displacement_impulse.scale(sign);
            body_travel[i] = self.reaction_before[i]
                .center_of_mass_velocity_world
                .scale(dt_s)
                .add(integrated_force.scale(self.properties[i].mass().recip()));
            // Apply the same full inverse inertia to the torque time integral.
            let integrated_torque =
                self.reaction_arms[i].cross(orientation.rotate_world_to_body(integrated_force));
            angular_travel[i] = self.reaction_before[i]
                .angular_velocity_body
                .scale(dt_s)
                .add(self.properties[i].angular_velocity_body_checked(integrated_torque)?);
            marker_travel[i] = body_travel[i]
                .add(orientation.rotate_body_to_world(angular_travel[i].cross(marker_arms[i])));
            reaction_travel[i] = body_travel[i].add(
                orientation.rotate_body_to_world(angular_travel[i].cross(self.reaction_arms[i])),
            );
        }
        let travel = finite(reaction_travel[0].sub(reaction_travel[1]).dot(self.normal))?;
        let marker_travel_normal = finite(marker_travel[0].sub(marker_travel[1]).dot(self.normal))?;
        let marker_velocity_residual = (relative_point(self.markers_before).dot(self.normal)
            - velocity_before)
            .abs()
            .max((relative_point(markers_after).dot(self.normal) - velocity_after).abs());
        let work = finite(impulse.dot(initial.scale(0.5).add(endpoint.scale(0.5))))?;
        let storage_change =
            finite(step.contact.reversible_energy_j - before_response.reversible_energy_j)?;
        let audit = RigidJkrAudit {
            normal_displacement_residual_m: travel - step.normal_increment_m,
            marker_displacement_residual_m: marker_travel_normal - travel,
            normal_velocity_residual_m_s: velocity_after - step.state().velocity_m_s(),
            marker_velocity_residual_m_s: marker_velocity_residual,
            work_energy_residual_j: phase.kinetic_change_j - work,
            normal_work_residual_j: phase.kinetic_change_j + step.normal_work_j,
            total_energy_residual_j: phase.kinetic_change_j
                + storage_change
                + step.released_energy_j,
            body_linear_momentum_residual_ns: linear_residual,
            body_angular_momentum_residual_nms: angular_residual,
            linear_momentum_change_ns: phase.linear_change_ns,
            angular_momentum_change_nms: phase.angular_change_nms,
        };
        let c = self.controls;
        if [
            audit.normal_displacement_residual_m,
            audit.marker_displacement_residual_m,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v.abs() > c.displacement_tolerance_m)
            || [audit.normal_velocity_residual_m_s, marker_velocity_residual]
                .into_iter()
                .any(|v| !v.is_finite() || v.abs() > c.velocity_tolerance_m_s)
            || [
                audit.work_energy_residual_j,
                audit.normal_work_residual_j,
                audit.total_energy_residual_j,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v.abs() > c.energy_tolerance_j)
            || linear_residual
                .into_iter()
                .chain([phase.linear_change_ns])
                .any(|v| !v.is_finite() || magnitude(v) > c.linear_momentum_tolerance_ns)
            || angular_residual
                .into_iter()
                .chain([phase.angular_change_nms])
                .any(|v| !v.is_finite() || magnitude(v) > c.angular_momentum_tolerance_nms)
            || body_travel
                .into_iter()
                .chain(angular_travel)
                .chain(marker_travel)
                .chain(reaction_travel)
                .any(|v| !v.is_finite())
        {
            return Err(RigidJkrError::Audit(Box::new(audit)));
        }
        poll(gate)?;
        Ok(RigidJkrCandidate {
            bodies_after: after,
            state_after: step.state().contact_state(),
            report: RigidJkrReport {
                motion: step,
                mobility_per_kg: self.mobility,
                force_impulse_on_a_world_ns: impulse,
                force_time_moment_on_a_world_ns2: time_moment,
                relative_normal_velocity_m_s: [velocity_before, velocity_after],
                relative_normal_displacement_m: travel,
                body_displacement_world_m: body_travel,
                angular_displacement_body_rad: angular_travel,
                marker_displacement_world_m: marker_travel,
                reaction_displacement_world_m: reaction_travel,
                kinetic_energy_change_j: phase.kinetic_change_j,
                impulse_work_j: work,
                normal_storage_change_j: storage_change,
                audit,
            },
        })
    }
}

/// Stage one complete trajectory without changing accepted native state.
#[allow(clippy::too_many_arguments)]
pub fn prepare_rigid_jkr_contact(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidJkrContact<'_>,
    previous: JkrState,
    dt_s: f64,
    controls: RigidJkrControls,
    gate: &CancelGate,
) -> Result<RigidJkrCandidate, RigidJkrError> {
    PreparedRigidJkrPair::new(bodies, properties, contact, controls, gate)?
        .advance(previous, dt_s, gate)
}

/// Publish both momenta and the original normal history only after admission.
/// Poses remain unchanged; the returned travel belongs to the surrounding pose
/// owner. An advanced source against stale marker geometry refuses continuation.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_jkr_contact(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidJkrContact<'_>,
    history: &mut JkrState,
    dt_s: f64,
    controls: RigidJkrControls,
    gate: &CancelGate,
) -> Result<RigidJkrReport, RigidJkrError> {
    let candidate =
        prepare_rigid_jkr_contact(bodies, properties, contact, *history, dt_s, controls, gate)?;
    poll(gate)?;
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    *history = candidate.state_after;
    Ok(candidate.report)
}

fn finite(value: f64) -> Result<f64, RigidJkrError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RigidJkrError::Invalid(
            "unrepresentable native JKR mechanics",
        ))
    }
}
