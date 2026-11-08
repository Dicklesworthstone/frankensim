//! Original JKR attachment and pull-off on two finite rigid bodies.
//!
//! Actual masses, full inertia tensors and reaction arms supply the scalar
//! normal mobility. The existing JKR trajectory integrates its original force,
//! impulse and first time moment. These restore every native momentum and the
//! fixed-Jacobian COM/marker travel, including trajectories with a reversal or
//! snap. A snap changes only contact storage; its release remains an explicit
//! outgoing energy account, never an invented impulse or automatically heat.
//! Declared constant COM forces and body torques drive that same original
//! trajectory through their actual normal acceleration. Their work uses native
//! COM/angular travel, including the contact force's first time moment.
//!
//! This frictionless normal model supplies no tangential force or free spin
//! couple; the full inertial cross response remains unconstrained.
//! Geometry markers authenticate approach against an immutable reference datum;
//! forces act at their common midpoint reaction plane. Poses and orientations
//! stay held during this phase. The surrounding pose owner must consume travel
//! and refresh geometry before continuation. No finite-rotation trajectory,
//! contact discovery, tangential law, support or additional event law is inferred.

use crate::{DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3, Wrench};
use fs_contact::normal_patch::adhesive::{
    JkrState, SourceJkrSphereContact,
    motion::{
        JkrNormalMotion, JkrNormalMotionBudget, JkrNormalMotionError, JkrNormalMotionState,
        JkrNormalMotionStep,
    },
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_tribo::partial_slip::TangentFrame;

use super::super::hertz::{RigidHertzError, points, relative_point, scaled_response};
use super::super::{
    ResultantImpulse, RigidContactBinding, RigidFrictionError, add_components, magnitude, poll,
    rounded_components, vector,
};

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
    /// Native kinetic change + original contact work - actual external work [J].
    pub normal_work_residual_j: f64,
    /// Native kinetic + contact storage + snap release - external work [J].
    pub total_energy_residual_j: f64,
    /// Actual minus commanded impulse, in original A/B order [N s].
    pub body_linear_momentum_residual_ns: [Vec3; 2],
    /// Actual minus commanded body-frame angular impulse [N m s].
    pub body_angular_momentum_residual_nms: [Vec3; 2],
    /// Complete pair linear-momentum change [N s].
    pub linear_momentum_change_ns: Vec3,
    /// Complete pair orbital-plus-spin change at the held geometry [N m s].
    pub angular_momentum_change_nms: Vec3,
    /// Pair momentum change minus the actual external force impulses [N s].
    pub linear_momentum_residual_ns: Vec3,
    /// Pair orbital-plus-spin change minus actual external wrench impulse [N m s].
    pub angular_momentum_residual_nms: Vec3,
}

/// Original events and physical travel for the surrounding native pose owner.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrReport {
    /// Source trajectory, including branch events and unpartitioned snap release.
    /// Fresh `advance` phases start at zero; `continue_motion` retains the
    /// original accepted source/mass/policy/time/step history.
    pub motion: JkrNormalMotionStep,
    /// Actual normal `J M^-1 J^T`, retaining both full inertia tensors [1/kg].
    pub mobility_per_kg: f64,
    /// Equal/opposite physical point impulse on side A [N s].
    pub force_impulse_on_a_world_ns: Vec3,
    /// Integral t times the force on A, with time local to this phase [N s²].
    pub force_time_moment_on_a_world_ns2: Vec3,
    /// Held world-frame COM forces and body-frame COM torques, in A/B order.
    pub external_wrenches: [Wrench; 2],
    /// Held scalar force derived from actual native normal acceleration and
    /// mobility; positive closes the gap [N]. It is not an imposed body load.
    pub applied_normal_force_n: f64,
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
    /// Combined impulses contracted with actual midpoint COM/angular rates [J].
    pub impulse_work_j: f64,
    /// Actual force·COM travel + body torque·angular travel [J]. Contact timing
    /// enters through its first time moment, including reversals and snaps.
    pub external_work_j: f64,
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

    /// Initialize a resumable original trajectory from the actual captured
    /// native geometry and normal rate. This creates no motion or branch event.
    /// Retain each accepted `report.motion.state()` for later continuation.
    pub fn initialize_motion(
        &self,
        previous: JkrState,
    ) -> Result<JkrNormalMotionState, RigidJkrError> {
        let geometry_residual = finite(previous.indentation_m() - self.approach_m)?;
        if geometry_residual.abs() > self.controls.geometry_tolerance_m {
            return Err(RigidJkrError::GeometryMismatch {
                residual_m: geometry_residual,
            });
        }
        let velocity = finite(relative_point(self.reaction_before)?.dot(self.normal))?;
        Ok(self.motion.initialize(previous, velocity)?)
    }

    /// Continue the original accepted clock and event history after the pose
    /// owner consumes travel and captures the refreshed native geometry.
    /// Actual approach/rate must agree within the original dimensional gates.
    /// Source, effective mass and numerical policy must match exactly; changed
    /// mobility is a refusal, not an implicit velocity/history reinitialization.
    pub fn continue_motion(
        &self,
        previous: JkrNormalMotionState,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<RigidJkrCandidate, RigidJkrError> {
        self.continue_motion_with_wrenches(previous, [Wrench::ZERO; 2], dt_s, gate)
    }

    /// Continue the same original trajectory with the new phase's held body
    /// wrenches. Loads may change between phases without changing the accepted
    /// source/mass/policy/time token or replaying a previous attachment/release.
    pub fn continue_motion_with_wrenches(
        &self,
        previous: JkrNormalMotionState,
        wrenches: [Wrench; 2],
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<RigidJkrCandidate, RigidJkrError> {
        self.advance_state(previous, wrenches, dt_s, gate)
    }

    /// Advance the explicit original branch at the ACTUAL initial body rate.
    /// Source tokens retain attachment history; body kinematics, not a cached
    /// effective mass or external scalar velocity, initialize every phase.
    pub fn advance(
        &self,
        previous: JkrState,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<RigidJkrCandidate, RigidJkrError> {
        self.advance_with_wrenches(previous, [Wrench::ZERO; 2], dt_s, gate)
    }

    /// Advance with held forces at the actual COMs and held body-frame torques,
    /// in original A/B order. Full native inertia determines their normal
    /// acceleration; unconstrained tangential and angular response is retained.
    /// Contact and external impulses are accumulated before either body changes.
    pub fn advance_with_wrenches(
        &self,
        previous: JkrState,
        wrenches: [Wrench; 2],
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<RigidJkrCandidate, RigidJkrError> {
        poll(gate)?;
        let state = self.initialize_motion(previous)?;
        self.advance_state(state, wrenches, dt_s, gate)
    }

    #[allow(clippy::too_many_lines)]
    fn advance_state(
        &self,
        state: JkrNormalMotionState,
        wrenches: [Wrench; 2],
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<RigidJkrCandidate, RigidJkrError> {
        poll(gate)?;
        let before_response = self.motion.response(state)?;
        let geometry_residual = finite(state.contact_state().indentation_m() - self.approach_m)?;
        if geometry_residual.abs() > self.controls.geometry_tolerance_m {
            return Err(RigidJkrError::GeometryMismatch {
                residual_m: geometry_residual,
            });
        }
        let initial = relative_point(self.reaction_before)?;
        let velocity_before = finite(initial.dot(self.normal))?;
        if (state.velocity_m_s() - velocity_before).abs() > self.controls.velocity_tolerance_m_s {
            return Err(RigidJkrError::Invalid(
                "original motion rate disagrees with actual native normal rate",
            ));
        }
        let mut free_acceleration = [Vec3::ZERO; 2];
        for i in 0..2 {
            wrenches[i].validate()?;
            let angular_acceleration =
                self.properties[i].angular_velocity_body_checked(wrenches[i].torque_body)?;
            free_acceleration[i] = wrenches[i]
                .force_world
                .scale(self.properties[i].mass().recip())
                .add(
                    self.before[i]
                        .pose()
                        .orientation()
                        .rotate_body_to_world(angular_acceleration.cross(self.reaction_arms[i])),
                );
        }
        let applied_normal_force = finite(
            free_acceleration[0]
                .sub(free_acceleration[1])
                .dot(self.normal)
                / self.mobility,
        )?;
        let step = self
            .motion
            .advance(state, applied_normal_force, dt_s, gate)?;
        let impulse = self.normal.scale(-step.normal_impulse_n_s);
        let time_moment = self.normal.scale(-step.normal_time_moment_n_s2);
        let mut after = self.before;
        let mut kinetic_change = Dd::from_f64(0.0);
        let mut work = Dd::from_f64(0.0);
        let mut linear_change = [Dd::from_f64(0.0); 3];
        let mut angular_change = [Dd::from_f64(0.0); 3];
        let mut external_linear = [Dd::from_f64(0.0); 3];
        let mut external_angular = [Dd::from_f64(0.0); 3];
        let mut linear_residual = [Vec3::ZERO; 2];
        let mut angular_residual = [Vec3::ZERO; 2];
        for (i, sign) in [(0, 1.0), (1, -1.0)] {
            poll(gate)?;
            let body = self.before[i];
            let orientation = body.pose().orientation();
            let external_impulse = wrenches[i].force_world.scale(dt_s);
            let external_torque = wrenches[i].torque_body.scale(dt_s);
            let mut resultant = ResultantImpulse::default();
            resultant.add(body, self.reaction_arms[i], impulse.scale(sign), Vec3::ZERO);
            add_components(&mut resultant.linear, external_impulse);
            // The real wrench torque is already in its declared body frame.
            add_components(&mut resultant.angular, external_torque);
            let receipt = body.apply_momentum_impulse(
                self.properties[i],
                rounded_components(resultant.linear),
                rounded_components(resultant.angular),
            )?;
            after[i] = receipt.state_after;
            linear_residual[i] = receipt.linear_residual_world;
            angular_residual[i] = receipt.angular_residual_body;
            kinetic_change = kinetic_change + Dd::from_f64(receipt.work.kinetic_energy_change_j);
            work = work + Dd::from_f64(receipt.work.impulse_work_j);
            add_components(&mut linear_change, receipt.linear_increment_world);
            add_components(&mut external_linear, external_impulse);
            add_components(
                &mut angular_change,
                body.pose()
                    .position_world()
                    .cross(receipt.linear_increment_world),
            );
            add_components(
                &mut angular_change,
                orientation.rotate_body_to_world(receipt.angular_increment_body),
            );
            add_components(
                &mut external_angular,
                body.pose().position_world().cross(external_impulse),
            );
            add_components(
                &mut external_angular,
                orientation.rotate_body_to_world(external_torque),
            );
        }
        let kinetic_change = finite(kinetic_change.to_f64())?;
        let work = finite(work.to_f64())?;
        let reaction_after = points(after, self.properties, self.reaction_arms)?;
        let marker_arms = [
            self.contact.binding.arm_a_body,
            self.contact.binding.arm_b_body,
        ];
        let markers_after = points(after, self.properties, marker_arms)?;
        let endpoint = relative_point(reaction_after)?;
        let velocity_after = finite(endpoint.dot(self.normal))?;
        // Integral (h-t) F(t) dt, not h/2 times the total impulse: branch
        // events and varying forces require the original first time moment.
        let displacement_impulse = impulse.scale(dt_s).sub(time_moment);
        let mut body_travel = [Vec3::ZERO; 2];
        let mut angular_travel = [Vec3::ZERO; 2];
        let mut marker_travel = [Vec3::ZERO; 2];
        let mut reaction_travel = [Vec3::ZERO; 2];
        let mut external_work = Dd::from_f64(0.0);
        for (i, sign) in [(0, 1.0), (1, -1.0)] {
            poll(gate)?;
            let orientation = self.before[i].pose().orientation();
            let integrated_contact = displacement_impulse.scale(sign);
            let integrated_force =
                integrated_contact.add(wrenches[i].force_world.scale(dt_s).scale(0.5 * dt_s));
            body_travel[i] = self.reaction_before[i]
                .center_of_mass_velocity_world
                .scale(dt_s)
                .add(integrated_force.scale(self.properties[i].mass().recip()));
            // Apply the same full inverse inertia to the torque time integral.
            let integrated_torque = self.reaction_arms[i]
                .cross(orientation.rotate_world_to_body(integrated_contact))
                .add(wrenches[i].torque_body.scale(dt_s).scale(0.5 * dt_s));
            angular_travel[i] = self.reaction_before[i]
                .angular_velocity_body
                .scale(dt_s)
                .add(self.properties[i].angular_velocity_body_checked(integrated_torque)?);
            marker_travel[i] = body_travel[i]
                .add(orientation.rotate_body_to_world(angular_travel[i].cross(marker_arms[i])));
            reaction_travel[i] = body_travel[i].add(
                orientation.rotate_body_to_world(angular_travel[i].cross(self.reaction_arms[i])),
            );
            external_work = external_work
                + Dd::from_f64(wrenches[i].force_world.dot(body_travel[i]))
                + Dd::from_f64(wrenches[i].torque_body.dot(angular_travel[i]));
        }
        let external_work = finite(external_work.to_f64())?;
        let travel = finite(reaction_travel[0].sub(reaction_travel[1]).dot(self.normal))?;
        let marker_travel_normal = finite(marker_travel[0].sub(marker_travel[1]).dot(self.normal))?;
        let marker_velocity_residual = (relative_point(self.markers_before)?.dot(self.normal)
            - velocity_before)
            .abs()
            .max((relative_point(markers_after)?.dot(self.normal) - velocity_after).abs());
        let storage_change =
            finite(step.contact.reversible_energy_j - before_response.reversible_energy_j)?;
        let linear_balance = rounded_components(core::array::from_fn(|i| {
            linear_change[i] - external_linear[i]
        }));
        let angular_balance = rounded_components(core::array::from_fn(|i| {
            angular_change[i] - external_angular[i]
        }));
        let audit = RigidJkrAudit {
            normal_displacement_residual_m: travel - step.normal_increment_m,
            marker_displacement_residual_m: marker_travel_normal - travel,
            normal_velocity_residual_m_s: velocity_after - step.state().velocity_m_s(),
            marker_velocity_residual_m_s: marker_velocity_residual,
            work_energy_residual_j: kinetic_change - work,
            normal_work_residual_j: kinetic_change + step.normal_work_j - external_work,
            total_energy_residual_j: kinetic_change + storage_change + step.released_energy_j
                - external_work,
            body_linear_momentum_residual_ns: linear_residual,
            body_angular_momentum_residual_nms: angular_residual,
            linear_momentum_change_ns: rounded_components(linear_change),
            angular_momentum_change_nms: rounded_components(angular_change),
            linear_momentum_residual_ns: linear_balance,
            angular_momentum_residual_nms: angular_balance,
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
                .chain([linear_balance])
                .any(|v| !v.is_finite() || magnitude(v) > c.linear_momentum_tolerance_ns)
            || angular_residual
                .into_iter()
                .chain([angular_balance])
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
                external_wrenches: wrenches,
                applied_normal_force_n: applied_normal_force,
                relative_normal_velocity_m_s: [velocity_before, velocity_after],
                relative_normal_displacement_m: travel,
                body_displacement_world_m: body_travel,
                angular_displacement_body_rad: angular_travel,
                marker_displacement_world_m: marker_travel,
                reaction_displacement_world_m: reaction_travel,
                kinetic_energy_change_j: kinetic_change,
                impulse_work_j: work,
                external_work_j: external_work,
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

/// Stage a trajectory under held actual body wrenches in original A/B order.
#[allow(clippy::too_many_arguments)]
pub fn prepare_rigid_jkr_contact_with_wrenches(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidJkrContact<'_>,
    previous: JkrState,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidJkrControls,
    gate: &CancelGate,
) -> Result<RigidJkrCandidate, RigidJkrError> {
    PreparedRigidJkrPair::new(bodies, properties, contact, controls, gate)?
        .advance_with_wrenches(previous, wrenches, dt_s, gate)
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
    apply_rigid_jkr_contact_with_wrenches(
        bodies,
        properties,
        contact,
        history,
        [Wrench::ZERO; 2],
        dt_s,
        controls,
        gate,
    )
}

/// Atomically publish loaded body momenta and source history after every check.
/// Wrenches act at the actual COMs, in original A/B order; poses remain held.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_jkr_contact_with_wrenches(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidJkrContact<'_>,
    history: &mut JkrState,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidJkrControls,
    gate: &CancelGate,
) -> Result<RigidJkrReport, RigidJkrError> {
    let candidate = prepare_rigid_jkr_contact_with_wrenches(
        bodies, properties, contact, *history, wrenches, dt_s, controls, gate,
    )?;
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
