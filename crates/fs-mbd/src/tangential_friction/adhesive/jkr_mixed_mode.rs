//! Native paired-body recoil with coupled full-stick JKR normal/shear storage.
//!
//! The original mixed-mode source equilibrates the actual disk and supplies
//! every effort, path tangent, work and sourced static-strength check. All nine
//! entries of the actual three-point-force mobility enter its joint solve.
//! Shear therefore changes normal force and area; normal and shear recoil are
//! not independent effective masses. Elastic unloading can return energy to
//! the bodies without being mislabeled friction heat.
//!
//! Collinear geometry markers and an immutable datum authenticate normal
//! approach. The two shear components are retained CONSTITUTIVE elastic shifts,
//! driven by actual velocities at the common midpoint reaction plane. They are
//! not raw tangential gaps between separated markers. This distinction keeps
//! all three forces work-conjugate while cancelling internal angular impulse.
//!
//! One phase holds the Jacobian, frame and mean contact effort. Its exact
//! quadratic discrete trajectory is screened by the original source owner;
//! the native poses stay held and their travel belongs to the surrounding pose
//! owner. History retains source, mobility, policy, clock, marker association
//! and tangent frame. A changed frame or mobility needs its own history
//! transfer. Torsion, sliding and radius instability refuse: this owner invents
//! no spin law, snap continuation, heat, contact migration or finite rotation.

/// Compatible endpoint poses with retained original native/source history.
pub mod geometry;
/// Original-source temperature continuation at unchanged native state.
pub mod temperature;

use crate::{DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3, Wrench};
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::adhesive::mixed_mode::{
    SourceJkrMixedMode,
    motion::{
        JkrMixedModeMotion, JkrMixedModeMotionBudget, JkrMixedModeMotionError,
        JkrMixedModeMotionState, JkrMixedModeMotionStep,
    },
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_tribo::partial_slip::TangentFrame;

use super::super::hertz::{RigidHertzError, point_response, points, relative_point, relative_spin};
use super::super::{
    ResultantImpulse, RigidContactBinding, RigidFrictionError, add_components, magnitude, poll,
    rounded_components, vector,
};

/// Original ordered interface and native material-point association.
#[derive(Clone, Copy, Debug)]
pub struct RigidJkrMixedModeContact<'a> {
    pub binding: RigidContactBinding,
    /// Positive A-minus-B rates increase [approach, elastic shear 1, shear 2].
    pub frame: TangentFrame,
    pub reference_approach_m: f64,
    pub source: &'a SourceJkrMixedMode,
    /// The same ordered source's measured/declared static junction strength.
    pub shear: &'a SourceAdhesiveShear,
}

/// Source work/solve budgets plus dimensional checks on actual native bodies.
#[derive(Clone, Copy, Debug)]
pub struct RigidJkrMixedModeControls {
    pub motion: JkrMixedModeMotionBudget,
    pub maximum_bodies: usize,
    pub geometry_tolerance_m: f64,
    pub displacement_tolerance_m: f64,
    pub velocity_tolerance_m_s: f64,
    /// Maximum unrepresented relative spin about the original normal [rad/s].
    pub angular_velocity_tolerance_rad_s: f64,
    pub linear_momentum_tolerance_ns: f64,
    pub angular_momentum_tolerance_nms: f64,
    pub energy_tolerance_j: f64,
}

/// Sealed source chronology plus the frame and material markers of elastic shear.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidJkrMixedModeState {
    motion: JkrMixedModeMotionState,
    frame: TangentFrame,
    marker_arms: [Vec3; 2],
    reference_approach_m: f64,
    endpoint: NativeEndpoint,
}

/// Accepted momenta remain at the phase geometry until a pose owner supplies
/// compatible endpoints. Pending COM travel is consumed only by that handoff.
#[derive(Clone, Copy, Debug, PartialEq)]
struct NativeEndpoint {
    body_indices: [usize; 2],
    properties: [MassProperties; 2],
    bodies: [RigidBodyState; 2],
    pending_com_travel_m: [Vec3; 2],
    // A sum of angular vectors is not a finite-rotation history. Retain whether
    // any pending phase moved so cancellation cannot erase that distinction.
    pending_angular_motion: [bool; 2],
}

impl RigidJkrMixedModeState {
    pub const fn motion_state(self) -> JkrMixedModeMotionState {
        self.motion
    }
    pub const fn position_m(self) -> [f64; 3] {
        self.motion.position_m()
    }
    pub const fn velocity_m_s(self) -> [f64; 3] {
        self.motion.velocity_m_s()
    }
    pub const fn time_s(self) -> f64 {
        self.motion.time_s()
    }
    pub const fn steps(self) -> u64 {
        self.motion.steps()
    }
}

/// Independent endpoint, action/reaction and whole-body energy checks.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrMixedModeAudit {
    pub displacement_residual_m: [f64; 3],
    pub velocity_residual_m_s: [f64; 3],
    pub marker_normal_displacement_residual_m: f64,
    pub marker_normal_velocity_residual_m_s: f64,
    /// Native kinetic change minus compensated resultant impulse work [J].
    pub impulse_work_residual_j: f64,
    /// Native kinetic change minus actual contact and external work [J].
    pub work_energy_residual_j: f64,
    /// Native signed contact work plus original conservative source work [J].
    pub contact_work_residual_j: f64,
    /// Native kinetic + original contact storage change - external work [J].
    pub total_energy_residual_j: f64,
    pub body_linear_momentum_residual_ns: [Vec3; 2],
    pub body_angular_momentum_residual_nms: [Vec3; 2],
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
}

/// Original coupled path and actual native motion, in original A/B order.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrMixedModeReport {
    pub motion: JkrMixedModeMotionStep,
    /// Full actual J M^-1 J^T for [approach, shear 1, shear 2] [1/kg].
    pub mobility_per_kg: [[f64; 3]; 3],
    pub force_impulse_on_a_world_ns: Vec3,
    pub external_wrenches: [Wrench; 2],
    pub driving_acceleration_m_s2: [f64; 3],
    pub relative_velocity_m_s: [[f64; 3]; 2],
    pub relative_spin_rad_s: [f64; 2],
    pub relative_displacement_m: [f64; 3],
    pub body_displacement_world_m: [Vec3; 2],
    /// Integrated held-body-frame rates, not a finite-rotation instruction.
    pub angular_displacement_body_rad: [Vec3; 2],
    pub marker_displacement_world_m: [Vec3; 2],
    pub reaction_displacement_world_m: [Vec3; 2],
    pub kinetic_energy_change_j: f64,
    pub impulse_work_j: f64,
    pub contact_impulse_work_j: f64,
    pub external_work_j: f64,
    pub contact_storage_change_j: f64,
    pub audit: RigidJkrMixedModeAudit,
}

/// Stage both native momenta and their original elastic history together.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrMixedModeCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub state_after: RigidJkrMixedModeState,
    pub report: RigidJkrMixedModeReport,
}

#[derive(Debug)]
pub enum RigidJkrMixedModeError {
    Invalid(&'static str),
    Motion(JkrMixedModeMotionError),
    Mechanics(RigidHertzError),
    StateMismatch,
    GeometryMismatch { residual_m: f64 },
    UnsupportedTorsion { before_rad_s: f64, after_rad_s: f64 },
    ParameterWorkBalance { residual_j: f64, tolerance_j: f64 },
    Audit(Box<RigidJkrMixedModeAudit>),
}
impl From<JkrMixedModeMotionError> for RigidJkrMixedModeError {
    fn from(error: JkrMixedModeMotionError) -> Self {
        Self::Motion(error)
    }
}
impl From<RigidHertzError> for RigidJkrMixedModeError {
    fn from(error: RigidHertzError) -> Self {
        Self::Mechanics(error)
    }
}
impl From<DynamicsError> for RigidJkrMixedModeError {
    fn from(error: DynamicsError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl From<RigidFrictionError> for RigidJkrMixedModeError {
    fn from(error: RigidFrictionError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl core::fmt::Display for RigidJkrMixedModeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native full-stick mixed-mode JKR: {self:?}")
    }
}
impl std::error::Error for RigidJkrMixedModeError {}

/// Captured physical pair with reusable original source solve storage.
/// Every trial starts at the same native endpoints; scratch reuse changes no
/// accepted source, body or elastic history. No caller effective mass is used.
pub struct PreparedRigidJkrMixedModePair<'a> {
    contact: RigidJkrMixedModeContact<'a>,
    before: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    marker_arms: [Vec3; 2],
    reaction_arms: [Vec3; 2],
    axes: [Vec3; 3],
    markers_before: [PointKinematics; 2],
    reaction_before: [PointKinematics; 2],
    approach_m: f64,
    mobility: [[f64; 3]; 3],
    motion: JkrMixedModeMotion<'a>,
    controls: RigidJkrMixedModeControls,
}

impl<'a> PreparedRigidJkrMixedModePair<'a> {
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contact: RigidJkrMixedModeContact<'a>,
        controls: RigidJkrMixedModeControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidJkrMixedModeError> {
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
                controls.angular_velocity_tolerance_rad_s,
                controls.linear_momentum_tolerance_ns,
                controls.angular_momentum_tolerance_nms,
                controls.energy_tolerance_j,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(RigidJkrMixedModeError::Invalid(
                "bounded distinct bodies and finite nonnegative controls required",
            ));
        }
        let before = [bodies[binding.body_a], bodies[binding.body_b]];
        let properties = [properties[binding.body_a], properties[binding.body_b]];
        let CapturedGeometry {
            marker_arms,
            reaction_arms,
            axes,
            markers_before,
            reaction_before,
            approach_m,
            mobility,
        } = capture_geometry(before, properties, contact, controls)?;
        let motion =
            JkrMixedModeMotion::new(contact.source, contact.shear, mobility, controls.motion)?;
        poll(gate)?;
        Ok(Self {
            contact,
            before,
            properties,
            marker_arms,
            reaction_arms,
            axes,
            markers_before,
            reaction_before,
            approach_m,
            mobility,
            motion,
            controls,
        })
    }

    pub const fn contact(&self) -> RigidJkrMixedModeContact<'a> {
        self.contact
    }
    pub const fn approach_m(&self) -> f64 {
        self.approach_m
    }
    pub const fn mobility_per_kg(&self) -> [[f64; 3]; 3] {
        self.mobility
    }
    pub const fn reaction_arms_body(&self) -> [Vec3; 2] {
        self.reaction_arms
    }

    /// Declare the original elastic shear preload. Normal approach and all
    /// three initial rates come from actual captured native kinematics.
    pub fn initialize(
        &self,
        elastic_shear_m: [f64; 2],
        gate: &CancelGate,
    ) -> Result<RigidJkrMixedModeState, RigidJkrMixedModeError> {
        poll(gate)?;
        let velocity = relative_point(self.reaction_before);
        let spin = finite(relative_spin(self.reaction_before).dot(self.axes[0]))?;
        self.admit_spin([spin; 2])?;
        let motion = self.motion.initialize(
            [self.approach_m, elastic_shear_m[0], elastic_shear_m[1]],
            self.axes.map(|axis| velocity.dot(axis)),
            gate,
        )?;
        Ok(RigidJkrMixedModeState {
            motion,
            frame: self.contact.frame,
            marker_arms: self.marker_arms,
            reference_approach_m: self.contact.reference_approach_m,
            endpoint: NativeEndpoint {
                body_indices: [self.contact.binding.body_a, self.contact.binding.body_b],
                properties: self.properties,
                bodies: self.before,
                pending_com_travel_m: [Vec3::ZERO; 2],
                pending_angular_motion: [false; 2],
            },
        })
    }

    pub fn advance(
        &mut self,
        previous: RigidJkrMixedModeState,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<RigidJkrMixedModeCandidate, RigidJkrMixedModeError> {
        self.advance_with_wrenches(previous, [Wrench::ZERO; 2], dt_s, gate)
    }

    /// One joint conservative path under held actual COM forces/body torques.
    /// Elastic assisting work is admitted; only the original full-stick domain
    /// and real static capacity may authorize the retained branch.
    #[allow(clippy::too_many_lines)]
    pub fn advance_with_wrenches(
        &mut self,
        previous: RigidJkrMixedModeState,
        wrenches: [Wrench; 2],
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<RigidJkrMixedModeCandidate, RigidJkrMixedModeError> {
        let (rates_before, spin_before) = self.admit_history(previous, gate)?;
        let old = previous.position_m();
        let initial = relative_point(self.reaction_before);
        let mut free = [Vec3::ZERO; 2];
        for i in 0..2 {
            wrenches[i].validate()?;
            let alpha =
                self.properties[i].angular_velocity_body_checked(wrenches[i].torque_body)?;
            free[i] = wrenches[i]
                .force_world
                .scale(self.properties[i].mass().recip())
                .add(
                    self.before[i]
                        .pose()
                        .orientation()
                        .rotate_body_to_world(alpha.cross(self.reaction_arms[i])),
                );
        }
        let acceleration = self.axes.map(|axis| free[0].sub(free[1]).dot(axis));
        let step = self
            .motion
            .advance(previous.motion, acceleration, dt_s, gate)?;
        let impulse = self
            .axes
            .into_iter()
            .zip(step.force_impulse_n_s)
            .fold(Vec3::ZERO, |sum, (axis, value)| sum.sub(axis.scale(value)));
        let mut after = self.before;
        let mut kinetic = Dd::from_f64(0.0);
        let mut work = Dd::from_f64(0.0);
        let mut linear = [Dd::from_f64(0.0); 3];
        let mut angular = [Dd::from_f64(0.0); 3];
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
            add_components(&mut resultant.angular, external_torque);
            let receipt = body.apply_momentum_impulse(
                self.properties[i],
                rounded_components(resultant.linear),
                rounded_components(resultant.angular),
            )?;
            after[i] = receipt.state_after;
            linear_residual[i] = receipt.linear_residual_world;
            angular_residual[i] = receipt.angular_residual_body;
            kinetic = kinetic + Dd::from_f64(receipt.work.kinetic_energy_change_j);
            work = work + Dd::from_f64(receipt.work.impulse_work_j);
            add_components(&mut linear, receipt.linear_increment_world);
            add_components(&mut external_linear, external_impulse);
            add_components(
                &mut angular,
                body.pose()
                    .position_world()
                    .cross(receipt.linear_increment_world),
            );
            add_components(
                &mut angular,
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
        let kinetic = finite(kinetic.to_f64())?;
        let work = finite(work.to_f64())?;
        let reaction_after = points(after, self.properties, self.reaction_arms)?;
        let markers_after = points(after, self.properties, self.marker_arms)?;
        let endpoint = relative_point(reaction_after);
        let rates_after = self.axes.map(|axis| endpoint.dot(axis));
        let spin_after = finite(relative_spin(reaction_after).dot(self.axes[0]))?;
        // The held mean effort and held external wrenches give affine angular
        // rates; both endpoints bound unrepresented torsion over this phase.
        self.admit_spin([spin_before, spin_after])?;
        let mut body_travel = [Vec3::ZERO; 2];
        let mut angular_travel = [Vec3::ZERO; 2];
        let mut marker_travel = [Vec3::ZERO; 2];
        let mut reaction_travel = [Vec3::ZERO; 2];
        let mut external_work = Dd::from_f64(0.0);
        for i in 0..2 {
            let orientation = self.before[i].pose().orientation();
            body_travel[i] = mean(
                self.reaction_before[i].center_of_mass_velocity_world,
                reaction_after[i].center_of_mass_velocity_world,
            )
            .scale(dt_s);
            angular_travel[i] = mean(
                self.reaction_before[i].angular_velocity_body,
                reaction_after[i].angular_velocity_body,
            )
            .scale(dt_s);
            marker_travel[i] = body_travel[i].add(
                orientation.rotate_body_to_world(angular_travel[i].cross(self.marker_arms[i])),
            );
            reaction_travel[i] = body_travel[i].add(
                orientation.rotate_body_to_world(angular_travel[i].cross(self.reaction_arms[i])),
            );
            external_work = external_work
                + Dd::from_f64(wrenches[i].force_world.dot(body_travel[i]))
                + Dd::from_f64(wrenches[i].torque_body.dot(angular_travel[i]));
        }
        let external_work = finite(external_work.to_f64())?;
        let relative_travel = reaction_travel[0].sub(reaction_travel[1]);
        let displacement = self.axes.map(|axis| relative_travel.dot(axis));
        let marker_velocity_residual = (relative_point(self.markers_before).dot(self.axes[0])
            - rates_before[0])
            .abs()
            .max((relative_point(markers_after).dot(self.axes[0]) - rates_after[0]).abs());
        let contact_work = finite(impulse.dot(mean(initial, endpoint)))?;
        let linear_balance =
            rounded_components(core::array::from_fn(|i| linear[i] - external_linear[i]));
        let angular_balance =
            rounded_components(core::array::from_fn(|i| angular[i] - external_angular[i]));
        let audit = RigidJkrMixedModeAudit {
            displacement_residual_m: core::array::from_fn(|i| {
                displacement[i] - (step.state().position_m()[i] - old[i])
            }),
            velocity_residual_m_s: core::array::from_fn(|i| {
                rates_after[i] - step.state().velocity_m_s()[i]
            }),
            marker_normal_displacement_residual_m: marker_travel[0]
                .sub(marker_travel[1])
                .dot(self.axes[0])
                - displacement[0],
            marker_normal_velocity_residual_m_s: marker_velocity_residual,
            impulse_work_residual_j: kinetic - work,
            work_energy_residual_j: kinetic - contact_work - external_work,
            contact_work_residual_j: contact_work + step.contact_work_j,
            total_energy_residual_j: kinetic + step.contact_work_j - external_work,
            body_linear_momentum_residual_ns: linear_residual,
            body_angular_momentum_residual_nms: angular_residual,
            linear_momentum_change_ns: rounded_components(linear),
            angular_momentum_change_nms: rounded_components(angular),
            linear_momentum_residual_ns: linear_balance,
            angular_momentum_residual_nms: angular_balance,
        };
        let c = self.controls;
        if audit
            .displacement_residual_m
            .into_iter()
            .chain([audit.marker_normal_displacement_residual_m])
            .any(|v| !v.is_finite() || v.abs() > c.displacement_tolerance_m)
            || audit
                .velocity_residual_m_s
                .into_iter()
                .chain([marker_velocity_residual])
                .any(|v| !v.is_finite() || v.abs() > c.velocity_tolerance_m_s)
            || [
                audit.impulse_work_residual_j,
                audit.work_energy_residual_j,
                audit.contact_work_residual_j,
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
            return Err(RigidJkrMixedModeError::Audit(Box::new(audit)));
        }
        let mut pending_com_travel_m = previous.endpoint.pending_com_travel_m;
        let mut pending_angular_motion = previous.endpoint.pending_angular_motion;
        for i in 0..2 {
            pending_com_travel_m[i] = pending_com_travel_m[i].add(body_travel[i]);
            pending_angular_motion[i] |= angular_travel[i] != Vec3::ZERO;
            if !pending_com_travel_m[i].is_finite() {
                return Err(RigidJkrMixedModeError::Invalid(
                    "unrepresentable pending COM travel",
                ));
            }
        }
        poll(gate)?;
        Ok(RigidJkrMixedModeCandidate {
            bodies_after: after,
            state_after: RigidJkrMixedModeState {
                motion: step.state(),
                endpoint: NativeEndpoint {
                    bodies: after,
                    pending_com_travel_m,
                    pending_angular_motion,
                    ..previous.endpoint
                },
                ..previous
            },
            report: RigidJkrMixedModeReport {
                contact_storage_change_j: step.contact_work_j,
                motion: step,
                mobility_per_kg: self.mobility,
                force_impulse_on_a_world_ns: impulse,
                external_wrenches: wrenches,
                driving_acceleration_m_s2: acceleration,
                relative_velocity_m_s: [rates_before, rates_after],
                relative_spin_rad_s: [spin_before, spin_after],
                relative_displacement_m: displacement,
                body_displacement_world_m: body_travel,
                angular_displacement_body_rad: angular_travel,
                marker_displacement_world_m: marker_travel,
                reaction_displacement_world_m: reaction_travel,
                kinetic_energy_change_j: kinetic,
                impulse_work_j: work,
                contact_impulse_work_j: contact_work,
                external_work_j: external_work,
                audit,
            },
        })
    }

    fn admit_history(
        &self,
        previous: RigidJkrMixedModeState,
        gate: &CancelGate,
    ) -> Result<([f64; 3], f64), RigidJkrMixedModeError> {
        poll(gate)?;
        if previous.frame != self.contact.frame
            || previous.marker_arms != self.marker_arms
            || previous.reference_approach_m != self.contact.reference_approach_m
            || previous.endpoint.properties != self.properties
            || previous.endpoint.body_indices
                != [self.contact.binding.body_a, self.contact.binding.body_b]
        {
            return Err(RigidJkrMixedModeError::StateMismatch);
        }
        let geometry_residual = finite(previous.position_m()[0] - self.approach_m)?;
        if geometry_residual.abs() > self.controls.geometry_tolerance_m {
            return Err(RigidJkrMixedModeError::GeometryMismatch {
                residual_m: geometry_residual,
            });
        }
        let relative = relative_point(self.reaction_before);
        let rates = self.axes.map(|axis| relative.dot(axis));
        if (0..3).any(|i| {
            !rates[i].is_finite()
                || (rates[i] - previous.velocity_m_s()[i]).abs()
                    > self.controls.velocity_tolerance_m_s
        }) {
            return Err(RigidJkrMixedModeError::StateMismatch);
        }
        let spin = finite(relative_spin(self.reaction_before).dot(self.axes[0]))?;
        self.admit_spin([spin; 2])?;
        if previous.endpoint.bodies != self.before {
            return Err(RigidJkrMixedModeError::StateMismatch);
        }
        Ok((rates, spin))
    }

    fn admit_spin(&self, spin: [f64; 2]) -> Result<(), RigidJkrMixedModeError> {
        if spin
            .into_iter()
            .any(|v| !v.is_finite() || v.abs() > self.controls.angular_velocity_tolerance_rad_s)
        {
            Err(RigidJkrMixedModeError::UnsupportedTorsion {
                before_rad_s: spin[0],
                after_rad_s: spin[1],
            })
        } else {
            Ok(())
        }
    }
}

/// Capture the real common reaction plane once, without preparing a duplicate
/// constitutive factor when only the geometry owner changes.
struct CapturedGeometry {
    marker_arms: [Vec3; 2],
    reaction_arms: [Vec3; 2],
    axes: [Vec3; 3],
    markers_before: [PointKinematics; 2],
    reaction_before: [PointKinematics; 2],
    approach_m: f64,
    mobility: [[f64; 3]; 3],
}

fn capture_geometry(
    before: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    contact: RigidJkrMixedModeContact<'_>,
    controls: RigidJkrMixedModeControls,
) -> Result<CapturedGeometry, RigidJkrMixedModeError> {
    let binding = contact.binding;
    let marker_arms = [binding.arm_a_body, binding.arm_b_body];
    let axes = [
        contact.frame.normal(),
        contact.frame.longitudinal(),
        contact.frame.lateral(),
    ]
    .map(vector);
    let markers_before = points(before, properties, marker_arms)?;
    let separation = markers_before[1]
        .point_world
        .sub(markers_before[0].point_world);
    let normal_separation = finite(separation.dot(axes[0]))?;
    if magnitude(separation.sub(axes[0].scale(normal_separation))) > controls.geometry_tolerance_m {
        return Err(RigidJkrMixedModeError::Invalid(
            "normal geometry markers must be collinear",
        ));
    }
    let approach_m = finite(contact.reference_approach_m - normal_separation)?;
    let reaction_arms = core::array::from_fn(|i| {
        marker_arms[i].add(
            before[i]
                .pose()
                .orientation()
                .rotate_world_to_body(separation.scale(if i == 0 { 0.5 } else { -0.5 })),
        )
    });
    let reaction_before = points(before, properties, reaction_arms)?;
    if magnitude(
        reaction_before[0]
            .point_world
            .sub(reaction_before[1].point_world),
    ) > controls.geometry_tolerance_m
    {
        return Err(RigidJkrMixedModeError::Invalid(
            "unrepresentable common reaction plane",
        ));
    }
    let mobility = point_response(before, properties, reaction_arms, axes)?;
    Ok(CapturedGeometry {
        marker_arms,
        reaction_arms,
        axes,
        markers_before,
        reaction_before,
        approach_m,
        mobility,
    })
}

/// Stage the complete native phase; wrenches are in original A/B order.
#[allow(clippy::too_many_arguments)]
pub fn prepare_rigid_jkr_mixed_mode(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidJkrMixedModeContact<'_>,
    previous: RigidJkrMixedModeState,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidJkrMixedModeControls,
    gate: &CancelGate,
) -> Result<RigidJkrMixedModeCandidate, RigidJkrMixedModeError> {
    PreparedRigidJkrMixedModePair::new(bodies, properties, contact, controls, gate)?
        .advance_with_wrenches(previous, wrenches, dt_s, gate)
}

/// Publish both body momenta and original elastic history after every gate.
/// Poses stay held; no fallible operation follows this publication boundary.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_jkr_mixed_mode(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidJkrMixedModeContact<'_>,
    history: &mut RigidJkrMixedModeState,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidJkrMixedModeControls,
    gate: &CancelGate,
) -> Result<RigidJkrMixedModeReport, RigidJkrMixedModeError> {
    let candidate = prepare_rigid_jkr_mixed_mode(
        bodies, properties, contact, *history, wrenches, dt_s, controls, gate,
    )?;
    poll(gate)?;
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    *history = candidate.state_after;
    Ok(candidate.report)
}

fn mean(a: Vec3, b: Vec3) -> Vec3 {
    a.scale(0.5).add(b.scale(0.5))
}
fn finite(value: f64) -> Result<f64, RigidJkrMixedModeError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RigidJkrMixedModeError::Invalid(
            "unrepresentable native mixed-mode JKR mechanics",
        ))
    }
}
