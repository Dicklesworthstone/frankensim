//! Native coupled normal/torsional motion on the original I/III JKR branch.
//!
//! Normal force and torque use the same radius-relaxed source potential. The
//! complete native two-coordinate mobility retains the coupling caused by
//! off-axis arms and full inertia tensors. Twist is an elastic constitutive
//! coordinate in radians, never a contact-radius-scaled displacement.
//!
//! A phase holds the actual Jacobian and applies one mean point force and free
//! couple, together with the declared COM forces and body torques. Poses stay
//! held. Sealed COM and angular travel belongs to the compatible geometry
//! continuation; changing a quaternion cannot reset the elastic twist. Shear,
//! relative rolling, slip, radial instability and detachment are not supplied
//! by this reversible attached-branch owner.

/// Compatible axial pose drift retaining the original source chronology.
pub mod geometry;

use crate::{DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3, Wrench};
use fs_contact::normal_patch::adhesive::mixed_mode::{
    SourceJkrMixedMode,
    torsion::motion::{
        JkrTorsionalMotion, JkrTorsionalMotionBudget, JkrTorsionalMotionError,
        JkrTorsionalMotionState, JkrTorsionalMotionStep,
    },
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_tribo::partial_slip::TangentFrame;

use super::super::hertz::{RigidHertzError, points, relative_point, relative_spin};
use super::super::{
    ResultantImpulse, RigidContactBinding, RigidFrictionError, add_components, magnitude, poll,
    rounded_components, vector,
};

/// Original ordered source and immutable material-marker normal datum.
#[derive(Clone, Copy, Debug)]
pub struct RigidJkrTorsionalContact<'a> {
    pub binding: RigidContactBinding,
    /// Positive A-minus-B normal closure and spin increase [approach, twist].
    pub frame: TangentFrame,
    pub reference_approach_m: f64,
    pub source: &'a SourceJkrMixedMode,
}

/// Original source solve limits and dimensional native acceptance gates.
#[derive(Clone, Copy, Debug)]
pub struct RigidJkrTorsionalControls {
    pub motion: JkrTorsionalMotionBudget,
    pub maximum_bodies: usize,
    pub geometry_tolerance_m: f64,
    pub displacement_tolerance_m: f64,
    pub twist_tolerance_rad: f64,
    /// Also bounds the unrepresented tangential point velocity.
    pub velocity_tolerance_m_s: f64,
    /// Also bounds the unrepresented relative rolling angular velocity.
    pub angular_velocity_tolerance_rad_s: f64,
    pub linear_momentum_tolerance_ns: f64,
    pub angular_momentum_tolerance_nms: f64,
    pub energy_tolerance_j: f64,
}

/// Original source chronology and its actual captured native endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidJkrTorsionalState {
    motion: JkrTorsionalMotionState,
    frame: TangentFrame,
    marker_arms: [Vec3; 2],
    reference_approach_m: f64,
    endpoint: NativeEndpoint,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct NativeEndpoint {
    body_indices: [usize; 2],
    properties: [MassProperties; 2],
    bodies: [RigidBodyState; 2],
    pending_com_travel_m: [Vec3; 2],
    pending_angular_displacement_body_rad: [Vec3; 2],
    // Sticky across pending phases: transverse motion cannot be erased by an
    // out-and-back angular excursion with zero net integrated angle.
    axial_drift_admissible: bool,
}

impl RigidJkrTorsionalState {
    pub const fn motion_state(self) -> JkrTorsionalMotionState {
        self.motion
    }
    /// [Normal approach in m, retained elastic twist in rad].
    pub const fn coordinates(self) -> [f64; 2] {
        self.motion.coordinates()
    }
    /// [Normal closure in m/s, relative normal spin in rad/s].
    pub const fn rates(self) -> [f64; 2] {
        self.motion.rates()
    }
    pub const fn time_s(self) -> f64 {
        self.motion.time_s()
    }
    pub const fn steps(self) -> u64 {
        self.motion.steps()
    }
}

/// Independent actual-body work, momentum and endpoint checks.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrTorsionalAudit {
    /// [m, rad], in the original source coordinate order.
    pub coordinate_residual: [f64; 2],
    /// [m/s, rad/s], in the original source coordinate order.
    pub rate_residual: [f64; 2],
    pub marker_normal_displacement_residual_m: f64,
    pub marker_normal_velocity_residual_m_s: f64,
    pub impulse_work_residual_j: f64,
    pub work_energy_residual_j: f64,
    pub contact_work_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub body_linear_momentum_residual_ns: [Vec3; 2],
    pub body_angular_momentum_residual_nms: [Vec3; 2],
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
}

/// Original source path and actual native motion, in binding A/B order.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrTorsionalReport {
    pub motion: JkrTorsionalMotionStep,
    /// Actual J M^-1 J^T. Units are [[1/kg,1/(kg m)],
    /// [1/(kg m),1/(kg m^2)]]; radians are dimensionless.
    pub mobility: [[f64; 2]; 2],
    pub force_impulse_on_a_world_ns: Vec3,
    pub angular_impulse_on_a_world_nms: Vec3,
    pub external_wrenches: [Wrench; 2],
    /// [m/s^2, rad/s^2], from actual held native wrenches.
    pub driving_acceleration: [f64; 2],
    /// Initial/final [normal closure m/s, normal spin rad/s].
    pub rates: [[f64; 2]; 2],
    /// Actual integrated [normal approach m, twist rad].
    pub coordinate_increment: [f64; 2],
    pub body_displacement_world_m: [Vec3; 2],
    /// Integrated held-body-frame rates, not a general finite-rotation path.
    pub angular_displacement_body_rad: [Vec3; 2],
    pub marker_displacement_world_m: [Vec3; 2],
    pub reaction_displacement_world_m: [Vec3; 2],
    pub kinetic_energy_change_j: f64,
    pub impulse_work_j: f64,
    pub contact_impulse_work_j: f64,
    pub external_work_j: f64,
    pub contact_storage_change_j: f64,
    pub audit: RigidJkrTorsionalAudit,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrTorsionalCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub state_after: RigidJkrTorsionalState,
    pub report: RigidJkrTorsionalReport,
}

#[derive(Debug)]
pub enum RigidJkrTorsionalError {
    Invalid(&'static str),
    Motion(JkrTorsionalMotionError),
    Mechanics(RigidHertzError),
    StateMismatch,
    GeometryMismatch { residual_m: f64 },
    UnsupportedShear { velocity_m_s: Vec3 },
    UnsupportedRolling { angular_velocity_rad_s: Vec3 },
    Audit(Box<RigidJkrTorsionalAudit>),
}
impl From<JkrTorsionalMotionError> for RigidJkrTorsionalError {
    fn from(error: JkrTorsionalMotionError) -> Self {
        Self::Motion(error)
    }
}
impl From<RigidHertzError> for RigidJkrTorsionalError {
    fn from(error: RigidHertzError) -> Self {
        Self::Mechanics(error)
    }
}
impl From<DynamicsError> for RigidJkrTorsionalError {
    fn from(error: DynamicsError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl From<RigidFrictionError> for RigidJkrTorsionalError {
    fn from(error: RigidFrictionError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl core::fmt::Display for RigidJkrTorsionalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native coupled normal/torsional JKR: {self:?}")
    }
}
impl std::error::Error for RigidJkrTorsionalError {}

/// Captured native geometry with reusable original constitutive storage.
pub struct PreparedRigidJkrTorsionalPair<'a> {
    contact: RigidJkrTorsionalContact<'a>,
    before: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    marker_arms: [Vec3; 2],
    reaction_arms: [Vec3; 2],
    axes: [Vec3; 3],
    markers_before: [PointKinematics; 2],
    reaction_before: [PointKinematics; 2],
    approach_m: f64,
    mobility: [[f64; 2]; 2],
    motion: JkrTorsionalMotion<'a>,
    controls: RigidJkrTorsionalControls,
}

impl<'a> PreparedRigidJkrTorsionalPair<'a> {
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contact: RigidJkrTorsionalContact<'a>,
        controls: RigidJkrTorsionalControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidJkrTorsionalError> {
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
                controls.twist_tolerance_rad,
                controls.velocity_tolerance_m_s,
                controls.angular_velocity_tolerance_rad_s,
                controls.linear_momentum_tolerance_ns,
                controls.angular_momentum_tolerance_nms,
                controls.energy_tolerance_j,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(RigidJkrTorsionalError::Invalid(
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
        let motion = JkrTorsionalMotion::new(contact.source, mobility, controls.motion)?;
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

    pub const fn contact(&self) -> RigidJkrTorsionalContact<'a> {
        self.contact
    }
    pub const fn approach_m(&self) -> f64 {
        self.approach_m
    }
    pub const fn mobility(&self) -> [[f64; 2]; 2] {
        self.mobility
    }
    pub const fn reaction_arms_body(&self) -> [Vec3; 2] {
        self.reaction_arms
    }

    /// Declare only the original elastic twist preload. Approach and both
    /// rates are authenticated by the actual native pair.
    pub fn initialize(
        &self,
        elastic_twist_rad: f64,
        gate: &CancelGate,
    ) -> Result<RigidJkrTorsionalState, RigidJkrTorsionalError> {
        poll(gate)?;
        self.admit_unrepresented(self.reaction_before)?;
        let motion = self.motion.initial_state(
            [self.approach_m, elastic_twist_rad],
            projected_rates(self.reaction_before, self.axes[0])?,
            gate,
        )?;
        Ok(RigidJkrTorsionalState {
            motion,
            frame: self.contact.frame,
            marker_arms: self.marker_arms,
            reference_approach_m: self.contact.reference_approach_m,
            endpoint: NativeEndpoint {
                body_indices: [self.contact.binding.body_a, self.contact.binding.body_b],
                properties: self.properties,
                bodies: self.before,
                pending_com_travel_m: [Vec3::ZERO; 2],
                pending_angular_displacement_body_rad: [Vec3::ZERO; 2],
                axial_drift_admissible: axial_drift_admissible(
                    self.before,
                    self.reaction_before,
                    self.axes[0],
                    self.controls,
                ),
            },
        })
    }

    pub fn advance(
        &self,
        previous: RigidJkrTorsionalState,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<RigidJkrTorsionalCandidate, RigidJkrTorsionalError> {
        self.advance_with_wrenches(previous, [Wrench::ZERO; 2], dt_s, gate)
    }

    /// Joint reversible recoil under held native COM forces and body torques.
    #[allow(clippy::too_many_lines)]
    pub fn advance_with_wrenches(
        &self,
        previous: RigidJkrTorsionalState,
        wrenches: [Wrench; 2],
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<RigidJkrTorsionalCandidate, RigidJkrTorsionalError> {
        let rates_before = self.admit_history(previous, gate)?;
        let old = previous.coordinates();
        let mut free_point = [Vec3::ZERO; 2];
        let mut free_spin = [Vec3::ZERO; 2];
        for i in 0..2 {
            wrenches[i].validate()?;
            let alpha =
                self.properties[i].angular_velocity_body_checked(wrenches[i].torque_body)?;
            let orientation = self.before[i].pose().orientation();
            free_spin[i] = orientation.rotate_body_to_world(alpha);
            free_point[i] = wrenches[i]
                .force_world
                .scale(self.properties[i].mass().recip())
                .add(orientation.rotate_body_to_world(alpha.cross(self.reaction_arms[i])));
        }
        let acceleration = [
            finite(free_point[0].sub(free_point[1]).dot(self.axes[0]))?,
            finite(free_spin[0].sub(free_spin[1]).dot(self.axes[0]))?,
        ];
        let step = self
            .motion
            .advance(previous.motion, dt_s, acceleration, gate)?;
        let impulse = self.axes[0].scale(-step.effort_impulse[0]);
        let couple = self.axes[0].scale(-step.effort_impulse[1]);
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
            resultant.add(
                body,
                self.reaction_arms[i],
                impulse.scale(sign),
                couple.scale(sign),
            );
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
        let rates_after = projected_rates(reaction_after, self.axes[0])?;
        // Held forces/couples give affine rates. Endpoint admission therefore
        // bounds every unrepresented shear and rolling rate of this phase.
        self.admit_unrepresented(reaction_after)?;
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
        // Cancel common body translation before forming small relative work.
        let displacement =
            core::array::from_fn(|i| dt_s * (0.5 * rates_before[i] + 0.5 * rates_after[i]));
        let marker_rates = [
            relative_point(self.markers_before)?.dot(self.axes[0]),
            relative_point(markers_after)?.dot(self.axes[0]),
        ];
        let marker_velocity_residual = (marker_rates[0] - rates_before[0])
            .abs()
            .max((marker_rates[1] - rates_after[0]).abs());
        let contact_work = finite(
            impulse.dot(mean(
                relative_point(self.reaction_before)?,
                relative_point(reaction_after)?,
            )) + couple.dot(mean(
                relative_spin(self.reaction_before),
                relative_spin(reaction_after),
            )),
        )?;
        let linear_balance =
            rounded_components(core::array::from_fn(|i| linear[i] - external_linear[i]));
        let angular_balance =
            rounded_components(core::array::from_fn(|i| angular[i] - external_angular[i]));
        let audit = RigidJkrTorsionalAudit {
            coordinate_residual: core::array::from_fn(|i| {
                displacement[i] - (step.state().coordinates()[i] - old[i])
            }),
            rate_residual: core::array::from_fn(|i| rates_after[i] - step.state().rates()[i]),
            marker_normal_displacement_residual_m: dt_s
                * (0.5 * marker_rates[0] + 0.5 * marker_rates[1])
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
        let coordinate_limits = [c.displacement_tolerance_m, c.twist_tolerance_rad];
        let rate_limits = [c.velocity_tolerance_m_s, c.angular_velocity_tolerance_rad_s];
        if (0..2).any(|i| {
            !audit.coordinate_residual[i].is_finite()
                || audit.coordinate_residual[i].abs() > coordinate_limits[i]
                || !audit.rate_residual[i].is_finite()
                || audit.rate_residual[i].abs() > rate_limits[i]
        }) || !audit.marker_normal_displacement_residual_m.is_finite()
            || audit.marker_normal_displacement_residual_m.abs() > c.displacement_tolerance_m
            || !marker_velocity_residual.is_finite()
            || marker_velocity_residual > c.velocity_tolerance_m_s
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
            return Err(RigidJkrTorsionalError::Audit(Box::new(audit)));
        }
        let mut pending_com_travel_m = previous.endpoint.pending_com_travel_m;
        let mut pending_angular_displacement_body_rad =
            previous.endpoint.pending_angular_displacement_body_rad;
        for i in 0..2 {
            pending_com_travel_m[i] = pending_com_travel_m[i].add(body_travel[i]);
            pending_angular_displacement_body_rad[i] =
                pending_angular_displacement_body_rad[i].add(angular_travel[i]);
            if !pending_com_travel_m[i].is_finite()
                || !pending_angular_displacement_body_rad[i].is_finite()
            {
                return Err(RigidJkrTorsionalError::Invalid(
                    "unrepresentable pending native travel",
                ));
            }
        }
        let axial_drift_admissible = previous.endpoint.axial_drift_admissible
            && axial_drift_admissible(self.before, self.reaction_before, self.axes[0], c)
            && axial_drift_admissible(after, reaction_after, self.axes[0], c);
        poll(gate)?;
        Ok(RigidJkrTorsionalCandidate {
            bodies_after: after,
            state_after: RigidJkrTorsionalState {
                motion: step.state(),
                endpoint: NativeEndpoint {
                    bodies: after,
                    pending_com_travel_m,
                    pending_angular_displacement_body_rad,
                    axial_drift_admissible,
                    ..previous.endpoint
                },
                ..previous
            },
            report: RigidJkrTorsionalReport {
                contact_storage_change_j: step.contact_work_j,
                motion: step,
                mobility: self.mobility,
                force_impulse_on_a_world_ns: impulse,
                angular_impulse_on_a_world_nms: couple,
                external_wrenches: wrenches,
                driving_acceleration: acceleration,
                rates: [rates_before, rates_after],
                coordinate_increment: displacement,
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
        previous: RigidJkrTorsionalState,
        gate: &CancelGate,
    ) -> Result<[f64; 2], RigidJkrTorsionalError> {
        poll(gate)?;
        if previous.frame != self.contact.frame
            || previous.marker_arms != self.marker_arms
            || previous.reference_approach_m != self.contact.reference_approach_m
            || previous.endpoint.properties != self.properties
            || previous.endpoint.body_indices
                != [self.contact.binding.body_a, self.contact.binding.body_b]
            || previous.endpoint.bodies != self.before
        {
            return Err(RigidJkrTorsionalError::StateMismatch);
        }
        let residual_m = finite(previous.coordinates()[0] - self.approach_m)?;
        if residual_m.abs() > self.controls.geometry_tolerance_m {
            return Err(RigidJkrTorsionalError::GeometryMismatch { residual_m });
        }
        let rates = projected_rates(self.reaction_before, self.axes[0])?;
        let limits = [
            self.controls.velocity_tolerance_m_s,
            self.controls.angular_velocity_tolerance_rad_s,
        ];
        if (0..2).any(|i| (rates[i] - previous.rates()[i]).abs() > limits[i]) {
            return Err(RigidJkrTorsionalError::StateMismatch);
        }
        self.admit_unrepresented(self.reaction_before)?;
        Ok(rates)
    }

    fn admit_unrepresented(
        &self,
        points: [PointKinematics; 2],
    ) -> Result<(), RigidJkrTorsionalError> {
        let velocity = relative_point(points)?;
        let shear = velocity.sub(self.axes[0].scale(velocity.dot(self.axes[0])));
        if !shear.is_finite() || magnitude(shear) > self.controls.velocity_tolerance_m_s {
            return Err(RigidJkrTorsionalError::UnsupportedShear {
                velocity_m_s: shear,
            });
        }
        let spin = relative_spin(points);
        let rolling = spin.sub(self.axes[0].scale(spin.dot(self.axes[0])));
        if !rolling.is_finite()
            || magnitude(rolling) > self.controls.angular_velocity_tolerance_rad_s
        {
            return Err(RigidJkrTorsionalError::UnsupportedRolling {
                angular_velocity_rad_s: rolling,
            });
        }
        Ok(())
    }
}

struct CapturedGeometry {
    marker_arms: [Vec3; 2],
    reaction_arms: [Vec3; 2],
    axes: [Vec3; 3],
    markers_before: [PointKinematics; 2],
    reaction_before: [PointKinematics; 2],
    approach_m: f64,
    mobility: [[f64; 2]; 2],
}

fn capture_geometry(
    before: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    contact: RigidJkrTorsionalContact<'_>,
    controls: RigidJkrTorsionalControls,
) -> Result<CapturedGeometry, RigidJkrTorsionalError> {
    let marker_arms = [contact.binding.arm_a_body, contact.binding.arm_b_body];
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
        return Err(RigidJkrTorsionalError::Invalid(
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
        return Err(RigidJkrTorsionalError::Invalid(
            "unrepresentable common reaction plane",
        ));
    }
    let mobility = normal_torsion_response(before, properties, reaction_arms, axes[0])?;
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

/// The first row is a point force; the second is a unit free couple. A couple
/// has no translational row and no changing contact-radius multiplier.
fn normal_torsion_response(
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    arms: [Vec3; 2],
    normal: Vec3,
) -> Result<[[f64; 2]; 2], RigidJkrTorsionalError> {
    let mut factor = [[0.0; 12]; 2];
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let direction = bodies[i].pose().orientation().rotate_world_to_body(normal);
        let translation = normal.scale(sign / properties[i].mass().sqrt());
        let torques = [arms[i].cross(direction).scale(sign), direction.scale(sign)];
        factor[0][6 * i..6 * i + 3].copy_from_slice(&[translation.x, translation.y, translation.z]);
        for j in 0..2 {
            let rotation = properties[i].angular_mobility_factor(torques[j]);
            factor[j][6 * i + 3..6 * i + 6].copy_from_slice(&[rotation.x, rotation.y, rotation.z]);
        }
    }
    let mut result = [[0.0; 2]; 2];
    for i in 0..2 {
        for j in 0..=i {
            result[i][j] = finite(factor[i].iter().zip(factor[j]).map(|(a, b)| a * b).sum())?;
            result[j][i] = result[i][j];
        }
    }
    Ok(result)
}

fn projected_rates(
    points: [PointKinematics; 2],
    normal: Vec3,
) -> Result<[f64; 2], RigidJkrTorsionalError> {
    Ok([
        finite(relative_point(points)?.dot(normal))?,
        finite(relative_spin(points).dot(normal))?,
    ])
}

/// Only fixed-axis principal spins support automatic finite-pose continuation.
/// This restriction does not discard any inertia entry of the held phase.
fn axial_drift_admissible(
    bodies: [RigidBodyState; 2],
    points: [PointKinematics; 2],
    normal: Vec3,
    _controls: RigidJkrTorsionalControls,
) -> bool {
    (0..2).all(|i| {
        let omega = points[i].angular_velocity_world;
        let momentum = bodies[i]
            .pose()
            .orientation()
            .rotate_body_to_world(bodies[i].angular_momentum_body());
        let transverse_omega = omega.sub(normal.scale(omega.dot(normal)));
        let transverse_momentum = momentum.sub(normal.scale(momentum.dot(normal)));
        // Native audit tolerances cannot authorize a different finite-rotation
        // path. Only representation roundoff may depart from the fixed axis.
        transverse_omega.is_finite()
            && transverse_momentum.is_finite()
            && magnitude(transverse_omega) <= 128.0 * f64::EPSILON * magnitude(omega)
            && magnitude(transverse_momentum) <= 128.0 * f64::EPSILON * magnitude(momentum)
    })
}

/// Stage the native body pair and its unchanged original source association.
#[allow(clippy::too_many_arguments)]
pub fn prepare_rigid_jkr_torsional(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidJkrTorsionalContact<'_>,
    previous: RigidJkrTorsionalState,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidJkrTorsionalControls,
    gate: &CancelGate,
) -> Result<RigidJkrTorsionalCandidate, RigidJkrTorsionalError> {
    PreparedRigidJkrTorsionalPair::new(bodies, properties, contact, controls, gate)?
        .advance_with_wrenches(previous, wrenches, dt_s, gate)
}

/// Atomically publish both momenta and the elastic state after all gates.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_jkr_torsional(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidJkrTorsionalContact<'_>,
    history: &mut RigidJkrTorsionalState,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidJkrTorsionalControls,
    gate: &CancelGate,
) -> Result<RigidJkrTorsionalReport, RigidJkrTorsionalError> {
    let candidate = prepare_rigid_jkr_torsional(
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
fn finite(value: f64) -> Result<f64, RigidJkrTorsionalError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RigidJkrTorsionalError::Invalid(
            "unrepresentable native normal/torsional JKR mechanics",
        ))
    }
}
