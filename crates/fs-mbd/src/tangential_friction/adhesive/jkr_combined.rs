//! Native normal, two-direction full-stick shear and twist on one JKR disk.
//!
//! All four efforts come from the same original radius-relaxed brittle source.
//! Actual masses, complete inertia tensors and one coincident reaction point
//! provide the full mixed-unit mobility. Elastic twist is an angle; changing
//! contact radius never rescales its coordinate or its physical free couple.
//!
//! This is a held-geometry momentum phase with held external COM forces and
//! body torques. It retains constitutive history and reports all integrated
//! native travel. [`geometry`] consumes compatible travel and commuting axial
//! spin without restarting elastic history; unrepresented finite attitudes and
//! transverse marker migration refuse. Relative rolling, slip, debonding and
//! detached trajectories are outside this attached full-stick source.
//! [`temperature`] continues the same loaded source at held geometry with
//! explicit normal, shear, torsional and surface parameter work.

/// Compatible actual poses and objective vector shear, with the source retained.
pub mod geometry;
/// Original material-temperature continuation of the loaded four-coordinate pair.
pub mod temperature;

use crate::{DynamicsError, MassProperties, PointKinematics, RigidBodyState, Vec3, Wrench};
use fs_contact::normal_patch::adhesive::mixed_mode::{
    SourceJkrMixedMode,
    combined::motion::{
        JkrCombinedMotion, JkrCombinedMotionBudget, JkrCombinedMotionError, JkrCombinedMotionState,
        JkrCombinedMotionStep,
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

/// Original ordered material source and immutable normal material-marker datum.
#[derive(Clone, Copy, Debug)]
pub struct RigidJkrCombinedContact<'a> {
    pub binding: RigidContactBinding,
    /// Positive A-minus-B point rates and normal spin increase the coordinates.
    pub frame: TangentFrame,
    pub reference_approach_m: f64,
    pub source: &'a SourceJkrMixedMode,
}

/// Original source solve limits and dimensional actual-body acceptance gates.
#[derive(Clone, Copy, Debug)]
pub struct RigidJkrCombinedControls {
    pub motion: JkrCombinedMotionBudget,
    pub maximum_bodies: usize,
    pub geometry_tolerance_m: f64,
    pub displacement_tolerance_m: f64,
    pub twist_tolerance_rad: f64,
    pub velocity_tolerance_m_s: f64,
    /// Also bounds unrepresented relative transverse rolling throughout a phase.
    pub angular_velocity_tolerance_rad_s: f64,
    pub linear_momentum_tolerance_ns: f64,
    pub angular_momentum_tolerance_nms: f64,
    pub energy_tolerance_j: f64,
}

/// Sealed original constitutive chronology and captured actual-body endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidJkrCombinedState {
    motion: JkrCombinedMotionState,
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
    /// Every held phase had individual principal-axis spin about the fixed normal.
    /// A later cancelling excursion cannot replace the unresolved attitude path.
    axial_drift_admissible: bool,
}

impl RigidJkrCombinedState {
    pub const fn motion_state(self) -> JkrCombinedMotionState {
        self.motion
    }
    /// [Normal approach m, elastic shear m, elastic shear m, elastic twist rad].
    pub const fn coordinates(self) -> [f64; 4] {
        self.motion.coordinates()
    }
    /// [Normal closure m/s, two tangent rates m/s, relative normal spin rad/s].
    pub const fn rates(self) -> [f64; 4] {
        self.motion.rates()
    }
    pub const fn time_s(self) -> f64 {
        self.motion.time_s()
    }
    pub const fn steps(self) -> u64 {
        self.motion.steps()
    }
    /// Accumulated held-geometry COM travel; no poses have been moved by this owner.
    pub const fn pending_com_travel_m(self) -> [Vec3; 2] {
        self.endpoint.pending_com_travel_m
    }
    /// Integrated held-body-frame rates, not a general finite-rotation path.
    pub const fn pending_angular_displacement_body_rad(self) -> [Vec3; 2] {
        self.endpoint.pending_angular_displacement_body_rad
    }
}

/// Independent actual-body work, momentum and source endpoint checks.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrCombinedAudit {
    /// [m,m,m,rad].
    pub coordinate_residual: [f64; 4],
    /// [m/s,m/s,m/s,rad/s].
    pub rate_residual: [f64; 4],
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

/// Original source path and actual motion, in binding A/B order.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrCombinedReport {
    pub motion: JkrCombinedMotionStep,
    /// J M^-1 J^T for [normal,tangent1,tangent2,free normal couple].
    /// The translational block has units 1/kg, mixed entries 1/(kg m),
    /// and the twist diagonal 1/(kg m^2).
    pub mobility: [[f64; 4]; 4],
    pub force_impulse_on_a_world_ns: Vec3,
    /// Free normal couple only; point-force arm torque is applied separately.
    pub angular_impulse_on_a_world_nms: Vec3,
    pub external_wrenches: [Wrench; 2],
    pub driving_acceleration: [f64; 4],
    pub rates: [[f64; 4]; 2],
    pub coordinate_increment: [f64; 4],
    pub body_displacement_world_m: [Vec3; 2],
    pub angular_displacement_body_rad: [Vec3; 2],
    pub marker_displacement_world_m: [Vec3; 2],
    pub reaction_displacement_world_m: [Vec3; 2],
    pub kinetic_energy_change_j: f64,
    pub impulse_work_j: f64,
    pub contact_impulse_work_j: f64,
    pub external_work_j: f64,
    pub contact_storage_change_j: f64,
    pub audit: RigidJkrCombinedAudit,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrCombinedCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub state_after: RigidJkrCombinedState,
    pub report: RigidJkrCombinedReport,
}

#[derive(Debug)]
pub enum RigidJkrCombinedError {
    Invalid(&'static str),
    Motion(JkrCombinedMotionError),
    Mechanics(RigidHertzError),
    StateMismatch,
    GeometryMismatch { residual_m: f64 },
    UnsupportedRolling { angular_velocity_rad_s: Vec3 },
    Audit(Box<RigidJkrCombinedAudit>),
}
impl From<JkrCombinedMotionError> for RigidJkrCombinedError {
    fn from(error: JkrCombinedMotionError) -> Self {
        Self::Motion(error)
    }
}
impl From<RigidHertzError> for RigidJkrCombinedError {
    fn from(error: RigidHertzError) -> Self {
        Self::Mechanics(error)
    }
}
impl From<DynamicsError> for RigidJkrCombinedError {
    fn from(error: DynamicsError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl From<RigidFrictionError> for RigidJkrCombinedError {
    fn from(error: RigidFrictionError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl core::fmt::Display for RigidJkrCombinedError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native combined normal/shear/twist JKR: {self:?}")
    }
}
impl std::error::Error for RigidJkrCombinedError {}

/// Captured actual geometry and the reusable original four-coordinate solver.
pub struct PreparedRigidJkrCombinedPair<'a> {
    contact: RigidJkrCombinedContact<'a>,
    before: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    marker_arms: [Vec3; 2],
    reaction_arms: [Vec3; 2],
    axes: [Vec3; 3],
    markers_before: [PointKinematics; 2],
    reaction_before: [PointKinematics; 2],
    approach_m: f64,
    mobility: [[f64; 4]; 4],
    motion: JkrCombinedMotion<'a>,
    controls: RigidJkrCombinedControls,
}

impl<'a> PreparedRigidJkrCombinedPair<'a> {
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contact: RigidJkrCombinedContact<'a>,
        controls: RigidJkrCombinedControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidJkrCombinedError> {
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
            return Err(RigidJkrCombinedError::Invalid(
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
        let motion = JkrCombinedMotion::new(contact.source, mobility, controls.motion)?;
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

    pub const fn contact(&self) -> RigidJkrCombinedContact<'a> {
        self.contact
    }
    pub const fn bodies(&self) -> &[RigidBodyState; 2] {
        &self.before
    }
    pub const fn properties(&self) -> &[MassProperties; 2] {
        &self.properties
    }
    pub const fn approach_m(&self) -> f64 {
        self.approach_m
    }
    pub const fn mobility(&self) -> [[f64; 4]; 4] {
        self.mobility
    }
    pub const fn reaction_arms_body(&self) -> [Vec3; 2] {
        self.reaction_arms
    }

    /// Declare [elastic shear1 m, shear2 m, twist rad] only. Actual material
    /// geometry supplies approach, and actual momenta supply all four rates.
    pub fn initialize(
        &self,
        elastic_preload: [f64; 3],
        gate: &CancelGate,
    ) -> Result<RigidJkrCombinedState, RigidJkrCombinedError> {
        poll(gate)?;
        self.admit_unrepresented(self.reaction_before)?;
        let motion = self.motion.initial_state(
            [
                self.approach_m,
                elastic_preload[0],
                elastic_preload[1],
                elastic_preload[2],
            ],
            projected_rates(self.reaction_before, self.axes)?,
            gate,
        )?;
        Ok(RigidJkrCombinedState {
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
                axial_drift_admissible: geometry::axial_drift_admissible(
                    self.before,
                    self.reaction_before,
                    self.axes[0],
                ),
            },
        })
    }

    pub fn advance(
        &self,
        previous: RigidJkrCombinedState,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<RigidJkrCombinedCandidate, RigidJkrCombinedError> {
        self.advance_with_wrenches(previous, [Wrench::ZERO; 2], dt_s, gate)
    }

    /// One coupled full-stick momentum phase under actual held native wrenches.
    #[allow(clippy::too_many_lines)]
    pub fn advance_with_wrenches(
        &self,
        previous: RigidJkrCombinedState,
        wrenches: [Wrench; 2],
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<RigidJkrCombinedCandidate, RigidJkrCombinedError> {
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
            finite(free_point[0].sub(free_point[1]).dot(self.axes[1]))?,
            finite(free_point[0].sub(free_point[1]).dot(self.axes[2]))?,
            finite(free_spin[0].sub(free_spin[1]).dot(self.axes[0]))?,
        ];
        let step = self
            .motion
            .advance(previous.motion, dt_s, acceleration, gate)?;
        let impulse = self.axes[0]
            .scale(-step.effort_impulse[0])
            .add(self.axes[1].scale(-step.effort_impulse[1]))
            .add(self.axes[2].scale(-step.effort_impulse[2]));
        let couple = self.axes[0].scale(-step.effort_impulse[3]);
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
        let rates_after = projected_rates(reaction_after, self.axes)?;
        // Rates are affine under this held resultant. Endpoint checks bound
        // unrepresented relative rolling throughout the whole physical phase.
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
        // Cancel common translation before forming small constitutive travel.
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
        let audit = RigidJkrCombinedAudit {
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
        let coordinate_limits = [
            c.displacement_tolerance_m,
            c.displacement_tolerance_m,
            c.displacement_tolerance_m,
            c.twist_tolerance_rad,
        ];
        let rate_limits = [
            c.velocity_tolerance_m_s,
            c.velocity_tolerance_m_s,
            c.velocity_tolerance_m_s,
            c.angular_velocity_tolerance_rad_s,
        ];
        if (0..4).any(|i| {
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
            return Err(RigidJkrCombinedError::Audit(Box::new(audit)));
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
                return Err(RigidJkrCombinedError::Invalid(
                    "unrepresentable pending native travel",
                ));
            }
        }
        let axial_drift_admissible = previous.endpoint.axial_drift_admissible
            && geometry::axial_drift_admissible(self.before, self.reaction_before, self.axes[0])
            && geometry::axial_drift_admissible(after, reaction_after, self.axes[0]);
        poll(gate)?;
        Ok(RigidJkrCombinedCandidate {
            bodies_after: after,
            state_after: RigidJkrCombinedState {
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
            report: RigidJkrCombinedReport {
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
        previous: RigidJkrCombinedState,
        gate: &CancelGate,
    ) -> Result<[f64; 4], RigidJkrCombinedError> {
        poll(gate)?;
        if previous.frame != self.contact.frame
            || previous.marker_arms != self.marker_arms
            || previous.reference_approach_m != self.contact.reference_approach_m
            || previous.endpoint.properties != self.properties
            || previous.endpoint.body_indices
                != [self.contact.binding.body_a, self.contact.binding.body_b]
            || previous.endpoint.bodies != self.before
        {
            return Err(RigidJkrCombinedError::StateMismatch);
        }
        let residual_m = finite(previous.coordinates()[0] - self.approach_m)?;
        if residual_m.abs() > self.controls.geometry_tolerance_m {
            return Err(RigidJkrCombinedError::GeometryMismatch { residual_m });
        }
        let rates = projected_rates(self.reaction_before, self.axes)?;
        let limits = [
            self.controls.velocity_tolerance_m_s,
            self.controls.velocity_tolerance_m_s,
            self.controls.velocity_tolerance_m_s,
            self.controls.angular_velocity_tolerance_rad_s,
        ];
        if (0..4).any(|i| (rates[i] - previous.rates()[i]).abs() > limits[i]) {
            return Err(RigidJkrCombinedError::StateMismatch);
        }
        self.admit_unrepresented(self.reaction_before)?;
        Ok(rates)
    }

    fn admit_unrepresented(
        &self,
        points: [PointKinematics; 2],
    ) -> Result<(), RigidJkrCombinedError> {
        let spin = relative_spin(points);
        let rolling = spin.sub(self.axes[0].scale(spin.dot(self.axes[0])));
        if !rolling.is_finite()
            || magnitude(rolling) > self.controls.angular_velocity_tolerance_rad_s
        {
            return Err(RigidJkrCombinedError::UnsupportedRolling {
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
    mobility: [[f64; 4]; 4],
}

/// The original ordered material markers and coincident reaction construction.
/// Pair-local arrays retain the caller's original binding indices separately.
fn capture_geometry(
    before: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    contact: RigidJkrCombinedContact<'_>,
    controls: RigidJkrCombinedControls,
) -> Result<CapturedGeometry, RigidJkrCombinedError> {
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
    let transverse_separation = separation.sub(axes[0].scale(normal_separation));
    if !transverse_separation.is_finite()
        || magnitude(transverse_separation) > controls.geometry_tolerance_m
    {
        return Err(RigidJkrCombinedError::Invalid(
            "normal material markers must be collinear",
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
    let reaction_gap = reaction_before[0]
        .point_world
        .sub(reaction_before[1].point_world);
    if !reaction_gap.is_finite() || magnitude(reaction_gap) > controls.geometry_tolerance_m {
        return Err(RigidJkrCombinedError::Invalid(
            "unrepresentable coincident reaction point",
        ));
    }
    let mobility = combined_response(before, properties, reaction_arms, axes)?;
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

/// Gram factorization retains every inertia and force/couple cross entry.
fn combined_response(
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    arms: [Vec3; 2],
    axes: [Vec3; 3],
) -> Result<[[f64; 4]; 4], RigidJkrCombinedError> {
    let mut factor = [[0.0; 12]; 4];
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let orientation = bodies[i].pose().orientation();
        for j in 0..4 {
            let direction = orientation.rotate_world_to_body(axes[if j == 3 { 0 } else { j }]);
            let torque = if j == 3 {
                direction.scale(sign)
            } else {
                let translation = axes[j].scale(sign / properties[i].mass().sqrt());
                factor[j][6 * i..6 * i + 3].copy_from_slice(&[
                    translation.x,
                    translation.y,
                    translation.z,
                ]);
                arms[i].cross(direction).scale(sign)
            };
            let rotation = properties[i].angular_mobility_factor(torque);
            factor[j][6 * i + 3..6 * i + 6].copy_from_slice(&[rotation.x, rotation.y, rotation.z]);
        }
    }
    let mut result = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..=i {
            result[i][j] = finite(factor[i].iter().zip(factor[j]).map(|(a, b)| a * b).sum())?;
            result[j][i] = result[i][j];
        }
    }
    Ok(result)
}

fn projected_rates(
    points: [PointKinematics; 2],
    axes: [Vec3; 3],
) -> Result<[f64; 4], RigidJkrCombinedError> {
    let point = relative_point(points)?;
    Ok([
        finite(point.dot(axes[0]))?,
        finite(point.dot(axes[1]))?,
        finite(point.dot(axes[2]))?,
        finite(relative_spin(points).dot(axes[0]))?,
    ])
}

/// Stage both native bodies and the unchanged original constitutive association.
#[allow(clippy::too_many_arguments)]
pub fn prepare_rigid_jkr_combined(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidJkrCombinedContact<'_>,
    previous: RigidJkrCombinedState,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidJkrCombinedControls,
    gate: &CancelGate,
) -> Result<RigidJkrCombinedCandidate, RigidJkrCombinedError> {
    PreparedRigidJkrCombinedPair::new(bodies, properties, contact, controls, gate)?
        .advance_with_wrenches(previous, wrenches, dt_s, gate)
}

/// Publish both momenta and retained four-coordinate history after the final gate.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_jkr_combined(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidJkrCombinedContact<'_>,
    history: &mut RigidJkrCombinedState,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidJkrCombinedControls,
    gate: &CancelGate,
) -> Result<RigidJkrCombinedReport, RigidJkrCombinedError> {
    let candidate = prepare_rigid_jkr_combined(
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
fn finite(value: f64) -> Result<f64, RigidJkrCombinedError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RigidJkrCombinedError::Invalid(
            "unrepresentable native combined JKR mechanics",
        ))
    }
}
