//! Accepted native travel, commuting axial attitudes and objective elastic shear.
//!
//! The fixed original normal and collinear material markers delimit this owner.
//! Every accepted angular phase must have remained an individual principal-axis
//! spin about that normal. The sealed integrated angles then determine both
//! actual attitudes and common material spin. Generic finite rotations, contact
//! migration and a changed normal require their own resolved geometry path.

use super::{
    NativeEndpoint, PointKinematics, PreparedRigidJkrCombinedPair, RigidBodyState,
    RigidJkrCombinedContact, RigidJkrCombinedError, RigidJkrCombinedState, Vec3, capture_geometry,
    finite, magnitude, points, poll, projected_rates, rounded_components,
};
use crate::{Pose, UnitQuaternion};
use fs_contact::normal_patch::adhesive::mixed_mode::combined::motion::{
    JkrCombinedMotionError, transport::JkrCombinedProjectionReport,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_tribo::{
    partial_slip::TangentFrame,
    tangential::transport::{TangentialHistoryRotation, TangentialTransportError},
};

/// Independent body defects. These are not impulses, work or pose corrections.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrCombinedGeometryAudit {
    pub com_displacement_residual_m: [Vec3; 2],
    /// Maximum quaternion-component discrepancy from the sealed axial path.
    /// The double-cover sign is immaterial; only representation roundoff passes.
    pub attitude_residual: [f64; 2],
    pub linear_momentum_residual_ns: [Vec3; 2],
    pub world_spin_momentum_residual_nms: [Vec3; 2],
    pub angular_momentum_residual_nms: Vec3,
    pub kinetic_energy_residual_j: f64,
    pub normal_approach_residual_m: f64,
}

/// Source storage is unchanged; native energy independently owns acceptance.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrCombinedGeometryReport {
    pub source: JkrCombinedProjectionReport,
    pub audit: RigidJkrCombinedGeometryAudit,
    pub mobility_before: [[f64; 4]; 4],
    pub mobility_after: [[f64; 4]; 4],
    pub rates: [[f64; 4]; 2],
    pub consumed_angular_displacement_body_rad: [Vec3; 2],
    /// Mean of the two retained signed axial angles, not a tangent gauge angle.
    pub common_spin_angle_rad: f64,
}

/// Original history, actual endpoint bodies and new projection adopted together.
pub struct RigidJkrCombinedGeometryContinuation<'a> {
    previous: RigidJkrCombinedState,
    pair: PreparedRigidJkrCombinedPair<'a>,
    state: RigidJkrCombinedState,
    report: RigidJkrCombinedGeometryReport,
}
impl<'a> RigidJkrCombinedGeometryContinuation<'a> {
    pub const fn previous_state(&self) -> RigidJkrCombinedState {
        self.previous
    }
    pub const fn pair(&self) -> &PreparedRigidJkrCombinedPair<'a> {
        &self.pair
    }
    pub const fn state(&self) -> RigidJkrCombinedState {
        self.state
    }
    pub const fn bodies(&self) -> &[RigidBodyState; 2] {
        &self.pair.before
    }
    pub const fn report(&self) -> &RigidJkrCombinedGeometryReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        PreparedRigidJkrCombinedPair<'a>,
        RigidJkrCombinedState,
        RigidJkrCombinedGeometryReport,
    ) {
        (self.pair, self.state, self.report)
    }
}

#[derive(Debug)]
pub enum RigidJkrCombinedGeometryError {
    Native(RigidJkrCombinedError),
    Frame(TangentialTransportError),
    /// Integrated angular vectors do not resolve the requested finite path.
    NoncommutingMotion,
    /// Only a tangent gauge change is supplied by this fixed-normal owner.
    ChangedNormal,
    Audit(Box<RigidJkrCombinedGeometryAudit>),
}
impl From<RigidJkrCombinedError> for RigidJkrCombinedGeometryError {
    fn from(value: RigidJkrCombinedError) -> Self {
        Self::Native(value)
    }
}
impl From<crate::DynamicsError> for RigidJkrCombinedGeometryError {
    fn from(value: crate::DynamicsError) -> Self {
        Self::Native(value.into())
    }
}
impl From<super::RigidFrictionError> for RigidJkrCombinedGeometryError {
    fn from(value: super::RigidFrictionError) -> Self {
        Self::Native(value.into())
    }
}
impl From<JkrCombinedMotionError> for RigidJkrCombinedGeometryError {
    fn from(value: JkrCombinedMotionError) -> Self {
        Self::Native(value.into())
    }
}
impl From<TangentialTransportError> for RigidJkrCombinedGeometryError {
    fn from(value: TangentialTransportError) -> Self {
        Self::Frame(value)
    }
}
impl core::fmt::Display for RigidJkrCombinedGeometryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "native combined JKR geometry continuation refused: {self:?}"
        )
    }
}
impl std::error::Error for RigidJkrCombinedGeometryError {}

impl<'a> PreparedRigidJkrCombinedPair<'a> {
    /// Turn the accepted COM travel and commuting axial spin into actual poses.
    /// The fixed world tangent gauge is retained; material shear follows the
    /// mean original body spin through the existing objective-history rotation.
    pub fn continue_axial_geometry(
        &self,
        previous: RigidJkrCombinedState,
        gate: &CancelGate,
    ) -> Result<RigidJkrCombinedGeometryContinuation<'a>, RigidJkrCombinedGeometryError> {
        poll(gate)?;
        let anchor = previous.endpoint;
        if !anchor.axial_drift_admissible {
            return Err(RigidJkrCombinedGeometryError::NoncommutingMotion);
        }
        let mut endpoint_bodies = anchor.bodies;
        for (i, body) in endpoint_bodies.iter_mut().enumerate() {
            poll(gate)?;
            let old = anchor.bodies[i];
            *body = RigidBodyState::new(
                Pose::new(
                    add_travel(old.pose().position_world(), anchor.pending_com_travel_m[i]),
                    axial_attitude(
                        old.pose().orientation(),
                        anchor.pending_angular_displacement_body_rad[i],
                    )?,
                )?,
                old.linear_momentum_world(),
                old.angular_momentum_body(),
            )?;
        }
        self.continue_at_geometry(previous, endpoint_bodies, self.contact.frame, gate)
    }

    /// Admit compatible caller-owned endpoint poses and an optional tangent gauge.
    ///
    /// The body pair is ordered A/B regardless of its original assembly indices.
    /// COM travel must consume the sealed accepted increments. Both attitudes
    /// must match their retained commuting axial paths; numerical audit tolerances
    /// cannot authorize a different finite rotation. Common material spin is
    /// derived from those angles, never inferred from the chosen tangent gauge.
    ///
    /// The normal, material markers and datum remain fixed. Transverse marker
    /// separation, unrepresented rolling, an unresolved orbital momentum change
    /// or a mismatch with the accepted source approach refuses. Full native
    /// momenta, kinetic energy and all four source stores are independently
    /// checked. No impulse, gap correction, marker migration, heat, work or time
    /// is added. The next momentum phase uses the actual new complete mobility
    /// and rates with the same source chronology and elastic shear/twist history.
    pub fn continue_at_geometry(
        &self,
        previous: RigidJkrCombinedState,
        endpoint_bodies: [RigidBodyState; 2],
        endpoint_frame: TangentFrame,
        gate: &CancelGate,
    ) -> Result<RigidJkrCombinedGeometryContinuation<'a>, RigidJkrCombinedGeometryError> {
        poll(gate)?;
        let anchor = previous.endpoint;
        if previous.frame != self.contact.frame
            || previous.marker_arms != self.marker_arms
            || previous.reference_approach_m != self.contact.reference_approach_m
            || anchor.properties != self.properties
            || anchor.body_indices != [self.contact.binding.body_a, self.contact.binding.body_b]
            || (0..2).any(|i| anchor.bodies[i].pose() != self.before[i].pose())
        {
            return Err(RigidJkrCombinedError::StateMismatch.into());
        }
        if endpoint_frame.normal() != self.contact.frame.normal() {
            return Err(RigidJkrCombinedGeometryError::ChangedNormal);
        }
        let old_points = points(anchor.bodies, self.properties, self.reaction_arms)
            .map_err(RigidJkrCombinedError::from)?;
        self.admit_unrepresented(old_points)?;
        if !anchor.axial_drift_admissible
            || !axial_drift_admissible(anchor.bodies, old_points, self.axes[0])
        {
            return Err(RigidJkrCombinedGeometryError::NoncommutingMotion);
        }
        let old_rates = projected_rates(old_points, self.axes)?;
        let rate_limits = [
            self.controls.velocity_tolerance_m_s,
            self.controls.velocity_tolerance_m_s,
            self.controls.velocity_tolerance_m_s,
            self.controls.angular_velocity_tolerance_rad_s,
        ];
        if (0..4).any(|i| (old_rates[i] - previous.rates()[i]).abs() > rate_limits[i]) {
            return Err(RigidJkrCombinedError::StateMismatch.into());
        }
        let mut audit = RigidJkrCombinedGeometryAudit {
            com_displacement_residual_m: [Vec3::ZERO; 2],
            attitude_residual: [0.0; 2],
            linear_momentum_residual_ns: [Vec3::ZERO; 2],
            world_spin_momentum_residual_nms: [Vec3::ZERO; 2],
            angular_momentum_residual_nms: Vec3::ZERO,
            kinetic_energy_residual_j: 0.0,
            normal_approach_residual_m: 0.0,
        };
        let mut axial_angle = [0.0; 2];
        let mut angular = [Dd::ZERO; 3];
        let mut kinetic = Dd::ZERO;
        for i in 0..2 {
            poll(gate)?;
            let old = anchor.bodies[i];
            let new = endpoint_bodies[i];
            let rotation = anchor.pending_angular_displacement_body_rad[i];
            let world_rotation = old.pose().orientation().rotate_body_to_world(rotation);
            if !axial(world_rotation, self.axes[0]) {
                return Err(RigidJkrCombinedGeometryError::NoncommutingMotion);
            }
            axial_angle[i] = finite(world_rotation.dot(self.axes[0]))?;
            audit.attitude_residual[i] = attitude_difference(
                new.pose().orientation(),
                axial_attitude(old.pose().orientation(), rotation)?,
            );
            let x0 = components(old.pose().position_world());
            let x1 = components(new.pose().position_world());
            let travel = components(anchor.pending_com_travel_m[i]);
            let displacement: [Dd; 3] =
                core::array::from_fn(|j| Dd::from_f64(x1[j]) - Dd::from_f64(x0[j]));
            audit.com_displacement_residual_m[i] = rounded_components(core::array::from_fn(|j| {
                displacement[j] - Dd::from_f64(travel[j])
            }));
            audit.linear_momentum_residual_ns[i] =
                new.linear_momentum_world().sub(old.linear_momentum_world());
            let spin0 = components(
                old.pose()
                    .orientation()
                    .rotate_body_to_world(old.angular_momentum_body()),
            );
            let spin1 = components(
                new.pose()
                    .orientation()
                    .rotate_body_to_world(new.angular_momentum_body()),
            );
            let spin: [Dd; 3] =
                core::array::from_fn(|j| Dd::from_f64(spin1[j]) - Dd::from_f64(spin0[j]));
            audit.world_spin_momentum_residual_nms[i] = rounded_components(spin);
            // Include any supplied momentum discrepancy in the orbital audit;
            // unchanged momenta reduce this to the accurate displacement cross p.
            let orbital = cross_displacement(displacement, old.linear_momentum_world());
            let changed_momentum = cross_displacement(
                components(new.pose().position_world()).map(Dd::from_f64),
                audit.linear_momentum_residual_ns[i],
            );
            for j in 0..3 {
                angular[j] = angular[j] + orbital[j] + changed_momentum[j] + spin[j];
            }
            kinetic =
                kinetic + Dd::from_f64(old.kinetic_energy_change_to(self.properties[i], new)?);
        }
        audit.angular_momentum_residual_nms = rounded_components(angular);
        audit.kinetic_energy_residual_j = finite(kinetic.to_f64())?;
        let contact = RigidJkrCombinedContact {
            frame: endpoint_frame,
            ..self.contact
        };
        let geometry = capture_geometry(endpoint_bodies, self.properties, contact, self.controls)?;
        self.admit_unrepresented(geometry.reaction_before)?;
        let rates = projected_rates(geometry.reaction_before, geometry.axes)?;
        audit.normal_approach_residual_m = geometry.approach_m - previous.coordinates()[0];
        let c = self.controls;
        if audit
            .com_displacement_residual_m
            .into_iter()
            .any(|v| !v.is_finite() || magnitude(v) > c.displacement_tolerance_m)
            || audit
                .attitude_residual
                .into_iter()
                .any(|v| !v.is_finite() || v > 128.0 * f64::EPSILON)
            || audit
                .linear_momentum_residual_ns
                .into_iter()
                .any(|v| !v.is_finite() || magnitude(v) > c.linear_momentum_tolerance_ns)
            || audit
                .world_spin_momentum_residual_nms
                .into_iter()
                .chain([audit.angular_momentum_residual_nms])
                .any(|v| !v.is_finite() || magnitude(v) > c.angular_momentum_tolerance_nms)
            || !audit.kinetic_energy_residual_j.is_finite()
            || audit.kinetic_energy_residual_j.abs() > c.energy_tolerance_j
            || !audit.normal_approach_residual_m.is_finite()
            || audit.normal_approach_residual_m.abs() > c.geometry_tolerance_m
        {
            return Err(RigidJkrCombinedGeometryError::Audit(Box::new(audit)));
        }
        let common_spin_angle_rad = finite(axial_angle[0].midpoint(axial_angle[1]))?;
        let rotation = TangentialHistoryRotation::between(
            previous.frame,
            endpoint_frame,
            common_spin_angle_rad,
        )?;
        let transfer = self.motion.rebind_projection(
            previous.motion,
            rotation,
            geometry.mobility,
            rates,
            gate,
        )?;
        let (motion, state, source) = transfer.into_parts();
        if !source.source_energy_residual_j.is_finite()
            || !source.absolute_storage_residual_j.is_finite()
            || source.source_energy_residual_j.abs() > c.energy_tolerance_j
            || source.absolute_storage_residual_j > c.energy_tolerance_j
        {
            return Err(JkrCombinedMotionError::EnergyBalance {
                residual_j: source
                    .source_energy_residual_j
                    .abs()
                    .max(source.absolute_storage_residual_j),
                tolerance_j: c.energy_tolerance_j,
            }
            .into());
        }
        let endpoint = NativeEndpoint {
            bodies: endpoint_bodies,
            pending_com_travel_m: [Vec3::ZERO; 2],
            pending_angular_displacement_body_rad: [Vec3::ZERO; 2],
            axial_drift_admissible: axial_drift_admissible(
                endpoint_bodies,
                geometry.reaction_before,
                geometry.axes[0],
            ),
            ..anchor
        };
        let pair = PreparedRigidJkrCombinedPair {
            contact,
            before: endpoint_bodies,
            properties: self.properties,
            marker_arms: geometry.marker_arms,
            reaction_arms: geometry.reaction_arms,
            axes: geometry.axes,
            markers_before: geometry.markers_before,
            reaction_before: geometry.reaction_before,
            approach_m: geometry.approach_m,
            mobility: geometry.mobility,
            motion,
            controls: c,
        };
        let state = RigidJkrCombinedState {
            motion: state,
            frame: endpoint_frame,
            endpoint,
            ..previous
        };
        let report = RigidJkrCombinedGeometryReport {
            source,
            audit,
            mobility_before: self.mobility,
            mobility_after: geometry.mobility,
            rates: [old_rates, rates],
            consumed_angular_displacement_body_rad: anchor.pending_angular_displacement_body_rad,
            common_spin_angle_rad,
        };
        poll(gate)?;
        Ok(RigidJkrCombinedGeometryContinuation {
            previous,
            pair,
            state,
            report,
        })
    }
}

/// Held momenta and angular rates are affine, so both phase endpoints determine
/// whether the complete phase stayed on the same individual principal axis.
pub(super) fn axial_drift_admissible(
    bodies: [RigidBodyState; 2],
    points: [PointKinematics; 2],
    normal: Vec3,
) -> bool {
    (0..2).all(|i| {
        axial(points[i].angular_velocity_world, normal)
            && axial(
                bodies[i]
                    .pose()
                    .orientation()
                    .rotate_body_to_world(bodies[i].angular_momentum_body()),
                normal,
            )
    })
}

fn axial(value: Vec3, normal: Vec3) -> bool {
    let scale = magnitude(value);
    let transverse = value.sub(normal.scale(value.dot(normal)));
    value.is_finite()
        && scale.is_finite()
        && transverse.is_finite()
        && magnitude(transverse) <= 128.0 * f64::EPSILON * scale
}
fn axial_attitude(
    orientation: UnitQuaternion,
    rotation: Vec3,
) -> Result<UnitQuaternion, crate::DynamicsError> {
    if rotation == Vec3::ZERO {
        Ok(orientation)
    } else {
        orientation.right_exp(rotation)
    }
}
fn attitude_difference(a: UnitQuaternion, b: UnitQuaternion) -> f64 {
    let a = a.components();
    let b = b.components();
    let same = (0..4).fold(0.0_f64, |r, i| r.max((a[i] - b[i]).abs()));
    let opposite = (0..4).fold(0.0_f64, |r, i| r.max((a[i] + b[i]).abs()));
    same.min(opposite)
}
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn add_travel(position: Vec3, travel: Vec3) -> Vec3 {
    rounded_components(core::array::from_fn(|i| {
        Dd::from_f64(components(position)[i]) + Dd::from_f64(components(travel)[i])
    }))
}
fn cross_displacement(d: [Dd; 3], momentum: Vec3) -> [Dd; 3] {
    let p = components(momentum).map(Dd::from_f64);
    [
        d[1] * p[2] - d[2] * p[1],
        d[2] * p[0] - d[0] * p[2],
        d[0] * p[1] - d[1] * p[0],
    ]
}
