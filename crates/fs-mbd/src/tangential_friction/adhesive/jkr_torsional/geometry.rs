//! Accepted axial spin and normal travel become actual native poses.
//!
//! This continuation is deliberately restricted to commuting spin about the
//! original fixed normal. It consumes sealed phase travel, retains the source
//! coordinates and clock, and rebinds to the actual new native mobility.
//! Arbitrary finite rotations require a separate resolved attitude path.

use super::{
    NativeEndpoint, PreparedRigidJkrTorsionalPair, RigidBodyState, RigidJkrTorsionalError,
    RigidJkrTorsionalState, Vec3, axial_drift_admissible, capture_geometry, finite, magnitude,
    points, poll, relative_point, relative_spin, rounded_components,
};
use crate::Pose;
use fs_contact::normal_patch::adhesive::mixed_mode::torsion::motion::{
    JkrTorsionalMotionError, JkrTorsionalProjectionReport,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// Actual geometry defects. No term authorizes an impulse or an energy repair.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrTorsionalGeometryAudit {
    pub com_displacement_residual_m: [Vec3; 2],
    pub linear_momentum_residual_ns: [Vec3; 2],
    pub world_spin_momentum_residual_nms: [Vec3; 2],
    pub angular_momentum_residual_nms: Vec3,
    pub kinetic_energy_residual_j: f64,
    pub normal_approach_residual_m: f64,
    /// Normal velocity [m/s] and normal-axis relative spin [rad/s].
    pub rate_residual: [f64; 2],
}

/// Repartition of reduced coordinates is diagnostic; native energy is unchanged.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrTorsionalGeometryReport {
    pub source: JkrTorsionalProjectionReport,
    pub audit: RigidJkrTorsionalGeometryAudit,
    pub mobility_before: [[f64; 2]; 2],
    pub mobility_after: [[f64; 2]; 2],
    pub consumed_angular_displacement_body_rad: [Vec3; 2],
}

/// New poses, source chronology and native projection form one private result.
pub struct RigidJkrTorsionalGeometryContinuation<'a> {
    previous: RigidJkrTorsionalState,
    pair: PreparedRigidJkrTorsionalPair<'a>,
    state: RigidJkrTorsionalState,
    report: RigidJkrTorsionalGeometryReport,
}
impl<'a> RigidJkrTorsionalGeometryContinuation<'a> {
    pub const fn previous_state(&self) -> RigidJkrTorsionalState {
        self.previous
    }
    pub const fn pair(&self) -> &PreparedRigidJkrTorsionalPair<'a> {
        &self.pair
    }
    pub const fn state(&self) -> RigidJkrTorsionalState {
        self.state
    }
    pub const fn bodies(&self) -> &[RigidBodyState; 2] {
        &self.pair.before
    }
    pub const fn report(&self) -> &RigidJkrTorsionalGeometryReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        PreparedRigidJkrTorsionalPair<'a>,
        RigidJkrTorsionalState,
        RigidJkrTorsionalGeometryReport,
    ) {
        (self.pair, self.state, self.report)
    }
}

#[derive(Debug)]
pub enum RigidJkrTorsionalGeometryError {
    Native(RigidJkrTorsionalError),
    /// A scalar twist and an integrated angular vector do not resolve this path.
    NoncommutingMotion,
    Audit(Box<RigidJkrTorsionalGeometryAudit>),
}
impl From<RigidJkrTorsionalError> for RigidJkrTorsionalGeometryError {
    fn from(value: RigidJkrTorsionalError) -> Self {
        Self::Native(value)
    }
}
impl From<crate::DynamicsError> for RigidJkrTorsionalGeometryError {
    fn from(value: crate::DynamicsError) -> Self {
        Self::Native(value.into())
    }
}
impl From<super::RigidFrictionError> for RigidJkrTorsionalGeometryError {
    fn from(value: super::RigidFrictionError) -> Self {
        Self::Native(value.into())
    }
}
impl From<JkrTorsionalMotionError> for RigidJkrTorsionalGeometryError {
    fn from(value: JkrTorsionalMotionError) -> Self {
        Self::Native(value.into())
    }
}
impl core::fmt::Display for RigidJkrTorsionalGeometryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "native normal/twist geometry continuation refused: {self:?}"
        )
    }
}
impl std::error::Error for RigidJkrTorsionalGeometryError {}

impl<'a> PreparedRigidJkrTorsionalPair<'a> {
    /// Consume accepted COM travel and commuting principal-axis spin.
    ///
    /// The actual source approach must equal the new marker separation. The
    /// original normal, marker association and datum stay fixed. Every pending
    /// phase must have had both individual angular velocity and spin momentum
    /// along that normal; cancelling a transverse excursion cannot admit it.
    ///
    /// The exponential map integrates the retained commuting spin exactly.
    /// It does not reinterpret an arbitrary integrated angular vector as a
    /// finite attitude. Linear momentum, world spin, total angular momentum
    /// and native kinetic energy must independently survive the pose change.
    /// Unresolved orbital drift refuses and calls for a refined physical phase.
    ///
    /// This handoff advances no source time and adds no work, heat, impulse,
    /// marker migration or gap correction. Both caller-owned bodies and
    /// accepted history remain untouched until the returned result is adopted.
    pub fn continue_axial_geometry(
        &self,
        previous: RigidJkrTorsionalState,
        gate: &CancelGate,
    ) -> Result<RigidJkrTorsionalGeometryContinuation<'a>, RigidJkrTorsionalGeometryError> {
        poll(gate)?;
        let anchor = previous.endpoint;
        if previous.frame != self.contact.frame
            || previous.marker_arms != self.marker_arms
            || previous.reference_approach_m != self.contact.reference_approach_m
            || anchor.properties != self.properties
            || anchor.body_indices != [self.contact.binding.body_a, self.contact.binding.body_b]
            || (0..2).any(|i| anchor.bodies[i].pose() != self.before[i].pose())
        {
            return Err(RigidJkrTorsionalError::StateMismatch.into());
        }
        let old_points = points(anchor.bodies, self.properties, self.reaction_arms)
            .map_err(RigidJkrTorsionalError::from)?;
        self.admit_unrepresented(old_points)?;
        if !anchor.axial_drift_admissible
            || !axial_drift_admissible(anchor.bodies, old_points, self.axes[0], self.controls)
        {
            return Err(RigidJkrTorsionalGeometryError::NoncommutingMotion);
        }
        let old_rates = [
            relative_point(old_points)?.dot(self.axes[0]),
            relative_spin(old_points).dot(self.axes[0]),
        ];
        if (old_rates[0] - previous.rates()[0]).abs() > self.controls.velocity_tolerance_m_s
            || (old_rates[1] - previous.rates()[1]).abs()
                > self.controls.angular_velocity_tolerance_rad_s
        {
            return Err(RigidJkrTorsionalError::StateMismatch.into());
        }
        let mut endpoint_bodies = anchor.bodies;
        let mut audit = RigidJkrTorsionalGeometryAudit {
            com_displacement_residual_m: [Vec3::ZERO; 2],
            linear_momentum_residual_ns: [Vec3::ZERO; 2],
            world_spin_momentum_residual_nms: [Vec3::ZERO; 2],
            angular_momentum_residual_nms: Vec3::ZERO,
            kinetic_energy_residual_j: 0.0,
            normal_approach_residual_m: 0.0,
            rate_residual: [0.0; 2],
        };
        let mut angular = [Dd::ZERO; 3];
        let mut kinetic = Dd::ZERO;
        for i in 0..2 {
            poll(gate)?;
            let old = anchor.bodies[i];
            let position = add_travel(old.pose().position_world(), anchor.pending_com_travel_m[i]);
            let rotation = anchor.pending_angular_displacement_body_rad[i];
            let orientation = if rotation == Vec3::ZERO {
                old.pose().orientation()
            } else {
                old.pose().orientation().right_exp(rotation)?
            };
            let new = RigidBodyState::new(
                Pose::new(position, orientation)?,
                old.linear_momentum_world(),
                old.angular_momentum_body(),
            )?;
            endpoint_bodies[i] = new;
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
            let orbital = cross_displacement(displacement, old.linear_momentum_world());
            for j in 0..3 {
                angular[j] = angular[j] + orbital[j] + spin[j];
            }
            kinetic =
                kinetic + Dd::from_f64(old.kinetic_energy_change_to(self.properties[i], new)?);
        }
        audit.angular_momentum_residual_nms = rounded_components(angular);
        audit.kinetic_energy_residual_j = finite(kinetic.to_f64())?;

        let geometry = capture_geometry(
            endpoint_bodies,
            self.properties,
            self.contact,
            self.controls,
        )?;
        self.admit_unrepresented(geometry.reaction_before)?;
        let rates = [
            relative_point(geometry.reaction_before)?.dot(self.axes[0]),
            relative_spin(geometry.reaction_before).dot(self.axes[0]),
        ];
        audit.normal_approach_residual_m = geometry.approach_m - previous.coordinates()[0];
        audit.rate_residual = core::array::from_fn(|i| rates[i] - previous.rates()[i]);
        let c = self.controls;
        if audit
            .com_displacement_residual_m
            .into_iter()
            .any(|v| !v.is_finite() || magnitude(v) > c.displacement_tolerance_m)
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
            || !audit.rate_residual[0].is_finite()
            || audit.rate_residual[0].abs() > c.velocity_tolerance_m_s
            || !audit.rate_residual[1].is_finite()
            || audit.rate_residual[1].abs() > c.angular_velocity_tolerance_rad_s
        {
            return Err(RigidJkrTorsionalGeometryError::Audit(Box::new(audit)));
        }

        let transfer =
            self.motion
                .rebind_projection(previous.motion, geometry.mobility, rates, gate)?;
        let (motion, state, source) = transfer.into_parts();
        let endpoint = NativeEndpoint {
            bodies: endpoint_bodies,
            pending_com_travel_m: [Vec3::ZERO; 2],
            pending_angular_displacement_body_rad: [Vec3::ZERO; 2],
            axial_drift_admissible: axial_drift_admissible(
                endpoint_bodies,
                geometry.reaction_before,
                self.axes[0],
                self.controls,
            ),
            ..anchor
        };
        let pair = PreparedRigidJkrTorsionalPair {
            contact: self.contact,
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
            controls: self.controls,
        };
        let state = RigidJkrTorsionalState {
            motion: state,
            endpoint,
            ..previous
        };
        let report = RigidJkrTorsionalGeometryReport {
            source,
            audit,
            mobility_before: self.mobility,
            mobility_after: geometry.mobility,
            consumed_angular_displacement_body_rad: anchor.pending_angular_displacement_body_rad,
        };
        poll(gate)?;
        Ok(RigidJkrTorsionalGeometryContinuation {
            previous,
            pair,
            state,
            report,
        })
    }
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
