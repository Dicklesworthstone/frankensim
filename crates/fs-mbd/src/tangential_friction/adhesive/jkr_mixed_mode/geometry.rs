//! Compatible pose endpoints with the original elastic history retained.
//!
//! This handoff consumes the accepted phase's sealed COM travel. It does not
//! integrate attitude, relocate material markers, repair a gap or supply moving
//! contact work. All such changes must already be compatible with the original
//! body's momentum, energy and normal approach. The next phase then uses the
//! actual new full mobility and rates, never a restarted constitutive history.

use super::{
    NativeEndpoint, PreparedRigidJkrMixedModePair, RigidBodyState, RigidJkrMixedModeContact,
    RigidJkrMixedModeError, RigidJkrMixedModeState, Vec3, capture_geometry, finite, magnitude,
    points, poll, relative_point, relative_spin, rounded_components,
};
use fs_contact::normal_patch::adhesive::mixed_mode::motion::{
    JkrMixedModeMotionError, transport::JkrMixedModeProjectionReport,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_tribo::{
    partial_slip::TangentFrame,
    tangential::transport::{TangentialHistoryRotation, TangentialTransportError},
};

/// Physical endpoint defects, not impulse, work, heat or position corrections.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrGeometryAudit {
    pub com_displacement_residual_m: [Vec3; 2],
    pub linear_momentum_residual_ns: [Vec3; 2],
    pub world_spin_momentum_residual_nms: [Vec3; 2],
    /// Includes the orbital change caused by moving each actual COM.
    pub angular_momentum_residual_nms: Vec3,
    pub kinetic_energy_residual_j: f64,
    pub normal_approach_residual_m: f64,
    pub relative_spin_rad_s: [f64; 2],
}

/// Same original source storage, with an actual new native kinematic projection.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrGeometryReport {
    pub source: JkrMixedModeProjectionReport,
    pub audit: RigidJkrGeometryAudit,
    pub mobility_before_per_kg: [[f64; 3]; 3],
    pub mobility_after_per_kg: [[f64; 3]; 3],
    pub relative_velocity_m_s: [[f64; 3]; 2],
}

/// Replacement captured pair and retained state, published together by caller.
pub struct RigidJkrGeometryContinuation<'a> {
    previous: RigidJkrMixedModeState,
    pair: PreparedRigidJkrMixedModePair<'a>,
    state: RigidJkrMixedModeState,
    report: RigidJkrGeometryReport,
}
impl<'a> RigidJkrGeometryContinuation<'a> {
    pub const fn previous_state(&self) -> RigidJkrMixedModeState {
        self.previous
    }
    pub const fn pair(&self) -> &PreparedRigidJkrMixedModePair<'a> {
        &self.pair
    }
    pub const fn state(&self) -> RigidJkrMixedModeState {
        self.state
    }
    pub const fn bodies(&self) -> &[RigidBodyState; 2] {
        &self.pair.before
    }
    pub const fn report(&self) -> &RigidJkrGeometryReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        PreparedRigidJkrMixedModePair<'a>,
        RigidJkrMixedModeState,
        RigidJkrGeometryReport,
    ) {
        (self.pair, self.state, self.report)
    }
}

#[derive(Debug)]
pub enum RigidJkrGeometryError {
    Native(RigidJkrMixedModeError),
    Frame(TangentialTransportError),
    Audit(Box<RigidJkrGeometryAudit>),
    StorageBalance { residual_j: f64, tolerance_j: f64 },
}
impl From<RigidJkrMixedModeError> for RigidJkrGeometryError {
    fn from(error: RigidJkrMixedModeError) -> Self {
        Self::Native(error)
    }
}
impl From<crate::DynamicsError> for RigidJkrGeometryError {
    fn from(error: crate::DynamicsError) -> Self {
        Self::Native(error.into())
    }
}
impl From<super::RigidFrictionError> for RigidJkrGeometryError {
    fn from(error: super::RigidFrictionError) -> Self {
        Self::Native(error.into())
    }
}
impl From<JkrMixedModeMotionError> for RigidJkrGeometryError {
    fn from(error: JkrMixedModeMotionError) -> Self {
        Self::Native(error.into())
    }
}
impl From<TangentialTransportError> for RigidJkrGeometryError {
    fn from(error: TangentialTransportError) -> Self {
        Self::Frame(error)
    }
}
impl core::fmt::Display for RigidJkrGeometryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native full-stick JKR geometry handoff: {self:?}")
    }
}
impl std::error::Error for RigidJkrGeometryError {}

impl<'a> PreparedRigidJkrMixedModePair<'a> {
    /// Accept caller-owned compatible poses, then retain the exact source clock.
    ///
    /// `endpoint_bodies` retain original A/B order, mass, body coordinates and
    /// momentum ownership. COM travel must equal the accepted phase travel;
    /// linear momentum, world spin, total orbital-plus-spin angular momentum
    /// and full kinetic energy must close under the original native gates.
    /// The original material markers and normal datum remain fixed. Actual
    /// normal approach must match the retained source coordinate; lateral
    /// marker separation, unowned corrections and relative torsion refuse.
    ///
    /// The existing tangent isometry transports elastic shear by minimum normal
    /// rotation plus the explicitly supplied common material spin. A tangent
    /// gauge change is not spin. Endpoint data do not resolve a general path
    /// of noncommuting finite rotations; the pose owner must refine that path.
    /// No pose integrator, sliding, marker migration, external impulse, work,
    /// heat or elapsed time is introduced here. Projected kinetic repartition
    /// is diagnostic only; actual native kinetic energy controls acceptance.
    pub fn continue_at_geometry(
        &self,
        previous: RigidJkrMixedModeState,
        endpoint_bodies: [RigidBodyState; 2],
        endpoint_frame: TangentFrame,
        common_spin_angle_rad: f64,
        gate: &CancelGate,
    ) -> Result<RigidJkrGeometryContinuation<'a>, RigidJkrGeometryError> {
        poll(gate)?;
        let anchor = previous.endpoint;
        if previous.frame != self.contact.frame
            || previous.marker_arms != self.marker_arms
            || previous.reference_approach_m != self.contact.reference_approach_m
            || anchor.properties != self.properties
            || anchor.body_indices != [self.contact.binding.body_a, self.contact.binding.body_b]
            || (0..2).any(|i| anchor.bodies[i].pose() != self.before[i].pose())
        {
            return Err(RigidJkrMixedModeError::StateMismatch.into());
        }
        if (0..2).any(|i| {
            !anchor.pending_angular_motion[i]
                && endpoint_bodies[i].pose().orientation() != anchor.bodies[i].pose().orientation()
        }) || (!anchor.pending_angular_motion.into_iter().any(|moved| moved)
            && common_spin_angle_rad != 0.0)
        {
            return Err(RigidJkrMixedModeError::Invalid(
                "endpoint attitude or common spin changed without accepted angular motion",
            )
            .into());
        }
        // A post-phase source coordinate need not equal the old held geometry.
        // Its actual endpoint rates come from the sealed accepted momenta.
        let old_points = points(anchor.bodies, self.properties, self.reaction_arms)
            .map_err(RigidJkrMixedModeError::from)?;
        let old_relative = relative_point(old_points)?;
        let old_rates = self.axes.map(|axis| old_relative.dot(axis));
        if (0..3).any(|i| {
            !old_rates[i].is_finite()
                || (old_rates[i] - previous.velocity_m_s()[i]).abs()
                    > self.controls.velocity_tolerance_m_s
        }) {
            return Err(RigidJkrMixedModeError::StateMismatch.into());
        }
        let old_spin = finite(relative_spin(old_points).dot(self.axes[0]))?;
        let contact = RigidJkrMixedModeContact {
            frame: endpoint_frame,
            ..self.contact
        };
        let geometry = capture_geometry(endpoint_bodies, self.properties, contact, self.controls)?;
        let new_relative = relative_point(geometry.reaction_before)?;
        let rates = geometry.axes.map(|axis| new_relative.dot(axis));
        let new_spin = finite(relative_spin(geometry.reaction_before).dot(geometry.axes[0]))?;
        self.admit_spin([old_spin, new_spin])?;

        let mut audit = RigidJkrGeometryAudit {
            com_displacement_residual_m: [Vec3::ZERO; 2],
            linear_momentum_residual_ns: [Vec3::ZERO; 2],
            world_spin_momentum_residual_nms: [Vec3::ZERO; 2],
            angular_momentum_residual_nms: Vec3::ZERO,
            kinetic_energy_residual_j: 0.0,
            normal_approach_residual_m: geometry.approach_m - previous.position_m()[0],
            relative_spin_rad_s: [old_spin, new_spin],
        };
        let mut angular = [Dd::from_f64(0.0); 3];
        let mut kinetic = Dd::from_f64(0.0);
        for i in 0..2 {
            let old = anchor.bodies[i];
            let new = endpoint_bodies[i];
            let x0 = components(old.pose().position_world());
            let x1 = components(new.pose().position_world());
            let travel = components(anchor.pending_com_travel_m[i]);
            audit.com_displacement_residual_m[i] = rounded_components(core::array::from_fn(|j| {
                Dd::from_f64(x1[j]) - Dd::from_f64(x0[j]) - Dd::from_f64(travel[j])
            }));
            audit.linear_momentum_residual_ns[i] =
                new.linear_momentum_world().sub(old.linear_momentum_world());
            let spin0 = old
                .pose()
                .orientation()
                .rotate_body_to_world(old.angular_momentum_body());
            let spin1 = new
                .pose()
                .orientation()
                .rotate_body_to_world(new.angular_momentum_body());
            audit.world_spin_momentum_residual_nms[i] = spin1.sub(spin0);
            let orbital0 = cross_dd(old.pose().position_world(), old.linear_momentum_world());
            let orbital1 = cross_dd(new.pose().position_world(), new.linear_momentum_world());
            let spin0 = components(spin0);
            let spin1 = components(spin1);
            for j in 0..3 {
                angular[j] = angular[j] + orbital1[j] - orbital0[j] + Dd::from_f64(spin1[j])
                    - Dd::from_f64(spin0[j]);
            }
            kinetic = kinetic + Dd::from_f64(new.kinetic_energy(self.properties[i])?)
                - Dd::from_f64(old.kinetic_energy(self.properties[i])?);
        }
        audit.angular_momentum_residual_nms = rounded_components(angular);
        audit.kinetic_energy_residual_j = kinetic.to_f64();
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
        {
            return Err(RigidJkrGeometryError::Audit(Box::new(audit)));
        }
        let rotation = TangentialHistoryRotation::between(
            previous.frame,
            endpoint_frame,
            common_spin_angle_rad,
        )?;
        let transfer = self.motion.rebind_native_projection(
            previous.motion,
            rotation,
            geometry.mobility,
            rates,
            gate,
        )?;
        let (motion, state, source) = transfer.into_parts();
        if !source.source_energy_residual_j.is_finite()
            || source.source_energy_residual_j.abs() > c.energy_tolerance_j
        {
            return Err(RigidJkrGeometryError::StorageBalance {
                residual_j: source.source_energy_residual_j,
                tolerance_j: c.energy_tolerance_j,
            });
        }
        let pair = PreparedRigidJkrMixedModePair {
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
        poll(gate)?;
        Ok(RigidJkrGeometryContinuation {
            previous,
            pair,
            state: RigidJkrMixedModeState {
                motion: state,
                frame: endpoint_frame,
                endpoint: NativeEndpoint {
                    bodies: endpoint_bodies,
                    pending_com_travel_m: [Vec3::ZERO; 2],
                    pending_angular_motion: [false; 2],
                    ..anchor
                },
                ..previous
            },
            report: RigidJkrGeometryReport {
                source,
                audit,
                mobility_before_per_kg: self.mobility,
                mobility_after_per_kg: geometry.mobility,
                relative_velocity_m_s: [old_rates, rates],
            },
        })
    }
}

fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn cross_dd(a: Vec3, b: Vec3) -> [Dd; 3] {
    let a = components(a).map(Dd::from_f64);
    let b = components(b).map(Dd::from_f64);
    core::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
