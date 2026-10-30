//! Ambient loads at the actual moving pose, with immutable reference inventory.
//!
//! Body forces already participate in the native potential and force assembly.
//! This view supplies their external impulse account without mistaking them for
//! internal elastic forces, or replacing body-frame attachment arms by a solid
//! mass centroid when fluid and solid densities differ.

use std::sync::Arc;

use fs_exec::CancelGate;
use fs_ga::Vec3;
use fs_material::fluid_transport::FluidTransport;

use super::{PreparedRodMotion, RodMotionPose, SpatialRodError, checkpoint};

/// External conservative body loads evaluated at one admitted native pose.
#[derive(Clone, Debug, PartialEq)]
pub struct RodBodyLoads {
    /// Solid weight and displaced-fluid force at each material node [N], world axes.
    pub force_n: Vec<[f64; 3]>,
    /// Pure body-load couple about each material node [N m], body axes.
    pub torque_body_nm: Vec<[f64; 3]>,
    /// Solid gravitational potential [J], relative to the original reference pose.
    pub gravity_potential_j: f64,
    /// Hydrostatic potential [J], using the same original reference pose.
    pub buoyancy_potential_j: f64,
}

impl<I> PreparedRodMotion<I> {
    /// Resolve the external part of the SAME native force and potential.
    /// Rotating strand offsets use the supplied dynamic frames. Fluid force is
    /// kept separate from inertia; no added mass or flow dynamics is implied.
    ///
    /// # Errors
    /// Invalid poses, moved stationary supports, arithmetic and cancellation refuse.
    pub fn body_loads(
        &self,
        pose: &RodMotionPose,
        gate: &CancelGate,
    ) -> Result<RodBodyLoads, SpatialRodError> {
        let rod = self.rod_at_pose(pose, gate)?;
        let source = self.source();
        let mut force_n = vec![[0.0; 3]; self.inertia.len()];
        if let Some(weight) = source.gravity_nodal_forces_n() {
            force_n.copy_from_slice(weight);
        }
        if let Some(fluid) = source.submerged_fluid_load() {
            for (total, buoyancy) in force_n.iter_mut().zip(fluid.nodal_forces_n()) {
                checkpoint(gate)?;
                for axis in 0..3 {
                    total[axis] += buoyancy[axis];
                }
            }
        }
        let mut torque_body_nm = vec![[0.0; 3]; self.inertia.len()];
        for (index, link) in source.tension_links().iter().enumerate() {
            checkpoint(gate)?;
            let Some(force) = source.link_body_load_force_n(index) else {
                continue;
            };
            for attachment in link.attachments() {
                let force_body = rod.frames[attachment.node]
                    .inverse()
                    .rotate(vector(force))
                    .map_err(crate::rod::group_error)?;
                let torque = vector(attachment.offset_body_m).cross(force_body);
                for (sum, value) in torque_body_nm[attachment.node]
                    .iter_mut()
                    .zip([torque.x, torque.y, torque.z])
                {
                    *sum += value;
                }
            }
        }
        if force_n
            .iter()
            .chain(&torque_body_nm)
            .flatten()
            .any(|v| !v.is_finite())
        {
            return Err(SpatialRodError::Invalid(
                "ambient rod wrench is unrepresentable",
            ));
        }
        let gravity_potential_j = source.gravity_potential_for(&rod)?;
        let buoyancy_potential_j = source.buoyancy_potential_for(&rod)?;
        checkpoint(gate)?;
        Ok(RodBodyLoads {
            force_n,
            torque_body_nm,
            gravity_potential_j,
            buoyancy_potential_j,
        })
    }
}

impl<I: Clone> PreparedRodMotion<I> {
    /// Replace an already bound hydrostatic snapshot through shared material
    /// transport. Reference volumes, mass, rotary inertia, source geometry and
    /// stationary supports are unchanged. This cannot flood a previously sealed
    /// void or add fluid mass to the specimen.
    ///
    /// A continuing time owner must evaluate the returned model at its actual
    /// pose and account for the changed potential before publishing it. This
    /// method alone does not advance time or transfer physical history.
    ///
    /// # Errors
    /// An absent fluid binding, incompatible fluid state or cancellation refuses.
    pub fn rebind_submerged_fluid(
        &self,
        fluid: &impl FluidTransport,
        gate: &CancelGate,
    ) -> Result<Self, SpatialRodError> {
        checkpoint(gate)?;
        let volume_model = self
            .source
            .submerged_fluid_load()
            .ok_or(SpatialRodError::Invalid(
                "fluid refresh requires an existing submerged-volume binding",
            ))?
            .volume_model();
        let source = self
            .source
            .with_submerged_fluid(fluid, volume_model, gate)?;
        let candidate = Self {
            source: Arc::new(source),
            inertia: self.inertia.clone(),
            fixed_position: self.fixed_position.clone(),
            fixed_frame: self.fixed_frame.clone(),
        };
        checkpoint(gate)?;
        Ok(candidate)
    }
}

fn vector(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}
