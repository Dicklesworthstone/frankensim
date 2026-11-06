//! Completed brittle cuts preserve the body's native physical motion.
//!
//! The original solid owner admits only failed, aligned, traction-free faces.
//! Its histories and geometric opening guards remain attached to the split
//! solid. This module duplicates the actual velocity field through the native
//! topology owner and reconstructs consistent inertia on the new nodal stars.
//! No impulse, elapsed time, artificial velocity or mass deletion is involved.

/// Simultaneous native normal impacts on the same released material faces.
pub mod contact;

use std::sync::Arc;

use fs_exec::Cx;

use super::{BrittleTetDynamicError, BrittleTetMotion};
use crate::brittle::solid::fracture::{
    BrittleTetFractureControls, BrittleTetFractureEnergy, BrittleTetFractureError,
};
use crate::brittle::solid::{BrittleTetError, checkpoint, energy_gate};
use crate::brittle::topology::{P1PhysicalInventory, TetFractureTopology};
use crate::inertia::P1ReferenceMass;
use crate::loads3::TetLoadGeometry3;
use crate::plastic3::SmallStrainTetError;

/// Accepted topology change with the original clock and constitutive history.
#[derive(Debug)]
pub struct BrittleTetFractureMotionStep<'a> {
    pub motion: BrittleTetMotion<'a>,
    /// Original bulk storage and cumulative fracture loss, without new work.
    pub energy: BrittleTetFractureEnergy,
    /// Independent consistent-P1 integrals on the original and split meshes.
    pub inventory_before: P1PhysicalInventory,
    pub inventory_after: P1PhysicalInventory,
    /// Ordered by the topology owner's `element_fragments()` identifiers.
    pub fragments: Vec<P1PhysicalInventory>,
    /// Retained old and reconstructed new native consistent-mass storage.
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Arithmetic difference in native kinetic storage; not physical loss.
    pub kinetic_defect_j: f64,
}

/// Physical release, topology, inertia or retained-energy admission failure.
#[derive(Clone, Debug, PartialEq)]
pub enum BrittleTetFractureMotionError {
    Fracture(BrittleTetFractureError),
    Motion(BrittleTetDynamicError),
}

impl core::fmt::Display for BrittleTetFractureMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "brittle fragment motion transfer refused: {self:?}")
    }
}

impl std::error::Error for BrittleTetFractureMotionError {}

impl From<BrittleTetFractureError> for BrittleTetFractureMotionError {
    fn from(error: BrittleTetFractureError) -> Self {
        Self::Fracture(error)
    }
}

impl From<BrittleTetDynamicError> for BrittleTetFractureMotionError {
    fn from(error: BrittleTetDynamicError) -> Self {
        Self::Motion(error)
    }
}

impl From<BrittleTetError> for BrittleTetFractureMotionError {
    fn from(error: BrittleTetError) -> Self {
        Self::Motion(error.into())
    }
}

impl<'a> BrittleTetMotion<'a> {
    /// Split fully failed, traction-free internal faces without a physical step.
    ///
    /// `target_geometry` must describe exactly the original topology owner's
    /// split mesh. Source bindings, element order and material history are
    /// retained by the original solid handoff. This method additionally retains
    /// the accepted reference-density snapshot even when later source queries
    /// would return another density. Mass is repartitioned among the new stars
    /// using the native consistent operator, never row-sum lumping.
    ///
    /// The topology owner independently reconstructs mass, first moment, linear
    /// and angular momentum, and kinetic energy on both meshes. Native kinetic
    /// storage must agree with those integrals within the explicit absolute
    /// `kinetic_tolerance_j`. No impulse or kinetic loss is attributed to a cut.
    /// Physical time and used/remaining step counts remain unchanged, including
    /// when the motion has exhausted its physical-step allowance.
    ///
    /// Subsequent native motion retains the solid's paired-face opening guard.
    /// Closure beyond its opening tolerance refuses: this entry does not add
    /// contact forces, arbitrary fragment collision detection or finite-rotation
    /// kinematics.
    /// Any refusal leaves the entire pre-cut body available for retry.
    pub fn split_fracture<'next>(
        &self,
        topology: &TetFractureTopology<'_>,
        target_geometry: &'next TetLoadGeometry3<'next>,
        controls: BrittleTetFractureControls,
        kinetic_tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetFractureMotionStep<'next>, BrittleTetFractureMotionError>
    where
        'a: 'next,
    {
        checkpoint(cx)?;
        if !kinetic_tolerance_j.is_finite() || kinetic_tolerance_j < 0.0 {
            return Err(BrittleTetDynamicError::Invalid(
                "finite nonnegative fracture-transfer kinetic tolerance required",
            )
            .into());
        }
        let fracture = self
            .solid
            .split_fracture(topology, target_geometry, controls, cx)?;
        let transfer = topology
            .transfer_motion(
                self.solid.solution.displacement_m.as_flattened(),
                self.velocity_mps.as_flattened(),
                &self.reference_density_kg_m3,
                cx,
            )
            .map_err(BrittleTetFractureError::Topology)?;
        if transfer.displacement_m != fracture.solid.solution.displacement_m.as_flattened() {
            return Err(BrittleTetDynamicError::Invalid(
                "solid and motion handoffs must retain the same displacement field",
            )
            .into());
        }
        let mass = P1ReferenceMass::new(
            target_geometry
                .tetrahedra()
                .iter()
                .zip(target_geometry.element_volumes_m3())
                .map(|(&nodes, &volume)| (nodes, volume)),
            &self.reference_density_kg_m3,
            transfer.velocity_m_s.len(),
            cx,
        )
        .map_err(SmallStrainTetError::from)
        .map_err(BrittleTetDynamicError::from)?;
        // Element order, volume and density are unchanged, so both native
        // mass constructors perform the identical scalar inventory reduction.
        if mass.reference_mass_kg() != self.mass.reference_mass_kg() {
            return Err(BrittleTetDynamicError::Invalid(
                "a fracture cut cannot change the native reference mass inventory",
            )
            .into());
        }
        let momentum = mass
            .action(&transfer.velocity_m_s, cx)
            .map_err(SmallStrainTetError::from)
            .map_err(BrittleTetDynamicError::from)?;
        let mut kinetic_after_j = 0.0;
        for (dof, (&velocity, &momentum)) in transfer.velocity_m_s.iter().zip(&momentum).enumerate()
        {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            kinetic_after_j = velocity.mul_add(momentum, kinetic_after_j);
        }
        kinetic_after_j *= 0.5;
        if !kinetic_after_j.is_finite() || kinetic_after_j < 0.0 {
            return Err(BrittleTetDynamicError::Invalid(
                "finite nonnegative native kinetic storage required after a cut",
            )
            .into());
        }
        let kinetic_defect_j = kinetic_after_j - self.kinetic_energy_j;
        for (account, defect) in [
            (
                "pre-fracture consistent kinetic inventory",
                transfer.before.kinetic_energy_j - self.kinetic_energy_j,
            ),
            (
                "post-fracture consistent kinetic inventory",
                kinetic_after_j - transfer.after.kinetic_energy_j,
            ),
            ("fracture kinetic-storage continuity", kinetic_defect_j),
        ] {
            energy_gate(account, defect, kinetic_tolerance_j)?;
        }
        let mut velocity_mps = Vec::with_capacity(target_geometry.nodes_m().len());
        for (index, velocity) in transfer.velocity_m_s.chunks_exact(3).enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            velocity_mps.push([velocity[0], velocity[1], velocity[2]]);
        }
        checkpoint(cx)?;
        Ok(BrittleTetFractureMotionStep {
            motion: BrittleTetMotion {
                solid: fracture.solid,
                velocity_mps,
                reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                mass: Arc::new(mass),
                kinetic_energy_j: kinetic_after_j,
                time_s: self.time_s,
                maximum_steps: self.maximum_steps,
                accepted_steps: self.accepted_steps,
            },
            energy: fracture.energy,
            inventory_before: transfer.before,
            inventory_after: transfer.after,
            fragments: transfer.fragments,
            kinetic_before_j: self.kinetic_energy_j,
            kinetic_after_j,
            kinetic_defect_j,
        })
    }
}
