//! Completed thermally driven tensile fractures retain native physical motion.
//!
//! Original thermal sources, reference expansion, mechanical history and total
//! displacement remain with the thermal solid. The shared completed-fracture
//! transfer reconstructs its original consistent mass on the split nodal stars,
//! preserving velocity, momentum, kinetic storage and the physical clock.
//! A connectivity release supplies no impulse, heat or fresh fracture energy.

/// Native simultaneous normal impacts on the retained released thermal faces.
pub mod contact;

use std::sync::Arc;

use fs_exec::Cx;

use super::ThermalBrittleTetMotion;
use crate::brittle::solid::checkpoint;
use crate::brittle::solid::dynamics::BrittleTetDynamicError;
use crate::brittle::solid::dynamics::fracture::{
    BrittleTetFractureMotionError, FractureMotionView, transfer_fracture_motion,
};
use crate::brittle::solid::fracture::{BrittleTetFractureControls, BrittleTetFractureEnergy};
use crate::brittle::topology::{P1PhysicalInventory, TetFractureTopology};
use crate::loads3::TetLoadGeometry3;

/// Accepted thermal fragment motion with unchanged material history and time.
#[derive(Debug)]
pub struct ThermalBrittleTetFractureMotionStep<'a> {
    pub motion: ThermalBrittleTetMotion<'a>,
    /// Original bulk storage and cumulative fracture loss; no new surface work.
    pub energy: BrittleTetFractureEnergy,
    pub inventory_before: P1PhysicalInventory,
    pub inventory_after: P1PhysicalInventory,
    /// Ordered by the native topology owner's element-fragment identifiers.
    pub fragments: Vec<P1PhysicalInventory>,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Arithmetic storage difference, never physical dissipation.
    pub kinetic_defect_j: f64,
}

impl<'a> ThermalBrittleTetMotion<'a> {
    /// Release completed, aligned, traction-free faces while retaining motion.
    ///
    /// The thermal solid first admits the target as the exact native split mesh.
    /// Both incident native histories must have completely failed on the face;
    /// source references, total strain, mechanical strain, physical stress and
    /// accumulated fracture loss remain unchanged.
    ///
    /// The original reference-density snapshot is repartitioned through the same
    /// consistent P1 mass transfer used by nonthermal fracture. A successor
    /// temperature's density cannot add matter or change momentum. Independent
    /// consistent-P1 inventories and native kinetic storage must agree within
    /// the explicit absolute kinetic tolerance. No physical time or step budget
    /// is consumed, even when the existing motion has no physical steps left.
    ///
    /// Subsequent thermal backward-Euler and midpoint endpoints retain the
    /// released matching-face opening guards. Closing endpoints require contact
    /// and refuse; this cut does not supply finite sliding or collision search.
    /// Every refusal leaves the complete original thermal motion available.
    pub fn split_fracture<'next>(
        &self,
        topology: &TetFractureTopology<'_>,
        target_geometry: &'next TetLoadGeometry3<'next>,
        controls: BrittleTetFractureControls,
        kinetic_tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetFractureMotionStep<'next>, BrittleTetFractureMotionError>
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
        let transfer = transfer_fracture_motion(
            FractureMotionView {
                displacement_m: &self.solid.solution.displacement_m,
                velocity_mps: &self.velocity_mps,
                reference_density_kg_m3: &self.reference_density_kg_m3,
                mass: &self.mass,
                kinetic_energy_j: self.kinetic_energy_j,
            },
            topology,
            target_geometry,
            &fracture.solid.solution.displacement_m,
            kinetic_tolerance_j,
            cx,
        )?;
        checkpoint(cx)?;
        Ok(ThermalBrittleTetFractureMotionStep {
            motion: ThermalBrittleTetMotion {
                solid: fracture.solid,
                velocity_mps: transfer.velocity_mps,
                reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                mass: Arc::new(transfer.mass),
                kinetic_energy_j: transfer.kinetic_after_j,
                time_s: self.time_s,
                maximum_steps: self.maximum_steps,
                accepted_steps: self.accepted_steps,
            },
            energy: fracture.energy,
            inventory_before: transfer.inventory_before,
            inventory_after: transfer.inventory_after,
            fragments: transfer.fragments,
            kinetic_before_j: self.kinetic_energy_j,
            kinetic_after_j: transfer.kinetic_after_j,
            kinetic_defect_j: transfer.kinetic_defect_j,
        })
    }
}
