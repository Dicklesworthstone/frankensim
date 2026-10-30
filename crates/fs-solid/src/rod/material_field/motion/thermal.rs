//! Continuing temperature fields on the original material-coordinate rod.
//!
//! `RodMaterialSource` remains the only owner of source selection, temperature
//! queries and integrated expansion. A refresh changes section coefficients and
//! thermal strains, while the original reference mass, geometry, attached links,
//! gravity, displaced volumes and stationary supports remain intact.

use std::sync::Arc;

use fs_exec::CancelGate;

use super::{PreparedRodMotion, RodMotionPose, SpatialRodError, admit_section, checkpoint};
use crate::rod::material_field::SpatialMaterialRod;
use crate::rod::material_field::thermal::source::RodMaterialSource;

/// Actual-pose native strain and section energy, excluding separate links/loads.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RodSegmentStrainState {
    /// Axial, shear, shear; intrinsic geometry removed, thermal strain retained.
    pub gamma: [f64; 3],
    /// Twist and two bending curvatures [rad/m], relative to reference geometry.
    pub curvature: [f64; 3],
    pub strain_energy_j: f64,
}

impl<I> PreparedRodMotion<I> {
    /// Observe each original section at one admitted pose through the SAME
    /// kinematic and energy owners used by analytical forces and motion.
    pub fn segment_strain_states(
        &self,
        pose: &RodMotionPose,
        gate: &CancelGate,
    ) -> Result<Vec<RodSegmentStrainState>, SpatialRodError> {
        let rod = self.rod_at_pose(pose, gate)?;
        let mut states = Vec::with_capacity(self.source.sections.len());
        for segment in 0..self.source.sections.len() {
            checkpoint(gate)?;
            let (gamma, curvature) = self.source.strains_for(&rod, segment)?;
            states.push(RodSegmentStrainState {
                gamma,
                curvature,
                strain_energy_j: self.source.segment_energy_for(&rod, segment)?,
            });
        }
        Ok(states)
    }
    /// Current section temperatures in segment/layer order [K].
    /// A mechanical source without an explicit `T` coordinate refuses; this
    /// operation does not invent a reference or ambient body temperature.
    pub fn section_temperatures_k(&self) -> Result<Vec<f64>, SpatialRodError> {
        temperatures(self.source())
    }

    /// Native conservative strain storage at the actual admitted dynamic pose.
    /// Includes attached-link storage and excludes gravity/hydrostatic potential.
    /// This is the original energy owner, not a new thermal energy definition.
    pub fn strain_energy_j(
        &self,
        pose: &RodMotionPose,
        gate: &CancelGate,
    ) -> Result<f64, SpatialRodError> {
        let rod = self.rod_at_pose(pose, gate)?;
        let energy = self.source.energy_for(&rod)?;
        checkpoint(gate)?;
        Ok(energy)
    }
}

impl<I: Clone> PreparedRodMotion<I> {
    /// Resolve the SAME source description at new local section temperatures.
    ///
    /// The complete original reference inventory and reference geometry must
    /// agree. Resolving the old temperatures must reproduce every current
    /// section, thermal path and strain bound exactly, so a caller cannot swap
    /// selected claims, stress-free temperatures, shear conventions or source
    /// families through this temperature-only operation. An initially unheated
    /// reference rod may acquire its original zero-strain thermal binding.
    ///
    /// The existing material source re-evaluates E/rho/G and integrates alpha
    /// from its ORIGINAL reference. Current rho(T) never changes closed-body
    /// mass or inertia. Geometry and reference displaced volumes stay frozen;
    /// attached tension links retain their separately bound material state.
    ///
    /// This model-level operation does not advance motion or supply heat. A
    /// continuing time owner must evaluate it at the actual pose and account
    /// for the changed strain storage before publishing the candidate.
    pub fn rebind_material_temperatures(
        &self,
        source: &RodMaterialSource<'_>,
        temperatures_k: &[f64],
        gate: &CancelGate,
    ) -> Result<Self, SpatialRodError> {
        checkpoint(gate)?;
        let reference = source.reference_rod();
        let same_geometry = match (
            self.source.reference_geometry(),
            reference.reference_geometry(),
        ) {
            (None, None) => true,
            (Some(a), Some(b)) => {
                a.length_m().to_bits() == b.length_m().to_bits()
                    && a.positions_m() == b.positions_m()
                    && a.frames() == b.frames()
            }
            _ => false,
        };
        if self.source.reference_mass() != reference.reference_mass()
            || self.source.reference_length_m().to_bits()
                != reference.reference_length_m().to_bits()
            || !same_geometry
        {
            return Err(SpatialRodError::Invalid(
                "thermal continuation must retain the exact original rod source inventory and reference geometry",
            ));
        }
        let old_temperatures = self.section_temperatures_k()?;
        let expected = source.resolve_at(&old_temperatures, gate)?;
        let same_response = match &self.source.thermal {
            Some(_) => {
                self.source.sections == expected.sections && self.source.thermal == expected.thermal
            }
            // The admitted reference carries its original selection-policy
            // receipts; future resolutions carry the pinned selection receipts.
            // Compare the exact original reference, never just sampled values.
            None => {
                old_temperatures == source.reference_temperatures_k()
                    && self.source.sections == reference.sections
            }
        };
        if !same_response {
            return Err(SpatialRodError::Invalid(
                "thermal continuation source must reproduce the current section claims, expansion reference and strain bounds",
            ));
        }
        let resolved = source.resolve_at(temperatures_k, gate)?;
        for section in &resolved.sections {
            checkpoint(gate)?;
            admit_section(section)?;
        }
        // Preserve the original retained support pose and every separate body
        // load/link owner. Only source-backed section response is refreshed.
        let mut updated = (*self.source).clone();
        updated.sections = resolved.sections;
        updated.stiffness = resolved.stiffness;
        updated.thermal = resolved.thermal;
        updated.rod.section = resolved.rod.section;
        let candidate = Self {
            source: Arc::new(updated),
            inertia: self.inertia.clone(),
            fixed_position: self.fixed_position.clone(),
            fixed_frame: self.fixed_frame.clone(),
        };
        checkpoint(gate)?;
        Ok(candidate)
    }
}

fn temperatures(source: &SpatialMaterialRod) -> Result<Vec<f64>, SpatialRodError> {
    source
        .sections()
        .iter()
        .flat_map(|section| section.elastic().layers())
        .map(|layer| {
            layer
                .material
                .query_point()
                .iter()
                .find_map(|(axis, temperature)| (axis == "T").then_some(*temperature))
                .filter(|temperature| temperature.is_finite() && *temperature > 0.0)
                .ok_or(SpatialRodError::Invalid(
                    "thermal continuation requires an explicit positive section temperature",
                ))
        })
        .collect()
}
