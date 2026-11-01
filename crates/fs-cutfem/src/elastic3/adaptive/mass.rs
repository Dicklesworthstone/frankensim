//! Consistent native-density inertia in independent adaptive master coordinates.

use std::ops::ControlFlow;

use super::AdaptiveElasticity3;
use crate::elastic3::{ElasticityError3, mass::MaterialMass3};

impl AdaptiveElasticity3 {
    /// Prepare `T^T M_cell T` once using the exact hanging-node reconstruction
    /// and homogeneous master clamps. No raw-field expansion is needed during
    /// later actions, diagonals or cell-occupancy pullbacks.
    ///
    /// `max_terms` counts upper-triangular scalar node-pair slots, including
    /// the current cell's constrained assembly map, before replication across
    /// the three vector components. Ghost penalties contribute no inertia.
    pub fn prepare_material_mass(
        &self,
        max_terms: usize,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<MaterialMass3, ElasticityError3> {
        MaterialMass3::prepare(
            &self.raw,
            Some(&self.rows),
            &self.fixed,
            max_terms,
            checkpoint,
        )
    }

    /// Require the original geometry, native material field and master space.
    pub fn validate_material_mass(&self, mass: &MaterialMass3) -> Result<(), ElasticityError3> {
        mass.validate_for(&self.raw, true, 3 * self.nodes.len())
    }
}
