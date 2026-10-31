//! Prepare the same native body-force columns in independent master coordinates.

use super::AdaptiveElasticity3;
use crate::elastic3::ElasticityError3;
use crate::elastic3::body::{BodyLoadSpace3, MaterialBodyLoad3};
use std::ops::ControlFlow;

impl AdaptiveElasticity3 {
    /// Integrate `rho0 N^T acceleration` once and apply the original hanging-node
    /// transpose to each cell column. Later occupancy actions/pullbacks operate
    /// directly on the independent master space without expanding raw fields.
    ///
    /// The scalar-term allowance includes the current cell's constrained
    /// assembly slots. Strong homogeneous master clamps carry no displacement
    /// work; their load entries are omitted exactly as in `body_load`.
    pub fn prepare_material_body_load(
        &self,
        acceleration: &dyn Fn([f64; 3]) -> [f64; 3],
        max_terms: usize,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<MaterialBodyLoad3, ElasticityError3> {
        MaterialBodyLoad3::prepare(
            &self.raw,
            Some(&self.rows),
            &self.fixed,
            BodyLoadSpace3::Adaptive,
            acceleration,
            max_terms,
            checkpoint,
        )
    }

    /// Refuse a stencil from another geometry, material field or master map.
    /// Wrapping this same operator in a preconditioner preserves the binding.
    pub fn validate_material_body_load(
        &self,
        load: &MaterialBodyLoad3,
    ) -> Result<(), ElasticityError3> {
        load.validate_for(&self.raw, BodyLoadSpace3::Adaptive, 3 * self.nodes.len())
    }
}
