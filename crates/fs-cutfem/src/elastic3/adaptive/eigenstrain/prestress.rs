//! Original Q1 reconstruction and transpose for frozen adaptive prestress.
use std::ops::ControlFlow;

use fs_solver::op::LinearOp;

use crate::elastic3::adaptive::AdaptiveElasticity3;
use crate::elastic3::eigenstrain::PreparedEigenstrain3;
use crate::elastic3::{ElasticityError3, PreparedPrestress3};

impl AdaptiveElasticity3 {
    /// Prepare the original signed geometric operator as T^T G T.
    ///
    /// Physical point stress uses the same homogeneous master reconstruction
    /// as elasticity, including B-bar projection of total and free strain.
    /// Selected Nitsche geometric flux terms are condensed with the same T.
    /// Original elastic penalties and ghost coefficients remain frozen in K_e.
    ///
    /// `max_points` bounds bulk plus selected surface stress evaluations.
    /// `max_entries` bounds upper-triangular scalar pairs, including the current
    /// cell's constrained map, before replication over xyz components.
    /// This state is physically admitted at retained points, not certified in
    /// equilibrium or throughout the continuum. No new geometry is integrated.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_eigenstrain_prestress(
        &self,
        prepared: &PreparedEigenstrain3,
        u: &[f64],
        maximum_kinematic_strain: f64,
        max_points: usize,
        max_entries: usize,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<PreparedPrestress3, ElasticityError3> {
        PreparedPrestress3::prepare(
            &self.raw,
            Some(&self.rows),
            &self.fixed,
            prepared,
            u,
            maximum_kinematic_strain,
            max_points,
            max_entries,
            |cell, field| {
                std::array::from_fn(|i| self.physical_coefficient(field, cell.nodes[i / 3], i % 3))
            },
            checkpoint,
        )
    }

    /// Require the original master space/material owner and exact frozen scales.
    /// A changed preload state needs new preparation even when scales match.
    pub fn validate_prestress(
        &self,
        prestress: &PreparedPrestress3,
    ) -> Result<(), ElasticityError3> {
        prestress.validate_for(&self.raw, true, self.n())
    }
}
