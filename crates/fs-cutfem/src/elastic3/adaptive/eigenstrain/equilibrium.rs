//! Thermal equilibrium through the original unmasked hanging-node transpose.
use std::ops::ControlFlow;

use fs_solver::op::LinearOp;

use crate::elastic3::ElasticityError3;
use crate::elastic3::adaptive::AdaptiveElasticity3;
use crate::elastic3::eigenstrain::PreparedEigenstrain3;
use crate::elastic3::equilibrium::{EquilibriumResponse3, admit_load, poll, response};
use crate::elastic3::surface::ReferenceLoad3;

impl AdaptiveElasticity3 {
    /// Complete thermal and external RHS in independent master coordinates.
    /// The original full T^T retains supported loads. Thermal columns enter
    /// this vector directly, without allocating a second physical thermal RHS.
    /// Reassemble after scale changes and use this same vector for the algebraic
    /// solve and [`Self::eigenstrain_equilibrium_response`].
    pub fn eigenstrain_equilibrium_load(
        &self,
        prepared: &PreparedEigenstrain3,
        load: ReferenceLoad3<'_>,
        prescribed: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        poll(&mut checkpoint)?;
        self.validate_eigenstrain(prepared)?;
        let mut rhs = self.equilibrium_load(load, prescribed, &mut checkpoint)?;
        self.raw.add_eigenstrain_full_load(
            prepared,
            &mut rhs,
            |node, component, value, out| self.scatter_full_force(node, component, value, out),
            &mut checkpoint,
        )?;
        Ok(rhs)
    }

    /// Full discrete thermal force balance in independent master coordinates.
    ///
    /// Reconstructs each physical cell with the original homogeneous clamps
    /// and applies the unmasked force transpose once. Support force/moment refer
    /// to masters, without duplicating hanging nodes. Total and elastic strain
    /// have separate retained-point admission, including both B-bar projections.
    /// Only the returned residual vector is allocated; there is no full
    /// displacement copy or bulk-response array.
    ///
    /// Supply the matching complete RHS, including supported external and
    /// eigenstrain loads and optional embedded lifting. Strong reactions have
    /// sign `T^T K_full T u - rhs_full`; embedded reactions remain separate.
    #[allow(clippy::too_many_arguments)]
    pub fn eigenstrain_equilibrium_response(
        &self,
        prepared: &PreparedEigenstrain3,
        u: &[f64],
        full_rhs: &[f64],
        maximum_kinematic_strain: f64,
        max_points: usize,
        moment_origin: [f64; 3],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<EquilibriumResponse3, ElasticityError3> {
        admit_load(full_rhs, self.n(), moment_origin, &mut checkpoint)?;
        self.admit_eigenstrain(
            prepared,
            u,
            maximum_kinematic_strain,
            max_points,
            &mut checkpoint,
        )?;
        let force = self.raw.uneliminated_force(
            self.n(),
            |node, component| self.physical_coefficient(u, node, component),
            |node, component, value, out| self.scatter_full_force(node, component, value, out),
            &mut checkpoint,
        )?;
        response(
            force,
            full_rhs,
            &self.nodes,
            &self.fixed,
            moment_origin,
            &mut checkpoint,
        )
    }
}
