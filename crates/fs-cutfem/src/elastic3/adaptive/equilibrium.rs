//! Strong-support reactions in the original independent master coordinates.
use super::*;
use crate::elastic3::equilibrium::{EquilibriumResponse3, admit_load, poll, response};
use crate::elastic3::surface::ReferenceLoad3;

impl AdaptiveElasticity3 {
    /// Complete body/surface/motion RHS, including strong master supports.
    /// The original hanging-node transpose retains supported entries here;
    /// the solver still masks them when enforcing homogeneous displacement.
    /// Use this SAME vector for solving and `equilibrium_response`. Existing
    /// masked/prepared loads cannot recover forces discarded at their supports.
    pub fn equilibrium_load(
        &self,
        load: ReferenceLoad3<'_>,
        prescribed: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        let raw = self
            .raw
            .equilibrium_load(load, prescribed, &mut checkpoint)?;
        let mut rhs = vec![0.0; self.n()];
        for node in 0..self.rows.len() {
            poll(&mut checkpoint)?;
            for component in 0..3 {
                self.scatter_full_force(node, component, raw[3 * node + component], &mut rhs);
            }
        }
        admit_load(&rhs, self.n(), [0.0; 3], &mut checkpoint)?;
        poll(&mut checkpoint)?;
        Ok(rhs)
    }

    /// `T^T K_full T u - rhs_full` with strong support rows retained.
    /// Input supports and hanging nodes use exactly the operator's reconstruction;
    /// force transfer uses its unmasked transpose. Residual entries, support forces
    /// and moment arms refer to independent master nodes, not duplicated slaves.
    ///
    /// `full_rhs` must include supported loads and any embedded motion lifting.
    /// Original material strain limits apply before force evaluation. The reported
    /// free residual is an equilibrium diagnostic, not a new convergence claim.
    /// Only the output vector is allocated; no full displacement copy is needed.
    pub fn equilibrium_response(
        &self,
        u: &[f64],
        full_rhs: &[f64],
        moment_origin: [f64; 3],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<EquilibriumResponse3, ElasticityError3> {
        admit_load(full_rhs, self.n(), moment_origin, &mut checkpoint)?;
        self.admit_material_domain(u, &mut checkpoint)?;
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

    pub(super) fn scatter_full_force(
        &self,
        node: usize,
        component: usize,
        value: f64,
        out: &mut [f64],
    ) {
        for &(master, weight) in &self.rows[node] {
            out[3 * master + component] += weight * value;
        }
    }
}
