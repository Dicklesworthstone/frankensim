//! Complete thermal forces and strong-support balance on the original operator.
use std::ops::ControlFlow;

use fs_solver::op::LinearOp;

use super::{PreparedEigenstrain3, poll};
use crate::elastic3::equilibrium::{EquilibriumResponse3, admit_load, response};
use crate::elastic3::surface::ReferenceLoad3;
use crate::elastic3::{CutElasticity3, ElasticityError3};

impl CutElasticity3 {
    /// Complete external, prescribed-motion and free-strain right-hand side.
    ///
    /// Retains ALL strongly supported entries of the original eigenstrain
    /// columns, including their embedded Nitsche contribution. The existing
    /// load owner supplies external body/surface forces and optional lifting.
    /// Neither external dead loads nor source density are stiffness-scaled.
    ///
    /// Use this same vector for `solve_controlled` and
    /// [`Self::eigenstrain_equilibrium_response`]. The solver alone masks the
    /// strong supports. Add independent nodal loads to the complete vector,
    /// including loads applied at supports, and reassemble after scale changes.
    /// Masked `eigenstrain_load` vectors cannot recover discarded support forces.
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
        self.add_eigenstrain_full_load(
            prepared,
            &mut rhs,
            |node, component, value, out| out[3 * node + component] += value,
            &mut checkpoint,
        )?;
        Ok(rhs)
    }

    // Accumulate into the caller's owned complete vector. The physical force
    // scatter is unmasked; adaptive callers use their original full T^T here.
    // No second global thermal vector or source/quadrature reconstruction.
    pub(in crate::elastic3) fn add_eigenstrain_full_load(
        &self,
        prepared: &PreparedEigenstrain3,
        rhs: &mut [f64],
        scatter: impl Fn(usize, usize, f64, &mut [f64]),
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        poll(checkpoint)?;
        self.validate_eigenstrain(prepared)?;
        if prepared.cells.len() != self.cells.len() {
            return Err(ElasticityError3::Invalid(
                "complete eigenstrain load cell alignment mismatch",
            ));
        }
        for (id, (cell, free)) in self.cells.iter().zip(&prepared.cells).enumerate() {
            poll(checkpoint)?;
            for (i, &force) in free.force.iter().enumerate() {
                scatter(cell.nodes[i / 3], i % 3, self.scales[id] * force, rhs);
            }
        }
        admit_load(rhs, rhs.len(), [0.0; 3], checkpoint)?;
        poll(checkpoint)?;
        Ok(())
    }

    /// Recover `K_full u - rhs_full` with thermal elastic-domain admission.
    ///
    /// The complete RHS must match this preparation, current scales, external
    /// loads and optional motion; obtain it from
    /// [`Self::eigenstrain_equilibrium_load`]. Fixed residual entries are forces
    /// exerted ON the solid by strong supports. Free entries expose imbalance.
    /// Moments use `(node - moment_origin) cross reaction`.
    ///
    /// Total small-deformation kinematics and actual/projected elastic strain
    /// are admitted separately at retained bulk points, exactly as in
    /// [`Self::eigenstrain_response`]. An arbitrary eigenstrain callback cannot
    /// certify between-point strain extrema. No bulk response array is allocated.
    /// Original B-bar, ghost, embedded-support and strong-clamp actions remain.
    /// Embedded reactions are reported separately by
    /// [`Self::eigenstrain_embedded_reaction`], without double counting them.
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
        let force = self.uneliminated_force(
            self.n(),
            |node, component| {
                if self.fixed[node] {
                    0.0
                } else {
                    u[3 * node + component]
                }
            },
            |node, component, value, out| out[3 * node + component] += value,
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
