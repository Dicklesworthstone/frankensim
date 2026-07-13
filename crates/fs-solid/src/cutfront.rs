//! Legacy fs-solid surface for canonical vector CutFEM elasticity.
//!
//! The assembly, stabilization, deterministic CG solve, topology metadata, and
//! error integration live in fs-cutfem. This module preserves the original
//! fs-solid field-based constructor for existing topology-optimization callers,
//! validates its material through fs-material, translates its historical
//! stabilization convention, and delegates every numerical operation.

use crate::SolidError;
use fs_cutfem::{
    CutElasticity as CanonicalCutElasticity, CutElasticitySolution, CutFemError, CutSdf, Quadtree,
};
use fs_material::IsotropicElastic;

pub use fs_cutfem::{BoundaryTraction, DesignBoxEdge, EdgeBand};

/// Compatibility validity bound attached to the legacy material card.
///
/// The old fs-solid constructor did not expose a strain-domain field. A finite,
/// positive value is required by the canonical frontend; this value is metadata
/// only and does not extend the small-strain validity claim.
const COMPATIBILITY_STRAIN_LIMIT: f64 = 1.0;
const SOLVER_TOL: f64 = 1e-12;
const SOLVER_MAX_ITERS: usize = 60_000;

/// A legacy fs-solid constructor for canonical CutFEM elasticity on
/// Ω = {φ < 0}.
///
/// New code should construct fs_cutfem::CutElasticity directly with an
/// IsotropicElastic material. This facade remains for existing callers that
/// provide Young's modulus and Poisson ratio as separate fields.
pub struct CutElasticity<'a> {
    /// Uniform or 2:1-balanced background quadtree.
    pub grid: &'a Quadtree,
    /// Certified negative-inside level set.
    pub sdf: &'a dyn CutSdf,
    /// Young's modulus.
    pub youngs: f64,
    /// Poisson ratio.
    pub poisson: f64,
    /// Historical Nitsche constant β, applied as β(λ+2μ)/h.
    ///
    /// The facade converts this to fs-cutfem's μ-scaled convention. When
    /// Nitsche interface data is active, the dimensionless translated
    /// coefficient must remain finite.
    pub nitsche_beta: f64,
    /// Historical ghost constant γ, applied as γ(λ+2μ)h per face.
    ///
    /// The facade converts this to fs-cutfem's μ-scaled convention. When
    /// ghost stabilization is active, the dimensionless translated
    /// coefficient must remain finite.
    pub ghost_gamma: f64,
    /// Certified cut-quadrature subdivision depth.
    pub quad_depth: u32,
    /// Optional zero-displacement clamp on design-box boundary nodes.
    pub clamp: Option<&'a dyn Fn(f64, f64) -> bool>,
    /// Optional legacy dead traction with uncertified support on every active
    /// design-box boundary edge. Prefer [`Self::solve_with_boundary_traction`]
    /// for checked named support.
    pub boundary_traction: Option<&'a dyn Fn(f64, f64) -> [f64; 2]>,
    /// Use a natural traction-free embedded interface instead of Nitsche
    /// displacement data.
    pub traction_free_interface: bool,
}

/// Canonical CutFEM solution returned through the legacy fs-solid name.
///
/// This preserves the established nodal, active-cell, and iteration accessors
/// while avoiding a second representation of canonical solver output.
pub type CutSolution = CutElasticitySolution;

impl CutElasticity<'_> {
    fn material(&self) -> Result<IsotropicElastic, SolidError> {
        IsotropicElastic::new(self.youngs, self.poisson, COMPATIBILITY_STRAIN_LIMIT).map_err(
            |error| SolidError::InvalidInput {
                what: format!("CutFEM compatibility material card refused E/nu: {error}"),
            },
        )
    }

    fn canonical<'a>(&'a self, material: &'a IsotropicElastic) -> CanonicalCutElasticity<'a> {
        let (lambda, mu) = material.lame();
        let legacy_stiffness_scale = (lambda + 2.0 * mu) / mu;
        CanonicalCutElasticity {
            grid: self.grid,
            sdf: self.sdf,
            material,
            nitsche_beta: self.nitsche_beta * legacy_stiffness_scale,
            ghost_gamma: self.ghost_gamma * legacy_stiffness_scale,
            quad_depth: self.quad_depth,
            clamp: self.clamp,
            boundary_traction: self.boundary_traction,
            traction_free_interface: self.traction_free_interface,
            solver_tol: SOLVER_TOL,
            solver_max_iters: SOLVER_MAX_ITERS,
        }
    }

    /// Assemble and solve −div σ(u) = f through fs-cutfem's canonical
    /// vector operator.
    ///
    /// # Errors
    /// Returns SolidError::SolveFailed when canonical CG misses its
    /// residual gate. Material, geometry, stabilization, callback, and
    /// unsupported-regime refusals return SolidError::InvalidInput.
    pub fn solve(
        &self,
        f: &dyn Fn(f64, f64) -> [f64; 2],
        g: &dyn Fn(f64, f64) -> [f64; 2],
    ) -> Result<CutSolution, SolidError> {
        let material = self.material()?;
        self.canonical(&material)
            .solve(f, g)
            .map_err(map_cutfem_error)
    }

    /// Assemble and solve through fs-cutfem with an explicit boundary-
    /// traction support descriptor.
    ///
    /// A checked [`BoundaryTraction::EdgeBand`] defines the applied traction
    /// as zero outside one named box-edge interval. This permits unrelated
    /// SDF crossings elsewhere on the same edge while retaining canonical
    /// refusal when the supported segment itself is cut.
    ///
    /// # Errors
    /// Returns [`SolidError::SolveFailed`] when canonical CG misses its
    /// residual gate. Material, geometry, stabilization, callback, support,
    /// and unsupported-regime refusals return [`SolidError::InvalidInput`].
    /// The legacy [`Self::boundary_traction`] field must be `None` so two load
    /// sources are never merged implicitly.
    pub fn solve_with_boundary_traction(
        &self,
        f: &dyn Fn(f64, f64) -> [f64; 2],
        g: &dyn Fn(f64, f64) -> [f64; 2],
        boundary_traction: BoundaryTraction<'_>,
    ) -> Result<CutSolution, SolidError> {
        let material = self.material()?;
        self.canonical(&material)
            .solve_with_boundary_traction(f, g, boundary_traction)
            .map_err(map_cutfem_error)
    }

    /// L2 and H1-seminorm errors through fs-cutfem's canonical integration.
    #[must_use]
    pub fn l2_h1_error(
        &self,
        solution: &CutSolution,
        exact: &dyn Fn(f64, f64) -> [f64; 2],
        exact_gradient: &dyn Fn(f64, f64) -> [[f64; 2]; 2],
    ) -> (f64, f64) {
        // Error integration depends only on grid, SDF, quadrature depth, and
        // the canonical solution. Use a fixed valid card so this infallible
        // legacy accessor does not invent a material-validation panic.
        let integration_material = IsotropicElastic {
            youngs: 1.0,
            poisson: 0.0,
            strain_limit: COMPATIBILITY_STRAIN_LIMIT,
        };
        self.canonical(&integration_material)
            .l2_h1_error(solution, exact, exact_gradient)
    }
}

fn map_cutfem_error(error: CutFemError) -> SolidError {
    match error {
        CutFemError::SolveNotConverged {
            iters,
            rel_residual,
        } => SolidError::SolveFailed {
            iters,
            rel_residual,
        },
        CutFemError::ConstraintCycle { node } => SolidError::InternalInvariant {
            what: format!("canonical fs-cutfem constraint graph cycled at node {node:?}"),
        },
        refusal => SolidError::InvalidInput {
            what: format!("canonical fs-cutfem refused the problem: {refusal}"),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constraint_cycle_is_not_mislabeled_as_caller_input() {
        assert!(matches!(
            map_cutfem_error(CutFemError::ConstraintCycle { node: (3, 5) }),
            SolidError::InternalInvariant { .. }
        ));
    }
}
