//! Source-bound proportional buckling for existing three-dimensional solids.
//!
//! Reuse the accepted reference stress and its original homogeneous incremental
//! supports. Solve `(K_e + lambda K_g) phi = 0`, not `(K_t + lambda K_g) phi = 0`:
//! the latter would count the accepted prestress twice and shift every load
//! factor by one. No material law, stress integration, or eigensolver is copied.
//!
//! A multiplier scales the COMPLETE reference stress with geometry, elastic
//! tangent and supports frozen. It is a linear buckling screen, not a nonlinear
//! load path. In particular it does not predict temperature-dependent stiffness,
//! changing supports, follower loads, or the material validity of the scaled
//! state. Fixed and proportional preloads must not be mixed in one reference.

use super::{BucklingResult, buckling_loads_cancellable};
use crate::SolidError;
use crate::linear3::{QuadraticPrestressedAssembly, TetElasticAssembly};
use crate::prestress::PrestressedTetAssembly;
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_sparse::Csr;

/// Critical proportional load factors tied to one accepted 3-D reference state.
/// The underlying free-DOF modes retain their geometric-work normalization and
/// independently recomputed residuals. Supports have zero incremental motion.
#[derive(Debug, Clone)]
pub struct SolidBucklingResult3 {
    reference_assembly_identity: ContentHash,
    node_count: usize,
    result: BucklingResult,
}

impl SolidBucklingResult3 {
    /// Exact accepted prestress assembly used to form the physical pencil.
    #[must_use]
    pub const fn reference_assembly_identity(&self) -> ContentHash {
        self.reference_assembly_identity
    }

    /// Load factors, free-DOF modes, residuals, mapping and iteration count.
    #[must_use]
    pub const fn result(&self) -> &BucklingResult {
        &self.result
    }

    /// Expand one mode to original node-major xyz components, retaining every
    /// supported node. No planar grouping or inference from the last free DOF.
    pub fn nodal_mode(&self, mode_index: usize) -> Result<Vec<[f64; 3]>, SolidError> {
        let mode = self
            .result
            .modes
            .get(mode_index)
            .ok_or_else(|| invalid("buckling mode index is outside the accepted result"))?;
        let mut full = vec![[0.0; 3]; self.node_count];
        for (dof, free) in self.result.dof_map.iter().enumerate() {
            if let Some(free) = free {
                full[dof / 3][dof % 3] = mode[*free];
            }
        }
        Ok(full)
    }
}

impl PrestressedTetAssembly {
    /// Buckling multipliers of the complete accepted P1 reference stress.
    ///
    /// Uses the original elastic stiffness K_e and signed geometric stiffness
    /// K_g. The already prestressed tangent K_e + K_g is deliberately not used
    /// as the unscaled matrix. Reference mass and geometry remain unchanged.
    /// Source-bound thermal, nonthermal eigenstrain, dead-load and prescribed-
    /// support assemblies all enter through the same existing stress owner.
    pub fn proportional_buckling(
        &self,
        modes: usize,
        maximum_iterations: usize,
        cx: &Cx<'_>,
    ) -> Result<SolidBucklingResult3, SolidError> {
        solve(
            self.pencil(),
            self.elastic_stiffness(),
            self.geometric_stiffness(),
            self.stress().nodal_internal_force_n().len(),
            self.identity(),
            modes,
            maximum_iterations,
            cx,
        )
    }
}

impl QuadraticPrestressedAssembly {
    /// The same proportional buckling screen for TET10, HEX8, HEX20 and HEX27.
    /// The accepted compatible or B-bar stress has already been integrated at
    /// the original stiffness points; no element-average stress is substituted.
    pub fn proportional_buckling(
        &self,
        modes: usize,
        maximum_iterations: usize,
        cx: &Cx<'_>,
    ) -> Result<SolidBucklingResult3, SolidError> {
        solve(
            self.pencil(),
            self.elastic_stiffness(),
            self.geometric_stiffness(),
            self.response().displacement_m().len(),
            self.identity(),
            modes,
            maximum_iterations,
            cx,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn solve(
    pencil: &TetElasticAssembly,
    elastic: &Csr,
    geometric: &Csr,
    node_count: usize,
    reference_assembly_identity: ContentHash,
    modes: usize,
    maximum_iterations: usize,
    cx: &Cx<'_>,
) -> Result<SolidBucklingResult3, SolidError> {
    cx.checkpoint().map_err(|_| SolidError::Cancelled)?;
    let full_dofs = node_count
        .checked_mul(3)
        .ok_or_else(|| invalid("three-dimensional buckling DOF count overflows"))?;
    if pencil.free_dofs.len() != elastic.nrows() {
        return Err(invalid(
            "three-dimensional buckling map differs from its elastic pencil",
        ));
    }
    let mut dof_map = vec![None; full_dofs];
    for (free, &full) in pencil.free_dofs.iter().enumerate() {
        if free % 256 == 0 {
            cx.checkpoint().map_err(|_| SolidError::Cancelled)?;
        }
        if full >= full_dofs || dof_map[full].is_some() {
            return Err(invalid("three-dimensional buckling map is not one-to-one"));
        }
        dof_map[full] = Some(free);
    }
    let result = buckling_loads_cancellable(
        elastic,
        geometric,
        &dof_map,
        modes,
        maximum_iterations,
        &|| cx.checkpoint().is_err(),
    )?;
    cx.checkpoint().map_err(|_| SolidError::Cancelled)?;
    Ok(SolidBucklingResult3 {
        reference_assembly_identity,
        node_count,
        result,
    })
}

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}
