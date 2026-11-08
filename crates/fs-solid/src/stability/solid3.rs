//! Source-bound proportional and fixed-preload buckling for 3-D solids.
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
//! state. For a fixed preload, pair compatible accepted base/reference states
//! through `incremental_buckling`: only their stress difference is scaled.

use super::{BucklingResult, buckling_loads_cancellable};
use crate::SolidError;
use crate::linear3::{QuadraticPrestressedAssembly, TetElasticAssembly};
use crate::prestress::PrestressedTetAssembly;
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_sparse::Csr;

/// Critical load factors tied to one accepted reference, or a base/reference pair.
/// The underlying free-DOF modes retain their geometric-work normalization and
/// independently recomputed residuals. Supports have zero incremental motion.
#[derive(Debug, Clone)]
pub struct SolidBucklingResult3 {
    base_assembly_identity: Option<ContentHash>,
    reference_assembly_identity: ContentHash,
    node_count: usize,
    result: BucklingResult,
}

impl SolidBucklingResult3 {
    /// Fixed base equilibrium for incremental loading; `None` for proportional
    /// loading of the complete reference stress. In incremental loading,
    /// lambda=0 denotes this base and lambda=1 denotes the reference state.
    #[must_use]
    pub const fn base_assembly_identity(&self) -> Option<ContentHash> {
        self.base_assembly_identity
    }

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

    /// Hold this accepted preload fixed and scale the load increment that
    /// reaches `reference` at lambda=1. Solve
    /// `(K_e + K_g,base + lambda (K_g,reference - K_g,base)) phi = 0`.
    ///
    /// Both states must retain the exact original geometry, material sources
    /// and coefficients, mesh frame and incremental support space. Thermal
    /// strains, dead forces and prescribed support values may differ. This
    /// models their affine load path with frozen stiffness, not a material or
    /// geometric update. The base tangent must be positive definite; the shared
    /// buckling owner returns `SolveFailed` for an unstable/singular base.
    pub fn incremental_buckling(
        &self,
        reference: &Self,
        modes: usize,
        maximum_iterations: usize,
        cx: &Cx<'_>,
    ) -> Result<SolidBucklingResult3, SolidError> {
        solve_increment(
            self.buckling_reference(),
            reference.buckling_reference(),
            modes,
            maximum_iterations,
            cx,
        )
    }

    fn buckling_reference(&self) -> BucklingReference<'_> {
        BucklingReference {
            pencil: self.pencil(),
            elastic: self.elastic_stiffness(),
            geometric: self.geometric_stiffness(),
            node_count: self.stress().nodal_internal_force_n().len(),
            identity: self.identity(),
            configuration: self.configuration_identity(),
        }
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

    /// Fixed-preload buckling for TET10, HEX8, HEX20 and HEX27. This state is
    /// lambda=0; `reference` is lambda=1. Only their geometric-stiffness
    /// difference scales, while this state's complete tangent stays fixed.
    ///
    /// Geometry, original/local material tensors and sources, quadrature,
    /// strain formulation and incremental supports must match. Loads and
    /// prescribed values may differ. An unstable/singular base refuses through
    /// the original SPD buckling owner; no negative stiffness is regularized.
    pub fn incremental_buckling(
        &self,
        reference: &Self,
        modes: usize,
        maximum_iterations: usize,
        cx: &Cx<'_>,
    ) -> Result<SolidBucklingResult3, SolidError> {
        solve_increment(
            self.buckling_reference(),
            reference.buckling_reference(),
            modes,
            maximum_iterations,
            cx,
        )
    }

    fn buckling_reference(&self) -> BucklingReference<'_> {
        BucklingReference {
            pencil: self.pencil(),
            elastic: self.elastic_stiffness(),
            geometric: self.geometric_stiffness(),
            node_count: self.response().displacement_m().len(),
            identity: self.identity(),
            configuration: self.configuration_identity(),
        }
    }
}

struct BucklingReference<'a> {
    pencil: &'a TetElasticAssembly,
    elastic: &'a Csr,
    geometric: &'a Csr,
    node_count: usize,
    identity: ContentHash,
    configuration: ContentHash,
}

fn solve_increment(
    base: BucklingReference<'_>,
    reference: BucklingReference<'_>,
    modes: usize,
    maximum_iterations: usize,
    cx: &Cx<'_>,
) -> Result<SolidBucklingResult3, SolidError> {
    cx.checkpoint().map_err(|_| SolidError::Cancelled)?;
    if base.configuration != reference.configuration
        || base.node_count != reference.node_count
        || base.pencil.free_dofs != reference.pencil.free_dofs
        || base.elastic != reference.elastic
    {
        return Err(invalid(
            "incremental buckling requires identical geometry, material sources/tensors, quadrature, formulation and incremental supports",
        ));
    }
    let increment = geometric_increment(base.geometric, reference.geometric, cx)?;
    // Cholesky admission of the complete base tangent precedes even the zero-
    // increment return in the original owner. A missing load direction must
    // never disguise an already unstable base as an unconditionally safe body.
    let mut result = solve(
        base.pencil,
        &base.pencil.stiffness,
        &increment,
        base.node_count,
        reference.identity,
        modes,
        maximum_iterations,
        cx,
    )?;
    result.base_assembly_identity = Some(base.identity);
    Ok(result)
}

// Merge canonical CSR rows in linear work; no dense intermediate or triplet
// sorting is needed to retain the signed reference-minus-base stress operator.
fn geometric_increment(base: &Csr, reference: &Csr, cx: &Cx<'_>) -> Result<Csr, SolidError> {
    let n = base.nrows();
    if base.ncols() != n || reference.nrows() != n || reference.ncols() != n {
        return Err(invalid("incremental geometric stiffness dimensions differ"));
    }
    let mut row_ptr = Vec::with_capacity(n + 1);
    let mut columns = Vec::new();
    let mut values = Vec::new();
    row_ptr.push(0);
    for row in 0..n {
        cx.checkpoint().map_err(|_| SolidError::Cancelled)?;
        let (a, av) = reference.row(row);
        let (b, bv) = base.row(row);
        let (mut i, mut j) = (0, 0);
        while i < a.len() || j < b.len() {
            let (column, value) = if j == b.len() || (i < a.len() && a[i] < b[j]) {
                let entry = (a[i], av[i]);
                i += 1;
                entry
            } else if i == a.len() || b[j] < a[i] {
                let entry = (b[j], -bv[j]);
                j += 1;
                entry
            } else {
                let entry = (a[i], av[i] - bv[j]);
                i += 1;
                j += 1;
                entry
            };
            if !value.is_finite() {
                return Err(invalid(
                    "incremental geometric stiffness is unrepresentable",
                ));
            }
            if value != 0.0 {
                columns.push(column);
                values.push(value);
            }
        }
        row_ptr.push(columns.len());
    }
    Csr::try_from_parts_with_checkpoint(n, n, row_ptr, columns, values, || {
        cx.checkpoint().map_err(|_| SolidError::Cancelled)
    })?
    .ok_or_else(|| invalid("incremental geometric stiffness is not canonical"))
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
        base_assembly_identity: None,
        reference_assembly_identity,
        node_count,
        result,
    })
}

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}
