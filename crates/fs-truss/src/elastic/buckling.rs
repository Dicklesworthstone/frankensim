//! Global initial-stress buckling of the native compatible planar truss.
//!
//! Two ordinary elastic solves establish the fixed baseline and the homogeneous
//! reference-load increment. The pencil is
//! (K_e + K_g(N_fixed) + lambda K_g(N_reference)) phi = 0.
//! Thermal/assembly strain, support settlement and explicit fixed nodal loads
//! stay fixed while the reference case's actual nodal loads are multiplied.
//!
//! For a bar, K_g = (N/L) t t^T on relative endpoint translations, with t the
//! transverse direction. Tension stiffens and compression softens. Geometry,
//! source E, selected areas and support flags remain frozen. This is the
//! small-displacement initial-stress approximation, not a finite-strain
//! equilibrium path, follower-load tangent, postbuckling or member Euler check.
//!
//! Sparse assembly uses fs-sparse. Positive-base admission, scaling, buckling
//! eigenpairs and original-pencil residual checks use fs-solid::stability,
//! including its bounded sparse/fs-spectral dispatch for large pencils.

use super::{
    ElasticControls, ElasticError, ElasticImposedState, ElasticResponse, checkpoint, invalid,
};
use crate::{GroundStructure, LayoutCase, TrussConstructionError};
use fs_exec::Cx;
use fs_solid::SolidError;
use fs_solid::stability::sparse::SparseBucklingLimits;
use fs_solid::stability::{BucklingResult, buckling_loads_cancellable};
use fs_sparse::Coo;

/// Hard cap on combined triplets before native deterministic COO assembly.
/// This bounds the existing assembler's unpolled sort and temporary storage.
pub const HARD_MAX_BUCKLING_STAGED_ENTRIES: usize = 8_000_000;

/// Fixed baseline data in the complete original member/node ordering.
#[derive(Debug, Clone, Copy, Default)]
pub struct BucklingLoadFamily<'a> {
    /// Fixed member eigenstrain and prescribed support displacement.
    pub imposed: ElasticImposedState<'a>,
    /// Optional fixed physical nodal loads [N], such as already computed gravity.
    /// None means zero. These loads are added to lambda * reference_case.loads().
    pub fixed_loads_n: Option<&'a [[f64; 2]]>,
}

/// Bounded equilibrium, sparse assembly and native spectral work.
#[derive(Debug, Clone, Copy)]
pub struct BucklingControls {
    /// Shared CG allowance for the fixed and reference-load equilibrium solves.
    pub elastic: ElasticControls,
    /// Number of lowest positive load factors requested, capped by free DOFs.
    pub modes: usize,
    /// Native buckling iteration allowance, from 1 through 1000000.
    pub maximum_eigen_iterations: usize,
    /// Combined nonzero contributions staged for both pencil matrices.
    /// Native sparse input/factor/work limits additionally remain in force.
    pub maximum_staged_entries: usize,
}

impl Default for BucklingControls {
    fn default() -> Self {
        Self {
            elastic: ElasticControls::default(),
            modes: 1,
            maximum_eigen_iterations: 200,
            maximum_staged_entries: 4_000_000,
        }
    }
}

/// No accepted buckling result is returned after either owner's refusal.
#[derive(Debug)]
pub enum BucklingError {
    /// Shape/material/load/control admission, unresolved equilibrium or cancellation.
    Elastic(ElasticError),
    /// Native positive-base, spectral, resource, residual or cancellation refusal.
    /// A refusal alone does not distinguish instability from numerical inability.
    Stability(SolidError),
}

impl From<ElasticError> for BucklingError {
    fn from(error: ElasticError) -> Self {
        Self::Elastic(error)
    }
}
impl From<TrussConstructionError> for BucklingError {
    fn from(error: TrussConstructionError) -> Self {
        Self::Elastic(error.into())
    }
}
impl From<SolidError> for BucklingError {
    fn from(error: SolidError) -> Self {
        Self::Stability(error)
    }
}
impl core::fmt::Display for BucklingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "truss initial-stress buckling: {self:?}")
    }
}
impl std::error::Error for BucklingError {}

/// Fixed/reference equilibrium states and residual-qualified numerical buckling.
#[derive(Debug, Clone)]
pub struct TrussBuckling {
    /// Actual fixed loads, fixed support motion and fixed eigenstrain response.
    pub baseline: ElasticResponse,
    /// Response to reference_case.loads() with zero support motion/eigenstrain.
    /// At multiplier lambda, displacements and forces superpose as
    /// baseline + lambda * reference_load; stored energies are quadratic.
    pub reference_load: ElasticResponse,
    /// Native critical multipliers and original-node/component DOF map.
    ///
    /// Empty loads means no positive factor was returned for the queried
    /// reference-load directions. It is not a nonlinear or complete-spectrum
    /// stability certificate. Modes obey homogeneous support perturbations.
    pub buckling: BucklingResult,
}

/// Solve initial-stress buckling with a fixed prestress baseline.
///
/// All original free coordinates enter the pencil, including coordinates that
/// the ordinary elastic solve can omit when their stiffness and load are zero.
/// An unloaded mechanism therefore cannot disappear from stability admission.
/// The native owner requires a positive-definite baseline K_e + K_g(N_fixed)
/// and refuses singular/indefinite baselines without regularizing springs.
///
/// Reference nodal loads scale together; eigenstrain, support motion, optional
/// fixed loads and material properties stay fixed. This distinction matters
/// for gravity and thermally restrained assemblies. Original source geometry
/// supplies lengths/directions; no deformed-geometry or follower-load stiffness
/// is inferred. Local compression-member Euler capacity remains independent.
///
/// Both primal solves use the existing matrix-free checked CG, sharing the
/// elastic iteration budget. Member assembly polls every 256 members. Native
/// COO sorts are bounded by maximum_staged_entries and checked before/after;
/// factorization and eigensolves retain their existing cancellation and limits.
/// Numerical eigenpairs must pass the native original-pencil residual gate.
///
/// # Errors
/// Invalid data, no free coordinates, exhausted work/storage, either unresolved
/// equilibrium, non-SPD baseline, native spectral refusal and cancellation
/// return no partial result. Positive roots do not establish missed-mode
/// exclusion, uniqueness, mesh convergence or a nonlinear collapse load.
#[allow(clippy::too_many_lines)] // Two native states, physical sparse pencil, native spectrum.
pub fn solve(
    ground: &GroundStructure,
    reference_case: &LayoutCase,
    areas_m2: &[f64],
    young_pa: &[f64],
    family: BucklingLoadFamily<'_>,
    controls: BucklingControls,
    cx: &Cx<'_>,
) -> Result<TrussBuckling, BucklingError> {
    checkpoint(cx)?;
    if reference_case.node_count() != ground.nodes().len() {
        return Err(invalid("truss buckling case", "must match original node ordering").into());
    }
    if controls.modes == 0
        || controls.maximum_eigen_iterations == 0
        || controls.maximum_eigen_iterations > 1_000_000
        || controls.maximum_staged_entries == 0
        || controls.maximum_staged_entries > HARD_MAX_BUCKLING_STAGED_ENTRIES
    {
        return Err(invalid(
            "truss buckling controls",
            "require positive mode/work allowances within their declared bounds",
        )
        .into());
    }
    let native_limits = SparseBucklingLimits::default();
    let full_dofs = 2 * ground.nodes().len();
    if full_dofs > native_limits.max_dofs.saturating_mul(8) {
        return Err(TrussConstructionError::WorkBudget {
            resource: "truss buckling original DOFs",
            limit: native_limits.max_dofs.saturating_mul(8),
            observed: full_dofs,
        }
        .into());
    }
    let mut dof_map = vec![None; full_dofs];
    let mut free_count = 0;
    for (node, supported) in reference_case.supported().iter().enumerate() {
        if node % 256 == 0 {
            checkpoint(cx)?;
        }
        for (axis, &fixed) in supported.iter().enumerate() {
            if !fixed {
                dof_map[2 * node + axis] = Some(free_count);
                free_count += 1;
            }
        }
    }
    if free_count == 0 {
        return Err(invalid(
            "truss buckling free DOFs",
            "require at least one free coordinate",
        )
        .into());
    }
    if free_count > native_limits.max_dofs {
        return Err(TrussConstructionError::WorkBudget {
            resource: "truss buckling free DOFs",
            limit: native_limits.max_dofs,
            observed: free_count,
        }
        .into());
    }
    let fixed_loads = family.fixed_loads_n.map_or_else(
        || vec![[0.0; 2]; ground.nodes().len()],
        <[[f64; 2]]>::to_vec,
    );
    let baseline_case = LayoutCase::try_new(
        reference_case.supported().to_vec(),
        fixed_loads,
        ground.nodes().len(),
    )?;
    let baseline = super::solve_with_imposed(
        ground,
        &baseline_case,
        areas_m2,
        young_pa,
        family.imposed,
        controls.elastic,
        cx,
    )?;
    let reference_load = super::solve(
        ground,
        reference_case,
        areas_m2,
        young_pa,
        ElasticControls {
            maximum_iterations: controls
                .elastic
                .maximum_iterations
                .saturating_sub(baseline.iterations),
            ..controls.elastic
        },
        cx,
    )?;
    let mut base = Coo::new(free_count, free_count);
    let mut reference = Coo::new(free_count, free_count);
    for (index, &(a, b)) in ground.members().iter().enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        if areas_m2[index] == 0.0 {
            continue;
        }
        let length = ground.lengths()[index];
        let stiffness = young_pa[index] * (areas_m2[index] / length);
        let fixed_geometric = baseline.member_force_n[index] / length;
        let reference_geometric = reference_load.member_force_n[index] / length;
        if !fixed_geometric.is_finite() || !reference_geometric.is_finite() {
            return Err(invalid("truss geometric stiffness", "N/L must remain finite").into());
        }
        let dx = (ground.nodes()[b][0] - ground.nodes()[a][0]) / length;
        let dy = (ground.nodes()[b][1] - ground.nodes()[a][1]) / length;
        let axial = [-dx, -dy, dx, dy];
        let transverse = [dy, -dx, -dy, dx];
        let ends = [2 * a, 2 * a + 1, 2 * b, 2 * b + 1];
        for (i, &row_end) in ends.iter().enumerate() {
            let Some(row) = dof_map[row_end] else {
                continue;
            };
            for (j, &column_end) in ends.iter().enumerate() {
                let Some(column) = dof_map[column_end] else {
                    continue;
                };
                // Multiply direction products first: transposed entries have
                // identical floating-point contributions and reduction order.
                let transverse_product = transverse[i] * transverse[j];
                let fixed =
                    stiffness * (axial[i] * axial[j]) + fixed_geometric * transverse_product;
                let scaled = reference_geometric * transverse_product;
                if !fixed.is_finite() || !scaled.is_finite() {
                    return Err(
                        invalid("truss buckling pencil", "entries must remain finite").into(),
                    );
                }
                let staged = base.staged()
                    + reference.staged()
                    + usize::from(fixed != 0.0)
                    + usize::from(scaled != 0.0);
                if staged > controls.maximum_staged_entries {
                    return Err(TrussConstructionError::WorkBudget {
                        resource: "truss buckling staged entries",
                        limit: controls.maximum_staged_entries,
                        observed: staged,
                    }
                    .into());
                }
                if fixed != 0.0 {
                    base.push(row, column, fixed);
                }
                if scaled != 0.0 {
                    reference.push(row, column, scaled);
                }
            }
        }
    }
    checkpoint(cx)?;
    let base_csr = base.assemble();
    drop(base);
    checkpoint(cx)?;
    let reference_csr = reference.assemble();
    drop(reference);
    checkpoint(cx)?;
    let stored = base_csr.nnz() + reference_csr.nnz();
    if stored > native_limits.max_matrix_entries {
        return Err(TrussConstructionError::WorkBudget {
            resource: "truss buckling stored entries",
            limit: native_limits.max_matrix_entries,
            observed: stored,
        }
        .into());
    }
    let buckling = buckling_loads_cancellable(
        &base_csr,
        &reference_csr,
        &dof_map,
        controls.modes,
        controls.maximum_eigen_iterations,
        &|| cx.checkpoint().is_err(),
    )?;
    checkpoint(cx)?;
    Ok(TrussBuckling {
        baseline,
        reference_load,
        buckling,
    })
}
