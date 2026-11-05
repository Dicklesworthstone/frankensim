//! Existing sparse direct inversion of the original consistent reference mass.
//!
//! No stiffness or density averaging enters this adapter. The exact element
//! kernel is assembled in original element order, then normalized by the
//! positive reference row sums. For positive P1 triangle density, the local
//! forms satisfy D/5 <= M <= D; assembly and fixed-coordinate elimination
//! preserve that bound. Scaling therefore changes coordinates without making
//! a light physical inventory a small pivot solely because another is heavy.

use super::{AxisInverse, MassFactor, MembraneError, MembranePatch, ScaledFactor, poll};
use fs_sparse::{
    Coo,
    direct::{DirectOrdering, LdltError, LdltOptions, SymbolicLdlt},
};
use std::sync::Arc;

/// Explicit allowances for one scalar sparse reference-mass preparation.
/// At most three distinct factors are retained; equal support patterns share
/// one. Limits on each factor are NOT a cap on whole-process RSS or on the
/// existing patch/nonlinear integration buffers.
#[derive(Clone, Copy, Debug)]
pub struct SparseConsistentMassBudget {
    /// Maximum free scalar nodes in each Cartesian component, admitted before
    /// assembly and the existing synchronous ordering/symbolic analysis.
    pub max_nodes_per_axis: usize,
    /// Maximum original element contributions before duplicate COO entries
    /// are assembled. Admitted before constructing the sparse matrix.
    pub max_staged_entries_per_axis: usize,
    /// Retained symbolic payload, checked AFTER original analysis and before
    /// numeric work. This does not cap symbolic-construction peak memory.
    pub max_symbolic_bytes_per_factor: usize,
    /// Original direct owner's numeric payload cap, covering factor, fronts,
    /// updates and numeric scratch before their admitted allocations.
    pub max_factor_bytes_per_factor: usize,
    /// Maximum nnz(L), including the unit diagonal. Checked on both symbolic
    /// and actual factors; it bounds work of each warm scalar triangular solve.
    pub max_factor_nonzeros_per_axis: usize,
    /// Existing owner's no-delay structural flop estimate, admitted before
    /// numeric factorization. This is an estimate, not a wall-clock deadline
    /// or a mid-factorization interruptible operation counter.
    pub max_estimated_factor_flops_per_factor: u64,
}

impl SparseConsistentMassBudget {
    pub(super) fn admit(self) -> Result<(), MembraneError> {
        if self.max_nodes_per_axis == 0
            || self.max_staged_entries_per_axis == 0
            || self.max_symbolic_bytes_per_factor == 0
            || self.max_factor_bytes_per_factor == 0
            || self.max_factor_nonzeros_per_axis == 0
            || self.max_estimated_factor_flops_per_factor == 0
        {
            return Err(MembraneError::Budget(
                "positive sparse membrane inertia allowances required",
            ));
        }
        Ok(())
    }
}

/// Storage and work of the prepared scalar factors, from their original
/// direct-solver statistics. These are counts/payloads, not throughput claims.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SparseMassStatistics {
    /// Number of independently prepared support patterns, at most three.
    pub distinct_factors: usize,
    /// Sum of input CSR nonzeros over distinct factors, including both halves.
    pub scalar_mass_nonzeros: usize,
    /// Sum of nnz(L) over distinct factors, including unit diagonals.
    pub stored_factor_nonzeros: usize,
    /// Sum of nnz(L) over active axes, counting shared factors per use.
    /// Each inverse action performs two triangular traversals and O(free DOFs)
    /// permutation, diagonal and physical-coordinate scaling work.
    pub factor_nonzeros_per_inverse: usize,
    /// Largest actual scalar frontal dimension among the prepared factors.
    pub maximum_front_nodes: usize,
    /// Sum of original numeric peak-payload bounds over distinct factors.
    /// It conservatively covers preparing one while retaining earlier factors;
    /// CSR, symbolic construction, allocator overhead and patch buffers are
    /// outside that numeric-owner account.
    pub summed_numeric_peak_bytes: usize,
}

pub(super) fn factor(
    patch: &MembranePatch,
    free: &[usize],
    budget: SparseConsistentMassBudget,
    ordering: DirectOrdering,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ScaledFactor, MembraneError> {
    poll(cancelled)?;
    let n = free.len();
    let mut index = vec![usize::MAX; patch.mesh.nodes.len()];
    let mut root_mass = Vec::with_capacity(n);
    for (row, &dof) in free.iter().enumerate() {
        poll(cancelled)?;
        let node = dof / 3;
        index[node] = row;
        root_mass.push(patch.nodal_mass_kg[node].sqrt());
    }
    let mut staged = 0usize;
    for element in &patch.elements {
        poll(cancelled)?;
        let free_corners = element
            .nodes
            .iter()
            .filter(|&&node| index[node] != usize::MAX)
            .count();
        staged = staged
            .checked_add(free_corners * free_corners)
            .ok_or(MembraneError::Budget(
                "sparse membrane mass contribution extent",
            ))?;
        if staged > budget.max_staged_entries_per_axis {
            return Err(MembraneError::Budget(
                "sparse membrane mass contribution allowance",
            ));
        }
    }
    let mut coo = Coo::new(n, n);
    for (element_index, element) in patch.elements.iter().enumerate() {
        poll(cancelled)?;
        let mass = patch.element_consistent_mass(element_index)?;
        for a in 0..3 {
            let row = index[element.nodes[a]];
            if row == usize::MAX {
                continue;
            }
            for b in 0..3 {
                let column = index[element.nodes[b]];
                if column != usize::MAX {
                    coo.push(row, column, mass[a][b]);
                }
            }
        }
    }
    drop(index);
    // Original deterministic COO assembly keeps element contribution order.
    // It is synchronous; its complete input extent was admitted above.
    poll(cancelled)?;
    let mut mass = coo.assemble();
    drop(coo);
    poll(cancelled)?;
    for row in 0..n {
        poll(cancelled)?;
        let entries = mass.row(row).0.len();
        for entry in 0..entries {
            let column = mass.row(row).0[entry];
            let lo = root_mass[row].min(root_mass[column]);
            let hi = root_mass[row].max(root_mass[column]);
            let value = &mut mass.row_values_mut(row)[entry];
            // Same physical sum and division order as the dense image. Both
            // symmetric entries use the smaller root first, including at
            // density contrasts where multiplying the two roots would fail.
            *value = (*value / lo) / hi;
            if !value.is_finite() {
                return Err(MembraneError::Invalid(
                    "unrepresentable scaled sparse membrane inertia",
                ));
            }
        }
    }
    let coordinates = if matches!(ordering, DirectOrdering::NestedDissection(_)) {
        let mut points = Vec::with_capacity(n);
        for &dof in free {
            poll(cancelled)?;
            let [x, y] = patch.mesh.nodes[dof / 3];
            points.push([x, y, 0.0]);
        }
        Some(points)
    } else {
        None
    };
    poll(cancelled)?;
    let analysis = SymbolicLdlt::analyze_with_coordinates(&mass, ordering, coordinates.as_deref());
    poll(cancelled)?;
    let symbolic = analysis.map_err(direct_error)?;
    drop(coordinates);
    if symbolic.symbolic_bytes() > budget.max_symbolic_bytes_per_factor {
        return Err(MembraneError::Budget(
            "retained sparse membrane symbolic payload",
        ));
    }
    if symbolic.nnz_l() > budget.max_factor_nonzeros_per_axis {
        return Err(MembraneError::Budget(
            "sparse membrane factor nonzero allowance",
        ));
    }
    if symbolic.estimated_flops() > budget.max_estimated_factor_flops_per_factor {
        return Err(MembraneError::Budget(
            "estimated sparse membrane factor work",
        ));
    }
    let options = LdltOptions {
        allow_2x2: false,
        max_factor_bytes: Some(budget.max_factor_bytes_per_factor),
        ..LdltOptions::default()
    };
    poll(cancelled)?;
    let numeric = symbolic.factor(&mass, &options);
    poll(cancelled)?;
    let factor = numeric.map_err(direct_error)?;
    let signs = factor.inertia();
    if signs.positive != n || signs.negative != 0 {
        return Err(MembraneError::Invalid(
            "sparse membrane inertia is not numerically positive definite",
        ));
    }
    if factor.stats().nnz_l > budget.max_factor_nonzeros_per_axis {
        return Err(MembraneError::Budget(
            "actual sparse membrane factor nonzero allowance",
        ));
    }
    poll(cancelled)?;
    Ok(ScaledFactor {
        root_mass,
        factor: MassFactor::Sparse(factor),
    })
}

fn direct_error(error: LdltError) -> MembraneError {
    match error {
        LdltError::FactorBudget { .. } => {
            MembraneError::Budget("sparse membrane numeric factor payload")
        }
        LdltError::InvalidOptions { .. } => {
            MembraneError::Invalid("invalid sparse membrane ordering or factor options")
        }
        _ => MembraneError::Invalid("original sparse membrane mass factorization refused"),
    }
}

pub(super) fn statistics(axes: &[AxisInverse]) -> Option<SparseMassStatistics> {
    let mut result = SparseMassStatistics::default();
    for (i, axis) in axes.iter().enumerate() {
        let MassFactor::Sparse(factor) = &axis.factor.factor else {
            return None;
        };
        let stats = factor.stats();
        result.factor_nonzeros_per_inverse = result
            .factor_nonzeros_per_inverse
            .saturating_add(stats.nnz_l);
        if axes[..i]
            .iter()
            .any(|previous| Arc::ptr_eq(&previous.factor, &axis.factor))
        {
            continue;
        }
        result.distinct_factors += 1;
        result.scalar_mass_nonzeros = result.scalar_mass_nonzeros.saturating_add(stats.nnz_a);
        result.stored_factor_nonzeros = result.stored_factor_nonzeros.saturating_add(stats.nnz_l);
        result.maximum_front_nodes = result.maximum_front_nodes.max(stats.max_front_dim);
        result.summed_numeric_peak_bytes = result
            .summed_numeric_peak_bytes
            .saturating_add(stats.peak_factor_bytes);
    }
    Some(result)
}
