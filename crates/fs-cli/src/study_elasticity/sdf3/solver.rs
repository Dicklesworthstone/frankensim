//! Native selection of the existing geometry-owned elasticity preparations.
//! Correction geometry supplies interpolation only; current fine-density cell
//! and ghost terms supply every Galerkin matrix. No stale-factor reuse.
use super::*;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveMultilevelOptions3;
use fs_solver::op::multilevel::{MultilevelBudget, MultilevelError};

#[derive(Debug, Clone, Copy, Default)]
pub(super) enum Policy {
    /// Preserve existing native optimization and enrichment arithmetic.
    #[default]
    Legacy,
    Multilevel {
        coarsest_level: u32,
        transfer_terms: usize,
        matrix_entries: usize,
        galerkin_products: usize,
        diagonal_contributions: usize,
    },
}
impl Policy {
    pub(super) fn validate(self, initial_level: u32) -> Result<()> {
        if let Self::Multilevel { coarsest_level, transfer_terms, matrix_entries,
            galerkin_products, diagonal_contributions } = self {
            if coarsest_level >= initial_level
                || !(1..=4_000_000).contains(&transfer_terms)
                || !(1..=8_000_000).contains(&matrix_entries)
                || !(1..=1_000_000_000).contains(&galerkin_products)
                || !(1..=100_000_000).contains(&diagonal_contributions)
            {
                return Err(fail("cli-study-sdf3-input",
                    "multilevel requires a coarsest level below the initial grid, 1..=4000000 transfer terms, 1..=8000000 matrix entries, 1..=1000000000 setup products and 1..=100000000 diagonal contributions"));
            }
        }
        Ok(())
    }

    /// Additional structural admission envelope, not measured allocator RSS.
    /// Include simultaneous retained transfer copies, sparse setup accumulators,
    /// factors and the correction geometries used throughout the computation.
    pub(super) fn memory_envelope(self, initial_level: u32) -> usize {
        match self {
            Self::Legacy => 0,
            Self::Multilevel { coarsest_level, transfer_terms, matrix_entries, .. } => {
                let cells: usize = (coarsest_level..initial_level).map(|level| 8usize.pow(level)).sum();
                192 * matrix_entries + 128 * transfer_terms + 65_536 * cells
            }
        }
    }

    fn options(self) -> Option<AdaptiveMultilevelOptions3> {
        match self {
            Self::Legacy => None,
            Self::Multilevel { transfer_terms, matrix_entries, galerkin_products,
                diagonal_contributions, .. } => Some(AdaptiveMultilevelOptions3 {
                max_transfer_terms: transfer_terms,
                max_diagonal_contributions: diagonal_contributions,
                hierarchy: MultilevelBudget {
                    max_fine_dofs: 50_000,
                    max_levels: 8,
                    max_transfer_entries: transfer_terms,
                    max_matrix_entries: matrix_entries,
                    max_galerkin_products: galerkin_products,
                    max_coarsest_dofs: 512,
                },
            }),
        }
    }

    pub(super) fn enrichment(self) -> GoalPreconditioner3 {
        match self.options() {
            Some(options) => GoalPreconditioner3::Multilevel { options },
            None => GoalPreconditioner3::TwoLevel {
                budget: TwoLevelBudget { max_fine_dofs: 50_000, ..Default::default() },
                max_diagonal_contributions: 100_000_000,
            },
        }
    }

    pub(super) fn wrap(self, operator: AdaptiveElasticity3, coarser: &[&AdaptiveElasticity3],
        checkpoint: impl FnMut() -> ControlFlow<()>)
        -> std::result::Result<AdaptiveSolveSpace3, GoalRefinementError3> {
        match self.options() {
            None => {
                if !coarser.is_empty() {
                    return Err(GoalRefinementError3::Invalid("unexpected correction geometry for legacy solver"));
                }
                Ok(AdaptiveSolveSpace3::jacobi(operator, 100_000_000))
            }
            Some(options) => AdaptiveSolveSpace3::multilevel(operator, coarser, options, checkpoint)
                .map_err(GoalRefinementError3::from),
        }
    }
}

/// Construct the immutable uniform ladder once, nearest first. The builder is
/// the SAME native domain/support/load-admission path used by the fine grid and
/// must charge these calls to the original shared quadrature allowance.
pub(super) fn correction_spaces(
    spec: &Spec,
    build: &mut impl FnMut(&Octree3, &mut dyn FnMut() -> ControlFlow<()>)
        -> std::result::Result<AdaptiveElasticity3, GoalRefinementError3>,
    checkpoint: &mut dyn FnMut() -> ControlFlow<()>,
) -> std::result::Result<Vec<AdaptiveElasticity3>, GoalRefinementError3> {
    let Policy::Multilevel { coarsest_level, .. } = spec.solver else { return Ok(Vec::new()); };
    let mut geometries = Vec::new();
    for level in (coarsest_level..spec.level).rev() {
        if checkpoint().is_break() { return Err(EvaluationStop::Cancelled.into()); }
        let tree = Octree3::uniform(level as u8, spec.max_level as u8, spec.leaves)
            .map_err(|_| GoalRefinementError3::Invalid("correction grid exceeds the study background envelope"))?;
        geometries.push(build(&tree, checkpoint)?);
    }
    if checkpoint().is_break() { return Err(EvaluationStop::Cancelled.into()); }
    Ok(geometries)
}

pub(super) fn setup_budget(error: &GoalRefinementError3) -> bool {
    matches!(error, GoalRefinementError3::Preconditioner(
        AdaptivePreconditionError3::Coarse(TwoLevelError::Budget(_))
            | AdaptivePreconditionError3::Hierarchy(MultilevelError::Budget(_))))
}

#[cfg(test)]
#[path = "solver_tests.rs"]
mod tests;
