//! Dyadic refinement of a projected design into a NEW, explicitly funded study.
//!
//! Bilinear prolongation retains the coarse field (up to interpolation rounding),
//! not a sampled binary mask or a remeshed boundary. The finer cut quadrature can
//! measure a different area, so the inherited area policy is applied again and
//! the resulting field is independently solved under every load. That field is
//! the fine study's baseline. A cross-grid objective change is NOT a descent step
//! or an error bound. Exact checkpoint continuation remains a separate operation.

use crate::robust_descent::{MultiLoadProjectedOptimizer, MultiLoadProjectedState};
use crate::{GridSdf, optimize::material_volume};
use fs_cutfem::{CutFemError, Quadtree};

/// Prolongated field and prescribed nodes on the next dyadic lattice.
#[derive(Debug, Clone)]
pub struct RefinedLevelSet {
    /// Bilinear prolongation; every old nodal bit is retained at its nested node.
    pub geometry: GridSdf,
    /// Sorted fixed nodes. A new node is fixed exactly when all coarse nodes
    /// with nonzero interpolation weights were fixed. Isolated fixed nodes do
    /// not freeze an adjacent cell, while fixed edges/regions remain fixed.
    pub fixed_nodes: Vec<(usize, f64)>,
}

/// Measured handoff between two separately baselined discretizations.
#[derive(Debug, Clone)]
pub struct ProjectedRefinementReport {
    pub coarse_level: u32,
    pub fine_level: u32,
    /// Completed updates and spent solves in the unchanged source study.
    pub source_updates: usize,
    pub source_solves_started: usize,
    /// Coarse-grid result; never used as the fine-grid descent threshold.
    pub coarse_endpoint: MultiLoadProjectedState,
    /// Fine-quadrature area BEFORE feasibility restoration, without a PDE solve.
    pub transferred_area: f64,
    /// Max nodal correction from area restoration. This is a field-value change,
    /// not a Hausdorff distance or a geometric error certificate.
    pub max_projection_change: f64,
    /// Fully solved fine-grid baseline AFTER area restoration.
    pub fine_baseline: MultiLoadProjectedState,
}

fn invalid(message: &'static str) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: message.into() }
}

/// Refine a dyadic lattice once without changing its bilinear representation.
///
/// Old nodes are injected bitwise. Edge midpoints use `f64::midpoint`; cell
/// centres use row-midpoints followed by their midpoint, fixing the rounding
/// order and avoiding overflowing sums. Fixed values must already match the
/// coarse field exactly; this operation does not repair prescribed data.
///
/// # Errors
/// Refuses non-dyadic/out-of-envelope lattices, non-finite fields, malformed
/// fixed-node lists, and inconsistent fixed values BEFORE allocating fine data.
pub fn prolongate_level_set(
    coarse: &GridSdf,
    fixed: &[(usize, f64)],
) -> Result<RefinedLevelSet, CutFemError> {
    let n = coarse.n();
    if !n.is_power_of_two() || !(2..=128).contains(&n) {
        return Err(invalid("refinement requires a dyadic coarse grid with 2..=128 cells per side"));
    }
    if coarse.nodes().iter().any(|v| !v.is_finite()) {
        return Err(invalid("refinement requires finite coarse nodal values"));
    }
    if fixed.len() > coarse.nodes().len() {
        return Err(invalid("refinement fixed-node count exceeds the coarse lattice"));
    }
    let mut previous = None;
    for &(index, value) in fixed {
        if index >= coarse.nodes().len() || !value.is_finite()
            || previous.is_some_and(|old| index <= old)
        {
            return Err(invalid("refinement fixed nodes must be finite, sorted, unique and in range"));
        }
        if coarse.nodes()[index].to_bits() != value.to_bits() {
            return Err(invalid("refinement fixed-node values do not match the retained field"));
        }
        previous = Some(index);
    }
    let mut pinned = vec![false; coarse.nodes().len()];
    for &(index, _) in fixed { pinned[index] = true; }
    let fine_n = 2 * n;
    let mut geometry = GridSdf::from_fn(fine_n, &|_, _| 0.0);
    let mut fixed_nodes = Vec::new();
    for j in 0..=fine_n {
        for i in 0..=fine_n {
            let (a, b, c, d) = (i / 2, i.div_ceil(2), j / 2, j.div_ceil(2));
            let v00 = coarse.node(a, c);
            let value = match (i % 2, j % 2) {
                (0, 0) => v00,
                (1, 0) => v00.midpoint(coarse.node(b, c)),
                (0, 1) => v00.midpoint(coarse.node(a, d)),
                _ => v00.midpoint(coarse.node(b, c))
                    .midpoint(coarse.node(a, d).midpoint(coarse.node(b, d))),
            };
            *geometry.node_mut(i, j) = value;
            let stride = n + 1;
            // For even coordinates some indices coincide; checking them again
            // does not introduce a zero-weight neighbour into fixed support.
            if [a + c * stride, b + c * stride, a + d * stride, b + d * stride]
                .iter().all(|&index| pinned[index])
            {
                fixed_nodes.push((i + j * (fine_n + 1), value));
            }
        }
    }
    Ok(RefinedLevelSet { geometry, fixed_nodes })
}

/// Fork a coarse endpoint into a one-level-finer, independently funded study.
///
/// `updates` and `max_solves` explicitly declare NEW work; the source study is
/// borrowed and is unchanged on both success and refusal. Loads/order, material,
/// objective, candidate controls, fixed regions and area/stress policies are
/// inherited. Area is restored before all baseline loads and optional stresses
/// are checked. Strict mode refuses an overstressed fine baseline; an explicitly
/// inherited restoration policy retains it as infeasible, never enlarging the
/// limit. No unassessed or merely interpolated displacement is published.
///
/// The new baseline, update ordinal and nucleation/AL schedule start a NEW study;
/// its checkpoint uses the existing format. No old candidate budget is refunded,
/// and no improvement is attributed to changing grids or restoring area.
/// Construction, like the existing baseline builder, is synchronous. A failed
/// construction returns no optimizer and may have spent baseline solve work.
///
/// # Errors
/// Refuses invalid new work, unsupported refinement, projection or physics.
pub fn refine_projected_study(
    coarse: &MultiLoadProjectedOptimizer,
    updates: usize,
    max_solves: usize,
) -> Result<(MultiLoadProjectedOptimizer, ProjectedRefinementReport), CutFemError> {
    let previous = coarse.settings();
    if !(1..=7).contains(&previous.level) || !(1..=10_000).contains(&updates)
        || max_solves < coarse.load_cases().len()
    {
        return Err(invalid("refinement needs coarse level 1..=7, 1..=10000 new updates and a complete new baseline solve allowance"));
    }
    let refined = prolongate_level_set(coarse.geometry(), coarse.fixed_nodes())?;
    let level = previous.level + 1;
    let transferred_area = material_volume(&Quadtree::uniform(level), &refined.geometry);
    if !(transferred_area.is_finite() && transferred_area > 0.0) {
        return Err(invalid("refined field has no finite positive numerical material area"));
    }
    let transferred = refined.geometry.clone();
    let settings = crate::OptimizeSettings { level, iterations: updates, ..previous };
    let controls = crate::robust_descent::MultiLoadProjectedSettings {
        max_solves, ..coarse.controls()
    };
    let mut fine = MultiLoadProjectedOptimizer::new(
        refined.geometry, coarse.load_cases(), settings, coarse.aggregate(),
        refined.fixed_nodes, coarse.projection_settings(), controls,
    )?;
    if let Some(limit) = coarse.stress_limit() {
        fine = match coarse.stress_restoration_reduction() {
            Some(reduction) => fine.with_stress_restoration(limit, reduction)?,
            None => fine.with_sampled_stress_limit(limit)?,
        };
    }
    let max_projection_change = transferred.nodes().iter().zip(fine.geometry().nodes())
        .map(|(a, b)| (a - b).abs()).fold(0.0_f64, f64::max);
    if !max_projection_change.is_finite() {
        return Err(invalid("refinement area-restoration displacement overflowed"));
    }
    let report = ProjectedRefinementReport {
        coarse_level: previous.level, fine_level: level,
        source_updates: coarse.next_iteration(), source_solves_started: coarse.solves_started(),
        coarse_endpoint: coarse.current(), transferred_area, max_projection_change,
        fine_baseline: fine.current(),
    };
    Ok((fine, report))
}

#[cfg(test)]
mod tests;
