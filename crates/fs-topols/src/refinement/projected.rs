//! A new, explicitly funded finer-grid study from an accepted volume-only design.
//! The source is borrowed. Prolongation and area restoration are preparation,
//! never an accepted descent step or an observed continuum error certificate.
use std::ops::ControlFlow;

use fs_cutfem::{CutFemError, Quadtree};
use crate::projected::{ProjectedOptimizer, ProjectedSetupStage};
use crate::{EvaluatedFinalState, OptimizeSettings, optimize::material_volume};
use super::prolongate_level_set;

/// Cooperative boundaries; one prolongation/area quadrature remains indivisible.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VolumeRefinementStage {
    Prepare,
    Prolongate,
    Area,
    Baseline(ProjectedSetupStage),
    Publish,
}

/// Measurements belong to different discretizations and are NOT an improvement.
#[derive(Debug, Clone, Copy)]
pub struct VolumeRefinementReport {
    pub coarse_level: u32,
    pub fine_level: u32,
    pub source_updates: usize,
    pub coarse_endpoint: EvaluatedFinalState,
    /// Fine-grid cut-quadrature area of the transferred field BEFORE projection.
    pub transferred_area: f64,
    /// Largest nodal field correction from projection, not a geometric distance.
    pub max_projection_change: f64,
    pub fine_baseline: EvaluatedFinalState,
}

/// Prolongate the actual accepted geometry, restore the SAME area policy, and
/// independently solve a one-level-finer baseline under the SAME load/material.
///
/// `updates` funds a NEW study. Its ordinal and AL/nucleation schedule restart;
/// no old work is refunded and no cross-grid compliance change is counted as
/// descent. Fixed values use the existing prolongator's nonzero-weight support.
/// Search controls and area tolerance are inherited without relaxation.
/// Cancellation returns no optimizer, including after the complete fine solve.
///
/// # Errors
/// Invalid new work, excessive level, unattainable area, or failed mechanics.
pub fn refine_controlled<B>(
    coarse: &ProjectedOptimizer,
    updates: usize,
    mut control: impl FnMut(VolumeRefinementStage) -> ControlFlow<B>,
) -> Result<ControlFlow<B, (ProjectedOptimizer, VolumeRefinementReport)>, CutFemError> {
    let prior = coarse.checkpoint().settings();
    if !(1..=7).contains(&prior.level) || !(1..=10_000).contains(&updates) {
        return Err(super::invalid("volume refinement requires level 1..=7 and 1..=10000 new updates"));
    }
    if let ControlFlow::Break(reason) = control(VolumeRefinementStage::Prepare) {
        return Ok(ControlFlow::Break(reason));
    }
    if let ControlFlow::Break(reason) = control(VolumeRefinementStage::Prolongate) {
        return Ok(ControlFlow::Break(reason));
    }
    let transfer = prolongate_level_set(coarse.checkpoint().geometry(), coarse.fixed_nodes())?;
    let level = prior.level + 1;
    if let ControlFlow::Break(reason) = control(VolumeRefinementStage::Area) {
        return Ok(ControlFlow::Break(reason));
    }
    let transferred_area = material_volume(&Quadtree::uniform(level), &transfer.geometry);
    if !(transferred_area.is_finite() && transferred_area > 0.0) {
        return Err(super::invalid("transferred design has no finite positive numerical area"));
    }
    let fine = match ProjectedOptimizer::new_controlled(
        &transfer.geometry, coarse.checkpoint().fixture(),
        OptimizeSettings { level, iterations: updates, ..prior },
        transfer.fixed_nodes, coarse.projection_settings(), coarse.controls(),
        |stage| control(VolumeRefinementStage::Baseline(stage)),
    )? {
        ControlFlow::Continue(fine) => fine,
        ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
    };
    let max_projection_change = transfer.geometry.nodes().iter()
        .zip(fine.checkpoint().geometry().nodes())
        .map(|(a, b)| (a - b).abs()).fold(0.0_f64, f64::max);
    if !max_projection_change.is_finite() {
        return Err(super::invalid("refinement projection correction overflowed"));
    }
    let report = VolumeRefinementReport {
        coarse_level: prior.level, fine_level: level,
        source_updates: coarse.checkpoint().next_iteration(),
        coarse_endpoint: coarse.current(), transferred_area, max_projection_change,
        fine_baseline: fine.current(),
    };
    if let ControlFlow::Break(reason) = control(VolumeRefinementStage::Publish) {
        return Ok(ControlFlow::Break(reason));
    }
    Ok(ControlFlow::Continue((fine, report)))
}

#[cfg(test)]
mod tests;
