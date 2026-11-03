//! Original quadratic volume maps over complete declared nodal paths.
//!
//! Boundary clearance does not admit a folded material map. This adapter checks
//! every original TET10 element over its complete closed material tetrahedron and
//! the complete affine coordinate path. The original reference map is admitted
//! once and bound by a borrow of the immutable surface topology. These local
//! orientation proofs do not establish global injectivity or prevent distinct
//! material regions from overlapping.

use super::QuadraticSolidContactSurface;
use fs_exec::Cx;
use fs_query::tetrahedron::quadratic::{
    QuadraticTetPathControls, QuadraticTetPathError, QuadraticTetPathVerdict,
    quadratic_tetrahedron_orientation,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum QuadraticSolidVolumePathVerdict {
    /// All original material maps stay positively oriented over the whole path.
    /// The minimum uses each element's own power-of-two normalization length;
    /// it is a dimensionless bound, not a physical determinant in cubic metres.
    Positive { minimum_scaled_determinant: f64 },
    /// Earliest unresolved time window. The retained material bounds identify
    /// the region that could not be admitted. This is not an asserted fold time.
    Unresolved {
        tetrahedron: usize,
        fraction_start: f64,
        fraction_end: f64,
        barycentric_bounds: [[f64; 2]; 4],
        scaled_determinant_bounds: [f64; 2],
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadraticSolidVolumePathReport {
    pub verdict: QuadraticSolidVolumePathVerdict,
    pub tested_tetrahedra: usize,
    /// Aggregate evaluated subdivision cells over the examined elements.
    pub examined_cells: usize,
}

#[derive(Debug)]
pub enum QuadraticSolidVolumeError {
    Invalid(&'static str),
    Query(QuadraticTetPathError),
    CellBudgetExhausted {
        maximum_cells: usize,
    },
    /// A reference mesh with unresolved regularity cannot authorize motion.
    ReferenceUnresolved {
        report: QuadraticSolidVolumePathReport,
    },
    Cancelled,
}

impl From<QuadraticTetPathError> for QuadraticSolidVolumeError {
    fn from(error: QuadraticTetPathError) -> Self {
        Self::Query(error)
    }
}
impl core::fmt::Display for QuadraticSolidVolumeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "quadratic solid volume path refused: {self:?}")
    }
}
impl std::error::Error for QuadraticSolidVolumeError {}

/// An admitted original reference mesh bound to its immutable native incidence.
/// Preparing this object copies neither coordinates nor constitutive histories.
#[derive(Debug)]
pub struct PreparedQuadraticSolidVolume<'s> {
    surface: &'s QuadraticSolidContactSurface,
    reference: QuadraticSolidVolumePathReport,
}

impl QuadraticSolidContactSurface {
    /// Admit the entire original reference map once before checking motion.
    /// Cell work is shared across the complete original element set. A frontier
    /// is held for only one element at a time; its bound remains unchanged.
    pub fn prepare_volume_paths(
        &self,
        controls: QuadraticTetPathControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedQuadraticSolidVolume<'_>, QuadraticSolidVolumeError> {
        let reference = query_volume_paths(
            self,
            &self.reference_nodes_m,
            &self.reference_nodes_m,
            controls,
            cx,
        )?;
        if !matches!(
            reference.verdict,
            QuadraticSolidVolumePathVerdict::Positive { .. }
        ) {
            return Err(QuadraticSolidVolumeError::ReferenceUnresolved { report: reference });
        }
        Ok(PreparedQuadraticSolidVolume {
            surface: self,
            reference,
        })
    }
}

impl PreparedQuadraticSolidVolume<'_> {
    pub fn reference_admission(&self) -> &QuadraticSolidVolumePathReport {
        &self.reference
    }

    /// Check every original node's declared affine motion. A successful positive
    /// verdict covers all material coordinates, including edges and vertices;
    /// accepted quadrature samples or positive endpoint corners cannot replace it.
    /// Limits apply to this complete trial and never grow during subdivision.
    pub fn affine_path_at_positions(
        &self,
        positions_start: &[[f64; 3]],
        positions_end: &[[f64; 3]],
        controls: QuadraticTetPathControls,
        cx: &Cx<'_>,
    ) -> Result<QuadraticSolidVolumePathReport, QuadraticSolidVolumeError> {
        query_volume_paths(self.surface, positions_start, positions_end, controls, cx)
    }
}

fn query_volume_paths(
    surface: &QuadraticSolidContactSurface,
    start: &[[f64; 3]],
    end: &[[f64; 3]],
    controls: QuadraticTetPathControls,
    cx: &Cx<'_>,
) -> Result<QuadraticSolidVolumePathReport, QuadraticSolidVolumeError> {
    checkpoint(cx)?;
    controls.validate()?;
    if start.len() != surface.reference_nodes_m.len()
        || end.len() != surface.reference_nodes_m.len()
        || surface.tetrahedra.is_empty()
    {
        return Err(QuadraticSolidVolumeError::Invalid(
            "complete matching original quadratic coordinate images required",
        ));
    }
    for (index, node) in start.iter().chain(end).enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        if node.iter().any(|value| !value.is_finite()) {
            return Err(QuadraticSolidVolumeError::Invalid(
                "finite quadratic path coordinates required",
            ));
        }
    }
    let mut report = QuadraticSolidVolumePathReport {
        verdict: QuadraticSolidVolumePathVerdict::Positive {
            minimum_scaled_determinant: f64::INFINITY,
        },
        tested_tetrahedra: 0,
        examined_cells: 0,
    };
    let mut minimum = f64::INFINITY;
    let mut unresolved: Option<QuadraticSolidVolumePathVerdict> = None;
    for (tetrahedron, nodes) in surface.tetrahedra.iter().enumerate() {
        checkpoint(cx)?;
        let remaining = controls.maximum_cells.saturating_sub(report.examined_cells);
        if remaining == 0 {
            return Err(QuadraticSolidVolumeError::CellBudgetExhausted {
                maximum_cells: controls.maximum_cells,
            });
        }
        let mut limits = controls;
        limits.maximum_cells = remaining;
        limits.maximum_active_cells = limits.maximum_active_cells.min(remaining);
        let element = quadratic_tetrahedron_orientation(
            nodes.map(|node| start[node]),
            nodes.map(|node| end[node]),
            limits,
            cx,
        )?;
        report.tested_tetrahedra += 1;
        report.examined_cells += element.examined_cells;
        match element.verdict {
            QuadraticTetPathVerdict::Positive {
                minimum_scaled_determinant,
            } => {
                minimum = minimum.min(minimum_scaled_determinant);
            }
            QuadraticTetPathVerdict::Unresolved {
                fraction_start,
                fraction_end,
                barycentric_bounds,
                scaled_determinant_bounds,
            } => {
                let earlier = match unresolved {
                    Some(QuadraticSolidVolumePathVerdict::Unresolved {
                        fraction_start: previous,
                        ..
                    }) => fraction_start < previous,
                    _ => true,
                };
                if earlier {
                    unresolved = Some(QuadraticSolidVolumePathVerdict::Unresolved {
                        tetrahedron,
                        fraction_start,
                        fraction_end,
                        barycentric_bounds,
                        scaled_determinant_bounds,
                    });
                }
                // No later element can precede time zero. Ascending element IDs
                // already select the canonical lowest index for this tie.
                if fraction_start == 0.0 {
                    break;
                }
            }
        }
    }
    report.verdict = unresolved.unwrap_or(QuadraticSolidVolumePathVerdict::Positive {
        minimum_scaled_determinant: minimum,
    });
    checkpoint(cx)?;
    Ok(report)
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), QuadraticSolidVolumeError> {
    cx.checkpoint()
        .map_err(|_| QuadraticSolidVolumeError::Cancelled)
}

#[cfg(test)]
mod tests;
