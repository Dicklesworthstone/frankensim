//! Authored keep-material and keep-void regions for the existing fixed-node path.
//!
//! Every positive-area cell intersecting an authored rectangle is covered. All
//! its corners are prescribed with the requested sign margin, so its exact
//! bilinear interpolant has that sign everywhere, not just at cell centres.
//! Coverage can extend beyond the requested rectangle by less than one cell per
//! side. This is discrete geometry authoring, not a manufacturing certificate.

use crate::GridSdf;
use fs_cutfem::CutFemError;
use std::convert::Infallible;
use std::ops::ControlFlow;

/// Sign that must survive all projected geometry updates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DesignPhase {
    /// Negative level-set values: retain material.
    Material,
    /// Positive level-set values: retain empty clearance.
    Void,
}

/// One axis-aligned rectangle in the normalized unit square.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DesignRegion {
    phase: DesignPhase,
    lower: [f64; 2],
    upper: [f64; 2],
    margin: f64,
}

impl DesignRegion {
    /// Declare positive-area bounds and a finite, strictly positive field margin.
    /// The margin has the same units as phi; it is NOT a geometric clearance.
    ///
    /// # Errors
    /// Refuses unordered/non-finite/out-of-domain bounds or an invalid margin.
    pub fn new(
        phase: DesignPhase, lower: [f64; 2], upper: [f64; 2], margin: f64,
    ) -> Result<Self, CutFemError> {
        if !(margin.is_finite() && margin > 0.0)
            || (0..2).any(|axis| !(lower[axis].is_finite() && upper[axis].is_finite()
                && lower[axis] >= 0.0 && lower[axis] < upper[axis] && upper[axis] <= 1.0))
        {
            return Err(invalid("design regions require ordered bounds in [0,1] and a finite positive phi margin"));
        }
        Ok(Self { phase, lower, upper, margin })
    }

    #[must_use]
    pub const fn phase(self) -> DesignPhase { self.phase }
    #[must_use]
    pub const fn lower(self) -> [f64; 2] { self.lower }
    #[must_use]
    pub const fn upper(self) -> [f64; 2] { self.upper }
    #[must_use]
    pub const fn margin(self) -> f64 { self.margin }
}

/// Authored input and the exact prescribed values consumed by volume projection.
#[derive(Debug, Clone)]
pub struct PreparedDesignRegions {
    /// Geometry after imposing regions, BEFORE area projection or physics.
    pub geometry: GridSdf,
    /// Union with the supplied fixed nodes, sorted and unique.
    pub fixed_nodes: Vec<(usize, f64)>,
    /// Per-request cell bounds [i_begin, j_begin, i_end, j_end], end exclusive.
    pub covered_cells: Vec<[usize; 4]>,
    /// Unique region-controlled nodes by phase; existing fixed nodes may overlap.
    pub material_nodes: usize,
    pub void_nodes: usize,
    /// Number of phi values actually changed from the input bits.
    pub changed_nodes: usize,
}

/// Cooperative rasterization boundaries. At most one lattice row per callback.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DesignRegionStage {
    Prepare,
    Rasterize { region: usize, row: usize },
    StageRow(usize),
    Publish,
}

fn invalid(message: impl Into<String>) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: message.into() }
}

/// Author non-design geometry without modifying the input or performing a solve.
///
/// Same-phase overlaps take the strongest margin, independent of request order.
/// Opposite phases sharing ANY required node refuse, even if their rectangles
/// do not overlap: this resolution cannot express both prescriptions. Existing
/// fixed values must match the input bits and satisfy every requested margin;
/// they are never replaced to resolve a conflict. Other nodes keep their bits.
///
/// Supply BOTH returned geometry and fixed nodes to the projected optimizer.
/// Its existing area projection, checkpoints and refinement then retain the
/// regions. Area/physics feasibility must be established by that owner, not here.
///
/// # Errors
/// Refuses malformed input, conflicting prescriptions or more than 64 regions.
pub fn prepare_design_regions(
    input: &GridSdf, fixed: &[(usize, f64)], regions: &[DesignRegion],
) -> Result<PreparedDesignRegions, CutFemError> {
    match prepare_design_regions_controlled(input, fixed, regions, |_| {
        ControlFlow::<Infallible>::Continue(())
    })? {
        ControlFlow::Continue(prepared) => Ok(prepared),
        ControlFlow::Break(never) => match never {},
    }
}

/// Controlled counterpart. Refusal and every interruption leave both inputs
/// unchanged and return no partially authored geometry. No wall-time bound is
/// claimed for validation, allocation, cloning, or an individual lattice row.
///
/// # Errors
/// Returns the same admission and conflict errors as [`prepare_design_regions`].
pub fn prepare_design_regions_controlled<B>(
    input: &GridSdf, fixed: &[(usize, f64)], regions: &[DesignRegion],
    mut control: impl FnMut(DesignRegionStage) -> ControlFlow<B>,
) -> Result<ControlFlow<B, PreparedDesignRegions>, CutFemError> {
    let n = input.n();
    if !n.is_power_of_two() || !(2..=256).contains(&n) || regions.len() > 64
        || fixed.len() > input.nodes().len()
        || input.nodes().iter().any(|value| !value.is_finite())
    {
        return Err(invalid("design regions require a finite dyadic 2..=256 lattice, at most 64 regions and an in-range fixed-node count"));
    }
    let mut previous = None;
    for &(index, value) in fixed {
        if index >= input.nodes().len() || !value.is_finite()
            || previous.is_some_and(|old| index <= old)
            || input.nodes()[index].to_bits() != value.to_bits()
        {
            return Err(invalid("existing fixed nodes must be sorted, unique, finite and bitwise equal to their input values"));
        }
        previous = Some(index);
    }
    if let ControlFlow::Break(reason) = control(DesignRegionStage::Prepare) {
        return Ok(ControlFlow::Break(reason));
    }
    let mut required: Vec<Option<(DesignPhase, f64)>> = vec![None; input.nodes().len()];
    let mut covered_cells = Vec::with_capacity(regions.len());
    #[allow(clippy::cast_precision_loss)]
    let scale = n as f64;
    for (ordinal, region) in regions.iter().enumerate() {
        // Multiplication by a dyadic grid extent is exact for these [0,1]
        // coordinates. Cover cells intersecting the rectangle's interior.
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let bounds = [
            (region.lower[0] * scale).floor() as usize,
            (region.lower[1] * scale).floor() as usize,
            (region.upper[0] * scale).ceil() as usize,
            (region.upper[1] * scale).ceil() as usize,
        ];
        covered_cells.push(bounds);
        for j in bounds[1]..=bounds[3] {
            if let ControlFlow::Break(reason) = control(DesignRegionStage::Rasterize { region: ordinal, row: j }) {
                return Ok(ControlFlow::Break(reason));
            }
            for i in bounds[0]..=bounds[2] {
                let index = i + j * (n + 1);
                match &mut required[index] {
                    Some((phase, margin)) => {
                        if *phase != region.phase {
                            return Err(invalid(format!("material and void regions conflict at lattice node ({i},{j}); separate them or refine the input lattice")));
                        }
                        *margin = (*margin).max(region.margin);
                    }
                    slot @ None => *slot = Some((region.phase, region.margin)),
                }
            }
        }
    }
    let mut prepared = PreparedDesignRegions {
        geometry: input.clone(), fixed_nodes: Vec::new(), covered_cells,
        material_nodes: 0, void_nodes: 0, changed_nodes: 0,
    };
    let mut cursor = 0;
    for j in 0..=n {
        if let ControlFlow::Break(reason) = control(DesignRegionStage::StageRow(j)) {
            return Ok(ControlFlow::Break(reason));
        }
        for i in 0..=n {
            let index = i + j * (n + 1);
            let original = input.nodes()[index];
            let pinned = fixed.get(cursor).is_some_and(|&(node, _)| node == index);
            let mut value = original;
            if let Some((phase, margin)) = required[index] {
                let satisfies = match phase {
                    DesignPhase::Material => {
                        prepared.material_nodes += 1;
                        value = original.min(-margin);
                        original <= -margin
                    }
                    DesignPhase::Void => {
                        prepared.void_nodes += 1;
                        value = original.max(margin);
                        original >= margin
                    }
                };
                if pinned && !satisfies {
                    return Err(invalid(format!("design region conflicts with the existing fixed value at node ({i},{j})")));
                }
            }
            if pinned { cursor += 1; value = original; }
            if pinned || required[index].is_some() { prepared.fixed_nodes.push((index, value)); }
            if value.to_bits() != original.to_bits() { prepared.changed_nodes += 1; }
            prepared.geometry.nodes_mut()[index] = value;
        }
    }
    if let ControlFlow::Break(reason) = control(DesignRegionStage::Publish) {
        return Ok(ControlFlow::Break(reason));
    }
    Ok(ControlFlow::Continue(prepared))
}

#[cfg(test)]
mod tests;
