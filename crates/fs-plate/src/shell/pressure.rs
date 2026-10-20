//! One oriented P1 normal trace for pressure work and displaced volume.
//!
//! A selection addresses the original midsurface triangles. It does not extrude
//! skins, infer which material bounds a fluid, certify a closed cavity, or add
//! an acoustic radiation law. Unselected structural triangles still belong to
//! the mechanical model. Positive and negative identify winding orientation,
//! not the material volume enclosed by the existing two-sided skin extrusion.

use super::{ShellMesh, profile::ProfileBudget};
use crate::PlateError;
use std::collections::BTreeSet;

pub use super::reduction::radiation::ShellFace;

/// One explicitly oriented original midsurface facet in a pressure port.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OrientedShellFacet {
    /// Original global triangle index; no copied or renumbered mesh.
    pub triangle: usize,
    /// Normal along or against that triangle's original winding.
    pub side: ShellFace,
}

/// Immutable physical pressure load and its source facet geometry.
#[derive(Clone, Debug)]
pub struct ShellPressurePort {
    facets: Vec<OrientedShellFacet>,
    geometry: Vec<([usize; 3], [f64; 3])>,
    physical_row: Vec<f64>,
    area_m2: f64,
}

fn bad(what: &'static str) -> PlateError { PlateError::BadConstraint { what } }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PlateError> {
    if cancelled() { Err(PlateError::Cancelled) } else { Ok(()) }
}

/// The original shell volume-port integral, shared by DKT and MITC consumers.
/// Each area vector already includes the selected orientation and facet area.
pub(crate) fn integrate_normal_displacement(
    facets: impl IntoIterator<Item = ([usize; 3], [f64; 3])>,
    mut displacement: impl FnMut(usize) -> [f64; 3],
) -> f64 {
    facets.into_iter().map(|(tri, area)| {
        tri.into_iter().map(|node| {
            let u = displacement(node);
            (area[0] * u[0] + area[1] * u[1] + area[2] * u[2]) / 3.0
        }).sum::<f64>()
    }).sum()
}

/// Complete scalar work for validating and projecting mode-major global rows.
pub fn projection_terms(nodes: usize, facets: usize, modes: usize) -> Option<usize> {
    nodes.checked_mul(6)?.checked_add(facets.checked_mul(12)?)?.checked_mul(modes)
}

impl ShellPressurePort {
    /// Integrate unit pressure on an explicit facet selection. A repeated facet
    /// refuses even with the opposite sign; absence explicitly means no load.
    /// The supplied geometry budget also bounds the full original node indexing.
    /// `max_feature_evaluations` is irrelevant to this explicit-mesh operation.
    pub fn new(mesh: &ShellMesh, facets: Vec<OrientedShellFacet>, budget: ProfileBudget,
        mut cancelled: impl FnMut() -> bool) -> Result<Self, PlateError>
    {
        poll(&mut cancelled)?;
        if mesh.nodes.len() < 3 || mesh.nodes.len() > budget.max_nodes
            || mesh.tris.is_empty() || mesh.tris.len() > budget.max_triangles
            || facets.len() > mesh.tris.len()
            || mesh.nodes.iter().flatten().any(|v| !v.is_finite()) {
            return Err(bad("shell pressure port exceeds its original geometry allowance"));
        }
        let width = mesh.nodes.len().checked_mul(6)
            .ok_or_else(|| bad("shell pressure row extent overflow"))?;
        let mut physical_row = vec![0.0; width];
        let mut geometry = Vec::with_capacity(facets.len());
        let mut seen = BTreeSet::new();
        let mut area_m2 = 0.0;
        for selection in &facets {
            poll(&mut cancelled)?;
            if !seen.insert(selection.triangle) {
                return Err(bad("shell pressure surface repeats an original facet"));
            }
            let facet = mesh.facet(selection.triangle)?;
            let triangle = mesh.tris[selection.triangle];
            let sign = match selection.side { ShellFace::Positive => 1.0, ShellFace::Negative => -1.0 };
            let area = facet.frame[2].map(|v| sign * facet.area_m2 * v);
            area_m2 += facet.area_m2;
            for &node in &triangle {
                for c in 0..3 { physical_row[6 * node + c] += area[c] / 3.0; }
            }
            geometry.push((triangle, area));
        }
        if !area_m2.is_finite() || physical_row.iter().any(|v| !v.is_finite()) {
            return Err(bad("shell pressure area or physical load overflows"));
        }
        Ok(Self { facets, geometry, physical_row, area_m2 })
    }

    /// Original facet identities and authored orientations.
    pub fn facets(&self) -> &[OrientedShellFacet] { &self.facets }
    /// Actual selected midsurface area [m²], independent of modal cancellation.
    /// An explicitly empty surface has zero area.
    pub const fn area_m2(&self) -> f64 { self.area_m2 }
    /// Unit-pressure load in original GLOBAL (u,v,w,rx,ry,rz) node order [m²].
    /// Dotting this same row with physical velocity gives signed volume flow.
    /// Pure midsurface normal traction supplies no thickness-offset moment.
    pub fn physical_row_m2(&self) -> &[f64] { &self.physical_row }

    /// Project the original oriented normal trace on global physical modes.
    /// Zero participation is retained: mechanical modes need not pump this fluid.
    pub fn project_modes(&self, modes: &[Vec<f64>], max_projection_terms: usize,
        mut cancelled: impl FnMut() -> bool) -> Result<Vec<f64>, PlateError>
    {
        poll(&mut cancelled)?;
        if projection_terms(self.physical_row.len() / 6, self.facets.len(), modes.len())
            .is_none_or(|n| n > max_projection_terms) {
            return Err(bad("shell pressure projection exceeds its declared work allowance"));
        }
        let mut shapes = Vec::with_capacity(modes.len());
        for mode in modes {
            poll(&mut cancelled)?;
            if mode.len() != self.physical_row.len() || mode.iter().any(|v| !v.is_finite()) {
                return Err(bad("shell pressure projection needs finite original physical modes"));
            }
            let value = integrate_normal_displacement(self.geometry.iter().copied(),
                |node| [mode[6 * node], mode[6 * node + 1], mode[6 * node + 2]]);
            if !value.is_finite() { return Err(bad("shell pressure projection overflows")); }
            shapes.push(value);
        }
        Ok(shapes)
    }
}

#[cfg(test)]
mod tests;
