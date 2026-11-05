//! Preserve the signed P1 motion of one plate eigensolve for spatial acoustics.
//!
//! A zero area integral is not a silent finite-wavelength mode. These samples
//! retain its separated positive and negative lobes in the original modal
//! coordinates. The existing Rayleigh owners supply propagation, interference
//! and receivers; no acoustic impedance or additional damping is introduced.

use super::{
    AcousticRealizeError, CompactBody, PlateChart, PlateChartRadiation, PlateMassModel,
    PlateMesh, PlateProjection, assemble_chart_model, harvest_radiators_with_modes,
};
use crate::pcm_wav::baffled::{MAX_RAYLEIGH_INPUTS, SurfaceSample};
use fs_math::det;

/// Complete positive-area surface quadrature limit, admitted before allocation.
/// Shape storage is additionally bounded by the existing 128-input Rayleigh cap.
pub const MAX_PLATE_SURFACE_SAMPLES: usize = 16_384;

/// Mechanical modes and their signed normal-velocity surface from ONE solve.
/// Public access is read-only so mode ordering and coordinate scaling stay
/// together until the caller consumes both with [`Self::into_parts`].
#[derive(Clone, Debug)]
pub struct PlateSurfaceRadiators {
    pub(crate) bodies: Vec<CompactBody>,
    pub(crate) surface: Vec<SurfaceSample>,
    pub(crate) maximum_cell_diameter_m: f64,
}

impl PlateSurfaceRadiators {
    /// The same coefficients, force projections and original modal coordinates
    /// returned by [`super::certified_chart_radiators_with_mass`].
    #[must_use]
    pub fn bodies(&self) -> &[CompactBody] {
        &self.bodies
    }

    /// Positive-area points on the chart's physical xy plane, with z=0.
    /// `mode_shape[k] * bodies[k].modal_velocity_m_s()` is that mode's normal
    /// velocity contribution. Shapes are not divided by mass or signed area.
    #[must_use]
    pub fn surface(&self) -> &[SurfaceSample] {
        &self.surface
    }

    /// Largest physical edge of any quadrature subtriangle [m]. It decreases
    /// by two per refinement level, independently of the structural mesh.
    /// This is geometric resolution, not an acoustic error certificate.
    #[must_use]
    pub const fn maximum_cell_diameter_m(&self) -> f64 {
        self.maximum_cell_diameter_m
    }

    /// Transfer the coupled coordinate definitions without copying either.
    #[must_use]
    pub fn into_parts(self) -> (Vec<CompactBody>, Vec<SurfaceSample>) {
        (self.bodies, self.surface)
    }
}

/// Retain the actual P1 modal surface while deriving the ordinary plate bodies.
/// Assembly, eigensolve, retained-mode order, damping and force normalization
/// are shared with [`super::certified_chart_radiators_with_mass`]. No second
/// solve, mode renormalization or silent-mode removal is performed.
///
/// `subdivisions` is the number of uniform **quartering levels**: zero uses
/// three positive degree-two points on each original triangle; one uses 12,
/// two uses 48, and so on. Each level halves the largest subtriangle edge.
/// This refines only acoustic quadrature of the SAME P1 displacement field.
/// It does not refine the structural eigenvectors or change their supports.
///
/// Geometry is the original flat chart at z=0, with positive normal velocity
/// toward +z. Receivers must use this same infinite-baffle frame and declared
/// medium. The receiver owner admits its wavelength, clearance, time/filter
/// and work limits; these samples alone supply no acoustic convergence bound,
/// finite-baffle scattering, fluid loading or higher-order plate displacement.
///
/// # Errors
/// Invalid refinement, more than [`MAX_PLATE_SURFACE_SAMPLES`] points or 128
/// requested modes, unrepresentable geometry/shapes, or the original complete
/// chart, force-footprint and modal admission failure. No partial image returns.
pub fn certified_chart_surface_radiators(
    chart: &PlateChart,
    options: &PlateChartRadiation,
    mass: PlateMassModel,
    subdivisions: usize,
) -> Result<PlateSurfaceRadiators, AcousticRealizeError> {
    let (divisions, points_per_triangle, samples) =
        surface_work(chart.mesh.tris.len(), options.n_modes, subdivisions).map_err(invalid)?;
    let model = assemble_chart_model(chart, options, mass)?;
    let projection = PlateProjection::with_force_weights(&chart.mesh, &options.unit_force_weights)?;
    let weights = triangle_weights(divisions, points_per_triangle);
    let (mut surface, maximum_cell_diameter_m) =
        sample_surface(&chart.mesh, &weights, divisions, options.n_modes, samples).map_err(invalid)?;
    let bodies = harvest_radiators_with_modes(
        &model,
        &projection,
        options.eigenvalue_window,
        options.n_modes,
        options.damping_ratio,
        None,
        |phi| append_mode(&chart.mesh, &model.dof_map, phi, &weights, &mut surface).map_err(invalid),
    )?;
    Ok(PlateSurfaceRadiators { bodies, surface, maximum_cell_diameter_m })
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

// All capacity arithmetic precedes rule, sample or modal allocation. There is
// no pairwise sample matrix: storage is O(points * requested_modes).
fn surface_work(
    triangles: usize, modes: usize, levels: usize,
) -> Result<(usize, usize, usize), &'static str> {
    if triangles == 0 || !(1..=MAX_RAYLEIGH_INPUTS).contains(&modes) {
        return Err("plate surface needs triangles and 1..=128 requested modes");
    }
    let shift = u32::try_from(levels).map_err(|_| "plate surface refinement exceeds sample budget")?;
    let divisions = 1usize.checked_shl(shift)
        .ok_or("plate surface refinement exceeds sample budget")?;
    let points = divisions.checked_mul(divisions).and_then(|n| n.checked_mul(3))
        .ok_or("plate surface refinement exceeds sample budget")?;
    let total = triangles.checked_mul(points)
        .filter(|&n| n <= MAX_PLATE_SURFACE_SAMPLES)
        .ok_or("plate surface refinement exceeds sample budget")?;
    Ok((divisions, points, total))
}

// Direct triangular lattice: the same uniform quartering as recursive cells,
// without allocating any intermediate cell trees. The three-point rule is
// positive and exact for P1 values and their quadratic products on every cell.
fn triangle_weights(divisions: usize, points: usize) -> Vec<[f64; 3]> {
    let inverse = (divisions as f64).recip();
    let vertex = |i: usize, j: usize| {
        [(divisions - i - j) as f64 * inverse, i as f64 * inverse, j as f64 * inverse]
    };
    let rule = [[2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0], [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0]];
    let mut weights = Vec::with_capacity(points);
    let mut cell = |vertices: [[f64; 3]; 3]| {
        for q in rule {
            weights.push(core::array::from_fn(|axis|
                (0..3).map(|i| q[i] * vertices[i][axis]).sum()));
        }
    };
    for i in 0..divisions {
        for j in 0..divisions - i {
            cell([vertex(i, j), vertex(i + 1, j), vertex(i, j + 1)]);
            if i + j + 1 < divisions {
                cell([vertex(i + 1, j), vertex(i + 1, j + 1), vertex(i, j + 1)]);
            }
        }
    }
    weights
}

// Mesh geometry/supports have already passed the original plate assembler.
fn sample_surface(
    mesh: &PlateMesh, weights: &[[f64; 3]], divisions: usize,
    modes: usize, samples: usize,
) -> Result<(Vec<SurfaceSample>, f64), &'static str> {
    let mut surface = Vec::with_capacity(samples);
    let mut maximum_cell_diameter_m = 0.0_f64;
    for tri in &mesh.tris {
        let p = tri.map(|node| mesh.nodes[node]);
        let area = 0.5 * ((p[1].0 - p[0].0) * (p[2].1 - p[0].1)
            - (p[2].0 - p[0].0) * (p[1].1 - p[0].1));
        let area_m2 = area / weights.len() as f64;
        if !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err("plate surface quadrature needs positive representable physical area");
        }
        for (a, b) in [(0, 1), (1, 2), (2, 0)] {
            let edge = det::hypot(p[b].0 - p[a].0, p[b].1 - p[a].1) / divisions as f64;
            if !edge.is_finite() || edge <= 0.0 {
                return Err("plate surface cell diameter is not representable");
            }
            maximum_cell_diameter_m = maximum_cell_diameter_m.max(edge);
        }
        for barycentric in weights {
            let x = (0..3).map(|i| barycentric[i] * p[i].0).sum::<f64>();
            let y = (0..3).map(|i| barycentric[i] * p[i].1).sum::<f64>();
            if !x.is_finite() || !y.is_finite() {
                return Err("plate surface quadrature position is not representable");
            }
            surface.push(SurfaceSample { position_m: [x, y, 0.0], area_m2,
                mode_shape: Vec::with_capacity(modes) });
        }
    }
    Ok((surface, maximum_cell_diameter_m))
}

// Retain the ORIGINAL coordinate scaling. Fetch each triangle's three nodal
// displacements once per mode; slope DOFs affect mechanics but are not a P1
// surface displacement. Eliminated transverse DOFs remain exactly at rest.
fn append_mode(
    mesh: &PlateMesh, dof_map: &[Option<usize>], phi: &[f64],
    weights: &[[f64; 3]], surface: &mut [SurfaceSample],
) -> Result<(), &'static str> {
    for (tri, points) in mesh.tris.iter().zip(surface.chunks_exact_mut(weights.len())) {
        let nodal = tri.map(|node| dof_map[3 * node].map_or(0.0, |i| phi[i]));
        for (sample, barycentric) in points.iter_mut().zip(weights) {
            let shape = (0..3).map(|i| barycentric[i] * nodal[i]).sum::<f64>();
            if !shape.is_finite() { return Err("plate surface modal shape is not representable"); }
            sample.mode_shape.push(shape);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
