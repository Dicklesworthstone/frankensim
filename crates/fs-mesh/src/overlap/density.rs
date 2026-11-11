//! Density-weighted moments on the existing tetrahedron/box intersections.

use super::{
    OverlapError, Point, clipped_tetrahedron, determinant, fan_det_noise, poll, position, sub,
};
use fs_exec::Cx;

mod moments;
#[cfg(test)]
mod tests;

/// Geometric volume and exact polynomial mass integrals of an affine density.
///
/// All coordinates and density samples use one consistent length/mass system.
/// The first and second moments are about the explicit caller-supplied origin.
/// Together they integrate affine momentum/enthalpy and quadratic kinetic or
/// angular moments without replacing the original density by a cell average.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DensityMoments {
    /// Unweighted geometric volume [length³]; zero for an empty intersection.
    pub volume: f64,
    /// Integral of the positive affine density [mass].
    pub mass: f64,
    /// Integral of rho * (x - moment_origin) [mass length].
    pub first: [f64; 3],
    /// Integral of rho * (x - moment_origin)(x - moment_origin)^T [mass length²].
    /// Symmetric entries share the same computed value.
    pub second: [[f64; 3]; 3],
}

/// Intersect an original tetrahedron and box, retaining its P1 density field.
///
/// Each density sample is positive mass per original geometric volume at the
/// corresponding input vertex. The same six-plane clipping used by
/// [super::tetrahedron_box_moments] carries that affine density through every
/// edge intersection. Positive simplex integration then includes the cubic
/// density/position product required for inertia and kinetic energy. No
/// centroid sample, density requery, or uniform-density approximation occurs.
///
/// Supplying transported rho_reference/J samples admits a deformed affine cell;
/// the caller owns that reference/current convention. This function integrates
/// the supplied state and does not infer thermodynamic expansion, material
/// mixing, phase change, or a conservative history update.
///
/// Bounds and cancellation are those of the original intersection owner. Work
/// and storage stay bounded by the same clipped faces; the field adds one
/// scalar per vertex, with no extra clipping or quadrature refinement. Moments
/// are exact polynomial integrals on measured f64 geometry, not certified
/// enclosures. Callers must audit partition coverage and transfer defects.
///
/// # Errors
/// Refuses invalid/nonpositive densities, a density ratio not representable
/// during normalized interpolation, unrepresentable positive mass/moments,
/// invalid geometry, the existing geometry limits, or cancellation. Empty
/// geometric intersections return zero volume and mass.
pub fn tetrahedron_box_density_moments(
    vertices: [[f64; 3]; 4],
    vertex_density: [f64; 4],
    lower: [f64; 3],
    upper: [f64; 3],
    moment_origin: [f64; 3],
    cx: &Cx<'_>,
) -> Result<DensityMoments, OverlapError> {
    poll(cx)?;
    if vertex_density
        .iter()
        .any(|rho| !rho.is_finite() || *rho <= 0.0)
    {
        return Err(OverlapError::InvalidInput);
    }
    let density_scale = vertex_density.iter().copied().fold(0.0_f64, f64::max);
    let relative = vertex_density.map(|rho| rho / density_scale);
    if relative.iter().any(|rho| *rho <= 0.0) {
        return Err(OverlapError::Unrepresentable);
    }
    let attributed =
        core::array::from_fn(|a| [vertices[a][0], vertices[a][1], vertices[a][2], relative[a]]);
    let Some(clipped) = clipped_tetrahedron(attributed, lower, upper, moment_origin, cx)? else {
        return Ok(DensityMoments::default());
    };
    let count = clipped.faces.iter().map(Vec::len).sum::<usize>();
    let center: [f64; 4] = core::array::from_fn(|i| {
        clipped.faces.iter().flatten().map(|p| p[i]).sum::<f64>() / count as f64
    });
    let center_position = position(center);
    let offset = sub(clipped.base, moment_origin);
    if offset.iter().any(|x| !x.is_finite()) {
        return Err(OverlapError::Unrepresentable);
    }
    let mut sums = [0.0; 14];
    let mut corrections = [0.0; 14];
    for face in &clipped.faces {
        poll(cx)?;
        for i in 1..face.len() - 1 {
            let vertices = [center, face[0], face[i], face[i + 1]];
            let points = vertices.map(position);
            let edges = [
                sub(points[1], center_position),
                sub(points[2], center_position),
                sub(points[3], center_position),
            ];
            let det = determinant(edges[0], edges[1], edges[2]);
            if det < 0.0 {
                // The error model of the fan determinant is in overlap.rs (fan_det_noise).
                if det >= -fan_det_noise(edges) {
                    continue;
                }
                return Err(OverlapError::Unrepresentable);
            }
            if det == 0.0 {
                continue;
            }
            let volume = ((det / 6.0 * clipped.scale) * clipped.scale) * clipped.scale;
            let relative_positions: [Point; 4] =
                points.map(|p| core::array::from_fn(|j| offset[j] + clipped.scale * p[j]));
            let values = moments::integrate(
                volume,
                relative_positions,
                vertices.map(|p| p[3]),
                density_scale,
            )
            .ok_or(OverlapError::Unrepresentable)?;
            for j in 0..14 {
                let adjusted = values[j] - corrections[j];
                let next = sums[j] + adjusted;
                corrections[j] = (next - sums[j]) - adjusted;
                sums[j] = next;
            }
        }
    }
    poll(cx)?;
    if sums.iter().any(|x| !x.is_finite())
        || sums[0] > 0.0 && (sums[1] <= 0.0 || (0..3).any(|j| sums[5 + 4 * j] <= 0.0))
    {
        return Err(OverlapError::Unrepresentable);
    }
    Ok(DensityMoments {
        volume: sums[0],
        mass: sums[1],
        first: [sums[2], sums[3], sums[4]],
        second: core::array::from_fn(|j| core::array::from_fn(|k| sums[5 + 3 * j + k])),
    })
}
