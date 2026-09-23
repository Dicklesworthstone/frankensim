//! Contact kinematics on the same skin used by shell radiation.
//!
//! A material point has displacement u + theta x d, where d is its signed
//! reference director offset. Transposing this SAME row maps a point force to
//! both translation and moment. Geometry remains the undeformed reference;
//! this is neither finite-rotation kinematics nor collision detection.
use super::{PlateError, ShellReduction, bad, dot};

/// Which side of the oriented midsurface carries the material point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShellFace {
    /// Along the area-weighted oriented vertex director.
    Positive,
    /// Opposite the area-weighted oriented vertex director.
    Negative,
}
impl ShellFace {
    fn sign(self) -> f64 { match self { Self::Positive => 1.0, Self::Negative => -1.0 } }
}

/// A finite-thickness reference point and its conjugate force/velocity row.
#[derive(Debug, Clone)]
pub struct ShellSurfacePort {
    /// Position in the original shell frame [m].
    pub position_m: [f64; 3],
    /// Unnormalized modal velocity participation [1/sqrt(kg)]. Dot with modal
    /// velocity for point velocity; multiply by force [N] for generalized force.
    pub weights: Vec<f64>,
}

// Shared by surface contacts and radiation: thickness and director geometry
// must never acquire two independently maintained definitions.
pub(super) fn directors(s: &ShellReduction, h: &[f64]) -> Result<Vec<[f64; 3]>, PlateError> {
    if h.len() != s.nodes || h.iter().any(|v| !v.is_finite() || *v <= 0.0) {
        return Err(bad("shell skin requires positive finite thickness at every source node"));
    }
    let mut out = vec![[0.0; 3]; s.nodes];
    for (f, tri) in s.triangles.iter().enumerate() {
        let actual = tri.iter().map(|&i| h[i] / 3.0).sum::<f64>();
        let expected = s.section_thicknesses[f];
        if !actual.is_finite() || (actual - expected).abs() > 1e-10 * expected {
            return Err(bad("shell skin thickness does not match the mechanical sections"));
        }
        for &i in tri { for c in 0..3 { out[i][c] += s.area_normals[f][c]; } }
    }
    for (i, d) in out.iter_mut().enumerate() {
        let norm = dot(*d, *d).sqrt();
        if !norm.is_finite() || norm <= 0.0 { return Err(bad("shell skin has an undefined vertex director")); }
        for v in d { *v = (*v / norm) * (0.5 * h[i]); }
    }
    Ok(out)
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}
impl ShellReduction {
    /// Project a physical point on the thickness-offset skin, retaining moments.
    /// Barycentrics use the ORIGINAL midsurface triangle's vertex order, on
    /// either side; the negative acoustic triangle has the reverse winding.
    /// Direction is a unit vector in the original shell frame. Thickness must
    /// match every source section, exactly as for `radiation_surface`.
    ///
    /// This cold query allocates O(nodes + modes) workspace. It does not locate
    /// contacts, move a point, renormalize shapes, or alter any mechanical state.
    ///
    /// # Errors
    /// Invalid facet/barycentrics/direction, incompatible thickness, undefined
    /// directors, or a nonfinite derived position/participation.
    pub fn surface_point_port(&self, nodal_thickness_m: &[f64], triangle: usize,
        barycentric: [f64; 3], face: ShellFace, direction: [f64; 3])
        -> Result<ShellSurfacePort, PlateError>
    {
        // Reuse the existing point-port admission (and its exact base row).
        let mut weights = self.point_port(triangle, barycentric, direction)?;
        let offsets = directors(self, nodal_thickness_m)?;
        let tri = self.triangles[triangle];
        let mut position_m = [0.0; 3];
        for a in 0..3 {
            let node = tri[a]; let offset = offsets[node].map(|v| v * face.sign());
            for c in 0..3 { position_m[c] += barycentric[a] * (self.reference_positions[node][c] + offset[c]); }
            for (mode, w) in weights.iter_mut().enumerate() {
                *w += barycentric[a] * dot(direction, cross(self.rotations[mode*self.nodes+node], offset));
            }
        }
        if position_m.iter().chain(weights.iter()).any(|v| !v.is_finite()) {
            return Err(bad("shell surface-point projection overflows"));
        }
        Ok(ShellSurfacePort { position_m, weights })
    }
}
