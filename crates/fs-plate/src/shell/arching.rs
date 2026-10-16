//! Explicit graph surfaces and independent graduation fields (IRP M1.2-W6).
//! Tensor cubic Hermite interpolation is C1 and reproduces bicubic data.
//! Nodal derivatives on a sampled grid are obtained from four-point Lagrange
//! polynomials; they are interpolation choices, not measured slope data.
//! A lift preserves topology and reference load state; it never unloads a scan.

use super::{ShellMesh, profile::ProfileBudget};
use crate::{PlateError, PlateMesh};
use fs_math::det;

fn bad(what: &'static str) -> PlateError {
    PlateError::BadSection { what }
}

/// The reference geometry's declared mechanical state, retained by every lift.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReferenceLoadState {
    /// Source explicitly identifies a stress-free geometry.
    Unloaded,
    /// Source identifies a loaded geometry. No stress-free shape is inferred.
    Loaded { description: String },
}

/// A section graph (coordinate [m], height [m], derivative [m/m]).
#[derive(Debug, Clone)]
pub struct Curve2 {
    coordinates_m: Vec<f64>,
    heights_m: Vec<f64>,
    slopes: Vec<f64>,
}

fn axis_valid(axis: &[f64], minimum: usize) -> bool {
    axis.len() >= minimum
        && axis.iter().all(|x| x.is_finite())
        && axis.windows(2).all(|p| p[1] > p[0])
}

fn cell(axis: &[f64], x: f64) -> Result<usize, PlateError> {
    if !x.is_finite() || x < axis[0] || x > axis[axis.len() - 1] {
        return Err(bad(
            "surface query is outside the declared finite parameter domain",
        ));
    }
    Ok(axis
        .partition_point(|a| *a <= x)
        .saturating_sub(1)
        .min(axis.len() - 2))
}

fn hermite(x0: f64, x1: f64, v0: f64, v1: f64, d0: f64, d1: f64, x: f64) -> (f64, f64) {
    let h = x1 - x0;
    let t = (x - x0) / h;
    let t2 = t * t;
    let t3 = t2 * t;
    let value = (2.0 * t3 - 3.0 * t2 + 1.0) * v0
        + (-2.0 * t3 + 3.0 * t2) * v1
        + (t3 - 2.0 * t2 + t) * h * d0
        + (t3 - t2) * h * d1;
    let derivative = ((6.0 * t2 - 6.0 * t) * v0 + (-6.0 * t2 + 6.0 * t) * v1) / h
        + (3.0 * t2 - 4.0 * t + 1.0) * d0
        + (3.0 * t2 - 2.0 * t) * d1;
    (value, derivative)
}

impl Curve2 {
    /// No slope or endpoint condition is defaulted.
    pub fn new(
        coordinates_m: Vec<f64>,
        heights_m: Vec<f64>,
        slopes: Vec<f64>,
    ) -> Result<Self, PlateError> {
        if !axis_valid(&coordinates_m, 2)
            || heights_m.len() != coordinates_m.len()
            || slopes.len() != coordinates_m.len()
            || heights_m.iter().chain(&slopes).any(|x| !x.is_finite())
        {
            return Err(bad(
                "section requires increasing coordinates, heights and explicit finite slopes",
            ));
        }
        Ok(Self {
            coordinates_m,
            heights_m,
            slopes,
        })
    }
    /// Height and analytic derivative of the declared Hermite section.
    pub fn sample(&self, x: f64) -> Result<(f64, f64), PlateError> {
        let i = cell(&self.coordinates_m, x)?;
        let result = hermite(
            self.coordinates_m[i],
            self.coordinates_m[i + 1],
            self.heights_m[i],
            self.heights_m[i + 1],
            self.slopes[i],
            self.slopes[i + 1],
            x,
        );
        if !result.0.is_finite() || !result.1.is_finite() {
            return Err(bad("section interpolation is unrepresentable"));
        }
        Ok(result)
    }
}

// Derivative at a knot of the local interpolating cubic. Four explicit
// samples, including at endpoints, make bicubic reproduction load-bearing.
fn knot_slope(axis: &[f64], values: &[f64], at: usize) -> f64 {
    let start = at.saturating_sub(1).min(axis.len() - 4);
    let x = axis[at];
    let mut sum = 0.0;
    for j in start..start + 4 {
        let mut derivative = 0.0;
        for k in start..start + 4 {
            if k == j {
                continue;
            }
            let mut term = 1.0 / (axis[j] - axis[k]);
            for m in start..start + 4 {
                if m != j && m != k {
                    term *= (x - axis[m]) / (axis[j] - axis[m]);
                }
            }
            derivative += term;
        }
        sum += derivative * values[j];
    }
    sum
}
fn sampled_curve(axis: &[f64], values: &[f64], x: f64) -> Result<(f64, f64), PlateError> {
    let i = cell(axis, x)?;
    Ok(hermite(
        axis[i],
        axis[i + 1],
        values[i],
        values[i + 1],
        knot_slope(axis, values, i),
        knot_slope(axis, values, i + 1),
        x,
    ))
}

/// A caller-sampled rectangular grid, row-major (y then x), all lengths in m.
#[derive(Debug, Clone)]
pub struct BicubicGrid {
    x_m: Vec<f64>,
    y_m: Vec<f64>,
    values_m: Vec<f64>,
}
impl BicubicGrid {
    /// At least four samples per axis; missing or nonfinite values refuse.
    pub fn new(x_m: Vec<f64>, y_m: Vec<f64>, values_m: Vec<f64>) -> Result<Self, PlateError> {
        if !axis_valid(&x_m, 4)
            || !axis_valid(&y_m, 4)
            || x_m.len().checked_mul(y_m.len()) != Some(values_m.len())
            || values_m.iter().any(|v| !v.is_finite())
        {
            return Err(bad(
                "bicubic grid needs four increasing samples per axis and every finite height",
            ));
        }
        Ok(Self { x_m, y_m, values_m })
    }
    /// (height, [dz/dx,dz/dy]); no extrapolation outside the supplied grid.
    pub fn sample(&self, p: [f64; 2]) -> Result<(f64, [f64; 2]), PlateError> {
        cell(&self.x_m, p[0])?;
        cell(&self.y_m, p[1])?;
        let mut heights = Vec::with_capacity(self.y_m.len());
        let mut slopes = Vec::with_capacity(self.y_m.len());
        for row in self.values_m.chunks_exact(self.x_m.len()) {
            let (z, dx) = sampled_curve(&self.x_m, row, p[0])?;
            heights.push(z);
            slopes.push(dx);
        }
        let (z, dy) = sampled_curve(&self.y_m, &heights, p[1])?;
        let (dx, _) = sampled_curve(&self.y_m, &slopes, p[1])?;
        if ![z, dx, dy].iter().all(|v| v.is_finite()) {
            return Err(bad("surface interpolation is unrepresentable"));
        }
        Ok((z, [dx, dy]))
    }
}

/// Geometry parameters only. Thickness is supplied through a separate field.
#[derive(Debug, Clone)]
pub enum ArchingParameters {
    /// Long coordinate is y, cross coordinate x. Every cross section includes
    /// x=0 and agrees there with the long section; four y stations are needed.
    LoftedSections {
        long: Curve2,
        cross: Vec<(f64, Curve2)>,
    },
    /// Sampled height, with explicit bicubic interpolation semantics above.
    HeightMap(BicubicGrid),
}
impl ArchingParameters {
    fn sections_per_query(&self) -> usize {
        match self {
            Self::HeightMap(grid) => grid.y_m.len(),
            Self::LoftedSections { cross, .. } => cross.len().saturating_add(1),
        }
    }
    fn validate(&self) -> Result<(), PlateError> {
        if let Self::LoftedSections { long, cross } = self {
            let stations: Vec<_> = cross.iter().map(|(s, _)| *s).collect();
            if !axis_valid(&stations, 4) {
                return Err(bad("loft requires four increasing cross-section stations"));
            }
            for (station, curve) in cross {
                let z = long.sample(*station)?.0;
                let middle = curve.sample(0.0)?.0;
                if (z - middle).abs() > 1e-12 * z.abs().max(middle.abs()).max(1.0) {
                    return Err(bad(
                        "cross section and long arch disagree at their intersection",
                    ));
                }
            }
        }
        Ok(())
    }
    fn sample(&self, p: [f64; 2]) -> Result<(f64, [f64; 2]), PlateError> {
        match self {
            Self::HeightMap(grid) => grid.sample(p),
            Self::LoftedSections { long, cross } => {
                let stations: Vec<_> = cross.iter().map(|(s, _)| *s).collect();
                let mut z = Vec::new();
                let mut dx = Vec::new();
                let mut centre = Vec::new();
                for (_, curve) in cross {
                    let sample = curve.sample(p[0])?;
                    z.push(sample.0);
                    dx.push(sample.1);
                    centre.push(curve.sample(0.0)?.0);
                }
                let (surface, dy) = sampled_curve(&stations, &z, p[1])?;
                let (slope_x, _) = sampled_curve(&stations, &dx, p[1])?;
                let (mid, mid_dy) = sampled_curve(&stations, &centre, p[1])?;
                let (long_z, long_dy) = long.sample(p[1])?;
                let result = (surface + long_z - mid, [slope_x, dy + long_dy - mid_dy]);
                if !result.0.is_finite() || result.1.iter().any(|x| !x.is_finite()) {
                    return Err(bad("loft interpolation is unrepresentable"));
                }
                Ok(result)
            }
        }
    }
}

/// Admitted geometric chart with explicit reference load state.
#[derive(Debug, Clone)]
pub struct ArchingSurface {
    parameters: ArchingParameters,
    reference_state: ReferenceLoadState,
}

/// Same connectivity and vertex ordering, analytic unit graph normals and the
/// original load-state declaration. Global self-intersection is not certified.
#[derive(Debug, Clone)]
pub struct LiftedSurface {
    pub mesh: ShellMesh,
    pub directors: Vec<[f64; 3]>,
    pub reference_state: ReferenceLoadState,
}
impl ArchingSurface {
    pub fn new(
        parameters: ArchingParameters,
        reference_state: ReferenceLoadState,
    ) -> Result<Self, PlateError> {
        parameters.validate()?;
        if matches!(&reference_state,ReferenceLoadState::Loaded{description} if description.trim().is_empty())
        {
            return Err(bad(
                "loaded reference geometry requires its source load-state description",
            ));
        }
        Ok(Self {
            parameters,
            reference_state,
        })
    }
    pub fn sample(&self, p: [f64; 2]) -> Result<(f64, [f64; 2]), PlateError> {
        self.parameters.sample(p)
    }
    /// Poll before allocation and per node, and retain every hole boundary.
    pub fn lift(
        &self,
        mesh: &PlateMesh,
        budget: ProfileBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LiftedSurface, PlateError> {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        let work = mesh
            .nodes
            .len()
            .checked_mul(self.parameters.sections_per_query())
            .ok_or(bad("arching evaluation count overflow"))?;
        if mesh.nodes.len() > budget.max_nodes
            || mesh.tris.len() > budget.max_triangles
            || work > budget.max_feature_evaluations
        {
            return Err(bad("arching lift exceeds explicit mesh/evaluation budgets"));
        }
        let mut nodes = Vec::with_capacity(mesh.nodes.len());
        let mut directors = Vec::with_capacity(mesh.nodes.len());
        for &(x, y) in &mesh.nodes {
            if cancelled() {
                return Err(PlateError::Cancelled);
            }
            let (z, [dx, dy]) = self.sample([x, y])?;
            let norm = det::sqrt(1.0 + dx * dx + dy * dy);
            if !z.is_finite() || !norm.is_finite() {
                return Err(bad("arching height or director is nonfinite"));
            }
            nodes.push([x, y, z]);
            directors.push([-dx / norm, -dy / norm, 1.0 / norm]);
        }
        Ok(LiftedSurface {
            mesh: ShellMesh::new(nodes, mesh.tris.clone())?,
            directors,
            reference_state: self.reference_state.clone(),
        })
    }
}

/// Independent graduation [m], with the same two interpolation choices.
#[derive(Debug, Clone)]
pub struct ThicknessField {
    parameters: ArchingParameters,
}
impl ThicknessField {
    pub fn new(parameters: ArchingParameters) -> Result<Self, PlateError> {
        parameters.validate()?;
        Ok(Self { parameters })
    }
    /// Refuse nonpositive interpolated thickness rather than clamp it.
    pub fn sample(&self, p: [f64; 2]) -> Result<f64, PlateError> {
        let h = self.parameters.sample(p)?.0;
        if !h.is_finite() || h <= 0.0 {
            return Err(bad("graduation must be positive and finite at every query"));
        }
        Ok(h)
    }
    pub fn sample_mesh(
        &self,
        mesh: &PlateMesh,
        maximum_evaluations: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, PlateError> {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        if mesh.nodes.len() > maximum_evaluations {
            return Err(bad("graduation exceeds explicit evaluation budget"));
        }
        mesh.nodes
            .iter()
            .map(|&(x, y)| {
                if cancelled() {
                    Err(PlateError::Cancelled)
                } else {
                    self.sample([x, y])
                }
            })
            .collect()
    }
}

#[cfg(test)]
mod tests;
