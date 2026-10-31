//! Source-owned P1 motion on the shell's original zero-thickness midsurface.
//!
//! This chart retains physical translations at each facet corner. Adjacent
//! folded facets have their own normals; they do not share a scalar velocity
//! obtained by averaging vertex normals. Rotations have zero offset here.
//! Structural thickness, inertia and stiffness remain in the source reduction.
use super::{PlateError, ShellReduction, bad, dot, reserve};
use std::collections::BTreeMap;

/// Final acoustic geometry/projection limits, independent of source mode count
/// and of the pressure-jump solver's own quadrature and dense-work limits.
#[derive(Clone, Copy, Debug)]
pub struct OpenScreenRadiationBudget {
    pub maximum_vertices: usize,
    pub maximum_facets: usize,
    /// Counts all three corners of every facet in every original source mode.
    pub maximum_corner_weights: usize,
}

/// An immutable acoustic subdivision of the original reference midsurface.
/// Rows retain every original mass-normalized source mode, in source order.
#[derive(Clone, Debug)]
pub struct ShellOpenScreenRadiation {
    vertices: Vec<[f64; 3]>,
    triangles: Vec<[usize; 3]>,
    weights: Vec<Vec<[f64; 3]>>,
    refinements: u32,
}

fn poll(cancelled: &mut dyn FnMut() -> bool) -> Result<(), PlateError> {
    if cancelled() {
        Err(PlateError::Cancelled)
    } else {
        Ok(())
    }
}

impl ShellReduction {
    /// Project THIS shell's own translations onto each original oriented facet.
    /// The result is a zero-thickness, two-sided open-screen source, not the
    /// extruded finite-thickness skin returned by `radiation_surface`.
    ///
    /// Each refinement splits every triangle into four conforming midpoint
    /// triangles, sharing new vertices on shared edges. The facet-local P1
    /// source is interpolated exactly, with no mode truncation, smoothing,
    /// area rescaling or prescribed zero velocity at acoustic free edges.
    /// This increases acoustic resolution while preserving the original
    /// piecewise-flat structural geometry and displacement space. It is NOT
    /// edge grading or structural/mode convergence. The downstream screen
    /// owner must still admit open manifold topology and its frequency band.
    ///
    /// # Errors
    /// Cancellation, overflow, allocation refusal or exceeded final geometry /
    /// complete projection limits. No partially refined chart is published.
    pub fn open_screen_radiation(
        &self,
        acoustic_refinements: u32,
        budget: OpenScreenRadiationBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ShellOpenScreenRadiation, PlateError> {
        poll(&mut cancelled)?;
        let modes = self.mode_count();
        let mut final_facets = self.triangles.len();
        for _ in 0..acoustic_refinements {
            poll(&mut cancelled)?;
            final_facets = final_facets
                .checked_mul(4)
                .ok_or_else(|| bad("open-screen refinement extent overflows"))?;
            if final_facets > budget.maximum_facets {
                return Err(bad("open-screen refinement exceeds its final facet limit"));
            }
        }
        if modes == 0
            || final_facets == 0
            || final_facets > budget.maximum_facets
            || self.nodes > budget.maximum_vertices
            || final_facets
                .checked_mul(modes)
                .and_then(|n| n.checked_mul(3))
                .is_none_or(|n| n > budget.maximum_corner_weights)
        {
            return Err(bad(
                "complete open-screen source exceeds its geometry/projection limits",
            ));
        }
        let mut vertices = reserve(self.nodes)?;
        vertices.extend_from_slice(&self.reference_positions);
        let mut triangles = reserve(self.triangles.len())?;
        triangles.extend_from_slice(&self.triangles);
        let mut weights = reserve(modes)?;
        for mode in 0..modes {
            poll(&mut cancelled)?;
            let mut row = reserve(self.triangles.len())?;
            for (facet, tri) in self.triangles.iter().enumerate() {
                if facet % 128 == 0 {
                    poll(&mut cancelled)?;
                }
                let area_normal = self.area_normals[facet];
                let length = dot(area_normal, area_normal).sqrt();
                if !length.is_finite() || length <= 0.0 {
                    return Err(bad("open-screen source normal is not representable"));
                }
                let normal = area_normal.map(|v| v / length);
                let corners =
                    tri.map(|node| dot(normal, self.translations[mode * self.nodes + node]));
                if corners.iter().any(|v| !v.is_finite()) {
                    return Err(bad("open-screen source projection overflows"));
                }
                row.push(corners);
            }
            weights.push(row);
        }
        for _ in 0..acoustic_refinements {
            poll(&mut cancelled)?;
            let mut edges = BTreeMap::<(usize, usize), usize>::new();
            let mut next = reserve(4 * triangles.len())?;
            for (facet, &[a, b, c]) in triangles.iter().enumerate() {
                if facet % 128 == 0 {
                    poll(&mut cancelled)?;
                }
                let mut midpoint = |i: usize, j: usize| -> Result<usize, PlateError> {
                    let key = (i.min(j), i.max(j));
                    if let Some(&node) = edges.get(&key) {
                        return Ok(node);
                    }
                    if vertices.len() >= budget.maximum_vertices {
                        return Err(bad("open-screen subdivision exceeds its vertex limit"));
                    }
                    let point = core::array::from_fn(|axis| {
                        0.5 * vertices[i][axis] + 0.5 * vertices[j][axis]
                    });
                    if point.iter().any(|v| !v.is_finite()) {
                        return Err(bad("open-screen midpoint is not representable"));
                    }
                    vertices
                        .try_reserve(1)
                        .map_err(|_| bad("open-screen vertex allocation refused"))?;
                    let node = vertices.len();
                    vertices.push(point);
                    edges.insert(key, node);
                    Ok(node)
                };
                let ab = midpoint(a, b)?;
                let bc = midpoint(b, c)?;
                let ca = midpoint(c, a)?;
                next.extend([[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]);
            }
            for row in &mut weights {
                poll(&mut cancelled)?;
                let mut refined = reserve(next.len())?;
                for (facet, &[a, b, c]) in row.iter().enumerate() {
                    if facet % 128 == 0 {
                        poll(&mut cancelled)?;
                    }
                    let (ab, bc, ca) = (0.5 * a + 0.5 * b, 0.5 * b + 0.5 * c, 0.5 * c + 0.5 * a);
                    refined.extend([[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]);
                }
                *row = refined;
            }
            triangles = next;
        }
        poll(&mut cancelled)?;
        Ok(ShellOpenScreenRadiation {
            vertices,
            triangles,
            weights,
            refinements: acoustic_refinements,
        })
    }
}

impl ShellOpenScreenRadiation {
    pub fn vertices_m(&self) -> &[[f64; 3]] {
        &self.vertices
    }
    pub fn triangles(&self) -> &[[usize; 3]] {
        &self.triangles
    }
    /// Mode-major facet-corner normal velocity [1/sqrt(kg)].
    pub fn normal_velocity_weights(&self) -> &[Vec<[f64; 3]>] {
        &self.weights
    }
    /// Uniform acoustic midpoint subdivisions, not an edge-grading receipt.
    pub fn acoustic_refinements(&self) -> u32 {
        self.refinements
    }
    /// Transfer the complete chart without copying its geometry or source rows.
    pub fn into_parts(self) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, Vec<Vec<[f64; 3]>>) {
        (self.vertices, self.triangles, self.weights)
    }
}

#[cfg(test)]
mod tests;
