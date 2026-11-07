//! Bounded linear density with the original directional face-flux owner.

use super::{
    DensityTransportCell, DensityTransportError, EdgeVolumeTransport, Result, TriMesh, admit_cell,
    allocate, poll, positive, positive_value,
};

#[derive(Clone, Copy, Default)]
pub(super) struct LinearCell {
    // Two vertex deviations from the mean, as linear functions of the three
    // neighboring density differences. The third deviation is minus their sum.
    variation: [[f64; 3]; 2],
}

pub(super) struct LinearReconstruction {
    cells: Vec<LinearCell>,
    // Volume-flux [centroid, variance] along each edge, for each sign interval.
    fractions: Vec<[[f64; 2]; 2]>,
}

fn invalid(message: &'static str) -> DensityTransportError {
    DensityTransportError::Invalid(message)
}

fn vertices(mesh: &TriMesh, cell: usize) -> Result<[[f64; 2]; 3]> {
    let mut points = [[0.0; 2]; 3];
    for (target, &vertex) in points.iter_mut().zip(&mesh.tris[cell]) {
        *target = *mesh
            .verts
            .get(vertex)
            .ok_or(invalid("original density vertex required"))?;
        if target.iter().any(|value| !value.is_finite()) {
            return Err(invalid("finite density reconstruction geometry required"));
        }
    }
    Ok(points)
}

fn centroid(points: [[f64; 2]; 3]) -> [f64; 2] {
    core::array::from_fn(|k| {
        points[0][k] + (points[1][k] - points[0][k]) / 3.0 + (points[2][k] - points[0][k]) / 3.0
    })
}

pub(super) fn neighbor(mesh: &TriMesh, cell: usize, edge: usize) -> usize {
    let (a, b) = mesh.edges[edge].tris;
    if a == cell { b } else { a }
}

impl LinearReconstruction {
    pub(super) fn prepare(
        mesh: &TriMesh,
        moments: &[[f64; 2]],
        rates: &[EdgeVolumeTransport],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let mut cells = allocate(mesh.tris.len(), LinearCell::default())?;
        for (cell, stencil) in cells.iter_mut().enumerate() {
            poll(cancelled)?;
            let points = vertices(mesh, cell)?;
            let center = centroid(points);
            let length = points
                .iter()
                .flat_map(|p| [p[0] - center[0], p[1] - center[1]])
                .fold(0.0_f64, |a, b| a.max(b.abs()));
            if !positive(length) || center.iter().any(|x| !x.is_finite()) {
                return Err(invalid("nondegenerate linear density geometry required"));
            }
            let mut vectors = [[0.0; 2]; 3];
            let mut normal = [0.0; 3];
            for (local, &(edge, _)) in mesh.tri_edges[cell].iter().enumerate() {
                let face = &mesh.edges[edge];
                let a = *mesh
                    .verts
                    .get(face.verts.0)
                    .ok_or(invalid("original density edge vertex required"))?;
                let b = *mesh
                    .verts
                    .get(face.verts.1)
                    .ok_or(invalid("original density edge vertex required"))?;
                if face.verts.0 == face.verts.1
                    || !mesh.tris[cell].contains(&face.verts.0)
                    || !mesh.tris[cell].contains(&face.verts.1)
                {
                    return Err(invalid("density face must belong to its original triangle"));
                }
                let other = neighbor(mesh, cell, edge);
                let sample = if other != usize::MAX {
                    centroid(vertices(mesh, other)?)
                } else if rates[edge].against_normal_m2_s() > 0.0 {
                    core::array::from_fn(|k| a[k] + 0.5 * (b[k] - a[k]))
                } else {
                    continue;
                };
                let d: [f64; 2] = core::array::from_fn(|k| (sample[k] - center[k]) / length);
                let square = d[0].mul_add(d[0], d[1] * d[1]);
                if !positive(square) {
                    return Err(invalid("finite distinct density stencil sites required"));
                }
                vectors[local] = d.map(|value| value / square);
                normal[0] += d[0] * vectors[local][0];
                normal[1] += d[0] * vectors[local][1];
                normal[2] += d[1] * vectors[local][1];
            }
            let determinant = normal[0].mul_add(normal[2], -normal[1] * normal[1]);
            let trace = normal[0] + normal[2];
            // A one-sided/collinear stencil cannot identify a planar slope.
            // Retain the original constant donor there; no ghost density is invented.
            if !determinant.is_finite() || determinant <= 64.0 * f64::EPSILON * trace * trace {
                continue;
            }
            for (vertex, row) in stencil.variation.iter_mut().enumerate() {
                let d: [f64; 2] =
                    core::array::from_fn(|k| (points[vertex][k] - center[k]) / length);
                let inverse = [
                    (d[0] * normal[2] - d[1] * normal[1]) / determinant,
                    (d[1] * normal[0] - d[0] * normal[1]) / determinant,
                ];
                for (weight, vector) in row.iter_mut().zip(&vectors) {
                    *weight = inverse[0].mul_add(vector[0], inverse[1] * vector[1]);
                    if !weight.is_finite() {
                        return Err(invalid("unrepresentable density reconstruction stencil"));
                    }
                }
            }
        }
        let mut fractions = allocate(mesh.edges.len(), [[0.5, 1.0 / 12.0]; 2])?;
        for (edge, target) in fractions.iter_mut().enumerate() {
            poll(cancelled)?;
            *target = flux_centroids(moments[edge]);
        }
        Ok(Self { cells, fractions })
    }

    pub(super) fn deviations(&self, cell: usize, differences: [f64; 3]) -> [f64; 3] {
        let mut deviation = [0.0; 3];
        for (value, weights) in deviation[..2].iter_mut().zip(&self.cells[cell].variation) {
            *value = weights[0].mul_add(
                differences[0],
                weights[1].mul_add(differences[1], weights[2] * differences[2]),
            );
        }
        deviation[2] = -deviation[0] - deviation[1];
        deviation
    }

    /// A specific value is a MASS mean. Shift the affine variation so its
    /// density-weighted triangle mean is zero: integral(rho*delta)/mass
    /// equals sum(rho_vertex*delta_vertex)/(4*sum(rho_vertex)).
    pub(super) fn specific_deviations(
        &self,
        cell: usize,
        differences: [f64; 3],
        density: [f64; 3],
    ) -> [f64; 3] {
        let deviation = self.deviations(cell, differences);
        let scale = density.into_iter().fold(0.0_f64, f64::max);
        let rho = density.map(|value| value / scale);
        let shift = (rho[0] - rho[2]).mul_add(deviation[0], (rho[1] - rho[2]) * deviation[1])
            / (4.0 * (rho[0] + rho[1] + rho[2]));
        deviation.map(|value| value - shift)
    }

    /// Covariance shifts the volume centroid to the mass centroid. Evaluating
    /// an affine specific trace here exactly integrates rho*s*u_n, a cubic,
    /// on this original sign interval; using the volume centroid loses that term.
    pub(super) fn mass_centroid(&self, edge: usize, direction: usize, rho: [f64; 2]) -> f64 {
        let [mean, variance] = self.fractions[edge][direction];
        let scale = rho[0].max(rho[1]);
        let [a, b] = rho.map(|value| value / scale);
        mean + (b - a) * variance / ((1.0 - mean) * a + mean * b)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn reconstruct_into<const N: usize>(
        &self,
        mesh: &TriMesh,
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        rates: &[EdgeVolumeTransport],
        vertex_density: &mut [[f64; 3]],
        vertex_bounds: &mut [[f64; 2]],
        face_density: &mut [[f64; 2]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        if cells.len() != self.cells.len()
            || vertex_density.len() != cells.len()
            || inflow.len() != self.fractions.len()
            || face_density.len() != inflow.len()
            || vertex_bounds.len() != mesh.verts.len()
        {
            return Err(invalid(
                "matching bounded linear density workspace required",
            ));
        }
        // One linear pass builds vertex-patch bounds. Face-neighbor bounds
        // alone can clip an exactly affine field at every triangle vertex.
        // The incident-cell patch preserves those slopes in the interior,
        // without searching a variable-valence stencil for every cell.
        for bound in vertex_bounds.iter_mut() {
            poll(cancelled)?;
            *bound = [f64::INFINITY, f64::NEG_INFINITY];
        }
        for (cell, state) in cells.iter().enumerate() {
            poll(cancelled)?;
            admit_cell(state)?;
            for &vertex in &mesh.tris[cell] {
                vertex_bounds[vertex][0] = vertex_bounds[vertex][0].min(state.density_kg_m3);
                vertex_bounds[vertex][1] = vertex_bounds[vertex][1].max(state.density_kg_m3);
            }
        }
        for (edge, exterior) in inflow.iter().enumerate() {
            poll(cancelled)?;
            if mesh.edges[edge].tris.1 != usize::MAX && exterior.is_some() {
                return Err(invalid("inflow states belong only to boundary edges"));
            }
            if let Some(state) = exterior {
                admit_cell(state)?;
                if rates[edge].against_normal_m2_s() > 0.0 {
                    for vertex in [mesh.edges[edge].verts.0, mesh.edges[edge].verts.1] {
                        vertex_bounds[vertex][0] =
                            vertex_bounds[vertex][0].min(state.density_kg_m3);
                        vertex_bounds[vertex][1] =
                            vertex_bounds[vertex][1].max(state.density_kg_m3);
                    }
                }
            } else if mesh.edges[edge].tris.1 == usize::MAX
                && rates[edge].against_normal_m2_s() > 0.0
            {
                return Err(invalid(
                    "every incoming boundary stream needs a physical donor state",
                ));
            }
        }
        for (cell, state) in cells.iter().enumerate() {
            poll(cancelled)?;
            admit_cell(state)?;
            let rho = state.density_kg_m3;
            let mut samples = [rho; 3];
            for (sample, &(edge, _)) in samples.iter_mut().zip(&mesh.tri_edges[cell]) {
                let other = neighbor(mesh, cell, edge);
                if other != usize::MAX {
                    *sample = cells[other].density_kg_m3;
                    if !positive(*sample) {
                        return Err(invalid("positive neighboring density required"));
                    }
                } else if rates[edge].against_normal_m2_s() > 0.0 {
                    *sample = inflow[edge]
                        .as_ref()
                        .expect("admitted incoming state")
                        .density_kg_m3;
                }
            }
            let bounds = mesh.tris[cell].map(|vertex| vertex_bounds[vertex]);
            let lower = bounds.iter().fold(rho, |value, bound| value.min(bound[0]));
            let upper = bounds.iter().fold(rho, |value, bound| value.max(bound[1]));
            let base = rho / upper;
            if lower == upper || base == 0.0 {
                vertex_density[cell] = [rho; 3];
                continue;
            }
            let differences = samples.map(|value| (value - rho) / upper);
            let deviation = self.deviations(cell, differences);
            let mut theta = 1.0_f64;
            for (&value, bound) in deviation.iter().zip(bounds) {
                if value > 0.0 {
                    theta = theta.min(((bound[1] - rho) / upper) / value);
                } else if value < 0.0 {
                    theta = theta.min(((rho - bound[0]) / upper) / -value);
                }
            }
            // Contract a limited slope very slightly towards its positive mean
            // so a large contrast cannot round a bound-touching vertex to zero.
            if theta < 1.0 {
                theta *= 1.0 - 64.0 * f64::EPSILON;
            }
            let reconstructed = deviation.map(|value| (theta * value).mul_add(upper, rho));
            vertex_density[cell] = if reconstructed.iter().all(|&value| positive(value)) {
                reconstructed
            } else {
                // An unrepresentable slope keeps the existing positive donor.
                [rho; 3]
            };
        }
        for (edge, target) in face_density.iter_mut().enumerate() {
            poll(cancelled)?;
            let face = &mesh.edges[edge];
            for (direction, donor) in [face.tris.0, face.tris.1].into_iter().enumerate() {
                target[direction] = if donor == usize::MAX {
                    inflow[edge]
                        .as_ref()
                        .map_or(cells[face.tris.0].density_kg_m3, |state| {
                            state.density_kg_m3
                        })
                } else {
                    let triangle = mesh.tris[donor];
                    let a = triangle
                        .iter()
                        .position(|&v| v == face.verts.0)
                        .ok_or(invalid("density face must belong to its original triangle"))?;
                    let b = triangle
                        .iter()
                        .position(|&v| v == face.verts.1)
                        .ok_or(invalid("density face must belong to its original triangle"))?;
                    let [first, second] = [vertex_density[donor][a], vertex_density[donor][b]];
                    if first == second {
                        first
                    } else {
                        let s = self.fractions[edge][direction][0];
                        positive_value(
                            (1.0 - s) * first + s * second,
                            "directional reconstructed density",
                        )?
                    }
                };
            }
        }
        Ok(())
    }
}

// Exact |linear normal velocity|-weighted centroid and variance on each sign
// interval. Variance also supplies the density-specific covariance term.
fn flux_centroids(moments: [f64; 2]) -> [[f64; 2]; 2] {
    let scale = moments[0].abs().max(moments[1].abs());
    if scale == 0.0 {
        return [[0.5, 1.0 / 12.0]; 2];
    }
    let mean = moments[0] / scale;
    let slope = 6.0 * (moments[1] / scale);
    let (left, right) = (mean - slope, mean + slope);
    let mut result = [[0.5, 1.0 / 12.0]; 2];
    if (left < 0.0 && right > 0.0) || (left > 0.0 && right < 0.0) {
        let split = left.abs() / (left.abs() + right.abs());
        let width = 1.0 - split;
        result[usize::from(left < 0.0)] = [split / 3.0, split * split / 18.0];
        result[usize::from(right < 0.0)] = [split + width * (2.0 / 3.0), width * width / 18.0];
    } else {
        let (a, b) = (left.abs(), right.abs());
        let direction = usize::from(left < 0.0 || right < 0.0);
        let weight = a / (a + b);
        result[direction] = [
            (a + 2.0 * b) / (3.0 * (a + b)),
            (1.0 + 2.0 * weight * (1.0 - weight)) / 18.0,
        ];
    }
    result
}
