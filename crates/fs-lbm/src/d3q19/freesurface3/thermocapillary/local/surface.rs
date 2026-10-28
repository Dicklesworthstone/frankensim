//! PLIC geometry and tangential derivatives of independently solved surface sigma.
//!
//! Surface samples live at cut-patch centroids. They must never be interpreted
//! as ambient liquid-cell samples or extrapolated again in the normal direction.

use super::*;

/// One volume-matching partial-cell patch, in lattice lengths and areas.
#[derive(Clone, Copy, Debug)]
pub(crate) struct LocalSurfacePatch3 {
    pub location: (u64, usize),
    pub liquid: [usize; 3],
    pub centroid: [f64; 3],
    pub outward_normal: [f64; 3],
    pub area: f64,
    pub volume: f64,
    // Retain the small cell-relative offset without subtracting large grid
    // coordinates when reconstructing an ambient source or clipped face.
    pub(super) offset: [f64; 3],
}

impl LocalSurfacePatch3 {
    /// Actual liquid-polyhedron boundary: six clipped cube faces plus the cut.
    /// This is an area budget, not a declaration of thermal contact topology.
    pub(crate) fn liquid_boundary_area(self) -> f64 {
        let normal = self.outward_normal;
        let offset = dot(normal, self.offset);
        let mut result = self.area;
        for axis in 0..3 {
            let tangent = match axis {
                0 => [1, 2],
                1 => [0, 2],
                _ => [0, 1],
            };
            for side in [-0.5, 0.5] {
                let corners = [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]];
                let mut clipped = [[0.0; 2]; 8];
                let mut count = 0;
                for i in 0..4 {
                    let p = corners[i];
                    let q = corners[(i + 1) % 4];
                    let distance = |r: [f64; 2]| {
                        normal[axis] * side + normal[tangent[0]] * r[0] + normal[tangent[1]] * r[1]
                            - offset
                    };
                    let (dp, dq) = (distance(p), distance(q));
                    if dp <= 0.0 {
                        clipped[count] = p;
                        count += 1;
                    }
                    if (dp < 0.0 && dq > 0.0) || (dp > 0.0 && dq < 0.0) {
                        let t = dp / (dp - dq);
                        clipped[count] = core::array::from_fn(|d| p[d] + t * (q[d] - p[d]));
                        count += 1;
                    }
                }
                let mut twice_area = 0.0;
                for i in 0..count {
                    let (p, q) = (clipped[i], clipped[(i + 1) % count]);
                    twice_area += p[0] * q[1] - p[1] * q[0];
                }
                result += 0.5 * twice_area.abs();
            }
        }
        result
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn tangent_basis(normal: [f64; 3]) -> [[f64; 3]; 2] {
    let axis = (0..3)
        .min_by(|&a, &b| normal[a].abs().total_cmp(&normal[b].abs()))
        .unwrap();
    let mut unit = [0.0; 3];
    unit[axis] = 1.0;
    let first = cross(normal, unit);
    let length = dot(first, first).sqrt();
    let first = first.map(|v| v / length);
    [first, cross(normal, first)]
}

impl FreeSurface3 {
    /// The same resolved geometry used by ambient and solved-surface traction.
    pub(crate) fn local_interface_geometry(
        &self,
        model: LocalSurface3,
        gate: &CancelGate,
    ) -> Result<Vec<LocalSurfacePatch3>, FreeSurfaceError3> {
        self.admit_local_surface_model(model)?;
        let mut patches = Vec::new();
        for (&key, tile) in &self.tiles {
            check_cancel(gate)?;
            for (lane, &cell) in tile.cells.iter().enumerate() {
                if cell != Cell3::Interface {
                    continue;
                }
                let (x, y, z) = Self::global_of(key, lane);
                let center = [x, y, z];
                let fill = self
                    .local_geometric_fill(center)
                    .ok_or(invalid(
                        "invalid liquid fill in local surface reconstruction",
                    ))?
                    .0;
                if fill == 0.0 || fill == 1.0 {
                    continue;
                }
                if fill < model.minimum_fill {
                    return Err(invalid(
                        "local surface cut cell is below its admitted liquid volume",
                    ));
                }
                let normal = self.local_height_normal(center).ok_or(invalid(
                    "local surface needs a complete well-conditioned monotone height stencil",
                ))?;
                let patch = cut_patch(normal, fill)
                    .ok_or(invalid("local surface plane has no representable cut area"))?;
                patches.push(LocalSurfacePatch3 {
                    location: (key, lane),
                    liquid: center.map(|v| v as usize),
                    centroid: core::array::from_fn(|d| center[d] as f64 + 0.5 + patch.centroid[d]),
                    outward_normal: normal,
                    area: patch.area,
                    volume: fill,
                    offset: patch.centroid,
                });
            }
        }
        if patches.is_empty() {
            return Err(invalid(
                "local surface has no resolved partial-cell patches",
            ));
        }
        Ok(patches)
    }

    /// Area-weighted tangent-plane fit of actual surface values on the local
    /// face-connected wet component. Opposing surface sheets cannot contribute.
    /// A rank-deficient or ill-conditioned two-dimensional fit refuses.
    fn local_surface_gradient(
        &self,
        patch: LocalSurfacePatch3,
        geometry: &BTreeMap<(u64, usize), LocalSurfacePatch3>,
        sigma: &BTreeMap<(u64, usize), f64>,
        maximum_condition: f64,
    ) -> Option<[f64; 3]> {
        let center = patch.liquid.map(|v| v as i64);
        let mut samples = [(0.0, (0, 0)); 27];
        for (i, sample) in samples.iter_mut().enumerate() {
            let offset = [
                (i % 3) as i64 - 1,
                (i / 3 % 3) as i64 - 1,
                (i / 9) as i64 - 1,
            ];
            *sample = self.local_geometric_fill(core::array::from_fn(|d| center[d] + offset[d]))?;
        }
        let connected = wet_component(&samples.map(|s| s.0), 13);
        let basis = tangent_basis(patch.outward_normal);
        let reference = *sigma.get(&patch.location)?;
        let (mut a, mut b, mut c) = (0.0, 0.0, 0.0);
        let mut rhs = [0.0; 2];
        for (i, &(_, location)) in samples.iter().enumerate() {
            if i == 13 || !connected[i] {
                continue;
            }
            let Some(neighbor) = geometry.get(&location) else {
                continue;
            };
            if dot(patch.outward_normal, neighbor.outward_normal) <= 0.0 {
                continue;
            }
            let delta = core::array::from_fn(|d| {
                neighbor.liquid[d] as f64 - patch.liquid[d] as f64 + neighbor.offset[d]
                    - patch.offset[d]
            });
            let tangent = basis.map(|axis| dot(axis, delta));
            let distance = tangent[0] * tangent[0] + tangent[1] * tangent[1];
            if distance <= 0.0 || !distance.is_finite() {
                return None;
            }
            let weight = neighbor.area / distance;
            let difference = sigma.get(&location)? - reference;
            a += weight * tangent[0] * tangent[0];
            b += weight * tangent[0] * tangent[1];
            c += weight * tangent[1] * tangent[1];
            for d in 0..2 {
                rhs[d] += weight * tangent[d] * difference;
            }
        }
        let norm = (a.abs() + b.abs()).max(c.abs() + b.abs());
        if norm <= 0.0 || !norm.is_finite() {
            return None;
        }
        let (a, b, c) = (a / norm, b / norm, c / norm);
        let determinant = a * c - b * b;
        let condition = (a.abs() + b.abs()).max(c.abs() + b.abs()) / determinant;
        if determinant <= 0.0 || !condition.is_finite() || condition > maximum_condition {
            return None;
        }
        let gradient = [
            (c * rhs[0] - b * rhs[1]) / (norm * determinant),
            (a * rhs[1] - b * rhs[0]) / (norm * determinant),
        ];
        let gradient: [f64; 3] =
            core::array::from_fn(|d| basis[0][d] * gradient[0] + basis[1][d] * gradient[1]);
        gradient.iter().all(|v| v.is_finite()).then_some(gradient)
    }

    /// Solved patch sigma drives both pressure and tangent stress. Bulk entries
    /// satisfy the native complete-field contract but never enter this fit.
    pub(crate) fn local_surface_marangoni_forces(
        &self,
        model: LocalSurface3,
        sigma: &mut BTreeMap<(u64, usize), f64>,
        gate: &CancelGate,
    ) -> Result<BTreeMap<(u64, usize), [f64; 3]>, FreeSurfaceError3> {
        self.admit_surface_tensions(sigma, gate)?;
        let geometry: BTreeMap<_, _> = self
            .local_interface_geometry(model, gate)?
            .into_iter()
            .map(|patch| (patch.location, patch))
            .collect();
        let mut forces = BTreeMap::new();
        let mut patches = BTreeMap::new();
        for (&location, &patch) in &geometry {
            check_cancel(gate)?;
            let gradient = self.local_surface_gradient(patch, &geometry, sigma, model.maximum_gradient_condition)
                .ok_or(invalid("local solved-surface sigma has insufficient or ill-conditioned tangent support"))?;
            let force = gradient.map(|v| v * patch.area / patch.volume);
            if force.iter().any(|v| !v.is_finite()) {
                return Err(invalid(
                    "local solved-surface traction is not representable",
                ));
            }
            if force != [0.0; 3] {
                forces.insert(location, force);
            }
            patches.insert(location, (patch.area, sigma[&location]));
        }
        self.extend_local_surface_tensions(patches, sigma, gate)?;
        Ok(forces)
    }
}
