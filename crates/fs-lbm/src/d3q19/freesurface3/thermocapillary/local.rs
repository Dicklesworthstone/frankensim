//! Resolved local thermocapillary traction for detached liquid interfaces.
//!
//! Raw VOF height columns share their geometry owner with native curvature.
//! A volume-matching plane cuts each partial cell; its area localizes the
//! tangential stress in that cell's occupied liquid volume. An affine fit of
//! connected, wet, cell-centered sigma samples supplies the ambient gradient
//! and an extrapolated surface coefficient. Independently solved surface values
//! instead use the separate tangent-plane fit at the same PLIC centroids; that
//! path neither samples bulk sigma nor extrapolates surface values again.

use super::super::curvature;
use super::*;

mod surface;
pub(crate) use surface::LocalSurfacePatch3;

/// Resolution envelope for local, detached-interface reconstruction.
#[derive(Clone, Copy, Debug)]
pub(crate) struct LocalSurface3 {
    /// Positive lower bound below one for every nonzero partial-cell fill.
    pub minimum_fill: f64,
    /// Upper bound, at least one, on the infinity-norm condition number of
    /// the ambient three-dimensional or surface-tangent two-dimensional fit.
    pub maximum_gradient_condition: f64,
}

#[derive(Clone, Copy)]
struct CutPatch {
    area: f64,
    centroid: [f64; 3],
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// Volume below a plane with sorted positive normal components, in [0,1]^3.
// Reflection confines alpha to the lower half. The piecewise convolution
// avoids subtracting nearly equal cubes when one normal component is small.
fn lower_cube_volume(sorted_lengths: [f64; 3], alpha: f64) -> f64 {
    let [a, b, c] = sorted_lengths;
    if a == 0.0 {
        return if b == 0.0 {
            alpha / c
        } else if alpha < b {
            (alpha / b) * (alpha / c) * 0.5
        } else {
            (alpha - 0.5 * b) / c
        };
    }
    if alpha < a {
        return (alpha / a) * (alpha / b) * (alpha / c) / 6.0;
    }
    if alpha < b {
        return (alpha * (alpha - a) + a * a / 3.0) / (2.0 * b * c);
    }
    let upper_corner = (a + b - alpha).max(0.0);
    let lower_corner = (alpha - c).max(0.0);
    (alpha - 0.5 * (a + b)) / c
        + ((upper_corner / a) * (upper_corner / b) * (upper_corner / c)
            - (lower_corner / a) * (lower_corner / b) * (lower_corner / c))
            / 6.0
}

/// Plane/cube intersection, with its offset fixed by the actual VOF volume.
/// All arrays are bounded by the twelve cube edges. Empty/full cells own no
/// patch, so interface closure ghosts cannot duplicate a surface force.
fn cut_patch(normal: [f64; 3], fill: f64) -> Option<CutPatch> {
    if !(0.0..1.0).contains(&fill) || fill == 0.0 {
        return None;
    }
    let lengths = normal.map(f64::abs);
    let total = lengths.iter().sum::<f64>();
    let mut ordered = lengths;
    ordered.sort_unstable_by(f64::total_cmp);
    let target = fill.min(1.0 - fill);
    let (mut lower, mut upper) = (0.0, 0.5 * total);
    // Fixed upper work bound; an unresolved tiny volume refuses below rather
    // than acquiring the area of the final representable bracket.
    for _ in 0..64 {
        let middle = 0.5 * (lower + upper);
        if middle == lower || middle == upper {
            break;
        }
        if lower_cube_volume(ordered, middle) < target {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    let alpha = 0.5 * (lower + upper);
    if (lower_cube_volume(ordered, alpha) - target).abs() > 1e-8 * target {
        return None;
    }
    let offset = if fill <= 0.5 { alpha } else { total - alpha } - 0.5 * total;
    let mut points = [[0.0; 3]; 12];
    let mut count = 0;
    let mut add = |point: [f64; 3]| {
        if points[..count]
            .iter()
            .any(|old| (0..3).map(|d| (old[d] - point[d]).powi(2)).sum::<f64>() < 1e-28)
        {
            return;
        }
        if count < points.len() {
            points[count] = point;
            count += 1;
        }
    };
    for axis in 0..3 {
        let tangent = match axis {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        };
        for a in [-0.5, 0.5] {
            for b in [-0.5, 0.5] {
                let mut start = [0.0; 3];
                start[axis] = -0.5;
                start[tangent[0]] = a;
                start[tangent[1]] = b;
                let mut end = start;
                end[axis] = 0.5;
                let first = dot(normal, start) - offset;
                let last = dot(normal, end) - offset;
                if first == 0.0 {
                    add(start);
                }
                if last == 0.0 {
                    add(end);
                }
                if (first < 0.0 && last > 0.0) || (first > 0.0 && last < 0.0) {
                    let fraction = first / (first - last);
                    start[axis] += fraction;
                    add(start);
                }
            }
        }
    }
    if count < 3 {
        return None;
    }
    let center: [f64; 3] =
        core::array::from_fn(|d| points[..count].iter().map(|p| p[d]).sum::<f64>() / count as f64);
    let mut axis = 0;
    for d in 1..3 {
        if lengths[d] > lengths[axis] {
            axis = d;
        }
    }
    let tangent = match axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    };
    // Angular order in a nonsingular Cartesian projection, without trig.
    points[..count].sort_unstable_by(|p, q| {
        let p = [
            p[tangent[0]] - center[tangent[0]],
            p[tangent[1]] - center[tangent[1]],
        ];
        let q = [
            q[tangent[0]] - center[tangent[0]],
            q[tangent[1]] - center[tangent[1]],
        ];
        let half = |r: [f64; 2]| r[1] < 0.0 || (r[1] == 0.0 && r[0] < 0.0);
        half(p)
            .cmp(&half(q))
            .then_with(|| (q[0] * p[1] - p[0] * q[1]).total_cmp(&0.0))
    });
    let mut area = 0.0;
    let mut centroid = [0.0; 3];
    for i in 0..count {
        let p = points[i];
        let q = points[(i + 1) % count];
        let u: [f64; 3] = core::array::from_fn(|d| p[d] - center[d]);
        let v: [f64; 3] = core::array::from_fn(|d| q[d] - center[d]);
        let cross = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let triangle = 0.5 * dot(cross, normal).abs();
        area += triangle;
        for d in 0..3 {
            centroid[d] += triangle * (center[d] + p[d] + q[d]) / 3.0;
        }
    }
    if !area.is_finite() || area <= 0.0 {
        return None;
    }
    for value in &mut centroid {
        *value /= area;
    }
    Some(CutPatch { area, centroid })
}

fn solve_gradient(
    matrix: [[f64; 3]; 3],
    rhs: [f64; 3],
    maximum_condition: f64,
) -> Option<[f64; 3]> {
    let norm = matrix
        .iter()
        .map(|row| row.iter().map(|x| x.abs()).sum::<f64>())
        .fold(0.0_f64, f64::max);
    if !norm.is_finite() || norm <= 0.0 {
        return None;
    }
    let m = matrix.map(|row| row.map(|x| x / norm));
    let [a, b, c] = m[0];
    let [_, d, e] = m[1];
    let f = m[2][2];
    let adj = [
        [d * f - e * e, c * e - b * f, b * e - c * d],
        [c * e - b * f, a * f - c * c, b * c - a * e],
        [b * e - c * d, b * c - a * e, a * d - b * b],
    ];
    let determinant = a * adj[0][0] + b * adj[0][1] + c * adj[0][2];
    if !determinant.is_finite() || determinant <= 0.0 {
        return None;
    }
    let condition = adj
        .iter()
        .map(|row| row.iter().map(|x| x.abs()).sum::<f64>() / determinant)
        .fold(0.0_f64, f64::max);
    if !condition.is_finite() || condition > maximum_condition {
        return None;
    }
    let rhs = rhs.map(|x| x / norm);
    let gradient = adj.map(|row| dot(row, rhs) / determinant);
    gradient.iter().all(|x| x.is_finite()).then_some(gradient)
}

// Cell-scale liquid connectivity, not a proof of subcell PLIC face overlap.
fn wet_component(fill: &[f64; 27], seed: usize) -> [bool; 27] {
    let mut connected = [false; 27];
    if fill[seed] <= 0.0 {
        return connected;
    }
    let mut queue = [0; 27];
    queue[0] = seed;
    connected[seed] = true;
    let (mut front, mut back) = (0, 1);
    while front < back {
        let index = queue[front];
        front += 1;
        let p = [index % 3, index / 3 % 3, index / 9];
        for axis in 0..3 {
            for shift in [-1, 1] {
                let v = p[axis] as i32 + shift;
                if !(0..3).contains(&v) {
                    continue;
                }
                let mut q = p;
                q[axis] = v as usize;
                let next = q[0] + 3 * q[1] + 9 * q[2];
                if !connected[next] && fill[next] > 0.0 {
                    connected[next] = true;
                    queue[back] = next;
                    back += 1;
                }
            }
        }
    }
    connected
}

impl FreeSurface3 {
    fn local_geometric_fill(&self, p: [i64; 3]) -> Option<(f64, (u64, usize))> {
        let (nx, ny, nz) = self.grid.tile_dims();
        let size = [nx * TILE, ny * TILE, nz * TILE];
        if (0..3).any(|d| p[d] < 0 || p[d] as usize >= size[d]) {
            return None;
        }
        let Some((slot, lane)) = self.grid.resolve_source(p[0], p[1], p[2]) else {
            // The sparse owner allocates a short gas halo, not every gas cell
            // in a wider height stencil. An inactive interior tile is empty
            // gas unless explicitly solid; it is not an unresolved wall.
            let q = p.map(|x| x as usize);
            let key = morton3(
                (q[0] / TILE) as u32,
                (q[1] / TILE) as u32,
                (q[2] / TILE) as u32,
            );
            if self.solid.contains(&key) {
                return None;
            }
            let lane = q[0] % TILE + TILE * (q[1] % TILE + TILE * (q[2] % TILE));
            return Some((0.0, (key, lane)));
        };
        let key = self.grid.active_keys()[slot];
        let tile = &self.tiles[&key];
        let fill = match tile.cells[lane] {
            Cell3::Gas => 0.0,
            Cell3::Fluid => 1.0,
            Cell3::Interface => {
                let rho = self.grid.populations(slot, lane).iter().sum::<f64>();
                if !rho.is_finite() || rho <= 0.0 {
                    return None;
                }
                tile.mass[lane] / rho
            }
        };
        if !fill.is_finite() || fill < 0.0 || fill > 1.0 + 64.0 * f64::EPSILON {
            return None;
        }
        Some((fill.min(1.0), (key, lane)))
    }

    fn local_height_normal(&self, center: [i64; 3]) -> Option<[f64; 3]> {
        let mut direction = [0.0; 3];
        for axis in 0..3 {
            let mut p = center;
            p[axis] -= 1;
            let below = self.local_geometric_fill(p)?.0;
            p[axis] += 2;
            direction[axis] = below - self.local_geometric_fill(p)?.0;
        }
        let mut axes = [0, 1, 2];
        // Directions subtract unit fills, whose population-density reductions
        // can move by a few ulps under a momentum-preserving force-frame change.
        // Resolve only that arithmetic ambiguity in coordinate order; choosing
        // a different height graph can otherwise rotate a curved-cell normal
        // by percent-scale amounts. Select against each remaining maximum,
        // rather than using a nontransitive tolerance-based sort comparator.
        for first in 0..axes.len() {
            let largest = axes[first..]
                .iter()
                .map(|&axis| direction[axis].abs())
                .fold(0.0_f64, f64::max);
            let selected = (first..axes.len())
                .filter(|&i| largest - direction[axes[i]].abs() <= 64.0 * f64::EPSILON)
                .min_by_key(|&i| axes[i])
                .expect("nonempty remaining height axes");
            axes.swap(first, selected);
        }
        for axis in axes {
            if direction[axis] == 0.0 {
                continue;
            }
            let mut orientation = [0.0; 3];
            orientation[axis] = direction[axis].signum();
            for half_depth in [3, 5, 7] {
                if let Some(geometry) =
                    curvature::height_geometry_with_depth(orientation, half_depth, |offset| {
                        self.local_geometric_fill(core::array::from_fn(|d| center[d] + offset[d]))
                            .map(|sample| sample.0)
                    })
                {
                    // Corner cuts can need longer columns even for a plane.
                    // The fixed maximum is fifteen axial cells; incomplete
                    // geometry has no smoothed-normal fallback.
                    // At a curved corner cut, the different Cartesian graphs
                    // sample slightly different surface points. Requiring an
                    // exactly largest normal component can reject every axis.
                    // Bound the graph metric by four instead (|n_axis| >= 1/2).
                    if geometry.normal[axis].abs() >= 0.5 {
                        return Some(geometry.normal);
                    }
                }
            }
        }
        None
    }

    fn local_sigma_gradient(
        &self,
        center: [i64; 3],
        sigma: &BTreeMap<(u64, usize), f64>,
        maximum_condition: f64,
    ) -> Option<[f64; 3]> {
        let mut samples = [(0.0, 0.0); 27];
        for (i, sample) in samples.iter_mut().enumerate() {
            let offset = [
                (i % 3) as i64 - 1,
                (i / 3 % 3) as i64 - 1,
                (i / 9) as i64 - 1,
            ];
            let (fill, location) =
                self.local_geometric_fill(core::array::from_fn(|d| center[d] + offset[d]))?;
            if fill > 0.0 {
                *sample = (fill, *sigma.get(&location)?);
            }
        }
        // Nearby liquid across a gas gap cannot supply a thermal gradient.
        // Only the face-connected wet component of the target is sampled.
        let connected = wet_component(&samples.map(|sample| sample.0), 13);
        let reference = samples[13].1;
        let mut matrix = [[0.0; 3]; 3];
        let mut rhs = [0.0; 3];
        for (i, &(fill, value)) in samples.iter().enumerate() {
            if !connected[i] || i == 13 {
                continue;
            }
            let offset = [
                (i % 3) as f64 - 1.0,
                (i / 3 % 3) as f64 - 1.0,
                (i / 9) as f64 - 1.0,
            ];
            let weight = fill / dot(offset, offset);
            for a in 0..3 {
                rhs[a] += weight * offset[a] * (value - reference);
                for b in 0..3 {
                    matrix[a][b] += weight * offset[a] * offset[b];
                }
            }
        }
        solve_gradient(matrix, rhs, maximum_condition)
    }

    /// Tangential traction on resolved local interface patches, including
    /// detached drops and surfaces that are not one global height graph.
    ///
    /// The native Guo owner receives (I-nn)grad(sigma)*area/fill; normal Laplace
    /// pressure remains exclusively in its existing boundary owner. Coefficients
    /// are reconstructed at PLIC centroids and extended to empty/full closure
    /// ghosts from neighboring partial patches. The input map is published only
    /// after complete admission; no field or lattice state changes on refusal.
    pub(crate) fn local_marangoni_forces(
        &self,
        model: LocalSurface3,
        sigma: &mut BTreeMap<(u64, usize), f64>,
        gate: &CancelGate,
    ) -> Result<BTreeMap<(u64, usize), [f64; 3]>, FreeSurfaceError3> {
        self.admit_surface_tensions(sigma, gate)?;
        let geometry: Vec<LocalSurfacePatch3> = self.local_interface_geometry(model, gate)?;
        let mut patches = BTreeMap::new();
        let mut forces = BTreeMap::new();
        for patch in geometry {
            check_cancel(gate)?;
            let center = patch.liquid.map(|v| v as i64);
            let gradient = self.local_sigma_gradient(center, sigma, model.maximum_gradient_condition).ok_or(invalid("local surface sigma gradient has insufficient or ill-conditioned connected wet support"))?;
            let surface_sigma = sigma[&patch.location] + dot(gradient, patch.offset);
            if !surface_sigma.is_finite()
                || surface_sigma <= 0.0
                || surface_sigma * Self::CURVATURE_LIMIT / CS2_3 >= 1.0
            {
                return Err(invalid(
                    "extrapolated surface tension exceeds native pressure-boundary domain",
                ));
            }
            let normal_gradient = dot(gradient, patch.outward_normal);
            let force = core::array::from_fn(|d| {
                (gradient[d] - normal_gradient * patch.outward_normal[d]) * patch.area
                    / patch.volume
            });
            if force.iter().any(|f| !f.is_finite()) {
                return Err(invalid("local surface traction is not representable"));
            }
            if force != [0.0; 3] {
                forces.insert(patch.location, force);
            }
            patches.insert(patch.location, (patch.area, surface_sigma));
        }
        self.extend_local_surface_tensions(patches, sigma, gate)?;
        Ok(forces)
    }

    fn admit_local_surface_model(&self, model: LocalSurface3) -> Result<(), FreeSurfaceError3> {
        if !model.minimum_fill.is_finite()
            || !(0.0..1.0).contains(&model.minimum_fill)
            || model.minimum_fill == 0.0
            || !model.maximum_gradient_condition.is_finite()
            || model.maximum_gradient_condition < 1.0
            || !self.solid.is_empty()
            || self.planar_contact.is_some()
            || self.planar_contact_field.is_some()
            || self.contact != ContactModel3::Neutral
            || self.carry != 0.0
        {
            return Err(invalid(
                "local thermocapillarity needs resolved liquid, neutral walls and bounded reconstruction conditioning",
            ));
        }
        Ok(())
    }

    fn extend_local_surface_tensions(
        &self,
        patches: BTreeMap<(u64, usize), (f64, f64)>,
        sigma: &mut BTreeMap<(u64, usize), f64>,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        let mut ghosts = Vec::new();
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
                    ghosts.push(((key, lane), center));
                }
            }
        }
        if patches.is_empty() {
            return Err(invalid(
                "local surface has no resolved partial-cell patches",
            ));
        }
        let mut extension = Vec::with_capacity(ghosts.len());
        for (location, center) in ghosts {
            check_cancel(gate)?;
            let mut samples = [(0.0, None); 27];
            for (i, sample) in samples.iter_mut().enumerate() {
                let offset = [
                    (i % 3) as i64 - 1,
                    (i / 3 % 3) as i64 - 1,
                    (i / 9) as i64 - 1,
                ];
                if let Some((fill, neighbor)) =
                    self.local_geometric_fill(core::array::from_fn(|d| center[d] + offset[d]))
                {
                    *sample = (fill, patches.get(&neighbor).copied());
                }
            }
            let occupied = samples[13].0 > 0.0;
            let seed = if occupied {
                13
            } else {
                samples
                    .iter()
                    .position(|sample| sample.1.is_some())
                    .ok_or(invalid(
                        "empty interface ghost has no neighboring resolved surface patch",
                    ))?
            };
            let connected = wet_component(&samples.map(|sample| sample.0), seed);
            let mut reference = None;
            let (mut weight, mut difference) = (0.0, 0.0);
            for (i, &(_, patch)) in samples.iter().enumerate() {
                let Some((area, value)) = patch else { continue };
                if !connected[i] {
                    if !occupied {
                        return Err(invalid(
                            "empty interface ghost touches disconnected surface components",
                        ));
                    }
                    continue;
                }
                let first = *reference.get_or_insert(value);
                let offset = [
                    (i % 3) as f64 - 1.0,
                    (i / 3 % 3) as f64 - 1.0,
                    (i / 9) as f64 - 1.0,
                ];
                let w = area / dot(offset, offset);
                weight += w;
                difference += w * (value - first);
            }
            let first = reference.ok_or(invalid(
                "interface closure ghost has no neighboring resolved surface patch",
            ))?;
            let value = first + difference / weight;
            if !value.is_finite() || value <= 0.0 || value * Self::CURVATURE_LIMIT / CS2_3 >= 1.0 {
                return Err(invalid(
                    "extended surface tension exceeds native pressure-boundary domain",
                ));
            }
            extension.push((location, value));
        }
        check_cancel(gate)?;
        for (location, (_, value)) in patches {
            sigma.insert(location, value);
        }
        for (location, value) in extension {
            sigma.insert(location, value);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const NODES: [f64; 4] = [
        -0.861_136_311_594_052_6,
        -0.339_981_043_584_856_3,
        0.339_981_043_584_856_3,
        0.861_136_311_594_052_6,
    ];
    const WEIGHTS: [f64; 4] = [
        0.347_854_845_137_453_8,
        0.652_145_154_862_546_1,
        0.652_145_154_862_546_1,
        0.347_854_845_137_453_8,
    ];

    fn model() -> LocalSurface3 {
        LocalSurface3 {
            minimum_fill: 1e-10,
            maximum_gradient_condition: 1e4,
        }
    }

    fn field(
        surface: &FreeSurface3,
        value: impl Fn([f64; 3]) -> f64,
    ) -> BTreeMap<(u64, usize), f64> {
        let mut result = BTreeMap::new();
        for (&key, tile) in &surface.tiles {
            for (lane, &cell) in tile.cells.iter().enumerate() {
                if cell != Cell3::Gas {
                    let (x, y, z) = FreeSurface3::global_of(key, lane);
                    result.insert(
                        (key, lane),
                        value([x as f64 + 0.5, y as f64 + 0.5, z as f64 + 0.5]),
                    );
                }
            }
        }
        result
    }

    fn assign_fill(surface: &mut FreeSurface3, fill: impl Fn([i64; 3]) -> f64) {
        for (&key, tile) in &mut surface.tiles {
            let slot = surface.grid.slot_of(key).unwrap();
            for lane in 0..TILE_CELLS {
                let (x, y, z) = FreeSurface3::global_of(key, lane);
                let fraction = fill([x, y, z]);
                let rho = surface.grid.populations(slot, lane).iter().sum::<f64>();
                tile.mass[lane] = fraction * rho;
                tile.cells[lane] = if fraction == 0.0 {
                    Cell3::Gas
                } else if fraction == 1.0 {
                    Cell3::Fluid
                } else {
                    Cell3::Interface
                };
            }
        }
        let mut closure = Vec::new();
        for (&key, tile) in &surface.tiles {
            for (lane, &cell) in tile.cells.iter().enumerate() {
                if cell == Cell3::Fluid
                    && (1..Q3).any(|q| {
                        matches!(
                            surface.neighbor_state(key, lane, q),
                            Neighbor::Active(_, _, Cell3::Gas)
                        )
                    })
                {
                    closure.push((key, lane));
                }
            }
        }
        for (key, lane) in closure {
            surface.tiles.get_mut(&key).unwrap().cells[lane] = Cell3::Interface;
        }
    }

    // Independent tensor Gauss integration of the exact sphere's z interval.
    // This does not call the height, PLIC-volume or traction implementation.
    fn sphere_fill(p: [f64; 3], radius: f64) -> f64 {
        let square = radius * radius;
        if p.iter()
            .map(|x| (x.abs() - 0.5).max(0.0).powi(2))
            .sum::<f64>()
            >= square
        {
            return 0.0;
        }
        if p.iter().map(|x| (x.abs() + 0.5).powi(2)).sum::<f64>() <= square {
            return 1.0;
        }
        let mut volume = 0.0;
        // Composite quadrature resolves the clipping contour inside a cell;
        // four unsplit nodes leave a planar-cut volume bias under grid refinement.
        for i in 0..32 {
            for j in 0..32 {
                let (a, b) = (i % 4, j % 4);
                let u = p[0] - 0.5 + (i as f64 / 4.0).floor() / 8.0 + (0.5 + 0.5 * NODES[a]) / 8.0;
                let v = p[1] - 0.5 + (j as f64 / 4.0).floor() / 8.0 + (0.5 + 0.5 * NODES[b]) / 8.0;
                let radial = square - u * u - v * v;
                if radial > 0.0 {
                    let height = radial.sqrt();
                    let interval = ((p[2] + 0.5).min(height) - (p[2] - 0.5).max(-height)).max(0.0);
                    volume += WEIGHTS[a] * WEIGHTS[b] * interval / 256.0;
                }
            }
        }
        volume.min(1.0)
    }

    fn drop_surface(radius: usize) -> (FreeSurface3, [f64; 3]) {
        let size = (2 * radius + 18).div_ceil(TILE) * TILE;
        let center = size / 2;
        let mut surface = FreeSurface3::new(
            size,
            size,
            size,
            0.8,
            [0.0; 3],
            0.01,
            ContactModel3::Neutral,
            |x, y, z| {
                [x, y, z]
                    .iter()
                    .map(|&p| (p as f64 - center as f64).powi(2))
                    .sum::<f64>()
                    < (radius * radius) as f64
            },
        )
        .unwrap();
        assign_fill(&mut surface, |p| {
            sphere_fill(p.map(|x| x as f64 - center as f64), radius as f64)
        });
        (surface, [center as f64 + 0.5; 3])
    }

    #[test]
    fn g1_local_plane_cut_area_and_centroid_match_exact_geometry() {
        for axis in 0..3 {
            for sign in [-1.0, 1.0] {
                let mut normal = [0.0; 3];
                normal[axis] = sign;
                for fill in [0.001, 0.2, 0.5, 0.8, 0.999] {
                    let patch = cut_patch(normal, fill).unwrap();
                    assert!((patch.area - 1.0).abs() < 1e-14);
                    assert!((patch.centroid[axis] - sign * (fill - 0.5)).abs() < 4e-15);
                    for d in 0..3 {
                        if d != axis {
                            assert!(patch.centroid[d].abs() < 1e-14);
                        }
                    }
                }
            }
        }
        let diagonal = [1.0 / 3.0_f64.sqrt(); 3];
        let patch = cut_patch(diagonal, 0.5).unwrap();
        // The centered body-diagonal plane cuts a regular hexagon.
        assert!((patch.area - 3.0_f64.sqrt() * 0.75).abs() < 1e-14);
        assert!(patch.centroid.iter().all(|x| x.abs() < 1e-14));
        let patch = cut_patch(diagonal, 1.0 / 6.0).unwrap();
        assert!((patch.area - 3.0_f64.sqrt() * 0.5).abs() < 2e-14);
        assert!(patch.centroid.iter().all(|x| (x + 1.0 / 6.0).abs() < 1e-14));
        for tiny in [1e-8, 1e-14, 1e-100] {
            let lengths = [tiny, tiny, 1.0];
            for fraction in [0.001, 0.1, 0.5] {
                let alpha = fraction + tiny;
                assert!((lower_cube_volume(lengths, alpha) - fraction).abs() < 2e-16);
            }
        }
        assert!(cut_patch(diagonal, 0.0).is_none());
        assert!(cut_patch(diagonal, 1.0).is_none());
    }

    #[test]
    fn g1_local_rotated_planes_reject_normal_only_linear_gradients() {
        for axis in 0..3 {
            let tangent = match axis {
                0 => [1, 2],
                1 => [0, 2],
                _ => [0, 1],
            };
            for sign in [-1.0, 1.0] {
                let mut surface = FreeSurface3::new(
                    24,
                    24,
                    24,
                    0.8,
                    [0.0; 3],
                    0.01,
                    ContactModel3::Neutral,
                    |_, _, z| z < 12,
                )
                .unwrap();
                let (hu, hv) = (0.2, -0.3);
                assign_fill(&mut surface, |p| {
                    let mut value = 0.0;
                    for a in 0..4 {
                        for b in 0..4 {
                            let u = p[tangent[0]] as f64 + 0.5 + 0.5 * NODES[a] - 12.5;
                            let v = p[tangent[1]] as f64 + 0.5 + 0.5 * NODES[b] - 12.5;
                            let below = (12.3 + hu * u + hv * v - p[axis] as f64).clamp(0.0, 1.0);
                            value += 0.25
                                * WEIGHTS[a]
                                * WEIGHTS[b]
                                * if sign > 0.0 { below } else { 1.0 - below };
                        }
                    }
                    value.clamp(0.0, 1.0)
                });
                let mut exact = [0.0; 3];
                exact[axis] = sign;
                exact[tangent[0]] = -sign * hu;
                exact[tangent[1]] = -sign * hv;
                let length = dot(exact, exact).sqrt();
                exact = exact.map(|n| n / length);
                let normal = surface.local_height_normal([12; 3]).unwrap();
                assert!((0..3).all(|d| (normal[d] - exact[d]).abs() < 4e-14));
                let sigma = field(&surface, |p| {
                    0.01 + 0.0001 * dot(exact, p.map(|x| x - 12.5))
                });
                let gradient = surface.local_sigma_gradient([12; 3], &sigma, 1e4).unwrap();
                let projected = core::array::from_fn::<_, 3, _>(|d| {
                    gradient[d] - dot(gradient, normal) * normal[d]
                });
                assert!(
                    projected.iter().all(|x| x.abs() < 2e-17),
                    "{axis}/{sign}: {projected:?}"
                );
                let uniform = field(&surface, |_| 0.01);
                assert_eq!(
                    surface.local_sigma_gradient([12; 3], &uniform, 1e4),
                    Some([0.0; 3])
                );
            }
        }
    }

    #[test]
    fn g1_local_detached_drop_uniform_sigma_and_closure_ghosts_add_no_force() {
        let gate = CancelGate::new_clock_free();
        let (mut surface, center) = drop_surface(16);
        let mut ghosts = Vec::new();
        for (&key, tile) in &surface.tiles {
            for (lane, &cell) in tile.cells.iter().enumerate() {
                if cell == Cell3::Interface
                    && surface.fill(key, lane) > 0.0
                    && surface.fill(key, lane) < 1.0
                {
                    for q in 1..Q3 {
                        if let Neighbor::Active(neighbor, index, Cell3::Gas) =
                            surface.neighbor_state(key, lane, q)
                        {
                            ghosts.push((neighbor, index));
                            break;
                        }
                    }
                    if !ghosts.is_empty() {
                        break;
                    }
                }
            }
            if !ghosts.is_empty() {
                break;
            }
        }
        assert!(!ghosts.is_empty());
        for &(key, lane) in &ghosts {
            surface.tiles.get_mut(&key).unwrap().cells[lane] = Cell3::Interface;
        }
        let mut sigma = field(&surface, |_| 0.01);
        assert!(
            surface
                .local_marangoni_forces(model(), &mut sigma, &gate)
                .unwrap()
                .is_empty()
        );
        assert!(sigma.values().all(|&x| x == 0.01));
        let mut sigma = field(&surface, |p| 0.01 + 0.00001 * (p[0] - center[0]));
        let forces = surface
            .local_marangoni_forces(model(), &mut sigma, &gate)
            .unwrap();
        assert!(!forces.is_empty());
        assert!(ghosts.iter().all(|location| !forces.contains_key(location)));
        assert!(forces.keys().all(|&(key, lane)| {
            let fill = surface.fill(key, lane);
            fill > 0.0 && fill < 1.0
        }));
    }

    #[test]
    fn g1_force_frame_rebind_preserves_roundoff_tied_local_height_geometry() {
        let gate = CancelGate::new_clock_free();
        let (original, center) = drop_surface(16);
        let geometry: BTreeMap<_, _> = original
            .local_interface_geometry(model(), &gate)
            .unwrap()
            .into_iter()
            .map(|patch| (patch.location, patch))
            .collect();
        let wet: Vec<_> = field(&original, |_| 0.0).into_keys().collect();
        // Cyclic source rotations exercise each Cartesian ranking. The larger
        // source exposes a second tied-column flip without advancing any time.
        for gradient in [
            [1e-6, -0.7e-6, 0.3e-6],
            [0.3e-6, 1e-6, -0.7e-6],
            [-0.7e-6, 0.3e-6, 1e-6],
            [1e-4, -0.7e-4, 0.3e-4],
        ] {
            let mut surface = original.clone();
            let mut sigma = field(&surface, |p| {
                0.01 + dot(gradient, core::array::from_fn(|d| p[d] - center[d]))
            });
            let traction = surface
                .local_marangoni_forces(model(), &mut sigma, &gate)
                .unwrap();
            assert!(!traction.is_empty());
            let forces = wet
                .iter()
                .map(|&location| {
                    (
                        location,
                        traction.get(&location).copied().unwrap_or([0.0; 3]),
                    )
                })
                .collect();
            surface.rebind_surface_tensions(sigma, &gate).unwrap();
            surface.rebind_force_field(forces, &gate).unwrap();
            assert_ne!(surface.grid().state_bits(), original.grid().state_bits());
            assert_eq!(
                surface.ledger_mass().to_bits(),
                original.ledger_mass().to_bits()
            );
            for &(key, lane) in &wet {
                assert_eq!(
                    surface.tiles[&key].cells[lane],
                    original.tiles[&key].cells[lane]
                );
                assert_eq!(
                    surface.tiles[&key].mass[lane].to_bits(),
                    original.tiles[&key].mass[lane].to_bits()
                );
                let (x, y, z) = FreeSurface3::global_of(key, lane);
                let before = original
                    .sample_at(x as usize, y as usize, z as usize)
                    .unwrap();
                let after = surface
                    .sample_at(x as usize, y as usize, z as usize)
                    .unwrap();
                assert!((after.density - before.density).abs() < 4e-15);
                assert!((0..3).all(|d| (after.velocity[d] - before.velocity[d]).abs() < 1e-14));
            }
            let after = surface.local_interface_geometry(model(), &gate).unwrap();
            assert_eq!(after.len(), geometry.len());
            for patch in after {
                let before = geometry[&patch.location];
                let normal_change = (0..3)
                    .map(|d| (patch.outward_normal[d] - before.outward_normal[d]).abs())
                    .fold(0.0_f64, f64::max);
                assert!(
                    normal_change <= 64.0 * f64::EPSILON,
                    "force-frame-only height switch at {:?}: {:?} -> {:?}",
                    patch.liquid,
                    before.outward_normal,
                    patch.outward_normal
                );
                assert!((patch.volume - before.volume).abs() <= 4.0 * f64::EPSILON);
                // Almost-full corner cuts amplify a one-ulp fill change in
                // area. This absolute cell-scale bound covers that observed
                // arithmetic sensitivity; normals retain the much tighter gate.
                assert!((patch.area - before.area).abs() <= 256.0 * f64::EPSILON);
            }
        }
    }

    #[test]
    fn g1_local_spherical_manufactured_tangential_force_refines() {
        let gate = CancelGate::new_clock_free();
        let (mut previous, mut previous_radial) = (f64::INFINITY, f64::INFINITY);
        for radius in [16, 32] {
            let (surface, center) = drop_surface(radius);
            let gradient = [
                0.0001 / radius as f64,
                -0.00006 / radius as f64,
                0.00004 / radius as f64,
            ];
            let mut sigma = field(&surface, |p| {
                0.01 + dot(gradient, core::array::from_fn(|d| p[d] - center[d]))
            });
            let forces = surface
                .local_marangoni_forces(model(), &mut sigma, &gate)
                .unwrap();
            let mut total = [0.0; 3];
            let mut area = 0.0;
            for (&(key, lane), &force) in &forces {
                let (x, y, z) = FreeSurface3::global_of(key, lane);
                let normal = surface.local_height_normal([x, y, z]).unwrap();
                assert!(
                    dot(force, normal).abs() < 2e-15,
                    "normal force leaked: {force:?}/{normal:?}"
                );
                let fill = surface.fill(key, lane);
                let patch = cut_patch(normal, fill).unwrap();
                area += patch.area;
                let point: [f64; 3] = core::array::from_fn(|d| {
                    [x, y, z][d] as f64 + 0.5 + patch.centroid[d] - center[d]
                });
                assert!((sigma[&(key, lane)] - (0.01 + dot(gradient, point))).abs() < 2e-16);
                for d in 0..3 {
                    total[d] += fill * force[d];
                }
            }
            // Independent surface integral: integral(I-nn)dA = (8pi R²/3) I.
            let scale = 8.0 * core::f64::consts::PI * (radius * radius) as f64 / 3.0;
            let relative = (0..3)
                .map(|d| (total[d] / (scale * gradient[d]) - 1.0).abs())
                .fold(0.0_f64, f64::max);
            let area_error =
                (area / (4.0 * core::f64::consts::PI * (radius * radius) as f64) - 1.0).abs();
            assert!(relative < 0.02 && area_error < 0.02);
            assert!(relative < previous * 0.8, "{relative} >= {previous}");
            previous = relative;

            // Radial sigma varies only normally at the exact spherical surface.
            // Its reconstructed tangential leakage must converge to zero.
            let amplitude = 0.0001;
            let mut radial_sigma = field(&surface, |p| {
                0.01 + amplitude * (0..3).map(|d| (p[d] - center[d]).powi(2)).sum::<f64>()
                    / (radius * radius) as f64
            });
            let radial_forces = surface
                .local_marangoni_forces(model(), &mut radial_sigma, &gate)
                .unwrap();
            let mut error_square = 0.0;
            for (&(key, lane), &force) in &radial_forces {
                let (x, y, z) = FreeSurface3::global_of(key, lane);
                let normal = surface.local_height_normal([x, y, z]).unwrap();
                let fill = surface.fill(key, lane);
                let patch = cut_patch(normal, fill).unwrap();
                error_square += dot(force, force) * fill * fill / patch.area;
            }
            let radial_error = (error_square / area).sqrt() / (2.0 * amplitude / radius as f64);
            eprintln!(
                "local sphere R={radius}: force error={relative:e}, area error={area_error:e}, radial leakage={radial_error:e}"
            );
            assert!(
                radial_error < 0.15 && radial_error < previous_radial * 0.8,
                "{radial_error} >= {previous_radial}"
            );
            previous_radial = radial_error;
        }
    }

    #[test]
    fn g0_local_source_and_reconstruction_refusals_preserve_candidates() {
        let gate = CancelGate::new_clock_free();
        let (surface, _) = drop_surface(16);
        let original = field(&surface, |_| 0.01);
        let mut missing = original.clone();
        missing.pop_first();
        let before = missing.clone();
        assert!(
            surface
                .local_marangoni_forces(model(), &mut missing, &gate)
                .is_err()
        );
        assert_eq!(missing, before);
        for value in [0.0, -0.01, CS2_3] {
            let mut invalid_source = original.clone();
            *invalid_source.first_entry().unwrap().get_mut() = value;
            let before = invalid_source.clone();
            assert!(
                surface
                    .local_marangoni_forces(model(), &mut invalid_source, &gate)
                    .is_err()
            );
            assert_eq!(invalid_source, before);
        }
        for model in [
            LocalSurface3 {
                minimum_fill: 0.9,
                ..model()
            },
            LocalSurface3 {
                maximum_gradient_condition: 1.0,
                ..model()
            },
        ] {
            let mut sigma = original.clone();
            assert!(
                surface
                    .local_marangoni_forces(model, &mut sigma, &gate)
                    .is_err()
            );
            assert_eq!(sigma, original);
        }
        assert!(
            solve_gradient([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0; 3]], [1.0; 3], 1e4).is_none()
        );
        // A sharp grid-face interface owns no strict partial-cell patches and
        // must refuse; it must never silently report zero Marangoni traction.
        let mut aligned = FreeSurface3::new(
            24,
            24,
            24,
            0.8,
            [0.0; 3],
            0.01,
            ContactModel3::Neutral,
            |_, _, z| z < 12,
        )
        .unwrap();
        assign_fill(&mut aligned, |p| if p[2] < 12 { 1.0 } else { 0.0 });
        let mut sigma = field(&aligned, |p| 0.01 + 0.00001 * p[0]);
        let before = sigma.clone();
        assert!(
            aligned
                .local_marangoni_forces(model(), &mut sigma, &gate)
                .is_err()
        );
        assert_eq!(sigma, before);
        let unresolved = FreeSurface3::new(
            16,
            16,
            16,
            0.8,
            [0.0; 3],
            0.01,
            ContactModel3::Neutral,
            |x, y, z| (6..10).contains(&x) && (6..10).contains(&y) && (6..10).contains(&z),
        )
        .unwrap();
        let mut sigma = field(&unresolved, |_| 0.01);
        let before = sigma.clone();
        assert!(
            unresolved
                .local_marangoni_forces(model(), &mut sigma, &gate)
                .is_err()
        );
        assert_eq!(sigma, before);
    }

    #[test]
    fn g1_local_cut_liquid_boundary_area_matches_independent_slabs() {
        for axis in 0..3 {
            for sign in [-1.0, 1.0] {
                for fill in [0.001, 0.2, 0.5, 0.8, 0.999] {
                    let mut normal = [0.0; 3];
                    normal[axis] = sign;
                    let cut = cut_patch(normal, fill).unwrap();
                    let patch = LocalSurfacePatch3 {
                        location: (0, 0),
                        liquid: [0; 3],
                        centroid: cut.centroid,
                        outward_normal: normal,
                        area: cut.area,
                        volume: fill,
                        offset: cut.centroid,
                    };
                    assert!((patch.liquid_boundary_area() - (2.0 + 4.0 * fill)).abs() < 3e-14);
                }
            }
        }
        let normal = [1.0 / 3.0_f64.sqrt(); 3];
        let cut = cut_patch(normal, 0.5).unwrap();
        let patch = LocalSurfacePatch3 {
            location: (0, 0),
            liquid: [0; 3],
            centroid: cut.centroid,
            outward_normal: normal,
            area: cut.area,
            volume: 0.5,
            offset: cut.centroid,
        };
        assert!((patch.liquid_boundary_area() - (3.0 + 0.75 * 3.0_f64.sqrt())).abs() < 3e-14);
    }

    #[test]
    fn g1_local_solved_surface_sigma_ignores_unrelated_bulk_gradient() {
        let gate = CancelGate::new_clock_free();
        let (surface, center) = drop_surface(16);
        let geometry = surface.local_interface_geometry(model(), &gate).unwrap();
        let mut sigma = field(&surface, |p| 0.01 + 0.00001 * (p[0] - center[0]));
        for patch in &geometry {
            sigma.insert(patch.location, 0.01);
        }
        assert!(
            surface
                .local_surface_marangoni_forces(model(), &mut sigma, &gate)
                .unwrap()
                .is_empty()
        );
        for patch in &geometry {
            assert_eq!(sigma[&patch.location], 0.01);
        }

        let gradient = [1e-6, -0.6e-6, 0.4e-6];
        let mut sigma = field(&surface, |p| 0.02 + 0.0001 * (p[1] - center[1]));
        for patch in &geometry {
            sigma.insert(
                patch.location,
                0.01 + dot(
                    gradient,
                    core::array::from_fn(|d| patch.centroid[d] - center[d]),
                ),
            );
        }
        let before = sigma.clone();
        let forces = surface
            .local_surface_marangoni_forces(model(), &mut sigma, &gate)
            .unwrap();
        let mut total = [0.0; 3];
        for patch in &geometry {
            assert_eq!(sigma[&patch.location], before[&patch.location]);
            let force = forces.get(&patch.location).copied().unwrap_or([0.0; 3]);
            assert!(dot(force, patch.outward_normal).abs() < 2e-16);
            for d in 0..3 {
                total[d] += force[d] * patch.volume;
            }
        }
        let scale = 8.0 * core::f64::consts::PI * 16.0_f64.powi(2) / 3.0;
        for d in 0..3 {
            assert!((total[d] / (scale * gradient[d]) - 1.0).abs() < 0.03);
        }

        let mut invalid = before;
        let original = invalid.clone();
        assert!(
            surface
                .local_surface_marangoni_forces(
                    LocalSurface3 {
                        maximum_gradient_condition: 1.0,
                        ..model()
                    },
                    &mut invalid,
                    &gate
                )
                .is_err()
        );
        assert_eq!(invalid, original);
    }
}
