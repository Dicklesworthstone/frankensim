//! Explicit planar bonds between independently triangulated P1 surfaces.
//!
//! Triangle intersections carry positive degree-two quadrature. The SAME
//! interpolation constructs separation, residual forces and tangent actions in
//! the existing cohesive solver. This exactly integrates an undamaged linear
//! law and preserves discrete virtual work even during nonlinear softening.
//! Nonlinear fracture still needs independent surface/load-path refinement.
//!
//! No proximity-based bond discovery occurs. Each supplied surface must be an
//! exterior, nonoverlapping triangulation of the same planar bonded area, with
//! disjoint displacement DOFs and the bodies on opposite sides. Coordinates
//! are projected to the explicitly admitted common plane within the supplied
//! relative tolerance; this is fixed small-strain reference geometry, not a
//! large-sliding contact or remeshing algorithm.

use std::collections::{BTreeMap, BTreeSet};

use super::{
    BilinearCohesive, CohesivePatch, Cx, Error, JumpTerm, MaterialStateId, P1Interfaces,
    ResolvedBilinearCohesive, SurfacePoint, TraceBudget, dot, face_without, finite, frame, length,
    poll, positive, subtract,
};

/// One explicitly bonded planar region with independent surface triangulations.
#[derive(Debug, Clone, Copy)]
pub struct NonmatchingSurface<'a> {
    pub id: u64,
    pub side_a: &'a [[usize; 3]],
    pub side_b: &'a [[usize; 3]],
    pub law: &'a BilinearCohesive,
}

/// The same region retaining its original ordered material/interface source.
#[derive(Debug, Clone, Copy)]
pub struct MaterialNonmatchingSurface<'a> {
    pub id: u64,
    pub side_a: &'a [[usize; 3]],
    pub side_b: &'a [[usize; 3]],
    pub source: &'a ResolvedBilinearCohesive,
    pub surface_materials: [&'a MaterialStateId; 2],
}

/// Explicit construction work, allocation and geometric admission limits.
#[derive(Debug, Clone, Copy)]
pub struct NonmatchingBudget {
    /// Existing mesh/DOF/quality limits; `maximum_interfaces` counts regions.
    pub mesh: TraceBudget,
    /// Total faces on both sides of all regions.
    pub maximum_surface_faces: usize,
    /// Candidate pairs with positive bounding-box overlap, including within-side
    /// overlap rejection. A deterministic sweep skips spatially disjoint pairs.
    pub maximum_pair_checks: usize,
    pub maximum_quadrature_points: usize,
    /// Relative plane, coverage and roundoff tolerance, in `(0, 1e-6]`.
    pub relative_geometry_tolerance: f64,
}

impl<'a> P1Interfaces<'a> {
    /// Bind independently meshed planar surfaces to the existing cohesive map.
    pub fn from_nonmatching(
        nodes_m: &[[f64; 3]],
        tetrahedra: &[[usize; 4]],
        free_dofs: &[usize],
        surfaces: &[NonmatchingSurface<'a>],
        budget: NonmatchingBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, Error> {
        Self::from_nonmatching_controlled(
            nodes_m,
            tetrahedra,
            free_dofs,
            surfaces,
            budget,
            &mut || cx.checkpoint().is_err(),
        )
    }

    /// Source-bound constructor with the same quadrature and work transfer.
    pub fn from_nonmatching_interfaces(
        nodes_m: &[[f64; 3]],
        tetrahedra: &[[usize; 4]],
        free_dofs: &[usize],
        surfaces: &[MaterialNonmatchingSurface<'a>],
        budget: NonmatchingBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, Error> {
        Self::from_nonmatching_interfaces_controlled(
            nodes_m,
            tetrahedra,
            free_dofs,
            surfaces,
            budget,
            &mut || cx.checkpoint().is_err(),
        )
    }

    /// Predicate form preserving the original source on every overlap point.
    pub fn from_nonmatching_interfaces_controlled<C: FnMut() -> bool>(
        nodes_m: &[[f64; 3]],
        tetrahedra: &[[usize; 4]],
        free_dofs: &[usize],
        surfaces: &[MaterialNonmatchingSurface<'a>],
        budget: NonmatchingBudget,
        cancelled: &mut C,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        if surfaces.len() > budget.mesh.maximum_interfaces {
            return Err(Error::Budget("nonmatching source interface count"));
        }
        let mut geometric = Vec::with_capacity(surfaces.len());
        let mut sources = BTreeMap::new();
        for surface in surfaces {
            poll(cancelled)?;
            surface
                .source
                .require_surface_materials(
                    surface.surface_materials[0],
                    surface.surface_materials[1],
                )
                .map_err(|error| Error::Source {
                    patch: surface.id,
                    error,
                })?;
            sources.insert(surface.id, surface.source);
            geometric.push(NonmatchingSurface {
                id: surface.id,
                side_a: surface.side_a,
                side_b: surface.side_b,
                law: surface.source.law(),
            });
        }
        let mut result = Self::from_nonmatching_controlled(
            nodes_m, tetrahedra, free_dofs, &geometric, budget, cancelled,
        )?;
        for ((full, reduced), point) in result
            .full_patches
            .iter_mut()
            .zip(&mut result.patches)
            .zip(&result.points)
        {
            poll(cancelled)?;
            full.source = Some(sources[&point.interface_id]);
            reduced.source = full.source;
        }
        Ok(result)
    }

    /// Predicate form. Cancellation is checked per face, tetrahedron, pair and
    /// overlap quadrature point; failed construction publishes no partial map.
    pub fn from_nonmatching_controlled<C: FnMut() -> bool>(
        nodes_m: &[[f64; 3]],
        tetrahedra: &[[usize; 4]],
        free_dofs: &[usize],
        surfaces: &[NonmatchingSurface<'a>],
        budget: NonmatchingBudget,
        cancelled: &mut C,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        positive(budget.relative_geometry_tolerance)?;
        positive(budget.mesh.minimum_face_quality)?;
        if budget.relative_geometry_tolerance > 1e-6 || budget.mesh.minimum_face_quality >= 1.0 {
            return Err(Error::Invalid(
                "nonmatching geometric tolerance or face quality",
            ));
        }
        if nodes_m.is_empty()
            || tetrahedra.is_empty()
            || surfaces.is_empty()
            || nodes_m.len() > budget.mesh.maximum_nodes
            || tetrahedra.len() > budget.mesh.maximum_tetrahedra
            || surfaces.len() > budget.mesh.maximum_interfaces
        {
            return Err(Error::Budget("nonmatching mesh or interface count"));
        }
        let full_dimension = nodes_m
            .len()
            .checked_mul(3)
            .ok_or(Error::Budget("P1 DOF overflow"))?;
        for chunk in nodes_m.chunks(256) {
            poll(cancelled)?;
            for node in chunk {
                for &x in node {
                    finite(x)?;
                }
            }
        }
        if free_dofs.len() > full_dimension {
            return Err(Error::Invalid("too many free P1 DOFs"));
        }
        let mut reduced_of = vec![None; full_dimension];
        for (i, &dof) in free_dofs.iter().enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            if dof >= full_dimension || (i > 0 && dof <= free_dofs[i - 1]) {
                return Err(Error::Invalid(
                    "P1 free DOFs must be distinct, sorted and in range",
                ));
            }
            reduced_of[dof] = Some(i);
        }
        let mut incidence = BTreeMap::<[usize; 3], (usize, usize)>::new();
        let mut ids = BTreeSet::new();
        let mut face_count = 0usize;
        for surface in surfaces {
            poll(cancelled)?;
            if surface.side_a.is_empty() || surface.side_b.is_empty() || !ids.insert(surface.id) {
                return Err(Error::Invalid("empty or repeated nonmatching interface"));
            }
            face_count = face_count
                .checked_add(surface.side_a.len())
                .and_then(|n| n.checked_add(surface.side_b.len()))
                .ok_or(Error::Budget("nonmatching face count overflow"))?;
            if face_count > budget.maximum_surface_faces {
                return Err(Error::Budget("nonmatching surface faces"));
            }
            let mut side_a_nodes = BTreeSet::new();
            for (side, faces) in [surface.side_a, surface.side_b].into_iter().enumerate() {
                for face in faces {
                    poll(cancelled)?;
                    let mut key = *face;
                    key.sort_unstable();
                    if key[2] >= nodes_m.len()
                        || key[0] == key[1]
                        || key[1] == key[2]
                        || incidence.insert(key, (0, 0)).is_some()
                    {
                        return Err(Error::Invalid("invalid or reused nonmatching face"));
                    }
                    for node in key {
                        if side == 0 {
                            side_a_nodes.insert(node);
                        } else if side_a_nodes.contains(&node) {
                            return Err(Error::Invalid(
                                "cohesive traces must have distinct vertices",
                            ));
                        }
                    }
                }
            }
        }
        for tet in tetrahedra {
            poll(cancelled)?;
            let mut sorted = *tet;
            sorted.sort_unstable();
            if sorted[3] >= nodes_m.len() || sorted.windows(2).any(|w| w[0] == w[1]) {
                return Err(Error::Invalid("invalid tetrahedral connectivity"));
            }
            for opposite in 0..4 {
                if let Some((count, vertex)) = incidence.get_mut(&face_without(*tet, opposite)) {
                    *count = count.saturating_add(1);
                    *vertex = tet[opposite];
                }
            }
        }
        if incidence.values().any(|&(count, _)| count != 1) {
            return Err(Error::Invalid(
                "cohesive face must have one incident tetrahedron",
            ));
        }
        let mut result = Self {
            full_dimension,
            free_dofs: free_dofs.to_vec(),
            reduced_of,
            full_patches: Vec::new(),
            patches: Vec::new(),
            points: Vec::new(),
        };
        let mut checks = 0usize;
        for surface in surfaces {
            poll(cancelled)?;
            let mut first_a = surface.side_a[0];
            first_a.sort_unstable();
            let mut first_b = surface.side_b[0];
            first_b.sort_unstable();
            let (axes, _) = frame(
                nodes_m,
                first_a,
                incidence[&first_a].1,
                incidence[&first_b].1,
                budget.mesh.minimum_face_quality,
            )?;
            let origin = nodes_m[first_a[0]];
            let scale = length(subtract(nodes_m[first_a[1]], origin))?
                .max(length(subtract(nodes_m[first_a[2]], origin))?);
            let plane_tolerance = positive(scale * budget.relative_geometry_tolerance)?;
            let mut sides = [Vec::new(), Vec::new()];
            for (side, faces) in [surface.side_a, surface.side_b].into_iter().enumerate() {
                for face in faces {
                    poll(cancelled)?;
                    let mut key = *face;
                    key.sort_unstable();
                    let distance =
                        finite(dot(axes[0], subtract(nodes_m[incidence[&key].1], origin)))?;
                    if (side == 0 && distance >= -plane_tolerance)
                        || (side == 1 && distance <= plane_tolerance)
                    {
                        return Err(Error::Invalid(
                            "bonded bodies must lie on opposite sides of the common plane",
                        ));
                    }
                    sides[side].push(Triangle::new(
                        nodes_m,
                        key,
                        origin,
                        axes,
                        plane_tolerance,
                        budget.mesh.minimum_face_quality,
                    )?);
                }
            }
            // Reject overlap within either trace: area sums alone could hide a
            // gap compensated by a doubled region elsewhere.
            for side in &sides {
                visit_pairs(
                    side,
                    None,
                    &mut checks,
                    budget.maximum_pair_checks,
                    cancelled,
                    &mut |i, j, _| {
                        let area = polygon_area(&intersection(&side[i], &side[j])?)?;
                        if area
                            > budget.relative_geometry_tolerance * side[i].area.min(side[j].area)
                        {
                            return Err(Error::Invalid("overlapping faces on one cohesive trace"));
                        }
                        Ok(())
                    },
                )?;
            }
            let mut covered_a = vec![0.0; sides[0].len()];
            let mut covered_b = vec![0.0; sides[1].len()];
            let point_start = result.points.len();
            visit_pairs(
                &sides[0],
                Some(&sides[1]),
                &mut checks,
                budget.maximum_pair_checks,
                cancelled,
                &mut |i, j, cancelled| {
                    let a = &sides[0][i];
                    let b = &sides[1][j];
                    let polygon = intersection(a, b)?;
                    for k in 1..polygon.len().saturating_sub(1) {
                        let vertices = [polygon[0], polygon[k], polygon[k + 1]];
                        let area = finite(
                            0.5 * cross2(
                                sub2(vertices[1], vertices[0]),
                                sub2(vertices[2], vertices[0]),
                            )
                            .abs(),
                        )?;
                        if area == 0.0 {
                            continue;
                        }
                        covered_a[i] = finite(covered_a[i] + area)?;
                        covered_b[j] = finite(covered_b[j] + area)?;
                        for q in 0..3 {
                            poll(cancelled)?;
                            if result.points.len() >= budget.maximum_quadrature_points {
                                return Err(Error::Budget("nonmatching quadrature points"));
                            }
                            let xy = core::array::from_fn(|c| {
                                (0..3)
                                    .map(|n| {
                                        vertices[n][c] * if n == q { 2.0 / 3.0 } else { 1.0 / 6.0 }
                                    })
                                    .sum()
                            });
                            let sa = a.shape(xy, budget.relative_geometry_tolerance)?;
                            let sb = b.shape(xy, budget.relative_geometry_tolerance)?;
                            let mut terms = Vec::with_capacity(18);
                            for (triangle, shape, sign) in [(a, sa, -1.0), (b, sb, 1.0)] {
                                for n in 0..3 {
                                    for c in 0..3 {
                                        if shape[n] != 0.0 {
                                            terms.push(JumpTerm {
                                                dof: 3 * triangle.nodes[n] + c,
                                                local: core::array::from_fn(|r| {
                                                    sign * shape[n] * axes[r][c]
                                                }),
                                            });
                                        }
                                    }
                                }
                            }
                            let reduced = terms
                                .iter()
                                .filter_map(|term| {
                                    result.reduced_of[term.dof].map(|dof| JumpTerm {
                                        dof,
                                        local: term.local,
                                    })
                                })
                                .collect();
                            let id = u64::try_from(result.points.len())
                                .map_err(|_| Error::Budget("nonmatching patch identity"))?;
                            let weight = positive(area / 3.0)?;
                            result.full_patches.push(CohesivePatch::new(
                                id,
                                weight,
                                surface.law,
                                terms,
                            )?);
                            result.patches.push(CohesivePatch::new(
                                id,
                                weight,
                                surface.law,
                                reduced,
                            )?);
                            let mut location_m = [0.0; 3];
                            for c in 0..3 {
                                location_m[c] =
                                    finite(origin[c] + xy[0] * axes[1][c] + xy[1] * axes[2][c])?;
                            }
                            result.points.push(SurfacePoint {
                                interface_id: surface.id,
                                point_in_face: result.points.len() - point_start,
                                location_m,
                                area_weight_m2: weight,
                                axes,
                                shape_values: sa,
                                side_b_shape_values: sb,
                                paired_nodes: [a.nodes, b.nodes],
                            });
                        }
                    }
                    Ok(())
                },
            )?;
            for (triangles, coverage) in sides.iter().zip([covered_a, covered_b]) {
                for (triangle, area) in triangles.iter().zip(coverage) {
                    poll(cancelled)?;
                    if (area - triangle.area).abs()
                        > budget.relative_geometry_tolerance * triangle.area
                    {
                        return Err(Error::Invalid(
                            "nonmatching traces do not cover the same bonded area",
                        ));
                    }
                }
            }
        }
        poll(cancelled)?;
        Ok(result)
    }
}

struct Triangle {
    nodes: [usize; 3],
    xy: [[f64; 2]; 3],
    area: f64,
    lower: [f64; 2],
    upper: [f64; 2],
}
impl Triangle {
    fn new(
        nodes_m: &[[f64; 3]],
        mut nodes: [usize; 3],
        origin: [f64; 3],
        axes: [[f64; 3]; 3],
        plane_tolerance: f64,
        quality: f64,
    ) -> Result<Self, Error> {
        let mut xy = [[0.0; 2]; 3];
        for i in 0..3 {
            let delta = subtract(nodes_m[nodes[i]], origin);
            if finite(dot(axes[0], delta))?.abs() > plane_tolerance {
                return Err(Error::Invalid(
                    "nonmatching interface is not planar within tolerance",
                ));
            }
            xy[i] = [finite(dot(axes[1], delta))?, finite(dot(axes[2], delta))?];
        }
        let mut doubled = finite(cross2(sub2(xy[1], xy[0]), sub2(xy[2], xy[0])))?;
        if doubled < 0.0 {
            nodes.swap(1, 2);
            xy.swap(1, 2);
            doubled = -doubled;
        }
        positive(doubled)?;
        let mut edge2 = 0.0_f64;
        for i in 0..3 {
            let edge = sub2(xy[(i + 1) % 3], xy[i]);
            edge2 = edge2.max(finite(edge[0] * edge[0] + edge[1] * edge[1])?);
        }
        if doubled / positive(edge2)? < quality {
            return Err(Error::Invalid("under-quality nonmatching face"));
        }
        let lower = core::array::from_fn(|c| xy.iter().map(|p| p[c]).fold(f64::INFINITY, f64::min));
        let upper =
            core::array::from_fn(|c| xy.iter().map(|p| p[c]).fold(f64::NEG_INFINITY, f64::max));
        Ok(Self {
            nodes,
            xy,
            area: positive(doubled / 2.0)?,
            lower,
            upper,
        })
    }
    fn shape(&self, p: [f64; 2], tolerance: f64) -> Result<[f64; 3], Error> {
        let v = sub2(p, self.xy[0]);
        let b = finite(cross2(v, sub2(self.xy[2], self.xy[0])) / (2.0 * self.area))?;
        let c = finite(cross2(sub2(self.xy[1], self.xy[0]), v) / (2.0 * self.area))?;
        let mut result = [1.0 - b - c, b, c];
        if result
            .iter()
            .any(|&x| x < -tolerance || x > 1.0 + tolerance)
        {
            return Err(Error::Invalid("unstable overlap barycentric coordinates"));
        }
        for x in &mut result {
            *x = x.max(0.0);
        }
        let total: f64 = result.iter().sum();
        Ok(result.map(|x| x / total))
    }
}

// Sweep along the wider bounding-box extent, with stable input-index ties.
// Active boxes can still be quadratic for long overlapping slivers, but exact
// clipping and the pair budget are spent only on boxes overlapping in 2D.
fn visit_pairs<C: FnMut() -> bool, V: FnMut(usize, usize, &mut C) -> Result<(), Error>>(
    first: &[Triangle],
    second: Option<&[Triangle]>,
    checks: &mut usize,
    maximum: usize,
    cancelled: &mut C,
    visit: &mut V,
) -> Result<(), Error> {
    let sides = [first, second.unwrap_or(&[])];
    let mut lower = [f64::INFINITY; 2];
    let mut upper = [f64::NEG_INFINITY; 2];
    let mut events = Vec::with_capacity(first.len() + sides[1].len());
    for (side, triangles) in sides.iter().enumerate() {
        for (index, triangle) in triangles.iter().enumerate() {
            poll(cancelled)?;
            events.push((side, index));
            for c in 0..2 {
                lower[c] = lower[c].min(triangle.lower[c]);
                upper[c] = upper[c].max(triangle.upper[c]);
            }
        }
    }
    let axis = usize::from(upper[1] - lower[1] > upper[0] - lower[0]);
    let other = 1 - axis;
    events.sort_by(|&(sa, ia), &(sb, ib)| {
        sides[sa][ia].lower[axis]
            .total_cmp(&sides[sb][ib].lower[axis])
            .then(sa.cmp(&sb))
            .then(ia.cmp(&ib))
    });
    let mut active = [Vec::<usize>::new(), Vec::<usize>::new()];
    for (side, index) in events {
        poll(cancelled)?;
        let triangle = &sides[side][index];
        let partner = if second.is_some() { 1 - side } else { 0 };
        // Compact active intervals in place, polling at bounded scan tiles.
        for s in 0..2 {
            let mut kept = 0;
            for read in 0..active[s].len() {
                if read % 256 == 0 {
                    poll(cancelled)?;
                }
                let candidate = active[s][read];
                if sides[s][candidate].upper[axis] > triangle.lower[axis] {
                    active[s][kept] = candidate;
                    kept += 1;
                }
            }
            active[s].truncate(kept);
        }
        for &candidate in &active[partner] {
            poll(cancelled)?;
            let b = &sides[partner][candidate];
            if triangle.upper[other] <= b.lower[other] || b.upper[other] <= triangle.lower[other] {
                continue;
            }
            admit_pair(checks, maximum, cancelled)?;
            if side == 0 {
                visit(index, candidate, cancelled)?;
            } else {
                visit(candidate, index, cancelled)?;
            }
        }
        active[side].push(index);
    }
    Ok(())
}

fn admit_pair<C: FnMut() -> bool>(
    checks: &mut usize,
    maximum: usize,
    cancelled: &mut C,
) -> Result<(), Error> {
    poll(cancelled)?;
    *checks = checks
        .checked_add(1)
        .ok_or(Error::Budget("nonmatching pair count overflow"))?;
    if *checks > maximum {
        return Err(Error::Budget("nonmatching triangle pair checks"));
    }
    Ok(())
}
fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}
fn cross2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn polygon_area(polygon: &[[f64; 2]]) -> Result<f64, Error> {
    let mut sum = 0.0;
    for i in 1..polygon.len().saturating_sub(1) {
        sum = finite(
            sum + 0.5
                * cross2(
                    sub2(polygon[i], polygon[0]),
                    sub2(polygon[i + 1], polygon[0]),
                )
                .abs(),
        )?;
    }
    Ok(sum)
}
fn intersection(a: &Triangle, b: &Triangle) -> Result<Vec<[f64; 2]>, Error> {
    if (0..2).any(|c| a.upper[c] <= b.lower[c] || b.upper[c] <= a.lower[c]) {
        return Ok(Vec::new());
    }
    let mut polygon = a.xy.to_vec();
    for i in 0..3 {
        if polygon.is_empty() {
            break;
        }
        let origin = b.xy[i];
        let edge = sub2(b.xy[(i + 1) % 3], origin);
        let mut clipped = Vec::with_capacity(8);
        let mut previous = polygon[polygon.len() - 1];
        let mut previous_distance = finite(cross2(edge, sub2(previous, origin)))?;
        for current in polygon {
            let distance = finite(cross2(edge, sub2(current, origin)))?;
            if (distance >= 0.0) != (previous_distance >= 0.0) {
                let t = finite(previous_distance / finite(previous_distance - distance)?)?;
                let point = [
                    finite(previous[0] + t * (current[0] - previous[0]))?,
                    finite(previous[1] + t * (current[1] - previous[1]))?,
                ];
                clipped.push(point);
            }
            if distance >= 0.0 {
                clipped.push(current);
            }
            previous = current;
            previous_distance = distance;
        }
        polygon = clipped;
    }
    Ok(polygon)
}

#[cfg(test)]
mod tests;
