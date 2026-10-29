//! A bounded-depth importance hierarchy for physical thermal mesh emitters.
//!
//! Internal nodes approximate incident luminance only to choose a proposal.
//! Actual radiance and leaf area densities remain exact. The same conditional
//! probabilities drive descent and the endpoint-support query used by MIS.

use std::collections::BTreeMap;

use fs_geom::Point3;
use fs_math::det;

use super::thermal::AdmittedThermalTriangle;
use super::{LightingError, TAU, clamp_unit_open};

/// Keep the exact solid-angle distribution and stream for small scenes.
const LINEAR_TRIANGLE_LIMIT: usize = 8;
/// A defensive branch probability keeps every positive leaf representable
/// through the maximum possible depth of this count-balanced tree. This is a
/// proposal choice, not a floor on physical emission or an angular cutoff.
const MIN_BRANCH_PROBABILITY: f64 = 1.0 / 65_536.0;

#[derive(Clone, Copy, Debug)]
struct Bounds {
    lo: [f64; 3],
    hi: [f64; 3],
}

impl Bounds {
    fn triangle(triangle: &AdmittedThermalTriangle) -> Self {
        let mut result = Self {
            lo: [f64::INFINITY; 3],
            hi: [f64::NEG_INFINITY; 3],
        };
        for point in triangle.source.vertices {
            for (axis, coordinate) in [point.x, point.y, point.z].into_iter().enumerate() {
                result.lo[axis] = result.lo[axis].min(coordinate);
                result.hi[axis] = result.hi[axis].max(coordinate);
            }
        }
        result
    }

    fn union(self, other: Self) -> Self {
        Self {
            lo: core::array::from_fn(|axis| self.lo[axis].min(other.lo[axis])),
            hi: core::array::from_fn(|axis| self.hi[axis].max(other.hi[axis])),
        }
    }

    fn center(self) -> [f64; 3] {
        core::array::from_fn(|axis| 0.5 * self.lo[axis] + 0.5 * self.hi[axis])
    }

    fn contains(self, point: Point3) -> bool {
        let coordinates = [point.x, point.y, point.z];
        let scale = coordinates
            .into_iter()
            .chain(self.lo)
            .chain(self.hi)
            .map(f64::abs)
            .fold(0.0_f64, f64::max);
        // The narrow phase admits rounded plane/barycentric coordinates.
        // This broader envelope includes that entire tolerance, including
        // shared edges, translated geometry, and exactly zero-thickness boxes.
        let tolerance = 1024.0 * f64::EPSILON * scale;
        coordinates.into_iter().enumerate().all(|(axis, value)| {
            value >= self.lo[axis] - tolerance && value <= self.hi[axis] + tolerance
        })
    }

    fn morton_key(self, domain: Self) -> u64 {
        let center = self.center();
        let mut key = 0_u64;
        for axis in 0..3 {
            let scale = domain.lo[axis].abs().max(domain.hi[axis].abs());
            let coordinate = if domain.lo[axis] == domain.hi[axis] || scale == 0.0 {
                0
            } else {
                let lo = domain.lo[axis] / scale;
                let extent = domain.hi[axis] / scale - lo;
                let fraction = ((center[axis] / scale - lo) / extent).clamp(0.0, 1.0);
                (fraction * ((1_u32 << 21) - 1) as f64) as u32
            };
            for bit in 0..21 {
                key |= u64::from((coordinate >> bit) & 1) << (3 * bit + axis);
            }
        }
        key
    }
}

#[derive(Clone, Copy, Debug)]
enum Kind {
    Leaf(usize),
    Branch { left: usize, right: usize },
}

#[derive(Clone, Copy, Debug)]
struct Node {
    kind: Kind,
    bounds: Bounds,
    /// Log of summed triangle area times cached visible luminance. The log
    /// representation prevents geometric scale alone overflowing admission.
    log_power: f64,
    power: f64,
    weight_cap: f64,
    log_weight_cap: f64,
    primitive_min: usize,
    primitive_max: usize,
}

impl Node {
    fn leaf(index: usize, triangle: &AdmittedThermalTriangle) -> Self {
        let weight_cap = TAU * triangle.luminance;
        Self {
            kind: Kind::Leaf(index),
            bounds: Bounds::triangle(triangle),
            log_power: det::ln(triangle.area) + det::ln(triangle.luminance),
            power: triangle.area * triangle.luminance,
            weight_cap,
            log_weight_cap: det::ln(weight_cap),
            primitive_min: triangle.source.primitive_index,
            primitive_max: triangle.source.primitive_index,
        }
    }

    fn weight(self, triangles: &[AdmittedThermalTriangle], origin: Point3) -> f64 {
        if let Kind::Leaf(index) = self.kind {
            return triangles[index].weight(origin);
        }
        let coordinates = [origin.x, origin.y, origin.z];
        if coordinates.iter().any(|value| !value.is_finite()) {
            return 0.0;
        }
        let center = self.bounds.center();
        let mut distance_squared = 0.0;
        for axis in 0..3 {
            let distance = coordinates[axis] - center[axis];
            let radius = 0.5 * self.bounds.hi[axis] - 0.5 * self.bounds.lo[axis];
            distance_squared += distance * distance + radius * radius;
        }
        if self.power > 0.0
            && self.power.is_finite()
            && distance_squared > 0.0
            && distance_squared.is_finite()
        {
            // Ordinary scene scales need no logarithms or exponentials in
            // internal-node queries. The logarithmic state is a fallback for
            // otherwise representable extreme physical sizes and powers.
            return (self.power / distance_squared)
                .max(f64::MIN_POSITIVE)
                .min(self.weight_cap);
        }
        let scale = coordinates
            .into_iter()
            .chain(self.bounds.lo)
            .chain(self.bounds.hi)
            .map(f64::abs)
            .fold(0.0_f64, f64::max);
        // A non-degenerate triangle guarantees a nonzero spatial scale.
        // Normalize before subtracting/squaring: finite metre coordinates
        // may otherwise overflow even when the resulting proposal is usable.
        let mut normalized_distance_squared = 0.0;
        for axis in 0..3 {
            let distance = coordinates[axis] / scale - center[axis] / scale;
            let radius = 0.5 * (self.bounds.hi[axis] / scale - self.bounds.lo[axis] / scale);
            normalized_distance_squared += distance * distance + radius * radius;
        }
        let log_weight =
            self.log_power - (2.0 * det::ln(scale) + det::ln(normalized_distance_squared));
        // The finite cap avoids a near-source singularity. It is a sampling
        // heuristic, not a bound or replacement for a cluster's radiance.
        if log_weight >= self.log_weight_cap {
            self.weight_cap
        } else {
            det::exp(log_weight)
                .max(f64::MIN_POSITIVE)
                .min(self.weight_cap)
        }
    }
}

#[derive(Debug)]
pub(super) struct ThermalLightTree {
    /// Preorder layout; each subtree occupies one contiguous node range.
    nodes: Vec<Node>,
    leaf_for_triangle: Vec<usize>,
}

impl ThermalLightTree {
    pub(super) fn try_new(
        triangles: &[AdmittedThermalTriangle],
        checkpoint: &mut impl FnMut() -> Result<(), LightingError>,
    ) -> Result<Option<Self>, LightingError> {
        let mut domain = None;
        let mut count = 0_usize;
        for triangle in triangles {
            checkpoint()?;
            if triangle.luminance == 0.0 {
                continue;
            }
            let bounds = Bounds::triangle(triangle);
            domain = Some(domain.map_or(bounds, |previous: Bounds| previous.union(bounds)));
            count += 1;
        }
        if count <= LINEAR_TRIANGLE_LIMIT {
            return Ok(None);
        }
        let domain = domain.ok_or(LightingError::EnvironmentTooLarge)?;
        // Ordered insertion is cancellable, unlike a single unbounded sort.
        // Geometry/state identities break Morton-cell ties independently of
        // caller order; identical same-primitive proposals remain equivalent.
        let mut ordered = BTreeMap::new();
        for (index, triangle) in triangles.iter().enumerate() {
            checkpoint()?;
            if triangle.luminance > 0.0 {
                let key = Bounds::triangle(triangle).morton_key(domain);
                ordered.insert((key, triangle.identity, index), index);
            }
        }
        let mut indices = Vec::new();
        indices
            .try_reserve_exact(count)
            .map_err(|_| LightingError::EnvironmentTooLarge)?;
        for index in ordered.into_values() {
            checkpoint()?;
            indices.push(index);
        }
        let node_count = count
            .checked_mul(2)
            .and_then(|value| value.checked_sub(1))
            .ok_or(LightingError::EnvironmentTooLarge)?;
        let mut nodes = Vec::new();
        nodes
            .try_reserve_exact(node_count)
            .map_err(|_| LightingError::EnvironmentTooLarge)?;
        let mut leaf_for_triangle = Vec::new();
        leaf_for_triangle
            .try_reserve_exact(triangles.len())
            .map_err(|_| LightingError::EnvironmentTooLarge)?;
        for _ in triangles {
            checkpoint()?;
            leaf_for_triangle.push(usize::MAX);
        }
        let mut result = Self {
            nodes,
            leaf_for_triangle,
        };
        result.build(&indices, triangles, checkpoint)?;
        checkpoint()?;
        Ok(Some(result))
    }

    fn build(
        &mut self,
        indices: &[usize],
        triangles: &[AdmittedThermalTriangle],
        checkpoint: &mut impl FnMut() -> Result<(), LightingError>,
    ) -> Result<usize, LightingError> {
        checkpoint()?;
        let index = self.nodes.len();
        self.nodes
            .push(Node::leaf(indices[0], &triangles[indices[0]]));
        if indices.len() == 1 {
            self.leaf_for_triangle[indices[0]] = index;
            return Ok(index);
        }
        let midpoint = indices.len() / 2;
        let left = self.build(&indices[..midpoint], triangles, checkpoint)?;
        let right = self.build(&indices[midpoint..], triangles, checkpoint)?;
        let a = self.nodes[left];
        let b = self.nodes[right];
        let maximum_power = a.log_power.max(b.log_power);
        let log_power = maximum_power + det::ln(1.0 + det::exp(-(a.log_power - b.log_power).abs()));
        self.nodes[index] = Node {
            kind: Kind::Branch { left, right },
            bounds: a.bounds.union(b.bounds),
            log_power,
            power: a.power + b.power,
            weight_cap: a.weight_cap.max(b.weight_cap),
            log_weight_cap: a.log_weight_cap.max(b.log_weight_cap),
            primitive_min: a.primitive_min.min(b.primitive_min),
            primitive_max: a.primitive_max.max(b.primitive_max),
        };
        Ok(index)
    }

    pub(super) fn weight(&self, triangles: &[AdmittedThermalTriangle], origin: Point3) -> f64 {
        self.nodes[0].weight(triangles, origin)
    }

    fn child_probabilities(
        &self,
        triangles: &[AdmittedThermalTriangle],
        origin: Point3,
        left: usize,
        right: usize,
    ) -> Option<[f64; 2]> {
        let a = self.nodes[left].weight(triangles, origin);
        let b = self.nodes[right].weight(triangles, origin);
        if a == 0.0 && b == 0.0 {
            return None;
        }
        if a == 0.0 {
            return Some([0.0, 1.0]);
        }
        if b == 0.0 {
            return Some([1.0, 0.0]);
        }
        let ratio = a.min(b) / a.max(b);
        let smaller = (ratio / (1.0 + ratio)).max(MIN_BRANCH_PROBABILITY);
        Some(if a <= b {
            [smaller, 1.0 - smaller]
        } else {
            [1.0 - smaller, smaller]
        })
    }

    /// One descent and a residual uniform variate for the triangle's area map.
    pub(super) fn sample(
        &self,
        triangles: &[AdmittedThermalTriangle],
        origin: Point3,
        mut u: f64,
    ) -> Option<(usize, f64)> {
        let mut node = 0;
        loop {
            match self.nodes[node].kind {
                Kind::Leaf(index) => return Some((index, u)),
                Kind::Branch { left, right } => {
                    let probabilities = self.child_probabilities(triangles, origin, left, right)?;
                    // Give the smaller probability the first interval so that
                    // subtracting it from one cannot erase its sampling support.
                    let first = usize::from(probabilities[0] > probabilities[1]);
                    let (choice, residual) = if u < probabilities[first] {
                        (first, u / probabilities[first])
                    } else {
                        let second = 1 - first;
                        (second, (u - probabilities[first]) / probabilities[second])
                    };
                    node = [left, right][choice];
                    u = clamp_unit_open(residual);
                }
            }
        }
    }

    pub(super) fn selection_probability(
        &self,
        triangles: &[AdmittedThermalTriangle],
        origin: Point3,
        triangle: usize,
    ) -> f64 {
        let Some(&target) = self.leaf_for_triangle.get(triangle) else {
            return 0.0;
        };
        if target == usize::MAX {
            return 0.0;
        }
        let mut node = 0;
        let mut probability = 1.0;
        while let Kind::Branch { left, right } = self.nodes[node].kind {
            let Some(probabilities) = self.child_probabilities(triangles, origin, left, right)
            else {
                return 0.0;
            };
            let choice = usize::from(target >= right);
            probability *= probabilities[choice];
            node = [left, right][choice];
        }
        probability
    }

    pub(super) fn mixture_pdf(
        &self,
        triangles: &[AdmittedThermalTriangle],
        primitive_index: usize,
        origin: Point3,
        point: Point3,
        tree_probability: f64,
    ) -> f64 {
        self.point_pdf(
            0,
            triangles,
            primitive_index,
            origin,
            point,
            tree_probability,
        )
    }

    fn point_pdf(
        &self,
        index: usize,
        triangles: &[AdmittedThermalTriangle],
        primitive_index: usize,
        origin: Point3,
        point: Point3,
        probability: f64,
    ) -> f64 {
        let node = self.nodes[index];
        if probability == 0.0
            || primitive_index < node.primitive_min
            || primitive_index > node.primitive_max
            || !node.bounds.contains(point)
        {
            return 0.0;
        }
        match node.kind {
            Kind::Leaf(triangle) => {
                probability * triangles[triangle].directional_pdf(origin, point)
            }
            Kind::Branch { left, right } => {
                let Some(probabilities) = self.child_probabilities(triangles, origin, left, right)
                else {
                    return 0.0;
                };
                self.point_pdf(
                    left,
                    triangles,
                    primitive_index,
                    origin,
                    point,
                    probability * probabilities[0],
                ) + self.point_pdf(
                    right,
                    triangles,
                    primitive_index,
                    origin,
                    point,
                    probability * probabilities[1],
                )
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lighting::{
        AdmittedLighting, EnvironmentMap, LightKind, LightSample, RectLight, ThermalTriangleLight,
    };
    use crate::spectral::lift_rgb;
    use crate::thermal_surface::GrayDiffuseThermalSurface;
    use fs_blake3::ContentHash;
    use fs_geom::Vec3;

    fn source(index: usize) -> ThermalTriangleLight {
        let x = (index % 8) as f64 - 4.0;
        let z = (index / 8) as f64 - 4.0;
        ThermalTriangleLight {
            vertices: [
                Point3::new(x, 2.0, z),
                Point3::new(x + 0.6, 2.0, z),
                Point3::new(x, 2.0, z + 0.6),
            ],
            primitive_index: index,
            surface: GrayDiffuseThermalSurface::try_new(
                2400.0 + (index % 7) as f64 * 50.0,
                0.7,
                ContentHash([73; 32]),
            )
            .unwrap(),
        }
    }

    fn admitted(sources: Vec<ThermalTriangleLight>) -> AdmittedLighting<'static> {
        AdmittedLighting::try_new_controlled(&[], None, sources, false, || Ok(())).unwrap()
    }

    #[test]
    fn g0_tree_selection_frequencies_and_full_mixture_densities_agree() {
        let rectangles = [RectLight {
            corner: Point3::new(-0.5, 3.0, -0.5),
            edge_u: Vec3::new(1.0, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 0.0, 1.0),
            prim: 1000,
            emission: (lift_rgb([1.0; 3]), 40.0),
        }];
        let environment = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let lights = AdmittedLighting::try_new_controlled(
            &rectangles,
            Some(&environment),
            (0..64).map(source).collect(),
            false,
            || Ok(()),
        )
        .unwrap();
        assert!(lights.thermal_tree.is_some());
        let origin = Point3::new(0.25, 0.0, 0.25);
        let entries = lights.diagnostics(origin).entries;
        assert!((entries.iter().map(|entry| entry.probability).sum::<f64>() - 1.0).abs() < 1e-14);
        let mut counts = [0_usize; 66];
        const SAMPLES: usize = 16_384;
        for i in 0..SAMPLES {
            let sample = lights
                .sample(
                    origin,
                    (i as f64 + 0.5) / SAMPLES as f64,
                    ((i * 97) % SAMPLES) as f64 / SAMPLES as f64,
                )
                .unwrap();
            match sample {
                LightSample::ThermalSurface(sample) => {
                    counts[sample.primitive_index] += 1;
                    let reverse =
                        lights.thermal_mixture_pdf(sample.primitive_index, origin, sample.point);
                    assert_eq!(sample.pdf_solid_angle.to_bits(), reverse.to_bits());
                    let entry = entries
                        .iter()
                        .find(|entry| {
                            entry.kind == LightKind::ThermalSurface(sample.primitive_index)
                        })
                        .unwrap();
                    let delta = sample.point.delta_from(origin);
                    let distance = delta.norm();
                    let expected = entry.probability * distance * distance
                        / (0.18 * sample.normal.dot(delta.scale(1.0 / distance)).abs());
                    assert!((reverse / expected - 1.0).abs() < 3e-13);
                }
                LightSample::Rectangle(sample) => {
                    counts[64] += 1;
                    let reverse = lights.rect_mixture_pdf(sample.light_index, origin, sample.point);
                    assert!((sample.pdf_solid_angle / reverse - 1.0).abs() < 2e-13);
                }
                LightSample::Environment(sample) => {
                    counts[65] += 1;
                    let reverse = lights
                        .environment_evaluation(origin, sample.direction)
                        .unwrap();
                    assert_eq!(
                        sample.pdf_solid_angle.to_bits(),
                        reverse.pdf_solid_angle.to_bits()
                    );
                }
            }
        }
        for entry in entries {
            let index = match entry.kind {
                LightKind::ThermalSurface(index) => index,
                LightKind::Rectangle(_) => 64,
                LightKind::Environment => 65,
            };
            assert!((counts[index] as f64 - entry.probability * SAMPLES as f64).abs() < 2.1);
        }
    }

    fn linear_support_pdf(
        lights: &AdmittedLighting<'_>,
        origin: Point3,
        point: Point3,
        prim: usize,
    ) -> f64 {
        lights
            .diagnostics(origin)
            .entries
            .into_iter()
            .filter_map(|entry| {
                let LightKind::ThermalSurface(index) = entry.kind else {
                    return None;
                };
                let triangle = lights.thermal_triangles[index];
                (triangle.source.primitive_index == prim)
                    .then(|| entry.probability * triangle.directional_pdf(origin, point))
            })
            .sum()
    }

    #[test]
    fn g0_tree_reverse_support_includes_overlaps_edges_and_translated_bounds() {
        let mut sources: Vec<_> = (0..32).map(source).collect();
        sources.push(sources[0]);
        let mut inner = sources[0];
        for point in &mut inner.vertices[1..] {
            point.x = sources[0].vertices[0].x + (point.x - sources[0].vertices[0].x) * 0.5;
            point.z = sources[0].vertices[0].z + (point.z - sources[0].vertices[0].z) * 0.5;
        }
        sources.push(inner);
        for translation in [0.0, 1.0e10] {
            let translated = sources
                .iter()
                .map(|source| {
                    let mut result = *source;
                    for point in &mut result.vertices {
                        point.x += translation;
                        point.z += translation;
                    }
                    result
                })
                .collect();
            let lights = admitted(translated);
            let origin = Point3::new(translation, 0.0, translation);
            for triangle in &lights.thermal_triangles {
                for point in triangle.source.vertices.into_iter().chain([
                    triangle.sample_point(0.25, 0.5),
                    triangle.sample_point(0.5, 0.0),
                ]) {
                    let expected =
                        linear_support_pdf(&lights, origin, point, triangle.source.primitive_index);
                    let actual =
                        lights.thermal_mixture_pdf(triangle.source.primitive_index, origin, point);
                    assert!(expected > 0.0);
                    assert!((actual / expected - 1.0).abs() < 5e-13);
                }
            }
        }
    }

    #[test]
    fn g5_tree_stream_is_independent_of_input_order_primitive_indices_and_winding() {
        let sources: Vec<_> = (0..64).map(source).collect();
        let original = admitted(sources.clone());
        let reordered = admitted(
            sources
                .into_iter()
                .rev()
                .enumerate()
                .map(|(index, mut source)| {
                    source.primitive_index = index + 1000;
                    source.vertices.swap(0, 2);
                    source
                })
                .collect(),
        );
        let origin = Point3::new(0.0, -0.5, 0.0);
        for i in 0..512 {
            let u1 = (f64::from(i) + 0.5) / 512.0;
            let u2 = f64::from((i * 73) % 512) / 512.0;
            let LightSample::ThermalSurface(a) = original.sample(origin, u1, u2).unwrap() else {
                unreachable!();
            };
            let LightSample::ThermalSurface(b) = reordered.sample(origin, u1, u2).unwrap() else {
                unreachable!();
            };
            assert_eq!(
                [a.point.x, a.point.y, a.point.z],
                [b.point.x, b.point.y, b.point.z]
            );
            assert_eq!(a.pdf_solid_angle.to_bits(), b.pdf_solid_angle.to_bits());
            assert_eq!(a.surface, b.surface);
        }
    }

    #[test]
    fn g0_tree_preserves_finite_support_at_extreme_power_and_geometry_scales() {
        for (scale, luminance, distance) in [(1e80, 1e200, 1e120), (1e-80, 1e-200, 1e-70)] {
            let triangles: Vec<_> = (0..32)
                .map(|index| {
                    let mut source = source(index);
                    for point in &mut source.vertices {
                        point.x *= scale;
                        point.y *= scale;
                        point.z *= scale;
                    }
                    AdmittedThermalTriangle::try_new(source, luminance).unwrap()
                })
                .collect();
            let tree = ThermalLightTree::try_new(&triangles, &mut || Ok(()))
                .unwrap()
                .unwrap();
            let origin = Point3::new(0.0, -distance, 0.0);
            let weight = tree.weight(&triangles, origin);
            assert!(weight.is_finite() && weight > 0.0);
            let sum: f64 = (0..triangles.len())
                .map(|index| tree.selection_probability(&triangles, origin, index))
                .sum();
            assert!((sum - 1.0).abs() < 1e-14);
            for index in 0..triangles.len() {
                let probability = tree.selection_probability(&triangles, origin, index);
                assert!(probability.is_finite() && probability > 0.0);
                let point = triangles[index].sample_point(0.25, 0.5);
                let expected = probability * triangles[index].directional_pdf(origin, point);
                let actual = tree.mixture_pdf(&triangles, index, origin, point, 1.0);
                assert!(expected.is_finite() && expected > 0.0);
                assert!((actual / expected - 1.0).abs() < 1e-13);
            }
        }
    }

    #[test]
    fn g0_tree_depth_and_spatial_query_work_scale_with_hierarchy() {
        let triangles: Vec<_> = (0..4096)
            .map(|index| AdmittedThermalTriangle::try_new(source(index), 1.0).unwrap())
            .collect();
        let tree = ThermalLightTree::try_new(&triangles, &mut || Ok(()))
            .unwrap()
            .unwrap();
        fn depth(tree: &ThermalLightTree, node: usize) -> usize {
            match tree.nodes[node].kind {
                Kind::Leaf(_) => 0,
                Kind::Branch { left, right } => 1 + depth(tree, left).max(depth(tree, right)),
            }
        }
        fn query_nodes(tree: &ThermalLightTree, node: usize, point: Point3) -> usize {
            if !tree.nodes[node].bounds.contains(point) {
                return 1;
            }
            match tree.nodes[node].kind {
                Kind::Leaf(_) => 1,
                Kind::Branch { left, right } => {
                    1 + query_nodes(tree, left, point) + query_nodes(tree, right, point)
                }
            }
        }
        assert_eq!(tree.nodes.len(), 2 * triangles.len() - 1);
        assert_eq!(depth(&tree, 0), 12);
        for index in [0, 127, 1024, 2047, 4095] {
            let visited = query_nodes(&tree, 0, triangles[index].sample_point(0.25, 0.5));
            assert!(
                visited <= 1 + 4 * depth(&tree, 0),
                "visited {visited} of {} nodes",
                tree.nodes.len()
            );
        }
    }

    #[test]
    fn g4_tree_admission_polls_during_hierarchy_construction() {
        let triangles: Vec<_> = (0..64)
            .map(|index| AdmittedThermalTriangle::try_new(source(index), 1.0).unwrap())
            .collect();
        let cancel_at = 4 * triangles.len() + 17;
        let mut polls = 0;
        let result = ThermalLightTree::try_new(&triangles, &mut || {
            polls += 1;
            if polls == cancel_at {
                Err(LightingError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(result, Err(LightingError::Cancelled)));
        assert_eq!(polls, cancel_at);
    }
}
