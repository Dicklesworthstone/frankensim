//! Generic axisymmetric tessellation with evaluated binary64 bounds and conversion adapters
//! (bead frankensim-b8bxd.4).
//!
//! This module provides:
//! - [`AxisymmetricTessellationConfig`]: configuration specifying error budget, angular sector bounds, and purpose;
//! - [`AxisymmetricTessellationReceipt`]: receipt carrying evaluated analytic bounds, measured topology, and input provenance;
//! - [`AxisymmetricMesh`]: tessellated triangle mesh with vertex positions, normals, and per-triangle feature IDs;
//! - [`AxisymmetricRenderMesh`] and [`AxisymmetricCollisionMesh`]: distinct domain-typed wrappers;
//! - [`tessellate_axisymmetric`]: tessellation with analytic sagitta admission, not directed-rounding certification.

use crate::axisymmetric::{AxisymmetricChart, MeridianPoint, MeridianSegment};
use fs_evidence::ProvenanceHash;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};

/// Purpose of the tessellated mesh artifact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TessellationPurpose {
    /// Visual rendering mesh.
    Rendering,
    /// Collision and contact detection mesh.
    Collision,
}

/// Configuration for axisymmetric tessellation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisymmetricTessellationConfig {
    /// Maximum tolerated Hausdorff error between the true revolved surface and the triangle mesh.
    pub max_hausdorff_error: f64,
    /// Minimum azimuthal sectors around the axis of revolution.
    pub min_azimuthal_sectors: u32,
    /// Maximum azimuthal sectors around the axis of revolution.
    pub max_azimuthal_sectors: u32,
    /// Minimum subdivisions per circular arc feature.
    pub min_arc_subdivisions: u32,
    /// Maximum subdivisions per circular arc feature.
    pub max_arc_subdivisions: u32,
    /// Destination domain purpose.
    pub purpose: TessellationPurpose,
}

impl AxisymmetricTessellationConfig {
    /// Construct default configuration for the given error budget and purpose.
    ///
    /// # Errors
    /// [`AxisymmetricTessellationError::InvalidBudget`] if `max_hausdorff_error` is non-positive or non-finite.
    pub fn new(
        max_hausdorff_error: f64,
        purpose: TessellationPurpose,
    ) -> Result<Self, AxisymmetricTessellationError> {
        let config = Self {
            max_hausdorff_error,
            min_azimuthal_sectors: 8,
            max_azimuthal_sectors: 4096,
            min_arc_subdivisions: 2,
            max_arc_subdivisions: 1024,
            purpose,
        };
        config.validate()?;
        Ok(config)
    }

    fn validate(self) -> Result<(), AxisymmetricTessellationError> {
        if !self.max_hausdorff_error.is_finite() || self.max_hausdorff_error <= 0.0 {
            return Err(AxisymmetricTessellationError::InvalidBudget {
                error: self.max_hausdorff_error,
            });
        }
        for (field, min, max, supported_minimum) in [
            (
                "azimuthal sectors",
                self.min_azimuthal_sectors,
                self.max_azimuthal_sectors,
                3,
            ),
            (
                "arc subdivisions",
                self.min_arc_subdivisions,
                self.max_arc_subdivisions,
                1,
            ),
        ] {
            if min < supported_minimum || max < min {
                return Err(AxisymmetricTessellationError::InvalidResolutionLimits {
                    field,
                    min,
                    max,
                    supported_minimum,
                });
            }
        }
        Ok(())
    }
}

/// Evidence receipt for an axisymmetric tessellation.
#[derive(Debug, Clone, PartialEq)]
pub struct AxisymmetricTessellationReceipt {
    /// Evaluated analytic azimuthal sagitta, not a directed-rounding enclosure.
    pub azimuthal_sagitta_bound: f64,
    /// Evaluated analytic meridional sagitta, not a directed-rounding enclosure.
    pub meridian_sagitta_bound: f64,
    /// Combined evaluated analytic bound; binary64 evaluation is uncertified.
    pub total_hausdorff_bound: f64,
    /// Maximum angular error of vertex normals within smooth features (radians).
    pub max_smooth_normal_angular_error: f64,
    /// Total vertex count.
    pub vertex_count: usize,
    /// Total triangle count.
    pub triangle_count: usize,
    /// Topological Euler characteristic `V - E + F`.
    pub euler_characteristic: i32,
    /// Whether the mesh is verified closed and watertight.
    pub is_watertight: bool,
    /// Whether edge consistency and analytic-feature face-normal checks pass.
    /// This is a binary64 observation, not a topology/orientation certificate.
    pub is_outward_oriented: bool,
    /// Tessellation purpose.
    pub purpose: TessellationPurpose,
    /// Provenance hash over the input chart and tessellation configuration.
    pub provenance: ProvenanceHash,
}

/// A tessellated triangle mesh produced from an axisymmetric chart.
#[derive(Debug, Clone)]
pub struct AxisymmetricMesh {
    /// 3D vertex positions.
    pub positions: Vec<Point3>,
    /// 3D outward unit vertex normals.
    pub normals: Vec<Vec3>,
    /// Triangle face connectivity (index triples into `positions`).
    pub triangles: Vec<[u32; 3]>,
    /// Source meridian feature index for each triangle.
    pub triangle_features: Vec<usize>,
    /// Evidence receipt for the tessellation.
    pub receipt: AxisymmetricTessellationReceipt,
}

/// Strongly-typed rendering mesh wrapper.
#[derive(Debug, Clone)]
pub struct AxisymmetricRenderMesh(pub AxisymmetricMesh);

/// Strongly-typed collision mesh wrapper.
#[derive(Debug, Clone)]
pub struct AxisymmetricCollisionMesh(pub AxisymmetricMesh);

/// Errors during axisymmetric tessellation.
#[derive(Debug, Clone, PartialEq)]
pub enum AxisymmetricTessellationError {
    /// Cooperative cancellation was observed.
    Cancelled,
    /// The requested error budget was non-positive or non-finite.
    InvalidBudget {
        /// Rejected budget value.
        error: f64,
    },
    /// A subdivision range is empty, inverted, or below its geometric minimum.
    InvalidResolutionLimits {
        /// The rejected subdivision range.
        field: &'static str,
        /// Requested minimum.
        min: u32,
        /// Requested maximum.
        max: u32,
        /// Smallest supported minimum for this range.
        supported_minimum: u32,
    },
    /// The budget could not be satisfied within the maximum sector/subdivision caps.
    BudgetInfeasible {
        /// Requested Hausdorff error.
        requested: f64,
        /// Evaluated bound of the capped discretization rejected by admission.
        min_achievable: f64,
    },
    /// Non-finite coordinate encountered during discretization.
    NonFiniteGeometry,
    /// The chart has no segments.
    EmptyChart,
}

impl core::fmt::Display for AxisymmetricTessellationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled => write!(f, "axisymmetric tessellation cancelled"),
            Self::InvalidBudget { error } => write!(
                f,
                "invalid tessellation budget {error}: must be strictly positive and finite"
            ),
            Self::InvalidResolutionLimits {
                field,
                min,
                max,
                supported_minimum,
            } => write!(
                f,
                "invalid {field} limits {min}..={max}: require {supported_minimum} <= min <= max"
            ),
            Self::BudgetInfeasible {
                requested,
                min_achievable,
            } => write!(
                f,
                "tessellation budget {requested} infeasible: evaluated capped bound is {min_achievable}"
            ),
            Self::NonFiniteGeometry => write!(f, "non-finite geometry in axisymmetric chart"),
            Self::EmptyChart => write!(f, "axisymmetric chart has no segments"),
        }
    }
}

impl core::error::Error for AxisymmetricTessellationError {}

/// Polyline vertex on the 2D meridian.
#[derive(Debug, Clone, Copy)]
struct ProfileVertex {
    radius: f64,
    axial: f64,
    normal_r: f64,
    normal_z: f64,
    feature_index: usize,
}

/// Tessellate a validated [`AxisymmetricChart`] with evaluated bounds and topology checks.
pub fn tessellate_axisymmetric(
    chart: &AxisymmetricChart,
    config: AxisymmetricTessellationConfig,
    cx: &Cx<'_>,
) -> Result<AxisymmetricMesh, AxisymmetricTessellationError> {
    config.validate()?;

    let segments = chart.segments();
    if segments.is_empty() {
        return Err(AxisymmetricTessellationError::EmptyChart);
    }

    // Step 1: Find maximum radius across all segments
    let max_radius = maximum_profile_radius(segments);
    if max_radius <= 0.0 || !max_radius.is_finite() {
        return Err(AxisymmetricTessellationError::NonFiniteGeometry);
    }

    // Step 2: Compute azimuthal sectors required to bound azimuthal sagitta
    // Sagitta s_theta = R_max * (1 - cos(pi / N_theta))
    // We target s_theta <= config.max_hausdorff_error / 2
    let target_sagitta_theta = config.max_hausdorff_error * 0.5;
    let cos_val = (1.0 - target_sagitta_theta / max_radius).clamp(-1.0, 1.0);
    let min_theta_angle = cos_val.acos();
    let needed_sectors = if min_theta_angle > 1e-12 {
        (core::f64::consts::PI / min_theta_angle).ceil() as u32
    } else {
        config.max_azimuthal_sectors
    };

    let n_theta = needed_sectors
        .max(config.min_azimuthal_sectors)
        .min(config.max_azimuthal_sectors);

    let actual_sagitta_theta =
        max_radius * (1.0 - (core::f64::consts::PI / f64::from(n_theta)).cos());

    // Step 3: Compute arc subdivisions for each arc feature
    let target_sagitta_arc = config.max_hausdorff_error - actual_sagitta_theta;
    if actual_sagitta_theta > config.max_hausdorff_error {
        return Err(AxisymmetricTessellationError::BudgetInfeasible {
            requested: config.max_hausdorff_error,
            min_achievable: actual_sagitta_theta,
        });
    }

    let (profile_vertices, max_meridian_sagitta, max_normal_angular_error) =
        discretize_meridian(segments, config, n_theta, target_sagitta_arc, cx)?;
    // Profile and revolution errors can align radially. Compose by the
    // triangle inequality, rather than assuming orthogonal displacements.
    let total_hausdorff_bound = actual_sagitta_theta + max_meridian_sagitta;
    if !total_hausdorff_bound.is_finite() {
        return Err(AxisymmetricTessellationError::NonFiniteGeometry);
    }
    if total_hausdorff_bound > config.max_hausdorff_error {
        return Err(AxisymmetricTessellationError::BudgetInfeasible {
            requested: config.max_hausdorff_error,
            min_achievable: total_hausdorff_bound,
        });
    }

    // Step 4: Revolve profile vertices into 3D positions and normals
    let (positions, normals, vertex_grid) = revolve_profile(&profile_vertices, n_theta);

    // Step 5: Construct triangles
    let (triangles, triangle_features) =
        triangulate_profile(&profile_vertices, &vertex_grid, n_theta);

    // Step 6: Verify topology (Euler characteristic and watertightness)
    let v_count = positions.len();
    let f_count = triangles.len();
    let (euler_characteristic, is_watertight, is_consistently_oriented) =
        mesh_topology(v_count, &triangles);
    let is_outward_oriented = is_consistently_oriented
        && all_faces_outward(chart, &positions, &triangles, &triangle_features);
    let provenance = tessellation_provenance(chart, config);

    let receipt = AxisymmetricTessellationReceipt {
        azimuthal_sagitta_bound: actual_sagitta_theta,
        meridian_sagitta_bound: max_meridian_sagitta,
        total_hausdorff_bound,
        max_smooth_normal_angular_error: max_normal_angular_error,
        vertex_count: v_count,
        triangle_count: f_count,
        euler_characteristic,
        is_watertight,
        is_outward_oriented,
        purpose: config.purpose,
        provenance,
    };

    Ok(AxisymmetricMesh {
        positions,
        normals,
        triangles,
        triangle_features,
        receipt,
    })
}

fn all_faces_outward(
    chart: &AxisymmetricChart,
    positions: &[Point3],
    triangles: &[[u32; 3]],
    features: &[usize],
) -> bool {
    !triangles.is_empty()
        && triangles.len() == features.len()
        && triangles.iter().zip(features).all(|(triangle, &feature)| {
            let [a, b, c] = triangle.map(|index| positions[index as usize]);
            let ab = Vec3::new(b.x - a.x, b.y - a.y, b.z - a.z);
            let ac = Vec3::new(c.x - a.x, c.y - a.y, c.z - a.z);
            let face = Vec3::new(
                ab.y * ac.z - ab.z * ac.y,
                ab.z * ac.x - ab.x * ac.z,
                ab.x * ac.y - ab.y * ac.x,
            );
            let center = Point3::new(
                (a.x + b.x + c.x) / 3.0,
                (a.y + b.y + c.y) / 3.0,
                (a.z + b.z + c.z) / 3.0,
            );
            // Average the generating meridian radii, not the Cartesian
            // centroid radius: the latter contracts across azimuthal sectors.
            let meridian_radius = (a.x.hypot(a.y) + b.x.hypot(b.y) + c.x.hypot(c.y)) / 3.0;
            let Some(segment) = chart.segments().get(feature).copied() else {
                return false;
            };
            let Some(normal) = feature_outward_normal(segment, center, meridian_radius) else {
                return false;
            };
            let alignment = face.dot(normal);
            alignment.is_finite() && alignment > 0.0
        })
}

fn feature_outward_normal(
    segment: MeridianSegment,
    center: Point3,
    meridian_radius: f64,
) -> Option<Vec3> {
    let radius = center.x.hypot(center.y);
    if !radius.is_finite() || radius <= 0.0 {
        return None;
    }
    let (radial, axial) = match segment {
        MeridianSegment::Line { start, end } => {
            (end.axial - start.axial, start.radius - end.radius)
        }
        MeridianSegment::Arc {
            center: origin,
            clockwise,
            ..
        } => {
            let sign = if clockwise { -1.0 } else { 1.0 };
            (
                sign * (meridian_radius - origin.radius),
                sign * (center.z - origin.axial),
            )
        }
    };
    Some(Vec3::new(
        radial * center.x / radius,
        radial * center.y / radius,
        axial,
    ))
}

fn tessellation_provenance(
    chart: &AxisymmetricChart,
    config: AxisymmetricTessellationConfig,
) -> ProvenanceHash {
    let mut bytes = b"fs-rep-frep/axisymmetric/tessellation/v1".to_vec();
    bytes.extend_from_slice(&crate::axisymmetric::profile_encoding(chart.segments()));
    bytes.extend_from_slice(&config.max_hausdorff_error.to_bits().to_le_bytes());
    for count in [
        config.min_azimuthal_sectors,
        config.max_azimuthal_sectors,
        config.min_arc_subdivisions,
        config.max_arc_subdivisions,
    ] {
        bytes.extend_from_slice(&count.to_le_bytes());
    }
    bytes.push(match config.purpose {
        TessellationPurpose::Rendering => 0,
        TessellationPurpose::Collision => 1,
    });
    ProvenanceHash::of_bytes(&bytes)
}

fn maximum_profile_radius(segments: &[MeridianSegment]) -> f64 {
    let mut max_radius = 0.0f64;
    for seg in segments {
        match *seg {
            MeridianSegment::Line { start, end } => {
                max_radius = max_radius.max(start.radius).max(end.radius);
            }
            MeridianSegment::Arc {
                start, end, center, ..
            } => {
                let arc_r = (start.radius - center.radius).hypot(start.axial - center.axial);
                max_radius = max_radius
                    .max(start.radius)
                    .max(end.radius)
                    .max(center.radius + arc_r);
            }
        }
    }
    max_radius
}

fn discretize_meridian(
    segments: &[MeridianSegment],
    config: AxisymmetricTessellationConfig,
    n_theta: u32,
    target_sagitta_arc: f64,
    cx: &Cx<'_>,
) -> Result<(Vec<ProfileVertex>, f64, f64), AxisymmetricTessellationError> {
    let mut profile_vertices = Vec::new();
    let mut max_meridian_sagitta = 0.0f64;
    let mut max_normal_angular_error = core::f64::consts::PI / f64::from(n_theta);
    for (feature_idx, seg) in segments.iter().copied().enumerate() {
        if cx.checkpoint().is_err() {
            return Err(AxisymmetricTessellationError::Cancelled);
        }
        match seg {
            MeridianSegment::Line { start, end } => {
                append_line_vertices(start, end, feature_idx, &mut profile_vertices);
            }
            MeridianSegment::Arc { .. } => {
                append_arc_vertices(
                    seg,
                    config,
                    target_sagitta_arc,
                    feature_idx,
                    &mut profile_vertices,
                    &mut max_meridian_sagitta,
                    &mut max_normal_angular_error,
                );
            }
        }
    }
    if profile_vertices.len() < 2 {
        return Err(AxisymmetricTessellationError::NonFiniteGeometry);
    }
    Ok((
        profile_vertices,
        max_meridian_sagitta,
        max_normal_angular_error,
    ))
}

fn append_line_vertices(
    start: MeridianPoint,
    end: MeridianPoint,
    feature_idx: usize,
    profile_vertices: &mut Vec<ProfileVertex>,
) {
    // Skip internal axis lines (r == 0).
    if start.radius < 1e-12 && end.radius < 1e-12 {
        return;
    }
    let dr = end.radius - start.radius;
    let dz = end.axial - start.axial;
    let len = dr.hypot(dz);
    if len <= 0.0 {
        return;
    }
    let nr = dz / len;
    let nz = -dr / len;
    if profile_vertices.is_empty()
        || (profile_vertices.last().unwrap().radius - start.radius).abs() > 1e-9
        || (profile_vertices.last().unwrap().axial - start.axial).abs() > 1e-9
    {
        profile_vertices.push(ProfileVertex {
            radius: start.radius,
            axial: start.axial,
            normal_r: nr,
            normal_z: nz,
            feature_index: feature_idx,
        });
    }
    profile_vertices.push(ProfileVertex {
        radius: end.radius,
        axial: end.axial,
        normal_r: nr,
        normal_z: nz,
        feature_index: feature_idx,
    });
}

fn append_arc_vertices(
    segment: MeridianSegment,
    config: AxisymmetricTessellationConfig,
    target_sagitta_arc: f64,
    feature_idx: usize,
    profile_vertices: &mut Vec<ProfileVertex>,
    max_meridian_sagitta: &mut f64,
    max_normal_angular_error: &mut f64,
) {
    let MeridianSegment::Arc {
        start,
        end,
        center,
        clockwise,
    } = segment
    else {
        unreachable!();
    };
    let arc_r = (start.radius - center.radius).hypot(start.axial - center.axial);
    let angle_start = (start.axial - center.axial).atan2(start.radius - center.radius);
    let angle_end = (end.axial - center.axial).atan2(end.radius - center.radius);
    let sweep = if clockwise {
        if angle_end <= angle_start {
            angle_start - angle_end
        } else {
            angle_start - angle_end + 2.0 * core::f64::consts::PI
        }
    } else {
        if angle_end >= angle_start {
            angle_end - angle_start
        } else {
            angle_end - angle_start + 2.0 * core::f64::consts::PI
        }
    };
    let cos_arc = (1.0 - target_sagitta_arc / arc_r).clamp(-1.0, 1.0);
    let max_arc_subdiv_angle = 2.0 * cos_arc.acos();
    let needed_arc_subdivs = if max_arc_subdiv_angle > 1e-12 {
        (sweep / max_arc_subdiv_angle).ceil() as u32
    } else {
        config.max_arc_subdivisions
    };
    let n_arc = needed_arc_subdivs
        .max(config.min_arc_subdivisions)
        .min(config.max_arc_subdivisions);
    let actual_arc_sagitta = arc_r * (1.0 - (sweep / (2.0 * f64::from(n_arc))).cos());
    *max_meridian_sagitta = max_meridian_sagitta.max(actual_arc_sagitta);
    *max_normal_angular_error = max_normal_angular_error.max(sweep / f64::from(n_arc));
    let d_angle = if clockwise {
        -sweep / f64::from(n_arc)
    } else {
        sweep / f64::from(n_arc)
    };
    for k in 0..=n_arc {
        let angle = angle_start + f64::from(k) * d_angle;
        let r = center.radius + arc_r * angle.cos();
        let z = center.axial + arc_r * angle.sin();
        let (nr, nz) = if clockwise {
            (-angle.cos(), -angle.sin())
        } else {
            (angle.cos(), angle.sin())
        };
        if k == 0 {
            if profile_vertices.is_empty()
                || (profile_vertices.last().unwrap().radius - r).abs() > 1e-9
                || (profile_vertices.last().unwrap().axial - z).abs() > 1e-9
            {
                profile_vertices.push(ProfileVertex {
                    radius: r,
                    axial: z,
                    normal_r: nr,
                    normal_z: nz,
                    feature_index: feature_idx,
                });
            }
        } else {
            profile_vertices.push(ProfileVertex {
                radius: r,
                axial: z,
                normal_r: nr,
                normal_z: nz,
                feature_index: feature_idx,
            });
        }
    }
}

fn revolve_profile(
    profile_vertices: &[ProfileVertex],
    n_theta: u32,
) -> (Vec<Point3>, Vec<Vec3>, Vec<Vec<u32>>) {
    let m = profile_vertices.len();
    let mut positions: Vec<Point3> = Vec::with_capacity(m * n_theta as usize);
    let mut normals: Vec<Vec3> = Vec::with_capacity(m * n_theta as usize);
    let mut vertex_grid: Vec<Vec<u32>> = vec![vec![0; n_theta as usize]; m];

    for (i, pv) in profile_vertices.iter().enumerate() {
        if pv.radius < 1e-12 {
            // Vertex on the axis: single vertex shared across all sectors
            let idx = positions.len() as u32;
            positions.push(Point3::new(0.0, 0.0, pv.axial));
            let nz = if pv.normal_z >= 0.0 { 1.0 } else { -1.0 };
            normals.push(Vec3::new(0.0, 0.0, nz));
            vertex_grid[i].fill(idx);
        } else {
            for (j, index) in vertex_grid[i].iter_mut().enumerate() {
                let theta = 2.0 * core::f64::consts::PI * j as f64 / f64::from(n_theta);
                let cos_t = theta.cos();
                let sin_t = theta.sin();
                let idx = positions.len() as u32;
                positions.push(Point3::new(pv.radius * cos_t, pv.radius * sin_t, pv.axial));
                normals.push(Vec3::new(
                    pv.normal_r * cos_t,
                    pv.normal_r * sin_t,
                    pv.normal_z,
                ));
                *index = idx;
            }
        }
    }

    (positions, normals, vertex_grid)
}

fn triangulate_profile(
    profile_vertices: &[ProfileVertex],
    vertex_grid: &[Vec<u32>],
    n_theta: u32,
) -> (Vec<[u32; 3]>, Vec<usize>) {
    let m = profile_vertices.len();
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    let mut triangle_features: Vec<usize> = Vec::new();

    for i in 0..m - 1 {
        // Each interval is generated by the feature that appended its end.
        let feat = profile_vertices[i + 1].feature_index;
        let r0 = profile_vertices[i].radius;
        let r1 = profile_vertices[i + 1].radius;

        for j in 0..n_theta as usize {
            let j_next = (j + 1) % n_theta as usize;
            let v00 = vertex_grid[i][j];
            let v10 = vertex_grid[i + 1][j];
            let v11 = vertex_grid[i + 1][j_next];
            let v01 = vertex_grid[i][j_next];

            if r0 < 1e-12 {
                // Top/bottom cap vertex fan (degenerate top edge)
                if v00 != v10 && v10 != v11 && v11 != v00 {
                    triangles.push([v00, v11, v10]);
                    triangle_features.push(feat);
                }
            } else if r1 < 1e-12 {
                // Top/bottom cap vertex fan (degenerate bottom edge)
                if v00 != v10 && v10 != v01 && v01 != v00 {
                    triangles.push([v00, v01, v10]);
                    triangle_features.push(feat);
                }
            } else {
                // Regular quad split into two oriented triangles
                triangles.push([v00, v11, v10]);
                triangle_features.push(feat);
                triangles.push([v00, v01, v11]);
                triangle_features.push(feat);
            }
        }
    }

    (triangles, triangle_features)
}

fn mesh_topology(v_count: usize, triangles: &[[u32; 3]]) -> (i32, bool, bool) {
    let f_count = triangles.len();

    let mut edges: std::collections::BTreeMap<[u32; 2], (u32, i32)> =
        std::collections::BTreeMap::new();
    for tri in triangles {
        for c in 0..3 {
            let (a, b) = (tri[c], tri[(c + 1) % 3]);
            if a == b {
                continue;
            }
            let key = if a < b { [a, b] } else { [b, a] };
            let entry = edges.entry(key).or_insert((0, 0));
            entry.0 += 1;
            entry.1 += if a < b { 1 } else { -1 };
        }
    }

    let e_count = edges.len();
    let mut is_watertight = true;
    let mut is_consistently_oriented = true;

    for (&_key, &(count, orient)) in &edges {
        if count != 2 {
            is_watertight = false;
        }
        if orient != 0 {
            is_consistently_oriented = false;
        }
    }

    let euler_characteristic = v_count as i32 - e_count as i32 + f_count as i32;

    (
        euler_characteristic,
        is_watertight,
        is_consistently_oriented,
    )
}
