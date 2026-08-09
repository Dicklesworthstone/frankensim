//! Deterministic conforming tetrahedral meshes for circularly filleted solid
//! cylinders.
//!
//! This is a reusable geometry primitive, not an Euler-disc special case. The
//! radial grid contains the exact cap/fillet tangent radius, the azimuthal grid
//! is periodic, and normalized axial layers follow the exact circular meridian
//! at every radial ring. The resulting boundary is a piecewise-planar
//! approximation whose chord error converges under radial/azimuthal refinement.

use std::collections::BTreeMap;

use fs_exec::Cx;

/// Parameterized geometry, resolution, and work envelope.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RoundedCylinderMeshSpec {
    /// Outer radius [m].
    pub outer_radius_m: f64,
    /// Total axial thickness [m].
    pub thickness_m: f64,
    /// Equal upper/lower circular fillet radius [m]. Zero means sharp rims.
    pub fillet_radius_m: f64,
    /// Radial intervals from the axis through the planar-cap region.
    pub core_radial_segments: u32,
    /// Radial intervals across the circular fillet. Must be zero for a sharp
    /// rim and positive for a nonzero fillet.
    pub fillet_radial_segments: u32,
    /// Periodic angular intervals around the axis.
    pub azimuthal_segments: u32,
    /// Axial intervals at every radial ring.
    pub axial_segments: u32,
    /// Maximum admitted vertices.
    pub maximum_vertices: usize,
    /// Maximum admitted tetrahedra.
    pub maximum_tetrahedra: usize,
}

impl RoundedCylinderMeshSpec {
    /// A bounded modal-analysis mesh resolution.
    #[must_use]
    pub const fn modal_default(
        outer_radius_m: f64,
        thickness_m: f64,
        fillet_radius_m: f64,
    ) -> Self {
        Self {
            outer_radius_m,
            thickness_m,
            fillet_radius_m,
            core_radial_segments: 6,
            fillet_radial_segments: if fillet_radius_m > 0.0 { 3 } else { 0 },
            azimuthal_segments: 24,
            axial_segments: 2,
            maximum_vertices: 100_000,
            maximum_tetrahedra: 600_000,
        }
    }
}

/// A boundary triangle surface with geometry needed by BEM and mode sampling.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundaryPanelMesh {
    /// Outward-oriented triangle vertex indices.
    pub triangles: Vec<[usize; 3]>,
    /// Triangle centroids [m].
    pub centroids_m: Vec<[f64; 3]>,
    /// Outward unit normals.
    pub normals: Vec<[f64; 3]>,
    /// Triangle areas [m^2].
    pub areas_m2: Vec<f64>,
}

/// Conforming volume mesh plus its derived exterior boundary.
#[derive(Debug, Clone, PartialEq)]
pub struct RoundedCylinderTetMesh {
    /// Vertex coordinates [m], centered at the geometric midplane/origin.
    pub nodes_m: Vec<[f64; 3]>,
    /// Conforming tetrahedral connectivity.
    pub tetrahedra: Vec<[usize; 4]>,
    /// Derived closed exterior triangle mesh.
    pub boundary: BoundaryPanelMesh,
    /// Maximum meridian circular-arc chord error [m].
    pub maximum_meridian_chord_error_m: f64,
    /// Maximum azimuthal chord error at the outer radius [m].
    pub maximum_azimuthal_chord_error_m: f64,
}

/// Typed refusal from rounded-cylinder volume meshing.
#[derive(Debug, Clone, PartialEq)]
pub enum RoundedCylinderMeshError {
    /// A named geometry scalar is non-finite or outside its admissible range.
    InvalidGeometry {
        /// Failed scalar/domain description.
        what: &'static str,
    },
    /// Resolution cannot produce a conforming periodic mesh.
    InvalidResolution {
        /// Failed resolution invariant.
        what: &'static str,
    },
    /// Checked size arithmetic overflowed.
    SizeOverflow,
    /// A declared count envelope was exceeded.
    BudgetExceeded {
        /// Bounded resource.
        what: &'static str,
        /// Requested count.
        requested: usize,
        /// Admitted count.
        maximum: usize,
    },
    /// Derived topology was non-manifold or internally inconsistent.
    Topology {
        /// Failed invariant.
        what: &'static str,
    },
    /// Cancellation was observed before publication.
    Cancelled,
}

impl core::fmt::Display for RoundedCylinderMeshError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidGeometry { what } => {
                write!(f, "FS-MESH-ROUNDED-CYLINDER-GEOMETRY: {what}")
            }
            Self::InvalidResolution { what } => {
                write!(f, "FS-MESH-ROUNDED-CYLINDER-RESOLUTION: {what}")
            }
            Self::SizeOverflow => write!(f, "FS-MESH-ROUNDED-CYLINDER-SIZE-OVERFLOW"),
            Self::BudgetExceeded {
                what,
                requested,
                maximum,
            } => write!(
                f,
                "FS-MESH-ROUNDED-CYLINDER-BUDGET: {what} requested {requested}, maximum {maximum}"
            ),
            Self::Topology { what } => {
                write!(f, "FS-MESH-ROUNDED-CYLINDER-TOPOLOGY: {what}")
            }
            Self::Cancelled => write!(f, "FS-MESH-ROUNDED-CYLINDER-CANCELLED"),
        }
    }
}

impl std::error::Error for RoundedCylinderMeshError {}

/// Build a deterministic body-fitted tetrahedral mesh.
///
/// # Errors
/// Returns [`RoundedCylinderMeshError`] for invalid geometry/resolution,
/// exceeded work bounds, topology failure, or cancellation.
pub fn rounded_cylinder_tet_mesh(
    spec: RoundedCylinderMeshSpec,
    cx: &Cx<'_>,
) -> Result<RoundedCylinderTetMesh, RoundedCylinderMeshError> {
    cx.checkpoint()
        .map_err(|_| RoundedCylinderMeshError::Cancelled)?;
    validate(spec)?;
    let radial = radial_coordinates(spec);
    let radial_count = radial.len() - 1;
    let azimuthal = usize::try_from(spec.azimuthal_segments)
        .map_err(|_| RoundedCylinderMeshError::SizeOverflow)?;
    let axial =
        usize::try_from(spec.axial_segments).map_err(|_| RoundedCylinderMeshError::SizeOverflow)?;
    let levels = axial
        .checked_add(1)
        .ok_or(RoundedCylinderMeshError::SizeOverflow)?;
    let nodes_per_level = radial_count
        .checked_mul(azimuthal)
        .and_then(|value| value.checked_add(1))
        .ok_or(RoundedCylinderMeshError::SizeOverflow)?;
    let vertex_count = nodes_per_level
        .checked_mul(levels)
        .ok_or(RoundedCylinderMeshError::SizeOverflow)?;
    let tets_per_sector_slab = radial_count
        .checked_sub(1)
        .and_then(|annuli| annuli.checked_mul(6))
        .and_then(|annulus_tets| annulus_tets.checked_add(3))
        .ok_or(RoundedCylinderMeshError::SizeOverflow)?;
    let tetrahedron_count = tets_per_sector_slab
        .checked_mul(azimuthal)
        .and_then(|value| value.checked_mul(axial))
        .ok_or(RoundedCylinderMeshError::SizeOverflow)?;
    check_budget("vertices", vertex_count, spec.maximum_vertices)?;
    check_budget("tetrahedra", tetrahedron_count, spec.maximum_tetrahedra)?;

    let mut nodes_m = Vec::with_capacity(vertex_count);
    for level in 0..levels {
        cx.checkpoint()
            .map_err(|_| RoundedCylinderMeshError::Cancelled)?;
        let eta = -1.0 + 2.0 * level as f64 / axial as f64;
        nodes_m.push([0.0, 0.0, eta * 0.5 * spec.thickness_m]);
        for &radius in radial.iter().skip(1) {
            let z = eta * half_height(spec, radius);
            for sector in 0..azimuthal {
                let theta = core::f64::consts::TAU * sector as f64 / azimuthal as f64;
                nodes_m.push([radius * theta.cos(), radius * theta.sin(), z]);
            }
        }
    }

    let node = |level: usize, ring: usize, sector: usize| -> usize {
        let offset = level * nodes_per_level;
        if ring == 0 {
            offset
        } else {
            offset + 1 + (ring - 1) * azimuthal + sector % azimuthal
        }
    };
    let mut tetrahedra = Vec::with_capacity(tetrahedron_count);
    for level in 0..axial {
        cx.checkpoint()
            .map_err(|_| RoundedCylinderMeshError::Cancelled)?;
        for sector in 0..azimuthal {
            let next = (sector + 1) % azimuthal;
            // Central triangular prism. Its outer rectangular face uses the
            // same lower-current to upper-next diagonal as the first annulus.
            let a = node(level, 0, sector);
            let b = node(level, 1, sector);
            let c = node(level, 1, next);
            let aa = node(level + 1, 0, sector);
            let bb = node(level + 1, 1, sector);
            let cc = node(level + 1, 1, next);
            tetrahedra.extend([[a, b, c, cc], [a, b, cc, bb], [a, aa, bb, cc]]);

            for inner_ring in 1..radial_count {
                let outer_ring = inner_ring + 1;
                let v000 = node(level, inner_ring, sector);
                let v100 = node(level, outer_ring, sector);
                let v110 = node(level, outer_ring, next);
                let v010 = node(level, inner_ring, next);
                let v001 = node(level + 1, inner_ring, sector);
                let v101 = node(level + 1, outer_ring, sector);
                let v111 = node(level + 1, outer_ring, next);
                let v011 = node(level + 1, inner_ring, next);
                tetrahedra.extend([
                    [v000, v100, v110, v111],
                    [v000, v110, v010, v111],
                    [v000, v010, v011, v111],
                    [v000, v011, v001, v111],
                    [v000, v001, v101, v111],
                    [v000, v101, v100, v111],
                ]);
            }
        }
    }
    if tetrahedra.len() != tetrahedron_count || nodes_m.len() != vertex_count {
        return Err(RoundedCylinderMeshError::Topology {
            what: "constructed counts disagree with preflight",
        });
    }
    let boundary = extract_boundary(&nodes_m, &tetrahedra, cx)?;
    let maximum_meridian_chord_error_m = if spec.fillet_radius_m == 0.0 {
        0.0
    } else {
        let angle = 0.5 * core::f64::consts::FRAC_PI_2 / f64::from(spec.fillet_radial_segments);
        spec.fillet_radius_m * (1.0 - angle.cos())
    };
    let maximum_azimuthal_chord_error_m = spec.outer_radius_m
        * (1.0 - (core::f64::consts::PI / f64::from(spec.azimuthal_segments)).cos());
    cx.checkpoint()
        .map_err(|_| RoundedCylinderMeshError::Cancelled)?;
    Ok(RoundedCylinderTetMesh {
        nodes_m,
        tetrahedra,
        boundary,
        maximum_meridian_chord_error_m,
        maximum_azimuthal_chord_error_m,
    })
}

fn validate(spec: RoundedCylinderMeshSpec) -> Result<(), RoundedCylinderMeshError> {
    if !(spec.outer_radius_m.is_finite() && spec.outer_radius_m > 0.0) {
        return Err(RoundedCylinderMeshError::InvalidGeometry {
            what: "outer_radius_m must be finite and positive",
        });
    }
    if !(spec.thickness_m.is_finite() && spec.thickness_m > 0.0) {
        return Err(RoundedCylinderMeshError::InvalidGeometry {
            what: "thickness_m must be finite and positive",
        });
    }
    let half = 0.5 * spec.thickness_m;
    if !(spec.fillet_radius_m.is_finite()
        && spec.fillet_radius_m >= 0.0
        && spec.fillet_radius_m < spec.outer_radius_m
        && spec.fillet_radius_m < half)
    {
        return Err(RoundedCylinderMeshError::InvalidGeometry {
            what: "fillet_radius_m must be finite and in [0,min(radius,thickness/2))",
        });
    }
    if spec.core_radial_segments == 0 {
        return Err(RoundedCylinderMeshError::InvalidResolution {
            what: "core_radial_segments must be positive",
        });
    }
    if (spec.fillet_radius_m == 0.0) != (spec.fillet_radial_segments == 0) {
        return Err(RoundedCylinderMeshError::InvalidResolution {
            what: "fillet_radial_segments must be zero exactly for a sharp rim",
        });
    }
    if spec.azimuthal_segments < 3 {
        return Err(RoundedCylinderMeshError::InvalidResolution {
            what: "azimuthal_segments must be at least three",
        });
    }
    if spec.axial_segments == 0 {
        return Err(RoundedCylinderMeshError::InvalidResolution {
            what: "axial_segments must be positive",
        });
    }
    Ok(())
}

fn check_budget(
    what: &'static str,
    requested: usize,
    maximum: usize,
) -> Result<(), RoundedCylinderMeshError> {
    if requested > maximum {
        return Err(RoundedCylinderMeshError::BudgetExceeded {
            what,
            requested,
            maximum,
        });
    }
    Ok(())
}

fn radial_coordinates(spec: RoundedCylinderMeshSpec) -> Vec<f64> {
    let mut radial = Vec::with_capacity(
        1 + spec.core_radial_segments as usize + spec.fillet_radial_segments as usize,
    );
    radial.push(0.0);
    let tangent = spec.outer_radius_m - spec.fillet_radius_m;
    for i in 1..=spec.core_radial_segments {
        radial.push(tangent * f64::from(i) / f64::from(spec.core_radial_segments));
    }
    for i in 1..=spec.fillet_radial_segments {
        radial.push(
            tangent + spec.fillet_radius_m * f64::from(i) / f64::from(spec.fillet_radial_segments),
        );
    }
    radial
}

fn half_height(spec: RoundedCylinderMeshSpec, radius: f64) -> f64 {
    let half = 0.5 * spec.thickness_m;
    let fillet = spec.fillet_radius_m;
    let tangent = spec.outer_radius_m - fillet;
    if fillet == 0.0 || radius <= tangent {
        half
    } else {
        let dr = radius - tangent;
        half - fillet + (fillet * fillet - dr * dr).max(0.0).sqrt()
    }
}

fn extract_boundary(
    nodes: &[[f64; 3]],
    tetrahedra: &[[usize; 4]],
    cx: &Cx<'_>,
) -> Result<BoundaryPanelMesh, RoundedCylinderMeshError> {
    let mut faces: BTreeMap<[usize; 3], ([usize; 3], u8)> = BTreeMap::new();
    for (index, tet) in tetrahedra.iter().enumerate() {
        if index.is_multiple_of(256) {
            cx.checkpoint()
                .map_err(|_| RoundedCylinderMeshError::Cancelled)?;
        }
        for opposite in 0..4 {
            let mut face = [0; 3];
            let mut cursor = 0;
            for (corner, &vertex) in tet.iter().enumerate() {
                if corner != opposite {
                    face[cursor] = vertex;
                    cursor += 1;
                }
            }
            if points_toward(nodes, face, tet[opposite]) {
                face.swap(1, 2);
            }
            let mut key = face;
            key.sort_unstable();
            let entry = faces.entry(key).or_insert((face, 0));
            entry.1 = entry
                .1
                .checked_add(1)
                .ok_or(RoundedCylinderMeshError::Topology {
                    what: "face incidence overflow",
                })?;
            if entry.1 > 2 {
                return Err(RoundedCylinderMeshError::Topology {
                    what: "non-manifold face has more than two incident tetrahedra",
                });
            }
        }
    }
    let triangles: Vec<[usize; 3]> = faces
        .values()
        .filter_map(|(face, count)| (*count == 1).then_some(*face))
        .collect();
    if triangles.is_empty() || faces.values().any(|(_, count)| *count == 0 || *count > 2) {
        return Err(RoundedCylinderMeshError::Topology {
            what: "derived boundary is empty or has invalid incidence",
        });
    }
    let mut centroids_m = Vec::with_capacity(triangles.len());
    let mut normals = Vec::with_capacity(triangles.len());
    let mut areas_m2 = Vec::with_capacity(triangles.len());
    for face in &triangles {
        let [a, b, c] = face.map(|vertex| nodes[vertex]);
        let cross = cross(sub(b, a), sub(c, a));
        let double_area = norm(cross);
        if !(double_area.is_finite() && double_area > 0.0) {
            return Err(RoundedCylinderMeshError::Topology {
                what: "boundary contains a degenerate triangle",
            });
        }
        centroids_m.push([
            (a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
            (a[2] + b[2] + c[2]) / 3.0,
        ]);
        normals.push([
            cross[0] / double_area,
            cross[1] / double_area,
            cross[2] / double_area,
        ]);
        areas_m2.push(0.5 * double_area);
    }
    Ok(BoundaryPanelMesh {
        triangles,
        centroids_m,
        normals,
        areas_m2,
    })
}

fn points_toward(nodes: &[[f64; 3]], face: [usize; 3], opposite: usize) -> bool {
    let [a, b, c] = face.map(|vertex| nodes[vertex]);
    dot(cross(sub(b, a), sub(c, a)), sub(nodes[opposite], a)) > 0.0
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_rep_mesh::TetComplex;

    fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 11,
                    kernel_id: 19,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn spec(core: u32, fillet: u32, azimuthal: u32) -> RoundedCylinderMeshSpec {
        RoundedCylinderMeshSpec {
            outer_radius_m: 0.038,
            thickness_m: 0.006,
            fillet_radius_m: 0.001,
            core_radial_segments: core,
            fillet_radial_segments: fillet,
            azimuthal_segments: azimuthal,
            axial_segments: 2,
            maximum_vertices: 1_000_000,
            maximum_tetrahedra: 6_000_000,
        }
    }

    #[test]
    fn g0_mesh_is_closed_oriented_and_chain_exact() {
        let mesh = with_cx(|cx| rounded_cylinder_tet_mesh(spec(3, 2, 12), cx)).unwrap();
        assert!(!mesh.boundary.triangles.is_empty());
        for ((centroid, normal), area) in mesh
            .boundary
            .centroids_m
            .iter()
            .zip(&mesh.boundary.normals)
            .zip(&mesh.boundary.areas_m2)
        {
            assert!(*area > 0.0);
            assert!((norm(*normal) - 1.0).abs() < 1.0e-12);
            // The solid is convex and centered, so every outward face normal
            // has nonnegative support at its centroid.
            assert!(dot(*centroid, *normal) >= -1.0e-14);
        }
        let tets_u32 = mesh
            .tetrahedra
            .iter()
            .map(|tet| tet.map(|v| u32::try_from(v).unwrap()))
            .collect();
        let complex = TetComplex::from_tets(mesh.nodes_m.len(), tets_u32);
        let vertex_values: Vec<i64> = (0..complex.vertex_count)
            .map(|index| (index as i64).wrapping_mul(17).wrapping_sub(9))
            .collect();
        let d1d0 = complex.d1().apply(&complex.d0().apply(&vertex_values));
        assert!(d1d0.iter().all(|value| *value == 0));
        let edge_values: Vec<i64> = (0..complex.edges.len())
            .map(|index| (index as i64).wrapping_mul(31).wrapping_add(5))
            .collect();
        let d2d1 = complex.d2().apply(&complex.d1().apply(&edge_values));
        assert!(d2d1.iter().all(|value| *value == 0));
    }

    #[test]
    fn g1_volume_converges_to_exact_filleted_solid() {
        let volume = |mesh: &RoundedCylinderTetMesh| {
            mesh.tetrahedra
                .iter()
                .map(|tet| {
                    let [a, b, c, d] = tet.map(|v| mesh.nodes_m[v]);
                    determinant([sub(b, a), sub(c, a), sub(d, a)]).abs() / 6.0
                })
                .sum::<f64>()
        };
        let exact = with_cx(|cx| {
            let chart = fs_rep_frep::AxisymmetricChart::squat_disc(
                0.038,
                0.006,
                fs_rep_frep::SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
            )
            .unwrap();
            chart.mass_properties(1.0, cx).unwrap().mass
        });
        let coarse = with_cx(|cx| rounded_cylinder_tet_mesh(spec(3, 2, 12), cx)).unwrap();
        let fine = with_cx(|cx| rounded_cylinder_tet_mesh(spec(8, 5, 36), cx)).unwrap();
        let coarse_error = (volume(&coarse) - exact).abs();
        let fine_error = (volume(&fine) - exact).abs();
        assert!(fine_error < coarse_error, "{fine_error} !< {coarse_error}");
        assert!(
            fine_error / exact < 0.01,
            "relative error {}",
            fine_error / exact
        );
    }

    fn determinant(columns: [[f64; 3]; 3]) -> f64 {
        dot(columns[0], cross(columns[1], columns[2]))
    }
}
