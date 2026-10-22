//! IRP M1.3-W4: solid parts from descriptions reach audited tet meshes.
//!
//! Fixtures are declared synthetic designs (named below), not specimens: a
//! violin-bridge-like outline with three holes, and a xylophone-like bar with
//! a smooth undercut and two cord holes across the width. Real profiles come
//! later from measured parts.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_feec::highorder::simplex::duffy_quadrature;
use fs_mesh::planar::{PlanarBudget, PlanarError};
use fs_mesh::solid_parts::{
    DepthProfile, HermiteKnot, OutlineSolidSpec, ProfileBarSpec, ReferenceState, SnappedTetMesh,
    SolidPart, SolidPartError, SolidPartOptions, ThicknessField, ThroughHole, outline_solid,
    profile_bar,
};
use fs_mesh::tet_improve::ImproveOptions;
use fs_mesh::{
    AuditedLabeledTetComplex, RecoveryOptions, RefinementOptions, RegionId, VolumetricError,
    VolumetricPolicy, volumetricize,
};
use fs_solid::linear3::tet10::{promote_p1_to_p2, tet10_point};
use std::collections::BTreeMap;
use std::time::Instant;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x5011D,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

// ---- Fixture: violin-bridge-like outline, declared synthetic design, in
// units of 2^-12 m (0.244 mm) so every coordinate is exactly representable
// and the dyadic section grid leaves it unchanged.

const BRIDGE_UNIT_M: f64 = 1.0 / 4096.0;
const BRIDGE_OUTLINE: [[f64; 2]; 19] = [
    [0.0, 0.0],
    [49.0, 0.0],
    [49.0, 10.0],
    [119.0, 10.0],
    [119.0, 0.0],
    [168.0, 0.0],
    [168.0, 16.0],
    [152.0, 33.0],
    [147.0, 57.0],
    [164.0, 82.0],
    [156.0, 111.0],
    [127.0, 127.0],
    [84.0, 133.0],
    [41.0, 127.0],
    [12.0, 111.0],
    [4.0, 82.0],
    [20.0, 57.0],
    [16.0, 33.0],
    [0.0, 16.0],
];
const BRIDGE_HEART: [[f64; 2]; 8] = [
    [74.0, 61.0],
    [84.0, 53.0],
    [94.0, 61.0],
    [96.0, 76.0],
    [90.0, 84.0],
    [84.0, 78.0],
    [78.0, 84.0],
    [72.0, 76.0],
];
const BRIDGE_KIDNEY: [[f64; 2]; 5] = [
    [37.0, 39.0],
    [53.0, 37.0],
    [59.0, 47.0],
    [53.0, 55.0],
    [41.0, 53.0],
];
const BRIDGE_THICKNESS_M: f64 = 4.2e-3;

fn scaled(ring: &[[f64; 2]], unit_m: f64) -> Vec<[f64; 2]> {
    ring.iter()
        .map(|p| [p[0] * unit_m, p[1] * unit_m])
        .collect()
}

fn bridge_spec_in(unit_m: f64, thickness: ThicknessField) -> OutlineSolidSpec {
    let right_kidney: Vec<[f64; 2]> = BRIDGE_KIDNEY.iter().map(|p| [168.0 - p[0], p[1]]).collect();
    OutlineSolidSpec {
        outline: scaled(&BRIDGE_OUTLINE, unit_m),
        holes: vec![
            scaled(&BRIDGE_HEART, unit_m),
            scaled(&BRIDGE_KIDNEY, unit_m),
            scaled(&right_kidney, unit_m),
        ],
        thickness,
        reference_state: Some(ReferenceState::Unloaded),
    }
}

fn bridge_spec(thickness: ThicknessField) -> OutlineSolidSpec {
    bridge_spec_in(BRIDGE_UNIT_M, thickness)
}

fn perimeter(ring: &[[f64; 2]]) -> f64 {
    let n = ring.len();
    (0..n)
        .map(|i| {
            let (a, b) = (ring[i], ring[(i + 1) % n]);
            (b[0] - a[0]).hypot(b[1] - a[1])
        })
        .sum()
}

fn bridge_options() -> SolidPartOptions {
    SolidPartOptions {
        chord_tolerance_m: 1e-4,
        maximum_edge_m: 2e-3,
        minimum_angle_deg: 20.0,
        // 2^-40 m: about 12 bits of midpoint headroom on a 0.3 m section.
        coordinate_grid_m: 1.0 / 1_099_511_627_776.0,
        planar_budget: PlanarBudget {
            max_vertices: 20_000,
            max_triangles: 40_000,
            max_refinement_steps: 20_000,
        },
        max_surface_triangles: 200_000,
    }
}

/// Shoelace area, computed here independently of the generator.
fn shoelace(ring: &[[f64; 2]]) -> f64 {
    let n = ring.len();
    (0..n)
        .map(|i| ring[i][0] * ring[(i + 1) % n][1] - ring[(i + 1) % n][0] * ring[i][1])
        .sum::<f64>()
        .abs()
        / 2.0
}

// ---- Fixture: xylophone-like undercut bar, declared synthetic design.

const BAR_LENGTH_M: f64 = 0.30;
const BAR_WIDTH_M: f64 = 0.038;
const BAR_DEPTH_M: f64 = 0.020;
const BAR_UNDERCUT_DEPTH_M: f64 = 0.010;
/// Cord holes at 0.224 L (the free-free fundamental's nodes of a uniform bar).
const BAR_HOLES: [ThroughHole; 2] = [
    ThroughHole {
        x_m: 0.0672,
        z_m: 0.010,
        radius_m: 0.0025,
    },
    ThroughHole {
        x_m: 0.2328,
        z_m: 0.010,
        radius_m: 0.0025,
    },
];
/// Declared quality floor of every rung: the census's own sliver line.
const MIN_DIHEDRAL_FLOOR_DEG: f64 = 5.0;

fn bar_spec() -> ProfileBarSpec {
    let knot = |x_m, depth_m| HermiteKnot {
        x_m,
        depth_m,
        slope: 0.0,
    };
    ProfileBarSpec {
        length_m: BAR_LENGTH_M,
        width_m: BAR_WIDTH_M,
        // Flat ends, then a C1 smoothstep undercut to the centre and back.
        depth_profile: DepthProfile {
            knots: vec![
                knot(0.0, BAR_DEPTH_M),
                knot(0.09, BAR_DEPTH_M),
                knot(0.15, BAR_UNDERCUT_DEPTH_M),
                knot(0.21, BAR_DEPTH_M),
                knot(BAR_LENGTH_M, BAR_DEPTH_M),
            ],
        },
        through_holes: BAR_HOLES.to_vec(),
        reference_state: Some(ReferenceState::Unloaded),
    }
}

fn bar_options() -> SolidPartOptions {
    SolidPartOptions {
        chord_tolerance_m: 1e-4,
        maximum_edge_m: 6e-3,
        ..bridge_options()
    }
}

// ---- Independent mesh measurements.

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

fn signed_volume(p: [[f64; 3]; 4]) -> f64 {
    dot(sub(p[1], p[0]), cross(sub(p[2], p[0]), sub(p[3], p[0]))) / 6.0
}

fn min_dihedral_deg(p: [[f64; 3]; 4]) -> f64 {
    let faces = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
    let normals: Vec<[f64; 3]> = faces
        .iter()
        .enumerate()
        .map(|(k, f)| {
            let n = cross(sub(p[f[1]], p[f[0]]), sub(p[f[2]], p[f[0]]));
            if dot(n, sub(p[k], p[f[0]])) > 0.0 {
                n.map(|v| -v)
            } else {
                n
            }
        })
        .collect();
    let mut best = 180.0f64;
    for k in 0..4 {
        for l in k + 1..4 {
            let c = dot(normals[k], normals[l])
                / (dot(normals[k], normals[k]) * dot(normals[l], normals[l])).sqrt();
            best = best.min(180.0 - c.clamp(-1.0, 1.0).acos().to_degrees());
        }
    }
    best
}

/// Neumaier-compensated running sum.
#[derive(Default)]
struct Sum {
    s: f64,
    c: f64,
}

impl Sum {
    fn add(&mut self, v: f64) {
        let t = self.s + v;
        self.c += if self.s.abs() >= v.abs() {
            (self.s - t) + v
        } else {
            (v - t) + self.s
        };
        self.s = t;
    }
    fn total(&self) -> f64 {
        self.s + self.c
    }
}

fn sum(values: impl Iterator<Item = f64>) -> f64 {
    let mut s = Sum::default();
    values.for_each(|v| s.add(v));
    s.total()
}

struct Topology {
    boundary_genus: i64,
    solid_euler: i64,
}

fn topology(tets: &[[u32; 4]]) -> Topology {
    let mut faces: BTreeMap<[u32; 3], u32> = BTreeMap::new();
    let mut edges = std::collections::BTreeSet::new();
    let mut vertices = std::collections::BTreeSet::new();
    for t in tets {
        for f in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
            let mut key = f.map(|i| t[i]);
            key.sort_unstable();
            *faces.entry(key).or_insert(0) += 1;
        }
        for (a, b) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
            edges.insert([t[a].min(t[b]), t[a].max(t[b])]);
        }
        vertices.extend(t.iter().copied());
    }
    let boundary: Vec<[u32; 3]> = faces
        .iter()
        .filter(|(_, c)| **c == 1)
        .map(|(f, _)| *f)
        .collect();
    let bv: std::collections::BTreeSet<u32> = boundary.iter().flatten().copied().collect();
    let be: std::collections::BTreeSet<[u32; 2]> = boundary
        .iter()
        .flat_map(|f| [[f[0], f[1]], [f[0], f[2]], [f[1], f[2]]])
        .collect();
    let chi = bv.len() as i64 - be.len() as i64 + boundary.len() as i64;
    Topology {
        boundary_genus: (2 - chi) / 2,
        solid_euler: vertices.len() as i64 - edges.len() as i64 + faces.len() as i64
            - tets.len() as i64,
    }
}

fn p1_volume(positions: &[[f64; 3]], tets: &[[u32; 4]]) -> f64 {
    let volumes: Vec<f64> = tets
        .iter()
        .map(|t| signed_volume(t.map(|v| positions[v as usize])))
        .collect();
    let positive = volumes.iter().filter(|v| **v > 0.0).count();
    assert!(
        positive == 0 || positive == volumes.len(),
        "tets are not consistently oriented"
    );
    sum(volumes.iter().map(|v| v.abs()))
}

/// Divergence-theorem volume of the tets' boundary faces (outward by the
/// opposite vertex), independent of the tet sum.
fn boundary_volume(positions: &[[f64; 3]], tets: &[[u32; 4]]) -> f64 {
    let local = |t: &[u32; 4]| {
        [
            ([t[1], t[2], t[3]], t[0]),
            ([t[0], t[2], t[3]], t[1]),
            ([t[0], t[1], t[3]], t[2]),
            ([t[0], t[1], t[2]], t[3]),
        ]
    };
    let sorted = |mut f: [u32; 3]| {
        f.sort_unstable();
        f
    };
    let mut count: BTreeMap<[u32; 3], u32> = BTreeMap::new();
    for t in tets {
        for (f, _) in local(t) {
            *count.entry(sorted(f)).or_insert(0) += 1;
        }
    }
    let mut total = Sum::default();
    for t in tets {
        for (f, opposite) in local(t) {
            if count[&sorted(f)] != 1 {
                continue;
            }
            let [a, b, c] = f.map(|v| positions[v as usize]);
            let n = cross(sub(b, a), sub(c, a));
            let n = if dot(n, sub(positions[opposite as usize], a)) > 0.0 {
                n.map(|x| -x)
            } else {
                n
            };
            let centroid: [f64; 3] = core::array::from_fn(|k| (a[k] + b[k] + c[k]) / 3.0);
            total.add(dot(centroid, n) / 6.0);
        }
    }
    total.total()
}

fn min_dihedral(positions: &[[f64; 3]], tets: &[[u32; 4]]) -> f64 {
    tets.iter()
        .map(|t| min_dihedral_deg(t.map(|v| positions[v as usize])))
        .fold(180.0, f64::min)
}

/// Exact volume of a quadratic tet mesh: `det J` is cubic in the
/// barycentrics and Duffy n = 3 integrates it exactly.
struct QuadraticVolume {
    /// Sum over every element whose sampled Jacobian stayed positive.
    volume: f64,
    /// Smallest sampled `det J / h^3` over those elements.
    min_scaled_jacobian: f64,
    /// Elements `tet10_point` refused (folded or below 1e-9 scaled).
    refused: Vec<usize>,
}

fn p2_volume(nodes: &[[f64; 3]], elements: &[[usize; 10]]) -> QuadraticVolume {
    let rule = duffy_quadrature(3);
    let mut min_scaled = f64::INFINITY;
    let mut volume = Sum::default();
    let mut refused = Vec::new();
    for (e, el) in elements.iter().enumerate() {
        let x: [[f64; 3]; 10] = core::array::from_fn(|a| nodes[el[a]]);
        let h = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            .iter()
            .map(|&(a, b)| dot(sub(x[a], x[b]), sub(x[a], x[b])).sqrt())
            .fold(0.0, f64::max);
        let points: Result<Vec<_>, _> = rule
            .iter()
            .enumerate()
            .map(|(q, (lambda, w))| tet10_point(x, *lambda, 1e-9, e, q).map(|p| (p, *w)))
            .collect();
        match points {
            Ok(points) => {
                for (point, w) in points {
                    min_scaled = min_scaled.min(point.determinant_m3 / (h * h * h));
                    volume.add(w * point.determinant_m3);
                }
            }
            Err(_) => refused.push(e),
        }
    }
    QuadraticVolume {
        volume: volume.total(),
        min_scaled_jacobian: min_scaled,
        refused,
    }
}

fn audited_volume(audited: &AuditedLabeledTetComplex) -> f64 {
    audited.witness().per_region_auditor[0].1
}

// ---- G1: bridge outline with holes.

#[test]
fn g1_bridge_outline_with_holes_extruded() {
    let spec = bridge_spec(ThicknessField::Uniform {
        thickness_m: BRIDGE_THICKNESS_M,
    });
    let area = shoelace(&spec.outline) - spec.holes.iter().map(|h| shoelace(h)).sum::<f64>();
    let expected = area * BRIDGE_THICKNESS_M;
    with_cx(|cx| {
        let started = Instant::now();
        let part = outline_solid(&spec, &bridge_options(), cx).expect("bridge surface");
        let receipt = part.receipt();
        let rel = |v: f64| ((v - expected) / expected).abs();
        assert!(
            rel(receipt.exact_volume_m3) <= 1e-12,
            "{}",
            receipt.to_json()
        );
        assert!(
            rel(receipt.polyhedral_volume_m3) <= 1e-12,
            "{}",
            receipt.to_json()
        );
        assert_eq!(
            receipt.max_quantization_m, 0.0,
            "on-grid input is unchanged"
        );
        assert_eq!(receipt.genus, 3, "three holes, genus three");
        let policy = VolumetricPolicy {
            length_unit: "m".to_string(),
            recovery: RecoveryOptions::default(),
            max_vertices: 100_000,
            max_tets: 2_000_000,
            refinement: None,
        };
        let audited = volumetricize(part.plc(RegionId(1)), policy, cx).expect("volumetricize");
        let volume = audited_volume(&audited);
        assert!(
            rel(volume) <= 1e-12,
            "tet volume {volume:e} vs {expected:e}"
        );
        let topo = topology(audited.labeled().tets());
        assert_eq!(topo.boundary_genus, 3, "hole count preserved in the tets");
        assert_eq!(topo.solid_euler, 1 - 3);
        let census = audited.labeled().quality();
        eprintln!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.3-W4\",\"test\":\"g1_bridge_outline_with_holes_extruded\",\
             \"inputs\":{{\"source\":\"parametric\",\"spec\":\"bridge outline 19 vertices, 3 holes, t={BRIDGE_THICKNESS_M:e} m\"}},\
             \"part\":{},\"observed\":{{\"shoelace_area_m2\":{area:e},\"expected_volume_m3\":{expected:e},\
             \"tet_volume_m3\":{volume:e},\"volume_rel_err\":{:e},\"genus\":{},\"census\":{}}},\
             \"tolerance\":1e-12,\"verdict\":\"pass\",\"elapsed_s\":{:.3}}}",
            receipt.to_json(),
            rel(volume),
            topo.boundary_genus,
            census.to_json(),
            started.elapsed().as_secs_f64()
        );

        // An affine thickness field (feet 4.5 mm thinning to the top) is
        // still piecewise planar: volume = t0 A + gx Mx + gy My.
        let affine = bridge_spec(ThicknessField::Affine {
            at_origin_m: 4.5e-3,
            gradient: [0.0, -0.1],
        });
        let moment_y = |ring: &[[f64; 2]]| {
            let n = ring.len();
            let s: f64 = (0..n)
                .map(|i| {
                    let (a, b) = (ring[i], ring[(i + 1) % n]);
                    (a[1] + b[1]) * (a[0] * b[1] - b[0] * a[1]) / 6.0
                })
                .sum();
            let signed: f64 = (0..n)
                .map(|i| ring[i][0] * ring[(i + 1) % n][1] - ring[(i + 1) % n][0] * ring[i][1])
                .sum();
            s * signed.signum()
        };
        let my = moment_y(&affine.outline) - affine.holes.iter().map(|h| moment_y(h)).sum::<f64>();
        let expected_affine = 4.5e-3 * area - 0.1 * my;
        let part = outline_solid(&affine, &bridge_options(), cx).expect("affine bridge");
        let r = part.receipt();
        for v in [r.exact_volume_m3, r.polyhedral_volume_m3] {
            assert!(
                ((v - expected_affine) / expected_affine).abs() <= 1e-12,
                "affine volume {v:e} vs {expected_affine:e}"
            );
        }
        let audited = volumetricize(
            part.plc(RegionId(1)),
            VolumetricPolicy {
                length_unit: "m".to_string(),
                recovery: RecoveryOptions::default(),
                max_vertices: 100_000,
                max_tets: 2_000_000,
                refinement: None,
            },
            cx,
        )
        .expect("volumetricize affine bridge");
        let v = audited_volume(&audited);
        assert!(((v - expected_affine) / expected_affine).abs() <= 1e-12);
        assert_eq!(topology(audited.labeled().tets()).boundary_genus, 3);

        // Off-grid input (the same design in decimal 0.244 mm steps) is
        // rounded to the grid; the volume moves by at most thickness x
        // boundary length x twice the receipted rounding.
        let decimal = bridge_spec_in(
            2.44e-4,
            ThicknessField::Uniform {
                thickness_m: BRIDGE_THICKNESS_M,
            },
        );
        let part = outline_solid(&decimal, &bridge_options(), cx).expect("off-grid bridge");
        let r = part.receipt();
        assert!(r.max_quantization_m > 0.0);
        assert!(r.max_quantization_m <= 0.5 * r.coordinate_grid_m);
        let boundary: f64 =
            perimeter(&decimal.outline) + decimal.holes.iter().map(|h| perimeter(h)).sum::<f64>();
        let bound = BRIDGE_THICKNESS_M * boundary * 2.0 * r.max_quantization_m;
        assert!(
            (r.polyhedral_volume_m3 - r.exact_volume_m3).abs() <= bound,
            "{} vs bound {bound:e}",
            r.to_json()
        );
    });
}

// ---- G1: the undercut bar ladder.

struct Rung {
    tets: usize,
    p1_volume: f64,
    /// `None` when any curved element was refused.
    p2_volume: Option<f64>,
    min_dihedral_deg: f64,
    min_scaled_jacobian: f64,
    max_snap_m: f64,
    curved_midsides: usize,
    /// Tets below 5 degrees, by how many of their vertices are on the boundary.
    slivers_by_boundary_vertices: [usize; 5],
    /// Refused curved elements and the largest P1 min dihedral among them.
    refused_curved: usize,
    refused_max_dihedral_deg: f64,
}

fn measure_rung(part: &SolidPart, mesh: &SnappedTetMesh, check_rule: bool, cx: &Cx<'_>) -> Rung {
    let topo = topology(&mesh.tets);
    assert_eq!(topo.boundary_genus, 2, "two holes, genus two");
    assert_eq!(topo.solid_euler, 1 - 2);
    let nodes: Vec<[f64; 3]> = mesh.positions.clone();
    let tets: Vec<[usize; 4]> = mesh.tets.iter().map(|t| t.map(|v| v as usize)).collect();
    let projection = part
        .midside_projection(mesh, cx)
        .expect("midside projection");
    let (p2_nodes, p2_elements) =
        promote_p1_to_p2(&nodes, &tets, Some(&projection), cx).expect("TET10 promotion");
    let curved = p2_volume(&p2_nodes, &p2_elements);
    let p1_volume = p1_volume(&mesh.positions, &mesh.tets);
    // The tets tile the domain: their sum is the boundary's enclosed volume
    // (an overlap or a gap left by the quality pass would show here).
    let enclosed = boundary_volume(&mesh.positions, &mesh.tets);
    assert!(
        ((p1_volume - enclosed) / enclosed).abs() <= 1e-12,
        "tet sum {p1_volume:e} vs enclosed {enclosed:e}"
    );
    if check_rule {
        // Straight-sided TET10 reproduces the P1 volume: the rule is exact.
        let (flat_nodes, flat_elements) =
            promote_p1_to_p2(&nodes, &tets, None, cx).expect("straight TET10 promotion");
        let flat = p2_volume(&flat_nodes, &flat_elements);
        assert!(flat.refused.is_empty());
        assert!(((flat.volume - p1_volume) / p1_volume).abs() <= 1e-12);
    }
    let on_boundary: std::collections::BTreeSet<u32> = mesh
        .boundary_faces
        .iter()
        .flat_map(|(f, _)| f.iter().copied())
        .collect();
    let mut slivers_by_boundary_vertices = [0usize; 5];
    for t in &mesh.tets {
        if min_dihedral_deg(t.map(|v| mesh.positions[v as usize])) < 5.0 {
            slivers_by_boundary_vertices[t.iter().filter(|v| on_boundary.contains(v)).count()] += 1;
        }
    }
    let refused_max_dihedral_deg = curved
        .refused
        .iter()
        .map(|&e| min_dihedral_deg(mesh.tets[e].map(|v| mesh.positions[v as usize])))
        .fold(0.0, f64::max);
    Rung {
        tets: mesh.tets.len(),
        p1_volume,
        p2_volume: curved.refused.is_empty().then_some(curved.volume),
        min_dihedral_deg: min_dihedral(&mesh.positions, &mesh.tets),
        min_scaled_jacobian: curved.min_scaled_jacobian,
        max_snap_m: mesh.max_snap_m,
        curved_midsides: projection.curved_edges(),
        slivers_by_boundary_vertices,
        refused_curved: curved.refused.len(),
        refused_max_dihedral_deg,
    }
}

#[test]
fn g1_undercut_bar_with_nodal_holes_from_parameters() {
    with_cx(|cx| {
        let started = Instant::now();
        let spec = bar_spec();
        let part = profile_bar(&spec, &bar_options(), cx).expect("bar surface");
        let receipt = part.receipt().clone();
        assert_eq!(receipt.genus, 2);
        let exact = receipt.exact_volume_m3;
        // Independent exact volume: Simpson is exact on each cubic piece.
        let h = |x: f64| spec.depth_profile.depth(x);
        let knots: Vec<f64> = spec.depth_profile.knots.iter().map(|k| k.x_m).collect();
        let integral: f64 = knots
            .windows(2)
            .map(|w| (w[1] - w[0]) / 6.0 * (h(w[0]) + 4.0 * h((w[0] + w[1]) / 2.0) + h(w[1])))
            .sum();
        let holes = 2.0 * std::f64::consts::PI * 0.0025 * 0.0025;
        assert!(((exact - BAR_WIDTH_M * (integral - holes)) / exact).abs() <= 1e-12);
        let policy = VolumetricPolicy {
            length_unit: "m".to_string(),
            recovery: RecoveryOptions::default(),
            max_vertices: 100_000,
            max_tets: 2_000_000,
            refinement: Some(RefinementOptions::default()),
        };
        let audited = volumetricize(part.plc(RegionId(7)), policy, cx).expect("volumetricize bar");
        let witnessed = audited_volume(&audited);
        assert!(
            ((witnessed - receipt.polyhedral_volume_m3) / witnessed).abs() <= 1e-9,
            "audited tets fill the PLC"
        );
        eprintln!(
            "bar part {} volumetricized in {:.2} s: census {}",
            receipt.to_json(),
            started.elapsed().as_secs_f64(),
            audited.labeled().quality().to_json()
        );
        // Quality pass after every snap, aimed at the declared floor: peel
        // boundary slivers (every boundary vertex stays on the exact surface),
        // then smooth interior vertices and flip interior faces. Volumes depend
        // on the boundary only; the floor is still checked below on what the
        // pass reached.
        let improve = Some(ImproveOptions {
            target_min_dihedral_deg: MIN_DIHEDRAL_FLOOR_DEG,
            max_passes: 32,
        });
        let mut mesh = part
            .snap(audited.labeled(), improve, cx)
            .expect("snap rung 0");
        eprintln!(
            "volumes: plc {:.15e} audited {:.15e} rung 0 {:.15e}",
            receipt.polyhedral_volume_m3,
            witnessed,
            p1_volume(&mesh.positions, &mesh.tets),
        );
        let mut rungs = Vec::new();
        for k in 0..3 {
            if k > 0 {
                mesh = part
                    .refine_snapped(&mesh, 4_000_000, improve, cx)
                    .expect("refine rung");
            }
            let t = Instant::now();
            let rung = measure_rung(&part, &mesh, k == 0, cx);
            eprintln!(
                "rung {k} quality pass: peeled {} then {}",
                mesh.peeled,
                mesh.improvement.map_or("none".to_string(), |s| s.to_json())
            );
            eprintln!(
                "rung {k}: tets {} p1 {:.15e} p2 {:?} min_dihedral {:.3} min_scaled_jac {:.3e} max_snap {:.3e} curved_midsides {} slivers<5deg by boundary vertices {:?} refused curved {} (largest P1 dihedral among them {:.3}) ({:.2} s)",
                rung.tets,
                rung.p1_volume,
                rung.p2_volume,
                rung.min_dihedral_deg,
                rung.min_scaled_jacobian,
                rung.max_snap_m,
                rung.curved_midsides,
                rung.slivers_by_boundary_vertices,
                rung.refused_curved,
                rung.refused_max_dihedral_deg,
                t.elapsed().as_secs_f64()
            );
            rungs.push(rung);
        }
        let err = |v: f64| ((exact - v) / exact).abs();
        let e1: Vec<f64> = rungs.iter().map(|r| err(r.p1_volume)).collect();
        // A rung with a refused curved element has no curved volume (NaN).
        let e2: Vec<f64> = rungs
            .iter()
            .map(|r| r.p2_volume.map_or(f64::NAN, err))
            .collect();
        let order = |e: &[f64]| [(e[0] / e[1]).log2(), (e[1] / e[2]).log2()];
        let (o1, o2) = (order(&e1), order(&e2));
        eprintln!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.3-W4\",\"test\":\"g1_undercut_bar_with_nodal_holes_from_parameters\",\
             \"inputs\":{{\"source\":\"parametric\",\"spec\":\"L={BAR_LENGTH_M} W={BAR_WIDTH_M} h=Hermite[{BAR_DEPTH_M},{BAR_UNDERCUT_DEPTH_M}] holes=2 r=0.0025\"}},\
             \"part\":{},\"rungs\":{{\"axis\":\"h\",\"refinement\":\"uniform 1->8\",\"tets\":{:?}}},\
             \"observed\":{{\"exact_volume_m3\":{exact:e},\"volume_rel_err_p1\":{e1:?},\"volume_rel_err_p2\":{e2:?},\
             \"min_dihedral_deg\":{:?},\"min_scaled_jacobian\":{:?},\"max_snap_m\":{:?},\"refused_curved\":{:?},\"slivers_by_boundary_vertices\":{:?},\"genus\":2}},\
             \"order\":{{\"p1\":{o1:?},\"p2\":{o2:?}}},\"tolerance\":{{\"p1\":\"2+-0.2\",\"p2\":\">=2.8\",\
             \"min_dihedral_floor_deg\":{MIN_DIHEDRAL_FLOOR_DEG}}},\"elapsed_s\":{:.3}}}",
            receipt.to_json(),
            rungs.iter().map(|r| r.tets).collect::<Vec<_>>(),
            rungs.iter().map(|r| r.min_dihedral_deg).collect::<Vec<_>>(),
            rungs
                .iter()
                .map(|r| r.min_scaled_jacobian)
                .collect::<Vec<_>>(),
            rungs.iter().map(|r| r.max_snap_m).collect::<Vec<_>>(),
            rungs.iter().map(|r| r.refused_curved).collect::<Vec<_>>(),
            rungs
                .iter()
                .map(|r| r.slivers_by_boundary_vertices)
                .collect::<Vec<_>>(),
            started.elapsed().as_secs_f64()
        );
        for (k, rung) in rungs.iter().enumerate() {
            assert_eq!(
                rung.refused_curved, 0,
                "rung {k} folded curved TET10 elements"
            );
        }
        assert!(
            (o1[1] - 2.0).abs() <= 0.2,
            "straight-sided volume order {o1:?} (errors {e1:?})"
        );
        assert!(
            o2[1] >= 2.8,
            "projected volume order {o2:?} (errors {e2:?})"
        );
        for (k, rung) in rungs.iter().enumerate() {
            assert!(
                rung.min_dihedral_deg >= MIN_DIHEDRAL_FLOOR_DEG,
                "rung {k} min dihedral {} below the declared floor",
                rung.min_dihedral_deg
            );
            assert!(rung.curved_midsides > 0);
        }
    });
}

// ---- G0: generator properties and refusals.

fn closed_and_oriented(part: &SolidPart) {
    let mut directed: BTreeMap<[u32; 2], i32> = BTreeMap::new();
    for t in part.triangles() {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *directed.entry([a, b]).or_insert(0) += 1;
        }
    }
    for ([a, b], count) in &directed {
        assert_eq!(*count, 1, "edge ({a},{b}) used twice in one direction");
        assert_eq!(directed.get(&[*b, *a]), Some(&1), "edge ({a},{b}) is open");
    }
    assert!(
        part.receipt().polyhedral_volume_m3 > 0.0,
        "outward orientation"
    );
}

fn bits(part: &SolidPart) -> Vec<[u64; 3]> {
    part.vertices()
        .iter()
        .map(|p| p.map(f64::to_bits))
        .collect()
}

#[test]
fn g0_profile_generator_is_deterministic_and_closed() {
    with_cx(|cx| {
        let a = profile_bar(&bar_spec(), &bar_options(), cx).expect("bar");
        let b = profile_bar(&bar_spec(), &bar_options(), cx).expect("bar again");
        assert_eq!(bits(&a), bits(&b));
        assert_eq!(a.triangles(), b.triangles());
        assert_eq!(a.triangle_patches(), b.triangle_patches());
        assert_eq!(a, b);
        closed_and_oriented(&a);
        let spec = bridge_spec(ThicknessField::Uniform {
            thickness_m: BRIDGE_THICKNESS_M,
        });
        let c = outline_solid(&spec, &bridge_options(), cx).expect("bridge");
        let d = outline_solid(&spec, &bridge_options(), cx).expect("bridge again");
        assert_eq!(bits(&c), bits(&d));
        assert_eq!(c, d);
        closed_and_oriented(&c);
        // Profiled-face vertices are grid-rounded samples of h(x) or planar
        // mesher splits of a chord: within the chord tolerance plus the grid.
        let profile = bar_spec().depth_profile;
        let slack = bar_options().chord_tolerance_m + bar_options().coordinate_grid_m;
        for (t, patch) in a.triangles().iter().zip(a.triangle_patches()) {
            if *patch == fs_mesh::solid_parts::SurfacePatch::DepthProfile {
                for &v in t {
                    let p = a.vertices()[v as usize];
                    let gap = (p[2] - profile.depth(p[0])).abs();
                    assert!(gap <= slack, "sagitta {gap:e}");
                }
            }
        }
    });
}

#[test]
fn g0_refuses_self_intersecting_outline_and_missing_reference_state() {
    with_cx(|cx| {
        let mut bow_tie = bridge_spec(ThicknessField::Uniform {
            thickness_m: BRIDGE_THICKNESS_M,
        });
        bow_tie.outline = vec![[0.0, 0.0], [0.04, 0.03], [0.04, 0.0], [0.0, 0.03]];
        bow_tie.holes.clear();
        assert!(matches!(
            outline_solid(&bow_tie, &bridge_options(), cx),
            Err(SolidPartError::Planar(PlanarError::SelfIntersection { .. }))
        ));
        let mut undeclared = bridge_spec(ThicknessField::Uniform {
            thickness_m: BRIDGE_THICKNESS_M,
        });
        undeclared.reference_state = None;
        assert_eq!(
            outline_solid(&undeclared, &bridge_options(), cx),
            Err(SolidPartError::ReferenceStateMissing)
        );
        let mut bar = bar_spec();
        bar.reference_state = None;
        assert_eq!(
            profile_bar(&bar, &bar_options(), cx),
            Err(SolidPartError::ReferenceStateMissing)
        );
        // A hole breaking through the undercut face, and one crossing the end.
        for hole in [
            ThroughHole {
                x_m: 0.15,
                z_m: 0.008,
                radius_m: 0.0025,
            },
            ThroughHole {
                x_m: 0.001,
                z_m: 0.010,
                radius_m: 0.0025,
            },
        ] {
            let mut bar = bar_spec();
            bar.through_holes.push(hole);
            assert!(matches!(
                profile_bar(&bar, &bar_options(), cx),
                Err(SolidPartError::InvalidSpec(_))
            ));
        }
        let mut thin = bridge_spec(ThicknessField::Affine {
            at_origin_m: 1e-3,
            gradient: [0.0, -0.1],
        });
        thin.holes.clear();
        assert!(matches!(
            outline_solid(&thin, &bridge_options(), cx),
            Err(SolidPartError::InvalidSpec(_))
        ));
    });
}

#[test]
fn g0_policy_cap_and_cap_plus_one() {
    with_cx(|cx| {
        let spec = bridge_spec(ThicknessField::Uniform {
            thickness_m: BRIDGE_THICKNESS_M,
        });
        let part = outline_solid(&spec, &bridge_options(), cx).expect("bridge");
        let triangles = part.triangles().len();
        let at_cap = SolidPartOptions {
            max_surface_triangles: triangles,
            ..bridge_options()
        };
        assert_eq!(outline_solid(&spec, &at_cap, cx).expect("at cap"), part);
        let over = SolidPartOptions {
            max_surface_triangles: triangles - 1,
            ..bridge_options()
        };
        assert!(matches!(
            outline_solid(&spec, &over, cx),
            Err(SolidPartError::Budget {
                what: "surface triangles",
                ..
            })
        ));
        let vertices = part.vertices().len();
        let policy = |max_vertices, max_tets| VolumetricPolicy {
            length_unit: "m".to_string(),
            recovery: RecoveryOptions::default(),
            max_vertices,
            max_tets,
            refinement: None,
        };
        let full = volumetricize(part.plc(RegionId(1)), policy(vertices, 2_000_000), cx)
            .expect("vertices exactly at cap");
        // The tet cap binds where carving first counts the retained tets
        // (flip repair can remove some afterwards).
        let carved = part
            .plc(RegionId(1))
            .admit(policy(vertices, 2_000_000), cx)
            .and_then(|admitted| admitted.recover(cx))
            .and_then(|recovered| recovered.carve_and_label(cx))
            .expect("carve")
            .tets()
            .len();
        assert!(matches!(
            volumetricize(part.plc(RegionId(1)), policy(vertices - 1, 2_000_000), cx),
            Err(VolumetricError::Budget {
                what: "vertices",
                ..
            })
        ));
        let at_tets = volumetricize(part.plc(RegionId(1)), policy(vertices, carved), cx)
            .expect("tets exactly at cap");
        assert!(at_tets.labeled().tets().len() <= carved);
        assert!(matches!(
            volumetricize(part.plc(RegionId(1)), policy(vertices, carved - 1), cx),
            Err(VolumetricError::Budget { what: "tets", .. })
        ));
        let tets = full.labeled().tets().len();
        // The snapped ladder's own cap.
        let mesh = part.snap(full.labeled(), None, cx).expect("snap");
        assert_eq!(mesh.max_snap_m, 0.0, "a polyhedral part needs no snap");
        let refined = part
            .refine_snapped(&mesh, 8 * tets, None, cx)
            .expect("refined exactly at cap");
        assert_eq!(refined.tets.len(), 8 * tets);
        assert!(matches!(
            part.refine_snapped(&mesh, 8 * tets - 1, None, cx),
            Err(SolidPartError::Budget {
                what: "refined tets",
                ..
            })
        ));
    });
}
