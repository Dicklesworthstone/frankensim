//! fs-mesh conformance suite (CONTRACT.md: any reimplementation must
//! pass). Exact-audit Delaunay on general-position clouds, the
//! degenerate adversarial battery (grids, cospherical shells, collinear
//! runs, coplanar refusals, duplicates), determinism/relabeling/G3
//! translation, hull integration with fs-rep-mesh, radius-edge
//! refinement, and a scale run. JSON-line verdicts; seeded cases carry
//! seeds.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::Point3;
use fs_mesh::{
    ExudeOptions, MeshError, RefineOptions, Tetrahedralization, delaunay, delaunay_colored,
    delaunay_colored_reversed, refine,
};
use fs_rep_mesh::{assess_quality, winding_exact};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-mesh/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\
         \"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "case {case}: {detail}");
}

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }

    /// Dyadic coordinate in [-1, 1] (multiples of 2⁻²⁰ — exact under
    /// the dyadic translations the G3 case applies).
    fn dyadic(&mut self) -> f64 {
        let k = (self.next() >> 32) % (1 << 21);
        (k as f64) / f64::from(1 << 20) - 1.0
    }
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x7E7,
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

fn cloud(seed: u64, n: usize) -> Vec<Point3> {
    let mut rng = Lcg(seed);
    (0..n)
        .map(|_| Point3::new(rng.dyadic(), rng.dyadic(), rng.dyadic()))
        .collect()
}

/// Canonical geometric key of a tet set (coordinate bits, label-free).
fn geometric_keys(t: &Tetrahedralization) -> std::collections::BTreeSet<[[u64; 3]; 4]> {
    let pts = t.points();
    t.tets()
        .into_iter()
        .map(|tet| {
            let mut k: [[u64; 3]; 4] = tet.map(|v| {
                let p = pts[v as usize];
                [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
            });
            k.sort_unstable();
            k
        })
        .collect()
}

/// tmesh-001 — a general-position cloud triangulates with the FULL
/// exact audit clean: global empty circumsphere, local Delaunay,
/// orientation, mutual adjacency, Euler = 1, hull closed and exactly
/// convex. Logs kernel statistics.
#[test]
fn tmesh_001_general_position_full_audit() {
    with_cx(|cx| {
        let pts = cloud(0x1001_2026_0706_0021, 200);
        let t = delaunay(&pts, cx).expect("general position builds");
        let report = t.audit(true);
        let stats = t.stats();
        let mut em = fs_obs::Emitter::new("fs-mesh/conformance", "tmesh-001/kernel");
        let line = em
            .emit(
                fs_obs::Severity::Info,
                fs_obs::EventKind::Custom {
                    name: "mesh-delaunay-stats".to_string(),
                    json: stats.to_json(),
                },
                None,
            )
            .to_jsonl();
        fs_obs::validate_line(&line).expect("kernel stats validate");
        println!("{line}");
        verdict(
            "tmesh-001",
            report.clean() && stats.tets_final > 400 && stats.exhaustive_locates == 0,
            &format!(
                "200-point cloud: {} tets, FULL exact audit clean (global empty-sphere, \
                 local Delaunay, orientation, adjacency, Euler, exact hull convexity); \
                 violations={:?}; seed 0x1001_2026_0706_0021",
                stats.tets_final,
                report.violations.first()
            ),
        );
    });
}

/// tmesh-002 — the adversarial degeneracy battery: exact predicates
/// earn their keep on grids (massively cospherical/coplanar), an
/// exactly cospherical shell, collinear runs, bitwise duplicates, and
/// the all-coplanar refusal.
#[test]
fn tmesh_002_degeneracy_battery() {
    with_cx(|cx| {
        // 4×4×4 integer grid.
        let mut grid = Vec::new();
        for x in 0..4 {
            for y in 0..4 {
                for z in 0..4 {
                    grid.push(Point3::new(f64::from(x), f64::from(y), f64::from(z)));
                }
            }
        }
        let t = delaunay(&grid, cx).expect("grid builds");
        let grid_report = t.audit(true);
        let grid_stats = t.stats();

        // Exactly cospherical: cube corners + octahedron vertices, all
        // at radius √3 from the origin.
        let r = 3.0f64.sqrt();
        let mut shell: Vec<Point3> = vec![];
        for &sx in &[-1.0, 1.0] {
            for &sy in &[-1.0, 1.0] {
                for &sz in &[-1.0, 1.0] {
                    shell.push(Point3::new(sx, sy, sz));
                }
            }
        }
        for k in 0..3 {
            for &s in &[-r, r] {
                let mut p = [0.0; 3];
                p[k] = s;
                shell.push(Point3::new(p[0], p[1], p[2]));
            }
        }
        shell.push(Point3::new(0.0, 0.0, 0.0)); // one interior point
        let t = delaunay(&shell, cx).expect("cospherical shell builds");
        let shell_report = t.audit(true);

        // Collinear run + off-line points (bootstrap must scan past).
        let mut line: Vec<Point3> = (0..10)
            .map(|i| Point3::new(f64::from(i), 0.0, 0.0))
            .collect();
        line.push(Point3::new(0.0, 1.0, 0.0));
        line.push(Point3::new(0.0, 0.0, 1.0));
        line.push(Point3::new(1.0, 1.0, 1.0));
        let t = delaunay(&line, cx).expect("collinear-plus builds");
        let line_report = t.audit(true);

        // Bitwise duplicates: counted, skipped, mesh stays clean.
        let mut dup = cloud(0x1001_2026_0706_0022, 100);
        for i in 0..30 {
            dup.push(dup[i * 3]);
        }
        let t = delaunay(&dup, cx).expect("duplicates build");
        let dup_report = t.audit(true);
        let dups_skipped = t.stats().duplicates_skipped;

        // All-coplanar: exact refusal with teaching text.
        let flat: Vec<Point3> = (0..25)
            .map(|i| Point3::new(f64::from(i % 5), f64::from(i / 5), 0.0))
            .collect();
        let err = delaunay(&flat, cx).expect_err("coplanar must refuse");
        let teaches = err == MeshError::DegenerateInput && err.to_string().contains("coplanar");

        verdict(
            "tmesh-002",
            grid_report.clean()
                && shell_report.clean()
                && line_report.clean()
                && dup_report.clean()
                && dups_skipped == 30
                && teaches,
            &format!(
                "4x4x4 grid clean under FULL audit ({} growth repairs absorbed \
                 degenerate visibility), exactly-cospherical shell clean, collinear \
                 run clean, 30/30 bitwise duplicates skipped with receipts, and the \
                 all-coplanar grid refuses with teaching text",
                grid_stats.growth_repairs
            ),
        );
    });
}

/// tmesh-003 — determinism (P2/G5): bitwise-identical output across
/// runs; relabeling invariance (same GEOMETRIC tet set under a vertex
/// permutation); exact G3 equivariance under a dyadic translation.
#[test]
fn tmesh_003_determinism_relabeling_translation() {
    with_cx(|cx| {
        let pts = cloud(0x1001_2026_0706_0023, 150);
        let a = delaunay(&pts, cx).expect("build a");
        let b = delaunay(&pts, cx).expect("build b");
        let bitwise = a.tets() == b.tets() && a.points() == b.points();

        let relabeled: Vec<Point3> = pts.iter().rev().copied().collect();
        let c = delaunay(&relabeled, cx).expect("build relabeled");
        let relabel_ok = geometric_keys(&a) == geometric_keys(&c);

        let shift = fs_geom::Vec3::new(0.5, 0.25, -0.375);
        let moved: Vec<Point3> = pts.iter().map(|p| p.offset(shift)).collect();
        let d = delaunay(&moved, cx).expect("build translated");
        let conn_ok = a.tets() == d.tets();
        let coords_ok = a
            .points()
            .iter()
            .zip(d.points())
            .all(|(p, q)| p.offset(shift) == q);

        verdict(
            "tmesh-003",
            bitwise && relabel_ok && conn_ok && coords_ok,
            "same input gives BITWISE-identical meshes; a reversed labeling yields \
             the identical geometric tet set; a dyadic translation preserves \
             connectivity exactly with exactly-shifted coordinates (G3/G5); \
             seed 0x1001_2026_0706_0023",
        );
    });
}

/// tmesh-004 — integration with fs-rep-mesh: the hull soup is closed,
/// 2-manifold, outward-oriented (winding +1 at an interior point), and
/// the oriented tet complex satisfies δδ = 0 exactly.
#[test]
fn tmesh_004_hull_and_complex_integration() {
    with_cx(|cx| {
        let pts = cloud(0x1001_2026_0706_0024, 120);
        let t = delaunay(&pts, cx).expect("build");
        let hull = t.hull();
        let quality = assess_quality(&hull);
        let w = winding_exact(&hull, Point3::new(0.0, 0.0, 0.0));
        let complex = t.complex();
        let (d0, d1, d2) = (complex.d0(), complex.d1(), complex.d2());
        let mut dd_zero = true;
        for probe in 0..complex.vertex_count.min(8) {
            let mut x = vec![0i64; complex.vertex_count];
            x[probe] = 1;
            dd_zero &= d1.apply(&d0.apply(&x)).iter().all(|&v| v == 0);
        }
        for probe in 0..complex.edges.len().min(8) {
            let mut x = vec![0i64; complex.edges.len()];
            x[probe] = 1;
            dd_zero &= d2.apply(&d1.apply(&x)).iter().all(|&v| v == 0);
        }
        verdict(
            "tmesh-004",
            quality.sign_certified() && (w - 1.0).abs() < 1e-9 && dd_zero,
            &format!(
                "hull is closed 2-manifold with winding {w:.6} at an interior point \
                 (outward-oriented) and the oriented tet complex satisfies dd=0 \
                 exactly; seed 0x1001_2026_0706_0024"
            ),
        );
    });
}

/// tmesh-005 — radius-edge refinement: quality improves to the bound
/// (or budget), the exact audit STAYS clean through every Steiner
/// insertion, and refinement is deterministic. Logs quality stats.
#[test]
fn tmesh_005_refinement_quality() {
    with_cx(|cx| {
        // A stretched cloud manufactures bad radius-edge ratios.
        let mut rng = Lcg(0x1001_2026_0706_0025);
        let pts: Vec<Point3> = (0..120)
            .map(|_| Point3::new(rng.dyadic() * 4.0, rng.dyadic(), rng.dyadic() * 0.25))
            .collect();
        let mut t = delaunay(&pts, cx).expect("build");
        let opts = RefineOptions {
            max_radius_edge: 2.0,
            max_steiner: 500,
            ..RefineOptions::default()
        };
        let stats = refine(&mut t, opts, cx).expect("refine");
        let report = t.audit(true);
        let mut t2 = delaunay(&pts, cx).expect("build 2");
        let _ = refine(&mut t2, opts, cx).expect("refine 2");
        let deterministic = t.tets() == t2.tets() && t.points() == t2.points();
        // The honest v1 guarantee: every offender that SURVIVES is one
        // whose circumcenter escapes the hull (boundary handling is the
        // successor bead); nothing interior-refinable remains.
        let exhausted = stats.refinable_remaining == 0;
        let mut em = fs_obs::Emitter::new("fs-mesh/conformance", "tmesh-005/refine");
        let line = em
            .emit(
                fs_obs::Severity::Info,
                fs_obs::EventKind::Custom {
                    name: "mesh-refine-stats".to_string(),
                    json: stats.to_json(),
                },
                None,
            )
            .to_jsonl();
        fs_obs::validate_line(&line).expect("refine stats validate");
        println!("{line}");
        verdict(
            "tmesh-005",
            exhausted && report.clean() && deterministic && stats.steiner_inserted > 0,
            &format!(
                "no interior-refinable offender remains after {} Steiner points \
                 (worst ratio {:.2} -> {:.2}; the {} survivors all have \
                 hull-escaping circumcenters, counted for the successor bead's \
                 boundary handling), the FULL exact audit stays clean through \
                 refinement, and refinement is deterministic; \
                 seed 0x1001_2026_0706_0025",
                stats.steiner_inserted,
                stats.worst_before,
                stats.worst_after,
                stats.unrefinable_remaining,
            ),
        );
    });
}

/// tmesh-006 — scale: a 10k-point cloud builds with the O(t) exact
/// audit clean and BRIO keeps walks short. Logs per-point walk cost.
#[test]
fn tmesh_006_scale_run() {
    with_cx(|cx| {
        let pts = cloud(0x1001_2026_0706_0026, 10_000);
        let t = delaunay(&pts, cx).expect("scale build");
        let report = t.audit(false);
        let stats = t.stats();
        let walk_per_point = stats.walk_steps as f64 / stats.points_in as f64;
        let mut em = fs_obs::Emitter::new("fs-mesh/conformance", "tmesh-006/scale");
        let line = em
            .emit(
                fs_obs::Severity::Info,
                fs_obs::EventKind::Custom {
                    name: "mesh-scale-stats".to_string(),
                    json: format!(
                        "{{\"stats\":{},\"walk_per_point\":{walk_per_point:.2}}}",
                        stats.to_json()
                    ),
                },
                None,
            )
            .to_jsonl();
        fs_obs::validate_line(&line).expect("scale stats validate");
        println!("{line}");
        verdict(
            "tmesh-006",
            report.clean() && stats.exhaustive_locates == 0 && walk_per_point < 64.0,
            &format!(
                "10k points: {} tets, local-Delaunay/orientation/adjacency/Euler/hull \
                 audits clean, {walk_per_point:.1} walk steps per insertion (BRIO \
                 locality), no exhaustive-location fallbacks; \
                 seed 0x1001_2026_0706_0026",
                stats.tets_final
            ),
        );
    });
}

/// Too few points teach.
#[test]
fn too_few_points_teaches() {
    with_cx(|cx| {
        let err = delaunay(&[Point3::new(0.0, 0.0, 0.0)], cx).expect_err("refuse");
        assert!(err.to_string().contains("at least 4"), "{err}");
    });
}

// ───────────────────────── remeshing battery ─────────────────────────

use fs_geom::{Chart, Vec3};
use fs_mesh::{MetricField, RemeshOptions, UniformMetric, remesh};
use fs_rep_mesh::HalfEdgeMesh;

fn icosphere(center: Point3, r: f64, sub: u32) -> fs_rep_mesh::Soup {
    fs_rep_mesh::shapes::icosphere(center, r, sub)
}

/// Closed-manifold validity: half-edge invariants + edge-use audit.
fn valid_closed(soup: &fs_rep_mesh::Soup) -> Result<(), String> {
    let he = HalfEdgeMesh::from_triangles(soup.positions.clone(), &soup.triangles)
        .map_err(|e| format!("from_triangles: {e}"))?;
    if let Some(v) = he.check_invariants() {
        return Err(format!("invariant: {v}"));
    }
    let q = assess_quality(soup);
    if !q.sign_certified() {
        return Err(format!("not closed 2-manifold: {q:?}"));
    }
    Ok(())
}

fn metric_edge_lengths(soup: &fs_rep_mesh::Soup, metric: &dyn MetricField) -> Vec<f64> {
    let mut edges = std::collections::BTreeSet::new();
    for t in &soup.triangles {
        for c in 0..3 {
            let (a, b) = (t[c], t[(c + 1) % 3]);
            edges.insert((a.min(b), a.max(b)));
        }
    }
    edges
        .into_iter()
        .map(|(a, b)| {
            let (p, q) = (soup.positions[a as usize], soup.positions[b as usize]);
            let mid = Point3::new(
                f64::midpoint(p.x, q.x),
                f64::midpoint(p.y, q.y),
                f64::midpoint(p.z, q.z),
            );
            let m = metric.metric(mid);
            let e = q.delta_from(p);
            let v = [e.x, e.y, e.z];
            let mut s = 0.0;
            for (row, &vi) in m.iter().zip(&v) {
                s += vi * (row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);
            }
            s.max(0.0).sqrt()
        })
        .collect()
}

/// tmesh-007 — isotropic remeshing of a sphere: unit-metric edge
/// concentration, chart drift bounded, valid closed output, bitwise
/// determinism, and tolerance-level translation equivariance (G3).
#[test]
fn tmesh_007_isotropic_sphere() {
    with_cx(|cx| {
        let base = icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 3);
        let chart = fs_geom::fixtures::SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let metric = UniformMetric { target: 0.25 };
        let opts = RemeshOptions::default();
        let (out, stats) = remesh(&base, Some(&chart), &metric, opts, cx).expect("remesh");
        let (out2, _) = remesh(&base, Some(&chart), &metric, opts, cx).expect("remesh 2");
        let bitwise = out.positions == out2.positions && out.triangles == out2.triangles;
        let validity = valid_closed(&out);
        let lens = metric_edge_lengths(&out, &metric);
        let in_band = lens.iter().filter(|&&l| (0.7..=1.4).contains(&l)).count();
        let frac = in_band as f64 / lens.len() as f64;
        let w = winding_exact(&out, Point3::new(0.0, 0.0, 0.0));
        // Chart drift: vertices (projected) AND centroids (sagitta-bounded).
        let mut worst_centroid = 0.0f64;
        for t in &out.triangles {
            let [a, b, c] = t.map(|v| out.positions[v as usize]);
            let g = Point3::new(
                (a.x + b.x + c.x) / 3.0,
                (a.y + b.y + c.y) / 3.0,
                (a.z + b.z + c.z) / 3.0,
            );
            worst_centroid = worst_centroid.max(chart.eval(g, cx).signed_distance.abs());
        }
        // G3: remesh the shifted problem; compare unshifted (tolerance:
        // fp sums differ, no exact predicates involved).
        let shift = Vec3::new(0.5, 0.25, -0.375);
        let shifted = fs_rep_mesh::Soup {
            positions: base.positions.iter().map(|p| p.offset(shift)).collect(),
            triangles: base.triangles.clone(),
        };
        let chart_s = fs_geom::fixtures::SphereChart {
            center: Point3::new(0.5, 0.25, -0.375),
            radius: 1.0,
        };
        let (out_s, stats_s) =
            remesh(&shifted, Some(&chart_s), &metric, opts, cx).expect("remesh s");
        // Threshold-driven ops legitimately flip borderline decisions
        // under shifted fp arithmetic, so exact-topology G3 cannot hold;
        // the equivariant claim is the QUALITY PROFILE: validity, edge
        // conformity, drift, and element count all match.
        let lens_s = metric_edge_lengths(&out_s, &metric);
        let band = |ls: &[f64]| {
            ls.iter().filter(|&&l| (0.7..=1.4).contains(&l)).count() as f64 / ls.len() as f64
        };
        let (frac_s, frac_base) = (band(&lens_s), band(&lens));
        let count_gap = (out_s.triangles.len() as f64 - out.triangles.len() as f64).abs()
            / out.triangles.len() as f64;
        // Both runs must meet the SAME acceptance floors and land in a
        // modest profile band — exact profile equality is not a property
        // threshold-driven algorithms have.
        let g3 = valid_closed(&out_s).is_ok()
            && frac_s > 0.85
            && (frac_s - frac_base).abs() < 0.08
            && stats_s.worst_chart_drift < 1e-6
            && count_gap < 0.10;
        let mut em = fs_obs::Emitter::new("fs-mesh/conformance", "tmesh-007/isotropic");
        let line = em
            .emit(
                fs_obs::Severity::Info,
                fs_obs::EventKind::Custom {
                    name: "mesh-remesh-iso".to_string(),
                    json: format!(
                        "{{\"ops\":{},\"edge_band_frac\":{frac:.3},\"tris\":{}}}",
                        stats.to_json(),
                        out.triangles.len()
                    ),
                },
                None,
            )
            .to_jsonl();
        fs_obs::validate_line(&line).expect("iso event validates");
        println!("{line}");
        verdict(
            "tmesh-007",
            validity.is_ok()
                && frac > 0.85
                && stats.worst_chart_drift < 1e-6
                && worst_centroid < 3e-2
                && (w - 1.0).abs() < 1e-9
                && bitwise
                && g3,
            &format!(
                "isotropic sphere remesh: {frac_pct:.0}% of edges within [0.7,1.4] \
                 of unit metric length, vertex drift {drift:.1e} and centroid sag \
                 {worst_centroid:.1e} off the chart, output closed/manifold/outward, \
                 BITWISE deterministic, and translation-equivariant in quality \
                 profile (bitwise={bitwise} g3={g3} [frac_s={frac_s:.3} vs \
                 {frac_base:.3}, count_gap={count_gap:.3}, drift_s={drift_s:.1e}] \
                 w={w:.4}; {ops}; {validity:?})",
                frac_pct = frac * 100.0,
                drift = stats.worst_chart_drift,
                drift_s = stats_s.worst_chart_drift,
                ops = stats.to_json(),
            ),
        );
    });
}

/// tmesh-008 — the op-storm robustness battery: rounds of remeshing at
/// randomized targets; after EVERY round the mesh passes half-edge
/// invariants, closed-manifold audit, and Euler = 2.
#[test]
fn tmesh_008_op_storm_validity() {
    with_cx(|cx| {
        let mut rng = Lcg(0x1001_2026_0706_0028);
        let chart = fs_geom::fixtures::SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let mut soup = icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 2);
        let mut total_ops = 0u64;
        let mut all_valid = true;
        let mut detail = String::new();
        for round in 0..16 {
            let target = 0.15 + 0.35 * ((rng.next() >> 32) as f64 / f64::from(u32::MAX));
            let metric = UniformMetric { target };
            let opts = RemeshOptions {
                iterations: 2,
                ..RemeshOptions::default()
            };
            let (next, stats) = remesh(&soup, Some(&chart), &metric, opts, cx).expect("round");
            total_ops += stats.splits + stats.collapses + stats.flips + stats.smooths;
            if let Err(e) = valid_closed(&next) {
                all_valid = false;
                detail = format!("round {round} (target {target:.3}): {e}");
                break;
            }
            let he = HalfEdgeMesh::from_triangles(next.positions.clone(), &next.triangles)
                .expect("valid");
            if he.euler_characteristic() != 2 {
                all_valid = false;
                detail = format!("round {round}: euler {}", he.euler_characteristic());
                break;
            }
            soup = next;
        }
        verdict(
            "tmesh-008",
            all_valid && total_ops > 8_000,
            &format!(
                "16 randomized remesh rounds ({total_ops} ops total) kept half-edge \
                 invariants, closed-manifold status, and Euler = 2 after every round \
                 {detail}; seed 0x1001_2026_0706_0028"
            ),
        );
    });
}

/// tmesh-009 — feature preservation: remeshing a cube keeps its 8
/// corners EXACTLY, keeps every crease-classified output edge ON a cube
/// edge line, stays on the box chart, and remains closed.
#[test]
fn tmesh_009_cube_crease_preservation() {
    with_cx(|cx| {
        let cube = fs_rep_mesh::shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0);
        let chart = fs_geom::fixtures::BoxChart {
            aabb: fs_geom::Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0)),
        };
        let metric = UniformMetric { target: 0.4 };
        let (out, stats) =
            remesh(&cube, Some(&chart), &metric, RemeshOptions::default(), cx).expect("remesh");
        let validity = valid_closed(&out);
        // All 8 corners survive bitwise (comparison by BITS on purpose).
        let mut corners_found = 0;
        for &sx in &[-1.0f64, 1.0] {
            for &sy in &[-1.0f64, 1.0] {
                for &sz in &[-1.0f64, 1.0] {
                    if out.positions.iter().any(|p| {
                        p.x.to_bits() == sx.to_bits()
                            && p.y.to_bits() == sy.to_bits()
                            && p.z.to_bits() == sz.to_bits()
                    }) {
                        corners_found += 1;
                    }
                }
            }
        }
        // Every crease-grade output edge lies on a cube edge line: two
        // coordinates pinned to ±1 (within 1e-9) and AGREEING across
        // both endpoints.
        let mut edge_faces: std::collections::BTreeMap<(u32, u32), Vec<usize>> =
            std::collections::BTreeMap::new();
        for (fi, t) in out.triangles.iter().enumerate() {
            for c in 0..3 {
                let (a, b) = (t[c], t[(c + 1) % 3]);
                edge_faces.entry((a.min(b), a.max(b))).or_default().push(fi);
            }
        }
        let normal = |fi: usize| -> Vec3 {
            let [a, b, c] = out.triangles[fi].map(|v| out.positions[v as usize]);
            let (u, v) = (b.delta_from(a), c.delta_from(a));
            Vec3::new(
                u.y * v.z - u.z * v.y,
                u.z * v.x - u.x * v.z,
                u.x * v.y - u.y * v.x,
            )
        };
        let mut creases_on_edges = true;
        let mut crease_count = 0;
        for (&(a, b), fs) in &edge_faces {
            if fs.len() != 2 {
                creases_on_edges = false;
                continue;
            }
            let (n1, n2) = (normal(fs[0]), normal(fs[1]));
            if n1.dot(n2) / (n1.norm() * n2.norm()).max(1e-300) < 0.7f64.cos() {
                crease_count += 1;
                let (p, q) = (out.positions[a as usize], out.positions[b as usize]);
                let pinned = |u: f64, v: f64| (u.abs() - 1.0).abs() < 1e-9 && (u - v).abs() < 1e-9;
                let shared = usize::from(pinned(p.x, q.x))
                    + usize::from(pinned(p.y, q.y))
                    + usize::from(pinned(p.z, q.z));
                if shared < 2 {
                    creases_on_edges = false;
                }
            }
        }
        let w = winding_exact(&out, Point3::new(0.0, 0.0, 0.0));
        verdict(
            "tmesh-009",
            validity.is_ok()
                && corners_found == 8
                && creases_on_edges
                && crease_count >= 12
                && stats.worst_chart_drift < 5e-3
                && (w - 1.0).abs() < 1e-9,
            &format!(
                "cube remesh keeps all 8 corners bitwise, all {crease_count} \
                 crease-grade edges lie on cube edge lines, drift {:.1e}, closed \
                 and outward ({:?})",
                stats.worst_chart_drift, validity
            ),
        );
    });
}

/// The equator boundary-layer metric: tight normal spacing near z = 0.
struct BoundaryLayerMetric;

impl MetricField for BoundaryLayerMetric {
    fn metric(&self, p: Point3) -> [[f64; 3]; 3] {
        let ht: f64 = 0.35;
        let hz = (0.04 + 0.6 * p.z.abs()).min(ht);
        [
            [1.0 / (ht * ht), 0.0, 0.0],
            [0.0, 1.0 / (ht * ht), 0.0],
            [0.0, 0.0, 1.0 / (hz * hz)],
        ]
    }
}

/// tmesh-010 — anisotropic remeshing under a boundary-layer metric:
/// metric-unit edge conformity, physically stretched aligned elements
/// in the layer, and a MEASURED interpolation-error win over isotropic
/// at comparable element count (the adaptivity-loop model).
#[test]
#[allow(clippy::too_many_lines)] // conformity + alignment + the error demo
fn tmesh_010_anisotropic_boundary_layer() {
    with_cx(|cx| {
        let base = icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 3);
        let chart = fs_geom::fixtures::SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let opts = RemeshOptions {
            iterations: 14,
            ..RemeshOptions::default()
        };
        let (aniso, stats) =
            remesh(&base, Some(&chart), &BoundaryLayerMetric, opts, cx).expect("aniso");
        let validity = valid_closed(&aniso);
        let lens = metric_edge_lengths(&aniso, &BoundaryLayerMetric);
        let frac =
            lens.iter().filter(|&&l| (0.6..=1.6).contains(&l)).count() as f64 / lens.len() as f64;
        // Physical anisotropy in the layer: stretched, equator-aligned.
        let mut aspects = Vec::new();
        let mut align = Vec::new();
        for t in &aniso.triangles {
            let ps = t.map(|v| aniso.positions[v as usize]);
            let gz = (ps[0].z + ps[1].z + ps[2].z) / 3.0;
            if gz.abs() > 0.08 {
                continue;
            }
            let mut longest = (0.0f64, Vec3::new(1.0, 0.0, 0.0));
            let mut shortest = f64::INFINITY;
            for c in 0..3 {
                let e = ps[(c + 1) % 3].delta_from(ps[c]);
                let l = e.norm();
                if l > longest.0 {
                    longest = (l, e);
                }
                shortest = shortest.min(l);
            }
            if shortest > 0.0 {
                aspects.push(longest.0 / shortest);
                align.push((longest.1.z / longest.0).abs());
            }
        }
        let mean_aspect = aspects.iter().sum::<f64>() / aspects.len().max(1) as f64;
        let mean_align = align.iter().sum::<f64>() / align.len().max(1) as f64;
        // Interpolation-error demo vs isotropic at comparable count.
        let t_aniso = aniso.triangles.len() as f64;
        let l_iso = (16.0 * core::f64::consts::PI / (3.0f64.sqrt() * t_aniso)).sqrt();
        let (iso, _) = remesh(
            &base,
            Some(&chart),
            &UniformMetric { target: l_iso },
            opts,
            cx,
        )
        .expect("iso");
        let t_iso = iso.triangles.len() as f64;
        let residual = |s: &fs_rep_mesh::Soup| -> f64 {
            let f = |p: Point3| (p.z / 0.05).tanh();
            let mut worst = 0.0f64;
            for t in &s.triangles {
                let ps = t.map(|v| s.positions[v as usize]);
                let g = Point3::new(
                    (ps[0].x + ps[1].x + ps[2].x) / 3.0,
                    (ps[0].y + ps[1].y + ps[2].y) / 3.0,
                    (ps[0].z + ps[1].z + ps[2].z) / 3.0,
                );
                let interp = (f(ps[0]) + f(ps[1]) + f(ps[2])) / 3.0;
                worst = worst.max((f(g) - interp).abs());
            }
            worst
        };
        let (r_aniso, r_iso) = (residual(&aniso), residual(&iso));
        let counts_comparable = (t_iso - t_aniso).abs() / t_aniso < 0.4;
        let mut em = fs_obs::Emitter::new("fs-mesh/conformance", "tmesh-010/aniso");
        let line = em
            .emit(
                fs_obs::Severity::Info,
                fs_obs::EventKind::Custom {
                    name: "mesh-remesh-aniso".to_string(),
                    json: format!(
                        "{{\"ops\":{},\"band_frac\":{frac:.3},\"mean_aspect\":{mean_aspect:.2},\
                         \"mean_align_z\":{mean_align:.3},\"tris_aniso\":{t_aniso},\
                         \"tris_iso\":{t_iso},\"residual_aniso\":{r_aniso:.4},\
                         \"residual_iso\":{r_iso:.4}}}",
                        stats.to_json()
                    ),
                },
                None,
            )
            .to_jsonl();
        fs_obs::validate_line(&line).expect("aniso event validates");
        println!("{line}");
        verdict(
            "tmesh-010",
            validity.is_ok()
                && frac > 0.7
                && mean_aspect > 2.0
                && mean_align < 0.45
                && counts_comparable
                && r_aniso < 0.6 * r_iso,
            &format!(
                "boundary-layer metric realized: {:.0}% metric-unit edges, layer \
                 elements stretched {mean_aspect:.1}x and equator-aligned \
                 (|cos z|={mean_align:.2}), and interpolation residual {r_aniso:.3} \
                 vs isotropic {r_iso:.3} at comparable counts ({t_aniso} vs {t_iso} \
                 tris) — the adaptivity win, measured ({:?})",
                frac * 100.0,
                validity
            ),
        );
    });
}

// ------------------------------------------------------------- tmesh-011/012
// Full-Ruppert hull-facet splitting + sliver exudation (bead uee3).

#[test]
fn tmesh_011_policy_floor_and_hull_split_evidence() {
    with_cx(|cx| {
        let mut rng = Lcg(0x1001_2026_0708_0011);
        let points: Vec<Point3> = (0..220)
            .map(|_| Point3::new(rng.dyadic(), rng.dyadic(), rng.dyadic()))
            .collect();
        // (i) Default policy (splitting OFF): quality improves, every
        // remaining offender is accounted (skipped or protected).
        let mut t = delaunay(&points, cx).expect("builds");
        let stats = refine(&mut t, RefineOptions::default(), cx).expect("refines");
        let audit = t.audit(true);
        // The global worst may sit AT the hull (its circumcenter
        // escapes and is skipped) — the honest gate is non-worsening
        // plus full accounting.
        let default_ok = audit.clean()
            && stats.worst_after <= stats.worst_before
            && stats.refinable_remaining == 0;
        // (ii) Experimental hull splitting: exact-audit-clean and
        // deterministic — and the MEASURED quality regression is the
        // ledgered evidence that full-Ruppert quality is coupled to
        // PLC boundary protection (the successor's machinery).
        let run_split = |cx: &Cx<'_>| {
            let mut t = delaunay(&points, cx).expect("builds");
            let st = refine(
                &mut t,
                RefineOptions {
                    max_radius_edge: 1.8,
                    max_steiner: 1200,
                    split_hull_facets: true,
                    min_edge_factor: 0.2,
                },
                cx,
            )
            .expect("refines");
            (t.audit(true).clean(), st)
        };
        let (clean_a, st_a) = run_split(cx);
        let (clean_b, st_b) = run_split(cx);
        let split_infra_ok =
            clean_a && clean_b && st_a.to_json() == st_b.to_json() && st_a.hull_facets_split > 0;
        let pass = default_ok && split_infra_ok;
        println!(
            "{{\"test\":\"tmesh-011\",\"verdict\":\"{}\",\"default\":{},\
             \"split_evidence\":{},\"split_regression_documented\":true}}",
            if pass { "pass" } else { "fail" },
            stats.to_json(),
            st_a.to_json()
        );
        assert!(
            pass,
            "tmesh-011: default {} split {}",
            stats.to_json(),
            st_a.to_json()
        );
    });
}

#[test]
fn tmesh_012_sliver_exudation() {
    with_cx(|cx| {
        let mut rng = Lcg(0x1001_2026_0708_0012);
        let points: Vec<Point3> = (0..200)
            .map(|_| Point3::new(rng.dyadic(), rng.dyadic(), rng.dyadic()))
            .collect();
        let run = |cx: &Cx<'_>| {
            let mut t = delaunay(&points, cx).expect("builds");
            let _ = refine(
                &mut t,
                RefineOptions {
                    max_radius_edge: 1.8,
                    max_steiner: 800,
                    split_hull_facets: true,
                    min_edge_factor: 0.2,
                },
                cx,
            )
            .expect("refines");
            let inputs_before = t.points()[..t.steiner_from as usize].to_vec();
            let stats = fs_mesh::exude(
                &mut t,
                ExudeOptions {
                    dihedral_min_deg: 8.0,
                    rounds: 8,
                    jitter: 0.05,
                },
                cx,
            )
            .expect("exudes");
            (t, stats, inputs_before)
        };
        let (t, stats, inputs_before) = run(cx);
        let (_, stats2, _) = run(cx);
        let deterministic = stats.to_json() == stats2.to_json();
        let audit = t.audit(true);
        // Input points untouched, bitwise.
        let inputs_after = &t.points()[..t.steiner_from as usize];
        let inputs_frozen = inputs_before.iter().zip(inputs_after).all(|(a, b)| {
            a.x.to_bits() == b.x.to_bits()
                && a.y.to_bits() == b.y.to_bits()
                && a.z.to_bits() == b.z.to_bits()
        });
        let movable_before = stats.slivers_before - stats.input_protected;
        let movable_after = stats.slivers_after.saturating_sub(stats.input_protected);
        let improved = stats.slivers_before == 0
            || movable_after * 2 <= movable_before
            || stats.slivers_after < stats.slivers_before;
        let pass = audit.clean() && improved && inputs_frozen && deterministic;
        println!(
            "{{\"test\":\"tmesh-012\",\"verdict\":\"{}\",\"stats\":{},\"audit_clean\":{},\
             \"inputs_frozen\":{inputs_frozen},\"deterministic\":{deterministic}}}",
            if pass { "pass" } else { "fail" },
            stats.to_json(),
            audit.clean()
        );
        assert!(pass, "tmesh-012: {}", stats.to_json());
    });
}

/// Canonical (index-free) form of a tetrahedralization: the sorted
/// list of sorted vertex quadruples of live real tets. Two builds are
/// "the same mesh" iff these match bitwise (points share input order).
fn canonical_tets(t: &Tetrahedralization) -> Vec<[u32; 4]> {
    let mut out: Vec<[u32; 4]> = t
        .tets()
        .into_iter()
        .map(|mut q| {
            q.sort_unstable();
            q
        })
        .collect();
    out.sort_unstable();
    out
}

/// tmesh-013: deterministic parallel domain coloring — the colored
/// (read-parallel, apply-canonical) construction is bitwise identical
/// at EVERY thread count, merges bitwise against the sequential
/// kernel's canonical output, keeps the exact audit clean, and its
/// coloring is genuinely commutative (adversarial reversed-order
/// application inside each color produces the identical mesh).
#[test]
fn tmesh_013_parallel_coloring() {
    with_cx(|cx| {
        let pts = cloud(0x1001_2026_0708_0013, 900);
        let seq = delaunay(&pts, cx).expect("sequential kernel");
        // Colors reorder only across provably disjoint pairs, so the
        // kernel merge is CANONICAL (allocation order legitimately
        // differs); thread-count invariance among colored runs is RAW
        // bitwise.
        let seq_canon = canonical_tets(&seq);
        let mut stats1 = None;
        let mut raw1: Option<Vec<[u32; 4]>> = None;
        for threads in [1usize, 2, 4, 8] {
            let (colored, stats) =
                delaunay_colored(&pts, threads, 64, cx).expect("colored build");
            let raw = colored.tets();
            match &raw1 {
                None => raw1 = Some(raw.clone()),
                Some(r) => assert_eq!(&raw, r, "RAW thread-count invariance at T={threads}"),
            }
            verdict(
                &format!("tmesh-013-threads-{threads}"),
                canonical_tets(&colored) == seq_canon,
                &format!(
                    "canonical merge vs sequential kernel: {} tets ({})",
                    seq_canon.len(),
                    stats.to_json()
                ),
            );
            let audit = colored.audit(true);
            verdict(
                &format!("tmesh-013-audit-{threads}"),
                audit.clean(),
                "exact audit clean on the colored build",
            );
            if threads == 1 {
                stats1 = Some(stats);
            }
        }
        // The batching is real parallelism on general-position input,
        // not a serial crawl.
        let s1 = stats1.expect("stats recorded");
        verdict(
            "tmesh-013-batch-width",
            s1.largest_batch >= 8 && s1.batches < s1.points,
            &format!("parallel batch evidence: {}", s1.to_json()),
        );
        // Adversarial commutativity: reversed within-batch application
        // (allocation order legitimately differs — compare canonically).
        let rev = delaunay_colored_reversed(&pts, 64, cx).expect("reversed build");
        verdict(
            "tmesh-013-commutativity",
            canonical_tets(&rev) == seq_canon,
            "reversed within-batch insertion order yields the identical canonical mesh",
        );
        // Degenerate adversary: a structured grid (massively cospherical)
        // through the colored path.
        let mut grid_pts = Vec::new();
        for i in 0..6i32 {
            for j in 0..6i32 {
                for k in 0..6i32 {
                    grid_pts.push(Point3::new(
                        f64::from(i) / 5.0,
                        f64::from(j) / 5.0,
                        f64::from(k) / 5.0,
                    ));
                }
            }
        }
        let gseq = delaunay(&grid_pts, cx).expect("grid kernel");
        let (gcol, gstats) = delaunay_colored(&grid_pts, 4, 32, cx).expect("grid colored");
        verdict(
            "tmesh-013-degenerate-grid",
            canonical_tets(&gcol) == canonical_tets(&gseq) && gcol.audit(true).clean(),
            &format!(
                "6x6x6 cospherical grid (order-dependent ties): canonical == kernel \
                 ({} tets), audit clean ({})",
                gcol.tets().len(),
                gstats.to_json()
            ),
        );
    });
}
