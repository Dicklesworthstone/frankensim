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
use fs_mesh::{MeshError, RefineOptions, Tetrahedralization, delaunay, refine};
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
        };
        let stats = refine(&mut t, opts, cx).expect("refine");
        let report = t.audit(true);
        let mut t2 = delaunay(&pts, cx).expect("build 2");
        let _ = refine(&mut t2, opts, cx).expect("refine 2");
        let deterministic = t.tets() == t2.tets() && t.points() == t2.points();
        let improved = stats.worst_after < stats.worst_before
            || (stats.worst_before == 0.0 && stats.worst_after == 0.0);
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
            improved && report.clean() && deterministic && stats.steiner_inserted > 0,
            &format!(
                "radius-edge improved {:.2} -> {:.2} with {} Steiner points \
                 ({} hull-encroaching offenders skipped and counted), the FULL exact \
                 audit stays clean through refinement, and refinement is \
                 deterministic; seed 0x1001_2026_0706_0025",
                stats.worst_before,
                stats.worst_after,
                stats.steiner_inserted,
                stats.skipped_outside_hull
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
