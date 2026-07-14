//! fs-rep-mesh conformance suite (CONTRACT.md: any reimplementation must
//! pass). Half-edge invariants under random edits, point-triangle
//! distance vs brute force, winding classification on nightmare soup,
//! dipole-vs-exact error, the repair battery with receipts, δδ = 0, and
//! watertight rays. JSON-line verdicts; seeded cases carry seeds.

use asupersync::types::Budget;
use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{
    Aabb, Chart, ChartSample, Differentiability, Point3, SamplingDomainError, TraceStepClaim, Vec3,
};
use fs_rep_mesh::{
    BracketCertificateError, BracketEvidenceIssue, ContourError, DC_MAX_CELLS_PER_AXIS, DcOptions,
    HalfEdgeMesh, MeshChart, TetComplex, WindingOctree, bracket_certificate, dual_contour,
    dual_contour_clipped, point_triangle_distance, ray_triangle_watertight, repair, shapes,
    winding_exact,
};
use std::sync::atomic::{AtomicU64, Ordering};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-rep-mesh/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\
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

    fn unit(&mut self) -> f64 {
        ((self.next() >> 11) as f64) / (1u64 << 53) as f64
    }

    fn below(&mut self, n: u64) -> u64 {
        (self.next() >> 32) % n
    }
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    with_gate_cx(&gate, f)
}

fn with_gate_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x9E54,
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

struct CountingPlane {
    x_offset: f64,
    evals: AtomicU64,
}

impl CountingPlane {
    fn new(x_offset: f64) -> Self {
        Self {
            x_offset,
            evals: AtomicU64::new(0),
        }
    }

    fn eval_count(&self) -> u64 {
        self.evals.load(Ordering::Relaxed)
    }
}

impl Chart for CountingPlane {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        self.evals.fetch_add(1, Ordering::Relaxed);
        let signed_distance = x.x - self.x_offset;
        ChartSample {
            signed_distance,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error: NumericalCertificate::enclosure(signed_distance, signed_distance),
        }
    }

    fn support(&self) -> Aabb {
        Aabb {
            min: Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
            max: Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
        }
    }

    fn name(&self) -> &'static str {
        "test/counting-plane"
    }

    fn differentiability(&self) -> Differentiability {
        Differentiability::C0
    }
}

#[derive(Clone, Copy)]
enum BracketAuthorityMode<'a> {
    NoClaim,
    Estimate,
    NoClaimEvidence,
    MalformedExact,
    RequestCancellation(&'a CancelGate),
}

struct BracketAuthorityPlane<'a> {
    mode: BracketAuthorityMode<'a>,
    evals: AtomicU64,
}

impl BracketAuthorityPlane<'_> {
    fn eval_count(&self) -> u64 {
        self.evals.load(Ordering::Relaxed)
    }
}

impl Chart for BracketAuthorityPlane<'_> {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let previous = self.evals.fetch_add(1, Ordering::Relaxed);
        if previous == 0
            && let BracketAuthorityMode::RequestCancellation(gate) = self.mode
        {
            gate.request();
        }
        let error = match self.mode {
            BracketAuthorityMode::Estimate => NumericalCertificate::estimate(x.x, x.x),
            BracketAuthorityMode::NoClaimEvidence => NumericalCertificate::no_claim(),
            BracketAuthorityMode::MalformedExact => NumericalCertificate {
                kind: NumericalKind::Exact,
                lo: x.x,
                hi: x.x + 1.0,
            },
            BracketAuthorityMode::NoClaim | BracketAuthorityMode::RequestCancellation(_) => {
                NumericalCertificate::enclosure(x.x, x.x)
            }
        };
        ChartSample {
            signed_distance: x.x,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error,
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        match self.mode {
            BracketAuthorityMode::NoClaim => TraceStepClaim::NoClaim,
            BracketAuthorityMode::Estimate
            | BracketAuthorityMode::NoClaimEvidence
            | BracketAuthorityMode::MalformedExact
            | BracketAuthorityMode::RequestCancellation(_) => TraceStepClaim::ExactDistance,
        }
    }

    fn name(&self) -> &'static str {
        "test/bracket-authority-plane"
    }
}

#[test]
fn rmesh_001_halfedge_invariants_survive_random_flip_batteries() {
    const SEED: u64 = 0x1001_2026_0706_0001;
    let soup = shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 2);
    let mut mesh =
        HalfEdgeMesh::from_triangles(soup.positions.clone(), &soup.triangles).expect("manifold");
    assert_eq!(mesh.euler_characteristic(), 2, "sphere: V - E + F = 2");
    let mut rng = Lcg(SEED);
    let mut flips_done = 0u32;
    for _ in 0..2_000 {
        let h = rng.below(mesh.half_edges.len() as u64) as u32;
        if mesh.flip_edge(h) {
            flips_done += 1;
        }
        if let Some(violation) = mesh.check_invariants() {
            verdict("rmesh-001", false, &format!("invariant broke: {violation}"));
        }
    }
    verdict(
        "rmesh-001",
        flips_done > 200 && mesh.euler_characteristic() == 2,
        &format!(
            "half-edge invariants held through {flips_done} random edge flips (seed {SEED:#x}); \
             Euler characteristic still 2"
        ),
    );
}

#[test]
fn rmesh_002_point_triangle_distance_matches_brute_force_and_chart_laws() {
    const SEED: u64 = 0x1002_2026_0706_D157;
    let mut rng = Lcg(SEED);
    // Distance vs dense barycentric sampling on random triangles.
    let mut worst = 0.0f64;
    for _ in 0..300 {
        let rp = |rng: &mut Lcg| {
            Point3::new(
                (rng.unit() - 0.5) * 4.0,
                (rng.unit() - 0.5) * 4.0,
                (rng.unit() - 0.5) * 4.0,
            )
        };
        let (a, b, c, p) = (rp(&mut rng), rp(&mut rng), rp(&mut rng), rp(&mut rng));
        let fast = point_triangle_distance(p, a, b, c);
        let mut brute = f64::INFINITY;
        let n = 60;
        for i in 0..=n {
            for j in 0..=(n - i) {
                let (u, v) = (f64::from(i) / f64::from(n), f64::from(j) / f64::from(n));
                let w = 1.0 - u - v;
                let q = Point3::new(
                    a.x * w + b.x * u + c.x * v,
                    a.y * w + b.y * u + c.y * v,
                    a.z * w + b.z * u + c.z * v,
                );
                brute = brute.min(p.delta_from(q).norm());
            }
        }
        worst = worst.max((fast - brute).abs());
        assert!(
            fast <= brute + 1e-12,
            "exact distance can never exceed a sampled distance"
        );
    }
    // Fixture observations on the icosphere: sd within the mesh's
    // approximation band of the analytic sphere; inside ⇔ sd < 0; measured
    // 1-Lipschitz behavior. rmesh-002c separately locks that raw soup does not
    // advertise this fixture observation as a theorem.
    let (band_ok, lip_ok) = with_cx(|cx| {
        let soup = shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 3);
        let chart = MeshChart::new(soup);
        let mut rng = Lcg(SEED ^ 0xC047);
        let mut band_ok = true;
        let mut lip_ok = true;
        for _ in 0..500 {
            let p = Point3::new(
                (rng.unit() - 0.5) * 3.0,
                (rng.unit() - 0.5) * 3.0,
                (rng.unit() - 0.5) * 3.0,
            );
            let s = chart.eval(p, cx);
            let analytic = p.delta_from(Point3::new(0.0, 0.0, 0.0)).norm() - 1.0;
            // Icosphere-3 chord error is < 6e-3.
            band_ok &= (s.signed_distance - analytic).abs() < 8e-3;
            band_ok &= chart.inside(p, cx) == (s.signed_distance < 0.0);
            let q = p.offset(Vec3::new(rng.unit() * 0.2, rng.unit() * 0.2, 0.0));
            let sq = chart.eval(q, cx);
            lip_ok &=
                (sq.signed_distance - s.signed_distance).abs() <= q.delta_from(p).norm() + 1e-9;
        }
        (band_ok, lip_ok)
    });
    verdict(
        "rmesh-002",
        worst < 0.02 && band_ok && lip_ok,
        &format!(
            "exact point-triangle distance under-approximates 1830-sample brute force by at \
             most sampling gap (worst {worst:.4}); chart tracks the analytic sphere and is \
             observed 1-Lipschitz on this fixture without advertising a generic theorem \
             (seed {SEED:#x})"
        ),
    );
}

#[test]
fn rmesh_002b_degenerate_triangles_yield_finite_correct_distances() {
    // Regression: on a degenerate triangle (repeated vertex or collinear
    // vertices) Ericson's edge-region parameter divides by |ab|²/|ac|²/|bc|²
    // or the area — all zero here — so e.g. the edge-AB branch's `d1/(d1-d3)`
    // was `0/0 = NaN`. The CONTRACT promises degenerate triangles yield
    // well-defined distances. The collapse is a segment/point, so the answer
    // is the nearest edge distance; verify it is finite AND matches brute force.
    const SEED: u64 = 0x1002_B026_0706_DE9E;
    let mut rng = Lcg(SEED);
    let rp = |rng: &mut Lcg| {
        Point3::new(
            (rng.unit() - 0.5) * 4.0,
            (rng.unit() - 0.5) * 4.0,
            (rng.unit() - 0.5) * 4.0,
        )
    };
    let brute = |p: Point3, a: Point3, b: Point3, c: Point3| {
        let n = 60;
        let mut best = f64::INFINITY;
        for i in 0..=n {
            for j in 0..=(n - i) {
                let (u, v) = (f64::from(i) / f64::from(n), f64::from(j) / f64::from(n));
                let w = 1.0 - u - v;
                let q = Point3::new(
                    a.x * w + b.x * u + c.x * v,
                    a.y * w + b.y * u + c.y * v,
                    a.z * w + b.z * u + c.z * v,
                );
                best = best.min(p.delta_from(q).norm());
            }
        }
        best
    };
    let mut worst = 0.0f64;
    let mut all_finite = true;
    for _ in 0..400 {
        let (a, c, p) = (rp(&mut rng), rp(&mut rng), rp(&mut rng));
        let d = Vec3::new(rng.unit() - 0.5, rng.unit() - 0.5, rng.unit() - 0.5);
        // Every family the guard must handle: repeated a=b, repeated b=c,
        // collinear (a, a+d, a+2d), and fully collapsed a=b=c.
        let cases = [
            (a, a, c),                                // a == b
            (a, c, c),                                // b == c
            (a, a.offset(d), a.offset(d.scale(2.0))), // collinear
            (a, a, a),                                // single point
        ];
        for (ta, tb, tc) in cases {
            let fast = point_triangle_distance(p, ta, tb, tc);
            all_finite &= fast.is_finite();
            let bf = brute(p, ta, tb, tc);
            worst = worst.max((fast - bf).abs());
            assert!(
                fast <= bf + 1e-12,
                "exact distance can never exceed a sampled distance (degenerate)"
            );
        }
    }
    // The exact repro that used to NaN: a == b collapses to segment a–c; the
    // closest point to (0.5,1,0) is (0.5,0,0), so the distance is exactly 1.
    let repro = point_triangle_distance(
        Point3::new(0.5, 1.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    );
    verdict(
        "rmesh-002b",
        all_finite && worst < 0.02 && (repro - 1.0).abs() < 1e-12,
        &format!(
            "degenerate triangles (repeated/collinear/collapsed) give finite distances that \
             match brute force (worst {worst:.4}); the a==b repro is {repro:.6} (want 1.0), \
             not NaN (seed {SEED:#x})"
        ),
    );
}

#[test]
fn rmesh_002c_raw_mesh_chart_never_promotes_clean_soup_authority() {
    with_cx(|cx| {
        let clean = shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0);
        let output = clean.clone();
        let chart = MeshChart::new(clean);

        assert_eq!(chart.trace_step_claim(), TraceStepClaim::NoClaim);
        let finite = chart.eval(Point3::new(2.0, 0.0, 0.0), cx);
        assert_eq!(finite.lipschitz, None);
        assert_eq!(finite.error.kind, NumericalKind::Estimate);
        assert!(finite.error.lo.is_finite() && finite.error.hi.is_finite());

        let nonfinite = chart.eval(Point3::new(f64::INFINITY, 0.0, 0.0), cx);
        assert_eq!(nonfinite.lipschitz, None);
        assert_eq!(nonfinite.error.kind, NumericalKind::NoClaim);

        assert_eq!(
            bracket_certificate(&chart, &output, 0.25, cx),
            Err(BracketCertificateError::UnsupportedTraceClaim {
                actual: TraceStepClaim::NoClaim,
            })
        );
    });
    verdict(
        "rmesh-002c",
        true,
        "a raw clean-looking closed soup remains TraceStepClaim::NoClaim with no Lipschitz bound; finite samples are Estimate and non-finite samples are NoClaim",
    );
}

#[test]
fn rmesh_003_winding_classifies_nightmare_soup() {
    const SEED: u64 = 0x1003_2026_0706_50FA;
    let mut rng = Lcg(SEED);
    // Nightmare corpus: icosphere with duplicates, degenerates, a flipped
    // patch, and a punched hole.
    let clean = shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 2);
    let soup = shapes::corrupt(clean, 6, 4, 40..48, Some(10));
    let octree = WindingOctree::build(&soup, 2.0);
    let mut correct = 0u32;
    let mut total = 0u32;
    for _ in 0..2_000 {
        // Sample away from the surface and the defects' shadow (|r-1| > 0.2).
        let dir = Vec3::new(rng.unit() - 0.5, rng.unit() - 0.5, rng.unit() - 0.5);
        let n = dir.norm().max(1e-9);
        let r = if rng.below(2) == 0 {
            0.3 + rng.unit() * 0.45
        } else {
            1.25 + rng.unit() * 1.0
        };
        let p = Point3::new(dir.x / n * r, dir.y / n * r, dir.z / n * r);
        let truly_inside = r < 1.0;
        if octree.inside(&soup, p) == truly_inside {
            correct += 1;
        }
        total += 1;
    }
    let rate = f64::from(correct) / f64::from(total);
    verdict(
        "rmesh-003",
        rate > 0.99,
        &format!(
            "winding classification on the nightmare soup (dups+degens+flipped patch+hole): \
             {correct}/{total} correct away from defects (seed {SEED:#x})"
        ),
    );
}

#[test]
fn rmesh_004_dipole_approximation_tracks_exact_within_declared_error() {
    const SEED: u64 = 0x1004_2026_0706_D1B0;
    let soup = shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 3);
    let octree = WindingOctree::build(&soup, 2.0);
    let mut rng = Lcg(SEED);
    let mut worst = 0.0f64;
    for _ in 0..400 {
        let p = Point3::new(
            (rng.unit() - 0.5) * 5.0,
            (rng.unit() - 0.5) * 5.0,
            (rng.unit() - 0.5) * 5.0,
        );
        // Skip the surface shell where winding is legitimately steep.
        let r = p.delta_from(Point3::new(0.0, 0.0, 0.0)).norm();
        if (r - 1.0).abs() < 0.1 {
            continue;
        }
        let exact = winding_exact(&soup, p);
        let approx = octree.winding(&soup, p);
        worst = worst.max((approx - exact).abs());
    }
    let mut em = fs_obs::Emitter::new("fs-rep-mesh/conformance", "rmesh-004/dipole");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "rep-mesh-dipole-error".to_string(),
                json: format!("{{\"worst_abs_error\":{worst:.6},\"beta\":2.0}}"),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("dipole event validates");
    println!("{line}");
    verdict(
        "rmesh-004",
        worst < 0.05,
        &format!(
            "dipole octree (beta=2) tracks exact winding within {worst:.4} off-surface \
             (seed {SEED:#x}; error ledgered)"
        ),
    );
}

#[test]
fn rmesh_005_repair_battery_heals_the_corpus_with_receipts() {
    let clean = shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 1);
    let n_clean = clean.triangles.len();
    let corrupted = shapes::corrupt(clean, 5, 3, 12..20, Some(11));
    let outcome = repair(corrupted, 8);
    // Receipts cover every defect class.
    let classes: std::collections::BTreeSet<&str> =
        outcome.receipts.iter().map(|r| r.defect).collect();
    let all_classes = [
        "boundary-hole",
        "degenerate-face",
        "duplicate-face",
        "flipped-patch",
    ]
    .iter()
    .all(|c| classes.contains(c));
    // Healed: face count restored, winding at center back to ~1, and the
    // mesh is manifold again (half-edge build succeeds).
    let healed_count = outcome.soup.triangles.len() == n_clean;
    let w = winding_exact(&outcome.soup, Point3::new(0.0, 0.0, 0.0));
    let winding_ok = (w - 1.0).abs() < 1e-6;
    let manifold =
        HalfEdgeMesh::from_triangles(outcome.soup.positions.clone(), &outcome.soup.triangles)
            .is_ok();
    let mut em = fs_obs::Emitter::new("fs-rep-mesh/conformance", "rmesh-005/repair");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "rep-mesh-repair-receipts".to_string(),
                json: format!("{{\"receipts\":{}}}", outcome.receipts_json()),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("receipts validate");
    println!("{line}");
    verdict(
        "rmesh-005",
        all_classes && healed_count && winding_ok && manifold,
        &format!(
            "dups+degens+flips+hole all repaired with receipts ({} actions); center winding \
             restored to {w:.6}; half-edge build succeeds again",
            outcome.receipts.len()
        ),
    );
}

#[test]
fn rmesh_006_incidence_satisfies_dd_zero_and_rays_are_watertight() {
    const SEED: u64 = 0x1006_2026_0706_DD00;
    // δδ = 0 exactly, on a 5-tet cube decomposition and a random fan.
    let five_tet = TetComplex::from_tets(
        8,
        vec![
            [0, 1, 2, 5],
            [0, 2, 3, 7],
            [0, 5, 7, 4],
            [2, 5, 7, 6],
            [0, 2, 5, 7],
        ],
    );
    let mut rng = Lcg(SEED);
    let mut fan = Vec::new();
    for i in 0..12u32 {
        fan.push([0u32, i + 1, i + 2, 14 + (rng.below(4) as u32)]);
    }
    let fan_complex = TetComplex::from_tets(20, fan);
    let mut dd_zero = true;
    for complex in [&five_tet, &fan_complex] {
        let (d0, d1, d2) = (complex.d0(), complex.d1(), complex.d2());
        for probe in 0..complex.vertex_count.min(6) {
            let mut x = vec![0i64; complex.vertex_count];
            x[probe] = 1;
            let dd = d1.apply(&d0.apply(&x));
            dd_zero &= dd.iter().all(|&v| v == 0);
        }
        for probe in 0..complex.edges.len().min(8) {
            let mut x = vec![0i64; complex.edges.len()];
            x[probe] = 1;
            let dd = d2.apply(&d1.apply(&x));
            dd_zero &= dd.iter().all(|&v| v == 0);
        }
    }
    // Watertight rays: through a cube, axis rays hit exactly twice
    // (entry+exit) even when passing through face diagonals shared by two
    // triangles; and the MeshChart raycast agrees with analytic t.
    let cube = shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0);
    let count_hits = |origin: Point3, dir: Vec3| -> usize {
        let mut hits = 0;
        for t in 0..cube.triangles.len() {
            let [a, b, c] = cube.tri(t);
            if ray_triangle_watertight(origin, dir, a, b, c).is_some() {
                hits += 1;
            }
        }
        hits
    };
    // Through the center: crosses the diagonal edge on both faces —
    // watertight means never a LEAK (>= entry + exit; exact-edge hits may
    // double-count, documented).
    let through_center = count_hits(Point3::new(-3.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
    let through_offcenter = count_hits(Point3::new(-3.0, 0.3, 0.17), Vec3::new(1.0, 0.0, 0.0));
    let no_leak = through_center >= 2 && through_offcenter == 2;
    let chart_ray_ok = {
        let chart = MeshChart::new(cube);
        let t = chart.raycast(Point3::new(-3.0, 0.3, 0.17), Vec3::new(1.0, 0.0, 0.0), 10.0);
        t.is_some_and(|t| (t - 2.0).abs() < 1e-9)
    };
    verdict(
        "rmesh-006",
        dd_zero && no_leak && chart_ray_ok,
        &format!(
            "d1∘d0 = 0 and d2∘d1 = 0 EXACTLY on 5-tet cube + random fan (seed {SEED:#x}); \
             axis rays never leak through shared edges (center {through_center} hits, \
             off-center {through_offcenter}); chart raycast hits at analytic t=2"
        ),
    );
}

#[test]
// One end-to-end scenario (build -> equivariance -> incremental -> downgrade);
// splitting would duplicate the expensive fixture builds.
#[allow(clippy::too_many_lines)]
fn rmesh_007_mesh_to_sdf_converter_is_honest_equivariant_and_incremental() {
    const SEED: u64 = 0x1007_2026_0706_C0DE;
    let (analytic_ok, equivariant, incremental_identical, downgrade_ok, samples) = with_cx(|cx| {
        // Measured accuracy: high-res icosphere -> SDF matches the analytic
        // sphere within the receipt's recorded estimate plus the mesh chord
        // band. Generic soup lacks the global validity certificate needed to
        // make that envelope rigorous.
        let ico = shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 3);
        let chart = MeshChart::new(ico.clone());
        let receipt = fs_rep_mesh::mesh_to_sdf(&chart, 0.08, cx).expect("convert");
        assert_eq!(
            receipt.numerical.kind,
            fs_evidence::NumericalKind::Estimate,
            "edge and aggregate-volume screens cannot certify generic soup"
        );
        assert_eq!(
            receipt.value.abstract_distance_kind(),
            fs_evidence::NumericalKind::Estimate
        );
        assert_eq!(
            receipt.qoi.to_bits(),
            receipt
                .value
                .abstract_distance_bound()
                .expect("sampled estimate carries a finite bound")
                .to_bits(),
            "mesh receipt uses the sampled payload's total abstract-distance estimate bound"
        );
        assert!(
            receipt
                .model
                .cards
                .contains(&"winding-sign-heuristic".to_string()),
            "clean-looking soup still names the uncertified sign model"
        );
        let mut rng = Lcg(SEED);
        let mut analytic_ok = true;
        for _ in 0..800 {
            let p = Point3::new(
                (rng.unit() - 0.5) * 2.6,
                (rng.unit() - 0.5) * 2.6,
                (rng.unit() - 0.5) * 2.6,
            );
            let sd = receipt.value.eval(p, cx).signed_distance;
            let analytic = p.delta_from(Point3::new(0.0, 0.0, 0.0)).norm() - 1.0;
            // measured envelope: recorded field estimate + icosphere-3 chord error.
            analytic_ok &= (sd - analytic).abs() <= receipt.qoi + 6e-3 + 1e-9;
        }
        // G3 translation equivariance: translate the mesh AND the queries;
        // the sampled fields must agree bitwise at corresponding points.
        let shift = fs_geom::Vec3::new(0.31, -0.17, 0.23);
        let moved = fs_rep_mesh::Soup {
            positions: ico.positions.iter().map(|p| p.offset(shift)).collect(),
            triangles: ico.triangles.clone(),
        };
        let moved_receipt =
            fs_rep_mesh::mesh_to_sdf(&MeshChart::new(moved), 0.08, cx).expect("convert moved");
        let mut equivariant = true;
        for _ in 0..200 {
            let p = Point3::new(
                (rng.unit() - 0.5) * 2.0,
                (rng.unit() - 0.5) * 2.0,
                (rng.unit() - 0.5) * 2.0,
            );
            let a = receipt.value.eval(p, cx).signed_distance;
            let b = moved_receipt
                .value
                .eval(p.offset(shift), cx)
                .signed_distance;
            // Grids are anchored to supports which translate with the mesh,
            // so samples align and values match to fp noise.
            equivariant &= (a - b).abs() < 1e-6;
        }
        assert!(
            receipt.certified().is_err(),
            "generic mesh receipt must not cross the rigorous-certification boundary"
        );
        // Incremental == full (G5): edit a vertex patch, refresh the dirty
        // box, compare bitwise against a full rebuild of the edited mesh.
        let mut edited = ico.clone();
        for p in edited.positions.iter_mut().take(12) {
            *p = Point3::new(p.x * 1.05, p.y * 1.05, p.z * 1.05);
        }
        let edited_chart = || MeshChart::new(edited.clone());
        let mut inc = fs_rep_mesh::IncrementalMeshSdf::build(MeshChart::new(ico.clone()), 0.08, cx)
            .expect("initial");
        // The dirty box: everything within reach of the moved vertices.
        // Distance fields have GLOBAL support in principle; for a 5% bump
        // the change is confined to the bump's distance cone — cover it
        // generously.
        inc.update(
            edited_chart(),
            fs_geom::Aabb::new(Point3::new(-2.0, -2.0, -2.0), Point3::new(2.0, 2.0, 2.0)),
            cx,
        )
        .expect("update");
        let full = fs_rep_mesh::mesh_to_sdf(&edited_chart(), 0.08, cx).expect("full rebuild");
        let mut incremental_identical = true;
        for _ in 0..400 {
            let p = Point3::new(
                (rng.unit() - 0.5) * 2.4,
                (rng.unit() - 0.5) * 2.4,
                (rng.unit() - 0.5) * 2.4,
            );
            let a = inc.sdf().eval(p, cx).signed_distance;
            let b = full.value.eval(p, cx).signed_distance;
            incremental_identical &= a.to_bits() == b.to_bits();
        }
        // Adversarial open/slivered soup also remains Estimate and names the
        // heuristic plus failed quality-screen diagnostics.
        let nasty = shapes::corrupt(
            shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 1),
            0,
            2,
            0..0,
            Some(5),
        );
        let nasty_receipt =
            fs_rep_mesh::mesh_to_sdf(&MeshChart::new(nasty), 0.15, cx).expect("soup builds");
        let downgrade_ok = nasty_receipt.numerical.kind == fs_evidence::NumericalKind::Estimate
            && nasty_receipt.value.abstract_distance_kind() == fs_evidence::NumericalKind::Estimate
            && nasty_receipt
                .model
                .cards
                .contains(&"winding-sign-heuristic".to_string());
        (
            analytic_ok,
            equivariant,
            incremental_identical,
            downgrade_ok,
            inc.last_update_samples,
        )
    });
    let mut em = fs_obs::Emitter::new("fs-rep-mesh/conformance", "rmesh-007/convert");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "rep-mesh-convert-stats".to_string(),
                json: format!("{{\"incremental_samples_refreshed\":{samples}}}"),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("convert stats validate");
    println!("{line}");
    verdict(
        "rmesh-007",
        analytic_ok && equivariant && incremental_identical && downgrade_ok,
        &format!(
            "mesh->SDF: analytic match within recorded estimate, translation-equivariant (G3), \
             incremental update bit-identical to full rebuild (G5, {samples} samples \
             refreshed), generic-soup payloads and receipts honestly capped at Estimate \
             (seed {SEED:#x})"
        ),
    );
}

#[test]
// One end-to-end scenario over shared reconstructions.
#[allow(clippy::too_many_lines)]
fn rmesh_008_dual_contouring_reconstructs_certifies_and_detects_bad_triangles() {
    use fs_geom::fixtures::{BoxChart, SphereChart};
    use fs_rep_mesh::{DcOptions, bracket_certificate, dual_contour};
    let (sphere_err, cert_ok, cert_margin, manifold_closed, winding_one, equivariant) =
        with_cx(|cx| {
            let sphere = SphereChart {
                center: Point3::new(0.05, -0.1, 0.02),
                radius: 1.0,
            };
            let (soup, stats) = dual_contour(&sphere, DcOptions::sharp(0.11), cx).expect("dc");
            assert!(stats.triangles > 100, "{}", stats.to_json());
            // Reconstruction accuracy: every DC vertex near the zero set.
            let mut worst = 0.0f64;
            for &v in &soup.positions {
                worst = worst.max((v.delta_from(sphere.center).norm() - 1.0).abs());
            }
            // Bracket certificate: proven within tolerance everywhere.
            let cert = bracket_certificate(&sphere, &soup, 0.2, cx).expect("exact-distance chart");
            let (cert_ok, margin) = match cert {
                Ok(report) => (true, report.worst_margin),
                Err(fails) => {
                    println!("bracket failures: {fails:?}");
                    (false, f64::NAN)
                }
            };
            // Manifold + closed + correctly oriented (winding at center = 1).
            let manifold =
                fs_rep_mesh::HalfEdgeMesh::from_triangles(soup.positions.clone(), &soup.triangles)
                    .is_ok();
            let closed = fs_rep_mesh::assess_quality(&soup).passes_basic_orientation_checks();
            let w = fs_rep_mesh::winding_exact(&soup, sphere.center);
            // G3 rigid motion: translate the chart; vertices translate.
            let shift = fs_geom::Vec3::new(0.5, 0.25, -0.375); // dyadic: exact fp
            let moved = SphereChart {
                center: sphere.center.offset(shift),
                radius: 1.0,
            };
            let (soup2, _) = dual_contour(&moved, DcOptions::sharp(0.11), cx).expect("dc moved");
            let mut equi = soup2.positions.len() == soup.positions.len();
            if equi {
                for (a, b) in soup.positions.iter().zip(&soup2.positions) {
                    let back = b.offset(shift.scale(-1.0));
                    equi &= back.delta_from(*a).norm() < 1e-9;
                }
            }
            (
                worst,
                cert_ok,
                margin,
                manifold && closed,
                (w - 1.0).abs() < 1e-9,
                equi,
            )
        });
    // Sharp features: box corners resolved by QEF, blurred by mass point.
    let (qef_corner_err, mass_corner_err, detects_bad) = with_cx(|cx| {
        use fs_rep_mesh::{DcOptions, Placement, bracket_certificate, dual_contour};
        let bx = BoxChart {
            aabb: fs_geom::Aabb::new(Point3::new(-0.8, -0.6, -0.7), Point3::new(0.7, 0.9, 0.6)),
        };
        let corner = Point3::new(0.7, 0.9, 0.6);
        let corner_dist = |soup: &fs_rep_mesh::Soup| {
            soup.positions
                .iter()
                .map(|v| v.delta_from(corner).norm())
                .fold(f64::INFINITY, f64::min)
        };
        let (qef, _) = dual_contour(&bx, DcOptions::sharp(0.1), cx).expect("qef");
        let (mass, _) = dual_contour(
            &bx,
            DcOptions {
                placement: Placement::MassPoint,
                ..DcOptions::sharp(0.1)
            },
            cx,
        )
        .expect("mass");
        // Seeded bad triangle: yank one vertex far off the surface; the
        // certificate must fail and LOCALIZE.
        let mut broken = qef.clone();
        broken.positions[7] = Point3::new(3.0, 3.0, 3.0);
        let verdict = bracket_certificate(&bx, &broken, 0.2, cx).expect("exact-distance chart");
        let detects = matches!(&verdict, Err(fails) if !fails.is_empty()
            && fails.iter().all(|f| f.proven_bound > f.tolerance));
        (corner_dist(&qef), corner_dist(&mass), detects)
    });
    let mut em = fs_obs::Emitter::new("fs-rep-mesh/conformance", "rmesh-008/dc");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "rep-mesh-dc-stats".to_string(),
                json: format!(
                    "{{\"sphere_worst_vertex_err\":{sphere_err:.5},\"cert_margin\":{cert_margin:.5},\
                     \"qef_corner_err\":{qef_corner_err:.5},\"mass_corner_err\":{mass_corner_err:.5}}}"
                ),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("dc stats validate");
    println!("{line}");
    verdict(
        "rmesh-008",
        sphere_err < 0.06
            && cert_ok
            && manifold_closed
            && winding_one
            && equivariant
            && qef_corner_err < 0.5 * mass_corner_err
            && detects_bad,
        &format!(
            "DC sphere vertices within {sphere_err:.4} of the zero set with the bracket \
             certificate proven (margin {cert_margin:.4}); output manifold, closed, and \
             outward-oriented; translation-equivariant (G3); QEF resolves the box corner \
             {qef_corner_err:.4} vs mass-point {mass_corner_err:.4}; a seeded bad triangle is \
             caught and localized"
        ),
    );
}

/// rmesh-009 — dual contouring has an explicit finite-domain boundary:
/// unresolved extended supports and invalid/excessive grids refuse before
/// chart evaluation, while the clipped API contours `chart ∩ clip` and is
/// translation-equivariant (G3).
#[test]
fn rmesh_009_dual_contour_sampling_domain_is_explicit_and_preflighted() {
    with_cx(|cx| {
        let clip = Aabb::new(Point3::new(-0.5, -0.5, -0.5), Point3::new(0.5, 0.5, 0.5));

        let unbounded = CountingPlane::new(0.0);
        assert!(matches!(
            dual_contour(&unbounded, DcOptions::sharp(0.125), cx),
            Err(ContourError::SamplingDomain(
                SamplingDomainError::UnboundedSupport { .. }
            ))
        ));
        assert_eq!(unbounded.eval_count(), 0, "support refusal precedes eval");

        let invalid = CountingPlane::new(0.0);
        assert!(matches!(
            dual_contour_clipped(&invalid, DcOptions::sharp(0.0), clip, cx),
            Err(ContourError::InvalidSpacing { .. })
        ));
        assert_eq!(invalid.eval_count(), 0, "invalid h precedes eval");

        let too_fine = CountingPlane::new(0.0);
        assert!(matches!(
            dual_contour_clipped(&too_fine, DcOptions::sharp(0.001), clip, cx),
            Err(ContourError::ResolutionTooFine {
                cap: DC_MAX_CELLS_PER_AXIS,
                ..
            })
        ));
        assert_eq!(too_fine.eval_count(), 0, "grid cap precedes eval");

        let source = CountingPlane::new(0.0);
        let (soup, stats) = dual_contour_clipped(&source, DcOptions::sharp(0.125), clip, cx)
            .expect("clipped half-box contours");
        assert!(stats.triangles > 0 && !soup.positions.is_empty());
        assert!(
            soup.positions
                .iter()
                .all(|p| clip.inflate(0.126).contains(*p)),
            "contour vertices stay in the clipped geometry's sampled halo"
        );

        let shift = Vec3::new(0.25, -0.125, 0.375);
        let moved_clip = Aabb::new(clip.min.offset(shift), clip.max.offset(shift));
        let (moved, _) = dual_contour_clipped(
            &CountingPlane::new(shift.x),
            DcOptions::sharp(0.125),
            moved_clip,
            cx,
        )
        .expect("translated clipped half-box contours");
        assert_eq!(soup.positions.len(), moved.positions.len());
        assert_eq!(soup.triangles, moved.triangles);
        for (a, b) in soup.positions.iter().zip(&moved.positions) {
            assert!(
                b.offset(shift.scale(-1.0)).delta_from(*a).norm() < 1e-9,
                "G3 contour vertex mismatch"
            );
        }
    });
    verdict(
        "rmesh-009",
        true,
        "dual contouring rejects unresolved extended support and invalid/excessive grids before evaluation; the clipped API contours source-intersection-clip and preserves translation equivariance (G3)",
    );
}

/// rmesh-010 — bracket authority is global and evidence-backed: a local
/// `Some(1)` cannot promote a no-claim/estimate chart, and the bracket
/// consumer observes cancellation requested by a non-polling chart.
#[test]
fn rmesh_010_bracket_authority_fails_closed_and_polls_directly() {
    let soup = fs_rep_mesh::Soup {
        positions: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.25, 0.0),
            Point3::new(0.0, 0.0, 0.25),
        ],
        triangles: vec![[0, 1, 2]],
    };

    with_cx(|cx| {
        let no_claim = BracketAuthorityPlane {
            mode: BracketAuthorityMode::NoClaim,
            evals: AtomicU64::new(0),
        };
        let error = bracket_certificate(&no_claim, &soup, 0.5, cx)
            .expect_err("local Some(1) must not upgrade NoClaim");
        assert_eq!(
            error,
            BracketCertificateError::UnsupportedTraceClaim {
                actual: TraceStepClaim::NoClaim,
            }
        );
        assert_eq!(
            no_claim.eval_count(),
            0,
            "trace-theorem refusal precedes evaluation"
        );

        let estimate = BracketAuthorityPlane {
            mode: BracketAuthorityMode::Estimate,
            evals: AtomicU64::new(0),
        };
        let error = bracket_certificate(&estimate, &soup, 0.5, cx)
            .expect_err("Estimate evidence must not become a verdict");
        assert!(matches!(
            error,
            BracketCertificateError::InvalidTraceEvidence {
                completed_evaluations: 1,
                issue: BracketEvidenceIssue::NonRigorous {
                    kind: NumericalKind::Estimate,
                },
                ..
            }
        ));
        assert_eq!(estimate.eval_count(), 1);

        let no_evidence = BracketAuthorityPlane {
            mode: BracketAuthorityMode::NoClaimEvidence,
            evals: AtomicU64::new(0),
        };
        let error = bracket_certificate(&no_evidence, &soup, 0.5, cx)
            .expect_err("NoClaim evidence must not become a verdict");
        assert!(matches!(
            error,
            BracketCertificateError::InvalidTraceEvidence {
                issue: BracketEvidenceIssue::NonRigorous {
                    kind: NumericalKind::NoClaim,
                },
                ..
            }
        ));

        let malformed = BracketAuthorityPlane {
            mode: BracketAuthorityMode::MalformedExact,
            evals: AtomicU64::new(0),
        };
        let error = bracket_certificate(&malformed, &soup, 0.5, cx)
            .expect_err("malformed Exact evidence must not become a verdict");
        assert!(matches!(
            error,
            BracketCertificateError::InvalidTraceEvidence {
                issue: BracketEvidenceIssue::MalformedExact,
                ..
            }
        ));
    });

    let gate = CancelGate::new();
    let cancelling = BracketAuthorityPlane {
        mode: BracketAuthorityMode::RequestCancellation(&gate),
        evals: AtomicU64::new(0),
    };
    with_gate_cx(&gate, |cx| {
        let error = bracket_certificate(&cancelling, &soup, 0.5, cx)
            .expect_err("bracket consumer must observe requested cancellation");
        assert_eq!(
            error,
            BracketCertificateError::Cancelled {
                completed_triangles: 0,
                completed_evaluations: 1,
            }
        );
    });
    assert_eq!(cancelling.eval_count(), 1);

    verdict(
        "rmesh-010",
        true,
        "NoClaim, Estimate, and malformed evidence cannot borrow authority from a local Lipschitz sample; direct post-evaluation polling reports cancellation progress",
    );
}
