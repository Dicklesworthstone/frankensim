//! fs-topo conformance suite (CONTRACT.md: any reimplementation must
//! pass). Manifold certificates on the defect zoo, exact
//! self-intersection proofs on the adversarial zoo, cubical Betti
//! numbers on the fixture solids, true 0-dim persistence with planted
//! features under noise, the stability theorem as a property test, and
//! chart-level topology verification with determinism and a ledgered
//! scale run. Verdicts use the canonical fs-obs schema; seeded cases carry
//! replay seeds.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Aabb, Chart, ChartSample, Differentiability, Point3, SamplingDomainError, Vec3};
use fs_rep_frep::{BoolOp, BoolStyle, Frep, FrepBuilder};
use fs_rep_mesh::{Soup, shapes};
use fs_topo::cubical::{
    MAX_VOXELIZE_CELLS, VoxelField, VoxelizeError, betti, count_persistent, persistence0,
    verify_topology, verify_topology_clipped, voxelize, voxelize_clipped,
};
use fs_topo::{IntersectKind, ManifoldDefect, manifold_certificate, self_intersection_certificate};
use std::sync::atomic::{AtomicU64, Ordering};

fn verdict(case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new("fs-topo/conformance", case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: "fs-topo/conformance".to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("topology verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("topology verdict must use the fs-obs wire schema");
    println!("{line}");
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
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x707,
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
    chart: Frep,
    evals: AtomicU64,
}

impl CountingPlane {
    fn new(x_offset: f64) -> Self {
        let mut builder = FrepBuilder::new();
        let plane = builder
            .half_space(Vec3::new(1.0, 0.0, 0.0), x_offset)
            .expect("valid plane");
        Self {
            chart: builder.finish(plane).expect("valid plane F-rep"),
            evals: AtomicU64::new(0),
        }
    }

    fn eval_count(&self) -> u64 {
        self.evals.load(Ordering::Relaxed)
    }
}

impl Chart for CountingPlane {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        self.evals.fetch_add(1, Ordering::Relaxed);
        self.chart.eval(x, cx)
    }

    fn support(&self) -> Aabb {
        self.chart.support()
    }

    fn name(&self) -> &'static str {
        "test/counting-plane"
    }

    fn differentiability(&self) -> Differentiability {
        self.chart.differentiability()
    }
}

struct BoundsCheckingChart {
    bounds: Aabb,
    evals: AtomicU64,
}

impl Chart for BoundsCheckingChart {
    fn eval(&self, point: Point3, _cx: &Cx<'_>) -> ChartSample {
        assert!(
            point.x.is_finite() && point.y.is_finite() && point.z.is_finite(),
            "voxel center must be finite: {point:?}"
        );
        assert!(
            self.bounds.contains(point),
            "voxel center escaped admitted bounds: {point:?} not in {:?}",
            self.bounds
        );
        self.evals.fetch_add(1, Ordering::Relaxed);
        ChartSample {
            signed_distance: point.x,
            gradient: None,
            lipschitz: Some(1.0),
            error: fs_evidence::NumericalCertificate::exact(point.x),
        }
    }

    fn support(&self) -> Aabb {
        self.bounds
    }

    fn name(&self) -> &'static str {
        "test/bounds-checking-chart"
    }
}

/// topo-001 — manifold certificates: clean fixtures certify; every
/// defect class is DETECTED and LOCALIZED on the zoo.
#[test]
fn topo_001_manifold_zoo() {
    let clean = shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 2);
    let cert = manifold_certificate(&clean, Some(Point3::new(0.0, 0.0, 0.0)));
    let clean_ok = cert.certified();

    // Punched hole → boundary edges localized.
    let holed = shapes::corrupt(clean.clone(), 0, 0, 0..0, Some(7));
    let hc = manifold_certificate(&holed, None);
    let hole_ok = !hc.closed
        && hc
            .defects
            .iter()
            .filter(|d| matches!(d, ManifoldDefect::BoundaryEdge { .. }))
            .count()
            == 3;

    // Duplicate face → non-manifold edges with use-count 3.
    let dup = shapes::corrupt(clean.clone(), 1, 0, 0..0, None);
    let dc = manifold_certificate(&dup, None);
    let dup_ok = !dc.manifold
        && dc
            .defects
            .iter()
            .any(|d| matches!(d, ManifoldDefect::NonManifoldEdge { uses: 3, .. }));

    // Flipped patch → misoriented edges.
    let flipped = shapes::corrupt(clean.clone(), 0, 0, 4..5, None);
    let fc = manifold_certificate(&flipped, None);
    let flip_ok = !fc.oriented
        && fc
            .defects
            .iter()
            .any(|d| matches!(d, ManifoldDefect::MisorientedEdge { .. }));

    // Degenerate face.
    let degen = shapes::corrupt(clean, 0, 2, 0..0, None);
    let gc = manifold_certificate(&degen, None);
    let degen_ok = gc
        .defects
        .iter()
        .any(|d| matches!(d, ManifoldDefect::DegenerateFace { .. }));

    verdict(
        "topo-001",
        clean_ok && hole_ok && dup_ok && flip_ok && degen_ok,
        "the clean icosphere certifies (manifold, closed, outward); a punched hole \
         localizes exactly 3 boundary edges; a duplicated face reads use-count 3; a \
         flipped patch reads misoriented edges; degenerate faces are named by index",
        0,
    );
}

/// topo-002 — self-intersection certificates: clean surfaces PROVEN
/// free; planted piercings, coincident patches detected and localized;
/// near-tangent surfaces do NOT false-FAIL; adjacency never flags.
#[test]
fn topo_002_self_intersection() {
    // Clean: proven free.
    let clean = shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 3);
    let rep = self_intersection_certificate(&clean);
    let clean_ok = rep.proven_free() && rep.pairs_tested > 0;

    // Spike a vertex through the far wall.
    let mut pierced = shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 2);
    let p0 = pierced.positions[0];
    pierced.positions[0] = Point3::new(-1.4 * p0.x, -1.4 * p0.y, -1.4 * p0.z);
    let rp = self_intersection_certificate(&pierced);
    let pierced_ok = !rp.proven_free()
        && rp
            .intersections
            .iter()
            .any(|&(_, _, k)| k == IntersectKind::Crossing);

    // Coincident patch: duplicate a face with FRESH vertices (no
    // shared indices → not adjacency-excused).
    let mut coincident = shapes::icosphere(Point3::new(0.0, 0.0, 0.0), 1.0, 2);
    let [a, b, c] = coincident.triangles[0];
    let base = coincident.positions.len() as u32;
    coincident.positions.push(coincident.positions[a as usize]);
    coincident.positions.push(coincident.positions[b as usize]);
    coincident.positions.push(coincident.positions[c as usize]);
    coincident.triangles.push([base, base + 1, base + 2]);
    let rc = self_intersection_certificate(&coincident);
    let coincident_ok = rc
        .intersections
        .iter()
        .any(|&(_, _, k)| k == IntersectKind::Touching);

    // Near-tangent (1e-4 gap): must NOT false-FAIL.
    let s1 = shapes::icosphere(Point3::new(-1.0 - 5e-5, 0.0, 0.0), 1.0, 2);
    let s2 = shapes::icosphere(Point3::new(1.0 + 5e-5, 0.0, 0.0), 1.0, 2);
    let mut merged = Soup {
        positions: s1.positions.clone(),
        triangles: s1.triangles.clone(),
    };
    let off = merged.positions.len() as u32;
    merged.positions.extend(s2.positions.iter().copied());
    merged
        .triangles
        .extend(s2.triangles.iter().map(|t| t.map(|v| v + off)));
    let rt = self_intersection_certificate(&merged);
    let tangent_ok = rt.proven_free();

    verdict(
        "topo-002",
        clean_ok && pierced_ok && coincident_ok && tangent_ok,
        &format!(
            "[clean={clean_ok} pierced={pierced_ok} coincident={coincident_ok} tangent={tangent_ok}] \
             the clean icosphere is PROVEN free ({} candidate pairs, exact narrow \
             phase), a planted spike reads Crossing with localization, an exactly \
             coincident patch reads Touching (the bounded conservative class), and \
             near-tangent spheres at 1e-4 separation do NOT false-FAIL",
            rep.pairs_tested
        ),
        0,
    );
}

fn frep_ball(c: Point3, r: f64) -> Frep {
    let mut b = FrepBuilder::new();
    let s = b.sphere(c, r).expect("s");
    b.finish(s).expect("frep")
}

/// topo-003 — cubical Betti numbers match the fixture zoo exactly:
/// ball (1,0,0), solid torus (1,1,0), hollow ball (1,0,1), two balls
/// (2,0,0).
#[test]
fn topo_003_cubical_betti() {
    with_cx(|cx| {
        let ball =
            verify_topology(&frep_ball(Point3::new(0.0, 0.0, 0.0), 1.0), 40, cx).expect("ball");
        let torus = {
            let mut b = FrepBuilder::new();
            let t = b.torus(Point3::new(0.0, 0.0, 0.0), 1.0, 0.35).expect("t");
            let f = b.finish(t).expect("frep");
            verify_topology(&f, 48, cx).expect("torus")
        };
        let hollow = {
            let mut b = FrepBuilder::new();
            let outer = b.sphere(Point3::new(0.0, 0.0, 0.0), 1.0).expect("o");
            let inner = b.sphere(Point3::new(0.0, 0.0, 0.0), 0.5).expect("i");
            let sh = b
                .boolean(BoolOp::Difference, BoolStyle::Hard, outer, inner)
                .expect("d");
            let f = b.finish(sh).expect("frep");
            verify_topology(&f, 48, cx).expect("hollow")
        };
        let two = {
            let mut b = FrepBuilder::new();
            let s1 = b.sphere(Point3::new(-1.5, 0.0, 0.0), 0.7).expect("s1");
            let s2 = b.sphere(Point3::new(1.5, 0.0, 0.0), 0.7).expect("s2");
            let u = b
                .boolean(BoolOp::Union, BoolStyle::Hard, s1, s2)
                .expect("u");
            let f = b.finish(u).expect("frep");
            verify_topology(&f, 48, cx).expect("two")
        };
        verdict(
            "topo-003",
            ball == (1, 0, 0) && torus == (1, 1, 0) && hollow == (1, 0, 1) && two == (2, 0, 0),
            &format!(
                "Betti triples read exactly: ball {ball:?} = (1,0,0), solid torus \
                 {torus:?} = (1,1,0) [the tunnel via Euler duality], hollow ball \
                 {hollow:?} = (1,0,1) [the cavity], two balls {two:?} = (2,0,0)"
            ),
            0,
        );
    });
}

#[test]
fn topo_007_diagonal_contact_is_not_a_phantom_tunnel() {
    // Two closed unit cubes touching only along an EDGE (2×2×1 grid, the
    // anti-diagonal filled). A union of CLOSED cubes is 26-connected, so this
    // is ONE contractible piece: true Betti (1,0,0). A 6-connected b0 made
    // b1 = b0 + b2 − χ manufacture a phantom tunnel (2,1,0) (bead rcvl).
    let edge = VoxelField {
        dims: [2, 2, 1],
        values: vec![-1.0, 1.0, 1.0, -1.0],
        h: 1.0,
    };
    let be = betti(&edge, 0.0);
    // Corner contact (2×2×2, opposite corners filled) — also one piece.
    let mut vals = vec![1.0; 8];
    vals[0] = -1.0; // (0,0,0)
    vals[7] = -1.0; // (1,1,1)
    let corner = VoxelField {
        dims: [2, 2, 2],
        values: vals,
        h: 1.0,
    };
    let bc = betti(&corner, 0.0);
    verdict(
        "topo-007",
        be == (1, 0, 0) && bc == (1, 0, 0),
        &format!(
            "diagonal voxel contacts are ONE 26-connected contractible component — \
             edge {be:?} and corner {bc:?} both = (1,0,0), no phantom tunnel"
        ),
        0,
    );
}

/// Two-well synthetic field with deterministic noise.
fn two_well_field(noise: f64, seed: u64) -> VoxelField {
    let n = 36u32;
    let c1 = [9.0, 18.0, 18.0];
    let c2 = [27.0, 18.0, 18.0];
    let mut rng = Lcg(seed);
    let mut values = Vec::with_capacity((n * n * n) as usize);
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let d = |c: [f64; 3]| -> f64 {
                    let dx = f64::from(x) - c[0];
                    let dy = f64::from(y) - c[1];
                    let dz = f64::from(z) - c[2];
                    (dx * dx + dy * dy + dz * dz).sqrt() / 18.0
                };
                let base = d(c1).min(d(c2) + 0.2);
                values.push(base + noise * (2.0 * rng.unit() - 1.0));
            }
        }
    }
    VoxelField {
        dims: [n, n, n],
        values,
        h: 1.0,
    }
}

/// topo-004 — 0-dim persistence: two planted wells stay exactly two
/// persistent features under noise (elder rule: the deeper well is the
/// essential class; the shallower dies at the saddle), while noise
/// bars stay short.
#[test]
fn topo_004_persistence_planted() {
    let field = two_well_field(0.05, 0x1001_2026_0707_0034);
    let bars = persistence0(&field);
    let persistent = count_persistent(&bars, 0.15);
    let essential: Vec<_> = bars.iter().filter(|b| b.death.is_infinite()).collect();
    // The saddle between wells sits near d/2 (wells 18 apart, scale 9
    // → saddle ≈ 1.0 + noise); the shallow well is born ≈ 0.2.
    let shallow = bars
        .iter()
        .filter(|b| b.death.is_finite() && b.persistence() > 0.15)
        .max_by(|a, b| {
            a.persistence()
                .partial_cmp(&b.persistence())
                .expect("finite")
        });
    let shallow_ok = shallow.is_some_and(|b| (b.birth - 0.2).abs() < 0.08);
    let noise_bars = bars
        .iter()
        .filter(|b| b.death.is_finite() && b.persistence() <= 0.1)
        .count();
    verdict(
        "topo-004",
        persistent == 2 && essential.len() == 1 && shallow_ok && noise_bars > 50,
        &format!(
            "exactly 2 persistent features at tau=0.15 (1 essential + the planted \
             shallow well born at {:.3} ~ 0.2), against {noise_bars} short noise \
             bars; seed 0x1001_2026_0707_0034",
            shallow.map_or(f64::NAN, |b| b.birth)
        ),
        0x1001_2026_0707_0034,
    );
}

/// topo-005 — the stability theorem as a property test: perturbing the
/// field by ≤ ε moves every surviving bar's endpoints by ≤ ε.
#[test]
fn topo_005_stability() {
    let eps = 0.02;
    let base = two_well_field(0.015, 0x1001_2026_0707_0035);
    let mut rng = Lcg(0x1001_2026_0707_0036);
    let mut perturbed = base.clone();
    for v in &mut perturbed.values {
        *v += eps * (2.0 * rng.unit() - 1.0);
    }
    let bars_a = persistence0(&base);
    let bars_b = persistence0(&perturbed);
    let survive = |bars: &[fs_topo::cubical::Bar]| -> Vec<fs_topo::cubical::Bar> {
        let mut v: Vec<_> = bars
            .iter()
            .copied()
            .filter(|b| b.persistence() > 2.0 * eps + 1e-9)
            .collect();
        v.sort_by(|a, b| a.birth.partial_cmp(&b.birth).expect("finite"));
        v
    };
    let (sa, sb) = (survive(&bars_a), survive(&bars_b));
    let counts_match = sa.len() == sb.len();
    let mut moved_ok = counts_match;
    if counts_match {
        for (a, b) in sa.iter().zip(&sb) {
            moved_ok &= (a.birth - b.birth).abs() <= eps + 1e-9;
            if a.death.is_finite() && b.death.is_finite() {
                moved_ok &= (a.death - b.death).abs() <= eps + 1e-9;
            } else {
                moved_ok &= a.death.is_infinite() == b.death.is_infinite();
            }
        }
    }
    verdict(
        "topo-005",
        counts_match && moved_ok,
        &format!(
            "under an eps={eps} perturbation, the {} bars with persistence > 2eps \
             survive with birth/death endpoints moved <= eps (the bottleneck \
             stability theorem, spot-verified); seed 0x1001_2026_0707_0035/36",
            sa.len()
        ),
        0x1001_2026_0707_0035,
    );
}

/// topo-006 — determinism + the ledgered scale run: persistence is
/// bitwise reproducible; a 96³ field's Betti + persistence timings go
/// to the ledger.
#[test]
fn topo_006_determinism_and_scale() {
    let field = two_well_field(0.02, 0x1001_2026_0707_0037);
    let a = persistence0(&field);
    let b = persistence0(&field);
    let bitwise = a.len() == b.len()
        && a.iter().zip(&b).all(|(x, y)| {
            x.birth.to_bits() == y.birth.to_bits() && x.death.to_bits() == y.death.to_bits()
        });
    // Scale run: 96³ ≈ 885k voxels.
    let n = 96u32;
    let mut rng = Lcg(0x1001_2026_0707_0038);
    let mut values = Vec::with_capacity((n * n * n) as usize);
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let d = f64::from(x.min(y).min(z)) / f64::from(n);
                values.push(d + 0.05 * rng.unit());
            }
        }
    }
    let big = VoxelField {
        dims: [n, n, n],
        values,
        h: 1.0,
    };
    let t0 = std::time::Instant::now();
    let bt = betti(&big, 0.3);
    let t_betti = t0.elapsed().as_millis();
    let t1 = std::time::Instant::now();
    let bars = persistence0(&big);
    let t_pers = t1.elapsed().as_millis();
    let mut em = fs_obs::Emitter::new("fs-topo/conformance", "topo-006/scale");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "topo-scale-ledger".to_string(),
                json: format!(
                    "{{\"voxels\":{},\"betti\":[{},{},{}],\"betti_ms\":{t_betti},\
                     \"bars\":{},\"persistence_ms\":{t_pers}}}",
                    n * n * n,
                    bt.0,
                    bt.1,
                    bt.2,
                    bars.len()
                ),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("scale ledger validates");
    println!("{line}");
    verdict(
        "topo-006",
        bitwise && bt.0 >= 1,
        &format!(
            "persistence is BITWISE reproducible; the 885k-voxel scale run ledgers \
             betti {bt:?} in {t_betti} ms and {} bars in {t_pers} ms (sequential v1; \
             chunked-parallel is the contract no-claim); seed 0x1001_2026_0707_0038",
            bars.len()
        ),
        0x1001_2026_0707_0038,
    );
}

/// topo-008 — chart voxelization/topology verification require an explicit
/// finite sampling domain. Invalid resolutions and excessive products refuse
/// before evaluation; clipped fields sample the geometric intersection and
/// preserve occupancy and Betti numbers under translation (G3).
#[test]
fn topo_008_chart_sampling_domain_is_explicit_and_preflighted() {
    with_cx(|cx| {
        let clip = Aabb::new(Point3::new(-0.5, -0.5, -0.5), Point3::new(0.5, 0.5, 0.5));

        let unbounded = CountingPlane::new(0.0);
        assert!(matches!(
            voxelize(&unbounded, 24, cx),
            Err(VoxelizeError::SamplingDomain(
                SamplingDomainError::UnboundedSupport { .. }
            ))
        ));
        assert_eq!(unbounded.eval_count(), 0, "support refusal precedes eval");

        let verify_unbounded = CountingPlane::new(0.0);
        assert!(matches!(
            verify_topology(&verify_unbounded, 24, cx),
            Err(VoxelizeError::SamplingDomain(
                SamplingDomainError::UnboundedSupport { .. }
            ))
        ));
        assert_eq!(verify_unbounded.eval_count(), 0);

        let invalid = CountingPlane::new(0.0);
        assert!(matches!(
            voxelize_clipped(&invalid, 0, clip, cx),
            Err(VoxelizeError::InvalidResolution { n: 0 })
        ));
        assert_eq!(invalid.eval_count(), 0, "zero resolution precedes eval");

        let too_large = CountingPlane::new(0.0);
        assert!(matches!(
            voxelize_clipped(&too_large, 101, clip, cx),
            Err(VoxelizeError::VoxelLimit {
                cap: MAX_VOXELIZE_CELLS,
                ..
            })
        ));
        assert_eq!(too_large.eval_count(), 0, "voxel cap precedes eval");

        let field = voxelize_clipped(&CountingPlane::new(0.0), 24, clip, cx)
            .expect("clipped half-box voxelizes");
        let topology = verify_topology_clipped(&CountingPlane::new(0.0), 24, clip, cx)
            .expect("clipped half-box topology");
        assert_eq!(topology, (1, 0, 0), "half-box is contractible");

        let shift = Vec3::new(0.25, -0.125, 0.375);
        let moved_clip = Aabb::new(clip.min.offset(shift), clip.max.offset(shift));
        let moved = voxelize_clipped(&CountingPlane::new(shift.x), 24, moved_clip, cx)
            .expect("translated clipped half-box voxelizes");
        assert_eq!(field.dims, moved.dims);
        assert_eq!(
            field
                .values
                .iter()
                .map(|value| *value < 0.0)
                .collect::<Vec<_>>(),
            moved
                .values
                .iter()
                .map(|value| *value < 0.0)
                .collect::<Vec<_>>(),
            "G3 occupancy mismatch"
        );
        assert_eq!(betti(&field, 0.0), betti(&moved, 0.0));

        let extreme_bounds = Aabb::new(
            Point3::new(-8.5e307, 1.79e308, -1.0),
            Point3::new(8.5e307, f64::MAX, 1.0),
        );
        let extreme = BoundsCheckingChart {
            bounds: extreme_bounds,
            evals: AtomicU64::new(0),
        };
        let extreme_field = voxelize(&extreme, 10, cx)
            .expect("ratio-first centers stay inside an extreme-aspect finite AABB");
        assert_eq!(extreme_field.dims, [10, 1, 1]);
        assert_eq!(extreme.evals.load(Ordering::Relaxed), 10);

        let threshold_bounds = Aabb::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.900_000_000_000_000_1, 0.1),
        );
        let threshold = BoundsCheckingChart {
            bounds: threshold_bounds,
            evals: AtomicU64::new(0),
        };
        let threshold_field = voxelize(&threshold, 10, cx)
            .expect("near-integer quotient preserves the maximum-width promise");
        assert_eq!(threshold_field.dims, [10, 10, 1]);
        for (span, dim) in [1.0, 0.900_000_000_000_000_1, 0.1]
            .into_iter()
            .zip(threshold_field.dims)
        {
            assert!(
                span / f64::from(dim) <= threshold_field.h,
                "realized cell width must not exceed reported h"
            );
        }
    });
    verdict(
        "topo-008",
        true,
        "voxelize and verify_topology reject unresolved extended support and invalid/excessive resolution before evaluation; paired clipped APIs sample source-intersection-clip and preserve occupancy and Betti numbers under translation (G3); ratio-first centers remain finite and inside extreme-aspect admitted domains, and realized widths stay within reported h across near-integer quotients",
        0,
    );
}

/// Vec3 kept in scope for fixture builders.
#[allow(dead_code)]
fn _unused(v: Vec3) -> f64 {
    v.norm()
}
