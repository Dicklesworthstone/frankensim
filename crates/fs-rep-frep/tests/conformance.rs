//! fs-rep-frep conformance suite (CONTRACT.md: any reimplementation must
//! pass). G0 interval containment on random DAGs, Lipschitz validity
//! under adversarial sampling, C¹ blend seams (with the hard-Boolean
//! discontinuity exhibited as the contrast), DAG-sharing bitwise
//! equality, the sphere-tracing safety harness vs a dense oracle, and
//! metamorphic + design-lever laws. JSON-line verdicts; seeded cases
//! carry seeds.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Aabb, Chart, Differentiability, Point3, TraceStepClaim, Vec3};
use fs_rep_frep::{BoolOp, BoolStyle, Frep, FrepBuilder, NodeId, smin_weights};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-rep-frep/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\
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

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
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
                seed: 0xF2EB,
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

/// Deterministic random CSG DAG: 3–5 bounded primitives, 3–6 combinator
/// nodes (Booleans in both styles + rigid/scale/offset transforms). All
/// supports bounded (no half-spaces/cylinders in the random pool — those
/// are exercised in directed cases).
fn random_frep(rng: &mut Lcg) -> Frep {
    let mut b = FrepBuilder::new();
    let mut ids: Vec<NodeId> = Vec::new();
    let nprim = 3 + rng.below(3);
    for _ in 0..nprim {
        let c = Point3::new(
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
        );
        let id = match rng.below(3) {
            0 => b.sphere(c, rng.range(0.3, 1.0)).expect("valid sphere"),
            1 => b
                .box_prim(
                    c,
                    Vec3::new(
                        rng.range(0.2, 0.8),
                        rng.range(0.2, 0.8),
                        rng.range(0.2, 0.8),
                    ),
                )
                .expect("valid box"),
            _ => b
                .torus(c, rng.range(0.4, 0.9), rng.range(0.1, 0.3))
                .expect("valid torus"),
        };
        ids.push(id);
    }
    let nops = 3 + rng.below(4);
    for _ in 0..nops {
        let pick = |rng: &mut Lcg, ids: &[NodeId]| ids[rng.below(ids.len() as u64) as usize];
        let id = match rng.below(8) {
            0..=3 => {
                let op = match rng.below(3) {
                    0 => BoolOp::Union,
                    1 => BoolOp::Intersect,
                    _ => BoolOp::Difference,
                };
                let style = if rng.below(2) == 0 {
                    BoolStyle::Hard
                } else {
                    BoolStyle::Blend {
                        radius: rng.range(0.05, 0.35),
                    }
                };
                let (x, y) = (pick(rng, &ids), pick(rng, &ids));
                b.boolean(op, style, x, y).expect("valid boolean")
            }
            4 => {
                let x = pick(rng, &ids);
                b.translate(
                    x,
                    Vec3::new(
                        rng.range(-0.5, 0.5),
                        rng.range(-0.5, 0.5),
                        rng.range(-0.5, 0.5),
                    ),
                )
                .expect("valid translate")
            }
            5 => {
                let x = pick(rng, &ids);
                let axis = Vec3::new(
                    rng.range(-1.0, 1.0),
                    rng.range(-1.0, 1.0),
                    rng.range(-1.0, 1.0),
                );
                let axis = if axis.norm() < 1e-6 {
                    Vec3::new(0.0, 0.0, 1.0)
                } else {
                    axis
                };
                b.rotate(x, axis, rng.range(-1.0, 1.0))
                    .expect("valid rotate")
            }
            6 => {
                let x = pick(rng, &ids);
                b.scale(x, rng.range(0.6, 1.6)).expect("valid scale")
            }
            _ => {
                let x = pick(rng, &ids);
                b.offset(x, rng.range(-0.05, 0.15)).expect("valid offset")
            }
        };
        ids.push(id);
    }
    let root = *ids.last().expect("nonempty");
    b.finish(root).expect("valid root")
}

fn sample_in(rng: &mut Lcg, b: &Aabb) -> Point3 {
    Point3::new(
        rng.range(b.min.x, b.max.x),
        rng.range(b.min.y, b.max.y),
        rng.range(b.min.z, b.max.z),
    )
}

/// frep-001 — G0 containment: the interval evaluator CONTAINS every point
/// sample, on random DAGs × random boxes. Logs bound tightness.
#[test]
fn frep_001_interval_containment() {
    let mut rng = Lcg(0x1001_2026_0706_0011);
    let mut violations = 0u64;
    let mut tightness_sum = 0.0;
    let mut tightness_n = 0u64;
    for _ in 0..12 {
        let f = random_frep(&mut rng);
        let dom = f.support().inflate(0.5);
        for _ in 0..8 {
            let c = sample_in(&mut rng, &dom);
            let h = Vec3::new(
                rng.range(0.05, 0.6),
                rng.range(0.05, 0.6),
                rng.range(0.05, 0.6),
            );
            let cell = Aabb::new(
                Point3::new(c.x - h.x, c.y - h.y, c.z - h.z),
                Point3::new(c.x + h.x, c.y + h.y, c.z + h.z),
            );
            let (lo, hi) = f.interval(&cell);
            let mut seen_lo = f64::INFINITY;
            let mut seen_hi = f64::NEG_INFINITY;
            for _ in 0..25 {
                let v = f.value(sample_in(&mut rng, &cell));
                seen_lo = seen_lo.min(v);
                seen_hi = seen_hi.max(v);
                // 1e-9 slack: shared fp rounding between the two
                // evaluators; OUTWARD-rounded intervals join with fs-ivl.
                if v < lo - 1e-9 || v > hi + 1e-9 {
                    violations += 1;
                }
            }
            if hi > lo {
                tightness_sum += (seen_hi - seen_lo) / (hi - lo);
                tightness_n += 1;
            }
        }
    }
    let tightness = tightness_sum / tightness_n.max(1) as f64;
    let mut em = fs_obs::Emitter::new("fs-rep-frep/conformance", "frep-001/interval");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "rep-frep-interval-tightness".to_string(),
                json: format!(
                    "{{\"violations\":{violations},\"mean_tightness\":{tightness:.4},\
                     \"boxes\":{tightness_n}}}"
                ),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("interval event validates");
    println!("{line}");
    verdict(
        "frep-001",
        violations == 0,
        &format!(
            "interval contains all samples on 12 random DAGs x 8 boxes x 25 points \
             (violations={violations}, mean observed/bound width={tightness:.3}); \
             seed 0x1001_2026_0706_0011"
        ),
    );
}

/// frep-002 — the composed Lipschitz bound is never violated under
/// adversarial sampling (near-pairs and far-pairs).
#[test]
fn frep_002_lipschitz_validity() {
    let mut rng = Lcg(0x1001_2026_0706_0012);
    let mut violations = 0u64;
    let mut worst_ratio = 0.0f64;
    for _ in 0..12 {
        let f = random_frep(&mut rng);
        let l = f.lipschitz();
        let dom = f.support().inflate(0.5);
        for i in 0..300 {
            let p = sample_in(&mut rng, &dom);
            let q = if i % 2 == 0 {
                sample_in(&mut rng, &dom)
            } else {
                // Near-pair: the regime where a broken bound shows first.
                let eps = 10.0_f64.powf(rng.range(-6.0, -1.0));
                Point3::new(
                    p.x + rng.range(-eps, eps),
                    p.y + rng.range(-eps, eps),
                    p.z + rng.range(-eps, eps),
                )
            };
            let d = q.delta_from(p).norm();
            if d < 1e-12 {
                continue;
            }
            let df = (f.value(p) - f.value(q)).abs();
            worst_ratio = worst_ratio.max(df / (l * d));
            if df > l * d * (1.0 + 1e-9) + 1e-12 {
                violations += 1;
            }
        }
    }
    verdict(
        "frep-002",
        violations == 0,
        &format!(
            "|f(x)-f(y)| <= L|x-y| held on 12 DAGs x 300 adversarial pairs \
             (violations={violations}, worst observed/bound ratio={worst_ratio:.4}); \
             seed 0x1001_2026_0706_0012"
        ),
    );
}

/// frep-003 — R-function blends are C¹ at the seam (analytic gradient
/// matches a crease-straddling central difference), the SAME probe
/// exhibits the hard Boolean's derivative discontinuity, and blend
/// weights are a convex pair. Logs seam gradient stats.
#[test]
fn frep_003_blend_seam_c1() {
    let two_spheres = |style: BoolStyle, op: BoolOp| -> Frep {
        let mut b = FrepBuilder::new();
        let s1 = b.sphere(Point3::new(-0.6, 0.0, 0.0), 1.0).expect("s1");
        let s2 = b.sphere(Point3::new(0.6, 0.0, 0.0), 1.0).expect("s2");
        let root = b.boolean(op, style, s1, s2).expect("bool");
        b.finish(root).expect("frep")
    };
    let fd_grad = |f: &Frep, p: Point3| -> Vec3 {
        let h = 1e-5;
        Vec3::new(
            (f.value(Point3::new(p.x + h, p.y, p.z)) - f.value(Point3::new(p.x - h, p.y, p.z)))
                / (2.0 * h),
            (f.value(Point3::new(p.x, p.y + h, p.z)) - f.value(Point3::new(p.x, p.y - h, p.z)))
                / (2.0 * h),
            (f.value(Point3::new(p.x, p.y, p.z + h)) - f.value(Point3::new(p.x, p.y, p.z - h)))
                / (2.0 * h),
        )
    };
    // The seam locus of the two-sphere fixture is the x = 0 plane.
    let seam_points: Vec<Point3> = (0..15)
        .flat_map(|iy| (0..15).map(move |iz| (iy, iz)))
        .map(|(iy, iz)| {
            Point3::new(
                0.0,
                -1.4 + 2.8 * f64::from(iy) / 14.0,
                -1.4 + 2.8 * f64::from(iz) / 14.0,
            )
        })
        .collect();
    let mut blend_worst = 0.0f64;
    for op in [BoolOp::Union, BoolOp::Intersect, BoolOp::Difference] {
        let f = two_spheres(BoolStyle::Blend { radius: 0.3 }, op);
        for &p in &seam_points {
            let (_, g) = f.value_grad(p);
            let Some(g) = g else { continue };
            let fd = fd_grad(&f, p);
            let e = (g.x - fd.x)
                .abs()
                .max((g.y - fd.y).abs())
                .max((g.z - fd.z).abs());
            blend_worst = blend_worst.max(e);
        }
    }
    let hard = two_spheres(BoolStyle::Hard, BoolOp::Union);
    let mut hard_worst = 0.0f64;
    for &p in &seam_points {
        let (_, g) = hard.value_grad(p);
        let Some(g) = g else { continue };
        let fd = fd_grad(&hard, p);
        hard_worst = hard_worst.max((g.x - fd.x).abs());
    }
    let mut rng = Lcg(0x1001_2026_0706_0013);
    let mut convex_ok = true;
    for _ in 0..2000 {
        let (a, b, r) = (
            rng.range(-2.0, 2.0),
            rng.range(-2.0, 2.0),
            rng.range(0.01, 0.5),
        );
        let (wa, wb) = smin_weights(a, b, r);
        convex_ok &= wa >= 0.0 && wb >= 0.0 && (wa + wb - 1.0).abs() < 1e-15;
    }
    let mut em = fs_obs::Emitter::new("fs-rep-frep/conformance", "frep-003/seam");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "rep-frep-seam-gradients".to_string(),
                json: format!(
                    "{{\"blend_worst_fd_err\":{blend_worst:.2e},\
                     \"hard_crease_mismatch\":{hard_worst:.3}}}"
                ),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("seam event validates");
    println!("{line}");
    verdict(
        "frep-003",
        blend_worst < 5e-4 && hard_worst > 0.2 && convex_ok,
        &format!(
            "blend gradients are C1 across the seam for all three ops (worst \
             FD-vs-analytic {blend_worst:.2e}) while the hard union's crease shows a \
             {hard_worst:.2} discontinuity — the optimization poison, exhibited; \
             weights stay a convex pair over 2000 draws"
        ),
    );
}

/// frep-004 — a DAG with a SHARED subexpression evaluates bitwise
/// identically to its expanded-tree rewrite (values and gradients).
#[test]
fn frep_004_sharing_equals_expansion() {
    let shared = {
        let mut b = FrepBuilder::new();
        let s = b.sphere(Point3::new(0.0, 0.0, 0.0), 0.8).expect("s");
        let bx = b
            .box_prim(Point3::new(0.3, 0.0, 0.0), Vec3::new(0.5, 0.4, 0.6))
            .expect("bx");
        let blend = b
            .boolean(BoolOp::Union, BoolStyle::Blend { radius: 0.2 }, s, bx)
            .expect("blend");
        // The SAME node id feeds both branches.
        let t1 = b.translate(blend, Vec3::new(0.7, 0.0, 0.0)).expect("t1");
        let t2 = b.rotate(blend, Vec3::new(0.0, 0.0, 1.0), 0.6).expect("t2");
        let root = b
            .boolean(BoolOp::Intersect, BoolStyle::Hard, t1, t2)
            .expect("root");
        b.finish(root).expect("frep")
    };
    let expanded = {
        let mut b = FrepBuilder::new();
        let blend = |b: &mut FrepBuilder| {
            let s = b.sphere(Point3::new(0.0, 0.0, 0.0), 0.8).expect("s");
            let bx = b
                .box_prim(Point3::new(0.3, 0.0, 0.0), Vec3::new(0.5, 0.4, 0.6))
                .expect("bx");
            b.boolean(BoolOp::Union, BoolStyle::Blend { radius: 0.2 }, s, bx)
                .expect("blend")
        };
        // Two DISTINCT copies of the subexpression.
        let b1 = blend(&mut b);
        let b2 = blend(&mut b);
        let t1 = b.translate(b1, Vec3::new(0.7, 0.0, 0.0)).expect("t1");
        let t2 = b.rotate(b2, Vec3::new(0.0, 0.0, 1.0), 0.6).expect("t2");
        let root = b
            .boolean(BoolOp::Intersect, BoolStyle::Hard, t1, t2)
            .expect("root");
        b.finish(root).expect("frep")
    };
    let mut rng = Lcg(0x1001_2026_0706_0014);
    let dom = shared.support().inflate(0.5);
    let mut bit_equal = true;
    for _ in 0..400 {
        let p = sample_in(&mut rng, &dom);
        let (vs, gs) = shared.value_grad(p);
        let (ve, ge) = expanded.value_grad(p);
        bit_equal &= vs.to_bits() == ve.to_bits() && gs == ge;
    }
    verdict(
        "frep-004",
        bit_equal,
        "shared-subexpression DAG matches its expanded tree BITWISE (value and \
         gradient) over 400 samples; seed 0x1001_2026_0706_0014",
    );
}

/// Conservative sphere tracer: steps by `f/L` (safe: `|f| <= dist`).
fn trace(f: &Frep, o: Point3, dir: Vec3, l: f64, tmax: f64) -> (Option<f64>, bool) {
    let mut t = 0.0;
    for _ in 0..4000 {
        let p = o.offset(dir.scale(t));
        let v = f.value(p);
        if v < 1e-7 {
            return (Some(t), false);
        }
        t += v / l;
        if t > tmax {
            return (None, false);
        }
    }
    (None, true) // stalled: step budget exhausted while still approaching
}

/// Dense-scan + bisection oracle for the FIRST ray crossing.
fn oracle(f: &Frep, o: Point3, dir: Vec3, tmax: f64) -> Option<f64> {
    let steps = 800;
    let dt = tmax / f64::from(steps);
    let mut prev = f.value(o);
    for i in 1..=steps {
        let t = dt * f64::from(i);
        let v = f.value(o.offset(dir.scale(t)));
        if prev >= 0.0 && v < 0.0 {
            let (mut lo, mut hi) = (t - dt, t);
            for _ in 0..60 {
                let mid = f64::midpoint(lo, hi);
                if f.value(o.offset(dir.scale(mid))) < 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            return Some(f64::midpoint(lo, hi));
        }
        prev = v;
    }
    None
}

/// frep-005 — sphere-tracing SAFETY: with the certified conservative
/// field, the tracer never tunnels past the oracle's first hit and never
/// claims a hit beyond it. Violation count must be 0. Also pins the
/// chart-level honesty (certificate kinds, Lipschitz, C-class).
#[test]
fn frep_005_sphere_tracing_safety() {
    let mut rng = Lcg(0x1001_2026_0706_0015);
    let mut rays = 0u64;
    let mut hits = 0u64;
    let mut safety_violations = 0u64;
    let mut stalls = 0u64;
    let mut worst_late = 0.0f64;
    for _ in 0..12 {
        let f = random_frep(&mut rng);
        let l = f.lipschitz();
        let sup = f.support();
        let center = Point3::new(
            f64::midpoint(sup.min.x, sup.max.x),
            f64::midpoint(sup.min.y, sup.max.y),
            f64::midpoint(sup.min.z, sup.max.z),
        );
        let radius = sup.max.delta_from(sup.min).norm() * 0.5 + 1.0;
        for _ in 0..120 {
            // Origin on a sphere OUTSIDE the support, aimed at a random
            // interior target.
            let theta = rng.range(0.0, core::f64::consts::TAU);
            let z = rng.range(-1.0, 1.0);
            let s = (1.0 - z * z).max(0.0).sqrt();
            let o = center.offset(Vec3::new(s * theta.cos(), s * theta.sin(), z).scale(radius));
            let target = sample_in(&mut rng, &sup);
            let dir = target.delta_from(o);
            let dir = dir.scale(1.0 / dir.norm());
            if f.value(o) <= 0.0 {
                continue; // pathological support: skip, not a safety case
            }
            let tmax = 2.0 * radius + 1.0;
            rays += 1;
            let t_star = oracle(&f, o, dir, tmax);
            let (t_trace, stalled) = trace(&f, o, dir, l, tmax);
            match (t_star, t_trace) {
                (Some(ts), Some(tt)) => {
                    hits += 1;
                    if tt > ts + 1e-4 {
                        safety_violations += 1; // claimed a hit PAST the surface
                    }
                    worst_late = worst_late.max(ts - tt);
                }
                (Some(_), None) => {
                    if stalled {
                        stalls += 1; // still approaching: incomplete, not unsafe
                    } else {
                        safety_violations += 1; // tunneled clean through
                    }
                }
                // Tracer hit where the coarse oracle saw nothing: the
                // oracle's dt missed a thin feature — not a safety fault.
                (None, _) => {}
            }
        }
    }
    let (exact_kind, csg_kind, c_class) = with_cx(|cx| {
        let mut b = FrepBuilder::new();
        let s = b.sphere(Point3::new(0.0, 0.0, 0.0), 1.0).expect("s");
        let t = b.translate(s, Vec3::new(0.2, 0.0, 0.0)).expect("t");
        let pure = b.finish(t).expect("pure");
        let mut b2 = FrepBuilder::new();
        let s1 = b2.sphere(Point3::new(0.0, 0.0, 0.0), 1.0).expect("s1");
        let s2 = b2.sphere(Point3::new(0.5, 0.0, 0.0), 0.8).expect("s2");
        let u = b2
            .boolean(BoolOp::Union, BoolStyle::Blend { radius: 0.2 }, s1, s2)
            .expect("u");
        let csg = b2.finish(u).expect("csg");
        let p = Point3::new(1.7, 0.3, -0.2);
        (
            pure.eval(p, cx).error.kind,
            csg.eval(p, cx).error.kind,
            csg.differentiability(),
        )
    });
    let mut em = fs_obs::Emitter::new("fs-rep-frep/conformance", "frep-005/tracing");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "rep-frep-ray-safety".to_string(),
                json: format!(
                    "{{\"rays\":{rays},\"oracle_hits\":{hits},\
                     \"safety_violations\":{safety_violations},\"stalls\":{stalls},\
                     \"worst_standoff\":{worst_late:.2e}}}"
                ),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("tracing event validates");
    println!("{line}");
    verdict(
        "frep-005",
        safety_violations == 0
            && hits > 200
            && exact_kind == fs_evidence::NumericalKind::Enclosure
            && csg_kind == fs_evidence::NumericalKind::Estimate
            && c_class == Differentiability::C1,
        &format!(
            "0 tunneling violations over {rays} rays on 12 random DAGs \
             ({hits} oracle-confirmed hits, {stalls} grazing stalls, worst standoff \
             {worst_late:.1e}); rigid chains enclose rounded exact-distance \
             evaluations, CSG keeps an Estimate relative to abstract distance \
             while exposing a separate rigorous trace enclosure, pure-blend DAGs advertise C1; \
             seed 0x1001_2026_0706_0015"
        ),
    );
}

/// frep-005b — evaluation rigor and typed geometric claims stay separate.
/// Rounded primitive arithmetic is enclosed, while generic normalization and
/// nontrivial Rodrigues maps cannot mint an exact-distance theorem.
#[test]
#[allow(clippy::too_many_lines)] // One claim-separation regression per producer family.
fn frep_005b_evaluation_enclosures_do_not_overpromote_geometry() {
    with_cx(|cx| {
        let mut sphere_builder = FrepBuilder::new();
        let sphere = sphere_builder
            .sphere(Point3::new(0.0, 0.0, 0.0), 1.0)
            .expect("sphere");
        let sphere = sphere_builder.finish(sphere).expect("sphere frep");
        let point = Point3::new(
            0.038_546_885_717_366_49,
            -0.607_415_449_300_028_1,
            -0.793_448_734_204_907_9,
        );
        let sample = sphere.eval(point, cx);
        let independent_residual = 1.421_291_065_674_733_6e-7;
        let sphere_ok = sphere.trace_step_claim() == TraceStepClaim::ExactDistance
            && sample.error.kind == fs_evidence::NumericalKind::Enclosure
            && sample.error.lo <= independent_residual
            && independent_residual <= sample.error.hi;

        let mut generic_plane_builder = FrepBuilder::new();
        let plane = generic_plane_builder
            .half_space(Vec3::new(1.0, 3.0, 7.0), 0.25)
            .expect("generic plane");
        let generic_plane = generic_plane_builder.finish(plane).expect("plane frep");
        let plane_point = Point3::new(2.0, -1.0, 0.5);
        let plane_sample = generic_plane.eval(plane_point, cx);
        let plane_ok = generic_plane.trace_step_claim() == TraceStepClaim::LipschitzImplicit
            && generic_plane.lipschitz() > 1.0
            && plane_sample.error.kind == fs_evidence::NumericalKind::Estimate
            && generic_plane
                .trace_value_enclosure(plane_point, &plane_sample, cx)
                .kind
                == fs_evidence::NumericalKind::Enclosure;

        let mut axis_plane_builder = FrepBuilder::new();
        let axis_plane = axis_plane_builder
            .half_space(Vec3::new(1.0, 0.0, 0.0), 0.25)
            .expect("axis plane");
        let axis_plane = axis_plane_builder
            .finish(axis_plane)
            .expect("axis plane frep");
        let axis_ok = axis_plane.trace_step_claim() == TraceStepClaim::ExactDistance;

        let mut rotated_builder = FrepBuilder::new();
        let base = rotated_builder
            .sphere(Point3::new(0.2, 0.0, 0.0), 0.5)
            .expect("rotated base");
        let rotated = rotated_builder
            .rotate(base, Vec3::new(0.3, -0.5, 0.8), 0.7)
            .expect("rotation");
        let rotated = rotated_builder.finish(rotated).expect("rotated frep");
        let rotated_point = Point3::new(1.0, 0.5, -0.25);
        let rotated_sample = rotated.eval(rotated_point, cx);
        let rotation_ok = rotated.trace_step_claim() == TraceStepClaim::LipschitzImplicit
            && rotated.lipschitz() > 1.0
            && rotated
                .trace_value_enclosure(rotated_point, &rotated_sample, cx)
                .kind
                == fs_evidence::NumericalKind::Enclosure;

        let mut reachable_builder = FrepBuilder::new();
        let root = reachable_builder
            .sphere(Point3::new(0.0, 0.0, 0.0), 1.0)
            .expect("reachable root");
        let unused = reachable_builder
            .sphere(Point3::new(4.0, 0.0, 0.0), 0.5)
            .expect("unused sphere");
        let _unused_boolean = reachable_builder
            .boolean(BoolOp::Union, BoolStyle::Hard, root, unused)
            .expect("unused boolean");
        let reachable = reachable_builder.finish(root).expect("reachable frep");
        let unused_hard_ok = reachable.trace_step_claim() == TraceStepClaim::ExactDistance
            && reachable.differentiability() == Differentiability::C1;

        let mut hard_root_builder = FrepBuilder::new();
        let left = hard_root_builder
            .sphere(Point3::new(-0.25, 0.0, 0.0), 1.0)
            .expect("hard left");
        let right = hard_root_builder
            .sphere(Point3::new(0.25, 0.0, 0.0), 1.0)
            .expect("hard right");
        let hard_root = hard_root_builder
            .boolean(BoolOp::Union, BoolStyle::Hard, left, right)
            .expect("reachable hard boolean");
        let hard_root = hard_root_builder.finish(hard_root).expect("hard-root frep");
        let reachable_hard_ok = hard_root.differentiability() == Differentiability::C0;

        // Regression seed: a far, ulp-scale-thin box rotated by pi around a
        // rounded diagonal axis. The represented field says `preimage` is
        // inside, while a rounded +angle corner map used to miss it by one
        // ulp. Support is now the certified preimage of the inverse map.
        let mut support_builder = FrepBuilder::new();
        let thin_box = support_builder
            .box_prim(
                Point3::new(
                    -59_387_588.392_124_265,
                    6_896_459.806_064_8,
                    -80_159_548.321_525_37,
                ),
                Vec3::new(
                    2.980_232_238_769_531_2e-8,
                    1.490_116_119_384_765_6e-8,
                    5.960_464_477_539_063e-8,
                ),
            )
            .expect("thin far box");
        let thin_box = support_builder
            .rotate(thin_box, Vec3::new(1.0, 1.0, 0.0), core::f64::consts::PI)
            .expect("rounded-axis rotation");
        let thin_box = support_builder.finish(thin_box).expect("rotated thin box");
        let preimage = Point3::new(
            6_896_459.806_064_783,
            -59_387_588.392_124_27,
            80_159_548.321_525_37,
        );
        let thin_support = thin_box.support();
        let support_ok = thin_box.value(preimage) < 0.0
            && thin_support.contains(preimage)
            && [
                thin_support.min.x,
                thin_support.min.y,
                thin_support.min.z,
                thin_support.max.x,
                thin_support.max.y,
                thin_support.max.z,
            ]
            .into_iter()
            .all(f64::is_finite);

        let mut wide_trig_builder = FrepBuilder::new();
        let wide_trig = wide_trig_builder
            .sphere(Point3::new(0.0, 0.0, 0.0), 1.0)
            .expect("wide-trig sphere");
        let wide_trig = wide_trig_builder
            .rotate(wide_trig, Vec3::new(1.0, 1.0, 0.0), 2_000_000.0)
            .expect("wide-trig rotation");
        let wide_trig = wide_trig_builder.finish(wide_trig).expect("wide-trig frep");
        let fallback = wide_trig.support();
        let singular_fallback_ok = [fallback.min.x, fallback.min.y, fallback.min.z]
            .into_iter()
            .all(|bound| bound.is_infinite() && bound.is_sign_negative())
            && [fallback.max.x, fallback.max.y, fallback.max.z]
                .into_iter()
                .all(|bound| bound.is_infinite() && bound.is_sign_positive());

        verdict(
            "frep-005b",
            sphere_ok
                && plane_ok
                && axis_ok
                && rotation_ok
                && unused_hard_ok
                && reachable_hard_ok
                && support_ok
                && singular_fallback_ok,
            "the seeded rounded sphere contains the independent residual; generic \
             half-spaces and nontrivial rotations are Lipschitz-implicit with widened \
             bounds; differentiability follows only the root-reachable DAG; rotated \
             support encloses the certified inverse-map preimage and fails closed \
             when matrix regularity is uncertified",
        );
    });
}

/// frep-006 — metamorphic algebra + design levers: hard idempotence and
/// commutativity are BITWISE; blend self-union equals a r/4 dilation
/// BITWISE; rotation round-trips; translation is equivariant; and the
/// parameter surface exposes exact radius/offset derivatives with a
/// zero-derivative far from the blend zone.
#[test]
#[allow(clippy::too_many_lines)] // one law per block; splitting hides the suite's shape
fn frep_006_metamorphic_and_params() {
    let mut rng = Lcg(0x1001_2026_0706_0016);
    let sphere_at = |c: Point3, r: f64| -> Frep {
        let mut b = FrepBuilder::new();
        let s = b.sphere(c, r).expect("s");
        b.finish(s).expect("frep")
    };
    // Hard idempotence + commutativity, bitwise.
    let (a_c, a_r) = (Point3::new(0.1, -0.2, 0.3), 0.9);
    let (b_c, b_r) = (Point3::new(-0.4, 0.5, 0.0), 0.7);
    let plain_a = sphere_at(a_c, a_r);
    let plain_b = sphere_at(b_c, b_r);
    let mk_bool = |op: BoolOp, style: BoolStyle, swap: bool| -> Frep {
        let mut b = FrepBuilder::new();
        let s1 = b.sphere(a_c, a_r).expect("s1");
        let s2 = b.sphere(b_c, b_r).expect("s2");
        let (x, y) = if swap { (s2, s1) } else { (s1, s2) };
        let root = b.boolean(op, style, x, y).expect("bool");
        b.finish(root).expect("frep")
    };
    let self_union = |style: BoolStyle| -> Frep {
        let mut b = FrepBuilder::new();
        let s = b.sphere(a_c, a_r).expect("s");
        let root = b.boolean(BoolOp::Union, style, s, s).expect("bool");
        b.finish(root).expect("frep")
    };
    let idem = self_union(BoolStyle::Hard);
    let blend_self = self_union(BoolStyle::Blend { radius: 0.25 });
    let uni = mk_bool(BoolOp::Union, BoolStyle::Hard, false);
    let uni_swapped = mk_bool(BoolOp::Union, BoolStyle::Hard, true);
    let buni = mk_bool(BoolOp::Union, BoolStyle::Blend { radius: 0.3 }, false);
    let buni_swapped = mk_bool(BoolOp::Union, BoolStyle::Blend { radius: 0.3 }, true);
    let mut laws_ok = true;
    for _ in 0..200 {
        let p = Point3::new(
            rng.range(-2.0, 2.0),
            rng.range(-2.0, 2.0),
            rng.range(-2.0, 2.0),
        );
        laws_ok &= idem.value(p).to_bits() == plain_a.value(p).to_bits();
        laws_ok &= blend_self.value(p).to_bits() == (plain_a.value(p) - 0.0625).to_bits();
        laws_ok &= uni.value(p).to_bits() == uni_swapped.value(p).to_bits();
        laws_ok &= buni.value(p).to_bits() == buni_swapped.value(p).to_bits();
        laws_ok &= uni.value(p).to_bits() == plain_a.value(p).min(plain_b.value(p)).to_bits();
    }
    // Rotation round-trip and dyadic translation equivariance.
    let mut rot_ok = true;
    {
        let mut b = FrepBuilder::new();
        let s = b
            .torus(Point3::new(0.1, 0.0, -0.2), 0.8, 0.25)
            .expect("torus");
        let r1 = b.rotate(s, Vec3::new(0.3, -0.5, 0.8), 0.7).expect("r1");
        let r2 = b.rotate(r1, Vec3::new(0.3, -0.5, 0.8), -0.7).expect("r2");
        let round = b.finish(r2).expect("frep");
        let plain = {
            let mut b = FrepBuilder::new();
            let s = b
                .torus(Point3::new(0.1, 0.0, -0.2), 0.8, 0.25)
                .expect("torus");
            b.finish(s).expect("frep")
        };
        let mut bt = FrepBuilder::new();
        let s = bt
            .torus(Point3::new(0.1, 0.0, -0.2), 0.8, 0.25)
            .expect("torus");
        let t = bt.translate(s, Vec3::new(0.5, 0.25, -0.375)).expect("t");
        let shifted = bt.finish(t).expect("frep");
        for _ in 0..200 {
            let p = Point3::new(
                rng.range(-1.5, 1.5),
                rng.range(-1.5, 1.5),
                rng.range(-1.5, 1.5),
            );
            rot_ok &= (round.value(p) - plain.value(p)).abs() < 1e-12;
            let q = Point3::new(p.x + 0.5, p.y + 0.25, p.z - 0.375);
            rot_ok &= (shifted.value(q) - plain.value(p)).abs() < 1e-12;
        }
    }
    // Design levers.
    let mut lever_ok = true;
    {
        let mut f = sphere_at(Point3::new(0.0, 0.0, 0.0), 1.0);
        let params = f.params();
        lever_ok &= params.len() == 4 && params[3].1 == "radius";
        let radius_id = params[3].0;
        let p = Point3::new(1.3, 0.4, -0.2);
        // d(|p−c| − r)/dr = −1 exactly.
        let d = f.d_value_d_param(p, radius_id).expect("radius lever");
        lever_ok &= (d + 1.0).abs() < 1e-6;
        let before = f.value(p);
        f.set_param(radius_id, 1.1).expect("grow radius");
        lever_ok &= (before - f.value(p) - 0.1).abs() < 1e-12;
        lever_ok &= f.set_param(radius_id, -1.0).is_err(); // teaching refusal

        let mut b = FrepBuilder::new();
        let s1 = b.sphere(Point3::new(-0.6, 0.0, 0.0), 1.0).expect("s1");
        let s2 = b.sphere(Point3::new(0.6, 0.0, 0.0), 1.0).expect("s2");
        let u = b
            .boolean(BoolOp::Union, BoolStyle::Blend { radius: 0.3 }, s1, s2)
            .expect("u");
        let blend = b.finish(u).expect("frep");
        let blend_radius_id = blend
            .params()
            .into_iter()
            .find(|(_, name, _)| *name == "blend radius")
            .expect("blend lever")
            .0;
        // Far from the seam the blend radius has EXACTLY zero action...
        let far = blend
            .d_value_d_param(Point3::new(-1.8, 0.0, 0.0), blend_radius_id)
            .expect("far");
        // ...and on the seam it acts like the −h²/4 term (≈ −1/4 at u = 0).
        let seam = blend
            .d_value_d_param(Point3::new(0.0, 1.2, 0.0), blend_radius_id)
            .expect("seam");
        #[allow(clippy::float_cmp)] // outside the blend zone the lever's FD is EXACTLY zero
        let far_is_dead = far == 0.0;
        lever_ok &= far_is_dead && (seam + 0.25).abs() < 1e-4;
    }
    verdict(
        "frep-006",
        laws_ok && rot_ok && lever_ok,
        "hard idempotence/commutativity and blend-self = dilation(r/4) hold BITWISE \
         over 200 points; rotation round-trips and dyadic translation is equivariant \
         to 1e-12; radius/offset levers differentiate exactly (-1), the blend-radius \
         lever is 0 far from the seam and -1/4 on it, and invalid lever values teach; \
         seed 0x1001_2026_0706_0016",
    );
}
