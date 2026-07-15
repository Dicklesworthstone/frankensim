//! Chart-backend conformance (beads qfx.2 + 8ll9; default-on).
//! Acceptance: ZERO missed intersections across the
//! adversarial ray battery (thin shells, grazing rays) vs the certified
//! root-finder oracle — the headline certificate test; NURBS Newton
//! matches analytic references; mixed-chart scenes agree across backend
//! kinds within tolerance; G3 frame invariance; ray rates measured and
//! ledgered (debug-build numbers; the perf GATE lives in perf-CI).
#![cfg(feature = "chart-backends")]

use core::sync::atomic::{AtomicUsize, Ordering};

use asupersync::types::Budget;
use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::{SphereChart, TorusChart};
use fs_geom::{Aabb, Chart, ChartSample, Point3, TraceStepClaim, Vec3};
use fs_math::eft::two_sum;
use fs_render::charts::{
    Backend, CHART_BACKEND_BIT_SEMANTICS_VERSION, NurbsRayError, Ray, SceneTraceError,
    TraceTermination, TriMesh, ray_intersect_nurbs, sphere_trace, trace_scene,
};
use fs_rep_frep::{BoolOp, BoolStyle, Frep, FrepBuilder};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-render/charts\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 9,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn unit(v: Vec3) -> Vec3 {
    v.scale(1.0 / v.norm())
}

fn exact_subtraction_bounds(lhs: f64, rhs: f64) -> (f64, f64) {
    let (rounded, tail) = two_sum(lhs, -rhs);
    if tail > 0.0 {
        (rounded, rounded.next_up())
    } else if tail < 0.0 {
        (rounded.next_down(), rounded)
    } else {
        (rounded, rounded)
    }
}

fn bounded_certificate(lo: f64, hi: f64, nominal: f64) -> NumericalCertificate {
    let lo = lo.min(nominal);
    let hi = hi.max(nominal);
    if lo == hi {
        NumericalCertificate::exact(lo)
    } else {
        NumericalCertificate::enclosure(lo, hi)
    }
}

fn lcg(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / (1u64 << 53) as f64
}

/// The THIN SHELL adversary: a sphere shell of thickness 2e-3 built as
/// an F-rep difference — the classic tunneling victim.
fn thin_shell() -> fs_rep_frep::Frep {
    thin_shell_with_thickness(2e-3)
}

fn thin_shell_with_thickness(thickness: f64) -> fs_rep_frep::Frep {
    let mut b = FrepBuilder::new();
    let outer = b.sphere(Point3::new(0.0, 0.0, 0.0), 1.0).expect("outer");
    let inner = b
        .sphere(Point3::new(0.0, 0.0, 0.0), 1.0 - thickness)
        .expect("inner");
    let shell = b
        .boolean(BoolOp::Difference, BoolStyle::Hard, outer, inner)
        .expect("shell");
    b.finish(shell).expect("frep")
}

/// Test chart for fail-closed audit states. The value/bound under test and its
/// point estimate are overridden together so unrelated certificate mismatch
/// does not mask the intended termination state.
struct SampleOverrideChart {
    inner: Frep,
    distance: f64,
    lipschitz: Option<f64>,
    trace_claim: TraceStepClaim,
    rigorous: bool,
}

impl Chart for SampleOverrideChart {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        let mut sample = self.inner.eval(x, cx);
        sample.signed_distance = self.distance;
        sample.gradient = None;
        sample.lipschitz = self.lipschitz;
        sample.error = NumericalCertificate::estimate(self.distance, self.distance);
        sample
    }

    fn support(&self) -> fs_geom::Aabb {
        self.inner.support()
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        self.trace_claim
    }

    fn trace_value_enclosure(
        &self,
        _x: Point3,
        _sample: &ChartSample,
        _cx: &Cx<'_>,
    ) -> NumericalCertificate {
        if self.rigorous && self.distance.is_finite() {
            NumericalCertificate::enclosure(self.distance.next_down(), self.distance.next_up())
        } else {
            NumericalCertificate::estimate(self.distance, self.distance)
        }
    }

    fn name(&self) -> &'static str {
        "trace-audit-override"
    }
}

/// Multiplies an F-rep field and its certified Lipschitz bound by the same
/// positive scale. The zero set is unchanged, but a marcher that silently
/// assumes `L = 1` can leap completely across the thin shell.
struct ScaledChart {
    inner: Frep,
    scale: f64,
    publish_bound: bool,
    trace_claim: TraceStepClaim,
}

impl Chart for ScaledChart {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        let mut sample = self.inner.eval(x, cx);
        sample.signed_distance *= self.scale;
        sample.gradient = sample.gradient.map(|gradient| gradient.scale(self.scale));
        sample.lipschitz = self.publish_bound.then(|| {
            let bound = sample.lipschitz.expect("F-rep publishes a bound") * self.scale;
            if bound.is_finite() {
                bound.next_up()
            } else {
                bound
            }
        });
        sample.error =
            NumericalCertificate::estimate(sample.signed_distance, sample.signed_distance);
        sample
    }

    fn support(&self) -> fs_geom::Aabb {
        self.inner.support()
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        self.trace_claim
    }

    fn trace_value_enclosure(
        &self,
        x: Point3,
        sample: &ChartSample,
        cx: &Cx<'_>,
    ) -> NumericalCertificate {
        if !self.scale.is_finite() || self.scale <= 0.0 {
            return NumericalCertificate::no_claim();
        }
        let inner_sample = self.inner.eval(x, cx);
        let inner = self.inner.trace_value_enclosure(x, &inner_sample, cx);
        if !matches!(inner.kind, NumericalKind::Exact | NumericalKind::Enclosure)
            || !inner.lo.is_finite()
            || !inner.hi.is_finite()
        {
            return NumericalCertificate::no_claim();
        }
        NumericalCertificate::enclosure(
            (inner.lo * self.scale)
                .next_down()
                .min(sample.signed_distance),
            (inner.hi * self.scale)
                .next_up()
                .max(sample.signed_distance),
        )
    }

    fn name(&self) -> &'static str {
        "scaled-thin-shell"
    }
}

/// Exact signed distance to an infinite x-aligned slab.
struct ExactSlabChart {
    lo: f64,
    hi: f64,
}

struct ConstantNoClaim;

impl Chart for ConstantNoClaim {
    fn eval(&self, _x: Point3, _cx: &Cx<'_>) -> ChartSample {
        ChartSample {
            signed_distance: 1.0,
            gradient: None,
            lipschitz: None,
            error: NumericalCertificate::estimate(1.0, 1.0),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "constant-no-claim"
    }
}

/// A valid Lipschitz-implicit field whose deliberately loose `L` makes the
/// normalized residual tiny one world unit from its zero set.
struct LooseLipschitzChart;

impl Chart for LooseLipschitzChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let value = 1e-12 * (x.x - 1.0);
        ChartSample {
            signed_distance: value,
            gradient: Some(Vec3::new(1e-12, 0.0, 0.0)),
            // The tight constant is 1e-12, so 1.0 is a valid but loose upper
            // bound and must not become geometric proximity evidence.
            lipschitz: Some(1.0),
            error: NumericalCertificate::estimate(value, value),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-2.0, -1.0, -1.0), Point3::new(2.0, 1.0, 1.0))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::LipschitzImplicit
    }

    fn trace_value_enclosure(
        &self,
        x: Point3,
        sample: &ChartSample,
        _cx: &Cx<'_>,
    ) -> NumericalCertificate {
        let delta = x.x - 1.0;
        NumericalCertificate::enclosure(
            (delta.next_down() * 1e-12)
                .next_down()
                .min(sample.signed_distance),
            (delta.next_up() * 1e-12)
                .next_up()
                .max(sample.signed_distance),
        )
    }

    fn name(&self) -> &'static str {
        "loose-lipschitz"
    }
}

/// An exact scalar plane evaluation that deliberately advertises only the
/// Lipschitz-implicit theorem. At `x = 0`, `[0, 0]` is still geometric zero-set
/// evidence even though nonzero residuals carry no proximity claim.
struct SingletonZeroImplicitChart;

impl Chart for SingletonZeroImplicitChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        ChartSample {
            signed_distance: x.x,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error: NumericalCertificate::estimate(x.x, x.x),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(0.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::LipschitzImplicit
    }

    fn trace_value_enclosure(
        &self,
        _x: Point3,
        sample: &ChartSample,
        _cx: &Cx<'_>,
    ) -> NumericalCertificate {
        NumericalCertificate::exact(sample.signed_distance)
    }

    fn name(&self) -> &'static str {
        "singleton-zero-implicit"
    }
}

/// A linear implicit field whose second evaluation requests cancellation. A
/// residual witness must observe that request before it can return a hit.
struct CancelOnWitnessChart<'a> {
    gate: &'a CancelGate,
    evaluations: AtomicUsize,
}

impl Chart for CancelOnWitnessChart<'_> {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        if self.evaluations.fetch_add(1, Ordering::SeqCst) == 1 {
            self.gate.request();
        }
        ChartSample {
            signed_distance: x.x,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error: NumericalCertificate::estimate(x.x, x.x),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(0.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::LipschitzImplicit
    }

    fn trace_value_enclosure(
        &self,
        _x: Point3,
        sample: &ChartSample,
        _cx: &Cx<'_>,
    ) -> NumericalCertificate {
        NumericalCertificate::exact(sample.signed_distance)
    }

    fn name(&self) -> &'static str {
        "cancel-on-witness"
    }
}

/// The current sample has a strict sign, while every nonnegative witness is
/// deliberately enclosed across zero. Containment remains rigorous, but it
/// cannot prove which side of the boundary the witness occupies.
struct StraddlingWitnessChart;

impl Chart for StraddlingWitnessChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        ChartSample {
            signed_distance: x.x,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error: NumericalCertificate::estimate(x.x, x.x),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(0.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::LipschitzImplicit
    }

    fn trace_value_enclosure(
        &self,
        _x: Point3,
        sample: &ChartSample,
        _cx: &Cx<'_>,
    ) -> NumericalCertificate {
        if sample.signed_distance < 0.0 {
            NumericalCertificate::exact(sample.signed_distance)
        } else {
            let radius = sample.signed_distance.abs().next_up();
            NumericalCertificate::enclosure(-radius, radius)
        }
    }

    fn name(&self) -> &'static str {
        "straddling-witness"
    }
}

/// An exact sphere deliberately advertised through the weaker implicit theorem
/// so a tangent ray exercises the no-sign-change boundary.
struct ImplicitSphereChart {
    inner: SphereChart,
}

impl Chart for ImplicitSphereChart {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        self.inner.eval(x, cx)
    }

    fn support(&self) -> Aabb {
        self.inner.support()
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::LipschitzImplicit
    }

    fn trace_value_enclosure(
        &self,
        x: Point3,
        sample: &ChartSample,
        cx: &Cx<'_>,
    ) -> NumericalCertificate {
        self.inner.trace_value_enclosure(x, sample, cx)
    }

    fn name(&self) -> &'static str {
        "implicit-sphere"
    }
}

struct CancelOnEvalChart<'a> {
    gate: &'a CancelGate,
}

impl Chart for CancelOnEvalChart<'_> {
    fn eval(&self, _x: Point3, _cx: &Cx<'_>) -> ChartSample {
        self.gate.request();
        ChartSample {
            signed_distance: 0.0,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error: NumericalCertificate::exact(0.0),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn name(&self) -> &'static str {
        "cancel-on-eval"
    }
}

impl Chart for ExactSlabChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let center = f64::midpoint(self.lo, self.hi);
        let half = 0.5 * (self.hi - self.lo);
        let distance = (x.x - center).abs() - half;
        let (lo, hi) = if x.x < self.lo {
            exact_subtraction_bounds(self.lo, x.x)
        } else if x.x > self.hi {
            exact_subtraction_bounds(x.x, self.hi)
        } else {
            let left = exact_subtraction_bounds(x.x, self.lo);
            let right = exact_subtraction_bounds(self.hi, x.x);
            let clearance = (left.0.min(right.0), left.1.min(right.1));
            (-clearance.1, -clearance.0)
        };
        ChartSample {
            signed_distance: distance,
            gradient: Some(Vec3::new((x.x - center).signum(), 0.0, 0.0)),
            lipschitz: Some(1.0),
            error: bounded_certificate(lo, hi, distance),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(
            Point3::new(self.lo, -1e6, -1e6),
            Point3::new(self.hi, 1e6, 1e6),
        )
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn name(&self) -> &'static str {
        "exact-thin-slab"
    }
}

/// Exact signed distance to the half-space x >= boundary.
struct ExactPlaneChart {
    boundary: f64,
    lipschitz: f64,
}

impl Chart for ExactPlaneChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let distance = self.boundary - x.x;
        let (lo, hi) = exact_subtraction_bounds(self.boundary, x.x);
        ChartSample {
            signed_distance: distance,
            gradient: Some(Vec3::new(-1.0, 0.0, 0.0)),
            lipschitz: Some(self.lipschitz),
            error: bounded_certificate(lo, hi, distance),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(
            Point3::new(self.boundary, -1e6, -1e6),
            Point3::new(1e6, 1e6, 1e6),
        )
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn name(&self) -> &'static str {
        "exact-plane"
    }
}

/// Deliberately wrong negative control: unit-Lipschitz marching ignores the
/// chart's published bound and therefore has no no-tunneling theorem.
fn naive_unit_bound_hits(chart: &dyn Chart, cx: &Cx<'_>, ray: &Ray, t_max: f64) -> bool {
    let mut t = 0.0;
    for _ in 0..128 {
        if t > t_max {
            return false;
        }
        let distance = chart.eval(ray.at(t), cx).signed_distance.abs();
        if distance <= 1e-6 {
            return true;
        }
        t += distance;
    }
    false
}

/// The certified ROOT-FINDER ORACLE: march at a fixed micro-step and bisect the
/// first sign-changing bracket — an independent earliest-root receipt.
fn oracle_first_hit(chart: &dyn Chart, cx: &Cx<'_>, ray: &Ray, t_max: f64) -> Option<f64> {
    let micro = 2e-4;
    let mut t = 0.0f64;
    let mut prev = chart.eval(ray.at(0.0), cx).signed_distance;
    while t < t_max {
        let previous_t = t;
        t = (t + micro).min(t_max);
        let d = chart.eval(ray.at(t), cx).signed_distance;
        if prev.signum() != d.signum() {
            let (mut lo, mut hi) = (previous_t, t);
            let mut flo = prev;
            for _ in 0..60 {
                let mid = f64::midpoint(lo, hi);
                let fmid = chart.eval(ray.at(mid), cx).signed_distance;
                if flo.signum() == fmid.signum() {
                    lo = mid;
                    flo = fmid;
                } else {
                    hi = mid;
                }
            }
            return Some(f64::midpoint(lo, hi));
        }
        prev = d;
    }
    None
}

#[test]
fn rb_001_zero_tunneling_headline() {
    with_cx(|cx| {
        assert_eq!(CHART_BACKEND_BIT_SEMANTICS_VERSION, 7);
        // Falsifier pairing: four thin-feature fields whose L > 1. The naive
        // unit-bound marcher tunnels in every case; the certified d/L path
        // approaches under a no-tunneling theorem and closes a short rigorous
        // sign bracket around the first crossing.
        let axial_ray = Ray {
            origin: Point3::new(0.0, 0.0, 3.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let fixtures = [(4e-3, 4.0), (2e-3, 8.0), (1e-3, 16.0), (5e-4, 32.0)];
        let mut naive_tunnels = 0usize;
        let mut certified_hits = 0usize;
        let mut residual_limits = 0usize;
        for (thickness, scale) in fixtures {
            let scaled = ScaledChart {
                inner: thin_shell_with_thickness(thickness),
                scale,
                publish_bound: true,
                trace_claim: TraceStepClaim::LipschitzImplicit,
            };
            let oracle_t = oracle_first_hit(&scaled, cx, &axial_ray, 6.0)
                .expect("the thin shell has a first crossing");
            if !naive_unit_bound_hits(&scaled, cx, &axial_ray, 6.0) {
                naive_tunnels += 1;
            }
            let (hit, audit) = sphere_trace(&scaled, cx, &axial_ray, 6.0, 1e-6, 1.0);
            let hit = hit.expect("a transverse short sign bracket certifies the first crossing");
            assert!(
                (hit.t - oracle_t).abs() <= 1e-6,
                "certified path must stay within tolerance of the earliest oracle root: {} vs {oracle_t}",
                hit.t
            );
            assert!(audit.certified, "published finite bound carries the claim");
            assert_eq!(audit.termination, TraceTermination::Hit);
            assert!(audit.certifies_hit());
            let normal = hit.normal.expect("scaled chart supplies a normal");
            assert!(
                (normal.norm() - 1.0).abs() <= 1e-12,
                "backend normal must be unit length, got {}",
                normal.norm()
            );
            certified_hits += 1;
        }
        // A chart that withholds the same bound may still use the preview
        // fallback, but the audit must refuse the certified headline.
        let unbounded = ScaledChart {
            inner: thin_shell(),
            scale: 8.0,
            publish_bound: false,
            trace_claim: TraceStepClaim::NoClaim,
        };
        let (_, unbounded_audit) = sphere_trace(&unbounded, cx, &axial_ray, 6.0, 1e-6, 1.0);
        assert!(
            !unbounded_audit.certified,
            "missing bound cannot mint a certificate"
        );
        assert_eq!(
            trace_scene(&[Backend::Chart(&unbounded)], cx, &axial_ray, 6.0, 1e-6),
            Err(SceneTraceError::UncertifiedTrace),
            "production scene tracing must preserve the uncertified refusal"
        );

        let shell = thin_shell();
        let mut state = 0xbeef_u64;
        let mut audited = 0usize;
        for k in 0..120 {
            // Adversarial battery: random origins on a far sphere aimed
            // near (but not always at) the shell — grazing rays included.
            let phi = lcg(&mut state) * std::f64::consts::TAU;
            let z = 2.0 * lcg(&mut state) - 1.0;
            let r = (1.0 - z * z).sqrt();
            let origin = Point3::new(3.0 * r * phi.cos(), 3.0 * r * phi.sin(), 3.0 * z);
            // Aim at a jittered point near the surface (grazing bias
            // every third ray).
            let graze = k % 3 == 0;
            let target_r = if graze { 0.999 } else { 0.6 * lcg(&mut state) };
            let tphi = lcg(&mut state) * std::f64::consts::TAU;
            let tz = 2.0 * lcg(&mut state) - 1.0;
            let tr = (1.0 - tz * tz).sqrt();
            let target = Point3::new(
                target_r * tr * tphi.cos(),
                target_r * tr * tphi.sin(),
                target_r * tz,
            );
            let ray = Ray {
                origin,
                dir: unit(target.delta_from(origin)),
            };
            let oracle = oracle_first_hit(&shell, cx, &ray, 6.0);
            let (traced, audit) = sphere_trace(&shell, cx, &ray, 6.0, 1e-6, 1.0);
            match (traced, audit.termination, oracle) {
                (Some(hit), TraceTermination::Hit, Some(oracle_t)) => {
                    assert!(
                        (hit.t - oracle_t).abs() <= 1e-6,
                        "ray {k}: bracketed hit must stay within tolerance of the earliest oracle root: {} vs {oracle_t}",
                        hit.t
                    );
                }
                (None, TraceTermination::ResidualLimit, None) => residual_limits += 1,
                (None, TraceTermination::Miss, None) => {}
                (traced, termination, oracle) => panic!(
                    "ray {k}: dishonest trace/oracle result: traced={} termination={termination:?} oracle={}",
                    traced.is_some(),
                    oracle.is_some()
                ),
            }
            // G0 step safety: no plain step ever exceeded |f|/L.
            assert!(
                audit.worst_step_ratio <= 1.0 + 1e-12,
                "step-safety property: {}",
                audit.worst_step_ratio
            );
            assert!(audit.certified);
            audited += 1;
        }
        println!(
            "{{\"metric\":\"tunneling-audit\",\"falsifiers\":{},\"naive_tunnels\":{naive_tunnels},\"certified_hits\":{certified_hits},\"residual_limits\":{residual_limits},\"oracle_rays\":{audited},\"geometric_false_hits\":0}}",
            fixtures.len()
        );
        assert_eq!(naive_tunnels, fixtures.len());
        assert_eq!(certified_hits, fixtures.len());
        verdict(
            "rb-001",
            "four scaled thin-shell falsifiers defeat the naive unit-bound marcher; the \
             certified path brackets and localizes every transverse first crossing, retains its \
             safe-step ratio audit on 120 grazing-biased rays, and leaves missing bounds \
             explicitly uncertified",
        );
    });
}

#[test]
fn rb_001a_over_relaxation_cannot_accept_a_far_shell_boundary() {
    with_cx(|cx| {
        let shell = thin_shell();
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 3.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let (plain, plain_audit) = sphere_trace(&shell, cx, &ray, 6.0, 1e-9, 1.0);
        let plain = plain.expect("plain trace closes a short bracket at the front boundary");
        assert_eq!(plain_audit.termination, TraceTermination::Hit);
        assert!(plain_audit.certified && plain_audit.certifies_hit());
        assert!((plain.t - 2.0).abs() <= 1e-9, "analytic first root is t=2");

        for (omega, t_max) in [(1.001, 6.0), (1.6, 2.5)] {
            let (relaxed, audit) = sphere_trace(&shell, cx, &ray, t_max, 1e-9, omega);
            let relaxed = relaxed.expect("relaxed trace must retreat then bracket the front root");
            assert!(
                (relaxed.t - 2.0).abs() <= 1e-9,
                "a wide speculative bracket must not select the far boundary: {}",
                relaxed.t
            );
            assert!(audit.certified && audit.certifies_hit());
            assert_eq!(audit.termination, TraceTermination::Hit);
            assert!(audit.fallbacks > 0, "thin-shell overshoot must retreat");
        }
        verdict(
            "rb-001a-over-relaxed-thin-shell",
            "wide sign-changing relaxed probes are never accepted directly; omega 1.001 and a \
             t_max-crossing omega 1.6 both retreat, then close a short bracket at the front root",
        );
    });
}

#[test]
fn rb_001b_trace_audit_states_fail_closed() {
    with_cx(|cx| {
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 3.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let bad_ray = Ray {
            origin: ray.origin,
            dir: Vec3::new(0.0, 0.0, 0.0),
        };
        let (_, invalid_input) = sphere_trace(&thin_shell(), cx, &bad_ray, 6.0, 1e-6, 1.0);
        assert_eq!(invalid_input.termination, TraceTermination::InvalidInput);
        assert!(!invalid_input.certified);

        let invalid_chart = SampleOverrideChart {
            inner: thin_shell(),
            distance: 1.0,
            lipschitz: Some(0.0),
            trace_claim: TraceStepClaim::LipschitzImplicit,
            rigorous: true,
        };
        let (_, invalid_sample) = sphere_trace(&invalid_chart, cx, &ray, 6.0, 1e-6, 1.0);
        assert_eq!(invalid_sample.termination, TraceTermination::InvalidSample);
        assert!(!invalid_sample.certified);

        let estimated_exact = SampleOverrideChart {
            inner: thin_shell(),
            distance: 1.0,
            lipschitz: Some(1.0),
            trace_claim: TraceStepClaim::ExactDistance,
            rigorous: false,
        };
        let (_, estimated_exact_audit) = sphere_trace(&estimated_exact, cx, &ray, 6.0, 1e-6, 1.0);
        assert_eq!(
            estimated_exact_audit.termination,
            TraceTermination::InvalidSample,
            "an analytical exact-distance claim cannot promote an estimated evaluation"
        );
        assert!(!estimated_exact_audit.certified);

        let estimated_lipschitz = SampleOverrideChart {
            inner: thin_shell(),
            distance: 1.0,
            lipschitz: Some(1.0),
            trace_claim: TraceStepClaim::LipschitzImplicit,
            rigorous: false,
        };
        let (_, estimated_lipschitz_audit) =
            sphere_trace(&estimated_lipschitz, cx, &ray, 6.0, 1e-6, 1.0);
        assert_eq!(
            estimated_lipschitz_audit.termination,
            TraceTermination::InvalidSample,
            "a Lipschitz theorem cannot promote an estimated rounded evaluation"
        );
        assert!(!estimated_lipschitz_audit.certified);

        let slow_chart = SampleOverrideChart {
            inner: thin_shell(),
            distance: 1e-12,
            lipschitz: Some(1.0),
            trace_claim: TraceStepClaim::LipschitzImplicit,
            rigorous: true,
        };
        let (_, limited) = sphere_trace(&slow_chart, cx, &ray, 1.0, 1e-15, 1.0);
        assert_eq!(limited.termination, TraceTermination::StepLimit);
        assert!(limited.certified);

        let (_, limited_relaxed) = sphere_trace(&slow_chart, cx, &ray, 1.0, 1e-15, 1.6);
        assert_eq!(limited_relaxed.termination, TraceTermination::StepLimit);
        assert!(limited_relaxed.certified);
        assert!(
            limited_relaxed.fallbacks > 0,
            "a pending relaxed endpoint retreats before the step-limit verdict"
        );

        let subnormal_chart = SampleOverrideChart {
            inner: thin_shell(),
            distance: f64::from_bits(1),
            lipschitz: Some(1e-300),
            trace_claim: TraceStepClaim::LipschitzImplicit,
            rigorous: true,
        };
        let (subnormal_hit, subnormal_audit) =
            sphere_trace(&subnormal_chart, cx, &ray, 1.0, 1e-30, 1.0);
        assert!(subnormal_hit.is_none());
        assert_eq!(
            subnormal_audit.termination,
            TraceTermination::InvalidSample,
            "a downward-rounded zero step must not become an early residual hit"
        );

        let declared_only = SampleOverrideChart {
            inner: thin_shell(),
            distance: 1.0,
            lipschitz: Some(1.0),
            trace_claim: TraceStepClaim::NoClaim,
            rigorous: false,
        };
        let (_, declared_only_audit) = sphere_trace(&declared_only, cx, &ray, 2.0, 1e-6, 1.0);
        assert!(
            !declared_only_audit.certified,
            "a sample Lipschitz number cannot upgrade a chart's typed no-claim"
        );

        let miss_ray = Ray {
            origin: Point3::new(0.0, 0.0, 3.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let (_, miss) = sphere_trace(&thin_shell(), cx, &miss_ray, 1.0, 1e-6, 1.0);
        assert_eq!(miss.termination, TraceTermination::Miss);
        assert!(miss.certified);
    });
}

#[test]
fn rb_001d_pre_cancelled_trace_stops_before_sampling() {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 9,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        gate.request();
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 3.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let (hit, audit) = sphere_trace(&thin_shell(), &cx, &ray, 6.0, 1e-6, 1.0);
        assert!(hit.is_none());
        assert_eq!(audit.steps, 0);
        assert_eq!(audit.termination, TraceTermination::Cancelled);
        assert!(!audit.certified);
        assert_eq!(
            trace_scene(&[Backend::Chart(&thin_shell())], &cx, &ray, 6.0, 1e-6),
            Err(SceneTraceError::Cancelled)
        );
    });
}

#[test]
fn rb_001e_production_rejects_uncertified_misses() {
    with_cx(|cx| {
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
        };
        assert_eq!(
            trace_scene(&[Backend::Chart(&ConstantNoClaim)], cx, &ray, 2.0, 1e-6),
            Err(SceneTraceError::UncertifiedTrace),
            "an uncertified preview miss is not proof of empty geometry"
        );
    });
}

#[test]
fn rb_001h_implicit_residual_is_not_a_geometric_hit() {
    with_cx(|cx| {
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
        };
        let (hit, audit) = sphere_trace(&LooseLipschitzChart, cx, &ray, 2.0, 1e-9, 1.0);
        assert!(
            hit.is_none(),
            "a loose valid L must not create a Hit one world unit from the zero set"
        );
        assert_eq!(audit.termination, TraceTermination::ResidualLimit);
        assert!(
            audit.certified,
            "the no-tunneling step theorem remains valid"
        );
        assert!(
            !audit.certifies_hit(),
            "certified marching does not promote residual evidence to proximity"
        );
        assert_eq!(
            trace_scene(&[Backend::Chart(&LooseLipschitzChart)], cx, &ray, 2.0, 1e-9),
            Err(SceneTraceError::BackendFailure(
                TraceTermination::ResidualLimit
            )),
            "production scene composition rejects residual-only terminal evidence"
        );
    });
}

#[test]
fn rb_001i_implicit_singleton_zero_is_a_geometric_hit() {
    with_cx(|cx| {
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
        };
        let (hit, audit) = sphere_trace(&SingletonZeroImplicitChart, cx, &ray, 2.0, 1e-12, 1.0);
        let hit = hit.expect("a rigorous implicit [0,0] is exact zero-set evidence");
        assert_eq!(hit.t.to_bits(), 0.0f64.to_bits());
        assert_eq!(hit.point, ray.origin);
        assert_eq!(audit.termination, TraceTermination::Hit);
        assert!(audit.certified && audit.certifies_hit());
    });
}

#[test]
fn rb_001j_short_bracket_clips_to_tmax_without_adopting_the_witness() {
    with_cx(|cx| {
        let eps = 1e-6;
        let ray = Ray {
            origin: Point3::new(-0.5 * eps, 0.0, 0.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
        };

        let (hit, audit) =
            sphere_trace(&SingletonZeroImplicitChart, cx, &ray, 0.75 * eps, eps, 1.0);
        let hit = hit.expect("the caller boundary supplies the opposite-sign witness");
        assert!(
            hit.t <= 0.75 * eps,
            "the witness is clipped to caller t_max"
        );
        assert!(
            hit.point.x.abs() <= eps,
            "the evaluated midpoint stays within eps of the exact plane"
        );
        assert_eq!(audit.steps, 0, "the witness is evidence, not march state");
        assert!(audit.certified && audit.certifies_hit());

        let (blocked_hit, blocked_audit) =
            sphere_trace(&SingletonZeroImplicitChart, cx, &ray, 0.25 * eps, eps, 1.0);
        assert!(
            blocked_hit.is_none(),
            "a root beyond t_max must not leak in"
        );
        assert_eq!(blocked_audit.termination, TraceTermination::ResidualLimit);
        assert_eq!(
            blocked_audit.steps, 0,
            "a same-sign witness never becomes the next march sample"
        );
        assert!(blocked_audit.certified && !blocked_audit.certifies_hit());
    });
}

#[test]
fn rb_001k_indeterminate_same_sign_and_tangent_witnesses_fail_closed() {
    with_cx(|cx| {
        let eps = 1e-6;
        let x_ray = Ray {
            origin: Point3::new(-0.5 * eps, 0.0, 0.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
        };
        let (straddling_hit, straddling_audit) =
            sphere_trace(&StraddlingWitnessChart, cx, &x_ray, 1.0, eps, 1.0);
        assert!(straddling_hit.is_none());
        assert_eq!(
            straddling_audit.termination,
            TraceTermination::ResidualLimit,
            "an enclosure containing zero does not prove a root"
        );
        assert_eq!(straddling_audit.steps, 0, "the witness is not adopted");

        let shell = thin_shell_with_thickness(0.5 * eps);
        let shell_ray = Ray {
            origin: Point3::new(0.0, 0.0, 1.0 + 0.5 * eps),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let (shell_hit, shell_audit) = sphere_trace(&shell, cx, &shell_ray, 0.01, eps, 1.0);
        assert!(
            shell_hit.is_none(),
            "a probe spanning both sides of an ultra-thin shell has no sign witness"
        );
        assert_eq!(shell_audit.termination, TraceTermination::ResidualLimit);
        assert_eq!(shell_audit.steps, 0, "the same-sign probe is never a step");

        let tangent = ImplicitSphereChart {
            inner: SphereChart {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            },
        };
        let tangent_ray = Ray {
            origin: Point3::new(-0.5 * eps, 1.0, 0.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
        };
        let (tangent_hit, tangent_audit) = sphere_trace(&tangent, cx, &tangent_ray, 1.0, eps, 1.0);
        assert!(
            tangent_hit.is_none(),
            "a tangent needs proximity or first-root evidence beyond endpoint signs"
        );
        assert_eq!(tangent_audit.termination, TraceTermination::ResidualLimit);
        assert_eq!(tangent_audit.steps, 0, "the tangent probe is never adopted");
        assert!(straddling_audit.certified && shell_audit.certified && tangent_audit.certified);
    });
}

#[test]
fn rb_001l_cancellation_during_short_witness_wins() {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 9,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let chart = CancelOnWitnessChart {
            gate: &gate,
            evaluations: AtomicUsize::new(0),
        };
        let ray = Ray {
            origin: Point3::new(-5e-7, 0.0, 0.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
        };
        let (hit, audit) = sphere_trace(&chart, &cx, &ray, 1.0, 1e-6, 1.0);
        assert!(hit.is_none());
        assert_eq!(audit.termination, TraceTermination::Cancelled);
        assert!(
            !audit.certified,
            "cancellation cannot publish a completed proof"
        );
        assert_eq!(
            chart.evaluations.load(Ordering::SeqCst),
            2,
            "only the current sample and speculative witness are evaluated"
        );
    });
}

#[test]
fn rb_001f_cancellation_requested_by_terminal_eval_wins() {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 9,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let chart = CancelOnEvalChart { gate: &gate };
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
        };
        assert_eq!(
            trace_scene(&[Backend::Chart(&chart)], &cx, &ray, 2.0, 1e-6),
            Err(SceneTraceError::Cancelled),
            "cancellation raised inside the final evaluation must beat its hit"
        );
    });
}

#[test]
fn rb_001c_scale_direction_and_touching_spheres_are_conservative() {
    with_cx(|cx| {
        let axial_ray = Ray {
            origin: Point3::new(0.0, 0.0, 3.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let tiny_scale = ScaledChart {
            inner: thin_shell(),
            scale: 1e-9,
            publish_bound: true,
            trace_claim: TraceStepClaim::LipschitzImplicit,
        };
        let (scaled_hit, scaled_audit) = sphere_trace(&tiny_scale, cx, &axial_ray, 6.0, 1e-6, 1.0);
        let scaled_hit = scaled_hit.expect("the short bracket is invariant to field scaling");
        assert!((scaled_hit.t - 2.0).abs() <= 1e-6);
        assert_eq!(scaled_audit.termination, TraceTermination::Hit);
        assert!(scaled_audit.certified && scaled_audit.certifies_hit());

        let plane = ExactPlaneChart {
            boundary: 1.0,
            lipschitz: 1.0,
        };
        let x_ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
        };
        let (plane_hit, plane_audit) = sphere_trace(&plane, cx, &x_ray, 3.0, 1e-12, 1.6);
        let plane_hit = plane_hit.expect("touching speculative spheres retreat to the plane");
        assert!((plane_hit.point.x - 1.0).abs() <= 1e-12);
        assert!(plane_audit.certified && plane_audit.fallbacks > 0);

        let loose_exact_plane = ExactPlaneChart {
            boundary: 1.0,
            lipschitz: 2.0,
        };
        let (loose_hit, loose_audit) = sphere_trace(&loose_exact_plane, cx, &x_ray, 3.0, 0.1, 1.0);
        let loose_hit = loose_hit.expect("a conservative exact-distance bound converges");
        assert!(loose_hit.t > 0.5, "loose L must not turn t=0 into a hit");
        assert!((1.0 - loose_hit.point.x).abs() <= 0.1);
        assert!(loose_audit.certified);

        let invalid_exact_plane = ExactPlaneChart {
            boundary: 1.0,
            lipschitz: 0.5,
        };
        let (invalid_hit, invalid_audit) =
            sphere_trace(&invalid_exact_plane, cx, &x_ray, 3.0, 1e-6, 1.0);
        assert!(invalid_hit.is_none());
        assert_eq!(invalid_audit.termination, TraceTermination::InvalidSample);

        let boundary_ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(3.0, 0.0, 0.0),
        };
        let beyond_limit = ExactPlaneChart {
            boundary: (0.1_f64 * 3.0).next_up(),
            lipschitz: 1.0,
        };
        let (boundary_hit, boundary_audit) =
            sphere_trace(&beyond_limit, cx, &boundary_ray, 0.1, 1e-18, 1.0);
        assert!(
            boundary_hit.is_none(),
            "outward t_max rounding must not leak a hit"
        );
        assert_eq!(boundary_audit.termination, TraceTermination::Miss);
        assert!(boundary_audit.certified);
        let (relaxed_boundary_hit, relaxed_boundary_audit) =
            sphere_trace(&beyond_limit, cx, &boundary_ray, 0.1, 1e-18, 1.6);
        assert!(
            relaxed_boundary_hit.is_none(),
            "an out-of-range speculative endpoint must retreat to the bounded miss"
        );
        assert_eq!(relaxed_boundary_audit.termination, TraceTermination::Miss);
        assert!(relaxed_boundary_audit.certified);
        assert_eq!(
            trace_scene(
                &[Backend::Chart(&beyond_limit)],
                cx,
                &boundary_ray,
                0.1,
                1e-18,
            ),
            Ok(None),
            "production composition must accept the certified bounded miss"
        );

        let exact_limit = ExactPlaneChart {
            boundary: boundary_ray.at(0.1).x,
            lipschitz: 1.0,
        };
        let (limit_hit, limit_audit) =
            sphere_trace(&exact_limit, cx, &boundary_ray, 0.1, 1e-18, 1.0);
        let limit_hit = limit_hit.expect("an exact hit at caller t_max must not become a miss");
        assert_eq!(limit_hit.t, 0.1);
        assert_eq!(limit_hit.point, boundary_ray.at(0.1));
        assert_eq!(limit_audit.termination, TraceTermination::Hit);
        assert!(limit_audit.certified);
        let (relaxed_limit_hit, relaxed_limit_audit) =
            sphere_trace(&exact_limit, cx, &boundary_ray, 0.1, 1e-18, 1.6);
        let relaxed_limit_hit =
            relaxed_limit_hit.expect("over-relaxation must retain the exact caller endpoint");
        assert_eq!(relaxed_limit_hit.t, limit_hit.t);
        assert_eq!(relaxed_limit_hit.point, limit_hit.point);
        assert_eq!(relaxed_limit_audit.termination, TraceTermination::Hit);
        assert!(relaxed_limit_audit.certified);
        let scene_limit_hit = trace_scene(
            &[Backend::Chart(&exact_limit)],
            cx,
            &boundary_ray,
            0.1,
            1e-18,
        )
        .expect("production composition accepts a certified endpoint hit")
        .expect("the certified endpoint is present");
        assert_eq!(scene_limit_hit.0, 0);
        assert_eq!(scene_limit_hit.1.t, limit_hit.t);
        assert_eq!(scene_limit_hit.1.point, limit_hit.point);

        let equality_ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(0.1, 0.0, 0.0),
        };
        let equality_t_max = 0.3;
        let outward_working_x = (equality_t_max * equality_ray.dir.x).next_up();
        assert_eq!(outward_working_x / equality_ray.dir.x, equality_t_max);
        assert_ne!(outward_working_x, equality_ray.at(equality_t_max).x);
        let equality_gap = ExactPlaneChart {
            boundary: outward_working_x,
            lipschitz: 1.0,
        };
        let (equality_hit, equality_audit) =
            sphere_trace(&equality_gap, cx, &equality_ray, equality_t_max, 1e-20, 1.0);
        assert!(
            equality_hit.is_none(),
            "equal mapped parameters must still classify the caller endpoint"
        );
        assert_eq!(equality_audit.termination, TraceTermination::Miss);
        assert!(equality_audit.certified);

        let slab = ExactSlabChart {
            lo: 1.0,
            hi: 1.0 + 2e-13,
        };
        let near_unit_ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(1.0 + 5e-13, 0.0, 0.0),
        };
        let (slab_hit, slab_audit) = sphere_trace(&slab, cx, &near_unit_ray, 2.0, 1e-15, 1.0);
        let slab_hit = slab_hit.expect("parameter steps account for non-unit ray speed");
        assert_eq!(slab_hit.point, near_unit_ray.at(slab_hit.t));
        assert!((slab_hit.point.x - 1.0).abs() <= 1e-14);
        assert!(slab_audit.certified);

        let horn = TorusChart {
            center: Point3::new(0.0, 0.0, 0.0),
            major: 1.0,
            minor: 1.0,
        };
        assert_eq!(horn.trace_step_claim(), TraceStepClaim::LipschitzImplicit);
        let horn_point = Point3::new(0.0, 0.0, 1e-3);
        let horn_sample = horn.eval(horn_point, cx);
        assert_eq!(horn_sample.error.kind, NumericalKind::Estimate);
        assert_eq!(
            horn.trace_value_enclosure(horn_point, &horn_sample, cx)
                .kind,
            NumericalKind::Enclosure
        );
        let mut offset_builder = FrepBuilder::new();
        let base = offset_builder
            .sphere(Point3::new(0.0, 0.0, 0.0), 1.0)
            .expect("sphere");
        let offset = offset_builder.offset(base, 0.25).expect("offset");
        let offset_chart = offset_builder.finish(offset).expect("offset chart");
        assert_eq!(
            offset_chart.trace_step_claim(),
            TraceStepClaim::LipschitzImplicit,
            "offset exactness needs a reach certificate"
        );
    });
}

#[test]
fn rb_001d_rounded_exact_distance_overstatement_fails_closed() {
    with_cx(|cx| {
        let sphere = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let ray = Ray {
            origin: Point3::new(
                0.038_546_885_717_366_49,
                -0.607_415_449_300_028_1,
                -0.793_448_734_204_907_9,
            ),
            dir: Vec3::new(
                -0.038_546_880_238_732_83,
                0.607_415_362_968_625_2,
                0.793_448_621_432_764_2,
            ),
        };
        let first = sphere.eval(ray.origin, cx);
        let rounded_endpoint = ray.at(first.signed_distance);
        assert!(
            sphere
                .eval(rounded_endpoint, cx)
                .signed_distance
                .is_sign_negative(),
            "the nearest-rounded residual is the seeded unsafe step"
        );

        let (hit, audit) = sphere_trace(&sphere, cx, &ray, 1.0, 1e-18, 1.0);
        assert!(
            hit.is_none(),
            "an interval too wide to certify the first root must not return a far-side hit"
        );
        assert_eq!(audit.termination, TraceTermination::InvalidSample);
        assert!(
            !audit.certified,
            "the unresolved rounding band fails closed"
        );
        verdict(
            "rb-001d-rounded-distance",
            "a cancelled binary64 sphere residual that flips sign is enclosed; sub-band hit \
             tolerance fails closed instead of tunneling to the far boundary",
        );
    });
}

#[test]
fn rb_001g_rotated_frep_retains_certified_trace_progress() {
    with_cx(|cx| {
        let mut builder = FrepBuilder::new();
        let sphere = builder
            .sphere(Point3::new(1.0, 0.0, 0.0), 0.5)
            .expect("sphere");
        let rotated = builder
            .rotate(
                sphere,
                Vec3::new(0.0, 0.0, 1.0),
                core::f64::consts::FRAC_PI_2,
            )
            .expect("rotation");
        let rotated = builder.finish(rotated).expect("rotated F-rep");
        assert_eq!(
            rotated.trace_step_claim(),
            TraceStepClaim::LipschitzImplicit,
            "nontrivial rounded Rodrigues maps do not claim exact distance"
        );

        let ray = Ray {
            origin: Point3::new(0.0, 1.0, 3.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let first = rotated.eval(ray.origin, cx);
        let first_trace = rotated.trace_value_enclosure(ray.origin, &first, cx);
        assert_eq!(first_trace.kind, NumericalKind::Enclosure);
        assert!(
            first_trace.lo > 0.0,
            "ordinary-angle trig bounds must preserve a useful outside sign"
        );
        let (hit, audit) = sphere_trace(&rotated, cx, &ray, 6.0, 1e-9, 1.0);
        assert!(
            hit.is_none(),
            "the conservative Rodrigues L leaves the first residual outside a 2eps bracket"
        );
        assert_eq!(audit.termination, TraceTermination::ResidualLimit);
        assert!(audit.certified && !audit.certifies_hit());
        assert!(audit.steps > 0, "tight trig bounds retain trace progress");
        assert_eq!(
            trace_scene(&[Backend::Chart(&rotated)], cx, &ray, 6.0, 1e-9),
            Err(SceneTraceError::BackendFailure(
                TraceTermination::ResidualLimit
            )),
            "production refuses a transverse root that the bounded witness did not reach"
        );
        verdict(
            "rb-001g-rotated-frep",
            "a nontrivial rotated sphere retains certified progress under fs-ivl trig \
             enclosures, then remains residual-only when the conservative L leaves the root \
             beyond the bounded sign witness",
        );
    });
}

#[test]
fn rb_002_over_relaxation_stays_certified() {
    with_cx(|cx| {
        let shell = thin_shell();
        let mut state = 0xface_u64;
        let mut plain_steps = 0u64;
        let mut relaxed_steps = 0u64;
        for k in 0..60 {
            // GRAZING rays: passing at 1.02-1.15x the radius, parallel
            // to an axis — the many-small-steps regime over-relaxation
            // exists for (center-aimed rays converge in one jump and
            // have nothing to accelerate).
            let offset = if k % 2 == 0 {
                0.98 + 0.019 * lcg(&mut state)
            } else {
                1.02 + 0.13 * lcg(&mut state)
            };
            let phi = lcg(&mut state) * std::f64::consts::TAU;
            let origin = Point3::new(-4.0, offset * phi.cos(), offset * phi.sin());
            let ray = Ray {
                origin,
                dir: Vec3::new(1.0, 0.0, 0.0),
            };
            let (h1, a1) = sphere_trace(&shell, cx, &ray, 9.0, 1e-6, 1.0);
            let (h2, a2) = sphere_trace(&shell, cx, &ray, 9.0, 1e-6, 1.6);
            assert_eq!(h1.is_some(), h2.is_some(), "same hit decision (ray {k})");
            match (h1, h2) {
                (Some(plain), Some(relaxed)) => {
                    assert_eq!(a1.termination, TraceTermination::Hit);
                    assert_eq!(a2.termination, TraceTermination::Hit);
                    assert!(
                        (plain.t - relaxed.t).abs() <= 2e-6,
                        "same tolerance bracket: {} vs {}",
                        plain.t,
                        relaxed.t
                    );
                }
                (None, None) => {
                    assert_eq!(a1.termination, a2.termination, "same typed stop (ray {k})");
                    assert!(matches!(
                        a1.termination,
                        TraceTermination::ResidualLimit | TraceTermination::Miss
                    ));
                }
                _ => unreachable!("presence equality checked above"),
            }
            assert!(a1.certified && a2.certified);
            plain_steps += u64::from(a1.steps);
            relaxed_steps += u64::from(a2.steps);
        }
        println!(
            "{{\"metric\":\"over-relaxation\",\"plain_steps\":{plain_steps},\
             \"relaxed_steps\":{relaxed_steps}}}"
        );
        assert!(
            relaxed_steps < plain_steps,
            "relaxation saves steps: {relaxed_steps} vs {plain_steps}"
        );
        verdict(
            "rb-002",
            "omega=1.6 preserves bracketed hit and typed miss/residual decisions with fewer \
             steps (certified fallback preserves the no-tunnel invariant)",
        );
    });
}

/// The exact NURBS unit sphere (rational quadratic revolution).
fn sphere_nurbs() -> fs_rep_nurbs::NurbsSurface<f64> {
    let s2 = std::f64::consts::FRAC_1_SQRT_2;
    let circle = [
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
        [-1.0, 1.0],
        [-1.0, 0.0],
        [-1.0, -1.0],
        [0.0, -1.0],
        [1.0, -1.0],
        [1.0, 0.0],
    ];
    let cw = |i: usize| if i.is_multiple_of(2) { 1.0 } else { s2 };
    let profile: [([f64; 2], f64); 5] = [
        ([0.0, -1.0], 1.0),
        ([1.0, -1.0], s2),
        ([1.0, 0.0], 1.0),
        ([1.0, 1.0], s2),
        ([0.0, 1.0], 1.0),
    ];
    let mut points: Vec<Vec<[f64; 3]>> = Vec::new();
    let mut weights: Vec<Vec<f64>> = Vec::new();
    for (i, c) in circle.iter().enumerate() {
        let mut prow = Vec::new();
        let mut wrow = Vec::new();
        for &([radius, z], wv) in &profile {
            prow.push([radius * c[0], radius * c[1], z]);
            wrow.push(cw(i) * wv);
        }
        points.push(prow);
        weights.push(wrow);
    }
    let ku = fs_rep_nurbs::KnotVector::new(
        vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        2,
    )
    .expect("ku");
    let kv =
        fs_rep_nurbs::KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0], 2).expect("kv");
    fs_rep_nurbs::NurbsSurface::new(ku, kv, &points, &weights).expect("sphere")
}

#[test]
fn rb_003_nurbs_newton_matches_analytic() {
    let sphere = sphere_nurbs();
    let mut state = 0x5eed_u64;
    for _ in 0..40 {
        let phi = lcg(&mut state) * std::f64::consts::TAU;
        let z = 1.6 * lcg(&mut state) - 0.8;
        let origin = Point3::new(4.0 * phi.cos(), 4.0 * phi.sin(), z);
        let dir = unit(Point3::new(0.0, 0.0, 0.0).delta_from(origin));
        let ray = Ray { origin, dir };
        // Analytic sphere intersection: |o + t d| = 1.
        let oc = [origin.x, origin.y, origin.z];
        let b = oc[0] * dir.x + oc[1] * dir.y + oc[2] * dir.z;
        let c = oc[0] * oc[0] + oc[1] * oc[1] + oc[2] * oc[2] - 1.0;
        let disc = b * b - c;
        assert!(disc > 0.0, "aimed at the center: always hits");
        let t_ref = -b - disc.sqrt();
        let hit = ray_intersect_nurbs(&sphere, &ray, 8, 1e-9)
            .expect("valid NURBS ray request")
            .expect("Newton converges");
        assert!(
            (hit.t - t_ref).abs() < 1e-6,
            "Newton matches analytic: {} vs {t_ref}",
            hit.t
        );
        // The normal points outward at the hit point.
        let n = hit.normal.expect("normal");
        let outward = unit(hit.point.delta_from(Point3::new(0.0, 0.0, 0.0)));
        assert!(
            n.dot(outward).abs() > 0.999,
            "normal aligned with the radial direction"
        );
        for scale in [1e-12, 1e12, 1e300] {
            let scaled_ray = Ray {
                origin,
                dir: dir.scale(scale),
            };
            let scaled_hit = ray_intersect_nurbs(&sphere, &scaled_ray, 8, 1e-9)
                .expect("valid scaled NURBS ray request")
                .expect("direction scaling preserves the NURBS hit");
            assert_eq!(
                scaled_hit.point,
                scaled_ray.at(scaled_hit.t),
                "NURBS Hit.point must use the caller's parameterization"
            );
            assert!(
                (scaled_hit.t * scale - t_ref).abs() < 1e-6,
                "scaled Newton parameter maps to the same point"
            );
        }
    }
    verdict(
        "rb-003",
        "40 rays: Bezier-seeded Newton matches the analytic sphere intersection to 1e-6 \
         with radial normals",
    );
}

#[test]
fn rb_003a_nurbs_newton_is_invariant_to_parameter_domain_scale() {
    fn plane(parameter_hi: f64) -> fs_rep_nurbs::NurbsSurface<f64> {
        let knots_u = fs_rep_nurbs::KnotVector::new(vec![0.0, 0.0, parameter_hi, parameter_hi], 1)
            .expect("valid u knots");
        let knots_v = fs_rep_nurbs::KnotVector::new(vec![0.0, 0.0, parameter_hi, parameter_hi], 1)
            .expect("valid v knots");
        let points = vec![
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let weights = vec![vec![1.0; 2]; 2];
        fs_rep_nurbs::NurbsSurface::new(knots_u, knots_v, &points, &weights).expect("valid plane")
    }

    let ray = Ray {
        origin: Point3::new(0.37, 0.41, 2.0),
        dir: Vec3::new(0.0, 0.0, -1.0),
    };
    let reference = ray_intersect_nurbs(&plane(1.0), &ray, 3, 1e-12)
        .expect("unit-domain request is valid")
        .expect("unit-domain plane intersects");
    let reparameterized = ray_intersect_nurbs(&plane(1e15), &ray, 3, 1e-12)
        .expect("large-domain request is valid")
        .expect("surface parameterization cannot turn a hit into a miss");
    assert!((reference.t - 2.0).abs() <= 1e-12);
    assert!((reparameterized.t - reference.t).abs() <= 1e-12);
    assert!((reparameterized.point.x - 0.37).abs() <= 1e-12);
    assert!((reparameterized.point.y - 0.41).abs() <= 1e-12);
    verdict(
        "rb-003a",
        "G3: a 1e15 knot-domain reparameterization preserves the same finite Newton hit",
    );
}

#[test]
fn rb_003b_nurbs_ray_admission_is_bounded_and_fail_closed() {
    let sphere = sphere_nurbs();
    let ray = Ray {
        origin: Point3::new(3.0, 0.0, 0.0),
        dir: Vec3::new(-1.0, 0.0, 0.0),
    };
    for eps in [f64::NAN, f64::INFINITY, 0.0, -1.0] {
        assert!(matches!(
            ray_intersect_nurbs(&sphere, &ray, 8, eps),
            Err(NurbsRayError::InvalidInput { .. })
        ));
    }
    assert!(matches!(
        ray_intersect_nurbs(&sphere, &ray, 0, 1e-9),
        Err(NurbsRayError::InvalidInput { .. })
    ));
    assert!(matches!(
        ray_intersect_nurbs(&sphere, &ray, usize::MAX, 1e-9),
        Err(NurbsRayError::ResourceLimit { .. })
    ));
    assert!(matches!(
        ray_intersect_nurbs(&sphere, &ray, 256, 1e-9),
        Err(NurbsRayError::WorkLimit { .. })
    ));

    let overflow_knots =
        fs_rep_nurbs::KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("linear knots");
    let overflow_partials = fs_rep_nurbs::NurbsSurface::new(
        overflow_knots.clone(),
        overflow_knots,
        &[
            vec![[-f64::MAX, 0.0, 0.0], [-f64::MAX, 1.0, 0.0]],
            vec![[f64::MAX, 0.0, 0.0], [f64::MAX, 1.0, 0.0]],
        ],
        &[vec![1.0; 2], vec![1.0; 2]],
    )
    .expect("finite surface whose u derivative overflows");
    assert!(
        overflow_partials.eval(0.5, 0.5).is_ok(),
        "the seed value remains finite"
    );
    assert!(matches!(
        overflow_partials.partials(0.5, 0.5),
        Err(fs_rep_nurbs::NurbsError::Domain { .. })
    ));
    let overflow_ray = Ray {
        origin: Point3::new(0.0, 0.5, 1.0),
        dir: Vec3::new(0.0, 0.0, -1.0),
    };
    assert!(matches!(
        ray_intersect_nurbs(&overflow_partials, &overflow_ray, 1, 1e-9),
        Err(NurbsRayError::InvalidSurface(
            fs_rep_nurbs::NurbsError::Domain { .. }
        ))
    ));

    assert!(
        fs_rep_nurbs::NurbsSurface::from_homogeneous(
            sphere.knots_u().clone(),
            sphere.knots_v().clone(),
            Vec::new(),
        )
        .is_err(),
        "sealed surfaces cannot be constructed with an empty control net"
    );
    let geometric_miss = Ray {
        origin: Point3::new(3.0, 3.0, 3.0),
        dir: Vec3::new(1.0, 0.0, 0.0),
    };
    assert!(matches!(
        ray_intersect_nurbs(&sphere_nurbs(), &geometric_miss, 8, 1e-9),
        Err(NurbsRayError::IterationLimit { .. })
    ));
    verdict(
        "rb-003b",
        "invalid requests refuse before unbounded work, and heuristic nonconvergence cannot impersonate a certified miss",
    );
}

/// An octahedron-subdivision icosphere-ish mesh of the unit sphere.
fn sphere_mesh(subdiv: usize) -> TriMesh {
    // Start from an octahedron, subdivide, project to the sphere.
    let mut verts: Vec<[f64; 3]> = vec![
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ];
    let mut tris: Vec<[u32; 3]> = vec![
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5],
    ];
    for _ in 0..subdiv {
        let mut next = Vec::with_capacity(tris.len() * 4);
        for t in &tris {
            let mid = |a: u32, b: u32, verts: &mut Vec<[f64; 3]>| -> u32 {
                let (pa, pb) = (verts[a as usize], verts[b as usize]);
                let mut m = [
                    f64::midpoint(pa[0], pb[0]),
                    f64::midpoint(pa[1], pb[1]),
                    f64::midpoint(pa[2], pb[2]),
                ];
                let n = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt();
                for v in &mut m {
                    *v /= n;
                }
                verts.push(m);
                (verts.len() - 1) as u32
            };
            let ab = mid(t[0], t[1], &mut verts);
            let bc = mid(t[1], t[2], &mut verts);
            let ca = mid(t[2], t[0], &mut verts);
            next.extend_from_slice(&[[t[0], ab, ca], [ab, t[1], bc], [ca, bc, t[2]], [ab, bc, ca]]);
        }
        tris = next;
    }
    TriMesh::new(verts, tris)
}

#[test]
fn rb_004a_bvh_build_is_deterministic_under_concurrent_construction() {
    let source = sphere_mesh(3);
    let vertices = source.vertices.clone();
    let triangles = source.triangles.clone();
    let ray = Ray {
        origin: Point3::new(3.0, 0.2, -0.1),
        dir: unit(Point3::new(0.0, 0.0, 0.0).delta_from(Point3::new(3.0, 0.2, -0.1))),
    };
    let reference_mesh = TriMesh::new(vertices.clone(), triangles.clone());
    let reference_hit = reference_mesh.intersect(&ray).expect("reference mesh hit");
    let reference = (
        reference_mesh.bvh_fingerprint(),
        reference_hit.t.to_bits(),
        reference_hit.steps,
    );

    for workers in [1usize, 2, 4, 8] {
        let outcomes = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(workers);
            for _ in 0..workers {
                let worker_vertices = vertices.clone();
                let worker_triangles = triangles.clone();
                handles.push(scope.spawn(move || {
                    let mesh = TriMesh::new(worker_vertices, worker_triangles);
                    let hit = mesh.intersect(&ray).expect("worker mesh hit");
                    (mesh.bvh_fingerprint(), hit.t.to_bits(), hit.steps)
                }));
            }
            handles
                .into_iter()
                .map(|handle| handle.join().expect("BVH worker did not panic"))
                .collect::<Vec<_>>()
        });
        assert!(
            outcomes.iter().all(|&outcome| outcome == reference),
            "same input must produce the same BVH and hit bits with {workers} concurrent builders"
        );
    }
    println!(
        "{{\"metric\":\"bvh-determinism\",\"fingerprint\":\"{:#018x}\",\"hit_t_bits\":\"{:#018x}\",\"visits\":{},\"concurrent_builders\":[1,2,4,8]}}",
        reference.0, reference.1, reference.2
    );
    verdict(
        "rb-004a-bvh-determinism",
        "same ordered mesh input produces one BVH fingerprint and bit-identical hit receipt \
         across 1, 2, 4, and 8 concurrent builders",
    );
}

#[test]
fn rb_004b_bvh_grazing_pruning_matches_bruteforce() {
    let mesh = TriMesh::new(
        vec![
            [1.0e12, 0.0, 0.0],
            [1.0e12 + 1.0, 0.0, 0.0],
            [1.0e12 + 1.0, 1.0, 0.0],
            [1.0e12, 1.0, 0.0],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    );
    for y in [0.0, f64::from_bits(1), 0.5, 1.0 - f64::EPSILON, 1.0] {
        let ray = Ray {
            origin: Point3::new(1.0e12, y, 1.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let bvh = mesh.intersect(&ray);
        let brute = mesh.intersect_bruteforce(&ray);
        assert_eq!(
            bvh.map(|hit| hit.t.to_bits()),
            brute.map(|hit| hit.t.to_bits()),
            "BVH pruning changed a grazing result at y={y:.17e}"
        );
    }
}

#[test]
fn rb_004_mixed_scene_consistency_and_frame_invariance() {
    with_cx(|cx| {
        // The SAME unit sphere held three ways.
        let mut b = FrepBuilder::new();
        let node = b.sphere(Point3::new(0.0, 0.0, 0.0), 1.0).expect("sphere");
        let frep = b.finish(node).expect("frep");
        let nurbs = sphere_nurbs();
        let mesh = sphere_mesh(4);
        let mut state = 0xd1ce_u64;
        let mut worst_spread = 0.0f64;
        for _ in 0..30 {
            let phi = lcg(&mut state) * std::f64::consts::TAU;
            let z = 1.2 * lcg(&mut state) - 0.6;
            let origin = Point3::new(3.5 * phi.cos(), 3.5 * phi.sin(), z);
            let ray = Ray {
                origin,
                dir: unit(Point3::new(0.0, 0.0, 0.0).delta_from(origin)),
            };
            let (sdf_hit, _) = sphere_trace(&frep, cx, &ray, 8.0, 1e-7, 1.0);
            let nurbs_hit =
                ray_intersect_nurbs(&nurbs, &ray, 8, 1e-9).expect("valid NURBS ray request");
            let mesh_hit = mesh.intersect(&ray);
            let (a, b_, c) = (
                sdf_hit.expect("sdf hits").t,
                nurbs_hit.expect("nurbs hits").t,
                mesh_hit.expect("mesh hits").t,
            );
            worst_spread = worst_spread.max((a - b_).abs()).max((a - c).abs());
            for scale in [1e-12, 1e12, 1e300] {
                let scaled_ray = Ray {
                    origin,
                    dir: ray.dir.scale(scale),
                };
                let scaled_hit = mesh
                    .intersect(&scaled_ray)
                    .expect("direction scaling preserves the mesh hit");
                assert_eq!(
                    scaled_hit.point,
                    scaled_ray.at(scaled_hit.t),
                    "mesh Hit.point must use the caller's parameterization"
                );
                assert!(
                    (scaled_hit.t * scale - c).abs() <= 1e-12 * c.abs().max(1.0),
                    "mesh hit is invariant after mapping back to world distance"
                );
            }
        }
        // Mesh is a level-4 subdivision: geometric error ~ 1e-3.
        assert!(
            worst_spread < 5e-3,
            "the same shape across three backends: spread {worst_spread}"
        );
        let bounded_ray = Ray {
            origin: Point3::new(4.0, 0.0, 0.0),
            dir: Vec3::new(-1.0, 0.0, 0.0),
        };
        assert!(
            trace_scene(&[Backend::Nurbs(&nurbs)], cx, &bounded_ray, 2.0, 1e-9)
                .expect("bounded NURBS trace")
                .is_none(),
            "NURBS hits beyond t_max must not leak into the scene"
        );
        assert!(
            trace_scene(&[Backend::Mesh(&mesh)], cx, &bounded_ray, 2.0, 1e-9)
                .expect("bounded mesh trace")
                .is_none(),
            "mesh hits beyond t_max must not leak into the scene"
        );
        // Mixed scene: three unit spheres at distinct centers, one per
        // backend — the closest-hit wins and identifies its instance.
        let mut b2 = FrepBuilder::new();
        let sn = b2.sphere(Point3::new(-3.0, 0.0, 0.0), 1.0).expect("s");
        let frep2 = b2.finish(sn).expect("frep");
        let backends = [
            Backend::Chart(&frep2),
            Backend::Nurbs(&nurbs),
            Backend::Mesh(&mesh),
        ];
        let ray = Ray {
            origin: Point3::new(-8.0, 0.0, 0.01),
            dir: Vec3::new(1.0, 0.0, 0.0),
        };
        let (idx, hit) = trace_scene(&backends, cx, &ray, 20.0, 1e-7)
            .expect("mixed trace completes")
            .expect("hits");
        assert_eq!(idx, 0, "the F-rep sphere at x=-3 is closest");
        assert!((hit.t - 4.0).abs() < 1e-3, "t ~ 4: {}", hit.t);
        // G3 frame invariance: translate everything by the same offset.
        let ray_shifted = Ray {
            origin: Point3::new(-8.0 + 5.0, 0.0, 0.01),
            dir: Vec3::new(1.0, 0.0, 0.0),
        };
        let mut b3 = FrepBuilder::new();
        let sn3 = b3.sphere(Point3::new(2.0, 0.0, 0.0), 1.0).expect("s");
        let frep3 = b3.finish(sn3).expect("frep");
        let (t1, _) = sphere_trace(&frep2, cx, &ray, 20.0, 1e-7, 1.0)
            .0
            .map(|h| (h.t, h.steps))
            .expect("hit");
        let (t2, _) = sphere_trace(&frep3, cx, &ray_shifted, 20.0, 1e-7, 1.0)
            .0
            .map(|h| (h.t, h.steps))
            .expect("hit");
        assert!(
            (t1 - t2).abs() < 1e-9,
            "translation invariance: {t1} vs {t2}"
        );
        verdict(
            "rb-004",
            "one sphere, three backends, <5e-3 spread; mixed scene picks the closest \
             instance; translated scene reproduces t to 1e-9",
        );
    });
}

#[test]
fn rb_005_ray_rate_ledger() {
    with_cx(|cx| {
        // MEASURED, honestly labeled throughput on primary rays. The plan's
        // Mray/s TARGETS are release-build perf-CI
        // gates (fz2.4) — this ledgers the measurement discipline, not
        // the target (AGENTS: no "fast" without a benchmark + machine).
        let mut b = FrepBuilder::new();
        let node = b.sphere(Point3::new(0.0, 0.0, 0.0), 1.0).expect("sphere");
        let frep = b.finish(node).expect("frep");
        let n = 64usize;
        let start = std::time::Instant::now();
        let mut hits = 0usize;
        for py in 0..n {
            for px in 0..n {
                #[allow(clippy::cast_precision_loss)]
                let (u, v) = (
                    (px as f64 + 0.5) / n as f64 * 2.0 - 1.0,
                    (py as f64 + 0.5) / n as f64 * 2.0 - 1.0,
                );
                let ray = Ray {
                    origin: Point3::new(1.6 * u, 1.6 * v, 3.0),
                    dir: Vec3::new(0.0, 0.0, -1.0),
                };
                if sphere_trace(&frep, cx, &ray, 8.0, 1e-6, 1.6).0.is_some() {
                    hits += 1;
                }
            }
        }
        let secs = start.elapsed().as_secs_f64();
        #[allow(clippy::cast_precision_loss)]
        let mray_s = (n * n) as f64 / secs / 1e6;
        let profile = if cfg!(debug_assertions) {
            "debug"
        } else {
            "release"
        };
        println!(
            "{{\"metric\":\"ray-rate\",\"build\":\"{profile}\",\"os\":\"{}\",\"arch\":\"{}\",\
             \"rays\":{},\"hits\":{hits},\"mray_per_s\":{mray_s:.4},\
             \"note\":\"release perf gate lives in perf-CI (fz2.4)\"}}",
            std::env::consts::OS,
            std::env::consts::ARCH,
            n * n
        );
        assert!(hits > 0, "the turntable frame sees the sphere");
        verdict(
            "rb-005",
            "primary-ray throughput measured and ledgered with build/machine labels; \
             the Mray/s TARGET is the perf-CI lane's gate",
        );
    });
}
