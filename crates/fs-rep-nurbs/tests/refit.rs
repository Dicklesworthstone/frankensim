//! SDF→NURBS refit conformance (the wqd.12 bead; runs under the
//! `nurbs-refit` feature). Acceptance: the NURBS→SDF→NURBS round-trip
//! recovers shape within retained sampled/analytic estimates (near-exactly on
//! unblended regions); Boolean-then-refit produces sampled sheaf
//! interface-agreement evidence on a CSG fixture (no continuum claim); seam continuity within tolerance
//! with exact G⁰; reports retain their sampled/no-certificate boundary;
//! localized fit residuals warn with locations instead of being silently
//! smoothed (without claiming that residual magnitude diagnoses thinness); the
//! patch-density budget knob trades fidelity monotonically.
#![cfg(feature = "nurbs-refit")]

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Chart, ChartSample, Point3, SheafComplex, SheafVerdict};
use fs_rep_nurbs::refit::{RefitConfig, refit_radial};
use fs_rep_nurbs::sdf::{Orientation, ShellSdf, ShellSdfChart};
use fs_rep_nurbs::{KnotVector, NurbsSurface};
use std::sync::atomic::{AtomicUsize, Ordering};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-rep-nurbs/refit\",\"case\":\"{case}\",\"verdict\":\"pass\",\
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
                seed: 11,
                kernel_id: 4,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

const S2: f64 = std::f64::consts::FRAC_1_SQRT_2;

/// Exact unit-sphere NURBS (the wqd.11 fixture, outward normals).
fn sphere_nurbs() -> NurbsSurface<f64> {
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
    let cw = |i: usize| if i.is_multiple_of(2) { 1.0 } else { S2 };
    let profile: [([f64; 2], f64); 5] = [
        ([0.0, -1.0], 1.0),
        ([1.0, -1.0], S2),
        ([1.0, 0.0], 1.0),
        ([1.0, 1.0], S2),
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
    let ku = KnotVector::new(
        vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        2,
    )
    .expect("ku");
    let kv = KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0], 2).expect("kv");
    NurbsSurface::new(ku, kv, &points, &weights).expect("sphere")
}

#[test]
fn rf_001_round_trip_through_the_real_converter() {
    // NURBS → SDF via the wqd.11 converter → NURBS via this bead.
    let shell =
        ShellSdf::new(vec![sphere_nurbs()], vec![None], Orientation::Outward).expect("shell");
    let sdf = |q: [f64; 3]| {
        let query = shell.distance(q, 5e-4, 300).expect("query");
        let sign = if (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt() < 1.0 {
            -1.0
        } else {
            1.0
        };
        sign * f64::midpoint(query.lower, query.upper)
    };
    let config = RefitConfig {
        samples_u: 24,
        samples_v: 24,
        ..RefitConfig::default()
    };
    let refit = refit_radial(&sdf, [0.0, 0.0, 0.0], 2.0, &config).expect("refit");
    // Recovery: fitted points sit on the unit sphere (unblended region —
    // near-exact at this density).
    let mut worst = 0.0f64;
    for a in 0..20 {
        for b in 1..20 {
            let (u, v) = (f64::from(a) / 20.0, f64::from(b) / 20.0);
            let p = refit.surface.eval(u, v).expect("eval");
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            worst = worst.max((r - 1.0).abs());
        }
    }
    assert!(worst < 5e-3, "round-trip radius recovery: {worst}");
    // This fixture is intentionally signed-distance-like, so its own test may
    // combine field residual and geometric probe spacing under that local
    // modeling assumption. The generic report keeps them separate.
    let fixture_unit_lipschitz_estimate =
        refit.report.spline_to_field_sampled + refit.report.spline_probe_spacing_estimate;
    assert!(
        fixture_unit_lipschitz_estimate < 8e-2,
        "fixture-conditional field-residual estimate closes: {fixture_unit_lipschitz_estimate}"
    );
    assert!(
        refit.report.warnings.is_empty(),
        "no retained sphere fit residual exceeds the warning threshold"
    );
    verdict(
        "rf-001",
        "NURBS->SDF->NURBS on the unit sphere: radius recovered to 5e-3; the \
         non-rigorous spline->field residual estimate closes",
    );
}

/// A chart adapter for an analytic CSG field (Booleans route through
/// F-rep — min/max ARE the Boolean policy).
struct CsgChart<F: Fn([f64; 3]) -> f64 + Send + Sync> {
    field: F,
    bound: f64,
}

impl<F: Fn([f64; 3]) -> f64 + Send + Sync> std::fmt::Debug for CsgChart<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("CsgChart")
    }
}

impl<F: Fn([f64; 3]) -> f64 + Send + Sync> Chart for CsgChart<F> {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let d = (self.field)([x.x, x.y, x.z]);
        ChartSample {
            signed_distance: d,
            gradient: None,
            lipschitz: Some(1.0),
            // This test adapter computes ordinary floating arithmetic; a
            // singleton Estimate records the observation without asserting
            // an exact-real or trusted-kernel proof.
            error: fs_evidence::NumericalCertificate::estimate(d, d),
        }
    }

    fn support(&self) -> fs_geom::Aabb {
        fs_geom::Aabb::new(
            Point3::new(-self.bound, -self.bound, -self.bound),
            Point3::new(self.bound, self.bound, self.bound),
        )
    }

    fn name(&self) -> &'static str {
        "test/csg"
    }
}

#[test]
fn rf_002_boolean_then_refit_sampled_interface_agreement() {
    with_cx(|cx| {
        // CSG union of two spheres (the F-rep Boolean), refit to NURBS,
        // then compare sampled interface evidence against the source field:
        // the refit chart and the CSG chart must agree at retained samples in
        // the shared surface band. This does not prove continuum watertightness.
        let union = |q: [f64; 3]| {
            let a = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt() - 1.0;
            let dx = q[0] - 0.55;
            let b = (dx * dx + q[1] * q[1] + q[2] * q[2]).sqrt() - 0.8;
            a.min(b)
        };
        let config = RefitConfig {
            nu: 16,
            nv: 16,
            samples_u: 48,
            samples_v: 48,
            lambda: 1e-5,
            ..RefitConfig::default()
        };
        let refit = refit_radial(&union, [0.2, 0.0, 0.0], 2.5, &config).expect("refit");
        // A minimum of exact sphere signed-distance functions is 1-Lipschitz,
        // so this fixture may combine the otherwise separate field-residual
        // and geometric probe-spacing estimates. This is test-local authority,
        // not authority supplied by the closure-based refit API.
        let conditional_field_estimate =
            refit.report.spline_to_field_sampled + refit.report.spline_probe_spacing_estimate;
        let sampled_tolerance_scale =
            conditional_field_estimate.max(refit.report.projected_target_to_spline_sampled);
        // Present both as charts and run the legacy-named sampled-interface
        // check.
        let refit_chart = ShellSdfChart::new(
            ShellSdf::new(vec![refit.surface], vec![None], Orientation::Outward).expect("shell"),
            1e-4,
            800,
            0.3,
        )
        .expect("admissible refit chart settings");
        let csg_chart = CsgChart {
            field: union,
            bound: 2.0,
        };
        let charts: Vec<&dyn Chart> = vec![&refit_chart, &csg_chart];
        let complex = SheafComplex::from_charts(&charts, cx).expect("finite patch supports");
        assert!(!complex.interfaces.is_empty(), "shared surface band found");
        let tol = 2.0 * sampled_tolerance_scale;
        let ev = complex.watertightness(tol);
        assert!(
            matches!(&ev.value, SheafVerdict::Unknown { .. }),
            "the measured NURBS distance chart must not be laundered into interval-verified seam authority: {:?}",
            ev.value
        );
        verdict(
            "rf-002",
            "CSG union -> refit retains a composite sampled residual, while the sheaf gate \
             correctly stays Unknown because the NURBS distance chart is measured rather \
             than outward-enclosed; no Hausdorff or continuum claim",
        );
    });
}

#[test]
fn rf_003_seam_g0_exact_g1_measured() {
    let sdf = |q: [f64; 3]| (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt() - 1.0;
    let refit = refit_radial(&sdf, [0.0, 0.0, 0.0], 2.0, &RefitConfig::default()).expect("refit");
    // G0: the tied control columns make the clamped endpoint evaluations
    // identical, not merely close at a nearby parameter.
    for b in 0..12 {
        let v = (f64::from(b) + 0.5) / 12.0;
        let p0 = refit.surface.eval(0.0, v).expect("eval");
        let p1 = refit.surface.eval(1.0, v).expect("eval");
        assert_eq!(p0, p1, "tied G0 seam endpoint at v={v}");
    }
    // G1: measured and small on a smooth field.
    assert!(
        refit.report.seam_g1_max < 1e-2,
        "seam tangent deviation: {}",
        refit.report.seam_g1_max
    );
    assert!(
        refit.report.seam_g1_excludes_v_endpoints,
        "the interior seam diagnostic must explicitly disclose its pole/endpoint exclusion"
    );
    verdict(
        "rf-003",
        "seam G0 exact by control tying; G1 deviation measured < 1e-2",
    );
}

#[test]
fn rf_004_localized_residual_warns_on_smoothed_spike() {
    // A sphere with a thin radial spike (capsule toward +x): far below
    // patch resolution at the default density.
    let spiky = |q: [f64; 3]| {
        let sphere = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt() - 1.0;
        // A capsule spur from (1,0,0) to (1.6,0,0), radius 0.12: an
        // azimuthal feature ~0.18 rad wide against a control spacing of
        // ~0.7 rad — below PATCH resolution, visible to the samples.
        let t = (q[0] - 1.0).clamp(0.0, 0.6);
        let spike = ((q[0] - 1.0 - t).powi(2) + q[1] * q[1] + q[2] * q[2]).sqrt() - 0.12;
        sphere.min(spike)
    };
    let refit = refit_radial(&spiky, [0.0, 0.0, 0.0], 2.2, &RefitConfig::default()).expect("refit");
    assert!(
        !refit.report.warnings.is_empty(),
        "a sub-resolution spike must leave a localized residual warning"
    );
    // The warnings localize to the spike (azimuth ~ 0, equator v ~ 0.5).
    let near_spike = refit
        .report
        .warnings
        .iter()
        .all(|w| (w.uv[0] < 0.1 || w.uv[0] > 0.9) && (w.uv[1] - 0.5).abs() < 0.15);
    assert!(
        near_spike,
        "residual warnings localized: {:?}",
        refit.report.warnings
    );
    // And the report says the fit did NOT follow the spike.
    assert!(
        refit.report.max_residual > 0.1,
        "the residual names the miss: {}",
        refit.report.max_residual
    );
    verdict(
        "rf-004",
        "sub-resolution spike produces localized fit-residual warnings without overdiagnosing their cause",
    );
}

#[test]
fn rf_005_patch_density_budget_knob() {
    let sdf = |q: [f64; 3]| {
        // A gently lobed star-shaped blob (needs real fitting power).
        let r = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
        if r < 1e-12 {
            return -1.0;
        }
        let bump = 0.15 * (3.0 * q[2] / r).sin() * (2.0 * q[0] / r).cos();
        r - (1.0 + bump)
    };
    let coarse = refit_radial(
        &sdf,
        [0.0, 0.0, 0.0],
        2.0,
        &RefitConfig {
            nu: 6,
            nv: 6,
            samples_u: 24,
            samples_v: 24,
            ..RefitConfig::default()
        },
    )
    .expect("coarse");
    let fine = refit_radial(
        &sdf,
        [0.0, 0.0, 0.0],
        2.0,
        &RefitConfig {
            nu: 16,
            nv: 16,
            samples_u: 48,
            samples_v: 48,
            ..RefitConfig::default()
        },
    )
    .expect("fine");
    assert!(
        fine.report.spline_to_field_sampled < coarse.report.spline_to_field_sampled,
        "more patches, better fidelity: fine {} vs coarse {}",
        fine.report.spline_to_field_sampled,
        coarse.report.spline_to_field_sampled
    );
    println!(
        "{{\"metric\":\"refit-budget-knob\",\"coarse_sampled\":{:.3e},\"fine_sampled\":{:.3e},\
         \"coarse_probe_spacing_estimate\":{:.3e},\"fine_probe_spacing_estimate\":{:.3e}}}",
        coarse.report.spline_to_field_sampled,
        fine.report.spline_to_field_sampled,
        coarse.report.spline_probe_spacing_estimate,
        fine.report.spline_probe_spacing_estimate
    );
    verdict(
        "rf-005",
        "the patch-density knob trades cost for fidelity monotonically (ledgered)",
    );
}

#[test]
fn rf_006_probe_work_is_admitted_before_calling_the_field() {
    let calls = AtomicUsize::new(0);
    let field = |_: [f64; 3]| {
        calls.fetch_add(1, Ordering::Relaxed);
        -1.0
    };
    let refusal = refit_radial(
        &field,
        [0.0, 0.0, 0.0],
        2.0,
        &RefitConfig {
            nu: 21,
            nv: 21,
            degree: 20,
            samples_u: 1,
            samples_v: 1,
            probe: 2_048,
            ..RefitConfig::default()
        },
    );
    assert!(
        refusal.is_err(),
        "billions of degree-sensitive probe operations must exceed the static legacy cap"
    );
    assert_eq!(
        calls.load(Ordering::Relaxed),
        0,
        "work admission must fail before invoking an arbitrary caller closure"
    );
    verdict(
        "rf-006",
        "degree-sensitive dense-probe work refuses before the first field call",
    );
}
