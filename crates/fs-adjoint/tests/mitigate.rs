//! Meshing-mitigation conformance (the bk0o.2 bead; runs under the
//! `diff-mitigations` feature). Acceptance: a gradient query routes to
//! the SDF/spline path when one exists; Hadamard boundary forms compute
//! shape gradients without mesh sensitivities (vs perturbation-resolve);
//! an unavoidable remesh yields an estimated-color gradient with a
//! discontinuity flag — never a silently-verified one; the flagged
//! discontinuity is REAL (the FD falsifier confirms it); router
//! tie-breaks are deterministic.
#![cfg(feature = "diff-mitigations")]

use std::collections::BTreeSet;

use fs_adjoint::mitigate::{GradientGrade, grade_ops, plan_gradient_route};
use fs_adjoint::transpose::fd_falsifier;
use fs_evidence::Color;
use fs_feec::kuhn_cube;
use fs_feec::whitney::element_geometry;
use fs_geom::{ConverterSpec, ErrorModel, MemoryCostOracle, RoutePlanError, RouteRequest, Router};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-adjoint/mitigate\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn edge(name: &str, from: &str, to: &str, cost: f64) -> ConverterSpec {
    ConverterSpec {
        from: from.to_string(),
        to: to.to_string(),
        name: name.to_string(),
        base_cost_s: cost,
        error: ErrorModel::AdditiveAbs(1e-6),
        certified: true,
    }
}

/// A router with two routes nurbs → sdf: a CHEAP mesh path (through a
/// non-differentiable remesh edge) and a pricier direct smooth path.
fn two_path_router() -> (Router, BTreeSet<String>) {
    let mut router = Router::new();
    router
        .register(edge("nurbs->mesh/remesh-v1", "nurbs", "mesh", 0.5))
        .expect("r");
    router
        .register(edge("mesh->sdf/rasterize-v1", "mesh", "sdf", 0.5))
        .expect("r");
    router
        .register(edge("nurbs->sdf/bezier-clip-v1", "nurbs", "sdf", 5.0))
        .expect("r");
    let mut nd = BTreeSet::new();
    nd.insert("nurbs->mesh/remesh-v1".to_string());
    (router, nd)
}

fn request() -> RouteRequest {
    request_with_cost(1e9)
}

fn request_with_cost(max_cost_s: f64) -> RouteRequest {
    RouteRequest {
        from: "nurbs".to_string(),
        to: "sdf".to_string(),
        scale: 1.0,
        max_abs_error: 1e-3,
        max_cost_s,
    }
}

#[test]
fn mg_001_gradient_queries_prefer_the_smooth_path() {
    let (router, nd) = two_path_router();
    let oracle = MemoryCostOracle::new();
    // WITHOUT a gradient request the cheap mesh path wins.
    let plain = router.plan(&request(), &oracle).expect("plain route");
    assert_eq!(
        plain.edges(),
        &[
            "nurbs->mesh/remesh-v1".to_string(),
            "mesh->sdf/rasterize-v1".to_string()
        ],
        "cost-only planning takes the cheap mesh path"
    );
    // WITH a gradient request the differentiability term dominates.
    let (plan, grade) =
        plan_gradient_route(&router, &request(), &oracle, &nd).expect("gradient route");
    assert_eq!(
        plan.edges(),
        &["nurbs->sdf/bezier-clip-v1".to_string()],
        "the smooth spline path is preferred under a gradient request"
    );
    assert!(
        matches!(grade, GradientGrade::Smooth { .. }),
        "smooth path, smooth grade: {grade:?}"
    );
    verdict(
        "mg-001",
        "same router, same costs: plain query takes the cheap remesh path, gradient \
         query routes onto the smooth spline edge",
    );
}

#[test]
fn mg_002_unavoidable_remesh_downgrades_never_fakes() {
    // A router where EVERY path crosses the remesh edge.
    let mut router = Router::new();
    router
        .register(edge("nurbs->mesh/remesh-v1", "nurbs", "mesh", 0.5))
        .expect("r");
    router
        .register(edge("mesh->sdf/rasterize-v1", "mesh", "sdf", 0.5))
        .expect("r");
    let mut nd = BTreeSet::new();
    nd.insert("nurbs->mesh/remesh-v1".to_string());
    let oracle = MemoryCostOracle::new();
    let (plan, grade) = plan_gradient_route(&router, &request_with_cost(2.0), &oracle, &nd)
        .expect("the real budget admits the unavoidable route");
    assert_eq!(plan.edges().len(), 2, "the only path is taken");
    match &grade {
        GradientGrade::EstimatedWithDiscontinuity {
            crossing, color, ..
        } => {
            assert_eq!(crossing, &vec!["nurbs->mesh/remesh-v1".to_string()]);
            match color {
                Color::Estimated {
                    estimator,
                    dispersion,
                } => {
                    assert!(estimator.contains("remesh"), "{estimator}");
                    assert!(
                        dispersion.is_infinite(),
                        "no spread claim across a topology event"
                    );
                }
                other => panic!("must be Estimated: {other:?}"),
            }
        }
        GradientGrade::Smooth { .. } => panic!("must downgrade, never fake: {grade:?}"),
    }
    // The tape-level twin: grade_ops applies the same rule.
    let g2 = grade_ops(&["convert/restrict", "mesh/remesh", "solver/spd"], &{
        let mut s = BTreeSet::new();
        s.insert("mesh/remesh".to_string());
        s
    });
    assert!(g2.discontinuity().is_some(), "op-level path flagged too");
    verdict(
        "mg-002",
        "when every path crosses the remesh, the gradient is Estimated with an \
         infinite-dispersion color and the edge named in the flag",
    );
}

#[test]
fn mg_003_hadamard_matches_perturbation_resolve() {
    // Volume shape gradient on a tet ball vs actually perturbing the
    // boundary and re-measuring volume: the mesh-free boundary form
    // agrees with perturbation-resolve without any mesh sensitivities.
    let (complex, positions) = kuhn_cube(3);
    let velocity = |p: [f64; 3]| -> [f64; 3] {
        // A smooth outward-ish field.
        [0.3 + 0.1 * p[1], -0.2 + 0.05 * p[0], 0.15 * p[2] + 0.05]
    };
    let hadamard = fs_adjoint::volume_shape_gradient(&complex, &positions, &velocity);
    // Perturbation-resolve: move EVERY vertex by t·V, re-measure volume.
    let volume_at = |t: f64| -> f64 {
        let moved: Vec<[f64; 3]> = positions
            .iter()
            .map(|p| {
                let v = velocity(*p);
                [p[0] + t * v[0], p[1] + t * v[1], p[2] + t * v[2]]
            })
            .collect();
        let geo = element_geometry(&complex, &moved);
        geo.vol_signed.iter().map(|v| v.abs()).sum()
    };
    let fd = fd_falsifier(
        &|x: &[f64]| volume_at(x[0]),
        &[0.0],
        &[1.0],
        hadamard,
        1e-5,
        1e-6,
    );
    assert!(
        fd.consistent,
        "Hadamard boundary form must match perturbation-resolve: {fd:?}"
    );
    verdict(
        "mg-003",
        "volume shape gradient by the Hadamard boundary form matches perturbation-\
         resolve through the FD falsifier (no mesh sensitivities computed)",
    );
}

#[test]
fn mg_004_the_flag_is_real_and_smooth_paths_check_out() {
    // A remesh-like objective: smooth in x plus a STEP when a "cell
    // budget" threshold crosses (the topology event). Across the step
    // the FD falsifier legitimately fails — which is exactly why the
    // flag exists; on the smooth side it passes.
    let objective = |x: &[f64]| -> f64 {
        let smooth = x[0] * x[0] + 0.5 * x[0];
        let cells = if x[0] > 0.5 { 9.0 } else { 8.0 }; // remesh event
        smooth + 0.01 * cells
    };
    // Smooth region: adjoint dd = 2x + 0.5 at x = 0.2.
    let smooth_dd = 2.0 * 0.2 + 0.5;
    let ok = fd_falsifier(&objective, &[0.2], &[1.0], smooth_dd, 1e-5, 1e-7);
    assert!(ok.consistent, "smooth side agrees: {ok:?}");
    // Straddling the event at x = 0.5: central FD sees the jump; the
    // smooth-formula adjoint sits inside a WIDENED band: the
    // Richardson self-error EXPLODES (FD at h and h/2 disagree by the
    // jump/h scale), so the conditioning-aware falsifier correctly
    // refuses to false-fire (bk0o.1's round-3 contract) — and that
    // explosion IS the measurable signature that the declared
    // discontinuity is real.
    assert!(
        (ok.fd_coarse - ok.fd_fine).abs() < 1e-6,
        "smooth side: negligible FD self-error"
    );
    let step_dd = 2.0 * 0.5 + 0.5;
    let at_step = fd_falsifier(&objective, &[0.5], &[1.0], step_dd, 1e-5, 1e-7);
    let self_error = (at_step.fd_coarse - at_step.fd_fine).abs();
    assert!(
        self_error > 100.0,
        "the topology event shows in the FD self-error: {self_error}"
    );
    assert!(
        at_step.consistent,
        "the falsifier does not false-fire on the ill-conditioned point: {at_step:?}"
    );
    verdict(
        "mg-004",
        "smooth side: tiny self-error, falsifier agrees; at the step the self-error \
         explodes (the real discontinuity) while the conditioning-aware falsifier \
         correctly declines to false-fire",
    );
}

#[test]
fn mg_005_deterministic_tie_break() {
    // Two EQUAL-fitness smooth paths: the router's choice must be
    // deterministic across repeated plans (the review-round-3 boundary).
    let mut router = Router::new();
    router
        .register(edge("nurbs->sdf/path-a", "nurbs", "sdf", 2.0))
        .expect("r");
    router
        .register(edge("nurbs->sdf/path-b", "nurbs", "sdf", 2.0))
        .expect("r");
    let oracle = MemoryCostOracle::new();
    let nd = BTreeSet::new();
    let first = plan_gradient_route(&router, &request(), &oracle, &nd)
        .expect("route")
        .0
        .edges()
        .to_vec();
    for _ in 0..5 {
        let again = plan_gradient_route(&router, &request(), &oracle, &nd)
            .expect("route")
            .0
            .edges()
            .to_vec();
        assert_eq!(again, first, "tie-break is deterministic");
    }
    verdict(
        "mg-005",
        "equal-fitness paths: five repeated plans pick the identical route",
    );
}

#[test]
fn mg_006_smooth_preference_is_lexicographic_not_a_fixed_cost_penalty() {
    let mut router = Router::new();
    router
        .register(edge("nurbs->sdf/remesh-v1", "nurbs", "sdf", 1.0))
        .expect("register remesh route");
    router
        .register(edge(
            "nurbs->sdf/smooth-expensive-v1",
            "nurbs",
            "sdf",
            2_000_000.0,
        ))
        .expect("register smooth route");
    let mut nd = BTreeSet::new();
    nd.insert("nurbs->sdf/remesh-v1".to_string());
    let oracle = MemoryCostOracle::new();

    let (plan, grade) = plan_gradient_route(&router, &request_with_cost(3_000_000.0), &oracle, &nd)
        .expect("the expensive smooth route is still within the real budget");
    assert_eq!(
        plan.edges(),
        &["nurbs->sdf/smooth-expensive-v1".to_string()]
    );
    assert!(matches!(grade, GradientGrade::Smooth { .. }));
    verdict(
        "mg-006",
        "an admissible smooth route wins lexicographically even when its real cost exceeds the old fixed penalty",
    );
}

#[test]
fn mg_007_smooth_budget_failure_falls_back_under_the_original_budget() {
    let mut router = Router::new();
    router
        .register(edge("nurbs->sdf/remesh-v1", "nurbs", "sdf", 1.0))
        .expect("register remesh route");
    router
        .register(edge("nurbs->sdf/smooth-v1", "nurbs", "sdf", 5.0))
        .expect("register smooth route");
    let mut nd = BTreeSet::new();
    nd.insert("nurbs->sdf/remesh-v1".to_string());
    let oracle = MemoryCostOracle::new();

    let (plan, grade) = plan_gradient_route(&router, &request_with_cost(2.0), &oracle, &nd)
        .expect("fallback must retain the request's real cost budget");
    assert_eq!(plan.edges(), &["nurbs->sdf/remesh-v1".to_string()]);
    assert!(matches!(
        grade,
        GradientGrade::EstimatedWithDiscontinuity { .. }
    ));
    verdict(
        "mg-007",
        "an over-budget smooth route falls back to the affordable remesh route with an explicit discontinuity grade",
    );
}

#[test]
fn mg_008_no_route_preserves_the_router_refusal() {
    let router = Router::new();
    let oracle = MemoryCostOracle::new();
    let error = plan_gradient_route(&router, &request_with_cost(2.0), &oracle, &BTreeSet::new())
        .expect_err("an empty graph has no route in either pass");
    match error {
        RoutePlanError::Infeasible(refusal) => {
            assert_eq!(refusal.binding, fs_geom::Binding::NoPath);
        }
        other => panic!("empty graph must produce an infeasibility, got {other:?}"),
    }
    verdict(
        "mg-008",
        "smooth-first planning preserves the original structured no-path refusal",
    );
}
