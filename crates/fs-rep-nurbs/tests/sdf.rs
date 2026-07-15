//! NURBS→SDF converter conformance (the wqd.11 bead; runs under the
//! `nurbs-sdf` feature). Acceptance: against exactly-representable
//! revolution fixtures (sphere/cylinder/torus as rational quadratics)
//! measured brackets are cross-checked against analytic distance everywhere
//! sampled; near-trim behavior widens honestly on adversarial trims;
//! sign assignment follows declared orientation with the unsigned
//! fallback named; G3 frame invariance; adaptive tiled generation with
//! throughput evidence.
#![cfg(feature = "nurbs-sdf")]

use asupersync::types::Budget;
use fs_evidence::NumericalKind;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Chart, Differentiability, Point3};
use fs_rep_nurbs::sdf::{Orientation, ShellSdf, ShellSdfChart, generate_tile};
use fs_rep_nurbs::{KnotVector, NurbsCurve, NurbsSurface, Rat, TrimLoop, TrimmedPatch};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-rep-nurbs/sdf\",\"case\":\"{case}\",\"verdict\":\"pass\",\
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
                seed: 7,
                kernel_id: 3,
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

/// The standard 9-point rational-quadratic full circle template:
/// (cos, sin) pairs at radius 1 with weights [1, s, 1, s, …].
const CIRCLE_XY: [[f64; 2]; 9] = [
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

fn circle_knots() -> KnotVector<f64> {
    KnotVector::new(
        vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        2,
    )
    .expect("circle knots")
}

fn circle_weight(i: usize) -> f64 {
    if i.is_multiple_of(2) { 1.0 } else { S2 }
}

/// Surface of revolution about z: profile points (radius_j, z_j) with
/// weights wv_j swept by the circle template (exactly representable —
/// the classic rational construction).
fn revolve(profile: &[([f64; 2], f64)], knots_v: KnotVector<f64>) -> NurbsSurface<f64> {
    let mut points: Vec<Vec<[f64; 3]>> = Vec::with_capacity(9);
    let mut weights: Vec<Vec<f64>> = Vec::with_capacity(9);
    for (i, c) in CIRCLE_XY.iter().enumerate() {
        let mut prow = Vec::with_capacity(profile.len());
        let mut wrow = Vec::with_capacity(profile.len());
        for &([radius, z], wv) in profile {
            prow.push([radius * c[0], radius * c[1], z]);
            wrow.push(circle_weight(i) * wv);
        }
        points.push(prow);
        weights.push(wrow);
    }
    NurbsSurface::new(circle_knots(), knots_v, &points, &weights).expect("revolution")
}

/// Exact unit sphere (half-circle profile revolved).
fn sphere() -> NurbsSurface<f64> {
    let half = KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0], 2).expect("half");
    // Profile runs SOUTH -> NORTH so du x dv points outward.
    revolve(
        &[
            ([0.0, -1.0], 1.0),
            ([1.0, -1.0], S2),
            ([1.0, 0.0], 1.0),
            ([1.0, 1.0], S2),
            ([0.0, 1.0], 1.0),
        ],
        half,
    )
}

fn collapsed_surface() -> NurbsSurface<f64> {
    let line = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line");
    let points = vec![vec![[0.0; 3]; 2]; 2];
    let weights = vec![vec![1.0; 2]; 2];
    NurbsSurface::new(line.clone(), line, &points, &weights).expect("collapsed surface")
}

/// Cylinder SIDE surface (radius r, z in [0, h]).
fn cylinder(r: f64, h: f64) -> NurbsSurface<f64> {
    let line = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line");
    revolve(&[([r, 0.0], 1.0), ([r, h], 1.0)], line)
}

/// Torus (major R, minor r): full-circle profile offset by R.
fn torus(big: f64, small: f64) -> NurbsSurface<f64> {
    let profile: Vec<([f64; 2], f64)> = CIRCLE_XY
        .iter()
        .enumerate()
        .map(|(i, c)| (([big + small * c[0], small * c[1]]), circle_weight(i)))
        .collect();
    revolve(&profile, circle_knots())
}

fn lcg(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / (1u64 << 53) as f64
}

fn sample_box(state: &mut u64, half: f64) -> [f64; 3] {
    [
        (lcg(state) - 0.5) * 2.0 * half,
        (lcg(state) - 0.5) * 2.0 * half,
        (lcg(state) - 0.5) * 2.0 * half,
    ]
}

#[test]
fn ns_001_containment_battery_vs_analytic() {
    type Case<'a> = (&'a ShellSdf, &'a dyn Fn([f64; 3]) -> f64, f64, &'a str);
    let sph = ShellSdf::new(vec![sphere()], vec![None], Orientation::Outward).expect("shell");
    let cyl =
        ShellSdf::new(vec![cylinder(0.7, 1.2)], vec![None], Orientation::Unknown).expect("shell");
    let tor =
        ShellSdf::new(vec![torus(1.0, 0.3)], vec![None], Orientation::Outward).expect("shell");
    let analytic_sphere =
        |q: [f64; 3]| ((q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt() - 1.0).abs();
    let analytic_cyl = |q: [f64; 3]| {
        let rho = (q[0] * q[0] + q[1] * q[1]).sqrt();
        let dz = (0.0f64).max(q[2] - 1.2).max(-q[2]);
        ((rho - 0.7).powi(2) + dz * dz).sqrt()
    };
    let analytic_torus = |q: [f64; 3]| {
        let rho = (q[0] * q[0] + q[1] * q[1]).sqrt();
        (((rho - 1.0).powi(2) + q[2] * q[2]).sqrt() - 0.3).abs()
    };
    let cases: [Case<'_>; 3] = [
        (&sph, &analytic_sphere, 1.6, "sphere"),
        (&cyl, &analytic_cyl, 1.8, "cylinder"),
        (&tor, &analytic_torus, 1.8, "torus"),
    ];
    let mut state = 0x5eed_cafe;
    let mut worst_width = 0.0f64;
    for (shell, analytic, half, name) in cases {
        for _ in 0..50 {
            let q = sample_box(&mut state, half);
            let query = shell.distance(q, 1e-6, 4000).expect("query");
            let truth = analytic(q);
            assert!(
                query.lower - 1e-7 <= truth && truth <= query.upper + 1e-7,
                "{name}: measured bracket [{}, {}] missed analytic {truth} at {q:?}",
                query.lower,
                query.upper
            );
            worst_width = worst_width.max(query.upper - query.lower);
        }
    }
    assert!(
        worst_width < 5e-3,
        "brackets converge: worst width {worst_width}"
    );
    verdict(
        "ns-001",
        "150 sampled points across sphere/cylinder/torus: measured brackets contain \
         the sampled analytic oracle; worst observed width < 5e-3 (no enclosure claim)",
    );
}

#[test]
fn ns_002_sign_follows_declared_orientation() {
    with_cx(|cx| {
        let signed = ShellSdfChart::new(
            ShellSdf::new(vec![sphere()], vec![None], Orientation::Outward).expect("shell"),
            1e-6,
            4000,
            0.5,
        )
        .expect("admissible signed chart settings");
        let inside = signed.eval(Point3::new(0.3, 0.1, -0.2), cx);
        assert!(inside.signed_distance < 0.0, "inside is negative");
        let center = signed.eval(Point3::new(0.0, 0.0, 0.0), cx);
        assert!(
            center.signed_distance < 0.0,
            "a singular-pole closest witness must use a regular oriented fallback at the sphere center: {center:?}"
        );
        let axial_inside = signed.eval(Point3::new(0.0, 0.0, -0.5), cx);
        assert!(
            axial_inside.signed_distance < 0.0,
            "inside axial query at a pole-degenerate closest parameter must remain negative"
        );
        let axial_outside = signed.eval(Point3::new(0.0, 0.0, -1.5), cx);
        assert!(
            axial_outside.signed_distance > 0.0,
            "outside axial query at a pole-degenerate closest parameter must remain positive"
        );
        let outside = signed.eval(Point3::new(1.5, 0.0, 0.0), cx);
        assert!(outside.signed_distance > 0.0, "outside is positive");
        assert!(
            (outside.signed_distance - 0.5).abs() < 1e-4,
            "sphere distance at 1.5: {}",
            outside.signed_distance
        );
        let g = outside.gradient.expect("gradient off the surface");
        assert!(g.x > 0.99, "gradient points outward: {g:?}");
        assert_eq!(signed.name(), "nurbs-sdf/estimated-signed");
        assert_eq!(
            signed.differentiability(),
            Differentiability::Unknown,
            "finite-budget witness switching grants no continuity claim"
        );
        assert_eq!(outside.error.kind, NumericalKind::Estimate);
        assert!(
            outside.lipschitz.is_none(),
            "measured field grants no Lipschitz authority"
        );
        assert_eq!(
            inside.lipschitz, None,
            "the measured field must not claim a certified Lipschitz constant"
        );
        // Unknown orientation: unsigned field, named as such.
        let unsigned = ShellSdfChart::new(
            ShellSdf::new(vec![sphere()], vec![None], Orientation::Unknown).expect("shell"),
            1e-6,
            4000,
            0.5,
        )
        .expect("admissible unsigned chart settings");
        let u_inside = unsigned.eval(Point3::new(0.3, 0.1, -0.2), cx);
        assert!(
            u_inside.signed_distance > 0.0,
            "no sign claim: non-negative"
        );
        assert_eq!(unsigned.name(), "nurbs-sdf/estimated-unsigned");
        let degenerate = ShellSdfChart::new(
            ShellSdf::new(vec![collapsed_surface()], vec![None], Orientation::Outward)
                .expect("degenerate shell remains diagnosable"),
            1e-6,
            8,
            0.5,
        )
        .expect("admissible degenerate chart settings");
        let no_sign = degenerate.eval(Point3::new(0.5, 0.0, 0.0), cx);
        assert_eq!(no_sign.error.kind, NumericalKind::NoClaim);
        assert!(
            no_sign.signed_distance.is_infinite(),
            "zero tangents/normals must fail closed instead of defaulting to outside"
        );
        verdict(
            "ns-002",
            "outward orientation signs correctly with outward gradients; unknown \
             orientation stays unsigned and says so in the chart name",
        );
    });
}

/// A closed degree-1 rational polyline loop in parameter space.
fn poly_loop(pts: &[[i64; 2]], scale_den: i64) -> TrimLoop {
    let n = pts.len();
    let mut knots = vec![Rat::int(0), Rat::int(0)];
    for k in 1..n {
        knots.push(Rat::new(k as i128, n as i128));
    }
    knots.push(Rat::int(1));
    knots.push(Rat::int(1));
    let kv = KnotVector::new(knots, 1).expect("polyline knots");
    let mut points: Vec<[Rat; 2]> = pts
        .iter()
        .map(|p| {
            [
                Rat::new(i128::from(p[0]), i128::from(scale_den)),
                Rat::new(i128::from(p[1]), i128::from(scale_den)),
            ]
        })
        .collect();
    points.push(points[0]);
    let weights = vec![Rat::int(1); points.len()];
    TrimLoop::new(NurbsCurve::new(kv, &points, &weights).expect("loop")).expect("closed")
}

#[test]
fn ns_003_trim_downgrade_is_honest() {
    // Keep only the NORTHERN parameter half (v in [1/2, 1] is north because
    // the sphere profile runs south to north): CCW square in (u, v).
    let keep_north = poly_loop(&[[0, 8], [16, 8], [16, 16], [0, 16]], 16);
    let trim = TrimmedPatch {
        loops: vec![keep_north],
        max_subdivision: 24,
    };
    let mut malformed_trim = trim.clone();
    malformed_trim.loops[0].curve.knots.knots.clear();
    assert!(
        ShellSdf::new(
            vec![sphere()],
            vec![Some(malformed_trim)],
            Orientation::Outward,
        )
        .is_err(),
        "shell ingestion must revalidate caller-mutable trim structure"
    );
    let shell =
        ShellSdf::new(vec![sphere()], vec![Some(trim)], Orientation::Outward).expect("shell");
    // A point nearest the SOUTH pole: its closest point is trimmed away.
    let south = shell
        .distance([0.0, 0.99, -0.99], 1e-6, 4000)
        .expect("query");
    assert!(
        south.trim_downgrade,
        "southern closest point is trimmed away"
    );
    assert!(
        south.upper.is_infinite(),
        "a trimmed-away point is not a kept-surface upper witness"
    );
    let alternate = ShellSdf::new(
        vec![sphere(), sphere()],
        vec![
            Some(TrimmedPatch {
                loops: vec![poly_loop(&[[0, 8], [16, 8], [16, 16], [0, 16]], 16)],
                max_subdivision: 24,
            }),
            None,
        ],
        Orientation::Outward,
    )
    .expect("shell with an untrimmed alternate");
    let alternate_south = alternate
        .distance([0.0, 0.99, -0.99], 1e-6, 4000)
        .expect("alternate query");
    assert!(
        !alternate_south.trim_downgrade && alternate_south.upper.is_finite(),
        "a valid kept-surface witness must outrank a closer trimmed-away diagnostic"
    );
    // A point nearest the NORTH pole: kept region, finite measured estimate.
    let north = shell
        .distance([0.0, 0.99, 0.99], 1e-6, 4000)
        .expect("query");
    assert!(!north.trim_downgrade, "northern closest point is kept");
    let expect = (0.99f64 * 0.99 * 2.0).sqrt() - 1.0;
    assert!(
        (north.upper - expect).abs() < 1e-4,
        "{} vs {expect}",
        north.upper
    );
    // The chart widens the trim-downgraded estimate instead of lying.
    with_cx(|cx| {
        let chart = ShellSdfChart::new(
            ShellSdf::new(
                vec![sphere()],
                vec![Some(TrimmedPatch {
                    loops: vec![poly_loop(&[[0, 8], [16, 8], [16, 16], [0, 16]], 16)],
                    max_subdivision: 24,
                })],
                Orientation::Outward,
            )
            .expect("shell"),
            1e-6,
            4000,
            0.5,
        )
        .expect("admissible trimmed chart settings");
        let s = chart.eval(Point3::new(0.0, 0.99, -0.99), cx);
        assert!(
            s.error.hi.is_infinite(),
            "downgraded estimate is widened, not asserted: {:?}",
            s.error
        );
        assert_eq!(s.error.kind, NumericalKind::NoClaim);
        assert!(
            s.signed_distance.is_infinite(),
            "a trimmed-away point must not survive as a finite nominal chart value"
        );
        let n = chart.eval(Point3::new(0.0, 0.99, 0.99), cx);
        assert!(
            n.error.hi.is_finite(),
            "kept region retains a finite estimate"
        );
        assert_eq!(n.error.kind, NumericalKind::Estimate);
    });
    verdict(
        "ns-003",
        "adversarial trim: closest-point-in-hole widens to an infinite upper estimate; \
         kept region retains a finite measured value near 0.4",
    );
}

#[test]
fn ns_004_frame_invariance_g3() {
    // Translate the whole shell (control points) and the query by the
    // same offset: the measured bracket should agree tightly.
    let base =
        ShellSdf::new(vec![torus(1.0, 0.3)], vec![None], Orientation::Outward).expect("shell");
    let mut moved_surface = torus(1.0, 0.3);
    let t = [13.0, -7.0, 3.5];
    for row in &mut moved_surface.cpw {
        for h in row.iter_mut() {
            let w = h[3];
            h[0] += t[0] * w;
            h[1] += t[1] * w;
            h[2] += t[2] * w;
        }
    }
    let moved =
        ShellSdf::new(vec![moved_surface], vec![None], Orientation::Outward).expect("shell");
    let mut state = 0xfeed_f00d;
    for _ in 0..25 {
        let q = sample_box(&mut state, 1.6);
        let a = base.distance(q, 1e-6, 4000).expect("base");
        let b = moved
            .distance([q[0] + t[0], q[1] + t[1], q[2] + t[2]], 1e-6, 4000)
            .expect("moved");
        assert!(
            (a.upper - b.upper).abs() < 1e-6,
            "translation invariance: {} vs {}",
            a.upper,
            b.upper
        );
    }
    verdict(
        "ns-004",
        "25 translated queries agree to 1e-6 (G3 metamorphic)",
    );
}

#[test]
fn ns_005_adaptive_tile_generation() {
    let chart = ShellSdfChart::new(
        ShellSdf::new(vec![sphere()], vec![None], Orientation::Outward).expect("shell"),
        1e-5,
        4000,
        0.5,
    )
    .expect("admissible adaptive chart settings");
    let aabb = fs_geom::Aabb::new(Point3::new(-1.4, -1.4, -1.4), Point3::new(1.4, 1.4, 1.4));
    let tile = generate_tile(&chart, &aabb, 8, 1e-5, 2000).expect("tile");
    assert_eq!(tile.values.len(), 512);
    assert_eq!(tile.downgraded, 0, "untrimmed shell never downgrades");
    assert_eq!(
        (
            tile.trim_downgraded,
            tile.sign_downgraded,
            tile.range_downgraded
        ),
        (0, 0, 0)
    );
    assert!(generate_tile(&chart, &aabb, 1, 1e-5, 1).is_err());
    assert!(generate_tile(&chart, &aabb, usize::MAX, 1e-5, 1).is_err());
    assert!(generate_tile(&chart, &aabb, 2, f64::NAN, 1).is_err());
    assert!(generate_tile(&chart, &aabb, 2, 1e-5, u32::MAX).is_err());
    assert!(
        ShellSdfChart::new(
            ShellSdf::new(vec![sphere()], vec![None], Orientation::Outward).expect("shell"),
            1e-5,
            1_048_576,
            0.0,
        )
        .is_err(),
        "chart construction must reject a configured query that can only fail its work ceiling"
    );
    let line = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("linear knots");
    let linear_surface = NurbsSurface::new(
        line.clone(),
        line,
        &vec![vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]; 2],
        &vec![vec![1.0; 2]; 2],
    )
    .expect("linear surface");
    let envelope_error = ShellSdfChart::new(
        ShellSdf::new(vec![linear_surface], vec![None], Orientation::Unknown)
            .expect("linear shell"),
        1e-5,
        1_048_576,
        0.0,
    )
    .expect_err("chart admission must share the delegated closest frontier envelope");
    assert!(
        envelope_error.to_string().contains("retained"),
        "closest frontier refusal must happen at chart construction: {envelope_error}"
    );
    assert!(
        ShellSdfChart::new(
            ShellSdf::new(vec![sphere()], vec![None], Orientation::Unknown).expect("shell"),
            f64::NAN,
            0,
            0.0,
        )
        .is_err(),
        "invalid chart settings must return a structured error, not panic"
    );
    let support_knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("support test knots");
    let extreme_surface = NurbsSurface::new(
        support_knots.clone(),
        support_knots,
        &vec![vec![[f64::MAX, 0.0, 0.0]; 2]; 2],
        &vec![vec![1.0; 2]; 2],
    )
    .expect("finite extreme surface");
    let extreme_shell = ShellSdf::new(vec![extreme_surface], vec![None], Orientation::Unknown)
        .expect("finite shell inputs remain diagnosable");
    assert!(
        extreme_shell.control_aabb(0.0).is_err(),
        "one-ULP outward support at f64::MAX must refuse infinity rather than admit an unbounded chart"
    );
    // Adaptive contract: near-surface cells are TIGHT; far cells may be
    // loose (that is the budget being spent where it matters).
    // The adaptive CONTRACT at the documented budget (2000 splits/cell):
    // near-band widths bounded and far cells strictly cheaper/looser.
    // Achieved widths are LEDGERED below; hull lower bounds near the
    // pole parameterization converge slowly, so the tight-claim budget
    // is a documented trade (2.6e-4 at 8000 splits, ~1e-3 at 2000).
    assert!(
        tile.worst_near_width < 5e-3,
        "near band bounded at the documented budget: {}",
        tile.worst_near_width
    );
    assert!(
        tile.worst_near_width < tile.worst_far_width,
        "effort concentrates near the surface: near {} vs far {}",
        tile.worst_near_width,
        tile.worst_far_width
    );
    // The center sample (inside) is negative; a corner is positive.
    let mid = tile.values[3 + 8 * 3 + 64 * 3]; // grid point (-0.2)^3, inside
    let corner = tile.values[0];
    assert!(corner > 0.0, "corner outside: {corner}");
    assert!(mid < 0.0, "center inside: {mid}");
    let mut far_surface = sphere();
    for row in &mut far_surface.cpw {
        for control in row {
            control[0] += 1.0e39 * control[3];
        }
    }
    let far_chart = ShellSdfChart::new(
        ShellSdf::new(vec![far_surface], vec![None], Orientation::Unknown)
            .expect("finite far shell"),
        1e-5,
        0,
        0.0,
    )
    .expect("admissible far-field chart settings");
    let small_box = fs_geom::Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0));
    let far_tile = generate_tile(&far_chart, &small_box, 2, 1e-5, 0)
        .expect("f64-finite distances outside f32 range downgrade");
    assert_eq!(far_tile.values, vec![f32::INFINITY; 8]);
    assert_eq!(far_tile.downgraded, 8);
    assert_eq!(far_tile.range_downgraded, 8);
    assert_eq!(far_tile.trim_downgraded, 0);
    assert_eq!(far_tile.sign_downgraded, 0);
    assert!(tile.worst_near_storage_error.is_finite());
    assert!(tile.worst_far_storage_error.is_finite());
    assert_eq!(far_tile.worst_near_storage_error, 0.0);
    assert_eq!(far_tile.worst_far_storage_error, f64::INFINITY);
    // Throughput evidence for the ledger line.
    println!(
        "{{\"metric\":\"nurbs-sdf-tile\",\"cells\":512,\"total_splits\":{},\
         \"near_width\":{:.2e},\"far_width\":{:.2e},\
         \"near_storage_error\":{:.2e},\"far_storage_error\":{:.2e}}}",
        tile.total_splits,
        tile.worst_near_width,
        tile.worst_far_width,
        tile.worst_near_storage_error,
        tile.worst_far_storage_error
    );
    verdict(
        "ns-005",
        "8^3 adaptive tile: near band tight, signs correct, splits ledgered",
    );
}
