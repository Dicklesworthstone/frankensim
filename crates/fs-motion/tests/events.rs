//! Validated-event battery (bead 6b8h). The Sev-0 rule under test:
//! no silent event misses, no fake finite certificates. Dense scans
//! FALSIFY the certified lane; they never prove it.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ga::Motor;
use fs_geom::Point3;
use fs_ivl::{Interval, TaylorModel1};
use fs_math::det;
use fs_motion::MotionError;
use fs_motion::analytic::{ScrewParams, screw_tube_with_derivative};
use fs_motion::events::{
    CertifiedEvent, CrossingDirection, EventScan, EventScanConfig, GuardFamily, GuardModel,
    PossibleReason, RootCountCertificate, ScanReceipt, ScanVerdict, enumerate_simultaneous,
    estimated_scan, plane_crossing_guard, scan_events,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x6B8,
                kernel_id: 11,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// A unit-rate rotation about z through the origin. Tests track the
/// point (1, 0, 0), so x(t) = cos t and y(t) = sin t analytically.
fn rotor_fixture() -> ScrewParams {
    ScrewParams {
        axis: [0.0, 0.0, 1.0],
        center: [0.0, 0.0, 0.0],
        omega: 1.0,
        axial_velocity: 0.0,
        base_pose: Motor::identity(),
    }
}

fn one_certified_event(window: Interval) -> EventScan {
    EventScan {
        certified: vec![CertifiedEvent {
            window,
            direction: CrossingDirection::Rising,
            guard_band: Interval::new(-1.0, 1.0),
        }],
        possible: Vec::new(),
        count: RootCountCertificate {
            confirmed: 1,
            possible_windows: 0,
            verdict: ScanVerdict::Complete,
        },
        receipt: ScanReceipt {
            intervals_examined: 1,
            max_depth: 0,
            excluded_leaves: 0,
            widest_leaf_guard_band: 2.0,
        },
    }
}

#[test]
fn ev_001_known_root_count_all_found_and_certified() {
    with_cx(|cx| {
        // Point at radius 1 rotating at ω = 1 for t ∈ [0, 6.2]:
        // x(t) = cos t. Guard x = 0.5 ⇒ cos t = 0.5 at t = π/3 (falling)
        // and t = 5π/3 (rising): exactly two roots in the span.
        let params = rotor_fixture();
        let (tube, rate) = screw_tube_with_derivative(&params, Interval::new(0.0, 6.2), 10, 8)
            .expect("tube builds");
        let family = plane_crossing_guard(
            &tube,
            &rate,
            Point3::new(1.0, 0.0, 0.0),
            [1.0, 0.0, 0.0],
            0.5,
        )
        .expect("guard builds");
        let scan = scan_events(
            &tube,
            &family,
            Interval::new(0.0, 6.2),
            &EventScanConfig::default(),
            cx,
        )
        .expect("scan completes");
        println!(
            "ev-001: confirmed {} possible {} verdict {:?} intervals {} depth {}",
            scan.count.confirmed,
            scan.count.possible_windows,
            scan.count.verdict,
            scan.receipt.intervals_examined,
            scan.receipt.max_depth
        );
        assert_eq!(scan.count.verdict, ScanVerdict::Complete);
        assert_eq!(scan.count.confirmed, 2, "exactly two crossings expected");
        assert_eq!(scan.count.possible_windows, 0);
        let t1 = std::f64::consts::FRAC_PI_3;
        let t2 = 5.0 * std::f64::consts::FRAC_PI_3;
        assert!(
            scan.certified[0].window.lo() <= t1 && t1 <= scan.certified[0].window.hi(),
            "first window {:?} must contain pi/3",
            scan.certified[0].window
        );
        assert_eq!(scan.certified[0].direction, CrossingDirection::Falling);
        assert!(
            scan.certified[1].window.lo() <= t2 && t2 <= scan.certified[1].window.hi(),
            "second window {:?} must contain 5pi/3",
            scan.certified[1].window
        );
        assert_eq!(scan.certified[1].direction, CrossingDirection::Rising);
    });
}

#[test]
fn ev_002_no_event_certificate_when_plane_out_of_reach() {
    with_cx(|cx| {
        let params = rotor_fixture();
        let (tube, rate) = screw_tube_with_derivative(&params, Interval::new(0.0, 6.2), 8, 6)
            .expect("tube builds");
        // The orbit has radius 1; the plane x = 1.5 is unreachable.
        let family = plane_crossing_guard(
            &tube,
            &rate,
            Point3::new(1.0, 0.0, 0.0),
            [1.0, 0.0, 0.0],
            1.5,
        )
        .expect("guard builds");
        let scan = scan_events(
            &tube,
            &family,
            Interval::new(0.0, 6.2),
            &EventScanConfig::default(),
            cx,
        )
        .expect("scan completes");
        println!(
            "ev-002: verdict {:?} excluded leaves {}",
            scan.count.verdict, scan.receipt.excluded_leaves
        );
        assert_eq!(scan.count.verdict, ScanVerdict::Complete);
        assert_eq!(scan.count.confirmed, 0);
        assert_eq!(scan.count.possible_windows, 0);
        assert!(scan.receipt.excluded_leaves >= 1);
    });
}

#[test]
fn ev_003_grazing_yields_unknown_not_a_fake_certificate() {
    with_cx(|cx| {
        // Tangency: plane exactly at the orbit radius. cos t = 1 at
        // t = 0 (endpoint) and the interior maximum at t = 2π grazes
        // the plane. A finite crossing certificate here would be the
        // Sev-0 failure; Unknown IS correct.
        let params = rotor_fixture();
        let (tube, rate) = screw_tube_with_derivative(&params, Interval::new(1.0, 7.0), 10, 8)
            .expect("tube builds");
        let family = plane_crossing_guard(
            &tube,
            &rate,
            Point3::new(1.0, 0.0, 0.0),
            [1.0, 0.0, 0.0],
            1.0,
        )
        .expect("guard builds");
        let scan = scan_events(
            &tube,
            &family,
            Interval::new(1.0, 7.0),
            &EventScanConfig {
                min_width: 1e-6,
                ..EventScanConfig::default()
            },
            cx,
        )
        .expect("scan completes");
        println!(
            "ev-003: confirmed {} possible {} verdict {:?}",
            scan.count.confirmed, scan.count.possible_windows, scan.count.verdict
        );
        assert_eq!(
            scan.count.confirmed, 0,
            "a tangency must never yield a certified crossing"
        );
        assert!(
            scan.count.possible_windows >= 1,
            "the grazing window must surface as a possible event"
        );
        assert!(
            scan.possible
                .iter()
                .any(|p| p.reason == PossibleReason::Grazing
                    && p.window.lo() <= 2.0 * std::f64::consts::PI
                    && 2.0 * std::f64::consts::PI <= p.window.hi()),
            "the grazing window must cover t = 2pi"
        );
        assert_eq!(scan.count.verdict, ScanVerdict::IncompleteUnknownWindows);
    });
}

#[test]
fn ev_004_dense_scan_falsifier_all_estimated_events_inside_windows() {
    with_cx(|cx| {
        // A screw with axial drift and an oblique plane: no analytic
        // count claimed; instead the 100x dense scan must never find a
        // sign change outside the certified/possible windows.
        let params = ScrewParams {
            axis: [0.0, 0.0, 1.0],
            center: [0.2, -0.1, 0.0],
            omega: 2.3,
            axial_velocity: 0.15,
            base_pose: Motor::rotor([1.0, 0.0, 0.0], 0.3),
        };
        let (tube, rate) = screw_tube_with_derivative(&params, Interval::new(0.0, 5.0), 10, 10)
            .expect("tube builds");
        let family = plane_crossing_guard(
            &tube,
            &rate,
            Point3::new(0.9, 0.1, -0.2),
            [0.3, -0.5, 0.8],
            0.25,
        )
        .expect("guard builds");
        let span = Interval::new(0.0, 5.0);
        let scan = scan_events(&tube, &family, span, &EventScanConfig::default(), cx)
            .expect("scan completes");
        let estimated = estimated_scan(&tube, &family, span, 4_000, cx).expect("dense scan");
        println!(
            "ev-004: certified {} possible {} estimated {}",
            scan.count.confirmed,
            scan.count.possible_windows,
            estimated.len()
        );
        for e in &estimated {
            let inside_certified = scan
                .certified
                .iter()
                .any(|c| c.window.lo() <= e.bracket.hi() && e.bracket.lo() <= c.window.hi());
            let inside_possible = scan
                .possible
                .iter()
                .any(|p| p.window.lo() <= e.bracket.hi() && e.bracket.lo() <= p.window.hi());
            assert!(
                inside_certified || inside_possible,
                "dense-scan sign change at {:?} escaped every window — Sev-0",
                e.bracket
            );
        }
        // Fixture-level cross-check only: this dense baseline happens to resolve
        // at least the certified crossings. It is not certificate evidence.
        assert!(scan.count.confirmed <= estimated.len() + scan.count.possible_windows);
    });
}

#[test]
fn ev_005_simultaneous_events_enumerate_admissible_orders() {
    with_cx(|cx| {
        // Two symmetric guards crossed by the same rotation within
        // overlapping windows: the order is genuinely undetermined at
        // window resolution, so both orders must be admissible.
        let params = rotor_fixture();
        let (tube, rate) = screw_tube_with_derivative(&params, Interval::new(0.0, 3.0), 10, 6)
            .expect("tube builds");
        let g1 = plane_crossing_guard(
            &tube,
            &rate,
            Point3::new(1.0, 0.0, 0.0),
            [1.0, 0.0, 0.0],
            0.0,
        )
        .expect("guard 1");
        // The same physical event observed through two guards — the
        // canonical simultaneous case: identical windows overlap by
        // construction, so the order is genuinely undetermined and
        // BOTH orders must be admissible. A third, disjoint crossing
        // must NOT join the group.
        let g2 = g1.clone();
        let g3 = plane_crossing_guard(
            &tube,
            &rate,
            Point3::new(1.0, 0.0, 0.0),
            [0.0, 1.0, 0.0],
            0.5,
        )
        .expect("guard 3");
        let config = EventScanConfig::default();
        let s1 = scan_events(&tube, &g1, Interval::new(0.0, 3.0), &config, cx).expect("scan 1");
        let s2 = scan_events(&tube, &g2, Interval::new(0.0, 3.0), &config, cx).expect("scan 2");
        let s3 = scan_events(&tube, &g3, Interval::new(0.0, 3.0), &config, cx).expect("scan 3");
        assert_eq!(s1.count.confirmed, 1, "x = 0 crosses once in [0, 3]");
        assert_eq!(
            s3.count.confirmed, 2,
            "y = 0.5 crosses twice in [0, 3] (pi/6 rising, 5pi/6 falling)"
        );
        let groups = enumerate_simultaneous(&[s1, s2, s3]);
        println!("ev-005: {} groups", groups.len());
        let multi: Vec<_> = groups.iter().filter(|g| g.members.len() > 1).collect();
        assert_eq!(multi.len(), 1, "exactly one simultaneous group expected");
        assert_eq!(multi[0].members.len(), 2, "the duplicated guard pair");
        let orders = multi[0]
            .admissible_orders
            .as_ref()
            .expect("small groups enumerate");
        assert_eq!(
            orders.len(),
            2,
            "both orders admissible — set-valued, not picked"
        );
        // The disjoint y-crossings (t ≈ π/6, 5π/6) stay their own groups.
        assert!(groups.iter().filter(|g| g.members.len() == 1).count() >= 2);
    });
}

#[test]
fn ev_006_budget_exhaustion_is_visible_never_silent() {
    with_cx(|cx| {
        let params = rotor_fixture();
        let (tube, rate) = screw_tube_with_derivative(&params, Interval::new(0.0, 6.2), 10, 8)
            .expect("tube builds");
        let family = plane_crossing_guard(
            &tube,
            &rate,
            Point3::new(1.0, 0.0, 0.0),
            [1.0, 0.0, 0.0],
            0.5,
        )
        .expect("guard builds");
        let scan = scan_events(
            &tube,
            &family,
            Interval::new(0.0, 6.2),
            &EventScanConfig {
                // Fewer than the tube's segment count: later segments
                // MUST drain into visible possible-event windows.
                max_subdivisions: 3,
                ..EventScanConfig::default()
            },
            cx,
        )
        .expect("scan completes structurally");
        println!(
            "ev-006: verdict {:?} confirmed {} possible {}",
            scan.count.verdict, scan.count.confirmed, scan.count.possible_windows
        );
        assert_eq!(scan.count.verdict, ScanVerdict::SubdivisionBudgetExhausted);
        assert!(
            scan.count.possible_windows >= 1,
            "unclassified remainders must surface as possible events"
        );
    });
}

#[test]
fn ev_007_scan_is_bitwise_deterministic() {
    with_cx(|cx| {
        let params = ScrewParams {
            axis: [0.0, 0.0, 1.0],
            center: [0.1, 0.2, 0.0],
            omega: 1.7,
            axial_velocity: -0.05,
            base_pose: Motor::translator(0.0, 0.3, 0.1),
        };
        let run = |cx: &Cx<'_>| {
            let (tube, rate) = screw_tube_with_derivative(&params, Interval::new(0.0, 4.0), 9, 7)
                .expect("tube builds");
            let family = plane_crossing_guard(
                &tube,
                &rate,
                Point3::new(0.8, -0.3, 0.0),
                [0.6, 0.8, 0.0],
                0.1,
            )
            .expect("guard builds");
            scan_events(
                &tube,
                &family,
                Interval::new(0.0, 4.0),
                &EventScanConfig::default(),
                cx,
            )
            .expect("scan completes")
        };
        let a = run(cx);
        let b = run(cx);
        assert_eq!(a.count.confirmed, b.count.confirmed);
        assert_eq!(a.count.possible_windows, b.count.possible_windows);
        assert_eq!(a.receipt.intervals_examined, b.receipt.intervals_examined);
        for (x, y) in a.certified.iter().zip(b.certified.iter()) {
            assert_eq!(x.window.lo().to_bits(), y.window.lo().to_bits());
            assert_eq!(x.window.hi().to_bits(), y.window.hi().to_bits());
        }
        println!(
            "ev-007: bitwise replay over {} intervals",
            a.receipt.intervals_examined
        );
    });
}

#[test]
fn ev_008_certified_windows_localize_the_analytic_roots_tightly() {
    with_cx(|cx| {
        // Localization quality gate: the certified windows for cos t =
        // 0.5 shrink to the resolution floor around the true roots and
        // the guard value at the true root is enclosed near zero.
        let params = rotor_fixture();
        let (tube, rate) = screw_tube_with_derivative(&params, Interval::new(0.0, 6.2), 10, 8)
            .expect("tube builds");
        let family = plane_crossing_guard(
            &tube,
            &rate,
            Point3::new(1.0, 0.0, 0.0),
            [1.0, 0.0, 0.0],
            0.5,
        )
        .expect("guard builds");
        let scan = scan_events(
            &tube,
            &family,
            Interval::new(0.0, 6.2),
            &EventScanConfig::default(),
            cx,
        )
        .expect("scan completes");
        assert_eq!(scan.count.confirmed, 2);
        for (window, truth) in scan.certified.iter().map(|c| c.window).zip([
            std::f64::consts::FRAC_PI_3,
            5.0 * std::f64::consts::FRAC_PI_3,
        ]) {
            let width = window.width();
            println!("ev-008: window width {width:e} around root {truth}");
            assert!(width < 1e-3, "window failed to localize: width {width:e}");
            // Cross-check with the independent pointwise value at the
            // analytic root: |cos(truth) − 0.5| ≈ 0 within rounding.
            let residual = (det::cos(truth) - 0.5).abs();
            assert!(residual < 1e-12);
        }
    });
}

#[test]
fn ev_009_guard_domains_must_match_tube_segments() {
    with_cx(|cx| {
        let (tube, _) = screw_tube_with_derivative(&rotor_fixture(), Interval::new(0.0, 1.0), 4, 1)
            .expect("tube builds");
        let wrong_domain = Interval::new(2.0, 3.0);
        let family = GuardFamily::new(vec![
            GuardModel::new(
                TaylorModel1::variable(wrong_domain, 1).expect("guard model"),
                TaylorModel1::constant(1.0, wrong_domain, 1).expect("derivative model"),
            )
            .expect("guard pair"),
        ])
        .expect("guard family");

        let certified_error = scan_events(
            &tube,
            &family,
            Interval::new(0.0, 1.0),
            &EventScanConfig::default(),
            cx,
        )
        .expect_err("misaligned certified scan must refuse");
        assert!(matches!(
            certified_error,
            MotionError::InvalidConfiguration {
                what: "each guard model domain must equal its tube segment domain"
            }
        ));

        let estimated_error = estimated_scan(&tube, &family, Interval::new(0.0, 1.0), 8, cx)
            .expect_err("misaligned estimated scan must refuse");
        assert!(matches!(
            estimated_error,
            MotionError::InvalidConfiguration {
                what: "each guard model domain must equal its tube segment domain"
            }
        ));
    });
}

#[test]
fn ev_010_receipt_includes_certified_leaf_guard_width() {
    with_cx(|cx| {
        let domain = Interval::new(0.0, 1.0);
        let (tube, _) =
            screw_tube_with_derivative(&rotor_fixture(), domain, 4, 1).expect("tube builds");
        let t = TaylorModel1::variable(domain, 1).expect("time model");
        let root = TaylorModel1::constant(0.5, domain, 1).expect("root model");
        let guard = t.try_sub(&root).expect("linear guard");
        let derivative = TaylorModel1::constant(1.0, domain, 1).expect("derivative");
        let family = GuardFamily::new(vec![
            GuardModel::new(guard, derivative).expect("guard pair"),
        ])
        .expect("guard family");

        let scan = scan_events(&tube, &family, domain, &EventScanConfig::default(), cx)
            .expect("scan completes");
        assert_eq!(scan.count.verdict, ScanVerdict::Complete);
        assert_eq!(scan.count.confirmed, 1);
        assert_eq!(scan.receipt.excluded_leaves, 0);
        assert!(
            scan.receipt.widest_leaf_guard_band >= 1.0,
            "the certified parent leaf spans guard values from -0.5 to 0.5"
        );
    });
}

#[test]
fn ev_011_zeno_budget_never_confirms_more_than_the_limit() {
    with_cx(|cx| {
        let params = rotor_fixture();
        let span = Interval::new(0.0, 6.2);
        let (tube, rate) = screw_tube_with_derivative(&params, span, 10, 8).expect("tube builds");
        let family = plane_crossing_guard(
            &tube,
            &rate,
            Point3::new(1.0, 0.0, 0.0),
            [1.0, 0.0, 0.0],
            0.5,
        )
        .expect("guard builds");
        let scan = scan_events(
            &tube,
            &family,
            span,
            &EventScanConfig {
                max_certified_events: 1,
                ..EventScanConfig::default()
            },
            cx,
        )
        .expect("scan returns a typed budget verdict");

        assert_eq!(scan.count.verdict, ScanVerdict::ZenoBudgetExceeded);
        assert_eq!(scan.count.confirmed, 1, "the configured cap is strict");
        assert!(
            !scan.possible.is_empty(),
            "the untouched tail stays visible"
        );
        let later_root = 5.0 * std::f64::consts::FRAC_PI_3;
        assert!(
            scan.possible
                .iter()
                .any(|event| event.window.lo() <= later_root && later_root <= event.window.hi()),
            "the unscanned second analytic crossing must remain inside a possible window"
        );
    });
}

#[test]
fn ev_012_chain_overlap_preserves_disjoint_pair_precedence() {
    let scans = [
        one_certified_event(Interval::new(0.0, 2.0)),
        one_certified_event(Interval::new(1.0, 3.0)),
        one_certified_event(Interval::new(2.5, 4.0)),
    ];
    let groups = enumerate_simultaneous(&scans);
    assert_eq!(groups.len(), 1, "overlap connectivity forms one group");
    let orders = groups[0]
        .admissible_orders
        .as_ref()
        .expect("three-member groups enumerate");
    assert_eq!(
        orders.len(),
        3,
        "only the three linear extensions are valid"
    );
    for order in orders {
        let first = order
            .iter()
            .position(|&member| member == 0)
            .expect("each permutation contains the first member");
        let third = order
            .iter()
            .position(|&member| member == 2)
            .expect("each permutation contains the third member");
        assert!(first < third, "disjoint [0, 2] must precede [2.5, 4]");
    }
}
