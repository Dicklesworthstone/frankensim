//! End-to-end battery: a deterministic PDHG cantilever iterate with an
//! advisory, endpoint-checked tropical load path from load to support.

use fs_evidence::Color;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_sparse::Csr;
use fs_truss::{
    LayoutCertificateProblem, LayoutCertificateRefusal, LayoutCertificateStatus, PdhgSettings,
};
use fs_truss_e2e::{
    LoadPathCertificateRefusal, LoadPathCertificateStatus, TrussError, analyze_load_path,
    certify_load_path, load_path_color_from_certificate, optimality_color_from_certificate,
    rescale_optimality_color, run_campaign,
};

fn with_gate_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x7A55,
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

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_gate_cx(&CancelGate::new(), f)
}

fn campaign(
    nx: usize,
    ny: usize,
    w: f64,
    h: f64,
    gap_tol: f64,
) -> Result<fs_truss_e2e::TrussReport, TrussError> {
    with_cx(|cx| run_campaign(nx, ny, w, h, gap_tol, cx))
}

#[test]
fn the_converged_truss_has_a_bounded_unique_load_path() {
    let report = campaign(4, 3, 4.0, 2.0, 1e-4).expect("valid tropical load path");
    // a real ground structure was optimized down to a sparse active set.
    assert!(report.num_members > report.num_active, "nothing was pruned");
    assert!(report.num_active > 0, "no active bars");
    assert!(report.total_volume > 0.0, "zero volume");
    // Solver convergence remains a diagnostic, independent of the outward
    // optimum certificate.
    assert!(
        report.solver_converged,
        "gap {} eq_res {}",
        report.gap, report.eq_residual
    );
    let Color::Verified { lo, hi } = report.optimality_color else {
        panic!("the repaired primal and checked dual must certify optimality");
    };
    assert!(lo.is_finite() && hi.is_finite());
    assert!(
        lo > 0.0,
        "the scaled solver dual must retain a useful bound"
    );
    assert!(lo <= hi, "inverted optimum interval [{lo}, {hi}]");
    // The advisory path is non-trivial and carries real rounded volume.
    assert!(
        report.critical_path.len() >= 2,
        "path too short: {:?}",
        report.critical_path
    );
    assert!(report.critical_path_volume > 0.0);
    assert!(report.bottleneck_member.is_some());
    assert!(
        report
            .critical_path
            .contains(&report.bottleneck_member.unwrap())
    );
    match (&report.load_path_status, &report.load_path_color) {
        (
            LoadPathCertificateStatus::Certified(path_certificate),
            Color::Verified {
                lo: path_lo,
                hi: path_hi,
            },
        ) => {
            assert!(path_lo.is_finite() && path_hi.is_finite() && *path_lo > 0.0);
            assert!(*path_lo <= report.critical_path_volume);
            assert!(report.critical_path_volume <= *path_hi);
            assert_eq!(path_certificate.analysis().members, report.critical_path);
            assert_eq!(
                path_certificate.analysis().bottleneck_member,
                report.bottleneck_member
            );
            assert_ne!(path_certificate.replay_golden(), 0);
        }
        (LoadPathCertificateStatus::Unavailable(_), Color::Estimated { dispersion, .. }) => {
            assert!(dispersion.is_infinite())
        }
        other => panic!("load-path color/status mismatch: {other:?}"),
    }
    // The selected path is a subset of the exactly feasible repaired design;
    // compare it to the certified whole-structure upper bound, not the nearby
    // rounded PDHG diagnostic.
    assert!(report.critical_path_volume <= hi + 1e-6);
    println!(
        "{{\"campaign\":\"trusspath\",\"members\":{},\"active\":{},\"volume\":{:.4},\"gap\":{:.2e},\
         \"eq_res\":{:.2e},\"iters\":{},\"path_len\":{},\"path_volume\":{:.4},\"bottleneck\":{:?}}}",
        report.num_members,
        report.num_active,
        report.total_volume,
        report.gap,
        report.eq_residual,
        report.iters,
        report.critical_path.len(),
        report.critical_path_volume,
        report.bottleneck_member,
    );
}

#[test]
fn the_campaign_is_deterministic() {
    let a = campaign(4, 3, 4.0, 2.0, 1e-4).expect("first run");
    let b = campaign(4, 3, 4.0, 2.0, 1e-4).expect("second run");
    assert_eq!(a.total_volume.to_bits(), b.total_volume.to_bits());
    assert_eq!(a.critical_path, b.critical_path);
    assert_eq!(a.bottleneck_member, b.bottleneck_member);
    assert_eq!(a.optimality_color, b.optimality_color);
    assert_eq!(a.load_path_color, b.load_path_color);
    assert_eq!(a.load_path_status, b.load_path_status);
}

#[test]
fn unavailable_certificate_never_promotes_finite_diagnostics() {
    let matrix = Csr::from_parts(1, 2, vec![0, 2], vec![0, 1], vec![1.0, -1.0]);
    let costs = [1.0, 1.0];
    let loads = [1.0];
    let problem = LayoutCertificateProblem::try_new(&matrix, &costs, &loads)
        .expect("well-formed paired fixture");
    let status = LayoutCertificateStatus::Unavailable(LayoutCertificateRefusal::RankDeficient {
        active_rows: 1,
        rank: 0,
    });
    let settings = PdhgSettings::default();
    with_cx(|cx| {
        assert!(matches!(
            optimality_color_from_certificate(
                &problem,
                &[0.0, 0.0],
                &[0.0],
                settings,
                &status,
                0.0,
                0.0,
                cx,
            )
            .expect("unavailable promotion fallback"),
            Color::Estimated {
                dispersion: 0.0,
                ..
            }
        ));
        assert!(matches!(
            optimality_color_from_certificate(
                &problem,
                &[0.0, 0.0],
                &[0.0],
                settings,
                &status,
                f64::NAN,
                0.0,
                cx,
            )
            .expect("non-finite diagnostic fallback"),
            Color::Estimated { dispersion, .. } if dispersion.is_infinite()
        ));
    });
}

#[test]
fn certificate_promotion_rejects_another_problem_and_rescales_outward() {
    let matrix = Csr::from_parts(1, 2, vec![0, 2], vec![0, 1], vec![1.0, -1.0]);
    let costs = [1.0, 1.0];
    let loads = [1.0];
    let other_loads = [2.0];
    let problem = LayoutCertificateProblem::try_new(&matrix, &costs, &loads)
        .expect("well-formed source problem");
    let other_problem = LayoutCertificateProblem::try_new(&matrix, &costs, &other_loads)
        .expect("well-formed distinct problem");
    let settings = PdhgSettings::default();
    with_cx(|cx| {
        let status = problem
            .certify_optimum(
                &[0.0, 0.0],
                &[0.0],
                settings,
                fs_truss::LayoutCertificateLimits::default(),
                cx,
            )
            .expect("source certificate attempt");
        assert!(matches!(status, LayoutCertificateStatus::Certified(_)));
        assert!(matches!(
            optimality_color_from_certificate(
                &other_problem,
                &[0.0, 0.0],
                &[0.0],
                settings,
                &status,
                0.0,
                0.0,
                cx,
            )
            .expect("context-mismatch preflight"),
            Color::Estimated { .. }
        ));
    });

    let scaled = rescale_optimality_color(&Color::Verified { lo: 1.0, hi: 2.0 }, 3.0);
    let Color::Verified { lo, hi } = scaled else {
        panic!("positive physical scaling must preserve Verified");
    };
    assert!(lo <= 1.0 / 3.0 && hi >= 2.0 / 3.0);
    assert!(matches!(
        rescale_optimality_color(&Color::Verified { lo: 1.0, hi: 2.0 }, 0.0),
        Color::Estimated { dispersion, .. } if dispersion.is_infinite()
    ));
}

#[test]
fn separated_member_intervals_mint_an_exact_path_receipt() {
    // Physical B = [[1, 0], [-1, 1]], followed by its negation for the
    // tension/compression split. q = [1, 1] is exactly equilibrated.
    let matrix = Csr::from_parts(
        2,
        4,
        vec![0, 2, 6],
        vec![0, 2, 0, 1, 2, 3],
        vec![1.0, -1.0, -1.0, 1.0, 1.0, -1.0],
    );
    let costs = [2.0, 1.0, 2.0, 1.0];
    let loads = [1.0, 0.0];
    let x = [1.0, 1.0, 0.0, 0.0];
    let y = [0.0, 0.0];
    let nodes = [[2.0, 0.0], [1.0, 0.0], [0.0, 0.0]];
    let members = [(0, 1), (1, 2)];
    let settings = PdhgSettings::default();
    let problem = LayoutCertificateProblem::try_new(&matrix, &costs, &loads)
        .expect("well-formed two-member fixture");

    with_cx(|cx| {
        let optimum = problem
            .certify_optimum(
                &x,
                &y,
                settings,
                fs_truss::LayoutCertificateLimits::default(),
                cx,
            )
            .expect("bounded optimum proof");
        let status = certify_load_path(
            &problem,
            &x,
            &y,
            settings,
            &optimum,
            &nodes,
            &members,
            0,
            &[2],
            cx,
        )
        .expect("bounded path proof");
        let LoadPathCertificateStatus::Certified(certificate) = &status else {
            panic!("strictly separated chain must certify: {status:?}");
        };
        assert_eq!(certificate.analysis().members, vec![0, 1]);
        assert_eq!(certificate.analysis().bottleneck_member, Some(0));
        assert!(certificate.path_weight_bounds().contains(3.0));
        assert!(certificate.member_weight_bounds()[0].contains(2.0));
        assert!(certificate.member_weight_bounds()[1].contains(1.0));
        assert!(certificate.active_threshold().hi() < 1.0);
        assert!(matches!(
            load_path_color_from_certificate(&status),
            Color::Verified { lo, hi } if lo <= 3.0 && hi >= 3.0
        ));
        assert!(
            certificate
                .verifies_for(
                    &problem,
                    &x,
                    &y,
                    settings,
                    &optimum,
                    &nodes,
                    &members,
                    0,
                    &[2],
                    cx,
                )
                .expect("exact receipt replay")
        );
        assert!(
            !certificate
                .verifies_for(
                    &problem,
                    &x,
                    &y,
                    settings,
                    &optimum,
                    &nodes,
                    &members,
                    0,
                    &[1],
                    cx,
                )
                .expect("altered endpoint must fail closed")
        );
    });
}

#[test]
fn tied_bottlenecks_and_direct_one_bar_paths_remain_estimated() {
    let matrix = Csr::from_parts(
        2,
        4,
        vec![0, 2, 6],
        vec![0, 2, 0, 1, 2, 3],
        vec![1.0, -1.0, -1.0, 1.0, 1.0, -1.0],
    );
    let costs = [1.0, 1.0, 1.0, 1.0];
    let loads = [1.0, 0.0];
    let x = [1.0, 1.0, 0.0, 0.0];
    let y = [0.0, 0.0];
    let settings = PdhgSettings::default();
    let problem = LayoutCertificateProblem::try_new(&matrix, &costs, &loads)
        .expect("well-formed tied fixture");
    with_cx(|cx| {
        let optimum = problem
            .certify_optimum(
                &x,
                &y,
                settings,
                fs_truss::LayoutCertificateLimits::default(),
                cx,
            )
            .expect("bounded optimum proof");
        let tied = certify_load_path(
            &problem,
            &x,
            &y,
            settings,
            &optimum,
            &[[2.0, 0.0], [1.0, 0.0], [0.0, 0.0]],
            &[(0, 1), (1, 2)],
            0,
            &[2],
            cx,
        )
        .expect("tied path proof attempt");
        assert!(matches!(
            tied,
            LoadPathCertificateStatus::Unavailable(
                LoadPathCertificateRefusal::BottleneckUnseparated
            )
        ));

        let direct_matrix = Csr::from_parts(1, 2, vec![0, 2], vec![0, 1], vec![1.0, -1.0]);
        let direct_costs = [1.0, 1.0];
        let direct_loads = [1.0];
        let direct_problem =
            LayoutCertificateProblem::try_new(&direct_matrix, &direct_costs, &direct_loads)
                .expect("well-formed direct fixture");
        let direct_status = direct_problem
            .certify_optimum(
                &[1.0, 0.0],
                &[0.0],
                settings,
                fs_truss::LayoutCertificateLimits::default(),
                cx,
            )
            .expect("direct optimum proof");
        let direct = certify_load_path(
            &direct_problem,
            &[1.0, 0.0],
            &[0.0],
            settings,
            &direct_status,
            &[[1.0, 0.0], [0.0, 0.0]],
            &[(0, 1)],
            0,
            &[1],
            cx,
        )
        .expect("direct path proof attempt");
        assert!(matches!(
            direct,
            LoadPathCertificateStatus::Unavailable(LoadPathCertificateRefusal::NoCompleteLoadPath)
        ));
    });
}

#[test]
fn equal_parallel_path_intervals_do_not_mint_a_unique_witness() {
    let matrix = Csr::from_parts(
        2,
        8,
        vec![0, 2, 6],
        vec![0, 4, 0, 1, 4, 5],
        vec![1.0, -1.0, -1.0, 1.0, 1.0, -1.0],
    );
    // Members 2 and 3 are equilibrium-neutral in this synthetic falsifier,
    // but their retained positive forces form the second geometric route.
    let costs = [2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0];
    let loads = [1.0, 0.0];
    let x = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    let y = [0.0, 0.0];
    let settings = PdhgSettings::default();
    let problem = LayoutCertificateProblem::try_new(&matrix, &costs, &loads)
        .expect("well-formed parallel-path fixture");
    let nodes = [[2.0, 0.0], [1.0, 1.0], [1.0, -1.0], [0.0, 0.0]];
    let members = [(0, 1), (1, 3), (0, 2), (2, 3)];
    with_cx(|cx| {
        let optimum = problem
            .certify_optimum(
                &x,
                &y,
                settings,
                fs_truss::LayoutCertificateLimits::default(),
                cx,
            )
            .expect("parallel-path optimum proof");
        let status = certify_load_path(
            &problem,
            &x,
            &y,
            settings,
            &optimum,
            &nodes,
            &members,
            0,
            &[3],
            cx,
        )
        .expect("parallel-path proof attempt");
        assert!(matches!(
            status,
            LoadPathCertificateStatus::Unavailable(
                LoadPathCertificateRefusal::CriticalPathUnseparated
            )
        ));
    });
}

#[test]
fn near_threshold_members_and_equal_distance_orientation_fail_closed() {
    // The third physical column is structurally zero, so its exact retained
    // force can probe threshold/orientation admission without changing the
    // independently equilibrated two-member chain.
    let matrix = Csr::from_parts(
        2,
        6,
        vec![0, 2, 6],
        vec![0, 3, 0, 1, 3, 4],
        vec![1.0, -1.0, -1.0, 1.0, 1.0, -1.0],
    );
    let costs = [2.0, 1.0, 1.0, 2.0, 1.0, 1.0];
    let loads = [1.0, 0.0];
    let y = [0.0, 0.0];
    let settings = PdhgSettings::default();
    let problem = LayoutCertificateProblem::try_new(&matrix, &costs, &loads)
        .expect("well-formed threshold fixture");
    let nodes = [[2.0, 0.0], [1.0, 0.0], [0.0, 0.0], [0.0, 1.0]];
    let members = [(0, 1), (1, 2), (1, 3)];
    with_cx(|cx| {
        let near_x = [1.0, 1.0, 1e-3, 0.0, 0.0, 0.0];
        let near_optimum = problem
            .certify_optimum(
                &near_x,
                &y,
                settings,
                fs_truss::LayoutCertificateLimits::default(),
                cx,
            )
            .expect("near-threshold optimum proof");
        let near = certify_load_path(
            &problem,
            &near_x,
            &y,
            settings,
            &near_optimum,
            &nodes,
            &members,
            0,
            &[2],
            cx,
        )
        .expect("near-threshold proof attempt");
        assert!(matches!(
            near,
            LoadPathCertificateStatus::Unavailable(
                LoadPathCertificateRefusal::ActiveSetUnseparated { member: 2 }
            )
        ));

        let equal_x = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
        let equal_optimum = problem
            .certify_optimum(
                &equal_x,
                &y,
                settings,
                fs_truss::LayoutCertificateLimits::default(),
                cx,
            )
            .expect("equal-distance optimum proof");
        let equal = certify_load_path(
            &problem,
            &equal_x,
            &y,
            settings,
            &equal_optimum,
            &nodes,
            &members,
            0,
            &[2],
            cx,
        )
        .expect("equal-distance proof attempt");
        assert!(matches!(
            equal,
            LoadPathCertificateStatus::Unavailable(
                LoadPathCertificateRefusal::OrientationUnseparated { member: 2 }
            )
        ));
    });
}

#[test]
fn cancelled_path_replay_publishes_no_certificate() {
    let matrix = Csr::from_parts(
        2,
        4,
        vec![0, 2, 6],
        vec![0, 2, 0, 1, 2, 3],
        vec![1.0, -1.0, -1.0, 1.0, 1.0, -1.0],
    );
    let costs = [2.0, 1.0, 2.0, 1.0];
    let loads = [1.0, 0.0];
    let x = [1.0, 1.0, 0.0, 0.0];
    let y = [0.0, 0.0];
    let settings = PdhgSettings::default();
    let problem = LayoutCertificateProblem::try_new(&matrix, &costs, &loads)
        .expect("well-formed cancellation fixture");
    let optimum = with_cx(|cx| {
        problem
            .certify_optimum(
                &x,
                &y,
                settings,
                fs_truss::LayoutCertificateLimits::default(),
                cx,
            )
            .expect("source optimum proof")
    });
    let gate = CancelGate::new();
    gate.request();
    let result = with_gate_cx(&gate, |cx| {
        certify_load_path(
            &problem,
            &x,
            &y,
            settings,
            &optimum,
            &[[2.0, 0.0], [1.0, 0.0], [0.0, 0.0]],
            &[(0, 1), (1, 2)],
            0,
            &[2],
            cx,
        )
    });
    assert!(matches!(
        result,
        Err(TrussError::Certificate(
            fs_truss::LayoutCertificateError::Cancelled { .. }
        ))
    ));
}

#[test]
fn invalid_or_unbounded_campaigns_refuse_before_ground_structure_work() {
    assert!(matches!(
        campaign(1, 2, 1.0, 1.0, 1e-4),
        Err(TrussError::InvalidInput {
            field: "grid dimensions",
            ..
        })
    ));
    assert!(matches!(
        campaign(17, 16, 1.0, 1.0, 1e-4),
        Err(TrussError::InvalidInput {
            field: "grid node count",
            ..
        })
    ));
    for (width, height, tolerance) in [
        (f64::NAN, 1.0, 1e-4),
        (1.0, f64::INFINITY, 1e-4),
        (1.0, 1.0, 0.0),
    ] {
        assert!(matches!(
            campaign(2, 2, width, height, tolerance),
            Err(TrussError::InvalidInput { .. })
        ));
    }
    assert!(matches!(
        campaign(2, 2, 0.01, 0.01, 1e-4),
        Err(TrussError::NoCandidateMembers)
    ));

    // 64 nodes is the exact cubic-preflight boundary. It reaches the later
    // candidate/solver budget; 65 nodes is refused before construction.
    let boundary = campaign(8, 8, 4.0, 2.0, 1e-4);
    assert!(
        matches!(boundary, Err(TrussError::WorkBudget { resource, .. }) if resource != "ground-structure triplet checks")
    );
    assert!(matches!(
        campaign(13, 5, 4.0, 2.0, 1e-4),
        Err(TrussError::WorkBudget {
            resource: "ground-structure triplet checks",
            ..
        })
    ));
}

#[test]
fn pre_cancelled_campaign_refuses_without_a_partial_report() {
    let gate = CancelGate::new();
    gate.request();
    let result = with_gate_cx(&gate, |cx| run_campaign(4, 3, 4.0, 2.0, 1e-4, cx));
    assert!(matches!(
        result,
        Err(TrussError::Construction(
            fs_truss::TrussConstructionError::Cancelled { .. }
        ))
    ));
}

#[test]
fn tight_tolerance_does_not_mislabel_the_iteration_cap() {
    let report = campaign(4, 3, 4.0, 2.0, f64::MIN_POSITIVE)
        .expect("bounded campaign still returns its final iterate");
    assert_eq!(report.iters, 60_000);
    assert!(!report.solver_converged);
}

#[test]
fn support_selection_is_index_based_even_below_the_old_coordinate_tolerance() {
    match campaign(2, 4, 1e-10, 1.0, 1e-4) {
        Ok(report) => {
            assert!(report.total_volume > 0.0);
            assert!(report.critical_path.len() >= 2);
        }
        Err(TrussError::NoCandidateMembers | TrussError::NoCompleteLoadPath) => {}
        Err(error) => panic!("unexpected narrow-grid refusal: {error}"),
    }
}

#[test]
fn path_analysis_excludes_disconnected_heavy_components_and_checks_endpoints() {
    let nodes = [[3.0, 0.0], [2.0, 0.0], [0.0, 0.0], [2.0, 2.0], [1.0, 2.0]];
    let members = [(0, 1), (1, 2), (3, 4)];
    let path = analyze_load_path(&nodes, &members, &[0, 1, 2], &[1.0, 2.0, 100.0], 0, &[2])
        .expect("the connected load-support chain survives filtering");
    assert_eq!(path.members, vec![0, 1]);
    assert_eq!(path.weight.to_bits(), 3.0_f64.to_bits());
    assert!(!path.members.contains(&2));

    assert!(matches!(
        analyze_load_path(&nodes, &members, &[0], &[1.0, 2.0, 100.0], 0, &[1]),
        Err(TrussError::NoCompleteLoadPath)
    ));
    assert!(matches!(
        analyze_load_path(
            &nodes,
            &members,
            &[0, 1],
            &[1.0, 2.0, 100.0],
            0,
            &[1, 2, 3, 4, 4, 4]
        ),
        Err(TrussError::InvalidLoadPath {
            reason: "support count must be within 1..=node count"
        })
    ));
}
