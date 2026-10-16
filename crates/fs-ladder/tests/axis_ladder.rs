//! G1 manufactured powers/ASME GCI, G3 liveness and G5 report identity.
use fs_ladder::axis::*;

fn policy() -> AnalysisPolicy {
    AnalysisPolicy {
        maximum_order: 16.0,
        order_tolerance: 0.02,
        relative_fit_tolerance: 0.001,
        asymptotic_tolerance: 0.02,
        minimum_refinement_ratio: 1.1,
    }
}
fn axis(model: ErrorModel) -> RefinableAxis {
    RefinableAxis {
        name: "h".into(),
        unit: "m".into(),
        model,
        group: None,
        coordinate: Coordinate::Direct,
        bound_targets: vec!["energy".into()],
    }
}
fn observable(value: f64) -> Observable {
    Observable {
        name: "energy".into(),
        value,
        unit: "J".into(),
    }
}
fn rung(h: f64, value: f64) -> RungRecord {
    RungRecord {
        requested: h,
        dimensions: vec![Dimension {
            name: "h".into(),
            value: h,
            unit: "m".into(),
        }],
        qoi: vec![observable(value)],
        bounds: vec![],
        regime: "interior".into(),
        solver_residual: 0.0,
        wall_s: 0.0,
    }
}
fn sensitivity() -> SensitivityFixture {
    SensitivityFixture {
        fixture: "G1: integral of a linear interpolant of x^2 on [0,1]".into(),
        realized: [1.0, 0.5],
        qoi: [observable(0.5), observable(0.375)],
        minimum_change: 1e-12,
    }
}
fn study(values: &[(f64, f64)], p: f64) -> AxisReport {
    evaluate_axis(
        &axis(ErrorModel::Algebraic { formal_order: p }),
        &values.iter().map(|(h, y)| rung(*h, *y)).collect::<Vec<_>>(),
        &sensitivity(),
        &policy(),
    )
}

#[test]
fn observed_order_recovers_manufactured_p() {
    for p in [1.0, 2.0, 4.0] {
        for hs in [
            vec![1.0, 0.5, 0.25, 0.125, 0.0625],
            vec![1.0, 0.5, 1.0 / 3.0, 1.0 / 6.0, 1.0 / 9.0],
        ] {
            let values: Vec<_> = hs
                .iter()
                .map(|h| (*h, 3.0 + 2.0 * fs_math::det::pow(*h, p)))
                .collect();
            let report = study(&values, p);
            let fit = report.qois[0].fit.as_ref().unwrap();
            assert_eq!(fit.points, hs.len());
            assert!((fit.observed_order.unwrap() - p).abs() <= 0.02, "{fit:?}");
            assert!(
                matches!(report.qois[0].verdict, AxisVerdict::Asymptotic { .. }),
                "{report:?}"
            );
        }
    }
}

#[test]
fn gci_matches_asme_hand_value() {
    // f=340+6250h^2; finest f=340.625, delta=1.875, r=2.
    // Fs=1.25, so GCI=1.25*(1.875/340.625)/(4-1).
    let report = study(&[(0.04, 350.0), (0.02, 342.5), (0.01, 340.625)], 2.0);
    let AxisVerdict::Asymptotic { gci: Some(gci), .. } = report.qois[0].verdict else {
        panic!("{report:?}");
    };
    let hand = 1.25 * (1.875 / 340.625) / 3.0;
    assert!(
        (gci / hand - 1.0).abs() <= 1e-12,
        "{gci:.17e} vs {hand:.17e}"
    );
    // Nonuniform finest ratio is 1.5, not an average of 2 and 1.5.
    let report = study(
        &[(1.0, 12.0), (0.5, 10.5), (1.0 / 3.0, 10.0 + 2.0 / 9.0)],
        2.0,
    );
    let AxisVerdict::Asymptotic { gci: Some(gci), .. } = report.qois[0].verdict else {
        panic!("{report:?}");
    };
    let hand = 1.25 * ((10.5 - (10.0 + 2.0 / 9.0)) / (10.0 + 2.0 / 9.0)) / (2.25 - 1.0);
    assert!(
        (gci / hand - 1.0).abs() <= 1e-12,
        "{gci:.17e} vs {hand:.17e}"
    );
}

#[test]
fn gci_refused_when_assumptions_fail() {
    let oscillatory = study(&[(1.0, 2.0), (0.5, 0.5), (0.25, 1.25), (0.125, 0.875)], 2.0);
    assert_eq!(oscillatory.qois[0].verdict, AxisVerdict::Oscillatory);
    let wrong_order = study(&[(1.0, 2.0), (0.5, 1.5), (0.25, 1.25), (0.125, 1.125)], 2.0);
    assert!(matches!(
        wrong_order.qois[0].verdict,
        AxisVerdict::PreAsymptotic { .. }
    ));
    let flat = study(&[(1.0, 1.0), (0.5, 1.0), (0.25, 1.0)], 2.0);
    assert_eq!(flat.qois[0].verdict, AxisVerdict::Unresolved);
    let zero_fine = study(&[(1.0, 15.0), (0.5, 3.0), (0.25, 0.0)], 2.0);
    assert!(matches!(
        zero_fine.qois[0].verdict,
        AxisVerdict::Asymptotic { gci: None, .. }
    ));
}

#[test]
fn eca_hoekstra_scattered_classification() {
    // The three finest rungs alone are a perfect p=2 sequence. Earlier
    // observations scatter; all-rung fitting must not discard that evidence.
    let report = study(
        &[
            (2.0, 1.8),
            (1.0, 3.0),
            (0.5, 1.25),
            (0.25, 1.0625),
            (0.125, 1.015625),
        ],
        2.0,
    );
    assert_eq!(report.qois[0].verdict, AxisVerdict::Scattered);
    assert_eq!(report.qois[0].fit.as_ref().unwrap().points, 5);
    // Monotonicity alone is insufficient: a poor power-law fit also scatters.
    let report = study(
        &[
            (1.0, 8.0),
            (0.5, 5.0),
            (0.25, 4.9),
            (0.125, 1.0),
            (0.0625, 0.8),
        ],
        2.0,
    );
    assert_eq!(report.qois[0].verdict, AxisVerdict::Scattered);
}

#[test]
fn two_rungs_refuse_three_accept() {
    let two = study(&[(1.0, 2.0), (0.5, 1.25)], 2.0);
    assert_eq!(two.qois[0].verdict, AxisVerdict::InsufficientRungs);
    assert!(!two.complete());
    assert!(study(&[(1.0, 2.0), (0.5, 1.25), (0.25, 1.0625)], 2.0).complete());
}

#[test]
fn axis_liveness_detects_clamped_axis() {
    let mut a = axis(ErrorModel::SpectralTail);
    a.name = "n_fft".into();
    a.unit = "samples".into();
    a.coordinate = Coordinate::Reciprocal;
    let rs: Vec<_> = [2048, 4096, 8192, 16384]
        .into_iter()
        .map(|request| {
            let actual = request.min(8192) as f64;
            RungRecord {
                requested: f64::from(request),
                dimensions: vec![Dimension {
                    name: a.name.clone(),
                    unit: a.unit.clone(),
                    value: actual,
                }],
                qoi: vec![observable(1.0 + 1.0 / actual)],
                bounds: vec![observable(1.0 / actual)],
                regime: "interior".into(),
                solver_residual: 0.0,
                wall_s: 0.0,
            }
        })
        .collect();
    let mut s = sensitivity();
    s.realized = [2048.0, 4096.0];
    let report = evaluate_axis(&a, &rs, &s, &policy());
    assert_eq!(report.dead_rungs, vec![3]);
    assert_eq!(report.live_rungs, 3);
    let report = evaluate_axis(&a, &rs[2..], &s, &policy());
    assert_eq!(report.live_rungs, 1);
    assert!(!report.complete());
}

#[test]
fn spectral_axis_reports_declared_bound() {
    for model in [
        ErrorModel::SpectralTail,
        ErrorModel::OperatorError,
        ErrorModel::ResidualBound,
        ErrorModel::Statistical,
    ] {
        let rs: Vec<_> = [1.0, 0.5, 0.25]
            .into_iter()
            .map(|h| {
                let mut r = rung(h, 2.0 + h);
                r.bounds = vec![observable(h)];
                r
            })
            .collect();
        let report = evaluate_axis(&axis(model), &rs, &sensitivity(), &policy());
        assert_eq!(report.qois[0].verdict, AxisVerdict::Bounded { bound: 0.25 });
        assert!(report.qois[0].fit.is_none());
        let mut missing = rs;
        missing[0].bounds.clear();
        assert!(
            !evaluate_axis(
                &axis(ErrorModel::SpectralTail),
                &missing,
                &sensitivity(),
                &policy()
            )
            .complete()
        );
    }
}

fn coupled_fixture() -> (Vec<RefinableAxis>, CoupledGroup, Vec<RungRecord>) {
    let mut a = axis(ErrorModel::Algebraic { formal_order: 2.0 });
    a.group = Some("duration".into());
    let mut b = a.clone();
    b.name = "n".into();
    b.unit = "samples".into();
    b.coordinate = Coordinate::Reciprocal;
    let group = CoupledGroup {
        name: "duration".into(),
        members: vec![("h".into(), 1.0), ("n".into(), 1.0)],
        invariant: 8.0,
        relative_tolerance: 1e-12,
    };
    let rs = [1.0, 0.5, 0.25]
        .into_iter()
        .map(|h| {
            let mut r = rung(h, 1.0 + h * h);
            r.dimensions.push(Dimension {
                name: "n".into(),
                unit: "samples".into(),
                value: 8.0 / h,
            });
            r
        })
        .collect();
    (vec![a, b], group, rs)
}

#[test]
fn coupled_group_refines_together_or_refuses() {
    let (axes, g, rs) = coupled_fixture();
    assert!(validate_group(&g, &axes, &rs).is_ok());
    let mut dead = rs.clone();
    dead[1].dimensions[1].value = 8.0;
    assert!(validate_group(&g, &axes, &dead).is_err());
    let mut moving = rs.clone();
    moving[1].dimensions[1].value = 17.0;
    assert!(validate_group(&g, &axes, &moving).is_err());
    assert!(validate_group(&g, &axes, &rs[..2]).is_err());
}

fn report() -> LadderReport {
    let (axes, g, rs) = coupled_fixture();
    let reports = axes
        .iter()
        .map(|a| evaluate_axis(a, &rs, &sensitivity(), &policy()))
        .collect();
    LadderReport::new(
        policy(),
        reports,
        vec![g],
        rs,
        EnsembleDecl::Deterministic {
            reason: "no stochastic source".into(),
        },
    )
    .unwrap()
}

#[test]
fn report_hash_is_deterministic() {
    let a = report();
    let b = report();
    assert_eq!(a.hash, b.hash);
    assert_eq!(a.canonical_bytes(), b.canonical_bytes());
    assert_eq!(a.joint.interactions.len(), 3);
    assert_eq!(a.joint.interactions[0].interaction, -(1.0 - 0.0625));
    let mut altered = a.clone();
    altered.axes[0].rungs[0].requested = 2.0;
    assert_ne!(a.canonical_bytes(), altered.canonical_bytes());
}

#[test]
fn invalid_data_regimes_and_insensitive_fixtures_refuse() {
    let rs = vec![rung(1.0, 2.0), rung(0.5, 1.25), rung(0.25, 1.0625)];
    let a = axis(ErrorModel::Algebraic { formal_order: 2.0 });
    let mut bad = rs.clone();
    bad[1].regime = "other register".into();
    let r = evaluate_axis(&a, &bad, &sensitivity(), &policy());
    assert_eq!(r.failed_rungs, vec![1]);
    assert!(!r.complete());
    let mut bad = rs.clone();
    bad[1].qoi[0].value = f64::NAN;
    assert!(!evaluate_axis(&a, &bad, &sensitivity(), &policy()).complete());
    let mut s = sensitivity();
    s.qoi[1].value = s.qoi[0].value;
    assert!(!evaluate_axis(&a, &rs, &s, &policy()).complete());
    assert!(!evaluate_axis(&axis(ErrorModel::Unknown), &rs, &sensitivity(), &policy()).complete());
    let mut bad = rs;
    bad[2].solver_residual = f64::NAN;
    assert!(!evaluate_axis(&a, &bad, &sensitivity(), &policy()).complete());
}

#[test]
fn legacy_evaluator_uses_all_rungs_without_scattered_gci() {
    let mut p = fs_ladder::ConvergencePlan::new("energy", 2.0);
    for (i, (h, q)) in [
        (2.0, 1.8),
        (1.0, 3.0),
        (0.5, 1.25),
        (0.25, 1.0625),
        (0.125, 1.015625),
    ]
    .into_iter()
    .enumerate()
    {
        p = p.with_rung(fs_ladder::MeshRung::new(
            i,
            format!("mesh-{i}"),
            h,
            "m",
            10,
            "energy",
            q,
            "J",
        ));
    }
    let r = fs_ladder::ConvergenceEvaluator::evaluate(&p);
    assert_eq!(r.admitted_rungs.len(), 5);
    assert!(r.discretization_error_gci.is_none());
    assert!(r.richardson_extrapolated_qoi.is_none());
}
