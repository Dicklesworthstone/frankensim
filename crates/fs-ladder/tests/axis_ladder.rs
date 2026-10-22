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

// ---- M5.5-W2 simplification-equivalence comparison ----

const Z_UNIT: &str = "Pa s/m^3";

fn solution(tag: u8, values: &[(f64, f64)], bounds: &[f64]) -> BoundedSolution {
    BoundedSolution {
        description: fs_blake3::ContentHash([tag; 32]),
        ladder: fs_blake3::ContentHash([tag + 1; 32]),
        unit: Z_UNIT.into(),
        values: values.to_vec(),
        bounds: bounds.to_vec(),
        sources: vec![BoundSource::Declared; values.len()],
    }
}
/// Nine points at x = 0..8 on one line, sampled on nested rungs
/// {0,4,8} < {0,2,4,6,8} < all.
fn domain9() -> EquivalenceDomain {
    EquivalenceDomain {
        description: "synthetic: nine points, three nested rungs".into(),
        rungs: vec![vec![0, 4, 8], vec![0, 2, 4, 6, 8], (0..9).collect()],
        lines: vec![(0..9).collect()],
        coordinate: (0..9).map(f64::from).collect(),
        ensemble: EnsembleDecl::Deterministic {
            reason: "passive, no stochastic source".into(),
        },
        limitations: vec![],
    }
}
/// Simplified side at 1 + 0i everywhere; the detailed side differs by `d[i]`.
fn pair(d: &[f64; 9], b: f64) -> (BoundedSolution, BoundedSolution) {
    let detailed: Vec<_> = d.iter().map(|x| (1.0 + x, 0.0)).collect();
    (
        solution(1, &detailed, &[b; 9]),
        solution(3, &[(1.0, 0.0); 9], &[b; 9]),
    )
}

#[test]
fn equivalence_verdict_three_way() {
    // |d| = 1 - (x - 3.3)^2 / 100 on x = 0..8: supremum 1 between samples.
    // Every rung's stencil sees u'' = -0.02, so the curvature terms are
    // 0.02 H^2 / 8 = 0.04, 0.01, 0.0025 for H = 4, 2, 1, and the sampled
    // maximum 0.9991 plus 0.0025 covers the supremum.
    let parabola: [f64; 9] =
        core::array::from_fn(|i| 1.0 - (i as f64 - 3.3) * (i as f64 - 3.3) / 100.0);
    let (a, b) = pair(&parabola, 0.0);
    let report = compare_equivalence(&a, &b, &domain9(), Metric::Absolute, 1.01).unwrap();
    for (got, want) in report.sampling_terms.iter().zip([0.04, 0.01, 0.0025]) {
        assert!((got - want).abs() < 1e-12, "{:?}", report.sampling_terms);
    }
    assert!((report.sampling_term - 0.0025).abs() < 1e-12);
    let EquivalenceVerdict::Equivalent { max_upper } = report.verdict else {
        panic!("{:?}", report.verdict);
    };
    assert!(max_upper >= 1.0 && (max_upper - 1.0016).abs() < 1e-12);
    assert!(report.authorizes_omission());
    // The sampled lower bound 0.9991 exceeds a 0.999 tolerance at x = 3.
    let report = compare_equivalence(&a, &b, &domain9(), Metric::Absolute, 0.999).unwrap();
    let EquivalenceVerdict::NotEquivalent { max_lower, at } = report.verdict else {
        panic!("{:?}", report.verdict);
    };
    assert!(at == 3 && (max_lower - 0.9991).abs() < 1e-12);
    assert!(!report.authorizes_omission());
    // Tolerance 1: lower 0.9991 <= 1 < 0.9991 + 0.0025 straddles at x = 3 only.
    let report = compare_equivalence(&a, &b, &domain9(), Metric::Absolute, 1.0).unwrap();
    assert_eq!(
        report.verdict,
        EquivalenceVerdict::Inconclusive {
            straddling: vec![3]
        }
    );
    // A spike only the finest rung sees makes the term grow (0, 0, 0.25):
    // an unresolved feature, so the sampling term is infinite.
    let mut spike = [0.0; 9];
    spike[5] = 1.0;
    let (a, b) = pair(&spike, 0.0);
    let report = compare_equivalence(&a, &b, &domain9(), Metric::Absolute, 10.0).unwrap();
    assert_eq!(report.sampling_terms, vec![0.0, 0.0, 0.25]);
    assert!(report.sampling_term.is_infinite());
    assert!(matches!(
        report.verdict,
        EquivalenceVerdict::Inconclusive { .. }
    ));
    // Level metric: 10x the modulus is 20 dB; a modulus interval reaching
    // zero has no finite upper bound.
    let mut tenfold = [0.0; 9];
    tenfold[0] = 9.0;
    let (a, b) = pair(&tenfold, 0.0);
    let report = compare_equivalence(&a, &b, &domain9(), Metric::LevelDb, 25.0).unwrap();
    assert!((report.upper[0] - 20.0).abs() < 1e-12 && (report.lower[0] - 20.0).abs() < 1e-12);
    assert!(report.upper[1] == 0.0 && report.lower[1] == 0.0);
    let (a, b) = pair(&tenfold, 1.5);
    let report = compare_equivalence(&a, &b, &domain9(), Metric::LevelDb, 25.0).unwrap();
    assert!(report.upper[1].is_infinite() && report.lower[1] == 0.0);
    // |a| in [8.5, 11.5], |b| in [0, 2.5]: at least 20 log10(8.5/2.5).
    let floor = 20.0 * (8.5_f64 / 2.5).log10();
    assert!(
        (report.lower[0] - floor).abs() < 1e-12,
        "{}",
        report.lower[0]
    );
}

#[test]
fn equivalence_requires_both_bounds() {
    let (a, b) = pair(&[0.0; 9], 0.0);
    let domain = domain9();
    assert!(compare_equivalence(&a, &b, &domain, Metric::Absolute, 1.0).is_ok());
    let mut short = b.clone();
    short.bounds.pop();
    assert!(compare_equivalence(&a, &short, &domain, Metric::Absolute, 1.0).is_err());
    let mut negative = b.clone();
    negative.bounds[3] = -1e-9;
    assert!(compare_equivalence(&a, &negative, &domain, Metric::Absolute, 1.0).is_err());
    let mut missing = b.clone();
    missing.bounds[3] = f64::NAN;
    assert!(compare_equivalence(&missing, &a, &domain, Metric::Absolute, 1.0).is_err());
    let mut other_unit = b.clone();
    other_unit.unit = "Pa/(m^3/s)".into();
    assert!(compare_equivalence(&a, &other_unit, &domain, Metric::Absolute, 1.0).is_err());
    let mut unsourced = b.clone();
    unsourced.sources.pop();
    assert!(compare_equivalence(&a, &unsourced, &domain, Metric::Absolute, 1.0).is_err());
    // Oscillating rungs have no GCI: their spread (2 - 1) is the bound, and
    // the point is marked as spread-bounded.
    let oscillating = ladder(&[(1.0, 2.0), (0.5, 1.0), (0.25, 1.5), (0.125, 1.2)]);
    let points = [("re".to_string(), "im".to_string())];
    let hash = fs_blake3::ContentHash([7; 32]);
    let spread = BoundedSolution::from_ladder(hash, &oscillating, Z_UNIT, &points).unwrap();
    assert_eq!(spread.sources, vec![BoundSource::Spread]);
    assert!((spread.bounds[0] - 1.0).abs() < 1e-12);
    // An observable a non-algebraic axis does not bound (ObservedOnly)
    // refuses when the solution is built, before any comparison.
    let refused = BoundedSolution::from_ladder(hash, &spectral_ladder(), Z_UNIT, &points);
    assert!(refused.unwrap_err().contains("im has no ladder bound"));
    let wrong = [("re".to_string(), "absent".to_string())];
    let manufactured = ladder(&[(1.0, 5.0), (0.5, 3.125), (0.25, 3.0 + 2.0 / 256.0)]);
    assert!(BoundedSolution::from_ladder(hash, &manufactured, Z_UNIT, &wrong).is_err());
    assert!(BoundedSolution::from_ladder(hash, &manufactured, "Pa", &points).is_err());
}

#[test]
fn domain_sampling_cap_three_rungs() {
    let (a, b) = pair(&[0.0; 9], 0.0);
    let mut domain = domain9();
    assert!(compare_equivalence(&a, &b, &domain, Metric::Absolute, 1.0).is_ok());
    domain.rungs.remove(0);
    assert!(compare_equivalence(&a, &b, &domain, Metric::Absolute, 1.0).is_err());
    let mut not_nested = domain9();
    not_nested.rungs[1] = vec![1, 2, 4, 6, 8];
    assert!(compare_equivalence(&a, &b, &not_nested, Metric::Absolute, 1.0).is_err());
    let mut omits = domain9();
    omits.rungs[2].pop();
    assert!(compare_equivalence(&a, &b, &omits, Metric::Absolute, 1.0).is_err());
    let refuses = |edit: &dyn Fn(&mut EquivalenceDomain)| {
        let mut bad = domain9();
        edit(&mut bad);
        compare_equivalence(&a, &b, &bad, Metric::Absolute, 1.0).is_err()
    };
    assert!(refuses(&|d| d.coordinate.truncate(8)));
    assert!(refuses(&|d| d.coordinate[4] = f64::NAN));
    assert!(refuses(&|d| d.coordinate.reverse()));
    assert!(refuses(&|d| d.lines = vec![(0..8).collect()]));
    assert!(refuses(&|d| d.lines.push(vec![3])));
    // The coarsest rung keeps one point of the first line: no stencil.
    assert!(refuses(
        &|d| d.lines = vec![vec![0, 1, 2, 3], (4..9).collect()]
    ));
    let mut seeded = domain9();
    seeded.ensemble = EnsembleDecl::Seeds {
        seeds: vec![1, 2],
        uncertainty: "bootstrap".into(),
    };
    assert!(compare_equivalence(&a, &b, &seeded, Metric::Absolute, 1.0).is_err());
    for tolerance in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(compare_equivalence(&a, &b, &domain9(), Metric::Absolute, tolerance).is_err());
    }
}

/// A one-axis report with observables `re` (the given rungs, `h` coarse to
/// fine) and `im` (flat at -4).
fn ladder(rungs: &[(f64, f64)]) -> LadderReport {
    let records: Vec<RungRecord> = rungs
        .iter()
        .map(|&(h, y)| {
            let mut r = rung(h, y);
            r.qoi = vec![
                Observable {
                    name: "re".into(),
                    value: y,
                    unit: Z_UNIT.into(),
                },
                Observable {
                    name: "im".into(),
                    value: -4.0,
                    unit: Z_UNIT.into(),
                },
            ];
            r
        })
        .collect();
    let report = evaluate_axis(
        &axis(ErrorModel::Algebraic { formal_order: 4.0 }),
        &records,
        &sensitivity(),
        &policy(),
    );
    LadderReport::new(
        policy(),
        vec![report],
        vec![],
        records,
        EnsembleDecl::Deterministic {
            reason: "no stochastic source".into(),
        },
    )
    .unwrap()
}

/// A spectral-tail axis that declares a bound for `re` only.
fn spectral_ladder() -> LadderReport {
    let mut spectral = axis(ErrorModel::SpectralTail);
    spectral.bound_targets = vec!["re".into()];
    let records: Vec<RungRecord> = [(1.0, 3.5), (0.5, 3.25), (0.25, 3.125)]
        .iter()
        .map(|&(h, y)| {
            let mut r = rung(h, y);
            let observable = |name: &str, value| Observable {
                name: name.into(),
                value,
                unit: Z_UNIT.into(),
            };
            r.qoi = vec![observable("re", y), observable("im", -4.0)];
            r.bounds = vec![observable("re", h)];
            r
        })
        .collect();
    let report = evaluate_axis(&spectral, &records, &sensitivity(), &policy());
    LadderReport::new(
        policy(),
        vec![report],
        vec![],
        records,
        EnsembleDecl::Deterministic {
            reason: "no stochastic source".into(),
        },
    )
    .unwrap()
}

#[test]
fn ladder_bounds_follow_axis_verdicts() {
    // y = 3 + 2 h^4: the finest rung's true error is 2 h^4, and the GCI bound
    // (Fs = 1.25 on an exact power law) is 1.25 times it. The flat imaginary
    // part is Unresolved and bounded at the flatness threshold.
    let rungs: Vec<(f64, f64)> = [1.0, 0.5, 0.25, 0.125]
        .iter()
        .map(|&h| (h, 3.0 + 2.0 * fs_math::det::pow(h, 4.0)))
        .collect();
    let report = ladder(&rungs);
    let (re, re_source) = qoi_bound(&report, "re").unwrap();
    let truth = 2.0 * fs_math::det::pow(0.125, 4.0);
    assert!((re / truth - 1.25).abs() < 1e-3, "{re} vs {truth}");
    assert_eq!(re_source, BoundSource::Gci);
    let (im, im_source) = qoi_bound(&report, "im").unwrap();
    assert!(im == 128.0 * f64::EPSILON * 4.0, "{im}");
    assert_eq!(im_source, BoundSource::Roundoff);
    assert!(qoi_bound(&report, "absent").is_none());
    // A declared spectral bound is the finest rung's: 0.25.
    let declared = qoi_bound(&spectral_ladder(), "re").unwrap();
    assert_eq!(declared, (0.25, BoundSource::Declared));
    let points = [("re".to_string(), "im".to_string())];
    let hash = fs_blake3::ContentHash([9; 32]);
    let s = BoundedSolution::from_ladder(hash, &report, Z_UNIT, &points).unwrap();
    assert_eq!(s.values, vec![(rungs[3].1, -4.0)]);
    assert!((s.bounds[0] - (re * re + im * im).sqrt()).abs() <= 1e-15 * s.bounds[0]);
    assert_eq!(s.ladder, report.hash);
    assert_eq!(s.sources, vec![BoundSource::Gci]);
}

#[test]
fn equivalence_report_identity_and_limitations() {
    let (a, b) = pair(&[0.01; 9], 0.0);
    let mut domain = domain9();
    domain.limitations = vec!["1-D chimney: evanescent near field omitted".into()];
    let one = compare_equivalence(&a, &b, &domain, Metric::Absolute, 1.0).unwrap();
    let two = compare_equivalence(&a, &b, &domain, Metric::Absolute, 1.0).unwrap();
    assert_eq!(one.hash, two.hash);
    assert!(matches!(one.verdict, EquivalenceVerdict::Equivalent { .. }));
    assert!(
        !one.authorizes_omission(),
        "a limitation withholds authorization"
    );
    assert_eq!(one.limitations, domain.limitations);
    let three = compare_equivalence(&a, &b, &domain, Metric::Absolute, 0.9).unwrap();
    assert_ne!(one.hash, three.hash);
    assert_eq!(one.descriptions, [a.description, b.description]);
    assert_eq!(one.ladders, [a.ladder, b.ladder]);
}
