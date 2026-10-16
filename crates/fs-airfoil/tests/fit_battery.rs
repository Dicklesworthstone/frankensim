//! E4.0b fit battery (bead wf-root-guzez.5.1.2). Per-item oracles: basis
//! partition-of-unity and linear reproduction (Greville), synthetic fit
//! round-trip, fail-closed constraint verification with falsifiers,
//! regime-continuity twins, convention-block refusals, caps at cap AND
//! cap±1, and a pinned fit golden.
//! Repro: cargo test -p fs-airfoil --test fit_battery

use fs_airfoil::fit::{
    BsplineAxis, DiffConstraint, DiffDirection, FitSample, ResidualSurface,
    verify_regime_continuity,
};
use fs_airfoil::table::{
    ANGLES_ID, AXES_ID, CoefficientTable, ConventionBlock, MOMENT_SIGNS_ID, RegimePatch,
    SurfaceKind,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-airfoil-fit\",\"case\":\"{case}\",{payload}}}");
}

fn axis(name: &'static str, lo: f64, hi: f64, n: usize) -> BsplineAxis {
    BsplineAxis {
        name,
        lo,
        hi,
        n_coef: n,
    }
}

fn good_block() -> ConventionBlock {
    ConventionBlock {
        axes_id: AXES_ID.into(),
        moment_signs_id: MOMENT_SIGNS_ID.into(),
        angles_id: ANGLES_ID.into(),
        reference_area_convention: "per-surface projected area (geometry-conventions-v1)".into(),
        wind_reference: "not-applicable".into(),
    }
}

#[test]
fn basis_partitions_unity_and_reproduces_linears() {
    let ax = axis("alpha_rad", -0.3, 0.4, 7);
    // Partition of unity across the whole domain (per-point oracle).
    for i in 0..=200 {
        let x = -0.3 + 0.7 * f64::from(i) / 200.0;
        let (_, b) = {
            // basis() is private; probe through a surface with unit coefs.
            let s = ResidualSurface {
                axes: [
                    ax.clone(),
                    axis("log10_re", 5.0, 7.0, 1),
                    axis("delta_rad", 0.0, 1.0, 1),
                ],
                coef: vec![1.0; 7],
                constraints: vec![],
            };
            (0, [s.eval([x, 6.0, 0.5]), 0.0, 0.0, 0.0])
        };
        assert!(
            (b[0] - 1.0).abs() < 1e-12,
            "partition of unity broke at x={x}: {}",
            b[0]
        );
    }
    // Linear reproduction via Greville abscissae: coef_i = 3 + 2·g_i.
    let coef: Vec<f64> = (0..7).map(|i| 3.0 + 2.0 * ax.greville(i)).collect();
    let s = ResidualSurface {
        axes: [
            ax.clone(),
            axis("log10_re", 5.0, 7.0, 1),
            axis("delta_rad", 0.0, 1.0, 1),
        ],
        coef,
        constraints: vec![],
    };
    for i in 0..=50 {
        let x = -0.3 + 0.7 * f64::from(i) / 50.0;
        let want = 3.0 + 2.0 * x;
        let got = s.eval([x, 6.0, 0.0]);
        assert!(
            (got - want).abs() < 1e-12,
            "linear reproduction at x={x}: {got} vs {want}"
        );
    }
    jlog(
        "basis",
        "\"partition_of_unity\":true,\"linear_reproduction\":true",
    );
}

#[test]
fn fit_round_trips_synthetic_data() {
    // Ground truth: a smooth deterministic function over (α, logRe, δ).
    let truth = |x: [f64; 3]| 0.4 * x[0] - 0.05 * (x[1] - 6.0) + 0.2 * x[0] * x[2] + 0.1;
    let axes = [
        axis("alpha_rad", -0.2, 0.3, 6),
        axis("log10_re", 5.0, 7.0, 4),
        axis("delta_rad", -0.5, 0.5, 4),
    ];
    let mut samples = Vec::new();
    for i in 0..12 {
        for j in 0..8 {
            for k in 0..8 {
                let x = [
                    -0.2 + 0.5 * f64::from(i) / 11.0,
                    5.0 + 2.0 * f64::from(j) / 7.0,
                    -0.5 + 1.0 * f64::from(k) / 7.0,
                ];
                samples.push(FitSample { x, y: truth(x) });
            }
        }
    }
    let s = ResidualSurface::fit(axes, vec![], &samples).expect("fit must succeed");
    // Round-trip oracle: evaluation matches truth on an OFF-GRID probe set.
    let mut worst = 0.0f64;
    for i in 0..7 {
        for j in 0..5 {
            let x = [
                -0.17 + 0.44 * f64::from(i) / 6.0,
                5.3 + 1.5 * f64::from(j) / 4.0,
                -0.31 + 0.6 * f64::from(i) / 6.0,
            ];
            worst = worst.max((s.eval(x) - truth(x)).abs());
        }
    }
    assert!(worst < 1e-8, "off-grid round-trip error {worst:e}");
    jlog(
        "round-trip",
        &format!("\"worst_offgrid\":{worst:e},\"n_coef\":{}", s.n_total()),
    );
}

#[test]
fn constraints_fail_closed_with_falsifier() {
    let axes = || {
        [
            axis("alpha_rad", 0.0, 1.0, 5),
            axis("log10_re", 5.0, 7.0, 1),
            axis("delta_rad", 0.0, 1.0, 1),
        ]
    };
    let nondecreasing = vec![DiffConstraint {
        axis: 0,
        direction: DiffDirection::NonDecreasing,
    }];
    // Monotone truth fits fine under the constraint.
    let mono: Vec<FitSample> = (0..24)
        .map(|i| {
            let x = f64::from(i) / 23.0;
            FitSample {
                x: [x, 6.0, 0.5],
                y: 0.8 * x,
            }
        })
        .collect();
    let ok = ResidualSurface::fit(axes(), nondecreasing.clone(), &mono);
    assert!(ok.is_ok(), "monotone data must pass the constraint");
    // FALSIFIER: non-monotone truth under the same constraint must refuse
    // with the typed code — never a silently projected fit.
    let bump: Vec<FitSample> = (0..24)
        .map(|i| {
            let x = f64::from(i) / 23.0;
            FitSample {
                x: [x, 6.0, 0.5],
                y: -(x - 0.5) * (x - 0.5),
            }
        })
        .collect();
    let refusal = ResidualSurface::fit(axes(), nondecreasing, &bump).unwrap_err();
    assert_eq!(refusal.code, "fit-constraint-violated");
    assert!(
        refusal.message.contains("alpha_rad"),
        "must name the axis: {}",
        refusal.message
    );
    // The SAME data fits fine without the wrong constraint (the physics-is-
    // not-monotone rule: constraints are per-regime declarations, not law).
    assert!(ResidualSurface::fit(axes(), vec![], &bump).is_ok());
    jlog(
        "constraints",
        "\"fail_closed\":true,\"falsifier\":\"executed\"",
    );
}

#[test]
fn insufficient_samples_at_cap_and_cap_minus_one() {
    let axes = || {
        [
            axis("alpha_rad", 0.0, 1.0, 4),
            axis("log10_re", 5.0, 7.0, 1),
            axis("delta_rad", 0.0, 1.0, 1),
        ]
    };
    let make = |m: usize| -> Vec<FitSample> {
        (0..m)
            .map(|i| {
                let x = f64::from(i as u32) / (m.max(2) - 1) as f64;
                FitSample {
                    x: [x, 6.0, 0.0],
                    y: x,
                }
            })
            .collect()
    };
    // Exactly n samples: admitted (cap).
    assert!(ResidualSurface::fit(axes(), vec![], &make(4)).is_ok());
    // n − 1 samples: refused (cap − 1).
    assert_eq!(
        ResidualSurface::fit(axes(), vec![], &make(3))
            .unwrap_err()
            .code,
        "insufficient-samples"
    );
    // Degenerate axis count refusals: n_coef = 3 is invalid, 4 valid, 1 valid.
    let bad = BsplineAxis {
        name: "alpha_rad",
        lo: 0.0,
        hi: 1.0,
        n_coef: 3,
    };
    assert_eq!(bad.admit().unwrap_err().code, "axis-coef-count-invalid");
    jlog(
        "caps",
        "\"samples\":\"n ok, n-1 refuses\",\"n_coef\":\"3 refuses, 4 ok\"",
    );
}

#[test]
fn regime_continuity_twins() {
    // Two constant patches: equal → pass; offset face → typed refusal.
    let patch = |lo: f64, hi: f64, value: f64| ResidualSurface {
        axes: [
            axis("alpha_rad", lo, hi, 4),
            axis("log10_re", 5.0, 7.0, 1),
            axis("delta_rad", 0.0, 1.0, 1),
        ],
        coef: vec![value; 4],
        constraints: vec![],
    };
    let a = patch(0.0, 0.3, 0.7);
    let b_good = patch(0.3, 1.0, 0.7);
    let worst = verify_regime_continuity(&a, &b_good, 1e-9).expect("equal faces must pass");
    assert!(worst < 1e-12);
    let b_bad = patch(0.3, 1.0, 0.7005);
    let refusal = verify_regime_continuity(&a, &b_bad, 1e-9).unwrap_err();
    assert_eq!(refusal.code, "regime-boundary-discontinuity");
    // Non-abutting patches are their own refusal.
    let b_gap = patch(0.4, 1.0, 0.7);
    assert_eq!(
        verify_regime_continuity(&a, &b_gap, 1e-9).unwrap_err().code,
        "regime-boundary-mismatch"
    );
    jlog("continuity", &format!("\"worst_equal_faces\":{worst:e}"));
}

#[test]
fn convention_block_and_table_gates() {
    let patch = ResidualSurface {
        axes: [
            axis("alpha_rad", -0.2, 0.3, 4),
            axis("log10_re", 5.0, 7.0, 1),
            axis("delta_rad", 0.0, 1.0, 1),
        ],
        coef: vec![0.01; 4],
        constraints: vec![],
    };
    let table = CoefficientTable {
        kind: SurfaceKind::Wing,
        channel: "cl-residual",
        dossier_record: "a2-ames-fullscale-1999".into(),
        conventions: good_block(),
        patches: vec![RegimePatch {
            regime: "attached",
            surface: patch.clone(),
        }],
    };
    assert!(table.validate(1e-9).is_ok());
    assert!((table.eval([0.1, 6.0, 0.2]).unwrap() - 0.01).abs() < 1e-12);
    // Fitted-domain refusal is DISTINCT from the global admission domain.
    assert_eq!(
        table.eval([0.9, 6.0, 0.2]).unwrap_err().code,
        "alpha-outside-table"
    );
    // Convention falsifiers: wrong axes id and empty wind reference.
    let mut wrong = table.clone();
    wrong.conventions.axes_id = "z-up-graphics".into();
    let refusal = wrong.validate(1e-9).unwrap_err();
    assert_eq!(refusal.code, "convention-block-mismatch");
    assert!(refusal.message.contains("axes_id"));
    let mut missing = table.clone();
    missing.conventions.wind_reference.clear();
    assert_eq!(
        missing.validate(1e-9).unwrap_err().code,
        "convention-block-missing"
    );
    // Provenance gate.
    let mut anon = table.clone();
    anon.dossier_record.clear();
    assert_eq!(anon.validate(1e-9).unwrap_err().code, "provenance-missing");
    jlog("table-gates", "\"convention_falsifiers\":\"executed\"");
}

#[test]
fn fit_golden_digest() {
    // Deterministic golden over the fitted coefficients of a fixed synthetic
    // problem (golden-bump protocol; measure-then-pin).
    let axes = [
        axis("alpha_rad", -0.2, 0.3, 5),
        axis("log10_re", 5.0, 7.0, 4),
        axis("delta_rad", 0.0, 1.0, 1),
    ];
    let samples: Vec<FitSample> = (0..120)
        .map(|i| {
            let t = f64::from(i) / 119.0;
            let x = [
                -0.2 + 0.5 * t,
                5.0 + 2.0 * f64::from((i * 7) % 120) / 119.0,
                0.5,
            ];
            FitSample {
                x,
                y: 0.3 * x[0] + 0.02 * (x[1] - 6.0),
            }
        })
        .collect();
    let s = ResidualSurface::fit(axes, vec![], &samples).unwrap();
    let mut payload = Vec::new();
    for c in &s.coef {
        payload.extend_from_slice(&c.to_bits().to_le_bytes());
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-airfoil.fit-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "f25f5e763c6322369cc0d72227b1807aad900e1b5af9025ca70cedfbe5ddc3cc",
        "fit golden moved — determinism regression or an intentional fit change \
         requiring the golden-bump protocol"
    );
}
