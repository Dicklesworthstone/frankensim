use super::*;
use fs_topopt::sdf3::stress::{StressError3, StressMeasure3, StressOptions3};

#[test]
fn thermal_reference_stress_adjoints_keep_relaxation_and_physical_scale_separate() {
    let op = adaptive();
    let field = op.prepare_eigenstrain(&free, 100_000, proceed).unwrap();
    let second = op
        .prepare_eigenstrain(
            &|cell, p| free(cell, p).map(|e| e.map(|v| -0.4 * v)),
            100_000,
            proceed,
        )
        .unwrap();
    let body = op
        .prepare_material_body_load(&|_| [0.0, -0.0002, 0.0], 100_000, proceed)
        .unwrap();
    let force = vec![0.0; op.n()];
    let cases = [
        CutEigenstrainCase3 {
            force: &force,
            weight: 0.2,
            eigenstrain: &field,
        },
        CutEigenstrainCase3 {
            force: &force,
            weight: 0.8,
            eigenstrain: &second,
        },
    ];
    let mut study =
        CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op, 100_000_000), 0.12, params())
            .with_material_body_load(body)
            .unwrap();
    let rho: Vec<_> = (0..study.cells())
        .map(|i| 0.45 + 0.01 * (i % 9) as f64)
        .collect();
    let incoming = study.operator().scales().to_vec();
    let mut callback = |_| proceed();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    for measure in [
        StressMeasure3::NormalizedAverage,
        StressMeasure3::SampledPeakBound,
    ] {
        let base = study
            .evaluate_eigenstrain_stress(
                &rho,
                &cases,
                policy(),
                Default::default(),
                measure,
                &mut control,
            )
            .unwrap();
        assert!(base.aggregate > 0.0);
        assert_eq!(study.operator().scales(), incoming);
        if measure == StressMeasure3::SampledPeakBound {
            assert!(base.aggregate >= base.sampled_physical_max);
            assert!(base.aggregate >= base.sampled_relaxed_max);
        }
        for direction in [
            vec![1.0; rho.len()],
            (0..rho.len()).map(|i| (i % 3) as f64 - 1.0).collect(),
        ] {
            let h = 1e-4;
            let a: Vec<_> = rho.iter().zip(&direction).map(|(r, d)| r + h * d).collect();
            let b: Vec<_> = rho.iter().zip(&direction).map(|(r, d)| r - h * d).collect();
            let plus = study
                .evaluate_eigenstrain_stress(
                    &a,
                    &cases,
                    policy(),
                    Default::default(),
                    measure,
                    &mut control,
                )
                .unwrap();
            let minus = study
                .evaluate_eigenstrain_stress(
                    &b,
                    &cases,
                    policy(),
                    Default::default(),
                    measure,
                    &mut control,
                )
                .unwrap();
            close(
                (plus.aggregate - minus.aggregate) / (2.0 * h),
                dot(&base.gradient, &direction),
                5e-4,
            );
        }
    }
}

#[test]
fn zero_field_preserves_cold_stress_and_total_point_budget_is_transactional() {
    let op = raw([2, 1, 2], 0.1, false);
    let field = op
        .prepare_eigenstrain(&|_, _| Ok([0.0; 6]), 100_000, proceed)
        .unwrap();
    let force = op.body_load(&|_| [0.0, -0.0001, 0.0], proceed).unwrap();
    let cases = [CutEigenstrainCase3 {
        force: &force,
        weight: 1.0,
        eigenstrain: &field,
    }];
    let loads = [LoadCase {
        force: &force,
        weight: 1.0,
    }];
    let mut study = CutDensityStudy3::new(op, 0.1, params());
    let rho = vec![0.55; study.cells()];
    let incoming = study.operator().scales().to_vec();
    let mut callback = |_| proceed();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let cold = study
        .evaluate_stress(&rho, &loads, Default::default(), &mut control)
        .unwrap();
    let thermal = study
        .evaluate_eigenstrain_stress(
            &rho,
            &cases,
            policy(),
            Default::default(),
            StressMeasure3::NormalizedAverage,
            &mut control,
        )
        .unwrap();
    close(cold.aggregate, thermal.aggregate, 1e-12);
    for (a, b) in cold.gradient.iter().zip(&thermal.gradient) {
        close(*a, *b, 1e-11);
    }
    assert_eq!(cold.displacements, thermal.displacements);
    let small = StressOptions3 {
        max_points: field.bulk_points(),
        ..Default::default()
    };
    assert!(matches!(
        study.evaluate_eigenstrain_stress(
            &rho,
            &[cases[0], cases[0]],
            policy(),
            small,
            StressMeasure3::NormalizedAverage,
            &mut control
        ),
        Err(StressError3::PointBudget)
    ));
    assert_eq!(study.operator().scales(), incoming);
}

#[test]
fn zero_weight_thermal_case_cannot_hide_elastic_domain_exit() {
    let op = raw([2, 1, 2], 0.01, true);
    let good = op.prepare_eigenstrain(&free, 100_000, proceed).unwrap();
    let bad = op
        .prepare_eigenstrain(&|_, _| Ok([0.2, 0.2, 0.2, 0.0, 0.0, 0.0]), 100_000, proceed)
        .unwrap();
    let force = vec![0.0; op.n()];
    let cases = [
        CutEigenstrainCase3 {
            force: &force,
            weight: 1.0,
            eigenstrain: &good,
        },
        CutEigenstrainCase3 {
            force: &force,
            weight: 0.0,
            eigenstrain: &bad,
        },
    ];
    let mut study = CutDensityStudy3::new(op, 0.1, params());
    let incoming = study.operator().scales().to_vec();
    let mut callback = |_| proceed();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    assert!(matches!(
        study.evaluate_eigenstrain_stress(
            &vec![0.6; study.cells()],
            &cases,
            policy(),
            Default::default(),
            StressMeasure3::NormalizedAverage,
            &mut control
        ),
        Err(StressError3::MaterialDomain(_))
    ));
    assert_eq!(study.operator().scales(), incoming);
}

#[test]
fn cancellation_inside_native_thermal_stress_preserves_incoming_scales() {
    let op = raw([2, 1, 2], 0.1, true);
    let field = op.prepare_eigenstrain(&free, 100_000, proceed).unwrap();
    let force = vec![0.0; op.n()];
    let cases = [CutEigenstrainCase3 {
        force: &force,
        weight: 1.0,
        eigenstrain: &field,
    }];
    let mut study = CutDensityStudy3::new(op, 0.1, params());
    let incoming = study.operator().scales().to_vec();
    let mut polls = 0;
    let mut callback = |progress: SolveProgress| {
        if progress.stage == "sdf3-eigenstrain-stress" {
            polls += 1;
        }
        if polls == 3 {
            ControlFlow::Break(())
        } else {
            proceed()
        }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    assert!(matches!(
        study.evaluate_eigenstrain_stress(
            &vec![0.6; study.cells()],
            &cases,
            policy(),
            Default::default(),
            StressMeasure3::SampledPeakBound,
            &mut control
        ),
        Err(StressError3::Evaluation(EvaluationStop::Cancelled))
    ));
    assert_eq!(study.operator().scales(), incoming);
}
