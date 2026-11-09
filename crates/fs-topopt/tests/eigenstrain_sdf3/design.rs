use super::*;
use fs_ascent::projected_al::ProjectedAlOptions;
use fs_topopt::sdf3::robust::{RobustDesignOptions3, RobustDesignStudy3, RobustOptions3};

#[test]
fn thermal_work_and_dilated_volume_follow_the_same_robust_density_maps() {
    let op = raw([2, 1, 2], 0.1, true);
    let field = op.prepare_eigenstrain(&free, 100_000, proceed).unwrap();
    let force = vec![0.0; op.n()];
    let cases = [CutEigenstrainCase3 {
        force: &force,
        weight: 1.0,
        eigenstrain: &field,
    }];
    let mut study = CutDensityStudy3::new(op, 0.12, params());
    let incoming = study.operator().scales().to_vec();
    let rho = vec![0.55; study.cells()];
    let options = RobustOptions3 {
        threshold_offset: 0.12,
    };
    let mut callback = |_| proceed();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let base = study
        .evaluate_eigenstrain_robust(&rho, &cases, policy(), options, &mut control)
        .unwrap();
    assert_eq!(study.operator().scales(), incoming);
    assert!(base.eroded.volume_fraction < base.nominal_volume_fraction);
    assert!(base.nominal_volume_fraction < base.dilated_volume_fraction);
    for i in 0..rho.len() {
        let h = 1e-4;
        let mut a = rho.clone();
        let mut b = rho.clone();
        a[i] += h;
        b[i] -= h;
        let plus = study
            .evaluate_eigenstrain_robust(&a, &cases, policy(), options, &mut control)
            .unwrap();
        let minus = study
            .evaluate_eigenstrain_robust(&b, &cases, policy(), options, &mut control)
            .unwrap();
        close(
            (plus.eroded.objective.compliance - minus.eroded.objective.compliance) / (2.0 * h),
            base.eroded.objective.gradient[i],
            4e-4,
        );
        close(
            (plus.dilated_volume_fraction - minus.dilated_volume_fraction) / (2.0 * h),
            base.dilated_volume_gradient[i],
            5e-5,
        );
    }
    assert_eq!(study.operator().scales(), incoming);
}

#[test]
fn projected_al_accepts_positive_thermal_work_gradients_and_retains_exact_fields() {
    let op = raw([1; 3], 1e-6, false);
    let field = op
        .prepare_eigenstrain(
            &|_, _| Ok([0.003, 0.0, 0.0, 0.0, 0.0, 0.0]),
            100_000,
            proceed,
        )
        .unwrap();
    let force = vec![0.0; op.n()];
    let cases = [CutEigenstrainCase3 {
        force: &force,
        weight: 1.0,
        eigenstrain: &field,
    }];
    let mut study = CutDensityStudy3::new(op, 0.0, params());
    let incoming = study.operator().scales().to_vec();
    let mut callback = |_| proceed();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let robustness = RobustOptions3 {
        threshold_offset: 0.05,
    };
    let baseline = study
        .evaluate_eigenstrain_robust(&[0.6], &cases, policy(), robustness, &mut control)
        .unwrap();
    assert!(baseline.eroded.objective.gradient[0] > 0.0);
    let options = RobustDesignOptions3 {
        robustness,
        volume_limit: 0.9,
        density_floor: 0.1,
        compliance_scale: baseline.eroded.objective.compliance,
        optimizer: ProjectedAlOptions {
            tolerance: 1e-6,
            max_evaluations: 200,
            ..Default::default()
        },
    };
    let accepted = {
        let mut session = RobustDesignStudy3::new_with_eigenstrain(
            &mut study,
            &cases,
            policy(),
            &[0.6],
            options,
            &mut control,
        )
        .unwrap();
        let before = session.optimizer_work();
        session.run(0).unwrap();
        assert_eq!(session.optimizer_work(), before);
        session.run(12).unwrap();
        assert!(session.point()[0] < 0.6);
        assert!(
            session.accepted().evaluation.eroded.objective.compliance
                < baseline.eroded.objective.compliance
        );
        assert!(session.best_feasible().is_some());
        session.accepted().clone()
    };
    assert_eq!(study.operator().scales(), incoming);
    let replay = study
        .evaluate_eigenstrain_robust(&accepted.rho, &cases, policy(), robustness, &mut control)
        .unwrap();
    assert_eq!(
        replay.eroded.objective.displacements,
        accepted.evaluation.eroded.objective.displacements
    );
    assert_eq!(
        replay.eroded.objective.gradient,
        accepted.evaluation.eroded.objective.gradient
    );
}
