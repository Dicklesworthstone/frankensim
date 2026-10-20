use super::*;
use crate::sdf3::robust::envelope::tests::fixture;
use crate::{SolveBudget, SolveProgress};
use std::ops::ControlFlow;
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance*a.abs().max(b.abs()).max(1e-10), "{a} != {b}");
}
fn policy() -> LoadSpace3 { LoadSpace3::new(2, 1000.0, 0.05).unwrap() }

#[test]
fn repeated_eigenvalues_have_a_smooth_objective_and_bad_grams_refuse() {
    let p = policy();
    let e = spectrum(vec![2.0, 0.0, 0.0, 2.0], &p).unwrap();
    close(e.smooth_maximum, 2.0+p.bias_bound(), 1e-14);
    assert_eq!(e.weights, vec![0.5; 2]);
    // Crossing eigenvalues must not turn into the derivative of ONE eigenmode.
    let h = 1e-5;
    let a = spectrum(vec![2.0+h, 0.0, 0.0, 2.0], &p).unwrap();
    let b = spectrum(vec![2.0-h, 0.0, 0.0, 2.0], &p).unwrap();
    close((a.smooth_maximum-b.smooth_maximum)/(2.0*h), 0.5, 1e-8);
    assert!(spectrum(vec![2.0, 0.5, 0.1, 1.0], &p).is_err());
    assert!(spectrum(vec![1.0, 2.0, 2.0, 1.0], &p).is_err());
    assert!(spectrum(vec![f64::NAN; 4], &p).is_err());
    for n in [0, MAX_LOAD_SPACE_DIMENSION+1] { assert!(LoadSpace3::new(n, 1.0, 0.1).is_err()); }
    assert!(LoadSpace3::new(2, 0.0, 0.1).is_err());
    assert!(LoadSpace3::new(2, 1.0, f64::INFINITY).is_err());
}

#[test]
fn simultaneous_combination_exceeds_individual_cases_and_replays_actual_equilibrium() {
    let (mut study, y, z) = fixture(2.0);
    let second: Vec<_> = y.iter().zip(&z).map(|(y,z)| y+0.2*z).collect();
    let rho = vec![0.45; study.cells()];
    let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let e = study.evaluate_robust_load_space(&rho, &[&y, &second], RobustOptions3::default(), &policy(), &mut control).unwrap();
    let basis_max = e.physical.eroded.objective.case_compliances.iter().copied().fold(0.0_f64, f64::max)/policy().reference();
    assert!(e.spectrum.maximum > 1.5*basis_max);
    close(e.worst_coefficients.iter().map(|q| q*q).sum(), 1.0, 1e-12);
    close(e.worst_compliance/policy().reference(), e.spectrum.maximum, 1e-10);
    assert!(e.spectrum.smooth_maximum >= e.spectrum.maximum);
    assert!(e.spectrum.smooth_maximum <= e.spectrum.maximum+policy().bias_bound()+1e-12);
    let direct = study.evaluate_robust(&rho, &[LoadCase { force: &e.worst_force, weight: 1.0 }],
        RobustOptions3::default(), &mut control).unwrap();
    close(direct.eroded.objective.compliance, e.worst_compliance, 1e-9);
    let error: f64 = direct.eroded.objective.displacements[0].iter().zip(&e.worst_displacement)
        .map(|(a,b)| (a-b)*(a-b)).sum();
    let norm: f64 = e.worst_displacement.iter().map(|u| u*u).sum();
    assert!(error <= 1e-16*norm);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn spectral_and_dilated_gradients_match_every_raw_coordinate() {
    for beta in [1.0, 4.0] {
        let (mut study, y, z) = fixture(beta);
        let second: Vec<_> = y.iter().zip(&z).map(|(y,z)| 0.4*y+z).collect();
        let rho: Vec<_> = (0..study.cells()).map(|i| 0.38+0.025*i as f64).collect();
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        let p = policy(); let robust = RobustOptions3::default(); let basis = [&y[..], &second[..]];
        let e = study.evaluate_robust_load_space(&rho, &basis, robust, &p, &mut control).unwrap();
        for i in 0..rho.len() {
            let h = 1e-4; let mut a = rho.clone(); let mut b = rho.clone(); a[i] += h; b[i] -= h;
            let plus = study.evaluate_robust_load_space(&a, &basis, robust, &p, &mut control).unwrap();
            let minus = study.evaluate_robust_load_space(&b, &basis, robust, &p, &mut control).unwrap();
            close((plus.spectrum.smooth_maximum-minus.spectrum.smooth_maximum)/(2.0*h), e.gradient[i], 4e-4);
            close((plus.physical.dilated_volume_fraction-minus.physical.dilated_volume_fraction)/(2.0*h),
                e.physical.dilated_volume_gradient[i], 3e-5);
        }
        let check = controlled_load_space_gradient_check(&mut study, &rho, &basis, robust, &p,
            GradientCheckOptions::default(), &mut control).unwrap();
        assert!(check.numerical.passed(), "{check:?}");
    }
}

#[test]
fn orthogonal_basis_rotation_preserves_the_load_set_objective_and_gradient() {
    let (mut study, y, z) = fixture(2.0);
    let s = 0.5_f64.sqrt();
    let a: Vec<_> = y.iter().zip(&z).map(|(y,z)| s*(y+z)).collect();
    let b: Vec<_> = y.iter().zip(&z).map(|(y,z)| s*(z-y)).collect();
    let rho = vec![0.45; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let p = policy(); let robust = RobustOptions3::default();
    let first = study.evaluate_robust_load_space(&rho, &[&y,&z], robust, &p, &mut control).unwrap();
    let rotated = study.evaluate_robust_load_space(&rho, &[&a,&b], robust, &p, &mut control).unwrap();
    close(first.spectrum.maximum, rotated.spectrum.maximum, 1e-9);
    close(first.spectrum.smooth_maximum, rotated.spectrum.smooth_maximum, 1e-9);
    for (a,b) in first.gradient.iter().zip(&rotated.gradient) { close(*a,*b,1e-8); }
    assert_eq!(first.physical.dilated_rho, rotated.physical.dilated_rho);
}

#[test]
fn dependent_opposite_columns_are_a_larger_load_set_not_cancelling_scenarios() {
    let (mut study, _, z) = fixture(1.0); let rho = vec![0.45; study.cells()];
    let negative: Vec<_> = z.iter().map(|v| -v).collect();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let one = study.evaluate_robust_load_space(&rho, &[&z], RobustOptions3::default(),
        &LoadSpace3::new(1, 1000.0, 0.05).unwrap(), &mut control).unwrap();
    let two = study.evaluate_robust_load_space(&rho, &[&z,&negative], RobustOptions3::default(), &policy(), &mut control).unwrap();
    close(two.spectrum.maximum, 2.0*one.spectrum.maximum, 1e-10);
    let zeros = vec![0.0; z.len()];
    let zero = study.evaluate_robust_load_space(&rho, &[&zeros,&zeros], RobustOptions3::default(), &policy(), &mut control).unwrap();
    assert_eq!(zero.spectrum.maximum, 0.0);
    close(zero.spectrum.smooth_maximum, policy().bias_bound(), 1e-14);
    assert!(zero.gradient.iter().all(|g| *g == 0.0));
}

#[test]
fn protected_material_and_interrupted_spectrum_preserve_nominal_state() {
    use crate::sdf3::PhysicalRegion3;
    let (study, y, z) = fixture(2.0);
    let mut regions = vec![PhysicalRegion3::Design; study.cells()];
    regions[0] = PhysicalRegion3::Solid; regions[1] = PhysicalRegion3::Void;
    let mut study = study.with_physical_regions(regions).unwrap(); let rho = vec![0.45; study.cells()];
    let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let e = study.evaluate_robust_load_space(&rho, &[&y,&z], RobustOptions3::default(), &policy(), &mut control).unwrap();
    for field in [&e.physical.eroded.projected_rho, &e.physical.nominal_rho, &e.physical.dilated_rho] {
        assert_eq!(field[0],1.0); assert_eq!(field[1],0.0);
    }
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-load-space-adjoint" { ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    assert_eq!(study.evaluate_robust_load_space(&rho, &[&y,&z], RobustOptions3::default(), &policy(), &mut control).unwrap_err(), EvaluationStop::Cancelled);
    assert!(control.work().linear_iterations > 0);
    assert_eq!(study.operator().scales(),scales);
}

#[test]
fn bad_basis_and_exhausted_shared_budget_publish_no_partial_load_space() {
    let (mut study, y, z) = fixture(1.0); let rho = vec![0.4; study.cells()];
    let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(study.evaluate_robust_load_space(&rho, &[&y], RobustOptions3::default(), &policy(), &mut control).is_err());
    assert_eq!(control.work().linear_solves, 0);
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() }, &mut poll);
    assert!(matches!(study.evaluate_robust_load_space(&rho, &[&y,&z], RobustOptions3::default(), &policy(), &mut control),
        Err(EvaluationStop::TotalBudget { .. })));
    assert_eq!(control.work().linear_iterations, 1);
    assert_eq!(study.operator().scales(), scales);
}
