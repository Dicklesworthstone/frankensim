//! G0/G3/G4: physical stress and exact curvature at original thermal equilibria.
use super::*;
use fs_topopt::elasticity::second_order::thermal::stress::ThermalStressOptions;

fn criterion(power: u32) -> ThermalStressOptions {
    ThermalStressOptions { allowable_pa: 1.0, power }
}

fn independent_moment(goal: &ThermalComplianceGoal<'_>, options: ThermalStressOptions) -> f64 {
    (0..goal.case_compliances().len()).map(|index| {
        goal.response(index).unwrap().elements.iter().map(|element| {
            let s = &element.effective_stress_mandel_pa;
            // Physical shear is Mandel shear / sqrt(2); compute the invariant
            // independently of the production stress moment and its derivatives.
            let vm = (0.5 * (
                (s[0] - s[1]).powi(2) + (s[1] - s[2]).powi(2) + (s[2] - s[0]).powi(2)
            ) + 1.5 * (s[3].powi(2) + s[4].powi(2) + s[5].powi(2))).sqrt();
            (vm / options.allowable_pa).powi(options.power as i32)
        }).sum::<f64>()
    }).sum()
}

#[test]
fn restrained_thermal_stress_remains_live_when_compliance_and_case_weight_are_zero() {
    let (elasticity, material) = one_tet(true);
    let thermal = state(&material, [0.0, 0.0, 0.01, 0.0, 0.0, 0.0], 0.1);
    let force = vec![0.0; elasticity.n()];
    let loads = [ThermalLoadCase {
        force: &force,
        weight: 0.0,
        strains: TetThermalStrainField::Uniform(&thermal),
    }];
    let m: f64 = 0.6;
    let goal = solve(&elasticity, &[m], &loads);
    assert_eq!(goal.value(), 0.0);
    assert_eq!(goal.gradient(), [0.0]);
    for p in [2, 3, 4, 8, 64] {
        let mut callback = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let stress = goal.stress_goal(criterion(p), &mut control).unwrap();
        let ratio = 1.2 * m;
        close(stress.maximum_ratio(), ratio, 1e-12);
        close(stress.value(), ratio.powi(p as i32), 1e-11);
        close(stress.gradient()[0], f64::from(p) * ratio.powi(p as i32) / m, 1e-10);
        assert_eq!(control.work().linear_solves, 1);
        let hv = stress.hessian_vector_product(&[0.2], &mut control).unwrap();
        close(
            hv[0],
            0.2 * f64::from(p * (p - 1)) * ratio.powi(p as i32) / (m * m),
            1e-9,
        );
        assert_eq!(control.work().linear_solves, 3);
        assert!(stress.adjoint_relative_residuals().iter().all(|r| *r < options().tolerance));
        assert!(std::ptr::eq(stress.source(), &goal));
    }
}

#[test]
fn physical_traction_stress_cancels_design_stiffness_and_thermal_load_derivatives() {
    let (elasticity, material) = one_tet(false);
    let thermal = state(&material, [0.0, 0.0, 0.01, 0.0, 0.0, 0.0], 0.1);
    let mut force = vec![0.0; elasticity.n()];
    force[11] = 0.02;
    let loads = [ThermalLoadCase {
        force: &force,
        weight: 1.0,
        strains: TetThermalStrainField::Uniform(&thermal),
    }];
    let goal = solve(&elasticity, &[0.4], &loads);
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let stress = goal.stress_goal(criterion(4), &mut control).unwrap();
    // For this free axial coordinate K_e=E/6 and strain=u_z, hence sigma=6F.
    // The nodal loading is not a uniform face traction.
    close(stress.maximum_ratio(), 6.0 * force[11], 1e-10);
    close(stress.gradient()[0], 0.0, 1e-7);
    close(stress.hessian_vector_product(&[0.3], &mut control).unwrap()[0], 0.0, 1e-7);
    let before = control.work();
    assert_eq!(stress.hessian_vector_product(&[0.0], &mut control).unwrap(), [0.0]);
    assert_eq!(before, control.work());
}

#[test]
fn anisotropic_multicase_stress_adjoint_and_hessian_match_original_response_perturbations() {
    let (_, _, elasticity, thermal, a, b) = cube();
    let loads = [
        ThermalLoadCase {
            force: &a,
            weight: 1.0,
            strains: TetThermalStrainField::PerElement(&thermal),
        },
        ThermalLoadCase {
            force: &b,
            weight: 0.0,
            strains: TetThermalStrainField::PerElement(&thermal),
        },
    ];
    let m: Vec<_> = (0..elasticity.cells()).map(|i| 0.5 + 0.04 * i as f64).collect();
    let direction: Vec<_> = (0..m.len()).map(|i| 0.1 * ((i % 3) as f64 - 1.0)).collect();
    let other: Vec<_> = (0..m.len()).map(|i| 0.02 * (i + 1) as f64).collect();
    let options = ThermalStressOptions { allowable_pa: 0.2, power: 4 };
    let goal = solve(&elasticity, &m, &loads);
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let stress = goal.stress_goal(options, &mut control).unwrap();
    close(stress.value(), independent_moment(&goal, options), 2e-12);
    assert_eq!(control.work().linear_solves, 2);
    let hv = stress.hessian_vector_product(&direction, &mut control).unwrap();
    assert_eq!(control.work().linear_solves, 6);
    let hw = stress.hessian_vector_product(&other, &mut control).unwrap();
    let step = 1e-5;
    for i in 0..m.len() {
        let mut plus = m.clone();
        let mut minus = m.clone();
        plus[i] += step;
        minus[i] -= step;
        let fd = (
            independent_moment(&solve(&elasticity, &plus, &loads), options)
                - independent_moment(&solve(&elasticity, &minus, &loads), options)
        ) / (2.0 * step);
        close(stress.gradient()[i], fd, 2e-5);
    }
    let plus: Vec<_> = m.iter().zip(&direction).map(|(m, v)| m + step * v).collect();
    let minus: Vec<_> = m.iter().zip(&direction).map(|(m, v)| m - step * v).collect();
    let gp = solve(&elasticity, &plus, &loads);
    let gm = solve(&elasticity, &minus, &loads);
    let sp = gp.stress_goal(options, &mut control).unwrap();
    let sm = gm.stress_goal(options, &mut control).unwrap();
    for i in 0..m.len() {
        close(hv[i], (sp.gradient()[i] - sm.gradient()[i]) / (2.0 * step), 3e-5);
    }
    close(
        direction.iter().zip(&hw).map(|(v, h)| v * h).sum(),
        other.iter().zip(&hv).map(|(v, h)| v * h).sum(),
        2e-8,
    );
    // Weight zero is still a physical experiment and contributes to the limit.
    let only = solve(&elasticity, &m, &loads[..1]);
    assert!(stress.value() > independent_moment(&only, options));
    assert_eq!(elasticity.moduli, vec![1.0; elasticity.cells()]);
}

#[test]
fn hydrostatic_stress_has_zero_moment_without_a_floor_and_mandel_shear_is_physical() {
    let (elasticity, material) = one_tet(true);
    let force = vec![0.0; elasticity.n()];
    for (free, expected) in [
        ([0.01, 0.01, 0.01, 0.0, 0.0, 0.0], 0.0),
        ([0.0, 0.0, 0.0, 0.01, 0.0, 0.0], 1.5_f64.sqrt() * 1.2 * 0.4),
    ] {
        let thermal = state(&material, free, 0.1);
        let goal = solve(&elasticity, &[0.4], &[ThermalLoadCase {
            force: &force,
            weight: 1.0,
            strains: TetThermalStrainField::Uniform(&thermal),
        }]);
        for power in [2, 4] {
            let mut callback = |_| ControlFlow::Continue(());
            let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
            let stress = goal.stress_goal(criterion(power), &mut control).unwrap();
            close(stress.maximum_ratio(), expected, 1e-12);
            let hv = stress.hessian_vector_product(&[0.1], &mut control).unwrap();
            if expected == 0.0 {
                assert_eq!(stress.value(), 0.0);
                assert_eq!(stress.gradient(), [0.0]);
                assert_eq!(hv, [0.0]);
            } else {
                close(
                    stress.gradient()[0],
                    f64::from(power) * expected.powi(power as i32) / 0.4,
                    1e-10,
                );
                close(
                    hv[0],
                    0.1 * f64::from(power * (power - 1)) * expected.powi(power as i32) / 0.16,
                    1e-10,
                );
            }
        }
    }
}

#[test]
fn stress_work_refusals_leave_the_original_thermal_point_usable() {
    let (_, _, elasticity, thermal, a, _) = cube();
    let goal = solve(&elasticity, &vec![0.6; elasticity.cells()], &[ThermalLoadCase {
        force: &a,
        weight: 1.0,
        strains: TetThermalStrainField::PerElement(&thermal),
    }]);
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    for options in [
        ThermalStressOptions { allowable_pa: 0.0, power: 4 },
        ThermalStressOptions { allowable_pa: 1.0, power: 1 },
        ThermalStressOptions { allowable_pa: f64::NAN, power: 4 },
        ThermalStressOptions { allowable_pa: 1.0, power: 65 },
    ] {
        assert!(matches!(
            goal.stress_goal(options, &mut control),
            Err(ComplianceError::InvalidInput(_))
        ));
    }
    assert_eq!(control.work().linear_solves, 0);
    let mut cancel = |_| ControlFlow::Break(());
    assert!(matches!(
        goal.stress_goal(
            criterion(4),
            &mut SolveControl::new(SolveBudget::default(), &mut cancel),
        ),
        Err(ComplianceError::Stopped(EvaluationStop::Cancelled))
    ));
    let mut none = SolveControl::new(
        SolveBudget { per_solve_iterations: 0, total_iterations: 0 },
        &mut callback,
    );
    assert!(matches!(
        goal.stress_goal(criterion(4), &mut none),
        Err(ComplianceError::Stopped(_))
    ));
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let stress = goal.stress_goal(criterion(4), &mut control).unwrap();
    let before = control.work();
    assert!(matches!(
        stress.hessian_vector_product(&[], &mut control),
        Err(ComplianceError::InvalidInput(_))
    ));
    assert!(matches!(
        stress.hessian_vector_product(&vec![f64::NAN; elasticity.cells()], &mut control),
        Err(ComplianceError::NonFinite(_))
    ));
    assert_eq!(before, control.work());
    assert!(stress.value().is_finite() && stress.value() > 0.0);
}
