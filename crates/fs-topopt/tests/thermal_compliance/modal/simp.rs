use super::*;
use fs_topopt::RobustPipeline;
use fs_topopt::simp_second_order::SimpThermalComplianceGoal;

#[test]
fn ordinary_and_robust_frequency_design_include_stress_equilibrium_and_reference_inertia() {
    for robust in [false, true] {
        let (mesh, positions, elasticity, thermal, a, b) = cube();
        let params = SimpParams {
            e_min: 0.03,
            penal: 2.7,
            beta: 2.3,
            eta: 0.5,
        };
        let ordinary = DesignPipeline {
            filter: DensityFilter::new(&mesh, &positions, 0.15),
            params,
        };
        let three_field = RobustPipeline {
            filter: DensityFilter::new(&mesh, &positions, 0.15),
            params,
            eta_offset: 0.08,
        };
        let reference = DesignPipeline {
            filter: DensityFilter::new(&mesh, &positions, 0.15),
            params: SimpParams {
                eta: if robust { 0.58 } else { 0.5 },
                ..params
            },
        };
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
        let rho: Vec<_> = (0..elasticity.cells())
            .map(|i| 0.42 + 0.03 * i as f64)
            .collect();
        let mut callback = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let goal = if robust {
            SimpThermalComplianceGoal::new_robust(
                &three_field,
                &elasticity,
                &rho,
                &loads,
                options(),
                &mut control,
            )
        } else {
            SimpThermalComplianceGoal::new(
                &ordinary,
                &elasticity,
                &rho,
                &loads,
                options(),
                &mut control,
            )
        }
        .unwrap();
        let count = elasticity.free().iter().filter(|f| **f).count();
        let before = control.work();
        let result = goal
            .prestressed_frequency_goal(1, modal_options(count), &mut control)
            .unwrap();
        assert_eq!(control.work().linear_solves, before.linear_solves + 2);
        assert_eq!(result.work(), control.work());
        assert!(std::ptr::eq(result.source(), &goal));
        assert!(std::ptr::eq(result.physical().source(), goal.physical()));
        assert_eq!(
            result.physical().mass_fractions(),
            goal.projected_densities()
        );
        if robust {
            assert!(
                goal.projected_densities()
                    .iter()
                    .zip(goal.volume_projected_densities())
                    .all(|(eroded, dilated)| eroded < dilated)
            );
        }

        // The independent first-order pipeline supplies the entire physical
        // point, without using the new design-map pullback.
        let mut original_value = |point: &[f64]| {
            let (_, projected, moduli) = reference.try_forward(point, &mut control).unwrap();
            let physical = solve(&elasticity, &moduli, &loads);
            frequency(&physical, 1, &projected, count).value()
        };
        close(result.value(), original_value(&rho), 2e-9);
        let h = 2e-5;
        for cell in 0..rho.len() {
            let mut plus = rho.clone();
            plus[cell] += h;
            let mut minus = rho.clone();
            minus[cell] -= h;
            let fd = (original_value(&plus) - original_value(&minus)) / (2.0 * h);
            close(result.gradient()[cell], fd, 1e-5);
        }
    }
}

#[test]
fn low_density_thermal_frequency_uses_actual_occupancy_without_a_hidden_void_mass_law() {
    let (elasticity, material) = one_tet(false);
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let positions = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    let pipeline = DesignPipeline {
        filter: DensityFilter::new(&mesh, &positions, 0.1),
        params: SimpParams {
            e_min: 0.05,
            penal: 3.0,
            beta: 0.0,
            eta: 0.5,
        },
    };
    let thermal = state(&material, [0., 0., 0.01, 0., 0., 0.], 0.1);
    let force = vec![0.0; elasticity.n()];
    let loads = [ThermalLoadCase {
        force: &force,
        weight: 1.0,
        strains: TetThermalStrainField::Uniform(&thermal),
    }];
    let r = 0.06_f64;
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let goal = SimpThermalComplianceGoal::new(
        &pipeline,
        &elasticity,
        &[r],
        &loads,
        options(),
        &mut control,
    )
    .unwrap();
    let result = goal
        .prestressed_frequency_goal(0, modal_options(3), &mut control)
        .unwrap();
    let m = 0.05 + 0.95 * r * r * r;
    let dm = 2.85 * r * r;
    let low = 100.0 * m / r;
    let high = 200.0 * m / r;
    let (value, w, z) = two_equal_low_modes(low, high, modal_options(3).aggregation_beta);
    close(result.value(), value, 2e-9);
    close(
        result.gradient()[0],
        (200.0 * w + 200.0 * z) * (dm * r - m) / (r * r),
        2e-8,
    );
    close(goal.mass_value_gradient(&mut control).unwrap().0, r, 1e-10);
    for (&actual, expected) in result
        .physical()
        .modal()
        .values
        .iter()
        .zip([low, low, high])
    {
        close(actual, expected, 2e-9);
    }
    // Cancelling the single raw pullback withholds the whole raw objective.
    let mut stop = |progress: fs_topopt::SolveProgress| {
        if progress.stage == "simp-thermal-frequency-pullback" {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        }
    };
    let mut cancelled = SolveControl::new(SolveBudget::default(), &mut stop);
    assert!(matches!(
        goal.prestressed_frequency_goal(0, modal_options(3), &mut cancelled),
        Err(GeneralizedError::Operator(ComplianceError::Stopped(
            EvaluationStop::Cancelled
        )))
    ));
    close(
        goal.prestressed_frequency_goal(0, modal_options(3), &mut control)
            .unwrap()
            .value(),
        result.value(),
        1e-12,
    );
}
