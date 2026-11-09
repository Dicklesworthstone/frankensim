use super::super::tests::{
    Fixture, alternating, budget, limits, material, near, reference, snapshot, temperatures,
};
use super::*;
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn controls() -> ThermalJetImplicitControls {
    ThermalJetImplicitControls {
        velocity_tolerance_m_s: 1e-11,
        cell_energy_tolerance_j: 1e-11,
        solver: NewtonKrylovConfig {
            absolute_tolerance: 1e-13,
            relative_tolerance: 1e-12,
            ..NewtonKrylovConfig::default()
        },
        ..ThermalJetImplicitControls::default()
    }
}

#[test]
fn g1_joint_imex_euler_refines_against_independent_continuous_coupled_ode() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([1e-5; 3], 2.0);
    let caloric = fixture.caloric(300.0);
    let conduction = fixture.conduction(&caloric);
    for heat_transfer in [None, Some(&conduction)] {
        let expected = reference(if heat_transfer.is_some() { 1e-5 } else { 0.0 });
        let mut errors = Vec::new();
        for dt in [0.004, 0.002, 0.001] {
            let mut jet = ThermalCarreauJet::new(
                alternating(&gate),
                &caloric,
                heat_transfer,
                &temperatures(),
                100.0,
                &gate,
            )
            .unwrap();
            let report = jet
                .advance_implicit(0.02, budget(dt), controls(), &gate)
                .unwrap();
            assert_eq!(report.method, ThermalJetImplicitMethod::ImexEuler);
            assert_eq!(report.outcome, SlenderJetOutcome::Completed);
            let actual = [
                1.0 - jet.flow().radius_squared_m2()[0],
                jet.flow().velocities_m_s().next().unwrap(),
                jet.cell_internal_energy_from_reference_j()[0],
                jet.cell_internal_energy_from_reference_j()[1],
                report.viscous_heat_j,
            ];
            errors.push(std::array::from_fn::<_, 5, _>(|i| {
                (actual[i] - expected[i]).abs()
            }));
            near(report.energy_balance_residual_j, 0.0, 1e-10);
            near(report.combined_energy_residual_j, 0.0, 3e-10);
            near(report.mass_balance_residual_kg, 0.0, 2e-13);
            assert!(report.implicit_kinetic_damping_j > 0.0);
            assert!(report.nonlinear_iterations > 0 && report.linear_iterations > 0);
            assert!(report.maximum_velocity_residual_m_s <= controls().velocity_tolerance_m_s);
            assert!(report.maximum_cell_energy_residual_j <= controls().cell_energy_tolerance_j);
            let heat = jet.flow().cell_viscous_heat().unwrap();
            near(heat.total_heat_j, report.viscous_heat_j, 0.0);
            near(
                heat.cell_heat_j.iter().sum(),
                report.internal_energy_change_j,
                1e-10,
            );
            if dt == 0.001 {
                // Independently solved discrete BE equations, distinct from
                // the continuous RK4 comparison used for temporal refinement.
                let values = if heat_transfer.is_some() {
                    [
                        0.01692920033903273,
                        0.08608013280238273,
                        301.57510646682056,
                        330.56028735454476,
                        0.10227986968705381,
                        0.0009555131056788664,
                    ]
                } else {
                    [
                        0.016929072218587263,
                        0.08607780868018494,
                        301.5438535539123,
                        330.58796786319414,
                        0.10228478726276659,
                        0.0009556235370834322,
                    ]
                };
                near(actual[0], values[0], 2e-10);
                near(actual[1], values[1], 2e-10);
                let temperatures: Vec<_> = jet.temperatures_k().collect();
                near(temperatures[0], values[2], 2e-7);
                near(temperatures[1], values[3], 2e-7);
                near(report.viscous_heat_j, values[4], 2e-9);
                near(report.implicit_kinetic_damping_j, values[5], 2e-9);
            }
        }
        for pair in errors.windows(2) {
            for i in 0..5 {
                assert!(pair[1][i] < 0.56 * pair[0][i], "errors={errors:?}");
            }
        }
    }
}

#[test]
fn g1_stiff_original_rheology_and_conduction_solve_together_beyond_explicit_diffusion_bound() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([1.0; 3], 2.0);
    let caloric = fixture.caloric(300.0);
    let conduction = fixture.conduction(&caloric);
    let mut p = material();
    p.rheology = ThermalCarreauYasuda::new(
        CarreauYasuda::new(80.0, 2.0, 10.0, 0.5, 2.0, [0.0, 1000.0]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::WilliamsLandelFerry {
                c1: 1.0,
                c2_k: 100.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    );
    let u: Vec<_> = (0..16)
        .map(|i| if i % 2 == 0 { 0.125 } else { -0.125 })
        .collect();
    let flow = PeriodicSlenderJet::new_carreau(p, 4.0, &[1.0; 16], &u, limits(), &gate).unwrap();
    let mut jet = ThermalCarreauJet::new(
        flow,
        &caloric,
        Some(&conduction),
        &temperatures(),
        100.0,
        &gate,
    )
    .unwrap();
    jet.advance(0.0, budget(0.02), &gate).unwrap();
    let explicit = jet.evaluate(false, &gate).unwrap();
    assert!(0.02 > 100.0 * stable_step(explicit.rate_per_s, budget(0.02).flow.cfl));
    let report = jet
        .advance_implicit(0.02, budget(0.02), controls(), &gate)
        .unwrap();
    assert_eq!(report.accepted_substeps, 1);
    let t: Vec<_> = jet.temperatures_k().collect();
    // Independent 3-variable BE root after native conservative transport.
    near(1.0 - jet.flow().radius_squared_m2()[0], 0.02, 2e-14);
    near(
        jet.flow().velocities_m_s().next().unwrap(),
        0.0011001480640860236,
        2e-10,
    );
    near(t[0], 315.3266870935633, 2e-6);
    near(t[1], 315.4707376991972, 2e-6);
    near(report.viscous_heat_j, 0.003425798271791962, 2e-9);
    near(report.implicit_kinetic_damping_j, 0.19290853317546924, 2e-9);
    near(report.combined_energy_residual_j, 0.0, 3e-10);
    near(
        report.internal_energy_change_j,
        report.viscous_heat_j,
        1e-10,
    );
    assert!(report.implicit_kinetic_damping_j > 50.0 * report.viscous_heat_j);
    // A stable coarse BE step damps this unresolved transient numerically;
    // none of that numerical loss may masquerade as physical heat.
    let sample = jet.flow().carreau_sample_at(0).unwrap().unwrap();
    assert_eq!(sample.temperature().temperature_k(), t[0]);
    assert_eq!(sample.law(), &p.rheology);
}

#[test]
fn g1_analytic_joint_jacobian_matches_original_residual_differences_in_both_cross_blocks() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([0.12, 0.2, 0.28], 2.0);
    let caloric = fixture.caloric(300.0);
    let conduction = fixture.conduction(&caloric);
    let n = 16;
    let radii: Vec<_> = (0..n)
        .map(|i| 1.0 + 0.03 * fs_math::det::sin(core::f64::consts::TAU * i as f64 / n as f64))
        .collect();
    let velocity: Vec<_> = (0..n)
        .map(|i| 0.125 + 0.05 * fs_math::det::cos(core::f64::consts::TAU * i as f64 / n as f64))
        .collect();
    let temperature: Vec<_> = (0..n).map(|i| 300.0 + 2.0 * i as f64).collect();
    let flow = PeriodicSlenderJet::new_carreau(material(), 4.0, &radii, &velocity, limits(), &gate)
        .unwrap();
    let mut jet = ThermalCarreauJet::new(
        flow,
        &caloric,
        Some(&conduction),
        &temperature,
        100.0,
        &gate,
    )
    .unwrap();
    jet.advance(0.0, budget(0.01), &gate).unwrap();
    let mut work = Workspace::new(n, &gate).unwrap();
    let mut evaluations = 0;
    let c = ThermalJetImplicitControls {
        specific_energy_scale_j_kg: 0.1,
        ..controls()
    };
    let problem = Problem {
        context: RefCell::new(StageContext {
            jet: &mut jet,
            work: &mut work,
            dt: 0.01,
            controls: c,
            evaluations: &mut evaluations,
            gate: &gate,
            failure: None,
            trial_refusal: None,
        }),
    };
    let x = vec![0.0; 2 * n];
    for active in [0, 1, 2] {
        let direction: Vec<_> = (0..2 * n)
            .map(|i| {
                if (active == 0 && i >= n) || (active == 1 && i < n) {
                    0.0
                } else {
                    0.2 * fs_math::det::sin(0.7 * (i + 1) as f64)
                }
            })
            .collect();
        let mut actual = vec![0.0; 2 * n];
        problem.jacobian_apply(&x, &direction, &mut actual);
        assert!(problem.context.borrow().failure.is_none());
        let h = 1e-4;
        let plus: Vec<_> = x.iter().zip(&direction).map(|(a, d)| a + h * d).collect();
        let minus: Vec<_> = x.iter().zip(&direction).map(|(a, d)| a - h * d).collect();
        let mut a = vec![0.0; 2 * n];
        let mut b = vec![0.0; 2 * n];
        problem.residual(&plus, &mut a);
        problem.residual(&minus, &mut b);
        let mut scale = 0.0_f64;
        let mut error = 0.0_f64;
        for i in 0..2 * n {
            let expected = (a[i] - b[i]) / (2.0 * h);
            scale = scale.max(expected.abs());
            error = error.max((actual[i] - expected).abs());
        }
        assert!(
            error < 5e-7 * scale,
            "active={active}, error={error}, scale={scale}"
        );
        if active == 0 {
            assert!(actual[n..].iter().any(|v| v.abs() > 1e-5));
        }
        if active == 1 {
            assert!(actual[..n].iter().any(|v| v.abs() > 1e-5));
        }
    }
    // An unsupported algebraic trial may be followed by an admitted trial:
    // original line-search backtracking must not poison the accepted state.
    let mut invalid = x.clone();
    invalid[n] = 1e6;
    let mut residual = vec![0.0; 2 * n];
    problem.residual(&invalid, &mut residual);
    assert!(residual.iter().all(|v| v.is_nan()));
    assert!(problem.context.borrow().failure.is_none());
    assert!(problem.context.borrow().trial_refusal.is_some());
    problem.residual(&x, &mut residual);
    assert!(residual.iter().all(|v| v.is_finite()));
    assert!(problem.context.borrow().trial_refusal.is_none());
}

#[test]
fn g4_implicit_work_source_and_cancellation_refusals_preserve_coupled_history() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([1e-5; 3], 2.0);
    let caloric = fixture.caloric(300.0);
    let mut jet = ThermalCarreauJet::new(
        alternating(&gate),
        &caloric,
        None,
        &temperatures(),
        100.0,
        &gate,
    )
    .unwrap();
    jet.advance_implicit(0.001, budget(0.001), controls(), &gate)
        .unwrap();
    let before = snapshot(&jet);
    let coefficients = jet.flow().parameters();
    let heat = jet.flow().cell_viscous_heat().unwrap().cell_heat_j.to_vec();
    let mut limited = budget(0.001);
    limited.flow.maximum_substeps = 1;
    assert!(matches!(
        jet.advance_implicit(0.003, limited, controls(), &gate),
        Err(ThermalJetImplicitError::Thermal(ThermalJetError::Flow(
            NonlinearJetError::SubstepBudget
        )))
    ));
    assert_eq!(snapshot(&jet), before);
    assert_eq!(jet.flow().parameters(), coefficients);
    assert_eq!(jet.flow().cell_viscous_heat().unwrap().cell_heat_j, heat);
    let mut limited = controls();
    limited.maximum_residual_evaluations = 1;
    assert!(matches!(
        jet.advance_implicit(0.001, budget(0.001), limited, &gate),
        Err(ThermalJetImplicitError::EvaluationBudget)
    ));
    assert_eq!(snapshot(&jet), before);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        jet.advance_implicit(0.001, budget(0.001), controls(), &cancelled),
        Err(ThermalJetImplicitError::Thermal(ThermalJetError::Flow(
            NonlinearJetError::Cancelled
        )))
    ));
    assert_eq!(snapshot(&jet), before);
    let mut hot = ThermalCarreauJet::new(
        alternating(&gate),
        &caloric,
        None,
        &[359.999; 16],
        100.0,
        &gate,
    )
    .unwrap();
    let before_hot = snapshot(&hot);
    assert!(
        hot.advance_implicit(0.02, budget(0.001), controls(), &gate)
            .is_err()
    );
    assert_eq!(snapshot(&hot), before_hot);
}
