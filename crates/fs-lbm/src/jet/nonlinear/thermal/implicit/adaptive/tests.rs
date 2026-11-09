//! G1/G4: heat-aware stiff resolution and whole-call transaction boundaries.

use super::super::super::tests::{
    Fixture, alternating, budget, limits, material, near, snapshot, temperatures,
};
use super::*;
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn nonlinear_controls() -> ThermalJetImplicitControls {
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

fn broad_local_controls(damping_j: f64) -> ThermalJetAdaptiveControls {
    let absolute = |value| ThermalJetAdaptiveTolerance {
        absolute: value,
        relative: 0.0,
    };
    ThermalJetAdaptiveControls {
        minimum_interval_s: 1e-12,
        maximum_refinements: 40,
        area_m2: absolute(1.0),
        radius_m: absolute(1.0),
        velocity_m_s: absolute(1.0),
        temperature_k: absolute(100.0),
        specific_energy_j_kg: absolute(1.0),
        cell_heat_j: absolute(1.0),
        total_heat_j: absolute(1.0),
        absolute_numerical_damping_j: damping_j,
        relative_numerical_damping: 0.0,
    }
}

#[test]
fn g1_adaptive_damping_budget_resolves_physical_heat_missed_by_stable_coarse_be() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([1.0; 3], 2.0);
    let caloric = fixture.caloric(300.0);
    let conduction = fixture.conduction(&caloric);
    let mut parameters = material();
    parameters.rheology = ThermalCarreauYasuda::new(
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
    let velocities: Vec<_> = (0..16)
        .map(|i| if i % 2 == 0 { 0.125 } else { -0.125 })
        .collect();
    let mut step_budget = budget(0.02);
    step_budget.energy_tolerance_j = 1e-8;
    let make_flow = || {
        PeriodicSlenderJet::new_carreau(parameters, 4.0, &[1.0; 16], &velocities, limits(), &gate)
            .unwrap()
    };
    let mut coarse = ThermalCarreauJet::new(
        make_flow(),
        &caloric,
        Some(&conduction),
        &temperatures(),
        100.0,
        &gate,
    )
    .unwrap();
    let coarse_report = coarse
        .advance_implicit(0.02, step_budget, nonlinear_controls(), &gate)
        .unwrap();
    assert_eq!(coarse_report.accepted_substeps, 1);
    assert!(coarse_report.viscous_heat_j < 0.004);
    assert!(coarse_report.implicit_kinetic_damping_j > 0.19);

    // An independent conservative reduced ODE, integrated by RK4 at 100k and
    // 200k intervals, agrees on these values. The equilibrium temperature also
    // follows from Cv(T)=1e-5*T and converting the initial alternating kinetic
    // energy to heat. No native constitutive or time-step code entered that oracle.
    let expected_heat = 0.196349540849355;
    let expected_temperature = 317.824637182205;
    let expected_area_amplitude = 0.00041857904457331;
    let mut errors = Vec::new();
    for damping_budget in [0.002, 0.0005] {
        let mut jet = ThermalCarreauJet::new(
            make_flow(),
            &caloric,
            Some(&conduction),
            &temperatures(),
            100.0,
            &gate,
        )
        .unwrap();
        let report = jet
            .advance_implicit_adaptive(
                0.02,
                step_budget,
                nonlinear_controls(),
                broad_local_controls(damping_budget),
                &gate,
            )
            .unwrap();
        let evolution = report.evolution;
        assert_eq!(evolution.method, ThermalJetImplicitMethod::ImexEuler);
        assert_eq!(evolution.outcome, SlenderJetOutcome::Completed);
        assert_eq!(evolution.end_time_s, 0.02);
        assert!(report.rejected_intervals > 0);
        assert!(report.accepted_intervals > 1 && report.accepted_intervals < 1000);
        assert_eq!(evolution.accepted_substeps, 2 * report.accepted_intervals);
        assert!(evolution.attempted_substeps >= 3 * report.accepted_intervals);
        assert!(evolution.attempted_substeps > evolution.accepted_substeps);
        assert_eq!(
            report.attempted_intervals,
            report.accepted_intervals + report.rejected_intervals
        );
        assert!(report.maximum_accepted_local.scaled_maximum <= 1.0);
        assert!(report.maximum_accepted_damping_ratio <= 1.0);
        assert!(report.minimum_accepted_interval_s > 0.0);
        assert!(report.maximum_accepted_interval_s <= step_budget.flow.maximum_step_s);
        near(report.numerical_damping_allowance_j, damping_budget, 0.0);
        assert!(evolution.implicit_kinetic_damping_j >= 0.0);
        assert!(evolution.implicit_kinetic_damping_j <= damping_budget);
        assert!(evolution.viscous_heat_j > 0.1943 && evolution.viscous_heat_j < expected_heat);
        near(evolution.energy_balance_residual_j, 0.0, 1e-8);
        near(evolution.combined_energy_residual_j, 0.0, 2e-8);
        near(evolution.mass_balance_residual_kg, 0.0, 2e-11);

        // Counting coarse/rejected heat or heating BE damping would fail this
        // accepted physical account, even with small endpoint discrepancies.
        let heat = jet.flow().cell_viscous_heat().unwrap();
        assert!(heat.cell_heat_j.iter().all(|q| q.is_finite() && *q >= 0.0));
        near(heat.total_heat_j, evolution.viscous_heat_j, 0.0);
        near(
            heat.cell_heat_j.iter().sum(),
            evolution.viscous_heat_j,
            1e-10,
        );
        near(
            evolution.internal_energy_change_j,
            evolution.viscous_heat_j,
            1e-8,
        );
        near(
            evolution.viscous_heat_j + evolution.implicit_kinetic_damping_j,
            expected_heat,
            2e-8,
        );
        let final_temperatures: Vec<_> = jet.temperatures_k().collect();
        for (cell, &temperature) in final_temperatures.iter().enumerate() {
            near(temperature, expected_temperature, 0.04);
            let sample = jet.flow().carreau_sample_at(cell).unwrap().unwrap();
            assert_eq!(sample.temperature().temperature_k(), temperature);
            assert_eq!(sample.law(), &parameters.rheology);
        }
        let area_amplitude = 1.0 - jet.flow().radius_squared_m2()[0];
        near(area_amplitude, expected_area_amplitude, 1e-5);
        errors.push([
            (evolution.viscous_heat_j - expected_heat).abs(),
            (final_temperatures[0] - expected_temperature).abs(),
            (area_amplitude - expected_area_amplitude).abs(),
        ]);
    }
    // Tightening only the physical damping allowance improves heat, temperature
    // and deformation even though endpoint tolerances were deliberately broad.
    for i in 0..3 {
        assert!(errors[1][i] < 0.4 * errors[0][i], "errors={errors:?}");
    }
}

#[test]
fn g4_adaptive_attempt_source_residual_and_cancellation_failures_are_atomic() {
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
    jet.advance_implicit(0.001, budget(0.001), nonlinear_controls(), &gate)
        .unwrap();
    let before = snapshot(&jet);
    let coefficients = jet.flow().parameters();
    let published_heat = jet.flow().cell_viscous_heat().unwrap().cell_heat_j.to_vec();
    let unchanged = |candidate: &ThermalCarreauJet<'_>| {
        assert_eq!(snapshot(candidate), before);
        assert_eq!(candidate.flow().parameters(), coefficients);
        assert_eq!(
            candidate.flow().cell_viscous_heat().unwrap().cell_heat_j,
            published_heat
        );
    };
    let adaptive = broad_local_controls(1.0);
    // Three native solves complete the first comparison; this budget then
    // refuses during a later comparison. Earlier internal acceptances must
    // remain unpublished together with their local heat.
    let mut limited_steps = budget(0.001);
    limited_steps.flow.maximum_substeps = 4;
    assert!(matches!(
        jet.advance_implicit_adaptive(0.003, limited_steps, nonlinear_controls(), adaptive, &gate,),
        Err(ThermalJetAdaptiveError::Implicit(
            ThermalJetImplicitError::Thermal(ThermalJetError::Flow(
                NonlinearJetError::SubstepBudget
            ))
        ))
    ));
    unchanged(&jet);
    let mut limited_newton = nonlinear_controls();
    limited_newton.maximum_residual_evaluations = 1;
    assert!(matches!(
        jet.advance_implicit_adaptive(0.001, budget(0.001), limited_newton, adaptive, &gate,),
        Err(ThermalJetAdaptiveError::Implicit(
            ThermalJetImplicitError::EvaluationBudget
        ))
    ));
    unchanged(&jet);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        jet.advance_implicit_adaptive(
            0.001,
            budget(0.001),
            nonlinear_controls(),
            adaptive,
            &cancelled,
        ),
        Err(ThermalJetAdaptiveError::Implicit(
            ThermalJetImplicitError::Thermal(ThermalJetError::Flow(NonlinearJetError::Cancelled))
        ))
    ));
    unchanged(&jet);
    let no_refinement = ThermalJetAdaptiveControls {
        absolute_numerical_damping_j: 0.0,
        maximum_refinements: 0,
        ..adaptive
    };
    assert!(matches!(
        jet.advance_implicit_adaptive(
            0.001,
            budget(0.001),
            nonlinear_controls(),
            no_refinement,
            &gate,
        ),
        Err(ThermalJetAdaptiveError::Resolution { .. })
    ));
    unchanged(&jet);

    let retry = jet
        .advance_implicit_adaptive(0.002, budget(0.001), nonlinear_controls(), adaptive, &gate)
        .unwrap();
    assert_eq!(retry.evolution.outcome, SlenderJetOutcome::Completed);
    near(jet.flow().time_s(), 0.003, 1e-16);
    near(
        retry.evolution.internal_energy_change_j,
        retry.evolution.viscous_heat_j,
        1e-10,
    );

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
        hot.advance_implicit_adaptive(0.02, budget(0.001), nonlinear_controls(), adaptive, &gate,)
            .is_err()
    );
    assert_eq!(snapshot(&hot), before_hot);
    assert_eq!(
        hot.flow().cell_viscous_heat().unwrap().cell_heat_j,
        &[0.0; 16]
    );
}

#[test]
fn g4_sub_ulp_positive_damping_cannot_spend_an_exhausted_allowance() {
    let controls = broad_local_controls(0.01);
    let spent = Ledger::default().append(0.0, 0.01);
    let remaining = spent.remaining(controls).unwrap();
    assert_eq!(remaining, 0.0);
    assert_eq!(0.01 + 1e-20, 0.01, "fixture must hide the increment in f64");
    assert!(ratio(1e-20, remaining).is_infinite());
    assert!(!spent.append(0.0, 1e-20).allowed(controls).unwrap());
    assert!(spent.append(0.0, 0.0).allowed(controls).unwrap());
    assert_eq!(
        ratio(0.0, remaining),
        0.0,
        "a quiet tail must remain admissible"
    );
}

#[test]
fn g3_local_indicators_ignore_baselines_and_preserve_opposing_cell_heat_errors() {
    let tolerance = ThermalJetAdaptiveTolerance {
        absolute: 0.0,
        relative: 0.1,
    };
    for baseline in [300.0, 1_099_511_627_776.0] {
        let (mut error, mut scaled) = (0.0, 0.0);
        component(
            Dd::from_f64(baseline),
            Dd::from_f64(baseline + 8.0),
            Dd::from_f64(baseline + 10.0),
            tolerance,
            &mut error,
            &mut scaled,
        )
        .unwrap();
        near(error, 2.0, 0.0);
        // Ten units of interval change give one unit of relative allowance.
        // A Kelvin datum or uniform translational baseline contributes none.
        near(scaled, 2.0, 1e-14);
    }

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
    jet.advance_implicit(0.0, budget(0.001), nonlinear_controls(), &gate)
        .unwrap();
    jet.evaluate_transport(false, &gate).unwrap();
    for i in 0..jet.accepted.len() {
        jet.flow.stage.a[i] = jet.flow.candidate.a[i];
        jet.flow.stage.q[i] = jet.flow.candidate.q[i];
        jet.stage[i] = jet.candidate[i];
    }
    let mut work = Workspace::new(jet.accepted.len(), &gate).unwrap();
    work.base
        .save(
            &jet.flow.candidate,
            &jet.candidate,
            &jet.workspace.profile,
            &gate,
        )
        .unwrap();
    work.coarse
        .save(&jet.flow.stage, &jet.stage, &jet.workspace.profile, &gate)
        .unwrap();
    // Exercise the actual cellwise indicator with nonnegative original heat.
    // The artificial comparison changes only its spatial distribution; the
    // aggregate signed discrepancy is zero, so it cannot be the error metric.
    jet.flow.spatial.start_viscous_heat(1.0, &gate).unwrap();
    let fine = jet.flow.spatial.trial_viscous_heat_j().unwrap().to_vec();
    work.coarse_heat.copy_from_slice(&fine);
    let difference = 0.25 * fine[0].min(fine[1]);
    assert!(difference > 0.0);
    work.coarse_heat[0] += difference;
    work.coarse_heat[1] -= difference;
    let coarse_total: f64 = work.coarse_heat.iter().sum();
    let fine_total: f64 = fine.iter().sum();
    near(coarse_total, fine_total, 1e-12);
    let controls = ThermalJetAdaptiveControls {
        total_heat_j: ThermalJetAdaptiveTolerance {
            absolute: difference,
            relative: 0.0,
        },
        ..broad_local_controls(1.0)
    };
    let indicators = work
        .indicators(&jet, coarse_total, fine_total, controls, &gate)
        .unwrap();
    near(indicators.cell_heat_j, difference, 1e-12 * difference);
    near(
        indicators.total_heat_j,
        2.0 * difference,
        1e-12 * difference,
    );
    near(indicators.scaled_maximum, 2.0, 1e-12);
}
