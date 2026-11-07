//! G1/G3/G4: adaptive SDIRK2 uses native stages, latent storage and accepted heat.

use super::*;
use fs_conduction::ThermalBc;
use fs_conduction::transient::enthalpy::EnthalpyError;
use fs_conduction::transient::enthalpy::adaptive::EnthalpySdirk2AcceptedTrial;
use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use fs_conduction::transient::enthalpy::sdirk::{EnthalpySdirk2Input, SDIRK2_GAMMA};

fn input(problem: ConductionProblem<'_>) -> EnthalpySdirk2Input<'_> {
    EnthalpySdirk2Input {
        problem,
        interfaces: None,
        radiation: EnthalpyRadiation::None,
        contacts: None,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} versus {expected:.16e}, tolerance {tolerance:e}"
    );
}

fn trace(trial: EnthalpySdirk2AcceptedTrial<'_>) -> Vec<u64> {
    let mut out = vec![trial.end_time_s.to_bits(), trial.error.normalized.to_bits()];
    for half in [trial.first_half, trial.second_half] {
        out.extend(bits(&half.first_stage.specific_enthalpy_j_kg));
        out.extend(bits(&half.second_stage.specific_enthalpy_j_kg));
        out.push(half.dirichlet_in_j.to_bits());
        out.push(half.stored_energy_change_j.to_bits());
    }
    out
}

#[test]
fn smooth_diffusion_refines_against_exact_mode_with_fewer_stages_than_backward_euler() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = sensible();
        let engine = prepare(cx, &mesh, &curve);
        let p = problem(&mesh, &boundary, &material, &source);
        let initial = [200.0, 250.0, 150.0, 200.0];
        // Independently, M^-1 K [0,1,-1,0] = 1.2 [0,1,-1,0] / s.
        let exact = 320.0 + 10.0 * (-1.2_f64 * 1.6).exp();
        let mut errors = Vec::new();
        let mut work = Vec::new();
        for tolerance in [1e-2, 1e-5] {
            let mut history = engine
                .adaptive_state(cx, &initial, 0.0, 1.0, config(tolerance))
                .unwrap();
            let result = history.advance_to_sdirk2(cx, input(p), 1.6, 10000).unwrap();
            assert_eq!(result.status, AdaptiveEnthalpyStatus::ReachedEnd);
            assert!(result.rejected > 0);
            assert_eq!(
                result.endpoint_solves,
                3 * (result.accepted + result.rejected)
            );
            errors.push((history.temperature()[1] - exact).abs());
            work.push(2 * result.endpoint_solves);
            close(result.integrated_net_input_j, 0.0, 1e-14);
            close(result.energy_residual_j, 0.0, 1e-8);
            assert!(history.last_endpoint().is_none());
            assert_eq!(
                bits(&history.last_sdirk2_step().unwrap().second_stage.temperature),
                bits(history.temperature())
            );
        }
        assert!(
            errors[1] < errors[0] / 20.0 && errors[1] < 1e-3,
            "{errors:?}"
        );
        let mut first_order = engine
            .adaptive_state(cx, &initial, 0.0, 1.0, config(1e-5))
            .unwrap();
        let reference = first_order.advance_to(cx, p, None, 1.6, 10000).unwrap();
        assert_eq!(reference.status, AdaptiveEnthalpyStatus::ReachedEnd);
        assert!(errors[1] < (first_order.temperature()[1] - exact).abs());
        assert!(
            work[1] < reference.endpoint_solves,
            "SDIRK implicit stages {}, backward Euler {}",
            work[1],
            reference.endpoint_solves
        );
    });
}

#[test]
fn latent_fraction_controls_error_even_when_plateau_temperatures_are_identical() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(1.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = latent();
        let engine = prepare(cx, &mesh, &curve);
        let mut policy = config(100.0);
        policy.enthalpy_absolute_tolerance_j_kg = 1000.0;
        policy.liquid_fraction_absolute_tolerance = 5e-7;
        let mut history = engine
            .adaptive_state(cx, &[80.0, 200.0, 200.0, 320.0], 0.0, 0.15, policy)
            .unwrap();
        let p = input(problem(&mesh, &boundary, &material, &source));
        let rejected = history.advance_to_sdirk2(cx, p, 0.15, 1).unwrap();
        assert_eq!((rejected.accepted, rejected.rejected), (0, 1));
        let error = rejected.last_error.unwrap();
        assert!(error.maximum_temperature_difference_k < policy.temperature_absolute_tolerance_k);
        assert!(error.maximum_enthalpy_difference_j_kg < policy.enthalpy_absolute_tolerance_j_kg);
        assert!(
            error.maximum_liquid_fraction_difference / 3.0
                > policy.liquid_fraction_absolute_tolerance
        );
        let report = history.advance_to_sdirk2(cx, p, 0.15, 1000).unwrap();
        assert_eq!(report.status, AdaptiveEnthalpyStatus::ReachedEnd);
        assert_eq!(&history.temperature()[1..3], &[350.0, 350.0]);
        // y=(T0-350,T3-350) obeys y'=[[-3,1],[1,-1]]y. Latent h1=h2
        // follows from the independently conserved sum, never from T=350 alone.
        let t = 0.15_f64;
        let a = 2.0_f64.sqrt() * t;
        let y0 = (-2.0 * t).exp() * (-10.0 * a.cosh() + 10.0 * 2.0_f64.sqrt() * a.sinh());
        let y3 = (-2.0 * t).exp() * (10.0 * a.cosh());
        for (&h, exact) in history.specific_enthalpy_j_kg().iter().zip([
            100.0 + 2.0 * y0,
            200.0 - y0 - y3,
            200.0 - y0 - y3,
            300.0 + 2.0 * y3,
        ]) {
            close(h, exact, 0.002);
        }
        close(report.energy_residual_j, 0.0, 1e-8);
    });
}

#[test]
fn accepted_heat_uses_both_native_stages_instead_of_endpoint_power() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "bath",
                |face| face.centroid[0] == 0.0,
                ThermalBc::dirichlet(320.0).unwrap(),
            )
            .unwrap()
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = sensible();
        let engine = prepare(cx, &mesh, &curve);
        let initial = [200.0, 400.0, 200.0, 200.0];
        let mut history = engine
            .adaptive_state(cx, &initial, 0.0, 0.5, config(1e-3))
            .unwrap();
        let mut heat = 0.0;
        let mut endpoint_power_heat = 0.0;
        let mut last_buffer = std::ptr::null();
        let result = history
            .advance_to_sdirk2_observed(
                cx,
                input(problem(&mesh, &boundary, &material, &source)),
                1.0,
                1000,
                |trial| {
                    for half in [trial.first_half, trial.second_half] {
                        // Vertex 1 alone is free: outward heat = k/6*(T1 - Tbath).
                        let first_rate = -0.5 * (half.first_stage.temperature[1] - 320.0);
                        let final_rate = -0.5 * (half.second_stage.temperature[1] - 320.0);
                        let integral = half.duration_s
                            * ((1.0 - SDIRK2_GAMMA) * first_rate + SDIRK2_GAMMA * final_rate);
                        close(half.dirichlet_in_j, integral, 1e-9);
                        heat += integral;
                        endpoint_power_heat += half.duration_s * final_rate;
                    }
                    last_buffer = trial
                        .second_half
                        .second_stage
                        .specific_enthalpy_j_kg
                        .as_ptr();
                    ControlFlow::Continue(())
                },
            )
            .unwrap();
        assert_eq!(result.status, AdaptiveEnthalpyStatus::ReachedEnd);
        close(result.integrated_net_input_j, heat, 1e-8);
        close(
            result.stored_energy_change_j,
            (history.specific_enthalpy_j_kg()[1] - 400.0) / 12.0,
            1e-9,
        );
        close(result.energy_residual_j, 0.0, 1e-8);
        assert!((heat - endpoint_power_heat).abs() > 0.01);
        assert_eq!(last_buffer, history.specific_enthalpy_j_kg().as_ptr());
        // Switching method retains the same accepted inventory and heat account.
        let switched = history
            .advance_to(
                cx,
                problem(&mesh, &boundary, &material, &source),
                None,
                1.1,
                1000,
            )
            .unwrap();
        assert_eq!(switched.status, AdaptiveEnthalpyStatus::ReachedEnd);
        assert!(history.last_sdirk2_step().is_none() && history.last_endpoint().is_some());
        close(switched.energy_residual_j, 0.0, 1e-8);
    });
}

#[test]
fn caloric_reference_offset_does_not_relax_relative_time_accuracy() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let mut outcomes = Vec::new();
        for datum in [0.0, 1_048_576.0] {
            let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
                ContentHash([0x6a; 32]),
                SolidLiquidPhase::Solid,
                [(0.0, 280.0), (1100.0, 500.0)]
                    .map(|(h, t)| EnthalpyPhaseKnot {
                        specific_enthalpy_j_kg: datum + h,
                        temperature_k: t,
                        liquid_mass_fraction: 0.0,
                        bulk_density_kg_m3: 2.0,
                    })
                    .to_vec(),
            )
            .unwrap();
            let engine = prepare(cx, &mesh, &curve);
            let mut policy = config(1e6);
            policy.enthalpy_absolute_tolerance_j_kg = 1e-6;
            policy.relative_tolerance = 1e-3;
            policy.endpoint.newton.absolute_tolerance = 1e-10;
            policy.endpoint.newton.relative_tolerance = 0.0;
            let initial = [200.0, 250.0, 150.0, 200.0].map(|h| datum + h);
            let mut history = engine
                .adaptive_state(cx, &initial, 0.0, 1.0, policy)
                .unwrap();
            let report = history
                .advance_to_sdirk2(
                    cx,
                    input(problem(&mesh, &boundary, &material, &source)),
                    1.6,
                    1000,
                )
                .unwrap();
            assert_eq!(report.status, AdaptiveEnthalpyStatus::ReachedEnd);
            assert!(report.rejected > 0);
            outcomes.push((
                history.temperature().to_vec(),
                history.accepted_steps(),
                history.rejected_steps(),
                history.next_step_s(),
            ));
        }
        assert_eq!(
            (outcomes[0].1, outcomes[0].2),
            (outcomes[1].1, outcomes[1].2)
        );
        for (&a, &b) in outcomes[0].0.iter().zip(&outcomes[1].0) {
            close(a, b, 1e-7);
        }
        close(outcomes[0].3, outcomes[1].3, 1e-7);
    });
}

#[test]
fn rejected_comparisons_and_trial_budget_resume_preserve_native_stage_history() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(10.0);
        let curve = sensible();
        let engine = prepare(cx, &mesh, &curve);
        let initial = [200.0, 250.0, 150.0, 200.0];
        let p = input(problem(&mesh, &boundary, &material, &source));
        let mut straight = engine
            .adaptive_state(cx, &initial, 0.0, 1.0, config(1e-5))
            .unwrap();
        let mut resumed = straight.clone();
        let mut expected = Vec::new();
        let full = straight
            .advance_to_sdirk2_observed(cx, p, 1.0, 1000, |trial| {
                expected.push(trace(trial));
                ControlFlow::Continue(())
            })
            .unwrap();
        let rejected = resumed
            .advance_to_sdirk2_observed(cx, p, 1.0, 1, |_| {
                panic!("a rejected SDIRK comparison must not publish native stages")
            })
            .unwrap();
        assert_eq!((rejected.accepted, rejected.rejected), (0, 1));
        assert_eq!(bits(resumed.specific_enthalpy_j_kg()), bits(&initial));
        assert_eq!(
            (resumed.time_s(), rejected.integrated_net_input_j),
            (0.0, 0.0)
        );
        let mut actual = Vec::new();
        for _ in 0..1000 {
            let report = resumed
                .advance_to_sdirk2_observed(cx, p, 1.0, 1, |trial| {
                    actual.push(trace(trial));
                    ControlFlow::Continue(())
                })
                .unwrap();
            if report.status == AdaptiveEnthalpyStatus::ReachedEnd {
                break;
            }
            assert_eq!(report.status, AdaptiveEnthalpyStatus::TrialLimit);
        }
        assert_eq!(resumed.time_s(), 1.0);
        assert_eq!(expected, actual);
        assert_eq!(
            bits(straight.specific_enthalpy_j_kg()),
            bits(resumed.specific_enthalpy_j_kg())
        );
        assert_eq!(
            straight.next_step_s().to_bits(),
            resumed.next_step_s().to_bits()
        );
        assert_eq!(straight.endpoint_solves(), resumed.endpoint_solves());
        let done = resumed.advance_to_sdirk2(cx, p, 1.0, 0).unwrap();
        assert_eq!(
            done.integrated_net_input_j.to_bits(),
            full.integrated_net_input_j.to_bits()
        );
        close(done.integrated_net_input_j, 10.0 / 6.0, 1e-12);
        close(done.energy_residual_j, 0.0, 1e-8);
    });
}

#[test]
fn minimum_step_and_observer_cancellation_never_publish_partial_stages() {
    support::with_gate(|gate, cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = sensible();
        let engine = prepare(cx, &mesh, &curve);
        let initial = [200.0, 250.0, 150.0, 200.0];
        let p = input(problem(&mesh, &boundary, &material, &source));
        let mut policy = config(1e-10);
        policy.minimum_step_s = 0.5;
        policy.maximum_step_s = 0.5;
        let mut refused = engine
            .adaptive_state(cx, &initial, 0.0, 0.5, policy)
            .unwrap();
        let report = refused
            .advance_to_sdirk2(cx, p, 0.5_f64.next_up(), 100)
            .unwrap();
        assert_eq!(report.status, AdaptiveEnthalpyStatus::MinimumStep);
        assert_eq!((report.accepted, report.rejected), (0, 1));
        assert_eq!(refused.time_s(), 0.0);
        assert_eq!(bits(refused.specific_enthalpy_j_kg()), bits(&initial));
        let mut history = engine
            .adaptive_state(cx, &initial, 0.0, 0.1, config(1e-3))
            .unwrap();
        let mut observed = None;
        let result = history.advance_to_sdirk2_observed(cx, p, 1.0, 100, |trial| {
            observed = Some((
                trial.end_time_s.to_bits(),
                bits(&trial.second_half.second_stage.specific_enthalpy_j_kg),
            ));
            gate.request();
            ControlFlow::Continue(())
        });
        assert!(matches!(
            result,
            Err(EnthalpyError::Conduction(
                fs_conduction::ConductionError::Cancelled { .. }
            ))
        ));
        let (time, field) = observed.unwrap();
        assert_eq!(history.time_s().to_bits(), time);
        assert_eq!(bits(history.specific_enthalpy_j_kg()), field);
        assert_eq!(history.accepted_steps(), 1);
        let before = history.endpoint_solves();
        assert!(history.advance_to_sdirk2(cx, p, 1.0, 1).is_err());
        assert_eq!(history.endpoint_solves(), before);
        assert_eq!(bits(history.specific_enthalpy_j_kg()), field);
        let cloned = history.clone();
        assert_eq!(
            bits(&cloned.last_sdirk2_step().unwrap().first_stage.temperature),
            bits(&history.last_sdirk2_step().unwrap().first_stage.temperature)
        );
    });
}
