//! G1/G3/G4: native SDIRK2 entropy across accepted adaptive history.

use super::*;
use fs_conduction::ThermalBc;
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

fn close(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual - expected).abs() <= allowance,
        "{actual:.16e} versus {expected:.16e}, allowance {allowance:e}"
    );
}

#[test]
fn audited_sdirk_keeps_original_fields_rejected_work_and_exact_resume() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder().finish().unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = sensible();
        let storage = prepare(cx, &mesh, &curve);
        let p = input(problem(&mesh, &boundary, &material, &source));
        let initial = [200.0, 250.0, 150.0, 200.0];
        let policy = config(1e-5);
        let mut plain = storage.adaptive_state(cx, &initial, 0.0, 1.0, policy).unwrap();
        let mut full = storage.adaptive_state_with_entropy(
            cx, &initial, 0.0, 1.0, policy, 1e-9,
        ).unwrap();
        let mut resumed = full.clone();

        let rejected = resumed.advance_to_sdirk2(cx, p, 1.6, 1).unwrap();
        assert_eq!((rejected.thermal.accepted, rejected.thermal.rejected), (0, 1));
        assert_eq!(rejected.entropy, EnthalpyEntropyBalance::default());
        assert_eq!(bits(resumed.history().specific_enthalpy_j_kg()), bits(&initial));

        let native = plain.advance_to_sdirk2(cx, p, 1.6, 10000).unwrap();
        let audited = full.advance_to_sdirk2(cx, p, 1.6, 10000).unwrap();
        assert_eq!(native, audited.thermal);
        let final_report = loop {
            let next = resumed.advance_to_sdirk2(cx, p, 1.6, 1).unwrap();
            if next.thermal.status == AdaptiveEnthalpyStatus::ReachedEnd {
                break next;
            }
            assert_eq!(next.thermal.status, AdaptiveEnthalpyStatus::TrialLimit);
        };
        assert_eq!(
            bits(plain.specific_enthalpy_j_kg()),
            bits(full.history().specific_enthalpy_j_kg())
        );
        assert_eq!(
            bits(resumed.history().specific_enthalpy_j_kg()),
            bits(full.history().specific_enthalpy_j_kg())
        );
        assert_eq!(final_report.entropy, audited.entropy);
        assert_eq!(resumed.history().endpoint_solves(), full.history().endpoint_solves());
        assert_eq!(
            bits(&resumed.history().last_sdirk2_step().unwrap().first_stage.temperature),
            bits(&full.history().last_sdirk2_step().unwrap().first_stage.temperature)
        );
        // Independently, each reference mass is 1/12 kg and Cp is 5 J/(kg K).
        let expected: f64 = full.history().temperature().iter().zip(initial)
            .map(|(&t, h)| {
                let old_t = 280.0 + h / 5.0;
                (5.0 / 12.0) * ((t - old_t) / old_t).ln_1p()
            })
            .sum();
        close(audited.entropy.thermal_storage_change_j_k, expected, 1e-11);
        assert!(audited.entropy.conduction_generation_j_k > 0.0);
        assert!(!audited.entropy.has_unclassified_inputs);
        close(audited.entropy.balance_residual_j_k, 0.0, 1e-9);
        assert!(full.history().last_endpoint().is_none());
    });
}

#[test]
fn accepted_bath_entropy_uses_both_stage_temperatures_and_original_storage() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "bath", |face| face.centroid[0] == 0.0,
                ThermalBc::dirichlet(320.0).unwrap(),
            ).unwrap()
            .adiabatic_remainder().finish().unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = sensible();
        let storage = prepare(cx, &mesh, &curve);
        let initial = [200.0, 400.0, 200.0, 200.0];
        let mut state = storage.adaptive_state_with_entropy(
            cx, &initial, 0.0, 0.5, config(1e-3), 1e-9,
        ).unwrap();
        let p = input(problem(&mesh, &boundary, &material, &source));
        let mut bath = 0.0;
        let mut generation = 0.0;
        let mut endpoint_only = 0.0;
        let mut observations = 0;
        let report = state.advance_to_sdirk2_observed(cx, p, 1.0, 1000, |trial| {
            observations += 1;
            for half in [trial.first_half, trial.second_half] {
                // One free vertex has conductance k/6=0.5 W/K to the bath.
                // The physical material/bath production is Q*(1/T-1/Tbath).
                let t1 = half.first_stage.temperature[1];
                let t2 = half.second_stage.temperature[1];
                let q1 = -0.5 * (t1 - 320.0);
                let q2 = -0.5 * (t2 - 320.0);
                let first = (1.0 - SDIRK2_GAMMA) * half.duration_s;
                let second = SDIRK2_GAMMA * half.duration_s;
                bath += (first * q1 + second * q2) / 320.0;
                generation += first * q1 * (1.0 / t1 - 1.0 / 320.0)
                    + second * q2 * (1.0 / t2 - 1.0 / 320.0);
                endpoint_only += half.duration_s * q2 / 320.0;
            }
            ControlFlow::Continue(())
        }).unwrap();
        assert_eq!(report.thermal.status, AdaptiveEnthalpyStatus::ReachedEnd);
        assert_eq!(observations, report.thermal.accepted);
        let exact_storage = (5.0 / 12.0) * (state.history().temperature()[1] / 360.0).ln();
        close(report.entropy.thermal_storage_change_j_k, exact_storage, 1e-11);
        close(report.entropy.dirichlet_entropy_input_j_k, bath, 1e-10);
        close(report.entropy.conduction_generation_j_k, generation, 1e-10);
        close(report.entropy.temporal_defect_j_k, exact_storage - bath - generation, 1e-10);
        close(report.entropy.balance_residual_j_k, 0.0, 1e-9);
        close(report.thermal.integrated_net_input_j / 320.0, bath, 1e-10);
        assert!((bath - endpoint_only).abs() > 1e-5);
        assert!(generation > 0.0);
    });
}

#[test]
fn latent_entropy_and_source_classification_survive_method_switches() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder().finish().unwrap();
        let material = ConductivityModel::isotropic_declared(1.0).unwrap();
        let curve = latent();
        let storage = prepare(cx, &mesh, &curve);
        let mut state = storage.adaptive_state_with_entropy(
            cx, &[200.0; 4], 0.0, 0.5, config(1e-3), 1e-9,
        ).unwrap();
        let heater = ScalarField::Uniform(6.0);
        let heating = state.advance_to_sdirk2(
            cx, input(problem(&mesh, &boundary, &material, &heater)), 1.0, 100,
        ).unwrap();
        // V=1/6 m³ and 6 W/m³ supply one joule at Tm=350 K.
        close(heating.entropy.thermal_storage_change_j_k, 1.0 / 350.0, 1e-11);
        close(heating.entropy.source_at_node_j_k, 1.0 / 350.0, 1e-11);
        assert!(heating.entropy.has_unclassified_inputs);
        assert!(state.history().temperature().iter().all(|&t| t == 350.0));
        let sink = ScalarField::Uniform(-6.0);
        let cooling = state.advance_to(
            cx, problem(&mesh, &boundary, &material, &sink), None, 2.0, 100,
        ).unwrap();
        close(cooling.entropy.thermal_storage_change_j_k, 0.0, 1e-11);
        close(cooling.entropy.source_at_node_j_k, 0.0, 1e-11);
        assert!(cooling.entropy.has_unclassified_inputs);
        assert!(state.history().last_endpoint().is_some());
        assert!(state.history().last_sdirk2_step().is_none());

        let again = state.advance_to_sdirk2(
            cx, input(problem(&mesh, &boundary, &material, &heater)), 3.0, 100,
        ).unwrap();
        close(again.entropy.thermal_storage_change_j_k, 1.0 / 350.0, 1e-11);
        close(again.entropy.source_at_node_j_k, 1.0 / 350.0, 1e-11);
        close(again.entropy.balance_residual_j_k, 0.0, 1e-9);
        assert!(again.entropy.has_unclassified_inputs);
        assert!(state.history().last_sdirk2_step().is_some());
        for &h in state.history().specific_enthalpy_j_kg() {
            close(h, 203.0, 1e-8);
        }
    });
}

#[test]
fn negative_native_stage_production_refuses_without_accepting_the_pair() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder().finish().unwrap();
        let c = 0.9999;
        let material = ConductivityModel::constant_tensor(
            [[1.0, c, 0.0], [c, 1.0, 0.0], [0.0, 0.0, 1.0]],
        ).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = sensible();
        let storage = prepare(cx, &mesh, &curve);
        let middle = (300.0_f64 * 480.0).sqrt();
        let initial = [middle, 300.0, 480.0, middle].map(|t| 5.0 * (t - 280.0));
        let p = input(problem(&mesh, &boundary, &material, &source));
        let policy = config(1.0);
        let mut plain = storage.adaptive_state(cx, &initial, 0.0, 0.001, policy).unwrap();
        let thermal = plain.advance_to_sdirk2(cx, p, 0.001, 100).unwrap();
        assert_eq!(thermal.status, AdaptiveEnthalpyStatus::ReachedEnd);
        let mut audited = storage.adaptive_state_with_entropy(
            cx, &initial, 0.0, 0.001, policy, 1e-10,
        ).unwrap();
        assert!(matches!(
            audited.advance_to_sdirk2(cx, p, 0.001, 100),
            Err(EnthalpyEntropyError::NegativeProduction {
                account: "SDIRK2 first-stage bulk/surface/contact", ..
            })
        ));
        assert_eq!(audited.history().time_s(), 0.0);
        assert_eq!(audited.history().accepted_steps(), 0);
        assert_eq!(bits(audited.history().specific_enthalpy_j_kg()), bits(&initial));
        assert_eq!(audited.entropy_balance(), EnthalpyEntropyBalance::default());
    });
}

#[test]
fn sdirk_observer_cancellation_keeps_one_complete_energy_and_entropy_pair() {
    support::with_gate(|gate, cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder().finish().unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = sensible();
        let storage = prepare(cx, &mesh, &curve);
        let p = input(problem(&mesh, &boundary, &material, &source));
        let initial = [200.0, 250.0, 150.0, 200.0];
        let mut state = storage.adaptive_state_with_entropy(
            cx, &initial, 0.0, 0.1, config(10.0), 1e-9,
        ).unwrap();
        assert!(matches!(
            state.advance_to_sdirk2_observed(cx, p, 1.0, 100, |_| {
                gate.request();
                ControlFlow::Continue(())
            }),
            Err(EnthalpyEntropyError::Physical(_))
        ));
        assert_eq!(state.history().accepted_steps(), 1);
        assert!(state.history().time_s() > 0.0);
        assert!(state.entropy_balance().thermal_storage_change_j_k > 0.0);
        assert!(state.history().last_sdirk2_step().is_some());
        let saved_h = bits(state.history().specific_enthalpy_j_kg());
        let saved_entropy = state.entropy_balance();
        assert!(state.advance_to_sdirk2(cx, p, 1.0, 100).is_err());
        assert_eq!(bits(state.history().specific_enthalpy_j_kg()), saved_h);
        assert_eq!(state.entropy_balance(), saved_entropy);
    });
}
