//! G1/G3/G4: original adaptive fields and accepted-only thermodynamic history.
#[path = "entropy/sdirk.rs"]
mod sdirk;

use super::*;
use fs_conduction::transient::enthalpy::entropy::{EnthalpyEntropyBalance, EnthalpyEntropyError};

#[test]
fn entropy_follows_only_accepted_fine_pairs_and_survives_budget_resume() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = sensible();
        let storage = prepare(cx, &mesh, &curve);
        let p = problem(&mesh, &boundary, &material, &source);
        let initial = [200.0, 250.0, 150.0, 200.0];
        let policy = config(1e-3);
        let mut plain = storage.adaptive_state(cx, &initial, 0.0, 1.0, policy).unwrap();
        let mut full = storage.adaptive_state_with_entropy(
            cx, &initial, 0.0, 1.0, policy, 1e-9,
        ).unwrap();
        let mut resumed = full.clone();
        let rejected = resumed.advance_to(cx, p, None, 1.0, 1).unwrap();
        assert_eq!(rejected.thermal.accepted, 0);
        assert_eq!(rejected.thermal.rejected, 1);
        assert_eq!(rejected.entropy, EnthalpyEntropyBalance::default());
        assert_eq!(bits(resumed.history().specific_enthalpy_j_kg()), bits(&initial));
        let native = plain.advance_to(cx, p, None, 1.0, 10000).unwrap();
        let report = full.advance_to(cx, p, None, 1.0, 10000).unwrap();
        let final_report = loop {
            let next = resumed.advance_to(cx, p, None, 1.0, 1).unwrap();
            if next.thermal.status == AdaptiveEnthalpyStatus::ReachedEnd { break next; }
            assert_eq!(next.thermal.status, AdaptiveEnthalpyStatus::TrialLimit);
        };
        assert_eq!(native, report.thermal);
        assert_eq!(bits(plain.specific_enthalpy_j_kg()), bits(full.history().specific_enthalpy_j_kg()));
        assert_eq!(bits(resumed.history().specific_enthalpy_j_kg()), bits(full.history().specific_enthalpy_j_kg()));
        assert_eq!(final_report.entropy, report.entropy);
        assert_eq!(resumed.history().endpoint_solves(), full.history().endpoint_solves());
        // Independently known rho*V/4=1/12 kg and Cp=5 J/(kg K).
        let exact: f64 = full.history().temperature().iter().zip(initial)
            .map(|(&t, h)| (5.0/12.0) * ((t-(280.0+h/5.0))/(280.0+h/5.0)).ln_1p())
            .sum();
        assert!((report.entropy.thermal_storage_change_j_k-exact).abs() < 1e-11);
        assert!(report.entropy.conduction_generation_j_k > 0.0);
        assert!(report.entropy.temporal_defect_j_k > 0.0);
        assert!(!report.entropy.has_unclassified_inputs);
        assert!(report.entropy.balance_residual_j_k.abs() < 1e-9);
    });
}

#[test]
fn isothermal_latent_storage_and_cancelling_sources_keep_the_original_account() {
    support::with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let material = ConductivityModel::isotropic_declared(1.0).unwrap();
        let curve = latent();
        let storage = prepare(cx, &mesh, &curve);
        let mut state = storage.adaptive_state_with_entropy(
            cx, &[200.0;4], 0.0, 0.5, config(1e-3), 1e-9,
        ).unwrap();
        let source = ScalarField::Uniform(6.0);
        let heating = state.advance_to(
            cx, problem(&mesh, &boundary, &material, &source), None, 1.0, 100,
        ).unwrap();
        // Volume 1/6 m³ gives exactly one joule; latent entropy is Q/Tm.
        assert!((heating.entropy.thermal_storage_change_j_k-1.0/350.0).abs() < 1e-11);
        assert!((heating.entropy.source_at_node_j_k-1.0/350.0).abs() < 1e-11);
        assert!(heating.entropy.has_unclassified_inputs);
        assert!(state.history().temperature().iter().all(|&t| t == 350.0));
        let sink = ScalarField::Uniform(-6.0);
        let cooling = state.advance_to(
            cx, problem(&mesh, &boundary, &material, &sink), None, 2.0, 100,
        ).unwrap();
        assert!(cooling.entropy.thermal_storage_change_j_k.abs() < 1e-11);
        assert!(cooling.entropy.source_at_node_j_k.abs() < 1e-11);
        // The source classification survives signed cancellation across calls.
        assert!(cooling.entropy.has_unclassified_inputs);
        assert!(cooling.entropy.balance_residual_j_k.abs() < 1e-9);
        assert!(state.history().specific_enthalpy_j_kg().iter().all(|h| (h-200.0).abs() < 1e-8));
    });
}

#[test]
fn observer_cancellation_retains_one_complete_physical_and_entropy_pair() {
    support::with_gate(|gate, cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let material = ConductivityModel::isotropic_declared(3.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = sensible();
        let storage = prepare(cx, &mesh, &curve);
        let p = problem(&mesh, &boundary, &material, &source);
        let initial = [200.0,250.0,150.0,200.0];
        let mut state = storage.adaptive_state_with_entropy(
            cx, &initial, 0.0, 0.1, config(10.0), 1e-9,
        ).unwrap();
        assert!(matches!(
            state.advance_to_observed(cx, p, None, 1.0, 100, |_| {
                gate.request();
                ControlFlow::Continue(())
            }),
            Err(EnthalpyEntropyError::Physical(_))
        ));
        assert_eq!(state.history().accepted_steps(), 1);
        assert!(state.history().time_s() > 0.0);
        assert!(state.entropy_balance().thermal_storage_change_j_k > 0.0);
        let saved = bits(state.history().specific_enthalpy_j_kg());
        let entropy = state.entropy_balance();
        assert!(state.advance_to(cx, p, None, 1.0, 100).is_err());
        assert_eq!(bits(state.history().specific_enthalpy_j_kg()), saved);
        assert_eq!(state.entropy_balance(), entropy);
    });
}
