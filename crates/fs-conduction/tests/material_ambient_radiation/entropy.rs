//! Original source-dependent emissivity and latent entropy share one heat load.
use super::*;
use fs_conduction::transient::enthalpy::entropy::EnthalpyEntropyError;
use fs_conduction::transient::enthalpy::EnthalpyError;

#[test]
fn material_emissivity_feedback_changes_latent_entropy_without_a_second_heat_owner() {
    let mesh = mesh();
    let boundary = boundary(&mesh);
    let card = card(Some("melt"));
    let context = context();
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let zero = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        source: &zero, element_materials: None };
    for datum in [0.0, 4096.0] {
        let curve = curve(&card, datum);
        for bath in [300.0_f64, 500.0] {
            let old = [datum + 500.0; 4];
            let dt = 0.2;
            // The exact implicit scalar phase-feedback equation on three
            // one-kilogram face vertices; apex h and every T stay unchanged.
            let coefficient = AREA * STEFAN_BOLTZMANN_W_M2_K4
                * (400.0_f64.powi(4) - bath.powi(4));
            let h = (1500.0 - dt * coefficient * 0.2)
                / (3.0 + dt * coefficient * 0.6 / LATENT);
            let heat = coefficient * (0.2 + 0.6 * h / LATENT);
            let expected_storage = 3.0 * (h - 500.0) / 400.0;
            with_cx(|cx| {
                let radiation = MaterialPointwiseAmbientRadiation::new(
                    cx, &mesh, &boundary, &[patch(&card, &curve, &context, bath)],
                ).unwrap();
                let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh,
                    &[ReferenceEnthalpyMaterial { curve: &curve, reference_density_kg_m3: RHO }],
                    &[0], budget(),
                ).unwrap();
                let original = engine.advance_with_radiation(
                    cx, problem, None, &old, dt, config(),
                    EnthalpyRadiation::MaterialAmbient(&radiation),
                ).unwrap();
                let (step, entropy) = engine.advance_with_radiation_and_entropy(
                    cx, problem, None, &old, dt, config(),
                    EnthalpyRadiation::MaterialAmbient(&radiation), 1e-10,
                ).unwrap();
                assert_eq!(step.specific_enthalpy_j_kg,
                    original.conduction().specific_enthalpy_j_kg);
                assert_eq!(step.temperature, [400.0; 4]);
                for &next in &step.specific_enthalpy_j_kg[..3] {
                    close(next - datum, h, 2e-9);
                }
                close(step.specific_enthalpy_j_kg[3], old[3], 2e-9);
                close(entropy.thermal_storage_change_j_k, expected_storage, 2e-11);
                close(entropy.radiation_entropy_input_j_k, -dt * heat / bath, 2e-11);
                close(entropy.radiation_generation_j_k,
                    dt * heat * (1.0 / bath - 1.0 / 400.0), 2e-11);
                close(entropy.temporal_defect_j_k, 0.0, 2e-11);
                close(entropy.algebraic_residual_j_k, 0.0, 2e-11);
                close(entropy.balance_residual_j_k, 0.0, 2e-11);
                close(engine.thermal_entropy_change_j_k(
                    cx, &old, &step.specific_enthalpy_j_kg,
                ).unwrap(), expected_storage, 2e-11);
                assert!(entropy.radiation_generation_j_k > 0.0);
                assert!(!entropy.has_unclassified_inputs);
                let report = &step.pointwise_radiation.as_ref().unwrap().patches[0];
                assert_eq!(report.material_states.len(), 3);
                for state in &report.material_states {
                    close(state.emissivity.value(), 0.2 + 0.6 * h / LATENT, 3e-12);
                    assert_eq!(state.emissivity.card_identity(), card.content_hash());
                }
                let (endpoint, audit) = engine.evaluate_endpoint_with_radiation_and_entropy(
                    cx, problem, None, &old, &step.specific_enthalpy_j_kg,
                    dt, config(), EnthalpyRadiation::MaterialAmbient(&radiation),
                ).unwrap();
                close(endpoint.energy_residual_j, 0.0, 1e-8);
                close(audit.radiation_entropy_input_j_k,
                    entropy.radiation_entropy_input_j_k, 1e-12);
                close(audit.radiation_generation_j_k, entropy.radiation_generation_j_k, 1e-12);
            });
            assert_eq!(old, [datum + 500.0; 4]);
        }
    }
}

#[test]
fn entropy_does_not_rebind_source_emissivity_to_an_equivalent_but_different_chart() {
    let mesh = mesh();
    let boundary = boundary(&mesh);
    let card = card(Some("melt"));
    let context = context();
    let curve = curve(&card, 0.0);
    let copied = curve.clone();
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let zero = ScalarField::Uniform(0.0);
    let old = [500.0; 4];
    with_cx(|cx| {
        let radiation = MaterialPointwiseAmbientRadiation::new(
            cx, &mesh, &boundary, &[patch(&card, &copied, &context, 300.0)],
        ).unwrap();
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            source: &zero, element_materials: None };
        assert!(matches!(engine.advance_with_radiation_and_entropy(
            cx, problem, None, &old, 0.2, config(),
            EnthalpyRadiation::MaterialAmbient(&radiation), 1e-10,
        ), Err(EnthalpyEntropyError::Physical(EnthalpyError::InvalidInput(_)))));
        assert!(matches!(engine.evaluate_endpoint_with_radiation_and_entropy(
            cx, problem, None, &old, &old, 0.2, config(),
            EnthalpyRadiation::MaterialAmbient(&radiation),
        ), Err(EnthalpyError::InvalidInput(_))));
    });
    assert_eq!(old, [500.0; 4]);
}
