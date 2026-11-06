//! G1/G3: original material emissivity, latent entropy and explicit open reservoirs.
use super::*;
use fs_conduction::transient::enthalpy::entropy::EnthalpyEntropyBalance;

#[test]
fn latent_material_reflection_preserves_native_states_and_datum_invariant_entropy() {
    let mesh = mesh();
    let boundary = boundary(&mesh);
    let context = context();
    let cards = [material_card("left"), material_card("right")];
    let enclosure = closed(&mesh, &cards, &context);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    let dt = 5.0;
    let mut baseline: Option<EnthalpyEntropyBalance> = None;
    for datum in [0.0, 4096.0] {
        let curves = [
            phase_curve(&cards[0], 400.0, datum),
            phase_curve(&cards[1], 350.0, datum),
        ];
        let declarations = bindings(&cards, &curves, &context);
        let old = [datum + 500.0; 8];
        with_cx(|cx| {
            let base = PointwiseEnclosureRadiation::new(
                cx,
                &mesh,
                &boundary,
                &enclosure,
                PointwiseEnclosureTolerance::default(),
            )
            .unwrap();
            let owner = MaterialPointwiseEnclosureRadiation::new(cx, &base, &declarations).unwrap();
            let engine = HeterogeneousEnthalpyBackwardEuler::new(
                cx,
                &mesh,
                &[
                    ReferenceEnthalpyMaterial {
                        curve: &curves[0],
                        reference_density_kg_m3: RHO,
                    },
                    ReferenceEnthalpyMaterial {
                        curve: &curves[1],
                        reference_density_kg_m3: RHO,
                    },
                ],
                &[0, 1],
                budget(),
            )
            .unwrap();
            let native = engine
                .advance_with_radiation(
                    cx,
                    problem,
                    None,
                    &old,
                    dt,
                    config(),
                    EnthalpyRadiation::MaterialEnclosure(&owner),
                )
                .unwrap();
            let (step, entropy) = engine
                .advance_with_radiation_and_entropy(
                    cx,
                    problem,
                    None,
                    &old,
                    dt,
                    config(),
                    EnthalpyRadiation::MaterialEnclosure(&owner),
                    1e-9,
                )
                .unwrap();
            assert_eq!(
                step.specific_enthalpy_j_kg,
                native.conduction().specific_enthalpy_j_kg
            );
            assert_eq!(step.temperature, native.conduction().temperature);
            assert_eq!(&step.temperature[..4], &[400.0; 4]);
            assert_eq!(&step.temperature[4..], &[350.0; 4]);
            let final_h: [f64; 8] = step.specific_enthalpy_j_kg.as_slice().try_into().unwrap();
            let (_, heat, bath) = oracle(&curves, &final_h, 1.0, 300.0);
            assert_eq!(bath, 0.0);
            let report = step.pointwise_enclosure_radiation.as_ref().unwrap();
            assert!(report.ambient.is_none());
            assert!(report.patches[0].material_mean_emissivity.unwrap() < 0.5);
            assert!(report.patches[1].material_mean_emissivity.unwrap() > 0.5);
            close(report.outward_heat_w, 0.0, 1e-9);
            close(step.stored_energy_change_j, 0.0, 1e-8);
            assert_eq!(entropy.radiation_entropy_input_j_k, 0.0);
            let generation = dt * heat[0] * (1.0 / 350.0 - 1.0 / 400.0);
            // Each original nodal mass is exactly 1 kg. All storage changes
            // remain latent; there is no sensible heat or added phase source.
            let storage: f64 = final_h
                .iter()
                .zip(old)
                .enumerate()
                .map(|(i, (h, old))| (h - old) / if i < 4 { 400.0 } else { 350.0 })
                .sum();
            close(entropy.radiation_generation_j_k, generation, 2e-10);
            close(entropy.thermal_storage_change_j_k, storage, 2e-10);
            close(entropy.thermal_storage_change_j_k, generation, 1e-9);
            close(entropy.temporal_defect_j_k, 0.0, 2e-10);
            close(entropy.balance_residual_j_k, 0.0, 1e-9);
            assert!(entropy.radiation_generation_j_k > 0.0);
            assert!(!entropy.has_unclassified_inputs);
            if let Some(first) = baseline {
                close(
                    entropy.thermal_storage_change_j_k,
                    first.thermal_storage_change_j_k,
                    1e-9,
                );
                close(
                    entropy.radiation_generation_j_k,
                    first.radiation_generation_j_k,
                    1e-9,
                );
            } else {
                baseline = Some(entropy);
            }
        });
        assert_eq!(old, [datum + 500.0; 8]);
    }
}

#[test]
fn nonuniform_material_enclosure_endpoint_uses_the_actual_open_bath_in_both_heat_directions() {
    let mesh = mesh();
    let boundary = boundary(&mesh);
    let context = context();
    let cards = [material_card("left"), material_card("right")];
    let curves = [
        phase_curve(&cards[0], 400.0, 0.0),
        phase_curve(&cards[1], 350.0, 0.0),
    ];
    let declarations = bindings(&cards, &curves, &context);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    let old = [500.0; 8];
    let candidate = [-30.0, 500.0, 1100.0, 500.0, -60.0, 200.0, 1050.0, 500.0];
    let dt = 0.2;
    with_cx(|cx| {
        let open = OpenDiffuseEnclosure::new(
            cx,
            surfaces(&mesh, &cards, &context),
            vec![vec![0.0, 0.4], vec![0.4, 0.0]],
            ViewFactorEvidence::Analytic {
                geometry: "reciprocal material entropy fixture".into(),
            },
            ViewFactorTolerance::default(),
        )
        .unwrap();
        let engine = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[
                ReferenceEnthalpyMaterial {
                    curve: &curves[0],
                    reference_density_kg_m3: RHO,
                },
                ReferenceEnthalpyMaterial {
                    curve: &curves[1],
                    reference_density_kg_m3: RHO,
                },
            ],
            &[0, 1],
            budget(),
        )
        .unwrap();
        for bath_temperature in [280.0, 520.0] {
            let base = PointwiseEnclosureRadiation::new_open(
                cx,
                &mesh,
                &boundary,
                &open,
                bath_temperature,
                PointwiseEnclosureTolerance::default(),
            )
            .unwrap();
            let owner = MaterialPointwiseEnclosureRadiation::new(cx, &base, &declarations).unwrap();
            let (endpoint, entropy) = engine
                .evaluate_endpoint_with_radiation_and_entropy(
                    cx,
                    problem,
                    None,
                    &old,
                    &candidate,
                    dt,
                    config(),
                    EnthalpyRadiation::MaterialEnclosure(&owner),
                )
                .unwrap();
            let (load, heat, bath) = oracle(&curves, &candidate, 0.4, bath_temperature);
            assert_eq!(bath > 0.0, bath_temperature == 520.0);
            let surface_entropy_rate: f64 = load
                .iter()
                .zip(&endpoint.temperature)
                .map(|(q, t)| q / t)
                .sum();
            close(endpoint.radiation_out_w, heat.iter().sum(), 1e-9);
            close(
                entropy.radiation_entropy_input_j_k,
                dt * bath / bath_temperature,
                2e-11,
            );
            close(
                entropy.radiation_generation_j_k,
                -dt * (surface_entropy_rate + bath / bath_temperature),
                2e-11,
            );
            close(entropy.balance_residual_j_k, 0.0, 1e-9);
            assert!(endpoint.residual_norm_j > 1.0);
            assert!(entropy.algebraic_residual_j_k.abs() > 0.1);
            assert!(!entropy.has_unclassified_inputs);
        }
    });
    assert_eq!(old, [500.0; 8]);
    assert_eq!(
        candidate,
        [-30.0, 500.0, 1100.0, 500.0, -60.0, 200.0, 1050.0, 500.0]
    );
}
