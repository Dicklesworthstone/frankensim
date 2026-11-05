//! Physical entropy on authentic SDIRK stages, against independent equations.
use super::*;
use fs_conduction::transient::enthalpy::entropy::EnthalpyEntropyError;

#[test]
fn uniform_robin_sdirk_entropy_uses_native_weights_and_signed_temporal_defect() {
    let mesh = tetrahedron();
    let chart = sensible_curve();
    let conductivity = ConductivityModel::isotropic_declared(3.0).unwrap();
    let zero = ScalarField::Uniform(0.0);
    let conductance = 2.0;
    let mut builder = ThermalBoundaryBuilder::new(&mesh);
    for (face, name) in mesh
        .boundary()
        .iter()
        .zip(["bath0", "bath1", "bath2", "bath3"])
    {
        let vertices = face.vertices;
        builder = builder
            .region(
                name,
                |f| f.vertices == vertices,
                ThermalBc::robin(conductance / (4.0 * face.area), 300.0).unwrap(),
            )
            .unwrap();
    }
    let boundary = builder.finish().unwrap();
    with_cx(|cx| {
        // Each face has g/4 integrated conductance; each lumped node has C/4.
        // Uniform states therefore satisfy one independent C*T'=g*(Tb-T).
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &chart, 12.0, budget()).unwrap();
        let p = input(problem(&mesh, &boundary, &conductivity, &zero));
        let initial = [1100.0; 4]; // 500 K.
        let saved = bits(&initial);
        let original = engine
            .advance_sdirk2(cx, [p; 2], &initial, 1.0, config())
            .unwrap();
        let (step, entropy) = engine
            .advance_sdirk2_with_entropy(cx, [p; 2], &initial, 1.0, config(), 1e-10)
            .unwrap();
        assert_eq!(
            bits(&step.first_stage.specific_enthalpy_j_kg),
            bits(&original.first_stage.specific_enthalpy_j_kg)
        );
        assert_eq!(
            bits(&step.second_stage.specific_enthalpy_j_kg),
            bits(&original.second_stage.specific_enthalpy_j_kg)
        );
        let capacity = 10.0;
        let a = SDIRK2_GAMMA;
        let w = 1.0 - a;
        let first = (capacity * 500.0 + a * conductance * 300.0) / (capacity + a * conductance);
        let carry = w * conductance * (300.0 - first);
        let second =
            (capacity * 500.0 + carry + a * conductance * 300.0) / (capacity + a * conductance);
        let q1 = conductance * (300.0 - first);
        let q2 = conductance * (300.0 - second);
        let state = capacity * (second / 500.0).ln();
        let bath = (w * q1 + a * q2) / 300.0;
        let generation =
            w * q1 * (1.0 / first - 1.0 / 300.0) + a * q2 * (1.0 / second - 1.0 / 300.0);
        let temporal = state - w * q1 / first - a * q2 / second;
        for &t in &step.first_stage.temperature {
            close(t, first, 1e-8);
        }
        for &t in &step.second_stage.temperature {
            close(t, second, 1e-8);
        }
        close(entropy.thermal_storage_change_j_k, state, 1e-10);
        close(entropy.robin_entropy_input_j_k, bath, 1e-10);
        close(entropy.robin_generation_j_k, generation, 1e-10);
        close(entropy.temporal_defect_j_k, temporal, 1e-10);
        assert!(entropy.temporal_defect_j_k < 0.0);
        assert!(entropy.robin_generation_j_k > 0.0);
        close(entropy.conduction_generation_j_k, 0.0, 1e-11);
        close(entropy.algebraic_residual_j_k, 0.0, 1e-10);
        close(entropy.balance_residual_j_k, 0.0, 1e-11);
        close(
            step.stored_energy_change_j,
            capacity * (second - 500.0),
            1e-9,
        );
        assert!(!entropy.has_unclassified_inputs);
        assert_eq!(bits(&initial), saved);
    });
}

#[test]
fn sdirk_latent_storage_and_film_entropy_use_original_history_and_source_quadrature() {
    use fs_conduction::transient::backward_euler::{ArealHeatCapacity, SurfaceHeatCapacity};
    let mesh = tetrahedron();
    let boundary = adiabatic(&mesh);
    let conductivity = ConductivityModel::isotropic_declared(3.0).unwrap();
    let sources = [ScalarField::Uniform(100.0), ScalarField::Uniform(-50.0)];
    for datum in [0.0, 1_048_576.0] {
        let chart = curve(&[
            (datum, 300.0, 0.0, 4.0),
            (datum + 100.0, 350.0, 0.0, 3.0),
            (datum + 300.0, 350.0, 1.0, 2.0),
            (datum + 600.0, 500.0, 1.0, 1.0),
        ]);
        with_cx(|cx| {
            let film = SurfaceHeatCapacity::uniform(0, ArealHeatCapacity::declared(12.0).unwrap());
            let engine =
                EnthalpyBackwardEuler::uniform(cx, &mesh, &chart, REFERENCE_DENSITY, budget())
                    .unwrap()
                    .with_surface_heat_capacity(cx, &[film])
                    .unwrap();
            let inputs =
                [0, 1].map(|i| input(problem(&mesh, &boundary, &conductivity, &sources[i])));
            let old = [datum + 150.0; 4];
            let saved = bits(&old);
            let mut policy = config();
            policy.newton.absolute_tolerance = 1e-10;
            policy.newton.relative_tolerance = 0.0;
            let (step, entropy) = engine
                .advance_sdirk2_with_entropy(cx, inputs, &old, 1.0, policy, 1e-9)
                .unwrap();
            let first_h = 150.0 + SDIRK2_GAMMA * 50.0;
            let final_h = 150.0 + (1.0 - SDIRK2_GAMMA) * 50.0 - SDIRK2_GAMMA * 25.0;
            let energy = VOLUME * ((1.0 - SDIRK2_GAMMA) * 100.0 - SDIRK2_GAMMA * 50.0);
            for i in 0..4 {
                close(
                    step.first_stage.specific_enthalpy_j_kg[i] - datum,
                    first_h,
                    1e-8,
                );
                close(
                    step.second_stage.specific_enthalpy_j_kg[i] - datum,
                    final_h,
                    1e-8,
                );
                assert_eq!(step.first_stage.temperature[i], 350.0);
                assert_eq!(step.second_stage.temperature[i], 350.0);
            }
            close(step.stored_energy_change_j, energy, 1e-9);
            close(entropy.thermal_storage_change_j_k, energy / 350.0, 1e-10);
            close(entropy.source_at_node_j_k, energy / 350.0, 1e-11);
            close(entropy.conduction_generation_j_k, 0.0, 1e-12);
            close(entropy.temporal_defect_j_k, 0.0, 1e-10);
            close(entropy.algebraic_residual_j_k, 0.0, 1e-10);
            close(entropy.balance_residual_j_k, 0.0, 1e-11);
            assert!(entropy.has_unclassified_inputs);
            assert_eq!(bits(&old), saved);
        });
    }
}

#[test]
fn heterogeneous_sdirk_entropy_retains_distinct_calorics_and_frozen_masses() {
    use fs_conduction::transient::enthalpy::heterogeneous::{
        HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
    };
    let mesh = ConductionMesh::new(
        TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [2.0, 0.0, 1.0],
        ],
    )
    .unwrap();
    let first = sensible_curve();
    let second = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x75; 32]),
        SolidLiquidPhase::Solid,
        [(0.0, 300.0), (800.0, 500.0)]
            .map(|(h, t)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 9.0,
            })
            .to_vec(),
    )
    .unwrap();
    let materials = [
        ReferenceEnthalpyMaterial {
            curve: &first,
            reference_density_kg_m3: 2.0,
        },
        ReferenceEnthalpyMaterial {
            curve: &second,
            reference_density_kg_m3: 6.0,
        },
    ];
    let boundary = adiabatic(&mesh);
    let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(-12.0);
    with_cx(|cx| {
        let engine = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &materials,
            &[0, 1],
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        let p = input(problem(&mesh, &boundary, &conductivity, &source));
        let (step, entropy) = engine
            .advance_sdirk2_with_entropy(
                cx,
                [p; 2],
                &[150.0, 150.0, 150.0, 150.0, 240.0, 240.0, 240.0, 240.0],
                0.5,
                config(),
                1e-10,
            )
            .unwrap();
        let a = SDIRK2_GAMMA;
        let t1 = [310.0 - a * 0.6, 360.0 - a * 0.25];
        let t2 = [309.4_f64, 359.75_f64];
        let state = (2.0 * VOLUME * 5.0) * (t2[0] / 310.0).ln()
            + (6.0 * VOLUME * 4.0) * (t2[1] / 360.0).ln();
        let at_node = -(1.0 - a) * (1.0 / t1[0] + 1.0 / t1[1]) - a * (1.0 / t2[0] + 1.0 / t2[1]);
        for i in 0..8 {
            close(
                step.second_stage.temperature[i],
                t2[usize::from(i >= 4)],
                1e-9,
            );
            close(
                engine.reference_nodal_masses_kg()[i],
                if i < 4 { 1.0 / 12.0 } else { 0.25 },
                1e-15,
            );
        }
        close(step.stored_energy_change_j, -2.0, 1e-9);
        close(entropy.thermal_storage_change_j_k, state, 1e-11);
        close(entropy.source_at_node_j_k, at_node, 1e-11);
        close(entropy.temporal_defect_j_k, state - at_node, 1e-11);
        close(entropy.balance_residual_j_k, 0.0, 1e-11);
        assert!(entropy.has_unclassified_inputs);
    });
}

#[test]
fn sdirk_stage_entropy_refuses_negative_spatial_production_despite_spd_material() {
    let mesh = tetrahedron();
    let boundary = adiabatic(&mesh);
    let chart = sensible_curve();
    let c = 0.9999;
    let conductivity =
        ConductivityModel::constant_tensor([[1.0, c, 0.0], [c, 1.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
    let zero = ScalarField::Uniform(0.0);
    let middle = (300.0_f64 * 480.0).sqrt();
    let old = [middle, 300.0, 480.0, middle].map(|t| 5.0 * (t - 280.0));
    let saved = bits(&old);
    with_cx(|cx| {
        let engine =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &chart, REFERENCE_DENSITY, budget()).unwrap();
        let p = input(problem(&mesh, &boundary, &conductivity, &zero));
        let original = engine
            .advance_sdirk2(cx, [p; 2], &old, 0.001, config())
            .unwrap();
        close(original.energy_residual_j, 0.0, 1e-9);
        assert!(matches!(
            engine.advance_sdirk2_with_entropy(cx, [p; 2], &old, 0.001, config(), 1e-10,),
            Err(EnthalpyEntropyError::NegativeProduction {
                account: "SDIRK2 first-stage bulk/surface/contact",
                ..
            })
        ));
        assert_eq!(bits(&old), saved);
    });
}

#[test]
fn coarse_prescribed_sdirk_history_refuses_entropy_without_counting_carry_as_heat() {
    let mesh = tetrahedron();
    let chart = sensible_curve();
    let conductivity = ConductivityModel::isotropic_declared(3.0).unwrap();
    let zero = ScalarField::Uniform(0.0);
    let boundaries = [330.0, 360.0].map(|t| {
        ThermalBoundaryBuilder::new(&mesh)
            .region("stage bath", |_| true, ThermalBc::dirichlet(t).unwrap())
            .unwrap()
            .finish()
            .unwrap()
    });
    with_cx(|cx| {
        let engine =
            EnthalpyBackwardEuler::uniform(cx, &mesh, &chart, REFERENCE_DENSITY, budget()).unwrap();
        let inputs = [0, 1].map(|i| input(problem(&mesh, &boundaries[i], &conductivity, &zero)));
        let old = [100.0; 4]; // 300 K.
        let saved = bits(&old);
        let original = engine
            .advance_sdirk2(cx, inputs, &old, 1.0, config())
            .unwrap();
        close(
            original.stored_energy_change_j,
            REFERENCE_DENSITY * VOLUME * 5.0 * 60.0,
            1e-10,
        );
        let r = (1.0 - SDIRK2_GAMMA) / SDIRK2_GAMMA;
        let expected = REFERENCE_DENSITY
            * VOLUME
            * 5.0
            * ((360.0_f64 / 300.0).ln() - r * 30.0 / 330.0 - (60.0 - r * 30.0) / 360.0);
        match engine.advance_sdirk2_with_entropy(cx, inputs, &old, 1.0, config(), 1e-10) {
            Err(EnthalpyEntropyError::NegativeProduction {
                account: "thermal storage minus inputs",
                value_j_k,
                ..
            }) => close(value_j_k, expected, 1e-12),
            other => panic!("expected coarse-step entropy refusal, got {other:?}"),
        }
        assert!(expected < 0.0);
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                engine.advance_sdirk2_with_entropy(cancelled, inputs, &old, 1.0, config(), 1e-10),
                Err(EnthalpyEntropyError::Physical(EnthalpyError::Conduction(
                    ConductionError::Cancelled { .. }
                )))
            ))
        });
        assert!(matches!(
            engine.advance_sdirk2_with_entropy(cx, inputs, &old, 1.0, config(), f64::NAN,),
            Err(EnthalpyEntropyError::Physical(EnthalpyError::InvalidInput(
                _
            )))
        ));
        assert_eq!(bits(&old), saved);
    });
}
