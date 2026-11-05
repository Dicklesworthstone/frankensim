//! Closed material/field/black-bath entropy on the native pointwise heat owner.
use super::*;
use fs_conduction::transient::enthalpy::entropy::EnthalpyEntropyError;
use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use fs_conduction::transient::enthalpy::sdirk::{EnthalpySdirk2Input, SDIRK2_GAMMA};

fn full_surface_rates(mesh: &ConductionMesh, t: &[f64], bath: f64) -> (f64, f64) {
    let mut heat = 0.0;
    let mut received = 0.0;
    for face in mesh.boundary() {
        let q = exact_face_load(face.vertices.map(|v| t[v as usize]), face.area, bath);
        for (i, &v) in face.vertices.iter().enumerate() {
            heat += q[i];
            received += q[i] / t[v as usize];
        }
    }
    (-heat / bath, heat / bath - received)
}

#[test]
fn pointwise_bath_entropy_retains_the_native_solution_in_both_heat_directions() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .remainder("exterior", ThermalBc::adiabatic()).unwrap().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let zero = ScalarField::Uniform(0.0);
    let curve = phase_curve(false);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        source: &zero, element_materials: None };
    let old = [CP * (400.0 - 200.0); 4];
    let dt = 0.05;
    let capacity = RHO * CP * 2.0_f64.sqrt() / 12.0;
    with_cx(|cx| {
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        for bath in [300.0, 500.0] {
            let patches = [patch(bath, 600.0)];
            let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
            let original = engine.advance_with_pointwise_ambient_radiation(
                cx, problem, None, &old, dt, config(), &radiation,
            ).unwrap();
            let (step, entropy) = engine.advance_with_radiation_and_entropy(
                cx, problem, None, &old, dt, config(),
                EnthalpyRadiation::PointwiseAmbient(&radiation), 1e-10,
            ).unwrap();
            assert_eq!(step.specific_enthalpy_j_kg, original.specific_enthalpy_j_kg);
            assert_eq!(step.newton, original.newton);
            let t = step.temperature[0];
            let (input, generation) = full_surface_rates(&mesh, &step.temperature, bath);
            let storage = capacity * (t / 400.0).ln();
            close(entropy.thermal_storage_change_j_k, storage, 1e-11);
            close(entropy.radiation_entropy_input_j_k, dt * input, 1e-11);
            close(entropy.radiation_generation_j_k, dt * generation, 1e-11);
            close(entropy.temporal_defect_j_k,
                storage - capacity * (t - 400.0) / t, 1e-11);
            close(entropy.algebraic_residual_j_k, 0.0, 1e-10);
            close(entropy.balance_residual_j_k, 0.0, 1e-11);
            close(engine.thermal_entropy_change_j_k(
                cx, &old, &step.specific_enthalpy_j_kg,
            ).unwrap(), storage, 1e-11);
            assert!(entropy.radiation_generation_j_k > 0.0);
            assert_eq!(entropy.robin_entropy_input_j_k, 0.0);
            assert_eq!(entropy.robin_generation_j_k, 0.0);
            assert!(!entropy.has_unclassified_inputs);
        }
    });
    assert_eq!(old, [2000.0; 4]);
}

#[test]
fn independent_baths_and_fixed_node_reactions_use_actual_nodal_radiation() {
    let mesh = regular_tet();
    let warm_face = mesh.boundary().iter().find(|face| face.centroid[2] > 0.0).unwrap().vertices;
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("fixed", |face| face.centroid[2] == 0.0, ThermalBc::dirichlet(400.0).unwrap())
        .unwrap()
        .region("warm", |face| face.vertices == warm_face, ThermalBc::adiabatic()).unwrap()
        .remainder("exterior", ThermalBc::adiabatic()).unwrap().finish().unwrap();
    let cold = patch(300.0, 600.0);
    let warm = AmbientRadiationPatch::new("warm", cold.emissivity().clone(), 600.0).unwrap();
    let patches = [cold, warm];
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let zero = ScalarField::Uniform(0.0);
    let curve = phase_curve(false);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        source: &zero, element_materials: None };
    let t = [400.0, 400.0, 400.0, 500.0];
    let h = t.map(|t| CP * (t - 200.0));
    let old = [1000.0; 4];
    let dt = 0.2;
    let mut nodal_heat = [0.0; 4];
    let mut bath_rate = 0.0;
    for face in mesh.boundary() {
        if face.centroid[2] == 0.0 { continue; }
        let bath = if face.vertices == warm_face { 600.0 } else { 300.0 };
        let q = exact_face_load(face.vertices.map(|v| t[v as usize]), face.area, bath);
        for (i, &v) in face.vertices.iter().enumerate() {
            nodal_heat[v as usize] += q[i];
            bath_rate += q[i] / bath;
        }
    }
    with_cx(|cx| {
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let (endpoint, entropy) = engine.evaluate_endpoint_with_radiation_and_entropy(
            cx, problem, None, &old, &h, dt, config(),
            EnthalpyRadiation::PointwiseAmbient(&radiation),
        ).unwrap();
        close(entropy.radiation_entropy_input_j_k, -dt * bath_rate, 2e-11);
        let at_node: f64 = nodal_heat.iter().zip(t).map(|(q, t)| q / t).sum();
        close(entropy.radiation_generation_j_k, dt * (bath_rate - at_node), 2e-11);
        close(entropy.dirichlet_entropy_input_j_k, dt * endpoint.dirichlet_in_w / 400.0, 2e-11);
        close(entropy.algebraic_residual_j_k,
            endpoint.residual_j.iter().zip(t).map(|(r, t)| r / t).sum(), 2e-11);
        close(entropy.balance_residual_j_k, 0.0, 2e-11);
        assert!(endpoint.residual_norm_j > 1.0, "supplied candidates are measured, not solved");
        assert_eq!(endpoint.free_dofs, 1);
    });
    assert_eq!(old, [1000.0; 4]);
    assert_eq!(h, [2000.0, 2000.0, 2000.0, 3000.0]);
}

#[test]
fn sdirk2_weights_both_actual_radiative_stage_baths_once() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .remainder("exterior", ThermalBc::adiabatic()).unwrap().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let zero = ScalarField::Uniform(0.0);
    let curve = phase_curve(false);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        source: &zero, element_materials: None };
    let patches1 = [patch(300.0, 600.0)];
    let patches2 = [patch(500.0, 600.0)];
    let old = [CP * (450.0 - 200.0); 4];
    let dt = 0.05;
    with_cx(|cx| {
        let first = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches1).unwrap();
        let second = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches2).unwrap();
        let inputs = [&first, &second].map(|radiation| EnthalpySdirk2Input {
            problem, interfaces: None, radiation: EnthalpyRadiation::PointwiseAmbient(radiation),
            contacts: None,
        });
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let original = engine.advance_sdirk2(cx, inputs, &old, dt, config()).unwrap();
        let (step, entropy) = engine.advance_sdirk2_with_entropy(
            cx, inputs, &old, dt, config(), 1e-10,
        ).unwrap();
        assert_eq!(step.first_stage.specific_enthalpy_j_kg,
            original.first_stage.specific_enthalpy_j_kg);
        assert_eq!(step.second_stage.specific_enthalpy_j_kg,
            original.second_stage.specific_enthalpy_j_kg);
        let a = full_surface_rates(&mesh, &step.first_stage.temperature, 300.0);
        let b = full_surface_rates(&mesh, &step.second_stage.temperature, 500.0);
        close(entropy.radiation_entropy_input_j_k,
            dt * ((1.0 - SDIRK2_GAMMA) * a.0 + SDIRK2_GAMMA * b.0), 2e-11);
        close(entropy.radiation_generation_j_k,
            dt * ((1.0 - SDIRK2_GAMMA) * a.1 + SDIRK2_GAMMA * b.1), 2e-11);
        let capacity = RHO * CP * 2.0_f64.sqrt() / 12.0;
        let ds = capacity * (step.second_stage.temperature[0] / 450.0).ln();
        close(entropy.thermal_storage_change_j_k, ds, 2e-11);
        close(entropy.balance_residual_j_k, 0.0, 2e-11);
        close(entropy.algebraic_residual_j_k, 0.0, 1e-10);
        assert_eq!(entropy.source_at_node_j_k, 0.0);
        assert_eq!(entropy.dirichlet_entropy_input_j_k, 0.0);
        assert!(entropy.radiation_generation_j_k > 0.0);
    });
    assert_eq!(old, [2500.0; 4]);
}

#[test]
fn signed_consistent_radiation_production_is_measured_and_refused() {
    let mesh = regular_tet();
    let face = &mesh.boundary()[0];
    let selected = face.vertices;
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("exterior", |face| face.vertices == selected, ThermalBc::adiabatic()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(1e-6).unwrap();
    let zero = ScalarField::Uniform(0.0);
    let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x82; 32]), SolidLiquidPhase::Solid,
        [(0.0, 200.0), (68000.0, 7000.0)].map(|(h, t)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h, temperature_k: t,
            liquid_mass_fraction: 0.0, bulk_density_kg_m3: RHO,
        }).to_vec(),
    ).unwrap();
    let patches = [patch(1000.0, 7000.0)];
    let mut t = [300.0; 4];
    t[selected[0] as usize] = 6000.0;
    let h = t.map(|t| CP * (t - 200.0));
    let dt = 1e-4;
    let q = exact_face_load([6000.0, 300.0, 300.0], face.area, 1000.0);
    let generation: f64 = q.iter().zip([6000.0, 300.0, 300.0])
        .map(|(q, t)| q * (1.0 / 1000.0 - 1.0 / t)).sum();
    assert!(generation < -100.0);
    with_cx(|cx| {
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            source: &zero, element_materials: None };
        let (_, entropy) = engine.evaluate_endpoint_with_radiation_and_entropy(
            cx, problem, None, &h, &h, dt, config(),
            EnthalpyRadiation::PointwiseAmbient(&radiation),
        ).unwrap();
        close(entropy.radiation_generation_j_k, dt * generation, 1e-12);
        close(entropy.balance_residual_j_k, 0.0, 1e-12);
        assert!(matches!(entropy.admit(1e-10),
            Err(EnthalpyEntropyError::NegativeProduction {
                account: "radiation material/field/bath", ..
            })));
    });
}

#[test]
fn near_equilibrium_radiation_is_not_replaced_by_photon_or_mean_surface_entropy() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .remainder("exterior", ThermalBc::adiabatic()).unwrap().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let zero = ScalarField::Uniform(0.0);
    let curve = phase_curve(false);
    let bath = 400.0;
    let t = bath + 1e-6;
    let h = [CP * (t - 200.0); 4];
    let patches = [patch(bath, 600.0)];
    with_cx(|cx| {
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            source: &zero, element_materials: None };
        let (endpoint, entropy) = engine.evaluate_endpoint_with_radiation_and_entropy(
            cx, problem, None, &h, &h, 1.0, config(),
            EnthalpyRadiation::PointwiseAmbient(&radiation),
        ).unwrap();
        let actual_t = endpoint.temperature[0];
        let q = EPSILON * STEFAN_BOLTZMANN_W_M2_K4 * 3.0_f64.sqrt()
            * (actual_t - bath) * (actual_t + bath) * (actual_t * actual_t + bath * bath);
        close(entropy.radiation_entropy_input_j_k, -q / bath, 2e-22);
        close(entropy.radiation_generation_j_k,
            (q / actual_t) * ((actual_t - bath) / bath), 2e-22);
        assert!(entropy.radiation_generation_j_k > 0.0);
        close(entropy.balance_residual_j_k, 0.0, 1e-14);
    });
}

#[test]
fn unqualified_models_invalid_support_and_cancelled_entropy_leave_history_intact() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .remainder("exterior", ThermalBc::adiabatic()).unwrap().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let zero = ScalarField::Uniform(0.0);
    let curve = phase_curve(false);
    let patches = [patch(300.0, 450.0)];
    let old = [2000.0; 4];
    with_cx(|cx| {
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            source: &zero, element_materials: None };
        assert!(matches!(engine.advance_with_radiation_and_entropy(
            cx, problem, None, &old, 0.1, config(),
            EnthalpyRadiation::Ambient { patches: &patches, config: radiation_config() }, 1e-10,
        ), Err(EnthalpyEntropyError::Physical(EnthalpyError::InvalidInput(_)))));
        assert!(engine.evaluate_endpoint_with_radiation_and_entropy(
            cx, problem, None, &old, &[3000.0; 4], 0.1, config(),
            EnthalpyRadiation::PointwiseAmbient(&radiation),
        ).is_err(), "500 K is outside the original 450 K surface support");
        assert!(engine.thermal_entropy_change_j_k(cx, &old[..3], &old).is_err());
        with_cancelled_cx(|cancelled| {
            assert!(matches!(engine.advance_with_radiation_and_entropy(
                cancelled, problem, None, &old, 0.1, config(),
                EnthalpyRadiation::PointwiseAmbient(&radiation), 1e-10,
            ), Err(EnthalpyEntropyError::Physical(EnthalpyError::Conduction(
                ConductionError::Cancelled { .. }
            )))));
            assert!(engine.thermal_entropy_change_j_k(cancelled, &old, &old).is_err());
        });
    });
    assert_eq!(old, [2000.0; 4]);
}
