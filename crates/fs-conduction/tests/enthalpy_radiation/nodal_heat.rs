//! G1/G3 combined signed source, local radiation and prescribed heat reactions.
use super::*;
use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use fs_conduction::PointwiseAmbientRadiation;
use fs_math::dd::Dd;

fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.16e} != {b:.16e}");
}

// Exact barycentric moments of N_i T_h^4, independent of native quadrature.
fn radiative_load(mesh: &ConductionMesh, t: [f64; 4], bath: f64) -> [f64; 4] {
    let mut out = [0.0; 4];
    for face in mesh.boundary() {
        let values = face.vertices.map(|v| t[v as usize]);
        for a in 0..=4 {
            for b in 0..=4 - a {
                let c = 4 - a - b;
                let monomial = values[0].powi(a) * values[1].powi(b) * values[2].powi(c);
                for (i, power) in [a, b, c].into_iter().enumerate() {
                    out[face.vertices[i] as usize] += EPSILON * STEFAN_BOLTZMANN_W_M2_K4
                        * face.area * f64::from(power + 1) * monomial / 105.0;
                }
            }
        }
        for vertex in face.vertices {
            out[vertex as usize] -= EPSILON * STEFAN_BOLTZMANN_W_M2_K4
                * face.area * bath.powi(4) / 3.0;
        }
    }
    out
}

#[test]
fn nonuniform_signed_heat_and_radiation_recover_an_independent_manufactured_endpoint() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .remainder("exterior", ThermalBc::robin(HTC, 300.0).unwrap()).unwrap().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let source = ScalarField::Uniform(7.0);
    let curve = phase_curve(false);
    let patches = [patch(300.0, 600.0)];
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        source: &source, element_materials: None };
    let t = [315.0, 350.0, 385.0, 420.0];
    let old = [420.0, 315.0, 390.0, 380.0].map(|v| CP * (v - 200.0));
    let expected = t.map(|v| CP * (v - 200.0));
    let volume = 2.0_f64.sqrt() / 12.0;
    let area = 3.0_f64.sqrt();
    let dt = 0.05;
    let radiation = radiative_load(&mesh, t, 300.0);
    let sum: f64 = t.iter().sum();
    let heat: [Dd; 4] = std::array::from_fn(|v| Dd::from_f64(
        RHO * volume / 4.0 * (expected[v] - old[v])
            + dt * (K * volume * (2.0 * t[v] - 0.5 * sum)
                + HTC * area / 24.0 * (2.0 * t[v] + sum - 1800.0)
                + radiation[v] - 7.0 * volume / 4.0)));
    assert!(heat.iter().any(|q| q.to_f64() < 0.0));
    assert!(heat.iter().any(|q| q.to_f64() > 0.0));
    with_cx(|cx| {
        let owner = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let step = engine.advance_with_nodal_heat_and_radiation(
            cx, problem, None, &old, dt, config(), &heat,
            EnthalpyRadiation::PointwiseAmbient(&owner),
        ).unwrap();
        for (&h, expected) in step.specific_enthalpy_j_kg.iter().zip(expected) {
            near(h, expected, 3e-8);
        }
        near(step.endpoint.energy.radiation_out_j.to_f64(),
            dt * radiation.iter().sum::<f64>(), 3e-8);
        near(step.endpoint.energy.robin_out_j.to_f64(),
            dt * HTC * area * (sum / 4.0 - 300.0), 3e-8);
        assert!(step.endpoint.pointwise_radiation.is_some());
        assert!(step.endpoint.pointwise_enclosure_radiation.is_none());
        assert!(step.endpoint.maximum_free_residual_j <= config().energy_tolerance_j);
        let mut changed = heat;
        changed[0] = changed[0] + Dd::from_f64(0.01);
        changed[1] = changed[1] - Dd::from_f64(0.01);
        let audit = engine.evaluate_endpoint_with_nodal_heat_and_radiation(
            cx, problem, None, &old, &step.specific_enthalpy_j_kg, dt, config(),
            &changed, EnthalpyRadiation::PointwiseAmbient(&owner),
        ).unwrap();
        assert!(audit.maximum_free_residual_j > 0.009);
        assert!(audit.energy.residual_j.to_f64().abs() <= config().energy_tolerance_j);
    });
}

#[test]
fn prescribed_node_reactions_include_signed_heat_and_actual_radiation_once() {
    let mesh = regular_tet();
    // Three prescribed faces fix all four vertices. The fourth natural trace
    // still radiates, and its power returns through the prescribed reactions.
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("clamp", |f| f.vertices.contains(&0), ThermalBc::dirichlet(350.0).unwrap()).unwrap()
        .remainder("exterior", ThermalBc::adiabatic()).unwrap().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let source = ScalarField::Uniform(7.0);
    let curve = phase_curve(false);
    let patches = [patch(300.0, 600.0)];
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        source: &source, element_materials: None };
    let old = [1000.0; 4];
    let heat = [1.0, -2.0, 3.0, 4.0].map(Dd::from_f64);
    let dt = 0.2;
    with_cx(|cx| {
        let owner = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let step = engine.advance_with_nodal_heat_and_radiation(cx, problem, None, &old,
            dt, config(), &heat, EnthalpyRadiation::PointwiseAmbient(&owner)).unwrap();
        let volume = 2.0_f64.sqrt() / 12.0;
        let stored = RHO * volume * 500.0;
        let emitted = dt * EPSILON * STEFAN_BOLTZMANN_W_M2_K4
            * 3.0_f64.sqrt() / 4.0 * (350.0_f64.powi(4) - 300.0_f64.powi(4));
        near(step.endpoint.energy.stored_energy_change_j.to_f64(), stored, 1e-9);
        near(step.endpoint.energy.radiation_out_j.to_f64(), emitted, 1e-9);
        near(step.endpoint.energy.dirichlet_in_j.to_f64(),
            stored - 6.0 - dt * 7.0 * volume + emitted, 1e-9);
        assert_eq!(step.endpoint.free_dofs, 0);
        assert_eq!(step.endpoint.residual_j, vec![0.0; 4]);
        assert!(step.endpoint.energy.residual_j.to_f64().abs() <= config().energy_tolerance_j);
    });
}

#[test]
fn original_radiation_identity_and_cancelled_steps_preserve_history() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .remainder("exterior", ThermalBc::adiabatic()).unwrap().finish().unwrap();
    let other = ThermalBoundaryBuilder::new(&mesh)
        .remainder("exterior", ThermalBc::adiabatic()).unwrap().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let source = ScalarField::Uniform(0.0);
    let curve = phase_curve(false); let patches = [patch(300.0, 600.0)];
    let old = [1000.0; 4]; let heat = [Dd::ZERO; 4];
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        source: &source, element_materials: None };
    with_cx(|cx| {
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let owner = PointwiseAmbientRadiation::new(cx, &mesh, &other, &patches).unwrap();
        assert!(engine.advance_with_nodal_heat_and_radiation(cx, problem, None,
            &old, 0.1, config(), &heat, EnthalpyRadiation::PointwiseAmbient(&owner)).is_err());
        assert!(engine.advance_with_nodal_heat_and_radiation(cx, problem, None,
            &old, 0.1, config(), &heat,
            EnthalpyRadiation::Ambient { patches: &patches, config: radiation_config() }).is_err());
        let owner = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        with_cancelled_cx(|cancelled| {
            assert!(engine.advance_with_nodal_heat_and_radiation(cancelled, problem, None,
                &old, 0.1, config(), &heat, EnthalpyRadiation::PointwiseAmbient(&owner)).is_err());
        });
        let retry = engine.advance_with_nodal_heat_and_radiation(cx, problem, None,
            &old, 0.1, config(), &heat, EnthalpyRadiation::PointwiseAmbient(&owner)).unwrap();
        for &h in &retry.specific_enthalpy_j_kg { near(h, 1000.0, 1e-10); }
        assert_eq!(old, [1000.0; 4]);
    });
}
