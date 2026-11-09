//! G1/G3 signed nodal heat against independent P1 modes and boundary work.

use super::*;
use fs_conduction::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
};
use fs_math::dd::Dd;

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 2e-9, "{actual} != {expected}");
}

#[test]
fn signed_interval_heat_is_additive_and_preserves_the_original_p1_operator() {
    let mesh = tetrahedron();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(12.0);
    let curve = sensible_curve();
    let old = [100.0; 4];
    let heat = [0.0, 0.7, -0.7, 0.0].map(Dd::from_f64);
    with_cx(|cx| {
        let native = EnthalpyBackwardEuler::uniform(
            cx, &mesh, &curve, REFERENCE_DENSITY, budget(),
        ).unwrap();
        let p = problem(&mesh, &boundary, &material, &source);
        let step = native.advance_with_nodal_heat(
            cx, p, None, &old, 0.25, config(), &heat,
        ).unwrap();
        // On this right tetrahedron K*[0,1,-1,0] = k*V*[0,1,-1,0].
        // Hence (rho*V/4 + dt*k*V/Cp)*dh = Q. The independent
        // uniform volumetric term adds dt*source/rho to every enthalpy.
        let mass = REFERENCE_DENSITY * VOLUME / 4.0;
        let amplitude = 0.7 / (mass + 0.25 * 2.0 * VOLUME / 5.0);
        let mean = 100.0 + 0.25 * 12.0 / REFERENCE_DENSITY;
        for (actual, expected) in step.specific_enthalpy_j_kg.iter().zip([
            mean, mean + amplitude, mean - amplitude, mean,
        ]) {
            close(*actual, expected);
        }
        close(step.endpoint.energy.nodal_heat_j.to_f64(), 0.0);
        close(step.endpoint.energy.volumetric_source_j.to_f64(), 0.5);
        close(step.endpoint.energy.stored_energy_change_j.to_f64(), 0.5);
        assert!(step.endpoint.maximum_free_residual_j <= config().energy_tolerance_j);
        assert!(step.endpoint.energy.residual_j.to_f64().abs() <= config().energy_tolerance_j);

        // Equal total source does not prove a coupled endpoint: changing only
        // the antisymmetric distribution leaves global energy unchanged.
        let changed = [0.0, 0.6, -0.6, 0.0].map(Dd::from_f64);
        let audit = native.evaluate_endpoint_with_nodal_heat(
            cx, p, None, &old, &step.specific_enthalpy_j_kg, 0.25, config(), &changed,
        ).unwrap();
        assert!(audit.energy.residual_j.to_f64().abs() <= config().energy_tolerance_j);
        assert!(audit.maximum_free_residual_j > 0.09);
        assert!(audit.residual_norm_j > 0.13);
        assert_eq!(old, [100.0; 4]);
    });
}

#[test]
fn prescribed_vertices_return_their_actual_supplied_heat_to_the_reservoir() {
    let mesh = tetrahedron();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("bath", |_| true, ThermalBc::dirichlet(300.0).unwrap())
        .unwrap().finish().unwrap();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(4.0);
    let curve = sensible_curve();
    let heat = [1.0, -2.0, 3.0, 4.0].map(Dd::from_f64);
    with_cx(|cx| {
        let native = EnthalpyBackwardEuler::uniform(
            cx, &mesh, &curve, REFERENCE_DENSITY, budget(),
        ).unwrap();
        let old = [150.0; 4];
        let step = native.advance_with_nodal_heat(
            cx, problem(&mesh, &boundary, &material, &source),
            None, &old, 0.5, config(), &heat,
        ).unwrap();
        let storage = REFERENCE_DENSITY * VOLUME * (100.0 - 150.0);
        close(step.endpoint.energy.stored_energy_change_j.to_f64(), storage);
        close(step.endpoint.energy.dirichlet_in_j.to_f64(), storage - 6.0 - 4.0 * VOLUME * 0.5);
        close(step.endpoint.energy.nodal_heat_j.to_f64(), 6.0);
        assert_eq!(step.endpoint.free_dofs, 0);
        assert_eq!(step.endpoint.residual_j, vec![0.0; 4]);
        assert!(step.endpoint.energy.residual_j.to_f64().abs() <= config().energy_tolerance_j);
        assert!(step.endpoint.temperature.iter().all(|t| (*t - 300.0).abs() < 1e-12));
        assert_eq!(old, [150.0; 4]);
    });
}

#[test]
fn heterogeneous_heat_keeps_each_original_latent_chart_and_reference_mass() {
    let mesh = ConductionMesh::new(
        TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
        vec![
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
            [2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0], [2.0, 0.0, 1.0],
        ],
    ).unwrap();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let solid = sensible_curve();
    let pcm = latent_curve();
    let materials = [
        ReferenceEnthalpyMaterial { curve: &solid, reference_density_kg_m3: 2.0 },
        ReferenceEnthalpyMaterial { curve: &pcm, reference_density_kg_m3: 6.0 },
    ];
    let old = [100.0, 100.0, 100.0, 100.0, 180.0, 180.0, 180.0, 180.0];
    let heat = [0.25, 0.25, 0.25, 0.25, -0.5, 0.5, -0.25, 0.25].map(Dd::from_f64);
    with_cx(|cx| {
        let native = HeterogeneousEnthalpyBackwardEuler::new(
            cx, &mesh, &materials, &[0, 1],
            EnthalpyBudget { max_vertices: 8, max_elements: 2 },
        ).unwrap();
        let p = problem(&mesh, &boundary, &material, &source);
        let step = native.advance_with_nodal_heat(
            cx, p, None, &old, 0.2, config(), &heat,
        ).unwrap();
        for vertex in 0..8 {
            let mass = if vertex < 4 { 1.0 / 12.0 } else { 0.25 };
            close(step.specific_enthalpy_j_kg[vertex], old[vertex] + heat[vertex].to_f64() / mass);
            if vertex >= 4 {
                close(step.endpoint.temperature[vertex], 350.0);
                close(step.liquid_mass_fraction[vertex],
                    (step.specific_enthalpy_j_kg[vertex] - 100.0) / 200.0);
            }
        }
        close(step.endpoint.energy.stored_energy_change_j.to_f64(), 1.0);
        let audit = native.evaluate_endpoint_with_nodal_heat(
            cx, p, None, &old, &step.specific_enthalpy_j_kg, 0.2, config(), &heat,
        ).unwrap();
        assert!(audit.maximum_free_residual_j <= config().energy_tolerance_j);
        assert!(audit.energy.residual_j.to_f64().abs() <= config().energy_tolerance_j);
        assert!(std::ptr::eq(native.phase_curve_at(7).unwrap(), &pcm));
    });
}

#[test]
fn invalid_heat_domain_and_cancellation_preserve_the_original_history() {
    let mesh = tetrahedron();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let curve = sensible_curve();
    let old = [100.0; 4];
    with_cx(|cx| {
        let native = EnthalpyBackwardEuler::uniform(
            cx, &mesh, &curve, REFERENCE_DENSITY, budget(),
        ).unwrap();
        let p = problem(&mesh, &boundary, &material, &source);
        for heat in [
            vec![Dd::ZERO; 3],
            vec![Dd::from_f64(f64::NAN); 4],
            vec![Dd::from_f64(-1e6); 4],
        ] {
            assert!(native.advance_with_nodal_heat(
                cx, p, None, &old, 0.1, config(), &heat,
            ).is_err());
        }
        assert!(native.evaluate_endpoint_with_nodal_heat(
            cx, p, None, &old, &[100.0; 3], 0.1, config(), &[Dd::ZERO; 4],
        ).is_err());
        with_cancelled_cx(|cancelled| {
            assert!(native.advance_with_nodal_heat(
                cancelled, p, None, &old, 0.1, config(), &[Dd::ZERO; 4],
            ).is_err());
        });
        let retry = native.advance_with_nodal_heat(
            cx, p, None, &old, 0.1, config(), &[Dd::ZERO; 4],
        ).unwrap();
        assert_eq!(retry.specific_enthalpy_j_kg, old);
        assert_eq!(old, [100.0; 4]);
    });
}
