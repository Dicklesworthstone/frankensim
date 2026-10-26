//! Physical endpoint checks reuse the primal residual, without derivative gates.

use super::*;
use fs_conduction::material::LiquidMassFractionConductivity;
use fs_conduction::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
};

fn right_tet() -> ConductionMesh {
    ConductionMesh::new(
        TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    )
    .unwrap()
}

#[test]
fn latent_nodal_defects_remain_visible_when_temperature_and_total_energy_agree() {
    let mesh = right_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let curve = phase_curve(true);
    let material = ConductivityModel::isotropic_declared(6.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    let old = [2000.0; 4];
    let candidate = [2000.0, 2001.0, 1999.0, 2000.0];
    with_cx(|cx| {
        let engine = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[ReferenceEnthalpyMaterial {
                curve: &curve,
                reference_density_kg_m3: 24.0,
            }],
            &[0],
            budget(),
        )
        .unwrap();
        let evaluated = engine
            .evaluate_endpoint(cx, problem, None, &old, &candidate, 1.0, config(), None)
            .unwrap();
        assert_eq!(evaluated.temperature, vec![350.0; 4]);
        assert_eq!(evaluated.free_dofs, 4);
        assert_eq!(evaluated.stored_energy_change_j, 0.0);
        assert_eq!(evaluated.energy_residual_j, 0.0);
        // V=1/6 and rho0=24 give one kg at each vertex. T is constant,
        // so the complete residual is precisely [0, 1, -1, 0] joules.
        assert!((evaluated.residual_norm_j - 2.0_f64.sqrt()).abs() < 1e-10);
        assert!(evaluated.residual_norm_j > config().energy_tolerance_j);
        // An enthalpy/temperature corner is a physical state even though
        // a derivative at the corner can be refused by an adjoint owner.
        let corner = [1000.0; 4];
        let evaluated = engine
            .evaluate_endpoint(cx, problem, None, &corner, &corner, 1.0, config(), None)
            .unwrap();
        assert!(evaluated.residual_norm_j < 1e-10);
    });
    assert_eq!(old, [2000.0; 4]);
    assert_eq!(candidate, [2000.0, 2001.0, 1999.0, 2000.0]);
}

#[test]
fn endpoint_uses_phase_scaled_transport_and_refuses_invalid_inputs_and_cancellation() {
    let mesh = right_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let curve = phase_curve(false);
    let material = ConductivityModel::isotropic_declared(6.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    let old = [1000.0; 4];
    let candidate = [1000.0, 1010.0, 990.0, 1000.0];
    with_cx(|cx| {
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 24.0, budget())
            .unwrap()
            .with_phase_conductivity(
                cx,
                vec![LiquidMassFractionConductivity::declared(2.0, 2.0).unwrap()],
            )
            .unwrap();
        let evaluated = engine
            .evaluate_endpoint(cx, problem, None, &old, &candidate, 1.0, config(), None)
            .unwrap();
        // h'=10 J/kg/K, m_i=1 kg. The right tet's k=12 spatial
        // residual at T=[300,301,299,300] is [0,2,-2,0] W.
        assert!((evaluated.residual_norm_j - 288.0_f64.sqrt()).abs() < 1e-10);
        assert_eq!(evaluated.stored_energy_change_j, 0.0);
        assert_eq!(evaluated.energy_residual_j, 0.0);
        assert!(
            engine
                .evaluate_endpoint(
                    cx,
                    problem,
                    None,
                    &old,
                    &candidate[..3],
                    1.0,
                    config(),
                    None,
                )
                .is_err()
        );
        let mut invalid_config = config();
        invalid_config.newton.relative_tolerance = f64::NAN;
        assert!(matches!(
            engine.evaluate_endpoint(
                cx,
                problem,
                None,
                &old,
                &candidate,
                1.0,
                invalid_config,
                None,
            ),
            Err(EnthalpyError::Newton(_))
        ));
        assert!(
            engine
                .evaluate_endpoint(cx, problem, None, &old, &[5000.0; 4], 1.0, config(), None,)
                .is_err()
        );
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                engine.evaluate_endpoint(
                    cancelled,
                    problem,
                    None,
                    &old,
                    &candidate,
                    1.0,
                    config(),
                    None,
                ),
                Err(EnthalpyError::Conduction(ConductionError::Cancelled { .. }))
            ));
        });
    });
}

#[test]
fn arbitrary_endpoint_separates_signed_radiation_from_convection_without_a_solve() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("exterior", |_| true, ThermalBc::robin(HTC, 300.0).unwrap())
        .unwrap()
        .finish()
        .unwrap();
    let curve = phase_curve(true);
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    with_cx(|cx| {
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        for ambient in [270.0_f64, 400.0] {
            let patches = [patch(ambient, 500.0)];
            let evaluated = engine
                .evaluate_endpoint(
                    cx,
                    problem,
                    None,
                    &[2000.0; 4],
                    &[2000.0; 4],
                    0.25,
                    config(),
                    Some(&patches),
                )
                .unwrap();
            let convection = 3.0_f64.sqrt() * HTC * 50.0;
            let radiation = 3.0_f64.sqrt()
                * EPSILON
                * STEFAN_BOLTZMANN_W_M2_K4
                * (350.0_f64.powi(4) - ambient.powi(4));
            assert!((evaluated.convective_out_w - convection).abs() < 1e-10);
            assert!((evaluated.radiation_out_w - radiation).abs() < 1e-9);
            let energy = 0.25 * (convection + radiation);
            assert!((evaluated.energy_residual_j - energy).abs() < 1e-9);
            // The regular tet has identical nodal boundary loads; four
            // equal entries give ||R||_2 = |sum R|/2.
            assert!((evaluated.residual_norm_j - energy.abs() / 2.0).abs() < 1e-9);
            assert_eq!(evaluated.stored_energy_change_j, 0.0);
            assert_eq!(evaluated.convective_robin_fluxes.len(), 1);
        }
    });
}

#[test]
fn radiation_only_endpoint_accounts_for_signed_heat_without_a_convection_port() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("exterior", |_| true, ThermalBc::adiabatic()).unwrap().finish().unwrap();
    let original = boundary.clone();
    let curve = phase_curve(true);
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        source: &source, element_materials: None };
    with_cx(|cx| {
        let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh,
            &[ReferenceEnthalpyMaterial { curve: &curve, reference_density_kg_m3: RHO }],
            &[0], budget()).unwrap();
        for ambient in [270.0_f64, 350.0, 400.0] {
            let patches = [patch(ambient, 500.0)];
            let got = engine.evaluate_endpoint(cx, problem, None, &[2000.0; 4], &[2000.0; 4],
                0.25, config(), Some(&patches)).unwrap();
            let heat = 3.0_f64.sqrt() * EPSILON * STEFAN_BOLTZMANN_W_M2_K4
                * (350.0_f64.powi(4) - ambient.powi(4));
            assert!((got.radiation_out_w - heat).abs() < 1e-9);
            assert!((got.energy_residual_j - 0.25 * heat).abs() < 1e-9);
            assert!((got.residual_norm_j - 0.125 * heat.abs()).abs() < 1e-9);
            assert_eq!(got.convective_out_w, 0.0);
            assert!(got.convective_robin_fluxes.is_empty());
            assert_eq!(got.stored_energy_change_j, 0.0);
        }
    });
    assert_eq!(boundary, original);
    assert!(!boundary.has_robin());
}

#[test]
fn radiation_only_primal_can_be_rechecked_and_perturbed_without_advancing_history() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("exterior", |_| true, ThermalBc::adiabatic()).unwrap().finish().unwrap();
    let curve = phase_curve(false);
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let source = ScalarField::Uniform(0.0);
    let patches = [patch(400.0, 500.0)];
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        source: &source, element_materials: None };
    let old = [1000.0; 4];
    with_cx(|cx| {
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let solved = engine.advance_with_ambient_radiation(cx, problem, None, &old, 0.1,
            config(), &patches, radiation_config()).unwrap();
        let got = engine.evaluate_endpoint(cx, problem, None, &old,
            &solved.conduction.specific_enthalpy_j_kg, 0.1, config(), Some(&patches)).unwrap();
        assert!(got.residual_norm_j <= solved.physical_residual_tolerance_j);
        assert!(got.energy_residual_j.abs() <= config().energy_tolerance_j);
        assert!(got.radiation_out_w < 0.0 && got.stored_energy_change_j > 0.0);
        assert!(got.convective_robin_fluxes.is_empty());
        let perturbed: Vec<_> = solved.conduction.specific_enthalpy_j_kg.iter().map(|h| h + 1.0).collect();
        let bad = engine.evaluate_endpoint(cx, problem, None, &old, &perturbed, 0.1,
            config(), Some(&patches)).unwrap();
        assert!(bad.residual_norm_j > 1e-3);
        assert!(bad.energy_residual_j.abs() > 1e-3);
    });
    assert_eq!(old, [1000.0; 4]);
}

#[test]
fn radiation_only_and_separate_air_faces_keep_independent_heat_accounts() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("exterior", |f| f.centroid[2] == 0.0, ThermalBc::adiabatic()).unwrap()
        .remainder("air", ThermalBc::robin(HTC, 300.0).unwrap()).unwrap().finish().unwrap();
    let curve = phase_curve(true);
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        source: &source, element_materials: None };
    with_cx(|cx| {
        let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let patches = [patch(400.0, 500.0)];
        let got = engine.evaluate_endpoint(cx, problem, None, &[2000.0; 4], &[2000.0; 4],
            0.25, config(), Some(&patches)).unwrap();
        let convection = 3.0 * 3.0_f64.sqrt() / 4.0 * HTC * 50.0;
        let radiation = 3.0_f64.sqrt() / 4.0 * EPSILON * STEFAN_BOLTZMANN_W_M2_K4
            * (350.0_f64.powi(4) - 400.0_f64.powi(4));
        assert!((got.convective_out_w - convection).abs() < 1e-10);
        assert!((got.radiation_out_w - radiation).abs() < 1e-9);
        assert!((got.energy_residual_j - 0.25 * (convection + radiation)).abs() < 1e-9);
        assert_eq!(got.convective_robin_fluxes.len(), 1);
        assert_eq!(got.convective_robin_fluxes[0].region, "air");
    });
}
