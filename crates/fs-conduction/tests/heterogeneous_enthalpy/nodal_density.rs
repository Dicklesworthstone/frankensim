//! G1/G3: graded reference storage in the original heterogeneous heat owner.
use super::*;
use fs_conduction::ThermalBc;

const DENSITY: [f64; 8] = [3.0, 15.0, 3.0, 3.0, 6.0, 18.0, 6.0, 6.0];
const MASSES: [f64; 8] = [0.225, 0.325, 0.225, 0.225, 0.35, 0.45, 0.35, 0.35];

fn materials<'c>(
    a: &'c EquilibriumEnthalpyPhaseCurve,
    b: &'c EquilibriumEnthalpyPhaseCurve,
) -> [ReferenceEnthalpyMaterial<'c>; 2] {
    [
        ReferenceEnthalpyMaterial {
            curve: a,
            reference_density_kg_m3: 2.0,
        },
        ReferenceEnthalpyMaterial {
            curve: b,
            reference_density_kg_m3: 6.0,
        },
    ]
}

#[test]
fn graded_multi_chart_storage_drives_sensible_and_latent_heating() {
    let mesh = mesh();
    let a = chart(0.0, 700.0);
    let b = chart(100.0, 700.0);
    let materials = materials(&a, &b);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    let old = [40.0, 40.0, 40.0, 40.0, 150.0, 150.0, 150.0, 150.0];
    let target = [50.0, 60.0, 70.0, 80.0, 180.0, 200.0, 220.0, 240.0];
    let temperature = [325.0, 330.0, 335.0, 340.0, 450.0, 450.0, 450.0, 450.0];
    let dt = 0.2;
    let mut q = [0.0; 8];
    // Independent right-tet stiffness k*V*(grad N_a dot grad N_b), and
    // inverse of the P1 source matrix V*(I+11^T)/20; no production assembly.
    for element in 0..2 {
        let start = 4 * element;
        let t = &temperature[start..start + 4];
        let conduction = [
            (3.0 * t[0] - t[1] - t[2] - t[3]) / 3.0,
            (t[1] - t[0]) / 3.0,
            (t[2] - t[0]) / 3.0,
            (t[3] - t[0]) / 3.0,
        ];
        let load: [f64; 4] = std::array::from_fn(|i| {
            let vertex = start + i;
            MASSES[vertex] * (target[vertex] - old[vertex]) / dt + conduction[i]
        });
        let sum: f64 = load.iter().sum();
        for i in 0..4 {
            q[start + i] = 120.0 * (load[i] - sum / 5.0);
        }
    }
    let source = ScalarField::nodal("graded two-chart heating", 8, q.to_vec()).unwrap();
    with_cx(|cx| {
        let stepper =
            HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh, &materials, &[0, 1], budget())
                .unwrap()
                .with_nodal_reference_density(cx, &DENSITY)
                .unwrap();
        let mass_before = stepper.reference_nodal_masses_kg().to_vec();
        for (actual, expected) in mass_before.iter().zip(MASSES) {
            assert!((actual - expected).abs() < 2e-15);
        }
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &conductivity,
            element_materials: None,
            source: &source,
        };
        let result = stepper
            .advance(cx, problem, None, &old, dt, config())
            .unwrap();
        for i in 0..8 {
            assert!((result.specific_enthalpy_j_kg[i] - target[i]).abs() < 2e-8);
            assert!((result.temperature[i] - temperature[i]).abs() < 2e-8);
        }
        let energy: f64 = (0..8).map(|i| MASSES[i] * (target[i] - old[i])).sum();
        assert!((result.stored_energy_change_j - energy).abs() < 1e-9);
        assert!(result.energy_residual_j.abs() < 1e-10);
        for i in 4..8 {
            assert!((result.liquid_mass_fraction[i] - (target[i] - 100.0) / 200.0).abs() < 1e-10);
        }
        assert_eq!(mass_before, stepper.reference_nodal_masses_kg());
        assert_eq!(
            stepper.nodal_reference_density_kg_m3(),
            Some(DENSITY.as_slice())
        );
        assert!(std::ptr::eq(stepper.phase_curve_at(0).unwrap(), &a));
        assert!(std::ptr::eq(stepper.phase_curve_at(7).unwrap(), &b));
        assert_eq!(stepper.element_material_ids(), &[0, 1]);
        assert_eq!(stepper.vertex_material_ids(), &[0, 0, 0, 0, 1, 1, 1, 1]);
        // A scalar source owner must refuse this declared graded inventory.
        assert!(stepper.reference_material(0).is_none());
        assert!(stepper.reference_material(1).is_none());
    });
}

#[test]
fn nodal_density_jvp_and_pullback_match_latent_storage_and_perturbed_solves() {
    let mesh = mesh();
    let a = chart(0.0, 700.0);
    let b = chart(100.0, 700.0);
    let materials = materials(&a, &b);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::uniform("latent deposition", 12.0).unwrap();
    let old = [160.0, 180.0, 200.0, 220.0, 170.0, 190.0, 210.0, 230.0];
    let seed = [0.7, -0.4, 0.2, 1.1, 0.3, -0.6, 0.9, -0.2];
    let direction = [0.4, -0.7, 0.2, 0.5, -0.3, 0.8, -0.6, 0.1];
    let dt = 0.5;
    with_cx(|cx| {
        let stepper =
            HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh, &materials, &[0, 1], budget())
                .unwrap()
                .with_nodal_reference_density(cx, &DENSITY)
                .unwrap();
        let problem = || ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &conductivity,
            element_materials: None,
            source: &source,
        };
        let response = stepper
            .linearize_step(cx, problem(), None, &old, dt, config())
            .unwrap();
        let gradient = response.pullback(cx, &seed, linear()).unwrap();
        let nodal = gradient.nodal_reference_density.as_ref().unwrap();
        let dr = response
            .apply_nodal_reference_density_jacobian(cx, &direction)
            .unwrap();
        close(
            dot(nodal, &direction),
            -dot(&gradient.transport.adjoint, &dr),
        );
        assert_eq!(
            response
                .transport()
                .temperature_pullback(cx, &seed)
                .unwrap(),
            vec![0.0; 8]
        );
        for i in 0..8 {
            let expected = old[i] + 0.25 / MASSES[i];
            assert!((response.primal().specific_enthalpy_j_kg[i] - expected).abs() < 1e-9);
            close(gradient.transport.previous_specific_enthalpy[i], seed[i]);
        }
        // Independent exact latent solution and P1 signed mass derivative.
        for element in 0..2 {
            let start = 4 * element;
            let sum: f64 = direction[start..start + 4].iter().sum();
            for i in start..start + 4 {
                let mass_direction = (sum + direction[i]) / 120.0;
                close(dr[i], mass_direction * 0.25 / MASSES[i]);
            }
            close(
                nodal[start..start + 4].iter().sum(),
                gradient.material_reference_density[element],
            );
        }
        let expected = five_point(|epsilon| {
            let density: Vec<_> = DENSITY
                .iter()
                .zip(direction)
                .map(|(rho, d)| rho + epsilon * d)
                .collect();
            let perturbed =
                HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh, &materials, &[0, 1], budget())
                    .unwrap()
                    .with_nodal_reference_density(cx, &density)
                    .unwrap();
            let result = perturbed
                .advance(cx, problem(), None, &old, dt, config())
                .unwrap();
            dot(&result.specific_enthalpy_j_kg, &seed)
        });
        close(dot(nodal, &direction), expected);
        let element = response
            .apply_reference_density_jacobian(cx, &[0.4, -0.7])
            .unwrap();
        let nodal_constant = response
            .apply_nodal_reference_density_jacobian(
                cx,
                &[0.4, 0.4, 0.4, 0.4, -0.7, -0.7, -0.7, -0.7],
            )
            .unwrap();
        for i in 0..8 {
            close(element[i], nodal_constant[i]);
        }
        assert!(
            response
                .apply_nodal_reference_density_jacobian(cx, &[1.0; 7])
                .is_err()
        );
        let mut invalid = direction;
        invalid[6] = f64::NAN;
        assert!(
            response
                .apply_nodal_reference_density_jacobian(cx, &invalid)
                .is_err()
        );
        with_cancelled_cx(|cancelled| {
            assert!(
                response
                    .apply_nodal_reference_density_jacobian(cancelled, &direction)
                    .is_err()
            );
            assert!(response.pullback(cancelled, &seed, linear()).is_err());
        });
    });
}

#[test]
fn prescribed_bath_uses_graded_mass_and_masks_nodal_density_residual_controls() {
    let mesh = mesh();
    let a = chart(0.0, 700.0);
    let b = chart(100.0, 700.0);
    let materials = materials(&a, &b);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("bath", |_| true, ThermalBc::dirichlet(425.0).unwrap())
        .unwrap()
        .finish()
        .unwrap();
    let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::uniform("generation", 4.0).unwrap();
    let old = [100.0; 8];
    let dt = 0.5;
    with_cx(|cx| {
        let stepper =
            HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh, &materials, &[0, 1], budget())
                .unwrap()
                .with_nodal_reference_density(cx, &DENSITY)
                .unwrap();
        let problem = || ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &conductivity,
            element_materials: None,
            source: &source,
        };
        let response = stepper
            .linearize_step(cx, problem(), None, &old, dt, config())
            .unwrap();
        let result = response.primal();
        let energy = 1.0 * (600.0 - 100.0) + 1.5 * (50.0 - 100.0);
        assert!((result.stored_energy_change_j - energy).abs() < 1e-10);
        assert!((result.dirichlet_in_w - (energy / dt - 4.0 / 3.0)).abs() < 1e-10);
        assert_eq!(
            response
                .apply_nodal_reference_density_jacobian(cx, &[1.0; 8])
                .unwrap(),
            vec![0.0; 8]
        );
        let gradient = response.pullback(cx, &[1.0; 8], linear()).unwrap();
        assert_eq!(gradient.nodal_reference_density, Some(vec![0.0; 8]));
        let next = stepper
            .advance(
                cx,
                problem(),
                None,
                &result.specific_enthalpy_j_kg,
                dt,
                config(),
            )
            .unwrap();
        assert_eq!(next.stored_energy_change_j, 0.0);
        assert!((next.dirichlet_in_w + 4.0 / 3.0).abs() < 1e-10);
    });
}

#[test]
fn constant_profiles_preserve_bits_but_do_not_claim_scalar_source_identity() {
    let mesh = mesh();
    let a = chart(0.0, 700.0);
    let b = chart(100.0, 700.0);
    let materials = materials(&a, &b);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::uniform("generation", 12.0).unwrap();
    let old = [50.0, 50.0, 50.0, 50.0, 180.0, 180.0, 180.0, 180.0];
    with_cx(|cx| {
        let prepare = || {
            HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh, &materials, &[0, 1], budget())
                .unwrap()
        };
        let original = prepare();
        assert!(std::ptr::eq(
            original.reference_material(0).unwrap().curve,
            &a
        ));
        assert_eq!(
            original
                .reference_material(1)
                .unwrap()
                .reference_density_kg_m3,
            6.0
        );
        assert!(original.nodal_reference_density_kg_m3().is_none());
        let mut profile = [2.0, 2.0, 2.0, 2.0, 6.0, 6.0, 6.0, 6.0];
        let nodal = prepare()
            .with_nodal_reference_density(cx, &profile)
            .unwrap();
        assert_eq!(
            original.reference_nodal_masses_kg(),
            nodal.reference_nodal_masses_kg()
        );
        assert!(nodal.reference_material(0).is_none());
        profile[0] = 9.0;
        assert_eq!(nodal.nodal_reference_density_kg_m3().unwrap()[0], 2.0);
        let problem = || ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &conductivity,
            element_materials: None,
            source: &source,
        };
        let a = original
            .linearize_step(cx, problem(), None, &old, 0.5, config())
            .unwrap();
        let b = nodal
            .linearize_step(cx, problem(), None, &old, 0.5, config())
            .unwrap();
        assert_eq!(
            a.primal().specific_enthalpy_j_kg,
            b.primal().specific_enthalpy_j_kg
        );
        assert_eq!(a.primal().temperature, b.primal().temperature);
        assert!(
            a.pullback(cx, &[1.0; 8], linear())
                .unwrap()
                .nodal_reference_density
                .is_none()
        );
        assert!(
            a.apply_nodal_reference_density_jacobian(cx, &[1.0; 8])
                .is_err()
        );
        for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut density = DENSITY;
            density[3] = value;
            assert!(
                prepare()
                    .with_nodal_reference_density(cx, &density)
                    .is_err()
            );
        }
        assert!(
            prepare()
                .with_nodal_reference_density(cx, &DENSITY[..7])
                .is_err()
        );
        let candidate = prepare();
        with_cancelled_cx(|cancelled| {
            assert!(
                candidate
                    .with_nodal_reference_density(cancelled, &DENSITY)
                    .is_err()
            );
        });
        assert_eq!(
            original
                .reference_material(1)
                .unwrap()
                .reference_density_kg_m3,
            6.0
        );
    });
}
