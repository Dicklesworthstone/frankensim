//! G1/G3: original graded inventory in loaded and constrained modal inertia.
use super::*;

const DENSITY: [[f64; 4]; 1] = [[3.0, 15.0, 3.0, 3.0]];

fn mesh(cx: &Cx<'_>) -> (Vec<[f64; 3]>, Vec<[usize; 10]>) {
    promote_p1_to_p2(
        &[[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        &[[0, 1, 2, 3]],
        None,
        cx,
    )
    .unwrap()
}

#[test]
fn g1_graded_loaded_corner_modes_follow_distribution_not_mean_density() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let (mu, bulk, stretch) = (2.0, 20.0, 1.2);
        let material = law(mu, bulk);
        let initial: Vec<_> = nodes
            .iter()
            .map(|x| [(stretch - 1.0) * x[0], 0.0, 0.0])
            .collect();
        let supports: Vec<_> = (0..3 * nodes.len())
            .filter(|&d| d != 3)
            .map(|d| (d, initial[d / 3][d % 3]))
            .collect();
        let mut p = problem(
            &nodes,
            &elements,
            &material,
            &supports,
            &[],
            VolumetricTreatment::OriginalMaterial,
        );
        p.settings.load_steps = 1;
        let pressure = [bulk * (stretch - 1.0); 4];
        let f = [stretch, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        // Integral (dN1/dX1)^2 = 1/10. The original constitutive tangent is
        // density-independent, while integral rho*N1^2 = (sum rho+2*rho1)/2520.
        let stiffness = material.tangent(&f).unwrap()[0][0] / 10.0;
        let mut eigenvalues = Vec::new();
        for density in [DENSITY, [[15.0, 3.0, 3.0, 3.0]]] {
            let spectrum = p
                .solve_vibration_from_with_reference_density(
                    &initial,
                    &pressure,
                    &density,
                    modal_settings(),
                    cx,
                )
                .unwrap();
            let mass = (density[0].iter().sum::<f64>() + 2.0 * density[0][1]) / 2520.0;
            assert_eq!(spectrum.free_dofs(), &[3]);
            assert_eq!(
                spectrum.reference_corner_density_kg_m3(),
                Some(density.as_slice())
            );
            near(spectrum.reference_density_kg_m3()[0], 6.0, 2e-13);
            near(spectrum.reference_mass_kg(), 1.0, 2e-13);
            near(spectrum.mass_kg()[0], mass, 2e-14);
            near(spectrum.stiffness_n_m()[0], stiffness, 2e-9 * stiffness);
            let mode = &spectrum.modes()[0].displacement;
            near(mode.lambda, stiffness / mass, 2e-9 * stiffness / mass);
            near(mode.phi[0] * mode.phi[0] * mass, 1.0, 2e-12);
            eigenvalues.push(mode.lambda);
        }
        near(eigenvalues[1] / eigenvalues[0], 1.8, 2e-10);
    });
}

#[test]
fn g1_graded_exact_pressure_modes_preserve_gauge_and_original_acceleration_loads() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let mu = 2.0;
        let material = law(mu, 20.0);
        let supports: Vec<_> = (0..3 * nodes.len())
            .filter(|d| ![0, 3, 7, 11].contains(d))
            .map(|d| (d, 0.0))
            .collect();
        let mut p = problem(
            &nodes,
            &elements,
            &material,
            &supports,
            &[],
            VolumetricTreatment::ExactIncompressible,
        );
        p.settings.load_steps = 1;
        for pin in [0, 1] {
            p.pressure_gauge = Some((pin, 0.0));
            let spectrum = if pin == 0 {
                p.solve_incompressible_vibration_with_reference_density(
                    &DENSITY,
                    modal_settings(),
                    cx,
                )
            } else {
                p.solve_incompressible_vibration_from_with_reference_density(
                    &vec![[0.0; 3]; nodes.len()],
                    &[0.0; 4],
                    &DENSITY,
                    modal_settings(),
                    cx,
                )
            }
            .unwrap();
            // The exact constraint leaves (N0-N1,-N2,-N3). Independent
            // simplex integrals give K=13*mu/15 and M=23/420 for this rho.
            let (stiffness, mass) = (13.0 * mu / 15.0, 23.0 / 420.0);
            assert_eq!(spectrum.dimension(), 1);
            assert_eq!(
                spectrum.reference_corner_density_kg_m3(),
                Some(DENSITY.as_slice())
            );
            near(spectrum.reference_mass_kg(), 1.0, 2e-13);
            let mode = &spectrum.modes()[0];
            near(mode.coordinates.lambda, stiffness / mass, 2e-8);
            let a = mode.displacement_shape[0];
            for &value in &mode.displacement_shape[1..] {
                near(value / a, -1.0, 2e-11);
            }
            near(a * a * mass, 1.0, 2e-11);
            assert_eq!(mode.pressure_shape[pin], 0.0);
            assert!(
                mode.constraint_relative_residual <= modal_settings().pressure_residual_tolerance
            );
            assert!(
                mode.momentum_relative_residual
                    <= modal_settings().spectrum.relative_residual_tolerance
            );
            let loads = spectrum
                .reference_acceleration_load_n(&vec![[0.0, 0.0, 1.0]; nodes.len()], cx)
                .unwrap();
            let (mut force, mut moment) = ([0.0; 3], [0.0; 3]);
            for &(dof, value) in &loads {
                force[dof % 3] += value;
                let x = nodes[dof / 3];
                if dof % 3 == 2 {
                    moment[0] += x[1] * value;
                    moment[1] -= x[0] * value;
                }
            }
            assert!(loads.iter().any(|&(d, f)| d == 2 && f < 0.0));
            near(force[2], 1.0, 2e-13);
            // rho=3+12*x has center of mass (13/40,9/40,9/40).
            near(moment[0], 9.0 / 40.0, 2e-13);
            near(moment[1], -13.0 / 40.0, 2e-13);
        }
    });
}

type Polynomial = Vec<([usize; 4], f64)>;

fn product(a: &Polynomial, b: &Polynomial) -> Polynomial {
    a.iter()
        .flat_map(|(p, c)| {
            b.iter()
                .map(move |(q, d)| (core::array::from_fn(|i| p[i] + q[i]), c * d))
        })
        .collect()
}

fn integral(poly: &Polynomial) -> f64 {
    let factorial = |n: usize| (1..=n).map(|i| i as f64).product::<f64>();
    poly.iter()
        .map(|(p, c)| {
            c * p.iter().map(|&i| factorial(i)).product::<f64>()
                / factorial(3 + p.iter().sum::<usize>())
        })
        .sum()
}

fn exact_mass(curvature: [f64; 3]) -> [[f64; 10]; 10] {
    let shapes: Vec<Polynomial> = (0..10)
        .map(|i| {
            let mut p = [0; 4];
            if i < 4 {
                p[i] = 1;
                let mut q = p;
                q[i] = 2;
                vec![(q, 2.0), (p, -1.0)]
            } else {
                let [a, b] = TET10_EDGES[i - 4];
                p[a] = 1;
                p[b] = 1;
                vec![(p, 4.0)]
            }
        })
        .collect();
    let density: Polynomial = DENSITY[0]
        .into_iter()
        .enumerate()
        .map(|(i, rho)| {
            let mut p = [0; 4];
            p[i] = 1;
            (p, rho)
        })
        .collect();
    // X_i=lambda_(i+1)+c_i*lambda_(i+1)^2 has diagonal Jacobian.
    // This exact polynomial integral does not use the production quadrature.
    let weight = curvature.iter().enumerate().fold(density, |poly, (i, c)| {
        let mut p = [0; 4];
        p[i + 1] = 1;
        product(&poly, &vec![([0; 4], 1.0), (p, 2.0 * c)])
    });
    core::array::from_fn(|a| {
        core::array::from_fn(|b| integral(&product(&product(&shapes[a], &shapes[b]), &weight)))
    })
}

#[test]
fn g1_graded_modal_mass_and_loading_match_every_affine_and_curved_coefficient() {
    with_cx(|cx| {
        let (affine, elements) = mesh(cx);
        let material = law(2.0, 20.0);
        let supports: Vec<_> = (0..30).filter(|&d| d != 3).map(|d| (d, 0.0)).collect();
        for curvature in [[0.0; 3], [0.4, 0.3, 0.2]] {
            let nodes: Vec<[f64; 3]> = affine
                .iter()
                .map(|x| core::array::from_fn(|i| x[i] + curvature[i] * x[i] * x[i]))
                .collect();
            let mut p = problem(
                &nodes,
                &elements,
                &material,
                &supports,
                &[],
                VolumetricTreatment::OriginalMaterial,
            );
            p.settings.load_steps = 1;
            if curvature != [0.0; 3] {
                p.quadrature_points_per_dir = 6;
            }
            let spectrum = p
                .solve_vibration_with_reference_density(&DENSITY, modal_settings(), cx)
                .unwrap();
            let expected = exact_mass(curvature);
            near(
                spectrum.reference_mass_kg(),
                expected.iter().flatten().sum(),
                3e-13,
            );
            near(spectrum.mass_kg()[0], expected[1][1], 2e-13);
            for b in 0..10 {
                let mut acceleration = [[0.0; 3]; 10];
                acceleration[elements[0][b]][0] = 1.0;
                let mut forces = [0.0; 30];
                for (dof, force) in spectrum
                    .reference_acceleration_load_n(&acceleration, cx)
                    .unwrap()
                {
                    forces[dof] = force;
                }
                for a in 0..10 {
                    near(forces[3 * elements[0][a]], expected[a][b], 2e-13);
                    assert_eq!(forces[3 * elements[0][a] + 1], 0.0);
                    assert_eq!(forces[3 * elements[0][a] + 2], 0.0);
                }
            }
            assert!(
                expected[0][4] < 0.0,
                "signed P2 support coupling must remain present"
            );
        }
    });
}

#[test]
fn g3_graded_modal_inventory_retains_scalar_bits_and_native_refusals() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(2.0, 20.0);
        let supports: Vec<_> = (0..30).filter(|&d| d != 3).map(|d| (d, 0.0)).collect();
        let p = problem(
            &nodes,
            &elements,
            &material,
            &supports,
            &[],
            VolumetricTreatment::OriginalMaterial,
        );
        let scalar = p.solve_vibration(&[6.0], modal_settings(), cx).unwrap();
        let uniform = p
            .solve_vibration_with_reference_density(&[[6.0; 4]], modal_settings(), cx)
            .unwrap();
        assert_eq!(scalar.mass_kg(), uniform.mass_kg());
        assert_eq!(scalar.stiffness_n_m(), uniform.stiffness_n_m());
        assert_eq!(scalar.reference_mass_kg(), uniform.reference_mass_kg());
        assert_eq!(
            scalar.modes()[0].displacement.lambda,
            uniform.modes()[0].displacement.lambda
        );
        let acceleration = [[1.0, -2.0, 3.0]; 10];
        assert_eq!(
            scalar
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap(),
            uniform
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap()
        );
        for invalid in [
            vec![],
            vec![[0.0, 1.0, 2.0, 3.0]],
            vec![[-1.0; 4]],
            vec![[f64::NAN; 4]],
            vec![[f64::INFINITY; 4]],
        ] {
            assert!(
                p.solve_vibration_with_reference_density(&invalid, modal_settings(), cx)
                    .is_err()
            );
        }
        let mut cap = modal_settings();
        cap.maximum_mass_assembly_bytes = 1;
        assert!(
            p.solve_vibration_with_reference_density(&DENSITY, cap, cx)
                .is_err()
        );
        assert!(scalar.reference_acceleration_load_n(&[], cx).is_err());
        assert!(
            scalar
                .reference_acceleration_load_n(&[[f64::NAN; 3]; 10], cx)
                .is_err()
        );
        let mut curved = nodes.clone();
        curved[elements[0][4]][0] += 0.02;
        let curved_problem = problem(
            &curved,
            &elements,
            &material,
            &supports,
            &[],
            VolumetricTreatment::OriginalMaterial,
        );
        assert!(
            curved_problem
                .solve_vibration_with_reference_density(&DENSITY, modal_settings(), cx)
                .is_err()
        );
        let gate = CancelGate::new();
        gate.request();
        with_gate(&gate, |cancelled| {
            assert!(matches!(
                p.solve_vibration_with_reference_density(&DENSITY, modal_settings(), cancelled),
                Err(
                    fs_solid::hyper3::mixed::vibration::MixedVibrationError::Mixed(
                        MixedUpError::Cancelled
                    )
                )
            ));
            assert!(
                scalar
                    .reference_acceleration_load_n(&acceleration, cancelled)
                    .is_err()
            );
        });
    });
}
