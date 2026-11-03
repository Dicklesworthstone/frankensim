//! Original graded inventory through the actual conservative finite-strain owner.
//! Mass entries use independent barycentric monomial integrals, not quadrature.
use super::*;
use fs_material::hyper::discrete_gradient::bulk::DiscreteBulkSettings;

const CORNERS: [[f64; 4]; 1] = [[3.0, 15.0, 3.0, 3.0]];
type Polynomial = Vec<([usize; 4], f64)>;

fn multiply(a: &Polynomial, b: &Polynomial) -> Polynomial {
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

fn measure(density: [f64; 4], curvature: [f64; 3]) -> Polynomial {
    let rho: Polynomial = density
        .into_iter()
        .enumerate()
        .map(|(i, value)| {
            let mut powers = [0; 4];
            powers[i] = 1;
            (powers, value)
        })
        .collect();
    // X_i = lambda_(i+1) + a_i*lambda_(i+1)^2, so its diagonal
    // Jacobian has determinant product(1+2*a_i*lambda_(i+1)).
    curvature.iter().enumerate().fold(rho, |poly, (i, a)| {
        let mut powers = [0; 4];
        powers[i + 1] = 1;
        multiply(&poly, &vec![([0; 4], 1.0), (powers, 2.0 * a)])
    })
}

fn exact_mass(curvature: [f64; 3]) -> [[f64; 10]; 10] {
    let shapes: Vec<Polynomial> = (0..10)
        .map(|i| {
            let mut powers = [0; 4];
            if i < 4 {
                powers[i] = 1;
                let mut quadratic = powers;
                quadratic[i] = 2;
                vec![(quadratic, 2.0), (powers, -1.0)]
            } else {
                let [a, b] = TET10_EDGES[i - 4];
                powers[a] = 1;
                powers[b] = 1;
                vec![(powers, 4.0)]
            }
        })
        .collect();
    let weight = measure(CORNERS[0], curvature);
    core::array::from_fn(|a| {
        core::array::from_fn(|b| integral(&multiply(&multiply(&shapes[a], &shapes[b]), &weight)))
    })
}

fn full_load(loads: &[(usize, f64)]) -> [[f64; 3]; 10] {
    let mut result = [[0.0; 3]; 10];
    for &(dof, force) in loads {
        result[dof / 3][dof % 3] += force;
    }
    result
}

fn exact_motion(
    mass: &[[f64; 10]; 10],
    element: &[usize; 10],
    nodes: &[[f64; 3]],
    u: &[[f64; 3]],
    v: &[[f64; 3]],
) -> (f64, [f64; 3], [f64; 3]) {
    let mut momentum = [[0.0; 3]; 10];
    for a in 0..10 {
        for b in 0..10 {
            for c in 0..3 {
                momentum[element[a]][c] += mass[a][b] * v[element[b]][c];
            }
        }
    }
    let (mut kinetic, mut linear, mut angular) = (0.0, [0.0; 3], [0.0; 3]);
    for i in 0..10 {
        let position = core::array::from_fn(|c| nodes[i][c] + u[i][c]);
        let spin = cross(position, momentum[i]);
        for c in 0..3 {
            kinetic += 0.5 * v[i][c] * momentum[i][c];
            linear[c] += momentum[i][c];
            angular[c] += spin[c];
        }
    }
    (kinetic, linear, angular)
}

#[test]
fn g1_graded_p2_mass_is_exact_on_affine_and_curved_reference_geometry() {
    with_cx(false, |cx| {
        let (affine, elements) = mesh(cx);
        let material = law();
        for curvature in [[0.0; 3], [0.4, 0.3, 0.2]] {
            let nodes: Vec<[f64; 3]> = affine
                .iter()
                .map(|x| core::array::from_fn(|i| x[i] + curvature[i] * x[i] * x[i]))
                .collect();
            let mut definition = definition(&nodes, &elements, &material, &[DENSITY]);
            if curvature != [0.0; 3] {
                for points in [4, 5] {
                    definition.quadrature_points_per_dir = points;
                    definition.maximum_quadrature_points = points.pow(3);
                    assert!(
                        definition
                            .prepare_with_reference_density(
                                &CORNERS,
                                &[[0.0; 3]; 10],
                                &[[0.0; 3]; 10],
                                &[],
                                cx,
                            )
                            .is_err()
                    );
                }
                definition.quadrature_points_per_dir = 6;
                definition.maximum_quadrature_points = 216;
            }
            let solid = definition
                .prepare_with_reference_density(&CORNERS, &[[0.0; 3]; 10], &[[0.0; 3]; 10], &[], cx)
                .unwrap();
            let expected = exact_mass(curvature);
            // Every original consistent-mass column, including signed corner
            // entries and all three physical component blocks, reaches M*a.
            for b in 0..10 {
                for c in 0..3 {
                    let mut acceleration = [[0.0; 3]; 10];
                    acceleration[elements[0][b]][c] = 1.0;
                    let actual = full_load(
                        &solid
                            .reference_acceleration_load_n(&acceleration, cx)
                            .unwrap(),
                    );
                    for a in 0..10 {
                        for d in 0..3 {
                            near(
                                actual[elements[0][a]][d],
                                if c == d { expected[a][b] } else { 0.0 },
                                3e-14,
                            );
                        }
                    }
                }
            }
            let mass = integral(&measure(CORNERS[0], curvature));
            let volume = integral(&measure([1.0; 4], curvature));
            near(solid.reference_mass_kg(), mass, 2e-14);
            near(solid.reference_density_kg_m3()[0], mass / volume, 3e-14);
            assert_eq!(
                solid.reference_corner_density_kg_m3(),
                Some(CORNERS.as_slice())
            );
            assert!(
                expected[0][4] < 0.0,
                "signed consistent entries must survive"
            );
        }
    });
}

#[test]
fn g1_graded_free_body_preserves_energy_and_both_momenta_during_deformation() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let material = law();
        let initial_u: Vec<[f64; 3]> = nodes
            .iter()
            .map(|x| {
                [
                    0.08 * x[0] + 0.48 * x[1] * x[2],
                    -0.06 * x[1] + 0.03 * x[2],
                    0.04 * x[2],
                ]
            })
            .collect();
        let initial_v: Vec<[f64; 3]> = nodes
            .iter()
            .zip(&initial_u)
            .map(|(x, u)| {
                let spin = cross(
                    [0.3, -0.2, 0.5],
                    core::array::from_fn(|i| x[i] + u[i] - 0.25),
                );
                [
                    0.2 + spin[0] + 0.1 * (x[0] - 0.25),
                    -0.15 + spin[1],
                    0.1 + spin[2],
                ]
            })
            .collect();
        let mass = exact_mass([0.0; 3]);
        let (kinetic0, linear0, angular0) =
            exact_motion(&mass, &elements[0], &nodes, &initial_u, &initial_v);
        let (averaged_kinetic, _, _) =
            independent_motion(&nodes, &elements[0], &initial_u, &initial_v);
        assert!(
            (kinetic0 - averaged_kinetic).abs() > 1e-4,
            "a cell mean must lose a real inertia moment"
        );
        let mut solid = definition(&nodes, &elements, &material, &[DENSITY])
            .prepare_with_reference_density(&CORNERS, &initial_u, &initial_v, &[], cx)
            .unwrap();
        near(solid.kinetic_energy_j(), kinetic0, 2e-14);
        let energy0 = independent_storage(&nodes, &elements[0], &initial_u, &material) + kinetic0;
        let original_mass = solid.reference_mass_kg().to_bits();
        for _ in 0..16 {
            let step = solid
                .advance(0.01, &[], &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
                .unwrap();
            check_work(&step);
            let next = &step.solid;
            let (kinetic, linear, angular) = exact_motion(
                &mass,
                &elements[0],
                &nodes,
                next.displacement_m(),
                next.velocity_mps(),
            );
            near(next.kinetic_energy_j(), kinetic, 2e-11);
            near(
                independent_storage(&nodes, &elements[0], next.displacement_m(), &material)
                    + kinetic,
                energy0,
                WORK_ALLOWANCE,
            );
            for c in 0..3 {
                near(linear[c], linear0[c], 1e-10);
                near(angular[c], angular0[c], 1e-10);
            }
            assert_eq!(next.reference_mass_kg().to_bits(), original_mass);
            assert_eq!(
                next.reference_corner_density_kg_m3(),
                Some(CORNERS.as_slice())
            );
            solid = step.solid;
        }
        let old_edge: [f64; 3] =
            core::array::from_fn(|c| nodes[1][c] + initial_u[1][c] - nodes[0][c] - initial_u[0][c]);
        let edge: [f64; 3] = core::array::from_fn(|c| {
            nodes[1][c] + solid.displacement_m()[1][c] - nodes[0][c] - solid.displacement_m()[0][c]
        });
        let squared = |v: [f64; 3]| v.into_iter().map(|x| x * x).sum::<f64>();
        assert!(
            (squared(edge) - squared(old_edge)).abs() > 1e-3,
            "the body must actually deform"
        );
    });
}

fn formulations() -> [ConservativeFormulation; 3] {
    [
        ConservativeFormulation::FullMaterial,
        ConservativeFormulation::ExactIncompressible {
            maximum_relative_volume_residual: 1e-9,
            maximum_initial_relative_volume_rate_per_s: 1e-10,
            maximum_point_volume_work_defect: 1e-12,
        },
        ConservativeFormulation::OriginalMaterial {
            bulk: DiscreteBulkSettings {
                minimum_volume_ratio: 0.25,
                maximum_volume_ratio: 2.5,
                maximum_iterations: 96,
                absolute_pressure_tolerance_pa: 1e-11,
                relative_pressure_tolerance: 1e-12,
                absolute_work_tolerance_pa: 1e-10,
                relative_work_tolerance: 1e-12,
            },
            maximum_relative_volume_residual: 1e-9,
            maximum_point_volume_work_defect: 1e-12,
        },
    ]
}

#[test]
fn g1_graded_acceleration_has_exact_motion_resultant_torque_and_support_reactions() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let material = law();
        let acceleration = [1.2, -9.81, 0.4];
        let velocity = [0.2, -0.1, 0.05];
        // Integral of rho*x divided by the unit total mass for rho=3+12*x.
        let centroid = [0.325, 0.225, 0.225];
        let expected_torque = cross(centroid, acceleration);
        for formulation in formulations() {
            let mut definition = definition(&nodes, &elements, &material, &[DENSITY]);
            definition.formulation = formulation;
            let mut solid = definition
                .prepare_with_reference_density(&CORNERS, &[[0.0; 3]; 10], &[velocity; 10], &[], cx)
                .unwrap();
            let loads = solid
                .reference_acceleration_load_n(&[acceleration; 10], cx)
                .unwrap();
            let full = full_load(&loads);
            let mut total = [0.0; 3];
            let mut torque = [0.0; 3];
            for (x, force) in nodes.iter().zip(full) {
                let moment = cross(*x, force);
                for c in 0..3 {
                    total[c] += force[c];
                    torque[c] += moment[c];
                }
            }
            for c in 0..3 {
                near(total[c], acceleration[c], 3e-14);
                near(torque[c], expected_torque[c], 3e-14);
            }
            // Negative corner shape integrals are physical consistent loads,
            // not negative density; deleting them destroys the exact moments.
            assert!(full[elements[0][0]][0] < 0.0);
            for interval in 1..=3 {
                let step = solid
                    .advance(0.02, &[], &loads, [0.0; 3], settings(), WORK_ALLOWANCE, cx)
                    .unwrap();
                check_work(&step);
                let time = f64::from(interval) * 0.02;
                for (u, v) in step
                    .solid
                    .displacement_m()
                    .iter()
                    .zip(step.solid.velocity_mps())
                {
                    for c in 0..3 {
                        near(
                            u[c],
                            time * velocity[c] + 0.5 * time * time * acceleration[c],
                            1e-10,
                        );
                        near(v[c], velocity[c] + time * acceleration[c], 1e-9);
                    }
                }
                near(step.solid.reference_mass_kg(), 1.0, 1e-14);
                near(step.solid.stored_energy_j(), 0.0, 1e-10);
                solid = step.solid;
            }
        }
        let solid = definition(&nodes, &elements, &material, &[DENSITY])
            .prepare_with_reference_density(&CORNERS, &[[0.0; 3]; 10], &[[0.0; 3]; 10], &[], cx)
            .unwrap();
        let loads = solid
            .reference_acceleration_load_n(&[acceleration; 10], cx)
            .unwrap();
        let supports: Vec<_> = (0..30).map(|dof| (dof, 0.0)).collect();
        let step = solid
            .advance(
                0.02,
                &supports,
                &loads,
                [0.0; 3],
                settings(),
                WORK_ALLOWANCE,
                cx,
            )
            .unwrap();
        let actual = full_load(&step.reactions_n);
        let expected = full_load(&loads);
        for i in 0..10 {
            for c in 0..3 {
                near(actual[i][c], -expected[i][c], 2e-12);
            }
        }
        check_work(&step);
        assert_eq!(step.energy.support_work_j, 0.0);
    });
}

#[test]
fn g0_scalar_bits_and_refused_graded_candidates_preserve_reference_inventory() {
    with_cx(false, |cx| {
        let (nodes, elements) = mesh(cx);
        let material = law();
        let definition = definition(&nodes, &elements, &material, &[DENSITY]);
        let velocity = [[0.2, -0.1, 0.05]; 10];
        let original = definition.prepare(&[[0.0; 3]; 10], &velocity, cx).unwrap();
        let constant = [[DENSITY; 4]];
        let graded = definition
            .prepare_with_reference_density(&constant, &[[0.0; 3]; 10], &velocity, &[], cx)
            .unwrap();
        assert!(original.reference_corner_density_kg_m3().is_none());
        assert_eq!(
            original.reference_density_kg_m3(),
            graded.reference_density_kg_m3()
        );
        assert_eq!(
            original.reference_mass_kg().to_bits(),
            graded.reference_mass_kg().to_bits()
        );
        assert_eq!(
            original.kinetic_energy_j().to_bits(),
            graded.kinetic_energy_j().to_bits()
        );
        let acceleration: Vec<_> = nodes
            .iter()
            .map(|x| [1.0 + x[0], -0.5 + 2.0 * x[1], 0.25 - x[2]])
            .collect();
        let loads = original
            .reference_acceleration_load_n(&acceleration, cx)
            .unwrap();
        assert_eq!(
            loads,
            graded
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap()
        );
        let old = original
            .advance(0.01, &[], &loads, [0.0; 3], settings(), WORK_ALLOWANCE, cx)
            .unwrap();
        let new = graded
            .advance(0.01, &[], &loads, [0.0; 3], settings(), WORK_ALLOWANCE, cx)
            .unwrap();
        assert_eq!(old.solid.displacement_m(), new.solid.displacement_m());
        assert_eq!(old.solid.velocity_mps(), new.solid.velocity_mps());
        assert_eq!(old.inertia, new.inertia);
        assert_eq!(old.energy, new.energy);
        for density in [
            &[][..],
            &[[1.0, 0.0, 2.0, 3.0]],
            &[[1.0, f64::NAN, 2.0, 3.0]],
            &[[1.0, f64::INFINITY, 2.0, 3.0]],
        ] {
            assert!(
                definition
                    .prepare_with_reference_density(density, &[[0.0; 3]; 10], &velocity, &[], cx,)
                    .is_err()
            );
        }
        let mut limited = definition;
        limited.maximum_mass_assembly_bytes = 0;
        assert!(
            limited
                .prepare_with_reference_density(&CORNERS, &[[0.0; 3]; 10], &velocity, &[], cx)
                .is_err()
        );
        assert!(new.solid.reference_acceleration_load_n(&[], cx).is_err());
        assert!(
            new.solid
                .reference_acceleration_load_n(&[[f64::NAN; 3]; 10], cx)
                .is_err()
        );
        with_cx(true, |cancelled| {
            assert!(
                definition
                    .prepare_with_reference_density(
                        &CORNERS,
                        &[[0.0; 3]; 10],
                        &velocity,
                        &[],
                        cancelled,
                    )
                    .is_err()
            );
            assert!(
                new.solid
                    .reference_acceleration_load_n(&acceleration, cancelled)
                    .is_err()
            );
            assert!(
                new.solid
                    .advance(
                        0.01,
                        &[],
                        &loads,
                        [0.0; 3],
                        settings(),
                        WORK_ALLOWANCE,
                        cancelled
                    )
                    .is_err()
            );
        });
        assert_eq!(
            new.solid.reference_corner_density_kg_m3(),
            Some(constant.as_slice())
        );
        assert_eq!(new.solid.displacement_m(), old.solid.displacement_m());
    });
}
