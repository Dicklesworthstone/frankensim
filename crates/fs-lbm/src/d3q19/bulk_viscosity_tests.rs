use super::{
    CollisionModel3, E3, Q3, W3, collide_cell3, equilibrium3, rheology_stress_components3,
};

fn source(velocity: [f64; 3], force: [f64; 3]) -> [f64; Q3] {
    let uf: f64 = velocity.iter().zip(force).map(|(u, f)| u * f).sum();
    core::array::from_fn(|q| {
        let e = [E3[q].0, E3[q].1, E3[q].2].map(f64::from);
        let eu: f64 = e.iter().zip(velocity).map(|(e, u)| e * u).sum();
        let ef: f64 = e.iter().zip(force).map(|(e, f)| e * f).sum();
        W3[q] * (3.0 * ef + 9.0 * eu * ef - 3.0 * uf)
    })
}

// Independent second-order Hermite quadrature: prescribed physical stress,
// with the held Guo half-source removed to obtain the stored populations.
fn state(rho: f64, velocity: [f64; 3], stress: [[f64; 3]; 3], force: [f64; 3]) -> [f64; Q3] {
    let equilibrium = equilibrium3(rho, velocity);
    let source = source(velocity, force);
    let trace = stress[0][0] + stress[1][1] + stress[2][2];
    core::array::from_fn(|q| {
        let e = [E3[q].0, E3[q].1, E3[q].2].map(f64::from);
        let mut contraction = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                contraction += e[i] * e[j] * stress[i][j];
            }
        }
        equilibrium[q] + 4.5 * W3[q] * (contraction - trace / 3.0) - 0.5 * source[q]
    })
}

#[test]
fn independent_bulk_relaxes_forced_trace_and_deviatoric_stress() {
    let (rho, velocity) = (1.3, [0.04, -0.03, 0.02]);
    let stress = [
        [0.003, -0.002, 0.001],
        [-0.002, -0.001, 0.0005],
        [0.001, 0.0005, 0.002],
    ];
    let trace = stress[0][0] + stress[1][1] + stress[2][2];
    for shear_tau in [0.56, 0.9, 1.7] {
        for bulk_tau in [0.6, 1.0, 1.8] {
            for force in [[0.0; 3], [0.001, -0.002, 0.0007]] {
                let model = CollisionModel3::CentralMomentBulk {
                    shear_rate: 1.0 / shear_tau,
                    bulk_rate: 1.0 / bulk_tau,
                    higher_order_rate: 1.3,
                };
                let initial = state(rho, velocity, stress, force);
                let post = collide_cell3(initial, model, force).unwrap();
                let equilibrium = equilibrium3(rho, velocity);
                let source = source(velocity, force);
                let mut observed = [[0.0; 3]; 3];
                let mut impulse = [0.0; 3];
                for q in 0..Q3 {
                    let e = [E3[q].0, E3[q].1, E3[q].2].map(f64::from);
                    let physical_neq = post[q] - equilibrium[q] - 0.5 * source[q];
                    for i in 0..3 {
                        impulse[i] += e[i] * (post[q] - initial[q]);
                        for j in 0..3 {
                            observed[i][j] += e[i] * e[j] * physical_neq;
                        }
                    }
                }
                assert!((post.iter().sum::<f64>() - rho).abs() < 5e-14);
                for i in 0..3 {
                    assert!((impulse[i] - force[i]).abs() < 5e-14);
                    for j in 0..3 {
                        let isotropic = if i == j { trace / 3.0 } else { 0.0 };
                        let expected = (1.0 - 1.0 / shear_tau) * (stress[i][j] - isotropic)
                            + (1.0 - 1.0 / bulk_tau) * isotropic;
                        assert!((observed[i][j] - expected).abs() < 5e-14);
                    }
                }
            }
        }
    }
}

#[test]
fn equal_bulk_and_shear_rates_preserve_legacy_collision_bits() {
    for shear_tau in [0.503, 0.8, 1.7] {
        for higher_order_rate in [0.6, 1.0, 1.8] {
            let legacy = CollisionModel3::CentralMoment {
                second_order_rate: 1.0 / shear_tau,
                higher_order_rate,
            };
            let independent = CollisionModel3::CentralMomentBulk {
                shear_rate: 1.0 / shear_tau,
                bulk_rate: 1.0 / shear_tau,
                higher_order_rate,
            };
            for force in [[0.0; 3], [0.0004, -0.0002, 0.0003]] {
                let mut populations = equilibrium3(1.3, [0.04, -0.03, 0.02]);
                for (q, value) in populations.iter_mut().enumerate() {
                    *value += (q as f64 - 9.0) * 2e-5;
                }
                assert_eq!(
                    collide_cell3(populations, independent, force)
                        .unwrap()
                        .map(f64::to_bits),
                    collide_cell3(populations, legacy, force)
                        .unwrap()
                        .map(f64::to_bits),
                );
            }
        }
    }
}

#[test]
fn material_strain_reconstruction_separates_bulk_and_shear_under_force() {
    for strain in [
        [[0.003, 0.0, 0.0], [0.0, 0.003, 0.0], [0.0, 0.0, 0.003]],
        [
            [0.003, -0.002, 0.001],
            [-0.002, -0.001, 0.0005],
            [0.001, 0.0005, 0.002],
        ],
    ] {
        let dilation = strain[0][0] + strain[1][1] + strain[2][2];
        let expected = (2.0_f64 * strain.iter().flatten().map(|s| s * s).sum::<f64>()).sqrt();
        for rho in [0.4, 2.5] {
            for shear_tau in [0.57, 1.3] {
                for bulk_tau in [0.6, 1.7] {
                    let stress = core::array::from_fn(|i| {
                        core::array::from_fn(|j| {
                            let isotropic = if i == j { dilation / 3.0 } else { 0.0 };
                            -2.0 * rho / 3.0
                                * (shear_tau * (strain[i][j] - isotropic) + bulk_tau * isotropic)
                        })
                    });
                    let model = CollisionModel3::CentralMomentBulk {
                        shear_rate: 1.0 / shear_tau,
                        bulk_rate: 1.0 / bulk_tau,
                        higher_order_rate: 1.2,
                    };
                    for force in [[0.0; 3], [0.001, -0.002, 0.0007]] {
                        let velocity = [0.04, -0.03, 0.02];
                        let populations = state(rho, velocity, stress, force);
                        let [numerator, fixed] =
                            rheology_stress_components3(&populations, rho, velocity, force, model);
                        assert!((fixed - (2.0_f64 / 3.0).sqrt() * dilation.abs()).abs() < 8e-14);
                        assert!(((numerator / shear_tau).hypot(fixed) - expected).abs() < 8e-14);
                    }
                }
            }
        }
    }
}

#[test]
fn independent_bulk_admits_only_finite_positive_relaxation() {
    for bulk_rate in [0.0, -0.1, 2.0, f64::INFINITY, f64::NAN] {
        let model = CollisionModel3::CentralMomentBulk {
            shear_rate: 1.25,
            bulk_rate,
            higher_order_rate: 1.0,
        };
        assert!(model.validate().is_err());
        assert!(collide_cell3(equilibrium3(1.0, [0.0; 3]), model, [0.0; 3]).is_err());
    }
}
