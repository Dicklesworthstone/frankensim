//! Graded reference inventory through the original prepared load and Newton owner.
use super::*;

#[test]
fn graded_gravity_drives_native_finite_stretch_pressure_and_supported_weight() {
    let nodes = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let tets = [[0, 1, 2, 3]];
    // rho0 = 3 + 12X: mass 1 kg, reference COM (0.325, 0.225, 0.225).
    let density = [3.0, 15.0, 3.0, 3.0];
    let law = material();
    let prescribed: Vec<_> = (0..12).filter(|&d| d != 3).map(|d| (d, 0.0)).collect();
    let p = HyperTetProblem {
        settings: HyperTetSettings {
            load_steps: 3,
            ..HyperTetSettings::default()
        },
        ..problem(&nodes, &tets, &law, &prescribed)
    };
    let stretch = 1.25_f64;
    let gradient = 3.0;
    let piola = 2.0 * (stretch - 1.0 / stretch) + 3.0 * stretch.ln() / stretch;
    with_cx(|cx| {
        let surface =
            ClosedTetPressureSurface::new(&nodes, &tets, Default::default(), 4, cx).unwrap();
        let prepared = p.prepare(cx).unwrap();
        for pressure in [0.0, 0.2] {
            // Exact cubic mass moments: row sum M_1*=0.325, M_11=0.15.
            // P_11/6 = a*0.325 - k*M_11*stretch - pressure/6.
            let a = ((piola + pressure) / 6.0 + gradient * 0.15 * stretch) / 0.325;
            let field = |x: [f64; 3]| {
                Ok(GravitySample3 {
                    potential_j_kg: 0.5 * gradient * x[0] * x[0] - a * x[0],
                    acceleration_m_s2: [a - gradient * x[0], 0.0, 0.0],
                    acceleration_gradient_s2: [[-gradient, 0.0, 0.0], [0.0; 3], [0.0; 3]],
                })
            };
            let solution = if pressure == 0.0 {
                p.solve_with_nodal_density_gravity(&density, &field, cx)
                    .unwrap()
            } else {
                p.solve_with_nodal_density_gravity_and_pressure(
                    &density,
                    &field,
                    &surface,
                    TetPressure::Uniform(pressure),
                    cx,
                )
                .unwrap()
            };
            near(solution.displacement_m[1][0], stretch - 1.0, 2e-10);
            let u: Vec<_> = solution.displacement_m.iter().flatten().copied().collect();
            let gravity = prepared
                .gravity_load_with_nodal_density(&density, &u, &field, cx)
                .unwrap();
            let direct = p
                .gravity_load_with_nodal_density(&density, &u, &field, cx)
                .unwrap();
            assert_eq!(gravity.force_n, direct.force_n);
            assert_eq!(gravity.force_tangent_n_m, direct.force_tangent_n_m);
            // J=1.25 changes current density to rho0/J; original mass remains 1.
            let total_weight: f64 = gravity.force_n.iter().step_by(3).sum();
            near(total_weight, a - gradient * 0.325 * stretch, 2e-10);
            let pressure_force = surface.force_n(&u, pressure, cx).unwrap();
            for i in 0..3 {
                let reaction: f64 = solution
                    .reactions_n
                    .iter()
                    .filter(|(d, _)| d % 3 == i)
                    .map(|(_, value)| value)
                    .sum();
                near(
                    reaction
                        + (0..4)
                            .map(|node| {
                                gravity.force_n[3 * node + i] + pressure_force[3 * node + i]
                            })
                            .sum::<f64>(),
                    0.0,
                    2e-10,
                );
            }
            // Node 1 carries the graded row moments M_1* = 0.325 and M_11 = 0.15. The
            // cell-mean density 6 kg/m3 holds the same 1 kg with moments 0.25 and 0.1,
            // and the graded node-1 force exceeds its force by 0.075 a - 0.05 k stretch
            // (-2.5e-3 N and +5.2e-3 N here). Both rules integrate these fields exactly.
            near(
                gravity.force_n[3],
                0.325 * a - 0.15 * gradient * stretch,
                2e-10,
            );
            let averaged = p.gravity_load(&[6.0], &u, &field, cx).unwrap();
            near(
                averaged.force_n[3],
                0.25 * a - 0.1 * gradient * stretch,
                2e-10,
            );
        }
    });
}

#[test]
fn uniform_nodal_public_gravity_preserves_original_samples_and_refuses_bad_inventory() {
    let mesh = block([1; 3], [1.0; 3]);
    let law = material();
    let p = problem(&mesh.nodes, &mesh.tets, &law, &[]);
    let u = vec![0.0; 3 * mesh.nodes.len()];
    let nodes = vec![7.3; mesh.nodes.len()];
    let cells = vec![7.3; mesh.tets.len()];
    with_cx(|cx| {
        let original = p.gravity_load(&cells, &u, &affine, cx).unwrap();
        let nodal = p
            .gravity_load_with_nodal_density(&nodes, &u, &affine, cx)
            .unwrap();
        assert_eq!(original.force_n, nodal.force_n);
        assert_eq!(original.force_tangent_n_m, nodal.force_tangent_n_m);
        assert_eq!(
            original.potential_energy_j.to_bits(),
            nodal.potential_energy_j.to_bits()
        );
        let mut bad = nodes.clone();
        bad[0] = 0.0;
        for rho in [&bad[..], &nodes[1..]] {
            assert!(
                p.gravity_load_with_nodal_density(rho, &u, &affine, cx)
                    .is_err()
            );
            assert!(
                p.solve_with_nodal_density_gravity(rho, &affine, cx)
                    .is_err()
            );
        }
    });
}
