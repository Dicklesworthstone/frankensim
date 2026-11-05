//! Full nodal loading, with an independent dense momentum oracle.
use super::*;

fn budget() -> ViscoDynamicPortBudget {
    ViscoDynamicPortBudget {
        maximum_ports: 4,
        maximum_payload_bytes: 1 << 20,
    }
}

#[test]
fn mutual_mean_impedance_matches_full_consistent_mass_and_material_solve() {
    with_cx(|cx, _| {
        let law = coupled();
        let free = [3, 6];
        let solid = system(&law, &free, cx);
        let dt = 0.002;
        let body = solid.dynamic_operator(dt, config(), cx).unwrap();
        let map = [1.0, -0.4, 0.3, 0.8];
        // Full nonsymmetric, rank-deficient, and zero direct matrices. None
        // may be diagonalized or inverted independently of the solid.
        for direct in [
            [90.0, 35.0, -15.0, 55.0],
            [40.0, -20.0, -20.0, 10.0],
            [0.0; 4],
        ] {
            let port = body
                .prepare_port_impedance(2, &map, &direct, budget(), cx)
                .unwrap();
            let mut initial = solid.state(cx).unwrap();
            let mut old = [0.0; 12];
            old[3] = 0.0003;
            old[6] = -0.0001;
            initial.step_internal_force(&old, 0.01, cx).unwrap();
            let mut velocity = [0.0; 12];
            velocity[3] = 0.002;
            velocity[6] = -0.001;
            let mut trajectory = body.start(initial, &velocity, cx).unwrap();
            let force = [
                0.0, 0.0, 0.0, 0.15, 0.0, 0.0, -0.08, 0.0, 0.0, 0.0, 0.0, 0.0,
            ];
            let history = [0.06, -0.035];

            // Assemble A+(h/2) B D B^T directly in physical free coordinates,
            // independent of the cached Schur/Woodbury inverse action.
            let mut held = trajectory.material_history().clone();
            let held = held.step_internal_force(&old, dt, cx).unwrap();
            let mut matrix = [0.0; 4];
            let mut rhs = [0.0; 2];
            for i in 0..2 {
                rhs[i] = (0..2)
                    .map(|j| solid.equilibrium().mass.get(i, j) * velocity[free[j]])
                    .sum::<f64>()
                    + 0.5
                        * dt
                        * (force[free[i]]
                            - held.reduced_mean_force_n[i]
                            - (0..2).map(|a| map[2 * i + a] * history[a]).sum::<f64>());
                for j in 0..2 {
                    let mut loading = 0.0;
                    for a in 0..2 {
                        for b in 0..2 {
                            loading += map[2 * i + a] * direct[2 * a + b] * map[2 * j + b];
                        }
                    }
                    matrix[2 * i + j] = body.effective_mass().get(i, j) + 0.5 * dt * loading;
                }
            }
            let det = matrix[0] * matrix[3] - matrix[1] * matrix[2];
            let mean = [
                (matrix[3] * rhs[0] - matrix[1] * rhs[1]) / det,
                (matrix[0] * rhs[1] - matrix[2] * rhs[0]) / det,
            ];
            let step = trajectory
                .advance_with_port_impedance(&port, &force, &history, cx)
                .unwrap();
            for i in 0..2 {
                close(
                    trajectory.velocity_mps()[free[i]],
                    2.0 * mean[i] - velocity[free[i]],
                    2e-11,
                );
                close(
                    trajectory.material_history().displacement_m()[free[i]],
                    old[free[i]] + dt * mean[i],
                    2e-12,
                );
            }
            let increment: Vec<_> = trajectory
                .material_history()
                .displacement_m()
                .iter()
                .zip(old)
                .map(|(a, b)| a - b)
                .collect();
            let mut nodal_boundary_work = 0.0;
            for i in 0..2 {
                for a in 0..2 {
                    nodal_boundary_work +=
                        map[2 * i + a] * step.port.mean_efforts[a] * increment[free[i]];
                }
            }
            close(step.port.boundary_work_j, nodal_boundary_work, 2e-12);
            close(
                step.port.applied_work_j,
                force.iter().zip(&increment).map(|(f, u)| f * u).sum(),
                2e-12,
            );
            close(
                step.body.external_work_j + step.port.boundary_work_j,
                step.port.applied_work_j,
                2e-11,
            );
            assert!(step.body.constitutive.step_dissipated_j >= 0.0);
            assert!(step.body.balance_residual_j.abs() < 1e-13);
        }
    });
}

#[test]
fn source_history_can_drive_a_resting_solid_through_zero_and_dark_ports() {
    with_cx(|cx, _| {
        let law = material(TensorProny::new_mandel(diagonal([1200.0; 6]), vec![], 0.02).unwrap());
        let solid = system(&law, &[3], cx);
        let dt = 0.03;
        let body = solid.dynamic_operator(dt, config(), cx).unwrap();
        let port = body
            .prepare_port_impedance(2, &[1.0, 0.0], &[0.0; 4], budget(), cx)
            .unwrap();
        let mut trajectory = body
            .start(solid.state(cx).unwrap(), &[0.0; 12], cx)
            .unwrap();
        let step = trajectory
            .advance_with_port_impedance(&port, &[0.0; 12], &[0.1, -90.0], cx)
            .unwrap();
        assert!(step.port.mean_flows[0] < 0.0);
        assert_eq!(step.port.mean_flows[1], 0.0);
        assert_eq!(step.port.mean_efforts, [0.1, -90.0]);
        assert!(step.port.boundary_work_j < 0.0);
        assert_eq!(step.port.applied_work_j, 0.0);
        let energy =
            trajectory.kinetic_energy_j() + trajectory.material_history().energy().stored_j;
        close(energy, -step.port.boundary_work_j, 2e-11);
        assert_eq!(step.body.constitutive.step_dissipated_j, 0.0);
        close(
            trajectory.velocity_mps()[3],
            -dt * 0.1 / body.effective_mass().get(0, 0),
            2e-12,
        );
    });
}

#[test]
fn zero_port_load_preserves_original_motion_and_wrong_owners_refuse_atomically() {
    with_cx(|cx, gate| {
        let law = coupled();
        let solid = system(&law, &[3, 6], cx);
        let body = solid.dynamic_operator(0.001, config(), cx).unwrap();
        let other = solid.dynamic_operator(0.001, config(), cx).unwrap();
        let map = [1.0, -0.4, 0.3, 0.8];
        let zero = body
            .prepare_port_impedance(2, &map, &[0.0; 4], budget(), cx)
            .unwrap();
        let wrong = other
            .prepare_port_impedance(2, &map, &[0.0; 4], budget(), cx)
            .unwrap();
        let mut bare = body
            .start(solid.state(cx).unwrap(), &[0.0; 12], cx)
            .unwrap();
        let mut loaded = bare.clone();
        let mut force = [0.0; 12];
        force[3] = 0.1;
        force[6] = -0.2;
        bare.advance(&force, cx).unwrap();
        let step = loaded
            .advance_with_port_impedance(&zero, &force, &[0.0; 2], cx)
            .unwrap();
        assert_eq!(loaded.velocity_mps(), bare.velocity_mps());
        assert_eq!(
            loaded.material_history().displacement_m(),
            bare.material_history().displacement_m()
        );
        assert_eq!(step.port.boundary_work_j, 0.0);
        let before = loaded.clone();
        assert!(
            loaded
                .advance_with_port_impedance(&wrong, &force, &[0.0; 2], cx)
                .is_err()
        );
        assert!(
            loaded
                .advance_with_port_impedance(&zero, &force, &[f64::NAN, 0.0], cx)
                .is_err()
        );
        assert!(
            body.prepare_port_impedance(usize::MAX, &[], &[], budget(), cx)
                .is_err()
        );
        assert!(
            body.prepare_port_impedance(
                2,
                &map,
                &[0.0; 4],
                ViscoDynamicPortBudget {
                    maximum_payload_bytes: 1,
                    ..budget()
                },
                cx
            )
            .is_err()
        );
        gate.request();
        assert!(
            loaded
                .advance_with_port_impedance(&zero, &force, &[0.0; 2], cx)
                .is_err()
        );
        assert_eq!(loaded.time_s(), before.time_s());
        assert_eq!(loaded.steps(), before.steps());
        assert_eq!(loaded.velocity_mps(), before.velocity_mps());
        assert_eq!(
            loaded.material_history().displacement_m(),
            before.material_history().displacement_m()
        );
    });
}
