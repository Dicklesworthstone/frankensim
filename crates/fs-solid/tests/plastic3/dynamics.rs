//! Independent rigid-body, one-mode and support-work oracles for native inertia.
use super::*;

const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}
fn problem<'a>(
    law: &'a J2Plasticity,
    prescribed: &'a [(usize, f64)],
    forces: &'a [(usize, f64)],
) -> SmallStrainTetProblem<'a, J2Plasticity> {
    SmallStrainTetProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        laws: SmallStrainLawField::Uniform(law),
        prescribed_m: prescribed,
        nodal_forces_n: forces,
        load_path: &[1.0],
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    }
}

#[test]
fn g1_inertia_free_translation_and_force_impulse_preserve_native_momentum() {
    with_cx(|cx| {
        let law = j2();
        let velocity = [[0.2, -0.1, 0.05]; 4];
        let dt = 0.002;
        for total_force in [0.0, 2.0] {
            let forces: Vec<_> = (0..4).map(|n| (3 * n, total_force / 4.0)).collect();
            let step = problem(&law, &[], &forces)
                .advance_backward_euler(
                    &[0.0; 12],
                    &velocity,
                    &[law.initial_state()],
                    &[6.0],
                    dt,
                    cx,
                )
                .unwrap();
            close(step.energy.reference_mass_kg, 1.0, 1e-14);
            for (v, u) in step.velocity_mps.iter().zip(&step.solution.displacement_m) {
                close(v[0], velocity[0][0] + dt * total_force, 1e-12);
                close(v[1], velocity[0][1], 1e-12);
                close(v[2], velocity[0][2], 1e-12);
                for c in 0..3 {
                    close(u[c], dt * v[c], 1e-15);
                }
            }
            assert!(
                step.solution.element_strain[0]
                    .iter()
                    .all(|s| s.abs() < 1e-14)
            );
            assert_eq!(step.solution.element_states, vec![law.initial_state()]);
            close(
                step.energy.kinetic_change_j + step.energy.numerical_loss_j,
                total_force * step.solution.displacement_m[0][0],
                1e-12,
            );
        }
    });
}

#[test]
fn g1_inertia_translated_continuation_never_discards_velocity_at_the_force_gate() {
    with_cx(|cx| {
        let law = j2();
        let native = problem(&law, &[], &[]);
        let mut old_u = [0.0; 12];
        for node in 0..4 {
            old_u[3 * node] = 1.0;
        }
        let velocity = [[1.0, 0.0, 0.0]; 4];
        let step = native
            .advance_backward_euler(&old_u, &velocity, &[law.initial_state()], &[6.0], 0.001, cx)
            .unwrap();
        for v in &step.velocity_mps {
            close(v[0], 1.0, 1e-10);
        }
        close(step.energy.kinetic_change_j, 0.0, 1e-10);

        // The requested 1e-11 m increment is close to absolute-coordinate
        // resolution at 1 m. Correct convergence or explicit nonconvergence
        // is admissible; accepting u_next=u_old and zero velocity is not.
        let tiny_velocity = [[1e-5, 0.0, 0.0]; 4];
        match native.advance_backward_euler(
            &old_u,
            &tiny_velocity,
            &[law.initial_state()],
            &[6.0],
            1e-6,
            cx,
        ) {
            Ok(step) => {
                for v in &step.velocity_mps {
                    close(v[0], 1e-5, 1e-12);
                }
            }
            Err(fs_solid::plastic3::SmallStrainTetError::NewtonStalled { .. }) => {}
            other => panic!("unexpected translated-motion result: {other:?}"),
        }
        assert_eq!(old_u[0], 1.0);
    });
}

#[test]
fn g1_inertia_stiff_rigid_body_cannot_stop_without_force() {
    with_cx(|cx| {
        // This exact small, stiff specimen passed the old stiffness*travel
        // gate at u_next=0: ||M*v_old/dt||/(maxdiag(K)*dt*v_old) < 1e-10.
        let nodes = NODES.map(|node| node.map(|x| x * 0.01));
        let law =
            J2Plasticity::new(IsotropicElastic::new(1.2e9, 0.2, 0.05).unwrap(), 3e6, 1e8).unwrap();
        let mut native = problem(&law, &[], &[]);
        native.nodes_m = &nodes;
        let velocity = [[0.1, 0.0, 0.0]; 4];
        match native.advance_backward_euler(
            &[0.0; 12],
            &velocity,
            &[law.initial_state()],
            &[1000.0],
            1.0,
            cx,
        ) {
            Ok(step) => {
                for (v, u) in step.velocity_mps.iter().zip(&step.solution.displacement_m) {
                    for c in 0..3 {
                        close(v[c], velocity[0][c], 1e-8);
                        close(u[c], velocity[0][c], 1e-8);
                    }
                }
            }
            // The real sparse solve may refuse this severe stiffness/mass
            // conditioning; it must not replace unresolved motion with rest.
            Err(
                fs_solid::plastic3::SmallStrainTetError::NewtonStalled { .. }
                | fs_solid::plastic3::SmallStrainTetError::LinearSolveFailed { .. },
            ) => {}
            other => panic!("unexpected stiff rigid-body result: {other:?}"),
        }
    });
}

#[test]
fn g1_inertia_unrelated_support_force_does_not_hide_free_momentum() {
    with_cx(|cx| {
        let law = j2();
        let fixed: Vec<_> = (0..12).filter(|d| d % 3 != 0).map(|d| (d, 0.0)).collect();
        // A y-support cancels this load without doing work. All x motion is
        // free and carries no applied force; its momentum must be retained.
        let forces = [(1, 1e20)];
        let velocity = [[0.1, 0.0, 0.0]; 4];
        let step = problem(&law, &fixed, &forces)
            .advance_backward_euler(
                &[0.0; 12],
                &velocity,
                &[law.initial_state()],
                &[6.0],
                0.01,
                cx,
            )
            .unwrap();
        for (v, u) in step.velocity_mps.iter().zip(&step.solution.displacement_m) {
            close(v[0], 0.1, 1e-10);
            close(u[0], 0.001, 1e-12);
        }
    });
}

#[test]
fn g1_inertia_uses_consistent_mass_and_shared_elastic_tangent() {
    with_cx(|cx| {
        let law = j2();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).map(|d| (d, 0.0)).collect();
        let mut velocity = [[0.0; 3]; 4];
        velocity[1][0] = 0.01;
        let dt = 0.01;
        let k = E * (1.0 - NU) / ((1.0 + NU) * (1.0 - 2.0 * NU)) / 6.0;
        for rho in [6.0, 12.0] {
            let m = rho / 60.0;
            let expected_u = m * dt * velocity[1][0] / (m + dt * dt * k);
            let step = problem(&law, &fixed, &[])
                .advance_backward_euler(
                    &[0.0; 12],
                    &velocity,
                    &[law.initial_state()],
                    &[rho],
                    dt,
                    cx,
                )
                .unwrap();
            close(step.solution.displacement_m[1][0], expected_u, 1e-14);
            close(step.velocity_mps[1][0], expected_u / dt, 1e-12);
            close(
                step.energy.kinetic_before_j,
                0.5 * m * velocity[1][0].powi(2),
                1e-15,
            );
            close(
                step.energy.inertial_work_j,
                -k * expected_u * expected_u,
                1e-14,
            );
            assert!(step.energy.numerical_loss_j > 0.0);
            assert_eq!(step.solution.element_states, vec![law.initial_state()]);
        }
    });
}

#[test]
fn g1_inertia_reactions_include_full_prescribed_support_work() {
    with_cx(|cx| {
        let law = j2();
        let displacement = 0.001;
        let fixed: Vec<_> = (0..12)
            .map(|d| (d, if d % 3 == 0 { displacement } else { 0.0 }))
            .collect();
        let step = problem(&law, &fixed, &[])
            .advance_backward_euler(
                &[0.0; 12],
                &[[0.0; 3]; 4],
                &[law.initial_state()],
                &[6.0],
                0.1,
                cx,
            )
            .unwrap();
        let work: f64 = fixed
            .iter()
            .zip(&step.solution.increments[0].reactions_n)
            .map(|((_, u), reaction)| u * reaction)
            .sum();
        close(work, 1e-4, 1e-14);
        close(step.energy.kinetic_after_j, 5e-5, 1e-14);
        close(step.energy.numerical_loss_j, 5e-5, 1e-14);
        close(work, step.energy.inertial_work_j, 1e-14);
    });
}

#[test]
fn g0_inertia_refuses_missing_inventory_invalid_time_and_load_paths() {
    with_cx(|cx| {
        let law = j2();
        let mut native = problem(&law, &[], &[]);
        let state = [law.initial_state()];
        for density in [&[][..], &[0.0][..], &[f64::NAN][..]] {
            assert!(
                native
                    .advance_backward_euler(&[0.0; 12], &[[0.0; 3]; 4], &state, density, 0.1, cx)
                    .is_err()
            );
        }
        for dt in [0.0, -1.0, f64::INFINITY, 1e-200, 1e200] {
            assert!(
                native
                    .advance_backward_euler(&[0.0; 12], &[[0.0; 3]; 4], &state, &[6.0], dt, cx)
                    .is_err()
            );
        }
        assert!(
            native
                .advance_backward_euler(&[0.0; 12], &[[0.0; 3]; 3], &state, &[6.0], 0.1, cx)
                .is_err()
        );
        native.load_path = &[0.0, 1.0];
        assert!(
            native
                .advance_backward_euler(&[0.0; 12], &[[0.0; 3]; 4], &state, &[6.0], 0.1, cx)
                .is_err()
        );
        assert_eq!(state, [law.initial_state()]);
    });
}
