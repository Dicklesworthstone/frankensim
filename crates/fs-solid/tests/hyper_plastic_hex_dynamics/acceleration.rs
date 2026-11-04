//! G1/G3: source density weights gravity and general reference acceleration
//! through the same mass that advances native Hex8 momentum.

use super::*;

fn graded_motion(
    formulation: FiniteHexFormulation,
    density_scale: f64,
    cx: &Cx<'_>,
) -> FiniteThermoplasticHexMotion<'static> {
    // On the unit cube rho = density_scale * (2 + 4 x). Total mass is 4*scale,
    // and its first x moment is (7/3)*scale; uniform density loses that moment.
    let density = [NODES.map(|p| density_scale * (2.0 + 4.0 * p[0]))];
    FiniteThermoplasticHexMotion::try_new_with_reference_density(
        initial(family(), formulation, cx),
        &density,
        &[[0.0; 3]; 8],
        cx,
    )
    .unwrap()
}

#[test]
fn graded_body_acceleration_retains_exact_nodal_force_and_moment() {
    with_cx(|cx| {
        let motion = graded_motion(FiniteHexFormulation::Standard, 1.0, cx);
        let gravity = [0.0, 0.0, -9.81];
        let loads = motion
            .reference_acceleration_load_n(&[gravity; 8], cx)
            .unwrap();
        assert_eq!(loads.len(), 8);
        let mut force = [0.0; 3];
        let mut torque = [0.0; 3];
        for &(dof, value) in &loads {
            let node = dof / 3;
            assert_eq!(dof % 3, 2);
            // Integral rho*N_a: (5/12) on x=0 and (7/12) on x=1.
            let nodal_mass = if NODES[node][0] == 0.0 {
                5.0 / 12.0
            } else {
                7.0 / 12.0
            };
            near(value, nodal_mass * gravity[2], 2e-13, 2e-14);
            force[2] += value;
            torque[0] += NODES[node][1] * value;
            torque[1] -= NODES[node][0] * value;
        }
        near(force[2], 4.0 * gravity[2], 3e-13, 2e-14);
        near(torque[0], 2.0 * gravity[2], 3e-13, 2e-14);
        near(torque[1], -(7.0 / 3.0) * gravity[2], 3e-13, 2e-14);

        // Spatial acceleration a_x=x must integrate rho*N_a*x, not nodal
        // lumped mass times a(node). The exact low/high x weights are 1/6,5/12.
        let acceleration = NODES.map(|p| [p[0], 0.0, 0.0]);
        let loads = motion
            .reference_acceleration_load_n(&acceleration, cx)
            .unwrap();
        assert_eq!(loads.len(), 8);
        for (dof, value) in loads {
            assert_eq!(dof % 3, 0);
            let expected = if NODES[dof / 3][0] == 0.0 {
                1.0 / 6.0
            } else {
                5.0 / 12.0
            };
            near(value, expected, 3e-14, 2e-14);
        }
    });
}

#[test]
fn graded_gravity_advances_uniformly_at_every_density_and_keeps_material_cold() {
    with_cx(|cx| {
        let acceleration = [0.2, -0.3, -0.5];
        let dt = 0.01;
        for formulation in [
            FiniteHexFormulation::Standard,
            FiniteHexFormulation::MeanDilatation,
        ] {
            for scale in [0.5, 2.0] {
                let motion = graded_motion(formulation, scale, cx);
                let loads = motion
                    .reference_acceleration_load_n(&[acceleration; 8], cx)
                    .unwrap();
                let step = motion
                    .advance(&[300.0], &[], &loads, dt, settings(), controls(), cx)
                    .unwrap();
                for (u, v) in step
                    .motion
                    .solid()
                    .displacement_m()
                    .iter()
                    .zip(step.motion.velocity_mps())
                {
                    for axis in 0..3 {
                        near(u[axis], acceleration[axis] * dt * dt, 2e-11, 2e-9);
                        near(v[axis], acceleration[axis] * dt, 2e-9, 2e-9);
                    }
                }
                near(step.motion.reference_mass_kg(), 4.0 * scale, 2e-13, 2e-14);
                near(step.motion.solid().free_energy_j(), 0.0, 1e-10, 0.0);
                assert_eq!(step.material_energy.plastic_dissipation_j, 0.0);
                assert_eq!(step.material_energy.creep_dissipation_j, 0.0);
                near(step.material_energy.caloric_heat_j, 0.0, 1e-10, 0.0);
                let kinetic = 0.5 * 4.0 * scale * dt * dt * norm2(acceleration);
                near(step.inertia.kinetic_after_j, kinetic, 2e-10, 2e-9);
                near(step.inertia.numerical_loss_j, kinetic, 2e-10, 2e-9);
                near(
                    step.material_energy.dead_load_endpoint_work_j,
                    2.0 * kinetic,
                    2e-10,
                    2e-9,
                );
                near(step.balance_error_j, 0.0, 2e-9, 0.0);
                // Dead reference loading is unchanged after the translated
                // state; it does not sample the material density again.
                assert_eq!(
                    step.motion
                        .reference_acceleration_load_n(&[acceleration; 8], cx)
                        .unwrap(),
                    loads
                );
                assert_eq!(motion.velocity_mps(), [[0.0; 3]; 8]);
            }
        }
    });
}

#[test]
fn reference_acceleration_rejects_bad_fields_and_cancellation_before_publication() {
    with_cx(|cx| {
        let motion = graded_motion(FiniteHexFormulation::Standard, 1.0, cx);
        assert!(
            motion
                .reference_acceleration_load_n(&[[0.0; 3]; 7], cx)
                .is_err()
        );
        for bad in [f64::NAN, f64::INFINITY] {
            let mut acceleration = [[0.0; 3]; 8];
            acceleration[4][1] = bad;
            assert!(
                motion
                    .reference_acceleration_load_n(&acceleration, cx)
                    .is_err()
            );
        }
        assert!(
            motion
                .reference_acceleration_load_n(&[[0.0; 3]; 8], cx)
                .unwrap()
                .is_empty()
        );
        let gate = CancelGate::new();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 7,
                    kernel_id: 178,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(
                motion
                    .reference_acceleration_load_n(&[[1.0; 3]; 8], &cancelled)
                    .is_err()
            );
        });
        assert_eq!(motion.velocity_mps(), [[0.0; 3]; 8]);
        near(motion.reference_mass_kg(), 4.0, 2e-13, 2e-14);
    });
}
