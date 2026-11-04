//! G1/G3/G4 graded reference inertia during native heating, flow and impacts.

use super::*;
use fs_solid::hyper3::plastic::thermal::{SolidImpulseControls, TetImpulsePoint};
use fs_solver::CheckedCgConfig;

const RHO: [f64; 4] = [3.0, 9.0, 3.0, 9.0];

fn impulse_controls() -> SolidImpulseControls {
    SolidImpulseControls {
        solve: CheckedCgConfig {
            tolerance: 1.0e-13,
            max_iterations: 64,
            max_corrections: 2,
        },
        maximum_workspace_components: 4096,
        linear_tolerance_kg_m_s: 1.0e-11,
        angular_tolerance_kg_m2_s: 1.0e-11,
        work_tolerance_j: 1.0e-11,
    }
}

fn heated(cx: &Cx<'_>) -> (FiniteThermoplasticTetSolid<'static>, Point) {
    let mut solid = initial(cx);
    let (mut point, _) = Point::at(0.0, 300.0, None);
    for (a, t) in [(0.20, 300.0), (0.20, 450.0)] {
        point = Point::at(a, t, Some(point)).0;
        solid = solid
            .advance(
                &[t],
                &prescribed_x(a),
                &[],
                [0.0; 3],
                0.05,
                settings(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
    }
    assert!(solid.element_states()[0].alpha > 0.0);
    (solid, point)
}

#[test]
fn g1_graded_density_drives_heated_plastic_momentum_and_gravity() {
    with_cx(|cx| {
        let (hot, old) = heated(cx);
        let history = hot.element_states().to_vec();
        let solid = hot.with_nodal_inertia(&RHO, &[[0.0; 3]; 4], cx).unwrap();
        assert_eq!(solid.element_states(), history);
        assert_eq!(solid.reference_nodal_density_kg_m3(), Some(RHO.as_slice()));
        assert_eq!(solid.reference_density_kg_m3(), Some([6.0].as_slice()));
        close(solid.reference_mass_kg().unwrap(), 1.0, 1.0e-14);

        // rho=3+6X+6Z on the unit reference simplex. Exact first moments
        // give nodal weights [9,11,9,11]/40, including restrained-node forces.
        let gravity = solid
            .reference_acceleration_load_n(&[[0.0, 0.0, -10.0]; 4], cx)
            .unwrap();
        assert_eq!(gravity.len(), 4);
        for ((dof, force), (expected_dof, expected_force)) in gravity.iter().zip([
            (2, -9.0 / 4.0),
            (5, -11.0 / 4.0),
            (8, -9.0 / 4.0),
            (11, -11.0 / 4.0),
        ]) {
            assert_eq!(*dof, expected_dof);
            close(*force, expected_force, 1.0e-14);
        }
        // The cubic rho*N_a*N_b integral, not cell-mean or row-sum inertia,
        // maps acceleration 60X in x to these exact full nodal loads.
        let mut acceleration = [[0.0; 3]; 4];
        acceleration[1][0] = 60.0;
        let loads = solid
            .reference_acceleration_load_n(&acceleration, cx)
            .unwrap();
        for ((dof, force), (expected_dof, expected_force)) in
            loads.iter().zip([(0, 3.0), (3, 7.0), (6, 3.0), (9, 3.5)])
        {
            assert_eq!(*dof, expected_dof);
            close(*force, expected_force, 1.0e-14);
        }

        let dt = 0.05;
        let (next, dissipation) = Point::at(0.25, 500.0, Some(old));
        let du = next.a.exp() - old.a.exp();
        let velocity = du / dt;
        let diagonal_mass = 7.0 / 60.0;
        let force = V * next.piola[0] + diagonal_mass * velocity / dt;
        let supports: Vec<_> = (0..12)
            .filter(|dof| *dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let mut refused = controls();
        refused.maximum_numerical_kinetic_loss_j = 0.0;
        assert!(
            solid
                .advance(
                    &[500.0],
                    &supports,
                    &[(3, force)],
                    [0.0; 3],
                    dt,
                    settings(),
                    refused,
                    cx,
                )
                .unwrap_err()
                .to_string()
                .contains("backward-Euler numerical kinetic loss")
        );
        assert_eq!(solid.element_states(), history);
        assert_eq!(solid.temperatures_k(), &[450.0]);
        assert_eq!(solid.velocity_mps(), Some([[0.0; 3]; 4].as_slice()));
        assert_eq!(solid.reference_nodal_density_kg_m3(), Some(RHO.as_slice()));

        let step = solid
            .advance(
                &[500.0],
                &supports,
                &[(3, force)],
                [0.0; 3],
                dt,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        close(
            step.solid.displacement_m()[1][0],
            next.a.exp() - 1.0,
            1.0e-9,
        );
        close(step.solid.velocity_mps().unwrap()[1][0], velocity, 1.0e-8);
        close(step.solid.element_states()[0].alpha, next.alpha, 1.0e-10);
        close(step.energy.plastic_dissipation_j, V * dissipation, 1.0e-6);
        let kinetic = 0.5 * diagonal_mass * velocity * velocity;
        close(step.inertia.unwrap().kinetic_after_j, kinetic, 1.0e-7);
        close(step.inertia.unwrap().numerical_loss_j, kinetic, 1.0e-7);
        close(step.energy.balance_error_j, 0.0, 1.0e-6);
        assert_eq!(
            step.solid.reference_nodal_density_kg_m3(),
            Some(RHO.as_slice())
        );
        close(step.solid.reference_mass_kg().unwrap(), 1.0, 1.0e-14);
        assert_eq!(
            step.solid
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap(),
            loads
        );
        assert_eq!(
            step.solid
                .reference_acceleration_load_n(&[[0.0, 0.0, -10.0]; 4], cx)
                .unwrap(),
            gravity
        );
    });
}

#[test]
fn g1_g3_density_redistribution_changes_thermal_impulse_and_retains_history() {
    with_cx(|cx| {
        let held: Vec<_> = (0..12).filter(|dof| *dof != 3).collect();
        let point = TetImpulsePoint {
            tetrahedron: 0,
            weights: [0.0, 1.0, 0.0, 0.0],
        };
        for (rho, diagonal_mass) in [(RHO, 7.0 / 60.0), ([9.0, 3.0, 9.0, 3.0], 1.0 / 12.0)] {
            let hot = heated(cx)
                .0
                .with_nodal_inertia(&rho, &[[0.0; 3]; 4], cx)
                .unwrap();
            close(hot.reference_mass_kg().unwrap(), 1.0, 1.0e-14);
            let mobility = hot
                .point_directional_inverse_mass_per_kg(
                    &held,
                    point,
                    [1.0, 0.0, 0.0],
                    impulse_controls(),
                    cx,
                )
                .unwrap();
            close(mobility.inverse_mass_per_kg, 1.0 / diagonal_mass, 1.0e-12);
            let impulse = 0.2;
            let event = hot
                .point_impulse_candidate(
                    &held,
                    point,
                    [impulse, 0.0, 0.0],
                    [0.0; 3],
                    impulse_controls(),
                    cx,
                )
                .unwrap();
            close(
                event.solid.velocity_mps().unwrap()[1][0],
                impulse / diagonal_mass,
                1.0e-12,
            );
            close(
                event.impulse.kinetic_change_j,
                0.5 * impulse * impulse / diagonal_mass,
                1.0e-12,
            );
            close(event.impulse.work_residual_j, 0.0, 1.0e-12);
            assert_eq!(event.solid.displacement_m(), hot.displacement_m());
            assert_eq!(event.solid.element_deformation(), hot.element_deformation());
            assert_eq!(event.solid.element_states(), hot.element_states());
            assert_eq!(
                event.solid.element_first_piola_pa(),
                hot.element_first_piola_pa()
            );
            assert_eq!(event.solid.temperatures_k(), hot.temperatures_k());
            assert_eq!(
                event.solid.free_energy_j().to_bits(),
                hot.free_energy_j().to_bits()
            );
            assert_eq!(
                event.solid.internal_energy_j().to_bits(),
                hot.internal_energy_j().to_bits()
            );
            assert_eq!(
                event.solid.reference_nodal_density_kg_m3(),
                Some(rho.as_slice())
            );
            assert_eq!(hot.velocity_mps(), Some([[0.0; 3]; 4].as_slice()));

            // Stop the kicked degree of freedom at its original position.
            // BE loses K numerically; the material receives no physical heat.
            let stopped = event
                .solid
                .advance(
                    &[450.0],
                    &prescribed_x(0.20),
                    &[],
                    [0.0; 3],
                    0.05,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            assert_eq!(
                stopped.solid.reference_nodal_density_kg_m3(),
                Some(rho.as_slice())
            );
            close(stopped.inertia.unwrap().kinetic_after_j, 0.0, 1.0e-20);
            close(
                stopped.inertia.unwrap().numerical_loss_j,
                event.impulse.kinetic_after_j,
                1.0e-12,
            );
            close(stopped.energy.caloric_heat_j, 0.0, 1.0e-6);
        }
    });
}

#[test]
fn g0_g3_g4_uniform_thermal_mass_keeps_bits_and_invalid_bindings_refuse() {
    with_cx(|cx| {
        let velocity = [[0.1, -0.2, 0.3]; 4];
        let scalar = initial(cx).with_inertia(&[6.0], &velocity, cx).unwrap();
        let nodal = initial(cx)
            .with_nodal_inertia(&[6.0; 4], &velocity, cx)
            .unwrap();
        let acceleration = [
            [1.0, 2.0, -3.0],
            [4.0, -5.0, 6.0],
            [7.0, 8.0, 9.0],
            [-1.0; 3],
        ];
        assert_eq!(
            scalar.reference_mass_kg().unwrap().to_bits(),
            nodal.reference_mass_kg().unwrap().to_bits()
        );
        assert_eq!(
            scalar
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap(),
            nodal
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap()
        );
        let held: Vec<_> = (0..12).filter(|dof| *dof != 3).collect();
        let a = scalar
            .impulse_candidate(&held, &[(3, 0.2)], [0.0; 3], impulse_controls(), cx)
            .unwrap();
        let b = nodal
            .impulse_candidate(&held, &[(3, 0.2)], [0.0; 3], impulse_controls(), cx)
            .unwrap();
        assert_eq!(a.solid.velocity_mps(), b.solid.velocity_mps());
        assert_eq!(
            a.impulse.kinetic_change_j.to_bits(),
            b.impulse.kinetic_change_j.to_bits()
        );
        assert_eq!(
            a.impulse.support_impulses_kg_m_s,
            b.impulse.support_impulses_kg_m_s
        );
        assert!(scalar.reference_nodal_density_kg_m3().is_none());
        assert_eq!(
            nodal.reference_nodal_density_kg_m3(),
            Some([6.0; 4].as_slice())
        );
        for rho in [
            &[][..],
            &[6.0; 3][..],
            &[6.0, 0.0, 6.0, 6.0][..],
            &[f64::NAN; 4][..],
        ] {
            assert!(initial(cx).with_nodal_inertia(rho, &velocity, cx).is_err());
        }
        assert!(
            initial(cx)
                .with_nodal_inertia(&RHO, &velocity[..3], cx)
                .is_err()
        );
        assert!(
            initial(cx)
                .with_nodal_inertia(&RHO, &[[f64::INFINITY; 3]; 4], cx)
                .is_err()
        );
        assert!(
            initial(cx)
                .with_inertia(&[6.0], &velocity, cx)
                .unwrap()
                .with_nodal_inertia(&RHO, &velocity, cx)
                .is_err()
        );
        assert!(
            initial(cx)
                .reference_acceleration_load_n(&acceleration, cx)
                .is_err()
        );
        assert!(nodal.reference_acceleration_load_n(&[], cx).is_err());
        assert!(
            nodal
                .reference_acceleration_load_n(&[[f64::NAN; 3]; 4], cx)
                .is_err()
        );

        let unbound = initial(cx);
        let gate = CancelGate::new_clock_free();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(
                unbound
                    .with_nodal_inertia(&RHO, &velocity, &cancelled)
                    .is_err()
            );
            assert!(
                nodal
                    .reference_acceleration_load_n(&acceleration, &cancelled)
                    .is_err()
            );
            assert!(
                nodal
                    .impulse_candidate(&held, &[(3, 0.2)], [0.0; 3], impulse_controls(), &cancelled)
                    .is_err()
            );
        });
        assert_eq!(nodal.velocity_mps(), Some(velocity.as_slice()));
        assert_eq!(
            nodal.reference_nodal_density_kg_m3(),
            Some([6.0; 4].as_slice())
        );
    });
}
