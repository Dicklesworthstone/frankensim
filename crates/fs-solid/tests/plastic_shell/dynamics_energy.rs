use super::*;
use fs_material::plastic::rate::LinearOverstressJ2;
use fs_solid::plastic_shell::dynamics::PlasticShellDensity;

fn triangle() -> ShellMesh {
    ShellMesh::new(
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}

fn bending_endpoint(mesh: &ShellMesh, curvature: f64) -> Vec<(usize, f64)> {
    mesh.nodes
        .iter()
        .flat_map(|&[x, _, _]| [0.0, 0.0, 0.5 * curvature * x * x, 0.0, -curvature * x, 0.0])
        .enumerate()
        .collect()
}

#[test]
fn bending_motion_closes_native_plastic_heat_kinetics_and_unloading_work() {
    with_cx(|_, cx| {
        let mesh = triangle();
        let material = section(metal());
        let owner = problem(&mesh, &material, &[], &[], &[1.0], 0.37);
        let density = 2.0;
        let dt = 0.25;
        let reference = owner
            .prepare_reference_inertia(PlasticShellDensity::Uniform(density), cx)
            .unwrap();
        let mut previous: Option<SmallStrainShellSolution<J2State>> = None;
        let mut velocity = vec![[0.0; 6]; mesh.nodes.len()];
        let mut previous_storage = 0.0;
        let mut previous_kinetic = 0.0;
        let mut cumulative_heat = 0.0;
        for factor in [1.0, 0.5, 0.0] {
            let curvature = 0.04 * factor;
            let prescribed = bending_endpoint(&mesh, curvature);
            let step = problem(&mesh, &material, &prescribed, &[], &[1.0], 0.37)
                .advance_backward_euler_with_energy(
                    previous.as_ref(),
                    &velocity,
                    &reference,
                    dt,
                    gates(),
                    cx,
                )
                .unwrap();
            let mut yield_heat = 0.0;
            for q in 0..3 {
                for (i, point) in material.points().iter().enumerate() {
                    let before = previous
                        .as_ref()
                        .map_or(0.0, |old| old.section_states[q].points[i].history.alpha);
                    let after = step.solution.section_states[q].points[i].history.alpha;
                    yield_heat += (0.5 / 3.0)
                        * point.weight_m
                        * material.layers()[point.layer].law.yield_stress
                        * (after - before);
                }
            }
            let e = &step.energy;
            let physical = &step.inertia.physical;
            let gauge = &step.inertia.drilling_gauge;
            near(e.physical_dissipation_j, yield_heat, 1e-12);
            near(e.initial.total_j(), previous_storage, 1e-11);
            near(physical.kinetic_before_j, previous_kinetic, 1e-14);
            near(gauge.kinetic_after_j, 0.0, 1e-25);
            near(e.balance_residual_j, 0.0, 1e-9);
            near(step.balance_residual_j(), 0.0, 1e-9);
            near(
                e.dead_load_endpoint_work_j
                    + e.support_endpoint_work_j
                    + e.transverse_residual_work_j,
                e.stored_change_j()
                    + e.physical_dissipation_j
                    + e.material_numerical_dissipation_j
                    + e.drilling_final_j
                    - e.drilling_initial_j
                    + e.drilling_numerical_loss_j
                    + physical.kinetic_change_j
                    + physical.numerical_loss_j
                    + gauge.kinetic_change_j
                    + gauge.numerical_loss_j,
                1e-9,
            );
            if previous.is_none() {
                // Only corner (1,0) moves. Exact P1 diagonal integral is A/6;
                // symmetric uniform density has no translation/director cross term.
                let expected = 0.5
                    * (density * H * 0.5 / 6.0)
                    * ((0.5 * curvature / dt).powi(2) + H * H / 12.0 * (curvature / dt).powi(2));
                near(physical.kinetic_after_j, expected, 1e-15);
                near(physical.numerical_loss_j, expected, 1e-15);
                let points = &step.solution.section_states[0].points;
                assert!(points.iter().any(|p| p.history.alpha > 0.0));
                assert!(points.iter().any(|p| p.history.alpha == 0.0));
            }
            cumulative_heat += e.physical_dissipation_j;
            previous_storage = e.final_state.total_j();
            previous_kinetic = physical.kinetic_after_j;
            velocity = step.velocity;
            previous = Some(step.solution);
        }
        assert!(cumulative_heat > 0.0);
        assert!(
            previous
                .unwrap()
                .section_states
                .iter()
                .any(|state| { state.points.iter().any(|p| p.history.alpha > 0.0) })
        );
    });
}

#[test]
fn drilling_pulse_keeps_gauge_kinetics_and_damping_out_of_material_heat() {
    with_cx(|_, cx| {
        let mesh = triangle();
        let material = section(metal());
        let angle = 0.01;
        let dt = 0.25;
        let density = 2.0;
        let prescribed: Vec<_> = (0..18)
            .map(|d| (d, if d == 5 { angle } else { 0.0 }))
            .collect();
        let owner = problem(&mesh, &material, &prescribed, &[], &[1.0], 0.37);
        let reference = owner
            .prepare_reference_inertia(PlasticShellDensity::Uniform(density), cx)
            .unwrap();
        let step = owner
            .advance_backward_euler_with_energy(None, &[[0.0; 6]; 3], &reference, dt, gates(), cx)
            .unwrap();
        // K_gauge = .5 * rho*h^3/12 * integral N0^2 dA * angular_rate^2.
        let expected = density * H.powi(3) * 0.5 * (angle / dt).powi(2) / 144.0;
        let gauge = &step.inertia.drilling_gauge;
        near(gauge.kinetic_after_j, expected, 1e-20);
        near(gauge.numerical_loss_j, expected, 1e-20);
        near(gauge.inertial_work_j, 2.0 * expected, 1e-20);
        near(step.inertia.physical.kinetic_after_j, 0.0, 0.0);
        near(step.inertia.physical.numerical_loss_j, 0.0, 0.0);
        near(step.energy.physical_dissipation_j, 0.0, 0.0);
        near(step.energy.material_numerical_dissipation_j, 0.0, 0.0);
        near(step.energy.final_state.total_j(), 0.0, 0.0);
        assert!(step.energy.drilling_final_j > 0.0);
        near(
            step.energy.support_endpoint_work_j,
            2.0 * (step.energy.drilling_final_j + expected),
            1e-16,
        );
        near(step.balance_residual_j(), 0.0, 1e-16);
    });
}

#[test]
fn timed_native_material_must_share_the_actual_inertial_interval() {
    with_cx(|_, cx| {
        let mesh = triangle();
        let duration = 0.25;
        let rate = LinearOverstressJ2::new(metal(), 10.0, duration).unwrap();
        // A borrowed law must preserve the same native interval metadata.
        let material = section(&rate);
        let prescribed = bending_endpoint(&mesh, 0.04);
        let owner = problem(&mesh, &material, &prescribed, &[], &[1.0], 0.37);
        let reference = owner
            .prepare_reference_inertia(PlasticShellDensity::Uniform(2.0), cx)
            .unwrap();
        let velocity = [[0.0; 6]; 3];
        let error = owner
            .advance_backward_euler_with_energy(
                None,
                &velocity,
                &reference,
                2.0 * duration,
                gates(),
                cx,
            )
            .unwrap_err();
        assert!(error.to_string().contains("durations must agree"));
        let accepted = owner
            .advance_backward_euler_with_energy(None, &velocity, &reference, duration, gates(), cx)
            .unwrap();
        assert!(accepted.energy.physical_dissipation_j > 0.0);
        near(accepted.balance_residual_j(), 0.0, 1e-9);
    });
}
