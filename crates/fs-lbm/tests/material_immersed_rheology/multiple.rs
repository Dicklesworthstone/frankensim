//! Simultaneous material viscosity, shared-cell weights and assembly rollback.

use super::*;

fn separated_bodies() -> [MovingPartialSaturationBody3; 2] {
    let mut bodies = [body(); 2];
    bodies[0].center = [2.4, 4.0, 4.0];
    bodies[1].center = [5.6, 4.0, 4.0];
    bodies[0].velocity = [0.007, -0.002, 0.001];
    bodies[1].velocity = [-0.004, 0.003, -0.002];
    bodies[1].mass = 1.1;
    bodies
}

#[test]
fn constant_material_matches_one_native_assembly_for_all_collision_families() {
    let reference = transport(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let shape = BodyFrameEllipsoid3::sphere(0.8).unwrap();
    let gate = CancelGate::new_clock_free();
    for family in 0..4 {
        let force = [2e-6, -1e-6, 5e-7];
        let mut actual = grid(model(1.37, family), force);
        let mut expected = grid(model(scaling.flow().tau, family), force);
        seed_velocity(&mut actual, |_| [0.024, -0.017, 0.011]);
        seed_velocity(&mut expected, |_| [0.024, -0.017, 0.011]);
        let mut a = separated_bodies();
        let mut b = a;
        let mut sa = core::array::from_fn(|_| BodyFrameSampler3::new(&actual, 3).unwrap());
        let mut sb = core::array::from_fn(|_| BodyFrameSampler3::new(&expected, 3).unwrap());
        let external = [
            [1e-4, -2e-5, 0.0, 0.0, 1e-5, 0.0],
            [-3e-5, 0.0, 2e-5, 1e-5, 0.0, 0.0],
        ];
        let mut queries = 0;
        let report = scaling
            .step_moving_bodies_with_material_rheology(
                &mut actual,
                &mut sa,
                [&shape as &dyn BodyFrameShape3; 2],
                &mut a,
                external,
                [0.3, -0.4, 0.2],
                &mut workspace(),
                CELLS,
                None,
                &gate,
                |point| {
                    queries += 1;
                    Ok(SingleRead {
                        point,
                        value: reference,
                        reads: Cell::new(0),
                    })
                },
            )
            .unwrap();
        let native = expected
            .step_with_moving_partial_saturation_bodies(
                &mut sb,
                [&shape as &dyn BodyFrameShape3; 2],
                &mut b,
                external,
                [0.3, -0.4, 0.2],
                &gate,
            )
            .unwrap();
        assert_eq!(queries, CELLS);
        assert_eq!(report.applied.len(), CELLS);
        assert_eq!(report.momentum, native);
        assert_eq!(a, b);
        assert_eq!(snapshot(&actual), snapshot(&expected));
        assert_eq!(actual.collision_model(), model(1.37, family));
        for sample in report.applied {
            assert_eq!(sample.material.reads.get(), 1);
            assert_eq!(sample.material.point, sample.point);
        }
    }
}

#[test]
fn shared_cell_material_weight_closes_two_independent_body_responses_and_fluid_momentum() {
    let reference = transport(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let gate = CancelGate::new_clock_free();
    let shape = BodyFrameEllipsoid3::sphere(0.33).unwrap();
    let force = [2e-6, 0.0, 0.0];
    let mut fluid = grid(model(1.37, 2), force);
    let speed = |cell: [usize; 3]| 0.001 * (cell[1] as f64 + 0.5 - 4.0).powi(2);
    seed_velocity(&mut fluid, |cell| [speed(cell), 0.0, 0.0]);
    let before = snapshot(&fluid);
    let mut bodies = [body(); 2];
    bodies[0].center = [3.5; 3];
    bodies[1].center = [4.2, 3.5, 3.5];
    bodies[1].mass = 0.9;
    let initial = bodies;
    let mut samplers: [BodyFrameSampler3; 2] =
        core::array::from_fn(|_| BodyFrameSampler3::new(&fluid, 8).unwrap());
    for a in 0..2 {
        samplers[a]
            .sample(bodies[a].center, bodies[a].orientation, &shape, &gate)
            .unwrap();
    }
    let fields = samplers.each_ref().map(|sampler| sampler.fractions());
    assert!(
        fields[0]
            .iter()
            .zip(fields[1])
            .any(|(&a, &b)| a > 0.0 && b > 0.0)
    );
    let viscosity = |cell: [usize; 3], rate: f64| carreau(rate, 7.0) * (1.0 + 0.1 * cell[0] as f64);
    let mut resistance = [0.0; 2];
    let mut drive = [0.0; 2];
    let mut isolated = [0.0; 2];
    for i in 0..CELLS {
        let cell = [i % N, i / N % N, i / (N * N)];
        let rate = (0.002 * (cell[1] as f64 + 0.5 - 4.0)).abs() / DT;
        let tau = 0.5 + 3.0 * viscosity(cell, rate) / 1000.0 * DT / (DX * DX);
        let total = fields[0][i] + fields[1][i];
        for a in 0..2 {
            let epsilon = fields[a][i];
            if epsilon == 0.0 {
                continue;
            }
            // These separated spheres occupy distinct subvolumes of a shared
            // voxel. Both centers and all occupied cells lie on the x axis
            // through y=z=3.5, so this physical drive has zero rigid torque.
            assert_eq!([cell[1], cell[2]], [3, 3]);
            let weight = epsilon * (tau - 0.5) / (1.0 - total + tau - 0.5);
            resistance[a] += weight;
            drive[a] += weight * speed(cell);
            isolated[a] += epsilon * (tau - 0.5) / (1.0 - epsilon + tau - 0.5);
        }
    }
    assert!(resistance[0] - isolated[0] > 1e-6);
    assert!(resistance[1] - isolated[1] > 1e-6);
    let expected =
        core::array::from_fn::<_, 2, _>(|a| drive[a] / (initial[a].mass + 0.5 * resistance[a]));
    let mut queries = 0;
    let report = scaling
        .step_moving_bodies_with_material_rheology(
            &mut fluid,
            &mut samplers,
            [&shape as &dyn BodyFrameShape3; 2],
            &mut bodies,
            [[0.0; 6]; 2],
            [0.0; 3],
            &mut workspace(),
            CELLS,
            None,
            &gate,
            |point| {
                queries += 1;
                Ok(transport(viscosity(point.cell, point.shear_rate_per_s)))
            },
        )
        .unwrap();
    assert_eq!(queries, CELLS);
    for a in 0..2 {
        close(bodies[a].velocity[0], expected[a], 3e-14);
        close(
            bodies[a].center[0],
            initial[a].center[0] + 0.5 * expected[a],
            2e-14,
        );
        close(
            report.momentum.coupling.bodies[a]
                .hydrodynamic
                .coupling_weight_sum,
            resistance[a],
            2e-12,
        );
        for value in bodies[a].velocity[1..]
            .iter()
            .copied()
            .chain(bodies[a].angular_velocity_world().unwrap())
        {
            close(value, 0.0, 3e-14);
        }
        close(report.momentum.bodies[a].work_residual, 0.0, 2e-14);
    }
    let after = snapshot(&fluid);
    let mut mass = 0.0;
    let mut momentum = [0.0; 3];
    for (old, new) in before.iter().zip(after) {
        for q in 0..Q3 {
            let delta = new[q] - old[q];
            mass += delta;
            let direction = [E3[q].0, E3[q].1, E3[q].2].map(f64::from);
            for axis in 0..3 {
                momentum[axis] += delta * direction[axis];
            }
        }
    }
    close(mass, 0.0, 3e-12);
    for axis in 0..3 {
        let body_change: f64 = (0..2)
            .map(|a| initial[a].mass * bodies[a].velocity[axis])
            .sum();
        close(
            momentum[axis] + body_change,
            CELLS as f64 * force[axis],
            4e-12,
        );
    }
}

#[test]
fn single_member_material_assembly_keeps_original_single_body_population_history() {
    let reference = transport(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let shape = BodyFrameEllipsoid3::sphere(1.13).unwrap();
    let mut assembly = grid(model(1.37, 3), [2e-6, -1e-6, 0.0]);
    let mut single = grid(model(1.37, 3), [2e-6, -1e-6, 0.0]);
    seed_velocity(&mut assembly, |cell| {
        [0.002 * (cell[1] as f64 - 4.0), 0.0, 0.0]
    });
    seed_velocity(&mut single, |cell| {
        [0.002 * (cell[1] as f64 - 4.0), 0.0, 0.0]
    });
    let mut bodies = [body()];
    let mut carrier = bodies[0];
    let mut samplers = [BodyFrameSampler3::new(&assembly, 3).unwrap()];
    let mut sampler = BodyFrameSampler3::new(&single, 3).unwrap();
    let gate = CancelGate::new_clock_free();
    let a = scaling
        .step_moving_bodies_with_material_rheology(
            &mut assembly,
            &mut samplers,
            [&shape as &dyn BodyFrameShape3],
            &mut bodies,
            [[0.0; 6]],
            [0.0; 3],
            &mut workspace(),
            CELLS,
            None,
            &gate,
            |point| Ok(transport(carreau(point.shear_rate_per_s, 4.0))),
        )
        .unwrap();
    let b = scaling
        .step_moving_body_with_material_rheology(
            &mut single,
            &mut sampler,
            &shape,
            &mut carrier,
            [0.0; 6],
            [0.0; 3],
            &mut workspace(),
            CELLS,
            None,
            &gate,
            |point| Ok(transport(carreau(point.shear_rate_per_s, 4.0))),
        )
        .unwrap();
    assert_eq!(bodies[0], carrier);
    assert_eq!(snapshot(&assembly), snapshot(&single));
    assert_eq!(a.applied, b.applied);
    assert_eq!(
        a.momentum.coupling.bodies[0].hydrodynamic,
        b.momentum.coupling.hydrodynamic
    );
    assert_eq!(a.momentum.coupling.forcing, b.momentum.coupling.forcing);
    assert_eq!(a.momentum.coupling.walls, b.momentum.coupling.walls);
    assert_eq!(a.momentum.coupling.periodic, b.momentum.coupling.periodic);
}

#[test]
fn late_assembly_geometry_cancellation_preserves_all_bodies_and_retries_once() {
    let reference = transport(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let mut fluid = grid(model(1.37, 1), [0.0; 3]);
    let mut fresh = grid(model(1.37, 1), [0.0; 3]);
    let saved = snapshot(&fluid);
    let mut bodies = [body(); 2];
    bodies[0].center = [2.0, 4.0, 4.0];
    bodies[1].center = [6.0, 4.0, 4.0];
    let initial = bodies;
    let mut expected = bodies;
    let mut samplers = core::array::from_fn(|_| BodyFrameSampler3::new(&fluid, 2).unwrap());
    let mut fresh_samplers = core::array::from_fn(|_| BodyFrameSampler3::new(&fresh, 2).unwrap());
    let shape = BodyFrameEllipsoid3::sphere(1.0).unwrap();
    let gate = CancelGate::new_clock_free();
    let cancelling = CancelOnEndpoint {
        gate: &gate,
        calls: Cell::new(0),
    };
    let error = scaling
        .step_moving_bodies_with_material_rheology(
            &mut fluid,
            &mut samplers,
            [&cancelling as &dyn BodyFrameShape3, &shape],
            &mut bodies,
            [[0.0; 6]; 2],
            [0.0; 3],
            &mut workspace(),
            CELLS,
            None,
            &gate,
            |_| Ok(reference),
        )
        .unwrap_err();
    assert!(matches!(
        error,
        MaterialImmersedError3::Flow(PartialSaturationError3::Cancelled)
    ));
    assert!(cancelling.calls.get() >= 129);
    assert_eq!(bodies, initial);
    assert_eq!(snapshot(&fluid), saved);
    let gate = CancelGate::new_clock_free();
    let retry = scaling
        .step_moving_bodies_with_material_rheology(
            &mut fluid,
            &mut samplers,
            [&shape as &dyn BodyFrameShape3; 2],
            &mut bodies,
            [[0.0; 6]; 2],
            [0.0; 3],
            &mut workspace(),
            CELLS,
            None,
            &gate,
            |_| Ok(reference),
        )
        .unwrap();
    let control = scaling
        .step_moving_bodies_with_material_rheology(
            &mut fresh,
            &mut fresh_samplers,
            [&shape as &dyn BodyFrameShape3; 2],
            &mut expected,
            [[0.0; 6]; 2],
            [0.0; 3],
            &mut workspace(),
            CELLS,
            None,
            &gate,
            |_| Ok(reference),
        )
        .unwrap();
    assert_eq!(bodies, expected);
    assert_eq!(snapshot(&fluid), snapshot(&fresh));
    assert_eq!(retry, control);
}
