//! G0/G1/G3/G4: density-consistent body loads on the native shell trajectory.
use super::*;

const DENSITY: [[f64; 3]; 2] = [[2.0, 7.0, 12.0], [9.0, 20.0, 14.0]];
// On each area-1/2, thickness-1/10 triangle,
// integral rho*h*N_i dA = (sum(rho) + rho_i)/240.
// Incident triangles deliberately disagree at their shared corners.
const NODE_MASS: [f64; 4] = [23.0 / 240.0, 80.0 / 240.0, 90.0 / 240.0, 63.0 / 240.0];
const SHELL_MASS: f64 = 16.0 / 15.0;

fn body() -> ShellRigidPayload {
    ShellRigidPayload {
        point: ShellLoadPoint {
            triangle: 0,
            barycentric: [0.2, 0.3, 0.5],
            offset_m: [0.14, -0.08, 0.23],
        },
        mass_kg: 2.3,
        inertia_com_kg_m2: [
            [0.006, 0.001, -0.0004],
            [0.001, 0.008, 0.0007],
            [-0.0004, 0.0007, 0.010],
        ],
    }
}

fn model(formulation: ShellMassFormulation, max_points: usize) -> MovingShellModel {
    MovingShellModel::new_with_payloads_and_reference_density(
        plate(),
        DENSITY.to_vec(),
        formulation,
        max_points,
        &[body()],
        || false,
    )
    .unwrap()
}

fn vertex(node: usize) -> ShellLoadPoint {
    let (triangle, barycentric) = match node {
        0 => (0, [1.0, 0.0, 0.0]),
        1 => (0, [0.0, 1.0, 0.0]),
        2 => (0, [0.0, 0.0, 1.0]),
        3 => (1, [0.0, 1.0, 0.0]),
        _ => panic!("fixture node"),
    };
    ShellLoadPoint {
        triangle,
        barycentric,
        offset_m: [0.0; 3],
    }
}

// Independent physical loads, without inspecting the assembled mass.
fn reference_wrenches(acceleration: [f64; 3]) -> Vec<ShellPointWrench> {
    let mut loads: Vec<_> = NODE_MASS
        .into_iter()
        .enumerate()
        .map(|(node, mass)| ShellPointWrench {
            point: vertex(node),
            force_n: acceleration.map(|a| mass * a),
            couple_nm: [0.0; 3],
        })
        .collect();
    loads.push(ShellPointWrench {
        point: body().point,
        force_n: acceleration.map(|a| body().mass_kg * a),
        couple_nm: [0.0; 3],
    });
    loads
}

#[test]
fn graded_body_loads_preserve_mass_centroid_payload_moment_and_potential() {
    let acceleration = [0.7, -1.2, -9.81];
    let payload = body();
    let payload_force = acceleration.map(|a| payload.mass_kg * a);
    let payload_couple = cross(payload.point.offset_m, payload_force);
    for formulation in [
        ShellMassFormulation::Lumped,
        ShellMassFormulation::ConsistentReference,
    ] {
        for supported in [false, true] {
            let mut p = plate();
            if supported {
                p.support = Some(ShellSupport {
                    node_indices: [2, 0, 1],
                    normal: [0.0, 0.0, 1.0],
                });
            }
            // Four distributed node forces and a payload need no point entries.
            let model = MovingShellModel::new_with_payloads_and_reference_density(
                p,
                DENSITY.to_vec(),
                formulation,
                1,
                &[payload],
                || false,
            )
            .unwrap();
            let load = model
                .project_with_body_acceleration(&[], acceleration, || false)
                .unwrap();
            let mut expected = vec![0.0; 24];
            let mut expected_moment = [0.0; 3];
            let mut payload_position = payload.point.offset_m;
            for (node, mass) in NODE_MASS.into_iter().enumerate() {
                let force = acceleration.map(|a| mass * a);
                let moment = cross(model.plate().nodes[node].position_m, force);
                for axis in 0..3 {
                    expected[6 * node + axis] = force[axis];
                    expected_moment[axis] += moment[axis];
                }
            }
            for (node, weight) in [0, 1, 2].into_iter().zip(payload.point.barycentric) {
                for axis in 0..3 {
                    expected[6 * node + axis] += weight * payload_force[axis];
                    expected[6 * node + 3 + axis] += weight * payload_couple[axis];
                    payload_position[axis] += weight * model.plate().nodes[node].position_m[axis];
                }
            }
            let orbital = cross(payload_position, payload_force);
            for axis in 0..3 {
                near(
                    load.resultant.linear[axis],
                    (SHELL_MASS + payload.mass_kg) * acceleration[axis],
                    2e-14,
                );
                near(
                    load.resultant.angular[axis],
                    expected_moment[axis] + orbital[axis],
                    3e-14,
                );
            }
            for (actual, expected) in load.full.iter().zip(&expected) {
                near(*actual, *expected, 5e-15);
            }
            for (&actual, &row) in load.reduced.iter().zip(&model.assembly().free_dofs) {
                near(actual, expected[row], 5e-15);
            }

            // Work is minus the variation of the independent shell/payload
            // gravitational potential, including director motion of the COM.
            let displacement: Vec<_> = (0..24).map(|i| 0.001 * (0.7 * i as f64).sin()).collect();
            let mut potential_change = 0.0;
            for (node, mass) in NODE_MASS.into_iter().enumerate() {
                let translation = std::array::from_fn(|a| displacement[6 * node + a]);
                potential_change -= mass * dot(acceleration, translation);
            }
            let mut com_displacement = [0.0; 3];
            for (node, weight) in [0, 1, 2].into_iter().zip(payload.point.barycentric) {
                let rotation = std::array::from_fn(|a| displacement[6 * node + 3 + a]);
                let arm_motion = cross(rotation, payload.point.offset_m);
                for axis in 0..3 {
                    com_displacement[axis] +=
                        weight * (displacement[6 * node + axis] + arm_motion[axis]);
                }
            }
            potential_change -= dot(payload_force, com_displacement);
            let work: f64 = load
                .full
                .iter()
                .zip(&displacement)
                .map(|(f, u)| f * u)
                .sum();
            near(work, -potential_change, 2e-17);

            let rigid_acceleration: Vec<_> = (0..24)
                .map(|i| if i % 6 < 3 { acceleration[i % 6] } else { 0.0 })
                .collect();
            let mass_force = model.assembly().full_mass.apply(&rigid_acceleration);
            for (actual, mass_force) in load.full.iter().zip(mass_force) {
                near(*actual, mass_force, 5e-15);
            }
            let point_load = wrench(point(1, [0.2, 0.3, 0.5]));
            let point_only = model.project(&[point_load], || false).unwrap();
            let combined = model
                .project_with_body_acceleration(&[point_load], acceleration, || false)
                .unwrap();
            for ((combined, body), point) in
                combined.full.iter().zip(&load.full).zip(point_only.full)
            {
                near(*combined, body + point, 2e-15);
            }
            assert!(
                model
                    .project_with_body_acceleration(&[point_load, point_load], acceleration, || {
                        false
                    })
                    .is_err()
            );
        }
    }
}

#[test]
fn equal_total_density_profiles_keep_distinct_centers_of_gravity() {
    let mut p = plate();
    p.nodes.truncate(3);
    p.triangles.truncate(1);
    for formulation in [
        ShellMassFormulation::Lumped,
        ShellMassFormulation::ConsistentReference,
    ] {
        let mut moments = Vec::new();
        for (density, expected_x_moment) in [
            ([1.0, 2.0, 3.0], 13.0 / 8.0),
            ([3.0, 2.0, 1.0], 41.0 / 24.0),
        ] {
            let model = MovingShellModel::new_with_reference_density(
                p.clone(),
                vec![density],
                formulation,
                1,
                || false,
            )
            .unwrap();
            let load = model
                .project_with_body_acceleration(&[], [0.0, 0.0, -10.0], || false)
                .unwrap();
            // Both bodies weigh 0.1 kg. Their barycentric centers are
            // [7,8,9]/24 and [9,8,7]/24, so only the y coordinate changes.
            near(load.resultant.linear[2], -1.0, 3e-16);
            near(load.resultant.angular[0], expected_x_moment, 5e-16);
            near(load.resultant.angular[1], 4.0 / 3.0, 5e-16);
            moments.push(load.resultant.angular[0]);
        }
        assert!(moments[1] > moments[0]);
    }
}

#[test]
fn graded_shell_and_offset_payload_fall_together_with_native_work_accounting() {
    let acceleration = [0.3, -0.8, -9.81];
    let total_mass = SHELL_MASS + body().mass_kg;
    let dt = config().dt_s;
    for formulation in [
        ShellMassFormulation::Lumped,
        ShellMassFormulation::ConsistentReference,
    ] {
        let model = model(formulation, 1);
        let zero = vec![0.0; model.assembly().mass.dimension()];
        let mut motion = model
            .start_with_body_acceleration(&zero, &zero, &[], acceleration, config(), || false)
            .unwrap();
        let mut work = 0.0;
        for _ in 0..12 {
            let report = motion
                .advance_with_body_acceleration(
                    motion.state().time_s + dt,
                    &[],
                    acceleration,
                    || false,
                )
                .unwrap();
            work += report.applied_work_j;
            let t = motion.state().time_s;
            for node in 0..4 {
                for axis in 0..3 {
                    near(
                        motion.state().acceleration[6 * node + axis],
                        acceleration[axis],
                        3e-10,
                    );
                    near(
                        motion.state().velocity[6 * node + axis],
                        acceleration[axis] * t,
                        3e-13,
                    );
                    near(
                        motion.state().displacement[6 * node + axis],
                        0.5 * acceleration[axis] * t * t,
                        3e-16,
                    );
                    near(motion.state().velocity[6 * node + 3 + axis], 0.0, 3e-13);
                    near(motion.state().displacement[6 * node + 3 + axis], 0.0, 3e-16);
                }
            }
            for axis in 0..3 {
                near(
                    report.applied_impulse.linear[axis],
                    total_mass * acceleration[axis] * dt,
                    2e-17,
                );
                near(report.momentum_residual.linear[axis], 0.0, 2e-14);
                near(report.momentum_residual.angular[axis], 0.0, 3e-14);
            }
            let kinetic = 0.5 * total_mass * dot(acceleration, acceleration) * t * t;
            near(motion.state().total_energy_j, kinetic, 3e-14);
            near(work, kinetic, 3e-14);
            near(report.energy_residual_j, 0.0, 3e-14);
            assert_eq!(report.support_work_j, 0.0);
        }
        // An independently supplied endpoint changes the load; acceleration is
        // neither inferred nor frozen at its initial value.
        let old_time = motion.state().time_s;
        let next_acceleration = [-0.4, 0.2, -3.0];
        let report = motion
            .advance_with_body_acceleration(old_time + dt, &[], next_acceleration, || false)
            .unwrap();
        work += report.applied_work_j;
        let expected_velocity = std::array::from_fn(|axis| {
            acceleration[axis] * old_time
                + 0.5 * dt * (acceleration[axis] + next_acceleration[axis])
        });
        for node in 0..4 {
            for axis in 0..3 {
                let expected_displacement = 0.5 * acceleration[axis] * old_time * old_time
                    + dt * acceleration[axis] * old_time
                    + 0.25 * dt * dt * (acceleration[axis] + next_acceleration[axis]);
                near(
                    motion.state().acceleration[6 * node + axis],
                    next_acceleration[axis],
                    3e-10,
                );
                near(
                    motion.state().velocity[6 * node + axis],
                    expected_velocity[axis],
                    3e-13,
                );
                near(
                    motion.state().displacement[6 * node + axis],
                    expected_displacement,
                    3e-16,
                );
            }
        }
        near(
            work,
            0.5 * total_mass * dot(expected_velocity, expected_velocity),
            3e-14,
        );
    }
}

#[test]
fn stationary_reactions_and_falling_driven_pins_use_the_complete_inventory() {
    let acceleration = [0.0, 0.0, -9.81];
    let mut p = plate();
    p.nodes.truncate(3);
    p.triangles.truncate(1);
    p.support = Some(ShellSupport {
        node_indices: [2, 0, 1],
        normal: [0.0, 0.0, 1.0],
    });
    let density = [2.0, 7.0, 12.0];
    let total_mass = 0.35 + body().mass_kg;
    for formulation in [
        ShellMassFormulation::Lumped,
        ShellMassFormulation::ConsistentReference,
    ] {
        let mut centered = body();
        centered.point.offset_m = [0.0, 0.0, 0.23];
        let stationary = MovingShellModel::new_with_payloads_and_reference_density(
            p.clone(),
            vec![density],
            formulation,
            1,
            &[centered],
            || false,
        )
        .unwrap();
        let zero = vec![0.0; stationary.assembly().mass.dimension()];
        let load = stationary
            .project_with_body_acceleration(&[], acceleration, || false)
            .unwrap();
        let mut motion = stationary
            .start_with_body_acceleration(&zero, &zero, &[], acceleration, config(), || false)
            .unwrap();
        let report = motion
            .advance_with_body_acceleration(config().dt_s, &[], acceleration, || false)
            .unwrap();
        for (reaction, force) in motion.support_reactions().iter().zip(&load.full) {
            near(*reaction, -force, 3e-15);
        }
        assert_eq!(motion.state().displacement, zero);
        assert_eq!(motion.state().velocity, zero);
        near(
            report.support_impulse.linear[2],
            -total_mass * acceleration[2] * config().dt_s,
            5e-18,
        );
        assert_eq!(report.applied_work_j, 0.0);

        // Offset COM moments cancel the corresponding prescribed-inertia term
        // when all three pins fall at the same acceleration as the body load.
        let falling = MovingShellModel::new_with_payloads_and_reference_density(
            p.clone(),
            vec![density],
            formulation,
            1,
            &[body()],
            || false,
        )
        .unwrap();
        let supports = ShellSupportMotion {
            acceleration_m_s2: [acceleration; 3],
            ..Default::default()
        };
        let mut motion = falling
            .start_driven_with_body_acceleration(
                &zero,
                &zero,
                supports,
                &[],
                acceleration,
                config(),
                || false,
            )
            .unwrap();
        let mut work = 0.0;
        for _ in 0..10 {
            let report = motion
                .advance_driven_with_body_acceleration(
                    motion.state().time_s + config().dt_s,
                    &[],
                    acceleration,
                    [acceleration; 3],
                    || false,
                )
                .unwrap();
            work += report.applied_work_j;
            let time = motion.state().time_s;
            for reaction in motion.support_reactions() {
                near(*reaction, 0.0, 5e-13);
            }
            for velocity in &motion.state().velocity {
                near(*velocity, 0.0, 3e-14);
            }
            for row in motion.support_motion().unwrap().velocity_m_s {
                for axis in 0..3 {
                    near(row[axis], acceleration[axis] * time, 2e-17);
                }
            }
            near(report.support_work_j, 0.0, 2e-17);
            near(report.energy_residual_j, 0.0, 2e-16);
            near(
                motion.state().total_energy_j,
                0.5 * total_mass * dot(acceleration, acceleration) * time * time,
                2e-16,
            );
        }
        near(work, motion.state().total_energy_j, 2e-16);
    }
}

#[test]
fn held_body_force_survives_impulses_and_matches_independent_physical_ports() {
    let acceleration = [0.3, -0.8, -9.81];
    let next_acceleration = [0.2, -0.3, -8.0];
    for formulation in [
        ShellMassFormulation::Lumped,
        ShellMassFormulation::ConsistentReference,
    ] {
        let mut p = plate();
        p.damping = DampingModel::Rayleigh {
            mass_proportional_per_s: 0.3,
            stiffness_proportional_s: 0.0001,
        };
        let model = MovingShellModel::new_with_payloads_and_reference_density(
            p,
            DENSITY.to_vec(),
            formulation,
            5,
            &[body()],
            || false,
        )
        .unwrap();
        let zero = vec![0.0; model.assembly().mass.dimension()];
        let mut native = model
            .start_with_body_acceleration(&zero, &zero, &[], acceleration, config(), || false)
            .unwrap();
        let mut ports = model
            .start(
                &zero,
                &zero,
                &reference_wrenches(acceleration),
                config(),
                || false,
            )
            .unwrap();
        let impulses = [ShellPointImpulse {
            point: vertex(0),
            linear_impulse_ns: [0.001, -0.002, 0.003],
            angular_impulse_nms: [0.00002, 0.0, -0.00001],
        }];
        assert_eq!(
            native.apply_impulses(&impulses, || false).unwrap(),
            ports.apply_impulses(&impulses, || false).unwrap(),
        );
        let state = native.state();
        let inertia = model.assembly().mass.apply(&state.acceleration);
        let elastic = model.assembly().stiffness.apply(&state.displacement);
        let viscous = model
            .assembly()
            .damping
            .as_ref()
            .unwrap()
            .apply(&state.velocity);
        let held = model
            .project(&reference_wrenches(acceleration), || false)
            .unwrap();
        for i in 0..inertia.len() {
            near(inertia[i] + elastic[i] + viscous[i], held.reduced[i], 3e-12);
        }
        let report = native
            .advance_with_body_acceleration(config().dt_s, &[], next_acceleration, || false)
            .unwrap();
        let reference = ports
            .advance(
                config().dt_s,
                &reference_wrenches(next_acceleration),
                || false,
            )
            .unwrap();
        for (a, b) in [
            (&native.state().displacement, &ports.state().displacement),
            (&native.state().velocity, &ports.state().velocity),
            (&native.state().acceleration, &ports.state().acceleration),
        ] {
            for (&a, &b) in a.iter().zip(b) {
                near(a, b, 3e-10);
            }
        }
        near(report.applied_work_j, reference.applied_work_j, 3e-16);
        near(report.energy_residual_j, 0.0, 3e-14);
        for axis in 0..3 {
            near(
                report.applied_impulse.linear[axis],
                (SHELL_MASS + body().mass_kg)
                    * 0.5
                    * config().dt_s
                    * (acceleration[axis] + next_acceleration[axis]),
                3e-17,
            );
        }
        // The original endpoint API explicitly removes body acceleration.
        // Its impulse still includes the retained previous body-force endpoint.
        let removed = native.advance(2.0 * config().dt_s, &[], || false).unwrap();
        for axis in 0..3 {
            near(
                removed.applied_impulse.linear[axis],
                (SHELL_MASS + body().mass_kg) * 0.5 * config().dt_s * next_acceleration[axis],
                3e-17,
            );
        }
    }
}

#[test]
fn zero_body_acceleration_preserves_bits_and_refused_endpoints_are_retryable() {
    let model = model(ShellMassFormulation::ConsistentReference, 1);
    let zero = vec![0.0; model.assembly().mass.dimension()];
    let loads = [wrench(point(1, [0.2, 0.3, 0.5]))];
    assert_eq!(
        model.project(&loads, || false).unwrap(),
        model
            .project_with_body_acceleration(&loads, [0.0; 3], || false)
            .unwrap(),
    );
    let mut original = model
        .start(&zero, &zero, &loads, config(), || false)
        .unwrap();
    let mut explicit_zero = model
        .start_with_body_acceleration(&zero, &zero, &loads, [0.0; 3], config(), || false)
        .unwrap();
    assert_eq!(original.state(), explicit_zero.state());
    assert_eq!(
        original.support_reactions(),
        explicit_zero.support_reactions()
    );
    assert_eq!(
        original.advance(config().dt_s, &loads, || false).unwrap(),
        explicit_zero
            .advance_with_body_acceleration(config().dt_s, &loads, [0.0; 3], || false)
            .unwrap(),
    );
    assert_eq!(original.state(), explicit_zero.state());

    let acceleration = [0.0, 0.0, -9.81];
    let next_acceleration = [0.2, -0.3, -8.0];
    let mut motion = model
        .start_with_body_acceleration(&zero, &zero, &loads, acceleration, config(), || false)
        .unwrap();
    let saved = motion.state().clone();
    let reactions = motion.support_reactions().to_vec();
    for bad in [
        [f64::NAN, 0.0, 0.0],
        [0.0, f64::INFINITY, 0.0],
        [f64::MAX; 3],
    ] {
        assert!(
            motion
                .advance_with_body_acceleration(config().dt_s, &loads, bad, || false)
                .is_err()
        );
        assert_eq!(motion.state(), &saved);
        assert_eq!(motion.support_reactions(), reactions);
    }
    let mut control = model
        .start_with_body_acceleration(&zero, &zero, &loads, acceleration, config(), || false)
        .unwrap();
    let mut total_polls = 0;
    let expected = control
        .advance_with_body_acceleration(config().dt_s, &loads, next_acceleration, || {
            total_polls += 1;
            false
        })
        .unwrap();
    for stop in [1, total_polls / 2, total_polls] {
        let mut polls = 0;
        assert!(matches!(
            motion.advance_with_body_acceleration(config().dt_s, &loads, next_acceleration, || {
                polls += 1;
                polls == stop
            }),
            Err(MovingLoadError::Cancelled),
        ));
        assert_eq!(motion.state(), &saved);
        assert_eq!(motion.support_reactions(), reactions);
    }
    assert_eq!(
        motion
            .advance_with_body_acceleration(config().dt_s, &loads, next_acceleration, || false)
            .unwrap(),
        expected,
    );
    assert_eq!(motion.state(), control.state());
    assert_eq!(motion.support_reactions(), control.support_reactions());
}
