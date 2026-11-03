//! Physical point inverse inertia without a contact law or surrogate solve.
use super::*;

fn payload() -> ShellRigidPayload {
    ShellRigidPayload {
        point: point(0, [0.2, 0.3, 0.5]),
        mass_kg: 0.7,
        inertia_com_kg_m2: [
            [0.002, 0.0001, 0.0],
            [0.0001, 0.003, -0.0002],
            [0.0, -0.0002, 0.004],
        ],
    }
}

fn vertex(node: usize) -> ShellLoadPoint {
    let (triangle, barycentric) = match node {
        0 => (0, [1.0, 0.0, 0.0]),
        1 => (0, [0.0, 1.0, 0.0]),
        2 => (0, [0.0, 0.0, 1.0]),
        3 => (1, [0.0, 1.0, 0.0]),
        _ => unreachable!(),
    };
    ShellLoadPoint {
        triangle,
        barycentric,
        offset_m: [0.0; 3],
    }
}

#[test]
fn bare_shell_matches_analytic_full_cross_point_blocks_and_redundant_nullspace() {
    let model = MovingShellModel::new(plate(), 2, || false).unwrap();
    let zero = vec![0.0; model.assembly().mass.dimension()];
    let mut motion = model.start(&zero, &zero, &[], config(), || false).unwrap();
    let mut second = point(1, [0.4, 0.2, 0.4]);
    second.offset_m = [-0.03, 0.04, -0.01];
    let points = [point(0, [0.2, 0.3, 0.5]), second];
    let response = motion.point_inverse_inertia(&points, || false).unwrap();
    assert_eq!(response.points(), points);
    assert_eq!(response.dimension(), 12);
    // Each half-square has mass 0.25 kg, split equally among its vertices.
    // The original uniform-thickness rotary inertia is m*h^2/12.
    let mass = [1.0 / 12.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 12.0];
    let weights = [[0.2, 0.3, 0.5, 0.0], [0.0, 0.4, 0.4, 0.2]];
    for column in 0..12 {
        let source = column / 6;
        let mut basis = [0.0; 6];
        basis[column % 6] = 1.0;
        let force = [basis[0], basis[1], basis[2]];
        let arm = cross(points[source].offset_m, force);
        let couple = std::array::from_fn(|axis| basis[axis + 3] + arm[axis]);
        for target in 0..2 {
            let inverse_mass: f64 = (0..4)
                .map(|n| weights[source][n] * weights[target][n] / mass[n])
                .sum();
            let inverse_rotary = inverse_mass * 12.0 / 0.1_f64.powi(2);
            let omega = couple.map(|value| inverse_rotary * value);
            let offset_velocity = cross(omega, points[target].offset_m);
            for axis in 0..3 {
                let expected = inverse_mass * force[axis] + offset_velocity[axis];
                near(
                    response.values()[(6 * target + axis) * 12 + column],
                    expected,
                    4e-12 * expected.abs().max(1.0),
                );
                near(
                    response.values()[(6 * target + 3 + axis) * 12 + column],
                    omega[axis],
                    4e-12 * omega[axis].abs().max(1.0),
                );
            }
        }
    }
    assert!(
        response.values()[6].abs() > 1e-4,
        "shared nodes must couple distinct material points"
    );
    let repeated = motion
        .point_inverse_inertia(&[points[0]; 2], || false)
        .unwrap();
    let opposite = [
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    assert_eq!(repeated.apply(&opposite).unwrap(), [[0.0; 6]; 2]);
}

#[test]
fn anisotropic_payload_response_drives_simultaneous_events_and_physical_work() {
    let body = payload();
    let model = MovingShellModel::new_with_payloads(plate(), 5, &[body], || false).unwrap();
    let zero = vec![0.0; model.assembly().mass.dimension()];
    let mut motion = model.start(&zero, &zero, &[], config(), || false).unwrap();
    let points = [vertex(0), vertex(1), vertex(2), vertex(3), body.point];
    let response = motion.point_inverse_inertia(&points, || false).unwrap();
    let dimension = response.dimension();
    // A finite payload links distinct vertices even though bare shell inertia
    // is nodally lumped. All reciprocal translation/rotation blocks survive.
    assert!(response.values()[6 * dimension].abs() > 1e-4);
    for row in 0..dimension {
        for col in 0..dimension {
            near(
                response.values()[row * dimension + col],
                response.values()[col * dimension + row],
                1e-9,
            );
        }
    }
    let mass = [1.0 / 12.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 12.0];
    let translation = [0.3, -0.1, 0.2];
    let mut rigid = [[0.0; 6]; 5];
    for node in 0..4 {
        for axis in 0..3 {
            rigid[node][axis] = mass[node] * translation[axis];
        }
    }
    for axis in 0..3 {
        rigid[4][axis] = body.mass_kg * translation[axis];
    }
    for value in response.apply(&rigid).unwrap() {
        for axis in 0..3 {
            near(value[axis], translation[axis], 2e-14);
            near(value[axis + 3], 0.0, 3e-13);
        }
    }

    let impulses = [
        [0.007, -0.003, 0.01, 0.0002, 0.0, -0.0001],
        [-0.002, 0.001, 0.003, 0.0, 0.0001, 0.0],
        [0.0; 6],
        [0.0; 6],
        [0.002, 0.001, -0.004, -0.0001, 0.0002, 0.00015],
    ];
    let predicted = response.apply(&impulses).unwrap();
    let kicks: Vec<_> = points
        .iter()
        .zip(impulses)
        .map(|(&point, j)| ShellPointImpulse {
            point,
            linear_impulse_ns: [j[0], j[1], j[2]],
            angular_impulse_nms: [j[3], j[4], j[5]],
        })
        .collect();
    let event = motion.apply_impulses(&kicks, || false).unwrap();
    for (&point, expected) in points.iter().zip(&predicted) {
        let actual = motion.observe(point).unwrap();
        for axis in 0..3 {
            near(actual.velocity_m_s[axis], expected[axis], 3e-15);
            near(
                actual.angular_velocity_rad_s[axis],
                expected[axis + 3],
                3e-13,
            );
        }
    }
    // Independently inventory bare nodal inertia plus the physical COM tensor.
    let v = &motion.state().velocity;
    let mut energy = 0.0;
    for node in 0..4 {
        let linear = std::array::from_fn(|axis| v[6 * node + axis]);
        let spin = std::array::from_fn(|axis| v[6 * node + 3 + axis]);
        energy +=
            0.5 * mass[node] * (dot(linear, linear) + 0.1_f64.powi(2) / 12.0 * dot(spin, spin));
    }
    let (mut com, mut omega) = ([0.0; 3], [0.0; 3]);
    for (node, weight) in [0, 1, 2].into_iter().zip(body.point.barycentric) {
        let spin = std::array::from_fn(|axis| v[6 * node + 3 + axis]);
        let offset = cross(spin, body.point.offset_m);
        for axis in 0..3 {
            com[axis] += weight * (v[6 * node + axis] + offset[axis]);
            omega[axis] += weight * spin[axis];
        }
    }
    energy += 0.5
        * (body.mass_kg * dot(com, com)
            + dot(omega, body.inertia_com_kg_m2.map(|row| dot(row, omega))));
    let work: f64 = 0.5
        * impulses
            .iter()
            .flatten()
            .zip(predicted.iter().flatten())
            .map(|(j, dv)| j * dv)
            .sum::<f64>();
    near(work, energy, 2e-16);
    near(event.kinetic_energy_change_j, energy, 2e-16);
    assert_eq!(
        motion.point_inverse_inertia(&points, || false).unwrap(),
        response
    );
}

#[test]
fn driven_constraints_budgets_and_cancellation_leave_motion_and_null_directions_intact() {
    let mut p = plate();
    p.nodes.truncate(3);
    p.triangles.truncate(1);
    p.support = Some(ShellSupport {
        node_indices: [2, 0, 1],
        normal: [0.0, 0.0, 1.0],
    });
    p.budget.max_matrix_entries = 18 * 18;
    let body = payload();
    let model = MovingShellModel::new_with_payloads(p, 4, &[body], || false).unwrap();
    let zero = vec![0.0; model.assembly().mass.dimension()];
    let support = ShellSupportMotion {
        velocity_m_s: [[0.05, -0.02, 0.01]; 3],
        ..Default::default()
    };
    let loads = [wrench(body.point)];
    let start = || {
        model
            .start_driven(&zero, &zero, support, &loads, config(), || false)
            .unwrap()
    };
    let mut motion = start();
    let before = motion.state().clone();
    let reactions = motion.support_reactions().to_vec();
    let points = [vertex(0), body.point];
    let mut polls = 0;
    let response = motion
        .point_inverse_inertia(&points, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut fixed = model
        .start(&zero, &zero, &loads, config(), || false)
        .unwrap();
    assert_eq!(
        fixed.point_inverse_inertia(&points, || false).unwrap(),
        response
    );
    for row in 0..12 {
        for axis in 0..3 {
            assert_eq!(response.values()[row * 12 + axis], 0.0);
            assert_eq!(response.values()[axis * 12 + row], 0.0);
        }
    }
    assert!(
        response.values()[6 * 12 + 6] > 0.0,
        "offset material points still move through the free rotations"
    );
    let mut invalid = points;
    invalid[1].triangle = usize::MAX;
    assert!(motion.point_inverse_inertia(&invalid, || false).is_err());
    // Matrix-entry budget refuses before a quadratic allocation, independently
    // of the still-admitted four-point input budget.
    assert!(
        motion
            .point_inverse_inertia(&[points[0]; 4], || false)
            .is_err()
    );
    assert!(
        motion
            .point_inverse_inertia(&[points[0]; 5], || false)
            .is_err()
    );
    let mut calls = 0;
    assert!(matches!(
        motion.point_inverse_inertia(&points, || {
            calls += 1;
            calls == polls
        }),
        Err(MovingLoadError::Cancelled)
    ));
    assert_eq!(motion.state(), &before);
    assert_eq!(motion.support_reactions(), reactions);
    assert_eq!(motion.support_motion(), Some(&support));
    assert!(response.apply(&[[0.0; 6]]).is_err());
    let mut invalid = [[0.0; 6]; 2];
    invalid[1][5] = f64::NAN;
    assert!(response.apply(&invalid).is_err());
    let empty = motion.point_inverse_inertia(&[], || false).unwrap();
    assert_eq!(empty.dimension(), 0);
    assert!(empty.apply(&[]).unwrap().is_empty());
    let kick = ShellPointImpulse {
        point: body.point,
        linear_impulse_ns: [0.007, -0.003, 0.01],
        angular_impulse_nms: [0.0002, -0.0001, 0.00015],
    };
    let mut control = start();
    assert_eq!(
        motion.apply_impulses(&[kick], || false).unwrap(),
        control.apply_impulses(&[kick], || false).unwrap()
    );
    assert_eq!(motion.state(), control.state());
    assert_eq!(motion.support_reactions(), control.support_reactions());
}
