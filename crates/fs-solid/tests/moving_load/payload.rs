//! G0/G1/G3/G4: finite payload inertia on the original shell trajectory.
use super::*;

fn payload() -> ShellRigidPayload {
    ShellRigidPayload {
        point: point(0, [0.2, 0.3, 0.5]),
        mass_kg: 2.3,
        inertia_com_kg_m2: [
            [0.006, 0.001, -0.0004],
            [0.001, 0.008, 0.0007],
            [-0.0004, 0.0007, 0.010],
        ],
    }
}

#[test]
fn finite_payload_matches_com_energy_and_reference_angular_momentum() {
    let mut p = plate();
    p.damping = DampingModel::Rayleigh {
        mass_proportional_per_s: 0.3,
        stiffness_proportional_s: 0.0001,
    };
    let bare = MovingShellModel::new(p.clone(), 2, || false).unwrap();
    let body = payload();
    let loaded = MovingShellModel::new_with_payloads(p, 2, &[body], || false).unwrap();
    assert_eq!(loaded.payloads(), &[body]);
    assert_eq!(loaded.assembly().full_stiffness, bare.assembly().full_stiffness);
    assert_eq!(loaded.assembly().full_damping, bare.assembly().full_damping);
    let n = bare.assembly().mass.dimension();
    let velocity: Vec<_> = (0..n).map(|i| 0.3 * (0.7 * i as f64).cos()).collect();
    let (mut com_velocity, mut omega, mut position) = ([0.; 3], [0.; 3], body.point.offset_m);
    // Independent material velocity at the COM, with no new inertia helper.
    for (node, weight) in bare.plate().triangles[body.point.triangle]
        .into_iter()
        .zip(body.point.barycentric)
    {
        let rotation = std::array::from_fn(|i| velocity[6 * node + 3 + i]);
        let spin = cross(rotation, body.point.offset_m);
        for axis in 0..3 {
            com_velocity[axis] += weight * (velocity[6 * node + axis] + spin[axis]);
            omega[axis] += weight * rotation[axis];
            position[axis] += weight * bare.plate().nodes[node].position_m[axis];
        }
    }
    let spin_momentum = body.inertia_com_kg_m2.map(|row| dot(row, omega));
    let linear = com_velocity.map(|v| body.mass_kg * v);
    let orbital = cross(position, linear);
    let expected_energy = 0.5 * (dot(com_velocity, linear) + dot(omega, spin_momentum));
    near(
        loaded.assembly().mass.quadratic_energy(&velocity)
            - bare.assembly().mass.quadratic_energy(&velocity),
        expected_energy,
        3e-15,
    );
    let actual: Vec<_> = loaded.assembly().mass.apply(&velocity)
        .into_iter()
        .zip(bare.assembly().mass.apply(&velocity))
        .map(|(a, b)| a - b)
        .collect();
    let (mut actual_linear, mut actual_angular) = ([0.; 3], [0.; 3]);
    for (node, momentum) in bare.plate().nodes.iter().zip(actual.chunks_exact(6)) {
        let orbit = cross(node.position_m, [momentum[0], momentum[1], momentum[2]]);
        for axis in 0..3 {
            actual_linear[axis] += momentum[axis];
            actual_angular[axis] += orbit[axis] + momentum[3 + axis];
        }
    }
    for axis in 0..3 {
        near(actual_linear[axis], linear[axis], 3e-15);
        near(actual_angular[axis], orbital[axis] + spin_momentum[axis], 6e-15);
    }
    // P1-consistent inertia couples different nodes; a lumped mass fails this.
    assert!(loaded.assembly().mass.values()[6] > 0.0);
    for i in 0..n {
        for j in 0..n {
            assert_eq!(
                loaded.assembly().mass.values()[i * n + j],
                loaded.assembly().mass.values()[j * n + i],
            );
        }
    }
}

#[test]
fn payload_force_drives_the_correct_native_acceleration_and_impulse() {
    let bare = MovingShellModel::new(plate(), 5, || false).unwrap();
    let body = payload();
    let model = MovingShellModel::new_with_payloads(plate(), 5, &[body], || false).unwrap();
    let n = model.assembly().mass.dimension();
    let acceleration = [0.3, -0.1, 0.2];
    let mut loads: Vec<_> = (0..4).map(|node| {
        let (triangle, barycentric) = match node {
            0 => (0, [1., 0., 0.]),
            1 => (0, [0., 1., 0.]),
            2 => (0, [0., 0., 1.]),
            _ => (1, [0., 1., 0.]),
        };
        let nodal_mass = bare.assembly().mass.values()[(6 * node) * n + 6 * node];
        ShellPointWrench {
            point: ShellLoadPoint { triangle, barycentric, offset_m: [0.; 3] },
            force_n: acceleration.map(|a| nodal_mass * a),
            couple_nm: [0.; 3],
        }
    }).collect();
    loads.push(ShellPointWrench {
        point: body.point,
        force_n: acceleration.map(|a| body.mass_kg * a),
        couple_nm: [0.; 3],
    });
    let zero = vec![0.; n];
    let mut motion = model.start(&zero, &zero, &loads, config(), || false).unwrap();
    for _ in 0..40 {
        let step = motion.advance(motion.state().time_s + config().dt_s, &loads, || false).unwrap();
        let t = motion.state().time_s;
        for node in 0..4 {
            for axis in 0..3 {
                near(motion.state().velocity[6 * node + axis], acceleration[axis] * t, 1e-14);
                near(
                    motion.state().displacement[6 * node + axis],
                    0.5 * acceleration[axis] * t * t,
                    1e-17,
                );
                near(motion.state().velocity[6 * node + 3 + axis], 0., 1e-14);
            }
        }
        // The authored unit square weighs 0.5 kg; the body contributes 2.3 kg.
        for axis in 0..3 {
            near(
                step.momentum_change.linear[axis],
                (0.5 + body.mass_kg) * acceleration[axis] * config().dt_s,
                1e-14,
            );
            near(step.momentum_residual.linear[axis], 0., 1e-14);
            near(step.momentum_residual.angular[axis], 0., 2e-14);
        }
        near(step.energy_residual_j, 0., 2e-15);
    }
}

#[test]
fn driven_pins_accelerate_payload_mass_and_supply_its_kinetic_energy() {
    let mut p = plate();
    p.nodes.truncate(3);
    p.triangles.truncate(1);
    p.support = Some(ShellSupport { node_indices: [2, 0, 1], normal: [0., 0., 1.] });
    let mut body = payload();
    body.point.offset_m = [0.; 3];
    body.inertia_com_kg_m2 = [[0.; 3]; 3];
    let model = MovingShellModel::new_with_payloads(p, 1, &[body], || false).unwrap();
    let zero = vec![0.; model.assembly().mass.dimension()];
    let acceleration = [0.3, -0.1, 0.2];
    let supports = ShellSupportMotion {
        acceleration_m_s2: [acceleration; 3],
        ..ShellSupportMotion::default()
    };
    let mut motion = model.start_driven(&zero, &zero, supports, &[], config(), || false).unwrap();
    let report = motion.advance_driven(config().dt_s, &[], [acceleration; 3], || false).unwrap();
    // The authored half-square shell weighs 0.25 kg.
    let total_mass = 0.25 + body.mass_kg;
    for axis in 0..3 {
        near(
            report.support_impulse.linear[axis],
            total_mass * acceleration[axis] * config().dt_s,
            1e-15,
        );
        near(report.momentum_residual.linear[axis], 0., 1e-15);
        near(report.momentum_residual.angular[axis], 0., 2e-15);
    }
    let kinetic = 0.5 * total_mass * dot(acceleration, acceleration) * config().dt_s.powi(2);
    near(motion.state().kinetic_energy_j, kinetic, 1e-20);
    near(report.support_work_j, kinetic, 1e-20);
    near(report.energy_residual_j, 0., 1e-20);
    for velocity in &motion.state().velocity {
        near(*velocity, 0., 1e-15);
    }
}

#[test]
fn payload_admission_empty_parity_and_atomic_refusal() {
    let bare = MovingShellModel::new(plate(), 1, || false).unwrap();
    let empty = MovingShellModel::new_with_payloads(plate(), 1, &[], || false).unwrap();
    assert_eq!(bare.assembly(), empty.assembly());
    assert!(empty.payloads().is_empty());
    let zero = vec![0.; bare.assembly().mass.dimension()];
    let load = wrench(point(0, [0.2, 0.3, 0.5]));
    let mut a = bare.start(&zero, &zero, &[load], config(), || false).unwrap();
    let mut b = empty.start(&zero, &zero, &[load], config(), || false).unwrap();
    assert_eq!(
        a.advance(config().dt_s, &[load], || false).unwrap(),
        b.advance(config().dt_s, &[load], || false).unwrap(),
    );
    assert_eq!(a.state(), b.state());
    let mut bad = [payload(); 5];
    bad[0].mass_kg = -1.;
    bad[1].inertia_com_kg_m2[0][1] += 0.001;
    bad[2].inertia_com_kg_m2 = [[1., 0., 0.], [0., 1., 0.], [0., 0., 3.]];
    bad[3].inertia_com_kg_m2 = [[1., 2., 0.], [2., 1., 0.], [0., 0., 1.]];
    bad[4].point.barycentric = [0.2, 0.3, 0.6];
    for body in bad {
        assert!(MovingShellModel::new_with_payloads(plate(), 1, &[body], || false).is_err());
    }
    assert!(MovingShellModel::new_with_payloads(plate(), 1, &[payload(); 2], || false).is_err());
    let mut calls = 0;
    assert!(matches!(
        MovingShellModel::new_with_payloads(plate(), 1, &[payload()], || {
            calls += 1;
            calls == 7
        }),
        Err(MovingLoadError::Cancelled)
    ));
    // Cancellation with non-diagonal payload inertia still retains old state.
    let model = MovingShellModel::new_with_payloads(plate(), 1, &[payload()], || false).unwrap();
    let mut motion = model.start(&zero, &zero, &[load], config(), || false).unwrap();
    let accepted = motion.state().clone();
    calls = 0;
    assert!(matches!(
        motion.advance(config().dt_s, &[load], || {
            calls += 1;
            calls == 6
        }),
        Err(MovingLoadError::Cancelled)
    ));
    assert_eq!(motion.state(), &accepted);
    let step = motion.advance(config().dt_s, &[load], || false).unwrap();
    let mut fresh = model.start(&zero, &zero, &[load], config(), || false).unwrap();
    assert_eq!(step, fresh.advance(config().dt_s, &[load], || false).unwrap());
    assert_eq!(motion.state(), fresh.state());
}
