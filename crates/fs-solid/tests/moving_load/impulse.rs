//! Native event impulses, checked against shell plus rigid-COM inventories.
use super::*;

fn body() -> ShellRigidPayload {
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

fn event(point: ShellLoadPoint) -> ShellPointImpulse {
    ShellPointImpulse {
        point,
        linear_impulse_ns: [0.007, -0.003, 0.01],
        angular_impulse_nms: [0.0002, -0.0001, 0.00015],
    }
}

fn full(model: &MovingShellModel, free: &[f64], pins: Option<&[[f64; 3]; 3]>) -> Vec<f64> {
    let mut result = vec![0.0; 6 * model.plate().nodes.len()];
    for (&index, &value) in model.assembly().free_dofs.iter().zip(free) {
        result[index] = value;
    }
    if let Some(pins) = pins {
        for (node, values) in model
            .plate()
            .support
            .unwrap()
            .node_indices
            .into_iter()
            .zip(pins)
        {
            result[6 * node..6 * node + 3].copy_from_slice(values);
        }
    }
    result
}

// Independent physical material-point interpolation; never use event receipts
// or the new payload matrix to construct the expected momentum and work.
fn material(
    p: &ShellPlate,
    point: ShellLoadPoint,
    velocity: &[f64],
) -> ([f64; 3], [f64; 3], [f64; 3]) {
    let (mut position, mut v, mut omega) = (point.offset_m, [0.0; 3], [0.0; 3]);
    for (node, weight) in p.triangles[point.triangle]
        .into_iter()
        .zip(point.barycentric)
    {
        let spin = std::array::from_fn(|axis| velocity[6 * node + 3 + axis]);
        let offset_velocity = cross(spin, point.offset_m);
        for axis in 0..3 {
            position[axis] += weight * p.nodes[node].position_m[axis];
            v[axis] += weight * (velocity[6 * node + axis] + offset_velocity[axis]);
            omega[axis] += weight * spin[axis];
        }
    }
    (position, v, omega)
}

fn inventory(
    bare: &MovingShellModel,
    body: ShellRigidPayload,
    velocity: &[f64],
) -> (ShellSpatialResultant, f64) {
    let mut result = ShellSpatialResultant::default();
    let shell_momentum = bare.assembly().full_mass.apply(velocity);
    for (node, momentum) in bare
        .plate()
        .nodes
        .iter()
        .zip(shell_momentum.chunks_exact(6))
    {
        let orbit = cross(node.position_m, [momentum[0], momentum[1], momentum[2]]);
        for axis in 0..3 {
            result.linear[axis] += momentum[axis];
            result.angular[axis] += orbit[axis] + momentum[axis + 3];
        }
    }
    let (position, v, omega) = material(bare.plate(), body.point, velocity);
    let linear = v.map(|value| body.mass_kg * value);
    let spin = body.inertia_com_kg_m2.map(|row| dot(row, omega));
    let orbit = cross(position, linear);
    for axis in 0..3 {
        result.linear[axis] += linear[axis];
        result.angular[axis] += orbit[axis] + spin[axis];
    }
    let energy = bare.assembly().full_mass.quadratic_energy(velocity)
        + 0.5 * (dot(linear, v) + dot(spin, omega));
    (result, energy)
}

fn held_force_balance(
    model: &MovingShellModel,
    motion: &MovingShellMotion<'_>,
    loads: &[ShellPointWrench],
) {
    let s = motion.support_motion();
    let q = full(
        model,
        &motion.state().displacement,
        s.map(|s| &s.displacement_m),
    );
    let v = full(model, &motion.state().velocity, s.map(|s| &s.velocity_m_s));
    let a = full(
        model,
        &motion.state().acceleration,
        s.map(|s| &s.acceleration_m_s2),
    );
    let inertia = model.assembly().full_mass.apply(&a);
    let elastic = model.assembly().full_stiffness.apply(&q);
    let viscous = model.assembly().full_damping.as_ref().map(|c| c.apply(&v));
    let force = model.project(loads, || false).unwrap().full;
    for i in 0..q.len() {
        near(
            inertia[i] + elastic[i] + viscous.as_ref().map_or(0.0, |c| c[i]) - force[i],
            motion.support_reactions()[i],
            3e-12,
        );
    }
}

#[test]
fn off_center_payload_event_preserves_native_momentum_work_and_force_continuation() {
    let mut p = plate();
    p.damping = DampingModel::Rayleigh {
        mass_proportional_per_s: 0.4,
        stiffness_proportional_s: 0.0002,
    };
    let body = body();
    let bare = MovingShellModel::new(p.clone(), 2, || false).unwrap();
    let model = MovingShellModel::new_with_payloads(p, 2, &[body], || false).unwrap();
    let n = model.assembly().mass.dimension();
    let q: Vec<_> = (0..n).map(|i| 1e-4 * (0.6 * i as f64).sin()).collect();
    let v: Vec<_> = (0..n).map(|i| 0.03 * (0.7 * i as f64).cos()).collect();
    let loads = [wrench(point(1, [0.3, 0.2, 0.5]))];
    let mut motion = model.start(&q, &v, &loads, config(), || false).unwrap();
    let before = motion.state().clone();
    let kick = event(point(0, [0.4, 0.35, 0.25]));
    let receipt = motion.apply_impulses(&[kick], || false).unwrap();
    assert_eq!(motion.state().time_s, before.time_s);
    assert_eq!(motion.state().step, before.step);
    assert_eq!(motion.state().displacement, before.displacement);
    assert_eq!(motion.state().strain_energy_j, before.strain_energy_j);
    let (old_momentum, old_energy) = inventory(&bare, body, &v);
    let (momentum, energy) = inventory(&bare, body, &motion.state().velocity);
    let (position, old_v, old_omega) = material(model.plate(), kick.point, &v);
    let (_, new_v, new_omega) = material(model.plate(), kick.point, &motion.state().velocity);
    let moment = cross(position, kick.linear_impulse_ns);
    for axis in 0..3 {
        near(
            momentum.linear[axis] - old_momentum.linear[axis],
            kick.linear_impulse_ns[axis],
            2e-15,
        );
        near(
            momentum.angular[axis] - old_momentum.angular[axis],
            moment[axis] + kick.angular_impulse_nms[axis],
            3e-15,
        );
        near(
            receipt.momentum_change.linear[axis],
            momentum.linear[axis] - old_momentum.linear[axis],
            2e-15,
        );
        near(
            receipt.momentum_change.angular[axis],
            momentum.angular[axis] - old_momentum.angular[axis],
            3e-15,
        );
        near(receipt.momentum_residual.linear[axis], 0.0, 2e-15);
        near(receipt.momentum_residual.angular[axis], 0.0, 3e-15);
    }
    let mean_v = std::array::from_fn(|i| 0.5 * (old_v[i] + new_v[i]));
    let mean_omega = std::array::from_fn(|i| 0.5 * (old_omega[i] + new_omega[i]));
    let work = dot(kick.linear_impulse_ns, mean_v) + dot(kick.angular_impulse_nms, mean_omega);
    near(receipt.applied_work_j, work, 2e-16);
    near(receipt.kinetic_energy_change_j, energy - old_energy, 2e-16);
    near(receipt.energy_residual_j, 0.0, 2e-16);
    near(motion.state().kinetic_energy_j, energy, 2e-16);
    assert_eq!(receipt.support_impulse, ShellSpatialResultant::default());
    assert_eq!(receipt.support_work_j, 0.0);
    held_force_balance(&model, &motion, &loads);
    let mut restarted = model
        .start(
            &motion.state().displacement,
            &motion.state().velocity,
            &loads,
            config(),
            || false,
        )
        .unwrap();
    assert_eq!(motion.state(), restarted.state());
    assert_eq!(
        motion.advance(config().dt_s, &loads, || false).unwrap(),
        restarted.advance(config().dt_s, &loads, || false).unwrap()
    );
    assert_eq!(motion.state(), restarted.state());
    let accepted = motion.state().clone();
    let again = motion.apply_impulses(&[kick], || false).unwrap();
    assert_eq!(again.time_s, accepted.time_s);
    assert_eq!(motion.state().time_s, accepted.time_s);
    assert_eq!(motion.state().step, accepted.step);
    assert_eq!(motion.state().displacement, accepted.displacement);
    held_force_balance(&model, &motion, &loads);
}

#[test]
fn prescribed_velocity_supplies_payload_event_reaction_and_its_actual_work() {
    for velocity in [[0.0; 3], [0.05, -0.02, 0.01]] {
        let mut p = plate();
        p.nodes.truncate(3);
        p.triangles.truncate(1);
        p.support = Some(ShellSupport {
            node_indices: [2, 0, 1],
            normal: [0.0, 0.0, 1.0],
        });
        p.damping = DampingModel::Rayleigh {
            mass_proportional_per_s: 0.4,
            stiffness_proportional_s: 0.0002,
        };
        let body = body();
        let bare = MovingShellModel::new(p.clone(), 1, || false).unwrap();
        let model = MovingShellModel::new_with_payloads(p, 1, &[body], || false).unwrap();
        let support = ShellSupportMotion {
            velocity_m_s: [velocity; 3],
            ..Default::default()
        };
        let zero = vec![0.0; model.assembly().mass.dimension()];
        let rotations: Vec<_> = (0..zero.len())
            .map(|i| 0.02 * (0.7 * i as f64).cos())
            .collect();
        let loads = [wrench(body.point)];
        let mut motion = model
            .start_driven(&zero, &rotations, support, &loads, config(), || false)
            .unwrap();
        let old_state = motion.state().clone();
        let old_full = full(&model, &rotations, Some(&support.velocity_m_s));
        let kick = event(body.point);
        let receipt = motion.apply_impulses(&[kick], || false).unwrap();
        let new_full = full(
            &model,
            &motion.state().velocity,
            Some(&support.velocity_m_s),
        );
        let (old_momentum, old_energy) = inventory(&bare, body, &old_full);
        let (momentum, energy) = inventory(&bare, body, &new_full);
        let (position, old_v, old_omega) = material(model.plate(), body.point, &old_full);
        let (_, new_v, new_omega) = material(model.plate(), body.point, &new_full);
        let reaction = std::array::from_fn(|i| {
            body.mass_kg * (new_v[i] - old_v[i]) - kick.linear_impulse_ns[i]
        });
        let midsurface = std::array::from_fn(|i| position[i] - body.point.offset_m[i]);
        let reaction_moment = cross(midsurface, reaction);
        for axis in 0..3 {
            near(receipt.support_impulse.linear[axis], reaction[axis], 2e-15);
            near(
                receipt.support_impulse.angular[axis],
                reaction_moment[axis],
                3e-15,
            );
            near(
                receipt.momentum_change.linear[axis],
                momentum.linear[axis] - old_momentum.linear[axis],
                2e-15,
            );
            near(
                receipt.momentum_change.angular[axis],
                momentum.angular[axis] - old_momentum.angular[axis],
                3e-15,
            );
            near(receipt.momentum_residual.linear[axis], 0.0, 2e-15);
            near(receipt.momentum_residual.angular[axis], 0.0, 3e-15);
        }
        let work = dot(
            kick.linear_impulse_ns,
            std::array::from_fn(|i| 0.5 * (old_v[i] + new_v[i])),
        ) + dot(
            kick.angular_impulse_nms,
            std::array::from_fn(|i| 0.5 * (old_omega[i] + new_omega[i])),
        );
        near(receipt.applied_work_j, work, 2e-16);
        near(receipt.support_work_j, dot(reaction, velocity), 2e-16);
        near(receipt.kinetic_energy_change_j, energy - old_energy, 2e-16);
        near(receipt.energy_residual_j, 0.0, 2e-16);
        assert_eq!(motion.support_motion(), Some(&support));
        assert_eq!(motion.state().displacement, old_state.displacement);
        assert_eq!(motion.state().time_s, old_state.time_s);
        assert_eq!(motion.state().step, old_state.step);
        held_force_balance(&model, &motion, &loads);

        // Zero free impulse: the pins receive all applied impulse and exchange
        // exactly opposite work even when they have a nonzero prescribed speed.
        let supported = ShellPointImpulse {
            point: ShellLoadPoint {
                triangle: 0,
                barycentric: [1.0, 0.0, 0.0],
                offset_m: [0.0; 3],
            },
            linear_impulse_ns: [0.01, -0.02, 0.03],
            angular_impulse_nms: [0.0; 3],
        };
        let before = motion.state().clone();
        let report = motion.apply_impulses(&[supported], || false).unwrap();
        assert_eq!(motion.state(), &before);
        for axis in 0..3 {
            near(
                report.support_impulse.linear[axis],
                -supported.linear_impulse_ns[axis],
                1e-17,
            );
        }
        assert_eq!(report.kinetic_energy_change_j, 0.0);
        near(
            report.support_work_j,
            -dot(supported.linear_impulse_ns, velocity),
            1e-17,
        );
        near(report.applied_work_j + report.support_work_j, 0.0, 1e-17);
    }
}

#[test]
fn event_refusal_late_cancel_and_unrepresentable_kicks_preserve_accepted_motion() {
    let mut p = plate();
    p.support = Some(ShellSupport {
        node_indices: [2, 0, 1],
        normal: [0.0, 0.0, 1.0],
    });
    let model = MovingShellModel::new_with_payloads(p, 1, &[body()], || false).unwrap();
    let zero = vec![0.0; model.assembly().mass.dimension()];
    let support = ShellSupportMotion {
        velocity_m_s: [[0.05, -0.02, 0.01]; 3],
        ..Default::default()
    };
    let loads = [wrench(body().point)];
    let start = || {
        model
            .start_driven(&zero, &zero, support, &loads, config(), || false)
            .unwrap()
    };
    let kick = event(body().point);
    let mut control = start();
    let mut polls = 0;
    let expected = control
        .apply_impulses(&[kick], || {
            polls += 1;
            false
        })
        .unwrap();
    let mut motion = start();
    let before = motion.state().clone();
    let reactions = motion.support_reactions().to_vec();
    let mut invalid = kick;
    invalid.linear_impulse_ns[1] = f64::NAN;
    assert!(motion.apply_impulses(&[invalid], || false).is_err());
    assert!(motion.apply_impulses(&[kick; 2], || false).is_err());
    let mut calls = 0;
    assert!(matches!(
        motion.apply_impulses(&[kick], || {
            calls += 1;
            calls == polls
        }),
        Err(MovingLoadError::Cancelled)
    ));
    assert_eq!(motion.state(), &before);
    assert_eq!(motion.support_reactions(), reactions);
    assert_eq!(motion.support_motion(), Some(&support));
    assert_eq!(motion.apply_impulses(&[kick], || false).unwrap(), expected);
    assert_eq!(motion.state(), control.state());
    assert_eq!(motion.support_reactions(), control.support_reactions());
    assert_eq!(
        motion
            .advance_driven(config().dt_s, &loads, support.acceleration_m_s2, || false)
            .unwrap(),
        control
            .advance_driven(config().dt_s, &loads, support.acceleration_m_s2, || false)
            .unwrap()
    );
    assert_eq!(motion.state(), control.state());

    let free = MovingShellModel::new(plate(), 1, || false).unwrap();
    let zero = vec![0.0; free.assembly().mass.dimension()];
    let mut huge = zero.clone();
    for node in 0..4 {
        huge[6 * node] = 1e12;
    }
    let mut motion = free.start(&zero, &huge, &[], config(), || false).unwrap();
    let before = motion.state().clone();
    let tiny = ShellPointImpulse {
        point: ShellLoadPoint {
            triangle: 0,
            barycentric: [1.0, 0.0, 0.0],
            offset_m: [0.0; 3],
        },
        linear_impulse_ns: [1e-200, 0.0, 0.0],
        angular_impulse_nms: [0.0; 3],
    };
    assert!(matches!(
        motion.apply_impulses(&[tiny], || false),
        Err(MovingLoadError::Invalid(_))
    ));
    assert_eq!(motion.state(), &before);
    let empty = motion.apply_impulses(&[], || false).unwrap();
    assert_eq!(motion.state(), &before);
    assert_eq!(empty.momentum_change, ShellSpatialResultant::default());
    assert_eq!(empty.kinetic_energy_change_j, 0.0);

    // A representable kick still rounds on addition. Its receipt must carry
    // the actual accepted inventory, not the ideal impulse passed to the LU.
    let mut background = zero.clone();
    background[0] = 1.0;
    let mut motion = free
        .start(&zero, &background, &[], config(), || false)
        .unwrap();
    let kick = ShellPointImpulse {
        linear_impulse_ns: [1e-4, 0.0, 0.0],
        ..tiny
    };
    let receipt = motion.apply_impulses(&[kick], || false).unwrap();
    let actual = free.assembly().mass.values()[0] * (motion.state().velocity[0] - 1.0);
    assert_ne!(actual, kick.linear_impulse_ns[0]);
    assert_eq!(receipt.momentum_change.linear[0], actual);
    assert_eq!(
        receipt.kinetic_energy_change_j,
        f64::midpoint(1.0, motion.state().velocity[0]) * actual
    );
}
