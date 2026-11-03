use super::*;

fn supported(damped: bool) -> MovingShellModel {
    let mut p = plate();
    // Deliberately neither node nor reduced-DOF order.
    p.support = Some(ShellSupport {
        node_indices: [2, 0, 1],
        normal: [0.0, 0.0, 1.0],
    });
    if damped {
        p.damping = DampingModel::Rayleigh {
            mass_proportional_per_s: 0.4,
            stiffness_proportional_s: 0.002,
        };
    }
    MovingShellModel::new(p, 4, || false).unwrap()
}

fn node_point(node: usize) -> ShellLoadPoint {
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

fn restrict(model: &MovingShellModel, full: &[f64]) -> Vec<f64> {
    model
        .assembly()
        .free_dofs
        .iter()
        .map(|&i| full[i])
        .collect()
}

fn prescribe(model: &MovingShellModel, q: &[f64], v: &[f64], a: &[f64]) -> ShellSupportMotion {
    let nodes = model.plate().support.unwrap().node_indices;
    let sample =
        |field: &[f64]| std::array::from_fn(|i| std::array::from_fn(|j| field[6 * nodes[i] + j]));
    ShellSupportMotion {
        displacement_m: sample(q),
        velocity_m_s: sample(v),
        acceleration_m_s2: sample(a),
    }
}

fn check_balance(step: &MovingShellStep, tolerance: f64) {
    for value in step
        .momentum_residual
        .linear
        .into_iter()
        .chain(step.momentum_residual.angular)
    {
        near(value, 0.0, tolerance);
    }
    near(step.energy_residual_j, 0.0, tolerance);
}

#[test]
fn driven_pins_reproduce_analytic_rigid_acceleration_and_their_actual_work() {
    let model = supported(false);
    let a = [0.3, -0.2, 0.6];
    let total_mass = 0.5;
    let free_mass = 1.0 / 12.0;
    let pin_mass = 5.0 / 12.0;
    let centroid = [1.5, -1.5, 0.5];
    let force = ShellPointWrench {
        point: node_point(3),
        force_n: a.map(|v| free_mass * v),
        couple_nm: [0.0; 3],
    };
    let support = ShellSupportMotion {
        acceleration_m_s2: [a; 3],
        ..Default::default()
    };
    let zero = vec![0.0; model.assembly().mass.dimension()];
    let mut motion = model
        .start_driven(&zero, &zero, support, &[force], config(), || false)
        .unwrap();
    let dt = config().dt_s;
    let mut time = 0.0;
    for _ in 0..128 {
        let midpoint_time = time + 0.5 * dt;
        time += dt;
        let step = motion
            .advance_driven(time, &[force], [a; 3], || false)
            .unwrap();
        let angular = cross(centroid, a.map(|x| dt * total_mass * x));
        for axis in 0..3 {
            near(
                step.momentum_change.linear[axis],
                dt * total_mass * a[axis],
                3e-16,
            );
            near(
                step.applied_impulse.linear[axis],
                dt * free_mass * a[axis],
                3e-17,
            );
            near(
                step.support_impulse.linear[axis],
                dt * pin_mass * a[axis],
                3e-16,
            );
            near(step.momentum_change.angular[axis], angular[axis], 8e-16);
        }
        near(
            step.support_work_j,
            dt * pin_mass * dot(a, a) * midpoint_time,
            2e-17,
        );
        near(
            step.applied_work_j,
            dt * free_mass * dot(a, a) * midpoint_time,
            2e-17,
        );
        near(
            motion.state().kinetic_energy_j,
            0.5 * total_mass * dot(a, a) * time * time,
            5e-17,
        );
        for node in 0..4 {
            let point = motion.observe(node_point(node)).unwrap();
            for axis in 0..3 {
                near(point.velocity_m_s[axis], a[axis] * time, 2e-14);
                near(
                    point.displacement_m[axis],
                    0.5 * a[axis] * time * time,
                    2e-16,
                );
            }
        }
        check_balance(&step, 2e-15);
    }
}

fn manufactured(time: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let q = (0..24)
        .map(|i| {
            let p = (i % 7) as f64 - 3.0;
            2e-5 * p + time * (0.002 * p + 0.5 * time * 0.3 * p)
        })
        .collect();
    let v = (0..24)
        .map(|i| {
            let p = (i % 7) as f64 - 3.0;
            0.002 * p + time * 0.3 * p
        })
        .collect();
    let a = (0..24).map(|i| 0.3 * ((i % 7) as f64 - 3.0)).collect();
    (q, v, a)
}

fn mechanical_force(model: &MovingShellModel, q: &[f64], v: &[f64], a: &[f64]) -> Vec<f64> {
    let assembly = model.assembly();
    let inertia = assembly.full_mass.apply(a);
    let elastic = assembly.full_stiffness.apply(q);
    let viscous = assembly.full_damping.as_ref().map(|c| c.apply(v));
    (0..q.len())
        .map(|i| inertia[i] + elastic[i] + viscous.as_ref().map_or(0.0, |c| c[i]))
        .collect()
}

fn free_loads(model: &MovingShellModel, force: &[f64]) -> Vec<ShellPointWrench> {
    let free = &model.assembly().free_dofs;
    let component = |i| if free.contains(&i) { force[i] } else { 0.0 };
    (0..4)
        .map(|node| ShellPointWrench {
            point: node_point(node),
            force_n: std::array::from_fn(|axis| component(6 * node + axis)),
            couple_nm: std::array::from_fn(|axis| component(6 * node + 3 + axis)),
        })
        .collect()
}

#[test]
fn nonuniform_prescribed_motion_keeps_full_elastic_viscous_and_reaction_coupling() {
    for damped in [false, true] {
        let model = supported(damped);
        let (q, v, a) = manufactured(0.0);
        let support = prescribe(&model, &q, &v, &a);
        let force = mechanical_force(&model, &q, &v, &a);
        let initial = free_loads(&model, &force);
        let mut motion = model
            .start_driven(
                &restrict(&model, &q),
                &restrict(&model, &v),
                support,
                &initial,
                config(),
                || false,
            )
            .unwrap();
        let mut time = 0.0;
        let mut support_work = 0.0;
        for _ in 0..128 {
            time += config().dt_s;
            let (q, v, a) = manufactured(time);
            let expected_support = prescribe(&model, &q, &v, &a);
            let force = mechanical_force(&model, &q, &v, &a);
            let loads = free_loads(&model, &force);
            let step = motion
                .advance_driven(time, &loads, expected_support.acceleration_m_s2, || false)
                .unwrap();
            for (free, &full) in model.assembly().free_dofs.iter().enumerate() {
                near(motion.state().displacement[free], q[full], 2e-16);
                near(motion.state().velocity[free], v[full], 3e-14);
                near(motion.state().acceleration[free], a[full], 2e-11);
            }
            for (full, &reaction) in motion.support_reactions().iter().enumerate() {
                let expected = if model.assembly().free_dofs.contains(&full) {
                    0.0
                } else {
                    force[full]
                };
                near(reaction, expected, 2e-13);
            }
            let energy = model.assembly().full_mass.quadratic_energy(&v)
                + model.assembly().full_stiffness.quadratic_energy(&q);
            near(motion.state().total_energy_j, energy, 2e-15);
            for node in 0..4 {
                let observed = motion.observe(node_point(node)).unwrap();
                for axis in 0..3 {
                    near(observed.displacement_m[axis], q[6 * node + axis], 2e-16);
                    near(observed.velocity_m_s[axis], v[6 * node + axis], 3e-14);
                }
            }
            support_work += step.support_work_j.abs();
            check_balance(&step, 3e-15);
        }
        assert!(
            support_work > 1e-7,
            "the external support must actually exchange work"
        );
    }
}

#[test]
fn changing_support_acceleration_retains_native_full_system_momentum_and_work() {
    let model = supported(true);
    let zero = vec![0.0; model.assembly().mass.dimension()];
    let mut motion = model
        .start_driven(
            &zero,
            &zero,
            ShellSupportMotion::default(),
            &[],
            config(),
            || false,
        )
        .unwrap();
    let mut time = 0.0;
    let mut work = 0.0;
    let mut reacted = 0.0;
    for index in 0..160 {
        time += config().dt_s;
        let phase = index as f64 * 0.07;
        let drive = [
            [0.0, 0.0, phase.sin()],
            [0.0, 0.0, -0.7 * phase.cos()],
            [0.0, 0.0, 0.3 * phase.sin()],
        ];
        let step = motion.advance_driven(time, &[], drive, || false).unwrap();
        work += step.support_work_j - step.damping_work_j;
        reacted += step.support_impulse.linear[2].abs();
        assert_eq!(step.applied_work_j, 0.0);
        check_balance(&step, 5e-15);
    }
    assert!(reacted > 1e-4);
    assert!(motion.state().total_energy_j > 1e-7);
    near(motion.state().total_energy_j, work, 3e-15);
}

#[test]
fn zero_prescribed_motion_keeps_the_stationary_support_path() {
    let model = supported(true);
    let zero = vec![0.0; model.assembly().mass.dimension()];
    let load = wrench(point(1, [0.2, 0.3, 0.5]));
    let mut fixed = model
        .start(&zero, &zero, &[load], config(), || false)
        .unwrap();
    let mut driven = model
        .start_driven(
            &zero,
            &zero,
            ShellSupportMotion::default(),
            &[load],
            config(),
            || false,
        )
        .unwrap();
    for _ in 0..32 {
        let time = fixed.state().time_s + config().dt_s;
        let a = fixed.advance(time, &[load], || false).unwrap();
        let b = driven
            .advance_driven(time, &[load], [[0.0; 3]; 3], || false)
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(fixed.state(), driven.state());
        assert_eq!(fixed.support_reactions(), driven.support_reactions());
    }
}

#[test]
fn driven_refusal_and_late_cancellation_preserve_both_free_and_prescribed_motion() {
    let model = supported(true);
    let (q, v, a) = manufactured(0.0);
    let support = prescribe(&model, &q, &v, &a);
    let loads = free_loads(&model, &mechanical_force(&model, &q, &v, &a));
    let start = || {
        model
            .start_driven(
                &restrict(&model, &q),
                &restrict(&model, &v),
                support,
                &loads,
                config(),
                || false,
            )
            .unwrap()
    };
    let mut control = start();
    let mut polls = 0;
    let expected = control
        .advance_driven(config().dt_s, &loads, support.acceleration_m_s2, || {
            polls += 1;
            false
        })
        .unwrap();
    for fail_at in [1, polls / 2, polls] {
        let mut actual = start();
        let before = actual.state().clone();
        let reactions = actual.support_reactions().to_vec();
        let mut calls = 0;
        assert!(
            actual
                .advance_driven(config().dt_s, &loads, support.acceleration_m_s2, || {
                    calls += 1;
                    calls == fail_at
                })
                .is_err()
        );
        assert_eq!(actual.state(), &before);
        assert_eq!(actual.support_motion(), Some(&support));
        assert_eq!(actual.support_reactions(), reactions);
        assert_eq!(
            actual
                .advance_driven(config().dt_s, &loads, support.acceleration_m_s2, || false)
                .unwrap(),
            expected
        );
        assert_eq!(actual.state(), control.state());
        assert_eq!(actual.support_motion(), control.support_motion());
    }
    let mut actual = start();
    let before = actual.state().clone();
    let mut invalid_acceleration = support.acceleration_m_s2;
    invalid_acceleration[2][1] = f64::NAN;
    assert!(
        actual
            .advance_driven(config().dt_s, &loads, invalid_acceleration, || false)
            .is_err()
    );
    assert!(actual.advance(config().dt_s, &loads, || false).is_err());
    assert!(
        actual
            .advance_driven(
                2.0 * config().dt_s,
                &loads,
                support.acceleration_m_s2,
                || false
            )
            .is_err()
    );
    assert_eq!(actual.state(), &before);
    assert_eq!(actual.support_motion(), Some(&support));
    let free = MovingShellModel::new(plate(), 1, || false).unwrap();
    let zero = vec![0.0; free.assembly().mass.dimension()];
    assert!(
        free.start_driven(&zero, &zero, support, &[], config(), || false)
            .is_err()
    );
}
