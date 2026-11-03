//! G0/G1/G3/G4: physical moving-wrench ports and native shell motion.
use fs_solid::moving_load::*;
use fs_solid::shell::*;
use fs_solid::shell_time::ShellTimeConfig;

#[path = "moving_load/driven.rs"]
mod driven;

fn plate() -> ShellPlate {
    ShellPlate {
        nodes: vec![
            [1.0, -2.0, 0.5],
            [2.0, -2.0, 0.5],
            [1.0, -1.0, 0.5],
            [2.0, -1.0, 0.5],
        ]
        .into_iter()
        .map(|position_m| ShellNode { position_m })
        .collect(),
        triangles: vec![[0, 1, 2], [1, 3, 2]],
        thickness_m: 0.1,
        material: ShellMaterial {
            youngs_modulus_pa: 1000.0,
            poisson_ratio: 0.25,
            density_kg_m3: 5.0,
        },
        identity: ShellIdentity {
            model_id: "authored-shell".into(),
            source_id: "explicit-synthetic-material".into(),
            state_id: "reference".into(),
        },
        support: None,
        damping: DampingModel::None,
        budget: AssemblyBudget::default(),
    }
}
fn point(triangle: usize, barycentric: [f64; 3]) -> ShellLoadPoint {
    ShellLoadPoint {
        triangle,
        barycentric,
        offset_m: [0.02, -0.01, 0.05],
    }
}
fn wrench(point: ShellLoadPoint) -> ShellPointWrench {
    ShellPointWrench {
        point,
        force_n: [0.3, -0.4, 0.7],
        couple_nm: [0.02, -0.01, 0.03],
    }
}
fn config() -> ShellTimeConfig {
    ShellTimeConfig {
        dt_s: 1e-4,
        duration_s: 0.1,
        max_steps: 1000,
        beta: 0.25,
        gamma: 0.5,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; error {:e}",
        a - b
    );
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}

#[test]
fn distributed_wrenches_keep_force_couple_and_reference_lever_arms() {
    let model = MovingShellModel::new(plate(), 4, || false).unwrap();
    let mut second = wrench(point(1, [0.2, 0.3, 0.5]));
    second.force_n = [-0.2, 0.8, -0.1];
    second.couple_nm = [-0.07, 0.01, 0.02];
    let loads = [wrench(point(0, [0.1, 0.6, 0.3])), second];
    let projection = model.project(&loads, || false).unwrap();
    let (mut force, mut moment) = ([0.0; 3], [0.0; 3]);
    for load in loads {
        let mut position = load.point.offset_m;
        for (node, weight) in model.plate().triangles[load.point.triangle]
            .into_iter()
            .zip(load.point.barycentric)
        {
            for axis in 0..3 {
                position[axis] += weight * model.plate().nodes[node].position_m[axis];
            }
        }
        let orbital = cross(position, load.force_n);
        for axis in 0..3 {
            force[axis] += load.force_n[axis];
            moment[axis] += orbital[axis] + load.couple_nm[axis];
        }
    }
    for axis in 0..3 {
        near(projection.resultant.linear[axis], force[axis], 3e-16);
        near(projection.resultant.angular[axis], moment[axis], 1e-15);
    }
    assert_eq!(projection.full, projection.reduced);
}

#[test]
fn rigid_translation_rotation_and_offset_material_velocity_are_reproduced() {
    let model = MovingShellModel::new(plate(), 2, || false).unwrap();
    let translation = [0.002, -0.003, 0.001];
    let rotation = [0.0001, 0.0002, -0.0003];
    let speed = [0.2, 0.3, -0.1];
    let omega = [-0.7, 0.4, 0.2];
    let (mut q, mut v) = (Vec::new(), Vec::new());
    for node in &model.plate().nodes {
        let arm = cross(rotation, node.position_m);
        let rotational_velocity = cross(omega, node.position_m);
        q.extend((0..3).map(|i| translation[i] + arm[i]));
        q.extend(rotation);
        v.extend((0..3).map(|i| speed[i] + rotational_velocity[i]));
        v.extend(omega);
    }
    for sample in [point(0, [0.3, 0.2, 0.5]), point(1, [0.0, 1.0, 0.0])] {
        let observed = model.observe(sample, &q, &v).unwrap();
        let arm = cross(rotation, observed.reference_position_m);
        let spin = cross(omega, observed.reference_position_m);
        for axis in 0..3 {
            near(
                observed.displacement_m[axis],
                translation[axis] + arm[axis],
                1e-18,
            );
            near(observed.velocity_m_s[axis], speed[axis] + spin[axis], 8e-16);
            near(observed.angular_velocity_rad_s[axis], omega[axis], 2e-16);
        }
        near(dot(observed.normal, observed.normal), 1.0, 4e-16);
    }
}

#[test]
fn moving_patch_power_is_material_work_not_the_convective_path_derivative() {
    let model = MovingShellModel::new(plate(), 2, || false).unwrap();
    let n = model.assembly().free_dofs.len();
    let q: Vec<_> = (0..n).map(|i| 1e-4 * ((i as f64) * 0.6).sin()).collect();
    let v: Vec<_> = (0..n).map(|i| 0.3 * ((i as f64) * 0.7).cos()).collect();
    let loads = [
        wrench(point(0, [0.1, 0.6, 0.3])),
        wrench(point(1, [0.2, 0.3, 0.5])),
    ];
    let projected = model.project(&loads, || false).unwrap();
    let generalized_power: f64 = projected.reduced.iter().zip(&v).map(|(f, v)| f * v).sum();
    let physical_power: f64 = loads
        .iter()
        .map(|load| {
            let sample = model.observe(load.point, &q, &v).unwrap();
            dot(load.force_n, sample.velocity_m_s)
                + dot(load.couple_nm, sample.angular_velocity_rad_s)
        })
        .sum();
    near(generalized_power, physical_power, 2e-16);
    let stationary = vec![0.0; n];
    let left = model.observe(loads[0].point, &q, &stationary).unwrap();
    let right = model.observe(loads[1].point, &q, &stationary).unwrap();
    assert_ne!(left.position_m, right.position_m);
    assert_eq!(left.velocity_m_s, [0.0; 3]);
    assert_eq!(right.velocity_m_s, [0.0; 3]);
    // Either original face gives the same interpolation at their common edge.
    let a = model
        .project(&[wrench(point(0, [0.0, 0.25, 0.75]))], || false)
        .unwrap();
    let b = model
        .project(&[wrench(point(1, [0.25, 0.0, 0.75]))], || false)
        .unwrap();
    assert_eq!(a, b);
}

#[test]
fn native_trajectory_carries_moving_force_and_moment_impulses_with_damping() {
    for damped in [false, true] {
        let mut p = plate();
        if damped {
            p.damping = DampingModel::Rayleigh {
                mass_proportional_per_s: 0.3,
                stiffness_proportional_s: 0.0001,
            };
        }
        let model = MovingShellModel::new(p, 2, || false).unwrap();
        let n = model.assembly().mass.dimension();
        let zeros = vec![0.0; n];
        let mut previous = wrench(point(0, [0.2, 0.3, 0.5]));
        let mut motion = model
            .start(&zeros, &zeros, &[previous], config(), || false)
            .unwrap();
        let mut expected_linear = [0.0; 3];
        let mut measured_linear = [0.0; 3];
        for sample in 1..=120 {
            let t = sample as f64 / 120.0;
            let mut next = wrench(point(
                usize::from(sample > 60),
                [0.2, 0.3 + 0.1 * t, 0.5 - 0.1 * t],
            ));
            next.force_n = [0.3 * (1.0 - t), -0.4 * (1.0 + t), 0.7];
            next.couple_nm = [0.02 * t, -0.01, 0.03 * (1.0 - t)];
            let before_velocity = motion.state().velocity.clone();
            let report = motion
                .advance(motion.state().time_s + config().dt_s, &[next], || false)
                .unwrap();
            // Independent native-mass inventory, not the new report's reduction.
            let delta: Vec<_> = motion
                .state()
                .velocity
                .iter()
                .zip(&before_velocity)
                .map(|(a, b)| a - b)
                .collect();
            let momentum = model.assembly().mass.apply(&delta);
            for axis in 0..3 {
                expected_linear[axis] +=
                    0.5 * config().dt_s * (previous.force_n[axis] + next.force_n[axis])
                        - report.damping_impulse.linear[axis];
                measured_linear[axis] += momentum.chunks_exact(6).map(|p| p[axis]).sum::<f64>();
                near(report.momentum_residual.linear[axis], 0.0, 2e-15);
                near(report.momentum_residual.angular[axis], 0.0, 5e-15);
            }
            near(report.energy_residual_j, 0.0, 2e-16);
            assert!(report.damping_work_j >= 0.0);
            if damped {
                assert!(report.damping_work_j > 0.0);
            }
            previous = next;
        }
        for axis in 0..3 {
            near(measured_linear[axis], expected_linear[axis], 1e-14);
        }
    }
}

#[test]
fn initially_equilibrated_acceleration_recovers_analytic_rigid_translation() {
    let model = MovingShellModel::new(plate(), 4, || false).unwrap();
    let n = model.assembly().mass.dimension();
    // Apply F_a = m_a * a directly through vertex ports. The exact native
    // shell rigid mode then has q=a*t²/2, v=a*t; no elastic approximation.
    let acceleration = [0.3, -0.1, 0.2];
    let mass = model.assembly().full_mass.values();
    let loads: Vec<_> = (0..4)
        .map(|node| {
            let (triangle, barycentric) = match node {
                0 => (0, [1., 0., 0.]),
                1 => (0, [0., 1., 0.]),
                2 => (0, [0., 0., 1.]),
                _ => (1, [0., 1., 0.]),
            };
            ShellPointWrench {
                point: ShellLoadPoint {
                    triangle,
                    barycentric,
                    offset_m: [0.; 3],
                },
                force_n: acceleration.map(|a| mass[(6 * node) * n + 6 * node] * a),
                couple_nm: [0.; 3],
            }
        })
        .collect();
    let zeros = vec![0.; n];
    let mut motion = model
        .start(&zeros, &zeros, &loads, config(), || false)
        .unwrap();
    for _ in 0..40 {
        motion
            .advance(motion.state().time_s + config().dt_s, &loads, || false)
            .unwrap();
        let t = motion.state().time_s;
        for node in 0..4 {
            for axis in 0..3 {
                near(
                    motion.state().velocity[6 * node + axis],
                    acceleration[axis] * t,
                    2e-17,
                );
                near(
                    motion.state().displacement[6 * node + axis],
                    0.5 * acceleration[axis] * t * t,
                    2e-20,
                );
            }
        }
    }
}

#[test]
fn supported_rows_receive_the_actual_load_and_reaction_without_spurious_motion() {
    let mut p = plate();
    p.nodes.truncate(3);
    p.triangles.truncate(1);
    p.support = Some(ShellSupport {
        node_indices: [0, 1, 2],
        normal: [0., 0., 1.],
    });
    let model = MovingShellModel::new(p, 1, || false).unwrap();
    let load = ShellPointWrench {
        point: ShellLoadPoint {
            triangle: 0,
            barycentric: [0.2, 0.3, 0.5],
            offset_m: [0.; 3],
        },
        force_n: [0.3, -0.4, 0.7],
        couple_nm: [0.; 3],
    };
    let projected = model.project(&[load], || false).unwrap();
    assert!(projected.reduced.iter().all(|f| *f == 0.));
    let zeros = vec![0.; model.assembly().mass.dimension()];
    let mut motion = model
        .start(&zeros, &zeros, &[load], config(), || false)
        .unwrap();
    let step = motion.advance(config().dt_s, &[load], || false).unwrap();
    for (force, reaction) in projected.full.iter().zip(motion.support_reactions()) {
        assert_eq!(*force, -reaction);
    }
    for axis in 0..3 {
        near(
            step.applied_impulse.linear[axis] + step.support_impulse.linear[axis],
            0.,
            1e-20,
        );
        near(
            step.applied_impulse.angular[axis] + step.support_impulse.angular[axis],
            0.,
            1e-19,
        );
    }
    assert_eq!(motion.state().displacement, zeros);
    assert_eq!(motion.state().velocity, zeros);
    assert_eq!(step.applied_work_j, 0.);
}

#[test]
fn invalid_cancelled_or_replayed_samples_preserve_the_complete_accepted_motion() {
    let model = MovingShellModel::new(plate(), 1, || false).unwrap();
    let n = model.assembly().mass.dimension();
    let zeros = vec![0.; n];
    let load = wrench(point(0, [0.2, 0.3, 0.5]));
    let mut motion = model
        .start(&zeros, &zeros, &[load], config(), || false)
        .unwrap();
    let original = motion.state().clone();
    let reactions = motion.support_reactions().to_vec();
    for time in [0., 2. * config().dt_s, f64::NAN] {
        assert!(motion.advance(time, &[load], || false).is_err());
    }
    assert!(
        motion
            .advance(config().dt_s, &[load, load], || false)
            .is_err()
    );
    let mut bad = load;
    bad.point.barycentric = [1.2, -0.2, 0.];
    assert!(motion.advance(config().dt_s, &[bad], || false).is_err());
    bad = load;
    bad.point.triangle = 2;
    assert!(motion.advance(config().dt_s, &[bad], || false).is_err());
    let mut polls = 0;
    assert!(matches!(
        motion.advance(config().dt_s, &[load], || {
            polls += 1;
            polls == 6
        }),
        Err(MovingLoadError::Cancelled)
    ));
    assert_eq!(motion.state(), &original);
    assert_eq!(motion.support_reactions(), reactions);
    let first = motion.advance(config().dt_s, &[load], || false).unwrap();
    let accepted = motion.state().clone();
    assert!(motion.advance(config().dt_s, &[load], || false).is_err());
    assert_eq!(motion.state(), &accepted);
    let mut independent = model
        .start(&zeros, &zeros, &[load], config(), || false)
        .unwrap();
    assert_eq!(
        independent
            .advance(config().dt_s, &[load], || false)
            .unwrap(),
        first
    );
    assert_eq!(motion.state(), independent.state());
}
