//! Physical Coulomb impulses, coupled offsets, prescribed work and atomicity.
use super::*;

fn site(id: u64, point: ShellLoadPoint, mu: f64) -> ShellFrictionalImpactSite {
    ShellFrictionalImpactSite {
        id,
        point,
        frame: [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        activation_tolerance_m: 0.0,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    }
}
fn planes(
    motion: &MovingShellMotion<'_>,
    sites: &[ShellFrictionalImpactSite],
    velocity_m_s: [f64; 3],
) -> Vec<ShellObstacleMotion> {
    sites
        .iter()
        .map(|site| ShellObstacleMotion {
            point_m: motion.observe(site.point).unwrap().position_m,
            velocity_m_s,
        })
        .collect()
}
fn vertex(node: usize) -> ShellLoadPoint {
    let (triangle, barycentric) = match node {
        0 => (0, [1.0, 0.0, 0.0]),
        1 => (0, [0.0, 1.0, 0.0]),
        2 => (0, [0.0, 0.0, 1.0]),
        _ => (1, [0.0, 1.0, 0.0]),
    };
    ShellLoadPoint {
        triangle,
        barycentric,
        offset_m: [0.0; 3],
    }
}

// Independent nodal masses/rotary inertias plus the actual physical COM tensor.
// No assembled inertia, inverse-mass operator or reported momentum is used.
fn inventory(
    motion: &MovingShellMotion<'_>,
    body: Option<ShellRigidPayload>,
) -> (ShellSpatialResultant, f64) {
    let p = plate();
    let mut momentum = ShellSpatialResultant::default();
    let mut energy = 0.0;
    for (node, mass) in [1.0 / 12.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 12.0]
        .into_iter()
        .enumerate()
    {
        let q = motion.observe(vertex(node)).unwrap();
        let linear = q.velocity_m_s.map(|v| mass * v);
        let spin = q
            .angular_velocity_rad_s
            .map(|w| mass * 0.1_f64.powi(2) / 12.0 * w);
        let orbit = cross(p.nodes[node].position_m, linear);
        for a in 0..3 {
            momentum.linear[a] += linear[a];
            momentum.angular[a] += orbit[a] + spin[a];
        }
        energy += 0.5 * (dot(q.velocity_m_s, linear) + dot(q.angular_velocity_rad_s, spin));
    }
    if let Some(body) = body {
        let q = motion.observe(body.point).unwrap();
        let linear = q.velocity_m_s.map(|v| body.mass_kg * v);
        let spin = body
            .inertia_com_kg_m2
            .map(|row| dot(row, q.angular_velocity_rad_s));
        let mut position = body.point.offset_m;
        for (node, w) in p.triangles[body.point.triangle]
            .into_iter()
            .zip(body.point.barycentric)
        {
            for a in 0..3 {
                position[a] += w * p.nodes[node].position_m[a];
            }
        }
        let orbit = cross(position, linear);
        for a in 0..3 {
            momentum.linear[a] += linear[a];
            momentum.angular[a] += orbit[a] + spin[a];
        }
        energy += 0.5 * (dot(q.velocity_m_s, linear) + dot(q.angular_velocity_rad_s, spin));
    }
    (momentum, energy)
}

#[test]
fn offset_coupling_resolves_stick_and_slip_reversal_from_actual_final_velocity() {
    let model = MovingShellModel::new(plate(), 1, || false).unwrap();
    let zero = vec![0.0; model.assembly().mass.dimension()];
    let mut point = vertex(0);
    point.offset_m = [0.05, 0.0, 0.03];
    // m=1/12, I=1/14400. W = (1/m)Id + (r^2 Id - r r^T)/I.
    // In normal-first order, W_nt=-21.6: a normal hit reverses initial x slip.
    for (mu, second_tangent) in [(0.05, 0.0), (0.3, 0.1)] {
        let mut velocity = zero.clone();
        velocity[..3].copy_from_slice(&[0.4, second_tangent, -1.0]);
        let mut motion = model
            .start(&zero, &velocity, &[], config(), || false)
            .unwrap();
        let mut sites = [site(7, point, mu)];
        if mu == 0.3 {
            sites[0].normal_law = NormalVelocityLaw::Target { speed: 0.0 };
        }
        let prepared = motion
            .prepare_frictional_impacts(&sites, Default::default(), || false)
            .unwrap();
        for (actual, expected) in prepared
            .inverse_effective_mass()
            .iter()
            .zip([48.0, -21.6, 0.0, -21.6, 24.96, 0.0, 0.0, 0.0, 60.96])
        {
            near(*actual, expected, 4e-13);
        }
        let obstacles = planes(&motion, &sites, [0.0; 3]);
        let before = motion.state().clone();
        let (initial_momentum, initial_energy) = inventory(&motion, None);
        let report = motion
            .apply_frictional_impacts(&prepared, &obstacles, || false)
            .unwrap();
        let p = if mu == 0.05 {
            [25.0 / 1173.0, 5.0 / 4692.0, 0.0]
        } else {
            [17.0 / 762.0, 5.0 / 1524.0, -second_tangent / 60.96]
        };
        for a in 0..3 {
            near(report.contact.impulses[0][a], p[a], 3e-11);
        }
        if mu == 0.05 {
            // Positive initial slip and POSITIVE friction impulse are correct:
            // the coupled final slip is negative. Initial-sign clipping fails.
            assert!(report.contact.impulses[0][1] > 0.0);
            near(
                report.relative_velocities_after_m_s[0][1],
                -66.0 / 1955.0,
                2e-9,
            );
        } else {
            assert!(p[1].hypot(p[2]) < mu * p[0]);
            for v in report.relative_velocities_after_m_s[0] {
                near(v, 0.0, 2e-9);
            }
        }
        let world = [p[1], p[2], p[0]];
        let torque = cross(point.offset_m, world);
        for a in 0..3 {
            near(
                motion.state().velocity[a],
                velocity[a] + 12.0 * world[a],
                2e-9,
            );
            near(motion.state().velocity[3 + a], 14400.0 * torque[a], 3e-8);
        }
        let (after_momentum, after_energy) = inventory(&motion, None);
        let event = report.impulse.unwrap();
        let point_position = [1.05, -2.0, 0.53];
        let angular = cross(point_position, world);
        for a in 0..3 {
            near(
                after_momentum.linear[a] - initial_momentum.linear[a],
                world[a],
                3e-11,
            );
            near(
                after_momentum.angular[a] - initial_momentum.angular[a],
                angular[a],
                6e-11,
            );
        }
        near(motion.state().kinetic_energy_j, after_energy, 1e-16);
        near(
            event.kinetic_energy_change_j,
            after_energy - initial_energy,
            1e-16,
        );
        near(
            report.dissipated_energy_j,
            initial_energy - after_energy,
            1e-11,
        );
        assert_eq!(event.support_impulse, ShellSpatialResultant::default());
        assert_eq!(report.obstacle_work_j, 0.0);
        assert_eq!(motion.state().displacement, before.displacement);
        assert_eq!(motion.state().time_s, before.time_s);
        assert_eq!(motion.state().step, before.step);
    }
}

#[test]
fn simultaneous_friction_keeps_cross_point_response_id_order_and_tangent_frame_invariance() {
    let model = MovingShellModel::new(plate(), 2, || false).unwrap();
    let zero = vec![0.0; model.assembly().mass.dimension()];
    let mut velocity = zero.clone();
    velocity[0] = 1.0;
    velocity[2] = -0.4;
    velocity[4] = 12.0;
    let sites = [0.05, -0.05].map(|r| {
        let mut p = vertex(0);
        p.offset_m = [r, 0.0, 0.0];
        site(if r > 0.0 { 20 } else { 10 }, p, 0.2)
    });
    let mut reference = None;
    for variant in 0..3 {
        let mut ordered = sites;
        if variant == 1 {
            ordered.reverse();
        }
        if variant == 2 {
            for site in &mut ordered {
                site.frame = [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]];
            }
        }
        let mut motion = model
            .start(&zero, &velocity, &[], config(), || false)
            .unwrap();
        let prepared = motion
            .prepare_frictional_impacts(&ordered, Default::default(), || false)
            .unwrap();
        // Shared node with opposite arms has genuine negative cross-normal mobility.
        near(prepared.inverse_effective_mass()[3], -24.0, 2e-13);
        let obstacles = planes(&motion, &ordered, [0.0; 3]);
        let before_energy = inventory(&motion, None).1;
        let report = motion
            .apply_frictional_impacts(&prepared, &obstacles, || false)
            .unwrap();
        for (i, site) in ordered.iter().enumerate() {
            let pn = if site.id == 20 { 0.025 } else { 1.0 / 120.0 };
            near(report.contact.impulses[i][0], pn, 2e-11);
            assert!(report.contact.impulses[i][0] > 0.0);
            let tangent_axis = if variant == 2 { 2 } else { 1 };
            near(
                report.contact.impulses[i][tangent_axis],
                if variant == 2 { 0.2 * pn } else { -0.2 * pn },
                2e-11,
            );
            near(report.relative_velocities_after_m_s[i][0], 0.0, 2e-9);
        }
        near(motion.state().velocity[0], 0.92, 3e-10);
        near(motion.state().velocity[2], 0.0, 3e-10);
        near(motion.state().velocity[4], 0.0, 3e-9);
        let after_energy = inventory(&motion, None).1;
        near(
            report.dissipated_energy_j,
            before_energy - after_energy,
            3e-11,
        );
        if let Some((state, loss)) = &reference {
            if variant == 1 {
                assert_eq!(motion.state(), state);
            } else {
                for (a, b) in motion.state().velocity.iter().zip(&state.velocity) {
                    near(*a, *b, 3e-9);
                }
            }
            near(report.dissipated_energy_j, *loss, 3e-11);
        } else {
            reference = Some((motion.state().clone(), report.dissipated_energy_j));
        }
    }
}

#[test]
fn payload_and_driven_supports_keep_coulomb_momentum_source_work_and_atomic_refusal() {
    let mut p = plate();
    p.support = Some(ShellSupport {
        node_indices: [2, 0, 1],
        normal: [0.0, 0.0, 1.0],
    });
    p.damping = DampingModel::Rayleigh {
        mass_proportional_per_s: 0.4,
        stiffness_proportional_s: 0.0002,
    };
    let mut point = point(1, [0.2, 0.5, 0.3]);
    point.offset_m = [0.08, -0.03, 0.04];
    let body = ShellRigidPayload {
        point,
        mass_kg: 0.7,
        inertia_com_kg_m2: [
            [0.002, 0.0001, 0.0],
            [0.0001, 0.003, -0.0002],
            [0.0, -0.0002, 0.004],
        ],
    };
    let model = MovingShellModel::new_with_payloads(p, 2, &[body], || false).unwrap();
    let support = ShellSupportMotion {
        velocity_m_s: [[0.03, -0.02, -0.1]; 3],
        ..Default::default()
    };
    let zero = vec![0.0; model.assembly().mass.dimension()];
    let velocity: Vec<_> = model
        .assembly()
        .free_dofs
        .iter()
        .map(|row| {
            if row % 6 < 3 {
                support.velocity_m_s[0][row % 6]
            } else {
                0.0
            }
        })
        .collect();
    let loads = [wrench(point)];
    let start = || {
        model
            .start_driven(&zero, &velocity, support, &loads, config(), || false)
            .unwrap()
    };
    let mut motion = start();
    let sites = [site(9, point, 0.25)];
    let prepared = motion
        .prepare_frictional_impacts(&sites, Default::default(), || false)
        .unwrap();
    let surface = [0.35, 0.12, 0.2];
    let obstacles = planes(&motion, &sites, surface);
    let before = motion.state().clone();
    let reactions = motion.support_reactions().to_vec();
    let (before_momentum, before_energy) = inventory(&motion, Some(body));
    let mut invalid = sites;
    invalid[0].frame[2] = invalid[0].frame[1];
    assert!(
        motion
            .prepare_frictional_impacts(&invalid, Default::default(), || false)
            .is_err()
    );
    assert!(
        motion
            .prepare_frictional_impacts(&[sites[0]; 2], Default::default(), || false)
            .is_err()
    );
    assert!(
        motion
            .apply_frictional_impacts(&prepared, &[], || false)
            .is_err()
    );
    let mut penetrating = obstacles.clone();
    penetrating[0].point_m[2] += 0.001;
    assert!(matches!(
        motion.apply_frictional_impacts(&prepared, &penetrating, || false),
        Err(ShellFrictionalImpactError::Contact(
            fs_contact::frictional_multicontact::FrictionalContactError::PenetratingContact { .. }
        ))
    ));
    let strict = motion
        .prepare_frictional_impacts(
            &sites,
            ShellFrictionalImpactControls {
                speed_tolerance_m_s: 0.0,
                work_tolerance_j: 0.0,
                ..Default::default()
            },
            || false,
        )
        .unwrap();
    assert!(matches!(
        motion.apply_frictional_impacts(&strict, &obstacles, || false),
        Err(ShellFrictionalImpactError::Audit { .. })
    ));
    let mut control = start();
    let mut polls = 0;
    let report = control
        .apply_frictional_impacts(&prepared, &obstacles, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        motion.apply_frictional_impacts(&prepared, &obstacles, || {
            calls += 1;
            calls == polls
        }),
        Err(ShellFrictionalImpactError::Cancelled)
    ));
    assert_eq!(motion.state(), &before);
    assert_eq!(motion.support_reactions(), reactions);
    assert_eq!(motion.support_motion(), Some(&support));
    assert_eq!(
        motion
            .apply_frictional_impacts(&prepared, &obstacles, || false)
            .unwrap(),
        report
    );
    assert_eq!(motion.state(), control.state());
    let (after_momentum, after_energy) = inventory(&motion, Some(body));
    let p = report.contact.impulses[0];
    let world = [p[1], p[2], p[0]];
    let position = [1.7 + 0.08, -1.2 - 0.03, 0.5 + 0.04];
    let applied_angular = cross(position, world);
    let reaction =
        std::array::from_fn(|a| after_momentum.linear[a] - before_momentum.linear[a] - world[a]);
    let event = report.impulse.unwrap();
    for a in 0..3 {
        near(event.support_impulse.linear[a], reaction[a], 4e-15);
        near(
            event.support_impulse.angular[a],
            after_momentum.angular[a] - before_momentum.angular[a] - applied_angular[a],
            8e-15,
        );
    }
    let obstacle_work = dot(world, surface);
    let support_work = dot(reaction, support.velocity_m_s[0]);
    near(report.obstacle_work_j, obstacle_work, 2e-16);
    near(event.support_work_j, support_work, 2e-16);
    assert!(obstacle_work.abs() > 1e-4 && support_work.abs() > 1e-5);
    assert!(world[0].abs() > 1e-4 && world[1].abs() > 1e-5);
    near(motion.state().kinetic_energy_j, after_energy, 3e-16);
    near(
        event.kinetic_energy_change_j,
        after_energy - before_energy,
        3e-16,
    );
    near(
        after_energy - before_energy,
        obstacle_work + support_work + report.relative_work_j,
        3e-16,
    );
    let actual = motion.observe(point).unwrap().velocity_m_s;
    let relative = [
        actual[2] - surface[2],
        actual[0] - surface[0],
        actual[1] - surface[1],
    ];
    for a in 0..3 {
        near(
            report.relative_velocities_after_m_s[0][a],
            relative[a],
            2e-16,
        );
    }
    let speed = relative[1].hypot(relative[2]);
    assert!(speed > 1e-3);
    near(
        p[1],
        -sites[0].friction_coefficient * p[0] * relative[1] / speed,
        2e-10,
    );
    near(
        p[2],
        -sites[0].friction_coefficient * p[0] * relative[2] / speed,
        2e-10,
    );
    let mut fresh = model
        .start_driven(
            &motion.state().displacement,
            &motion.state().velocity,
            support,
            &loads,
            config(),
            || false,
        )
        .unwrap();
    assert_eq!(motion.state(), fresh.state());
    assert_eq!(motion.support_reactions(), fresh.support_reactions());
    assert_eq!(
        motion
            .advance_driven(config().dt_s, &loads, support.acceleration_m_s2, || false)
            .unwrap(),
        fresh
            .advance_driven(config().dt_s, &loads, support.acceleration_m_s2, || false)
            .unwrap()
    );
}
