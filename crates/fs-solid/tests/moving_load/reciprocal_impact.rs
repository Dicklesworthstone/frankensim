//! Two original finite shell inertias, actual reciprocal work and atomic motion.
use super::*;

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
fn site(id: u64, point: ShellLoadPoint, mu: f64, restitution: f64) -> ShellReciprocalImpactSite {
    ShellReciprocalImpactSite {
        first: ShellFrictionalImpactSite {
            id,
            point,
            frame: [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            activation_tolerance_m: 0.0,
            friction_coefficient: mu,
            normal_law: NormalVelocityLaw::NewtonRestitution {
                coefficient: restitution,
            },
        },
        second: point,
    }
}
fn model(scale: f64) -> MovingShellModel {
    let mut p = plate();
    p.material.density_kg_m3 *= scale;
    MovingShellModel::new(p, 2, || false).unwrap()
}
fn inventory(
    motion: &MovingShellMotion<'_>,
    scale: f64,
    body: Option<ShellRigidPayload>,
) -> (ShellSpatialResultant, f64) {
    let p = plate();
    let mut result = ShellSpatialResultant::default();
    let mut energy = 0.0;
    for (node, mass) in [1.0 / 12.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 12.0]
        .into_iter()
        .enumerate()
    {
        let q = motion.observe(vertex(node)).unwrap();
        let linear = q.velocity_m_s.map(|v| mass * scale * v);
        let spin = q.angular_velocity_rad_s.map(|w| mass * scale / 1200.0 * w);
        let orbit = cross(p.nodes[node].position_m, linear);
        for a in 0..3 {
            result.linear[a] += linear[a];
            result.angular[a] += orbit[a] + spin[a];
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
            result.linear[a] += linear[a];
            result.angular[a] += orbit[a] + spin[a];
        }
        energy += 0.5 * (dot(q.velocity_m_s, linear) + dot(q.angular_velocity_rad_s, spin));
    }
    (result, energy)
}

#[test]
fn unequal_finite_shells_recoil_with_analytic_restitution_and_no_hidden_reference_couple() {
    let a = model(1.0);
    let b = model(3.0);
    let zero = vec![0.0; a.assembly().mass.dimension()];
    let mut va = zero.clone();
    va[2] = -1.0;
    let mut vb = zero.clone();
    vb[2] = 0.2;
    for e in [0.0, 0.5, 1.0] {
        let mut first = a.start(&zero, &va, &[], config(), || false).unwrap();
        let mut second = b.start(&zero, &vb, &[], config(), || false).unwrap();
        let prepared = first
            .prepare_reciprocal_impacts(
                &mut second,
                &[site(1, vertex(0), 0.0, e)],
                Default::default(),
                || false,
            )
            .unwrap();
        near(prepared.inverse_effective_mass()[0], 16.0, 3e-14);
        let before = [inventory(&first, 1.0, None), inventory(&second, 3.0, None)];
        let report = first
            .apply_reciprocal_impacts(&mut second, &prepared, || false)
            .unwrap();
        // Reduced mass = (1/12)*(1/4)/(1/12+1/4) = 1/16 kg.
        let impulse = 0.075 * (1.0 + e);
        near(report.contact.impulses[0][0], impulse, 3e-13);
        near(first.state().velocity[2], -0.1 + 0.9 * e, 3e-12);
        near(second.state().velocity[2], -0.1 - 0.3 * e, 3e-12);
        near(report.relative_velocities_after_m_s[0][0], 1.2 * e, 3e-12);
        near(report.dissipated_energy_j, 0.045 * (1.0 - e * e), 3e-13);
        let after = [inventory(&first, 1.0, None), inventory(&second, 3.0, None)];
        for axis in 0..3 {
            near(
                after[0].0.linear[axis] + after[1].0.linear[axis],
                before[0].0.linear[axis] + before[1].0.linear[axis],
                4e-16,
            );
            near(
                after[0].0.angular[axis] + after[1].0.angular[axis],
                before[0].0.angular[axis] + before[1].0.angular[axis],
                8e-16,
            );
        }
        near(
            report.kinetic_energy_change_j,
            after[0].1 + after[1].1 - before[0].1 - before[1].1,
            2e-16,
        );
        assert_eq!(report.support_work_j, 0.0);
        for event in report.impulses.unwrap() {
            assert_eq!(event.support_impulse, ShellSpatialResultant::default());
        }
    }
    // Generous geometric admission must not hide a finite angular couple.
    let mut first = a.start(&zero, &va, &[], config(), || false).unwrap();
    let mut second = b.start(&zero, &vb, &[], config(), || false).unwrap();
    let initial = [first.state().clone(), second.state().clone()];
    let mut offset = site(2, vertex(0), 0.0, 0.0);
    offset.second.offset_m[0] = 1e-4;
    let prepared = first
        .prepare_reciprocal_impacts(
            &mut second,
            &[offset],
            ShellReciprocalImpactControls {
                reference_point_tolerance_m: 1e-3,
                tangential_position_tolerance_m: 1e-3,
                ..Default::default()
            },
            || false,
        )
        .unwrap();
    assert!(matches!(
        first.apply_reciprocal_impacts(&mut second, &prepared, || false),
        Err(ShellFrictionalImpactError::Audit {
            quantity: "reference offset contact couple",
            ..
        })
    ));
    assert_eq!(first.state(), &initial[0]);
    assert_eq!(second.state(), &initial[1]);
}

#[test]
fn simultaneous_offset_friction_transfers_both_spin_inventories_with_id_and_frame_invariance() {
    let a = model(1.0);
    let b = model(2.0);
    let zero = vec![0.0; a.assembly().mass.dimension()];
    let mut velocity = zero.clone();
    velocity[0] = 1.0;
    velocity[2] = -0.4;
    velocity[4] = 12.0;
    let sites = [0.05, -0.05].map(|r| {
        let mut point = vertex(0);
        point.offset_m = [r, 0.0, 0.0];
        site(if r > 0.0 { 20 } else { 10 }, point, 0.2, 0.0)
    });
    let mut reference: Option<([fs_solid::shell_time::DynamicState; 2], f64)> = None;
    for variant in 0..3 {
        let mut ordered = sites;
        if variant == 1 {
            ordered.reverse();
        }
        if variant == 2 {
            for site in &mut ordered {
                site.first.frame = [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]];
            }
        }
        let mut first = a.start(&zero, &velocity, &[], config(), || false).unwrap();
        let mut second = b.start(&zero, &zero, &[], config(), || false).unwrap();
        let before = [inventory(&first, 1.0, None), inventory(&second, 2.0, None)];
        let prepared = first
            .prepare_reciprocal_impacts(&mut second, &ordered, Default::default(), || false)
            .unwrap();
        // W_relative=(1+1/2)*[[48,-24],[-24,48]] on the normal block.
        near(prepared.inverse_effective_mass()[0], 72.0, 3e-13);
        near(prepared.inverse_effective_mass()[3], -36.0, 3e-13);
        let report = first
            .apply_reciprocal_impacts(&mut second, &prepared, || false)
            .unwrap();
        for (i, site) in ordered.iter().enumerate() {
            let pn = if site.first.id == 20 {
                1.0 / 60.0
            } else {
                1.0 / 180.0
            };
            near(report.contact.impulses[i][0], pn, 2e-11);
            let axis = if variant == 2 { 2 } else { 1 };
            near(
                report.contact.impulses[i][axis],
                if variant == 2 { 0.2 * pn } else { -0.2 * pn },
                2e-11,
            );
            near(report.relative_velocities_after_m_s[i][0], 0.0, 2e-9);
        }
        near(first.state().velocity[0], 71.0 / 75.0, 3e-10);
        near(second.state().velocity[0], 2.0 / 75.0, 3e-10);
        for motion in [&first, &second] {
            near(motion.state().velocity[2], -2.0 / 15.0, 3e-10);
            near(motion.state().velocity[4], 4.0, 3e-9);
        }
        let after = [inventory(&first, 1.0, None), inventory(&second, 2.0, None)];
        for axis in 0..3 {
            near(
                after[0].0.linear[axis] + after[1].0.linear[axis],
                before[0].0.linear[axis] + before[1].0.linear[axis],
                3e-15,
            );
            near(
                after[0].0.angular[axis] + after[1].0.angular[axis],
                before[0].0.angular[axis] + before[1].0.angular[axis],
                5e-15,
            );
        }
        near(
            report.dissipated_energy_j,
            before[0].1 + before[1].1 - after[0].1 - after[1].1,
            2e-14,
        );
        if let Some((states, loss)) = &reference {
            if variant == 1 {
                assert_eq!([first.state(), second.state()], [&states[0], &states[1]]);
            } else {
                for (motion, state) in [&first, &second].into_iter().zip(states) {
                    for (x, y) in motion.state().velocity.iter().zip(&state.velocity) {
                        near(*x, *y, 3e-9);
                    }
                }
            }
            near(report.dissipated_energy_j, *loss, 2e-14);
        } else {
            reference = Some((
                [first.state().clone(), second.state().clone()],
                report.dissipated_energy_j,
            ));
        }
    }
}

#[test]
fn deformed_driven_payload_pair_keeps_both_support_work_and_late_refusals_atomic() {
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
    let bodies = [
        body,
        ShellRigidPayload {
            mass_kg: 2.0 * body.mass_kg,
            inertia_com_kg_m2: body.inertia_com_kg_m2.map(|r| r.map(|v| 2.0 * v)),
            ..body
        },
    ];
    let models = [1.0, 2.0].map(|scale| {
        let mut p = plate();
        p.material.density_kg_m3 *= scale;
        p.material.youngs_modulus_pa *= scale;
        p.support = Some(ShellSupport {
            node_indices: [2, 0, 1],
            normal: [0.0, 0.0, 1.0],
        });
        // Stiffness damping permits an independent uniform translational boost.
        p.damping = DampingModel::Rayleigh {
            mass_proportional_per_s: 0.0,
            stiffness_proportional_s: 0.0002,
        };
        MovingShellModel::new_with_payloads(
            p,
            2,
            &[bodies[if scale == 1.0 { 0 } else { 1 }]],
            || false,
        )
        .unwrap()
    });
    let base = [0.03, -0.02, -0.1];
    let boost = [0.2, -0.1, 0.15];
    let dt = config().dt_s;
    let loads = [
        vec![wrench(point)],
        vec![ShellPointWrench {
            force_n: wrench(point).force_n.map(|f| 2.0 * f),
            couple_nm: wrench(point).couple_nm.map(|f| 2.0 * f),
            point,
        }],
    ];
    let start = |index: usize| {
        let shift: [f64; 3] = std::array::from_fn(|a| {
            if index == 0 {
                if a == 2 { 2e-10 } else { 0.0 }
            } else {
                -dt * boost[a]
            }
        });
        let v = std::array::from_fn(|a| base[a] + if index == 1 { boost[a] } else { 0.0 });
        let support = ShellSupportMotion {
            displacement_m: [[0.001 + shift[0], -0.0005 + shift[1], 0.0003 + shift[2]]; 3],
            velocity_m_s: [v; 3],
            ..Default::default()
        };
        let q: Vec<_> = models[index]
            .assembly()
            .free_dofs
            .iter()
            .map(|r| {
                if r % 6 < 3 {
                    [0.001, -0.0005, 0.0003][r % 6] + 0.0002 + shift[r % 6]
                } else {
                    0.002 * (0.7 * (*r as f64)).cos()
                }
            })
            .collect();
        let velocity: Vec<_> = models[index]
            .assembly()
            .free_dofs
            .iter()
            .map(|r| if r % 6 < 3 { v[r % 6] } else { 0.0 })
            .collect();
        let mut motion = models[index]
            .start_driven(&q, &velocity, support, &loads[index], config(), || false)
            .unwrap();
        motion
            .advance_driven(dt, &loads[index], [[0.0; 3]; 3], || false)
            .unwrap();
        if index == 0 {
            motion
                .apply_impulses(
                    &[ShellPointImpulse {
                        point,
                        linear_impulse_ns: [-0.005, -0.003, -0.04],
                        angular_impulse_nms: [0.0; 3],
                    }],
                    || false,
                )
                .unwrap();
        }
        motion
    };
    let mut first = start(0);
    let mut second = start(1);
    assert_eq!(first.state().time_s, dt);
    assert_eq!(second.state().time_s, dt);
    assert!(first.state().displacement.iter().any(|v| v.abs() > 1e-4));
    let mut sites = [site(3, point, 0.25, 0.0)];
    sites[0].first.activation_tolerance_m = 1e-8;
    let prepared = first
        .prepare_reciprocal_impacts(&mut second, &sites, Default::default(), || false)
        .unwrap();
    let before = [
        inventory(&first, 1.0, Some(bodies[0])),
        inventory(&second, 2.0, Some(bodies[1])),
    ];
    let initial = [first.state().clone(), second.state().clone()];
    let reactions = [
        first.support_reactions().to_vec(),
        second.support_reactions().to_vec(),
    ];
    let supports = [
        *first.support_motion().unwrap(),
        *second.support_motion().unwrap(),
    ];
    let strict = first
        .prepare_reciprocal_impacts(
            &mut second,
            &sites,
            ShellReciprocalImpactControls {
                impact: ShellFrictionalImpactControls {
                    speed_tolerance_m_s: 0.0,
                    work_tolerance_j: 0.0,
                    ..Default::default()
                },
                ..Default::default()
            },
            || false,
        )
        .unwrap();
    assert!(matches!(
        first.apply_reciprocal_impacts(&mut second, &strict, || false),
        Err(ShellFrictionalImpactError::Audit { .. })
    ));
    let mut control_a = start(0);
    let mut control_b = start(1);
    let mut polls = 0;
    let report = control_a
        .apply_reciprocal_impacts(&mut control_b, &prepared, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        first.apply_reciprocal_impacts(&mut second, &prepared, || {
            calls += 1;
            calls == polls
        }),
        Err(ShellFrictionalImpactError::Cancelled)
    ));
    for (i, motion) in [&first, &second].into_iter().enumerate() {
        assert_eq!(motion.state(), &initial[i]);
        assert_eq!(motion.support_reactions(), reactions[i]);
        assert_eq!(motion.support_motion(), Some(&supports[i]));
    }
    assert_eq!(
        first
            .apply_reciprocal_impacts(&mut second, &prepared, || false)
            .unwrap(),
        report
    );
    assert_eq!(first.state(), control_a.state());
    assert_eq!(second.state(), control_b.state());
    let after = [
        inventory(&first, 1.0, Some(bodies[0])),
        inventory(&second, 2.0, Some(bodies[1])),
    ];
    let p = report.contact.impulses[0];
    let world = [p[1], p[2], p[0]];
    let angular = cross([1.78, -1.23, 0.54], world);
    let mut support_work = 0.0;
    for i in 0..2 {
        let sign = if i == 0 { 1.0 } else { -1.0 };
        let reaction =
            std::array::from_fn(|a| after[i].0.linear[a] - before[i].0.linear[a] - sign * world[a]);
        let event = report.impulses.unwrap()[i];
        for a in 0..3 {
            near(event.support_impulse.linear[a], reaction[a], 8e-15);
            near(
                event.support_impulse.angular[a],
                after[i].0.angular[a] - before[i].0.angular[a] - sign * angular[a],
                2e-14,
            );
        }
        let work = dot(reaction, supports[i].velocity_m_s[0]);
        support_work += work;
        near(event.support_work_j, work, 3e-16);
        assert!(work.abs() > 1e-5);
        near(
            event.kinetic_energy_change_j,
            after[i].1 - before[i].1,
            4e-16,
        );
    }
    assert!(support_work.abs() > 1e-5);
    near(report.support_work_j, support_work, 4e-16);
    near(
        report.kinetic_energy_change_j,
        after[0].1 + after[1].1 - before[0].1 - before[1].1,
        5e-16,
    );
    near(
        report.kinetic_energy_change_j,
        support_work + report.relative_work_j,
        5e-16,
    );
    for (index, motion) in [&mut first, &mut second].into_iter().enumerate() {
        let mut fresh = models[index]
            .start_driven(
                &motion.state().displacement,
                &motion.state().velocity,
                supports[index],
                &loads[index],
                config(),
                || false,
            )
            .unwrap();
        assert_eq!(motion.state().acceleration, fresh.state().acceleration);
        assert_eq!(motion.support_reactions(), fresh.support_reactions());
        motion
            .advance_driven(2.0 * dt, &loads[index], [[0.0; 3]; 3], || false)
            .unwrap();
        fresh
            .advance_driven(dt, &loads[index], [[0.0; 3]; 3], || false)
            .unwrap();
        assert_eq!(motion.state().velocity, fresh.state().velocity);
        assert_eq!(motion.state().displacement, fresh.state().displacement);
        assert_eq!(motion.support_reactions(), fresh.support_reactions());
    }

    // The first native trial succeeds, but the finite massive second body's
    // tiny recoil cannot be represented atop its accepted speed. Neither commits.
    let small = model(1.0);
    let huge = model(2.0_f64.powi(60));
    let zero = vec![0.0; small.assembly().mass.dimension()];
    let mut va = zero.clone();
    let mut vb = zero.clone();
    va[0] = 2.0_f64.powi(40) - 1.0;
    vb[0] = 2.0_f64.powi(40);
    let mut a = small.start(&zero, &va, &[], config(), || false).unwrap();
    let mut b = huge.start(&zero, &vb, &[], config(), || false).unwrap();
    let mut probe = small.start(&zero, &va, &[], config(), || false).unwrap();
    probe
        .apply_impulses(
            &[ShellPointImpulse {
                point: vertex(0),
                linear_impulse_ns: [1.0 / 12.0, 0.0, 0.0],
                angular_impulse_nms: [0.0; 3],
            }],
            || false,
        )
        .unwrap();
    let mut pair = site(1, vertex(0), 0.0, 0.0);
    pair.first.frame = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let prepared = a
        .prepare_reciprocal_impacts(&mut b, &[pair], Default::default(), || false)
        .unwrap();
    let old = [a.state().clone(), b.state().clone()];
    assert!(matches!(
        a.apply_reciprocal_impacts(&mut b, &prepared, || false),
        Err(ShellFrictionalImpactError::Shell(MovingLoadError::Invalid(
            "shell event impulse is not representable in the accepted velocity"
        )))
    ));
    assert_eq!(a.state(), &old[0]);
    assert_eq!(b.state(), &old[1]);
}
