//! Native coupled shell collisions, source work, and atomic measured acceptance.
use super::*;
use fs_contact::multicontact::{MultiImpactDeterminacy, MultiImpactError};

fn site(point: ShellLoadPoint, restitution: f64) -> ShellNormalImpactSite {
    ShellNormalImpactSite {
        point,
        normal: [0.0, 0.0, 1.0],
        activation_tolerance_m: 0.0,
        restitution,
    }
}
fn planes(
    motion: &MovingShellMotion<'_>,
    sites: &[ShellNormalImpactSite],
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
fn supported() -> (MovingShellModel, ShellRigidPayload) {
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
    let mut point = point(0, [0.2, 0.3, 0.5]);
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
    (
        MovingShellModel::new_with_payloads(p, 2, &[body], || false).unwrap(),
        body,
    )
}
fn support() -> ShellSupportMotion {
    ShellSupportMotion {
        velocity_m_s: [[0.03, -0.02, -0.1]; 3],
        ..Default::default()
    }
}
fn driven(model: &MovingShellModel, body: ShellRigidPayload) -> MovingShellMotion<'_> {
    let zero = vec![0.0; model.assembly().mass.dimension()];
    model
        .start_driven(
            &zero,
            &zero,
            support(),
            &[wrench(body.point)],
            config(),
            || false,
        )
        .unwrap()
}

#[test]
fn simultaneous_offset_contacts_stop_an_initially_separating_point_and_keep_native_restitution() {
    let model = MovingShellModel::new(plate(), 2, || false).unwrap();
    let zero = vec![0.0; model.assembly().mass.dimension()];
    let mut velocity = zero.clone();
    velocity[2] = -0.4;
    velocity[4] = 12.0;
    let points = [0.05, -0.05].map(|arm| ShellLoadPoint {
        triangle: 0,
        barycentric: [1.0, 0.0, 0.0],
        offset_m: [arm, 0.0, 0.0],
    });
    // Bare node 0: m=1/12 kg, I=m*h^2/12=1/14400 kg m^2.
    // Two opposite 5 cm arms give W=[[48,-24],[-24,48]].
    // Initial normal speeds are [-1, +0.2]: the separating row must participate.
    for restitution in [0.0, 0.5, 1.0] {
        let mut motion = model
            .start(&zero, &velocity, &[], config(), || false)
            .unwrap();
        let sites = points.map(|point| site(point, restitution));
        let prepared = motion
            .prepare_normal_impacts(&sites, ShellNormalImpactControls::default(), || false)
            .unwrap();
        for (actual, expected) in prepared
            .inverse_effective_mass()
            .iter()
            .zip([48.0, -24.0, -24.0, 48.0])
        {
            near(*actual, expected, 3e-13);
        }
        let obstacles = planes(&motion, &sites, [0.0; 3]);
        let before = motion.state().clone();
        let result = motion.apply_normal_impacts(&prepared, &obstacles, || false);
        if restitution == 1.0 {
            // These local Newton targets add net energy in the assembled event.
            // The native owner refuses; no hidden restitution reduction occurs.
            assert!(matches!(
                result,
                Err(ShellNormalImpactError::Contact(
                    MultiImpactError::EnergyIncrease { .. }
                ))
            ));
            assert_eq!(motion.state(), &before);
            continue;
        }
        let report = result.unwrap();
        let p = [
            (48.0 * (1.0 + restitution) - 24.0 * 0.2) / 1728.0,
            (24.0 * (1.0 + restitution) - 48.0 * 0.2) / 1728.0,
        ];
        assert_eq!(report.contact.active_indices, [0, 1]);
        assert!(report.contact.normal_impulses[1] > 0.0);
        assert!(matches!(
            report.contact.determinacy,
            MultiImpactDeterminacy::UniqueImpulses
        ));
        for i in 0..2 {
            near(report.contact.normal_impulses[i], p[i], 3e-16);
            near(
                report.relative_normal_speeds_after_m_s[i],
                if i == 0 { restitution } else { 0.0 },
                4e-14,
            );
        }
        let final_translation = 0.5 * restitution;
        let final_rotation = -10.0 * restitution;
        near(motion.state().velocity[2], final_translation, 2e-14);
        near(motion.state().velocity[4], final_rotation, 3e-13);
        let expected_energy =
            0.5 * final_translation.powi(2) / 12.0 + 0.5 * final_rotation.powi(2) / 14400.0;
        let initial_energy = 7.0 / 600.0;
        near(motion.state().kinetic_energy_j, expected_energy, 2e-16);
        near(
            report.dissipated_energy_j,
            initial_energy - expected_energy,
            2e-16,
        );
        let impulse = report.impulse.unwrap();
        near(impulse.momentum_change.linear[2], p[0] + p[1], 2e-16);
        near(
            impulse.momentum_change.angular[0],
            -2.0 * (p[0] + p[1]),
            3e-16,
        );
        near(
            impulse.momentum_change.angular[1],
            -(p[0] + p[1]) + 0.05 * (p[1] - p[0]),
            3e-16,
        );
        assert_eq!(impulse.support_impulse, ShellSpatialResultant::default());
        near(report.work_balance_residual_j, 0.0, 2e-16);
        assert_eq!(motion.state().displacement, before.displacement);
        assert_eq!(motion.state().time_s, before.time_s);
        assert_eq!(motion.state().step, before.step);
    }
}

#[test]
fn off_center_tensor_payload_impact_keeps_plane_and_driven_pin_work_separate() {
    let (model, body) = supported();
    let mut motion = driven(&model, body);
    let sites = [site(body.point, 0.25)];
    let prepared = motion
        .prepare_normal_impacts(&sites, ShellNormalImpactControls::default(), || false)
        .unwrap();
    // Condense the bare nodal rotary inertias analytically to the COM rotation.
    // I_eff = I_COM + m*(r^2 Id-r r^T) + I_node/sum(N_i^2)*Id.
    let r = body.point.offset_m;
    let bare_rotary = (1.0 / 14400.0) / body.point.barycentric.iter().map(|n| n * n).sum::<f64>();
    let inertia: [[f64; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            body.inertia_com_kg_m2[i][j]
                + body.mass_kg * (if i == j { dot(r, r) } else { 0.0 } - r[i] * r[j])
                + if i == j { bare_rotary } else { 0.0 }
        })
    });
    let cofactors = [
        cross(inertia[1], inertia[2]),
        cross(inertia[2], inertia[0]),
        cross(inertia[0], inertia[1]),
    ];
    let determinant = dot(inertia[0], cofactors[0]);
    let torque = cross(r, [0.0, 0.0, 1.0]);
    let angular_response = std::array::from_fn(|axis| {
        (0..3).map(|i| torque[i] * cofactors[i][axis]).sum::<f64>() / determinant
    });
    let mobility = dot(torque, angular_response);
    near(prepared.inverse_effective_mass()[0], mobility, 2e-14);
    let obstacle_velocity = [0.01, 0.04, 0.2];
    let obstacles = planes(&motion, &sites, obstacle_velocity);
    let before = motion.state().clone();
    let before_com = support().velocity_m_s[0];
    let magnitude = 0.3 * 1.25 / mobility;
    let omega = angular_response.map(|value| magnitude * value);
    let delta_com = cross(omega, r);
    let after_com = std::array::from_fn(|i| before_com[i] + delta_com[i]);
    let report = motion
        .apply_normal_impacts(&prepared, &obstacles, || false)
        .unwrap();
    let sampled = motion.observe(body.point).unwrap();
    for axis in 0..3 {
        near(sampled.velocity_m_s[axis], after_com[axis], 3e-15);
        near(sampled.angular_velocity_rad_s[axis], omega[axis], 2e-14);
    }
    near(report.contact.normal_impulses[0], magnitude, 3e-15);
    near(report.relative_normal_speeds_after_m_s[0], 0.075, 3e-15);
    let reaction =
        std::array::from_fn(|i| body.mass_kg * delta_com[i] - if i == 2 { magnitude } else { 0.0 });
    let mid = [1.3, -1.5, 0.5];
    let reaction_moment = cross(mid, reaction);
    let impulse = report.impulse.unwrap();
    for axis in 0..3 {
        near(impulse.support_impulse.linear[axis], reaction[axis], 3e-15);
        near(
            impulse.support_impulse.angular[axis],
            reaction_moment[axis],
            5e-15,
        );
        near(impulse.momentum_residual.linear[axis], 0.0, 3e-15);
        near(impulse.momentum_residual.angular[axis], 0.0, 5e-15);
    }
    let plane_work = magnitude * obstacle_velocity[2];
    let pin_work = dot(reaction, before_com);
    near(report.obstacle_work_j, plane_work, 2e-16);
    near(impulse.support_work_j, pin_work, 2e-16);
    assert!(plane_work > 1e-4 && pin_work.abs() > 1e-4);

    // Independent bare shell + physical COM tensor energy inventory.
    let bare_translation_energy = 0.5 * 0.25 * dot(before_com, before_com);
    let mut after_energy = bare_translation_energy;
    for spin in motion.state().velocity.chunks_exact(3) {
        after_energy += 0.5 / 14400.0 * spin.iter().map(|v| v * v).sum::<f64>();
    }
    after_energy += 0.5 * body.mass_kg * dot(after_com, after_com)
        + 0.5 * dot(omega, body.inertia_com_kg_m2.map(|row| dot(row, omega)));
    let before_energy = bare_translation_energy + 0.5 * body.mass_kg * dot(before_com, before_com);
    let loss = 0.5 * magnitude * (0.3 - 0.075);
    near(motion.state().kinetic_energy_j, after_energy, 2e-16);
    near(
        impulse.kinetic_energy_change_j,
        after_energy - before_energy,
        2e-16,
    );
    near(
        after_energy - before_energy,
        plane_work + pin_work - loss,
        2e-16,
    );
    near(report.dissipated_energy_j, loss, 2e-16);
    near(report.work_balance_residual_j, 0.0, 2e-16);
    assert_eq!(motion.support_motion(), Some(&support()));
    assert_eq!(motion.state().displacement, before.displacement);
    let loads = [wrench(body.point)];
    let mut fresh = model
        .start_driven(
            &motion.state().displacement,
            &motion.state().velocity,
            support(),
            &loads,
            config(),
            || false,
        )
        .unwrap();
    assert_eq!(motion.state(), fresh.state());
    assert_eq!(motion.support_reactions(), fresh.support_reactions());
    assert_eq!(
        motion
            .advance_driven(config().dt_s, &loads, support().acceleration_m_s2, || false)
            .unwrap(),
        fresh
            .advance_driven(config().dt_s, &loads, support().acceleration_m_s2, || false)
            .unwrap()
    );
}

#[test]
fn invalid_geometry_ambiguous_sources_and_late_audits_do_not_publish_partial_collisions() {
    let (model, body) = supported();
    let mut motion = driven(&model, body);
    let sites = [site(body.point, 0.0)];
    let prepared = motion
        .prepare_normal_impacts(&sites, ShellNormalImpactControls::default(), || false)
        .unwrap();
    let obstacles = planes(&motion, &sites, [0.0, 0.0, 0.2]);
    let before = motion.state().clone();
    let reactions = motion.support_reactions().to_vec();
    let mut invalid = sites;
    invalid[0].normal = [0.0, 0.0, 2.0];
    assert!(
        motion
            .prepare_normal_impacts(&invalid, ShellNormalImpactControls::default(), || false)
            .is_err()
    );
    assert!(
        motion
            .prepare_normal_impacts(&[sites[0]; 3], ShellNormalImpactControls::default(), || {
                false
            })
            .is_err()
    );
    assert!(
        motion
            .apply_normal_impacts(&prepared, &[], || false)
            .is_err()
    );
    let mut penetrating = obstacles.clone();
    penetrating[0].point_m[2] += 0.001;
    assert!(matches!(
        motion.apply_normal_impacts(&prepared, &penetrating, || false),
        Err(ShellNormalImpactError::Contact(
            MultiImpactError::InvalidContact { .. }
        ))
    ));
    let (other_model, other_body) = supported();
    let mut other = driven(&other_model, other_body);
    let other_before = other.state().clone();
    assert!(
        other
            .apply_normal_impacts(&prepared, &obstacles, || false)
            .is_err()
    );
    assert_eq!(other.state(), &other_before);

    let duplicated = motion
        .prepare_normal_impacts(&[sites[0]; 2], ShellNormalImpactControls::default(), || {
            false
        })
        .unwrap();
    let mut different_sources = vec![obstacles[0]; 2];
    different_sources[1].velocity_m_s[2] += 0.01;
    assert!(matches!(
        motion.apply_normal_impacts(&duplicated, &different_sources, || false),
        Err(ShellNormalImpactError::Contact(
            MultiImpactError::UnresolvedWork { .. }
        ))
    ));
    let strict = motion
        .prepare_normal_impacts(
            &sites,
            ShellNormalImpactControls {
                speed_tolerance_m_s: 0.0,
                work_tolerance_j: 0.0,
                ..Default::default()
            },
            || false,
        )
        .unwrap();
    assert!(matches!(
        motion.apply_normal_impacts(&strict, &obstacles, || false),
        Err(ShellNormalImpactError::Audit { .. })
    ));
    assert_eq!(motion.state(), &before);
    assert_eq!(motion.support_reactions(), reactions);
    assert_eq!(motion.support_motion(), Some(&support()));

    let mut control = driven(&model, body);
    let mut polls = 0;
    let expected = control
        .apply_normal_impacts(&prepared, &obstacles, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        motion.apply_normal_impacts(&prepared, &obstacles, || {
            calls += 1;
            calls == polls
        }),
        Err(ShellNormalImpactError::Cancelled)
    ));
    assert_eq!(motion.state(), &before);
    assert_eq!(motion.support_reactions(), reactions);
    assert_eq!(motion.support_motion(), Some(&support()));
    assert_eq!(
        motion
            .apply_normal_impacts(&prepared, &obstacles, || false)
            .unwrap(),
        expected
    );
    assert_eq!(motion.state(), control.state());
    assert_eq!(motion.support_reactions(), control.support_reactions());

    let mut redundant_motion = driven(&model, body);
    let redundant = redundant_motion
        .apply_normal_impacts(&duplicated, &[obstacles[0]; 2], || false)
        .unwrap();
    assert!(matches!(
        redundant.contact.determinacy,
        MultiImpactDeterminacy::UniqueVelocity(_)
    ));
    near(
        redundant.contact.normal_impulses.iter().sum::<f64>(),
        expected.contact.normal_impulses[0],
        1e-12,
    );
    near(
        redundant.dissipated_energy_j,
        expected.dissipated_energy_j,
        1e-12,
    );
    for (a, b) in redundant_motion
        .state()
        .velocity
        .iter()
        .zip(&control.state().velocity)
    {
        near(*a, *b, 1e-11);
    }
    let accepted = motion.state().clone();
    let resting = motion
        .apply_normal_impacts(&prepared, &obstacles, || false)
        .unwrap();
    assert!(resting.impulse.is_none());
    assert_eq!(motion.state(), &accepted);
    let empty = motion
        .prepare_normal_impacts(&[], ShellNormalImpactControls::default(), || false)
        .unwrap();
    assert!(
        motion
            .apply_normal_impacts(&empty, &[], || false)
            .unwrap()
            .impulse
            .is_none()
    );
    assert_eq!(motion.state(), &accepted);
}
