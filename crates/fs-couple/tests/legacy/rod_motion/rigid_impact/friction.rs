//! G0/G1/G3/G4: real Coulomb impacts, native spatial inertia and continuing motion.

use super::*;
use fs_contact::frictional_multicontact::FrictionalContactConfig;
use fs_couple::rod_motion::impulse::RodPointImpulse;
use fs_couple::rod_motion::rigid_impact::friction::{
    RodRigidFrictionalContact, RodRigidFrictionalControls, RodRigidFrictionalError,
    RodRigidFrictionalReport,
};

fn bounds() -> RodRigidFrictionalControls {
    let mut impact = impact_controls();
    impact.velocity_tolerance_m_s = 2e-9;
    impact.impulse.work_tolerance_j = 2e-9;
    RodRigidFrictionalControls {
        impact,
        contact: FrictionalContactConfig {
            velocity_tolerance: 2e-11,
            impulse_tolerance: 2e-11,
            complementarity_tolerance: 2e-11,
            energy_tolerance: 2e-11,
            ..FrictionalContactConfig::default()
        },
    }
}
fn contact(mu: f64) -> RodRigidFrictionalContact {
    RodRigidFrictionalContact {
        geometry: witness(1, 0.35),
        tangent_world: [0.0, 1.0, 0.0],
        friction_coefficient: mu,
    }
}
fn audit(report: &RodRigidFrictionalReport, mu: f64) {
    let p = report.contact.impulses[0];
    let v = report.relative_velocity_after_m_s;
    assert!(p[0] >= -2e-11 && v[0] >= -2e-9);
    assert!(p[1].hypot(p[2]) <= mu * p[0] + 2e-9);
    let slip = v[1].hypot(v[2]);
    if slip > 2e-9 {
        near(p[1], -mu * p[0] * v[1] / slip, 2e-9);
        near(p[2], -mu * p[0] * v[2] / slip, 2e-9);
    }
    vector_near(report.endpoint_velocity_residual_m_s, [0.0; 3], 2e-9);
    vector_near(report.linear_momentum_residual_n_s, [0.0; 3], 2e-11);
    vector_near(report.angular_momentum_residual_nm_s, [0.0; 3], 2e-11);
    near(report.energy_residual_j, 0.0, 2e-9);
    near(report.work_energy_residual_j, 0.0, 2e-9);
    assert!(report.contact.dissipated_energy >= 0.0);
}

#[test]
fn glancing_finite_masses_stick_or_slide_and_continue_on_both_native_owners() {
    let gate = CancelGate::new();
    let rod = free_rod([-0.6, 0.8, -0.6]);
    let mass = inertia(2.3);
    let rigid = body(mass, [0.8, 0.0, 0.0], [1.4, 0.0, 0.0]);
    let rod_mass = rod.model().inertia()[1].mass_kg();
    let inverse_mass = 1.0 / rod_mass + 1.0 / mass.mass();
    let normal = 2.0 * 1.35 / inverse_mass;
    for mu in [0.0, 0.1, 0.7] {
        let (next, finite, report) = rod
            .apply_rigid_frictional_impact(rigid, mass, contact(mu), bounds(), &gate)
            .unwrap();
        let tangent = (1.0 / inverse_mass).min(mu * normal);
        let expected = [normal, -0.8 * tangent, 0.6 * tangent];
        vector_near(report.impulse_world_n_s, expected, 2e-9);
        vector_near(
            next.velocities_m_s().unwrap()[1],
            add([-0.6, 0.8, -0.6], expected.map(|x| x / rod_mass)),
            2e-9,
        );
        vector_near(
            array(finite.center_of_mass_velocity_world(mass).unwrap()),
            sub([1.4, 0.0, 0.0], expected.map(|x| x / mass.mass())),
            2e-9,
        );
        let loss = 0.5 * (1.0 - 0.35_f64.powi(2)) * 4.0 / inverse_mass + tangent
            - 0.5 * inverse_mass * tangent * tangent;
        near(report.contact.dissipated_energy, loss, 2e-9);
        near(report.kinetic_energy_change_j, -loss, 2e-9);
        if mu == 0.7 {
            near(report.relative_velocity_after_m_s[1], 0.0, 2e-9);
            near(report.relative_velocity_after_m_s[2], 0.0, 2e-9);
        } else {
            assert!(report.relative_velocity_after_m_s[1] > 0.0);
        }
        audit(&report, mu);
        assert_eq!(next.pose(), rod.pose());
        assert_eq!(finite.pose(), rigid.pose());
        assert_eq!(next.time_s(), rod.time_s());
        assert_eq!(next.accepted_steps(), rod.accepted_steps());
        assert!(std::ptr::eq(next.model(), rod.model()));
        assert_eq!(
            next.linear_momenta_kg_m_s()[0],
            rod.linear_momenta_kg_m_s()[0]
        );
        assert_eq!(next.body_spins_kg_m2_s(), rod.body_spins_kg_m2_s());
        let continued = next
            .step_candidate(1e-6, &[RodNodalWrench::default(); 2], controls(1e-6), &gate)
            .unwrap()
            .0;
        assert_eq!(continued.accepted_steps(), rod.accepted_steps() + 1);
        let body_step = fs_mbd::RigidBodyIntegrator::new(fs_mbd::Gravity::ZERO)
            .step(finite, mass, fs_mbd::Wrench::ZERO, 1e-6)
            .unwrap();
        assert_eq!(
            body_step.state_after.linear_momentum_world(),
            finite.linear_momentum_world()
        );
        assert_ne!(
            body_step.state_after.pose().position_world(),
            finite.pose().position_world()
        );
        if mu == 0.0 {
            let original = rod
                .apply_rigid_normal_impact(rigid, mass, witness(1, 0.35), impact_controls(), &gate)
                .unwrap();
            vector_near(
                next.linear_momenta_kg_m_s()[1],
                original.0.linear_momenta_kg_m_s()[1],
                2e-9,
            );
            vector_near(
                array(finite.linear_momentum_world()),
                array(original.1.linear_momentum_world()),
                2e-9,
            );
            near(
                report.contact.dissipated_energy,
                original.2.contact.dissipated_energy,
                2e-9,
            );
        }
    }
}

#[test]
fn eccentric_off_centre_impulses_match_independent_spatial_mass_and_rotate_covariantly() {
    let gate = CancelGate::new();
    let model = Arc::new(
        eccentric_source()
            .into_spatial_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let original = model.inertia()[1];
    let center = original.center_of_mass_body_m();
    let tensor = original.inertia_com_body_kg_m2();
    let spatial = fs_mbd::articulated::SpatialInertia::new(
        original.mass_kg(),
        Vec3::new(center[0], center[1], center[2]),
        Mat3 {
            m: std::array::from_fn(|c| tensor[c / 3][c % 3]),
        },
    )
    .unwrap()
    .matrix()
    .unwrap();
    let rod_factor = fs_la::factor::lu(&spatial.m, 6).unwrap();
    let body_tensor = Mat3 {
        m: [0.07, 0.01, -0.004, 0.01, 0.09, 0.006, -0.004, 0.006, 0.08],
    };
    let body_factor = fs_la::factor::lu(&body_tensor.m, 3).unwrap();
    let mass = MassProperties::from_inertia_tensor(2.4, BodyVec::ZERO, body_tensor).unwrap();
    let rod_arm = [0.03, -0.06, 0.04];
    let body_arm = [0.1, -0.08, 0.06];
    // Independent six-by-six origin-coordinate solve, including the spatial
    // translational/rotational blocks, rather than a central-inertia formula.
    let mut reference_rod = [[0.0; 3]; 3];
    let mut reference_body = [[0.0; 3]; 3];
    for j in 0..3 {
        let axis = std::array::from_fn(|c| if c == j { 1.0 } else { 0.0 });
        let torque = cross(rod_arm, axis);
        let mut response = [torque[0], torque[1], torque[2], axis[0], axis[1], axis[2]];
        rod_factor.solve(&mut response);
        let velocity = add(
            [response[3], response[4], response[5]],
            cross([response[0], response[1], response[2]], rod_arm),
        );
        let mut omega = cross(body_arm, axis);
        body_factor.solve(&mut omega);
        let body_velocity = add(axis.map(|x| x / mass.mass()), cross(omega, body_arm));
        for i in 0..3 {
            reference_rod[i][j] = velocity[i];
            reference_body[i][j] = body_velocity[i];
        }
    }
    assert!((reference_rod[0][1] + reference_body[0][1]).abs() > 1e-3);
    let rod_v = [-0.2, 1.1, -0.4];
    let rod_w = [0.2, -0.1, 0.3];
    let body_v = [1.0, -0.2, 0.15];
    let body_w = [-0.3, 0.1, 0.2];
    let q = So3::exp(So3Tangent::new(Vec3::new(0.4, -0.2, 0.3))).unwrap();
    let mut reference_impulse = None;
    let mut reference_loss = None;
    for (q, shift, boost, hint) in [
        (So3::identity(), [0.0; 3], [0.0; 3], [0.0, 1.0, 0.0]),
        (q, [0.7, -0.2, 0.1], [0.4, -0.3, 0.2], [0.0, 1.0, 0.0]),
        (q, [0.7, -0.2, 0.1], [0.4, -0.3, 0.2], [0.0, 0.6, 0.8]),
    ] {
        let mut pose = model.initial_pose();
        for node in 0..2 {
            pose.positions_m[node] = add(rotate(q, pose.positions_m[node]), shift);
            pose.frames[node] = q;
        }
        let rod = SpatialRodMotion::new_spatial(
            model.clone(),
            pose,
            &[add(rotate(q, rod_v), boost); 2],
            &[rod_w; 2],
            0.0,
            1e-12,
            &gate,
        )
        .unwrap();
        // Positive roundoff-sized gap avoids a sign ambiguity after a rotated
        // coordinate transform. The actual separated-point moment is audited.
        let position = add(
            rotate(
                q,
                sub(
                    add([0.8, 0.0, 0.0], rod_arm),
                    add(body_arm, [1e-13, 0.0, 0.0]),
                ),
            ),
            shift,
        );
        let rigid = RigidBodyState::new(
            Pose::new(vector(position), orientation(q)).unwrap(),
            vector(add(rotate(q, body_v), boost)).scale(mass.mass()),
            mass.angular_momentum_body_checked(vector(body_w)).unwrap(),
        )
        .unwrap();
        let witness = RodRigidFrictionalContact {
            geometry: RodRigidNormalContact {
                point_body_m: rod_arm,
                rigid_arm_body_m: body_arm,
                normal_world: rotate(q, [1.0, 0.0, 0.0]),
                ..witness(1, 0.35)
            },
            tangent_world: rotate(q, hint),
            friction_coefficient: 0.15,
        };
        let (next, finite, report) = rod
            .apply_rigid_frictional_impact(rigid, mass, witness, bounds(), &gate)
            .unwrap();
        for i in 0..3 {
            let left = rotate(q.inverse(), report.frame_world[i]);
            for j in 0..3 {
                let right = rotate(q.inverse(), report.frame_world[j]);
                let quadratic = |matrix: [[f64; 3]; 3]| {
                    (0..3)
                        .map(|a| {
                            (0..3)
                                .map(|b| left[a] * matrix[a][b] * right[b])
                                .sum::<f64>()
                        })
                        .sum::<f64>()
                };
                near(
                    report.rod_inverse_mass_kg_inverse[i][j],
                    quadratic(reference_rod),
                    2e-11,
                );
                near(
                    report.rigid_inverse_mass_kg_inverse[i][j],
                    quadratic(reference_body),
                    2e-11,
                );
            }
        }
        let p = rotate(q.inverse(), report.impulse_world_n_s);
        let torque = cross(rod_arm, p);
        let mut response = [torque[0], torque[1], torque[2], p[0], p[1], p[2]];
        rod_factor.solve(&mut response);
        vector_near(
            next.angular_velocities_body_rad_s().unwrap()[1],
            add(rod_w, [response[0], response[1], response[2]]),
            2e-9,
        );
        vector_near(
            next.velocities_m_s().unwrap()[1],
            add(
                rotate(q, add(rod_v, [response[3], response[4], response[5]])),
                boost,
            ),
            2e-9,
        );
        let mut angular = cross(body_arm, p);
        body_factor.solve(&mut angular);
        vector_near(
            array(
                mass.angular_velocity_body_checked(finite.angular_momentum_body())
                    .unwrap(),
            ),
            sub(body_w, angular),
            2e-9,
        );
        assert!(angular.iter().any(|x| x.abs() > 1e-3));
        assert!(report.contact.impulses[0][1].hypot(report.contact.impulses[0][2]) > 1e-3);
        let normal_only = rod
            .apply_rigid_normal_impact(rigid, mass, witness.geometry, impact_controls(), &gate)
            .unwrap();
        let extra_spin = finite
            .angular_momentum_body()
            .sub(normal_only.1.angular_momentum_body());
        assert!(
            extra_spin.norm_squared() > 1e-8,
            "glancing friction must transfer additional spin"
        );
        if let Some(previous) = reference_impulse {
            vector_near(p, previous, 2e-9);
        }
        if let Some(previous) = reference_loss {
            near(report.contact.dissipated_energy, previous, 2e-9);
        }
        reference_impulse = Some(p);
        reference_loss = Some(report.contact.dissipated_energy);
        audit(&report, 0.15);
    }
}

#[test]
fn stationary_supports_return_the_missing_linear_and_spin_reactions() {
    let gate = CancelGate::new();
    for lock_frame in [false, true] {
        let source = SpatialMaterialRod::try_straight(0.8, vec![section(1.0)])
            .unwrap()
            .with_supports(
                &[RodSupport {
                    node: 1,
                    position_m: [Some(0.8), None, None],
                    frame: lock_frame.then_some(So3::identity()),
                }],
                &gate,
            )
            .unwrap();
        let model = Arc::new(
            source
                .into_motion_model(RodRootBoundary::Free, &gate)
                .unwrap(),
        );
        let rod = RodMotion::new(
            model.clone(),
            model.initial_pose(),
            &[[0.0; 3]; 2],
            &[[0.0; 3]; 2],
            0.0,
            &gate,
        )
        .unwrap();
        let mass = inertia(2.0);
        let rigid = body(mass, [0.8, 0.07, 0.0], [1.0, -1.0, 0.0]);
        let witness = RodRigidFrictionalContact {
            geometry: RodRigidNormalContact {
                point_body_m: [0.0, 0.07, 0.0],
                ..witness(1, 0.5)
            },
            ..contact(0.2)
        };
        let (next, _, report) = rod
            .apply_rigid_frictional_impact(rigid, mass, witness, bounds(), &gate)
            .unwrap();
        let rotational = if lock_frame {
            0.0
        } else {
            0.07 * 0.07 / model.inertia()[1].principal_kg_m2()[2]
        };
        let jn = 1.5 / (rotational + 1.0 / mass.mass());
        let jt = -(0.2 * jn).min(1.0 / (1.0 / model.inertia()[1].mass_kg() + 1.0 / mass.mass()));
        vector_near(report.impulse_world_n_s, [jn, jt, 0.0], 2e-9);
        near(next.linear_momenta_kg_m_s()[1][0], 0.0, 2e-11);
        near(next.linear_momenta_kg_m_s()[1][1], jt, 2e-9);
        near(
            next.body_spins_kg_m2_s()[1][2],
            if lock_frame { 0.0 } else { -0.07 * jn },
            2e-9,
        );
        let support = &report.rod_impulse.as_ref().unwrap().support_impulses[1];
        vector_near(support.linear_n_s, [-jn, 0.0, 0.0], 2e-9);
        vector_near(
            support.angular_nm_s,
            [0.0, 0.0, if lock_frame { 0.07 * jn } else { 0.0 }],
            2e-9,
        );
        audit(&report, 0.2);
    }
}

#[test]
fn inactive_cancelled_and_late_refused_friction_events_leave_both_states_unchanged() {
    let gate = CancelGate::new();
    let rod = free_rod([-0.6, 0.8, -0.6]);
    let mass = inertia(2.3);
    let rigid = body(mass, [0.8, 0.0, 0.0], [1.4, 0.0, 0.0]);
    let original = rod.observe([0.0; 3]).unwrap();
    for other in [
        body(mass, [0.79, 0.0, 0.0], [1.4, 0.0, 0.0]),
        body(mass, [0.8, 0.0, 0.0], [-1.4, 0.0, 0.0]),
    ] {
        let (next, finite, report) = rod
            .apply_rigid_frictional_impact(other, mass, contact(0.1), bounds(), &gate)
            .unwrap();
        assert_eq!(report.impulse_world_n_s, [0.0; 3]);
        assert!(report.rod_impulse.is_none() && report.rigid_impulse.is_none());
        assert_eq!(next.linear_momenta_kg_m_s(), rod.linear_momenta_kg_m_s());
        assert_eq!(next.body_spins_kg_m2_s(), rod.body_spins_kg_m2_s());
        assert_eq!(next.pose(), rod.pose());
        assert_eq!(finite, other);
    }
    let bad_axis = RodRigidFrictionalContact {
        tangent_world: [1.0, 0.0, 0.0],
        ..contact(0.1)
    };
    assert!(
        rod.apply_rigid_frictional_impact(rigid, mass, bad_axis, bounds(), &gate)
            .is_err()
    );
    assert!(matches!(
        rod.apply_rigid_frictional_impact(rigid, mass, contact(-0.1), bounds(), &gate),
        Err(RodRigidFrictionalError::Contact(_))
    ));
    let separated = body(mass, [0.8 - 1e-9, 0.0, 0.0], [1.4, 0.0, 0.0]);
    assert!(matches!(
        rod.apply_rigid_frictional_impact(separated, mass, contact(0.1), bounds(), &gate),
        Err(RodRigidFrictionalError::Impact(
            RodRigidImpactError::Invalid(_)
        ))
    ));
    let crooked = body(mass, [0.8, 1e-4, 0.0], [1.4, 0.0, 0.0]);
    let mut loose_geometry = bounds();
    loose_geometry.impact.geometry_tolerance_m = 1e-3;
    // The law and both native impulses can succeed, but their actual arm
    // separation creates an inadmissible couple. No hidden midpoint repairs it.
    assert!(matches!(
        rod.apply_rigid_frictional_impact(crooked, mass, contact(0.1), loose_geometry, &gate),
        Err(RodRigidFrictionalError::Impact(
            RodRigidImpactError::Balance { .. }
        ))
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        rod.apply_rigid_frictional_impact(rigid, mass, contact(0.1), bounds(), &cancelled)
            .is_err()
    );
    let mut exhausted = bounds();
    exhausted.contact.maximum_evaluations = 0;
    assert!(
        rod.apply_rigid_frictional_impact(rigid, mass, contact(0.1), exhausted, &gate)
            .is_err()
    );
    assert_eq!(
        rod.observe([0.0; 3]).unwrap().kinetic_energy_j,
        original.kinetic_energy_j
    );
    let first = rod
        .apply_rigid_frictional_impact(rigid, mass, contact(0.1), bounds(), &gate)
        .unwrap();
    let retry = rod
        .apply_rigid_frictional_impact(rigid, mass, contact(0.1), bounds(), &gate)
        .unwrap();
    assert_eq!(
        first.0.linear_momenta_kg_m_s(),
        retry.0.linear_momenta_kg_m_s()
    );
    assert_eq!(first.1, retry.1);
}

fn friction_controls() -> RodRigidFrictionalControls {
    RodRigidFrictionalControls {
        impact: impact_controls(),
        contact: FrictionalContactConfig {
            maximum_contacts: 1,
            velocity_tolerance: 1e-12,
            impulse_tolerance: 1e-12,
            complementarity_tolerance: 1e-12,
            energy_tolerance: 1e-12,
            ..Default::default()
        },
    }
}
fn friction_witness(mu: f64, restitution: f64) -> RodRigidFrictionalContact {
    RodRigidFrictionalContact {
        geometry: witness(1, restitution),
        tangent_world: [0.0, 1.0, 0.0],
        friction_coefficient: mu,
    }
}
fn matrix_vector(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| (0..3).map(|j| m[i][j] * v[j]).sum())
}

#[test]
fn full_spatial_mass_cross_terms_match_independent_dense_solve_in_rotated_frames() {
    let gate = CancelGate::new();
    let model = Arc::new(
        eccentric_source()
            .into_spatial_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let node = model.inertia()[1];
    let center = node.center_of_mass_body_m();
    let tensor = node.inertia_com_body_kg_m2();
    let rod_mass = fs_mbd::articulated::SpatialInertia::new(
        node.mass_kg(),
        Vec3::new(center[0], center[1], center[2]),
        Mat3 {
            m: std::array::from_fn(|c| tensor[c / 3][c % 3]),
        },
    )
    .unwrap()
    .matrix()
    .unwrap();
    let rod_factor = fs_la::factor::lu(&rod_mass.m, 6).unwrap();
    let body_tensor = Mat3 {
        m: [0.07, 0.01, -0.004, 0.01, 0.09, 0.006, -0.004, 0.006, 0.08],
    };
    let mass = MassProperties::from_inertia_tensor(2.4, BodyVec::ZERO, body_tensor).unwrap();
    let body_factor = fs_la::factor::lu(&body_tensor.m, 3).unwrap();
    let rod_arm = [0.03, -0.06, 0.04];
    let body_arm = [0.1, -0.08, 0.06];
    let mut response = [[0.0; 6]; 3];
    let mut mobility = [[0.0; 3]; 3];
    for j in 0..3 {
        let direction = std::array::from_fn(|c| f64::from(c == j));
        let torque = cross(rod_arm, direction);
        response[j] = [
            torque[0],
            torque[1],
            torque[2],
            direction[0],
            direction[1],
            direction[2],
        ];
        rod_factor.solve(&mut response[j]);
        let rod_flow = add(
            [response[j][3], response[j][4], response[j][5]],
            cross([response[j][0], response[j][1], response[j][2]], rod_arm),
        );
        let mut body_omega = cross(body_arm, direction);
        body_factor.solve(&mut body_omega);
        let body_flow = add(
            direction.map(|x| x / mass.mass()),
            cross(body_omega, body_arm),
        );
        for i in 0..3 {
            mobility[i][j] = rod_flow[i] + body_flow[i];
        }
    }
    assert!(
        mobility[0][1].abs() > 1e-3,
        "fixture must exercise coupled directions"
    );
    let impulse = [1.0, -0.03, 0.02];
    let relative = matrix_vector(mobility, impulse).map(|x| -x);
    assert!(relative[0] < 0.0);
    let rod_omega = [0.2, -0.1, 0.3];
    let body_omega = [-0.3, 0.1, 0.2];
    let body_velocity = [0.1, 0.2, -0.05];
    let rod_velocity = sub(
        add(relative, add(body_velocity, cross(body_omega, body_arm))),
        cross(rod_omega, rod_arm),
    );
    let jump: [f64; 6] = std::array::from_fn(|i| (0..3).map(|j| response[j][i] * impulse[j]).sum());
    for (q, translation, boost) in [
        (So3::identity(), [0.0; 3], [0.0; 3]),
        (
            So3::exp(So3Tangent::new(Vec3::new(0.4, -0.2, 0.3))).unwrap(),
            [0.7, -0.2, 0.1],
            [0.4, -0.3, 0.2],
        ),
    ] {
        for angle in [0.0_f64, 0.7] {
            let basis = [
                [1.0, 0.0, 0.0],
                [0.0, angle.cos(), angle.sin()],
                [0.0, -angle.sin(), angle.cos()],
            ];
            let mut pose = model.initial_pose();
            for i in 0..2 {
                pose.positions_m[i] = add(rotate(q, pose.positions_m[i]), translation);
                pose.frames[i] = q;
            }
            let rod = SpatialRodMotion::new_spatial(
                model.clone(),
                pose,
                &[add(rotate(q, rod_velocity), boost); 2],
                &[rod_omega; 2],
                0.0,
                1e-12,
                &gate,
            )
            .unwrap();
            let position = add(
                rotate(q, sub(add([0.8, 0.0, 0.0], rod_arm), body_arm)),
                translation,
            );
            let rigid = RigidBodyState::new(
                Pose::new(vector(position), orientation(q)).unwrap(),
                vector(add(rotate(q, body_velocity), boost)).scale(mass.mass()),
                mass.angular_momentum_body_checked(vector(body_omega))
                    .unwrap(),
            )
            .unwrap();
            let contact = RodRigidFrictionalContact {
                geometry: RodRigidNormalContact {
                    point_body_m: rod_arm,
                    rigid_arm_body_m: body_arm,
                    normal_world: rotate(q, basis[0]),
                    ..witness(1, 0.0)
                },
                tangent_world: rotate(q, basis[1]),
                friction_coefficient: 0.15,
            };
            let (next, finite, report) = rod
                .apply_rigid_frictional_impact(rigid, mass, contact, friction_controls(), &gate)
                .unwrap();
            for i in 0..3 {
                for j in 0..3 {
                    let expected: f64 = (0..3)
                        .flat_map(|a| {
                            (0..3).map(move |b| basis[i][a] * mobility[a][b] * basis[j][b])
                        })
                        .sum();
                    near(
                        report.rod_inverse_mass_kg_inverse[i][j]
                            + report.rigid_inverse_mass_kg_inverse[i][j],
                        expected,
                        3e-11,
                    );
                }
            }
            let world_impulse = report
                .frame_world
                .iter()
                .zip(report.contact.impulses[0])
                .fold([0.0; 3], |s, (axis, j)| add(s, axis.map(|x| x * j)));
            vector_near(world_impulse, rotate(q, impulse), 4e-11);
            vector_near(report.relative_velocity_after_m_s, [0.0; 3], 4e-11);
            vector_near(
                next.velocities_m_s().unwrap()[1],
                add(
                    rotate(q, add(rod_velocity, [jump[3], jump[4], jump[5]])),
                    boost,
                ),
                4e-11,
            );
            vector_near(
                next.angular_velocities_body_rad_s().unwrap()[1],
                add(rod_omega, [jump[0], jump[1], jump[2]]),
                4e-11,
            );
            vector_near(
                array(finite.linear_momentum_world()),
                sub(array(rigid.linear_momentum_world()), rotate(q, impulse)),
                4e-11,
            );
            vector_near(
                array(finite.angular_momentum_body()),
                sub(
                    array(rigid.angular_momentum_body()),
                    cross(body_arm, impulse),
                ),
                4e-11,
            );
            near(report.energy_residual_j, 0.0, 4e-11);
        }
    }
}

#[test]
fn small_frictional_loss_survives_large_common_translation() {
    let gate = CancelGate::new();
    let common = 1_073_741_824.0;
    let base = free_rod([common, common, 0.0]);
    let m = base.model().inertia()[1].mass_kg();
    let mut bounds = friction_controls();
    bounds.impact.impulse.work_tolerance_j = 1e-6;
    bounds.impact.velocity_tolerance_m_s = 1e-6;
    let (rod, _) = base
        .apply_point_impulses(
            &[RodPointImpulse {
                node: 1,
                point_body_m: [0.0; 3],
                linear_world_n_s: [-1.0, 1.0, 0.0],
                couple_world_nm_s: [0.0; 3],
            }],
            bounds.impact.impulse,
            &gate,
        )
        .unwrap();
    let mass = inertia(m);
    let rigid = RigidBodyState::new(
        Pose::new(BodyVec::new(0.8, 0.0, 0.0), UnitQuaternion::IDENTITY).unwrap(),
        BodyVec::new(m * common + 1.0, m * common - 1.0, 0.0),
        BodyVec::ZERO,
    )
    .unwrap();
    let (next, finite, report) = rod
        .apply_rigid_frictional_impact(rigid, mass, friction_witness(0.25, 0.0), bounds, &gate)
        .unwrap();
    near(report.contact.dissipated_energy, 1.4375 / m, 2e-10);
    near(report.kinetic_energy_change_j, -1.4375 / m, 2e-10);
    near(report.work_energy_residual_j, 0.0, 2e-10);
    let before =
        rod.observe([0.0; 3]).unwrap().kinetic_energy_j + rigid.kinetic_energy(mass).unwrap();
    let after =
        next.observe([0.0; 3]).unwrap().kinetic_energy_j + finite.kinetic_energy(mass).unwrap();
    assert!(((after - before) - report.kinetic_energy_change_j).abs() > 0.1);
}

#[test]
fn rotational_approach_cannot_be_hidden_as_noop_by_common_translation() {
    let gate = CancelGate::new();
    let base = free_rod([0.0; 3]);
    let model = Arc::new(base.model().clone());
    let rod = RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[[1e20, 0.0, 0.0]; 2],
        &[[0.0, 0.0, 1.0]; 2],
        0.0,
        &gate,
    )
    .unwrap();
    let mass = inertia(model.inertia()[1].mass_kg());
    let rigid = body(mass, [0.8, 0.1, 0.0], [1e20, 0.0, 0.0]);
    let mut contact = friction_witness(0.2, 0.4);
    contact.geometry.point_body_m = [0.0, 0.1, 0.0];
    assert!(matches!(
        rod.apply_rigid_frictional_impact(rigid, mass, contact, friction_controls(), &gate),
        Err(RodRigidFrictionalError::Impact(RodRigidImpactError::Rod(
            RodMotionError::ImpulseBalance { .. }
        )))
    ));
}
