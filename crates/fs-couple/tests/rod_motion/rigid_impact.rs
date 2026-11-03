//! G0/G1/G3/G4: actual momentum exchange with a finite rigid participant.

#[path = "rigid_impact/friction.rs"]
mod friction;

#[path = "rigid_impact/multiple.rs"]
mod multiple;

use super::*;
use fs_couple::rod_motion::SpatialRodMotion;
use fs_couple::rod_motion::impulse::RodImpulseControls;
use fs_couple::rod_motion::rigid_impact::{
    RodRigidImpactControls, RodRigidImpactError, RodRigidNormalContact,
};
use fs_ga::Mat3;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3 as BodyVec};

fn impact_controls() -> RodRigidImpactControls {
    RodRigidImpactControls {
        impulse: RodImpulseControls {
            balance_origin_m: [-0.3, 0.7, -0.5],
            linear_tolerance_n_s: 2e-11,
            angular_tolerance_nm_s: 2e-11,
            work_tolerance_j: 2e-11,
        },
        geometry_tolerance_m: 1e-11,
        velocity_tolerance_m_s: 2e-11,
    }
}
fn witness(node: usize, restitution: f64) -> RodRigidNormalContact {
    RodRigidNormalContact {
        node,
        point_body_m: [0.0; 3],
        rigid_arm_body_m: [0.0; 3],
        normal_world: [1.0, 0.0, 0.0],
        activation_tolerance_m: 1e-8,
        restitution,
    }
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|c| a[c] + b[c])
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|c| a[c] - b[c])
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn rotate(q: So3, x: [f64; 3]) -> [f64; 3] {
    let y = q.rotate(Vec3::new(x[0], x[1], x[2])).unwrap();
    [y.x, y.y, y.z]
}
fn vector(x: [f64; 3]) -> BodyVec {
    BodyVec::new(x[0], x[1], x[2])
}
fn array(x: BodyVec) -> [f64; 3] {
    [x.x, x.y, x.z]
}
fn vector_near(a: [f64; 3], b: [f64; 3], tol: f64) {
    for c in 0..3 {
        near(a[c], b[c], tol);
    }
}
fn orientation(q: So3) -> UnitQuaternion {
    let q = q.as_quat();
    UnitQuaternion::new(q.w, q.x, q.y, q.z).unwrap()
}
fn inertia(mass: f64) -> MassProperties {
    MassProperties::new(mass, BodyVec::ZERO, BodyVec::new(0.06, 0.08, 0.1)).unwrap()
}
fn body(mass: MassProperties, position: [f64; 3], velocity: [f64; 3]) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(vector(position), UnitQuaternion::IDENTITY).unwrap(),
        vector(velocity).scale(mass.mass()),
        BodyVec::ZERO,
    )
    .unwrap()
}
fn free_rod(velocity: [f64; 3]) -> RodMotion {
    let gate = CancelGate::new();
    let model = Arc::new(
        SpatialMaterialRod::try_straight(0.8, vec![section(1.0)])
            .unwrap()
            .into_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[velocity; 2],
        &[[0.0; 3]; 2],
        2.0,
        &gate,
    )
    .unwrap()
}

#[test]
fn collinear_finite_masses_match_independent_collision_and_continue() {
    let gate = CancelGate::new();
    let rod = free_rod([-0.6, 0.0, 0.0]);
    let mass = inertia(2.3);
    let rigid = body(mass, [0.8, 0.0, 0.0], [1.7, 0.0, 0.0]);
    let rod_mass = rod.model().inertia()[1].mass_kg();
    for restitution in [0.0, 0.35, 1.0] {
        let (next, finite, report) = rod
            .apply_rigid_normal_impact(
                rigid,
                mass,
                witness(1, restitution),
                impact_controls(),
                &gate,
            )
            .unwrap();
        let reduced_mass = rod_mass * mass.mass() / (rod_mass + mass.mass());
        let impulse = reduced_mass * (1.0 + restitution) * 2.3;
        near(report.contact.normal_impulse, impulse, 2e-12);
        near(
            next.velocities_m_s().unwrap()[1][0],
            -0.6 + impulse / rod_mass,
            2e-12,
        );
        near(
            finite.center_of_mass_velocity_world(mass).unwrap().x,
            1.7 - impulse / mass.mass(),
            2e-12,
        );
        let loss = 0.5 * reduced_mass * (1.0 - restitution * restitution) * 2.3 * 2.3;
        near(report.contact.dissipated_energy, loss, 2e-12);
        let initial_energy =
            rod.observe([0.0; 3]).unwrap().kinetic_energy_j + rigid.kinetic_energy(mass).unwrap();
        let final_energy =
            next.observe([0.0; 3]).unwrap().kinetic_energy_j + finite.kinetic_energy(mass).unwrap();
        near(initial_energy - final_energy, loss, 3e-12);
        vector_near(report.linear_momentum_residual_n_s, [0.0; 3], 2e-12);
        vector_near(report.angular_momentum_residual_nm_s, [0.0; 3], 2e-12);
        near(report.energy_residual_j, 0.0, 3e-12);
        assert_eq!(next.pose(), rod.pose());
        assert_eq!(finite.pose(), rigid.pose());
        assert_eq!(next.time_s(), rod.time_s());
        assert_eq!(next.accepted_steps(), rod.accepted_steps());
        assert!(std::ptr::eq(next.model(), rod.model()));
        assert_eq!(
            next.linear_momenta_kg_m_s()[0],
            rod.linear_momenta_kg_m_s()[0]
        );
        let continued = next
            .step_candidate(1e-5, &[RodNodalWrench::default(); 2], controls(1e-5), &gate)
            .unwrap()
            .0;
        assert_eq!(continued.accepted_steps(), rod.accepted_steps() + 1);
        assert!(continued.time_s() > rod.time_s());
    }
}

fn eccentric_source() -> SpatialMaterialRod {
    let source = section(1.0);
    let link = MaterialTensionLink::try_new(
        source.elastic().clone(),
        [
            TensionLinkAttachment {
                node: 0,
                offset_body_m: [0.05, 0.08, -0.03],
            },
            TensionLinkAttachment {
                node: 1,
                offset_body_m: [0.05, 0.08, -0.03],
            },
        ],
        0.9,
        1,
        0.2,
    )
    .unwrap();
    SpatialMaterialRod::try_straight(0.8, vec![source])
        .unwrap()
        .with_tension_links(vec![link], &CancelGate::new())
        .unwrap()
}

#[test]
fn eccentric_rotated_attachments_match_spatial_mass_solve_and_objectivity() {
    let gate = CancelGate::new();
    let model = Arc::new(
        eccentric_source()
            .into_spatial_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let node_inertia = model.inertia()[1];
    let center = node_inertia.center_of_mass_body_m();
    let tensor = node_inertia.inertia_com_body_kg_m2();
    let spatial = fs_mbd::articulated::SpatialInertia::new(
        node_inertia.mass_kg(),
        Vec3::new(center[0], center[1], center[2]),
        Mat3 {
            m: std::array::from_fn(|c| tensor[c / 3][c % 3]),
        },
    )
    .unwrap()
    .matrix()
    .unwrap();
    let rod_arm = [0.03, -0.06, 0.04];
    let rigid_arm = [0.1, -0.08, 0.06];
    let rod_torque = cross(rod_arm, [1.0, 0.0, 0.0]);
    let mut rod_response = [rod_torque[0], rod_torque[1], rod_torque[2], 1.0, 0.0, 0.0];
    fs_la::factor::lu(&spatial.m, 6)
        .unwrap()
        .solve(&mut rod_response);
    let rod_response_omega = [rod_response[0], rod_response[1], rod_response[2]];
    let rod_response_velocity = [rod_response[3], rod_response[4], rod_response[5]];
    let exact_rod_mobility = add(rod_response_velocity, cross(rod_response_omega, rod_arm))[0];
    let body_tensor = Mat3 {
        m: [0.07, 0.01, -0.004, 0.01, 0.09, 0.006, -0.004, 0.006, 0.08],
    };
    let mass = MassProperties::from_inertia_tensor(2.4, BodyVec::ZERO, body_tensor).unwrap();
    let mut body_response_omega = cross(rigid_arm, [1.0, 0.0, 0.0]);
    fs_la::factor::lu(&body_tensor.m, 3)
        .unwrap()
        .solve(&mut body_response_omega);
    let exact_body_mobility = 1.0 / mass.mass() + cross(body_response_omega, rigid_arm)[0];
    let rod_v = [-0.2, 0.1, -0.1];
    let rod_w = [0.2, -0.1, 0.3];
    let body_v = [1.0, -0.2, 0.15];
    let body_w = [-0.3, 0.1, 0.2];
    let speed = sub(
        add(rod_v, cross(rod_w, rod_arm)),
        add(body_v, cross(body_w, rigid_arm)),
    )[0];
    let expected_impulse = -1.4 * speed / (exact_rod_mobility + exact_body_mobility);
    let mut previous_loss = None;
    for (q, translation, boost) in [
        (So3::identity(), [0.0; 3], [0.0; 3]),
        (
            So3::exp(So3Tangent::new(Vec3::new(0.4, -0.2, 0.3))).unwrap(),
            [0.7, -0.2, 0.1],
            [0.4, -0.3, 0.2],
        ),
    ] {
        let mut pose = model.initial_pose();
        for node in 0..2 {
            pose.positions_m[node] = add(rotate(q, pose.positions_m[node]), translation);
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
        let position = add(
            rotate(
                q,
                sub(
                    add([0.8, 0.0, 0.0], rod_arm),
                    add(rigid_arm, [1e-9, 0.0, 0.0]),
                ),
            ),
            translation,
        );
        let rigid = RigidBodyState::new(
            Pose::new(vector(position), orientation(q)).unwrap(),
            vector(add(rotate(q, body_v), boost)).scale(mass.mass()),
            mass.angular_momentum_body_checked(vector(body_w)).unwrap(),
        )
        .unwrap();
        let contact = RodRigidNormalContact {
            point_body_m: rod_arm,
            rigid_arm_body_m: rigid_arm,
            normal_world: rotate(q, [1.0, 0.0, 0.0]),
            ..witness(1, 0.4)
        };
        let (next, finite, report) = rod
            .apply_rigid_normal_impact(rigid, mass, contact, impact_controls(), &gate)
            .unwrap();
        near(
            report.rod_inverse_mass_kg_inverse,
            exact_rod_mobility,
            2e-11,
        );
        near(
            report.rigid_inverse_mass_kg_inverse,
            exact_body_mobility,
            2e-12,
        );
        near(report.contact.normal_impulse, expected_impulse, 2e-12);
        vector_near(
            next.angular_velocities_body_rad_s().unwrap()[1],
            add(rod_w, rod_response_omega.map(|x| x * expected_impulse)),
            2e-11,
        );
        vector_near(
            next.velocities_m_s().unwrap()[1],
            add(
                rotate(
                    q,
                    add(rod_v, rod_response_velocity.map(|x| x * expected_impulse)),
                ),
                boost,
            ),
            2e-12,
        );
        vector_near(
            array(
                mass.angular_velocity_body_checked(finite.angular_momentum_body())
                    .unwrap(),
            ),
            sub(body_w, body_response_omega.map(|x| x * expected_impulse)),
            2e-12,
        );
        let before = rod
            .observe(impact_controls().impulse.balance_origin_m)
            .unwrap();
        let after = next
            .observe(impact_controls().impulse.balance_origin_m)
            .unwrap();
        let delta_p = sub(
            array(finite.linear_momentum_world()),
            array(rigid.linear_momentum_world()),
        );
        vector_near(
            add(
                sub(after.linear_momentum_kg_m_s, before.linear_momentum_kg_m_s),
                delta_p,
            ),
            [0.0; 3],
            3e-12,
        );
        let body_spin_change = rotate(
            q,
            sub(
                array(finite.angular_momentum_body()),
                array(rigid.angular_momentum_body()),
            ),
        );
        let body_moment_change = add(
            cross(
                sub(position, impact_controls().impulse.balance_origin_m),
                delta_p,
            ),
            body_spin_change,
        );
        vector_near(
            add(
                sub(
                    after.angular_momentum_kg_m2_s,
                    before.angular_momentum_kg_m2_s,
                ),
                body_moment_change,
            ),
            [0.0; 3],
            4e-12,
        );
        if let Some(loss) = previous_loss {
            near(report.contact.dissipated_energy, loss, 2e-12);
        }
        previous_loss = Some(report.contact.dissipated_energy);
    }
}

#[test]
fn supported_node_rotates_and_returns_the_missing_linear_reaction() {
    let gate = CancelGate::new();
    for clamped_frame in [false, true] {
        let source = SpatialMaterialRod::try_straight(0.8, vec![section(1.0)])
            .unwrap()
            .with_supports(
                &[RodSupport {
                    node: 1,
                    position_m: [Some(0.8), None, None],
                    frame: clamped_frame.then_some(So3::identity()),
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
        let rigid = body(mass, [0.8, 0.07, 0.0], [1.0, 0.0, 0.0]);
        let contact = RodRigidNormalContact {
            point_body_m: [0.0, 0.07, 0.0],
            ..witness(1, 0.5)
        };
        let (next, finite, report) = rod
            .apply_rigid_normal_impact(rigid, mass, contact, impact_controls(), &gate)
            .unwrap();
        let i_z = model.inertia()[1].principal_kg_m2()[2];
        let mobility = if clamped_frame {
            0.0
        } else {
            0.07 * 0.07 / i_z
        };
        let impulse = 1.5 / (mobility + 1.0 / mass.mass());
        near(report.rod_inverse_mass_kg_inverse, mobility, 2e-12);
        near(report.contact.normal_impulse, impulse, 2e-12);
        assert_eq!(next.linear_momenta_kg_m_s(), rod.linear_momenta_kg_m_s());
        let expected_spin = if clamped_frame { 0.0 } else { -0.07 * impulse };
        near(next.body_spins_kg_m2_s()[1][2], expected_spin, 2e-13);
        let support = &report.rod_impulse.as_ref().unwrap().support_impulses[1];
        vector_near(support.linear_n_s, [-impulse, 0.0, 0.0], 2e-12);
        let expected_couple = if clamped_frame { 0.07 * impulse } else { 0.0 };
        near(support.angular_nm_s[2], expected_couple, 2e-13);
        near(
            finite.linear_momentum_world().x - rigid.linear_momentum_world().x,
            support.linear_n_s[0],
            2e-12,
        );
        near(report.energy_residual_j, 0.0, 3e-12);
    }
}

#[test]
fn large_common_translation_preserves_small_collision_loss() {
    let gate = CancelGate::new();
    let base = free_rod([0.0; 3]);
    let m = base.model().inertia()[1].mass_kg();
    let common = 1_073_741_824.0;
    let rod = free_rod([common - 1.0 / m, 0.0, 0.0]);
    let mass = inertia(m);
    let rigid = RigidBodyState::new(
        Pose::new(BodyVec::new(0.8, 0.0, 0.0), UnitQuaternion::IDENTITY).unwrap(),
        BodyVec::new(m * common + 1.0, 0.0, 0.0),
        BodyVec::ZERO,
    )
    .unwrap();
    let mut bounds = impact_controls();
    bounds.impulse.work_tolerance_j = 1e-6;
    bounds.velocity_tolerance_m_s = 1e-6;
    let (next, finite, report) = rod
        .apply_rigid_normal_impact(rigid, mass, witness(1, 0.0), bounds, &gate)
        .unwrap();
    assert!(report.contact.dissipated_energy > 0.1);
    near(
        report.kinetic_energy_change_j,
        -report.contact.dissipated_energy,
        2e-10,
    );
    near(report.work_energy_residual_j, 0.0, 2e-10);
    let before =
        rod.observe([0.0; 3]).unwrap().kinetic_energy_j + rigid.kinetic_energy(mass).unwrap();
    let after =
        next.observe([0.0; 3]).unwrap().kinetic_energy_j + finite.kinetic_energy(mass).unwrap();
    assert!(
        ((after - before) - report.kinetic_energy_change_j).abs()
            > 100.0 * report.contact.dissipated_energy,
        "absolute endpoint energy subtraction cannot resolve this collision"
    );
}

#[test]
fn inactive_separating_and_refused_events_preserve_both_participants() {
    let gate = CancelGate::new();
    let rod = free_rod([-0.6, 0.0, 0.0]);
    let mass = inertia(2.3);
    let rigid = body(mass, [0.8, 0.0, 0.0], [1.7, 0.0, 0.0]);
    for other in [
        body(mass, [0.79, 0.0, 0.0], [1.7, 0.0, 0.0]),
        body(mass, [0.8, 0.0, 0.0], [-1.7, 0.0, 0.0]),
    ] {
        let (next, finite, report) = rod
            .apply_rigid_normal_impact(other, mass, witness(1, 0.4), impact_controls(), &gate)
            .unwrap();
        assert!(!report.contact.impacted);
        assert!(report.rod_impulse.is_none() && report.rigid_impulse.is_none());
        assert_eq!(next.linear_momenta_kg_m_s(), rod.linear_momenta_kg_m_s());
        assert_eq!(next.body_spins_kg_m2_s(), rod.body_spins_kg_m2_s());
        assert_eq!(next.pose(), rod.pose());
        assert_eq!(finite, other);
    }
    let before = rod.observe([0.0; 3]).unwrap();
    let penetrating = body(mass, [0.80001, 0.0, 0.0], [1.7, 0.0, 0.0]);
    assert!(matches!(
        rod.apply_rigid_normal_impact(penetrating, mass, witness(1, 0.4), impact_controls(), &gate),
        Err(RodRigidImpactError::Contact(_))
    ));
    let crooked = body(mass, [0.8, 1e-4, 0.0], [1.7, 0.0, 0.0]);
    let mut bounds = impact_controls();
    bounds.geometry_tolerance_m = 1e-3;
    // Both native momentum jumps stage successfully, but their actual moment
    // fails the joint angular-momentum gate. A valid retry sees the old states.
    assert!(matches!(
        rod.apply_rigid_normal_impact(crooked, mass, witness(1, 0.4), bounds, &gate),
        Err(RodRigidImpactError::Balance { .. })
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        rod.apply_rigid_normal_impact(rigid, mass, witness(1, 0.4), bounds, &cancelled)
            .is_err()
    );
    assert_eq!(
        rod.observe([0.0; 3]).unwrap().kinetic_energy_j,
        before.kinetic_energy_j
    );
    let first = rod
        .apply_rigid_normal_impact(rigid, mass, witness(1, 0.4), impact_controls(), &gate)
        .unwrap();
    let retry = rod
        .apply_rigid_normal_impact(rigid, mass, witness(1, 0.4), impact_controls(), &gate)
        .unwrap();
    assert_eq!(
        first.0.linear_momenta_kg_m_s(),
        retry.0.linear_momenta_kg_m_s()
    );
    assert_eq!(first.1, retry.1);
}

#[test]
fn rotational_approach_cannot_hide_unresolved_linear_impulse_in_large_drift() {
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
    let contact = RodRigidNormalContact {
        point_body_m: [0.0, 0.1, 0.0],
        ..witness(1, 0.4)
    };
    assert!(matches!(
        rod.apply_rigid_normal_impact(rigid, mass, contact, impact_controls(), &gate),
        Err(RodRigidImpactError::Rod(
            RodMotionError::ImpulseBalance { .. }
        ))
    ));
    assert_eq!(
        rod.body_spins_kg_m2_s(),
        base.model()
            .inertia()
            .iter()
            .map(|i| [0.0, 0.0, i.principal_kg_m2()[2]])
            .collect::<Vec<_>>()
    );
}
